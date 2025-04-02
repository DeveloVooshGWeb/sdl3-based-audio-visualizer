#include "FFTProcessor.hpp"

FFTProcessor::FFTProcessor(size_t sr, uint16_t ffts, uint16_t bs)
{
	_sampleRate = sr;
	_fftSize = ffts;
	_binSize = _fftSize / 2.0;
	_bands = bs;
	_prevBands = (double*)malloc(_bands * sizeof(double));
	for (int i = 0; i < _bands; i++)
	{
		_prevBands[i] = 0.0;
	}
	_splineWeights = (float*)malloc(_binSize * 4 * sizeof(float));
	_curBands = (double*)malloc(_bands * sizeof(double));
}

FFTProcessor::~FFTProcessor()
{
	free(_prevBands);
	free(_splineWeights);
}

double FFTProcessor::lerp(double a, double b, double f)
{
	return a * (1.0 - f) + (b * f);
}

double FFTProcessor::logint(double a, double b, double f)
{
	return a * pow(b / a, f);
}

double FFTProcessor::clamp(double a, double b, double c)
{
	return max(min(a, c), b);
}

double FFTProcessor::toDb(double a)
{
	if (a > 0.0) return 20.0 * log10(a);
	return DB_MIN;
}

void FFTProcessor::assign(double* bandIndices, double* bandWidths, float t, double emaValue)
{
	_bandIndices = bandIndices;
	_bandWidths = bandWidths;
	//_generateCatmullRomWeights(t);
	_emaValue = emaValue;
}

void FFTProcessor::_generateCatmullRomWeights(float t)
{
	float matrix[4][4] = {
		{ 0, -t, 2 * t, -t },
		{ 1, 0, t - 3, 2 - t },
		{ 0, t, 3 - (2 * t), t - 2 },
		{ 0, 0, -t, t }
	};

	for (uint16_t i = 0; i < _bands; i++)
	{
		double u = _bandIndices[i] - floor(_bandIndices[i]);
		double row[4] = { 1, u, u * u, u * u * u };
		for (size_t j = 0; j < 4; j++)
		{
			float sum = 0;
			for (size_t k = 0; k < 4; k++) sum += row[k] * matrix[j][k];
			_splineWeights[(i * 4) + j] = sum;
		}
	}
}

void FFTProcessor::pipe(double* bins, double* output)
{
	_bins = bins;
	//_transformToDb();
	//_smoothingExponentialMovingAverage();
	_genericFunction(output);
	//_interpolateCatmullRom();
	_smoothingExponentialMovingAverage(output);
}

void FFTProcessor::_genericFunction(double* output)
{
	for (uint16_t j = 0; j < _bands; j++)
	{
		_curBands[j] = 0.0;
	}
	for (uint16_t i = 0; i < _binSize; i++)
	{
		double mag = _bins[i];
		double freq = (double)i * (double)_sampleRate / (double)_fftSize;
		for (uint16_t j = 0; j < _bands; j++)
		{
			double idx = _bandIndices[j];
			double las = 19050;
			if (j < _bands - 1) las = _bandIndices[j + 1];
			if (freq >= idx && freq < las)
			{
				if (mag > _curBands[j]) _curBands[j] = mag;
			}
		}
	}
}

void FFTProcessor::_transformToDb()
{
	for (size_t i = 0; i < _binSize; i++)
	{
		_bins[i] = toDb(_bins[i]);
	}
}

void FFTProcessor::_smoothingExponentialMovingAverage(double* output)
{
	for (size_t i = 0; i < _bands; i++)
	{
		output[i] = _emaValue * _curBands[i] + (1.0 - _emaValue) * _prevBands[i];
		_prevBands[i] = output[i];
	}
}

#ifdef _MSC_VER
#define WAV_FORCE_INLINE __forceinline
#else
#define WAV_FORCE_INLINE __attribute__((always_inline)) inline
#endif

/*
static WAV_FORCE_INLINE double horizontal_sum(__m128 vec)
{
	auto low = vec;
	auto high = _mm_shuffle_ps(low, low, _MM_SHUFFLE(3, 2, 3, 2));  // high[0] = low[2], high[1] = low[3]
	low = _mm_add_ps(high, low);                                    // (h[0] + l[0]) (h[1] + l[1])
	high = _mm_movehdup_ps(low);                                    // high[0] = low[1]
	return _mm_cvtss_f32(_mm_add_ss(high, low));
}

static WAV_FORCE_INLINE double horizontal_sum(__m256 vec)
{
	return horizontal_sum(_mm_add_ps(_mm256_extractf128_ps(vec, 1), _mm256_castps256_ps128(vec)));
}
*/
/*
static WAV_FORCE_INLINE double horizontal_sum(__m128 vec) {
	auto low = vec;
	auto high = _mm_shuffle_ps(low, low, _MM_SHUFFLE(3, 2, 3, 2));  // high[0] = low[2], high[1] = low[3]
	low = _mm_add_ps(high, low);                                     // (h[0] + l[0]) (h[1] + l[1])
	high = _mm_movehdup_ps(low);                                     // high[0] = low[1]
	return _mm_cvtss_f32(_mm_add_ss(high, low));
}

static WAV_FORCE_INLINE double horizontal_sum(__m256d vec) {
	__m128d low = _mm256_castpd256_pd128(vec);
	__m128d high = _mm256_extractf128_pd(vec, 1);
	__m128d sum = _mm_add_pd(low, high);
	sum = _mm_hadd_pd(sum, sum);
	return _mm_cvtsd_f64(sum);
}
*/

static WAV_FORCE_INLINE double horizontal_sum(__m256d vec) {
	__m128d low = _mm256_castpd256_pd128(vec);
	__m128d high = _mm256_extractf128_pd(vec, 1);
	__m128d sum = _mm_add_pd(low, high);
	sum = _mm_hadd_pd(sum, sum);
	return _mm_cvtsd_f64(sum);
}

double getKnotInterval(double alpha, double p0, double p1, double y0, double y1)
{
	double p2 = p0 - p1;
	double y2 = y0 - y1;
	return pow(p2 * p2 + y2 * y2, 0.5 * alpha);
}

static double lerp(double a, double b, double f)
{
	return a * (1.0 - f) + (b * f);
}

static double catmullRom(double p0, double p1, double p2, double p3, double t, double alpha = 0.5)
{
	double t0 = 0.0;
	double t1 = 0.5;
	double t2 = 1.0;
	double t3 = 1.5;
	//t = lerp(t1, t2, t);
	double u = lerp(t1, t2, t);
	double A1 = (t1 - u) / (t1 - t0) * p0 + (u - t0) / (t1 - t0) * p1;
	double A2 = (t2 - u) / (t2 - t1) * p1 + (u - t1) / (t2 - t1) * p2;
	double A3 = (t3 - u) / (t3 - t2) * p2 + (u - t2) / (t3 - t2) * p3;
	double B1 = (t2 - u) / (t2 - t0) * A1 + (u - t0) / (t2 - t0) * A2;
	double B2 = (t3 - u) / (t3 - t1) * A2 + (u - t1) / (t3 - t1) * A3;
	//double A1 = lerp(p0, p1, (u - t0) / (t1 - t0));
	//double A2 = lerp(p1, p2, (u - t1) / (t2 - t1));
	//double A3 = lerp(p2, p3, (u - t2) / (t3 - t2));
	//double B1 = lerp(A1, A2, (u - t0) / (t2 - t0));
	//double B2 = lerp(A2, A3, (u - t1) / (t3 - t1));
	double C = lerp(B1, B2, (u - t1) / (t2 - t1));
	return C;
}

void FFTProcessor::_interpolateCatmullRom()
{
	double idx;
	uint16_t idxI;
	double t;
	double len;
	double sum;
	double points[] = {0.0, 0.0, 0.0, 0.0};
	for (uint16_t i = 0; i < _bands; i++)
	{
		_curBands[i] = 0.0;
		idx = _bandIndices[i];
		len = _bandWidths[i];
		t = idx - floor(idx);
		//sum = 0.0;
		for (uint16_t j = 0; j < len; j++, idx += 1.0)
		{
			idxI = (uint16_t)idx;
			if (idxI < (int)_bandIndices[0]) {
				points[0] = 0.0;
			}
			else
			{
				points[0] = _bins[idxI - 1];
			}
			points[1] = _bins[idxI];
			/*
			if (idxI >= _binSize - 2)
			{
				if (idxI >= _binSize - 1)
				{
					points[2] = _bins[idxI];
				}
				else
				{
					points[2] = _bins[idxI + 1];
				}
				points[3] = points[2];
			}
			else
			{
				points[2] = points[1];
				points[3] = points[2];
			}*/
			points[2] = _bins[min(_binSize - 1, idxI + 1)];
			points[3] = _bins[min(_binSize - 1, idxI + 2)];
			double rom = catmullRom(points[0], points[1], points[2], points[3], t, 0.5);
			if (rom > _curBands[i]) _curBands[i] = rom;
		}
		//_curBands[i] = sum / len;
	}
	/*
	constexpr auto step = sizeof(__m128) / sizeof(float);
	const auto sse_stop = (intmax_t)_binSize - 2;
	const auto bands = (intmax_t)_bands;
	for (intmax_t i = 0, k = 0, l = 0; i < bands; ++i)
	{
		auto vecsum = _mm_setzero_ps();
		auto count = (intmax_t)_bandWidths[i];
		for (intmax_t j = 0; j < count; ++j, ++k, l += step)
		{
			auto index = (intmax_t)_bandIndices[k];
			if ((index >= 1) && (index < sse_stop))
			{
				const float bin = static_cast<const float>(_bins[index - 1]);
				vecsum = _mm_fmadd_ps(_mm_loadu_ps(&bin), _mm_load_ps(&_splineWeights[l]), vecsum);
			}
			else
			{
				
				const auto start = index - 1;
				const auto stop = std::min(index + 3, (intmax_t)_binSize);

				auto sum = _mm_setzero_ps();
				for (auto m = std::max(start, (intmax_t)0); m < stop; ++m)
				{
					const float bin = static_cast<const float>(_bins[m]);
					sum = _mm_fmadd_ss(_mm_load_ss(&bin), _mm_load_ss(&_splineWeights[l + (m - start)]), sum);
				}
				vecsum = _mm_add_ss(vecsum, sum);
			}
		}
		output[i] = horizontal_sum(vecsum) / count;
	}
	*/
	/*
	constexpr auto step = sizeof(__m128) / sizeof(float);  // 4 floats per __m128
	const auto sse_stop = static_cast<intmax_t>(_binSize) - 2;
	const auto bands = static_cast<intmax_t>(_bands);

	for (intmax_t i = 0, k = 0, l = 0; i < bands; ++i) {
		__m256d vecsum = _mm256_setzero_pd();  // Use double-precision AVX
		auto count = static_cast<intmax_t>(_bandWidths[i]);

		for (intmax_t j = 0; j < count; ++j, ++k, l += step) {
			auto index = static_cast<intmax_t>(_bandIndices[k]);

			if (index >= 1 && index < sse_stop) {
				// Load double from _bins, convert to float for weights
				__m256d bin_vec = _mm256_set1_pd(_bins[index - 1]);
				__m128 weights = _mm_load_ps(&_splineWeights[l]);
				__m256d weights_d = _mm256_cvtps_pd(weights);  // Convert weights to double
				vecsum = _mm256_fmadd_pd(bin_vec, weights_d, vecsum);
			}
			else {
				// Edge case: handle boundary conditions
				const auto start = index - 1;
				const auto stop = std::min(index + 3, static_cast<intmax_t>(_binSize));

				__m256d sum = _mm256_setzero_pd();
				for (auto m = std::max(start, static_cast<intmax_t>(0)); m < stop; ++m) {
					__m256d bin = _mm256_set1_pd(_bins[m]);
					__m128 weight = _mm_load_ss(&_splineWeights[l + (m - start)]);
					__m256d weight_d = _mm256_cvtps_pd(weight);  // Convert weight to double
					sum = _mm256_fmadd_pd(bin, weight_d, sum);
				}
				vecsum = _mm256_add_pd(vecsum, sum);
			}
		}
		output[i] = horizontal_sum(vecsum) / count;
	}
	*/
	/*
	const auto sse_stop = static_cast<intmax_t>(_binSize) - 2;
	const auto bands = static_cast<intmax_t>(_bands);
	for (intmax_t i = 0, k = 0, l = 0; i < bands; ++i) {
		float sum = 0.0f;
		auto count = static_cast<intmax_t>(_bandWidths[i]);
		for (intmax_t j = 0; j < count; ++j, ++k, l += 4) {
			auto index = static_cast<intmax_t>(_bandIndices[k]);
			if (index >= 1 && index < sse_stop) {
				// Valid index: use the same bin for all 4 weights
				const float bin = static_cast<const float>(_bins[index - 1]);
				for (int m = 0; m < 4; ++m) {
					sum += bin * _splineWeights[l + m];
				}
			}
			else {
				// Handle edge cases: each weight corresponds to a different bin
				const auto start = index - 1;
				const auto stop = std::min(index + 3, static_cast<intmax_t>(_binSize));
				for (auto m = std::max(start, static_cast<intmax_t>(0)); m < stop; ++m) {
					const float bin = static_cast<const float>(_bins[m]);
					const intmax_t weight_idx = l + (m - start);
					sum += bin * _splineWeights[weight_idx];
				}
			}
		}
		output[i] = sum / count;
	}
	*/
	/*
	constexpr int step = 4;  // 4 weights per bin (Catmull-Rom)
	const auto binSize = static_cast<intmax_t>(_binSize);
	const auto bands = static_cast<intmax_t>(_bands);

	for (intmax_t i = 0, k = 0, l = 0; i < bands; ++i) {
		__m256d sum_vec = _mm256_setzero_pd();
		const auto count = static_cast<intmax_t>(_bandWidths[i]);

		for (intmax_t j = 0; j < count; ++j, ++k, l += step) {
			const intmax_t index = static_cast<intmax_t>(_bandIndices[k]);
			const intmax_t start = index - 1;  // Catmull-Rom uses [i-1, i, i+1, i+2]

			// Load 4 weights (precomputed Catmull-Rom coefficients)
			__m256d weights = _mm256_cvtps_pd(_mm_load_ps(&_splineWeights[l]));

			// Load 4 bins: _bins[start], _bins[start+1], _bins[start+2], _bins[start+3]
			__m256d bins = _mm256_set_pd(
				(start + 3 < binSize) ? _bins[start + 3] : 0.0,
				(start + 2 < binSize) ? _bins[start + 2] : 0.0,
				(start + 1 < binSize) ? _bins[start + 1] : 0.0,
				(start >= 0) ? _bins[start] : 0.0
			);

			// sum_vec += bins * weights (Catmull-Rom interpolation)
			sum_vec = _mm256_fmadd_pd(bins, weights, sum_vec);
		}

		output[i] = horizontal_sum(sum_vec) / count;
	}
	*/

	/*
	double sum;
	for (size_t i = 0, k = 0, l = 0; i < _bands; i++)
	{
		sum = 0.0;
		uint16_t bw = _bandWidths[i];
		for (uint16_t j = 0; j < bw; j++, k++, l += 4)
		{
			uint16_t idx = (uint16_t)(_bandIndices[k]);
			if (idx >= 1 && idx < _binSize - 2)
			{
				for (size_t m = 0; m < 4; m++)
				{
					sum += _bins[idx - 1 + m] * _splineWeights[l + m];
				}
			}
			else
			{
				uint16_t beg = idx - 1;
				uint16_t end = min(idx + 3, (int)_binSize);
				for (size_t m = max((int)beg, 0); m < end; m++)
				{
					sum += _bins[m] * _splineWeights[l + (m - beg)];
				}
			}
		}
		output[i] = sum / bw;
	}
	*/
	
	/*
	unsigned long long step = sizeof(__m128) / sizeof(float);
	intmax_t sse_stop = (intmax_t)_binSize - 2;
	intmax_t bands = _bands;
	for (intmax_t i = 0, k = 0, l = 0; i < bands; ++i)
	{
		auto vecsum = _mm_setzero_ps();
		auto count = (intmax_t)_bandWidths[i];
		for (intmax_t j = 0; j < count; ++j, ++k, l += step)
		{
			auto index = (intmax_t)_bandIndices[k];
			if ((index >= 1) && (index < sse_stop))
				vecsum = _mm_fmadd_ps(_mm_loadu_ps(reinterpret_cast<const float*>(&_bins[index - 1])), _mm_load_ps(&_splineWeights[l]), vecsum);
			else
			{
				const auto start = index - 1;
				const auto stop = min(index + 3, (intmax_t)_binSize);

				auto sum = _mm_setzero_ps();
				for (auto m = max(start, (intmax_t)0); m < stop; ++m)
					sum = _mm_fmadd_ss(_mm_load_ss(reinterpret_cast<const float*>(&_bins[m])), _mm_load_ss(&_splineWeights[l + (m - start)]), sum);
				vecsum = _mm_add_ss(vecsum, sum);
			}
		}
		output[i] = horizontal_sum(vecsum) / count;
	}
	*/
}

uint16_t FFTProcessor::getBinSize()
{
	return _binSize;
}