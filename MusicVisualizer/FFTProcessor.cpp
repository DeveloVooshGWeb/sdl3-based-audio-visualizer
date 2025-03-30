#include "FFTProcessor.hpp"

FFTProcessor::FFTProcessor(size_t sr, uint16_t ffts, uint16_t bs)
{
	_sampleRate = sr;
	_fftSize = ffts;
	_binSize = _fftSize / 2.0;
	_bands = bs;
	_prevBands = (double*)malloc(_bands * sizeof(double));
	_splineWeights = (float*)malloc(_binSize * 4 * sizeof(float));
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

void FFTProcessor::assign(double* bandIndices, uint16_t* bandWidths, float t, double emaValue)
{
	_bandIndices = bandIndices;
	_bandWidths = bandWidths;
	_generateCatmullRomWeights(t);
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
	
	_interpolateCatmullRom(output);
	//_smoothingExponentialMovingAverage(output);
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
		output[i] = _emaValue * output[i] + (1.0 - _emaValue) * _prevBands[i];
		_prevBands[i] = output[i];
	}
}

#ifdef _MSC_VER
#define WAV_FORCE_INLINE __forceinline
#else
#define WAV_FORCE_INLINE __attribute__((always_inline)) inline
#endif

static WAV_FORCE_INLINE float horizontal_sum(__m128 vec)
{
	auto low = vec;
	auto high = _mm_shuffle_ps(low, low, _MM_SHUFFLE(3, 2, 3, 2));  // high[0] = low[2], high[1] = low[3]
	low = _mm_add_ps(high, low);                                    // (h[0] + l[0]) (h[1] + l[1])
	high = _mm_movehdup_ps(low);                                    // high[0] = low[1]
	return _mm_cvtss_f32(_mm_add_ss(high, low));
}

static WAV_FORCE_INLINE float horizontal_sum(__m256 vec)
{
	return horizontal_sum(_mm_add_ps(_mm256_extractf128_ps(vec, 1), _mm256_castps256_ps128(vec)));
}


void FFTProcessor::_interpolateCatmullRom(double* output)
{
	double sum;
	for (size_t i = 0, k = 0, l = 0; i < _bands; i++)
	{
		sum = 0;
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
		output[i] = toDb(sum / bw);
	}
	
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