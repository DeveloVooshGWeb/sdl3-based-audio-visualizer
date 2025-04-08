#include "FFTProcessor.hpp"

FFTProcessor::FFTProcessor(size_t sr, uint16_t ffts, uint16_t bs, uint16_t fStart, uint16_t fEnd)
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
	_freqStart = fStart;
	_freqEnd = fEnd;
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

void FFTProcessor::assign(double* bandFreqs, float t, double emaValue)
{
	_bandFreqs = bandFreqs;
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
		double u = _bandFreqs[i] - floor(_bandFreqs[i]);
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
	//_genericFunction(output);
	//_interpolateCatmullRom();
	_extractHighestMagnitudes(_curBands);
	//_interpolateCatmullRom(_curBands);
	_interpolateLinear(_curBands);
	//copy(_curBands, _curBands + _bands, output);
	_smoothingExponentialMovingAverage(output);
}

void FFTProcessor::_extractHighestMagnitudes(double* output)
{
	double coeff = _sampleRate / _fftSize;
	double idx, len, mag;
	int cIdx, pIdx = -1;
	for (uint16_t i = 0; i < _bands; i++)
	{
		output[i] = -DBL_MAX;
		idx = _bandFreqs[i] / coeff;
		cIdx = (int)idx;
		if (cIdx == pIdx) continue;
		len = (_bandFreqs[i + 1] / coeff) - idx;
		mag = *max_element(_bins + cIdx, _bins + (int)(idx + len));
		if (mag > output[i]) output[i] = mag;
		//cout << "DEBUG INFO: " << i << " | " << _curBands[i] << endl;
		pIdx = cIdx;
	}
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
			double idx = _bandFreqs[j];
			double las = 19050;
			if (j < _bands - 1) las = _bandFreqs[j + 1];
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

double getKnotInterval(double a, double* p0, double* p1)
{
	double x = p1[0] - p0[0];
	double y = p1[1] - p0[1];
	return pow(x * x + y * y, 0.5 * a);
}

static double lerp(double a, double b, double f)
{
	return a * (1.0 - f) + (b * f);
}

static double catmullRom(double* p0, double* p1, double* p2, double* p3, double x, double a = 0.5)
{
	double t0 = 0.0;
	double t1 = getKnotInterval(a, p0, p1);
	double t2 = getKnotInterval(a, p1, p2) + t1;
	double t3 = getKnotInterval(a, p2, p3) + t2;
	//t = lerp(t1, t2, t);
	double t = (x - p0[0]) / (p3[0] - p0[0]);
	//cout << "t: " << t << endl;
	double u = lerp(t1, t2, t);
	double A1 = (t1 - u) / (t1 - t0) * p0[1] + (u - t0) / (t1 - t0) * p1[1];
	double A2 = (t2 - u) / (t2 - t1) * p1[1] + (u - t1) / (t2 - t1) * p2[1];
	double A3 = (t3 - u) / (t3 - t2) * p2[1] + (u - t2) / (t3 - t2) * p3[1];
	double B1 = (t2 - u) / (t2 - t0) * A1 + (u - t0) / (t2 - t0) * A2;
	double B2 = (t3 - u) / (t3 - t1) * A2 + (u - t1) / (t3 - t1) * A3;
	//double A1 = lerp(p0, p1, (u - t0) / (t1 - t0));
	//double A2 = lerp(p1, p2, (u - t1) / (t2 - t1));
	//double A3 = lerp(p2, p3, (u - t2) / (t3 - t2));
	//double B1 = lerp(A1, A2, (u - t0) / (t2 - t0));
	//double B2 = lerp(A2, A3, (u - t1) / (t3 - t1));
	double C = (t2 - u) / (t2 - t1) * B1 + (u - t1) / (t2 - t1) * B2;
	return C;
}

void FFTProcessor::_interpolateCatmullRom(double* output)
{
	double mag, lastFreq, lastMag, curFreq, cRealFreq;
	uint16_t lastU;
	uint16_t u = 0;
	double points[4][2] = { { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 }, { 0.0, 0.0 } };
	uint16_t* realFreq = (uint16_t*)malloc(_bands * sizeof(uint16_t));
	double* realMag = (double*)malloc(_bands * sizeof(double));
	uint16_t realSz = 0;
	for (uint16_t i = 0; i < _bands; i++)
	{
		mag = output[i];
		if (mag <= -DBL_MAX) continue;
		realFreq[realSz] = _bandFreqs[i];
		realMag[realSz] = mag;
		realSz++;
	}
	for (uint16_t i = 0; i < _bands; i++)
	{
		
		mag = output[i];
		curFreq = _bandFreqs[i];
		cRealFreq = realFreq[u];
		if (u == 0)
		{
			points[0][0] = -realFreq[0];
			points[0][1] = 0.0;
		}
		else
		{
			points[0][0] = realFreq[u - 1];
			points[0][1] = realMag[u - 1];
		}
		cRealFreq = realFreq[u];
		points[1][0] = cRealFreq;
		points[1][1] = realMag[u];
		lastU = u < realSz - 1 ? u + 1 : u;
		points[2][0] = realFreq[lastU];
		points[2][1] = realMag[lastU];
		points[3][0] = realFreq[u < realSz - 2 ? u + 2 : lastU];
		points[3][1] = realMag[u < realSz - 2 ? u + 2 : lastU];
		double curvePoint = catmullRom(points[0], points[1], points[2], points[3], curFreq, 0.5);
		output[i] = curvePoint;
		if (i > 0.0 && mag > -DBL_MAX)
		{
			u++;
		}
	}
	free(realFreq);
	free(realMag);
}

void FFTProcessor::_interpolateLinear(double* output)
{
	double x0, x1, y0, y1, curFreq, mag;
	uint16_t u = 0;
	uint16_t* realFreq = (uint16_t*)malloc(_bands * sizeof(uint16_t));
	double* realMag = (double*)malloc(_bands * sizeof(double));
	uint16_t realSz = 0;
	for (uint16_t i = 0; i < _bands; i++)
	{
		mag = output[i];
		if (mag <= -DBL_MAX) continue;
		realFreq[realSz] = _bandFreqs[i];
		realMag[realSz] = mag;
		realSz++;
	}
	for (uint16_t i = 0; i < _bands; i++)
	{
		mag = output[i];
		if (mag <= -DBL_MAX)
		{
			curFreq = _bandFreqs[i];
			double z = curFreq;
			double w = floor(z);
			x0 = realFreq[u];
			x1 = realFreq[min(realSz - 1, u + 1)];
			y0 = realMag[u];
			y1 = realMag[min(realSz - 1, u + 1)];
			output[i] = y0 + ((y1 - y0) * ((z - x0) / (x1 - x0)));
		}
		else if (i > 0)
		{
			u++;
		}
	}
	free(realFreq);
	free(realMag);
}

uint16_t FFTProcessor::getBinSize()
{
	return _binSize;
}