#ifndef FFTProcessor_hpp
#define FFTProcessor_hpp

#include <stdio.h>
#include <vector>
#include "fftw3.h"
#include <algorithm>
#include <xmmintrin.h>
#include <immintrin.h>

using namespace std;

class FFTProcessor
{
public:
	FFTProcessor(size_t sr, uint16_t ffts, uint16_t bs);
	~FFTProcessor();

	double lerp(double a, double b, double f);
	double logint(double a, double b, double f);
	double clamp(double a, double b, double c);
	double toDb(double a);
	void assign(double* bandIndices, uint16_t* bandWidths, float t, double emaValue);
	void pipe(double* bins, double* output);
	uint16_t getBinSize();

	const double DB_MIN = 20.0f * log10(numeric_limits<float>::min());
	
private:
	void _generateCatmullRomWeights(float t);
	void _transformToDb();
	void _smoothingExponentialMovingAverage(double* output);
	void _interpolateCatmullRom(double* output);

	size_t _sampleRate;
	uint16_t _fftSize;
	uint16_t _binSize;
	uint16_t _bands;
	double* _prevBands;
	double* _bins;
	double* _bandIndices;
	uint16_t* _bandWidths;
	float* _splineWeights;
	double _emaValue;
};

#endif