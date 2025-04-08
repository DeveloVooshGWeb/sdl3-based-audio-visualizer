#ifndef FFTProcessor_hpp
#define FFTProcessor_hpp

#include <stdio.h>
#include <vector>
#include "fftw3.h"
#include <algorithm>
#include <xmmintrin.h>
#include <immintrin.h>
#include <iostream>

using namespace std;

typedef struct {
	double* in;
	fftw_complex* out;
	double* magnitude;
	fftw_plan p;
	int startIdx;
	int specSize;
} streamData;

class FFTProcessor
{
public:
	FFTProcessor(size_t sr, uint16_t ffts, uint16_t bs, uint16_t fStart, uint16_t fEnd);
	~FFTProcessor();

	double lerp(double a, double b, double f);
	double logint(double a, double b, double f);
	double clamp(double a, double b, double c);
	double toDb(double a);
	void assign(double* bandFreqs, float t, double emaValue);
	void pipe(double* bins, double* output);
	uint16_t getBinSize();

	const double DB_MIN = 20.0f * log10(numeric_limits<float>::min());
	
private:
	void _generateCatmullRomWeights(float t);
	void _transformToDb();
	void _smoothingExponentialMovingAverage(double* output);
	void _interpolateCatmullRom(double* output);
	void _extractHighestMagnitudes(double* output);
	void _genericFunction(double* output);
	void _interpolateLinear(double* output);

	size_t _sampleRate;
	uint16_t _fftSize;
	uint16_t _binSize;
	uint16_t _bands;
	double* _prevBands;
	double* _bins;
	double* _bandFreqs;
	float* _splineWeights;
	double _emaValue;
	double* _curBands;
	uint16_t _freqStart;
	uint16_t _freqEnd;
};

#endif