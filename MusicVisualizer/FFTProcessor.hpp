#ifndef FFTProcessor_hpp
#define FFTProcessor_hpp

#include <stdio.h>
#include <vector>
#include "fftw3.h"
#include <string>

using namespace std;

class FFTProcessor
{
public:
	FFTProcessor(size_t sr, uint16_t ffts, uint16_t bs);
	~FFTProcessor();

	float lerp(float a, float b, float f);
	double logint(double a, double b, double f);
	double clamp(double a, double b, double c);
	double toDb(double a);


private:
	size_t sampleRate;
	uint16_t fftSize;
	uint16_t bands;
};

#endif