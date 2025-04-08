#ifndef Sound_hpp
#define Sound_hpp

#include "SDL3/SDL.h"
#include <stdio.h>
#include <vector>
#include "fftw3.h"
#include <string>

typedef struct {
	double* in;
	double* out;
	fftw_plan p;
	int startIdx;
	int specSize;
} streamData;

class Sound
{
public:

};

#endif