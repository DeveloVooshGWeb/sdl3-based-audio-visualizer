#ifndef Exporter_hpp
#define Exporter_hpp

#include <iostream>
#include "SDL3/SDL.h"
#include "fftw3.h"
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <vector>
#include "MP4Encoder.hpp"

using namespace std;

class Exporter
{
	public:
		Exporter();
		~Exporter();

	private:
		MP4Encoder _encoder;
};

#endif