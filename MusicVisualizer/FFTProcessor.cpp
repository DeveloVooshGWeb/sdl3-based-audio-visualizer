#include "FFTProcessor.hpp"

FFTProcessor::FFTProcessor(size_t sr, uint16_t ffts, uint16_t bs)
{
	sampleRate = sr;
	fftSize = ffts;
	bands = bs;
}

FFTProcessor::~FFTProcessor()
{
}

