#ifndef Player_hpp
#define Player_hpp
#define _USE_MATH_DEFINES
#define FREQ_START 24
#define FREQ_END 16000
#define FFT_SIZE 1024

#define BAR_WIDTH 6.0
#define BAR_SPACE 3.0

#define MINIMP3_IMPLEMENTATION
//#define MINIMP3_ALLOW_MONO_STEREO_TRANSITION

#include "SDL3/SDL.h"
#include <stdio.h>
#include <vector>
#include "fftw3.h"
#include <string>
#include <iostream>
#include <math.h>
#include <algorithm>
#include "FFTProcessor.hpp"

using namespace std;

//double freq_bin[] = { 19.0, 140.0, 400.0, 2600.0, 5200.0, nyquist }

class Player
{
public:
	Player(SDL_Window* win, SDL_Renderer* ren);
	~Player();

	static bool mustCall;
	static int len;
	static float toLog;

	//static Player& instance;
	static streamData* specData;
	static void SDLCALL onAudioData(void* userdata, const SDL_AudioSpec* spec, float* buffer, int buflen);
	static SDL_AudioSpec audioSpec;
	static int frameCount;
	static int currentFrameSize;
	static int bytesElapsed;
	uint16_t bands = 64;

	void init();
	void eventCall(SDL_Event* event);
	void update(double delta);
	void draw(double delta);
	void clean();
	void close();
	void playSound();
	void stopSound();
	bool loadMP3(string path);
	bool loadWAV(string path);
	
	bool running() { return isRunning; };

private:
	FFTProcessor* processor;
	SDL_Window* window;
	SDL_Renderer* renderer;
	//vector<double> magnitudes;
	vector<SDL_FRect> rectangles;
	SDL_AudioDeviceID audioDevice;
	//double* freqBin;
	double* bins;
	double* bandFreqs;
	//double* magnitudes;
	//double* decibelsBuf;
	//double* decibels;
	double* spectrum;
	//double* normalSpectrum;
	//int* sizes;
	SDL_AudioStream* audioStream;
	//int16_t* audBufferTmp;
	void* audBuffer;
	int audSize;
	double sampleRatio;
	double sampleRate;

	bool isRunning = true;
};

#endif