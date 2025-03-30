#ifndef Player_hpp
#define Player_hpp
#define _USE_MATH_DEFINES
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

using namespace std;

//double freq_bin[] = { 19.0, 140.0, 400.0, 2600.0, 5200.0, nyquist }

typedef struct {
	double* in;
	fftw_complex* out;
	double* magnitude;
	fftw_plan p;
	int startIdx;
	int specSize;
} streamData;

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
	int bands = 42;

	void init();
	void eventCall(SDL_Event event);
	void update(double delta);
	void draw(double delta);
	void clean();
	void close();
	void playSound();
	void stopSound();
	bool loadMP3(string path);
	bool loadWAV(string path);
	float lerp(float a, float b, float f);
	double logint(double a, double b, double f);
	double clamp(double a, double b, double c);
	double toDb(double a);
	const double DB_MIN = 20.0f * std::log10(std::numeric_limits<float>::min());
	
	bool running() { return isRunning; };

private:
	bool bully2 = false;
	bool isRunning = true;
	SDL_Window* window;
	SDL_Renderer* renderer;
	//vector<double> magnitudes;
	vector<SDL_FRect> rectangles;
	SDL_AudioDeviceID audioDevice;
	//double* freqBin;
	double* bandIndices;
	short* bandWidths;
	//double* magnitudes;
	double* decibelsBuf;
	double* decibels;
	//double* normalSpectrum;
	//int* sizes;
	SDL_AudioStream* audioStream;
	//int16_t* audBufferTmp;
	void* audBuffer;
	int audSize;
	double binSpacing;
	
	double sampleRatio;
};

#endif