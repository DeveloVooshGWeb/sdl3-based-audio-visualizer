#ifndef Exporter_hpp
#define Exporter_hpp
#define _USE_MATH_DEFINES
#define FREQ_START 24
#define FREQ_END 16000
#define FFT_SIZE 1024

#define BAR_WIDTH 6.0
#define BAR_SPACE 3.0

#define VBR 8000000
#define ABR 192000
#define GOP_SIZE 12

#define MINIMP3_IMPLEMENTATION

#include "MP4Encoder.hpp"
#include <iostream>
#include "SDL3/SDL.h"
#include <math.h>
#include <algorithm>
#include <stdio.h>
#include <vector>
#include "FFTProcessor.hpp"

using namespace std;

class Exporter
{
public:
	Exporter(SDL_Window* win, SDL_Renderer* ren, int width, int height, int FPS);
	~Exporter();

	streamData* specData = (streamData*)malloc(sizeof(streamData));
	SDL_AudioSpec audioSpec = SDL_AudioSpec();
	int frameCount = 0;
	int currentFrameSize = 0;
	uint16_t bands = 64;

	void init();
	void eventCall(SDL_Event* event);
	void update(double delta);
	void draw(double delta);
	void clean();
	void close();
	//void playSound();
	//void stopSound();
	bool loadMP3(string path);
	bool loadWAV(string path);

	bool running() { return isRunning; };

private:
	MP4Encoder* encoder;
	FFTProcessor* processor;
	SDL_Window* window;
	SDL_Renderer* renderer;
	vector<SDL_FRect> rectangles;
	SDL_AudioDeviceID audioDevice;
	double* bins;
	double* bandFreqs;
	double* spectrum;
	SDL_AudioStream* audioStream;
	void* audBuffer;
	int audSize;
	double sampleRatio;
	double sampleRate;
	int curFrame;
	void* audFrame;
	size_t totalSamples;
	int WINDOW_WIDTH;
	int WINDOW_HEIGHT;
	int FPS;
	SDL_Rect scrRect;
	double timeBase;
	size_t prevIdx;

	bool isRunning = true;
};

#endif