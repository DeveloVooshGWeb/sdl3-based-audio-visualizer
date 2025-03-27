#include "Player.hpp"

#include "minimp3_ex.h"

#define FREQ_START 34
#define FREQ_END 19044
#define FFT_SIZE 5120

streamData* Player::specData = (streamData*)malloc(sizeof(streamData));
int Player::len = 0;
SDL_AudioSpec Player::audioSpec = SDL_AudioSpec();
float* Player::audioFrames = (float*)malloc(sizeof(float) * FFT_SIZE * 8);
float* Player::cBuf = (float*)malloc(sizeof(float) * 4096);
int Player::currentFrameSize = 0;
int Player::frameCount = 0;
bool Player::mustCall = false;
float Player::toLog = 0.0;
int Player::bytesElapsed = 0;


Player::Player(SDL_Window* win, SDL_Renderer* ren)
{
	window = win;
	renderer = ren;
}

Player::~Player()
{
}

void SDLCALL Player::onAudioData(void* userdata, const SDL_AudioSpec* spec, float* buffer, int buflen)
{
	if (bytesElapsed < FFT_SIZE)
	{
		int prevBytesElapsed = 0 + bytesElapsed;
		int actualSize = buflen / 2 / sizeof(float);
		bytesElapsed += actualSize;
		if (bytesElapsed > FFT_SIZE)
		{
			actualSize = (bytesElapsed - FFT_SIZE);
			bytesElapsed = FFT_SIZE;
		}
		for (int i = 0; i < actualSize; i++)
		{
			// Hann Function
			double multiplier = 0.5 * (1 - cos(2 * M_PI * (i + prevBytesElapsed) / (FFT_SIZE - 1)));
			specData->in[i + prevBytesElapsed] = multiplier * ((buffer[i * 2] + buffer[i * 2 + 1]) / 2.0);
		}
	}
}

bool Player::loadMP3(string path)
{
	mp3dec_ex_t dec;
	if (mp3dec_ex_open(&dec, path.c_str(), MP3D_SEEK_TO_SAMPLE))
	{
		SDL_Log("Sound failed to load: %s", SDL_GetError());
		isRunning = false;
		return false;
	}
	audioSpec.channels = dec.info.channels;
	audioSpec.freq = dec.info.hz;
	audioSpec.format = SDL_AUDIO_S16LE;
	audSize = dec.samples * sizeof(int16_t);
	audBuffer = malloc(audSize);
	size_t bytes_read = mp3dec_ex_read(&dec, (int16_t*)audBuffer, dec.samples);
	if (bytes_read != dec.samples)
	{
		SDL_Log("Warning in decoding.");
		if (dec.last_error)
		{
			SDL_Log("Sound failed to decode: %s", SDL_GetError());
			isRunning = false;
			return false;
		}
	}
	return true;
}

bool Player::loadWAV(string path)
{
	Uint32 len = 0;
	if (!SDL_LoadWAV(path.c_str(), &audioSpec, (Uint8**)&audBuffer, &len)) {
		SDL_Log("Sound failed to load: %s", SDL_GetError());
		isRunning = false;
		return false;
	}
	audSize = len;
	return true;
}

void Player::init()
{
	specData->in = (double*)malloc(sizeof(double) * FFT_SIZE);
	specData->out = (fftw_complex*)malloc(sizeof(fftw_complex) * FFT_SIZE); //(double*)malloc(sizeof(double) * FFT_SIZE);
	specData->magnitude = (double*)malloc(sizeof(double) * FFT_SIZE);
	if (specData->in == NULL || specData->out == NULL)
	{
		SDL_Log("Could not allocate spectrum data.");
		isRunning = false;
		return;
	}
	specData->p = fftw_plan_dft_r2c_1d(FFT_SIZE, specData->in, specData->out, FFTW_MEASURE);
	// Load MP3 to buffer;
	if (!loadMP3("samplestostudy.mp3")) return;
	//if (!loadWAV("sunsetfmastered.wav")) return;
	// Load WAV to buffer
	sampleRatio = FFT_SIZE / (double)audioSpec.freq;
	specData->startIdx = ceil(sampleRatio * FREQ_START);
	specData->specSize = min(
		ceil(sampleRatio * FREQ_END),
		FFT_SIZE / 2.0
	) - specData->startIdx;
	// Create audio stream
	audioStream = SDL_OpenAudioDeviceStream(SDL_AUDIO_DEVICE_DEFAULT_PLAYBACK, &audioSpec, NULL, NULL);
	if (!SDL_PutAudioStreamData(audioStream, audBuffer, audSize))
	{
		SDL_Log("Audio stream failed to write: %s", SDL_GetError());
		isRunning = false;
		return;
	}
	if (!SDL_FlushAudioStream(audioStream))
	{
		SDL_Log("Audio stream failed to flush: %s", SDL_GetError());
		isRunning = false;
		return;
	}
	// Play audio
	playSound();
	if (!SDL_SetAudioPostmixCallback(SDL_GetAudioStreamDevice(audioStream), onAudioData, NULL))
	{
		SDL_Log("Audio stream failed to set callback: %s", SDL_GetError());
		isRunning = false;
		return;
	}
	// Create pools along with frequency bins for each bucket
	freqBin = (double*)malloc((buckets + 1) * sizeof(double));
	spectrum = (double*)malloc(buckets * sizeof(double));
	//normalSpectrum = (double*)malloc((FFT_SIZE / 2) * sizeof(double));
	//sizes = (int*)malloc(buckets * sizeof(int));
	binSpacing = (FREQ_END - FREQ_START) / buckets;
	// Iterate frequency bins
	double logMin = log10(FREQ_START);
	double logMax = log10(FREQ_END);
	for (int i = 0; i < buckets; i++)
	{
		rectangles.push_back(SDL_FRect());
		spectrum[i] = 0; //1.7E-308;
		double t = (double)i / buckets;
		// Logarithmic spacing
		freqBin[i] = pow(10, logMin + t * (logMax - logMin));
	}
	freqBin[buckets] = FREQ_END;
	
}

void Player::playSound()
{
	SDL_ResumeAudioStreamDevice(audioStream);
}

void Player::stopSound()
{
	SDL_PauseAudioStreamDevice(audioStream);
}

void Player::eventCall(SDL_Event event)
{
	switch (event.type) {
		case SDL_EVENT_QUIT:
			SDL_Log("SDL3 event quit");
			isRunning = false;
			break;
	}
}

float Player::lerp(float a, float b, float f)
{
	return a * (1.0 - f) + (b * f);
}

void Player::update(double delta)
{
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
	SDL_RenderClear(renderer);
	SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
	if (bytesElapsed == FFT_SIZE)
	{
		// Execute FFT
		fftw_execute(specData->p);
		string toPrint = "";
		// Prepare rendering
		
		// Get overall peak values
		//double peakmax = 1.7E-308;
		//int max_index = -1;
		for (int i = 0; i < buckets; i++)
		{
			spectrum[i] = 0;
			//sizes[i] = 0;
		}
		int minIdx = 0;
		int maxIdx = 0;
		float maxM = 0;
		int lastBucket = 0;
		// Calculate maximum magnitudes for each bucket
		for (int i = 0; i < FFT_SIZE / 2.0; i++)
		{
			double real = specData->out[i][0];
			double imag = specData->out[i][1];
			double magnitude = sqrt(real * real + imag * imag);

			double freq = i * (double)audioSpec.freq / FFT_SIZE;
			//cout << freq << endl;
			for (int j = 0; j < buckets; j++)
			{
				if (freq >= freqBin[j] && freq < freqBin[j + 1])
				{
					if (magnitude > spectrum[j])
					{
						spectrum[j] = magnitude;
						break;
					}
				}
			}
		}
		bytesElapsed = 0;
	}
	// Render visualizer bars
	for (int i = 0; i < rectangles.size(); i++)
	{
		rectangles[i].x = (1920.0 / 2.0 - (BAR_WIDTH * BAR_SPACE * (rectangles.size() / 2.0))) + BAR_WIDTH * BAR_SPACE * i;
		rectangles[i].y = 960;
		rectangles[i].w = BAR_WIDTH;
		// Create a multiplier magnitude value
		double f = (10.0 * log10(spectrum[i]) / 48.0); //i * ((log10(FREQ_END) - log10(FREQ_START)) / (buckets))) / 24.0; //(i + 1))) / 32.0; //(10.0 * log10(maxFreqs[i] * (buckets-i))) / 24.0;
		float h = f * -180.0;
		rectangles[i].h = lerp(rectangles[i].h, h > -4.0 ? -4.0 : h, delta * 12.0);
		SDL_RenderFillRect(renderer, &rectangles[i]);
	}
	// Render
	SDL_RenderPresent(renderer);
	
}

void Player::draw(double delta)
{
}

void Player::clean()
{
	bytesElapsed = 0;
	SDL_DestroyAudioStream(audioStream);
	fftw_destroy_plan(specData->p);
	//fftw_free(specData->in);
	//fftw_free(specData->out);
	free(specData);
	free(audBuffer);
	free(freqBin);
	free(spectrum);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
}

void Player::close()
{
	isRunning = false;
}