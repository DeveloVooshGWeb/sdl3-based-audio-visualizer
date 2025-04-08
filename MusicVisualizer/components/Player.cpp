#include "Player.hpp"

#include "../libs/minimp3_ex.h"

streamData* Player::specData = (streamData*)malloc(sizeof(streamData));
int Player::len = 0;
SDL_AudioSpec Player::audioSpec = SDL_AudioSpec();
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
			actualSize = bytesElapsed - (bytesElapsed - FFT_SIZE);
			bytesElapsed = FFT_SIZE;
		}
		//cout << actualSize << endl;
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
	specData->p = fftw_plan_dft_r2c_1d(FFT_SIZE, specData->in, specData->out, FFTW_ESTIMATE);
	// Load MP3 to buffer;
	if (!loadMP3("springcmastered.mp3")) return;
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
	sampleRate = audioSpec.freq;
	processor = new FFTProcessor(sampleRate, FFT_SIZE, bands, FREQ_START, FREQ_END);
	bins = (double*)malloc(processor->getBinSize() * sizeof(double));
	bandFreqs = (double*)malloc((bands+1) * sizeof(double));
	spectrum = (double*)malloc(bands * sizeof(double));
	// Iterate frequency bins
	//double logMin = log10(FREQ_START);
	//double logMax = log10(FREQ_END);
	double lastBin = processor->getBinSize() - 1;
	double freqStart = processor->clamp(FREQ_START * FFT_SIZE / sampleRate, 1.0, lastBin);
	double freqEnd = processor->clamp(FREQ_END * FFT_SIZE / sampleRate, 1.0, lastBin);
	for (int i = 0; i < bands; i++)
	{
		rectangles.push_back(SDL_FRect());
		bandFreqs[i] = processor->clamp(processor->logint(FREQ_START, FREQ_END, (double)i / bands), FREQ_START, FREQ_END);
		
		double t = (double)i / bands;
		// Logarithmic spacing
		//bandFreqs[i] = pow(10, logMin + t * (logMax - logMin));
		//cout << bandFreqs[i] << endl;
	}
	bandFreqs[bands] = FREQ_END;
	processor->assign(bandFreqs, 0.5, 0.55);
}

void Player::playSound()
{
	SDL_ResumeAudioStreamDevice(audioStream);
}

void Player::stopSound()
{
	SDL_PauseAudioStreamDevice(audioStream);
}

void Player::eventCall(SDL_Event* event)
{
	switch (event->type) {
		case SDL_EVENT_QUIT:
			SDL_Log("SDL3 event quit");
			isRunning = false;
			break;
	}
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
		// Calculate magnitudes
		for (size_t i = 0; i < processor->getBinSize(); i++)
		{
			double real = specData->out[i][0];
			double imag = specData->out[i][1];
			double magnitude = sqrt(real * real + imag * imag);
			bins[i] = magnitude;
			//bins[i] = processor->toDb(magnitude);
		}
		processor->pipe(bins, spectrum);
		bytesElapsed = 0;
	}
	// Render visualizer bars
	for (int i = 0; i < rectangles.size(); i++)
	{
		rectangles[i].x = (1920.0 / 2.0 - (BAR_WIDTH * BAR_SPACE * (rectangles.size() / 2.0))) + BAR_WIDTH * BAR_SPACE * i;
		rectangles[i].y = 960;
		rectangles[i].w = BAR_WIDTH;
		// Create a multiplier magnitude value
		double f = 10.0 * log10(spectrum[i]*(i*0.01+1)) / 24.0; //24.0; //i * ((log10(FREQ_END) - log10(FREQ_START)) / (buckets))) / 24.0; //(i + 1))) / 32.0; //(10.0 * log10(maxFreqs[i] * (buckets-i))) / 24.0;
		//cout << f << endl;
		float h = f * -180.0;
		h = isnan(h) ? -4.0 : h;
		rectangles[i].h = h > -4.0 ? -4.0 : h; //lerp(rectangles[i].h, h > -4.0 ? -4.0 : h, delta * 12.0);
		SDL_RenderFillRect(renderer, &rectangles[i]);
		if (i > 0)
		{
			//SDL_RenderLine(renderer, rectangles[i - 1].x, rectangles[i - 1].y + rectangles[i - 1].h, rectangles[i].x, rectangles[i].y + rectangles[i].h);
		}
		//SDL_RenderPoint(renderer, rectangles[i].x, rectangles[i].y + rectangles[i].h);
	}
	for (int i = 0; i < FFT_SIZE/2.0; i++)
	{
		double pain = processor->toDb(bins[i]);
		//SDL_RenderPoint(renderer, 1920.0/2.0 - (i-(FFT_SIZE / 2.0/2.0))*1.0, pain+256.0);
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
	//free(freqBin);
	free(processor);
	free(bandFreqs);
	free(spectrum);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
}

void Player::close()
{
	isRunning = false;
}