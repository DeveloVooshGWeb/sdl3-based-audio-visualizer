#include "Exporter.hpp"

#include "../libs/minimp3_ex.h"

Exporter::Exporter(SDL_Window* win, SDL_Renderer* ren, int width, int height, int framerate)
{
	window = win;
	renderer = ren;
	WINDOW_WIDTH = width;
	WINDOW_HEIGHT = height;
	FPS = framerate;
}

Exporter::~Exporter()
{
}

bool Exporter::loadMP3(string path)
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

bool Exporter::loadWAV(string path)
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


void Exporter::eventCall(SDL_Event* event)
{
	switch (event->type) {
	case SDL_EVENT_QUIT:
		SDL_Log("SDL3 event quit");
		isRunning = false;
		break;
	}
}

void Exporter::init()
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
	if (!loadMP3("springcmastered.mp3")) return;
	sampleRatio = FFT_SIZE / (double)audioSpec.freq;
	specData->startIdx = ceil(sampleRatio * FREQ_START);
	specData->specSize = min(
		ceil(sampleRatio * FREQ_END),
		FFT_SIZE / 2.0
	) - specData->startIdx;
	// Create pools along with frequency bins for each bucket
	sampleRate = audioSpec.freq;
	processor = new FFTProcessor(sampleRate, FFT_SIZE, bands, FREQ_START, FREQ_END);
	bins = (double*)malloc(processor->getBinSize() * sizeof(double));
	bandFreqs = (double*)malloc((bands + 1) * sizeof(double));
	spectrum = (double*)malloc(bands * sizeof(double));
	double lastBin = processor->getBinSize() - 1;
	double freqStart = processor->clamp(FREQ_START * FFT_SIZE / sampleRate, 1.0, lastBin);
	double freqEnd = processor->clamp(FREQ_END * FFT_SIZE / sampleRate, 1.0, lastBin);
	for (int i = 0; i < bands; i++)
	{
		rectangles.push_back(SDL_FRect());
		bandFreqs[i] = processor->clamp(processor->logint(FREQ_START, FREQ_END, (double)i / bands), FREQ_START, FREQ_END);

		double t = (double)i / bands;
	}
	bandFreqs[bands] = FREQ_END;
	processor->assign(bandFreqs, 0.5, 0.55);
	totalSamples = audSize / sizeof(int16_t) / 2.0;
	int fllen = 4;
	string* fl = new string[fllen]{ "-preset", "ultrafast", "-crf", "0"};
	encoder = new MP4Encoder(MP4Data(WINDOW_WIDTH, WINDOW_HEIGHT, FPS, (int)sampleRate, VBR, ABR, GOP_SIZE, AV_PIX_FMT_YUV420P, AV_SAMPLE_FMT_S16, string("exporttest.mp4"), fl, fllen));
	timeBase = 1.0 / FPS;
	scrRect = SDL_Rect();
	scrRect.w = WINDOW_WIDTH;
	scrRect.h = WINDOW_HEIGHT;
	encoder->init_audio_write(AV_SAMPLE_FMT_S16);
	encoder->init_video_write(AV_PIX_FMT_RGB24);
}

void Exporter::update(double delta)
{
	// RENDERING
	size_t idx = curFrame * timeBase * sampleRate;
	if (idx >= totalSamples)
	{
		// Finalize
		size_t diff = totalSamples - prevIdx;
		if (diff > 0)
		{
			encoder->write_audio_samples((uint8_t*)audBuffer + prevIdx, diff * sizeof(int16_t), diff);
		}
		encoder->finalize();
		cout << "Finished rendering!";
		isRunning = false;
		return;
	}
	cout << "Rendering frame: " << curFrame << endl;
	if (audioSpec.format != SDL_AUDIO_S16)
	{
		cout << "Not a signed int16." << endl;
		isRunning = false;
		return;
	}
	memcpy(audFrame, (char*)audBuffer + idx, FFT_SIZE * 2 * sizeof(int16_t));
	int16_t* audBufferI16 = static_cast<int16_t*>(audBuffer);
	double d1, d2;
	// Extract frames for FFT function
	for (size_t i = 0; i < FFT_SIZE; i++)
	{
		// Hann Function
		double multiplier = 0.5 * (1 - cos(2 * M_PI * i / (FFT_SIZE - 1)));
		d1 = audBufferI16[i * 2] / 32767.0;
		d2 = audBufferI16[i * 2 + 1] / 32767.0;
		specData->in[i] = multiplier * ((d1 + d2) / 2.0);
	}
	// Perform FFT
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
	SDL_RenderClear(renderer);
	SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
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
	// Render visualizer bars
	for (int i = 0; i < rectangles.size(); i++)
	{
		rectangles[i].x = (1920.0 / 2.0 - (BAR_WIDTH * BAR_SPACE * (rectangles.size() / 2.0))) + BAR_WIDTH * BAR_SPACE * i;
		rectangles[i].y = 960;
		rectangles[i].w = BAR_WIDTH;
		// Create a multiplier magnitude value
		double f = 10.0 * log10(spectrum[i] * (i * 0.01 + 1)) / 24.0; //24.0; //i * ((log10(FREQ_END) - log10(FREQ_START)) / (buckets))) / 24.0; //(i + 1))) / 32.0; //(10.0 * log10(maxFreqs[i] * (buckets-i))) / 24.0;
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
	for (int i = 0; i < FFT_SIZE / 2.0; i++)
	{
		double pain = processor->toDb(bins[i]);
		//SDL_RenderPoint(renderer, 1920.0/2.0 - (i-(FFT_SIZE / 2.0/2.0))*1.0, pain+256.0);
	}
	// Render (Export instead of presenting to screen)
	SDL_Surface* surface = SDL_RenderReadPixels(renderer, &scrRect);
	SDL_PixelFormat pixFmt = SDL_PIXELFORMAT_RGB24;
	if (surface && surface->format != pixFmt) {
		SDL_Surface* converted = SDL_ConvertSurface(surface, pixFmt);
		SDL_DestroySurface(surface);
		surface = converted;
	}
	// Render to FFMPEG
	if (surface) {
		encoder->write_image_matrix((uint8_t**)surface->pixels, 3);
		SDL_DestroySurface(surface);
	}
	// Write Audio Samples
	if (idx > 0)
	{
		encoder->write_audio_samples(((uint8_t*)audBuffer) + prevIdx, (idx - prevIdx) * sizeof(int16_t), idx - prevIdx);
	}
	prevIdx = idx;
	curFrame++;
}

void Exporter::draw(double delta)
{
}

void Exporter::clean()
{
	free(encoder);
	SDL_DestroyAudioStream(audioStream);
	fftw_destroy_plan(specData->p);
	free(specData);
	free(audBuffer);
	free(processor);
	free(bandFreqs);
	free(spectrum);
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
}

void Exporter::close()
{
	isRunning = false;
}