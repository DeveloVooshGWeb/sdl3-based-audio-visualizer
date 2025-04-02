#include "SDL3/SDL.h"
//#include "SDL3/SDL_main.h"
#include "Player.hpp"

static SDL_Window* window = NULL;
static SDL_Renderer* renderer = NULL;

Uint64 ticks = 0;

#define WINDOW_WIDTH 1920
#define WINDOW_HEIGHT 1080
#define FRAMERATE 60

const float fps = 1000 / FRAMERATE;

using namespace std;

int main() {
	int result = SDL_Init(SDL_INIT_AUDIO | SDL_INIT_VIDEO);
	if (result < 0) {
		SDL_Log("SDL_Init error: %s", SDL_GetError());
		return -2;
	}
	if (!SDL_CreateWindowAndRenderer("Music Visualizer", WINDOW_WIDTH, WINDOW_HEIGHT, 0, &window, &renderer))
	{
		return -1;
	}
	Player player = Player(window, renderer);
	player.init();
	SDL_Event* event;
	event = new SDL_Event;
	while (player.running()) {
		while (SDL_PollEvent(event)) {
			player.eventCall(event);
		}
		double delta = (double)((SDL_GetTicks() - ticks)) / 1000;
		player.update(delta);
		player.draw(delta);
		ticks = SDL_GetTicks();
		SDL_Delay(fps);
	}
	player.clean();

	SDL_Quit();
	return 0;
}

