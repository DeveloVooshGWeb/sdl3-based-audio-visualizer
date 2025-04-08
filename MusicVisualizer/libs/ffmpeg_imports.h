#pragma once

#define _CRTDBG_MAP_ALLOC  
#include <stdlib.h>  
#include <crtdbg.h> 

extern "C" {
	#include <libavcodec/avcodec.h>
	#include <libavutil/opt.h>
	#include <libavutil/imgutils.h>
	#include <libswscale/swscale.h>
	#include <libavformat/avformat.h>
	#include <libavutil/timestamp.h>
	#include <libavutil/mathematics.h>
	#include <libavutil/channel_layout.h>
	#include <libavutil/avassert.h>
	#include <libavfilter/avfilter.h>
	#include <libavfilter/buffersink.h>
	#include <libavfilter/buffersrc.h>
	#include <libswresample/swresample.h>
	#include <libavutil/audio_fifo.h>
}

#ifdef av_err2str
#undef av_err2str
#include <string>
	av_always_inline std::string av_err2string(int errnum) {
		char str[AV_ERROR_MAX_STRING_SIZE];
		return av_make_error_string(str, AV_ERROR_MAX_STRING_SIZE, errnum);
	}
#define av_err2str(err) av_err2string(err).c_str()
#endif  // av_err2str

#ifdef av_ts2timestr
#undef av_ts2timestr
#include <string>
	av_always_inline std::string av_ts2timestring(int64_t ts, AVRational* tb) {
		char str[AV_TS_MAX_STRING_SIZE];
		return av_ts_make_time_string(str, ts, tb);
	}
#define av_ts2timestr(ts, tb) av_ts2timestring((ts), (tb)).c_str()
#endif  // av_ts2timestr

#ifdef av_ts2str
#undef av_ts2str
#include <string>
	av_always_inline std::string av_ts2str(int64_t a) {
		char str[AV_TS_MAX_STRING_SIZE];
		return av_ts_make_string(str, a);
	}
#define av_ts2str(a) av_ts2str((a)).c_str()
#endif  // av_ts2str

#ifdef AV_CHANNEL_LAYOUT_STEREO
#undef AV_CHANNEL_LAYOUT_STEREO
#define AV_CHANNEL_LAYOUT_STEREO { AV_CHANNEL_ORDER_NATIVE, (2), { ((1ULL << AV_CHAN_FRONT_LEFT) | (1ULL << AV_CHAN_FRONT_RIGHT)) }, 0 }
#endif

// a wrapper around a single output AVStream
typedef struct OutputStream {
	AVStream* st;
	AVCodecContext* enc;

	/* pts of the next frame that will be generated */
	int64_t next_pts;
	int samples_count;

	AVFrame* frame;
	AVFrame* tmp_frame;

	AVPacket* tmp_pkt;

	struct SwsContext* sws_ctx;
	struct SwrContext* swr_ctx;
} OutputStream;