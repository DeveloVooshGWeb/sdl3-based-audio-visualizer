#ifndef MP4Encoder_hpp
#define MP4Encoder_hpp

#include "../libs/ffmpeg_imports.h"
#include <iostream>
#include <string>

#define SCALE_FLAGS SWS_BICUBIC

using namespace std;

typedef struct MP4Data {
	int width;
	int height;
	int framerate;
	int samplerate;
	int vbitrate;
	int abitrate;
	int gop_size;
	AVPixelFormat pix_fmt;
	AVSampleFormat smp_fmt;
	string fpath;
	string* flags;
	size_t flag_len;
	MP4Data() {};
	MP4Data(int w, int h, int fr, int sr, int vbr, int abr, int gsz, AVPixelFormat pfmt, AVSampleFormat sfmt, string fp, string* fl, size_t fllen) : width(w), height(h), framerate(fr), samplerate(sr), vbitrate(vbr), abitrate(abr), gop_size(gsz), pix_fmt(pfmt), smp_fmt(sfmt), fpath(fp), flags(fl), flag_len(fllen) {};
} MP4Data;

class MP4Encoder
{
public:
	MP4Encoder(MP4Data data);
	~MP4Encoder();
	
	//void log_packet(const AVFormatContext* fmt_ctx, const AVPacket* pkt);
	void init_video_write(AVPixelFormat pfmt);
	int write_image_matrix(uint8_t** data, uint8_t channels);
	void init_audio_write(AVSampleFormat sfmt);
	int write_audio_samples(uint8_t* data, size_t sz, int64_t samples);
	bool isWorking();
	void finalize();

private:
	bool _working = false;
	MP4Data _data;
	OutputStream _video_st;
	OutputStream _audio_st;
	const AVOutputFormat* _fmt;
	AVFormatContext* _oc;
	const AVCodec* _acodec;
	const AVCodec* _vcodec;
	AVDictionary* _opt = NULL;
	int _img_pixels;

	void _log_packet(const AVFormatContext* fmt_ctx, const AVPacket* pkt);

	void _add_stream(OutputStream* ost, const AVCodec** codec, enum AVCodecID codec_id);

	int _write_frame(AVCodecContext* c, AVStream* st, AVFrame* frame, AVPacket* pkt);
	AVFrame* _alloc_audio_frame(enum AVSampleFormat sample_fmt, const AVChannelLayout* channel_layout, int sample_rate, int nb_samples);
	void _open_audio(const AVCodec* codec, OutputStream* ost);

	AVFrame* _alloc_frame(enum AVPixelFormat pix_fmt, int width, int height);
	void _open_video(const AVCodec* codec, OutputStream* ost);
	void _close_stream(OutputStream* ost);

	void _finalize();
};

#endif
