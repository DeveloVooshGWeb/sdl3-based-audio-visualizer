#include "MP4Encoder.hpp"

/*
MP4Encoder::MP4Encoder(int width, int height, int framerate, int bitrate, AVPixelFormat pix_fmt, string fpath, char** flags, size_t flagLen)
{
	//_data = MP4Data();
	_data.width = width;
	_data.height = height;
	_data.framerate = framerate;
	_data.bitrate = bitrate;
	_data.pix_fmt = pix_fmt;
	_data.fpath = fpath;
	_data.flags = flags;
	_data.flagLen = flagLen;
	*/
MP4Encoder::MP4Encoder(MP4Data data)
{
	_data = data;
	for (size_t i = 0; i + 1 < _data.flag_len; i += 2)
	{
		av_dict_set(&_opt, _data.flags[i].c_str(), _data.flags[i + 1].c_str(), 0);
	}
	avformat_alloc_output_context2(&_oc, NULL, NULL, _data.fpath.c_str());
	if (!_oc) {
		cout << "Could not deduce output format from file extension: using MPEG." << endl;
		avformat_alloc_output_context2(&_oc, NULL, "mpeg", _data.fpath.c_str());
	}
	if (!_oc)
	{
		cout << "Failed!";
		return;
	}
	_fmt = _oc->oformat;
	if (_fmt->video_codec != AV_CODEC_ID_NONE)
	{
		_add_stream(&_video_st, &_vcodec, _fmt->video_codec);
	}
	else
	{
		return;
	}
	if (_fmt->audio_codec != AV_CODEC_ID_NONE)
	{
		_add_stream(&_audio_st, &_acodec, _fmt->audio_codec);
	}
	else
	{
		return;
	}
	_open_video(_vcodec, &_video_st);
	_open_audio(_acodec, &_audio_st);
	int r;
	r = avio_open(&_oc->pb, _data.fpath.c_str(), AVIO_FLAG_WRITE);
	if (r < 0)
	{
		cout << "Could not open: " << _data.fpath << ", " << av_err2str(r) << endl;
		return;
	}
	r = avformat_write_header(_oc, &_opt);
	if (r < 0)
	{
		cout << "Error occured when opening output file: " << _data.fpath << ", " << av_err2str(r) << endl;
		return;
	}
	_img_pixels = _data.width * _data.height;
	_working = true;
}

MP4Encoder::~MP4Encoder()
{
}

void MP4Encoder::_log_packet(const AVFormatContext* fmt_ctx, const AVPacket* pkt)
{
	AVRational* time_base = &fmt_ctx->streams[pkt->stream_index]->time_base;

	printf("pts:%s pts_time:%s dts:%s dts_time:%s duration:%s duration_time:%s stream_index:%d\n",
		av_ts2str(pkt->pts), av_ts2timestr(pkt->pts, time_base),
		av_ts2str(pkt->dts), av_ts2timestr(pkt->dts, time_base),
		av_ts2str(pkt->duration), av_ts2timestr(pkt->duration, time_base),
		pkt->stream_index);
}

bool MP4Encoder::isWorking()
{
	return _working;
}

void MP4Encoder::_add_stream(OutputStream* ost, const AVCodec** codec, enum AVCodecID codec_id)
{
	AVCodecContext* c;
	//int i;
	*codec = avcodec_find_encoder(codec_id);
	if (!(*codec))
	{
		cout << "Could not find encoder for " << avcodec_get_name(codec_id) << endl;
		_working = false;
		return;
	}

	ost->tmp_pkt = av_packet_alloc();
	if (!ost->tmp_pkt)
	{
		cout << "Could not allocate AVPacket" << endl;
		_working = false;
		return;
	}

	ost->st = avformat_new_stream(_oc, NULL);
	if (!ost->st)
	{
		cout << "Could not allocate stream" << endl;
		return;
	}
	ost->st->id = _oc->nb_streams - 1;
	c = avcodec_alloc_context3(*codec);
	if (!c)
	{
		cout << "Could not allocate an encoding context" << endl;
		return;
	}
	ost->enc = c;
	AVChannelLayout stereo = AV_CHANNEL_LAYOUT_STEREO;
	switch ((*codec)->type)
	{
		case AVMEDIA_TYPE_AUDIO:
			c->sample_fmt = _data.smp_fmt; //(*codec)->sample_fmts ? (*codec)->sample_fmts[0] : AV_SAMPLE_FMT_FLTP;
			c->bit_rate = _data.abitrate;
			c->sample_rate = _data.samplerate;
			/*
			if ((*codec)->supported_samplerates)
			{
				//avcodec_get_supported_config()
				c->sample_rate = (*codec)->supported_samplerates[0];
				for (i = 0; (*codec)->supported_samplerates[i]; i++)
				{
					if ((*codec)->supported_samplerates[i] == 44100)
					{
						c->sample_rate = 44100;
					}
				}
			}
			*/
			av_channel_layout_copy(&c->ch_layout, &stereo);
			//c->ch_layout = AV_CHANNEL_LAYOUT_STEREO;
			ost->st->time_base = AVRational{ 1, c->sample_rate };
			break;

		case AVMEDIA_TYPE_VIDEO:
			c->codec_id = codec_id;
			c->bit_rate = _data.vbitrate;
			c->width = _data.width;
			c->height = _data.height;
			ost->st->time_base = AVRational{ 1, _data.framerate };
			c->time_base = ost->st->time_base;

			c->gop_size = _data.gop_size;
			c->pix_fmt = _data.pix_fmt; //AV_PIX_FMT_YUVA420P;
			if (c->codec_id == AV_CODEC_ID_MPEG2VIDEO)
			{
				c->max_b_frames = 2;
			}
			if (c->codec_id == AV_CODEC_ID_MPEG1VIDEO)
			{
				c->mb_decision = 2;
			}
			break;

		default:
			break;
	}
}

int MP4Encoder::_write_frame(AVCodecContext* c, AVStream* st, AVFrame* frame, AVPacket* pkt)
{
	int r;
	r = avcodec_send_frame(c, frame);
	if (r < 0)
	{
		cout << "Error sending a frame to the encoder!, " << r << endl;
		_working = false;
		return 1;
	}
	while (r >= 0)
	{
		r = avcodec_receive_packet(c, pkt);
		if (r == AVERROR(EAGAIN) || r == AVERROR_EOF)
		{
			break;
		}
		else if (r < 0)
		{
			cout << "Error encoding a frame!, " << av_err2str(r) << endl;
			_working = false;
			return 1;
		}
		av_packet_rescale_ts(pkt, c->time_base, st->time_base);
		pkt->stream_index = st->index;

		//log_packet(_oc, pkt);
		r = av_interleaved_write_frame(_oc, pkt);
		if (r < 0)
		{
			cout << "Error while writing output packet!, " << av_err2str(r) << endl;
			_working = false;
			return 1;
		}
	}
	return r == AVERROR_EOF ? 1 : 0;
}

AVFrame* MP4Encoder::_alloc_audio_frame(enum AVSampleFormat sample_fmt, const AVChannelLayout* channel_layout, int sample_rate, int nb_samples)
{
	AVFrame* frame = av_frame_alloc();
	if (!frame)
	{
		cout << "Error allocating an audio frame" << endl;
		_working = false;
		return NULL;
	}
	frame->format = sample_fmt;
	av_channel_layout_copy(&frame->ch_layout, channel_layout);
	frame->sample_rate = sample_rate;
	frame->nb_samples = nb_samples;

	if (nb_samples)
	{
		if (av_frame_get_buffer(frame, 0) < 0)
		{
			cout << "Error allocating an audio buffer" << endl;
			_working = false;
			return NULL;
		}
	}

	return frame;
}

void MP4Encoder::_open_audio(const AVCodec* codec, OutputStream* ost)
{
	AVCodecContext* c;
	int nb_samples;
	int r;
	AVDictionary* opt = NULL;

	c = ost->enc;

	av_dict_copy(&opt, _opt, 0);
	r = avcodec_open2(c, codec, &opt);
	av_dict_free(&opt);
	if (r < 0)
	{
		cout << "Could not open audio codec: " << av_err2str(r) << endl;
		_working = false;
		return;
	}

	if (c->codec->capabilities && AV_CODEC_CAP_VARIABLE_FRAME_SIZE)
	{
		nb_samples = 10000;
	}
	else
	{
		nb_samples = c->frame_size;
	}

	ost->frame = _alloc_audio_frame(c->sample_fmt, &c->ch_layout, c->sample_rate, nb_samples);
	
	r = avcodec_parameters_from_context(ost->st->codecpar, c);
	if (r < 0)
	{
		cout << "Could not copy the stream parameters" << endl;
		_working = false;
		return;
	}

	ost->swr_ctx = swr_alloc();
	if (!ost->swr_ctx)
	{
		cout << "Could not allocate resampler context" << endl;
		_working = false;
		return;
	}

	av_opt_set_chlayout(ost->swr_ctx, "in_chlayout", &c->ch_layout, 0);
	av_opt_set_int(ost->swr_ctx, "in_sample_rate", c->sample_rate, 0);
	av_opt_set_sample_fmt(ost->swr_ctx, "in_sample_fmt", AV_SAMPLE_FMT_S16, 0); //AV_SAMPLE_FMT_S16, 0);
	av_opt_set_chlayout(ost->swr_ctx, "out_chlayout", &c->ch_layout, 0);
	av_opt_set_int(ost->swr_ctx, "out_sample_rate", c->sample_rate, 0);
	av_opt_set_sample_fmt(ost->swr_ctx, "out_sample_fmt", c->sample_fmt, 0);

	if ((r = swr_init(ost->swr_ctx)) < 0)
	{
		cout << "Failed to initialize the resampling context" << endl;
		_working = false;
		return;
	}
}

void MP4Encoder::_open_video(const AVCodec* codec, OutputStream* ost)
{
	int r;
	AVCodecContext* c = ost->enc;
	AVDictionary* opt = NULL;

	av_dict_copy(&opt, _opt, 0);

	r = avcodec_open2(c, codec, &opt);
	av_dict_free(&opt);
	if (r < 0)
	{
		cout << "Could not open video codec: " << av_err2str(r) << endl;
		_working = false;
		return;
	}

	ost->frame = _alloc_frame(c->pix_fmt, c->width, c->height);
	if (!ost->frame)
	{
		cout << "Could not allocate video frame" << endl;
		_working = false;
		return;
	}

	ost->tmp_frame = NULL;
	if (c->pix_fmt != AV_PIX_FMT_YUV420P)
	{
		ost->tmp_frame = _alloc_frame(AV_PIX_FMT_YUV420P, c->width, c->height);
		if (!ost->tmp_frame)
		{
			cout << "Could not allocate temporary video frame" << endl;
			_working = false;
			return;
		}
	}

	r = avcodec_parameters_from_context(ost->st->codecpar, c);
	if (r < 0)
	{
		cout << "Could not copy the stream parameters" << endl;
		_working = false;
		return;
	}
}

void MP4Encoder::init_video_write(AVPixelFormat pfmt)
{
	OutputStream* ost = &_video_st;
	AVCodecContext* c = ost->enc;

	if (!ost->sws_ctx)
	{
		ost->sws_ctx = sws_getContext(c->width, c->height, pfmt, c->width, c->height, c->pix_fmt, SCALE_FLAGS, NULL, NULL, NULL);
		if (!ost->sws_ctx)
		{
			cout << "Could not initialize the conversion context" << endl;
			_working = false;
			return;
		}
	}
}

int MP4Encoder::write_image_matrix(uint8_t** data, uint8_t channels)
{
	OutputStream* ost = &_video_st;
	AVCodecContext* c = ost->enc;

	// (1.0 - stream duration)
	if (av_compare_ts(ost->next_pts, c->time_base, 1.0, AVRational{ 1, 1 }) > 0)
	{
		return 1;
	}

	if (av_frame_make_writable(ost->frame) < 0)
	{
		_working = false;
		return 1;
	}

	for (uint8_t i = 0; i < channels; i++)
	{
		for (uint8_t j = 0; j < _img_pixels; j++)
		{
			ost->tmp_frame->data[i][j] = data[i][j];
		}
	}

	sws_scale(ost->sws_ctx, (const uint8_t* const*)ost->tmp_frame->data, ost->tmp_frame->linesize, 0, c->height, ost->frame->data, ost->frame->linesize);

	ost->frame->pts = ost->next_pts++;
	
	return _write_frame(c, ost->st, ost->frame, ost->tmp_pkt);
}

void MP4Encoder::init_audio_write(AVSampleFormat sfmt)
{
	OutputStream* ost = &_audio_st;
	AVCodecContext* c = ost->enc;
	ost->tmp_frame = _alloc_audio_frame(sfmt, &c->ch_layout, c->sample_rate, ost->frame->nb_samples);
}

int MP4Encoder::write_audio_samples(uint8_t* data, size_t sz, int64_t samples)
{
	OutputStream* ost = &_audio_st;
	AVCodecContext* c = ost->enc;
	int r;
	int dst_nb_samples;
	AVFrame* frame = ost->tmp_frame;
	uint8_t* q = (uint8_t*)frame->data[0];
	for (size_t i = 0; i < sz; i++)
	{
		q[i] = data[i];
	}
	frame->pts = ost->next_pts;
	ost->next_pts += samples;

	if (frame)
	{
		dst_nb_samples = swr_get_delay(ost->swr_ctx, c->sample_rate) + samples;
		av_assert0(dst_nb_samples == samples);
		r = av_frame_make_writable(ost->frame);
		if (r < 0)
		{
			cout << "Could not make audio frame writable" << endl;
			_working = false;
			return 1;
		}

		r = swr_convert(ost->swr_ctx, ost->frame->data, dst_nb_samples, (const uint8_t**)frame->data, samples);
		if (r < 0)
		{
			cout << "Error while converting" << endl;
			_working = false;
			return 1;
		}

		frame = ost->frame;

		frame->pts = av_rescale_q(ost->samples_count, AVRational{ 1, c->sample_rate }, c->time_base);
		ost->samples_count += dst_nb_samples;
	}

	return _write_frame(c, ost->st, frame, ost->tmp_pkt);
}

void MP4Encoder::_close_stream(OutputStream* ost)
{
	avcodec_free_context(&ost->enc);
	av_frame_free(&ost->frame);
	if (ost->tmp_frame)	av_frame_free(&ost->tmp_frame);
	av_packet_free(&ost->tmp_pkt);
	if (ost->sws_ctx) sws_freeContext(ost->sws_ctx);
	swr_free(&ost->swr_ctx);
}

void MP4Encoder::finalize()
{
	_finalize();
}

void MP4Encoder::_finalize()
{
	av_write_trailer(_oc);
	_close_stream(&_video_st);
	_close_stream(&_audio_st);
	avio_closep(&_oc->pb);
	avformat_free_context(_oc);
	MP4Encoder::~MP4Encoder();
}

AVFrame* MP4Encoder::_alloc_frame(enum AVPixelFormat pix_fmt, int width, int height)
{
	AVFrame* frame;
	int r;

	frame = av_frame_alloc();
	if (!frame)
	{
		cout << "Frame allocation failed" << endl;
		_working = false;
		return NULL;
	}
	
	frame->format = pix_fmt;
	frame->width = width;
	frame->height = height;

	r = av_frame_get_buffer(frame, 0);
	if (r < 0)
	{
		cout << "Could not allocate frame data." << endl;
		_working = false;
		return NULL;
	}

	return frame;
}