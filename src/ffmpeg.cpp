#include "ffmpeg.h"

bool FFMpegDecoder :: bInit = false;

#define SampleFloatFrom8BitInt 0.0078125f
#define SampleFloatFrom16BitInt 3.051757812500000e-05f
#define SampleFloatFrom32BitInt 4.656612873077393e-10f

FFMpegDecoder :: FFMpegDecoder(const char *filename)
  : frameBuffer(1024)
{
  fmt_ctx = NULL;
  audio_dec_ctx = NULL;
  audio_stream = NULL;
  audio_stream_idx = -1;
  frame = NULL;
  audio_frame_count = 0;
  error = 0;

  if(!bInit) av_register_all();
  bInit = true;

  if (avformat_open_input(&fmt_ctx, filename, NULL, NULL) < 0) {
    fprintf(stderr, "Could not open source file %s\n", filename);
    error = -1;
  }
  if (avformat_find_stream_info(fmt_ctx, NULL) < 0) {
    fprintf(stderr, "Could not find stream information\n");
    error = -2;
  }

  int ret;

  ret = av_find_best_stream(fmt_ctx, AVMEDIA_TYPE_AUDIO, -1, -1, NULL, 0);
  if(ret < 0) {
    fprintf(stderr, "Could not find %s stream in input file '%s'\n",
            av_get_media_type_string(AVMEDIA_TYPE_AUDIO), filename);
    error = ret;
    return;
  }

  /* find decoder for the stream */
  audio_stream_idx = ret;
  AVStream *st = fmt_ctx->streams[audio_stream_idx];
  AVCodecContext *dec_ctx = st->codec;
  AVCodec *dec = avcodec_find_decoder(dec_ctx->codec_id);
  if (!dec) {
    fprintf(stderr, "Failed to find %s codec\n",
            av_get_media_type_string(AVMEDIA_TYPE_AUDIO));
    error = AVERROR(EINVAL);
    return;
  }

  /* Init the decoders, with or without reference counting */
  AVDictionary *opts = NULL;
  if ((ret = avcodec_open2(dec_ctx, dec, &opts)) < 0) {
    fprintf(stderr, "Failed to open %s codec\n",
            av_get_media_type_string(AVMEDIA_TYPE_AUDIO));
    error = ret;
    return;
  }

 
  audio_stream = fmt_ctx->streams[audio_stream_idx];
  if(!audio_stream) {
    fprintf(stderr, "Could not find audio stream in the input, aborting\n");
    error = -3;
    return;
  }

  audio_dec_ctx = audio_stream->codec;
  Fs = audio_dec_ctx->sample_rate;
  channels = audio_dec_ctx->channels;
  AVRational base = audio_stream->time_base;
  duration = ((audio_stream->duration * Fs * base.num) / (base.den));

  frame = av_frame_alloc();
  if(!frame) {
    fprintf(stderr, "Could not allocate frame\n");
    error = AVERROR(ENOMEM);
    return;
  }

  /* initialize packet, set data to NULL, let the demuxer fill it */
  av_init_packet(&pkt);
  pkt.data = NULL;
  pkt.size = 0;
}

FFMpegDecoder :: ~FFMpegDecoder()
{
  close();
}

bool FFMpegDecoder :: isError()
{
  return (error != 0);
}

int FFMpegDecoder :: decode_packet(int *got_frame)
{
  int ret = 0;
  int decoded = pkt.size;
  *got_frame = 0;
  if (pkt.stream_index == audio_stream_idx) {
    ret = avcodec_decode_audio4(audio_dec_ctx, frame, got_frame, &pkt);
    if(ret < 0) {
      fprintf(stderr, "Error decoding audio frame (%d)\n",ret);
      error = ret;
      return error;
    }
    decoded = FFMIN(ret, pkt.size);
    if (*got_frame) {
      enum AVSampleFormat fmt = audio_dec_ctx->sample_fmt;
      int totSamples = frame->nb_samples * channels;

      frameBuffer.grow(totSamples);
      float *buf = frameBuffer.getWriteBuf();
      if(av_sample_fmt_is_planar(fmt)) {
        int c,k;
        uint8_t *in;
        switch(fmt) {
        case AV_SAMPLE_FMT_U8P:
          for(c=0; c<channels; c++) {
            in = frame->extended_data[c];
            for(k=c; k<totSamples; k+=channels) {
              buf[k] = (*(uint8_t *)in - 0x80) * SampleFloatFrom8BitInt;
              in++;
            }
          }
          break;
        case AV_SAMPLE_FMT_S16P:
          for(c=0; c<channels; c++) {
            in = frame->extended_data[c];
            for(k=c; k<totSamples; k+=channels) {
              buf[k] = (*(uint16_t *)in) * SampleFloatFrom16BitInt;
              in+=2;
            }
          }
          break;
        case AV_SAMPLE_FMT_S32P:
          for(c=0; c<channels; c++) {
            in = frame->extended_data[c];
            for(k=c; k<totSamples; k+=channels) {
              buf[k] = (*(uint32_t *)in) * SampleFloatFrom32BitInt;
              in+=4;
            }
          }
          break;
        case AV_SAMPLE_FMT_FLTP:
          for(c=0; c<channels; c++) {
            in = frame->extended_data[c];
            for(k=c; k<totSamples; k+=channels) {
              buf[k] = (*(float*)in);
              in+=4;
            }
          }
          break;
        case AV_SAMPLE_FMT_DBLP:
          for(c=0; c<channels; c++) {
            in = frame->extended_data[c];
            for(k=c; k<totSamples; k+=channels) {
              buf[k] = (float)(*(double*)in);
              in+=8;
            }
          }
          break;
        }
      } else {
        int c,k;
        uint8_t *in = frame->extended_data[0];
        switch(fmt) {
        case AV_SAMPLE_FMT_U8:
          for(k=0; k<totSamples; k++) {
            buf[k] = (*(uint8_t *)in - 0x80) * SampleFloatFrom8BitInt;
            in++;
          }
          break;
        case AV_SAMPLE_FMT_S16:
          for(k=0; k<totSamples; k++) {
            buf[k] = (*(uint16_t *)in) * SampleFloatFrom16BitInt;
            in+=2;
          }
          break;
        case AV_SAMPLE_FMT_S32:
          for(k=0; k<totSamples; k++) {
            buf[k] = (*(uint32_t *)in) * SampleFloatFrom32BitInt;
            in+=4;
          }
          break;
        case AV_SAMPLE_FMT_FLT:
          for(k=0; k<totSamples; k++) {
            buf[k] = (*(float*)in);
            in+=4;
          }
          break;
        case AV_SAMPLE_FMT_DBL:
          for(k=0; k<totSamples; k++) {
            buf[k] = (float)(*(double*)in);
            in+=8;
          }
          break;
        }
      }
      frameBuffer.writePos += totSamples;
    }
  }
  return decoded;
}

long FFMpegDecoder :: read(float *buf, long n)
{
  int got_frame;
  long nRead = 0;
  while(n > 0 && av_read_frame(fmt_ctx, &pkt) >= 0) {
    AVPacket orig_pkt = pkt;
    do {
      int ret = decode_packet(&got_frame);
      if(ret < 0)
        break;
      pkt.data += ret;
      pkt.size -= ret;
    } while (pkt.size > 0);
    av_free_packet(&orig_pkt);
    long nToRead = min(n,frameBuffer.nReadable());
    frameBuffer.read(buf,nToRead);
    n -= nToRead/channels;
    nRead += nToRead/channels;
  }
  pkt.data = NULL;
  pkt.size = 0;
  do {
    decode_packet(&got_frame);
    long nToRead = min(n,frameBuffer.nReadable());
    frameBuffer.read(buf,nToRead);
    n -= nToRead/channels;
    nRead += nToRead/channels;
  } while (got_frame && n > 0);

  return nRead;
}

bool FFMpegDecoder :: close()
{
  avcodec_close(audio_dec_ctx);
  avformat_close_input(&fmt_ctx);
  av_frame_free(&frame);
}


int FFMpegDecoder :: getSampleRate()
{
  return Fs;
}

int FFMpegDecoder :: getChannels()
{
  return channels;
}

long long int FFMpegDecoder :: getFrames()
{
  return duration;
}
