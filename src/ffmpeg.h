#ifndef FFMPEG_H
#define FFMPEG_H

#include "import.h"
#include "buffer.h"

using namespace _sbsms_;

extern "C" {
#include <libavutil/imgutils.h>
#include <libavutil/samplefmt.h>
#include <libavutil/timestamp.h>
#include <libavformat/avformat.h>
}

class FFMpegDecoder 
  : public AudioDecoder
{
public:
  FFMpegDecoder(const char *filename);
  ~FFMpegDecoder();
  long read(float *buf, long block_size);
  int getChannels();
  int getSampleRate();
  long long int getFrames();
  bool isError();
  bool close();

  static int open_codec_context(int *stream_idx,
                                AVFormatContext *fmt_ctx,
                                enum AVMediaType type);
 protected:

  int decode_packet(int *got_frame);
  ArrayRingBuffer<float> frameBuffer;
  AVFormatContext *fmt_ctx;
  AVCodecContext *audio_dec_ctx;
  AVStream *audio_stream;
  int audio_stream_idx;
  AVFrame *frame;
  AVPacket pkt;
  int audio_frame_count;
  int error;
  static bool bInit;
  int channels;
  int Fs;
  long long int duration;
};

#endif
