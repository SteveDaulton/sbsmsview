#ifndef CONVERT_H
#define CONVERT_H

#include "sbsms.h"
#include "import.h"

using namespace _sbsms_;


long audio_convert_from(float *to, long off1, audio *from, long off2, long n, int channels = 2);
long audio_convert_to(audio *to, long off1, float *from, long off2, long n, int channels = 2);

class SBSMSInterfaceDecoderImp;

class SBSMSInterfaceDecoderImp;
class SBSMSInterfaceDecoder : public SBSMSInterfaceSliding {
public:
  SBSMSInterfaceDecoder(Slide *rateSlide, Slide *pitchSlide, bool bPitchInputReference,
                        int channels, const SampleCountType &samples, long preSamples,
                        SBSMSQuality *quality, AudioDecoder *decoder, float pitch);
  virtual ~SBSMSInterfaceDecoder();
  long samples(audio *buf, long n);
  
protected:
  SBSMSInterfaceDecoderImp *imp;
};


typedef bool (*progress_cb)(float progress, const char *msg, void *data);

bool sbsms_convert(const char *filenameIn, const char *filenameOut, bool bAnalyze, bool bSynthesize, progress_cb progressCB, void *data, float rate0, float rate1, float pitch0, float pitch1, float volume);
SBSMS *createCachedSBSMS(const char *filenameIn);

#endif 
