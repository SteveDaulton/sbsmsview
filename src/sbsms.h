// -*- mode: c++ -*-
#ifndef SBSMS_INCLUDE
#define SBSMS_INCLUDE

#include <stdio.h>
#include <vector>

namespace _sbsms_ {

typedef float t_fft[2];
typedef t_fft audio;
typedef long long int SampleCountType;
typedef long long int TimeType;
typedef unsigned int TrackIndexType;

enum {
  maxBands = 10,
  numQualityParams = 52
};


enum {
  NumGrainTypes = 3
};

struct SBSMSQualityParams {
  int bands;
  int H;
  int N[maxBands];
  int N0[maxBands];
  int N1[maxBands];
  int N2[maxBands];
  int res[maxBands];
};

class SBSMSQuality {
 public:
  SBSMSQuality(const SBSMSQualityParams *params);
  SBSMSQualityParams params;
  long getFrameSize();
  long getMaxPresamples();
  int getBandGrainsPerFrame(int band);
};

extern const SBSMSQualityParams SBSMSQualityStandard;

struct SBSMSFrame {
  float ratio0;
  float ratio1;
  audio *buf;
  long size;
};

typedef long (*SBSMSResampleCB)(void *cbData, SBSMSFrame *frame);

class SBSMSInterface {
 public:
  virtual ~SBSMSInterface() {}
  virtual long samples(audio *buf, long n) { return 0; }
  virtual float getStretch(float t)=0;
  virtual float getPitch(float t)=0;
  virtual long getPresamples()=0;
  virtual SampleCountType getSamplesToInput()=0;
  virtual SampleCountType getSamplesToOutput()=0;
};

struct SBSMSRenderChunk {
  int band;
  SampleCountType samplePos;
  TimeType time;
  float chunkSize;
  float stepSize;
  float chunkPos;
  int length;
  float f0;
  float f1;
};

class TrackPoint;
class Debugger {
public:
  virtual bool shouldAssign(int band, TrackPoint *tp, TrackIndexType index)=0;
  virtual bool shouldScale(int band, TrackPoint *tp, TrackIndexType index, int c)=0;
  virtual bool shouldOnset(int band, TrackPoint *tp, TrackIndexType index)=0;
};

class Track;
class SBSMSRenderer {
 public:
  virtual ~SBSMSRenderer() {}
  virtual void startFrame() {}
  virtual void startTime(const SBSMSRenderChunk &i) {}
  virtual void render(const SBSMSRenderChunk &i, Track *t, Debugger *dbg) {}
  virtual void endTime(const SBSMSRenderChunk &i) {}
  virtual void endFrame() {}
  virtual void end(const SampleCountType &samples) {}
};

enum SBSMSError {
  SBSMSErrorNone = 0,
  SBSMSErrorInvalidRate
};

class SBSMSImp;
class Cache;

class SBSMS {
 public:
  SBSMS(int channels, SBSMSQuality *quality, bool bSynthesize, bool bCache);
  SBSMS(SBSMS *src, bool bSynthesize);
  ~SBSMS();

  long read(SBSMSInterface *iface, audio *buf, long n);
  void addRenderer(SBSMSRenderer *renderer);
  void removeRenderer(SBSMSRenderer *renderer);
  long renderFrame(SBSMSInterface *iface);
  long renderFrameFromCache(SBSMSInterface *iface);
  long synthFromCache(SBSMSInterface *iface, audio *buf, long n, Debugger *dbg);
  long getInputFrameSize();
  void reset(bool bFlushInput);
  void seek(SBSMSInterface *iface, SampleCountType samplePos);
  void setTotalSamples(const SampleCountType &samples);
  SampleCountType getTotalSamples();
  SampleCountType getSamplePos();
  void setLeftPos(SampleCountType pos);
  void setRightPos(SampleCountType pos);

  //blargh
  Cache *getCache(int band);
  int getChannels();
  SBSMSQuality *getQuality();

  SBSMSError getError();
  friend class SBSMSImp;
 protected:
  SBSMSImp *imp;
};

enum SlideType {
  SlideIdentity = 0,
  SlideConstant,
  SlideLinearInputRate,
  SlideLinearOutputRate,
  SlideLinearInputStretch,
  SlideLinearOutputStretch,
  SlideGeometricInput,
  SlideGeometricOutput
};

class SlideImp;

class Slide {
 public:
  Slide(SlideType slideType, float rate0 = 1.0f, float rate1 = 1.0f, const SampleCountType &n = 0);
  ~Slide();
  float getTotalStretch();
  float getStretchedTime(float t);
  float getRate(float t);
  float getStretch(float t);
  float getRate();
  float getStretch();
  void step();
 protected:
  SlideImp *imp;
};

 
class SBSMSInterfaceSlidingImp;

class SBSMSInterfaceSliding : public SBSMSInterface {
public:
  SBSMSInterfaceSliding(Slide *rateSlide, 
                        Slide *pitchSlide, 
                        bool bPitchReferenceInput, 
                        const SampleCountType &samplesToInput, 
                        long preSamples,
                        SBSMSQuality *quality);
  virtual ~SBSMSInterfaceSliding();
  virtual float getStretch(float t);
  virtual float getPitch(float t);
  virtual long getPresamples();
  virtual SampleCountType getSamplesToInput();
  virtual SampleCountType getSamplesToOutput();

  friend class SBSMSInterfaceSlidingImp;
protected:
  SBSMSInterfaceSlidingImp *imp;
};



class SBSMSInterfaceVariableRateImp;

class SBSMSInterfaceVariableRate : public SBSMSInterface {
public:
  SBSMSInterfaceVariableRate(SampleCountType samplesToInput);
  virtual ~SBSMSInterfaceVariableRate();
  virtual float getStretch(float t);
  virtual float getPitch(float t);
  void setRate(float rate);
  void setPitch(float pitch);
  virtual long getPresamples();
  virtual SampleCountType getSamplesToInput();
  virtual SampleCountType getSamplesToOutput();

  friend class SBSMSInterfaceVariableRateImp;
protected:
  SBSMSInterfaceVariableRateImp *imp;
};

class ResamplerImp;

class Resampler {
 public:
  Resampler(SBSMSResampleCB func, void *data, SlideType slideType = SlideConstant);
  ~Resampler();
  long read(audio *audioOut, long frames);
  void reset();
  long samplesInOutput();

 protected:
  ResamplerImp *imp;
};

}

#endif
