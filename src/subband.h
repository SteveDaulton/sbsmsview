// -*- mode: c++ -*-
#ifndef SUBBAND_H
#define SUBBAND_H

#include "real.h"
#include "buffer.h"
#include "sms.h"
#include <stdio.h>
#include "config.h"
#ifdef MULTITHREADED
#include "pthread.h"
#endif

namespace _sbsms_ {

enum {
  NDownSample = 256,
  SDownSample = 4,
  subBufSize = 512,
  hSub = NDownSample/(2*SDownSample)
};

class SubBand {
 public:
  SubBand(SubBand *parent, int band, int channels, SBSMSQuality *quality, bool bAnalyze, bool bSynthesize, SubBand *source, TrackStitcher *stitcher, bool bCreateCache);
  ~SubBand();

  void seek(long framePos, SampleCountType samplePos);
  void reset(bool flushInput);
  void writingComplete();
  long getDrop(float stretch);
  long assignFromCache(float stretch, float pitch);
  long assignFromCache();
  long assignFromCacheInit(bool bSet);
  long synthFromCache(audio *buf, long n, float stretch, float pitch, Debugger *dbg);
  long synthFromCacheInit(long n, float *stretch, Debugger *dbg);
  long synth(long n, float stretch, float pitch, Debugger *dbg);
  void stepSynthGrain(int n);
  
  long write(audio *buf, long n, float stretch, float pitch);
  long read(audio *buf, long n);
  long getInputFrameSize();
  bool isDone();
  void setLeftPos(SampleCountType pos);
  void setRightPos(SampleCountType pos);

  // UGH!
  Cache *getCache(int band);

  bool writeInit();
  long analyzeInit(int,bool,long n=0);
  long extractInit(bool);
  long markInit(bool);
  long assignInit(bool);
  long trialInit(bool);
  long adjustInit(bool);
  long scoreInit(bool);
  long propagateInit(bool);
  long trimInit(bool);
  long renderInit(bool);
  long readInit();
  long renderFromCache();
  long renderFromCacheInit(bool);
  SampleCountType getSamplePos();

  void analyze(int);
  void extract();
  void mark();
  void assign();
  void trial();
  void trialStart();
  void trialTrial();
  void trialEnd();
  void adjust();
  void score();
  void propagate();
  void trim();
  void advance();
  void render();

  void assignStart();
  void assignInit();
  void assignFind();
  bool assignConnect();
  void assignStep();
  void splitMerge();
  void addRenderer(SBSMSRenderer *);
  void removeRenderer(SBSMSRenderer *);
  long renderSynchronous();
  void renderComplete(const SampleCountType &samples);
  void process(bool);

  void stepAnalyzeFrame(int);
  void stepExtractFrame();
  void stepMarkFrame();
  void stepAssignFrame();
  void stepTrialFrame();
  void stepAdjustFrame();
  void stepInitFrame();
  void stepScoreFrame();
  void stepPropagateFrame();
  void stepTrimFrame();
  void stepRenderFrame();
  void stepReadFrame();

#ifdef MULTITHREADED
  pthread_mutex_t dataMutex;
  pthread_mutex_t grainMutex[NumGrainTypes];
#endif
  friend class SBSMSImp;

 protected:
  void init();  
  long getFramesAtFront(int);
  void readSubSamples();
  void setStretch(float stretch);
  void setPitch(float pitch);

  int nMarkLatency;
  int nAssignLatency;
  int nTrialLatency;
  int nAdjustLatency;
  int nScoreLatency;
  int nPropagateLatency;
  int nTrimLatency;
  int nRenderLatency;
  int nRenderLatencyOriginal;
  int nWriteSlack;
  int nExtractSlack;
  int nAnalyzeSlack;
  int nMarkSlack;
  int nAssignSlack;
  int nTrialSlack;
  int nAdjustSlack;
  int nScoreSlack;
  int nPropagateSlack;
  int nTrimSlack;
  int nRenderSlack;

  list<SBSMSRenderer*> renderers;
  RingBuffer<float> stretchRender;
  RingBuffer<float> pitchRender;
  int inputFrameSize;
  RingBuffer<int> outputFrameSize;
  float totalSizef;
  SBSMSQuality *quality;
  int channels;
  int N;
  int h;
  int band;
  long nReadFromOutputFrame;
  long nToWriteForGrain;
  long res;
  long resMask;
  long nGrainsPerFrame;
  long nToDrop1;
  long nToDrop2;
  long nToStart;
  long nToPrepad;
  bool bAnalyze;
  bool bSynthesize;


  bool bBackwards;
  bool bSynthStarted;
  long synthFramePos;
  long bSynthGrainStart;
  long nGrainsSynthed;


  long nGrainsToAnalyze[NumGrainTypes];
  long nGrainsToExtract;
  long nGrainsToMark;
  long nGrainsToAssign;
  long nGrainsToAdvance;
  long nGrainsToTrial;
  long nGrainsToAdjust;
  long nGrainsToScore;
  long nGrainsToPropagate;
  long nGrainsToTrim;
  long nGrainsToRender;
  long nGrainsWritten;
  long nGrainsMarked;
  long nGrainsAssigned;
  long nGrainsTrialed;
  long nGrainsAdjusted;
  long nGrainsScored;
  long nGrainsPropagated;
  long nGrainsTrimmed;
  long nGrainsAdvanced;
  long nGrainsRendered;
  long nGrainsRead;

  long nFramesAnalyzed[NumGrainTypes];
  long nFramesExtracted;
  long nFramesMarked;
  long nFramesAssigned;
  long nFramesTrialed;
  long nFramesAdjusted;
  long nFramesScored;
  long nFramesPropagated;
  long nFramesTrimmed;
  long nFramesRendered;
  long nFramesRead;

  SubBand *parent;
  SubBand *sub;
  SampleBufBase *outMixer;
  SynthRenderer *synthRenderer;
  SMS *sms;
  SampleBuf *samplesSubIn;
  SampleBuf *samplesSubOut;
  GrainBuf *grains[NumGrainTypes];
  GrainBuf *analyzedGrains[NumGrainTypes];
  GrainBuf *grainsIn;
  GrainAllocator *downSampledGrainAllocator;
  SampleCountType samplePos;
  long nToDrop;
  bool bWritingComplete;
  int M;
  int N1;
  int N2;

  bool bCreateCache;
  Cache *cache;
  CacheRenderer *cacheRenderer;
};

}

#endif
