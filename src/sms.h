// -*- mode: c++ -*-
#ifndef SMS_H
#define SMS_H

#include "config.h"
#ifdef MULTITHREADED
#include "pthread.h"
#endif
#include "sbsms.h"
#include "track.h"
#include "grain.h"
#include "buffer.h"
#include <queue>
#include <list>
#include <map>
using namespace std;

namespace _sbsms_ {

enum {
  minTrialBand = 0
};

class GenFactory 
{
public:
  virtual Gen *create()=0;
};

class DefaultGenFactory : public GenFactory 
{
public:
  DefaultGenFactory(int channels);
  int channels;
  Gen *create();
};


class TrackSynthesizer 
{
public:
  virtual void synth( const SBSMSRenderChunk &i, audio *in, audio *out, Track *t, Debugger *dbg)=0;
};


class TrackStitcher : public TrackSynthesizer
{
public:
  TrackStitcher(GenFactory *genFactory, int channels);
  void synth( const SBSMSRenderChunk &i, audio *in, audio *out, Track *t, Debugger *dbg);
  void reset();
  int channels;
  map< TrackPoint*, Gen* > stitch;
  GenFactory *genFactory;

};

class SynthRenderer : public SBSMSRenderer, public SampleBufBase {
 public:
  SynthRenderer(int channels, int h, TrackSynthesizer *synth);
  ~SynthRenderer();
  void startTime( const SBSMSRenderChunk &i);
  void render( const SBSMSRenderChunk &i, Track *t, Debugger *dbg);
  void endTime( const SBSMSRenderChunk &i);
  long read(audio *out, long n);
  void reset(bool bFlush);

 protected:
  TrackSynthesizer *synth;
  int channels;
  ArrayRingBuffer<audio> *sines;
#ifdef MULTITHREADED
  pthread_mutex_t bufferMutex;
#endif
};

class Cache {
public:
  Cache(int channels);
  ~Cache();
  int channels;
  long frames;
  int kLo;
  int kHi;
  int minK;
  int maxK;
  int h;
  SampleBuf wave;
  vector<TrackPoint*> trackPoints;
  vector<TrackIndexType> trackIndex;
  vector<long> indexAtTime;
  vector<int> nTracksAtTime;
  vector<float> meritOn[2];
  vector<float*> mag1Cache[2];
  vector<float*> mag2Cache[2];
  vector<float*> magTrialCache[2];
  vector< list<int> > cuts[2];
};

class CacheRenderer : public SBSMSRenderer {
public:
  CacheRenderer(Cache *cache, int channels);
  void startTime( const SBSMSRenderChunk &i);
  void render( const SBSMSRenderChunk &i, Track *t, Debugger *dbg);
  void endTime( const SBSMSRenderChunk &i);

  Cache *cache;
  int channels;
  int nTracks;
};

class SMS {
 public:
  SMS(SMS *lo, int N, int band, int bandMax, int h, int res, int N1, int N2, int channels, audio *peak2);
  ~SMS();
  void render( float stretch, float f0, float f1, list<SBSMSRenderer*> &renderers);
  void add(grain *g1, grain *g2, grain *gT, Cache *cache);
  void assignStart(long offset);
  void assignInit(long offset);
  void assignFind(long offset);
  bool assignConnect(long offset,  bool bLastDitch);
  void start(long offset);
  void splitMerge();
  void mark(long offset);
  void advance();
  void trial();
  void trialStart();
  void trialEnd();
  void adjust(Cache *cache);
  void adjustInit(ArrayRingBuffer<float> **trialRingBuf,
                  GrainBuf *trialGrainBuf);
  void adjust(GrainBuf *trialGrainBuf,
              queue<float*> *magQueue,
              int minCutSep,
              float **_mag1,
              float **_dmag1,
              audio **x1,
              const TimeType &time,
              Slice *slice,
              Cache *cache);
  void trim(Cache *cache);
  void score(AnalysisTrack *t);
  void score();
  int getTrimLatency();
  int getRenderLatency();
  void propagateScores(long offset);
  void propagateScoresDown(long offset);
  void prepad(audio *buf, long n);
  int getTrialLatency();
  void reset();
  void seek(const TimeType &time);
  bool isPastLeft();
  bool isPastRight();
  void pruneTracks();
  void setLeftPos(SampleCountType pos);
  void setRightPos(SampleCountType pos);
  SampleCountType getSamplePos();
  long synthInit( int n, bool bBackwards, float stretch);
  int synthTracks( int n, SBSMSRenderer *renderer, bool bBackwards, float stretch, float pitch, Debugger *dbg);
  void stepSynth( int n);
  bool assignTrackPointsFromCache( Cache *cache);
  bool assignTrackPointsFromCache( bool bBackwards, int offset, Cache *cache, Debugger *dbg);
  

  
 protected:
  float cost(AnalysisTrackPoint *tp0, AnalysisTrackPoint *tp1);
  float costOffsetMatch(float f, float m2, AnalysisTrackPoint *tp1);
  float costOffsetMatch(AnalysisTrackPoint *tp0, AnalysisTrackPoint *tp1);
  void connect(AnalysisTrackPoint *tp0, AnalysisTrackPoint *tp1, int ilo);
  void mark(long offset, long offsetlo);
  AnalysisTrackPoint *nearestForward(AnalysisTrackPoint **begin, AnalysisTrackPoint *tp0, float *minCost2, float maxCost2, float maxDF2, float dMCoeff2);
  AnalysisTrackPoint *nearestReverse(AnalysisTrackPoint **begin, AnalysisTrackPoint *tp0, float *minCost2, float maxCost2, float maxDF2, float dMCoeff2);
  AnalysisTrackPoint *nearestForwardIncludeOwned(AnalysisTrackPoint **begin, AnalysisTrackPoint *tp0, float *minCost2, float maxCost2, float maxDF2, float dMCoeff2);
  AnalysisTrackPoint *nearestReverseIncludeOwned(AnalysisTrackPoint **begin, AnalysisTrackPoint *tp0, float *minCost2, float maxCost2, float maxDF2, float dMCoeff2);
  AnalysisTrackPoint *nearestForwardMerge(AnalysisTrackPoint **begin, AnalysisTrackPoint *tp0, float *minCost2, float maxCost2, float maxDF, float dMCoeff2); 
  AnalysisTrackPoint *nearestReverseMerge(AnalysisTrackPoint **begin, AnalysisTrackPoint *tp0, float *minCost2, float maxCost2, float maxDF, float dMCoeff2); 
  float interp2(int k, int ko1, float kf);
  float findExtremum(float *mag, int k);
  float findMagnitude(float *mag, float x);
  void calcmags(float *mag, audio *x);
  void calcDT(float *dt, float *mag2, audio *x2, audio *xT);
  int findCut(float *dmag, int k, int maxK);
  AnalysisTrack *createTrack( TrackPoint *tp, const TimeType &time, AnalysisTrack *precursor);
  Track *createTrack( TrackPoint *tp, const TimeType &time, bool bStitch, TrackIndexType index);
  void returnTrackIndex( Track *t);
  void insertTrackForRender( Track *t);
  Track *getTrack( TrackIndexType index);
  void assignTrackPointFromCache(TrackPoint *tp,  bool bBackwards, TrackIndexType index);

  float calcEnergy(Slice *slice, float **mag2, float **dec2, bool bK);
  list<AnalysisTrackPoint*> ended;
  list<AnalysisTrackPoint*> started;
  list<AnalysisTrackPoint*> connected;
  int minTrackSize; 
  int peakWidth1;
  int peakWidth2;
  int minCutSep;
  int minK;
  int maxK;
  float peakThresh;
  float maxCost2;
  float maxDF;
  float dMCoeff2;
  float dNCoeff2;
  float maxCost2Match;
  float maxDFMatch;
  float dMCoeff2Match;
  float maxCost2OffsetMatch;
  float maxDFOffsetMatch;
  float dMCoeff2OffsetMatch;
  float maxCost2SplitMerge;
  float maxDFSplitMerge;
  float dMCoeff2SplitMerge;
  float maxFMatchM;
  float minFMatchL;
  float minFLo;
  float maxFHi;
  float minFMid;
  float maxFMid;
  int kStart;
  int kEnd;
  int kLo;
  int kHi;
  float mNorm;
  float localFavorRatio;
  queue<Slice*> adjustSliceQueue;
  RingBuffer<Slice*> sliceBuffer;
  Slice* sliceM0;
  Slice* sliceL0;
  Slice* sliceH0;
  Slice* sliceM1;
  Slice* sliceL1;
  Slice* sliceM2;
  Slice* sliceH1;
  Slice* sliceH_1;
  float *magTot;
  audio* x10[2];
  audio* x11[2];
  float* dmag1[2];
  float* mag1[2];
  float *mag2[2];
  float *mak2[2];
  float *dT[2];
  audio* x2[2];
  audio* xT[2];
  float* dec2[2];
  float* dek2[2];
  float *peak20;
  float *peak2N;
  int N;
  int Nover2;
  SMS *lo;
  SMS *hi;
  map<TrackIndexType, Track *> liveTracks;
  queue<float*> magQueue[2];
  audio *trialBuf;
  GrainBuf *trialGrainBuf;
  list<AnalysisTrack*> assignTracks;
  list<Track*> renderTracks;
  TimeType addtime;
  TimeType assigntime;
  TimeType adjusttime;
  TimeType trialtime;
  TimeType scoretime;
  TimeType trimtime;
  TimeType synthtime;
  double h2cum;
  int channels;
  long res;
  long resMask;
  int h;
  int dtmax;
  int N1;
  int N2;
  float M;
  double h1;
  int band;  
  float crossingPenalty;
  friend class SubBand;

#ifdef MULTITHREADED
  pthread_mutex_t sliceMutex;
  pthread_mutex_t magMutex[2];
  pthread_mutex_t renderMutex;
  pthread_mutex_t trialMutex;
  pthread_mutex_t trackMutex;
#endif
  
  float bandMeritOn[2];
  float bandMeritOff[2];
  float bandM[2];

  RingBuffer<float> meritBuffer[2];
  RingBuffer<float> meritLoBuffer[2];
  RingBuffer<float> meritHiBuffer[2];

  RingBuffer<float> meritOffBuffer[2];
  RingBuffer<float> meritOffLoBuffer[2];
  RingBuffer<float> meritOffHiBuffer[2];

  RingBuffer<float> mBuffer[2];
  RingBuffer<float> mLoBuffer[2];
  RingBuffer<float> mHiBuffer[2];
  bool bScorePropagated;
  int onsetLength;

  SampleCountType samplePos;
  bool bAssignDone;
  double synthOffset;
  double grainPos;
  double grainLength;
  double samplePosCum;

  TimeType beginTime;
  SampleCountType leftPos;
  SampleCountType rightPos;
  
  long trackSerialNum;
};

}

#endif
