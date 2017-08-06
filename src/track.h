// -*- mode: c++ -*-
#ifndef TRACK_H
#define TRACK_H

#include "trackpoint.h"
#include "config.h"
#include "utils.h"

#ifdef MULTITHREADED
#include "pthread.h"
#endif

#include <vector>
using namespace std;

namespace _sbsms_ {

enum {
  trackIndexNone = 0
};


enum {
  TrackStart = 0x1,
  TrackTailStart = 0x2,
  TrackStitchStart = 0x4,
  TrackEnd = 0x8,
  TrackTailEnd = 0x10,
  TrackStitchEnd = 0x20,
  TrackSingleton = 0x40,
  TrackOnset = 0x80,
  TrackOffset = 0x100,
  TrackDummy = 0x200,
  TrackMultiple = 0x400
};

#define PhShift 5
#define WShift 21
#define Ph1 65535
#define WPI 536870912
#define W2PI 1073741824
#define W2PIMask 1073741823
#define WScale 1.708913188941079e8f
#define MScale 4.656683928435187e-10f

enum SynthMode {
  synthModeOutput = 0,
  synthModeTrial
};

class SMS;


class Gen {
public:
  virtual ~Gen() { }
  float w;
  float m1[2];
  float ph1[2];
  float m[2];
  float ph[2];
  int channels;
  virtual void gen(float dw, float *dm, audio *in, audio *out, int n)=0;
};

class Osc1 : public Gen
{
public:
  virtual ~Osc1() {}
  void gen(float dw, float *dm, audio *in, audio *out, int n);
};

class Osc2 : public Gen
{
public:
  virtual ~Osc2() {}
  void gen(float dw, float *dm, audio *in, audio *out, int n);
};



class Track {
public:
  Track(int band, int h, int N2, TrackIndexType index, TrackPoint *p, const TimeType &time, const TimeType &beginTime, bool bStitch);
  ~Track();

  TrackIndexType getIndex();
  bool isFirst(const TimeType &synthtime);
  bool isLast(const TimeType &synthtime);
  TimeType getStart();
  TimeType getEnd();

  virtual void push_back(TrackPoint *p);
  virtual void pop_back();
  virtual void updateM(const TimeType &time);
  virtual void step(const TimeType &time);
  virtual void absorb();
  virtual TrackPoint *cacheTrackPoint(const TimeType &time);
  virtual void endTrack(bool bStitch);


  int band;
  TrackIndexType index;
  TimeType start;
  TimeType first;
  TimeType end;
  TimeType last;
  bool tailStart;
  bool tailEnd;
  bool bStitchStart;
  bool bStitchEnd;

  TrackPoint * &back();
  TrackPoint *getTrackPoint(const TimeType &time);
  TimeType size();  
  void prune(const TimeType &time);
  
  friend class SMS;
  friend class SynthRenderer;
  friend class TrackPoint;
  friend class TrackStitcher;
  friend class CacheRenderer;
 protected:
  Gen *gen;
  vector<TrackPoint*> point;
  int h;
  int N2;
  bool bEnd;
  bool bEnded;
  bool bRender;

};

class AnalysisTrack : public Track 
{
public:
  AnalysisTrack(int band, int h, int N2, TrackIndexType index, TrackPoint *p, const TimeType &time, AnalysisTrack *precursor, int channels);
  virtual ~AnalysisTrack();

  virtual void push_back(TrackPoint *p);
  virtual void updateM(const TimeType &time);
  virtual void step(const TimeType &time);
  virtual void absorb();
  virtual TrackPoint *cacheTrackPoint(const TimeType &time);
  virtual void endTrack(bool bStitchEnd);
  void synth(audio *out, const TimeType &synthtime, int n, float f0, float f1);
  void trimOnset(const TimeType &time, int n, int c);
  void trimOffset(const TimeType &time, int c);
  void trim(const TimeType &time, int c);
  void dummy();
  AnalysisTrackPoint * &back();
  
  AnalysisTrack *precursor;
  bool bDummy;
  int channels;
  bool bOn[2];
  TimeType offsetTime[2];
  vector<float> ss[2];
  vector<float> meritsOn[2];
  vector<float> meritsOff[2];
  Average<float> meanS[2];
  Average<float> meanM[2];
  Average<float> meanJ[2];
  Average<float> meandt[2];
  Diff<float> diffdt2[2];
};

}

#endif
