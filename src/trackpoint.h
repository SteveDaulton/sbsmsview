// -*- mode: c++ -*-
#ifndef TRACKPOINT_H
#define TRACKPOINT_H

#include "sbsms.h"

namespace _sbsms_ {

enum { TrackPointNoCont = 65535 };

class Track;
class AnalysisTrack;
class Slice;

class TrackPointBare {
public:
  float f;
  float m[2];
  float ph[2];
};

struct Merit {
  float e0;
  float e1;
  float emean;
  float j;
  float jmean;
  float dbs;
  float total;
};

class TrackPoint : public TrackPointBare {
public:
  unsigned int flags;
  float s[2];
  float f2[2];
  float dt[2]; // meh?
  float dph[2];
  //float meritOn[2];
  Merit meritOn[2];
  float meritOff[2];
  bool bScale[2];
};

class AnalysisTrackPoint : public TrackPoint
{
public:
  AnalysisTrackPoint(Slice *slice, float *peak, float x, float **mag2, float **dT, audio **gx, int N, int band, int channels);
  ~AnalysisTrackPoint();
  void destroy();
  void absorb();
  void updateM();
  TrackPoint *makeTrackPoint();
protected:
  TrackPoint *cache;
  AnalysisTrackPoint *pp;
  AnalysisTrackPoint *pn;
  AnalysisTrackPoint *dupcont;
  AnalysisTrackPoint *cont;
  AnalysisTrackPoint *dup[3];
  AnalysisTrack *owner;
  Slice *slice;
  float *peak;
  float x01;
  float y01[2];
  float phSynth[2];
  float xtp2;
  float xtn2;
  int refCount;
  int channels;
  float x;
  float dtmax;
  float mTot2;
  float y[2];
  float contF;
  float m2[2];
  bool bConnected;
  bool bConnect;
  bool bDelete;
  bool bOwned;
  bool bMarked;
  bool bSplit;
  bool bMerge;
  bool bOnset[2];
  bool bOffset[2];

  int order; // used in track crossing stuffs
  friend class Slice;
  friend class SMS;
  friend class Track;
  friend class AnalysisTrack;
};

class Slice {
public:
  Slice(int band, const TimeType &time);
  ~Slice();
  void remove(AnalysisTrackPoint *tp);
  AnalysisTrackPoint *bottom;
  AnalysisTrackPoint *top;
  int band;
  TimeType time;
};

}

#endif
