#include "track.h"
#include "real.h"
#include "dBTable.h"
#include "utils.h"
#include "synthTable.h"
#include <assert.h>
#include <algorithm>
using namespace std;

namespace _sbsms_ {


void Osc2 :: gen(float dw, float *dm, audio *in, audio *out, int n)
{
  int iph0 = lrintf(ph[0] * WScale);
  if(iph0>=W2PI) iph0 -= W2PI;
  else if(iph0<0) iph0 += W2PI;
  
  int iph1 = lrintf(ph[1] * WScale);
  if(iph1>=W2PI) iph1 -= W2PI;
  else if(iph1<0) iph1 += W2PI;
  
  int iw = lrintf(w * WScale);
  int idw = lrintf(dw * WScale);
  audio *end = out + n;

  while(out != end) {
    if(iw < WPI) {
      long f = (iph0>>PhShift)&Ph1;
      long i = iph0>>WShift;
      ((float*)out)[0] += m[0] * (float)(synthTable1[i] + f * synthTable2[i]);
      f = (iph1>>PhShift)&Ph1;
      i = iph1>>WShift;
      ((float*)out)[1] += m[1] * (float)(synthTable1[i] + f * synthTable2[i]);
    }
    iph0 += iw;
    iph1 += iw;
    iw += idw;
    iph0 &= W2PIMask;
    iph1 &= W2PIMask;
    m[0] += dm[0];
    m[1] += dm[1];
    out++;
  }

  w = (float)iw / WScale;
  ph[0] = (float)iph0 / WScale;
  ph[1] = (float)iph1 / WScale;
}

void Osc1 :: gen(float dw, float *dm, audio *in, audio *out, int n)
{
  int iph = lrintf(ph[0] * WScale);
  if(iph>=W2PI) iph -= W2PI;
  else if(iph<0) iph += W2PI;

  int iw = lrintf(w * WScale);
  int idw = lrintf(dw * WScale);
  audio *end = out + n;

  while(out != end) {
    if(iw < WPI) {
      long f = (iph>>PhShift)&Ph1;
      long i = iph>>WShift;
      ((float*)out)[0] += m[0] * (float)(synthTable1[i] + f * synthTable2[i]);
    }
    iph += iw;
    iw += idw;
    iph &= W2PIMask;
    m[0] += dm[0];
    out++;
  }

  w = (float)iw / WScale;
  ph[0] = (float)iph / WScale;

}

Track :: Track(int band, int h, int N2, TrackIndexType index, TrackPoint *p, const TimeType &time, const TimeType &beginTime, bool bStitchStart)
{
  gen = NULL;
  this->band =  band;
  this->N2 = N2;
  this->h = h;
  this->index = index;
  bEnd = false;
  bEnded = false;
  bRender = false;
  first = time;
  start = time;
  bStitchEnd = false;
  if(bStitchStart) {
    this->bStitchStart = true;
    tailStart = false;
  } else {
    this->bStitchStart = false;
    if(start > beginTime) {
      tailStart = true;
      start--;
    } else {
      tailStart = false;
    }
  }
  point.push_back(p);
  end = time;
  last = time;
}

TrackIndexType Track :: getIndex()
{
  return index;
}

bool Track :: isFirst(const TimeType &time)
{
  return (time == first);
}

bool Track :: isLast(const TimeType &time)
{
  return (time == last);
}

TimeType Track :: getStart()
{
  return start;
}

TimeType Track :: getEnd()
{
  return end;
}

TimeType Track :: size()
{
  return point.size();
}

TrackPoint * &Track :: back()
{
  return point.back(); 
}

void Track :: prune(const TimeType &time)
{
  while(point.size() && last > time) {
    if(end > 0 && (last == end)) {
      if(tailEnd) delete point.back();
      last--;
      end = -1;
      tailEnd = 0;
    }
    if(tailStart && (last == start)) {
      delete point.back();
    }
    point.pop_back();
    last--;
  }
}

TrackPoint *Track :: getTrackPoint(const TimeType &time)
{
  if(time >= first)
  return point[time - first];
  else return NULL;
}

void Track :: updateM(const TimeType &time)
{
}


void Track :: step(const TimeType &time)
{
  if(time > first && time < last) {
    TrackPoint *tp = point[time-first];
    point[time-first] = NULL;
  }
}

void Track :: push_back(TrackPoint *p)
{
  point.push_back(p);
  last++;
  end++;
}

void Track :: pop_back()
{
  last--;
  if(bEnded) {
    end -= 2;
  } else {
    end--;
  }
}

void Track :: endTrack(bool bStitchEnd)
{
  if(bStitchEnd) {
    this->bStitchEnd = true;
    tailEnd = false;
  } else {
    this->bStitchEnd = false;
    tailEnd = true;
    end++;
  }

  bEnded = true;
}

void AnalysisTrack :: endTrack(bool bStitchEnd)
{
  Track::endTrack(bStitchEnd);
  if(bStitchEnd) {
    point.back()->flags |= (TrackEnd | TrackStitchEnd);
  } else {
    point.back()->flags |= (TrackEnd | TrackTailEnd);
  }
}

TrackPoint *Track :: cacheTrackPoint(const TimeType &time)
{
  TrackPoint *tp = getTrackPoint(time);
  return tp;
}

void Track :: absorb()
{
}

Track :: ~Track()
{
  if(tailEnd && gen) delete gen;
}

/* For analysis from audio */
AnalysisTrack :: AnalysisTrack(int band, int h, int N2, TrackIndexType index, TrackPoint *p, const TimeType &time, AnalysisTrack *precursor, int channels) : Track(band,h,N2,index,p,time,0,precursor!=NULL)
{
  this->precursor = precursor;
  this->channels = channels;
  bDummy = false;

  meandt[0].setSize(max(2,min(6,N2/h/4)));
  meandt[1].setSize(max(2,min(6,N2/h/4)));
  meanS[0].setSize(max(1,min(5,N2/h/4)));
  meanS[1].setSize(max(1,min(5,N2/h/4)));
  meanJ[0].setSize(max(1,min(4,N2/h/4)));
  meanJ[1].setSize(max(1,min(4,N2/h/4)));
  meanM[0].setSize(1);
  meanM[1].setSize(1);
  
  offsetTime[0] = time;
  offsetTime[1] = time;
  bOn[0] = true;
  bOn[1] = true;

  ((AnalysisTrackPoint*)p)->refCount++;
  ((AnalysisTrackPoint*)p)->owner = this;

  if(tailStart) {
    p->flags |= (TrackStart | TrackTailStart);
  } else {
    p->flags |= (TrackStart | TrackStitchStart);
  }

}

AnalysisTrack :: ~AnalysisTrack() 
{
  for(vector<TrackPoint*>::iterator i = point.begin();
      i != point.end();
      ++i) {
    TrackPoint *tp = (*i);
    if(tp) {
      ((AnalysisTrackPoint*)tp)->destroy();
    }
  }
}

void AnalysisTrack :: step(const TimeType &time)
{
  if(time >= first && time <= last) {
    AnalysisTrackPoint *tp = (AnalysisTrackPoint*) point[time-first];
    tp->destroy();
    point[time-first] = NULL;
  }
}

void AnalysisTrack :: push_back(TrackPoint *p)
{
  point.push_back(p);
  ((AnalysisTrackPoint*)p)->owner = this;
  ((AnalysisTrackPoint*)p)->refCount++;
  last++;
  end++;
}
          
void AnalysisTrack :: updateM(const TimeType &time)
{
  if(time == first && time == start) {
    AnalysisTrackPoint *tp0 = (AnalysisTrackPoint*)getTrackPoint(time);
    tp0->updateM();
  }
  if(time < last) {
    AnalysisTrackPoint *tp1 = (AnalysisTrackPoint*)getTrackPoint(time+1);
    tp1->updateM();
  }
}

void AnalysisTrack :: absorb()
{
  for(vector<TrackPoint*>::iterator i = point.begin();
      i != point.end();
      ++i) {
    AnalysisTrackPoint *tp = (AnalysisTrackPoint*)(*i);
    tp->absorb();
  }
}

// if already exists, this is a shared trackpoint, and the last call
// was on the precursor track
TrackPoint *AnalysisTrack :: cacheTrackPoint(const TimeType &time)
{
  AnalysisTrackPoint *tp = (AnalysisTrackPoint*)getTrackPoint(time);
  if(!tp) return NULL;
  TrackPoint *cache;

  if(tp->cache) { 
    cache = tp->cache;
    cache->flags += TrackMultiple;
  } else {
    cache = new TrackPoint();
    tp->cache = cache;
    cache->flags = 0;
    cache->f = tp->f;
    cache->m[0] = tp->m[0];
    cache->m[1] = tp->m[1];
    cache->ph[0] = tp->ph[0];
    cache->ph[1] = tp->ph[1];
    cache->s[0] = tp->s[0];
    cache->s[1] = tp->s[1];
    cache->dph[0] = tp->dph[0];
    cache->dph[1] = tp->dph[1];
    cache->dt[0] = tp->dt[0];
    cache->dt[1] = tp->dt[1];
    cache->meritOn[0] = tp->meritOn[0];
    cache->meritOn[1] = tp->meritOn[1];
    cache->meritOff[0] = tp->meritOff[0];
    cache->meritOff[1] = tp->meritOff[1];
    cache->bScale[0] = tp->bScale[0];
    cache->bScale[1] = tp->bScale[1];
  }
  if((tailStart && start == time - 1) || (!tailStart && start == time)) {
    cache->flags |= TrackStart;
    if(bStitchStart) {
      cache->flags |= TrackStitchStart;
    }
    if(tailStart) {
      cache->flags |= TrackTailStart;
    }
  }
  if(bEnded && ((tailEnd && end == time + 1) || (!tailEnd && end == time))) {
    cache->flags |= TrackEnd;
    if(bStitchEnd) {
      cache->flags |= TrackStitchEnd;
    }
    if(tailEnd) {
      cache->flags |= TrackTailEnd;
    }
  }
  
  if(size() == 1) {
    if(!(tailEnd && !tailStart)) {
      abort();
    }
    cache->flags |= TrackSingleton;
  }

  if(tp->bOnset[0]) {
    cache->flags |= TrackOnset;
  }
  if(tp->bOffset[0]) {
    cache->flags |= TrackOffset;
  }
  if(bDummy) {
    cache->flags |= TrackDummy;
  }

  return cache;
}


void AnalysisTrack :: synth(audio *out,
                            const TimeType &time,
                            int n,
                            float f0,
                            float f1)
{

  float m0[2], m1[2];
  float w0, w1;
  bool bTailStart;
  bool bTailEnd;
  AnalysisTrackPoint *tp0;
  AnalysisTrackPoint *tp1;

  if(time >= end) return;
  if(time < last) {
    tp1 = (AnalysisTrackPoint*)getTrackPoint(time+1);
    w1 = tp1->f * f1;
    for(int c=0; c<channels; c++) {
      m1[c] = tp1->m[c];
    }
    bTailStart = false;
    bTailEnd = false;
  } else {
    bTailStart = false;
    bTailEnd = (last != end);
  }
  if(time >= first) {
    tp0 = (AnalysisTrackPoint*)getTrackPoint(time);
    w0 = tp0->f * f0;
    for(int c=0; c<channels; c++) {
      m0[c] = tp0->m[c];
    }
  } else {
    bTailStart = true;
  }

  if(bTailEnd) {
    float dm[2];
    long iph[2];
    int fall = min(n,w0==0.0f?384:min(384,(int)lrintf(PI * 4.0f / w0)));
    for(int c=0; c<channels; c++) { 
      dm[c] = m0[c] / fall; 
      iph[c] = lrintf(tp0->phSynth[c] * WScale);
      if(iph[c]>=W2PI) iph[c] -= W2PI;
    }
    float w = w0;
    audio *end = out + fall;
    long iw = lrintf(w * WScale);
    if(channels == 1) {
      while(out != end) {
        if(iw < WPI) {
          long f = (iph[0]>>PhShift)&Ph1;
          long i = iph[0]>>WShift;
          ((float*)out)[0] += m0[0] * (float)(synthTable1[i] + f * synthTable2[i]);
        }
        out++;
        m0[0] -= dm[0];
        iph[0] += iw;
        iph[0] &= W2PIMask;
      }
    } else {
      while(out != end) {
        if(iw < WPI) {
          long f = (iph[0]>>PhShift)&Ph1;
          long i = iph[0]>>WShift;
          ((float*)out)[0] += m0[0] * (float)(synthTable1[i] + f * synthTable2[i]);
          f = (iph[1]>>PhShift)&Ph1;
          i = iph[1]>>WShift;
          ((float*)out)[1] += m0[1] * (float)(synthTable1[i] + f * synthTable2[i]);
        }
        out++;
        m0[0] -= dm[0];
        m0[1] -= dm[1];
        iph[0] += iw;
        iph[0] &= W2PIMask;
        iph[1] += iw;
        iph[1] &= W2PIMask;
      }
    }
  }

  if(bTailStart) {
    float dm[2];
    long iph[2];
    tp1->phSynth[0] = tp1->ph[0];
    tp1->phSynth[1] = tp1->ph[1];
    int rise = min(n,w1==0.0f?384:min(384,(int)lrintf(PI * 3.0f / w1)));
    float w = w1;
    out += n;
    long iw = lrintf(w * WScale);
    for(int c=0; c<channels; c++) { 
      dm[c] = m1[c] / rise; 
      iph[c] = lrintf(tp1->phSynth[c] * WScale);
      if(iph[c]>=W2PI) iph[c] -= W2PI;
    }
    audio *end = out-rise;
    if(channels == 1) {
      while(out != end) {
        out--;
        m1[0] -= dm[0];
        iph[0] -= iw;
        if(iph[0]<0) iph[0] += W2PI;
        if(iw < WPI) {
          long f = (iph[0]>>PhShift)&Ph1;
          long i = iph[0]>>WShift;
          ((float*)out)[0] += m1[0] * (float)(synthTable1[i] + f * synthTable2[i]);
        }
      }
    } else {
      while(out != end) {
        out--;
        m1[0] -= dm[0];
        m1[1] -= dm[1];
        iph[0] -= iw;
        if(iph[0]<0) iph[0] += W2PI;
        iph[1] -= iw;
        if(iph[1]<0) iph[1] += W2PI;
        if(iw < WPI) {
          long f = (iph[0]>>PhShift)&Ph1;
          long i = iph[0]>>WShift;
          ((float*)out)[0] += m1[0] * (float)(synthTable1[i] + f * synthTable2[i]);
          f = (iph[1]>>PhShift)&Ph1;
          i = iph[1]>>WShift;
          ((float*)out)[1] += m1[1] * (float)(synthTable1[i] + f * synthTable2[i]);
        }
      }
    }
  }

  if(!(bTailStart || bTailEnd)) {
    float dw = (w1 - w0) / n;
    float w = w0 + 0.5f * dw;
    float dm[2];
    long iph[2];
    for(int c=0; c<channels; c++) { 
      dm[c] = (m1[c] - m0[c]) / n;
      iph[c] = lrintf(tp0->phSynth[c] * WScale);
      if(iph[c]>=W2PI) iph[c] -= W2PI;
    }
    long iw = lrintf(w * WScale);
    long idw = lrintf(dw * WScale);
    audio *end = out + n;
    if(channels == 1) {
      while(out != end) {
        if(iw < WPI) {
          long f = (iph[0]>>PhShift)&Ph1;
          long i = iph[0]>>WShift;
          ((float*)out)[0] += m0[0] * (float)(synthTable1[i] + f * synthTable2[i]);
        }
        iph[0] += iw;
        iw += idw;
        iph[0] &= W2PIMask;
        m0[0] += dm[0];
        out++;
      }
      tp1->phSynth[0] = (float)iph[0] / WScale;
    } else {
      while(out != end) {
        if(iw < WPI) {
          long f = (iph[0]>>PhShift)&Ph1;
          long i = iph[0]>>WShift;
          ((float*)out)[0] += m0[0] * (float)(synthTable1[i] + f * synthTable2[i]);
          f = (iph[1]>>PhShift)&Ph1;
          i = iph[1]>>WShift;
          ((float*)out)[1] += m0[1] * (float)(synthTable1[i] + f * synthTable2[i]);
        }
        iph[0] += iw;
        iph[0] &= W2PIMask;
        iph[1] += iw;
        iph[1] &= W2PIMask;
        iw += idw;
        m0[0] += dm[0];
        m0[1] += dm[1];
        out++;
      }
      tp1->phSynth[0] = (float)iph[0] / WScale;
      tp1->phSynth[1] = (float)iph[1] / WScale;
    }
  }
}

void AnalysisTrack :: dummy()
{
  bDummy = true;
}


void AnalysisTrack :: trimOnset(const TimeType &onset, int n, int c) 
{
  int k0 = max(0,(int)(onset-first-n));
  float m = 1.0f;
  if(k0 == 0) {
    m = 0.0f;
  } else {
    for(int k=k0; k<onset-first; k++) {
      TrackPoint *tp = point[k];
      if(tp->m[c] < m) {
        m = min(m,tp->m[c]);
        k0 = k;
      }
    }
  }
  
  for(int k=k0; k<onset-first; k++) {
    TrackPoint *tp = point[k];
    //tp->m[c] = m;
    tp->bScale[c] = true;
  }

  AnalysisTrackPoint *tp = (AnalysisTrackPoint*) point[onset-first];
  tp->bOnset[c] = true;
  /*  
  if(!bOn[c]) {
    for(TimeType time = onset-1; time >= offsetTime[c]; time--) {
      TrackPoint *tp = point[time-first];
      tp->m[c] = 0.0f;
    }
  }
  */
  bOn[c] = true;

  //bOnset = true;
  //point.erase(point.begin(),point.begin() + (time - first));
  //first = time;
  //  assert(tailStart);
  //start = max((TimeType)0,first - 1);
}

void AnalysisTrack :: trimOffset(const TimeType &time, int c) 
{
  ((AnalysisTrackPoint*)point[time-first])->bOffset[c] = true;
  bOn[c] = false;
}


void AnalysisTrack :: trim(const TimeType &time, int c) 
{
  if(!bOn[c]) {
    //point[time-first]->m[c] = 0.0f;
    point[time-first]->bScale[c] = true;
  }
}

AnalysisTrackPoint * &AnalysisTrack :: back()
{
  return (AnalysisTrackPoint*&)point.back(); 
}

}
