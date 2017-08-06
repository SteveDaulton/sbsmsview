#include "sms.h"
#include "real.h"
#include "utils.h"
#include "dBTable.h"
#include <stdlib.h>
#include <math.h>
#include <set>
#include <assert.h>
#include "synthTable.h"
#include <algorithm>

using namespace std;


#define COND2 (0 && ((band == 4 && assigntime >= 2120 && assigntime <= 2167) || (band == 5 && assigntime >= 2120 && assigntime <= 2167 )))

#define COND0 (0 && ((band == 6 && assigntime >= 907 && assigntime <= 908) || (band == 5 && assigntime >= 1815 && assigntime <= 1816 )))

#define COND (0 && ((band == 3 && assigntime >= 1183 && assigntime <= 1184) || (band == 2 && assigntime >= 2367 && assigntime <= 2368 )))

#define COND1 (0 && ((i.band == 3 && i.time >= 3584 && i.time <= 3586) || (i.band == 4 && i.time >= 3584 && i.time <= 3586 )))

// XXX check for singleton double stitch tracks

namespace _sbsms_ {


bool TrackPointSortFunc(TrackPoint *tp0, TrackPoint *tp1)
{
  return tp0->f < tp1->f;
}

/* Renderers */
// Cache Renderer
// Synth Renderer
// File Renderer


DefaultGenFactory :: DefaultGenFactory(int channels) : channels(channels)
{
}

Gen *DefaultGenFactory :: create()
{
  if(channels == 1) {
    return new Osc1();
  } else {
    return new Osc2();
  }
}

TrackStitcher :: TrackStitcher(GenFactory *genFactory, int channels)
{
  this->genFactory = genFactory;
  this->channels = channels;
}

void TrackStitcher :: reset()
{
  for(map<TrackPoint*, Gen*>::iterator j = stitch.begin(); j != stitch.end(); ++j) {
    delete j->second;
  }
  stitch.clear();
}

// map key's are tp pointers so this works iff all tps are persistent
void TrackStitcher :: synth(const SBSMSRenderChunk &i, audio *in, audio *out, Track *t, Debugger *dbg)
{
  Gen *gen = NULL;
  TrackPoint *tp0 = i.time>=t->first?t->getTrackPoint(i.time):NULL;
  TrackPoint *tp1 = i.time<t->last?t->getTrackPoint(i.time+1):NULL;
  bool bDummy = (tp0?tp0->flags&TrackDummy:false) || (tp1?tp1->flags&TrackDummy:false);
  if(bDummy) return;

  bool bStart = (i.time == t->start);
  bool bEnd = (t->bEnded && i.time + 1 == t->end);
  if(i.chunkPos == 0.0f && bStart) {
    if(t->bStitchStart) {
      map<TrackPoint*, Gen*>::iterator j = stitch.find(tp0);
      if(j != stitch.end()) {
        gen = j->second;
        stitch.erase(j);
      }
      if(!gen) {
        gen = genFactory->create();
        for(int c=0; c<channels; c++) {
          gen->m[c] = tp0->m[c];
          gen->ph[c] = tp0->ph[c];
          if(dbg && dbg->shouldScale(i.band,tp0,t->index,c)) gen->m[c] *= tp0->s[c];
        }
      }
    } else if(t->tailStart) {
      gen = genFactory->create();
      for(int c=0; c<channels; c++) {
        gen->ph[c] = tp1->ph[c];
        gen->m[c] = 0.0f;
      }
    } else {
      gen = genFactory->create();
      for(int c=0; c<channels; c++) {
        gen->ph[c] = tp0->ph[c];
        gen->m[c] = tp0->m[c];
        if(dbg && dbg->shouldScale(i.band,tp0,t->index,c)) gen->m[c] *= tp0->s[c];
      }
    }
    t->gen = gen;
  } else {
    if(t->gen) {
      gen = t->gen;
    } else {
      gen = genFactory->create();
      t->gen = gen;
    }
  }

  float h2 = fabsf(i.chunkSize);
  int n = i.length;
  
  float w0;
  float w1;
  float m0[2];
  float m1[2];
  float dm[2];

  if(t->end > t->last+1) {
    abort();
  }

  if(tp0) {
    for(int c=0; c<channels; c++) {
      m0[c] = tp0->m[c];
      if(dbg && dbg->shouldScale(i.band,tp0,t->index,c)) m0[c] *= tp0->s[c];
    }
    w0 = i.f0 * tp0->f;
  }
  if(tp1) {
    for(int c=0; c<channels; c++) {
      m1[c] = tp1->m[c];
      if(dbg && dbg->shouldScale(i.band,tp1,t->index,c)) m1[c] *= tp1->s[c];
    }
    w1 = i.f1 * tp1->f;
  }
  
  bool bNewOnset = tp1 && (tp1->flags & TrackOnset) && dbg->shouldOnset(i.band,tp1,t->index);
  if(0 && bNewOnset && tp0 && tp1) {
    float dp = tp1->ph - tp0->ph;
    float dp0 = 0.5f*i.stepSize*(w0 + w1);
    float dw = canonPI(dp - dp0)/i.chunkSize;
    w0 += dw;
    w1 += dw;
  }
  
  //bNewOnset = false;
  bool bOnset = (bStart && t->tailStart) || bNewOnset;
  bool bOffset = (bEnd && t->tailEnd) || bNewOnset;

  
  if(i.time < t->end && !bOnset && !bOffset) {
    float dw;
    if(h2 == 0.0f) {
      dw = 0.0f;
      for(int c=0; c<channels; c++) { dm[c] = 0.0f; }
      gen->w = w0 + i.chunkPos * (w1 - w0);
    } else {
      if(i.chunkPos == 0.0f) {
        gen->w = w0;
      } else {
        gen->w = w0 + i.chunkPos * (w1 - w0);
      }
      if(n < 0) {
        float dt = (float)max(-n,(int)lrintf(i.chunkPos * h2));
        for(int c=0; c<channels; c++) { dm[c] = (gen->m[c] - m0[c]) / dt; }
        dw = (gen->w - w0) / dt;
        gen->w -= 0.5f * dw;
      } else {
        float dt = (float)max(n,(int)lrintf(h2 - i.chunkPos * h2));
        for(int c=0; c<channels; c++) { dm[c] = (m1[c] - gen->m[c]) / dt; }
        dw = (w1 - gen->w) / dt;
        gen->w += 0.5f * dw;
      }
    }
    if(n < 0) {
      for(int c=0; c<channels; c++) { dm[c] = -dm[c]; }
      gen->gen(-dw,dm,in,out,-n);
    } else {
      gen->gen(dw,dm,in,out,n);
    }   
  }

  
  if(bOffset) {
    if(i.chunkPos == 0.0f) {
      gen->w = w0;
    }
    if(n < 0) {
      int iStart;
      if(h2 == 0.0f) {
        iStart = 0;
        for(int c=0; c<channels; c++) { dm[c] = 0.0f; }
      } else {
        iStart = lrintf(h2 * i.chunkPos);
        int rise = w0==0.0f?384:min(384,(int)lrintf(PI * 3.0f / w0));
        rise = min(rise,max(-n,iStart));
        for(int c=0; c<channels; c++) { dm[c] = (m0[c] - gen->m[c]) / rise; }
        iStart = max(0,iStart-rise);
      }
      if(iStart < -n) {
        out += iStart;
        in += iStart;
        gen->gen(0,dm,in,out,-n-iStart);
      }
    } else {
      int fall;
      if(i.chunkPos == 0.0f) {
        for(int c=0; c<channels; c++) { gen->m1[c] = gen->m[c]; }
        for(int c=0; c<channels; c++) { gen->ph1[c] = gen->ph[c]; }
      } else {
        for(int c=0; c<channels; c++) { swap(gen->m[c],gen->m1[c]); }
        for(int c=0; c<channels; c++) { swap(gen->ph[c],gen->ph1[c]); }
      }
      if(h2 == 0.0f) {
        for(int c=0; c<channels; c++) { dm[c] = 0.0f; }
        fall = n;
      } else {
        int nFrame = lrintf(h2 - h2 * i.chunkPos);
        fall = w0==0.0f?384:min(384,(int)lrintf(PI * 3.0f / w0));
        fall = min(max(0,(int)lrintf(fall-h2*i.chunkPos)),max(n,nFrame));
        for(int c=0; c<channels; c++) { dm[c] = -gen->m[c] / fall; }
        fall = min(n,fall);
      }
      gen->gen(0,dm,in,out,fall);
      for(int c=0; c<channels; c++) { swap(gen->m[c],gen->m1[c]); }
      for(int c=0; c<channels; c++) { swap(gen->ph[c],gen->ph1[c]); }
    }
  } 

  if(bOnset) {
    if(i.chunkPos == 0.0f) {
      gen->w = w1;
    }
    if(n < 0) {
      int iStart;
      int rise;
      if(h2 == 0.0f) {
        iStart = 0;
        for(int c=0; c<channels; c++) { dm[c] = 0.0f; }
        rise = -n;
      } else {
        iStart = lrintf(h2 - h2 * i.chunkPos);
        rise = w1==0.0f?384:min(384,(int)lrintf(PI * 3.0f / w1));
        rise = min(max(0,(int)lrintf(rise - h2 + h2*i.chunkPos)),max(-n,iStart));
        for(int c=0; c<channels; c++) { dm[c] = -gen->m[c] / rise; }
        rise = min(-n,rise);
      }
      gen->gen(0,dm,in,out,rise);
    } else {
      int iStart;
      if(h2 == 0.0f) {
        iStart = 0;
        for(int c=0; c<channels; c++) { dm[c] = 0.0f; }
      } else {
        iStart = lrintf(h2 - h2 * i.chunkPos);
        int rise = w1==0.0f?384:min(384,(int)lrintf(PI * 3.0f / w1));
        // XXX
        //rise = 384;
        rise = min(rise,max(n,iStart));
        iStart = max(0,iStart-rise);
        if(i.chunkPos == 0.0f) {
          for(int c=0; c<channels; c++) { 
            gen->ph[c] = canon2PI(tp1->ph[c] - w1 * rise);
            gen->m[c] = 0.0f;
          }
        }
        for(int c=0; c<channels; c++) { dm[c] = (m1[c]-gen->m[c]) / rise; }
      }
      if(iStart < n) {
        out += iStart;
        in += iStart;
        gen->gen(0,dm,in,out,n-iStart);
      }
    }
  } 
  

  if(bEnd) {
    if(t->tailEnd) {
      //delete gen;
    } else {
      stitch[tp1] = gen;
    }
  }
}


Cache :: Cache(int channels) : wave(0)
{
  this->channels = channels;
  this->frames = 0;
}

// Make sure not to delete duplicates which exists due to stitching
Cache :: ~Cache() 
{
  for(unsigned int i=0; i<trackPoints.size(); i++) {
    TrackPoint *tp = trackPoints.at(i);
    if(tp->flags >= TrackMultiple) {
      tp->flags -= TrackMultiple;
    } else {
      delete tp;
    }
  }
  for(int c=0; c<channels; c++) {
    for(unsigned int i=0; i<mag1Cache[c].size(); i++) {
      delete magTrialCache[c].at(i);
    }
    for(unsigned int i=0; i<mag1Cache[c].size(); i++) {
      delete mag1Cache[c].at(i);
    }
    for(unsigned int i=0; i<mag2Cache[c].size(); i++) {
      delete mag2Cache[c].at(i);
    }
  }
}

CacheRenderer :: CacheRenderer(Cache *cache, int channels) 
{
  this->channels = channels;
  this->cache = cache;
}

void CacheRenderer :: startTime(const SBSMSRenderChunk &i)
 { 
  nTracks = 0;
  cache->indexAtTime.push_back(cache->trackPoints.size());
}

void CacheRenderer :: render(const SBSMSRenderChunk &i, Track *t, Debugger *dbg)
{
  if(i.time >= t->first && i.time <= t->last) {
    TrackPoint *tp = ((Track*)t)->cacheTrackPoint(i.time);
    cache->trackPoints.push_back(tp);
    cache->trackIndex.push_back(t->index);
    nTracks++;
  }
}

void CacheRenderer :: endTime(const SBSMSRenderChunk &i) 
{
  cache->nTracksAtTime.push_back(nTracks);
}

SynthRenderer :: SynthRenderer(int channels, int h, TrackSynthesizer *synth)
{
  this->channels = channels;
  this->synth = synth;
  sines = new ArrayRingBuffer<audio>(0);
#ifdef MULTITHREADED
  pthread_mutex_init(&bufferMutex,NULL);
#endif
}

SynthRenderer :: ~SynthRenderer()
{
  delete sines;
}

void SynthRenderer :: startTime(const SBSMSRenderChunk &i)
{
  int n = abs(i.length);
  sines->grow(n);
}

void SynthRenderer :: render(const SBSMSRenderChunk &i, Track *t, Debugger *dbg)
{
  this->synth->synth(i,NULL,sines->getWriteBuf(),t,dbg);
}

void SynthRenderer :: endTime(const SBSMSRenderChunk &i)
{
  sines->writePos += i.length;
}

long SynthRenderer :: read(audio *out, long n)
{
  long ret;
#ifdef MULTITHREADED
  pthread_mutex_lock(&bufferMutex);
#endif
  ret = sines->read(out, n);
#ifdef MULTITHREADED
  pthread_mutex_unlock(&bufferMutex);
#endif
  return ret;
}


void SynthRenderer :: reset(bool flushInput)
{
  sines->clear();
}

/* End Renderers */


void SMS :: insertTrackForRender(Track *t)
{
  list<Track*>::reverse_iterator tt0 = renderTracks.rbegin();
  while(tt0 != renderTracks.rend()) {
    Track *t0 = *tt0;
    if(t->start >= t0->start) {
      break;
    }
    tt0++;
  }
  renderTracks.insert(tt0.base(),t);        
  t->bRender = true;
}

// For Analysis
AnalysisTrack *SMS :: createTrack(TrackPoint *tp, const TimeType &time, AnalysisTrack *precursor) 
{
  TrackIndexType index = trackSerialNum++;
  AnalysisTrack *t = new AnalysisTrack(band,h,N2,index,tp,time,precursor,channels);
  assignTracks.push_back(t);
  liveTracks[index] = t;
  return t;
}

// For Cache
Track *SMS :: createTrack(TrackPoint *tp, const TimeType &time, bool bStitch, TrackIndexType index) 
{
  Track *t = new Track(band,h,N2,index,tp,time,beginTime,bStitch);
  liveTracks[index] = t;
  return t;
}

void SMS :: returnTrackIndex(Track *t)
{
  liveTracks.erase(t->index);
}

Track *SMS :: getTrack(TrackIndexType index)
{
  if(liveTracks.find(index) == liveTracks.end()) {
    return NULL;
  } else {
    return liveTracks[index];
  }
}

void SMS :: assignTrackPointFromCache(TrackPoint *tp, bool bBackwards, TrackIndexType index)
{
  Track *t = getTrack(index);
  bool bStart;
  bool bStitchStart;
  bool bTailStart;
  bool bEnd;
  bool bStitchEnd;
  bool bTailEnd;
  if(bBackwards) {
    bStart = (tp->flags & TrackEnd);
    bTailStart = (tp->flags & TrackTailEnd);
    bStitchStart = (tp->flags & TrackStitchEnd);
    bEnd = (tp->flags & TrackStart);
    bTailEnd = (tp->flags & TrackTailStart);
    bStitchEnd = (tp->flags & TrackStitchStart);
  } else {
    bStart = (tp->flags & TrackStart);
    bTailStart = (tp->flags & TrackTailStart);
    bStitchStart = (tp->flags & TrackStitchStart);
    bEnd = (tp->flags & TrackEnd);
    bTailEnd = (tp->flags & TrackTailEnd);
    bStitchEnd = (tp->flags & TrackStitchEnd);
  }
  if(t) {
    t->push_back(tp);
    if(bEnd) {
      if(t->bEnded) {
        abort();
      }      
      t->endTrack(bStitchEnd);
    }
  } else {
    t = createTrack(tp,assigntime,!(bStart&&bTailStart),index);
    if(tp->flags & TrackSingleton) {
      t->endTrack(false);
    }
    insertTrackForRender(t);  
  }
}

bool SMS :: assignTrackPointsFromCache(Cache *cache)
{
  long time = (long)addtime;
  if((long)cache->indexAtTime.size() > time) {
    long trackPointIndex = cache->indexAtTime.at(time);
    int nTracks = cache->nTracksAtTime.at(time);
    for(int k=0;k<nTracks;k++) {
      TrackPoint *tp = cache->trackPoints.at(trackPointIndex+k);
      TrackIndexType index = cache->trackIndex.at(trackPointIndex+k);
      assignTrackPointFromCache(tp,false,index);
    }
    assigntime++;
    addtime++;
    return true;
  } else {
    return false;
  }
}

bool SMS :: assignTrackPointsFromCache(bool bBackwards, int offset, Cache *cache, Debugger *dbg)
{
  long time = (long)addtime + offset;
  if((long)cache->indexAtTime.size() > time && time >= 0) {
    long trackPointIndex = cache->indexAtTime.at(time);
    int nTracks = cache->nTracksAtTime.at(time);
    for(int k=0;k<nTracks;k++) {
      TrackPoint *tp = cache->trackPoints.at(trackPointIndex+k);
      TrackIndexType index = cache->trackIndex.at(trackPointIndex+k);
      if(!dbg || dbg->shouldAssign(band,tp,index))
        assignTrackPointFromCache(tp,bBackwards,index);
    }
    for(list<Track*>::iterator tt=renderTracks.begin(); 
        tt != renderTracks.end();
        tt++) {
      Track *t = (*tt);      
      if(!t->bEnded && t->last < assigntime) {
         t->endTrack(false);
      }
    }
    assigntime++;
    return true;
  } else {
    if(time < 0) {
      samplePos = 0;
      samplePosCum = 0.0;
    } else if((long)cache->indexAtTime.size() <= time) {
      samplePos = (SampleCountType)cache->indexAtTime.size() * (SampleCountType)(h * M);
      samplePosCum = 0.0;
    }
    return false;
  }
}

SMS :: SMS(SMS *lo, int N, int band, int bandMax, int h, int res, int N1, int N2, int channels, audio *peak2)
{
  this->lo = lo;
  if(lo) lo->hi = this;
  hi = NULL;
  this->band = band;
  this->h = h;
  this->h1 = (double)(h<<band);
  this->res = res;
  this->resMask = res - 1;
  this->channels = channels;
  this->N = N;
  this->Nover2 = N/2;
  this->N1 = N1;
  this->N2 = N2;
  float pad2 = (float)N/(float)N2;
  float pad1 = (float)N/(float)N1;

  dtmax = min(N2/4,2*h);

  onsetLength = N2/h/2;
  bScorePropagated = false;
  crossingPenalty = 0.25f;
  M = (float)(1<<band);
  peakThresh = 1e-8f;
  //float maxDF2 = square(0.005f * (float)h) / M;
  float maxDF2 = square(0.00075f * (float)h);
  maxDF = sqrt(maxDF2);
  maxCost2 = 1.5f * maxDF2;
  dMCoeff2 = 0.0002f * maxDF2;

  maxDFMatch = .06f / M;
  float maxDF2Match = square(maxDFMatch);
  dMCoeff2Match = 0.0009f * maxDF2Match;
  maxCost2Match = .3f * maxDF2Match;

  maxDFOffsetMatch = .12f / M;
  float maxDF2OffsetMatch = square(maxDFOffsetMatch);
  dMCoeff2OffsetMatch = 0.001f * maxDF2OffsetMatch;
  maxCost2OffsetMatch = .6f * maxDF2OffsetMatch;

  float maxDF2SplitMerge = square(0.001f * (float)h) / M;
  maxDFSplitMerge = sqrt(maxDF2SplitMerge);
  maxCost2SplitMerge = 1.0f * maxDF2SplitMerge;
  dMCoeff2SplitMerge = 0.006f * maxDF2SplitMerge;
  
  maxDFSplitMerge = maxDF;
  maxCost2SplitMerge = maxCost2;
  dMCoeff2SplitMerge = dMCoeff2;

  int peakWidth0 = lrintf(pad1*(float)N*0.0055f) + 1;
  peakWidth1 = lrintf(pad1*(float)N*0.0055f) + 1;
  peakWidth2 = lrintf(pad2*(float)N*0.0055f) + 1;

  //int peakWidth0 = lrintf(pad1*(float)N*0.011f) + 1;
  //peakWidth1 = lrintf(pad1*(float)N*0.011f) + 1;
  //peakWidth2 = lrintf(pad2*(float)N*0.011f) + 1;

  minTrackSize = max(384/(h<<band),N2/h/2);
  
  //minTrackSize = 1;

  // trackpoints between kLo and kHi
  minCutSep = max((int)lrintf(0.008f * (float)N),peakWidth0);
  if(band==bandMax) kLo = 1;
  else kLo = max(1L,lrintf(floor(0.5f*(float)N/(float)lo->N*(float)lo->kHi-maxDFMatch*M/TWOPI*(float)N)));
  if(band==0) kHi = Nover2;
  else kHi = max(1L,lrintf(0.4785f * N)-peakWidth0*2);
  kStart = max(1,kLo-peakWidth0);
  kEnd = min(Nover2-1,kHi+peakWidth0*2);
  float kNorm = TWOPI / (float)(M * N);
  maxFHi = (float)kHi * kNorm + maxDF;
  minFLo = (float)kLo * kNorm - maxDF;
  if(lo) maxFMatchM = (float)lo->kHi * TWOPI / (float)(lo->N * M * 2) + maxDFMatch + maxDF;
  else maxFMatchM = 0.0f;

  minFMatchL = (float)kLo * kNorm - maxDFMatch;
  if(lo) maxFMid = (float)lo->kHi * TWOPI / (float)(lo->N * M * 2) + maxDF;
  else maxFMid = 0.0f;
  if(lo) lo->minFMid = (float)kLo * kNorm - lo->maxDF;
  if(lo && lo->lo) {
    minK = max(1L,(lrintf(0.25f * (float)N / (float)lo->lo->N * (float)lo->lo->kHi + peakWidth0)));
  } else {
    minK = 1;
  }
  maxK = min(kEnd,kHi + peakWidth0);
  beginTime = 0;

  localFavorRatio = 1.1f;
  mNorm = MScale * MScale * 16.061113032124002f * pad2 / square((float)N);
  //mNorm = MScale * MScale * 22.345537015550892f * pad2 / square((float)N);
  

  leftPos = 0;
  rightPos = 0;
  //bAssignDone = false;
  addtime = 0;
  assigntime = 0;
  trialtime = 0;
  adjusttime = 0;
  synthtime = 0;
  scoretime = 0;
  trimtime = 0;
  samplePos = 0;
  synthOffset = 0.0;
  grainPos = 0.0;
  h2cum = 0.0;

  magTot = (float*)malloc((Nover2+1)*sizeof(float));  
  trialBuf = new audio[h*res];
  trialGrainBuf = new GrainBuf(N,h,N1,hannpoisson);
  for(int c=0; c<channels; c++) {
    dmag1[c] = (float*)malloc(N*sizeof(float));
    mag1[c] = (float*)malloc((Nover2+1)*sizeof(float));
    x11[c] = (audio*)malloc(N*sizeof(audio));
    x10[c] = (audio*)malloc(N*sizeof(audio));
    x2[c] = (audio*)malloc(N*sizeof(audio));    
    xT[c] = (audio*)malloc(N*sizeof(audio));    
    dT[c] = (float*)malloc((Nover2+1)*sizeof(float));
    mag2[c] = (float*)malloc((Nover2+1)*sizeof(float));
    dec2[c] = (float*)malloc(N*sizeof(float));

    mak2[c] = (float*)malloc((Nover2+1)*sizeof(float));
    dek2[c] = (float*)malloc(N*sizeof(float));
#ifdef MULTITHREADED
    pthread_mutex_init(&magMutex[c],NULL);
#endif
  }
#ifdef MULTITHREADED
  pthread_mutex_init(&renderMutex,NULL);
  pthread_mutex_init(&trackMutex,NULL);
  pthread_mutex_init(&sliceMutex,NULL);
  pthread_mutex_init(&trialMutex,NULL);
#endif
  peak20 = (float*)calloc(2*N,sizeof(float));
  peak2N = peak20 + N;
  for(int k=-Nover2;k<=Nover2;k++) {
    peak2N[k] = norm2(peak2[(k+N)%N]);
  }
}

SMS :: ~SMS()
{
  reset();
  for(int c=0;c<channels;c++) {
    while(!magQueue[c].empty()) {
      delete magQueue[c].front();
      magQueue[c].pop();
    }
    free(x10[c]);
    free(x11[c]);
    free(x2[c]);
    free(dec2[c]);
    free(mag2[c]);
    free(dmag1[c]);
    free(mag1[c]);
  }
  set<Slice*> slices;
  while(!adjustSliceQueue.empty()) {
    slices.insert(adjustSliceQueue.front());
    adjustSliceQueue.pop();
  }
  for(long k=sliceBuffer.readPos; k<sliceBuffer.writePos; k++) {
    slices.insert(sliceBuffer.read(k));
  }
  for(set<Slice*>::iterator i = slices.begin();
      i != slices.end();
      ++i) {
    Slice *s = *i;
    AnalysisTrackPoint *tp = s->bottom;
    delete s;
    while(tp) {
      AnalysisTrackPoint *tpn = tp->pn;
      if(!tp->owner) tp->destroy();
      tp = tpn;
    }
  }
  free(peak20);
  delete trialGrainBuf;
}

void SMS :: reset()
{
  liveTracks.clear();
  set<Track*> tracks;
  for(list<AnalysisTrack*>::iterator i=assignTracks.begin(); 
      i != assignTracks.end(); 
      ) {
    list<AnalysisTrack*>::iterator eraseMe = i;
    tracks.insert(*i);
    i++;
    assignTracks.erase(eraseMe);
  }
  for(list<Track*>::iterator i=renderTracks.begin(); 
      i != renderTracks.end(); 
      ) {
    list<Track*>::iterator eraseMe = i;
    tracks.insert(*i);
    i++;
    renderTracks.erase(eraseMe);
  }
  for(set<Track*>::iterator i=tracks.begin(); 
      i != tracks.end(); 
      ++i ) {
    returnTrackIndex(*i);
    delete *i;
  }
  samplePos = 0;
  samplePosCum = 0.0;
  h2cum = 0.0;
  synthOffset = 0.0;
  grainPos = 0.0f;
  grainLength = h1;
  addtime = 0;
  assigntime = 0;
  trialtime = 0;
  adjusttime = 0;
  scoretime = 0;
  trimtime = 0;
  synthtime = 0;
  trackSerialNum = 0;
}

void SMS :: seek(const TimeType &time)
{
  samplePos = ((SampleCountType)time) * (SampleCountType)(h * M);
  addtime = time;
  assigntime = time;
  synthtime = time;
  beginTime = time;
}

void SMS :: trialStart()
{
  if(band >= minTrialBand) {
    memset(trialBuf,0,h*res*sizeof(audio));
  }
}

void SMS :: trialEnd()
{
  if(band < minTrialBand) return;
#ifdef MULTITHREADED
  pthread_mutex_lock(&trialMutex[c]);
#endif
  trialGrainBuf->write(trialBuf,h*res);
#ifdef MULTITHREADED
  pthread_mutex_unlock(&trialMutex[c]);
#endif
}

void SMS :: trial()
{
#ifdef MULTITHREADED
  pthread_mutex_lock(&trackMutex);
#endif
  for(list<Track*>::iterator tt = renderTracks.begin(); 
      tt != renderTracks.end();
      ++tt) {
    AnalysisTrack *t = (AnalysisTrack*)(*tt);
    if(trialtime >= t->start) {
      if(trialtime > t->last) continue;
      t->updateM(trialtime);
      if(hi && hi->band >= minTrialBand) {
        float f = 0.5f*M;
        t->synth(hi->trialBuf,trialtime,h<<1,f,f);
      }
      if(lo && lo->band >= minTrialBand) {
        float f = 2.0f*M;
        t->synth(lo->trialBuf+(trialtime&(res*lo->res-1))*(h>>1),trialtime,h>>1,f,f);
      }
      if(band >= minTrialBand) {
        float f = M;
        t->synth(trialBuf+(trialtime&resMask)*h,trialtime,h,f,f);
      }
    } else {
      break;
    }
  }
#ifdef MULTITHREADED
  pthread_mutex_unlock(&trackMutex);
#endif
  trialtime++;
}

void SMS :: adjust(Cache *cache)
{
  Slice* slice;
#ifdef MULTITHREADED
  pthread_mutex_lock(&sliceMutex);
#endif
  slice = adjustSliceQueue.front(); adjustSliceQueue.pop();
#ifdef MULTITHREADED
  pthread_mutex_unlock(&sliceMutex);
#endif

  if(band >= minTrialBand) {
    adjust(trialGrainBuf,magQueue,minCutSep,mag1,dmag1,x11,adjusttime,slice,cache);
  }
  delete slice;
  adjusttime++;
}

int SMS :: findCut(float *dmag, int k0, int maxK) 
{
  int k;
  for(k = max(1,k0); k <= maxK; k++) {
    float dd0 = dmag[k+1] - dmag[k];
    if(dd0 > 0.0f) {
      float d02 = square(dmag[k+1] + dmag[k]);
      if(dd0 * square(dmag[k] + dmag[k-1]) > (dmag[k] - dmag[k-1]) * d02
         &&
         dd0 * square(dmag[k+2] + dmag[k+1]) > (dmag[k+2] - dmag[k+1]) * d02) {
        break;
      }
    }
  }
  return k;
}

void SMS :: adjust(GrainBuf *trialGrainBuf,
                   queue<float*> *magQueue,
                   int minCutSep,
                   float **_mag1,
                   float **_dmag1,
                   audio **x1,
                   const TimeType &time,
                   Slice *slice,
                   Cache *cache)
{
  grain *g = trialGrainBuf->read(trialGrainBuf->readPos);
  g->analyze();
  for(int c=0; c<channels; c++) {    
    AnalysisTrackPoint *p = slice->bottom;
    if(c == 0) {
      c2even(g->x, x1[0], N);
    } else {
      c2odd(g->x, x1[1], N);
    }
    float *mag1 = _mag1[c];
    calcmags(mag1, x1[c]);

#ifdef MULTITHREADED
    pthread_mutex_lock(&magMutex[c]);
#endif
    float *mag0 = magQueue[c].front(); magQueue[c].pop();
#ifdef MULTITHREADED
    pthread_mutex_unlock(&magMutex[c]);
#endif
    list<int> cuts;
    if(p) {
      // mag0 or mag1
      float *dmag = _dmag1[c];
      int k3 = min(Nover2,maxK+2);
      dmag[0] = mag1[0];
      for(int k=max(1,minK-2); k<k3; k++) {
        dmag[k] = mag1[k] - mag1[k-1];
      }
      int k = minK;
      while(true) {
        k = findCut(dmag,k+1,maxK);
        if(k >= maxK) {
          break;
        } else  {
          cuts.push_back(k);
        }
      }
      bool bDone = false;
      while(!bDone) {
        bDone = true;
        for(list<int>::iterator i = cuts.begin();
            i != cuts.end();
            ++i) {
          int k0 = *i;
          list<int>::iterator ibad = cuts.end();
          list<int>::iterator i2 = i;
          ++i2;
          float maxY = 0.0f;
          for(;
              i2 != cuts.end();
              ++i2) {
            int k2 = *i2;
            if(k2 - k0 >= minCutSep) break;
            float y = mag0[k2] * mag1[k2];
            if(y >= maxY) {
              maxY = y;
              ibad = i2;
            }
            k0 = k2;
          }
          if(ibad != cuts.end()) {
            if(mag0[*i] * mag1[*i] > maxY) {
              ibad = i;
            }
            cuts.erase(ibad);
            bDone = false;
            break;
          }
        }
      }
      cuts.push_front(minK);
      cuts.push_back(maxK);

      list<int>::iterator i = cuts.begin();
      while(p) {
        int k0 = *i;
        ++i;
        if(i == cuts.end()) break;
        int k2 = *i;
        if(p->x > k2) continue;
        float m0 = 0.0f;
        float m1 = 0.0f;
        for(int k=k0;k<=k2;k++) {
          m0 += mag0[k];
          m1 += mag1[k];
        }
        float s = (m1>m0?sqrt(m0/m1):1.0f);
        //float s = (m1>0.0f)?sqrt(m0/m1):1.0f;
        while(p && p->x <= k2) {
          p->s[c] = s;
          //p->m *= s;
          p = p->pn;
        }
      }
    }

    if(cache) {
      float *mag1dup = (float*)malloc((Nover2+1)*sizeof(float));
      memcpy(mag1dup,mag1,(Nover2+1)*sizeof(float));
      cache->magTrialCache[c].push_back(mag1dup);
      cache->cuts[c].push_back(cuts);
    }

    free(mag0);
  }
  trialGrainBuf->advance(1);
}

void SMS :: render(float stretch, float f0, float f1, list<SBSMSRenderer*> &renderers)
{
  TimeType time = synthtime;
  double h2 = stretch * h1;
  h2cum += h2;
  int h2i = lrint(h2cum);
  h2cum -= h2i; 

  SBSMSRenderChunk chunk;
  chunk.band = band;
  chunk.samplePos = samplePos;
  chunk.time = time;
  chunk.stepSize = h1;
  chunk.chunkSize = h2i;
  chunk.chunkPos = 0;
  chunk.length = h2i;
  chunk.f0 = f0;
  chunk.f1 = f1;

  for(list<SBSMSRenderer*>::iterator i = renderers.begin(); i != renderers.end(); ++i) {
    SBSMSRenderer *renderer = *i;
    renderer->startTime(chunk);
  }
#ifdef MULTITHREADED
  pthread_mutex_lock(&trackMutex);
#endif
  for(list<Track*>::iterator tt = renderTracks.begin(); 
      tt != renderTracks.end();) {
    Track *t = (*tt);
    if(t->bEnded && time > t->last) {
      list<Track*>::iterator eraseMe = tt;
      ++tt;
      renderTracks.erase(eraseMe);
      returnTrackIndex(t);
      delete t;
    } else if(time >= t->start) {
      if(time <= t->last) {
        for(list<SBSMSRenderer*>::iterator i = renderers.begin(); i != renderers.end(); ++i) {
          SBSMSRenderer *renderer = *i;
          renderer->render(chunk,t,NULL);
          if(t->end > t->last+1) {
            abort();
          }
        }
        t->step(time);
      }
      ++tt;
    } else {
      break;
    }
  }
#ifdef MULTITHREADED
  pthread_mutex_unlock(&trackMutex);
#endif  
  for(list<SBSMSRenderer*>::iterator i = renderers.begin(); i != renderers.end(); ++i) {
    SBSMSRenderer *renderer = *i;  
    renderer->endTime(chunk);
  }
  samplePos += h2i;
  synthtime++;
}


AnalysisTrackPoint *SMS :: nearestForward(AnalysisTrackPoint **begin, AnalysisTrackPoint *tp0, float *minCost2, float maxCost2, float maxDF, float dMCoeff2) 
{
  *minCost2 = TrackPointNoCont;
  float minF = tp0->f - maxDF;
  float maxF = tp0->f + maxDF;
  float maxDF2 = square(maxDF);
  while((*begin) && (*begin)->f < minF) {
    (*begin) = (*begin)->pn;
  }
  AnalysisTrackPoint *mintp1 = NULL;
  for(AnalysisTrackPoint *tp1 = (*begin);
      tp1;
      tp1 = tp1->pn) {
    if(0&&COND) {
      printf("%g %g %d %g %g\n",tp0->f,tp1->f,tp1->bOwned,minF,maxF);
    }


    if(tp1->bOwned) continue;
    float df2 = square(tp1->f - tp0->f);
    if(df2 > maxDF2) break;
    float dM2 = dB2Approx(tp1->mTot2,tp0->mTot2);
    float cost2 = (df2+dMCoeff2*dM2);
    if(cost2 > maxCost2) continue;
    if(cost2 < (*minCost2)) {
      (*minCost2) = cost2;
      mintp1 = tp1;
    }
  }
  return mintp1;
}

AnalysisTrackPoint *SMS :: nearestReverse(AnalysisTrackPoint **begin, AnalysisTrackPoint *tp0, float *minCost2, float maxCost2, float maxDF, float dMCoeff2) 
{
  *minCost2 = TrackPointNoCont;
  float minF = tp0->f - maxDF;
  float maxF = tp0->f + maxDF;
  float maxDF2 = square(maxDF);
  while((*begin) && (*begin)->f > maxF) {
    (*begin) = (*begin)->pp;
  }
  AnalysisTrackPoint *mintp1 = NULL;
  for(AnalysisTrackPoint *tp1 = (*begin);
      tp1;
      tp1 = tp1->pp) {
    if(COND) printf("%d %g\n",tp1->bOwned,tp1->f);             
    if(tp1->bOwned) continue;
    float df2 = square(tp1->f - tp0->f);
    if(df2 > maxDF2) break;
    float dM2 = dB2Approx(tp1->mTot2,tp0->mTot2);
    float cost2 = (df2+dMCoeff2*dM2);
    if(cost2 > maxCost2) continue;
    if(cost2 < (*minCost2)) {
      (*minCost2) = cost2;
      mintp1 = tp1;
    }
  }
  return mintp1;
}

AnalysisTrackPoint *SMS :: nearestForwardIncludeOwned(AnalysisTrackPoint **begin, AnalysisTrackPoint *tp0, float *minCost2, float maxCost2, float maxDF, float dMCoeff2) 
{
  *minCost2 = TrackPointNoCont;
  float minF = tp0->f - maxDF;
  float maxF = tp0->f + maxDF;
  float maxDF2 = square(maxDF);
  while((*begin) && (*begin)->f < minF) {
    (*begin) = (*begin)->pn;
  }
  AnalysisTrackPoint *mintp1 = NULL;
  for(AnalysisTrackPoint *tp1 = (*begin);
      tp1;
      tp1 = tp1->pn) {
    float df2 = square(tp1->f - tp0->f);
    if(df2 > maxDF2) break;
    float dM2 = dB2Approx(tp1->mTot2,tp0->mTot2);
    float cost2 = (df2+dMCoeff2*dM2);
    if(cost2 > maxCost2) continue;
    if(cost2 < (*minCost2)) {
      (*minCost2) = cost2;
      mintp1 = tp1;
    }
  }
  return mintp1;
}

AnalysisTrackPoint *SMS :: nearestReverseIncludeOwned(AnalysisTrackPoint **begin, AnalysisTrackPoint *tp0, float *minCost2, float maxCost2, float maxDF, float dMCoeff2) 
{
  *minCost2 = TrackPointNoCont;
  float minF = tp0->f - maxDF;
  float maxF = tp0->f + maxDF;
  float maxDF2 = square(maxDF);
  while((*begin) && (*begin)->f > maxF) {
    (*begin) = (*begin)->pp;
  }
  AnalysisTrackPoint *mintp1 = NULL;
  for(AnalysisTrackPoint *tp1 = (*begin);
      tp1;
      tp1 = tp1->pp) {
    float df2 = square(tp1->f - tp0->f);
    if(df2 > maxDF2) break;
    float dM2 = dB2Approx(tp1->mTot2,tp0->mTot2);
    float cost2 = (df2+dMCoeff2*dM2);
    if(cost2 > maxCost2) continue;
    if(cost2 < (*minCost2)) {
      (*minCost2) = cost2;
      mintp1 = tp1;
    }
  }
  return mintp1;
}

/*
 Note that the tp0 track may acquire a trackpoint from a different band
 and the stitched track may begin with trackpoints from the tp0 band
 */

void SMS :: connect(AnalysisTrackPoint *tp0, AnalysisTrackPoint *tp1, int ilo)
{
  if(COND) {
    printf("connect!  %lld %lld %lld %d %d %g %g %p\n",assigntime,tp0->slice->time,tp1->slice->time,band,tp1->slice->band,tp0->f,tp1->f,tp0->dupcont);
    }

  tp0->cont = tp1;
  connected.push_back(tp0);
  
  TimeType time = assigntime;
  if(tp0->slice->band == tp1->slice->band) {
#ifdef MULTITHREADED
    pthread_mutex_lock(&trackMutex);
#endif    
    tp0->owner->push_back(tp1);
#ifdef MULTITHREADED
    pthread_mutex_unlock(&trackMutex);
#endif
  } else if(tp0->slice->band < tp1->slice->band) {
    AnalysisTrack *precursor = tp0->owner;

    if(ilo == 1) {
#ifdef MULTITHREADED
      pthread_mutex_lock(&trackMutex);
#endif
      precursor->push_back(tp1);
      precursor->endTrack(true);
      TimeType time = precursor->end/res;
#ifdef MULTITHREADED
      pthread_mutex_unlock(&trackMutex);
#endif
#ifdef MULTITHREADED
      pthread_mutex_lock(&lo->trackMutex);
#endif
      AnalysisTrack *t = lo->createTrack(tp1,time,precursor);
#ifdef MULTITHREADED
      pthread_mutex_unlock(&lo->trackMutex);
#endif
    } else {
#ifdef MULTITHREADED
      pthread_mutex_lock(&trackMutex);
#endif
      TimeType time = precursor->end/res;
      precursor->endTrack(true);
      AnalysisTrackPoint *last = (AnalysisTrackPoint*)precursor->back();
#ifdef MULTITHREADED
      pthread_mutex_unlock(&trackMutex);
#endif
#ifdef MULTITHREADED
      pthread_mutex_lock(&lo->trackMutex);
#endif
      AnalysisTrack *t = lo->createTrack(last,time,precursor);
      t->push_back(tp1);
#ifdef MULTITHREADED
      pthread_mutex_unlock(&lo->trackMutex);
#endif
      last->owner = precursor;
    }
  } else {
    AnalysisTrack *precursor = tp0->owner;
#ifdef MULTITHREADED
    pthread_mutex_lock(&trackMutex);
#endif
    precursor->push_back(tp1);
    precursor->endTrack(true);
    TimeType time = precursor->end*hi->res;
#ifdef MULTITHREADED
    pthread_mutex_unlock(&trackMutex);
#endif
#ifdef MULTITHREADED
    pthread_mutex_lock(&hi->trackMutex);
#endif
    hi->createTrack(tp1,time,precursor);
#ifdef MULTITHREADED
    pthread_mutex_unlock(&hi->trackMutex);
#endif
  }

  tp0->bConnected = true;
  tp1->bConnected = true;
  tp0->bOwned = true;
  tp1->bOwned = true;
  if(tp0->dupcont) {
    AnalysisTrackPoint *dup = tp0->dupcont;
    if(!dup->owner) {
      dup->bOwned = true;
      dup->bDelete = true;
    }
  }

  // unused dups before or above the conncetion will be deleted
  // So if a res==2 downstitch is made, M0->L1
  // M0-M1-M2-M3
  // L0    L1   
  // M1 = L1->dup[0]; M2 = L1->dup[1]; M2 = L1->dup[2] will be deleted
  // 
  // If a res==2 upstitch is made, M0->H1
  // H0-H1-H2
  // M0    M1
  //
  // If a res stitch is made
  for(int d=0;d<3;d++) {
    AnalysisTrackPoint *dup = tp1->dup[d];
    if(dup && !dup->owner && (d<2 || dup->slice->band < tp1->slice->band)) {
      dup->bOwned = true;
      dup->bDelete = true;
    }
  }
  AnalysisTrackPoint *dup2 = tp0->dup[2];
  if(dup2 && dup2 != tp1 && !dup2->owner) {
    dup2->bOwned = true;
    dup2->bDelete = true;
  }

  if(COND) {
    printf("connect %d %d %lld %g %g %lld %lld %p \n",tp0->slice->band,tp1->slice->band,assigntime, tp0->f, tp1->f, tp0->owner->start, tp0->owner->end,tp0->dupcont);

  }

  
}

void SMS :: mark(long offset)
{
  mark(offset,0);
  if(offset&resMask) {
    mark(offset,1);
  }
}

void SMS :: mark(long offset, long offsetlo)
{
  if(!lo) return;
#ifdef MULTITHREADED
  pthread_mutex_lock(&lo->sliceMutex);
#endif
  Slice *sliceL1 = lo->sliceBuffer.read(lo->sliceBuffer.readPos+offset/res+offsetlo);
#ifdef MULTITHREADED
  pthread_mutex_unlock(&lo->sliceMutex);
#endif
#ifdef MULTITHREADED
  pthread_mutex_lock(&sliceMutex);
#endif
  Slice *sliceM1 = sliceBuffer.read(sliceBuffer.readPos+offset);
#ifdef MULTITHREADED
  pthread_mutex_unlock(&sliceMutex);
#endif
  bool b0 = !(offset&resMask);
  bool bDone = false;
  bool bLastDitch = false;
  /*
  float maxDFMatch;
  float maxCost2Match;
  float dMCoeff2Match;
  if(b0) {
    maxDFMatch = this->maxDFMatch;
    dMCoeff2Match = 0.002f * square(maxDFMatch);
    maxCost2Match = .5f * square(maxDFMatch);
  } else {
    maxDFMatch = this->maxDFMatch + maxDF;
    dMCoeff2Match = 0.002f * square(maxDFMatch);
    maxCost2Match = .5f * square(maxDFMatch);
  }
  */

  while(!bDone) {
    int nToCont = 0;
    int nCont = 0;
    AnalysisTrackPoint *rbegin = NULL;
    //M0   M1      M0   M1
    //     |   or         \
    //     L1               L1
    AnalysisTrackPoint *begin = sliceL1->bottom;
    for(AnalysisTrackPoint *tp = sliceM1->bottom;
        tp;
        tp = tp->pn) {
      if(tp->bMarked) continue;
      if(tp->f > maxFMatchM) {
        break;
      } else {
        rbegin = tp;
      }
      float F;
      tp->cont = nearestForward(&begin,tp,&F,maxCost2Match,maxDFMatch,dMCoeff2Match);
      if(tp->cont) nToCont++;
    }
    if(sliceL1) {
      for(AnalysisTrackPoint *tp = sliceL1->top;
          tp;
          tp = tp->pp) {
        if(tp->f < minFLo) break;
        float F;
        tp->cont = nearestReverse(&rbegin,tp,&F,maxCost2Match,maxDFMatch,dMCoeff2Match);
      }
    }

    for(AnalysisTrackPoint *tp0 = sliceM1->bottom;
        tp0;
        tp0 = tp0->pn) {
      if(tp0->bMarked) continue;
      if(tp0->f > maxFMatchM) {
        break;
      }
      AnalysisTrackPoint *tp1 = tp0->cont;
      if(tp1) {
        if(bLastDitch || tp1->cont == tp0) {
          nCont++;
          bool bAlreadyMarked = false;
          if(b0) {
            if(tp1->dup[1] || tp0->dup[1]) {
              bAlreadyMarked = true;
            }
          } else {
            if(tp1->dup[2-2*offsetlo] || tp0->dup[2*offsetlo]) {
              bAlreadyMarked = true;
            }
          }
          if(!bAlreadyMarked) {
            if(b0) {
              tp1->dup[1] = tp0;
              tp0->dup[1] = tp1;
            } else {
              tp1->dup[2-2*offsetlo] = tp0;
              tp0->dup[2*offsetlo] = tp1;
            }
            TimeType t = assigntime + offset;
            if(band ==5 && t >= 888 && t <= 894) {
              printf("mark %lld %d %d %g %g %g %g %ld %ld %d \n",t,tp0->slice->band, tp1->slice->band, tp0->f, tp1->f, tp0->ph[0], tp1->ph[0], offset, offsetlo,b0);
            }
          } else {
            TimeType t = assigntime + offset;
            if(band ==5 && t >= 888 && t <= 894) {
              printf("no mark %lld %d %d %g %g %ld %ld %d \n",t,tp0->slice->band, tp1->slice->band, tp0->f, tp1->f, offset, offsetlo,b0);
            }
          }
          tp0->bMarked = true;
        }
      }
    }
    bDone = (nToCont == nCont);
    bLastDitch = (!bDone && nCont==0);
  }
}

void SMS :: assignStart(long offset)
{
  connected.clear();
  //bAssignDone = false;
#ifdef MULTITHREADED
  pthread_mutex_lock(&sliceMutex);
#endif
  sliceM0 = sliceBuffer.read(sliceBuffer.readPos+offset);
  sliceM1 = sliceBuffer.read(sliceBuffer.readPos+offset+1);
  if(res == 2) {
    sliceM2 = sliceBuffer.read(sliceBuffer.readPos+offset+2);
  } else {
    sliceM2 = NULL;
  }
#ifdef MULTITHREADED
  pthread_mutex_unlock(&sliceMutex);	
#endif
  for(AnalysisTrackPoint *tp = sliceM0->bottom;
      tp;
      tp = tp->pn) {
    if(!tp->owner->bEnded) {
      tp->owner->bEnd = true;
      tp->bConnected = false;
      tp->bOwned = false;
    } else {
      tp->bConnected = true;
      tp->bOwned = true;
    }
  }
#ifdef MULTITHREADED
  if(hi) pthread_mutex_lock(&hi->sliceMutex);
#endif
  sliceH_1 = hi&&hi->res==2&&(offset>0)?hi->sliceBuffer.read(hi->sliceBuffer.readPos+(offset)*hi->res-1):NULL;  
  sliceH0 = hi?hi->sliceBuffer.read(hi->sliceBuffer.readPos+(offset)*hi->res):NULL;  
  sliceH1 = hi?hi->sliceBuffer.read(hi->sliceBuffer.readPos+(offset+1)*hi->res):NULL;  
#ifdef MULTITHREADED 
  if(hi) pthread_mutex_unlock(&hi->sliceMutex);
#endif
#ifdef MULTITHREADED
  if(lo) pthread_mutex_lock(&lo->sliceMutex);
#endif
  sliceL0 = lo?lo->sliceBuffer.read(lo->sliceBuffer.readPos+offset/res):NULL;
  sliceL1 = lo?lo->sliceBuffer.read(lo->sliceBuffer.readPos+offset/res+1):NULL;
#ifdef MULTITHREADED
	if(lo) pthread_mutex_unlock(&lo->sliceMutex);
#endif
}

void SMS :: assignInit(long offset)
{
  for(AnalysisTrackPoint *tp = sliceM1->bottom;
      tp;
      tp = tp->pn) {
    tp->cont = NULL;
    tp->contF = TrackPointNoCont;
  }
  if(sliceM2) {
    for(AnalysisTrackPoint *tp = sliceM2->bottom;
        tp;
        tp = tp->pn) {
      tp->cont = NULL;
      tp->contF = TrackPointNoCont;
    }
  }
}

void SMS :: assignFind(long offset)
{
  //if(bAssignDone) return;

  if(COND) printf("assignFind %d %lld\n",band,assigntime);
  Slice *sliceM0 = this->sliceM0;
  Slice *sliceM1 = this->sliceM1;
  Slice *sliceM2 = this->sliceM2;
  Slice *sliceL1 = this->sliceL1;
  Slice *sliceH1 = this->sliceH1;
  AnalysisTrackPoint *begin;

  // M0 <- M1
  begin = sliceM0->bottom;
  for(AnalysisTrackPoint *tp = sliceM1->bottom;
      tp;
      tp = tp->pn) {
    if(tp->bOwned) continue;
    float F;	
    tp->bConnect = false;
    AnalysisTrackPoint *minM = nearestForward(&begin,tp,&F,maxCost2,maxDF,dMCoeff2);
    if(minM && F < tp->contF) {
      tp->cont = minM;
      tp->contF = F;
    }
  }
  if(sliceL1) {
    // M0 <- L1
    AnalysisTrackPoint *rbegin = sliceM0->top;
    for(AnalysisTrackPoint *tp = sliceL1->top;
        tp;
        tp = tp->pp) {
      if(tp->bOwned) continue;
      if(tp->f < minFLo) break;
      float F;
      AnalysisTrackPoint *minL = nearestReverse(&rbegin,tp,&F,maxCost2,maxDF,dMCoeff2);
      if(minL) {
        F *= localFavorRatio;
        if(F < tp->contF) {
          tp->cont = minL;
          tp->contF = F;
        }
      }
      if(COND) {
        printf("m0<-l1 %lld %g %g %g\n",assigntime,tp->f,minL?minL->f:0,tp->cont?tp->cont->f:0);
      }
    }
  }
  begin = sliceM0->bottom;
  if(sliceH1) {
    // M0 <- H1
    for(AnalysisTrackPoint *tp = sliceH1->bottom;
        tp;
        tp = tp->pn) {
      if(tp->bOwned) continue;
      if(tp->f > maxFHi) break;
      float F;
      AnalysisTrackPoint *minH = nearestForward(&begin,tp,&F,maxCost2,maxDF,dMCoeff2);
      if(minH) {
        F *= localFavorRatio;
        if(F < tp->contF) {
          tp->cont = minH;
          tp->contF = F;
        }
      }
    }
  }
  if(sliceM2 && !(offset&resMask)) {
    // M1 <- M2
    // include tp from M1 that are owned
    begin = sliceM1->bottom;
    for(AnalysisTrackPoint *tp = sliceM2->bottom;
        tp;
        tp = tp->pn) {
      if(tp->bOwned) continue;
      float F;	
      tp->bConnect = false;
      AnalysisTrackPoint *minM = nearestForwardIncludeOwned(&begin,tp,&F,maxCost2,maxDF,dMCoeff2);
      if(COND) {
        printf("m1<-m2 %lld %g %g\n",assigntime,tp->f,minM?minM->f:0);
      }
      if(minM) {
        tp->cont = minM;
        tp->contF = F;
      }
    }
    // XXX
    if(sliceL1) {
      // M1 <- L1
      AnalysisTrackPoint *rbegin = sliceM1->top;
      for(AnalysisTrackPoint *tp = sliceL1->top;
          tp;
          tp = tp->pp) {
        if(tp->bOwned) continue;
        if(tp->f < minFLo) break;
        float F;	
        AnalysisTrackPoint *minM = nearestReverseIncludeOwned(&rbegin,tp,&F,maxCost2,maxDF,dMCoeff2);
        if(COND) {
          printf("m1<-l1 %lld %g %g\n",assigntime,tp->f,minM?minM->f:0);
        }
        if(minM) {
          F *= localFavorRatio;
          if(F < tp->contF) {
            tp->cont = minM;
            tp->contF = F;
          }
        }
      }
    }
  }
}

/* There are two cases where dups may not be found.  
If M0->H1 where H res == 2, then H1 tp1->cont may be in H(0.5), with (f,m) approximately interpolating tp0 in M0 and tp1 in H1.  This is not strictly a dup, and it can't straightforwardly be predicted ahead of time whether this will be the case.  
Similarly, if M0->L1 and M res == 2, then there should be a tp1->cont in M1 which interpolates M0 and L1.  

So, if
*/
bool SMS :: assignConnect(long offset, bool bLastDitch)
{
  //if(bAssignDone) return false;
  if(COND) printf("assignConnect %d %lld\n",band,assigntime);
  Slice *sliceM0 = this->sliceM0;
  Slice *sliceM1 = this->sliceM1;
  Slice *sliceL1 = this->sliceL1;
  Slice *sliceH1 = this->sliceH1;
  int nToCont = 0;
  int nCont = 0;
  bool b0 = !(offset&resMask);
  int ilo;
  if(res == 2 && b0) {
    ilo = 0;
  } else {
    ilo = 1;
  }
  AnalysisTrackPoint *beginM1 = sliceM1->bottom;
 AnalysisTrackPoint *beginH1;
  if(sliceH1) beginH1 = sliceH1->bottom;
  // M0 -> M1,L1,H1
  for(AnalysisTrackPoint *tp = sliceM0->bottom;
      tp;
      tp = tp->pn) {
    if(tp->bOwned) continue;
    float FM1 = TrackPointNoCont;
    float FL1 = TrackPointNoCont;
    float FH1 = TrackPointNoCont;
    // M0 -> M1
    // M0 -> L1
    // M0 -> H1
    AnalysisTrackPoint *minM1 = nearestForward(&beginM1,tp,&FM1,maxCost2,maxDF,dMCoeff2);
    AnalysisTrackPoint *minL1 = NULL;
    if(sliceL1 && tp->f < maxFMid) {
      AnalysisTrackPoint *rbeginL1 = sliceL1->top;
      minL1 = nearestReverse(&rbeginL1,tp,&FL1,maxCost2,maxDF,dMCoeff2);
      FL1 *= localFavorRatio;
    }
    AnalysisTrackPoint *minH1 = NULL;
    if(sliceH1 && tp->f > minFMid) {
      minH1 = nearestForward(&beginH1,tp,&FH1,maxCost2,maxDF,dMCoeff2);
      if(COND) {
        //printf("%g %g\n",tp->f,minH1?minH1->f:NULL);
      }

      FH1 *= localFavorRatio;
    }
    if(minM1 &&
       ((FM1<=FH1 && FM1<=FL1)
        ||(minL1 && FL1<=FH1 && FL1<=FM1 && minL1->dup[ilo] == minM1)
        ||(minH1 && FH1<=FL1 && FH1<=FM1 && minH1->dup[1] == minM1))) {
      if(ilo == 1 && minL1 && minL1->dup[1] == minM1) {
        tp->dupcont = minL1;
      } else if(minH1 && minH1->dup[1] == minM1) {
        tp->dupcont = minH1;
      } else {
        tp->dupcont = NULL;
      }
      tp->contF = FM1;
      tp->cont = minM1;
      nToCont++;
    } else if(minL1 && FL1<=FM1 && FL1<=FH1) {
      if(minM1 && minL1->dup[ilo] == minM1) {
        tp->dupcont = minM1;
      } else {
        tp->dupcont = NULL;
      }
      tp->contF = FL1;
      tp->cont = minL1;
      nToCont++;
    } else if(minH1 && FH1<=FM1 && FH1<=FL1) {
      if(minM1 && minH1->dup[1] == minM1) {
        tp->dupcont = minM1;
      } else {
        tp->dupcont = NULL;
      }
      tp->contF = FH1;
      tp->cont = minH1;
      nToCont++;
    } else {
      tp->cont = NULL;
    }
    
    // XXX YYY
    if(0 && tp->cont) {
      // H0          H0.5           H1
      // M0                         M1
      //
      //           tp1->cont        tp1
      // tp                                           
      AnalysisTrackPoint *tp1 = tp->cont;
      if(hi && hi->res == 2 && tp1->slice->band < tp->slice->band && tp1->cont && tp1->cont->slice->band == tp1->slice->band && (tp1->cont->flags & TrackStart)) {
        float cost2 = costOffsetMatch(0.5f*(tp->f + tp1->f), 0.5f*(tp->mTot2 + tp1->mTot2), tp1->cont);        
        if(COND) 
          printf("%lld %d %g %g %g %g\n",assigntime,band,tp->f,tp1->f,tp1->cont->f,cost2/maxCost2Match);

        if(cost2 < maxCost2OffsetMatch) {
          bool bMatch = true;
          if(tp->dup[2]) {
            float cost2Dup = costOffsetMatch(tp,tp->dup[2]);
            if(cost2 > cost2Dup) {
              bMatch = false;
            }
            if(COND) 
              printf("%g %g\n",cost2Dup,cost2);
          }
          if(bMatch && tp1->cont->dup[0]) {
            float cost2Dup = costOffsetMatch(tp1->cont,tp1->cont->dup[0]);
            if(cost2 > cost2Dup) {
              bMatch = false;
            }
            if(COND) 
              printf("%g %g\n",cost2Dup,cost2);
          }


          if(bMatch) {
            if(tp->dup[2]) tp->dup[2]->dup[0] = NULL;
            if(tp1->cont->dup[0]) tp1->cont->dup[0]->dup[2] = NULL;          
            tp->dup[2] = tp1->cont;
            tp1->cont->dup[0] = tp;
            tp1->cont = tp;
          }
        }
      }
      // M0          M1             M2 
      //                            L1
      //
      // tp  tp1->dup[1]->cont   tp1->dup[1]
      //                           tp1
      else if(res == 2 && ilo == 0 && tp1->cont == tp && tp1->slice->band > tp->slice->band && tp1->dup[1] && tp1->dup[1]->cont && tp1->dup[1]->slice->band == tp1->dup[1]->cont->slice->band) {
        float cost2 = costOffsetMatch(0.5f*(tp->f + tp1->f), 0.5f*(tp->mTot2 + tp1->mTot2), tp1->dup[1]->cont);
        if(cost2 < maxCost2OffsetMatch) {
          bool bMatch = true;
          if(tp1->dup[1]->cont->dup[2]) {
            float cost2Dup = costOffsetMatch(tp1->dup[1]->cont,tp1->dup[1]->cont->dup[2]);
            if(cost2 > cost2Dup) {
              bMatch = false;
            }
          }
          if(bMatch && tp1->dup[0]) {
            float cost2Dup = costOffsetMatch(tp1,tp1->dup[0]);
            if(cost2 > cost2Dup) {
              bMatch = false;
            }
          }
          if(bMatch) {
            if(tp1->dup[1]->cont->dup[2]) tp1->dup[1]->cont->dup[2]->dup[0] = NULL;
            if(tp1->dup[0]) tp1->dup[0]->dup[2] = NULL;
            tp1->dup[0] = tp1->dup[1]->cont;
            tp1->dup[0]->dup[2] = tp1;
          }
        }
      }
    }
    
    if(COND) {
      printf("connect?  %lld %d %g %g\n",assigntime,band,tp->f,tp->contF);

      if(tp->cont) {
        printf("\t %d %g %p\n",tp->cont->slice->band,tp->cont->f,tp->dupcont);

        for(int d=0;d<3;d++) {
          AnalysisTrackPoint *dup = tp->cont->dup[d];
          if(dup) {
            printf("dup %d %d %d %g\n",d,tp->cont->slice->band,dup->slice->band,dup->f);
          }
        }


        if(tp->cont->cont) {
          printf("\t %d %g %g %p\n",tp->cont->cont->slice->band,tp->cont->cont->f,tp->cont->cont->contF,tp->cont->cont->dupcont);
        }
      }
    }
  }
  for(AnalysisTrackPoint *tp0 = sliceM0->bottom;
      tp0;
      tp0 = tp0->pn) {

    if(tp0->bOwned) continue;
    tp0->bConnect = false;
    AnalysisTrackPoint *tp1 = tp0->cont;
    TimeType time = assigntime;
    if(tp1 && !tp1->bOwned &&
       (bLastDitch ||
        (tp1->cont == tp0) ||
        ((tp1->cont && tp0->contF <= tp1->cont->contF) &&
         ((tp1->cont->dup[0] == tp0) ||
          (tp1->cont->dup[1] == tp0))))) {
      tp1->cont = tp0;
      tp0->bConnect = true;
      tp1->bConnect = true;
    }
  }
  for(AnalysisTrackPoint *tp0 = sliceM0->bottom;
      tp0;
      tp0 = tp0->pn) {
    if(tp0->bOwned) continue;
    AnalysisTrackPoint *tp1 = tp0->cont;
    if(tp0->bConnect && tp1 && !tp1->bOwned && tp1->bConnect && tp1->cont == tp0) {
      AnalysisTrackPoint *dupcont = tp0->dupcont;
      if(dupcont && dupcont->bConnect) {
        if(!tp1->bConnected && !dupcont->bConnected) {
          if(!tp0->bConnected && (dupcont->cont == NULL || tp0->contF <= dupcont->cont->contF)) {
            nCont++;
            connect(tp0,tp1,ilo);
            tp0->owner->bEnd = false;
            dupcont->bConnect = false;
          } else if(dupcont->cont && !dupcont->cont->bConnected) {
            nCont++;
            connect(dupcont->cont,dupcont,ilo);
            dupcont->cont->owner->bEnd = false;
            tp1->bConnect = false;
          }
        }
      } else if(!tp0->bConnected && !tp1->bConnected) {
        nCont++;
        connect(tp0,tp1,ilo);
        tp0->owner->bEnd = false;
      }
    }
  }
  bool bAssignDone = (nToCont == nCont || bLastDitch);
  return !(bAssignDone || nCont == 0);
}

extern bool bBlob;
void SMS :: start(long offset)
{
  started.clear();
  ended.clear();
#ifdef MULTITHREADED
  pthread_mutex_lock(&trackMutex);
#endif
  for(list<AnalysisTrack*>::iterator tt = assignTracks.begin(); 
      tt != assignTracks.end(); ) {
    AnalysisTrack *t = (*tt);
    bool bKeep;  
    if(t->bEnded) {
      bKeep = ((!t->bRender) && t->size() > 1 && (t->bStitchStart || t->bStitchEnd || t->size() >= minTrackSize));
      if(assigntime > t->last) {
        list<AnalysisTrack*>::iterator eraseMe = tt;
        ++tt;
        assignTracks.erase(eraseMe);
      } else {
        ++tt;
      }
    } else if(t->bEnd) {
      bKeep = (t->size() >= minTrackSize || t->bStitchStart);
      if(bKeep) {
        list<AnalysisTrack*>::iterator eraseMe = tt;
        ++tt;
        // will check in splitMerge stage to see if this track should be continued
        bKeep = !t->bRender;
        t->endTrack(false);
        if(t->size() >= 2) {
          ended.push_back((AnalysisTrackPoint*)t->back());
        } else {
          assignTracks.erase(eraseMe);
        }
      } else {
        list<AnalysisTrack*>::iterator eraseMe = tt;
        ++tt;
        assignTracks.erase(eraseMe);
        if(COND2) {
          
          printf("absorb %d %lld \n",band,assigntime);
          bBlob = true;
        }
        if(t->size() >= 2) {
          t->dummy();
          t->endTrack(false);
          t->absorb();
        bBlob = false;

          bKeep = !t->bRender;
        } else {
          returnTrackIndex(t);
          t->absorb();
        bBlob = false;

          delete t;
          continue;
        }
      }
    } else {
      bKeep = ((!t->bRender) && (t->size() >= minTrackSize || (t->bStitchStart && t->size() > 1)));
      ++tt;
    }
    if(bKeep) {
      //if(band == 3) printf("%lld %g\n",assigntime,t->point[assigntime-t->first-1]->f);
      insertTrackForRender(t);
    }
  }
#ifdef MULTITHREADED
  pthread_mutex_unlock(&trackMutex);
#endif
#ifdef MULTITHREADED
  pthread_mutex_lock(&sliceMutex);
#endif
  Slice *sliceM0 = sliceBuffer.read(sliceBuffer.readPos+offset);
  adjustSliceQueue.push(sliceM0);
#ifdef MULTITHREADED
  pthread_mutex_unlock(&sliceMutex);	
#endif
  for(AnalysisTrackPoint *tp = sliceM0->bottom;
      tp;) {
    AnalysisTrackPoint *tpn = tp->pn;
    if(tp->bOwned) {

      if(tp->bDelete) {
        tp->destroy();
      }
    } else {
      AnalysisTrack *t = createTrack(tp,assigntime,NULL);
     started.push_back(tp);
      for(int d=0;d<2;d++) {
        AnalysisTrackPoint *dup = tp->dup[d];
        if(COND) {
           printf("start? %lld %d %g %p %p %d\n",assigntime,band,tp->f,dup,dup?dup->owner:0,dup?dup->refCount:0);
        }
        if(dup && !dup->owner) {
          dup->destroy();
        }
      }
    }
    tp = tpn;
  }
  assigntime++;
}

float SMS :: costOffsetMatch(AnalysisTrackPoint *tp0, AnalysisTrackPoint *tp1)
{
  return costOffsetMatch(tp0->f,tp0->mTot2,tp1);
}

float SMS :: costOffsetMatch(float f, float m2, AnalysisTrackPoint *tp1)
{
  float maxDF2 = square(maxDFOffsetMatch);
  float df2 = square(tp1->f - f);
  if(df2 > maxDF2) return TrackPointNoCont; 
  float dM2 = dB2Approx(tp1->mTot2,m2);
  return (df2+dMCoeff2OffsetMatch*dM2);        
}

float SMS :: cost(AnalysisTrackPoint *tp0, AnalysisTrackPoint *tp1)
{
  float maxDF2 = square(maxDF);
  float df2 = square(tp1->f - tp0->f);
  if(df2 > maxDF2) return TrackPointNoCont; 
  float dM2 = dB2Approx(tp1->mTot2,tp0->mTot2);
  return (df2+dMCoeff2*dM2);        
}

/* 
The tps must be the start of a track 
*/
AnalysisTrackPoint *SMS :: nearestForwardMerge(AnalysisTrackPoint **begin, AnalysisTrackPoint *tp0, float *minCost2, float maxCost2, float maxDF, float dMCoeff2) 
{
  *minCost2 = TrackPointNoCont;
  float minF = tp0->f - maxDF;
  float maxF = tp0->f + maxDF;
  float maxDF2 = square(maxDF);
  while((*begin) && (*begin)->f < minF) {
    (*begin) = (*begin)->pn;
  }
  AnalysisTrackPoint *mintp1 = NULL;
  for(AnalysisTrackPoint *tp1 = (*begin);
      tp1;
      tp1 = tp1->pn) {
    if(!((tp1->flags & TrackStart) && (tp1->flags & TrackTailStart) && !tp1->owner->bEnded)) continue;
    float df2 = square(tp1->f - tp0->f);
    if(df2 > maxDF2) break;
    float dM2 = dB2Approx(tp1->mTot2,tp0->mTot2);
    float cost2 = (df2+dMCoeff2*dM2);
    if(cost2 > maxCost2) continue;
    if(cost2 < (*minCost2)) {
      (*minCost2) = cost2;
      mintp1 = tp1;
    }
  }
  return mintp1;
}

AnalysisTrackPoint *SMS :: nearestReverseMerge(AnalysisTrackPoint **begin, AnalysisTrackPoint *tp0, float *minCost2, float maxCost2, float maxDF, float dMCoeff2)
{
  *minCost2 = TrackPointNoCont;
  float minF = tp0->f - maxDF;
  float maxF = tp0->f + maxDF;
  float maxDF2 = square(maxDF);
  while((*begin) && (*begin)->f > maxF) {
    (*begin) = (*begin)->pp;
  }
  AnalysisTrackPoint *mintp1 = NULL;
  for(AnalysisTrackPoint *tp1 = (*begin);
      tp1;
      tp1 = tp1->pp) {
    if(!((tp1->flags & TrackStart) && (tp1->flags & TrackTailStart) && !tp1->owner->bEnded)) continue;
    float df2 = square(tp1->f - tp0->f);
    if(df2 > maxDF2) break;
    float dM2 = dB2Approx(tp1->mTot2,tp0->mTot2);
    float cost2 = (df2+dMCoeff2*dM2);
    if(cost2 > maxCost2) continue;
    if(cost2 < (*minCost2)) {
      (*minCost2) = cost2;
      mintp1 = tp1;
    }
  }
  return mintp1;
}

/*
  Look M0<->M0  M0<->L0  M0<->L1 M0<->H0 M0<->H(0.5)
 */


void SMS :: splitMerge()
{
  //printf("%d %lld \n",band,assigntime);
  TimeType time = assigntime - 1;
#ifdef MULTITHREADED
  pthread_mutex_lock(&trackMutex);
#endif
  /*
  // attempt to untangle crossings  
  vector<AnalysisTrackPoint*> left;
  vector<AnalysisTrackPoint*> right;
  for(list<AnalysisTrackPoint*>::iterator i = connected.begin(); 
      i != connected.end();
      ++i) {
    AnalysisTrackPoint *tp0 = *i;
    AnalysisTrackPoint *tp1 = tp0->cont;
    tp1->cont = tp0;
    if(tp0->owner->band == tp1->owner->band) {
      tp0->contF = cost(tp0,tp1);
      left.push_back(tp0);
      right.push_back(tp1);
    }
  }
  
  sort(right.begin(),right.end(),TrackPointSortFunc);
  for(int i=0; i<left.size(); i++) {
    left[i]->order = i;
  }
  for(int i=0; i<right.size(); i++) {
    right[i]->order = i;
  }
  bool bSwap;
  float maxDF2 = square(maxDF);
  do {
    bSwap = false;
    float maxChange = crossingPenalty;
    pair<int,int> maxChangeSwap;
    for(int i=0; i<left.size(); i++) {
      AnalysisTrackPoint *tp0 = left[i];
      int ro = tp0->cont->order;
      for(int j=0; j<ro; j++) {
        AnalysisTrackPoint *tp1 = right[j];
        int ii = tp1->cont->order;
        if(ii > i) {
          if(tp0->slice->band == tp1->slice->band) {
            float change = tp0->contF / cost(tp0,tp1);
            if(change > maxChange) {
              maxChange = change;
              maxChangeSwap = pair<int,int>(i,ii);
              bSwap = true;
            }
          }
        }
      }
    }
    if(bSwap) {
      int i = maxChangeSwap.first;
      int ii = maxChangeSwap.second;
      AnalysisTrackPoint *tp0 = left[i];
      AnalysisTrackPoint *tp1 = left[ii];
      //printf("swapping %d %d %d %d\n",i,ii,tp0->cont->order,tp1->cont->order);
      //printf("swapping %d %lld %g %g %g %g\n",band,assigntime,tp0->f,tp0->cont->f,tp1->f,tp1->cont->f);
      swap(tp0->cont,tp1->cont);
      tp0->cont->cont = tp0;
      tp1->cont->cont = tp1;
      tp0->contF = cost(tp0,tp0->cont);
      tp1->contF = cost(tp1,tp1->cont);
      swap(tp0->owner->back(),tp1->owner->back());
      swap(tp0->owner,tp1->owner);
    }
  } while(bSwap);
  */

  AnalysisTrackPoint *rbeginL0 = sliceL0?sliceL0->top:NULL;
  AnalysisTrackPoint *beginM0 = sliceM0->bottom;
  AnalysisTrackPoint *beginH0 = sliceH0?sliceH0->bottom:NULL;
  AnalysisTrackPoint *beginH_1 = sliceH_1?sliceH_1->bottom:NULL;

  for(list<AnalysisTrackPoint*>::iterator i = ended.begin();
      i != ended.end();
      ++i) {
    AnalysisTrackPoint *tp = *i;
    float F, FL, FH;
    tp->cont = nearestForwardMerge(&beginM0,tp,&F,maxCost2SplitMerge,maxDFSplitMerge,dMCoeff2SplitMerge);
    AnalysisTrackPoint *minL = nearestReverseMerge(&rbeginL0,tp,&FL,maxCost2SplitMerge,maxDFSplitMerge,dMCoeff2SplitMerge);
    if(minL) {
      FL *= localFavorRatio;
      if(FL < F) {
        tp->cont = minL;
        F = FL;
      }
    }
    AnalysisTrackPoint *minH = nearestForwardMerge(&beginH0,tp,&FH,maxCost2SplitMerge,maxDFSplitMerge,dMCoeff2SplitMerge);
    if(minH) {
      FH *= localFavorRatio;
      if(FH < F) {
        tp->cont = minH;
        F = FH;
      }
    }

    AnalysisTrackPoint *minH_1 = nearestForwardMerge(&beginH_1,tp,&FH,maxCost2SplitMerge,maxDFSplitMerge,dMCoeff2SplitMerge);
    if(minH_1) {
      FH *= localFavorRatio;
      if(FH < F) {
        tp->cont = minH_1;
        F = FH;
      }
    }
    AnalysisTrack *t = tp->owner;
    if(1&&tp->cont) {
      AnalysisTrackPoint *tp1 = tp->cont;
      AnalysisTrack *t1 = tp1->owner;
      if(1 && tp1->owner->band == tp->owner->band) {
        assignTracks.remove(t1);
        returnTrackIndex(t1);
        tp1->flags ^= (tp1->flags & (TrackStart | TrackTailStart | TrackStitchStart));
        t->point.pop_back();
        t->end -= 2;
        t->last--;
        for(int c=0; c<channels; c++) {
          tp1->m2[c] += tp->m2[c];
        }
        tp1->mTot2 += tp->mTot2;
        for(vector<TrackPoint*>::iterator i = t1->point.begin(); i != t1->point.end(); ++i) {
          t->push_back(*i);
        }
        t1->point.clear();
        delete t1;
        t->bEnded = false;
        t->bEnd = false;
        tp->destroy();
      } else if(1 && tp1->owner->band < tp->owner->band) {        
        // M0 -> H0 or M0 -> H_1
        if(hi->res == 2) {
          if(t1->first & 1 && t1->size() >= 2) {
            t1->point.erase(t1->point.begin());
            t1->start += 2;
            t1->first++;
            tp1->destroy();
            AnalysisTrackPoint *tp2 = (AnalysisTrackPoint*)t1->point.front();
            tp2->refCount++;
            tp2->flags |= (TrackStart | TrackEnd | TrackStitchStart);
            for(int c=0; c<channels; c++) {
              tp2->m2[c] += tp->m2[c];
            }
            tp2->mTot2 += tp->mTot2;
            tp->destroy();
            t->point.pop_back();
            t->point.push_back(tp2);
            t->tailEnd = false;
            t->bStitchEnd = true;
            t->end--;
            t1->tailStart = false;
            t1->bStitchStart = true;
            t1->precursor = t;
          } else {
            t1->start++;
            tp1->refCount++;
            tp1->flags ^= (tp1->flags & (TrackTailStart));
            tp1->flags |= (TrackStart | TrackEnd | TrackStitchStart);
            for(int c=0; c<channels; c++) {
              tp1->m2[c] += tp->m2[c];
            }
            tp1->mTot2 += tp->mTot2;
            tp->destroy();
            t->point.pop_back();
            t->point.push_back(tp1);
            t->tailEnd = false;
            t->bStitchEnd = true;
            t->end--;
            t1->tailStart = false;
            t1->bStitchStart = true;
            t1->precursor = t;
          }
        } else {
          // M0 -> H0
          t1->start++;
          tp1->refCount++;
          tp1->flags ^= (tp1->flags & (TrackTailStart));
          tp1->flags |= (TrackStart | TrackEnd) | TrackStitchStart;
          tp->destroy();
          t->point.pop_back();
          t->point.push_back(tp1);
          t->tailEnd = false;
          t->bStitchEnd = true;
          t->end--;
          t1->tailStart = false;
          t1->bStitchStart = true;
          t1->precursor = t;
        }
      } else if(1) {
        // M0 -> L0
        if(res == 2) {
          if(t->last & 1) {
            t1->start++;
            tp1->refCount++;
            tp1->flags ^= (tp1->flags & (TrackTailStart));
            tp1->flags |= (TrackStart | TrackEnd | TrackStitchStart);
            t->point.pop_back();
            AnalysisTrackPoint *tp0 = t->back();
            t->point.pop_back();
            for(int c=0; c<channels; c++) {
              tp1->m2[c] += tp0->m2[c];
            }
            tp1->mTot2 += tp0->mTot2;
            tp0->destroy();
            tp->destroy();
            t1->tailStart = false;
            t1->bStitchStart = true;
            if(t->size() == 0 && t->precursor) {
              assignTracks.remove(t);
              renderTracks.remove(t);
              t->precursor->point.pop_back();
              t->precursor->point.push_back(tp1);
              delete t;         
              t1->precursor = t->precursor;     
            } else {
              t->point.push_back(tp1);
              t->end -= 2;
              t->last--;
              t->tailEnd = false;
              t->bStitchEnd = true;
              t1->precursor = t;
            }
          } else {
            t1->start++;
            tp1->refCount++;
            tp1->flags ^= (tp1->flags & (TrackTailStart));
            tp1->flags |= (TrackStart | TrackEnd | TrackStitchStart);
            t->point.pop_back();
            for(int c=0; c<channels; c++) {
              tp1->m2[c] += tp->m2[c];
            }
            tp1->mTot2 += tp->mTot2;
            tp->destroy();
            t->point.push_back(tp1);
            t->end--;
            t->tailEnd = false;
            t->bStitchEnd = true;
            t1->tailStart = false;
            t1->bStitchStart = true;
            t1->precursor = t;
          }
        } else {
          t1->start++;
          tp1->refCount++;
          tp1->flags ^= (tp1->flags & (TrackTailStart));
          tp1->flags |= (TrackStart | TrackEnd | TrackStitchStart);
          t->point.pop_back();
          for(int c=0; c<channels; c++) {
            tp1->m2[c] += tp->m2[c];
          }
          tp1->mTot2 += tp->mTot2;
          tp->destroy();
          t->point.push_back(tp1);
          t->end--;
          t->tailEnd = false;
          t->bStitchEnd = true;
          t1->tailStart = false;
          t1->bStitchStart = true;
          t1->precursor = t;
        }
      }    
    } else {
      // XXX keep iterators so this can be done O(1)
      assignTracks.remove(t);
    }
  }

#ifdef MULTITHREADED
  pthread_mutex_unlock(&trackMutex);
#endif

}

void SMS :: advance()
{
#ifdef MULTITHREADED
  pthread_mutex_lock(&sliceMutex);
#endif
  sliceBuffer.advance(1);
#ifdef MULTITHREADED
  pthread_mutex_unlock(&sliceMutex);
#endif
}

int SMS :: getTrimLatency()
{
  return 3;
}

int SMS :: getRenderLatency()
{
  return 1+onsetLength;
}

// trim needs to determine peaks in the merit score, so it requires a lookahead of at least 2
void SMS :: trim(Cache *cache)
{
  float meritOn[3][2];
  float meritOff[3][2];
  float wBand = 0.5f;
  float meritOnThresh = 9.0f;
  //float wBand = 2.0f;
  //float meritOnThresh = 2.0f;
  float meritOffThresh = 5.0f;

  /*
  if(band <= 4) {
    meritOnThresh = 1e6;
    meritOffThresh = 1e6;
  }
  */
    
  for(int c=0; c<channels; c++) {
    for(int i=0; i<3; i++) {
      float totOn = -meritBuffer[c].get(i);
      totOn += meritHiBuffer[c].get(i);   
      totOn += meritLoBuffer[c].get(i);   

      float totOff = -meritOffBuffer[c].get(i);
      totOff += meritOffHiBuffer[c].get(i);   
      totOff += meritOffLoBuffer[c].get(i);   
      
      float m = -mBuffer[c].get(i);
      m += mHiBuffer[c].get(i);
      m += mLoBuffer[c].get(i);

      meritOn[i][c] = wBand * (m>0.0f?totOn/m:0.0f);
      meritOff[i][c] = wBand * (m>0.0f?totOff/m:0.0f);
    }
    if(cache) cache->meritOn[c].push_back(meritOn[0][c]);
  }

                                
  
  //printf("%d %g %g %g\n",band,trimtime*h1,meritOn[0][0],meritOff[0][0]);

#ifdef MULTITHREADED
  pthread_mutex_lock(&trackMutex);
#endif
  for(list<Track*>::iterator tt = renderTracks.begin(); 
      tt != renderTracks.end();
      ) {
    AnalysisTrack *t = (AnalysisTrack*)(*tt);
    list<Track*>::iterator eraseMe = tt;
    ++tt;

    //AnalysisTrackPoint *tp = (AnalysisTrackPoint*)t->point[trimtime-t->first];
    //if(band ==3 && trimtime == 35 && trimtime >= t->first && trimtime <= t->last)  printf("trime %g\n",tp->f);
    if(t->bEnded && trimtime == t->last) {
      /*
      TimeType time = trimtime;
      TimeType first = max(t->first,time-onsetLength);
      bool bOff[2] = {true,true};
      bool bContinue = true;
      while(bContinue && time >= first) {
        bContinue = false;
        TrackPoint *tp = t->getTrackPoint(time);
        for(int c=0; c<channels; c++) {
          if(bOff[c] && tp->dt[c] > dtmax) {
            tp->m[c] = 0.0f;
            bContinue = true;
          } else {
            bOff[c] = false;
          }
        }
        time--;
      }
      */
    } else if(trimtime >= t->first && trimtime < t->last) {

      bool bOnsetThresh[2] = {false,false};
      bool bOffsetThresh[2] = {false,false};
      float onsetTotal[3] = {0.0f, 0.0f, 0.0f};
      float offsetTotal[3] = {0.0f, 0.0f, 0.0f};


      for(int c=0; c<channels; c++) {
        long pos = trimtime - t->first;
        float meritOnTrack0 = t->meritsOn[c][pos];
        float meritOffTrack0 = t->meritsOff[c][pos];
        if(trimtime + 1 <= t->last) {
          float meritOnTrack1 = t->meritsOn[c][pos + 1];
          float meritOffTrack1 = t->meritsOff[c][pos + 1];
          if(trimtime + 2 <= t->last) {
            float meritOnTrack2 = t->meritsOn[c][pos + 2];
            float meritOffTrack2 = t->meritsOff[c][pos + 2];
            float meritOnTrackThresh = meritOnThresh;
            if(!t->bOn[c]) meritOnTrackThresh *= 0.75f;            
            if(meritOnTrack0 + meritOnTrack1 + meritOnTrack2 >= meritOnTrackThresh) {
              bOnsetThresh[c] = true;
              onsetTotal[0] += meritOnTrack0 + meritOn[0][c];
              onsetTotal[1] += meritOnTrack1 + meritOn[1][c];
              onsetTotal[2] += meritOnTrack2 + meritOn[2][c];
            } else if(meritOffTrack0 + meritOffTrack1 + meritOffTrack2 >= meritOffThresh) {
              bOffsetThresh[c] = true;
              offsetTotal[0] += meritOffTrack0 + meritOff[0][c];
              offsetTotal[1] += meritOffTrack1 + meritOff[1][c];
              offsetTotal[2] += meritOffTrack2 + meritOff[2][c];
            }
          }
        }
      }
      
      if(channels == 1) {
        bool bChange = false;
        if(bOnsetThresh[0] && 
           onsetTotal[1] > onsetTotal[0] && onsetTotal[1] > onsetTotal[2]) {
          t->trimOnset(trimtime+1,onsetLength,0);
          bChange = true;
        } else if(bOffsetThresh[0] && 
                  offsetTotal[1] > offsetTotal[0] && offsetTotal[1] > offsetTotal[2]) {
          t->trimOffset(trimtime+1,0);
          bChange = true;
        }
        if(!bChange) {
          t->trim(trimtime+1,0);
        }
      } else {
        bool bChange = false;
        if((bOnsetThresh[0] && bOnsetThresh[1]) && 
           onsetTotal[1] > onsetTotal[0] && onsetTotal[1] > onsetTotal[2]) {
          t->trimOnset(trimtime+1,onsetLength,0);
          t->trimOnset(trimtime+1,onsetLength,1);
          bChange = true;
        } else if((bOffsetThresh[0] && bOffsetThresh[1]) && 
                  offsetTotal[1] > offsetTotal[0] && offsetTotal[1] > offsetTotal[2]) {
          t->trimOffset(trimtime+1,0);
          t->trimOffset(trimtime+1,1);
          bChange = true;
        }
        if(!bChange) {
          t->trim(trimtime+1,0);
          t->trim(trimtime+1,1);
        }
      }
    }
  }
#ifdef MULTITHREADED
  pthread_mutex_unlock(&trackMutex);
#endif  

  for(int c=0; c<channels; c++) {
    meritBuffer[c].advance(1);
    meritLoBuffer[c].advance(1);
    meritHiBuffer[c].advance(1);
    meritOffBuffer[c].advance(1);
    meritOffLoBuffer[c].advance(1);
    meritOffHiBuffer[c].advance(1);
    mBuffer[c].advance(1);
    mLoBuffer[c].advance(1);
    mHiBuffer[c].advance(1);
  }

  trimtime++;
}

// meritHi contains the sum of merits for bands above and including this band
void SMS :: propagateScoresDown(long offset)
{
  //printf("down %d %lld %ld\n",band,scoretime,offset);          
  if(lo) {
    if(res == 2) {
      long offset2 = offset>>1;
      if(offset & 1) {
        for(int c=0; c<channels; c++) {
          float &loMeritHi0 = lo->meritHiBuffer[c].get(offset2);
          float merit = meritHiBuffer[c].get(offset);   
          loMeritHi0 += 0.25f * merit;
          float &loMeritHi1 = lo->meritHiBuffer[c].get(offset2 + 1);
          loMeritHi1 += 0.25f * merit;

          float &loMeritOffHi0 = lo->meritOffHiBuffer[c].get(offset2);
          float meritOff = meritOffHiBuffer[c].get(offset);   
          loMeritOffHi0 += 0.25f * meritOff;
          float &loMeritOffHi1 = lo->meritOffHiBuffer[c].get(offset2 + 1);
          loMeritOffHi1 += 0.25f * meritOff;

          float &loMHi0 = lo->mHiBuffer[c].get(offset2);
          float m = mHiBuffer[c].get(offset);   
          float &loMHi1 = lo->mHiBuffer[c].get(offset2 + 1);
          loMHi0 += 0.25f * m;
          loMHi1 += 0.25f * m;
        }
        lo->propagateScoresDown(offset2);
      } else {
        for(int c=0; c<channels; c++) {
          float &loMeritHi0 = lo->meritHiBuffer[c].get(offset2);
          float merit = meritHiBuffer[c].get(offset);   

          float &loMeritOffHi0 = lo->meritOffHiBuffer[c].get(offset2);
          float meritOff = meritOffHiBuffer[c].get(offset);   

          float &loMHi0 = lo->mHiBuffer[c].get(offset2);
          float m = mHiBuffer[c].get(offset);   

          if(bScorePropagated) {
            loMeritHi0 += 0.5f * merit;
            loMeritOffHi0 += 0.5f * meritOff;
            loMHi0 += 0.5f * m;
          } else {
            loMeritHi0 += 0.75f * merit;
            loMeritOffHi0 += 0.75f * meritOff;
            loMHi0 += 0.75f * m;
          }          
        }
        bScorePropagated = true;
      }
    } else {
      for(int c=0; c<channels; c++) {
        float &loMeritHi = lo->meritHiBuffer[c].get(offset);
        float merit = meritHiBuffer[c].get(offset);
        loMeritHi += merit;

        float &loMeritOffHi = lo->meritOffHiBuffer[c].get(offset);
        float meritOff = meritOffHiBuffer[c].get(offset);
        loMeritOffHi += meritOff;

        float &loMHi = lo->mHiBuffer[c].get(offset);
        float m = mHiBuffer[c].get(offset);
        loMHi += m;
      }
      lo->propagateScoresDown(offset);
    }
  }

}

// meritLo contains the sum of merits for bands below and including this band
void SMS :: propagateScores(long offset)
{
  //printf("up %d %lld %ld\n",band,scoretime,offset);
  if(lo) {
    if(res == 2) {
      long offset2 = offset>>1;
      if(offset & 1) {
        for(int c=0; c<channels; c++) {
          float &meritLo = meritLoBuffer[c].get(offset);
          float merit = lo->meritLoBuffer[c].get(offset2);
          meritLo += 0.5f * merit;
          merit = lo->meritLoBuffer[c].get(offset2 + 1);
          meritLo += 0.5f * merit;

          float &meritOffLo = meritOffLoBuffer[c].get(offset);
          float meritOff = lo->meritOffLoBuffer[c].get(offset2);
          meritOffLo += 0.5f * meritOff;
          meritOff = lo->meritOffLoBuffer[c].get(offset2 + 1);
          meritOffLo += 0.5f * meritOff;

          float &mLo = mLoBuffer[c].get(offset);
          float m = lo->mLoBuffer[c].get(offset2);
          mLo += 0.5f * m;
          m = lo->mLoBuffer[c].get(offset2 + 1);
          mLo += 0.5f * m;
        }
      } else {
        for(int c=0; c<channels; c++) {
          float &meritLo = meritLoBuffer[c].get(offset);
          float merit = lo->meritLoBuffer[c].get(offset2);
          meritLo += merit;

          float &meritOffLo = meritOffLoBuffer[c].get(offset);
          float meritOff = lo->meritOffLoBuffer[c].get(offset2);
          meritOffLo += meritOff;

          float &mLo = mLoBuffer[c].get(offset);
          float m = lo->mLoBuffer[c].get(offset2);
          mLo += m;
        }
      }
    } else {
      for(int c=0; c<channels; c++) {
        float &meritLo = meritLoBuffer[c].get(offset);
        float merit = lo->meritLoBuffer[c].get(offset);
        meritLo += merit; 

        float &meritOffLo = meritOffLoBuffer[c].get(offset);
        float meritOff = lo->meritOffLoBuffer[c].get(offset);
        meritOffLo += meritOff; 

        float &mLo = mLoBuffer[c].get(offset);
        float m = lo->mLoBuffer[c].get(offset);
        mLo += m; 
      }
    }
  }
  // propagate down
  if(!hi) {
    propagateScoresDown(offset);
  }
}

void SMS :: score(AnalysisTrack *t)
{
  AnalysisTrackPoint *tp = (AnalysisTrackPoint*)t->getTrackPoint(scoretime);
  if(scoretime+1<=t->last) {
    AnalysisTrackPoint *tp1 = (AnalysisTrackPoint*)t->getTrackPoint(scoretime+1);
    for(int c=0; c<channels; c++) {
      float dp = tp1->ph[c] - tp->ph[c];
      float dp0 = 0.5f*h1*(tp->f + tp1->f);
      float dph = canonPI(dp - dp0);
      tp1->dph[c] = dph;
    }
  }
  float costOn[2];
  float meritOn[2];
  float costOff[2];
  float meritOff[2];
 
  // Expected dt at onset
  float E0 = N2/8.0f;
  float E0Off = -E0;
  // Expected D(dt) at onset
  float E1 = -min(E0,h*0.75f);
  float E1Off = E1;
  // Expected running average before onset
  float Emean = N2/6.0f;
  float EmeanOff = -N2/16.0f;

  
  float W0 = 1.5f / square(E0);
  float W1 = 1.5f / square(E1);
  float Wmean = 1.5f / square(Emean);
  float WmeanOff = 1.5f / square(EmeanOff);
  float WJ = 0.4f;
  float WDbS = 0.3f;
  float WMeanJ = 0.5f;

  float wStereo = 0.2f;
  float meritStereo = 0.2f;
  float wStereoOff = 0.1f;
  float meritStereoOff = 0.1f;
  float s[2];

  for(int c=0; c<channels; c++) {
    float dt = tp->dt[c];
    s[c] = scoretime==t->first?1.0f:t->meanS[c].get();
    float j;                                    
    float js;
    if(scoretime==t->first) {
      j = 0.0f;
      js = 0.0f;
      t->meanS[c].push(0.0f);
    } else {
      j = dBApprox(t->meanM[c].get(),tp->m[c]);
      js = t->meanS[c].get();
      float dbs = dBApprox(tp->s[c],1.0f);
      t->meanS[c].push(dbs);
      if(t->meanM[c].get() > tp->m[c]) j = -j;
    }
    t->diffdt2[c].push(dt);
    t->meandt[c].push(dt);
    
    float e0 =  max(-1.0f, 1.0f - W0 * square(dt  - E0));
    tp->meritOn[c].e0 = e0;

    float e1 = max(0.0f, 1.0f - W1 * square(t->diffdt2[c].get() - E1));
    tp->meritOn[c].e1 = e1;

    float emean = max(0.0f, 1.0f - Wmean * square(t->meandt[c].get() - Emean));
    tp->meritOn[c].emean = emean;

    // 0th order
    costOn[c] = e0;
    // 1st order
    costOn[c] += e1;
    // mean
    costOn[c] += emean;


    tp->meritOn[c].jmean = max(-1.0f, min(4.0f,WMeanJ * t->meanJ[c].get()));
    tp->meritOn[c].dbs = max(-1.0f, min(2.0f, WDbS * js));
    tp->meritOn[c].j = max(-1.0f, min(3.0f, WJ * j));

    costOn[c] += tp->meritOn[c].jmean;
    costOn[c] += tp->meritOn[c].dbs;
    costOn[c] += tp->meritOn[c].j;

    // 0th order
    costOff[c] = max(-1.0f, 1.0f - W0 * square(dt - E0Off));
    // 1st order
    costOff[c] += max(-0.5f, 1.0f - W1 * square(t->diffdt2[c].get() - E1Off));
    // mean
    costOff[c] += max(0.0f, 1.0f - WmeanOff * square(t->meandt[c].get() - EmeanOff));
    
    //costOff[c] += max(-1.0f, min(4.0f, -WMeanJ * t->meanJ[c].get()));
    //costOff[c] += max(-1.0f, min(2.0f, -WDbS * js));
    //costOff[c] += max(-1.0f, min(3.0f, -WJ * j));

    t->meanM[c].push(tp->m[c]);
    t->meanJ[c].push(j);

  }

  if(channels == 1) {
    meritOn[0] = costOn[0];
    meritOff[0] = costOff[0];
    tp->meritOn[0].total = meritOn[0];
  } else {
    meritOn[0] = 0.8f * costOn[0] + 0.2f * costOn[1];
    meritOn[1] = 0.8f * costOn[1] + 0.2f * costOn[0];
    tp->meritOn[0].total = meritOn[0];
    tp->meritOn[1].total = meritOn[1];
    meritOff[0] = 0.8f * costOff[0] + 0.2f * costOff[1];
    meritOff[1] = 0.8f * costOff[1] + 0.2f* costOff[0];
  }

  for(int c=0; c<channels; c++) {
    //tp->meritOn[c] = meritOn[c];
    tp->meritOff[c] = meritOff[c];
    t->ss[c].push_back(s[c]);
    t->meritsOn[c].push_back(meritOn[c]);
    t->meritsOff[c].push_back(meritOff[c]);
    bandMeritOn[c] += tp->m[c] * meritOn[c];
    bandMeritOff[c] += tp->m[c] * meritOff[c];
    bandM[c] += tp->m[c];
  }
}

void SMS :: score()
{
  for(int c=0; c<channels; c++) {
    bandMeritOn[c] = 0.0f;
    bandMeritOff[c] = 0.0f;
    bandM[c] = 0.0f;
  }
#ifdef MULTITHREADED
  pthread_mutex_lock(&trackMutex);
#endif
  for(list<Track*>::iterator tt = renderTracks.begin(); 
      tt != renderTracks.end();) {
    AnalysisTrack *t = (AnalysisTrack*)(*tt);
    if(scoretime >= t->first && scoretime <= t->last) {
      score(t);
    }
    ++tt;
  }
#ifdef MULTITHREADED
  pthread_mutex_unlock(&trackMutex);
#endif  

  for(int c=0; c<channels; c++) {

    meritBuffer[c].write(bandMeritOn[c]);
    meritLoBuffer[c].write(bandMeritOn[c]);
    meritHiBuffer[c].write(bandMeritOn[c]);

    meritOffBuffer[c].write(bandMeritOff[c]);
    meritOffLoBuffer[c].write(bandMeritOff[c]);
    meritOffHiBuffer[c].write(bandMeritOff[c]);

    mBuffer[c].write(bandM[c]);
    mLoBuffer[c].write(bandM[c]);
    mHiBuffer[c].write(bandM[c]);
  }

  scoretime++;
}

void SMS :: add(grain *g1, grain *g2, grain *gT, Cache *cache)
{
  if(channels == 1) {
    float *mag1;
    if(band >= minTrialBand) {
      c2even(g1->x, x10[0], N);
      mag1 = (float*)malloc((Nover2+1)*sizeof(float));
      calcmags(mag1, x10[0]);
#ifdef MULTITHREADED
      pthread_mutex_lock(&magMutex[0]);
#endif
      magQueue[0].push(mag1);
#ifdef MULTITHREADED
      pthread_mutex_unlock(&magMutex[0]);
#endif
      if(cache) {
        float *mag1dup = (float*)malloc((Nover2+1)*sizeof(float));
        memcpy(mag1dup,mag1,(Nover2+1)*sizeof(float));
        cache->mag1Cache[0].push_back(mag1dup);
      }
    }
    c2even(g2->x, x2[0], N);
    calcmags(mag2[0],x2[0]);

    c2even(gT->x, xT[0], N);
    calcDT(dT[0],mag2[0],x2[0],xT[0]);
    //c2even(gT->x, xT[0], N);
    //calcDT(dT[0],mag1,x10[0],xT[0]);
    
    for(int k=0; k<=Nover2; k++) {
      magTot[k] = mag2[0][k];
    }
    if(cache) {
      float *mag2dup = (float*)malloc((Nover2+1)*sizeof(float));
      memcpy(mag2dup,mag2[0],(Nover2+1)*sizeof(float));
      cache->mag2Cache[0].push_back(mag2dup);
    }
  } else {
    if(band >= minTrialBand) {
      c2even(g1->x, x10[0], N);
      c2odd(g1->x, x10[1], N);

      float *mag1[2];
      for(int c=0; c<2; c++) {
        mag1[c] = (float*)malloc((Nover2+1)*sizeof(float));
        calcmags(mag1[c], x10[c]);
#ifdef MULTITHREADED
        pthread_mutex_lock(&magMutex[c]);
#endif
        magQueue[c].push(mag1[c]);
#ifdef MULTITHREADED
        pthread_mutex_unlock(&magMutex[c]);
#endif
        if(cache) {
          float *mag1dup = (float*)malloc((Nover2+1)*sizeof(float));
          memcpy(mag1dup,mag1[c],(Nover2+1)*sizeof(float));
          cache->mag1Cache[c].push_back(mag1dup);
        }
      }
    }
    c2even(g2->x, x2[0], N);
    c2odd(g2->x, x2[1], N);
    
    calcmags(mag2[0],x2[0]);
    calcmags(mag2[1],x2[1]);

    c2even(gT->x, xT[0], N);
    calcDT(dT[0],mag2[0],x2[0],xT[0]);
    //c2even(gT->x, xT[0], N);
    //calcDT(dT[0],mag1[0],x10[0],xT[0]);
    
    c2odd(gT->x, xT[1], N);
    calcDT(dT[1],mag2[1],x2[1],xT[1]);
    //c2even(gT->x, xT[1], N);
    //calcDT(dT[1],mag1[1],x10[1],xT[1]);

    for(int k=0; k<=Nover2; k++) {
      magTot[k] = mag2[0][k] + mag2[1][k];
    }

    if(cache) {
      for(int c=0; c<2; c++) {
        float *mag2dup = (float*)malloc((Nover2+1)*sizeof(float));
        memcpy(mag2dup,mag2[c],(Nover2+1)*sizeof(float));
        cache->mag2Cache[c].push_back(mag2dup);
      }
    }
  }

  float magmax = magTot[0];
  for(int k=1;k<=kEnd;k++) {
    if(magmax < magTot[k]) magmax = magTot[k];
  }
  for(int c=0; c<channels; c++) {
    for(int k=0; k<=Nover2; k++) {
      mak2[c][k] = mag2[c][k] * k;
    }
  }
  float peakmin = magmax * peakThresh;

  float xt2 = 1.0f;
  bool bTroughN1 = false;
  bool bTroughN2 = false;
  float x0 = 1.0f;
  float y0[2];
  float x1 = 0.0f;
  float y1[2];
  bool bX0 = !lo;
  bool bX1 = false;
  AnalysisTrackPoint *prev = NULL;

  for(int c=0; c<channels; c++) {
    y0[c] = mag2[c][0];
    y1[c] = 0.0f;
  }
  Slice *slice = new Slice(band,addtime);

  for(int k=1; k<=kEnd; k++) {
    if(magTot[k] > peakmin && magTot[k] > magTot[k-1] && magTot[k] >= magTot[k+1]) {
      if(k < kLo) {
        x0 = findExtremum(magTot,k);
        for(int c=0; c<channels; c++) { y0[c] = findMagnitude(mag2[c],x0); }
        bX0 = true;
      } else if(k > kHi) {
        if(!bX1) {
          x1 = findExtremum(magTot,k);
          for(int c=0; c<channels; c++) { y1[c] = findMagnitude(mag2[c],x1); }
          if(prev) {
            prev->x01 = x1;
            for(int c=0; c<channels; c++) { prev->y01[c] = y1[c]; }
          }
          bX1 = true;
        }
      } else {
        float x = findExtremum(magTot,k); 
        AnalysisTrackPoint *p = new AnalysisTrackPoint(slice,peak2N,x,mag2,dT,x2,N,band,channels);
        if(prev) {
          prev->pn = p;
          p->pp = prev;
        } else {
          slice->bottom = p;
        }
        slice->top = p;
        prev = p;
        p->xtn2 = (float)maxK;
        bTroughN1 = true;
        bTroughN2 = true;
         p->xtp2 = xt2;
        p->x01 = x0;
        for(int c=0; c<channels; c++) { prev->y01[c] = y0[c]; }
      }
    } else if(magTot[k] <= magTot[k-1] && magTot[k] <= magTot[k+1]) {
      xt2 = findExtremum(magTot,k);
      xt2 = max(1.0f,xt2);
      xt2 = min((float)kEnd,xt2);
      if(bTroughN2) {
        prev->xtn2 = xt2;
        bTroughN2 = false;
      }
    }
  }
  if(bTroughN2) {
    prev->xtn2 = (float)kEnd;
  }
  if(!bX1 && !hi) {
    x1 = (float)kEnd;
    for(int c=0; c<channels; c++) { y1[c] = magTot[kEnd]; }
    bX1 = true;
  }
  for(int c=0; c<channels; c++) {
    float *dec = dec2[c];
    memset(dec,0,(Nover2+1)*sizeof(float));
    float *dek = dek2[c];
    memset(dek,0,(Nover2+1)*sizeof(float));
    if(bX0 && prev) {
      int k1 = lrintf(x0);
      int ko1 = k1 > x0 ? -1 : 1;
      float kf1 = k1 > x0 ? k1 - x0 : x0 - k1; 
      int k3 = min(kEnd,k1+peakWidth2);
      for(int k=lrintf(slice->bottom->xtp2);k<=k3;k++) {
        float m = interp2(k-k1,ko1,kf1);
        dec[k] += m * y0[c];
        m *= x0;
        dek[k] += m * y0[c];
      }
    }
    if(bX1 && prev) {
      int k1 = lrintf(x1);
      int ko1 = k1 > x1 ? -1 : 1;
      float kf1 = k1 > x1 ? k1 - x1 : x1 - k1; 
      int k3 = lrintf(slice->top->xtn2);
      for(int k=max(0,k1-peakWidth2);k<=k3;k++) {
        float m = interp2(k-k1,ko1,kf1);
        dec[k] += m * y1[c];
        m *= x1;
        dek[k] += m * y1[c];
      }
    }
  }

  float m2max = calcEnergy(slice,mag2,dec2,false);
  calcEnergy(slice,mak2,dek2,true);

  float m2min = m2max * peakThresh;
  for(AnalysisTrackPoint *p = slice->bottom;
      p; ) {
    AnalysisTrackPoint *pn = p->pn;
    bool bKeep = false;
    for(int c=0; c<channels; c++) {
      if(p->m2[c] > m2min) {
        bKeep = true;
        //printf("%g %g\n",p->f,TWOPI * (p->f2[c]/p->m2[c]) / (N<<band));
        p->f = TWOPI * (p->f2[c]/p->m2[c]) / (N<<band);
        //if(p->pp && p->f < p->pp->f) abort();
      } else {
        if(p->m2[c] < 0) { p->m2[c] = 0; }
      }
    }
    p->mTot2 = p->m2[0] + p->m2[1];
    
    if(1 && ((band == 5 && addtime >= 2986 && addtime <= 2986)))
      //if(1 && ((band == 4 && addtime >= 360 && addtime <= 370) || (band == 5 && addtime >= 360 && addtime <= 370 ))) 
    {
      printf("%d %lld %g %g\n",band,addtime,p->f,sqrt(p->m2[0])/MScale,sqrt(p->m2[1])/MScale);
    }

    if(!bKeep) {
      delete p;
    }
    p = pn;
  }

#ifdef MULTITHREADED
  pthread_mutex_lock(&sliceMutex);
#endif
  sliceBuffer.write(slice);
#ifdef MULTITHREADED
  pthread_mutex_unlock(&sliceMutex);
#endif
  addtime++;
}

float SMS :: calcEnergy(Slice *slice, float **mag2, float **dec2, bool bK)
{
  for(AnalysisTrackPoint *p = slice->bottom;
      p;
      p = p->pn) {
    
    for(int c=0; c<channels; c++) {
      float *dec = dec2[c];
      float *mag = mag2[c];
      float m2 = 0.0f;
      int k1 = lrintf(p->x);
      int ko1 = k1 > p->x ? -1 : 1;
      float kf1 = k1 > p->x ? k1 - p->x : p->x - k1; 
      int k0 = lrintf(p->xtp2);
      float kf0 = (k0 > p->xtp2 ? k0 - p->xtp2 : p->xtp2 - k0);
      int k2 = lrintf(p->xtn2);
      float kf2 = (k2 > p->xtn2 ? k2 - p->xtn2 : p->xtn2 - k2);
      if(k0 < p->xtp2) {
        m2 += (mag[k0] + mag[k0+1]) * 0.5f * (1.0f - kf0) + 0.5f * mag[k0+1];
        int i = k0 - k1;
        float m = interp2(i,ko1,kf1) * p->y[c];
        if(bK) m *= k0;
        m = min(m,mag[k0]) * 0.5f * (1.0f + kf0);
        m2 += m;
        dec[k0] += m;
        m = interp2(i+1,ko1,kf1) * p->y[c];
        if(bK) m *= (k0+1);
        m = min(m,mag[k0+1]) * 0.5f * kf0;
        m2 += m;
        dec[k0+1] += m;
        m = interp2(i-1,ko1,kf1) * p->y[c];
        if(bK) m *= (k0-1);
        m = min(m,mag[k0-1]);
        m2 += m;
        dec[k0-1] += m;
      } else {
        m2 += (mag[k0] + mag[k0-1]) * 0.5f * kf0 + 0.5f * mag[k0] + mag[k0+1];
        int i = k0 - k1;
        float m = interp2(i,ko1,kf1) * p->y[c];
        if(bK) m *= k0;
        m = min(m,mag[k0]) * 0.5f * (1.0f - kf0);
        m2 += m;
        dec[k0] += m;
        m = interp2(i-1,ko1,kf1) * p->y[c];
        if(bK) m *= (k0-1);
        m = min(m,mag[k0-1]) * 0.5f * (2.0f - kf0);
        m2 += m;
        dec[k0-1] += m;
      }
      if(k2 < p->xtn2) {
        m2 += mag[k2-1] + 0.5f * mag[k2] + 0.5f * kf2 * (mag[k2] + mag[k2+1]);
        int i = k2 - k1;
        float m = interp2(i,ko1,kf1) * p->y[c];
        if(bK) m *= k2;
        m = min(m,mag[k2]) * 0.5f * (1.0f - kf2);
        m2 += m;
        dec[k2] += m;
        m = interp2(i+1,ko1,kf1) * p->y[c];
        if(bK) m *= (k2+1);
        m = min(m,mag[k2+1]) * 0.5f * (2.0f - kf2);
        m2 += m;
        dec[k2+1] += m;
      } else {
        m2 += (mag[k2-1] + mag[k2]) * (1.0f - kf2) * 0.5f + 0.5f * mag[k2-1];
        int i = k2 - k1;
        float m = interp2(i,ko1,kf1) * p->y[c];
        if(bK) m *= k2;
        m = min(m,mag[k2]) * 0.5f * (1.0f + kf2);
        m2 += m;
        dec[k2] += m;
        m = interp2(i-1,ko1,kf1) * p->y[c];
        if(bK) m *= (k2-1);
        m = min(m,mag[k2-1]) * 0.5f * kf2;
        m2 += m;
        dec[k2-1] += m;
        m = interp2(i+1,ko1,kf1) * p->y[c];
        if(bK) m *= (k2+1);
        m = min(m,mag[k2+1]);
        m2 += m;
        dec[k2+1] += m;
      }    
      for(int k=k0+2;k<k2-1;k++) {
        m2 += mag[k];
      }
      if(k0 + 1 == k2 - 1) {
        m2 -= mag[k0+1];
      }
      for(int k=max(0,k1-peakWidth2);k<k0-1;k++) {
        float m = interp2(k-k1,ko1,kf1) * p->y[c];
        if(bK) m *= k;
        m = min(m,mag[k]);
        m2 += m;
        dec[k] += m;
      }
      int k3 = min(kEnd,k1+peakWidth2);
      for(int k=k2+2;k<=k3;k++) {
        float m = interp2(k-k1,ko1,kf1) * p->y[c];
        if(bK) m *= k;
        m = min(m,mag[k]);
        m2 += m;
        dec[k] += m;
      }
      if(bK) p->f2[c] = m2;
      else p->m2[c] = m2;
    }
  }

  float m2max = 0.0f;
  for(AnalysisTrackPoint *p = slice->bottom;
      p;
      p = p->pn) {

    for(int c=0; c<channels; c++) {    
      float *mag = mag2[c];
      float *dec = dec2[c];
      
      int k1 = lrintf(p->x);
      int ko1 = k1 > p->x ? -1 : 1;
      float kf1 = k1 > p->x ? k1 - p->x : p->x - k1; 
      int k0 = lrintf(p->xtp2);
      float kf0 = (k0 > p->xtp2 ? k0 - p->xtp2 : p->xtp2 - k0);
      int k2 = lrintf(p->xtn2);
      float kf2 = (k2 > p->xtn2 ? k2 - p->xtn2 : p->xtn2 - k2);
      float mdec = 0.0f;
      if(k0 < p->xtp2) {
        mdec += (dec[k0] + dec[k0+1]) * 0.5f * (1.0f - kf0) + 0.5f * dec[k0+1];
      } else {
        mdec += (dec[k0] + dec[k0-1]) * 0.5f * kf0 + 0.5f * dec[k0] + dec[k0+1];
      }
      if(k2 < p->xtn2) {
        mdec += dec[k2-1] + 0.5f * dec[k2] + 0.5f * kf2 * (dec[k2] + dec[k2+1]);
      } else {
        mdec += (dec[k2-1] + dec[k2]) * (1.0f - kf2) * 0.5f + 0.5f * dec[k2-1];
      }
      for(int k=k0+2;k<k2-1;k++) {
        mdec += dec[k];
      }
      if(k0 + 1 == k2 - 1) {
        mdec -= dec[k0+1];
      }
      if(bK) {
        p->f2[c] -= mdec;
        p->f2[c] *= mNorm;
      } else {
        p->m2[c] -= mdec;
        p->m2[c] *= mNorm;
        if(p->m2[c] > m2max) {
          m2max = p->m2[c];
        }
      }
    }
  }
  return m2max;
}

void SMS :: prepad(audio *buf, long n)
{
  if(band >= minTrialBand) {
    trialGrainBuf->write(buf,n);
  }
}

int SMS :: getTrialLatency()
{ 
  return minTrackSize;
}

float SMS :: interp2(int k, int ko, float kf) 
{
  return (1.0f-kf)*peak2N[k] + kf*peak2N[k+ko];
}

float SMS :: findExtremum(float *mag, int k)
{
  float y0 = mag[k-1];
  float y1 = mag[k];
  float y2 = mag[k+1];
  float d = (y0 + y2 - y1 - y1);
  float x = (d==0.0f?k:k + 0.5f * (y0 - y2) / d);
  return x;
}

float SMS :: findMagnitude(float *mag, float x)
{
  int ki = lrintf(x);
  float kf = ki<x?x-ki:ki-x;
  int ki1 = ki<x?ki+1:ki-1;
  return ((1.0f-kf)*mag[ki] + kf*mag[ki1]);
}

void SMS :: calcmags(float *mag, audio *x) {
  for(int k=0;k<=Nover2;k++) {
    mag[k] = norm2(x[k]);
  }
}

void SMS :: calcDT(float *dt, float *mag2, audio *x2, audio *xT)
{
  for(int k=0;k<=Nover2;k++) {
    dt[k] = mag2[k] == 0?0:((x2[k][0] * xT[k][0] + x2[k][1] * xT[k][1]) / mag2[k]);
  }
}

// From sbsmsampler
long SMS :: synthInit(int n, bool bBackwards, float stretch)
{
  if(stretch == 0.0f) return n;
  if(bBackwards) stretch = -stretch;
  bool bBackwardsInGrain = (stretch < 0.0f);
  if(bBackwards != bBackwardsInGrain) {
    n = max(0,min(n,(int)lrint(((samplePos - leftPos) * fabsf(stretch)))));
  } else {
    n = max(0,min(n,(int)lrint(((rightPos - samplePos) * fabsf(stretch)))));
  }
  if(bBackwardsInGrain) {
    int n0 = (int)lrint(ceil(synthOffset * fabsf(stretch)));
    n = min(n, max(0, n0));
  } else {
    int n0 = (int)lrint(ceil((h1 - synthOffset) * fabsf(stretch)));
    n = min(n, max(0, n0));
  }

  return n;
}

int SMS :: synthTracks(int n, SBSMSRenderer *renderer, bool bBackwards, float stretch, float pitch, Debugger *dbg)
{
#ifdef MULTITHREADED
  pthread_mutex_lock(&trackMutex);
#endif
  bool bBackwardsInGrain = (bBackwards != (stretch < 0.0f));
  int nInGrain = bBackwardsInGrain?-n:n;
  float stretchAbs = fabsf(stretch);
  float h2 = (stretch == 0.0f?0.0f:(float)(stretch * h1));
  
  SBSMSRenderChunk chunk;
  chunk.band = band;
  chunk.samplePos = samplePos;
  chunk.time = synthtime;
  chunk.stepSize = h1;
  chunk.chunkSize = h2;
  chunk.chunkPos = grainPos;
  chunk.length = nInGrain;
  chunk.f0 = pitch;
  chunk.f1 = pitch;

  renderer->startTime(chunk);

  for(list<Track*>::iterator tt=renderTracks.begin(); tt != renderTracks.end();) {
    Track *t = (*tt);      
    if(synthtime >= t->getStart()) {
      renderer->render(chunk,t,dbg);
      tt++;
    } else {
      break;
    }
  }
#ifdef MULTITHREADED
  pthread_mutex_unlock(&trackMutex);
#endif
  renderer->endTime(chunk);

  if(stretch != 0.0f) {
    synthOffset += (double)nInGrain / stretchAbs;
    grainPos += (float)nInGrain / (float)(stretchAbs * grainLength);
    if(M == 1) {
      double n1 = (double)nInGrain / stretch;
      samplePosCum += n1;
      int n1i = lrint(samplePosCum);
      samplePosCum -= (double)n1i;
      samplePos += n1i;
      if(samplePos < leftPos) samplePos = leftPos;
      else if(samplePos > rightPos) samplePos = rightPos;
    }
  }

  if(n && synthOffset <= 0.0) {
    return -1;
  } else if(n && synthOffset >= h1) {
    return 1;
  } else {
    return 0;
  }
}

void SMS :: stepSynth(int n)
{
  addtime += n;
  if(n == 0) {
    grainPos = 0.0f;
    grainLength = h1 + synthOffset;
    pruneTracks();
    assigntime--;
  } else {
    synthOffset -= h1;
    grainPos = 0.0f;
    grainLength = h1 - synthOffset;
    synthtime++;
    for(list<Track*>::iterator tt=renderTracks.begin(); 
        tt != renderTracks.end();) {
      Track *t = (*tt);      
      if(t->bEnded && synthtime >= t->getEnd()) {
        list<Track*>::iterator eraseMe = tt;
        tt++;
        renderTracks.erase(eraseMe);
        returnTrackIndex(t);
        delete t;
      } else {
        tt++;
      }
    }
  }
}

void SMS :: pruneTracks()
{
#ifdef MULTITHREADED
  pthread_mutex_lock(&trackMutex);
#endif
    for(list<Track*>::iterator tt=renderTracks.begin(); 
        tt != renderTracks.end();) {
      Track *t = (*tt);      
    if(synthtime <= t->getStart()) {
      list<Track*>::iterator eraseMe = tt;
      tt++;
      renderTracks.erase(eraseMe);
      returnTrackIndex(t);
      delete t;
    } else {
      t->prune(synthtime);
      tt++;
    }
  }
#ifdef MULTITHREADED
  pthread_mutex_unlock(&trackMutex);
#endif
}

SampleCountType SMS :: getSamplePos()
{
  return samplePos;
}

bool SMS :: isPastLeft()
{
  return samplePos <= leftPos;
}

bool SMS :: isPastRight()
{
  return samplePos >= rightPos;
}

void SMS :: setLeftPos(SampleCountType pos)
{
  leftPos = pos;
}

void SMS :: setRightPos(SampleCountType pos)
{
  rightPos = pos;
}

}
