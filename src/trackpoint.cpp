#include <math.h>
#include "real.h"
#include "utils.h"
#include "trackpoint.h"
#include "track.h"

namespace _sbsms_ {


void AnalysisTrackPoint :: destroy()
{
  refCount--;
  if(refCount <= 0) {
    delete this;
  }
}

AnalysisTrackPoint :: AnalysisTrackPoint(Slice *slice, float *peak, float x, float **mag2, float **dT, audio **X2, int N, int band, int channels)
{
  flags = 0;
  s[0] = 1.0f;
  s[1] = 1.0f;
  dph[0] = 0.0f;
  dph[1] = 0.0f;
  bScale[0]= false;
  bScale[1]= false;
  cache = NULL;
  bOnset[0] = false;
  bOnset[1] = false;
  bOffset[0] = false;
  bOffset[1] = false;
  
  refCount = 0;
  for(int d=0;d<3;d++) {
    dup[d] = NULL; 
  }
  for(int c=0; c<2; c++) { 
    y01[c] = 0.0f;
    m2[c] = 0.0f;
  }
  this->channels = channels;
  pp = NULL;
  pn = NULL;
  bConnect = false;
  bConnected = false;
  bDelete = false;
  bOwned = false;
  bMarked = false;
  bSplit = false;
  bMerge = false;
  owner = NULL;
  this->slice = slice;
  this->peak = peak;
  int ki = lrintf(x);
  int ki1;
  float kf;
  this->x = x;
  if(ki<x) {
    ki1 = ki + 1;
    kf = x - ki;
  } else {
    ki1 = ki - 1;
    kf = ki - x;
  }

  f = TWOPI*x/(float)(N*(1<<band));
  for(int c=0; c<channels; c++) {
    y[c] = ((1.0f-kf)*mag2[c][ki] + kf*mag2[c][ki1]);
    dt[c] = ((1.0f-kf)*dT[c][ki] + kf*dT[c][ki1]);
    audio *gx = X2[c];
    float norm0 = square(gx[ki][0]) + square(gx[ki][1]); 
    float ph0;
    if(norm0 > 0.0f) {
      ph0 = atan2(gx[ki][1],gx[ki][0]);
    } else {
      ph0 = 0.0f;
    }
    float ph1;
    float norm1 = square(gx[ki1][0]) + square(gx[ki1][1]);
    if(norm1 > 0.0f) {
      ph1 = atan2(gx[ki1][1],gx[ki1][0]);
    } else { 
      ph1 = 0.0f;
    }
    ph0 += (float)(ki&1)*PI;
    ph1 += (float)(ki1&1)*PI;
    if(kf < 0.5f) {
      ph1 = ph0 + canonPI(ph1 - ph0);
    } else {
      ph0 = ph1 + canonPI(ph0 - ph1);
    }
    ph[c] = canon2PI((1.0f-kf)*ph0 + kf*ph1);
    phSynth[c] = ph[c];
  }
}

bool bBlob = false;
void AnalysisTrackPoint :: absorb()
{
  if(bBlob) printf("absorb! %g\n",f);
  if(pp && pn) {
    if(x - pp->x < pn->x - x) {
      if((x - pp->x) / x < 0.05f) {
        for(int c=0; c<channels; c++) pp->m2[c] += m2[c];
        if(bBlob) printf("absorb 0 %g %g %g %g\n",f,pp->f,sqrt(m2[0])/MScale,sqrt(pp->m2[0])/MScale);
      }
    } else {
      if((pn->x - x) / x < 0.05f) {
        for(int c=0; c<channels; c++) pn->m2[c] += m2[c];
        if(bBlob) printf("absorb 1 %g %g %g %g\n",f,pn->f,sqrt(m2[0])/MScale,sqrt(pn->m2[0])/MScale);
      }
    }
  } else if(pp) {
    if((x - pp->x) / x < 0.05f) {
      for(int c=0; c<channels; c++) pp->m2[c] += m2[c];
      if(bBlob) printf("absorb 2 %g %g %g %g\n",f,pp->f,sqrt(m2[0])/MScale,sqrt(pp->m2[0])/MScale);
    }
  } else if(pn) {
    if((pn->x - x) / x < 0.05f) {
      for(int c=0; c<channels; c++) pn->m2[c] += m2[c];
      if(bBlob) printf("absorb 3 %g %g %g %g\n",f,pn->f,sqrt(m2[0])/MScale,sqrt(pn->m2[0])/MScale);
    }
  }
}

void AnalysisTrackPoint :: updateM()
{
  for(int c=0; c<channels; c++) {
    m[c] = (m2[c]>0.0f?sqrt(m2[c]):0.0f);
  }
}

AnalysisTrackPoint :: ~AnalysisTrackPoint()
{
  for(int d=0;d<3;d++) {
    if(dup[d]) {
      dup[d]->dup[2-d] = NULL;
    }
  }
  if(slice) slice->remove(this);
  if(pp && pn) {
    pp->pn = pn;
    pn->pp = pp;
  } else if(pp) {
    pp->pn = NULL;  
  } else if(pn) {
    pn->pp = NULL;
  }
}


Slice :: Slice(int band, const TimeType &time)
{
  this->band = band;
  this->time = time;
  bottom = NULL;
  top = NULL;
}

void Slice :: remove(AnalysisTrackPoint *tp)
{
  if(tp == top) {
    top = top->pp;
  }
  if(tp == bottom) {
    bottom = bottom->pn;
  }
}

Slice :: ~Slice()
{
  for(AnalysisTrackPoint *tp = bottom;
      tp;
      tp = tp->pn) {
    tp->slice = NULL;
  }
}

}
