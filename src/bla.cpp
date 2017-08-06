
    // if the start of a new track has a large forward time offset (dt)
    // then delay the track onse
    if(bStart && t->bStitchStart) {
      t->state = tp0->state;
    }
    if(tp1) {
      if(bStart && t->tailStart) {
        tp1->dt[c] / tp1->dtmax;
        if(tp1->dt[c] > dtmax) {
        t->bDelayed[c] = true;
      }
    } else if(t->bDelayed[c]) {
      // if the track onset has been delayed
      // then continue delaying while the dt remains large
      if(tp1 && tp1->dt[c] > dtmax) {
        //printf("%lld %d %g %g\n",i.time,i.band,tp1->f,tp1->dt[c]);
      } 
      // when the onset occurs, reset the phase
      else {
        gen->ph[c] = tp0->ph[c];
        t->bDelayed[c] = false;
      }
    } else {
      // if the track onset has occured
      if(tp1) {
        // but the track magnitude has consistently dropped off
        // then commence checking for a new onset
        if(tp1->m[c] < 0.2f * t->average1[c].get() && 
           t->average0[c].get() < 0.8f * t->average1[c].get()) {
          t->bDelayed[c] = true;
        } 
        // or if the dt is large again
        else if(tp1->dt[c] > dtmax) {
          // the previous track onset was spurious
          t->bDelayed[c] = true;
        }
      }
    }
    // When delaying a track onset, just set the trackpoint magnitude to 0
    // continue synthesizing both channels in case only one channel is delayed
    // 
    if(t->bDelayed[c]) {
      tp1->bDelayed[c] = true;
      m1[c] = 0.0f;
      tp1->m[c] = 0.0f;
    }
  }









      TrackPoint *tp = mintp;
      Cache *cache = sbsms->getCache(minband);
      mintime--;

      while(mintime >= 0) {
        long index = cache->indexAtTime[channel][mintime];
        int nTracks = cache->nTracksAtTime[channel][mintime];
        long trackIndex = cache->trackIndex[channel][mintime];
        TrackPoint *prev = NULL;
        for(int k = 0; k < nTracks; k++) {
          TrackIndexType trackIndex = cache->trackIndex[channel][index+k];
          if(trackIndex == minIndex) {
            prev = cache->trackPoints[channel][index+k];
            break;
          }
        }
        if(prev) {
          tp = prev;
        } else {
          break;
        }
        mintime--;
      }

    int top = bands-1;
    int bot = 0;
    for(int band=0; band<bands; band++) {
      if(view->endF >= topF[band]) {
        top = band;
        break;
      }
    }
    for(int band=bands-1; band>=0; band--) {
      if(view->startF <= botF[band]) {
        bot = band;
        break;
      }
    }

    printf("%d %d %g %g\n",top,bot,view->endF,view->startF);
    if(view->startF == botF[bot] && view->endF == topF[top]) {
      if(top == bot) {
        view->endF = view->startF + (view->endF - view->startF) / 2.0;
      } else {
        view->endF = topF[top+1];
      }
    } else {
      if(top == bot) {
        rangeY /= 4;
        double newStartF = view->startF + rangeY;
        view->startF = min(maxF,newStartF);
        double newEndF = view->endF - rangeY;
        view->endF = max(0.0,newEndF);
      } else {
        view->endF = topF[top];
        view->startF = botF[bot];
      }
    }

class ZoomScrollBar : public ScrollBar
{
public:
  int lastMouseTransPos;
  int dragStartMouseTransPos;
  double dragStartLength;

  ZoomScrollBar(bool bVertical) : ScrollBar(bVertical) {}

  void mouseDown(const MouseEvent& e)
  {
    ScrollBar::mouseDown(e);
    lastMouseTransPos = isVertical() ? e.x : -e.y;
    dragStartMouseTransPos = lastMouseTransPos;
    dragStartLength = getCurrentRange().getLength();
  }


  void mouseDrag (const MouseEvent& e)
  {
    ScrollBar::mouseDrag(e);
    const int mouseTransPos = isVertical() ? e.x : -e.y;
    if(lastMouseTransPos != mouseTransPos) {
      const int deltaPixels = mouseTransPos - dragStartMouseTransPos;
      double midRange = getCurrentRange().getStart() + 0.5 * getCurrentRange().getLength();
      double newLength = pow(1.25,((double)deltaPixels)/getHeight()) * dragStartLength;
      //setCurrentRange(midRange-0.5*newLength,newLength);
      setCurrentRange(getCurrentRange().getStart(),newLength);
    }
    lastMouseTransPos = mouseTransPos;
  }
};

void Track :: synth(float *out,
                    const TimeType &time,
                    int n,
                    float f0,
                    float f1,
                    int c)
{

  float m0, m1;
  float w0, w1;
  float ph0, ph1;
  bool bTailStart;
  bool bTailEnd;
  if(time >= end) return;
  if(time < last) {
    TrackPoint *tp1 = getTrackPoint(time+1);
    w1 = tp1->f * f1;
    m1 = tp1->m;
    ph1 = tp1->ph;
    if(bMerge && time + 1 == last) {
      m1 = 0.0f;
    }
    bTailStart = false;
    bTailEnd = false;
  } else {
    bTailStart = false;
    bTailEnd = (last != end);
  }
  if(time >= first) {
    TrackPoint *tp0 = getTrackPoint(time);
    w0 = tp0->f * f0;
    m0 = tp0->m;
    ph0 = tp0->ph;
    if(bSplit && time == first) {
      m0 = 0.0f;
    }
  } else {
    bTailStart = true;
  }

  if(bTailEnd) {
    int fall = min(n,w0==0.0f?384:min(384,(int)lrintf(PI * 4.0f / w0)));
    float dm = m0 / (float)fall;
    float w = w0;
    float *out2 = out;
    float *end = out + fall;
    long iph = lrintf(ph * WScale);
    if(iph>=W2PI) iph -= W2PI;
    long iw = lrintf(w * WScale);
    while(out2 != end) {
      if(iw < WPI) {
        long f = (iph>>PhShift)&Ph1;
        long i = iph>>WShift;
        *out2 += m0 * (float)(synthTable1[i] + f * synthTable2[i]);
      }
      out2++;
      m0 -= dm;
      iph += iw;
      iph &= W2PIMask;
    }
  }

  if(bTailStart) {
    int rise = min(n,w1==0.0f?384:min(384,(int)lrintf(PI * 3.0f / w1)));
    float dm = m1 / (float)rise;
    float w = w1;
    out += n;
    float *end = out-rise;
    long iph = lrintf(ph1 * WScale);
    iph &= W2PIMask;
    long iw = lrintf(w * WScale);
    while(out != end) {
      out--;
      m1 -= dm;
      iph -= iw;
      if(iph<0) iph += W2PI;
      if(iw < WPI) {
        long f = (iph>>PhShift)&Ph1;
        long i = iph>>WShift;
        *out += m1 * (float)(synthTable1[i] + f * synthTable2[i]);
      }
    }
  }

  if(!(bTailStart || bTailEnd)) {
    float dw = (w1 - w0) / n;
    float w = w0 + 0.5f * dw;
    float dm = (m1 - m0) / n;
    long iph = lrintf(ph0 * WScale);
    if(iph>=W2PI) iph -= W2PI;
    long iw = lrintf(w * WScale);
    long idw = lrintf(dw * WScale);

    float *end = out + n;
    while(out != end) {
      if(iw < WPI) {
        long f = (iph>>PhShift)&Ph1;
        long i = iph>>WShift;
        *out += m0 * (float)(synthTable1[i] + f * synthTable2[i]);
      }
      iph += iw;
      iw += idw;
      iph &= W2PIMask;
      m0 += dm;
      out++;
    }
  }
}
