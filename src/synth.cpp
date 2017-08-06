

void VoiceSynthesizer :: synth(int c, const SBSMSRenderChunk &i, float *in, float *out, Track *t)

=0;synth(int c,
                               float *in,
                               float *out,
                               countType synthtime,
                               float h2,
                               float offset,
                               int n,
                               float fScale0,
                               float fScale1,
                               Track *t)
{
  float m0, m1;
  float w0, w1;
  float ph0, ph1;
  bool bTailStart;
  bool bTailEnd;
  Gen *gen;
  if(time < last) {
    TrackPoint *tp1 = getTrackPoint(time+1);
    w1 = tp1->f * f1;
    m1 = tp1->m;
    ph1 = tp1->ph;
    if(time + 1 == last) {
      
    }
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
    gen = t
    w0 = tp0->f * f0;
    m0 = tp0->m;
    ph0 = tp0->ph;
    if(bSplit && time == first) {
      m0 = 0.0f;
    }
  } else {
    bTailStart = true;
  }


  if(!tp0 || !tp1) return false;

  bool bStart = (synthtime == t->start);
  bool bEnd = (synthtime+1 == t->end);

  GenParams params;
  SampleSynthesizer *s = sampleSynth;
  bool bNewGen = false;

  if(!t->gen || synthMode != s->synthMode) {
    if(t->gen) delete t->gen;
    bNewGen = true;
    synthMode = s->synthMode;
    if(synthMode == SynthModeOsc) {
      t->gen = new Osc(c);
    } else if(synthMode == SynthModeFilter) {
      t->gen = new Bandpass(c);
    } else if(synthMode == SynthModeDelay) {
      t->gen = new Delay(c);
    } else if(synthMode == SynthModeDecimate) {
      t->gen = new Decimator(c);
    } else if(synthMode == SynthModeGranulate) {
      t->gen = new Granulator(c,granuGrain);
    } else if(synthMode == SynthModeDWGS) {
      t->gen = new DWGS(c);
    }
  }
  
  if(offset == 0.0f) {
    t->w00 = tp0->f;
    t->w11 = tp1->f;
    
    if(bStart) {
      Track *precursor = t->getPrecursor();
      if(precursor) {
        t->state = precursor->stateDescendant;
      }
    }
    t->w00 *= t->mScale;
    t->w11 *= t->mScale;
  }

  Gen *gen = t->gen;

  float m0 = t->mScale * tp0->y;
  float m1 = t->mScale * tp1->y;

  bool bThresh;
  if(m0 >= s->mThresh || m1 >= s->mThresh || gen->m >= s->mThresh) {
    bThresh = true;
  } else {
    gen->m = 0.0f;
    bThresh = false;
  }

  float w0 = t->w00 * fScale0 * gen->x;
  float w1 = t->w11 * fScale1 * gen->x;

  float nf = (float)abs(n);

  
  float &pAM = gen->pAM;
  float &pFM = gen->pFM;

  if(bStart) {
    pAM = 1.0f;
    pFM = 0.0f;
  }

  float w01 = 0.5f * (w0 + w1);
  float w2 = (log(w01) - s->modPivot);
  float sm = getSidebandMod(w01,c);  
  
  float dpAM;
  float dpFM;
  if(s->fAM1Mode == 0) {
    dpAM = nf * s->fAM0 * exp(s->fAM1 * w2);
  } else {
    dpAM = nf * s->fAM0 * exp(s->fAM1 * sm);
  }
  if(s->fAM1Mode == 0) {
    dpFM = nf * s->fFM0 * exp(s->fFM1 * w2);
  } else {
    dpFM = nf * s->fFM0 * exp(s->fFM1 * sm);
  }

  float tAM0 = triangle(pAM);
  pAM += dpAM;
  while(pAM > 4.0f) pAM -= 4.0f;
  float tAM1 = triangle(pAM);

  float tFM0 = triangle(pFM);
  pFM += dpFM;
  while(pFM > 4.0f) pFM -= 4.0f;
  float tFM1 = triangle(pFM);

  if(bThresh) {
    float mm;
    if(s->mAM1Mode == 0) {
      mm = max(0.0f,(s->mAM0 * (1.0f + s->mAM1 * w2)));
    } else {
      mm = max(0.0f,(s->mAM0 * (1.0f + s->mAM1 * sm)));
    }
    m0 *= max(0.0f,tAM0 * mm + 1.0f);
    m1 *= max(0.0f,tAM1 * mm + 1.0f);    
    
    float se = (1.0f - s->sidebandEnv) + s->sidebandEnv * getSidebandEnv(w01,c);
    m0 *= se;
    m1 *= se;
    
    float dwFM;
    if(s->mFM1Mode == 0) {
      dwFM = max(-1.0f,min(1.0f, s->mFM1 * w2));
    } else {
      dwFM = max(-1.0f,min(1.0f, s->mFM1 * sm));
    }
    w0 = (w0 + tFM0 * (s->mFM0 + dwFM * w0));
    w1 = (w1 + tFM1 * (s->mFM0 + dwFM * w1));
    
    w0 = min(6.0f,max(0.001f,w0));
    w1 = min(6.0f,max(0.001f,w1));
    
    m0 = min(1.0f,max(0.0f,m0));
    m1 = min(1.0f,max(0.0f,m1));
    
    h2 = fabsf(h2);
    
    gen->state = t->state;
    
    float dm;
    
    //printf("%g %g\n",w01,sm);
    
    int dist;
    if(s->dist1Mode == 0) {
      dist = max(0,min((int)DistMax,(int)lrintf(s->dist0 * (1.0f + min(1.0f, s->dist1 * w2)))));
    } else {
      dist = max(0,min((int)DistMax,(int)lrintf(s->dist0 * (1.0f + min(1.0f, s->dist1 * sm)))));
    }
    params.tab1 = distSynthTable.getDistSynthTable1(dist);
    params.tab2 = distSynthTable.getDistSynthTable2(dist);
    
    if(s->Q1Mode == 0) {
      params.Q = s->Q0 * exp(4.0f * s->Q1 * w2);
    } else {
      params.Q = s->Q0 * exp(4.0f * s->Q1 * sm);
    }
    
    if(s->decBits1Mode == 0) {
      params.bits = min(32.0f,max(4.0f,s->decBits0 * (1.0f + s->decBits1 * w2)));
    } else {
      params.bits = min(32.0f,max(4.0f,s->decBits0 * (1.0f + s->decBits1 * sm)));
    }
    
    if(s->combFB1Mode == 0) {
      float F = 1.0f / (1.0f - s->combFB0);
      F = max(1.0f, F * expf(2.0f * s->combFB1 * w2));
      params.fb = F / (1.0f + F);
    } else {
      float F = 1.0f / (1.0f - s->combFB0);
      F = max(1.0f, F * expf(2.0f * s->combFB1 * sm));
      params.fb = F / (1.0f + F);
    }
    
    params.granMode = s->granMode;
    if(s->granRate1Mode == 0) {
      params.granRate = s->granRate0 * (1.0f + s->granRate1 * w2);
    } else {
      params.granRate = s->granRate0 * (1.0f + s->granRate1 * sm);
    }
    params.granRate = max(0.0f,min(1.0f,params.granRate));
    if(s->granSmooth1Mode == 0) {
      params.granSmooth = s->granSmooth0 * (1.0f + s->granSmooth1 * w2);
    } else {
      params.granSmooth = s->granSmooth0 * (1.0f + s->granSmooth1 * sm);
    }
    params.granSmooth = max(0.0f,min(1.0f,params.granSmooth));
    if(s->dwgsDecay1Mode == 0) {
      params.c1 = s->dwgsDecay0 * exp(2.0f * s->dwgsDecay1 * w2);
    } else {
      params.c1 = s->dwgsDecay0 * exp(2.0f * s->dwgsDecay1 * sm);
    }
    if(s->dwgsLopass1Mode == 0) {
      params.c3 = s->dwgsLopass0 * exp(2.0f * s->dwgsLopass1 * w2);
    } else {
      params.c3 = s->dwgsLopass0 * exp(2.0f * s->dwgsLopass1 * sm);
    }
    if(s->dwgsStringPos1Mode == 0) {
      params.inpos = max(0.0f,min(1.0f,s->dwgsStringPos0 * (1.0f + s->dwgsStringPos1 * w2)));
    } else {
      params.inpos = max(0.0f,min(1.0f,s->dwgsStringPos0 * (1.0f + s->dwgsStringPos1 * sm)));
    }
    
    gen->setParams(&params);
    
    if(gen->m < 0.0f) gen->m = 0.0f;
    //  if(t->index == 74) blob();
    
    //printf("%d %g %g %d %d %d %d\n",n,h2,offset,bStart,t->tailStart,bEnd,t->tailEnd);
    if(bStart && t->tailStart) {
      if(offset == 0.0f) {
        gen->w = w1;
      }
      if(n < 0) {
      int iStart;
      int rise;
      if(h2 == 0.0f) {
        iStart = 0;
        dm = 0.0f;
        rise = -n;
      } else {
        iStart = lrintf(h2 - h2 * offset);
        rise = min(max(0,(int)lrintf(t->rise - h2 + h2*offset)),max(-n,iStart));
        dm = gen->m / rise;
        rise = min(-n,rise);
      }
      gen->gen(0,-dm,in,out,rise);
      } else {
        int iStart;
        if(h2 == 0.0f) {
          iStart = 0;
          dm = 0.0f;
        } else {
          iStart = lrintf(h2 - h2 * offset);
          int rise = min(t->rise,max(n,iStart));
          dm = (m1 - gen->m) / rise;
          iStart = max(0,iStart-rise);
          if(synthMode == SynthModeOsc && offset == 0.0f) {
            gen->state.ph = canon(gen->state.ph - w0 * rise);
          }
        }
        if(iStart < n) {
          out += iStart;
          in += iStart;
          gen->gen(0,dm,in,out,n-iStart);
        }
      }
    } else if(bEnd && t->tailEnd) {
      if(offset == 0.0f) {
        gen->w = w0;
      }
      if(n < 0) {
        int iStart;
        if(h2 == 0.0f) {
          iStart = 0;
          dm = 0.0f;
        } else {
          iStart = lrintf(h2 * offset);
          int rise = min(t->fall,max(-n,iStart));
          dm = (m0 - gen->m) / rise;
          iStart = max(0,iStart-rise);
        }
        if(iStart < -n) {
          out += iStart;
          in += iStart;
          gen->gen(0,dm,in,out,-n-iStart);
        }
      } else {
        int fall;
        if(h2 == 0.0f) {
          dm = 0.0f;
          fall = n;
        } else {
          int nFrame = lrintf(h2 - h2 * offset);
          fall = min(max(0,(int)lrintf(t->fall-h2*offset)),max(n,nFrame));
          dm = gen->m / (float)fall;
          fall = min(n,fall);
        }
        gen->gen(0,-dm,in,out,fall);
      }
    } else {
      float dw;
      if(h2 == 0.0f) {
        dw = 0.0f;
        dm = 0.0f;
        gen->w = w0 + offset * (w1 - w0);
      } else {
        if(offset == 0.0f) {
          gen->w = w0;
        } else {
          gen->w = w0 + offset * (w1 - w0);
        }
        if(n < 0) {
          float dt = (float)max(-n,(int)lrintf(offset * h2));
          dm = (gen->m - m0) / dt;
          dw = (gen->w - w0) / dt;
          gen->w -= 0.5f * dw;
        } else {
          float dt = (float)max(n,(int)lrintf(h2 - offset * h2));
        dm = (m1 - gen->m) / dt;
        dw = (w1 - gen->w) / dt;
        gen->w += 0.5f * dw;
        }
      }
      //    if(t->index == 74) blob();
      
      if(n < 0) {
        gen->gen(-dw,-dm,in,out,-n);
      } else {
        gen->gen(dw,dm,in,out,n);
      }   
    } 
  }

  if(bEnd) {
    if(n >= 0 && offset == 0.0f) {
      Track *descendant = t->getDescendant();
      if(descendant && descendant->M < t->M) {
        t->stateDescendant = t->state;
      }
    }
    t->state = gen->state;
  } else if(bStart && t->tailStart) {
    t->state = gen->state;
  } else {
    t->state = gen->state;
    if(n >= 0) {
      Track *descendant = t->getDescendant();
      if(descendant && descendant->M > t->M) {
        t->stateDescendant = t->state;
      }
    }
  }
  return true;
}
