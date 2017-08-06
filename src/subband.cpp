#include "subband.h"
#include "real.h"
#include "sbsms.h"
#include "utils.h"
#include <algorithm>
using namespace std;

namespace _sbsms_ {

Cache *SubBand :: getCache(int band)
{
  if(band == 0) {
    return cache;
  } else {
    return sub->getCache(band-1);
  }
}

SubBand :: SubBand(SubBand *parent, int band, int channels, SBSMSQuality *quality, bool bAnalyze, bool bSynthesize, SubBand *source, TrackStitcher *stitcher, bool bCreateCache)
{
  if(band<quality->params.bands-1) {
    sub = new SubBand(this,band+1,channels,quality,bAnalyze,bSynthesize,source?source->sub:NULL,stitcher,bCreateCache);
  } else {
    sub = NULL;
  }
  this->quality = quality;
  this->channels = channels;
  this->parent = parent;
  this->band = band;
  this->M = (1<<band);
  int M_MAX = 1<<(quality->params.bands-1);
  this->N = quality->params.N[band];
  this->N1 = quality->params.N1[band];  
  this->N2 = quality->params.N2[band];
  this->res = quality->params.res[band];
  this->resMask = res - 1;
  this->bAnalyze = bAnalyze;
  this->bSynthesize = bSynthesize;
  nGrainsPerFrame = res;
  if(sub) nGrainsPerFrame *= sub->nGrainsPerFrame;
  inputFrameSize = M_MAX*quality->params.H;
  h = inputFrameSize / (M*nGrainsPerFrame);

  if(sub && bAnalyze) {
    samplesSubIn = new SampleBuf(NDownSample/2);
    grainsIn = new GrainBuf(NDownSample, NDownSample/SDownSample, NDownSample, hann);
    downSampledGrainAllocator = new GrainAllocator(NDownSample/2,NDownSample/2,hann);
  } else {
    samplesSubIn = NULL;
    grainsIn = NULL;
    downSampledGrainAllocator = NULL;
  }
  if(band >= minTrialBand) {
    grains[0] = new GrainBuf(N, h, N1, hannpoisson);
  } else {
    grains[0] = NULL;
  }
  grains[1] = new GrainBuf(N, h, N2, hannpoisson);
  // XXX
  grains[2] = new GrainBuf(N, h, N2, hannpoissonT);
  if(band >= minTrialBand) {
    analyzedGrains[0] = new GrainBuf(N, h, N1, hannpoisson);
  } else {
    analyzedGrains[0] = NULL;
  }
  analyzedGrains[1] = new GrainBuf(N, h, N2, hannpoisson);
  analyzedGrains[2] = new GrainBuf(N, h, N2, hannpoissonT);

#ifdef MULTITHREADED
  pthread_mutex_init(&dataMutex, NULL);
  for(int i=0; i<NumGrainTypes; i++) {
    pthread_mutex_init(&grainMutex[i], NULL);
  }
#endif
  sms = new SMS(sub?sub->sms:NULL,N,band,quality->params.bands-1,h,res,N1,N2,channels,analyzedGrains[1]->getWindowFFT());


  if(bCreateCache) {
    this->bCreateCache = true;
    this->cache = new Cache(channels);
    this->cacheRenderer = new CacheRenderer(cache,channels);
    cache->kLo = sms->kLo;
    cache->kHi = sms->kHi;
    cache->minK = sms->minK;
    cache->maxK = sms->maxK;
    cache->h = sms->h;
    renderers.push_back(cacheRenderer);
  } else {
    this->bCreateCache = false;
    this->cache = source?source->cache:NULL;
    this->cacheRenderer = NULL;
  }

  nTrialLatency = sms->getTrialLatency() / nGrainsPerFrame + 1;
  if(sms->getTrialLatency() % nGrainsPerFrame) nTrialLatency++;
  int nAdjustLatencyGrains = N1/(2*h);
  nAdjustLatency = nAdjustLatencyGrains / nGrainsPerFrame + 1;
  if(nAdjustLatencyGrains % nGrainsPerFrame) nAdjustLatency++;
  if(sub) nTrialLatency = max(nTrialLatency,sub->nTrialLatency);
  if(sub) nAdjustLatency = max(nAdjustLatency,sub->nAdjustLatency);

  nTrimLatency = sms->getTrimLatency() / nGrainsPerFrame + 1;
  if(sms->getTrialLatency() % nGrainsPerFrame) nTrimLatency++;
  if(sub) nTrimLatency = max(nTrimLatency,sub->nTrimLatency);

  nScoreLatency = 1;
  nPropagateLatency = 2;
  nMarkLatency = 1;
  nAssignLatency = 1;
  if(bAnalyze) {
    nRenderLatencyOriginal = sms->getRenderLatency();
  } else {
    nRenderLatencyOriginal = 1;
  }
  if(band==0) {
    SubBand *s = sub;
    while(s) {
      s->nTrialLatency = nTrialLatency;
      s->nAdjustLatency = nAdjustLatency;
      s->nTrimLatency = nTrimLatency;
      s = s->sub;
    }
  }

#ifdef MULTITHREADED
  nWriteSlack = 6;
  nAnalyzeSlack = 6;
  nExtractSlack = 6;
  nMarkSlack = 6;
  nAssignSlack = 6;
  nTrialSlack = 6;
  nAdjustSlack = 6;
  nScoreSlack = 6;
  nPropagateSlack = 6;
  nTrimSlack = 6;
  nRenderSlack = 6;
#else
  nWriteSlack = 2;
  nAnalyzeSlack = 2;
  nExtractSlack = 2;
  nMarkSlack = 2;
  nAssignSlack = 2;
  nAdjustSlack = 2;
  nTrialSlack = 2;
  nScoreSlack = 2;
  nPropagateSlack = 2;
  nTrimSlack = 2;
  nRenderSlack = 2;
#endif

  synthRenderer = NULL;
  outMixer = NULL;

  if(bSynthesize) {
    synthRenderer = new SynthRenderer(channels,M*h,stitcher);
    renderers.push_back(synthRenderer);
    if(sub) {
      samplesSubOut = new SampleBuf(0);
      outMixer = new Mixer(synthRenderer,samplesSubOut);
    } else {
      outMixer = synthRenderer;
    }
  }       
  init();
}

SubBand :: ~SubBand() 
{
  for(int i=0; i<NumGrainTypes; i++) {
    if(grains[i]) {
      delete grains[i];
    }
    if(analyzedGrains[i]) delete analyzedGrains[i];
  }
  delete sms;
  if(sub) {
    delete sub;
    if(bAnalyze) {
      delete grainsIn;
      delete samplesSubIn;
      delete downSampledGrainAllocator;
    }
    if(bSynthesize) {
      delete samplesSubOut;
      delete outMixer;
    }
  }
  if(cacheRenderer) {
    delete cacheRenderer;
    delete cache;
  }

  if(bSynthesize) delete synthRenderer;
}

void SubBand :: addRenderer(SBSMSRenderer *renderer)
{
  if(sub) sub->addRenderer(renderer);
  renderers.push_back(renderer);
}

void SubBand :: removeRenderer(SBSMSRenderer *renderer)
{
  if(sub) sub->removeRenderer(renderer);
  renderers.remove(renderer);
}

void SubBand :: setStretch(float stretch)
{
#ifdef MULTITHREADED
  pthread_mutex_lock(&dataMutex);
#endif
  if(!parent) {
    float oFrameSizef = (stretch==0.0f?1.0f:stretch)*(float)inputFrameSize;
    totalSizef += oFrameSizef;
    long oFrameSizei = lrintf(totalSizef);
    totalSizef -= oFrameSizei;
    outputFrameSize.write(oFrameSizei);
  }
  stretchRender.write(stretch);
#ifdef MULTITHREADED
  pthread_mutex_unlock(&dataMutex);
#endif
  if(sub) sub->setStretch(stretch);
}

void SubBand :: setPitch(float f)
{
  if(sub) sub->setPitch(f);
#ifdef MULTITHREADED
    pthread_mutex_lock(&dataMutex);
#endif
    pitchRender.write(f);
#ifdef MULTITHREADED
    pthread_mutex_unlock(&dataMutex);
#endif    
}

void SubBand :: stepAnalyzeFrame(int i)
{
  if(sub) sub->stepAnalyzeFrame(i);
  nFramesAnalyzed[i]++;
}

void SubBand :: stepExtractFrame()
{
  if(sub) sub->stepExtractFrame();
  nFramesExtracted++;
}

void SubBand :: stepMarkFrame()
{
  if(sub) sub->stepMarkFrame();
  nFramesMarked++;
}

void SubBand :: stepAssignFrame()
{
  if(sub) sub->stepAssignFrame();
  nFramesAssigned++;
}

void SubBand :: stepTrialFrame()
{
  if(sub) sub->stepTrialFrame();
  nFramesTrialed++;
}

void SubBand :: stepAdjustFrame()
{
  if(sub) sub->stepAdjustFrame();
  nFramesAdjusted++;
}

void SubBand :: stepScoreFrame()
{
  if(sub) sub->stepScoreFrame();
  nFramesScored++;
}

void SubBand :: stepPropagateFrame()
{
  if(sub) sub->stepPropagateFrame();
  nFramesPropagated++;
}

void SubBand :: stepTrimFrame()
{
  if(sub) sub->stepTrimFrame();
  nFramesTrimmed++;
}

void SubBand :: stepRenderFrame()
{
  if(sub) sub->stepRenderFrame();
#ifdef MULTITHREADED
  pthread_mutex_lock(&dataMutex);
#endif
  stretchRender.advance(1);
  pitchRender.advance(1);
#ifdef MULTITHREADED
  pthread_mutex_unlock(&dataMutex);
#endif
  nFramesRendered++;
}

void SubBand :: stepReadFrame()
{
  if(sub) sub->stepReadFrame();
  nFramesRead++;
  if(bCreateCache) cache->frames++;
}

bool SubBand :: writeInit()
{
  long n = getFramesAtFront(0);
  for(int i=1; i<NumGrainTypes; i++) {
    n = min(n,getFramesAtFront(i));
  }
  return (n <= nWriteSlack);
}

long SubBand :: readInit()
{
  long n = max(0L,min(1L,min(n,nFramesRendered-nFramesRead)));
  if(sub) n = min(n,sub->readInit());
  return n;
}

long SubBand :: analyzeInit(int i, bool bSet, long n)
{
  if(!parent) {
    n = getFramesAtFront(i);
    n = max(0L,min(1L,min(n,nAnalyzeSlack-(long)(nFramesAnalyzed[i]-nFramesExtracted))));
  }
  if(bSet) {
    nGrainsToAnalyze[i] = n * nGrainsPerFrame;
    if(sub) {
      sub->analyzeInit(i,bSet,n);
    }
  }
  return n;
}

long SubBand :: extractInit(bool bSet)
{
  long n;
  if(sub) n = res*sub->extractInit(bSet);
  if(!sub) {
    n = max(0L,min(1L,nExtractSlack+nMarkLatency-(long)(nFramesExtracted-nFramesMarked)));
    for(int i=0; i<NumGrainTypes; i++) {
      n = max(0L,min(1L,min(n,(long)(nFramesAnalyzed[i]-nFramesExtracted))));
    }
  }
  if(bSet) {
    nGrainsToExtract = n;
  }
  return n;
}

long SubBand :: markInit(bool bSet)
{
  long n;
  if(sub) n = res*sub->markInit(bSet);
  if(!sub) n = max(0L,min(1L,min((long)(nFramesExtracted-nFramesMarked)-nMarkLatency,
                                 nMarkSlack+nAssignLatency-(long)(nFramesMarked-nFramesAssigned))));
  if(bSet) {
    nGrainsToMark = n;
  }
  return n;
}

long SubBand :: assignInit(bool bSet)
{
  long n;
  if(sub) {
    n = res * sub->assignInit(bSet);
  } else {
    n = max(0L,min(1L,min((long)(nFramesMarked-nFramesAssigned)-nAssignLatency,
                          nAssignSlack+nTrialLatency-(long)(nFramesAssigned-nFramesTrialed))));
  }
  if(bSet) {
    nGrainsToAdvance = n;
    nGrainsToAssign = n;
    if(n) {
      if(nFramesAssigned==0) {
        sms->start(0);
      }      
    }
  }
  return n;
}

long SubBand :: trialInit(bool bSet)
{
  long n;
  if(sub) {
    n = res * sub->trialInit(bSet);
  } else {
    n = max(0L,min(1L,min((long)(nFramesAssigned-nFramesTrialed)-nTrialLatency,
                          nTrialSlack+nAdjustLatency-(long)(nFramesTrialed-nFramesAdjusted))));
  }
  if(bSet) {
    nGrainsToTrial = n;
    nGrainsTrialed = 0;
  }
  return n;
}

long SubBand :: adjustInit(bool bSet)
{
  long n;
  if(sub) {
    n = res * sub->adjustInit(bSet);
  } else {
    n = 1;
    n = min(n,(long)(nFramesTrialed-nFramesAdjusted-nAdjustLatency));
    n = min(n,nAdjustSlack+nScoreLatency-(long)(nFramesAdjusted-nFramesScored));
    n = max(0L,n);
  }
  if(bSet) {
    nGrainsToAdjust = n;
    nGrainsAdjusted = 0;
  }
  return n;
}


long SubBand :: scoreInit(bool bSet)
{
  long n;
  if(sub) {
    n = res * sub->scoreInit(bSet);
  } else {
    n = 1;
    n = min(n,(long)(nFramesAdjusted-nFramesScored-nScoreLatency));
    n = min(n,nScoreSlack+nPropagateLatency-(long)(nFramesScored-nFramesPropagated));
    n = max(0L,n);
  }
  if(bSet) {
    nGrainsToScore = n;
    nGrainsScored = 0;
  }
  return n;
}

long SubBand :: propagateInit(bool bSet)
{
  long n;
  if(sub) {
    n = res * sub->propagateInit(bSet);
  } else {
    n = 1;
    n = min(n,(long)(nFramesScored-nFramesPropagated-nPropagateLatency));
    n = min(n,nPropagateSlack+nTrimLatency-(long)(nFramesPropagated-nFramesTrimmed));
    n = max(0L,n);
  }
  if(bSet) {
    nGrainsToPropagate = n;
    if(n) {
      if(nFramesPropagated == 0) {
        sms->propagateScores(0);
        nGrainsPropagated++;
      }
    }
  }
  return n;
}

long SubBand :: trimInit(bool bSet)
{
  long n;
  if(sub) {
    n = res * sub->trimInit(bSet);
  } else {
    n = 1;
    n = min(n,(long)(nFramesPropagated-nFramesTrimmed-nTrimLatency));
    n = min(n,nTrimSlack+nRenderLatency-(long)(nFramesTrimmed-nFramesRendered));
    n = max(0L,n);
  }
  if(bSet) {
    nGrainsToTrim = n;
    nGrainsTrimmed = 0;
  }
  return n;
}

long SubBand :: renderInit(bool bSet) 
{
  long n;
  if(sub) {
    n = res * sub->renderInit(bSet);
  } else {
    n = max(0L,min(1L,min((long)(nFramesTrimmed-nFramesRendered)-nRenderLatency,
                          nRenderSlack-(long)(nFramesRendered-nFramesRead))));
  }
  if(bSet) {
    nGrainsRendered = 0;
    nGrainsToRender = n;
  }
  return n;
}

void SubBand :: analyze(int i)
{
  if(sub) sub->analyze(i);
  if(grains[i]) {
    vector<grain*> gV;
#ifdef MULTITHREADED
    pthread_mutex_lock(&grainMutex[i]);
#endif
    for(int k=grains[i]->readPos;k<grains[i]->readPos+nGrainsToAnalyze[i];k++) {
      grain *g = grains[i]->read(k);
      gV.push_back(g);
    }
#ifdef MULTITHREADED
    pthread_mutex_unlock(&grainMutex[i]);
#endif

    for(int k=0;k<nGrainsToAnalyze[i];k++) {
      gV[k]->analyze();
    }

#ifdef MULTITHREADED
    pthread_mutex_lock(&grainMutex[i]);
#endif
    for(int k=0;k<nGrainsToAnalyze[i];k++) {
      analyzedGrains[i]->write(gV[k]);
    }
    grains[i]->advance(nGrainsToAnalyze[i]);
#ifdef MULTITHREADED
    pthread_mutex_unlock(&grainMutex[i]);
#endif
  }
}

void SubBand :: extract()
{
  if(sub) sub->extract();
  vector<grain*> gV[NumGrainTypes];

  for(int i=0; i<NumGrainTypes; i++) {
    if(grains[i]) {
#ifdef MULTITHREADED
      pthread_mutex_lock(&grainMutex[i]);
#endif    
      for(int k=analyzedGrains[i]->readPos;
          k<analyzedGrains[i]->readPos+nGrainsToExtract;
          k++) {
        grain *g = analyzedGrains[i]->read(k);
        gV[i].push_back(g);
      }
#ifdef MULTITHREADED
      pthread_mutex_unlock(&grainMutex[i]);
#endif
    }
  }

  for(int k=0;k<nGrainsToExtract;k++) {
    grain *g1 = (grains[0]?gV[0][k]:NULL);
    grain *g2 = gV[1][k];
    grain *gT = gV[2][k];
    sms->add(g1,g2,gT,cache);
  }

  for(int i=0; i<NumGrainTypes; i++) {
    if(grains[i]) {
#ifdef MULTITHREADED
      pthread_mutex_lock(&grainMutex[i]);
#endif
      analyzedGrains[i]->advance(nGrainsToExtract);
#ifdef MULTITHREADED
      pthread_mutex_unlock(&grainMutex[i]);
#endif
    }
  }
}

void SubBand :: mark()
{
  long ntodo = parent?1:nGrainsToMark;
  long ndone = 0;
  while(ndone<ntodo) {  
    sms->mark(nGrainsMarked);
    if(nGrainsMarked&resMask || res==1) {
      if(sub) sub->mark();
    }
    ndone++;
    nGrainsMarked++;
  }
}

void SubBand :: assign() 
{
  for(int ndone=0; ndone<nGrainsToAssign; ndone++) {
    assignStart();
    bool bCont = true;
    while(bCont) {
      assignInit();
      assignFind();
      bCont = assignConnect();
    }
    assignStep();
    splitMerge();
  }
}

void SubBand :: assignStart()
{
  if(sub && !(nGrainsAssigned&resMask)) sub->assignStart();
  sms->assignStart(nGrainsAssigned);
}

void SubBand :: assignInit()
{
  if(sub) sub->assignInit();
  sms->assignInit(nGrainsAssigned);
}

void SubBand :: assignFind()
{
  if(sub) sub->assignFind();
  sms->assignFind(nGrainsAssigned);
}

bool SubBand :: assignConnect()
{
  bool bCont = false;
  if(sub) {
    if(sub->assignConnect()) {
      bCont = true;
    }
  }  
  if(sms->assignConnect(nGrainsAssigned,false)) {
    bCont = true;
  }
  return bCont;
}

void SubBand :: assignStep()
{
  sms->assignConnect(nGrainsAssigned,true);
  if(sub && !((nGrainsAssigned+1)&resMask)) {
    sub->assignStep();
  }
  sms->start(nGrainsAssigned+1);
}

void SubBand :: splitMerge()
{
  nGrainsAssigned++;
  if(sub && !(nGrainsAssigned&resMask)) {
    sub->splitMerge();
  }
  sms->splitMerge();
}

void SubBand :: advance()
{
  long ntodo = parent?1:nGrainsToAdvance;
  long ndone = 0;
  while(ndone<ntodo) {
    if(sub && !(nGrainsAdvanced&resMask)) {
      sub->advance();
    }
    sms->advance();
    nGrainsMarked--;
    nGrainsAssigned--;
    nGrainsAdvanced++;
    ndone++;
  }
}

void SubBand :: trial()
{
  for(int i=0; i<nGrainsToTrial; i++) {  
    trialStart();
    trialTrial();
    trialEnd();
  }
}

void SubBand :: trialStart()
{
  if(!(nGrainsTrialed&resMask)) {
    if(sub) sub->trialStart();
    sms->trialStart();
  }
}

void SubBand :: trialTrial()
{
  if(sub && !(nGrainsTrialed&resMask)) {
    sub->trialTrial();
  }
  sms->trial();
}

void SubBand :: trialEnd()
{
  nGrainsTrialed++;
  if(!(nGrainsTrialed&resMask)) {
    if(sub) sub->trialEnd();
    sms->trialEnd();
  }
}

void SubBand :: adjust()
{
  long ntodo = parent?1:nGrainsToAdjust;
  long ndone = 0;
  while(ndone<ntodo) {  
    if(!(nGrainsAdjusted&resMask)) {
      if(sub) sub->adjust();
    }
    sms->adjust(cache);
    ndone++;
    nGrainsAdjusted++;
  }
}

void SubBand :: score()
{
  long ntodo = parent?1:nGrainsToScore;
  long ndone = 0;
  while(ndone<ntodo) {  
    if(!(nGrainsScored&resMask)) {
      if(sub) sub->score();
    }
    sms->score();
    ndone++;
    nGrainsScored++;
  }
}

void SubBand :: propagate()
{
  long ntodo = parent?1:nGrainsToPropagate;
  long ndone = 0;
  while(ndone<ntodo) {  
    if(!((nGrainsPropagated+1)&resMask)) {
      if(sub) sub->propagate();
    }
    sms->propagateScores(nGrainsPropagated);
    ndone++;
    nGrainsPropagated++;
  }
}

void SubBand :: trim()
{
  long ntodo = parent?1:nGrainsToTrim;
  long ndone = 0;
  while(ndone<ntodo) {  
    if(!(nGrainsTrimmed&resMask)) {
      if(sub) sub->trim();
    }
    sms->trim(cache);
    ndone++;
    nGrainsPropagated--;
    nGrainsTrimmed++;
  }
}

void SubBand :: readSubSamples()
{
  if(sub) sub->readSubSamples();
  if(sub) {
    audio fromSub[subBufSize];
    long nFromSub = 0;
    do {
      nFromSub = sub->outMixer->read(fromSub,subBufSize);
      samplesSubOut->write(fromSub, nFromSub);
    } while(nFromSub>0);
  }
}

long SubBand :: read(audio *buf, long n) 
{
  long nRead = 0;
  long nToRead = n;
  readSubSamples();
  while(nToRead && nRead < n && outputFrameSize.nReadable()) {
    long nToReadFromOutputFrame = outputFrameSize.read();
    nToRead = min(n-nRead,nToReadFromOutputFrame-nReadFromOutputFrame);
    nToRead = outMixer->read(buf+nRead, nToRead);
    nReadFromOutputFrame += nToRead;
    nRead += nToRead;
    if(nReadFromOutputFrame == nToReadFromOutputFrame) {
      nReadFromOutputFrame = 0;
      outputFrameSize.advance(1);
      stepReadFrame();
    }
  }
  return nRead;
}


long SubBand :: renderSynchronous() 
{
  for(list<SBSMSRenderer*>::iterator i = renderers.begin(); i != renderers.end(); ++i) {
    SBSMSRenderer *renderer = *i;
    renderer->startFrame();
  }
  renderInit(true);
  render();
  stepRenderFrame();
  for(list<SBSMSRenderer*>::iterator i = renderers.begin(); i != renderers.end(); ++i) {
    SBSMSRenderer *renderer = *i;
    renderer->endFrame();
  }
  long samples = outputFrameSize.read();
  outputFrameSize.advance(1);
  stepReadFrame();
  return samples;
}

void SubBand :: render()
{

#ifdef MULTITHREADED
  pthread_mutex_lock(&dataMutex);
#endif
  float stretch = stretchRender.read();
  float f0 = pitchRender.read(pitchRender.readPos);
  float f1;
  if(pitchRender.nReadable()>=2) {
    f1 = pitchRender.read(pitchRender.readPos+1);
  } else {
    f1 = f0;
  }
#ifdef MULTITHREADED
  pthread_mutex_unlock(&dataMutex);
#endif
  long ntodo = parent?1:nGrainsToRender;
  long ndone = 0;
  long nRenderedTotal = 0;

  while(ndone<ntodo) {
    if(nGrainsRendered%res==0)
      if(sub) sub->render();
    float df = (f1-f0)/(float)nGrainsToRender;
    sms->render(stretch,f0+nGrainsRendered*df,f0+(nGrainsRendered+1)*df,renderers);
    int pos = samplePos%inputFrameSize;
    if(pos && !bSynthStarted) {
      nToDrop = lrintf(pos * stretch);
    }
    bSynthStarted = true;
    nGrainsRendered++;
    ndone++;
  }
}

void SubBand :: renderComplete(const SampleCountType &samples)
{
  for(list<SBSMSRenderer*>::iterator i = renderers.begin(); i != renderers.end(); ++i) {
    SBSMSRenderer *renderer = *i;
    renderer->end(samples);
  }
}

long SubBand :: write(audio *inBuf, long n, float stretch, float pitch)
{
  long nWritten = 0;

  if(cache) {    
    if(nToStart) {
      if(nToStart >= n) {
        nToStart -= min(nToStart, n);
      } else {
        cache->wave.write(inBuf+nToStart,n-nToStart);
        nToStart = 0;
      }
    } else {
      cache->wave.write(inBuf,n);
    }
  }

  while(nWritten<n) {
    long nToWrite = min(nToWriteForGrain,n-nWritten);
    if(nToDrop2) {
      nToWrite = min(nToDrop2,nToWrite);
      nToDrop2 -= nToWrite;
      nToDrop1 -= nToWrite;
    } else {
      if(nToDrop1) {
        nToWrite = min(nToDrop1,nToWrite);
        nToDrop1 -= nToWrite;
      } else {
        if(nToPrepad) {
          nToWrite = min(nToPrepad,nToWrite);
          sms->prepad(inBuf+nWritten, nToWrite);
          nToPrepad -= nToWrite;
        }
#ifdef MULTITHREADED
        pthread_mutex_lock(&grainMutex[0]);
#endif      
        if(grains[0]) {
          grains[0]->write(inBuf+nWritten, nToWrite);
        }
#ifdef MULTITHREADED
        pthread_mutex_unlock(&grainMutex[0]);
#endif
      }
#ifdef MULTITHREADED
      pthread_mutex_lock(&grainMutex[1]);
#endif      
      grains[1]->write(inBuf+nWritten, nToWrite);
#ifdef MULTITHREADED
      pthread_mutex_unlock(&grainMutex[1]);
#endif

#ifdef MULTITHREADED
      pthread_mutex_lock(&grainMutex[2]);
#endif      
      grains[2]->write(inBuf+nWritten, nToWrite);
#ifdef MULTITHREADED
      pthread_mutex_unlock(&grainMutex[2]);
#endif
    }
    nWritten += nToWrite;
    nToWriteForGrain -= nToWrite;
    if(nToWriteForGrain == 0) {
      nToWriteForGrain = h;
      if(!parent) {
        if(nGrainsWritten == 0) {
          setStretch(stretch);
          setPitch(pitch);
        }
        nGrainsWritten++;
        if(nGrainsWritten == nGrainsPerFrame) {
          nGrainsWritten = 0;
        }
      }
    }
  }

  if(sub) {
    grainsIn->write(inBuf, n);
    long nGrainsRead = 0;
    for(int k=grainsIn->readPos;k<grainsIn->writePos;k++) {
      grain *g = grainsIn->read(k); g->analyze();
      grain *gdown = downSampledGrainAllocator->create();
      g->downsample(gdown);
      samplesSubIn->write(gdown, hSub);
      downSampledGrainAllocator->forget(gdown);
      nGrainsRead++;
    }
    grainsIn->advance(nGrainsRead);  
    long nWriteToSub = samplesSubIn->nReadable();
    audio *subBuf = samplesSubIn->getReadBuf();
    nWriteToSub = sub->write(subBuf,nWriteToSub,stretch,pitch);
    samplesSubIn->advance(nWriteToSub);
  }
  return n;
}

bool SubBand :: isDone()
{
  if(!nFramesRead || nFramesRead < nFramesAssigned)
    return false;
  return true;
}

void SubBand :: process(bool bRender)
{
  for(int i=0; i<NumGrainTypes; i++) {
    if(analyzeInit(i,true)) {
      analyze(i);
      stepAnalyzeFrame(i);
    }
  }

  if(extractInit(true)) {
    extract();
    stepExtractFrame();
  }

  if(markInit(true)) {
    mark();
    stepMarkFrame();
  }
  
  if(assignInit(true)) {
    assign();
    advance();
    stepAssignFrame();
  }
  
  if(trialInit(true)) {
    trial();
    stepTrialFrame();
  }
  
  if(adjustInit(true)) {
    adjust();
    stepAdjustFrame();
  }

  if(scoreInit(true)) {
    score();
    stepScoreFrame();
  }

  if(propagateInit(true)) {
    propagate();
    stepPropagateFrame();
  }

  if(trimInit(true)) {
    trim();
    stepTrimFrame();
  }
  
  if(bRender) {
    if(renderInit(true)) {
      render();
      stepRenderFrame();
    }
  }
}

long SubBand :: getFramesAtFront(int i)
{
  long n = 65536;
#ifdef MULTITHREADED
  pthread_mutex_lock(&grainMutex[i]);
#endif
  if(grains[i]) {
    n = grains[i]->nReadable() / nGrainsPerFrame;
  }
#ifdef MULTITHREADED
  pthread_mutex_unlock(&grainMutex[i]);
#endif
  if(sub) n = min(n,sub->getFramesAtFront(i));
  return n;
}

long SubBand :: getInputFrameSize()
{
  return inputFrameSize;
}




/* From sbsmsampler */


void SubBand :: reset(bool flushInput) {

  if(sub) sub->reset(flushInput);
  init();
  outputFrameSize.clear();
  bSynthStarted = false;
  synthFramePos = 0;
  bSynthGrainStart = true;
  nGrainsSynthed = 0;
  stretchRender.clear();
  pitchRender.clear();
  for(int i=0; i<NumGrainTypes; i++) {
    if(grains[i]) {
      grains[i]->clear();
      analyzedGrains[i]->clear();
    }
  }
  if(sub) {
    if(bAnalyze) {
      grainsIn->clear();
      samplesSubIn->clear();
    }
    if(bSynthesize) {
      samplesSubOut->clear();
    }
  }
  sms->reset();
  if(synthRenderer) synthRenderer->reset(flushInput);
}


void SubBand :: init()
{
  nToDrop1 = (quality->getMaxPresamples()/M - N1/2);
  nToDrop2 = (quality->getMaxPresamples()/M - N2/2);
  nToStart = quality->getMaxPresamples()/M;
  nReadFromOutputFrame = 0;
  nToWriteForGrain = (quality->getMaxPresamples()/M + N2/2);
  nRenderLatency = nRenderLatencyOriginal;
  nToPrepad = N1/2;
  for(int i=0; i<NumGrainTypes; i++) {
    nFramesAnalyzed[i] = 0;
  }
  nFramesExtracted = 0;
  nFramesMarked = 0;
  nFramesAssigned = 0;
  nFramesTrialed = 0;
  nFramesRendered = 0;
  nGrainsMarked = 0;
  nGrainsAssigned = 0;
  nGrainsTrialed = 0;
  nGrainsAdvanced = 0;
  nGrainsWritten = 0;
  nFramesAdjusted = 0;
  nGrainsAdjusted = 0;
  nFramesScored = 0;
  nGrainsScored = 0;
  nFramesPropagated = 0;
  nGrainsPropagated = 0;
  nFramesTrimmed = 0;
  nGrainsTrimmed = 0;
  samplePos = 0;
  nToDrop = 0;
  nFramesRead = 0;
  totalSizef = 0.0;
  bWritingComplete = false;
}


void SubBand :: seek(long framePos, SampleCountType samplePos) {
  if(sub) sub->seek(framePos,samplePos);
  this->samplePos = samplePos;
  sms->seek(framePos*nGrainsPerFrame);
  nFramesAssigned = framePos;
  nFramesRendered = framePos;
  synthFramePos = framePos;
  nGrainsSynthed = framePos*nGrainsPerFrame;
  nFramesRead = framePos;
}

long SubBand :: renderFromCache() 
{
  for(list<SBSMSRenderer*>::iterator i = renderers.begin(); i != renderers.end(); ++i) {
    SBSMSRenderer *renderer = *i;
    renderer->startFrame();
  }
  renderFromCacheInit(true);
  render();
  stepRenderFrame();

  for(list<SBSMSRenderer*>::iterator i = renderers.begin(); i != renderers.end(); ++i) {
    SBSMSRenderer *renderer = *i;
    renderer->endFrame();
  }
  long samples = outputFrameSize.read();
  outputFrameSize.advance(1);
  stepReadFrame();
  return samples;
}

long SubBand :: renderFromCacheInit(bool bSet) 
{
  long n;
  if(sub) {
    n = res * sub->renderFromCacheInit(bSet);
  } else {
    n = max(0L,min(1L,min((long)(nFramesAssigned-nFramesRendered)-nRenderLatency,
                          nRenderSlack-(long)(nFramesRendered-nFramesRead))));
  }
  if(bSet) {
    nGrainsRendered = 0;
    nGrainsToRender = n;
  }
  return n;
}

long SubBand :: assignFromCache(float stretch, float pitch)
{
  if(assignFromCacheInit(false)) {
    setPitch(pitch);
    setStretch(stretch);
    assignFromCacheInit(true);
    assignFromCache();
    stepAssignFrame();
    return inputFrameSize;
  } else {
    return 0;
  }
}

long SubBand :: assignFromCache()
{
  long ntodo = parent?1:nGrainsToAssign;
  long ndone = 0;

  while(ndone<ntodo) {
    if(nGrainsAssigned%res==0)
      if(sub) sub->assignFromCache();
    sms->assignTrackPointsFromCache(cache);
    nGrainsAssigned++;
    ndone++;
  }
  return ndone;
}

long SubBand :: assignFromCacheInit(bool bSet)
{
  long n = 1;
  if(sub) n = res*sub->assignFromCacheInit(bSet);
  if(!sub) {
    n = min(n,max(0L,min(1L,(long)(cache->frames-nFramesAssigned))));
  }
  if(bSet) {
    nGrainsAssigned = 0;
    nGrainsToAssign = n;
  }
  return n;
}

void SubBand :: writingComplete()
{
  if(sub) sub->writingComplete();
  nRenderLatency = 0;
  bWritingComplete = true;
}


long SubBand :: synthFromCache(audio *buf, long n, float stretch, float pitch, Debugger *dbg) 
{
  if(!n) return 0;
  long nSynth;
  float stretch2 = stretch;
  nSynth = synthFromCacheInit(n,&stretch2,dbg);
  this->synth(nSynth,stretch2,pitch,dbg);
  nReadFromOutputFrame += nSynth;
  while(nReadFromOutputFrame >= inputFrameSize) {
    nReadFromOutputFrame -= inputFrameSize;
  }
  long nRead = 0;
  long nToRead = -1;
  readSubSamples();
  while(nToRead && nRead < n) {
    nToRead = n - nRead;
    nToRead = outMixer->read(buf+nRead, nToRead);
    nRead += nToRead;
  }
  return nRead;
}


long SubBand :: getDrop(float stretch)
{
  if(!bSynthStarted) {
    bool bBackwards = (stretch<0.0f);
    long pos = (samplePos%inputFrameSize);
    if(bBackwards) {
      nToDrop = inputFrameSize - pos;
    } else {
      nToDrop = pos;
    }
  }
  return nToDrop;
}

long SubBand :: synthFromCacheInit(long n, float *stretch, Debugger *dbg)
{
  if(!parent && 
     (sms->isPastLeft() && *stretch < 0.0f ||
      sms->isPastRight() && *stretch >= 0.0f)) return 0;

  if(!bSynthStarted) {
    bBackwards = (*stretch<0.0f);
    if(!parent) {
      int pos = samplePos%inputFrameSize;
      nToDrop = 0;
      if(bBackwards) {
        if(pos) {
          synthFramePos++;
          sms->seek(synthFramePos*nGrainsPerFrame);
          nToDrop = inputFrameSize - pos;
        }
      } else {
        if(pos) {
          nToDrop = pos;
        }
      }
    }
    if(nToDrop) {
      // no inputs
      //prepadInputs(nToDrop);
    }
    // blaa
    sms->assignTrackPointsFromCache(bBackwards,0,cache,dbg);
    bSynthStarted = true;
  }
  int ism = synthFramePos;
  if(bSynthGrainStart) {
    bool bBackwardsNext = (*stretch<0.0f);
    if(bBackwards != bBackwardsNext) {
      sms->pruneTracks();
    }
    bBackwards = bBackwardsNext;
    bSynthGrainStart = false;
    if(bBackwards && nGrainsSynthed%nGrainsPerFrame == 0) {
      ism = ism - 1;
    }
    // XXXX fafefawef
    if(!sms->assignTrackPointsFromCache(bBackwards,bBackwards?-1:1,cache,dbg)) {
      n = 0;
    }
  }
  if(!parent && nToDrop) {
    n = min(n,nToDrop);
    if(bBackwards) {
      *stretch = -1.0f;
    } else {
      *stretch = 1.0f;
    }
  }
  if(!parent) {
    long n0 = sms->synthInit(n,bBackwards,*stretch);
    n = min(n,n0);
  }
  if(n == 0) bSynthGrainStart = true;
  if(sub) sub->synthFromCacheInit(n,stretch,dbg);
  return n;
}

long SubBand :: synth(long n, float stretch, float pitch, Debugger *dbg) {
  if(sub) sub->synth(n,stretch,pitch,dbg);
  int grainCompleted = sms->synthTracks(n,synthRenderer,bBackwards,stretch,pitch,dbg);
  if(!parent && grainCompleted) {
    stepSynthGrain(bBackwards?-grainCompleted:grainCompleted);
  }
  return n;
}

void SubBand :: stepSynthGrain(int n)
{
  int n0;
  if(bBackwards) {
    if(n == -1) {
      nGrainsSynthed--;
      n0 = -1;
    } else {
      n0 = 0;
    }
  } else {
    if(n == 1) {
      nGrainsSynthed++;
      n0 = 1;
    } else {
      n0 = 0;
    } 
  }
  if(sub && nGrainsSynthed%res == 0) sub->stepSynthGrain(n);
  sms->stepSynth(n0);
  synthFramePos = nGrainsSynthed / nGrainsPerFrame;
  bSynthGrainStart = true;
}


void SubBand :: setLeftPos(SampleCountType pos)
{
  sms->setLeftPos(pos);
}

void SubBand :: setRightPos(SampleCountType pos)
{
  sms->setRightPos(pos);
}

SampleCountType SubBand :: getSamplePos()
{
  return sms->getSamplePos();
}


}
