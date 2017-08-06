#include "config.h"
#include "sbsms.h"
#include "real.h"
#include "subband.h"
#ifdef MULTITHREADED
#include <pthread.h>
#endif
#include <stdlib.h>
#include <string.h>
#include <algorithm>
using namespace std;

namespace _sbsms_ {

const SBSMSQualityParams SBSMSQualityStandard = {
  8,3,
  {512,512,512,384,384,384,384,384,0,0},
  {168,144,128,96,64,36,24,14,0,0},
  //{384,288,256,168,128,84,52,28,0,0},
  {168,144,128,96,64,36,24,14,0,0},
  {512,448,360,288,192,108,84,44,0,0},
  //{512,456,400,320,256,180,108,48,0,0},
  {1,1,2,1,1,2,1,1,0,0}
};


SBSMSQuality :: SBSMSQuality(const SBSMSQualityParams *params)
{
  this->params = *params;
}

long SBSMSQuality :: getFrameSize()
{
  return (1<<(params.bands-1)) * params.H;
}

long SBSMSQuality :: getMaxPresamples()
{
  long prepad = 0;
  for(int i=0; i<params.bands; i++) {
    prepad = max(prepad,(long)((1<<i) * (params.N2[i]>>1)));
  }
  prepad += ((1<<(params.bands-1)) - 1) * (NDownSample>>1);
  long framesize = (1<<(params.bands-1)) * params.H;
  long frames = prepad / framesize;
  if(prepad%framesize) frames++;
  frames++;
  prepad = frames * framesize;
  return prepad;
}

int SBSMSQuality :: getBandGrainsPerFrame(int band)
{
  int g = 1;
  for(int i=band; i<params.bands; i++) {
    g *= params.res[i];
  }
  return g;
}

#ifdef MULTITHREADED
class ThreadInterface;
#endif

class SBSMSImp {
public:
  SBSMSImp(SBSMSImp *src, bool bSynthesize);
  SBSMSImp(int channels, SBSMSQuality *quality, bool bSynthesize, bool bCache);
  ~SBSMSImp();
  inline long read(SBSMSInterface *iface, audio *buf, long n);
  inline void addRenderer(SBSMSRenderer *renderer);
  inline void removeRenderer(SBSMSRenderer *renderer);
  inline long renderFrame(SBSMSInterface *iface);
  inline long renderFrameFromCache(SBSMSInterface *iface);
  inline long synthFromCache(SBSMSInterface *iface, audio *buf, long n, Debugger *dbg);
  inline void reset(bool bFlushInput);
  inline void setTotalSamples(const SampleCountType &samples);
  inline SampleCountType getSamplePos();
  inline SampleCountType getTotalSamples();
  inline void setLeftPos(SampleCountType pos);
  inline void setRightPos(SampleCountType pos);
  inline void seek(SBSMSInterface *iface,SampleCountType samplePos);  
  inline Cache *getCache(int band);
  inline int getChannels();
  inline SBSMSQuality *getQuality();

  SubBand *top;  
#ifdef MULTITHREADED
  friend class ThreadInterface;
  ThreadInterface *threadInterface;
#endif
  friend class SBSMS;
protected:
  float getInputTime(SBSMSInterface *iface);
  long write(SBSMSInterface *);
  long assignFromCache(SBSMSInterface *iface);
  SBSMSError error;
  SBSMSImp *source;
  long nPrepad;
  long nPrepadDone;
  long nPresamplesDone;
  SampleCountType totalSamples;
  SampleCountType nSamplesInputed;
  SampleCountType nSamplesOutputed;
  int channels;
  SBSMSQuality *quality;
  audio *ina;
  GenFactory *genFactory;
  TrackStitcher *stitcher;
};

#ifdef MULTITHREADED

struct channel_thread_data {
  int c;
  ThreadInterface *threadInterface;
};

struct analyze_thread_data {
  int i;
  ThreadInterface *threadInterface;
};

class ThreadInterface {
public:
  friend class SBSMSImp;
  ThreadInterface(SBSMSImp *sbsms, bool bSynthesize);  
  ~ThreadInterface();
  void signalReadWrite();
  void signalAnalyze();
  void signalExtract(int c);
  void signalMark(int c);
  void signalAssign(int c);
  void signalTrial(int c);
  void signalAdjust();
  void signalRender(int c);
  void waitReadWrite();
  void waitAnalyze(int i);
  void waitExtract(int c);
  void waitAssign(int c);
  void waitTrial(int c);
  void waitAdjust();
  void waitRender(int c);
  SubBand *top;
  int channels;
  pthread_mutex_t readWriteMutex;
  pthread_cond_t readWriteCond;
  pthread_t analyzeThread[NumGrainTypes];
  pthread_mutex_t analyzeMutex[NumGrainTypes];
  pthread_cond_t analyzeCond[NumGrainTypes];
  pthread_t extractThread[2];
  pthread_mutex_t extractMutex[2];
  pthread_cond_t extractCond[2];
  pthread_t assignThread[2];
  pthread_mutex_t assignMutex[2];
  pthread_cond_t assignCond[2];
  pthread_t trialThread[2];
  pthread_mutex_t trialMutex[2];
  pthread_cond_t trialCond[2];
  pthread_t adjustThread;
  pthread_mutex_t adjustMutex;
  pthread_cond_t adjustCond;
  bool bRenderThread;
  pthread_t renderThread[2];
  pthread_mutex_t renderMutex[2];
  pthread_cond_t renderCond[2];
  channel_thread_data channelData[2];
  analyze_thread_data analyzeData[NumGrainTypes];
  bool bActive;
};

void *analyzeThreadCB(void *data) {
  analyze_thread_data *analyzeData = (analyze_thread_data*)data;
  ThreadInterface *threadInterface = analyzeData->threadInterface;
  SubBand *top = threadInterface->top;
  int i = analyzeData->i;
  int channels = threadInterface->channels;
  while(threadInterface->bActive) {
    threadInterface->waitAnalyze(i);
    if(top->analyzeInit(i,true)) {
      top->analyze(i);
      top->stepAnalyzeFrame(i);
      threadInterface->signalReadWrite();
      for(int c=0; c<channels; c++) {
        threadInterface->signalExtract(c);
      }
    }
  }
  pthread_exit(NULL);
  return NULL;
}

void *extractThreadCB(void *data) {
  channel_thread_data *channelData = (channel_thread_data*)data;
  ThreadInterface *threadInterface = channelData->threadInterface;
  SubBand *top = threadInterface->top;
  int c = channelData->c;
  while(threadInterface->bActive) {
    threadInterface->waitExtract(c);
    if(top->extractInit(c,true)) {
      top->extract(c);
      top->stepExtractFrame(c);
      threadInterface->signalAnalyze();
      threadInterface->signalMark(c);        
    }
  }
  pthread_exit(NULL);
  return NULL;
}

void *assignThreadCB(void *data) {
  channel_thread_data *channelData = (channel_thread_data*)data;
  ThreadInterface *threadInterface = channelData->threadInterface;
  SubBand *top = threadInterface->top;
  int c = channelData->c;
  while(threadInterface->bActive) {
    threadInterface->waitAssign(c);
    if(top->markInit(c,true)) {
      top->mark(c);
      top->stepMarkFrame(c);      
      threadInterface->signalExtract(c);
    }
    if(top->assignInit(c,true)) {
      top->assign(c);
      top->advance(c);
      top->stepAssignFrame(c);
      threadInterface->signalTrial(c);
    }
  }
  pthread_exit(NULL);
  return NULL;
}

void *trialThreadCB(void *data) {
  channel_thread_data *channelData = (channel_thread_data*)data;
  ThreadInterface *threadInterface = channelData->threadInterface;
  SubBand *top = threadInterface->top;
  int c = channelData->c;
  while(threadInterface->bActive) {
    threadInterface->waitTrial(c);
    if(top->trialInit(c,true)) {
      top->trial(c);
      top->stepTrialFrame(c);
      threadInterface->signalAssign(c);
      threadInterface->signalAdjust();
    }
  }
  pthread_exit(NULL);
  return NULL;
}

void *adjustThreadCB(void *data) {
  ThreadInterface *threadInterface = (ThreadInterface*)data;
  int channels = threadInterface->channels;
  SubBand *top = threadInterface->top;
  while(threadInterface->bActive) {
    threadInterface->waitAdjust();
    if(top->adjustInit(true)) {
      top->adjust();
      top->stepAdjustFrame();
      for(int c=0; c<channels; c++) {
        threadInterface->signalTrial(c);
      }
      if(threadInterface->bRenderThread) {
        for(int c=0; c<channels; c++) {
          threadInterface->signalRender(c);
        }
      } else {
        threadInterface->signalReadWrite();
      }
    }
  }
  pthread_exit(NULL);
  return NULL;
}

void *renderThreadCB(void *data) {
  channel_thread_data *channelData = (channel_thread_data*)data;
  ThreadInterface *threadInterface = channelData->threadInterface;
  SubBand *top = threadInterface->top;
  int c = channelData->c;
  while(threadInterface->bActive) {
    threadInterface->waitRender(c);
    if(top->renderInit(c,true)) {
      top->render(c);
      top->stepRenderFrame(c);
      threadInterface->signalAdjust();
      threadInterface->signalReadWrite();
    }
  }
  pthread_exit(NULL);
  return NULL;
}

ThreadInterface :: ThreadInterface(SBSMSImp *sbsms, bool bSynthesize) 
{
  this->top = sbsms->top;
  this->channels = sbsms->channels;
  bActive = true;
  pthread_cond_init(&readWriteCond, NULL);
  pthread_mutex_init(&readWriteMutex, NULL);
  if(bSynthesize) {
    bRenderThread = true;
  } else {
    bRenderThread = false;
  }
  for(int i=0; i<NumGrainTypes; i++) {
    analyzeData[i].i = i;
    analyzeData[i].threadInterface = this;
    pthread_cond_init(&analyzeCond[i], NULL);
    pthread_mutex_init(&analyzeMutex[i], NULL);
  }
  for(int c=0; c<channels; c++) {
    channelData[c].c = c;
    channelData[c].threadInterface = this;
    pthread_cond_init(&extractCond[c], NULL);
    pthread_mutex_init(&extractMutex[c], NULL);
    pthread_cond_init(&assignCond[c], NULL);    
    pthread_mutex_init(&assignMutex[c], NULL);
    pthread_cond_init(&trialCond[c], NULL);    
    pthread_mutex_init(&trialMutex[c], NULL);
    if(bRenderThread) {
      pthread_cond_init(&renderCond[c], NULL);
      pthread_mutex_init(&renderMutex[c], NULL);
    }
  }
  for(int i=0; i<NumGrainTypes; i++) {
    pthread_create(&analyzeThread[i], NULL, analyzeThreadCB, (void*)&analyzeData[i]);
  }
  for(int c=0; c<channels; c++) {
    pthread_create(&extractThread[c], NULL, extractThreadCB, (void*)&channelData[c]);
    pthread_create(&assignThread[c], NULL, assignThreadCB, (void*)&channelData[c]);
    pthread_create(&trialThread[c], NULL, trialThreadCB, (void*)&channelData[c]);
    if(bRenderThread) {
      pthread_create(&renderThread[c], NULL, renderThreadCB, (void*)&channelData[c]);
    }
  }
  pthread_cond_init(&adjustCond, NULL);    
  pthread_mutex_init(&adjustMutex, NULL);
  pthread_create(&adjustThread, NULL, adjustThreadCB, this);
}

ThreadInterface :: ~ThreadInterface() 
{
  bActive = false;
  for(int i=0; i<NumGrainTypes; i++) {
    pthread_mutex_lock(&analyzeMutex[i]);
    pthread_cond_broadcast(&analyzeCond[i]);
    pthread_mutex_unlock(&analyzeMutex[i]);
    pthread_join(analyzeThread[i],NULL);
  }
  for(int c=0; c<channels; c++) {
    pthread_mutex_lock(&extractMutex[c]);
    pthread_cond_broadcast(&extractCond[c]);
    pthread_mutex_unlock(&extractMutex[c]);
    pthread_join(extractThread[c],NULL);
    pthread_mutex_lock(&assignMutex[c]);
    pthread_cond_broadcast(&assignCond[c]);
    pthread_mutex_unlock(&assignMutex[c]);
    pthread_join(assignThread[c],NULL);

    pthread_mutex_lock(&trialMutex[c]);
    pthread_cond_broadcast(&trialCond[c]);
    pthread_mutex_unlock(&trialMutex[c]);
    pthread_join(trialThread[c],NULL);
    pthread_mutex_lock(&adjustMutex);
    pthread_cond_broadcast(&adjustCond);
    pthread_mutex_unlock(&adjustMutex);
    pthread_join(adjustThread,NULL);

    if(bRenderThread) {
      pthread_mutex_lock(&renderMutex[c]);
      pthread_cond_broadcast(&renderCond[c]);
      pthread_mutex_unlock(&renderMutex[c]);
      pthread_join(renderThread[c],NULL);
    }
  }
}

void ThreadInterface :: signalReadWrite() 
{
  pthread_mutex_lock(&readWriteMutex);
  bool bReady;
  if(bRenderThread) {
    bReady = (top->writeInit() || top->readInit());
  } else {
    if(top->writeInit()) {
      bReady = true;  
    } else {
      bReady = true;
      for(int c=0; c<channels; c++) {
        if(!top->renderInit(c,false)) {
          bReady = false;
          break;
        }
      }
    }
  }
  if(bReady) {
    pthread_cond_broadcast(&readWriteCond);
  }
  pthread_mutex_unlock(&readWriteMutex);
}

void ThreadInterface :: signalAnalyze() 
{
  for(int i=0; i<NumGrainTypes; i++) {
    pthread_mutex_lock(&analyzeMutex[i]);
    if(top->analyzeInit(i,false)) {
      pthread_cond_broadcast(&analyzeCond[i]);
    }
    pthread_mutex_unlock(&analyzeMutex[i]);
  }
}

void ThreadInterface :: signalExtract(int c) {
  pthread_mutex_lock(&extractMutex[c]);
  if(top->extractInit(c,false)) {
    pthread_cond_broadcast(&extractCond[c]);
  }
  pthread_mutex_unlock(&extractMutex[c]);
}

void ThreadInterface :: signalMark(int c) {
  pthread_mutex_lock(&assignMutex[c]);
  if(top->markInit(c,false)) {
    pthread_cond_broadcast(&assignCond[c]);
  }
  pthread_mutex_unlock(&assignMutex[c]);
}

void ThreadInterface :: signalAssign(int c) {
  pthread_mutex_lock(&assignMutex[c]);
  if(top->assignInit(c,false)) {
    pthread_cond_broadcast(&assignCond[c]);
  }
  pthread_mutex_unlock(&assignMutex[c]);
}

void ThreadInterface :: signalTrial(int c) {
  pthread_mutex_lock(&trialMutex[c]);
  if(top->trialInit(c,false)) {
    pthread_cond_broadcast(&trialCond[c]);
  }
  pthread_mutex_unlock(&trialMutex[c]);
}

void ThreadInterface :: signalAdjust() {
  pthread_mutex_lock(&adjustMutex);
  if(top->adjustInit(false)) {
    pthread_cond_broadcast(&adjustCond);
  }
  pthread_mutex_unlock(&adjustMutex);
}

void ThreadInterface :: signalRender(int c) {
  pthread_mutex_lock(&renderMutex[c]);
  if(top->renderInit(c,false)) {
    pthread_cond_broadcast(&renderCond[c]);
  }
  pthread_mutex_unlock(&renderMutex[c]);
}

void ThreadInterface :: waitReadWrite() {
  pthread_mutex_lock(&readWriteMutex);
  bool bReady;
  if(bRenderThread) {
    bReady = (top->writeInit() || top->readInit());
  } else {
    if(top->writeInit()) {
      bReady = true;  
    } else {
      bReady = true;
      for(int c=0; c<channels; c++) {
        if(!top->renderInit(c,false)) {
          bReady = false;
          break;
        }
      }
    }
  }
  if(!bReady) {
    pthread_cond_wait(&readWriteCond,&readWriteMutex);
  }
  pthread_mutex_unlock(&readWriteMutex);
}

void ThreadInterface :: waitAnalyze(int i) {
  pthread_mutex_lock(&analyzeMutex[i]);
  if(!top->analyzeInit(i,false)) {
    pthread_cond_wait(&analyzeCond[i],&analyzeMutex[i]);
  }
  pthread_mutex_unlock(&analyzeMutex[i]);
}

void ThreadInterface :: waitExtract(int c) {
  pthread_mutex_lock(&extractMutex[c]);
  if(!top->extractInit(c,false)) {
    pthread_cond_wait(&extractCond[c],&extractMutex[c]);
  }
  pthread_mutex_unlock(&extractMutex[c]);
}

void ThreadInterface :: waitAssign(int c) {
  pthread_mutex_lock(&assignMutex[c]);
  if(!top->markInit(c,false) && !top->assignInit(c,false)) {
    pthread_cond_wait(&assignCond[c],&assignMutex[c]);
  }
  pthread_mutex_unlock(&assignMutex[c]);
}

void ThreadInterface :: waitTrial(int c) {
  pthread_mutex_lock(&trialMutex[c]);
  if(!top->trialInit(c,false)) {
    pthread_cond_wait(&trialCond[c],&trialMutex[c]);
  }
  pthread_mutex_unlock(&trialMutex[c]);
}

void ThreadInterface :: waitAdjust() {
  pthread_mutex_lock(&adjustMutex);
  if(!top->adjustInit(false)) {
    pthread_cond_wait(&adjustCond,&adjustMutex);
  }
  pthread_mutex_unlock(&adjustMutex);
}

void ThreadInterface :: waitRender(int c) {
  pthread_mutex_lock(&renderMutex[c]);
  if(!top->renderInit(c,false)) {
    pthread_cond_wait(&renderCond[c],&renderMutex[c]);
  }
  pthread_mutex_unlock(&renderMutex[c]);
}

#endif

void SBSMS :: reset(bool bFlushInput) { imp->reset(bFlushInput); }
void SBSMSImp :: reset(bool bFlushInput)
{
  nSamplesInputed = 0;
  nSamplesOutputed = 0;
  nPrepadDone = 0;
  nPresamplesDone = 0;
  top->reset(bFlushInput);
  if(stitcher) stitcher->reset();
}

SBSMS :: SBSMS(int channels, SBSMSQuality *quality, bool bSynthesize, bool bCache)
{ imp = new SBSMSImp(channels,quality,bSynthesize,bCache); }
SBSMSImp :: SBSMSImp(int channels, SBSMSQuality *quality, bool bSynthesize, bool bCache)
{
  this->source = NULL;
  this->channels = channels;
  this->quality = new SBSMSQuality(&quality->params);
  error = SBSMSErrorNone;
  genFactory = new DefaultGenFactory(channels);
  stitcher = new TrackStitcher(genFactory,channels);
  SubBand *source = NULL;
  top = new SubBand(NULL,0,channels,quality,true,bSynthesize,source,stitcher,bCache);
  ina = (audio*)malloc(quality->getFrameSize()*sizeof(audio));
  nPrepad = quality->getMaxPresamples();
  reset(false);
#ifdef MULTITHREADED
  threadInterface = new ThreadInterface(this,bSynthesize);
#endif
}


SBSMS :: SBSMS(SBSMS *source, bool bSynthesize) { imp = new SBSMSImp(source->imp, bSynthesize); }
SBSMSImp :: SBSMSImp(SBSMSImp *source, bool bSynthesize)
{
  this->source = source;
  this->channels = source->channels;
  this->quality = new SBSMSQuality(&source->quality->params);
  setTotalSamples(source->getTotalSamples());
  error = SBSMSErrorNone;
  genFactory = new DefaultGenFactory(channels);
  stitcher = new TrackStitcher(genFactory,channels);
  top = new SubBand(NULL,0,channels,quality,false,bSynthesize,source->top,stitcher,false);  
  setLeftPos(0);
  setRightPos(getTotalSamples());
  ina = NULL;
  reset(false);
#ifdef MULTITHREADED
  threadInterface = new ThreadInterface(this,false,false);
#endif
}

SBSMS :: ~SBSMS() { delete imp; }
SBSMSImp :: ~SBSMSImp()
{
#ifdef MULTITHREADED
  if(threadInterface) delete threadInterface;
#endif
  delete genFactory;
  delete stitcher;
  if(top) delete top;
  if(ina) free(ina);
  if(quality) delete quality;
}

void SBSMS :: addRenderer(SBSMSRenderer *renderer) { imp->addRenderer(renderer); }
void SBSMSImp :: addRenderer(SBSMSRenderer *renderer)
{
  top->addRenderer(renderer);
}

void SBSMS :: removeRenderer(SBSMSRenderer *renderer) { imp->removeRenderer(renderer); }
void SBSMSImp :: removeRenderer(SBSMSRenderer *renderer)
{
  top->removeRenderer(renderer);
}

SBSMSError SBSMS :: getError()
{
  return imp->error;
}

float SBSMSImp :: getInputTime(SBSMSInterface *iface)
{
  return (float)nSamplesInputed / (float)iface->getSamplesToInput();
}

long SBSMSImp :: assignFromCache(SBSMSInterface *iface)
{
  long nWrite = 0;
  float t = getInputTime(iface);
  float stretch = iface->getStretch(t);
  float pitch = iface->getPitch(t);

  if(source) {
    nWrite = top->assignFromCache(stretch,pitch);
    nSamplesInputed += nWrite;
    if(nWrite == 0) {
      top->writingComplete();
    }
  }
  return nWrite;
}

long SBSMS :: renderFrameFromCache(SBSMSInterface *iface) { return imp->renderFrameFromCache(iface); }
long SBSMSImp :: renderFrameFromCache(SBSMSInterface *iface)
{
  long nRendered = 0;
  while(!nRendered && !top->isDone()) {
    assignFromCache(iface);
    if(top->renderFromCacheInit(false)) {
      nRendered = top->renderFromCache();
    } else {
      nRendered = 0;
    }
    if(nSamplesOutputed >= iface->getSamplesToOutput()) {
      top->renderComplete(iface->getSamplesToOutput());
    }
    nSamplesOutputed += nRendered;
  }
  return nRendered;
}

long SBSMS :: synthFromCache(SBSMSInterface *iface, audio *buf, long n, Debugger *dbg) { return imp->synthFromCache(iface,buf,n,dbg); }
long SBSMSImp :: synthFromCache(SBSMSInterface *iface, audio *buf, long n, Debugger *dbg)
{
  long nReadTotal = 0;
  float t = getInputTime(iface);
  float stretch = iface->getStretch(t);
  float pitch = iface->getPitch(t);

  long nRead = -1;
  while(nRead && nReadTotal < n) {
    long nToDrop = top->getDrop(stretch);
    if(nToDrop) {
      audio nullBuf[512];
      nRead = min(512L,top->nToDrop);
      nRead = top->synthFromCache(nullBuf,nRead,stretch,pitch,dbg);
      top->nToDrop -= nRead;
    } else {
      nRead = n - nReadTotal;
      nRead = top->synthFromCache(buf+nReadTotal,nRead,stretch,pitch,dbg);
      nReadTotal += nRead;
    }
  }

  nSamplesInputed += nReadTotal;

  return stretch<0.0f?-nReadTotal:nReadTotal;
}

long SBSMSImp :: write(SBSMSInterface *iface)
{
  long nWrite = 0;
  float t = getInputTime(iface);
  float stretch = iface->getStretch(t);
  float pitch = iface->getPitch(t);
  long nPresamples = iface->getPresamples();
  if(nPrepadDone < nPrepad - nPresamples) {
    stretch = 1.0f;
    nWrite = min(quality->getFrameSize(),nPrepad - nPresamples - nPrepadDone);
    memset(ina,0,nWrite*sizeof(audio));
    nPrepadDone += nWrite;
  } else if(nPresamplesDone < nPresamples) {
    stretch = 1.0f;
    nWrite = min(quality->getFrameSize(),nPresamples - nPresamplesDone);
    nWrite = iface->samples(ina,nWrite);
    if(nWrite == 0) {
      nWrite = quality->getFrameSize();
      memset(ina,0,nWrite*sizeof(audio));
    } else {
      nPresamplesDone += nWrite;
    }
  } else {
    nWrite = iface->samples(ina,quality->getFrameSize());
    nSamplesInputed += nWrite;
    if(nWrite == 0) {
      nWrite = quality->getFrameSize();
      memset(ina,0,nWrite*sizeof(audio));
    }
  }
  nWrite = top->write(ina, nWrite, stretch, pitch);

  return nWrite;
}

long SBSMS :: read(SBSMSInterface *iface, audio *buf, long n) { return imp->read(iface,buf,n); }
long SBSMSImp :: read(SBSMSInterface *iface, audio *buf, long n)
{
  long nReadTotal = 0;
  while(nReadTotal < n) {
    long nRead;
    nRead = n - nReadTotal;
    nRead = top->read(buf+nReadTotal,nRead);
    nReadTotal += nRead;
    if(nRead) {
#ifdef MULTITHREADED
      if(threadInterface->bRenderThread) {
        for(int c=0; c<channels; c++) {
          threadInterface->signalRender(c);
        }
      }
#endif
    } else {
#ifdef MULTITHREADED
      threadInterface->waitReadWrite();
#endif
      if(top->writeInit()) {
        write(iface);
#ifdef MULTITHREADED
        threadInterface->signalAnalyze();
#endif
      }
    }
#ifdef MULTITHREADED
    if(!threadInterface->bRenderThread) {
      for(int c=0; c<channels; c++) {     
        threadInterface->signalRender(c);
      }
    }
#else
    top->process(true);
#endif
    nSamplesOutputed += nRead;
  }
  return nReadTotal;
}

long SBSMS :: renderFrame(SBSMSInterface *iface) { return imp->renderFrame(iface); }
long SBSMSImp :: renderFrame(SBSMSInterface *iface)
{
  long nRendered = 0;
  while(!nRendered && !top->isDone()) {
    bool bReady = true;
    if(top->renderInit(false)) {
      nRendered = top->renderSynchronous();
    } else {
      nRendered = 0;
    }
    if(nRendered) {
#ifdef MULTITHREADED
      threadInterface->signalAdjust();
#endif
    } else {
#ifdef MULTITHREADED
      threadInterface->waitReadWrite();  
#endif      
      if(top->writeInit()) {
        write(iface);
      }
      
#ifdef MULTITHREADED
      threadInterface->signalAnalyze();
#endif
    }
#ifdef MULTITHREADED
#else
    top->process(false);
#endif
    if(nSamplesOutputed >= iface->getSamplesToOutput()) {
      top->renderComplete(iface->getSamplesToOutput());
    }
    nSamplesOutputed += nRendered;
  }
  return nRendered;
}

long SBSMS :: getInputFrameSize()
{
  return imp->top->getInputFrameSize();
}

SampleCountType SBSMS :: getSamplePos() { return imp->getSamplePos(); }
SampleCountType SBSMSImp :: getSamplePos()
{
  return top->getSamplePos();
}

SampleCountType SBSMS :: getTotalSamples() { return imp->getTotalSamples(); }
SampleCountType SBSMSImp :: getTotalSamples()
{
  return totalSamples;
}

void SBSMS :: setLeftPos(SampleCountType pos) { imp->setLeftPos(pos); }
void SBSMSImp :: setLeftPos(SampleCountType pos) { top->setLeftPos(pos); }
void SBSMS :: setRightPos(SampleCountType pos) { imp->setRightPos(pos); }
void SBSMSImp :: setRightPos(SampleCountType pos) { top->setRightPos(pos); }

void SBSMS :: setTotalSamples(const SampleCountType &samples) { imp->setTotalSamples(samples); }
void SBSMSImp :: setTotalSamples(const SampleCountType &samples)
{
  totalSamples = samples;
}


void SBSMS :: seek(SBSMSInterface *iface, SampleCountType samplePos) { imp->seek(iface,samplePos); }
void SBSMSImp :: seek(SBSMSInterface *iface, SampleCountType samplePos)
{
  long i = (long) (samplePos / top->getInputFrameSize());
  //if(fpIn) FSEEK(fpIn,frameByteOffset->at(i),SEEK_SET);  
  //nSamplesInputed = frameSampleOffset->at(i);
  nSamplesInputed = i * top->getInputFrameSize();
  float t = getInputTime(iface);
  float stretch = iface->getStretch(t);
  nSamplesOutputed = (SampleCountType)((double)nSamplesInputed * (double)stretch);
  top->seek(i,samplePos);
}

int SBSMS :: getChannels() { return imp->getChannels(); }
int SBSMSImp :: getChannels() { return channels; }

SBSMSQuality *SBSMS :: getQuality() { return imp->getQuality(); }
SBSMSQuality *SBSMSImp :: getQuality() { return quality; }


Cache *SBSMS :: getCache(int band) { return imp->getCache(band); }

Cache *SBSMSImp :: getCache(int band)
{
  return top->getCache(band);
}

class SBSMSInterfaceVariableRateImp {
public:
  SBSMSInterfaceVariableRateImp(SampleCountType samplesToInput);
  inline SampleCountType getSamplesToInput();
  inline SampleCountType getSamplesToOutput();
  friend class SBSMSInterfaceVariableRate;
protected:
  double stretch;
  float pitch;
  SampleCountType samplesToInput;
};

SBSMSInterfaceVariableRate :: SBSMSInterfaceVariableRate(SampleCountType samplesToInput)
{
  imp = new SBSMSInterfaceVariableRateImp(samplesToInput);
}

SBSMSInterfaceVariableRate :: ~SBSMSInterfaceVariableRate()
{
  delete imp;
}

SBSMSInterfaceVariableRateImp :: SBSMSInterfaceVariableRateImp(SampleCountType samplesToInput)
{
  this->samplesToInput = samplesToInput;
  this->stretch = 1.0;
  this->pitch = 1.0;
}

float SBSMSInterfaceVariableRate :: getStretch(float t)
{
  return (float)imp->stretch;
}

float SBSMSInterfaceVariableRate :: getPitch(float t)
{
  return imp->pitch;
}

void SBSMSInterfaceVariableRate :: setRate(float rate)
{
  imp->stretch = (rate == 0.0f?0.0f:1.0 / (double)rate);
}

void SBSMSInterfaceVariableRate :: setPitch(float pitch)
{
  imp->pitch = pitch;
}

long SBSMSInterfaceVariableRate :: getPresamples() 
{
  return 0L;
}

SampleCountType SBSMSInterfaceVariableRate :: getSamplesToInput() { return imp->getSamplesToInput(); }
SampleCountType SBSMSInterfaceVariableRateImp :: getSamplesToInput()
{
  return samplesToInput;
}

SampleCountType SBSMSInterfaceVariableRate :: getSamplesToOutput() { return imp->getSamplesToOutput(); }
SampleCountType SBSMSInterfaceVariableRateImp :: getSamplesToOutput()
{
  return (SampleCountType)((double)samplesToInput * stretch);
}

}
