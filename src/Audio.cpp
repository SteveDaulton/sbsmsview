#include "Audio.h"
#include <algorithm>
#include "sms.h"
using namespace std;

void *writeThreadCB(void *data) 
{
  SBSMSAudioSource *source = (SBSMSAudioSource*)data;
  while(source->isPlaying()) {
    source->write();
  }
  source->writingComplete();
  pthread_exit(NULL);
  return NULL;
}


SBSMSAudioSource :: SBSMSAudioSource(SBSMS *sbsmsSrc, _sbsms_::Debugger *debugger)
{
  dbg = debugger;
  sbsms = new SBSMS(sbsmsSrc,true);
  iface = new SBSMSInterfaceVariableRate(sbsms->getTotalSamples());
  sbsms->reset(false);
  sbsms->seek(iface,0);
  channels = sbsms->getChannels();
  pthread_mutex_init(&sbsmsMutex, NULL);
  pthread_mutex_init(&playMutex, NULL);
  rb = NULL;
  fbuf = NULL;
  outBuf = NULL;
  abuf = NULL;
  frameSize = sbsms->getInputFrameSize();
  Fs = 44100.0f;
  denv = 1.0f / (0.001f * (float)Fs);
  bPlaying = false;
  bDonePlaying = false;
  bPlayedToEnd = false;
  bWriteThread = false;
  setRate(1.0f);
  
}

void SBSMSAudioSource :: writingComplete()
{
  rb->flush();
  rb->writingComplete();
  pthread_mutex_lock(&sbsmsMutex);
  sbsms->reset(false);
  pthread_mutex_unlock(&sbsmsMutex);
}

long SBSMSAudioSource :: write()
{
  SampleCountType nextPos;
  pthread_mutex_lock(&sbsmsMutex);
  long nWrite;
  if(bCache) {
    int end = min(rightPos,sbsms->getTotalSamples());
    nWrite = min(frameSize,(int)(end-currPos));
    fprintf(stderr,"%ld %ld %ld\n",end,nWrite,currPos);

    memcpy(abuf,sbsms->getCache(0)->wave.getReadBuf()+currPos,nWrite*sizeof(audio));
  } else {
    nWrite = sbsms->synthFromCache(iface,abuf,frameSize,dbg);
  }
  if(!nWrite) rb->writingComplete();
  if(bCache) {
    currPos += nWrite;
  } else {
    nextPos = sbsms->getSamplePos();
  }
  pthread_mutex_unlock(&sbsmsMutex);
  
  if(nWrite) {
    if(channels==1) {
      for(int k=0;k<nWrite;k++) {	
        fbuf[k] = env*abuf[k][0];
        env += denv;
        if(env > 1.0f) env = 1.0f;
      }
    } else if(channels==2) {
      for(int k=0;k<nWrite;k++) {
        int k2 = k<<1;
        fbuf[k2] = env*abuf[k][0];
        fbuf[k2+1] = env*abuf[k][1];
        env += denv;
        if(env > 1.0f) env = 1.0f;
      }
    }
    rb->write(fbuf,nWrite);
  } else {
    bWriting = false;
    bPlayedToEnd = true;
  }

  float rate = getRate() * 44100.0f / Fs;
  if(!bCache) {
    if(rate != 0.0f) {
      currPos = nextPos - (int)((float)rb->getSamplesQueued() * rate);
      currPos = max(leftPos,min(rightPos,currPos));
    }
  }

  return nWrite;
}

long SBSMSAudioSource :: read(float *buf, long n, bool bFlush)
{
  n = rb->read(buf,n,bFlush);
  return n;
}

bool SBSMSAudioSource :: play(bool bCache)
{
  if(isPlaying()) return false;
  pthread_mutex_lock(&playMutex);
  setNextReadPosition(leftPos);
  env = 0.0f;
  bDonePlaying = false;
  bPlayedToEnd = false;
  bWriting = true;
  this->bCache = bCache;
  while(write() && !rb->isFull());
  bPlaying = true;
  pthread_create(&writeThread, NULL, writeThreadCB, (void*)this);
  bWriteThread = true;
  pthread_mutex_unlock(&playMutex);
  return true;
}

bool SBSMSAudioSource :: stop()
{
  if(!isPlaying()) return false;
  rb->readingComplete();
  pthread_mutex_lock(&playMutex);
  bDonePlaying = true;
  bPlaying = false;
  bWriting = false;
  if(bWriteThread) {
    pthread_join(writeThread,NULL);
  }
  pthread_mutex_unlock(&playMutex);
  return true;
}

bool SBSMSAudioSource :: isPlaying()
{  
  return bPlaying;
}

bool SBSMSAudioSource :: isWriting()
{
  return bWriting;
}

bool SBSMSAudioSource :: isDonePlaying()
{
  return bDonePlaying;
}

bool SBSMSAudioSource :: isPlayedToEnd()
{
  return bDonePlaying && bPlayedToEnd;
}

void SBSMSAudioSource :: prepareToPlay(int blockSize, double sampleRate)
{
  if(rb) delete rb; rb = new AudioBuffer(6*blockSize,channels);
  bPlaying = false;
  bWriting = false;
  bDonePlaying = false;
  bPlayedToEnd = false;
  currPos = 0;
  if(abuf) free(abuf); abuf = (audio*)calloc(blockSize,sizeof(audio));
  if(fbuf) free(fbuf); fbuf = (float*)calloc(blockSize*2,sizeof(float));
  if(outBuf) free(outBuf); outBuf = (float*)calloc(blockSize*2,sizeof(float));

}

void SBSMSAudioSource :: releaseResources()
{
  bWriting = false;
  bPlaying = false;
  bDonePlaying = false;
  if(bWriteThread) {
    pthread_join(writeThread,NULL);
  }
  if(rb) delete rb; rb = NULL;
  if(abuf) free(abuf); abuf = NULL;
  if(outBuf) free(outBuf); outBuf = NULL;
}

float SBSMSAudioSource :: getRate()  const 
{
  return rate;
}

void SBSMSAudioSource :: setRate(float rate)
{
  this->rate = rate;
  iface->setRate(rate);
}

SampleCountType SBSMSAudioSource :: getCurrPos()
{
  return currPos;
}

void SBSMSAudioSource :: setLeftPos(SampleCountType pos)
{
  leftPos = pos;
  pthread_mutex_lock(&sbsmsMutex);
  sbsms->setLeftPos(leftPos);
  pthread_mutex_unlock(&sbsmsMutex);
}

void SBSMSAudioSource :: setRightPos(SampleCountType pos)
{
  rightPos = pos;
  pthread_mutex_lock(&sbsmsMutex);
  sbsms->setRightPos(rightPos);
  pthread_mutex_unlock(&sbsmsMutex);
}

void SBSMSAudioSource :: setNextReadPosition (int64 newPosition) 
{
  pthread_mutex_lock(&sbsmsMutex);
  sbsms->reset(false);
  sbsms->seek(iface,newPosition);
  currPos = newPosition;
  pthread_mutex_unlock(&sbsmsMutex);
}

 
int64 SBSMSAudioSource :: getNextReadPosition ()  const
{
  int64 ret;
  pthread_mutex_lock((pthread_mutex_t*)&sbsmsMutex);
  ret = sbsms->getSamplePos();
  pthread_mutex_unlock((pthread_mutex_t*)&sbsmsMutex);
  return ret;
}

int64 SBSMSAudioSource :: getTotalLength() const
{
  int64 ret;
  pthread_mutex_lock((pthread_mutex_t*)&sbsmsMutex);
  ret = (int64)lrintf(sbsms->getTotalSamples()/rate);
  pthread_mutex_unlock((pthread_mutex_t*)&sbsmsMutex);
}

bool SBSMSAudioSource :: isLooping() const  
{
  return false;
}
 
void SBSMSAudioSource :: setLooping(bool shouldLoop) 
{
}


void SBSMSAudioSource :: getNextAudioBlock(const AudioSourceChannelInfo& bufferToFill)
{
  bufferToFill.clearActiveBufferRegion();
  unsigned long pos = 0;
  long ret = isPlaying()?1:0;
  long chunk;
  bool bStop = false;

  while(pos<bufferToFill.numSamples && ret) {
    int chunk = bufferToFill.numSamples-pos;
    ret = rb->read(outBuf+pos*channels,chunk);
    pos += ret;
    if(!ret) bStop = true;
  }

  
  if(channels == 1)  {
    for(int c=0; c<bufferToFill.buffer->getNumChannels(); c++) {
      float *buf = bufferToFill.buffer->getWritePointer(c);
      memcpy(buf+bufferToFill.startSample,outBuf,pos*sizeof(float));
    }
  } else {
    float *buf0 = bufferToFill.buffer->getWritePointer(0);
    float *buf1 = bufferToFill.buffer->getWritePointer(1);
    for(int k=0; k<pos; k++) {
      int k2 = k<<1;
      buf0[k+bufferToFill.startSample] = outBuf[k2];
      buf1[k+bufferToFill.startSample] = outBuf[k2+1];
    }
  }

  if(bStop) {
    stop();
  }
}
