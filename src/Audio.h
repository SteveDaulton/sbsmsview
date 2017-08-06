#ifndef SBSMSAUDIO_H
#define SBSMSAUDIO_H

#include "audiobuffer.h"
#include <pthread.h>
#include "../JuceLibraryCode/JuceHeader.h"
#include "sbsms.h"
using namespace _sbsms_;

class SBSMSAudioSource : 
public PositionableAudioSource
{
 public:
  SBSMSAudioSource(SBSMS *sbsmsSrc, _sbsms_::Debugger *debugger);

  bool play(bool bWave);
  bool stop();
  bool isPlaying();
  bool isWriting();
  bool isDonePlaying();
  bool isPlayedToEnd();
  float getRate() const;
  void setRate(float rate);
  float getPitch() const;
  void setPitch(float pitch);
  void setLeftPos(SampleCountType pos);
  void setRightPos(SampleCountType pos);
  SampleCountType getCurrPos();

  void prepareToPlay(int samplesPerBlockExpected, double sampleRate) override;
  void releaseResources() override;
  void getNextAudioBlock(const AudioSourceChannelInfo& bufferToFill) override;
  void setNextReadPosition (int64 newPosition) override;
  int64 getNextReadPosition() const override; 
  int64 getTotalLength () const override;
  bool isLooping() const override;
  void setLooping(bool shouldLoop) override;

  long read(float *buf, long n, bool bFlush);
  long write();
  void writingComplete();

 protected:

  bool bCache;
  int frameSize;
  int channels;
  AudioBuffer *rb;
  SBSMS *sbsms;
  SBSMSInterfaceVariableRate *iface;
  float env;
  float denv;
  bool bWriteThread;
  bool bOpen;
  bool bPlaying;
  bool bDonePlaying;
  bool bPlayedToEnd;
  bool bWriting;
  float Fs;
  float rate;
  pthread_t writeThread;
  audio *abuf;
  float *fbuf;
  float *outBuf;

  SampleCountType leftPos;
  SampleCountType currPos;
  SampleCountType rightPos;
  pthread_mutex_t sbsmsMutex;
  pthread_mutex_t playMutex;

  _sbsms_::Debugger *dbg;

};

#endif
