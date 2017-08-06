#ifndef TFVIEW_H
#define TFVIEW_H

#include "trackpoint.h"
#include "sbsms.h"
#include "Audio.h"

#include <set>
#include <map>
#include <vector>

using namespace std;
using namespace _sbsms_;
#include "../JuceLibraryCode/JuceHeader.h"

#include "OpenGL.h"


class TimeFrequencyView;


enum {
  TimeAxisMSM = 0,
  TimeAxisGrains,
  TimeAxisSamples,
  FrequencyAxisHz,
  FrequencyAxis2PI,
  FrequencyAxisNotes
  
};


class OpenGLChildComponent : public Component {
 public:
  virtual ~OpenGLChildComponent() {}
  virtual void renderOpenGL()=0;
  virtual void openGLContextClosing()=0;
  virtual void newOpenGLContextCreated(OpenGLContext *openGlContext)=0;

  OpenGLContext *openGlContext;
};

class WaveDisplay : public OpenGLChildComponent {
public:
  WaveDisplay(TimeFrequencyView *view);
  ~WaveDisplay() {}
  void renderOpenGL();
  void openGLContextClosing();
  void newOpenGLContextCreated(OpenGLContext *openGlContext);
  
  void paint(Graphics &g) override;

  TimeFrequencyView *tfView;

  OpenGLShaderProgram shaderLines;

  GLuint positionAttribId;
  GLuint pvMatrixUniformId;
  GLuint colorUniformId;

  GLuint minMaxVAOId;
  GLuint minMaxVBOId;
  GLuint rmsVAOId;
  GLuint rmsVBOId;

};

class TrackView 
: public OpenGLChildComponent
{
 public:
  TrackView(TimeFrequencyView *tfView);
  ~TrackView() {}
  void paint(Graphics &g) override;

  // this is not an openglrenderer, called from parent openglrenderer
  void renderOpenGL();
  void newOpenGLContextCreated(OpenGLContext *openGlContext); 
  void openGLContextClosing();

  TimeFrequencyView *tfView;
  OpenGLShaderProgram shaderLines;

  GLuint positionAttribId;
  GLuint pvMatrixUniformId;
  GLuint colorUniformId;

  GLuint trackVAOId;
  GLuint trackVBOId;
  GLuint trackScaleVAOId;
  GLuint trackScaleVBOId;
  GLuint trackDtVAOId;
  GLuint trackDtVBOId;
  GLuint trackDPhVAOId;
  GLuint trackDPhVBOId;
  GLuint trackMeritVAOId;
  GLuint trackMeritVBOId;
  GLuint totalMeritVAOId;
  GLuint totalMeritVBOId;
  GLuint trackOnsetVAOId;
  GLuint trackOnsetVBOId;
  GLuint trackOffsetVAOId;
  GLuint trackOffsetVBOId;
};

class SpectrumView 
: public OpenGLChildComponent
{
 public:
  SpectrumView(TimeFrequencyView *tfView, int which, float *color);
  ~SpectrumView() {}
  void paint(Graphics &g) override;

  // this is not an openglrenderer, called from parent openglrenderer
  void renderOpenGL();
  void newOpenGLContextCreated(OpenGLContext *openGlContext); 
  void openGLContextClosing();

  TimeFrequencyView *tfView;
  int which;
  OpenGLShaderProgram shaderLines;

  GLuint positionAttribId;
  GLuint pvMatrixUniformId;
  GLuint colorUniformId;

  GLuint graphVAOId[10];
  GLuint graphVBOId[10];

  GLuint cutsVAOId[10];
  GLuint cutsVBOId[10];

  float color[4];
};


class StatusListener {
public:
  virtual void statusChanged(SampleCountType pos, float f, int command)=0;
  virtual void statusChanged(const Rectangle<float> &r, int command)=0;
};

class StatusView : public Component
{
public:
  int band[2];
  int k[2];
  int time[2];
  float freq;
  int numBands;
  SampleCountType pos;

  StatusView();
  void setStatus(float freq, int *band, int *k, int *time, int numBands, SampleCountType pos);
  void paint(Graphics &g) override;
};

class TimeFrequencyRange
{
public:
  // view
  SampleCountType startSample;
  SampleCountType endSample;
  double startF;
  double endF;

  // play selection
  SampleCountType leftPos;
  SampleCountType currPos;
  SampleCountType rightPos;
  int frameSize;

  double xScale;
  double yScale;
};

class TimeFrequencyOverlay : public Component,
                             public TimeFrequencyRange
{
 public:
  TimeFrequencyOverlay();
  void resized() override;
  void rescale();
  void paint(Graphics& g) override;
  void addStatusListener(StatusListener *listener);
  void mouseExit (const MouseEvent& e) override;
  void mouseMove (const MouseEvent& e) override;
  void mouseUp (const MouseEvent& e) override;
  void mouseDown (const MouseEvent& e) override;
  void mouseDrag (const MouseEvent& e) override;
  vector<StatusListener*> listeners;
  int dragging;
  SampleCountType anchorPos;
  Point<float> selectAnchorPoint;
  Rectangle<float> selectionRect;
};

struct position {
  float x,y;
};

struct vertex 
{
  float x,y;
  unsigned char r,g,b,a;
};


struct pixel {
  unsigned char r,g,b,a;
};

struct texture {
  int width;
  int height;
  int offsetX;
  int offsetY;
  pixel *data;
};


struct texQuad {
  texture tex;
  GLuint texId;
  GLuint vaoId;
  GLuint vboId;
};

enum {
  colorMapSize = 1024
};

struct track {
  vector<TrackPoint*> trackpoints;
  vector<vertex> vertices;
  TimeType start;
};

class TimeFrequencyView 
: public Component,
  public ScrollBar::Listener,
  public KeyListener,
  public SBSMSRenderer,
  private OpenGLRenderer,
  public Slider::Listener,
  public Timer,
  public StatusListener,
  public _sbsms_::Debugger
{
public:
  TimeFrequencyView(SBSMS *sbsms);
  ~TimeFrequencyView();

  // opaque
  void paint(Graphics& g);  
  void render(const SBSMSRenderChunk &i, Track *t, Debugger *dbg);
	void resized();
  
  void sliderValueChanged(Slider * slider);
  void mouseWheelMove (const MouseEvent&, const MouseWheelDetails& d) override;
  void statusChanged(SampleCountType pos, float y, int command) override;
  void statusChanged(const Rectangle<float> &r, int command) override;

  bool shouldAssign(int band, TrackPoint *tp, TrackIndexType index);
  bool shouldScale(int band, TrackPoint *tp, TrackIndexType index, int c);
  bool shouldOnset(int band, TrackPoint *tp, TrackIndexType index);
  void selectAll();

  void markTrackOnset();
  void markTrackOffset();
  void outputSelectedTrack();
  void saveSelection();

  void setView(int pos0, int pos1, float f0, float f1);
  void zoomT(bool bIn);
  void zoomF(bool bIn);
  void setHScrollRange();
  void setVScrollRange();
  bool keyPressed(const KeyPress &key, Component * originatingComponent);
  void scrollBarMoved(ScrollBar *scrollBarThatHasMoved, double newRangeStart);
  ScrollBar *hScroll;
  ScrollBar *vScroll;
	SBSMS *sbsms;

  TrackPoint *selectedPoint;
  int selectedPointBand;
  TimeType selectedPointTime;
  TrackIndexType selectedPointTrackIndex;

  SpectrumView *spectrum2;
  SpectrumView *spectrum1;
  SpectrumView *spectrumTrial;
  TrackView *trackView;

  pixel colorMap[colorMapSize];
  void initColorMaps();
	
  bool bRenderTracks;
  bool bRenderSpectrum;

  int channel;
  bool bEnableBand[10];
  float botF[10];
  float topF[10];
  int bands;
  
  GLuint vboId[2][10];
  GLuint vaoId[2][10];

  OpenGLContext openGLContext;
  OpenGLShaderProgram shaderLines;
  OpenGLShaderProgram shaderTex;

  GLuint positionAttribId_Lines;
  GLuint colorAttribId_Lines;
  GLuint pvMatrixUniformId_Lines;
  GLuint colorScaleUniformId_Lines;

  GLuint positionAttribId_Tex;
  GLuint textureCoordAttribId_Tex;
  GLuint pvMatrixUniformId_Tex;
  GLuint textureMapUniformId_Tex;
  
  pthread_mutex_t glMutex;
  void renderOpenGL();
  map<TrackIndexType, track > trackVBO[2][10];
  void newOpenGLContextCreated(); 
  void openGLContextClosing();
  map<TrackIndexType,unsigned long> trackStartIndex[2][10];
  map<TrackIndexType,TimeType> trackStartTime[2][10];
  long selectedTrackOffset;
  long selectedTrackSize;
  long selectedTrackIndex;
  long selectedTrackTime;
  int selectedTrackBand;
  bool bSelectedTrackLocked;
  set<TimeType> trackOnsets;
  set<TimeType> trackOffsets;
  
  set<TrackIndexType> selected[10];
  vector<int> selectedOffsets[10];
  vector<int> selectedCounts[10];
  vector<vertex> allVBO[2][10];
	vector<int> allCounts[2][10];
	vector<int> allOffsets[2][10];
  vector<texQuad> spectrogram[3][2][10];

  int timeAxisMode;
  int freqAxisMode;
  int scrollSize;
  int xAxisHeight;
  int yAxisWidth;
  int controlsHeight;
  int spectrumWidth;
  int trackViewHeight;
  int waveHeight;
  Slider *rateCtrl;

  WaveDisplay *waveDisplay;
  TimeFrequencyOverlay *view;
  StatusView *status;
  
	SampleCountType totSamples;
  double minF;
  double maxF;
  float cursorF;
  int cursorTime;
  SampleCountType cursorSample;

  int whichgram;
  bool bSelectAll;
  bool bAdjust;
  bool bOnset;
  bool bForceAdjust;

  BufferingAudioSource *audioSource;
  AudioDeviceManager audioDeviceManager;
  AudioSourcePlayer audioSourcePlayer;
  SBSMSAudioSource *sbsmsAudioSource;
  TimeSliceThread *audioThread;

  void timerCallback() override;
};



#endif

