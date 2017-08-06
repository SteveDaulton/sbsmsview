#include "TimeFrequencyView.h"
#include "convert.h"
#include "track.h"
#include <algorithm>
#include <gl-matrix.h>
#include "sms.h"
#include "pcm.h"
#include <pthread.h>

class ASCII {
public:
static const int a = 65;
static const int s = 83;
};

const char vertexShaderLinesUniformColor[] = {
  "attribute vec4 position;\n"
  "\n"
  "uniform mat4 pvMatrix;\n"
  "\n"
  "void main()\n"
  "{\n"
  "    gl_Position = pvMatrix * position;\n"
  "}\n"
};

const char fragmentShaderLinesUniformColor[] = {
  "uniform vec4 color;\n"
  "void main()\n"
  "{\n"
  "    gl_FragColor = color;"
  "}\n"
};

const char vertexShaderLines[] = {
  "attribute vec4 position;\n"
  "attribute vec4 sourceColour;\n"
  "\n"
  "uniform mat4 pvMatrix;\n"
  "\n"
  "varying vec4 destinationColour;\n"
  "\n"
  "void main()\n"
  "{\n"
  "    destinationColour = sourceColour;\n"
  "    gl_Position = pvMatrix * position;\n"
  "}\n"
};

const char fragmentShaderLines[] = {
  "uniform vec4 colorScale;\n"
#if JUCE_OPENGL_ES
  "varying lowp vec4 destinationColour;\n"
#else
  "varying vec4 destinationColour;\n"
#endif
  "\n"
  "void main()\n"
  "{\n"
  "    gl_FragColor = colorScale * destinationColour;"
  "}\n"
};



const char vertexShaderTex[] = {
  "attribute vec4 position;\n"
  "attribute vec2 textureCoordIn;\n"
  "\n"
  "uniform mat4 pvMatrix;\n"
  "\n"
  "varying vec2 textureCoordOut;\n"
  "\n"
  "void main()\n"
  "{\n"
  "    textureCoordOut = textureCoordIn;\n"
  "    gl_Position = pvMatrix * position;\n"
  "}\n"
};

const char fragmentShaderTex[] = {
#if JUCE_OPENGL_ES
  "varying lowp vec2 textureCoordOut;\n"
#else
  "varying vec2 textureCoordOut;\n"
#endif
  "\n"
  "uniform sampler2D textureMap;\n"
  "\n"
  "void main()\n"
  "{\n"
  "    gl_FragColor = texture2D (textureMap, textureCoordOut);\n"
  "}\n"
};

WaveDisplay :: WaveDisplay(TimeFrequencyView *tfView) : tfView(tfView), shaderLines(tfView->openGLContext)
{
}

void WaveDisplay :: paint(Graphics &g)
{
  pthread_mutex_lock(&tfView->glMutex);
  if(tfView->cursorSample >= 0) {
    int w = getWidth();
    int h = getHeight();
    TimeFrequencyOverlay *view = tfView->view;
    float x = view->xScale * (tfView->cursorSample - view->startSample);
    g.setColour(Colours::white.withAlpha(0.5f));
    g.drawLine(x,0,x,h);
  }
  pthread_mutex_unlock(&tfView->glMutex);
}

void WaveDisplay :: newOpenGLContextCreated(OpenGLContext *openGlContext)
{
  this->openGlContext = openGlContext;
  shaderLines.addVertexShader(OpenGLHelpers::translateVertexShaderToV3(vertexShaderLinesUniformColor));
  shaderLines.addFragmentShader(OpenGLHelpers::translateFragmentShaderToV3(fragmentShaderLinesUniformColor));
  shaderLines.link();
  
  positionAttribId = glGetAttribLocation(shaderLines.getProgramID(), "position");
  pvMatrixUniformId = glGetUniformLocation(shaderLines.getProgramID(), "pvMatrix");
  colorUniformId = glGetUniformLocation(shaderLines.getProgramID(), "color");

  glGenBuffers(1, &minMaxVBOId);
  glGenVertexArrays(1, &minMaxVAOId);
  glBindVertexArray(minMaxVAOId);
  glBindBuffer(GL_ARRAY_BUFFER, minMaxVBOId);
  glEnableVertexAttribArray(positionAttribId);
  glVertexAttribPointer(positionAttribId, 2, GL_FLOAT, GL_FALSE, sizeof(position), 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glDisableVertexAttribArray(positionAttribId);

  glGenBuffers(1, &rmsVBOId);
  glGenVertexArrays(1, &rmsVAOId);
  glBindVertexArray(rmsVAOId);
  glBindBuffer(GL_ARRAY_BUFFER, rmsVBOId);
  glEnableVertexAttribArray(positionAttribId);
  glVertexAttribPointer(positionAttribId, 2, GL_FLOAT, GL_FALSE, sizeof(position), 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glDisableVertexAttribArray(positionAttribId);
}

void WaveDisplay :: renderOpenGL()
{
  audio *buf = tfView->sbsms->getCache(0)->wave.getReadBuf();

  //printf("%ld\n",tfView->sbsms->getCache(0)->wave.nReadable());

  const float desktopScale = (float) openGlContext->getRenderingScale();
  Rectangle<int> r = getBounds();
  int x0 = desktopScale * r.getX();
  int y0 = desktopScale * (getParentComponent()->getHeight()-(r.getHeight()+r.getY()));
  int w = desktopScale * r.getWidth();
  int h = desktopScale * r.getHeight();
  glViewport(x0,y0,w,h);
  glScissor(x0,y0,w,h);

  TimeFrequencyOverlay *view = tfView->view;
  SampleCountType len = view->endSample - view->startSample;
  int channel = tfView->channel;

  shaderLines.use();
  glLineWidth(1);
  
  mat4_t pMatrix = mat4_create(NULL);
  mat4_ortho(0,w,0,h,0.1,1.0,pMatrix);
  mat4_t pvMatrix = mat4_create(NULL);
  mat4_t vMatrix = mat4_create(NULL);
  vec3_t v = vec3_create(NULL);
  mat4_identity(vMatrix);
  mat4_multiply(pMatrix,vMatrix,pvMatrix);
  v[0] = 0.0f;
  v[1] = 0.0f;
  v[2] = -0.5f;
  mat4_translate(pvMatrix,v,NULL);
  glUniformMatrix4fv(pvMatrixUniformId, 1, GL_FALSE, (const GLfloat*)pvMatrix);
  

  float samplesPerPixel = (float)len / (float)w;
  float mid = 0.5f * h;
    
  if(samplesPerPixel < 1.0f) {
    
  } else {
    vector<position> rmsPos;
    vector<position> minMaxPos;  

    int start = view->startSample;
    float endGroup = start + samplesPerPixel;
    int end = lrintf(endGroup);
    int lastMinY;
    int lastMaxY;
    for(int x=0; x<w; x++) {
      float theMin = 1.0f;
      float theMax = -1.0f;
      float sum = 0.0f;

      for(int k=start; k<end; k++) {
        float a = buf[k][channel];
        if(a > theMax) {
          theMax = a;
        }
        if(a < theMin) {
          theMin = a;
        }
        sum += a*a;
      }
      int n = end-start;
      float rms = sqrt(sum / n);
      start = end;      
      endGroup += samplesPerPixel;
      end = lrintf(endGroup);
      int minY = lrintf(theMin * mid + mid);
      int maxY = lrintf(theMax * mid + mid);
      int rmsPY = lrintf(rms * mid + mid);
      int rmsNY = lrintf(-rms * mid + mid);
      if(x > 0) {
        if(minY > lastMaxY) {
          minY = lastMaxY + 1;
        }
        if(maxY < lastMinY) {
          maxY = lastMinY - 1;
        }
      }      
      lastMinY = minY;
      lastMaxY = maxY;
      if(rmsPY > maxY-1) {
        rmsPY = maxY-1;
      }
      if(rmsNY < minY+1) {
        rmsNY = minY+1;
      }


      position p;
      p.x = x;

      p.y = minY;
      minMaxPos.push_back(p);
      p.y = maxY;
      minMaxPos.push_back(p);

      p.y = rmsNY;
      rmsPos.push_back(p);
      p.y = rmsPY;
      rmsPos.push_back(p);

      //printf("%d %g %g %g %d %d %d %d\n",x,theMin,theMax,rms,minY,maxY,rmsNY,rmsPY);
      
    }
    
    
    float color[4] = {0.6f,0.6f,1.0f,1.0f};
    glUniform4fv(colorUniformId, 1, color);
    glBindVertexArray(minMaxVAOId);
    glBindBuffer(GL_ARRAY_BUFFER, minMaxVBOId);
    glBufferData(GL_ARRAY_BUFFER, sizeof(position) * minMaxPos.size(), minMaxPos.data(), GL_STREAM_DRAW);
    glDrawArrays(GL_LINES, 0, minMaxPos.size());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    v[0] = 0.0f;
    v[1] = 0.0f;
    v[2] = 0.25f;
    mat4_translate(pvMatrix,v,NULL);
    glUniformMatrix4fv(pvMatrixUniformId, 1, GL_FALSE, (const GLfloat*)pvMatrix);

    float colorRms[4] = {0.4f,0.4f,0.8f,1.0f};
    glUniform4fv(colorUniformId, 1, colorRms);
    glBindVertexArray(rmsVAOId);
    glBindBuffer(GL_ARRAY_BUFFER, rmsVBOId);
    glBufferData(GL_ARRAY_BUFFER, sizeof(position) * rmsPos.size(), rmsPos.data(), GL_STREAM_DRAW);
    glDrawArrays(GL_LINES, 0, rmsPos.size());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
  }
  
}

void WaveDisplay :: openGLContextClosing()
{
}

bool TimeFrequencyView :: shouldAssign(int band, TrackPoint *tp, TrackIndexType index) 
{
  return (bSelectAll != (selected[band].find(index) != selected[band].end())) && bEnableBand[band];
}

bool TimeFrequencyView :: shouldScale(int band, TrackPoint *tp, TrackIndexType index, int c) 
{
  return (bAdjust && tp->bScale[c]) || bForceAdjust;
}


bool TimeFrequencyView :: shouldOnset(int band, TrackPoint *tp, TrackIndexType index)
{
  return bOnset;
}

TrackView :: TrackView(TimeFrequencyView *tfView)
  : shaderLines(tfView->openGLContext)
{
  this->tfView = tfView;
}

void TrackView :: paint(Graphics &g)
{
  pthread_mutex_lock(&tfView->glMutex);
  if(tfView->cursorSample >= 0) {
    int w = getWidth();
    int h = getHeight();
    TimeFrequencyOverlay *view = tfView->view;
    float x = view->xScale * (tfView->cursorSample - view->startSample);
    g.setColour(Colours::white.withAlpha(0.5f));
    g.drawLine(x,0,x,h);
  }
  pthread_mutex_unlock(&tfView->glMutex);
}


void TrackView :: newOpenGLContextCreated(OpenGLContext *openGlContext) 
{
  this->openGlContext = openGlContext;
  shaderLines.addVertexShader(OpenGLHelpers::translateVertexShaderToV3(vertexShaderLinesUniformColor));
  shaderLines.addFragmentShader(OpenGLHelpers::translateFragmentShaderToV3(fragmentShaderLinesUniformColor));
  shaderLines.link();

  positionAttribId = glGetAttribLocation(shaderLines.getProgramID(), "position");
  pvMatrixUniformId = glGetUniformLocation(shaderLines.getProgramID(), "pvMatrix");
  colorUniformId = glGetUniformLocation(shaderLines.getProgramID(), "color");

  glGenBuffers(1, &trackVBOId);
  glGenVertexArrays(1, &trackVAOId);
  glBindVertexArray(trackVAOId);
  glBindBuffer(GL_ARRAY_BUFFER, trackVBOId);
  glEnableVertexAttribArray(positionAttribId);
  glVertexAttribPointer(positionAttribId, 2, GL_FLOAT, GL_FALSE, sizeof(position), 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glDisableVertexAttribArray(positionAttribId);

  glGenBuffers(1, &trackScaleVBOId);
  glGenVertexArrays(1, &trackScaleVAOId);
  glBindVertexArray(trackScaleVAOId);
  glBindBuffer(GL_ARRAY_BUFFER, trackScaleVBOId);
  glEnableVertexAttribArray(positionAttribId);
  glVertexAttribPointer(positionAttribId, 2, GL_FLOAT, GL_FALSE, sizeof(position), 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glDisableVertexAttribArray(positionAttribId);

  glGenBuffers(1, &trackDtVBOId);
  glGenVertexArrays(1, &trackDtVAOId);
  glBindVertexArray(trackDtVAOId);
  glBindBuffer(GL_ARRAY_BUFFER, trackDtVBOId);
  glEnableVertexAttribArray(positionAttribId);
  glVertexAttribPointer(positionAttribId, 2, GL_FLOAT, GL_FALSE, sizeof(position), 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glDisableVertexAttribArray(positionAttribId);

  glGenBuffers(1, &trackDPhVBOId);
  glGenVertexArrays(1, &trackDPhVAOId);
  glBindVertexArray(trackDPhVAOId);
  glBindBuffer(GL_ARRAY_BUFFER, trackDPhVBOId);
  glEnableVertexAttribArray(positionAttribId);
  glVertexAttribPointer(positionAttribId, 2, GL_FLOAT, GL_FALSE, sizeof(position), 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glDisableVertexAttribArray(positionAttribId);

  glGenBuffers(1, &trackMeritVBOId);
  glGenVertexArrays(1, &trackMeritVAOId);
  glBindVertexArray(trackMeritVAOId);
  glBindBuffer(GL_ARRAY_BUFFER, trackMeritVBOId);
  glEnableVertexAttribArray(positionAttribId);
  glVertexAttribPointer(positionAttribId, 2, GL_FLOAT, GL_FALSE, sizeof(position), 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glDisableVertexAttribArray(positionAttribId);

  glGenBuffers(1, &totalMeritVBOId);
  glGenVertexArrays(1, &totalMeritVAOId);
  glBindVertexArray(totalMeritVAOId);
  glBindBuffer(GL_ARRAY_BUFFER, totalMeritVBOId);
  glEnableVertexAttribArray(positionAttribId);
  glVertexAttribPointer(positionAttribId, 2, GL_FLOAT, GL_FALSE, sizeof(position), 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glDisableVertexAttribArray(positionAttribId);

  glGenBuffers(1, &trackOnsetVBOId);
  glGenVertexArrays(1, &trackOnsetVAOId);
  glBindVertexArray(trackOnsetVAOId);
  glBindBuffer(GL_ARRAY_BUFFER, trackOnsetVBOId);
  glEnableVertexAttribArray(positionAttribId);
  glVertexAttribPointer(positionAttribId, 2, GL_FLOAT, GL_FALSE, sizeof(position), 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glDisableVertexAttribArray(positionAttribId);

  glGenBuffers(1, &trackOffsetVBOId);
  glGenVertexArrays(1, &trackOffsetVAOId);
  glBindVertexArray(trackOffsetVAOId);
  glBindBuffer(GL_ARRAY_BUFFER, trackOffsetVBOId);
  glEnableVertexAttribArray(positionAttribId);
  glVertexAttribPointer(positionAttribId, 2, GL_FLOAT, GL_FALSE, sizeof(position), 0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);
  glDisableVertexAttribArray(positionAttribId);
}

void TrackView :: openGLContextClosing()
{
  
}


void TrackView :: renderOpenGL()
{
  if(tfView->cursorTime < 0 || tfView->cursorF < 0.0f) {
    return;
  }

  const float desktopScale = (float) openGlContext->getRenderingScale();
  Rectangle<int> r = getBounds();
  int x0 = desktopScale * r.getX();
  int y0 = desktopScale * (getParentComponent()->getHeight()-(r.getHeight()+r.getY()));
  int w = desktopScale * r.getWidth();
  int h = desktopScale * r.getHeight();
  glViewport(x0,y0,w,h);
  glScissor(x0,y0,w,h);

  TimeFrequencyOverlay *view = tfView->view;
  int channel = tfView->channel;

  shaderLines.use();
  glLineWidth(1);
  

  mat4_t pMatrix = mat4_create(NULL);
  mat4_ortho(0,w,0,h,0.1,1.0,pMatrix);
  mat4_t pvMatrix = mat4_create(NULL);
  mat4_t vMatrix = mat4_create(NULL);
  vec3_t v = vec3_create(NULL);
  mat4_identity(vMatrix);

  v[0] = desktopScale * view->xScale;
  v[1] = 1.0f;
  v[2] = 1.0f;
  mat4_scale(vMatrix,v,NULL);

  mat4_multiply(pMatrix,vMatrix,pvMatrix);

  v[0] = 0.0f;
  v[1] = 0.0f;
  v[2] = -0.5f;
  mat4_translate(pvMatrix,v,NULL);
  glUniformMatrix4fv(pvMatrixUniformId, 1, GL_FALSE, (const GLfloat*)pvMatrix);

  if(tfView->selectedTrackSize) {
    // current selected track
    vector<position> trackPos;
    vector<position> trackScalePos;
    vector<position> trackDtPos;
    vector<position> trackDPhPos;
    vector<position> trackMeritPos;
    vector<position> totalMeritPos;
    track &tr = tfView->trackVBO[channel][tfView->selectedTrackBand][tfView->selectedTrackIndex];
    int k=0;
    for(int k=0; k<tr.vertices.size(); k++) {
      vertex &v = tr.vertices[k];
      TrackPoint *tp = tr.trackpoints[k];
      position pos;
      Cache *cache = tfView->sbsms->getCache(tfView->selectedTrackBand);
      pos.x = (v.x - view->startSample);
      pos.y = (1.0f + 0.33f * log10f(tp->m[channel]/MScale + 1e-8)) * h;
      trackPos.push_back(pos);
      pos.y = (1.0f + 0.33f * log10f(tp->s[channel] * tp->m[channel]/MScale + 1e-8)) * h;
      trackScalePos.push_back(pos);
      pos.y = 0.5f * h + 0.5f * h / PI * tp->dph[channel];
      trackDPhPos.push_back(pos);
      pos.y = 0.5f * h + 0.125f * h * tp->dt[channel] / cache->h;
      trackDtPos.push_back(pos);
      pos.y = h * max(0.0f,0.3f+log10f(tp->meritOn[channel].total)) * 0.3f;
      trackMeritPos.push_back(pos);
      pos.y = h * max(0.0f,0.3f+log10f(cache->meritOn[channel][tr.start+k])) * 0.3f;
      totalMeritPos.push_back(pos);
    }
    
    // magnitude
    float color[4] = {1.0f,1.0f,1.0f,1.0f};
    glUniform4fv(colorUniformId, 1, color);
    
    glBindVertexArray(trackVAOId);
    glBindBuffer(GL_ARRAY_BUFFER, trackVBOId);
    glBufferData(GL_ARRAY_BUFFER, sizeof(position) * trackPos.size(), trackPos.data(), GL_STREAM_DRAW);
    glDrawArrays(GL_LINE_STRIP, 0, trackPos.size());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);


    // magnitude scaled
    float colorScale[4] = {1.90f,0.2f,0.4f,1.0f};
    glUniform4fv(colorUniformId, 1, colorScale);
    
    glBindVertexArray(trackScaleVAOId);
    glBindBuffer(GL_ARRAY_BUFFER, trackScaleVBOId);
    glBufferData(GL_ARRAY_BUFFER, sizeof(position) * trackScalePos.size(), trackScalePos.data(), GL_STREAM_DRAW);
    glDrawArrays(GL_LINE_STRIP, 0, trackScalePos.size());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    
    // dt
    float colordt[4] = {0.5f,1.0f,1.0f,1.0f};
    glUniform4fv(colorUniformId, 1, colordt);
    glBindVertexArray(trackDtVAOId);
    glBindBuffer(GL_ARRAY_BUFFER, trackDtVBOId);
    glBufferData(GL_ARRAY_BUFFER, sizeof(position) * trackDtPos.size(), trackDtPos.data(), GL_STREAM_DRAW);
    glDrawArrays(GL_LINE_STRIP, 0, trackDtPos.size());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    // dph
    float colordph[4] = {0.0f,0.4f,1.0f,1.0f};
    glUniform4fv(colorUniformId, 1, colordph);
    glBindVertexArray(trackDPhVAOId);
    glBindBuffer(GL_ARRAY_BUFFER, trackDPhVBOId);
    glBufferData(GL_ARRAY_BUFFER, sizeof(position) * trackDPhPos.size(), trackDPhPos.data(), GL_STREAM_DRAW);
    glDrawArrays(GL_LINE_STRIP, 0, trackDPhPos.size());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    
    // track merit
    float colormerit[4] = {1.0f,1.0f,0.5f,1.0f};
    glUniform4fv(colorUniformId, 1, colormerit);
    glBindVertexArray(trackMeritVAOId);
    glBindBuffer(GL_ARRAY_BUFFER, trackMeritVBOId);
    glBufferData(GL_ARRAY_BUFFER, sizeof(position) * trackMeritPos.size(), trackMeritPos.data(), GL_STREAM_DRAW);
    glDrawArrays(GL_LINE_STRIP, 0, trackMeritPos.size());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);


    // total merit
    float colorTotalMerit[4] = {1.0f,0.7f,0.4f,1.0f};
    glUniform4fv(colorUniformId, 1, colorTotalMerit);
    glBindVertexArray(totalMeritVAOId);
    glBindBuffer(GL_ARRAY_BUFFER, totalMeritVBOId);
    glBufferData(GL_ARRAY_BUFFER, sizeof(position) * totalMeritPos.size(), totalMeritPos.data(), GL_STREAM_DRAW);
    glDrawArrays(GL_LINE_STRIP, 0, totalMeritPos.size());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    
    int grainSize = tfView->sbsms->getInputFrameSize() / tfView->sbsms->getQuality()->getBandGrainsPerFrame(tfView->selectedTrackBand);
    
    // onsets
    float onsetColor[4] = {0.0f,1.0f,1.0f,1.0f};
    glUniform4fv(colorUniformId, 1, onsetColor);
    vector<position> trackOnsets;
    for(set<TimeType>::iterator i = tfView->trackOnsets.begin(); i != tfView->trackOnsets.end(); ++i) {
      position pos;
      pos.x = (*i * grainSize - view->startSample) * view->xScale;
      pos.y = 0;
      trackOnsets.push_back(pos);
      pos.y = h;
      trackOnsets.push_back(pos);
    }
    glBindVertexArray(trackOnsetVAOId);
    glBindBuffer(GL_ARRAY_BUFFER, trackOnsetVBOId);
    glBufferData(GL_ARRAY_BUFFER, sizeof(position) * trackOnsets.size(), trackOnsets.data(), GL_STREAM_DRAW);
    glDrawArrays(GL_LINES, 0, trackOnsets.size());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    
    
    // offsets
    float offsetColor[4] = {1.0f,1.0f,0.0f,1.0f};
    glUniform4fv(colorUniformId, 1, offsetColor);
    vector<position> trackOffsets;
    for(set<TimeType>::iterator i = tfView->trackOffsets.begin(); i != tfView->trackOffsets.end(); ++i) {
      position pos;
      pos.x = (*i * grainSize - view->startSample) * view->xScale;
      pos.y = 0;
      trackOffsets.push_back(pos);
      pos.y = h;
      trackOffsets.push_back(pos);
    }
    glBindVertexArray(trackOffsetVAOId);
    glBindBuffer(GL_ARRAY_BUFFER, trackOffsetVBOId);
    glBufferData(GL_ARRAY_BUFFER, sizeof(position) * trackOffsets.size(), trackOffsets.data(), GL_STREAM_DRAW);
    glDrawArrays(GL_LINES, 0, trackOffsets.size());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
  }
}

SpectrumView :: SpectrumView(TimeFrequencyView *tfView, int which, float *color)
  : shaderLines(tfView->openGLContext)
{
  this->which = which;
  this->tfView = tfView;
  memcpy(this->color,color,4*sizeof(float));
}

void SpectrumView :: newOpenGLContextCreated(OpenGLContext *openGlContext) 
{
  this->openGlContext = openGlContext;
  shaderLines.addVertexShader(OpenGLHelpers::translateVertexShaderToV3(vertexShaderLinesUniformColor));
  shaderLines.addFragmentShader(OpenGLHelpers::translateFragmentShaderToV3(fragmentShaderLinesUniformColor));
  shaderLines.link();

  positionAttribId = glGetAttribLocation(shaderLines.getProgramID(), "position");
  pvMatrixUniformId = glGetUniformLocation(shaderLines.getProgramID(), "pvMatrix");
  colorUniformId = glGetUniformLocation(shaderLines.getProgramID(), "color");

  SBSMS *sbsms = tfView->sbsms;
  int bands = sbsms->getQuality()->params.bands;
  glGenBuffers(bands, graphVBOId);
  glGenVertexArrays(bands, graphVAOId);
  if(which == 0) {
    glGenBuffers(bands, cutsVBOId);
    glGenVertexArrays(bands, cutsVAOId);
  }

  // spectrum graphs
  for(int band=0; band<bands; band++) {
    glBindVertexArray(graphVAOId[band]);
    glBindBuffer(GL_ARRAY_BUFFER, graphVBOId[band]);
    
    glEnableVertexAttribArray(positionAttribId);
    glVertexAttribPointer(positionAttribId, 2, GL_FLOAT, GL_FALSE, sizeof(position), 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    glDisableVertexAttribArray(positionAttribId);

    if(which == 0) {
      glBindVertexArray(cutsVAOId[band]);
      glBindBuffer(GL_ARRAY_BUFFER, cutsVBOId[band]);
      glEnableVertexAttribArray(positionAttribId);
      glVertexAttribPointer(positionAttribId, 2, GL_FLOAT, GL_FALSE, sizeof(position), 0);
      glBindBuffer(GL_ARRAY_BUFFER, 0);
      glBindVertexArray(0);
      glDisableVertexAttribArray(positionAttribId);
    }
  }  
}

void SpectrumView :: openGLContextClosing()
{
  
}


void SpectrumView :: renderOpenGL()
{
  if(tfView->cursorTime < 0 || tfView->cursorF < 0.0f) {
    return;
  }

  const float desktopScale = (float) openGlContext->getRenderingScale();
  Rectangle<int> r = getBounds();
  int x0 = desktopScale * r.getX();
  int y0 = desktopScale * (getParentComponent()->getHeight()-(r.getHeight()+r.getY()));
  int w = desktopScale * r.getWidth();
  int h = desktopScale * r.getHeight();
  glViewport(x0,y0,w,h);
  glScissor(x0,y0,w,h);

  TimeFrequencyOverlay *view = tfView->view;
  SBSMS *sbsms = tfView->sbsms;
  int bands = sbsms->getQuality()->params.bands;
  int channel = tfView->channel;


  glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
  glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ZERO);

  shaderLines.use();
  glLineWidth(1);

  mat4_t pMatrix = mat4_create(NULL);
  mat4_ortho(0,w,0,h,0.1,1.0,pMatrix);
  mat4_t pvMatrix = mat4_create(NULL);
  mat4_t vMatrix = mat4_create(NULL);
  vec3_t v = vec3_create(NULL);
  mat4_identity(vMatrix);

  v[0] = 1.0f;
  v[1] = desktopScale * view->yScale;
  v[2] = 1.0f;
  mat4_scale(vMatrix,v,NULL);

  mat4_multiply(pMatrix,vMatrix,pvMatrix);
  v[0] = 0.0f;
  v[1] = 0.0f;
  v[2] = -0.5f;
  mat4_translate(pvMatrix,v,NULL);
  glUniformMatrix4fv(pvMatrixUniformId, 1, GL_FALSE, (const GLfloat*)pvMatrix);
  glUniform4fv(colorUniformId, 1, color);

  int bandtime = tfView->cursorTime;
  for(int band=0; band<bands; band++) {
    if(view->startF > tfView->topF[band]) break;
    if(tfView->bEnableBand[band] && view->endF >= tfView->botF[band]) {
      Cache *cache = sbsms->getCache(band);
      float *mag;
      if(which == 2) {
        mag = cache->mag2Cache[channel][bandtime];
      } else if(which == 1) {
        mag = cache->mag1Cache[channel].size()?cache->mag1Cache[channel][bandtime]:NULL;
      } else {
        mag = cache->magTrialCache[channel].size()?cache->magTrialCache[channel][bandtime]:NULL;
      }
      position graphData[512];
      if(mag) {
        float scale = -0.10f * w;        
        int kLo = cache->kLo;
        int kHi = cache->kHi;
        int minK = cache->minK;
        int maxK = cache->maxK;
        
        float y = (tfView->botF[band] - view->startF);
        float yStep = (tfView->topF[band] - tfView->botF[band])/ (kHi - kLo - 1);
        y += (minK - kLo) * yStep;
        for(int k=minK; k<maxK; k++) {
          float x0 = min((float)w,max(0.0f,scale * mag[k]));
          graphData[k].x = x0;
          graphData[k].y = y;
          y += yStep;          
        }
      
        int size = (maxK-minK);
        glBindVertexArray(graphVAOId[band]);
        glBindBuffer(GL_ARRAY_BUFFER, graphVBOId[band]);
        glBufferData(GL_ARRAY_BUFFER, sizeof(position) * size, graphData+minK, GL_STREAM_DRAW);
        glDrawArrays(GL_LINE_STRIP, 0, size);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);

        
        if(which == 0) {
          float color[4] = {1.0f,0.7f,0.8f,0.8f};
          glUniform4fv(colorUniformId, 1, color);
          vector<position> lines;
          list<int> &cuts = cache->cuts[channel][bandtime];
          float y0 = (tfView->botF[band] - view->startF);
          for(list<int>::iterator i = cuts.begin(); i != cuts.end(); ++i) {
            int k = *i;
            position pos;
            pos.x = w;
            pos.y = y0 + yStep * (k - kLo);            
            lines.push_back(pos);
            pos.x = min((float)w,max(0.0f,scale * mag[k]));
            lines.push_back(pos);
            //printf("%d %d %d %d %d %g\n",band,bandtime,minK,k,maxK,pos.x);
          }
          glBindVertexArray(cutsVAOId[band]);
          glBindBuffer(GL_ARRAY_BUFFER, cutsVBOId[band]);
          glBufferData(GL_ARRAY_BUFFER, sizeof(position) * lines.size(), lines.data(), GL_STREAM_DRAW);
          glDrawArrays(GL_LINES, 0, lines.size());
          glBindBuffer(GL_ARRAY_BUFFER, 0);
          glBindVertexArray(0);
        }
      }
    }
    bandtime /= sbsms->getQuality()->params.res[band];      
  }
}

void SpectrumView :: paint(Graphics &g)
{
  pthread_mutex_lock(&tfView->glMutex);
  if(tfView->cursorF > 0) {
    int w = getWidth();
    int h = getHeight();
    TimeFrequencyOverlay *view = tfView->view;
    float y = h-view->yScale * (tfView->cursorF - view->startF);
    g.setColour(Colours::white.withAlpha(0.5f));
    g.drawLine(0,y,w,y);
  }
  pthread_mutex_unlock(&tfView->glMutex);
}


enum {
  CommandNone = 0,
  CommandDeselect,
  CommandSelectNew,
  CommandSelectAppend
};

enum {
  None = 0,
  Left = 1,
  Right = 2,
  Both = 3
};

StatusView :: StatusView() 
{
  numBands = 0;
  pos = 0;
}
  
void StatusView :: setStatus(float freq, int *band, int *k, int *time, int numBands,SampleCountType pos) 
{
  this->freq = freq;
  this->pos = pos;
  this->numBands = numBands;
  for(int i=0; i<numBands; i++) {
    this->band[i] = band[i];
    this->time[i] = time[i];
    this->k[i] = k[i];
  }
}

void StatusView :: paint(Graphics &g) 
{
  if(numBands) {
    g.setColour(Colours::white);
    String str;
    const int textHeight = 16;
    int w = getWidth();
    int y = 2;
    str = String(freq,1) + String(" Hz, k") + String(band[0]) + String("=") +String(k[0]);
    if(numBands == 2) {
      str += String("k") + String(band[1]) + String("=") +String(k[1]);
    }
    g.drawText(str,2,y,w-3,textHeight,Justification::centredLeft,false);
    y += textHeight;
    str = String("t") + String(band[0]) + String("=") + String(time[0]);
    if(numBands == 2) {
      str += String(", t") + String(band[1]) + String("=") + String(time[1]);
    }
    g.drawText(str,2,y,w-3,textHeight,Justification::centredLeft,false);
    y += textHeight;
    str = String("t = ") + String(pos);
    g.drawText(str,2,y,w-3,textHeight,Justification::centredLeft,false);
  }
}
    
TimeFrequencyOverlay :: TimeFrequencyOverlay() 
{
  dragging = None;
  leftPos = -1;
  rightPos = -1;
  currPos = -1;
  startSample = -1;
  endSample = -1;
  selectionRect *= 0.0f;
}

void TimeFrequencyOverlay :: addStatusListener(StatusListener *listener)
{
  listeners.push_back(listener);
}

void TimeFrequencyOverlay :: resized() 
{
  rescale();
}

void TimeFrequencyOverlay :: rescale() 
{
  if(endSample > startSample) 
    xScale = (double)getWidth() / (double)(endSample - startSample);
  if(endF > startF)
    yScale = (double)getHeight() / (double)(endF - startF);
}

void TimeFrequencyOverlay :: paint(Graphics& g) 
{
  int h = getHeight();
  int w = getWidth();

  // draw cursors
  int x0;
  if(leftPos >= 0) {
    x0 = lrintf(xScale * (leftPos - startSample));
    if(x0 >= 0 && x0 < w) {
      g.setColour(Colours::white.withAlpha(0.5f));
      g.drawLine(x0,0,x0,h);
    }
  }

  if(rightPos >= 0) {
    int x1 = lrintf(xScale * (rightPos - startSample));
    if(x1 >= 0 && x1 < w) {
      g.setColour(Colours::white.withAlpha(0.5f));
      g.drawLine(x1,0,x1,h);
    }
    if(x0 < w && x1 >= 0) {
      if(x0 < 0) x0 = 0;
      if(x1 > w) x1 = w;
      g.setColour(Colours::white.withAlpha(0.15f));
      g.fillRect(x0,0,x1-x0+1,h);
    }
  }

  // play cursor
  if(currPos >= 0) {
    int xc = lrintf(xScale * (currPos - startSample));
    if(xc >= 0 && xc < w) {
      g.setColour(Colours::white.withAlpha(0.5f));
      g.drawLine(xc,0,xc,h);
    }
  }

  // selection 
  if(!selectionRect.isEmpty()) {
    g.setColour(Colours::green.withAlpha(0.25f));
    g.fillRect(selectionRect);
  }
}

const float mouseSelectDist = 2.0f;

void TimeFrequencyOverlay :: mouseMove(const MouseEvent& e)
{
  if(fabsf(e.x - xScale * (leftPos - startSample)) < mouseSelectDist) {
    setMouseCursor(MouseCursor::LeftRightResizeCursor);
  } else if(rightPos > 0 && fabsf(e.x - xScale * (rightPos - startSample)) < mouseSelectDist) {
    setMouseCursor(MouseCursor::LeftRightResizeCursor);
  } else {
    setMouseCursor(MouseCursor::NormalCursor);
  }
  float f = endF - e.y / yScale;
  SampleCountType pos = startSample + lrintf(e.x / xScale);
  for(vector<StatusListener*>::iterator i = listeners.begin(); i != listeners.end(); ++i) {
    (*i)->statusChanged(pos,f,CommandNone);
  }
}

void TimeFrequencyOverlay :: mouseExit(const MouseEvent& e)
{
  for(vector<StatusListener*>::iterator i = listeners.begin(); i != listeners.end(); ++i) {
    (*i)->statusChanged(-1,-1.0f,CommandNone);
  }
}

void TimeFrequencyOverlay :: mouseDown (const MouseEvent& e)
{
  if(e.mods.isCommandDown()) {
    selectAnchorPoint = e.position;
  } else {
    if(fabsf(e.x - xScale * (leftPos - startSample)) < mouseSelectDist) {
      setMouseCursor(MouseCursor::LeftRightResizeCursor);
      dragging = Left;
    } else if(rightPos > 0 && fabsf(e.x - xScale * (rightPos - startSample)) < mouseSelectDist) {
      setMouseCursor(MouseCursor::LeftRightResizeCursor);
      dragging = Right;
    } else {
      setMouseCursor(MouseCursor::NormalCursor);
      dragging = Both;
      anchorPos = startSample + lrintf(e.x / xScale);
      leftPos = anchorPos;
      rightPos = -1;
      repaint();
    }
  }
}

void TimeFrequencyOverlay :: mouseUp (const MouseEvent& e)
{
  if(e.mouseWasClicked()) {
    float f = endF - e.y / yScale;
    SampleCountType pos = startSample + lrintf(e.x / xScale);
    if(pos >= startSample && pos < endSample) {
      int command;
      if(e.mods.isCommandDown()) {
        command = CommandSelectAppend;
      } else if(e.mods.isCtrlDown()) {
        command = CommandSelectNew;
      } else {
        command = CommandNone;
      }
      for(vector<StatusListener*>::iterator i = listeners.begin(); i != listeners.end(); ++i) {
        (*i)->statusChanged(pos,f,command);
      }
    }
  }
  selectionRect *= 0.0f;
  dragging = None;
}

void TimeFrequencyOverlay :: mouseDrag (const MouseEvent& e)
{
  SampleCountType pos = startSample + lrintf(e.x / xScale);
  pos = max(startSample,pos);
  pos = min(endSample,pos);
  if(e.mods.isCommandDown()) {
    Point<float> selectPoint = e.position;
    selectionRect = Rectangle<float>(selectAnchorPoint,selectPoint);
    Rectangle<float> scaledRect = selectionRect/Point<float>(xScale,-yScale) + Point<float>(startSample,endF);
    int command;
    if(e.mods.isShiftDown()) {
      command = CommandSelectAppend;
    } else {
      command = CommandSelectNew;
    }
    for(vector<StatusListener*>::iterator i = listeners.begin(); i != listeners.end(); ++i) {
      (*i)->statusChanged(scaledRect,command);
    }
    repaint();
  } else {
    if(dragging == Left) {
      if(pos < rightPos) leftPos = pos;
      repaint();
    } else if(dragging == Right) {
      if(pos > leftPos) rightPos = pos;
      repaint();
    } else if(dragging == Both) {
      if(pos > anchorPos) {
        leftPos = anchorPos;
        rightPos = pos;
      } else {
        leftPos = pos;
        rightPos = anchorPos;
      }    
      repaint();
    }
  }
}

void TimeFrequencyView :: initColorMaps() 
{
  float r, g, b;
  int i;
  for (i=0; i<colorMapSize; i++) {
    float value = float(i)/colorMapSize;
    
    const int gsteps = 4;
    float gradient[gsteps + 1][3] = {
      {float(0.75), float(0.75), float(0.75)},    // lt gray
      {float(0.30), float(0.60), float(1.00)},    // lt blue
      {float(0.90), float(0.10), float(0.90)},    // violet
      {float(1.00), float(0.00), float(0.00)},    // red
      {float(1.00), float(1.00), float(1.00)}     // white
    };                        
          
    int left = int (value * gsteps);
    int right = (left == gsteps ? gsteps : left + 1);
    
    float rweight = (value * gsteps) - left;
    float lweight = 1.0 - rweight;
    
    r = (gradient[left][0] * lweight) + (gradient[right][0] * rweight);
    g = (gradient[left][1] * lweight) + (gradient[right][1] * rweight);
    b = (gradient[left][2] * lweight) + (gradient[right][2] * rweight);
  
    colorMap[i].r = (unsigned char) (255 * r);
    colorMap[i].g = (unsigned char) (255 * g);
    colorMap[i].b = (unsigned char) (255 * b);
    colorMap[i].a = (unsigned char) (255);
  }
}



// guranteed time >= t->start && time <= t->last
void TimeFrequencyView :: render(const SBSMSRenderChunk &i, Track *t, Debugger *dbg)
{
  TrackPoint *tp0 = i.time >= t->first?t->getTrackPoint(i.time):NULL;
  if(!tp0) return;
  float f0 = tp0->f;
  float m0[2] = {tp0->m[0], tp0->m[1] };
  bool bDummy = (tp0?tp0->flags&TrackDummy:false);
 
  for(int c=0; c<sbsms->getChannels(); c++) {
    float k = max(0.0f,min(1.0f,(m0[c]==0.0f?0.0f:(0.25f*log10f(m0[c]/MScale)+1.1f))));
    unsigned char col = 32 + lrintf(223.0f * k);
    vertex v;
    v.x = i.samplePos;
    v.y = f0;
    
    if(bDummy) {
      v.r = col/2;
      v.g = col;
      v.b = 0;
      v.a = col;
    } else {
      v.r = col;
      v.g = col/2;
      v.b = col;
      v.a = col;
    }

    if(tp0->flags & TrackOnset) {
      v.g = v.a;
    }

    if(tp0->flags & TrackOffset) {
      v.g = 0;
    }
    
    track &tr = trackVBO[c][i.band][t->index];
    if(tr.vertices.empty()) {
      tr.start = i.time;
    }
    tr.vertices.push_back(v);
    tr.trackpoints.push_back(tp0);
  }
}

//constructor
TimeFrequencyView :: TimeFrequencyView(SBSMS *sbsmsSrc) 
  : shaderLines(openGLContext),
    shaderTex(openGLContext)
{
  hScroll = new ScrollBar(false);
  hScroll->addListener(this);
  hScroll->setAutoHide(false);
  vScroll = new ScrollBar(true);
  vScroll->addListener(this);
  vScroll->setAutoHide(false);
  rateCtrl = new Slider();
  addAndMakeVisible(rateCtrl);
  rateCtrl->setColour(Slider::thumbColourId,Colours::white);
  rateCtrl->setColour(Slider::trackColourId,Colours::grey);
  rateCtrl->setColour(Slider::rotarySliderOutlineColourId,Colours::grey);
  rateCtrl->setColour(Slider::rotarySliderFillColourId,Colours::white);
   
  rateCtrl->setRange(0.0, 2.0f);
  rateCtrl->setSliderStyle(Slider::RotaryHorizontalVerticalDrag);
  rateCtrl->addListener(this);
	
  initColorMaps();
  setOpaque(true);
  trackViewHeight = 100;
  scrollSize = 12;
  xAxisHeight = 20;
  yAxisWidth = 40;
  controlsHeight = 74;
  spectrumWidth = 256;
  waveHeight = 128;
  cursorTime = -1;
  cursorSample = -1;
  cursorF = -1.0f;
  selectedTrackOffset = -1;
  selectedTrackSize = 0;
  bSelectedTrackLocked = false;
  selectedPoint = NULL;
  whichgram = 2;
  view = new TimeFrequencyOverlay();
  status = new StatusView();
  view->addStatusListener(this);
  bSelectAll = true;
  bAdjust = false;
  bOnset = true;
  bForceAdjust = false;
  timeAxisMode = TimeAxisSamples;
  freqAxisMode = FrequencyAxisHz;

  float color0[4] = {0.5f,0.5f,1.0f,1.0f};
  float color1[4] = {1.0f,1.0f,0.5f,1.0f};
  float color2[4] = {0.5f,1.0f,1.0f,1.0f};

  spectrum1 = new SpectrumView(this, 1, color1);
  spectrum2 = new SpectrumView(this, 2, color2);
  spectrumTrial = new SpectrumView(this, 0, color0);
  trackView = new TrackView(this);
  waveDisplay = new WaveDisplay(this);

  addAndMakeVisible(waveDisplay);
  addAndMakeVisible(hScroll);
  addAndMakeVisible(vScroll);
  addAndMakeVisible(view);
  addAndMakeVisible(status);
  addAndMakeVisible(spectrum2);
  addAndMakeVisible(spectrum1);
  addAndMakeVisible(spectrumTrial);
  addAndMakeVisible(trackView);

  pthread_mutex_init(&glMutex,NULL);
  openGLContext.setRenderer (this);
  openGLContext.attachTo (*this);
  
  this->sbsms = new SBSMS(sbsmsSrc,false);

  minF = 0.0f;
  maxF = M_PI;

  sbsms->addRenderer(this);
  totSamples = sbsms->getTotalSamples();
  view->startSample = 0;
  view->endSample = totSamples;
  view->leftPos = 0;
  view->rightPos = -1;
  view->frameSize = sbsms->getInputFrameSize();
  view->startF = minF;
  view->endF = maxF;

  bRenderSpectrum = true;
  bRenderTracks = true;


  SBSMSInterfaceVariableRate iface(totSamples);

  bands = sbsms->getQuality()->params.bands;
  channel = 0;
  for(int band=0; band<bands; band++) {
    bEnableBand[band] = true;
  }

  sbsms->reset(true);
  sbsms->seek(&iface,view->startSample);
  while(sbsms->renderFrameFromCache(&iface));
 
  // bitmap
  for(int which=0; which<3; which++) {
    const float a = 1.2*colorMapSize;
    const float b = colorMapSize / 8.0f;
    for(int c=0; c<sbsms->getChannels(); c++) {
      float bandF = TWOPI;
      for(int band=0; band<bands; band++) {
        Cache *cache = sbsms->getCache(band);
        
        vector<float*> &mags = which==0?cache->magTrialCache[c]:which==1?cache->mag1Cache[c]:cache->mag2Cache[c];

        int N = sbsms->getQuality()->params.N[band];
        int N2 = which==2?sbsms->getQuality()->params.N2[band]:sbsms->getQuality()->params.N1[band];

        float scale = 2.0f * 16.061113032124002f / (N2 * N);

        int kLo = cache->kLo;
        int kHi = cache->kHi;
        int minK = min(cache->minK,kLo);
        int maxK = max(cache->maxK,kHi);
        
        botF[band] = bandF * (float)cache->kLo / (float)N;
        topF[band] = bandF * ((float)cache->kHi-1) / (float)N;
        
        for(int t=0; t<mags.size(); t += 1024) {
          texQuad quad;
          texture &tex = quad.tex;
          tex.offsetX = t;
          tex.offsetY = 0;
          tex.width = kHi - kLo;
          tex.height = min(1024UL,mags.size()-t);
          tex.data = new pixel[tex.width*tex.height];
          pixel *data = tex.data;
          memset(data,0,sizeof(pixel)*tex.width*tex.height);
          
          for(int i=0; i<tex.height; i++) {
            float *mag = mags[t+i];
            for(int k=minK; k<=maxK; k++) {          
              mag[k] = log10f(1e-15f+mag[k]*scale);
            }
            
            for(int k=kLo; k<kHi; k++) {          
              int m = lrintf(a+b*mag[k]);
              m = min(colorMapSize-1,max(0,m));
              *data = colorMap[m];
              data++;
            }
          }
          spectrogram[which][c][band].push_back(quad);
        }
        bandF /= 2.0f;
      }
    }
  }

  audioDeviceManager.initialise(2, 2, 0, true, String::empty, 0);  
  sbsmsAudioSource = new SBSMSAudioSource(sbsms,this);
  audioSourcePlayer.setSource(sbsmsAudioSource);
  audioDeviceManager.addAudioCallback(&audioSourcePlayer);

  hScroll->setRangeLimits(0,totSamples,dontSendNotification);
  hScroll->setCurrentRange(view->startSample,view->endSample,dontSendNotification);
  vScroll->setRangeLimits(-maxF,-minF,dontSendNotification);
  vScroll->setCurrentRange(-view->endF,view->endF-view->startF,dontSendNotification);

  rateCtrl->setValue(1.0f);
}


TimeFrequencyView :: ~TimeFrequencyView()
{
  openGLContext.detach();
  delete hScroll;
  delete vScroll;
  delete view;
  delete status;
  delete trackView;
  delete spectrum1;
  delete spectrum2;
  delete spectrumTrial;
  delete rateCtrl;
  audioSourcePlayer.setSource(NULL);
  audioDeviceManager.removeAudioCallback(&audioSourcePlayer);
  delete sbsmsAudioSource;
  delete sbsms;
}

void TimeFrequencyView :: selectAll()
{
  pthread_mutex_lock(&glMutex);
  bSelectAll = true;
  for(int band=0; band<bands; band++) {
    selected[band].clear();
    selectedOffsets[band].clear();
    selectedCounts[band].clear();
  }
  trackOnsets.clear();
  trackOffsets.clear();
  bSelectedTrackLocked = false;
  pthread_mutex_unlock(&glMutex);
  repaint();
}


void TimeFrequencyView :: statusChanged(const Rectangle<float> &r, int command)
{
  pthread_mutex_lock(&glMutex);
  for(int band=0; band<bands; band++) {
    if(command == CommandSelectNew || command == CommandSelectAppend) {
      trackOnsets.clear();
      trackOffsets.clear();
    }
    if(command == CommandSelectNew) {
      selected[band].clear();
      selectedOffsets[band].clear();
      selectedCounts[band].clear();
    }
    if(bEnableBand[band]) {
      if(r.getY() >= botF[band] && r.getBottom() <= topF[band]) {
        Cache *cache = sbsms->getCache(band);
        int grainSize = sbsms->getInputFrameSize() / sbsms->getQuality()->getBandGrainsPerFrame(band);
        int time0 = max(0,(int)lrintf(r.getX() / grainSize));
        int time1 = min((int)cache->nTracksAtTime.size()-1,
                        (int)lrintf(r.getRight() / grainSize));
        for(int time=time0; time<=time1; time++) {        
          int nTracks = cache->nTracksAtTime[time];
          long cacheIndex = cache->indexAtTime[time];
          for(int k=0; k<nTracks; k++) {
            TrackPoint *tp = cache->trackPoints[cacheIndex+k];
            if(tp->f > r.getBottom() && tp->f < r.getY()) {
              int trackIndex = cache->trackIndex[cacheIndex+k];
              if(selected[band].insert(trackIndex).second) {
                unsigned int glIndex = trackStartIndex[channel][band][trackIndex];
                selectedOffsets[band].push_back(allOffsets[channel][band][glIndex]);
                selectedCounts[band].push_back(allCounts[channel][band][glIndex]);
              }
            }
          }
        }
      }
    }
  }
  pthread_mutex_unlock(&glMutex);
}


void TimeFrequencyView :: statusChanged(SampleCountType pos, float f, int command)
{
  pthread_mutex_lock(&glMutex);
  if(pos < 0 || f < 0.0f) {
    cursorTime = -1;
    cursorSample = -1;
    cursorF = -1.0f;
    status->setStatus(f,NULL,NULL,NULL,0,0);
    selectedPoint = NULL;
    if(!bSelectedTrackLocked) {
      selectedTrackSize = 0;
      selectedTrackOffset = -1;
    }
  } else {
    cursorF = f;
    cursorSample = pos;

    int grainSize = sbsms->getInputFrameSize() / sbsms->getQuality()->getBandGrainsPerFrame(0);
    int time = pos / grainSize;
    if((pos % grainSize) > grainSize/2) time++;
    cursorTime = time;

    int whichBands[2];
    int whichK[2];
    int whichTimes[2];
    int nBands = 0;
    TrackPoint *mintp = NULL;
    float mindf = TWOPI;
    int minband = -1;
    int mintime = -1;
    TrackIndexType minTrackIndex = 0;
    float bandF = TWOPI;
    for(int band=0; band<bands; band++) {
      if(bEnableBand[band]) {
        if(f >= botF[band] && f <= topF[band]) {
          Cache *cache = sbsms->getCache(band);
          int nTracks = cache->nTracksAtTime[time];
          long index = cache->indexAtTime[time];
          for(int i = 0; i < nTracks; i++) {
            TrackPoint *tp = cache->trackPoints[index+i];
            float df = fabsf(f - tp->f);
            if(df < mindf) {
              mindf = df;
              mintp = tp;
              minband = band;
              minTrackIndex = cache->trackIndex[index+i];
              mintime = time;
            }
          }
          int k = lrintf(f / bandF * sbsms->getQuality()->params.N[band]);
          whichBands[nBands] = band;
          whichK[nBands] = k;
          whichTimes[nBands] = time;
          nBands++;
        }
      }
      bandF /= 2.0f;
      time /= sbsms->getQuality()->params.res[band];
    }

    // selected track
    if(mintp) {
      if(command == CommandSelectNew || !bSelectedTrackLocked) {
        unsigned int glIndex = trackStartIndex[channel][minband][minTrackIndex];
        selectedTrackIndex = minTrackIndex;
        selectedTrackTime = mintime;
        selectedTrackSize = allCounts[channel][minband][glIndex];
        selectedTrackOffset = allOffsets[channel][minband][glIndex];
        selectedTrackBand = minband;
      }

      selectedPoint = mintp;
      selectedPointBand = minband;
      selectedPointTime = mintime;
      selectedPointTrackIndex = minTrackIndex;

      if(command == CommandSelectAppend || command == CommandSelectNew) {
        if(command == CommandSelectNew) {
          trackOnsets.clear();
          trackOffsets.clear();
          for(int band=0; band<bands; band++) {
            selected[band].clear();
            selectedOffsets[band].clear();
            selectedCounts[band].clear();
          }
        }
        if(selected[minband].insert(minTrackIndex).second) {
          unsigned int glIndex = trackStartIndex[channel][minband][minTrackIndex];
          selectedOffsets[minband].push_back(allOffsets[channel][minband][glIndex]);
          selectedCounts[minband].push_back(allCounts[channel][minband][glIndex]);
        } else {
          selected[minband].erase(minTrackIndex);
          unsigned int glIndex = trackStartIndex[channel][minband][minTrackIndex];
          selectedOffsets[minband].erase(remove(selectedOffsets[minband].begin(),
                                                selectedOffsets[minband].end(),
                                                allOffsets[channel][minband][glIndex]),
                                         selectedOffsets[minband].end());
          selectedCounts[minband].erase(remove(selectedCounts[minband].begin(),
                                               selectedCounts[minband].end(),
                                               allCounts[channel][minband][glIndex]),
                                        selectedCounts[minband].end()
                                        );
        }
        if(command == CommandSelectNew) {
          bSelectedTrackLocked = true;
        }
      }
    }
    status->setStatus(f*44100.0f/TWOPI,whichBands,whichK,whichTimes,nBands,pos);
  }
  pthread_mutex_unlock(&glMutex);
  repaint();
}

void TimeFrequencyView :: paint(Graphics &g)
{
  pthread_mutex_lock(&glMutex);
  int h = getHeight()-scrollSize-xAxisHeight-controlsHeight-trackViewHeight-waveHeight;
  int w = getWidth()-scrollSize-yAxisWidth-spectrumWidth;
  int x0 = yAxisWidth + spectrumWidth;
  int y = getHeight()-xAxisHeight-controlsHeight-trackViewHeight;
  int yAxisTop = scrollSize + waveHeight;

  g.setColour(Colours::white);
  // time axis
  const int tickTextWidth = 80;
  const int tickTextHeight = 32;
  const int tickTextMinSep = 10;
  int numTicks = lrintf(w/(tickTextWidth+tickTextMinSep)) + 1;
  SampleCountType rangeT = view->endSample - view->startSample;
  int lasttime = -1;
  for(int i=0; i<numTicks; i++) {
    SampleCountType pos = view->startSample + (i * rangeT) / (numTicks);
    int time = pos / view->frameSize;
    if(lasttime==time) continue;
    if(pos % view->frameSize > view->frameSize/2) time++;
    int x = x0 + lrintf(view->xScale * ((time * view->frameSize) - view->startSample));
    g.drawLine(x,y,x,y+2);

    String str;
    if(timeAxisMode == TimeAxisMSM) {
      float t = pos/44100.0f;
      int min = lrintf(floor(t/60.0f));
      if(min > 0) {
        str = String(min) + String(":");
      }
      t -= 60.0f * min;
      int sec = lrintf(floor(t));
      str += String(sec) + String(".");
      t -= sec;
      int ms = lrintf(t*1000.0f);
      str += String(ms).paddedLeft('0',3);
    } else if(timeAxisMode == TimeAxisSamples) {
      str = String(time*view->frameSize);
    } else if(timeAxisMode == TimeAxisGrains) {
      str = String(time);
    }
    g.drawText(str,x-tickTextWidth/2,y+4,tickTextWidth,tickTextHeight,Justification::centredTop,false);
    lasttime = time;
  }

  // frequency axis
  numTicks = lrintf(h/tickTextHeight) + 1;
  float rangeF = view->endF - view->startF;
  for(int i=0; i<numTicks; i++) {
    float f = view->startF + i * rangeF / numTicks;
    String str;    
    if(freqAxisMode == FrequencyAxisHz) {
      int freq = lrintf(f * 44100.0/TWOPI);
      if(freq >= 1000) {
        int thousands = freq/1000;
        int hundreds = lrintf(freq*.01)%10;
        str = String(thousands) + String(".") + String(hundreds) + String("k");
        f = (thousands * 1000 + hundreds * 100) / 44100.0f * TWOPI;
      } else {
        f = freq / 44100.0f * TWOPI;
        str = String(freq);
      }
    } else if(freqAxisMode == FrequencyAxis2PI) {
      str = String(f,3);
      f = str.getFloatValue();
    }
    if(f >= view->startF && f <= view->endF) {
      g.drawLine(x0-2,y,x0,y);
      y = yAxisTop + lrintf(view->yScale * (view->endF - f));
      g.drawText(str,x0-tickTextWidth-4,y-tickTextHeight/2,tickTextWidth,tickTextHeight,Justification::centredRight,false);
    }
  }

  // settings
  String str = bForceAdjust?String("Scaling forced"):(bAdjust?String("Scaling on") : String("Scaling off"));
  g.drawText(str,rateCtrl->getBounds().getRight()+8,getHeight()-controlsHeight+3,80,tickTextHeight,Justification::centredLeft,false);

  str = String("Spectrogram ") + String(whichgram) + String("\nchannel ") + String(channel) + String("\nonset ") + (bOnset?String("on"):String("off"));
  g.drawFittedText(str,rateCtrl->getBounds().getRight()+84,getHeight()-controlsHeight+3,96,tickTextHeight,Justification::centredLeft,2);

  // trackpoint data
  if(selectedPoint) {
    str = String("f=") + String(selectedPoint->f);
    str += String("\nm=") + String(selectedPoint->m[channel]/MScale) + String(" * ") + String(selectedPoint->s[channel]);
    str += String("\ndt=") + String(selectedPoint->dt[channel]) + String(" / ") + String(sbsms->getCache(selectedPointBand)->h) + String(" / ") + String(sbsms->getQuality()->params.N2[selectedPointBand]);
    str += String("\nband=") + String(selectedPointBand);
    g.drawFittedText(str,rateCtrl->getBounds().getRight()+200,getHeight()-controlsHeight+3,200,tickTextHeight,Justification::topLeft,6);


    str = String("merit=") + String(selectedPoint->meritOn[channel].total,4) + " / " + String(selectedPoint->meritOff[channel],4);
    if(selectedPoint->bScale[channel]) {
      str += String(" scaled");
    }
    str += String("\n(e0,e1,mean)=(") + String(selectedPoint->meritOn[channel].e0);
    str += String(", ") + String(selectedPoint->meritOn[channel].e1);
    str += String(", ") + String(selectedPoint->meritOn[channel].emean);
    str += String(")");
    str += String("\n(j,jmean,dbs)=(") + String(selectedPoint->meritOn[channel].j);
    str += String(", ") + String(selectedPoint->meritOn[channel].jmean);
    str += String(", ") + String(selectedPoint->meritOn[channel].dbs);
    str += String(")");
    str += String("\ndph=") + String(selectedPoint->dph[channel]);
    str += String(", ph=") + String(selectedPoint->ph[channel]);
    g.drawFittedText(str,rateCtrl->getBounds().getRight()+400,getHeight()-controlsHeight+3,330,tickTextHeight,Justification::topLeft,6);
    
  }


  pthread_mutex_unlock(&glMutex);
}
    
void TimeFrequencyView :: renderOpenGL()
{
  pthread_mutex_lock(&glMutex);
  OpenGLHelpers::clear (Colours::black);
  
  const float desktopScale = (float) openGLContext.getRenderingScale();
  int w = desktopScale * (getWidth()-scrollSize-yAxisWidth-spectrumWidth);
  int h = desktopScale * (getHeight()-scrollSize-xAxisHeight-controlsHeight-trackViewHeight-waveHeight);
  int x0 = desktopScale * (yAxisWidth+spectrumWidth);
  int y0 = desktopScale * (xAxisHeight+controlsHeight+trackViewHeight);
  glViewport(x0, y0, w, h);
  glScissor(x0, y0, w, h);

  glEnable(GL_SCISSOR_TEST);
  glEnable (GL_DEPTH_TEST);
  glDepthFunc (GL_LESS);
  glEnable (GL_BLEND);
  glBlendEquationSeparate(GL_FUNC_ADD, GL_FUNC_ADD);
  glBlendFuncSeparate(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, GL_ONE, GL_ZERO);

  mat4_t pMatrix = mat4_create(NULL);
  mat4_ortho(0,w,0,h,0.1,1.0,pMatrix);
  mat4_t pvMatrix = mat4_create(NULL);
  mat4_t vMatrix = mat4_create(NULL);
  vec3_t v = vec3_create(NULL);
  mat4_identity(vMatrix);
  
  v[0] = desktopScale * view->xScale;
  v[1] = desktopScale * view->yScale;
  v[2] = 1.0f;
  mat4_scale(vMatrix,v,NULL);
  
  v[0] = -view->startSample;
  v[1] = -view->startF;
  v[2] = 0.0f;
  mat4_translate(vMatrix,v,NULL);
  
  mat4_multiply(pMatrix,vMatrix,pvMatrix);

  v[0] = 0.0f;
  v[1] = 0.0f;
  v[2] = -0.8f;
  mat4_translate(pvMatrix,v,NULL);

  if(bRenderSpectrum) {
  shaderTex.use();
  glUniformMatrix4fv(pvMatrixUniformId_Tex, 1, GL_FALSE, (const GLfloat*)pvMatrix);
  glUniform1i(textureMapUniformId_Tex,0);
  

  for(int band=0; band<bands; band++) {
    if(bEnableBand[band])
      for(int i=0; i<spectrogram[whichgram][channel][band].size(); i++) {
        texQuad &quad = spectrogram[whichgram][channel][band][i];
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D,quad.texId);
        glBindVertexArray(quad.vaoId);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        glBindTexture(GL_TEXTURE_2D,0);
        glBindVertexArray(0);
      }
  }

  }


  if(bRenderTracks) {
  shaderLines.use();
  glLineWidth(1);

  float colorScale[4] = {1.0f,1.0f,1.0f,1.0f};
  glUniform4fv(colorScaleUniformId_Lines, 1, colorScale);

  v[0] = 0.0f;
  v[1] = 0.0f;
  v[2] = 0.5f;
  mat4_translate(pvMatrix,v,NULL);

  glUniformMatrix4fv(pvMatrixUniformId_Tex, 1, GL_FALSE, (const GLfloat*)pvMatrix);
;

 // all tracks
 for(int band=0; band<bands; band++) {
   if(bEnableBand[band]) {
     glBindVertexArray(vaoId[channel][band]);
     glMultiDrawArrays(GL_LINE_STRIP, allOffsets[channel][band].data(), allCounts[channel][band].data(), allCounts[channel][band].size());
     glBindVertexArray(0);
   }
 }

  // all points
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
  float scaleF = 30.0f;
  float scaleT = 1.0f;
  for(int band=0; band<bands; band++) {
    if(bEnableBand[band]) {
      int pointSize = max(1,
                          min(8,
                              min((int)lrintf(w/(scaleT*(view->endSample-view->startSample)/(view->frameSize/sbsms->getQuality()->getBandGrainsPerFrame(band)))),
                                  (int)lrintf(h/(scaleF*(view->endF-view->startF))))));
      glPointSize(pointSize);      
      glBindVertexArray(vaoId[channel][band]);
      glMultiDrawArrays(GL_POINTS, allOffsets[channel][band].data(), allCounts[channel][band].data(), allCounts[channel][band].size());
      glBindVertexArray(0);
    }
    scaleF *= 2.0f;
  }

 // selected tracks
  float colorScaleSelected[4] = {1.25f,1.5f,0.75f,1.0f};
  glUniform4fv(colorScaleUniformId_Lines, 1, colorScaleSelected);
  for(int band=0; band<bands; band++) {
    if(bEnableBand[band]) {
      glLineWidth(8);
      glBindVertexArray(vaoId[channel][band]);
      glMultiDrawArrays(GL_LINE_STRIP, selectedOffsets[band].data(), selectedCounts[band].data(), selectedCounts[band].size());
      glBindVertexArray(0);
    }
  }

  // current selected track
  if(selectedTrackSize) {
    if(bEnableBand[selectedTrackBand]) {
      float colorScaleCurrent[4] = {1.5f,2.0f,1.0f,1.0f};
      glUniform4fv(colorScaleUniformId_Lines, 1, colorScaleCurrent);
      glLineWidth(6);
      glBindVertexArray(vaoId[channel][selectedTrackBand]);
      glDrawArrays(GL_LINE_STRIP, selectedTrackOffset, selectedTrackSize);
      glBindVertexArray(0);
    }
  }

  
  // selected point
  if(selectedPoint && bEnableBand[selectedPointBand]) {
      float scaleF = 20.0f;
      float scaleT = 0.4f;
      int pointSize = max(2,
                          min(12,
                              min((int)lrintf(w/(scaleT*(view->endSample-view->startSample)/(view->frameSize/sbsms->getQuality()->getBandGrainsPerFrame(selectedPointBand)))),
                                  (int)lrintf(h/(scaleF*(view->endF-view->startF))))));
      TimeType startTime = trackStartTime[channel][selectedPointBand][selectedPointTrackIndex];
      glPointSize(pointSize);
      glBindVertexArray(vaoId[channel][selectedPointBand]);
      glDrawArrays(GL_POINTS, selectedTrackOffset+selectedPointTime-startTime, 1);
      glBindVertexArray(0);
  }
  
  }

  spectrum1->renderOpenGL();
  spectrum2->renderOpenGL();
  spectrumTrial->renderOpenGL();
  trackView->renderOpenGL();
  waveDisplay->renderOpenGL();

  GLenum error = glGetError();
  if(error) {
    printf("%d\n",error);
    abort();
  }

  glUseProgram(0);
  glDisable(GL_POINT_SMOOTH);
  glDisable(GL_BLEND);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_SCISSOR_TEST);

  pthread_mutex_unlock(&glMutex); 
}


void TimeFrequencyView :: newOpenGLContextCreated() 
{
  shaderLines.addVertexShader(OpenGLHelpers::translateVertexShaderToV3(vertexShaderLines));
  shaderLines.addFragmentShader(OpenGLHelpers::translateFragmentShaderToV3(fragmentShaderLines));
  shaderLines.link();

  shaderTex.addVertexShader(OpenGLHelpers::translateVertexShaderToV3(vertexShaderTex));
  shaderTex.addFragmentShader(OpenGLHelpers::translateFragmentShaderToV3(fragmentShaderTex));
  shaderTex.link();

  positionAttribId_Lines = glGetAttribLocation(shaderLines.getProgramID(), "position");
  colorAttribId_Lines = glGetAttribLocation(shaderLines.getProgramID(), "sourceColour");
  pvMatrixUniformId_Lines = glGetUniformLocation(shaderLines.getProgramID(), "pvMatrix");
  colorScaleUniformId_Lines = glGetUniformLocation(shaderLines.getProgramID(), "colorScale");

  positionAttribId_Tex = glGetAttribLocation(shaderTex.getProgramID(), "position");
  textureCoordAttribId_Tex = glGetAttribLocation(shaderTex.getProgramID(), "textureCoordIn");
  pvMatrixUniformId_Tex = glGetUniformLocation(shaderTex.getProgramID(), "pvMatrix");
textureMapUniformId_Tex = glGetUniformLocation(shaderTex.getProgramID(), "textureMap");

  // tracks
  for(int c=0; c<sbsms->getChannels(); c++) {
    glGenBuffers(bands, vboId[c]);
    glGenVertexArrays(bands, vaoId[c]);

    for(int band=0; band<bands; band++) {
      int start = 0;  
      unsigned int ntracks = 0;
      for(map<TrackIndexType, track >::iterator i = trackVBO[c][band].begin(); i != trackVBO[c][band].end(); ++i) {
        allVBO[c][band].insert(allVBO[c][band].end(),i->second.vertices.begin(),i->second.vertices.end());
        allCounts[c][band].push_back(i->second.vertices.size());
        allOffsets[c][band].push_back(start);
        trackStartIndex[c][band][i->first] = ntracks;
        trackStartTime[c][band][i->first] = i->second.start;
        ntracks++;
        start += i->second.vertices.size();
      }
      
      glBindVertexArray(vaoId[c][band]);
      glBindBuffer(GL_ARRAY_BUFFER, vboId[c][band]);
      glBufferData(GL_ARRAY_BUFFER, sizeof(vertex) * allVBO[c][band].size(), allVBO[c][band].data(), GL_STATIC_DRAW);
      glEnableVertexAttribArray(positionAttribId_Lines);
      glEnableVertexAttribArray(colorAttribId_Lines);
      glVertexAttribPointer(positionAttribId_Lines, 2, GL_FLOAT, GL_FALSE, sizeof(vertex), 0);
      glVertexAttribPointer(colorAttribId_Lines, 4, GL_UNSIGNED_BYTE, GL_TRUE, sizeof(vertex),   (void*)(sizeof( float ) * 2));
      glBindBuffer(GL_ARRAY_BUFFER, 0);
      glBindVertexArray(0);
      glDisableVertexAttribArray(positionAttribId_Lines);
      glDisableVertexAttribArray(colorAttribId_Lines);

    }
  }
  
  // spectogram textures
  for(int which=0; which<3; which++) {
    for(int c=0; c<sbsms->getChannels(); c++) {
      for(int band=0; band<bands; band++) {
        for(int i=0; i<spectrogram[which][c][band].size(); i++) {
          texQuad &quad = spectrogram[which][c][band][i];
          float q = sbsms->getInputFrameSize() / sbsms->getQuality()->getBandGrainsPerFrame(band);
          texture &tex = quad.tex;
          float x0 = tex.offsetX * q;
          float x1 = x0 + tex.height * q;
          float coords[16] = {x0, botF[band], 0.0f,0.0f,
                              x0, topF[band], 1.0f,0.0f,
                              x1, botF[band], 0.0f,1.0f,
                              x1, topF[band], 1.0f,1.0f};
        
          
          glGenTextures(1,&quad.texId);
          glGenVertexArrays(1,&quad.vaoId);
          glGenBuffers(1,&quad.vboId);

          glBindTexture(GL_TEXTURE_2D,quad.texId);      
          glPixelStorei (GL_UNPACK_ALIGNMENT, 4);
          glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tex.width, tex.height, 0, GL_RGBA, GL_UNSIGNED_BYTE, tex.data);
          
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
          glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
          
          glBindTexture(GL_TEXTURE_2D, 0);
          glBindVertexArray(quad.vaoId);
          glBindBuffer(GL_ARRAY_BUFFER, quad.vboId);
          glEnableVertexAttribArray(positionAttribId_Tex);
          glEnableVertexAttribArray(textureCoordAttribId_Tex);
          glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 16, coords, GL_STATIC_DRAW);
          glVertexAttribPointer(positionAttribId_Tex, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), 0);
          glVertexAttribPointer(textureCoordAttribId_Tex, 2, GL_FLOAT, GL_FALSE, 4*sizeof(float), (void*)(sizeof( float ) * 2));
          glBindBuffer(GL_ARRAY_BUFFER, 0);
          glBindVertexArray(0);
          glDisableVertexAttribArray(positionAttribId_Tex);
          glDisableVertexAttribArray(textureCoordAttribId_Tex);
        }
      }
    }
  }
  spectrum1->newOpenGLContextCreated(&openGLContext);
  spectrum2->newOpenGLContextCreated(&openGLContext);
  spectrumTrial->newOpenGLContextCreated(&openGLContext);
  trackView->newOpenGLContextCreated(&openGLContext);
  waveDisplay->newOpenGLContextCreated(&openGLContext);
}

void TimeFrequencyView :: openGLContextClosing()
{
  spectrum1->openGLContextClosing();
  spectrum2->openGLContextClosing();
  spectrumTrial->openGLContextClosing();
  trackView->openGLContextClosing();
  waveDisplay->openGLContextClosing();
}

void TimeFrequencyView :: markTrackOnset()
{
  if(selectedTrackSize) {
    int grainSize = sbsms->getInputFrameSize() / sbsms->getQuality()->getBandGrainsPerFrame(selectedTrackBand);
    TimeType t = cursorSample / grainSize;
    if(cursorSample - t*grainSize > grainSize/2) {
      t++;
    }
    if(!trackOnsets.insert(t).second) {
      trackOnsets.erase(t);
    }
  }
}

void TimeFrequencyView :: markTrackOffset()
{
  if(selectedTrackSize) {
    int grainSize = sbsms->getInputFrameSize() / sbsms->getQuality()->getBandGrainsPerFrame(selectedTrackBand);
    TimeType t = cursorSample / grainSize;
    if(cursorSample - t*grainSize > grainSize/2) {
      t++;
    }
    if(!trackOffsets.insert(t).second) {
      trackOffsets.erase(t);
    }
  }
}

void TimeFrequencyView :: outputSelectedTrack()
{
  if(selectedTrackSize) {
    FILE * fp = fopen("traintracks","a");
    track &tr = trackVBO[channel][selectedTrackBand][selectedTrackIndex];
    int grainSize = sbsms->getInputFrameSize() / sbsms->getQuality()->getBandGrainsPerFrame(selectedTrackBand);
    long t = tr.vertices[0].x / grainSize;
    for(vector<TrackPoint*>::iterator i = tr.trackpoints.begin(); i != tr.trackpoints.end(); ++i) {
      TrackPoint *tp = *i;
      int s;
      if(trackOnsets.find(t) != trackOnsets.end()) {
        s = 1;
      } else if(trackOffsets.find(t) != trackOffsets.end()) {
        s = -1;
      } else {
        s = 0;
      }
      printf("%d %d %g %g %g\n",s,selectedTrackBand,tp->f,tp->m[0]/MScale,tp->dt[0]);
      fprintf(fp,"%d %d %g %g %g\n",s,selectedTrackBand,tp->f,tp->m[0]/MScale,tp->dt[0]);
      t++;
    }
    printf("\n");
    fprintf(fp,"\n");
    fclose(fp);
  }  
}

void TimeFrequencyView :: saveSelection()
{
  if(!(view->leftPos >= 0 && view->rightPos <= totSamples && view->leftPos < view->rightPos)) return;
  
  File lastFile = File::getCurrentWorkingDirectory().getChildFile("out.wav").getFullPathName();
  FileChooser myChooser ("Save file...",
                         lastFile,
                         "*.wav",
                         false);
  if(myChooser.browseForFileToSave(false)) {
    File file = myChooser.getResult();
    const char *out = file.getFullPathName().getCharPointer();
    int channels = sbsms->getChannels();
    PcmWriter writer(out, view->rightPos-view->leftPos, 44100, channels);

    SBSMS sbsms(this->sbsms,true);
    SBSMSInterfaceVariableRate iface(sbsms.getTotalSamples());
    sbsms.reset(false);
    sbsms.seek(&iface,view->leftPos);
    sbsms.setLeftPos(view->leftPos);
    sbsms.setRightPos(view->rightPos);
    int  frameSize = sbsms.getInputFrameSize();
    audio *abuf = new audio[frameSize];
    float *fbuf = new float[frameSize<<1];
    long nWrite = -1;
    while(nWrite) {
      nWrite = sbsms.synthFromCache(&iface,abuf,frameSize,this);
      if(nWrite) {
        if(channels==1) {
          for(int k=0;k<nWrite;k++) {	
            fbuf[k] = abuf[k][0];
          }
        } else if(channels==2) {
          for(int k=0;k<nWrite;k++) {
            int k2 = k<<1;
            fbuf[k2] = abuf[k][0];
            fbuf[k2+1] = abuf[k][1];
          }
        }
        writer.write(fbuf,nWrite);
      }      
    }
    delete [] abuf;
    delete fbuf;    
    writer.close();
  }
}

bool TimeFrequencyView :: keyPressed(const KeyPress &key, Component * originatingComponent)
{
  switch(key.getTextCharacter()) {
  case 'h':
    // zoom in t
    zoomT(true);
    setHScrollRange();
    repaint();
    break;
  case 'g':
    // zoom out t
    zoomT(false);
    setHScrollRange();
    repaint();
    break;
  case 'b':
    // zoom in f
    zoomF(true);
    setVScrollRange();
    repaint();
    break;
  case 'v':
    // zoom out f
    zoomF(false);
    setVScrollRange();
    repaint();
    break;
  case 'o':
    bOnset = !bOnset;
    repaint();
    break;
  case 's':
    bRenderSpectrum = !bRenderSpectrum;
    repaint();
    break;
  case 't':
    bRenderTracks = !bRenderTracks;
    repaint();
    break;
  case 'w':
    whichgram = (whichgram + 1)%3;
    repaint();
    break;
  case 'c':
    if(sbsms && sbsms->getChannels() == 2) {
      channel = 1 - channel;
      // invalidate current track
      pthread_mutex_lock(&glMutex);
      selectedTrackSize = 0;
      pthread_mutex_unlock(&glMutex);
      statusChanged(cursorSample,cursorF,CommandNone);
      repaint();
    }
    break;
  case '0':
    bEnableBand[0] = !bEnableBand[0];
    repaint();
    break;
  case '1':
    bEnableBand[1] = !bEnableBand[1];
    repaint();
    break;
  case '2':
    bEnableBand[2] = !bEnableBand[2];
    repaint();
    break;
  case '3':
    bEnableBand[3] = !bEnableBand[3];
    repaint();
    break;
  case '4':
    bEnableBand[4] = !bEnableBand[4];
    repaint();
    break;
  case '5':
    bEnableBand[5] = !bEnableBand[5];
    repaint();
    break;
  case '6':
    bEnableBand[6] = !bEnableBand[6];
    repaint();
    break;
  case '7':
    bEnableBand[7] = !bEnableBand[7];
    repaint();
    break;
  case 'i':
    bSelectAll = !bSelectAll;
    break;
  case 'q':
    if(bAdjust) {
      bAdjust = false;
      bForceAdjust = true;
    } else if(bForceAdjust) {
      bForceAdjust = false;
      bAdjust = false;
    } else {
      bAdjust = true;
      bForceAdjust = false;
    }
    repaint();
    break;
  case 'z':
    markTrackOnset();
    break;
  case 'x':
    markTrackOffset();
    break;
  case ']':
    if(timeAxisMode == TimeAxisMSM)  {
      timeAxisMode = TimeAxisGrains;
    } else if(timeAxisMode == TimeAxisGrains) {
      timeAxisMode = TimeAxisSamples;
    } else {
      timeAxisMode = TimeAxisMSM;
    }
    repaint();
    break;
  case '[':
    if(freqAxisMode == FrequencyAxisHz)  {
      freqAxisMode = FrequencyAxis2PI;
    } else {
      freqAxisMode = FrequencyAxisHz;
    }
    repaint();
    break;
  }


  if(key.isKeyCode(ASCII::a)) {
    if(key.getModifiers().isCommandDown()) {
      selectAll();
    }
  } else if(key.isKeyCode(ASCII::s)) {
    if(key.getModifiers().isCommandDown()) {
      saveSelection();
    } 
  } else if(key.isKeyCode(KeyPress::returnKey)) {
    outputSelectedTrack();
  } else if(key.isKeyCode(KeyPress::spaceKey)) {
    fprintf(stderr,"play %d\n",key.getModifiers().isCommandDown());
    if(sbsmsAudioSource->isPlaying()) {
      sbsmsAudioSource->stop();
      view->currPos = -1;
      view->repaint();
    } else {
      startTimer(50);
      sbsmsAudioSource->setLeftPos(view->leftPos);
      sbsmsAudioSource->setRightPos(view->rightPos==-1?totSamples:view->rightPos);
      sbsmsAudioSource->play(key.getModifiers().isCommandDown());
    }
  } else if(key.isKeyCode(KeyPress::leftKey)) {
    pthread_mutex_lock(&glMutex);
    SampleCountType range = view->endSample - view->startSample;
    view->startSample = max((SampleCountType)0,view->startSample-range/10);
    view->endSample = view->startSample + range;
    pthread_mutex_unlock(&glMutex);
    setHScrollRange();
    repaint();
  } else if(key.isKeyCode(KeyPress::rightKey)) {
    pthread_mutex_lock(&glMutex);
    SampleCountType range = view->endSample - view->startSample;
    view->endSample = min((SampleCountType)totSamples,view->endSample+range/10);
    view->startSample = view->endSample - range;
    pthread_mutex_unlock(&glMutex);
    setHScrollRange();
    repaint();
  } else if(key.isKeyCode(KeyPress::downKey)) {
    pthread_mutex_lock(&glMutex);
    double range = view->endF - view->startF;
    view->startF = max(0.0,view->startF-range/10.0);
    view->endF = view->startF + range;
    pthread_mutex_unlock(&glMutex);
    setVScrollRange();
    repaint();
  } else if(key.isKeyCode(KeyPress::upKey)) {
    pthread_mutex_lock(&glMutex);
    double range = view->endF - view->startF;
    view->endF = min(maxF,view->endF+range/10.0);
    view->startF = view->endF - range;
    pthread_mutex_unlock(&glMutex);
    setVScrollRange();
    repaint();
  } else if(key.isKeyCode(KeyPress::homeKey)) {
    view->leftPos = 0;
    view->rightPos = -1;
    repaint();
  }
}

void TimeFrequencyView :: mouseWheelMove (const MouseEvent& e, const MouseWheelDetails& d)
{
  float scaleX = 0.1f;
  float scaleY = 0.1f;
  float deltaX = d.isReversed?d.deltaX:-d.deltaX;
  if(deltaX > 0) {
    pthread_mutex_lock(&glMutex);
    SampleCountType range = view->endSample - view->startSample;
    view->startSample = max((SampleCountType)0,view->startSample - lrintf(scaleX * deltaX * range));
    view->endSample = view->startSample + range;
    pthread_mutex_unlock(&glMutex);
    setHScrollRange();
    repaint();
  } else if(deltaX < 0) {
    pthread_mutex_lock(&glMutex);
    SampleCountType range = view->endSample - view->startSample;
    view->endSample = min((SampleCountType)totSamples,view->endSample - lrintf(scaleX * deltaX * range));
    view->startSample = view->endSample - range;
    pthread_mutex_unlock(&glMutex);
    setHScrollRange();
    repaint();
  }

  if(d.deltaY < 0) {
    pthread_mutex_lock(&glMutex);
    float range = view->endF - view->startF;
    view->startF = max(0.0,view->startF + scaleY * d.deltaY * range);
    view->endF = view->startF + range;
    pthread_mutex_unlock(&glMutex);
    setVScrollRange();
    repaint();
  } else if(d.deltaY > 0) {
    pthread_mutex_lock(&glMutex);
    float range = view->endF - view->startF;
    view->endF = min(maxF,view->endF + scaleY * d.deltaY * range);
    view->startF = view->endF - range;
    pthread_mutex_unlock(&glMutex);
    setVScrollRange();
    repaint();
  }
}

void TimeFrequencyView :: zoomT(bool bIn)
{
  pthread_mutex_lock(&glMutex);
  SampleCountType range = view->endSample - view->startSample;
  if(bIn) {
    range /= 4;
    SampleCountType newStartSample = view->leftPos - range;
    view->startSample = max((SampleCountType)0,newStartSample);
    SampleCountType newEndSample = view->leftPos + range;
    view->endSample = min(totSamples,newEndSample);
    } else {
    range /= 2;
    SampleCountType newStartSample = view->startSample - range;
    view->startSample = max((SampleCountType)0,newStartSample);
    SampleCountType newEndSample = view->endSample + range;
    view->endSample = min(totSamples,newEndSample);
  }
  pthread_mutex_unlock(&glMutex);
}

void TimeFrequencyView :: zoomF(bool bIn)
{
  pthread_mutex_lock(&glMutex);
  double rangeY = view->endF - view->startF;
  if(bIn) {
    rangeY /= 2;
    double newEndF = view->endF - rangeY;
    view->endF = min(maxF,newEndF);   
  } else {
    rangeY /= 2;
    double newStartF = view->startF - rangeY;
    view->startF = max(0.0,newStartF);
    double newEndF = view->endF + rangeY;
    view->endF = min(maxF,newEndF);   
  }
  pthread_mutex_unlock(&glMutex);
}

void TimeFrequencyView :: setView(int pos0, int pos1, float f0, float f1)
{
  pthread_mutex_lock(&glMutex);
  view->startSample = pos0 * view->frameSize;
  view->endSample = pos1 * view->frameSize;
  view->startF = f0 * TWOPI / 44100.0f;
  view->endF = f1 * TWOPI / 44100.0f;
  view->rescale();
  pthread_mutex_unlock(&glMutex);
  repaint();
}

void TimeFrequencyView :: resized()
{
  int w = getWidth()-scrollSize-yAxisWidth-spectrumWidth;
  int h = getHeight()-scrollSize-xAxisHeight-trackViewHeight-controlsHeight-waveHeight;

  hScroll->setBounds(yAxisWidth+spectrumWidth,waveHeight,w,scrollSize);
  vScroll->setBounds(getWidth()-scrollSize,waveHeight+scrollSize-1,scrollSize,h);
  view->setBounds(yAxisWidth+spectrumWidth,scrollSize+waveHeight,w,h);
  spectrum1->setBounds(0,scrollSize+waveHeight,spectrumWidth,h);
  spectrum2->setBounds(0,scrollSize+waveHeight,spectrumWidth,h);
  spectrumTrial->setBounds(0,scrollSize+waveHeight,spectrumWidth,h);
  trackView->setBounds(yAxisWidth+spectrumWidth,getHeight()-controlsHeight-trackViewHeight,w,trackViewHeight);
  waveDisplay->setBounds(yAxisWidth+spectrumWidth,0,w,waveHeight);
  status->setBounds(2,getHeight()-controlsHeight + 2,spectrumWidth,controlsHeight - 4);
  rateCtrl->setBounds(spectrumWidth + 12,getHeight()-controlsHeight + 2,80,controlsHeight - 4);
}


void TimeFrequencyView :: scrollBarMoved (ScrollBar *scrollBar, double newRangeStart)
{
  pthread_mutex_lock(&glMutex);
  if(scrollBar == hScroll) {
    view->startSample = (SampleCountType)newRangeStart;
    view->endSample = (SampleCountType)scrollBar->getCurrentRange().getEnd();
  } else {
    view->endF = -newRangeStart;
    view->startF = -scrollBar->getCurrentRange().getEnd();
  }
  view->rescale();
  pthread_mutex_unlock(&glMutex);
}

void TimeFrequencyView :: setHScrollRange()
{
  pthread_mutex_lock(&glMutex);
  view->rescale();
  pthread_mutex_unlock(&glMutex);
  hScroll->setCurrentRange(view->startSample,view->endSample-view->startSample);
}

void TimeFrequencyView :: setVScrollRange()
{
  pthread_mutex_lock(&glMutex);
  view->rescale();
  pthread_mutex_unlock(&glMutex);
  vScroll->setCurrentRange(-view->endF,view->endF-view->startF);
}

void TimeFrequencyView :: sliderValueChanged(Slider * slider)
{
  sbsmsAudioSource->setRate(slider->getValue());
}

void TimeFrequencyView :: timerCallback()
{
  if(!sbsmsAudioSource->isPlaying()) {
    stopTimer();
    view->currPos = -1;
  }
  view->currPos = sbsmsAudioSource->getCurrPos();
  view->repaint();
}

