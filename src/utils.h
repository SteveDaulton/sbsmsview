// -*- mode: c++ -*-
#ifndef UTILS_H
#define UTILS_H

#include "real.h"
#include "config.h"
#include "sbsms.h"
#include <algorithm>
#include <queue>
using namespace std;

namespace _sbsms_ {

#define ONEOVERTWOPI 0.15915494309189533576888376337251f
#define PI 3.1415926535897932384626433832795f
#define TWOPI 6.28318530717958647692528676655900576f

inline void c2even(audio *eo, audio *even, int N)
{
  int Nover2 = N/2;
  even[0][0] = eo[0][0];
  even[0][1] = 0.0f;
  for(int k=1;k<=Nover2;k++) {
    int Nk = N-k;
    even[k][0] = 0.5f*(eo[k][0] + eo[Nk][0]);
    even[k][1] = 0.5f*(eo[k][1] - eo[Nk][1]);
  }
}

inline void c2odd(audio *eo, audio *odd, int N)
{
  int Nover2 = N/2;
  odd[0][0] = eo[0][1];
  odd[0][1] = 0.0f;
  for(int k=1;k<=Nover2;k++) {
    int Nk = N-k;
    odd[k][0] = 0.5f*(eo[k][1] + eo[Nk][1]);
    odd[k][1] = 0.5f*(eo[Nk][0] - eo[k][0]);
  }
}

inline float canonPI(float ph) 
{
  ph -= TWOPI * lrintf(ph * ONEOVERTWOPI);
  if(ph < -PI) ph += TWOPI;
  else if(ph >= PI) ph -= TWOPI;
  return ph;
}

inline float canon2PI(float ph) 
{
  ph -= TWOPI * lrintf(ph * ONEOVERTWOPI);
  if(ph < 0.0f) ph += TWOPI;
  if(ph >= TWOPI) ph -= TWOPI;
  return ph;
}

inline float square(float x)
{ 
  return x*x;
}

inline float norm2(t_fft x)
{
  return square(x[0]) + square(x[1]);
}

template<typename T>
class Diff {
public:
  Diff() : d(0), t1(0) {}
  T d;
  T t1;
  void push(T t) {
    d = t - t1;
    t1 = t;
  }
  
  T const get() const { return d; }
};

template<typename T>
class Average {
public:
  Average() : N(1), a(0), f(0) {}

  void setSize(int N) { this->N = N; }
  void push(T t) {
    q.push(t);
    a += t;
    if(q.size() <= N) {
      f = (T)1/(T)q.size();
    } else {
      a-= q.front();
      q.pop();
    }
  }
  
  T const get() const {
    return f * a;
  }

  queue<T> q;
  int N;
  T a;
  T f;
};

}

#endif
