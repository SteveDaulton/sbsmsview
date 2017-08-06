// -*- mode: c++ -*-
#ifndef DBTABLE_H
#define DBTABLE_H

#include "real.h"

namespace _sbsms_ {

enum { dBTableSize = 4096, dBTableScale = dBTableSize - 1 };

extern float dBTable[dBTableSize];
extern float dB2Table[dBTableSize];

inline float dB2Approx(float x, float y) 
{
  if(x < y) {
    return dB2Table[lrintf(x/y*dBTableScale)];
  } else if(x == 0.0f) {
    return 0;
  } else {
    return dB2Table[lrintf(y/x*dBTableScale)];
  }
}

inline float dBApprox(float x, float y) 
{
  if(x < y) {
    return dBTable[lrintf(x/y*dBTableScale)];
  } else if(x == 0.0f) {
    return 0;
  } else {
    return dBTable[lrintf(y/x*dBTableScale)];
  }
}

}

#endif
