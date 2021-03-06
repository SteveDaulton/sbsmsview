#include <iostream>
#include "import.h"
#include "pcm.h"
#include "mp3.h"
#include "ffmpeg.h"
#include <map>
#include <string>

using namespace std;

string lower(string strToConvert)
{
  for(unsigned int i=0;i<strToConvert.length();i++) {
    strToConvert[i] = tolower(strToConvert[i]);
  }
  return strToConvert;
}
 
AudioDecoder *import(const char *filename)
{
  string fname(filename);
  size_t i = fname.find(".");
  AudioDecoder *decoder = NULL;

  if(i) {
    string ext = fname.substr(i+1);
    string extl = lower(ext);
    if(!extl.compare("wav") || !extl.compare("aif") || !extl.compare("aiff")) {
      decoder = new PcmReader(filename);
    } else if (!extl.compare("mp3")) {
      decoder = new MP3Reader(filename);
    } else if (!extl.compare("aac") || !extl.compare("m4a")) {
      decoder = new FFMpegDecoder(filename);
    } else {
      perror("Error importing file");
      return NULL;
    }
  }
  if(decoder->isError()) {
    abort();
  }
  return decoder;
}
