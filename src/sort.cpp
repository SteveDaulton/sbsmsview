#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <utility>

using namespace std;

class AnalysisTrackPoint {
public:
  float f;
  int order;
  AnalysisTrackPoint *cont;
};


bool TrackPointSortFunc(AnalysisTrackPoint *tp0, AnalysisTrackPoint *tp1)
{
  return tp0->f < tp1->f;
}

int main(int argc, char **argv)
{
  int N = atoi(argv[1]);
  
  AnalysisTrackPoint *tp0s = new AnalysisTrackPoint[N];
  AnalysisTrackPoint *tp1s = new AnalysisTrackPoint[N];
  vector<AnalysisTrackPoint*> left;
  vector<AnalysisTrackPoint*> right;

  for(int i=0; i<N; i++) {
    tp0s[i].f = (float)i/(float)N;
    tp1s[i].f = (float)rand()/(float)RAND_MAX;
    tp0s[i].cont = &tp1s[i];
    tp1s[i].cont = &tp0s[i];
    left.push_back(&tp0s[i]);
    right.push_back(&tp1s[i]);
  }

  vector<pair<AnalysisTrackPoint*,AnalysisTrackPoint*> > crossings;
  sort(right.begin(),right.end(),TrackPointSortFunc);
  for(int i=0; i<left.size(); i++) {
    left[i]->order = i;
  }
  for(int i=0; i<right.size(); i++) {
    right[i]->order = i;
  }
  for(int i=0; i<left.size(); i++) {
    float f0 = left[i]->f;
    float f1 = left[i]->cont->f;
    printf("%g %g\n",f0,f1);
  }
  exit(0);
  for(int i=0; i<left.size(); i++) {
    int ro = left[i]->cont->order;
    printf("%d %d\n",i,ro);
    for(int j=0; j<ro; j++) {
      if(right[j]->cont->order > i) {
        pair<AnalysisTrackPoint*,AnalysisTrackPoint*> crossing(left[i],right[j]->cont);
        crossings.push_back(crossing);
        printf("%d %d %g %g %g %g\n",i,j,left[i]->f,left[i]->cont->f,right[j]->cont->f,right[j]->f);
      }
    }
  }
  printf("\n");
}
