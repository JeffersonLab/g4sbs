#ifndef G4SBSECaloutput_h 
#define G4SBSECaloutput_h 1

//#include <set>
#include <vector>
//#include "G4ThreeVector.h"

using namespace std;

class G4SBSECaloutput {
public:
  G4SBSECaloutput();
  ~G4SBSECaloutput();

  void Clear();
  
  //output of the ECal detector simulation consists of hits in individual PMTs. "Hits" in GEANT4 are really individual tracking steps. 
  //In our end of event routine, we will collect all tracking steps (GEANT4 "hits") occuring within a certain time window on the same PMT into "logical hits"
  double timewindow, threshold; //Number of photo-electrons within timewindow must exceed threshold

  //"Hit" variables:
  int nhits_ECal; //number of "logical hits" in the ECal detector.
  vector<int> PMTnumber;
  vector<int> row;
  vector<int> col;
  vector<int> plane;
  vector<double> xcell,ycell,zcell, xgcell, ygcell, zgcell;
  //vector<int> NumPhotons;
  vector<int> NumPhotoelectrons;
  vector<double> Time_avg; //average global time of a hit
  vector<double> Time_rms; //std. deviation of photon arrival times in this hit.
  vector<double> Time_min; //arrival time of earliest photon in the hit
  vector<double> Time_max; //arrival time of latest photon in the hit
  
  //G4SBSParticleOutput ParticleHistory;
  
  
};

#endif
