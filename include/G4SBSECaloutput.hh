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
  int ntimebins;

  vector< vector<double> > NPE_vs_time;
  double gatewidth;

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

  //Add bookkeeping indices for "original", "primary", and "SD boundary" tracks:
  vector<int> otridx,ptridx,sdtridx;
  //G4SBSParticleOutput ParticleHistory;
  //"Part" keeps track of all unique particles detected in the ECal
  int npart_ECAL;   // Number of optical photons detected in this hit
  //vector<int> ihit; // hit index associated with this particle 
  vector<int> part_PMT; // PMT number for this particle track
  vector<int> trid; // Track ID of this optical photon
  vector<double> E; // Energy of particle
  vector<double> t; // Global time of particle
  vector<bool> detected; // Was this photon detected?
  
};

#endif
