#ifndef G4SBSRICHoutput_h 
#define G4SBSRICHoutput_h 1

//#include <set>
#include <vector>
#include "G4SBSParticleOutput.hh"
//#include "G4ThreeVector.h"

using namespace std;

class G4SBSRICHoutput {
public:
  G4SBSRICHoutput();
  ~G4SBSRICHoutput();

  void Clear();
  
  //output of the RICH detector simulation consists of hits in individual PMTs. "Hits" in GEANT4 are really individual tracking steps. 
  //In our end of event routine, we will collect all tracking steps (GEANT4 "hits") occuring within a certain time window on the same PMT into "logical hits"
  double timewindow, threshold; //Number of photo-electrons within timewindow must exceed threshold

  //"Hit" variables:
  int nhits_RICH; //number of "logical hits" in the RICH detector.
  vector<int> PMTnumber;
  vector<int> row;
  vector<int> col;
  //vector<double> xpmt, ypmt, zpmt;
  //vector<int> NumPhotons;
  vector<int> NumPhotoelectrons;
  vector<double> Time_avg; //average global time of a hit
  vector<double> Time_rms; //std. deviation of photon arrival times in this hit.
  vector<double> Time_min; //Earliest photon arrival time in this hit
  vector<double> Time_max; //Latest photon arrival time in this hit
  vector<int> mTrackNo; //Track of the mother particle that produced this hit. Note: THIS IS NOT THE GEANT4 TRACK ID, BUT RATHER THE INDEX IN THE ARRAY THAT 
                          //IS DEFINED BELOW!
  //Average position of photons in this hit (local coordinates):
  vector<double> xhit;
  vector<double> yhit;
  vector<double> zhit;
  //Average momentum direction of photons in this hit (local coordinates):
  vector<double> pxhit;
  vector<double> pyhit;
  vector<double> pzhit; 

  //Average emission vertex of photons in this hit (global coordinates)
  vector<double> pvx;
  vector<double> pvy; 
  vector<double> pvz;
  //Average momentum ***direction*** of photons in this hit when emitted (global coordinates) (ignore photon energy here)
  vector<double> ppx;
  vector<double> ppy; 
  vector<double> ppz;
  //Average emission angle of photons in this hit (wrt mother track):
  //vector<double> thetaC; //how do we figure this out?
  vector<int>    volume_flag; //flag corresponding to production volume: 1=aerogel, 2=gas, 3=lucite, 0=other

  //Add PMT coordinates, both "local" to the the detector mother volume and "global" to the Hall
  vector<double> xpmt,ypmt,zpmt;
  vector<double> xgpmt,ygpmt,zgpmt;

  // //"Track" variables:
  // int ntracks_RICH; //number of tracks producing RICH hits.
  // vector<int> mPID; //Particle ID of tracks producing RICH hits. (according to PDG encoding scheme).
  // vector<int> mTID; //Track ID number of tracks producing RICH hits
  // vector<int> mMID; //Mother Track ID number of tracks producing RICH hits: for primary particles, this is zero!
  // //vertex position of production of track producing RICH hits:
  // vector<double> mvx; 
  // vector<double> mvy;
  // vector<double> mvz;
  // //Initial momentum of particle at vertex:
  // vector<double> mpx;
  // vector<double> mpy;
  // vector<double> mpz;
  // //Difficult (and costly) to keep track of mother particle momentum at each Ckov photon emission vertex, but perhaps we can 
  // //calculate some kind of average. Let's defer this question for now.
  
  G4SBSParticleOutput ParticleHistory;

};

#endif
