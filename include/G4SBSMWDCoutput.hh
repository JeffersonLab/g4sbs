#ifndef G4SBSMWDCoutput_h
#define G4SBSMWDCoutput_h 1

#include <vector>
#include "G4SBSParticleOutput.hh"
using namespace std;

class G4SBSMWDCoutput {
public:
  G4SBSMWDCoutput();
  ~G4SBSMWDCoutput();

  void Clear();

  double timewindow, threshold;

  int nhits_MWDC;
  vector<int> plane, strip;
  vector<double> x,y,z,t,trms,tmin,tmax;
  vector<double> polx,poly,polz;
  
  vector<double> dx,dy;
  vector<double> tx,ty,txp,typ;
  vector<int> trid,mid,pid;
  vector<double> vx,vy,vz;
  vector<double> p,edep,beta;
  
  G4SBSParticleOutput ParticleHistory;

};

#endif
