#ifndef G4SBSGEMoutput_h
#define G4SBSGEMoutput_h 1

#include <vector>
#include "G4SBSParticleOutput.hh"
using namespace std;

class G4SBSGEMoutput {
public:
  G4SBSGEMoutput();
  ~G4SBSGEMoutput();

  void Clear();

  double timewindow, threshold;

  int nhits_GEM;
  vector<int> plane, strip;
  vector<double> x,y,z,t,trms,tmin,tmax;
  vector<double> polx,poly,polz;
  vector<double> dx,dy;
  vector<double> tx,ty; 
  vector<double> xin,yin,zin;
  vector<double> xout,yout,zout;
  vector<double> txp,typ;
  vector<double> xg,yg,zg;
  vector<int> trid,mid,pid;
  vector<double> vx,vy,vz;
  vector<double> p,edep,beta;
  //Add bookkeeping indices for "original", "primary", and "SD boundary" tracks:
  vector<int> otridx,ptridx,sdtridx;
  
  G4SBSParticleOutput ParticleHistory;

};

#endif
