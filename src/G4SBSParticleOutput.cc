#include "G4SBSParticleOutput.hh"

G4SBSParticleOutput::G4SBSParticleOutput(){
  npart=0;
  Clear();
}

G4SBSParticleOutput::~G4SBSParticleOutput(){;}

void G4SBSParticleOutput::Clear(){
  npart=0;
  PID.clear();
  MID.clear();
  TID.clear();
  nbounce.clear();
  vx.clear();
  vy.clear();
  vz.clear();
  px.clear();
  py.clear();
  pz.clear();
  
  hitindex.clear();
  
}
