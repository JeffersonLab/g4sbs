#include "G4SBSGEMoutput.hh"
#include "G4SystemOfUnits.hh"

G4SBSGEMoutput::G4SBSGEMoutput(){
  nhits_GEM = 0;
  timewindow = 1000.0*ns;
  threshold = 0.0*eV;
  Clear();
}

G4SBSGEMoutput::~G4SBSGEMoutput(){
  ;
}

void G4SBSGEMoutput::Clear(){
  nhits_GEM = 0;
  timewindow = 1000.0*ns;
  threshold = 0.0*eV;
  plane.clear();
  strip.clear();
  x.clear();
  y.clear();
  z.clear();
  polx.clear();
  poly.clear();
  polz.clear();
  t.clear();
  trms.clear();
  tmin.clear();
  tmax.clear();
  dx.clear();
  dy.clear();
  tx.clear();
  ty.clear();
  xin.clear();
  yin.clear();
  zin.clear();
  xout.clear();
  yout.clear();
  zout.clear();
  txp.clear();
  typ.clear();
  //EFuchey 2017-01-30: add hit position in global frame
  xg.clear();
  yg.clear();
  zg.clear();
  trid.clear();
  mid.clear();
  pid.clear();
  vx.clear();
  vy.clear();
  vz.clear();
  p.clear();
  edep.clear();
  beta.clear();
  otridx.clear();
  ptridx.clear();
  sdtridx.clear();

  ParticleHistory.Clear();
}
