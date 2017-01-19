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
  //EFuchey 2017-01-13: add track position at the entrance and exit of the ionizable gas layer
  // x_in.clear();
  // y_in.clear();
  // z_in.clear();
  // x_out.clear();
  // y_out.clear();
  // z_out.clear();
  dx.clear();
  dy.clear();
  tx.clear();
  ty.clear();
  txp.clear();
  typ.clear();
  trid.clear();
  mid.clear();
  pid.clear();
  vx.clear();
  vy.clear();
  vz.clear();
  p.clear();
  edep.clear();
  beta.clear();

  ParticleHistory.Clear();
}
