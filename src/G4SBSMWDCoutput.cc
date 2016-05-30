#include "G4SBSMWDCoutput.hh"
#include "G4SystemOfUnits.hh"

G4SBSMWDCoutput::G4SBSMWDCoutput(){
  nhits_MWDC = 0;
  timewindow = 1000.0*ns;
  threshold = 0.0*eV;
  Clear();
}

G4SBSMWDCoutput::~G4SBSMWDCoutput(){
  ;
}

void G4SBSMWDCoutput::Clear(){
  nhits_MWDC = 0;
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
