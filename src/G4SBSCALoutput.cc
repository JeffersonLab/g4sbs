#include "G4SBSCALoutput.hh"
#include "G4SystemOfUnits.hh"

G4SBSCALoutput::G4SBSCALoutput(){
  timewindow = 1000.0*ns;
  threshold = 0.0*eV;
  keeppart = true;
  Clear();
}

G4SBSCALoutput::~G4SBSCALoutput(){
  ;
}

void G4SBSCALoutput::Clear(){
  nhits_CAL = 0;
  //timewindow = 1000.0*ns;
  //threshold = 0.0*eV;

  Esum = 0.0;
  
  row.clear();
  col.clear();
  plane.clear();
  wire.clear();
  cell.clear();
  xcell.clear();
  ycell.clear();
  zcell.clear();
  xcellg.clear();
  ycellg.clear();
  zcellg.clear();
  xhit.clear();
  yhit.clear();
  zhit.clear();
  sumedep.clear();
  tavg.clear();
  trms.clear();
  tmin.clear();
  tmax.clear();

  npart_CAL = 0;
  px.clear();
  py.clear();
  pz.clear();
  x.clear();
  y.clear();
  z.clear();
  t.clear();
  E.clear();
  dt.clear();
  L.clear();
  // trms.clear();
  vx.clear();
  vy.clear();
  vz.clear();
  trid.clear();
  mid.clear();
  pid.clear();
  p.clear();
  edep.clear();

  ParticleHistory.Clear();
}
