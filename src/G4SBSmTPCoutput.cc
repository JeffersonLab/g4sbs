#include "G4SBSmTPCoutput.hh"
#include "G4SystemOfUnits.hh"

G4SBSmTPCoutput::G4SBSmTPCoutput(){
  timewindow = 1000.0*ns;
  threshold = 0.0*eV;
  // keeppart = true;
  Clear();
}

G4SBSmTPCoutput::~G4SBSmTPCoutput(){
  ;
}

void G4SBSmTPCoutput::Clear(){
  nhits_mTPC = 0;
  //timewindow = 1000.0*ns;
  //threshold = 0.0*eV;

  Esum = 0.0;
  
  // row.clear();
  // col.clear();
  // plane.clear();
  // wire.clear();
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
  mid.clear();
  pid.clear();
  trid.clear();
  tmax.clear();
  // Ehit = 0.0;
  // pxhit = 0.0;
  // pyhit = 0.0;
  // pzhit = 0.0;
  hitL.clear();
 
  npart_mTPC = 0;
  px.clear();
  py.clear();
  pz.clear();
  px_v.clear();
  py_v.clear();
  pz_v.clear();
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
  trid_.clear();
  mid_.clear();
  pid_.clear();
  p.clear();
  edep.clear();
  ztravel.clear();
  nstrips.clear();

  ParticleHistory.Clear();
}
