#include "G4SBSECaloutput.hh"
#include "G4SystemOfUnits.hh"

G4SBSECaloutput::G4SBSECaloutput(){
  //Set sensible default values for threshold and time window:

  timewindow = 250.0*ns; //Also declared in G4SBSEventAction.cc -> G4SBSECaloutput section
  threshold  = 0.5; //single photo-electron threshold!
  

  Clear();
}

G4SBSECaloutput::~G4SBSECaloutput(){
  ;
}

void G4SBSECaloutput::Clear(){
  nhits_ECal = 0;
  //ntracks_ECal = 0;
  
  PMTnumber.clear();
  row.clear();
  col.clear();
  plane.clear();
  xcell.clear();
  ycell.clear();
  zcell.clear();
  xgcell.clear();
  ygcell.clear();
  zgcell.clear();
  //NumPhotons.clear();
  NumPhotoelectrons.clear();
  Time_avg.clear();
  Time_rms.clear();
  Time_min.clear();
  Time_max.clear();

  otridx.clear();
  ptridx.clear();
  sdtridx.clear();

  // Clear individual particle list
  npart_ECAL = 0;
  part_PMT.clear();
  E.clear();
  t.clear();
  trid.clear();
  detected.clear();
  //ParticleHistory.Clear();
  /* 
  mTrackNo.clear();

  xhit.clear();
  yhit.clear();
  zhit.clear();
  
  pxhit.clear();
  pyhit.clear();
  pzhit.clear();

  pvx.clear();
  pvy.clear();
  pvz.clear();

  ppx.clear();
  ppy.clear();
  ppz.clear();
  */
  //  thetaC.clear();
  // volume_flag.clear();
  
  // mPID.clear();
  // mTID.clear();
  // mMID.clear();
  // mvx.clear();
  // mvy.clear();
  // mvz.clear();
  // mpx.clear();
  // mpy.clear();
  // mpz.clear();

}
