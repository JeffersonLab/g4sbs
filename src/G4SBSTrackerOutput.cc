#include "G4SBSTrackerOutput.hh"

G4SBSTrackerOutput::G4SBSTrackerOutput(){
  Clear();
}

G4SBSTrackerOutput::~G4SBSTrackerOutput(){
  ;
}

void G4SBSTrackerOutput::Clear(){
  ntracks = 0;
  //TrackerID.clear();
  TrackTID.clear();
  TrackPID.clear();
  TrackMID.clear();
  
  NumHits.clear();
  NumPlanes.clear();
  NDF.clear();
  Chi2fit.clear();
  Chi2true.clear();

  TrackX.clear();
  TrackY.clear();
  TrackXp.clear();
  TrackYp.clear();
  TrackT.clear();
  TrackP.clear();

  TrackSx.clear();
  TrackSy.clear();
  TrackSz.clear();
  
  TrackXfit.clear();
  TrackYfit.clear();
  TrackXpfit.clear();
  TrackYpfit.clear();

  otridx.clear();
  ptridx.clear();
  sdtridx.clear();
}
