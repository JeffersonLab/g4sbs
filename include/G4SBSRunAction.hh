
#ifndef G4SBSRunAction_h
#define G4SBSRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

class G4Timer;
class G4Run;
class G4SBSIO;
class G4SBSTrackingAction;
class G4SBSSteppingAction;

class G4SBSRunAction : public G4UserRunAction
{
public:
  G4SBSRunAction();
  ~G4SBSRunAction();
  
public:
  void BeginOfRunAction(const G4Run* aRun);
  void EndOfRunAction(const G4Run* aRun);
  
  void SetIO( G4SBSIO *io ){ fIO = io; }
  void SetTrackingAction( G4SBSTrackingAction *trkact ){ ftrkact = trkact; }
  void SetSteppingAction( G4SBSSteppingAction *stepact ){ fstepact = stepact; }
  
  void SetNtries( int n ){ Ntries = n; }
  G4int GetNtries(){ return Ntries; }

private:
  G4int Ntries; //Keep track of total number of tries to throw an event during the run (efficiency of MC generation needed for correct normalization).
  
  G4Timer* timer;
  
  G4SBSIO *fIO;
  G4SBSTrackingAction *ftrkact;
  G4SBSSteppingAction *fstepact;
};

#endif

