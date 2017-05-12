// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class SteppingMessenger
// Control of max time for tracking
// 20/05/13 JRMA adapted from SBS equivalent, under construction
//
#ifndef SteppingMessenger_h
#define SteppingMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class SteppingAction;
class G4UIdirectory;
class G4UIcmdWithADoubleAndUnit;

class SteppingMessenger: public G4UImessenger
{
public:
  SteppingMessenger(SteppingAction*);
  ~SteppingMessenger();
  void SetNewValue(G4UIcommand*, G4String);  
private:
  SteppingAction* SBSAction;
  G4UIdirectory* fStepDir;
  G4UIcmdWithADoubleAndUnit* SetMaxTimeCmd;

};

#endif

