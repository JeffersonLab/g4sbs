// SBS a Geant-4 Based Model of Hall-A Experiments with 11 GeV
// J.R.M Annand, University of Glasgow
// Class EventActionMessenger
// Online control of Event Generator via keyboard
// 20/06/09 JRMA adapted from A2 equivalent, under construction

#ifndef EventActionMessenger_h
#define EventActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class EventAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithoutParameter;

//----------------------------------------------------------------------------
class EventActionMessenger: public G4UImessenger
{
public:
  EventActionMessenger(EventAction*);
  ~EventActionMessenger();
  void SetNewValue(G4UIcommand*, G4String);  
private:
  EventAction*     feventAction;
  G4UIdirectory*        feventDir;   
  G4UIcmdWithAString*   fDrawCmd;
  G4UIcmdWithAString*   fOutFileCmd;
  G4UIcmdWithAString*   fHitDrawCmd;
  G4UIcmdWithAnInteger* fPrintCmd;     
};

#endif
