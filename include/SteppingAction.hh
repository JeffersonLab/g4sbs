// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class SteppingAction
// Max time to continue tracking
// 20/05/13 JRMA adapted from SBS equivalent, under construction

#ifndef SteppingAction_h
#define SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class DetectorConstruction;
class EventAction;
class SteppingMessenger;


class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(DetectorConstruction*, EventAction*);
  ~SteppingAction();
  void UserSteppingAction(const G4Step*);
  G4double GetMaxTime(){ return fMaxTime; }
  void SetMaxTime(G4double t){ fMaxTime = t; }
private:
  DetectorConstruction* detector;
  EventAction*          eventaction;
  G4double fMaxTime;
  SteppingMessenger* fSteppingMessenger;
};


#endif
