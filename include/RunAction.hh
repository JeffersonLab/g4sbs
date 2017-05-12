// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class RunAction
// MC run control
// 20/05/13 JRMA adapted from SBS equivalent, under construction
//
#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "EventAction.hh"

class G4Run;

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
   ~RunAction();
  public:
    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);
  private:
  EventAction *fEventAction;
};

#endif

