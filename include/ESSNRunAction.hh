// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class ESSNRunAction
// MC run control
// 20/05/13 JRMA adapted from SBS equivalent, under construction
//
#ifndef ESSNRunAction_h
#define ESSNRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "ESSNEventAction.hh"

class G4Run;

class ESSNRunAction : public G4UserRunAction
{
  public:
    ESSNRunAction();
   ~ESSNRunAction();
  public:
    void BeginOfRunAction(const G4Run*);
    void   EndOfRunAction(const G4Run*);
  private:
  ESSNEventAction *fEventAction;
};

#endif

