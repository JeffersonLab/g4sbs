// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class SteppingAction
// Max time to continue tracking
// 20/05/13 JRMA adapted from SBS equivalent, under construction

#include "SteppingAction.hh"
#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "G4Track.hh"
#include "G4Gamma.hh"
#include "G4Proton.hh"
#include "G4SteppingManager.hh"
#include "SteppingMessenger.hh"

SteppingAction::SteppingAction(DetectorConstruction* det,
                                         EventAction* evt)
:detector(det), eventaction(evt)                                         
{ 
  fMaxTime = 1000.0;         // track up to 1000 ns
  fSteppingMessenger = new SteppingMessenger(this);
}

SteppingAction::~SteppingAction()
{ }

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  //  return;
  G4Track* track = aStep->GetTrack();
  // abort track if over max time
  if(track->GetGlobalTime() > fMaxTime*ns)
    track->SetTrackStatus(fStopAndKill);
  // bug in phot process, can't get rid of gamma with energy 1.2E-5MeV
  // goes into infinite loop!
  /*
  if( ( track->GetDefinition()->GetParticleName() ==
	G4Gamma::Gamma()->GetParticleName()) && 
      ( track->GetKineticEnergy()/MeV < 1E-4 ) )
    track->SetTrackStatus(fStopAndKill);
  if( ( track->GetDefinition()->GetParticleName() ==
	G4Proton::Proton()->GetParticleName()) && 
      ( track->GetKineticEnergy()/MeV < 1E-1 ) )
    track->SetTrackStatus(fStopAndKill);
  */
}





