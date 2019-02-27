#include "G4SBSSteppingAction.hh"
//#include "G4SBSSteppingActionMessenger.hh"

#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SteppingManager.hh"

#include "G4StepPoint.hh"
#include "G4Material.hh"
#include "G4VProcess.hh"


G4SBSSteppingAction::G4SBSSteppingAction()
:drawFlag(false)
{
 // new G4SBSSteppingActionMessenger(this);
}

// void G4SBSSteppingAction::UserSteppingAction(const G4Step*)

//Montgomery Nov 2018, added stepping action functionality for option to kill any tracks enterring GEM and readout pcb of mtpc, if want to preserve kryptonite functionality from GEMC - I do not think this is the best way to do this. Temporary - must check!
void G4SBSSteppingAction::UserSteppingAction(const G4Step* theStep)
{


  // // get volume of the current step
  // G4Track* theTrack = theStep->GetTrack();

  // G4StepPoint *prestep = theStep->GetPreStepPoint();
  // G4Material *mat = prestep->GetMaterial();
  // G4StepPoint *poststep = theStep->GetPostStepPoint();
  // G4Material *mat2 = prestep->GetMaterial();
  // if(mat->GetName() != mat2->GetName())  theTrack->SetTrackStatus(fStopAndKill);



  // G4VPhysicalVolume* volume = theStep->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  // // G4VPhysicalVolume* volume = theStep->GetPostStepPoint()->GetTouchableHandle()->GetVolume();
  // if(volume->GetName()=="mTPCouterwall2_phys"
  //    || volume->GetName()=="mTPCReadoutDisc_phys"
  //    || volume->GetName()=="mTPCReadoutGEMGap_phys"
  //    || volume->GetName()=="mTPCGEMSurf1_phys"
  //    || volume->GetName()=="mTPCGEMDielec_phys"
  //    || volume->GetName()=="mTPCGEMSurf2_phys"
  //    || volume->GetName()=="mTPCGEMGap_phys")
  //   {
  //     theTrack->SetTrackStatus(fStopAndKill);
  //   }

  //fKillTrackAndSecondaries - kill current track and also associated secondaries
  //fSuspend - suspend processing of current track and push it and its secondaries to the stack
  // mTPC in GEMC, anything passing readout GEMS and discs were killed as they were made of kryptonite
  //GetMaterial refers to pre-step point, so hopefully any info already stored in SD is recorded?
  // if(theTrack->GetMaterial()->GetName() == "BonusPCB"){ theTrack->SetTrackStatus(fStopAndKill); }
  // if(track->GetMaterial()->GetName() == "BonusPCB"){ track->SetTrackStatus(fStopAndKillSecondaries); }
  // if(track->GetMaterial()->GetName() == "BonusPCB"){ track->SetTrackStatus(fSuspend); }



  // G4StepPoint *prestep = s->GetPreStepPoint();
  // G4Material *mat = prestep->GetMaterial();
    // G4cout << "Before step, material = " << mat->GetName() << G4endl;
    // G4StepPoint *poststep = s->GetPostStepPoint();
    // mat = poststep->GetMaterial();
    // G4cout << "After step, material = " << mat->GetName() << G4endl;
  // G4MaterialPropertiesTable *mpt = mat->GetMaterialPropertiesTable();

  // if( mpt != NULL ){
  //   //mpt->DumpTable();
  // }
  
  // G4cout << poststep->GetProcessDefinedStep()->GetProcessName() << G4endl;

  //  G4cout << process->GetProcessName() << G4endl;

}


