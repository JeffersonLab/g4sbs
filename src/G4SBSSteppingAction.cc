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
#include "G4LogicalVolume.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Track.hh"
#include "G4SBSTrackInformation.hh"
#include "G4SBSDetectorConstruction.hh"


G4SBSSteppingAction::G4SBSSteppingAction()
:drawFlag(false)
{
 // new G4SBSSteppingActionMessenger(this);
}
/*
//<<<<<<< HEAD
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
  =======
*/

void G4SBSSteppingAction::Initialize( G4SBSDetectorConstruction *fdc ){
  fSDboundaryVolumes = fdc->SDboundaryVolumes; 
}

void G4SBSSteppingAction::UserSteppingAction(const G4Step *aStep)
{

  // G4cout << "Begin User Stepping Action..." << G4endl;
  
  G4Track *theTrack = aStep->GetTrack();
  if( theTrack->GetTrackStatus() != fAlive ){ return; }
  
  G4String lvname_prestep = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName();
  G4String lvname_poststep = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetName();

  //Here we are interested in the cases when a track ENTERS a "SD boundary volume":
  //So we want the pre-step logical volume to NOT be in the list, and the post-step logical volume to be in the list.
  //Moreover, we want the pre-step logical volume to be the mother volume of the post-step
  // logical volume!
  // std::pair<std::map<G4String,set<G4String> >::iterator,bool> pretest = fSDboundaryVolumes.find( lvname_prestep );
  // std::pair<std::map<G4String,set<G4String> >::iterator,bool> posttest = fSDboundaryVolumes.find( lvname_poststep );

  G4bool prestepismother = ( aStep->GetPostStepPoint()->GetPhysicalVolume()->GetMotherLogical() == 
			     aStep->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume() );
  
  G4bool SDentry = ( lvname_prestep != lvname_poststep && fSDboundaryVolumes.find( lvname_poststep ) != fSDboundaryVolumes.end() );
  
  if( SDentry && prestepismother ){
    // G4cout << "Entered SD boundary volume, prestep volume = " << lvname_prestep
    // 	   << ", poststep volume = " << lvname_poststep << G4endl;
    
    G4SBSTrackInformation *theTrackInfo = (G4SBSTrackInformation*) theTrack->GetUserInformation();

    //Set the track information for all SDs contained within this boundary volume:
    for( set<G4String>::iterator it = fSDboundaryVolumes[lvname_poststep].begin(); it != fSDboundaryVolumes[lvname_poststep].end(); ++it ){
      //SetTrackSDInformation takes care of checking whether the SD is already in the list of SDs
      //associated with this track:
      theTrackInfo->SetTrackSDInformation( *it, theTrack );
    }
    
    //Again, need to copy since parent track info gets destroyed when track goes out of scope:
    //But is this actually true? I don't think so.
    //Here we are actually just modifying the existing track info object of the current track.
    //Let's comment the following lines out and see if it still works:
    //G4SBSTrackInformation *theNewTrackInfo = new G4SBSTrackInformation(theTrackInfo);

    //    theTrack->SetUserInformation(theNewTrackInfo);
  }
}


