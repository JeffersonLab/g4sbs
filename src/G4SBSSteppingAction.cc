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
///  new G4SBSSteppingActionMessenger(this);
}

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


