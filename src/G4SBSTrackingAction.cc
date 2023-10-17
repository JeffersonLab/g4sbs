//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file runAndEvent/RE01/src/RE01TrackingAction.cc
/// \brief Implementation of the RE01TrackingAction class
//
//

#include "G4SBSTrackingAction.hh"
//#include "G4SBSTrajectory.hh"
#include "G4SBSTrackInformation.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4Trajectory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
G4SBSTrackingAction::G4SBSTrackingAction()
:G4UserTrackingAction()
{;}

void G4SBSTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // When a new track is created, ask: 
  // 1. Is this a primary particle? If so, create the User information
  // 2. If it is a secondary, was it created in a logical volume of interest (target, analyzer)?
  //    a) IF SO: set "original"
  //    b) IF NOT: do nothing:

  //G4cout << "Starting tracking for track ID = " << aTrack->GetTrackID() << ", parent track ID = " << aTrack->GetParentID() << G4endl;
  
  G4SBSTrackInformation *trackInfo;

  // aTrack->SetUserInformation( trackInfo );

  G4String lvname = aTrack->GetVolume()->GetLogicalVolume()->GetName();

  
  
  if( aTrack->GetParentID() == 0 ){ //Primary particle: 
    trackInfo = new G4SBSTrackInformation( aTrack ); //sets 

    //Declare a pointer to track and assign it to the function argument,
    //as a way to get around compiler error due to the fact that the argument is "const".
    G4Track *theTrack = (G4Track*) aTrack;
    theTrack->SetUserInformation( trackInfo );
    
    // aTrack->SetUserInformation( trackInfo );
    // G4cout << "New primary track created track ID = " << aTrack->GetTrackID() << G4endl
    // 	   << " Parent ID = " << aTrack->GetParentID() << G4endl
    // 	   << " LV name = \"" << lvname << "\"" << G4endl
    // 	   << " Position = " << aTrack->GetPosition()/CLHEP::cm << G4endl
    // 	   << " Momentum = " << aTrack->GetMomentum()/CLHEP::GeV << G4endl
    // 	   << " Polarization = " << aTrack->GetPolarization() << G4endl
    // 	   << " Total Energy = " << aTrack->GetTotalEnergy()/CLHEP::GeV << G4endl
    // 	   << " Global time = " << aTrack->GetGlobalTime() << G4endl << G4endl;
    
  } else { //secondary: grab the track information that will (presumably) have already been copied from the parent track:
    trackInfo = (G4SBSTrackInformation*) aTrack->GetUserInformation(); //this preserves the primary track info. The "origin track info" 

    // G4cout << "This is a secondary: (PID, parent PID) of track " << aTrack->GetTrackID() << " = (" << aTrack->GetDefinition()->GetPDGEncoding() << ", "
    // 	   << trackInfo->GetOriginalParentPID()->GetPDGEncoding() << ")" << G4endl;
    
  }

  // Whether this is a primary or secondary track, overwrite the "original track" info and set the
  // appropriate tracking status IFF this track was produced in one of the "target" or "analyzer" volumes
  
  if( fTargetVolumes.find( lvname ) != fTargetVolumes.end() ){
    trackInfo->SetTrackingStatus( 1 );
    trackInfo->SetOriginalTrackInformation( aTrack );
    // if( aTrack->GetParentID() != 0 ){
    //   // G4cout << "New secondary track created in target, track ID = " << aTrack->GetTrackID() << G4endl
    //   // 	     << " Parent ID = " << aTrack->GetParentID() << G4endl
    //   // 	     << " PID = " << aTrack->GetDefinition()->GetPDGEncoding() << G4endl
    //   // 	     << " Parent PID = " << trackInfo->GetOriginalParentPID()->GetPDGEncoding() << G4endl
    //   // 	     << " LV name = \"" << lvname << "\"" << G4endl
    //   // 	     << " Position = " << aTrack->GetPosition()/CLHEP::cm << G4endl
    //   // 	     << " Momentum = " << aTrack->GetMomentum()/CLHEP::GeV << G4endl
    //   // 	     << " Polarization = " << aTrack->GetPolarization() << G4endl
    //   // 	     << " Total Energy = " << aTrack->GetTotalEnergy()/CLHEP::GeV << G4endl
    //   // 	     << " Global time = " << aTrack->GetGlobalTime() << G4endl << G4endl;
    // }
  }

  if( fAnalyzerVolumes.find( lvname ) != fAnalyzerVolumes.end() ){
    trackInfo->SetTrackingStatus( 2 );
    trackInfo->SetOriginalTrackInformation( aTrack );
    //if( aTrack->GetParentID() != 0 ){
      // G4cout << "New secondary track created in analyzer, track ID = " << aTrack->GetTrackID() << G4endl
      // 	     << " Parent ID = " << aTrack->GetParentID() << G4endl
      // 	     << " PID = " << aTrack->GetDefinition()->GetPDGEncoding() << G4endl
      // 	     << " LV name = \"" << lvname << "\"" << G4endl
      // 	     << " Position = " << aTrack->GetPosition() << G4endl
      // 	     << " Momentum = " << aTrack->GetMomentum() << G4endl
      // 	     << " Polarization = " << aTrack->GetPolarization() << G4endl
      // 	     << " Total Energy = " << aTrack->GetTotalEnergy()/CLHEP::GeV << G4endl
      // 	     << " Global time = " << aTrack->GetGlobalTime() << G4endl << G4endl;
    //}
  }
  
  //We need to use the copy constructor to create a new instance of the track info because when the parent
  //track gets deleted, its UserTrackInformation object gets deleted too:
  //Actually, I'm not sure this is true at the beginning of tracking a new track:
  //G4SBSTrackInformation *trackInfoNew = new G4SBSTrackInformation(trackInfo);
  
  //After appropriately modifying trackInfo, set the usertrackinformation:
  //if( aTrack->GetParentID() == 0 ) aTrack->SetUserInformation( trackInfo );
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void G4SBSTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  //Information that we want to copy from all primaries to all secondaries:
  //1. Primary vertex particle information
  //2. "Source" track information
  //3. "Origin" track information
  //4. Parent Particle ID: 

  //G4cout << "Finishing up tracking for track ID = " << aTrack->GetTrackID() << G4endl;
  
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if(secondaries)
  {
    G4SBSTrackInformation* info = 
      (G4SBSTrackInformation*)(aTrack->GetUserInformation());

    //aTrack here represents the primary, and we want to copy the PID of THIS track to the parent PID of the secondaries
    G4ParticleDefinition *ParentPID = aTrack->GetDefinition();
    
    size_t nSeco = secondaries->size();
    if(nSeco>0)
    {
      for(size_t i=0;i<nSeco;i++)
      {
	// This default behavior copied from example RE01 is acceptable
	// it copies the track information to all secondaries.
	// The preusertrackingaction handles modifications to "OriginalTrack" when the secondaries are
	// produced in volumes of interest ("targets" or "analyzers")
        G4SBSTrackInformation* infoNew = new G4SBSTrackInformation(info); 

	//if( (*secondaries)[i]->GetParentID() == 0 ) G4cout << "Secondary particle with parent ID of zero, TID == " << (*secondaries)[i]->GetTrackID() << G4endl;

	//We set the generic "Parent PID" regardless of the volume where the secondary is produced
	infoNew->SetParentPID( ParentPID );
	//We need to add a check HERE so that we only modify the "parent PID" information if we are in a
	// "target" or "analyzer" volume. But will this work? 
	//Check logical volume:

	G4String lvname = (*secondaries)[i]->GetVolume()->GetLogicalVolume()->GetName();

	if( fTargetVolumes.find( lvname ) != fTargetVolumes.end() ||
	    fAnalyzerVolumes.find( lvname ) != fAnalyzerVolumes.end() ){
	  infoNew->SetOriginalParentPID( ParentPID ); //this copies the Particle ID of the immediate parent particle of each secondary to the "original track" parent PID info
	}
	//Copy track information to all secondaries:
	
        (*secondaries)[i]->SetUserInformation(infoNew); //now the question is, does this nuke the information from the PreUserTrackingAction? -->NO
	// G4cout << "Copying user track information to secondary track number " << i << " of parent track "
	//        << (*secondaries)[i]->GetParentID() << ", parent PID = " << ParentPID->GetPDGEncoding() << ", "
	//        << ", child PID = " << (*secondaries)[i]->GetDefinition()->GetPDGEncoding() << G4endl;
      }
    }
  }

  //G4cout << "... done. " << G4endl;
}

void G4SBSTrackingAction::Initialize( G4SBSDetectorConstruction *fdc ){ //Copy the information ONCE at beginning of run:
  SetTargetVolumes( fdc->GetTargetVolumes() );
  SetAnalyzerVolumes( fdc->GetAnalyzerVolumes() );
}


