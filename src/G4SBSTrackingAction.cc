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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
G4SBSTrackingAction::G4SBSTrackingAction()
:G4UserTrackingAction()
{;}

void G4SBSTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  // Create trajectory only for track in tracking region
  G4SBSTrackInformation* trackInfo = 
    (G4SBSTrackInformation*)(aTrack->GetUserInformation());

  // if(trackInfo->GetTrackingStatus() > 0)
  // {
  //   fpTrackingManager->SetStoreTrajectory(true);
  //   fpTrackingManager->SetTrajectory(new G4SBSTrajectory(aTrack));
  // }
  // else
  // { fpTrackingManager->SetStoreTrajectory(false); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 
void G4SBSTrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
  G4TrackVector* secondaries = fpTrackingManager->GimmeSecondaries();
  if(secondaries)
  {
    G4SBSTrackInformation* info = 
      (G4SBSTrackInformation*)(aTrack->GetUserInformation());
    size_t nSeco = secondaries->size();
    if(nSeco>0)
    {
      for(size_t i=0;i<nSeco;i++)
      { 
        G4SBSTrackInformation* infoNew = new G4SBSTrackInformation(info);
        (*secondaries)[i]->SetUserInformation(infoNew);
      }
    }
  }
}


