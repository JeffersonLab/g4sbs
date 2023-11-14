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
/// \file runAndEvent/G4SBS/src/G4SBSTrackInformation.cc
/// \brief Implementation of the G4SBSTrackInformation class
//
//
// Using RE01 track information class as a template and modifying for G4SBS as appropriate

#include "G4SBSTrackInformation.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"


G4ThreadLocal G4Allocator<G4SBSTrackInformation> *
aTrackInformationAllocator = 0;

//default constructor: set everything to zero, clear out "SD maps"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4SBSTrackInformation::G4SBSTrackInformation()
  : G4VUserTrackInformation()
{
  fTrackingStatus = 0;
  
  fOriginalTrackID = 0;
  fOriginalParentID = 0;
  fOriginalDefinition = 0;
  fOriginalParentPID = 0;
  fOriginalPosition = G4ThreeVector(0.,0.,0.);
  fOriginalMomentum = G4ThreeVector(0.,0.,0.);
  fOriginalPolarization = G4ThreeVector(0.,0.,0.);
  fOriginalEnergy = 0.;
  fOriginalTime = 0.;

  fParentPID = 0;
  
  fPrimaryTrackID = 0;
  fPrimaryDefinition = 0;
  fPrimaryPosition = G4ThreeVector(0.,0.,0.);
  fPrimaryMomentum = G4ThreeVector(0.,0.,0.);
  fPrimaryPolarization = G4ThreeVector(0.,0.,0.);
  fPrimaryEnergy = 0.;
  fPrimaryTime = 0.;

  //fNbounce = 0;
  //fPIDbounce.clear();
  
  fSDlist.clear();
  fSDDefinition.clear();
  fSDParentPID.clear();
  fSDTrackID.clear();
  fSDPosition.clear();
  fSDMomentum.clear();
  fSDPolarization.clear();
  fSDEnergy.clear();
  fSDTime.clear();
  fSDVertexPosition.clear();
  fSDVertexDirection.clear();
  fSDVertexKineticEnergy.clear();

  //fSDNbounce.clear();
  //fSDPIDbounce.clear();
}

//This G4Track based constructor will typically get invoked only when new tracks are created:
//Should we clear the SD map info here? Maybe...
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4SBSTrackInformation::G4SBSTrackInformation(const G4Track* aTrack)
  : G4VUserTrackInformation()
{
  
  fTrackingStatus = 0;
  
  fOriginalTrackID = aTrack->GetTrackID();
  fOriginalParentID = aTrack->GetParentID();
  fOriginalDefinition = aTrack->GetDefinition();
  fOriginalParentPID = aTrack->GetDefinition(); //When we invoke the G4Track-based constructor we initially set the Original Parent PID and parent pid to the PID of the track itself. 
  fOriginalPosition = aTrack->GetPosition();
  fOriginalMomentum = aTrack->GetMomentum();
  fOriginalPolarization = aTrack->GetPolarization();
  fOriginalEnergy = aTrack->GetTotalEnergy();
  fOriginalTime = aTrack->GetGlobalTime();

  fParentPID = aTrack->GetDefinition(); 
  
  if( aTrack->GetParentID() == 0 ){ //Primary particle:
    fPrimaryTrackID = aTrack->GetTrackID();
    fPrimaryDefinition = aTrack->GetDefinition();
    fPrimaryPosition = aTrack->GetPosition();
    fPrimaryMomentum = aTrack->GetMomentum();
    fPrimaryPolarization = aTrack->GetPolarization();
    fPrimaryEnergy = aTrack->GetTotalEnergy();
    fPrimaryTime   = aTrack->GetGlobalTime();
    //fNbounce = 0;
    //fPIDbounce.clear();
    //fPIDbounce.push_back( fPrimaryDefinition->GetPDGEncoding() );
  }
  
  fSDlist.clear();
  fSDDefinition.clear();
  fSDTrackID.clear();
  fSDParentID.clear();
  fSDParentPID.clear();
  fSDPosition.clear();
  fSDMomentum.clear();
  fSDPolarization.clear();
  fSDEnergy.clear();
  fSDTime.clear();
  fSDVertexPosition.clear();
  fSDVertexDirection.clear();
  fSDVertexKineticEnergy.clear();
  //fSDNbounce.clear();
  //fSDPIDbounce.clear();
 
}

//COPY constructor: just copy everything
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4SBSTrackInformation
::G4SBSTrackInformation(const G4SBSTrackInformation* aTrackInfo)
  : G4VUserTrackInformation()
{
  fOriginalTrackID = aTrackInfo->fOriginalTrackID;
  fOriginalParentID = aTrackInfo->fOriginalParentID;
  fOriginalParentPID = aTrackInfo->fOriginalParentPID;
  fOriginalDefinition = aTrackInfo->fOriginalDefinition;
  fOriginalPosition = aTrackInfo->fOriginalPosition;
  fOriginalMomentum = aTrackInfo->fOriginalMomentum;
  fOriginalPolarization = aTrackInfo->fOriginalPolarization;
  fOriginalEnergy = aTrackInfo->fOriginalEnergy;
  fOriginalTime = aTrackInfo->fOriginalTime;

  fParentPID = aTrackInfo->fParentPID;
  
  fTrackingStatus = aTrackInfo->fTrackingStatus;

  fPrimaryTrackID = aTrackInfo->fPrimaryTrackID;
  fPrimaryDefinition = aTrackInfo->fPrimaryDefinition;
  fPrimaryPosition = aTrackInfo->fPrimaryPosition;
  fPrimaryMomentum = aTrackInfo->fPrimaryMomentum;
  fPrimaryPolarization = aTrackInfo->fPrimaryPolarization;
  fPrimaryEnergy = aTrackInfo->fPrimaryEnergy;
  fPrimaryTime = aTrackInfo->fPrimaryTime;

  //fNbounce = aTrackInfo->fNbounce;
  //fPIDbounce = aTrackInfo->fPIDbounce;
  
  fSDlist     = aTrackInfo->fSDlist;
  fSDDefinition = aTrackInfo->fSDDefinition;
  fSDTrackID    = aTrackInfo->fSDTrackID;
  fSDParentID   = aTrackInfo->fSDParentID;
  fSDParentPID = aTrackInfo->fSDParentPID;
  fSDPosition = aTrackInfo->fSDPosition;
  fSDMomentum = aTrackInfo->fSDMomentum;
  fSDPolarization = aTrackInfo->fSDPolarization;
  fSDEnergy = aTrackInfo->fSDEnergy;
  fSDTime = aTrackInfo->fSDTime;
  fSDVertexPosition = aTrackInfo->fSDVertexPosition;
  fSDVertexDirection = aTrackInfo->fSDVertexDirection;
  fSDVertexKineticEnergy = aTrackInfo->fSDVertexKineticEnergy;
  //fSDNbounce = aTrackInfo->fSDNbounce;
  //fSDPIDbounce = aTrackInfo->fSDPIDbounce;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4SBSTrackInformation::~G4SBSTrackInformation()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// Assignment operator: Just copy everything
G4SBSTrackInformation& G4SBSTrackInformation
::operator =(const G4SBSTrackInformation& aTrackInfo)
{
  fOriginalTrackID = aTrackInfo.fOriginalTrackID;
  fOriginalParentID = aTrackInfo.fOriginalParentID;
  fOriginalParentPID = aTrackInfo.fOriginalParentPID;
  fOriginalDefinition = aTrackInfo.fOriginalDefinition;
  fOriginalPosition = aTrackInfo.fOriginalPosition;
  fOriginalMomentum = aTrackInfo.fOriginalMomentum;
  fOriginalPolarization = aTrackInfo.fOriginalPolarization;
  fOriginalEnergy = aTrackInfo.fOriginalEnergy;
  fOriginalTime = aTrackInfo.fOriginalTime;

  fParentPID = aTrackInfo.fParentPID;
  
  fTrackingStatus = aTrackInfo.fTrackingStatus;

  //  fNbounce = aTrackInfo.fNbounce;
  //fPIDbounce = aTrackInfo.fPIDbounce;
  
  fPrimaryTrackID = aTrackInfo.fPrimaryTrackID;
  fPrimaryDefinition = aTrackInfo.fPrimaryDefinition;
  fPrimaryPosition = aTrackInfo.fPrimaryPosition;
  fPrimaryMomentum = aTrackInfo.fPrimaryMomentum;
  fPrimaryPolarization = aTrackInfo.fPrimaryPolarization;
  fPrimaryEnergy = aTrackInfo.fPrimaryEnergy;
  fPrimaryTime = aTrackInfo.fPrimaryTime;

  fSDlist     = aTrackInfo.fSDlist;
  fSDDefinition = aTrackInfo.fSDDefinition;
  fSDParentPID = aTrackInfo.fSDParentPID;
  fSDTrackID    = aTrackInfo.fSDTrackID;
  fSDParentID   = aTrackInfo.fSDParentID;
  fSDPosition = aTrackInfo.fSDPosition;
  fSDMomentum = aTrackInfo.fSDMomentum;
  fSDPolarization = aTrackInfo.fSDPolarization;
  fSDEnergy = aTrackInfo.fSDEnergy;
  fSDTime = aTrackInfo.fSDTime;
  fSDVertexPosition = aTrackInfo.fSDVertexPosition;
  fSDVertexDirection = aTrackInfo.fSDVertexDirection;
  fSDVertexKineticEnergy = aTrackInfo.fSDVertexKineticEnergy;
  //fSDNbounce = aTrackInfo.fSDNbounce;
  return *this;
}

//Utility methods:
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4SBSTrackInformation::SetPrimaryTrackInformation(const G4Track* aTrack)
{
  fPrimaryTrackID = aTrack->GetTrackID();
  fPrimaryDefinition = aTrack->GetDefinition();
  fPrimaryPosition = aTrack->GetPosition();
  fPrimaryMomentum = aTrack->GetMomentum();
  fPrimaryPolarization = aTrack->GetPolarization();
  fPrimaryEnergy = aTrack->GetTotalEnergy();
  fPrimaryTime = aTrack->GetGlobalTime();
}

void G4SBSTrackInformation::SetOriginalTrackInformation(const G4Track* aTrack)
{
  fOriginalTrackID = aTrack->GetTrackID();
  fOriginalParentID = aTrack->GetParentID();
  //fOriginalParentPID = aTrack->GetParentPID();
  fOriginalDefinition = aTrack->GetDefinition();
  if( aTrack->GetParentID() == 0 ) fOriginalParentPID = aTrack->GetDefinition(); //if this is a primary particle, set "parent PID" to the PID of the particle itself, to avoid edge cases with erroneous parent PID presumably left over from previous events and/or tracks
  fOriginalPosition = aTrack->GetPosition();
  fOriginalMomentum = aTrack->GetMomentum();
  fOriginalPolarization = aTrack->GetPolarization();
  fOriginalEnergy = aTrack->GetTotalEnergy();
  fOriginalTime = aTrack->GetGlobalTime();
}

//This will typically get invoked in G4SBS stepping action, when the track enters a relevant SD volume:
void G4SBSTrackInformation::SetTrackSDInformation(G4String SDname, const G4Track *aTrack ){

  std::pair<std::set<G4String>::iterator,bool> newSD = fSDlist.insert( SDname ); 
  if( newSD.second ){ //Only do anything if this is the first instance of this track crossing this SD boundary:
    fSDDefinition[SDname] = aTrack->GetDefinition();
    //fSDParentPID[SDname] = aTrack->GetDefinition();
    fSDTrackID[SDname] = aTrack->GetTrackID();
    fSDParentID[SDname] = aTrack->GetParentID();
    fSDPosition[SDname] = aTrack->GetPosition();
    fSDMomentum[SDname] = aTrack->GetMomentum();
    fSDPolarization[SDname] = aTrack->GetPolarization();
    fSDEnergy[SDname] = aTrack->GetTotalEnergy();
    fSDTime[SDname] = aTrack->GetGlobalTime();
    fSDVertexPosition[SDname] = aTrack->GetVertexPosition();
    fSDVertexDirection[SDname] = aTrack->GetVertexMomentumDirection();
    fSDVertexKineticEnergy[SDname] = aTrack->GetVertexKineticEnergy();
    
    G4SBSTrackInformation *info = (G4SBSTrackInformation *) ( aTrack->GetUserInformation() );

    fSDParentPID[SDname] = info->GetParentPID(); //If this is a primary track, GetOriginalParentPID() should return the primary track PID. If this is a secondary, then GetParentPID() should return the PID of the immediate parent of the track.

    //fSDNbounce[SDname] = info->GetNbounce();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4SBSTrackInformation::Print() const
{
  G4cout 
    << "Primary track ID " << fPrimaryTrackID << " (" 
    << fPrimaryDefinition->GetParticleName() << ","
    << fPrimaryEnergy/GeV << "[GeV]) at " << fPrimaryPosition << G4endl;
  G4cout
    << "Original primary track ID " << fOriginalTrackID << " (" 
    << fOriginalDefinition->GetParticleName() << ","
    << fOriginalEnergy/GeV << "[GeV])" << G4endl;
}

