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
/// \file runAndEvent/RE01/include/RE01TrackInformation.hh
/// \brief Definition of the RE01TrackInformation class
//
//
//
//Use RE01 track information and Tracking Action classes as template: modify for G4SBS as appropriate.

#ifndef G4SBSTrackInformation_h
#define G4SBSTrackInformation_h 1

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4Allocator.hh"
#include "G4VUserTrackInformation.hh"
#include "G4String.hh"
#include <map>
#include <set>
#include <vector>

using namespace std;

class G4SBSTrackInformation : public G4VUserTrackInformation 
{
public:
  G4SBSTrackInformation();
  G4SBSTrackInformation(const G4Track* aTrack);
  G4SBSTrackInformation(const G4SBSTrackInformation* aTrackInfo);
  virtual ~G4SBSTrackInformation();
   
  inline void *operator new(size_t);
  inline void operator delete(void *aTrackInfo);

  G4SBSTrackInformation& operator =(const G4SBSTrackInformation& right);

  void SetOriginalTrackInformation(const G4Track* aTrack);
  void SetPrimaryTrackInformation(const G4Track* aTrack);
  void SetTrackSDInformation(G4String SDname, const G4Track *aTrack);
  virtual void Print() const;

public:
  inline G4int GetTrackingStatus() const {return fTrackingStatus;}
  inline void  SetTrackingStatus(G4int i) {fTrackingStatus = i;}
  //inline G4int GetSourceTrackID() const {return fSourceTrackID;}
  // inline void  SetSuspendedStepID(G4int i) {fSuspendedStepID = i;}
  // inline G4int GetSuspendedStepID() const {return fSuspendedStepID;}

  inline G4int GetOriginalTrackID() const { return fOriginalTrackID; }
  inline G4int GetOriginalParentID() const { return fOriginalParentID; }
  inline G4ParticleDefinition *GetOriginalDefinition() const { return fOriginalDefinition; }
  inline G4ThreeVector GetOriginalPosition() const { return fOriginalPosition; }
  inline G4ThreeVector GetOriginalMomentum() const { return fOriginalMomentum; }
  inline G4ThreeVector GetOriginalPolarization() const { return fOriginalPolarization; }
  inline G4double GetOriginalEnergy() const { return fOriginalEnergy; }
  inline G4double GetOriginalTime() const { return fOriginalTime; }

  inline G4int GetPrimaryTrackID() const { return fPrimaryTrackID; }
  inline G4ParticleDefinition *GetPrimaryDefinition() const { return fPrimaryDefinition; }
  inline G4ThreeVector GetPrimaryPosition() const { return fPrimaryPosition; }
  inline G4ThreeVector GetPrimaryMomentum() const { return fPrimaryMomentum; }
  inline G4ThreeVector GetPrimaryPolarization() const { return fPrimaryPolarization; }
  inline G4double GetPrimaryEnergy() const { return fPrimaryEnergy; }
  inline G4double GetPrimaryTime() const { return fPrimaryTime; }

private:
  //THIS NOTATION (COPIED FROM EXAMPLE RE01) SUCKS, BTW. "Original" and "Source" sound the same
  //and totally interchangeable and ambiguous.
  //THESE COMMENTS SUCK TOO
  //What do we want to use the track information for in G4SBS anyway?
  // 1. Associate hits in sensitive volumes with primary and/or "source" particles.
  //    This we would usually take to mean the primary particles generated at the vertex.
  //    Considering the four types of sensitive detectors currently in g4sbs, this information would be
  //    mostly redundant for GEMs, but will be useful for RICH, Cal, and "ECal" sensitivity types, which only
  //    get translated into "hit" level information after the production of large numbers of secondaries
  // 2. For simulating recoil nucleon polarimetry. In this case, we are interested not only in the primary particles
  //    generated at the event vertex, but also secondary particles produced in polarization-analyzing
  //    nucleon-nucleus scattering reactions, including neutrons.
  // 3. For beam background studies: here we are often interested in locating the sources of backgrounds
  // For case (1), we are not only interested in the origin information of the track, but also where
  // it first enters the relevant detection volume. For example, in the BigBite detector stack, we are
  // usually looking for a single primary electron track producing signals in all the detectors
  // We want to be sufficiently generic to preserve all relevant information for all "hits" in sensitive detectors.
  // We also want to be efficient in terms of memory usage.
  // So this class, together with G4SBSSteppingAction, should basically be interested in recording the
  // properties of particles when the beginning (or end) of the track occurs within a handful of situations:
  //  1) secondaries produced in the TARGET (for beam background simulations) or ANALYZERS (for recoil nucleon polarimetry)
  //  2) particles ENTERING the "region" of sensitive detectors (trackers or calorimeters), this will likely require this class to interact with "G4SBSSteppingAction"
  //Strategy: For each hit in each sensitive detector, we want to record the "source track" ID and the "primary track" ID. To do so, we need to define a logical volume for each sensitive detector, which could be but is not necesarily the
  // sensitive volume, upon entry to which we record the source track information.
  // For sensitive detectors whose logical volumes are placed many times, e.g., calorimeters, typically the relevant
  // volume boundary crossing for source track information recording is a mother volume enclosing all
  // individual segments of the detector.

  //Is this the right place to store a map of sensitive detectors traversed by the track?
  //Maybe. Except this gets deleted at the end of an event. So should we store Trajectories for
  //Tracks crossing sensitive volume boundaries? Maybe. But the standard G4Trajectory may not be adequate for this.
  //What we actually want is for the source track information to propagate down to the SD hit information.
  //So a map of sensitive detectors entered and the track information as they cross the relevant boundaries
  //seems appropriate:
  
  // Information of the primary track at the primary vertex:
  // As written, this is actually the track information of every track when it is produced.
  G4int                 fOriginalTrackID;  // Track ID of original particle
  G4int                 fOriginalParentID; //IF applicable
  G4ParticleDefinition* fOriginalDefinition;
  G4ThreeVector         fOriginalPosition;
  G4ThreeVector         fOriginalMomentum;
  G4ThreeVector         fOriginalPolarization;
  G4double              fOriginalEnergy;
  G4double              fOriginalTime;

  // This is the primary track responsible for the current track at the primary vertex:
  G4int                 fPrimaryTrackID;  // Track ID of primary particle
  G4ParticleDefinition* fPrimaryDefinition;
  G4ThreeVector         fPrimaryPosition;
  G4ThreeVector         fPrimaryMomentum;
  G4ThreeVector         fPrimaryPolarization;
  G4double              fPrimaryEnergy;
  G4double              fPrimaryTime;

public:
  set<G4String>                       fSDlist; //list of SD names associated with this track
  map<G4String,G4ParticleDefinition*> fSDDefinition; //Particle ID of track when entering SD
  map<G4String,G4int>                 fSDTrackID;    //Track ID of track entering SD
  map<G4String,G4int>                 fSDParentID;   //Parent ID of track entering SD
  map<G4String,G4ThreeVector>         fSDPosition; //Global position of track when entering SD.
  map<G4String,G4ThreeVector>         fSDMomentum; //Global three-momentum of track when entering SD.
  map<G4String,G4ThreeVector>         fSDPolarization; //Global Polarization of track when entering SD
  map<G4String,G4ThreeVector>         fSDVertexPosition; //Vertex position of track that enters SD boundary volume
  map<G4String,G4ThreeVector>         fSDVertexDirection; //Momentum direction at vertex of track that enters SD boundary volume
  map<G4String,G4double>              fSDVertexKineticEnergy; //Kinetic energy at vertex of track entering SD volume
  map<G4String,G4double>              fSDEnergy; //Total energy of track when entering SD
  map<G4String,G4double>              fSDTime;   //Global time of track when entering SD

  // //"Source" track information means information about a track as it enters the "region of interest" of a detector
  // G4int                 fSourceTrackID;
  // G4ParticleDefinition* fSourceDefinition;
  // G4ThreeVector         fSourcePosition;
  // G4ThreeVector         fSourceMomentum;
  // G4ThreeVector         fSourcePolarization;
  // G4double              fSourceEnergy;
  // G4double              fSourceTime;
  // //Source track information when it was produced:
  // G4ThreeVector         fSourceVertexPosition;
  // G4ThreeVector         fSourceVertexMomentum;
  // G4ThreeVector         fSourceVertexEnergy;
  // G4ThreeVector         fSourceVertexTime; 
  // //G4int                 fSuspendedStepID;

  //"Tracking Status": Flag to indicate if particle was produced in a "region of interest"
  // 1: TARGET
  // 2: ANALYZER
  // 0: all other volumes
private:
  G4int                 fTrackingStatus; 
};

extern G4ThreadLocal
 G4Allocator<G4SBSTrackInformation> * aTrackInformationAllocator;

inline void* G4SBSTrackInformation::operator new(size_t)
{
  if(!aTrackInformationAllocator)
    aTrackInformationAllocator = new G4Allocator<G4SBSTrackInformation>;
  return (void*)aTrackInformationAllocator->MallocSingle();
}

inline void G4SBSTrackInformation::operator delete(void *aTrackInfo)
{ aTrackInformationAllocator->FreeSingle((G4SBSTrackInformation*)aTrackInfo);}

#endif
