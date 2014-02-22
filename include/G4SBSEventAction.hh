
#ifndef G4SBSEventAction_h
#define G4SBSEventAction_h 1

#include "TBuffer.h"
#include "TString.h"
#include "TMatrixTBase.h"

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "G4SBSRICHHit.hh"
#include "G4SBSRICHoutput.hh"

#define __MAXGEM 100

class G4Event;
class G4SBSIO;
class G4SBSEventGen;

class G4SBSEventAction : public G4UserEventAction
{
public:
  G4SBSEventAction();
  virtual ~G4SBSEventAction();
  
public:
  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);
  
  void SetIO( G4SBSIO *io ){ fIO = io; }
  void SetEvGen( G4SBSEventGen *gen ){ fevgen = gen; }
  void SetGEMRes( double r ){ fGEMres = r; }
  
  void LoadSigmas(const char *filename);
  
  void MapTracks(const G4Event *);
  void FillRICHData( const G4Event*, G4SBSRICHHitsCollection*, G4SBSRICHoutput &);

private:
  //Hit collection IDs:
  G4int gemCollID, hcalCollID, bbcalCollID, RICHCollID;
  
  //maps associating trajectory index with track ID, parent track ID and particle ID:
  //For Track IDs, we only need to use vector, because there is a one-to-one mapping between 
  //trajectory number and track ID. 
  //For mother track IDs, we use map, because many tracks can have the same mother. 
  //This map will be an associative mapping between ***TRACK*** ID and mother ID, as opposed to 
  //Parent ID and mother ID:
  map<G4int,G4int> TrajectoryIndex;
  map<G4int,G4int> MotherTrackIDs;

  double fGEMres;
  
  G4SBSIO *fIO;
  G4SBSEventGen *fevgen;
  
  double fGEMsigma[__MAXGEM];
  
public:
};

#endif

    
