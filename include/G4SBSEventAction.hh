#ifndef G4SBSEventAction_h
#define G4SBSEventAction_h 1

#include "TBuffer.h"
#include "TString.h"
#include "TMatrixTBase.h"

#include "G4UserEventAction.hh"
#include "globals.hh"

#include "G4SBSRICHHit.hh"
#include "G4SBSRICHoutput.hh"
#include "G4SBSCalHit.hh"
#include "G4SBSCALoutput.hh"
#include "G4SBSECalHit.hh"
#include "G4SBSECaloutput.hh"
#include "G4SBSGEMHit.hh"
#include "G4SBSGEMoutput.hh"
#include "G4SBSTrackerOutput.hh"
#include "G4SBSParticleOutput.hh"
#include "sbstypes.hh"
#include <set> 

using namespace std;

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
  void SetTreeFlag( G4int f ){ fTreeFlag = f; }

  void LoadSigmas(const char *filename);
  
  void MapTracks(const G4Event *);
  void FillGEMData( const G4Event*, G4SBSGEMHitsCollection*, G4SBSGEMoutput & );
  void FillCalData( const G4Event*, G4SBSCalHitsCollection*, G4SBSCALoutput & );
  void FillRICHData( const G4Event*, G4SBSRICHHitsCollection*, G4SBSRICHoutput & );
  void FillTrackData( G4SBSGEMoutput, G4SBSTrackerOutput & );
  void FillECalData( G4SBSECalHitsCollection*, G4SBSECaloutput & );

  //map<G4String, G4VSensitiveDetector*> SDlist; //List of all sensitive detectors in the run. 
  set<G4String> SDlist;
  map<G4String, SDet_t> SDtype;
  //map<G4String, Arm_t> SDarm;
  void SetEventStatusEvery(G4int n) { fEventStatusEvery = n; };

private:
  //Hit collection IDs:
  G4int gemCollID, hcalCollID, bbcalCollID, RICHCollID, ECalCollID;

  //maps associating trajectory index with track ID, parent track ID and particle ID
  //For mother track IDs, we use map, because many tracks can have the same mother. 
  //This map will be an associative mapping between ***TRACK*** ID and mother ID, as opposed to 
  //Parent ID and mother ID:
  map<G4int,G4int> TrajectoryIndex;
  map<G4int,G4int> MotherTrackIDs;

  double fGEMres;
  
  G4SBSIO *fIO;
  G4SBSEventGen *fevgen;
  
  double fGEMsigma[__MAXGEM];

  G4int fTreeFlag;
  G4int fEventStatusEvery; //< Print event status at every N-entries
  
public:
};

#endif
  
