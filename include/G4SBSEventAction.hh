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

#include "G4SBSBeamDiffuserSD.hh"
#include "G4SBSBDHit.hh"
#include "G4SBSBDoutput.hh"

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

  //Although these functions don't directly modify the "SDTrackOutput", we pass them by const reference
  //to avoid the overhead of copying:
  void FillGEMData( const G4Event*, G4SBSGEMHitsCollection*, G4SBSGEMoutput &, G4SBSSDTrackOutput & );
  void FillCalData( const G4Event*, G4SBSCalHitsCollection*, G4SBSCALoutput &, G4SBSSDTrackOutput & );
  void FillRICHData( const G4Event*, G4SBSRICHHitsCollection*, G4SBSRICHoutput &, G4SBSSDTrackOutput & );
  void FillTrackData( G4SBSGEMoutput, G4SBSTrackerOutput & );
  void FillECalData( G4SBSECalHitsCollection*, G4SBSECaloutput &, G4SBSSDTrackOutput & );
  void FillBDData(const G4Event *evt,G4SBSBDHitsCollection *hc,G4SBSBDoutput &out); // for the Beam Diffuser (BD)

  //map<G4String, G4VSensitiveDetector*> SDlist; //List of all sensitive detectors in the run. 
  set<G4String> SDlist;
  map<G4String, G4SBS::SDet_t> SDtype;
  //map<G4String, G4SBS::Arm_t> SDarm;
  void SetEventStatusEvery(G4int n) { fEventStatusEvery = n; };

private:
  //Hit collection IDs:
  G4int gemCollID, hcalCollID, bbcalCollID, RICHCollID, ECalCollID;
  // G4int bdCollID; // for the beam diffuser (bd)  

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

  G4int fhistogram_index;
  
public:
};

#endif
  
