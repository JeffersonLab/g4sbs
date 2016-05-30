#ifndef G4SBSMWDCSD_h
#define G4SBSMWDCSD_h 1

#include "G4SBSMWDCHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4SBSDetMap.hh"
#include "sbstypes.hh"

#include <map>

using namespace std;

class G4SBSMWDCSD : public G4VSensitiveDetector
{

public:
  G4SBSMWDCSD(G4String name, G4String colname);
  ~G4SBSMWDCSD();

  void Initialize(G4HCofThisEvent*HCE);
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  void EndOfEvent(G4HCofThisEvent*HCE);
  void clear();
  void DrawAll();
  void PrintAll();

  void SetZoffset(double z){ fZoffset = z; }
  
  // map<G4String, int> MWDCTrackerIDs; //Map to associate logical volume names with tracking modules
  // map<int, Arm_t>    MWDCArmIDs;

  G4SBSDetMap detmap;

private:
  G4SBSMWDCHitsCollection *hitCollection;
  double fZoffset;

  
};




#endif

