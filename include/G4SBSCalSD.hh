#ifndef G4SBSCalSD_h
#define G4SBSCalSD_h 1

#include "G4SBSCalHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4SBSDetMap.hh"

using namespace std;

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;

class G4SBSCalSD : public G4VSensitiveDetector
{

public:
  G4SBSCalSD(G4String name, G4String colname);
  ~G4SBSCalSD();

  void Initialize(G4HCofThisEvent*HCE);
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
  void EndOfEvent(G4HCofThisEvent*HCE);
  void clear();
  void DrawAll();
  void PrintAll();

  G4SBSDetMap detmap;

  G4int GetNTimeBins(){ return fNTimeBins; }
  G4double GetTimeWindow(){ return fHitTimeWindow; }
  G4double GetEnergyThreshold(){ return fEnergyThreshold; }
  
  void SetNTimeBins( G4int N ){ fNTimeBins = N; }
  void SetHitTimeWindow( G4double T ){ fHitTimeWindow = T; } 
  void SetEnergyThreshold( G4double Emin){ fEnergyThreshold = Emin; } 
  
private:

  G4int    fNTimeBins;       //Number of time bins to record energy deposition after the start of the hit
  G4double fHitTimeWindow;   //Maximum time after the first tracking step of the hit before starting a new hit. 
  G4double fEnergyThreshold; //Threshold on the minimum summed energy deposition of a hit to record to the output.
  G4SBSCalHitsCollection *hitCollection;
};




#endif

