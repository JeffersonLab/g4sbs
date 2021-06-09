#ifndef G4SBSmTPCSD_h
#define G4SBSmTPCSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4SBSmTPCHit.hh"
#include "G4SBSDetMap.hh"
#include "G4Step.hh"

class G4SBSmTPCSD : public G4VSensitiveDetector
{
public:
  G4SBSmTPCSD( G4String name, G4String collname );
  ~G4SBSmTPCSD();

  void Initialize(G4HCofThisEvent *HC);
  G4bool ProcessHits( G4Step *s, G4TouchableHistory *h );
  void EndOfEvent( G4HCofThisEvent *HC );
  
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
  G4SBSmTPCHitsCollection *hitCollection;
  // from CAL SD output set up
  G4int    fNTimeBins;       //Number of time bins to record energy deposition after the start of the hit
  G4double fHitTimeWindow;   //Maximum time after the first tracking step of the hit before starting a new hit. 
  G4double fEnergyThreshold; //Threshold on the minimum summed energy deposition of a hit to record to the output.

};

#endif
