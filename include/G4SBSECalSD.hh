#ifndef G4SBSECalSD_h
#define G4SBSECalSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4SBSECalHit.hh"
#include "G4SBSDetMap.hh"
#include "G4Step.hh"
#include "G4SBSSDTrackOutput.hh"

class G4SBSECalSD : public G4VSensitiveDetector
{
public:
  G4SBSECalSD( G4String name, G4String collname );
  ~G4SBSECalSD();

  void Initialize(G4HCofThisEvent *HC);
  G4bool ProcessHits( G4Step *s, G4TouchableHistory *h );
  void EndOfEvent( G4HCofThisEvent *HC );
  
  void clear();
  void DrawAll();
  void PrintAll();

  G4SBSDetMap detmap;

  // *****
  G4int GetNTimeBins(){ return fNTimeBins; }
  G4double GetTimeWindow(){ return fHitTimeWindow; }
  G4double GetPEThreshold(){ return fPEThreshold; }

  void SetNTimeBins( G4int N ){ fNTimeBins = N; }
  void SetHitTimeWindow( G4double T ){ fHitTimeWindow = T; } 
  void SetPEThreshold( G4double Emin){ fPEThreshold = Emin; }
  // *****
  
  G4SBSSDTrackOutput SDtracks;
  
private:

  G4int    fNTimeBins;       //Number of time bins to record energy deposition after the start of the hit
  G4double fHitTimeWindow;   //Maximum time after the first tracking step of the hit before starting a new hit. 
  G4double fPEThreshold; //Threshold on the minimum summed energy deposition of a hit to record to the output.
  G4SBSECalHitsCollection *hitCollection;
};

#endif
