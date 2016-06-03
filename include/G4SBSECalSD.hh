#ifndef G4SBSECalSD_h
#define G4SBSECalSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4SBSECalHit.hh"
#include "G4SBSDetMap.hh"
#include "G4Step.hh"

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

private:
  G4SBSECalHitsCollection *hitCollection;
  G4String fSDname;
  void SetName(G4String name) {fSDname = name;} 
  G4String GetSDName() {return fSDname;}
};

#endif
