#ifndef G4SBSCalSD_h
#define G4SBSCalSD_h 1

#include "G4SBSCalHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"

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

  private:
      G4SBSCalHitsCollection *hitCollection;
};




#endif

