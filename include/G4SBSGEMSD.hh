#ifndef G4SBSGEMSD_h
#define G4SBSGEMSD_h 1

#include "G4SBSGEMHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"

class G4SBSGEMSD : public G4VSensitiveDetector
{

  public:
      G4SBSGEMSD(G4String name, G4String colname);
      ~G4SBSGEMSD();

      void Initialize(G4HCofThisEvent*HCE);
      G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
      void EndOfEvent(G4HCofThisEvent*HCE);
      void clear();
      void DrawAll();
      void PrintAll();

      void SetZoffset(double z){ fZoffset = z; }

  private:
      G4SBSGEMHitsCollection *hitCollection;
      double fZoffset;
};




#endif

