#ifndef G4SBSRICHSD_h
#define G4SBSRICHSD_h 1

#include "G4VSensitiveDetector.hh"
#include "G4SBSRICHHit.hh"
#include "G4SBSDetMap.hh"
#include "G4Step.hh"
#include "G4SBSSDTrackOutput.hh"

class G4SBSRICHSD : public G4VSensitiveDetector
{
public:
  G4SBSRICHSD( G4String name, G4String collname ); //why "colname"? Is this for HitsCollection?
  ~G4SBSRICHSD();

  void Initialize(G4HCofThisEvent *HC);
  G4bool ProcessHits( G4Step *s, G4TouchableHistory *h );
  void EndOfEvent( G4HCofThisEvent *HC );
  
  void clear();
  void DrawAll();
  void PrintAll();

  G4SBSDetMap detmap;

  G4SBSSDTrackOutput SDtracks;

private:
  G4SBSRICHHitsCollection *hitCollection;

};

#endif
