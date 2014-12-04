#ifndef G4SBSCalSD_h
#define G4SBSCalSD_h 1

#include "G4SBSCalHit.hh"

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include <map>
#include <set>
#include <vector>

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

  map<G4String, map<int, int> > RowMap; //Mapping of SD names to row and column maps (Calorimeters with different kinds of layouts)
  map<G4String, map<int, int> > ColMap; //Mapping of SD names to row and column maps (Calorimeters with different layouts)
  map<G4String, map<int, int> > DepthMap; //Number of levels up (or down?) the volume hierarchy the volume containing the row, col, cell, xcell, ycell info is
  map<G4String, map<int, double> > XMap; //Mapping of "X" coordinate of rows and columns to SD names.
  map<G4String, map<int, double> > YMap; //Mapping of "Y" coordinate of rows and columns to SD names.

private:

  G4SBSCalHitsCollection *hitCollection;
};




#endif

