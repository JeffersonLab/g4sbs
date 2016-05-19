#ifndef __G4SBSTrackerBuilder_hh
#define __G4SBSTrackerBuilder_hh


#include "G4SBSComponent.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include "sbstypes.hh"

using namespace std;

class G4SBSTrackerBuilder: public G4SBSComponent {
public:
  G4SBSTrackerBuilder(G4SBSDetectorConstruction *);
  ~G4SBSTrackerBuilder();

  void BuildComponent(G4LogicalVolume *);
  void BuildComponent(G4LogicalVolume *, G4RotationMatrix *, G4ThreeVector, unsigned int, vector<double>, vector<double>, vector<double>, G4String );
	

private:
};

#endif//__G4SBSTrackerBuilder_hh
