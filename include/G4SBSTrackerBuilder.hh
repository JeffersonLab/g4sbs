#ifndef __G4SBSTrackerBuilder_hh
#define __G4SBSTrackerBuilder_hh

#include "G4SBSComponent.hh"

class G4SBSTrackerBuilder: public G4SBSComponent {
    public:
	G4SBSTrackerBuilder(G4SBSDetectorConstruction *);
	~G4SBSTrackerBuilder();

	void BuildComponent(G4LogicalVolume *);
	void BuildComponent(G4LogicalVolume *, G4RotationMatrix *, G4ThreeVector, G4int, vector<double>, vector<double>, vector<double>);
	

    private:
}

#endif//__G4SBSTrackerBuilder_hh
