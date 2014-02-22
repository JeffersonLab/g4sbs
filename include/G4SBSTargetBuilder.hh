#ifndef __G4SBSTargetBuilder_hh
#define __G4SBSTargetBuilder_hh

#include "G4SBSComponent.hh"

class G4SBSTargetBuilder: public G4SBSComponent {
    public:
	G4SBSTargetBuilder(G4DetectorConstruction *);
	~G4SBSTargetBuilder();

	void BuildComponent(G4LogicalVolume *);

	void SetTarget(Targ_t t){fTargType = t;}
	void SetTargLen(double len){ fTargLen = len;}
	void SetTargDen(double den){ fTargDen = den;}

    private:
	double fTargLen;
	double fTargDen;

	Targ_t fTargType;
}

#endif//__G4SBSTargetBuilder_hh
