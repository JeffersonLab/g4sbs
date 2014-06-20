#ifndef __G4SBSBeamlineBuilder_hh
#define __G4SBSBeamlineBuilder_hh

#include "G4SBSComponent.hh"

class G4LogicalVolume;
class GSBSDetectorConstruction;

class G4SBSBeamlineBuilder: public G4SBSComponent {
    public:
	G4SBSBeamlineBuilder(G4SBSDetectorConstruction *);
	~G4SBSBeamlineBuilder();

	void BuildComponent(G4LogicalVolume *);

    private:
	void MakeGEpLead(G4LogicalVolume *);
	void MakeGEnLead(G4LogicalVolume *);
	void MakeGEnClamp(G4LogicalVolume *);
  void MakeSIDISLead( G4LogicalVolume * );

};

#endif//__G4SBSBeamlineBuilder_hh
