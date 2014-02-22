#include "G4SBSComponent.hh"

G4SBSComponent::G4SBSComponent(G4SBSDetectorConstruction *dc):
    fDetCon(dc){
	;
}

G4SBSComponent::~G4SBSComponent();

G4SBSComponent::GetMaterial(G4String mat){
    return fDetCon->GetMaterial(mat);
}
