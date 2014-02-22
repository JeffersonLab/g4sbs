#include "G4SBSComponent.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4SBSDetectorConstruction.hh"

G4SBSComponent::G4SBSComponent(G4SBSDetectorConstruction *dc):
    fDetCon(dc){
	;
}

G4SBSComponent::~G4SBSComponent(){;}

G4Material *G4SBSComponent::GetMaterial(G4String mat){
    return fDetCon->GetMaterial(mat);
}

G4OpticalSurface *G4SBSComponent::GetOpticalSurface(G4String sur){
    return fDetCon->GetOpticalSurface(sur);
}
