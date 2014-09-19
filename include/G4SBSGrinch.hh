#ifndef G4SBSGrinch_hh
#define G4SBSGrinch_hh 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "G4Colour.hh"
#include "G4Polyhedra.hh"
#include "G4Cons.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ElementTable.hh"
#include "G4UnitsTable.hh"
#include "G4Isotope.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4PVPlacement.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4PVParameterised.hh"
#include "G4Sphere.hh"
#include "G4RotationMatrix.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4VSolid.hh"
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"

#include "G4GenericTrap.hh"
#include "G4TessellatedSolid.hh"

#include <vector>
#include <map>
#include "G4UnitsTable.hh"


#include "G4SDManager.hh"
#include "G4PSPassageCellCurrent.hh"
#include "G4SDParticleWithEnergyFilter.hh"

#include "G4SBSComponent.hh"

class G4LogicalVolume;

class G4SBSGrinch : public G4SBSComponent {
	public:
		G4SBSGrinch( G4SBSDetectorConstruction *dc);
		~G4SBSGrinch();

	public:
		void  BuildComponent(G4LogicalVolume *);

                void SetZOffset(G4double off){ fDetOffset = off ;}

	private:
		G4VSolid* ConstructSimple(const G4String& aName, const G4String& aShape, const G4ThreeVector& aFullSize);
		G4LogicalVolume* Hall_log;
		G4Colour GetColor(const G4String& aColor);

	private:
		G4LogicalVolume* GC_Tank_log;
		G4VPhysicalVolume* Tank_phys;

                G4double fDetOffset;
}
;

#endif /*G4SBSGrinch_hh*/
