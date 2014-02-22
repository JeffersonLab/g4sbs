#ifndef __G4SBSComponent_hh
#define __G4SBSComponent_hh

/*  Class to break up SBS representation into
    different sections such as target, BigBite
    arm, etc.  

    Should be pretty self contained. Gets the
    main detector constructor to which it can
    talk to for information, materials, etc.

    Seamus Riordan
    Feb 22, 2014
*/

#include <map>
#include "G4String.hh"

class G4Material;
class G4SensitiveDetectorManager;
class G4LogicalVolume;
class G4OpticalSurface;
class G4SBSDetectorConstruction;

class G4SBSComponent {
    public:
	G4SBSComponent(G4SBSDetectorConstruction * );
	virtual ~G4SBSComponent();
	virtual void BuildComponent(G4LogicalVolume *) = 0;

    protected:
	G4SBSDetectorConstruction *fDetCon;
	G4Material *GetMaterial(G4String);
	G4OpticalSurface*GetOpticalSurface(G4String);
};


#endif//__G4SBSComponent_hh
