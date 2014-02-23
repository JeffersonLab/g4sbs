#include "G4SBSBeamlineBuilder.hh"

#include "G4SBSHArmBuilder.hh"
#include "G4SBSDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"

G4SBSBeamlineBuilder::G4SBSBeamlineBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
}

G4SBSBeamlineBuilder::~G4SBSBeamlineBuilder(){;}

void G4SBSBeamlineBuilder::BuildComponent(G4LogicalVolume *worldlog){
    Targ_t targtype = fDetCon->fTargType;
    
    //Material definition moved to "ConstructMaterials":


    //G4Material* aluminum = new G4Material("Aluminum", 13., 26.98*g/mole, density=2.7*g/cm3);
    // G4Material *aluminum = GetMaterial("Aluminum");
    // G4Material* vacuum = GetMaterial("Vacuum");

    double swallrad = 1.143*m/2;
    double beamheight = 10.0*12*2.54*cm; // 10 feet off the ground

    // Stainless
    G4double ent_len = 10*m;
    G4double ent_rin = 31.75*mm;
    G4double ent_rou = ent_rin+0.120*mm;

    G4Tubs *ent_tube = new G4Tubs("ent_tube", ent_rin, ent_rou, ent_len/2, 0.*deg, 360.*deg );
    G4Tubs *ent_vac  = new G4Tubs("ent_vac", 0.0, ent_rin, ent_len/2, 0.*deg, 360.*deg );

    G4LogicalVolume *entLog = new G4LogicalVolume(ent_tube, GetMaterial("Stainless"), "ent_log", 0, 0, 0);
    G4LogicalVolume *entvacLog = new G4LogicalVolume(ent_vac, GetMaterial("Vacuum"), "entvac_log", 0, 0, 0);


    if( targtype == kH2 || targtype == k3He || targtype == kNeutTarg ){
	// gas target -  1.5m in air
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-1.5*m), entLog, "ent_phys", worldlog, false,0);
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-1.5*m), entvacLog, "entvac_phys", worldlog,false,0);

	// Add in Be window
	G4double winthick = 0.0127*cm;

	G4Tubs *ent_win = new G4Tubs("ent_win", 0.0, ent_rin, winthick/2, 0.*deg, 360.*deg );
	G4LogicalVolume *ent_winlog = new G4LogicalVolume(ent_win, GetMaterial("Beryllium"), "entwin_log", 0, 0, 0);
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, ent_len/2-winthick/2), ent_winlog, "entwin_phys", entvacLog,false,0);
	ent_winlog->SetVisAttributes(new G4VisAttributes(G4Colour(0.7,1.0,0.0)));
    } else {
	// Cryotarget - up against the chamber wall
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad), entLog, "ent_phys", worldlog, false,0);
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad), entvacLog, "entvac_phys", worldlog,false,0);
    }

    // Aluminum
    /*
    int nsec = 24;
    //  Definition taken from HAPLOG 2722 by Juliette, but offset by 31.54 cm
    G4double exit_z[]   = {206*cm, 234.01*cm, 234.02*cm, 253.02*cm, 253.03*cm, 268.26*cm, 268.27*cm,305.29*cm, 305.30*cm,328.71*cm, 328.72*cm, 356.33*cm,356.34*cm, 378.7*cm,378.71*cm, 473.16*cm,473.17*cm, 503.64*cm,503.65*cm, 609.84*cm,609.85*cm, 1161.02*cm, 1161.03*cm,2725.66*cm };
    G4double exit_zero[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    G4double exit_rin[] = {4.128*cm, 4.128*cm, 4.445*cm, 4.445*cm,4.763*cm, 4.763*cm, 5.08*cm,5.08*cm, 6.35*cm, 6.35*cm, 7.62*cm, 7.62*cm,10.16*cm, 10.16*cm,10.478*cm, 10.478*cm,12.7*cm, 12.7*cm,15.24*cm,15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };
    G4double exit_rou[] = {4.432*cm, 4.432*cm, 4.75*cm, 4.75*cm,5.067*cm, 5.067*cm, 5.385*cm,5.385*cm, 6.655*cm, 6.655*cm, 7.925*cm, 7.925*cm, 10.478*cm,10.478*cm,  10.795*cm, 10.795*cm, 13.018*cm, 13.018*cm,15.558*cm, 15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };
    */

    int nsec = 7;
    //  Definition taken from GEN_10M.opc by Bogdan to z = 5.92.  2mm thickness assumed
    G4double exit_z[]   = { 162.2*cm, 592.2*cm, 609.84*cm,609.85*cm, 1161.02*cm, 1161.03*cm,2725.66*cm };
    G4double exit_zero[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    G4double exit_rin[] = { 4.8*cm, 14.8*cm,15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };
    G4double exit_rou[] = { 5.0*cm, 15.0*cm,15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };


    G4Polycone *ext_cone = new G4Polycone("ext_tube", 0.0*deg, 360.0*deg, nsec, exit_z, exit_rin, exit_rou);
    G4Polycone *ext_vac  = new G4Polycone("ext_vac ", 0.0*deg, 360.0*deg, nsec, exit_z, exit_zero, exit_rin);

    G4LogicalVolume *extLog = new G4LogicalVolume(ext_cone, GetMaterial("Aluminum"), "ext_log", 0, 0, 0);
    G4LogicalVolume *extvacLog = new G4LogicalVolume(ext_vac, GetMaterial("Vacuum"), "extvac_log", 0, 0, 0);

    new G4PVPlacement(0,G4ThreeVector(), extLog, "ext_phys", worldlog, false,0);
    new G4PVPlacement(0,G4ThreeVector(), extvacLog, "extvac_phys", worldlog,false,0);

    // Seal this up if we have a gas target
    if( fDetCon->fTargType == kH2 || fDetCon->fTargType == k3He || fDetCon->fTargType == kNeutTarg ){
	// Add in exit Al window

	double extwin_thick = 5.0e-4*cm;

	G4Tubs *extwin = new G4Tubs("ext_win", 0.0, exit_rin[0], extwin_thick/2, 0.*deg, 360.*deg );
	G4LogicalVolume *ext_winlog = new G4LogicalVolume(extwin, GetMaterial("Aluminum"), "entwin_log", 0, 0, 0);
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, exit_z[0] - extwin_thick/2), ext_winlog, "extwin_phys", worldlog,false,0);

	ext_winlog->SetVisAttributes(new G4VisAttributes(G4Colour(0.6,0.6,0.6)));
    }


    double floorthick = 1.0*m;
    G4Tubs *floor_tube = new G4Tubs("floor_tube", 0.0, 30*m, floorthick/2, 0.*deg, 360.*deg );

    G4RotationMatrix *floorrm = new G4RotationMatrix;
    floorrm->rotateX(90*deg);

    G4LogicalVolume *floorLog = new G4LogicalVolume(floor_tube, GetMaterial("Concrete"), "floor_log", 0, 0, 0);
    new G4PVPlacement(floorrm, G4ThreeVector(0.0, -floorthick/2 - beamheight, 0.0), floorLog, "floor_phys", worldlog, false, 0);


    extvacLog->SetVisAttributes(G4VisAttributes::Invisible);
    entvacLog->SetVisAttributes(G4VisAttributes::Invisible);

    G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));

    extLog->SetVisAttributes(pipeVisAtt);
    entLog->SetVisAttributes(pipeVisAtt);


    /*    G4VisAttributes *floorVisAtt= new G4VisAttributes(G4Colour(0.9,0.9,0.9));
	  floorLog->SetVisAttributes(floorVisAtt); */
    floorLog->SetVisAttributes(G4VisAttributes::Invisible);

    return;


    return;

}
