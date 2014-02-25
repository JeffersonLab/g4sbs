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
#include "G4Box.hh"
#include "G4Polycone.hh"

G4SBSBeamlineBuilder::G4SBSBeamlineBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
    assert(dc);
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

    if( fDetCon->fExpType == kGEp && fDetCon->fLeadOption == 1 ){
	MakeGEpLead(worldlog);
    }

    if( fDetCon->fExpType == kNeutronExp && fDetCon->fTargType != kLD2 ){
	MakeGEnClamp(worldlog);
    }

    if( fDetCon->fExpType == kNeutronExp && fDetCon->fTargType != kLD2 && fDetCon->fLeadOption == 1){
	MakeGEnLead(worldlog);
    }

    return;

}


void G4SBSBeamlineBuilder::MakeGEpLead(G4LogicalVolume *worldlog){
    double maxrad = 25*cm;

    // Lead from scattering chamber to exit pipe in TargetBuilder
    
    // Lead in magnet
    int nsec = 4;
    //  Definition taken from GEN_10M.opc by Bogdan to z = 5.92.  2mm thickness assumed
    G4double exit_z[4]   = {130.0*cm, 162.2*cm, 592.2*cm, 609.84*cm};
    G4double exit_rou[4] = {7.0*cm,  7.0*cm, 17.0*cm ,18.00*cm};
    G4double exit_rin[4] = {0.0*cm,  0.0*cm, 0.0*cm, 0.0*cm };


    // 160 -> 310 cm  box in the magnet
    
    double leadstart = 160*cm;
    double leadend   = 310*cm;
    double magleadlen = leadend-leadstart;

    G4Box  *leadbox = new G4Box( "leadbox",  maxrad, 15.0*cm, magleadlen/2 );
    G4Polycone *ext_cone = new G4Polycone("hollowing_tube", 0.0*deg, 360.0*deg, nsec, exit_z, exit_rin, exit_rou);

    G4SubtractionSolid *leadinmag = new G4SubtractionSolid("lead_w_hole", leadbox, ext_cone, 0, G4ThreeVector(0.0, 0.0, -magleadlen/2 - leadstart ) );

    double cbsize = 50*cm;
    G4Box *leadclip = new G4Box("leadclip_beam", cbsize, cbsize, cbsize);
    G4RotationMatrix *cliprm = new G4RotationMatrix();
    double ang48d48 = fDetCon->fHArmBuilder->f48D48ang;
    cliprm->rotateY( -ang48d48 );

    // Cut away side that interferes with magnet
    leadinmag = new G4SubtractionSolid("lead_w_hole_cut", leadinmag, leadclip, cliprm, 
	    G4ThreeVector( 12.0*cm + cbsize, 0.0, -magleadlen/2 ) );

    G4LogicalVolume *leadinmag_log = new G4LogicalVolume( leadinmag, GetMaterial("Lead"), "leadinmag", 0, 0, 0 );

    new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, leadstart + magleadlen/2), leadinmag_log, "leadinmag_phys", worldlog,false,0);

    

    // Lead from magnet on
    // 311 cm -> 592 cm
    leadstart = 311*cm;
    leadend   = 592*cm;
    magleadlen = leadend-leadstart;

    G4Tubs *leadtube= new G4Tubs( "leadtube",  0*cm, maxrad, magleadlen/2, 0.*deg, 360*deg );
    G4SubtractionSolid *leadafter = new G4SubtractionSolid("lead_after", leadtube, ext_cone, 0, G4ThreeVector(0.0, 0.0, -leadstart-magleadlen/2 ) );

    G4LogicalVolume *leadafter_log = new G4LogicalVolume( leadafter, GetMaterial("Lead"), "leadafter_log", 0, 0, 0 );

    new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, leadstart + magleadlen/2), leadafter_log, "leadafter_phys", worldlog,false,0);






    G4VisAttributes *leadVisAtt= new G4VisAttributes(G4Colour(0.15,0.15,0.15));
    leadinmag_log->SetVisAttributes(leadVisAtt);
    leadafter_log->SetVisAttributes(leadVisAtt);

}

void G4SBSBeamlineBuilder::MakeGEnClamp(G4LogicalVolume *worldlog){
    int nsec = 2;
    //  Definition taken from GEN_10M.opc by Bogdan to z = 5.92.  2mm thickness assumed
    G4double shield_z[]   = { 2.5*m, 5.35*m };
    G4double shield_rin[] = { 8.12*cm, 14.32*cm};
    G4double shield_rou[] = { 10.11*cm, 16.33*cm };

    G4Polycone *shield_cone1 = new G4Polycone("shield_cone1", 0.0*deg, 360.0*deg, nsec, shield_z, shield_rin, shield_rou);
    G4LogicalVolume *shield_cone1_log = new G4LogicalVolume( shield_cone1, GetMaterial("Lead"), "shield_cone1_log", 0, 0, 0 );
    new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, 0.0), shield_cone1_log, "shield_cone1_phys", worldlog,false,0);
}


void G4SBSBeamlineBuilder::MakeGEnLead(G4LogicalVolume *worldlog){

    int nsec = 2;
    G4double clamp1_z[]   = { 162.2*cm, 228.0*cm};
    G4double clamp1_rin[] = { 5.0*cm, 10.5*cm};
    G4double clamp1_rou[] = { 25.0*cm, 25.0*cm};

    G4double clamp2_z[]   = { 2.45*m, 2.85*m,  };
    G4double clamp2_rin[] = { 11.00*cm, 12.0*cm };
    G4double clamp2_rou[] = { 25.0*cm, 25.0*cm};

    G4double clamp3_z[]   = { 4.4*m, 5.90*m,  };
    G4double clamp3_rin[] = { 16.0*cm, 17.00*cm };
    G4double clamp3_rou[] = { 25.0*cm, 25.0*cm};

    G4Polycone *clamp_cone1 = new G4Polycone("clamp_cone1", 0.0*deg, 360.0*deg, nsec, clamp1_z, clamp1_rin, clamp1_rou);
    G4LogicalVolume *clamp_cone1_log = new G4LogicalVolume( clamp_cone1, GetMaterial("Lead"), "clamp_cone1_log", 0, 0, 0 );
    new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, 0.0), clamp_cone1_log, "clamp_cone1_phys", worldlog,false,0);

    G4Polycone *clamp_cone2 = new G4Polycone("clamp_cone2", 0.0*deg, 360.0*deg, nsec, clamp2_z, clamp2_rin, clamp2_rou);
    G4LogicalVolume *clamp_cone2_log = new G4LogicalVolume( clamp_cone2, GetMaterial("Lead"), "clamp_cone2_log", 0, 0, 0 );
    new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, 0.0), clamp_cone2_log, "clamp_cone2_phys", worldlog,false,0);
    
    G4Polycone *clamp_cone3 = new G4Polycone("clamp_cone3", 0.0*deg, 360.0*deg, nsec, clamp3_z, clamp3_rin, clamp3_rou);
    G4LogicalVolume *clamp_cone3_log = new G4LogicalVolume( clamp_cone3, GetMaterial("Lead"), "clamp_cone3_log", 0, 0, 0 );
    new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, 0.0), clamp_cone3_log, "clamp_cone3_phys", worldlog,false,0);


    // 290 -> 435 cm  box in the magnet
    G4double shield_z[]   = { 2.855*m, 4.395*m };
    G4double shield_rin[] = { 0.0, 0.0 };
    G4double shield_rou[] = { 10.50*cm, 17.*cm };

    G4Polycone *shield_cone1 = new G4Polycone("shield_cone1", 0.0*deg, 360.0*deg, nsec, shield_z, shield_rin, shield_rou);
    
    double leadstart = 290*cm;
    double leadend   = 435*cm;
    double magleadlen = leadend-leadstart;

    G4Box  *leadbox = new G4Box( "leadbox",  25*cm, 15.0*cm, magleadlen/2 );
    G4Polycone *ext_cone = new G4Polycone("hollowing_tube", 0.0*deg, 360.0*deg, nsec, shield_z, shield_rin, shield_rou);

    G4SubtractionSolid *leadinmag = new G4SubtractionSolid("lead_w_hole", leadbox, ext_cone, 0, G4ThreeVector(0.0, 0.0, -magleadlen/2 - leadstart ) );

    G4LogicalVolume *leadinmag_log = new G4LogicalVolume( leadinmag, GetMaterial("Lead"), "leadinmag", 0, 0, 0 );

    new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, leadstart + magleadlen/2), leadinmag_log, "leadinmag_phys", worldlog,false,0);

    ///////////  around opening of 48D48 ///////////////////////////////////////////////

    double gapwidth = 22*cm;
    double gapheight= 70*cm;

    double shieldblock3_height = (fDetCon->fHArmBuilder->f48D48depth - gapheight)/2;

    G4Box *extblocklead2 = new G4Box("extblocklead2", 10*cm, 65*cm, 10*cm );
    G4Box *shieldblock3 = new G4Box("shieldblock3", 17*cm/2, shieldblock3_height/2, 10*cm  );
    G4RotationMatrix *windowshieldrm = new G4RotationMatrix();
    windowshieldrm->rotateY(-fDetCon->fHArmBuilder->f48D48ang);

    G4LogicalVolume *windowshield_log = new G4LogicalVolume(extblocklead2, GetMaterial("Lead"),"windowshield_log");
    G4LogicalVolume *shieldblock3_log = new G4LogicalVolume(shieldblock3, GetMaterial("Lead"),"shieldblock3_log");

    G4RotationMatrix *iwindowshieldrm = new G4RotationMatrix( windowshieldrm->inverse() );
    new G4PVPlacement(windowshieldrm, (*iwindowshieldrm)*G4ThreeVector(19*cm, 0.0, fDetCon->fHArmBuilder->f48D48dist-15*cm), windowshield_log, "windowshield_phys", worldlog, false, 0);
    new G4PVPlacement(windowshieldrm, (*iwindowshieldrm)*G4ThreeVector(-19*cm, 0.0, fDetCon->fHArmBuilder->f48D48dist-15*cm), windowshield_log, "windowshield_phys2", worldlog, false, 0);

    new G4PVPlacement(windowshieldrm, (*iwindowshieldrm)*G4ThreeVector(0*cm, gapheight/2+shieldblock3_height/2, fDetCon->fHArmBuilder->f48D48dist-15*cm), shieldblock3_log, "windowshield_phys3", worldlog, false, 0);
    new G4PVPlacement(windowshieldrm, (*iwindowshieldrm)*G4ThreeVector(0*cm, -gapheight/2-shieldblock3_height/2, fDetCon->fHArmBuilder->f48D48dist-15*cm), shieldblock3_log, "windowshield_phys4", worldlog, false, 0);

    






    G4VisAttributes *leadVisAtt= new G4VisAttributes(G4Colour(0.15,0.15,0.15));
    clamp_cone1_log->SetVisAttributes(leadVisAtt);
    clamp_cone2_log->SetVisAttributes(leadVisAtt);
    clamp_cone3_log->SetVisAttributes(leadVisAtt);
    leadinmag_log->SetVisAttributes(leadVisAtt);
    windowshield_log->SetVisAttributes(leadVisAtt);
    shieldblock3_log->SetVisAttributes(leadVisAtt);

}







