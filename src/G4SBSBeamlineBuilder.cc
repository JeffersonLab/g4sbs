#include "G4SBSBeamlineBuilder.hh"

#include "G4SBSHArmBuilder.hh"
#include "G4SBSTargetBuilder.hh"
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
  double swallrad_inner = 1.041/2.0*m; 
  double beamheight = 10.0*12*2.54*cm; // 10 feet off the ground

  // Stainless
  G4double ent_len = 10*m;
  G4double ent_rin = 31.75*mm;
  G4double ent_rou = ent_rin+0.120*mm;

  G4Tubs *ent_tube = new G4Tubs("ent_tube", ent_rin, ent_rou, ent_len/2, 0.*deg, 360.*deg );
  G4Tubs *ent_vac  = new G4Tubs("ent_vac", 0.0, ent_rin, ent_len/2, 0.*deg, 360.*deg );

  //We want to subtract this cylinder from the entry tube/pipe:
  G4Tubs *cut_cylinder = new G4Tubs("cut_cylinder", 0.0, swallrad, 1.0*m, 0.0*deg, 360.0*deg );

  G4RotationMatrix *cut_cylinder_rot = new G4RotationMatrix;
  cut_cylinder_rot->rotateX( -90.0*deg );

  G4SubtractionSolid *ent_tube_cut = new G4SubtractionSolid( "ent_tube_cut", ent_tube, cut_cylinder, cut_cylinder_rot, 
							     G4ThreeVector( 0.0, 0.0, ent_len/2.0 + swallrad_inner ) );
  G4SubtractionSolid *ent_vac_cut = new G4SubtractionSolid( "ent_vac_cut", ent_vac, cut_cylinder, cut_cylinder_rot, 
							    G4ThreeVector( 0.0, 0.0, ent_len/2.0 + swallrad_inner ) );

  G4LogicalVolume *entLog = new G4LogicalVolume(ent_tube, GetMaterial("Stainless"), "ent_log", 0, 0, 0);
  G4LogicalVolume *entvacLog = new G4LogicalVolume(ent_vac, GetMaterial("Vacuum"), "entvac_log", 0, 0, 0);

  G4LogicalVolume *entLog_cut = new G4LogicalVolume(ent_tube_cut, GetMaterial("Stainless"), "ent_log_cut", 0, 0, 0);
  G4LogicalVolume *entvacLog_cut = new G4LogicalVolume(ent_vac_cut, GetMaterial("Vacuum"), "entvac_log_cut", 0, 0, 0);

  if( targtype == kH2 || targtype == k3He || targtype == kNeutTarg ){
    //if( fDetCon->fTargetBuilder->GetSchamFlag() != 1 ){
    // gas target -  1.5m in air
    new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-1.5*m), entLog, "ent_phys", worldlog, false,0);
    new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-1.5*m), entvacLog, "entvac_phys", worldlog,false,0);
    
    // Add in Be window if no scattering chamber is to be defined:
    if( fDetCon->fTargetBuilder->GetSchamFlag() != 1 ){
      G4double winthick = 0.0127*cm;
    
      G4Tubs *ent_win = new G4Tubs("ent_win", 0.0, ent_rin, winthick/2, 0.*deg, 360.*deg );
      G4LogicalVolume *ent_winlog = new G4LogicalVolume(ent_win, GetMaterial("Beryllium"), "entwin_log", 0, 0, 0);
      new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, ent_len/2-winthick/2), ent_winlog, "entwin_phys", entvacLog,false,0);
      ent_winlog->SetVisAttributes(new G4VisAttributes(G4Colour(0.7,1.0,0.0)));
    } // else {
    //   //Don't add window: we want the beam to interact with the target first. Butt up against the outer edge of the scattering chamber:
    //   new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner), entLog_cut, "ent_phys", worldlog, false,0);
    //   new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner), entvacLog_cut, "entvac_phys", worldlog,false,0);
    // }
  } else {
    // Cryotarget - up against the chamber wall
    new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner), entLog_cut, "ent_phys", worldlog, false,0);
    new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner), entvacLog_cut, "entvac_phys", worldlog,false,0);
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

  entvacLog_cut->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));

  extLog->SetVisAttributes(pipeVisAtt);
  entLog->SetVisAttributes(pipeVisAtt);
  entLog_cut->SetVisAttributes(pipeVisAtt);

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

  if( fDetCon->fExpType == kSIDISExp && fDetCon->fLeadOption == 1 ){
    MakeSIDISLead(worldlog);
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

  double leadstart = 290*cm;
  double leadend   = 435*cm;
  double magleadlen = leadend-leadstart;

  G4Box  *leadbox = new G4Box( "leadbox",  25*cm, 15.0*cm, magleadlen/2 );
  G4Polycone *ext_cone = new G4Polycone("hollowing_tube", 0.0*deg, 360.0*deg, nsec, shield_z, shield_rin, shield_rou);

  G4SubtractionSolid *leadinmag = new G4SubtractionSolid("lead_w_hole", leadbox, ext_cone, 0, G4ThreeVector(0.0, 0.0, -magleadlen/2 - leadstart ) );

  G4LogicalVolume *leadinmag_log = new G4LogicalVolume( leadinmag, GetMaterial("Lead"), "leadinmag", 0, 0, 0 );

  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, leadstart + magleadlen/2), leadinmag_log, "leadinmag_phys", worldlog,false,0);

  ///////////  around opening of 48D48 ///////////////////////////////////////////////

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



void G4SBSBeamlineBuilder::MakeSIDISLead( G4LogicalVolume *worldlog ){
  //Let's fill the beam slot of the magnet with lead: 
  //Height = 31 cm
  //Width = magnet width / 2 - 35 cm
  //depth = magnet depth:
  G4double Beamslot_lead_width = (fDetCon->fHArmBuilder->f48D48width)/2.0-35.0*cm;
  G4double Beamslot_lead_height = 31.0*cm;
  G4double Beamslot_lead_depth = fDetCon->fHArmBuilder->f48D48depth;

  G4Box *Beamslot_lead_box = new G4Box("Beamslot_lead_box", Beamslot_lead_width/2.0, Beamslot_lead_height/2.0, Beamslot_lead_depth/2.0 );
  //G4LogicalVolume *Beamslot_lead_log = new G4LogicalVolume( 

  G4double SBSang = fDetCon->fHArmBuilder->f48D48ang;

  G4RotationMatrix *Beamslot_lead_rm = new G4RotationMatrix;
  Beamslot_lead_rm->rotateY( -SBSang );

  G4ThreeVector SBS_zaxis( sin(SBSang), 0.0, cos(SBSang) );
  G4ThreeVector SBS_yaxis( 0.0, 1.0, 0.0 );
  G4ThreeVector SBS_xaxis( cos(SBSang), 0.0, -sin(SBSang) );

  G4double Beamslot_lead_xoffset = 35.0*cm + Beamslot_lead_width/2.0;

  G4ThreeVector Beamslot_lead_position = (fDetCon->fHArmBuilder->f48D48dist + (fDetCon->fHArmBuilder->f48D48depth)/2.0) * SBS_zaxis - Beamslot_lead_xoffset * SBS_xaxis;

  //Define a subtraction cone for the beam slot in the SBS magnet:
  // G4double zbeampipe[2] = {162.2*cm, 592.2*cm};
  // G4double rinbeampipe[2] = {0.0*cm, 0.0*cm};

  G4double beampipe_subtraction_cone_dz = ((592.2 - 162.2)/2.0)*cm;
  G4double beampipe_subtraction_cone_zpos = ((592.2+162.2)/2.0)*cm;

  G4ThreeVector beampipe_subtraction_cone_position( 0.0, 0.0, beampipe_subtraction_cone_zpos );

  G4Cons *beampipe_subtraction_cone = new G4Cons( "beampipe_subtraction_cone", 0.0*cm, 5.1*cm, 0.0*cm, 15.1*cm, beampipe_subtraction_cone_dz, 0.0*deg, 360.0*deg );

  G4RotationMatrix *Beamslot_lead_rm_inv = new G4RotationMatrix;
  Beamslot_lead_rm_inv->rotateY( SBSang );

  G4ThreeVector beampipe_beamslot_relative_position = beampipe_subtraction_cone_position - Beamslot_lead_position;

  G4ThreeVector beampipe_beamslot_relative_position_local( beampipe_beamslot_relative_position.dot( SBS_xaxis ), 
							   beampipe_beamslot_relative_position.dot( SBS_yaxis ), 
							   beampipe_beamslot_relative_position.dot( SBS_zaxis ) );

  //The subtraction that we want to perform is Beamslot lead box - beampipe cone:
  G4SubtractionSolid *Beamslot_lead_with_hole = new G4SubtractionSolid( "Beamslot_lead_with_hole", Beamslot_lead_box, beampipe_subtraction_cone, Beamslot_lead_rm_inv, beampipe_beamslot_relative_position_local );

  G4LogicalVolume *Beamslot_lead_log = new G4LogicalVolume( Beamslot_lead_with_hole, GetMaterial("Lead"), "Beamslot_lead_log" );
  G4PVPlacement *Beamslot_lead_pv = new G4PVPlacement( Beamslot_lead_rm, Beamslot_lead_position, Beamslot_lead_log, "Beamslot_lead_pv", worldlog, 0, false, 0 );

  //Let's assume a thickness of 5 cm (2 inches lead, whose inner dimensions follow those of the downstream beamline:
  // int nz = 2;
  // G4double zlead[2] = {162.2*cm, 592.2*cm};
  // G4double rinlead[2] = {5.1*cm, 15.1*cm};
  // G4double routlead[2] = {10.1*cm, 20.1*cm};

  //  G4double zstart = Beamslot_lead_position.z();
  G4double zstart = 162.2*cm;
  G4double rinstart = 5.1*cm + (zstart-162.2*cm)/(2.0*beampipe_subtraction_cone_dz) * 10.0*cm;
  G4double zend = 592.2*cm;
  G4double rinend = 15.1*cm;

  G4Cons *leadcone = new G4Cons("leadcone", rinstart, rinstart+5.0*cm, rinend, rinend + 5.0*cm, (zend-zstart)/2.0, 0.0*deg, 360.0*deg );

  G4ThreeVector leadcone_global_position(0.0, 0.0, (zstart+zend)/2.0 );

  // G4ThreeVector leadcone_relative_position = leadcone_global_position - Beamslot_lead_position;
  // G4ThreeVector leadcone_relative_position_local( leadcone_relative_position.dot( SBS_xaxis ), 
  // 						  leadcone_relative_position.dot( SBS_yaxis ),
  // 						  leadcone_relative_position.dot( SBS_zaxis ) );
  G4ThreeVector cutbox_relative_position = Beamslot_lead_position - leadcone_global_position;

  G4SubtractionSolid *leadcone_cut = new G4SubtractionSolid( "leadcone_cut", leadcone, Beamslot_lead_box, Beamslot_lead_rm, cutbox_relative_position );
  
  G4LogicalVolume *leadcone_cut_log = new G4LogicalVolume( leadcone_cut, GetMaterial("Lead"), "leadcone_cut_log" );
  G4PVPlacement *leadcone_cut_pv = new G4PVPlacement( 0, leadcone_global_position, leadcone_cut_log, "leadcone_cut_pv", worldlog, 0, false, 0 );

  // G4Polycone *SIDISlead_cone = new G4Polycone( "SIDISlead_cone", 0.0*deg, 360.0*deg, nz, zlead, rinlead, routlead );
  // G4LogicalVolume *SIDISlead_log = new G4LogicalVolume( SIDISlead_cone, GetMaterial("Lead"), "SIDISlead_log", 0, 0, 0 );
  // G4PVPlacement *SIDISlead_pv = new G4PVPlacement( 0, G4ThreeVector(), SIDISlead_log, "SIDISlead_pv", worldlog, false, 0 );
  
  //We also want to put some lead and/or Iron shielding, i.e., a "collimator" in front of the SBS magnet gap:

  double SBScollwidth = 469.9*mm;
  double SBScollheight = 1219.2*mm;
  double SBScolldepth = 10.0*cm;
  
  double coilspace = 214.5*mm + 20.63*mm;
  double SBS_coll_R = fDetCon->fHArmBuilder->f48D48dist - coilspace - SBScolldepth/2.0 - 5.0*cm;

  double SBS_coll_gapwidth = 50.0*cm*SBS_coll_R/(fDetCon->fHArmBuilder->fRICHdist - 0.5*m) + fDetCon->fTargetBuilder->GetTargLen()*sin( SBSang );
  double SBS_coll_gapheight = 200.0*cm*SBS_coll_R/(fDetCon->fHArmBuilder->fRICHdist - 0.5*m);

  G4Box *SBScoll = new G4Box("SBScoll", 1.5*SBScollwidth, 1.5*SBScollheight, SBScolldepth/2.0 );
  G4Box *SBScoll_hole = new G4Box("SBScoll_hole", SBS_coll_gapwidth/2.0, SBS_coll_gapheight/2.0, SBScolldepth/2.0+1.0*cm );

  G4Cons *SBScoll_cutcone = new G4Cons("SBScoll_cutcone", 0.0*cm, rinstart+5.0*cm, 0.0, rinend + 5.0*cm, (zend-zstart)/2.0, 0.0*deg, 360.0*deg );

  

  G4ThreeVector SBScoll_pos( SBS_coll_R*sin(SBSang), 0.0, SBS_coll_R*cos(SBSang) );
  G4ThreeVector cutcone_relative_pos = leadcone_global_position - SBScoll_pos;

  G4ThreeVector cutcone_relative_pos_local( cutcone_relative_pos.dot(SBS_xaxis),
					    cutcone_relative_pos.dot(SBS_yaxis),
					    cutcone_relative_pos.dot(SBS_zaxis) );

  G4SubtractionSolid *SBS_collimator = new G4SubtractionSolid( "SBS_collimator", SBScoll, SBScoll_hole );
  G4SubtractionSolid *SBS_collimator_beamcut = new G4SubtractionSolid("SBS_collimator_beamcut", SBS_collimator, SBScoll_cutcone, Beamslot_lead_rm_inv, 
								      cutcone_relative_pos_local );

  G4LogicalVolume *SBS_collimator_log = new G4LogicalVolume( SBS_collimator_beamcut, GetMaterial("Lead"), "SBS_collimator_log" );
    
  new G4PVPlacement( Beamslot_lead_rm, G4ThreeVector( SBS_coll_R*sin(SBSang), 0.0, SBS_coll_R*cos(SBSang) ), SBS_collimator_log, "SBS_collimator_phys", worldlog, 
		     0, false, 0 );

  G4VisAttributes *leadVisAtt= new G4VisAttributes(G4Colour(0.25,0.25,0.25));
  //SIDISlead_log->SetVisAttributes(leadVisAtt);

  Beamslot_lead_log->SetVisAttributes(leadVisAtt);
  leadcone_cut_log->SetVisAttributes(leadVisAtt);
  SBS_collimator_log->SetVisAttributes(leadVisAtt);
}




