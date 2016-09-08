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
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4RotationMatrix.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


G4SBSBeamlineBuilder::G4SBSBeamlineBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
  assert(dc);
}

G4SBSBeamlineBuilder::~G4SBSBeamlineBuilder(){;}

void G4SBSBeamlineBuilder::BuildComponent(G4LogicalVolume *worldlog){
  Targ_t targtype = fDetCon->fTargType;

  double beamheight = 10.0*12*2.54*cm; // 10 feet off the ground
  
  if( fDetCon->fExpType == kGEp && (targtype == kLH2 || targtype == kLD2) ){
    
    MakeGEpBeamline(worldlog);
    
  } else {
    
    double swallrad = 1.143*m/2;
    double swallrad_inner = 1.041/2.0*m; 
    
    
    // Stainless
    G4double ent_len = 10*m;
    //ent_len = ent_len+1.1*m;// for background studies;
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
      new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-1.5*m), entLog, "ent_phys", worldlog, false,0);// -- Nominal beamline --
      new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-1.5*m), entvacLog, "entvac_phys", worldlog,false,0);// -- Nominal beamline --
      //new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-40.0*cm), entLog, "ent_phys", worldlog, false,0);// -- Extended beamline for background studies (2016/09/07)
      //new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-40.0*cm), entvacLog, "entvac_phys", worldlog,false,0);// -- Extended beamline for background studies (2016/09/07)
      
      // Add in Be window if no scattering chamber is to be defined:
      if( fDetCon->fTargetBuilder->GetSchamFlag() != 1 ){
	G4double winthick = 0.0127*cm;
	
	G4Tubs *ent_win = new G4Tubs("ent_win", 0.0, ent_rin, winthick/2, 0.*deg, 360.*deg );
	G4LogicalVolume *ent_winlog = new G4LogicalVolume(ent_win, GetMaterial("Beryllium"), "entwin_log", 0, 0, 0);
	
	// my cancel Be window for GEp experiment 09/29/2014 
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, ent_len/2-winthick/2), ent_winlog, "entwin_phys", entvacLog,false,0);	// => uncommented on 2016/09/07 for background studies
	// my cancel Be window for GEp experiment 09/29/2014
	/* Note from  2016/09/07: */ 
	/* Is that normal that this was commented ? My guess would be not. */
      	
	ent_winlog->SetVisAttributes(new G4VisAttributes(G4Colour(0.7,1.0,0.0)));
      } else {
	//Don't add window: we want the beam to interact with the target first. Butt up against the outer edge of the scattering chamber:
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner), entLog_cut, "ent_phys", worldlog, false,0);
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner), entvacLog_cut, "entvac_phys", worldlog,false,0);
      }
    } else {
      // Cryotarget - up against the chamber wall
      
      //*****6/17
      new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner), entLog_cut, "ent_phys", worldlog, false,0);
      new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner), entvacLog_cut, "entvac_phys", worldlog,false,0);
      //*****6/17
      
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
    G4double exit_z[]   = { 162.2*cm, 592.2*cm, 609.84*cm,609.85*cm, 1161.02*cm, 1161.03*cm,2725.66*cm };// -- Nominal beamline --
    //G4double exit_z[]   = { 37.2*cm, 592.2*cm, 609.84*cm,609.85*cm, 1161.02*cm, 1161.03*cm,2725.66*cm };// -- Extended beamline for background studies (2016/09/07)
    //G4double exit_z_vac[] = { 162.2*cm, 592.2*cm, 610.24*cm,610.35*cm, 1161.52*cm, 1161.53*cm,2726.46*cm };
    
    G4double exit_zero[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    //G4double exit_rin[] = { 4.8*cm, 14.8*cm,15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };
    //G4double exit_rou[] = { 5.0*cm, 15.0*cm,15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };

    G4double exit_rin[] = { 6.065*2.54*cm/2., 14.8*cm,15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };// -- Nominal beamline --
    G4double exit_rou[] = { (6.065/2.0+0.28)*2.54*cm, 15.0*cm,15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };// -- Nominal beamline --
    //G4double exit_rin[] = { 5.64*cm, 14.8*cm, 15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };// -- Extended beamline for background studies (2016/09/07)
    //G4double exit_rou[] = { (5.64+0.28*2.54)*cm, 15.0*cm,15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };// -- Extended beamline for background studies (2016/09/07)
    
    
    G4Polycone *ext_cone = new G4Polycone("ext_cone", 0.0*deg, 360.0*deg, nsec, exit_z, exit_rin, exit_rou);
    //G4Polycone *ext_vac  = new G4Polycone("ext_vac ", 0.0*deg, 360.0*deg, nsec, exit_z_vac, exit_zero, exit_rin);
    G4Polycone *ext_vac  = new G4Polycone("ext_vac ", 0.0*deg, 360.0*deg, nsec, exit_z, exit_zero, exit_rin);
    
    G4LogicalVolume *extLog = new G4LogicalVolume(ext_cone, GetMaterial("Aluminum"), "ext_log", 0, 0, 0);
    G4LogicalVolume *extvacLog = new G4LogicalVolume(ext_vac, GetMaterial("Vacuum"), "extvac_log", 0, 0, 0);
    
    //*****6/17
    new G4PVPlacement(0,G4ThreeVector(), extLog, "ext_phys", worldlog, false,0);
    new G4PVPlacement(0,G4ThreeVector(), extvacLog, "extvac_phys", worldlog,false,0);
    //*****6/17
    
    G4VisAttributes *extVisAtt= new G4VisAttributes(G4Colour(0.9,0.1,0.9));
    extLog->SetVisAttributes(extVisAtt);
    
    
    // Seal this up if we have a gas target
    if( fDetCon->fTargType == kH2 || fDetCon->fTargType == k3He || fDetCon->fTargType == kNeutTarg ){
      // Add in exit Al window
      
      double extwin_thick = 5.0e-4*cm;
      
      G4Tubs *extwin = new G4Tubs("ext_win", 0.0, exit_rin[0], extwin_thick/2, 0.*deg, 360.*deg );
      G4LogicalVolume *ext_winlog = new G4LogicalVolume(extwin, GetMaterial("Aluminum"), "extwin_log", 0, 0, 0);
      
      
      ///*//  cancel Al window to vacuum pipe 10/02/14
      new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, exit_z[0] - extwin_thick/2), ext_winlog, "extwin_phys", worldlog,false,0);
      //*/  // //cancel Al window to vacuum pipe 10/02/14
      

      ext_winlog->SetVisAttributes(new G4VisAttributes(G4Colour(0.5,0.2,0.6)));
    }
    extvacLog->SetVisAttributes(G4VisAttributes::Invisible);
    entvacLog->SetVisAttributes(G4VisAttributes::Invisible);
    
    entvacLog_cut->SetVisAttributes(G4VisAttributes::Invisible);
    
    G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
    
    extLog->SetVisAttributes(pipeVisAtt);
    entLog->SetVisAttributes(pipeVisAtt);
    entLog_cut->SetVisAttributes(pipeVisAtt);
  }

  double floorthick = 1.0*m;
  G4Tubs *floor_tube = new G4Tubs("floor_tube", 0.0, 30*m, floorthick/2, 0.*deg, 360.*deg );

  G4RotationMatrix *floorrm = new G4RotationMatrix;
  floorrm->rotateX(90*deg);

  G4LogicalVolume *floorLog = new G4LogicalVolume(floor_tube, GetMaterial("Concrete"), "floor_log", 0, 0, 0);
  new G4PVPlacement(floorrm, G4ThreeVector(0.0, -floorthick/2 - beamheight, 0.0), floorLog, "floor_phys", worldlog, false, 0);

  /*    G4VisAttributes *floorVisAtt= new G4VisAttributes(G4Colour(0.9,0.9,0.9));
	floorLog->SetVisAttributes(floorVisAtt); */
  floorLog->SetVisAttributes(G4VisAttributes::Invisible);

  if( fDetCon->fExpType == kGEp && fDetCon->fLeadOption == 1 ){
    MakeGEpLead(worldlog);
  }

  if( fDetCon->fExpType == kNeutronExp && fDetCon->fTargType != kLD2 ){ //ought to be in Harm builder?
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


// GEp Beamline Construction --- following Sergey's Fortran code
void G4SBSBeamlineBuilder::MakeGEpBeamline(G4LogicalVolume *worldlog) {

  //Define visualization attributes here:
  G4VisAttributes *ironColor= new G4VisAttributes(G4Colour(0.3,0.3,0.3));
  G4VisAttributes *AlColor= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
  G4VisAttributes *Vacuum_visatt = new G4VisAttributes(G4Colour(0.1, 0.5, 0.9 ) );
  //Vacuum_visatt->SetForceSolid(true);
  Vacuum_visatt->SetVisibility(false);
  G4VisAttributes *SteelColor = new G4VisAttributes( G4Colour( 0.75, 0.75, 0.75 ) );
  G4VisAttributes *CopperColor = new G4VisAttributes( G4Colour( 0.7, 0.3, 0.3 ) );
  
  G4double inch = 2.54*cm;
  
  G4double TargetCenter_zoffset = 6.50*inch;
  G4double ScatChamberRadius = 23.8*inch;

  //Need to make an upstream beamline: 
  G4double upstream_beampipe_zstart = -10.0*m;
  G4double upstream_beampipe_zstop = 0.0;
  G4double upstream_beampipe_rin = 31.75*mm;
  G4double upstream_beampipe_rout = upstream_beampipe_rin + 0.12*mm;

  G4Tubs *upstream_beampipe = new G4Tubs("upstream_beampipe", upstream_beampipe_rin, upstream_beampipe_rout, (upstream_beampipe_zstop - upstream_beampipe_zstart)/2.0, 0.*deg, 360.*deg );
  G4Tubs *upstream_beampipe_vac = new G4Tubs("upstream_beampipe_vac", 0.0, upstream_beampipe_rin, (upstream_beampipe_zstop - upstream_beampipe_zstart)/2.0, 0.*deg, 360.*deg );

  G4Tubs *cut_cylinder = new G4Tubs("cut_cylinder", 0.0, ScatChamberRadius, 0.4*m, 0.0*deg, 360.0*deg );

  G4RotationMatrix *rot_temp = new G4RotationMatrix;
  rot_temp->rotateX(90.0*deg);
  
  G4ThreeVector pos_temp(0,0,0.5*(upstream_beampipe_zstop-upstream_beampipe_zstart)-TargetCenter_zoffset);

  G4SubtractionSolid *upstream_beampipe_cut = new G4SubtractionSolid( "upstream_beampipe_cut", upstream_beampipe, cut_cylinder, rot_temp, pos_temp );
  G4SubtractionSolid *upstream_beampipe_vac_cut = new G4SubtractionSolid("upstream_beampipe_vac_cut", upstream_beampipe_vac, cut_cylinder, rot_temp, pos_temp );

  G4LogicalVolume *upstream_beampipe_log = new G4LogicalVolume(upstream_beampipe_cut, GetMaterial("Stainless_Steel"), "upstream_beampipe_log" );
  G4LogicalVolume *upstream_beampipe_vac_log = new G4LogicalVolume(upstream_beampipe_vac_cut, GetMaterial("Vacuum"), "upstream_beampipe_vac_log" );
  
  upstream_beampipe_log->SetVisAttributes( SteelColor );
  upstream_beampipe_vac_log->SetVisAttributes( Vacuum_visatt );

  pos_temp.set( 0, 0, 0.5*(upstream_beampipe_zstop+upstream_beampipe_zstart) );

  new G4PVPlacement( 0, pos_temp, upstream_beampipe_log, "upstream_beampipe_phys", worldlog, false, 0 );
  new G4PVPlacement( 0, pos_temp, upstream_beampipe_vac_log, "upstream_beampipe_vac_phys", worldlog, false, 0 );


  G4double z_formed_bellows = 52.440*inch - TargetCenter_zoffset; //relative to "target center"? or "origin"?
  G4double z_spool_piece = 58.44*inch - TargetCenter_zoffset;
  G4double z_conic_vacline_weldment = 62.8*inch - TargetCenter_zoffset;
  G4double z_outer_magnetic = 71.782*inch - TargetCenter_zoffset;
  G4double z_inner_magnetic = 73.782*inch - TargetCenter_zoffset;
  G4double z_welded_bellows = 201.632*inch - TargetCenter_zoffset;

  G4double X=0.0, Y=0.0, Z=0.0;
  G4ThreeVector zero(0.0, 0.0, 0.0);
  
  // Conic vacuum line weldment upstream flange:

  G4double Rin, Rout, Thick;
  G4double Rin1, Rout1, Rin2, Rout2;
  Rin = 3.517*inch/2.0;
  Rout = 6.75*inch/2.0;
  
  Thick = 0.84*inch;
  Rin1 = Rin - Thick/2.0*tan( 1.5*deg );
  Rout1 = Rout;
  Rin2 = Rin + Thick/2.0*tan( 1.5*deg );
  Rout2 = Rout;
  
  //G4Cons *CVLW_Flange1 = new G4Cons( "CVLW_Flange1", Rin1, Rout1, Rin2, Rout2, Thick/2.0, 0.0, twopi );
  G4Tubs *CVLW_Flange1 = new G4Tubs("CVLW_Flange1", Rin, Rout, Thick/2.0, 0.0, twopi );
  
  // Fill with vacuum
  Rin = 0.0;
  Rout = 3.517*inch/2.0;

  Rout1 = Rout - Thick/2.0 * tan( 1.5*deg );
  Rout2 = Rout + Thick/2.0 * tan( 1.5*deg );
  
  //G4Cons *CVLW_Flange1_vac = new G4Cons("CVLW_Flange1_vac", Rin, Rout1, Rin, Rout2, Thick/2.0, 0.0, twopi );
  G4Tubs *CVLW_Flange1_vac = new G4Tubs("CVLW_Flange1_vac", Rin, Rout, Thick/2.0, 0.0, twopi );
  
  G4LogicalVolume *CVLW_Flange1_log = new G4LogicalVolume( CVLW_Flange1, GetMaterial("Stainless_Steel"), "CVLW_Flange1_log" );
  G4LogicalVolume *CVLW_Flange1_vac_log = new G4LogicalVolume( CVLW_Flange1_vac, GetMaterial("Vacuum"), "CVLW_Flange1_vac_log" );

  CVLW_Flange1_log->SetVisAttributes( SteelColor );
  CVLW_Flange1_vac_log->SetVisAttributes( Vacuum_visatt );
  
  // Then place the vacuum inside the Iron Tube
  Z = z_conic_vacline_weldment + Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange1_vac_log, "CVLW_Flange1_vac_phys", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange1_log, "CVLW_Flange1_phys", worldlog, false, 0 );

  //conic vacuum line weldment:
  //Thick = 3.50*m;
  Thick = 138.83*inch - 1.12*inch - 0.84*inch;
  Rin1 = 3.517*inch/2.0;
  Rout1 = Rin1 + 0.125*inch;
  Rin2 = 10.734*inch/2.0;
  Rout2 = Rin2 + 0.125*inch;
  
  G4Cons *CVLW = new G4Cons( "CVLW", Rin1, Rout1, Rin2, Rout2, Thick/2.0, 0.0, twopi );   
  // Fill with vacuum				

  Rin1 = 0.0;
  Rout1 = 3.517*inch/2.0;
  Rin2 = 0.0;
  Rout2 = 10.734*inch/2.0;
  
  G4Cons *CVLW_vac = new G4Cons( "CVLW_vac", Rin1, Rout1, Rin2, Rout2, Thick/2.0, 0.0, twopi );
  // Convert into logical volumes
  G4LogicalVolume *CVLW_log = new G4LogicalVolume( CVLW, GetMaterial("Stainless_Steel"), "CVLW_log" );
  G4LogicalVolume *CVLW_vac_log = new G4LogicalVolume( CVLW_vac, GetMaterial("Vacuum"), "CVLW_vac_log" );

  CVLW_log->SetVisAttributes( SteelColor );
  CVLW_vac_log->SetVisAttributes( Vacuum_visatt );
  // Then place the vacuum inside the Iron Cone
  //Z = (159.51 + 2.13)*cm + pDz;
  Z = z_conic_vacline_weldment + 0.84*inch + Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_vac_log, "CVLW_vac_phys", worldlog, false, 0);
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_log, "CVLW_phys", worldlog, false, 0 );

  // Flange 2:
  Rin = 10.734/2.0*inch;
  Rout = 14.0/2.0*inch;
  Thick = 1.12*inch;
  Rin1 = Rin - Thick/2.0 * tan(1.5*deg);
  Rout1 = Rout;
  Rin2 = Rin + Thick/2.0 * tan(1.5*deg);
  Rout2 = Rout;
  //G4Cons *CVLW_Flange2 = new G4Cons("CVLW_Flange2", Rin1, Rout1, Rin2, Rout2, Thick/2.0, 0.0, twopi );
  G4Tubs *CVLW_Flange2 = new G4Tubs("CVLW_Flange2", Rin, Rout, Thick/2.0, 0.0, twopi );
  // Fill with vacuum
  Rin = 0.0;
  Rout1 = Rin1;
  Rout2 = Rin2;
  Rout = 10.734*inch/2.0;
  
  //G4Cons *CVLW_Flange2_vac = new G4Cons("CVLW_Flange2_vac", Rin, Rout1, Rin, Rout2, Thick/2.0, 0.0, twopi );
  G4Tubs *CVLW_Flange2_vac = new G4Tubs( "CVLW_Flange2_vac", Rin, Rout, Thick/2.0, 0.0, twopi );
  // Convert into logical volumes
  // G4LogicalVolume *FLN2_log = new G4LogicalVolume( FLN2_tube, GetMaterial("Iron"), "FLN2_log" );
  // G4LogicalVolume *FVL2_log = new G4LogicalVolume( FVL2_tube, GetMaterial("Vacuum"), "FVL2_log");
  G4LogicalVolume *CVLW_Flange2_log = new G4LogicalVolume( CVLW_Flange2, GetMaterial("Stainless_Steel"), "CVLW_Flange2_log" );
  G4LogicalVolume *CVLW_Flange2_vac_log = new G4LogicalVolume( CVLW_Flange2_vac, GetMaterial("Vacuum"), "CVLW_Flange2_vac_log" );

  CVLW_Flange2_log->SetVisAttributes(SteelColor );
  CVLW_Flange2_vac_log->SetVisAttributes( Vacuum_visatt );

  // Then place the vacuum inside the Iron Tube
  Z = z_welded_bellows - Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange2_vac_log, "CVLW_Flange2_vac_phys", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange2_log, "CVLW_Flange2_phys", worldlog, false, 0 );

  //Next: "Welded bellows"
  G4double dz_welded_bellows = 207.144*inch - z_welded_bellows - TargetCenter_zoffset; // = =5.512 inches
  
  Rin = 11.750/2.0*inch;
  Rout = 14.0/2.0*inch;
  Thick = 1.12*inch;
  
  G4Tubs *WB_Flange = new G4Tubs( "WB_Flange", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *WB_Flange_log = new G4LogicalVolume( WB_Flange, GetMaterial("Stainless_Steel"), "WB_Flange_log" );

  WB_Flange_log->SetVisAttributes( SteelColor );
  
  Z = z_welded_bellows + Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Flange_log, "WB_Flange1_phys", worldlog, false, 0 );

  Z = z_welded_bellows + dz_welded_bellows - Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Flange_log, "WB_Flange2_phys", worldlog, false, 1 );
  
  Rout = Rin + 0.125*inch;
  Thick = dz_welded_bellows - 2*1.12*inch;
  G4Tubs *WB_Bellows = new G4Tubs( "WB_Bellows", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *WB_Bellows_log = new G4LogicalVolume(WB_Bellows, GetMaterial("Stainless_Steel"), "WB_Bellows_log" );

  WB_Bellows_log->SetVisAttributes( SteelColor );
  
  Z = z_welded_bellows + 1.12*inch + Thick/2.0;
  
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Bellows_log, "WB_Bellows_phys", worldlog, false, 0 );
  
  Rin = 0.0;
  Rout = 11.750/2.0*inch;
  Thick = dz_welded_bellows;
  G4Tubs *WB_Vacuum = new G4Tubs( "WB_Vacuum", Rin, Rout, Thick/2.0, 0.0, twopi );

  G4LogicalVolume *WB_Vacuum_log = new G4LogicalVolume(WB_Vacuum, GetMaterial("Vacuum"), "WB_Vacuum_log" );

  WB_Vacuum_log->SetVisAttributes( Vacuum_visatt );
  
  Z = z_welded_bellows + dz_welded_bellows/2.0;

  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Vacuum_log, "WB_Vacuum_phys", worldlog, false, 0 );
  
  // // Here a bellow and we assign wall of 0.03 cm
  // tRmin = (0.5*27.62)*cm;
  // tRmax = (0.5*27.62 + 0.03)*cm;
  // tDzz  = 0.5*(4.237*2.54 - 2.84)*cm;
  // G4Tubs *TBL8_tube = new G4Tubs("TBL8_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  // // Fill with vacuum
  // tRmin = 0.0*cm;
  // tRmax = (0.5*27.62)*cm;
  // G4Tubs *TVL8_tube = new G4Tubs("TVL8_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  // // Convert into logical volumes
  // G4LogicalVolume *TBL8_log = new G4LogicalVolume( TBL8_tube, GetMaterial("Iron"), "TBL8_log" );
  // G4LogicalVolume *TVL8_log = new G4LogicalVolume( TVL8_tube, GetMaterial("Vacuum"), "TVL8_log" );
  // // Then place the vacuum inside the Iron Tube
  // Z = (159.51 + 352.636 - 2.84*0.5 + 4.237*2.54*0.5)*cm;
  // new G4PVPlacement( 0, zero, TVL8_log, "Bellow_Vac", TBL8_log, false, 0 );
  // new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TBL8_log, "Bellow_Iron", worldlog, false, 0 );

  //return;
  
  // EXTEND VACUUM LINE by using Maduka geometry
  // ============================================
  G4double tRmin, tRmax, tDzz, pDz, pRmax1, pRmax2, tSPhi, tDphi, pRmin1, pRmin2;
  tSPhi = 0.0;
  tDphi = twopi;
  
  tRmin = 0.5*12.0*2.54*cm; 
  tRmax = 13.0*2.54*0.5*cm;
  tDzz  = 0.5*41.0*2.54*cm;
  G4Tubs *TBL9_tube = new G4Tubs("TBL9_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  // Fill with vacuum
  tRmin = 0.0*cm;
  tRmax = 0.5*12.0*2.54*cm;
  G4Tubs *TVL9_tube = new G4Tubs("TVL9_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  // Convert into logical volumes
  G4LogicalVolume *TBL9_log = new G4LogicalVolume( TBL9_tube, GetMaterial("Aluminum"), "TBL9_log" );
  G4LogicalVolume *TVL9_log = new G4LogicalVolume( TVL9_tube, GetMaterial("Vacuum"), "TVL9_log" );

  TBL9_log->SetVisAttributes( AlColor );
  TVL9_log->SetVisAttributes( Vacuum_visatt );
  
  // Then place the vacuum inside the Al Tube
  //Z = (159.51 + 352.636 - 2.84*0.5 + 4.237*2.54 + 41.0*2.54*0.5)*cm;
  Z = 207.144*inch + tDzz - TargetCenter_zoffset;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TVL9_log, "Extended_Vac1", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TBL9_log, "Extended_Al1", worldlog, false, 0 );

  tRmin = 0.5*24.0*2.54*cm;
  tRmax = 25.0*2.54*0.5*cm;
  tDzz  = 0.5*217.0*2.54*cm;
  G4Tubs *TML9_tube = new G4Tubs( "TML9_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  // Fill with vacuum
  tRmin = 0.0*cm;
  tRmax = 0.5*24.0*2.54*cm;
  G4Tubs *TMV9_tube = new G4Tubs("TMV9_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  // Convert into logical volumes
  G4LogicalVolume *TML9_log = new G4LogicalVolume( TML9_tube, GetMaterial("Aluminum"), "TML9_log" );
  G4LogicalVolume *TMV9_log = new G4LogicalVolume( TMV9_tube, GetMaterial("Vacuum"), "TMV9_log" );

  TML9_log->SetVisAttributes( AlColor );
  TMV9_log->SetVisAttributes( Vacuum_visatt );
  // Then place vacuum inside of Al tube
  //Z = (159.51 + 352.636 - 2.84*0.5 + 4.237*2.54 + 41.0*2.54 + 0.5*217.0*2.54)*cm;
  Z = 207.144*inch + 41.0*inch + tDzz - TargetCenter_zoffset;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TMV9_log, "Extended_Vac2", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TML9_log, "Extended_Al2", worldlog, false, 0 );

  // For CPU speed, extend vacuum all the way to the edge of the "world" volume, so that we don't track beam electrons in air beyond interesting region.
  G4double Zstop = 50.0*m;
  G4double Zstart = Z + tDzz;
  G4double Zwidth = (Zstop-Zstart);
  G4Tubs *FakeVacuumExtension = new G4Tubs( "FakeVacuumExtension", tRmin, tRmax, Zwidth/2.0, tSPhi, tDphi );
  G4LogicalVolume *FakeVacuumExtension_log = new G4LogicalVolume( FakeVacuumExtension, GetMaterial("Vacuum"), "FakeVacuumExtension_log" );
  FakeVacuumExtension_log->SetVisAttributes( Vacuum_visatt );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0.5*(Zstop+Zstart)), FakeVacuumExtension_log, "FakeVacuumExtension_phys", worldlog,false,0);

  //-----------------------------------------------------
  //       magnetic tubes

  //Inner Magnetic 1:
  Rin1 = 4.354*inch/2.0;
  Rout1 = Rin1 + 0.25*inch;
  Rin2 = 6.848*inch/2.0;
  Rout2 = Rin2 + 0.25*inch;
  Thick = 47.625*inch;

  G4Cons *IM1 = new G4Cons( "IM1", Rin1, Rout1, Rin2, Rout2, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *IM1_log = new G4LogicalVolume( IM1, GetMaterial("Iron"), "IM1_log" );

  IM1_log->SetVisAttributes( ironColor );
  
  Z = z_inner_magnetic + Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), IM1_log, "IM1_phys", worldlog, false, 0 );
  
  Rin1 = 8.230*inch/2.0;
  Rout1 = Rin1 + 0.25*inch;
  Rin2 = 10.971*inch/2.0;
  Rout2 = Rin2 + 0.25*inch;
  Thick = 52.347*inch;
  
  G4Cons *IM2 = new G4Cons( "IM2", Rin1, Rout1, Rin2, Rout2, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *IM2_log = new G4LogicalVolume( IM2, GetMaterial("Iron"), "IM2_log" );

  IM2_log->SetVisAttributes( ironColor );
  
  Z = z_inner_magnetic + 74.00*inch + Thick/2.0;

  new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), IM2_log, "IM2_phys", worldlog, false, 0 );
  
  G4double OMthick = 1.625*inch;
  G4double OMspace = 0.375*inch;

  G4double zmin = z_outer_magnetic;
  G4double zmax = zmin + 26.0*OMthick + 25.0*OMspace;
  
  G4double  Rin_min = 5.5*inch/2.0;
  G4double  Rin_max = 8.178*inch/2.0;
  for( G4int i=0; i<26; i++ ){
    char cname[100];
    sprintf(cname,"OM1_ring%d", i);
    G4String name = cname;
    
    G4double zstart = zmin + i*(OMthick + OMspace);
    G4double zstop = zstart + OMthick;
    G4double Rin_start = Rin_min + (zstart-zmin)/(zmax-zmin)*(Rin_max - Rin_min);
    G4double Rout_start = Rin_start + 0.5*inch;
    G4double Rin_stop = Rin_min + (zstop-zmin)/(zmax-zmin)*(Rin_max - Rin_min);
    G4double Rout_stop = Rin_stop + 0.5*inch;

    G4Cons *ring = new G4Cons( name,Rin_start, Rout_start, Rin_stop, Rout_stop, OMthick/2.0, 0.0, twopi );

    name += "_log";
    G4LogicalVolume *ring_log = new G4LogicalVolume( ring, GetMaterial("Iron"), name );

    ring_log->SetVisAttributes( ironColor );
    
    name = cname;
    name += "_phys";
    
    Z = 0.5*(zstart + zstop);
    new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), ring_log, name, worldlog, false, 0 );
    
  }

  zmin = z_outer_magnetic + 74.0*inch;
  zmax = zmin + 27.0*OMthick + 26.0*OMspace;

  Rin_min = 9.349*inch/2.0;
  Rin_max = 12.156*inch/2.0;
  for(G4int i=0; i<27; i++){
    char cname[100];
    sprintf(cname,"OM2_ring%d", i );
    G4String name = cname;

    G4double zstart = zmin + i*(OMthick+OMspace);
    G4double zstop = zstart + OMthick;
    G4double Rin_start = Rin_min + (zstart-zmin)/(zmax-zmin)*(Rin_max-Rin_min);
    G4double Rout_start = Rin_start + 0.5*inch;
    G4double Rin_stop = Rin_min + (zstop-zmin)/(zmax-zmin)*(Rin_max-Rin_min);
    G4double Rout_stop = Rin_stop + 0.5*inch;

    G4Cons *ring = new G4Cons("name", Rin_start, Rout_start, Rin_stop, Rout_stop, OMthick/2.0, 0.0, twopi );
    name += "_log";
    G4LogicalVolume *ring_log = new G4LogicalVolume( ring, GetMaterial("Iron"), name );

    ring_log->SetVisAttributes( ironColor );
    name = cname;
    name += "_phys";
    Z = 0.5*(zstart + zstop);
    new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), ring_log, name, worldlog, false, 0 );
  }
  
  G4double dz_spool_piece = z_conic_vacline_weldment - z_spool_piece;

  //Make Spool piece vacuum:
  Rout = 3.76*inch/2.0;
  Rin = 0.0;
  Thick = dz_spool_piece;
  
  G4Tubs *SpoolPiece_vac = new G4Tubs("SpoolPiece_vac", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *SpoolPiece_vac_log = new G4LogicalVolume( SpoolPiece_vac, GetMaterial("Vacuum"), "SpoolPiece_vac_log" );

  SpoolPiece_vac_log->SetVisAttributes( Vacuum_visatt );
  
  Z = z_spool_piece + Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), SpoolPiece_vac_log, "SpoolPiece_vac_phys", worldlog, false, 0 );

  Rin = 3.76*inch/2.0;
  Rout = 6.00*inch/2.0;
  Thick = 0.84*inch;
  
  G4Tubs *SpoolPiece_Flange1 = new G4Tubs("SpoolPiece_Flange1", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *SpoolPiece_Flange1_log = new G4LogicalVolume( SpoolPiece_Flange1, GetMaterial("Stainless_Steel"), "SpoolPiece_Flange1_log" );

  SpoolPiece_Flange1_log->SetVisAttributes( SteelColor );
  
  Z = z_spool_piece + Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), SpoolPiece_Flange1_log, "SpoolPiece_Flange1_phys", worldlog, false, 0 );

  Rout = 6.75*inch/2.0;
  Thick = 0.84*inch;
  
  G4Tubs *SpoolPiece_Flange2 = new G4Tubs("SpoolPiece_Flange2", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *SpoolPiece_Flange2_log = new G4LogicalVolume( SpoolPiece_Flange2, GetMaterial("Stainless_Steel"), "SpoolPiece_Flange2_log" );

  SpoolPiece_Flange2_log->SetVisAttributes( SteelColor );

  Z = z_conic_vacline_weldment - Thick/2.0;

  new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), SpoolPiece_Flange2_log, "SpoolPiece_Flange2_phys", worldlog, false, 0 );

  Rout = 4.0*inch/2.0;
  Thick = dz_spool_piece - 2.0*0.84*inch;

  G4Tubs *SpoolPiece_tube = new G4Tubs("SpoolPiece_tube", Rin, Rout, Thick/2.0, 0.0, twopi );

  G4LogicalVolume *SpoolPiece_tube_log = new G4LogicalVolume( SpoolPiece_tube, GetMaterial("Stainless_Steel"), "SpoolPiece_tube_log" );

  SpoolPiece_tube_log->SetVisAttributes( SteelColor );
  
  Z = z_spool_piece + dz_spool_piece/2.0;

  new G4PVPlacement( 0,  G4ThreeVector( X, Y, Z ), SpoolPiece_tube_log, "SpoolPiece_tube_phys", worldlog, false, 0 );
 
  //Last but not least: formed bellows! defer to tomorrow...

  G4double dz_formed_bellows = 6.00*inch;
  Rin = 0.0;
  Rout = 3.81*inch/2.0;
  Thick = dz_formed_bellows;
  //define vacuum volume for formed bellows
  G4Tubs *FormedBellows_vac = new G4Tubs("FormedBellows_vac", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *FormedBellows_vac_log = new G4LogicalVolume( FormedBellows_vac, GetMaterial("Vacuum"), "FormedBellows_vac_log" );

  FormedBellows_vac_log->SetVisAttributes( Vacuum_visatt );
  
  Z = z_formed_bellows + Thick/2.0;

  new G4PVPlacement(  0,  G4ThreeVector( X, Y, Z ), FormedBellows_vac_log, "FormedBellows_vac_phys", worldlog, false, 0 );

  Rin = 3.81*inch/2.0;
  Rout = 6.00*inch/2.0;
  Thick = 0.84*inch;

  //Flanges for formed bellows:
  G4Tubs *FormedBellows_Flange = new G4Tubs("FormedBellows_Flange", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *FormedBellows_Flange_log = new G4LogicalVolume( FormedBellows_Flange, GetMaterial("Stainless_Steel"), "FormedBellows_Flange_log" );

  FormedBellows_Flange_log->SetVisAttributes( SteelColor );
  
  Z = z_formed_bellows + Thick/2.0;

  new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), FormedBellows_Flange_log, "FormedBellows_Flange1_phys", worldlog, false, 0 );

  Z = z_formed_bellows + dz_formed_bellows - Thick/2.0;

  new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), FormedBellows_Flange_log, "FormedBellows_Flange2_phys", worldlog, false, 1 );

  //Tube for formed bellows:

  Rout = Rin + 0.125*inch; //This is just a guess!!
  Thick = dz_formed_bellows - 2.0*0.84*inch;

  G4Tubs *FormedBellows_tube = new G4Tubs( "FormedBellows_tube", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *FormedBellows_tube_log = new G4LogicalVolume( FormedBellows_tube, GetMaterial("Stainless_Steel"), "FormedBellows_tube_log" );

  FormedBellows_tube_log->SetVisAttributes( SteelColor );
  
  Z = z_formed_bellows + dz_formed_bellows/2.0;

  new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), FormedBellows_tube_log, "FormedBellows_tube_phys", worldlog, false, 0 );

  //Two more "Iron" tubes to connect Snout to "formed bellows"
  G4double dz_iron_tubes = z_formed_bellows - 49.56*inch + TargetCenter_zoffset;

  Thick = dz_iron_tubes/2.0;
  Rin = 5.0*cm;
  Rout = 7.0*cm;

  G4Tubs *IronTube1 = new G4Tubs("IronTube1", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *IronTube1_log = new G4LogicalVolume( IronTube1, GetMaterial("Iron"), "IronTube1_log" );
  IronTube1_log->SetVisAttributes( ironColor );
  
  G4Tubs *IronTube1_vac = new G4Tubs("IronTube1_vac", 0.0, Rin, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *IronTube1_vac_log = new G4LogicalVolume( IronTube1_vac, GetMaterial("Vacuum"), "IronTube1_vac_log" );

  IronTube1_vac_log->SetVisAttributes( Vacuum_visatt );

  Z = 49.56*inch + Thick/2.0 - TargetCenter_zoffset;

  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), IronTube1_log, "IronTube1_phys", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), IronTube1_vac_log, "IronTube1_vac_phys", worldlog, false, 0 );

  Rin = 2.415*inch;
  Rout = 2.5*inch;
  
  G4Tubs *IronTube2 = new G4Tubs("IronTube2", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4Tubs *IronTube2_vac = new G4Tubs("IronTube2_vac", 0.0, Rin, Thick/2.0, 0.0, twopi );

  G4LogicalVolume *IronTube2_log = new G4LogicalVolume( IronTube2, GetMaterial("Iron"), "IronTube2_log" );
  G4LogicalVolume *IronTube2_vac_log = new G4LogicalVolume( IronTube2_vac, GetMaterial("Vacuum"), "IronTube2_vac_log" );

  IronTube2_log->SetVisAttributes(ironColor);
  IronTube2_vac_log->SetVisAttributes(Vacuum_visatt);

  Z += Thick;

  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), IronTube2_log, "IronTube2_phys", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), IronTube2_vac_log, "IronTube2_vac_phys", worldlog, false, 0 );


  //Next, corrector magnets:
  G4double UpstreamCoilThickY = 1.68*inch;
  G4double UpstreamCoilThickX = 3.46*inch;
  //G4double UpstreamCoilWidth = 3.46*inch;
  G4double UpstreamCoilHeight = 8.17*inch;
  G4double UpstreamCoilDepth = 6.60*inch;
  G4double UpstreamCoilWidth = 7.56*inch;

  G4Box *UpstreamCoil_outer = new G4Box("UpstreamCoil_outer", UpstreamCoilThickX/2.0, (UpstreamCoilHeight+2.0*UpstreamCoilThickY)/2.0, (UpstreamCoilDepth + 2.0*UpstreamCoilThickY)/2.0 );
  G4Box *UpstreamCoil_inner = new G4Box("UpstreamCoil_inner", UpstreamCoilThickX/2.0 + cm, UpstreamCoilHeight/2.0, UpstreamCoilDepth/2.0 );

  G4SubtractionSolid *UpstreamCoil = new G4SubtractionSolid( "UpstreamCoil", UpstreamCoil_outer, UpstreamCoil_inner, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *UpstreamCoil_log = new G4LogicalVolume(UpstreamCoil, GetMaterial("Copper"), "UpstreamCoil_log" );

  UpstreamCoil_log->SetVisAttributes( CopperColor );

  Z = z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY;
  X = (UpstreamCoilWidth+UpstreamCoilThickX)/2.0;
  Y = 0.0;

  //two placements of upstream coil:
  
  new G4PVPlacement( 0, G4ThreeVector(-X, Y, Z), UpstreamCoil_log, "UpstreamCoil_phys_right", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), UpstreamCoil_log, "UpstreamCoil_phys_left", worldlog, false, 1 );

  G4double UpstreamPoleDepth = 6.3*inch;
  G4double UpstreamPoleWidth = 4.02*inch;
  G4double UpstreamPoleHeight = 7.87*inch;
  //Next, make poles:
  G4Box *UpstreamPole = new G4Box( "UpstreamPole", UpstreamPoleWidth/2.0, UpstreamPoleHeight/2.0, UpstreamPoleDepth/2.0 );
  G4LogicalVolume *UpstreamPole_log = new G4LogicalVolume( UpstreamPole, GetMaterial("Iron"), "UpstreamPole_log" );
  UpstreamPole_log->SetVisAttributes( ironColor );
  //two placements of upstream poles:

  new G4PVPlacement( 0, G4ThreeVector(-X, Y, Z), UpstreamPole_log, "UpstreamPole_phys_right", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), UpstreamPole_log, "UpstreamPole_phys_left", worldlog, false, 1 );

  //Next, make surrounding yoke:
  G4double YokeTopPiece_Width = 15.04*inch;
  G4double YokeTopPiece_Height = 3.94*inch;
  G4double YokeTopPiece_Depth = 6.30*inch;

  G4Box *YokeTopPiece = new G4Box("YokeTopPiece", YokeTopPiece_Width/2.0, YokeTopPiece_Height/2.0, YokeTopPiece_Depth/2.0 );
  G4LogicalVolume *YokeTopPiece_log = new G4LogicalVolume( YokeTopPiece, GetMaterial("Iron"), "YokeTopPiece_log" );

  YokeTopPiece_log->SetVisAttributes( ironColor );
  
  X = 0.0;
  Y = (11.81*inch + YokeTopPiece_Height)/2.0;

  //two placements of yoke top piece (top and bottom symmetric):
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeTopPiece_log, "UpstreamYokeTop_phys", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(X,-Y,Z), YokeTopPiece_log, "UpstreamYokeBottom_phys", worldlog, false, 1 );

  G4double YokeLeftPiece_Width = 2.76*inch;
  G4double YokeLeftPiece_Height = 11.81*inch + 2.0*YokeTopPiece_Height;
  G4double YokeLeftPiece_Depth = YokeTopPiece_Depth;

  G4Box *YokeLeftPiece = new G4Box("YokeLeftPiece", YokeLeftPiece_Width/2.0, YokeLeftPiece_Height/2.0, YokeLeftPiece_Depth/2.0 );
  G4LogicalVolume *YokeLeftPiece_log = new G4LogicalVolume( YokeLeftPiece, GetMaterial("Iron"), "YokeLeftPiece_log" );
  YokeLeftPiece_log->SetVisAttributes(ironColor );
  
  X = 7.52*inch + YokeLeftPiece_Width/2.0;
  Y = 0.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeLeftPiece_log, "UpstreamYokeLeftPiece_phys", worldlog, false, 0 );

  G4double YokeRightNotchAngle = 18.43*deg;
  G4double YokeRightWidthFinal = YokeLeftPiece_Width;
  G4double YokeRightZFinal = YokeLeftPiece_Depth - 0.39*inch;
  G4double YokeRightWidthInitial = YokeRightWidthFinal - YokeRightZFinal*tan(YokeRightNotchAngle );

  //I *think* this is correct:
  G4Trap *YokeRight_trap = new G4Trap( "YokeRight_trap", YokeRightZFinal/2.0, atan( (YokeRightWidthFinal-YokeRightWidthInitial)/2.0/YokeRightZFinal ), 180.0*deg,
				       YokeLeftPiece_Height/2.0, YokeRightWidthInitial/2.0, YokeRightWidthInitial/2.0, 0.0,
				       YokeLeftPiece_Height/2.0, YokeRightWidthFinal/2.0, YokeRightWidthFinal/2.0, 0.0 ); 

  G4Box *YokeRight_box = new G4Box( "YokeRight_box", YokeRightWidthFinal/2.0, YokeLeftPiece_Height/2.0, 0.39*inch/2.0 );

  X = 0.5*(YokeRightWidthFinal/2.0 + YokeRightWidthInitial/2.0) - YokeRightWidthFinal/2.0;
  
  G4UnionSolid *YokeRightPiece = new G4UnionSolid("YokeRightPiece", YokeRight_trap, YokeRight_box, 0, G4ThreeVector( X, 0, (YokeRightZFinal+0.39*inch)/2.0 ) );
  G4LogicalVolume *YokeRightPiece_log = new G4LogicalVolume(YokeRightPiece, GetMaterial("Iron"), "YokeRightPiece_log" );

  YokeRightPiece_log->SetVisAttributes(ironColor);
  
  X = -7.52*inch - 0.5*(YokeRightWidthFinal/2.0 + YokeRightWidthInitial/2.0);
  Y = 0.0;
  Z = z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeRightPiece_log, "UpstreamYokeRightPiece_phys", worldlog, false, 0 );

  //Downstream Corrector:
  G4double DownstreamTotalWidth = 17.58*inch + 2.0*2.76*inch;
  G4double DownstreamTotalHeight = 20.16*inch + 2.0*3.94*inch;
  G4double DownstreamYokeDepth = 15.75*inch;

  G4double DownstreamYokeGapWidth = 17.58*inch;
  G4double DownstreamYokeGapHeight = 20.16*inch;
  G4Box *DownstreamYoke_box = new G4Box("DownstreamYoke_box", DownstreamTotalWidth/2.0, DownstreamTotalHeight/2.0, DownstreamYokeDepth/2.0 );
  G4Box *DownstreamYoke_gap = new G4Box("DownstreamYoke_gap", DownstreamYokeGapWidth/2.0, DownstreamYokeGapHeight/2.0, DownstreamYokeDepth/2.0+cm );
  G4SubtractionSolid *DownstreamYoke = new G4SubtractionSolid( "DownstreamYoke", DownstreamYoke_box, DownstreamYoke_gap, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *DownstreamYoke_log = new G4LogicalVolume( DownstreamYoke, GetMaterial("Iron"), "DownstreamYoke_log" );

  DownstreamYoke_log->SetVisAttributes( ironColor );

  X = 0.0; Y = 0.0;
  Z = z_formed_bellows + 76.09*inch + 1.71*inch + DownstreamYokeDepth/2.0;

  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), DownstreamYoke_log, "DownstreamYoke_phys", worldlog, false, 0 );

  G4double DS_coil_depth = 8.91*inch;
  G4double DS_coil_height = 12.04*inch;
  G4double DS_coil_ThickX = 2.90*inch;
  G4double DS_coil_ThickY = 1.68*inch;

  G4Box *DS_coil_outer = new G4Box( "DS_coil_outer", DS_coil_ThickX/2.0, (DS_coil_height + 2.0*DS_coil_ThickY)/2.0, (DS_coil_depth + 2.0*DS_coil_ThickY)/2.0 );
  G4Box *DS_coil_inner = new G4Box( "DS_coil_inner", DS_coil_ThickX/2.0+cm, DS_coil_height/2.0, DS_coil_depth/2.0 );

  G4SubtractionSolid *DS_coil = new G4SubtractionSolid( "DS_coil", DS_coil_outer, DS_coil_inner, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *DS_coil_log = new G4LogicalVolume( DS_coil, GetMaterial("Copper"), "DS_coil_log" );
  DS_coil_log->SetVisAttributes(CopperColor );
  
  X = 11.67*inch/2.0 + DS_coil_ThickX/2.0;
  Y = 0.0;
  Z = z_formed_bellows + 76.09*inch + DS_coil_ThickY + DS_coil_depth/2.0;

  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), DS_coil_log, "DS_coil_phys_left", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(-X,Y,Z), DS_coil_log, "DS_coil_phys_right", worldlog, false, 1 );

  //Now just need poles:
  G4double DSpole_depth = 8.76*inch;
  G4double DSpole_width = (17.58-11.00)*inch/2.0;
  G4double DSpole_height = 11.81*inch;

  G4Box *DSpole = new G4Box("DSpole", DSpole_width/2.0, DSpole_height/2.0, DSpole_depth/2.0 );
  G4LogicalVolume *DSpole_log = new G4LogicalVolume(DSpole, GetMaterial("Iron"), "DSpole_log" );

  DSpole_log->SetVisAttributes(ironColor);
  
  X = (17.58+11.00)*inch/4.0;
  Y = 0.0;
  //two placements of poles:
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), DSpole_log, "DSpole_phys_left", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(-X,Y,Z), DSpole_log, "DSpole_phys_right", worldlog, false, 1 );

  // VISUALS
  
  // CVLW_Flange1_log->SetVisAttributes( ironColor );
  // CVLW_log->SetVisAttributes( ironColor );
  // CVLW_Flange2_log->SetVisAttributes( ironColor );
  // WB_Flange_log->SetVisAttributes( ironColor );
  // WB_Bellows_log->SetVisAttributes( ironColor );
  // //TBL8_log->SetVisAttributes( ironColor );
  // TBM1_log->SetVisAttributes( ironColor );
  // TBM2_log->SetVisAttributes( ironColor );
  // TBM3_log->SetVisAttributes( ironColor );
  // TBM4_log->SetVisAttributes( ironColor );
  // TBT1_log->SetVisAttributes( ironColor );
  // TBT2_log->SetVisAttributes( ironColor );

  
  // TBL9_log->SetVisAttributes( AlColor );
  // TML9_log->SetVisAttributes( AlColor );

  // // Vacuum
  // FVL1_log->SetVisAttributes( G4VisAttributes::Invisible );
  // FVL2_log->SetVisAttributes( G4VisAttributes::Invisible );
  // FVL3_log->SetVisAttributes( G4VisAttributes::Invisible );
  // FVL5_log->SetVisAttributes( G4VisAttributes::Invisible );
  // FVL6_log->SetVisAttributes( G4VisAttributes::Invisible );
  // FVL7_log->SetVisAttributes( G4VisAttributes::Invisible );
  // TVB1_log->SetVisAttributes( G4VisAttributes::Invisible );
  // TVL8_log->SetVisAttributes( G4VisAttributes::Invisible );
  // TVL9_log->SetVisAttributes( G4VisAttributes::Invisible );
  // TMV9_log->SetVisAttributes( G4VisAttributes::Invisible );
  // TTV1_log->SetVisAttributes( G4VisAttributes::Invisible );
  // TTV2_log->SetVisAttributes( G4VisAttributes::Invisible );
}

//  Here is lead shield of beam line for GEp

void G4SBSBeamlineBuilder::MakeGEpLead(G4LogicalVolume *worldlog){

  G4VisAttributes *lead_visatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
  
  G4double inch = 2.54*cm;
  G4double TargetCenter_zoffset = 6.50*inch;

  G4double z_outer_magnetic = 182.33*cm - TargetCenter_zoffset;
  
  G4double zstart_lead1 = 170.0*cm;
  G4double z_formed_bellows = 133.2*cm - TargetCenter_zoffset;
  G4double zstop_lead1 = z_formed_bellows + 75.0*inch;

  //G4cout << "zmag, zstart, zstop = " << z_outer_magnetic << ", " << zstart_lead1 << ", " << zstop_lead1 << G4endl;
  
  G4double Rin1 = 9.3*cm;
  G4double Rout1 = Rin1 + 5.0*cm;
  G4double Rin2 = Rin1 + (zstop_lead1 - zstart_lead1)*tan(1.5*deg );
  G4double Rout2 = Rin2 + 5.0*cm;

  G4Cons *leadcone1 = new G4Cons("leadcone1", Rin1, Rout1, Rin2, Rout2, (zstop_lead1-zstart_lead1)/2.0, 90.0*deg, 180.0*deg );

  G4double width = 0.5*fDetCon->fHArmBuilder->f48D48width;
  G4double angle = fDetCon->fHArmBuilder->f48D48ang;
  G4double dist = fDetCon->fHArmBuilder->f48D48dist;
  G4double depth = fDetCon->fHArmBuilder->f48D48depth;
  
  //We want to subtract the overlap between the cone and the SBS magnet.
  G4Box *cutbox1 = new G4Box( "cutbox1", width/2.0, Rout2+cm, depth/2.0 );
  G4Box *slottemp = new G4Box( "slottemp", width/2.0, 15.5*cm, depth/2.0 + cm );
  G4RotationMatrix *rot_temp = new G4RotationMatrix;

  rot_temp->rotateY(angle);

  G4SubtractionSolid *boxwithslot = new G4SubtractionSolid( "boxwithslot", cutbox1, slottemp, 0, G4ThreeVector(35.0*cm,0,0) );
  //G4LogicalVolume *boxwithslot_log = new G4LogicalVolume( boxwithslot, GetMaterial("Air"), "boxwithslot_log" );
  
  G4double Rbox = dist + depth/2.0;
  
  G4ThreeVector pos_box_withslot( -Rbox*sin(angle) + width/2.0*cos(angle), 0, Rbox*cos(angle) + width/2.0*sin(angle) );

  //new G4PVPlacement( rot_temp, pos_box_withslot, boxwithslot_log, "boxwithslot_phys", worldlog, false, 0 );
  
  G4ThreeVector pos(0,0,0.5*(zstart_lead1+zstop_lead1));
  
  G4ThreeVector posrel_boxwithslot = pos_box_withslot - pos;
  
  G4SubtractionSolid *leadcone1_cut = new G4SubtractionSolid("leadcone1_cut", leadcone1, boxwithslot, rot_temp, posrel_boxwithslot );

  G4LogicalVolume *leadcone1_log = new G4LogicalVolume(leadcone1_cut, GetMaterial("Lead"), "leadcone1_log" );

  leadcone1_log->SetVisAttributes( lead_visatt );
  
  
  //new G4PVPlacement( 0, pos, leadcone1_log, "leadcone1_phys", worldlog, false, 0 );
  
  G4double zsections[3] = {z_outer_magnetic + 74.0*inch,
			   201.632*inch - TargetCenter_zoffset,
			   207.144*inch - TargetCenter_zoffset + 40.0*inch };
  G4double Rin_sections[3] = { 15.0*inch/2.0 + (zsections[0]-zsections[1])*tan(1.5*deg),
			       15.0*inch/2.0,
			       15.0*inch/2.0 };
  G4double Rout_sections[3] = {Rin_sections[0] + 5.*cm,
			       Rin_sections[1] + 5.*cm,
			       Rin_sections[2] + 5.*cm };

  G4Polycone *leadshield2 = new G4Polycone( "leadshield2", 90.0*deg, 180.0*deg, 3, zsections, Rin_sections, Rout_sections );
  G4LogicalVolume *leadshield2_log = new G4LogicalVolume( leadshield2, GetMaterial("Lead"), "leadshield2_log" );

  leadshield2_log->SetVisAttributes( lead_visatt );
  
  //new G4PVPlacement( 0, G4ThreeVector(), leadshield2_log, "leadshield2_phys", worldlog, false, 0 );

  ////// New geometry with vertical wall(s): 

  G4ThreeVector zaxis_temp( -sin(16.9*deg), 0.0, cos(16.9*deg) );
  G4ThreeVector yaxis_temp( 0,1,0);
  G4ThreeVector xaxis_temp = (yaxis_temp.cross(zaxis_temp)).unit();
  
  G4ThreeVector frontcorner_pos = 1.6*m*zaxis_temp;

  G4double zstart_lead_wall1 = z_outer_magnetic + 15*cm;
  G4double zstop_lead_wall1 = zstart_lead_wall1 + 1.25*m;

  G4Box *lead_wall1 = new G4Box("lead_wall1", 5.0*cm/2.0, 31.0*cm/2.0, 1.25*m/2.0 );
  G4LogicalVolume *lead_wall1_log = new G4LogicalVolume( lead_wall1, GetMaterial("Lead"), "lead_wall1_log" );

  G4double xtemp = -( 5.5*inch/2.0 + 1.5*inch + (1.25/2.0+0.15)*m*tan(1.5*deg) + 2.5*cm/cos(1.5*deg) );
  
  rot_temp = new G4RotationMatrix;

  rot_temp->rotateY( 1.5*deg );
  
  new G4PVPlacement( rot_temp, G4ThreeVector( xtemp, 0.0, zstart_lead_wall1 + 0.5*1.25*m ), lead_wall1_log, "lead_wall1_phys", worldlog, false, 0 );

  G4cout << "Lead wall A (x,y,z) = (" << xtemp/cm << ", " << 0.0 << ", " << (zstart_lead_wall1 + 0.5*1.25*m)/cm << ")" << G4endl;
  
  lead_wall1_log->SetVisAttributes( lead_visatt );

  G4double zstart_lead_wall2 = z_formed_bellows + 76.09*inch + 1.71*inch + 15.75*inch + 1.0*inch;
  G4double zstop_lead_wall2 = 207.144*inch - TargetCenter_zoffset + 40.0*inch;

  G4cout << "Lead wall B zstart - zstop = " << (zstop_lead_wall2 - zstart_lead_wall2)/cm << G4endl;
  
  G4double zpos_lead_wall2 = 0.5*(zstart_lead_wall2 + zstop_lead_wall2 );
  //we want x position to have x = 
  G4double xpos_lead_wall2 = -(8.0*inch + 2.5*cm + (zpos_lead_wall2 - 201.632*inch + TargetCenter_zoffset )*tan(1.5*deg));

  G4Box *lead_wall2 = new G4Box("lead_wall2", 5.0*cm/2.0, 24.0*inch/2.0, 0.5*(zstop_lead_wall2 - zstart_lead_wall2) );

  G4LogicalVolume *lead_wall2_log = new G4LogicalVolume( lead_wall2, GetMaterial("Lead"), "lead_wall2_log" );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY( 1.5*deg );

  new G4PVPlacement( rot_temp, G4ThreeVector( xpos_lead_wall2, 0, zpos_lead_wall2 ), lead_wall2_log, "lead_wall2_phys", worldlog, false, 0 );

  G4cout << "Lead wall B (x,y,z) = (" << xpos_lead_wall2/cm << ", " << 0.0 << ", " << zpos_lead_wall2/cm << ")" << G4endl;
  
  lead_wall2_log->SetVisAttributes( lead_visatt );
  
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
  //Notch depth 25 cm
  G4double Beamslot_notch_depth = 25.0*cm;
  
  G4double Beamslot_lead_width = (fDetCon->fHArmBuilder->f48D48width)/2.0-35.0*cm;
  G4double Beamslot_lead_height = 31.0*cm;
  G4double Beamslot_lead_depth = fDetCon->fHArmBuilder->f48D48depth;

  G4Box *Beamslot_lead_box = new G4Box("Beamslot_lead_box", Beamslot_lead_width/2.0, Beamslot_lead_height/2.0, Beamslot_lead_depth/2.0 );
  //G4LogicalVolume *Beamslot_lead_log = new G4LogicalVolume( 

  G4Box *Beamslot_cut_box = new G4Box( "Beamslot_cut_box", 50.0*cm/2.0, Beamslot_lead_height/2.0 + mm, Beamslot_notch_depth );

  G4RotationMatrix *notch_rot = new G4RotationMatrix;
  notch_rot->rotateY( 45.0*deg );

  G4double notch_angle = 45.0*deg;
  
  G4SubtractionSolid *Beamslot_lead_box_cut = new G4SubtractionSolid( "Beamslot_lead_box_cut", Beamslot_lead_box, Beamslot_cut_box, notch_rot,
								      G4ThreeVector( Beamslot_lead_width/2.0,
										     0.0,
										     -Beamslot_lead_depth/2.0 ) );
  
  G4double SBSang = fDetCon->fHArmBuilder->f48D48ang;

  G4RotationMatrix *Beamslot_lead_rm = new G4RotationMatrix;
  Beamslot_lead_rm->rotateY( SBSang );

  G4ThreeVector SBS_zaxis( -sin(SBSang), 0.0, cos(SBSang) );
  G4ThreeVector SBS_yaxis( 0.0, 1.0, 0.0 );
  G4ThreeVector SBS_xaxis( cos(SBSang), 0.0, sin(SBSang) );

  G4double Beamslot_lead_xoffset = 35.0*cm + Beamslot_lead_width/2.0;

  G4ThreeVector Beamslot_lead_position = (fDetCon->fHArmBuilder->f48D48dist + (fDetCon->fHArmBuilder->f48D48depth)/2.0) * SBS_zaxis + Beamslot_lead_xoffset * SBS_xaxis;

  //Define a subtraction cone for the beam slot in the SBS magnet:
  // G4double zbeampipe[2] = {162.2*cm, 592.2*cm};
  // G4double rinbeampipe[2] = {0.0*cm, 0.0*cm};

  G4double beampipe_subtraction_cone_dz = ((592.2 - 162.2)/2.0)*cm;
  G4double beampipe_subtraction_cone_zpos = ((592.2+162.2)/2.0)*cm;

  G4ThreeVector beampipe_subtraction_cone_position( 0.0, 0.0, beampipe_subtraction_cone_zpos );

  G4Cons *beampipe_subtraction_cone = new G4Cons( "beampipe_subtraction_cone", 0.0*cm, 5.1*cm, 0.0*cm, 15.1*cm, beampipe_subtraction_cone_dz, 0.0*deg, 360.0*deg );

  G4RotationMatrix *Beamslot_lead_rm_inv = new G4RotationMatrix;
  Beamslot_lead_rm_inv->rotateY( -SBSang );

  G4ThreeVector beampipe_beamslot_relative_position = beampipe_subtraction_cone_position - Beamslot_lead_position;

  G4ThreeVector beampipe_beamslot_relative_position_local( beampipe_beamslot_relative_position.dot( SBS_xaxis ), 
							   beampipe_beamslot_relative_position.dot( SBS_yaxis ), 
							   beampipe_beamslot_relative_position.dot( SBS_zaxis ) );

  //The subtraction that we want to perform is Beamslot lead box - beampipe cone:
  G4SubtractionSolid *Beamslot_lead_with_hole = new G4SubtractionSolid( "Beamslot_lead_with_hole", Beamslot_lead_box_cut, beampipe_subtraction_cone, Beamslot_lead_rm_inv, beampipe_beamslot_relative_position_local );

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

  G4SubtractionSolid *leadcone_cut = new G4SubtractionSolid( "leadcone_cut", leadcone, Beamslot_lead_box_cut, Beamslot_lead_rm, cutbox_relative_position );
  
  G4LogicalVolume *leadcone_cut_log = new G4LogicalVolume( leadcone_cut, GetMaterial("Lead"), "leadcone_cut_log" );
  G4PVPlacement *leadcone_cut_pv = new G4PVPlacement( 0, leadcone_global_position, leadcone_cut_log, "leadcone_cut_pv", worldlog, 0, false, 0 );

  // G4Polycone *SIDISlead_cone = new G4Polycone( "SIDISlead_cone", 0.0*deg, 360.0*deg, nz, zlead, rinlead, routlead );
  // G4LogicalVolume *SIDISlead_log = new G4LogicalVolume( SIDISlead_cone, GetMaterial("Lead"), "SIDISlead_log", 0, 0, 0 );
  // G4PVPlacement *SIDISlead_pv = new G4PVPlacement( 0, G4ThreeVector(), SIDISlead_log, "SIDISlead_pv", worldlog, false, 0 );
  
  //We also want to put some lead and/or Fe shielding, i.e., a "collimator" in front of the SBS magnet gap:

  // double SBScollwidth = 469.9*mm;
  // double SBScollheight = 1219.2*mm;
  // double SBScolldepth = 10.0*cm;
  
  // double coilspace = 214.5*mm + 20.63*mm;
  // double SBS_coll_R = fDetCon->fHArmBuilder->f48D48dist - coilspace - SBScolldepth/2.0 - 5.0*cm;

  // // double SBS_coll_gapwidth = 50.0*cm*SBS_coll_R/(fDetCon->fHArmBuilder->fRICHdist - 0.5*m) + fDetCon->fTargetBuilder->GetTargLen()*sin( SBSang );
  // // double SBS_coll_gapheight = 200.0*cm*SBS_coll_R/(fDetCon->fHArmBuilder->fRICHdist - 0.5*m)
  // G4double SBS_coll_gapwidth =
  

  // G4Box *SBScoll = new G4Box("SBScoll", 1.5*SBScollwidth, 1.5*SBScollheight, SBScolldepth/2.0 );
  // G4Box *SBScoll_hole = new G4Box("SBScoll_hole", SBS_coll_gapwidth/2.0, SBS_coll_gapheight/2.0, SBScolldepth/2.0+1.0*cm );

  // G4Cons *SBScoll_cutcone = new G4Cons("SBScoll_cutcone", 0.0*cm, rinstart+5.0*cm, 0.0, rinend + 5.0*cm, (zend-zstart)/2.0, 0.0*deg, 360.0*deg );

  

  // G4ThreeVector SBScoll_pos( SBS_coll_R*sin(SBSang), 0.0, SBS_coll_R*cos(SBSang) );
  // G4ThreeVector cutcone_relative_pos = leadcone_global_position - SBScoll_pos;

  // G4ThreeVector cutcone_relative_pos_local( cutcone_relative_pos.dot(SBS_xaxis),
  // 					    cutcone_relative_pos.dot(SBS_yaxis),
  // 					    cutcone_relative_pos.dot(SBS_zaxis) );

  // G4SubtractionSolid *SBS_collimator = new G4SubtractionSolid( "SBS_collimator", SBScoll, SBScoll_hole );
  // G4SubtractionSolid *SBS_collimator_beamcut = new G4SubtractionSolid("SBS_collimator_beamcut", SBS_collimator, SBScoll_cutcone, Beamslot_lead_rm_inv, 
  // 								      cutcone_relative_pos_local );

  // G4LogicalVolume *SBS_collimator_log = new G4LogicalVolume( SBS_collimator_beamcut, GetMaterial("Lead"), "SBS_collimator_log" );
    
  // new G4PVPlacement( Beamslot_lead_rm, G4ThreeVector( SBS_coll_R*sin(SBSang), 0.0, SBS_coll_R*cos(SBSang) ), SBS_collimator_log, "SBS_collimator_phys", worldlog, 
  // 		     0, false, 0 );

  G4VisAttributes *leadVisAtt= new G4VisAttributes(G4Colour(0.25,0.25,0.25));
  //SIDISlead_log->SetVisAttributes(leadVisAtt);

  Beamslot_lead_log->SetVisAttributes(leadVisAtt);
  leadcone_cut_log->SetVisAttributes(leadVisAtt);
  //SBS_collimator_log->SetVisAttributes(leadVisAtt);
}

//++++++++++++++++++++++++++++++++++++++++++  correcting magnets ++++ 09/23/2014

/////////////////////////////////////////////////////////////////////
//  09/16/2014  implement of correcting magnet ob beam line 

//void G4SBSBeamlineBuilder::MakeGEpCorrMag(G4LogicalVolume *worldlog){
// }


///////////////////////////////////////////////////////////////////////////



