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
#include "G4RotationMatrix.hh"

#include "G4SystemOfUnits.hh"


G4SBSBeamlineBuilder::G4SBSBeamlineBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
  assert(dc);
}

G4SBSBeamlineBuilder::~G4SBSBeamlineBuilder(){;}

void G4SBSBeamlineBuilder::BuildComponent(G4LogicalVolume *worldlog){
  Targ_t targtype = fDetCon->fTargType;

  if( fDetCon->fExpType == kGEp ){
    /////////////////////////////////////////////////////////////////////
    //  09/16/2014  implement of correcting magnet ob beam line 
    // AJP 02/23/2015: Need experiment-dependent positioning of correction magnets: This configuration is technically only applicable to GEp 12 GeV^2 configuration! 
    //
    // 
    //  Entry iron tube to shield beam from stray field
    //
   
    // *****6/17

   //  G4double tRmin = 5.*cm;
   //  G4double tRmax = 7.*cm;
   //  G4double tDzz = 20.*cm/2;
   //  G4double tSPhi = 0.*deg;
   //  G4double tDphi = 360.*deg;


   //  G4Tubs *iron_ent_tube = new G4Tubs("iron_ent_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);

   //  G4LogicalVolume *iron_ent_tube_log = new G4LogicalVolume(iron_ent_tube, GetMaterial("Iron"), "iron_ent_tube", 0, 0, 0 );

   // /new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, 110.936*cm + tDzz), iron_ent_tube_log, "iron_ent_tube_phys", worldlog,false,0);

   //  //  next entrance correction magnet

   //  G4Box* box_1 = new G4Box("EnMag_1",50.*cm/2., 50.*cm/2., 16.*cm/2.);
   //  G4Box* box_2 = new G4Box("EnMag_2",20.*cm/2., 30.*cm/2., 17.*cm/2.); // aperture to pass beam
   //  G4SubtractionSolid* EnMag = new G4SubtractionSolid("EnMag", box_1, box_2);   
   //  G4LogicalVolume * EnMag_log = new G4LogicalVolume(EnMag , GetMaterial("Iron"), "EnMag", 0, 0, 0); 
   //  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, 134.935*cm+8.*cm), EnMag_log, "EnMag_phys", worldlog, false, 0);

   //  // exit correction magnet 

   //  G4Box* box_3 = new G4Box("ExtMag_1",50.*cm/2., 54.*cm/2., 40.*cm/2.);
   //  G4Box* box_4 = new G4Box("ExtMag_2",20.*cm/2., 40.*cm/2., 41.*cm/2.); // aperture to pass beam
   //  G4SubtractionSolid* ExtMag = new G4SubtractionSolid("ExtMag", box_3, box_4);   
   //  G4LogicalVolume * ExtMag_log = new G4LogicalVolume(ExtMag , GetMaterial("Iron"), "ExtMag", 0, 0, 0); 
   //  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, 351.435*cm+20.*cm), ExtMag_log, "ExtMag_phys", worldlog, false, 0);

   //  // iron conical tube on the beamline inside SBS magnet split

   //  G4double tcRmin1 = 7.4*cm;
   //  G4double tcRmax1 = 8.68*cm;
   //  G4double tcRmin2 = 10.55*cm;
   //  G4double tcRmax2 = 11.82*cm;
   //  G4double tcDzz = 120.*cm/2;
   //  G4double tcSPhi = 0.*deg;
   //  G4double tcDphi = 360.*deg;
  
   //  G4Cons *iron_con_tube = new G4Cons("iron_con_tube", tcRmin1, tcRmax1,tcRmin2, tcRmax2, tcDzz, tcSPhi, tcDphi);

   //  G4LogicalVolume *iron_con_tube_log = new G4LogicalVolume(iron_con_tube, GetMaterial("Iron"), "iron_con_tube", 0, 0, 0 );

   //  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, 170.944*cm + tcDzz), iron_con_tube_log, "iron_con_tube_phys", worldlog,false,0);

   //  G4VisAttributes *McorrVisAtt= new G4VisAttributes(G4Colour(0.9,0.9,0.9));
   //  iron_con_tube_log->SetVisAttributes(McorrVisAtt);
   //  ExtMag_log->SetVisAttributes(McorrVisAtt);
   //  EnMag_log->SetVisAttributes(McorrVisAtt);
   //  iron_ent_tube_log->SetVisAttributes(McorrVisAtt);

    // *****6/17
    
    // Code Imported from BeamLine.F (Sergey) - 6/18/2015

    MakeGEpBeamline(worldlog);

  }

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

      /*  // my cancel Be window for GEp experiment 09/29/2014
	  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, ent_len/2-winthick/2), ent_winlog, "entwin_phys", entvacLog,false,0);
      */  // my cancel Be window for GEp experiment 09/29/2014

      ent_winlog->SetVisAttributes(new G4VisAttributes(G4Colour(0.7,1.0,0.0)));
    } // else {
    //   //Don't add window: we want the beam to interact with the target first. Butt up against the outer edge of the scattering chamber:
    //   new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner), entLog_cut, "ent_phys", worldlog, false,0);
    //new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner), entvacLog_cut, "entvac_phys", worldlog,false,0);
    // }
  } else {
    // Cryotarget - up against the chamber wall

    //*****6/17
    //new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner), entLog_cut, "ent_phys", worldlog, false,0);
    //new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner), entvacLog_cut, "entvac_phys", worldlog,false,0);
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
  G4double exit_z[]   = { 162.2*cm, 592.2*cm, 609.84*cm,609.85*cm, 1161.02*cm, 1161.03*cm,2725.66*cm };
  G4double exit_z_vac[] = { 162.5*cm, 592.5*cm, 610.24*cm,610.35*cm, 1161.52*cm, 1161.53*cm,2726.46*cm };

  G4double exit_zero[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  G4double exit_rin[] = { 4.8*cm, 14.8*cm,15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };
  G4double exit_rou[] = { 5.0*cm, 15.0*cm,15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };


  G4Polycone *ext_cone = new G4Polycone("ext_cone", 0.0*deg, 360.0*deg, nsec, exit_z, exit_rin, exit_rou);
  G4Polycone *ext_vac  = new G4Polycone("ext_vac ", 0.0*deg, 360.0*deg, nsec, exit_z_vac, exit_zero, exit_rin);

  G4LogicalVolume *extLog = new G4LogicalVolume(ext_cone, GetMaterial("Aluminum"), "ext_log", 0, 0, 0);
  G4LogicalVolume *extvacLog = new G4LogicalVolume(ext_vac, GetMaterial("Vacuum"), "extvac_log", 0, 0, 0);

  //*****6/17
  //new G4PVPlacement(0,G4ThreeVector(), extLog, "ext_phys", worldlog, false,0);
  //new G4PVPlacement(0,G4ThreeVector(), extvacLog, "extvac_phys", worldlog,false,0);
  //*****6/17

  G4VisAttributes *extVisAtt= new G4VisAttributes(G4Colour(0.9,0.1,0.9));
  extLog->SetVisAttributes(extVisAtt);


  // Seal this up if we have a gas target
  if( fDetCon->fTargType == kH2 || fDetCon->fTargType == k3He || fDetCon->fTargType == kNeutTarg ){
    // Add in exit Al window

    double extwin_thick = 5.0e-4*cm;

    G4Tubs *extwin = new G4Tubs("ext_win", 0.0, exit_rin[0], extwin_thick/2, 0.*deg, 360.*deg );
    G4LogicalVolume *ext_winlog = new G4LogicalVolume(extwin, GetMaterial("Aluminum"), "entwin_log", 0, 0, 0);


    /*  cancel Al window to vacuum pipe 10/02/14
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, exit_z[0] - extwin_thick/2), ext_winlog, "extwin_phys", worldlog,false,0);
    */   //cancel Al window to vacuum pipe 10/02/14


    ext_winlog->SetVisAttributes(new G4VisAttributes(G4Colour(0.5,0.2,0.6)));
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

  G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.2,0.6,0.2));

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

// GEp Beamline Construction --- following Sergey's Fortran code
void G4SBSBeamlineBuilder::MakeGEpBeamline(G4LogicalVolume *worldlog) {

    G4double tRmin = 0.0*cm;
    G4double tRmax = 0.5*17.14*cm;
    G4double tDzz  = 0.5*2.13*cm;
    G4double tSPhi = 0.0*deg;
    G4double tDphi = 360.0*deg;

    G4double X=0, Y=0, Z=0;
    G4ThreeVector zero(0.0, 0.0, 0.0);

    // BeamLine-tube main with flanges
    // Flange 1:
    G4Tubs *FLN1_tube = new G4Tubs("FLN1_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);

    // Bore out middle, fill with vaccum
    tRmax = 0.5*8.93*cm;
    tDzz += 0.01*cm;
    G4Tubs *FVL1_tube = new G4Tubs("FVL1_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
    G4SubtractionSolid *FLN1_sub = new G4SubtractionSolid( "FLN1_sub", FLN1_tube, FVL1_tube, 0, zero );
    
    // Convert into logical volumes
    G4LogicalVolume *FLN1_log = new G4LogicalVolume( FLN1_sub, GetMaterial("Iron"), "FLN1_log" );
    G4LogicalVolume *FVL1_log = new G4LogicalVolume( FVL1_tube, GetMaterial("Vacuum"), "FVL1_log");
    // Then place the vacuum inside the Iron Tube
    Z = (159.51 + 0.5*2.13)*cm;
    new G4PVPlacement( 0, zero, FVL1_log, "Flange1_Vac", FLN1_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), FLN1_log, "Flange1_Iron", worldlog, false, 0 );

    // Conical tube for vacuum
    G4double pDz    = (352.636 - 2.13*0.5 - 2.84*0.5)*0.5*cm;
    G4double pRmin1 = 0.0*cm;
    G4double pRmax1 = (0.5*8.93 + 0.318)*cm;
    G4double pRmin2 = 0.0*cm;
    G4double pRmax2 = (0.5*27.62 + 0.318)*cm;

    G4Cons *TB01_cons = new G4Cons( "TB01_cons", pRmin1, pRmax1, pRmin2, pRmax2, pDz, tSPhi, tDphi );
    
    // Bore out middle, fill with vacuum
    pRmax1 = 0.5*8.93*cm;
    pRmax2 = 0.5*27.62*cm;
    pDz += 0.01*cm;
    G4Cons *TVB1_cons = new G4Cons( "TVB1_cons", pRmin1, pRmax1, pRmin2, pRmax2, pDz, tSPhi, tDphi );
    G4SubtractionSolid *TB01_sub = new G4SubtractionSolid( "TB01_sub", TB01_cons, TVB1_cons, 0, zero );

    // Convert into logical volumes
    G4LogicalVolume *TB01_log = new G4LogicalVolume( TB01_sub, GetMaterial("Iron"), "TB01_log" );
    G4LogicalVolume *TVB1_log = new G4LogicalVolume( TVB1_cons, GetMaterial("Vacuum"), "TVB1_log" );
    // Then place the vacuum inside the Iron Cone
    Z = (159.51 + 2.13)*cm + pDz;
    new G4PVPlacement( 0, zero, TVB1_log, "Cone1_Vac", TB01_log, false, 0);
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TB01_log, "Cone1_Iron", worldlog, false, 0 );

    // Flange 2:
    tRmin = 0.0*cm;
    tRmax = 0.5*35.56*cm;
    tDzz  = 0.5*2.84*cm;
    G4Tubs *FLN2_tube = new G4Tubs("FLN2_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);

    // Bore out middle, fill with vaccum
    tRmax = 0.5*27.62*cm;
    tDzz += 0.01*cm;
    G4Tubs *FVL2_tube = new G4Tubs("FVL2_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
    G4SubtractionSolid *FLN2_sub = new G4SubtractionSolid( "FLN2_sub", FLN2_tube, FVL2_tube, 0, zero );

    // Convert into logical volumes
    G4LogicalVolume *FLN2_log = new G4LogicalVolume( FLN2_sub, GetMaterial("Iron"), "FLN2_log" );
    G4LogicalVolume *FVL2_log = new G4LogicalVolume( FVL2_tube, GetMaterial("Vacuum"), "FVL2_log");
    // Then place the vacuum inside the Iron Tube
    Z = (159.51 + 352.636 - 2.84*0.5)*cm;
    new G4PVPlacement( 0, zero, FVL2_log, "Flange2_Vac", FLN2_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), FLN2_log, "Flange2_Iron", worldlog, false, 0 );

    //last second big flange the same an placed at distance 4.237"
    // ??? I am assuming this comment means "use flange2 again, place at distance of 4.237 inches "
    Z = (159.51+352.636-2.84*0.5+4.237*2.54)*cm;
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), FLN2_log, "Flange7_Iron", worldlog, false, 0 );

    // Here a bellow and we assign wall of 0.03 cm
    tRmin = 0.0*cm;
    tRmax = (0.5*27.62 + 0.03)*cm;
    tDzz  = 0.5*(4.237*2.54 - 2.84)*cm;

    G4Tubs *TBL8_tube = new G4Tubs("TBL8_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);

    // Bore out middle, fill with vaccum
    tRmax = (0.5*27.62)*cm;
    tDzz += 0.01*cm;
    G4Tubs *TVL8_tube = new G4Tubs("TVL8_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
    G4SubtractionSolid *TBL8_sub = new G4SubtractionSolid( "TBL8_sub", TBL8_tube, TVL8_tube, 0, zero );

    // Convert into logical volumes
    G4LogicalVolume *TBL8_log = new G4LogicalVolume( TBL8_sub, GetMaterial("Iron"), "TBL8_log" );
    G4LogicalVolume *TVL8_log = new G4LogicalVolume( TVL8_tube, GetMaterial("Vacuum"), "TVL8_log" );
    // Then place the vacuum inside the Iron Tube
    Z = (159.51 + 352.636 - 2.84*0.5 + 4.237*2.54*0.5)*cm;
    new G4PVPlacement( 0, zero, TVL8_log, "Bellow_Vac", TBL8_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TBL8_log, "Bellow_Iron", worldlog, false, 0 );

    // EXTEND VACUUM LINE by using Maduka geometry
    // ============================================
    tRmin = 0.0*cm;
    tRmax = 13.0*2.54*0.5*cm;
    tDzz  = 0.5*41.0*2.54*cm;

    G4Tubs *TBL9_tube = new G4Tubs("TBL9_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);

    // Bore out middle, fill with vaccum
    tRmax = 0.5*12.0*2.54*cm;
    tDzz += 0.01*cm;
    G4Tubs *TVL9_tube = new G4Tubs("TVL9_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
    G4SubtractionSolid *TBL9_sub = new G4SubtractionSolid( "TBL9_sub", TBL9_tube, TVL9_tube, 0, zero );

    // Convert into logical volumes
    G4LogicalVolume *TBL9_log = new G4LogicalVolume( TBL9_sub, GetMaterial("Aluminum"), "TBL9_log" );
    G4LogicalVolume *TVL9_log = new G4LogicalVolume( TVL9_tube, GetMaterial("Vacuum"), "TVL9_log" );
    // Then place the vacuum inside the Al Tube
    Z = (159.51 + 352.636 - 2.84*0.5 + 4.237*2.54 + 41.0*2.54*0.5)*cm;
    new G4PVPlacement( 0, zero, TVL9_log, "Extended_Vac1", TBL9_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TBL9_log, "Extended_Al1", worldlog, false, 0 );

    tRmin = 0.0*cm;
    tRmax = 25.0*2.54*0.5*cm;
    tDzz  = 0.5*217.0*2.54*cm;
    G4Tubs *TML9_tube = new G4Tubs( "TML9_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);

    // Bore out middle, fill with vaccum
    tRmax = 0.5*24.0*2.54*cm;
    tDzz += 0.01*cm;
    G4Tubs *TMV9_tube = new G4Tubs("TMV9_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
    G4SubtractionSolid *TML9_sub = new G4SubtractionSolid( "TML9_sub", TML9_tube, TMV9_tube, 0, zero );

    // Convert into logical volumes
    G4LogicalVolume *TML9_log = new G4LogicalVolume( TML9_sub, GetMaterial("Aluminum"), "TML9_log" );
    G4LogicalVolume *TMV9_log = new G4LogicalVolume( TMV9_tube, GetMaterial("Vacuum"), "TMV9_log" );
    // Then place vac inside of Al tube
    Z = (159.51 + 352.636 - 2.84*0.5 + 4.237*2.54 + 41.0*2.54 + 0.5*217.0*2.54)*cm;
    new G4PVPlacement( 0, zero, TMV9_log, "Extended_Vac2", TML9_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TML9_log, "Extended_Al2", worldlog, false, 0 );

    //-----------------------------------------------------
    //       magnetic tubes 
    pDz    = 120.9*0.5*cm;
    pRmin1 = 0.5*11.06*cm;
    pRmax1 = (0.5*11.06 + 0.635)*cm;
    pRmin2 = 0.5*17.4*cm;
    pRmax2 = (0.5*17.4 + 0.635)*cm;

    G4Cons *TBM1_cons = new G4Cons( "TBM1_cons", pRmin1, pRmax1, pRmin2, pRmax2, pDz, tSPhi, tDphi );
    G4LogicalVolume *TBM1_log = new G4LogicalVolume( TBM1_cons, GetMaterial("Iron"), "TBM1_log" );
    Z = (187.41+120.9*0.5)*cm;
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TBM1_log, "Magnetic_Tube1", worldlog, false, 0 );

    pDz    = 132.96*0.5*cm;
    pRmin1 = 0.5*20.9*cm;
    pRmax1 = (0.5*20.9 + 0.635)*cm;
    pRmin2 = 0.5*27.87*cm;
    pRmax2 = (0.5*27.87 + 0.635)*cm;

    G4Cons *TBM2_cons = new G4Cons( "TBM2_cons", pRmin1, pRmax1, pRmin2, pRmax2, pDz, tSPhi, tDphi );
    G4LogicalVolume *TBM2_log = new G4LogicalVolume( TBM2_cons, GetMaterial("Iron"), "TBM2_log" );
    Z = (187.41 + 187.96 + 132.96*0.5)*cm;
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TBM2_log, "Magnetic_Tube2", worldlog, false, 0 );

    //*       magnetic shield from rings  thickness reduced by factor 0.819
    //*       thickness 12.7 mm * 0.819 = 10.4 mm , we need reduce outer radius
    //*       by (12.7-10.4)/2 = 1.15 mm and increase internal by 1.15 mm
    pDz    = 131.13*0.5*cm;
    pRmin1 = (0.5*13.91 + 0.115)*cm;
    pRmax1 = (0.5*13.91 + 1.27 - 0.115)*cm;
    pRmin2 = (0.5*20.77 + 0.115)*cm;
    pRmax2 = (0.5*20.77 + 1.27 - 0.115)*cm;

    G4Cons *TBM3_cons = new G4Cons( "TBM3_cons", pRmin1, pRmax1, pRmin2, pRmax2, pDz, tSPhi, tDphi );
    G4LogicalVolume *TBM3_log = new G4LogicalVolume( TBM3_cons, GetMaterial("Iron"), "TBM3_log" );
    Z = (182.33 + 131.13*0.5)*cm;
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TBM3_log, "Magnetic_Tube3", worldlog, false, 0 );

    pDz    = 132.08*0.5*cm;
    pRmin1 = (0.5*23.75 + 0.115)*cm;
    pRmax1 = (0.5*23.75 + 1.27 - 0.115)*cm;
    pRmin2 = (0.5*30.88 + 0.115)*cm;
    pRmax2 = (0.5*30.88 + 1.27 - 0.115)*cm;

    G4Cons *TBM4_cons = new G4Cons( "TBM4_cons", pRmin1, pRmax1, pRmin2, pRmax2, pDz, tSPhi, tDphi );
    G4LogicalVolume *TBM4_log = new G4LogicalVolume( TBM4_cons, GetMaterial("Iron"), "TBM4_log" );
    Z = (182.33 + 187.96 + 132.08*0.5)*cm;
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TBM4_log, "Magnetic_Tube4", worldlog, false, 0 );

    //*------------------------------------------------------
    //*        flanges close to the scattered chamber
    tRmin = 0.0*cm;
    tRmax = 0.5*6.0*2.54*cm;
    tDzz  = 0.5*0.84*2.54*cm;
    G4Tubs *FLN3_tube = new G4Tubs("FLN3_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);

    tRmax = 0.5*3.81*2.54*cm;
    tDzz += 0.01*cm;
    G4Tubs *FVL3_tube = new G4Tubs("FVL3_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
    G4SubtractionSolid *FLN3_sub = new G4SubtractionSolid( "FLN3_sub", FLN3_tube, FVL3_tube, 0, zero );
    G4LogicalVolume *FLN3_log = new G4LogicalVolume( FLN3_sub, GetMaterial("Iron"), "FLN3_log" );
    G4LogicalVolume *FVL3_log = new G4LogicalVolume( FVL3_tube, GetMaterial("Vacuum"), "FVL3_log");

    Z = (133.2 + 0.5*0.84*2.54)*cm;
    new G4PVPlacement( 0, zero, FVL3_log, "Flange3_Vac", FLN3_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), FLN3_log, "Flange3_Iron", worldlog, false, 0 );


    tRmin = 0.0*cm;
    tRmax = 0.5*3.81*2.54*cm;
    tDzz  = 0.5*0.84*2.54*cm;
    G4Tubs *FVL5_tube = new G4Tubs("FVL5_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
    G4LogicalVolume *FVL5_log = new G4LogicalVolume( FVL5_tube, GetMaterial("Vacuum"), "FVL5_log" );
    Z = (148.4 - 0.5*0.84*2.54)*cm;
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), FVL5_log, "Flange5_Vac", worldlog, false, 0 );
    
    G4Tubs *FVL6_tube = new G4Tubs("FVL6_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
    G4LogicalVolume *FVL6_log = new G4LogicalVolume( FVL6_tube, GetMaterial("Vacuum"), "FVL6_log" );
    Z = (148.4 - 0.5*0.84*2.54 + 0.84*2.54)*cm;
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), FVL6_log, "Flange5_Vac", worldlog, false, 0 );


    tRmin = 0.0*cm;
    tRmax = 0.5*17.14*cm;
    tDzz  = 0.5*2.13*cm;
    G4Tubs *FLN7_tube = new G4Tubs("FLN7_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);

    tRmax = 0.5*8.93*cm;
    tDzz += 0.01*cm;
    G4Tubs *FVL7_tube = new G4Tubs("FVL7_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
    G4SubtractionSolid *FLN7_sub = new G4SubtractionSolid( "FLN7_sub", FLN7_tube, FVL7_tube, 0, zero );
    G4LogicalVolume *FLN7_log = new G4LogicalVolume( FLN7_sub, GetMaterial("Iron"), "FLN7_log" );
    G4LogicalVolume *FVL7_log = new G4LogicalVolume( FVL7_tube, GetMaterial("Vacuum"), "FVL7_log");

    Z = (159.51 + 2.13*0.5 - 2.13)*cm;
    new G4PVPlacement( 0, zero, FVL7_log, "Flange7_Vac", FLN7_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), FLN7_log, "Flange7_Iron", worldlog, false, 0 );

    //*      we need to extend vacuum after FLN flange by using Maduka geometry
    //*      of vacuum pipe

    //* Describing tubes between flanges  
    tRmin = 0.0*cm;
    tRmax = 4.00*2.54*0.5*cm;
    tDzz  = (3.25*2.54 - 0.84*2.54)*0.5*cm;
    G4Tubs *TBT1_tube = new G4Tubs("TBT1_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);

    tRmax = 3.76*2.54*0.5*cm;
    tDzz += 0.01*cm;
    G4Tubs *TTV1_tube = new G4Tubs("TTV1_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
    G4SubtractionSolid *TBT1_sub = new G4SubtractionSolid( "TBT1_sub", TBT1_tube, TTV1_tube, 0, zero );
    G4LogicalVolume *TBT1_log = new G4LogicalVolume( TBT1_sub, GetMaterial("Iron"), "TBT1_log" );
    G4LogicalVolume *TTV1_log = new G4LogicalVolume( TTV1_tube, GetMaterial("Vacuum"), "TTV1_log");

    Z = (148.4 + 0.84*2.54 + (3.25*2.54 - 0.84*2.54)*0.5)*cm;
    new G4PVPlacement( 0, zero, TTV1_log, "TTV1_Vac", TBT1_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TBT1_log, "TBT1_Iron", worldlog, false, 0 ); // between large and small flanges
  
    tRmin = 0.0*cm;
    tRmax = (3.81*2.54*0.5 + 0.03)*cm;
    tDzz  = (6.0 - 0.84 - 0.84)*2.54*0.5*cm;
    G4Tubs *TBT2_tube = new G4Tubs("TBT2_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);

    tRmax = 3.81*2.54*0.5*cm;
    tDzz += 0.01*cm;
    G4Tubs *TTV2_tube = new G4Tubs("TTV2_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
    G4SubtractionSolid *TBT2_sub = new G4SubtractionSolid( "TBT2_sub", TBT2_tube, TTV2_tube, 0, zero );
    G4LogicalVolume *TBT2_log = new G4LogicalVolume( TBT2_sub, GetMaterial("Iron"), "TBT2_log" );
    G4LogicalVolume *TTV2_log = new G4LogicalVolume( TTV2_tube, GetMaterial("Vacuum"), "TTV2_log");

    Z = (133.2 + 0.84*2.54)*cm + tDzz;
    new G4PVPlacement( 0, zero, TTV2_log, "TTV2_Vac", TBT2_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TBT2_log, "TBT2_Iron", worldlog, false, 0 );

    //C------BeamLine-Vacuum----------------------------------------------------------
    //C-------------------------------------------------------------------------------


    //*  ------ Shield of Beam Line tube------------
    G4double deg2rad = 3.1415926/180.0;

    G4double BLLength = 600.0;                                //cm
    G4double BLDist = 53.73*2.54;                             //cm
    G4double BLREntMin = 2.415*2.54;                          //cm  
    G4double BLThick = 1.0;                                   //cm
    G4double BLAngle = 1.5*(3.141592/180.0);                  //rad
    G4double BLRExitMin = BLREntMin + BLLength*tan(BLAngle);  //cm
    G4double BLREntMax = BLREntMin + BLThick;                 //cm
    G4double BLRExitMax = BLRExitMin + BLThick;               //cm

    G4double th1 = (90.0+1.5)*deg2rad;
    G4double ph1 = 0.0*deg2rad;
    G4double th2 = 90.0*deg2rad;
    G4double ph2 = 90.0*deg2rad;
    G4double th3 = 1.5*deg2rad;
    G4double ph3 = 0.0*deg2rad;

    G4ThreeVector u( sin(th1)*cos(ph1), sin(th1)*sin(ph1), cos(ph1) );
    G4ThreeVector v( sin(th2)*cos(ph2), sin(th2)*sin(ph2), cos(ph2) );
    G4ThreeVector w( sin(th3)*cos(ph3), sin(th3)*sin(ph3), cos(ph3) );
    G4RotationMatrix *rotm1 = new G4RotationMatrix(u,v,w);

    G4double  dx1 = 2.5*cm;
    G4double  dx2 = 2.5*cm;
    G4double  dy1 = (BLREntMax + 5.0 + 1.83 + 0.4 + 5.0)*cm;
    G4double  dy2 = (BLRExitMax + 5.0 + 5.0)*cm;
    G4double  dz  = (0.5*(BLLength - 70.0) - 7.0)*cm;

    G4Trd *BLSH = new G4Trd( "BLSH", dx1, dx2, dy1, dy2, dz );
    G4LogicalVolume *BLSH_log = new G4LogicalVolume( BLSH, GetMaterial("Lead"), "BLSH_log" );
   
    X = (dy1+dy2)*0.5 + 10.0*cm;
    Z = (BLDist + 0.5*(BLLength+70.0) + 7.0)*cm;
    //new G4PVPlacement( rotm1, G4ThreeVector(X, Y, Z), BLSH_log, "BLSH_Lead", worldlog, false, 0 );

    // VISUALS
    G4VisAttributes *ironColor= new G4VisAttributes(G4Colour(0.52,0.47,0.47));
    FLN1_log->SetVisAttributes( ironColor );
    FLN2_log->SetVisAttributes( ironColor );
    FLN3_log->SetVisAttributes( ironColor );
    FLN7_log->SetVisAttributes( ironColor );
    TB01_log->SetVisAttributes( ironColor );
    TBL8_log->SetVisAttributes( ironColor );
    TBM1_log->SetVisAttributes( ironColor );
    TBM2_log->SetVisAttributes( ironColor );
    TBM3_log->SetVisAttributes( ironColor );
    TBM4_log->SetVisAttributes( ironColor );
    TBT1_log->SetVisAttributes( ironColor );
    TBT2_log->SetVisAttributes( ironColor );

    G4VisAttributes *AlColor= new G4VisAttributes(G4Colour(0.4,0.4,0.4));
    TBL9_log->SetVisAttributes( AlColor );
    TML9_log->SetVisAttributes( AlColor );

    G4VisAttributes *leadColor= new G4VisAttributes(G4Colour(0.9,0.9,0.9));
    BLSH_log->SetVisAttributes( leadColor );

    // Vacuum
    FVL1_log->SetVisAttributes( G4VisAttributes::Invisible );
    FVL2_log->SetVisAttributes( G4VisAttributes::Invisible );
    FVL3_log->SetVisAttributes( G4VisAttributes::Invisible );
    FVL5_log->SetVisAttributes( G4VisAttributes::Invisible );
    FVL6_log->SetVisAttributes( G4VisAttributes::Invisible );
    FVL7_log->SetVisAttributes( G4VisAttributes::Invisible );
    TVB1_log->SetVisAttributes( G4VisAttributes::Invisible );
    TVL8_log->SetVisAttributes( G4VisAttributes::Invisible );
    TVL9_log->SetVisAttributes( G4VisAttributes::Invisible );
    TMV9_log->SetVisAttributes( G4VisAttributes::Invisible );
    TTV1_log->SetVisAttributes( G4VisAttributes::Invisible );
    TTV2_log->SetVisAttributes( G4VisAttributes::Invisible );
}

//  Here is lead shield of beam line for GEp

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

  /*    remove magnet conical 09/30/2014 
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, leadstart + magleadlen/2), leadinmag_log, "leadinmag_phys", worldlog,false,0);

  */   // 09/30/2014 we will shield by lead bloks around beam line


  // Lead from magnet on
  // 311 cm -> 592 cm
  leadstart = 311*cm;
  leadend   = 592*cm;
  magleadlen = leadend-leadstart;

  G4Tubs *leadtube= new G4Tubs( "leadtube",  0*cm, maxrad, magleadlen/2, 0.*deg, 360*deg );
  G4SubtractionSolid *leadafter = new G4SubtractionSolid("lead_after", leadtube, ext_cone, 0, G4ThreeVector(0.0, 0.0, -leadstart-magleadlen/2 ) );

  G4LogicalVolume *leadafter_log = new G4LogicalVolume( leadafter, GetMaterial("Lead"), "leadafter_log", 0, 0, 0 );

  /*   cancell 09/30/2014  we will shield by lead bloks around beam line
       new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, leadstart + magleadlen/2), leadafter_log, "leadafter_phys", worldlog,false,0);
  */   // cancell 09/30/2014


  G4VisAttributes *leadVisAtt= new G4VisAttributes(G4Colour(0.9,0.9,0.9));
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
  //Notch depth 25 cm
  G4double Beamslot_notch_depth = 25.0*cm;
  
  G4double Beamslot_lead_width = (fDetCon->fHArmBuilder->f48D48width)/2.0-35.0*cm;
  G4double Beamslot_lead_height = 31.0*cm;
  G4double Beamslot_lead_depth = fDetCon->fHArmBuilder->f48D48depth;

  G4Box *Beamslot_lead_box = new G4Box("Beamslot_lead_box", Beamslot_lead_width/2.0, Beamslot_lead_height/2.0, Beamslot_lead_depth/2.0 );
  //G4LogicalVolume *Beamslot_lead_log = new G4LogicalVolume( 

  G4Box *Beamslot_cut_box = new G4Box( "Beamslot_cut_box", 50.0*cm/2.0, Beamslot_lead_height/2.0 + mm, Beamslot_notch_depth );

  G4RotationMatrix *notch_rot = new G4RotationMatrix;
  notch_rot->rotateY( -45.0*deg );

  G4double notch_angle = 45.0*deg;
  
  G4SubtractionSolid *Beamslot_lead_box_cut = new G4SubtractionSolid( "Beamslot_lead_box_cut", Beamslot_lead_box, Beamslot_cut_box, notch_rot,
								      G4ThreeVector( -Beamslot_lead_width/2.0,
										     0.0,
										     -Beamslot_lead_depth/2.0 ) );
  
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



