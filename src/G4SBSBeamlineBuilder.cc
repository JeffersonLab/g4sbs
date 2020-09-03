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
#include "G4SBSCalSD.hh"

#include "G4SBSBDParameterisation.hh"
#include "G4SBSBeamDiffuserSD.hh"

G4SBSBeamlineBuilder::G4SBSBeamlineBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
  assert(dc);
}

G4SBSBeamlineBuilder::~G4SBSBeamlineBuilder(){;}

void G4SBSBeamlineBuilder::BuildComponent(G4LogicalVolume *worldlog){
  
  double beamheight = 10.0*12*2.54*cm; // 10 feet off the ground
 
  G4bool bdEnable = fDetCon->GetBeamDiffuserEnable();  // is the beam diffuser enabled? 
 
  // EFuchey 2017/03/29: organized better this with a switch instead of an endless chain of if...  else...
  //if( (fDetCon->fTargType == kLH2 || fDetCon->fTargType == kLD2) ){
  switch(fDetCon->fExpType){
  case(G4SBS::kGEp):
  case(G4SBS::kGEPpositron):
    printf("GEp experiment: forcing beamline configuration 1 \n");
    fDetCon->fBeamlineConf = 1;
    MakeGEpBeamline(worldlog);
    if(fDetCon->fLeadOption == 1){
      MakeGEpLead(worldlog);
    }
    break;
  case(G4SBS::kGMN):// GMn
    MakeGMnBeamline(worldlog);
    if(fDetCon->fLeadOption == 1){
      MakeGMnLead(worldlog);
    }
    break;
  case(G4SBS::kGEnRP):// GEnRP
    MakeGMnBeamline(worldlog);
    if(fDetCon->fLeadOption == 1){
      MakeGMnLead(worldlog);
      MakeGEnRPLead(worldlog);
    }
    break;
  case(G4SBS::kGEN):// GEn
    printf("GEn experiment: forcing beamline configuration 2 \n");
    fDetCon->fBeamlineConf = 2;
    Make3HeBeamline(worldlog);
    MakeGEnClamp(worldlog);
    if(bdEnable) MakeBeamDiffuser(worldlog);
    if(fDetCon->fLeadOption == 1){
      MakeGEnLead(worldlog);
    }
    break;
  case(G4SBS::kSIDISExp):// SIDIS
    Make3HeBeamline(worldlog);
    if(fDetCon->fLeadOption == 1){
      MakeSIDISLead(worldlog);
    }
    break;
  case(G4SBS::kGEMHCtest):// Hall C GEM test
    MakeGMnBeamline(worldlog);
    break;  
  default:
    MakeDefaultBeamline(worldlog);
    break;
  }
  
  double floorthick = 1.0*m;
  G4Tubs *floor_tube = new G4Tubs("floor_tube", 0.0, 30*m, floorthick/2, 0.*deg, 360.*deg );

  G4RotationMatrix *floorrm = new G4RotationMatrix;
  floorrm->rotateX(90*deg);

  G4LogicalVolume *floorLog = new G4LogicalVolume(floor_tube, GetMaterial("Concrete"), "floor_log", 0, 0, 0);
  new G4PVPlacement(floorrm, G4ThreeVector(0.0, -floorthick/2 - beamheight, 0.0), floorLog, "floor_phys", worldlog, false, 0);

  // G4VisAttributes *floorVisAtt= new G4VisAttributes(G4Colour(0.9,0.1,0.9));
  // 	floorLog->SetVisAttributes(floorVisAtt); 
  
  floorLog->SetVisAttributes(G4VisAttributes::Invisible);

  /*
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
  */
  return;

}

/*
// Entrance Beam line construction
void G4SBSBeamlineBuilder::MakeEntranceBeamline(G4LogicalVolume *worldlog){
  G4VisAttributes *Vacuum_visatt = new G4VisAttributes(G4Colour(0.1, 0.5, 0.9 ) );
  //Vacuum_visatt->SetForceSolid(true);
  Vacuum_visatt->SetVisibility(false);
  G4VisAttributes *SteelColor = new G4VisAttributes( G4Colour( 0.75, 0.75, 0.75 ) );
  
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
}
*/

void G4SBSBeamlineBuilder::MakeCommonExitBeamline(G4LogicalVolume *worldlog) {
  bool ChkOverlaps = false;
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
  
  G4double z_formed_bellows = 52.440*inch - TargetCenter_zoffset; //relative to "target center"? or "origin"?
  G4double z_spool_piece = 58.44*inch - TargetCenter_zoffset;
  if(fDetCon->fBeamlineConf>2)z_spool_piece = 27.903*inch;
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
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange1_vac_log, "CVLW_Flange1_vac_phys", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange1_log, "CVLW_Flange1_phys", worldlog, false, 0 , ChkOverlaps );

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
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_vac_log, "CVLW_vac_phys", worldlog, false, 0 , ChkOverlaps);
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_log, "CVLW_phys", worldlog, false, 0 , ChkOverlaps );

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
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange2_vac_log, "CVLW_Flange2_vac_phys", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange2_log, "CVLW_Flange2_phys", worldlog, false, 0 , ChkOverlaps );

  //Next: "Welded bellows"
  G4double dz_welded_bellows = 207.144*inch - z_welded_bellows - TargetCenter_zoffset; // = =5.512 inches
  
  Rin = 11.750/2.0*inch;
  Rout = 14.0/2.0*inch;
  Thick = 1.12*inch;
  
  G4Tubs *WB_Flange = new G4Tubs( "WB_Flange", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *WB_Flange_log = new G4LogicalVolume( WB_Flange, GetMaterial("Stainless_Steel"), "WB_Flange_log" );

  WB_Flange_log->SetVisAttributes( SteelColor );
    
  Z = z_welded_bellows + Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Flange_log, "WB_Flange1_phys", worldlog, false, 0 , ChkOverlaps );

  Z = z_welded_bellows + dz_welded_bellows - Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Flange_log, "WB_Flange2_phys", worldlog, false, 1 , ChkOverlaps );
  
  Rout = Rin + 0.125*inch;
  Thick = dz_welded_bellows - 2*1.12*inch;
  G4Tubs *WB_Bellows = new G4Tubs( "WB_Bellows", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *WB_Bellows_log = new G4LogicalVolume(WB_Bellows, GetMaterial("Stainless_Steel"), "WB_Bellows_log" );

  WB_Bellows_log->SetVisAttributes( SteelColor );
  
  Z = z_welded_bellows + 1.12*inch + Thick/2.0;
  
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Bellows_log, "WB_Bellows_phys", worldlog, false, 0 , ChkOverlaps );
  
  Rin = 0.0;
  Rout = 11.750/2.0*inch;
  Thick = dz_welded_bellows;
  G4Tubs *WB_Vacuum = new G4Tubs( "WB_Vacuum", Rin, Rout, Thick/2.0, 0.0, twopi );

  G4LogicalVolume *WB_Vacuum_log = new G4LogicalVolume(WB_Vacuum, GetMaterial("Vacuum"), "WB_Vacuum_log" );

  WB_Vacuum_log->SetVisAttributes( Vacuum_visatt );
  
  Z = z_welded_bellows + dz_welded_bellows/2.0;

  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Vacuum_log, "WB_Vacuum_phys", worldlog, false, 0 , ChkOverlaps );
  
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
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TVL9_log, "Extended_Vac1", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TBL9_log, "Extended_Al1", worldlog, false, 0 , ChkOverlaps );

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
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TMV9_log, "Extended_Vac2", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TML9_log, "Extended_Al2", worldlog, false, 0 , ChkOverlaps );

  // For CPU speed, extend vacuum all the way to the edge of the "world" volume, so that we don't track beam electrons in air beyond interesting region.
  G4double Zstop  = 30.0*m;
  G4double Zstart = Z + tDzz;

  G4double Zwidth = (Zstop-Zstart);
  G4Tubs *FakeVacuumExtension = new G4Tubs( "FakeVacuumExtension", tRmin, tRmax, Zwidth/2.0, tSPhi, tDphi );
  G4LogicalVolume *FakeVacuumExtension_log = new G4LogicalVolume( FakeVacuumExtension, GetMaterial("Vacuum"), "FakeVacuumExtension_log" );
  FakeVacuumExtension_log->SetVisAttributes( Vacuum_visatt );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0.5*(Zstop+Zstart)), FakeVacuumExtension_log, "FakeVacuumExtension_phys", worldlog,false,0 , ChkOverlaps);

  //-----------------------------------------------------
  //       magnetic tubes
    
  // Inner and outer Magnetic shieldings: 
  // arrays of variables to parameterize the different geometries with the beamline config flag
  G4int Ndivs = 1;// number of segments with shielding
  //G4double Rin_array[6];// radii for inner shielding elements
  //G4double Zin_array[6];// z for inner shielding elements
  //G4int Nrings_out[3];// number of outer elements per segments
  //G4double Rout_array[6];// radii for outer shielding elements
  //G4double Zout_array[6];// z for outer shielding elements
  std::vector<G4double> Rin_array;// radii for inner shielding elements
  std::vector<G4double> Zin_array;// z for inner shielding elements
  std::vector<G4int> Nrings_out;// number of outer elements per segments
  std::vector<G4double> Rout_array;// radii for inner shielding elements
  std::vector<G4double> Zout_array;// z for inner shielding elements

  G4double OMthick = 1.625*inch;
  G4double OMspace = 0.375*inch;

  switch(fDetCon->fBeamlineConf){
  case(1):// reminder: beamline config 1 = GEp
    Ndivs = 2;
    
    Rin_array.push_back( 4.354*inch/2.0 );
    Rin_array.push_back( 6.848*inch/2.0 );
    Rin_array.push_back( 8.230*inch/2.0 );
    Rin_array.push_back( 10.971*inch/2.0 );
    
    Zin_array.push_back( z_inner_magnetic );
    Zin_array.push_back( z_inner_magnetic + 47.625*inch );
    Zin_array.push_back( z_inner_magnetic + 74.00*inch );
    Zin_array.push_back( z_inner_magnetic + (74.00 + 52.347)*inch );
    
    Nrings_out.push_back( 26 );
    Nrings_out.push_back( 27 );
    
    Rout_array.push_back( 5.5*inch/2.0 );
    Rout_array.push_back( 8.178*inch/2.0 );
    Rout_array.push_back( 9.349*inch/2.0 );
    Rout_array.push_back( 12.156*inch/2.0 );
    
    Zout_array.push_back( z_outer_magnetic );
    Zout_array.push_back( z_outer_magnetic + 26.0*OMthick + 25.0*OMspace );
    Zout_array.push_back( z_outer_magnetic + 74.0*inch );
    Zout_array.push_back( z_outer_magnetic + 74.0*inch + 27.0*OMthick + 26.0*OMspace );
    
    //G4double Rin_array_tmp1[6] = {4.354*inch/2.0, 6.848*inch/2.0, 8.230*inch/2.0, 10.971*inch/2.0, 0.0, 0.0};
    //G4double Zin_array_tmp1[6] = {z_inner_magnetic, z_inner_magnetic + 47.625*inch, z_inner_magnetic + 74.00*inch, z_inner_magnetic + (74.00 + 52.347)*inch, 0.0, 0.0};
    //G4int Nrings_out_tmp1[3] = {26, 27, 0};
    //G4double Rout_array_tmp1[6] = {5.5*inch/2.0, 8.178*inch/2.0, 9.349*inch/2.0, 12.156*inch/2.0, 0.0, 0.0};
    //G4double Zout_array_tmp1[6] = {z_outer_magnetic, z_outer_magnetic + 26.0*OMthick + 25.0*OMspace, z_outer_magnetic + 74.0*inch, z_outer_magnetic + 74.0*inch + 27.0*OMthick + 26.0*OMspace, 0.0, 0.0};
    break;
  // case(2):// reminder: beamline config 2 = GEn, SIDIS
  // Ndivs = 3;
  //   break;
  case(3):// reminder: beamline config 3 = GMn, all Q^2
    Ndivs = 3;
    
    Rin_array.push_back( 3.7745*inch/2.0 );
    Rin_array.push_back( 4.307*inch/2.0 );
    Rin_array.push_back( 5.264*inch/2.0 );
    Rin_array.push_back( 7.905*inch/2.0 );
    Rin_array.push_back( 9.310*inch/2.0 );
    Rin_array.push_back( 10.930*inch/2.0 );
    
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 +11.62 - 1.625)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 +11.62 + 14.38 + 2.0)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 - 2.0)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 + 22.38 + 2.0)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 138.83 - 1.12 - 1.14)*inch );
    
    Nrings_out.push_back( 6 );
    Nrings_out.push_back( 27 );
    Nrings_out.push_back( 17 );
    
    Rout_array.push_back( (3.7745/2+0.38)*inch );
    Rout_array.push_back( (4.392/2+0.38)*inch );
    Rout_array.push_back( (5.158/2+0.38)*inch );
    Rout_array.push_back( (8.012/2+0.38)*inch );
    Rout_array.push_back( (9.203/2+0.38)*inch );
    Rout_array.push_back( (10.923/2+0.38)*inch );
    
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 + 22.38)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 138.83 - 1.12 - 1.14)*inch );
    
    //G4double Rin_array_tmp3[6] = {3.774*inch/2.0, 4.307*inch/2.0, 5.264*inch/2.0, 7.905*inch/2.0, 9.310*inch/2.0, 10.930*inch/2.0};
    //G4double Zin_array_tmp3[6] = {z_conic_vacline_weldment + (0.84 + 0.14)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62 - 1.625)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 2.0)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 - 2.0)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 + 22.38 + 2.0)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 138.83 - 1.12 - 1.14)*inch};
    //G4int Nrings_out_tmp3[3] = {6, 27, 17};
    //G4double Rout_array_tmp3[6] = {(3.774/2+0.38)*inch, (4.392/2+0.38)*inch, (5.158/2+0.38)*inch, (8.012/2+0.38)*inch, (9.203/2+0.38)*inch, (10.923/2+0.38)*inch};
    //G4double Zout_array_tmp3[6] = {z_conic_vacline_weldment + (0.84 + 0.14)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 + 22.38)*inch, z_conic_vacline_weldment + (0.84 + 138.83 - 1.12 - 1.14)*inch};
    break;
  case(4):// reminder: beamline config 4 = GMn, Q^2 = 13.5 GeV^2 (+calibrations)
    Ndivs = 2;
    
    Rin_array.push_back( 3.7745*inch/2.0 );
    Rin_array.push_back( 6.096*inch/2.0 );
    Rin_array.push_back( 7.074*inch/2.0 );
    Rin_array.push_back( 9.609*inch/2.0 );
    
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 - 2.0)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38 + 2.0)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38 + 51.62 - 2.0)*inch );
    
    Nrings_out.push_back( 23 );
    Nrings_out.push_back( 26 );
    
    Rout_array.push_back( (3.7745/2+0.38)*inch );
    Rout_array.push_back( (6.202/2+0.38)*inch );
    Rout_array.push_back( (6.968/2+0.38)*inch );
    Rout_array.push_back( (9.715/2+0.38)*inch );
    
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14)*inch  );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 45.63)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38 + 51.62)*inch );
    
    //G4double Rin_array_tmp4[6] = {3.774*inch/2.0, 6.096*inch/2.0, 7.074*inch/2.0, 9.609*inch/2.0, 0.0, 0.0};
    //G4double Zin_array_tmp4[6] = {z_conic_vacline_weldment + (0.84 + 0.14)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 - 2.0)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38 + 2.0)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38 + 51.62 - 2.0)*inch, 0.0, 0.0};
    //G4int Nrings_out_tmp4[3] = {23, 26, 0};
    //G4double Rout_array_tmp4[6] = {(3.774/2+0.38)*inch, (6.202/2+0.38)*inch, (6.968/2+0.38)*inch, (9.715/2+0.38)*inch, 0.0, 0.0};
    //G4double Zout_array_tmp4[6] = {z_conic_vacline_weldment + (0.84 + 0.14)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 45.63)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38)*inch, z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38 + 51.62)*inch, 0.0, 0.0};    
    break;
  default:// default: all elements built
    Ndivs = 1;
    Rin_array.push_back( 3.7745*inch/2.0 );
    Rin_array.push_back( 10.923*inch/2.0 );

    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 138.83 - 1.12 - 1.14)*inch );

    Nrings_out.push_back( 68 );
    
    Rout_array.push_back( (3.7745/2+0.38)*inch );
    Rout_array.push_back( (10.923/2+0.38)*inch );

    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 138.83 - 1.12 - 1.14)*inch );
    
    //G4double Rin_array_tmp0[6] = {3.774*inch/2.0, 10.923*inch/2.0, 0.0, 0.0, 0.0, 0.0};
    //G4double Zin_array_tmp0[6] = {z_conic_vacline_weldment + (0.84 + 0.14)*inch, z_conic_vacline_weldment + (0.84 + 138.83 - 1.12 - 1.14)*inch, 0.0, 0.0, 0.0, 0.0};
    //G4int Nrings_out_tmp0[3] = {68, 0, 0};
    //G4double Rout_array_tmp0[6] = {(3.774/2+0.38)*inch, (10.923/2+0.38)*inch, 0.0, 0.0, 0.0, 0.0};
    //G4double Zout_array_tmp0[6] = {z_conic_vacline_weldment + (0.84 + 0.14)*inch, z_conic_vacline_weldment + (0.84 + 138.83 - 1.12 - 1.14)*inch, 0.0, 0.0, 0.0, 0.0};
    break;
  }

  // Building beamline shielding:
  for(G4int i = 0; i<Ndivs; i++){
    // Building beamline shielding: inner elements
    Rin1 = Rin_array[2*i];
    Rout1 = Rin1 + 0.25*inch;
    Rin2 = Rin_array[2*i+1];
    Rout2 = Rin2 + 0.25*inch;
    Thick = Zin_array[2*i+1]-Zin_array[2*i];
    
    char cname[100];
    sprintf(cname,"IM_%d", i);
    G4String name = cname;
     
    G4Cons *IM_ = new G4Cons( name, Rin1, Rout1, Rin2, Rout2, Thick/2.0, 0.0, twopi );
    name += "_log";
    G4LogicalVolume *IM__log = new G4LogicalVolume( IM_, GetMaterial("Iron"), name );
    
    IM__log->SetVisAttributes( ironColor );
        
    Z = (Zin_array[2*i]+Zin_array[2*i+1])/2.0;
    
    name = cname;
    name += "_phys";
    new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), IM__log, name, worldlog, false, 0 , ChkOverlaps );
    
    // Building beamline shielding: outer elements
    G4double zmin = Zout_array[2*i];
    G4double zmax = Zout_array[2*i+1];
    
    G4double Rin_min = Rout_array[2*i];
    G4double Rin_max = Rout_array[2*i+1];
    for( G4int j=0; j<Nrings_out[i]; j++ ){
      char cname[100];
      sprintf(cname,"OM_%d_ring%d", i, j);
      G4String name = cname;
      
      G4double zstart = zmin + j*(OMthick + OMspace);
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
      new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), ring_log, name, worldlog, false, 0  , ChkOverlaps);
    }
  }
  
  if(fDetCon->fBeamlineConf!=2){
    G4double dz_spool_piece = z_conic_vacline_weldment - z_spool_piece;
    
    //Make Spool piece vacuum:
    Rout = 3.76*inch/2.0;
    Rin = 0.0;
    Thick = dz_spool_piece;
    
    G4Tubs *SpoolPiece_vac = new G4Tubs("SpoolPiece_vac", Rin, Rout, Thick/2.0, 0.0, twopi );
    G4LogicalVolume *SpoolPiece_vac_log = new G4LogicalVolume( SpoolPiece_vac, GetMaterial("Vacuum"), "SpoolPiece_vac_log" );
    
    SpoolPiece_vac_log->SetVisAttributes( Vacuum_visatt );
    
    Z = z_spool_piece + Thick/2.0;
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), SpoolPiece_vac_log, "SpoolPiece_vac_phys", worldlog, false, 0 , ChkOverlaps );
    
    Rin = 3.76*inch/2.0;
    Rout = 6.00*inch/2.0;
    Thick = 0.84*inch;
    
    G4Tubs *SpoolPiece_Flange1 = new G4Tubs("SpoolPiece_Flange1", Rin, Rout, Thick/2.0, 0.0, twopi );
    G4LogicalVolume *SpoolPiece_Flange1_log = new G4LogicalVolume( SpoolPiece_Flange1, GetMaterial("Stainless_Steel"), "SpoolPiece_Flange1_log" );
    
    SpoolPiece_Flange1_log->SetVisAttributes( SteelColor );
    
    Z = z_spool_piece + Thick/2.0;
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), SpoolPiece_Flange1_log, "SpoolPiece_Flange1_phys", worldlog, false, 0 , ChkOverlaps );
    
    Rout = 6.75*inch/2.0;
    Thick = 0.84*inch;
    
    G4Tubs *SpoolPiece_Flange2 = new G4Tubs("SpoolPiece_Flange2", Rin, Rout, Thick/2.0, 0.0, twopi );
    G4LogicalVolume *SpoolPiece_Flange2_log = new G4LogicalVolume( SpoolPiece_Flange2, GetMaterial("Stainless_Steel"), "SpoolPiece_Flange2_log" );
    
    SpoolPiece_Flange2_log->SetVisAttributes( SteelColor );
    
    Z = z_conic_vacline_weldment - Thick/2.0;
    
    new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), SpoolPiece_Flange2_log, "SpoolPiece_Flange2_phys", worldlog, false, 0 , ChkOverlaps );
    
    Rout = 4.0*inch/2.0;
    Thick = dz_spool_piece - 2.0*0.84*inch;
    
    G4Tubs *SpoolPiece_tube = new G4Tubs("SpoolPiece_tube", Rin, Rout, Thick/2.0, 0.0, twopi );
    
    G4LogicalVolume *SpoolPiece_tube_log = new G4LogicalVolume( SpoolPiece_tube, GetMaterial("Stainless_Steel"), "SpoolPiece_tube_log" );
    
    SpoolPiece_tube_log->SetVisAttributes( SteelColor );
    
    Z = z_spool_piece + dz_spool_piece/2.0;
    
    new G4PVPlacement( 0,  G4ThreeVector( X, Y, Z ), SpoolPiece_tube_log, "SpoolPiece_tube_phys", worldlog, false, 0 , ChkOverlaps );
    
    if(fDetCon->fBeamlineConf>1){
      Rin1 = 4.00*inch/2.0+0.02*inch;
      Rout1 = Rin1 + 0.25*inch;
      //Thick = dz_spool_piece - 2.0*0.84*inch;
      
      G4Tubs *IM0 = new G4Tubs( "IM0", Rin1, Rout1, Thick/2.0, 0.0, twopi );
      G4LogicalVolume *IM0_log = new G4LogicalVolume( IM0, GetMaterial("Iron"), "IM0_log" );
      
      IM0_log->SetVisAttributes( ironColor );
      
      //Z = z_spool_piece + dz_spool_piece/2.0;
      new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), IM0_log, "IM0_phys", worldlog, false, 0 , ChkOverlaps );
      
      G4double zmin = z_spool_piece+0.84*inch+0.44*inch;
      G4double zmax = zmin + 13.0*OMthick + 12.0*OMspace;
  
      G4double Rin_min = 5.3*inch/2.0;
      for( G4int i=0; i<13; i++ ){
	char cname[100];
	sprintf(cname,"OM0_ring%d", i);
	G4String name = cname;
	
	G4double zstart = zmin + i*(OMthick + OMspace);
	G4double zstop = zstart + OMthick;
	
	G4Tubs *ring = new G4Tubs( name,Rin_min, Rin_min+0.5*inch, OMthick/2.0, 0.0, twopi );
	
	name += "_log";
	G4LogicalVolume *ring_log = new G4LogicalVolume( ring, GetMaterial("Iron"), name );
	
	ring_log->SetVisAttributes( ironColor );
	
	name = cname;
	name += "_phys";
	
	Z = 0.5*(zstart + zstop);
	new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), ring_log, name, worldlog, false, 0 , ChkOverlaps );
	
      }
      
    }
  }
    
  if(fDetCon->fBeamlineConf==1){
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
    
    new G4PVPlacement(  0,  G4ThreeVector( X, Y, Z ), FormedBellows_vac_log, "FormedBellows_vac_phys", worldlog, false, 0 , ChkOverlaps );
    
    Rin = 3.81*inch/2.0;
    Rout = 6.00*inch/2.0;
    Thick = 0.84*inch;
    
    //Flanges for formed bellows:
    G4Tubs *FormedBellows_Flange = new G4Tubs("FormedBellows_Flange", Rin, Rout, Thick/2.0, 0.0, twopi );
    G4LogicalVolume *FormedBellows_Flange_log = new G4LogicalVolume( FormedBellows_Flange, GetMaterial("Stainless_Steel"), "FormedBellows_Flange_log" );

    FormedBellows_Flange_log->SetVisAttributes( SteelColor );
    
    Z = z_formed_bellows + Thick/2.0;
    
    new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), FormedBellows_Flange_log, "FormedBellows_Flange1_phys", worldlog, false, 0 , ChkOverlaps );
    
    Z = z_formed_bellows + dz_formed_bellows - Thick/2.0;
    
    new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), FormedBellows_Flange_log, "FormedBellows_Flange2_phys", worldlog, false, 1 , ChkOverlaps );
    
    //Tube for formed bellows:
    
    Rout = Rin + 0.125*inch; //This is just a guess!!
    Thick = dz_formed_bellows - 2.0*0.84*inch;
    
    G4Tubs *FormedBellows_tube = new G4Tubs( "FormedBellows_tube", Rin, Rout, Thick/2.0, 0.0, twopi );
    G4LogicalVolume *FormedBellows_tube_log = new G4LogicalVolume( FormedBellows_tube, GetMaterial("Stainless_Steel"), "FormedBellows_tube_log" );

    FormedBellows_tube_log->SetVisAttributes( SteelColor );
  
    Z = z_formed_bellows + dz_formed_bellows/2.0;

    new G4PVPlacement( 0, G4ThreeVector( X, Y, Z ), FormedBellows_tube_log, "FormedBellows_tube_phys", worldlog, false, 0 , ChkOverlaps );

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

    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), IronTube1_log, "IronTube1_phys", worldlog, false, 0 , ChkOverlaps );
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), IronTube1_vac_log, "IronTube1_vac_phys", worldlog, false, 0 , ChkOverlaps );

    Rin = 2.415*inch;
    Rout = 2.5*inch;
  
    G4Tubs *IronTube2 = new G4Tubs("IronTube2", Rin, Rout, Thick/2.0, 0.0, twopi );
    G4Tubs *IronTube2_vac = new G4Tubs("IronTube2_vac", 0.0, Rin, Thick/2.0, 0.0, twopi );

    G4LogicalVolume *IronTube2_log = new G4LogicalVolume( IronTube2, GetMaterial("Iron"), "IronTube2_log" );
    G4LogicalVolume *IronTube2_vac_log = new G4LogicalVolume( IronTube2_vac, GetMaterial("Vacuum"), "IronTube2_vac_log" );

    IronTube2_log->SetVisAttributes(ironColor);
    IronTube2_vac_log->SetVisAttributes(Vacuum_visatt);

    Z += Thick;

    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), IronTube2_log, "IronTube2_phys", worldlog, false, 0 , ChkOverlaps );
    new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), IronTube2_vac_log, "IronTube2_vac_phys", worldlog, false, 0 , ChkOverlaps );
  }
  
  //Next, corrector magnets:
  // Define some dimensions that are going to be useful to define the distances
  G4double UpstreamCoilThickY = 1.68*inch;
  G4double UpstreamCoilThickX = 3.46*inch;
  //G4double UpstreamCoilWidth = 3.46*inch;
  G4double UpstreamCoilHeight = 8.17*inch;
  G4double UpstreamCoilDepth = 6.60*inch;
  G4double UpstreamCoilWidth = 7.56*inch;
    
  G4double YokeTopPiece_Width = 15.04*inch;
  G4double YokeTopPiece_Height = 3.94*inch;
  G4double YokeTopPiece_Depth = 6.30*inch;
  
  G4double YokeLeftPiece_Width = 2.76*inch;
  G4double YokeLeftPiece_Height = 11.81*inch + 2.0*YokeTopPiece_Height;
  G4double YokeLeftPiece_Depth = YokeTopPiece_Depth;
  
  G4double YokeRightNotchAngle = 18.43*deg;
  G4double YokeRightWidthFinal = YokeLeftPiece_Width;
  G4double YokeRightZFinal = YokeLeftPiece_Depth - 0.39*inch;
  G4double YokeRightWidthInitial = YokeRightWidthFinal - YokeRightZFinal*tan(YokeRightNotchAngle );
  
  G4double DownstreamTotalWidth = 17.58*inch + 2.0*2.76*inch;
  G4double DownstreamTotalHeight = 20.16*inch + 2.0*3.94*inch;
  G4double DownstreamYokeDepth = 15.75*inch;

  G4double DS_coil_depth = 8.91*inch;
  G4double DS_coil_height = 12.04*inch;
  G4double DS_coil_ThickX = 2.90*inch;
  G4double DS_coil_ThickY = 1.68*inch;
  
  // Z Array to change easily z values with beamline configuration;
  // Right now, it looks like X and Y do NOT need to change depending on the configuration; only Z does
  std::vector<G4double> z_Magnets_array;
  switch(fDetCon->fBeamlineConf){
  case(1):// reminder: beamline config 1 = GEp
    z_Magnets_array.push_back( z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY );
    z_Magnets_array.push_back( z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0 );
    z_Magnets_array.push_back( z_formed_bellows + 76.09*inch + 1.71*inch + DownstreamYokeDepth/2.0 );
    z_Magnets_array.push_back( z_formed_bellows + 76.09*inch + DS_coil_ThickY + DS_coil_depth/2.0 );
    
    //z_Magnets_array = {z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY, z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0, z_formed_bellows + 76.09*inch + 1.71*inch + DownstreamYokeDepth/2.0, z_formed_bellows + 76.09*inch + DS_coil_ThickY + DS_coil_depth/2.0};
     break;
  // case(2):// reminder: beamline config 2 = GEn, SIDIS ?
  //   Z = z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY;
  //   break;
  case(3):// reminder: beamline config 3 = GMn 
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 15.94)*inch + UpstreamCoilDepth/2.0 );
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 15.94 + 8.3 - 6.47)*inch - UpstreamCoilThickY + YokeRightZFinal/2.0 );
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 85.78)*inch + DownstreamYokeDepth/2.0 );
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 85.78 - 1.71)*inch + DS_coil_ThickY + DS_coil_depth/2.0 );
    
    //z_Magnets_array = {z_conic_vacline_weldment + (0.84 + 0.14 + 15.94)*inch + UpstreamCoilDepth/2.0, z_conic_vacline_weldment + (0.84 + 0.14 + 15.94 + 8.3 - 6.47)*inch - UpstreamCoilThickY + YokeRightZFinal/2.0, z_conic_vacline_weldment + (0.84 + 0.14 + 85.78)*inch + DownstreamYokeDepth/2.0, z_conic_vacline_weldment + (0.84 + 0.14 + 85.78 - 1.71)*inch + DS_coil_ThickY + DS_coil_depth/2.0};    
    break;
  case(4):// reminder: beamline config 3 = GMn, Q^2 = 13.5 GeV^2
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 49.47)*inch + UpstreamCoilDepth/2.0 );
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 49.47 + 8.3 - 6.47)*inch - UpstreamCoilThickY + YokeRightZFinal/2.0  );
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 118.34)*inch + DownstreamYokeDepth/2.0 );
    z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 118.34 - 1.71)*inch + DS_coil_ThickY + DS_coil_depth/2.0 );
    
    //z_Magnets_array = {z_conic_vacline_weldment + (0.84 + 0.14 + 49.47)*inch + UpstreamCoilDepth/2.0, z_conic_vacline_weldment + (0.84 + 0.14 + 49.47 + 8.3 - 6.47)*inch - UpstreamCoilThickY + YokeRightZFinal/2.0, z_conic_vacline_weldment + (0.84 + 0.14 + 118.34)*inch + DownstreamYokeDepth/2.0, z_conic_vacline_weldment + (0.84 + 0.14 + 118.34 - 1.71)*inch + DS_coil_ThickY + DS_coil_depth/2.0};
    break;
  default:
    //sent back: who cares...
    z_Magnets_array.push_back( -10.0*m + 0.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY );
    z_Magnets_array.push_back( -10.0*m + 2.3*inch + YokeRightZFinal/2.0 );
    z_Magnets_array.push_back( -10.0*m + 15.24*inch + 1.71*inch + DownstreamYokeDepth/2.0 );
    z_Magnets_array.push_back( -10.0*m + 15.24*inch + DS_coil_ThickY + DS_coil_depth/2.0 );
    //z_Magnets_array = {z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY, z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0, z_formed_bellows + 76.09*inch + 1.71*inch + DownstreamYokeDepth/2.0, z_formed_bellows + 76.09*inch + DS_coil_ThickY + DS_coil_depth/2.0};
    break;
  }
  
  G4Box *UpstreamCoil_outer = new G4Box("UpstreamCoil_outer", UpstreamCoilThickX/2.0, (UpstreamCoilHeight+2.0*UpstreamCoilThickY)/2.0, (UpstreamCoilDepth + 2.0*UpstreamCoilThickY)/2.0 );
  G4Box *UpstreamCoil_inner = new G4Box("UpstreamCoil_inner", UpstreamCoilThickX/2.0 + cm, UpstreamCoilHeight/2.0, UpstreamCoilDepth/2.0 );

  G4SubtractionSolid *UpstreamCoil = new G4SubtractionSolid( "UpstreamCoil", UpstreamCoil_outer, UpstreamCoil_inner, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *UpstreamCoil_log = new G4LogicalVolume(UpstreamCoil, GetMaterial("Copper"), "UpstreamCoil_log" );

  UpstreamCoil_log->SetVisAttributes( CopperColor );

  Z = z_Magnets_array[0];//z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY;
  X = (UpstreamCoilWidth+UpstreamCoilThickX)/2.0;
  Y = 0.0;

  //two placements of upstream coil:
  
  new G4PVPlacement( 0, G4ThreeVector(-X, Y, Z), UpstreamCoil_log, "UpstreamCoil_phys_right", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), UpstreamCoil_log, "UpstreamCoil_phys_left", worldlog, false, 1  , ChkOverlaps);

  G4double UpstreamPoleDepth = 6.3*inch;
  G4double UpstreamPoleWidth = 4.02*inch;
  G4double UpstreamPoleHeight = 7.87*inch;
  //Next, make poles:
  G4Box *UpstreamPole = new G4Box( "UpstreamPole", UpstreamPoleWidth/2.0, UpstreamPoleHeight/2.0, UpstreamPoleDepth/2.0 );
  G4LogicalVolume *UpstreamPole_log = new G4LogicalVolume( UpstreamPole, GetMaterial("Iron"), "UpstreamPole_log" );
  UpstreamPole_log->SetVisAttributes( ironColor );
  //two placements of upstream poles:

  new G4PVPlacement( 0, G4ThreeVector(-X, Y, Z), UpstreamPole_log, "UpstreamPole_phys_right", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), UpstreamPole_log, "UpstreamPole_phys_left", worldlog, false, 1 , ChkOverlaps );

  //Next, make surrounding yoke:
  // G4double YokeTopPiece_Width = 15.04*inch;
  // G4double YokeTopPiece_Height = 3.94*inch;
  // G4double YokeTopPiece_Depth = 6.30*inch;

  G4Box *YokeTopPiece = new G4Box("YokeTopPiece", YokeTopPiece_Width/2.0, YokeTopPiece_Height/2.0, YokeTopPiece_Depth/2.0 );
  G4LogicalVolume *YokeTopPiece_log = new G4LogicalVolume( YokeTopPiece, GetMaterial("Iron"), "YokeTopPiece_log" );

  YokeTopPiece_log->SetVisAttributes( ironColor );
  
  X = 0.0;
  Y = (11.81*inch + YokeTopPiece_Height)/2.0;

  //two placements of yoke top piece (top and bottom symmetric):
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeTopPiece_log, "UpstreamYokeTop_phys", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X,-Y,Z), YokeTopPiece_log, "UpstreamYokeBottom_phys", worldlog, false, 1 , ChkOverlaps );

  // G4double YokeLeftPiece_Width = 2.76*inch;
  // G4double YokeLeftPiece_Height = 11.81*inch + 2.0*YokeTopPiece_Height;
  // G4double YokeLeftPiece_Depth = YokeTopPiece_Depth;
  
  G4Box *YokeLeftPiece = new G4Box("YokeLeftPiece", YokeLeftPiece_Width/2.0, YokeLeftPiece_Height/2.0, YokeLeftPiece_Depth/2.0 );
  G4LogicalVolume *YokeLeftPiece_log = new G4LogicalVolume( YokeLeftPiece, GetMaterial("Iron"), "YokeLeftPiece_log" );
  YokeLeftPiece_log->SetVisAttributes(ironColor );
  
  X = 7.52*inch + YokeLeftPiece_Width/2.0;
  Y = 0.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeLeftPiece_log, "UpstreamYokeLeftPiece_phys", worldlog, false, 0 , ChkOverlaps );

  // G4double YokeRightNotchAngle = 18.43*deg;
  // G4double YokeRightWidthFinal = YokeLeftPiece_Width;
  // G4double YokeRightZFinal = YokeLeftPiece_Depth - 0.39*inch;
  // G4double YokeRightWidthInitial = YokeRightWidthFinal - YokeRightZFinal*tan(YokeRightNotchAngle );

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
  Z = z_Magnets_array[1];//z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeRightPiece_log, "UpstreamYokeRightPiece_phys", worldlog, false, 0 , ChkOverlaps );

  //Downstream Corrector:
  // G4double DownstreamTotalWidth = 17.58*inch + 2.0*2.76*inch;
  // G4double DownstreamTotalHeight = 20.16*inch + 2.0*3.94*inch;
  // G4double DownstreamYokeDepth = 15.75*inch;

  G4double DownstreamYokeGapWidth = 17.58*inch;
  G4double DownstreamYokeGapHeight = 20.16*inch;
  G4Box *DownstreamYoke_box = new G4Box("DownstreamYoke_box", DownstreamTotalWidth/2.0, DownstreamTotalHeight/2.0, DownstreamYokeDepth/2.0 );
  G4Box *DownstreamYoke_gap = new G4Box("DownstreamYoke_gap", DownstreamYokeGapWidth/2.0, DownstreamYokeGapHeight/2.0, DownstreamYokeDepth/2.0+cm );
  G4SubtractionSolid *DownstreamYoke = new G4SubtractionSolid( "DownstreamYoke", DownstreamYoke_box, DownstreamYoke_gap, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *DownstreamYoke_log = new G4LogicalVolume( DownstreamYoke, GetMaterial("Iron"), "DownstreamYoke_log" );

  DownstreamYoke_log->SetVisAttributes( ironColor );

  X = 0.0; Y = 0.0;
  Z = z_Magnets_array[2];//z_formed_bellows + 76.09*inch + 1.71*inch + DownstreamYokeDepth/2.0;

  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), DownstreamYoke_log, "DownstreamYoke_phys", worldlog, false, 0 , ChkOverlaps );

  // G4double DS_coil_depth = 8.91*inch;
  // G4double DS_coil_height = 12.04*inch;
  // G4double DS_coil_ThickX = 2.90*inch;
  // G4double DS_coil_ThickY = 1.68*inch;
  
  G4Box *DS_coil_outer = new G4Box( "DS_coil_outer", DS_coil_ThickX/2.0, (DS_coil_height + 2.0*DS_coil_ThickY)/2.0, (DS_coil_depth + 2.0*DS_coil_ThickY)/2.0 );
  G4Box *DS_coil_inner = new G4Box( "DS_coil_inner", DS_coil_ThickX/2.0+cm, DS_coil_height/2.0, DS_coil_depth/2.0 );

  G4SubtractionSolid *DS_coil = new G4SubtractionSolid( "DS_coil", DS_coil_outer, DS_coil_inner, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *DS_coil_log = new G4LogicalVolume( DS_coil, GetMaterial("Copper"), "DS_coil_log" );
  DS_coil_log->SetVisAttributes(CopperColor );
  
  X = 11.67*inch/2.0 + DS_coil_ThickX/2.0;
  Y = 0.0;
  Z = z_Magnets_array[3];//z_formed_bellows + 76.09*inch + DS_coil_ThickY + DS_coil_depth/2.0;
  
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), DS_coil_log, "DS_coil_phys_left", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(-X,Y,Z), DS_coil_log, "DS_coil_phys_right", worldlog, false, 1 , ChkOverlaps );

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
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), DSpole_log, "DSpole_phys_left", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(-X,Y,Z), DSpole_log, "DSpole_phys_right", worldlog, false, 1 , ChkOverlaps );

  
  if(fDetCon->fBLneutronDet){//TO-DO: set the possibility to deactivate it.
    // EFuchey: 2018/05/29: add a small dummy detector to study the neutron production by the shielding.
    // 
    // double x_blndet[8] = {0.3*m, 0.6*m, 1.6*m, 0.0*m, -0.7*m, 1.0*m, 2.2*m, -2.2*m};
    // double y_blndet[8] = {0.0*m, 0.0*m, 0.0*m, 0.6*m,  0.0*m, 0.0*m, 0.0*m,  0.0*m};
    // double z_blndet[8] = {1.5*m, 2.5*m, 5.0*m, 2.5*m, +0.7*m, 0.0*m, 2.2*m,  2.2*m};
    //
    // G4double ElecX = 2.0*cm;
    // G4double ElecY = 2.0*cm;
    // G4double ElecZ = 2.0*cm;

    double x_blndet = 3.0*m;
    double y_blndet = 0.0*m;
    double z_blndet = 2.5*m;
    
    G4double ElecX = 5.0*cm;
    G4double ElecY = 100.0*cm;
    G4double ElecZ = 100.0*cm;
    
    G4Box *Electronics = new G4Box( "Electronics" , ElecX/2.0, ElecY/2.0, ElecZ/2.0);
    G4LogicalVolume *Electronics_log = new G4LogicalVolume( Electronics , GetMaterial("Silicon"), "Electronics_log" );
    Electronics_log->SetVisAttributes(G4VisAttributes::Invisible);
    G4String GEMElectronicsname = "BLneutronDet";
    G4String  GEMElectronicscollname = "BLneutronDet";
    G4SBSCalSD *GEMElecSD = NULL;
    
    switch(fDetCon->fExpType){
    case(G4SBS::kGEp):
      GEMElectronicsname += "GEp";
      GEMElectronicscollname += "GEp";
      break;
    case(G4SBS::kGMN):// GMn
      //case(kGEN): //
      GEMElectronicsname += "GMn";
      GEMElectronicscollname += "GMn";
      break;
    case(G4SBS::kGEnRP):// GEnRP
      GEMElectronicsname += "GMn";
      GEMElectronicscollname += "GMn";
      break;
    default:
      
      break;
    }
    
    //for(int i_blndet = 0; i_blndet<8; i_blndet++){
    if( !( (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(GEMElectronicsname) )){
      G4cout << "Adding GEM electronics Sensitive Detector to SDman..." << G4endl;
      GEMElecSD = new G4SBSCalSD( GEMElectronicsname, GEMElectronicscollname );
      fDetCon->fSDman->AddNewDetector(GEMElecSD);
      (fDetCon->SDlist).insert(GEMElectronicsname);
      fDetCon->SDtype[GEMElectronicsname] = G4SBS::kCAL;
      (GEMElecSD->detmap).depth = 0;
    }
    Electronics_log->SetSensitiveDetector( GEMElecSD );
    
    if( (fDetCon->StepLimiterList).find( GEMElectronicsname ) != (fDetCon->StepLimiterList).end() ){
      Electronics_log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
    }
    
    // Place the electronics in our hut:
    // new G4PVPlacement( 0, G4ThreeVector(0.0, -ShieldMotherY/2.0 + GPlateY2 + ElecY/2.0, ShieldMotherZ/2.0 - GPlateZ1 - ElecZ/2.0),
    // 		     Electronics_log, "Electronics", ShieldLog, false, 0);
    // new G4PVPlacement( 0, G4ThreeVector(x_blndet[i_blndet], y_blndet[i_blndet], z_blndet[i_blndet]),
    // 		       Electronics_log, "GMn_Electronics", worldlog, false, i_blndet);
    new G4PVPlacement( 0, G4ThreeVector(x_blndet, y_blndet, z_blndet),
		       Electronics_log, "GMn_Electronics", worldlog, false, 0 , ChkOverlaps);
  }
  //}
    
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


// GEp Beamline Construction --- following Sergey's Fortran code
void G4SBSBeamlineBuilder::MakeGEpBeamline(G4LogicalVolume *worldlog) {
  bool ChkOverlaps = false;
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

  new G4PVPlacement( 0, pos_temp, upstream_beampipe_log, "upstream_beampipe_phys", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, pos_temp, upstream_beampipe_vac_log, "upstream_beampipe_vac_phys", worldlog, false, 0 , ChkOverlaps );
  
  MakeCommonExitBeamline(worldlog);
  /*
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
  */
}

// This is the beam line for GMn
void G4SBSBeamlineBuilder::MakeGMnBeamline(G4LogicalVolume *worldlog){
  bool ChkOverlaps = false;
  double swallrad = 1.143*m/2;
  double swallrad_inner = 1.041/2.0*m; 
  //EFuchey: 2017/02/14: change parameters for Standard scat chamber:
  double sc_entbeampipeflange_dist = 25.375*2.54*cm;// entrance pipe flange distance from hall center
  double sc_exbeampipeflange_dist = 27.903*2.54*cm;// exit pipe flange distance from hall center
  
  // Stainless
  G4double ent_len = 10*m;
  //ent_len = ent_len+1.1*m;// for background studies;
  G4double ent_rin = 31.75*mm;
  G4double ent_rou = ent_rin+0.120*mm;
  
  G4Tubs *ent_tube = new G4Tubs("ent_tube", ent_rin, ent_rou, ent_len/2, 0.*deg, 360.*deg );
  G4Tubs *ent_vac  = new G4Tubs("ent_vac", 0.0, ent_rin, ent_len/2, 0.*deg, 360.*deg );
  
  G4LogicalVolume *entLog = new G4LogicalVolume(ent_tube, GetMaterial("Stainless"), "ent_log", 0, 0, 0);
  G4LogicalVolume *entvacLog = new G4LogicalVolume(ent_vac, GetMaterial("Vacuum"), "entvac_log", 0, 0, 0);
  
  //We want to subtract this cylinder from the entry tube/pipe: 
  // NOT for GMn, because the bealine connects with 
  // G4Tubs *cut_cylinder = new G4Tubs("cut_cylinder", 0.0, swallrad, 1.0*m, 0.0*deg, 360.0*deg );
  // G4RotationMatrix *cut_cylinder_rot = new G4RotationMatrix;
  // cut_cylinder_rot->rotateX( -90.0*deg );
  // G4SubtractionSolid *ent_tube_cut = new G4SubtractionSolid( "ent_tube_cut", ent_tube, cut_cylinder, cut_cylinder_rot, G4ThreeVector( 0.0, 0.0, ent_len/2.0 + swallrad_inner ) );
  // G4SubtractionSolid *ent_vac_cut = new G4SubtractionSolid( "ent_vac_cut", ent_vac, cut_cylinder, cut_cylinder_rot, G4ThreeVector( 0.0, 0.0, ent_len/2.0 + swallrad_inner ) );
  // G4LogicalVolume *entLog_cut = new G4LogicalVolume(ent_tube_cut, GetMaterial("Stainless"), "ent_log_cut", 0, 0, 0);
  // G4LogicalVolume *entvacLog_cut = new G4LogicalVolume(ent_vac_cut, GetMaterial("Vacuum"), "entvac_log_cut", 0, 0, 0);
  
  // EFuchey: 2017/02/14
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-sc_entbeampipeflange_dist), entLog, "ent_phys", worldlog, false,0 , ChkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-sc_entbeampipeflange_dist), entvacLog, "entvac_phys", worldlog,false,0 , ChkOverlaps);
  //}
   
  MakeCommonExitBeamline(worldlog);
  
  /*
  // EFuchey: 2017/02/14: add the possibility to change the first parameters for the beam line polycone 
  // Default set of values;
  //double z0 = sc_exbeampipeflange_dist, rin_0 = 6.20*cm, rout_0 = (6.20+0.28*2.54)*cm;
  
  int nsec = 7;
  G4double exit_z[]   = { sc_exbeampipeflange_dist, 592.2*cm, 609.84*cm,609.85*cm, 1161.02*cm, 1161.03*cm,2725.66*cm };// -- Extended beamline for background studies (2016/09/07)

  G4double exit_zero[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  G4double exit_rin[] = { 6.20*cm, 14.8*cm, 15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };// -- Extended beamline for background studies (2016/09/07)
  G4double exit_rou[] = { (6.20+0.28*2.54)*cm, 15.0*cm,15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };// -- Extended beamline for background studies (2016/09/07)
  
  G4Polycone *ext_cone = new G4Polycone("ext_cone", 0.0*deg, 360.0*deg, nsec, exit_z, exit_rin, exit_rou);
  G4Polycone *ext_vac  = new G4Polycone("ext_vac ", 0.0*deg, 360.0*deg, nsec, exit_z, exit_zero, exit_rin);
  
  G4LogicalVolume *extLog = new G4LogicalVolume(ext_cone, GetMaterial("Aluminum"), "ext_log", 0, 0, 0);
  G4LogicalVolume *extvacLog = new G4LogicalVolume(ext_vac, GetMaterial("Vacuum"), "extvac_log", 0, 0, 0);
  
  new G4PVPlacement(0,G4ThreeVector(), extLog, "ext_phys", worldlog, false,0);
  new G4PVPlacement(0,G4ThreeVector(), extvacLog, "extvac_phys", worldlog,false,0);
  
  G4VisAttributes *extVisAtt= new G4VisAttributes(G4Colour(0.9,0.9,0.9));
  extLog->SetVisAttributes(extVisAtt);
  // extvacLog->SetVisAttributes(G4VisAttributes::Invisible);
  // extLog->SetVisAttributes(pipeVisAtt);
  */
  
  G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
  entvacLog->SetVisAttributes(G4VisAttributes::Invisible);
  entLog->SetVisAttributes(pipeVisAtt);
  
  //entvacLog_cut->SetVisAttributes(G4VisAttributes::Invisible);
  //entLog_cut->SetVisAttributes(pipeVisAtt);
}


// This is the beam line for 3He
void G4SBSBeamlineBuilder::Make3HeBeamline(G4LogicalVolume *worldlog){// for GEn, A1n, SIDIS
  double swallrad = 1.143*m/2;
  double swallrad_inner = 1.041/2.0*m; 
  G4double inch = 2.54*cm;
  bool ChkOverlaps = false;
  
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

  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-40.0*cm), entLog, "ent_phys", worldlog, false,0);// -- Extended beamline for background studies (2016/09/07)
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-40.0*cm), entvacLog, "entvac_phys", worldlog,false,0);// -- Extended beamline for background studies (2016/09/07)

  ////Changed 40.0 cm to 13.9 inch from information at Autodesk Viewer

  
  // Add in Be window if no scattering chamber is to be defined:
  if( fDetCon->fTargetBuilder->GetSchamFlag() != 1 ){
    G4double winthick = 0.0127*cm;    
    G4Tubs *ent_win = new G4Tubs("ent_win", 0.0, ent_rin, winthick/2, 0.*deg, 360.*deg );
    G4LogicalVolume *ent_winlog = new G4LogicalVolume(ent_win, GetMaterial("Beryllium"), "entwin_log", 0, 0, 0);
    
    // my cancel Be window for GEp experiment 09/29/2014 
    new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, ent_len/2-winthick/2), ent_winlog, "entwin_phys", entvacLog,false,0);	// => uncommented on 2016/09/07 for background studies
    // Note from  2016/09/07: 
    // Is that normal that this was commented ? My guess would be not.
    
    ent_winlog->SetVisAttributes(new G4VisAttributes(G4Colour(0.7,1.0,0.0)));
  }
  
  // EFuchey: 2017/02/14: add the possibility to change the first parameters for the beam line polycone 
  // // Default set of values;
  // double z0 = 37.2*cm;
  // double rin_0 = 5.64*cm; 
  // double rout_0 = rin_0+(0.28*2.54)*cm;
  
  int nsec = 7;
  //  Definition taken from GEN_10M.opc by Bogdan to z = 5.92.  2mm thickness assumed
  
  /*
  G4double exit_z[]   = { 37.2*cm, 592.2*cm, 609.84*cm,609.85*cm, 1161.02*cm, 1161.03*cm,2725.66*cm };// -- Extended beamline for background studies (2016/09/07)
  G4double exit_zero[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  
  G4double exit_rin[] = { 5.64*cm, 14.8*cm, 15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };// -- Extended beamline for background studies (2016/09/07)
  G4double exit_rou[] = { (5.64+0.28*2.54)*cm, 15.0*cm,15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };// -- Extended beamline for background studies (2016/09/07)
  
  G4Polycone *ext_cone = new G4Polycone("ext_cone", 0.0*deg, 360.0*deg, nsec, exit_z, exit_rin, exit_rou);
  G4Polycone *ext_vac  = new G4Polycone("ext_vac ", 0.0*deg, 360.0*deg, nsec, exit_z, exit_zero, exit_rin);
  
  G4LogicalVolume *extLog = new G4LogicalVolume(ext_cone, GetMaterial("Aluminum"), "ext_log", 0, 0, 0);
  G4LogicalVolume *extvacLog = new G4LogicalVolume(ext_vac, GetMaterial("Vacuum"), "extvac_log", 0, 0, 0);
  
  new G4PVPlacement(0,G4ThreeVector(), extLog, "ext_phys", worldlog, false,0);
  new G4PVPlacement(0,G4ThreeVector(), extvacLog, "extvac_phys", worldlog,false,0);
  
  G4VisAttributes *extVisAtt= new G4VisAttributes(G4Colour(0.9,0.1,0.9));
  extLog->SetVisAttributes(extVisAtt);
    
  double extwin_thick = 5.0e-4*cm;
  
  G4Tubs *extwin = new G4Tubs("ext_win", 0.0, exit_rin[0], extwin_thick/2, 0.*deg, 360.*deg );
  G4LogicalVolume *ext_winlog = new G4LogicalVolume(extwin, GetMaterial("Aluminum"), "extwin_log", 0, 0, 0);
    
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, exit_z[0] - extwin_thick/2), ext_winlog, "extwin_phys", worldlog,false,0);
  
  ext_winlog->SetVisAttributes(new G4VisAttributes(G4Colour(0.5,0.2,0.6)));
  
  extvacLog->SetVisAttributes(G4VisAttributes::Invisible);
  entvacLog->SetVisAttributes(G4VisAttributes::Invisible);
    
  entvacLog_cut->SetVisAttributes(G4VisAttributes::Invisible);
    
  G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
    
  extLog->SetVisAttributes(pipeVisAtt);
  entLog->SetVisAttributes(pipeVisAtt);
  entLog_cut->SetVisAttributes(pipeVisAtt);
  */

  //MakeCommonExitBeamline(worldlog);
  
 






  ////BEGIN EXIT BEAMLINE UPDATE FOR SIDIS - S.SEEDS - MOST RECENT UPDATE: 7.24.20
  ////In progress - change visuals, and verify P1initPlacement_z
  ////Opted to leave out 80/20 rails and related fixtures as a first approximation
  ////Corrected initial placement with help from D. Flay.
  ////Corrected side shields - material to be aluminum per R Wines. 7.24.20

  //Section One
  //This section details all components of the exit beamline from the target chamber to the first cone and shielding including all simple cylinders. All labels numbered by proximity to target chamber.

  //General Specifications
  G4double P1tubeTh = 0.035*inch;
  G4double P1ringL = 0.125/2*inch;
  // G4double initPlacement_z = 15.14*inch-0.2*inch; //Offset from upstream beampipe with cad data and with autodesk viewer
  G4double initPlacement_z = 26.74*inch; //Offset from upstream beampipe with cad data and with autodesk viewer
  G4VisAttributes *Aluminum = new G4VisAttributes(G4Colour(0.3,0.3,1.0));
  G4VisAttributes *Iron = new G4VisAttributes(G4Colour(0.3,0.3,0.3));
  G4VisAttributes *LeadColor = new G4VisAttributes(G4Colour(0.4,0.4,0.4));
  G4VisAttributes *DebugRed = new G4VisAttributes(G4Colour(1.0,0.,0.));

  //G4Tubs *P1ringA = new G4Tubs("P1ringA", P1ringA_rin, P1ringA_rin+P1tubeTh, P1ringA_L, 0.*deg, 360.*deg);

  //Ring A - Most proximal ring to target chamber
  G4double P1ringA_L = 0.11/2*inch;
  G4double P1ringA_rin = 2.93/2*inch;

  G4Tubs *P1ringA = new G4Tubs("P1ringA", P1ringA_rin, P1ringA_rin+P1tubeTh, P1ringA_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringALog = new G4LogicalVolume(P1ringA, GetMaterial("Aluminum"), "P1ringA_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1ringA_L), P1ringALog, "P1ringALog_pv", worldlog, false, 0 , ChkOverlaps );

  //Disk - Empty space in cad file ostensibly left open to house Be window containing vacuum through exit beamline. Will have to check this.
  G4double P1disk_L = 0.015/2*inch;

  G4Tubs *P1disk = new G4Tubs("P1disk", 0.0, P1ringA_rin+P1tubeTh, P1disk_L, 0.*deg, 360.*deg);
 
  G4LogicalVolume *P1diskLog = new G4LogicalVolume(P1disk, GetMaterial("Beryllium"), "P1disk_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+P1disk_L), P1diskLog, "P1diskLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube A - First of three ascending radius tubes. Mates on distant end with ring.
  G4double P1tubeA_L = 1.171/2*inch;

  G4Tubs *P1tubeA = new G4Tubs("P1tubeA", P1ringA_rin, P1ringA_rin+P1tubeTh, P1tubeA_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeALog = new G4LogicalVolume(P1tubeA, GetMaterial("Aluminum"), "P1tubeA_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+P1tubeA_L), P1tubeALog, "P1tubeALog_pv", worldlog, false, 0 , ChkOverlaps );
 
  //Ring Bin - Small increase in outer radius
  G4double P1ringBin_rou = 3.08/2*inch;

  G4Tubs *P1ringBin = new G4Tubs("P1ringBin", P1ringA_rin, P1ringBin_rou, P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringBinLog = new G4LogicalVolume(P1ringBin, GetMaterial("Aluminum"), "P1ringBin_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+P1ringL), P1ringBinLog, "P1ringBinLog_pv", worldlog, false, 0 , ChkOverlaps );
  
  //Ring Bou - Large increase in outer radius to match mating tube
  G4double P1ringBou_rou = 3.75/2*inch;

  G4Tubs *P1ringBou = new G4Tubs("P1ringBou", P1ringA_rin, P1ringBou_rou, P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringBouLog = new G4LogicalVolume(P1ringBou, GetMaterial("Aluminum"), "P1ringBou_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+2*P1ringL+P1ringL), P1ringBouLog, "P1ringBouLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube B - Second of three ascending radius tubes. Mates on both ends with rings.
  G4double P1tubeB_L = 2.377/2*inch;
  G4double P1tubeB_rin = 3.68/2*inch;

  G4Tubs *P1tubeB = new G4Tubs("P1tubeB", P1tubeB_rin, P1tubeB_rin+P1tubeTh, P1tubeB_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeBLog = new G4LogicalVolume(P1tubeB, GetMaterial("Aluminum"), "P1tubeB_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+4*P1ringL+P1tubeB_L), P1tubeBLog, "P1tubeBLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring Cin - Small increase in outer radius
  G4double P1ringCin_rou = 3.83/2*inch;

  G4Tubs *P1ringCin = new G4Tubs("P1ringCin", P1tubeB_rin, P1ringCin_rou, P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringCinLog = new G4LogicalVolume(P1ringCin, GetMaterial("Aluminum"), "P1ringCin_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+4*P1ringL+2*P1tubeB_L+P1ringL), P1ringCinLog, "P1ringCinLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring Cou - Large increase in outer radius to match mating tube
  G4double P1ringCou_rou = 4.5/2*inch;

  G4Tubs *P1ringCou = new G4Tubs("P1ringCou", P1tubeB_rin, P1ringCou_rou, P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringCouLog = new G4LogicalVolume(P1ringCou, GetMaterial("Aluminum"), "P1ringCou_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+6*P1ringL+2*P1tubeB_L+P1ringL), P1ringCouLog, "P1ringCouLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube C - Third of three ascending radius tubes. Mates on both ends with rings.
  G4double P1tubeC_L = 2.187/2*inch;
  G4double P1tubeC_rin = 4.43/2*inch;

  G4Tubs *P1tubeC = new G4Tubs("P1tubeC", P1tubeC_rin, P1tubeC_rin+P1tubeTh, P1tubeC_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeCLog = new G4LogicalVolume(P1tubeC, GetMaterial("Aluminum"), "P1tubeC_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+8*P1ringL+2*P1tubeB_L+P1tubeC_L), P1tubeCLog, "P1tubeCLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring D in and out twice the length of previous rings - remaining elements deviate from general specifications above
  //Ring Din
  G4double P1ringDin_rou = 4.57/2*inch;

  G4Tubs *P1ringDin = new G4Tubs("P1ringDin", P1tubeC_rin, P1ringDin_rou, 2*P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringDinLog = new G4LogicalVolume(P1ringDin, GetMaterial("Aluminum"), "P1ringDin_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+8*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1ringL), P1ringDinLog, "P1ringDinLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring Dou
  G4double P1ringDou_rou = 6.0/2*inch;

  G4Tubs *P1ringDou = new G4Tubs("P1ringCou", P1tubeC_rin, P1ringDou_rou, 2*P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringDouLog = new G4LogicalVolume(P1ringDou, GetMaterial("Aluminum"), "P1ringDou_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+12*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1ringL), P1ringDouLog, "P1ringDouLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube D
  G4double P1tubeD_L = 19.219/2*inch;
  G4double P1tubeD_rin = 5.75/2*inch;
  G4double P1tubeD_rou = 6.0/2*inch;

  G4Tubs *P1tubeD = new G4Tubs("P1tubeD", P1tubeD_rin, P1tubeD_rou, P1tubeD_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeDLog = new G4LogicalVolume(P1tubeD, GetMaterial("Aluminum"), "P1tubeD_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+P1tubeD_L), P1tubeDLog, "P1tubeDLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring E
  G4double P1ringE_L = 0.5/2*inch;
  G4double P1ringE_rou = 6.299/2*inch;

  G4Tubs *P1ringE = new G4Tubs("P1ringE", P1tubeD_rin, P1ringE_rou, P1ringE_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringELog = new G4LogicalVolume(P1ringE, GetMaterial("Aluminum"), "P1ringE_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1tubeD_L+P1ringE_L), P1ringELog, "P1ringELog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube E
  G4double P1tubeE_L = 2.571/2*inch;
  G4double P1tubeE_rin = 5.906/2*inch;

  G4Tubs *P1tubeE = new G4Tubs("P1tubeD", P1tubeE_rin, P1ringE_rou, P1tubeE_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeELog = new G4LogicalVolume(P1tubeE, GetMaterial("Aluminum"), "P1tubeE_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1tubeD_L+2*P1ringE_L+P1tubeE_L), P1tubeELog, "P1tubeELog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring F
  G4double P1ringF_L = 0.866/2*inch;
  G4double P1ringF_rou = 7.97/2*inch;

  G4Tubs *P1ringF = new G4Tubs("P1ringF", P1tubeE_rin, P1ringF_rou, P1ringF_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringFLog = new G4LogicalVolume(P1ringF, GetMaterial("Aluminum"), "P1ringF_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1tubeD_L+2*P1ringE_L+2*P1tubeE_L+P1ringF_L), P1ringFLog, "P1ringFLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring G
  G4double P1ringG_rin = 3.87/2*inch;

  G4Tubs *P1ringG = new G4Tubs("P1ringG", P1ringG_rin, P1ringF_rou, P1ringF_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringGLog = new G4LogicalVolume(P1ringG, GetMaterial("Aluminum"), "P1ringG_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1tubeD_L+2*P1ringE_L+2*P1tubeE_L+2*P1ringF_L+P1ringF_L), P1ringGLog, "P1ringGLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube F
  G4double P1tubeF_L = 4.28/2*inch;
  G4double P1tubeF_rou = 4.0/2*inch;

  G4Tubs *P1tubeF = new G4Tubs("P1tubeF", P1ringG_rin, P1tubeF_rou, P1tubeF_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeFLog = new G4LogicalVolume(P1tubeF, GetMaterial("Aluminum"), "P1tubeF_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1tubeD_L+2*P1ringE_L+2*P1tubeE_L+4*P1ringF_L+P1tubeF_L), P1tubeFLog, "P1tubeFLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring H
  G4double P1ringH_L = 0.42/2*inch;
  G4double P1ringH_rou = 6.7/2*inch;

  G4Tubs *P1ringH = new G4Tubs("P1ringH", P1ringG_rin, P1ringH_rou, P1ringH_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringHLog = new G4LogicalVolume(P1ringH, GetMaterial("Aluminum"), "P1ringH_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1tubeD_L+2*P1ringE_L+2*P1tubeE_L+4*P1ringF_L+2*P1tubeF_L+P1ringH_L), P1ringHLog, "P1ringHLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring I
  G4double P1ringI_rin = 3.92/2*inch;

  G4Tubs *P1ringI = new G4Tubs("P1ringI", P1ringI_rin, P1ringH_rou, P1ringH_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringILog = new G4LogicalVolume(P1ringI, GetMaterial("Aluminum"), "P1ringI_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1tubeD_L+2*P1ringE_L+2*P1tubeE_L+4*P1ringF_L+2*P1tubeF_L+2*P1ringH_L+P1ringH_L), P1ringILog, "P1ringILog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring J
  //This ring mates with section two cone at 3.571 inches diameter at entrance of cone and difference is ignored.
  G4double P1ringJ_L = 0.84/2*inch;
  G4double P1ringJ_rin = 3.5/2*inch;

  G4Tubs *P1ringJ = new G4Tubs("P1ringJ", P1ringJ_rin, P1ringH_rou, P1ringJ_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringJLog = new G4LogicalVolume(P1ringJ, GetMaterial("Aluminum"), "P1ringJ_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1tubeD_L+2*P1ringE_L+2*P1tubeE_L+4*P1ringF_L+2*P1tubeF_L+4*P1ringH_L+P1ringJ_L), P1ringJLog, "P1ringJLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Vacuum fill solids - adding an object for each new rin

  //Tube A Vacuum
  G4Tubs *P1tubeA_vac = new G4Tubs("P1tubeA_vac", 0.0, P1ringA_rin, P1tubeA_L+2*P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeA_vacLog = new G4LogicalVolume(P1tubeA_vac, GetMaterial("Vacuum"), "P1tubeA_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1tubeA_L+2*P1ringL), P1tubeA_vacLog, "P1tubeA_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube B Vacuum
  G4Tubs *P1tubeB_vac = new G4Tubs("P1tubeB_vac", 0.0, P1tubeB_rin, P1tubeB_L+2*P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeB_vacLog = new G4LogicalVolume(P1tubeB_vac, GetMaterial("Vacuum"), "P1tubeB_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1tubeA_L+4*P1ringL+P1tubeB_L+2*P1ringL), P1tubeB_vacLog, "P1tubeB_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube C Vacuum
  G4Tubs *P1tubeC_vac = new G4Tubs("P1tubeC_vac", 0.0, P1tubeC_rin, P1tubeC_L+4*P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeC_vacLog = new G4LogicalVolume(P1tubeC_vac, GetMaterial("Vacuum"), "P1tubeC_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1tubeA_L+8*P1ringL+2*P1tubeB_L+P1tubeC_L+4*P1ringL), P1tubeC_vacLog, "P1tubeC_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube D Vacuum
  G4Tubs *P1tubeD_vac = new G4Tubs("P1tubeD_vac", 0.0, P1tubeD_rin, P1tubeD_L+P1ringE_L, 0.*deg, 360.*deg);
  
  G4LogicalVolume *P1tubeD_vacLog = new G4LogicalVolume(P1tubeD_vac, GetMaterial("Vacuum"), "P1tubeD_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+P1tubeD_L+P1ringE_L), P1tubeD_vacLog, "P1tubeD_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube E Vacuum
  G4Tubs *P1tubeE_vac = new G4Tubs("P1tubeE_vac", 0.0, P1tubeE_rin, P1tubeE_L+P1ringF_L, 0.*deg, 360.*deg);
  
  G4LogicalVolume *P1tubeE_vacLog = new G4LogicalVolume(P1tubeE_vac, GetMaterial("Vacuum"), "P1tubeE_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1tubeD_L+2*P1ringE_L+P1tubeE_L+P1ringF_L), P1tubeE_vacLog, "P1tubeE_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube F Vacuum
  G4Tubs *P1tubeF_vac = new G4Tubs("P1tubeF_vac", 0.0, P1ringG_rin, P1ringF_L+P1tubeF_L+P1ringH_L, 0.*deg, 360.*deg);
  
  G4LogicalVolume *P1tubeF_vacLog = new G4LogicalVolume(P1tubeF_vac, GetMaterial("Vacuum"), "P1tubeF_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1tubeD_L+2*P1ringE_L+2*P1tubeE_L+2*P1ringF_L+P1ringF_L+P1tubeF_L+P1ringH_L), P1tubeF_vacLog, "P1tubeF_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring I Vacuum
  G4Tubs *P1ringI_vac = new G4Tubs("P1ringI_vac", 0.0, P1ringI_rin, P1ringH_L, 0.*deg, 360.*deg);
  
  G4LogicalVolume *P1ringI_vacLog = new G4LogicalVolume(P1ringI_vac, GetMaterial("Vacuum"), "P1ringI_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1tubeD_L+2*P1ringE_L+2*P1tubeE_L+2*P1ringF_L+2*P1ringF_L+2*P1tubeF_L+2*P1ringH_L+P1ringH_L), P1ringI_vacLog, "P1ringI_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring J Vacuum
  G4Tubs *P1ringJ_vac = new G4Tubs("P1ringJ_vac", 0.0, P1ringG_rin, P1ringJ_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringJ_vacLog = new G4LogicalVolume(P1ringJ_vac, GetMaterial("Vacuum"), "P1ringJ_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1tubeD_L+2*P1ringE_L+2*P1tubeE_L+2*P1ringF_L+2*P1ringF_L+2*P1tubeF_L+4*P1ringH_L+P1ringJ_L), P1ringJ_vacLog, "P1ringJ_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  /*
  G4LogicalVolume *P1ringALog = new G4LogicalVolume(P1ringA, GetMaterial("Aluminum"), "P1ringA_log", 0, 0, 0);
  G4LogicalVolume *P1diskLog = new G4LogicalVolume(P1disk, GetMaterial("Aluminum"), "P1disk_log", 0, 0, 0);
  G4LogicalVolume *P1tubeALog = new G4LogicalVolume(P1tubeA, GetMaterial("Aluminum"), "P1tubeA_log", 0, 0, 0);
  G4LogicalVolume *P1ringBinLog = new G4LogicalVolume(P1ringBin, GetMaterial("Aluminum"), "P1ringBin_log", 0, 0, 0);
  G4LogicalVolume *P1ringBouLog = new G4LogicalVolume(P1ringBou, GetMaterial("Aluminum"), "P1ringBou_log", 0, 0, 0);
  G4LogicalVolume *P1tubeBLog = new G4LogicalVolume(P1tubeB, GetMaterial("Aluminum"), "P1tubeB_log", 0, 0, 0);
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z), P1ringALog, "P1ringALog_pv", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1ringA_L), P1diskLog, "P1diskLog_pv", worldlog, false, 0 , ChkOverlaps );

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z), P1ringALog, "P1ringALog_pv", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1ringA_L), P1diskLog, "P1diskLog_pv", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1ringA_L+P1disk_L), P1tubeALog, "P1tubeALog_pv", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1ringA_L+P1disk_L+P1tubeA_L), P1ringBinLog, "P1ringBinLog_pv", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1ringA_L+P1disk_L+P1tubeA_L+P1ringL), P1ringBouLog, "P1ringBouLog_pv", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1ringA_L+P1disk_L+P1tubeA_L+2*P1ringL), P1tubeBLog, "P1tubeBLog_pv", worldlog, false, 0 , ChkOverlaps );

  P1ringALog->SetVisAttributes( Aluminum);
  P1diskLog->SetVisAttributes( Aluminum);
  P1tubeALog->SetVisAttributes( Aluminum);
  P1ringBinLog->SetVisAttributes( Aluminum);
  P1ringBouLog->SetVisAttributes( Aluminum);
  P1tubeBLog->SetVisAttributes( Aluminum);
  */

  //Placement of Elements
  //Distance to origin (presumably the target... must scour code)
  //G4double initPlacement_z = ent_len/2+40.0*cm;
  //Positions of cylinders in proximal order to target chamber

  //Too clunky - try something else. Create vectors which contain positions in proximal order from target - something like the following.
   
  /*
  G4double P1SolidLengths[] = { P1ringA_L, P1disk_L, P1tubeA_L, P1ringL, P1ringL, P1tubeB_L, P1ringL, P1ringL, P1tubeC_L, 2*P1ringL, 2*P1ringL, P1tubeC_L, P1ringE_L, P1tubeE_L, P1ringF_L, P1ringF_L, P1tubeF_L, P1ringH_L, p1ringH_L, P1ringJ_L};
  G4double P1VacLengths[] = { P1tubeA_L+2*P1ringL, P1tubeB_L+2*P1ringL, P1tubeC_L+4*P1ringL, P1tubeD_L+2*P1ringE_L, P1tubeE_L+2*P1ringF_L, P1ringG_L+P1tubeF_L+P1ringH_L, P1ringH_L, P1ringJ_L};

  G4double P1Pos[] = { 0.0, P1ringA_L};
  G4double P1PosVac[] = { P1ringA_L+P1disk_L, P1tubeA_L+2*P1ringL}
  int index = 2;

  while (index < P1SolidLengths.size()){
    P1Pos.push_back(P1Pos[index-1]+P1SolidLengths[index]);
    if(index < P1VacLengths.size()){
      P1PosVac.push_back(P1PosVac[index-1]+P1VacLengths[index]);
    }
    index++;
  }

  //Now place all elements with these values
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[0]), P1ringALog, "P1ringALog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[1]), P1diskLog, "P1diskLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[2]), P1tubeALog, "P1tubeALog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[3]), P1ringBinLog, "P1ringBinLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[4]), P1ringBouLog, "P1ringBouLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[5]), P1tubeBLog, "P1tubeBLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[6]), P1ringCinLog, "P1ringCinLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[7]), P1ringCouLog, "P1ringCouLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[8]), P1tubeCLog, "P1tubeCLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[9]), P1ringDinLog, "P1ringDinLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[10]), P1ringDouLog, "P1ringDouLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[11]), P1tubeDLog, "P1tubeDLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[12]), P1ringELog, "P1ringELog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[13]), P1tubeELog, "P1tubeELog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[14]), P1ringFLog, "P1ringFLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[15]), P1ringGLog, "P1ringGLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[16]), P1tubeFLog, "P1tubeFLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[17]), P1ringHLog, "P1ringHLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[18]), P1ringILog, "P1ringILog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1Pos[19]), P1ringJLog, "P1ringJLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1PosVac[0]), P1tubeA_vacLog, "P1tubeA_vacLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1PosVac[1]), P1tubeB_vacLog, "P1tubeB_vacLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1PosVac[2]), P1tubeC_vacLog, "P1tubeC_vacLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1PosVac[3]), P1tubeD_vacLog, "P1tubeD_vacLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1PosVac[4]), P1tubeE_vacLog, "P1tubeE_vacLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1PosVac[5]), P1tubeF_vacLog, "P1tubeF_vacLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1PosVac[6]), P1ringI_vacLog, "P1ringI_vacLog_pv", worldlog, 0, false, 0 );
  new G4PVPlacement( 0.0, G4ThreeVector( 0.0, 0.0, initPlacement_z+P1PosVac[7]), P1ringJ_vacLog, "P1ringJ_vacLog_pv", worldlog, 0, false, 0 );
  */

  
  //Visuals
  //G4VisAttributes *Aluminum = new G4VisAttributes( G4Colour( 0.3, 0.3, 1.0));
  G4VisAttributes * WireFrameVisAtt = new G4VisAttributes(G4Colour(0.,1.,1.));
  WireFrameVisAtt->SetForceWireframe(true);
  //P1tubeBLog->SetVisAttributes( WireFrameVisAtt );
  

  P1ringALog->SetVisAttributes( Aluminum);
  P1diskLog->SetVisAttributes( Aluminum);
  P1tubeALog->SetVisAttributes( Aluminum);
  P1ringBinLog->SetVisAttributes( Aluminum);
  P1ringBouLog->SetVisAttributes( Aluminum);
  P1tubeBLog->SetVisAttributes( Aluminum);
  P1ringCinLog->SetVisAttributes( Aluminum);
  P1ringCouLog->SetVisAttributes( Aluminum);
  P1tubeCLog->SetVisAttributes( Aluminum);
  P1ringDinLog->SetVisAttributes( Aluminum);
  P1ringDouLog->SetVisAttributes( Aluminum);
  P1tubeDLog->SetVisAttributes( Aluminum);
  P1ringELog->SetVisAttributes( Aluminum);
  P1tubeELog->SetVisAttributes( Aluminum);
  P1ringFLog->SetVisAttributes( Aluminum);
  P1ringGLog->SetVisAttributes( Aluminum);
  P1tubeFLog->SetVisAttributes( Aluminum);
  P1ringHLog->SetVisAttributes( Aluminum);
  P1ringILog->SetVisAttributes( Aluminum);
  P1ringJLog->SetVisAttributes( Aluminum);
  P1tubeA_vacLog->SetVisAttributes( G4VisAttributes::Invisible);
  P1tubeB_vacLog->SetVisAttributes( G4VisAttributes::Invisible);
  P1tubeC_vacLog->SetVisAttributes( G4VisAttributes::Invisible);
  P1tubeD_vacLog->SetVisAttributes( G4VisAttributes::Invisible);
  P1tubeE_vacLog->SetVisAttributes( G4VisAttributes::Invisible);
  P1tubeF_vacLog->SetVisAttributes( G4VisAttributes::Invisible);
  P1ringI_vacLog->SetVisAttributes( G4VisAttributes::Invisible);
  P1ringJ_vacLog->SetVisAttributes( G4VisAttributes::Invisible);
  
  //Will need to confirm materials throughout section one


  //SECTION TWO 

  //This section includes all components up to and including the first corrector magnet.
  //Specifications - ordered first from inside to out then from target-proximal to distant
  
  //G4double P2initPlacement_z = initPlacement_z+P1tubeA_L+8*P1ringL+P1tubeB_L+P1tubeC_L+P1tubeD_L+P1ringE_L+P1tubeE_L+P1ringF_L+P1ringF_L+P1tubeF_L+2*P1ringH_L+P1ringJ_L;
  
  G4double P2initPlacement_z = initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1tubeD_L+2*P1ringE_L+2*P1tubeE_L+4*P1ringF_L+2*P1tubeF_L+4*P1ringH_L+2*P1ringJ_L;
  G4double P2_offset1 = -0.309*inch; //Long inner cone is seated within P1ringJ
  G4double P2_offset2 = 0.141*inch; //First shielding cone is separated from P1ringJ

  //Inner cone
  G4double P2coneA_rin1 = 3.517/2*inch;
  G4double P2coneA_rou1 = 3.767/2*inch;
  //G4double P2coneA_rin2 = 5.446*inch; //Section one rin2
  //G4double P2coneA_rou2 = 5.696*inch; //Section one rou2
  G4double P2coneA_rin2 = 10.734/2*inch;
  G4double P2coneA_rou2 = 10.984/2*inch;
  //G4double P2coneA_L = 45.55*inch; //Section one length only
  G4double P2coneA_L = 137.8/2*inch; //subtracting mating section with P5 and including volume with first ring there -> -0.620*inch

  G4Cons *P2coneA = new G4Cons("P2coneA", P2coneA_rin1, P2coneA_rou1, P2coneA_rin2, P2coneA_rou2, P2coneA_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P2coneALog = new G4LogicalVolume(P2coneA, GetMaterial("Aluminum"), "P2coneA_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset1+P2coneA_L), P2coneALog, "P2coneALog_pv", worldlog, false, 0 , ChkOverlaps);

  G4VisAttributes *DebugGreen = new G4VisAttributes(G4Colour(0.,1.,0.));
  //P2coneALog->SetVisAttributes(DebugGreen);
  
  P2coneALog->SetVisAttributes( Aluminum);

  //Inner Cone Vacuum
  G4Cons *P2coneA_vac = new G4Cons("P2coneA_vac", 0.0, P2coneA_rin1, 0.0, P2coneA_rin2, P2coneA_L, 0.*deg, 360.*deg);
  
  G4LogicalVolume *P2coneA_vacLog = new G4LogicalVolume(P2coneA_vac, GetMaterial("Vacuum"), "P2coneA_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset1+P2coneA_L), P2coneA_vacLog, "P2coneA_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Middle cone - staged, wrapped around inner cone
  G4double P2coneB_rin1 = 3.517/2*inch;
  G4double P2coneB_rou1 = 3.767/2*inch;
  G4double P2coneB_rin2 = 5.903/2*inch;
  G4double P2coneB_rou2 = 6.150/2*inch;
  G4double P2coneB_L = 33.625/2*inch;
  
  G4Cons *P2coneB = new G4Cons("P2coneB", P2coneB_rin1, P2coneB_rou1, P2coneB_rin2, P2coneB_rou2, P2coneB_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P2coneBLog = new G4LogicalVolume(P2coneB, GetMaterial("Aluminum"), "P2coneB_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+P2coneB_L), P2coneBLog, "P2coneBLog_pv", worldlog, false, 0 , ChkOverlaps);

  //P2coneBLog->SetVisAttributes( DebugRed );
  P2coneBLog->SetVisAttributes( Aluminum);

  //Rings General Specifications, 18 in total
  G4double P2ringTh = 0.5*inch;
  G4double P2ringr0 = 1.895*inch;
  G4double P2ringL = 1.625/2*inch;
  G4double P2DAngle = 1.5*deg;
  G4double P2ringSep = 0.375*inch;

  //Ring 1
  G4double P2ring1_rin1 = P2ringr0;
  G4double P2ring1_rin2 = P2ringr0+2*P2ringL*tan(P2DAngle);

  G4Cons *P2ring1 = new G4Cons("P2ring1", P2ring1_rin1, P2ring1_rin1+P2ringTh, P2ring1_rin2, P2ring1_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring1Log = new G4LogicalVolume(P2ring1, GetMaterial("Iron"), "P2ring1_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+P2ringL), P2ring1Log, "P2ring1Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring1Log->SetVisAttributes(Iron);

  //Before this mark, updated and working 6.15

  //Ring 2
  G4double P2ring2_rin1 = P2ringr0+(2*P2ringL+P2ringSep)*tan(P2DAngle);
  G4double P2ring2_rin2 = P2ringr0+(4*P2ringL+P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring2 = new G4Cons("P2ring2", P2ring2_rin1, P2ring2_rin1+P2ringTh, P2ring2_rin2, P2ring2_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring2Log = new G4LogicalVolume(P2ring2, GetMaterial("Iron"), "P2ring2_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+2*P2ringL+P2ringSep+P2ringL), P2ring2Log, "P2ring2Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring2Log->SetVisAttributes(Iron);

  //Ring 3
  G4double P2ring3_rin1 = P2ringr0+(4*P2ringL+2*P2ringSep)*tan(P2DAngle);
  G4double P2ring3_rin2 = P2ringr0+(6*P2ringL+2*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring3 = new G4Cons("P2ring3", P2ring3_rin1, P2ring3_rin1+P2ringTh, P2ring3_rin2, P2ring3_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring3Log = new G4LogicalVolume(P2ring3, GetMaterial("Iron"), "P2ring3_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+4*P2ringL+2*P2ringSep+P2ringL), P2ring3Log, "P2ring3Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring3Log->SetVisAttributes(Iron);

  //Ring 4
  G4double P2ring4_rin1 = P2ringr0+(6*P2ringL+3*P2ringSep)*tan(P2DAngle);
  G4double P2ring4_rin2 = P2ringr0+(8*P2ringL+3*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring4 = new G4Cons("P2ring4", P2ring4_rin1, P2ring4_rin1+P2ringTh, P2ring4_rin2, P2ring4_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring4Log = new G4LogicalVolume(P2ring4, GetMaterial("Iron"), "P2ring4_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+6*P2ringL+3*P2ringSep+P2ringL), P2ring4Log, "P2ring4Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring4Log->SetVisAttributes(Iron);

  //Ring 5
  G4double P2ring5_rin1 = P2ringr0+(8*P2ringL+4*P2ringSep)*tan(P2DAngle);
  G4double P2ring5_rin2 = P2ringr0+(10*P2ringL+4*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring5 = new G4Cons("P2ring5", P2ring5_rin1, P2ring5_rin1+P2ringTh, P2ring5_rin2, P2ring5_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring5Log = new G4LogicalVolume(P2ring5, GetMaterial("Iron"), "P2ring5_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+8*P2ringL+4*P2ringSep+P2ringL), P2ring5Log, "P2ring5Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring5Log->SetVisAttributes(Iron);

  //Ring 6
  G4double P2ring6_rin1 = P2ringr0+(10*P2ringL+5*P2ringSep)*tan(P2DAngle);
  G4double P2ring6_rin2 = P2ringr0+(12*P2ringL+5*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring6 = new G4Cons("P2ring6", P2ring6_rin1, P2ring6_rin1+P2ringTh, P2ring6_rin2, P2ring6_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring6Log = new G4LogicalVolume(P2ring6, GetMaterial("Iron"), "P2ring6_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+10*P2ringL+5*P2ringSep+P2ringL), P2ring6Log, "P2ring6Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring6Log->SetVisAttributes(Iron);

  //Ring 7
  G4double P2ring7_rin1 = P2ringr0+(12*P2ringL+6*P2ringSep)*tan(P2DAngle);
  G4double P2ring7_rin2 = P2ringr0+(14*P2ringL+6*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring7 = new G4Cons("P2ring7", P2ring7_rin1, P2ring7_rin1+P2ringTh, P2ring7_rin2, P2ring7_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring7Log = new G4LogicalVolume(P2ring7, GetMaterial("Iron"), "P2ring7_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+12*P2ringL+6*P2ringSep+P2ringL), P2ring7Log, "P2ring7Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring7Log->SetVisAttributes(Iron);

  //Ring 8
  G4double P2ring8_rin1 = P2ringr0+(14*P2ringL+7*P2ringSep)*tan(P2DAngle);
  G4double P2ring8_rin2 = P2ringr0+(16*P2ringL+7*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring8 = new G4Cons("P2ring8", P2ring8_rin1, P2ring8_rin1+P2ringTh, P2ring8_rin2, P2ring8_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring8Log = new G4LogicalVolume(P2ring8, GetMaterial("Iron"), "P2ring8_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+14*P2ringL+7*P2ringSep+P2ringL), P2ring8Log, "P2ring8Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring8Log->SetVisAttributes(Iron);

  //Ring 9
  G4double P2ring9_rin1 = P2ringr0+(16*P2ringL+8*P2ringSep)*tan(P2DAngle);
  G4double P2ring9_rin2 = P2ringr0+(18*P2ringL+8*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring9 = new G4Cons("P2ring9", P2ring9_rin1, P2ring9_rin1+P2ringTh, P2ring9_rin2, P2ring9_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring9Log = new G4LogicalVolume(P2ring9, GetMaterial("Iron"), "P2ring9_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+16*P2ringL+8*P2ringSep+P2ringL), P2ring9Log, "P2ring9Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring9Log->SetVisAttributes(Iron);

  //Ring 10
  G4double P2ring10_rin1 = P2ringr0+(18*P2ringL+9*P2ringSep)*tan(P2DAngle);
  G4double P2ring10_rin2 = P2ringr0+(20*P2ringL+9*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring10 = new G4Cons("P2ring10", P2ring10_rin1, P2ring10_rin1+P2ringTh, P2ring10_rin2, P2ring10_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring10Log = new G4LogicalVolume(P2ring10, GetMaterial("Iron"), "P2ring10_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+18*P2ringL+9*P2ringSep+P2ringL), P2ring10Log, "P2ring10Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring10Log->SetVisAttributes(Iron);

  //Ring 11
  G4double P2ring11_rin1 = P2ringr0+(20*P2ringL+10*P2ringSep)*tan(P2DAngle);
  G4double P2ring11_rin2 = P2ringr0+(22*P2ringL+10*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring11 = new G4Cons("P2ring11", P2ring11_rin1, P2ring11_rin1+P2ringTh, P2ring11_rin2, P2ring11_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring11Log = new G4LogicalVolume(P2ring11, GetMaterial("Iron"), "P2ring11_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+20*P2ringL+10*P2ringSep+P2ringL), P2ring11Log, "P2ring11Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring11Log->SetVisAttributes(Iron);

  //Ring 12
  G4double P2ring12_rin1 = P2ringr0+(22*P2ringL+11*P2ringSep)*tan(P2DAngle);
  G4double P2ring12_rin2 = P2ringr0+(24*P2ringL+11*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring12 = new G4Cons("P2ring12", P2ring12_rin1, P2ring12_rin1+P2ringTh, P2ring12_rin2, P2ring12_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring12Log = new G4LogicalVolume(P2ring12, GetMaterial("Iron"), "P2ring12_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+22*P2ringL+11*P2ringSep+P2ringL), P2ring12Log, "P2ring12Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring12Log->SetVisAttributes(Iron);

  //Ring 13
  G4double P2ring13_rin1 = P2ringr0+(24*P2ringL+12*P2ringSep)*tan(P2DAngle);
  G4double P2ring13_rin2 = P2ringr0+(26*P2ringL+12*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring13 = new G4Cons("P2ring13", P2ring13_rin1, P2ring13_rin1+P2ringTh, P2ring13_rin2, P2ring13_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring13Log = new G4LogicalVolume(P2ring13, GetMaterial("Iron"), "P2ring13_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+24*P2ringL+12*P2ringSep+P2ringL), P2ring13Log, "P2ring13Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring13Log->SetVisAttributes(Iron);

  //Ring 14
  G4double P2ring14_rin1 = P2ringr0+(26*P2ringL+13*P2ringSep)*tan(P2DAngle);
  G4double P2ring14_rin2 = P2ringr0+(28*P2ringL+13*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring14 = new G4Cons("P2ring14", P2ring14_rin1, P2ring14_rin1+P2ringTh, P2ring14_rin2, P2ring14_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring14Log = new G4LogicalVolume(P2ring14, GetMaterial("Iron"), "P2ring14_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+26*P2ringL+13*P2ringSep+P2ringL), P2ring14Log, "P2ring14Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring14Log->SetVisAttributes(Iron);

  //Ring 15
  G4double P2ring15_rin1 = P2ringr0+(28*P2ringL+14*P2ringSep)*tan(P2DAngle);
  G4double P2ring15_rin2 = P2ringr0+(30*P2ringL+14*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring15 = new G4Cons("P2ring15", P2ring15_rin1, P2ring15_rin1+P2ringTh, P2ring15_rin2, P2ring15_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring15Log = new G4LogicalVolume(P2ring15, GetMaterial("Iron"), "P2ring15_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+28*P2ringL+14*P2ringSep+P2ringL), P2ring15Log, "P2ring15Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring15Log->SetVisAttributes(Iron);

  //Ring 16
  G4double P2ring16_rin1 = P2ringr0+(30*P2ringL+15*P2ringSep)*tan(P2DAngle);
  G4double P2ring16_rin2 = P2ringr0+(32*P2ringL+15*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring16 = new G4Cons("P2ring16", P2ring16_rin1, P2ring16_rin1+P2ringTh, P2ring16_rin2, P2ring16_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring16Log = new G4LogicalVolume(P2ring16, GetMaterial("Iron"), "P2ring16_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+30*P2ringL+15*P2ringSep+P2ringL), P2ring16Log, "P2ring16Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring16Log->SetVisAttributes(Iron);

  //Ring 17
  G4double P2ring17_rin1 = P2ringr0+(32*P2ringL+16*P2ringSep)*tan(P2DAngle);
  G4double P2ring17_rin2 = P2ringr0+(34*P2ringL+16*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring17 = new G4Cons("P2ring17", P2ring17_rin1, P2ring17_rin1+P2ringTh, P2ring17_rin2, P2ring17_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring17Log = new G4LogicalVolume(P2ring17, GetMaterial("Iron"), "P2ring17_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+32*P2ringL+16*P2ringSep+P2ringL), P2ring17Log, "P2ring17Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring17Log->SetVisAttributes(Iron);

  //Ring 18
  G4double P2ring18_rin1 = P2ringr0+(34*P2ringL+17*P2ringSep)*tan(P2DAngle);
  G4double P2ring18_rin2 = P2ringr0+(36*P2ringL+17*P2ringSep)*tan(P2DAngle);

  G4Cons *P2ring18 = new G4Cons("P2ring18", P2ring18_rin1, P2ring18_rin1+P2ringTh, P2ring18_rin2, P2ring18_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P2ring18Log = new G4LogicalVolume(P2ring18, GetMaterial("Iron"), "P2ring18_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset2+34*P2ringL+17*P2ringSep+P2ringL), P2ring18Log, "P2ring18Log_pv", worldlog, false, 0 , ChkOverlaps);

  P2ring18Log->SetVisAttributes(Iron);


  //Clunky - more efficient method with loop, but don't know how make file parses vector of G4 objects - will need to look into this.
  /*
  G4double P2ring_rin1[] = {};
  G4double P2ring_rin2[] = {};
  
  G4Cons P2ring[] = {};

  G4LogicalVolume *P2ringLog[] = {};

  G4PVPlacement *P2ringLogpv[] = {};

  
  for (int i=0, i<18){

  P2ring_rin1[i] = P2ringr0+(2*i*P2ringL+i*P2ringSep)*tan(P2DAngle);
  P2ring_rin2[i] = P2ringr0+(2*(i+1)*P2ringL+i*P2ringSep)*tan(P2DAngle);

  string P2label1text = "P2ring" + to_string([i]);
  string P2label2text = "P2ringLog" + to_string([i]);
  string P2label3text = "P2ringLogPV" + to_string([i]);

  P2ring[i] = new G4Cons( to_string(P2label1text) , P2ring_rin1[i], P2ring_rin1[i]+P2ringTh, P2ring_rin2[i], P2ring_rin2[i]+P2ringTh, P2ringL, 0.*deg, 360.*deg);

  P2ringLog[i] = new G4LogicalVolume( P2label1text, GetMaterial("Iron"), to_string(P2label2text), 0, 0, 0);

  P2ringLogpv[i] = new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset+2*i*P2ringL+i*P2ringSep+P1ringL), P2label2text, to_string(P2label3text), worldlog, false, 0 , ChkOverlaps);

  }
  */

  
  //Side shields - NOT shields per material update 7.24.20
  //Specifications
  G4double sideshieldTh = 0.25*inch;
  G4double sideshieldA = 1.06*deg;
  G4double P2sideshieldL = 37.119/2*inch;
  G4double P2sideshieldW1 = 4.762/2*inch;
  G4double P2sideshieldW2 = 6.137/2*inch;
  G4double P2sideshield_xoffset = 8.282/2*inch+sideshieldTh;
  G4double P2sideshield_zoffset = 0.116*inch;

  //Trapezoid class 
  G4Trd *P2sideshield1 = new G4Trd("P2sideshield1", P2sideshieldW1, P2sideshieldW2, sideshieldTh, sideshieldTh, P2sideshieldL);
  G4Trd *P2sideshield2 = new G4Trd("P2sideshield1", P2sideshieldW1, P2sideshieldW2, sideshieldTh, sideshieldTh, P2sideshieldL);

  G4LogicalVolume *P2sideshield1_log = new G4LogicalVolume( P2sideshield1, GetMaterial("Aluminum"), "P2sideshield1_log" );
  G4LogicalVolume *P2sideshield2_log = new G4LogicalVolume( P2sideshield2, GetMaterial("Aluminum"), "P2sideshield2_log" );
  
  G4RotationMatrix *P2rot1_temp = new G4RotationMatrix;
  P2rot1_temp->rotateZ(+90*deg);
  P2rot1_temp->rotateX(-sideshieldA);
  G4RotationMatrix *P2rot2_temp = new G4RotationMatrix;
  P2rot2_temp->rotateZ(-90*deg);
  P2rot2_temp->rotateX(-sideshieldA);

  new G4PVPlacement( P2rot1_temp, G4ThreeVector(-P2sideshield_xoffset-P2sideshieldL*sin(sideshieldA), 0, P2initPlacement_z+P2sideshield_zoffset+P2sideshieldL*cos(sideshieldA)), P2sideshield1_log, "P2sideshield1_log", worldlog, false, 0, ChkOverlaps);
  P2sideshield1_log->SetVisAttributes(Aluminum);
 
  new G4PVPlacement( P2rot2_temp, G4ThreeVector(P2sideshield_xoffset+P2sideshieldL*sin(sideshieldA), 0, P2initPlacement_z+P2sideshield_zoffset+P2sideshieldL*cos(sideshieldA)), P2sideshield2_log, "P2sideshield2_log", worldlog, false, 0, ChkOverlaps);
  P2sideshield2_log->SetVisAttributes(Aluminum);

 
  /*

  //G4Trap *P2sideshield1 = new G4Trap("P2sideshield1", P2sideshieldL, P2sideshieldTh, P2sideshieldW2, P2sideshieldW1);
  G4Trap *P2sideshield2 = new G4Trap("P2sideshield2", P2sideshieldL, P2sideshieldTh, P2sideshieldW2, P2sideshieldW1);
  
  G4Trap *P2sideshield1 = new G4Trap("P2sideshield1", P2sideshieldTh, P2sideshieldL, P2sideshieldW1/50*inch, P2sideshieldW2/50*inch); //Not sure why these factors of 50 are necessary yet, these values are guesses
  //G4Trap *P2sideshield12 = new G4Trap("P2sideshield12", P2sideshieldTh, P2sideshieldL, P2sideshieldW2, P2sideshieldW1);
 
  //G4Trap *P2sideshield1 = new G4Trap("P2sideshield1", 20.*inch, 20.*inch, 5.*inch, 5.*inch, P2sideshieldL);

  G4LogicalVolume *P2sideshield1_log = new G4LogicalVolume( P2sideshield1, GetMaterial("Aluminum"), "P2sideshield1_log" );
  G4LogicalVolume *P2sideshield2_log = new G4LogicalVolume( P2sideshield2, GetMaterial("Aluminum"), "P2sideshield2_log" );

  //G4RotationMatrix *P2rot1_temp = new G4RotationMatrix;
  //P2rot1_temp->rotateZ(+90*deg);
  //P2rot1_temp->rotateX(+2*deg);
  //G4RotationMatrix *P2rot2_temp = new G4RotationMatrix;
  //P2rot2_temp->rotateZ(+90*deg);
  //P2rot2_temp->rotateX(-2*deg);
  G4RotationMatrix *P2rot1_temp = new G4RotationMatrix;
  P2rot1_temp->rotateZ(+90.*deg);
  P2rot1_temp->rotateX(-91.5*deg);
  G4RotationMatrix *P2rot2_temp = new G4RotationMatrix;
  P2rot2_temp->rotateZ(-90.*deg);
  P2rot2_temp->rotateX(-88.5*deg);
  G4RotationMatrix *P2rot3_temp = new G4RotationMatrix;
  P2rot3_temp->rotateZ(+90.*deg);
  P2rot3_temp->rotateX(-88.5*deg);
  G4RotationMatrix *P2rot4_temp = new G4RotationMatrix;
  P2rot4_temp->rotateZ(-90.*deg);
  P2rot4_temp->rotateX(-91.5*deg);

  
  //Assuming aluminum sideshields - Assembling a single trapezoidal shield from two G4Traps where three sides of the plane are at right angles and the last is defined by W1 and W2. The y offsets are guesses at this stage and will follow when the factors of 50 in W1 and W2 are understood.

  //new G4PVPlacement( P2rot2_temp, G4ThreeVector(-(P2sideshield_xoffset), 0, P2initPlacement_z+P2_offset2+P2sideshieldL/2), P2sideshield2_log, "P2sideshield2_log", worldlog, false, 0, ChkOverlaps);
  //P2sideshield2_log->SetVisAttributes(Aluminum);
  new G4PVPlacement( P2rot1_temp, G4ThreeVector(-P2sideshield_xoffset, -P2sideshieldW1/4, P2initPlacement_z+P2_offset2+P2sideshieldL/2), P2sideshield1_log, "P2sideshield1_log", worldlog, false, 0, ChkOverlaps);
  P2sideshield1_log->SetVisAttributes(Aluminum);
  new G4PVPlacement( P2rot2_temp, G4ThreeVector(-P2sideshield_xoffset, P2sideshieldW1/4, P2initPlacement_z+P2_offset2+P2sideshieldL/2), P2sideshield1_log, "P2sideshield1_log", worldlog, false, 0, ChkOverlaps);
  P2sideshield1_log->SetVisAttributes(Aluminum);
 
  //new G4PVPlacement( P2rot1_temp, G4ThreeVector(P2sideshield_xoffset, 0, P2initPlacement_z+P2_offset2+P2sideshieldL/2), P2sideshield1_log, "P2sideshield1_log", worldlog, false, 0, ChkOverlaps);
  //P2sideshield1_log->SetVisAttributes(Aluminum);
  new G4PVPlacement( P2rot3_temp, G4ThreeVector(P2sideshield_xoffset, -P2sideshieldW1/4, P2initPlacement_z+P2_offset2+P2sideshieldL/2), P2sideshield1_log, "P2sideshield1_log", worldlog, false, 0, ChkOverlaps);
  P2sideshield1_log->SetVisAttributes(Aluminum);
  new G4PVPlacement( P2rot4_temp, G4ThreeVector(P2sideshield_xoffset, P2sideshieldW1/4, P2initPlacement_z+P2_offset2+P2sideshieldL/2), P2sideshield1_log, "P2sideshield1_log", worldlog, false, 0, ChkOverlaps);
  P2sideshield1_log->SetVisAttributes(Aluminum);
 
  */

  //First Corrector Magnet and Fixtures

  //Visuals


  //SECTION THREE
  //This section includes all components from the first corrector magnet up to and including the second corrector magnet.
  //Specifications - ordered first from inside to out then from target-proximal to distant
  
  //dz 52.0*inch from beginning of first middle ring to the beginning of second, 72.0*inch from beginning of second middle ring to the beginning of third

  //G4double P2initPlacement_z = initPlacement_z+2*P1ringA_L+2*P1disk_L+2*P1tubeA_L+16*P1ringL+2*P1tubeB_L+2*P1tubeC_L+2*P1tubeD_L+2*P1ringE_L+2*P1tubeE_L+4*P1ringF_L+2*P1tubeF_L+4*P1ringH_L+2*P1ringJ_L;

  G4double P2P3displacement = 50.141*inch;
  G4double P3initPlacement_z = P2initPlacement_z+P2P3displacement;

  //Inner cone included in section one code

  //Rings General Specifications, 27 in total
  G4double P3ringTh = 0.5*inch;
  //G4double P3ringr0 = 1.895*inch; //Measured - ends up inside first cone geometry
  //G4double P3ringr0 = 3.837*inch;
  G4double P3ringr0 = P2ringr0+(P2P3displacement-P2_offset2)*tan(P2DAngle);

  G4double P3ringL = 1.625/2*inch;
  G4double P3ringSep = 0.375*inch;
  //P3DAngle = P2DAngle

  G4double P3_offset = 2*P3ringL+P3ringSep;

  //Middle cone - staged, wrapping the inner cone in the central section
  G4double P3coneB_rin1 = 3.257*inch;
  G4double P3coneB_rou1 = 3.507*inch;
  G4double P3coneB_rin2 = 4.556*inch;
  G4double P3coneB_rou2 = 4.807*inch;
  G4double P3coneB_L = 49.625/2*inch;
  
  G4Cons *P3coneB = new G4Cons("P3coneB", P3coneB_rin1, P3coneB_rou1, P3coneB_rin2, P3coneB_rou2, P3coneB_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P3coneBLog = new G4LogicalVolume(P3coneB, GetMaterial("Aluminum"), "P3coneB_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+P3_offset+P3coneB_L), P3coneBLog, "P3coneBLog_pv", worldlog, false, 0 , ChkOverlaps);

  //P3coneBLog->SetVisAttributes( DebugRed );
  P3coneBLog->SetVisAttributes(Aluminum);


  //Ring 1
  G4double P3ring1_rin1 = P3ringr0;
  G4double P3ring1_rin2 = P3ringr0+2*P3ringL*tan(P2DAngle);

  G4Cons *P3ring1 = new G4Cons("P3ring1", P3ring1_rin1, P3ring1_rin1+P3ringTh, P3ring1_rin2, P3ring1_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring1Log = new G4LogicalVolume(P3ring1, GetMaterial("Iron"), "P3ring1_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+P3ringL), P3ring1Log, "P3ring1Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring1Log->SetVisAttributes(Iron);

  //Ring 2
  G4double P3ring2_rin1 = P3ringr0+P3_offset*tan(P2DAngle);
  G4double P3ring2_rin2 = P3ringr0+(P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring2 = new G4Cons("P3ring2", P3ring2_rin1, P3ring2_rin1+P3ringTh, P3ring2_rin2, P3ring2_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring2Log = new G4LogicalVolume(P3ring2, GetMaterial("Iron"), "P3ring2_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+P3_offset+P3ringL), P3ring2Log, "P3ring2Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring2Log->SetVisAttributes(Iron);

  //Ring 3
  G4double P3ring3_rin1 = P3ringr0+2*P3_offset*tan(P2DAngle);
  G4double P3ring3_rin2 = P3ringr0+(2*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring3 = new G4Cons("P3ring3", P3ring3_rin1, P3ring3_rin1+P3ringTh, P3ring3_rin2, P3ring3_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring3Log = new G4LogicalVolume(P3ring3, GetMaterial("Iron"), "P3ring3_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+2*P3_offset+P3ringL), P3ring3Log, "P3ring3Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring3Log->SetVisAttributes(Iron);

  //Ring 4
  G4double P3ring4_rin1 = P3ringr0+3*P3_offset*tan(P2DAngle);
  G4double P3ring4_rin2 = P3ringr0+(3*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring4 = new G4Cons("P3ring4", P3ring4_rin1, P3ring4_rin1+P3ringTh, P3ring4_rin2, P3ring4_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring4Log = new G4LogicalVolume(P3ring4, GetMaterial("Iron"), "P3ring4_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+3*P3_offset+P3ringL), P3ring4Log, "P3ring4Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring4Log->SetVisAttributes(Iron);

  //Ring 5
  G4double P3ring5_rin1 = P3ringr0+4*P3_offset*tan(P2DAngle);
  G4double P3ring5_rin2 = P3ringr0+(4*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring5 = new G4Cons("P3ring5", P3ring5_rin1, P3ring5_rin1+P3ringTh, P3ring5_rin2, P3ring5_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring5Log = new G4LogicalVolume(P3ring5, GetMaterial("Iron"), "P3ring5_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+4*P3_offset+P3ringL), P3ring5Log, "P3ring5Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring5Log->SetVisAttributes(Iron);

  //Ring 6
  G4double P3ring6_rin1 = P3ringr0+5*P3_offset*tan(P2DAngle);
  G4double P3ring6_rin2 = P3ringr0+(5*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring6 = new G4Cons("P3ring6", P3ring6_rin1, P3ring6_rin1+P3ringTh, P3ring6_rin2, P3ring6_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring6Log = new G4LogicalVolume(P3ring6, GetMaterial("Iron"), "P3ring6_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+5*P3_offset+P3ringL), P3ring6Log, "P3ring6Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring6Log->SetVisAttributes(Iron);

  //Ring 7
  G4double P3ring7_rin1 = P3ringr0+6*P3_offset*tan(P2DAngle);
  G4double P3ring7_rin2 = P3ringr0+(6*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring7 = new G4Cons("P3ring7", P3ring7_rin1, P3ring7_rin1+P3ringTh, P3ring7_rin2, P3ring7_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring7Log = new G4LogicalVolume(P3ring7, GetMaterial("Iron"), "P3ring7_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+6*P3_offset+P3ringL), P3ring7Log, "P3ring7Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring7Log->SetVisAttributes(Iron);

  //Ring 8
  G4double P3ring8_rin1 = P3ringr0+7*P3_offset*tan(P2DAngle);
  G4double P3ring8_rin2 = P3ringr0+(7*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring8 = new G4Cons("P3ring8", P3ring8_rin1, P3ring8_rin1+P3ringTh, P3ring8_rin2, P3ring8_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring8Log = new G4LogicalVolume(P3ring8, GetMaterial("Iron"), "P3ring8_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+7*P3_offset+P3ringL), P3ring8Log, "P3ring8Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring8Log->SetVisAttributes(Iron);

  //Ring 9
  G4double P3ring9_rin1 = P3ringr0+8*P3_offset*tan(P2DAngle);
  G4double P3ring9_rin2 = P3ringr0+(8*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring9 = new G4Cons("P3ring9", P3ring9_rin1, P3ring9_rin1+P3ringTh, P3ring9_rin2, P3ring9_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring9Log = new G4LogicalVolume(P3ring9, GetMaterial("Iron"), "P3ring9_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+8*P3_offset+P3ringL), P3ring9Log, "P3ring9Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring9Log->SetVisAttributes(Iron);

  //Ring 10
  G4double P3ring10_rin1 = P3ringr0+9*P3_offset*tan(P2DAngle);
  G4double P3ring10_rin2 = P3ringr0+(9*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring10 = new G4Cons("P3ring10", P3ring10_rin1, P3ring10_rin1+P3ringTh, P3ring10_rin2, P3ring10_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring10Log = new G4LogicalVolume(P3ring10, GetMaterial("Iron"), "P3ring10_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+9*P3_offset+P3ringL), P3ring10Log, "P3ring10Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring10Log->SetVisAttributes(Iron);

  //Ring 11
  G4double P3ring11_rin1 = P3ringr0+10*P3_offset*tan(P2DAngle);
  G4double P3ring11_rin2 = P3ringr0+(10*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring11 = new G4Cons("P3ring11", P3ring11_rin1, P3ring11_rin1+P3ringTh, P3ring11_rin2, P3ring11_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring11Log = new G4LogicalVolume(P3ring11, GetMaterial("Iron"), "P3ring11_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+10*P3_offset+P3ringL), P3ring11Log, "P3ring11Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring11Log->SetVisAttributes(Iron);

  //Ring 12
  G4double P3ring12_rin1 = P3ringr0+11*P3_offset*tan(P2DAngle);
  G4double P3ring12_rin2 = P3ringr0+(11*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring12 = new G4Cons("P3ring12", P3ring12_rin1, P3ring12_rin1+P3ringTh, P3ring12_rin2, P3ring12_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring12Log = new G4LogicalVolume(P3ring12, GetMaterial("Iron"), "P3ring12_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+11*P3_offset+P3ringL), P3ring12Log, "P3ring12Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring12Log->SetVisAttributes(Iron);

  //Ring 13
  G4double P3ring13_rin1 = P3ringr0+12*P3_offset*tan(P2DAngle);
  G4double P3ring13_rin2 = P3ringr0+(12*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring13 = new G4Cons("P3ring13", P3ring13_rin1, P3ring13_rin1+P3ringTh, P3ring13_rin2, P3ring13_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring13Log = new G4LogicalVolume(P3ring13, GetMaterial("Iron"), "P3ring13_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+12*P3_offset+P3ringL), P3ring13Log, "P3ring13Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring13Log->SetVisAttributes(Iron);

  //Ring 14
  G4double P3ring14_rin1 = P3ringr0+13*P3_offset*tan(P2DAngle);
  G4double P3ring14_rin2 = P3ringr0+(13*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring14 = new G4Cons("P3ring14", P3ring14_rin1, P3ring14_rin1+P3ringTh, P3ring14_rin2, P3ring14_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring14Log = new G4LogicalVolume(P3ring14, GetMaterial("Iron"), "P3ring14_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+13*P3_offset+P3ringL), P3ring14Log, "P3ring14Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring14Log->SetVisAttributes(Iron);

  //Ring 15
  G4double P3ring15_rin1 = P3ringr0+14*P3_offset*tan(P2DAngle);
  G4double P3ring15_rin2 = P3ringr0+(14*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring15 = new G4Cons("P3ring15", P3ring15_rin1, P3ring15_rin1+P3ringTh, P3ring15_rin2, P3ring15_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring15Log = new G4LogicalVolume(P3ring15, GetMaterial("Iron"), "P3ring15_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+14*P3_offset+P3ringL), P3ring15Log, "P3ring15Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring15Log->SetVisAttributes(Iron);

  //Ring 16
  G4double P3ring16_rin1 = P3ringr0+15*P3_offset*tan(P2DAngle);
  G4double P3ring16_rin2 = P3ringr0+(15*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring16 = new G4Cons("P3ring16", P3ring16_rin1, P3ring16_rin1+P3ringTh, P3ring16_rin2, P3ring16_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring16Log = new G4LogicalVolume(P3ring16, GetMaterial("Iron"), "P3ring16_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+15*P3_offset+P3ringL), P3ring16Log, "P3ring16Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring16Log->SetVisAttributes(Iron);

  //Ring 17
  G4double P3ring17_rin1 = P3ringr0+16*P3_offset*tan(P2DAngle);
  G4double P3ring17_rin2 = P3ringr0+(16*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring17 = new G4Cons("P3ring17", P3ring17_rin1, P3ring17_rin1+P3ringTh, P3ring17_rin2, P3ring17_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring17Log = new G4LogicalVolume(P3ring17, GetMaterial("Iron"), "P3ring17_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+16*P3_offset+P3ringL), P3ring17Log, "P3ring17Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring17Log->SetVisAttributes(Iron);

  //Ring 18
  G4double P3ring18_rin1 = P3ringr0+17*P3_offset*tan(P2DAngle);
  G4double P3ring18_rin2 = P3ringr0+(17*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring18 = new G4Cons("P3ring18", P3ring18_rin1, P3ring18_rin1+P3ringTh, P3ring18_rin2, P3ring18_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring18Log = new G4LogicalVolume(P3ring18, GetMaterial("Iron"), "P3ring18_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+17*P3_offset+P3ringL), P3ring18Log, "P3ring18Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring18Log->SetVisAttributes(Iron);

  //Ring 19
  G4double P3ring19_rin1 = P3ringr0+18*P3_offset*tan(P2DAngle);
  G4double P3ring19_rin2 = P3ringr0+(18*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring19 = new G4Cons("P3ring19", P3ring19_rin1, P3ring19_rin1+P3ringTh, P3ring19_rin2, P3ring19_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring19Log = new G4LogicalVolume(P3ring19, GetMaterial("Iron"), "P3ring19_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+18*P3_offset+P3ringL), P3ring19Log, "P3ring19Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring19Log->SetVisAttributes(Iron);

  //Ring 20
  G4double P3ring20_rin1 = P3ringr0+19*P3_offset*tan(P2DAngle);
  G4double P3ring20_rin2 = P3ringr0+(19*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring20 = new G4Cons("P3ring20", P3ring20_rin1, P3ring20_rin1+P3ringTh, P3ring20_rin2, P3ring20_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring20Log = new G4LogicalVolume(P3ring20, GetMaterial("Iron"), "P3ring20_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+19*P3_offset+P3ringL), P3ring20Log, "P3ring20Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring20Log->SetVisAttributes(Iron);

  //Ring 21
  G4double P3ring21_rin1 = P3ringr0+20*P3_offset*tan(P2DAngle);
  G4double P3ring21_rin2 = P3ringr0+(20*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring21 = new G4Cons("P3ring21", P3ring21_rin1, P3ring21_rin1+P3ringTh, P3ring21_rin2, P3ring21_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring21Log = new G4LogicalVolume(P3ring21, GetMaterial("Iron"), "P3ring21_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+20*P3_offset+P3ringL), P3ring21Log, "P3ring21Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring21Log->SetVisAttributes(Iron);

  //Ring 22
  G4double P3ring22_rin1 = P3ringr0+21*P3_offset*tan(P2DAngle);
  G4double P3ring22_rin2 = P3ringr0+(21*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring22 = new G4Cons("P3ring22", P3ring22_rin1, P3ring22_rin1+P3ringTh, P3ring22_rin2, P3ring22_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring22Log = new G4LogicalVolume(P3ring22, GetMaterial("Iron"), "P3ring22_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+21*P3_offset+P3ringL), P3ring22Log, "P3ring22Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring22Log->SetVisAttributes(Iron);

  //Ring 23
  G4double P3ring23_rin1 = P3ringr0+22*P3_offset*tan(P2DAngle);
  G4double P3ring23_rin2 = P3ringr0+(22*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring23 = new G4Cons("P3ring23", P3ring23_rin1, P3ring23_rin1+P3ringTh, P3ring23_rin2, P3ring23_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring23Log = new G4LogicalVolume(P3ring23, GetMaterial("Iron"), "P3ring23_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+22*P3_offset+P3ringL), P3ring23Log, "P3ring23Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring23Log->SetVisAttributes(Iron);

  //Ring 24
  G4double P3ring24_rin1 = P3ringr0+23*P3_offset*tan(P2DAngle);
  G4double P3ring24_rin2 = P3ringr0+(23*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring24 = new G4Cons("P3ring24", P3ring24_rin1, P3ring24_rin1+P3ringTh, P3ring24_rin2, P3ring24_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring24Log = new G4LogicalVolume(P3ring24, GetMaterial("Iron"), "P3ring24_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+23*P3_offset+P3ringL), P3ring24Log, "P3ring24Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring24Log->SetVisAttributes(Iron);

  //Ring 25
  G4double P3ring25_rin1 = P3ringr0+24*P3_offset*tan(P2DAngle);
  G4double P3ring25_rin2 = P3ringr0+(24*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring25 = new G4Cons("P3ring25", P3ring25_rin1, P3ring25_rin1+P3ringTh, P3ring25_rin2, P3ring25_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring25Log = new G4LogicalVolume(P3ring25, GetMaterial("Iron"), "P3ring25_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+24*P3_offset+P3ringL), P3ring25Log, "P3ring25Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring25Log->SetVisAttributes(Iron);

  //Ring 26
  G4double P3ring26_rin1 = P3ringr0+25*P3_offset*tan(P2DAngle);
  G4double P3ring26_rin2 = P3ringr0+(25*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring26 = new G4Cons("P3ring26", P3ring26_rin1, P3ring26_rin1+P3ringTh, P3ring26_rin2, P3ring26_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring26Log = new G4LogicalVolume(P3ring26, GetMaterial("Iron"), "P3ring26_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+25*P3_offset+P3ringL), P3ring26Log, "P3ring26Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring26Log->SetVisAttributes(Iron);

  //Ring 27
  G4double P3ring27_rin1 = P3ringr0+26*P3_offset*tan(P2DAngle);
  G4double P3ring27_rin2 = P3ringr0+(26*P3_offset+P3ringL)*tan(P2DAngle);

  G4Cons *P3ring27 = new G4Cons("P3ring27", P3ring27_rin1, P3ring27_rin1+P3ringTh, P3ring27_rin2, P3ring27_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P3ring27Log = new G4LogicalVolume(P3ring27, GetMaterial("Iron"), "P3ring27_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+26*P3_offset+P3ringL), P3ring27Log, "P3ring27Log_pv", worldlog, false, 0 , ChkOverlaps);

  P3ring27Log->SetVisAttributes(Iron);

  //Still clunky - more efficient method with loop, but don't know how make file parses vector of G4 objects
  /*
  G4double P3ring_rin1[] = {};
  G4double P3ring_rin2[] = {};
  
  G4Cons P3ring[] = {};

  G4LogicalVolume *P3ringLog[] = {};

  G4PVPlacement *P3ringLogpv[] = {};

  
  for (int i=0, i<27){

  P3ring_rin1[i] = P3ringr0+i*P3_offset*tan(P2DAngle);
  P3ring_rin2[i] = P3ringr0+(i*P3_offset+P3ringL)*tan(P2DAngle);

  string P3label1text = "P3ring" + to_string([i]);
  string P3label2text = "P3ringLog" + to_string([i]);
  string P3label3text = "P3ringLogPV" + to_string([i]);

  P3ring[i] = new G4Cons( to_string(P3label1text) , P3ring_rin1[i], P3ring_rin1[i]+P3ringTh, P3ring_rin2[i], P3ring_rin2[i]+P3ringTh, P3ringL, 0.*deg, 360.*deg);

  P3ringLog[i] = new G4LogicalVolume( P3label1text, GetMaterial("Iron"), to_string(label2text), 0, 0, 0);

  P3ringLogpv[i] = new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+i*P3_offset+P3ringL), P3label2text, to_string(label3text), worldlog, false, 0 , ChkOverlaps);

  }
  */

  //Updated before this mark 6.16.20


  
  //P3 Side shields - NOT shields per material update 7.24.20
  //Specifications
  G4double P3sideshieldL = 57.813/2*inch;
  G4double P3sideshieldW1 = 6.547/2*inch;
  G4double P3sideshieldW2 = 8.688/2*inch;
  G4double P3sideshield_xoffset = 10.067/2*inch+sideshieldTh;
  G4double P3sideshield_zoffset = 48.308*inch;

  //Trapezoid class 
  G4Trd *P3sideshield1 = new G4Trd("P3sideshield1", P3sideshieldW1, P3sideshieldW2, sideshieldTh, sideshieldTh, P3sideshieldL);
  G4Trd *P3sideshield2 = new G4Trd("P3sideshield1", P3sideshieldW1, P3sideshieldW2, sideshieldTh, sideshieldTh, P3sideshieldL);

  G4LogicalVolume *P3sideshield1_log = new G4LogicalVolume( P3sideshield1, GetMaterial("Aluminum"), "P3sideshield1_log" );
  G4LogicalVolume *P3sideshield2_log = new G4LogicalVolume( P3sideshield2, GetMaterial("Aluminum"), "P3sideshield2_log" );
  
  G4RotationMatrix *P3rot1_temp = new G4RotationMatrix;
  P3rot1_temp->rotateZ(+90*deg);
  P3rot1_temp->rotateX(-sideshieldA);
  G4RotationMatrix *P3rot2_temp = new G4RotationMatrix;
  P3rot2_temp->rotateZ(-90*deg);
  P3rot2_temp->rotateX(-sideshieldA);

  new G4PVPlacement( P3rot1_temp, G4ThreeVector(-P3sideshield_xoffset-P3sideshieldL*sin(sideshieldA), 0, P2initPlacement_z+P3sideshield_zoffset+P3sideshieldL*cos(sideshieldA)), P3sideshield1_log, "P3sideshield1_log", worldlog, false, 0, ChkOverlaps);
  P3sideshield1_log->SetVisAttributes(Aluminum);
 
  new G4PVPlacement( P3rot2_temp, G4ThreeVector(P3sideshield_xoffset+P3sideshieldL*sin(sideshieldA), 0, P2initPlacement_z+P3sideshield_zoffset+P3sideshieldL*cos(sideshieldA)), P3sideshield2_log, "P3sideshield2_log", worldlog, false, 0, ChkOverlaps);
  P3sideshield2_log->SetVisAttributes(Aluminum);


  /*


  //G4Trap *P3sideshield1 = new G4Trap("P3sideshield1", P3sideshieldL, P3sideshieldTh, P3sideshieldW2, P3sideshieldW1);
  G4Trap *P3sideshield2 = new G4Trap("P3sideshield2", P3sideshieldL, P3sideshieldTh, P3sideshieldW2, P3sideshieldW1);
  G4Trap *P3sideshield1 = new G4Trap("P3sideshield1", P3sideshieldTh, P3sideshieldL, P3sideshieldW1/50*inch, P3sideshieldW2/50*inch); //Not sure why these factors of 50 are necessary yet, both values are guesses.

  G4LogicalVolume *P3sideshield1_log = new G4LogicalVolume( P3sideshield1, GetMaterial("Aluminum"), "P3sideshield1_log" );
  G4LogicalVolume *P3sideshield2_log = new G4LogicalVolume( P3sideshield2, GetMaterial("Aluminum"), "P3sideshield2_log" );

 G4RotationMatrix *P3rot1_temp = new G4RotationMatrix;
  P3rot1_temp->rotateZ(+90*deg);
  P3rot1_temp->rotateX(-91.5*deg);
  G4RotationMatrix *P3rot2_temp = new G4RotationMatrix;
  P3rot2_temp->rotateZ(-90*deg);
  P3rot2_temp->rotateX(-88.5*deg);
  G4RotationMatrix *P3rot3_temp = new G4RotationMatrix;
  P3rot3_temp->rotateZ(+90*deg);
  P3rot3_temp->rotateX(-88.5*deg);
  G4RotationMatrix *P3rot4_temp = new G4RotationMatrix;
  P3rot4_temp->rotateZ(-90*deg);
  P3rot4_temp->rotateX(-91.5*deg);

  
  //Assuming aluminum sideshields - Assembling a single trapezoidal shield from two G4Traps where three sides of the plane are at right angles and the last is defined by W1 and W2. The y offsets are guesses at this stage and will follow when the factors of 50 in W1 and W2 are understood.

  new G4PVPlacement( P3rot1_temp, G4ThreeVector(-P3sideshield_xoffset, -P3sideshieldW1/4, P2initPlacement_z+P3sideshield_zoffset+P3sideshieldL/2), P3sideshield1_log, "P3sideshield1_log", worldlog, false, 0, ChkOverlaps);
  P3sideshield1_log->SetVisAttributes(Aluminum);
  new G4PVPlacement( P3rot2_temp, G4ThreeVector(-P3sideshield_xoffset, P3sideshieldW1/4, P2initPlacement_z+P3sideshield_zoffset+P3sideshieldL/2), P3sideshield1_log, "P3sideshield1_log", worldlog, false, 0, ChkOverlaps);
  P3sideshield1_log->SetVisAttributes(Aluminum);
 
  new G4PVPlacement( P3rot3_temp, G4ThreeVector(P3sideshield_xoffset, -P3sideshieldW1/4, P2initPlacement_z+P3sideshield_zoffset+P3sideshieldL/2), P3sideshield1_log, "P3sideshield1_log", worldlog, false, 0, ChkOverlaps);
  P3sideshield1_log->SetVisAttributes(Aluminum);
  new G4PVPlacement( P3rot4_temp, G4ThreeVector(P3sideshield_xoffset, P3sideshieldW1/4, P2initPlacement_z+P3sideshield_zoffset+P3sideshieldL/2), P3sideshield1_log, "P3sideshield1_log", worldlog, false, 0, ChkOverlaps);
  P3sideshield1_log->SetVisAttributes(Aluminum);
 
  */


  //Second corrector magnet

  //Visuals


  

  //SECTION FOUR
  //This section includes all components from the second corrector magnet.
  //Specifications - ordered first from inside to out then from target-proximal to distant

  G4double P2P4displacement = 122.141*inch;
  G4double P4initPlacement_z = P2initPlacement_z+P2P4displacement;

  //Inner cone included in section one code

  //Rings General Specifications, 27 in total
  G4double P4ringTh = 0.5*inch;
  //G4double P4ringr0 = 5.722*inch;
  G4double P4ringr0 = P2ringr0+(P2P4displacement-P2_offset2)*tan(P2DAngle);
  G4double P4ringL = 1.625/2*inch;
  //P4DAngle = P2DAngle
  G4double P4ringSep = 0.375*inch;
  G4double P4_offset = 2*P4ringL+P4ringSep;

  //Middle cone - staged, wrapping the inner cone in the central section
  G4double P4coneB_rin1 = 5.142*inch;
  G4double P4coneB_rou1 = 5.392*inch;
  G4double P4coneB_rin2 = 5.446*inch;
  G4double P4coneB_rou2 = 5.696*inch;
  G4double P4coneB_L = 11.591/2*inch;
  
  G4Cons *P4coneB = new G4Cons("P4coneB", P4coneB_rin1, P4coneB_rou1, P4coneB_rin2, P4coneB_rou2, P4coneB_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P4coneBLog = new G4LogicalVolume(P4coneB, GetMaterial("Aluminum"), "P4coneB_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P4initPlacement_z+P4_offset+P4coneB_L), P4coneBLog, "P4coneBLog_pv", worldlog, false, 0 , ChkOverlaps);

  //P4coneBLog->SetVisAttributes( DebugRed );
  P4coneBLog->SetVisAttributes(Aluminum);

  //Ring 1
  G4double P4ring1_rin1 = P4ringr0;
  G4double P4ring1_rin2 = P4ringr0+2*P4ringL*tan(P2DAngle);

  G4Cons *P4ring1 = new G4Cons("P4ring1", P4ring1_rin1, P4ring1_rin1+P4ringTh, P4ring1_rin2, P4ring1_rin2+P4ringTh, P4ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P4ring1Log = new G4LogicalVolume(P4ring1, GetMaterial("Iron"), "P4ring1_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P4initPlacement_z+P4ringL), P4ring1Log, "P4ring1Log_pv", worldlog, false, 0 , ChkOverlaps);

  P4ring1Log->SetVisAttributes(Iron);

  //Ring 2
  G4double P4ring2_rin1 = P4ringr0+P4_offset*tan(P2DAngle);
  G4double P4ring2_rin2 = P4ringr0+(P4_offset+P4ringL)*tan(P2DAngle);

  G4Cons *P4ring2 = new G4Cons("P4ring2", P4ring2_rin1, P4ring2_rin1+P4ringTh, P4ring2_rin2, P4ring2_rin2+P4ringTh, P4ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P4ring2Log = new G4LogicalVolume(P4ring2, GetMaterial("Iron"), "P4ring2_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P4initPlacement_z+P4_offset+P4ringL), P4ring2Log, "P4ring2Log_pv", worldlog, false, 0 , ChkOverlaps);

  P4ring2Log->SetVisAttributes(Iron);

  //Ring 3
  G4double P4ring3_rin1 = P4ringr0+2*P4_offset*tan(P2DAngle);
  G4double P4ring3_rin2 = P4ringr0+(2*P4_offset+P4ringL)*tan(P2DAngle);

  G4Cons *P4ring3 = new G4Cons("P4ring3", P4ring3_rin1, P4ring3_rin1+P4ringTh, P4ring3_rin2, P4ring3_rin2+P4ringTh, P4ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P4ring3Log = new G4LogicalVolume(P4ring3, GetMaterial("Iron"), "P4ring3_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P4initPlacement_z+2*P4_offset+P4ringL), P4ring3Log, "P4ring3Log_pv", worldlog, false, 0 , ChkOverlaps);

  P4ring3Log->SetVisAttributes(Iron);

  //Ring 4
  G4double P4ring4_rin1 = P4ringr0+3*P4_offset*tan(P2DAngle);
  G4double P4ring4_rin2 = P4ringr0+(3*P4_offset+P4ringL)*tan(P2DAngle);

  G4Cons *P4ring4 = new G4Cons("P4ring4", P4ring4_rin1, P4ring4_rin1+P4ringTh, P4ring4_rin2, P4ring4_rin2+P4ringTh, P4ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P4ring4Log = new G4LogicalVolume(P4ring4, GetMaterial("Iron"), "P4ring4_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P4initPlacement_z+3*P4_offset+P4ringL), P4ring4Log, "P4ring4Log_pv", worldlog, false, 0 , ChkOverlaps);

  P4ring4Log->SetVisAttributes(Iron);

  //Ring 5
  G4double P4ring5_rin1 = P4ringr0+4*P4_offset*tan(P2DAngle);
  G4double P4ring5_rin2 = P4ringr0+(4*P4_offset+P4ringL)*tan(P2DAngle);

  G4Cons *P4ring5 = new G4Cons("P4ring5", P4ring5_rin1, P4ring5_rin1+P4ringTh, P4ring5_rin2, P4ring5_rin2+P4ringTh, P4ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P4ring5Log = new G4LogicalVolume(P4ring5, GetMaterial("Iron"), "P4ring5_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P4initPlacement_z+4*P4_offset+P4ringL), P4ring5Log, "P4ring5Log_pv", worldlog, false, 0 , ChkOverlaps);

  P4ring5Log->SetVisAttributes(Iron);

  //Ring 6
  G4double P4ring6_rin1 = P4ringr0+5*P4_offset*tan(P2DAngle);
  G4double P4ring6_rin2 = P4ringr0+(5*P4_offset+P4ringL)*tan(P2DAngle);

  G4Cons *P4ring6 = new G4Cons("P4ring6", P4ring6_rin1, P4ring6_rin1+P4ringTh, P4ring6_rin2, P4ring6_rin2+P4ringTh, P4ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P4ring6Log = new G4LogicalVolume(P4ring6, GetMaterial("Iron"), "P4ring6_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P4initPlacement_z+5*P4_offset+P4ringL), P4ring6Log, "P4ring6Log_pv", worldlog, false, 0 , ChkOverlaps);

  P4ring6Log->SetVisAttributes(Iron);

  //Ring 7
  G4double P4ring7_rin1 = P4ringr0+6*P4_offset*tan(P2DAngle);
  G4double P4ring7_rin2 = P4ringr0+(6*P4_offset+P4ringL)*tan(P2DAngle);

  G4Cons *P4ring7 = new G4Cons("P4ring7", P4ring7_rin1, P4ring7_rin1+P4ringTh, P4ring7_rin2, P4ring7_rin2+P4ringTh, P4ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P4ring7Log = new G4LogicalVolume(P4ring7, GetMaterial("Iron"), "P4ring7_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P4initPlacement_z+6*P4_offset+P4ringL), P4ring7Log, "P4ring7Log_pv", worldlog, false, 0 , ChkOverlaps);

  P4ring7Log->SetVisAttributes(Iron);


  //P4 Side shields - NOT shields per material update 7.24.20
  //Specifications
  G4double P4sideshieldL = 13.498/2*inch;
  G4double P4sideshieldW1 = 9.280/2*inch;
  G4double P4sideshieldW2 = 9.780/2*inch;
  G4double P4sideshield_xoffset = 12.8/2*inch+sideshieldTh;
  G4double P4sideshield_zoffset = 122.113*inch;
 
  //Trapezoid class
  G4Trd *P4sideshield1 = new G4Trd("P4sideshield1", P4sideshieldW1, P4sideshieldW2, sideshieldTh, sideshieldTh, P4sideshieldL);
  G4Trd *P4sideshield2 = new G4Trd("P4sideshield1", P4sideshieldW1, P4sideshieldW2, sideshieldTh, sideshieldTh, P4sideshieldL);

  G4LogicalVolume *P4sideshield1_log = new G4LogicalVolume( P4sideshield1, GetMaterial("Aluminum"), "P4sideshield1_log" );
  G4LogicalVolume *P4sideshield2_log = new G4LogicalVolume( P4sideshield2, GetMaterial("Aluminum"), "P4sideshield2_log" );
  
  G4RotationMatrix *P4rot1_temp = new G4RotationMatrix;
  P4rot1_temp->rotateZ(+90*deg);
  P4rot1_temp->rotateX(-sideshieldA);
  G4RotationMatrix *P4rot2_temp = new G4RotationMatrix;
  P4rot2_temp->rotateZ(-90*deg);
  P4rot2_temp->rotateX(-sideshieldA);

  new G4PVPlacement( P4rot1_temp, G4ThreeVector(-P4sideshield_xoffset-P4sideshieldL*sin(sideshieldA), 0, P2initPlacement_z+P4sideshield_zoffset+P4sideshieldL*cos(sideshieldA)), P4sideshield1_log, "P4sideshield1_log", worldlog, false, 0, ChkOverlaps);
  P4sideshield1_log->SetVisAttributes(Aluminum);
 
  new G4PVPlacement( P4rot2_temp, G4ThreeVector(P4sideshield_xoffset+P4sideshieldL*sin(sideshieldA), 0, P2initPlacement_z+P4sideshield_zoffset+P4sideshieldL*cos(sideshieldA)), P4sideshield2_log, "P4sideshield2_log", worldlog, false, 0, ChkOverlaps);
  P4sideshield2_log->SetVisAttributes(Aluminum);


  //P5 Endcap

  //General Specifications
  G4double P5initPlacement_z = P2initPlacement_z + 136.872*inch;
  
  //Ring A
  G4double P5ringA_rin = 10.734/2*inch-0.6208*inch*tan(P2DAngle); //Backing out the inner radius of the cone at interface to define ring there as an approximation.
  G4double P5ringA_rou = 13.940/2*inch;
  G4double P5ringA_L = 1.120*inch;

  G4Tubs *P5ringA = new G4Tubs("P5ringA", P5ringA_rin, P5ringA_rou, P5ringA_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringALog = new G4LogicalVolume(P5ringA, GetMaterial("Aluminum"), "P5ringA_log", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5initPlacement_z+P5ringA_L), P5ringALog, "P5ringALog_pv", worldlog, false, 0, ChkOverlaps);

  P5ringALog->SetVisAttributes(Aluminum);

  //Ring A Vacuum
  G4Tubs *P5ringA_vac = new G4Tubs("P5ringA_vac", 0.0, P5ringA_rin, P5ringA_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringA_vacLog = new G4LogicalVolume(P5ringA_vac, GetMaterial("Vacuum"), "P5ringA_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P5initPlacement_z+P5ringA_L), P5ringA_vacLog, "P5ringA_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  P5ringA_vacLog->SetVisAttributes( G4VisAttributes::Invisible);

  //Ring B
  G4double P5ringB_rin = 11.750/2*inch; 
  G4double P5ringB_rou = 14.0/2*inch;
  G4double P5ringB_L = 1.120*inch;

  G4Tubs *P5ringB = new G4Tubs("P5ringB", P5ringB_rin, P5ringB_rou, P5ringB_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringBLog = new G4LogicalVolume(P5ringB, GetMaterial("Aluminum"), "P5ringB_log", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5initPlacement_z+2*P5ringA_L+P5ringB_L), P5ringBLog, "P5ringBLog_pv", worldlog, false, 0, ChkOverlaps);

  P5ringBLog->SetVisAttributes(Aluminum);

  //Ring B Vacuum
  G4Tubs *P5ringB_vac = new G4Tubs("P5ringB_vac", 0.0, P5ringB_rin, P5ringB_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringB_vacLog = new G4LogicalVolume(P5ringB_vac, GetMaterial("Vacuum"), "P5ringB_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P5initPlacement_z+2*P5ringB_L+P5ringB_L), P5ringB_vacLog, "P5ringB_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  P5ringB_vacLog->SetVisAttributes( G4VisAttributes::Invisible);

  //Ring C
  G4double P5ringC_rin = 11.750/2*inch; 
  G4double P5ringC_rou = 12.0/2*inch;
  G4double P5ringC_L = 0.405*inch;

  G4Tubs *P5ringC = new G4Tubs("P5ringC", P5ringC_rin, P5ringC_rou, P5ringC_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringCLog = new G4LogicalVolume(P5ringC, GetMaterial("Aluminum"), "P5ringC_log", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5initPlacement_z+2*P5ringA_L+2*P5ringB_L+P5ringC_L), P5ringCLog, "P5ringCLog_pv", worldlog, false, 0, ChkOverlaps);

  P5ringCLog->SetVisAttributes(Aluminum);

  //Ring C Vacuum
  G4Tubs *P5ringC_vac = new G4Tubs("P5ringC_vac", 0.0, P5ringC_rin, P5ringC_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringC_vacLog = new G4LogicalVolume(P5ringC_vac, GetMaterial("Vacuum"), "P5ringC_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P5initPlacement_z+2*P5ringA_L+2*P5ringB_L+P5ringC_L), P5ringC_vacLog, "P5ringC_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  P5ringC_vacLog->SetVisAttributes( G4VisAttributes::Invisible);

  //Ring D
  G4double P5ringD_rin = 12.710/2*inch; 
  G4double P5ringD_rou = 12.750/2*inch;
  G4double P5ringD_L = 2.487*inch;

  G4Tubs *P5ringD = new G4Tubs("P5ringD", P5ringD_rin, P5ringD_rou, P5ringD_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringDLog = new G4LogicalVolume(P5ringD, GetMaterial("Aluminum"), "P5ringD_log", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5initPlacement_z+2*P5ringA_L+2*P5ringB_L+2*P5ringC_L+P5ringD_L), P5ringDLog, "P5ringDLog_pv", worldlog, false, 0, ChkOverlaps);

  P5ringDLog->SetVisAttributes(Aluminum);

  //Ring D Vacuum
  G4Tubs *P5ringD_vac = new G4Tubs("P5ringD_vac", 0.0, P5ringD_rin, P5ringD_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringD_vacLog = new G4LogicalVolume(P5ringD_vac, GetMaterial("Vacuum"), "P5ringD_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P5initPlacement_z+2*P5ringA_L+2*P5ringB_L+2*P5ringC_L+P5ringD_L), P5ringD_vacLog, "P5ringD_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  P5ringD_vacLog->SetVisAttributes( G4VisAttributes::Invisible);

  //Ring E
  G4double P5ringE_rin = 11.750/2*inch; 
  G4double P5ringE_rou = 12.0/2*inch;
  G4double P5ringE_L = 0.380*inch;

  G4Tubs *P5ringE = new G4Tubs("P5ringE", P5ringE_rin, P5ringE_rou, P5ringE_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringELog = new G4LogicalVolume(P5ringE, GetMaterial("Aluminum"), "P5ringE_log", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5initPlacement_z+2*P5ringA_L+2*P5ringB_L+2*P5ringC_L+2*P5ringD_L+P5ringE_L), P5ringELog, "P5ringELog_pv", worldlog, false, 0, ChkOverlaps);

  P5ringELog->SetVisAttributes(Aluminum);

  //Ring E Vacuum
  G4Tubs *P5ringE_vac = new G4Tubs("P5ringE_vac", 0.0, P5ringE_rin, P5ringE_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringE_vacLog = new G4LogicalVolume(P5ringE_vac, GetMaterial("Vacuum"), "P5ringE_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P5initPlacement_z+2*P5ringA_L+2*P5ringB_L+2*P5ringC_L+2*P5ringD_L+P5ringE_L), P5ringE_vacLog, "P5ringE_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  P5ringE_vacLog->SetVisAttributes( G4VisAttributes::Invisible);

  //Ring F
  G4double P5ringF_rin = 11.750/2*inch; 
  G4double P5ringF_rou = 13.960/2*inch;
  G4double P5ringF_L = 1.120*inch;

  G4Tubs *P5ringF = new G4Tubs("P5ringF", P5ringF_rin, P5ringF_rou, P5ringF_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringFLog = new G4LogicalVolume(P5ringF, GetMaterial("Aluminum"), "P5ringF_log", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5initPlacement_z+2*P5ringA_L+2*P5ringB_L+2*P5ringC_L+2*P5ringD_L+2*P5ringE_L+P5ringF_L), P5ringFLog, "P5ringFLog_pv", worldlog, false, 0, ChkOverlaps);

  P5ringFLog->SetVisAttributes(Aluminum);

  //Ring F Vacuum
  G4Tubs *P5ringF_vac = new G4Tubs("P5ringF_vac", 0.0, P5ringF_rin, P5ringF_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringF_vacLog = new G4LogicalVolume(P5ringF_vac, GetMaterial("Vacuum"), "P5ringF_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P5initPlacement_z+2*P5ringA_L+2*P5ringB_L+2*P5ringC_L+2*P5ringD_L+2*P5ringE_L+P5ringF_L), P5ringF_vacLog, "P5ringF_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  P5ringF_vacLog->SetVisAttributes( G4VisAttributes::Invisible);


  //Welded bellows and extended beamline out to dump from commonexitbeamline placed in location as described there. Will need to check position relative to target and mating with modified beamline geometry described from P1 - P4 at P4 end. 

  //Visuals
  G4VisAttributes *ironColor= new G4VisAttributes(G4Colour(0.3,0.3,0.3));
  //G4VisAttributes *AlColor= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
  G4VisAttributes *Vacuum_visatt = new G4VisAttributes(G4Colour(0.1, 0.5, 0.9 ) );
  Vacuum_visatt->SetVisibility(false);
  G4VisAttributes *SteelColor = new G4VisAttributes( G4Colour( 0.75, 0.75, 0.75 ) );
  G4VisAttributes *CopperColor = new G4VisAttributes( G4Colour( 0.7, 0.3, 0.3 ) );
  
  G4double TargetCenter_zoffset = 6.50*inch;
  
  G4double z_formed_bellows = 52.440*inch - TargetCenter_zoffset; //relative to "target center"? or "origin"?
  G4double z_spool_piece = 58.44*inch - TargetCenter_zoffset;
  if(fDetCon->fBeamlineConf>2)z_spool_piece = 27.903*inch;
  G4double z_conic_vacline_weldment = 62.8*inch - TargetCenter_zoffset;
  G4double z_outer_magnetic = 71.782*inch - TargetCenter_zoffset;
  G4double z_inner_magnetic = 73.782*inch - TargetCenter_zoffset;
  G4double z_welded_bellows = 201.632*inch - TargetCenter_zoffset;

  G4double X=0.0, Y=0.0, Z=0.0;
  G4ThreeVector zero(0.0, 0.0, 0.0);
  
  G4double dz_welded_bellows = 207.144*inch - z_welded_bellows - TargetCenter_zoffset; // = =5.512 inches
  
  G4double Rin = 11.750/2.0*inch;
  G4double Rout = 14.0/2.0*inch;
  G4double Thick = 1.12*inch;
  
  G4Tubs *WB_Flange = new G4Tubs( "WB_Flange", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *WB_Flange_log = new G4LogicalVolume( WB_Flange, GetMaterial("Stainless_Steel"), "WB_Flange_log" );

  WB_Flange_log->SetVisAttributes( SteelColor );
    
  Z = z_welded_bellows + Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Flange_log, "WB_Flange1_phys", worldlog, false, 0 , ChkOverlaps );

  Z = z_welded_bellows + dz_welded_bellows - Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Flange_log, "WB_Flange2_phys", worldlog, false, 1 , ChkOverlaps );
  
  Rout = Rin + 0.125*inch;
  Thick = dz_welded_bellows - 2*1.12*inch;
  G4Tubs *WB_Bellows = new G4Tubs( "WB_Bellows", Rin, Rout, Thick/2.0, 0.0, twopi );
  G4LogicalVolume *WB_Bellows_log = new G4LogicalVolume(WB_Bellows, GetMaterial("Stainless_Steel"), "WB_Bellows_log" );

  WB_Bellows_log->SetVisAttributes( SteelColor );
  
  Z = z_welded_bellows + 1.12*inch + Thick/2.0;
  
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Bellows_log, "WB_Bellows_phys", worldlog, false, 0 , ChkOverlaps );
  
  Rin = 0.0;
  Rout = 11.750/2.0*inch;
  Thick = dz_welded_bellows;
  G4Tubs *WB_Vacuum = new G4Tubs( "WB_Vacuum", Rin, Rout, Thick/2.0, 0.0, twopi );

  G4LogicalVolume *WB_Vacuum_log = new G4LogicalVolume(WB_Vacuum, GetMaterial("Vacuum"), "WB_Vacuum_log" );

  WB_Vacuum_log->SetVisAttributes( Vacuum_visatt );
  
  Z = z_welded_bellows + dz_welded_bellows/2.0;

  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Vacuum_log, "WB_Vacuum_phys", worldlog, false, 0 , ChkOverlaps );
  
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

  TBL9_log->SetVisAttributes( Aluminum );
  TVL9_log->SetVisAttributes( Vacuum_visatt );
  
  // Then place the vacuum inside the Al Tube
  //Z = (159.51 + 352.636 - 2.84*0.5 + 4.237*2.54 + 41.0*2.54*0.5)*cm;
  Z = 207.144*inch + tDzz - TargetCenter_zoffset;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TVL9_log, "Extended_Vac1", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TBL9_log, "Extended_Al1", worldlog, false, 0 , ChkOverlaps );

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

  TML9_log->SetVisAttributes( Aluminum );
  TMV9_log->SetVisAttributes( Vacuum_visatt );
  // Then place vacuum inside of Al tube
  //Z = (159.51 + 352.636 - 2.84*0.5 + 4.237*2.54 + 41.0*2.54 + 0.5*217.0*2.54)*cm;
  Z = 207.144*inch + 41.0*inch + tDzz - TargetCenter_zoffset;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TMV9_log, "Extended_Vac2", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), TML9_log, "Extended_Al2", worldlog, false, 0 , ChkOverlaps );

  // For CPU speed, extend vacuum all the way to the edge of the "world" volume, so that we don't track beam electrons in air beyond interesting region.
  G4double Zstop = 30.0*m;
  G4double Zstart = Z + tDzz;

  std::cout << "**************** z = " << 0.5*(Zstop+Zstart)/cm << " cm" << std::endl;  
  G4double dz_df = 15.0*cm; 
  G4VisAttributes *Vacuum_visatt_df = new G4VisAttributes();
  Vacuum_visatt_df->SetColour( G4Colour::White() ); 
  Vacuum_visatt_df->SetForceWireframe(true); 
 
  G4double Zwidth = (Zstop-Zstart) - dz_df;
  G4Tubs *FakeVacuumExtension = new G4Tubs( "FakeVacuumExtension", tRmin, tRmax, Zwidth/2.0, tSPhi, tDphi );
  G4LogicalVolume *FakeVacuumExtension_log = new G4LogicalVolume( FakeVacuumExtension, GetMaterial("Vacuum"), "FakeVacuumExtension_log" );
  FakeVacuumExtension_log->SetVisAttributes( Vacuum_visatt_df );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0.5*(Zstop+Zstart)+dz_df), FakeVacuumExtension_log, "FakeVacuumExtension_phys", worldlog,false,0 , ChkOverlaps);





  
  //Ripping both corrector magnets code from common exit beamline. Foregoing placement vector Z component - will define manually with measurements from cad file.




  //Next, corrector magnets:
  // Define some dimensions that are going to be useful to define the distances

  //G4double TargetCenter_zoffset = 6.50*inch;
  
  //G4double z_formed_bellows = 52.440*inch - TargetCenter_zoffset; //relative to "target center"? or "origin"?

  z_formed_bellows = P2initPlacement_z+39.616*inch;
  G4double Bellows1L = 6.299/2*inch;
  G4double Bellows2L = 15.748/2*inch;

  G4double UpstreamCoilThickY = 1.68*inch;
  G4double UpstreamCoilThickX = 3.46*inch;
  //G4double UpstreamCoilWidth = 3.46*inch;
  G4double UpstreamCoilHeight = 8.17*inch;
  G4double UpstreamCoilDepth = 6.60*inch;
  //G4double UpstreamCoilDepth2 = 6.30*inch;
  G4double UpstreamCoilWidth = 7.56*inch;
    
  G4double YokeTopPiece_Width = 15.04*inch;
  G4double YokeTopPiece_Height = 3.94*inch;
  G4double YokeTopPiece_Depth = 6.30*inch;
  
  G4double YokeLeftPiece_Width = 2.76*inch;
  G4double YokeLeftPiece_Height = 11.81*inch + 2.0*YokeTopPiece_Height;
  G4double YokeLeftPiece_Depth = YokeTopPiece_Depth;
  
  G4double YokeRightNotchAngle = 18.43*deg;
  G4double YokeRightWidthFinal = YokeLeftPiece_Width;
  G4double YokeRightZFinal = YokeLeftPiece_Depth - 0.39*inch;
  G4double YokeRightWidthInitial = YokeRightWidthFinal - YokeRightZFinal*tan(YokeRightNotchAngle );
  
  G4double DownstreamTotalWidth = 17.58*inch + 2.0*2.76*inch;
  G4double DownstreamTotalHeight = 20.16*inch + 2.0*3.94*inch;
  G4double DownstreamYokeDepth = 15.75*inch;

  G4double DS_coil_depth = 8.91*inch;
  G4double DS_coil_height = 12.04*inch;
  G4double DS_coil_ThickX = 2.90*inch;
  G4double DS_coil_ThickY = 1.68*inch;
  
  //G4Box *UpstreamCoil_outerS = new G4Box("UpstreamCoil_outerS", UpstreamCoilThickX/2.0, (UpstreamCoilHeight +2.0*UpstreamCoilThickY)/2.0, (YokeTopPiece_Depth + 2.0*UpstreamCoilThickY)/2.0);
  G4Box *UpstreamCoil_outer = new G4Box("UpstreamCoil_outer", UpstreamCoilThickX/2.0, (UpstreamCoilHeight+2.0*UpstreamCoilThickY)/2.0, (UpstreamCoilDepth + 2.0*UpstreamCoilThickY)/2.0 );
  G4Box *UpstreamCoil_inner = new G4Box("UpstreamCoil_inner", UpstreamCoilThickX/2.0 + cm, UpstreamCoilHeight/2.0, UpstreamCoilDepth/2.0 );

  G4SubtractionSolid *UpstreamCoil = new G4SubtractionSolid( "UpstreamCoil", UpstreamCoil_outer, UpstreamCoil_inner, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *UpstreamCoil_log = new G4LogicalVolume(UpstreamCoil, GetMaterial("Copper"), "UpstreamCoil_log" );

  UpstreamCoil_log->SetVisAttributes( CopperColor );

  //Z = z_Magnets_array[0];//z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY;  
  //Z = z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY;
  //Z = 0.0;

  Z = z_formed_bellows+Bellows1L;
  X = (UpstreamCoilWidth+UpstreamCoilThickX)/2.0;
  Y = 0.0;

  //two placements of upstream coil:
  
  //new G4PVPlacement( 0, G4ThreeVector(-X, Y, Z), UpstreamCoil_log, "UpstreamCoil_phys_right", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), UpstreamCoil_log, "UpstreamCoil_phys_left", worldlog, false, 1  , ChkOverlaps);

  G4double UpstreamPoleDepth = 6.3*inch;
  G4double UpstreamPoleWidth = 4.02*inch;
  G4double UpstreamPoleHeight = 7.87*inch;
  //Next, make poles:
  G4Box *UpstreamPole = new G4Box( "UpstreamPole", UpstreamPoleWidth/2.0, UpstreamPoleHeight/2.0, UpstreamPoleDepth/2.0 );
  G4LogicalVolume *UpstreamPole_log = new G4LogicalVolume( UpstreamPole, GetMaterial("Iron"), "UpstreamPole_log" );
  UpstreamPole_log->SetVisAttributes( ironColor );
  //two placements of upstream poles:

  new G4PVPlacement( 0, G4ThreeVector(-X, Y, Z), UpstreamPole_log, "UpstreamPole_phys_right", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), UpstreamPole_log, "UpstreamPole_phys_left", worldlog, false, 1 , ChkOverlaps );

  //Next, make surrounding yoke:
  // G4double YokeTopPiece_Width = 15.04*inch;
  // G4double YokeTopPiece_Height = 3.94*inch;
  // G4double YokeTopPiece_Depth = 6.30*inch;

  G4Box *YokeTopPiece = new G4Box("YokeTopPiece", YokeTopPiece_Width/2.0, YokeTopPiece_Height/2.0, YokeTopPiece_Depth/2.0 );
  G4LogicalVolume *YokeTopPiece_log = new G4LogicalVolume( YokeTopPiece, GetMaterial("Iron"), "YokeTopPiece_log" );

  YokeTopPiece_log->SetVisAttributes( ironColor );
  
  X = 0.0;
  Y = (11.81*inch + YokeTopPiece_Height)/2.0;

  //two placements of yoke top piece (top and bottom symmetric):
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeTopPiece_log, "UpstreamYokeTop_phys", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X,-Y,Z), YokeTopPiece_log, "UpstreamYokeBottom_phys", worldlog, false, 1 , ChkOverlaps );

  // G4double YokeLeftPiece_Width = 2.76*inch;
  // G4double YokeLeftPiece_Height = 11.81*inch + 2.0*YokeTopPiece_Height;
  // G4double YokeLeftPiece_Depth = YokeTopPiece_Depth;
  
  G4Box *YokeLeftPiece = new G4Box("YokeLeftPiece", YokeLeftPiece_Width/2.0, YokeLeftPiece_Height/2.0, YokeLeftPiece_Depth/2.0 );
  G4LogicalVolume *YokeLeftPiece_log = new G4LogicalVolume( YokeLeftPiece, GetMaterial("Iron"), "YokeLeftPiece_log" );
  YokeLeftPiece_log->SetVisAttributes(ironColor );
  
  X = 7.52*inch + YokeLeftPiece_Width/2.0;
  Y = 0.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeLeftPiece_log, "UpstreamYokeLeftPiece_phys", worldlog, false, 0 , ChkOverlaps );

  // G4double YokeRightNotchAngle = 18.43*deg;
  // G4double YokeRightWidthFinal = YokeLeftPiece_Width;
  // G4double YokeRightZFinal = YokeLeftPiece_Depth - 0.39*inch;
  // G4double YokeRightWidthInitial = YokeRightWidthFinal - YokeRightZFinal*tan(YokeRightNotchAngle );

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
  //Z = z_Magnets_array[1];//z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0;
  //Z = z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0;
  Z = z_formed_bellows + Bellows1L;

  //Z = 0.0;

  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeRightPiece_log, "UpstreamYokeRightPiece_phys", worldlog, false, 0 , ChkOverlaps );

  //Downstream Corrector:
  // G4double DownstreamTotalWidth = 17.58*inch + 2.0*2.76*inch;
  // G4double DownstreamTotalHeight = 20.16*inch + 2.0*3.94*inch;
  // G4double DownstreamYokeDepth = 15.75*inch;

  G4double DownstreamYokeGapWidth = 17.58*inch;
  G4double DownstreamYokeGapHeight = 20.16*inch;
  G4Box *DownstreamYoke_box = new G4Box("DownstreamYoke_box", DownstreamTotalWidth/2.0, DownstreamTotalHeight/2.0, DownstreamYokeDepth/2.0 );
  G4Box *DownstreamYoke_gap = new G4Box("DownstreamYoke_gap", DownstreamYokeGapWidth/2.0, DownstreamYokeGapHeight/2.0, DownstreamYokeDepth/2.0+cm );
  G4SubtractionSolid *DownstreamYoke = new G4SubtractionSolid( "DownstreamYoke", DownstreamYoke_box, DownstreamYoke_gap, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *DownstreamYoke_log = new G4LogicalVolume( DownstreamYoke, GetMaterial("Iron"), "DownstreamYoke_log" );

  DownstreamYoke_log->SetVisAttributes( ironColor );

  z_formed_bellows = P2initPlacement_z + 108.479*inch;

  X = 0.0; Y = 0.0;
  //Z = z_Magnets_array[2];//z_formed_bellows + 76.09*inch + 1.71*inch + DownstreamYokeDepth/2.0;
  //Z = z_formed_bellows + 76.09*inch + 1.71*inch + DownstreamYokeDepth/2.0;

  Z = z_formed_bellows+Bellows2L;
  //Z = 0.0;


  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), DownstreamYoke_log, "DownstreamYoke_phys", worldlog, false, 0 , ChkOverlaps );

  // G4double DS_coil_depth = 8.91*inch;
  // G4double DS_coil_height = 12.04*inch;
  // G4double DS_coil_ThickX = 2.90*inch;
  // G4double DS_coil_ThickY = 1.68*inch;
  
  G4Box *DS_coil_outer = new G4Box( "DS_coil_outer", DS_coil_ThickX/2.0, (DS_coil_height + 2.0*DS_coil_ThickY)/2.0, (DS_coil_depth + 2.0*DS_coil_ThickY)/2.0 );
  G4Box *DS_coil_inner = new G4Box( "DS_coil_inner", DS_coil_ThickX/2.0+cm, DS_coil_height/2.0, DS_coil_depth/2.0 );

  G4SubtractionSolid *DS_coil = new G4SubtractionSolid( "DS_coil", DS_coil_outer, DS_coil_inner, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *DS_coil_log = new G4LogicalVolume( DS_coil, GetMaterial("Copper"), "DS_coil_log" );
  DS_coil_log->SetVisAttributes(CopperColor );
  
  X = 11.67*inch/2.0 + DS_coil_ThickX/2.0;
  Y = 0.0;
  //Z = z_Magnets_array[3];//z_formed_bellows + 76.09*inch + DS_coil_ThickY + DS_coil_depth/2.0;
  //Z = z_formed_bellows + 76.09*inch + DS_coil_ThickY + DS_coil_depth/2.0;
  
  Z = z_formed_bellows+Bellows2L;

  G4double DSCoil_offset_z =(DownstreamYokeDepth-DS_coil_height)/2+1.757*inch;

  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z-DSCoil_offset_z), DS_coil_log, "DS_coil_phys_left", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(-X,Y,Z-DSCoil_offset_z), DS_coil_log, "DS_coil_phys_right", worldlog, false, 1 , ChkOverlaps );

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
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z-DSCoil_offset_z), DSpole_log, "DSpole_phys_left", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(-X,Y,Z-DSCoil_offset_z), DSpole_log, "DSpole_phys_right", worldlog, false, 1 , ChkOverlaps );


 

}




















/*
//Lead shielding for GEn/SIDIS
void G4SBSBeamlineBuilder::MakeSIDISLead(G4LogicalVolume *worldlog){
  bool leadring = true;
  bool checkoverlaps = false;

  G4VisAttributes* LeadColor = new G4VisAttributes(G4Colour(0.4,0.4,0.4));
  G4VisAttributes* AlColor = new G4VisAttributes(G4Colour(0.75,0.75,0.75));
  
  G4double inch = 2.54*cm;
  G4double TargetCenter_zoffset = 6.50*inch;
  G4double z_conic_vacline_weldment = 62.8*inch - TargetCenter_zoffset;
  // G4double rout_OM0 = 3.15*inch;
  
  // Shielding for Moller electrons: at the source :)
  // Ring of material (Pb 2 cm ? or Al 12 cm ?), to stop electrons up to 100 MeV, between 6 and 12 degrees
  G4double th_ringshield = 0.5*inch;
  G4double z1_ringshield = 1.625*inch;
  G4double z2_ringshield = z1_ringshield+th_ringshield;
  G4double rin_ringshield = z2_ringshield*tan(6.0*deg);
  G4double rout_ringshield = z2_ringshield*tan(12.*deg);
  
  //G4cout << "rin_ringshield (mm) " << rin_ringshield << endl;
  
  G4Tubs* ringshield = new G4Tubs("ringshield", rin_ringshield, rout_ringshield, th_ringshield, 
				  0.0, 360.0*deg);
  
  G4LogicalVolume *ringshield_log = new G4LogicalVolume( ringshield, GetMaterial("Lead"), "ringshield_log" );

  //G4RotationMatrix* rot_temp = new G4RotationMatrix;
  //  rot_temp->rotateZ(+135.0*deg);
  
  if(leadring)new G4PVPlacement( rot_temp, G4ThreeVector( 0, 0, z1_ringshield+th_ringshield/2.0 ), ringshield_log, "ringshield_phys", worldlog, false, 0, checkoverlaps );
  ringshield_log->SetVisAttributes(LeadColor);
  
*/

/*
  // Shielding for Scattering chamber:
  // 
  G4double mindist_SCshield = 6.125*inch;
  G4double th_SCshield = 4.0*inch;
  G4double h_SCshield = 12.0*inch;//18.0*inch;
  G4double w1_SCshield = z1_ringshield*tan(25.1*deg)-mindist_SCshield;
  G4double w2_SCshield = (z1_ringshield+th_SCshield)*tan(25.1*deg)-mindist_SCshield;
  // G4double rin_ringshield = z1_ringshield*sin(6.0*deg);
  // G4double rout_ringshield = z2_ringshield*sin(12.0*deg);
  */
/*
  //Side Shields, exit beamline
  G4double z_SCshield = z1_ringshield+th_SCshield/2.0;
  G4double x_SCshield = mindist_SCshield+(w1_SCshield+w2_SCshield)/4.0;
  
  G4Trap* SCshield = new G4Trap("SCshield", h_SCshield, th_SCshield, w2_SCshield, 
   				w1_SCshield);
  
  G4LogicalVolume *SCshield_log = new G4LogicalVolume( SCshield, GetMaterial("Lead"), "SCshield_log" );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateX(+90*deg);
  
  //if(lead)
  new G4PVPlacement( rot_temp, G4ThreeVector( x_SCshield, 0, z_SCshield ), SCshield_log, "SCshield_phys", worldlog, false, 0, checkoverlaps );
  SCshield_log->SetVisAttributes(LeadColor);
  
  // Shielding for Spool piece.
  // 
  G4double z1_spoolshield = z1_ringshield+th_SCshield;
  G4double z2_spoolshield = //z1_spoolshield+1.89*m;
    z_conic_vacline_weldment + (0.84 + 0.14 + 11.62)*inch;
  
  G4double z_spoolshield = (z2_spoolshield+z1_spoolshield)/2.0;
  G4double L_spoolshield = z2_spoolshield-z1_spoolshield;
  G4double th_spoolshield = 2.0*inch;
  G4double d_spoolshield = mindist_SCshield + th_spoolshield/2.0;
  G4double H_spoolshield = h_SCshield;//12.0*inch;//4*rout_OM0;
  
  G4Box *spoolshield = new G4Box("spoolshield", th_spoolshield/2.0, H_spoolshield/2.0, L_spoolshield/2.0 );

  G4LogicalVolume *spoolshield_log = new G4LogicalVolume( spoolshield, GetMaterial("Lead"), "spoolshield_log" );

  rot_temp = new G4RotationMatrix;
  
  new G4PVPlacement( rot_temp, G4ThreeVector( d_spoolshield, 0, z_spoolshield ), spoolshield_log, "spoolshield_phys", worldlog, false, 0, checkoverlaps );
  spoolshield_log->SetVisAttributes(LeadColor);
  
  // Beamline shielding : between before 1st corrector magnets (BL4 only)
  //
  G4double th_BLshield1 = 2.0*inch;
  G4double L_BLshield1 = 34.0*inch;
  G4double h_BLshield1 = h_SCshield;
  
  G4double z_BLshield1 = z_conic_vacline_weldment + (0.84 + 0.14 + 45.62 + 14.38*0.65 - 34.0/2.0)*inch;
  G4double x_BLshield1 = (12.0)*inch;
  
  x_BLshield1+= th_BLshield1/2.0;
  
  G4Box *BLshield1 = new G4Box("BLshield1", th_BLshield1/2.0, h_BLshield1/2.0, L_BLshield1/2.0 );
  
  G4LogicalVolume *BLshield1_log = new G4LogicalVolume( BLshield1, GetMaterial("Lead"), "BLshield1_log" );
  
  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY(-1.5*deg);
  
  if(fDetCon->fBeamlineConf==4){
    //if(lead)
    new G4PVPlacement( rot_temp, G4ThreeVector( x_BLshield1, 0, z_BLshield1 ), BLshield1_log, "BLshield1_phys", worldlog, false, 0, checkoverlaps );
    BLshield1_log->SetVisAttributes(LeadColor);
  }

*/

// This is the "default" beam line (for C16)
void G4SBSBeamlineBuilder::MakeDefaultBeamline(G4LogicalVolume *worldlog){// Old beam line...
  G4SBS::Targ_t targtype = fDetCon->fTargType;

  double swallrad = 1.143*m/2;
  double swallrad_inner = 1.041/2.0*m; 
  //EFuchey: 2017/02/14: change parameters for Standard scat chamber:
  double sc_entbeampipeflange_dist = 25.375*2.54*cm;// entrance pipe flange distance from hall center
  double sc_exbeampipeflange_dist = 27.903*2.54*cm;// exit pipe flange distance from hall center
  
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
  
  //Don't add window: we want the beam to interact with the target first. Butt up against the outer edge of the scattering chamber:
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner), entLog_cut, "ent_phys", worldlog, false,0);
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner), entvacLog_cut, "entvac_phys", worldlog,false,0);

  // EFuchey: 2017/02/14: add the possibility to change the first parameters for the beam line polycone 
  // Default set of values;
  // double z0 = sc_exbeampipeflange_dist, rin_0 = 6.20*cm, rout_0 = (6.20+0.28*2.54)*cm;
  
  // if( fDetCon->fTargType == kH2 || fDetCon->fTargType == k3He || fDetCon->fTargType == kNeutTarg ){
  //   z0 = 37.2*cm;
  //   rin_0 = 5.64*cm;
  //   rout_0 = (5.64+0.28*2.54)*cm;
  // }  
  
  int nsec = 7;
  
  G4double exit_z[]   = { 162.2*cm, 592.2*cm, 609.84*cm,609.85*cm, 1161.02*cm, 1161.03*cm,2725.66*cm };
  G4double exit_z_vac[] = { 162.2*cm, 592.2*cm, 610.24*cm,610.35*cm, 1161.52*cm, 1161.53*cm,2726.46*cm };
  
  G4double exit_zero[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  //G4double exit_rin[] = { 4.8*cm, 14.8*cm,15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };
  //G4double exit_rou[] = { 5.0*cm, 15.0*cm,15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };
  
  G4double exit_rin[] = { 6.065*2.54*cm/2., 14.8*cm,15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };
  G4double exit_rou[] = { (6.065/2.0+0.28)*2.54*cm, 15.0*cm,15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };
  
  G4Polycone *ext_cone = new G4Polycone("ext_cone", 0.0*deg, 360.0*deg, nsec, exit_z, exit_rin, exit_rou);
  //G4Polycone *ext_vac  = new G4Polycone("ext_vac ", 0.0*deg, 360.0*deg, nsec, exit_z_vac, exit_zero, exit_rin);
  G4Polycone *ext_vac  = new G4Polycone("ext_vac ", 0.0*deg, 360.0*deg, nsec, exit_z, exit_zero, exit_rin);
  
  G4LogicalVolume *extLog = new G4LogicalVolume(ext_cone, GetMaterial("Aluminum"), "ext_log", 0, 0, 0);
  G4LogicalVolume *extvacLog = new G4LogicalVolume(ext_vac, GetMaterial("Vacuum"), "extvac_log", 0, 0, 0);
  
  new G4PVPlacement(0,G4ThreeVector(), extLog, "ext_phys", worldlog, false,0);
  new G4PVPlacement(0,G4ThreeVector(), extvacLog, "extvac_phys", worldlog,false,0);
    
  G4VisAttributes *extVisAtt= new G4VisAttributes(G4Colour(0.9,0.9,0.9));
  extLog->SetVisAttributes(extVisAtt);
  
  extvacLog->SetVisAttributes(G4VisAttributes::Invisible);
  entvacLog->SetVisAttributes(G4VisAttributes::Invisible);
    
  entvacLog_cut->SetVisAttributes(G4VisAttributes::Invisible);
    
  G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
    
  extLog->SetVisAttributes(pipeVisAtt);
  entLog->SetVisAttributes(pipeVisAtt);
  entLog_cut->SetVisAttributes(pipeVisAtt);
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

//lead shielding for GMn
void G4SBSBeamlineBuilder::MakeGMnLead(G4LogicalVolume *worldlog){
  bool leadring = true;
  bool checkoverlaps = false;

  G4VisAttributes* LeadColor = new G4VisAttributes(G4Colour(0.4,0.4,0.4));
  G4VisAttributes* AlColor = new G4VisAttributes(G4Colour(0.75,0.75,0.75));
  
  G4double inch = 2.54*cm;
  G4double TargetCenter_zoffset = 6.50*inch;
  G4double z_conic_vacline_weldment = 62.8*inch - TargetCenter_zoffset;
  G4double rout_OM0 = 3.15*inch;
  
  // Shielding for Moller electrons: at the source :)
  // Ring of material (Pb 2 cm ? or Al 12 cm ?), to stop electrons up to 100 MeV, between 6 and 12 degrees
  G4double th_ringshield = 1.5*inch;
  G4double z1_ringshield = 27.12*inch;
  G4double z2_ringshield = z1_ringshield+th_ringshield;
  G4double rin_ringshield = z2_ringshield*tan(6.0*deg);
  G4double rout_ringshield = z2_ringshield*tan(12.*deg);
  
  //G4cout << "rin_ringshield (mm) " << rin_ringshield << endl;
  
  G4Tubs* ringshield = new G4Tubs("ringshield", rin_ringshield, rout_ringshield, th_ringshield/2.0, 
				  0.0, 270.0*deg);
  
  G4LogicalVolume *ringshield_log = new G4LogicalVolume( ringshield, GetMaterial("Lead"), "ringshield_log" );

  G4RotationMatrix* rot_temp = new G4RotationMatrix;
  rot_temp->rotateZ(+135.0*deg);
  
  if(leadring)new G4PVPlacement( rot_temp, G4ThreeVector( 0, 0, z1_ringshield+th_ringshield/2.0 ), ringshield_log, "ringshield_phys", worldlog, false, 0, checkoverlaps );
  ringshield_log->SetVisAttributes(LeadColor);
  
  // Shielding for Scattering chamber:
  // 
  G4double mindist_SCshield = 6.125*inch;
  G4double th_SCshield = 4.0*inch;
  G4double h_SCshield = 12.0*inch;//18.0*inch;
  G4double w1_SCshield = z1_ringshield*tan(25.1*deg)-mindist_SCshield;
  G4double w2_SCshield = (z1_ringshield+th_SCshield)*tan(25.1*deg)-mindist_SCshield;
  // G4double rin_ringshield = z1_ringshield*sin(6.0*deg);
  // G4double rout_ringshield = z2_ringshield*sin(12.0*deg);
  
  G4double z_SCshield = z1_ringshield+th_SCshield/2.0;
  G4double x_SCshield = mindist_SCshield+(w1_SCshield+w2_SCshield)/4.0;
  
  //G4Trap* SCshield = new G4Trap("SCshield", w1_SCshield/2.0, w2_SCshield/2.0, 
  //h_SCshield/2.0, h_SCshield/2.0, th_SCshield/2.0);
  G4Trap* SCshield = new G4Trap("SCshield", h_SCshield, th_SCshield, w2_SCshield, 
   				w1_SCshield);
  
  //G4Box* SCshield = new G4Box("SCshield", w1_SCshield/2.0, h_SCshield/2.0, th_SCshield/2.0);//
  
  G4LogicalVolume *SCshield_log = new G4LogicalVolume( SCshield, GetMaterial("Lead"), "SCshield_log" );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateX(+90*deg);
  
  //if(lead)
  new G4PVPlacement( rot_temp, G4ThreeVector( x_SCshield, 0, z_SCshield ), SCshield_log, "SCshield_phys", worldlog, false, 0, checkoverlaps );
  SCshield_log->SetVisAttributes(LeadColor);
  
  // Shielding for Spool piece.
  // 
  G4double z1_spoolshield = z1_ringshield+th_SCshield;
  G4double z2_spoolshield = //z1_spoolshield+1.89*m;
    z_conic_vacline_weldment + (0.84 + 0.14 + 11.62)*inch;
  
  G4double z_spoolshield = (z2_spoolshield+z1_spoolshield)/2.0;
  G4double L_spoolshield = z2_spoolshield-z1_spoolshield;
  G4double th_spoolshield = 2.0*inch;
  G4double d_spoolshield = mindist_SCshield + th_spoolshield/2.0;
  G4double H_spoolshield = h_SCshield;//12.0*inch;//4*rout_OM0;
  
  G4Box *spoolshield = new G4Box("spoolshield", th_spoolshield/2.0, H_spoolshield/2.0, L_spoolshield/2.0 );

  G4LogicalVolume *spoolshield_log = new G4LogicalVolume( spoolshield, GetMaterial("Lead"), "spoolshield_log" );

  rot_temp = new G4RotationMatrix;
  
  new G4PVPlacement( rot_temp, G4ThreeVector( d_spoolshield, 0, z_spoolshield ), spoolshield_log, "spoolshield_phys", worldlog, false, 0, checkoverlaps );
  spoolshield_log->SetVisAttributes(LeadColor);
  
  // Beamline shielding : between before 1st corrector magnets (BL4 only)
  //
  G4double th_BLshield1 = 2.0*inch;
  G4double L_BLshield1 = 34.0*inch;
  G4double h_BLshield1 = h_SCshield;
  
  G4double z_BLshield1 = z_conic_vacline_weldment + (0.84 + 0.14 + 45.62 + 14.38*0.65 - 34.0/2.0)*inch;
  G4double x_BLshield1 = (12.0)*inch;
  
  x_BLshield1+= th_BLshield1/2.0;
  
  G4Box *BLshield1 = new G4Box("BLshield1", th_BLshield1/2.0, h_BLshield1/2.0, L_BLshield1/2.0 );
  
  G4LogicalVolume *BLshield1_log = new G4LogicalVolume( BLshield1, GetMaterial("Lead"), "BLshield1_log" );
  
  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY(-1.5*deg);
  
  if(fDetCon->fBeamlineConf==4){
    //if(lead)
    new G4PVPlacement( rot_temp, G4ThreeVector( x_BLshield1, 0, z_BLshield1 ), BLshield1_log, "BLshield1_phys", worldlog, false, 0, checkoverlaps );
    BLshield1_log->SetVisAttributes(LeadColor);
  }
  
 



  /*
  // Beamline shielding : between corrector magnets
  //
  G4double th_BLshield2 = 2.0*inch;
  G4double L_BLshield2 = 55.0*inch;
  G4double h_BLshield2 = 12.0*inch;
  
  G4double z_BLshield2 = z_conic_vacline_weldment + (0.84 + 0.14 + 11.62 + 14.38 + 53.62/2.0)*inch;
  G4double x_BLshield2 = (3.6725+1.5)*inch;
  
  if(fDetCon->fBeamlineConf==4){
    z_BLshield2 = z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38 + 51.62/2.0)*inch;
    x_BLshield2 = (4.55075+1.5)*inch;
  }
  
  x_BLshield2+= th_BLshield2/2.0;
  
  G4Box *BLshield2 = new G4Box("BLshield2", th_BLshield2/2.0, h_BLshield2/2.0, L_BLshield2/2.0 );

  G4LogicalVolume *BLshield2_log = new G4LogicalVolume( BLshield2, GetMaterial("Lead"), "BLshield2_log" );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY(-1.5*deg);
  
  if(lead)new G4PVPlacement( rot_temp, G4ThreeVector( x_BLshield2, 0, z_BLshield2 ), BLshield2_log, "BLshield2_phys", worldlog, false, 0 );
  BLshield2_log->SetVisAttributes(LeadColor);
  
  // Beamline shielding : after 2nd corrector magnets (BL3 only)
  //
  G4double th_BLshield3 = 2.0*inch;
  G4double L_BLshield3 = 32.0*inch;
  G4double h_BLshield3 = 12.0*inch;
  
  G4double z_BLshield3 = z_conic_vacline_weldment + (0.84 + 0.14 + 120.5)*inch;
  G4double x_BLshield3 = (5.43825+1.5)*inch;
  
  x_BLshield3+= th_BLshield3/2.0;
  
  G4Box *BLshield3 = new G4Box("BLshield3", th_BLshield3/2.0, h_BLshield3/2.0, L_BLshield3/2.0 );
  
  G4LogicalVolume *BLshield3_log = new G4LogicalVolume( BLshield3, GetMaterial("Lead"), "BLshield3_log" );
  
  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY(-1.5*deg);
  
  if(fDetCon->fBeamlineConf==3){
    if(lead)new G4PVPlacement( rot_temp, G4ThreeVector( x_BLshield3, 0, z_BLshield3 ), BLshield3_log, "BLshield3_phys", worldlog, false, 0 );
    BLshield3_log->SetVisAttributes(LeadColor);
  }
  */
  
  
  
  /*
  //"Wood" shielding ? -> Polyethylene...
  
  G4double L_woodshield = 3.0*m;
  G4double h_woodshield = 4.0*m;
  G4double th1_woodshield = 0.525*m;//z1_ringshield*tan(25.1*deg)-rout_ringshield;
  G4double th2_woodshield = 0.9*m;//z2_ringshield*tan(25.1*deg)-rout_ringshield;
  // G4double rin_ringshield = z1_ringshield*sin(6.0*deg);
  // G4double rout_ringshield = z2_ringshield*sin(12.0*deg);
  
  G4double z_woodshield = 4.3*m*cos(16.7*deg);
  G4double x_woodshield = 4.3*m*sin(16.7*deg);

  G4Trap* sideshield = new G4Trap("sideshield", h_woodshield, L_woodshield, th2_woodshield, 
				  th1_woodshield);
  
  G4LogicalVolume *sideshield_log = new G4LogicalVolume( sideshield, GetMaterial("Air"), "sideshield_log" );  

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY(-25.75*deg);
  rot_temp->rotateX(+90*deg);
  
  new G4PVPlacement( rot_temp, G4ThreeVector( x_woodshield, 0, z_woodshield ), sideshield_log, "sideshield_phys", worldlog, false, 0 );
  sideshield_log->SetVisAttributes( G4VisAttributes::Invisible );
  
  G4Trap* woodshield = new G4Trap("woodshield", h_woodshield, L_woodshield, th2_woodshield-2.5*cm, 
				  th1_woodshield-2.5*cm);

  G4LogicalVolume *woodshield_log = new G4LogicalVolume( woodshield, GetMaterial("Polyethylene"), "woodshield_log" );  

  rot_temp = new G4RotationMatrix;
  
  //new G4PVPlacement( rot_temp, G4ThreeVector( -2.5*cm, 0, 0), woodshield_log, "woodshield_phys", sideshield_log, false, 0 );
  woodshield_log->SetVisAttributes( G4Colour(0.2,0.9,0.75));

  G4Box* leadblanket = new G4Box("leadblanket", 1.0*cm/2.0, L_woodshield/2.0-1.0*mm, h_woodshield/2.0-1.0*mm); 

  G4LogicalVolume *leadblanket_log = new G4LogicalVolume( leadblanket, GetMaterial("Lead"), "leadblanket_log" ); 
  leadblanket_log->SetVisAttributes(LeadColor);
  
  rot_temp = new G4RotationMatrix;
  rot_temp->rotateZ(-7.05*deg);
  
  //new G4PVPlacement( rot_temp, G4ThreeVector( +0.335*m, 0, 0 ), leadblanket_log, "leadblanket_phys", sideshield_log, false, 0 );
  */
  
  G4double L_sideshield = 2.0*m;
  G4double h_sideshield = 1.3*m;
  G4double th_sideshield = 0.50*m;
  // To be reoptimized

  G4double z_sideshield = z_conic_vacline_weldment + (0.84 + 0.14 + 45.63 + 14.38 + (51.62+22.38)/2.0)*inch;
  G4double x_sideshield = 40.0*cm+th_sideshield/2.0;
  
  G4Box* sideshield = new G4Box("sideshield", th_sideshield/2.0, h_sideshield/2.0, L_sideshield/2.0);
  
  G4LogicalVolume *sideshield_log = new G4LogicalVolume( sideshield, GetMaterial("Air"), "sideshield_log" );  

  rot_temp = new G4RotationMatrix;
  
  new G4PVPlacement( rot_temp, G4ThreeVector( x_sideshield, 0, z_sideshield ), sideshield_log, "sideshield_phys", worldlog, false, 0, checkoverlaps );
  sideshield_log->SetVisAttributes( G4VisAttributes::Invisible );
  
  G4double th_Alshield = 4.0*inch;
  G4double th_SSshield = 1.0*inch;
    
  G4Box* Alshield = new G4Box("Alshield", th_Alshield/2.0, h_sideshield/2.0, (L_sideshield-50.0*cm)/2.0);
  
  G4LogicalVolume *Alshield_log = new G4LogicalVolume( Alshield, GetMaterial("Aluminum"), "Alshield_log" );  

  rot_temp = new G4RotationMatrix;
  
  if(!leadring)new G4PVPlacement( rot_temp, G4ThreeVector( -th_sideshield/2.0+th_Alshield/2.0, 0, -25.0*cm), Alshield_log, "Alshield_phys", sideshield_log, false, 0, checkoverlaps );
  Alshield_log->SetVisAttributes(AlColor);

  G4Box* leadblanket = new G4Box("leadblanket", th_SSshield/2.0, h_sideshield/2.0, (L_sideshield-50.0*cm)/2.0); 

  G4LogicalVolume *leadblanket_log = new G4LogicalVolume( leadblanket, GetMaterial("Lead"), "leadblanket_log" ); 
  leadblanket_log->SetVisAttributes(LeadColor);
  
  rot_temp = new G4RotationMatrix;
  
  if(!leadring)new G4PVPlacement( rot_temp, G4ThreeVector( -th_sideshield/2.0+th_Alshield+th_SSshield/2.0, 0, -25.0*cm ), leadblanket_log, "leadblanket_phys", sideshield_log, false, 0, checkoverlaps );
  /**/
}


void G4SBSBeamlineBuilder::MakeGEnRPLead(G4LogicalVolume *worldlog){

  bool chkovrlps = false;
  G4VisAttributes* LeadColor = new G4VisAttributes(G4Colour(0.4,0.4,0.4));

  G4double l_leadwall1 = 5.0*12.0*2.54*cm;
  G4double h_leadwall1 = 42.0*2.54*cm;
  G4double w_leadwall1 = 2*2.54*cm;

  G4Box *leadwall1 = new G4Box("leadwall1", w_leadwall1/2.0, h_leadwall1/2.0, l_leadwall1/2.0);
  G4LogicalVolume *leadwall1_log = new G4LogicalVolume(leadwall1, GetMaterial("Lead"), "leadwall1_log");

  G4double X1 = -17.0*2.54*cm; 
  G4double Y1 = 0.0;
  G4double Z1 = 2.5*76.325*2.54*cm + 1.5*2.54*cm;  // guesstimate 
  
  new G4PVPlacement( 0, G4ThreeVector(X1,Y1,Z1), leadwall1_log, "leadwall1_phys", worldlog, false, 0, chkovrlps );
  leadwall1_log->SetVisAttributes(LeadColor);
  
  G4double l_leadwall2 = 3.0*12.0*2.54*cm;
  G4double h_leadwall2 = h_leadwall1;
  G4double w_leadwall2 = w_leadwall1;

  G4Box *leadwall2 = new G4Box("leadwall2", w_leadwall2/2.0, h_leadwall2/2.0, l_leadwall2/2.0);
  G4LogicalVolume *leadwall2_log = new G4LogicalVolume(leadwall2, GetMaterial("Lead"), "leadwall2_log");

  G4double X2 = -21.0*2.54*cm; 
  G4double Y2 = 0.0;
  G4double Z2 = Z1 + l_leadwall1/2.0 + l_leadwall2/2.0 - 1.*2.54*cm;  

  new G4PVPlacement( 0, G4ThreeVector(X2,Y2,Z2), leadwall2_log, "leadwall2_phys", worldlog, false, 0, chkovrlps );
  leadwall2_log->SetVisAttributes(LeadColor);
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

void G4SBSBeamlineBuilder::MakeToyBeamline(G4LogicalVolume *motherlog){ //This is the "toy" beamline to go with the "toy" scattering chamber:

  // Don't do anything yet, just make the code compile;
}

void G4SBSBeamlineBuilder::MakeBeamDiffuser(G4LogicalVolume *logicMother){
   // A beam diffuser that sits right in front of the beam dump
   // Added by D. Flay (JLab) in Aug 2020  
      
   G4cout << "[G4SBSBeamlineBuilder]: Adding the Beam Diffuser to the beam line..." << G4endl; 
   
   G4double inch        = 25.4*mm;

   // A case for diffuser
   // - made of vacuum 
   // - allows placement of the volume in same mother as the calorimeter
   //   (can't have two replicas or parameterised volumes in same mother...)  
   // G4double diffCase_x  = 12.*inch;
   // G4double diffCase_y  = 6.*inch;   
   // G4double diffCase_z  = 15.*cm;    
   // G4VSolid *diffCaseS  = new G4Box("diffCase",diffCase_x/2.,diffCase_y/2.,diffCase_z/2.); 

   // make the case match the FakeVacuumExtension dimensions 
   G4double dcRmin     = 0.0*cm;
   G4double dcRmax     = 12*inch;
   G4double dcLen      = 15.*cm;
   G4double dcStartPhi = 0*deg;  
   G4double dcDPhi     = 360*deg;  
   G4VSolid *diffCaseS = new G4Tubs("diffCase",dcRmin,dcRmax,dcLen/2.,dcStartPhi,dcDPhi);

   G4LogicalVolume *diffCaseLV = new G4LogicalVolume(diffCaseS,GetMaterial("Vacuum"),"diffCase");

   // where to place the diffuser 
   // note: the (x,y) center of the diffuser plates is centered on this logical volume 
   // double ft   = 12.*inch; 
   // double beamHeight = 10.0*ft; // 10 feet off the ground
   // double floorThick = 1.0*m;   // floor is 1 m thick 
   G4double xd = 0.;
   G4double yd = 0; // beamHeight + 0.5*floorThick; // do we need this?
   G4double zd = 11.75*m;  // FakeVacuumExtension ends at xx m  
   G4ThreeVector P_case = G4ThreeVector(xd,yd,zd);

   bool checkOverlaps = true;

   new G4PVPlacement(0,                // no rotation
	             P_case,           // location in mother volume 
	             diffCaseLV,       // its logical volume                         
	             "diffCase",       // its name
	             logicMother,      // its mother  volume
	             false,            // no boolean operation
	             0,                // copy number
	             checkOverlaps);   // checking overlaps 

   G4VisAttributes *visCase = new G4VisAttributes();
   visCase->SetForceWireframe();

   // diffCaseLV->SetVisAttributes(G4VisAttributes::GetInvisible());
   diffCaseLV->SetVisAttributes(visCase);

   // parameterised build of the diffuser
   // build first plate (same for Hall A or C)  
   G4double r_min    = 2.*inch;
   G4double r_max    = 5.*inch;
   G4double thk      = 0.125*inch;
   G4double startPhi = 255.*deg;
   G4double dPhi     = 30.*deg;

   // choose the origin of the device (where the first plate starts, relative to the mother volume)  
   G4ThreeVector P0 = G4ThreeVector(0,0,0);

   G4VisAttributes *vis = new G4VisAttributes();
   vis->SetColour( G4Colour::Blue() );

   // first plate 
   G4VSolid *plateSolid     = new G4Tubs("plate",r_min,r_max,thk/2.,startPhi,dPhi);
   G4LogicalVolume *plateLV = new G4LogicalVolume(plateSolid,GetMaterial("Aluminum"),"plateLV");
   plateLV->SetVisAttributes(vis);

   // parameterisation (Hall A)
   char Hall   = 'A';
   int NPlanes = 15; 
   if(Hall=='C') NPlanes = 16;  
   G4VPVParameterisation *plateParam = new G4SBSBDParameterisation(Hall,P0);
   // placement
   new G4PVParameterised("BeamDiffuser",plateLV,diffCaseLV,kZAxis,NPlanes,plateParam);

   // Attach sensitive detector (SD) functionality; follow the GEM example

   // name of SD and the hitCollection  
   G4String bdSDname = "BD";  // FIXME: is this ok, or do we need directory structure like the GEMs? 
   // We have to remove all the directory structure from the 
   // Hits Collection name or else GEANT4 SDmanager routines will not handle correctly.
   G4String bdSDname_nopath = bdSDname;
   bdSDname_nopath.remove(0,bdSDname.last('/')+1);
   G4String bdColName = bdSDname_nopath; 
   bdColName += "HitsCollection";

   // check to see if this SD exists already; if not, create a new SD object and append to the list of SDs  
   G4SBSBeamDiffuserSD *bdSD; 
   if( !(bdSD = (G4SBSBeamDiffuserSD *)fDetCon->fSDman->FindSensitiveDetector(bdSDname)) ){
      G4cout << "[G4SBSBeamlineBuilder]: Adding Beam Diffuser SD functionality..." << G4endl; 
      bdSD = new G4SBSBeamDiffuserSD(bdSDname,bdColName);
      plateLV->SetSensitiveDetector(bdSD);  
      fDetCon->fSDman->AddNewDetector(bdSD);
      (fDetCon->SDlist).insert(bdSDname); 
      fDetCon->SDtype[bdSDname] = G4SBS::kBD; 
   }

   G4cout << "[G4SBSBeamlineBuilder]: --> Done." << G4endl;
 
}
