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
#include "G4Sphere.hh"
#include "G4Polycone.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4RotationMatrix.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4SBSCalSD.hh"

#include "G4SBSBDParameterisation.hh"
#include "G4SBSBeamDiffuserSD.hh"
#include "TString.h"

G4SBSBeamlineBuilder::G4SBSBeamlineBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
  assert(dc);
}

G4SBSBeamlineBuilder::~G4SBSBeamlineBuilder(){;}

void G4SBSBeamlineBuilder::BuildComponent(G4LogicalVolume *worldlog){
  
  double beamheight = 10.0*12*2.54*cm; // 10 feet off the ground
 
  // EFuchey 2017/03/29: organized better this with a switch instead of an endless chain of if...  else...
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
  case(G4SBS::kGEp_BB):
    fDetCon->fBeamlineConf = 3;
    MakeGMnBeamline(worldlog);
    //if(fDetCon->fLeadOption == 1){
    MakeGMnLead(worldlog);
    //}
    break;
  case(G4SBS::kGEnRP):// GEnRP
    fDetCon->fBeamlineConf = 3;
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
    if(fDetCon->fLeadOption == 1){
      //Sseeds - leaving option in, but removing obsolete geometry
      MakeGEnLead(worldlog);
    }
    break;
  case(G4SBS::kSIDISExp):// SIDIS
    fDetCon->fBeamlineConf = 2;
    Make3HeBeamline(worldlog);
    if(fDetCon->fLeadOption == 1){
      //Sseeds - leaving option in, but removing obsolete geometry
      MakeSIDISLead(worldlog);
    }
    break;
  case(G4SBS::kALL):// ALL
    fDetCon->fBeamlineConf = 2;
    Make3HeBeamline(worldlog);
    if(fDetCon->fLeadOption == 1){
      //Sseeds - leaving option in, but removing obsolete geometry
      MakeALLLead(worldlog);
    }
    break;
  case(G4SBS::kGEMHCtest):// Hall C GEM test
    fDetCon->fBeamlineConf = 3;
    MakeGMnBeamline(worldlog);
    break;  
  case(G4SBS::kTDIS):
    fDetCon->fBeamlineConf = 2;
    MakeTDISBeamline(worldlog);
    break;  
  case(G4SBS::kNDVCS):
    fDetCon->fBeamlineConf = 2;
    MakeTDISBeamline(worldlog);
    break;  
  case(G4SBS::kMTPConly):
    //fDetCon->fBeamlineConf = 2;
    //MakeTDISBeamline(worldlog);
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
  
  floorLog->SetVisAttributes(G4VisAttributes::GetInvisible());

  return;

}

void G4SBSBeamlineBuilder::MakeCommonExitBeamline(G4LogicalVolume *worldlog) {
  bool ChkOverlaps = false;
  //Define visualization attributes here:
  G4VisAttributes *ironColor= new G4VisAttributes(G4Colour(0.3,0.3,0.3));
  G4VisAttributes *AlColor= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
  G4VisAttributes *Vacuum_visatt = new G4VisAttributes(G4Colour(0.1, 0.5, 0.9 ) );
  Vacuum_visatt->SetVisibility(false);
  G4VisAttributes *SteelColor = new G4VisAttributes( G4Colour( 0.75, 0.75, 0.75 ) );
  G4VisAttributes *CopperColor = new G4VisAttributes( G4Colour( 0.7, 0.3, 0.3 ) );
  
  G4double inch = 2.54*cm;
  
  G4double TargetCenter_zoffset = 0.0*inch;  //SSeeds. Working to eliminate this offset here
  
  G4double z_formed_bellows = 52.440*inch - TargetCenter_zoffset; //relative to "target center"? or "origin"?
  G4double z_spool_piece = 100.44*inch - TargetCenter_zoffset;
  if(fDetCon->fBeamlineConf>2)z_spool_piece = 27.903*inch;
  G4double z_conic_vacline_weldment = 63.331*inch - TargetCenter_zoffset; 
  G4double z_outer_magnetic = 71.782*inch - TargetCenter_zoffset;
  G4double z_inner_magnetic = 73.782*inch - TargetCenter_zoffset;
  G4double z_welded_bellows = 201.132*inch; 

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
  //Z = z_conic_vacline_weldment + Thick/2.0;
  //new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange1_vac_log, "CVLW_Flange1_vac_phys", worldlog, false, 0 , ChkOverlaps );
  //new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange1_log, "CVLW_Flange1_phys", worldlog, false, 0 , ChkOverlaps );

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
  //Z = z_conic_vacline_weldment + 0.84*inch + Thick/2.0;
  Z = z_conic_vacline_weldment + Thick/2.0 + 0.31*inch;  //SSJT, added 0.31*inch to preserve ring dims from previous build

  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_vac_log, "CVLW_vac_phys", worldlog, false, 0 , ChkOverlaps);
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_log, "CVLW_phys", worldlog, false, 0 , ChkOverlaps );

  //Section 6 Bellows (S6) - Leading into target to midpipe section.

  //Cylinder 1 (C1)
  Rin = 5.875*inch;
  Rout = 7.0*inch;
  Thick = 1.64*inch;
  G4Tubs *S6_C1 = new G4Tubs("S6_C1",Rin,Rout,Thick/2.0,0.0,twopi);
  G4Tubs *S6_C1_vac = new G4Tubs("S6_C1_vac",0.0,Rin,Thick/2.0,0.0,twopi);

  G4LogicalVolume *S6_C1_log = new G4LogicalVolume(S6_C1,GetMaterial("Stainless_Steel"),"S2_C1_log" );
  G4LogicalVolume *S6_C1_vaclog = new G4LogicalVolume(S6_C1_vac,GetMaterial("Vacuum"),"S2_C1_vaclog" );

  S6_C1_log->SetVisAttributes(SteelColor);
  
  Z = 200.512*inch+Thick/2.0; //Direct JT measurement to start

  new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S6_C1_log,"S6_C1_phys",worldlog,false,0,ChkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S6_C1_vaclog,"S6_C1_vac",worldlog,false,0,ChkOverlaps);

  //Cylinder 2 (C2)
  Z += Thick/2.0;
  Rin = 5.88*inch;
  Rout = 6.0*inch;
  Thick = 3.267*inch;
  G4Tubs *S6_C2 = new G4Tubs("S6_C2",Rin,Rout,Thick/2.0,0.0,twopi);
  G4Tubs *S6_C2_vac = new G4Tubs("S6_C2_vac",0.0,Rin,Thick/2.0,0.0,twopi);

  G4LogicalVolume *S6_C2_log = new G4LogicalVolume(S6_C2,GetMaterial("Stainless_Steel"),"S6_C2_log" );
  G4LogicalVolume *S6_C2_vaclog = new G4LogicalVolume(S6_C2_vac,GetMaterial("Vacuum"),"S6_C2_vaclog" );

  S6_C2_log->SetVisAttributes(SteelColor);
  
  Z += Thick/2.0;

  new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S6_C2_log,"S6_C1_phys",worldlog,false,0,ChkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S6_C2_vaclog,"S6_C1_vac",worldlog,false,0,ChkOverlaps);

  //Cylinder 3 (C3)
  Z += Thick/2.0;
  Rin = 5.88*inch;
  Rout = 8.25*inch;
  Thick = 1.766*inch;
  G4Tubs *S6_C3 = new G4Tubs("S6_C3",Rin,Rout,Thick/2.0,0.0,twopi);
  G4Tubs *S6_C3_vac = new G4Tubs("S6_C3_vac",0.0,Rin,Thick/2.0,0.0,twopi);

  G4LogicalVolume *S6_C3_log = new G4LogicalVolume(S6_C3,GetMaterial("Stainless_Steel"),"S6_C3_log" );
  G4LogicalVolume *S6_C3_vaclog = new G4LogicalVolume(S6_C3_vac,GetMaterial("Vacuum"),"S6_C3_vaclog" );

  S6_C3_log->SetVisAttributes(SteelColor);
  
  Z += Thick/2.0;

  new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S6_C3_log,"S6_C3_phys",worldlog,false,0,ChkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S6_C3_vaclog,"S6_C3_vac",worldlog,false,0,ChkOverlaps);


  //-----------------------------------------------------
  //       magnetic tubes
    
  // Inner and outer Magnetic shieldings: 
  // arrays of variables to parameterize the different geometries with the beamline config flag
  G4int Ndivs = 1;// number of segments with shielding
  G4double MPA = 1.5*deg;
  G4double MPTh = 0.25*inch;
  G4double MPinit_disp = 9.42*inch;
  G4double SSTh = 2.5*inch;
  std::vector<G4double> Rin_array;// radii for inner shielding elements
  std::vector<G4double> Zin_array;// z for inner shielding elements
  std::vector<G4int> Nrings_out;// number of outer elements per segments
  std::vector<G4double> Rout_array;// radii for inner shielding elements
  std::vector<G4double> Zout_array;// z for inner shielding elements
  std::vector<G4double> W1_array;// first width of side shields
  std::vector<G4double> W2_array;// last width of side shields
  std::vector<G4double> MPl_array;// Mounting plate lengths
  std::vector<G4double> MPxdisp_array;// x displacement for each shield and mounting plate section
  std::vector<G4double> SSxdisp_array;
  std::vector<G4double> MPzmin_array;

  G4double OMthick = 1.625*inch;
  G4double OMspace = 0.375*inch;

  //SSeeds will likely need to eliminate these options from the commonexitbeamline function and keep only option 3. After GEp and SIDIS layouts are known, will evaluate. Jan 2021
  /*
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
    
    //Rin_array.push_back( 3.7745*inch/2.0 );
    //Rin_array.push_back( 4.307*inch/2.0 );
    //Rin_array.push_back( 5.264*inch/2.0 );
    //Rin_array.push_back( 7.905*inch/2.0 );
    //Rin_array.push_back( 9.310*inch/2.0 );
    //Rin_array.push_back( 10.930*inch/2.0 );
    

    
    //Rin_array.push_back( 2.527*inch );
    //Rin_array.push_back( 2.58*inch );
    //Rin_array.push_back( 2.632*inch );
    //Rin_array.push_back( 2.685*inch );
    //Rin_array.push_back( 2.685*inch );
    //Rin_array.push_back( 2.737*inch );
    
    
    //SSeeds - updating with direct measurements to accomodate shift in z. Would be better to loop over calculation for ring radii. Will update.
    Rin_array.push_back( 1.895*inch );
    Rin_array.push_back( 2.157*inch );
    Rin_array.push_back( 2.628*inch );
    Rin_array.push_back( 3.928*inch );
    Rin_array.push_back( 4.619*inch );
    Rin_array.push_back( 5.446*inch );
    
    
    
    //Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14)*inch );
    //Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 +11.62 - 1.625)*inch );
    //Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 +11.62 + 14.38 + 2.0)*inch );
    //Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 - 2.0)*inch );
    //Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 + 22.38 + 2.0)*inch );
    //Zin_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 138.83 - 1.12 - 1.14)*inch );
    
    
    Zin_array.push_back( z_conic_vacline_weldment + 0.451*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62 - 1.625)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62 + 14.38 + 2.0)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62  + 14.38 + 53.62 - 2.0)*inch );
    Zin_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62  + 14.38 + 53.62 + 22.38 + 2.0)*inch );
    //Zin_array.push_back( z_conic_vacline_weldment + (0.451 + 138.83 - 1.12 - 1.14)*inch );
    Zin_array.push_back((167.782+31.591)*inch ); //SSeeds JT direct measure

    //Trapezoidal mounting plate shorter end length
    W1_array.push_back(5.646*inch);
    W1_array.push_back(6.575*inch);
    W1_array.push_back(9.435*inch);
    
    //Trapezoidal mounting plate longer end length
    W2_array.push_back(6.16*inch);
    W2_array.push_back(8.75*inch);
    W2_array.push_back(10.68*inch);

    //Overall length of mounting plates and shield supports
    MPl_array.push_back(13.875*inch);
    MPl_array.push_back(58.75*inch);
    MPl_array.push_back(33.631*inch);

    //Starting x disp of each mounting plate section
    //MPxdisp_array.push_back(6.521*inch/2.0+MPTh/2.0);
    MPxdisp_array.push_back(9.420*inch/2.0+MPTh/2.0);
    MPxdisp_array.push_back(10.332*inch/2.0+MPTh/2.0);
    MPxdisp_array.push_back(13.198*inch/2.0+MPTh/2.0);

    //Starting x disp of each shield support section
    SSxdisp_array.push_back(6.521*inch/2.0+MPTh/2.0);
    SSxdisp_array.push_back(7.812*inch/2.0+MPTh/2.0);
    SSxdisp_array.push_back(11.864*inch/2.0+MPTh/2.0);


    MPzmin_array.push_back( z_conic_vacline_weldment + 0.451*inch );
    MPzmin_array.push_back(88.375*inch);
    MPzmin_array.push_back(z_conic_vacline_weldment + (0.451 + 11.62  + 14.38 + 53.62 + 22.38)*inch );
    
    Nrings_out.push_back( 6 );
    Nrings_out.push_back( 27 );
    Nrings_out.push_back( 17 );
    
    
    //Rout_array.push_back( (3.7745/2+0.38)*inch );
    //Rout_array.push_back( (4.392/2+0.38)*inch );
    //Rout_array.push_back( (5.158/2+0.38)*inch );
    //Rout_array.push_back( (8.012/2+0.38)*inch );
    //Rout_array.push_back( (9.203/2+0.38)*inch );
    //Rout_array.push_back( (10.923/2+0.38)*inch );
    

    //SSeeds - updating with direct measurements to accomodate shift in z. Would be better to loop over calculation for ring radii. Will update.
    Rout_array.push_back( 2.527*inch );
    Rout_array.push_back( 2.832*inch );
    Rout_array.push_back( 3.208*inch );
    Rout_array.push_back( 4.613*inch );
    Rout_array.push_back( 5.198*inch );
    Rout_array.push_back( 6.079*inch );
    
    
    //Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14)*inch );
    //Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62)*inch );
    //Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38)*inch );
    //Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62)*inch );
    //Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 11.62  + 14.38 + 53.62 + 22.38)*inch );
    //Zout_array.push_back( z_conic_vacline_weldment + (0.84 + 138.83 - 1.12 - 1.14)*inch );
    
    
    Zout_array.push_back( z_conic_vacline_weldment + 0.451*inch);
    Zout_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62  + 14.38)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62  + 14.38 + 53.62)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62  + 14.38 + 53.62 + 22.38)*inch );
    Zout_array.push_back( z_conic_vacline_weldment + (0.451 + 138.83 - 1.12 - 1.0)*inch );
    
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
*/

  //SSeeds - updating with direct measurements to accomodate shift in z. Would be better to loop over calculation for ring radii. Will update.
  Ndivs = 3;  

  Rin_array.push_back( 1.895*inch );
  Rin_array.push_back( 2.157*inch );
  Rin_array.push_back( 2.628*inch );
  Rin_array.push_back( 3.928*inch );
  Rin_array.push_back( 4.619*inch );
  Rin_array.push_back( 5.446*inch );
    
  Zin_array.push_back( z_conic_vacline_weldment + 0.451*inch );
  Zin_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62 - 1.625)*inch );
  Zin_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62 + 14.38 + 2.0)*inch );
  Zin_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62  + 14.38 + 53.62 - 2.0)*inch );
  Zin_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62  + 14.38 + 53.62 + 22.38 + 2.0)*inch );
  //Zin_array.push_back( z_conic_vacline_weldment + (0.451 + 138.83 - 1.12 - 1.14)*inch );
  Zin_array.push_back((167.782+31.591)*inch ); //SSeeds JT direct measure

  //Trapezoidal mounting plate shorter end length
  W1_array.push_back(5.646*inch);
  W1_array.push_back(6.575*inch);
  W1_array.push_back(9.435*inch);
    
  //Trapezoidal mounting plate longer end length
  W2_array.push_back(6.16*inch);
  W2_array.push_back(8.75*inch);
  W2_array.push_back(10.68*inch);

  //Overall length of mounting plates and shield supports
  MPl_array.push_back(13.875*inch);
  MPl_array.push_back(58.75*inch);
  MPl_array.push_back(33.631*inch);

  //Starting x disp of each mounting plate section
  //MPxdisp_array.push_back(6.521*inch/2.0+MPTh/2.0);
  MPxdisp_array.push_back(9.420*inch/2.0+MPTh/2.0);
  MPxdisp_array.push_back(10.332*inch/2.0+MPTh/2.0);
  MPxdisp_array.push_back(13.198*inch/2.0+MPTh/2.0);

  //Starting x disp of each shield support section
  SSxdisp_array.push_back(6.521*inch/2.0+MPTh/2.0);
  SSxdisp_array.push_back(7.812*inch/2.0+MPTh/2.0);
  SSxdisp_array.push_back(11.864*inch/2.0+MPTh/2.0);


  MPzmin_array.push_back( z_conic_vacline_weldment + 0.451*inch );
  MPzmin_array.push_back(88.375*inch);
  MPzmin_array.push_back(z_conic_vacline_weldment + (0.451 + 11.62  + 14.38 + 53.62 + 22.38)*inch );
    
  Nrings_out.push_back( 6 );
  Nrings_out.push_back( 27 );
  Nrings_out.push_back( 17 );

  //SSeeds - updating with direct measurements to accomodate shift in z. Would be better to loop over calculation for ring radii. Will update.
  Rout_array.push_back( 2.527*inch );
  Rout_array.push_back( 2.832*inch );
  Rout_array.push_back( 3.208*inch );
  Rout_array.push_back( 4.613*inch );
  Rout_array.push_back( 5.198*inch );
  Rout_array.push_back( 6.079*inch );
    
  Zout_array.push_back( z_conic_vacline_weldment + 0.451*inch);
  Zout_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62)*inch );
  Zout_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62  + 14.38)*inch );
  Zout_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62  + 14.38 + 53.62)*inch );
  Zout_array.push_back( z_conic_vacline_weldment + (0.451 + 11.62  + 14.38 + 53.62 + 22.38)*inch );
  Zout_array.push_back( z_conic_vacline_weldment + (0.451 + 138.83 - 1.12 - 1.0)*inch );
  
  // Building beamline mounting plates:
  for(G4int i = 0; i<Ndivs; i++){
    
    //Build Shields and Mounts
    //Specifications
    G4double SMPL = MPl_array[i]/2.0;
    G4double SMPW1 = W1_array[i]/2.0;
    G4double SMPW2 = W2_array[i]/2.0;
    G4double SMP_xoffset = MPxdisp_array[i];
    G4double SS_xoffset = SSxdisp_array[i];
    
    //Symmetric mounting plates - construct using trapezoid class 
    G4Trd *SMP1 = new G4Trd(Form("SMP1_%d",i), SMPW1, SMPW2, MPTh/2.0, MPTh/2.0, SMPL);
    G4Trd *SMP2 = new G4Trd(Form("SMP2_%d",i), SMPW1, SMPW2, MPTh/2.0, MPTh/2.0, SMPL);
    G4Box *SSS = new G4Box(Form("SSS_%d",i),SSTh/2.0,MPTh/2.0,SMPL);
    
    G4LogicalVolume *SMP1_log = new G4LogicalVolume( SMP1, GetMaterial("Aluminum"), Form("SMP1_%d_log",i) );
    G4LogicalVolume *SMP2_log = new G4LogicalVolume( SMP2, GetMaterial("Aluminum"), Form("SMP2_%d_log",i) );
    G4LogicalVolume *SSS_log = new G4LogicalVolume( SSS,GetMaterial("Stainless_Steel"),Form("SSS_%d_log",i));

    SMP1_log->SetVisAttributes( AlColor );
    SMP2_log->SetVisAttributes( AlColor );
    SSS_log->SetVisAttributes( SteelColor );

    G4RotationMatrix *Srot1_temp = new G4RotationMatrix;
    Srot1_temp->rotateZ(+90*deg);
    Srot1_temp->rotateX(-MPA);
    G4RotationMatrix *Srot2_temp = new G4RotationMatrix;
    Srot2_temp->rotateZ(-90*deg);
    Srot2_temp->rotateX(-MPA);
    G4RotationMatrix *S2rot1_temp = new G4RotationMatrix;
    S2rot1_temp->rotateZ(+135*deg);
    S2rot1_temp->rotateX(-MPA);
    G4RotationMatrix *S2rot2_temp = new G4RotationMatrix;
    S2rot2_temp->rotateZ(-135*deg);
    S2rot2_temp->rotateX(-MPA);
    G4RotationMatrix *S2rot3_temp = new G4RotationMatrix;
    S2rot3_temp->rotateZ(+45*deg);
    S2rot3_temp->rotateX(-MPA);
    G4RotationMatrix *S2rot4_temp = new G4RotationMatrix;
    S2rot4_temp->rotateZ(-45*deg);
    S2rot4_temp->rotateX(-MPA);

    
    //Place both mounting plates
    new G4PVPlacement( Srot1_temp, G4ThreeVector(-SMP_xoffset-SMPL*sin(MPA), 0, MPzmin_array[i]+SMPL*cos(MPA)), SMP1_log, Form("SMP1_%d_phys",i), worldlog, false, 0, ChkOverlaps);
    SMP1_log->SetVisAttributes(AlColor);
    new G4PVPlacement( Srot2_temp, G4ThreeVector(SMP_xoffset+SMPL*sin(MPA), 0, MPzmin_array[i]+SMPL*cos(MPA)), SMP2_log, Form("SMP2_%d_phys",i), worldlog, false, 0, ChkOverlaps);
    SMP2_log->SetVisAttributes(AlColor);

    //Place shield supports
    G4double SSX = (SS_xoffset+SMPL*sin(MPA))*sin(45*deg)+SSTh/4.0*sin(45*deg);
    G4double SSY = (SS_xoffset+SMPL*sin(MPA))*sin(45*deg)-SSTh/4.0*sin(45*deg);
    
    new G4PVPlacement(S2rot1_temp,G4ThreeVector(-SSX,SSY,MPzmin_array[i]+SMPL*cos(MPA)),SSS_log,Form("SSS1_phys%d",i),worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(S2rot3_temp,G4ThreeVector(-SSX,-SSY,MPzmin_array[i]+SMPL*cos(MPA)),SSS_log,Form("SSS2_phys%d",i),worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(S2rot4_temp,G4ThreeVector(SSX,-SSY,MPzmin_array[i]+SMPL*cos(MPA)),SSS_log,Form("SSS3_phys%d",i),worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(S2rot2_temp,G4ThreeVector(SSX,SSY,MPzmin_array[i]+SMPL*cos(MPA)),SSS_log,Form("SSS4_phys%d",i),worldlog,false,0,ChkOverlaps);

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
    G4LogicalVolume *IM__log = new G4LogicalVolume( IM_, GetMaterial("Iron"), name ); //Solid inner cones
    
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
      //ring_log->SetVisAttributes( G4Colour::Green() );

      
      
      name = cname;
      name += "_phys";
      
      Z = 0.5*(zstart + zstop);
      new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), ring_log, name, worldlog, false, 0  , ChkOverlaps); //Conical rings
    }
  }
  
  if(fDetCon->fBeamlineConf!=2){
    //G4double dz_spool_piece = z_conic_vacline_weldment - z_spool_piece;

    //SSeeds - replacing target-proximal bellows for GMn 2020. Nominally section 1 (S1) and section 2 (S2). First element closest to target.

    //SECTION 1
    //Build cylinder 1 (C1)
    Rin = 1.969*inch;
    Rout = 2.165*inch;
    Thick = 3.966*inch;
    
    G4Tubs *S1_C1 = new G4Tubs("S1_C1",Rin,Rout,Thick/2.0,0.0,twopi);
    G4Tubs *S1_C1_vac = new G4Tubs("S1_C1_vac",0.0,Rin,Thick/2.0,0.0,twopi);
    
    G4LogicalVolume *S1_C1_log = new G4LogicalVolume(S1_C1,GetMaterial("Aluminum"),"S1_C1_log"); //Only aluminum element in S1
    G4LogicalVolume *S1_C1_vaclog = new G4LogicalVolume(S1_C1_vac,GetMaterial("Vacuum"),"S1_C1_vaclog");
    
    S1_C1_log->SetVisAttributes(AlColor);

    Z = 23.130*inch+Thick/2.0; //Initial placement from JT and adding length to central placment of cylinder
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C1_log,"S1_C1_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C1_vaclog,"S1_C1_vacPhys",worldlog,false,0,ChkOverlaps);
    
    //Build cylinder 2 (C2)
    Z += Thick/2.0; //Shift z placement by the length of the first cylinder element
    Rin = 1.880*inch;
    Rout = 2.985*inch;
    Thick = 1.607*inch;
    
    G4Tubs *S1_C2 = new G4Tubs("S1_C2",Rin,Rout,Thick/2.0,0.0,twopi);
    G4Tubs *S1_C2_vac = new G4Tubs("S1_C2_vac",0.0,Rin,Thick/2.0,0.0,twopi);
    
    G4LogicalVolume *S1_C2_log = new G4LogicalVolume(S1_C2,GetMaterial("Stainless_Steel"),"S1_C2_log");
    G4LogicalVolume *S1_C2_vaclog = new G4LogicalVolume(S1_C2_vac,GetMaterial("Vacuum"),"S1_C2_vaclog");
    
    S1_C2_log->SetVisAttributes( SteelColor);

    Z += Thick/2.0; //Shift by half of next cylinder as Thick is redefined
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C2_log,"S1_C2_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C2_vaclog,"S1_C2_vacPhys",worldlog,false,0,ChkOverlaps);
       
    //Build cylinder 3 (C3)
    Z += Thick/2.0; 
    Rin = 1.880*inch;
    Rout = 2.0*inch;
    Thick = 2.255*inch;   
    
    G4Tubs *S1_C3 = new G4Tubs("S1_C3",Rin,Rout,Thick/2.0,0.0,twopi);
    G4Tubs *S1_C3_vac = new G4Tubs("S1_C3_vac",0.0,Rin,Thick/2.0,0.0,twopi);
    
    G4LogicalVolume *S1_C3_log = new G4LogicalVolume(S1_C3,GetMaterial("Stainless_Steel"),"S1_C3_log");
    G4LogicalVolume *S1_C3_vaclog = new G4LogicalVolume(S1_C3_vac,GetMaterial("Vacuum"),"S1_C3_vaclog");
    
    S1_C3_log->SetVisAttributes(SteelColor);

    Z += Thick/2.0; 
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C3_log,"S1_C3_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C3_vaclog,"S1_C3_vacPhys",worldlog,false,0,ChkOverlaps);
       
    //Build cylinder 4 (C4)
    Z += Thick/2.0; 
    Rin = 1.905*inch;
    Rout = 2.985*inch;
    Thick = 0.889*inch;
    
    G4Tubs *S1_C4 = new G4Tubs("S1_C4",Rin,Rout,Thick/2.0,0.0,twopi);
    G4Tubs *S1_C4_vac = new G4Tubs("S1_C4_vac",0.0,Rin,Thick/2.0,0.0,twopi);
    
    G4LogicalVolume *S1_C4_log = new G4LogicalVolume(S1_C4,GetMaterial("Stainless_Steel"),"S1_C4_log");
    G4LogicalVolume *S1_C4_vaclog = new G4LogicalVolume(S1_C4_vac,GetMaterial("Vacuum"),"S1_C4_vaclog");
    
    S1_C4_log->SetVisAttributes(SteelColor);

    Z += Thick/2.0; 
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C4_log,"S1_C4_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C4_vaclog,"S1_C4_vacPhys",worldlog,false,0,ChkOverlaps);
       
    //Build cylinder 5 (C5)
    Z += Thick/2.0;
    Rin = 1.969*inch;
    Rout = 4.375*inch;
    Thick = 0.571*inch;
    
    G4Tubs *S1_C5 = new G4Tubs("S1_C5",Rin,Rout,Thick/2.0,0.0,twopi);
    G4Tubs *S1_C5_vac = new G4Tubs("S1_C5_vac",0.0,Rin,Thick/2.0,0.0,twopi);
    
    G4LogicalVolume *S1_C5_log = new G4LogicalVolume(S1_C5,GetMaterial("Stainless_Steel"),"S1_C5_log");
    G4LogicalVolume *S1_C5_vaclog = new G4LogicalVolume(S1_C5_vac,GetMaterial("Vacuum"),"S1_C5_vaclog");
    
    S1_C5_log->SetVisAttributes(SteelColor);

    Z += Thick/2.0; 
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C5_log,"S1_C5_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C5_vaclog,"S1_C5_vacPhys",worldlog,false,0,ChkOverlaps);
       
    //Build Box 1 (B1)
    Z += Thick/2.0; 
    G4double Lx = 7.756*inch;
    G4double Ly = 13.110*inch;
    G4double Box_th = 0.0985*inch; //Hollow steel box thickness
    Thick = 2.126*inch;
    
    G4Box *S1_B1 = new G4Box("S1_B1",Lx/2.0,Ly/2.0,Thick/2.0);
    G4Box *S1_B1_vac = new G4Box("S1_B1_vac",Lx/2.0-Box_th,Ly/2.0-Box_th,Thick/2.0-Box_th);
    G4Tubs *S1_B1_hole = new G4Tubs("S1_B1_hole",0,Rin,Thick/2.0,0.0,twopi);
    G4Tubs *S1_B1_vac2 = new G4Tubs("S1_B1_hole",0,Rin,Box_th,0.0,twopi);
    
    G4SubtractionSolid *S1_B1_sub1 = new G4SubtractionSolid("S1_B1_sub1",S1_B1,S1_B1_vac,0,G4ThreeVector(0,0,0));
    G4SubtractionSolid *S1_B1_sub2 = new G4SubtractionSolid("S1_B1_sub2",S1_B1_sub1,S1_B1_hole,0,G4ThreeVector(0,Ly/4.0,0));  

    G4LogicalVolume *S1_B1_sublog = new G4LogicalVolume(S1_B1_sub2,GetMaterial("Stainless_Steel"),"S1_B1_log");
    G4LogicalVolume *S1_B1_vaclog = new G4LogicalVolume(S1_B1_vac,GetMaterial("Vacuum"),"S1_B1_vaclog");
    G4LogicalVolume *S1_B1_vac2log = new G4LogicalVolume(S1_B1_vac2,GetMaterial("Vacuum"),"S1_B1_vac2log");
    
    S1_B1_sublog->SetVisAttributes(SteelColor);

    Z += Thick/2.0; 
    new G4PVPlacement(0,G4ThreeVector(X,-Ly/4.0,Z),S1_B1_sublog,"S1_B1_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,-Ly/4.0,Z),S1_B1_vaclog,"S1_B1_vacPhys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z-Thick/2.0+Box_th/2.0),S1_B1_vac2log,"S1_B1_vac2PhysUS",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z+Thick/2.0-Box_th/2.0),S1_B1_vac2log,"S1_B1_vac2PhysDS",worldlog,false,0,ChkOverlaps);
    
    //Build cylinder 6 (C6)
    Z += Thick/2.0; 
    Rin = 1.969*inch;
    Rout = 4.375*inch;
    Thick = 0.571*inch;
    
    G4Tubs *S1_C6 = new G4Tubs("S1_C6",Rin,Rout,Thick/2.0,0.0,twopi);
    G4Tubs *S1_C6_vac = new G4Tubs("S1_C6_vac",0.0,Rin,Thick/2.0,0.0,twopi);
    
    G4LogicalVolume *S1_C6_log = new G4LogicalVolume(S1_C6,GetMaterial("Stainless_Steel"),"S1_C6_log");
    G4LogicalVolume *S1_C6_vaclog = new G4LogicalVolume(S1_C6_vac,GetMaterial("Vacuum"),"S1_C6_vaclog");
    
    S1_C6_log->SetVisAttributes( SteelColor);

    Z += Thick/2.0; 
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C6_log,"S1_C6_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C6_vaclog,"S1_C6_vacPhys",worldlog,false,0,ChkOverlaps);
       
    //Build cylinder 7 (C7)
    Z += Thick/2.0; 
    Rin = 1.960*inch;
    Rout = 2.985*inch;
    Thick = 0.821*inch;
    
    G4Tubs *S1_C7 = new G4Tubs("S1_C7",Rin,Rout,Thick/2.0,0.0,twopi);
    G4Tubs *S1_C7_vac = new G4Tubs("S1_C7_vac",0.0,Rin,Thick/2.0,0.0,twopi);
    
    G4LogicalVolume *S1_C7_log = new G4LogicalVolume(S1_C7,GetMaterial("Stainless_Steel"),"S1_C7_log");
    G4LogicalVolume *S1_C7_vaclog = new G4LogicalVolume(S1_C7_vac,GetMaterial("Vacuum"),"S1_C7_vaclog");
    
    S1_C7_log->SetVisAttributes( SteelColor);

    Z += Thick/2.0; 
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C7_log,"S1_C7_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C7_vaclog,"S1_C7_vacPhys",worldlog,false,0,ChkOverlaps);
       
    //Build cylinder 8 (C8)
    Z += Thick/2.0; 
    Rin = 2.335*inch;
    Rout = 2.375*inch;
    Thick = 4.380*inch;
    
    G4Tubs *S1_C8 = new G4Tubs("S1_C8",Rin,Rout,Thick/2.0,0.0,twopi);
    G4Tubs *S1_C8_vac = new G4Tubs("S1_C8_vac",0.0,Rin,Thick/2.0,0.0,twopi);
    
    G4LogicalVolume *S1_C8_log = new G4LogicalVolume(S1_C8,GetMaterial("Stainless_Steel"),"S1_C8_log");
    G4LogicalVolume *S1_C8_vaclog = new G4LogicalVolume(S1_C8_vac,GetMaterial("Vacuum"),"S1_C8_vaclog");
    
    S1_C8_log->SetVisAttributes( SteelColor);

    Z += Thick/2.0; 
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C8_log,"S1_C8_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C8_vaclog,"S1_C8_vacPhys",worldlog,false,0,ChkOverlaps);
       
    //Build cylinder 9 (C9)
    Z += Thick/2.0; 
    Rin = 1.960*inch;
    Rout = 2.985*inch;
    Thick = 1.7*inch;
    
    G4Tubs *S1_C9 = new G4Tubs("S1_C9",Rin,Rout,Thick/2.0,0.0,twopi);
    G4Tubs *S1_C9_vac = new G4Tubs("S1_C9_vac",0.0,Rin,Thick/2.0,0.0,twopi);
    
    G4LogicalVolume *S1_C9_log = new G4LogicalVolume(S1_C9,GetMaterial("Stainless_Steel"),"S1_C9_log");
    G4LogicalVolume *S1_C9_vaclog = new G4LogicalVolume(S1_C9_vac,GetMaterial("Vacuum"),"S1_C9_vaclog");
    
    S1_C9_log->SetVisAttributes( SteelColor);

    Z += Thick/2.0; 
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C9_log,"S1_C9_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C9_vaclog,"S1_C9_vacPhys",worldlog,false,0,ChkOverlaps);
       
    //Build cylinder 10 (C10)
    Z += Thick/2.0; 
    Rin = 1.880*inch;
    Rout = 2.0*inch;
    Thick = 1.081*inch;
    
    G4Tubs *S1_C10 = new G4Tubs("S1_C10",Rin,Rout,Thick/2.0,0.0,twopi);
    G4Tubs *S1_C10_vac = new G4Tubs("S1_C10_vac",0.0,Rin,Thick/2.0,0.0,twopi);
    
    G4LogicalVolume *S1_C10_log = new G4LogicalVolume(S1_C10,GetMaterial("Stainless_Steel"),"S1_C10_log");
    G4LogicalVolume *S1_C10_vaclog = new G4LogicalVolume(S1_C10_vac,GetMaterial("Vacuum"),"S1_C10_vaclog");
    
    S1_C10_log->SetVisAttributes( SteelColor);

    Z += Thick/2.0; 
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C10_log,"S1_C10_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S1_C10_vaclog,"S1_C10_vacPhys",worldlog,false,0,ChkOverlaps);


    //SECTION 2 (S2)
    //This section is layered. Will build in order of proximity to target, then proximity to z axis.

    //Build cylinder 1 (C1)
    Z += Thick/2.0;
    Rin = 1.880*inch;
    Rout = 2.0*inch;
    Thick = 18.843*inch;
    
    G4Tubs *S2_C1 = new G4Tubs("S2_C1",Rin,Rout,Thick/2.0,0.0,twopi);
    G4Tubs *S2_C1_vac = new G4Tubs("S2_C1_vac",0.0,Rin,Thick/2.0,0.0,twopi);
    
    G4LogicalVolume *S2_C1_log = new G4LogicalVolume(S2_C1,GetMaterial("Stainless_Steel"),"S2_C1_log");
    G4LogicalVolume *S2_C1_vaclog = new G4LogicalVolume(S2_C1_vac,GetMaterial("Vacuum"),"S2_C1_vaclog");
    
    S2_C1_log->SetVisAttributes(SteelColor);

    Z += Thick/2.0; 
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S2_C1_log,"S2_C1_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S2_C1_vaclog,"S2_C1_vacPhys",worldlog,false,0,ChkOverlaps);
    
    //Build cylinder 2 (C2)
    Z -= Thick/2.0; //Starting back where cylinder 1 begins
    Rin = 2.02*inch;
    Rout = 2.27*inch;
    Thick = 18.5*inch;
    
    G4Tubs *S2_C2 = new G4Tubs("S2_C2",Rin,Rout,Thick/2.0,0.0,twopi);
    G4Tubs *S2_C2_vac = new G4Tubs("S2_C2_vac",0.0,Rin,Thick/2.0,0.0,twopi);
    
    G4LogicalVolume *S2_C2_log = new G4LogicalVolume(S2_C2,GetMaterial("Iron"),"S2_C2_log");
    G4LogicalVolume *S2_C2_vaclog = new G4LogicalVolume(S2_C2_vac,GetMaterial("Vacuum"),"S2_C2_vaclog");
    
    S2_C2_log->SetVisAttributes(ironColor);

    Z += Thick/2.0; 
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S2_C2_log,"S2_C2_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S2_C2_vaclog,"S2_C2_vacPhys",worldlog,false,0,ChkOverlaps);

    //Build Rings
    Z -= Thick/2.0; //Starting again back where cylinder 1 begins
    G4double initZ = Z + 0.438*inch;  //Adding offset to first ring
    G4double S2ringTh = 0.5*inch;
    G4double S2ringr = 2.651*inch;
    G4double S2ringL = 1.625*inch;
    G4double S2ringSep = 0.375*inch;
    G4double ringPlac;
    
    for (int i=0; i<9; i++){
      ringPlac = initZ+i*S2ringL+i*S2ringSep+S2ringL/2.0;
      
      G4Tubs *S2ring = new G4Tubs(Form("S2ring%d",i+1),S2ringr,S2ringr+S2ringTh,S2ringL/2.0,0.0,twopi);
      
      G4LogicalVolume *S2ringLog = new G4LogicalVolume(S2ring,GetMaterial("Iron"),Form("S2ring%d_log",i+1));
      
      new G4PVPlacement(0,G4ThreeVector(X,Y,ringPlac),S2ringLog,Form("S2ring%dLog_pv",i+1),worldlog,false,0,ChkOverlaps);
      
      S2ringLog->SetVisAttributes(ironColor);
    }
    
    //Build Shields and Mounts
    G4double B_Lx = 0.25*inch;
    G4double B1_Ly = 5.848*inch;
    G4double B2_Ly = 2.5*inch;
    G4double B_Lz = 17.625*inch;
    
    G4RotationMatrix *S2rot1_temp = new G4RotationMatrix;
    S2rot1_temp->rotateZ(45*deg);
    G4RotationMatrix *S2rot2_temp = new G4RotationMatrix;
    S2rot2_temp->rotateZ(-45*deg);
    
    G4double B1offset = 9.620*inch/2.0;
    G4double B2offset = 6.802*inch/2.0;
    
    G4Box *S2_B1 = new G4Box("S2_B1",B_Lx/2.0,B1_Ly/2.0,B_Lz/2.0);
    G4Box *S2_B2 = new G4Box("S2_B2",B_Lx/2.0,B2_Ly/2.0,B_Lz/2.0);
    
    G4LogicalVolume *S2_B1_log = new G4LogicalVolume(S2_B1,GetMaterial("Aluminum"),"S2_B1_log");
    G4LogicalVolume *S2_B2_log = new G4LogicalVolume(S2_B2,GetMaterial("Stainless_Steel"),"S2_B2_log");
    
    //Larger aluminum mounts
    new G4PVPlacement(0,G4ThreeVector(B1offset+B_Lx/2.0,Y,initZ+B_Lz/2.0),S2_B1_log,"S2_B1_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(-B1offset-B_Lx/2.0,Y,initZ+B_Lz/2.0),S2_B1_log,"S2_B1_phys",worldlog,false,0,ChkOverlaps);
    
    X = (B2offset+B_Lx/2.0)*sin(45*deg);
    Y = (B2offset+B_Lx/2.0)*sin(45*deg);
    
    //Smaller steel shields
    new G4PVPlacement(S2rot1_temp,G4ThreeVector(-X,Y,initZ+B_Lz/2.0),S2_B2_log,"S2_B1_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(S2rot2_temp,G4ThreeVector(-X,-Y,initZ+B_Lz/2.0),S2_B2_log,"S2_B1_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(S2rot1_temp,G4ThreeVector(X,-Y,initZ+B_Lz/2.0),S2_B2_log,"S2_B1_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(S2rot2_temp,G4ThreeVector(X,Y,initZ+B_Lz/2.0),S2_B2_log,"S2_B1_phys",worldlog,false,0,ChkOverlaps);
    
    S2_B1_log->SetVisAttributes(AlColor);
    S2_B2_log->SetVisAttributes(SteelColor);

    X = 0.0;
    Y = 0.0;
   
    //Build cylinder 3 (C3)
    Z += 18.843*inch; //Shifting Z by length of the first cylinder
    Rin = 1.880*inch;
    Rout = 3.375*inch;
    Thick = 1.701*inch/2.0;
    
    G4Tubs *S2_C3 = new G4Tubs("S2_C3",Rin,Rout,Thick/2.0,0.0,twopi);
    G4Tubs *S2_C3_vac = new G4Tubs("S2_C3_vac",0.0,Rin,Thick/2.0,0.0,twopi);
    
    G4LogicalVolume *S2_C3_log = new G4LogicalVolume(S2_C3,GetMaterial("Stainless_Steel"),"S2_C3_log");
    G4LogicalVolume *S2_C3_vaclog = new G4LogicalVolume(S2_C3_vac,GetMaterial("Vacuum"),"S2_C3_vaclog");
    
    S2_C3_log->SetVisAttributes(SteelColor);

    Z += Thick/2.0; 
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S2_C3_log,"S2_C3_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S2_C3_vaclog,"S2_C3_vacPhys",worldlog,false,0,ChkOverlaps);
    
    //Build cylinder 4 (C4)
    Z += Thick/2.0;
    Rin = 1.750*inch;
    Rout = 3.375*inch;
    Thick = 1.701*inch/2.0;
    
    G4Tubs *S2_C4 = new G4Tubs("S2_C4",Rin,Rout,Thick/2.0,0.0,twopi);
    G4Tubs *S2_C4_vac = new G4Tubs("S2_C4_vac",0.0,Rin,Thick/2.0,0.0,twopi);
    
    G4LogicalVolume *S2_C4_log = new G4LogicalVolume(S2_C4,GetMaterial("Stainless_Steel"),"S2_C4_log");
    G4LogicalVolume *S2_C4_vaclog = new G4LogicalVolume(S2_C4_vac,GetMaterial("Vacuum"),"S2_C4_vaclog");
    
    S2_C4_log->SetVisAttributes(SteelColor);

    Z += Thick/2.0; 
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S2_C4_log,"S2_C4_phys",worldlog,false,0,ChkOverlaps);
    new G4PVPlacement(0,G4ThreeVector(X,Y,Z),S2_C4_vaclog,"S2_C4_vacPhys",worldlog,false,0,ChkOverlaps);
    
  }
    
  if(fDetCon->fBeamlineConf==1){
    
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

  /*
  switch(fDetCon->fBeamlineConf){
  case(1):// reminder: beamline config 1 = GEp
    z_Magnets_array.push_back( z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY );
    //z_Magnets_array.push_back( z_formed_bellows + UpstreamCoilDepth/2.0 + UpstreamCoilThickY );
    z_Magnets_array.push_back( z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0 );
    z_Magnets_array.push_back( z_formed_bellows + 76.09*inch + 1.71*inch + DownstreamYokeDepth/2.0 );
    z_Magnets_array.push_back( z_formed_bellows + 76.09*inch + DS_coil_ThickY + DS_coil_depth/2.0 );
    
    //z_Magnets_array = {z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY, z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0, z_formed_bellows + 76.09*inch + 1.71*inch + DownstreamYokeDepth/2.0, z_formed_bellows + 76.09*inch + DS_coil_ThickY + DS_coil_depth/2.0};
     break;
  // case(2):// reminder: beamline config 2 = GEn, SIDIS ?
  //   Z = z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY;
  //   break;
  case(3):// reminder: beamline config 3 = GMn
    
    //z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 15.94)*inch + UpstreamCoilDepth/2.0 );
    //z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 15.94 + 8.3 - 6.47)*inch - UpstreamCoilThickY + YokeRightZFinal/2.0 );
    //z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 85.78)*inch + DownstreamYokeDepth/2.0 );
    //z_Magnets_array.push_back( z_conic_vacline_weldment + (0.84 + 0.14 + 85.78 - 1.71)*inch + DS_coil_ThickY + DS_coil_depth/2.0 );
    
    //SSeeds 2020 update
    z_Magnets_array.push_back( (79.724)*inch + UpstreamCoilDepth/2.0 );
    z_Magnets_array.push_back( (79.724 + 8.3 - 6.47)*inch - UpstreamCoilThickY + YokeRightZFinal/2.0 );
    z_Magnets_array.push_back( (149.559)*inch + DownstreamYokeDepth/2.0 );
    z_Magnets_array.push_back( (149.559 - 1.71)*inch + DS_coil_ThickY + DS_coil_depth/2.0 );
    
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
*/
  //SSeeds 2020 update
  z_Magnets_array.push_back( (79.724)*inch + UpstreamCoilDepth/2.0 );
  z_Magnets_array.push_back( (79.724 + 8.3 - 6.47)*inch - UpstreamCoilThickY + YokeRightZFinal/2.0 );
  z_Magnets_array.push_back( (149.559)*inch + DownstreamYokeDepth/2.0 );
  z_Magnets_array.push_back( (149.559 - 1.71)*inch + DS_coil_ThickY + DS_coil_depth/2.0 );
    
  G4Box *UpstreamCoil_outer = new G4Box("UpstreamCoil_outer", UpstreamCoilThickX/2.0, (UpstreamCoilHeight+2.0*UpstreamCoilThickY)/2.0, (UpstreamCoilDepth + 2.0*UpstreamCoilThickY)/2.0 );
  G4Box *UpstreamCoil_inner = new G4Box("UpstreamCoil_inner", UpstreamCoilThickX/2.0 + cm, UpstreamCoilHeight/2.0, UpstreamCoilDepth/2.0 );

  G4SubtractionSolid *UpstreamCoil = new G4SubtractionSolid( "UpstreamCoil", UpstreamCoil_outer, UpstreamCoil_inner, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *UpstreamCoil_log = new G4LogicalVolume(UpstreamCoil, GetMaterial("Copper"), "UpstreamCoil_log" );

  UpstreamCoil_log->SetVisAttributes( CopperColor );

  Z = z_Magnets_array[0]-0.15*inch;//z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY; (- (UpstreamCoilDepth-UpstreamPoleDepth)/2 - very convoluted dimensions here)
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
  Z = z_Magnets_array[1]-0.15*inch;//z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0;
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

  //Sseeds - Not sure if this is still relevant
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
    Electronics_log->SetVisAttributes(G4VisAttributes::GetInvisible());
    G4String GEMElectronicsname = "BLneutronDet";
    G4String  GEMElectronicscollname = "BLneutronDet";
    G4SBSCalSD *GEMElecSD = NULL;

    GEMElectronicsname += "GMn";
    GEMElectronicscollname += "GMn";
    
    /*
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
    */

    
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

    //SSJT - commented as this object does not appear in jt file
    
    //new G4PVPlacement( 0, G4ThreeVector(x_blndet, y_blndet, z_blndet), Electronics_log, "GMn_Electronics", worldlog, false, 0 , ChkOverlaps);
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
  // FVL1_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // FVL2_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // FVL3_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // FVL5_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // FVL6_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // FVL7_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // TVB1_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // TVL8_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // TVL9_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // TMV9_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // TTV1_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // TTV2_log->SetVisAttributes( G4VisAttributes::GetInvisible() );

  /*
  //SSeeds 12.17.20 - test to see where common exit beamline connects with target to midpipe section. Test ring marks beginning of target to midpipe section according to JT file - Dec 2020
  G4double P1testRing2_rin = 14.75*inch; 
  G4double P1testRing2_rou = 15.0*inch; //CJT 13.0 - making larger for debug
  G4double P1testRing2_L = 0.187/10*inch;
  G4double P1placement = 43.097*inch; //Beginning S.2
  //G4double P2placement = 23.75*inch; //Radius of scattering chamber
  G4double P2placement = 23.130*inch; //Start of first cylinder, beginning S.1
  //G4double P3placement = 63.331*inch; //Beginning S.3 (and conical beampipe section)
  G4double P3placement = 63.641*inch; //Beginning S.3 (and conical beampipe section)
  G4double P4placement = 89.782*inch; //Beginning S.4 rings
  G4double P5placement = 165.782*inch; //Beginning S.5 rings
  //G4double P6placement = 201.132*inch; //Beginning exit bellows
  //G4double P6placement = 167.782*inch+31.591*inch; //Beginning exit bellows
  G4double P6placement = 200.512*inch; 
  G4double P7placement = 212.37*inch; //Beginning target to midpipe (green)
  G4double P8placement = 0.0*inch; //Target Center

  G4double P9placement = 63.782*inch; //Checking ring placement
  G4double P10placement = 75.407*inch; //End of S.3 rings
  G4double P11placement = 79.724*inch; //Face of first corrector magnet
  G4double P12placement = 149.559*inch; //Face of second corrector magnet
  

  G4Tubs *P1testRing2 = new G4Tubs("P1testRing2", P1testRing2_rin, P1testRing2_rou, P1testRing2_L, 0.*deg, 360.*deg);
  G4Tubs *P1testRing3 = new G4Tubs("P1testRing3", P1testRing2_rin, P1testRing2_rou, P1testRing2_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1testRing2Log = new G4LogicalVolume(P1testRing2, GetMaterial("Air"), "P1testRing2_log", 0, 0, 0);
  G4LogicalVolume *P1testRing3Log = new G4LogicalVolume(P1testRing3, GetMaterial("Air"), "P1testRing3_log", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P1placement), P1testRing2Log, "P1testRing2Log_pv", worldlog, false, 0, ChkOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P2placement), P1testRing2Log, "P1testRing2Log_pv", worldlog, false, 0, ChkOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P3placement), P1testRing2Log, "P1testRing2Log_pv", worldlog, false, 0, ChkOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P4placement), P1testRing2Log, "P1testRing2Log_pv", worldlog, false, 0, ChkOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5placement), P1testRing2Log, "P1testRing2Log_pv", worldlog, false, 0, ChkOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P6placement), P1testRing2Log, "P1testRing2Log_pv", worldlog, false, 0, ChkOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P7placement), P1testRing2Log, "P1testRing2Log_pv", worldlog, false, 0, ChkOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P8placement), P1testRing2Log, "P1testRing2Log_pv", worldlog, false, 0, ChkOverlaps);
  
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P9placement), P1testRing3Log, "P1testRing3Log_pv", worldlog, false, 0, ChkOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P10placement), P1testRing3Log, "P1testRing3Log_pv", worldlog, false, 0, ChkOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P11placement), P1testRing3Log, "P1testRing3Log_pv", worldlog, false, 0, ChkOverlaps);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P12placement), P1testRing3Log, "P1testRing3Log_pv", worldlog, false, 0, ChkOverlaps);

  P1testRing2Log->SetVisAttributes( G4Colour::Red()); //Debug

  P1testRing3Log->SetVisAttributes( G4Colour::Green()); //Debug
  */
  
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
  
  G4double TargetCenter_zoffset = 6.5*inch;
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
  
  //MakeCommonExitBeamline(worldlog);

  //////SSeeds Jan. 2021 - Start, Temporary awaiting hall dimensions from engineering team.//////

  
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
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange1_vac_log, "CVLW_Flange1_vac_phys", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange1_log, "CVLW_Flange1_phys", worldlog, false, 0 , ChkOverlaps );

  //conic vacuum line weldment:
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
  G4Tubs *CVLW_Flange2 = new G4Tubs("CVLW_Flange2", Rin, Rout, Thick/2.0, 0.0, twopi );
  // Fill with vacuum
  Rin = 0.0;
  Rout1 = Rin1;
  Rout2 = Rin2;
  Rout = 10.734*inch/2.0;
  
  //G4Cons *CVLW_Flange2_vac = new G4Cons("CVLW_Flange2_vac", Rin, Rout1, Rin, Rout2, Thick/2.0, 0.0, twopi );
  G4Tubs *CVLW_Flange2_vac = new G4Tubs( "CVLW_Flange2_vac", Rin, Rout, Thick/2.0, 0.0, twopi );

  G4LogicalVolume *CVLW_Flange2_log = new G4LogicalVolume( CVLW_Flange2, GetMaterial("Stainless_Steel"), "CVLW_Flange2_log" );
  G4LogicalVolume *CVLW_Flange2_vac_log = new G4LogicalVolume( CVLW_Flange2_vac, GetMaterial("Vacuum"), "CVLW_Flange2_vac_log" );

  CVLW_Flange2_log->SetVisAttributes(SteelColor );
  CVLW_Flange2_vac_log->SetVisAttributes( Vacuum_visatt );

  // Then place the vacuum inside the Iron Tube
  Z = z_welded_bellows - Thick/2.0;
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange2_vac_log, "CVLW_Flange2_vac_phys", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), CVLW_Flange2_log, "CVLW_Flange2_phys", worldlog, false, 0 , ChkOverlaps );

  //Next: "Welded bellows"
  //G4double dz_welded_bellows = 207.144*inch - z_welded_bellows - TargetCenter_zoffset; // = =5.512 inches
  G4double dz_welded_bellows = 212.37*inch - TargetCenter_zoffset;
  
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
  //Thick = 212.37*inch - Z + 6.5*inch; //SSeeds 2021, Temporary extension of vacuum  
  G4Tubs *WB_Vacuum = new G4Tubs( "WB_Vacuum", Rin, Rout, Thick/2.0, 0.0, twopi );

  G4LogicalVolume *WB_Vacuum_log = new G4LogicalVolume(WB_Vacuum, GetMaterial("Vacuum"), "WB_Vacuum_log" );

  //WB_Vacuum_log->SetVisAttributes( Vacuum_visatt );
  
  Z = z_welded_bellows + dz_welded_bellows/2.0;

  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), WB_Vacuum_log, "WB_Vacuum_phys", worldlog, false, 0 , ChkOverlaps );


  G4double P7placement = 212.37*inch; //Beginning target to midpipe (green)
  

  //Magnetic Tubes
  G4int Ndivs;// number of segments with shielding

  std::vector<G4double> Rin_array;// radii for inner shielding elements
  std::vector<G4double> Zin_array;// z for inner shielding elements
  std::vector<G4int> Nrings_out;// number of outer elements per segments
  std::vector<G4double> Rout_array;// radii for inner shielding elements
  std::vector<G4double> Zout_array;// z for inner shielding elements

  G4double OMthick = 1.625*inch;
  G4double OMspace = 0.375*inch;


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
      
    
  //Target-proximal Bellows and Flanges
  
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

  Thick = dz_iron_tubes/2.0; //Sseeds 2021 - extensions of proximal tubing
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
  
  //Next, corrector magnets:
  //Define some dimensions that are going to be useful to define the distances
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

  std::vector<G4double> z_Magnets_array;

  z_Magnets_array.push_back( z_formed_bellows + 6.47*inch + UpstreamCoilDepth/2.0 + UpstreamCoilThickY );
  z_Magnets_array.push_back( z_formed_bellows + 8.3*inch + YokeRightZFinal/2.0 );
  z_Magnets_array.push_back( z_formed_bellows + 76.09*inch + 1.71*inch + DownstreamYokeDepth/2.0 );
  z_Magnets_array.push_back( z_formed_bellows + 76.09*inch + DS_coil_ThickY + DS_coil_depth/2.0 );

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

  G4Box *YokeTopPiece = new G4Box("YokeTopPiece", YokeTopPiece_Width/2.0, YokeTopPiece_Height/2.0, YokeTopPiece_Depth/2.0 );
  G4LogicalVolume *YokeTopPiece_log = new G4LogicalVolume( YokeTopPiece, GetMaterial("Iron"), "YokeTopPiece_log" );

  YokeTopPiece_log->SetVisAttributes( ironColor );
  
  X = 0.0;
  Y = (11.81*inch + YokeTopPiece_Height)/2.0;

  //two placements of yoke top piece (top and bottom symmetric):
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeTopPiece_log, "UpstreamYokeTop_phys", worldlog, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X,-Y,Z), YokeTopPiece_log, "UpstreamYokeBottom_phys", worldlog, false, 1 , ChkOverlaps );
  
  G4Box *YokeLeftPiece = new G4Box("YokeLeftPiece", YokeLeftPiece_Width/2.0, YokeLeftPiece_Height/2.0, YokeLeftPiece_Depth/2.0 );
  G4LogicalVolume *YokeLeftPiece_log = new G4LogicalVolume( YokeLeftPiece, GetMaterial("Iron"), "YokeLeftPiece_log" );
  YokeLeftPiece_log->SetVisAttributes(ironColor );
  
  X = 7.52*inch + YokeLeftPiece_Width/2.0;
  Y = 0.0;
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeLeftPiece_log, "UpstreamYokeLeftPiece_phys", worldlog, false, 0 , ChkOverlaps );

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

    double x_blndet = 3.0*m;
    double y_blndet = 0.0*m;
    double z_blndet = 2.5*m;
    
    G4double ElecX = 5.0*cm;
    G4double ElecY = 100.0*cm;
    G4double ElecZ = 100.0*cm;
    
    G4Box *Electronics = new G4Box( "Electronics" , ElecX/2.0, ElecY/2.0, ElecZ/2.0);
    G4LogicalVolume *Electronics_log = new G4LogicalVolume( Electronics , GetMaterial("Silicon"), "Electronics_log" );
    Electronics_log->SetVisAttributes(G4VisAttributes::GetInvisible());
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

    new G4PVPlacement( 0, G4ThreeVector(x_blndet, y_blndet, z_blndet), Electronics_log, "GMn_Electronics", worldlog, false, 0 , ChkOverlaps);
  }


  //////SSeeds Jan. 2021 - End, Temporary awaiting hall dimensions from engineering team.//////
  
  MakeBeamExit(worldlog,TargetCenter_zoffset); // Added by D Flay (Sept 2020) 
}

// This is the beam line for GMn
void G4SBSBeamlineBuilder::MakeGMnBeamline(G4LogicalVolume *worldlog){
  bool ChkOverlaps = false;
  double swallrad = 1.143*m/2;
  double swallrad_inner = 1.041/2.0*m; 
  //EFuchey: 2017/02/14: change parameters for Standard scat chamber:
  double sc_entbeampipeflange_dist = 25.375*2.54*cm;// entrance pipe flange distance from hall center. Re-verified by SSeeds for 2020 JT
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

  // EFuchey: 2017/02/14
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-sc_entbeampipeflange_dist), entLog, "ent_phys", worldlog, false,0 , ChkOverlaps);
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-sc_entbeampipeflange_dist), entvacLog, "entvac_phys", worldlog,false,0 , ChkOverlaps);
  //}
   
  MakeCommonExitBeamline(worldlog);

  // Added by D Flay (Sept 2020) 
  G4double inch = 2.54*cm; 
  MakeBeamExit(worldlog,0.0*inch); // account for offset of 6.5" in MakeCommonExitBeamline   
  
  G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
  entvacLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  entLog->SetVisAttributes(pipeVisAtt);
  
  //entvacLog_cut->SetVisAttributes(G4VisAttributes::GetInvisible());
  //entLog_cut->SetVisAttributes(pipeVisAtt);
}


// This is the beam line for 3He
void G4SBSBeamlineBuilder::Make3HeBeamline(G4LogicalVolume *worldlog){  // for GEn, A1n, SIDIS

  //===== UPSTREAM - PIPE =====//
 
  // beam collimator to protect the target (D. Flay study) 
  bool enableBC_dnstr     = fDetCon->GetBeamCollimatorEnable_dnstr();
  bool enableBC_upstr     = fDetCon->GetBeamCollimatorEnable_upstr();
 
  //General Parameters
  G4double inch = 2.54*cm;
  bool ChkOverlaps = false;
  G4VisAttributes *Aluminum = new G4VisAttributes(G4Colour(0.3,0.3,1.0));
  G4VisAttributes *Iron = new G4VisAttributes(G4Colour(0.3,0.3,0.3));
  G4VisAttributes *LeadColor = new G4VisAttributes(G4Colour(0.4,0.4,0.4));
  G4VisAttributes *DebugRed = new G4VisAttributes(G4Colour(1.0,0.,0.));
  G4VisAttributes *Beryllium = new G4VisAttributes(G4Colour(1.0,0.,0.));
  G4VisAttributes *SteelColor = new G4VisAttributes(G4Colour(0.75,0.75,0.75));
  
  //Section Zero
  //This section details all components of the enterance beamline from the target chamber upstream.
  
  //CJT taking distance measurements from the target cylinder, not the end. The radius of the hemisphere in CJT is 0.41". Must subtract this from the total target length for calculation of offsets.
  G4double targetEndOffset_z = 22.8/2.0*inch; //(23.62"-2*0.41")/2
  //G4double P0initPlacement_z = -(23.62/2*inch) - 2.66*inch;  //From updated CJT file. -(half target length) - offset to beampipe flange
  G4double P0initPlacement_z = -(targetEndOffset_z) - 6.395*inch;
  
  //Ring 0A - Be window housing flange. Most proximal piece to target. Bolted to beampipe flange.
  G4double P0ringA_L = 0.510/2.0*inch;
  G4double P0ringA_rin = 1.46/2.0*inch; //This is the inner radius of the proximal upstream pipe per dwg no 67507-0023.
  G4double P0ringA_rou = 1.37*inch; //CJT

  G4Tubs *P0ringA = new G4Tubs("P0ringA", P0ringA_rin, P0ringA_rou, P0ringA_L, 0.0*deg, 360.0*deg);

  G4double P0disk_cut_L = 0.020/2.0*inch; //0.030" (Nominal cut)  - 0.010" (Thickness of window). Ignoring any Be beyond the inner radius of the ring and setting it flush with the face.
  G4double P0disk_cut_rou = 2.105/2.0*inch;

  G4Tubs *P0disk_cut = new G4Tubs("P0disk_cut", 0.0, P0disk_cut_rou, P0disk_cut_L, 0.0*deg, 360.0*deg );
  
  G4RotationMatrix *P0disk_cut_rot = new G4RotationMatrix; //No rotation necessary
  
  G4SubtractionSolid *P0ringA_cut = new G4SubtractionSolid( "P0ringA_cut", P0ringA, P0disk_cut, P0disk_cut_rot, G4ThreeVector( 0.0, 0.0, P0ringA_L - P0disk_cut_L));

  G4LogicalVolume *P0ringA_cutLog = new G4LogicalVolume(P0ringA_cut, GetMaterial("Aluminum"), "P0ringA_cut_log", 0, 0, 0);

  //Place ring 0A
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P0initPlacement_z-P0ringA_L+2*P0ringA_L), P0ringA_cutLog, "P0ringA_cutLog_pv", worldlog, false, 0 , ChkOverlaps ); //Moving forward by one length to attach to proximal beampipe flange

  //Be window 0A - implementing Be "Dome" window upstream from dwg. no. 67507-0023
  G4double P0shell_w = 0.010*inch;
  G4double P0shell_r = 2.85*inch;
  G4double P0dome_r = 1.46/2.0*inch;
  G4double P0dome_vd = 0.096*inch; //max depth into beampipe of dome
  G4double P0dome_th = atan(P0dome_r/(P0shell_r-P0dome_vd));
  
  G4RotationMatrix *P0dome_rot = new G4RotationMatrix;
  P0dome_rot->rotateX(-180.0*deg);
  
  //Construct sphere with theta constraint to produce dome window
  G4Sphere *P0domeA = new G4Sphere("P0domeA", P0shell_r, P0shell_r+P0shell_w, 0.0*deg, 360.0*deg, 0.0*deg, P0dome_th);

  G4LogicalVolume *P0domeALog = new G4LogicalVolume(P0domeA, GetMaterial("Beryllium"), "P0dome_log", 0, 0, 0);

  fDetCon->InsertTargetVolume(P0domeALog->GetName());

  //Place Be window 0A dome
  new G4PVPlacement( P0dome_rot, G4ThreeVector( 0.0, 0.0, P0initPlacement_z+P0shell_r-2.0*P0disk_cut_L-P0dome_vd+2.0*P0ringA_L), P0domeALog, "P0domeALog_pv", worldlog, false, 0, ChkOverlaps ); //Moving forward by one length to attach to proximal beampipe flange
  
  //Flange 0A - must add back the old P0 flange as the actual beampipe flange - flange length same as Be housing flange
  G4double P0flange_rin = 1.46/2.0*inch; //CJT file confirms larger inner radius for pipe.
  G4double P0flange_rou = 1.37*inch; //CJT

  G4Tubs *P0flange = new G4Tubs("P0flange", P0flange_rin, P0flange_rou, P0ringA_L, 0.0*deg, 360*deg);

  G4LogicalVolume *P0flangeLog = new G4LogicalVolume(P0flange, GetMaterial("Aluminum"), "P0flange_log", 0, 0, 0);

  //Place flange 0A
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P0initPlacement_z-P0ringA_L), P0flangeLog, "P0flange_log_pv", worldlog, false, 0 , ChkOverlaps );
  //Flange vacuum added to Vacuum 0A
  
  //Tube 0A
  G4double P0tubeA_L = 3.427/2.0*inch; //CJT
  G4double P0tubeA_rou = 0.846*inch; //CJT

  G4Tubs *P0tubeA = new G4Tubs("P0tubeA", P0ringA_rin, P0tubeA_rou, P0tubeA_L, 0.0*deg, 360.0*deg);

  G4LogicalVolume *P0tubeALog = new G4LogicalVolume(P0tubeA, GetMaterial("Aluminum"), "P0tubeA_log", 0, 0, 0);

  //Place tube 0A
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P0initPlacement_z-2.0*P0ringA_L-P0tubeA_L), P0tubeALog, "P0tubeALog_pv", worldlog, false, 0 , ChkOverlaps );

  //Vacuum 0A
  G4Tubs *P0tubeA_vac = new G4Tubs("P0tubeA_vac", 0.0, P0ringA_rin, 3.0*P0ringA_L+P0tubeA_L-P0disk_cut_L-P0shell_w, 0.*deg, 360.*deg);

  //need to subtract the vacuum volume walled off by the Be dome
  G4Sphere *P0sphere_cut = new G4Sphere("P0sphere_cut", 0.0, P0shell_r+P0shell_w, 0.0*deg, 360.0*deg, 0.0*deg, 360.0*deg);

  G4SubtractionSolid *P0tubeA_winvac = new G4SubtractionSolid("P0tubeA_winvac", P0tubeA_vac, P0sphere_cut, 0, G4ThreeVector(0.0, 0.0, (P0ringA_L+P0tubeA_L-P0disk_cut_L)+(P0shell_r+P0shell_w)-P0dome_vd));

  G4LogicalVolume *P0tubeA_winvacLog = new G4LogicalVolume(P0tubeA_winvac, GetMaterial("Vacuum"), "P0tubeA_winvac_log", 0, 0, 0);

  // edit by D Flay to make it easier to read 
  G4double P0_vac_a_z   = P0initPlacement_z-2.0*P0disk_cut_L-P0shell_w-(P0ringA_L+P0tubeA_L-P0disk_cut_L)+2.0*P0ringA_L; 
  //Place vacuum 0A
  new G4PVPlacement(0, 
	// G4ThreeVector( 0.0, 0.0, P0initPlacement_z-2.0*P0disk_cut_L-P0shell_w-(P0ringA_L+P0tubeA_L-P0disk_cut_L)+2.0*P0ringA_L), 
	G4ThreeVector(0.0,0.0,P0_vac_a_z), 
	P0tubeA_winvacLog,
	"P0tubeA_winvacLog_pv",
	worldlog,
	false,
	0,
	ChkOverlaps );
  

  //Tube 0B
  //G4double P0tubeB_L = 15.303/2.0*inch; //CJT
  G4double P0tubeB_L = 15.303/2.0*inch; //CJT
  G4double P0tubeB_rin = 1.46/2.0*inch; //CJT
  G4double P0tubeB_rou = 0.846*inch; //CJT

  G4Tubs *P0tubeB = new G4Tubs("P0tubeB", P0tubeB_rin, P0tubeB_rou, P0tubeB_L, 0.0*deg, 360.0*deg);

  G4LogicalVolume *P0tubeBLog = new G4LogicalVolume(P0tubeB, GetMaterial("Aluminum"), "P0tubeB_log", 0, 0, 0);

  // edit by D Flay to make it easier to read 
  G4double P0_b_z     = P0initPlacement_z-2.0*P0ringA_L-2.0*P0tubeA_L-P0tubeB_L;  
  G4double P0_vac_b_z = P0_b_z; // for diagnostics   
  //Place tube 0B
  new G4PVPlacement(0,
                    // G4ThreeVector(0.0,0.0,P0initPlacement_z-2.0*P0ringA_L-2.0*P0tubeA_L-P0tubeB_L),
                    G4ThreeVector(0.0,0.0,P0_b_z),
		    P0tubeBLog,
		    "P0tubeBLog_pv",
		    worldlog,
		    false,
		    0,
		    ChkOverlaps);
 
  //Place vacuum 0B
  G4Tubs *P0tubeB_vac             = new G4Tubs("P0tubeB_vac", 0.0, P0tubeB_rin, P0tubeB_L, 0.*deg, 360.*deg);
  G4LogicalVolume *P0tubeB_vacLog = new G4LogicalVolume(P0tubeB_vac, GetMaterial("Vacuum"), "P0tubeB_vac_log", 0, 0, 0);

  // in the event we need to divide up the vacuum 
  G4double lengthBC_dnstr = fDetCon->GetBeamCollimatorL_dnstr();
  G4double z0_bc = P0_vac_b_z;                            // center of original vacuum insert 
  G4double z1_bc = fDetCon->GetBeamCollimatorZ_dnstr();   // center of beam collimator
  G4double L0_bc = 2.*P0tubeB_L; 
  G4double L1_bc = lengthBC_dnstr; 
  // lengths of upstream and downstream elements
  G4double Lu_bc = 0.5*(L0_bc-L1_bc) + ( fabs(z0_bc)-fabs(z1_bc) ); // upsream
  G4double Ld_bc = 0.5*(L0_bc-L1_bc) - ( fabs(z0_bc)-fabs(z1_bc) ); // downstream
  // z coordinates of upstream and downstream elements 
  G4double zu_bc = z1_bc - 0.5*L1_bc - 0.5*Lu_bc;     // upstream
  G4double zd_bc = z1_bc + 0.5*L1_bc + 0.5*Ld_bc;     // downstream

  G4VisAttributes *vis_vac_wire = new G4VisAttributes();
  vis_vac_wire->SetForceWireframe(true); 
  // vis_vac_wire->SetColour( G4Colour::White() ); 
 
  G4Tubs *P0tubeB_upstr_vac             = new G4Tubs("P0tubeB_upstr_vac",0,P0tubeB_rin,Lu_bc/2.,0*deg,360.*deg); 
  G4LogicalVolume *P0tubeB_upstr_vacLog = new G4LogicalVolume(P0tubeB_upstr_vac,GetMaterial("Vacuum"),"P0tubeB_upstr_vac_log",0,0,0);
  P0tubeB_upstr_vacLog->SetVisAttributes(vis_vac_wire); 
 
  G4Tubs *P0tubeB_dnstr_vac             = new G4Tubs("P0tubeB_upstr_vac",0,P0tubeB_rin,Ld_bc/2.,0*deg,360.*deg); 
  G4LogicalVolume *P0tubeB_dnstr_vacLog = new G4LogicalVolume(P0tubeB_dnstr_vac,GetMaterial("Vacuum"),"P0tubeB_dnstr_vac_log",0,0,0);
  P0tubeB_dnstr_vacLog->SetVisAttributes(vis_vac_wire); 

  if(!enableBC_dnstr){
     //Vacuum 0B
     new G4PVPlacement(0,
	   // G4ThreeVector(0.0,0.0,P0initPlacement_z-2.0*P0ringA_L-2.0*P0tubeA_L-P0tubeB_L), 
	   G4ThreeVector(0.0,0.0,P0_vac_b_z), 
	   P0tubeB_vacLog,
	   "P0tubeB_vacLog_pv",
	   worldlog,
	   false,
	   0,
	   ChkOverlaps);
  }else{
     // need to modify things for the placement of the beam collimator
     std::cout << "******** [DOWNSTREAM] BEAM LINE COLLIMATOR INSERTED FOR GEn TARGET PROTECTION [for study only] ********" << std::endl;
     std::cout << "******** Sub-dividing the beam pipe vacuum so the collimator does not overlap: " << std::endl;
     std::cout << "******** upstream vac: z = " << zu_bc/cm << " cm, len = " << Lu_bc/cm << " cm" << std::endl;
     std::cout << "******** Beam collimator length: " << lengthBC_dnstr/cm << " cm" << std::endl; 
     new G4PVPlacement(0,
                       G4ThreeVector(0.0,0.0,zu_bc), 
                       P0tubeB_upstr_vacLog,
                       "P0tubeB_upstr_vacLog_pv",
                       worldlog,
                       false,
                       0,
                       ChkOverlaps);
     std::cout << "******** downstream vac: z = " << zd_bc/cm << " cm, len = " << Ld_bc/cm << " cm" << std::endl;
     new G4PVPlacement(0,
                       G4ThreeVector(0.0,0.0,zd_bc), 
                       P0tubeB_dnstr_vacLog,
                       "P0tubeB_dnstr_vacLog_pv",
                       worldlog,
                       false,
                       0,
                       ChkOverlaps);
  }

  //EPAF: add the Cu radiator *under certain conditions only*
  // i.e the radiator is set to use, but it's distance is above
  if(fDetCon->fTargetBuilder->UseRad()){
    G4double zrad = fDetCon->fTargetBuilder->RadZoffset()+fDetCon->fTargetBuilder->GetTargLen()/2.0;
    //G4cout << "  " << targetEndOffset_z+0.1*cm << " <? " << fDetCon->fTargetBuilder->RadZoffset() << " <? " <<  targetEndOffset_z+P0tubeA_L*2.0 << G4endl;
    G4cout << "  " << P0_vac_a_z+P0tubeA_L << " <? " << -zrad << " <? " << P0_vac_a_z-P0tubeA_L << G4endl;
    
    if(P0_vac_a_z-P0tubeA_L<-zrad && -zrad<P0_vac_a_z+P0tubeA_L){
      cout << zrad << " " << P0_vac_a_z << " => " << -zrad-P0_vac_a_z << " " << P0tubeA_L << endl;
      fDetCon->fTargetBuilder->BuildRadiator( P0tubeA_winvacLog, 0, G4ThreeVector(0., 0., -zrad-P0_vac_a_z) );
    }else if(!enableBC_dnstr){
      //G4cout << "  " << targetEndOffset_z+P0tubeA_L << " <? " << fDetCon->fTargetBuilder->RadZoffset() << " <? " << targetEndOffset_z+P0tubeA_L*2.0+P0tubeB_L*2.0 << G4endl;
      G4cout << "  " << P0_vac_b_z+P0tubeB_L << " <? " << -zrad << " <? " << P0_vac_b_z-P0tubeB_L << G4endl;
      
      if(P0_vac_b_z-P0tubeB_L<-zrad && -zrad<P0_vac_b_z+P0tubeB_L){
	cout << zrad << " " << P0_vac_b_z << " => " << -zrad-P0_vac_b_z << " " << P0tubeB_L << endl;
	fDetCon->fTargetBuilder->BuildRadiator( P0tubeB_vacLog, 0, G4ThreeVector(0., 0.,-zrad-P0_vac_b_z) );
      }
    }
  }
  
  //Tube 0C
  G4double P0tubeC_L = 3.427/2.0*inch; //CJT
  G4double P0tubeC_rin = 1.46/2.0*inch; //CJT
  G4double P0tubeC_rou = 0.846*inch; //CJT

  G4Tubs *P0tubeC = new G4Tubs("P0tubeC", P0tubeC_rin, P0tubeC_rou, P0tubeC_L, 0.0*deg, 360.0*deg);

  G4LogicalVolume *P0tubeCLog = new G4LogicalVolume(P0tubeC, GetMaterial("Aluminum"), "P0tubeC_log", 0, 0, 0);

  //Place tube 0C
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P0initPlacement_z-2.0*P0ringA_L-2.0*P0tubeA_L-2.0*P0tubeB_L-P0tubeC_L), P0tubeCLog, "P0tubeCLog_pv", worldlog, false, 0 , ChkOverlaps );
  
  //Ring 0B
  G4double P0ringB_L = 0.51*inch; //CJT - Doubling length to encompass two flanges here
  G4double P0ringB_rou = 2.74/2.0*inch; //CJT

  G4Tubs *P0ringB = new G4Tubs("P0ringB", P0tubeC_rin, P0ringB_rou, P0ringB_L, 0.0*deg, 360.0*deg);

  G4LogicalVolume *P0ringBLog = new G4LogicalVolume(P0ringB, GetMaterial("Aluminum"), "P0ringB_log", 0, 0, 0);

  //Place ring 0B
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P0initPlacement_z-2.0*P0ringA_L-2.0*P0tubeA_L-2.0*P0tubeB_L-2.0*P0tubeC_L-P0ringB_L), P0ringBLog, "P0ringBLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Vacuum 0C
  G4Tubs *P0tubeC_vac = new G4Tubs("P0tubeC_vac", 0.0, P0tubeC_rin, P0tubeC_L+P0ringB_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P0tubeC_vacLog = new G4LogicalVolume(P0tubeC_vac, GetMaterial("Vacuum"), "P0tubeC_vac_log", 0, 0, 0);

  //Place vacuum 0C
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P0initPlacement_z-2.0*P0ringA_L-2.0*P0tubeA_L-2.0*P0tubeB_L-(P0tubeC_L+P0ringB_L)), P0tubeC_vacLog, "P0tubeC_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube 0D
  //G4double P0tubeD_L = 37.623/2.0*inch; //CJT
  //G4double P0tubeD_L = 117.623/2.0*inch; //Extended beamline for beam studies (4 m extension)
  G4double P0tubeD_L = 433.633/2.0*inch; //Extended beamline for beam studies (10 m extension)

  G4double P0tubeD_rin = 0.685*inch; //CJT
  G4double P0tubeD_rou = 0.75*inch; //CJT

  G4Tubs *P0tubeD = new G4Tubs("P0tubeD", P0tubeD_rin, P0tubeD_rou, P0tubeD_L, 0.0*deg, 360.0*deg);

  G4LogicalVolume *P0tubeDLog = new G4LogicalVolume(P0tubeD, GetMaterial("Aluminum"), "P0tubeD_log", 0, 0, 0);

  //Place tube 0D
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P0initPlacement_z-2.0*P0ringA_L-2.0*P0tubeA_L-2.0*P0tubeB_L-2.0*P0tubeC_L-2.0*P0ringB_L-P0tubeD_L), P0tubeDLog, "P0tubeDLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring 0C
  G4double P0ringC_L = 0.51*inch; //CJT - Doubling length to encompass two flanges here
  G4double P0ringC_rou = 2.74/2.0*inch; //CJT

  G4Tubs *P0ringC = new G4Tubs("P0ringC", P0tubeD_rin, P0ringC_rou, P0ringC_L, 0.0*deg, 360.0*deg);

  G4LogicalVolume *P0ringCLog = new G4LogicalVolume(P0ringC, GetMaterial("Aluminum"), "P0ringC_log", 0, 0, 0);

  //Place ring 0C
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P0initPlacement_z-2.0*P0ringA_L-2.0*P0tubeA_L-2.0*P0tubeB_L-2.0*P0tubeC_L-2.0*P0ringB_L-2.0*P0tubeD_L-P0ringC_L), P0ringCLog, "P0ringCLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Vacuum 0D
  G4Tubs *P0tubeD_vac = new G4Tubs("P0tubeD_vac", 0.0, P0tubeD_rin, P0tubeD_L+P0ringC_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P0tubeD_vacLog = new G4LogicalVolume(P0tubeD_vac, GetMaterial("Vacuum"), "P0tubeD_vac_log", 0, 0, 0);

  //Place vacuum 0D
  // new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P0initPlacement_z-2.0*P0ringA_L-2.0*P0tubeA_L-2.0*P0tubeB_L-2.0*(P0tubeC_L+P0ringB_L)-(P0tubeD_L+P0ringC_L)), P0tubeD_vacLog, "P0tubeD_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  /*
  //SSeeds - Verification that upstream beampipe terminus at -12.0 m from target center
  G4double testRing2_rin = 14.75*inch; 
  G4double testRing2_rou = 15.0*inch; //CJT 13.0 - making larger for debug
  G4double testRing2_L = 0.187/10*inch;
  G4double testPlacement = -472.4*inch; //Beginning S.2
  G4Tubs *testRing2 = new G4Tubs("testRing2", testRing2_rin, testRing2_rou, testRing2_L, 0.*deg, 360.*deg);
  G4LogicalVolume *testRing2Log = new G4LogicalVolume(testRing2, GetMaterial("Air"), "testRing2_log", 0, 0, 0);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, testPlacement), testRing2Log, "testRing2Log_pv", worldlog, false, 0, ChkOverlaps);
  testRing2Log->SetVisAttributes( G4Colour::Red()); //Debug
  */

  // placement of P0 vacuum elements 
  G4double P0_vac_c_z   = P0initPlacement_z-2.0*P0ringA_L-2.0*P0tubeA_L-2.0*P0tubeB_L-(P0tubeC_L+P0ringB_L); 
  G4double P0_vac_d_z   = P0initPlacement_z-2.0*P0ringA_L-2.0*P0tubeA_L-2.0*P0tubeB_L-2.0*(P0tubeC_L+P0ringB_L)-(P0tubeD_L+P0ringC_L);
  // length of P0 vacuum elements (NOTE: Sebastian uses half-lengths! 
  G4double P0_vac_a_len = 2.*P0tubeA_L; 
  G4double P0_vac_b_len = 2.*P0tubeB_L; 
  G4double P0_vac_c_len = 2.*P0tubeC_L;  
  G4double P0_vac_d_len = 2.*(P0tubeD_L+P0ringC_L); // NOTE the difference here! 

  std::vector<G4double> P0_vac_z;
  P0_vac_z.push_back(P0_vac_a_z);   
  P0_vac_z.push_back(P0_vac_b_z);   
  P0_vac_z.push_back(P0_vac_c_z);   
  P0_vac_z.push_back(P0_vac_d_z);   
  std::vector<G4double> P0_vac_len;
  P0_vac_len.push_back(P0_vac_a_len);   
  P0_vac_len.push_back(P0_vac_b_len);   
  P0_vac_len.push_back(P0_vac_c_len);   
  P0_vac_len.push_back(P0_vac_d_len); 
  std::vector<char> P0_label;
  P0_label.push_back('A');  
  P0_label.push_back('B');  
  P0_label.push_back('C');  
  P0_label.push_back('D');  
  const int NP0 = P0_vac_z.size();
 
  char msg[200];  
  // printing info to screen for diagnostics
  if(enableBC_dnstr || enableBC_upstr){ 
     std::cout << "P0 Part Details: " << std::endl;
     for(int i=0;i<NP0;i++){
        sprintf(msg,"tube %c: z = %.3lf cm, len = %.3lf cm",P0_label[i],P0_vac_z[i]/cm,P0_vac_len[i]/cm); 
        std::cout << msg << std::endl;
     }
  }

  // in the event we need to divide up the vacuum 
  G4double lengthBC_upstr = fDetCon->GetBeamCollimatorL_upstr();
  G4double z0_bc_2 = P0_vac_d_z;                      // center of original vacuum insert 
  G4double z1_bc_2 = fDetCon->GetBeamCollimatorZ_upstr();  // center of beam collimator
  G4double L0_bc_2 = P0_vac_d_len; 
  G4double L1_bc_2 = lengthBC_upstr; 
  // lengths of upstream and downstream elements
  G4double Lu_bc_2 = 0.5*(L0_bc_2-L1_bc_2) + ( fabs(z0_bc_2)-fabs(z1_bc_2) ); // upsream
  G4double Ld_bc_2 = 0.5*(L0_bc_2-L1_bc_2) - ( fabs(z0_bc_2)-fabs(z1_bc_2) ); // downstream
  // z coordinates of upstream and downstream elements 
  G4double zu_bc_2 = z1_bc_2 - 0.5*L1_bc_2 - 0.5*Lu_bc_2;     // upstream
  G4double zd_bc_2 = z1_bc_2 + 0.5*L1_bc_2 + 0.5*Ld_bc_2;     // downstream
 
  G4Tubs *P0tubeD_upstr_vac             = new G4Tubs("P0tubeD_upstr_vac",0,P0tubeD_rin,Lu_bc_2/2.,0*deg,360.*deg); 
  G4LogicalVolume *P0tubeD_upstr_vacLog = new G4LogicalVolume(P0tubeD_upstr_vac,GetMaterial("Vacuum"),"P0tubeD_upstr_vac_log",0,0,0);
  P0tubeD_upstr_vacLog->SetVisAttributes(vis_vac_wire); 
 
  G4Tubs *P0tubeD_dnstr_vac             = new G4Tubs("P0tubeD_upstr_vac",0,P0tubeD_rin,Ld_bc_2/2.,0*deg,360.*deg); 
  G4LogicalVolume *P0tubeD_dnstr_vacLog = new G4LogicalVolume(P0tubeD_dnstr_vac,GetMaterial("Vacuum"),"P0tubeD_dnstr_vac_log",0,0,0);
  P0tubeD_dnstr_vacLog->SetVisAttributes(vis_vac_wire); 

  if(!enableBC_upstr){
     //Vacuum 0D
     new G4PVPlacement(0,
	   G4ThreeVector(0.0,0.0,P0_vac_d_z), 
	   P0tubeD_vacLog,
	   "P0tubeD_vacLog_pv",
	   worldlog,
	   false,
	   0,
	   ChkOverlaps);
  }else{
     // need to modify things for the placement of the beam collimator
     std::cout << "******** [UPSTREAM] BEAM LINE COLLIMATOR INSERTED FOR GEn TARGET PROTECTION [for study only] ********" << std::endl;
     std::cout << "******** Sub-dividing the beam pipe vacuum so the collimator does not overlap: " << std::endl;
     std::cout << "******** upstream vac: z = " << zu_bc_2/cm << " cm, len = " << Lu_bc_2/cm << " cm" << std::endl;
     std::cout << "******** Beam collimator length: " << lengthBC_upstr/cm << " cm" << std::endl; 
     new G4PVPlacement(0,
                       G4ThreeVector(0.0,0.0,zu_bc_2), 
                       P0tubeD_upstr_vacLog,
                       "P0tubeD_upstr_vacLog_pv",
                       worldlog,
                       false,
                       0,
                       ChkOverlaps);
     std::cout << "******** downstream vac: z = " << zd_bc_2/cm << " cm, len = " << Ld_bc_2/cm << " cm" << std::endl;
     new G4PVPlacement(0,
                       G4ThreeVector(0.0,0.0,zd_bc_2), 
                       P0tubeD_dnstr_vacLog,
                       "P0tubeD_dnstr_vacLog_pv",
                       worldlog,
                       false,
                       0,
                       ChkOverlaps);
  }      
   
  //Section zero visual attributes
  P0ringA_cutLog->SetVisAttributes( Aluminum);
  P0tubeALog->SetVisAttributes( Aluminum);
  P0tubeBLog->SetVisAttributes( Aluminum);
  P0tubeCLog->SetVisAttributes( Aluminum);
  P0ringBLog->SetVisAttributes( Aluminum);
  P0domeALog->SetVisAttributes( Beryllium);
  P0flangeLog->SetVisAttributes( Aluminum);
  P0tubeDLog->SetVisAttributes( Aluminum);
  P0ringCLog->SetVisAttributes( Aluminum);

  //===== UPSTREAM - PIPE - END =====// UPDATED 10.27.20


  //===== DOWNSTREAM - PIPE - UPSTREAM FLANGES AND WELDMENT =====//
  
  //Section One
  //This section details all components of the exit beamline from the target chamber to the first cone and shielding including all simple cylinders. All labels numbered by proximity to target chamber.

  //G4double targetEndOffset_z = 22.8/2.0*inch;
  //G4double initPlacement_z = (23.62/2.0*inch)+6.267*inch; //CJT -> half target length + target to tube 1A
  G4double initPlacement_z = (targetEndOffset_z)+6.267*inch; //CJT -> half target length + target to tube 1A
  
  //General Specifications
  G4double P1tubeTh = 0.035*inch;
  G4double P1ringL = 0.125/2*inch;

  //Ring A - Most proximal ring to target chamber
  G4double P1ringA_L = 0.11/2*inch;
  G4double P1ringA_rin = 2.93/2*inch;

  G4Tubs *P1ringA = new G4Tubs("P1ringA", P1ringA_rin, P1ringA_rin+P1tubeTh, P1ringA_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringALog = new G4LogicalVolume(P1ringA, GetMaterial("Aluminum"), "P1ringA_log", 0, 0, 0);

  //Place ring A
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, initPlacement_z-P1ringA_L-0.020*inch), P1ringALog, "P1ringALog_pv", worldlog, false, 0 , ChkOverlaps ); //Moving closer to target by half length of ring and thickness of Be window
 
  //Constructing Be "Dome" window upstream from dwg no 67507-0023. sseeds 10.7.20
  G4double P1shell_w = 0.020*inch;
  G4double P1shell_r = 7.136*inch; //Determined from dwg
  G4double P1dome_r = 2.930/2.0*inch;
  G4double P1dome_vd = 0.152*inch; //max depth into beampipe of dome
  G4double P1dome_th = atan(P1dome_r/(P1shell_r-P1dome_vd));
  G4double P1placement; //Make dynamic G4double to keep track of position for each component
  
  G4RotationMatrix *P1dome_rot = new G4RotationMatrix; //Unused
  
  //Construct sphere with theta constraint to produce dome window
  G4Sphere *P1domeA = new G4Sphere("P1domeA", P1shell_r, P1shell_r+P1shell_w, 0.0*deg, 360.0*deg, 0.0*deg, P1dome_th);

  G4LogicalVolume *P1domeALog = new G4LogicalVolume(P1domeA, GetMaterial("Beryllium"), "P1dome_log", 0, 0, 0);
  fDetCon->InsertTargetVolume(P1domeALog->GetName());

  //Place Be window dome
  new G4PVPlacement( P1dome_rot, G4ThreeVector( 0.0, 0.0, initPlacement_z-(P1shell_r+P1shell_w)+P1dome_vd), P1domeALog, "P1domeALog_pv", worldlog, false, 0, ChkOverlaps );
  
  //Tube A - First of three ascending radius tubes. Mates on distant end with ring.
  G4double P1tubeA_L = 1.145/2*inch; //CJT

  P1placement = initPlacement_z;

  G4Tubs *P1tubeA = new G4Tubs("P1tubeA", P1ringA_rin, P1ringA_rin+P1tubeTh, P1tubeA_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeALog = new G4LogicalVolume(P1tubeA, GetMaterial("Aluminum"), "P1tubeA_log", 0, 0, 0);

  //Place tube A
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1tubeA_L), P1tubeALog, "P1tubeALog_pv", worldlog, false, 0 , ChkOverlaps );
 
  //Ring Bin - Small increase in outer radius
  G4double P1ringBin_rou = 3.08/2*inch;
  P1placement += 2.0*P1tubeA_L;

  G4Tubs *P1ringBin = new G4Tubs("P1ringBin", P1ringA_rin, P1ringBin_rou, P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringBinLog = new G4LogicalVolume(P1ringBin, GetMaterial("Aluminum"), "P1ringBin_log", 0, 0, 0);

  //Place ring Bin
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1ringL), P1ringBinLog, "P1ringBinLog_pv", worldlog, false, 0 , ChkOverlaps );
  
  //Ring Bou - Large increase in outer radius to match mating tube
  G4double P1ringBou_rou = 3.75/2*inch;
  P1placement += 2.0*P1ringL;

  G4Tubs *P1ringBou = new G4Tubs("P1ringBou", P1ringA_rin, P1ringBou_rou, P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringBouLog = new G4LogicalVolume(P1ringBou, GetMaterial("Aluminum"), "P1ringBou_log", 0, 0, 0);

  //Place ring Bou
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1ringL), P1ringBouLog, "P1ringBouLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube B - Second of three ascending radius tubes. Mates on both ends with rings.
  G4double P1tubeB_L = 2.440/2*inch; //CJT
  G4double P1tubeB_rin = 3.68/2*inch;
  P1placement += 2.0*P1ringL;

  G4Tubs *P1tubeB = new G4Tubs("P1tubeB", P1tubeB_rin, P1tubeB_rin+P1tubeTh, P1tubeB_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeBLog = new G4LogicalVolume(P1tubeB, GetMaterial("Aluminum"), "P1tubeB_log", 0, 0, 0);

  //Place tube B
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1tubeB_L), P1tubeBLog, "P1tubeBLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring Cin - Small increase in outer radius
  G4double P1ringCin_rou = 3.83/2*inch;
  P1placement += 2.0*P1tubeB_L;

  G4Tubs *P1ringCin = new G4Tubs("P1ringCin", P1tubeB_rin, P1ringCin_rou, P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringCinLog = new G4LogicalVolume(P1ringCin, GetMaterial("Aluminum"), "P1ringCin_log", 0, 0, 0);

  //Place ring Cin
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1ringL), P1ringCinLog, "P1ringCinLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring Cou - Large increase in outer radius to match mating tube
  G4double P1ringCou_rou = 4.5/2*inch;
  P1placement += 2.0*P1ringL;

  G4Tubs *P1ringCou = new G4Tubs("P1ringCou", P1tubeB_rin, P1ringCou_rou, P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringCouLog = new G4LogicalVolume(P1ringCou, GetMaterial("Aluminum"), "P1ringCou_log", 0, 0, 0);

  //Place ring Cou
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1ringL), P1ringCouLog, "P1ringCouLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube C - Third of three ascending radius tubes. Mates on both ends with rings.
  G4double P1tubeC_L = 2.25/2*inch; //CJT
  G4double P1tubeC_rin = 4.43/2*inch;
  P1placement += 2.0*P1ringL;

  G4Tubs *P1tubeC = new G4Tubs("P1tubeC", P1tubeC_rin, P1tubeC_rin+P1tubeTh, P1tubeC_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeCLog = new G4LogicalVolume(P1tubeC, GetMaterial("Aluminum"), "P1tubeC_log", 0, 0, 0);

  //Place tube C
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1tubeC_L), P1tubeCLog, "P1tubeCLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring D in and out twice the length of previous rings - remaining elements deviate from general specifications above
  //Ring Din
  G4double P1ringDin_rou = 4.57/2*inch;
  P1placement += 2.0*P1tubeC_L;

  G4Tubs *P1ringDin = new G4Tubs("P1ringDin", P1tubeC_rin, P1ringDin_rou, 2*P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringDinLog = new G4LogicalVolume(P1ringDin, GetMaterial("Aluminum"), "P1ringDin_log", 0, 0, 0);

  //Place ring Din
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+2.0*P1ringL), P1ringDinLog, "P1ringDinLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring Dou
  G4double P1ringDou_rou = 6.0/2*inch;
  P1placement += 4.0*P1ringL;

  G4Tubs *P1ringDou = new G4Tubs("P1ringCou", P1tubeC_rin, P1ringDou_rou, 2*P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringDouLog = new G4LogicalVolume(P1ringDou, GetMaterial("Aluminum"), "P1ringDou_log", 0, 0, 0);

  //Place ring Dou
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+2.0*P1ringL), P1ringDouLog, "P1ringDouLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube D
  G4double P1tubeD_L = 19.344/2*inch; 
  G4double P1tubeD_rin = 2.875*inch; //CJT
  G4double P1tubeD_rou = 6.0/2*inch;
  P1placement += 4.0*P1ringL;

  G4Tubs *P1tubeD = new G4Tubs("P1tubeD", P1tubeD_rin, P1tubeD_rou, P1tubeD_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeDLog = new G4LogicalVolume(P1tubeD, GetMaterial("Aluminum"), "P1tubeD_log", 0, 0, 0);

  //Place tube D
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1tubeD_L), P1tubeDLog, "P1tubeDLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring E  //CJT combines ring E and tube E by color, but each differs in Rin
  G4double P1ringE_L = 0.5/2*inch;
  G4double P1ringE_rou = 3.150*inch; //CJT
  P1placement += 2.0*P1tubeD_L;

  G4Tubs *P1ringE = new G4Tubs("P1ringE", P1tubeD_rin, P1ringE_rou, P1ringE_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringELog = new G4LogicalVolume(P1ringE, GetMaterial("Aluminum"), "P1ringE_log", 0, 0, 0);

  //Place ring E
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1ringE_L), P1ringELog, "P1ringELog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube E
  G4double P1tubeE_L = 2.571/2*inch;
  G4double P1tubeE_rin = 2.953*inch; //CJT
  P1placement += 2.0*P1ringE_L;

  G4Tubs *P1tubeE = new G4Tubs("P1tubeD", P1tubeE_rin, P1ringE_rou, P1tubeE_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeELog = new G4LogicalVolume(P1tubeE, GetMaterial("Aluminum"), "P1tubeE_log", 0, 0, 0);

  //Place tube E
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1tubeE_L), P1tubeELog, "P1tubeELog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring F
  G4double P1ringF_L = 0.866/2*inch;
  G4double P1ringF_rou = 3.996*inch; //CJT
  P1placement += 2.0*P1tubeE_L;

  G4Tubs *P1ringF = new G4Tubs("P1ringF", P1tubeE_rin, P1ringF_rou, P1ringF_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringFLog = new G4LogicalVolume(P1ringF, GetMaterial("Aluminum"), "P1ringF_log", 0, 0, 0);

  //Place ring F
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1ringF_L), P1ringFLog, "P1ringFLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring G
  G4double P1ringG_L = 0.52/2*inch; //CJT
  G4double P1ringG_rin = 2.312*inch; //Slight modification in geometry here to simplify. Actual 2.46*inch by CJT.
  P1placement += 2.0*P1ringF_L;

  G4Tubs *P1ringG = new G4Tubs("P1ringG", P1ringG_rin, P1ringF_rou, P1ringG_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringGLog = new G4LogicalVolume(P1ringG, GetMaterial("Aluminum"), "P1ringG_log", 0, 0, 0);

  //Place ring G
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1ringG_L), P1ringGLog, "P1ringGLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube F
  G4double P1tubeF_L = 13.89/2*inch; //CJT
  G4double P1tubeF_rin = 2.312*inch; //CJT
  G4double P1tubeF_rou = 2.5*inch; //CJT
  P1placement += 2.0*P1ringG_L;

  G4Tubs *P1tubeF = new G4Tubs("P1tubeF", P1tubeF_rin, P1tubeF_rou, P1tubeF_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeFLog = new G4LogicalVolume(P1tubeF, GetMaterial("Aluminum"), "P1tubeF_log", 0, 0, 0);

  //Place tube F
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1tubeF_L), P1tubeFLog, "P1tubeFLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring H
  G4double P1ringH_L = 0.42/2*inch;
  G4double P1ringH_rou = 3.36*inch; //CJT
  P1placement += 2.0*P1tubeF_L;

  G4Tubs *P1ringH = new G4Tubs("P1ringH", P1tubeF_rin, P1ringH_rou, P1ringH_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringHLog = new G4LogicalVolume(P1ringH, GetMaterial("Aluminum"), "P1ringH_log", 0, 0, 0);

  //Place ring H
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1ringH_L), P1ringHLog, "P1ringHLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring I
  G4double P1ringI_rin = 3.92/2*inch;
  P1placement += 2.0*P1ringH_L;

  G4Tubs *P1ringI = new G4Tubs("P1ringI", P1ringI_rin, P1ringH_rou, P1ringH_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringILog = new G4LogicalVolume(P1ringI, GetMaterial("Aluminum"), "P1ringI_log", 0, 0, 0);

  //Place ring I
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1ringH_L), P1ringILog, "P1ringILog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring J
  //This ring mates with section two cone at 3.571 inches diameter at entrance of cone and difference is ignored.
  G4double P1ringJ_L = 0.695/2*inch; //CJT
  G4double P1ringJ_rin = 3.5/2*inch;
  P1placement += 2.0*P1ringH_L;

  G4Tubs *P1ringJ = new G4Tubs("P1ringJ", P1ringJ_rin, P1ringH_rou, P1ringJ_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringJLog = new G4LogicalVolume(P1ringJ, GetMaterial("Aluminum"), "P1ringJ_log", 0, 0, 0);

  //Place ring J
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1ringJ_L), P1ringJLog, "P1ringJLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Vacuum fill solids - adding an object for each new rin only

  //Tube A vacuum
  P1placement = initPlacement_z; //Resetting placement double
  G4Tubs *P1tubeA_vac = new G4Tubs("P1tubeA_vac", 0.0, P1ringA_rin, P1tubeA_L+2*P1ringL, 0.*deg, 360.*deg);

  //Subract the colume walled off by the Be dome
  G4Sphere *P1sphere_cut = new G4Sphere("P0sphere_cut", 0.0, P1shell_r+P1shell_w, 0.0*deg, 360.0*deg, 0.0*deg, 360.0*deg);

  G4SubtractionSolid *P1tubeA_winvac = new G4SubtractionSolid("P1tubeA_winvac", P1tubeA_vac, P1sphere_cut, 0, G4ThreeVector(0.0, 0.0, -(P1tubeA_L+2*P1ringL)-(P1shell_r+P1shell_w)+P1dome_vd));

  G4LogicalVolume *P1tubeA_winvacLog = new G4LogicalVolume(P1tubeA_winvac, GetMaterial("Vacuum"), "P1tubeA_winvac_log", 0, 0, 0);

  //Place tube A vacuum
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1tubeA_L+2.0*P1ringL+P1shell_w), P1tubeA_winvacLog, "P1tubeA_winvacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube B vacuum
  P1placement += 2.0*P1tubeA_L+4.0*P1ringL;
  G4Tubs *P1tubeB_vac = new G4Tubs("P1tubeB_vac", 0.0, P1tubeB_rin, P1tubeB_L+2.0*P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeB_vacLog = new G4LogicalVolume(P1tubeB_vac, GetMaterial("Vacuum"), "P1tubeB_vac_log", 0, 0, 0);

  //Place tube B vacuum
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1tubeB_L+2.0*P1ringL), P1tubeB_vacLog, "P1tubeB_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube C vacuum
  P1placement += 2.0*P1tubeB_L+4.0*P1ringL;
  G4Tubs *P1tubeC_vac = new G4Tubs("P1tubeC_vac", 0.0, P1tubeC_rin, P1tubeC_L+4.0*P1ringL, 0.*deg, 360.*deg);

  G4LogicalVolume *P1tubeC_vacLog = new G4LogicalVolume(P1tubeC_vac, GetMaterial("Vacuum"), "P1tubeC_vac_log", 0, 0, 0);

  //Place tube C vacuum
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1tubeC_L+4.0*P1ringL), P1tubeC_vacLog, "P1tubeC_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube D vacuum
  P1placement += 2.0*P1tubeC_L+8.0*P1ringL;
  G4Tubs *P1tubeD_vac = new G4Tubs("P1tubeD_vac", 0.0, P1tubeD_rin, P1tubeD_L+P1ringE_L, 0.*deg, 360.*deg);
  
  G4LogicalVolume *P1tubeD_vacLog = new G4LogicalVolume(P1tubeD_vac, GetMaterial("Vacuum"), "P1tubeD_vac_log", 0, 0, 0);

  //Place tube D vacuum
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1tubeD_L+P1ringE_L), P1tubeD_vacLog, "P1tubeD_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube E vacuum
  P1placement += 2.0*P1tubeD_L+2.0*P1ringE_L;
  G4Tubs *P1tubeE_vac = new G4Tubs("P1tubeE_vac", 0.0, P1tubeE_rin, P1tubeE_L+P1ringF_L, 0.*deg, 360.*deg);
  
  G4LogicalVolume *P1tubeE_vacLog = new G4LogicalVolume(P1tubeE_vac, GetMaterial("Vacuum"), "P1tubeE_vac_log", 0, 0, 0);

  //Place tube E vacuum
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1tubeE_L+P1ringF_L), P1tubeE_vacLog, "P1tubeE_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Tube F Vacuum
  P1placement += 2.0*P1tubeE_L+2.0*P1ringF_L;
  G4Tubs *P1tubeF_vac = new G4Tubs("P1tubeF_vac", 0.0, P1ringG_rin, P1ringG_L+P1tubeF_L+P1ringH_L, 0.*deg, 360.*deg);
  
  G4LogicalVolume *P1tubeF_vacLog = new G4LogicalVolume(P1tubeF_vac, GetMaterial("Vacuum"), "P1tubeF_vac_log", 0, 0, 0);

  //Place tube F vacuum
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1ringG_L+P1tubeF_L+P1ringH_L), P1tubeF_vacLog, "P1tubeF_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring I Vacuum
  P1placement += 2.0*P1ringG_L+2.0*P1tubeF_L+2.0*P1ringH_L;
  G4Tubs *P1ringI_vac = new G4Tubs("P1ringI_vac", 0.0, P1ringI_rin, P1ringH_L, 0.*deg, 360.*deg);
  
  G4LogicalVolume *P1ringI_vacLog = new G4LogicalVolume(P1ringI_vac, GetMaterial("Vacuum"), "P1ringI_vac_log", 0, 0, 0);

  //Place ring I vacuum
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1ringH_L), P1ringI_vacLog, "P1ringI_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring J Vacuum
  P1placement += 2.0*P1ringH_L;
  G4Tubs *P1ringJ_vac = new G4Tubs("P1ringJ_vac", 0.0, P1ringG_rin, P1ringJ_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P1ringJ_vacLog = new G4LogicalVolume(P1ringJ_vac, GetMaterial("Vacuum"), "P1ringJ_vac_log", 0, 0, 0);

  //Place ring J vacuum
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P1placement+P1ringJ_L), P1ringJ_vacLog, "P1ringJ_vacLog_pv", worldlog, false, 0 , ChkOverlaps );
  
  //Visuals
  G4VisAttributes * WireFrameVisAtt = new G4VisAttributes(G4Colour(0.,1.,1.));  //Keeping for debugging purposes
  WireFrameVisAtt->SetForceWireframe(true);
  //P1tubeBLog->SetVisAttributes( WireFrameVisAtt );
  //P1tubeA_winvacLog->SetVisAttributes( G4VisAttributes::GetInvisible());

  P1ringALog->SetVisAttributes( Aluminum);
  P1domeALog->SetVisAttributes( Beryllium);
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


  //===== DOWNSTREAM - PIPE - UPSTREAM FLANGES AND WELDMENT - END =====//

  //===== DOWNSTREAM - PIPE - CONES =====//
  //===== DOWNSTREAM - MAG SHIELDING ====//
  
  //SECTION TWO 
  //This section includes all components up to and including the first corrector magnet.
  //Specifications - ordered first from inside to out then from target-proximal to distant

  //===== DOWNSTREAM - PIPE - US MIDDLE CONE =====//
  
  //General
  //G4double targetEndOffset_z = 22.8/2.0*inch;
  G4double P2initPlacement_z = (targetEndOffset_z)+51.934*inch; //CJT -> Half target length + target to tube 1A + length of first weldment and flanges 
  G4double P2_offset = 0.451*inch; //CJT -> Gap between start of inner cone and outer-cone/mag-shielding
  G4VisAttributes *DebugGreen = new G4VisAttributes(G4Colour(0.0,1.0,0.0)); //Debug option left in
  
  //Inner cone - runs the length of the section, through both corrector magnets to the endcap section of weldments and flanges before exit beamline
  G4double P2coneA_rin1 = 3.517/2*inch;
  G4double P2coneA_rou1 = 3.767/2*inch;
  G4double P2coneA_rin2 = 10.734/2*inch;
  G4double P2coneA_rou2 = 10.984/2*inch;
  G4double P2coneA_L = 137.8/2*inch; //CJT

  G4Cons *P2coneA = new G4Cons("P2coneA", P2coneA_rin1, P2coneA_rou1, P2coneA_rin2, P2coneA_rou2, P2coneA_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P2coneALog = new G4LogicalVolume(P2coneA, GetMaterial("Aluminum"), "P2coneA_log", 0, 0, 0);

  //Place the inner cone
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2coneA_L), P2coneALog, "P2coneALog_pv", worldlog, false, 0 , ChkOverlaps);
  
  P2coneALog->SetVisAttributes( Aluminum);

  //Inner Cone Vacuum
  G4Cons *P2coneA_vac = new G4Cons("P2coneA_vac", 0.0, P2coneA_rin1, 0.0, P2coneA_rin2, P2coneA_L, 0.*deg, 360.*deg);
  
  G4LogicalVolume *P2coneA_vacLog = new G4LogicalVolume(P2coneA_vac, GetMaterial("Vacuum"), "P2coneA_vac_log", 0, 0, 0);

  //Place inner cone vacuum
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2coneA_L), P2coneA_vacLog, "P2coneA_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Middle cone -> staged, wrapped around inner cone; section from US start of cone to first corrector magnet
  G4double P2coneB_rin1 = 3.517/2*inch;
  G4double P2coneB_rou1 = 3.767/2*inch;
  G4double P2coneB_rin2 = 5.278/2*inch; //CJT 
  G4double P2coneB_rou2 = 5.528/2*inch; //CJT
  G4double P2coneB_L = 33.625/2*inch; //CJT
  
  G4Cons *P2coneB = new G4Cons("P2coneB", P2coneB_rin1, P2coneB_rou1, P2coneB_rin2, P2coneB_rou2, P2coneB_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P2coneBLog = new G4LogicalVolume(P2coneB, GetMaterial("Aluminum"), "P2coneB_log", 0, 0, 0);

  //Place middle cone
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2initPlacement_z+P2_offset+P2coneB_L), P2coneBLog, "P2coneBLog_pv", worldlog, false, 0 , ChkOverlaps);

  P2coneBLog->SetVisAttributes( Aluminum);

  //===== DOWNSTREAM - PIPE - US MIDDLE CONE - END =====//

  //===== MAG SHIELD, RINGS - US SECTION =====//
  
  //Rings General Specifications, 18 in total
  //G4double P2_offset = 0.451*inch; //DUPE FOR BREAKOUT
  G4double P2ringTh = 0.5*inch;
  G4double P2ringr0 = 1.895*inch;
  G4double P2ringL = 1.625/2*inch;
  G4double P2DAngle = 1.5*deg;
  G4double P2ringSep = 0.375*inch;
  G4double P2placement; //Dynamic G4double to keep track of position
  G4double P2ring_rin1; //US radius of conic section
  G4double P2ring_rin2; //DS radius of conic section


  for (int i=0; i<18; i++){
    P2ring_rin1 = P2ringr0+(2*i*P2ringL+i*P2ringSep)*tan(P2DAngle);
    P2ring_rin2 = P2ringr0+(2*(i+1)*P2ringL+i*P2ringSep)*tan(P2DAngle);
    P2placement = P2initPlacement_z+P2_offset+2*i*P2ringL+i*P2ringSep+P2ringL;

    G4Cons *P2ring = new G4Cons(Form("P2ring%d",i+1), P2ring_rin1, P2ring_rin1+P2ringTh, P2ring_rin2, P2ring_rin2+P2ringTh, P2ringL, 0.*deg, 360.*deg);
    
    G4LogicalVolume *P2ringLog = new G4LogicalVolume(P2ring, GetMaterial("Iron"), Form("P2ring%d_log",i+1), 0, 0, 0);

    new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P2placement), P2ringLog, Form("P2ring%dLog_pv",i+1), worldlog, false, 0 , ChkOverlaps);

    P2ringLog->SetVisAttributes(Iron);
  }

  //===== MAG SHIELD, RINGS - US SECTION - END =====//
  
  //===== MOUNTING PLATES - US SECTION =====//
  
  //Mounting plates - Aluminum plates NOT radiation shields per material update 7.24.20
  //Specifications
  G4double sideshieldTh = 0.25*inch;
  G4double sideshieldA = 1.06*deg; //CJT
  G4double P2sideshieldL = 37.125/2*inch; //CJT
  G4double P2sideshieldW1 = 4.762/2*inch;
  G4double P2sideshieldW2 = 6.137/2*inch;
  G4double P2sideshield_xoffset = 8.282/2*inch+sideshieldTh;

  //Symmetric mounting plates - construct using trapezoid class 
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

  //Place both mounting plates
  new G4PVPlacement( P2rot1_temp, G4ThreeVector(-P2sideshield_xoffset-P2sideshieldL*sin(sideshieldA), 0, P2initPlacement_z+P2_offset+P2sideshieldL*cos(sideshieldA)), P2sideshield1_log, "P2sideshield1_log", worldlog, false, 0, ChkOverlaps);
  P2sideshield1_log->SetVisAttributes(Aluminum);
  new G4PVPlacement( P2rot2_temp, G4ThreeVector(P2sideshield_xoffset+P2sideshieldL*sin(sideshieldA), 0, P2initPlacement_z+P2_offset+P2sideshieldL*cos(sideshieldA)), P2sideshield2_log, "P2sideshield2_log", worldlog, false, 0, ChkOverlaps);
  P2sideshield2_log->SetVisAttributes(Aluminum);

  //===== MOUNTING PLATES - US SECTION - END =====//

  //SECTION THREE
  //This section includes all components from the first corrector magnet up to and including the second corrector magnet.

  //===== DOWNSTREAM - PIPE - MIDDLE MIDDLE CONE =====//
  
  //General specifications - ordered first from inside to out then from target-proximal to distant
  //G4double targetEndOffset_z = 22.8/2.0*inch;
  G4double P3initPlacement_z = (targetEndOffset_z)+102.385*inch; //Direct CJT measure

  //Inner cone included in section one code

  G4double P3ringL = 1.625/2*inch;
  G4double P3ringSep = 0.375*inch;
  
  G4double P3_offset = 2*P3ringL+P3ringSep;

  //Middle cone - staged, wrapping the inner cone in the central section
  G4double P3coneB_rin1 = 3.257*inch;
  G4double P3coneB_rou1 = 3.507*inch;
  G4double P3coneB_rin2 = 4.556*inch;
  G4double P3coneB_rou2 = 4.807*inch;
  G4double P3coneB_L = 49.625/2*inch;
  
  G4Cons *P3coneB = new G4Cons("P3coneB", P3coneB_rin1, P3coneB_rou1, P3coneB_rin2, P3coneB_rou2, P3coneB_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P3coneBLog = new G4LogicalVolume(P3coneB, GetMaterial("Aluminum"), "P3coneB_log", 0, 0, 0);

  //Place middle cone
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3initPlacement_z+P3_offset+P3coneB_L), P3coneBLog, "P3coneBLog_pv", worldlog, false, 0 , ChkOverlaps);

  //P3coneBLog->SetVisAttributes( DebugRed );  //Debug, left in
  P3coneBLog->SetVisAttributes(Aluminum);

  //===== DOWNSTREAM - PIPE - MIDDLE MIDDLE CONE - END =====//

  //===== MAG SHIELD, RINGS - MIDDLE SECTION =====//
  
  //Rings General Specifications, 27 in total
  //G4double P3initPlacement_z = (23.62/2.0*inch)+102.385*inch; //Direct CJT measure //DUPE FOR BREAKOUT
  
  G4double P3ringTh = 0.5*inch;
  G4double P2P3displacement = 50.451*inch; //CJT distance between start P2 and start P3 elements
  G4double P3ringr0 = P2ringr0+(P2P3displacement-P2_offset)*tan(P2DAngle);

  //G4double P3ringL = 1.625/2*inch; //DUPE FOR BREAKOUT
  //G4double P3ringSep = 0.375*inch; //DUPE FOR BREAKOUT

  G4double P3ring_rin1;
  G4double P3ring_rin2;
  G4double P3placement;
  
  for (int i=0; i<27; i++){
    P3ring_rin1 = P3ringr0+(2*i*P3ringL+i*P3ringSep)*tan(P2DAngle);
    P3ring_rin2 = P3ringr0+(2*(i+1)*P3ringL+i*P3ringSep)*tan(P2DAngle);
    P3placement = P3initPlacement_z+2*i*P3ringL+i*P3ringSep+P3ringL;

    G4Cons *P3ring = new G4Cons(Form("P3ring%d",i+1), P3ring_rin1, P3ring_rin1+P3ringTh, P3ring_rin2, P3ring_rin2+P3ringTh, P3ringL, 0.*deg, 360.*deg);
    
    G4LogicalVolume *P3ringLog = new G4LogicalVolume(P3ring, GetMaterial("Iron"), Form("P3ring%d_log",i+1), 0, 0, 0);

    new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P3placement), P3ringLog, Form("P3ring%dLog_pv",i+1), worldlog, false, 0 , ChkOverlaps);

    P3ringLog->SetVisAttributes(Iron);
  }

  //===== MAG SHIELD, RINGS - MIDDLE SECTION - END =====//

  //===== MOUNTING PLATES - MIDDLE SECTION =====//
  
  //P3 Side shields - NOT shields per material update 7.24.20
  //Specifications
  G4double P3sideshieldL = 57.750/2*inch; //CJT
  G4double P3sideshieldW1 = 6.547/2*inch;
  G4double P3sideshieldW2 = 8.688/2*inch;
  G4double P3sideshield_xoffset = 10.067/2*inch+sideshieldTh;

  //Symmetric mounting plates - construct using trapezoid class 
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

  //Place both mounting plates
  new G4PVPlacement( P3rot1_temp, G4ThreeVector(-P3sideshield_xoffset-P3sideshieldL*sin(sideshieldA), 0, P3initPlacement_z-P3_offset+P3sideshieldL*cos(sideshieldA)), P3sideshield1_log, "P3sideshield1_log", worldlog, false, 0, ChkOverlaps);
  P3sideshield1_log->SetVisAttributes(Aluminum);
  new G4PVPlacement( P3rot2_temp, G4ThreeVector(P3sideshield_xoffset+P3sideshieldL*sin(sideshieldA), 0, P3initPlacement_z-P3_offset+P3sideshieldL*cos(sideshieldA)), P3sideshield2_log, "P3sideshield2_log", worldlog, false, 0, ChkOverlaps);
  P3sideshield2_log->SetVisAttributes(Aluminum);

  //===== MOUNTING PLATES - MIDDLE SECTION - END =====//
  
  //===== DOWNSTREAM - PIPE - DS MIDDLE CONE =====//

  //SECTION FOUR
  //This section includes all components from the second corrector magnet.
  //Specifications - ordered first from inside to out then from target-proximal to distant

  //G4double P2P4displacement = 122.141*inch;
  G4double P2P4displacement = 122.451*inch;  //CJT
  //G4double targetEndOffset_z = 22.8/2.0*inch;
  G4double P4initPlacement_z = (targetEndOffset_z)+174.385*inch; //Direct CJT measure
  
  //Inner cone included in section one code

  G4double P4ringL = 1.625/2*inch;
  G4double P4ringSep = 0.375*inch;
  G4double P4_offset = 2*P4ringL+P4ringSep;

  //DS middle cone - staged, wrapping the inner cone in the central section
  G4double P4coneB_rin1 = 5.142*inch;
  G4double P4coneB_rou1 = 5.392*inch;
  G4double P4coneB_rin2 = 5.446*inch;
  G4double P4coneB_rou2 = 5.696*inch;
  G4double P4coneB_L = 11.591/2*inch;
  
  G4Cons *P4coneB = new G4Cons("P4coneB", P4coneB_rin1, P4coneB_rou1, P4coneB_rin2, P4coneB_rou2, P4coneB_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P4coneBLog = new G4LogicalVolume(P4coneB, GetMaterial("Aluminum"), "P4coneB_log", 0, 0, 0);

  //Place DS middle cone
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P4initPlacement_z+P4_offset+P4coneB_L), P4coneBLog, "P4coneBLog_pv", worldlog, false, 0 , ChkOverlaps);

  //P4coneBLog->SetVisAttributes( DebugRed ); //Debug, keep in
  P4coneBLog->SetVisAttributes(Aluminum);

  //===== DOWNSTREAM - PIPE - DS MIDDLE CONE - END =====//

  //===== MAG SHIELD, RINGS - DS SECTION - END =====//
  
  //Rings General Specifications
  //G4double P2P4displacement = 122.451*inch;  //CJT //DUPE FOR BREAKOUT
  //G4double P4initPlacement_z = (23.62/2.0*inch)+174.385*inch; //Direct CJT measure //DUPE FOR BREAKOUT
  G4double P4ringTh = 0.5*inch;
  G4double P4ringr0 = P2ringr0+(P2P4displacement-P2_offset)*tan(P2DAngle);
  //G4double P4ringL = 1.625/2*inch; //DUPE FOR BREAKOUT
  //G4double P4ringSep = 0.375*inch; //DUPE FOR BREAKOUT
  G4double P4ring_rin1;
  G4double P4ring_rin2;
  G4double P4placement;

  for (int i=0; i<7; i++){
    P4ring_rin1 = P4ringr0+(2*i*P4ringL+i*P4ringSep)*tan(P2DAngle);
    P4ring_rin2 = P4ringr0+(2*(i+1)*P4ringL+i*P4ringSep)*tan(P2DAngle);
    P4placement = P4initPlacement_z+2*i*P4ringL+i*P4ringSep+P4ringL;

    G4Cons *P4ring = new G4Cons(Form("P4ring%d",i+1), P4ring_rin1, P4ring_rin1+P4ringTh, P4ring_rin2, P4ring_rin2+P4ringTh, P4ringL, 0.*deg, 360.*deg);
    
    G4LogicalVolume *P4ringLog = new G4LogicalVolume(P4ring, GetMaterial("Iron"), Form("P4ring%d_log",i+1), 0, 0, 0);

    new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P4placement), P4ringLog, Form("P4ring%dLog_pv",i+1), worldlog, false, 0 , ChkOverlaps);

    P4ringLog->SetVisAttributes(Iron);
  }

  //===== MAG SHIELD, RINGS - DS SECTION - END =====//

  //===== MOUNTING PLATES - DS SECTION =====//
  
  //P4 Side shields - NOT shields per material update 7.24.20
  //Specifications
  G4double P4sideshieldL = 13.625/2*inch;  //CJT
  G4double P4sideshieldW1 = 9.280/2*inch;
  G4double P4sideshieldW2 = 9.780/2*inch;
  G4double P4sideshield_xoffset = 12.8/2*inch+sideshieldTh;
 
  //Symmetric mounting plates - construct using trapezoid class 
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

  //Place both mounting plates
  new G4PVPlacement( P4rot1_temp, G4ThreeVector(-P4sideshield_xoffset-P4sideshieldL*sin(sideshieldA), 0, P4initPlacement_z+P4sideshieldL*cos(sideshieldA)), P4sideshield1_log, "P4sideshield1_log", worldlog, false, 0, ChkOverlaps);
  P4sideshield1_log->SetVisAttributes(Aluminum);
  new G4PVPlacement( P4rot2_temp, G4ThreeVector(P4sideshield_xoffset+P4sideshieldL*sin(sideshieldA), 0, P4initPlacement_z+P4sideshieldL*cos(sideshieldA)), P4sideshield2_log, "P4sideshield2_log", worldlog, false, 0, ChkOverlaps);
  P4sideshield2_log->SetVisAttributes(Aluminum);

  //===== MOUNTING PLATES - DS SECTION - END =====//
  
  //P5 Endcap
  
  //General Specifications
  //G4double targetEndOffset_z = 22.8/2.0*inch;
  G4double P5initPlacement_z = (targetEndOffset_z)+189.735*inch; //Direct CJT measure
  G4double P5placement = P5initPlacement_z;
  
  //Ring A
  G4double P5ringA_rin = 10.734/2*inch; //Same as inner cone rin at wider end
  G4double P5ringA_rou = 13.940/2*inch;
  G4double P5ringA_L = 0.52/2*inch; //CJT

  G4Tubs *P5ringA = new G4Tubs("P5ringA", P5ringA_rin, P5ringA_rou, P5ringA_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringALog = new G4LogicalVolume(P5ringA, GetMaterial("Aluminum"), "P5ringA_log", 0, 0, 0);

  //Place ring A
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5placement+P5ringA_L), P5ringALog, "P5ringALog_pv", worldlog, false, 0, ChkOverlaps);
  P5ringALog->SetVisAttributes(Aluminum);

  //Ring A Vacuum
  G4Tubs *P5ringA_vac = new G4Tubs("P5ringA_vac", 0.0, P5ringA_rin, P5ringA_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringA_vacLog = new G4LogicalVolume(P5ringA_vac, GetMaterial("Vacuum"), "P5ringA_vac_log", 0, 0, 0);

  //Place ring A vacuum
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P5placement+P5ringA_L), P5ringA_vacLog, "P5ringA_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring B
  G4double P5ringB_rin = 5.88*inch; //CJT
  G4double P5ringB_rou = 14.0/2*inch; //CJT
  G4double P5ringB_L = 1.120/2*inch;
  P5placement += 2.0*P5ringA_L;

  G4Tubs *P5ringB = new G4Tubs("P5ringB", P5ringB_rin, P5ringB_rou, P5ringB_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringBLog = new G4LogicalVolume(P5ringB, GetMaterial("Aluminum"), "P5ringB_log", 0, 0, 0);

  //Place ring B
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5placement+P5ringB_L), P5ringBLog, "P5ringBLog_pv", worldlog, false, 0, ChkOverlaps);
  P5ringBLog->SetVisAttributes(Aluminum);

  //Ring B Vacuum
  G4Tubs *P5ringB_vac = new G4Tubs("P5ringB_vac", 0.0, P5ringB_rin, P5ringB_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringB_vacLog = new G4LogicalVolume(P5ringB_vac, GetMaterial("Vacuum"), "P5ringB_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P5placement+P5ringB_L), P5ringB_vacLog, "P5ringB_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //Ring C - Simplifying geometry of this odd weldment
  //G4double P5ringC_rin = 11.750/2*inch; 
  G4double P5ringC_rin = 5.88*inch; //CJT
  //G4double P5ringC_rou = 12.0/2*inch; //CJT
  G4double P5ringC_rou = 7.0*inch; //CJT
  //G4double P5ringC_L = 0.405*inch;
  //G4double P5ringC_L = 0.75/2*inch; //CJT
  G4double P5ringC_L = 3.267/2*inch; //CJT
  P5placement += 2.0*P5ringB_L;

  G4Tubs *P5ringC = new G4Tubs("P5ringC", P5ringC_rin, P5ringC_rou, P5ringC_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringCLog = new G4LogicalVolume(P5ringC, GetMaterial("Aluminum"), "P5ringC_log", 0, 0, 0);

  //Place ring C
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5placement+P5ringC_L), P5ringCLog, "P5ringCLog_pv", worldlog, false, 0, ChkOverlaps);
  P5ringCLog->SetVisAttributes(Aluminum);

  //Ring C Vacuum
  G4Tubs *P5ringC_vac = new G4Tubs("P5ringC_vac", 0.0, P5ringC_rin, P5ringC_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringC_vacLog = new G4LogicalVolume(P5ringC_vac, GetMaterial("Vacuum"), "P5ringC_vac_log", 0, 0, 0);

  //Place ring C vacuum
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P5placement+P5ringC_L), P5ringC_vacLog, "P5ringC_vacLog_pv", worldlog, false, 0 , ChkOverlaps );
  
  //Ring D - Combining two flanges, simplifying rin (slightly smaller 5.875" CJT actual for second half)
  //G4double P5ringD_rin = 12.710/2*inch; 
  G4double P5ringD_rin = 5.88*inch; //CJT
  //G4double P5ringD_rou = 12.750/2*inch;
  G4double P5ringD_rou = 8.25*inch; //CJT
  //G4double P5ringD_L = 2.487*inch;
  //G4double P5ringD_L = 2.260/2*inch; //CJT
  G4double P5ringD_L = 1.120/2*inch; //CJT - ignoring the adjacent flange as it is included in target to midpipe geometry within exit beamline 
  P5placement += 2.0*P5ringC_L;

  G4Tubs *P5ringD = new G4Tubs("P5ringD", P5ringD_rin, P5ringD_rou, P5ringD_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringDLog = new G4LogicalVolume(P5ringD, GetMaterial("Aluminum"), "P5ringD_log", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5placement+P5ringD_L), P5ringDLog, "P5ringDLog_pv", worldlog, false, 0, ChkOverlaps);
  P5ringDLog->SetVisAttributes(Aluminum);
  //P5ringDLog->SetVisAttributes(G4Colour::Green());
  
  //Ring D Vacuum
  G4Tubs *P5ringD_vac = new G4Tubs("P5ringD_vac", 0.0, P5ringD_rin, P5ringD_L+0.015*inch, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringD_vacLog = new G4LogicalVolume(P5ringD_vac, GetMaterial("Vacuum"), "P5ringD_vac_log", 0, 0, 0);

  //Place ring D vacuum
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P5placement+P5ringD_L), P5ringD_vacLog, "P5ringD_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  //P5ringD_vacLog->SetVisAttributes( Beryllium);

  /*
  
  //Ring E  - simplifying the small articulation in this weldment
  //G4double P5ringE_rin = 11.750/2*inch; 
  G4double P5ringE_rin = 5.88*inch; //CJT
  //G4double P5ringE_rou = 12.0/2*inch;
  G4double P5ringE_rou = 6.5*inch; //CJT
  //G4double P5ringE_L = 0.380*inch;
  G4double P5ringE_L = 4.071/2*inch; //CJT
  P5placement += 2.0*P5ringD_L;

  G4Tubs *P5ringE = new G4Tubs("P5ringE", P5ringE_rin, P5ringE_rou, P5ringE_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringELog = new G4LogicalVolume(P5ringE, GetMaterial("Aluminum"), "P5ringE_log", 0, 0, 0);

  //Place ring E
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5placement+P5ringE_L), P5ringELog, "P5ringELog_pv", worldlog, false, 0, ChkOverlaps);

  //P5ringELog->SetVisAttributes(Aluminum);
  P5ringELog->SetVisAttributes(Beryllium);

  //Ring E Vacuum
  G4Tubs *P5ringE_vac = new G4Tubs("P5ringE_vac", 0.0, P5ringE_rin, P5ringE_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringE_vacLog = new G4LogicalVolume(P5ringE_vac, GetMaterial("Vacuum"), "P5ringE_vac_log", 0, 0, 0);

  //Place ring E vacuum
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P5placement+P5ringE_L), P5ringE_vacLog, "P5ringE_vacLog_pv", worldlog, false, 0 , ChkOverlaps );

  */
  /*
  //TEST RING - should exist in exit beamline. Must delete after check.
  //CJT - total distance from DS end of target to wide flange which begins the exit beamline -> 200.973 inches

  //G4double targetEndOffset_z = 22.8/2.0*inch;
  G4double P5testRing_rin = 5.88*inch; 
  G4double P5testRing_rou = 15.0*inch; //CJT 13.0 - making larger for debug
  G4double P5testRing_L = 0.187/2*inch;
  P5placement = (targetEndOffset_z)+200.973*inch; //Direct measurment from CJT plus half length of target

  G4Tubs *P5testRing = new G4Tubs("P5testRing", P5testRing_rin, P5testRing_rou, P5testRing_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5testRingLog = new G4LogicalVolume(P5testRing, GetMaterial("Air"), "P5testRing_log", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5placement), P5testRingLog, "P5testRingLog_pv", worldlog, false, 0, ChkOverlaps);

  P5testRingLog->SetVisAttributes(G4Colour::Green()); //Debug

  */
  /*
  //TEST RING - should exist in exit beamline. Must delete after check.
  //CJT - total distance from DS end of target to wide flange which begins the exit beamline -> 200.973 inches
  
  G4double P5testRing2_rin = 5.88*inch; 
  G4double P5testRing2_rou = 15.0*inch; //CJT 13.0 - making larger for debug
  G4double P5testRing2_L = 0.187/2*inch;
  P5placement = 212.37*inch; //Direct measurment from CJT plus half length of target

  G4Tubs *P5testRing2 = new G4Tubs("P5testRing2", P5testRing2_rin, P5testRing2_rou, P5testRing2_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5testRing2Log = new G4LogicalVolume(P5testRing2, GetMaterial("Air"), "P5testRing2_log", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5placement), P5testRing2Log, "P5testRing2Log_pv", worldlog, false, 0, ChkOverlaps);

  P5testRing2Log->SetVisAttributes( G4Colour::Green()); //Debug
  */
  //P5ringE_vacLog->SetVisAttributes( G4VisAttributes::GetInvisible());
  /*
  //Ring F
  G4double P5ringF_rin = 11.750/2*inch; 
  G4double P5ringF_rou = 13.960/2*inch;
  //G4double P5ringF_L = 1.120*inch;
  G4double P5ringF_L = 1.120*inch+0.28*inch; //TEMPORARY - adding enough to mate with beam dump pipe

  G4Tubs *P5ringF = new G4Tubs("P5ringF", P5ringF_rin, P5ringF_rou, P5ringF_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringFLog = new G4LogicalVolume(P5ringF, GetMaterial("Aluminum"), "P5ringF_log", 0, 0, 0);

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, P5initPlacement_z+2*P5ringA_L+2*P5ringB_L+2*P5ringC_L+2*P5ringD_L+2*P5ringE_L+P5ringF_L), P5ringFLog, "P5ringFLog_pv", worldlog, false, 0, ChkOverlaps);

  P5ringFLog->SetVisAttributes(Aluminum);

  //Ring F Vacuum
  G4Tubs *P5ringF_vac = new G4Tubs("P5ringF_vac", 0.0, P5ringF_rin, P5ringF_L, 0.*deg, 360.*deg);

  G4LogicalVolume *P5ringF_vacLog = new G4LogicalVolume(P5ringF_vac, GetMaterial("Vacuum"), "P5ringF_vac_log", 0, 0, 0);

  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, P5initPlacement_z+2*P5ringA_L+2*P5ringB_L+2*P5ringC_L+2*P5ringD_L+2*P5ringE_L+P5ringF_L), P5ringF_vacLog, "P5ringF_vacLog_pv", worldlog, false, 0 , ChkOverlaps );
  */
  /*
  //Visuals
  G4VisAttributes *ironColor= new G4VisAttributes(G4Colour(0.3,0.3,0.3));
  //G4VisAttributes *AlColor= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
  //G4VisAttributes *Vacuum_visatt = new G4VisAttributes(G4Colour(0.1, 0.5, 0.9 ) );
  G4VisAttributes *Vacuum_visatt = new G4VisAttributes(G4Colour(0.0, 1.0, 0.0 ) );
  //Vacuum_visatt->SetVisibility(false);  //To determine overlaps, expect green color for vacuum, sseeds
  G4VisAttributes *CopperColor = new G4VisAttributes( G4Colour( 0.7, 0.3, 0.3 ) );
  */

  //===== START PASSING FUNCTIONS - EXIT-BEAMLINE/CORRECTOR-MAGNETS =====//  11.1.20
  
  // Exit beam line piping: use this instead of the commented out section below.  Added by D Flay (Sept 2020)
  //G4double TargetCenter_zoffset = 6.50*inch;
  G4double TargetCenter_zoffset = 0.0*inch;
  
  MakeBeamExit(worldlog,TargetCenter_zoffset);  

  //G4double corrMag1_offset_z = 39.616*inch;  //Offset distance between beginning of conical beamline and first corrector magnet for GEn
  G4double corrMag1_offset_z = 39.925*inch;  //CJT
  G4double corrMag_dz = 68.863*inch;  //CJT distance between corrector magnet 1 and 2 for GEn
  
  MakeCorrectorMagnets(worldlog, P2initPlacement_z+corrMag1_offset_z, corrMag_dz);
  
}





//sseeds Oct 2020 modified from EFuchey CommonExitBeamline code
void G4SBSBeamlineBuilder::MakeCorrectorMagnets(G4LogicalVolume *logicMother, G4double z0, G4double dz){

  //Units, placement vars, and visuals
  bool ChkOverlaps = false;
  G4double inch = 2.54*cm;
  G4VisAttributes *ironColor= new G4VisAttributes(G4Colour(0.3,0.3,0.3));
  G4VisAttributes *CopperColor = new G4VisAttributes(G4Colour(0.7,0.3,0.3));
    //P1tubeA_winvacLog->SetVisAttributes( G4VisAttributes::GetInvisible());
  
  G4double X=0.0, Y=0.0, Z=0.0;
  G4ThreeVector zero(0.0, 0.0, 0.0);

  //Dimension variables for bellows, yoke, and coils
  G4double Rin = 11.750/2.0*inch;
  G4double Rout = 14.0/2.0*inch;
  G4double Thick = 1.12*inch;

  //Bellows dims
  G4double z_formed_bellows = z0;
  G4double Bellows1L = 6.299/2*inch;
  G4double Bellows2L = 15.748/2*inch;

  //Coil dims
  G4double UpstreamCoilThickY = 1.68*inch;
  G4double UpstreamCoilThickX = 3.46*inch;
  G4double UpstreamCoilHeight = 8.17*inch;
  G4double UpstreamCoilDepth = 6.60*inch;
  G4double UpstreamCoilWidth = 7.56*inch;

  G4double DS_coil_depth = 8.91*inch;
  G4double DS_coil_height = 12.04*inch;
  G4double DS_coil_ThickX = 2.90*inch;
  G4double DS_coil_ThickY = 1.68*inch;
  
  //Yoke dims
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

  G4double DownstreamYokeGapWidth = 17.58*inch;
  G4double DownstreamYokeGapHeight = 20.16*inch;

  //Pole dims
  G4double UpstreamPoleDepth = 6.3*inch;
  G4double UpstreamPoleWidth = 4.02*inch;
  G4double UpstreamPoleHeight = 7.87*inch;

  G4double DSpole_depth = 8.76*inch;
  G4double DSpole_width = (17.58-11.00)*inch/2.0;
  G4double DSpole_height = 11.81*inch;

  //Building geometry for upstream coils
  G4Box *UpstreamCoil_outer = new G4Box("UpstreamCoil_outer", UpstreamCoilThickX/2.0, (UpstreamCoilHeight+2.0*UpstreamCoilThickY)/2.0, (UpstreamCoilDepth + 2.0*UpstreamCoilThickY)/2.0 );
  G4Box *UpstreamCoil_inner = new G4Box("UpstreamCoil_inner", UpstreamCoilThickX/2.0 + cm, UpstreamCoilHeight/2.0, UpstreamCoilDepth/2.0 );

  G4SubtractionSolid *UpstreamCoil = new G4SubtractionSolid( "UpstreamCoil", UpstreamCoil_outer, UpstreamCoil_inner, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *UpstreamCoil_log = new G4LogicalVolume(UpstreamCoil, GetMaterial("Copper"), "UpstreamCoil_log" );

  UpstreamCoil_log->SetVisAttributes( CopperColor );
  //UpstreamCoil_log->SetVisAttributes( G4VisAttributes::GetInvisible()); //Debug, keep in
  
  Z = z_formed_bellows+Bellows1L;
  X = (UpstreamCoilWidth+UpstreamCoilThickX)/2.0;
  Y = 0.0;

  //Two placements of upstream coil
  new G4PVPlacement( 0, G4ThreeVector(-X, Y, Z), UpstreamCoil_log, "UpstreamCoil_phys_right", logicMother, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), UpstreamCoil_log, "UpstreamCoil_phys_left", logicMother, false, 1  , ChkOverlaps);

  //Building geometry for upstream poles
  G4Box *UpstreamPole = new G4Box( "UpstreamPole", UpstreamPoleWidth/2.0, UpstreamPoleHeight/2.0, UpstreamPoleDepth/2.0 );
  G4LogicalVolume *UpstreamPole_log = new G4LogicalVolume( UpstreamPole, GetMaterial("Iron"), "UpstreamPole_log" );
  UpstreamPole_log->SetVisAttributes( ironColor );

  //Two placements of upstream poles
  new G4PVPlacement( 0, G4ThreeVector(-X, Y, Z), UpstreamPole_log, "UpstreamPole_phys_right", logicMother, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X, Y, Z), UpstreamPole_log, "UpstreamPole_phys_left", logicMother, false, 1 , ChkOverlaps );

  //Building geometry for yoke top
  G4Box *YokeTopPiece = new G4Box("YokeTopPiece", YokeTopPiece_Width/2.0, YokeTopPiece_Height/2.0, YokeTopPiece_Depth/2.0 );
  G4LogicalVolume *YokeTopPiece_log = new G4LogicalVolume( YokeTopPiece, GetMaterial("Iron"), "YokeTopPiece_log" );

  YokeTopPiece_log->SetVisAttributes( ironColor );
  //YokeTopPiece_log->SetVisAttributes( G4VisAttributes::GetInvisible()); //Debug, keep in
  
  X = 0.0;
  Y = (11.81*inch + YokeTopPiece_Height)/2.0;

  //Two placements of yoke top piece (top and bottom symmetric)
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeTopPiece_log, "UpstreamYokeTop_phys", logicMother, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(X,-Y,Z), YokeTopPiece_log, "UpstreamYokeBottom_phys", logicMother, false, 1 , ChkOverlaps );

  //Building geometry for yoke left
  G4Box *YokeLeftPiece = new G4Box("YokeLeftPiece", YokeLeftPiece_Width/2.0, YokeLeftPiece_Height/2.0, YokeLeftPiece_Depth/2.0 );
  G4LogicalVolume *YokeLeftPiece_log = new G4LogicalVolume( YokeLeftPiece, GetMaterial("Iron"), "YokeLeftPiece_log" );
  YokeLeftPiece_log->SetVisAttributes(ironColor );
  //YokeLeftPiece_log->SetVisAttributes( G4VisAttributes::GetInvisible()); //Debug, keep in
  
  X = 7.52*inch + YokeLeftPiece_Width/2.0;
  Y = 0.0;

  //Placement of yoke left
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeLeftPiece_log, "UpstreamYokeLeftPiece_phys", logicMother, false, 0 , ChkOverlaps );

  //Building geometry for yoke right
  G4Trap *YokeRight_trap = new G4Trap( "YokeRight_trap", YokeRightZFinal/2.0, atan( (YokeRightWidthFinal-YokeRightWidthInitial)/2.0/YokeRightZFinal ), 180.0*deg, YokeLeftPiece_Height/2.0, YokeRightWidthInitial/2.0, YokeRightWidthInitial/2.0, 0.0, YokeLeftPiece_Height/2.0, YokeRightWidthFinal/2.0, YokeRightWidthFinal/2.0, 0.0 ); 

  G4Box *YokeRight_box = new G4Box( "YokeRight_box", YokeRightWidthFinal/2.0, YokeLeftPiece_Height/2.0, 0.39*inch/2.0 );

  X = 0.5*(YokeRightWidthFinal/2.0 + YokeRightWidthInitial/2.0) - YokeRightWidthFinal/2.0;
  
  G4UnionSolid *YokeRightPiece = new G4UnionSolid("YokeRightPiece", YokeRight_trap, YokeRight_box, 0, G4ThreeVector( X, 0, (YokeRightZFinal+0.39*inch)/2.0 ) );
  G4LogicalVolume *YokeRightPiece_log = new G4LogicalVolume(YokeRightPiece, GetMaterial("Iron"), "YokeRightPiece_log" );

  YokeRightPiece_log->SetVisAttributes(ironColor);
  //YokeRightPiece_log->SetVisAttributes( G4VisAttributes::GetInvisible()); //Debug, keep in

  X = -7.52*inch - 0.5*(YokeRightWidthFinal/2.0 + YokeRightWidthInitial/2.0);
  Y = 0.0;
  Z = z_formed_bellows + Bellows1L;

  //Placement of yoke right
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), YokeRightPiece_log, "UpstreamYokeRightPiece_phys", logicMother, false, 0 , ChkOverlaps );

  //Building geometry for downstream (DS) yoke
  G4Box *DownstreamYoke_box = new G4Box("DownstreamYoke_box", DownstreamTotalWidth/2.0, DownstreamTotalHeight/2.0, DownstreamYokeDepth/2.0 );
  G4Box *DownstreamYoke_gap = new G4Box("DownstreamYoke_gap", DownstreamYokeGapWidth/2.0, DownstreamYokeGapHeight/2.0, DownstreamYokeDepth/2.0+cm );
  G4SubtractionSolid *DownstreamYoke = new G4SubtractionSolid( "DownstreamYoke", DownstreamYoke_box, DownstreamYoke_gap, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *DownstreamYoke_log = new G4LogicalVolume( DownstreamYoke, GetMaterial("Iron"), "DownstreamYoke_log" );

  DownstreamYoke_log->SetVisAttributes( ironColor );
  //DownstreamYoke_log->SetVisAttributes( G4VisAttributes::GetInvisible()); //Debug, keep in

  
  //Setting DS Corrector position from input
  z_formed_bellows = z0 + dz;
  X = 0.0; Y = 0.0;
  Z = z_formed_bellows+Bellows2L;

  //Placement of DS yoke
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z), DownstreamYoke_log, "DownstreamYoke_phys", logicMother, false, 0 , ChkOverlaps );
  
  //Building geometry for DS coils
  G4Box *DS_coil_outer = new G4Box( "DS_coil_outer", DS_coil_ThickX/2.0, (DS_coil_height + 2.0*DS_coil_ThickY)/2.0, (DS_coil_depth + 2.0*DS_coil_ThickY)/2.0 );
  G4Box *DS_coil_inner = new G4Box( "DS_coil_inner", DS_coil_ThickX/2.0+cm, DS_coil_height/2.0, DS_coil_depth/2.0 );

  G4SubtractionSolid *DS_coil = new G4SubtractionSolid( "DS_coil", DS_coil_outer, DS_coil_inner, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *DS_coil_log = new G4LogicalVolume( DS_coil, GetMaterial("Copper"), "DS_coil_log" );
  DS_coil_log->SetVisAttributes(CopperColor );
  
  X = 11.67*inch/2.0 + DS_coil_ThickX/2.0;
  Y = 0.0;
  Z = z_formed_bellows+Bellows2L;

  G4double DSCoil_offset_z =(DownstreamYokeDepth-DS_coil_height)/2+1.757*inch;

  //Placement of DS yoke
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z-DSCoil_offset_z), DS_coil_log, "DS_coil_phys_left", logicMother, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(-X,Y,Z-DSCoil_offset_z), DS_coil_log, "DS_coil_phys_right", logicMother, false, 1 , ChkOverlaps );

  //Building geometry for DS poles
  G4Box *DSpole = new G4Box("DSpole", DSpole_width/2.0, DSpole_height/2.0, DSpole_depth/2.0 );
  G4LogicalVolume *DSpole_log = new G4LogicalVolume(DSpole, GetMaterial("Iron"), "DSpole_log" );

  DSpole_log->SetVisAttributes(ironColor);
  
  X = (17.58+11.00)*inch/4.0;
  Y = 0.0;
  
  //Two placements of poles:
  new G4PVPlacement( 0, G4ThreeVector(X,Y,Z-DSCoil_offset_z), DSpole_log, "DSpole_phys_left", logicMother, false, 0 , ChkOverlaps );
  new G4PVPlacement( 0, G4ThreeVector(-X,Y,Z-DSCoil_offset_z), DSpole_log, "DSpole_phys_right", logicMother, false, 1 , ChkOverlaps );

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

void G4SBSBeamlineBuilder::MakeTDISBeamline(G4LogicalVolume *worldlog){// Old beam line...                                                                                                              
  G4bool fOvLap = true;
  //MakeDefaultBeamline(worldlog);
  //add stuff
  // double Wcollimator_length = 35.*mm;
  // double Wcollimator_bore = 7.*mm;
  // double Wcollimator_diameter = 50.*mm;
  // double BeWindowThickness = 20*um;
  
  // G4Tubs *Wcollimator_sol = new G4Tubs("Wcollimator_sol", Wcollimator_bore, Wcollimator_diameter, Wcollimator_length/2, 0.*deg, 360.*deg );
  // G4Tubs *Wcollimatorhole_sol  = new G4Tubs("Wcollimatorhole_sol", 0.0, Wcollimator_bore, Wcollimator_length/2, 0.*deg, 360.*deg );

  // G4LogicalVolume *Wcollimator_log = new G4LogicalVolume(Wcollimator_sol, GetMaterial("TargetBeamCollimator_Material"), "Wcollimatorhole_log", 0, 0, 0);
  // G4LogicalVolume *Wcollimatorhole_log = new G4LogicalVolume(Wcollimatorhole_sol, GetMaterial("Vacuum"), "Wcollimatorhole_log", 0, 0, 0);
  
  // // new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -250.*mm-BeWindowThickness-Wcollimator_length/2.), Wcollimator_log, "Wcollimator_phys", worldlog, false,0);
  // // new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -250.*mm-BeWindowThickness-Wcollimator_length/2.), Wcollimatorhole_log, "Wcollimatorhole_phys", worldlog, false,0);
  
  // G4Tubs *BeWindow_sol = new G4Tubs("BeWindow_sol", 0, Wcollimator_diameter, BeWindowThickness/2, 0.*deg, 360.*deg );
  
  // G4LogicalVolume *BeWindow_log = new G4LogicalVolume(BeWindow_sol, GetMaterial("Beryllium"), "BeWindow_log", 0, 0, 0);

  // new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -250.*mm-BeWindowThickness/2.), BeWindow_log, "BeWindow_phys0", worldlog, false,0);  
  // new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -350.*mm+BeWindowThickness/2.), BeWindow_log, "BeWindow_phys1", worldlog, false,0);  
  
  MakeCommonExitBeamline(worldlog);  
  
  //Downstream beam pipe
  int nsec = 6;
  //  Definition taken from GEN_10M.opc by Bogdan to z = 5.92.  2mm thickness assumed
  G4double exit_z[]   = { 527.*cm, 609.84*cm,609.85*cm, 1161.02*cm, 1161.03*cm,2725.66*cm };
  G4double exit_z_vac[] = { 527.*cm, 610.24*cm,610.35*cm, 1161.52*cm, 1161.53*cm,2726.46*cm };

  G4double exit_zero[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  G4double exit_rin[] = { 14.8*cm,15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };
  G4double exit_rou[] = { 15.0*cm,15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };


  G4Polycone *ext_cone = new G4Polycone("ext_cone", 0.0*deg, 360.0*deg, nsec, exit_z, exit_rin, exit_rou);
  G4Polycone *ext_vac  = new G4Polycone("ext_vac ", 0.0*deg, 360.0*deg, nsec, exit_z_vac, exit_zero, exit_rin);

  G4LogicalVolume *extLog =
    new G4LogicalVolume(ext_cone, GetMaterial("Aluminum"), "ext_log", 0, 0, 0);
  G4LogicalVolume *extvacLog =
    new G4LogicalVolume(ext_vac, GetMaterial("Vacuum"), "extvac_log", 0, 0, 0);

  new G4PVPlacement(0,G4ThreeVector(0, 0, 0.*m), extLog, "ext_phys", worldlog, false,0,fOvLap);
  new G4PVPlacement(0,G4ThreeVector(0, 0, 0.*m), extvacLog, "extvac_phys", worldlog,false,0,fOvLap);

  extvacLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.2,0.6,0.2));
  extLog->SetVisAttributes(pipeVisAtt);
  
  
}


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
  
  extvacLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  entvacLog->SetVisAttributes(G4VisAttributes::GetInvisible());
    
  entvacLog_cut->SetVisAttributes(G4VisAttributes::GetInvisible());
    
  G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
    
  extLog->SetVisAttributes(pipeVisAtt);
  entLog->SetVisAttributes(pipeVisAtt);
  entLog_cut->SetVisAttributes(pipeVisAtt);
}

//  Here is lead shield of beam line for GEp

void G4SBSBeamlineBuilder::MakeGEpLead(G4LogicalVolume *worldlog){

  G4VisAttributes *lead_visatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
  
  G4double inch = 2.54*cm;
  G4double TargetCenter_zoffset = 6.50*inch; //Remove offset - GEp self contained. SSeeds 2021
  //G4double TargetCenter_zoffset = 0.0*inch;

  G4double z_outer_magnetic = 182.33*cm - TargetCenter_zoffset;
  //G4double z_outer_magnetic = 182.33*cm - 6.50*inch;
  
  G4double zstart_lead1 = 170.0*cm; //Remove offset - GEp self contained. SSeeds 2021
  //G4double zstart_lead1 = 170.0*cm+6.50*inch;
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

  ///////////////////// New geometry with vertical wall(s): 

  G4ThreeVector zaxis_temp( -sin(16.9*deg), 0.0, cos(16.9*deg) );
  G4ThreeVector yaxis_temp( 0,1,0);
  G4ThreeVector xaxis_temp = (yaxis_temp.cross(zaxis_temp)).unit();
  
  G4ThreeVector frontcorner_pos = 1.6*m*zaxis_temp;

  G4double temp_shift = 4.5*inch; //Temporary extension of lead wall closest to target 

  G4double zstart_lead_wall1 = z_outer_magnetic + 15*cm - temp_shift; //SSeeds 2021 - Temporary shift pending JT
  G4double zstop_lead_wall1 = zstart_lead_wall1 + 1.25*m + temp_shift; //SSeeds 2021 - Temporary extension pending JT

  //G4Box *lead_wall1 = new G4Box("lead_wall1", 5.0*cm/2.0, 31.0*cm/2.0, 1.25*m/2.0);
  G4Box *lead_wall1 = new G4Box("lead_wall1", 5.0*cm/2.0, 31.0*cm/2.0, 1.25*m/2.0  + temp_shift/2.0);  //SSeeds 2021 - Temporary extension pending JT
  G4LogicalVolume *lead_wall1_log = new G4LogicalVolume( lead_wall1, GetMaterial("Lead"), "lead_wall1_log" );

  G4double xtemp = -( 5.5*inch/2.0 + 1.5*inch + (1.25/2.0+0.15)*m*tan(1.5*deg) + 2.5*cm/cos(1.5*deg) );
  
  rot_temp = new G4RotationMatrix;

  rot_temp->rotateY( 1.5*deg );
  
  //new G4PVPlacement( rot_temp, G4ThreeVector( xtemp, 0.0, zstart_lead_wall1 + 0.5*1.25*m ), lead_wall1_log, "lead_wall1_phys", worldlog, false, 0 );

  G4cout << "Lead wall A (x,y,z) = (" << xtemp/cm << ", " << 0.0 << ", " << (zstart_lead_wall1 + 0.5*1.25*m)/cm << ")" << G4endl;
  
  lead_wall1_log->SetVisAttributes( lead_visatt );

  //first lead wall is lead_wall2 adjacent to iron rings,lead_wall1 is not being used

  G4double zstart_lead_wall2 = z_formed_bellows + 76.09*inch + 1.71*inch + 15.75*inch + 1.0*inch;
  //G4double zstop_lead_wall2 = 207.144*inch - TargetCenter_zoffset + 40.0*inch;
  G4double zstop_lead_wall2 = 207.144*inch - TargetCenter_zoffset + 12.5*inch; //SSeeds 2021 - Temporary extension pending JT

  G4cout << "Lead wall B zstart - zstop = " << (zstop_lead_wall2 - zstart_lead_wall2)/cm << G4endl;
  
  G4double zpos_lead_wall2 = 0.5*(zstart_lead_wall2 + zstop_lead_wall2 - (10.0*inch));
  //zpos offset by 10 inches to fully shield length of beam pipe from line of sight of GEMs
  //we want x position to have x = 
  G4double xpos_lead_wall2 = -(8.0*inch + 2.5*cm + (zpos_lead_wall2 - 201.632*inch + TargetCenter_zoffset )*tan(1.5*deg) + 2.0*inch);
  //xpos_lead_wall2 offset by 2 inches in negative x direction such that wall 3 doesnt clip into the beam pipe, this is done assuming walls 2 and 3 are coupled and offset by 4.5 inches from eachother as given by the model from jlab
  G4Box *lead_wall2 = new G4Box("lead_wall2", 5.0*cm/2.0, 24.0*inch/2.0, 60.0*inch/2.0 );

  G4LogicalVolume *lead_wall2_log = new G4LogicalVolume( lead_wall2, GetMaterial("Lead"), "lead_wall2_log" );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY( 1.5*deg );

  new G4PVPlacement( rot_temp, G4ThreeVector( xpos_lead_wall2, 0, zpos_lead_wall2 ), lead_wall2_log, "lead_wall2_phys", worldlog, false, 0 );

  G4cout << "Lead wall B (x,y,z) = (" << xpos_lead_wall2/cm << ", " << 0.0 << ", " << zpos_lead_wall2/cm << ")" << G4endl;
  
  lead_wall2_log->SetVisAttributes( lead_visatt );
  
  //second lead wall downstream from lead_wall2

  G4Box *lead_wall3 = new G4Box("lead_wall3", 5.0*cm/2.0, 24.0*inch/2.0, 36.0*inch/2.0 );

  G4LogicalVolume *lead_wall3_log = new G4LogicalVolume( lead_wall3, GetMaterial("Lead"), "lead_wall3_log" );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY( 1.5*deg );

  new G4PVPlacement( rot_temp, G4ThreeVector( xpos_lead_wall2-(4.5*inch), 0, zpos_lead_wall2+(30.0*inch)+(18.0*inch)), lead_wall3_log, "lead_wall3_phys", worldlog, false, 0 );

  lead_wall3_log->SetVisAttributes( lead_visatt );

}

//lead shielding for GMn
void G4SBSBeamlineBuilder::MakeGMnLead(G4LogicalVolume *worldlog){

  //SSeeds 1.14.21 - Commenting as obsolete. Will leave here for potential beam studies.
  /**/
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
  G4double w1_SCshield = z1_ringshield*tan(20.0*deg)-mindist_SCshield;
  G4double w2_SCshield = (z1_ringshield+th_SCshield)*tan(20.0*deg)-mindist_SCshield;
  // G4double rin_ringshield = z1_ringshield*sin(6.0*deg);
  // G4double rout_ringshield = z2_ringshield*sin(12.0*deg);
  
  G4double z_SCshield = z1_ringshield+th_SCshield/2.0+0.5*inch;
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
  G4double z1_spoolshield = z1_ringshield+0.5*inch+th_SCshield;
  G4double z2_spoolshield = //z1_spoolshield+1.89*m;
    //z_conic_vacline_weldment + (0.84 + 0.14 + 11.62)*inch;
    z1_spoolshield+36.5*inch;
  //cout << "old length " << (z_conic_vacline_weldment + (0.84 + 0.14 + 11.62)*inch-z1_ringshield-th_SCshield)/inch << endl;
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
  G4double L_BLshield1 = 28.0*inch;//34.0*inch;
  G4double h_BLshield1 = h_SCshield;
  
  G4double z_BLshield1 = z_conic_vacline_weldment + (0.84 + 0.14 + 45.62 + 14.38*0.65 - 34.0/2.0)*inch;
  G4double x_BLshield1 = (12.0)*inch;
  
  x_BLshield1+= th_BLshield1/2.0;
  
  G4Box *BLshield1 = new G4Box("BLshield1", th_BLshield1/2.0, h_BLshield1/2.0, L_BLshield1/2.0 );
  
  G4LogicalVolume *BLshield1_log = new G4LogicalVolume( BLshield1, GetMaterial("Lead"), "BLshield1_log" );
  
  rot_temp = new G4RotationMatrix;
  //rot_temp->rotateY(-1.5*deg);
  
  if(fDetCon->fLeadOption){
    G4cout << "Adding downstream lead plate" << G4endl;
    //if(lead)
    new G4PVPlacement( rot_temp, G4ThreeVector( x_BLshield1, 0, z_BLshield1 ), BLshield1_log, "BLshield1_phys", worldlog, false, 0, checkoverlaps );
    BLshield1_log->SetVisAttributes(LeadColor);
  }
  
  /**/



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
  sideshield_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  
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

  //SSeeds 1.14.21 - Commenting as obsolete
  /*
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
  sideshield_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  
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
  

  */
  
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
  G4double Z2 = Z1 + l_leadwall1/2.0 + l_leadwall2/2.0 - 1*2.54*cm;  

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
  //Add geometry here for GEn lead shielding. None in current build - Jan 2021.
  
}

void G4SBSBeamlineBuilder::MakeSIDISLead( G4LogicalVolume *worldlog ){
  //Add geometry here for SIDIS lead shielding. None in current build - Jan 2021.


}

void G4SBSBeamlineBuilder::MakeALLLead( G4LogicalVolume *worldlog ){
  //Add geometry here for ALL lead shielding.
  
  //adding lead shielding for WAPP2:
  G4cout << " ****************** Adding beamline shielding for SIDIS/WAPP **************" << G4endl;
  bool checkoverlaps = true;
  G4VisAttributes *visShield = new G4VisAttributes();
  visShield->SetColour( G4Colour(0.4,0.4,0.4) );
  // first for radiator: add a 3.75 cm thick, 18cm tall, 25cm long slab of lead, 
  // 15 cm away from BL and as close as possible from tgt collimator A
  G4double radshieldthick = 3.75*cm;
  G4double radshieldlength = 30.0*cm;
  G4double radshieldheight = 20.0*cm;
  G4double radshield_zpos = -30*cm;
  G4double radshield_xposl = 10*cm;
  G4double radshield_xposr = -10*cm;
  G4Box *radshield_solid = new G4Box("radshield_solid", radshieldthick/2.0, radshieldheight/2.0, radshieldlength/2.0);
  
  G4LogicalVolume *radshield_log = new G4LogicalVolume( radshield_solid, GetMaterial("Lead"), "radshield_log" );
  radshield_log->SetVisAttributes(visShield);
  
  //G4RotationMatrix rot_temp = new G4RotationMatrix;

  //new G4PVPlacement( 0, G4ThreeVector( radshield_xposl, 0, radshield_zpos ), radshield_log, "radshield_phys_left", worldlog, false, 0, checkoverlaps );
  
  //do another one for right GEMs
  //new G4PVPlacement( 0, G4ThreeVector( radshield_xposr, 0, radshield_zpos ), radshield_log, "radshield_phys_right", worldlog, false, 0, checkoverlaps );
    
  // second for beamline: add a 5 cm thick, 30cm tall,  200 cm long slab of lead, 
  // 30 cm away from BL and center around z=200cm
  G4double blshieldthick = 5.0*cm;
  G4double blshieldlength = 240.0*cm;
  G4double blshieldheight = 40.0*cm;
  G4double blshield_xpos = 30*cm;
  G4double blshield_zpos = 238*cm;
    
  G4Box *blshield_solid = new G4Box("blshield_solid", blshieldthick/2.0, blshieldheight/2.0, blshieldlength/2.0);
  
  G4LogicalVolume *blshield_log = new G4LogicalVolume( blshield_solid, GetMaterial("Lead"), "blshield_log" );
  
  new G4PVPlacement( 0, G4ThreeVector( blshield_xpos, 0, blshield_zpos ), blshield_log, "blshield_phys", worldlog, false, 0, checkoverlaps );
  blshield_log->SetVisAttributes(visShield);
  
}


void G4SBSBeamlineBuilder::MakeToyBeamline(G4LogicalVolume *motherlog){
  //Add toy beamline and scattering chamber for checking extreme forward angles and effects of other detection parameters. Nothing implemented - Jan 2021
  
}

void G4SBSBeamlineBuilder::MakeBeamDump(G4LogicalVolume *logicMother,G4double dz){
   // build the Hall A beam dump
   // dz = global offset to be applied to all downstream components  
   // Added by D. Flay (JLab) in Sept 2020 

   G4double inch = 2.54*cm; 
   // location of beam diffuser front face relative to target pivot (from Ron Lassiter Sept 2020)
   // target pivot to upstream pipe conical flange:                         1052.4797 inches 
   // upstream pipe conical flange to ISO wall weldment:                    207.1108 inches
   // upstream ISO wall weldment to upstream face of diffuser:              24.56 inches
   // upstream face of beam diffuser to upstream face of downstream flange: 17.3943 inches
   G4double dz3   = -1.5*cm;                   // FIXME: fudge factor to be flush against mid pipe to dump part (see MakeBeamExit, dz3 term)  
   G4double z_us  = dz + dz3 + 1052.4797*inch; 
   G4double z_iso = z_us + 207.1108*inch; 
   G4double z_bd  = z_iso + 24.56*inch; 
   G4double z_ds  = z_bd + 17.3943*inch; 
   // CheckZPos(logicMother,z_bd);
   MakeBeamDump_UpstreamPipe(logicMother,z_us);
   MakeBeamDump_ISOWallWeldment(logicMother,z_iso);

   if( fDetCon->GetBeamDiffuserEnable() ) MakeBeamDump_Diffuser(logicMother,z_bd);
   MakeBeamDump_DownstreamPipe(logicMother,z_ds);
}

void G4SBSBeamlineBuilder::CheckZPos(G4LogicalVolume *logicMother,G4double z0){
   // a dummy function to check positioning
   // z0 = position of DOWNSTREAM face of this part.  All components are spaced relative to this point 

   G4double inch = 2.54*cm; 

   std::cout << "[G4SBSBeamlineBuilder::CheckZPos]: Downstream face of part is at z = " << z0/m << " m" << std::endl; 

   G4double xl = 20.*inch;
   G4double yl = 120.*inch;
   G4double zl = 1.*inch;
   G4Box *solidBox = new G4Box("solidBox",xl/2.,yl/2.,zl/2.); 
   
   G4LogicalVolume *boxLV = new G4LogicalVolume(solidBox,GetMaterial("Aluminum"),"boxLV");

   G4double z = z0 - zl/2.;  // back face is at z0  
   G4ThreeVector P = G4ThreeVector(0,0,z);

   new G4PVPlacement(0,                // no rotation
	             P,                // location in mother volume 
	             boxLV,            // its logical volume                         
	             "testBox_PHY",    // its name
	             logicMother,      // its mother  volume
	             false,            // no boolean operation
	             0,                // copy number
	             true);            // checking overlaps 

}

void G4SBSBeamlineBuilder::MakeBeamDump_Diffuser(G4LogicalVolume *logicMother,G4double z0){
   // A beam diffuser that sits right in front of the beam dump
   // z0 = position of upstream face of this part.  All components are spaced relative to this point 
   // Added by D. Flay (JLab) in Aug 2020  
  
   G4double inch = 25.4*mm;

   G4double BDLength = 0.;

   char Hall = 'A';
   if(Hall=='A') BDLength = 11.44*cm;
   if(Hall=='C') BDLength = 11.30*cm;  

   G4double r_min    = 0.;
   G4double r_max    = 25.*inch;
   G4double len      = BDLength;
   G4double startPhi = 0.*deg;
   G4double dPhi     = 360.*deg;
   G4Tubs *diffCaseS = new G4Tubs("diffCase",r_min,r_max,len/2.,startPhi,dPhi);
   G4LogicalVolume *diffCaseLV = new G4LogicalVolume(diffCaseS,GetMaterial("Vacuum"),"diffCase"); // its name

   // place the diffuser 
   // note: the (x,y) center of the diffuser plates is centered on this logical volume 
   // z0 = location of FRONT FACE of the beam diffuser
   // zz = location of CENTER of the beam diffuser CASE (coincides with the center of the BD)  
   G4double zz = z0 + 0.5*BDLength;
   G4ThreeVector P_case = G4ThreeVector(0,0,zz); 

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
   r_min        = 17.67*inch;
   r_max        = 24.*inch;
   dPhi         = 38.*deg;
   startPhi     = 270.*deg - dPhi/2.;
   G4double thk = 0.125*inch;

   // choose the origin of the device (where the first plate starts, relative to the mother volume) 
   zz = -BDLength/2.;  
   G4ThreeVector P0 = G4ThreeVector(0,0,zz);

   G4VisAttributes *vis = new G4VisAttributes();
   vis->SetColour( G4Colour::Blue() );

   // first plate 
   G4VSolid *plateSolid     = new G4Tubs("plate",r_min,r_max,thk/2.,startPhi,dPhi);
   G4LogicalVolume *plateLV = new G4LogicalVolume(plateSolid,GetMaterial("Aluminum"),"plateLV");
   plateLV->SetVisAttributes(vis);

   // parameterisation
   int NPlanes=0;
   if(Hall=='A') NPlanes = 15;   
   if(Hall=='C') NPlanes = 16;  
   G4VPVParameterisation *plateParam = new G4SBSBDParameterisation(Hall,P0);
   // // placement
   new G4PVParameterised("BeamDiffuser",plateLV,diffCaseLV,kZAxis,NPlanes,plateParam);

   // Attach sensitive detector (SD) functionality; follow the GEM example
   //G4bool bdEnable = fDetCon->GetBeamDiffuserEnable();  // is the beam diffuser enabled? 

   // name of SD and the hitCollection  
   G4String bdSDname = "BD";  // FIXME: is this ok, or do we need directory structure like the GEMs? 
   // We have to remove all the directory structure from the 
   // Hits Collection name or else GEANT4 SDmanager routines will not handle correctly.
   G4String bdSDname_nopath = bdSDname;
   //bdSDname_nopath.remove(0,bdSDname.last('/')+1);
   bdSDname_nopath.erase(0,bdSDname.find_last_of('/')+1);
   G4String bdColName = bdSDname_nopath; 
   bdColName += "HitsCollection";

   G4SBSBeamDiffuserSD *bdSD = nullptr; 
   //if(bdEnable){
   if( !(bdSD = (G4SBSBeamDiffuserSD *)fDetCon->fSDman->FindSensitiveDetector(bdSDname)) ){
     // check to see if this SD exists already; if not, create a new SD object and append to the list of SDs  
     // G4cout << "[G4SBSBeamlineBuilder]: Adding Beam Diffuser SD functionality..." << G4endl;
     G4cout << "Adding Beam Diffuser sensitive detector to SDman..." << G4endl;
     bdSD = new G4SBSBeamDiffuserSD(bdSDname,bdColName);
     plateLV->SetSensitiveDetector(bdSD);  
     fDetCon->fSDman->AddNewDetector(bdSD);
     (fDetCon->SDlist).insert(bdSDname); 
     fDetCon->SDtype[bdSDname] = G4SBS::kBD; 
     // G4cout << "[G4SBSBeamlineBuilder]: --> Done." << G4endl;
   }
      // }
 
}

void G4SBSBeamlineBuilder::MakeBeamDump_ISOWallWeldment(G4LogicalVolume *logicMother,G4double z0){
   // Hall A Beam Dump: ISO Wall Weldment 
   // z0 = position of upstream face of this part.  All components are spaced relative to this point 
   // Drawings: JL0015694, JL0015725, JL0016212 
   // Added by D. Flay (JLab) in Sept 2020  

   G4double inch     = 25.4*mm;
   G4double startPhi = 0.*deg;
   G4double dPhi     = 360.*deg;

   // vacuum tube hub [drawing JL0016212] 
   // - outer component  
   G4double r_min_vth   = 0.5*16.*inch; 
   G4double r_max_vth   = 0.5*20.*inch;
   G4double len_vth     = 2.75*inch; 
   G4Tubs *solidVTH_cyl = new G4Tubs("solidVTH_cyl",r_min_vth,r_max_vth,len_vth/2.,startPhi,dPhi);
   // - cut component   
   G4double r_min_vth_cc = 0.5*16.5*inch; 
   G4double r_max_vth_cc = 0.5*25.0*inch;
   G4double len_vth_cc   = 2.50*inch; 
   G4Tubs *solidVTH_cc = new G4Tubs("solidVTH_cc",r_min_vth_cc,r_max_vth_cc,len_vth_cc/2.,startPhi,dPhi);
   // subtraction 
   G4ThreeVector Pcc = G4ThreeVector(0,0,len_vth_cc-len_vth);
   G4SubtractionSolid *solidVTH = new G4SubtractionSolid("solidVTH",solidVTH_cyl,solidVTH_cc,0,Pcc);

   // wall [drawings JL0015694, JL0015725]  
   G4double x_len = 77.0*inch;   // from JL0015694  
   G4double y_len = 117.0*inch;  // from JL0015694
   G4double z_len = 0.25*inch;   // from JL0015725

   // solid box 
   G4Box *solidBox = new G4Box("isoBox",x_len/2.,y_len/2.,z_len/2.);

   // cut a circular hole 
   G4double r_min    = 0.*mm;
   G4double r_max    = 0.5*16.*inch;
   G4double len      = z_len + 5.*inch; // make sure it cuts through  
   G4Tubs *solidTube = new G4Tubs("isoTube",r_min,r_max,len/2.,startPhi,dPhi);

   // cut the hole in the box
   // - relative to bottom of object, vertical distance is 74.48 inches in drawing JL0015725, 
   //   where the plate is 114.25 inches tall.  Note that this drawing is the inner portion of JL0015694.
   //   However, the height of the whole wall is actually 117 inches (JL0015694), 
   //   so we add 1.375 inches to get 75.86 inches from the bottom of the wall. 
   // - this means we have a distance from the center of the wall: 
   G4double yc     = y_len/2. - (y_len - 75.86*inch); 
   G4ThreeVector P = G4ThreeVector(0.*inch,yc,z_len/2.);
   G4SubtractionSolid *solidWall = new G4SubtractionSolid("solidWall",solidBox,solidTube,0,P);

   // now union the solidVTH and solidWall
   G4double zw = len_vth/2. + z_len/2.;  
   G4ThreeVector Pw = G4ThreeVector(0,-yc,zw); 
   G4UnionSolid *solidISO = new G4UnionSolid("solidISO",solidVTH,solidWall,0,Pw);

   // visualization
   G4VisAttributes *vis = new G4VisAttributes();
   vis->SetColour( G4Colour::Red() );
   // vis->SetForceWireframe();

   // logical volume
   // WARNING: Drawing says Aluminum 6062.  It seems like there is only 6061-T6 alloys.  We use the latter.
   G4LogicalVolume *isoWallLV = new G4LogicalVolume(solidISO,GetMaterial("Aluminum_6061"),"beamDump_isoWall_LV");
   isoWallLV->SetVisAttributes(vis);

   // placement
   bool checkOverlaps = true;  
   G4double z = z0 - 0.5*len_vth; // move center of weldment forward so front face of wall is at z0
   G4ThreeVector P_wall = G4ThreeVector(0,0,z);
   new G4PVPlacement(0,                        // no rotation
                     P_wall,                   // location in mother volume 
                     isoWallLV,                // its logical volume                         
                     "beamDump_isoWall_PHY",   // its name
                     logicMother,              // its mother  volume
                     true,                     // boolean operation? 
                     0,                        // copy number
                     checkOverlaps);          // checking overlaps    

}

void G4SBSBeamlineBuilder::MakeBeamDump_UpstreamPipe(G4LogicalVolume *logicMother,G4double z0){
   // Hall A Beam Dump: Pipe upstream of ISO Weldment
   // z0 = position of upstream face of this part.  All components are spaced relative to this point 
   // Drawing: JL0009934-C-VAC SPOOL REGION UPPER LEVEL
   // Added by D. Flay (JLab) in Sept 2020  

   G4double inch         = 25.4*mm;
   G4double TOTAL_LENGTH = 196.91*inch;
   G4double startPhi     = 0.*deg;
   G4double dPhi         = 360.*deg;

   // large conical tube [item 1]
   G4double delta      = 0.005*inch;  // FIXME: arbitrary!   
   G4double r_min1_lgc = 0.5*23.54*inch;
   G4double r_max1_lgc = r_min1_lgc + delta;
   G4double r_min2_lgc = 0.5*37.54*inch;
   G4double r_max2_lgc = r_min2_lgc + delta;
   G4double len_lgc    =  9.95*inch;
   G4Cons *solidConeLG = new G4Cons("solidConeLG",r_min1_lgc,r_max1_lgc,r_min2_lgc,r_max2_lgc,len_lgc/2.,startPhi,dPhi);

   // large conical tube [item 1, vacuum] 
   G4Cons *solidConeLG_vac = new G4Cons("solidConeLG_vac",0,r_min1_lgc,0,r_min2_lgc,len_lgc/2.,startPhi,dPhi);

   // small conical tube [item 3] 
   G4double r_min1_smc = 0.5*11.99*inch;
   G4double r_max1_smc = r_min1_smc + delta;
   G4double r_min2_smc = 0.5*23.49*inch;
   G4double r_max2_smc = r_min2_smc + delta;
   G4double len_smc    =  6.13*inch;
   G4Cons *solidConeSM = new G4Cons("solidConeSM",r_min1_smc,r_max1_smc,r_min2_smc,r_max2_smc,len_smc/2.,startPhi,dPhi);

   // small conical tube [item 3, vacuum] 
   G4Cons *solidConeSM_vac = new G4Cons("solidConeSM_vac",0,r_min1_smc,0,r_min2_smc,len_smc/2.,startPhi,dPhi);

   // vacuum window tube [item 4] 
   G4double r_min_vwt = 0.5*12.12*inch;
   G4double r_max_vwt = 0.5*12.50*inch;
   G4double len_vwt   = 1.88*inch;
   G4Tubs *solidVWTube = new G4Tubs("solidVWTube",r_min_vwt,r_max_vwt,len_vwt/2.,startPhi,dPhi);

   // vacuum window tube [item 4, vacuum] 
   G4Tubs *solidVWTube_vac = new G4Tubs("solidVWTube_vac",0,r_min_vwt,len_vwt/2.,startPhi,dPhi);

   // main tube [item 6] 
   G4double delta2   = 0.4775*inch;          // fudge factor to get the length to match TOTAL_LENGTH 
   G4double delta3   = 0.005*inch;           // FIXME: This is arbitrary!  
   G4double r_max_m6 = 0.5*24.00*inch;          
   G4double r_min_m6 = r_max_m6 - delta3;      
   G4double len_m6   = 39.85*inch - delta2;  // ESTIMATE: total length of item 6 is 85.36*inch; split into two parts 
   G4Tubs *solidTubeM6 = new G4Tubs("solidTubeM6",r_min_m6,r_max_m6,len_m6/2.,startPhi,dPhi);

   // main tube [item 6, vacuum] 
   G4Tubs *solidTubeM6_vac = new G4Tubs("solidTubeM6_vac",0,r_min_m6,len_m6/2.,startPhi,dPhi);

   // main tube [item 6b, ESTIMATE]  
   G4double r_max_m6b = 0.5*24.00*inch;            
   G4double r_min_m6b = r_max_m6b - delta3;      
   G4double len_m6b   = 45.51*inch - delta2;   // derived number.   
   G4Tubs *solidTubeM6b = new G4Tubs("solidTubeM6b",r_min_m6b,r_max_m6b,len_m6b/2.,startPhi,dPhi);

   // main tube [item 6b, vacuum] 
   G4Tubs *solidTubeM6b_vac = new G4Tubs("solidTubeM6b_vac",0,r_min_m6b,len_m6b/2.,startPhi,dPhi);

   // main tube [item 2]  
   G4double r_max_m2 = 0.5*24.00*inch;      
   G4double r_min_m2 = r_max_m2 - delta3;  
   G4double len_m2   = 52.69*inch - delta2;    // ESTIMATE: total length of item 2 is 91.38*inch; split into two parts 
   G4Tubs *solidTubeM2 = new G4Tubs("solidTubeM2",r_min_m2,r_max_m2,len_m2/2.,startPhi,dPhi);

   // main tube [item 2, vacuum] 
   G4Tubs *solidTubeM2_vac = new G4Tubs("solidTubeM2_vac",0,r_min_m2,len_m2/2.,startPhi,dPhi);

   // main tube [item 2b, ESTIMATE]  
   G4double r_max_m2b = 0.5*24.00*inch;      
   G4double r_min_m2b = r_max_m2b - delta3;  
   G4double len_m2b   = 38.69*inch - delta2;   // derived number.   
   G4Tubs *solidTubeM2b = new G4Tubs("solidTubeM2b",r_min_m2b,r_max_m2b,len_m2b/2.,startPhi,dPhi);

   // main tube [item 2b, vacuum] 
   G4Tubs *solidTubeM2b_vac = new G4Tubs("solidTubeM2b_vac",0,r_min_m2b,len_m2b/2.,startPhi,dPhi);

   // large flange [item 13, drawing JL0012786] 
   G4double r_min_lgf  = 0.5*37.50*inch;
   G4double r_max_lgf  = 0.5*46.00*inch;
   G4double len_lgf    = 0.500*inch;
   G4Tubs *solidLGF = new G4Tubs("solidLGF",r_min_lgf,r_max_lgf,len_lgf/2.,startPhi,dPhi);

   // large flange [item 13, vacuum]  
   G4Tubs *solidLGF_vac = new G4Tubs("solidLGF_vac",0,r_min_lgf,len_lgf/2.,startPhi,dPhi);

   // flange with o-ring [item 9, drawing JL0029536] 
   G4double r_min_for = 0.5*23.25*inch;
   G4double r_max_for = 0.5*29.53*inch;
   G4double len_for   = 1.00*inch;
   G4Tubs *solidFORing = new G4Tubs("solidFORing",r_min_for,r_max_for,len_for/2.,startPhi,dPhi);

   // flange with o-ring [item 9, vacuum] 
   G4Tubs *solidFORing_vac = new G4Tubs("solidFORing_vac",0,r_min_for,len_for/2.,startPhi,dPhi);

   // aperture plate flange [item 14, drawing JL0058855]
   // simplified approach: single material of aluminum 
   G4double r_min_apf = 0.5*3.50*inch;
   G4double r_max_apf = 0.5*30.43*inch;
   G4double len_apf   = 1.50*inch;
   G4Tubs *solidAPF = new G4Tubs("solidAPF",r_min_apf,r_max_apf,len_apf/2.,startPhi,dPhi);

   // vacuum insert into aperture plate flange 
   G4Tubs *solidAPF_vac = new G4Tubs("solidAPF_vac",0,r_min_apf,len_apf/2.,startPhi,dPhi);

   // vacuum tube support ring [item 12] 
   G4double r_min_vtsr = 0.5*24.01*inch;
   G4double r_max_vtsr = 0.5*27.01*inch;
   G4double len_vtsr   = 0.25*inch;
   G4Tubs *solidVTSRing = new G4Tubs("solidVTSRing",r_min_vtsr,r_max_vtsr,len_vtsr/2.,startPhi,dPhi);

   // vacuum tube support ring [item 12, vacuum] 
   G4Tubs *solidVTSRing_vac = new G4Tubs("solidVTSRing_vac",0,r_min_vtsr,len_vtsr/2.,startPhi,dPhi);

   // vacuum step-down flange [item 5, drawing JL0009940] 
   G4double r_min_vsdf = 0.5*12.25*inch;
   G4double r_max_vsdf = 0.5*16.00*inch;
   G4double len_vsdf   = 0.620*inch;
   G4Tubs *solidVSDF = new G4Tubs("solidVSDF",r_min_vsdf,r_max_vsdf,len_vsdf/2.,startPhi,dPhi);

   // vacuum step-down flange [item 5, vacuum] 
   G4Tubs *solidVSDF_vac = new G4Tubs("solidVSDF_vac",0,r_min_vsdf,len_vsdf/2.,startPhi,dPhi);

      // union solid
   // - start with large flange + cone [origin is center of LGF]  
   G4double zz = 0.5*len_lgf + 0.5*len_lgc;
   G4ThreeVector P_0 = G4ThreeVector(0,0,-zz);
   G4UnionSolid *upstrPipe = new G4UnionSolid("lgf_lgc",solidLGF,solidConeLG,0,P_0);
   // - attach m2b 
   zz  = 0.5*len_lgf + len_lgc + 0.5*len_m2b;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b",upstrPipe,solidTubeM2b,0,P_0);
   // - vacuum tube support ring [item 12]  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + 0.5*len_vtsr;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr",upstrPipe,solidVTSRing,0,P_0);
   // - attach m2 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + 0.5*len_m2;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2",upstrPipe,solidTubeM2,0,P_0);
   // - attach aperture plate flange [item 14] 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + 0.5*len_apf;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf",upstrPipe,solidAPF,0,P_0);
   // - attach flange with o-ring [item 9] 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + 0.5*len_for;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for",upstrPipe,solidFORing,0,P_0);
   // - attach m6b 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + 0.5*len_m6b;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b",upstrPipe,solidTubeM6b,0,P_0);
   // - vacuum tube support ring [item 12]  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + 0.5*len_vtsr;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr",upstrPipe,solidVTSRing,0,P_0);
   // - attach m6 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + 0.5*len_m6;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr_m6",upstrPipe,solidTubeM6,0,P_0);
   // - attach small cone 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + len_m6
       + 0.5*len_smc;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr_m6_smc",upstrPipe,solidConeSM,0,P_0);
   // - attach vacuum window tube  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + len_m6
       + len_smc + 0.5*len_vwt;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr_m6_smc_vwt",upstrPipe,solidVWTube,0,P_0);
   // - attach vacuum window stepdown flange  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + len_m6
       + len_smc + len_vwt + 0.5*len_vsdf;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe = new G4UnionSolid("upstreamPipe",upstrPipe,solidVSDF,0,P_0);

   // this is a check  
   G4double TOT_LEN_SUM = len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b
                        + len_vtsr + len_m6 + len_smc + len_vwt + len_vsdf;

   if(TOT_LEN_SUM!=TOTAL_LENGTH){
      std::cout << "[G4SBSBeamlineBuilder::MakeBeamDump_UpstreamPipe]: TOTAL_LENGTH = " << TOTAL_LENGTH/inch
                << " inches, sum of parts = " << TOT_LEN_SUM/inch << " inches" << std::endl;
      exit(1);
   }

   // visualization
   G4VisAttributes *vis = new G4VisAttributes();
   vis->SetColour( G4Colour::Green() );
   // vis->SetForceWireframe(true); 

   // logical volume
   G4LogicalVolume *tubeLV = new G4LogicalVolume(upstrPipe,GetMaterial("Aluminum_5052"),"beamDump_usPipe_LV");
   tubeLV->SetVisAttributes(vis);

   // placement
   bool checkOverlaps = true;  
   G4double z = z0 + 0.5*len_lgf;  // upstream face is at z0 
   G4ThreeVector P = G4ThreeVector(0.,0.,z);
   G4RotationMatrix *rm = new G4RotationMatrix();
   rm->rotateY(180*deg);
   new G4PVPlacement(rm,                       // rotation
                     P,                        // location in mother volume 
                     tubeLV,                   // its logical volume                         
                     "beamDump_usPipe_PHY",    // its name
                     logicMother,              // its mother  volume
                     false,                    // boolean operation? 
                     0,                        // copy number
                     checkOverlaps);          // checking overlaps  


   // union solid [vacuum] 
   // - start with large flange + cone [origin is center of LGF]  
   zz = 0.5*len_lgf + 0.5*len_lgc;
   P_0 = G4ThreeVector(0,0,-zz);
   G4UnionSolid *upstrPipe_vac = new G4UnionSolid("lgf_lgc",solidLGF_vac,solidConeLG_vac,0,P_0);
   // - attach m2b 
   zz  = 0.5*len_lgf + len_lgc + 0.5*len_m2b;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b",upstrPipe_vac,solidTubeM2b_vac,0,P_0);
   // - vacuum tube support ring [item 12]  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + 0.5*len_vtsr;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr",upstrPipe_vac,solidVTSRing_vac,0,P_0);
   // - attach m2 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + 0.5*len_m2;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2",upstrPipe_vac,solidTubeM2_vac,0,P_0);
   // - attach aperture plate flange [item 14] 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + 0.5*len_apf;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf",upstrPipe_vac,solidAPF_vac,0,P_0);
   // - attach flange with o-ring [item 9] 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + 0.5*len_for;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for",upstrPipe_vac,solidFORing_vac,0,P_0);
   // - attach m6b 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + 0.5*len_m6b;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b",upstrPipe_vac,solidTubeM6b_vac,0,P_0);
   // - vacuum tube support ring [item 12]  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + 0.5*len_vtsr;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr",upstrPipe_vac,solidVTSRing_vac,0,P_0);
   // - attach m6 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + 0.5*len_m6;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr_m6",upstrPipe_vac,solidTubeM6_vac,0,P_0);
   // - attach small cone 
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + len_m6
       + 0.5*len_smc;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr_m6_smc",upstrPipe_vac,solidConeSM_vac,0,P_0);
   // - attach vacuum window tube  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + len_m6
       + len_smc + 0.5*len_vwt;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("lgf_lgc_m2b_vtsr_m2_apf_for_m6b_vtsr_m6_smc_vwt",upstrPipe_vac,solidVWTube_vac,0,P_0);
   // - attach vacuum window stepdown flange  
   zz  = 0.5*len_lgf + len_lgc + len_m2b + len_vtsr + len_m2 + len_apf + len_for + len_m6b + len_vtsr + len_m6
       + len_smc + len_vwt + 0.5*len_vsdf;
   P_0 = G4ThreeVector(0,0,-zz);
   upstrPipe_vac = new G4UnionSolid("upstreamPipe_vac",upstrPipe_vac,solidVSDF_vac,0,P_0);

   // visualization
   G4VisAttributes *vis_vac = new G4VisAttributes();
   vis_vac->SetForceWireframe(true);

   // logical volume
   std::string name;
   G4double a=0,density=0;
   G4LogicalVolume *tubeLV_vac = new G4LogicalVolume(upstrPipe_vac,GetMaterial("Vacuum"),"beamDump_usPipe_LV");
   tubeLV_vac->SetVisAttributes(vis_vac);

   new G4PVPlacement(rm,                       // no rotation
                     P,                        // location in mother volume 
                     tubeLV_vac,               // its logical volume                         
                     "beamDump_usPipeVac_PHY", // its name
                     logicMother,              // its mother  volume
                     false,                    // boolean operation? 
                     0,                        // copy number
                     checkOverlaps);          // checking overlaps 

}

void G4SBSBeamlineBuilder::MakeBeamDump_DownstreamPipe(G4LogicalVolume *logicMother,G4double z0){
   // Hall A Beam Dump: Pipe downstream of ISO Weldment
   // z0 = position of upstream face of this part.  All components are spaced relative to this point 
   // Drawing: JL0011756_27020E0145 [modified]
   // Added by D. Flay (JLab) in Sept 2020  

   G4double inch     = 25.4*mm;
   G4double startPhi = 0.*deg;
   G4double dPhi     = 360.*deg;
   G4double TOTAL_LENGTH = 328.04*inch;

   // main tube
   G4double r_min_m  = 0.5*12.11*inch;   // FIXME: This is arbitrary! 
   G4double r_max_m  = 0.5*12.12*inch;
   G4double len_m    = 327.66*inch;
   G4Tubs *solidTubeM = new G4Tubs("solidTube",r_min_m,r_max_m,len_m/2.,startPhi,dPhi);

   // vacuum insert 
   G4Tubs *solidVacuumInsert = new G4Tubs("solidVacuumInsert",0.,r_min_m,TOTAL_LENGTH/2.,startPhi,dPhi);

   // bookend flange [upstream]  
   G4double r_min_us = 0.5*12.11*inch;
   G4double r_max_us = 0.5*14.00*inch;   // FIXME: This is arbitrary!
   G4double len_us   = 0.188*inch;
   G4Tubs *solidTubeUS = new G4Tubs("solidTubeUS",r_min_us,r_max_us,len_us/2.,startPhi,dPhi);

   // bookend flange [downstream]  
   G4double r_min_ds = 0.5*12.11*inch;
   G4double r_max_ds = 0.5*14.00*inch;   // FIXME: This is arbitrary!
   G4double len_ds   = 0.188*inch;
   G4Tubs *solidTubeDS = new G4Tubs("solidTubeDS",r_min_ds,r_max_ds,len_ds/2.,startPhi,dPhi);

   // union solid 
   // - start with upstream flange and main tube 
   G4double zz = 0.5*len_us + 0.5*len_m;
   G4ThreeVector P_0 = G4ThreeVector(0,0,-zz);
   G4UnionSolid *dwnstrPipe = new G4UnionSolid("us_m",solidTubeUS,solidTubeM,0,P_0);
   // - downstream flange 
   zz = 0.5*len_us + len_m + 0.5*len_ds;
   P_0 = G4ThreeVector(0,0,-zz);
   dwnstrPipe = new G4UnionSolid("dwnstrPipe",dwnstrPipe,solidTubeDS,0,P_0);

   // visualization
   G4VisAttributes *vis = new G4VisAttributes();
   vis->SetColour( G4Colour::Yellow() );
   // vis->SetForceWireframe(); 

   // logical volume
   G4LogicalVolume *tubeLV = new G4LogicalVolume(dwnstrPipe,GetMaterial("Aluminum"),"beamDump_dsPipe_LV");
   tubeLV->SetVisAttributes(vis);

   // placement
   bool checkOverlaps = true;
   G4double delta = 5.0*cm; // FIXME: arbitrary!  
   G4double Z = z0 + TOTAL_LENGTH + 0.5*len_us; // upstream face is at z0
   G4ThreeVector P = G4ThreeVector(0,0,Z);
   new G4PVPlacement(0,                        // no rotation
                     P,                        // location in mother volume 
                     tubeLV,                   // its logical volume                         
                     "beamDump_dsPipe_PHY",    // its name
                     logicMother,              // its mother  volume
                     false,                    // boolean operation? 
                     0,                        // copy number
                     checkOverlaps);          // checking overlaps  
   
   // vacuum insert 
   G4VisAttributes *visV = new G4VisAttributes();
   visV->SetForceWireframe();

   std::string name;
   G4double z=0,a=0,density=0;

   // logical volume
   G4LogicalVolume *vacLV = new G4LogicalVolume(solidVacuumInsert,GetMaterial("Vacuum"),"vacuum_dsPipe_LV");
   vacLV->SetVisAttributes(visV);

   Z = z0 + TOTAL_LENGTH/2.;
   P = G4ThreeVector(0,0,Z);
   new G4PVPlacement(0,                    // no rotation
                     P,                    // location in mother volume 
                     vacLV,                // its logical volume                         
                     "vacuum_dsPipe_PHY",  // its name
                     logicMother,          // its mother  volume
                     false,                // boolean operation? 
                     0,                    // copy number
                     checkOverlaps);       // checking overlaps   

}

void G4SBSBeamlineBuilder::MakeBeamExit(G4LogicalVolume *logicMother,G4double dz){
   // Build the Hall A exit beam line  
   // Added by D. Flay (JLab) in Sept 2020
   // Distances from drawing A00000-02-08-0300, A00000-02-08-0700
   // previously: dz set to 7 inches to avoid overlaps with P5ringD_vacLog_pv 

   G4double inch   = 2.54*cm; 
   //G4double dz2    = 1.5*cm; // FIXME: still an overlap despite 6.5" from dz??  adjusting here for simplicity...
   G4double dz2    = 0*mm;            
   G4double z_tmp  = 207.179*inch + dz + dz2;  
   G4double z_mpd  = 749.4997*inch + dz2; // derived, based on numbers from Ron Lassiter (see MakeBeamDump)  
   // CheckZPos(logicMother,z_bd); 
   MakeBeamExit_TargetToMidPipe(logicMother,z_tmp);  
   MakeBeamExit_MidPipeToDump(logicMother,z_mpd); 
   if(fDetCon->GetBeamDumpEnable()) MakeBeamDump(logicMother,dz2); 
}

void G4SBSBeamlineBuilder::MakeBeamExit_TargetToMidPipe(G4LogicalVolume *logicMother,G4double z0){
   // SBS exit beam pipe.  This is immediately upstream of the mid pipe
   // z0 = position of upstream face of this part
   // Added by D. Flay (JLab) in Sept 2020
   // Drawings: 
   // - A00000-02-08-0900
   // - A00000-02-08-0901
   // - 04-65620-E-38500-07_rev 
   // Total length as built here: 438.57*inch (matches first drawing listed)  

   G4double inch     = 2.54*cm; 
   G4double startPhi = 0.*deg; 
   G4double dPhi     = 360.*deg;
   G4double len_fl   = 0.188*inch; // flange thickness  
   G4double wall     = 0.060*inch; // from drawings.  The "wave" part of the corrugation tends to be 0.5"; this seems to be too thick to use  
   G4double TOTAL_LENGTH = 0;  

   // visualization 
   // G4VisAttributes *AlColor = new G4VisAttributes( G4Colour(0.3,0.3,1.0) );
   G4VisAttributes *AlColor = new G4VisAttributes( G4Colour::Green() );
   // G4VisAttributes *vis_vac = new G4VisAttributes( G4Colour(0.1,0.5,0.9) );
   G4VisAttributes *vis_vac = new G4VisAttributes();
  
   // pipe [part 11] 
   G4double r_min_11          = 0.5*39.88*inch;   
   G4double r_max_11          = 0.5*(39.88*inch + 2.*wall); // use wall thickness // 0.5*42.00*inch;  // based on 0.060 x 3 corrugation 
   G4double len_11            = 216.191*inch;    // derived from drawing 
   G4Tubs *solidTube11        = new G4Tubs("solidTube11",r_min_11,r_max_11,len_11/2.,startPhi,dPhi); 
   TOTAL_LENGTH += len_11;
   // vacuum insert  
   G4Tubs *solidTube11_vac    = new G4Tubs("solidTube11_vac",0,r_min_11,len_11/2.,startPhi,dPhi); 

   // flange [part 9] 
   G4double r_min_09          = 0.5*24.00*inch;  
   G4double r_max_09          = 0.5*40.00*inch; 
   G4double len_09            = len_fl;      
   G4Tubs *solidTube09        = new G4Tubs("solidTube09",r_min_09,r_max_09,len_09/2.,startPhi,dPhi); 
   TOTAL_LENGTH += len_09;
   // vacuum insert  
   G4Tubs *solidTube09_vac    = new G4Tubs("solidTube09_vac",0,r_min_09,len_09/2.,startPhi,dPhi); 

   // flange [part 10] 
   G4double r_min_10          = 0.5*36.00*inch;  
   G4double r_max_10          = 0.5*43.00*inch; 
   G4double len_10            = len_fl;      
   G4Tubs *solidTube10        = new G4Tubs("solidTube10",r_min_10,r_max_10,len_10/2.,startPhi,dPhi); 
   TOTAL_LENGTH += len_10; 
   // vacuum insert 
   G4Tubs *solidTube10_vac    = new G4Tubs("solidTube10_vac",0,r_min_10,len_10/2.,startPhi,dPhi); 

   // pipe [part 8] 
   G4double r_min_08          = 0.5*25.88*inch; // estimate to make r_max = 26*inch // 0.5*24.00*inch;   
   G4double r_max_08          = 0.5*(25.88*inch + 2.*wall); // use wall thickness // 0.5*29.334*inch;  // based on 0.5 x 2-2/3 corrugation  
   G4double len_08            = 216.624*inch;     // derived from drawing 
   G4Tubs *solidTube08        = new G4Tubs("solidTube08",r_min_08,r_max_08,len_08/2.,startPhi,dPhi); 
   TOTAL_LENGTH += len_08;
   // vacuum insert 
   G4Tubs *solidTube08_vac    = new G4Tubs("solidTube08_vac",0,r_min_08,len_08/2.,startPhi,dPhi); 

   // flange [part 6] 
   G4double r_min_06          = 0.5*11.75*inch;  
   G4double r_max_06          = 0.5*26.00*inch; 
   G4double len_06            = len_fl;  
   G4Tubs *solidTube06        = new G4Tubs("solidTube06",r_min_06,r_max_06,len_06/2.,startPhi,dPhi); 
   TOTAL_LENGTH += len_06;
   // vacuum insert  
   G4Tubs *solidTube06_vac    = new G4Tubs("solidTube06_vac",0,r_min_06,len_06/2.,startPhi,dPhi); 

   // these two parts related to part 06 are named by D. Flay since there is no "official" 
   // name, but should be included here 
   // 12-in flange [part 6a] 
   // - Drawing A00000-02-08-0900
   G4double r_min_06a         = 0.5*11.28*inch;  
   G4double r_max_06a         = 0.5*16.50*inch; 
   G4double len_06a           = 1.12*inch;   
   G4Tubs *solidTube06a       = new G4Tubs("solidTube06a",r_min_06a,r_max_06a,len_06a/2.,startPhi,dPhi); 
   TOTAL_LENGTH += len_06a;
   // vacuum insert  
   G4Tubs *solidTube06a_vac   = new G4Tubs("solidTube06a_vac",0,r_min_06a,len_06a/2.,startPhi,dPhi); 

   // transfer tube [part 6b] 
   // - Drawing A00000-02-08-0900
   G4double r_min_06b         = 0.5*11.75*inch;  
   G4double r_max_06b         = r_min_06b + 0.5*inch; // FIXME: arbitrary  
   G4double len_06b           = 5.191*inch - len_06a;        
   G4Tubs *solidTube06b       = new G4Tubs("solidTube06b",r_min_06b,r_max_06b,len_06b/2.,startPhi,dPhi); 
   TOTAL_LENGTH += len_06b;
   // vacuum insert  
   G4Tubs *solidTube06b_vac   = new G4Tubs("solidTube06b_vac",0,r_min_06b,len_06b/2.,startPhi,dPhi); 

   // NOTE: these are NOT part of the setup! 
   // // pipe [part 5] 
   // G4double r_min_05          = 0.5*12.00*inch;   
   // G4double r_max_05          = 0.5*17.334*inch; // based on 0.5 x 2-2/3 corrugation  
   // G4double len_05            = 41.8125*inch; 
   // G4Tubs *solidTube05        = new G4Tubs("solidTube05",r_min_05,r_max_05,len_05/2.,startPhi,dPhi); 
   // TOTAL_LENGTH += len_05;
   // // vacuum insert
   // G4Tubs *solidTube05_vac    = new G4Tubs("solidTube05_vac",0,r_min_05,len_05/2.,startPhi,dPhi); 

   // // flange [part 7] 
   // G4double r_min_07          = 0.5*10.00*inch;  
   // G4double r_max_07          = 0.5*13.00*inch; 
   // G4double len_07            = 0.125*inch;      
   // G4Tubs *solidTube07        = new G4Tubs("solidTube07",r_min_07,r_max_07,len_07/2.,startPhi,dPhi); 
   // TOTAL_LENGTH += len_07;
   // // vacuum insert  
   // G4Tubs *solidTube07_vac    = new G4Tubs("solidTube07_vac",0,r_min_07,len_07/2.,startPhi,dPhi); 
  
   // // pipe [part 3]  
   // G4double r_min_03          = 0.5*10.00*inch;   
   // G4double r_max_03          = 0.5*10.25*inch;   
   // G4double len_03            = 12.00*inch; 
   // G4Tubs *solidTube03        = new G4Tubs("solidTube03",r_min_03,r_max_03,len_03/2.,startPhi,dPhi); 
   // TOTAL_LENGTH += len_03;
   // // vacuum insert
   // G4Tubs *solidTube03_vac    = new G4Tubs("solidTube03_vac",0,r_min_03,len_03/2.,startPhi,dPhi); 

   // // flange [part 4] 
   // G4double r_min_04          = 0.5*8.25*inch;  
   // G4double r_max_04          = 0.5*10.00*inch; 
   // G4double len_04            = 0.125*inch;      
   // G4Tubs *solidTube04        = new G4Tubs("solidTube04",r_min_04,r_max_04,len_04/2.,startPhi,dPhi); 
   // TOTAL_LENGTH += len_04;
   // // vacuum insert  
   // G4Tubs *solidTube04_vac    = new G4Tubs("solidTube04_vac",0,r_min_04,len_04/2.,startPhi,dPhi); 

   // // pipe [part 2]  
   // G4double r_min_02          = 0.5*8.00*inch;   
   // G4double r_max_02          = 0.5*8.25*inch;   
   // G4double len_02            = 36.875*inch; 
   // G4Tubs *solidTube02        = new G4Tubs("solidTube02",r_min_02,r_max_02,len_02/2.,startPhi,dPhi); 
   // TOTAL_LENGTH += len_02;
   // // vacuum insert
   // G4Tubs *solidTube02_vac    = new G4Tubs("solidTube02_vac",0,r_min_02,len_02/2.,startPhi,dPhi); 

   // // flange [part 1] 
   // G4double r_min_01          = 0.5*8.25*inch;  
   // G4double r_max_01          = 0.5*10.00*inch; 
   // G4double len_01            = 0.125*inch;      
   // G4Tubs *solidTube01        = new G4Tubs("solidTube01",r_min_01,r_max_01,len_01/2.,startPhi,dPhi); 
   // TOTAL_LENGTH += len_01;
   // // vacuum insert  
   // G4Tubs *solidTube01_vac    = new G4Tubs("solidTube01_vac",0,r_min_01,len_01/2.,startPhi,dPhi);

   // union: put it all together
   // // - start with 01 and 02  
   // G4double zp = 0.5*len_01 + 0.5*len_02;
   // G4ThreeVector P = G4ThreeVector(0,0,zp); 
   // G4UnionSolid *tgtToMidPipe = new G4UnionSolid("t01",solidTube01,solidTube02,0,P);
   // // - add 04 
   // zp = 0.5*len_01 + len_02 + 0.5*len_04;
   // P  = G4ThreeVector(0,0,zp); 
   // tgtToMidPipe = new G4UnionSolid("t04",tgtToMidPipe,solidTube04,0,P);  
   // // - add 03 
   // zp = 0.5*len_01 + len_02 + len_04 + 0.5*len_03;
   // P  = G4ThreeVector(0,0,zp); 
   // tgtToMidPipe = new G4UnionSolid("t03",tgtToMidPipe,solidTube03,0,P);  
   // // - add 07 
   // zp = 0.5*len_01 + len_02 + len_04 + len_03 + 0.5*len_07;
   // P  = G4ThreeVector(0,0,zp); 
   // tgtToMidPipe = new G4UnionSolid("t07",tgtToMidPipe,solidTube07,0,P);  
   // // - add 05 
   // zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + 0.5*len_05;
   // P  = G4ThreeVector(0,0,zp); 
   // tgtToMidPipe = new G4UnionSolid("t05",tgtToMidPipe,solidTube05,0,P); 
   // start with 06a and 06b 
   G4double zp = 0.5*len_06a + 0.5*len_06b;
   G4ThreeVector P  = G4ThreeVector(0,0,zp); 
   G4UnionSolid *tgtToMidPipe = new G4UnionSolid("t06ab",solidTube06a,solidTube06b,0,P);  
   // - add 06 
   zp = 0.5*len_06a + len_06b + 0.5*len_06;
   P  = G4ThreeVector(0,0,zp); 
   tgtToMidPipe = new G4UnionSolid("t06",tgtToMidPipe,solidTube06,0,P);  
   // - add 08 
   zp = 0.5*len_06a + len_06b + len_06 + 0.5*len_08;
   // zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + 0.5*len_08;
   P  = G4ThreeVector(0,0,zp); 
   tgtToMidPipe = new G4UnionSolid("t08",tgtToMidPipe,solidTube08,0,P);  
   // - add 09 
   zp = 0.5*len_06a + len_06b + len_06 + len_08 + 0.5*len_09;
   // zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + len_08 + 0.5*len_09;
   P  = G4ThreeVector(0,0,zp); 
   tgtToMidPipe = new G4UnionSolid("t09",tgtToMidPipe,solidTube09,0,P);  
   // - add 11 
   zp = 0.5*len_06a + len_06b + len_06 + len_08 + len_09 + 0.5*len_11;
   // zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + len_08 + len_09 + 0.5*len_11;
   P  = G4ThreeVector(0,0,zp); 
   tgtToMidPipe = new G4UnionSolid("t11",tgtToMidPipe,solidTube11,0,P);  
   // - add 10 
   zp = 0.5*len_06a + len_06b + len_06 + len_08 + len_09 + len_11 + 0.5*len_10;
   // zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + len_08 + len_09 + len_11 + 0.5*len_10;
   P  = G4ThreeVector(0,0,zp); 
   tgtToMidPipe = new G4UnionSolid("tgtToMidPipe",tgtToMidPipe,solidTube10,0,P); 

   G4LogicalVolume *tgtMP_LV = new G4LogicalVolume(tgtToMidPipe,GetMaterial("Aluminum"),"tgtMP_LV"); 
   tgtMP_LV->SetVisAttributes(AlColor);

   G4double z = z0 + 0.5*len_06a;  // upstream face at z0
   G4ThreeVector Pz = G4ThreeVector(0,0,z);
   
   new G4PVPlacement(0,                // no rotation
                     Pz,               // location in mother volume 
                     tgtMP_LV,         // its logical volume                         
                     "tgtMidPipe_PHY", // its name
                     logicMother,      // its mother  volume
                     true,             // boolean operation? 
                     0,                // copy number
                     true);            // checking overlaps   
   
   // union: put it all together [vacuum] 
   // // - start with 01 and 02  
   // zp = 0.5*len_01 + 0.5*len_02;
   // P = G4ThreeVector(0,0,zp); 
   // G4UnionSolid *tgtToMidPipe_vac = new G4UnionSolid("t01v",solidTube01_vac,solidTube02_vac,0,P);
   // // - add 04 
   // zp = 0.5*len_01 + len_02 + 0.5*len_04;
   // P  = G4ThreeVector(0,0,zp); 
   // tgtToMidPipe_vac = new G4UnionSolid("t04v",tgtToMidPipe_vac,solidTube04_vac,0,P);  
   // // - add 03 
   // zp = 0.5*len_01 + len_02 + len_04 + 0.5*len_03;
   // P  = G4ThreeVector(0,0,zp); 
   // tgtToMidPipe_vac = new G4UnionSolid("t03v",tgtToMidPipe_vac,solidTube03_vac,0,P);  
   // // - add 07 
   // zp = 0.5*len_01 + len_02 + len_04 + len_03 + 0.5*len_07;
   // P  = G4ThreeVector(0,0,zp); 
   // tgtToMidPipe_vac = new G4UnionSolid("t07v",tgtToMidPipe_vac,solidTube07_vac,0,P);  
   // // - add 05 
   // zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + 0.5*len_05;
   // P  = G4ThreeVector(0,0,zp); 
   // tgtToMidPipe_vac = new G4UnionSolid("t05v",tgtToMidPipe_vac,solidTube05_vac,0,P); 
   // start with 06a and 06b 
   zp = 0.5*len_06a + 0.5*len_06b;  
   P  = G4ThreeVector(0,0,zp); 
   G4UnionSolid *tgtToMidPipe_vac = new G4UnionSolid("t06abv",solidTube06a_vac,solidTube06b_vac,0,P);  
   // - add 06 
   zp = 0.5*len_06a + len_06b + 0.5*len_06;
   P  = G4ThreeVector(0,0,zp); 
   tgtToMidPipe_vac = new G4UnionSolid("t06v",tgtToMidPipe_vac,solidTube06_vac,0,P);  
   // - add 08 
   zp = 0.5*len_06a + len_06b + len_06 + 0.5*len_08;
   // zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + 0.5*len_08;
   P  = G4ThreeVector(0,0,zp); 
   tgtToMidPipe_vac = new G4UnionSolid("t08v",tgtToMidPipe_vac,solidTube08_vac,0,P);  
   // - add 09 
   zp = 0.5*len_06a + len_06b + len_06 + len_08 + 0.5*len_09;
   // zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + len_08 + 0.5*len_09;
   P  = G4ThreeVector(0,0,zp); 
   tgtToMidPipe_vac = new G4UnionSolid("t09v",tgtToMidPipe_vac,solidTube09_vac,0,P);  
   // - add 11 
   zp = 0.5*len_06a + len_06b + len_06 + len_08 + len_09 + 0.5*len_11;
   // zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + len_08 + len_09 + 0.5*len_11;
   P  = G4ThreeVector(0,0,zp); 
   tgtToMidPipe_vac = new G4UnionSolid("t11v",tgtToMidPipe_vac,solidTube11_vac,0,P);  
   // - add 10 
   zp = 0.5*len_06a + len_06b + len_06 + len_08 + len_09 + len_11 + 0.5*len_10;
   // zp = 0.5*len_01 + len_02 + len_04 + len_03 + + len_07 + len_05 + len_06 + len_08 + len_09 + len_11 + 0.5*len_10;
   P  = G4ThreeVector(0,0,zp); 
   tgtToMidPipe_vac = new G4UnionSolid("tgtToMidPipe",tgtToMidPipe_vac,solidTube10_vac,0,P); 

   G4LogicalVolume *tgtMP_vac_LV = new G4LogicalVolume(tgtToMidPipe_vac,GetMaterial("Vacuum"),"tgtMP_vac_LV"); 
   tgtMP_vac_LV->SetVisAttributes(vis_vac);
   
   new G4PVPlacement(0,                    // no rotation
                     Pz,                   // location in mother volume 
                     tgtMP_vac_LV,         // its logical volume                         
                     "tgtMidPipe_vac_PHY", // its name
                     logicMother,          // its mother  volume
                     true,                 // boolean operation? 
                     0,                    // copy number
                     true);                // checking overlaps   
   
} 

void G4SBSBeamlineBuilder::MakeBeamExit_MidPipeToDump(G4LogicalVolume *logicMother,G4double z0){
   // SBS exit beam pipe.  This is immediately upstream of the beam dump
   // z0 = position of upstream face of this part
   // Drawings: 
   // - Full assembly: A00000-02-08-0300
   // - Mid pipe to dump: A00000-02-08-700, A00000-02-02-0001 [rev] 
   // Added by D. Flay (JLab) in Sept 2020 

   G4double inch     = 2.54*cm;
   G4double startPhi = 0.*deg; 
   G4double dPhi     = 360.*deg;  
   G4double wall     = 0.060*inch;  

   // visualization 
   // G4VisAttributes *AlColor = new G4VisAttributes( G4Colour(0.3,0.3,1.0) );
   G4VisAttributes *AlColor = new G4VisAttributes( G4Colour(1.0,0.0,0.0) );
   // G4VisAttributes *vis_vac = new G4VisAttributes( G4Colour(0.1,0.5,0.9) );
   G4VisAttributes *vis_vac = new G4VisAttributes();
  
   // from drawing A00000-02-02-0001 [rev]

   G4double delta  = 103.75*inch; // 96.62*inch; // 235.0*mm;     // FIXME: Arbitrary fudge factor to make everything connect from tgtMidPipe to dump  

   G4double r_min  = 0.5*42.88*inch; // 0.5*36.00*inch; // 0.5*42.94*inch;  
   G4double r_max  = 0.5*(42.88*inch + 2.*wall); // 0.5*43.00*inch; 
   G4double len    = 302.98*inch + delta;
   G4Tubs *solidMP = new G4Tubs("solidMP",r_min,r_max,len/2.,startPhi,dPhi); 
    
   G4LogicalVolume *midPipeLV = new G4LogicalVolume(solidMP,GetMaterial("Aluminum"),"midPipeLV");
   midPipeLV->SetVisAttributes(AlColor); 

   // placement 
   G4double Z = z0 + len/2. - delta; // place upstream face at z0 
   G4ThreeVector P = G4ThreeVector(0,0,Z);
   new G4PVPlacement(0,                        // no rotation
                     P,                        // location in mother volume 
                     midPipeLV,                // its logical volume                         
                     "midPipe_toDump_PHY",     // its name
                     logicMother,              // its mother  volume
                     false,                    // boolean operation? 
                     0,                        // copy number
                     true);                    // checking overlaps  

   // fill with vacuum 
   G4Tubs *solidMP_vac            = new G4Tubs("solidMP_vac",0,r_min,len/2.,startPhi,dPhi);
   G4LogicalVolume *midPipeLV_vac = new G4LogicalVolume(solidMP_vac,GetMaterial("Vacuum"),"midPipeLV_vac"); 
   midPipeLV_vac->SetVisAttributes(vis_vac);  
    
   new G4PVPlacement(0,                        // no rotation
                     P,                        // location in mother volume 
                     midPipeLV_vac,            // its logical volume                         
                     "midPipe_toDump_vac_PHY", // its name
                     logicMother,              // its mother  volume
                     false,                    // boolean operation? 
                     0,                        // copy number
                     true);                    // checking overlaps     

}
