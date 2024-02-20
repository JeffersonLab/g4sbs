#include "G4SBSECal.hh"

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
#include "G4Trd.hh"
#include "G4TwoVector.hh"
#include "G4GenericTrap.hh"
#include "G4Polycone.hh"

#include "G4SBSCDet.hh"
#include "G4SBSCalSD.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4SBSRICHSD.hh"
#include "G4SBSECalSD.hh"
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>

#include "sbstypes.hh"

#include "G4PhysicalConstants.hh"

#include "TString.h"

// To supress errors with TString, system of units should be included last
#include "G4SystemOfUnits.hh"

using namespace std;

G4SBSECal::G4SBSECal(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
  fAng = 29.0*deg;
  fDist = 4.9*m;

  fnzsegments_leadglass_ECAL = 1;
  fnzsegments_leadglass_C16 = 1;
  
  assert(fDetCon);
}

G4SBSECal::~G4SBSECal(){;}

void G4SBSECal::BuildComponent(G4LogicalVolume *worldlog){
  G4SBS::Exp_t exptype = fDetCon->fExpType;
  
  if( exptype == G4SBS::kGEp || exptype == G4SBS::kGEPpositron ) //Subsystems unique to the GEp experiment include FPP and BigCal:
    {
      //MakeBigCal( worldlog );
      MakeECal_new( worldlog );
    }
  if( exptype == G4SBS::kC16 ) 
    {
      MakeC16( worldlog );
    }
}


G4LogicalVolume* G4SBSECal::MakeSuperModule( G4double SMWidth, 
					     G4double SMHeight,
					     G4double TiWallLength)
{
  // Make Super Module "shell", which dimensions are parameterized with length and width of the blocks
  // First parameter: Width of the super module; second: length of the super module Ti wall
  
  G4double inch = 2.54*cm;
  
  G4double SMLength = 21.73*inch;
  //G4double SMWidth = 5.11*inch;
  //G4double SMHeight = 5.00*inch;
  
  //G4double TiWallLength = 14.76*inch;
  G4double TiWallThick = 0.032*inch;

  G4double FrontPlateThick = 0.25*inch;
  
  G4double ClampingBarThick = 0.125*inch;
  G4double ClampingBarWidth = 0.78*inch;

  // Mother Volume
  G4Box* SM_mother_box_0 = 
    new G4Box("SM_mother_box_0", SMWidth/2.0, SMHeight/2.0, SMLength/2.0);
  G4Box* SM_mother_box_cut = 
    new G4Box("SM_mother_box_cut", SMWidth/2.0-TiWallThick, SMHeight/2.0+1.0*mm, SMLength/2.0);
  
  G4SubtractionSolid* SM_mother_box = 
    new G4SubtractionSolid("SM_mother_box", SM_mother_box_0, SM_mother_box_cut, 0, G4ThreeVector(0.0, 0.0, -FrontPlateThick-ClampingBarThick));
  
  G4LogicalVolume* SM_mother_log = 
    new G4LogicalVolume(SM_mother_box, GetMaterial("Special_Air"), "SM_mother_log");
  
  // Titanium side walls
  
  G4Box* TiSideWall_box = 
    new G4Box("TiSideWall_box", TiWallThick/2.0, SMHeight/2.0, TiWallLength/2.0);
  G4LogicalVolume* TiSideWall_log = 
    new G4LogicalVolume(TiSideWall_box, GetMaterial("Titanium"), "TiSideWall_log");
  
  G4RotationMatrix* temp_rot = new G4RotationMatrix();
  
  G4ThreeVector temp_pos(0, 0, 0);
  temp_pos.setY(0.0);
  temp_pos.setZ( (SMLength-TiWallLength)/2.0 );
  for(int i = 0; i<2; i++){
    temp_pos.setX( (SMWidth-TiWallThick)/2.0*pow(-1, i) );
    new G4PVPlacement(temp_rot, temp_pos, TiSideWall_log, "TiSideWall_phys", SM_mother_log, false, 0 );
  }
  
  // Front plate
  G4Box* FrontPlate_box = 
    new G4Box("FrontPlate_box", SMWidth/2.0-TiWallThick, SMHeight/2.0, FrontPlateThick/2.0);
  G4LogicalVolume* FrontPlate_log = 
    new G4LogicalVolume(FrontPlate_box, GetMaterial("Aluminum"), "FrontPlate_log");
  
  temp_pos.setX(0.0);
  temp_pos.setY(0.0);
  temp_pos.setZ( (SMLength-FrontPlateThick)/2.0 );
  new G4PVPlacement(temp_rot, temp_pos, FrontPlate_log, "FrontPlate_phys", SM_mother_log, false, 0 );
  
  // Clamping bar
  G4Box* ClampingBar_box = 
    new G4Box("ClampingBar_box", ClampingBarWidth/2.0, ClampingBarWidth/2.0, ClampingBarThick/2.0);
  G4LogicalVolume* ClampingBar_log = 
    new G4LogicalVolume(ClampingBar_box, GetMaterial("Aluminum"), "ClampingBar_log");
  
  temp_pos.setZ( (SMLength-2.0*FrontPlateThick-ClampingBarThick)/2.0 );
  for(int i = 0; i<3; i++){
    for(int j = 0; j<3; j++){
      temp_pos.setX(SMWidth/3.0*(i-1));
      temp_pos.setY(SMHeight/3.0*(j-1));
      new G4PVPlacement(temp_rot, temp_pos, ClampingBar_log, "ClampingBar_phys", SM_mother_log, false, 0 );
    }
  }
  
  //Visualization:
  // G4VisAttributes *mother_visatt = new G4VisAttributes( G4Colour( 1, 1, 1 ) );
  // mother_visatt->SetForceWireframe(true);
  // SM_mother_log->SetVisAttributes(mother_visatt);
  SM_mother_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  // G4VisAttributes *Ti_visatt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0.7 ) );
  // Ti_visatt->SetForceWireframe(true);
  // TiSideWall_log->SetVisAttributes(Ti_visatt);
  TiSideWall_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  // G4VisAttributes *Al_visatt = new G4VisAttributes( G4Colour( 0.7, 0.7, 0.7 ) );
  // Al_visatt->SetForceWireframe(true);
  // FrontPlate_log->SetVisAttributes(Al_visatt);
  // ClampingBar_log->SetVisAttributes(Al_visatt);
  FrontPlate_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  ClampingBar_log->SetVisAttributes(G4VisAttributes::GetInvisible());
  
  return(SM_mother_log);
}
/* KIP CUTOFF//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void G4SBSECal::MakeECal_new(G4LogicalVolume *motherlog){
  // Define the inch
  G4double inch = 2.54*cm;

  // Pointer to SDmanager, used frequently in this routine
  G4SDManager *sdman = fDetCon->fSDman;

  G4double width_42 = 4.2*cm;
  G4double width_40 = 4.0*cm;
  
  G4double depth_42 = 34.3*cm;
  G4double depth_40 = 40.0*cm;
  
  G4double depth_leadglass = 45.0*cm;
  G4double depth_ecal_pmt = 0.3*cm;
  G4double depth_lightguide_short = 20.0*cm;
  G4double diam_lightguide = 2.5*cm;
  G4double radius_ecal_pmt = 1.25*cm;
  // G4double depth_ecal_frontplate = 2.54*cm;
  G4double depth_CH2 = 20.0*cm; //This goes directly in front of CDET:
  G4double depth_CDET = 45.0*cm; // CDET starts 45 cm in front of ECAL:
  
  //Define "Earm" box a bit wider than total ECAL area:
  G4double width_earm = 150.0*cm;
  G4double height_earm = 340.0*cm;
  G4double depth_earm = depth_CH2 + depth_CDET + depth_ecal_pmt + depth_leadglass + depth_lightguide_short; // 
  
  // Rotations, offsets, and distances from Target:
  G4RotationMatrix *bbrm = new G4RotationMatrix;
  bbrm->rotateY(-fAng);
  
  G4Box *earm_mother_box = new G4Box( "earm_mother_box", width_earm/2.0, height_earm/2.0, (depth_earm+1.0*mm)/2.0 );
  G4LogicalVolume *earm_mother_log = new G4LogicalVolume( earm_mother_box, GetMaterial("Air"), "earm_mother_log");

  //If "earm_mother_log" is in the step limiter list, make ECAL a total-absorber with "kCAL" sensitivity:
  if( (fDetCon->StepLimiterList).find( "earm_mother_log" ) != (fDetCon->StepLimiterList).end() ){
    earm_mother_log->SetUserLimits( new G4UserLimits( 0.0, 0.0, 0.0, DBL_MAX, DBL_MAX ) );

    G4String sdname = "Earm/ECAL_box";
    G4String collname = "ECAL_boxHitsCollection";
    G4SBSCalSD *earm_mother_SD = NULL;
    if( !( (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(sdname) ) ){
      G4cout << "Adding ECAL_box sensitive detector to SDman..." << G4endl;
      earm_mother_SD = new G4SBSCalSD( sdname, collname );
      fDetCon->fSDman->AddNewDetector( earm_mother_SD );
      (fDetCon->SDlist).insert( sdname );
      fDetCon->SDtype[sdname] = G4SBS::kCAL;
      (earm_mother_SD->detmap).depth = 0;
      
      earm_mother_log->SetSensitiveDetector( earm_mother_SD );
    }
  }
  
  //fDist should now be interpreted to mean the distance from the origin to the ****FRONT**** of the ECAL lead-glass:
  //G4double zback_ECAL = depth_earm/2.0 - depth_ecal_pmt;
  G4double zfront_ECAL = depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass;
  G4double R_Earm = fDist - zfront_ECAL;
  
  G4ThreeVector pos_ECAL( R_Earm*sin(fAng), 0.0, R_Earm*cos(fAng) );
  
  new G4PVPlacement( bbrm, pos_ECAL, earm_mother_log, "earm_mother_phys", motherlog, false, 0 );
  
  //Blocks
  G4double mylar_thick = 0.001*inch;
  G4double air_thick = mylar_thick;
  
  G4double copper_thick = 0.005*inch;
  G4double al_thick = 0.50*mm;
  G4VisAttributes* hc_visAtt = new G4VisAttributes();
  
  G4double hcf_thick = copper_thick;
  const char* hcf_mat_name = "Copper";
  hc_visAtt->SetColour(1.0, 0.5, 0.0);
  hc_visAtt->SetForceWireframe(true);
  // G4double hcf_thick = al_thick;
  // const char* hcf_mat_name = "Aluminum";
  // hc_visAtt->SetColour(0.7, 0.7, 0.7);
  
  G4Tubs *LightGuide_42 = new G4Tubs("LightGuide_42", 0.0, diam_lightguide/2.0, (depth_lightguide_short+depth_40-depth_42)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LightGuide_42_log = new G4LogicalVolume( LightGuide_42, GetMaterial("Pyrex_Glass"), "LightGuide_42_log" );
  
  G4Tubs *LightGuide_40 = new G4Tubs("LightGuide_40", 0.0, diam_lightguide/2.0, (depth_lightguide_short)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LightGuide_40_log = new G4LogicalVolume( LightGuide_40, GetMaterial("Pyrex_Glass"), "LightGuide_40_log" );

  G4Tubs *LGWrap_42 = new G4Tubs( "LGWrap_42", diam_lightguide/2.0+air_thick, diam_lightguide/2.0 + air_thick + mylar_thick, (depth_lightguide_short+depth_40-depth_42)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LGWrap_42_log = new G4LogicalVolume( LGWrap_42, GetMaterial("Mylar"), "LGWrap_42_log" );

  G4Tubs *LGWrap_40 = new G4Tubs( "LGWrap_40", diam_lightguide/2.0+air_thick, diam_lightguide/2.0 + air_thick + mylar_thick, (depth_lightguide_short)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LGWrap_40_log = new G4LogicalVolume( LGWrap_40, GetMaterial("Mylar"), "LGWrap_40_log" );
  
  G4Tubs *LG42 = new G4Tubs("LG42", 0.0, diam_lightguide/2.0 + air_thick + mylar_thick+0.1*mm, (depth_lightguide_short+depth_40-depth_42)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LG42_log = new G4LogicalVolume( LG42, GetMaterial("Special_Air"), "LG42_log" );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LightGuide_42_log, "LightGuide_42_phys", LG42_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LGWrap_42_log, "LGWrap_42_phys", LG42_log, false, 0 );

  G4Tubs *LG40 = new G4Tubs("LG40", 0.0, diam_lightguide/2.0 + air_thick + mylar_thick+0.1*mm, (depth_lightguide_short)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LG40_log = new G4LogicalVolume( LG40, GetMaterial("Special_Air"), "LG40_log" );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LightGuide_40_log, "LightGuide_40_phys", LG40_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LGWrap_40_log, "LGWrap_40_phys", LG40_log, false, 0 );
  
  width_42 = width_42 + 2*( hcf_thick+mylar_thick+air_thick );
  width_40 = width_40 + 2*( hcf_thick+mylar_thick+air_thick );
  
  G4Box *Module_42 = new G4Box( "Module_42", width_42/2.0, width_42/2.0, depth_42/2.0 );
  G4LogicalVolume *Module_42_log = new G4LogicalVolume( Module_42, GetMaterial("Special_Air"), "Module_42_log" );

  G4Box *Module_40 = new G4Box( "Module_40", width_40/2.0, width_40/2.0, depth_40/2.0 );
  G4LogicalVolume *Module_40_log = new G4LogicalVolume( Module_40, GetMaterial("Special_Air"), "Module_40_log" );

  G4Box *Module_42_k = new G4Box( "Module_42_k", width_42/2.0, width_42/2.0, depth_42/2.0-0.25*mm );
  G4Box *Module_40_k = new G4Box( "Module_42_k", width_40/2.0, width_40/2.0, depth_40/2.0-0.25*mm );
  
  G4Box *hc_42 = new G4Box( "hc_42", width_42/2.0 - hcf_thick, width_42/2.0 - hcf_thick, depth_42/2.0 );
  G4SubtractionSolid *hc_foil_42 = new G4SubtractionSolid( "hc_foil_42", Module_42_k, hc_42, 0, G4ThreeVector( 0, 0, 0 ) );
  G4LogicalVolume *hc_foil_42_log = new G4LogicalVolume( hc_foil_42, GetMaterial(hcf_mat_name), "hc_foil_42_log" );
  hc_foil_42_log->SetVisAttributes(hc_visAtt);

  G4Box *hc_40 = new G4Box( "hc_40", width_40/2.0 - hcf_thick, width_40/2.0 - hcf_thick, depth_40/2.0 );
  G4SubtractionSolid *hc_foil_40 = new G4SubtractionSolid( "hc_foil_40", Module_40_k, hc_40, 0, G4ThreeVector( 0, 0, 0 ) );
  G4LogicalVolume *hc_foil_40_log = new G4LogicalVolume( hc_foil_40, GetMaterial(hcf_mat_name), "hc_foil_40_log" );
  hc_foil_40_log->SetVisAttributes(hc_visAtt);

  //Next, we want to make a subtraction solid for the mylar:
  G4Box *Mylar_42 = new G4Box( "Mylar_42", width_42/2.0 - hcf_thick - mylar_thick, width_42/2.0 - hcf_thick - mylar_thick, depth_42/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_42 = new G4SubtractionSolid( "Mylar_wrap_42", hc_42, Mylar_42, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm ) );
  G4LogicalVolume *Mylar_wrap_42_log = new G4LogicalVolume( Mylar_wrap_42, GetMaterial("Mylar"), "Mylar_wrap_42_log" );
  
  G4Box *Mylar_40 = new G4Box( "Mylar_40", width_40/2.0 - hcf_thick - mylar_thick, width_40/2.0 - hcf_thick - mylar_thick, depth_40/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_40 = new G4SubtractionSolid( "Mylar_wrap_40", hc_40, Mylar_40, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm ) );
  G4LogicalVolume *Mylar_wrap_40_log = new G4LogicalVolume( Mylar_wrap_40, GetMaterial("Mylar"), "Mylar_wrap_40_log" );
  
  ////// Define Sensitive Detector for lead-glass:
  G4String ECalTF1SDname = "Earm/ECalTF1";
  G4String ECalTF1collname = "ECalTF1HitsCollection";
  G4SBSCalSD *ECalTF1SD = NULL;
    
  if( !( ECalTF1SD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(ECalTF1SDname) ) ){
    G4cout << "Adding ECal TF1 Sensitive Detector to SDman..." << G4endl;
    ECalTF1SD = new G4SBSCalSD( ECalTF1SDname, ECalTF1collname );
    fDetCon->fSDman->AddNewDetector( ECalTF1SD );
    (fDetCon->SDlist).insert(ECalTF1SDname);
    fDetCon->SDtype[ECalTF1SDname] = G4SBS::kCAL;
    //fDetCon->SDarm[ECalTF1SDname] = G4SBS::kEarm;

    (ECalTF1SD->detmap).depth = 1;

    fDetCon->SetThresholdTimeWindowAndNTimeBins( ECalTF1SDname, 0.0*MeV, 100.0*ns, 25 );
  }

  fDetCon->InsertSDboundaryVolume( earm_mother_log->GetName(), ECalTF1SDname );

  if( fDetCon->GetC16Segmentation() <= 0 ){
    G4Box *LeadGlass_42 = new G4Box("LeadGlass_42", width_42/2.0 - hcf_thick - mylar_thick - air_thick, width_42/2.0 - hcf_thick - mylar_thick - air_thick, (depth_42 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_42_log = new G4LogicalVolume( LeadGlass_42, GetMaterial("TF1"), "LeadGlass_42_log" );

    G4Box *LeadGlass_40 = new G4Box("LeadGlass_40", width_40/2.0 - hcf_thick - mylar_thick - air_thick, width_40/2.0 - hcf_thick - mylar_thick - air_thick, (depth_40 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_40_log = new G4LogicalVolume( LeadGlass_40, GetMaterial("TF1"), "LeadGlass_40_log" );
    
    //Assign "kCAL" sensitivity to the lead-glass:
    LeadGlass_42_log->SetSensitiveDetector( ECalTF1SD );
    LeadGlass_40_log->SetSensitiveDetector( ECalTF1SD );

    G4VisAttributes *TF1visatt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0 ) );
    LeadGlass_42_log->SetVisAttributes( TF1visatt );
    LeadGlass_40_log->SetVisAttributes( TF1visatt );
    
    //Positioning of lead-glass in module:
    // z + Lz/2 - m/2 - a/2 = Lz/2 --> z = m/2 + a/2
    //lead-glass:
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, (mylar_thick + air_thick)/2.0 ), LeadGlass_42_log, "LeadGlass_42_phys", Module_42_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, (mylar_thick + air_thick)/2.0 ), LeadGlass_40_log, "LeadGlass_40_phys", Module_40_log, false, 0 );
  } else {
    G4int nsegments = fDetCon->GetC16Segmentation();
    G4double segthick = fDetCon->GetSegmentThickC16();

    //Segment the lead-glass longitudinally so that different optical properties can be defined:
    G4int nseg_42 = G4int( (depth_42 - mylar_thick - air_thick)/segthick );
    if( nseg_42 > nsegments ){
      nseg_42 = nsegments;
    }
    G4double remainder_42 = (depth_42 - mylar_thick - air_thick) - nseg_42 * segthick; //remainder is always non-negative!

    if( nseg_42 == 0 ){
      nseg_42 = 1;
      segthick = (depth_42 - mylar_thick - air_thick);
      remainder_42 = 0.0;
    }

    G4Box *segment_42 = new G4Box( "segment_42", (width_42 - 2.*(mylar_thick+air_thick) )/2.0, (width_42 - 2.*(mylar_thick+air_thick) )/2.0, segthick/2.0 );
    G4Box *lastsegment_42 = new G4Box( "lastsegment_42", (width_42 - 2.*(mylar_thick+air_thick) )/2.0, (width_42 - 2.*(mylar_thick+air_thick) )/2.0, (segthick+remainder_42)/2.0 );

    G4double ztemp = -depth_42/2.0 + mylar_thick + air_thick;
    
    for( int seg42=0; seg42<nseg_42; seg42++ ){
      TString material_name, lv_name, pv_name;
      material_name.Form( "TF1_anneal_ECAL_z%d", seg42 );
      lv_name.Form("segment_42_log_%d", seg42 );
      pv_name.Form("segment_42_phys_%d", seg42 );
      G4LogicalVolume *Seg42_log;
      if( seg42 + 1 < nseg_42 || nseg_42 == 1 ){
	Seg42_log = new G4LogicalVolume( segment_42, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += segthick/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg42_log, pv_name.Data(), Module_42_log, false, 0 );
	ztemp += segthick/2.0;
      } else {
	Seg42_log = new G4LogicalVolume( lastsegment_42, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += (segthick+remainder_42)/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg42_log, pv_name.Data(), Module_42_log, false, 0 );
	ztemp += (segthick+remainder_42)/2.0;
      }
      Seg42_log->SetSensitiveDetector( ECalTF1SD );

      G4VisAttributes *Segment_VisAtt = new G4VisAttributes( G4Colour( 0.8*(seg42/15.0)+0.20, 0.8*(seg42/15.0)+0.20, 0.0 ) );
      Seg42_log->SetVisAttributes( Segment_VisAtt );
    }


    //Segment the lead-glass longitudinally so that different optical properties can be defined:
    G4int nseg_40 = G4int( (depth_40 - mylar_thick - air_thick)/segthick );
    if( nseg_40 > nsegments ){
      nseg_40 = nsegments;
    }
    G4double remainder_40 = (depth_40 - mylar_thick - air_thick) - nseg_40 * segthick; //remainder is always non-negative!

    if( nseg_40 == 0 ){
      nseg_40 = 1;
      segthick = (depth_40 - mylar_thick - air_thick);
      remainder_40 = 0.0;
    }
    
    G4Box *segment_40 = new G4Box( "segment_40", (width_40 - 2.*(mylar_thick+air_thick) )/2.0, (width_40 - 2.*(mylar_thick+air_thick) )/2.0, segthick/2.0 );
    G4Box *lastsegment_40 = new G4Box( "lastsegment_40", (width_40 - 2.*(mylar_thick+air_thick) )/2.0, (width_40 - 2.*(mylar_thick+air_thick) )/2.0, (segthick+remainder_40)/2.0 );

    ztemp = -depth_40/2.0 + mylar_thick + air_thick;
    
    for( int seg40=0; seg40<nseg_40; seg40++ ){
      TString material_name, lv_name, pv_name;
      material_name.Form( "TF1_anneal_ECAL_z%d", seg40 );
      lv_name.Form("segment_40_log_%d", seg40 );
      pv_name.Form("segment_40_phys_%d", seg40 );
      G4LogicalVolume *Seg40_log;
      if( seg40 + 1 < nseg_40 || nseg_40 == 1 ){
	Seg40_log = new G4LogicalVolume( segment_40, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += segthick/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg40_log, pv_name.Data(), Module_40_log, false, 0 );
	ztemp += segthick/2.0;
      } else {
	Seg40_log = new G4LogicalVolume( lastsegment_40, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += (segthick+remainder_40)/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg40_log, pv_name.Data(), Module_40_log, false, 0 );
	ztemp += (segthick+remainder_40)/2.0;
      }
      Seg40_log->SetSensitiveDetector( ECalTF1SD );

      G4VisAttributes *Segment_VisAtt = new G4VisAttributes( G4Colour( 0.8*(seg40/15.0)+0.20, 0.8*(seg40/15.0)+0.20, 0.0 ) );
      Seg40_log->SetVisAttributes( Segment_VisAtt );
    }
  }
    
  //Place lead-glass and mylar wrap inside module:
  
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), hc_foil_42_log, "hc_foil_42_phys", Module_42_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), hc_foil_40_log, "hc_foil_40_phys", Module_40_log, false, 0 );
  
  //mylar:
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), Mylar_wrap_42_log, "Mylar_wrap_42_phys", Module_42_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), Mylar_wrap_40_log, "Mylar_wrap_40_phys", Module_40_log, false, 0 );

  new G4LogicalSkinSurface( "Mylar_skin_42", Mylar_wrap_42_log, GetOpticalSurface("Mirrsurf") );
  new G4LogicalSkinSurface( "Mylar_skin_40", Mylar_wrap_40_log, GetOpticalSurface("Mirrsurf") );

  new G4LogicalSkinSurface( "lgwrap_42", LGWrap_42_log, GetOpticalSurface("Mirrsurf") );
  new G4LogicalSkinSurface( "lgwrap_40", LGWrap_40_log, GetOpticalSurface("Mirrsurf") );
  
  //Next: PMTs:
  G4Tubs *ecal_PMT = new G4Tubs( "ecal_PMT", 0.0, radius_ecal_pmt, depth_ecal_pmt/2.0, 0.0, twopi );
  G4LogicalVolume *ecal_PMT_log = new G4LogicalVolume( ecal_PMT, GetMaterial("Photocathode_material_ecal"), "ecal_PMT_log" );

  G4String sdname = "Earm/ECAL";
  G4String collname = "ECALHitsCollection";

  G4SBSECalSD *ECalSD = NULL;

  if( !( (G4SBSECalSD*) sdman->FindSensitiveDetector( sdname ) ) ){
    G4cout << "Adding ECAL PMT sensitive detector to SDman..." << G4endl;
    ECalSD = new G4SBSECalSD( sdname, collname );
    sdman->AddNewDetector( ECalSD );
    (fDetCon->SDlist).insert( sdname );
    fDetCon->SDtype[sdname] = G4SBS::kECAL;
    (ECalSD->detmap).depth = 0;
  }

  ecal_PMT_log->SetSensitiveDetector( ECalSD );
  fDetCon->InsertSDboundaryVolume( earm_mother_log->GetName(), sdname );
  
  // Place the blocks
  G4int NrowsSM_40 = 10; 
  G4int NrowsSM_42 = 15; 
  
  G4int NrowsSM_40 = 0; 
  G4int NrowsSM_42 = 24; 
  */

  //G4int NcolsSM_40[10] = {3, 5, 6, 7, 8, 8, 9, 9, 9, 9};// from bottom to top
  //G4int NcolsSM_42[15] = {9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 7, 6, 5, 3};// from bottom to top
  /*
  G4int NcolsSM_42[24] = {4, 6, 7, 7, 8, 8, 
			  9, 9, 9, 9, 9, 9, 
			  9, 9, 9, 9, 9, 9, 
			  8, 8, 7, 7, 6, 4};// from bottom to top
  */
  //G4double yfp_start_40[10] = {-62.0*cm, -62.0*cm, -62.0*cm, -62.0*cm, -62.0*cm, -58.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm};//start of the block edge

  //G4double yfp_start_42[15] = {-54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -58.0*cm, -62.0*cm, -62.0*cm, -62.0*cm, -62.0*cm, -62.0*cm};//start of the block edge
  /*  
  G4double yfp_start_42[24] = {-54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, 
			       -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -58.0*cm,
			       -62.0*cm, -62.0*cm, -62.0*cm, -62.0*cm, -62.0*cm};//start of the block edge
  */

/* KIP CUTOFF///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  G4double xfpstart = -153.8*cm;
  G4int copy_nb = 0;
  G4double X_block, Y_block;
  
  // Make Super Modules for 42 mm and 40 mm wide blocks
  G4double BlockSpace_42 = 1.69*inch;
  G4double BlockFirst_42 = 0.84*inch;
  
  G4double SMWidth_42_1 = 5.06*inch;
  G4double SMWidth_42_2 = 5.11*inch;
  G4double SMHeight_42 = 5.00*inch;

  G4double TiWallLength_42_1 = 14.83*inch;
  G4double TiWallLength_42_2 = 15.83*inch;
  G4double TiWallThick = 0.032*inch;

  G4double ClampingBarThick = 0.125*inch;
  G4double ClampingBarWidth = 0.78*inch;

  G4double SpacerLength_42_1 = 0.43*inch;
  G4double SpacerLength_42_2 = 1.43*inch;
  G4double SpacerOutDiam = 1.5*inch;
  
  G4double FrontFlangeThick = 0.25*inch;
  G4double BackFlangeThick = 0.5*inch;
  
  G4double BackFlangeHoleSize = 1.37*inch;
  G4double PMTFlangeHoleSize = 1.42*inch;
  
  G4double StandoffLength_42 = 3.00*inch;
  G4double StandoffDiam = 0.375*inch;

  G4double BlockSpace_40 = 1.61*inch;
  G4double BlockFirst_40 = 0.81*inch;
  
  G4double SMWidth_40_1 = 4.84*inch;
  G4double SMWidth_40_2 = 4.89*inch;
  G4double SMHeight_40 = 4.78*inch;

  G4double TiWallLength_40_1 = 17.075*inch;
  G4double TiWallLength_40_2 = 18.075*inch;

  G4double SpacerLength_40_1 = 0.43*inch;
  G4double SpacerLength_40_2 = 1.43*inch;
  
  G4double StandoffLength_40 = 3.00*inch;
  
  //Titanium Walls
  G4Box* TiSideWall_42_1_box = 
    new G4Box("TiSideWall_42_1_box", TiWallThick/2.0, SMHeight_42/2.0, TiWallLength_42_1/2.0);
  G4LogicalVolume* TiSideWall_42_1_log = 
    new G4LogicalVolume(TiSideWall_42_1_box, GetMaterial("Titanium"), "TiSideWall_42_1_log");

  G4Box* TiSideWall_42_2_box = 
    new G4Box("TiSideWall_42_2_box", TiWallThick/2.0, SMHeight_42/2.0, TiWallLength_42_2/2.0);
  G4LogicalVolume* TiSideWall_42_2_log = 
    new G4LogicalVolume(TiSideWall_42_2_box, GetMaterial("Titanium"), "TiSideWall_42_2_log");
  
  G4Box* TiSideWall_40_1_box = 
    new G4Box("TiSideWall_40_1_box", TiWallThick/2.0, SMHeight_40/2.0, TiWallLength_40_1/2.0);
  G4LogicalVolume* TiSideWall_40_1_log = 
    new G4LogicalVolume(TiSideWall_40_1_box, GetMaterial("Titanium"), "TiSideWall_40_1_log");

  G4Box* TiSideWall_40_2_box = 
    new G4Box("TiSideWall_40_2_box", TiWallThick/2.0, SMHeight_40/2.0, TiWallLength_40_2/2.0);
  G4LogicalVolume* TiSideWall_40_2_log = 
    new G4LogicalVolume(TiSideWall_40_2_box, GetMaterial("Titanium"), "TiSideWall_40_2_log");

  // Flanges:
  // Front
  G4Box* FrontFlange_42_1_box = 
    new G4Box("FrontFlange_42_1_box", SMWidth_42_1/2.0, SMHeight_42/2.0, FrontFlangeThick/2.0);
  G4Box* FrontFlange_42_2_box = 
    new G4Box("FrontFlange_42_2_box", SMWidth_42_2/2.0, SMHeight_42/2.0, FrontFlangeThick/2.0);
  
  G4Box* ClampingBar_box = 
    new G4Box("ClampingBar_box", ClampingBarWidth/2.0, ClampingBarWidth/2.0, ClampingBarThick/2.0);
  
  G4UnionSolid* FrontFlange_42_1_solid = 
    new G4UnionSolid("FrontFlange_42_1_solid", FrontFlange_42_1_box, ClampingBar_box, 0, 
		     G4ThreeVector(0.0, 0.0, (FrontFlangeThick+ClampingBarThick)/2.0));
  
  G4UnionSolid* FrontFlange_42_2_solid = 
    new G4UnionSolid("FrontFlange_42_2_solid", FrontFlange_42_2_box, ClampingBar_box, 0, 
		     G4ThreeVector(0.0, 0.0, (FrontFlangeThick+ClampingBarThick)/2.0));
  
  G4Box* FrontFlange_40_1_box = 
    new G4Box("FrontFlange_40_1_box", SMWidth_40_1/2.0, SMHeight_40/2.0, FrontFlangeThick/2.0);
  G4Box* FrontFlange_40_2_box = 
    new G4Box("FrontFlange_40_2_box", SMWidth_40_2/2.0, SMHeight_40/2.0, FrontFlangeThick/2.0);
  
  G4UnionSolid* FrontFlange_40_1_solid = 
    new G4UnionSolid("FrontFlange_40_1_solid", FrontFlange_40_1_box, ClampingBar_box, 0, 
		     G4ThreeVector(0.0, 0.0, (FrontFlangeThick+ClampingBarThick)/2.0));
  
  G4UnionSolid* FrontFlange_40_2_solid = 
    new G4UnionSolid("FrontFlange_40_2_solid", FrontFlange_40_2_box, ClampingBar_box, 0, 
		     G4ThreeVector(0.0, 0.0, (FrontFlangeThick+ClampingBarThick)/2.0));

  // Back
  G4Box* BackFlange_42_1_box = 
    new G4Box("BackFlange_42_1_box", SMWidth_42_1/2.0, SMHeight_42/2.0, BackFlangeThick/2.0);
  
  G4Box* BackFlange_42_2_box = 
    new G4Box("BackFlange_42_2_box", SMWidth_42_2/2.0, SMHeight_42/2.0, BackFlangeThick/2.0);
  
  G4Tubs *BackFlange_hole = new G4Tubs("BackFlange_hole", 0.0, BackFlangeHoleSize/2.0, BackFlangeThick, 0.0*deg, 360.0*deg );
  
  
  G4SubtractionSolid* BackFlange_42_1_solid = 
    new G4SubtractionSolid("BackFlange_42_1_solid", BackFlange_42_1_box, BackFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  
  G4SubtractionSolid* BackFlange_42_2_solid = 
    new G4SubtractionSolid("BackFlange_42_2_solid", BackFlange_42_2_box, BackFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  
  G4Box* BackFlange_40_1_box = 
    new G4Box("BackFlange_40_1_box", SMWidth_40_1/2.0, SMHeight_40/2.0, BackFlangeThick/2.0);
  
  G4Box* BackFlange_40_2_box = 
    new G4Box("BackFlange_40_2_box", SMWidth_40_2/2.0, SMHeight_40/2.0, BackFlangeThick/2.0);
  
  G4SubtractionSolid* BackFlange_40_1_solid = 
    new G4SubtractionSolid("BackFlange_40_1_solid", BackFlange_40_1_box, BackFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  
  G4SubtractionSolid* BackFlange_40_2_solid = 
    new G4SubtractionSolid("BackFlange_40_2_solid", BackFlange_40_2_box, BackFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  
  // PMT
  G4Tubs *PMTFlange_hole = new G4Tubs("PMTFlange_hole", 0.0, PMTFlangeHoleSize/2.0, FrontFlangeThick, 0.0*deg, 360.0*deg );
  
  G4SubtractionSolid* PMTFlange_42_1_solid = 
    new G4SubtractionSolid("PMTFlange_42_1_solid", FrontFlange_42_1_box, PMTFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  
  G4SubtractionSolid* PMTFlange_42_2_solid = 
    new G4SubtractionSolid("PMTFlange_42_2_solid", FrontFlange_42_2_box, PMTFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  
  G4SubtractionSolid* PMTFlange_40_1_solid = 
    new G4SubtractionSolid("PMTFlange_40_1_solid", FrontFlange_40_1_box, PMTFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  
  G4SubtractionSolid* PMTFlange_40_2_solid = 
    new G4SubtractionSolid("PMTFlange_40_2_solid", FrontFlange_40_2_box, PMTFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  
  // Flanges boolean solids
  for(int i = 0; i<3; i++){
    for(int j = 0; j<3; j++){
      FrontFlange_42_1_solid = 
	new G4UnionSolid("FrontFlange_42_1_solid", FrontFlange_42_1_solid, ClampingBar_box, 0, 
			 G4ThreeVector(BlockSpace_42*(i-1),BlockSpace_42*(j-1), (FrontFlangeThick+ClampingBarThick)/2.0));
      FrontFlange_42_2_solid = 
	new G4UnionSolid("FrontFlange_42_1_solid", FrontFlange_42_2_solid, ClampingBar_box, 0, 
			 G4ThreeVector(BlockSpace_42*(i-1),BlockSpace_42*(j-1), (FrontFlangeThick+ClampingBarThick)/2.0));
      
      
       BackFlange_42_1_solid = new G4SubtractionSolid("BackFlange_42_1_solid", BackFlange_42_1_solid, BackFlange_hole, 0, G4ThreeVector(BlockSpace_42*(i-1),BlockSpace_42*(j-1),0.0));
      
       BackFlange_42_2_solid = new G4SubtractionSolid("BackFlange_42_2_solid", BackFlange_42_2_solid, BackFlange_hole, 0, G4ThreeVector(BlockSpace_42*(i-1),BlockSpace_42*(j-1),0.0));
      
      PMTFlange_42_1_solid = new G4SubtractionSolid("PMTFlange_42_1_solid", PMTFlange_42_1_solid, PMTFlange_hole, 0, G4ThreeVector(BlockSpace_42*(i-1),BlockSpace_42*(j-1),0.0));
      
      PMTFlange_42_2_solid = new G4SubtractionSolid("PMTFlange_42_2_solid", PMTFlange_42_2_solid, PMTFlange_hole, 0, G4ThreeVector(BlockSpace_42*(i-1),BlockSpace_42*(j-1),0.0));

      FrontFlange_40_1_solid = 
	new G4UnionSolid("FrontFlange_40_1_solid", FrontFlange_40_1_solid, ClampingBar_box, 0, 
			 G4ThreeVector(BlockSpace_40*(i-1),BlockSpace_40*(j-1), (FrontFlangeThick+ClampingBarThick)/2.0));
      FrontFlange_40_2_solid = 
	new G4UnionSolid("FrontFlange_40_1_solid", FrontFlange_40_2_solid, ClampingBar_box, 0, 
			 G4ThreeVector(BlockSpace_40*(i-1),BlockSpace_40*(j-1), (FrontFlangeThick+ClampingBarThick)/2.0));
      
      
       BackFlange_40_1_solid = new G4SubtractionSolid("BackFlange_40_1_solid", BackFlange_40_1_solid, BackFlange_hole, 0, G4ThreeVector(BlockSpace_40*(i-1),BlockSpace_40*(j-1),0.0));
      
       BackFlange_40_2_solid = new G4SubtractionSolid("BackFlange_40_2_solid", BackFlange_40_2_solid, BackFlange_hole, 0, G4ThreeVector(BlockSpace_40*(i-1),BlockSpace_40*(j-1),0.0));
      
      PMTFlange_40_1_solid = new G4SubtractionSolid("PMTFlange_40_1_solid", PMTFlange_40_1_solid, PMTFlange_hole, 0, G4ThreeVector(BlockSpace_40*(i-1),BlockSpace_40*(j-1),0.0));
      
      PMTFlange_40_2_solid = new G4SubtractionSolid("PMTFlange_40_2_solid", PMTFlange_40_2_solid, PMTFlange_hole, 0, G4ThreeVector(BlockSpace_40*(i-1),BlockSpace_40*(j-1),0.0));
    }
  }
  
  G4LogicalVolume* FrontFlange_42_1_log = 
    new G4LogicalVolume(FrontFlange_42_1_solid, GetMaterial("Aluminum"), "FrontFlange_42_1_log");
  G4LogicalVolume* FrontFlange_42_2_log = 
    new G4LogicalVolume(FrontFlange_42_2_solid, GetMaterial("Aluminum"), "FrontFlange_42_2_log");
  
  G4LogicalVolume* BackFlange_42_1_log = 
    new G4LogicalVolume(BackFlange_42_1_solid, GetMaterial("Aluminum"), "BackFlange_42_1_log");
  G4LogicalVolume* BackFlange_42_2_log = 
    new G4LogicalVolume(BackFlange_42_2_solid, GetMaterial("Aluminum"), "BackFlange_42_2_log");
  
  G4LogicalVolume* PMTFlange_42_1_log = 
    new G4LogicalVolume(PMTFlange_42_1_solid, GetMaterial("Aluminum"), "PMTFlange_42_1_log");
  G4LogicalVolume* PMTFlange_42_2_log = 
    new G4LogicalVolume(PMTFlange_42_2_solid, GetMaterial("Aluminum"), "PMTFlange_42_2_log");

  G4LogicalVolume* FrontFlange_40_1_log = 
    new G4LogicalVolume(FrontFlange_40_1_solid, GetMaterial("Aluminum"), "FrontFlange_40_1_log");
  G4LogicalVolume* FrontFlange_40_2_log = 
    new G4LogicalVolume(FrontFlange_40_2_solid, GetMaterial("Aluminum"), "FrontFlange_40_2_log");
  
  G4LogicalVolume* BackFlange_40_1_log = 
    new G4LogicalVolume(BackFlange_40_1_solid, GetMaterial("Aluminum"), "BackFlange_40_1_log");
  G4LogicalVolume* BackFlange_40_2_log = 
    new G4LogicalVolume(BackFlange_40_2_solid, GetMaterial("Aluminum"), "BackFlange_40_2_log");
  
  G4LogicalVolume* PMTFlange_40_1_log = 
    new G4LogicalVolume(PMTFlange_40_1_solid, GetMaterial("Aluminum"), "PMTFlange_40_1_log");
  G4LogicalVolume* PMTFlange_40_2_log = 
    new G4LogicalVolume(PMTFlange_40_2_solid, GetMaterial("Aluminum"), "PMTFlange_40_2_log");

  // Spacers
  G4Tubs *Spacer_42_1_solid = new G4Tubs("Spacer_42_1_solid", BackFlangeHoleSize/2.0, SpacerOutDiam/2.0, SpacerLength_42_1/2.0, 0.0*deg, 360.0*deg );
  G4Tubs *Spacer_42_2_solid = new G4Tubs("Spacer_42_2_solid", BackFlangeHoleSize/2.0, SpacerOutDiam/2.0, SpacerLength_42_2/2.0, 0.0*deg, 360.0*deg );
   
  G4LogicalVolume* Spacer_42_1_log = 
    new G4LogicalVolume(Spacer_42_1_solid, 
			GetMaterial("Aluminum"), "Spacer_42_1_log");
  G4LogicalVolume* Spacer_42_2_log = 
    new G4LogicalVolume(Spacer_42_2_solid, 
			GetMaterial("Aluminum"), "Spacer_42_2_log");
  
  G4Tubs *Spacer_40_1_solid = new G4Tubs("Spacer_40_1_solid", BackFlangeHoleSize/2.0, SpacerOutDiam/2.0, SpacerLength_40_1/2.0, 0.0*deg, 360.0*deg );
  G4Tubs *Spacer_40_2_solid = new G4Tubs("Spacer_40_2_solid", BackFlangeHoleSize/2.0, SpacerOutDiam/2.0, SpacerLength_40_2/2.0, 0.0*deg, 360.0*deg );
   
  G4LogicalVolume* Spacer_40_1_log = 
    new G4LogicalVolume(Spacer_40_1_solid, 
			GetMaterial("Aluminum"), "Spacer_40_1_log");
  G4LogicalVolume* Spacer_40_2_log = 
    new G4LogicalVolume(Spacer_40_2_solid, 
			GetMaterial("Aluminum"), "Spacer_40_2_log");
  
  //Standoffs
  G4Tubs *Standoff_42_solid = new G4Tubs("Standoff_42_solid", 0, StandoffDiam/2.0, StandoffLength_42/2.0, 0.0*deg, 360.0*deg );
  
  G4LogicalVolume* Standoff_42_log = 
    new G4LogicalVolume(Standoff_42_solid, 
			GetMaterial("Aluminum"), "Standoff_42_log");
  
  G4Tubs *Standoff_40_solid = new G4Tubs("Standoff_40_solid", 0, StandoffDiam/2.0, StandoffLength_40/2.0, 0.0*deg, 360.0*deg );
  
  G4LogicalVolume* Standoff_40_log = 
    new G4LogicalVolume(Standoff_40_solid, 
			GetMaterial("Aluminum"), "Standoff_40_log");

  ofstream mapfile("database/ecal_gep_blockmap.txt");

  TString  output_line = "";
  output_line.Form("#%16s, %16s, %16s, %16s, %16s", "Cell", "Row", "Column", "Xcenter (cm)", "Ycenter (cm)" ); 

  mapfile << output_line << endl;
  
  G4int SM_num = 0;
  X_block = xfpstart+BlockFirst_40;
  for(int i_ = 0; i_<NrowsSM_40*3; i_++){
    Y_block = yfp_start_40[i_/3]+TiWallThick+BlockFirst_40;
    for(int j_ = 0; j_<NcolsSM_40[i_/3]*3; j_++){
      //printf("i_ = %d, j = %d, i_/3, copy_nb = %d, X_block = %f, Y_Block = %f \n", i_, j_, i_/3, copy_nb, X_block, Y_block);
      G4ThreeVector modpos( Y_block, X_block, zfront_ECAL + depth_40/2.0 );
      // depth_earm/2.0 - depth_ecal_pmt - depth_leadglass/2.0 );
      new G4PVPlacement( 0, modpos, Module_40_log, "Module_40_phys", earm_mother_log, false, copy_nb );

      output_line.Form( "%16d, %16d, %16d, %16.3f, %16.3f", copy_nb, i_, j_, Y_block/cm, X_block/cm );
      mapfile << output_line << endl;
      
      (ECalTF1SD->detmap).Row[copy_nb] = i_;
      (ECalTF1SD->detmap).Col[copy_nb] = j_;
      (ECalTF1SD->detmap).LocalCoord[copy_nb] = modpos;
      
      G4ThreeVector pmtpos( modpos.x(), modpos.y(), 
			    modpos.z() + depth_40/2.0 + LG40->GetZHalfLength()*2 + depth_ecal_pmt/2.0 );
      new G4PVPlacement( 0, pmtpos, ecal_PMT_log, "ecal_PMT_phys", earm_mother_log, false, copy_nb );
      
      (ECalSD->detmap).Row[copy_nb] = i_;
      (ECalSD->detmap).Col[copy_nb] = j_;
      (ECalSD->detmap).LocalCoord[copy_nb] = pmtpos;
      
      //Add light-guide with mylar wrap:
      G4ThreeVector LGpos( modpos.x(), modpos.y(), modpos.z() + depth_40/2.0 + LG40->GetZHalfLength() );
      new G4PVPlacement( 0, LGpos, LG40_log, "LG40_phys", earm_mother_log, false, copy_nb );
      
      // //EFuchey 2017-01-12: Need to make sensitive the three volumes above, to measure their dose.
      // // shall be temporary, and not end in the repo...
      // (ECalLGSD->detmap).Row[copy_nb] = global_row;
      // (ECalLGSD->detmap).Col[copy_nb] = col;
      // (ECalLGSD->detmap).LocalCoord[copy_nb] = modpos;
      
      // Placing the supermodule
      // Spacers...
      G4double z0_SM = zfront_ECAL - FrontFlangeThick - ClampingBarThick;
      G4RotationMatrix* rm_temp = new G4RotationMatrix();
      G4ThreeVector pos_temp(Y_block, X_block, z0_SM + TiWallLength_40_1-BackFlangeThick-SpacerLength_40_1/2.0);
      new G4PVPlacement(rm_temp, pos_temp, Spacer_40_1_log, "Spacer_40_1_phys", earm_mother_log, false, 0 );
      
      if(i_%3==1 && j_%3==1){
	// Front Flange
	pos_temp.setZ(z0_SM + FrontFlangeThick/2.0);
	new G4PVPlacement(rm_temp, pos_temp, FrontFlange_40_1_log, "FrontFlange_40_1_phys", earm_mother_log, false, 0 );
	
	// Back Flange
	pos_temp.setZ(z0_SM + TiWallLength_40_1 - BackFlangeThick/2.0);
	new G4PVPlacement(rm_temp, pos_temp, BackFlange_40_1_log, "BackFlange_40_1_phys", earm_mother_log, false, 0 );
	
	// PMT Flange
	pos_temp.setZ(z0_SM + TiWallLength_40_1 + StandoffLength_40 + FrontFlangeThick/2.0);
	new G4PVPlacement(rm_temp, pos_temp, PMTFlange_40_1_log, "PMTFlange_40_1_phys", earm_mother_log, false, 0 );
	
	// Ti walls and Standoffs
	for(int k = 0; k<2; k++){
	  pos_temp.setZ(z0_SM + TiWallLength_40_1/2.0);
	  pos_temp.setY(X_block);
	  pos_temp.setX(Y_block+(k-0.5)*(SMWidth_40_1+TiWallThick/2.0));
	  new G4PVPlacement(rm_temp, pos_temp, TiSideWall_40_1_log, "TiSideWall_40_1_phys", 
			    earm_mother_log, false, 0 );
	  for(int l = 0; l<2; l++){
	    pos_temp.setX(Y_block+BlockSpace_40*(l-0.5));
	    pos_temp.setY(X_block+(SMHeight_40-StandoffDiam)*(k-0.5));
	    pos_temp.setZ(z0_SM + TiWallLength_40_1 + StandoffLength_40/2.0);
	    new G4PVPlacement(rm_temp, pos_temp, Standoff_40_log, "Standoff_40_phys", 
			      earm_mother_log, false, 0 );
	  }
	}
	SM_num++;
	//if(NcolsSM_40[i_/3]%2==0 && j_>(NcolsSM_40[i_/3]-1)*3)SM_num++;
      }
      
      Y_block+= BlockSpace_40;
      if(j_%3==2)Y_block+= 2*BlockFirst_40+2*TiWallThick-BlockSpace_40;
      copy_nb++;
    }
    X_block+= BlockSpace_40;
  }
  
  X_block+= BlockFirst_42+BlockFirst_40-BlockSpace_40;
  for(int i_ = 0; i_<NrowsSM_42*3; i_++){
    Y_block = yfp_start_42[i_/3]+TiWallThick+BlockFirst_42;
    for(int j_ = 0; j_<NcolsSM_42[i_/3]*3; j_++){
      //printf("i_ = %d, j = %d, i_/3, copy_nb = %d, X_block = %f, Y_Block = %f \n", i_, j_, i_/3, copy_nb, X_block, Y_block);
      G4ThreeVector modpos( Y_block, X_block, zfront_ECAL + depth_42/2.0 );
      // depth_earm/2.0 - depth_ecal_pmt - depth_leadglass/2.0 );
      new G4PVPlacement( 0, modpos, Module_42_log, "Module_42_phys", earm_mother_log, false, copy_nb );

      output_line.Form( "%16d, %16d, %16d, %16.3f, %16.3f", copy_nb, i_ + NrowsSM_40*3, j_, Y_block/cm, X_block/cm );
      mapfile << output_line << endl;
      
      (ECalTF1SD->detmap).Row[copy_nb] = i_ + 3*NrowsSM_40;
      (ECalTF1SD->detmap).Col[copy_nb] = j_;
      (ECalTF1SD->detmap).LocalCoord[copy_nb] = modpos;

      G4ThreeVector pmtpos( modpos.x(), modpos.y(),  
			    modpos.z() + depth_42/2.0 + LG42->GetZHalfLength()*2 + depth_ecal_pmt/2.0 );
      new G4PVPlacement( 0, pmtpos, ecal_PMT_log, "ecal_PMT_phys", earm_mother_log, false, copy_nb );

      (ECalSD->detmap).Row[copy_nb] = i_ + 3*NrowsSM_40;
      (ECalSD->detmap).Col[copy_nb] = j_;
      (ECalSD->detmap).LocalCoord[copy_nb] = pmtpos;

      //Add light-guide with mylar wrap:
      G4ThreeVector LGpos( modpos.x(), modpos.y(), modpos.z() + depth_42/2.0 + LG42->GetZHalfLength() );
      new G4PVPlacement( 0, LGpos, LG42_log, "LG42_phys", earm_mother_log, false, copy_nb );

      // //EFuchey 2017-01-12: Need to make sensitive the three volumes above, to measure their dose.
      // // shall be temporary, and not end in the repo...
      // (ECalLGSD->detmap).Row[copy_nb] = global_row;
      // (ECalLGSD->detmap).Col[copy_nb] = col;
      // (ECalLGSD->detmap).LocalCoord[copy_nb] = modpos;
      
      // Placing the supermodule
      // Spacers...
      G4double z0_SM = zfront_ECAL - FrontFlangeThick - ClampingBarThick;
      G4RotationMatrix* rm_temp = new G4RotationMatrix();
      G4ThreeVector pos_temp(Y_block, X_block, z0_SM + TiWallLength_42_1-BackFlangeThick-SpacerLength_42_1/2.0);
      new G4PVPlacement(rm_temp, pos_temp, Spacer_42_1_log, "Spacer_42_1_phys", earm_mother_log, false, 0 );
      
      if(i_%3==1 && j_%3==1){
	// Front Flange
	pos_temp.setZ(z0_SM + FrontFlangeThick/2.0);
	new G4PVPlacement(rm_temp, pos_temp, FrontFlange_42_1_log, "FrontFlange_42_1_phys", earm_mother_log, false, 0 );
	
	// Back Flange
	pos_temp.setZ(z0_SM + TiWallLength_42_1 - BackFlangeThick/2.0);
	new G4PVPlacement(rm_temp, pos_temp, BackFlange_42_1_log, "BackFlange_42_1_phys", earm_mother_log, false, 0 );
	
	// PMT Flange
	pos_temp.setZ(z0_SM + TiWallLength_42_1 + StandoffLength_42 + FrontFlangeThick/2.0);
	new G4PVPlacement(rm_temp, pos_temp, PMTFlange_42_1_log, "PMTFlange_42_1_phys", earm_mother_log, false, 0 );
	
	// Ti walls and Standoffs
	for(int k = 0; k<2; k++){
	  pos_temp.setZ(z0_SM + TiWallLength_42_1/2.0);
	  pos_temp.setY(X_block);
	  pos_temp.setX(Y_block+(k-0.5)*(SMWidth_42_1+TiWallThick/2.0));
	  new G4PVPlacement(rm_temp, pos_temp, TiSideWall_42_1_log, "TiSideWall_42_1_phys", 
			    earm_mother_log, false, 0 );
	  for(int l = 0; l<2; l++){
	    pos_temp.setX(Y_block+BlockSpace_42*(l-0.5));
	    pos_temp.setY(X_block+(SMHeight_42-StandoffDiam)*(k-0.5));
	    pos_temp.setZ(z0_SM + TiWallLength_42_1 + StandoffLength_42/2.0);
	    new G4PVPlacement(rm_temp, pos_temp, Standoff_42_log, "Standoff_42_phys", 
			      earm_mother_log, false, 0 );
	  }
	}
	SM_num++;
	//if(NcolsSM_42[i_/3]%2==0 && j_>(NcolsSM_42[i_/3]-1)*3)SM_num++;
      }
      Y_block+= BlockSpace_42;
      if(j_%3==2)Y_block+= 2*BlockFirst_42+2*TiWallThick-BlockSpace_42;
      copy_nb++;
    }
    X_block+= BlockSpace_42;
    if(i_==2)X_block+= 2.0*mm;// to leave space for the steel plate instered there
  }

  mapfile.close();
  
  // Build CDet:
  if(fDetCon->fExpType==G4SBS::kGEp || fDetCon->fExpType == G4SBS::kGEPpositron ){
    //Next: CH2 filter:
    G4Box *CH2_filter = new G4Box( "CH2_filter", width_earm/2.0, height_earm/2.0, depth_CH2/2.0 );
    G4LogicalVolume *CH2_filter_log = new G4LogicalVolume( CH2_filter, GetMaterial("Polyethylene"), "CH2_filter_log" );
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, zfront_ECAL - depth_CDET - depth_CH2/2.0 ), CH2_filter_log, "CH2_filter_phys", earm_mother_log, false, 0 );
    
    G4double z0_CDET = -depth_earm/2.0 + depth_CH2;
    //G4double R0_CDET = R_Earm - depth_leadglass - depth_CDET;
    //G4double R0_CDET = fDist - depth_leadglass - depth_CDET;
    G4double R0_CDET = fDist - depth_CDET;
    
    G4SBSCDet* CDet = new G4SBSCDet(fDetCon);
    CDet->SetArmName("Earm");
    CDet->SetR0(R0_CDET);
    CDet->SetZ0(z0_CDET);
    CDet->SetPlanesHOffset(0.0);
    CDet->SetPlanesInterDistance(1.0*cm);
    CDet->BuildComponent( earm_mother_log );
    //MakeCDET( earm_mother_log );
    
    G4VisAttributes *CH2_visatt = new G4VisAttributes( G4Colour( 0, 0.6, 0.6 ) );
    CH2_visatt->SetForceWireframe(true);
    
    CH2_filter_log->SetVisAttributes(CH2_visatt);
  }

  //Visualization:
  G4VisAttributes *mother_visatt = new G4VisAttributes( G4Colour( 1, 1, 1 ) );
  mother_visatt->SetForceWireframe(true);
  earm_mother_log->SetVisAttributes(mother_visatt);
   
  //earm_mother_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  Module_42_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  Module_40_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  
  G4VisAttributes *Mylarvisatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
  Mylarvisatt->SetForceWireframe(true);
  Mylar_wrap_42_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  Mylar_wrap_40_log->SetVisAttributes( G4VisAttributes::GetInvisible() );

  G4VisAttributes *Ti_visatt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0.7 ) );
  Ti_visatt->SetForceWireframe(true);
  TiSideWall_42_1_log->SetVisAttributes(Ti_visatt);
  TiSideWall_42_2_log->SetVisAttributes(Ti_visatt);
  TiSideWall_40_1_log->SetVisAttributes(Ti_visatt);
  TiSideWall_40_2_log->SetVisAttributes(Ti_visatt);
  // TiSideWall_log->SetVisAttributes(G4VisAttributes::GetInvisible());//
  
  G4VisAttributes *Al_visatt = new G4VisAttributes( G4Colour( 0.7, 0.7, 0.7 ) );
  Al_visatt->SetForceWireframe(true);
  FrontFlange_42_1_log->SetVisAttributes(Al_visatt);
  FrontFlange_42_2_log->SetVisAttributes(Al_visatt);
  FrontFlange_40_1_log->SetVisAttributes(Al_visatt);
  FrontFlange_40_2_log->SetVisAttributes(Al_visatt);
  BackFlange_42_1_log->SetVisAttributes(Al_visatt);
  BackFlange_42_2_log->SetVisAttributes(Al_visatt);
  BackFlange_40_1_log->SetVisAttributes(Al_visatt);
  BackFlange_40_2_log->SetVisAttributes(Al_visatt);
  PMTFlange_42_1_log->SetVisAttributes( G4Colour( 0.0, 0.7, 0.0 ) );
  PMTFlange_42_2_log->SetVisAttributes( G4Colour( 0.7, 0.0, 0.0 ) );
  PMTFlange_40_1_log->SetVisAttributes( G4Colour( 0.0, 0.7, 0.0 ) );
  PMTFlange_40_2_log->SetVisAttributes( G4Colour( 0.7, 0.0, 0.0 ) );
  Spacer_42_1_log->SetVisAttributes( G4VisAttributes::GetInvisible() ); //G4Colour( 0.0, 0.7, 0.0 ));
  Spacer_42_2_log->SetVisAttributes( G4VisAttributes::GetInvisible() ); //G4Colour( 0.0, 0.7, 0.0 ));
  Spacer_40_1_log->SetVisAttributes( G4VisAttributes::GetInvisible() ); //G4Colour( 0.0, 0.7, 0.0 ));
  Spacer_40_2_log->SetVisAttributes( G4VisAttributes::GetInvisible() ); //G4Colour( 0.0, 0.7, 0.0 ));
  Standoff_42_log->SetVisAttributes( G4VisAttributes::GetInvisible() ); //Al_visatt);
  Standoff_40_log->SetVisAttributes( G4VisAttributes::GetInvisible() ); //Al_visatt);
  
  
  G4VisAttributes *ECALpmtvisatt = new G4VisAttributes( G4Colour( 0, 0, 1 ) );
  ecal_PMT_log->SetVisAttributes( ECALpmtvisatt );
}


*/
//KIP CUTOFF///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Workspace for Kip updating the ECal Geometry, this is an edited version of what is above for the MakeECal_new code, if messed up, uncomment where its labelled KIP CUTOFF and comment out everything in the KIP workspace below
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////KIP
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////KIP
void G4SBSECal::MakeECal_new(G4LogicalVolume *motherlog){
  // Define the inch
  G4double inch = 2.54*cm;

  // Pointer to SDmanager, used frequently in this routine
  G4SDManager *sdman = fDetCon->fSDman;

  G4double width_42 = 4.25*cm;
  G4double width_40 = 4.0*cm;
  
  G4double depth_42 = 34.3*cm;
  G4double depth_40 = 40.0*cm;
  
  G4double depth_leadglass = 45.0*cm;
  G4double depth_ecal_pmt = 0.3*cm;
  G4double depth_lightguide_short = 20.0*cm;
  G4double diam_lightguide = 2.5*cm;
  G4double radius_ecal_pmt = 1.25*cm;
  G4double absorber_depth = 7.0*2.54*cm;
  // G4double depth_ecal_frontplate = 2.54*cm;
  G4double depth_CH2 = 20.0*cm; //This goes directly in front of CDET:
  G4double depth_CDET = 45.0*cm; // CDET starts 45 cm in front of ECAL:
  
  //Define "Earm" box a bit wider than total ECAL area:
  G4double width_earm = 150.0*cm;
  G4double height_earm = 340.0*cm;
  G4double depth_earm = depth_CH2 + depth_CDET + depth_ecal_pmt + depth_leadglass + depth_lightguide_short + absorber_depth; //   

  // Rotations, offsets, and distances from Target:
  G4RotationMatrix *bbrm = new G4RotationMatrix;
  bbrm->rotateY(-fAng);
  
  G4Box *earm_mother_box = new G4Box( "earm_mother_box", width_earm/2.0, height_earm/2.0, (depth_earm+1.0*mm)/2.0 );
  G4LogicalVolume *earm_mother_log = new G4LogicalVolume( earm_mother_box, GetMaterial("Air"), "earm_mother_log");

  //If "earm_mother_log" is in the step limiter list, make ECAL a total-absorber with "kCAL" sensitivity:
  if( (fDetCon->StepLimiterList).find( "earm_mother_log" ) != (fDetCon->StepLimiterList).end() ){
    earm_mother_log->SetUserLimits( new G4UserLimits( 0.0, 0.0, 0.0, DBL_MAX, DBL_MAX ) );

    G4String sdname = "Earm/ECAL_box";
    G4String collname = "ECAL_boxHitsCollection";
    G4SBSCalSD *earm_mother_SD = NULL;
    if( !( (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(sdname) ) ){
      G4cout << "Adding ECAL_box sensitive detector to SDman..." << G4endl;
      earm_mother_SD = new G4SBSCalSD( sdname, collname );
      fDetCon->fSDman->AddNewDetector( earm_mother_SD );
      (fDetCon->SDlist).insert( sdname );
      fDetCon->SDtype[sdname] = G4SBS::kCAL;
      (earm_mother_SD->detmap).depth = 0;
      
      earm_mother_log->SetSensitiveDetector( earm_mother_SD );
    }
  }
  
  //fDist should now be interpreted to mean the distance from the origin to the ****FRONT**** of the ECAL lead-glass:
  //G4double zback_ECAL = depth_earm/2.0 - depth_ecal_pmt;
  G4double zfront_ECAL = depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass;
  //?since depth_earm has been modified to account for the absorber, does this(zfront_ECAL) need to be modifified as well?
  G4double R_Earm = fDist - zfront_ECAL;
  
  G4ThreeVector pos_ECAL( R_Earm*sin(fAng), 0.0, R_Earm*cos(fAng) );
  
  new G4PVPlacement( bbrm, pos_ECAL, earm_mother_log, "earm_mother_phys", motherlog, false, 0 );
  
  //Blocks
  G4double mylar_thick = 0.001*inch;
  G4double air_thick = mylar_thick;   

  G4double copper_thick = 0.005*inch;
  G4double al_thick = 0.50*mm;
  G4VisAttributes* hc_visAtt = new G4VisAttributes();
  
  G4double hcf_thick = copper_thick;
  const char* hcf_mat_name = "Copper";
  hc_visAtt->SetColour(1.0, 0.5, 0.0);
  hc_visAtt->SetForceWireframe(true);
  // G4double hcf_thick = al_thick;
  // const char* hcf_mat_name = "Aluminum";
  // hc_visAtt->SetColour(0.7, 0.7, 0.7);
  

  //Creating Light Guides

  G4Tubs *LightGuide_42 = new G4Tubs("LightGuide_42", 0.0, diam_lightguide/2.0, (depth_lightguide_short+depth_40-depth_42)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LightGuide_42_log = new G4LogicalVolume( LightGuide_42, GetMaterial("Pyrex_Glass"), "LightGuide_42_log" );
  //creates solid lightguide volume   

  G4Tubs *LGWrap_42 = new G4Tubs( "LGWrap_42", diam_lightguide/2.0+air_thick, diam_lightguide/2.0 + air_thick + mylar_thick, (depth_lightguide_short+depth_40-depth_42)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LGWrap_42_log = new G4LogicalVolume( LGWrap_42, GetMaterial("Mylar"), "LGWrap_42_log" );
  //creates cylindrical shell of mylar 

  G4Tubs *LG42 = new G4Tubs("LG42", 0.0, diam_lightguide/2.0 + air_thick + mylar_thick+0.1*mm, (depth_lightguide_short+depth_40-depth_42)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LG42_log = new G4LogicalVolume( LG42, GetMaterial("Special_Air"), "LG42_log" );
  //creates logical volume to contain both lg and mylar

  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LightGuide_42_log, "LightGuide_42_phys", LG42_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LGWrap_42_log, "LGWrap_42_phys", LG42_log, false, 0 );
  //places lg and mylar in physical volume

  width_42 = width_42 + 2*( hcf_thick + mylar_thick + air_thick );
  //this should come out to be 4.25 cm
  
  G4Box *Module_42 = new G4Box( "Module_42", width_42/2.0, width_42/2.0, depth_42/2.0 );
  G4LogicalVolume *Module_42_log = new G4LogicalVolume( Module_42, GetMaterial("Special_Air"), "Module_42_log" );

  G4Box *Module_42_k = new G4Box( "Module_42_k", width_42/2.0, width_42/2.0, depth_42/2.0-0.25*mm );
  //why the 0.25 subtraction?

  G4Box *hc_42 = new G4Box( "hc_42", width_42/2.0 - hcf_thick, width_42/2.0 - hcf_thick, depth_42/2.0 );
  G4SubtractionSolid *hc_foil_42 = new G4SubtractionSolid( "hc_foil_42", Module_42_k, hc_42, 0, G4ThreeVector( 0, 0, 0 ) );
  G4LogicalVolume *hc_foil_42_log = new G4LogicalVolume( hc_foil_42, GetMaterial(hcf_mat_name), "hc_foil_42_log" );
  hc_foil_42_log->SetVisAttributes(hc_visAtt);
  //defining the volume of the the hc foil

  G4Box *Mylar_42 = new G4Box( "Mylar_42", width_42/2.0 - hcf_thick - mylar_thick, width_42/2.0 - hcf_thick - mylar_thick, depth_42/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_42 = new G4SubtractionSolid( "Mylar_wrap_42", hc_42, Mylar_42, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm ) );
  G4LogicalVolume *Mylar_wrap_42_log = new G4LogicalVolume( Mylar_wrap_42, GetMaterial("Mylar"), "Mylar_wrap_42_log" );
  //defining the volume of the mylar wrapping


  ////// Define Sensitive Detector for lead-glass:
  G4String ECalTF1SDname = "Earm/ECalTF1";
  G4String ECalTF1collname = "ECalTF1HitsCollection";
  G4SBSCalSD *ECalTF1SD = NULL;
    
  if( !( ECalTF1SD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(ECalTF1SDname) ) ){
    G4cout << "Adding ECal TF1 Sensitive Detector to SDman..." << G4endl;
    ECalTF1SD = new G4SBSCalSD( ECalTF1SDname, ECalTF1collname );
    fDetCon->fSDman->AddNewDetector( ECalTF1SD );
    (fDetCon->SDlist).insert(ECalTF1SDname);
    fDetCon->SDtype[ECalTF1SDname] = G4SBS::kCAL;
    //fDetCon->SDarm[ECalTF1SDname] = G4SBS::kEarm;

    (ECalTF1SD->detmap).depth = 1;

    fDetCon->SetThresholdTimeWindowAndNTimeBins( ECalTF1SDname, 0.0*MeV, 100.0*ns, 25 );
  }

  fDetCon->InsertSDboundaryVolume( earm_mother_log->GetName(), ECalTF1SDname );

  if( fDetCon->GetC16Segmentation() <= 0 ){
    G4Box *LeadGlass_42 = new G4Box("LeadGlass_42", width_42/2.0 - hcf_thick - mylar_thick - air_thick, width_42/2.0 - hcf_thick - mylar_thick - air_thick, (depth_42 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_42_log = new G4LogicalVolume( LeadGlass_42, GetMaterial("TF1"), "LeadGlass_42_log" );

    //Assign "kCAL" sensitivity to the lead-glass:
    LeadGlass_42_log->SetSensitiveDetector( ECalTF1SD );

    G4VisAttributes *TF1visatt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0 ) );
    LeadGlass_42_log->SetVisAttributes( TF1visatt );
    
    //Positioning of lead-glass in module:
    // z + Lz/2 - m/2 - a/2 = Lz/2 --> z = m/2 + a/2
    //lead-glass:
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, (mylar_thick + air_thick)/2.0 ), LeadGlass_42_log, "LeadGlass_42_phys", Module_42_log, false, 0 );

  } else {
    G4int nsegments = fDetCon->GetC16Segmentation();
    G4double segthick = fDetCon->GetSegmentThickC16();

    //Segment the lead-glass longitudinally so that different optical properties can be defined:
    G4int nseg_42 = G4int( (depth_42 - mylar_thick - air_thick)/segthick );
    if( nseg_42 > nsegments ){
      nseg_42 = nsegments;
    }
    G4double remainder_42 = (depth_42 - mylar_thick - air_thick) - nseg_42 * segthick; //remainder is always non-negative!

    if( nseg_42 == 0 ){
      nseg_42 = 1;
      segthick = (depth_42 - mylar_thick - air_thick);
      remainder_42 = 0.0;
    }

    G4Box *segment_42 = new G4Box( "segment_42", (width_42 - 2.*(mylar_thick+air_thick) )/2.0, (width_42 - 2.*(mylar_thick+air_thick) )/2.0, segthick/2.0 );
    G4Box *lastsegment_42 = new G4Box( "lastsegment_42", (width_42 - 2.*(mylar_thick+air_thick) )/2.0, (width_42 - 2.*(mylar_thick+air_thick) )/2.0, (segthick+remainder_42)/2.0 );

    G4double ztemp = -depth_42/2.0 + mylar_thick + air_thick;
    
    for( int seg42=0; seg42<nseg_42; seg42++ ){
      TString material_name, lv_name, pv_name;
      material_name.Form( "TF1_anneal_ECAL_z%d", seg42 );
      lv_name.Form("segment_42_log_%d", seg42 );
      pv_name.Form("segment_42_phys_%d", seg42 );
      G4LogicalVolume *Seg42_log;
      if( seg42 + 1 < nseg_42 || nseg_42 == 1 ){
	Seg42_log = new G4LogicalVolume( segment_42, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += segthick/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg42_log, pv_name.Data(), Module_42_log, false, 0 );
	ztemp += segthick/2.0;
      } else {
	Seg42_log = new G4LogicalVolume( lastsegment_42, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += (segthick+remainder_42)/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg42_log, pv_name.Data(), Module_42_log, false, 0 );
	ztemp += (segthick+remainder_42)/2.0;
      }
      Seg42_log->SetSensitiveDetector( ECalTF1SD );

      G4VisAttributes *Segment_VisAtt = new G4VisAttributes( G4Colour( 0.8*(seg42/15.0)+0.20, 0.8*(seg42/15.0)+0.20, 0.0 ) );
      Seg42_log->SetVisAttributes( Segment_VisAtt );
    }   
  }
  /////End of Sensitivity coding for lead-glass





  //Place lead-glass and mylar wrap inside module:
  
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), hc_foil_42_log, "hc_foil_42_phys", Module_42_log, false, 0 );
  
  //mylar:
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), Mylar_wrap_42_log, "Mylar_wrap_42_phys", Module_42_log, false, 0 );

  new G4LogicalSkinSurface( "Mylar_skin_42", Mylar_wrap_42_log, GetOpticalSurface("Mirrsurf") );

  new G4LogicalSkinSurface( "lgwrap_42", LGWrap_42_log, GetOpticalSurface("Mirrsurf") );
  
  //Next: PMTs:
  G4Tubs *ecal_PMT = new G4Tubs( "ecal_PMT", 0.0, radius_ecal_pmt, depth_ecal_pmt/2.0, 0.0, twopi );
  G4LogicalVolume *ecal_PMT_log = new G4LogicalVolume( ecal_PMT, GetMaterial("Photocathode_material_ecal"), "ecal_PMT_log" );

  //sensitivity for pmts
  G4String sdname = "Earm/ECAL";
  G4String collname = "ECALHitsCollection";

  G4SBSECalSD *ECalSD = NULL;

  if( !( (G4SBSECalSD*) sdman->FindSensitiveDetector( sdname ) ) ){
    G4cout << "Adding ECAL PMT sensitive detector to SDman..." << G4endl;
    ECalSD = new G4SBSECalSD( sdname, collname );
    sdman->AddNewDetector( ECalSD );
    (fDetCon->SDlist).insert( sdname );
    fDetCon->SDtype[sdname] = G4SBS::kECAL;
    (ECalSD->detmap).depth = 0;
  }

  ecal_PMT_log->SetSensitiveDetector( ECalSD );
  fDetCon->InsertSDboundaryVolume( earm_mother_log->GetName(), sdname );
  //end of sensitivity for pmts

  
  //decided doing a subtraction for creating the surrounding inactive pb-glass would be better than dealing with overlaps due to placing the SM's inside of this new parent volume
  //create box of inactive pb-glass to surround ECal

  
  //using depth of SM such that none of the inactive pb-glass covers the front face of the pb-glass in ECal
  
  /*
  G4Box *inactivePbGls = new G4Box( "inactivePbGls", inactiveBox_width/2.0, inactiveBox_height/2.0, depth_42/2.0);
  G4LogicalVolume *inactivePbGls_log = new G4LogicalVolume( inactivePbGls, GetMaterial("TF1_dead"), "inactivePbGls_log" );
  new G4PVPlacement( 0, G4ThreeVector( 0 , 0 , zfront_ECAL + (depth_42/2.0)), inactivePbGls_log, "inactivePbGls_phys", earm_mother_log, false, 0 );
  */

  // Place the super modules 
   
  G4int NrowsSM_42 = 23; 
 
  G4int Nblocks_per_rowSM_42[23] = {4, 6, 7, 7, 8, 8, 
			  9, 9, 9, 9, 9, 9, 
			  9, 9, 9, 9, 9, 9, 
			  8, 8, 7, 7, 6};// from bottom to top

  // G4double yfp_start_40[10] = {-62.0*cm, -62.0*cm, -62.0*cm, -62.0*cm, -62.0*cm, -58.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm};//start of the block edge

  // G4double yfp_start_42[15] = {-54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -58.0*cm, -62.0*cm, -62.0*cm, -62.0*cm, -62.0*cm, -62.0*cm};//start of the block edge
  /*
  G4double yfp_start_42[23] = {-61.88*cm, -57.92*cm, -61.88*cm, -57.92*cm, -61.88*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, 
			       -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -54.0*cm, -61.79*cm,
			       -58.16*cm, -61.95*cm, -58.32*cm};//start of the block edge's right side, there are 14 starting at -54 from 12 rows of 9 SM's and 2 rows of 8 SM's alligned at their right side
  //this allignment is assuming the 12 rows of 9 SM's should start at -54 cm, this will likely need to be shifted
  */

  //shifted to be centered on beam line based on center in JT model
  G4double yfp_start_42[23] = {-58.69*cm, -54.73*cm, -58.69*cm, -54.73*cm, -58.69*cm, -50.81*cm, -50.81*cm, -50.81*cm, -50.81*cm, -50.81*cm, 
			       -50.81*cm, -50.81*cm, -50.81*cm, -50.81*cm, -50.81*cm, -50.81*cm, -50.81*cm, -50.81*cm, -50.81*cm, -58.60*cm,
			       -54.97*cm, -58.76*cm, -55.13*cm};// from bottom to top

  //xfp and yfp are following the analysis coordinate system, not the same as g4 coord system
  
  //xfpstart based on center of ECal JT model
  G4double xfpstart = -147.22*cm;
  G4int copy_nb = 0;
  G4double X_block, Y_block;

  //the width(or height) of a SM is defined elsewhere as 5.06in, however everything seems to be properly alligned only when using 5.04in
  G4double width42 = 5.134*2.54*cm;
  G4double height42 = 5.07*2.54*cm;
  G4double depthInactive = 15.75*2.54*cm;
  G4double inactiveBox_width = 54.42*2.54*cm;
  G4double inactiveBox_height = 125.12*2.54*cm;

  G4Box *inactivePbGls = new G4Box( "inactivePbGls", inactiveBox_width/2.0, width42/2.0, depthInactive/2.0 );

  G4LogicalVolume *inactivePbGls_log = new G4LogicalVolume( inactivePbGls, GetMaterial("TF1_dead"), "inactivePbGls_log" );
  
  //Make subtraction shell of inactive pb glass surrounding the Supermodules in ECal

  //this is extremely messy just because im currently not sure how to do this with a for loop, but trying to be methodical, numbering follows the order of the above arrays used to describe the rows of SM's going from top to bottom

  //width of # of supermodules to be subtracted from box of inactive pbglass, these are just the numbers needed from Nblocks_per_rowSM_42
  /* 
 G4double cutout4wide = 4*width42;
  G4double cutout6wide = 6*width42;
  G4double cutout7wide = 7*width42;
  G4double cutout8wide = 8*width42;
  G4double cutout9wide = 9*width42;

  */

  //z offset between front face of SM's and front face of inactive pbglass
  G4double zOffset42 = 0.5*2.54*cm;

  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[0] + (inactiveBox_width/2.0) , xfpstart - (1*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactivePbGls_log, "inactivePbGls_phys", earm_mother_log, false, 1 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[0] + (inactiveBox_width/2.0) , xfpstart + (47*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactivePbGls_log, "inactivePbGls_phys", earm_mother_log, false, 1 );
  
  //left side of inactive pb glass

  G4double widthGlsL0 = 29.45*2.54*cm;
  G4double widthGlsL1 = 16.46*2.54*cm;
  G4double widthGlsL2 = 12.92*2.54*cm;
  G4double widthGlsL3 = 11.43*2.54*cm;
  G4double widthGlsL4 = 8.61*2.54*cm;
  G4double widthGlsL5 = 5.47*2.54*cm;
  G4double widthGlsL6_17 = 1.0*2.54*cm;
  G4double widthGlsL18 = 5.77*2.54*cm;
  G4double widthGlsL19 = 8.75*2.54*cm;
  G4double widthGlsL20 = 13.29*2.54*cm;
  G4double widthGlsL21 = 15.1*2.54*cm;
  G4double widthGlsL22 = 18.76*2.54*cm;

  G4Box *inactiveGlsL0 = new G4Box( "inactiveGlsL0",( widthGlsL0 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsL1 = new G4Box( "inactiveGlsL1",( widthGlsL1 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsL2 = new G4Box( "inactiveGlsL2",( widthGlsL2 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsL3 = new G4Box( "inactiveGlsL3",( widthGlsL3 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsL4 = new G4Box( "inactiveGlsL4",( widthGlsL4 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsL5 = new G4Box( "inactiveGlsL5",( widthGlsL5 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsL6_17 = new G4Box( "inactiveGlsL6_17",( widthGlsL6_17 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsL18 = new G4Box( "inactiveGlsL18",( widthGlsL18 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsL19 = new G4Box( "inactiveGlsL19",( widthGlsL19 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsL20 = new G4Box( "inactiveGlsL20",( widthGlsL20 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsL21 = new G4Box( "inactiveGlsL21",( widthGlsL21 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsL22 = new G4Box( "inactiveGlsL22",( widthGlsL22 )/2.0, height42/2.0, depthInactive/2.0 );
  
  G4LogicalVolume *inactiveGlsL0_log = new G4LogicalVolume( inactiveGlsL0, GetMaterial("TF1_dead"), "inactiveGlsL0_log" );
  G4LogicalVolume *inactiveGlsL1_log = new G4LogicalVolume( inactiveGlsL1, GetMaterial("TF1_dead"), "inactiveGlsL1_log" );
  G4LogicalVolume *inactiveGlsL2_log = new G4LogicalVolume( inactiveGlsL2, GetMaterial("TF1_dead"), "inactiveGlsL2_log" );
  G4LogicalVolume *inactiveGlsL3_log = new G4LogicalVolume( inactiveGlsL3, GetMaterial("TF1_dead"), "inactiveGlsL3_log" );
  G4LogicalVolume *inactiveGlsL4_log = new G4LogicalVolume( inactiveGlsL4, GetMaterial("TF1_dead"), "inactiveGlsL4_log" );
  G4LogicalVolume *inactiveGlsL5_log = new G4LogicalVolume( inactiveGlsL5, GetMaterial("TF1_dead"), "inactiveGlsL5_log" );
  G4LogicalVolume *inactiveGlsL6_17_log = new G4LogicalVolume( inactiveGlsL6_17, GetMaterial("TF1_dead"), "inactiveGlsL6_17_log" );
  G4LogicalVolume *inactiveGlsL18_log = new G4LogicalVolume( inactiveGlsL18, GetMaterial("TF1_dead"), "inactiveGlsL18_log" );
  G4LogicalVolume *inactiveGlsL19_log = new G4LogicalVolume( inactiveGlsL19, GetMaterial("TF1_dead"), "inactiveGlsL19_log" );
  G4LogicalVolume *inactiveGlsL20_log = new G4LogicalVolume( inactiveGlsL20, GetMaterial("TF1_dead"), "inactiveGlsL20_log" );
  G4LogicalVolume *inactiveGlsL21_log = new G4LogicalVolume( inactiveGlsL21, GetMaterial("TF1_dead"), "inactiveGlsL21_log" );
  G4LogicalVolume *inactiveGlsL22_log = new G4LogicalVolume( inactiveGlsL22, GetMaterial("TF1_dead"), "inactiveGlsL22_log" );
  
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[0] + width42*Nblocks_per_rowSM_42[0] + widthGlsL0/2.0 , xfpstart + (1*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsL0_log, "inactiveGlsL0_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[1] + width42*Nblocks_per_rowSM_42[1] + widthGlsL1/2.0 , xfpstart + (3*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsL1_log, "inactiveGlsL1_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[2] + width42*Nblocks_per_rowSM_42[2] + widthGlsL2/2.0 , xfpstart + (5*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsL2_log, "inactiveGlsL2_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[3] + width42*Nblocks_per_rowSM_42[3] + widthGlsL3/2.0 , xfpstart + (7*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsL3_log, "inactiveGlsL3_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[4] + width42*Nblocks_per_rowSM_42[4] + widthGlsL4/2.0 , xfpstart + (9*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsL4_log, "inactiveGlsL4_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[5] + width42*Nblocks_per_rowSM_42[5] + widthGlsL5/2.0 , xfpstart + (11*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsL5_log, "inactiveGlsL5_phys", earm_mother_log, false, 0 );

  for( int i = 6; i < 18; i++ ){

    new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[i] + width42*Nblocks_per_rowSM_42[i] + widthGlsL6_17/2.0 , xfpstart + (((i*2)+1)*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsL6_17_log, "inactiveGlsL6_17_phys", earm_mother_log, false, 12 );

  }
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[18] + width42*Nblocks_per_rowSM_42[18] + widthGlsL18/2.0 , xfpstart + (37*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsL18_log, "inactiveGlsL18_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[19] + width42*Nblocks_per_rowSM_42[19] + widthGlsL19/2.0 , xfpstart + (39*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsL19_log, "inactiveGlsL19_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[20] + width42*Nblocks_per_rowSM_42[20] + widthGlsL20/2.0 , xfpstart + (41*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsL20_log, "inactiveGlsL20_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[21] + width42*Nblocks_per_rowSM_42[21] + widthGlsL21/2.0 , xfpstart + (43*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsL21_log, "inactiveGlsL21_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[22] + width42*Nblocks_per_rowSM_42[22] + widthGlsL22/2.0 , xfpstart + (45*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsL22_log, "inactiveGlsL22_phys", earm_mother_log, false, 0 ); 

  //right side of inactive pb glass

  G4double widthGlsR0 = 1.0*2.54*cm;
  G4double widthGlsR1 = 2.56*2.54*cm;
  G4double widthGlsR2 = 1.0*2.54*cm;
  G4double widthGlsR3 = 2.56*2.54*cm;
  G4double widthGlsR4 = 1.0*2.54*cm;
  G4double widthGlsR5 = 4.13*2.54*cm;
  G4double widthGlsR6_17 = 4.13*2.54*cm;
  G4double widthGlsR18 = 4.13*2.54*cm;
  G4double widthGlsR19 = 1.0*2.54*cm;
  G4double widthGlsR20 = 2.56*2.54*cm;
  G4double widthGlsR21 = 1.0*2.54*cm;
  G4double widthGlsR22 = 2.56*2.54*cm;

  G4Box *inactiveGlsR0 = new G4Box( "inactiveGlsR0",( widthGlsR0 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsR1 = new G4Box( "inactiveGlsR1",( widthGlsR1 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsR2 = new G4Box( "inactiveGlsR2",( widthGlsR2 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsR3 = new G4Box( "inactiveGlsR3",( widthGlsR3 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsR4 = new G4Box( "inactiveGlsR4",( widthGlsR4 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsR5 = new G4Box( "inactiveGlsR5",( widthGlsR5 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsR6_17 = new G4Box( "inactiveGlsR6_17",( widthGlsR6_17 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsR18 = new G4Box( "inactiveGlsR18",( widthGlsR18 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsR19 = new G4Box( "inactiveGlsR19",( widthGlsR19 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsR20 = new G4Box( "inactiveGlsR20",( widthGlsR20 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsR21 = new G4Box( "inactiveGlsR21",( widthGlsR21 )/2.0, height42/2.0, depthInactive/2.0 );
  G4Box *inactiveGlsR22 = new G4Box( "inactiveGlsR22",( widthGlsR22 )/2.0, height42/2.0, depthInactive/2.0 );
  
  G4LogicalVolume *inactiveGlsR0_log = new G4LogicalVolume( inactiveGlsR0, GetMaterial("TF1_dead"), "inactiveGlsR0_log" );
  G4LogicalVolume *inactiveGlsR1_log = new G4LogicalVolume( inactiveGlsR1, GetMaterial("TF1_dead"), "inactiveGlsR1_log" );
  G4LogicalVolume *inactiveGlsR2_log = new G4LogicalVolume( inactiveGlsR2, GetMaterial("TF1_dead"), "inactiveGlsR2_log" );
  G4LogicalVolume *inactiveGlsR3_log = new G4LogicalVolume( inactiveGlsR3, GetMaterial("TF1_dead"), "inactiveGlsR3_log" );
  G4LogicalVolume *inactiveGlsR4_log = new G4LogicalVolume( inactiveGlsR4, GetMaterial("TF1_dead"), "inactiveGlsR4_log" );
  G4LogicalVolume *inactiveGlsR5_log = new G4LogicalVolume( inactiveGlsR5, GetMaterial("TF1_dead"), "inactiveGlsR5_log" );
  G4LogicalVolume *inactiveGlsR6_17_log = new G4LogicalVolume( inactiveGlsR6_17, GetMaterial("TF1_dead"), "inactiveGlsR6_17_log" );
  G4LogicalVolume *inactiveGlsR18_log = new G4LogicalVolume( inactiveGlsR18, GetMaterial("TF1_dead"), "inactiveGlsR18_log" );
  G4LogicalVolume *inactiveGlsR19_log = new G4LogicalVolume( inactiveGlsR19, GetMaterial("TF1_dead"), "inactiveGlsR19_log" );
  G4LogicalVolume *inactiveGlsR20_log = new G4LogicalVolume( inactiveGlsR20, GetMaterial("TF1_dead"), "inactiveGlsR20_log" );
  G4LogicalVolume *inactiveGlsR21_log = new G4LogicalVolume( inactiveGlsR21, GetMaterial("TF1_dead"), "inactiveGlsR21_log" );
  G4LogicalVolume *inactiveGlsR22_log = new G4LogicalVolume( inactiveGlsR22, GetMaterial("TF1_dead"), "inactiveGlsR22_log" );

  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[0] - widthGlsR0/2.0 , xfpstart + (1*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsR0_log, "inactiveGlsR0_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[1] - widthGlsR1/2.0 , xfpstart + (3*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsR1_log, "inactiveGlsR1_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[2] - widthGlsR2/2.0 , xfpstart + (5*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsR2_log, "inactiveGlsR2_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[3] - widthGlsR3/2.0 , xfpstart + (7*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsR3_log, "inactiveGlsR3_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[4] - widthGlsR4/2.0 , xfpstart + (9*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsR4_log, "inactiveGlsR4_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[5] - widthGlsR5/2.0 , xfpstart + (11*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsR5_log, "inactiveGlsR5_phys", earm_mother_log, false, 0 );

  for( int i = 6; i < 18; i++ ){

    new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[i] - widthGlsR6_17/2.0 , xfpstart + (((i*2)+1)*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsR6_17_log, "inactiveGlsR6_17_phys", earm_mother_log, false, 12 );

  }
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[18] - widthGlsR18/2.0 , xfpstart + (37*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsR18_log, "inactiveGlsR18_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[19] - widthGlsR19/2.0 , xfpstart + (39*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsR19_log, "inactiveGlsR19_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[20] - widthGlsR20/2.0 , xfpstart + (41*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsR20_log, "inactiveGlsR20_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[21] - widthGlsR21/2.0 , xfpstart + (43*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsR21_log, "inactiveGlsR21_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( yfp_start_42[22] - widthGlsR22/2.0 , xfpstart + (45*height42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveGlsR22_log, "inactiveGlsR22_phys", earm_mother_log, false, 0 ); 

  G4VisAttributes *tf1dead_visatt = new G4VisAttributes( G4Colour( 0, 0.6, 0.5 ) );

  inactiveGlsL0_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsL1_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsL2_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsL3_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsL4_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsL5_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsL6_17_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsL18_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsL19_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsL20_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsL21_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsL22_log->SetVisAttributes(tf1dead_visatt);

  inactiveGlsR0_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsR1_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsR2_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsR3_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsR4_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsR5_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsR6_17_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsR18_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsR19_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsR20_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsR21_log->SetVisAttributes(tf1dead_visatt);
  inactiveGlsR22_log->SetVisAttributes(tf1dead_visatt);

  

  /*

  //boxes to be cut out, making them larger than the inactivePbGls to avoid issues of subtracting at coincident planes 
  G4Box *PbGlsCutout4wide = new G4Box( "PbGlsCutout4wide", cutout4wide/2.0, width42/2.0 + 1.0*cm, depthInactive/2.0 + 1.0*cm );
  G4Box *PbGlsCutout6wide = new G4Box( "PbGlsCutout6wide", cutout6wide/2.0, width42/2.0 + 1.0*cm, depthInactive/2.0 + 1.0*cm );
  G4Box *PbGlsCutout7wide = new G4Box( "PbGlsCutout7wide", cutout7wide/2.0, width42/2.0 + 1.0*cm, depthInactive/2.0 + 1.0*cm );
  G4Box *PbGlsCutout8wide = new G4Box( "PbGlsCutout8wide", cutout8wide/2.0, width42/2.0 + 1.0*cm, depthInactive/2.0 + 1.0*cm );
  G4Box *PbGlsCutout9wide = new G4Box( "PbGlsCutout9wide", cutout9wide/2.0, width42/2.0 + 1.0*cm, depthInactive/2.0 + 1.0*cm );

  //subtraction of cutouts from inactive pbglass and empty space centerd to match placement of SM rows within the mother volume
  //note:rows 6 through 17 have 9 SM's and are alligned
  G4SubtractionSolid *inactiveShell0 = new G4SubtractionSolid( "inactiveShell0", inactivePbGls, PbGlsCutout4wide, 0, G4ThreeVector( yfp_start_42[0] + (width42*Nblocks_per_rowSM_42[0]/2.0), 0, 0 ) );
  G4SubtractionSolid *inactiveShell1 = new G4SubtractionSolid( "inactiveShell1", inactivePbGls, PbGlsCutout6wide, 0, G4ThreeVector( yfp_start_42[1] + (width42*Nblocks_per_rowSM_42[1]/2.0), 0, 0 ) );
  G4SubtractionSolid *inactiveShell2 = new G4SubtractionSolid( "inactiveShell2", inactivePbGls, PbGlsCutout7wide, 0, G4ThreeVector( yfp_start_42[2] + (width42*Nblocks_per_rowSM_42[2]/2.0), 0, 0 ) );
  G4SubtractionSolid *inactiveShell3 = new G4SubtractionSolid( "inactiveShell3", inactivePbGls, PbGlsCutout7wide, 0, G4ThreeVector( yfp_start_42[3] + (width42*Nblocks_per_rowSM_42[3]/2.0), 0, 0 ) );
  G4SubtractionSolid *inactiveShell4 = new G4SubtractionSolid( "inactiveShell4", inactivePbGls, PbGlsCutout8wide, 0, G4ThreeVector( yfp_start_42[4] + (width42*Nblocks_per_rowSM_42[4]/2.0), 0, 0 ) );
  G4SubtractionSolid *inactiveShell5 = new G4SubtractionSolid( "inactiveShell5", inactivePbGls, PbGlsCutout8wide, 0, G4ThreeVector( yfp_start_42[5] + (width42*Nblocks_per_rowSM_42[5]/2.0), 0, 0 ) );
  G4SubtractionSolid *inactiveShell6_17 = new G4SubtractionSolid( "inactiveShell6_17", inactivePbGls, PbGlsCutout9wide, 0, G4ThreeVector( yfp_start_42[6] + (width42*Nblocks_per_rowSM_42[6]/2.0), 0, 0 ) );
  G4SubtractionSolid *inactiveShell18 = new G4SubtractionSolid( "inactiveShell18", inactivePbGls, PbGlsCutout8wide, 0, G4ThreeVector( yfp_start_42[18] + (width42*Nblocks_per_rowSM_42[18]/2.0), 0, 0 ) );
  G4SubtractionSolid *inactiveShell19 = new G4SubtractionSolid( "inactiveShell19", inactivePbGls, PbGlsCutout8wide, 0, G4ThreeVector( yfp_start_42[19] + (width42*Nblocks_per_rowSM_42[19]/2.0), 0, 0 ) );
  G4SubtractionSolid *inactiveShell20 = new G4SubtractionSolid( "inactiveShell20", inactivePbGls, PbGlsCutout7wide, 0, G4ThreeVector( yfp_start_42[20] + (width42*Nblocks_per_rowSM_42[20]/2.0), 0, 0 ) );
  G4SubtractionSolid *inactiveShell21 = new G4SubtractionSolid( "inactiveShell21", inactivePbGls, PbGlsCutout7wide, 0, G4ThreeVector( yfp_start_42[21] + (width42*Nblocks_per_rowSM_42[21]/2.0), 0, 0 ) );
  G4SubtractionSolid *inactiveShell22 = new G4SubtractionSolid( "inactiveShell22", inactivePbGls, PbGlsCutout6wide, 0, G4ThreeVector( yfp_start_42[22] + (width42*Nblocks_per_rowSM_42[22]/2.0), 0, 0 ) );
  
  G4LogicalVolume *inactiveShell0_log = new G4LogicalVolume( inactiveShell0, GetMaterial("TF1_dead"), "inactiveShell0_log" );
  G4LogicalVolume *inactiveShell1_log = new G4LogicalVolume( inactiveShell1, GetMaterial("TF1_dead"), "inactiveShell1_log" );
  G4LogicalVolume *inactiveShell2_log = new G4LogicalVolume( inactiveShell2, GetMaterial("TF1_dead"), "inactiveShell2_log" );
  G4LogicalVolume *inactiveShell3_log = new G4LogicalVolume( inactiveShell3, GetMaterial("TF1_dead"), "inactiveShell3_log" );
  G4LogicalVolume *inactiveShell4_log = new G4LogicalVolume( inactiveShell4, GetMaterial("TF1_dead"), "inactiveShell4_log" );
  G4LogicalVolume *inactiveShell5_log = new G4LogicalVolume( inactiveShell5, GetMaterial("TF1_dead"), "inactiveShell5_log" );
  G4LogicalVolume *inactiveShell6_17_log = new G4LogicalVolume( inactiveShell6_17, GetMaterial("TF1_dead"), "inactiveShell6_17_log" );
  G4LogicalVolume *inactiveShell18_log = new G4LogicalVolume( inactiveShell18, GetMaterial("TF1_dead"), "inactiveShell18_log" );
  G4LogicalVolume *inactiveShell19_log = new G4LogicalVolume( inactiveShell19, GetMaterial("TF1_dead"), "inactiveShell19_log" );
  G4LogicalVolume *inactiveShell20_log = new G4LogicalVolume( inactiveShell20, GetMaterial("TF1_dead"), "inactiveShell20_log" );
  G4LogicalVolume *inactiveShell21_log = new G4LogicalVolume( inactiveShell21, GetMaterial("TF1_dead"), "inactiveShell21_log" );
  G4LogicalVolume *inactiveShell22_log = new G4LogicalVolume( inactiveShell22, GetMaterial("TF1_dead"), "inactiveShell22_log" );

  //z offset between front face of SM's and front face of inactive pbglass
  G4double zOffset42 = 0.5*2.54*cm;

  //starting at bottom and going up 
  new G4PVPlacement( 0, G4ThreeVector( 0 , xfpstart - (1*width42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactivePbGls_log, "inactivePbGls_phys", earm_mother_log, false, 1 );
  new G4PVPlacement( 0, G4ThreeVector( 0 , xfpstart + (1*width42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveShell0_log, "inactiveShell0_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0 , xfpstart + (3*width42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveShell1_log, "inactiveShell1_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0 , xfpstart + (5*width42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveShell2_log, "inactiveShell2_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0 , xfpstart + (7*width42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveShell3_log, "inactiveShell3_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0 , xfpstart + (9*width42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveShell4_log, "inactiveShell4_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0 , xfpstart + (11*width42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveShell5_log, "inactiveShell5_phys", earm_mother_log, false, 0 );

  for( int i = 6; i < 18; i++ ){
    new G4PVPlacement( 0, G4ThreeVector( 0 , xfpstart + ((1+(i*2))*width42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveShell6_17_log, "inactiveShell6_17_phys", earm_mother_log, false, 12 );
  }

  new G4PVPlacement( 0, G4ThreeVector( 0 , xfpstart + (37*width42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveShell18_log, "inactiveShell18_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0 , xfpstart + (39*width42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveShell19_log, "inactiveShell19_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0 , xfpstart + (41*width42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveShell20_log, "inactiveShell20_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0 , xfpstart + (43*width42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveShell21_log, "inactiveShell21_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0 , xfpstart + (45*width42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactiveShell22_log, "inactiveShell22_phys", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0 , xfpstart + (47*width42/2.0) , zfront_ECAL + (depthInactive/2.0) - zOffset42), inactivePbGls_log, "inactivePbGls_phys", earm_mother_log, false, 1 );

  //color the inactive pbglass
  G4VisAttributes *tf1dead_visatt = new G4VisAttributes( G4Colour( 0, 0.6, 0.5 ) );

  inactivePbGls_log->SetVisAttributes(tf1dead_visatt);
  inactiveShell0_log->SetVisAttributes(tf1dead_visatt);
  inactiveShell1_log->SetVisAttributes(tf1dead_visatt);
  inactiveShell2_log->SetVisAttributes(tf1dead_visatt);
  inactiveShell3_log->SetVisAttributes(tf1dead_visatt);
  inactiveShell4_log->SetVisAttributes(tf1dead_visatt);
  inactiveShell5_log->SetVisAttributes(tf1dead_visatt);
  inactiveShell6_17_log->SetVisAttributes(tf1dead_visatt);
  inactiveShell18_log->SetVisAttributes(tf1dead_visatt);
  inactiveShell19_log->SetVisAttributes(tf1dead_visatt);
  inactiveShell20_log->SetVisAttributes(tf1dead_visatt);
  inactiveShell21_log->SetVisAttributes(tf1dead_visatt);
  inactiveShell22_log->SetVisAttributes(tf1dead_visatt);

  */

  // Make Super Modules for 42 mm and 40 mm wide blocks
  // G4double BlockSpace_42 = 1.68*inch;
  // G4double BlockFirst_42 = 0.84*inch;
  G4double BlockSpace_42 = 1.69*inch;
  G4double BlockFirst_42 = 0.845*inch;
  
  // G4double SMWidth_42_1 = 5.06*inch;
  G4double SMWidth_42_1 = 5.07*inch;
  G4double SMWidth_42_2 = 5.12*inch;
  G4double SMHeight_42 = 5.07*inch;

  G4double TiWallLength_42_1 = 14.83*inch;
  G4double TiWallLength_42_2 = 15.83*inch;
  G4double TiWallThick = 0.032*inch;

  G4double ClampingBarThick = 0.125*inch;
  G4double ClampingBarWidth = 0.78*inch;

  G4double SpacerLength_42_1 = 0.43*inch;
  G4double SpacerLength_42_2 = 1.43*inch;
  G4double SpacerOutDiam = 1.5*inch;
  
  G4double FrontFlangeThick = 0.25*inch;
  G4double BackFlangeThick = 0.5*inch;
  
  G4double BackFlangeHoleSize = 1.37*inch;
  G4double PMTFlangeHoleSize = 1.42*inch;
  
  G4double StandoffLength_42 = 3.00*inch;
  G4double StandoffDiam = 0.375*inch;

  //using only 42 for quick testing, should change instances of 40 elsewhere to 42 in order to get rid of unneccesary definitions but this is faster for now
  G4double BlockSpace_40 = BlockSpace_42;
  G4double BlockFirst_40 = BlockFirst_42;
  
  G4double SMWidth_40_1 = SMWidth_42_1;
  G4double SMWidth_40_2 = SMWidth_42_2;
  G4double SMHeight_40 = SMHeight_42;

  G4double TiWallLength_40_1 = TiWallLength_42_1;
  G4double TiWallLength_40_2 = TiWallLength_42_2;

  G4double SpacerLength_40_1 = SpacerLength_42_1;
  G4double SpacerLength_40_2 = SpacerLength_42_2;
  
  G4double StandoffLength_40 = StandoffLength_42;
  
  //Titanium Walls
  G4Box* TiSideWall_42_1_box = 
    new G4Box("TiSideWall_42_1_box", TiWallThick/2.0, SMHeight_42/2.0, TiWallLength_42_1/2.0);
  G4LogicalVolume* TiSideWall_42_1_log = 
    new G4LogicalVolume(TiSideWall_42_1_box, GetMaterial("Titanium"), "TiSideWall_42_1_log");

  G4Box* TiSideWall_42_2_box = 
    new G4Box("TiSideWall_42_2_box", TiWallThick/2.0, SMHeight_42/2.0, TiWallLength_42_2/2.0);
  G4LogicalVolume* TiSideWall_42_2_log = 
    new G4LogicalVolume(TiSideWall_42_2_box, GetMaterial("Titanium"), "TiSideWall_42_2_log");
  /* 
  G4Box* TiSideWall_40_1_box = 
    new G4Box("TiSideWall_40_1_box", TiWallThick/2.0, SMHeight_40/2.0, TiWallLength_40_1/2.0);
  G4LogicalVolume* TiSideWall_40_1_log = 
    new G4LogicalVolume(TiSideWall_40_1_box, GetMaterial("Titanium"), "TiSideWall_40_1_log");

  G4Box* TiSideWall_40_2_box = 
    new G4Box("TiSideWall_40_2_box", TiWallThick/2.0, SMHeight_40/2.0, TiWallLength_40_2/2.0);
  G4LogicalVolume* TiSideWall_40_2_log = 
    new G4LogicalVolume(TiSideWall_40_2_box, GetMaterial("Titanium"), "TiSideWall_40_2_log");
  */
  // Flanges:
  // Front
  G4Box* FrontFlange_42_1_box = 
    new G4Box("FrontFlange_42_1_box", SMWidth_42_1/2.0, SMHeight_42/2.0, FrontFlangeThick/2.0);
  G4Box* FrontFlange_42_2_box = 
    new G4Box("FrontFlange_42_2_box", SMWidth_42_2/2.0, SMHeight_42/2.0, FrontFlangeThick/2.0);
  
  G4Box* ClampingBar_box = 
    new G4Box("ClampingBar_box", ClampingBarWidth/2.0, ClampingBarWidth/2.0, ClampingBarThick/2.0);
  
  G4UnionSolid* FrontFlange_42_1_solid = 
    new G4UnionSolid("FrontFlange_42_1_solid", FrontFlange_42_1_box, ClampingBar_box, 0, 
		     G4ThreeVector(0.0, 0.0, (FrontFlangeThick+ClampingBarThick)/2.0));
  
  G4UnionSolid* FrontFlange_42_2_solid = 
    new G4UnionSolid("FrontFlange_42_2_solid", FrontFlange_42_2_box, ClampingBar_box, 0, 
		     G4ThreeVector(0.0, 0.0, (FrontFlangeThick+ClampingBarThick)/2.0));
  /*
  G4Box* FrontFlange_40_1_box = 
    new G4Box("FrontFlange_40_1_box", SMWidth_40_1/2.0, SMHeight_40/2.0, FrontFlangeThick/2.0);
  G4Box* FrontFlange_40_2_box = 
    new G4Box("FrontFlange_40_2_box", SMWidth_40_2/2.0, SMHeight_40/2.0, FrontFlangeThick/2.0);
  
  G4UnionSolid* FrontFlange_40_1_solid = 
    new G4UnionSolid("FrontFlange_40_1_solid", FrontFlange_40_1_box, ClampingBar_box, 0, 
		     G4ThreeVector(0.0, 0.0, (FrontFlangeThick+ClampingBarThick)/2.0));
  
  G4UnionSolid* FrontFlange_40_2_solid = 
    new G4UnionSolid("FrontFlange_40_2_solid", FrontFlange_40_2_box, ClampingBar_box, 0, 
		     G4ThreeVector(0.0, 0.0, (FrontFlangeThick+ClampingBarThick)/2.0));
  */
  // Back
  G4Box* BackFlange_42_1_box = 
    new G4Box("BackFlange_42_1_box", SMWidth_42_1/2.0, SMHeight_42/2.0, BackFlangeThick/2.0);
  
  G4Box* BackFlange_42_2_box = 
    new G4Box("BackFlange_42_2_box", SMWidth_42_2/2.0, SMHeight_42/2.0, BackFlangeThick/2.0);
  
  G4Tubs *BackFlange_hole = new G4Tubs("BackFlange_hole", 0.0, BackFlangeHoleSize/2.0, BackFlangeThick, 0.0*deg, 360.0*deg );
  
  
  G4SubtractionSolid* BackFlange_42_1_solid = 
    new G4SubtractionSolid("BackFlange_42_1_solid", BackFlange_42_1_box, BackFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  
  G4SubtractionSolid* BackFlange_42_2_solid = 
    new G4SubtractionSolid("BackFlange_42_2_solid", BackFlange_42_2_box, BackFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  /*
  G4Box* BackFlange_40_1_box = 
    new G4Box("BackFlange_40_1_box", SMWidth_40_1/2.0, SMHeight_40/2.0, BackFlangeThick/2.0);
  
  G4Box* BackFlange_40_2_box = 
    new G4Box("BackFlange_40_2_box", SMWidth_40_2/2.0, SMHeight_40/2.0, BackFlangeThick/2.0);
  
  G4SubtractionSolid* BackFlange_40_1_solid = 
    new G4SubtractionSolid("BackFlange_40_1_solid", BackFlange_40_1_box, BackFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  
  G4SubtractionSolid* BackFlange_40_2_solid = 
    new G4SubtractionSolid("BackFlange_40_2_solid", BackFlange_40_2_box, BackFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  */
  // PMT
  G4Tubs *PMTFlange_hole = new G4Tubs("PMTFlange_hole", 0.0, PMTFlangeHoleSize/2.0, FrontFlangeThick, 0.0*deg, 360.0*deg );
  
  G4SubtractionSolid* PMTFlange_42_1_solid = 
    new G4SubtractionSolid("PMTFlange_42_1_solid", FrontFlange_42_1_box, PMTFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  
  G4SubtractionSolid* PMTFlange_42_2_solid = 
    new G4SubtractionSolid("PMTFlange_42_2_solid", FrontFlange_42_2_box, PMTFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  /*
  G4SubtractionSolid* PMTFlange_40_1_solid = 
    new G4SubtractionSolid("PMTFlange_40_1_solid", FrontFlange_40_1_box, PMTFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  
  G4SubtractionSolid* PMTFlange_40_2_solid = 
    new G4SubtractionSolid("PMTFlange_40_2_solid", FrontFlange_40_2_box, PMTFlange_hole, 0, G4ThreeVector(0.0, 0.0, 0.0));
  */
  // Flanges boolean solids
  for(int i = 0; i<3; i++){
    for(int j = 0; j<3; j++){
      FrontFlange_42_1_solid = 
	new G4UnionSolid("FrontFlange_42_1_solid", FrontFlange_42_1_solid, ClampingBar_box, 0, 
			 G4ThreeVector(BlockSpace_42*(i-1),BlockSpace_42*(j-1), (FrontFlangeThick+ClampingBarThick)/2.0));
      FrontFlange_42_2_solid = 
	new G4UnionSolid("FrontFlange_42_1_solid", FrontFlange_42_2_solid, ClampingBar_box, 0, 
			 G4ThreeVector(BlockSpace_42*(i-1),BlockSpace_42*(j-1), (FrontFlangeThick+ClampingBarThick)/2.0));
      
      
       BackFlange_42_1_solid = new G4SubtractionSolid("BackFlange_42_1_solid", BackFlange_42_1_solid, BackFlange_hole, 0, G4ThreeVector(BlockSpace_42*(i-1),BlockSpace_42*(j-1),0.0));
      
       BackFlange_42_2_solid = new G4SubtractionSolid("BackFlange_42_2_solid", BackFlange_42_2_solid, BackFlange_hole, 0, G4ThreeVector(BlockSpace_42*(i-1),BlockSpace_42*(j-1),0.0));
      
      PMTFlange_42_1_solid = new G4SubtractionSolid("PMTFlange_42_1_solid", PMTFlange_42_1_solid, PMTFlange_hole, 0, G4ThreeVector(BlockSpace_42*(i-1),BlockSpace_42*(j-1),0.0));
      
      PMTFlange_42_2_solid = new G4SubtractionSolid("PMTFlange_42_2_solid", PMTFlange_42_2_solid, PMTFlange_hole, 0, G4ThreeVector(BlockSpace_42*(i-1),BlockSpace_42*(j-1),0.0));
      /*
      FrontFlange_40_1_solid = 
	new G4UnionSolid("FrontFlange_40_1_solid", FrontFlange_40_1_solid, ClampingBar_box, 0, 
			 G4ThreeVector(BlockSpace_40*(i-1),BlockSpace_40*(j-1), (FrontFlangeThick+ClampingBarThick)/2.0));
      FrontFlange_40_2_solid = 
	new G4UnionSolid("FrontFlange_40_1_solid", FrontFlange_40_2_solid, ClampingBar_box, 0, 
			 G4ThreeVector(BlockSpace_40*(i-1),BlockSpace_40*(j-1), (FrontFlangeThick+ClampingBarThick)/2.0));
      
      
       BackFlange_40_1_solid = new G4SubtractionSolid("BackFlange_40_1_solid", BackFlange_40_1_solid, BackFlange_hole, 0, G4ThreeVector(BlockSpace_40*(i-1),BlockSpace_40*(j-1),0.0));
      
       BackFlange_40_2_solid = new G4SubtractionSolid("BackFlange_40_2_solid", BackFlange_40_2_solid, BackFlange_hole, 0, G4ThreeVector(BlockSpace_40*(i-1),BlockSpace_40*(j-1),0.0));
      
      PMTFlange_40_1_solid = new G4SubtractionSolid("PMTFlange_40_1_solid", PMTFlange_40_1_solid, PMTFlange_hole, 0, G4ThreeVector(BlockSpace_40*(i-1),BlockSpace_40*(j-1),0.0));
      
      PMTFlange_40_2_solid = new G4SubtractionSolid("PMTFlange_40_2_solid", PMTFlange_40_2_solid, PMTFlange_hole, 0, G4ThreeVector(BlockSpace_40*(i-1),BlockSpace_40*(j-1),0.0));
      */
      }
  }
  
  G4LogicalVolume* FrontFlange_42_1_log = 
    new G4LogicalVolume(FrontFlange_42_1_solid, GetMaterial("Aluminum"), "FrontFlange_42_1_log");
  G4LogicalVolume* FrontFlange_42_2_log = 
    new G4LogicalVolume(FrontFlange_42_2_solid, GetMaterial("Aluminum"), "FrontFlange_42_2_log");
  
  G4LogicalVolume* BackFlange_42_1_log = 
    new G4LogicalVolume(BackFlange_42_1_solid, GetMaterial("Aluminum"), "BackFlange_42_1_log");
  G4LogicalVolume* BackFlange_42_2_log = 
    new G4LogicalVolume(BackFlange_42_2_solid, GetMaterial("Aluminum"), "BackFlange_42_2_log");
  
  G4LogicalVolume* PMTFlange_42_1_log = 
    new G4LogicalVolume(PMTFlange_42_1_solid, GetMaterial("Aluminum"), "PMTFlange_42_1_log");
  G4LogicalVolume* PMTFlange_42_2_log = 
    new G4LogicalVolume(PMTFlange_42_2_solid, GetMaterial("Aluminum"), "PMTFlange_42_2_log");
  /*
  G4LogicalVolume* FrontFlange_40_1_log = 
    new G4LogicalVolume(FrontFlange_40_1_solid, GetMaterial("Aluminum"), "FrontFlange_40_1_log");
  G4LogicalVolume* FrontFlange_40_2_log = 
    new G4LogicalVolume(FrontFlange_40_2_solid, GetMaterial("Aluminum"), "FrontFlange_40_2_log");
  
  G4LogicalVolume* BackFlange_40_1_log = 
    new G4LogicalVolume(BackFlange_40_1_solid, GetMaterial("Aluminum"), "BackFlange_40_1_log");
  G4LogicalVolume* BackFlange_40_2_log = 
    new G4LogicalVolume(BackFlange_40_2_solid, GetMaterial("Aluminum"), "BackFlange_40_2_log");
  
  G4LogicalVolume* PMTFlange_40_1_log = 
    new G4LogicalVolume(PMTFlange_40_1_solid, GetMaterial("Aluminum"), "PMTFlange_40_1_log");
  G4LogicalVolume* PMTFlange_40_2_log = 
    new G4LogicalVolume(PMTFlange_40_2_solid, GetMaterial("Aluminum"), "PMTFlange_40_2_log");
  */
  // Spacers
  G4Tubs *Spacer_42_1_solid = new G4Tubs("Spacer_42_1_solid", BackFlangeHoleSize/2.0, SpacerOutDiam/2.0, SpacerLength_42_1/2.0, 0.0*deg, 360.0*deg );
  G4Tubs *Spacer_42_2_solid = new G4Tubs("Spacer_42_2_solid", BackFlangeHoleSize/2.0, SpacerOutDiam/2.0, SpacerLength_42_2/2.0, 0.0*deg, 360.0*deg );
   
  G4LogicalVolume* Spacer_42_1_log = 
    new G4LogicalVolume(Spacer_42_1_solid, 
			GetMaterial("Aluminum"), "Spacer_42_1_log");
  G4LogicalVolume* Spacer_42_2_log = 
    new G4LogicalVolume(Spacer_42_2_solid, 
			GetMaterial("Aluminum"), "Spacer_42_2_log");
  /*
  G4Tubs *Spacer_40_1_solid = new G4Tubs("Spacer_40_1_solid", BackFlangeHoleSize/2.0, SpacerOutDiam/2.0, SpacerLength_40_1/2.0, 0.0*deg, 360.0*deg );
  G4Tubs *Spacer_40_2_solid = new G4Tubs("Spacer_40_2_solid", BackFlangeHoleSize/2.0, SpacerOutDiam/2.0, SpacerLength_40_2/2.0, 0.0*deg, 360.0*deg );
   
  G4LogicalVolume* Spacer_40_1_log = 
    new G4LogicalVolume(Spacer_40_1_solid, 
			GetMaterial("Aluminum"), "Spacer_40_1_log");
  G4LogicalVolume* Spacer_40_2_log = 
    new G4LogicalVolume(Spacer_40_2_solid, 
			GetMaterial("Aluminum"), "Spacer_40_2_log");
  */
  //Standoffs
  G4Tubs *Standoff_42_solid = new G4Tubs("Standoff_42_solid", 0, StandoffDiam/2.0, StandoffLength_42/2.0, 0.0*deg, 360.0*deg );
  
  G4LogicalVolume* Standoff_42_log = 
    new G4LogicalVolume(Standoff_42_solid, 
			GetMaterial("Aluminum"), "Standoff_42_log");
  /*
  G4Tubs *Standoff_40_solid = new G4Tubs("Standoff_40_solid", 0, StandoffDiam/2.0, StandoffLength_40/2.0, 0.0*deg, 360.0*deg );
  
  G4LogicalVolume* Standoff_40_log = 
    new G4LogicalVolume(Standoff_40_solid, 
			GetMaterial("Aluminum"), "Standoff_40_log");
  */
  ofstream mapfile("database/ecal_gep_blockmap.txt");

  TString  output_line = "";
  output_line.Form("#%16s, %16s, %16s, %16s, %16s", "Cell", "Row", "Column", "Xcenter (cm)", "Ycenter (cm)" ); 

  mapfile << output_line << endl;
  

  //replace all 40's with 42 here, original 42's placement currently commented out
   
  G4int SM_num = 0;
  X_block = xfpstart+BlockFirst_42;
  for(int i_ = 0; i_<NrowsSM_42*3; i_++){
    Y_block = yfp_start_42[i_/3]+TiWallThick+BlockFirst_42;
    for(int j_ = 0; j_<Nblocks_per_rowSM_42[i_/3]*3; j_++){
      //printf("i_ = %d, j = %d, i_/3, copy_nb = %d, X_block = %f, Y_Block = %f \n", i_, j_, i_/3, copy_nb, X_block, Y_block);
      G4ThreeVector modpos( Y_block, X_block, zfront_ECAL + depth_42/2.0 );

      // depth_earm/2.0 - depth_ecal_pmt - depth_leadglass/2.0 );
      new G4PVPlacement( 0, modpos, Module_42_log, "Module_42_phys", earm_mother_log, false, copy_nb );

      output_line.Form( "%16d, %16d, %16d, %16.3f, %16.3f", copy_nb, i_, j_, Y_block/cm, X_block/cm );
      mapfile << output_line << endl;
      
      (ECalTF1SD->detmap).Row[copy_nb] = i_;
      (ECalTF1SD->detmap).Col[copy_nb] = j_;
      (ECalTF1SD->detmap).LocalCoord[copy_nb] = modpos;
      
      G4ThreeVector pmtpos( modpos.x(), modpos.y(), 
			    modpos.z() + depth_42/2.0 + LG42->GetZHalfLength()*2 + depth_ecal_pmt/2.0 );
      new G4PVPlacement( 0, pmtpos, ecal_PMT_log, "ecal_PMT_phys", earm_mother_log, false, copy_nb );
      
      (ECalSD->detmap).Row[copy_nb] = i_;
      (ECalSD->detmap).Col[copy_nb] = j_;
      (ECalSD->detmap).LocalCoord[copy_nb] = pmtpos;
      
      //Add light-guide with mylar wrap:
      G4ThreeVector LGpos( modpos.x(), modpos.y(), modpos.z() + depth_42/2.0 + LG42->GetZHalfLength() );
      new G4PVPlacement( 0, LGpos, LG42_log, "LG42_phys", earm_mother_log, false, copy_nb );
      
      // //EFuchey 2017-01-12: Need to make sensitive the three volumes above, to measure their dose.
      // // shall be temporary, and not end in the repo...
      // (ECalLGSD->detmap).Row[copy_nb] = global_row;
      // (ECalLGSD->detmap).Col[copy_nb] = col;
      // (ECalLGSD->detmap).LocalCoord[copy_nb] = modpos;
      
      // Placing the supermodule
      // Spacers...
      G4double z0_SM = zfront_ECAL - FrontFlangeThick - ClampingBarThick;
      G4RotationMatrix* rm_temp = new G4RotationMatrix();
      G4ThreeVector pos_temp(Y_block, X_block, z0_SM + TiWallLength_42_1-BackFlangeThick-SpacerLength_42_1/2.0);
      new G4PVPlacement(rm_temp, pos_temp, Spacer_42_1_log, "Spacer_42_1_phys", earm_mother_log, false, 0 );
      
      if(i_%3==1 && j_%3==1){
	// Front Flange
	pos_temp.setZ(z0_SM + FrontFlangeThick/2.0);
	new G4PVPlacement(rm_temp, pos_temp, FrontFlange_42_1_log, "FrontFlange_42_1_phys", earm_mother_log, false, 0 );
	
	// Back Flange
	pos_temp.setZ(z0_SM + TiWallLength_42_1 - BackFlangeThick/2.0);
	new G4PVPlacement(rm_temp, pos_temp, BackFlange_42_1_log, "BackFlange_42_1_phys", earm_mother_log, false, 0 );
	
	// PMT Flange
	pos_temp.setZ(z0_SM + TiWallLength_42_1 + StandoffLength_42 + FrontFlangeThick/2.0);
	new G4PVPlacement(rm_temp, pos_temp, PMTFlange_42_1_log, "PMTFlange_42_1_phys", earm_mother_log, false, 0 );
	
	// Ti walls and Standoffs
	for(int k = 0; k<2; k++){
	  pos_temp.setZ(z0_SM + TiWallLength_42_1/2.0);
	  pos_temp.setY(X_block);
	  //pos_temp.setX(Y_block+(k-0.5)*(SMWidth_42_1+TiWallThick/2.0));
	  pos_temp.setX(Y_block+(k-0.5)*((3*BlockSpace_42)+(2*TiWallThick)-TiWallThick));
	  new G4PVPlacement(rm_temp, pos_temp, TiSideWall_42_1_log, "TiSideWall_42_1_phys", 
			    earm_mother_log, false, 0 );
	  for(int l = 0; l<2; l++){
	    pos_temp.setX(Y_block+BlockSpace_42*(l-0.5));
	    pos_temp.setY(X_block+(SMHeight_42-StandoffDiam)*(k-0.5));
	    pos_temp.setZ(z0_SM + TiWallLength_42_1 + StandoffLength_42/2.0);
	    new G4PVPlacement(rm_temp, pos_temp, Standoff_42_log, "Standoff_42_phys", 
			      earm_mother_log, false, 0 );
	  }
	}
	SM_num++;
	//if(NcolsSM_40[i_/3]%2==0 && j_>(NcolsSM_40[i_/3]-1)*3)SM_num++;
      }
      
      Y_block+= BlockSpace_42;
      // if(j_%3==2)Y_block+= 2*BlockFirst_42+2*TiWallThick-BlockSpace_42;
      if(j_%3==2)Y_block+= 2*TiWallThick;
      copy_nb++;
    }
    X_block+= BlockSpace_42;
  }
  
  /*

  //arranging and placing the 42's
  X_block+= BlockFirst_42+BlockFirst_40-BlockSpace_40;
  for(int i_ = 0; i_<NrowsSM_42*3; i_++){
    Y_block = yfp_start_42[i_/3]+TiWallThick+BlockFirst_42;
    for(int j_ = 0; j_<Nblocks_per_rowSM_42[i_/3]*3; j_++){
      //printf("i_ = %d, j = %d, i_/3, copy_nb = %d, X_block = %f, Y_Block = %f \n", i_, j_, i_/3, copy_nb, X_block, Y_block);
      G4ThreeVector modpos( Y_block, X_block, zfront_ECAL + depth_42/2.0 );
      // depth_earm/2.0 - depth_ecal_pmt - depth_leadglass/2.0 );
      new G4PVPlacement( 0, modpos, Module_42_log, "Module_42_phys", earm_mother_log, false, copy_nb );

      output_line.Form( "%16d, %16d, %16d, %16.3f, %16.3f", copy_nb, i_ + NrowsSM_40*3, j_, Y_block/cm, X_block/cm );
      mapfile << output_line << endl;
      
      (ECalTF1SD->detmap).Row[copy_nb] = i_ + 3*NrowsSM_40;
      (ECalTF1SD->detmap).Col[copy_nb] = j_;
      (ECalTF1SD->detmap).LocalCoord[copy_nb] = modpos;

      G4ThreeVector pmtpos( modpos.x(), modpos.y(),  
			    modpos.z() + depth_42/2.0 + LG42->GetZHalfLength()*2 + depth_ecal_pmt/2.0 );
      new G4PVPlacement( 0, pmtpos, ecal_PMT_log, "ecal_PMT_phys", earm_mother_log, false, copy_nb );

      (ECalSD->detmap).Row[copy_nb] = i_ + 3*NrowsSM_40;
      (ECalSD->detmap).Col[copy_nb] = j_;
      (ECalSD->detmap).LocalCoord[copy_nb] = pmtpos;

      //Add light-guide with mylar wrap:
      G4ThreeVector LGpos( modpos.x(), modpos.y(), modpos.z() + depth_42/2.0 + LG42->GetZHalfLength() );
      new G4PVPlacement( 0, LGpos, LG42_log, "LG42_phys", earm_mother_log, false, copy_nb );

      // //EFuchey 2017-01-12: Need to make sensitive the three volumes above, to measure their dose.
      // // shall be temporary, and not end in the repo...
      // (ECalLGSD->detmap).Row[copy_nb] = global_row;
      // (ECalLGSD->detmap).Col[copy_nb] = col;
      // (ECalLGSD->detmap).LocalCoord[copy_nb] = modpos;
      
      // Placing the supermodule
      // Spacers...
      G4double z0_SM = zfront_ECAL - FrontFlangeThick - ClampingBarThick;
      G4RotationMatrix* rm_temp = new G4RotationMatrix();
      G4ThreeVector pos_temp(Y_block, X_block, z0_SM + TiWallLength_42_1-BackFlangeThick-SpacerLength_42_1/2.0);
      new G4PVPlacement(rm_temp, pos_temp, Spacer_42_1_log, "Spacer_42_1_phys", earm_mother_log, false, 0 );
      
      if(i_%3==1 && j_%3==1){
	// Front Flange
	pos_temp.setZ(z0_SM + FrontFlangeThick/2.0);
	new G4PVPlacement(rm_temp, pos_temp, FrontFlange_42_1_log, "FrontFlange_42_1_phys", earm_mother_log, false, 0 );
	
	// Back Flange
	pos_temp.setZ(z0_SM + TiWallLength_42_1 - BackFlangeThick/2.0);
	new G4PVPlacement(rm_temp, pos_temp, BackFlange_42_1_log, "BackFlange_42_1_phys", earm_mother_log, false, 0 );
	
	// PMT Flange
	pos_temp.setZ(z0_SM + TiWallLength_42_1 + StandoffLength_42 + FrontFlangeThick/2.0);
	new G4PVPlacement(rm_temp, pos_temp, PMTFlange_42_1_log, "PMTFlange_42_1_phys", earm_mother_log, false, 0 );
	
	// Ti walls and Standoffs
	for(int k = 0; k<2; k++){
	  pos_temp.setZ(z0_SM + TiWallLength_42_1/2.0);
	  pos_temp.setY(X_block);
	  pos_temp.setX(Y_block+(k-0.5)*(SMWidth_42_1+TiWallThick/2.0));
	  new G4PVPlacement(rm_temp, pos_temp, TiSideWall_42_1_log, "TiSideWall_42_1_phys", 
			    earm_mother_log, false, 0 );
	  for(int l = 0; l<2; l++){
	    pos_temp.setX(Y_block+BlockSpace_42*(l-0.5));
	    pos_temp.setY(X_block+(SMHeight_42-StandoffDiam)*(k-0.5));
	    pos_temp.setZ(z0_SM + TiWallLength_42_1 + StandoffLength_42/2.0);
	    new G4PVPlacement(rm_temp, pos_temp, Standoff_42_log, "Standoff_42_phys", 
			      earm_mother_log, false, 0 );
	  }
	}
	SM_num++;
	//if(NcolsSM_42[i_/3]%2==0 && j_>(NcolsSM_42[i_/3]-1)*3)SM_num++;
      }
      Y_block+= BlockSpace_42;
      if(j_%3==2)Y_block+= 2*BlockFirst_42+2*TiWallThick-BlockSpace_42;
      copy_nb++;
    }
    X_block+= BlockSpace_42;
    if(i_==2)X_block+= 2.0*mm;// to leave space for the steel plate instered there
  }
  */
  mapfile.close();
  
  // Build CDet:
  if(fDetCon->fExpType==G4SBS::kGEp || fDetCon->fExpType == G4SBS::kGEPpositron ){
    //Next: CH2 filter:
    G4Box *CH2_filter = new G4Box( "CH2_filter", width_earm/2.0, height_earm/2.0, depth_CH2/2.0 );
    G4LogicalVolume *CH2_filter_log = new G4LogicalVolume( CH2_filter, GetMaterial("Polyethylene"), "CH2_filter_log" );
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, zfront_ECAL - depth_CDET - depth_CH2/2.0 ), CH2_filter_log, "CH2_filter_phys", earm_mother_log, false, 0 );
    
    G4double z0_CDET = -depth_earm/2.0 + depth_CH2;
    //G4double R0_CDET = R_Earm - depth_leadglass - depth_CDET;
    //G4double R0_CDET = fDist - depth_leadglass - depth_CDET;
    G4double R0_CDET = fDist - depth_CDET;
    
    G4SBSCDet* CDet = new G4SBSCDet(fDetCon);
    CDet->SetArmName("Earm");
    CDet->SetR0(R0_CDET);
    CDet->SetZ0(z0_CDET);
    CDet->SetPlanesHOffset(0.0);
    CDet->SetPlanesInterDistance(1.0*cm);
    CDet->BuildComponent( earm_mother_log );
    //MakeCDET( earm_mother_log );
    
    G4VisAttributes *CH2_visatt = new G4VisAttributes( G4Colour( 0, 0.6, 0.6 ) );
    CH2_visatt->SetForceWireframe(true);
    
    CH2_filter_log->SetVisAttributes(CH2_visatt);
  }

   //Create heat shielding in between ECal supermodule array and CDet

  G4double HSCent_width = 8.908*2.54*cm;
  G4double HSSide_width = 36.2*2.54*cm;
  G4double HSAll_height = 24.058*2.54*cm;
  G4double HSAll_depth = 9.944*2.54*cm;
  G4double HSAll_vert_space = 0.29*2.54*cm;
  G4double HSAll_hrz_space = 0.19*2.54*cm;
  G4double HSAll_ystart = 2.5*HSAll_height + 2.0*HSAll_vert_space;
  G4double HSSide_xblockcenter = HSCent_width/2.0 + HSAll_hrz_space + HSSide_width/2.0;
  
  G4Box *HeatShield_cent = new G4Box( "HeatShield_cent", HSCent_width/2.0, HSAll_height/2.0, HSAll_depth/2.0 );
  G4LogicalVolume *HeatShield_cent_log = new G4LogicalVolume( HeatShield_cent, GetMaterial("SiO2_C16"), "HeatShield_cent_log" );

  G4Box *HeatShield_side = new G4Box( "HeatShield_side", HSSide_width/2.0, HSAll_height/2.0, HSAll_depth/2.0 );
  G4LogicalVolume *HeatShield_side_log = new G4LogicalVolume( HeatShield_side, GetMaterial("SiO2_C16"), "HeatShield_side_log" );

  for(int i = 0; i < 5; i++){ 
    new G4PVPlacement( 0, G4ThreeVector( 0, HSAll_ystart - (HSAll_height/2.0) - i*(HSAll_height + HSAll_vert_space), zfront_ECAL - 2.23*2.54*cm - HSAll_depth/2.0 ), HeatShield_cent_log, "HeatShield_cent_phys", earm_mother_log, false, 5 );
    new G4PVPlacement( 0, G4ThreeVector( HSSide_xblockcenter, HSAll_ystart - (HSAll_height/2.0) - i*(HSAll_height + HSAll_vert_space), zfront_ECAL - 2.23*2.54*cm - HSAll_depth/2.0 ), HeatShield_side_log, "HeatShield_side_phys_left", earm_mother_log, false, 5 );
    new G4PVPlacement( 0, G4ThreeVector( -HSSide_xblockcenter, HSAll_ystart - (HSAll_height/2.0) - i*(HSAll_height + HSAll_vert_space), zfront_ECAL - 2.23*2.54*cm - HSAll_depth/2.0 ), HeatShield_side_log, "HeatShield_side_phys_right", earm_mother_log, false, 5 ); 
  }
  //color for heat shielding
  G4VisAttributes *HeatShield_visatt = new G4VisAttributes( G4Colour( 0.9, 0.7, 0.9 ) );
  HeatShield_cent_log->SetVisAttributes(HeatShield_visatt);
  HeatShield_side_log->SetVisAttributes(HeatShield_visatt);

  //Create absorbers in front of cdet

  G4double absorb_width = 48*2.54*cm;
  G4double absorb_height = 120*2.54*cm;
  G4double absorb_depth = 6*2.54*cm;

  G4Box *absorber = new G4Box( "absorber", absorb_width/2.0, absorb_height/2.0, absorb_depth/2.0 );
  G4LogicalVolume *absorber_log = new G4LogicalVolume( absorber, GetMaterial("Polyethylene"), "absorber_log" );
  new G4PVPlacement( 0, G4ThreeVector( 0 , 0 , zfront_ECAL - 28.45*2.54*cm - (absorb_depth/2.0)), absorber_log, "absorber_phys", earm_mother_log, false, 0); 

  //color for absorber
  G4VisAttributes *absorber_visatt = new G4VisAttributes( G4Colour( 0, 1, 1 ) );
  absorber_log->SetVisAttributes(absorber_visatt);

  //Visualization:
  G4VisAttributes *mother_visatt = new G4VisAttributes( G4Colour( 1, 1, 1 ) );
  mother_visatt->SetForceWireframe(true);
  earm_mother_log->SetVisAttributes(mother_visatt);
   
  //earm_mother_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  Module_42_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // Module_40_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  
  G4VisAttributes *Mylarvisatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
  Mylarvisatt->SetForceWireframe(true);
  Mylar_wrap_42_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  //Mylar_wrap_40_log->SetVisAttributes( G4VisAttributes::GetInvisible() );

  G4VisAttributes *Ti_visatt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0.7 ) );
  Ti_visatt->SetForceWireframe(true);
  TiSideWall_42_1_log->SetVisAttributes(Ti_visatt);
  TiSideWall_42_2_log->SetVisAttributes(Ti_visatt);
  //TiSideWall_40_1_log->SetVisAttributes(Ti_visatt);
  //TiSideWall_40_2_log->SetVisAttributes(Ti_visatt);
  // TiSideWall_log->SetVisAttributes(G4VisAttributes::GetInvisible());//
  
  G4VisAttributes *Al_visatt = new G4VisAttributes( G4Colour( 0.7, 0.7, 0.7 ) );
  Al_visatt->SetForceWireframe(true);
  FrontFlange_42_1_log->SetVisAttributes(Al_visatt);
  FrontFlange_42_2_log->SetVisAttributes(Al_visatt);
  //FrontFlange_40_1_log->SetVisAttributes(Al_visatt);
  //FrontFlange_40_2_log->SetVisAttributes(Al_visatt);
  BackFlange_42_1_log->SetVisAttributes(Al_visatt);
  BackFlange_42_2_log->SetVisAttributes(Al_visatt);
  //BackFlange_40_1_log->SetVisAttributes(Al_visatt);
  //BackFlange_40_2_log->SetVisAttributes(Al_visatt);
  PMTFlange_42_1_log->SetVisAttributes( G4Colour( 0.0, 0.7, 0.0 ) );
  PMTFlange_42_2_log->SetVisAttributes( G4Colour( 0.7, 0.0, 0.0 ) );
  //PMTFlange_40_1_log->SetVisAttributes( G4Colour( 0.0, 0.7, 0.0 ) );
  //PMTFlange_40_2_log->SetVisAttributes( G4Colour( 0.7, 0.0, 0.0 ) );
  Spacer_42_1_log->SetVisAttributes( G4VisAttributes::GetInvisible() ); //G4Colour( 0.0, 0.7, 0.0 ));
  Spacer_42_2_log->SetVisAttributes( G4VisAttributes::GetInvisible() ); //G4Colour( 0.0, 0.7, 0.0 ));
  //Spacer_40_1_log->SetVisAttributes( G4VisAttributes::GetInvisible() ); //G4Colour( 0.0, 0.7, 0.0 ));
  //Spacer_40_2_log->SetVisAttributes( G4VisAttributes::GetInvisible() ); //G4Colour( 0.0, 0.7, 0.0 ));
  Standoff_42_log->SetVisAttributes( G4VisAttributes::GetInvisible() ); //Al_visatt);
  //Standoff_40_log->SetVisAttributes( G4VisAttributes::GetInvisible() ); //Al_visatt);
  
  
  G4VisAttributes *ECALpmtvisatt = new G4VisAttributes( G4Colour( 0, 0, 1 ) );
  ecal_PMT_log->SetVisAttributes( ECALpmtvisatt );
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////KIP
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////KIP

void G4SBSECal::MakeC16( G4LogicalVolume *motherlog ){
  printf("C16 at %f deg\n", fAng/deg);

  G4SDManager *sdman = fDetCon->fSDman;

  // Module specs
  G4double width_42 = 4.2*cm;
  G4double depth_leadglass = 34.3*cm;
  G4int nrows = 4;
  G4int ncols = 4;

  // Wave-guide specs
  G4double depth_WG  = 15.0*cm;
  G4double radius_WG = 1.25*cm;
  
  // PMT specs
  G4double depth_ecal_pmt = 0.3*cm;
  G4double radius_ecal_pmt = 1.25*cm;

  // We should be using aluminized foil instead of mylar.
  G4double alum_thick = 0.001*2.54*cm;
  G4double air_thick = alum_thick;

  // Make a C16 Mother Volume to house all the modules+WG+PMTs:
  G4double C16_width  = nrows*width_42;
  G4double C16_depth  = depth_leadglass + depth_WG + depth_ecal_pmt;

  // Geometric information & C16 placing
  G4double r_C16 = fDist + C16_depth/2.0;      // Distance from target to front of C16
  G4double angle_C16 = fAng;                    // C16 BB angle
  G4ThreeVector R_C16( r_C16*sin( angle_C16 ), 0.0, r_C16*cos( angle_C16 ) );

  G4RotationMatrix *bbrm_C16 = new G4RotationMatrix;
  bbrm_C16->rotateY( -angle_C16 );

  G4Box *C16_Box = new G4Box( "C16_Box", C16_width/2.0, C16_width/2.0, C16_depth/2.0 );
  G4LogicalVolume *C16_Log = new G4LogicalVolume( C16_Box, GetMaterial("Special_Air"), "C16_Log" );
  new G4PVPlacement( bbrm_C16, R_C16, C16_Log, "C16_Phys", motherlog, false, 0 );

  // Getting the PMTs set up - use same PMT as ECal, same sensitivity
  G4Tubs *ecal_PMT = new G4Tubs( "ecal_PMT", 0.0, radius_ecal_pmt, depth_ecal_pmt/2.0, 0.0, twopi );
  G4LogicalVolume *ecal_PMT_log = new G4LogicalVolume( ecal_PMT, GetMaterial("Photocathode_material_ecal"), "ecal_PMT_log" );

  G4String C16SDname = "Earm/C16";
  G4String collname = "C16HitsCollection";
  G4SBSECalSD *C16SD = NULL;

  if( !((G4SBSECalSD*) sdman->FindSensitiveDetector( C16SDname )) ){
    G4cout << "Adding C16 PMT sensitive detector to SDman..." << G4endl;
    C16SD = new G4SBSECalSD( C16SDname, collname );
    sdman->AddNewDetector( C16SD );
    (fDetCon->SDlist).insert( C16SDname );
    fDetCon->SDtype[C16SDname] = G4SBS::kECAL;
    (C16SD->detmap).depth = 0;
  }

  ecal_PMT_log->SetSensitiveDetector( C16SD );

  fDetCon->InsertSDboundaryVolume( C16_Log->GetName(), C16SDname );
  
  // Getting the Wave Guides set up, *****need to UPDATE material properties*****
  G4Tubs *C16_WG = new G4Tubs( "C16_WG", 0.0, radius_WG, depth_WG/2.0, 0.0, twopi );
  G4LogicalVolume *C16_WG_Log = new G4LogicalVolume( C16_WG, GetMaterial("Pyrex_Glass"), "C16_WG_Log" );
  //new G4LogicalSkinSurface( "C16_WG_Skin", C16_WG_Log, GetOpticalSurface("osWLSToAir") );

  // Make a trapzoidal aluminum piece that starts at the end of the TF1 and ends
  // half way down the WG
  G4double  dx1 = width_42 / 2.0;
  G4double  dx2 = (radius_WG + 2*alum_thick + 2*air_thick);
  G4double  dy1 = width_42 / 2.0;
  G4double  dy2 = (radius_WG + 2*alum_thick + 2*air_thick);
  G4double  dz  = depth_WG / 4.0;
  G4Trd *Al_endpiece = new G4Trd( "Al_endpiece", dx1, dx2, dy1, dy2, dz );

  // Make a subtraction in order to get Trapezoidal Al foil
  dx1 = (width_42 - 2*alum_thick - 2*air_thick) / 2.0;
  dx2 = radius_WG;
  dy1 = (width_42 - 2*alum_thick - 2*air_thick) / 2.0;
  dy2 = radius_WG;
  dz = (depth_WG + 0.1*cm) / 4.0;
  G4Trd *Al_wrap_endpiece_sub = new G4Trd( "Al_wrap_endpiece_sub", dx1, dx2, dy1, dy2, dz );
  G4SubtractionSolid *Al_wrap_endpiece = new G4SubtractionSolid( "Al_wrap_endpiece", Al_endpiece, Al_wrap_endpiece_sub, 0, G4ThreeVector(0.0,0.0,0.0) );
  G4LogicalVolume *Al_wrap_endpiece_log = new G4LogicalVolume( Al_wrap_endpiece, GetMaterial("Aluminum"), "Al_wrap_endpiece_log" );
  new G4LogicalSkinSurface( "Al_endpiece_skin", Al_wrap_endpiece_log, GetOpticalSurface("Mirrsurf") );

  // There is a 10cm Foam Glass wrap around everything except the PMT area
  G4double foam_thick = 10.0*cm;
  G4double foam_XY = C16_width + 2*foam_thick;
  G4double foam_Z  = C16_depth + foam_thick - depth_ecal_pmt;
  G4Box *foam_box = new G4Box( "foam_box", foam_XY/2.0, foam_XY/2.0, foam_Z/2.0 );

  // Cut out the inside, but keep 10cm on one longitudinal face
  G4Box *foam_sub = new G4Box( "foam_sub", C16_width/2.0, C16_width/2.0, foam_Z/2.0 );
  G4SubtractionSolid *foam_wrap = new G4SubtractionSolid( "foam_wrap", foam_box, foam_sub, 0, G4ThreeVector(0.0,0.0,foam_thick) );
  G4LogicalVolume *foam_wrap_log = new G4LogicalVolume( foam_wrap, GetMaterial("SiO2_C16"), "foam_wrap_log" );
  G4ThreeVector R_TF1_Foam( (depth_ecal_pmt/2.0+foam_thick/2.0)*sin( angle_C16 ), 0.0, (depth_ecal_pmt/2.0+foam_thick/2.0)*cos( angle_C16 ) );
  G4ThreeVector R_Foam = R_C16 - R_TF1_Foam;
  new G4PVPlacement(bbrm_C16, R_Foam, foam_wrap_log, "Foam_phys", motherlog, false, 0 ,true);

  // There is a small Aluminium plate 5" upstream of the lead-glass face
  G4double Al_depth = 0.25*2.54*cm;
  G4double Al_width = 6.0*2.54*cm;
  G4double Distfrom_TF1 = 5.0*2.54*cm + C16_depth/2.0;
  G4ThreeVector R_TF1_Al( Distfrom_TF1*sin( angle_C16 ), 0.0, Distfrom_TF1*cos( angle_C16 ) );
  G4ThreeVector R_Al = R_C16 - R_TF1_Al;

  G4Box *Al_Plate_Box = new G4Box( "Al_Plate_Box", Al_width/2.0, Al_width/2.0, Al_depth/2.0 );
  G4LogicalVolume *Al_Plate_Log = new G4LogicalVolume( Al_Plate_Box, GetMaterial("Aluminum"), "Al_Plate_Log" );
  new G4PVPlacement( bbrm_C16, R_Al, Al_Plate_Log, "Al_Plate", motherlog, false, 0 );

  // VISUALS
  G4VisAttributes *TF1visatt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0 ) );
  G4VisAttributes *Alvisatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
  G4VisAttributes *ECALpmtvisatt = new G4VisAttributes( G4Colour( 0, 0, 1 ) );
  G4VisAttributes *C16WG_visatt = new G4VisAttributes( G4Colour(0.54, 0.53, 0.79) );
  G4VisAttributes *Foam_visatt = new G4VisAttributes( G4Colour( 0.0, 0.6, 0.6) );
  // C16 Mother:
  C16_Log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // Al Plate & Foil:
  Alvisatt->SetForceWireframe(true);
  Al_Plate_Log->SetVisAttributes( Alvisatt );
  Al_wrap_endpiece_log->SetVisAttributes( Alvisatt );
  // PMTs:
  ecal_PMT_log->SetVisAttributes( ECALpmtvisatt );
  // WaveGuides:
  C16_WG_Log->SetVisAttributes( C16WG_visatt );
  // Foam Wrap
  Foam_visatt->SetForceWireframe(true);
  foam_wrap_log->SetVisAttributes( Foam_visatt );

  /////////////////////////////////////////////////////////////////////////////////
  //   There are two options to build the TF1 ( /g4sbs/segmentC16 int )          //
  //     /g4sbs/segmentC16 N segments the TF1 into N "equal" sections, used       //
  //     to analyze dose rate as a function of longitudinal dimension of module  //
  //                                                                             //
  //     /g4sbs/segmentC16 0 builds the normal ECal modules                      //
  /////////////////////////////////////////////////////////////////////////////////

  if( fDetCon->GetC16Segmentation() <= 0 ) { //Default case: room temperature, no radiation damage:
    // Make a C16 Module  
    G4Box *Module_42 = new G4Box( "Module_42", width_42/2.0, width_42/2.0, depth_leadglass/2.0 );
    G4LogicalVolume *Module_42_log = new G4LogicalVolume( Module_42, GetMaterial("Special_Air"), "Module_42_log" );
    // Next, we want to make a subtraction solid for the Al:
    G4Box *Al_42 = new G4Box( "Al_42", (width_42 - alum_thick)/2.0, (width_42 - alum_thick)/2.0, depth_leadglass/2.0 + 1.0*cm );
    //
    G4SubtractionSolid *Al_wrap_42 = new G4SubtractionSolid( "Al_wrap_42", Module_42, Al_42, 0, G4ThreeVector( 0, 0, alum_thick + 1.0*cm ) );
    G4LogicalVolume *Al_wrap_42_log = new G4LogicalVolume( Al_wrap_42, GetMaterial("Aluminum"), "Al_wrap_42_log" );
    // Make lead-glass
    G4Box *LeadGlass_42 = new G4Box("LeadGlass_42", (width_42 - alum_thick - air_thick)/2.0, (width_42 - alum_thick - air_thick)/2.0, (depth_leadglass - alum_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_42_log = new G4LogicalVolume( LeadGlass_42, GetMaterial("TF1"), "LeadGlass_42_log" );

    // Define Sensitive Detector for lead-glass of type kCAL
    G4String C16TF1SDname = "Earm/C16TF1";
    G4String C16TF1collname = "C16TF1HitsCollection";
    G4SBSCalSD *C16TF1SD = NULL;
    
    if( !( C16TF1SD = (G4SBSCalSD*) sdman->FindSensitiveDetector(C16TF1SDname) ) ){
      G4cout << "Adding C16 TF1 Sensitive Detector to SDman..." << G4endl;
      C16TF1SD = new G4SBSCalSD( C16TF1SDname, C16TF1collname );
      fDetCon->fSDman->AddNewDetector( C16TF1SD );
      (fDetCon->SDlist).insert( C16TF1SDname );
      fDetCon->SDtype[C16TF1SDname] = G4SBS::kCAL;
      (C16TF1SD->detmap).depth = 1;

      fDetCon->SetThresholdTimeWindowAndNTimeBins( C16TF1SDname, 0.0*MeV, 100.0*ns, 25 );
    }
    // Assign "kCAL" sensitivity to the lead-glass:
    LeadGlass_42_log->SetSensitiveDetector( C16TF1SD );

    fDetCon->InsertSDboundaryVolume( C16_Log->GetName(), C16TF1SDname );
    
    // Place lead-glass and Al wrap inside module:
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, (alum_thick + air_thick)/2.0 ), LeadGlass_42_log, "LeadGlass_42_phys", Module_42_log, false, 0 );
    // Al:
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), Al_wrap_42_log, "Al_wrap_42_phys", Module_42_log, false, 0 );
    new G4LogicalSkinSurface( "Al_skin_42", Al_wrap_42_log, GetOpticalSurface("Mirrsurf") );

    // Construct C16
    G4int copyID = 0;
    for( G4int i=0; i<nrows; i++ ) {
      for( G4int j=0; j<ncols; j++ ) {
	G4double tempx = C16_width/2.0 - width_42/2.0 - j*width_42;
	G4double tempy = C16_width/2.0 - width_42/2.0 - i*width_42;
	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy, -depth_WG/2.0 - depth_ecal_pmt/2.0), Module_42_log, "C16_Module", C16_Log, false, copyID );
	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy, -C16_depth/2.0 + depth_leadglass + depth_WG/4.0), Al_wrap_endpiece_log, "Aluminum_wrap_endpiece_phys", C16_Log, false, copyID );
	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy, -depth_ecal_pmt/2.0 + depth_leadglass/2.0), C16_WG_Log, "C16_WG", C16_Log, false, copyID );
	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy,  depth_leadglass/2.0 + depth_WG/2.0), ecal_PMT_log, "C16_PMT", C16_Log, false, copyID );

	(C16TF1SD->detmap).Row[copyID] = i;
	(C16TF1SD->detmap).Col[copyID] = j;
	(C16TF1SD->detmap).LocalCoord[copyID] = G4ThreeVector( tempx, tempy, 0.0 );

	(C16SD->detmap).Row[copyID] = i;
	(C16SD->detmap).Col[copyID] = j;
	(C16SD->detmap).LocalCoord[copyID] = G4ThreeVector( tempx, tempy, 0.0 );
	
	copyID++;
      }
    }
    // Set Visuals
    Module_42_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
    LeadGlass_42_log->SetVisAttributes( TF1visatt );
    Al_wrap_42_log->SetVisAttributes( Alvisatt );  
  } else {
    // The strategy is to place the aluminum wrap within C16_Log (the mother volume), then iteratively
    // place TF1 modules in the longitudinal direction. Therefore, we can define different material properties
    // to each TF1 "segment" within a C16 Module.

    // Make a C16 Module  
    G4Box *Module_42 = new G4Box( "Module_42", width_42/2.0, width_42/2.0, depth_leadglass/2.0 );
    G4LogicalVolume *Module_42_log = new G4LogicalVolume( Module_42, GetMaterial("Special_Air"), "Module_42_log" );
    
    // Next, we want to make a subtraction solid for the Al:
    G4Box *Al_42 = new G4Box( "Al_42", (width_42 - alum_thick)/2.0, (width_42 - alum_thick)/2.0, depth_leadglass/2.0 + 1.0*cm );
    G4SubtractionSolid *Al_wrap_42 = new G4SubtractionSolid( "Al_wrap_42", Module_42, Al_42, 0, G4ThreeVector( 0, 0, alum_thick + 1.0*cm ) );
    G4LogicalVolume *Al_wrap_42_log = new G4LogicalVolume( Al_wrap_42, GetMaterial("Aluminum"), "Al_wrap_42_log" );
    new G4LogicalSkinSurface( "Al_skin_42", Al_wrap_42_log, GetOpticalSurface("Mirrsurf") );

    // Make TF1 - Logical Volume will be defined iteratively below in order to change material 
    // properties based on segmentation 
    //G4double Nsegments = 10.0;
    //    G4double Nsegments = G4double(fDetCon->GetC16Segmentation());
    //G4double segment_depth = (depth_leadglass - alum_thick - air_thick) / Nsegments;
    //Note: Because segment thickness is always the same for the detector and the material definition,
    // as long as nsegments * segmentthick > depth of lead glass, the material definition will always exist!
    G4double segment_depth = fDetCon->GetSegmentThickC16();

    G4int Nsegments = G4int( (depth_leadglass - alum_thick - air_thick ) / segment_depth ); //truncate the remainder.
    
    G4double remainder = (depth_leadglass - alum_thick - air_thick ) - Nsegments * segment_depth; //remainder is always non-negative!

    if( Nsegments > fDetCon->GetC16Segmentation() ){
      Nsegments = fDetCon->GetC16Segmentation();
      remainder = (depth_leadglass - alum_thick - air_thick ) - Nsegments * segment_depth;
    }
    
    if( Nsegments == 0 ){
      Nsegments = 1;
      segment_depth = (depth_leadglass - alum_thick - air_thick );
      remainder = 0.0;
    }
    
    G4Box *Segments_TF1 = new G4Box( "Segments_TF1", (width_42 - 2.0*alum_thick - 2.0*air_thick)/2.0, 
					   (width_42 - 2.0*alum_thick - 2.0*air_thick)/2.0, segment_depth/2.0 );

    //Always add the remainder to the last segment, since optical properties will be varying less with z deeper in the glass anyway!
    G4Box *LastSegment_TF1 = new G4Box( "LastSegment_TF1",
					(width_42 - 2.0*alum_thick - 2.0*air_thick)/2.0, 
					(width_42 - 2.0*alum_thick - 2.0*air_thick)/2.0,
					(segment_depth + remainder)/2.0 );  
    // TF1 is a Sensitive Detector of type kCAL. Sensitivity will be assigned to a LV
    // within the loop below:
    G4String C16TF1SDname = "Earm/C16TF1";
    G4String C16TF1collname = "C16TF1HitsCollection";
    G4SBSCalSD *C16TF1SD = NULL;
    
    if( !( C16TF1SD = (G4SBSCalSD*) sdman->FindSensitiveDetector(C16TF1SDname) ) ){
      G4cout << "Adding C16 TF1 Segmented Sensitive Detector to SDman..." << G4endl;
      C16TF1SD = new G4SBSCalSD( C16TF1SDname, C16TF1collname );
      fDetCon->fSDman->AddNewDetector( C16TF1SD );
      (fDetCon->SDlist).insert( C16TF1SDname );
      fDetCon->SDtype[C16TF1SDname] = G4SBS::kCAL;
      (C16TF1SD->detmap).depth = 0;

      fDetCon->SetThresholdTimeWindowAndNTimeBins( C16TF1SDname, 0.0*MeV, 100.0*ns, 25 );
    }

    fDetCon->InsertSDboundaryVolume( C16_Log->GetName(), C16TF1SDname );
    
    G4int cell_number = 0 ;    // cell #
    G4int TF1_number = 0 ;     // TF1 identifier
    for( G4int i = 0; i < nrows; i++ ) {
      for( G4int j = 0; j < ncols; j++ ) {
	G4double tempx = C16_width/2.0 - width_42/2.0 - j*width_42;
	G4double tempy = C16_width/2.0 - width_42/2.0 - i*width_42;
	
	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy, -C16_depth/2.0 + depth_leadglass/2.0), 
			   Al_wrap_42_log, "Aluminum_wrap_phys", C16_Log, false, cell_number );

	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy, -C16_depth/2.0 + depth_leadglass + depth_WG/4.0), 
			   Al_wrap_endpiece_log, "Aluminum_wrap_endpiece_phys", C16_Log, false, cell_number );

	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy, C16_depth/2.0 - depth_ecal_pmt - depth_WG/2.0), 
			   C16_WG_Log, "C16_WG_phys", C16_Log, false, cell_number );

	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy, C16_depth/2.0 - depth_ecal_pmt/2.0), 
			   ecal_PMT_log, "C16_PMT_phys", C16_Log, false, cell_number );

	(C16SD->detmap).Row[cell_number] = i;
	(C16SD->detmap).Col[cell_number] = j;
	(C16SD->detmap).LocalCoord[cell_number] = G4ThreeVector( tempx, tempy, 0.0 );
	
	cell_number++;

	for( G4int planeN = 0; planeN < Nsegments; planeN++ ) {
	  // ostringstream temp, temp1, temp2;
	  // temp << planeN;// segment #
	  // temp1 << i;    // row #
	  // temp2 << j;    // col #
	  // G4String tempstring  = temp.str();
	  // G4String tempstring1 = temp1.str();
	  // G4String tempstring2 = temp2.str();
	  // G4String seg_TF1_name_log  = "TF1_log_seg_"  + tempstring + "_row_" + tempstring1 + "_col_" + tempstring2;
	  // G4String seg_TF1_name_phys = "TF1_phys_seg_" + tempstring + "_row_" + tempstring1 + "_col_" + tempstring2;
	  // G4String seg_TF1_material = "TF1_anneal_" + tempstring;

	  TString material_name, lv_name, pv_name;
	  material_name.Form( "TF1_anneal_C16_row%d_col%d_z%d", i+1, j+1, planeN );

	  lv_name.Form( "TF1_log_z%d", planeN );
	  pv_name.Form( "TF1_pv_z%d", planeN );
	  
	  // Assign each TF1 segment a kCAL Sensitivity
	  G4LogicalVolume *Segments_TF1_log;
	  if( planeN + 1 < Nsegments || Nsegments == 1 ){
	    Segments_TF1_log = new G4LogicalVolume( Segments_TF1, GetMaterial(material_name.Data()), lv_name.Data() );
	  } else {
	    Segments_TF1_log = new G4LogicalVolume( LastSegment_TF1, GetMaterial(material_name.Data()), lv_name.Data() );
	  }
	    
	  Segments_TF1_log->SetSensitiveDetector( C16TF1SD );

	  G4VisAttributes *Segment_VisAtt = new G4VisAttributes( G4Colour( 0.8*(planeN/15.0)+0.20, 0.8*(planeN/15.0)+0.20, 0.0 ) );
	  Segments_TF1_log->SetVisAttributes( Segment_VisAtt );

	  // Place the TF1 segments longitudinally down the module
	  // Therefore, seg_0 corresponds to the face of C16
	  G4double tempz;
	  if( planeN + 1 < Nsegments ){
	    tempz = -C16_depth/2.0 + air_thick + alum_thick + (planeN + 0.5)*segment_depth;
	  } else {
	    tempz = -C16_depth/2.0 + air_thick + alum_thick + (Nsegments-1)*segment_depth +
	      0.5* ( segment_depth + remainder );
	  }
	  //G4double tempz = -C16_depth/2.0 + segment_depth/2.0 + air_thick + alum_thick + planeN*segment_depth;
	    
	  new G4PVPlacement(0, G4ThreeVector(tempx, tempy, tempz), Segments_TF1_log, pv_name.Data(), C16_Log, false, TF1_number );
	  
	  // Record useful information using Detmap:
	  (C16TF1SD->detmap).Row[TF1_number] = i;
	  (C16TF1SD->detmap).Col[TF1_number] = j;
	  (C16TF1SD->detmap).LocalCoord[TF1_number] = G4ThreeVector(tempx,tempy,tempz);
	  (C16TF1SD->detmap).Plane[TF1_number] = planeN;
	  TF1_number++;
	}
      }
    }
    // Set Visuals
    Module_42_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
    Al_wrap_42_log->SetVisAttributes( Alvisatt );
  }
}

void G4SBSECal::MakeBigCal(G4LogicalVolume *motherlog){
  printf("BigCal at %f deg\n", fAng/deg);

  // Pointer to SDmanager, used frequently in this routine
  G4SDManager *sdman = fDetCon->fSDman;

  //Parameters of ECAL for GEP:
  G4double xfpmin_cal = -160.0*cm;
  G4double xfpmax_cal = 160.0*cm;

  //number of available blocks:
  G4int nblocks_42 = 1000;
  G4int nblocks_40 = 700;
  G4int nblocks_38 = 300;

  //assume for now that all lead-glass blocks are 40 cm deep:
  G4double width_42 = 4.2*cm;
  G4double width_40 = 4.0*cm;
  G4double width_38 = 3.8*cm;

  G4double depth_42 = 34.3*cm;
  G4double depth_40 = 40.0*cm;
  G4double depth_38 = 45.0*cm;
  
  G4double depth_leadglass = 45.0*cm;
  G4double depth_ecal_pmt = 0.3*cm;
  G4double depth_lightguide_short = 15.0*cm;
  G4double radius_ecal_pmt = 1.25*cm;
  G4double depth_ecal_frontplate = 2.54*cm;
  G4double depth_CH2 = 20.0*cm; //This goes directly in front of CDET:
  G4double depth_CDET = 45.0*cm; // CDET starts 45 cm in front of ECAL:
  
  //Define "Earm" box a bit wider than total ECAL area:
  G4double width_earm = 150.0*cm;
  G4double height_earm = 340.0*cm;
  G4double depth_earm = depth_CH2 + depth_CDET + depth_ecal_pmt + depth_leadglass + depth_lightguide_short; // 
  
  G4int nSuperRows = 20;
  //These are the "yfp" coordinates of the box lower and upper edges. +yfp
  //Note: 0-9 and 10-19 are mirror images of each other!
  G4double yfpmin_super_rows[20] = { -56.0*cm,
				   -64.0*cm, -64.0*cm, -64.0*cm, -64.0*cm, -64.0*cm, -64.0*cm,
				   -56.0*cm,
				   -48.0*cm, -48.0*cm, -48.0*cm, -48.0*cm,
				   -56.0*cm,
				   -64.0*cm, -64.0*cm, -64.0*cm, -64.0*cm, -64.0*cm, -64.0*cm,
				   -56.0*cm };

  G4double yfpmax_super_rows[20] = { -24.0*cm, -8.0*cm, 8.0*cm, 24.0*cm, 32.0*cm,
				     48.0*cm, 48.0*cm, 56.0*cm, 56.0*cm, 56.0*cm,
				     56.0*cm, 56.0*cm, 56.0*cm, 48.0*cm, 48.0*cm,
				     32.0*cm, 24.0*cm, 8.0*cm, -8.0*cm, -24.0*cm };

  //G4double ycalo_min = -160.0*cm; //This can be adjusted if necessary!
  G4double ycalo_min = -164.0*cm;
  
  G4double xcalomin_super_rows[20], xcalomax_super_rows[20], width_super_rows[20];
  G4int ncol42_row[20], ncol40_row[20], ncol38_row[20];
  G4int nblock42_super_row[20], nblock40_super_row[20], nblock38_super_row[20];
  for(G4int row=0; row<20; row++){//
    xcalomin_super_rows[row] = yfpmin_super_rows[row];
    xcalomax_super_rows[row] = yfpmax_super_rows[row];
    width_super_rows[row] = xcalomax_super_rows[row] - xcalomin_super_rows[row];
    ncol42_row[row] = G4int( width_super_rows[row]/width_42 ) + 1; //truncate and add 1:
    ncol40_row[row] = G4int( width_super_rows[row]/width_40 ) + 1;
    ncol38_row[row] = G4int( width_super_rows[row]/width_38 ) + 1;

    //compute number of blocks required to fill a "super row":
    //Note that this is always a multiple of 4
    nblock42_super_row[row] = 4*ncol42_row[row];
    nblock40_super_row[row] = 4*ncol40_row[row];
    nblock38_super_row[row] = 4*ncol38_row[row];
  }
				     
  // Rotations, offsets, and distances from Target:
  G4RotationMatrix *bbrm = new G4RotationMatrix;
  bbrm->rotateY(-fAng);

  G4Box *earm_mother_box = new G4Box( "earm_mother_box", width_earm/2.0, height_earm/2.0, (depth_earm+1.0*mm)/2.0 );
  G4LogicalVolume *earm_mother_log = new G4LogicalVolume( earm_mother_box, GetMaterial("Air"), "earm_mother_log");

  //If "earm_mother_log" is in the step limiter list, make ECAL a total-absorber with "kCAL" sensitivity:
  if( (fDetCon->StepLimiterList).find( "earm_mother_log" ) != (fDetCon->StepLimiterList).end() ){
    earm_mother_log->SetUserLimits( new G4UserLimits( 0.0, 0.0, 0.0, DBL_MAX, DBL_MAX ) );

    G4String sdname = "Earm/ECAL_box";
    G4String collname = "ECAL_boxHitsCollection";
    G4SBSCalSD *earm_mother_SD = NULL;
    if( !( (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(sdname) ) ){
      G4cout << "Adding ECAL_box sensitive detector to SDman..." << G4endl;
      earm_mother_SD = new G4SBSCalSD( sdname, collname );
      fDetCon->fSDman->AddNewDetector( earm_mother_SD );
      (fDetCon->SDlist).insert( sdname );
      fDetCon->SDtype[sdname] = G4SBS::kCAL;
      (earm_mother_SD->detmap).depth = 0;
   
      earm_mother_log->SetSensitiveDetector( earm_mother_SD );
    }
  }

  //fDist should now be interpreted to mean the distance from the origin to the ****FRONT**** of the ECAL lead-glass:
  //G4double zback_ECAL = depth_earm/2.0 - depth_ecal_pmt;
  G4double zfront_ECAL = depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass;
  G4double R_Earm = fDist - zfront_ECAL;

  G4ThreeVector pos_ECAL( R_Earm*sin(fAng), 0.0, R_Earm*cos(fAng) );
  
  new G4PVPlacement( bbrm, pos_ECAL, earm_mother_log, "earm_mother_phys", motherlog, false, 0 );

  //assume lead-glass is surrounded by 1-mil (.001") thickness of mylar:
  //assume lead-glass is also surrounded by another mil of air:
  G4double mylar_thick = 0.001*2.54*cm;
  G4double air_thick = mylar_thick;
  
  G4double copper_thick = 0.25*mm;
  G4double al_thick = 0.50*mm;
  G4VisAttributes* hc_visAtt = new G4VisAttributes();

  G4double hcf_thick = copper_thick;
  const char* hcf_mat_name = "Copper";
  hc_visAtt->SetColour(1.0, 0.5, 0.0);
  // G4double hcf_thick = al_thick;
  // const char* hcf_mat_name = "Aluminum";
  // hc_visAtt->SetColour(0.7, 0.7, 0.7);
  // G4double hcf_thick = 0.0;
  // const char* hcf_mat_name = "Special_Air";
  // hc_visAtt->SetVisibility(0);
  
  //EFuchey 2017-01-11: Declaring sensitive detector for light guide 
  // shall be temporary, and not end in the repo...
  // G4String ECalLGSDname = "Earm/ECalLG";
  // G4String ECalLGcollname = "ECalLGHitsCollection";
  // G4SBSCalSD *ECalLGSD = NULL;
  // if( !( ECalLGSD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(ECalLGSDname) ) ){
  //   G4cout << "Adding ECal light guide Sensitive Detector to SDman..." << G4endl;
  //   ECalLGSD = new G4SBSCalSD( ECalLGSDname, ECalLGcollname );
  //   fDetCon->fSDman->AddNewDetector( ECalLGSD );
  //   (fDetCon->SDlist).insert(ECalLGSDname);
  //   fDetCon->SDtype[ECalLGSDname] = kCAL;
  //   (ECalLGSD->detmap).depth = 1;//?????
  // }
  
  //Now place things in ECAL:
  //Start with the lead-glass modules, PMTs and light guides:
  G4Tubs *LightGuide_42 = new G4Tubs("LightGuide_42", 0.0, 2.5*cm/2.0, (depth_lightguide_short+depth_38-depth_42)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LightGuide_42_log = new G4LogicalVolume( LightGuide_42, GetMaterial("Pyrex_Glass"), "LightGuide_42_log" );
  
  G4Tubs *LightGuide_40 = new G4Tubs("LightGuide_40", 0.0, 2.5*cm/2.0, (depth_lightguide_short+depth_38-depth_40)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LightGuide_40_log = new G4LogicalVolume( LightGuide_40, GetMaterial("Pyrex_Glass"), "LightGuide_40_log" );

  G4Tubs *LightGuide_38 = new G4Tubs("LightGuide_38", 0.0, 2.5*cm/2.0, (depth_lightguide_short+depth_38-depth_38)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LightGuide_38_log = new G4LogicalVolume( LightGuide_38, GetMaterial("Pyrex_Glass"), "LightGuide_38_log" );

  //EFuchey 2017-01-11: Need to make sensitive the three volumes above, to measure their dose.
  // shall be temporary, and not end in the repo...
  // LightGuide_42_log->SetSensitiveDetector( ECalLGSD );
  // LightGuide_40_log->SetSensitiveDetector( ECalLGSD );
  // LightGuide_38_log->SetSensitiveDetector( ECalLGSD );
  
  G4Tubs *LGWrap_42 = new G4Tubs( "LGWrap_42", 2.5*cm/2.0+air_thick, 2.5*cm/2.0 + air_thick + mylar_thick, (depth_lightguide_short+depth_38-depth_42)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LGWrap_42_log = new G4LogicalVolume( LGWrap_42, GetMaterial("Mylar"), "LGWrap_42_log" );

  G4Tubs *LGWrap_40 = new G4Tubs( "LGWrap_40", 2.5*cm/2.0+air_thick, 2.5*cm/2.0 + air_thick + mylar_thick, (depth_lightguide_short+depth_38-depth_40)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LGWrap_40_log = new G4LogicalVolume( LGWrap_40, GetMaterial("Mylar"), "LGWrap_40_log" );

  G4Tubs *LGWrap_38 = new G4Tubs( "LGWrap_38", 2.5*cm/2.0+air_thick, 2.5*cm/2.0 + air_thick + mylar_thick, (depth_lightguide_short+depth_38-depth_38)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LGWrap_38_log = new G4LogicalVolume( LGWrap_38, GetMaterial("Mylar"), "LGWrap_38_log" );
  
  G4Tubs *LG42 = new G4Tubs("LG42", 0.0, 2.5*cm/2.0 + air_thick + mylar_thick+0.1*mm, (depth_lightguide_short+depth_38-depth_42)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LG42_log = new G4LogicalVolume( LG42, GetMaterial("Special_Air"), "LG42_log" );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LightGuide_42_log, "LightGuide_42_phys", LG42_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LGWrap_42_log, "LGWrap_42_phys", LG42_log, false, 0 );

  G4Tubs *LG40 = new G4Tubs("LG40", 0.0, 2.5*cm/2.0 + air_thick + mylar_thick+0.1*mm, (depth_lightguide_short+depth_38-depth_40)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LG40_log = new G4LogicalVolume( LG40, GetMaterial("Special_Air"), "LG40_log" );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LightGuide_40_log, "LightGuide_40_phys", LG40_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LGWrap_40_log, "LGWrap_40_phys", LG40_log, false, 0 );
  
  G4Tubs *LG38 = new G4Tubs("LG38", 0.0, 2.5*cm/2.0 + air_thick + mylar_thick+0.1*mm, (depth_lightguide_short+depth_38-depth_38)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LG38_log = new G4LogicalVolume( LG38, GetMaterial("Special_Air"), "LG38_log" );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LightGuide_38_log, "LightGuide_38_phys", LG38_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LGWrap_38_log, "LGWrap_38_phys", LG38_log, false, 0 );

  G4Box *Module_42 = new G4Box( "Module_42", width_42/2.0, width_42/2.0, depth_42/2.0 );
  G4LogicalVolume *Module_42_log = new G4LogicalVolume( Module_42, GetMaterial("Special_Air"), "Module_42_log" );

  G4Box *Module_40 = new G4Box( "Module_40", width_40/2.0, width_40/2.0, depth_40/2.0 );
  G4LogicalVolume *Module_40_log = new G4LogicalVolume( Module_40, GetMaterial("Special_Air"), "Module_40_log" );

  G4Box *Module_38 = new G4Box( "Module_38", width_38/2.0, width_38/2.0, depth_38/2.0 );
  G4LogicalVolume *Module_38_log = new G4LogicalVolume( Module_38, GetMaterial("Special_Air"), "Module_38_log" );
  
  /*
  //Next, we want to make a subtraction solid for the mylar:
  G4Box *Mylar_42 = new G4Box( "Mylar_42", width_42/2.0 - mylar_thick, width_42/2.0 - mylar_thick, depth_42/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_42 = new G4SubtractionSolid( "Mylar_wrap_42", Module_42, Mylar_42, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm ) );
  G4LogicalVolume *Mylar_wrap_42_log = new G4LogicalVolume( Mylar_wrap_42, GetMaterial("Mylar"), "Mylar_wrap_42_log" );
  
  G4Box *Mylar_40 = new G4Box( "Mylar_40", width_40/2.0 - mylar_thick, width_40/2.0 - mylar_thick, depth_40/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_40 = new G4SubtractionSolid( "Mylar_wrap_40", Module_40, Mylar_40, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm ) );
  G4LogicalVolume *Mylar_wrap_40_log = new G4LogicalVolume( Mylar_wrap_40, GetMaterial("Mylar"), "Mylar_wrap_40_log" );

  G4Box *Mylar_38 = new G4Box( "Mylar_38", width_38/2.0 - mylar_thick, width_38/2.0 - mylar_thick, depth_38/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_38 = new G4SubtractionSolid( "Mylar_wrap_38", Module_38, Mylar_38, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm) );
  G4LogicalVolume *Mylar_wrap_38_log = new G4LogicalVolume( Mylar_wrap_38, GetMaterial("Mylar"), "Mylar_wrap_38_log" );
  */
  G4Box *Module_42_k = new G4Box( "Module_42_k", width_42/2.0, width_42/2.0, depth_42/2.0-0.25*mm );
  G4Box *Module_40_k = new G4Box( "Module_42_k", width_40/2.0, width_40/2.0, depth_40/2.0-0.25*mm );
  G4Box *Module_38_k = new G4Box( "Module_42_k", width_38/2.0, width_38/2.0, depth_38/2.0-0.25*mm );
  
  // Solid for heat conducting foils: ouside of the wrapping, supposedly
  G4Box *hc_42 = new G4Box( "hc_42", width_42/2.0 - hcf_thick, width_42/2.0 - hcf_thick, depth_42/2.0 );
  G4SubtractionSolid *hc_foil_42 = new G4SubtractionSolid( "hc_foil_42", Module_42_k, hc_42, 0, G4ThreeVector( 0, 0, 0 ) );
  G4LogicalVolume *hc_foil_42_log = new G4LogicalVolume( hc_foil_42, GetMaterial(hcf_mat_name), "hc_foil_42_log" );
  hc_foil_42_log->SetVisAttributes(hc_visAtt);

  G4Box *hc_40 = new G4Box( "hc_40", width_40/2.0 - hcf_thick, width_40/2.0 - hcf_thick, depth_40/2.0 );
  G4SubtractionSolid *hc_foil_40 = new G4SubtractionSolid( "hc_foil_40", Module_40_k, hc_40, 0, G4ThreeVector( 0, 0, 0 ) );
  G4LogicalVolume *hc_foil_40_log = new G4LogicalVolume( hc_foil_40, GetMaterial(hcf_mat_name), "hc_foil_40_log" );
  hc_foil_40_log->SetVisAttributes(hc_visAtt);
  
  G4Box *hc_38 = new G4Box( "hc_38", width_38/2.0 - hcf_thick, width_38/2.0 - hcf_thick, depth_38/2.0 );
  G4SubtractionSolid *hc_foil_38 = new G4SubtractionSolid( "hc_foil_38", Module_38_k, hc_38, 0, G4ThreeVector( 0, 0, 0 ) );
  G4LogicalVolume *hc_foil_38_log = new G4LogicalVolume( hc_foil_38, GetMaterial(hcf_mat_name), "hc_foil_38_log" );
  hc_foil_38_log->SetVisAttributes(hc_visAtt);
  
  //Next, we want to make a subtraction solid for the mylar:
  G4Box *Mylar_42 = new G4Box( "Mylar_42", width_42/2.0 - hcf_thick - mylar_thick, width_42/2.0 - hcf_thick - mylar_thick, depth_42/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_42 = new G4SubtractionSolid( "Mylar_wrap_42", hc_42, Mylar_42, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm ) );
  G4LogicalVolume *Mylar_wrap_42_log = new G4LogicalVolume( Mylar_wrap_42, GetMaterial("Mylar"), "Mylar_wrap_42_log" );
  
  G4Box *Mylar_40 = new G4Box( "Mylar_40", width_40/2.0 - hcf_thick - mylar_thick, width_40/2.0 - hcf_thick - mylar_thick, depth_40/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_40 = new G4SubtractionSolid( "Mylar_wrap_40", hc_40, Mylar_40, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm ) );
  G4LogicalVolume *Mylar_wrap_40_log = new G4LogicalVolume( Mylar_wrap_40, GetMaterial("Mylar"), "Mylar_wrap_40_log" );

  G4Box *Mylar_38 = new G4Box( "Mylar_38", width_38/2.0 - hcf_thick - mylar_thick, width_38/2.0 - hcf_thick - mylar_thick, depth_38/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_38 = new G4SubtractionSolid( "Mylar_wrap_38", hc_38, Mylar_38, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm) );
  G4LogicalVolume *Mylar_wrap_38_log = new G4LogicalVolume( Mylar_wrap_38, GetMaterial("Mylar"), "Mylar_wrap_38_log" );
  
  
  ////// Define Sensitive Detector for lead-glass:
  G4String ECalTF1SDname = "Earm/ECalTF1";
  G4String ECalTF1collname = "ECalTF1HitsCollection";
  G4SBSCalSD *ECalTF1SD = NULL;
    
  if( !( ECalTF1SD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(ECalTF1SDname) ) ){
    G4cout << "Adding ECal TF1 Sensitive Detector to SDman..." << G4endl;
    ECalTF1SD = new G4SBSCalSD( ECalTF1SDname, ECalTF1collname );
    fDetCon->fSDman->AddNewDetector( ECalTF1SD );
    (fDetCon->SDlist).insert(ECalTF1SDname);
    fDetCon->SDtype[ECalTF1SDname] = G4SBS::kCAL;
    //fDetCon->SDarm[ECalTF1SDname] = G4SBS::kEarm;

    (ECalTF1SD->detmap).depth = 1;

    fDetCon->SetThresholdTimeWindowAndNTimeBins( ECalTF1SDname, 0.0*MeV, 100.0*ns, 25 );
  }

  fDetCon->InsertSDboundaryVolume( earm_mother_log->GetName(), ECalTF1SDname );
  
  //Make lead-glass and place in modules:
  
  if( fDetCon->GetC16Segmentation() <= 0 ){
    /*
    // EFuchey:2017/03/03 Was that correct anyhow ??? Don't think so...
    // block section shall be: module_section - 2*mylar_thick - 2*airthick ( - 2*hcf_thick, but I added that)
    G4Box *LeadGlass_42 = new G4Box("LeadGlass_42", width_42/2.0 - mylar_thick - air_thick, width_42/2.0 - mylar_thick - air_thick, (depth_42 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_42_log = new G4LogicalVolume( LeadGlass_42, GetMaterial("TF1"), "LeadGlass_42_log" );

    G4Box *LeadGlass_40 = new G4Box("LeadGlass_40", width_40/2.0 - mylar_thick - air_thick, width_40/2.0 - mylar_thick - air_thick, (depth_40 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_40_log = new G4LogicalVolume( LeadGlass_40, GetMaterial("TF1"), "LeadGlass_40_log" );

    G4Box *LeadGlass_38 = new G4Box("LeadGlass_38", width_38/2.0 - mylar_thick - air_thick, width_38/2.0 - mylar_thick - air_thick, (depth_38 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_38_log = new G4LogicalVolume( LeadGlass_38, GetMaterial("TF1"), "LeadGlass_38_log" );
    */
    G4Box *LeadGlass_42 = new G4Box("LeadGlass_42", width_42/2.0 - hcf_thick - mylar_thick - air_thick, width_42/2.0 - hcf_thick - mylar_thick - air_thick, (depth_42 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_42_log = new G4LogicalVolume( LeadGlass_42, GetMaterial("TF1"), "LeadGlass_42_log" );

    G4Box *LeadGlass_40 = new G4Box("LeadGlass_40", width_40/2.0 - hcf_thick - mylar_thick - air_thick, width_40/2.0 - hcf_thick - mylar_thick - air_thick, (depth_40 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_40_log = new G4LogicalVolume( LeadGlass_40, GetMaterial("TF1"), "LeadGlass_40_log" );

    G4Box *LeadGlass_38 = new G4Box("LeadGlass_38", width_38/2.0 - hcf_thick - mylar_thick - air_thick, width_38/2.0 - hcf_thick - mylar_thick - air_thick, (depth_38 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_38_log = new G4LogicalVolume( LeadGlass_38, GetMaterial("TF1"), "LeadGlass_38_log" );
    
    
    //Assign "kCAL" sensitivity to the lead-glass:
    LeadGlass_42_log->SetSensitiveDetector( ECalTF1SD );
    LeadGlass_40_log->SetSensitiveDetector( ECalTF1SD );
    LeadGlass_38_log->SetSensitiveDetector( ECalTF1SD );

    G4VisAttributes *TF1visatt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0 ) );
    LeadGlass_42_log->SetVisAttributes( TF1visatt );
    LeadGlass_40_log->SetVisAttributes( TF1visatt );
    LeadGlass_38_log->SetVisAttributes( TF1visatt );
    
    //Positioning of lead-glass in module:
    // z + Lz/2 - m/2 - a/2 = Lz/2 --> z = m/2 + a/2
    //lead-glass:
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, (mylar_thick + air_thick)/2.0 ), LeadGlass_42_log, "LeadGlass_42_phys", Module_42_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, (mylar_thick + air_thick)/2.0 ), LeadGlass_40_log, "LeadGlass_40_phys", Module_40_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, (mylar_thick + air_thick)/2.0 ), LeadGlass_38_log, "LeadGlass_38_phys", Module_38_log, false, 0 );
  } else {
    G4int nsegments = fDetCon->GetC16Segmentation();
    G4double segthick = fDetCon->GetSegmentThickC16();

    //Segment the lead-glass longitudinally so that different optical properties can be defined:
    G4int nseg_42 = G4int( (depth_42 - mylar_thick - air_thick)/segthick );
    if( nseg_42 > nsegments ){
      nseg_42 = nsegments;
    }
    G4double remainder_42 = (depth_42 - mylar_thick - air_thick) - nseg_42 * segthick; //remainder is always non-negative!

    if( nseg_42 == 0 ){
      nseg_42 = 1;
      segthick = (depth_42 - mylar_thick - air_thick);
      remainder_42 = 0.0;
    }

    G4Box *segment_42 = new G4Box( "segment_42", (width_42 - 2.*(mylar_thick+air_thick) )/2.0, (width_42 - 2.*(mylar_thick+air_thick) )/2.0, segthick/2.0 );
    G4Box *lastsegment_42 = new G4Box( "lastsegment_42", (width_42 - 2.*(mylar_thick+air_thick) )/2.0, (width_42 - 2.*(mylar_thick+air_thick) )/2.0, (segthick+remainder_42)/2.0 );

    G4double ztemp = -depth_42/2.0 + mylar_thick + air_thick;
    
    for( int seg42=0; seg42<nseg_42; seg42++ ){
      TString material_name, lv_name, pv_name;
      material_name.Form( "TF1_anneal_ECAL_z%d", seg42 );
      lv_name.Form("segment_42_log_%d", seg42 );
      pv_name.Form("segment_42_phys_%d", seg42 );
      G4LogicalVolume *Seg42_log;
      if( seg42 + 1 < nseg_42 || nseg_42 == 1 ){
	Seg42_log = new G4LogicalVolume( segment_42, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += segthick/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg42_log, pv_name.Data(), Module_42_log, false, 0 );
	ztemp += segthick/2.0;
      } else {
	Seg42_log = new G4LogicalVolume( lastsegment_42, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += (segthick+remainder_42)/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg42_log, pv_name.Data(), Module_42_log, false, 0 );
	ztemp += (segthick+remainder_42)/2.0;
      }
      Seg42_log->SetSensitiveDetector( ECalTF1SD );

      G4VisAttributes *Segment_VisAtt = new G4VisAttributes( G4Colour( 0.8*(seg42/15.0)+0.20, 0.8*(seg42/15.0)+0.20, 0.0 ) );
      Seg42_log->SetVisAttributes( Segment_VisAtt );
    }


    //Segment the lead-glass longitudinally so that different optical properties can be defined:
    G4int nseg_40 = G4int( (depth_40 - mylar_thick - air_thick)/segthick );
    if( nseg_40 > nsegments ){
      nseg_40 = nsegments;
    }
    G4double remainder_40 = (depth_40 - mylar_thick - air_thick) - nseg_40 * segthick; //remainder is always non-negative!

    if( nseg_40 == 0 ){
      nseg_40 = 1;
      segthick = (depth_40 - mylar_thick - air_thick);
      remainder_40 = 0.0;
    }

    G4Box *segment_40 = new G4Box( "segment_40", (width_40 - 2.*(mylar_thick+air_thick) )/2.0, (width_40 - 2.*(mylar_thick+air_thick) )/2.0, segthick/2.0 );
    G4Box *lastsegment_40 = new G4Box( "lastsegment_40", (width_40 - 2.*(mylar_thick+air_thick) )/2.0, (width_40 - 2.*(mylar_thick+air_thick) )/2.0, (segthick+remainder_40)/2.0 );

    ztemp = -depth_40/2.0 + mylar_thick + air_thick;
    
    for( int seg40=0; seg40<nseg_40; seg40++ ){
      TString material_name, lv_name, pv_name;
      material_name.Form( "TF1_anneal_ECAL_z%d", seg40 );
      lv_name.Form("segment_40_log_%d", seg40 );
      pv_name.Form("segment_40_phys_%d", seg40 );
      G4LogicalVolume *Seg40_log;
      if( seg40 + 1 < nseg_40 || nseg_40 == 1 ){
	Seg40_log = new G4LogicalVolume( segment_40, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += segthick/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg40_log, pv_name.Data(), Module_40_log, false, 0 );
	ztemp += segthick/2.0;
      } else {
	Seg40_log = new G4LogicalVolume( lastsegment_40, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += (segthick+remainder_40)/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg40_log, pv_name.Data(), Module_40_log, false, 0 );
	ztemp += (segthick+remainder_40)/2.0;
      }
      Seg40_log->SetSensitiveDetector( ECalTF1SD );

      G4VisAttributes *Segment_VisAtt = new G4VisAttributes( G4Colour( 0.8*(seg40/15.0)+0.20, 0.8*(seg40/15.0)+0.20, 0.0 ) );
      Seg40_log->SetVisAttributes( Segment_VisAtt );
    }

    //Segment the lead-glass longitudinally so that different optical properties can be defined:
    G4int nseg_38 = G4int( (depth_38 - mylar_thick - air_thick)/segthick );
    if( nseg_38 > nsegments ){
      nseg_38 = nsegments;
    }
    G4double remainder_38 = (depth_38 - mylar_thick - air_thick) - nseg_38 * segthick; //remainder is always non-negative!

    if( nseg_38 == 0 ){
      nseg_38 = 1;
      segthick = (depth_38 - mylar_thick - air_thick);
      remainder_38 = 0.0;
    }

    G4Box *segment_38 = new G4Box( "segment_38", (width_38 - 2.*(mylar_thick+air_thick) )/2.0, (width_38 - 2.*(mylar_thick+air_thick) )/2.0, segthick/2.0 );
    G4Box *lastsegment_38 = new G4Box( "lastsegment_38", (width_38 - 2.*(mylar_thick+air_thick) )/2.0, (width_38 - 2.*(mylar_thick+air_thick) )/2.0, (segthick+remainder_38)/2.0 );

    ztemp = -depth_38/2.0 + mylar_thick + air_thick;
    
    for( int seg38=0; seg38<nseg_38; seg38++ ){
      TString material_name, lv_name, pv_name;
      material_name.Form( "TF1_anneal_ECAL_z%d", seg38 );
      lv_name.Form("segment_38_log_%d", seg38 );
      pv_name.Form("segment_38_phys_%d", seg38 );
      G4LogicalVolume *Seg38_log;
      if( seg38 + 1 < nseg_38 || nseg_38 == 1 ){
	Seg38_log = new G4LogicalVolume( segment_38, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += segthick/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg38_log, pv_name.Data(), Module_38_log, false, 0 );
	ztemp += segthick/2.0;
      } else {
	Seg38_log = new G4LogicalVolume( lastsegment_38, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += (segthick+remainder_38)/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg38_log, pv_name.Data(), Module_38_log, false, 0 );
	ztemp += (segthick+remainder_38)/2.0;
      }
      Seg38_log->SetSensitiveDetector( ECalTF1SD );

      G4VisAttributes *Segment_VisAtt = new G4VisAttributes( G4Colour( 0.8*(seg38/15.0)+0.20, 0.8*(seg38/15.0)+0.20, 0.0 ) );
      Seg38_log->SetVisAttributes( Segment_VisAtt );
    }
    
  }
    
  //Place lead-glass and mylar wrap inside module:
  
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), hc_foil_42_log, "hc_foil_42_phys", Module_42_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), hc_foil_40_log, "hc_foil_40_phys", Module_40_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), hc_foil_38_log, "hc_foil_38_phys", Module_38_log, false, 0 );
  
  //mylar:
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), Mylar_wrap_42_log, "Mylar_wrap_42_phys", Module_42_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), Mylar_wrap_40_log, "Mylar_wrap_40_phys", Module_40_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), Mylar_wrap_38_log, "Mylar_wrap_38_phys", Module_38_log, false, 0 );

  new G4LogicalSkinSurface( "Mylar_skin_42", Mylar_wrap_42_log, GetOpticalSurface("Mirrsurf") );
  new G4LogicalSkinSurface( "Mylar_skin_40", Mylar_wrap_40_log, GetOpticalSurface("Mirrsurf") );
  new G4LogicalSkinSurface( "Mylar_skin_38", Mylar_wrap_38_log, GetOpticalSurface("Mirrsurf") );

  new G4LogicalSkinSurface( "lgwrap_42", LGWrap_42_log, GetOpticalSurface("Mirrsurf") );
  new G4LogicalSkinSurface( "lgwrap_40", LGWrap_40_log, GetOpticalSurface("Mirrsurf") );
  new G4LogicalSkinSurface( "lgwrap_38", LGWrap_38_log, GetOpticalSurface("Mirrsurf") );

  

  //Filling order:
  //G4int index_super_rows[20] = {0,19,1,18,2,17,3,16,4,15,5,14,6,13,7,12,8,11,9,10};
  //Now, start placing the modules:
  //Approach: Put larger blocks on bottom, smaller blocks on top? 

  //Next: PMTs:
  G4Tubs *ecal_PMT = new G4Tubs( "ecal_PMT", 0.0, radius_ecal_pmt, depth_ecal_pmt/2.0, 0.0, twopi );
  G4LogicalVolume *ecal_PMT_log = new G4LogicalVolume( ecal_PMT, GetMaterial("Photocathode_material_ecal"), "ecal_PMT_log" );

  G4String sdname = "Earm/ECAL";
  G4String collname = "ECALHitsCollection";

  G4SBSECalSD *ECalSD = NULL;

  if( !( (G4SBSECalSD*) sdman->FindSensitiveDetector( sdname ) ) ){
    G4cout << "Adding ECAL PMT sensitive detector to SDman..." << G4endl;
    ECalSD = new G4SBSECalSD( sdname, collname );
    sdman->AddNewDetector( ECalSD );
    (fDetCon->SDlist).insert( sdname );
    fDetCon->SDtype[sdname] = G4SBS::kECAL;
    (ECalSD->detmap).depth = 0;
  }

  fDetCon->InsertSDboundaryVolume( earm_mother_log->GetName(), sdname );
  
  ecal_PMT_log->SetSensitiveDetector( ECalSD );
  
  int lastrow42 = 0;
  int nused42 = 0, nused40=0, nused38=0;
  int icell=0, global_row=0, global_col=0;

  int nrows42 = 0, nrows40 = 0, nrows38 = 0;

  G4double ysum = ycalo_min;

  G4VisAttributes *Alvisatt = new G4VisAttributes( G4Colour( 0.4, 0.4, 0.4 ) );

  //Fill out bottom with Al:
  G4Box *bottom_Al = new G4Box( "bottom_Al", width_earm/2.0, (ycalo_min + height_earm/2.0)/2.0, depth_leadglass/2.0 );
  G4LogicalVolume *bottom_Al_log = new G4LogicalVolume( bottom_Al, GetMaterial("Al"), "bottom_Al_log" );
  new G4PVPlacement( 0, G4ThreeVector( 0, (ycalo_min - height_earm/2.0)/2.0, depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 ), bottom_Al_log, "bottom_Al_phys", earm_mother_log, false, 0 ); 

  // bottom_Al_log->SetVisAttributes( Alvisatt );
  bottom_Al_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  for( int super_row=0; super_row<20; super_row++ ){
    for( int sub_row=0; sub_row<4; sub_row++ ){
      global_row = sub_row + 4*super_row;
      int ncol;

      G4double xlow_row, xhigh_row;
      //check whether there are enough 4.2 cm blocks to fill this row:
      if( nused42 + ncol42_row[super_row] <= nblocks_42 ){
	//Fill this row with 4.2-cm blocks:
	ncol = ncol42_row[super_row];
	for( int col=0; col<ncol; col++ ){
	  G4ThreeVector modpos( xcalomin_super_rows[super_row] + (col + 0.5*(1+sub_row%2))*width_42,
				ysum + 0.5*width_42,
				zfront_ECAL + depth_42/2.0 );
				// depth_earm/2.0 - depth_ecal_pmt - depth_leadglass/2.0 );
	  new G4PVPlacement( 0, modpos, Module_42_log, "Module_42_phys", earm_mother_log, false, icell );

	  (ECalTF1SD->detmap).Row[icell] = global_row;
	  (ECalTF1SD->detmap).Col[icell] = col;
	  (ECalTF1SD->detmap).LocalCoord[icell] = modpos;

	  G4ThreeVector pmtpos( modpos.x(), modpos.y(), depth_earm/2.0 - depth_ecal_pmt/2.0 );
	  new G4PVPlacement( 0, pmtpos, ecal_PMT_log, "ecal_PMT_phys", earm_mother_log, false, icell );

	  (ECalSD->detmap).Row[icell] = global_row;
	  (ECalSD->detmap).Col[icell] = col;
	  (ECalSD->detmap).LocalCoord[icell] = pmtpos;

	  //Add light-guide with mylar wrap:
	  G4ThreeVector LGpos( modpos.x(), modpos.y(), modpos.z() + depth_42/2.0 + LG42->GetZHalfLength() );
	  new G4PVPlacement( 0, LGpos, LG42_log, "LG42_phys", earm_mother_log, false, icell );
	  // new G4PVPlacement( 0, LGpos, LightGuide_42_log, "LightGuide_42_phys", earm_mother_log, false, icell );
	  // new G4PVPlacement( 0, LGpos, LGWrap_42_log, "LGWrap_42_phys", earm_mother_log, false, icell );

	  // //EFuchey 2017-01-12: Need to make sensitive the three volumes above, to measure their dose.
	  // // shall be temporary, and not end in the repo...
	  // (ECalLGSD->detmap).Row[icell] = global_row;
	  // (ECalLGSD->detmap).Col[icell] = col;
	  // (ECalLGSD->detmap).LocalCoord[icell] = modpos;
	  
	  if( col == 0 ) xlow_row = modpos.x() - 0.5*width_42;
	  if( col+1 == ncol ) xhigh_row = modpos.x() + 0.5*width_42;
	  
	  icell++;
	}

	ysum += width_42;
	
	nused42 += ncol42_row[super_row];
	nrows42++;
	
	char prefix[255];
	
	//Fill out row with Al:
	if( xlow_row > -width_earm/2.0 ){
	  sprintf( prefix, "Albox1_row%d", super_row );
	  G4String boxname = prefix;
	  
	  G4Box *Al_filler1 = new G4Box( boxname, (xlow_row + width_earm/2.0)/2.0, width_42/2.0, depth_leadglass/2.0 );
	  
	  G4String logname = boxname + "log";
	  
	  G4LogicalVolume *Al_filler1_log = new G4LogicalVolume( Al_filler1, GetMaterial("Al"), logname );
	  
	  //Al_filler1_log->SetVisAttributes( Alvisatt );
	  Al_filler1_log->SetVisAttributes( G4VisAttributes::GetInvisible() );

	  G4ThreeVector pos( (xlow_row - width_earm/2.0)/2.0,
			     ysum - 0.5*width_42,
			     depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 );

	  G4String physname = boxname + "phys";
	  
	  new G4PVPlacement( 0, pos, Al_filler1_log, physname, earm_mother_log, false, 0 );
	}
	if( xhigh_row < width_earm/2.0 ){
	  sprintf( prefix, "Albox2_row%d", super_row );
	  G4String boxname = prefix;
	  
	  G4Box *Al_filler2 = new G4Box( boxname, (-xhigh_row + width_earm/2.0)/2.0, width_42/2.0, depth_leadglass/2.0 );
	  
	  G4String logname = boxname + "log";
	  
	  G4LogicalVolume *Al_filler2_log = new G4LogicalVolume( Al_filler2, GetMaterial("Al"), logname );
	  //Al_filler2_log->SetVisAttributes( Alvisatt );

	  Al_filler2_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
	  
	  //G4ThreeVector pos( (xhigh_row + width_earm/2.0)/2.0, modpos.y(), modpos.z() );
	  G4ThreeVector pos( (xhigh_row + width_earm/2.0)/2.0,
			     ysum - 0.5*width_42,
			     depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 );
	  
	  G4String physname = boxname + "phys";
	  
	  new G4PVPlacement( 0, pos, Al_filler2_log, physname, earm_mother_log, false, 0 );
	}
	
      } else if( nused40 + ncol40_row[super_row] <= nblocks_40 ){
	ncol = ncol40_row[super_row];
	for( int col=0; col<ncol; col++ ){
	  G4ThreeVector modpos( xcalomin_super_rows[super_row] + (col + 0.5*(1+sub_row%2))*width_40,
				ysum + 0.5*width_40,
				zfront_ECAL + depth_40/2.0 );
				//depth_earm/2.0 - depth_ecal_pmt - depth_leadglass/2.0 );
	  new G4PVPlacement( 0, modpos, Module_40_log, "Module_40_phys", earm_mother_log, false, icell );

	  (ECalTF1SD->detmap).Row[icell] = global_row;
	  (ECalTF1SD->detmap).Col[icell] = col;
	  (ECalTF1SD->detmap).LocalCoord[icell] = modpos;

	  G4ThreeVector pmtpos( modpos.x(), modpos.y(), depth_earm/2.0 - depth_ecal_pmt/2.0 );
	  new G4PVPlacement( 0, pmtpos, ecal_PMT_log, "ecal_PMT_phys", earm_mother_log, false, icell );

	  (ECalSD->detmap).Row[icell] = global_row;
	  (ECalSD->detmap).Col[icell] = col;
	  (ECalSD->detmap).LocalCoord[icell] = pmtpos;

	  //Add light-guide with mylar wrap:
	  G4ThreeVector LGpos( modpos.x(), modpos.y(), modpos.z() + depth_40/2.0 + LightGuide_40->GetZHalfLength() );
	  new G4PVPlacement( 0, LGpos, LG40_log, "LG40_phys", earm_mother_log, false, icell );
	  // new G4PVPlacement( 0, LGpos, LightGuide_40_log, "LightGuide_40_phys", earm_mother_log, false, icell );
	  // new G4PVPlacement( 0, LGpos, LGWrap_40_log, "LGWrap_40_phys", earm_mother_log, false, icell );
	  
	  //EFuchey 2017-01-12: Need to make sensitive the three volumes above, to measure their dose.
	  // shall be temporary, and not end in the repo...
	  // (ECalLGSD->detmap).Row[icell] = global_row;
	  // (ECalLGSD->detmap).Col[icell] = col;
	  // (ECalLGSD->detmap).LocalCoord[icell] = modpos;
	  
	  if( col == 0 ) xlow_row = modpos.x() - 0.5*width_40;
	  if( col+1 == ncol ) xhigh_row = modpos.x() + 0.5*width_40;
	  
	  icell++;
	}
	ysum += width_40;
	nused40 += ncol40_row[super_row];
	nrows40++;

	char prefix[255];
	//Fill out row with Al:
	if( xlow_row > -width_earm/2.0 ){
	  sprintf( prefix, "Albox1_row%d", super_row );
	  G4String boxname = prefix;
	  
	  G4Box *Al_filler1 = new G4Box( boxname, (xlow_row + width_earm/2.0)/2.0, width_40/2.0, depth_leadglass/2.0 );
	  
	  G4String logname = boxname + "log";
	  
	  G4LogicalVolume *Al_filler1_log = new G4LogicalVolume( Al_filler1, GetMaterial("Al"), logname );
	  //Al_filler1_log->SetVisAttributes( Alvisatt );

	  Al_filler1_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
	  
	  G4ThreeVector pos( (xlow_row - width_earm/2.0)/2.0,
			     ysum - 0.5*width_40,
			     depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 );

	  G4String physname = boxname + "phys";
	  
	  new G4PVPlacement( 0, pos, Al_filler1_log, physname, earm_mother_log, false, 0 );
	}
	if( xhigh_row < width_earm/2.0 ){
	  sprintf( prefix, "Albox2_row%d", super_row );
	  G4String boxname = prefix;
	  
	  G4Box *Al_filler2 = new G4Box( boxname, (-xhigh_row + width_earm/2.0)/2.0, width_40/2.0, depth_leadglass/2.0 );
	  
	  G4String logname = boxname + "log";
	  
	  G4LogicalVolume *Al_filler2_log = new G4LogicalVolume( Al_filler2, GetMaterial("Al"), logname );
	  //Al_filler2_log->SetVisAttributes( Alvisatt );

	  Al_filler2_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
	  
	  //G4ThreeVector pos( (xhigh_row + width_earm/2.0)/2.0, modpos.y(), modpos.z() );
	  G4ThreeVector pos( (xhigh_row + width_earm/2.0)/2.0,
			     ysum - 0.5*width_40,
			     depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 );
	  
	  G4String physname = boxname + "phys";
	  
	  new G4PVPlacement( 0, pos, Al_filler2_log, physname, earm_mother_log, false, 0 );
	}
	
      } else if( nused38 + ncol38_row[super_row] <= nblocks_38 ){
	ncol = ncol38_row[super_row];
	for( int col=0; col<ncol; col++ ){
	  G4ThreeVector modpos( xcalomin_super_rows[super_row] + (col + 0.5*(1+sub_row%2))*width_38,
				ysum + 0.5*width_38,
				zfront_ECAL + depth_38/2.0 );
				//depth_earm/2.0 - depth_ecal_pmt - depth_leadglass/2.0 );
	  new G4PVPlacement( 0, modpos, Module_38_log, "Module_38_phys", earm_mother_log, false, icell );

	  (ECalTF1SD->detmap).Row[icell] = global_row;
	  (ECalTF1SD->detmap).Col[icell] = col;
	  (ECalTF1SD->detmap).LocalCoord[icell] = modpos;

	  G4ThreeVector pmtpos( modpos.x(), modpos.y(), depth_earm/2.0 - depth_ecal_pmt/2.0 );
	  new G4PVPlacement( 0, pmtpos, ecal_PMT_log, "ecal_PMT_phys", earm_mother_log, false, icell );

	  (ECalSD->detmap).Row[icell] = global_row;
	  (ECalSD->detmap).Col[icell] = col;
	  (ECalSD->detmap).LocalCoord[icell] = pmtpos;

	   //Add light-guide with mylar wrap:
	  G4ThreeVector LGpos( modpos.x(), modpos.y(), modpos.z() + depth_38/2.0 + LightGuide_38->GetZHalfLength() );
	  new G4PVPlacement( 0, LGpos, LG38_log, "LG38_phys", earm_mother_log, false, icell );
	  // new G4PVPlacement( 0, LGpos, LightGuide_38_log, "LightGuide_38_phys", earm_mother_log, false, icell );
	  // new G4PVPlacement( 0, LGpos, LGWrap_38_log, "LGWrap_38_phys", earm_mother_log, false, icell );

	  //EFuchey 2017-01-12: Need to make sensitive the three volumes above, to measure their dose.
	  // shall be temporary, and not end in the repo...
	  // (ECalLGSD->detmap).Row[icell] = global_row;
	  // (ECalLGSD->detmap).Col[icell] = col;
	  // (ECalLGSD->detmap).LocalCoord[icell] = modpos;
	  
	  if( col == 0 ) xlow_row = modpos.x() - 0.5*width_38;
	  if( col+1 == ncol ) xhigh_row = modpos.x() + 0.5*width_38;
	  
	  icell++;
	}
	ysum += width_38;
	nused38 += ncol38_row[super_row];
	nrows38++;

	char prefix[255];
	//Fill out row with Al:
	if( xlow_row > -width_earm/2.0 ){
	  sprintf( prefix, "Albox1_row%d", super_row );
	  G4String boxname = prefix;
	  
	  G4Box *Al_filler1 = new G4Box( boxname, (xlow_row + width_earm/2.0)/2.0, width_38/2.0, depth_leadglass/2.0 );
	  
	  G4String logname = boxname + "log";
	  
	  G4LogicalVolume *Al_filler1_log = new G4LogicalVolume( Al_filler1, GetMaterial("Al"), logname );
	  //Al_filler1_log->SetVisAttributes( Alvisatt );

	  Al_filler1_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
	  
	  G4ThreeVector pos( (xlow_row - width_earm/2.0)/2.0,
			     ysum - 0.5*width_38,
			     depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 );

	  G4String physname = boxname + "phys";
	  
	  new G4PVPlacement( 0, pos, Al_filler1_log, physname, earm_mother_log, false, 0 );
	}
	if( xhigh_row < width_earm/2.0 ){
	  sprintf( prefix, "Albox2_row%d", super_row );
	  G4String boxname = prefix;
	  
	  G4Box *Al_filler2 = new G4Box( boxname, (-xhigh_row + width_earm/2.0)/2.0, width_38/2.0, depth_leadglass/2.0 );
	  
	  G4String logname = boxname + "log";
	  
	  G4LogicalVolume *Al_filler2_log = new G4LogicalVolume( Al_filler2, GetMaterial("Al"), logname );
	  //Al_filler2_log->SetVisAttributes( Alvisatt );

	  Al_filler2_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
	  
	  //G4ThreeVector pos( (xhigh_row + width_earm/2.0)/2.0, modpos.y(), modpos.z() );
	  G4ThreeVector pos( (xhigh_row + width_earm/2.0)/2.0,
			     ysum - 0.5*width_38,
			     depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 );
	  
	  G4String physname = boxname + "phys";
	  
	  new G4PVPlacement( 0, pos, Al_filler2_log, physname, earm_mother_log, false, 0 );
	}
	
      }
    }
  }

  //Fill out top with Al:
  G4Box *top_Al = new G4Box( "top_Al", width_earm/2.0, (-ysum + height_earm/2.0)/2.0, depth_leadglass/2.0 );
  G4LogicalVolume *top_Al_log = new G4LogicalVolume( top_Al, GetMaterial("Al"), "top_Al_log" );
  new G4PVPlacement( 0, G4ThreeVector( 0, (ysum + height_earm/2.0)/2.0, depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 ), top_Al_log, "top_Al_phys", earm_mother_log, false, 0 ); 

  //  top_Al_log->SetVisAttributes( Alvisatt );
  top_Al_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  
  //Next: Put front Aluminum plate in front of ECAL (make wireframe):
  G4Box *ECAL_FrontPlate = new G4Box( "ECAL_FrontPlate", width_earm/2.0, height_earm/2.0, depth_ecal_frontplate/2.0 );
  G4LogicalVolume *ECAL_FrontPlate_log = new G4LogicalVolume( ECAL_FrontPlate, GetMaterial("Al"), "ECAL_FrontPlate_log" );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, zfront_ECAL - depth_ecal_frontplate/2.0 ), ECAL_FrontPlate_log, "ECAL_FrontPlate_phys", earm_mother_log, false, 0 );


  if(fDetCon->fExpType==G4SBS::kGEp){
    //Next: CH2 filter:
    G4Box *CH2_filter = new G4Box( "CH2_filter", width_earm/2.0, height_earm/2.0, depth_CH2/2.0 );
    G4LogicalVolume *CH2_filter_log = new G4LogicalVolume( CH2_filter, GetMaterial("Polyethylene"), "CH2_filter_log" );
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, zfront_ECAL - depth_CDET - depth_CH2/2.0 ), CH2_filter_log, "CH2_filter_phys", earm_mother_log, false, 0 );
    
    G4double z0_CDET = -depth_earm/2.0 + depth_CH2;
    //G4double R0_CDET = R_Earm - depth_leadglass - depth_CDET;
    //G4double R0_CDET = fDist - depth_leadglass - depth_CDET;
    G4double R0_CDET = fDist - depth_CDET;
    
    G4SBSCDet* CDet = new G4SBSCDet(fDetCon);
    CDet->SetR0(R0_CDET);
    CDet->SetZ0(z0_CDET);
    CDet->SetPlanesHOffset(0.0);
    CDet->BuildComponent( earm_mother_log );
    //MakeCDET( earm_mother_log );
    
    G4VisAttributes *CH2_visatt = new G4VisAttributes( G4Colour( 0, 0.6, 0.6 ) );
    CH2_visatt->SetForceWireframe(true);
    
    CH2_filter_log->SetVisAttributes(CH2_visatt);
    
  }
  
  //Visualization:
  G4VisAttributes *FrontPlate_visatt = new G4VisAttributes( G4Colour( 0.7, 0.7, 0.7 ) );
  FrontPlate_visatt->SetForceWireframe(true);
  ECAL_FrontPlate_log->SetVisAttributes( FrontPlate_visatt );

  earm_mother_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  Module_42_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  Module_40_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  Module_38_log->SetVisAttributes( G4VisAttributes::GetInvisible() );

  G4VisAttributes *Mylarvisatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
  Mylarvisatt->SetForceWireframe(true);
  Mylar_wrap_42_log->SetVisAttributes( Mylarvisatt );
  Mylar_wrap_40_log->SetVisAttributes( Mylarvisatt );
  Mylar_wrap_38_log->SetVisAttributes( Mylarvisatt );

  G4VisAttributes *ECALpmtvisatt = new G4VisAttributes( G4Colour( 0, 0, 1 ) );
  ecal_PMT_log->SetVisAttributes( ECALpmtvisatt );
  
}
