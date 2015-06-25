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
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4TwoVector.hh"
#include "G4GenericTrap.hh"
#include "G4Polycone.hh"


#include "G4SBSHArmBuilder.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

G4SBSTargetBuilder::G4SBSTargetBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
  assert(fDetCon);
  fTargLen = 60.0*cm;
  fTargType = kH2;
  fTargDen = 10.5*atmosphere/(300*kelvin*k_Boltzmann);

  fTargPos = G4ThreeVector( 0, 0, 0 );
  fTargDir = G4ThreeVector( 0, 0, 1 );
  
  fSchamFlag = 0;
}

G4SBSTargetBuilder::~G4SBSTargetBuilder(){;}

void G4SBSTargetBuilder::BuildComponent(G4LogicalVolume *worldlog){
  fTargType = fDetCon->fTargType;

  if( fTargType == kLH2 || fTargType == kLD2 ){
    BuildCryoTarget( worldlog ); //The cryotarget is placed entirely inside the scattering chamber vacuum:
  } else {
    BuildGasTarget( worldlog );
  }

  return;

}

void G4SBSTargetBuilder::BuildCryoTarget(G4LogicalVolume *worldlog){

  
  //New version of buildcryotarget updated with scattering chamber for GEP:
  //Start with vacuum snout:
  G4double inch = 2.54*cm;

  //In the following dimensions, "left" and "right" are as viewed from downstream!!!!
  G4double SnoutCenterPlate_width = 9.425*inch;
  G4double SnoutLeftPlate_width = 34.364*inch;
  G4double SnoutRightPlate_width = 27.831*inch;

  G4double SnoutHeight = 37.75*inch;
  G4double SnoutThick = 1.0*inch;

  //Window dimensions:
  G4double SnoutLeftWindow_Rbend_corners = 5.750*inch;
  G4double SnoutLeftWindow_Width = 2.0*9.625*inch;
  G4double SnoutLeftWindow_Height = 2.0*15.875*inch;
  G4double SnoutLeftOffset = -0.744*inch;

  G4double SnoutRightWindow_Rbend_corners = 5.000*inch;
  G4double SnoutRightWindow_Width = 2.0*8.438*inch;
  G4double SnoutRightWindow_Height = 2.0*9.563*inch;
  G4double SnoutRightOffset = -1.4595*inch;

  //Aluminium Dimensions - just need to be bigger than window
  G4double LeftAlHeight = SnoutLeftWindow_Height + 2*inch;
  G4double LeftAlWidth = SnoutLeftWindow_Width + 2*inch;
  G4double LeftAlThick = 0.032*inch;

  G4double RightAlHeight = SnoutRightWindow_Height + 2*inch;
  G4double RightAlWidth = SnoutRightWindow_Width + 2*inch;
  G4double RightAlThick = 0.020*inch;
    
  G4double SnoutLeftWindow_xcenter = 16.438*inch;
  G4double SnoutRightWindow_xcenter = 12.456*inch;

  G4double SnoutLeftWindowAngle = 27.5*deg;
  G4double SnoutRightWindowAngle = 22.0*deg;

  G4double SnoutUpstreamHoleDiameter = 4.870*inch; //All the way through:
  G4double SnoutDownstreamHoleDiameter = 5.010*inch; //To a depth of .26 inches from the front:
  G4double SnoutDownstreamHoleDepth = 0.260*inch;

  G4double SnoutBeamHole_xcenter = 4.591*inch;

  G4double SnoutCenterPlate_xpos = SnoutCenterPlate_width/2.0 - SnoutBeamHole_xcenter;

  //G4ThreeVector TargetCenter_GEP(0,0,16.5*cm);
  G4ThreeVector TargetCenter_GEP(0,0,0);

  //Assume for the moment that this is the distance to the window on the proton side along the line at 22 deg from the target center, which is 16.5 cm downstream of the Hall A origin:
  G4double SnoutRightWindow_R = 41.74*inch;
  G4double SnoutRightWindow_zintercept = SnoutRightWindow_R/cos(SnoutRightWindowAngle) + TargetCenter_GEP.z();
  //The x coordinate of the front corner along the proton side:
  G4double SnoutFrontPlate_xintercept = SnoutCenterPlate_width - SnoutBeamHole_xcenter; //4.834 inches = 12.2784 cm = 
  //The z coordinate of the FRONT of the center plate:
  G4double SnoutFrontPlate_zfront = SnoutRightWindow_zintercept - SnoutFrontPlate_xintercept * tan( SnoutRightWindowAngle );

  //Now we need to make the solids; Start with the center plate:
  G4Box *SnoutCenterPlate_box = new G4Box("SnoutCenterPlate_box", SnoutCenterPlate_width/2.0, SnoutHeight/2.0, SnoutThick/2.0 );

  //Make cylindrical cutouts for beam through-hole:
  G4Tubs *SnoutBeamHoleInner = new G4Tubs("SnoutBeamHoleInner", 0.0, SnoutUpstreamHoleDiameter/2.0, SnoutThick/2.0 + 1.0*mm, 0.0, twopi );
  G4Tubs *SnoutBeamHoleOuter = new G4Tubs("SnoutBeamHoleOuter", 0.0, SnoutDownstreamHoleDiameter/2.0, SnoutDownstreamHoleDepth/2.0 + 1.0*mm, 0.0, twopi );

  //Subtract inner hole:
  G4SubtractionSolid *SnoutCenterPlate_minus_hole1 = new G4SubtractionSolid( "SnoutCenterPlate_minus_hole1", SnoutCenterPlate_box, SnoutBeamHoleInner, 0,
									     G4ThreeVector( -SnoutCenterPlate_width/2.0 + SnoutBeamHole_xcenter, 0, 0 ) );

  //Subtract outer hole:
  //z coordinate should be such that z - (D/2 + 1 mm) + D = T/2 --> z = T/2 - D/2 + 1 mm
  G4SubtractionSolid *SnoutCenterPlate_minus_hole2 = new G4SubtractionSolid( "SnoutCenterPlate_minus_hole2", SnoutCenterPlate_minus_hole1, SnoutBeamHoleOuter, 0,
									     G4ThreeVector( -SnoutCenterPlate_width/2.0 + SnoutBeamHole_xcenter, 0, SnoutThick/2.0 - SnoutDownstreamHoleDepth/2.0 + 1.0*mm ) );

  G4LogicalVolume *SnoutCenterPlate_log = new G4LogicalVolume( SnoutCenterPlate_minus_hole2, GetMaterial("Stainless_Steel"), "SnoutCenterPlate_log" );
  // SnoutCenterPlate_log is placed near Snout Window placement

  //x coordinate should be located such that the beam goes through the center of the hole in the center plate: 


  // **********************************************************************************************
  // Sergey's Fortran Code ------------------------
  // **********************************************************************************************

  // Snout Center Plate
  G4double BLPlateX = SnoutCenterPlate_width;
  G4double BLPlateY = SnoutHeight;
  G4double BLPlateZ = SnoutThick;
  //G4double BLPlateDist = 42.13*inch;  //distance from origin to face of center plate
  G4double BLPlateDist = SnoutFrontPlate_zfront - SnoutThick;

  // Snout Left Plate
  G4double EAPlateX = SnoutLeftPlate_width;
  G4double EAPlateAng = SnoutLeftWindowAngle; 

  // Snout Right Plate
  G4double SBPlateX = SnoutRightPlate_width;
  G4double SBPlateAng = SnoutRightWindowAngle;

  // Scattering Chamber
  G4double ScChRmin = 13.5*inch;
  G4double ScChRmax = 14.5*inch;
  G4double ScChHeight = 37.75*inch; //same as SnoutHeight
  G4double ScChAngle1 = 35.6*deg;   //used to define Scattering Chamber walls, but commented out in Sergey's code
  G4double ScChAngle2 = 42.7*deg;   //used to define Scattering Chamber walls, but commented out in Sergey's code
  
  // Tube bored out of Snout Center Vacuum
  G4double VIV2Length = 0.5*(BLPlateDist - ScChRmax + 5.0*cm);
  G4double VIAngle1 = 27.0*deg;
  G4double VIAngle2 = 31.5*deg;

  G4double BLTube1Rmin = 5.00*cm;
  G4double BLTube1Rmax = 7.00*cm;
  G4double BLTube1Z = 6.5*cm;

  // Tubes defined upstream from Center Snout Plate
  G4double BLTube2Rmin = 5.00*cm;
  G4double BLTube2Rmax = 7.00*cm;
  G4double BLTube2Z = 5.3*inch; 
 
  G4double BLTube3Rmin = 2.415*inch;
  G4double BLTube3Rmax = 2.50*inch;
  G4double BLTube3Z = 5.3*inch; 

  // Convenient Vectors
  G4ThreeVector zero( 0.0, 0.0, 0.0 );      // used for SubtractionSolid
  G4ThreeVector place( 0.0, 0.0, 0.0 );     // used as a placement vector for PVPlacement & SubtractionSolid
  G4ThreeVector unionize( 0.0, 0.0, 0.0 );  // used for UnionSolid

  // Define Rotations
  G4RotationMatrix *ScatChambRot = new G4RotationMatrix;
  ScatChambRot->rotateX( 90.0*deg );

  G4RotationMatrix *SnoutRot2 = new G4RotationMatrix;
  SnoutRot2->rotateY( -SBPlateAng );

  G4RotationMatrix *SnoutRot3 = new G4RotationMatrix;
  SnoutRot3->rotateY( EAPlateAng );

  G4RotationMatrix *RightWedgeRot = new G4RotationMatrix;
  RightWedgeRot->rotateX( 90.0*deg );
  RightWedgeRot->rotateZ( 90.0*deg );

  G4RotationMatrix *LeftWedgeRot = new G4RotationMatrix;
  LeftWedgeRot->rotateX( 90.0*deg );
  LeftWedgeRot->rotateZ( 90.0*deg );
  LeftWedgeRot->rotateZ( EAPlateAng );

  // Define Solids for Scattering Chamber(1) & Snout(3), then union into one Solid

  // Creating the Scattering Chamber 
  G4double tRmin = 0.0*cm;
  G4double tRmax = ScChRmax;       //cm
  G4double tDzz  = 0.5*ScChHeight; //cm
  G4double tSPhi = 0.0*deg;
  G4double tDphi = 360.0*deg;
  G4Tubs *SCCH_tube = new G4Tubs( "SCCH_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);

  // Chamber Wall ----Commented out in Sergey's Code

  G4ThreeVector n1(cos(SBPlateAng), 0.0, -sin(SBPlateAng) );
  G4ThreeVector n2(-sin(SBPlateAng), 0.0, -cos(SBPlateAng) );
  G4ThreeVector xhat(1,0,0);

  G4ThreeVector x1(0.0, 0.0, BLPlateDist); //  plate center face
  G4ThreeVector x2(BLPlateX/2.0, 0.0, SnoutThick); 
  G4ThreeVector x3 = (SnoutRightPlate_width/2.0)*n1;
  G4ThreeVector x4 = SnoutThick*n2;
  G4ThreeVector RightCenter = x1+x2+x3+x4;
  G4double s_temp = (-n1.cross(RightCenter - x1)).mag();
  G4double s_temp1 = (-n1.cross(xhat)).mag();
  G4double S_Right = s_temp / s_temp1;

  G4double t_temp = ((x1-RightCenter).cross(xhat)).mag();
  G4double t_temp1 = (-n1.cross(xhat)).mag();
  G4double T_Right = t_temp / t_temp1;

  G4ThreeVector n1_EA(-cos(EAPlateAng), 0.0, -sin(EAPlateAng) );
  G4ThreeVector n2_EA(sin(EAPlateAng), 0.0, -cos(EAPlateAng) );
  G4ThreeVector x2_EA(-BLPlateX/2.0, 0.0, SnoutThick);
  
  G4ThreeVector LeftCenter = x1 + x2_EA + (SnoutLeftPlate_width/2.0)*n1_EA + SnoutThick*n2_EA;
  G4double s_temp2 = (-n1_EA.cross(LeftCenter-x1)).mag();
  G4double s_temp3 = (-n1_EA.cross(xhat)).mag();
  G4double S_Left = s_temp2 / s_temp3;

  G4double t_temp2 = ((x1-LeftCenter).cross(xhat)).mag();
  G4double t_temp3 = (-n1_EA.cross(xhat)).mag();
  G4double T_Left = t_temp2 / t_temp3;

  G4double TrapOffset = (S_Right-S_Left)/2.0;

  // Creating the Snout - need three Trapezoids (Left, Center, Right) which will be combined via UnionSolid
  // Volume 1 ( RIGHT SIDE ):
  G4double dx1 = 0.5*SnoutRightPlate_width - (BLPlateDist - ScChRmax)*tan(VIAngle1) + 1.47*cm; // 1.47cm yields same dx1 from Sergey's code, and avoids a negative length declaration
  G4double dx2 = 0.5*SnoutRightPlate_width - 0.5*(SnoutRightPlate_width*0.5 - T_Right);
  G4double dy1 = 0.5*ScChHeight;
  G4double dy2 = 0.5*ScChHeight;
  G4double dz = 0.5*(BLPlateDist - ScChRmax);
  G4Trd *VIV1_trap = new G4Trd( "VIV1_trap", dx1, dx2, dy1, dy2, dz );
  
  // Volume 2 ( CENTER ):
  dx1 = 0.5*BLPlateX + (BLPlateDist - ScChRmax)*tan(VIAngle2 - SBPlateAng);
  //dx2 = 0.5*BLPlateX;
  dx2 = 0.5*S_Right+0.5*S_Left;
  dy1 = 0.5*ScChHeight;
  dy2 = 0.5*ScChHeight;
  dz = 0.5*(BLPlateDist - ScChRmax + 5.0*cm);
  G4Trd *VIV2_trap = new G4Trd( "VIV2_trap", dx1, dx2, dy1, dy2, dz );
 
  // Remove cylindrical section from Center Snout piece, replace with Iron tube filled with vacuum
  // This piece is downstream from center snout plate, and placed inside VIV2 (center trapezoid)
  tRmin = 0.0*cm;
  tRmax = BLTube1Rmax;
  tDzz  = 0.5*BLTube1Z;
  G4Tubs *TUB1_sub_tube = new G4Tubs( "TUB1_sub_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  place.set( 0.0, 0.0, VIV2Length - 0.5*BLTube1Z );
  G4SubtractionSolid *BoreHoleInVIV2 = new G4SubtractionSolid( "BoreHoleInVIV2", VIV2_trap, TUB1_sub_tube, 0, place );

  // Bore out center of Iron Tube, fill with vacuum
  tRmax = BLTube1Rmin;
  tDzz += 0.1*cm;
  G4Tubs *TUB1_subtract = new G4Tubs( "TUB1_subtract", tRmin, tRmax, tDzz, tSPhi, tDphi);
  G4SubtractionSolid *TUB1_sub = new G4SubtractionSolid( "TUB1_sub", TUB1_sub_tube, TUB1_subtract, 0, zero );

  tDzz -= 0.1*cm;
  G4Tubs *TUV1_tube = new G4Tubs( "TUV1_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  G4LogicalVolume *TUV1_log = new G4LogicalVolume( TUV1_tube, GetMaterial("Vacuum"), "TUV1_log" );
  G4LogicalVolume *TUB1_log = new G4LogicalVolume( TUB1_sub, GetMaterial("Iron"), "TUB1_log" ); 
  new G4PVPlacement( 0, zero, TUV1_log, "TUB1_Vac", TUB1_log, false, 0 );
  // Place TUB1_log inside the snout - done once the "unionized" Snout LogicalVolume is defined...

  // Volume 3 ( LEFT SIDE ):
  dx1 = 0.5*SnoutLeftPlate_width - (BLPlateDist - ScChRmax)*tan(VIAngle2) + 1.2*cm; // 1.2cm yields same dx1 from Sergey's code, and avoids a negative length declaration
  dx2 = 0.5*SnoutLeftPlate_width - 0.5*(SnoutLeftPlate_width*0.5 - T_Left);
  dy1 = 0.5*ScChHeight;
  dy2 = 0.5*ScChHeight;
  dz = 0.5*(BLPlateDist - ScChRmax);
  G4Trd *VIV3_trap = new G4Trd( "VIV3_trap", dx1, dx2, dy1, dy2, dz );

  // // "Unionize" Volumes 1-3 and Scattering Chamber
  // // Center piece with Right (looking in -z direction)
  // unionize.set( 0.5*(BLPlateX + SnoutRightPlate_width*cos(SBPlateAng) - (BLPlateDist - ScChRmax)*sin(SBPlateAng)), 0.0, 0.0);
  // unionize.setZ( 0.5*(BLPlateDist - ScChRmax + 5.0*cm) - 0.5*(SnoutRightPlate_width*sin(SBPlateAng) + (BLPlateDist - ScChRmax)*cos(SBPlateAng))  ) ;
  // G4UnionSolid *RightSide = new G4UnionSolid( "RightSide", BoreHoleInVIV2, VIV1_trap, SnoutRot2, unionize );

  // //"RightSide" with Left
  // unionize.setX( 0.5*(-BLPlateX - SnoutLeftPlate_width*cos(EAPlateAng) + (BLPlateDist - ScChRmax)*sin(EAPlateAng)) );
  // unionize.setZ( 0.5*(BLPlateDist - ScChRmax + 5.0*cm) - 0.5*(SnoutLeftPlate_width*sin(EAPlateAng) + (BLPlateDist - ScChRmax)*cos(EAPlateAng)) ); 
  // G4UnionSolid *Snout_Union = new G4UnionSolid( "Snout_Union", RightSide, VIV3_trap, SnoutRot3, unionize );

  // // Snout_Union with Scattering Chamber
  // unionize.set(0.0, 0.0, -(BLPlateDist - 0.5*(BLPlateDist - ScChRmax + 5.0*cm)) );
  // G4UnionSolid *Snout_ScChamb_Union = new G4UnionSolid( "Snout_ScChamb_Union", Snout_Union, SCCH_tube, ScatChambRot, unionize );
  // G4LogicalVolume *Snout_ScChamb = new G4LogicalVolume( Snout_ScChamb_Union, GetMaterial("Vacuum"), "Snout_ScChamb" );


  // Start with Scattering Chamber and Center Plate
  G4RotationMatrix *ChamberRot = new G4RotationMatrix;
  ChamberRot->rotateX( 90.0*deg );
  // Rotation makes y' = z now, all unionize vectors will be constructed accordingly
  unionize.set( TrapOffset, BLPlateDist - 0.5*(BLPlateDist - ScChRmax + 5.0*cm), 0.0);
  G4UnionSolid *ScatPlusCenter = new G4UnionSolid( "ScatPlusCenter", SCCH_tube, BoreHoleInVIV2, ChamberRot, unionize );

  // Add on Left and Right Trapezoids
  // Right:
  G4RotationMatrix *RightTrapRot = new G4RotationMatrix;
  RightTrapRot->rotateX( 90.0*deg );
  RightTrapRot->rotateY( -SBPlateAng );
  unionize.setX( 0.5*(BLPlateX + SBPlateX*cos(SBPlateAng) - (BLPlateDist - ScChRmax + 8.0*(0.5*SBPlateX-T_Right))*sin(SBPlateAng)) );
  unionize.setY( BLPlateDist - 0.5*(SBPlateX*sin(SBPlateAng) + (BLPlateDist - ScChRmax - (1/4.0)*(0.5*SBPlateX-T_Right))*cos(SBPlateAng)) );
  unionize.setZ( 0.0 );
  G4UnionSolid *ScatPlusCenterPlusRight = new G4UnionSolid( "ScatPlusCenterPlusRight", ScatPlusCenter, VIV1_trap, RightTrapRot, unionize );
  
  // Left:
  G4RotationMatrix *LeftTrapRot = new G4RotationMatrix;
  LeftTrapRot->rotateX( 90.0*deg );
  LeftTrapRot->rotateY( EAPlateAng );
  // unionize.setX( 0.5*(-BLPlateX - EAPlateX*cos(EAPlateAng) + (BLPlateDist - ScChRmax + 8.0*(0.5*EAPlateX-T_Left))*sin(EAPlateAng)) );
  //unionize.setY( BLPlateDist - 0.5*(EAPlateX*sin(EAPlateAng) + (BLPlateDist - ScChRmax - (1/2.0)*(0.5*EAPlateX-T_Left) )*cos(EAPlateAng)) ); //(1/4.0)*(0.5*EAPlateX-T_Left)
  unionize.setX( 0.5*(-BLPlateX - EAPlateX*cos(EAPlateAng) + (BLPlateDist - ScChRmax + 6.0*(0.5*EAPlateX-T_Left))*sin(EAPlateAng)) );
  unionize.setY( BLPlateDist - 0.5*(EAPlateX*sin(EAPlateAng) + (BLPlateDist - ScChRmax - (1/2.0)*(0.5*EAPlateX-T_Left) )*cos(EAPlateAng)) );
  unionize.setZ( 0.0 );
  G4UnionSolid *Snout_ScChamb_Union = new G4UnionSolid( "Snout_ScChamb_Union", ScatPlusCenterPlusRight, VIV3_trap, LeftTrapRot, unionize );

  G4RotationMatrix *ChamberRot2 = new G4RotationMatrix;
  ChamberRot2->rotateX( -90.0*deg );

  G4LogicalVolume *Snout_ScChamb = new G4LogicalVolume( Snout_ScChamb_Union, GetMaterial("Vacuum"), "Snout_ScChamb" );
  place.set( SnoutCenterPlate_xpos, 0.0, 0.0 ); 
  
  // Place Unionized Scattering Chamber & Snout
  new G4PVPlacement( ChamberRot2, place, Snout_ScChamb, "SnoutandChamber", worldlog, false, 0 );

  // Place Iron tube VIV2 inside CENTER piece
  place.set( 0.0, BLPlateDist - 0.5*BLTube1Z, 0.0);
  new G4PVPlacement( ChamberRot, place, TUB1_log, "TUB1", Snout_ScChamb, false, 0 ); 


  // Making Tube2 (outside of Snout - upstream from snout center plate)
  tRmin = 0.0*cm;
  tRmax = BLTube2Rmax;
  tDzz  = 0.5*BLTube2Z;
  G4Tubs *TUB2_sub_tube = new G4Tubs( "TUB2_sub_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);

  tRmax = BLTube2Rmin;
  tDzz += 0.1*cm;
  G4Tubs *TUB2_subtract = new G4Tubs( "TUB2_subtract", tRmin, tRmax, tDzz, tSPhi, tDphi);
  G4SubtractionSolid *TUB2_sub = new G4SubtractionSolid( "TUB2_sub", TUB2_sub_tube, TUB2_subtract, 0, zero );
  tDzz -= 0.1*cm;
  G4Tubs *TUV2_tube = new G4Tubs( "TUV2_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  G4LogicalVolume *TUV2_log = new G4LogicalVolume( TUV2_tube, GetMaterial("Vacuum"), "TUV2_log" );
  G4LogicalVolume *TUB2_log = new G4LogicalVolume( TUB2_sub, GetMaterial("Iron"), "TUB2_log" ); 
  new G4PVPlacement( 0, zero, TUV2_log, "TUB2_Vac", TUB2_log, false, 0 );
  place.set(0.0, 0.0, BLPlateDist + SnoutThick + 0.5*BLTube2Z);
  new G4PVPlacement( 0, place, TUB2_log, "TUB2", worldlog, false, 0 );

  // Making Tube3
  tRmin = 0.0*cm;
  tRmax = BLTube3Rmax;
  tDzz  = 0.5*BLTube3Z;
  G4Tubs *TUB3_sub_tube = new G4Tubs( "TUB3_sub_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);

  tRmax = BLTube3Rmin;
  tDzz += 0.1*cm;
  G4Tubs *TUB3_subtract = new G4Tubs( "TUB3_subtract", tRmin, tRmax, tDzz, tSPhi, tDphi);
  G4SubtractionSolid *TUB3_sub = new G4SubtractionSolid( "TUB3_sub", TUB3_sub_tube, TUB3_subtract, 0, zero );
  tDzz -= 0.1*cm;
  G4Tubs *TUV3_tube = new G4Tubs( "TUV3_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);
  G4LogicalVolume *TUV3_log = new G4LogicalVolume( TUV3_tube, GetMaterial("Vacuum"), "TUV3_log" );
  G4LogicalVolume *TUB3_log = new G4LogicalVolume( TUB3_sub, GetMaterial("Iron"), "TUB3_log" ); 
  new G4PVPlacement( 0, zero, TUV3_log, "TUB3_Vac", TUB3_log, false, 0 );
  place.set(0.0, 0.0, BLPlateDist + SnoutThick + BLTube2Z + 0.5*BLTube3Z);
  new G4PVPlacement( 0, place, TUB3_log, "TUB3", worldlog, false, 0 );

  // SNOUT PLATES
  // Snout Center Plate (defined at the top) 
  place.set( 0.0, 0.0, BLPlateDist + 0.5*SnoutThick );
  //new G4PVPlacement( 0, place, SnoutCenterPlate_log, "SnoutCenterPlatePhys", worldlog, false, 0 );

  // Dummy variable for SubtractionSolid - need to avoid overlapping surfaces or else undefined behavior should be expected
  G4double subtract = 10*cm; 

  // RIGHT SIDE WINDOW CUT
  // Define a Mother Volume that will house Aluminum + Stainless_Steel Window
  // G4Box *SnoutRightMother_Box = new G4Box( "SnoutRightMother_Box", SnoutRightPlate_width/2.0, SnoutHeight/2.0, (SnoutThick+RightAlThick)/2.0 );
  // G4LogicalVolume *SnoutRightMother = new G4LogicalVolume( SnoutRightMother_Box, GetMaterial("Vacuum"), "SnoutRightMother" );
  
  // // Make Aluminium Window 
  // G4Box *RightAlBox = new G4Box( "RightAlBox", RightAlWidth/2.0, RightAlHeight/2.0, RightAlThick/2.0 );
  // G4LogicalVolume *RightAl_Log = new G4LogicalVolume( RightAlBox, GetMaterial("Aluminum"), "RightAl_Log" );

  // // Define Stainless_Steel box which will be cut by BuildSnoutWindows routine
  // G4Box *SnoutRightPlateBox = new G4Box( "SnoutRightPlateBox", SnoutRightPlate_width/2.0, SnoutHeight/2.0, SnoutThick/2.0 ); 
  // G4LogicalVolume *SnoutRightWindow = 
  //   BuildSnoutWindows( SnoutRightPlateBox, SnoutRightWindow_Width, SnoutRightWindow_Height,
  // 		       subtract, SnoutRightWindow_Rbend_corners, SnoutRightOffset );

  // // Place Al & Window in SnoutRightMother, then place SnoutRightMother
  // new G4PVPlacement( 0, G4ThreeVector(SnoutRightOffset, 0.0, SnoutThick/2.0), RightAl_Log, "RightAluminiumWindow", SnoutRightMother, false, 0);
  // new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, -RightAlThick/2.0), SnoutRightWindow, "SnoutRightWindowPhys", SnoutRightMother, false, 0 );
  // place.setX( SnoutCenterPlate_width/2.0 + (SnoutRightPlate_width/2.0)*cos(SBPlateAng) + ((SnoutThick+RightAlThick)/2.0)*sin(SBPlateAng) );
  // place.setY( 0.0 );
  // place.setZ( BLPlateDist + 0.5*SnoutThick - SnoutRightPlate_width/2.0*sin(SBPlateAng) - SnoutThick/2.0 + ((SnoutThick+RightAlThick)/2.0)*cos(SBPlateAng));
  // new G4PVPlacement( SnoutRot2, place, SnoutRightMother, "SnoutRightMotherPhys", worldlog, false, 0 );
  
  // // LEFT SIDE WINDOW CUT
  // G4Box *SnoutLeftMother_Box = new G4Box( "SnoutLeftMother_Box",  SnoutLeftPlate_width/2.0, SnoutHeight/2.0, (SnoutThick+LeftAlThick)/2.0 );
  // G4LogicalVolume *SnoutLeftMother = new G4LogicalVolume( SnoutLeftMother_Box, GetMaterial("Vacuum"), "SnoutLeftMother" );

  // // Make Aluminium Window 
  // G4Box *LeftAlBox = new G4Box( "LeftAlBox", LeftAlWidth/2.0, LeftAlHeight/2.0, LeftAlThick/2.0 );
  // G4LogicalVolume *LeftAl_Log = new G4LogicalVolume( LeftAlBox, GetMaterial("Aluminum"), "LeftAl_Log" );

  // // Define Stainless_Steel box which will be cut using BuildSnoutWindows routine
  // G4Box *SnoutLeftPlateBox = new G4Box( "SnoutLeftPlateBox", SnoutLeftPlate_width/2.0, SnoutHeight/2.0, SnoutThick/2.0 ); 
  // G4LogicalVolume *SnoutLeftWindow = 
  //   BuildSnoutWindows( SnoutLeftPlateBox, SnoutLeftWindow_Width, SnoutLeftWindow_Height, subtract, SnoutLeftWindow_Rbend_corners, SnoutLeftOffset );
 
  // // Place Al & Window in SnoutLeftMother, and place SnoutLeftMother 
  // new G4PVPlacement( 0, G4ThreeVector(SnoutLeftOffset, 0.0, SnoutThick/2.0), LeftAl_Log, "LeftAluminiumWindow", SnoutLeftMother, false, 0);
  // new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, -LeftAlThick/2.0), SnoutLeftWindow, "SnoutLeftWindowPhys", SnoutLeftMother, false, 0 );
  // place.setX( -SnoutCenterPlate_width/2.0 - (SnoutLeftPlate_width/2.0)*cos(EAPlateAng) - ((SnoutThick + LeftAlThick)/2.0)*sin(EAPlateAng) );
  // place.setY( 0.0 );
  // place.setZ( BLPlateDist + 0.5*SnoutThick - SnoutLeftPlate_width/2.0*sin(EAPlateAng) - SnoutThick/2.0 + ((SnoutThick+LeftAlThick)/2.0)*cos(EAPlateAng));
  // new G4PVPlacement( SnoutRot3, place, SnoutLeftMother, "SnoutLeftMotherPhys", worldlog, false, 0 );

  // // Triangular Wedges place inbetween Left/Right and Center Snout Plates
  // // RIGHT SIDE:
  // tRmin = 0.0*cm;
  // tRmax = SnoutThick;       
  // tDzz  = 0.5*ScChHeight; 
  // tSPhi = 0.0*deg;
  // tDphi = SBPlateAng;
  // G4Tubs *SBWedgeTub = new G4Tubs( "SBWedgeTub", tRmin, tRmax, tDzz, tSPhi, tDphi);
  // G4LogicalVolume *SBWedgeLog = new G4LogicalVolume( SBWedgeTub, GetMaterial("Stainless_Steel"), "SBWedgeLog" );
  // place.set( BLPlateX/2.0, 0.0, BLPlateDist );
  // new G4PVPlacement( RightWedgeRot, place, SBWedgeLog, "SBWedge", worldlog, false, 0 );

  // // LEFT SIDE:
  // tSPhi = 0.0*deg;
  // tDphi = EAPlateAng;
  // G4Tubs *EAWedgeTub = new G4Tubs( "EAWedgeTub", tRmin, tRmax, tDzz, tSPhi, tDphi);
  // G4LogicalVolume *EAWedgeLog = new G4LogicalVolume( EAWedgeTub, GetMaterial("Stainless_Steel"), "EAWedgeLog" );
  // place.set( -BLPlateX/2.0, 0.0, BLPlateDist );
  // new G4PVPlacement( LeftWedgeRot, place, EAWedgeLog, "EAWedge", worldlog, false, 0 );

  // 
  // Making Windows - then make union between windows and center plate
  // Right Side:
  G4Box *SnoutRightPlateBox = new G4Box( "SnoutRightPlateBox", SnoutRightPlate_width/2.0, SnoutHeight/2.0, SnoutThick/2.0 );
  G4Box *SnoutRightPlateWindowSub = new G4Box( "SnoutRightPlateWindowSub", SnoutRightWindow_Width/2.0, SnoutRightWindow_Height/2.0, subtract );
  place.set( SnoutRightOffset, 0.0, 0.0 );
  G4SubtractionSolid *SnoutRightPlateNoWindow = new G4SubtractionSolid( "SnoutRightPlateNoWindow", SnoutRightPlateBox, SnoutRightPlateWindowSub, 0, place );

  // Make a Union between SnoutRightPlateNoWindow & Snout Center Plate
  unionize.setX( SnoutCenterPlate_width/2.0 + (SnoutRightPlate_width/2.0)*cos(SBPlateAng) - (SnoutThick/2.0)*sin(SBPlateAng) );
  unionize.setY( 0.0 );
  unionize.setZ( SnoutThick/2.0 - SnoutRightPlate_width/2.0*sin(SBPlateAng) - (SnoutThick/2.0)*cos(SBPlateAng) );
  G4UnionSolid *RightPlusCenterUnion = new G4UnionSolid( "RightPlusCenterUnion", SnoutCenterPlate_minus_hole2, SnoutRightPlateNoWindow,
							 SnoutRot2, unionize );
	
  // Left Side:
  G4Box *SnoutLeftPlateBox = new G4Box( "SnoutLeftPlateBox", SnoutLeftPlate_width/2.0, SnoutHeight/2.0, SnoutThick/2.0 );
  G4Box *SnoutLeftPlateWindowSub = new G4Box( "SnoutLeftPlateWindowSub", SnoutLeftWindow_Width/2.0, SnoutLeftWindow_Height/2.0, subtract );
  place.set( SnoutLeftOffset, 0.0, 0.0 );
  G4SubtractionSolid *SnoutLeftPlateNoWindow = new G4SubtractionSolid( "SnoutLeftPlateNoWindow", SnoutLeftPlateBox, SnoutLeftPlateWindowSub, 0, place );
  unionize.setX( -SnoutCenterPlate_width/2.0 - (SnoutLeftPlate_width/2.0)*cos(EAPlateAng) + (SnoutThick/2.0)*sin(EAPlateAng) );
  unionize.setY( 0.0 );
  unionize.setZ( SnoutThick/2.0 - SnoutLeftPlate_width/2.0*sin(EAPlateAng) - (SnoutThick/2.0)*cos(EAPlateAng) );
  G4UnionSolid *PlateUnion = new G4UnionSolid( "PlateUnion", RightPlusCenterUnion, SnoutLeftPlateNoWindow,
								 SnoutRot3, unionize );
  G4LogicalVolume *PlateUnionLog = new G4LogicalVolume( PlateUnion, GetMaterial("Stainless_Steel"), "PlateUnionLog" );
  
  // Make Circular Corners which will be placed in Left/Right side Windows 
  // Make a square, then remove a cylinder from the square in order to achieve the desired effect
  // Left:
  G4Box *LeftCornerBox = new G4Box( "LeftCornerBox", SnoutLeftWindow_Rbend_corners/2.0, SnoutLeftWindow_Rbend_corners/2.0, SnoutThick/2.0 );
  tRmin = 0.0*cm;
  tRmax = SnoutLeftWindow_Rbend_corners;
  tDzz = subtract;
  tSPhi = 0.0*deg;
  tDphi = 360.0*deg;
  G4Tubs *LeftTubSub = new G4Tubs( "LeftTubSub", tRmin, tRmax, tDzz, tSPhi, tDphi); 
  place.set( -tRmax/2.0, -tRmax/2.0, 0.0 );
  // Top Left of Left Window
  G4SubtractionSolid *LeftCorner = new G4SubtractionSolid( "LeftCorner", LeftCornerBox, LeftTubSub, 0, place ); 
  G4LogicalVolume *LeftCornerLog = new G4LogicalVolume( LeftCorner, GetMaterial("Stainless_Steel"), "LeftCornerLog" );

  G4ThreeVector SnoutLeftZaxis( -sin( SnoutLeftWindowAngle ), 0, cos( SnoutLeftWindowAngle ) );
  G4ThreeVector SnoutLeftYaxis( 0, 1, 0 );
  G4ThreeVector SnoutLeftXaxis( cos( SnoutLeftWindowAngle ), 0, sin( SnoutLeftWindowAngle ) );
  
  // Right:
  G4Box *RightCornerBox = new G4Box( "RightCornerBox", SnoutRightWindow_Rbend_corners/2.0, SnoutRightWindow_Rbend_corners/2.0, SnoutThick/2.0 );
  tRmin = 0.0*cm;
  tRmax = SnoutRightWindow_Rbend_corners;
  tDzz = subtract;
  tSPhi = 0.0*deg;
  tDphi = 360.0*deg;
  G4Tubs *RightTubSub = new G4Tubs( "RightTubSub", tRmin, tRmax, tDzz, tSPhi, tDphi); 
  place.set( -tRmax/2.0, -tRmax/2.0, 0.0 );
  // Top Left of Right Window
  G4SubtractionSolid *RightCorner = new G4SubtractionSolid( "RightCorner", RightCornerBox, RightTubSub, 0, place ); 
  G4LogicalVolume *RightCornerLog = new G4LogicalVolume( RightCorner, GetMaterial("Stainless_Steel"), "RightCornerLog" );

  G4ThreeVector SnoutRightZaxis( sin( SnoutRightWindowAngle ), 0, cos( SnoutRightWindowAngle ) );
  G4ThreeVector SnoutRightYaxis( 0, 1, 0 );
  G4ThreeVector SnoutRightXaxis( cos( SnoutRightWindowAngle ), 0, -sin( SnoutRightWindowAngle ) );
  
  

  // Placing Corners and Plates	
  // Plates
  place.set( SnoutCenterPlate_xpos, 0.0, BLPlateDist + 0.5*SnoutThick );

  new G4PVPlacement( 0, place, PlateUnionLog, "TestPhys", worldlog, false, 0 );

  //Now get ready to place rounded corners of window cutouts:
  G4ThreeVector FrontRightCornerPos = place + G4ThreeVector( SnoutCenterPlate_width/2.0, 0, SnoutThick/2.0 );
  
  G4ThreeVector RightWindowCutoutCenterPos = FrontRightCornerPos - SnoutThick/2.0 * SnoutRightZaxis + SnoutRightWindow_xcenter * SnoutRightXaxis;

  //Make a box filled with vacuum for the cutout:
  G4Box *RightWindowCutoutVacuum = new G4Box("RightWindowCutoutVacuum", SnoutRightWindow_Width/2.0, SnoutRightWindow_Height/2.0, SnoutThick/2.0 );
  G4LogicalVolume *RightWindowCutoutVacuum_log = new G4LogicalVolume( RightWindowCutoutVacuum, GetMaterial("Vacuum"), "RightWindowCutoutVacuum_log" );
  G4RotationMatrix *rottemp = new G4RotationMatrix;
  rottemp->rotateY( -SnoutRightWindowAngle );
  
  new G4PVPlacement( rottemp, RightWindowCutoutCenterPos, RightWindowCutoutVacuum_log, "RightWindowCutoutVacuum_phys", worldlog, false, 0 );

  //now place corners inside Vacuum volume:
  G4double dxtemp = SnoutRightWindow_Width/2.0 - SnoutRightWindow_Rbend_corners/2.0;
  G4double dytemp = SnoutRightWindow_Height/2.0 - SnoutRightWindow_Rbend_corners/2.0;

  new G4PVPlacement( 0, G4ThreeVector( dxtemp, dytemp, 0 ), RightCornerLog, "RightCornerPhys1", RightWindowCutoutVacuum_log, false, 0 );
  
  rottemp = new G4RotationMatrix;
  rottemp->rotateZ( 90.0*deg );

  new G4PVPlacement( rottemp, G4ThreeVector( dxtemp, -dytemp, 0 ), RightCornerLog, "RightCornerPhys2", RightWindowCutoutVacuum_log, false, 1 );

  rottemp = new G4RotationMatrix;
  rottemp->rotateZ( 180.0*deg );

  new G4PVPlacement( rottemp, G4ThreeVector( -dxtemp, -dytemp, 0 ), RightCornerLog, "RightCornerPhys3", RightWindowCutoutVacuum_log, false, 2 );

  rottemp = new G4RotationMatrix;
  rottemp->rotateZ( -90.0*deg );

  new G4PVPlacement( rottemp, G4ThreeVector( -dxtemp, dytemp, 0 ), RightCornerLog, "RightCornerPhys4", RightWindowCutoutVacuum_log, false, 3 );

  G4ThreeVector FrontLeftCornerPos = place + G4ThreeVector( -SnoutCenterPlate_width/2.0, 0, SnoutThick/2.0 );

  G4ThreeVector LeftWindowCutoutCenterPos = FrontLeftCornerPos - SnoutThick/2.0 * SnoutLeftZaxis - (SnoutLeftPlate_width - SnoutLeftWindow_xcenter)*SnoutLeftXaxis;

  G4Box *LeftWindowCutoutVacuum = new G4Box("LeftWindowCutoutVacuum", SnoutLeftWindow_Width/2.0, SnoutLeftWindow_Height/2.0, SnoutThick/2.0 );
  G4LogicalVolume *LeftWindowCutoutVacuum_log = new G4LogicalVolume( LeftWindowCutoutVacuum, GetMaterial("Vacuum"), "LeftWindowCutoutVacuum_log" );
  rottemp = new G4RotationMatrix;
  rottemp->rotateY( SnoutLeftWindowAngle );

  new G4PVPlacement( rottemp, LeftWindowCutoutCenterPos, LeftWindowCutoutVacuum_log, "LeftWindowCutoutVacuum_phys", worldlog, false, 0 );

  dxtemp = SnoutLeftWindow_Width/2.0 - SnoutLeftWindow_Rbend_corners/2.0;
  dytemp = SnoutLeftWindow_Height/2.0 - SnoutLeftWindow_Rbend_corners/2.0;

  new G4PVPlacement( 0, G4ThreeVector( dxtemp, dytemp, 0 ), LeftCornerLog, "LeftCornerPhys1", LeftWindowCutoutVacuum_log, false, 0 );

  rottemp = new G4RotationMatrix;
  rottemp->rotateZ( 90.0*deg );

  new G4PVPlacement( rottemp, G4ThreeVector( dxtemp, -dytemp, 0 ), LeftCornerLog, "LeftCornerPhys2", LeftWindowCutoutVacuum_log, false, 1 );

  rottemp = new G4RotationMatrix;
  rottemp->rotateZ( 180.0*deg );

  new G4PVPlacement( rottemp, G4ThreeVector( -dxtemp, -dytemp, 0 ), LeftCornerLog, "LeftCornerPhys3", LeftWindowCutoutVacuum_log, false, 2 );

  rottemp = new G4RotationMatrix;
  rottemp->rotateZ( -90.0*deg );

  new G4PVPlacement( rottemp, G4ThreeVector( -dxtemp, dytemp, 0 ), LeftCornerLog, "LeftCornerPhys4", LeftWindowCutoutVacuum_log, false, 3 );  

  //Now let's make a cryotarget:
  G4double Rcell = 4.0*cm;
  G4double uthick = 0.1*mm;
  G4double dthick = 0.15*mm;
  G4double sthick = 0.2*mm;

  G4Tubs *TargetCell = new G4Tubs("TargetCell", 0, Rcell, fTargLen/2.0, 0, twopi );

  G4LogicalVolume *TargetCell_log = new G4LogicalVolume( TargetCell, GetMaterial("LH2"), "TargetCell_log" );

  G4Tubs *TargetWall = new G4Tubs("TargetWall", Rcell, Rcell + sthick, fTargLen/2.0, 0, twopi );

  G4LogicalVolume *TargetWall_log = new G4LogicalVolume( TargetWall, GetMaterial("Al"), "TargetWall_log" );

  G4Tubs *UpstreamWindow = new G4Tubs("UpstreamWindow", 0, Rcell + sthick, uthick/2.0, 0, twopi );
  G4Tubs *DownstreamWindow = new G4Tubs("DownstreamWindow", 0, Rcell + sthick, dthick/2.0, 0, twopi );

  G4LogicalVolume *uwindow_log = new G4LogicalVolume( UpstreamWindow, GetMaterial("Al"), "uwindow_log" );
  G4LogicalVolume *dwindow_log = new G4LogicalVolume( DownstreamWindow, GetMaterial("Al"), "dwindow_log" );

  // G4ThreeVector targ_pos(0,0,
  //Now place everything:
  //Need to fix this later: Union solid defining vacuum chamber needs to be defined with the cylinder as the first solid so that we can place the
  //target as a daughter volume at the origin!
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), TargetCell_log, "TargetCell_phys", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), TargetWall_log, "TargetWall_phys", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,-(fTargLen+uthick)/2.0), uwindow_log, "uwindow_phys", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,(fTargLen+dthick)/2.0 ), dwindow_log, "dwindow_phys", worldlog, false, 0 );
  

  // VISUALS
  // Mother Volumes

  // Vacuum
  TUV1_log->SetVisAttributes( G4VisAttributes::Invisible );
  TUV2_log->SetVisAttributes( G4VisAttributes::Invisible );
  TUV3_log->SetVisAttributes( G4VisAttributes::Invisible );

  G4VisAttributes *Snout_VisAtt = new G4VisAttributes( G4Colour( 0.6, 0.55, 0.65 ) );
  PlateUnionLog->SetVisAttributes( Snout_VisAtt );
  LeftCornerLog->SetVisAttributes( Snout_VisAtt );
  RightCornerLog->SetVisAttributes( Snout_VisAtt );
  RightWindowCutoutVacuum_log->SetVisAttributes( G4VisAttributes::Invisible );
  LeftWindowCutoutVacuum_log->SetVisAttributes( G4VisAttributes::Invisible );

  G4VisAttributes *Snout_VisAttWire = new G4VisAttributes( G4Colour( 0.6, 0.55, 0.65 ) );
  Snout_VisAttWire->SetForceWireframe(true);
  Snout_ScChamb->SetVisAttributes( Snout_VisAttWire );

  G4VisAttributes *AlColor= new G4VisAttributes(G4Colour(0.4,0.4,0.4));
  //LeftAl_Log->SetVisAttributes( AlColor );
  //RightAl_Log->SetVisAttributes( AlColor );

  G4VisAttributes *ironColor= new G4VisAttributes(G4Colour(0.52,0.47,0.47));
  TUB1_log->SetVisAttributes( ironColor );
  TUB2_log->SetVisAttributes( ironColor );
  TUB3_log->SetVisAttributes( ironColor );
}

G4LogicalVolume *G4SBSTargetBuilder::BuildSnoutWindows( G4Box *Plate, G4double x, G4double y, 
							G4double z,   G4double r, G4double offset ) {
  // Subtract out 4 circles, 2 rectangles - all have larger depths than Plate (Left/Right Snout Windows)
  // Use circles to cut corners, creating the desired "bend radius", then remove remaining material with rectangles.
  // x,y,z are conventional global coordinates
  // r is bend radius of Snout Windows
  // Windows are offset from center by some small amount according to CAD drawings

  G4ThreeVector tempvec( offset, 0.0, 0.0 );

  G4Box *box1 = new G4Box( "box1", (x-2*r)/2.0, y/2.0, z/2.0 );
  G4SubtractionSolid *WindowCut1 = new G4SubtractionSolid( "WindowCut1", Plate, box1, 0, tempvec );

  G4Box *box2 = new G4Box( "box2", x/2.0, (y-2*r)/2.0, z/2.0 );
  G4SubtractionSolid *WindowCut2 = new G4SubtractionSolid( "WindowCut2", WindowCut1, box2, 0, tempvec );

  tempvec.set( x/2.0 - r + offset, y/2.0 - r, 0.0 );
  G4Tubs *CircleCut = new G4Tubs( "CircleCut", 0.0, r, z, 0.0, twopi );
  G4SubtractionSolid *WindowCut3 = new G4SubtractionSolid( "WindowCut3", WindowCut2, CircleCut, 0, tempvec );

  tempvec.set( 0.5*x - r + offset, -0.5*y + r, 0.0 );
  G4SubtractionSolid *WindowCut4 = new G4SubtractionSolid( "WindowCut4", WindowCut3, CircleCut, 0, tempvec );

  tempvec.set( -0.5*x + r + offset, 0.5*y - r, 0.0 );
  G4SubtractionSolid *WindowCut5 = new G4SubtractionSolid( "WindowCut5", WindowCut4, CircleCut, 0, tempvec );

  tempvec.set( -0.5*x + r + offset, -0.5*y + r, 0.0 );
  G4SubtractionSolid *WindowCut6 = new G4SubtractionSolid( "WindowCut6", WindowCut5, CircleCut, 0, tempvec );
  G4LogicalVolume *WindowLog = new G4LogicalVolume( WindowCut6, GetMaterial("Stainless_Steel"), "WindowLog" );
  return WindowLog;
}

// void G4SBSTargetBuilder::BuildCryoTarget(G4LogicalVolume *worldlog){
//   //////////////////////////////////////////////////////////////////

//   double snoutclear = 0.45*m;
//   double snout_r = fDetCon->fHArmBuilder->f48D48dist - snoutclear; 

//   G4double entpipe_rin = 31.75*mm;
//   G4double entpipe_rout = entpipe_rin + 0.12*mm;

//   G4double extpipe_rin = 48.00*mm;
  
//   G4double extpipestart = 1.62*m;
//   // 1.62m is where the main exit pipe starts
//   if( snout_r >= extpipestart - 2.0*cm ) snout_r = extpipestart - 2.0*cm;

//   G4double extpipe_len;

//   double sheight = 1.2*m;
  
//   double snoutwallthick   = 1*cm;
//   double swinthick    = 0.38*mm;
//   double swallrad     = 1.143*m/2;
//   double swallrad_in  = 1.041*m/2;
  
//   double hcal_ang_min = -55*deg;
//   //    double hcal_ang_max = -7*deg;
//   //    double hcal_ang_max = -10.5*deg;
//   double hcal_ang_max =  45*deg;
//   double hcal_win_h = 0.4*m;

//   double bb_ang_min = 18*deg;
//   double bb_ang_max = 80*deg;
//   double bb_win_h = 0.5*m;

//   if( bb_ang_min < hcal_ang_max ) bb_ang_min = hcal_ang_max + (swallrad-swallrad_in)/swallrad/4;

//   extpipe_len = extpipestart -  snout_r;

//   G4Tubs *swall = new G4Tubs("scham_wall", swallrad_in, swallrad, sheight/2, 0.*deg, 360.*deg );
  
//   // Cut out for windows

//   G4Tubs *swall_hcalcut = new G4Tubs("scham_wall_hcalcut", swallrad_in-2*cm, swallrad+2*cm, hcal_win_h, hcal_ang_min, hcal_ang_max-hcal_ang_min);
//   G4Tubs *swall_bbcut = new G4Tubs("scham_wall_bbcut", swallrad_in-2*cm, swallrad+2*cm, bb_win_h, bb_ang_min, bb_ang_max-bb_ang_min);

//   G4SubtractionSolid *swallcut = new G4SubtractionSolid("swallcut1", swall, swall_hcalcut);
//   swallcut = new G4SubtractionSolid("swallcut2", swallcut, swall_bbcut);

//   G4Tubs *swall_hcalwin = new G4Tubs("scham_wall_hcalwin", swallrad_in, swallrad_in+swinthick, hcal_win_h, hcal_ang_min, hcal_ang_max-hcal_ang_min);
//   G4Tubs *swall_bbwin = new G4Tubs("scham_wall_bbwin", swallrad_in, swallrad_in+swinthick, bb_win_h, bb_ang_min, bb_ang_max-bb_ang_min);

//   ////    exit pipe
//   //
//   G4Tubs *exttube = new G4Tubs("exitpipetube", extpipe_rin, extpipe_rin+0.120*cm, extpipe_len/2, 0.*deg, 360.*deg );
//   G4Tubs *extvactube = new G4Tubs("exitpipetube_vac", 0.0, extpipe_rin, extpipe_len, 0.*deg, 360.*deg );

//   G4LogicalVolume *extpipe_log = new G4LogicalVolume(exttube, GetMaterial("Aluminum"),"extpipe_log");
//   G4LogicalVolume *extvac_log = new G4LogicalVolume(extvactube, GetMaterial("Vacuum"),"extvac_log");

//   ///////////////////////// SHIELDING //////////////////////////////////////////////////////////////////////////
//   //  FIXME  This should be moved to HArmBuilder
//   double shieldrad = 25.*cm;
//   double gapwidth = 22*cm;
//   double gapheight= 70*cm;
//   /*
//     double shieldlen = 57*cm;

//     G4Tubs *exttubelead = new G4Tubs("exitpipeleadtube", extpipe_rin+0.120*cm, shieldrad, shieldlen/2, 0.*deg, 360.*deg );
//     double cbsize = shieldlen;
//     G4Box *leadclip = new G4Box("leadclip", cbsize, cbsize, cbsize);
//     G4RotationMatrix *cliprm = new G4RotationMatrix();
//     double ang48d48 = fDetCon->fHArmBuilder->f48D48ang;
//     cliprm->rotateY( -ang48d48 );
//     G4SubtractionSolid *extshield = new G4SubtractionSolid("extshield", exttubelead, leadclip, cliprm, 
//     G4ThreeVector(0.0, cbsize*cm*sin(ang48d48) - shieldrad, (cbsize + shieldlen/2 )*cos(ang48d48) - shieldrad*sin(ang48d48) ) );
//     // snout clear 
//     cliprm = new G4RotationMatrix();
//     cliprm->rotateY( hcal_ang_max );
//     double pipeclear = 3.*cm;
//     extshield = new G4SubtractionSolid("extshield2", exttubelead, leadclip, cliprm, 
//     G4ThreeVector( pipeclear + cbsize - sin(hcal_ang_max)*cbsize, 0.0 ,-extpipe_len/2 + cos(hcal_ang_max)*cbsize ) );
//     G4LogicalVolume *extshield_log = new G4LogicalVolume(extshield, GetMaterial("Lead"),"extshield_log");
//   */

//   double shieldlen2 = 29*cm;

//   G4Tubs *exttubelead2 = new G4Tubs("exitpipelead2tube", extpipe_rin+0.120*cm, shieldrad, shieldlen2/2, 0.*deg, 360.*deg );
//   // Mate this with 130cm tall 20x20cm block oriented so it is in the face of the 48d48 magnet at 16.9 deg
//   G4Box *extblocklead2 = new G4Box("extblocklead2", 10*cm, 65*cm, 10*cm );
//   G4RotationMatrix *windowshieldrm = new G4RotationMatrix();
//   windowshieldrm->rotateY(-16.9*deg);
//   G4UnionSolid *exttube_windowshield2 = new G4UnionSolid("exttube_windowshield2", exttubelead2, extblocklead2, windowshieldrm,
// 							 G4ThreeVector(19.0*cm, 0.0, 3.0*cm) );

//   // We now also need blocks for the top and bottom of the window
//   // Window is 70*cm high, but full magnet gap height is 48in
//   double shieldblock3_height = (fDetCon->fHArmBuilder->f48D48depth - gapheight)/2;

//   G4Box *shieldblock3 = new G4Box("shieldblock3", 22*cm/2, shieldblock3_height/2, 10*cm  );

//   G4LogicalVolume *extshield2_log = new G4LogicalVolume(exttube_windowshield2, GetMaterial("Lead"),"extshield2_log");
//   G4LogicalVolume *windowshield_log = new G4LogicalVolume(extblocklead2, GetMaterial("Lead"),"windowshield_log");
//   G4LogicalVolume *shieldblock3_log = new G4LogicalVolume(shieldblock3, GetMaterial("Lead"),"shieldblock3_log");

//   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//   //  Place exit pipe tube

//   //G4Tubs *swall_enthole = new G4Tubs("scham_wall_enthole", 0.0, entpipe_rin, 5.0*cm, 0.*deg, 360.*deg );
//   G4Tubs *swall_enthole = new G4Tubs("scham_wall_enthole", 0.0, entpipe_rout, 5.0*cm, 0.*deg, 360.*deg );
//   G4Tubs *swall_exthole = new G4Tubs("scham_wall_exthole", 0.0, extpipe_rin, 2.0*cm, 0.*deg, 360.*deg );
 
//   G4RotationMatrix *chamholerot = new G4RotationMatrix;
//   chamholerot->rotateY(90.0*deg);

//   //  Cut holes in the scattering chamber
//   G4SubtractionSolid* swall_holes = new G4SubtractionSolid("swall_enthole", swallcut, swall_enthole, chamholerot, G4ThreeVector(-(swallrad+swallrad_in)/2, 0.0, 0.0) );
//   //    swall_holes = new G4SubtractionSolid("swall_holes", swall_holes, swall_exthole, chamholerot, G4ThreeVector((swallrad+swallrad_in)/2, 0.0, 0.0) );

//   //sc_entry_hole is the actual entry port of the scattering chamber
//   G4IntersectionSolid* sc_entry_hole = new G4IntersectionSolid("sc_entry_hole",swallcut, swall_enthole, chamholerot, 
// 							       G4ThreeVector(-(swallrad+swallrad_in)/2, 0.0, 0.0) );

//   G4LogicalVolume *sc_entry_hole_vacuum_log = new G4LogicalVolume( sc_entry_hole, GetMaterial("Vacuum"), "sc_entry_hole_vacuum_log");

//   G4LogicalVolume *swall_log = new G4LogicalVolume(swall_holes, GetMaterial("Aluminum"),"scham_wall_log");

//   G4LogicalVolume *sc_hcalwin_log = new G4LogicalVolume(swall_hcalwin, GetMaterial("Aluminum"),"sc_hcalwin_log");
//   G4LogicalVolume *sc_bbwin_log = new G4LogicalVolume(swall_bbwin, GetMaterial("Aluminum"),"sc_bbwin_log");

//   G4RotationMatrix *schamrot = new G4RotationMatrix;
//   schamrot->rotateX(-90.0*deg);
//   schamrot->rotateZ(-90.0*deg);

//   G4Tubs *chamber_inner = new G4Tubs("chamber_inner", 0.0, swallrad_in,  sheight/2, 0*deg, 360*deg );
//   G4LogicalVolume* chamber_inner_log = new G4LogicalVolume(chamber_inner, GetMaterial("Vacuum"), "cham_inner_log");

//   // Top and bottom
//   G4Tubs *sc_topbottom = new G4Tubs("scham_topbottom", 0.0, swallrad, (swallrad-swallrad_in)/2, 0.*deg, 360.*deg );
//   G4LogicalVolume* sc_topbottom_log = new G4LogicalVolume(sc_topbottom, GetMaterial("Aluminum"), "scham_topbottom_log");

//   //  SNOUT ////////////////////////////////////////////////
    
//   // 0.4m is to give clearance for clamps

//   /*
//     if( fTargType == kLH2 ){
//     // GEp kinematic - 48d48 is 1.6m away
//     snout_r = 1.2*m;
//     }
//     if(fTargType == kLD2 ){
//     // GMn kinematic - 48d48 is 3 m away
//     snout_r = 2.2*m;
//     }
//   */

//   double snoutang_min = hcal_ang_min - (swallrad-swallrad_in)/swallrad/4;
//   double snoutang_max = hcal_ang_max + (swallrad-swallrad_in)/swallrad/4;

//   G4Tubs *snoutbase= new G4Tubs("snoutbase", swallrad, snout_r, hcal_win_h+(swallrad-swallrad_in)/2, -snoutang_max, snoutang_max-snoutang_min );
    
//   G4Tubs *snouthollow= new G4Tubs("snouthollow", swallrad_in, snout_r-snoutwallthick, hcal_win_h, -hcal_ang_max, hcal_ang_max-hcal_ang_min );

//   // Window is nominall 22cm across
//   double hcalwinstart = fDetCon->fHArmBuilder->f48D48ang - gapwidth*1.1/snout_r/2;
//   double hcalwinstop  = fDetCon->fHArmBuilder->f48D48ang + gapwidth*1.1/snout_r/2;

//   G4Tubs *snoutwindowcut= new G4Tubs("snoutwindowcut", snout_r-snoutwallthick-2*cm, snout_r+2*cm, gapheight/2, hcalwinstart, hcalwinstop-hcalwinstart );
//   G4Tubs *snoutwindow  = new G4Tubs("snoutwindow", snout_r-snoutwallthick, snout_r-snoutwallthick+swinthick, gapheight/2, hcalwinstart, hcalwinstop-hcalwinstart );

//   G4SubtractionSolid *snoutsub = new G4SubtractionSolid("snoutsub", snoutbase, snouthollow );
//   snoutsub = new G4SubtractionSolid("snoutsub_beamhole", snoutsub, swall_exthole, chamholerot, G4ThreeVector(snout_r, 0.0, 0.0));
//   snoutsub = new G4SubtractionSolid("snoutsub_beamhole_window", snoutsub, snoutwindowcut, 0, G4ThreeVector(0.0, 0.0, 0.0));
//   //    swall_holes = new G4SubtractionSolid("swall_holes", swall_holes, swall_exthole, chamholerot, G4ThreeVector((swallrad+swallrad_in)/2, 0.0, 0.0) );
//   //


//   G4LogicalVolume* snout_log = new G4LogicalVolume(snoutsub, GetMaterial("Aluminum"), "snout_log");
//   G4LogicalVolume* snoutvacuum_log = new G4LogicalVolume(snouthollow, GetMaterial("Vacuum"), "snoutvacuum_log");
//   G4LogicalVolume* snoutwindow_log = new G4LogicalVolume(snoutwindow, GetMaterial("Aluminum"), "snoutwindow_log");


//   G4RotationMatrix *rm_snout = new G4RotationMatrix();
//   rm_snout->rotateY(90*deg);
//   rm_snout->rotateX(90*deg);

//   //////////////////////////////////////////////////////////

  
//   // Scattering chamber
//   new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), swall_log,
//   		    "scham_wall_phys", worldlog, false, 0);
//   new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), chamber_inner_log,
//   		    "chamber_inner_phys", worldlog, false, 0);

//   new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), sc_bbwin_log,
//   		    "sc_bbwin_phys", worldlog, false, 0);

//   new G4PVPlacement(schamrot, G4ThreeVector(0.0, sheight/2.0 + (swallrad-swallrad_in)/2, 0.0), sc_topbottom_log,
//   		    "scham_top_phys", worldlog, false, 0);
//   new G4PVPlacement(schamrot, G4ThreeVector(0.0, -sheight/2.0 - (swallrad-swallrad_in)/2, 0.0), sc_topbottom_log,
//   		    "scham_bot_phys", worldlog, false, 0);

//   new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), sc_entry_hole_vacuum_log, "sc_entry_hole_vacuum_phys",
//   		    worldlog, false, 0 );

  
//   // **** 6/17/2015 commented out by RFO while importing Sergey's BeamlineBuilder code
//   // Also, should these lines be in BeamlineBuilder??
//   //new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, extpipestart-extpipe_len/2), extpipe_log, "extpipe_phys", worldlog, false, 0);
//   //new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, extpipestart-extpipe_len/2), extvac_log, "extvacpipe_phys", worldlog, false, 0);
//   // **** 6/17/2015


//   if( fDetCon->fExpType == kGEp && fDetCon->fLeadOption == 1 ){
//     //	    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, swallrad+shieldlen/2), extshield_log, "extshield_phys", worldlog, false, 0);
//     new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 129*cm+shieldlen2/2), extshield2_log, "extshield2_phys", worldlog, false, 0);
      
//     G4RotationMatrix *iwindowshieldrm = new G4RotationMatrix( windowshieldrm->inverse() );
//     new G4PVPlacement(windowshieldrm, (*iwindowshieldrm)*G4ThreeVector(19*cm, 0.0, fDetCon->fHArmBuilder->f48D48dist-15*cm), windowshield_log, "windowshield_phys2", worldlog, false, 0);

//     new G4PVPlacement(windowshieldrm, (*iwindowshieldrm)*G4ThreeVector(-2.5*cm, gapheight/2+shieldblock3_height/2, fDetCon->fHArmBuilder->f48D48dist-15*cm), shieldblock3_log, "windowshield_phys3", worldlog, false, 0);
//     new G4PVPlacement(windowshieldrm, (*iwindowshieldrm)*G4ThreeVector(-2.5*cm, -gapheight/2-shieldblock3_height/2, fDetCon->fHArmBuilder->f48D48dist-15*cm), shieldblock3_log, "windowshield_phys4", worldlog, false, 0);
//   }

//   /*
//     new G4PVPlacement(targrot, G4ThreeVector(fTargLen/2.0+downcapthick/2.0, 0.0, 0.0), targ_dcap_log,
//     "targ_dcap_phys", chamber_inner_log, false, 0);
//     new G4PVPlacement(targrot, G4ThreeVector(-fTargLen/2.0-upcapthick/2.0, 0.0, 0.0), targ_ucap_log,
//     "targ_ucap_phys", chamber_inner_log, false, 0);
//   */

//   new G4PVPlacement(rm_snout, G4ThreeVector(0,0,0), snout_log, "snout_phys", worldlog, false, 0);
//   new G4PVPlacement(rm_snout, G4ThreeVector(0,0,0), snoutvacuum_log, "snoutvacuum_phys", worldlog, false, 0);
//   new G4PVPlacement(rm_snout, G4ThreeVector(0,0,0), snoutwindow_log, "snoutwindow_phys", worldlog, false, 0);

//   double wallthick   = 635*um;
//   double upcapthick    = 100*um;
//   double downcapthick  = 125*um;
//   double targconeang = 15.*deg;
//   double cellupradius = 2.0*cm;

//   double celldownradius = fTargLen*sin(targconeang);
//   double cellconelen = fTargLen*cos(targconeang);
//   double cellconeang = atan((celldownradius-cellupradius)/cellconelen);

//   // Aluminum shell sphere
//   G4Sphere *shellsph = new G4Sphere("shellsph", 0, fTargLen, 0, 360.0*deg, 0, targconeang);
//   // Aluminum shell cone
//   G4Cons *shellcon = new G4Cons("shellcon", 0.0, cellupradius, 0.0, celldownradius,  cellconelen/2, 0.0, 360.0*deg);
//   // Union 
//   G4UnionSolid *cryoshell = new G4UnionSolid("cryoshell", shellcon, shellsph, 0, G4ThreeVector(0,0,-cellconelen/2));

//   double cryoupradius =  cellupradius - wallthick/cos(cellconeang);
//   double cryoconelen = cellconelen - upcapthick - downcapthick; 

//   double cryodownradius = cryoupradius + cryoconelen*tan(cellconeang);

//   double cryooffset = cryoconelen/2+upcapthick- cellconelen/2;

//   double cryoang = 14.95*deg;

//   // Cryo sphere
//   G4Sphere *cryosph = new G4Sphere("cryosph", 0, fTargLen-downcapthick, 0*deg, 360.0*deg, 0, cryoang);
//   // Cryo cone
//   G4Cons *cryocon = new G4Cons("shellcon", 0.0, cryoupradius, 0.0, cryodownradius,  cryoconelen/2,  0.0*deg, 360.0*deg);
//   // Union 
//   G4UnionSolid *cryovol1 = new G4UnionSolid("cryovol1", cryocon, cryosph, 0, G4ThreeVector(0.0,0.0,-cryoconelen/2.-upcapthick ));

//   double trimboxsize = 50.0*cm;
//   G4Box *cryotrimbox = new G4Box("cryotrimbox", trimboxsize, trimboxsize, trimboxsize);

//   G4SubtractionSolid *cryovol = new G4SubtractionSolid("cryovol", cryovol1, cryotrimbox, 0, G4ThreeVector(0.0,0.0,-trimboxsize-cryoconelen/2 ));

//   /*
//     targ_tube = new G4Tubs("targ_tube", cellradius-wallthick, cellradius, fTargLen/2.0, 0.*deg, 360.*deg );
//     G4Tubs *targ_ucap = new G4Tubs("targ_ucap", 0.0, cellradius, upcapthick/2.0, 0.*deg, 360.*deg );
//     G4Tubs *targ_dcap = new G4Tubs("targ_dcap", 0.0, cellradius, downcapthick/2.0, 0.*deg, 360.*deg );
//     targ_tube_log = new G4LogicalVolume(targ_tube, Aluminum,"targ_tube_log");
//   */

//   G4LogicalVolume *targ_tube_log = new G4LogicalVolume(cryoshell, GetMaterial("Aluminum"),"targ_tube_log");

//   G4RotationMatrix *targrot = new G4RotationMatrix;
//   targrot->rotateY(-90.0*deg);

//   new G4PVPlacement(targrot, G4ThreeVector(0.0, 0.0, 0.0), targ_tube_log,
//   		      "targ_tube_phys", chamber_inner_log, false, 0);

//   /**/

//   //  G4Tubs *cryo_tube = new G4Tubs("cryo_tube", 0.0, cellradius-wallthick,fTargLen/2.0, 0.*deg, 360.*deg );
//   G4LogicalVolume* cryo_tube_log = NULL;


//   if( fTargType == kLH2 ){
//     cryo_tube_log = new G4LogicalVolume(cryovol, GetMaterial("LH2"), "cryo_tube_log");
//   }
//   if( fTargType == kLD2 ){
//     cryo_tube_log = new G4LogicalVolume(cryovol, GetMaterial("LD2"), "cryo_tube_log");
//   }

//   //    cryo_tube_log = new G4LogicalVolume(cryovol, GetMaterial("Vacuum"), "cryo_tube_vacuum_log");

  
//   new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, cryooffset), cryo_tube_log,
//   		    "cryo_tube_phys", targ_tube_log, false, 0);
  

//   G4VisAttributes *cryoVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
//   if( cryo_tube_log ){
//     cryo_tube_log->SetVisAttributes(cryoVisAtt);
//   }

//   //  Vis attributes
//   chamber_inner_log->SetVisAttributes(G4VisAttributes::Invisible);
//   snoutvacuum_log->SetVisAttributes(G4VisAttributes::Invisible);
//   G4VisAttributes * schamVisAtt
//     = new G4VisAttributes(G4Colour(0.7,0.7,1.0));
//   schamVisAtt->SetForceWireframe(true);
//   //      = new G4VisAttributes(G4VisAttributes::Invisible);
//   swall_log->SetVisAttributes(schamVisAtt);
//   sc_topbottom_log->SetVisAttributes(schamVisAtt);
//   snout_log->SetVisAttributes(schamVisAtt);

//   sc_entry_hole_vacuum_log->SetVisAttributes(schamVisAtt);

//   G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
//   extpipe_log->SetVisAttributes(pipeVisAtt);
//   extvac_log->SetVisAttributes( G4VisAttributes::Invisible );

//   G4VisAttributes *winVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
//   sc_hcalwin_log->SetVisAttributes(winVisAtt);
//   sc_bbwin_log->SetVisAttributes(winVisAtt);
//   snoutwindow_log->SetVisAttributes(winVisAtt);

//   G4VisAttributes *leadVisAtt = new G4VisAttributes(G4Colour(0.15,0.15,0.15));
//   //    extshield_log->SetVisAttributes(leadVisAtt);
//   extshield2_log->SetVisAttributes(leadVisAtt);
//   windowshield_log->SetVisAttributes(leadVisAtt);
//   shieldblock3_log->SetVisAttributes(leadVisAtt);
  
//   return;
// }

void G4SBSTargetBuilder::BuildGasTarget(G4LogicalVolume *worldlog){

  //Desired minimum scattering angle for SBS = 5 deg (for SIDIS @10 deg):
  double sbs_scattering_angle_min = 5.0*deg;
  double sc_exitpipe_radius_inner = 48.0*mm;
  double sc_exitpipe_radius_outer = 50.0*mm;
  double sc_entrypipe_radius_inner = 31.75*mm;
  double sc_entrypipe_radius_outer = sc_entrypipe_radius_inner+0.12*mm;
  double sc_winthick = 0.38*mm;

  double zstart_sc = -1.5*m;
  double zend_sc   = 162.2*cm;

  double zpos_sc = (zstart_sc + zend_sc)/2.0;
  double dz_sc = (zend_sc - zstart_sc)/2.0 - sc_winthick;

  //Material definition was moved to ConstructMaterials();
  //--------- Glass target cell -------------------------------

  double wallthick = 1.61*mm;
  double capthick  = 0.126*mm;
  double cellradius    = 0.75*2.54*cm/2.0;

  double sc_radius_inner = sc_exitpipe_radius_outer;
  double sc_radius_outer = sc_radius_inner + sc_winthick;

  G4Tubs *sc_wall = new G4Tubs("sc_tube", sc_radius_inner, sc_radius_outer, dz_sc, 0.*deg, 360.*deg );
  G4Tubs *sc_vacuum = new G4Tubs("sc_vac", 0.0, sc_radius_inner, dz_sc, 0.*deg, 360.*deg );

  G4Tubs *sc_cap_upstream = new G4Tubs("sc_cap_upstream", sc_entrypipe_radius_inner, sc_radius_outer, sc_winthick/2.0, 0.*deg, 360.*deg );
  G4Tubs *sc_cap_downstream = new G4Tubs("sc_cap_downstream", sc_exitpipe_radius_inner, sc_radius_outer, sc_winthick/2.0, 0.*deg, 360.*deg );
  
  G4LogicalVolume *sc_wall_log = new G4LogicalVolume( sc_wall, GetMaterial("Aluminum"), "sc_wall_log" );
  G4LogicalVolume *sc_vacuum_log = new G4LogicalVolume( sc_vacuum, GetMaterial("Vacuum"), "sc_vacuum_log" );
  G4LogicalVolume *sc_cap_upstream_log = new G4LogicalVolume( sc_cap_upstream, GetMaterial("Aluminum"), "sc_cap_upstream_log" );
  G4LogicalVolume *sc_cap_downstream_log = new G4LogicalVolume( sc_cap_downstream, GetMaterial("Aluminum"), "sc_cap_downstream_log" );
  
  G4Tubs *targ_tube = new G4Tubs("targ_tube", cellradius-wallthick, cellradius, fTargLen/2.0, 0.*deg, 360.*deg );
  G4Tubs *targ_cap = new G4Tubs("targ_cap", 0.0, cellradius, capthick/2.0, 0.*deg, 360.*deg );

  G4LogicalVolume* targ_tube_log = new G4LogicalVolume(targ_tube, GetMaterial("GE180"),"targ_tube_log");
  G4LogicalVolume* targ_cap_log = new G4LogicalVolume(targ_cap, GetMaterial("GE180"),"targ_cap_log");

  // gas
  G4Tubs *gas_tube = new G4Tubs("gas_tube", 0.0, cellradius-wallthick,fTargLen/2.0, 0.*deg, 360.*deg );
  G4LogicalVolume* gas_tube_log = NULL;


  if( fTargType == kH2 || fTargType == kNeutTarg ){
    gas_tube_log = new G4LogicalVolume(gas_tube, GetMaterial("refH2"), "gas_tube_log");
  }
  if( fTargType == k3He ){
    gas_tube_log = new G4LogicalVolume(gas_tube, GetMaterial("pol3He"), "gas_tube_log");
  }

  G4LogicalVolume *motherlog = worldlog;
  double target_zpos = 0.0;
  
  if( fSchamFlag == 1 ){
    motherlog = sc_vacuum_log;
    target_zpos = -zpos_sc;
  }

  //if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, target_zpos), targ_tube_log,
		    "targ_tube_phys", motherlog, false, 0);
  
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, target_zpos+fTargLen/2.0+capthick/2.0), targ_cap_log,
		    "targ_cap_phys1", motherlog, false, 0);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, target_zpos-fTargLen/2.0-capthick/2.0), targ_cap_log,
		    "targ_cap_phys2", motherlog, false, 1);
  
  assert(gas_tube_log);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, target_zpos), gas_tube_log,
		    "gas_tube_phys", motherlog, false, 0);
  
  //Place scattering chamber:
  if( fSchamFlag == 1 ){
    new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, zpos_sc), sc_wall_log, "sc_wall_phys", worldlog, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, zpos_sc-dz_sc-sc_winthick/2.0), sc_cap_upstream_log, "sc_cap_upstream_phys", worldlog, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, zpos_sc+dz_sc+sc_winthick/2.0), sc_cap_downstream_log, "sc_cap_downstream_phys", worldlog, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, zpos_sc), sc_vacuum_log, "sc_vacuum_phys", worldlog, false, 0 );
  }
  
  sc_vacuum_log->SetVisAttributes( G4VisAttributes::Invisible );

  //Visualization attributes:
  //ScatteringChamber_log->SetVisAttributes( G4VisAttributes::Invisible );

  G4VisAttributes *sc_wall_visatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
  sc_wall_visatt->SetForceWireframe(true);

  sc_wall_log->SetVisAttributes( sc_wall_visatt );
  sc_cap_upstream_log->SetVisAttributes( sc_wall_visatt );
  sc_cap_upstream_log->SetVisAttributes( sc_wall_visatt );

  // G4VisAttributes *sc_exit_pipe_visatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
  // sc_exit_pipe_log->SetVisAttributes( sc_exit_pipe_visatt );

  // G4VisAttributes *sc_win_visatt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0.0 ) );
  // sc_window_sbs_log->SetVisAttributes( sc_win_visatt );
  // //sc_window_bb_log->SetVisAttributes( sc_win_visatt );

  G4VisAttributes *tgt_cell_visatt = new G4VisAttributes( G4Colour( 1.0, 1.0, 1.0 ) );
  tgt_cell_visatt->SetForceWireframe(true);

  targ_cap_log->SetVisAttributes( tgt_cell_visatt );
  targ_tube_log->SetVisAttributes( tgt_cell_visatt );

  G4VisAttributes *tgt_gas_visatt = new G4VisAttributes( G4Colour( 0.0, 1.0, 1.0 ) );
  gas_tube_log->SetVisAttributes( tgt_gas_visatt );

}
