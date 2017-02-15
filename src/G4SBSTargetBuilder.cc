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
#include "G4ExtrudedSolid.hh"

#include "G4SBSCalSD.hh"

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

  fFlux = false;
  
  fSchamFlag = 0;
}

G4SBSTargetBuilder::~G4SBSTargetBuilder(){;}

void G4SBSTargetBuilder::BuildComponent(G4LogicalVolume *worldlog){
  fTargType = fDetCon->fTargType;
  
  // EFuchey 2017/02/10: organized better this with a switch instead of an endless chain of if...  else...
  if( (fTargType == kLH2 || fTargType == kLD2) ){
    switch(fDetCon->fExpType){
    case(kGEp):
      BuildGEpScatCham( worldlog );
      break;
    case(kC16):
      BuildC16ScatCham( worldlog );
      break;
    default:
      //BuildGEpScatCham( worldlog );
      BuildStandardScatCham( worldlog );
      break;
    }
  } else {
    BuildGasTarget( worldlog );
  }
  
  return;

}

// EFuchey: 2017/02/10: Making a standard function to build the cryotarget itself.
// The code for building C16 and GEp are indeed almost identical.
void G4SBSTargetBuilder::BuildStandardCryoTarget(G4LogicalVolume *motherlog, 
						 G4RotationMatrix *rot_targ, G4ThreeVector targ_offset){
  // Now let's make a cryotarget:
  G4double Rcell = 4.0*cm;
  G4double uthick = 0.1*mm;
  G4double dthick = 0.15*mm;
  G4double sthick = 0.2*mm;
  
  G4Tubs *TargetMother_solid = new G4Tubs( "TargetMother_solid", 0, Rcell + sthick, (fTargLen+uthick+dthick)/2.0, 0.0, twopi );
  G4LogicalVolume *TargetMother_log = new G4LogicalVolume( TargetMother_solid, GetMaterial("Vacuum"), "TargetMother_log" );
  
  G4Tubs *TargetCell = new G4Tubs( "TargetCell", 0, Rcell, fTargLen/2.0, 0, twopi );
  
  G4LogicalVolume *TargetCell_log;
  
  if( fTargType == kLH2 ){
    TargetCell_log = new G4LogicalVolume( TargetCell, GetMaterial("LH2"), "TargetCell_log" );
  } else {
    TargetCell_log = new G4LogicalVolume( TargetCell, GetMaterial("LD2"), "TargetCell_log" );
  }
  
  G4Tubs *TargetWall = new G4Tubs("TargetWall", Rcell, Rcell + sthick, fTargLen/2.0, 0, twopi );
  
  G4LogicalVolume *TargetWall_log = new G4LogicalVolume( TargetWall, GetMaterial("Al"), "TargetWall_log" );
  
  G4Tubs *UpstreamWindow = new G4Tubs("UpstreamWindow", 0, Rcell + sthick, uthick/2.0, 0, twopi );
  G4Tubs *DownstreamWindow = new G4Tubs("DownstreamWindow", 0, Rcell + sthick, dthick/2.0, 0, twopi );
  
  G4LogicalVolume *uwindow_log = new G4LogicalVolume( UpstreamWindow, GetMaterial("Al"), "uwindow_log" );
  G4LogicalVolume *dwindow_log = new G4LogicalVolume( DownstreamWindow, GetMaterial("Al"), "dwindow_log" );
  
  // Now place everything:
  // Need to fix this later: Union solid defining vacuum chamber 
  // needs to be defined with the cylinder as the first solid 
  // so that we can place the target as a daughter volume at the origin!
  
  G4double ztemp = -(fTargLen+uthick+dthick)/2.0;
  // Place upstream window:
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+uthick/2.0), uwindow_log, "uwindow_phys", TargetMother_log, false, 0 );
  // Place target and side walls:
  ztemp += uthick;
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+fTargLen/2.0), TargetCell_log, "TargetCell_phys", TargetMother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+fTargLen/2.0), TargetWall_log, "TargetWall_phys", TargetMother_log, false, 0 );
  ztemp += fTargLen;
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+dthick/2.0), dwindow_log, "dwindow_phys", TargetMother_log, false, 0 );
  
  G4double targ_zcenter = (uthick-dthick)/2.0; //position of target center relative to target mother volume
   
  //Compute position of target relative to scattering chamber:
  //The target center should be at
  //G4RotationMatrix *rot_temp = new G4RotationMatrix;
  //rot_temp = new G4RotationMatrix();
  //rot_temp->rotateX(90.0*deg);

  //for z of target center to be at zero, 
  G4double temp = targ_offset.y();
  targ_offset.setY(temp+targ_zcenter);
  
  new G4PVPlacement( rot_targ, targ_offset, TargetMother_log, "TargetMother_phys", motherlog, false, 0 );
  
  if( fFlux ){ //Make a sphere to compute particle flux:
    G4Sphere *fsph = new G4Sphere( "fsph", 1.5*fTargLen/2.0, 1.5*fTargLen/2.0+cm, 0.0*deg, 360.*deg,
				   0.*deg, 150.*deg );
    G4LogicalVolume *fsph_log = new G4LogicalVolume( fsph, GetMaterial("Air"), "fsph_log" );
    new G4PVPlacement( rot_targ, targ_offset, fsph_log, "fsph_phys", motherlog, false, 0 );
    
    G4String FluxSDname = "FLUX";
    G4String Fluxcollname = "FLUXHitsCollection";
    G4SBSCalSD *FluxSD = NULL;
    if( !( FluxSD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(FluxSDname) ) ){
      G4cout << "Adding FLUX SD to SDman..." << G4endl;
      FluxSD = new G4SBSCalSD( FluxSDname, Fluxcollname );
      fDetCon->fSDman->AddNewDetector( FluxSD );
      (fDetCon->SDlist).insert( FluxSDname );
      fDetCon->SDtype[FluxSDname] = kCAL;
      
      (FluxSD->detmap).depth = 0;
    }
    fsph_log->SetSensitiveDetector( FluxSD );
    fsph_log->SetUserLimits(  new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }

  
  G4VisAttributes *Targ_visatt = new G4VisAttributes( G4Colour( 0.1, 0.05, 0.9 ) );
  TargetCell_log->SetVisAttributes( Targ_visatt );

  G4VisAttributes *TargWall_visatt = new G4VisAttributes( G4Colour( 0.9, .05, 0.1 ) );
  TargWall_visatt->SetForceWireframe( true );
  TargetWall_log->SetVisAttributes( TargWall_visatt );
  uwindow_log->SetVisAttributes( TargWall_visatt );
  dwindow_log->SetVisAttributes( TargWall_visatt );
  TargetMother_log->SetVisAttributes( G4VisAttributes::Invisible );
}

//This function is meant to build the "Standard" scattering chamber for GMn
void G4SBSTargetBuilder::BuildStandardScatCham(G4LogicalVolume *worldlog ){
  G4double inch = 2.54*cm;
  
  G4LogicalVolume *logicScatChamberTank =0;
  G4LogicalVolume *logicScatChamberExitFlangePlate =0;
  G4LogicalVolume *logicScatChamberFrontClamshell =0;
  G4LogicalVolume *logicScatChamberBackClamshell =0;
  G4LogicalVolume *logicScatChamber =0;

  // Scattering chamber tank:
  // basic volume:
  G4double SCHeight = 44.75*inch;
  G4double SCRadius = 20.0*inch;
  G4double SCTankThickness = 2.5*inch;
  G4double SCTankRadius = SCRadius+SCTankThickness;
  G4double SCTankHeight = SCHeight;
  G4double SCOffset = 3.75*inch;
  
  G4Tubs* solidSCTank_0 = 
    new G4Tubs("SCTank_0", SCRadius, SCTankRadius, 0.5*SCTankHeight, 0.0*deg, 360.0*deg);
	       
  
  // exit flange:
  G4double SCExitFlangePlateHLength = 22.5*sin(25.5*atan(1)/45.0)*inch;
  G4double SCExitFlangePlateHeight = 11.0*inch;
  G4double SCExitFlangePlateThick = 1.25*inch;
  G4double SCExitFlangeHAngleApert = atan(SCExitFlangePlateHLength/(SCTankRadius+SCExitFlangePlateThick));
  G4double SCExitFlangeMaxRad = SCExitFlangePlateHLength/sin(SCExitFlangeHAngleApert);
  
  G4Tubs* solidSCExitFlangetubs = 
    new G4Tubs("SCExFlange_tubs", SCRadius, SCExitFlangeMaxRad, 
	       0.5*SCExitFlangePlateHeight, 0.0, 2.0*SCExitFlangeHAngleApert);
  
  G4RotationMatrix* rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(180.0*deg+SCExitFlangeHAngleApert);
  
  G4UnionSolid* solidSCTank_0_exft = 
    new G4UnionSolid("solidSCTank_0_exft", solidSCTank_0, solidSCExitFlangetubs,
		     rot_temp, G4ThreeVector(0,0,SCOffset));
  
  G4Box* ExitFlangeHeadCut = new G4Box("ExitFlangeHeadCut", 0.5*m, 0.5*m, 0.5*m); 
  
  
  G4SubtractionSolid* solidSCTank_0_exf = 
    new G4SubtractionSolid("solidSCTank_0_exf", solidSCTank_0_exft, ExitFlangeHeadCut,
			   0, G4ThreeVector(-SCTankRadius-0.5*m,0,0));
  
  // exit flange hole:
  G4double SCExitFlangeHoleHeight = 7.85*inch;
  G4double SCExitFlangeHoleAngleApert = 38.25*deg;
  
  G4Tubs* solidSCExFH = 
    new G4Tubs("SCExFH", SCRadius-1.0*cm, SCTankRadius+1.5*inch,
	       0.5*SCExitFlangeHoleHeight, 0.0, SCExitFlangeHoleAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(180.0*deg+SCExitFlangeHoleAngleApert*0.5);
  
  G4SubtractionSolid* solidSCTank_0_exfh = 
    new G4SubtractionSolid("solidSCTank_0_exfh", solidSCTank_0_exf, solidSCExFH,
			   rot_temp, G4ThreeVector(0,0,SCOffset));
  
  // windows holes: 
  G4double SCWindowHeight = 18.0*inch;
  G4double SCWindowAngleApert = 149.0*deg;
  G4double SCWindowAngleOffset = 11.0*deg;
  
  G4Tubs* solidSCWindow = 
    new G4Tubs("SCWindow", SCRadius-1.0*cm, SCTankRadius+1.0*cm, 0.5*SCWindowHeight, 0.0, SCWindowAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(90.0*deg+SCWindowAngleApert*0.5-SCWindowAngleOffset);
  
  G4SubtractionSolid* solidSCTank_0_wf = 
    new G4SubtractionSolid("solidSCTank_0_wf", solidSCTank_0_exfh, solidSCWindow,
			   rot_temp, G4ThreeVector(0,0,SCOffset));
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(-90.0*deg+SCWindowAngleApert*0.5+SCWindowAngleOffset);
  
  G4SubtractionSolid* solidSCTank_0_wb = 
    new G4SubtractionSolid("solidSCTank_0_wb", solidSCTank_0_wf, solidSCWindow,
  			   rot_temp, G4ThreeVector(0,0,SCOffset));
  
  //Solid Scattering chamber tank
  G4VSolid* solidScatChamberTank = solidSCTank_0_wb;
  
  //Logic scat chamber tank
  logicScatChamberTank = 
    new G4LogicalVolume(solidScatChamberTank, GetMaterial("Aluminum"), "ScatChamberTank_log");
  
  //Scattering chamber tank placement:
  G4RotationMatrix* rotSC = new G4RotationMatrix();
  rotSC->rotateX(-90.0*deg);
  
  G4ThreeVector* SCPlacement = new G4ThreeVector(0,0,-SCOffset);
  SCPlacement->rotateX(90*deg);
  
  new G4PVPlacement(rotSC, *SCPlacement, logicScatChamberTank, "ScatChamberTankPhys", worldlog, false, 0);
  
  //Exit Flange Plate
  G4Box* solidSCExitFlangePlate = 
    new G4Box("SCExitFlangePlate_sol", SCExitFlangePlateThick/2.0, 
	      SCExitFlangePlateHeight/2.0, SCExitFlangePlateHLength); 
  
  logicScatChamberExitFlangePlate = 
    new G4LogicalVolume(solidSCExitFlangePlate, GetMaterial("Aluminum"), "SCExitFlangePlate_log");

  new G4PVPlacement(0, G4ThreeVector(-SCTankRadius-SCExitFlangePlateThick/2.0,0,0), 
		    logicScatChamberExitFlangePlate, "SCExitFlangePlate", worldlog, false, 0); 
  
  //Front and back Clamshells...
  //Basic solid: 
  G4double SCClamHeight = 20.0*inch;
  G4double SCClamThick = 1.25*inch;
  G4double SCClamAngleApert = 151.0*deg;
  
  G4Tubs* solidSCClamshell_0 = 
    new G4Tubs("solidSCClamshell_0", SCTankRadius, SCTankRadius+SCClamThick, 
	       0.5*SCClamHeight, 0.0, SCClamAngleApert);
  
  // Front Clamshell:
  G4double SCFrontClamOuterRadius = SCTankRadius+SCClamThick;
  G4double SCBeamExitAngleOffset = 64.5*deg;
  
  // Exit Beam Pipe
  G4Tubs* solidExitBeamPipe = new G4Tubs("solidExitBeamPipe", 
					 0.0, 55.0*mm, 2.100*inch, 0.0, 360.0*deg);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateY(-90.0*deg+SCBeamExitAngleOffset);
  
  G4UnionSolid* solidSCFrontClam_0_ebp =
    new G4UnionSolid("solidSCClamshell_0_ebp", solidSCClamshell_0, 
		     solidExitBeamPipe, rot_temp, 
		     G4ThreeVector((SCFrontClamOuterRadius+2.053*inch)*sin(90.0*deg-SCBeamExitAngleOffset),
				   (SCFrontClamOuterRadius+2.053*inch)*cos(90.0*deg-SCBeamExitAngleOffset), 
				   0.0)
		     );
  
  G4Tubs* solidExitBeamPipeFlange = new G4Tubs("solidExitBeamPipeFlange", 
					       0.0, 2.985*inch, 0.3925*inch, 0.0, 360.0*deg);
  
  G4UnionSolid* solidSCFrontClam_0_ebpf =
    new G4UnionSolid("solidSCClamshell_0_ebpf", solidSCFrontClam_0_ebp, 
		     solidExitBeamPipeFlange, rot_temp, 
		     G4ThreeVector((SCFrontClamOuterRadius+3.7605*inch)*sin(90.0*deg-SCBeamExitAngleOffset), 
				   (SCFrontClamOuterRadius+3.7605*inch)*cos(90.0*deg-SCBeamExitAngleOffset), 
				   0.0));
  
  G4Tubs* solidExitBeamPipeHole = new G4Tubs("solidBackViewPipeHole", 
					     0.0, 50.0*mm, 4.0*inch, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidSCFrontClam_0_ebph =
    new G4SubtractionSolid("solidSCClamshell_0_ebph", solidSCFrontClam_0_ebpf, 
			   solidExitBeamPipeHole, rot_temp, 
			   G4ThreeVector(SCFrontClamOuterRadius*sin(90.0*deg-SCBeamExitAngleOffset),
					 SCFrontClamOuterRadius*cos(90.0*deg-SCBeamExitAngleOffset), 
					 0.0)
			   );
  
  // ICI
  
  G4VSolid* solidSCFrontClamShell = solidSCFrontClam_0_ebph;
  
  logicScatChamberFrontClamshell = 
    new G4LogicalVolume(solidSCFrontClamShell, GetMaterial("Aluminum"), "SCFrontClamshell_log");
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateZ(90.0*deg+SCClamAngleApert*0.5-SCWindowAngleOffset);
  
  new G4PVPlacement(rot_temp, G4ThreeVector(0,0,0), logicScatChamberFrontClamshell, 
  		    "SCFrontClamshell", worldlog, false, 0);
  
  //Back Clamshell:
  G4double SCBackClamThick = 0.80*inch;
  G4double SCBackClamHeight = 12.0*inch;
  G4double SCBackClamAngleApert = 145.7*deg;
  G4double SCBackClamAngleOffset = -2.5*deg;
  G4double SCBackClamOuterRadius = SCTankRadius+SCClamThick+SCBackClamThick;
  G4double SCBeamEntranceAngleOffset = SCBackClamAngleOffset-84*deg;
  G4double SCBackViewPipeAngleOffset = SCBackClamAngleOffset-39.9*deg;
  
  G4Tubs* solidSCAddBackClam = 
    new G4Tubs("solidSCAddBackClam", SCTankRadius+SCClamThick-0.5*cm, SCBackClamOuterRadius, 
	       0.5*SCBackClamHeight, 0.0, SCBackClamAngleApert);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateZ(SCBackClamAngleOffset);
  
  G4UnionSolid* solidSCBackClam_0_abc = 
    new G4UnionSolid("solidSCClamshell_0_abc", solidSCClamshell_0, solidSCAddBackClam,
		     rot_temp, G4ThreeVector(0,0,0));
  
  //Entrance beam pipe
  G4Tubs* solidEntranceBeamPipe = 
    new G4Tubs("solidEntranceBeamPipe", 0.0, 2.25*inch, 0.825*inch, 0.0, 360.0*deg);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateY(-90.0*deg-SCBeamEntranceAngleOffset);
  
  G4UnionSolid* solidSCBackClam_0_ebp =
    new G4UnionSolid("solidSCClamshell_0_ebp", solidSCBackClam_0_abc, solidEntranceBeamPipe, rot_temp, 
		     G4ThreeVector(SCBackClamOuterRadius*sin(90.0*deg+SCBeamEntranceAngleOffset),
				   SCBackClamOuterRadius*cos(90.0*deg+SCBeamEntranceAngleOffset),
				   0.0)
		     );
  
  G4Tubs* solidEntranceBeamPipeHole = 
    new G4Tubs("solidEntranceBeamPipeHole", 0.0, 1.0*inch, 4.0*inch, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidSCBackClam_0_ebph =
    new G4SubtractionSolid("solidSCClamshell_0_ebph", solidSCBackClam_0_ebp, 
			   solidEntranceBeamPipeHole, rot_temp, 
			   G4ThreeVector(SCBackClamOuterRadius*
					 sin(90.0*deg+SCBeamEntranceAngleOffset),
					 SCBackClamOuterRadius*
					 cos(90.0*deg+SCBeamEntranceAngleOffset),
					 0.0)
			   );
  
  //Back view pipe
  G4Tubs* solidBackViewPipe = new G4Tubs("solidBackViewPipe", 
					 0.0, 2.12*inch, 1.555*inch, 0.0, 360.0*deg);
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateY(-90.0*deg-SCBackViewPipeAngleOffset);
  
  G4UnionSolid* solidSCBackClam_0_bvp =
    new G4UnionSolid("solidSCClamshell_0_bvp", solidSCBackClam_0_ebph, 
		     solidBackViewPipe, rot_temp, 
		     G4ThreeVector((SCBackClamOuterRadius+1.501*inch)*sin(90.0*deg+SCBackViewPipeAngleOffset),
				   (SCBackClamOuterRadius+1.501*inch)*cos(90.0*deg+SCBackViewPipeAngleOffset), 
				   0.0)
		     );
  
  G4Tubs* solidBackViewPipeFlange = new G4Tubs("solidBackViewPipeFlange", 
					       0.0, 3.0*inch, 0.5*inch, 0.0, 360.0*deg);
  
  G4UnionSolid* solidSCBackClam_0_bvpf =
    new G4UnionSolid("solidSCClamshell_0_bvpf", solidSCBackClam_0_bvp, 
		     solidBackViewPipeFlange, rot_temp, 
		     G4ThreeVector((SCBackClamOuterRadius+2.556*inch)*sin(90.0*deg+SCBackViewPipeAngleOffset), 
				   (SCBackClamOuterRadius+2.556*inch)*cos(90.0*deg+SCBackViewPipeAngleOffset), 
				   0.0));
  
  G4Tubs* solidBackViewPipeHole = new G4Tubs("solidBackViewPipeHole", 
					     0.0, 2.0*inch, 4.0*inch, 0.0, 360.0*deg);
  
  G4SubtractionSolid* solidSCBackClam_0_bvph =
    new G4SubtractionSolid("solidSCClamshell_0_bvph", solidSCBackClam_0_bvpf, 
			   solidBackViewPipeHole, rot_temp, 
			   G4ThreeVector(SCBackClamOuterRadius*sin(90.0*deg+SCBackViewPipeAngleOffset),
					 SCBackClamOuterRadius*cos(90.0*deg+SCBackViewPipeAngleOffset), 
					 0.0)
			   );
  
  
  G4VSolid* solidSCBackClamShell = solidSCBackClam_0_bvph;
    
  logicScatChamberBackClamshell = new G4LogicalVolume(solidSCBackClamShell, 
						      GetMaterial("Aluminum"), 
						      "SCBackClamshell_log");
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  rot_temp->rotateZ(-90.0*deg+SCClamAngleApert*0.5+SCWindowAngleOffset);
  
  new G4PVPlacement(rot_temp, G4ThreeVector(0,0,0), logicScatChamberBackClamshell, 
  		    "SCBackClamshell", worldlog, false, 0);
  
  // Scattering chamber volume
  //
  G4Tubs* solidScatChamber = new G4Tubs("SC", 0.0, SCRadius, 0.5* SCHeight, 0.0*deg, 360.0*deg);
			       
  
  logicScatChamber = new G4LogicalVolume(solidScatChamber, GetMaterial("Vacuum"), "ScatChamber_log");
  
  new G4PVPlacement(rotSC, *SCPlacement, logicScatChamber, "ScatChamberPhys", worldlog, false, 0);
  
  //G4double YOffset_std = 0.0;
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(90.0*deg);
  
  //Call BuildStandardCryoTarget HERE !
  BuildStandardCryoTarget(logicScatChamber, rot_temp, G4ThreeVector(0, 0, SCOffset));
  
  G4VisAttributes* colourGrey= new G4VisAttributes(G4Colour(0.5,0.5,0.5)); 
  colourGrey->SetVisibility(true);
  colourGrey->SetForceWireframe(true);
  
  G4VisAttributes* Invisible = new G4VisAttributes(G4Colour(0.,0.,0.)); 
  Invisible->SetVisibility(false);

  G4VisAttributes* colourCyan= new G4VisAttributes(G4Colour(0.,1.,1.)); 
  colourCyan->SetVisibility(true);
  //colourCyan->SetForceWireframe(true);
  
  logicScatChamberTank->SetVisAttributes(colourGrey);
  logicScatChamberFrontClamshell->SetVisAttributes(colourCyan);
  logicScatChamberBackClamshell->SetVisAttributes(colourGrey);
  logicScatChamberExitFlangePlate->SetVisAttributes(colourCyan);
  logicScatChamber->SetVisAttributes(Invisible);
  
}

// This function will replaced BuildCryoTarget
// It is identical except it uses the function BuildStandardCryoTarget.
// It is also used for GEp only now, and not for GMn
void G4SBSTargetBuilder::BuildGEpScatCham(G4LogicalVolume *worldlog ){
  // GEP scattering chamber:
  
  // This shall go with Target construction now
  // if( fFlux ){ //Make a sphere to compute particle flux:
  //   G4Sphere *fsph = new G4Sphere( "fsph", 1.5*fTargLen/2.0, 1.5*fTargLen/2.0+cm, 0.0*deg, 360.*deg,
  // 				   0.*deg, 150.*deg );
  //   G4LogicalVolume *fsph_log = new G4LogicalVolume( fsph, GetMaterial("Air"), "fsph_log" );

  //   fsph_log->SetUserLimits(  new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
    
  //   new G4PVPlacement( 0, G4ThreeVector(0,0,0), fsph_log, "fsph_phys", worldlog, false, 0 );

  //   G4String FluxSDname = "FLUX";
  //   G4String Fluxcollname = "FLUXHitsCollection";
  //   G4SBSCalSD *FluxSD = NULL;
  //   if( !( FluxSD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(FluxSDname) ) ){
  //     G4cout << "Adding FLUX SD to SDman..." << G4endl;
  //     FluxSD = new G4SBSCalSD( FluxSDname, Fluxcollname );
  //     fDetCon->fSDman->AddNewDetector( FluxSD );
  //     (fDetCon->SDlist).insert( FluxSDname );
  //     fDetCon->SDtype[FluxSDname] = kCAL;

  //     (FluxSD->detmap).depth = 0;
  //   }
  //   fsph_log->SetSensitiveDetector( FluxSD );
  // }
  
  //Start with vacuum snout:
  G4double inch = 2.54*cm;

  //In the following dimensions, "left" and "right" are as viewed from downstream!!!!
  G4double SnoutBeamPlate_Width = 9.425*inch;
  G4double SnoutEarmPlate_Width = 34.364*inch;
  G4double SnoutHarmPlate_Width = 27.831*inch;

  G4double Snout_Height = 37.75*inch;
  G4double Snout_Thick = 1.0*inch;

  //Window dimensions:
  G4double SnoutEarmWindow_Rbend_corners = 5.750*inch;
  G4double SnoutEarmWindow_Width = 2.0*9.625*inch;
  G4double SnoutEarmWindow_Height = 2.0*15.875*inch;
  //G4double SnoutEarmOffset = -0.744*inch;

  G4double SnoutHarmWindow_Rbend_corners = 5.000*inch;
  G4double SnoutHarmWindow_Width = 2.0*8.438*inch;
  G4double SnoutHarmWindow_Height = 2.0*9.562*inch;
  //G4double SnoutRightOffset = -1.4595*inch;

  //Aluminium Dimensions - just need to be bigger than window
  G4double EarmWindowThick = 0.032*inch;
  G4double HarmWindowThick = 0.020*inch;
    
  G4double SnoutEarmWindow_xcenter = SnoutEarmPlate_Width/2.0 - 16.438*inch;
  G4double SnoutHarmWindow_xcenter = -SnoutHarmPlate_Width/2.0 + 15.375*inch;
  G4double SnoutEarmWindowAngle = 27.5*deg;
  G4double SnoutHarmWindowAngle = 22.0*deg;

  G4double SnoutUpstreamHoleDiameter = 4.870*inch; //All the way through:
  G4double SnoutDownstreamHoleDiameter = 5.010*inch; //To a depth of .26 inches from the front:
  G4double SnoutDownstreamHoleDepth = 0.260*inch;

  G4double SnoutBeamPlate_xcenter = SnoutBeamPlate_Width/2.0 - 4.591*inch;

  // x coord. relative to center plate! To get coordinate in hall, we want center hole to be positioned at x = 0:
  G4double SnoutBeamPlate_xcoord = -SnoutBeamPlate_xcenter; // 4.591 - w/2 = -.1215 inch
  G4double SnoutBeamPlate_zcoord = 48.56*inch; //distance to upstream edge of snout beam plate
  //If we keep the target center as the origin of Hall A for detector positioning, then we will have to offset everything else in the
  //target chamber construction
  G4double TargetCenter_zoffset = 6.50*inch; //offset of target center wrt scattering chamber:
  
  //Make Boxes for Snout plates:
  G4Box *SnoutBeamPlate_Box = new G4Box("SnoutBeamPlate_Box", SnoutBeamPlate_Width/2.0, Snout_Height/2.0, Snout_Thick/2.0 );
  G4Box *SnoutEarmPlate_Box = new G4Box("SnoutEarmPlate_Box", SnoutEarmPlate_Width/2.0, Snout_Height/2.0, Snout_Thick/2.0 );
  G4Box *SnoutHarmPlate_Box = new G4Box("SnoutHarmPlate_Box", SnoutHarmPlate_Width/2.0, Snout_Height/2.0, Snout_Thick/2.0 );

  //Make cylinder cutouts for center plate:
  G4Tubs *SnoutBeamPlate_ThroughHole = new G4Tubs( "SnoutBeamPlate_ThroughHole", 0.0, SnoutUpstreamHoleDiameter/2.0, Snout_Thick/2.0 + mm, 0.0, twopi );
  G4Tubs *SnoutBeamPlate_CounterBore = new G4Tubs( "SnoutBeamPlate_CounterBore", 0.0, SnoutDownstreamHoleDiameter/2.0, SnoutDownstreamHoleDepth/2.0 + mm, 0.0, twopi );

  //Make subtraction solid(s) for beam plate:
  G4SubtractionSolid *SnoutBeamPlate_cut1 = new G4SubtractionSolid( "SnoutBeamPlate_cut1", SnoutBeamPlate_Box, SnoutBeamPlate_ThroughHole, 0, G4ThreeVector( SnoutBeamPlate_xcenter, 0, 0 ) );
  // z - (d/2+mm) = t/2 - d --> z = t/2 - d + d/2 + mm = t/2 - d/2 + mm
  G4double zoff_cbore = Snout_Thick/2.0 - SnoutDownstreamHoleDepth/2.0 + mm;
    
  G4SubtractionSolid *SnoutBeamPlate_cut2 = new G4SubtractionSolid( "SnoutBeamPlate_cut2", SnoutBeamPlate_cut1, SnoutBeamPlate_CounterBore, 0, G4ThreeVector( SnoutBeamPlate_xcenter, 0, zoff_cbore ) );

  //Make Window Box cutouts for Earm and Harm plates:
  G4Box *EarmWindowCutout_box = new G4Box("EarmWindowCutout_box", SnoutEarmWindow_Width/2.0, SnoutEarmWindow_Height/2.0, Snout_Thick/2.0+mm );
  G4SubtractionSolid *EarmPlate_cut = new G4SubtractionSolid( "EarmPlate_cut", SnoutEarmPlate_Box, EarmWindowCutout_box, 0, G4ThreeVector( SnoutEarmWindow_xcenter, 0, 0 ) );
							      
  G4Box *HarmWindowCutout_box = new G4Box("HarmWindowCutout_box", SnoutHarmWindow_Width/2.0, SnoutHarmWindow_Height/2.0, Snout_Thick/2.0+mm );
  G4SubtractionSolid *HarmPlate_cut = new G4SubtractionSolid( "HarmPlate_cut", SnoutHarmPlate_Box, HarmWindowCutout_box, 0, G4ThreeVector( SnoutHarmWindow_xcenter, 0, 0 ) );

  //Make rounded corner pieces for later placement (subtraction of box and cylinder):
  G4Box *EarmWindowCorner_box = new G4Box("EarmWindowCorner_box", SnoutEarmWindow_Rbend_corners/2.0, SnoutEarmWindow_Rbend_corners/2.0, Snout_Thick/2.0 );

  G4Tubs *EarmWindowCorner_cyl = new G4Tubs("EarmWindowCorner_cyl", 0.0, SnoutEarmWindow_Rbend_corners, Snout_Thick/2.0 + mm, 0.0, twopi );

  G4SubtractionSolid *EarmWindowCorner = new G4SubtractionSolid( "EarmWindowCorner", EarmWindowCorner_box, EarmWindowCorner_cyl, 0, G4ThreeVector( SnoutEarmWindow_Rbend_corners/2.0, SnoutEarmWindow_Rbend_corners/2.0, 0.0 ) );
  G4LogicalVolume *EarmWindowCorner_log = new G4LogicalVolume( EarmWindowCorner, GetMaterial("Stainless_Steel"), "EarmWindowCorner_log" );
  
  G4Box *HarmWindowCorner_box = new G4Box("HarmWindowCorner_box", SnoutHarmWindow_Rbend_corners/2.0, SnoutHarmWindow_Rbend_corners/2.0, Snout_Thick/2.0 );

  G4Tubs *HarmWindowCorner_cyl = new G4Tubs("HarmWindowCorner_cyl", 0.0, SnoutHarmWindow_Rbend_corners, Snout_Thick/2.0 + mm, 0.0, twopi );

  G4SubtractionSolid *HarmWindowCorner = new G4SubtractionSolid( "HarmWindowCorner", HarmWindowCorner_box, HarmWindowCorner_cyl, 0, G4ThreeVector( SnoutHarmWindow_Rbend_corners/2.0, SnoutHarmWindow_Rbend_corners/2.0, 0.0 ) );
  G4LogicalVolume *HarmWindowCorner_log = new G4LogicalVolume( HarmWindowCorner, GetMaterial("Stainless_Steel"), "HarmWindowCorner_log" );
  
  //Now we want to put together the three plates that make the Snout:
  G4ThreeVector Harm_zaxis( -sin(SnoutHarmWindowAngle), 0, cos(SnoutHarmWindowAngle ) );
  G4ThreeVector Harm_yaxis(0,1,0);
  G4ThreeVector Harm_xaxis( cos(SnoutHarmWindowAngle), 0, sin(SnoutHarmWindowAngle) );

  G4ThreeVector FrontRightCorner_pos_local( -SnoutBeamPlate_Width/2.0, 0, Snout_Thick/2.0 );

  G4ThreeVector HarmPlate_pos_relative = FrontRightCorner_pos_local - Snout_Thick/2.0 * Harm_zaxis - SnoutHarmPlate_Width/2.0 * Harm_xaxis;

  G4ThreeVector Earm_zaxis( sin(SnoutEarmWindowAngle), 0, cos(SnoutEarmWindowAngle) );
  G4ThreeVector Earm_yaxis(0,1,0);
  G4ThreeVector Earm_xaxis( cos(SnoutEarmWindowAngle), 0, -sin(SnoutEarmWindowAngle) );
  
  G4ThreeVector FrontLeftCorner_pos_local( SnoutBeamPlate_Width/2.0, 0, Snout_Thick/2.0 );
  G4ThreeVector EarmPlate_pos_relative = FrontLeftCorner_pos_local - Snout_Thick/2.0 * Earm_zaxis + SnoutEarmPlate_Width/2.0 * Earm_xaxis;

  G4RotationMatrix *rot_harm_window = new G4RotationMatrix;
  rot_harm_window->rotateY( SnoutHarmWindowAngle );

  G4UnionSolid *BeamHarmUnion = new G4UnionSolid( "BeamHarmUnion", SnoutBeamPlate_cut2, HarmPlate_cut, rot_harm_window, HarmPlate_pos_relative );
  G4RotationMatrix *rot_earm_window = new G4RotationMatrix;
  rot_earm_window->rotateY( -SnoutEarmWindowAngle );
  G4UnionSolid *Snout_solid = new G4UnionSolid("Snout_solid", BeamHarmUnion, EarmPlate_cut, rot_earm_window, EarmPlate_pos_relative );
  
  G4LogicalVolume *Snout_log = new G4LogicalVolume( Snout_solid, GetMaterial("Stainless_Steel"), "Snout_log" );

  G4ThreeVector Snout_position_global( SnoutBeamPlate_xcoord, 0, SnoutBeamPlate_zcoord - TargetCenter_zoffset + Snout_Thick/2.0 );

  new G4PVPlacement( 0, Snout_position_global, Snout_log, "Snout_phys", worldlog, false, 0 );

  //Fill the cutouts with vacuum, and then with steel corner pieces:
  G4LogicalVolume *SnoutHarmWindowCutout_log = new G4LogicalVolume( HarmWindowCutout_box, GetMaterial("Vacuum"), "SnoutHarmWindowCutout_log" );
  G4LogicalVolume *SnoutEarmWindowCutout_log = new G4LogicalVolume( EarmWindowCutout_box, GetMaterial("Vacuum"), "SnoutEarmWindowCutout_log" );

  G4RotationMatrix *rot_temp = new G4RotationMatrix;

  G4double xtemp, ytemp;
  xtemp = SnoutHarmWindow_Width/2.0 - SnoutHarmWindow_Rbend_corners/2.0;
  ytemp = SnoutHarmWindow_Height/2.0 - SnoutHarmWindow_Rbend_corners/2.0;
  
  new G4PVPlacement( 0, G4ThreeVector( -xtemp, -ytemp, 0 ), HarmWindowCorner_log, "HarmWindowCorner_phys_bottom_right", SnoutHarmWindowCutout_log, false, 0 );

  rot_temp->rotateZ( 90.0*deg ); //clockwise as viewed from downstream, ccw as viewed from upstream!
  
  new G4PVPlacement( rot_temp, G4ThreeVector( -xtemp, ytemp, 0 ), HarmWindowCorner_log, "HarmWindowCorner_phys_top_right", SnoutHarmWindowCutout_log, false, 1 );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateZ( -90.0*deg );

  new G4PVPlacement( rot_temp, G4ThreeVector( xtemp, -ytemp, 0 ), HarmWindowCorner_log, "HarmWindowCorner_phys_bottom_left", SnoutHarmWindowCutout_log, false, 2 );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateZ( 180.0*deg );

  new G4PVPlacement( rot_temp, G4ThreeVector( xtemp, ytemp, 0 ), HarmWindowCorner_log, "HarmWindowCorner_phys_top_left", SnoutHarmWindowCutout_log, false, 3 );

  //Finally, place the window cutout globally:

  G4ThreeVector HarmCutout_pos_global = Snout_position_global + FrontRightCorner_pos_local - Snout_Thick/2.0 * Harm_zaxis - ( SnoutHarmPlate_Width/2.0 - SnoutHarmWindow_xcenter ) * Harm_xaxis;

  new G4PVPlacement( rot_harm_window, HarmCutout_pos_global, SnoutHarmWindowCutout_log, "SnoutHarmWindowCutout_phys", worldlog, false, 0 );

  xtemp = SnoutEarmWindow_Width/2.0 - SnoutEarmWindow_Rbend_corners/2.0;
  ytemp = SnoutEarmWindow_Height/2.0 - SnoutEarmWindow_Rbend_corners/2.0;

  new G4PVPlacement( 0, G4ThreeVector(-xtemp,-ytemp,0), EarmWindowCorner_log, "EarmWindowCorner_phys_bottom_right", SnoutEarmWindowCutout_log, false, 0 );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateZ(90.0*deg );

  new G4PVPlacement( rot_temp, G4ThreeVector( -xtemp,ytemp,0), EarmWindowCorner_log, "EarmWindowCorner_phys_top_right", SnoutEarmWindowCutout_log, false, 1 );
  
  rot_temp = new G4RotationMatrix;
  rot_temp->rotateZ(-90.0*deg );

  new G4PVPlacement( rot_temp, G4ThreeVector( xtemp,-ytemp,0), EarmWindowCorner_log, "EarmWindowCorner_phys_bottom_left", SnoutEarmWindowCutout_log, false, 2 );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateZ(180.0*deg );

  new G4PVPlacement( rot_temp, G4ThreeVector( xtemp,ytemp,0), EarmWindowCorner_log, "EarmWindowCorner_phys_top_left", SnoutEarmWindowCutout_log, false, 3 );
  
  G4ThreeVector EarmCutout_pos_global = Snout_position_global + FrontLeftCorner_pos_local - Snout_Thick/2.0 * Earm_zaxis + (SnoutEarmPlate_Width/2.0 + SnoutEarmWindow_xcenter) * Earm_xaxis;

  new G4PVPlacement( rot_earm_window, EarmCutout_pos_global, SnoutEarmWindowCutout_log, "SnoutEarmWindowCutout_phys", worldlog, false, 0 );

  //What's next? Define scattering chamber vacuum volume:
  G4double ScatChamberRadius = 23.80*inch;
  G4double ScatChamberHeight = Snout_Height;

  G4Tubs *ScatChamber_solid = new G4Tubs("ScatChamber_solid", 0, ScatChamberRadius, ScatChamberHeight/2.0, 0.0, twopi );
  G4LogicalVolume *ScatChamber_log = new G4LogicalVolume( ScatChamber_solid, GetMaterial("Vacuum"), "ScatChamber_log" );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateX( 90.0*deg );

  new G4PVPlacement( rot_temp, G4ThreeVector(0,0,-TargetCenter_zoffset), ScatChamber_log, "ScatChamber_phys", worldlog, false, 0 );
  
  rot_temp = new G4RotationMatrix;
  rot_temp->rotateX( -90.0*deg );
  
  //Call BuildStandardCryoTarget HERE !
  BuildStandardCryoTarget(ScatChamber_log, rot_temp, G4ThreeVector(0, -TargetCenter_zoffset, 0));
  
  /*
  //HERE
  //Now let's make a cryotarget:
  G4double Rcell = 4.0*cm;
  G4double uthick = 0.1*mm;
  G4double dthick = 0.15*mm;
  G4double sthick = 0.2*mm;

  G4Tubs *TargetMother_solid = new G4Tubs("TargetMother_solid", 0, Rcell + sthick, (fTargLen+uthick+dthick)/2.0, 0.0, twopi );
  G4LogicalVolume *TargetMother_log = new G4LogicalVolume(TargetMother_solid, GetMaterial("Vacuum"), "TargetMother_log" );
  
  G4Tubs *TargetCell = new G4Tubs("TargetCell", 0, Rcell, fTargLen/2.0, 0, twopi );

  G4LogicalVolume *TargetCell_log;

  if( fTargType == kLH2 ){
    TargetCell_log = new G4LogicalVolume( TargetCell, GetMaterial("LH2"), "TargetCell_log" );
  } else {
    TargetCell_log = new G4LogicalVolume( TargetCell, GetMaterial("LD2"), "TargetCell_log" );
  }
  
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
  G4double ztemp = -(fTargLen+uthick+dthick)/2.0;
  //place upstream window:
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+uthick/2.0), uwindow_log, "uwindow_phys", TargetMother_log, false, 0 );
  //place target and side walls:
  ztemp += uthick;
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+fTargLen/2.0), TargetCell_log, "TargetCell_phys", TargetMother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+fTargLen/2.0), TargetWall_log, "TargetWall_phys", TargetMother_log, false, 0 );
  ztemp += fTargLen;
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+dthick/2.0), dwindow_log, "dwindow_phys", TargetMother_log, false, 0 );

  G4double targ_zcenter = (uthick-dthick)/2.0; //position of target center relative to target mother volume
  
  //Compute position of target relative to scattering chamber:
  //The target center should be at
  rot_temp = new G4RotationMatrix;
  rot_temp->rotateX( -90.0*deg );

  //for z of target center to be at zero, 
  
  new G4PVPlacement( rot_temp, G4ThreeVector(0, -(TargetCenter_zoffset-targ_zcenter) ,0), TargetMother_log, "TargetMother_phys", ScatChamber_log, false, 0 );
  */
  
  //Next: Snout vacuum geometry (complicated!)

  //Angles of snout side walls relative to (scattering chamber) origin:
  G4double phisnout_harm = (90.0 - 54.88)*deg; //35.12 degrees
  G4double phisnout_earm = (90.0 - 46.82)*deg; //43.18 degrees

  //The global coordinates of the center of the Earm snout plate:
  G4ThreeVector EarmPlanePos_global = Snout_position_global + EarmPlate_pos_relative - Snout_Thick/2.0 * Earm_zaxis;
  G4ThreeVector HarmPlanePos_global = Snout_position_global + HarmPlate_pos_relative - Snout_Thick/2.0 * Harm_zaxis;

  G4ThreeVector ScatChamberPos_global = G4ThreeVector( 0, 0, -TargetCenter_zoffset );
  
  //The Snout vacuum volume will be the INTERSECTION of a polycone and an extruded solid:
  std::vector<G4TwoVector> snout_polygon;

  xtemp = ScatChamberRadius * sin( phisnout_earm );
  ytemp = ScatChamberRadius * cos( phisnout_earm ) - TargetCenter_zoffset;

  snout_polygon.push_back( G4TwoVector( xtemp, ytemp ) );

  //G4double distance_schbr_origin_to_earm_plane = (EarmPlanePos_global - ScatChamberPos_global).mag(); 

  //Intersection with a line and a plane:
  G4ThreeVector nhat_snout_earm_side( sin( phisnout_earm ), 0, cos( phisnout_earm ) );
  G4ThreeVector nhat_snout_harm_side( -sin( phisnout_harm ), 0, cos( phisnout_harm ) );

  // (r0 + nhat * s - rplane ) dot nplane = 0
  // s = (rplane - r0) dot nplane / (nhat dot nplane )

  G4double s_temp = (EarmPlanePos_global - ScatChamberPos_global).dot( Earm_zaxis )/(nhat_snout_earm_side.dot( Earm_zaxis ) );
  
  G4ThreeVector Snout_earm_side_intercept = ScatChamberPos_global + s_temp * nhat_snout_earm_side;

  xtemp = Snout_earm_side_intercept.x();
  ytemp = Snout_earm_side_intercept.z();

  snout_polygon.push_back( G4TwoVector( xtemp, ytemp ) );

  //inner corner position of snout on earm side:
  G4ThreeVector Snout_earm_corner_pos_rel = Snout_position_global + FrontLeftCorner_pos_local - Snout_Thick*Earm_zaxis + Snout_Thick*tan( SnoutEarmWindowAngle/2.0 ) * Earm_xaxis;
 
  xtemp = Snout_earm_corner_pos_rel.x();
  ytemp = Snout_earm_corner_pos_rel.z();

  snout_polygon.push_back( G4TwoVector( xtemp, ytemp ) );
  
  G4ThreeVector Snout_harm_corner_pos_rel = Snout_position_global + FrontRightCorner_pos_local - Snout_Thick*Harm_zaxis - Snout_Thick*tan( SnoutHarmWindowAngle/2.0 ) * Harm_xaxis;

  xtemp = Snout_harm_corner_pos_rel.x();
  ytemp = Snout_harm_corner_pos_rel.z();

  snout_polygon.push_back( G4TwoVector( xtemp, ytemp ) );
  
  s_temp = (HarmPlanePos_global - ScatChamberPos_global).dot( Harm_zaxis ) / (nhat_snout_harm_side.dot( Harm_zaxis ) );

  G4ThreeVector Snout_harm_side_intercept = ScatChamberPos_global + s_temp * nhat_snout_harm_side;

  xtemp = Snout_harm_side_intercept.x();
  ytemp = Snout_harm_side_intercept.z();

  snout_polygon.push_back( G4TwoVector( xtemp, ytemp ) );

  //last point:
  xtemp = -ScatChamberRadius * sin( phisnout_harm );
  ytemp =  ScatChamberRadius * cos( phisnout_harm ) - TargetCenter_zoffset;

  snout_polygon.push_back( G4TwoVector( xtemp, ytemp ) );

  G4ExtrudedSolid *Snout_pgon = new G4ExtrudedSolid( "Snout_pgon", snout_polygon, Snout_Height/2.0,
						     G4TwoVector(0,0), 1.0,
						     G4TwoVector(0,0), 1.0 );

  //This extruded solid has the x axis to beam left, the yaxis along the beam direction, and the z axis in the -y direction. This means we want to rotate it by 90 degrees ccw about x when looking in the +x direction
  rot_temp = new G4RotationMatrix;
  rot_temp->rotateX( -90.0*deg );

  // G4LogicalVolume *Snout_pgon_log = new G4LogicalVolume( Snout_pgon, GetMaterial("Vacuum"), "Snout_pgon_log" );

  // //Positioning: The origin of this pgon is the same as the origin of the hall, so should be 0,0,0:
  // new G4PVPlacement( rot_temp, G4ThreeVector(0,0,0), Snout_pgon_log, "Snout_pgon_phys", worldlog, false, 0 );

  //Next: make the polycone:
  
  G4double SnoutTaperAngle = 38.81*deg;
  G4double SnoutTaper_zpos = 14.88*inch;

  const G4int nzplanes = 4;
  G4double zplanes[nzplanes] = {-Snout_Height/2.0,
				-(ScatChamberRadius-SnoutTaper_zpos)*tan(SnoutTaperAngle),
				(ScatChamberRadius-SnoutTaper_zpos)*tan(SnoutTaperAngle),
				Snout_Height/2.0 };
  G4double rplanes[nzplanes] = { ScatChamberRadius + (zplanes[1]-zplanes[0])/tan(SnoutTaperAngle),
				 ScatChamberRadius,
				 ScatChamberRadius,
				 ScatChamberRadius + (zplanes[3]-zplanes[2])/tan(SnoutTaperAngle) };

  G4double routplanes[nzplanes] = { ScatChamberRadius + (zplanes[1]-zplanes[0])/tan(SnoutTaperAngle),
				    ScatChamberRadius + (zplanes[1]-zplanes[0])/tan(SnoutTaperAngle),
				    ScatChamberRadius + (zplanes[3]-zplanes[2])/tan(SnoutTaperAngle),
				    ScatChamberRadius + (zplanes[3]-zplanes[2])/tan(SnoutTaperAngle) };

  for( G4int pl=0; pl<nzplanes; ++pl ){
    routplanes[pl] += 40.0*cm;
  }
  
  //We want the z axis of the polycone to coincide with the z axis of the polygon, so to have +x to beam left and +y along beam direction, we need to have +z be vertically down. In this case, the phi start would be
  G4double phistart = 46.82*deg;
  G4double dphi = phisnout_earm + phisnout_harm;

  G4Polycone *SnoutTaper = new G4Polycone( "SnoutTaper", phistart, dphi, nzplanes, zplanes, rplanes, routplanes );
  
  //The polycone has its central axis along the vertical line through the "origin", whereas the polygon has its central axis
  // along the vertical line through the target center;
  //Therefore, we have to offset the position of the polygon in y by +targ z center:
  G4IntersectionSolid *SnoutTaperPgon_intersect = new G4IntersectionSolid( "SnoutTaperPgon_intersect", SnoutTaper, Snout_pgon, 0, G4ThreeVector(0, TargetCenter_zoffset, 0) );

  G4LogicalVolume *SnoutVacuum_log = new G4LogicalVolume( SnoutTaperPgon_intersect, GetMaterial("Vacuum"), "SnoutVacuum_log" );

  new G4PVPlacement( rot_temp, G4ThreeVector(0,0,-TargetCenter_zoffset), SnoutVacuum_log, "SnoutVacuum_phys", worldlog, false, 0 );

  //Iron Tube inside the Snout vacuum volume:
  G4double IronTube_Rmin = 5.0*cm;
  G4double IronTube_Rmax = 7.0*cm;
  G4double IronTube_Thick = 6.5*cm;

  G4Tubs *IronTube = new G4Tubs("IronTube", IronTube_Rmin, IronTube_Rmax, IronTube_Thick/2.0, 0.0, twopi );
  G4LogicalVolume *IronTube_log = new G4LogicalVolume( IronTube, GetMaterial("Iron"), "IronTube_log" );

  IronTube_log->SetVisAttributes( new G4VisAttributes(G4Colour(0.3,0.3,0.3) ) );
  
  //relative to the origin of the snout polycone, the iron tube is at y = snout z - half thick:
  G4ThreeVector IronTube_pos_rel( 0, SnoutBeamPlate_zcoord - IronTube_Thick/2.0, 0 );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateX( 90.0*deg );

  new G4PVPlacement( rot_temp, IronTube_pos_rel, IronTube_log, "IronTube_phys", SnoutVacuum_log, false, 0 );
  

  //Finally, put windows and bolt plates:
  
  G4double HarmFlangeWidth = 2.0*11.44*inch;
  G4double HarmFlangeHeight = 2.0*12.56*inch;
  G4double HarmFlangeThick = 0.980*inch;

  G4Box *HarmAlWindow_box = new G4Box("HarmAlWindow_box", HarmFlangeWidth/2.0, HarmFlangeHeight/2.0, HarmWindowThick/2.0 );
  G4LogicalVolume *HarmWindow_log = new G4LogicalVolume( HarmAlWindow_box, GetMaterial("Aluminum"), "HarmWindow_log" );
  //Figure out the placement of the Harm window:
  G4ThreeVector Harm_window_pos = Snout_position_global + FrontRightCorner_pos_local + (-SnoutHarmPlate_Width/2.0 + SnoutHarmWindow_xcenter) * Harm_xaxis + (HarmWindowThick/2.0) * Harm_zaxis;

  new G4PVPlacement( rot_harm_window, Harm_window_pos, HarmWindow_log, "HarmWindow_phys", worldlog, false, 0 );
  
  G4Box *HarmFlange_box = new G4Box("HarmFlange_box", HarmFlangeWidth/2.0, HarmFlangeHeight/2.0, HarmFlangeThick/2.0 );

  G4SubtractionSolid *HarmFlange_cut = new G4SubtractionSolid( "HarmFlange_cut", HarmFlange_box, HarmWindowCutout_box, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *HarmFlange_log = new G4LogicalVolume( HarmFlange_cut, GetMaterial("Aluminum"), "HarmFlange_log" );

  G4ThreeVector HarmFlange_pos = Harm_window_pos + (HarmWindowThick + HarmFlangeThick)/2.0 * Harm_zaxis;

  new G4PVPlacement( rot_harm_window, HarmFlange_pos, HarmFlange_log, "HarmFlange_phys", worldlog, false, 0 );

  //Create a new logical volume of aluminum instead of steel for the rounded corners of the flange:
  G4LogicalVolume *HarmFlangeCorner_log = new G4LogicalVolume( HarmWindowCorner, GetMaterial("Aluminum"), "HarmFlangeCorner_log" );
  
  G4ThreeVector pos_temp;
  
  xtemp = SnoutHarmWindow_Width/2.0 - SnoutHarmWindow_Rbend_corners/2.0;
  ytemp = SnoutHarmWindow_Height/2.0 - SnoutHarmWindow_Rbend_corners/2.0;

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY( SnoutHarmWindowAngle );
  
  pos_temp = HarmFlange_pos - xtemp * Harm_xaxis - ytemp * Harm_yaxis;
  
  new G4PVPlacement( rot_temp, pos_temp, HarmFlangeCorner_log, "HarmFlangeCorner_phys_bottom_right", worldlog, false, 0 );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY( SnoutHarmWindowAngle );
  rot_temp->rotateZ( 90.0*deg );
  pos_temp = HarmFlange_pos - xtemp * Harm_xaxis + ytemp * Harm_yaxis;

  new G4PVPlacement( rot_temp, pos_temp, HarmFlangeCorner_log, "HarmFlangeCorner_phys_top_right", worldlog, false, 1 );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY( SnoutHarmWindowAngle );
  rot_temp->rotateZ( -90.0*deg );
  pos_temp = HarmFlange_pos + xtemp * Harm_xaxis - ytemp * Harm_yaxis;

  new G4PVPlacement( rot_temp, pos_temp, HarmFlangeCorner_log, "HarmFlangeCorner_phys_bottom_left", worldlog, false, 2 );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY( SnoutHarmWindowAngle );
  rot_temp->rotateZ( 180.0*deg );
  pos_temp = HarmFlange_pos + xtemp * Harm_xaxis + ytemp * Harm_yaxis;

  new G4PVPlacement( rot_temp, pos_temp, HarmFlangeCorner_log, "HarmFlangeCorner_phys_top_left", worldlog, false, 3 );

  G4double EarmFlangeWidth = 2.0*12.62*inch;
  G4double EarmFlangeHeight = Snout_Height;
  G4double EarmFlangeThick = 0.980*inch;

  G4Box *EarmAlWindow_box = new G4Box( "EarmAlWindow_box", EarmFlangeWidth/2.0, EarmFlangeHeight/2.0, EarmWindowThick/2.0 );

  G4LogicalVolume *EarmWindow_log = new G4LogicalVolume( EarmAlWindow_box, GetMaterial("Aluminum"), "EarmWindow_log" );
  G4ThreeVector Earm_window_pos = Snout_position_global + FrontLeftCorner_pos_local + (SnoutEarmPlate_Width/2.0 + SnoutEarmWindow_xcenter) * Earm_xaxis + EarmWindowThick/2.0 * Earm_zaxis;

  new G4PVPlacement( rot_earm_window, Earm_window_pos, EarmWindow_log, "EarmWindow_phys", worldlog, false, 0 );
  
  //Next: make flange:

  G4Box *EarmFlangeBox = new G4Box( "EarmFlangeBox", EarmFlangeWidth/2.0, EarmFlangeHeight/2.0, EarmFlangeThick/2.0 );
  G4SubtractionSolid *EarmFlange_cut = new G4SubtractionSolid( "EarmFlange_cut", EarmFlangeBox, EarmWindowCutout_box, 0, G4ThreeVector(0,0,0));

  G4LogicalVolume *EarmFlange_log = new G4LogicalVolume( EarmFlange_cut, GetMaterial("Aluminum"), "EarmFlange_log" );

  G4ThreeVector EarmFlange_pos = Earm_window_pos + (EarmWindowThick + EarmFlangeThick)/2.0 * Earm_zaxis;

  new G4PVPlacement( rot_earm_window, EarmFlange_pos, EarmFlange_log, "EarmFlange_phys", worldlog, false, 0 );

  G4LogicalVolume *EarmFlangeCorner_log = new G4LogicalVolume( EarmWindowCorner, GetMaterial("Aluminum"), "EarmFlangeCorner_log" );

  xtemp = SnoutEarmWindow_Width/2.0 - SnoutEarmWindow_Rbend_corners/2.0;
  ytemp = SnoutEarmWindow_Height/2.0 - SnoutEarmWindow_Rbend_corners/2.0;

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY( -SnoutEarmWindowAngle );
  
  pos_temp = EarmFlange_pos - xtemp * Earm_xaxis - ytemp * Earm_yaxis;
  
  new G4PVPlacement( rot_temp, pos_temp, EarmFlangeCorner_log, "EarmFlangeCorner_phys_bottom_right", worldlog, false, 0 );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY( -SnoutEarmWindowAngle );
  rot_temp->rotateZ( 90.0*deg );
  pos_temp = EarmFlange_pos - xtemp * Earm_xaxis + ytemp * Earm_yaxis;
  
  new G4PVPlacement( rot_temp, pos_temp, EarmFlangeCorner_log, "EarmFlangeCorner_phys_top_right", worldlog, false, 1 );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY( -SnoutEarmWindowAngle );
  rot_temp->rotateZ( -90.0*deg );
  pos_temp = EarmFlange_pos + xtemp * Earm_xaxis - ytemp * Earm_yaxis;
  
  new G4PVPlacement( rot_temp, pos_temp, EarmFlangeCorner_log, "EarmFlangeCorner_phys_bottom_left", worldlog, false, 2 );

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateY( -SnoutEarmWindowAngle );
  rot_temp->rotateZ( 180.0*deg );
  pos_temp = EarmFlange_pos + xtemp * Earm_xaxis + ytemp * Earm_yaxis;
  
  new G4PVPlacement( rot_temp, pos_temp, EarmFlangeCorner_log, "EarmFlangeCorner_phys_top_left", worldlog, false, 3 );
  
  // VISUALS
  // Mother Volumes

  // Vacuum

  G4VisAttributes *Snout_VisAtt = new G4VisAttributes( G4Colour( 0.6, 0.55, 0.65 ) );
  Snout_log->SetVisAttributes( Snout_VisAtt );
  HarmWindowCorner_log->SetVisAttributes( Snout_VisAtt );
  EarmWindowCorner_log->SetVisAttributes( Snout_VisAtt );

  G4VisAttributes *Flange_VisAtt = new G4VisAttributes( G4Colour( 0.15, 0.8, 0.9 ) );
  HarmFlange_log->SetVisAttributes( Flange_VisAtt );
  HarmFlangeCorner_log->SetVisAttributes( Flange_VisAtt );

  EarmFlange_log->SetVisAttributes( Flange_VisAtt );
  EarmFlangeCorner_log->SetVisAttributes( Flange_VisAtt );
  // PlateUnionLog->SetVisAttributes( Snout_VisAtt );
  // LeftCornerLog->SetVisAttributes( Snout_VisAtt );
  // RightCornerLog->SetVisAttributes( Snout_VisAtt );
  // RightWindowCutoutVacuum_log->SetVisAttributes( G4VisAttributes::Invisible );
  // LeftWindowCutoutVacuum_log->SetVisAttributes( G4VisAttributes::Invisible );

  SnoutHarmWindowCutout_log->SetVisAttributes( G4VisAttributes::Invisible );
  SnoutEarmWindowCutout_log->SetVisAttributes( G4VisAttributes::Invisible );
  
  G4VisAttributes *Snout_VisAttWire = new G4VisAttributes( G4Colour( 0.6, 0.55, 0.65 ) );
  Snout_VisAttWire->SetForceWireframe(true);

  G4VisAttributes *Window_visatt = new G4VisAttributes( G4Colour( 0.9, 0.8, 0.15 ) );
  Window_visatt->SetForceWireframe(true);
  HarmWindow_log->SetVisAttributes( Window_visatt );
  EarmWindow_log->SetVisAttributes( Window_visatt );
  //Snout_ScChamb->SetVisAttributes( Snout_VisAttWire );

  ScatChamber_log->SetVisAttributes( Snout_VisAttWire );
  SnoutVacuum_log->SetVisAttributes( Snout_VisAttWire );
  
  // G4VisAttributes *Targ_visatt = new G4VisAttributes( G4Colour( 0.1, 0.05, 0.9 ) );
  // TargetCell_log->SetVisAttributes( Targ_visatt );
  // G4VisAttributes *TargWall_visatt = new G4VisAttributes( G4Colour( 0.9, .05, 0.1 ) );
  // TargWall_visatt->SetForceWireframe( true );
  // TargetWall_log->SetVisAttributes( TargWall_visatt );
  // uwindow_log->SetVisAttributes( TargWall_visatt );
  // dwindow_log->SetVisAttributes( TargWall_visatt );
  // TargetMother_log->SetVisAttributes( G4VisAttributes::Invisible );
  
  G4VisAttributes *AlColor= new G4VisAttributes(G4Colour(0.4,0.4,0.4));
  //LeftAl_Log->SetVisAttributes( AlColor );
  //RightAl_Log->SetVisAttributes( AlColor );

  G4VisAttributes *ironColor= new G4VisAttributes(G4Colour(0.52,0.47,0.47));
  // TUB1_log->SetVisAttributes( ironColor );
  // TUB2_log->SetVisAttributes( ironColor );
  // TUB3_log->SetVisAttributes( ironColor );
}

// This function has replaced BuildC16CryoTarget
// It is identical except it uses the function BuildStandardCryoTarget.
void G4SBSTargetBuilder::BuildC16ScatCham(G4LogicalVolume *worldlog ){
  G4RotationMatrix *rot_temp;
  
  G4double inch = 2.54*cm;
  
  G4double entpipe_rin = 31.75*mm;
  G4double entpipe_rout = entpipe_rin + 0.12*mm;

  double sheight = 1.2*m;
  double swinthick    = 0.38*mm;
  double swallrad     = 1.143*m/2.0;
  double swallrad_in  = 1.041*m/2.0;
  
  double hcal_ang_min = -55.0*deg;
  double hcal_ang_max =  45.0*deg;
  double hcal_win_h = 15.2*cm;

  double bb_ang_min = 18.0*deg;
  double bb_ang_max = 80.0*deg;
  double bb_win_h = 0.5*m;

  //G4double extpipe_rin = 48.00*mm;
  G4double extpipe_rin = 6.0*inch/2.0;
  G4double extpipestart = 1.622*m;
  
  //  if( bb_ang_min < hcal_ang_max ) bb_ang_min = hcal_ang_max + (swallrad-swallrad_in)/swallrad/4;
  G4double dvcs_ang_min = 8.0*deg;
  G4double dvcs_ang_max = (90.0+52.0)*deg; //142 deg

  G4double angthick_dvcs_snout = 2.5*deg; 
  G4double Rin_dvcs_snout = 29.315*inch;
  G4double Rout_dvcs_snout = 30.315*inch;

  G4double extpipe_len = extpipestart - Rin_dvcs_snout;
  
  G4double height_dvcs_snout = 12.5*inch;
  G4double Rin_dvcs_beampipe = 6.065*inch/2.0;
  G4double Rout_dvcs_beampipe = Rin_dvcs_beampipe + 0.28*inch;
  
  // Make the chamber
  G4Tubs *swall = new G4Tubs("scham_wall", swallrad_in, swallrad, sheight/2.0, 0.0*deg, 360.0*deg );
  G4Tubs *dvcs_snout_cut = new G4Tubs("dvcs_snout_cut", swallrad_in - 1.0*cm, swallrad + 1.0*cm, height_dvcs_snout/2.0, dvcs_ang_min, dvcs_ang_max-dvcs_ang_min );
  
  // Cut out for windows
  // G4Tubs *swall_hcalcut = new G4Tubs("scham_wall_hcalcut", swallrad_in-2*cm, swallrad+2*cm, hcal_win_h, hcal_ang_min, hcal_ang_max-hcal_ang_min);
  // G4Tubs *swall_bbcut = new G4Tubs("scham_wall_bbcut", swallrad_in-2*cm, swallrad+2*cm, bb_win_h, bb_ang_min, bb_ang_max-bb_ang_min);

  //  G4Tubs *bigdvcscut = new G4Tubs("bigdvcscut", 

  rot_temp = new G4RotationMatrix;
  rot_temp->rotateX(-90.0*deg);
  rot_temp->rotateY(180.0*deg);
  
  G4SubtractionSolid *swallcut1 = new G4SubtractionSolid( "swallcut1", swall, dvcs_snout_cut );

  G4Tubs *upstream_beam_hole = new G4Tubs("upstream_beam_hole", 0.0, entpipe_rin, 5.0*inch, 0, 360.*deg);
  G4SubtractionSolid *swallcut2 = new G4SubtractionSolid( "swallcut2", swallcut1, upstream_beam_hole, rot_temp, G4ThreeVector(0,-0.5*(swallrad+swallrad_in),0) );
  
  G4Tubs *dvcs_snout = new G4Tubs("dvcs_snout", Rin_dvcs_snout, Rout_dvcs_snout, height_dvcs_snout/2.0, dvcs_ang_min - angthick_dvcs_snout, dvcs_ang_max - dvcs_ang_min + 2.*angthick_dvcs_snout );

  G4Tubs *dvcs_beamleft_cut = new G4Tubs("dvcs_beamleft_cut", 0.0, Rout_dvcs_snout + cm, 5.0*inch/2.0, 100.6*deg, 38.8*deg );

  G4Tubs *dvcs_beamright_cut = new G4Tubs("dvcs_beamright_cut", 0.0, Rout_dvcs_snout + cm, 5.0*inch/2.0, 10.6*deg, 47.6*deg );

  G4SubtractionSolid *dvcs_snout_cut1 = new G4SubtractionSolid( "dvcs_snout_cut1", dvcs_snout, dvcs_beamleft_cut );
  G4SubtractionSolid *dvcs_snout_cut2 = new G4SubtractionSolid( "dvcs_snout_cut2", dvcs_snout_cut1, dvcs_beamright_cut );

  G4Tubs *dvcs_snout_vacuum = new G4Tubs("dvcs_snout_vacuum", swallrad, Rout_dvcs_snout, height_dvcs_snout/2.0, dvcs_ang_min - angthick_dvcs_snout, dvcs_ang_max - dvcs_ang_min + 2.*angthick_dvcs_snout );

  G4LogicalVolume *dvcs_snout_vacuum_log = new G4LogicalVolume(dvcs_snout_vacuum, GetMaterial("Vacuum"), "dvcs_snout_vacuum_log" );

  G4Tubs *dvcs_snout_beamhole = new G4Tubs("dvcs_snout_beamhole", 0.0, Rin_dvcs_beampipe, 5.0*inch, 0.0, 360.0*deg );

  
  
  G4SubtractionSolid *dvcs_snout_cut3 = new G4SubtractionSolid( "dvcs_snout_cut3", dvcs_snout_cut2, dvcs_snout_beamhole, rot_temp, G4ThreeVector( 0.0, 0.5*(Rin_dvcs_snout+Rout_dvcs_snout), 0.0) );

  G4LogicalVolume *dvcs_snout_log = new G4LogicalVolume(dvcs_snout_cut3, GetMaterial("Aluminum"), "dvcs_snout_log");

  new G4PVPlacement( 0, G4ThreeVector(0,0,0), dvcs_snout_log, "dvcs_snout_phys", dvcs_snout_vacuum_log, false, 0 );

  new G4PVPlacement( rot_temp, G4ThreeVector(0,0,0), dvcs_snout_vacuum_log, "dvcs_snout_vacuum_phys", worldlog, false, 0);

  G4double dvcs_win_thick = (30.331-30.315)*inch;

  G4Tubs *dvcs_beamleft_window = new G4Tubs("dvcs_beamleft_window", Rout_dvcs_snout, Rout_dvcs_snout + dvcs_win_thick, 7.328*inch/2.0, (100.6-2.23)*deg, (38.8+4.46)*deg );
  G4LogicalVolume *dvcs_beamleft_window_log = new G4LogicalVolume( dvcs_beamleft_window, GetMaterial("Aluminum"), "dvcs_beamleft_window_log" );

  new G4PVPlacement( rot_temp, G4ThreeVector(0,0,0), dvcs_beamleft_window_log, "dvcs_beamleft_window_phys", worldlog, false, 0 );

  G4Tubs *dvcs_beamright_window = new G4Tubs("dvcs_beamright_window", Rout_dvcs_snout, Rout_dvcs_snout + dvcs_win_thick, 7.328*inch/2.0, 8.0*deg, (47.6+4.4)*deg );
  G4LogicalVolume *dvcs_beamright_window_log = new G4LogicalVolume( dvcs_beamright_window, GetMaterial("Aluminum"), "dvcs_beamright_window_log" );

  new G4PVPlacement( rot_temp, G4ThreeVector(0,0,0), dvcs_beamright_window_log, "dvcs_beamright_window_phys", worldlog, false, 0 );
  
  G4VisAttributes *window_visatt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0.0 ) );
  window_visatt->SetForceWireframe(true);
  dvcs_beamleft_window_log->SetVisAttributes( window_visatt );
  dvcs_beamright_window_log->SetVisAttributes( window_visatt );
  
  G4Tubs *scham_vacuum = new G4Tubs("scham_vacuum", 0.0, swallrad, sheight/2.0, 0.0, 360.0*deg );
  G4LogicalVolume *scham_vacuum_log = new G4LogicalVolume( scham_vacuum, GetMaterial("Vacuum"), "scham_vacuum_log" );
  G4LogicalVolume *scham_wall_log = new G4LogicalVolume( swallcut2, GetMaterial("Aluminum"), "scham_wall_log" );

  new G4PVPlacement( 0, G4ThreeVector(0,0,0), scham_wall_log, "scham_wall_phys", scham_vacuum_log, false, 0 );
  new G4PVPlacement( rot_temp, G4ThreeVector(0,0,0), scham_vacuum_log, "scham_vacuum_phys", worldlog, false, 0 );
  
  //swallcut = new G4SubtractionSolid("swallcut2", swallcut, swall_bbcut);

  // G4Tubs *swall_hcalwin = new G4Tubs("scham_wall_hcalwin", swallrad_in, swallrad_in+swinthick, hcal_win_h, hcal_ang_min, hcal_ang_max-hcal_ang_min);
  // G4Tubs *swall_bbwin = new G4Tubs("scham_wall_bbwin", swallrad_in, swallrad_in+swinthick, bb_win_h, bb_ang_min, bb_ang_max-bb_ang_min)
  // Exit pipe, prepping to bore out holes in scattering chamber
  
  // G4Tubs *exttube = new G4Tubs("exitpipetube", extpipe_rin, extpipe_rin+0.28*inch, extpipe_len/2, 0.*deg, 360.*deg );
  // G4Tubs *extvactube = new G4Tubs("exitpipetube_vac", 0.0, extpipe_rin, extpipe_len, 0.*deg, 360.*deg );

  // G4LogicalVolume *extpipe_log = new G4LogicalVolume(exttube, GetMaterial("Aluminum"),"extpipe_log");
  // G4LogicalVolume *extvac_log = new G4LogicalVolume(extvactube, GetMaterial("Vacuum"),"extvac_log");

  // double tempZ = 5.0*cm + 1*cm;
  // G4Tubs *swall_enthole = new G4Tubs("scham_wall_enthole", 0.0, entpipe_rout, tempZ, 0.0*deg, 360.0*deg );
  // G4Tubs *swall_exthole = new G4Tubs("scham_wall_exthole", 0.0, extpipe_rin, tempZ, 0.*deg, 360.*deg );

  // G4RotationMatrix *chamholerot = new G4RotationMatrix;
  // chamholerot->rotateY(90.0*deg);

  //  Cut holes in the scattering chamber and HCal aluminum window 
  // G4SubtractionSolid* swall_holes = new G4SubtractionSolid("swall_enthole", swallcut, swall_enthole, chamholerot, G4ThreeVector(-(swallrad+swallrad_in)/2.0, 0.0, 0.0) );
  // swall_holes = new G4SubtractionSolid("swall_holes", swall_holes, swall_exthole, chamholerot, G4ThreeVector((swallrad+swallrad_in)/2, 0.0, 0.0) );

  // G4SubtractionSolid *swall_hcalwin_cut = new G4SubtractionSolid("swall_hcalwin_cut", swall_hcalwin, swall_exthole, chamholerot, G4ThreeVector((swallrad+swallrad_in-5*cm)/2.0,0,0) );

  // // Fill the Entry / Exit holes with vacuum, and place them 
  // tempZ = 5.1*cm;
  // swall_enthole = new G4Tubs("scham_wall_enthole", 0.0, entpipe_rout, tempZ/2.0, 0.0*deg, 360.0*deg );
  // G4LogicalVolume *sc_entry_hole_vacuum_log = new G4LogicalVolume( swall_enthole, GetMaterial("Vacuum"), "sc_entry_hole_vacuum_log");
  // new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -(swallrad+swallrad_in)/2.0), sc_entry_hole_vacuum_log, "sc_entry_hole_vacuum_phys",
  // 		    worldlog, false, 0 );

  // swall_exthole = new G4Tubs("scham_wall_exthole", 0.0, extpipe_rin, tempZ/2.0, 0.0*deg, 360.0*deg );
  // G4LogicalVolume *sc_exit_hole_vacuum_log = new G4LogicalVolume( swall_exthole, GetMaterial("Vacuum"), "sc_exit_hole_vacuum_log" );
  // new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, (swallrad+swallrad_in)/2.0), sc_exit_hole_vacuum_log, "sc_exit_hole_vacuum_phys",
  // 		    worldlog, false, 0 );

  // // Turn all subtraction solids into Logical Volumes of the appropriate material.
  // G4LogicalVolume *swall_log = new G4LogicalVolume(swall_holes, GetMaterial("Aluminum"),"scham_wall_log");
  // G4LogicalVolume *sc_hcalwin_log = new G4LogicalVolume(swall_hcalwin_cut, GetMaterial("Aluminum"),"sc_hcalwin_log");
  // G4LogicalVolume *sc_bbwin_log = new G4LogicalVolume(swall_bbwin, GetMaterial("Aluminum"),"sc_bbwin_log");

  // G4RotationMatrix *schamrot = new G4RotationMatrix;
  // schamrot->rotateX(-90.0*deg);
  // schamrot->rotateZ(-90.0*deg);

  // // Fill Scattering Chamber with Vacuum
  // G4Tubs *chamber_inner = new G4Tubs("chamber_inner", 0.0, swallrad_in,  sheight/2, 0.0*deg, 360.0*deg );
  // G4LogicalVolume* chamber_inner_log = new G4LogicalVolume(chamber_inner, GetMaterial("Vacuum"), "cham_inner_log");

  // Top and bottom
  G4Tubs *sc_topbottom = new G4Tubs("scham_topbottom", 0.0, swallrad, (swallrad-swallrad_in)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume* sc_topbottom_log = new G4LogicalVolume(sc_topbottom, GetMaterial("Aluminum"), "scham_topbottom_log");

  //////////////////////////////////////////////////////////

  // // Scattering chamber
  // new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), swall_log,
  // 		    "scham_wall_phys", worldlog, false, 0);

  // new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), chamber_inner_log,
  // 		    "chamber_inner_phys", worldlog, false, 0);

  // new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), sc_hcalwin_log,
  // 		    "sc_hcalwin_phys", worldlog, false, 0);

  new G4PVPlacement(rot_temp, G4ThreeVector(0.0, sheight/2.0 + (swallrad-swallrad_in)/2, 0.0), sc_topbottom_log,
  		    "scham_top_phys", worldlog, false, 0);

  new G4PVPlacement(rot_temp, G4ThreeVector(0.0, -sheight/2.0 - (swallrad-swallrad_in)/2, 0.0), sc_topbottom_log,
  		    "scham_bot_phys", worldlog, false, 0);

  // new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, extpipestart-extpipe_len/2), extpipe_log, "extpipe_phys", worldlog, false, 0);
  // new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, extpipestart-extpipe_len/2), extvac_log, "extvacpipe_phys", worldlog, false, 0);

  //Make exit beam pipe and vacuum:
  G4Tubs *exit_pipe = new G4Tubs( "exit_pipe", Rin_dvcs_beampipe, Rout_dvcs_beampipe, extpipe_len/2.0, 0.0, 360.*deg );
  G4Tubs *exit_vacuum = new G4Tubs( "exit_vacuum", 0.0, Rin_dvcs_beampipe, extpipe_len/2.0, 0.0, 360.*deg );

  //Cut the exit pipe and exit vacuum using dvcs vacuum snout:
  G4double z0_exitpipe = 0.5*(extpipestart + Rin_dvcs_snout);

  G4RotationMatrix rinv = rot_temp->inverse();
  
  G4SubtractionSolid *exit_pipe_cut = new G4SubtractionSolid( "exit_pipe_cut", exit_pipe, dvcs_snout_vacuum, &rinv, G4ThreeVector(0,0,-z0_exitpipe) );

  G4SubtractionSolid *exit_vacuum_cut = new G4SubtractionSolid( "exit_vacuum_cut", exit_vacuum, dvcs_snout_vacuum, &rinv, G4ThreeVector(0,0,-z0_exitpipe ) );

  G4LogicalVolume *exit_pipe_log = new G4LogicalVolume( exit_pipe_cut, GetMaterial("Aluminum"), "exit_pipe_log" );
  G4LogicalVolume *exit_vacuum_log = new G4LogicalVolume( exit_vacuum_cut, GetMaterial("Vacuum"), "exit_vacuum_log" );
  new G4PVPlacement( 0, G4ThreeVector(0,0,z0_exitpipe ), exit_pipe_log, "exit_pipe_phys", worldlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,z0_exitpipe ), exit_vacuum_log, "exit_vacuum_phys", worldlog, false, 0 );
  
  rot_temp = new G4RotationMatrix();
  rot_temp->rotateX(+90.0*deg);
  
  //Call BuildStandardCryoTarget HERE !
  BuildStandardCryoTarget(scham_vacuum_log, rot_temp, G4ThreeVector(0, +0.025*mm, 0));
  /*
  // Now let's make a cryotarget:
  G4double Rcell = 4.0*cm;
  G4double uthick = 0.1*mm;
  G4double dthick = 0.15*mm;
  G4double sthick = 0.2*mm;

  G4Tubs *TargetMother_solid = new G4Tubs( "TargetMother_solid", 0, Rcell + sthick, (fTargLen+uthick+dthick)/2.0, 0.0, twopi );
  G4LogicalVolume *TargetMother_log = new G4LogicalVolume( TargetMother_solid, GetMaterial("Vacuum"), "TargetMother_log" );
  
  G4Tubs *TargetCell = new G4Tubs( "TargetCell", 0, Rcell, fTargLen/2.0, 0, twopi );

  G4LogicalVolume *TargetCell_log;

  if( fTargType == kLH2 ){
    TargetCell_log = new G4LogicalVolume( TargetCell, GetMaterial("LH2"), "TargetCell_log" );
  } else {
    TargetCell_log = new G4LogicalVolume( TargetCell, GetMaterial("LD2"), "TargetCell_log" );
  }

  G4Tubs *TargetWall = new G4Tubs("TargetWall", Rcell, Rcell + sthick, fTargLen/2.0, 0, twopi );

  G4LogicalVolume *TargetWall_log = new G4LogicalVolume( TargetWall, GetMaterial("Al"), "TargetWall_log" );

  G4Tubs *UpstreamWindow = new G4Tubs("UpstreamWindow", 0, Rcell + sthick, uthick/2.0, 0, twopi );
  G4Tubs *DownstreamWindow = new G4Tubs("DownstreamWindow", 0, Rcell + sthick, dthick/2.0, 0, twopi );

  G4LogicalVolume *uwindow_log = new G4LogicalVolume( UpstreamWindow, GetMaterial("Al"), "uwindow_log" );
  G4LogicalVolume *dwindow_log = new G4LogicalVolume( DownstreamWindow, GetMaterial("Al"), "dwindow_log" );

  // Now place everything:

  G4double ztemp = -(fTargLen+uthick+dthick)/2.0;
  // Place upstream window:
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+uthick/2.0), uwindow_log, "uwindow_phys", TargetMother_log, false, 0 );
  // Place target and side walls:
  ztemp += uthick;
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+fTargLen/2.0), TargetCell_log, "TargetCell_phys", TargetMother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+fTargLen/2.0), TargetWall_log, "TargetWall_phys", TargetMother_log, false, 0 );
  ztemp += fTargLen;
  new G4PVPlacement( 0, G4ThreeVector(0,0,ztemp+dthick/2.0), dwindow_log, "dwindow_phys", TargetMother_log, false, 0 );
  
  // Place target at origin of the Scattering Chamber 
  G4RotationMatrix *targrot = new G4RotationMatrix;
  targrot->rotateX(90.0*deg);
  new G4PVPlacement( targrot, G4ThreeVector(0.0, 0.0,0.0), TargetMother_log, "TargetMother_phys", scham_vacuum_log, false, 0 );

  if( fFlux ){ //Make a sphere to compute particle flux:
    G4Sphere *fsph = new G4Sphere( "fsph", 1.5*fTargLen/2.0, 1.5*fTargLen/2.0+cm, 0.0*deg, 360.*deg,
    0.*deg, 150.*deg );
    G4LogicalVolume *fsph_log = new G4LogicalVolume( fsph, GetMaterial("Air"), "fsph_log" );
    new G4PVPlacement( targrot, G4ThreeVector(0,0,0), fsph_log, "fsph_phys", scham_vacuum_log, false, 0 );
    
    G4String FluxSDname = "FLUX";
    G4String Fluxcollname = "FLUXHitsCollection";
    G4SBSCalSD *FluxSD = NULL;
    if( !( FluxSD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(FluxSDname) ) ){
    G4cout << "Adding FLUX SD to SDman..." << G4endl;
    FluxSD = new G4SBSCalSD( FluxSDname, Fluxcollname );
    fDetCon->fSDman->AddNewDetector( FluxSD );
    (fDetCon->SDlist).insert( FluxSDname );
    fDetCon->SDtype[FluxSDname] = kCAL;
    
    (FluxSD->detmap).depth = 0;
    }
    fsph_log->SetSensitiveDetector( FluxSD );
    fsph_log->SetUserLimits(  new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
    }
  */
  
  //  VISUALS
  G4VisAttributes * schamVisAtt = new G4VisAttributes(G4Colour(0.7,0.7,1.0));
  schamVisAtt->SetForceWireframe(true);
  scham_wall_log->SetVisAttributes(schamVisAtt);
  scham_vacuum_log->SetVisAttributes(G4VisAttributes::Invisible);
  sc_topbottom_log->SetVisAttributes(schamVisAtt);
  dvcs_snout_vacuum_log->SetVisAttributes( schamVisAtt );
  // chamber_inner_log->SetVisAttributes(G4VisAttributes::Invisible);
  // sc_entry_hole_vacuum_log->SetVisAttributes( G4VisAttributes::Invisible );
  // sc_exit_hole_vacuum_log->SetVisAttributes( G4VisAttributes::Invisible );

  // G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
  // extpipe_log->SetVisAttributes(pipeVisAtt);
  // extvac_log->SetVisAttributes( G4VisAttributes::Invisible );

  // G4VisAttributes *winVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  // sc_hcalwin_log->SetVisAttributes(winVisAtt);
}

void G4SBSTargetBuilder::BuildGasTarget(G4LogicalVolume *worldlog){

  if( fFlux ){ //Make a sphere to compute particle flux:
    G4Sphere *fsph = new G4Sphere( "fsph", 1.5*fTargLen/2.0, 1.5*fTargLen/2.0+cm, 0.0*deg, 360.*deg,
				   0.*deg, 150.*deg );
    G4LogicalVolume *fsph_log = new G4LogicalVolume( fsph, GetMaterial("Air"), "fsph_log" );
    new G4PVPlacement( 0, G4ThreeVector(0,0,0), fsph_log, "fsph_phys", worldlog, false, 0 );

    fsph_log->SetUserLimits(  new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
    
    G4String FluxSDname = "FLUX";
    G4String Fluxcollname = "FLUXHitsCollection";
    G4SBSCalSD *FluxSD = NULL;
    if( !( FluxSD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(FluxSDname) ) ){
      G4cout << "Adding FLUX SD to SDman..." << G4endl;
      FluxSD = new G4SBSCalSD( FluxSDname, Fluxcollname );
      fDetCon->fSDman->AddNewDetector( FluxSD );
      (fDetCon->SDlist).insert( FluxSDname );
      fDetCon->SDtype[FluxSDname] = kCAL;

      (FluxSD->detmap).depth = 0;
    }
    fsph_log->SetSensitiveDetector( FluxSD );
  }
  
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
