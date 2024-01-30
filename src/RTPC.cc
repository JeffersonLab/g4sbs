// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class RTPC
// Geometry and materials of the RTPC
// 03/05/14 JRMA
// 11/01/16 JRMA Update more "realistic" geometry

#include "RTPC.hh"
#include "DetectorConstruction.hh"
#include "G4Material.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"
#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_SpinEqRhs.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ClassicalRK4.hh"
#include "G4NystromRK4.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
//
#include "G4CutTubs.hh"
#include "G4Cons.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"
#include "G4FieldManager.hh"
#include "G4MagneticField.hh"
#include "TOSCAField3D.hh"
#include "TOSCAField2D.hh"
#include "BonusBField.hh"
#include "Field3D.hh"
#include "ArraySD.hh"
#include "BeamLine.hh"
//#include "TH3D.h"
enum {EUniField, ENonUniField, EToscaField3D, EBonusBField, EToscaField2D};

//----------------------------------------------------------------------------
RTPC::RTPC(DetectorConstruction* rectagg)
{
  fRtag = rectagg;
  fMaw = NULL;
  //  fPWT=NULL;
  fIsInteractive=1;
  fIsSrcPb = false;
  fIsOverlapVol = false;
  fNtarget = 1;
  fNgas = 2;
  fXmin = fBmin = fBmax = fBz = 0.0;
  fTx = fTy = fTz = fRst = fTst = fZst = fRHe1 = fRHe1a = fRHe2 = fRw = fRbl =
    fTbl = fSbl = fZbl =
    fRbaf = fTbaf = fRG = fTG = fTend1 = fTend2 = fTsh = fZsh = 0;
  fNwI = fNwO = fNbl = 0;
  fShI = fShO = fShZI = fShZO = fShTh = 0.0;
  fFieldMap = NULL;
  fBFieldType = EUniField;
  fBScaleFac = 1.0;
  fMst = 0; // 0 = kapton, 1 = Al
  fOvHe = 0.0;
  fWMat = 0;
}

//----------------------------------------------------------------------------
RTPC::~RTPC()
{
 
}

//---------------------------------------------------------------------------
G4int RTPC::ReadParameters(G4String file)
{
  //
  // Initialise magnet setup
  //
  DefaultInit();
  char* keylist[] = { (char *)"Bfield-Box:", (char *)"Target-Dim:", 
		      (char *)"Uni-Field:", (char *)"Gas-Volume:", 
		      (char *)"End-Foil:", (char *)"Target-Mat:",
		      (char *)"Gas-Mat:", (char *)"BL-Dim:",
		      (char *)"GEM-Dim:",(char *)"Non-Uni-Field:",
		      (char *)"Tosca-Field:",(char *)"Baffle-Dim:",
		      (char *)"GEM-Pixels:", (char*)"Outer-Chamber:",
		      (char*)"Be-Window:", (char*)"Beam-Colli:",
		      (char*)"Bonus-Field:", (char*)"Tosca-2DField:",NULL };
  enum { EBfieldBox, ETargDim, EMagField, EGasVol, EEndFoil, ETargMat, EGasMat,
	 EBLDim, EGEMDim, ENonUniField, EToscaField, EBaffleDim, EGEMPixels,
	 EOutCham, EBeWind, EBeamColli, EBonusField, ETosca2Field,ENULL };
  G4int ikey, iread, ierr;
  ierr = 0;
  char line[256];
  char delim[64];
  FILE* pdata;
  if( (pdata = fopen(file.data(),"r")) == NULL ){
    printf("Error opening source parameter file: %s\n",file.data());
    return -1;
  }
  while( fgets(line,256,pdata) ){
    if( line[0] == '#' ) continue;       // comment #
    sscanf(line,"%s",delim);
    for(ikey=0; ikey<ENULL; ikey++)
      if(!strcmp(keylist[ikey],delim)) break;
    switch( ikey ){
    default:
      printf("Unknown setup key\n");
      return -1;
      break;
    case EBfieldBox:
      // 1/2 dimensions of magnetic field region
      iread = sscanf(line,"%*s%lf%lf%lf",&fTx,&fTy,&fTz);
      if( iread != 3 ) ierr++;
      break;
    case ETargDim:
      // dimensions target straw
      iread = sscanf(line,"%*s%lf%lf%lf%lf",&fRst,&fZst,&fTst,&fMst);
      if( iread < 3 ) ierr++;
      break;
    case EMagField:
      // Z component of uniform magnetic field
      iread = sscanf(line,"%*s%lf",&fBz);
      if( iread != 1 ) ierr++;
      fBFieldType = EUniField;
      break;
    case EGasVol:
      // Radii of active gas region
      iread =
	sscanf(line,"%*s%lf%lf%lf%d%lf%d%d%lf",
	       &fRHe1,&fRHe1a,&fRHe2,&fIsSep,&fRw,&fNwI,&fNwO,&fOvHe);
      if( iread < 8 ) ierr++;
      break;
    case EOutCham:
      // Outer Al shell, inner/outer radius, inner/outer length, thickness
      iread =
	sscanf(line,"%*s%lf%lf%lf%lf%lf",&fShI,&fShO,&fShZI,&fShZO,&fShTh);
      if( iread != 5 ) ierr++;
      break;
    case EBeWind:
      // Be Window
      iread =
	sscanf(line,"%*s%lf%d",&fWTh,&fWMat);
      if( iread < 1 ) ierr++;
      break;
    case EBeamColli:
      // Beam collimators
      iread =
	sscanf(line,"%*s%lf%lf%lf%lf",&fBmClen1,&fBmCr1,&fBmClen2,&fBmCr2);
      if( iread != 4 ) ierr++;
      break;
    case EEndFoil:
      // End cover
      iread = sscanf(line,"%*s%lf%lf",&fTend1,&fTend2);
      if( iread != 2 ) ierr++;
      break;
    case ETargMat:
      // Target material and density
      iread = sscanf(line,"%*s%d%lf",&fNtarget,&fTDens);
      if( iread != 2 ) ierr++;
      break;
    case EGasMat:
      // He density
      iread = sscanf(line,"%*s%d%lf",&fNgas,&fHeDens);
      if( iread != 2 ) ierr++;
      break;
    case EBLDim:
      // Beam line radii
      //iread = sscanf(line,"%*s%lf%lf%d%lf%lf",&fRblI,&fRblO,&fNbl,&fTsh,&fZsh);
      iread = sscanf(line,"%*s%lf%lf%d%lf%lf",&fRbl,&fTbl,&fNbl,&fSbl,&fZbl);
      if( iread != 5) ierr++;
      break;
    case EGEMDim:
      // Be interior
      iread = sscanf(line,"%*s%lf%lf",&fRG,&fTG);
      if( iread != 2 ) ierr++;
      break;
    case ENonUniField:
      // Non uniform magnetic field
      iread = sscanf(line,"%*s%lf%lf%lf",&fXmin,&fBmin,&fBmax);
      if( iread != 3 ) ierr++;
      fBFieldType = ENonUniField;
      break;
    case EToscaField:
      // Tosca magnetic field
      fFieldMap = new char[64];
      iread = sscanf(line,"%*s%lf%lf%lf%s%lf",
		     &fTXoff,&fTYoff,&fTZoff,fFieldMap,&fBScaleFac);
      if( iread != 5 ) ierr++;
      fBFieldType = EToscaField3D;
      break;
    case ETosca2Field:
      // Tosca magnetic field
      fFieldMap = new char[64];
      iread = sscanf(line,"%*s%lf%lf%lf%s%lf",
		     &fTXoff,&fTYoff,&fTZoff,fFieldMap,&fBScaleFac);
      if( iread != 5 ) ierr++;
      fBFieldType = EToscaField2D;
      break;
    case EBonusField:
      // Bonus magnetic field
      fFieldMap = new char[64];
      iread = sscanf(line,"%*s%lf%lf%lf%s%lf",
		     &fTXoff,&fTYoff,&fTZoff,fFieldMap,&fBScaleFac);
      if( iread != 5 ) ierr++;
      fBFieldType = EBonusBField;
      break;
    case EBaffleDim:
      // Baffle dimensions
      iread = sscanf(line,"%*s%lf%lf",&fRbaf,&fTbaf);
      if( iread != 2 ) ierr++;
      break;
    case EGEMPixels:
      // no. pixels in z,phi for GEM pads
      iread = sscanf(line,"%*s%d%d",&fZPixG,&fPhiPixG);
      if( iread != 2 ) ierr++;
      break;
    }
  }
  if( ierr ) printf("Error detected in procedure RTPC::ReadParameters: %s\n",
		    line);
  return ierr;
}

//----------------------------------------------------------------------------
G4VPhysicalVolume* RTPC::Construct(G4LogicalVolume* maw){
  // Main tank build
  // Then build view ports and target
  //
  fMaw = maw;             // save the mother volume
  fZHe = fZst + fOvHe;
  //
  //G4Box* Bfield = new G4Box("Bfield", fTx,fTy,fTz);
  G4Tubs* Bfield = new G4Tubs("Bfield", 0.0,fTy,fTz,0.0,360*deg);
  //
  fLBfield =
    new G4LogicalVolume(Bfield, fRtag->GetAir(),"LBfield",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),fLBfield,"PBfield",fMaw,false,0,
		    fIsOverlapVol);
  //fLBfield->SetVisAttributes (G4VisAttributes::Invisible);
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  G4double parm[5];
  parm[0] = fZst;
  parm[1] = fRHe1;
  parm[2] = fRHe2;
  parm[3] = (G4double)fZPixG;
  parm[4] = (G4double)fPhiPixG;
  fArrSD = new ArraySD("ArraySD", 9, parm );
  SDman->AddNewDetector(fArrSD);
  //
  ConEndCap();
  ConShell();  
  ConTarget();
  ConGas();
  ConWire();
  ConGEM();
  ConBeamLine();
  Beamline* bl = new Beamline(fRtag);
  bl->BuildComponent(fMaw);  
  SetMagField();
  return NULL;
}
//----------------------------------------------------------------------------
void RTPC::ConGEM()
{
  // Cylindrical triple GEM chamber surrounding He gas volume of RTPC
  //
  G4double rg = fRHe2;       // start at outer radius of He volume
  G4double tgkap = 0.05*mm;  // 50um kapton GEM foils
  G4double tgkap1 = 0.2*mm;  // 200um kapton readout foil
  G4double tgcu = 0.005*mm;  // 5 um Cu cladding of foils
  G4double tgap = 2.0*mm;    // radial spacing between foils
  G4double rginner = rg;     // inner GEM radius
  G4double rgouter;
  // Inner GEM
  G4Tubs* G1a = new G4Tubs("G1a",rg,rg+tgcu,fZHe,0.0,360*deg);
  G4LogicalVolume* LG1a =
    new G4LogicalVolume(G1a,fRtag->GetCu(),"LG1a",0,0,0);
  rg += tgcu;
  G4Tubs* G1b = new G4Tubs("G1b",rg,rg+tgkap,fZHe,0.0,360*deg);
  G4LogicalVolume* LG1b =
    new G4LogicalVolume(G1b,fRtag->GetKapton(),"LG1b",0,0,0);
  rg += tgkap;
  G4Tubs* G1c = new G4Tubs("G1c",rg,rg+tgcu,fZHe,0.0,360*deg);
  G4LogicalVolume* LG1c =
    new G4LogicalVolume(G1c,fRtag->GetCu(),"LG1c",0,0,0);
  LG1a->SetVisAttributes(G4VisAttributes::Invisible);
  LG1b->SetVisAttributes(G4VisAttributes::Invisible);
  LG1c->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
  // Middle GEM
  rg = fRHe2 + tgap;
  G4Tubs* G2a = new G4Tubs("G2a",rg,rg+tgcu,fZHe,0.0,360*deg);
  G4LogicalVolume* LG2a =
    new G4LogicalVolume(G2a,fRtag->GetCu(),"LG2a",0,0,0);
  rg += tgcu;
  G4Tubs* G2b = new G4Tubs("G2b",rg,rg+tgkap,fZHe,0.0,360*deg);
  G4LogicalVolume* LG2b =
    new G4LogicalVolume(G2b,fRtag->GetKapton(),"LG2b",0,0,0);
  rg += tgkap;
  G4Tubs* G2c = new G4Tubs("G2c",rg,rg+tgcu,fZHe,0.0,360*deg);
  G4LogicalVolume* LG2c =
    new G4LogicalVolume(G2c,fRtag->GetCu(),"LG2c",0,0,0);
  LG2a->SetVisAttributes(G4VisAttributes::Invisible);
  LG2b->SetVisAttributes(G4VisAttributes::Invisible);
  LG2c->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
  // Outer GEM
  rg = fRHe2 + 2*tgap;
  G4Tubs* G3a = new G4Tubs("G3a",rg,rg+tgcu,fZHe,0.0,360*deg);
  G4LogicalVolume* LG3a =
    new G4LogicalVolume(G3a,fRtag->GetCu(),"LG3a",0,0,0);
  rg += tgcu;
  G4Tubs* G3b = new G4Tubs("G3b",rg,rg+tgkap,fZHe,0.0,360*deg);
  G4LogicalVolume* LG3b =
    new G4LogicalVolume(G3b,fRtag->GetKapton(),"LG3b",0,0,0);
  rg += tgkap;
  G4Tubs* G3c = new G4Tubs("G3c",rg,rg+tgcu,fZHe,0.0,360*deg);
  G4LogicalVolume* LG3c =
    new G4LogicalVolume(G3c,fRtag->GetCu(),"LG3c",0,0,0);
  LG3a->SetVisAttributes(G4VisAttributes::Invisible);
  LG3b->SetVisAttributes(G4VisAttributes::Invisible);
  LG3c->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
  // Redaout Foil
  rg = fRHe2 + 3*tgap;
  G4Tubs* G4a = new G4Tubs("G4a",rg,rg+tgcu,fZHe,0.0,360*deg);
  G4LogicalVolume* LG4a =
    new G4LogicalVolume(G4a,fRtag->GetCu(),"LG4a",0,0,0);
  rg += tgcu;
  G4Tubs* G4b = new G4Tubs("G4b",rg,rg+tgkap1,fZHe,0.0,360*deg);
  G4LogicalVolume* LG4b =
    new G4LogicalVolume(G4b,fRtag->GetKapton(),"LG4b",0,0,0);
  rg += tgkap1;
  G4Tubs* G4c = new G4Tubs("G4c",rg,rg+tgcu,fZHe,0.0,360*deg);
  rgouter = rg + tgcu;
  //
  G4LogicalVolume* LG4c =
    new G4LogicalVolume(G4c,fRtag->GetCu(),"LG4c",0,0,0);
  LG4a->SetVisAttributes(G4VisAttributes::Invisible);
  LG4b->SetVisAttributes(G4VisAttributes::Invisible);
  LG4c->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
  //
  // Volume to hold all GEM components...fill with RTPC gas
  G4Tubs* Gall = new G4Tubs("Gall",rginner,rgouter,fZHe,0.0,360*deg);
  G4LogicalVolume* LGall = 
    new G4LogicalVolume(Gall,fMgas,"LGall",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
  		    LGall,"PGall",fLrtpc,false,2,fIsOverlapVol);
  LGall->SetVisAttributes(G4VisAttributes::Invisible);
  LGall->SetSensitiveDetector( fArrSD );

  //
  // Place all layers of GEM foil and readout in GEM volume
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG1a,"PG1a",LGall,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG1b,"PG1b",LGall,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG1c,"PG1c",LGall,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG2a,"PG2a",LGall,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG2b,"PG2b",LGall,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG2c,"PG2b",LGall,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG3a,"PG3a",LGall,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG3b,"PG3b",LGall,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG3c,"PG3c",LGall,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG4a,"PG4a",LGall,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG4b,"PG4b",LGall,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG4c,"PG4c",LGall,false,0,fIsOverlapVol);
}
//----------------------------------------------------------------------------
void RTPC::ConTarget()
{
  // Construct straw target at centre of RTPC
  // Target tube
  G4Tubs* Straw = new G4Tubs("Straw",fRst-fTst,fRst,fZHe,0.0,360*deg);
  // Target gas
  G4Tubs* Straw1 = new G4Tubs("Straw1",0.0,fRst-fTst,fZst,0.0,360*deg);
  // Target window
  G4Tubs* WindT = new G4Tubs("WindT",0.0,fRst-fTst,fWTh,0.0,360*deg);
  G4Material* targWall;
  if( fMst == 0 )
    targWall = fRtag->GetKapton();
  else
    targWall = fRtag->GetAl();
  G4LogicalVolume* LStraw =
    new G4LogicalVolume(Straw, targWall,"LStraw",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),LStraw,"PStraw",fLrtpc,false,0,
		    fIsOverlapVol);
  G4Material* target;
  if(fNtarget == 1) target = fRtag->GetCoH2(fTDens);
  else target = fRtag->GetCoD2(fTDens);
  G4cout << target << endl;
  G4LogicalVolume* LStraw1 =
    new G4LogicalVolume(Straw1, target, "LStraw1 ",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,-fOvHe),
		    LStraw1,"PStraw1",fLrtpc,false,3,fIsOverlapVol);
  G4LogicalVolume* LWindT =
	new G4LogicalVolume(WindT,fRtag->GetBe(),"LWindT",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,-fZst+2*fWTh),
		    LWindT,"PWindT1",LStraw1,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,fZst-2*fWTh),
		    LWindT,"PWindT2",LStraw1,false,0,fIsOverlapVol);
  LStraw->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
  LStraw1->SetVisAttributes(new G4VisAttributes(G4Colour::Blue()));
  LStraw1->SetSensitiveDetector( fArrSD );
}
//----------------------------------------------------------------------------
void RTPC::ConWire()
{
  // Construct electric field wires
  //
  G4Tubs* WireI = new G4Tubs("WireI",0.0,fRw,fZHe-2.275,0.0,360*deg);
  G4Tubs* WireO = new G4Tubs("WireO",0.0,fRw,fZHe,0.0,360*deg);
    // Field Wires....fNw > 0
  G4LogicalVolume* LWireI;
  G4LogicalVolume* LWireO;
  G4double r,dth;
  char wnm[32];
  if(fNwI > 0){
    LWireI = new G4LogicalVolume(WireI, fRtag->GetW(),"LWireI",0,0,0);
    LWireI->SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
    r = fRHe1+fRw;
    dth = 360*deg/fNwI;
    for(G4int i=0; i<fNwI; i++){
      G4double th = dth * (i-1);
      G4double x = r*cos(th);
      G4double y = r*sin(th);
      sprintf(wnm,"PwI%d",i);
      new G4PVPlacement(0,G4ThreeVector(x,y,0),
			LWireI,wnm,fLHe2,false,0,fIsOverlapVol);
    }
  }
  if(fNwO > 0){
    LWireO = new G4LogicalVolume(WireO, fRtag->GetW(),"LWireO",0,0,0);
    LWireO->SetVisAttributes(new G4VisAttributes(G4Colour::Green()));
    r = fRHe1a+fRw;
    dth = 360*deg/fNwO;
    for(G4int i=0; i<fNwO; i++){
      G4double th = dth * (i-1);
      G4double x = r*cos(th);
      G4double y = r*sin(th);
      sprintf(wnm,"PwO%d",i);
      new G4PVPlacement(0,G4ThreeVector(x,y,0),
			LWireO,wnm,fLHe2,false,0,fIsOverlapVol);
    }
  }
}
//----------------------------------------------------------------------------
void RTPC::ConGas()
{
  // RTPC ionising gas
  G4double rHe2a = 90.0*mm;
  G4double rHe2b = 50.0*mm;
  G4double zHe2a = fZHe - 2.275*mm;
  G4double zHe2b = fZHe - 4.0*mm;
  //
  G4Tubs* He2 = new G4Tubs("He2",rHe2a,fRHe2,fZHe,0.0,360*deg);
  G4Tubs* He2a = new G4Tubs("He2a",rHe2b,rHe2a,zHe2a,0.0,360*deg);
  G4Tubs* He2b = new G4Tubs("He2b",fRst,rHe2b,zHe2b,0.0,360*deg);
  G4UnionSolid* uHe2a = new G4UnionSolid("He2a",He2,He2a);
  G4UnionSolid* uHe2b = new G4UnionSolid("He2b",uHe2a,He2b);
  //
  //G4Material* gas;
  if(fNgas == 1) fMgas = fRtag->GetCoH2(fHeDens);
  else fMgas = fRtag->GetHe(fHeDens);
  G4cout << fMgas << endl;
  fLHe2 = 
    new G4LogicalVolume(uHe2b, fMgas,"LHe2",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    fLHe2,"PHe2",fLrtpc,false,0,fIsOverlapVol);
  // Gas invisible and sensitive volume
  //fLHe2->SetVisAttributes(G4VisAttributes::Invisible);
  fLHe2->SetSensitiveDetector( fArrSD );

}
//----------------------------------------------------------------------------
void RTPC::ConBeamLine()
{
  // Construct electron beam line
  //
  //G4double Zbl = 0.5*(fTz - fZrtpc);
  //G4double ZblOff = -fTz + Zbl;
  G4double Zbli = 0.5*(fZbl - fZrtpc);
  G4double Zble = Zbli/fNbl;
  //G4double ZblOff = fZrtpc + Zble;
  //G4double Zbl2 = 2*Zbl;
  //fZbl = Zbl/fNbl - fTbl;

  G4int nbl = 3*fNbl;
  G4Tubs** Bl = new G4Tubs*[nbl];
  G4LogicalVolume** LBl = new G4LogicalVolume*[nbl];
  char name[64];
  G4Tubs* BlinI = new G4Tubs("BLinI",0.0,fRbl,Zbli,0.0,360*deg);
  G4Tubs* BlinO = new G4Tubs("BLinO",0.0,fRbl+fTbl,Zbli,0.0,360*deg);
  G4Tubs* Colli1 = new G4Tubs("Colli1",fBmCr1,fRbl,fBmClen1,0.0,360*deg);
  G4Tubs* Colli2 = new G4Tubs("Colli2",fBmCr2,fRbl,fBmClen2,0.0,360*deg);
  G4double rbl = fRbl;
  for(G4int ibl=0; ibl<nbl; ibl+=3){
    sprintf(name,"BlO-%d",ibl);
    Bl[ibl] = new G4Tubs(name,0.0,rbl+fTbl,Zble,0.0,360*deg);
    sprintf(name,"BlI-%d",ibl);
    Bl[ibl+1] = new G4Tubs(name,0.0,rbl,Zble,0.0,360*deg);
    sprintf(name,"BlE-%d",ibl);
    Bl[ibl+2] = new G4Tubs(name,rbl+fTbl,rbl+fTbl+fSbl,fTbl,0.0,360*deg);
    rbl += fSbl;
  }
  // Create vacuum and Al beam pipe logical volumes
  // Put vacuum inside Al. 
  G4LogicalVolume* LBlinI = 
    new G4LogicalVolume(BlinI,fRtag->GetVac(),"LBlinI",0,0,0);
  G4LogicalVolume* LBlinO =
    new G4LogicalVolume(BlinO,fRtag->GetAl(),"LBlinO",0,0,0);
  G4LogicalVolume* LColli1 =
    new G4LogicalVolume(Colli1,fRtag->GetW(),"LColli1",0,0,0);
  G4LogicalVolume* LColli2 =
    new G4LogicalVolume(Colli2,fRtag->GetW(),"LColli2",0,0,0);
  //
  new G4PVPlacement(0,G4ThreeVector(0,0,0),LColli1,"PColli1",LBlinI,false,
		    0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),LBlinI,"PBlin",LBlinO,false,
		    0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,+Zbli-fBmClen2),LColli2,
		    "PColli2",LBlinI,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,-(fZrtpc+Zbli)),LBlinO,"PBlU",fLBfield,
		    false,0,fIsOverlapVol);
  //
  for(G4int ibl=0; ibl<nbl; ibl+=3){
    sprintf(name,"LBlI-%d",ibl);
    LBl[ibl+1] = new G4LogicalVolume(Bl[ibl+1],fRtag->GetVac(),name,0,0,0);
    sprintf(name,"LBlO-%d",ibl);
    LBl[ibl] = new G4LogicalVolume(Bl[ibl],fRtag->GetAl(),name,0,0,0);
    sprintf(name,"LBlE-%d",ibl+1);
    LBl[ibl+2] = new G4LogicalVolume(Bl[ibl+2],fRtag->GetAl(),name,0,0,0);
    sprintf(name,"PB-%d",ibl);
    new G4PVPlacement(0,G4ThreeVector(0,0,0),LBl[ibl+1],name,LBl[ibl],false,
		    0,fIsOverlapVol);
  }
  // Place beam-line elements
  //G4double zz = fZst + fZbl;
  //G4double zz = fShZO + fZbl;
  G4double zz = fZrtpc + Zble;
  for(G4int ibl=0; ibl<nbl; ibl+=3){
    G4LogicalVolume* lbf = fLBfield;
    //if(ibl >= nbl/2) lbf = fMaw;
    sprintf(name,"PBl-%d",ibl);
    new G4PVPlacement(0,G4ThreeVector(0,0,zz),LBl[ibl],name,lbf,false,
		      0,fIsOverlapVol);
    //zz = zz + fZbl + fTbl;
    sprintf(name,"PBlE-%d",ibl);
    if(ibl+3 < nbl)
      new G4PVPlacement(0,G4ThreeVector(0,0,zz+Zble-fTbl),LBl[ibl+2],name,lbf,false,
		      0,fIsOverlapVol);
    zz = zz + 2*Zble;
  }
  // Optional downstream shield
  /*
  G4Tubs* Shield = NULL;
  G4LogicalVolume* LShield = NULL;
  if(fTsh){
    Shield = new G4Tubs("Shield",fRbl+fTbl,fRGtot,fTsh,0.0,360*deg);
    LShield =  new G4LogicalVolume(Shield,fRtag->GetAl(),"Shield",0,0,0);
    new G4PVPlacement(0,G4ThreeVector(0,0,-ZblOff+fZsh),LShield,"PShield",
		      fLBfield,false,0,fIsOverlapVol);
  }
  */
  LColli1->SetVisAttributes(new G4VisAttributes(G4Colour::Blue()));
  LColli2->SetVisAttributes(new G4VisAttributes(G4Colour::Blue()));
}
//----------------------------------------------------------------------------
void RTPC::ConShell()
{
  // Construct outer shell of RTPC
  //
  G4Tubs* ShI = new G4Tubs("ShI",fShI,fShI+fShTh,fShZI,0.0,360*deg);         
  G4Tubs* ShIE = new G4Tubs("ShIE",fRst1,fShI,fShTh,0.0,360*deg);         
  G4Tubs* ShO = new G4Tubs("ShO",fShO,fShO+fShTh,fShZO,0.0,360*deg);         
  G4Tubs* ShOE = new G4Tubs("ShOE",fRst1,fShO,fShTh,0.0,360*deg);         
  G4LogicalVolume* LShI =
    new G4LogicalVolume(ShI,fRtag->GetAl(),"LShI",0,0,0);
  G4LogicalVolume* LShIE =
    new G4LogicalVolume(ShIE,fRtag->GetAl(),"LShIE",0,0,0);
  G4LogicalVolume* LShO =
    new G4LogicalVolume(ShO,fRtag->GetAl(),"LShO",0,0,0);
  G4LogicalVolume* LShOE =
    new G4LogicalVolume(ShOE,fRtag->GetAl(),"LShOE",0,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),LShI,"PShI",fLrtpc,false,0,
		    fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,+fShZI-fShTh),LShIE,"PShIE1",fLrtpc,
		    false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,-fShZI+fShTh),LShIE,"PShIE2",fLrtpc,
		    false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),LShO,"PShO",fLrtpc,false,0,
		    fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,+fShZO-fShTh),LShOE,"PShOE1",fLrtpc,
		    false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,-fShZO+fShTh),LShOE,"PShOE2",fLrtpc,
		    false,0,fIsOverlapVol);
  LShI->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
  LShIE->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
  LShO->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
  LShOE->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
}
//----------------------------------------------------------------------------
void RTPC::ConEndCap()
{
  G4double sphrad = 1000*mm;
  G4double spht = 1*mm;
  //G4double sphth = 8.62693*deg;
  //G4double sphth = 8.5*deg; // trunc angle..do not protrude above outer ring
  G4double sphth = asin(fRHe2/(sphrad+spht));
  G4double Rcollar = fRst + 5;
  G4double sphth1 = asin(Rcollar/(sphrad+spht));
  G4double wrth = 5.0*mm;
  G4double Rcan = fRHe2 + 4*fRG;
  //G4double tcan = 3*mm;
  G4double zcan = fZHe + 2*wrth;
  fRGtot = fShO + fShTh;
  G4double zext = 20.0*mm;
  fZrtpc = fZHe + 2*zext;
  fRst1 = fRst + 1*mm;
  G4RotationMatrix* rotW = new G4RotationMatrix();
  rotW->rotateX(180*deg);
  G4Sphere* SphW =
    new G4Sphere("SphW",sphrad,sphrad+spht,0.0,360*deg,sphth1,sphth-sphth1);
  G4Tubs* rtpc = new G4Tubs("rtpc",0,fRGtot,fZrtpc,0,360*deg);
  fLrtpc =
    new G4LogicalVolume(rtpc,fRtag->GetVac(),"Lrtpc",0,0,0);
  G4Tubs* r1 = new G4Tubs("wr1",fRHe1-wrth,fRHe1+wrth,wrth,0,360*deg);
  G4Tubs* r2 = new G4Tubs("wr2",fRHe1a-wrth,fRHe1a+wrth,wrth,0,360*deg);
  G4Tubs* r3 = new G4Tubs("wr3",fRHe2-wrth,fRHe2,wrth,0,360*deg);
  G4Tubs* r4 = new G4Tubs("wr4",fRst,Rcollar,wrth,0,360*deg);
  G4Tubs* r5 = new G4Tubs("wr5",fRst-fTst,fRst1,zext,0.0,360*deg);
  G4Tubs* r6 = new G4Tubs("wr6",0.0,fRst,3*zext,0,360*deg);
  G4Tubs* r7 = new G4Tubs("wr7",fRHe2,Rcan+fTG,zcan,0,360*deg);
  G4Tubs* r8 = new G4Tubs("wr8",fRHe2,Rcan,zcan-fTG,0,360*deg);
  G4Tubs* r9 = new G4Tubs("wr9",fRst1,fRbl+fTbl,fTbl,0,360*deg);
  G4RotationMatrix* rm = new G4RotationMatrix();
  //G4ThreeVector pm1 = G4ThreeVector(0,0,0);
  G4ThreeVector pm2(0.0,0.0,-992.5*mm);
  G4RotationMatrix irm = rm->invert();
  //G4Transform3D tm1 = G4Transform3D(rm,pm1);
  G4Transform3D tm2 = G4Transform3D(irm,pm2);
  G4ThreeVector pm3(0.0,0.0,3.0);
  G4Transform3D tm3 = G4Transform3D(irm,pm3);
  G4ThreeVector pm4(0.0,0.0,-(zext-wrth));
  G4Transform3D tm4 = G4Transform3D(irm,pm4);
  G4UnionSolid* mus1 = new G4UnionSolid("EC",r1,SphW,tm2);
  G4UnionSolid* mus2 = new G4UnionSolid("EC1",mus1,r2);
  G4UnionSolid* mus3 = new G4UnionSolid("EC2",mus2,r3);
  G4UnionSolid* mus4 = new G4UnionSolid("EC3",mus3,r4,tm3);
  G4UnionSolid* mus5 = new G4UnionSolid("EC4",mus4,r5,tm4);
  G4SubtractionSolid* mus6 = new G4SubtractionSolid("EC5",mus5,r6);
  G4SubtractionSolid* mus7 = new G4SubtractionSolid("EC6",r7,r8);
  G4LogicalVolume* LSphW =
    new G4LogicalVolume(mus6,fRtag->GetAl(),"LSphW",0,0,0);
  G4LogicalVolume* Lcan =
    new G4LogicalVolume(mus7,fRtag->GetAl(),"Lcan",0,0,0);
  G4LogicalVolume* Lend =
    new G4LogicalVolume(r9,fRtag->GetAl(),"Lend",0,0,0);
  //
  new G4PVPlacement(0,G4ThreeVector(0,0,-(fZHe+wrth)),
		    LSphW,"PEC1",fLrtpc,false,0,fIsOverlapVol);
  new G4PVPlacement(rotW,G4ThreeVector(0,0,+(fZHe+wrth)),
  		    LSphW,"PEC2",fLrtpc,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
  		    Lcan,"Pcan",fLrtpc,false,0,fIsOverlapVol);
  new G4PVPlacement(0,G4ThreeVector(0,0,+fZrtpc-fTbl),
  		    Lend,"Pend",fLrtpc,false,0,fIsOverlapVol);  
  new G4PVPlacement(0,G4ThreeVector(0,0,-fZrtpc+fTbl),
  		    Lend,"Pend",fLrtpc,false,0,fIsOverlapVol);  
  //new G4PVPlacement(0,G4ThreeVector(0,370,fOvHe*mm),
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
  		    fLrtpc,"Prtpc",fLBfield,false,0,fIsOverlapVol);
  fLrtpc->SetVisAttributes (G4VisAttributes::Invisible);
}
//----------------------------------------------------------------------------
void RTPC::SetMagField( )
{
  // 2v4 use default stepper
  // original caused chord length error
  G4FieldManager  *magFieldMgr;
  G4MagneticField* magField;
  G4double pos[3];
  G4double b[3];
  pos[0] = pos[1] = pos[2] = 0.0;
  switch(fBFieldType){
  case EUniField:
  default:
    magField = new G4UniformMagField( G4ThreeVector(0,0,fBz*tesla) );
    break;
  case ENonUniField:
    magField = new Field3D(fXmin, fBmin, fBmax);
    break;
  case EToscaField3D:
    magField = new TOSCAField3D(fFieldMap, fTXoff, fTYoff, fTZoff, fBScaleFac);
    break;
  case EBonusBField:
    magField = new BonusBField(fFieldMap, fTXoff, fTYoff, fTZoff, fBScaleFac);
    /*
    for(G4int i=0; i<10; i++){
      pos[0] = 0;
      pos[1] = i*0.6;
      pos[2] = i*0.6;
      magField->GetFieldValue(pos,b);
      printf("%g %g %g %g %g\n",pos[1],pos[2],b[0],b[1],b[2]);
      //GetFieldStrength(Double_t* pos, Double_t b)
    }
    */
    break;
  case EToscaField2D:
    magField = new TOSCAField2D(fFieldMap, fTXoff, fTYoff, fTZoff, fBScaleFac);
    for(G4int i=0;i<100;i++){
      pos[0] = 0.0;
      pos[1] = 0.0;
      pos[2] = -1500 + i*30;
      magField->GetFieldValue(pos,b);
      printf("%g, %g, %g, %g\n",pos[2],b[0]/tesla,b[1]/tesla,b[2]/tesla);
    }
    for(G4int i=0;i<100;i++){
      pos[0] = -500 + i*10;
      pos[1] = 0.0;
      pos[2] = 0.0;
      magField->GetFieldValue(pos,b);
      printf("%g, %g, %g, %g\n",pos[0],b[0]/tesla,b[1]/tesla,b[2]/tesla);
    }
    break;
  }
  magField->GetFieldValue(pos,b);
  printf("Value of Magnetic Field at (0,0.01,0) is Bx:%g  By:%g Bz:%g\n",
	 b[0],b[1],b[2]);
  //
  //  G4Mag_UsualEqRhs* Equation = new G4Mag_UsualEqRhs(magField);
  //  G4MagIntegratorStepper* Stepper = new G4HelixImplicitEuler(Equation);
  G4double minStep = 0.10*mm;
  //  G4ChordFinder* chordF = new G4ChordFinder(magField,minStep,Stepper);
  //
  magFieldMgr = new G4FieldManager(magField);
  magFieldMgr->SetDetectorField(magField);
  magFieldMgr->CreateChordFinder(magField);
  //  magFieldMgr->SetChordFinder(chordF);
  magFieldMgr->GetChordFinder()->SetDeltaChord(minStep);

  fLBfield->SetFieldManager(magFieldMgr,true);
}
G4int* RTPC::GetNhits(){ return fArrSD->GetNhits(); }
G4int* RTPC::GetHitID(){ return fArrSD->GetHitID(); }
G4int* RTPC::GetHits(){ return fArrSD->GetHits(); }
