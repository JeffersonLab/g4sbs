#include "G4SBSDetectorConstruction.hh"
#include "G4SBSGEMSD.hh"
#include "G4SBSCalSD.hh"

#include "G4UserLimits.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ProductionCuts.hh"
#include "G4ExtrudedSolid.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4GenericTrap.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4TwoVector.hh"
#include "G4RotationMatrix.hh"

#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4SBSBigBiteField.hh"
#include "G4SBS48D48Field.hh"
#include "G4FieldManager.hh"

#include "G4MagIntegratorStepper.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4ChordFinder.hh"
#include "G4SimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"

#include <vector>


G4SBSDetectorConstruction::G4SBSDetectorConstruction()
{
    fBBang  = 40.0*deg;
    fBBdist = 1.5*m;

    fBBCaldist = 0.8*m;

    f48D48ang  = 39.4*deg;
    f48D48dist = 2.8*m;
    fHCALdist  = 17.0*m;

    fTargLen = 60.0*cm;

    fTargType = kH2;
    fTargDen = 10.5*atmosphere/(300*kelvin*k_Boltzmann);

    fCerDepth = 60.0*cm;
    fCerDist  =  7.0*cm;

    fGEMDist  = 70.0*cm;

    fbbfield = new G4SBSBigBiteField( fBBdist, NULL );
    
    f48d48field = NULL;

    fGEMOption = 1;

    fTotalAbs = true;
}



G4SBSDetectorConstruction::~G4SBSDetectorConstruction()
{;}

G4VPhysicalVolume* G4SBSDetectorConstruction::Construct(){
    // Just nothing so we don't step on toes
    // ConstructAll is where the real magic happens
  //G4double a, iz, z, density;
  G4double density;
  G4String name, symbol;
  G4int nel;

  G4Element *elN = new G4Element("Nitrogen", "N", 7, 14.007*g/mole );
  G4Element *elO = new G4Element("Oxygen", "O", 8, 16.000*g/mole );

  density = 1e-9*mg/cm3;
  G4Material* Air = new G4Material(name="BlandAir", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  G4Box *WorldBox= new G4Box("WorldBox",50*m, 50*m, 50*m);
  G4LogicalVolume *WorldLog=new G4LogicalVolume(WorldBox,Air,
						  "WorldLogical", 0, 0, 0);
  G4PVPlacement *WorldPhys=new G4PVPlacement(0,G4ThreeVector(),
					       "WorldPhysical",
					       WorldLog,
					       0,false,0);
  return WorldPhys;
}

G4VPhysicalVolume* G4SBSDetectorConstruction::ConstructAll()
{
  G4cout << "\nG4SBSDetectorConstruction....\n" << G4endl;
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  //--------- Material definition ---------
  
  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;
  
  a = 126.9*g/mole;
  G4Element* elI = new G4Element(name="Iodine",   symbol="I",  iz=53., a);
  a = 132.9*g/mole;
  G4Element* elCs= new G4Element(name="Cesium",   symbol="Cs", iz=55., a);

  G4Element *elH = new G4Element("Hydrogen", "H", 1, 1.007*g/mole );
  G4Element *elC = new G4Element("Carbon", "C", 6, 12.011*g/mole );
  G4Element *elN = new G4Element("Nitrogen", "N", 7, 14.007*g/mole );
  G4Element *elO = new G4Element("Oxygen", "O", 8, 16.000*g/mole );
  G4Element *elF = new G4Element("Fluorine", "F", 9, 18.998*g/mole );
  G4Element *elAl = new G4Element("Aluminum", "Al", 13, 26.982*g/mole );
  G4Element *elSi = new G4Element("Silicon", "Si", 14, 28.086*g/mole );
  G4Element *elCa = new G4Element("Calcium", "Ca", 20, 40.078*g/mole );
  G4Element *elSr = new G4Element("Strontium", "Sr", 38, 87.62*g/mole );
  G4Element *elBa = new G4Element("Barium", "Ba", 56, 137.327*g/mole );

   G4Element* elCl  = new G4Element("Chlorine",  "Cl", z=17, a=   35.453*g/mole);
   G4Element* elAr  = new G4Element("Argon",     "Ar", z=18, a=    39.95*g/mole);

  G4Material* Lead  = new G4Material(name="Lead", z=82., a=208.0*g/mole, density=11.34*g/cm3);

  G4Material* Aluminum = new G4Material(name="Aluminum", z=13., a=26.98*g/mole, density=2.7*g/cm3);


  density = 4.51*g/cm3;
  G4Material* CsI = new G4Material(name="CsI", density, nel = 2);
  CsI->AddElement(elI, .5);
  CsI->AddElement(elCs,.5);
  a = 4.0*g/mole;
  density = 0.1786e-03*g/cm3;
  
  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Fer", z=26., a, density);
  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  double bigden = 1e9*g/cm3;

  // Cell Glass - GE180 Aluminosilicate Glass
  // SiO2 60.3%
  G4Material* SiO2 = new G4Material("SiO2", 2.2*g/cm3, 2 );
  SiO2->AddElement(elSi, 1);
  SiO2->AddElement(elO, 2);
  // BaO  18.2%
  G4Material* BaO = new G4Material("BaO", bigden, 2 );
  BaO->AddElement(elBa, 1);
  BaO->AddElement(elO, 1);
  // Al2O3 14.3%
  G4Material* Al2O3 = new G4Material("Al2O3", bigden, 2 );
  Al2O3->AddElement(elAl, 2);
  Al2O3->AddElement(elO, 3);
  // CaO   6.5%
  G4Material* CaO = new G4Material("CaO", bigden, 2 );
  CaO->AddElement(elCa, 1);
  CaO->AddElement(elO, 1);
  // SrO   0.25%
  G4Material* SrO = new G4Material("SrO", bigden, 2 );
  SrO->AddElement(elSr, 1);
  SrO->AddElement(elO, 1);
  //
  density = 1.19*g/cm3;
  G4Material* Acrylic = new G4Material(name="Acrylic", density, nel=3);
  Acrylic->AddElement(elC, 5);
  Acrylic->AddElement(elH, 8);
  Acrylic->AddElement(elO, 2);
  // Density 2.76 g/cm^3
  // Index of Refraction 1.536
  G4Material* GE180 = new G4Material("GE180", 2.76*g/cm3, 5);
  GE180->AddMaterial(SiO2, 0.6039);
  GE180->AddMaterial(BaO, 0.1829);
  GE180->AddMaterial(Al2O3, 0.1439);
  GE180->AddMaterial(CaO, 0.0659);
  GE180->AddMaterial(SrO, 0.0034);

  //--------- GEM Materials  ---------
  // (stolen from GEMC)

  G4Material* NOMEX_pure = new G4Material("NOMEX_pure", density = 1.38*g/cm3, 5);
  NOMEX_pure -> AddElement(elH,0.04);
  NOMEX_pure -> AddElement(elC,0.54);
  NOMEX_pure -> AddElement(elN,0.09);
  NOMEX_pure -> AddElement(elO,0.10);
  NOMEX_pure -> AddElement(elCl,0.23);

  G4Material* NOMEX = new G4Material("NOMEX",density = 0.04*g/cm3, 2);
  NOMEX -> AddMaterial(NOMEX_pure,0.45);
  NOMEX -> AddMaterial(Air,0.55);

  G4Material* NEMAG10 = new G4Material("NEMAG10", 1.70*g/cm3, nel=4);
  NEMAG10 -> AddElement(elSi, 1);
  NEMAG10 -> AddElement(elO , 2);
  NEMAG10 -> AddElement(elC , 3);
  NEMAG10 -> AddElement(elH , 3);

  G4double density_Ar = 1.7823*mg/cm3 ;
  G4Material* Argon = new G4Material("Argon"  , density_Ar, nel=1);
  Argon->AddElement(elAr, 1);

  G4double density_CO2 = 1.977*mg/cm3;
  G4Material* CO2 = new G4Material("CO2", density_CO2, nel=2);
  CO2->AddElement(elC, 1);
  CO2->AddElement(elO, 2);

  // 1.5 Atmosphere C4F8O for cerkenkov
  G4double density_C4F8O = 9.64*mg/cm3; // density at 1ATM
  G4Material* C4F8O = new G4Material("C4F8O", density_C4F8O*1.5, nel=3);
  C4F8O->AddElement(elC, 4);
  C4F8O->AddElement(elF, 8);
  C4F8O->AddElement(elO, 1);

  G4double density_ArCO2 = .7*density_Ar + .3*density_CO2;
  // Use ArCO2
  G4Material *GEMgas= new G4Material("GEMgas", density_ArCO2, nel=2);
  GEMgas->AddMaterial(Argon, 0.7*density_Ar/density_ArCO2) ;
  GEMgas->AddMaterial(CO2, 0.3*density_CO2/density_ArCO2) ;

  G4Material *Copper= new G4Material("Copper", z=29, a=   63.55*g/mole, density = 8.96*g/cm3);

  G4Material *Kapton = new G4Material("Kapton",   density = 1.42*g/cm3, nel=4);
  Kapton->AddElement(elH, 0.026362);
  Kapton->AddElement(elC, 0.691133);
  Kapton->AddElement(elN, 0.073270);
  Kapton->AddElement(elO, 0.209235);



  //--------- G4VSolid, G4LogicalVolume, G4VPhysicalVolume  ---------
  
  //--------------
  // World:
  //--------------
  G4Box *WorldBox= new G4Box("WorldBox",20*m, 20*m, 20*m);
  G4LogicalVolume *WorldLog=new G4LogicalVolume(WorldBox,Air,
						  "WorldLogical", 0, 0, 0);
  G4PVPlacement *WorldPhys=new G4PVPlacement(0,G4ThreeVector(),
					       "WorldPhysical",
					       WorldLog,
					       0,false,0);

  //-----------------------------
  //  BigBite

  printf("BigBite at %f deg\n", fBBang/deg);

  G4RotationMatrix *bbrm = new G4RotationMatrix;
  bbrm->rotateY(fBBang);

  G4RotationMatrix *bbykrm = new G4RotationMatrix;
  bbykrm->rotateX(90.0*deg);

  double motherdepth = 600.0*cm;

  double bbmagwidth  = 1670.0*mm;
  double bbmagheight = 486.0*mm;

  double eps = 0.1*mm;

  // Mother box
  G4Box *bbmotherBox= new G4Box("bbmotherBox", bbmagwidth/2.0 + eps, 250*cm, motherdepth/2.0);

  // We need to account for offsets so we can fit BigBite and detectors in without running into
  // the target
  // Need 70 cm clearance from front of the spectrometer to front of the mother volume

  double clear = 70.0*cm;

  double motherr = fBBdist + motherdepth/2.0 - clear;
  



  std::vector<G4TwoVector> bbyokepts;
  bbyokepts.push_back( G4TwoVector( bbmagheight, 0.0*mm));
  bbyokepts.push_back( G4TwoVector( bbmagheight, 464.0*mm));
  bbyokepts.push_back( G4TwoVector( -bbmagheight, 805.0*mm));
  bbyokepts.push_back( G4TwoVector( -bbmagheight, 0.0*mm));
  bbyokepts.push_back( G4TwoVector( bbmagheight, 0.0*mm));
  bbyokepts.push_back( G4TwoVector( bbmagheight, 464.0*mm));
  bbyokepts.push_back( G4TwoVector( -bbmagheight, 805.0*mm));
  bbyokepts.push_back( G4TwoVector( -bbmagheight, 0.0*mm));

  G4GenericTrap *bbyokeTrap = new G4GenericTrap("bbyokeTrap",
	  bbmagwidth/2.0, bbyokepts );


  std::vector<G4TwoVector> bbairpts;
  bbairpts.push_back( G4TwoVector(  bbmagheight-133.1*mm, 0.0*mm - eps));
  bbairpts.push_back( G4TwoVector(  bbmagheight, 464.1*mm));
  bbairpts.push_back( G4TwoVector( -bbmagheight, 805.1*mm));
  bbairpts.push_back( G4TwoVector( -bbmagheight, 0.0*mm - eps));

  bbairpts.push_back( G4TwoVector(  bbmagheight-133.1*mm, 0.0*mm - eps));
  bbairpts.push_back( G4TwoVector(  bbmagheight, 464.1*mm));
  bbairpts.push_back( G4TwoVector( -bbmagheight, 805.1*mm));
  bbairpts.push_back( G4TwoVector( -bbmagheight, 0.0*mm - eps));

  double gapsize = 250.0*mm;

  G4GenericTrap *bbairTrap = new G4GenericTrap("bbairTrap",
	  gapsize/2.0, bbairpts );



  double coilsize = 320.0*mm;
  double coilwidth = 90.0*mm;

  std::vector<G4TwoVector> bbcoilpts;
  bbcoilpts.push_back( G4TwoVector(  bbmagheight-133.0*mm+coilsize, 0.0*mm-coilsize));
  bbcoilpts.push_back( G4TwoVector(  bbmagheight+coilsize, 464.0*mm+coilsize*(1.0-tan(20.0*deg))));
  bbcoilpts.push_back( G4TwoVector( -bbmagheight-coilsize, 805.0*mm+coilsize*(1.0+tan(20.0*deg))));
  bbcoilpts.push_back( G4TwoVector( -bbmagheight-coilsize, 0.0*mm-coilsize));

  bbcoilpts.push_back( G4TwoVector(  bbmagheight-133.0*mm+coilsize, 0.0*mm-coilsize));
  bbcoilpts.push_back( G4TwoVector(  bbmagheight+coilsize, 464.0*mm+coilsize*(1.0-tan(20.0*deg))));
  bbcoilpts.push_back( G4TwoVector( -bbmagheight-coilsize, 805.0*mm+coilsize*(1.0+tan(20.0*deg))));
  bbcoilpts.push_back( G4TwoVector( -bbmagheight-coilsize, 0.0*mm-coilsize));


  G4GenericTrap *bbcoilTrap = new G4GenericTrap("bbyokeTrap",
	 coilwidth, bbcoilpts );

  double coilrot = 40.0*mrad;

  G4RotationMatrix *coilrot1 = new G4RotationMatrix;
  coilrot1->rotateX(-coilrot );
  G4RotationMatrix *coilrot2 = new G4RotationMatrix;
  coilrot2->rotateX( coilrot );

  G4ThreeVector zcoil1 = G4ThreeVector( 0.0, 0.0, coilwidth+gapsize/2.0+(805.0*mm/2.0+coilsize)*tan(coilrot));
  G4ThreeVector zcoil2 = G4ThreeVector( 0.0, 0.0, -coilwidth-gapsize/2.0-(805.0*mm/2.0+coilsize)*tan(coilrot));
  
  // "Guitar" cuts
  std::vector<G4TwoVector> bbleftcutpts;
  bbleftcutpts.push_back( G4TwoVector(  0.0, ((bbmagwidth+gapsize)/2.0+coilwidth)/2.0 ));
  bbleftcutpts.push_back( G4TwoVector(  bbmagheight,  bbmagwidth/2.0));
  bbleftcutpts.push_back( G4TwoVector(  0.0,  bbmagwidth/2.0 + 1.0*cm));
  bbleftcutpts.push_back( G4TwoVector(  -bbmagheight*2.0/3.0,  bbmagwidth/2.0));

  bbleftcutpts.push_back( G4TwoVector(  0.0, ((bbmagwidth+gapsize)/2.0+coilwidth)/2.0 ));
  bbleftcutpts.push_back( G4TwoVector(  bbmagheight,  bbmagwidth/2.0));
  bbleftcutpts.push_back( G4TwoVector(  0.0,  bbmagwidth/2.0 + 1.0*cm));
  bbleftcutpts.push_back( G4TwoVector(  -bbmagheight*2.0/3.0,  bbmagwidth/2.0));


  G4GenericTrap *bbleftcutTrap = new G4GenericTrap("bbleftcutTrap",
	  2010.1*mm, bbleftcutpts );

  G4RotationMatrix *leftcutrot = new G4RotationMatrix;
  leftcutrot->rotateX( -90*deg );

  G4RotationMatrix *leftcutrot2 = new G4RotationMatrix;
  leftcutrot2->rotateZ( -90*deg );

  std::vector<G4TwoVector> bbrightcutpts;
  bbrightcutpts.push_back( G4TwoVector(  0.0, -((bbmagwidth+gapsize)/2.0+coilwidth)/2.0 ));
  bbrightcutpts.push_back( G4TwoVector(  bbmagheight,  -bbmagwidth/2.0));
  bbrightcutpts.push_back( G4TwoVector(  0.0, -( bbmagwidth/2.0 + 1.0*cm)));
  bbrightcutpts.push_back( G4TwoVector(  -bbmagheight*2.0/3.0,  -bbmagwidth/2.0));

  bbrightcutpts.push_back( G4TwoVector(  0.0, -((bbmagwidth+gapsize)/2.0+coilwidth)/2.0 ));
  bbrightcutpts.push_back( G4TwoVector(  bbmagheight,  -bbmagwidth/2.0));
  bbrightcutpts.push_back( G4TwoVector(  0.0, -( bbmagwidth/2.0 + 1.0*cm)));
  bbrightcutpts.push_back( G4TwoVector(  -bbmagheight*2.0/3.0,  -bbmagwidth/2.0));

  G4GenericTrap *bbrightcutTrap = new G4GenericTrap("bbrightcutTrap",
	 2010.1*mm, bbrightcutpts );

  G4RotationMatrix *rightcutrot = new G4RotationMatrix;
  rightcutrot->rotateX( 90*deg );

  G4RotationMatrix *rightcutrot2 = new G4RotationMatrix;
  rightcutrot2->rotateZ( -90*deg );

  G4UnionSolid *fullyoke;

  fullyoke = new G4UnionSolid("yoke_coil1", bbyokeTrap, bbcoilTrap, coilrot1, zcoil1 );
  fullyoke = new G4UnionSolid("yoke_coils", fullyoke, bbcoilTrap, coilrot2, zcoil2 );

  double topheight = 66.0*cm;
  G4Box *yoketopbox = new G4Box("yoketopbox", topheight/2.0, 464.0*mm/2.0, bbmagwidth/2.0);
  G4ThreeVector topboxpos = G4ThreeVector( topheight/2.0 + bbmagheight, 464.0*mm/2.0, 0.0);

  double bottomheight = 82.0*cm;
  G4Box *yokebottombox = new G4Box("yokebotbox", bottomheight/2.0, 805.0*mm/2.0, bbmagwidth/2.0);
  G4ThreeVector bottomboxpos = G4ThreeVector( -bottomheight/2.0 - bbmagheight, 805.0*mm/2.0, 0.0);

  fullyoke = new G4UnionSolid("yoke_top", fullyoke, yoketopbox, 0, topboxpos );
  fullyoke = new G4UnionSolid("yokefull", fullyoke, yokebottombox, 0, bottomboxpos );

  G4SubtractionSolid* yokewgap = new G4SubtractionSolid("yoke_with_gap", fullyoke, bbairTrap);
  yokewgap = new G4SubtractionSolid("yoke_with_gap_lcut", yokewgap, bbleftcutTrap, leftcutrot, G4ThreeVector());
  yokewgap = new G4SubtractionSolid("yoke_with_gap_lrcut", yokewgap, bbrightcutTrap, leftcutrot, G4ThreeVector());

  G4LogicalVolume *bbyokewgapLog=new G4LogicalVolume(yokewgap, Fe,
						  "bbyokewgapLog", 0, 0, 0);

  if( fTotalAbs ){
      bbyokewgapLog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }

  G4RotationMatrix *yokerm = new G4RotationMatrix;
  yokerm->rotateY(90.0*deg);
  yokerm->rotateZ(-90.0*deg);
  yokerm->rotateX(180.0*deg);

  // Cut mother box
  G4SubtractionSolid* bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxL", bbmotherBox, bbleftcutTrap, leftcutrot2, G4ThreeVector(-10*eps, 0.0, -motherdepth/2.0+clear));
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR", bbmothercutBox, bbrightcutTrap, rightcutrot2, G4ThreeVector(10*eps, 0.0, -motherdepth/2.0+clear));

  G4Box *frontboxcut = new G4Box("frontboxcut",(bbmagwidth-gapsize)/4.0-coilwidth, 250*cm + eps, clear/2);
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR_fL", bbmothercutBox, frontboxcut, 0, G4ThreeVector( ((bbmagwidth+gapsize)/4.0+coilwidth)+ 10*eps +5*cm, 0, -motherdepth/2.0+clear/2 - 10*eps ));
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR_fLR", bbmothercutBox, frontboxcut, 0, G4ThreeVector( -((bbmagwidth+gapsize)/4.0+coilwidth)-10*eps -5*cm, 0, -motherdepth/2.0+clear/2 - 10*eps ));

  

  G4Box *bottomboxcut = new G4Box("bottomboxcut",bbmagwidth+eps, 250*cm/2, motherdepth+eps);
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR_fLR_floor", bbmothercutBox, bottomboxcut, 0, G4ThreeVector( 0.0, -1.25*m*2-20*cm, 0.0));


  //   Make logical volume for mother vol and place
  G4LogicalVolume *bbmotherLog=new G4LogicalVolume(bbmothercutBox,Air,
						  "bbmotherLog", 0, 0, 0);

  new G4PVPlacement(bbrm, G4ThreeVector(-motherr*sin(fBBang), 0.0, motherr*cos(fBBang)),
	  bbmotherLog, "bbmotherPhys", WorldLog, 0,false,0);


  new G4PVPlacement(yokerm,G4ThreeVector(0.0, 0.0, -motherdepth/2.0+clear),
	  bbyokewgapLog, "bbyokewgapPhysical", bbmotherLog, 0,false,0);

  //  Bigbite field log volume
  G4LogicalVolume *bbfieldLog=new G4LogicalVolume(bbairTrap, Air,
						  "bbfieldLog", 0, 0, 0);

  fbbfield->SetRM( bbrm );
  G4FieldManager *bbfm = new G4FieldManager(fbbfield);

//  G4EqMagElectricField* fequation= new G4EqMagElectricField(fbbfield); 
  G4Mag_UsualEqRhs* fequation= new G4Mag_UsualEqRhs(fbbfield); 
  G4MagIntegratorStepper *stepper = new G4ExplicitEuler(fequation, 8);
//  G4MagIntegratorStepper *stepper = new G4ImplicitEuler(fequation, 8);
//  G4MagIntegratorStepper *stepper = new G4CashKarpRKF45(fequation);
//  G4MagIntegratorStepper *stepper = new G4SimpleRunge(fequation, 8);
    
//  G4MagInt_Driver *intgrDriver = new G4MagInt_Driver(100.0*um, stepper, stepper->GetNumberOfVariables() );
 // G4ChordFinder *chordfinder = new G4ChordFinder(fbbfield, 1.0*nm, stepper);
  new G4ChordFinder(fbbfield, 1.0*nm, stepper);
//  bbfm->SetChordFinder(chordfinder);
// bbfm->GetChordFinder()->SetDeltaChord(1.0*um);

  /*
  bbfm->SetMinimumEpsilonStep( 1e-6 );
  bbfm->SetMaximumEpsilonStep( 1e-5 );
  */

  bbmotherLog->SetFieldManager(bbfm,true);

  new G4PVPlacement(0, G4ThreeVector(), bbfieldLog, "bbfieldPhysical", bbyokewgapLog, 0,false,0);

  //--------- BigBite Detector Volumes ------------
  //
  // Mother volume is at 10 degrees with BigBite

  double detboxang    = 10.0*deg;
  double detboxheight = 2.5*m;
  double detboxdepth  = 4.0*m;
  double detboxplace  = 0.8*m; // From midplane pivot

  G4RotationMatrix *bbdetrot = new G4RotationMatrix();
  // Y is "down"
  bbdetrot->rotateZ(180.0*deg);
  bbdetrot->rotateX(-detboxang);

  // 325mm is about half the depth of the field volume
  // this is in mother box coords
  double midplanez    = -motherdepth/2.0+clear+325.0*mm;

  G4Box *bbdetbox = new G4Box("bbdetbox", bbmagwidth/2.0, detboxheight/2.0, detboxdepth/2.0);
  G4LogicalVolume *bbdetLog=new G4LogicalVolume(bbdetbox, Air,
						  "bbdetLog", 0, 0, 0);
  new G4PVPlacement(bbdetrot, G4ThreeVector(0.0, 
	               (detboxplace+detboxdepth/2.0)*sin(detboxang),
	               (detboxplace+detboxdepth/2.0)*cos(detboxang)+midplanez),
	  bbdetLog, "bbdetPhysical", bbmotherLog, 0,false,0);

  //  Just interested in the GEMs for now

  double detoffset = 0.05*m -detboxdepth/2.0;

  const int nplane = 24;

  int i;
  int ngem = 0;
  double gemdsep;
  double *gemz;
  double *gemw;
  double *gemh;



  switch( fGEMOption ){
      case 1:
	  ngem = 4;
	  gemdsep = 0.05*m;
	  break;
      case 2:
	  ngem = 5;
	  gemdsep = 0.10*m;
	  break;
      case 3:
	  ngem = 3;
	  gemdsep = 0.35*m;
	  break;
      default:
	  ngem = 4;
	  gemdsep = 0.05*m;
	  break;
  }
  //
  // GEM option 1
  double gemz_opt1[] = { 0.0*cm, gemdsep, fGEMDist, fGEMDist+gemdsep};
  double gemw_opt1[] = { 40.0*cm, 40.0*cm, 50.0*cm, 50.0*cm };
  double gemh_opt1[] = { 150.0*cm, 150.0*cm, 200.0*cm, 200.0*cm };
  
  // GEM option 2
  gemdsep = 0.15*m;
  double gemz_opt2[] = { 0.0*cm, gemdsep, 2.0*gemdsep, 3.0*gemdsep, fGEMDist};
  double gemw_opt2[] = { 40.0*cm, 40.0*cm, 40.0*cm, 40.0*cm, 50.0*cm };
  double gemh_opt2[] = { 150.0*cm, 150.0*cm, 150.0*cm, 150.0*cm, 200.0*cm };
  
  // GEM option 3
  double gemz_opt3[] = { 0.0*cm, gemdsep, gemdsep*2.0};
  double gemw_opt3[] = { 40.0*cm, 50.0*cm, 50.0*cm};
  double gemh_opt3[] = { 150.0*cm, 200.0*cm, 200.0*cm};

  gemz = new double[ngem];
  gemw = new double[ngem];
  gemh = new double[ngem];

  for( i = 0; i < ngem; i++ ){
      if( fGEMOption == 1 ){
	  for( i = 0; i < ngem; i++ ){
	      gemz[i] = gemz_opt1[i];
	      gemw[i] = gemw_opt1[i];
	      gemh[i] = gemh_opt1[i];
	  }
      }
      if( fGEMOption == 2 ){
	  for( i = 0; i < ngem; i++ ){
	      gemz[i] = gemz_opt2[i];
	      gemw[i] = gemw_opt2[i];
	      gemh[i] = gemh_opt2[i];
	  }
      }
      if( fGEMOption == 3 ){
	  for( i = 0; i < ngem; i++ ){
	      gemz[i] = gemz_opt3[i];
	      gemw[i] = gemw_opt3[i];
	      gemh[i] = gemh_opt3[i];
	  }
      }
  }

  double gempz[] = {
      120.0*um, 3.0*mm, 120.0*um, 5.0*um, 50.0*um, 3.0*mm, // cover + cathode
      5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM0
      5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM1
      5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM2
      10.0*um, 50.0*um, 180.0*um, 120.0*um, 3.0*mm, 120.0*um // readout + honeycomb
  };

  G4Material *gempm[] = {
      NEMAG10, NOMEX, NEMAG10, Copper, Kapton, GEMgas,
      Copper, Kapton, Copper, GEMgas,
      Copper, Kapton, Copper, GEMgas,
      Copper, Kapton, Copper, GEMgas,
      Copper, Kapton, NEMAG10, NEMAG10, NOMEX, NEMAG10
  };

  int gidx, gpidx; 

  double gemmaxw = 0.0;
  double gemmaxh = 0.0;

  for( gidx = 0; gidx < ngem; gidx++ ){
      if( gemw[gidx] > gemmaxw ) gemmaxw = gemw[gidx];
      if( gemh[gidx] > gemmaxh ) gemmaxh = gemh[gidx];
  }

  double gempzsum = 0.0;
  for( gpidx = 0; gpidx < nplane; gpidx++ ){
      gempzsum += gempz[gpidx];
  }

  char gpname[20][50][3][255];

//  G4Box **gbox = new G4Box* [ngem*nplane];


  G4Box *gpbox;
  G4LogicalVolume *gplog;

  G4String GEMSDname = "G4SBS/GEM";
  G4String GEMcolname = "GEMcol";
  G4SBSGEMSD* GEMSD;

  if( !(GEMSD = (G4SBSGEMSD*) SDman->FindSensitiveDetector(GEMSDname)) ){
      GEMSD = new G4SBSGEMSD( GEMSDname, GEMcolname );
      SDman->AddNewDetector(GEMSD);
  }

  char gemname[10][255], gemboxname[10][255], gemlogname[10][255];

  for(  gidx = 0; gidx < ngem; gidx++ ){
      sprintf(gemboxname[gidx], "gembox_%02d", gidx);
      sprintf(gemlogname[gidx], "gemlog_%02d", gidx);

//      G4Box *gembox = new G4Box("gembox", gemmaxw/2.0, gemmaxh/2.0, gempzsum/2.0 );
      G4Box *gembox = new G4Box(gemboxname[gidx], gemw[gidx]/2.0, gemh[gidx]/2.0, gempzsum/2.0 );
      G4LogicalVolume *gemlog = new G4LogicalVolume(gembox, Air, gemlogname[gidx], 0, 0, 0);

      double zrun = 0.0;
      for( gpidx = 0; gpidx < nplane; gpidx++ ){
	  sprintf(gpname[gidx][gpidx][0], "gemplane_%02d_%03d_box", gidx, gpidx );
	  sprintf(gpname[gidx][gpidx][1], "gemplane_%02d_%03d_log", gidx, gpidx );
	  sprintf(gpname[gidx][gpidx][2], "gemplane_%02d_%03d_phy", gidx, gpidx );

	  zrun += gempz[gpidx]/2.0;
	  gpbox = new G4Box( gpname[gidx][gpidx][0], gemw[gidx]/2.0, gemh[gidx]/2.0, gempz[gpidx]/2.0);
	  gplog = new G4LogicalVolume(gpbox, gempm[gpidx], gpname[gidx][gpidx][1], 0, 0, 0);
	          new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, zrun-gempzsum/2.0), gplog, 
		  gpname[gidx][gpidx][2], gemlog, 0, 0, false );
	  zrun += gempz[gpidx]/2.0;

	  // Assign sensitive volume
	  if( gpidx == 5 ){
	      gplog->SetSensitiveDetector(GEMSD);
	  }
      }

      sprintf(gemname[gidx], "gemphys_%02d", gidx);
      new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, gemz[gidx]+detoffset), gemlog, 
	      gemname[gidx], bbdetLog, true, gidx+1, false );
  }
 
  //--------- 48D48 -------------------------------

  Make48D48(WorldLog, f48D48dist + 1219.2*mm/2 );


  ConstructTarget(WorldLog);

  //--------- Calorimeters -------------------------------
  
  // BigBite Preshower

  double psheight = 27*8.5*cm;
  double pswidth  = 2.0*37.0*cm;
  double psdepth  = 8.5*cm;

  G4Box *bbpsbox = new G4Box("bbpsbox", pswidth/2.0, psheight/2.0, psdepth/2.0 );

  G4LogicalVolume* bbpslog = new G4LogicalVolume(bbpsbox, Air, "bbpslog");
  
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth/2.0), bbpslog,
	      "bbpsphys", bbdetLog, false, 0, false);

  // BigBite Shower
  
  double calheight = 27*8.5*cm;
  double calwidth  = 7*8.5*cm;
  double caldepth  = 37.0*cm;

  G4Box *bbcalbox = new G4Box("bbcalbox", calwidth/2.0, calheight/2.0, caldepth/2.0 );
  G4LogicalVolume* bbcallog = new G4LogicalVolume(bbcalbox, Lead, "bbcallog");
  

  
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth+caldepth/2.0+5.0*cm), bbcallog,
	      "bbcalphys", bbdetLog, false, 0, false);

  G4String BBCalSDname = "G4SBS/BBCal";
  G4String BBCalcolname = "BBCalcol";
  G4SBSCalSD* BBCalSD;

  if( !(BBCalSD = (G4SBSCalSD*) SDman->FindSensitiveDetector(BBCalSDname)) ){
      BBCalSD = new G4SBSCalSD( BBCalSDname, BBCalcolname );
      SDman->AddNewDetector(BBCalSD);
  }

  bbcallog->SetSensitiveDetector(BBCalSD);
  bbcallog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );

  // HCAL 
  double hcalheight = 330.0*cm;
  double hcalwidth  = 165.0*cm;
  double hcaldepth  = 101.0*cm;
  double hcalr = fHCALdist+hcaldepth/2.0;

  G4RotationMatrix *hcalrm = new G4RotationMatrix;
  hcalrm->rotateY(-f48D48ang);

  G4Box *hcalbox = new G4Box("hcalbox", hcalwidth/2.0, hcalheight/2.0, hcaldepth/2.0 );
  G4LogicalVolume* hcallog = new G4LogicalVolume(hcalbox, Lead, "hcallog");
     
  new G4PVPlacement(hcalrm, G4ThreeVector(hcalr*sin(f48D48ang), 0.0, hcalr*cos(f48D48ang) ), hcallog,
	      "hcalphys", WorldLog, false, 0, false);

  G4String HCALSDname = "G4SBS/HCAL";
  G4String HCALcolname = "HCALcol";
  G4SBSCalSD* HCalSD;

  if( !(HCalSD = (G4SBSCalSD*) SDman->FindSensitiveDetector(HCALSDname)) ){
      HCalSD = new G4SBSCalSD( HCALSDname, HCALcolname );
      SDman->AddNewDetector(HCalSD);
  }

  SDman->AddNewDetector(HCalSD);
  hcallog->SetSensitiveDetector(HCalSD);
  hcallog->SetUserLimits(  new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );

  //--------- BigBite Cerenkov ------------------------------

//  double cer_mirrorthick = 0.635*mm;
  double cer_mirrorthick = 3.00*mm;
  //double cer_winthick_in   = 1.0*mm;
  double cer_winthick_in   = 0.1*mm;
  double cer_winthick_out  = 0.2*mm;

  double cer_width  =  50.0*cm;
  double cer_height = 200.0*cm;

//  G4Box *cer_winbox = new G4Box("cer_winbox", cer_width/2.0, cer_height/2.0, cer_winthick/2.0 );
  G4Box *cer_winbox_in = new G4Box("cer_winbox_in", cer_width/2.0, cer_height/2.0, cer_winthick_in/2.0 );
  G4Box *cer_winbox_out = new G4Box("cer_winbox_out", cer_width/2.0, cer_height/2.0, cer_winthick_out/2.0 );
  G4Box *cer_mirbox = new G4Box("cer_mirbox", cer_width/2.0, cer_height/2.0, cer_mirrorthick/2.0 );
  G4Box *cer_gasbox = new G4Box("cer_gasbox", cer_width/2.0, cer_height/2.0, fCerDepth/2.0 );

  G4LogicalVolume* cer_winlog_in = new G4LogicalVolume(cer_winbox_in, Aluminum, "cer_winlog_in");
  G4LogicalVolume* cer_winlog_out = new G4LogicalVolume(cer_winbox_out, Aluminum, "cer_winlog_out");
//  G4LogicalVolume* cer_mirlog = new G4LogicalVolume(cer_mirbox, SiO2, "cer_mirlog");
  G4LogicalVolume* cer_mirlog = new G4LogicalVolume(cer_mirbox, Acrylic, "cer_mirlog");
  G4LogicalVolume* cer_gaslog = new G4LogicalVolume(cer_gasbox, C4F8O, "cer_gaslog");

  double thisz = detoffset+fCerDist;

  /*
   *  FIXME  - NO CERENKOV */
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, thisz + cer_winthick_in/2.0 ), cer_winlog_in, "cerwin1", bbdetLog, false, 0, false);
  thisz += cer_winthick_in;
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, thisz + fCerDepth/2.0 ), cer_gaslog, "cergas", bbdetLog, false, 0, false);
  thisz += fCerDepth;
 // new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, thisz + cer_mirrorthick/2.0 ), cer_mirlog, "cermir", bbdetLog, false, 0, false);
  if( fCerDepth > 20.0*cm ){
      new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, fCerDepth/2.0-20.0*cm ), cer_mirlog, "cermir", cer_gaslog, false, 0, false);
  } else {
      new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, fCerDepth/2.0 ), cer_mirlog, "cermir", cer_gaslog, false, 0, false);
  }
//  thisz += cer_mirrorthick;
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, thisz + cer_winthick_out/2.0 ), cer_winlog_out, "cerwin2", bbdetLog, false, 0, false);


  ConstructBeamline(WorldLog);

   /* */
 
  //--------- Reference boxes -------------------------------

  /*
  G4Box *upbox = new G4Box("upbox", 5*cm, 5*cm, 5*cm);
  G4LogicalVolume *upLog=new G4LogicalVolume(upbox, Air,
						  "upLog", 0, 0, 0);
  G4PVPlacement *upPhys=new G4PVPlacement(0,G4ThreeVector(0.0, 4.0*m, 0.0*m),
	  				    upLog,
					       "upPhys",
					       WorldLog,
					       0,false,0);


  G4Box *dbox = new G4Box("dbox", 5*cm, 5*cm, 5*cm);
  G4LogicalVolume *dLog=new G4LogicalVolume(dbox, Air,
						  "dLog", 0, 0, 0);
  G4PVPlacement *dPhys=new G4PVPlacement(0,G4ThreeVector(0.0, 0.0*m, 4.0*m),
	  				    dLog,
					       "dPhys",
					       WorldLog,
					       0,false,0);

  G4Box *xbox = new G4Box("xbox", 5*cm, 5*cm, 5*cm);
  G4LogicalVolume *xLog=new G4LogicalVolume(xbox, Air,
						  "xLog", 0, 0, 0);
  G4PVPlacement *xPhys=new G4PVPlacement(0,G4ThreeVector(4.0*m, 0.0*m, 0.0*m),
	  				    xLog,
					       "xPhys",
					       WorldLog,
					       0,false,0);
  
					       */


  //--------- Visualization attributes -------------------------------
  WorldLog->SetVisAttributes(G4VisAttributes::Invisible);
  bbdetLog->SetVisAttributes(G4VisAttributes::Invisible);
  bbfieldLog->SetVisAttributes(G4VisAttributes::Invisible);
  bbmotherLog->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes * yokeVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
//  yokeVisAtt->SetForceWireframe(true);
  bbyokewgapLog->SetVisAttributes(yokeVisAtt);

  G4VisAttributes * dVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  dVisAtt->SetForceWireframe(true);
  //dLog->SetVisAttributes(dVisAtt);

  G4VisAttributes * xVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  xVisAtt->SetForceWireframe(true);
  //xLog->SetVisAttributes(xVisAtt);

  G4VisAttributes * alVisAtt
    = new G4VisAttributes(G4Colour(0.1,0.1,0.1));
  cer_winlog_in->SetVisAttributes(alVisAtt);
  cer_winlog_out->SetVisAttributes(alVisAtt);

  G4VisAttributes * gasVisAtt
    = new G4VisAttributes(G4Colour(0.6,0.6,1.0));
  gasVisAtt->SetForceWireframe(true);
  cer_gaslog->SetVisAttributes(gasVisAtt);

  G4VisAttributes * psVisAtt
    = new G4VisAttributes(G4Colour(0.3,0.9,0.3));
  psVisAtt->SetForceWireframe(true);
  bbpslog->SetVisAttributes(psVisAtt);

  G4VisAttributes * bbcalVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.6,0.0));
  bbcallog->SetVisAttributes(bbcalVisAtt);
  hcallog->SetVisAttributes(bbcalVisAtt);


  /*
  hcallog->SetVisAttributes(G4VisAttributes::Invisible);
 big48d48Log->SetVisAttributes(G4VisAttributes::Invisible);
 */


  //bigfieldLog->SetVisAttributes(dVisAtt);
  //bbfieldLog->SetVisAttributes(dVisAtt);
  /*
  
  G4VisAttributes * calorimeterBoxVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  calorimeterBoxVisAtt->SetForceWireframe(true);
  calorimeterLog->SetVisAttributes(calorimeterBoxVisAtt);
  
  G4VisAttributes * crystalVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  crystalVisAtt->SetForceWireframe(true);
  theCrystalLog->SetVisAttributes(crystalVisAtt);

  G4VisAttributes * hadCaloBoxVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,1.0));
  hadCaloBoxVisAtt->SetForceWireframe(true);
  hadCaloLog->SetVisAttributes(hadCaloBoxVisAtt);
  
  G4VisAttributes * towerVisAtt
    = new G4VisAttributes(G4Colour(0.5,0.0,1.0));
  towerVisAtt->SetForceWireframe(true);
  theTowerLog->SetVisAttributes(towerVisAtt);
  
  //------------------------------------------------------------------
  
  */
  
  //-----------------------
  // Returns the pointer to
  // the physical world:
  //-----------------------
  return WorldPhys;
}


G4VPhysicalVolume* G4SBSDetectorConstruction::ConstructAllGEp()
{
  G4cout << "\nG4SBSDetectorConstruction....\n" << G4endl;
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  //--------- Material definition ---------
  
  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;
  
  a = 126.9*g/mole;
  G4Element* elI = new G4Element(name="Iodine",   symbol="I",  iz=53., a);
  a = 132.9*g/mole;
  G4Element* elCs= new G4Element(name="Cesium",   symbol="Cs", iz=55., a);

  G4Element *elH = new G4Element("Hydrogen", "H", 1, 1.007*g/mole );
  G4Element *elC = new G4Element("Carbon", "C", 6, 12.011*g/mole );
  G4Element *elN = new G4Element("Nitrogen", "N", 7, 14.007*g/mole );
  G4Element *elO = new G4Element("Oxygen", "O", 8, 16.000*g/mole );
  G4Element *elF = new G4Element("Fluorine", "F", 9, 18.998*g/mole );
  G4Element *elAl = new G4Element("Aluminum", "Al", 13, 26.982*g/mole );
  G4Element *elSi = new G4Element("Silicon", "Si", 14, 28.086*g/mole );
  G4Element *elCa = new G4Element("Calcium", "Ca", 20, 40.078*g/mole );
  G4Element *elSr = new G4Element("Strontium", "Sr", 38, 87.62*g/mole );
  G4Element *elBa = new G4Element("Barium", "Ba", 56, 137.327*g/mole );

   G4Element* elCl  = new G4Element("Chlorine",  "Cl", z=17, a=   35.453*g/mole);
     G4Element* elAr  = new G4Element("Argon",     "Ar", z=18, a=    39.95*g/mole);

  G4Material* Lead  = new G4Material(name="Lead", z=82., a=208.0*g/mole, density=11.34*g/cm3);

  density = 4.51*g/cm3;
  G4Material* CsI = new G4Material(name="CsI", density, nel = 2);
  CsI->AddElement(elI, .5);
  CsI->AddElement(elCs,.5);
  a = 4.0*g/mole;
  density = 0.1786e-03*g/cm3;
  
  /*
  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Fer", z=26., a, density);
  */
  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  double bigden = 1e9*g/cm3;

  // Cell Glass - GE180 Aluminosilicate Glass
  // SiO2 60.3%
  G4Material* SiO2 = new G4Material("SiO2", 2.2*g/cm3, 2 );
  SiO2->AddElement(elSi, 1);
  SiO2->AddElement(elO, 2);
  // BaO  18.2%
  G4Material* BaO = new G4Material("BaO", bigden, 2 );
  BaO->AddElement(elBa, 1);
  BaO->AddElement(elO, 1);
  // Al2O3 14.3%
  G4Material* Al2O3 = new G4Material("Al2O3", bigden, 2 );
  Al2O3->AddElement(elAl, 2);
  Al2O3->AddElement(elO, 3);
  // CaO   6.5%
  G4Material* CaO = new G4Material("CaO", bigden, 2 );
  CaO->AddElement(elCa, 1);
  CaO->AddElement(elO, 1);
  // SrO   0.25%
  G4Material* SrO = new G4Material("SrO", bigden, 2 );
  SrO->AddElement(elSr, 1);
  SrO->AddElement(elO, 1);
  //
  density = 1.19*g/cm3;
  G4Material* Acrylic = new G4Material(name="Acrylic", density, nel=3);
  Acrylic->AddElement(elC, 5);
  Acrylic->AddElement(elH, 8);
  Acrylic->AddElement(elO, 2);
  // Density 2.76 g/cm^3
  // Index of Refraction 1.536
  G4Material* GE180 = new G4Material("GE180", 2.76*g/cm3, 5);
  GE180->AddMaterial(SiO2, 0.6039);
  GE180->AddMaterial(BaO, 0.1829);
  GE180->AddMaterial(Al2O3, 0.1439);
  GE180->AddMaterial(CaO, 0.0659);
  GE180->AddMaterial(SrO, 0.0034);

  //--------- GEM Materials  ---------
  // (stolen from GEMC)

  G4Material* NOMEX_pure = new G4Material("NOMEX_pure", density = 1.38*g/cm3, 5);
  NOMEX_pure -> AddElement(elH,0.04);
  NOMEX_pure -> AddElement(elC,0.54);
  NOMEX_pure -> AddElement(elN,0.09);
  NOMEX_pure -> AddElement(elO,0.10);
  NOMEX_pure -> AddElement(elCl,0.23);

  G4Material* NOMEX = new G4Material("NOMEX",density = 0.04*g/cm3, 2);
  NOMEX -> AddMaterial(NOMEX_pure,0.45);
  NOMEX -> AddMaterial(Air,0.55);

  G4Material* NEMAG10 = new G4Material("NEMAG10", 1.70*g/cm3, nel=4);
  NEMAG10 -> AddElement(elSi, 1);
  NEMAG10 -> AddElement(elO , 2);
  NEMAG10 -> AddElement(elC , 3);
  NEMAG10 -> AddElement(elH , 3);

  G4double density_Ar = 1.7823*mg/cm3 ;
  G4Material* Argon = new G4Material("Argon"  , density_Ar, nel=1);
  Argon->AddElement(elAr, 1);

  G4double density_CO2 = 1.977*mg/cm3;
  G4Material* CO2 = new G4Material("CO2", density_CO2, nel=2);
  CO2->AddElement(elC, 1);
  CO2->AddElement(elO, 2);

  G4double density_CH2 = 0.95*g/cm3;
  G4Material* CH2 = new G4Material("CH2", density_CH2, nel=2);
  CH2->AddElement(elC, 1);
  CH2->AddElement(elH, 2);

  G4double density_CH = 0.95*g/cm3;
  G4Material* CH = new G4Material("CH", density_CH, nel=2);
  CH->AddElement(elC, 1);
  CH->AddElement(elH, 1);

  // 1.5 Atmosphere C4F8O for cerkenkov
  G4double density_C4F8O = 9.64*mg/cm3; // density at 1ATM
  G4Material* C4F8O = new G4Material("C4F8O", density_C4F8O*1.5, nel=3);
  C4F8O->AddElement(elC, 4);
  C4F8O->AddElement(elF, 8);
  C4F8O->AddElement(elO, 1);

  G4double density_ArCO2 = .7*density_Ar + .3*density_CO2;
  // Use ArCO2
  G4Material *GEMgas= new G4Material("GEMgas", density_ArCO2, nel=2);
  GEMgas->AddMaterial(Argon, 0.7*density_Ar/density_ArCO2) ;
  GEMgas->AddMaterial(CO2, 0.3*density_CO2/density_ArCO2) ;

  G4Material *Copper= new G4Material("Copper", z=29, a=   63.55*g/mole, density = 8.96*g/cm3);

  G4Material *Kapton = new G4Material("Kapton",   density = 1.42*g/cm3, nel=4);
  Kapton->AddElement(elH, 0.026362);
  Kapton->AddElement(elC, 0.691133);
  Kapton->AddElement(elN, 0.073270);
  Kapton->AddElement(elO, 0.209235);



  //--------- G4VSolid, G4LogicalVolume, G4VPhysicalVolume  ---------
  
  //--------------
  // World:
  //--------------
  G4Box *WorldBox= new G4Box("WorldBox",20*m, 20*m, 20*m);
  G4LogicalVolume *WorldLog=new G4LogicalVolume(WorldBox,Air,
						  "WorldLogical", 0, 0, 0);
  G4PVPlacement *WorldPhys=new G4PVPlacement(0,G4ThreeVector(),
					       "WorldPhysical",
					       WorldLog,
					       0,false,0);

  //-----------------------------
  //  BigBite

  printf("BigCal at %f deg\n", fBBang/deg);

  G4RotationMatrix *bbrm = new G4RotationMatrix;
  bbrm->rotateY(fBBang);

  // Ecal will act as BBcal detector

  double bigcalheight = (24*4.5+32*4.0)*cm;
  double bigcalwidth  = 44.10*2.54*cm;
  double bigcaldepth  = 15.75*2.54*cm;
  double bbr = fBBdist+bigcaldepth/2.0;

  double CH2depth = 15.0*cm;
  double CHdepth  = 6.0*cm;

  G4Box *bigcalbox = new G4Box("bigcalbox", bigcalwidth/2.0, bigcalheight/2.0, bigcaldepth/2.0 );
  G4LogicalVolume* bigcallog = new G4LogicalVolume(bigcalbox, Lead, "bigcallog");

  G4Box *CH2box = new G4Box("ch2box", bigcalwidth/2.0, bigcalheight/2.0, CH2depth/2.0 );
  G4LogicalVolume* ch2boxlog = new G4LogicalVolume(CH2box, CH2, "ch2log");
  G4Box *CHbox = new G4Box("chbox", bigcalwidth/2.0, bigcalheight/2.0, CHdepth/2.0 );
  G4LogicalVolume* chboxlog = new G4LogicalVolume(CHbox, CH, "chlog");

  double ch2r = bbr - bigcaldepth/2.0 - CHdepth - CH2depth/2.0;
  double chr = bbr - bigcaldepth/2.0 - CHdepth/2.0;

  new G4PVPlacement(bbrm, G4ThreeVector(bbr*sin(-fBBang), 0.0, bbr*cos(-fBBang) ), bigcallog,
	      "bigcalphys", WorldLog, false, 0, false);
  new G4PVPlacement(bbrm, G4ThreeVector(ch2r*sin(-fBBang), 0.0, ch2r*cos(-fBBang) ), ch2boxlog,
	      "ch2boxphys", WorldLog, false, 0, false);
  new G4PVPlacement(bbrm, G4ThreeVector(chr*sin(-fBBang), 0.0, chr*cos(-fBBang) ), chboxlog,
	      "chboxphys", WorldLog, false, 0, false);

  G4String BBCALSDname = "G4SBS/BBCal";
  G4String BBCALcolname = "BBCalcol";
  G4SBSCalSD* BBCalSD;

  if( !(BBCalSD = (G4SBSCalSD*) SDman->FindSensitiveDetector(BBCALSDname)) ){
      BBCalSD = new G4SBSCalSD( BBCALSDname, BBCALcolname );
      SDman->AddNewDetector(BBCalSD);
  }

  SDman->AddNewDetector(BBCalSD);
  bigcallog->SetSensitiveDetector(BBCalSD);
  bigcallog->SetUserLimits(  new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );

  /////////////////////////////////////////////////////////
  G4RotationMatrix *hcalrm = new G4RotationMatrix;
  hcalrm->rotateY(-f48D48ang);

  double sbsboxpitch = 5.0*deg;

  G4RotationMatrix *gemrm = new G4RotationMatrix;
  gemrm->rotateY(-f48D48ang);
  gemrm->rotateX( sbsboxpitch);

  // SBS box

  double sbsdepth  = 3.0*m;
  double sbswidth  = 2.0*m;
  double sbsheight = 2.1*m;

  double sbsr = fHCALdist-4.106*m + sbsheight*sin(sbsboxpitch)/2+sbsdepth/2;

  G4Box *sbsbox = new G4Box("sbsbox", sbswidth/2.0, sbsheight/2.0, sbsdepth/2.0 );
  G4LogicalVolume* sbslog = new G4LogicalVolume(sbsbox, Air, "sbslog");
  new G4PVPlacement(gemrm, G4ThreeVector(sbsr*sin(f48D48ang), (sbsr-f48D48dist)*sin(sbsboxpitch), sbsr*cos(f48D48ang) ), sbslog,
	      "sbsphys", WorldLog, false, 0, false);

  //  6 GEMs in the front tracker

  double detoffset = 0.05*m - sbsdepth/2.0;

  const int nplane = 24;

  int i;
  int ngem = 0;
  //double gemdsep;
  double *gemz;
  double *gemw;
  double *gemh;

  ngem = 14;
//  gemdsep = 0.05*m;

  //double pairspac  = 0.30*m;
  //
  // GEM option 1
  
  gemz = new double[ngem];
  gemw = new double[ngem];
  gemh = new double[ngem];

  for( i = 0; i < ngem; i++ ){
      if( i < 6 ){
	  gemz[i] = ((double) i)*9*cm;
	  gemw[i] = 40.0*cm;
	  gemh[i] = 150.0*cm;
      } else if( i < 10 ) {
	  gemz[i] = ((double) i-6)*16*cm + 1.2*m;
//	  gemz[i] = pairspac*((i-6)/2) + (i%2)*gemdsep + 1.2*m;
	  gemw[i] = 50.0*cm;
	  gemh[i] = 200.0*cm;
      } else {
	  gemz[i] = ((double) i-10)*16.*cm + 2.316*m;
//	  gemz[i] = pairspac*((i-10)/2) + (i%2)*gemdsep + 2.316*m;
	  gemw[i] = 50.0*cm;
	  gemh[i] = 200.0*cm;
      }

      //printf("i = %d  z = %f\n", i, gemz[i]/m);
  }


  double gempz[] = {
      120.0*um, 3.0*mm, 120.0*um, 5.0*um, 50.0*um, 3.0*mm, // cover + cathode
      5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM0
      5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM1
      5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM2
      10.0*um, 50.0*um, 180.0*um, 120.0*um, 3.0*mm, 120.0*um // readout + honeycomb
  };

  G4Material *gempm[] = {
      NEMAG10, NOMEX, NEMAG10, Copper, Kapton, GEMgas,
      Copper, Kapton, Copper, GEMgas,
      Copper, Kapton, Copper, GEMgas,
      Copper, Kapton, Copper, GEMgas,
      Copper, Kapton, NEMAG10, NEMAG10, NOMEX, NEMAG10
  };

  int gidx, gpidx; 

  double gemmaxw = 0.0;
  double gemmaxh = 0.0;

  for( gidx = 0; gidx < ngem; gidx++ ){
      if( gemw[gidx] > gemmaxw ) gemmaxw = gemw[gidx];
      if( gemh[gidx] > gemmaxh ) gemmaxh = gemh[gidx];
  }

  double gempzsum = 0.0;
  for( gpidx = 0; gpidx < nplane; gpidx++ ){
      gempzsum += gempz[gpidx];
  }

  char gpname[20][50][3][255];

//  G4Box **gbox = new G4Box* [ngem*nplane];


  G4Box *gpbox;
  G4LogicalVolume *gplog;

  G4String GEMSDname = "G4SBS/GEM";
  G4String GEMcolname = "GEMcol";
  G4SBSGEMSD* GEMSD;

  if( !(GEMSD = (G4SBSGEMSD*) SDman->FindSensitiveDetector(GEMSDname)) ){
      GEMSD = new G4SBSGEMSD( GEMSDname, GEMcolname );
      SDman->AddNewDetector(GEMSD);
  }

  char gemname[100][255], gemboxname[100][255], gemlogname[100][255];

  for(  gidx = 0; gidx < ngem; gidx++ ){
      sprintf(gemboxname[gidx], "gembox_%02d", gidx);
      sprintf(gemlogname[gidx], "gemlog_%02d", gidx);

//      G4Box *gembox = new G4Box("gembox", gemmaxw/2.0, gemmaxh/2.0, gempzsum/2.0 );
      G4Box *gembox = new G4Box(gemboxname[gidx], gemw[gidx]/2.0, gemh[gidx]/2.0, gempzsum/2.0 );
      G4LogicalVolume *gemlog = new G4LogicalVolume(gembox, Air, gemlogname[gidx], 0, 0, 0);

      double zrun = 0.0;
      for( gpidx = 0; gpidx < nplane; gpidx++ ){
	  sprintf(gpname[gidx][gpidx][0], "gemplane_%02d_%03d_box", gidx, gpidx );
	  sprintf(gpname[gidx][gpidx][1], "gemplane_%02d_%03d_log", gidx, gpidx );
	  sprintf(gpname[gidx][gpidx][2], "gemplane_%02d_%03d_phy", gidx, gpidx );

	  zrun += gempz[gpidx]/2.0;
	  gpbox = new G4Box( gpname[gidx][gpidx][0], gemw[gidx]/2.0, gemh[gidx]/2.0, gempz[gpidx]/2.0);
	  gplog = new G4LogicalVolume(gpbox, gempm[gpidx], gpname[gidx][gpidx][1], 0, 0, 0);
	          new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, zrun-gempzsum/2.0), gplog, 
		  gpname[gidx][gpidx][2], gemlog, 0, 0, false );
	  zrun += gempz[gpidx]/2.0;

	  // Assign sensitive volume
//	  if( gpidx == 5 && gidx < 6 ){ // just make first six tracking for now
	  if( gpidx == 5 ){
	      gplog->SetSensitiveDetector(GEMSD);
	  }
      }

      sprintf(gemname[gidx], "gemphys_%02d", gidx);
      new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, gemz[gidx]+detoffset), gemlog, 
	      gemname[gidx], sbslog, true, gidx+1, false );
  }
 
  ////   CH2 blocks

  double anaheight = 200.0*cm;
  double anawidth  = 44.0*2.54*cm;
  double anadepth  = 22.0*2.54*cm;
  
  G4Box *anabox = new G4Box("anabox", anawidth/2.0, anaheight/2.0, anadepth/2.0 );
  G4LogicalVolume* analog = new G4LogicalVolume(anabox, CH2, "analog");

     new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset +  58.53*cm + anadepth/2.0 ), analog,
	      "anaphys1", sbslog, false, 0, false);
     new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset + 170.3*cm + anadepth/2.0 ), analog,
	      "anaphys1", sbslog, false, 0, false);

  //  HCAL
  
  double hcalheight = 330.0*cm;
  double hcalwidth  = 165.0*cm;
  double hcaldepth  = 101.0*cm;

  G4Box *hcalbox = new G4Box("hcalbox", hcalwidth/2.0, hcalheight/2.0, hcaldepth/2.0 );
  G4LogicalVolume* hcallog = new G4LogicalVolume(hcalbox, Lead, "hcallog");

  /*
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset + 2.97*m + hcaldepth/2.0 ), hcallog,
	      "hcalphys", sbslog, false, 0, false);
	      */

  double HCALvertoffset = 49.7*cm;
  new G4PVPlacement(hcalrm, G4ThreeVector( (fHCALdist + hcaldepth/2.0)*sin(f48D48ang), HCALvertoffset, (fHCALdist + hcaldepth/2.0)*cos(f48D48ang) ), hcallog,
	      "hcalphys", WorldLog, false, 0, false);

  G4String HCALSDname = "G4SBS/HCAL";
  G4String HCALcolname = "HCALcol";
  G4SBSCalSD* HCalSD;

  if( !(HCalSD = (G4SBSCalSD*) SDman->FindSensitiveDetector(HCALSDname)) ){
      HCalSD = new G4SBSCalSD( HCALSDname, HCALcolname );
      SDman->AddNewDetector(HCalSD);
  }

  SDman->AddNewDetector(HCalSD);
  hcallog->SetSensitiveDetector(HCalSD);
  hcallog->SetUserLimits(  new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );



  //--------- 48D48 -------------------------------
  Make48D48(WorldLog, f48D48dist + 1219.2*mm/2 );

  ConstructTarget(WorldLog);
  ConstructBeamline(WorldLog);


  //--------- Visualization attributes -------------------------------
  WorldLog->SetVisAttributes(G4VisAttributes::Invisible);
  sbslog->SetVisAttributes(G4VisAttributes::Invisible);


  G4VisAttributes * dVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  dVisAtt->SetForceWireframe(true);

  G4VisAttributes * xVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  xVisAtt->SetForceWireframe(true);


  G4VisAttributes * bbcalVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.6,0.0));
  bigcallog->SetVisAttributes(bbcalVisAtt);
  hcallog->SetVisAttributes(bbcalVisAtt);

  G4VisAttributes * ch2VisAtt
    = new G4VisAttributes(G4Colour(1.0,0.7,0.0));
  G4VisAttributes * chVisAtt
    = new G4VisAttributes(G4Colour(1.0,1.0,0.0));

  ch2boxlog->SetVisAttributes(ch2VisAtt);
  chboxlog->SetVisAttributes(chVisAtt);

  G4VisAttributes * anaVisAtt
    = new G4VisAttributes(G4Colour(0.6,0.6,0.0));
  analog->SetVisAttributes(anaVisAtt);


  return WorldPhys;
}

void G4SBSDetectorConstruction::ConstructTarget( G4LogicalVolume *worldlog ){

  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;
  
  a = 126.9*g/mole;
  G4Element* elI = new G4Element(name="Iodine",   symbol="I",  iz=53., a);
  a = 132.9*g/mole;
  G4Element* elCs= new G4Element(name="Cesium",   symbol="Cs", iz=55., a);

  G4Element *elH = new G4Element("Hydrogen", "H", 1, 1.007*g/mole );
  G4Element *elD = new G4Element("Deuterium", "D", 1, 2.014*g/mole );
  G4Element *el3He = new G4Element("Helium3", "3He", 2, 3.016*g/mole );
  G4Element *elC = new G4Element("Carbon", "C", 6, 12.011*g/mole );
  G4Element *elN = new G4Element("Nitrogen", "N", 7, 14.007*g/mole );
  G4Element *elO = new G4Element("Oxygen", "O", 8, 16.000*g/mole );
  G4Element *elAl = new G4Element("Aluminum", "Al", 13, 26.982*g/mole );
  G4Element *elSi = new G4Element("Silicon", "Si", 14, 28.086*g/mole );
  G4Element *elCa = new G4Element("Calcium", "Ca", 20, 40.078*g/mole );
  G4Element *elSr = new G4Element("Strontium", "Sr", 38, 87.62*g/mole );
  G4Element *elBa = new G4Element("Barium", "Ba", 56, 137.327*g/mole );

  G4Material* Vacuum = new G4Material(name="vacuum", z=1., a=1.0*g/mole, density=1e-9*g/cm3);

  G4Material* Aluminum = new G4Material(name="Aluminum", z=13., a=26.98*g/mole, density=2.7*g/cm3);


  density = 4.51*g/cm3;
  G4Material* CsI = new G4Material(name="CsI", density, nel = 2);
  CsI->AddElement(elI, .5);
  CsI->AddElement(elCs,.5);
  
  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  double bigden = 1e9*g/cm3;

  // Cell Glass - GE180 Aluminosilicate Glass
  // SiO2 60.3%
  G4Material* SiO2 = new G4Material("SiO2", 2.2*g/cm3, 2 );
  SiO2->AddElement(elSi, 1);
  SiO2->AddElement(elO, 2);
  // BaO  18.2%
  G4Material* BaO = new G4Material("BaO", bigden, 2 );
  BaO->AddElement(elBa, 1);
  BaO->AddElement(elO, 1);
  // Al2O3 14.3%
  G4Material* Al2O3 = new G4Material("Al2O3", bigden, 2 );
  Al2O3->AddElement(elAl, 2);
  Al2O3->AddElement(elO, 3);
  // CaO   6.5%
  G4Material* CaO = new G4Material("CaO", bigden, 2 );
  CaO->AddElement(elCa, 1);
  CaO->AddElement(elO, 1);
  // SrO   0.25%
  G4Material* SrO = new G4Material("SrO", bigden, 2 );
  SrO->AddElement(elSr, 1);
  SrO->AddElement(elO, 1);
  //
  density = 1.19*g/cm3;
  G4Material* Acrylic = new G4Material(name="Acrylic", density, nel=3);
  Acrylic->AddElement(elC, 5);
  Acrylic->AddElement(elH, 8);
  Acrylic->AddElement(elO, 2);
  // Density 2.76 g/cm^3
  // Index of Refraction 1.536
  G4Material* GE180 = new G4Material("GE180", 2.76*g/cm3, 5);
  GE180->AddMaterial(SiO2, 0.6039);
  GE180->AddMaterial(BaO, 0.1829);
  GE180->AddMaterial(Al2O3, 0.1439);
  GE180->AddMaterial(CaO, 0.0659);
  GE180->AddMaterial(SrO, 0.0034);


  double gasden = fTargDen;
//  gasden = 10.5*atmosphere/(300*kelvin*k_Boltzmann);
  G4Material *refH2 = new G4Material("refH2", gasden, 1 );
  refH2->AddElement(elH, 1);

 // gasden = 10.5*atmosphere/(300*kelvin*k_Boltzmann);
  G4Material *refN2 = new G4Material("refN2", gasden, 1 );
  refN2->AddElement(elN, 1);

 // gasden = 10.77*atmosphere/(300*kelvin*k_Boltzmann);
  G4Material *pol3He = new G4Material("pol3He", gasden, 1 );
  pol3He->AddElement(el3He, 1);

  double LH2den = 0.071*g/cm3;
  G4Material *LH2mat = new G4Material("LH2", LH2den, 1 );
  LH2mat->AddElement(elH, 1);

  double LD2den = 162.4*kg/m3;
  G4Material *LD2mat = new G4Material("LD2", LD2den, 1 );
  LD2mat->AddElement(elD, 1);


  //--------- Glass target cell -------------------------------

  double wallthick = 1.61*mm;
  double capthick  = 0.126*mm;
  double cellradius    = 0.75*2.54*cm/2.0;

  G4Tubs *targ_tube = new G4Tubs("targ_tube", cellradius-wallthick, cellradius, fTargLen/2.0, 0.*deg, 360.*deg );
  G4Tubs *targ_cap = new G4Tubs("targ_cap", 0.0, cellradius, capthick/2.0, 0.*deg, 360.*deg );

  G4LogicalVolume* targ_tube_log = new G4LogicalVolume(targ_tube, GE180,"targ_tube_log");
  G4LogicalVolume* targ_cap_log = new G4LogicalVolume(targ_cap, GE180,"targ_cap_log");

  /* FIXME
   * */
  if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
      new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), targ_tube_log,
	      "targ_tube_phys", worldlog, false, 0);

      new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, fTargLen/2.0+capthick/2.0), targ_cap_log,
	      "targ_cap_phys1", worldlog, false, 0);
      new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -fTargLen/2.0-capthick/2.0), targ_cap_log,
	      "targ_cap_phys2", worldlog, false, 0);
  }

  
  /**/


  // gas
  G4Tubs *gas_tube = new G4Tubs("gas_tube", 0.0, cellradius-wallthick,fTargLen/2.0, 0.*deg, 360.*deg );
  G4LogicalVolume* gas_tube_log;


  if( fTargType == kH2 || fTargType == kNeutTarg ){
      gas_tube_log = new G4LogicalVolume(gas_tube, refH2, "gas_tube_log");
  }
  if( fTargType == k3He ){
      gas_tube_log = new G4LogicalVolume(gas_tube, pol3He, "gas_tube_log");
  }

  /*
   * FIXME*/
  if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
      new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), gas_tube_log,
	      "gas_tube_phys", worldlog, false, 0);
  }
  /**/
  
  //--------- Cryo target cell -------------------------------

  wallthick   = 635*um;
  double upcapthick    = 100*um;
  double downcapthick  = 125*um;
  double targconeang = 15.*deg;
  double cellupradius = 2.0*cm;

  double celldownradius = fTargLen*sin(targconeang);
  double cellconelen = fTargLen*cos(targconeang);
  double cellconeang = atan((celldownradius-cellupradius)/cellconelen);

  // Aluminum shell sphere
  G4Sphere *shellsph = new G4Sphere("shellsph", 0, fTargLen, 0, 360.0*deg, 0, targconeang);
  // Aluminum shell cone
  G4Cons *shellcon = new G4Cons("shellcon", 0.0, cellupradius, 0.0, celldownradius,  cellconelen/2, 0.0, 360.0*deg);
  // Union 
  G4UnionSolid *cryoshell = new G4UnionSolid("cryoshell", shellcon, shellsph, 0, G4ThreeVector(0,0,-cellconelen/2));

  double cryoupradius =  cellupradius - wallthick/cos(cellconeang);
  double cryoconelen = cellconelen - upcapthick - downcapthick; 

  double cryodownradius = cryoupradius + cryoconelen*tan(cellconeang);

  double cryooffset = cryoconelen/2+upcapthick- cellconelen/2;

  double cryoang = 14.95*deg;
  
  // Cryo sphere
  G4Sphere *cryosph = new G4Sphere("cryosph", 0, fTargLen-downcapthick, 0*deg, 360.0*deg, 0, cryoang);
  // Cryo cone
  G4Cons *cryocon = new G4Cons("shellcon", 0.0, cryoupradius, 0.0, cryodownradius,  cryoconelen/2,  0.0*deg, 360.0*deg);
  // Union 
  G4UnionSolid *cryovol1 = new G4UnionSolid("cryovol1", cryocon, cryosph, 0, G4ThreeVector(0.0,0.0,-cryoconelen/2.-upcapthick ));

  double trimboxsize = 50.0*cm;
  G4Box *cryotrimbox = new G4Box("cryotrimbox", trimboxsize, trimboxsize, trimboxsize);

  G4SubtractionSolid *cryovol = new G4SubtractionSolid("cryovol", cryovol1, cryotrimbox, 0, G4ThreeVector(0.0,0.0,-trimboxsize-cryoconelen/2 ));

  /*
  targ_tube = new G4Tubs("targ_tube", cellradius-wallthick, cellradius, fTargLen/2.0, 0.*deg, 360.*deg );
  G4Tubs *targ_ucap = new G4Tubs("targ_ucap", 0.0, cellradius, upcapthick/2.0, 0.*deg, 360.*deg );
  G4Tubs *targ_dcap = new G4Tubs("targ_dcap", 0.0, cellradius, downcapthick/2.0, 0.*deg, 360.*deg );
  targ_tube_log = new G4LogicalVolume(targ_tube, Aluminum,"targ_tube_log");
  */

  targ_tube_log = new G4LogicalVolume(cryoshell, Aluminum,"targ_tube_log");

  //////////////////////////////////////////////////////////////////

  G4double entpipe_rin = 31.75*mm;
  G4double extpipe_rin = 41.28*mm;

  G4double extpipestart = 2.06*m;
  // 2.06m is where the main exit pipe starts
  G4double extpipe_len;

  double sheight = 1.2*m;

  double swallthick   = 0.38*mm;
  double swallrad     = 1.143*m/2;
  double swallrad_in  = 1.041*m/2;

  double hcal_ang_min = -55*deg;
  double hcal_ang_max = -7*deg;
  double hcal_win_h = 0.4*m;

  double bb_ang_min = 18*deg;
  double bb_ang_max = 80*deg;
  double bb_win_h = 0.5*m;

  if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
      // Gas target
      extpipe_len = extpipestart - 1.0*m;
  } else {
      // Cryotarget
      extpipe_len = extpipestart -  swallrad;
  }

  G4Tubs *swall = new G4Tubs("scham_wall", swallrad_in, swallrad, sheight/2, 0.*deg, 360.*deg );

  // Cut out for windows
  
  G4Tubs *swall_hcalcut = new G4Tubs("scham_wall_hcalcut", swallrad_in-2*cm, swallrad+2*cm, hcal_win_h, hcal_ang_min, hcal_ang_max-hcal_ang_min);
  G4Tubs *swall_bbcut = new G4Tubs("scham_wall_bbcut", swallrad_in-2*cm, swallrad+2*cm, bb_win_h, bb_ang_min, bb_ang_max-bb_ang_min);

  G4SubtractionSolid *swallcut = new G4SubtractionSolid("swallcut1", swall, swall_hcalcut);
  swallcut = new G4SubtractionSolid("swallcut2", swallcut, swall_bbcut);

  G4Tubs *swall_hcalwin = new G4Tubs("scham_wall_hcalwin", swallrad_in, swallrad_in+swallthick, hcal_win_h, hcal_ang_min, hcal_ang_max-hcal_ang_min);
  G4Tubs *swall_bbwin = new G4Tubs("scham_wall_bbwin", swallrad_in, swallrad_in+swallthick, bb_win_h, bb_ang_min, bb_ang_max-bb_ang_min);

  G4Tubs *exttube = new G4Tubs("exitpipetube", extpipe_rin, extpipe_rin+0.120*cm, extpipe_len/2, 0.*deg, 360.*deg );
  G4Tubs *extvactube = new G4Tubs("exitpipetube_vac", 0.0, extpipe_rin, extpipe_len, 0.*deg, 360.*deg );
  G4LogicalVolume *extpipe_log = new G4LogicalVolume(exttube, Aluminum,"extpipe_log");
  G4LogicalVolume *extvac_log = new G4LogicalVolume(extvactube, Vacuum,"extvac_log");


  if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
      // Add in exit Al window
      
      double extwin_thick = 5.0e-4*cm;
      
      // 5cm is exit pipe opening
      G4Tubs *extwin = new G4Tubs("ent_win", 0.0, 5.0*cm, extwin_thick/2, 0.*deg, 360.*deg );
      G4LogicalVolume *ext_winlog = new G4LogicalVolume(extwin, Aluminum, "extwin_log", 0, 0, 0);
      //  place at beginning of pipe
      new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, 162.2*cm - extwin_thick/2), ext_winlog, "extwin_phys", worldlog,false,0);

      ext_winlog->SetVisAttributes(new G4VisAttributes(G4Colour(0.6,0.6,0.6)));
  } else {
      // Scattering chamber parts

      new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, extpipestart-extpipe_len/2), extpipe_log, "extpipe_phys", worldlog, false, 0);
      new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, extpipestart-extpipe_len/2), extvac_log, "extvacpipe_phys", worldlog, false, 0);
  }


  //  Place exit pipe tube

  G4Tubs *swall_enthole = new G4Tubs("scham_wall_enthole", 0.0, entpipe_rin, 20.0*cm, 0.*deg, 360.*deg );
  G4Tubs *swall_exthole = new G4Tubs("scham_wall_exthole", 0.0, extpipe_rin, 20.0*cm, 0.*deg, 360.*deg );

  G4RotationMatrix *chamholerot = new G4RotationMatrix;
  chamholerot->rotateY(90.0*deg);

  //  Cut holes in the scattering chamber
  G4SubtractionSolid* swall_holes = new G4SubtractionSolid("swall_enthole", swallcut, swall_enthole, chamholerot, G4ThreeVector(-(swallrad+swallrad_in)/2, 0.0, 0.0) );
  swall_holes = new G4SubtractionSolid("swall_holes", swall_holes, swall_exthole, chamholerot, G4ThreeVector((swallrad+swallrad_in)/2, 0.0, 0.0) );

  G4LogicalVolume *swall_log = new G4LogicalVolume(swall_holes, Aluminum,"scham_wall_log");

  G4LogicalVolume *sc_hcalwin_log = new G4LogicalVolume(swall_hcalwin, Aluminum,"sc_hcalwin_log");
  G4LogicalVolume *sc_bbwin_log = new G4LogicalVolume(swall_bbwin, Aluminum,"sc_bbwin_log");

  G4RotationMatrix *schamrot = new G4RotationMatrix;
  schamrot->rotateX(-90.0*deg);
  schamrot->rotateZ(-90.0*deg);

  G4RotationMatrix *targrot = new G4RotationMatrix;
  targrot->rotateY(-90.0*deg);

  G4Tubs *chamber_inner = new G4Tubs("chamber_inner", 0.0, swallrad_in,  sheight/2, 0*deg, 360*deg );
  G4LogicalVolume* chamber_inner_log = new G4LogicalVolume(chamber_inner, Vacuum, "cham_inner_log");

  // Top and bottom
  G4Tubs *sc_topbottom = new G4Tubs("scham_topbottom", 0.0, swallrad, (swallrad-swallrad_in)/2, 0.*deg, 360.*deg );
  G4LogicalVolume* sc_topbottom_log = new G4LogicalVolume(sc_topbottom, Aluminum, "scham_topbottom_log");

  /*
   * FIXME*/
  if( fTargType == kLH2 || fTargType == kLD2 ){
      new G4PVPlacement(targrot, G4ThreeVector(0.0, 0.0, 0.0), targ_tube_log,
	      "targ_tube_phys", chamber_inner_log, false, 0);

      /*
      new G4PVPlacement(targrot, G4ThreeVector(fTargLen/2.0+downcapthick/2.0, 0.0, 0.0), targ_dcap_log,
	      "targ_dcap_phys", chamber_inner_log, false, 0);
      new G4PVPlacement(targrot, G4ThreeVector(-fTargLen/2.0-upcapthick/2.0, 0.0, 0.0), targ_ucap_log,
	      "targ_ucap_phys", chamber_inner_log, false, 0);
	      */

      // Scattering chamber
	  new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), swall_log,
		  "scham_wall_phys", worldlog, false, 0);
	  new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), chamber_inner_log,
		  "chamber_inner_phys", worldlog, false, 0);

	  new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), sc_hcalwin_log,
		  "sc_hcalwin_phys", worldlog, false, 0);
	  new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), sc_bbwin_log,
		  "sc_bbwin_phys", worldlog, false, 0);

	  new G4PVPlacement(schamrot, G4ThreeVector(0.0, sheight/2.0 + (swallrad-swallrad_in)/2, 0.0), sc_topbottom_log,
		  "scham_top_phys", worldlog, false, 0);
	  new G4PVPlacement(schamrot, G4ThreeVector(0.0, -sheight/2.0 - (swallrad-swallrad_in)/2, 0.0), sc_topbottom_log,
		  "scham_bot_phys", worldlog, false, 0);
  }
  /**/

//  G4Tubs *cryo_tube = new G4Tubs("cryo_tube", 0.0, cellradius-wallthick,fTargLen/2.0, 0.*deg, 360.*deg );
  G4LogicalVolume* cryo_tube_log = NULL;


  if( fTargType == kLH2 ){
      cryo_tube_log = new G4LogicalVolume(cryovol, LH2mat, "cryo_tube_log");
  }
  if( fTargType == kLD2 ){
      cryo_tube_log = new G4LogicalVolume(cryovol, LD2mat, "cryo_tube_log");
  }

  /*
   * FIXME */
  if( fTargType == kLD2 || fTargType == kLH2 ){
	  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, cryooffset), cryo_tube_log,
	  "cryo_tube_phys", targ_tube_log, false, 0);
  }
 /*  */


  //  Vis attributes
  chamber_inner_log->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes * schamVisAtt
     = new G4VisAttributes(G4Colour(0.7,0.7,1.0));
//      = new G4VisAttributes(G4VisAttributes::Invisible);
  swall_log->SetVisAttributes(schamVisAtt);
  sc_topbottom_log->SetVisAttributes(schamVisAtt);

  G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
  extpipe_log->SetVisAttributes(pipeVisAtt);
  extvac_log->SetVisAttributes( G4VisAttributes::Invisible );

  G4VisAttributes *cryoVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  if( cryo_tube_log ){
      cryo_tube_log->SetVisAttributes(cryoVisAtt);
  }

  G4VisAttributes *winVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  sc_hcalwin_log->SetVisAttributes(winVisAtt);
  sc_bbwin_log->SetVisAttributes(winVisAtt);

  return;
}

void G4SBSDetectorConstruction::ConstructBeamline( G4LogicalVolume *worldlog ){

    G4double density;


    G4Material* aluminum = new G4Material("Aluminum", 13., 26.98*g/mole, density=2.7*g/cm3);
    G4Material* vacuum = new G4Material("vacuum", 1., 1.0*g/mole, density=1e-9*g/cm3);

    G4Element *elO = new G4Element("Oxygen", "O", 8, 16.000*g/mole );
    G4Element *elSi = new G4Element("Silicon", "Si", 14, 28.086*g/mole );
    G4Element *elAl = new G4Element("Aluminum", "Al", 13, 26.982*g/mole );
    G4Element *elCa = new G4Element("Calcium", "Ca", 20, 40.078*g/mole );
    G4Element *elFe = new G4Element("Iron", "Fe", 26, 55.850*g/mole );
    G4Element *elNa = new G4Element("Sodium", "Na", 11, 22.99*g/mole );
    G4Element* elCr  = new G4Element("Chromium","Cr",24.,52.0*g/mole);
    G4Element* elMn   =  new G4Element("Manganese","Mn", 25.,54.94*g/mole);
    G4Element* elNi  = new G4Element("Nickel","Ni",28.,58.70*g/mole);

    density = 2.5*g/cm3;
    G4Material *Concrete = new G4Material("Concrete",density,6);
    Concrete->AddElement(elO, 0.52);
    Concrete->AddElement(elSi, 0.325);
    Concrete->AddElement(elCa, 0.06);
    Concrete->AddElement(elNa, 0.015);
    Concrete->AddElement(elFe, 0.04);
    Concrete->AddElement(elAl, 0.04);

    density = 8.02*g/cm3 ;
    G4Material *stainless = new G4Material("Stainless steel",density,5);
    stainless->AddElement(elMn, 0.02);
    stainless->AddElement(elSi, 0.01);
    stainless->AddElement(elCr, 0.19);
    stainless->AddElement(elNi, 0.10);
    stainless->AddElement(elFe, 0.68);

    density = 8.02*g/cm3 ;
    G4Material* matBe = new G4Material("Berylium", 4., 9.012*g/mole, density=1.85*g/cm3);

    double swallrad = 1.143*m/2;
    double beamheight = 10.0*12*2.54*cm; // 10 feet off the ground

    // Stainless
    G4double ent_len = 10*m;
    G4double ent_rin = 31.75*mm;
    G4double ent_rou = ent_rin+0.120*mm;

    G4Tubs *ent_tube = new G4Tubs("ent_tube", ent_rin, ent_rou, ent_len/2, 0.*deg, 360.*deg );
    G4Tubs *ent_vac  = new G4Tubs("ent_vac", 0.0, ent_rin, ent_len/2, 0.*deg, 360.*deg );

    G4LogicalVolume *entLog = new G4LogicalVolume(ent_tube, stainless, "ent_log", 0, 0, 0);
    G4LogicalVolume *entvacLog = new G4LogicalVolume(ent_vac, vacuum, "entvac_log", 0, 0, 0);


    if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
	// gas target -  1.5m in air
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-1.5*m), entLog, "ent_phys", worldlog, false,0);
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-1.5*m), entvacLog, "entvac_phys", worldlog,false,0);

	// Add in Be window
	G4double winthick = 0.0127*cm;

	G4Tubs *ent_win = new G4Tubs("ent_win", 0.0, ent_rin, winthick/2, 0.*deg, 360.*deg );
	G4LogicalVolume *ent_winlog = new G4LogicalVolume(ent_win, matBe, "entwin_log", 0, 0, 0);
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

    G4LogicalVolume *extLog = new G4LogicalVolume(ext_cone, aluminum, "ext_log", 0, 0, 0);
    G4LogicalVolume *extvacLog = new G4LogicalVolume(ext_vac, vacuum, "extvac_log", 0, 0, 0);

    new G4PVPlacement(0,G4ThreeVector(), extLog, "ext_phys", worldlog, false,0);
    new G4PVPlacement(0,G4ThreeVector(), extvacLog, "extvac_phys", worldlog,false,0);


    double floorthick = 1.0*m;
    G4Tubs *floor_tube = new G4Tubs("floor_tube", 0.0, 30*m, floorthick/2, 0.*deg, 360.*deg );

    G4RotationMatrix *floorrm = new G4RotationMatrix;
    floorrm->rotateX(90*deg);

    G4LogicalVolume *floorLog = new G4LogicalVolume(floor_tube, Concrete, "floor_log", 0, 0, 0);
    new G4PVPlacement(floorrm, G4ThreeVector(0.0, -floorthick/2 - beamheight, 0.0), floorLog, "floor_phys", worldlog, false, 0);


    extvacLog->SetVisAttributes(G4VisAttributes::Invisible);
    entvacLog->SetVisAttributes(G4VisAttributes::Invisible);

    G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));

    extLog->SetVisAttributes(pipeVisAtt);
    entLog->SetVisAttributes(pipeVisAtt);


/*    G4VisAttributes *floorVisAtt= new G4VisAttributes(G4Colour(0.9,0.9,0.9));
    floorLog->SetVisAttributes(floorVisAtt); */
    floorLog->SetVisAttributes(G4VisAttributes::Invisible);

    return;
}


void G4SBSDetectorConstruction::Make48D48( G4LogicalVolume *worldlog, double r48d48 ){
    int nel;
    double z;
    G4String name;

  G4Element *elN = new G4Element("Nitrogen", "N", 7, 14.007*g/mole );
  G4Element *elO = new G4Element("Oxygen", "O", 8, 16.000*g/mole );


    double a = 55.85*g/mole;
  double density = 7.87*g/cm3;
  G4Material* Fe = new G4Material(name="Fer", z=26., a, density);
  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);



  double bigcoilwidth = 214.5*mm;
  double bigcoilheight = 263.7*mm;

  double bigwidth = 2324.1*mm;
  double bigheight = 3721.1*mm;
  double bigdepth = 1219.2*mm;

  double notchdepth = 25*cm;

  G4Box *biggap  = new G4Box("biggap",  469.9*mm/2+0.1*mm, 187.*cm/2.-bigcoilheight,  bigdepth/2+0.1*mm);

  std::vector<G4TwoVector> bigpoly;
  bigpoly.push_back( G4TwoVector(-bigwidth/2.0,  bigdepth/2.0 ));
  bigpoly.push_back( G4TwoVector(-bigwidth/2.0, -bigdepth/2.0 ));
  bigpoly.push_back( G4TwoVector( bigwidth/2.0, -bigdepth/2.0 ));
  bigpoly.push_back( G4TwoVector( bigwidth/2.0, bigdepth/2.0 - notchdepth*sqrt(2.0) ));
  bigpoly.push_back( G4TwoVector( bigwidth/2.0 - notchdepth*sqrt(2.0), bigdepth/2.0  ));

  G4ExtrudedSolid *bigbox_ext = new G4ExtrudedSolid("bigbox_ext", bigpoly, bigheight/2.0, G4TwoVector(), 1.0, G4TwoVector(), 1.0);
  G4RotationMatrix *bigboxrm = new G4RotationMatrix;
  bigboxrm->rotateY(-f48D48ang);
  bigboxrm->rotateX( -90.*deg);
  bigboxrm->rotateZ( 180.*deg);

  G4RotationMatrix *bigboxaddrm = new G4RotationMatrix;
  bigboxaddrm->rotateZ( -180.*deg);
  bigboxaddrm->rotateX( 90.*deg);


  G4Box *bclampgap  = new G4Box("bclampgap",  23.*cm, 65.*cm,  15.5*cm);
  G4Box *fclampgap  = new G4Box("fclampgap",  11.*cm, 35.*cm,  12.*cm/2.);

  G4SubtractionSolid* bigbase = new G4SubtractionSolid("bigbase", bigbox_ext, biggap, bigboxaddrm, G4ThreeVector());

  double coilgapwidth = (60.*cm - bigcoilheight)*2;
  double coilgapheight = 160*cm-bigcoilheight;

  G4Box *bigcoilbase = new G4Box("bigcoilbase", (bigcoilheight+coilgapwidth/2.0)/2.0, bigcoilheight+coilgapheight/2, bigcoilwidth/2.0);
  G4Box *bigcoilgap = new G4Box("bigcoilgap", coilgapwidth/4.0+1.0*mm, coilgapheight/2, bigcoilwidth/2.0+0.1*mm);

//  double coilspace = 6.63*mm;
  double coilspace = 20.63*mm;



  std::vector<G4TwoVector> woundpoly;
  woundpoly.push_back( G4TwoVector(0.0,  -bigdepth/2.0 -coilspace ));
  woundpoly.push_back( G4TwoVector(0.0, -bigdepth/2.0 -coilspace-bigcoilwidth));
  woundpoly.push_back( G4TwoVector( bigwidth/2.0+coilspace+bigcoilwidth, -bigdepth/2.0 -coilspace-bigcoilwidth));
  woundpoly.push_back( G4TwoVector( bigwidth/2.0+coilspace+bigcoilwidth, bigdepth/2.0 - notchdepth*sqrt(2.0) ));

  woundpoly.push_back( G4TwoVector( bigwidth/2.0+coilspace+bigcoilwidth  - 2.0*(coilspace+bigcoilwidth)*sin(pi/8.)*sin(pi/8.) , 
	                               bigdepth/2.0 - notchdepth*sqrt(2.0) + 2.0*(coilspace+bigcoilwidth)*sin(pi/8)*cos(pi/8) ));

  woundpoly.push_back( G4TwoVector( bigwidth/2.0- notchdepth*sqrt(2.0) + 2.0*(coilspace+bigcoilwidth)*sin(pi/8)*cos(pi/8) , 
	                               bigdepth/2.0 + coilspace+bigcoilwidth  - 2.0*(coilspace+bigcoilwidth)*sin(pi/8)*sin(pi/8) ));

  woundpoly.push_back( G4TwoVector( bigwidth/2.0 - notchdepth*sqrt(2.0), bigdepth/2.0 +coilspace+bigcoilwidth ));

  ////
  woundpoly.push_back( G4TwoVector(0.0,  bigdepth/2.0 +coilspace+bigcoilwidth));
  woundpoly.push_back( G4TwoVector(0.0,  bigdepth/2.0 +coilspace ));

  woundpoly.push_back( G4TwoVector( bigwidth/2.0 - notchdepth*sqrt(2.0), bigdepth/2.0 +coilspace ));

  // arc here
  woundpoly.push_back( G4TwoVector( bigwidth/2.0+coilspace  - 2.0*(coilspace)*sin(pi/8.)*sin(pi/8.) , 
	                               bigdepth/2.0 - notchdepth*sqrt(2.0) + 2.0*(coilspace)*sin(pi/8)*cos(pi/8) ));

  woundpoly.push_back( G4TwoVector( bigwidth/2.0- notchdepth*sqrt(2.0) + 2.0*(coilspace)*sin(pi/8)*cos(pi/8) , 
	                               bigdepth/2.0 + coilspace - 2.0*(coilspace)*sin(pi/8)*sin(pi/8) ));

  woundpoly.push_back( G4TwoVector( bigwidth/2.0+coilspace, bigdepth/2.0 - notchdepth*sqrt(2.0) ));

  woundpoly.push_back( G4TwoVector( bigwidth/2.0+coilspace, bigdepth/2.0 +coilspace - notchdepth*sqrt(2.0) ));
  woundpoly.push_back( G4TwoVector( bigwidth/2.0+coilspace, -bigdepth/2.0 -coilspace));

  G4ExtrudedSolid *woundcoil_ext = new G4ExtrudedSolid("woundcoil_ext", woundpoly, bigcoilheight/2.0, G4TwoVector(), 1.0, G4TwoVector(), 1.0);

  // Pull out left side of gap
  G4SubtractionSolid *bigcoil = new G4SubtractionSolid("bigcoil", bigcoilbase, bigcoilgap, 0, G4ThreeVector(-coilgapwidth/4.0, 0.0, 0.0) );

//  double coilfrontback = 150*cm;
  double coilfrontback = bigdepth+2.0*coilspace+bigcoilwidth;

  G4Box *bigcoilthr = new G4Box("bigcoilthr", bigcoilwidth,  bigcoilheight/2,  coilfrontback/2.0 );

  // Sum together coils
  

  // Sum together base iron plus coils

  G4UnionSolid* big48d48;

  G4Box *bigbeamslot = new G4Box("bigbeamslot",  bigwidth/2, 15.5*cm, 2.0*m ); // Height is roughly beam pipe outer radius at 3m
 
  big48d48 = new G4UnionSolid("big48d48_1", bigbase, bigcoilthr, bigboxaddrm, 
	  G4ThreeVector(0.0, 0.0, (coilgapheight+bigcoilheight)/2.0));
  big48d48 = new G4UnionSolid("big48d48_2", big48d48, bigcoilthr, bigboxaddrm, 
	  G4ThreeVector(0.0, 0.0, -(coilgapheight+bigcoilheight)/2.0));

  big48d48 = new G4UnionSolid("big48d48_3", big48d48, bigcoil, bigboxaddrm, 
	  G4ThreeVector(-(bigcoilheight+coilgapwidth/2.0)/2.0-1.0*mm, coilfrontback/2.0, 0.0));
  big48d48 = new G4UnionSolid("big48d48_4", big48d48, bigcoil, bigboxaddrm, 
	  G4ThreeVector(-(bigcoilheight+coilgapwidth/2.0)/2.0-1.0*mm, -coilfrontback/2.0, 0.0));

  big48d48 = new G4UnionSolid("big48d48_5", big48d48, woundcoil_ext, 0, 
	  G4ThreeVector( 1.0*mm, 0.0,  coilgapheight/2.+bigcoilheight/2.0));
  big48d48 = new G4UnionSolid("big48d48_6", big48d48, woundcoil_ext, 0, 
	  G4ThreeVector( 1.0*mm, 0.0,  -coilgapheight/2.-bigcoilheight/2.0));


  //  Cut out slot - from magnet center to inside of cut is ~35cm
  G4SubtractionSolid *big48d48_wslot = new G4SubtractionSolid("big48d48_5", big48d48, bigbeamslot, bigboxaddrm, 
	  G4ThreeVector( bigwidth/2+35*cm, 0.0, 0.0) );

  G4LogicalVolume *big48d48Log=new G4LogicalVolume(big48d48_wslot, Fe,
						  "b48d48Log", 0, 0, 0);

  if( fTotalAbs ){
      big48d48Log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }

  G4RotationMatrix *bigrm = new G4RotationMatrix;
  bigrm->rotateY(-f48D48ang);

  /*
  new G4PVPlacement(bigrm, 
	  G4ThreeVector(-(f48D48dist+bigdepth/2.0)*sin(-f48D48ang), 0.0, (f48D48dist+bigdepth/2.0)*cos(-f48D48ang)),
	  		    big48d48Log, "big48d48Physical", worldlog, 0,false,0);
  */


  new G4PVPlacement(bigboxrm, 
//	  G4ThreeVector(-(f48D48dist+bigdepth/2.0)*sin(-f48D48ang), 0.0, (f48D48dist+bigdepth/2.0)*cos(-f48D48ang)),
	  G4ThreeVector(-r48d48*sin(-f48D48ang), 0.0, r48d48*cos(-f48D48ang)),
	  		    big48d48Log, "big48d48Physical", worldlog, 0,false,0);

  // Associate magnetic field with gap

  G4LogicalVolume *bigfieldLog=new G4LogicalVolume(biggap, Air,
						  "bigfieldLog", 0, 0, 0);

  // use uniform field for now with 48D48

  double fieldValue = 1.4*tesla;
  G4UniformMagField* magField
            = new G4UniformMagField(G4ThreeVector(fieldValue*cos(f48D48ang), 0.0, -fieldValue*sin(f48D48ang)));

  G4FieldManager *bigfm = NULL;

  if( f48d48field ){
      f48d48field->SetRM( bigrm );
      if( f48d48field->GoodParams() ){
	  bigfm = new G4FieldManager(f48d48field);
	  G4Mag_UsualEqRhs* fequation= new G4Mag_UsualEqRhs(fbbfield); 
	  G4MagIntegratorStepper *stepper = new G4ExplicitEuler(fequation, 8);
	  new G4ChordFinder(fbbfield, 1.0*nm, stepper);
	  worldlog->SetFieldManager(bigfm,true);
      }
  } else {
      bigfm = new G4FieldManager(magField);
      bigfm->SetDetectorField(magField);
      bigfm->CreateChordFinder(magField);
      bigfieldLog->SetFieldManager(bigfm,true);
  }

  assert(bigfm);

  new G4PVPlacement(bigrm, 
//	  G4ThreeVector(-(f48D48dist+bigdepth/2.0)*sin(-f48D48ang), 0.0, (f48D48dist+bigdepth/2.0)*cos(-f48D48ang)),
	  G4ThreeVector(-r48d48*sin(-f48D48ang), 0.0, r48d48*cos(-f48D48ang)),
	  		    bigfieldLog, "bigfieldPhysical", worldlog, 0,false,0);

  // Clamps
  
  // The positioning and acceptance gaps were taken directly from CAD
  // The position widths for the beam pipe holes are fudged around so
  // that it doesn't interfere with the beam pipe
  //

  double clampdepth = 10.*cm;
  double clampoffset = 35*cm;

  G4Box *frontclampbase  = new G4Box("frontclampbase", (150+115)*cm/2, bigheight/2,  clampdepth/2);
  G4SubtractionSolid *frontclamp = new G4SubtractionSolid("frontclamp1", frontclampbase, fclampgap, 0,
	 G4ThreeVector( clampoffset, 0,0 ) );

// G4Box *frontclampbeamhole  = new G4Box("frontclampbeamhole", 20.*cm/2, 20.*cm/2,  clampdepth/2+2*cm);
//  frontclamp = new G4SubtractionSolid("frontclamp2", frontclamp, frontclampbeamhole, 0, G4ThreeVector(-55*cm, 0, 0) );
//  frontclamp = new G4SubtractionSolid("frontclamp3", frontclamp, frontclampbeamhole, 0, G4ThreeVector(-28*cm, 0, 0) );

  G4Box *frontclampbeamhole  = new G4Box("frontclampbeamhole", 50.*cm/2, 14.*cm/2,  clampdepth/2+2*cm);
  frontclamp = new G4SubtractionSolid("frontclamp2", frontclamp, frontclampbeamhole, 0, G4ThreeVector(-55*cm+clampoffset, 0, 0) );

  G4Box *frontclampecalhole  = new G4Box("frontclampecalhole", 100.*cm/2, 236.*cm/4.,  clampdepth/2+2*cm);
  frontclamp = new G4SubtractionSolid("frontclamp3", frontclamp, frontclampecalhole, 0, G4ThreeVector(-120*cm+clampoffset, 0, 0) );

  G4LogicalVolume *frontclampLog=new G4LogicalVolume(frontclamp, Fe, "frontclampLog", 0, 0, 0);
  if( fTotalAbs ){
      frontclampLog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }

  double backclampheight = 290.0*cm;

  G4Box *backclampbase  = new G4Box("backclampbase", (150+115)*cm/2, backclampheight/2,  clampdepth/2);

  double backclampaddheight = 55.0*cm;
  G4Box *backclampadd  = new G4Box("backclampadd", (150+115)*cm/2, backclampaddheight/2,  clampdepth);
  G4UnionSolid *backclampfull = new G4UnionSolid("backclampfull1", backclampbase, backclampadd, 0,
	  G4ThreeVector( 0.0, backclampheight/2.0 - backclampaddheight/2.0, clampdepth/2.0-0.1*mm ));
  backclampfull = new G4UnionSolid("backclampfull1", backclampfull, backclampadd, 0,
	  G4ThreeVector( 0.0, -backclampheight/2.0 + backclampaddheight/2.0, clampdepth/2.0 ));

  G4SubtractionSolid *backclamp = new G4SubtractionSolid("backclamp1", backclampfull, bclampgap, 0, 
	 G4ThreeVector(clampoffset, 0, 0) );

  G4Box *backclampbeamhole  = new G4Box("backclampbeamhole", 83.*cm/2., 16.*cm/2.,  clampdepth/2+2*cm);
  backclamp = new G4SubtractionSolid("backclamp2", backclamp, backclampbeamhole, 0, G4ThreeVector(-128*cm+clampoffset, 0, 0) );

  G4LogicalVolume *backclampLog=new G4LogicalVolume(backclamp, Fe, "backclampLog", 0, 0, 0);

  if( fTotalAbs ){
      backclampLog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }

  double frontclampz = -100*cm + clampdepth/2.0;

  new G4PVPlacement(bigrm, 
//	  G4ThreeVector(-(f48D48dist+bigdepth/2.0+frontclampz)*sin(-f48D48ang)-cos(-f48D48ang)*clampoffset, 0.0, (f48D48dist+bigdepth/2.0+frontclampz)*cos(-f48D48ang)-sin(-f48D48ang)*clampoffset),
	  G4ThreeVector(-(r48d48+frontclampz)*sin(-f48D48ang)-cos(-f48D48ang)*clampoffset, 0.0, (r48d48+frontclampz)*cos(-f48D48ang)-sin(-f48D48ang)*clampoffset),
	  		    frontclampLog, "frontclampPhysical", worldlog, 0,false,0);


  /*
   *  No more back clamp in GEN-10M
  double backclampz  =  100*cm - clampdepth/2.0;
  new G4PVPlacement(bigrm, 
	  G4ThreeVector(-(r48d48+backclampz)*sin(-f48D48ang)-cos(-f48D48ang)*clampoffset, 0.0, (r48d48+backclampz)*cos(-f48D48ang)-sin(-f48D48ang)*clampoffset),
	  		    backclampLog, "backclampPhysical", worldlog, 0,false,0);
			    */

  bigfieldLog->SetVisAttributes(G4VisAttributes::Invisible);
//  backclampLog->SetVisAttributes(G4VisAttributes::Invisible);
//  frontclampLog->SetVisAttributes(G4VisAttributes::Invisible);
}



void G4SBSDetectorConstruction::Set48D48Field(int n){
    switch(n){
	case 1:
	    f48d48field = new G4SBS48D48Field( f48D48dist, NULL );
	    break;
	case 0:
	    if( f48d48field ){ delete f48d48field; }
	    f48d48field = NULL;
	    break;
	default:
	    break;
    }
    return;
}










