#include "G4SBSDetectorConstruction.hh"
#include "G4SBSGEMSD.hh"
#include "G4SBSCalSD.hh"

#include "G4UserLimits.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4ProductionCuts.hh"
#include "G4ElementTable.hh"
#include "G4Box.hh"
#include "G4GenericTrap.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
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

    fGEMOption = 1;
}

G4SBSDetectorConstruction::~G4SBSDetectorConstruction()
{;}

G4VPhysicalVolume* G4SBSDetectorConstruction::Construct(){
    // Just nothing so we don't step on toes
    // ConstructAll is where the real magic happens
  G4double a, iz, z, density;
  G4String name, symbol;
  G4int nel;

  G4Element *elN = new G4Element("Nitrogen", "N", 7, 14.007*g/mole );
  G4Element *elO = new G4Element("Oxygen", "O", 8, 16.000*g/mole );

  density = 1e-9*mg/cm3;
  G4Material* Air = new G4Material(name="BlandAir", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);

  G4Box *WorldBox= new G4Box("WorldBox",20*m, 20*m, 20*m);
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
  G4Element *elD = new G4Element("Deuterium", "D", 1, 2.014*g/mole );
  G4Element *el3He = new G4Element("Helium3", "3He", 2, 3.016*g/mole );
  G4Element *elB = new G4Element("Boron", "B", 5, 10.811*g/mole );
  G4Element *elC = new G4Element("Carbon", "C", 6, 12.011*g/mole );
  G4Element *elN = new G4Element("Nitrogen", "N", 7, 14.007*g/mole );
  G4Element *elO = new G4Element("Oxygen", "O", 8, 16.000*g/mole );
  G4Element *elF = new G4Element("Fluorine", "F", 9, 18.998*g/mole );
  G4Element *elMg = new G4Element("Magnesium", "Mg", 12, 24.305*g/mole );
  G4Element *elAl = new G4Element("Aluminum", "Al", 13, 26.982*g/mole );
  G4Element *elSi = new G4Element("Silicon", "Si", 14, 28.086*g/mole );
  G4Element *elK = new G4Element("Potassium", "K", 19, 39.098*g/mole );
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
  G4Material* He  = new G4Material(name="He", z=2., a, density);
  
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


  //--------- Target Materials  ---------

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

  // Mother box
  G4Box *bbmotherBox= new G4Box("bbmotherBox",100*cm, 250*cm, motherdepth/2.0);
  G4LogicalVolume *bbmotherLog=new G4LogicalVolume(bbmotherBox,Air,
						  "bbmotherLog", 0, 0, 0);

  // We need to account for offsets so we can fit BigBite and detectors in without running into
  // the target
  // Need 70 cm clearance from front of the spectrometer to front of the mother volume

  double clear = 70.0*cm;

  double motherr = fBBdist + motherdepth/2.0 - clear;
  
  G4PVPlacement *bbmotherPhys=new G4PVPlacement(bbrm,
	  G4ThreeVector(-motherr*sin(fBBang), 0.0, motherr*cos(fBBang)),
					       bbmotherLog, "bbmotherPhys",
					       WorldLog, 0,false,0);

  double bbmagwidth  = 1670.0*mm;
  double bbmagheight = 486.0*mm;

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
  bbairpts.push_back( G4TwoVector(  (bbmagheight-133.1)*mm, 0.0*mm));
  bbairpts.push_back( G4TwoVector(  bbmagheight, 464.1*mm));
  bbairpts.push_back( G4TwoVector( -bbmagheight, 805.1*mm));
  bbairpts.push_back( G4TwoVector( -bbmagheight, 0.0*mm));

  bbairpts.push_back( G4TwoVector(  (bbmagheight-133.1)*mm, 0.0*mm));
  bbairpts.push_back( G4TwoVector(  bbmagheight, 464.1*mm));
  bbairpts.push_back( G4TwoVector( -bbmagheight, 805.1*mm));
  bbairpts.push_back( G4TwoVector( -bbmagheight, 0.0*mm));

  double gapsize = 250.0*mm;

  G4GenericTrap *bbairTrap = new G4GenericTrap("bbairTrap",
	  gapsize/2.0, bbairpts );


  double coilsize = 320.0*mm;
  double coilwidth = 90.0*mm;

  std::vector<G4TwoVector> bbcoilpts;
  bbcoilpts.push_back( G4TwoVector(  (bbmagheight-133.0)*mm+coilsize, 0.0*mm-coilsize));
  bbcoilpts.push_back( G4TwoVector(  bbmagheight+coilsize, 464.0*mm+coilsize*(1.0-tan(20.0*deg))));
  bbcoilpts.push_back( G4TwoVector( -bbmagheight-coilsize, 805.0*mm+coilsize*(1.0+tan(20.0*deg))));
  bbcoilpts.push_back( G4TwoVector( -bbmagheight-coilsize, 0.0*mm-coilsize));

  bbcoilpts.push_back( G4TwoVector(  (bbmagheight-133.0)*mm+coilsize, 0.0*mm-coilsize));
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

  G4LogicalVolume *bbyokewgapLog=new G4LogicalVolume(yokewgap, Fe,
						  "bbyokewgapLog", 0, 0, 0);

  bbyokewgapLog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );

  G4RotationMatrix *yokerm = new G4RotationMatrix;
  yokerm->rotateY(90.0*deg);
  yokerm->rotateZ(-90.0*deg);
  yokerm->rotateX(180.0*deg);

  G4PVPlacement *bbyokewgapPhys=new G4PVPlacement(yokerm,G4ThreeVector(0.0, 0.0, -motherdepth/2.0+clear),
	  				    bbyokewgapLog,
					       "bbyokewgapPhysical",
					       bbmotherLog,
					       0,false,0);

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
  G4ChordFinder *chordfinder = new G4ChordFinder(fbbfield, 1.0*nm, stepper);
//  bbfm->SetChordFinder(chordfinder);
// bbfm->GetChordFinder()->SetDeltaChord(1.0*um);

  /*
  bbfm->SetMinimumEpsilonStep( 1e-6 );
  bbfm->SetMaximumEpsilonStep( 1e-5 );
  */

  bbmotherLog->SetFieldManager(bbfm,true);

  G4PVPlacement *bbfieldPhys= new G4PVPlacement(0, G4ThreeVector(),
	  		    bbfieldLog, "bbfieldPhysical", bbyokewgapLog, 0,false,0);

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
  G4PVPlacement *bbdetPhys= new G4PVPlacement(bbdetrot, 
	  G4ThreeVector(0.0, 
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

  int gidx, gpidx, lidx; 

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

  char gpname[ngem][nplane][3][255];

//  G4Box **gbox = new G4Box* [ngem*nplane];


  G4Box *gpbox;
  G4LogicalVolume *gplog;
  G4PVPlacement *gpphy;

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
	  gpphy = new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, zrun-gempzsum/2.0), gplog, 
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

  double bigcoilwidth = 231.0*mm;
  double bigcoilheight = 323.0*mm;

  double bigdepth = 1219.2*mm;

  G4Box *bigbox  = new G4Box("bigbox", 2324.1*mm/2, 3721.1*mm/2,  bigdepth/2);
  G4Box *biggap  = new G4Box("biggap",  469.9*mm/2+0.1*mm, 1219.2*mm/2+0.1*mm,  bigdepth/2+0.1*mm);

  G4SubtractionSolid* bigbase = new G4SubtractionSolid("bigbase", bigbox, biggap);
  //G4Box* bigbase = bigbox;

  G4Box *bigcoilbase = new G4Box("bigcoilbase", bigcoilheight+825.5/2.0*mm, 1866.9*mm/2, bigcoilwidth/2.0);
  G4Box *bigcoilgap = new G4Box("bigcoilgap", (825.5/2.0)*mm+1*mm, 1219.2*mm/2+1*mm, bigcoilwidth/2.0+1*mm);

  G4SubtractionSolid* bigcoil = new G4SubtractionSolid("bigcoil", bigcoilbase, bigcoilgap);

  G4Box *bigcoilthr = new G4Box("bigcoilthr", bigcoilwidth,  bigcoilheight/2,  1416.0*mm/2.0+bigcoilwidth );

  // Sum together base iron plus coils

  G4UnionSolid* big48d48;
 
  big48d48 = new G4UnionSolid("big48d48_1", bigbase, bigcoilthr, 0, 
	  G4ThreeVector(0.0, (1219.2*mm+bigcoilheight)/2.0, 0.0));
  big48d48 = new G4UnionSolid("big48d48_2", big48d48, bigcoilthr, 0, 
	  G4ThreeVector(0.0, -(1219.2*mm+bigcoilheight)/2.0, 0.0));
  big48d48 = new G4UnionSolid("big48d48_3", big48d48, bigcoil, 0, 
	  G4ThreeVector(0.0, 0.0, (1416.0*mm+bigcoilwidth)/2.0));
  big48d48 = new G4UnionSolid("big48d48_4", big48d48, bigcoil, 0, 
	  G4ThreeVector(0.0, 0.0, -(1416.0*mm+bigcoilwidth)/2.0));

  G4LogicalVolume *big48d48Log=new G4LogicalVolume(big48d48, Fe,
						  "b48d48Log", 0, 0, 0);

  big48d48Log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );

  G4RotationMatrix *bigrm = new G4RotationMatrix;
  bigrm->rotateY(-f48D48ang);

  G4PVPlacement *big48d48Phys= new G4PVPlacement(bigrm, 
	  G4ThreeVector(-(f48D48dist+bigdepth/2.0)*sin(-f48D48ang), 0.0, (f48D48dist+bigdepth/2.0)*cos(-f48D48ang)),
	  		    big48d48Log, "big48d48Physical", WorldLog, 0,false,0);

  // Associate magnetic field with gap

  G4LogicalVolume *bigfieldLog=new G4LogicalVolume(biggap, Air,
						  "bigfieldLog", 0, 0, 0);

  // use uniform field for now with 48D48

  double fieldValue = 1.4*tesla;
  G4UniformMagField* magField
            = new G4UniformMagField(G4ThreeVector(fieldValue*cos(f48D48ang), 0.0, -fieldValue*sin(f48D48ang)));

  G4FieldManager *bigfm = new G4FieldManager(magField);

  bigfm->SetDetectorField(magField);
  bigfm->CreateChordFinder(magField);

  bigfieldLog->SetFieldManager(bigfm,true);


  G4PVPlacement *bigfieldPhys= new G4PVPlacement(bigrm, 
	  G4ThreeVector(-(f48D48dist+bigdepth/2.0)*sin(-f48D48ang), 0.0, (f48D48dist+bigdepth/2.0)*cos(-f48D48ang)),
	  		    bigfieldLog, "bigfieldPhysical", WorldLog, 0,false,0);

  //--------- Glass target cell -------------------------------

  double wallthick = 1.61*mm;
  double capthick  = 0.126*mm;
  double cellradius    = 0.75*2.54*cm/2.0;

  G4Tubs *targ_tube = new G4Tubs("targ_tube", cellradius-wallthick, cellradius, fTargLen/2.0, 0.*deg, 360.*deg );
  G4Tubs *targ_cap = new G4Tubs("targ_cap", 0.0, cellradius, capthick/2.0, 0.*deg, 360.*deg );

  G4LogicalVolume* targ_tube_log = new G4LogicalVolume(targ_tube, GE180,"targ_tube_log");
  G4LogicalVolume* targ_cap_log = new G4LogicalVolume(targ_cap, GE180,"targ_cap_log");

  G4VPhysicalVolume* targ_cap_phys;

  /* FIXME
   * */
  if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
      G4VPhysicalVolume* targ_tube_phys
      = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), targ_tube_log,
	      "targ_tube_phys", WorldLog, false, 0);

      targ_cap_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, fTargLen/2.0+capthick/2.0), targ_cap_log,
	      "targ_cap_phys1", WorldLog, false, 0);
      targ_cap_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -fTargLen/2.0-capthick/2.0), targ_cap_log,
	      "targ_cap_phys2", WorldLog, false, 0);
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
      G4VPhysicalVolume* gas_tube_phys
	  = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), gas_tube_log,
		  "gas_tube_phys", WorldLog, false, 0);
  }
  /**/
  
  //--------- Cryo target cell -------------------------------

  wallthick   = 178*um;
  double upcapthick    = 71*um;
  double downcapthick  = 102*um;
  cellradius  = 63.5*mm/2.0;

  double swallthick   = 0.38*mm;
  double swallrad     = 1.037*m/2;

  targ_tube = new G4Tubs("targ_tube", cellradius-wallthick, cellradius, fTargLen/2.0, 0.*deg, 360.*deg );
  G4Tubs *targ_ucap = new G4Tubs("targ_ucap", 0.0, cellradius, upcapthick/2.0, 0.*deg, 360.*deg );
  G4Tubs *targ_dcap = new G4Tubs("targ_dcap", 0.0, cellradius, downcapthick/2.0, 0.*deg, 360.*deg );

  targ_tube_log = new G4LogicalVolume(targ_tube, Aluminum,"targ_tube_log");
  G4LogicalVolume* targ_ucap_log = new G4LogicalVolume(targ_cap, Aluminum,"targ_ucap_log");
  G4LogicalVolume* targ_dcap_log = new G4LogicalVolume(targ_cap, Aluminum,"targ_dcap_log");

  G4VPhysicalVolume* targ_ucap_phys, *targ_dcap_phys;

  G4Tubs *swall = new G4Tubs("scham_wall", swallrad, swallrad+swallthick, 0.75*m, 0.*deg, 360.*deg );
  G4LogicalVolume *swall_log = new G4LogicalVolume(swall, Aluminum,"scham_wall_log");

  G4RotationMatrix *schamrot = new G4RotationMatrix;
  schamrot->rotateX(90.0*deg);

  /*
   * FIXME*/
  if( fTargType == kLH2 || fTargType == kLD2 ){
      G4VPhysicalVolume* targ_tube_phys
      = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), targ_tube_log,
	      "targ_tube_phys", WorldLog, false, 0);

      targ_cap_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, fTargLen/2.0+downcapthick/2.0), targ_dcap_log,
	      "targ_dcap_phys", WorldLog, false, 0);
      targ_cap_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -fTargLen/2.0-upcapthick/2.0), targ_ucap_log,
	      "targ_ucap_phys", WorldLog, false, 0);

      // Scattering chamber
      G4VPhysicalVolume* swall_phys
	  = new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), swall_log,
		  "scham_wall_phys", WorldLog, false, 0);
  }
  /**/

  G4Tubs *cryo_tube = new G4Tubs("cryo_tube", 0.0, cellradius-wallthick,fTargLen/2.0, 0.*deg, 360.*deg );
  G4LogicalVolume* cryo_tube_log;


  if( fTargType == kLH2 ){
      cryo_tube_log = new G4LogicalVolume(cryo_tube, LH2mat, "cryo_tube_log");
  }
  if( fTargType == kLD2 ){
      cryo_tube_log = new G4LogicalVolume(cryo_tube, LD2mat, "cryo_tube_log");
  }

  /*
   * FIXME */
  if( fTargType == kLD2 || fTargType == kLH2 ){
      G4VPhysicalVolume* cryo_tube_phys
	  = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), cryo_tube_log,
		  "cryo_tube_phys", WorldLog, false, 0);
  }
 /*  */




  //--------- Calorimeters -------------------------------
  
  // BigBite Preshower

  double psheight = 27*8.5*cm;
  double pswidth  = 2.0*37.0*cm;
  double psdepth  = 8.5*cm;

  G4Box *bbpsbox = new G4Box("bbpsbox", pswidth/2.0, psheight/2.0, psdepth/2.0 );

  G4LogicalVolume* bbpslog = new G4LogicalVolume(bbpsbox, Air, "bbpslog");
  G4VPhysicalVolume* bbpsphys
      = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth/2.0), bbpslog,
	      "bbpsphys", bbdetLog, false, 0, false);

  // BigBite Shower
  
  double calheight = 27*8.5*cm;
  double calwidth  = 7*8.5*cm;
  double caldepth  = 37.0*cm;

  G4Box *bbcalbox = new G4Box("bbcalbox", calwidth/2.0, calheight/2.0, caldepth/2.0 );
  G4LogicalVolume* bbcallog = new G4LogicalVolume(bbcalbox, Lead, "bbcallog");
  

  G4VPhysicalVolume* bbcalphys							       // 5cm clearance
      = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth+caldepth/2.0+5.0*cm), bbcallog,
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
  G4VPhysicalVolume* hcalphys
      = new G4PVPlacement(hcalrm, G4ThreeVector(hcalr*sin(f48D48ang), 0.0, hcalr*cos(f48D48ang) ), hcallog,
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
  bigfieldLog->SetVisAttributes(G4VisAttributes::Invisible);

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

  G4VisAttributes * schamVisAtt
    = new G4VisAttributes(G4Colour(0.7,0.7,1.0));
  swall_log->SetVisAttributes(schamVisAtt);

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
  G4Element *elD = new G4Element("Deuterium", "D", 1, 2.014*g/mole );
  G4Element *el3He = new G4Element("Helium3", "3He", 2, 3.016*g/mole );
  G4Element *elB = new G4Element("Boron", "B", 5, 10.811*g/mole );
  G4Element *elC = new G4Element("Carbon", "C", 6, 12.011*g/mole );
  G4Element *elN = new G4Element("Nitrogen", "N", 7, 14.007*g/mole );
  G4Element *elO = new G4Element("Oxygen", "O", 8, 16.000*g/mole );
  G4Element *elF = new G4Element("Fluorine", "F", 9, 18.998*g/mole );
  G4Element *elMg = new G4Element("Magnesium", "Mg", 12, 24.305*g/mole );
  G4Element *elAl = new G4Element("Aluminum", "Al", 13, 26.982*g/mole );
  G4Element *elSi = new G4Element("Silicon", "Si", 14, 28.086*g/mole );
  G4Element *elK = new G4Element("Potassium", "K", 19, 39.098*g/mole );
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
  G4Material* He  = new G4Material(name="He", z=2., a, density);
  
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

  G4double density_CH2 = 0.95*g/cm3;
  G4Material* CH2 = new G4Material("CH2", density_CH2, nel=2);
  CH2->AddElement(elC, 1);
  CH2->AddElement(elH, 2);

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


  //--------- Target Materials  ---------

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

  // Ecal will act as Hcal detector

  double bigcalheight = (24*4.5+32*4.0)*cm;
  double bigcalwidth  = 44.10*2.54*cm;
  double bigcaldepth  = 15.75*2.54*cm;
  double bbr = fBBdist+bigcaldepth/2.0;

  G4Box *bigcalbox = new G4Box("bigcalbox", bigcalwidth/2.0, bigcalheight/2.0, bigcaldepth/2.0 );
  G4LogicalVolume* bigcallog = new G4LogicalVolume(bigcalbox, Lead, "bigcallog");
  G4VPhysicalVolume* bigcalphys
      = new G4PVPlacement(bbrm, G4ThreeVector(bbr*sin(-fBBang), 0.0, bbr*cos(-fBBang) ), bigcallog,
	      "bigcalphys", WorldLog, false, 0, false);

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

  // SBS box

  double sbsdepth  = 6.0*m;
  double sbswidth  = 2.0*m;
  double sbsheight = 3.5*m;

  double sbsr = fHCALdist+sbsdepth/2;
  G4Box *sbsbox = new G4Box("sbsbox", sbswidth/2.0, sbsheight/2.0, sbsdepth/2.0 );
  G4LogicalVolume* sbslog = new G4LogicalVolume(sbsbox, Air, "sbslog");
  G4VPhysicalVolume* sbsphys
      = new G4PVPlacement(hcalrm, G4ThreeVector(sbsr*sin(f48D48ang), 0.0, sbsr*cos(f48D48ang) ), sbslog,
	      "sbsphys", WorldLog, false, 0, false);

  //  6 GEMs in the front tracker

  double detoffset = 0.05*m - sbsdepth/2.0;

  const int nplane = 24;

  int i;
  int ngem = 0;
  double gemdsep;
  double *gemz;
  double *gemw;
  double *gemh;

  ngem = 14;
  gemdsep = 0.05*m;

  double pairspac  = 0.30*m;
  //
  // GEM option 1
  
  gemz = new double[ngem];
  gemw = new double[ngem];
  gemh = new double[ngem];

  for( i = 0; i < ngem; i++ ){
      if( i < 6 ){
	  gemz[i] = pairspac*(i/2) + (i%2)*gemdsep;
	  gemw[i] = 40.0*cm;
	  gemh[i] = 150.0*cm;
      } else if( i < 10 ) {
	  gemz[i] = pairspac*((i-6)/2) + (i%2)*gemdsep + 1.32*m + 7*cm;
	  gemw[i] = 50.0*cm;
	  gemh[i] = 200.0*cm;
      } else {
	  gemz[i] = pairspac*((i-10)/2) + (i%2)*gemdsep + 2.43*m + 7*cm;
	  gemw[i] = 50.0*cm;
	  gemh[i] = 200.0*cm;
      }

      printf("i = %d  z = %f\n", i, gemz[i]/m);
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

  int gidx, gpidx, lidx; 

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

  char gpname[ngem][nplane][3][255];

//  G4Box **gbox = new G4Box* [ngem*nplane];


  G4Box *gpbox;
  G4LogicalVolume *gplog;
  G4PVPlacement *gpphy;

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
	  gpphy = new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, zrun-gempzsum/2.0), gplog, 
		  gpname[gidx][gpidx][2], gemlog, 0, 0, false );
	  zrun += gempz[gpidx]/2.0;

	  // Assign sensitive volume
	  if( gpidx == 5 && gidx < 6 ){ // just make first six tracking for now
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

  G4VPhysicalVolume* anaphys;
    anaphys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset +  73.7*cm + anadepth/2.0 ), analog,
	      "anaphys1", sbslog, false, 0, false);
    anaphys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset + 185.4*cm + anadepth/2.0 ), analog,
	      "anaphys1", sbslog, false, 0, false);

  //  HCAL
  
  double hcalheight = 330.0*cm;
  double hcalwidth  = 165.0*cm;
  double hcaldepth  = 101.0*cm;

  G4Box *hcalbox = new G4Box("hcalbox", hcalwidth/2.0, hcalheight/2.0, hcaldepth/2.0 );
  G4LogicalVolume* hcallog = new G4LogicalVolume(hcalbox, Lead, "hcallog");
  G4VPhysicalVolume* hcalphys
      = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset + 2.97*m + hcaldepth/2.0 ), hcallog,
	      "hcalphys", sbslog, false, 0, false);

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

  double bigcoilwidth = 231.0*mm;
  double bigcoilheight = 323.0*mm;

  double bigdepth = 1219.2*mm;

  G4Box *bigbox  = new G4Box("bigbox", 2324.1*mm/2, 3721.1*mm/2,  bigdepth/2);
  G4Box *biggap  = new G4Box("biggap",  469.9*mm/2+0.1*mm, 1219.2*mm/2+0.1*mm,  bigdepth/2+0.1*mm);

  G4SubtractionSolid* bigbase = new G4SubtractionSolid("bigbase", bigbox, biggap);
  //G4Box* bigbase = bigbox;

  G4Box *bigcoilbase = new G4Box("bigcoilbase", bigcoilheight+825.5/2.0*mm, 1866.9*mm/2, bigcoilwidth/2.0);
  G4Box *bigcoilgap = new G4Box("bigcoilgap", (825.5/2.0)*mm+1*mm, 1219.2*mm/2+1*mm, bigcoilwidth/2.0+1*mm);

  G4SubtractionSolid* bigcoil = new G4SubtractionSolid("bigcoil", bigcoilbase, bigcoilgap);

  G4Box *bigcoilthr = new G4Box("bigcoilthr", bigcoilwidth,  bigcoilheight/2,  1416.0*mm/2.0+bigcoilwidth );

  // Sum together base iron plus coils

  G4UnionSolid* big48d48;
 
  big48d48 = new G4UnionSolid("big48d48_1", bigbase, bigcoilthr, 0, 
	  G4ThreeVector(0.0, (1219.2*mm+bigcoilheight)/2.0, 0.0));
  big48d48 = new G4UnionSolid("big48d48_2", big48d48, bigcoilthr, 0, 
	  G4ThreeVector(0.0, -(1219.2*mm+bigcoilheight)/2.0, 0.0));
  big48d48 = new G4UnionSolid("big48d48_3", big48d48, bigcoil, 0, 
	  G4ThreeVector(0.0, 0.0, (1416.0*mm+bigcoilwidth)/2.0));
  big48d48 = new G4UnionSolid("big48d48_4", big48d48, bigcoil, 0, 
	  G4ThreeVector(0.0, 0.0, -(1416.0*mm+bigcoilwidth)/2.0));

  G4LogicalVolume *big48d48Log=new G4LogicalVolume(big48d48, Fe,
						  "b48d48Log", 0, 0, 0);

  big48d48Log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );

  G4RotationMatrix *bigrm = new G4RotationMatrix;
  bigrm->rotateY(-f48D48ang);

  G4PVPlacement *big48d48Phys= new G4PVPlacement(bigrm, 
	  G4ThreeVector(-(f48D48dist+bigdepth/2.0)*sin(-f48D48ang), 0.0, (f48D48dist+bigdepth/2.0)*cos(-f48D48ang)),
	  		    big48d48Log, "big48d48Physical", WorldLog, 0,false,0);

  // Associate magnetic field with gap

  G4LogicalVolume *bigfieldLog=new G4LogicalVolume(biggap, Air,
						  "bigfieldLog", 0, 0, 0);

  // use uniform field for now with 48D48

  double fieldValue = 1.4*tesla;
  G4UniformMagField* magField
            = new G4UniformMagField(G4ThreeVector(fieldValue*cos(f48D48ang), 0.0, -fieldValue*sin(f48D48ang)));

  G4FieldManager *bigfm = new G4FieldManager(magField);

  bigfm->SetDetectorField(magField);
  bigfm->CreateChordFinder(magField);

  bigfieldLog->SetFieldManager(bigfm,true);


  G4PVPlacement *bigfieldPhys= new G4PVPlacement(bigrm, 
	  G4ThreeVector(-(f48D48dist+bigdepth/2.0)*sin(-f48D48ang), 0.0, (f48D48dist+bigdepth/2.0)*cos(-f48D48ang)),
	  		    bigfieldLog, "bigfieldPhysical", WorldLog, 0,false,0);

  //--------- Glass target cell -------------------------------

  double wallthick = 1.61*mm;
  double capthick  = 0.126*mm;
  double cellradius    = 0.75*2.54*cm/2.0;

  G4Tubs *targ_tube = new G4Tubs("targ_tube", cellradius-wallthick, cellradius, fTargLen/2.0, 0.*deg, 360.*deg );
  G4Tubs *targ_cap = new G4Tubs("targ_cap", 0.0, cellradius, capthick/2.0, 0.*deg, 360.*deg );

  G4LogicalVolume* targ_tube_log = new G4LogicalVolume(targ_tube, GE180,"targ_tube_log");
  G4LogicalVolume* targ_cap_log = new G4LogicalVolume(targ_cap, GE180,"targ_cap_log");

  G4VPhysicalVolume* targ_cap_phys;

  if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
      G4VPhysicalVolume* targ_tube_phys
      = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), targ_tube_log,
	      "targ_tube_phys", WorldLog, false, 0);

      targ_cap_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, fTargLen/2.0+capthick/2.0), targ_cap_log,
	      "targ_cap_phys1", WorldLog, false, 0);
      targ_cap_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -fTargLen/2.0-capthick/2.0), targ_cap_log,
	      "targ_cap_phys2", WorldLog, false, 0);
  }



  // gas
  G4Tubs *gas_tube = new G4Tubs("gas_tube", 0.0, cellradius-wallthick,fTargLen/2.0, 0.*deg, 360.*deg );
  G4LogicalVolume* gas_tube_log;


  if( fTargType == kH2 || fTargType == kNeutTarg ){
      gas_tube_log = new G4LogicalVolume(gas_tube, refH2, "gas_tube_log");
  }
  if( fTargType == k3He ){
      gas_tube_log = new G4LogicalVolume(gas_tube, pol3He, "gas_tube_log");
  }

  if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
      G4VPhysicalVolume* gas_tube_phys
	  = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), gas_tube_log,
		  "gas_tube_phys", WorldLog, false, 0);
  }
  
  //--------- Cryo target cell -------------------------------

  /*
  wallthick   = 178*um;
  double upcapthick    = 71*um;
  double downcapthick  = 102*um;
  cellradius  = 63.5*mm/2.0;
  */
  // From CDR
  wallthick   = 100*um;
  double upcapthick    = 50*um;
  double downcapthick  = 50*um;

  double swallthick   = 0.38*mm;
  double swallrad     = 1.037*m/2;

  targ_tube = new G4Tubs("targ_tube", cellradius-wallthick, cellradius, fTargLen/2.0, 0.*deg, 360.*deg );
  G4Tubs *targ_ucap = new G4Tubs("targ_ucap", 0.0, cellradius, upcapthick/2.0, 0.*deg, 360.*deg );
  G4Tubs *targ_dcap = new G4Tubs("targ_dcap", 0.0, cellradius, downcapthick/2.0, 0.*deg, 360.*deg );

  targ_tube_log = new G4LogicalVolume(targ_tube, Aluminum,"targ_tube_log");
  G4LogicalVolume* targ_ucap_log = new G4LogicalVolume(targ_cap, Aluminum,"targ_ucap_log");
  G4LogicalVolume* targ_dcap_log = new G4LogicalVolume(targ_cap, Aluminum,"targ_dcap_log");

  G4VPhysicalVolume* targ_ucap_phys, *targ_dcap_phys;

  G4Tubs *swall = new G4Tubs("scham_wall", swallrad, swallrad+swallthick, 0.75*m, 0.*deg, 360.*deg );
  G4LogicalVolume *swall_log = new G4LogicalVolume(swall, Aluminum,"scham_wall_log");

  G4RotationMatrix *schamrot = new G4RotationMatrix;
  schamrot->rotateX(90.0*deg);

  if( fTargType == kLH2 || fTargType == kLD2 ){
      G4VPhysicalVolume* targ_tube_phys
      = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), targ_tube_log,
	      "targ_tube_phys", WorldLog, false, 0);

      targ_cap_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, fTargLen/2.0+downcapthick/2.0), targ_dcap_log,
	      "targ_dcap_phys", WorldLog, false, 0);
      targ_cap_phys = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -fTargLen/2.0-upcapthick/2.0), targ_ucap_log,
	      "targ_ucap_phys", WorldLog, false, 0);

      // Scattering chamber
      G4VPhysicalVolume* swall_phys
	  = new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), swall_log,
		  "scham_wall_phys", WorldLog, false, 0);
  }

  G4Tubs *cryo_tube = new G4Tubs("cryo_tube", 0.0, cellradius-wallthick,fTargLen/2.0, 0.*deg, 360.*deg );
  G4LogicalVolume* cryo_tube_log;


  if( fTargType == kLH2 ){
      cryo_tube_log = new G4LogicalVolume(cryo_tube, LH2mat, "cryo_tube_log");
  }
  if( fTargType == kLD2 ){
      cryo_tube_log = new G4LogicalVolume(cryo_tube, LD2mat, "cryo_tube_log");
  }

  if( fTargType == kLD2 || fTargType == kLH2 ){
      G4VPhysicalVolume* cryo_tube_phys
	  = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), cryo_tube_log,
		  "cryo_tube_phys", WorldLog, false, 0);
  }





  //--------- Visualization attributes -------------------------------
  WorldLog->SetVisAttributes(G4VisAttributes::Invisible);
  sbslog->SetVisAttributes(G4VisAttributes::Invisible);
  bigfieldLog->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes * yokeVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));

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

  G4VisAttributes * anaVisAtt
    = new G4VisAttributes(G4Colour(0.6,0.6,0.0));
  analog->SetVisAttributes(anaVisAtt);

  G4VisAttributes * schamVisAtt
    = new G4VisAttributes(G4Colour(0.7,0.7,1.0));
  swall_log->SetVisAttributes(schamVisAtt);

  return WorldPhys;
}


