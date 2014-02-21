#include "G4SBSDetectorConstruction.hh"
#include "G4SBSGEMSD.hh"
#include "G4SBSCalSD.hh"
#include "G4SBSRICHSD.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UserLimits.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4MaterialPropertiesTable.hh" 
#include "G4NistManager.hh"
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
#include "G4IntersectionSolid.hh"
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
#include "G4FieldManager.hh"

#include "G4MagIntegratorStepper.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4ChordFinder.hh"
#include "G4SimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include <vector>
#include <map>
//#include <pair>

using namespace std;

G4SBSDetectorConstruction::G4SBSDetectorConstruction()
{
    fBBang  = 40.0*deg;
    fBBdist = 1.5*m;

    fBBCaldist = 0.8*m;

    f48D48ang  = 39.4*deg;
    f48D48dist = 2.8*m;

    f48D48_uniform_bfield = 1.4*tesla;
    f48D48_fieldclamp_config = 1; //0 = No field clamps. 1 = GEp (default). 2 = BigBite experiments:

    fHCALdist  = 17.0*m;

    fRICHdist  = 15.0*m;

    fTargLen = 60.0*cm;

    fTargType = kH2;
    fTargDen = 10.5*atmosphere/(300*kelvin*k_Boltzmann);

    fCerDepth = 60.0*cm;
    fCerDist  =  7.0*cm;

    fGEMDist  = 70.0*cm;

    fbbfield = new G4SBSBigBiteField( fBBdist, NULL );

    fGEMOption = 1;

    fTotalAbs = true;

    fExpType = kNeutronExp;

    ConstructMaterials(); //Now we want to construct all materials at the beginning, so that the physics tables can get built properly!!!
    
}

G4SBSDetectorConstruction::~G4SBSDetectorConstruction()
{;}

G4VPhysicalVolume* G4SBSDetectorConstruction::Construct(){
    // Just nothing so we don't step on toes
    // ConstructAll is where the real magic happens
  //G4double a, iz, z, density;
  //Moved all material definitions to ConstructMaterials()

  if( MaterialsMap.empty() ) ConstructMaterials();

  G4Material *Mtemp = GetMaterial("BlandAir");

  G4Box *WorldBox= new G4Box("WorldBox",50*m, 50*m, 50*m);
  G4LogicalVolume *WorldLog=new G4LogicalVolume(WorldBox,Mtemp,
						  "WorldLogical", 0, 0, 0);
  G4PVPlacement *WorldPhys=new G4PVPlacement(0,G4ThreeVector(),
					       "WorldPhysical",
					       WorldLog,
					       0,false,0);
  return WorldPhys;
}

void G4SBSDetectorConstruction::ConstructMaterials(){

  G4NistManager *man = G4NistManager::Instance();

  G4MaterialPropertiesTable *MPT_temp; //pointer to hold material optical properties:

  G4double fractionmass;

  G4double density;
  G4String name, symbol;
  G4int nel, natoms;
  
  map<G4String,G4Material*>::iterator itest;
  //pair<G4String, G4Material>  ptest;

  G4Element *elN = new G4Element("Nitrogen", "N", 7, 14.007*g/mole );
  G4Element *elO = new G4Element("Oxygen", "O", 8, 16.000*g/mole );
  
  density = 1e-9*mg/cm3;
  G4Material* BlandAir = new G4Material(name="BlandAir", density, nel=2);
  BlandAir->AddElement(elN, .7);
  BlandAir->AddElement(elO, .3);

  MaterialsMap["BlandAir"] = BlandAir;

  G4double a, iz, z;
  
  //Seamus's element definitions:
  a = 126.9*g/mole;
  G4Element* elI = new G4Element(name="Iodine",   symbol="I",  iz=53., a);
  a = 132.9*g/mole;
  G4Element* elCs= new G4Element(name="Cesium",   symbol="Cs", iz=55., a);
  
  G4Element *elH = new G4Element("Hydrogen", "H", 1, 1.007*g/mole );
  G4Element *elD = new G4Element("Deuterium", "D", 1, 2.014*g/mole );
  G4Element *el3He = new G4Element("Helium3", "3He", 2, 3.016*g/mole );
  G4Element *elC = new G4Element("Carbon", "C", 6, 12.011*g/mole );
  //G4Element *elN = new G4Element("Nitrogen", "N", 7, 14.007*g/mole );
  //G4Element *elO = new G4Element("Oxygen", "O", 8, 16.000*g/mole );
  G4Element *elF = new G4Element("Fluorine", "F", 9, 18.998*g/mole );
  G4Element *elNa = new G4Element("Sodium", "Na", 11, 22.99*g/mole );
  G4Element *elAl = new G4Element("Aluminum", "Al", 13, 26.982*g/mole );
  G4Element *elSi = new G4Element("Silicon", "Si", 14, 28.086*g/mole );
  G4Element *elCa = new G4Element("Calcium", "Ca", 20, 40.078*g/mole );
  G4Element *elFe = new G4Element("Iron", "Fe", 26, 55.850*g/mole );
  G4Element *elSr = new G4Element("Strontium", "Sr", 38, 87.62*g/mole );
  G4Element *elBa = new G4Element("Barium", "Ba", 56, 137.327*g/mole );
  
  G4Element* elCl  = new G4Element("Chlorine",  "Cl", z=17, a=   35.453*g/mole);
  G4Element* elAr  = new G4Element("Argon",     "Ar", z=18, a=    39.95*g/mole);

  G4Element* elCr  = new G4Element("Chromium","Cr",24.,52.0*g/mole);
  G4Element* elMn   =  new G4Element("Manganese","Mn", 25.,54.94*g/mole);
  G4Element* elNi  = new G4Element("Nickel","Ni",28.,58.70*g/mole);
  
  MaterialsMap["Vacuum"] = new G4Material(name="Vacuum", z=1., a=1.0*g/mole, density=1e-9*g/cm3);

  //G4Material* Lead  = new G4Material(name="Lead", z=82., a=208.0*g/mole, density=11.34*g/cm3);
  //G4Material* Aluminum = new G4Material(name="Aluminum", z=13., a=26.98*g/mole, density=2.7*g/cm3);

  if( MaterialsMap.find("Lead") == MaterialsMap.end() ) MaterialsMap["Lead"] = new G4Material(name="Lead", z=82., a=208.0*g/mole, density=11.34*g/cm3);
  if( MaterialsMap.find("Aluminum") == MaterialsMap.end() ) MaterialsMap["Aluminum"] = new G4Material(name="Aluminum", z=13., a=26.98*g/mole, density=2.7*g/cm3);
    
  density = 4.51*g/cm3;
  G4Material* CsI = new G4Material(name="CsI", density, nel = 2);
  CsI->AddElement(elI, .5);
  CsI->AddElement(elCs,.5);

  MaterialsMap["CsI"] = CsI;
 
  // a = 4.0*g/mole;
  // density = 0.1786e-03*g/cm3;
  

  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  MaterialsMap["Fer"] = new G4Material(name="Fer", z=26., a, density);

  density = 1.29e-03*g/cm3;
  G4Material* Air = new G4Material(name="Air", density, nel=2);
  Air->AddElement(elN, .7);
  Air->AddElement(elO, .3);
    
  MaterialsMap["Air"] = Air;

  double bigden = 1e9*g/cm3;

  // Cell Glass - GE180 Aluminosilicate Glass
  //Changed names of materials in this composition since the G4Materials used here are only used to make up GE180:
  // SiO2 60.3%
  G4Material* SiO2 = new G4Material("GE180_SiO2", 2.2*g/cm3, 2 );
  SiO2->AddElement(elSi, 1);
  SiO2->AddElement(elO, 2);
  MaterialsMap["GE180_SiO2"] = SiO2;
  
  // BaO  18.2%
  G4Material* BaO = new G4Material("GE180_BaO", bigden, 2 );
  BaO->AddElement(elBa, 1);
  BaO->AddElement(elO, 1);
  MaterialsMap["GE180_BaO"] = BaO;
  // Al2O3 14.3%
  G4Material* Al2O3 = new G4Material("GE180_Al2O3", bigden, 2 );
  Al2O3->AddElement(elAl, 2);
  Al2O3->AddElement(elO, 3);
  MaterialsMap["GE180_Al2O3"] = Al2O3;
  // CaO   6.5%
  G4Material* CaO = new G4Material("GE180_CaO", bigden, 2 );
  CaO->AddElement(elCa, 1);
  CaO->AddElement(elO, 1);
  MaterialsMap["GE180_CaO"] = CaO;
  // SrO   0.25%
  G4Material* SrO = new G4Material("GE180_SrO", bigden, 2 );
  SrO->AddElement(elSr, 1);
  SrO->AddElement(elO, 1);
  MaterialsMap["GE180_SrO"] = SrO;
  
  // Density 2.76 g/cm^3
  // Index of Refraction 1.536
  G4Material* GE180 = new G4Material("GE180", 2.76*g/cm3, 5);
  GE180->AddMaterial(SiO2, 0.6039);
  GE180->AddMaterial(BaO, 0.1829);
  GE180->AddMaterial(Al2O3, 0.1439);
  GE180->AddMaterial(CaO, 0.0659);
  GE180->AddMaterial(SrO, 0.0034);
  MaterialsMap["GE180"] = GE180;
  //
  density = 1.19*g/cm3;
  G4Material* Acrylic = new G4Material(name="Acrylic", density, nel=3);
  Acrylic->AddElement(elC, 5);
  Acrylic->AddElement(elH, 8);
  Acrylic->AddElement(elO, 2);
  MaterialsMap["Acrylic"] = Acrylic;
  
  
  //--------- GEM Materials  ---------
  // (stolen from GEMC)

  G4Material* NOMEX_pure = new G4Material("NOMEX_pure", density = 1.38*g/cm3, 5);
  NOMEX_pure -> AddElement(elH,0.04);
  NOMEX_pure -> AddElement(elC,0.54);
  NOMEX_pure -> AddElement(elN,0.09);
  NOMEX_pure -> AddElement(elO,0.10);
  NOMEX_pure -> AddElement(elCl,0.23);
  MaterialsMap["NOMEX_pure"] = NOMEX_pure;

  G4Material* NOMEX = new G4Material("NOMEX",density = 0.04*g/cm3, 2);
  NOMEX -> AddMaterial(NOMEX_pure,0.45);
  NOMEX -> AddMaterial(Air,0.55);
  MaterialsMap["NOMEX"] = NOMEX;

  G4Material* NEMAG10 = new G4Material("NEMAG10", 1.70*g/cm3, nel=4);
  NEMAG10 -> AddElement(elSi, 1);
  NEMAG10 -> AddElement(elO , 2);
  NEMAG10 -> AddElement(elC , 3);
  NEMAG10 -> AddElement(elH , 3);
  MaterialsMap["NEMAG10"] = NEMAG10;

  G4double density_Ar = 1.7823*mg/cm3 ;
  G4Material* Argon = new G4Material("Argon"  , density_Ar, nel=1);
  Argon->AddElement(elAr, 1);
  MaterialsMap["Argon"] = Argon;

  G4double density_CO2 = 1.977*mg/cm3;
  G4Material* CO2 = new G4Material("CO2", density_CO2, nel=2);
  CO2->AddElement(elC, 1);
  CO2->AddElement(elO, 2);
  MaterialsMap["CO2"] = CO2;

  // 1.5 Atmosphere C4F8O for cerkenkov
  G4double density_C4F8O = 9.64*mg/cm3; // density at 1ATM
  G4Material* C4F8O = new G4Material("C4F8O", density_C4F8O*1.5, nel=3);
  C4F8O->AddElement(elC, 4);
  C4F8O->AddElement(elF, 8);
  C4F8O->AddElement(elO, 1);
  MaterialsMap["C4F8O"] = C4F8O;

  G4double density_ArCO2 = .7*density_Ar + .3*density_CO2;
  // Use ArCO2
  G4Material *GEMgas= new G4Material("GEMgas", density_ArCO2, nel=2);
  GEMgas->AddMaterial(Argon, 0.7*density_Ar/density_ArCO2) ;
  GEMgas->AddMaterial(CO2, 0.3*density_CO2/density_ArCO2) ;

  MaterialsMap["GEMgas"] = GEMgas;

  G4Material *Copper= new G4Material("Copper", z=29, a=   63.55*g/mole, density = 8.96*g/cm3);
  MaterialsMap["Copper"] = Copper;

  G4Material *Kapton = new G4Material("Kapton",   density = 1.42*g/cm3, nel=4);
  Kapton->AddElement(elH, 0.026362);
  Kapton->AddElement(elC, 0.691133);
  Kapton->AddElement(elN, 0.073270);
  Kapton->AddElement(elO, 0.209235);

  MaterialsMap["Kapton"] = Kapton;

  //CH2 for FPP analyzers for GEP:
  G4double density_CH2 = 0.95*g/cm3;
  G4Material* CH2 = new G4Material("CH2", density_CH2, nel=2);
  CH2->AddElement(elC, 1);
  CH2->AddElement(elH, 2);

  MaterialsMap["CH2"] = CH2;

  G4double density_CH = 0.95*g/cm3;
  G4Material* CH = new G4Material("CH", density_CH, nel=2);
  CH->AddElement(elC, 1);
  CH->AddElement(elH, 1);

  MaterialsMap["CH"] = CH;

  //Target materials:
  double gasden = fTargDen;
  //  gasden = 10.5*atmosphere/(300*kelvin*k_Boltzmann);
  G4Material *refH2 = new G4Material("refH2", gasden, 1 );
  refH2->AddElement(elH, 1);
  
  MaterialsMap["refH2"] = refH2;

  // gasden = 10.5*atmosphere/(300*kelvin*k_Boltzmann);
  G4Material *refN2 = new G4Material("refN2", gasden, 1 );
  refN2->AddElement(elN, 1);
  
  MaterialsMap["refN2"] = refN2;

  // gasden = 10.77*atmosphere/(300*kelvin*k_Boltzmann);
  G4Material *pol3He = new G4Material("pol3He", gasden, 1 );
  pol3He->AddElement(el3He, 1);

  MaterialsMap["pol3He"] = pol3He;

  double LH2den = 0.071*g/cm3;
  G4Material *LH2mat = new G4Material("LH2", LH2den, 1 );
  LH2mat->AddElement(elH, 1);

  MaterialsMap["LH2"] = LH2mat;

  double LD2den = 162.4*kg/m3;
  G4Material *LD2mat = new G4Material("LD2", LD2den, 1 );
  LD2mat->AddElement(elD, 1);

  MaterialsMap["LD2"] = LD2mat;

  //Beamline materials:
  density = 2.5*g/cm3;
  G4Material *Concrete = new G4Material("Concrete",density,6);
  Concrete->AddElement(elO, 0.52);
  Concrete->AddElement(elSi, 0.325);
  Concrete->AddElement(elCa, 0.06);
  Concrete->AddElement(elNa, 0.015);
  Concrete->AddElement(elFe, 0.04);
  Concrete->AddElement(elAl, 0.04);
  
  MaterialsMap["Concrete"] = Concrete;

  density = 8.02*g/cm3 ;
  G4Material *stainless = new G4Material("Stainless steel",density,5);
  stainless->AddElement(elMn, 0.02);
  stainless->AddElement(elSi, 0.01);
  stainless->AddElement(elCr, 0.19);
  stainless->AddElement(elNi, 0.10);
  stainless->AddElement(elFe, 0.68);

  MaterialsMap["Stainless"] = stainless;

  density = 8.02*g/cm3 ;
  G4Material* matBe = new G4Material("Berylium", 4., 9.012*g/mole, density=1.85*g/cm3);

  MaterialsMap["Beryllium"] = matBe;

  //RICH Materials:
  // Materials for SBS RICH: 
  // List of needed materials:
  // 1. Aerogel tiles
  // 2. Spherical mirror: Aluminum metal "skin surface" and carbon fiber composite structural backing
  // 3. PMT photocathodes and light collection cones: Specify quantum efficiency, index of refraction, etc. 
  // 4. Radiator gas: C4F10 or other?
  // 5. Entry and exit windows: Aluminum & UVT lucite:

  G4Element *H  = man->FindOrBuildElement("H");
  G4Element *Si = man->FindOrBuildElement("Si");
  G4Element *O  = man->FindOrBuildElement("O");
  
  G4Element *C  = man->FindOrBuildElement("C");
  G4Element *Al = man->FindOrBuildElement("Al");
  G4Element *F  = man->FindOrBuildElement("F");

  //Add sodium, potassium, and Boron for Borosilicate glass windows of PMTs:
  //Including the borosilicate glass in the PMT windows is important to measure the contribution of neutron backgrounds and the signal rate induced by the B(n,alpha) reaction.
  G4Element *B = man->FindOrBuildElement("B");
  G4Element *Na = man->FindOrBuildElement("Na");
  G4Element *K  = man->FindOrBuildElement("K");
  
  G4Element *Sb = man->FindOrBuildElement("Sb");
  G4Element *Cs = man->FindOrBuildElement("Cs");
  
  //Define AIR for RICH mother volume:
  G4Material *RICH_Air = man->FindOrBuildMaterial("G4_AIR");
  MaterialsMap["RICH_Air"] = RICH_Air;

  //Define tedlar for gaps between aerogel tiles. These will have no optical properties and therefore optical photons 
  //crossing tile boundaries will be killed.
  G4Material *Tedlar = man->FindOrBuildMaterial("G4_POLYVINYLIDENE_FLUORIDE");
  MaterialsMap["Tedlar"] = Tedlar;

  G4double den_lucite = 1.19*g/cm3;
  
  // Composition of "lucite" according to NIST:
  G4Material *UVT_Lucite = new G4Material("UVT_Lucite", den_lucite, nel=3 );
  UVT_Lucite->AddElement( H, fractionmass=0.080538 );
  UVT_Lucite->AddElement( C, fractionmass=0.599848 );
  UVT_Lucite->AddElement( O, fractionmass=0.319614 );

  MPT_temp = new G4MaterialPropertiesTable();
  //MPT_temp->AddConstProperty("RINDEX", 1.5); //Neglect dispersion in lucite for now, index of refraction is sufficient information.

  const G4int nentries_lucite = 21;
  
  //Lucite is strongly absorptive below 300 nm: transmission data from P. Carter Cal Tech thesis.
  G4double Ephoton_lucite[nentries_lucite] = {1.719888896*eV, 1.81647209*eV, 1.924545279*eV, 2.054623669*eV, 2.293705949*eV, 2.598432583*eV, 2.900093346*eV, 3.145984951*eV, 
					      3.465763144*eV, 3.704958425*eV, 3.973356238*eV, 4.049879001*eV, 4.122678976*eV, 4.184217105*eV, 4.219202221*eV, 4.276408832*eV, 
					      4.320338518*eV, 4.36519548*eV, 4.410978*eV, 4.521654214*eV, 4.521654214*eV};
  G4double Abslength_lucite[nentries_lucite] = {4.359093727*cm, 4.478062096*cm, 4.603526015*cm, 4.603526015*cm, 4.478062096*cm, 4.359093727*cm, 4.035983645*cm, 3.938294782*cm, 
						3.51091739*cm, 3.164087254*cm, 2.727358613*cm, 2.319695564*cm, 1.529808225*cm, 0.873567226*cm, 0.532172776*cm, 0.28887992*cm, 
						0.148639788*cm, 0.093529994*cm, 0.07329661*cm, 0.046324745*cm, 0.046324745*cm };

  G4double Ephoton_rindex_lucite[2] = {1.77*eV, 6.2*eV}; 
  G4double Rindex_lucite[2] = {1.5, 1.5};
  
  MPT_temp->AddProperty( "ABSLENGTH", Ephoton_lucite, Abslength_lucite, nentries_lucite );
  MPT_temp->AddProperty( "RINDEX", Ephoton_rindex_lucite, Rindex_lucite, 2 );

  UVT_Lucite->SetMaterialPropertiesTable( MPT_temp );
  // G4cout << "Material optical properties for material UVT_Lucite" << G4endl;
  // MPT_temp->DumpTable();

  MaterialsMap["UVT_Lucite"] = UVT_Lucite;

  G4Material *Steel = man->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  MaterialsMap["Steel"] = Steel;

  //Need to define "QuartzWindow" *before* Aerogel because the Aerogel Material properties definition refers to Quartz Window:
  
  // For the time being, we treat the mirror as being a single shell of Carbon, with a uniform thickness corresponding 
  // to 0.01 X0:

  G4double den_composite = 2.21*g/cm3;

  G4Material *MirrorComposite = new G4Material( "MirrorComposite", den_composite, nel=1 );
  MirrorComposite->AddElement( C, fractionmass = 1.0 );
  MaterialsMap["MirrorComposite"] = MirrorComposite;

  //Next we will need to define Aluminum for the RICH entry and exit windows and other parts of the containment volume.
  //Let's use the standard NIST database parameters for Al "G4_Al";
  G4Material *RICHAluminum = man->FindOrBuildMaterial("G4_Al");
  MaterialsMap["RICHAluminum"] = RICHAluminum;

  G4double den_pmtwindow = 2.20*g/cm3;
  
  G4Material *QuartzWindow = new G4Material( "QuartzWindow", den_pmtwindow, nel=2 );
  QuartzWindow->AddElement( O, natoms=2 );
  QuartzWindow->AddElement( Si, natoms=1 );

  //Define refractive index and absorption length for quartz PMT windows:
  const G4int nentries_quartz = 51;

  G4double Ephoton_quartz[nentries_quartz] = { 1.77120301*eV, 1.796872619*eV, 1.823297216*eV, 1.850510608*eV, 1.878548647*eV, 1.907449396*eV, 1.937253292*eV, 1.968003345*eV, 
					       1.999745334*eV, 2.032528044*eV, 2.066403512*eV, 2.1014273*eV, 2.137658805*eV, 2.175161591*eV, 2.214003763*eV, 2.254258377*eV, 
					       2.296003902*eV, 2.33932473*eV, 2.384311744*eV, 2.431062955*eV, 2.479684214*eV, 2.530290015*eV, 2.58300439*eV, 2.63796193*eV, 
					       2.695308928*eV, 2.755204682*eV, 2.817822971*eV, 2.883353737*eV, 2.952005017*eV, 3.024005139*eV, 3.099605268*eV, 3.179082326*eV, 
					       3.262742387*eV, 3.350924614*eV, 3.444005853*eV, 3.54240602*eV, 3.646594433*eV, 3.757097294*eV, 3.874506585*eV, 3.999490668*eV, 
					       4.132807024*eV, 4.275317611*eV, 4.428007525*eV, 4.592007804*eV, 4.768623489*eV, 4.959368428*eV, 5.16600878*eV, 5.390617857*eV, 
					       5.635645941*eV, 5.904010034*eV, 6.199210536*eV };

  //Refractive index data from "Refractiveindex.info"
  G4double Rindex_quartz[nentries_quartz] = { 1.455292466, 1.455524071, 1.455763571, 1.456011496, 1.456268423, 1.456534974, 1.456811819, 1.457099689, 
					      1.457399374, 1.457711733, 1.458037702, 1.4583783, 1.458734641, 1.459107942, 1.459499536, 1.459910886, 
					      1.460343603, 1.460799458, 1.461280408, 1.461788618, 1.462326487, 1.462896682, 1.463502175, 1.464146283, 
					      1.464832722, 1.465565665, 1.466349815, 1.467190482, 1.46809369, 1.469066293, 1.470116119, 1.471252144, 
					      1.472484709, 1.473825777, 1.475289258, 1.476891413, 1.478651361, 1.48059172, 1.482739429, 1.485126813, 
					      1.487792976, 1.490785646, 1.494163661, 1.498000361, 1.502388312, 1.507446007, 1.513327606, 1.520237459, 
					      1.528452449, 1.53835762, 1.550505538 };
  
  //Typical absorption length for fused silica from Thorlabs.com (uncoated UV fused silica, includes fresnel reflections)
  G4double Abslength_quartz[nentries_quartz] = {15.65792444*cm, 15.78550788*cm, 15.7794917*cm, 15.60910249*cm, 15.72664954*cm, 15.72488912*cm, 15.57290011*cm, 15.68021339*cm, 
						15.73546266*cm, 15.55685196*cm, 15.55490625*cm, 15.63907251*cm, 15.48113765*cm, 15.54074565*cm, 15.39638598*cm, 15.50169846*cm, 
						15.4950396*cm, 15.36125979*cm, 15.41113687*cm, 15.33874196*cm, 15.24165927*cm, 15.25602267*cm, 15.23330157*cm, 15.14071666*cm, 
						15.13642486*cm, 15.06590584*cm, 15.05023293*cm, 14.99006002*cm, 14.91826095*cm, 14.79500397*cm, 14.80590431*cm, 14.66396966*cm, 
						14.57959363*cm, 14.47194788*cm, 14.40952367*cm, 14.20967861*cm, 14.11981056*cm, 13.98888512*cm, 13.79714319*cm, 13.85187177*cm, 
						13.61931079*cm, 13.27721911*cm, 12.75893298*cm, 12.51276543*cm, 12.13753406*cm, 11.84114847*cm, 10.96932156*cm, 10.13778162*cm, 
						9.92333989*cm, 9.282597031*cm, 7.443241349*cm };

  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("RINDEX", Ephoton_quartz, Rindex_quartz, nentries_quartz);
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_quartz, Abslength_quartz, nentries_quartz);
  
  QuartzWindow->SetMaterialPropertiesTable( MPT_temp );
  MaterialsMap["QuartzWindow"] = QuartzWindow;

  //Next define borosilicate "UV glass" for the PMT window. For chemical composition, use data for Pyrex Corning 7740 borosilicate glass from PDG, but assign same optical properties as quartz above.

  G4double den_UVglass = 2.23*g/cm3;
  
  G4Material *UVglass = new G4Material( "UVglass", den_UVglass, nel=6 );
  UVglass->AddElement( B, fractionmass=0.040061 );
  UVglass->AddElement( O, fractionmass=0.539564 );
  UVglass->AddElement( Na, fractionmass=0.028191 );
  UVglass->AddElement( Al, fractionmass=0.011644 );
  UVglass->AddElement( Si, fractionmass=0.377220 );
  UVglass->AddElement( K,  fractionmass=0.003321 );
  UVglass->SetMaterialPropertiesTable( MPT_temp );

  MaterialsMap["UVglass"] = UVglass;

  G4double den_C4F10 = 10.62*mg/cm3; //Wow, that really is a heavy gas, heavier than HTCC mirror substrates made of Rohacell31!
  //We also need to define the radiator gas C4F10.
  G4Material *C4F10_gas = new G4Material( "C4F10_gas", den_C4F10, nel=2 );
  
  C4F10_gas->AddElement( C, natoms=4 );
  C4F10_gas->AddElement( F, natoms=10 );

  MPT_temp = new G4MaterialPropertiesTable();
  //MPT_temp->AddConstProperty( "RINDEX", 1.00137 ); //assume a constant refractive index for the gas radiator for now. 
  //MPT_temp->AddConstProperty( "ABSLENGTH", 100.0*m );
  const G4int nentries_C4F10 = 2;
  G4double Ephoton_C4F10[nentries_C4F10] = { 1.77*eV, 6.20*eV };
  G4double Rindex_C4F10[nentries_C4F10] = {1.00137, 1.00137};
  G4double Abslength_C4F10[nentries_C4F10] = {100.0*m, 100.0*m};

  MPT_temp->AddProperty("RINDEX", Ephoton_C4F10, Rindex_C4F10, nentries_C4F10 );
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_C4F10, Abslength_C4F10, nentries_C4F10 );

  C4F10_gas->SetMaterialPropertiesTable( MPT_temp );
  MaterialsMap["C4F10_gas"] = C4F10_gas;
  
  //Quantum efficiency for PMT photocathode ("typical", from XP1911/UV data sheet):
  const G4int nentries_QE = 40;
  
  G4double Ephoton_QE[nentries_QE] = { 1.939182049*eV, 1.963579994*eV, 1.999929682*eV, 2.0346944*eV, 2.067643295*eV, 2.098539267*eV, 2.127141474*eV, 2.153230717*eV, 
				       2.179967872*eV, 2.224858112*eV, 2.264313372*eV, 2.293899315*eV, 2.339752453*eV, 2.359406493*eV, 2.403825804*eV, 2.432973578*eV, 
				       2.484622785*eV, 2.538512459*eV, 2.599590307*eV, 2.673824751*eV, 2.768699841*eV, 2.858865976*eV, 2.986530569*eV, 3.105387967*eV, 
				       3.168432345*eV, 3.279402652*eV, 3.349799677*eV, 3.526419346*eV, 3.674006182*eV, 3.953218421*eV, 4.139635948*eV, 4.385211047*eV, 
				       4.570741196*eV, 4.871970386*eV, 5.177284871*eV, 5.437857883*eV, 5.611675974*eV, 5.845211394*eV, 5.969437475*eV, 6.234402273*eV }; 
  
  G4double PMT_QuantumEfficiency[nentries_QE] = { 0.001154539, 0.001948441, 0.003714689, 0.006298763, 0.009646797, 0.013657231, 0.018010571, 0.022819293, 
				     0.028693173, 0.038099138, 0.048532161, 0.057843653, 0.075581257, 0.084283662, 0.105012092, 0.116629941, 
				     0.132736748, 0.149971254, 0.16338945, 0.176043552, 0.193934134, 0.213040691, 0.227782388, 0.229627292, 
				     0.223657738, 0.215914559, 0.213826088, 0.213229177, 0.200888412, 0.17403965, 0.161019419, 0.142756599, 
				     0.126474694, 0.104423377, 0.080171292, 0.057628786, 0.041016469, 0.024094955, 0.012166848, 0.004610016 };

  //Define another material for photocathode: 
  G4double den_photocathode = 2.57*g/cm3;
  G4Material *Photocathode_material = new G4Material( "Photocathode_material", den_photocathode, nel=3 );
  
  Photocathode_material->AddElement( Sb, natoms=1 );
  Photocathode_material->AddElement( K, natoms=2 );
  Photocathode_material->AddElement( Cs, natoms=1 );
  
  G4double Ephot_Rcathode[2] = {1.77*eV, 6.20*eV};
  G4double Rcathode[2] = {0.0, 0.0};

  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("EFFICIENCY", Ephoton_QE, PMT_QuantumEfficiency, nentries_QE );
  MPT_temp->AddProperty("RINDEX", Ephoton_quartz, Rindex_quartz, nentries_quartz );
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_quartz, Abslength_quartz, nentries_quartz );
  //MPT_temp->AddProperty("REFLECTIVITY", Ephot_Rcathode, Rcathode, 2 );

  Photocathode_material->SetMaterialPropertiesTable( MPT_temp );
  
  MaterialsMap["Photocathode_material"] = Photocathode_material;

  // G4cout << "Material Properties for C4F10_gas:" << G4endl;
  // MPT_temp->DumpTable();

  G4double RefractiveIndexAerogel = 1.0304; //At 633 nm; what is the dispersion?
  G4double DensityAerogel = ( (pow(RefractiveIndexAerogel,2)-1.0)/0.424 )*g/cm3; //about 0.146 g/cm^3
  
  G4Material *Aerogel = new G4Material("Aerogel", DensityAerogel, nel=3 );
  
  //These values are from NIM A 433, 396 (1999). The chemical composition does not exactly match that of the 
  //HERMES RICH aerogel, but the density and refractive index do. Should be okay for simulation purposes:
  Aerogel->AddElement( O,  fractionmass = 0.543192 );
  Aerogel->AddElement( Si, fractionmass = 0.453451 );
  Aerogel->AddElement( H,  fractionmass = 0.003357 );
  
  //Next we must define the optical properties of aerogel: most of this information is taken from NIM A 40, 338 (2000):
  // And also Paul Carter's Ph.D. thesis from Cal. Tech. 
  // The properties of interest for our purposes are index of refraction and Rayleigh scattering length.

  //One important property of aerogel that significantly affects the light yield is the relatively short scattering length. Absorption length is 
  //much longer and for our purposes can be neglected. Let's define the scattering length according to the formula:
  // T = A exp(-C t/ lambda^4 ), where T is the transmission. From NIM A 40, 338, we find that the measured average value of A is 0.964
  // and the measured average value of Ct is 0.0094 um^4. 

  G4double A_aerogel = 0.964; //for a thickness of 1 cm. This implies that A = e^{-1 cm/Labs}; ln A = -1 cm / Labs --> Labs = -1 cm / ln(A)
  G4double Ct_aerogel = 0.0094; //microns^4

  G4double t_aerogel = 0.5*(1.125+1.0); //nominal average thickness of aerogel tile in cm
  
  //Since Ct is in units of um^4, let's convert t to microns to compute C in um^3:
  G4double C_aerogel = Ct_aerogel / (t_aerogel*1.0e4); //This quantity has units of um^3!

  //The scattering length at a given lambda is found from L = lambda^4 / C:
  
  G4double lmin = 200.0, lmax = 700.0; //nm
  G4double lstep = 10.0; //nm
  const G4int nsteps = G4int( (lmax-lmin)/lstep + 0.5 ) + 1;

  G4double hbarc_eV_nm = 197.3269718; //eV nm
  
  
  G4double Ephoton_aerogel[nsteps];
  G4double Rindex_aerogel[nsteps];
  G4double Rayleigh_aerogel[nsteps];

  G4bool inrange = true;

  G4double Ephoton_633 = ( twopi * hbarc_eV_nm / 633.0 )*eV;
  G4double nquartz_633 = ( ( QuartzWindow->GetMaterialPropertiesTable() )->GetProperty("RINDEX") )->GetValue( Ephoton_633, inrange );

  for(int i=0; i<nsteps; i++){
    G4double ltemp = lmin + i*lstep; //nm
    Ephoton_aerogel[i] = (twopi * hbarc_eV_nm / ltemp )*eV;
    //Since C is in um^3, we need to compute lambda^4 in um^4. The result will be in microns, which we then convert to cm:
    Rayleigh_aerogel[i] = (pow( ltemp/1000.0, 4 )/C_aerogel / 10000.0 )*cm; //scattering length in cm
    //From the references above, the mean refractive index at 633 nm is 1.0304;
    // Assume that the dispersion follows that of fused silica, so scale n at each wavelength by the ratio of 
    // fused silica at that wavelength to fused silica at 633 nm. 
    G4double nquartz_temp = ( ( QuartzWindow->GetMaterialPropertiesTable() )->GetProperty("RINDEX") )->GetValue( Ephoton_aerogel[i], inrange );
    Rindex_aerogel[i] = 1.0304 * nquartz_temp / nquartz_633;
  }
    
  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("RINDEX", Ephoton_aerogel, Rindex_aerogel, nsteps );
  MPT_temp->AddProperty("RAYLEIGH", Ephoton_aerogel, Rayleigh_aerogel, nsteps );
  //MPT_temp->AddConstProperty("ABSLENGTH", 10.0*m );

  Aerogel->SetMaterialPropertiesTable( MPT_temp );
  MaterialsMap["Aerogel"] = Aerogel;

  //  G4cout << "Material properties for Aerogel:" << G4endl; 
  // MPT_temp->DumpTable();

  //Define the reflectivity of the mirror surface using the "Logical Skin surface":
  const G4int nentries_mirr = 15;
  
  G4double Ephoton_mirr[nentries_mirr] = { 1.787025848*eV, 1.851389327*eV, 2.027907076*eV, 2.206874769*eV, 2.468929557*eV, 
					   2.822981327*eV, 3.15509551*eV, 3.548346967*eV, 3.916325599*eV, 4.348919321*eV, 
					   4.496873938*eV, 4.750862572*eV, 5.176123788*eV, 5.633059855*eV, 6.178512519*eV };  
  
  G4double Reflectivity_mirr[nentries_mirr] = { 0.867162, 0.872027, 0.879324, 0.882973, 0.884189, 0.884189, 0.882973, 
						0.878108, 0.858649, 0.841622, 0.823378, 0.765, 0.687162, 0.619054, 0.557027 };
  
  /// Optical surfaces (for mirrors, etc.):
  //First, define a default optical surface for general use; Aluminum with 100% reflectivity:
  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddConstProperty( "REFLECTIVITY", 1.0 );
  G4OpticalSurface *DefaultOpticalSurface = new G4OpticalSurface("MirrorDefault");
 
  DefaultOpticalSurface->SetType( dielectric_metal );
  DefaultOpticalSurface->SetFinish( polished );
  DefaultOpticalSurface->SetModel( glisur );
  DefaultOpticalSurface->SetMaterialPropertiesTable(MPT_temp);
 
  OpticalSurfacesMap["MirrorDefault"] = DefaultOpticalSurface;
  
  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty( "REFLECTIVITY", Ephoton_mirr, Reflectivity_mirr, nentries_mirr );
  
  G4OpticalSurface *Mirrsurf = new G4OpticalSurface("Mirrsurf");

  Mirrsurf->SetType( dielectric_metal );
  Mirrsurf->SetFinish( polished );
  Mirrsurf->SetModel( glisur );
  
  Mirrsurf->SetMaterialPropertiesTable( MPT_temp );

  

  //G4cout << "Material properties for Mirrsurf:" << G4endl;
  
}

G4Material *G4SBSDetectorConstruction::GetMaterial(G4String name){
  
  if( MaterialsMap.empty() ) ConstructMaterials();

  map<G4String, G4Material*>::iterator it = MaterialsMap.find( name );

  if( it != MaterialsMap.end() ){ 
    return MaterialsMap[name];
  } else {
    it = MaterialsMap.begin();
    return it->second;
  }
}

G4OpticalSurface *G4SBSDetectorConstruction::GetOpticalSurface( G4String name ){
  if( OpticalSurfacesMap.empty() ) ConstructMaterials();
  map<G4String, G4OpticalSurface*>::iterator it = OpticalSurfacesMap.find( name );
  if( it != OpticalSurfacesMap.end() ){
    return OpticalSurfacesMap[name];
  } else {
    it = OpticalSurfacesMap.begin();
    return it->second;
  }
}

G4VPhysicalVolume* G4SBSDetectorConstruction::ConstructAll()
{
  G4cout << "\nG4SBSDetectorConstruction....\n" << G4endl;
 
  if( MaterialsMap.empty() ) ConstructMaterials();
 
  SDman = G4SDManager::GetSDMpointer();
  //--------- Material definition was moved to ConstructMaterials()---------
  //--------- G4VSolid, G4LogicalVolume, G4VPhysicalVolume  ---------
  
  //--------------
  // World:
  //--------------
  G4Box *WorldBox= new G4Box("WorldBox",20*m, 20*m, 20*m);
  G4LogicalVolume *WorldLog=new G4LogicalVolume(WorldBox,GetMaterial("Air"),
						  "WorldLogical", 0, 0, 0);
  G4PVPlacement *WorldPhys=new G4PVPlacement(0,G4ThreeVector(),
					       "WorldPhysical",
					       WorldLog,
					       0,false,0);

  // In the new version of ConstructAll(), now we call individual modular subsystem creation routines depending on fExpType:

  //All three types of experiments have a target:
  ConstructTarget(WorldLog);
  
  //All three types of experiments have a beam line:
  ConstructBeamline(WorldLog);

  // All three types of experiments have a 48D48 magnet:
  //--------- 48D48 -------------------------------
  Make48D48(WorldLog, f48D48dist + 1219.2*mm/2 );

  //--------------- HCAL --------------------------
  //All the experiments use HCAL:

  G4double HCAL_vertical_offset = 0.0*cm; //Neutron/SIDIS experiments have no vertical offset for HCAL (in Neutron case because it is detecting neutrons, which don't bend in a magnetic field, and in SIDIS case because we are detecting +/- charged hadrons simultaneously, want to have symmetric acceptance).
  if( fExpType == kGEp ) HCAL_vertical_offset = 49.7*cm; //A number like this, which represents a positioning offset, shouldn't be hard-coded!

  MakeHCAL( WorldLog, HCAL_vertical_offset );

  
  //  The neutron experiments and the SIDIS experiment use BigBite:
  //------------ BigBite: -----------------------------------------------------
  if( fExpType == kNeutronExp || fExpType == kSIDISExp ) 
    {
      MakeBigBite( WorldLog );
    }
  //The SIDIS experiment uses a  RICH in SBS:
  //--------- RICH (experimental): -------------------------
  if( fExpType == kSIDISExp ) //SIDIS experiment requires a RICH detector and a tracker for SBS: 
    {
      //Let's make a simple tracker: 6 planes of GEMs, equally spaced in z, separation in z between planes of 10 cm. Then total length of tracker is ~60 cm + about 1.6 cm
      G4double SBStracker_dist = fRICHdist - 0.5*m;
      G4ThreeVector SBStracker_pos( SBStracker_dist * sin( f48D48ang ), 0.0, SBStracker_dist * cos( f48D48ang ) );
      
      G4RotationMatrix *SBStracker_rot_I = new G4RotationMatrix(G4RotationMatrix::IDENTITY);
      
      G4RotationMatrix *SBStracker_rot = new G4RotationMatrix;
      SBStracker_rot->rotateY( -f48D48ang );

      G4Box *SBStracker_box = new G4Box("SBStracker_box", 26.0*cm, 102.0*cm, 30.0*cm );
      
      G4LogicalVolume *SBStracker_log = new G4LogicalVolume( SBStracker_box, GetMaterial("Air"), "SBStracker_log" );
      
      new G4PVPlacement( SBStracker_rot, SBStracker_pos, SBStracker_log, "SBStracker_phys", WorldLog, false, 0 );

      int ngems_SBStracker = 6;
      vector<double> zplanes_SBStracker, wplanes_SBStracker, hplanes_SBStracker;
      
      G4double zspacing_SBStracker = 10.0*cm;
      G4double zoffset_SBStracker = -25.0*cm;

      for(int i=0; i<ngems_SBStracker; i++ ){
	zplanes_SBStracker.push_back( zoffset_SBStracker + i*zspacing_SBStracker );
	wplanes_SBStracker.push_back( 50.0*cm );
	hplanes_SBStracker.push_back( 200.0*cm );
      }

      MakeTracker( SBStracker_log, SBStracker_rot_I, G4ThreeVector(0,0,0), 
		   ngems_SBStracker, zplanes_SBStracker, wplanes_SBStracker, hplanes_SBStracker );
      MakeRICH( WorldLog );

      SBStracker_log->SetVisAttributes(G4VisAttributes::Invisible);
    }
  //---------------------------------------------------------
  if( fExpType == kGEp ) //Subsystems unique to the GEp experiment include FPP and BigCal:
    {
      //Let's make a box and then put the FPP in it:
      //Define the rotation matrix for the FPP (pitch angle of 5 deg relative to vertical): 
      G4double sbsboxpitch = 5.0*deg;
      G4RotationMatrix *SBS_FPP_rm = new G4RotationMatrix;
      SBS_FPP_rm->rotateY( -f48D48ang );
      SBS_FPP_rm->rotateX( sbsboxpitch );

      //FPP box: 
      double sbsdepth  = 3.0*m;
      double sbswidth  = 2.0*m;
      double sbsheight = 2.1*m;
      
      double sbsr = fHCALdist-4.106*m + sbsheight*sin(sbsboxpitch)/2+sbsdepth/2;

      G4Box *sbsbox = new G4Box("sbsbox", sbswidth/2.0, sbsheight/2.0, sbsdepth/2.0 );
      G4LogicalVolume* sbslog = new G4LogicalVolume(sbsbox, GetMaterial("Air"), "sbslog");

      sbslog->SetVisAttributes( G4VisAttributes::Invisible );
      //Now position and orient the FPP "box":
      new G4PVPlacement(SBS_FPP_rm, G4ThreeVector(sbsr*sin(f48D48ang), (sbsr-f48D48dist)*sin(sbsboxpitch), sbsr*cos(f48D48ang) ), sbslog,
	      "sbsphys", WorldLog, false, 0, false);

      G4RotationMatrix *rot_I = new G4RotationMatrix;
      
      double detoffset = 0.05*m - sbsdepth/2.0;

      MakeFPP( sbslog, rot_I, G4ThreeVector( 0.0, 0.0, detoffset) );
      MakeBigCal( WorldLog );
    }
  

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
  

  

  G4VisAttributes * dVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  dVisAtt->SetForceWireframe(true);
  //dLog->SetVisAttributes(dVisAtt);

  G4VisAttributes * xVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  xVisAtt->SetForceWireframe(true);
  //xLog->SetVisAttributes(xVisAtt);

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

void G4SBSDetectorConstruction::ConstructTarget( G4LogicalVolume *worldlog ){

  //Material definition was moved to ConstructMaterials();
  //--------- Glass target cell -------------------------------

  double wallthick = 1.61*mm;
  double capthick  = 0.126*mm;
  double cellradius    = 0.75*2.54*cm/2.0;

  G4Tubs *targ_tube = new G4Tubs("targ_tube", cellradius-wallthick, cellradius, fTargLen/2.0, 0.*deg, 360.*deg );
  G4Tubs *targ_cap = new G4Tubs("targ_cap", 0.0, cellradius, capthick/2.0, 0.*deg, 360.*deg );

  G4LogicalVolume* targ_tube_log = new G4LogicalVolume(targ_tube, GetMaterial("GE180"),"targ_tube_log");
  G4LogicalVolume* targ_cap_log = new G4LogicalVolume(targ_cap, GetMaterial("GE180"),"targ_cap_log");

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
    gas_tube_log = new G4LogicalVolume(gas_tube, GetMaterial("refH2"), "gas_tube_log");
  }
  if( fTargType == k3He ){
    gas_tube_log = new G4LogicalVolume(gas_tube, GetMaterial("pol3He"), "gas_tube_log");
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

  targ_tube_log = new G4LogicalVolume(cryoshell, GetMaterial("Aluminum"),"targ_tube_log");

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
  G4LogicalVolume *extpipe_log = new G4LogicalVolume(exttube, GetMaterial("Aluminum"),"extpipe_log");
  G4LogicalVolume *extvac_log = new G4LogicalVolume(extvactube, GetMaterial("Vacuum"),"extvac_log");

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, extpipestart-extpipe_len/2), extpipe_log, "extpipe_phys", worldlog, false, 0);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, extpipestart-extpipe_len/2), extvac_log, "extvacpipe_phys", worldlog, false, 0);

  if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
      // Add in exit Al window
      
      double extwin_thick = 5.0e-4*cm;
      
      G4Tubs *extwin = new G4Tubs("ent_win", 0.0, extpipe_rin, extwin_thick/2, 0.*deg, 360.*deg );
      G4LogicalVolume *ext_winlog = new G4LogicalVolume(extwin, GetMaterial("Aluminum"), "entwin_log", 0, 0, 0);
      new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -extpipe_len/2 + extwin_thick/2), ext_winlog, "extwin_phys", extvac_log,false,0);

      ext_winlog->SetVisAttributes(new G4VisAttributes(G4Colour(0.6,0.6,0.6)));
  }


  //  Place exit pipe tube

  G4Tubs *swall_enthole = new G4Tubs("scham_wall_enthole", 0.0, entpipe_rin, 20.0*cm, 0.*deg, 360.*deg );
  G4Tubs *swall_exthole = new G4Tubs("scham_wall_exthole", 0.0, extpipe_rin, 20.0*cm, 0.*deg, 360.*deg );

  G4RotationMatrix *chamholerot = new G4RotationMatrix;
  chamholerot->rotateY(90.0*deg);

  //  Cut holes in the scattering chamber
  G4SubtractionSolid* swall_holes = new G4SubtractionSolid("swall_enthole", swallcut, swall_enthole, chamholerot, G4ThreeVector(-(swallrad+swallrad_in)/2, 0.0, 0.0) );
  swall_holes = new G4SubtractionSolid("swall_holes", swall_holes, swall_exthole, chamholerot, G4ThreeVector((swallrad+swallrad_in)/2, 0.0, 0.0) );

  G4LogicalVolume *swall_log = new G4LogicalVolume(swall_holes, GetMaterial("Aluminum"),"scham_wall_log");

  G4LogicalVolume *sc_hcalwin_log = new G4LogicalVolume(swall_hcalwin, GetMaterial("Aluminum"),"sc_hcalwin_log");
  G4LogicalVolume *sc_bbwin_log = new G4LogicalVolume(swall_bbwin, GetMaterial("Aluminum"),"sc_bbwin_log");

  G4RotationMatrix *schamrot = new G4RotationMatrix;
  schamrot->rotateX(-90.0*deg);
  schamrot->rotateZ(-90.0*deg);

  G4RotationMatrix *targrot = new G4RotationMatrix;
  targrot->rotateY(-90.0*deg);

  G4Tubs *chamber_inner = new G4Tubs("chamber_inner", 0.0, swallrad_in,  sheight/2, 0*deg, 360*deg );
  G4LogicalVolume* chamber_inner_log = new G4LogicalVolume(chamber_inner, GetMaterial("Vacuum"), "cham_inner_log");

  // Top and bottom
  G4Tubs *sc_topbottom = new G4Tubs("scham_topbottom", 0.0, swallrad, (swallrad-swallrad_in)/2, 0.*deg, 360.*deg );
  G4LogicalVolume* sc_topbottom_log = new G4LogicalVolume(sc_topbottom, GetMaterial("Aluminum"), "scham_topbottom_log");

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
    cryo_tube_log = new G4LogicalVolume(cryovol, GetMaterial("LH2mat"), "cryo_tube_log");
  }
  if( fTargType == kLD2 ){
    cryo_tube_log = new G4LogicalVolume(cryovol, GetMaterial("LD2mat"), "cryo_tube_log");
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

  //Material definition moved to "ConstructMaterials":
  
  
  //G4Material* aluminum = new G4Material("Aluminum", 13., 26.98*g/mole, density=2.7*g/cm3);
  // G4Material *aluminum = GetMaterial("Aluminum");
  // G4Material* vacuum = GetMaterial("Vacuum");

    double swallrad = 1.143*m/2;
    double beamheight = 10.0*12*2.54*cm; // 10 feet off the ground

    // Stainless
    G4double ent_len = 10*m;
    G4double ent_rin = 31.75*mm;
    G4double ent_rou = ent_rin+0.120*mm;

    G4Tubs *ent_tube = new G4Tubs("ent_tube", ent_rin, ent_rou, ent_len/2, 0.*deg, 360.*deg );
    G4Tubs *ent_vac  = new G4Tubs("ent_vac", 0.0, ent_rin, ent_len/2, 0.*deg, 360.*deg );

    G4LogicalVolume *entLog = new G4LogicalVolume(ent_tube, GetMaterial("Stainless"), "ent_log", 0, 0, 0);
    G4LogicalVolume *entvacLog = new G4LogicalVolume(ent_vac, GetMaterial("Vacuum"), "entvac_log", 0, 0, 0);


    if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
	// gas target -  1.5m in air
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-1.5*m), entLog, "ent_phys", worldlog, false,0);
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-1.5*m), entvacLog, "entvac_phys", worldlog,false,0);

	// Add in Be window
	G4double winthick = 0.0127*cm;

	G4Tubs *ent_win = new G4Tubs("ent_win", 0.0, ent_rin, winthick/2, 0.*deg, 360.*deg );
	G4LogicalVolume *ent_winlog = new G4LogicalVolume(ent_win, GetMaterial("Beryllium"), "entwin_log", 0, 0, 0);
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, ent_len/2-winthick/2), ent_winlog, "entwin_phys", entvacLog,false,0);
	ent_winlog->SetVisAttributes(new G4VisAttributes(G4Colour(0.7,1.0,0.0)));
    } else {
	// Cryotarget - up against the chamber wall
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad), entLog, "ent_phys", worldlog, false,0);
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad), entvacLog, "entvac_phys", worldlog,false,0);
    }

    // Aluminum
    int nsec = 24;
    //  Definition taken from HAPLOG 2722 by Juliette, but offset by 31.54 cm
    G4double exit_z[]   = {206*cm, 234.01*cm, 234.02*cm, 253.02*cm, 253.03*cm, 268.26*cm, 268.27*cm,305.29*cm, 305.30*cm,328.71*cm, 328.72*cm, 356.33*cm,356.34*cm, 378.7*cm,378.71*cm, 473.16*cm,473.17*cm, 503.64*cm,503.65*cm, 609.84*cm,609.85*cm, 1161.02*cm, 1161.03*cm,2725.66*cm };
    G4double exit_zero[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    G4double exit_rin[] = {4.128*cm, 4.128*cm, 4.445*cm, 4.445*cm,4.763*cm, 4.763*cm, 5.08*cm,5.08*cm, 6.35*cm, 6.35*cm, 7.62*cm, 7.62*cm,10.16*cm, 10.16*cm,10.478*cm, 10.478*cm,12.7*cm, 12.7*cm,15.24*cm,15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };
    G4double exit_rou[] = {4.432*cm, 4.432*cm, 4.75*cm, 4.75*cm,5.067*cm, 5.067*cm, 5.385*cm,5.385*cm, 6.655*cm, 6.655*cm, 7.925*cm, 7.925*cm, 10.478*cm,10.478*cm,  10.795*cm, 10.795*cm, 13.018*cm, 13.018*cm,15.558*cm, 15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };

    G4Polycone *ext_cone = new G4Polycone("ext_tube", 0.0*deg, 360.0*deg, nsec, exit_z, exit_rin, exit_rou);
    G4Polycone *ext_vac  = new G4Polycone("ext_vac ", 0.0*deg, 360.0*deg, nsec, exit_z, exit_zero, exit_rin);

    G4LogicalVolume *extLog = new G4LogicalVolume(ext_cone, GetMaterial("Aluminum"), "ext_log", 0, 0, 0);
    G4LogicalVolume *extvacLog = new G4LogicalVolume(ext_vac, GetMaterial("Vacuum"), "extvac_log", 0, 0, 0);

    new G4PVPlacement(0,G4ThreeVector(), extLog, "ext_phys", worldlog, false,0);
    new G4PVPlacement(0,G4ThreeVector(), extvacLog, "extvac_phys", worldlog,false,0);


    double floorthick = 1.0*m;
    G4Tubs *floor_tube = new G4Tubs("floor_tube", 0.0, 30*m, floorthick/2, 0.*deg, 360.*deg );

    G4RotationMatrix *floorrm = new G4RotationMatrix;
    floorrm->rotateX(90*deg);

    G4LogicalVolume *floorLog = new G4LogicalVolume(floor_tube, GetMaterial("Concrete"), "floor_log", 0, 0, 0);
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

  //moved definition of clamp gaps to field clamp method:
  // G4Box *bclampgap  = new G4Box("bclampgap",  23.*cm, 65.*cm,  12.*cm/2.);
  // G4Box *fclampgap  = new G4Box("fclampgap",  11.*cm, 35.*cm,  12.*cm/2.);

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

  G4Box *bigbeamslot = new G4Box("bigbeamslot",  bigwidth/2, 15.5*cm/2.0, 2.0*m ); // Height is roughly beam pipe outer radius at 3m
 
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

  G4LogicalVolume *big48d48Log=new G4LogicalVolume(big48d48_wslot, GetMaterial("Fer"),
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

  G4LogicalVolume *bigfieldLog=new G4LogicalVolume(biggap, GetMaterial("Air"),
						  "bigfieldLog", 0, 0, 0);

  // use uniform field for now with 48D48

  double fieldValue = f48D48_uniform_bfield;
  G4UniformMagField* magField
            = new G4UniformMagField(G4ThreeVector(fieldValue*cos(f48D48ang), 0.0, -fieldValue*sin(f48D48ang)));

  G4FieldManager *bigfm = new G4FieldManager(magField);

  bigfm->SetDetectorField(magField);
  bigfm->CreateChordFinder(magField);

  bigfieldLog->SetFieldManager(bigfm,true);


  new G4PVPlacement(bigrm, 
//	  G4ThreeVector(-(f48D48dist+bigdepth/2.0)*sin(-f48D48ang), 0.0, (f48D48dist+bigdepth/2.0)*cos(-f48D48ang)),
	  G4ThreeVector(-r48d48*sin(-f48D48ang), 0.0, r48d48*cos(-f48D48ang)),
	  		    bigfieldLog, "bigfieldPhysical", worldlog, 0,false,0);


  
  // Clamps
  
  // The positioning and acceptance gaps were taken directly from CAD
  // The position widths for the beam pipe holes are fudged around so
  // that it doesn't interfere with the beam pipe
  //
  if( f48D48_fieldclamp_config != 0 ){

    MakeSBSFieldClamps(worldlog);

  }

  bigfieldLog->SetVisAttributes(G4VisAttributes::Invisible);
//  backclampLog->SetVisAttributes(G4VisAttributes::Invisible);
//  frontclampLog->SetVisAttributes(G4VisAttributes::Invisible);
}

void G4SBSDetectorConstruction::MakeSBSFieldClamps( G4LogicalVolume *motherlog ){
  double clampdepth = 10.*cm;
  double clampoffset = 35*cm;
  double bigheight = 3721.1*mm;

  G4Box *bclampgap  = new G4Box("bclampgap",  23.*cm, 65.*cm,  12.*cm/2.);
  G4Box *fclampgap  = new G4Box("fclampgap",  11.*cm, 35.*cm,  12.*cm/2.);

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

  G4LogicalVolume *frontclampLog=new G4LogicalVolume(frontclamp, GetMaterial("Fer"), "frontclampLog", 0, 0, 0);
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

  G4LogicalVolume *backclampLog=new G4LogicalVolume(backclamp, GetMaterial("Fer"), "backclampLog", 0, 0, 0);

  if( fTotalAbs ){
      backclampLog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }

  double frontclampz = -100*cm + clampdepth/2.0;
  double backclampz  =  100*cm - clampdepth/2.0;
  
  G4RotationMatrix *rot = new G4RotationMatrix;
  rot->rotateY( -f48D48ang );

  double r48d48 = f48D48dist + 1219.2*mm/2.0;

  new G4PVPlacement(rot, 
		    //	  G4ThreeVector(-(f48D48dist+bigdepth/2.0+frontclampz)*sin(-f48D48ang)-cos(-f48D48ang)*clampoffset, 0.0, (f48D48dist+bigdepth/2.0+frontclampz)*cos(-f48D48ang)-sin(-f48D48ang)*clampoffset),
		    G4ThreeVector(-(r48d48+frontclampz)*sin(-f48D48ang)-cos(-f48D48ang)*clampoffset, 0.0, (r48d48+frontclampz)*cos(-f48D48ang)-sin(-f48D48ang)*clampoffset),
		    frontclampLog, "frontclampPhysical", motherlog, 0,false,0);
  new G4PVPlacement(rot, 
		    //	  G4ThreeVector(-(f48D48dist+bigdepth/2.0+backclampz)*sin(-f48D48ang)-cos(-f48D48ang)*clampoffset, 0.0, (f48D48dist+bigdepth/2.0+backclampz)*cos(-f48D48ang)-sin(-f48D48ang)*clampoffset),
		    G4ThreeVector(-(r48d48+backclampz)*sin(-f48D48ang)-cos(-f48D48ang)*clampoffset, 0.0, (r48d48+backclampz)*cos(-f48D48ang)-sin(-f48D48ang)*clampoffset),
		    backclampLog, "backclampPhysical", motherlog, 0,false,0);
}

void G4SBSDetectorConstruction::MakeBigBite( G4LogicalVolume *motherlog ){

  //G4SDManager* SDman = G4SDManager::GetSDMpointer();

  //Lines of code used to build BigBite moved to their own method:
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

  G4LogicalVolume *bbyokewgapLog=new G4LogicalVolume(yokewgap, GetMaterial("Fer"),
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
  G4LogicalVolume *bbmotherLog=new G4LogicalVolume(bbmothercutBox,GetMaterial("Air"),
						  "bbmotherLog", 0, 0, 0);

  new G4PVPlacement(bbrm, G4ThreeVector(-motherr*sin(fBBang), 0.0, motherr*cos(fBBang)),
	  bbmotherLog, "bbmotherPhys", motherlog, 0,false,0);


  new G4PVPlacement(yokerm,G4ThreeVector(0.0, 0.0, -motherdepth/2.0+clear),
	  bbyokewgapLog, "bbyokewgapPhysical", bbmotherLog, 0,false,0);

  //  Bigbite field log volume
  G4LogicalVolume *bbfieldLog=new G4LogicalVolume(bbairTrap, GetMaterial("Air"),
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
  G4LogicalVolume *bbdetLog=new G4LogicalVolume(bbdetbox, GetMaterial("Air"),
						"bbdetLog", 0, 0, 0);
  new G4PVPlacement(bbdetrot, G4ThreeVector(0.0, 
					    (detboxplace+detboxdepth/2.0)*sin(detboxang),
					    (detboxplace+detboxdepth/2.0)*cos(detboxang)+midplanez),
		    bbdetLog, "bbdetPhysical", bbmotherLog, 0,false,0);
  
  //  Just interested in the GEMs for now:

  double detoffset = 0.05*m -detboxdepth/2.0; //z offset of GEM plane positions within BB detector volume: "global" GEM plane z = detoffset + gemz[plane]

  int i;
  int ngem = 0;
  double gemdsep;
  // double *gemz;
  // double *gemw;
  // double *gemh;
  vector<double> gemz, gemw, gemh;

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

  gemz.resize(ngem);
  gemw.resize(ngem);
  gemh.resize(ngem);
  //
  // GEM option 1
  double gemz_opt1[] = { 0.0*cm, gemdsep, fGEMDist, fGEMDist+gemdsep};
  double gemw_opt1[] = { 40.0*cm, 40.0*cm, 50.0*cm, 50.0*cm };
  double gemh_opt1[] = { 150.0*cm, 150.0*cm, 200.0*cm, 200.0*cm };

  // GEM option 2
  gemdsep = 0.15*m; //Is this intended?
  double gemz_opt2[] = { 0.0*cm, gemdsep, 2.0*gemdsep, 3.0*gemdsep, fGEMDist};
  double gemw_opt2[] = { 40.0*cm, 40.0*cm, 40.0*cm, 40.0*cm, 50.0*cm };
  double gemh_opt2[] = { 150.0*cm, 150.0*cm, 150.0*cm, 150.0*cm, 200.0*cm };
  
  // GEM option 3
  double gemz_opt3[] = { 0.0*cm, gemdsep, gemdsep*2.0};
  double gemw_opt3[] = { 40.0*cm, 50.0*cm, 50.0*cm};
  double gemh_opt3[] = { 150.0*cm, 200.0*cm, 200.0*cm};

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

  G4RotationMatrix *rot_identity = new G4RotationMatrix;

  //This routine creates and positions GEM planes in bbdetLog:

  //------------------------------------------- BigBite GEMs: ----------------------------------------//
  MakeTracker( bbdetLog, rot_identity, G4ThreeVector( 0.0, 0.0, detoffset ), ngem, gemz, gemw, gemh );
  //----- Note: Lines of code that are common to the construction of all individual GEM planes/modules were moved to MakeTracker() -----// 
  //----- All we do here in MakeBigBite() is define the number of planes, their z positions, and their transverse dimensions ------//

  // BigBite Preshower 
  // AJP: Why make a preshower box if it's just going to be full of air and not a sensitive detector?

  double psheight = 27*8.5*cm;
  double pswidth  = 2.0*37.0*cm;
  double psdepth  = 8.5*cm;

  G4Box *bbpsbox = new G4Box("bbpsbox", pswidth/2.0, psheight/2.0, psdepth/2.0 );

  G4LogicalVolume* bbpslog = new G4LogicalVolume(bbpsbox, GetMaterial("Air"), "bbpslog");
  
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth/2.0), bbpslog,
	      "bbpsphys", bbdetLog, false, 0, false);

  // BigBite Shower
  
  double calheight = 27*8.5*cm;
  double calwidth  = 7*8.5*cm;
  double caldepth  = 37.0*cm;

  G4Box *bbcalbox = new G4Box("bbcalbox", calwidth/2.0, calheight/2.0, caldepth/2.0 );
  G4LogicalVolume* bbcallog = new G4LogicalVolume(bbcalbox, GetMaterial("Lead"), "bbcallog");
  
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
  
  G4LogicalVolume* cer_winlog_in = new G4LogicalVolume(cer_winbox_in, GetMaterial("Aluminum"), "cer_winlog_in");
  G4LogicalVolume* cer_winlog_out = new G4LogicalVolume(cer_winbox_out, GetMaterial("Aluminum"), "cer_winlog_out");
  //  G4LogicalVolume* cer_mirlog = new G4LogicalVolume(cer_mirbox, SiO2, "cer_mirlog");
  G4LogicalVolume* cer_mirlog = new G4LogicalVolume(cer_mirbox, GetMaterial("Acrylic"), "cer_mirlog");
  G4LogicalVolume* cer_gaslog = new G4LogicalVolume(cer_gasbox, GetMaterial("C4F8O"), "cer_gaslog");
  
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
  
  //--------- Visualization attributes -------------------------------
  bbdetLog->SetVisAttributes(G4VisAttributes::Invisible);
  bbfieldLog->SetVisAttributes(G4VisAttributes::Invisible);
  bbmotherLog->SetVisAttributes(G4VisAttributes::Invisible);
  
  G4VisAttributes * yokeVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  //  yokeVisAtt->SetForceWireframe(true);
  bbyokewgapLog->SetVisAttributes(yokeVisAtt);
  
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

}

void G4SBSDetectorConstruction::MakeHCAL( G4LogicalVolume *motherlog, G4double VerticalOffset=0.0*cm ){
  //G4SDManager* SDman = G4SDManager::GetSDMpointer();

  // HCAL 
  double hcalheight = 330.0*cm;
  double hcalwidth  = 165.0*cm;
  double hcaldepth  = 101.0*cm;
  double hcalr = fHCALdist+hcaldepth/2.0;

  G4RotationMatrix *hcalrm = new G4RotationMatrix;
  hcalrm->rotateY(-f48D48ang);

  G4Box *hcalbox = new G4Box("hcalbox", hcalwidth/2.0, hcalheight/2.0, hcaldepth/2.0 );
  G4LogicalVolume* hcallog = new G4LogicalVolume(hcalbox, GetMaterial("Lead"), "hcallog");
     
  new G4PVPlacement(hcalrm, G4ThreeVector(hcalr*sin(f48D48ang), VerticalOffset, hcalr*cos(f48D48ang) ), hcallog,
	      "hcalphys", motherlog, false, 0, false);

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

  G4VisAttributes * hcalVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.6,0.0));

  hcallog->SetVisAttributes(hcalVisAtt);
}

void G4SBSDetectorConstruction::MakeRICH( G4LogicalVolume *motherlog ){

  //*********************************************************************************************************************************//
  //                                  BEGIN GEOMETRY DEFINITION FOR SBS RICH COUNTER                                                 //
  //*********************************************************************************************************************************//
  
  
  //First, define a RICH box that will serve as the mother volume for the positioning of the RICH components relative to each other and 
  //as the containment volume for the C4F10 radiator gas:
  G4double RICHbox_dx=165.0*cm, RICHbox_dy=283.0*cm, RICHbox_dz=127.0*cm;

  G4Box *RICHbox = new G4Box( "RICHbox", RICHbox_dx/2.0, RICHbox_dy/2.0, RICHbox_dz/2.0 );
  G4LogicalVolume *RICHbox_log = new G4LogicalVolume( RICHbox, GetMaterial("C4F10_gas"), "RICHbox_log" );
  
  //We also want to define exterior walls of the RICH box as 1-inch thick aluminum: we will do a subtraction solid of RICHbox and RICHwalls 
  //as well as the entry windows...
  G4double RICHwall_dx = RICHbox_dx - 2.0*2.54*cm;
  G4double RICHwall_dy = RICHbox_dy - 2.0*2.54*cm;
  G4double RICHwall_dz = RICHbox_dz - 2.0*2.54*cm;
  G4Box *RICHwall = new G4Box( "RICHwall", RICHwall_dx/2.0, RICHwall_dy/2.0, RICHwall_dz/2.0 );

  G4SubtractionSolid *RICHbox_wall = new G4SubtractionSolid("RICHbox_wall", RICHbox, RICHwall );

  // The stacking of the tiles is defined by three parameters, nx, ny, and nz, the number of tiles along each dimension. x is assumed to be the horizontal direction, y the 
  // vertical direction, and z the nominal direction of particle motion. The design constraint is that nx*ny*nz <= 850:
  G4int nx_aero=5, ny_aero=17, nz_aero=5;

  //define tile dimensions: 
  G4double aero_dx = (11.4/2.0)*cm, aero_dy=(11.4/2.0)*cm, aero_dz=(1.13/2.0)*cm;

  G4Box *Aerogel_tile_solid = new G4Box("Aerogel_tile", aero_dx, aero_dy, aero_dz );  
  G4LogicalVolume *Aerogel_tile_log = new G4LogicalVolume( Aerogel_tile_solid, GetMaterial("Aerogel"), "Aerogel_tile_log" );

  //Assume 1 mil gap between tiles, filled with tedlar

  //Decide where we want center of aerogel coordinate system to be located. 
  //Assume that the geometric center of the combined aerogel box is at the origin of RICHbox
  
  G4double tilegap = 0.00254*cm; //gap between tiles in x and y:

  G4double Lx_aero = 2.0*nx_aero*aero_dx + tilegap*(nx_aero-1);
  G4double Ly_aero = 2.0*ny_aero*aero_dy + tilegap*(ny_aero-1);
  G4double Lz_aero = 2.0*nz_aero*aero_dz;

  G4double RICH_entrywindow_thick = 0.1*cm;
  G4double RICH_entrywindow_dz = RICH_entrywindow_thick/2.0;
  G4double RICH_entrywindow_dx = Lx_aero/2.0;
  G4double RICH_entrywindow_dy = Ly_aero/2.0;

  //1 mm-thick aluminum entry window for RICH, assumed to have same dimensions as aerogel tiles. 
  G4Box *RICH_entrywindow = new G4Box("RICH_entrywindow", RICH_entrywindow_dx, RICH_entrywindow_dy, RICH_entrywindow_dz );
  G4LogicalVolume *RICH_entrywindow_log = new G4LogicalVolume( RICH_entrywindow, GetMaterial("RICHAluminum"), "RICH_entrywindow_log" );
  
  G4double Aero_exit_dx = Lx_aero/2.0;
  G4double Aero_exit_dy = Ly_aero/2.0;
  G4double Aero_exit_dz = (0.32/2.0)*cm;
  //3.2 mm-thick UVT-lucite exit window for aerogel:
  G4Box *Aero_exit = new G4Box( "Aero_exit", Aero_exit_dx, Aero_exit_dy, Aero_exit_dz );
  G4LogicalVolume *Aero_exitwindow = new G4LogicalVolume( Aero_exit, GetMaterial("UVT_Lucite"), "Aero_exitwindow" );

  //1 mm-thick aluminum exit window for RICH, dimensions are up to us, but let's start with HERMES case:
  G4double RICH_exitwindow_thick = 0.1*cm;
  G4double RICH_exitwindow_dx = (59.0/2.0)*cm;
  G4double RICH_exitwindow_dy = (257.0/2.0)*cm;
  G4double RICH_exitwindow_dz = RICH_exitwindow_thick/2.0;

  G4Box *RICH_exit = new G4Box( "RICH_exit", RICH_exitwindow_dx, RICH_exitwindow_dy, RICH_exitwindow_dz );
  G4LogicalVolume *RICH_exitwindow = new G4LogicalVolume( RICH_exit, GetMaterial("RICHAluminum"), "RICH_exitwindow" );

  G4double aero_xoffset = 7.191*cm;

  G4double z0_entrywindow = -RICHbox_dz/2.0 + RICH_entrywindow_dz;
  G4double y0_entrywindow = 0.0;
  G4double x0_entrywindow = -RICHbox_dx/2.0 + aero_xoffset + Lx_aero/2.0;

  //Position entry and exit windows inside RICHbox:
  G4VPhysicalVolume *RICH_entrywindow_pv = new G4PVPlacement( 0, 
							      G4ThreeVector( x0_entrywindow, y0_entrywindow, z0_entrywindow), 
							      RICH_entrywindow_log, 
							      "RICH_entrywindow_pv", 
							      RICHbox_log, 
							      false, 
							      0 );

							      
  G4double x0_aeroexit = x0_entrywindow;
  G4double y0_aeroexit = 0.0;
  G4double z0_aeroexit = z0_entrywindow + RICH_entrywindow_dz + Lz_aero + Aero_exit_dz;

  //This is aerogel exit window.
  G4VPhysicalVolume *Aero_exitwindow_pv = new G4PVPlacement( 0, 
							     G4ThreeVector( x0_aeroexit, y0_aeroexit, z0_aeroexit ),
							     Aero_exitwindow, 
							     "Aero_exitwindow_pv",
							     RICHbox_log, 
							     false,
							     0 );

  G4double x0_RICHexit = -RICHbox_dx/2.0 + aero_xoffset + RICH_exitwindow_dx;
  G4double y0_RICHexit = 0.0;
  G4double z0_RICHexit = RICHbox_dz/2.0 - RICH_exitwindow_dz;

  G4VPhysicalVolume *RICH_exitwindow_pv = new G4PVPlacement( 0, 
							     G4ThreeVector( x0_RICHexit, y0_RICHexit, z0_RICHexit ),
							     RICH_exitwindow, 
							     "RICH_exitwindow_pv",
							     RICHbox_log, 
							     false, 
							     0 );

  //We need to define cutouts from the 1"-thick Aluminum box for the entry and exit windows:
  G4Box *RICHentry_cutout = new G4Box( "RICHentry_cutout", RICH_entrywindow_dx, RICH_entrywindow_dy, 10.0*cm );
  G4SubtractionSolid *RICHbox_wall_entrycut = new G4SubtractionSolid( "RICHbox_wall_entrycut", RICHbox_wall, RICHentry_cutout, 0, G4ThreeVector( x0_entrywindow, y0_entrywindow, -RICHwall_dz/2.0 ) );

  G4Box *RICHexit_cutout = new G4Box( "RICHexit_cutout", RICH_exitwindow_dx, RICH_exitwindow_dy, 10.0*cm );
  G4SubtractionSolid *RICHbox_wall_entryexitcut = new G4SubtractionSolid( "RICHbox_wall_entryexitcut", RICHbox_wall_entrycut, RICHexit_cutout, 0, G4ThreeVector( x0_RICHexit, y0_RICHexit, RICHwall_dz/2.0 ) );

  G4LogicalVolume *RICH_container_walls = new G4LogicalVolume( RICHbox_wall_entryexitcut, GetMaterial("RICHAluminum"), "RICH_container_walls" );
  
  G4VPhysicalVolume *RICH_container_walls_placement = new G4PVPlacement( 0, G4ThreeVector(0,0,0), RICH_container_walls, "RICH_container_walls_placement", RICHbox_log, false, 0 );

  //We also need to define the Tedlar spacers: Let us define horizontal and vertical spacers. Our convention will be that the 
  //vertical spacer fills the corner region:
  G4double horizontal_spacer_dx = aero_dx, horizontal_spacer_dy = tilegap/2.0, horizontal_spacer_dz=Lz_aero/2.0;
  G4double vertical_spacer_dx = tilegap/2.0, vertical_spacer_dy = Ly_aero/2.0, vertical_spacer_dz=Lz_aero/2.0;

  G4Box *Horizontal_spacer = new G4Box( "Horizontal_spacer", horizontal_spacer_dx, horizontal_spacer_dy, horizontal_spacer_dz );
  G4Box *Vertical_spacer = new G4Box( "Vertical_spacer", vertical_spacer_dx, vertical_spacer_dy, vertical_spacer_dz );
  
  G4LogicalVolume *Horizontal_spacer_log = new G4LogicalVolume( Horizontal_spacer, GetMaterial("Tedlar"), "Horizontal_spacer_log" );
  G4LogicalVolume *Vertical_spacer_log = new G4LogicalVolume( Vertical_spacer, GetMaterial("Tedlar"), "Vertical_spacer_log" );

  G4int icopy = 1;
  G4int icopy_vertical_spacer = 1;
  G4int icopy_horizontal_spacer = 1;
 
  G4String tilename;

  //Next: position aerogel tiles: 
  for(G4int ix=0; ix<nx_aero; ix++ ){
    
    G4double xtemp = -RICHbox_dx/2.0 + aero_xoffset + (ix+0.5)*2.0*aero_dx + ix*tilegap;
    
    for(G4int iy=0; iy<ny_aero; iy++ ){
    
      G4double ytemp = -Ly_aero/2.0 + (iy+0.5)*2.0*aero_dy + iy*tilegap;
      
      for(G4int iz=0; iz<nz_aero; iz++ ){
	//compute center coordinates of aerogel tiles and position in RICHbox:
	
	G4double ztemp = -RICHbox_dz/2.0 + RICH_entrywindow_thick + (iz+0.5)*2.0*aero_dz;
	
	tilename = "Aerogel_tile_pv_";
	
	char ccopy[20];
	
	sprintf( ccopy, "%d", icopy );

	tilename += ccopy;

	new G4PVPlacement( 0, 
			   G4ThreeVector( xtemp, ytemp, ztemp ), 
			   Aerogel_tile_log, 
			   tilename,
			   RICHbox_log,
			   false,
			   icopy++ );
	
      }
      if( iy > 0 && iy+1 < ny_aero ){
	G4double xspacer = xtemp;
	G4double yspacer = ytemp + aero_dy + tilegap/2.0;
	G4double zspacer = -RICHbox_dz/2.0 + RICH_entrywindow_thick + Lz_aero/2.0;

	new G4PVPlacement( 0, 
			   G4ThreeVector( xspacer, yspacer, zspacer ), 
			   Horizontal_spacer_log, 
			   "Horizontal_spacer_pv",
			   RICHbox_log, 
			   false, 
			   icopy_horizontal_spacer++ );
	
      }
    } 
    //position vertical tedlar spacers:
    if( ix>0 && ix+1 < nx_aero ){
      //vertical spacer position is equal to tile position + half tile width + half gap width:
      G4double xspacer = xtemp + aero_dx + tilegap/2.0;
      G4double yspacer = 0.0;
      G4double zspacer = -RICHbox_dz/2.0 + RICH_entrywindow_thick + Lz_aero/2.0;
      
      new G4PVPlacement( 0, 
			 G4ThreeVector( xspacer, yspacer, zspacer ),
			 Vertical_spacer_log, 
			 "Vertical_spacer_pv",
			 RICHbox_log, 
			 false,
			 icopy_vertical_spacer++ );

    }
  }
  
  //Next, let's try to define the mirror. For this, we can probably use a "spherical shell section" without resorting to 
  //solid operations. Alternatively we could use polycone.
 
  G4double MirrorCenter_x = -RICHbox_dx/2.0 + 136.403*cm;
  G4double MirrorCenter_y = 0.0*cm;
  G4double MirrorCenter_z = -RICHbox_dz/2.0 - 103.112*cm;

  G4double MirrorShell_thick = 0.1932*cm; //This corresponds to 1% X0 of graphite.

  G4double MirrorRadius = 220.0*cm;
  //G4double MirrorRadius_XZproject = 180.190*cm;

  //This is relative to RICHbox
  G4double Mirror_xmin = -RICHbox_dx/2.0 + 5.540*cm;
  G4double Mirror_xmax = 79.505*cm + Mirror_xmin;

  G4double Mirror_zmin = MirrorCenter_z + 123.890*cm;
  G4double Mirror_zmax = RICHbox_dz/2.0;

  G4double Mirror_ymin = -126.2*cm;
  G4double Mirror_ymax = 126.2*cm;

  //Let's make a spherical shell and a "cut box"
  
  //The spherical shell should go from 0 to 90 degrees in theta and +/- 90 degrees in phi, and then we form the intersection 
  // with the cut box defining the minimum and maximum z planes:
  //Include the entire forward hemisphere for simplicity:
  G4Sphere *RICH_mirror_shell = new G4Sphere( "RICH_mirror_shell", MirrorRadius, MirrorRadius + MirrorShell_thick, 0.0, twopi, 0.0, halfpi );
  G4Box *RICH_mirror_cutbox = new G4Box( "RICH_mirror_cutbox", (Mirror_xmax-Mirror_xmin)/2.0, (Mirror_ymax-Mirror_ymin)/2.0, (Mirror_zmax-Mirror_zmin)/2.0 );

  //Now, we want to make the intersection of the two solids, so we need to express the coordinates of the center of the box in the mirror
  // coordinate system. 
  G4ThreeVector MirrorCenterCoords( MirrorCenter_x, MirrorCenter_y, MirrorCenter_z );
  G4ThreeVector BoxCenterCoords( 0.5*(Mirror_xmin+Mirror_xmax), 0.0, 0.5*(Mirror_zmin+Mirror_zmax) );
  G4ThreeVector RelativeCoords = BoxCenterCoords - MirrorCenterCoords;
  
  G4IntersectionSolid *Mirror_solid = new G4IntersectionSolid( "Mirror_solid", RICH_mirror_shell, RICH_mirror_cutbox, 0, RelativeCoords );
  G4LogicalVolume *Mirror_log = new G4LogicalVolume( Mirror_solid, GetMaterial("MirrorComposite"), "Mirror_log" );

  G4VPhysicalVolume *Mirror_pv = new G4PVPlacement( 0, 
						    MirrorCenterCoords, 
						    Mirror_log, 
						    "Mirror_pv",
						    RICHbox_log, 
						    false,
						    0 );

  G4LogicalSkinSurface *Mirrskin = new G4LogicalSkinSurface( "Mirrskin", Mirror_log, GetOpticalSurface("Mirrsurf") );

  //What is left? We've done the mirror, the aerogel, the gas box. All that remains is the PMTs and the structure of the containment vessel. Let's start with the PMTs: 

  ////////////////////////////////////////////////////////////////////////
  //                         !!!PMTS!!!                                 //
  ////////////////////////////////////////////////////////////////////////


  //Define the PMT windows as 1 mm-thick discs of "UVglass":
  G4Tubs *PMTwindow = new G4Tubs( "PMTwindow", 0.0*cm, (1.86/2.0)*cm, 0.05*cm, 0.0, twopi ); 
  //Define the PMT photocathode as a thin disc of 15 mm
  G4Tubs *PMTcathode = new G4Tubs( "PMTcathode", 0.0*cm, (1.50/2.0)*cm, 0.025*cm, 0.0, twopi );
  G4Tubs *PMTtube    = new G4Tubs( "PMTtube", (1.66/2.0)*cm, (1.86/2.0)*cm, (8.7/2.0)*cm, 0.0, twopi );

  //"Quartz window" is a different, sealed window that separates the PMT from the C4F10 environment.
  G4Tubs *PMTQuartzWindow = new G4Tubs( "PMTQuartzWindow", 0.0*cm, (2.13/2.0)*cm, 0.15*cm, 0.0, twopi );
  //CollectionCone is a light-collecting cone that increases the effective collection efficiency:
  G4Cons *CollectionCone = new G4Cons( "CollectionCone", 0.75*cm, (2.13/2.0)*cm, (2.13/2.0)*cm, (2.13/2.0)*cm, 0.75*cm, 0.0, twopi );

  G4double PMT_total_length = 10.65*cm;
  G4double PMT_max_radius = 1.065*cm;

  G4LogicalVolume *PMTwindow_log  = new G4LogicalVolume( PMTwindow, GetMaterial("UVglass"), "PMTwindow_log" );
  G4LogicalVolume *PMTcathode_log = new G4LogicalVolume( PMTcathode, GetMaterial("Photocathode_material"), "PMTcathode_log" );

  //PMTcathode_log is the sensitive detector for the RICH:

  //  G4SDManager *SDman = G4SDManager::GetSDMpointer();

  G4String RICHSDname = "G4SBS/RICH";
  G4String RICHcollname = "RICHcoll";
  G4SBSRICHSD *RICHSD = NULL;

  if( !( (G4SBSRICHSD*) SDman->FindSensitiveDetector(RICHSDname) ) ){
    G4cout << "Adding RICH sensitive detector to SDman..." << G4endl;
    RICHSD = new G4SBSRICHSD( RICHSDname, RICHcollname );
    SDman->AddNewDetector( RICHSD );
  }
  //}

  PMTcathode_log->SetSensitiveDetector( RICHSD ); //This assigns the sensitive detector type "RICHSD" to the logical volume PMTcathode!
  //We make this a hollow cylinder with length and radius approximately equal to that of the PMT housing, made of steel 
  //to approximate the material shielding the PMT.
  G4LogicalVolume *PMTtube_log    = new G4LogicalVolume( PMTtube, GetMaterial("Steel"), "PMTtube_log" ); 
  G4LogicalVolume *PMTquartzwindow_log = new G4LogicalVolume( PMTQuartzWindow, GetMaterial("QuartzWindow"), "PMTQuartzWindow_log" );
  G4LogicalVolume *CollectionCone_log = new G4LogicalVolume( CollectionCone, GetMaterial("Steel"), "CollectionCone_log" );
  //Define a logical skin surface for the collection cone and assign it the same reflectivity as the mirror:
  G4LogicalSkinSurface *Coneskin = new G4LogicalSkinSurface( "Coneskin", CollectionCone_log, GetOpticalSurface("Mirrsurf") );

  //Within the RICHbox, each PMT assembly unit is rotated so that its z-axis makes an angle of 50 degrees with respect to the 
  //local z axis of the RICHbox. Therefore, we rotate by an angle of 
  G4double PMT_rotation_angle = 50.0*degree;
  G4RotationMatrix *rot_PMT = new G4RotationMatrix;
  rot_PMT->rotateY( PMT_rotation_angle );

  G4int icopy_PMT_assembly = 0;
  
  G4double xfp = 119.350*cm - RICHbox_dx/2.0;
  G4double yfp = 0.0;
  G4double zfp = 42.521*cm - RICHbox_dz/2.0;

  G4ThreeVector focalpoint_position( xfp, yfp, zfp );
  
  G4ThreeVector PMT_zaxis( -sin(PMT_rotation_angle), 0.0, cos(PMT_rotation_angle) );
  G4ThreeVector PMT_yaxis( 0, 1, 0 );
  G4ThreeVector PMT_xaxis( (PMT_yaxis.cross( PMT_zaxis ) ).unit() );
  
  G4double ymin_PMT = -72.5376*cm, ymax_PMT = 72.5376*cm;
  G4double xmin_PMT[2] = { -29.083*cm, -30.24632*cm };
  G4double xmax_PMT[2] = { 29.083*cm, 30.24632*cm };
  G4int nrows_PMT[2] = {26, 27};

  for( G4int icol=0; icol<=72; icol++){
    G4int evenoddcol = icol%2;
    for( G4int irow=0; irow<nrows_PMT[evenoddcol]; irow++ ){
      G4double xtemp = xmin_PMT[evenoddcol] + irow * ( xmax_PMT[evenoddcol] - xmin_PMT[evenoddcol] )/( G4double(nrows_PMT[evenoddcol]-1) );
      G4double ytemp = ymin_PMT + icol*(ymax_PMT-ymin_PMT)/( 72.0 );

      G4ThreeVector PMT_position = focalpoint_position - PMT_zaxis * PMT_total_length/2.0 + xtemp * PMT_xaxis + ytemp * PMT_yaxis;
      
      //Place PMT components inside RICHbox.
      G4ThreeVector Pos_temp;
      //Steel tube (mainly for visualization and shielding
      G4double ztube = -PMT_total_length/2.0 + 4.35*cm;
      Pos_temp = PMT_position + ztube * PMT_zaxis;
      G4VPhysicalVolume *PMTtube_pv = new G4PVPlacement( rot_PMT, Pos_temp, PMTtube_log, "PMTtube_pv", RICHbox_log, false, icopy_PMT_assembly );
      //Photocathode (this is the sensitive part!!):
      G4double zcathode = ztube + 4.375*cm;
      Pos_temp = PMT_position + zcathode * PMT_zaxis;
      G4VPhysicalVolume *PMTcathode_pv = new G4PVPlacement( rot_PMT, Pos_temp, PMTcathode_log, "PMTcathode_pv", RICHbox_log, false, icopy_PMT_assembly );
      //UV-glass PMT window:
      G4double zwindow = zcathode + 0.075*cm;
      Pos_temp = PMT_position + zwindow * PMT_zaxis;
      G4VPhysicalVolume *PMTwindow_pv = new G4PVPlacement( rot_PMT, Pos_temp, PMTwindow_log, "PMTwindow_pv", RICHbox_log, false, icopy_PMT_assembly );
      //Quartz window between PMT and gas:
      G4double zquartz = zwindow + 0.2*cm;
      Pos_temp = PMT_position + zquartz * PMT_zaxis;
      G4VPhysicalVolume *PMTquartzwindow_pv = new G4PVPlacement( rot_PMT, Pos_temp, PMTquartzwindow_log, "PMTquartzwindow_pv", RICHbox_log, false, icopy_PMT_assembly );
      //Light collection cone:
      G4double zcone = zquartz + 0.9*cm;
      Pos_temp = PMT_position + zcone * PMT_zaxis;
      G4VPhysicalVolume *CollectionCone_pv = new G4PVPlacement( rot_PMT, Pos_temp, CollectionCone_log, "CollectionCone_pv", RICHbox_log, false, icopy_PMT_assembly );

      // G4VPhysicalVolume *PMT_placement = new G4PVPlacement( rot_PMT, 
      // 							   PMT_position, 
      // 							   PMT_assembly, 
      // 							   "PMT_placement", 
      // 							   RICHbox_log, 
      // 							   false, 
      // 							   icopy_PMT_assembly++ );
      icopy_PMT_assembly++;

      
    }
  }
      

  ////////////////////////////////////////////////////////////////////////
  //                         !!!END OF PMTS!!!                                 //
  ////////////////////////////////////////////////////////////////////////

  //At the end we have to define the translation and rotation to apply for correct positioning of RICHbox:
  //Rotation is easy, same as HCAL, we rotate about the Y axis by -f48D48ang:

  G4RotationMatrix *rot_RICH = new G4RotationMatrix;
  rot_RICH->rotateY( -f48D48ang );

  //We want the center of the RICH entry window to be located at a distance equal to fRICHdist along the line at angle f48D48ang from the origin. For this condition to be satisfied, the center of the RICH box must be offset from this line:
  G4ThreeVector RICHcoord_global( fRICHdist*sin( f48D48ang ), 0.0, fRICHdist*cos( f48D48ang ) );

  G4ThreeVector RICH_zaxis( RICHcoord_global.unit() );
  G4ThreeVector RICH_yaxis( 0.0, 1.0, 0.0 );
  G4ThreeVector RICH_xaxis( (RICH_yaxis.cross( RICH_zaxis )).unit() );

  //RICH center coordinates

  G4ThreeVector RICH_centercoord_global = RICHcoord_global - x0_entrywindow * RICH_xaxis - y0_entrywindow * RICH_yaxis - z0_entrywindow * RICH_zaxis;

  G4VPhysicalVolume *RICHbox_pv = new G4PVPlacement( rot_RICH, 
						     RICH_centercoord_global,
						     RICHbox_log,
						     "RICHbox_pv",
						     motherlog, 
						     false, 
						     0 );


  G4VisAttributes *RICHbox_vis = new G4VisAttributes( G4Colour( 1.0, 1.0, 1.0 ) );
  RICHbox_vis->SetForceWireframe(true);

  RICHbox_log->SetVisAttributes( RICHbox_vis ); 
  
  //Set color and transparency for RICH windows (Aluminum)
  G4VisAttributes *RICHwindow_visatt = new G4VisAttributes( G4Colour( 0.75,0.75,0.75) );
  
  //RICH entry and exit windows are not inherently interesting, so we force them to wireframe:
  RICHwindow_visatt->SetForceWireframe(true);

  RICH_exitwindow->SetVisAttributes( RICHwindow_visatt );
  RICH_entrywindow_log->SetVisAttributes( RICHwindow_visatt );
  
  //Set aerogel exit window to a magenta color (equal parts red and blue) and also wireframe:
  G4VisAttributes *Lucitewindow_visatt = new G4VisAttributes( G4Colour( 1.0,0.0,1.0 ) );
  Lucitewindow_visatt->SetForceWireframe(true); 

  Aero_exitwindow->SetVisAttributes( Lucitewindow_visatt );
  
  G4VisAttributes *aero_tile_visatt = new G4VisAttributes( G4Colour( 0.0, 0.8, 0.8 ) );
  Aerogel_tile_log->SetVisAttributes( aero_tile_visatt );

  G4VisAttributes *tedlar_vis = new G4VisAttributes( G4Colour(0.3,0.3,0.3) );
  Vertical_spacer_log->SetVisAttributes( tedlar_vis );
  Horizontal_spacer_log->SetVisAttributes( tedlar_vis );

  G4VisAttributes *mirror_vis = new G4VisAttributes( G4Colour(0.0, 0.5, 1.0) );
  Mirror_log->SetVisAttributes( mirror_vis );

  //PMT_assembly->SetVisAttributes( G4VisAttributes::GetInvisible() );

  //  G4VisAttributes for PMT assemblies:
  
  G4VisAttributes *PMTtube_vis = new G4VisAttributes( G4Colour( 0.4, 0.4, 0.4 ) );
  PMTtube_vis->SetForceLineSegmentsPerCircle( 24 );
  PMTtube_log->SetVisAttributes( PMTtube_vis );

  G4VisAttributes *PMTwindow_vis = new G4VisAttributes( G4Colour::Cyan() );
  PMTwindow_vis->SetForceLineSegmentsPerCircle( 24 );
  PMTwindow_vis->SetForceWireframe( true );
  PMTwindow_log->SetVisAttributes( PMTwindow_vis );

  //  PMTcathode_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  G4VisAttributes *PMTcathode_vis = new G4VisAttributes( G4Colour::Blue() );
  PMTcathode_vis->SetForceLineSegmentsPerCircle( 24 );
  PMTcathode_log->SetVisAttributes( PMTcathode_vis );
  
  G4VisAttributes *PMTquartzwindow_vis = new G4VisAttributes( G4Colour::Green() );
  PMTquartzwindow_vis->SetForceLineSegmentsPerCircle( 24 );
  PMTquartzwindow_vis->SetForceWireframe( true );
  PMTquartzwindow_log->SetVisAttributes( PMTquartzwindow_vis );

  G4VisAttributes *CollectionCone_vis = new G4VisAttributes( G4Colour::Red() );
  CollectionCone_vis->SetForceLineSegmentsPerCircle( 24 );
  CollectionCone_log->SetVisAttributes( CollectionCone_vis );
  
  G4VisAttributes *RICHwalls_vis = new G4VisAttributes( G4Colour::Gray() );
  RICHwalls_vis->SetForceWireframe( true );
  RICH_container_walls->SetVisAttributes( RICHwalls_vis );
  

}

//This routine allows us to flexibly position GEM modules without code duplication:
void G4SBSDetectorConstruction::MakeTracker(G4LogicalVolume *Mother, G4RotationMatrix *rot, G4ThreeVector pos, G4int nplanes, vector<double> zplanes, vector<double> wplanes, vector<double> hplanes) 
{
  //This routine will create and position a GEM tracker consisting of nplanes planes centered at position pos oriented with rotation rot wrt logical volume Mother. 
  //The list of z coordinates, widths and heights of the planes are passed as arguments:
  
  G4String MotherName = Mother->GetName();

  //  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  
  //Define sensitive detector, if not already done:
  G4String GEMSDname = "G4SBS/GEM";
  G4String GEMcolname = "GEMcol";
  G4SBSGEMSD* GEMSD;

  if( !(GEMSD = (G4SBSGEMSD*) SDman->FindSensitiveDetector(GEMSDname)) ){
      GEMSD = new G4SBSGEMSD( GEMSDname, GEMcolname );
      SDman->AddNewDetector(GEMSD);
  }

  if( !( nplanes > 0 && zplanes.size() == nplanes && wplanes.size() == nplanes && hplanes.size() == nplanes ) ){
    G4cout << "MakeTracker called with invalid arguments, doing nothing..." << G4endl;
    return;
  }

  //Define z extent of various GEM layer components:

  const int nlayers = 24;

  double gempz[nlayers] = {
    120.0*um, 3.0*mm, 120.0*um, 5.0*um, 50.0*um, 3.0*mm, // cover + cathode
    5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM0
    5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM1
    5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM2
    10.0*um, 50.0*um, 180.0*um, 120.0*um, 3.0*mm, 120.0*um // readout + honeycomb
  };
  
  //Define GEM materials:
  
  G4Material *gempm[nlayers] = {
    GetMaterial("NEMAG10"), GetMaterial("NOMEX"), GetMaterial("NEMAG10"), GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("GEMgas"),
    GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("Copper"), GetMaterial("GEMgas"),
    GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("Copper"), GetMaterial("GEMgas"),
    GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("Copper"), GetMaterial("GEMgas"),
    GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("NEMAG10"), GetMaterial("NEMAG10"), GetMaterial("NOMEX"), GetMaterial("NEMAG10")
  };
  
  int gidx, gpidx;
  
  //Determine the largest transverse dimensions of a plane in this tracker:
  double gemmaxw = 0.0;
  double gemmaxh = 0.0;

  for( gidx = 0; gidx < nplanes; gidx++ ){
    if( wplanes[gidx] > gemmaxw ) gemmaxw = wplanes[gidx];
    if( hplanes[gidx] > gemmaxh ) gemmaxh = hplanes[gidx];
  }

  //Determine the total z extent of a plane in this tracker:
  double gempzsum = 0.0;
  for( gpidx = 0; gpidx < nlayers; gpidx++ ){
      gempzsum += gempz[gpidx];
  }
  
  // char gpname[20][50][3][255];
  // char gemname[10][255], gemboxname[10][255], gemlogname[10][255];
  char cgidx[255];
  char cgpidx[255];

  G4Box *gpbox;
  G4LogicalVolume *gplog;

  G4VisAttributes *gemvisatt = new G4VisAttributes( G4Colour( 1.0, 0.0, 0.0 ) );
  
  gemvisatt->SetForceWireframe(true);

  G4VisAttributes *gemsdvisatt = new G4VisAttributes( G4Colour( 1.0, 0.0, 0.0 ) );
  
  for( gidx = 0; gidx < nplanes; gidx++ ){
    sprintf( cgidx, "%02d", gidx );

    // sprintf(gemboxname[gidx], "gembox_%02d", gidx);
    // sprintf(gemlogname[gidx], "gemlog_%02d", gidx);
    G4String gemboxname = MotherName + G4String("_gembox_");
    gemboxname += cgidx;
    G4String gemlogname = MotherName + G4String("_gemlog_");
    gemlogname += cgidx;
    
    G4Box *gembox = new G4Box( gemboxname, wplanes[gidx]/2.0, hplanes[gidx]/2.0, gempzsum/2.0 );
    G4LogicalVolume *gemlog = new G4LogicalVolume( gembox, GetMaterial("Air"), gemlogname, 0, 0, 0 );

    gemlog->SetVisAttributes( gemvisatt );

    double ztemp = 0.0;
    for( gpidx = 0; gpidx < nlayers; gpidx++ ){
      sprintf( cgpidx, "_%02d_%03d_", gidx, gpidx );
      // sprintf(gpname[gidx][gpidx][0], "gemplane_%02d_%03d_box", gidx, gpidx );
      // sprintf(gpname[gidx][gpidx][1], "gemplane_%02d_%03d_log", gidx, gpidx );
      // sprintf(gpname[gidx][gpidx][2], "gemplane_%02d_%03d_phy", gidx, gpidx );
      G4String gempboxname = MotherName + G4String("_gemplane") + cgpidx + G4String("box");
      G4String gemplogname = MotherName + G4String("_gemplane") + cgpidx + G4String("log");
      G4String gempphysname = MotherName + G4String("_gemplane") + cgpidx + G4String("phy");
      
      ztemp += gempz[gpidx]/2.0;
      
      gpbox = new G4Box( gempboxname, wplanes[gidx]/2.0, hplanes[gidx]/2.0, gempz[gpidx]/2.0 );
      gplog = new G4LogicalVolume( gpbox, gempm[gpidx], gemplogname, 0, 0, 0 );
      
      //gplog->SetVisAttributes( 

      new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, ztemp - gempzsum/2.0 ), gplog, gempphysname, gemlog, 0, 0, false ); //??
      
      ztemp += gempz[gpidx]/2.0;

      //Assign sensitive volume: why 5?
      if( gpidx == 5 ){
	gplog->SetSensitiveDetector(GEMSD);
	gplog->SetVisAttributes( gemsdvisatt );
      } else {
	gplog->SetVisAttributes( G4VisAttributes::Invisible );
      }
    }
    
    G4String gemname = MotherName + G4String("_gemphys_") + cgidx;
    //Now place the fully constructed GEM plane in the mother logical volume:
    G4ThreeVector plane_pos = pos + G4ThreeVector( 0.0, 0.0, zplanes[gidx] );
    new G4PVPlacement( rot, plane_pos, gemlog, gemname, Mother, true, gidx+1, false );
  }
}

void G4SBSDetectorConstruction::MakeBigCal( G4LogicalVolume *Mother ){

  //G4SDManager* SDman = G4SDManager::GetSDMpointer();

  //-----------------------------
  //  BigCal: Currently just a box with dimensions and sensitivity
  
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
  G4LogicalVolume* bigcallog = new G4LogicalVolume(bigcalbox, GetMaterial("Lead"), "bigcallog");

  G4Box *CH2box = new G4Box("ch2box", bigcalwidth/2.0, bigcalheight/2.0, CH2depth/2.0 );
  G4LogicalVolume* ch2boxlog = new G4LogicalVolume(CH2box, GetMaterial("CH2"), "ch2log");
  G4Box *CHbox = new G4Box("chbox", bigcalwidth/2.0, bigcalheight/2.0, CHdepth/2.0 );
  G4LogicalVolume* chboxlog = new G4LogicalVolume(CHbox, GetMaterial("CH"), "chlog");

  double ch2r = bbr - bigcaldepth/2.0 - CHdepth - CH2depth/2.0;
  double chr = bbr - bigcaldepth/2.0 - CHdepth/2.0;

  new G4PVPlacement(bbrm, G4ThreeVector(bbr*sin(-fBBang), 0.0, bbr*cos(-fBBang) ), bigcallog,
	      "bigcalphys", Mother, false, 0, false);
  new G4PVPlacement(bbrm, G4ThreeVector(ch2r*sin(-fBBang), 0.0, ch2r*cos(-fBBang) ), ch2boxlog,
	      "ch2boxphys", Mother, false, 0, false);
  new G4PVPlacement(bbrm, G4ThreeVector(chr*sin(-fBBang), 0.0, chr*cos(-fBBang) ), chboxlog,
	      "chboxphys", Mother, false, 0, false);
  
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

}

void G4SBSDetectorConstruction::MakeFPP( G4LogicalVolume *Mother, G4RotationMatrix *rot, G4ThreeVector pos ){
  //FPP consists of GEM tracker interspersed w/CH2 analyzers:

  int i;
  int ngem = 0;

  vector<double> gemz, gemw, gemh;

  ngem = 14; 
  gemz.resize(ngem);
  gemw.resize(ngem);
  gemh.resize(ngem);
  
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

  MakeTracker( Mother, rot, pos, ngem, gemz, gemw, gemh );

  //CH2 analyzers:
  
  double anaheight = 200.0*cm;
  double anawidth  = 44.0*2.54*cm;
  double anadepth  = 22.0*2.54*cm;
  
  G4Box *anabox = new G4Box("anabox", anawidth/2.0, anaheight/2.0, anadepth/2.0 );
  G4LogicalVolume* analog = new G4LogicalVolume(anabox, GetMaterial("CH2"), "analog");

  G4ThreeVector Ana1_pos = pos + G4ThreeVector( 0.0, 0.0, 58.53*cm + anadepth/2.0 );
  G4ThreeVector Ana2_pos = pos + G4ThreeVector( 0.0, 0.0, 170.3*cm + anadepth/2.0 );
  
  new G4PVPlacement(0, Ana1_pos, analog,
		    "anaphys1", Mother, false, 0, false);
  new G4PVPlacement(0, Ana2_pos, analog,
		    "anaphys1", Mother, false, 0, false);
  
}
