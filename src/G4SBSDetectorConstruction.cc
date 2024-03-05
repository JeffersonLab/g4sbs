#include "G4SBSDetectorConstruction.hh"

#include "G4PhysicalConstants.hh"
#include "G4UserLimits.hh"
//#include "G4String.hh"

#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4MaterialPropertiesTable.hh" 
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Element.hh"
#include "G4ProductionCuts.hh"
#include "G4ElementTable.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RotationMatrix.hh"

#include "G4MagneticField.hh"
#include "G4SBSGlobalField.hh"
#include "G4SBSBigBiteField.hh"
#include "G4SBSToscaField.hh"
#include "G4SBSConstantField.hh"
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

#include "G4SBSBeamlineBuilder.hh"
#include "G4SBSTargetBuilder.hh"
#include "G4SBSEArmBuilder.hh"
#include "G4SBSHArmBuilder.hh"
#include "G4SBSECal.hh"
#include "G4SBSCDet.hh"

//#include "RTPC.hh"

#include "G4Mag_SpinEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4DormandPrince745.hh"

#include "G4SBSCalSD.hh"
#include "G4SBSGEMSD.hh"
#include "G4SBSECalSD.hh"

#include "TSpline.h"

#include <vector>
#include <map>
#include <algorithm>
//#include <pair>

// To supress errors with TString, system of units should be included last
#include "G4SystemOfUnits.hh"

using namespace std;

G4SBSDetectorConstruction::G4SBSDetectorConstruction()
{
  f48D48_uniform_bfield = 1.4*tesla;

  fFieldScale_SBS = 1.0;
  fFieldScale_BB = 1.0;
  
  fbbfield = NULL;
  f48d48field = NULL;
  //fIO = NULL;

  fTotalAbs = false;

  fExpType = G4SBS::kGMN;

  // By default, don't check detectors for overlaps
  fCheckOverlap = false;

  //flags controlling ECAL thermal annealing model:
  fSegmentC16 = 0; //default to no segmentation!
  fSegmentThickC16 = 4.0*cm; //default thickness of 4 cm for lead-glass longitudinal segmentation.
  fDoseRateC16 = 0.0; //Default radiation dose rate of ZERO (no rad. damage!)
  
  //ConstructMaterials(); //Now we want to construct all materials at the beginning, so that the physics tables can get built properly (do this in "Construct()", not here)!!!

  fTargetBuilder   = new G4SBSTargetBuilder(this);
  fBeamlineBuilder = new G4SBSBeamlineBuilder(this);
  fEArmBuilder     = new G4SBSEArmBuilder(this);
  fHArmBuilder     = new G4SBSHArmBuilder(this);
  fECal            = new G4SBSECal(this);
  fCDet            = new G4SBSCDet(this);
  //fRTPC = 0;
  fmTPC = new G4SBSmTPC(this);
  
  fHArmBuilder->fFieldStrength = f48D48_uniform_bfield;

  fGlobalField = new G4SBSGlobalField();

  fUseGlobalField = false;

  fBeamlineConf = 3;
  fLeadOption = 1;
  fBLneutronDet = false;
    
  SDlist.clear();
  SDtype.clear();
  StepLimiterList.clear();

  fCDetOption = 1;

  fGEMflip = false;
  //    TrackerIDnumber = 0;
  //TrackerArm.clear();

  //mtpc glabal variables default values
  fmTPCHeGasFraction = 0.7;//fraction of 1
  fmTPCCH4GasFraction = 0.3;//fraction of 1
  fmTPCGasTemp = 296.15;//K
  fmTPCGasPressure = 1.0;//atm

  // D. Flay (8/25/20) 
  // beam diffuser switch 
  fBeamDumpEnable            = true;
  fBeamDiffuserEnable        = false;
  // taret collimators 
  fGEnTgtCollimatorEnable    = true;  
  fGEnTgtCollimatorAEnable   = true;  
  fGEnTgtCollimatorBEnable   = true;  
  fGEnTgtCollimatorCEnable   = true;  
  
  // D. Flay (9/29/20) 
  // GEn 3He target angular misalignment
  fGEnTgtDRX = 0.;  
  fGEnTgtDRY = 0.;  
  fGEnTgtDRZ = 0.;

  // D. Flay (12/9/20) 
  // GEn 3He target as a sensitive detector 
  fGEnTgtSDEnable = false;  

  // D. Flay (10/15/20) 
  // Ion chamber (testing) 
  fIonChamberEnable = false; 
  fIonChamberX      = 0;  
  fIonChamberY      = 0;  
  fIonChamberZ      = 0;  
  fIonChamberRX     = 0;  
  fIonChamberRY     = 0;  
  fIonChamberRZ     = 0;  

  // D. Flay (11/5/20) 
  // beam collimator (testing) 
  fBeamCollimatorEnable_upstr = false; 
  fBeamCollimatorX_upstr      = 0;  
  fBeamCollimatorY_upstr      = 0;  
  fBeamCollimatorZ_upstr      = -145.0*cm;  
  fBeamCollimatorL_upstr      = 4.0*cm;  
  fBeamCollimatorDmin_upstr   = 15*mm;  
  fBeamCollimatorDmax_upstr   = 30*mm;  

  fBeamCollimatorEnable_dnstr = false; 
  fBeamCollimatorX_dnstr      = 0;  
  fBeamCollimatorY_dnstr      = 0;  
  fBeamCollimatorZ_dnstr      = -60.0*cm;  
  fBeamCollimatorL_dnstr      = 4.0*cm;  
  fBeamCollimatorDmin_dnstr   = 15*mm;  
  fBeamCollimatorDmax_dnstr   = 30*mm;  

  //Default Aluminum shielding around GEMs to false:
  fGEMuseAlshield = false;
  fGEMAlShieldThick = 50.0*um;
  fGEMAirGapThick = 0.0*um;
  
}

G4SBSDetectorConstruction::~G4SBSDetectorConstruction()
{;}

G4VPhysicalVolume* G4SBSDetectorConstruction::Construct(){
  // Just nothing so we don't step on toes
  // ConstructAll is where the real magic happens
  //G4double a, iz, z, density;
  //Moved all material definitions to ConstructMaterials()

  if( fMaterialsMap.empty() ) ConstructMaterials();

  G4Material *Mtemp = GetMaterial("BlandAir");

  G4Box *WorldBox= new G4Box("WorldBox",50*m, 50*m, 50*m);
  //G4Box *WorldBox= new G4Box("WorldBox",28*m, 28*m, 28*m);
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

  G4double GC_Gas_Pressure= 1.0*atmosphere;

  const G4int nEntries=25;
  G4double PhotonWaveLength[nEntries]={
    650*nm, 600*nm, 550*nm, 500*nm, 450*nm, 400*nm, 350*nm,
    350*nm, 330*nm, 310*nm, 290*nm, 270*nm, 250*nm, 230*nm, 210*nm,
    210*nm, 208*nm, 206*nm, 204*nm, 202*nm, 200*nm,
    199*nm, 195*nm, 190*nm, 185*nm
  };
  G4double nm_lambda;
  G4double PhotonEnergy[nEntries];
  G4double Std_RefractiveIndex[nEntries];
  for ( int i = 0; i < nEntries; ++i ) {
    PhotonEnergy[i]=(1240*nm/PhotonWaveLength[i])*eV;
    Std_RefractiveIndex[i]=1.00;
  }

  G4MaterialPropertiesTable* Std_MPT = new G4MaterialPropertiesTable();
  Std_MPT->AddProperty("RINDEX", PhotonEnergy, Std_RefractiveIndex, nEntries);



  map<G4String,G4Material*>::iterator itest;
  //pair<G4String, G4Material>  ptest;

  //Let's Fix all these element definitions to use NIST standardized versions:
  
  // G4Element *elN = new G4Element("Nitrogen", "N", 7, 14.007*g/mole );
  // G4Element *elO = new G4Element("Oxygen", "O", 8, 16.000*g/mole );

  G4Element *elN = man->FindOrBuildElement("N");
  G4Element *elO = man->FindOrBuildElement("O");

  //"BlandAir" == vacuum
  density = 1e-9*mg/cm3;
  G4Material* BlandAir = new G4Material(name="BlandAir", density, nel=2);
  BlandAir->AddElement(elN, .7);
  BlandAir->AddElement(elO, .3);

  fMaterialsMap["BlandAir"] = BlandAir;

  G4double a, iz, z, abundance;
  G4int Z, N, A, ncomponents;
  
  //Seamus's element definitions:
  // a = 126.9*g/mole;
  // G4Element* elI = new G4Element(name="Iodine",   symbol="I",  iz=53., a);
  // a = 132.9*g/mole;
  // G4Element* elCs= new G4Element(name="Cesium",   symbol="Cs", iz=55., a);

  G4Element *elI = man->FindOrBuildElement("I");
  G4Element *elCs = man->FindOrBuildElement("Cs");
  
  // G4Element *elH = new G4Element("Hydrogen", "H", 1, 1.007*g/mole );
  // G4Element *elD = new G4Element("Deuterium", "D", 1, 2.014*g/mole );
  // G4Element *el3He = new G4Element("Helium3", "3He", 2, 3.016*g/mole );
  // G4Element *el4He = new G4Element("Helium4", "4He", 2, 4.0026*g/mole );
  // G4Element *elC = new G4Element("Carbon", "C", 6, 12.011*g/mole ) ;
  // G4Element *elF = new G4Element("Fluorine", "F", 9, 18.998*g/mole );
  // G4Element *elNa = new G4Element("Sodium", "Na", 11, 22.99*g/mole );
  // G4Element *elAl = new G4Element("Aluminum", "Al", 13, 26.982*g/mole );
  // G4Element *elSi = new G4Element("Silicon", "Si", 14, 28.086*g/mole );
  // G4Element *elCa = new G4Element("Calcium", "Ca", 20, 40.078*g/mole );
  // G4Element *elFe = new G4Element("Iron", "Fe", 26, 55.850*g/mole );
  // G4Element *elSr = new G4Element("Strontium", "Sr", 38, 87.62*g/mole );
  // G4Element *elBa = new G4Element("Barium", "Ba", 56, 137.327*g/mole );
  //NOTE: "N" in "G4Isotope" refers to the number of nucleons, NOT the number of neutrons, as is the usual notation in nuclear physics
  //We don't need to supply the molar mass, because G4NistManager does it for us!
  G4Isotope *iso_2H = new G4Isotope( "H2", Z=1, N=2 );
  G4Isotope *iso_3H = new G4Isotope( "H3", Z=1, N=3 );
  G4Isotope *iso_3He = new G4Isotope( "He3", Z=2, N=3 );
  
  G4Element *elH = man->FindOrBuildElement("H"); //Define hydrogen with natural isotopic abundances
  G4Element *elD = new G4Element( "Deuterium", "2H", ncomponents=1 ); //Define pure deuterium:
  elD->AddIsotope( iso_2H, abundance=100.0*perCent);
  G4Element *el3He = new G4Element("Helium3", "3He", ncomponents=1 ); //Define isotopically pure Helium-3 
  el3He->AddIsotope( iso_3He, abundance=100.0*perCent );
  G4Element *el3H = new G4Element("Tritium", "3H", ncomponents=1 ); //Define isotopically pure Tritium
  el3H->AddIsotope( iso_3H, abundance=100.0*perCent );
  G4Element *el4He = man->FindOrBuildElement("He");
  G4Element *elC = man->FindOrBuildElement("C");
  G4Element *elF = man->FindOrBuildElement("F");
  G4Element *elNa = man->FindOrBuildElement("Na");
  G4Element *elAl = man->FindOrBuildElement("Al");
  G4Element *elSi = man->FindOrBuildElement("Si");
  G4Element *elCa = man->FindOrBuildElement("Ca");
  G4Element *elFe = man->FindOrBuildElement("Fe");
  G4Element *elSr = man->FindOrBuildElement("Sr");
  G4Element *elBa = man->FindOrBuildElement("Ba");

  G4Element *elZn = man->FindOrBuildElement("Zn");
  G4Element *elTi = man->FindOrBuildElement("Ti");
  G4Element *elCu = man->FindOrBuildElement("Cu");
  G4Element *elMg = man->FindOrBuildElement("Mg");
  
  G4Element* elCl  = man->FindOrBuildElement("Cl");
  G4Element* elAr  = man->FindOrBuildElement("Ar");

  G4Element* elCr  = man->FindOrBuildElement("Cr");
  G4Element* elMn   =  man->FindOrBuildElement("Mn");
  G4Element* elNi  = man->FindOrBuildElement("Ni");
  G4Element *elP  = man->FindOrBuildElement("P");
  G4Element *elS  = man->FindOrBuildElement("S");
  //G4Element *elW = man->FindOrBuildElement("W");

  G4Material *Vacuum =new G4Material(name="Vacuum", z=1., a=1.0*g/mole, density=1e-9*g/cm3);
  //Vacuum->SetMaterialPropertiesTable(Std_MPT);
  fMaterialsMap["Vacuum"] = Vacuum;

  if( fMaterialsMap.find("Lead") == fMaterialsMap.end() ){ 
    fMaterialsMap["Lead"] = man->FindOrBuildMaterial( "G4_Pb" );
  }
  if( fMaterialsMap.find("Aluminum") == fMaterialsMap.end() ){ 
    //fMaterialsMap["Aluminum"] = new G4Material(name="Aluminum", z=13., a=26.98*g/mole, density=2.7*g/cm3);
    fMaterialsMap["Aluminum"] = man->FindOrBuildMaterial( "G4_Al" );
    fMaterialsMap["Aluminum"]->SetName("Aluminum");
  }
  if( fMaterialsMap.find("Silicon") == fMaterialsMap.end() ){ 
    //fMaterialsMap["Silicon"] = new G4Material(name="Silicon", z=14., 28.086*g/mole, density=2.33*g/cm3);
    fMaterialsMap["Silicon"] = man->FindOrBuildMaterial( "G4_Si" );
  }
  density = 4.51*g/cm3;
  // G4Material* CsI = new G4Material(name="CsI", density, nel = 2);
  // CsI->AddElement(elI, .5);
  // CsI->AddElement(elCs,.5);
  //Use standardized versions wherever possible:
  G4Material *CsI = man->FindOrBuildMaterial("G4_CESIUM_IODIDE");

  fMaterialsMap["CsI"] = CsI;

  a = 55.85*g/mole;
  density = 7.87*g/cm3;
  //fMaterialsMap["Fer"] = new G4Material(name="Fer", z=26., a, density);

  fMaterialsMap["Iron"] = man->FindOrBuildMaterial("G4_Fe");
  fMaterialsMap["Fer"] = fMaterialsMap["Iron"];
  
  density = 1.29e-03*g/cm3; //This air is at a temp. of 0 deg. C 
  // G4Material* Air = new G4Material(name="Air", density, nel=2);
  // Air->AddElement(elN, .7);
  // Air->AddElement(elO, .3);
  //Use standardized G4_AIR, which corresponds to 20 degrees C and atmospheric pressure

  G4Material *Air = man->FindOrBuildMaterial("G4_AIR");

  //Air->SetMaterialPropertiesTable(Std_MPT);

  fMaterialsMap["Air"] = Air;

  // Montgomery July 2018 gold for mtpc
  G4Element* elAu  = man->FindOrBuildElement("Au");
  G4Material *Au = man->FindOrBuildMaterial("G4_Au");
  fMaterialsMap["Au"] = Au;
  // bonus electronics readout board
  double BonusPCBDen = 1.27*g/cm3;
  G4Material* BonusPCB = new G4Material("BonusPCB", BonusPCBDen, 4 );
  BonusPCB->AddElement(elH, 37);
  BonusPCB->AddElement(elC, 24);
  BonusPCB->AddElement(elO, 6);
  BonusPCB->AddElement(elN, 2);
  fMaterialsMap["BonusPCB"] = BonusPCB;


  G4Material *G4_polystyrene = man->FindOrBuildMaterial( "G4_POLYSTYRENE" );
  fMaterialsMap["POLYSTYRENE"] = G4_polystyrene;
    
  double bigden = 1e9*g/cm3;

  // Cell Glass - GE180 Aluminosilicate Glass
  //Changed names of materials in this composition since the G4Materials used here are only used to make up GE180:
  // SiO2 60.3%
  G4Material* SiO2 = new G4Material("GE180_SiO2", 2.2*g/cm3, 2 );
  SiO2->AddElement(elSi, 1);
  SiO2->AddElement(elO, 2);
  fMaterialsMap["GE180_SiO2"] = SiO2;

  // BaO  18.2%
  G4Material* BaO = new G4Material("GE180_BaO", bigden, 2 );
  BaO->AddElement(elBa, 1);
  BaO->AddElement(elO, 1);
  fMaterialsMap["GE180_BaO"] = BaO;
  // Al2O3 14.3%
  G4Material* Al2O3 = new G4Material("GE180_Al2O3", bigden, 2 );
  Al2O3->AddElement(elAl, 2);
  Al2O3->AddElement(elO, 3);
  fMaterialsMap["GE180_Al2O3"] = Al2O3;
  // CaO   6.5%
  G4Material* CaO = new G4Material("GE180_CaO", bigden, 2 );
  CaO->AddElement(elCa, 1);
  CaO->AddElement(elO, 1);
  fMaterialsMap["GE180_CaO"] = CaO;
  // SrO   0.25%
  G4Material* SrO = new G4Material("GE180_SrO", bigden, 2 );
  SrO->AddElement(elSr, 1);
  SrO->AddElement(elO, 1);
  fMaterialsMap["GE180_SrO"] = SrO;

  // Density 2.76 g/cm^3
  // Index of Refraction 1.536
  G4Material* GE180 = new G4Material("GE180", 2.76*g/cm3, 5);
  GE180->AddMaterial(SiO2, 0.6039);
  GE180->AddMaterial(BaO, 0.1829);
  GE180->AddMaterial(Al2O3, 0.1439);
  GE180->AddMaterial(CaO, 0.0659);
  GE180->AddMaterial(SrO, 0.0034);
  fMaterialsMap["GE180"] = GE180;
  //
  density = 1.19*g/cm3;
  G4Material* Acrylic = new G4Material(name="Acrylic", density, nel=3);
  Acrylic->AddElement(elC, 5);
  Acrylic->AddElement(elH, 8);
  Acrylic->AddElement(elO, 2);
  fMaterialsMap["Acrylic"] = Acrylic;


  //--------- GEM Materials  ---------
  // (stolen from GEMC)

  G4Material* NOMEX_pure = new G4Material("NOMEX_pure", density = 1.38*g/cm3, 5);
  NOMEX_pure -> AddElement(elH,0.04);
  NOMEX_pure -> AddElement(elC,0.54);
  NOMEX_pure -> AddElement(elN,0.09);
  NOMEX_pure -> AddElement(elO,0.10);
  NOMEX_pure -> AddElement(elCl,0.23);
  fMaterialsMap["NOMEX_pure"] = NOMEX_pure;

  G4Material* NOMEX = new G4Material("NOMEX",density = 0.04*g/cm3, 2);
  NOMEX -> AddMaterial(NOMEX_pure,0.45);
  NOMEX -> AddMaterial(Air,0.55);
  fMaterialsMap["NOMEX"] = NOMEX;

  G4Material* NEMAG10 = new G4Material("NEMAG10", 1.70*g/cm3, nel=4);
  NEMAG10 -> AddElement(elSi, 1);
  NEMAG10 -> AddElement(elO , 2);
  NEMAG10 -> AddElement(elC , 3);
  NEMAG10 -> AddElement(elH , 3);
  fMaterialsMap["NEMAG10"] = NEMAG10;

  G4double density_Ar = 1.7823*mg/cm3 ;
  G4Material* Argon = new G4Material("Argon"  , density_Ar, nel=1);
  Argon->AddElement(elAr, 1);
  fMaterialsMap["Argon"] = Argon;

  G4double density_CO2 = 1.977*mg/cm3;
  G4Material* CO2 = new G4Material("CO2", density_CO2, nel=2);
  CO2->AddElement(elC, 1);
  CO2->AddElement(elO, 2);


  G4double CO2_RefractiveIndex[nEntries];
  for ( int i = 0; i < nEntries; ++i ) {
    CO2_RefractiveIndex[i]=1.00045;
    CO2_RefractiveIndex[i]=1+(CO2_RefractiveIndex[i]-1)*GC_Gas_Pressure/atmosphere;
  }
  G4MaterialPropertiesTable* MPCO2 = new G4MaterialPropertiesTable();
  MPCO2->AddProperty("RINDEX",PhotonEnergy, CO2_RefractiveIndex, nEntries);

  if( fMaterialsListOpticalPhotonDisabled.find( "CO2" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    CO2->SetMaterialPropertiesTable(MPCO2);
  }

  fMaterialsMap["CO2"] = CO2;

  // 1.5 Atmosphere C4F8O for Cherenkov => changed to 1 atm : 2017/03/13
  G4double density_C4F8O = 9.64*mg/cm3; // density at 1ATM
  //G4Material* C4F8O = new G4Material("C4F8O", density_C4F8O*1.5, nel=3);
  G4Material* C4F8O = new G4Material("C4F8O", density_C4F8O, nel=3);// 
  C4F8O->AddElement(elC, 4);
  C4F8O->AddElement(elF, 8);
  C4F8O->AddElement(elO, 1);

  G4double C4F8O_RefractiveIndex[nEntries]; //at 1.atm
  G4double C4F8O_ABSLENGTH[nEntries]; //at 1.atm
  for ( int i = 0; i < nEntries; ++i ) {
    C4F8O_RefractiveIndex[i]=1.00135;
    C4F8O_ABSLENGTH[i]=7424.75*cm;
    C4F8O_RefractiveIndex[i]=1+(C4F8O_RefractiveIndex[i]-1)*GC_Gas_Pressure/atmosphere;
  }
  G4MaterialPropertiesTable* MPC4F8O = new G4MaterialPropertiesTable();
  MPC4F8O->AddProperty("RINDEX",PhotonEnergy, C4F8O_RefractiveIndex, nEntries);
  MPC4F8O->AddProperty("ABSLENGTH",PhotonEnergy, C4F8O_ABSLENGTH, nEntries);

  if( fMaterialsListOpticalPhotonDisabled.find( "C4F8O" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    C4F8O->SetMaterialPropertiesTable(MPC4F8O);
  }
  
  fMaterialsMap["C4F8O"] = C4F8O;

  G4double density_ArCO2 = .7*density_Ar + .3*density_CO2;
  // Use ArCO2
  G4Material *GEMgas= new G4Material("GEMgas", density_ArCO2, nel=2);
  GEMgas->AddMaterial(Argon, 0.7*density_Ar/density_ArCO2) ;
  GEMgas->AddMaterial(CO2, 0.3*density_CO2/density_ArCO2) ;

  fMaterialsMap["GEMgas"] = GEMgas;

  //G4Material *Copper= new G4Material("Copper", z=29, a=   63.55*g/mole, density = 8.96*g/cm3);
  G4Material *Copper = man->FindOrBuildMaterial( "G4_Cu" );
  fMaterialsMap["Copper"] = Copper;

  G4Material *Kapton = new G4Material("Kapton",   density = 1.42*g/cm3, nel=4);
  Kapton->AddElement(elH, 0.026362);
  Kapton->AddElement(elC, 0.691133);
  Kapton->AddElement(elN, 0.073270);
  Kapton->AddElement(elO, 0.209235);

  fMaterialsMap["Kapton"] = Kapton;

  //CH2 for FPP analyzers for GEP:
  G4double density_CH2 = 0.95*g/cm3;
  G4Material* CH2 = new G4Material("CH2", density_CH2, nel=2);
  CH2->AddElement(elC, 1);
  CH2->AddElement(elH, 2);

  fMaterialsMap["CH2"] = CH2;

  G4double density_CH = 0.95*g/cm3;
  G4Material* CH = new G4Material("CH", density_CH, nel=2);
  CH->AddElement(elC, 1);
  CH->AddElement(elH, 1);

  //adding a random comment here to test git push issues
  
  fMaterialsMap["CH"] = CH;

  //Target materials: these are hard-coded in two different places, leading to inconsistent calculations for luminosity, etc;
  //It would be better to figure out a sensible way for the G4SBSDetectorConstruction and the G4SBSEventGen classes to talk to each other
  //By the time /g4sbs/run is invoked, all the geometry should be frozen; therefore, the messenger should invoke the SetTargDen method of G4SBSEventGen at
  // that time, and grabbing the density of the relevant material from the fdetcon->fTargetBuilder based on the user's target choice
  //This would eliminate inconsistencies between the actual material density used in GEANT4 and the density assumed in calculating luminosity and
  //normalizing to a rate:
  double gasden = 10.5*atmosphere*(1.0079*2*g/Avogadro)/(300*kelvin*k_Boltzmann);
  // G4Material *refH2 = new G4Material("refH2", gasden, 1 );
  // double H2den = 1.2*atmosphere*(1.00794*2*g/Avogadro)/(90.0*kelvin*k_Boltzmann);
  gasden = 4.0*atmosphere*(1.00794*2*g/Avogadro)/(293.0*kelvin*k_Boltzmann);
  G4Material *refH2 = new G4Material("refH2", gasden, 1 );
  
  refH2->AddElement(elH, 1);

  fMaterialsMap["refH2"] = refH2;

  //gasden = 1.0*atmosphere*(2.0141*2*g/Avogadro)/(77*kelvin*k_Boltzmann);
  //G4Material *refD2 = new G4Material("refD2", gasden, 1 );
  //double D2den = 1.2*atmosphere*(2.0141*2*g/Avogadro)/(90.0*kelvin*k_Boltzmann);
  gasden = 4.0*atmosphere*(2.0141*2*g/Avogadro)/(293.0*kelvin*k_Boltzmann);
  G4Material *refD2 = new G4Material("refD2", gasden, 1 );
  refD2->AddElement(elD, 1);

  fMaterialsMap["refD2"] = refD2;

  /*
  // TDIS: ATM HARD CODE HIGHER PRESSURE BUT NEED TO MAKE THIS A MESSENGER OPTION
  // Redundant with above
  gasden = 4.0*atmosphere*(1.0079*2*g/Avogadro)/(300.0*kelvin*k_Boltzmann);
  G4Material *mTPCH2 = new G4Material("mTPCH2", gasden, 1 );
  mTPCH2->AddElement(elH, 1);
  fMaterialsMap["mTPCH2"] = mTPCH2;

  gasden = 4.0*atmosphere*(2.0141*2*g/Avogadro)/(300.0*kelvin*k_Boltzmann);
  G4Material *mTPCD2 = new G4Material("mTPCD2", gasden, 1 );
  mTPCD2->AddElement(elD, 1);
  fMaterialsMap["mTPCD2"] = mTPCD2;
  //
  */
  gasden = 10.5*atmosphere*(14.0067*2*g/Avogadro)/(300*kelvin*k_Boltzmann);
  G4Material *refN2 = new G4Material("refN2", gasden, 1 );
  refN2->AddElement(elN, 1);

  fMaterialsMap["refN2"] = refN2;

  gasden = 10.77*atmosphere*(3.016*g/Avogadro)/(300*kelvin*k_Boltzmann);
  G4Material *pol3He = new G4Material("pol3He", gasden, 1 );
  pol3He->AddElement(el3He, 1);

  fMaterialsMap["pol3He"] = pol3He;

  double LH2den = 0.071*g/cm3;
  G4Material *LH2mat = new G4Material("LH2", LH2den, 1 );
  LH2mat->AddElement(elH, 1);

  fMaterialsMap["LH2"] = LH2mat;

  double LD2den = 162.4*kg/m3;
  G4Material *LD2mat = new G4Material("LD2", LD2den, 1 );
  LD2mat->AddElement(elD, 1);
  
  fMaterialsMap["LD2"] = LD2mat;
  
  //TPC/future targets? material
  //G4double density_4He = 0.1*atmosphere*(4.0026*g/Avogadro)/(90*kelvin*k_Boltzmann);
  G4double density_4He = 1.0*atmosphere*(4.0026*g/Avogadro)/(293.0*kelvin*k_Boltzmann);
  G4cout << "Density of Helium (g/cm3): " << density_4He/(g/cm3) << G4endl;
  G4Material *ref4He = new G4Material("ref4He", density_4He, 1 );
  ref4He->AddElement(el4He, 1);

  fMaterialsMap["ref4He"] = ref4He;
  
  // GEn polarized 3He materials (D Flay, July 2020)  
  // - includes materials for the full geometry
  // - most materials already defined above; adding carbon steel, 
  //   which is used for the magnetic shield, and Ultem which is for the target ladder  

  // AISI 1008 carbon steel
  // - details from http://www.iron-foundry.com/AISI-1008-SAE-UNS-G10080-Carbon-Steel-Foundry.html
  // - NOTE: may throw a warning because this doesn't add to 100% (adds to 99.8%) 
  G4Material *Carbon_Steel_1008 = new G4Material("Carbon_Steel_1008",7.872*g/cm3,5);
  //Carbon_Steel_1008->AddElement(elFe,0.9931);
  //Carbon_Steel_1008->AddElement(elMn,0.0030);
  //shouldn't be a big difference, but it will stop throwing warnings...
  Carbon_Steel_1008->AddElement(elFe,0.9941);
  Carbon_Steel_1008->AddElement(elMn,0.0040);
  Carbon_Steel_1008->AddElement(elC ,0.0010);
  Carbon_Steel_1008->AddElement(elS ,0.0005);
  Carbon_Steel_1008->AddElement(elP ,0.0004);
  fMaterialsMap["Carbon_Steel_1008"] = Carbon_Steel_1008;

  // Beam line exit materials (D Flay, Sept 2020) 

  // Aluminum 5052 alloy 
  // - details from https://unitedaluminum.com/united-aluminum-alloy-5052/ 
  G4Material *Aluminum_5052 = new G4Material("Aluminum_5052",2.68*g/cm3,8);
  Aluminum_5052->AddElement(elAl,0.9635);
  Aluminum_5052->AddElement(elSi,0.0025);
  Aluminum_5052->AddElement(elFe,0.0040);
  Aluminum_5052->AddElement(elCu,0.0010);
  Aluminum_5052->AddElement(elMn,0.0010);
  Aluminum_5052->AddElement(elMg,0.0250);
  Aluminum_5052->AddElement(elCr,0.0020);
  Aluminum_5052->AddElement(elZn,0.0010);
  fMaterialsMap["Aluminum_5052"] = Aluminum_5052; 

  // Aluminum 6061 alloy 
  // - details from https://unitedaluminum.com/united-aluminum-alloy-6061/
  G4Material *Aluminum_6061 = new G4Material("Aluminum_6061",2.70*g/cm3,9);
  //Aluminum_6061->AddElement(elAl,0.9635);
  //shouldn't be a big difference, but it will stop throwing warnings...
  Aluminum_6061->AddElement(elAl,0.9668);
  Aluminum_6061->AddElement(elSi,0.0060);
  Aluminum_6061->AddElement(elFe,0.0070);
  Aluminum_6061->AddElement(elCu,0.0028);
  Aluminum_6061->AddElement(elMn,0.0015);
  Aluminum_6061->AddElement(elMg,0.0100);
  Aluminum_6061->AddElement(elCr,0.0019);
  Aluminum_6061->AddElement(elZn,0.0025);
  Aluminum_6061->AddElement(elTi,0.0015);
  fMaterialsMap["Aluminum_6061"] = Aluminum_6061; 

  // Molybdenum.  Possibly use for GEn 3He target collimators?
  // density = 10.22 g/cm^3  
  fMaterialsMap["Molybdenum"] = man->FindOrBuildMaterial("G4_Mo"); 

  // D Flay's mock ion chamber material 
  // Nitrogen.  Take this from DF notes on ion chambers. 
  // This density is from a typical LHC device 
  G4double gasden_icN2 = 1.08*atmosphere*(14.0067*2*g/Avogadro)/(300*kelvin*k_Boltzmann);
  G4Material *icN2 = new G4Material("icN2",gasden_icN2,1);
  icN2->AddElement(elN,1);
  fMaterialsMap["GEnTarget_ionChamber_N2"] = icN2; 

  // Ultem (polyetherimide plastic, similar to PEEK)
  // - details from http://www.polymerprocessing.com/polymers/PEI.html
  G4Material *Ultem = new G4Material("Ultem",1.27*g/cm3,nel=4);
  Ultem->AddElement(elC,37); 
  Ultem->AddElement(elH,24);
  Ultem->AddElement(elO,6 );
  Ultem->AddElement(elN,2 ); 
  fMaterialsMap["Ultem"] = Ultem; 
  
  //G4double density_CH4 = 0.1*atmosphere*((12.0107+4*1.0079)*g/Avogadro)/(90*kelvin*k_Boltzmann);
  G4double density_CH4 = 1.0*atmosphere*((12.0107+4*1.0079)*g/Avogadro)/(293.0*kelvin*k_Boltzmann);
  G4Material *CH4 = new G4Material("CH4", density_CH4, nel=2 );
  CH4->AddElement(elC, 1);
  CH4->AddElement(elH, 4);

  G4cout << "H2 density (g/cm3) = " << refH2->GetDensity()/(g/cm3) 
	 << ", D2 density (g/cm3) = " << refD2->GetDensity()/(g/cm3) 
	 << ", 4He density (g/cm3) = " << ref4He->GetDensity()/(g/cm3) 
	 << ", CH4 density (g/cm3) = " << CH4->GetDensity()/(g/cm3) 
	 << G4endl;
  
  fMaterialsMap["CH4"] = CH4;

  G4double density_TPCgas = 0.7*density_4He+0.3*density_CH4;
  G4Material *TPCgas= new G4Material("TPCgas", density_TPCgas, nel=2);
  TPCgas->AddMaterial(ref4He, 0.7*density_4He/density_TPCgas) ;
  TPCgas->AddMaterial(CH4, 0.3*density_CH4/density_TPCgas) ;

  fMaterialsMap["TPCgas"] = TPCgas;


  // Adjustable mTPC gas
  //mtpc glabal variables default values
  fmTPCHeGasFraction = 0.7;//fraction of 1
  fmTPCCH4GasFraction = 0.3;//fraction of 1
  fmTPCGasTemp = 296.15;//K
  fmTPCGasPressure = 1.0;//atm
  G4double density_mTPCgas = fmTPCHeGasFraction*density_4He+fmTPCCH4GasFraction*density_CH4;
  G4Material *mTPCgas= new G4Material("mTPCgas", density_mTPCgas, nel=2);
  mTPCgas->AddMaterial(ref4He, fmTPCHeGasFraction*density_4He/density_mTPCgas) ;
  mTPCgas->AddMaterial(CH4, fmTPCCH4GasFraction*density_CH4/density_mTPCgas) ;
  fMaterialsMap["mTPCgas"] = mTPCgas;
  // G4cout <<  "density_4He " << ref4He->GetDensity()/(g/cm3) << G4endl;
  // G4cout <<  "density_mTPCgas " << mTPCgas->GetDensity()/(g/cm3) << G4endl;

  // // Adjustable mTPC gas
  // //mtpc glabal variables default values
  // fmTPCHeGasFraction = 0.7;//fraction of 1
  // fmTPCCH4GasFraction = 0.3;//fraction of 1
  // fmTPCGasTemp = 296.15;//K
  // fmTPCGasPressure = 1.0;//atm
  // //He4
  // G4double density_mTPC_4He = fmTPCGasPressure*atmosphere*(4.0026*g/Avogadro)/(fmTPCGasTemp*kelvin*k_Boltzmann);
  // G4Material *refmTPC4He = new G4Material("refmTPC4He", density_mTPC_4He, 1 );
  // refmTPC4He->AddElement(el4He, 1);
  // G4cout << "density of mtpc he gas " << density_mTPC_4He << G4endl;
  // fMaterialsMap["refmTPC4He"] = refmTPC4He;
  // //CH4
  // G4double density_mTPC_CH4 = fmTPCGasPressure*atmosphere*((12.0107+4*1.0079)*g/Avogadro)/(fmTPCGasTemp*kelvin*k_Boltzmann);
  // G4Material *mTPCCH4 = new G4Material("mTPCCH4", density_mTPC_CH4, nel=2 );
  // mTPCCH4->AddElement(elC, 1);
  // mTPCCH4->AddElement(elH, 4);
  // G4cout << "density of mtpc ch4 gas " << density_mTPC_CH4 << G4endl;
  // fMaterialsMap["mTPCCH4"] = mTPCCH4;
  // // He4 and CH4 mix
  // G4double density_mTPCgas = fmTPCHeGasFraction*density_mTPC_4He + fmTPCCH4GasFraction*density_mTPC_CH4;
  // G4Material *mTPCgas= new G4Material("mTPCgas", density_mTPCgas, nel=2);
  // mTPCgas->AddMaterial(refmTPC4He, fmTPCHeGasFraction*density_mTPC_4He/density_mTPCgas) ;
  // mTPCgas->AddMaterial(mTPCCH4, fmTPCCH4GasFraction*density_mTPC_CH4/density_mTPCgas) ;
  // fMaterialsMap["mTPCgas"] = mTPCgas;
  // G4cout << "density of mtpc gas " << density_mTPCgas << G4endl;


  //Beamline materials:
  density = 2.5*g/cm3;
  G4Material *Concrete = new G4Material("Concrete",density,6);
  Concrete->AddElement(elO, 0.52);
  Concrete->AddElement(elSi, 0.325);
  Concrete->AddElement(elCa, 0.06);
  Concrete->AddElement(elNa, 0.015);
  Concrete->AddElement(elFe, 0.04);
  Concrete->AddElement(elAl, 0.04);

  fMaterialsMap["Concrete"] = Concrete;

  density = 8.02*g/cm3 ;
  G4Material *stainless = new G4Material("Stainless steel",density,5);
  stainless->AddElement(elMn, 0.02);
  stainless->AddElement(elSi, 0.01);
  stainless->AddElement(elCr, 0.19);
  stainless->AddElement(elNi, 0.10);
  stainless->AddElement(elFe, 0.68);

  fMaterialsMap["Stainless"] = stainless;

  density = 8.02*g/cm3 ;
  //G4Material* matBe = new G4Material("Berylium", 4., 9.012*g/mole, density=1.85*g/cm3);
  G4Material *matBe = man->FindOrBuildMaterial( "G4_Be" );
  fMaterialsMap["Beryllium"] = matBe;

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

  G4Element *S  = man->FindOrBuildElement("S"); //sulfur
  
  //Define tedlar for gaps between aerogel tiles. These will have no optical properties and therefore optical photons 
  //crossing tile boundaries will be killed.
  G4Material *Tedlar = man->FindOrBuildMaterial("G4_POLYVINYLIDENE_FLUORIDE");
  fMaterialsMap["Tedlar"] = Tedlar;

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
  G4double Ephoton_lucite[nentries_lucite] = 
    {   1.719888896*eV, 1.81647209*eV,  1.924545279*eV, 2.054623669*eV, 2.293705949*eV, 
	2.598432583*eV, 2.900093346*eV, 3.145984951*eV, 3.465763144*eV, 3.704958425*eV, 
	3.973356238*eV, 4.049879001*eV, 4.122678976*eV, 4.184217105*eV, 4.219202221*eV, 
	4.276408832*eV, 4.320338518*eV, 4.36519548*eV,  4.410978*eV,    4.521654214*eV, 
	4.521654214*eV
    };

  G4double Abslength_lucite[nentries_lucite] = 
    {   4.359093727*cm, 4.478062096*cm, 4.603526015*cm, 4.603526015*cm, 4.478062096*cm, 
	4.359093727*cm, 4.035983645*cm, 3.938294782*cm, 3.51091739*cm,  3.164087254*cm, 
	2.727358613*cm, 2.319695564*cm, 1.529808225*cm, 0.873567226*cm, 0.532172776*cm, 
	0.28887992*cm,  0.148639788*cm, 0.093529994*cm, 0.07329661*cm,  0.046324745*cm, 
	0.046324745*cm };

  G4double Ephoton_rindex_lucite[2] = {1.77*eV, 6.2*eV}; 
  G4double Rindex_lucite[2] = {1.5, 1.5};

  MPT_temp->AddProperty( "ABSLENGTH", Ephoton_lucite, Abslength_lucite, nentries_lucite );
  MPT_temp->AddProperty( "RINDEX", Ephoton_rindex_lucite, Rindex_lucite, 2 );

  if( fMaterialsListOpticalPhotonDisabled.find( "UVT_Lucite" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    UVT_Lucite->SetMaterialPropertiesTable( MPT_temp );
  }
  // G4cout << "Material optical properties for material UVT_Lucite" << G4endl;
  // MPT_temp->DumpTable();

  fMaterialsMap["UVT_Lucite"] = UVT_Lucite;


  // Grinch Lucite

  G4Material* Lucite=man->FindOrBuildMaterial("G4_LUCITE");
  G4double Lucite_RefractiveIndex[nEntries];
  for ( int i = 0; i < nEntries; ++i ) {
    Lucite_RefractiveIndex[i]=1.5;
  }
  G4MaterialPropertiesTable* MPLucite = new G4MaterialPropertiesTable();
  MPLucite->AddProperty("RINDEX",PhotonEnergy, Lucite_RefractiveIndex, nEntries);
  if( fMaterialsListOpticalPhotonDisabled.find( "Lucite" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    Lucite->SetMaterialPropertiesTable(MPLucite);
  }
  Lucite->SetName("Lucite");
  fMaterialsMap["Lucite"] = Lucite;


  G4Material* Glass=new G4Material("Glass",1.397*g/cm3,3);
  Glass->AddElement(H,8);
  Glass->AddElement(C,10);
  Glass->AddElement(O,4);
  fMaterialsMap["Glass"] = Glass;

  G4Material* Polycarbonate=man->FindOrBuildMaterial("G4_POLYCARBONATE");
  fMaterialsMap["PolyCarbonate"] = Polycarbonate;
 
  
  G4Material* matAl=man->FindOrBuildMaterial("G4_Al");
  matAl->SetName("Al");
  fMaterialsMap["Al"] = matAl;

  G4Element *Fe = man->FindOrBuildElement("Fe");
  G4Element *Cr = man->FindOrBuildElement("Cr");
  G4Element *Ni = man->FindOrBuildElement("Ni");
  G4Element *Mn = man->FindOrBuildElement("Mn");
  G4Element *P = man->FindOrBuildElement("P");
  //G4Element *S = man->FindOrBuildElement("S");
  G4Element *Mo = man->FindOrBuildElement("Mo");

  G4Material* Carbon_Steel = new G4Material("Carbon_Steel", density=7.85*g/cm3, 7);
  Carbon_Steel->AddElement(Fe, 73.4*perCent);
  Carbon_Steel->AddElement(Cr, 8.*perCent);
  Carbon_Steel->AddElement(Ni, 18.*perCent);
  Carbon_Steel->AddElement(C, 0.1*perCent);
  Carbon_Steel->AddElement(Mn, 0.41*perCent);
  Carbon_Steel->AddElement(P, 0.04*perCent);
  Carbon_Steel->AddElement(S, 0.05*perCent);
  fMaterialsMap["Carbon_Steel"] = Carbon_Steel;

  // G4Material* Mylar=new G4Material("Mylar", density= 1.397*g/cm3, 3);
  // Mylar->AddElement(C, 10);
  // Mylar->AddElement(H, 8);
  // Mylar->AddElement(O, 4);
  // fMaterialsMap["Mylar"] = Mylar;

  G4Material* mu_metal=new G4Material("mu-metal", density= 8.250*g/cm3, 6);
  mu_metal->AddElement(C, 0.02*perCent);
  mu_metal->AddElement(Ni, 80.00*perCent);
  mu_metal->AddElement(Mn, 0.50*perCent);
  mu_metal->AddElement(Mo, 4.20*perCent);
  mu_metal->AddElement(Si, 0.35*perCent);
  mu_metal->AddElement(Fe, 14.93*perCent);
  fMaterialsMap["mu-metal"] = mu_metal;



  G4Material *Steel = man->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  fMaterialsMap["Steel"] = Steel;
  fMaterialsMap["Stainless_Steel"] = Steel;

  G4Material *Carbon = man->FindOrBuildMaterial("G4_GRAPHITE");
  fMaterialsMap["Carbon"] = Carbon;
  
  
  //Need to define "QuartzWindow" *before* Aerogel because 
  //the Aerogel Material properties definition refers to Quartz Window:

  // For the time being, we treat the mirror as being a single shell 
  // of Carbon, with a uniform thickness corresponding to 0.01 X0:

  G4double den_composite = 2.21*g/cm3;

  G4Material *MirrorComposite = new G4Material( "MirrorComposite", den_composite, nel=1 );
  MirrorComposite->AddElement( C, fractionmass = 1.0 );
  fMaterialsMap["MirrorComposite"] = MirrorComposite;

  //Next we will need to define Aluminum for the RICH entry and exit windows and other parts of the containment volume.
  //Let's use the standard NIST database parameters for Al "G4_Al";
  G4Material *RICHAluminum = man->FindOrBuildMaterial("G4_Al");
  fMaterialsMap["RICHAluminum"] = RICHAluminum;

  G4double den_pmtwindow = 2.20*g/cm3;

  G4Material *QuartzWindow = new G4Material( "QuartzWindow", den_pmtwindow, nel=2 );
  QuartzWindow->AddElement( O, natoms=2 );
  QuartzWindow->AddElement( Si, natoms=1 );

  //Define refractive index and absorption length for quartz PMT windows:
  const G4int nentries_quartz = 51;

  G4double Ephoton_quartz[nentries_quartz] = { 
    1.77120301*eV, 1.796872619*eV, 1.823297216*eV, 1.850510608*eV, 1.878548647*eV, 
    1.907449396*eV, 1.937253292*eV, 1.968003345*eV, 1.999745334*eV, 2.032528044*eV, 
    2.066403512*eV, 2.1014273*eV, 2.137658805*eV, 2.175161591*eV, 2.214003763*eV, 
    2.254258377*eV, 2.296003902*eV, 2.33932473*eV, 2.384311744*eV, 2.431062955*eV, 
    2.479684214*eV, 2.530290015*eV, 2.58300439*eV, 2.63796193*eV, 2.695308928*eV, 
    2.755204682*eV, 2.817822971*eV, 2.883353737*eV, 2.952005017*eV, 3.024005139*eV, 
    3.099605268*eV, 3.179082326*eV, 3.262742387*eV, 3.350924614*eV, 3.444005853*eV, 
    3.54240602*eV, 3.646594433*eV, 3.757097294*eV, 3.874506585*eV, 3.999490668*eV, 
    4.132807024*eV, 4.275317611*eV, 4.428007525*eV, 4.592007804*eV, 4.768623489*eV, 
    4.959368428*eV, 5.16600878*eV, 5.390617857*eV, 5.635645941*eV, 5.904010034*eV, 
    6.199210536*eV 
  };

  //Refractive index data from "Refractiveindex.info"
  G4double Rindex_quartz[nentries_quartz] = { 
    1.455292466, 1.455524071, 1.455763571, 1.456011496, 1.456268423, 1.456534974, 1.456811819, 1.457099689, 
    1.457399374, 1.457711733, 1.458037702, 1.4583783, 1.458734641, 1.459107942, 1.459499536, 1.459910886, 
    1.460343603, 1.460799458, 1.461280408, 1.461788618, 1.462326487, 1.462896682, 1.463502175, 1.464146283, 
    1.464832722, 1.465565665, 1.466349815, 1.467190482, 1.46809369, 1.469066293, 1.470116119, 1.471252144, 
    1.472484709, 1.473825777, 1.475289258, 1.476891413, 1.478651361, 1.48059172, 1.482739429, 1.485126813, 
    1.487792976, 1.490785646, 1.494163661, 1.498000361, 1.502388312, 1.507446007, 1.513327606, 1.520237459, 
    1.528452449, 1.53835762, 1.550505538 
  };

  //Typical absorption length for fused silica from Thorlabs.com 
  //(uncoated UV fused silica, includes fresnel reflections)
  // G4double Abslength_quartz[nentries_quartz] = {
  //   15.65792444*cm, 15.78550788*cm, 15.7794917*cm, 15.60910249*cm, 15.72664954*cm, 15.72488912*cm, 15.57290011*cm, 15.68021339*cm, 
  //   15.73546266*cm, 15.55685196*cm, 15.55490625*cm, 15.63907251*cm, 15.48113765*cm, 15.54074565*cm, 15.39638598*cm, 15.50169846*cm, 
  //   15.4950396*cm, 15.36125979*cm, 15.41113687*cm, 15.33874196*cm, 15.24165927*cm, 15.25602267*cm, 15.23330157*cm, 15.14071666*cm, 
  //   15.13642486*cm, 15.06590584*cm, 15.05023293*cm, 14.99006002*cm, 14.91826095*cm, 14.79500397*cm, 14.80590431*cm, 14.66396966*cm, 
  //   14.57959363*cm, 14.47194788*cm, 14.40952367*cm, 14.20967861*cm, 14.11981056*cm, 13.98888512*cm, 13.79714319*cm, 13.85187177*cm, 
  //   13.61931079*cm, 13.27721911*cm, 12.75893298*cm, 12.51276543*cm, 12.13753406*cm, 11.84114847*cm, 10.96932156*cm, 10.13778162*cm, 
  //   9.92333989*cm, 9.282597031*cm, 7.443241349*cm 
  // };
  //Internal transmittance data for UV grade fused silica from Melles-Griot (bulk absorption is actually much smaller than implied by external transmittance data; including this leads to "double-counting" absorption that is actually Fresnel reflections):
  const G4int Nabs_quartz = 10;
  G4double Ephoton_abs_quartz[Nabs_quartz] = {1.7491*eV,1.8793*eV,2.0012*eV,2.1734*eV,2.2993*eV,2.4195*eV,4.0372*eV,4.8901*eV,6.0646*eV,6.4897*eV};
  G4double Abslength_quartz[Nabs_quartz] = {522.51*cm,2666.17*cm,1633.49*cm,756.50*cm,418.96*cm,233.58*cm,23.46*cm,22.58*cm,15.76*cm,14.74*cm};
  
  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("RINDEX", Ephoton_quartz, Rindex_quartz, nentries_quartz);
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_abs_quartz, Abslength_quartz, Nabs_quartz);

  if( fMaterialsListOpticalPhotonDisabled.find( "QuartzWindow" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    QuartzWindow->SetMaterialPropertiesTable( MPT_temp );
  }
  fMaterialsMap["QuartzWindow"] = QuartzWindow;

  //Next define borosilicate "UV glass" for the PMT window. For chemical 
  //composition, use data for Pyrex Corning 7740 borosilicate glass from 
  //PDG, but assign same optical properties as quartz above.

  G4double den_UVglass = 2.23*g/cm3;

  G4Material *UVglass = new G4Material( "UVglass", den_UVglass, nel=6 );
  UVglass->AddElement( B, fractionmass=0.040061 );
  UVglass->AddElement( O, fractionmass=0.539564 );
  UVglass->AddElement( Na, fractionmass=0.028191 );
  UVglass->AddElement( Al, fractionmass=0.011644 );
  UVglass->AddElement( Si, fractionmass=0.377220 );
  UVglass->AddElement( K,  fractionmass=0.003321 );

  if( fMaterialsListOpticalPhotonDisabled.find( "UVglass" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    UVglass->SetMaterialPropertiesTable( MPT_temp );
  }

  fMaterialsMap["UVglass"] = UVglass;

  G4double den_C4F10 = 10.62*mg/cm3; //Wow, that really is a heavy gas,
  //We also need to define the radiator gas C4F10.
  G4Material *C4F10_gas = new G4Material( "C4F10_gas", den_C4F10, nel=2 );

  C4F10_gas->AddElement( C, natoms=4 );
  C4F10_gas->AddElement( F, natoms=10 );

  MPT_temp = new G4MaterialPropertiesTable();
  //MPT_temp->AddConstProperty( "RINDEX", 1.00137 ); //assume a constant refractive index for the gas radiator for now. 
  //MPT_temp->AddConstProperty( "ABSLENGTH", 100.0*m ); //For some reason, defining a constant refractive index does not seem to work.
  const G4int nentries_C4F10 = 2;
  G4double Ephoton_C4F10[nentries_C4F10] = { 1.77*eV, 6.20*eV };
  G4double Rindex_C4F10[nentries_C4F10] = {1.00137, 1.00137};
  G4double Abslength_C4F10[nentries_C4F10] = {100.0*m, 100.0*m};

  MPT_temp->AddProperty("RINDEX", Ephoton_C4F10, Rindex_C4F10, nentries_C4F10 );
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_C4F10, Abslength_C4F10, nentries_C4F10 );

  if( fMaterialsListOpticalPhotonDisabled.find( "C4F10_gas" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    C4F10_gas->SetMaterialPropertiesTable( MPT_temp );
  }
  
  fMaterialsMap["C4F10_gas"] = C4F10_gas;
  
  //EFuchey: 2019/08/15: Add C4F8
  G4double den_C4F8 = 8.82*mg/cm3; 
  G4Material *C4F8_gas = new G4Material( "C4F8_gas", den_C4F8, nel=2 );
  
  C4F8_gas->AddElement( C, natoms=4 );
  C4F8_gas->AddElement( F, natoms=8 );

  MPT_temp = new G4MaterialPropertiesTable();
  const G4int nentries_C4F8 = 2;
  G4double Ephoton_C4F8[nentries_C4F8] = { 1.95*eV, 5.00*eV };
  G4double Rindex_C4F8[nentries_C4F8] = {1.00132, 1.00132};
  G4double Abslength_C4F8[nentries_C4F8] = {100.0*m, 100.0*m};

  MPT_temp->AddProperty("RINDEX", Ephoton_C4F8, Rindex_C4F8, nentries_C4F8 );
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_C4F8, Abslength_C4F8, nentries_C4F8 );

  if( fMaterialsListOpticalPhotonDisabled.find( "C4F8_gas" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    C4F8_gas->SetMaterialPropertiesTable( MPT_temp );
  }
  fMaterialsMap["C4F8_gas"] = C4F8_gas;
   
  
  //EFuchey: 2017/07/24: Add CF4, just in case we want to study it as an option for GRINCH
  G4double den_CF4 = 3.72*mg/cm3; 
  G4Material *CF4_gas = new G4Material( "CF4_gas", den_CF4, nel=2 );
  
  CF4_gas->AddElement( C, natoms=1 );
  CF4_gas->AddElement( F, natoms=4 );

  MPT_temp = new G4MaterialPropertiesTable();
  const G4int nentries_CF4 = 2;
  G4double Ephoton_CF4[nentries_CF4] = { 1.95*eV, 5.00*eV };
  G4double Rindex_CF4[nentries_CF4] = {1.00048, 1.00048};
  G4double Abslength_CF4[nentries_CF4] = {100.0*m, 100.0*m};

  MPT_temp->AddProperty("RINDEX", Ephoton_CF4, Rindex_CF4, nentries_CF4 );
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_CF4, Abslength_CF4, nentries_CF4 );

  if( fMaterialsListOpticalPhotonDisabled.find( "CF4_gas" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    CF4_gas->SetMaterialPropertiesTable( MPT_temp );
  }
  fMaterialsMap["CF4_gas"] = CF4_gas;
  //SF6 = 146.05 g/mol, P = rho RT / Mmol, rho = Mmol*P/RT;
  G4double Mmol_SF6 = 146.05*gram;
  G4double density_SF6 = Mmol_SF6 * atmosphere / (Avogadro * k_Boltzmann * 293.0*kelvin );

  G4cout << "Defining Sulfur Hexafluoride, density = " << density_SF6/g*cm3 << " g/cm^3" << G4endl;

  G4Material *mat_SF6 = new G4Material( "SF6_gas", density_SF6, nel=2 );

  mat_SF6->AddElement( S, natoms=1 );
  mat_SF6->AddElement( F, natoms=6 );
  
  const G4int nentries_SF6 = 2;
  G4double Ephoton_SF6[nentries_SF6] = { 1240.0/700.0*eV, 1240.0/200.0*eV };
  G4double Rindex_SF6[nentries_SF6] = { 1.000783, 1.000783 }; //neglecting dispersion for now
  G4double Abslength_SF6[nentries_SF6] = {10.0*m, 10.0*m };
  
  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("RINDEX", Ephoton_SF6, Rindex_SF6, nentries_SF6 );
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_SF6, Abslength_SF6, nentries_SF6 );

  if( fMaterialsListOpticalPhotonDisabled.find( "SF6_gas" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    mat_SF6->SetMaterialPropertiesTable( MPT_temp );
  }
  fMaterialsMap["SF6_gas"] = mat_SF6;
 
  
  G4double Ephoton_RICH_air[nentries_C4F10] = { 1.77*eV, 6.20*eV };
  G4double Rindex_RICH_air[nentries_C4F10] = { 1.000277, 1.000277 };
  G4double Abslength_RICH_air[nentries_C4F10] = {100.0*m, 100.0*m };

  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("RINDEX", Ephoton_RICH_air, Rindex_RICH_air, nentries_C4F10 );
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_RICH_air, Abslength_RICH_air, nentries_C4F10 );

  // G4Material *RICH_air = man->FindOrBuildMaterial("G4_AIR");
  G4Material *RICH_air = new G4Material( "RICH_air", 1.20479*kg/m3, 4 );
  RICH_air->AddElement( elC, fractionmass = 0.000124 );
  RICH_air->AddElement( elN, fractionmass = 0.755268 );
  RICH_air->AddElement( elO, fractionmass = 0.231781 );
  RICH_air->AddElement( elAr, fractionmass = 0.012827 );
  
  if( fMaterialsListOpticalPhotonDisabled.find( "RICH_air" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    RICH_air->SetMaterialPropertiesTable( MPT_temp );
  }

  fMaterialsMap["RICH_air"] = RICH_air;
  
  //Quantum efficiency for PMT photocathode ("typical", from XP1911/UV data sheet):
  const G4int nentries_QE_RICH = 40;

  G4double Ephoton_QE_RICH[nentries_QE_RICH] = { 
    1.939182049*eV, 1.963579994*eV, 1.999929682*eV, 2.0346944*eV, 2.067643295*eV, 2.098539267*eV, 2.127141474*eV, 2.153230717*eV, 
    2.179967872*eV, 2.224858112*eV, 2.264313372*eV, 2.293899315*eV, 2.339752453*eV, 2.359406493*eV, 2.403825804*eV, 2.432973578*eV, 
    2.484622785*eV, 2.538512459*eV, 2.599590307*eV, 2.673824751*eV, 2.768699841*eV, 2.858865976*eV, 2.986530569*eV, 3.105387967*eV, 
    3.168432345*eV, 3.279402652*eV, 3.349799677*eV, 3.526419346*eV, 3.674006182*eV, 3.953218421*eV, 4.139635948*eV, 4.385211047*eV, 
    4.570741196*eV, 4.871970386*eV, 5.177284871*eV, 5.437857883*eV, 5.611675974*eV, 5.845211394*eV, 5.969437475*eV, 6.234402273*eV }; 

  G4double PMT_QuantumEfficiency_RICH[nentries_QE_RICH] = { 
    0.001154539, 0.001948441, 0.003714689, 0.006298763, 0.009646797, 0.013657231, 0.018010571, 0.022819293, 
    0.028693173, 0.038099138, 0.048532161, 0.057843653, 0.075581257, 0.084283662, 0.105012092, 0.116629941, 
    0.132736748, 0.149971254, 0.16338945, 0.176043552, 0.193934134, 0.213040691, 0.227782388, 0.229627292, 
    0.223657738, 0.215914559, 0.213826088, 0.213229177, 0.200888412, 0.17403965, 0.161019419, 0.142756599, 
    0.126474694, 0.104423377, 0.080171292, 0.057628786, 0.041016469, 0.024094955, 0.012166848, 0.004610016 };

  const G4int nentries_QE_RICH_NIM = 16;

  G4double Ephoton_QE_RICH_NIM[nentries_QE_RICH_NIM] = {1.8321513*eV, 1.95954488*eV, 2.066666667*eV, 2.227011494*eV, 2.464228935*eV, 2.672413793*eV,
							2.980769231*eV, 3.304904051*eV, 3.655660377*eV, 4.057591623*eV, 4.454022989*eV, 4.98392283*eV,
							5.516014235*eV, 6.031128405*eV, 6.623931624*eV, 7.276995305*eV };
  
  G4double PMT_QuantumEfficiency_RICH_NIM[nentries_QE_RICH_NIM] = {0.0057971, 0.0217391, 0.0434783, 0.0913043,
								   0.162319, 0.214493, 0.26087, 0.288406,
								   0.294203, 0.282609, 0.256522, 0.215942,
								   0.169565, 0.12029, 0.0710145, 0.015942 };

  //Define another material for photocathode: 
  G4double den_photocathode = 2.57*g/cm3;
  //EFuchey 2017-07-24: declare a different photocathod material for RICH and GRINCH 
  //(so we may set them a different quantum efficiency )
  G4Material *Photocathode_material_RICH = new G4Material( "Photocathode_material_RICH", den_photocathode, nel=3 );

  Photocathode_material_RICH->AddElement( Sb, natoms=1 );
  Photocathode_material_RICH->AddElement( K, natoms=2 );
  Photocathode_material_RICH->AddElement( Cs, natoms=1 );

  //G4double Ephot_Rcathode[2] = {1.77*eV, 6.20*eV};
  //G4double Rcathode[2] = {0.0, 0.0};

  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("EFFICIENCY", Ephoton_QE_RICH_NIM, PMT_QuantumEfficiency_RICH_NIM, nentries_QE_RICH_NIM );
  MPT_temp->AddProperty("RINDEX", Ephoton_quartz, Rindex_quartz, nentries_quartz );
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_abs_quartz, Abslength_quartz, nentries_quartz );
  //MPT_temp->AddProperty("REFLECTIVITY", Ephot_Rcathode, Rcathode, 2 );

  if( fMaterialsListOpticalPhotonDisabled.find( "Photocathode_material_RICH" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    Photocathode_material_RICH->SetMaterialPropertiesTable( MPT_temp );
  }

  fMaterialsMap["Photocathode_material_RICH"] = Photocathode_material_RICH;
  
  const G4int nentries_QE_GRINCH = 50;
  
  G4double Ephoton_QE_GRINCH[nentries_QE_GRINCH] = {
    2.00*eV, 2.10*eV, 2.15*eV, 2.20*eV, 2.25*eV, 2.30*eV, 2.35*eV, 2.40*eV, 2.45*eV, 2.50*eV, 
    2.55*eV, 2.60*eV, 2.65*eV, 2.70*eV, 2.75*eV, 2.80*eV, 2.85*eV, 2.90*eV, 2.95*eV, 3.00*eV, 
    3.05*eV, 3.10*eV, 3.15*eV, 3.20*eV, 3.25*eV, 3.30*eV, 3.35*eV, 3.40*eV, 3.45*eV, 3.50*eV, 
    3.55*eV, 3.60*eV, 3.65*eV, 3.70*eV, 3.75*eV, 3.80*eV, 3.85*eV, 3.90*eV, 3.95*eV, 4.00*eV, 
    4.05*eV, 4.10*eV, 4.15*eV, 4.20*eV, 4.25*eV, 4.30*eV, 4.35*eV, 4.40*eV, 4.45*eV, 4.50*eV};

  G4double PMT_QuantumEfficiency_GRINCH[nentries_QE_GRINCH] = {
    0.004004980, 0.010049800, 0.021369829, 0.037179902, 0.054958921, 
    0.073774945, 0.092877242, 0.111676268, 0.129724628, 0.146699023, 
    0.162383195, 0.176651856, 0.189455605, 0.200806842, 0.210766665, 
    0.219432753, 0.226928250, 0.233391630, 0.238967550, 0.243798700, 
    0.248018636, 0.251745604, 0.255077357, 0.258086955, 0.260819559, 
    0.263290215, 0.265482625, 0.267348907, 0.268810346, 0.269759136, 
    0.270061107, 0.269559448, 0.268079408, 0.265434003, 0.261430698, 
    0.255879083, 0.248599542, 0.239432909, 0.228251110, 0.214968800, 
    0.199555986, 0.182051641, 0.162578307, 0.141357686, 0.118727222, 
    0.095157675, 0.071271677, 0.047863288, 0.025918527, 0.006636911};
  
  G4Material *Photocathode_material_GRINCH = new G4Material( "Photocathode_material_GRINCH", den_photocathode, nel=3 );
  
  Photocathode_material_GRINCH->AddElement( Sb, natoms=1 );
  Photocathode_material_GRINCH->AddElement( K, natoms=2 );
  Photocathode_material_GRINCH->AddElement( Cs, natoms=1 );

  //G4double Ephot_Rcathode[2] = {1.77*eV, 6.20*eV};
  //G4double Rcathode[2] = {0.0, 0.0};

  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("EFFICIENCY", Ephoton_QE_GRINCH, PMT_QuantumEfficiency_GRINCH, nentries_QE_GRINCH );
  MPT_temp->AddProperty("RINDEX", Ephoton_quartz, Rindex_quartz, nentries_quartz );
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_abs_quartz, Abslength_quartz, nentries_quartz );
  //MPT_temp->AddProperty("REFLECTIVITY", Ephot_Rcathode, Rcathode, 2 );

  if( fMaterialsListOpticalPhotonDisabled.find( "Photocathode_material_GRINCH" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    Photocathode_material_GRINCH->SetMaterialPropertiesTable( MPT_temp );
  }

  fMaterialsMap["Photocathode_material_GRINCH"] = Photocathode_material_GRINCH;

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

  //G4double A_aerogel = 0.964; //for a thickness of 1 cm. This implies that A = e^{-1 cm/Labs}; ln A = -1 cm / Labs --> Labs = -1 cm / ln(A)
  G4double Ct_aerogel = 0.0094; //microns^4

  G4double t_aerogel = 0.5*(1.125+1.0); //nominal average thickness of aerogel tile in cm

  //Since Ct is in units of um^4, let's convert t to microns to compute C in um^3:
  G4double C_aerogel = Ct_aerogel / (t_aerogel*1.0e4); //This quantity has units of um^3!

  //The scattering length at a given lambda is found from L = lambda^4 / C:

  G4double lmin = 200.0, lmax = 700.0; //nm
  G4double lstep = 10.0; //nm
  const G4int nsteps = G4int( (lmax-lmin)/lstep + 0.5 ) + 1;

  G4double hbarc_eV_nm = 197.3269718; //eV nm


  G4double *Ephoton_aerogel  = new G4double [nsteps];
  G4double *Rindex_aerogel   = new G4double [nsteps];
  G4double *Rayleigh_aerogel = new G4double [nsteps];

  G4bool inrange = true;

  G4double Ephoton_633 = ( twopi * hbarc_eV_nm / 633.0 )*eV;
  //G4double nquartz_633 = ( ( QuartzWindow->GetMaterialPropertiesTable() )->GetProperty("RINDEX") )->GetValue( Ephoton_633, inrange );
  G4double nquartz_633 = 1.542599196;
  G4double nair_633 = 1.00027653;

  G4double Bratio = (1.0304-nair_633)/(nquartz_633 - nair_633);

  for(int i=0; i<nsteps; i++){
    G4double ltemp = lmin + i*lstep; //nm

    G4double nquartz = sqrt( 1.28604141 + 1.07044083*pow(ltemp/1000.0,2)/(pow(ltemp/1000.0,2)-1.00585997e-2) + 1.10202242*pow(ltemp/1000.0,2)/(pow(ltemp/1000.0,2)-100) );
    G4double nair = 1 + 5792105E-8/(238.0185-pow(ltemp/1000.0,-2)) + 167917E-8/(57.362-pow(ltemp/1000.0,-2));
	  
    int idx = nsteps - (i+1);
	  
	
    Ephoton_aerogel[idx] = (twopi * hbarc_eV_nm / ltemp )*eV;
    //Since C is in um^3, we need to compute lambda^4 in um^4. The result will be in microns, which we then convert to cm:
    Rayleigh_aerogel[idx] = (pow( ltemp/1000.0, 4 )/C_aerogel / 10000.0 )*cm; //scattering length in cm
    //From the references above, the mean refractive index at 633 nm is 1.0304;
    // Assume that the dispersion follows that of fused silica, so scale n at each wavelength by the ratio of 
    // fused silica at that wavelength to fused silica at 633 nm. 
    //G4double nquartz_temp = ( ( QuartzWindow->GetMaterialPropertiesTable() )->GetProperty("RINDEX") )->GetValue( Ephoton_aerogel[i], inrange );
    Rindex_aerogel[idx] = (1.0-Bratio)*nair + Bratio*nquartz;
  }


  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("RINDEX", Ephoton_aerogel, Rindex_aerogel, nsteps );
  MPT_temp->AddProperty("RAYLEIGH", Ephoton_aerogel, Rayleigh_aerogel, nsteps );
  //MPT_temp->AddConstProperty("ABSLENGTH", 10.0*m );

  if( fMaterialsListOpticalPhotonDisabled.find( "Aerogel" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    Aerogel->SetMaterialPropertiesTable( MPT_temp );
  }
  fMaterialsMap["Aerogel"] = Aerogel;
  
  // Grinch Quartz

  G4Material* Quartz=new G4Material("Quartz", density= 2.200*g/cm3, 2);
  Quartz->AddElement(Si, 1);
  Quartz->AddElement(O, 2);
  G4double rindex_Quartz[nEntries];
  G4double absl_Quartz[nEntries];
  for ( int i = 0; i < nEntries; ++i ) {
    nm_lambda=PhotonWaveLength[i]/nm;
    rindex_Quartz[i]=1.736-0.001396*nm_lambda+2.326e-06*nm_lambda*nm_lambda-1.281e-09*nm_lambda*nm_lambda*nm_lambda;
    absl_Quartz[i]=(-93.08+0.07151*nm_lambda+0.001799*nm_lambda*nm_lambda-1.678e-06*nm_lambda*nm_lambda*nm_lambda)*m;
    if ( absl_Quartz[i]<0 ) {
      absl_Quartz[i]=0;
    }
  }
  G4MaterialPropertiesTable* mptQuartz=new G4MaterialPropertiesTable();
  mptQuartz->AddProperty("RINDEX",PhotonEnergy,rindex_Quartz,nEntries);
  mptQuartz->AddProperty("ABSLENGTH",PhotonEnergy,absl_Quartz,nEntries);

  if( fMaterialsListOpticalPhotonDisabled.find( "Quartz" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    Quartz->SetMaterialPropertiesTable(mptQuartz);
  }
  fMaterialsMap["Quartz"] = Quartz;


  //Define the reflectivity of the mirror surface using the "Logical Skin surface":
  const G4int nentries_mirr = 15;

  G4double Ephoton_mirr[nentries_mirr] = { 1.787025848*eV, 1.851389327*eV, 2.027907076*eV, 2.206874769*eV, 2.468929557*eV, 
					   2.822981327*eV, 3.15509551*eV, 3.548346967*eV, 3.916325599*eV, 4.348919321*eV, 
					   4.496873938*eV, 4.750862572*eV, 5.176123788*eV, 5.633059855*eV, 6.178512519*eV };  

  G4double Reflectivity_mirr[nentries_mirr] = { 0.867162, 0.872027, 0.879324, 0.882973, 0.884189, 0.884189, 0.882973, 
						0.878108, 0.858649, 0.841622, 0.823378, 0.765, 0.687162, 0.619054, 0.557027 };

  // Optical surfaces (for mirrors, etc.):
  // First, define a default optical surface for general use; Aluminum with 100% reflectivity:
  MPT_temp = new G4MaterialPropertiesTable();
  //MPT_temp->AddConstProperty( "REFLECTIVITY", 1.0 );
  G4double Ephoton_mirrdefault[2] = { (1240./1000.)*eV, (1240./100.)*eV};
  G4double Reflectivity_mirrdefault[2] = { 1.0, 1.0 };
  G4OpticalSurface *DefaultOpticalSurface = new G4OpticalSurface("MirrorDefault");

  MPT_temp->AddProperty("REFLECTIVITY", Ephoton_mirrdefault, Reflectivity_mirrdefault, 2 );
  
  DefaultOpticalSurface->SetType( dielectric_metal );
  DefaultOpticalSurface->SetFinish( polished );
  DefaultOpticalSurface->SetModel( glisur );
  DefaultOpticalSurface->SetMaterialPropertiesTable(MPT_temp);

  fOpticalSurfacesMap["MirrorDefault"] = DefaultOpticalSurface;

  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty( "REFLECTIVITY", Ephoton_mirr, Reflectivity_mirr, nentries_mirr );

  G4OpticalSurface *Mirrsurf = new G4OpticalSurface("Mirrsurf");

  Mirrsurf->SetType( dielectric_metal );
  Mirrsurf->SetFinish( polished );
  Mirrsurf->SetModel( glisur );

  Mirrsurf->SetMaterialPropertiesTable( MPT_temp );

  fOpticalSurfacesMap["Mirrsurf"] = Mirrsurf;

  //G4cout << "Material properties for Mirrsurf:" << G4endl;

    
  //********************************************************************
  //************                 ECAL             **********************
  //********************************************************************  
  // Additional elements needed for the materials:
  G4Element* elK = man->FindOrBuildElement("K");
  G4Element* elAs = man->FindOrBuildElement("As");
  G4Element* elPb = man->FindOrBuildElement("Pb");
  G4Element* elBe = man->FindOrBuildElement("Be");
  G4Element* elW = man->FindOrBuildElement("W");
  
  G4Material* Titanium = man->FindOrBuildMaterial("G4_Ti");
  fMaterialsMap["Titanium"] = Titanium;
  // if( fMaterialsMap.find("Titanium") == fMaterialsMap.end() ){ 
  //   fMaterialsMap["Titanium"] = new G4Material(name="Titanium", z=22., a=47.867*g/mole, density=4.54*g/cm3);
  // }
  
  // Materials necessary to build TF1 aka lead-glass
  G4Material* PbO = new G4Material("TF1_PbO", bigden, 2);
  PbO->AddElement(elPb, 1);
  PbO->AddElement(elO, 1);
  fMaterialsMap["TF1_PbO"] = PbO;

  G4Material* K2O = new G4Material("TF1_K2O", bigden, 2);
  K2O->AddElement(elK, 2);
  K2O->AddElement(elO, 1);
  fMaterialsMap["TF1_K2O"] = K2O;

  G4Material* As2O3 = new G4Material("TF1_As2O3", bigden, 2);
  As2O3->AddElement(elAs, 2);
  As2O3->AddElement(elO, 3);
  fMaterialsMap["TF1_As2O3"] = As2O3;
  
  // Simulating annealing: http://hallaweb.jlab.org/12GeV/SuperBigBite/SBS-minutes/2014/Sergey_Abrahamyan_LGAnnealing_2014.pdf
  // const G4int nentries_annealing_model=50;

  // G4double Ephoton_annealing_model[nentries_annealing_model] = {
  //   1.513*eV, 1.531*eV, 1.548*eV, 1.569*eV, 1.590*eV,
  //   1.609*eV, 1.632*eV, 1.653*eV, 1.676*eV, 1.698*eV,
  //   1.724*eV, 1.746*eV, 1.773*eV, 1.797*eV, 1.825*eV,
  //   1.849*eV, 1.878*eV, 1.909*eV, 1.935*eV, 1.965*eV,
  //   2.001*eV, 2.033*eV, 2.069*eV, 2.100*eV, 2.135*eV,
  //   2.171*eV, 2.213*eV, 2.259*eV, 2.297*eV, 2.337*eV,
  //   2.386*eV, 2.431*eV, 2.481*eV, 2.530*eV, 2.587*eV,
  //   2.638*eV, 2.697*eV, 2.752*eV, 2.820*eV, 2.887*eV,
  //   2.954*eV, 3.024*eV, 3.102*eV, 3.179*eV, 3.261*eV,
  //   3.351*eV, 3.447*eV, 3.543*eV, 3.655*eV, 3.752*eV};

  // G4double Absorption_avg[nentries_annealing_model] = {
  //   197.05*cm, 197.62*cm, 197.05*cm, 197.01*cm, 197.48*cm, 
  //   196.86*cm, 197.39*cm, 196.74*cm, 197.25*cm, 196.61*cm, 
  //   196.55*cm, 197.08*cm, 196.45*cm, 196.40*cm, 196.94*cm, 
  //   196.88*cm, 196.85*cm, 195.69*cm, 196.25*cm, 196.21*cm, 
  //   196.16*cm, 194.94*cm, 192.57*cm, 189.58*cm, 184.89*cm, 
  //   179.64*cm, 174.91*cm, 167.96*cm, 160.95*cm, 154.04*cm, 
  //   148.19*cm, 140.62*cm, 137.08*cm, 135.03*cm, 130.38*cm, 
  //   124.69*cm, 116.54*cm, 106.33*cm,  95.16*cm,  83.84*cm, 
  //   74.26*cm,  70.30*cm,  60.21*cm,  43.32*cm,  33.15*cm, 
  //   20.07*cm,   9.22*cm,   5.19*cm,   1.73*cm,   0.58*cm}; 

  //AJRP Jan. 10, 2016: Sergey's model:

  const G4int nentries_atilde = 31;

  G4double Ephoton_atilde[nentries_atilde] =
    { 1.905921267*eV, 2.004354622*eV, 2.032730235*eV, 2.067855684*eV, 2.10114292*eV, 2.138694327*eV,
      2.174317107*eV, 2.21455858*eV, 2.25631769*eV, 2.296002918*eV, 2.34092119*eV, 2.383666502*eV,
      2.43211172*eV, 2.47828998*eV, 2.530700045*eV, 2.585380093*eV, 2.637618826*eV, 2.697071068*eV,
      2.753970493*eV, 2.818848091*eV, 2.881060974*eV, 2.946081949*eV, 3.020448926*eV, 3.091991741*eV,
      3.181036915*eV, 3.260489442*eV, 3.351822418*eV, 3.440153808*eV, 3.550739785*eV, 3.650021782*eV,
      3.75501552*eV };   

  G4double atilde[nentries_atilde] =
    { 198.544*cm, 199.029*cm, 197.573*cm, 195.146*cm, 191.748*cm, 186.893*cm, 182.524*cm, 176.699*cm, 170.388*cm, 164.078*cm,
      156.796*cm, 150*cm, 143.204*cm, 139.806*cm, 138.35*cm, 133.981*cm, 127.67*cm, 119.417*cm, 109.223*cm, 98.0583*cm,
      85.9223*cm, 75.2427*cm, 71.3592*cm, 61.6505*cm, 44.6602*cm, 34.466*cm, 20.8738*cm, 9.70874*cm, 5.82524*cm, 1.94175*cm, 0.970874*cm };

  TSpline3 *spline_atilde = new TSpline3( "spline_atilde", Ephoton_atilde, atilde, nentries_atilde ); 
  
  const G4int nentries_btilde = 33;

  G4double Ephoton_btilde[nentries_btilde] =
    { 1.90414474*eV, 1.934990559*eV, 1.964154186*eV, 1.994213564*eV, 2.03094239*eV, 2.063097508*eV, 2.099357157*eV, 2.136910476*eV,
      2.172538068*eV, 2.212780233*eV, 2.251005237*eV, 2.294239978*eV, 2.339163636*eV, 2.381921981*eV, 2.430381335*eV, 2.480858377*eV,
      2.533476625*eV, 2.579056398*eV, 2.635970367*eV, 2.70052486*eV, 2.752381143*eV, 2.817298211*eV, 2.879548933*eV, 2.950680798*eV,
      3.025408552*eV, 3.090643078*eV, 3.179772544*eV, 3.251904593*eV, 3.350726486*eV, 3.439151972*eV, 3.541106824*eV, 3.64002501*eV,
      3.744617113*eV };

  G4double btilde[nentries_btilde] =
    { 29.1935*cm, 28.871*cm, 28.0645*cm, 27.2581*cm, 25.8065*cm, 24.1935*cm, 22.2581*cm, 20.4839*cm,
      18.7097*cm, 17.0968*cm, 15.3226*cm, 14.0323*cm, 12.5806*cm, 11.6129*cm, 10.1613*cm, 9.03226*cm,
      8.3871*cm, 7.58065*cm, 6.77419*cm, 6.12903*cm, 5.96774*cm, 5.32258*cm, 4.83871*cm, 4.51613*cm,
      4.03226*cm, 4.03226*cm, 3.70968*cm, 3.54839*cm, 3.22581*cm, 3.22581*cm, 3.22581*cm, 2.41935*cm,
      2.58065*cm };

  TSpline3 *spline_btilde = new TSpline3( "spline_btilde", Ephoton_btilde, btilde, nentries_btilde );
  
  // const G4int nentries_Cz0 = 47;

  // G4double z_Cz0[nentries_Cz0] =
  //   { 0.794393*cm, 1.58879*cm, 2.38318*cm, 3.2243*cm, 3.97196*cm, 4.81308*cm, 5.60748*cm, 6.40187*cm, 7.19626*cm, 8.03738*cm,
  //     8.83178*cm, 9.62617*cm, 10.4206*cm, 11.215*cm, 12.0093*cm, 12.8037*cm, 13.5981*cm, 14.3925*cm, 15.1869*cm, 16.028*cm,
  //     16.8224*cm, 17.6168*cm, 18.4112*cm, 19.2056*cm, 20*cm, 20.7477*cm, 21.5888*cm, 22.4299*cm, 23.2243*cm, 24.0187*cm, 24.8131*cm,
  //     25.6542*cm, 26.4486*cm, 27.1963*cm, 27.9907*cm, 28.785*cm, 29.6262*cm, 30.4206*cm, 31.215*cm, 32.0093*cm, 32.8505*cm, 33.5981*cm,
  //     36.028*cm, 36.8224*cm, 37.6636*cm, 38.4579*cm, 40*cm };

  // G4double Cz0[nentries_Cz0] =
  //   { 0.759136, 0.750884, 0.72613, 0.691749, 0.649116, 0.607859, 0.565226, 0.521218, 0.482711,
  //     0.446955, 0.413949, 0.385069, 0.358939, 0.33556, 0.313556, 0.294303, 0.276424, 0.259921,
  //     0.247544, 0.232417, 0.218664, 0.206287, 0.195285, 0.185658, 0.174656, 0.166405, 0.159528,
  //     0.154028, 0.147151, 0.143026, 0.1389, 0.134774, 0.133399, 0.129273, 0.126523, 0.123772,
  //     0.121022, 0.118271, 0.115521, 0.11277, 0.111395, 0.11002, 0.104519, 0.101768, 0.101768,
  //     0.100393, 0.0990177 };

  const G4double C0_80krad = 0.76; //"Rad damage at Z = 0 at 410 nm for 80 krad radiation dose"
  // For C16, we'll define dose rates as a function of z individually by row and column
  // For ECAL, we'll define a universal dose rate profile vs z, averaged over the surface of the calorimeter.

  const G4int nbins_z_ECAL_dose = 20;
  G4double z_ECAL_dose[nbins_z_ECAL_dose] = { 1.0*cm, 3.0*cm, 5.0*cm, 7.0*cm, 9.0*cm,
					      11.0*cm, 13.0*cm, 15.0*cm, 17.0*cm, 19.0*cm,
					      21.0*cm, 23.0*cm, 25.0*cm, 27.0*cm, 29.0*cm,
					      31.0*cm, 33.0*cm, 35.0*cm, 37.0*cm, 39.0*cm };
  G4double DoseRate_ECAL_vs_z[nbins_z_ECAL_dose] =
    { 0.0987343, 0.080827, 0.0721091, 0.064972, 0.0559033, 
      0.0489642, 0.0380918, 0.0306929, 0.0245048, 0.0211042, 
      0.0155789, 0.0140177, 0.0110763, 0.00985617, 0.00788075, 
      0.00623784, 0.00566375, 0.0045522, 0.00399444, 0.00439388 }; //in krad/h

  G4double z_DoseRate_C16[17];
  G4double DoseRate_C16_vs_z[17][4][4];
  for( int i=0; i<4; i++){
    for( int j=0; j<4; j++ ){
      for( int k=0; k<16; k++ ){
	DoseRate_C16_vs_z[k][i][j] = 0.0;
      }
    }
  }
  TString line;
  ifstream infile("database/C16_doserate.txt");

  while( infile && !infile.eof() ){
    int nbinstemp;
    double ztemp, ratetemp;
    int row, col;
    if( infile >> row >> col >> nbinstemp ){
      //cout << "row, col, nbins = " << row << ", " << col << ", " << nbinstemp << endl;
      for( int bin=0; bin<nbinstemp; bin++ ){
	infile >> ztemp >> ratetemp;
	//cout << "z, rate = " << ztemp << ", " << ratetemp << endl;
	if( bin<17 ){
	  z_DoseRate_C16[bin] = ztemp*cm;
	  DoseRate_C16_vs_z[bin][row-1][col-1] = ratetemp;
	}
      }
    } else break;
  }
  
  
  TSpline3 *spline_DoseRate_ECAL = new TSpline3( "spline_DoseRate_ECAL", z_ECAL_dose, DoseRate_ECAL_vs_z, nbins_z_ECAL_dose );
  
  //TSpline3* spline_Cz0 = new TSpline3( "Rad_damage_profile_80krad", z_Cz0, Cz0, nentries_Cz0 ); 
  
  //Annealing lifetime is given in hours by:
  // tau = tau_0 * exp( -T / T0 ) (where T is absolute temperature in Kelvin).
  G4double AnnealingLifetime_T0 = 43.5; //Kelvin
  G4double AnnealingLifetime_tau0 = 1.674e5; //h

  //Temperature profile is linear with "z"
  const G4double Temp_front_C16 = 250.0; //Degrees C
  const G4double Temp_back_C16 = 190.0; //Degrees C

  const G4double Temp_front_ECAL = 225.0; //Degrees C
  const G4double Temp_back_ECAL = 185.0; //Degrees C
  
  // Absorption increase with temperature,
  // parametrized as absorption in 40 cm = Amin + C * T^2
  // absorption in 40 cm is related to absorption length by:
  // I(40)/I0 = exp( -40 cm/Labs ) = 1 - A40cm;
  // -40 cm/L = ln( 1 - A40cm );
  // L = -40 cm/ln( 1 - A40cm );
  // Ratio L/L0 = -ln( 1 - A40cm_min )/-ln
  const G4double Abs40cm_min = 0.336;
  const G4double Abs40cm_T2coeff = 6.06e-6;
  
  // Values come from old GSTAR code written by K.Shestermanov 
  const G4int nentries_ecal_QE = 37;

  G4double Ephoton_ECAL_QE[nentries_ecal_QE] = {
    1.91*eV, 1.98*eV, 2.07*eV, 2.11*eV, 2.16*eV,
    2.20*eV, 2.25*eV, 2.31*eV, 2.37*eV, 2.41*eV,
    2.52*eV, 2.60*eV, 2.71*eV, 2.77*eV, 2.84*eV,
    2.90*eV, 2.97*eV, 3.03*eV, 3.07*eV, 3.09*eV,
    3.11*eV, 3.13*eV, 3.15*eV, 3.16*eV, 3.18*eV,
    3.20*eV, 3.22*eV, 3.25*eV, 3.27*eV, 3.29*eV,
    3.31*eV, 3.33*eV, 3.36*eV, 3.42*eV, 3.54*eV,
    3.81*eV, 3.95*eV}; 

  G4double Rindex_TF1[nentries_ecal_QE] = {
    1.58546, 1.59412, 1.60393, 1.60788, 1.61250,
    1.61597, 1.62003, 1.62455, 1.62873, 1.63133,
    1.63786, 1.64208, 1.64728, 1.64985, 1.65265,
    1.65489, 1.65732, 1.65928, 1.66052, 1.66112,
    1.66171, 1.66229, 1.66286, 1.66314, 1.66369,
    1.66423, 1.66476, 1.66554, 1.66605, 1.66655,
    1.66704, 1.66752, 1.66823, 1.66958, 1.67209,
    1.67690, 1.67903};  

  G4double Abslength_TF1[nentries_ecal_QE] = {
    300.0000*cm, 300.0000*cm, 300.0000*cm, 300.0000*cm, 300.0000*cm,
    300.0000*cm, 300.0000*cm, 300.0000*cm, 293.8310*cm, 284.1470*cm,
    264.2000*cm, 236.9510*cm, 214.8510*cm, 190.8860*cm, 172.2030*cm,
    156.7000*cm, 137.0000*cm, 116.9000*cm, 102.2920*cm, 92.5750*cm,
    85.2866*cm, 80.0355*cm, 70.8671*cm, 63.0679*cm, 56.3423*cm,
    50.5331*cm, 45.2678*cm, 40.6448*cm, 36.3173*cm, 32.5037*cm,
    28.8809*cm, 25.3716*cm, 22.0407*cm, 17.8329*cm, 14.7174*cm,
    10.7150*cm,  4.88581*cm}; 

  //Density 3.86 g/cm^3
  //TF1 as defined by GSTAR code - uses arrays of size nentries_ecal_QE only
  G4Material* TF1 = new G4Material("TF1", 3.86*g/cm3, 4);
  TF1->AddMaterial(PbO, 0.512);
  TF1->AddMaterial(SiO2, 0.413);
  TF1->AddMaterial(K2O, 0.070);
  TF1->AddMaterial(As2O3, 0.005);

  //lead-glass with no optical properties:
  G4Material *TF1_dead = new G4Material("TF1_dead", 3.86 * g/cm3, 1 );
  TF1_dead->AddMaterial( TF1, 1. );

  //This is an approximation of the BigBite preshower lead-glass, which is more dense:
  //chemical composition is not correct, but higher density is needed to get the
  //preshower energy deposition right:
  G4Material *TF5 = new G4Material("TF5", 4.78*g/cm3, 4);
  TF5->AddMaterial(PbO, 0.512);
  TF5->AddMaterial(SiO2, 0.413);
  TF5->AddMaterial(K2O, 0.070);
  TF5->AddMaterial(As2O3, 0.005);
  
  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("RINDEX", Ephoton_ECAL_QE, Rindex_TF1, nentries_ecal_QE );
  //MPT_temp->AddProperty("ABSLENGTH", Ephoton_ECAL_QE, Abslength_TF1, nentries_ecal_QE );
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_atilde, atilde, nentries_atilde );

  if( fMaterialsListOpticalPhotonDisabled.find( "TF1" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    TF1->SetMaterialPropertiesTable( MPT_temp );
  }
  if( fMaterialsListOpticalPhotonDisabled.find( "TF5" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    TF5->SetMaterialPropertiesTable( MPT_temp );
  }
  fMaterialsMap["TF1"] = TF1; //Default TF1: no temperature increase, no rad. damage.
  fMaterialsMap["TF5"] = TF5;
  fMaterialsMap["TF1_dead"] = TF1_dead;
  
  //For C16 and/or ECAL, we need the following configurations of lead-glass:
  //C16: Elevated temp, no rad. damage: z-dependent temperature --> z-dependent absorption
  //C16: Elevated temp, rad. damage: z-dependent temperature and z-dependent radiation dose rate/thermal annealing rate.
  //ECAL: Default: room temp., no rad damage (only need one material property here)
  //ECAL: Equilibrium state between elevated temperature and rad. damage rate/thermal annealing rate
  //We will need to make spline interpolations of all the various curves: atilde, btilde, Cz0
  //
  const G4int Ntemp = 40;

  //G4cout << "Number of C16 segments = " << fSegmentC16 << G4endl;
  
  for( G4int segment=0; segment<fSegmentC16; segment++ ){
    //define nsegments different lead-glass material properties: 
    G4double zsegment = (segment + 0.5) * fSegmentThickC16;
    //Get temperature, dose rate at z=0 and annealing rate of the segment:
    G4double Temp_z_ECAL = Temp_front_ECAL + (Temp_back_ECAL - Temp_front_ECAL)/(G4double(fSegmentC16*fSegmentThickC16)) * zsegment;
    G4double Temp_z_C16 = Temp_front_C16 + (Temp_back_C16 - Temp_front_C16)/(G4double(fSegmentC16*fSegmentThickC16)) * zsegment;
    
    // G4double AbsIncreaseFactor_ECAL = log( 1.0 - Abs40cm_min )/log( 1.0 - (Abs40cm_min + Abs40cm_T2coeff * pow( max(Temp_z_ECAL-20.0,0.0), 2 ) ) );
    // G4double AbsIncreaseFactor_C16 = log( 1.0 - Abs40cm_min )/log( 1.0 - (Abs40cm_min + Abs40cm_T2coeff * pow( max(Temp_z_C16-20.0,0.0), 2 ) ) );
    
    G4double tau_annealing_ECAL = AnnealingLifetime_tau0 * exp( -( (Temp_z_ECAL + 273.15)/AnnealingLifetime_T0 ) ); //in hours!
    G4double tau_annealing_C16 = AnnealingLifetime_tau0 * exp( -( (Temp_z_C16 + 273.15)/AnnealingLifetime_T0 ) ); //in hours!
    
    //This computes the equilibrium radiation damage profile as a function of z!
    double zmin = spline_DoseRate_ECAL->GetXmin();
    double zmax = spline_DoseRate_ECAL->GetXmax();
    double zeval = zsegment >= zmin ? zsegment : zmin;
    zeval = zsegment <= zmax ? zsegment : zmax;

    //Re-interpret user dose rate parameter as an overall scale factor for the dose rate.
    G4double Cz_eq_ECAL = C0_80krad * spline_DoseRate_ECAL->Eval( zeval ) * fDoseRateC16 / 80.0 * tau_annealing_ECAL;
    
    //G4double Cz_eq_ECAL = spline_Cz0->Eval( zeval ) * fDoseRateC16 * tau_annealing_ECAL / 80.0; //dose rate is assumed to be given in krad/hour!
    //G4double Cz_eq_C16 = spline_Cz0->Eval( zeval ) * fDoseRateC16 * tau_annealing_C16 / 80.0; //dose rate is assume to be given in krad/hour!

    G4double Cz_eq_C16[4][4];
    for( int row=0; row<4; row++ ){
      for( int col=0; col<4; col++ ){
	G4double rate_temp[17];
	for( int bintemp=0; bintemp<17; bintemp++ ){
	  rate_temp[bintemp] = DoseRate_C16_vs_z[bintemp][row][col];
	}
	TSpline3 stemp( "splinetemp", z_DoseRate_C16, rate_temp, 17 );
	zeval = zsegment >= z_DoseRate_C16[0] ? zsegment : z_DoseRate_C16[0];
	zeval = zsegment <= z_DoseRate_C16[16] ? zsegment : z_DoseRate_C16[16];
	
	Cz_eq_C16[row][col] = C0_80krad * stemp.Eval( zeval ) * fDoseRateC16 / 80.0 * tau_annealing_C16;

	// G4cout << "Row, Col, zeval, doserate, Ceq(z) = " << row+1 << ", " << col+1 << ", " << zeval/cm
	//        << ", " << stemp.Eval( zeval ) * fDoseRateC16 << ", " << Cz_eq_C16[row][col] << G4endl;
      }
    }
    
    G4double Ephot_min = 1.91*eV;
    G4double Ephot_max = 3.74*eV;

    G4double Ephoton_abslength[Ntemp+1];
    G4double abslength_ECAL[Ntemp+1];
    G4double abslength_C16[Ntemp+1][4][4];

    //   G4cout << "Ntemp = " << Ntemp << G4endl;
    
    for( G4int iE=0; iE<=Ntemp; iE++ ){
      Ephoton_abslength[iE] = Ephot_min + (Ephot_max-Ephot_min)/G4double(Ntemp) * iE;
      abslength_ECAL[iE] = 1.0 / (1.0/spline_atilde->Eval( Ephoton_abslength[iE] ) + Cz_eq_ECAL / spline_btilde->Eval( Ephoton_abslength[iE] ) );
      //abslength_C16[iE] = 1.0 / (1.0/spline_atilde->Eval( Ephoton_abslength[iE] ) + Cz_eq_C16 / spline_btilde->Eval( Ephoton_abslength[iE] ) );

      for( int row=0; row<4; row++ ){
	for( int col=0; col<4; col++ ){
	  abslength_C16[iE][row][col] = 1.0 / (1.0/spline_atilde->Eval( Ephoton_abslength[iE] ) + Cz_eq_C16[row][col] / spline_btilde->Eval( Ephoton_abslength[iE] ) );
	}
      }
      
      // G4cout << "z, Ephoton, lambda, Labs(ECAL), Labs(C16), atilde = " << zsegment/cm << ", "
      // 	     << Ephoton_abslength[iE] / eV << ", "
      // 	     << twopi * hbarc / Ephoton_abslength[iE] / nm << ", "
      // 	     << abslength_ECAL[iE]/cm << ", " << abslength_C16[iE]/cm << ", " << spline_atilde->Eval( Ephoton_abslength[iE] )/cm << G4endl;
    }

    spline_atilde->Delete();
    spline_btilde->Delete();
    spline_DoseRate_ECAL->Delete();
    
    //Next: define new materials!
    TString matname;
    matname.Form( "TF1_anneal_ECAL_z%d", segment );
    G4Material *mat_temp = new G4Material( matname.Data(), 3.86*g/cm3, 1 );
    mat_temp->AddMaterial( TF1, 1.0 );

    MPT_temp = new G4MaterialPropertiesTable();
    MPT_temp->AddProperty("RINDEX", Ephoton_ECAL_QE, Rindex_TF1, nentries_ecal_QE );
    MPT_temp->AddProperty("ABSLENGTH", Ephoton_abslength, abslength_ECAL, Ntemp+1 );

    if( fMaterialsListOpticalPhotonDisabled.find( matname.Data() ) == fMaterialsListOpticalPhotonDisabled.end() ){
      mat_temp->SetMaterialPropertiesTable( MPT_temp );
    }
    fMaterialsMap[matname.Data()] = mat_temp;

    for( int row=0; row<4; row++ ){
      for( int col=0; col<4; col++ ){
	
	matname.Form( "TF1_anneal_C16_row%d_col%d_z%d", row+1, col+1, segment );
	mat_temp = new G4Material( matname.Data(), 3.86*g/cm3, 1 );
	mat_temp->AddMaterial( TF1, 1.0 );

	MPT_temp = new G4MaterialPropertiesTable();
	MPT_temp->AddProperty("RINDEX", Ephoton_ECAL_QE, Rindex_TF1, nentries_ecal_QE );

	G4double abslength_temp[Ntemp+1];
	for( int iE=0; iE<Ntemp+1; iE++ ){
	  abslength_temp[iE] = abslength_C16[iE][row][col];
	  // G4cout << "Row " << row+1 << " Column " << col+1 << " Ephoton = " << Ephoton_abslength[iE]/eV
	  // 	 << " eV, Abs. length = " << abslength_temp[iE]/cm << " cm, atilde = "
	  // 	 << spline_atilde->Eval( Ephoton_abslength[iE] )/cm << G4endl;
	}
	
	MPT_temp->AddProperty("ABSLENGTH", Ephoton_abslength, abslength_temp, Ntemp+1 );

	if( fMaterialsListOpticalPhotonDisabled.find( matname.Data() ) == fMaterialsListOpticalPhotonDisabled.end() ){
	  mat_temp->SetMaterialPropertiesTable( MPT_temp );
	}
	fMaterialsMap[matname.Data()] = mat_temp;
      }
    } 
  }

  //spline
  
  //****  TF1 implementing annealing model  ****
  // G4Material* TF1_anneal = new G4Material("TF1_anneal", 3.86*g/cm3, 4);
  // TF1_anneal->AddMaterial(PbO, 0.512);
  // TF1_anneal->AddMaterial(SiO2, 0.413);
  // TF1_anneal->AddMaterial(K2O, 0.070);
  // TF1_anneal->AddMaterial(As2O3, 0.005);

  // MPT_temp = new G4MaterialPropertiesTable();
  // MPT_temp->AddProperty("RINDEX", Ephoton_ECAL_QE, Rindex_TF1, nentries_ecal_QE );
  // MPT_temp->AddProperty("ABSLENGTH", Ephoton_annealing_model, Absorption_avg, nentries_annealing_model );

  // TF1_anneal->SetMaterialPropertiesTable( MPT_temp );
  // fMaterialsMap["TF1_anneal"] = TF1_anneal;

  G4Material *Mylar = man->FindOrBuildMaterial("G4_MYLAR");
  fMaterialsMap["Mylar"] = Mylar;

  //Photocathode & window material/optical properties!
  G4double den_pmtwindow_ecal = 2.20*g/cm3;

  G4Material *QuartzWindow_ECal = new G4Material( "QuartzWindow_ECal", den_pmtwindow_ecal, nel=2 );
  QuartzWindow_ECal->AddElement( O, natoms=2 );
  QuartzWindow_ECal->AddElement( Si, natoms=1 );

  //Rindex comes from GSTAR code 
  G4double Rindex_quartz_ecal[nentries_ecal_QE] = { 
    1.41935, 1.42989, 1.44182, 1.44661, 1.45221,
    1.45640, 1.46131, 1.46677, 1.47181, 1.47495,
    1.48280, 1.48788, 1.49412, 1.49720, 1.50055,
    1.50322, 1.50614, 1.50847, 1.50995, 1.51067,
    1.51137, 1.51206, 1.51274, 1.51307, 1.51373,
    1.51437, 1.51501, 1.51593, 1.51654, 1.51713,
    1.51771, 1.51828, 1.51912, 1.52073, 1.52370,
    1.52938, 1.53188};

  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("RINDEX", Ephoton_ECAL_QE, Rindex_quartz_ecal, nentries_ecal_QE);
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_ECAL_QE, Abslength_TF1, nentries_ecal_QE); //Do we need this?? 

  if( fMaterialsListOpticalPhotonDisabled.find( "QuartzWindow_ECal" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    QuartzWindow_ECal->SetMaterialPropertiesTable( MPT_temp );
  }
  fMaterialsMap["QuartzWindow_ECal"] = QuartzWindow_ECal;

  //Photocathode material - RIndex same as QuartzWindow_ECal, define a QE from GSTAR code
  G4double PMT_ECAL_QE[nentries_ecal_QE] = {
    0.004, 0.016, 0.031, 0.043, 0.052,
    0.059, 0.069, 0.083, 0.100, 0.110,
    0.131, 0.150, 0.171, 0.184, 0.191,
    0.196, 0.200, 0.198, 0.195, 0.193,
    0.190, 0.188, 0.185, 0.183, 0.180,
    0.177, 0.173, 0.169, 0.163, 0.160,
    0.155, 0.150, 0.139, 0.123, 0.084,
    0.023, 0.003 }; 

  G4double den_photocathode_ecal = 2.57*g/cm3;
  G4Material *Photocathode_material_ecal = new G4Material( "Photocathode_material_ecal", den_photocathode_ecal, nel=3 );
  Photocathode_material_ecal->AddElement( Sb, natoms=1 );
  Photocathode_material_ecal->AddElement( K, natoms=2 );
  Photocathode_material_ecal->AddElement( Cs, natoms=1 );

  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("EFFICIENCY", Ephoton_ECAL_QE, PMT_ECAL_QE, nentries_ecal_QE ); 
  MPT_temp->AddProperty("RINDEX", Ephoton_ECAL_QE, Rindex_quartz_ecal, nentries_ecal_QE );

  if( fMaterialsListOpticalPhotonDisabled.find( "Photocathode_material_ecal" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    Photocathode_material_ecal->SetMaterialPropertiesTable( MPT_temp );
  }
  fMaterialsMap["Photocathode_material_ecal"] = Photocathode_material_ecal;

  //Define optical properties for AIR:
  //G4Material *Special_Air = man->FindOrBuildMaterial("G4_AIR");
  G4double airdensity = 1.20479*kg/m3;

  G4Material *Special_Air = new G4Material( "Special_Air", airdensity, 4 );
  Special_Air->AddElement( elC, fractionmass = 0.000124 );
  Special_Air->AddElement( elN, fractionmass = 0.755268 );
  Special_Air->AddElement( elO, fractionmass = 0.231781 );
  Special_Air->AddElement( elAr, fractionmass = 0.012827 );
  
  MPT_temp = new G4MaterialPropertiesTable();

  G4double Rindex_air[nentries_ecal_QE] = {
    1.0003, 1.0003, 1.0003, 1.0003, 1.0003,
    1.0003, 1.0003, 1.0003, 1.0003, 1.0003,
    1.0003, 1.0003, 1.0003, 1.0003, 1.0003,
    1.0003, 1.0003, 1.0003, 1.0003, 1.0003,
    1.0003, 1.0003, 1.0003, 1.0003, 1.0003,
    1.0003, 1.0003, 1.0003, 1.0003, 1.0003,
    1.0003, 1.0003, 1.0003, 1.0003, 1.0003,
    1.0003, 1.0003};  

  G4double Abslength_air[nentries_ecal_QE] = {
    1000.0*cm, 1000.0*cm, 1000.0*cm, 1000.0*cm, 1000.0*cm,
    1000.0*cm, 1000.0*cm, 1000.0*cm, 1000.0*cm, 1000.0*cm,
    1000.0*cm, 1000.0*cm, 1000.0*cm, 1000.0*cm, 1000.0*cm,
    1000.0*cm, 1000.0*cm, 1000.0*cm, 1000.0*cm, 1000.0*cm,
    1000.0*cm, 1000.0*cm, 1000.0*cm, 1000.0*cm, 1000.0*cm,
    1000.0*cm, 1000.0*cm, 1000.0*cm, 1000.0*cm, 1000.0*cm,
    1000.0*cm, 1000.0*cm, 1000.0*cm, 1000.0*cm, 1000.0*cm,
    1000.0*cm, 1000.0*cm};

  MPT_temp->AddProperty("RINDEX", Ephoton_ECAL_QE, Rindex_air, nentries_ecal_QE );
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_ECAL_QE, Abslength_air, nentries_ecal_QE );

  if( fMaterialsListOpticalPhotonDisabled.find( "Special_Air" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    Special_Air->SetMaterialPropertiesTable( MPT_temp );
  }
  fMaterialsMap["Special_Air"] = Special_Air;

  //Poly "filter"
  G4Material *Polyethylene = man->FindOrBuildMaterial("G4_POLYETHYLENE");
  fMaterialsMap["Polyethylene"] = Polyethylene;


  //   ************************
  //   *          CDet        *
  //   ************************
  G4Material *PLASTIC_SC_VINYLTOLUENE = man->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    

  const G4int nentries_CDET_scint = 46;

  G4double Ephoton_CDET[nentries_CDET_scint] = {
    2.40675006*eV, 2.445105532*eV, 2.488158883*eV, 2.522044473*eV, 2.555653572*eV,
    2.588921153*eV, 2.603982726*eV, 2.6346434*eV, 2.656802618*eV, 2.667359638*eV,
    2.696833226*eV, 2.709078112*eV, 2.722809302*eV, 2.745084471*eV, 2.762023812*eV,
    2.795087475*eV, 2.806774486*eV, 2.81855964*eV, 2.831937605*eV, 2.843935469*eV,
    2.863654459*eV, 2.888297428*eV, 2.905489513*eV, 2.921289742*eV, 2.938878028*eV,
    2.953425474*eV, 2.969759939*eV, 2.976340141*eV, 2.98627609*eV, 2.996271359*eV,
    3.011394009*eV, 3.019858358*eV, 3.030087486*eV, 3.038657397*eV, 3.050740039*eV,
    3.06117225*eV, 3.073442588*eV, 3.084030912*eV, 3.107225739*eV, 3.12531095*eV,
    3.154686087*eV, 3.195988284*eV, 3.234491656*eV, 3.265965379*eV, 3.326666918*eV,
    3.413361452*eV };

  G4double RelativeYield_CDET[nentries_CDET_scint] = {
    0.0300926, 0.0486111, 0.0752315, 0.0983796, 0.131944,
    0.168981, 0.190972, 0.234954, 0.283565, 0.313657,
    0.380787, 0.413194, 0.446759, 0.491898, 0.525463,
    0.599537, 0.645833, 0.738426, 0.840278, 0.90625,
    0.961806, 0.990741, 0.974537, 0.943287, 0.893519,
    0.861111, 0.805556, 0.77662, 0.721065, 0.678241,
    0.597222, 0.549769, 0.491898, 0.456019, 0.398148,
    0.356481, 0.297454, 0.25463, 0.201389, 0.173611,
    0.136574, 0.0983796, 0.0787037, 0.0648148, 0.0451389,
    0.03125 };

  G4double Ephoton_rindex_CDET[2] = { 2.40675006*eV, 3.413361452*eV };
  G4double Rindex_CDET[2] = { 1.58, 1.58 };
  G4double AbsLength_CDET[2] = {3.80*m, 3.80*m};
    
  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("SCINTILLATIONCOMPONENT1", Ephoton_CDET, RelativeYield_CDET, nentries_CDET_scint );
  MPT_temp->AddConstProperty("SCINTILLATIONTIMECONSTANT1", 2.1*ns);
  MPT_temp->AddConstProperty("SCINTILLATIONYIELD", 0.64*17400.0/MeV);
  MPT_temp->AddConstProperty("RESOLUTIONSCALE", 1.0 );
  MPT_temp->AddProperty("RINDEX", Ephoton_rindex_CDET, Rindex_CDET, 2 );
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_rindex_CDET, AbsLength_CDET, 2 );

  if( fMaterialsListOpticalPhotonDisabled.find( "CDET_BC408" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    PLASTIC_SC_VINYLTOLUENE->SetMaterialPropertiesTable(MPT_temp);
  }
    
  fMaterialsMap["CDET_BC408"] = PLASTIC_SC_VINYLTOLUENE;
    
  //Specs come from Hamamatsu Datasheet H8711 maPMT
  const G4int nentries_CDet = 37;
  G4double den_photocathode_CDet = 2.57*g/cm3;
  G4Material *Photocathode_CDet = new G4Material( "Photocathode_CDet", den_photocathode_CDet, nel=3 );
  Photocathode_CDet->AddElement( Sb, natoms=1 );
  Photocathode_CDet->AddElement( K, natoms=2 );
  Photocathode_CDet->AddElement( Cs, natoms=1 );

  G4double EPhoton_CDet[nentries_CDet] = {
    1.71482*eV, 1.73907*eV, 1.75771*eV, 1.77037*eV, 1.77676*eV, 
    1.78644*eV, 1.79295*eV, 1.80281*eV, 1.80612*eV, 1.82963*eV, 
    1.85028*eV, 1.85727*eV, 1.87497*eV, 1.87855*eV, 1.88575*eV, 
    1.89301*eV, 1.89666*eV, 1.94158*eV, 1.95703*eV, 1.98067*eV, 
    2.00081*eV, 2.02553*eV, 2.05515*eV, 2.10349*eV, 2.21728*eV, 
    2.30563*eV, 2.38385*eV, 2.47376*eV, 2.57073*eV, 2.79735*eV, 
    3.09675*eV, 3.44364*eV, 3.81785*eV, 4.07062*eV, 4.26474*eV, 
    4.43789*eV, 4.58261*eV }; 
  G4double PMT_CDet_QE[nentries_CDet] = {
    1.773e-05, 3.434e-05, 5.187e-05, 6.984e-05, 9.368e-05, 
    0.0001138, 0.0001321, 0.0001570, 0.0001759, 0.0003521, 
    0.0005321, 0.0007219, 0.0009634, 0.0012201, 0.0014003, 
    0.0016252, 0.0019037, 0.0037228, 0.0055455, 0.0075018, 
    0.0099054, 0.0119881, 0.0158992, 0.0215097, 0.0437984, 
    0.0688255, 0.0940591, 0.1233840, 0.1532860, 0.1928390, 
    0.2018930, 0.2076390, 0.1883090, 0.1606140, 0.1258940, 
    0.0886464, 0.0468604 }; 
  //Reused from ECal
  G4double Rindex_CDet[nentries_CDet] = { 
    1.41935, 1.42989, 1.44182, 1.44661, 1.45221,
    1.45640, 1.46131, 1.46677, 1.47181, 1.47495,
    1.48280, 1.48788, 1.49412, 1.49720, 1.50055,
    1.50322, 1.50614, 1.50847, 1.50995, 1.51067,
    1.51137, 1.51206, 1.51274, 1.51307, 1.51373,
    1.51437, 1.51501, 1.51593, 1.51654, 1.51713,
    1.51771, 1.51828, 1.51912, 1.52073, 1.52370,
    1.52938, 1.53188 };

  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("EFFICIENCY", EPhoton_CDet, PMT_CDet_QE, nentries_CDet);
  MPT_temp->AddProperty("RINDEX", EPhoton_CDet, Rindex_CDet, nentries_CDet);

  if( fMaterialsListOpticalPhotonDisabled.find( "Photocathode_CDet" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    Photocathode_CDet->SetMaterialPropertiesTable( MPT_temp );
  }
  fMaterialsMap["Photocathode_CDet"] = Photocathode_CDet;

  ///////////// TO DO: Make WLSFiber with actual cladding!!! //////////////////////////////

  G4Material *BCF_92 = man->FindOrBuildMaterial( "G4_POLYSTYRENE" );

  //Set up optical properties for BCF_92:
      
  MPT_temp = new G4MaterialPropertiesTable();

  G4double Ephoton_BCF92_rindex[2] = {2.0*eV, 3.5*eV};
  G4double Rindex_BCF_92[2] = {1.60, 1.60};
  G4double Abslength_BCF_92[2] = {3.5*m, 3.5*m};
  MPT_temp->AddProperty( "RINDEX", Ephoton_BCF92_rindex, Rindex_BCF_92, 2 );
  MPT_temp->AddProperty( "ABSLENGTH", Ephoton_BCF92_rindex, Abslength_BCF_92, 2 );
  MPT_temp->AddConstProperty( "WLSTIMECONSTANT", 2.7*ns );
    
  const G4int nentries_BCF92_abs = 27;

  G4double Ephoton_BCF92_abs[nentries_BCF92_abs] = {
    2.630707893*eV, 2.664734002*eV, 2.686117734*eV, 2.71333491*eV, 2.738300242*eV,
    2.75709789*eV, 2.777106798*eV, 2.803271443*eV, 2.82594757*eV, 2.844966538*eV,
    2.866276065*eV, 2.906722434*eV, 2.950487623*eV, 3.015832968*eV, 3.050184651*eV,
    3.072345879*eV, 3.099605268*eV, 3.110406153*eV, 3.144492904*eV, 3.176817825*eV,
    3.20982255*eV, 3.217527656*eV, 3.259302965*eV, 3.303523516*eV, 3.325399923*eV,
    3.391382895*eV, 3.443709311*eV };

  G4double WLSabslength_BCF92[nentries_BCF92_abs] = {
    10.28081206*cm, 3.185612619*cm, 1.803324028*cm, 0.865410817*cm, 0.508661494*cm,
    0.317450506*cm, 0.21767341*cm, 0.151335806*cm, 0.105467803*cm, 0.086143145*cm,
    0.072386575*cm, 0.054623942*cm, 0.041499166*cm, 0.021491391*cm, 0.011369803*cm,
    0.013554007*cm, 0.02312288*cm, 0.025719555*cm, 0.032357641*cm, 0.039288583*cm,
    0.048017391*cm, 0.050098717*cm, 0.065481686*cm, 0.085430618*cm, 0.095633564*cm,
    0.133619777*cm, 0.176473376*cm };

  const G4int nentries_BCF92_emission = 32;

  G4double Ephoton_BCF92_emission[nentries_BCF92_emission] = {
    2.673600455*eV, 2.650668218*eV, 2.628984229*eV, 2.620410753*eV, 2.606807804*eV,
    2.600056427*eV, 2.595017157*eV, 2.58250403*eV, 2.575883146*eV, 2.568471124*eV,
    2.552964507*eV, 2.536839201*eV, 2.52408797*eV, 2.498190812*eV, 2.479684214*eV,
    2.45918933*eV, 2.439035253*eV, 2.429443326*eV, 2.402566631*eV, 2.375567827*eV,
    2.338914317*eV, 2.310000796*eV, 2.279201125*eV, 2.252992604*eV, 2.233575468*eV,
    2.218774127*eV, 2.19813402*eV, 2.180237337*eV, 2.154512292*eV, 2.12207618*eV,
    2.087884362*eV, 2.069597591*eV };
  
  G4double BCF92_emission_relative[nentries_BCF92_emission] = {
    0, 0.042618, 0.129376, 0.196347, 0.394216,
    0.506849, 0.598174, 0.715373, 0.805175, 0.890411,
    0.968037, 0.99239, 0.993912, 0.987823, 0.971081,
    0.929985, 0.855403, 0.803653, 0.703196, 0.596651,
    0.48554, 0.392694, 0.316591, 0.266362, 0.2207,
    0.190259, 0.165906, 0.136986, 0.112633, 0.0791476,
    0.0669711, 0.0502283 };

  //this needs to be sorted in increasing order of energy; since the above is sorted in descending order,
  //we don't need any fancy sort algorithm, just reverse the order:
  G4double Ephoton_BCF92_emission_sorted[nentries_BCF92_emission];
  G4double BCF92_emission_relative_sorted[nentries_BCF92_emission];
  for( int i=0; i<nentries_BCF92_emission; i++ ){
    Ephoton_BCF92_emission_sorted[i] = Ephoton_BCF92_emission[nentries_BCF92_emission-i-1];
    BCF92_emission_relative_sorted[i] = BCF92_emission_relative[nentries_BCF92_emission-i-1];
  }
    
  MPT_temp->AddProperty("WLSABSLENGTH", Ephoton_BCF92_abs, WLSabslength_BCF92, nentries_BCF92_abs );
  MPT_temp->AddProperty("WLSCOMPONENT", Ephoton_BCF92_emission_sorted, BCF92_emission_relative_sorted, nentries_BCF92_emission );

  if( fMaterialsListOpticalPhotonDisabled.find( "BCF_92" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    BCF_92->SetMaterialPropertiesTable( MPT_temp );
  }
    
  fMaterialsMap["BCF_92"] = BCF_92;

  G4double cdet_acrylic_density = 1.19*g/cm3;
    
  G4Material *CDET_Acrylic = new G4Material( "CDET_Acrylic", cdet_acrylic_density, 3 );
  CDET_Acrylic->AddElement( elH, fractionmass = 0.080538 );
  CDET_Acrylic->AddElement( elC, fractionmass = 0.599848 );
  CDET_Acrylic->AddElement( elO, fractionmass = 0.319614 );

  MPT_temp = new G4MaterialPropertiesTable();

  G4double ephoton_acrylic[2] = {2.0*eV, 3.5*eV};
  G4double Rindex_acrylic[2] = {1.49, 1.49 };
  G4double abslength_acrylic[2] = {3.5*m, 3.5*m};
    
  MPT_temp->AddProperty("RINDEX", ephoton_acrylic, Rindex_acrylic, 2 );
  MPT_temp->AddProperty("ABSLENGTH", ephoton_acrylic, abslength_acrylic, 2 );

  if( fMaterialsListOpticalPhotonDisabled.find( "CDET_Acrylic" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    CDET_Acrylic->SetMaterialPropertiesTable(MPT_temp);
  }

  fMaterialsMap["CDET_Acrylic"] = CDET_Acrylic;

  //  *****************************
  //  *    BB timing hodoscope    *
  //  *****************************
  G4double d_PolyVinylToluene = 0.57*g/cm3;
  G4double d_Anthracene = 1.28*g/cm3;
  G4double d_BBHodo_Scinti = 1.023*g/cm3;
  G4Material* PolyVinylToluene = new G4Material( "PolyVinylToluene", d_PolyVinylToluene, 2 );
  PolyVinylToluene->AddElement(elC, fractionmass = 0.91471);
  PolyVinylToluene->AddElement(elH, fractionmass = 0.08529);
  
  G4Material* Anthracene = new G4Material( "Anthracene", d_Anthracene, 2 );
  Anthracene->AddElement(elC, fractionmass = 0.943447);
  Anthracene->AddElement(elH, fractionmass = 0.056553);
  
  G4Material *BBHodo_Scinti = new G4Material( "BBHodo_Scinti", d_BBHodo_Scinti,  2 );
  BBHodo_Scinti->AddMaterial(PolyVinylToluene, fractionmass = 0.36);
  BBHodo_Scinti->AddMaterial(Anthracene, fractionmass = 0.64);
  fMaterialsMap["BBHodo_Scinti"] = BBHodo_Scinti;
  
  //  *****************************
  //  *     Preshower/Shower      * 
  //  *****************************
  //Photonis XP3312B - http://www.qsl.net/k0ff/01%20Manuals/PMT/Photonis/XP3312B.pdf
  const G4int nentries_BB = 35;
  G4double den_photocathode_BB = 2.57*g/cm3;
  G4Material *Photocathode_BB = new G4Material( "Photocathode_BB", den_photocathode_BB, nel=3 );
  Photocathode_BB->AddElement( Sb, natoms=1 );
  Photocathode_BB->AddElement( K, natoms=2 );
  Photocathode_BB->AddElement( Cs, natoms=1 );

  G4double EPhoton_BB[nentries_BB] = {
    1.92794*eV, 1.97312*eV, 1.99651*eV, 2.01562*eV, 2.02532*eV, 
    2.03511*eV, 2.04499*eV, 2.05247*eV, 2.06000*eV, 2.06759*eV, 
    2.13847*eV, 2.20281*eV, 2.25290*eV, 2.34709*eV, 2.41439*eV, 
    2.48200*eV, 2.63747*eV, 2.75843*eV, 3.09276*eV, 3.22296*eV, 
    3.33795*eV, 3.51930*eV, 3.61758*eV, 3.71329*eV, 3.86675*eV, 
    3.94834*eV, 4.06263*eV, 4.14257*eV, 4.17337*eV, 4.18373*eV, 
    4.20463*eV, 4.23636*eV, 4.26858*eV, 4.30128*eV, 4.33448*eV };
    
  G4double PMT_BB_QE[nentries_BB] = {
    0.001987, 0.003958, 0.00600, 0.00811, 0.01017, 
    0.012283, 0.014383, 0.01644, 0.01836, 0.02083, 
    0.042897, 0.065273, 0.09065, 0.11700, 0.14572, 
    0.174573, 0.211272, 0.24593, 0.30927, 0.28516, 
    0.267383, 0.245646, 0.21668, 0.18511, 0.15323, 
    0.117895, 0.080874, 0.03693, 0.03343, 0.02942, 
    0.025184, 0.021118, 0.01717, 0.01294, 0.00862 };
  //Reused from ECal
  G4double Rindex_BB[nentries_BB] = { 
    1.41935, 1.42989, 1.44182, 1.44661, 1.45221,
    1.45640, 1.46131, 1.46677, 1.47181, 1.47495,
    1.48280, 1.48788, 1.49412, 1.49720, 1.50055,
    1.50322, 1.50614, 1.50847, 1.50995, 1.51067,
    1.51137, 1.51206, 1.51274, 1.51307, 1.51373,
    1.51437, 1.51501, 1.51593, 1.51654, 1.51713,
    1.51771, 1.51828, 1.51912, 1.52073, 1.52370 };

  MPT_temp = new G4MaterialPropertiesTable();
  MPT_temp->AddProperty("EFFICIENCY", EPhoton_BB, PMT_BB_QE, nentries_BB ); 
  MPT_temp->AddProperty("RINDEX", EPhoton_BB, Rindex_BB, nentries_BB );

  if( fMaterialsListOpticalPhotonDisabled.find( "Photocathode_BB" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    Photocathode_BB->SetMaterialPropertiesTable( MPT_temp );
  }
  fMaterialsMap["Photocathode_BB"] = Photocathode_BB;

  //   ************************
  //   *          HCAL        *
  //   ************************
  // Everything has been "taken" from Vahe's HCaloMaterials.cc code 
  //Elements:
  //G4Element *Gd = new G4Element("Gadolinium", "Gd" , z=64.0 , a=157.25*g/mole);
  G4Element *Gd = man->FindOrBuildElement( "Gd" );
  
  //- EJ-232
  G4Material* EJ232 = new G4Material("EJ232",density=1.02*g/cm3, nel=2,kStateSolid,293.15*kelvin,1.0*atmosphere);
  EJ232->AddElement(H , natoms=11);
  EJ232->AddElement(C , natoms=10);

  //for G4 version 11, everything here needs to be re-ordered in increasing order of photon energy, unfortunately:
  
  // -- EJ232 optical properties 
  const G4int nEntriesEJ232 = 66;
  G4double PhotonWaveLength_HC[nEntriesEJ232] = 
    {
      235, 240, 245, 250, 255, 260, 265, 270, 275, 280,
      285, 290, 295, 300, 305, 310, 315, 320, 325, 330,
      335, 340, 345, 350, 355, 360, 365, 370, 375, 380,
      385, 390, 395, 400, 405, 410, 415, 420, 425, 430,
      435, 440, 445, 450, 455, 460, 465, 470, 475, 480,
      485, 490, 495, 500, 505, 510, 515, 520, 525, 530,
      535, 540, 545, 550, 555, 560 
    };

  G4double ABSL_EJ232[nEntriesEJ232] = 
    { 
      3.8000e-03, 3.6000e-03, 3.3000e-03, 3.1000e-03, 2.7000e-03, 2.4000e-03, 2.0000e-03, 1.6000e-03, 1.4000e-03, 1.1000e-03,
      1.0000e-03, 9.0000e-04, 8.0000e-04, 6.9729e-04, 6.9946e-04, 6.9972e-04, 6.9983e-04, 7.9983e-04, 9.9978e-04, 1.1973e-02,
      1.7947e-02, 2.9874e-02, 6.1538e-02, 2.2268e-01, 4.6133e+00, 8.9015e+00, 1.4659e+01, 4.2145e+01, 1.2060e+02, 1.5600e+02,
      2.2000e+02, 3.4100e+02, 5.8400e+02, 1.1200e+03, 1.6000e+03, 1.8700e+03, 2.0000e+03, 2.8600e+03, 3.5300e+03, 3.8900e+03,
      4.0500e+03, 4.1400e+03, 4.2000e+03, 4.2200e+03, 4.2400e+03, 4.2600e+03, 4.2800e+03, 4.3000e+03, 4.3300e+03, 4.3500e+03,
      4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03,
      4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03
    };

  G4double FAST_EJ232[nEntriesEJ232] = 
    {
      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00,
      0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.01,
      0.02, 0.05, 0.12, 0.36, 0.52, 0.52, 0.71, 0.92, 1.00, 0.88,
      0.83, 0.84, 0.77, 0.63, 0.53, 0.45, 0.39, 0.34, 0.30, 0.24,
      0.19, 0.14, 0.11, 0.08, 0.06, 0.05, 0.04, 0.03, 0.02, 0.02,
      0.02, 0.02, 0.02, 0.01, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00,
      0.00, 0.00, 0.00, 0.00, 0.00,0.00
    };

  G4double ABSL_EJ232_sorted[nEntriesEJ232];
  G4double FAST_EJ232_sorted[nEntriesEJ232];
  
  G4double PhotonEnergyEJ232[nEntriesEJ232];
  G4double RefractiveIndexEJ232[nEntriesEJ232];
  for(int ii = 0; ii < nEntriesEJ232; ii++) {
    PhotonEnergyEJ232[ii] = 1240.0/(PhotonWaveLength_HC[nEntriesEJ232-ii-1]) *eV;
    RefractiveIndexEJ232[ii] = 1.58;
    ABSL_EJ232_sorted[ii]           = ABSL_EJ232[nEntriesEJ232-ii-1]*mm;
    FAST_EJ232_sorted[ii] = FAST_EJ232[nEntriesEJ232-ii-1];
  }

  G4MaterialPropertiesTable* MPT_EJ232 = new G4MaterialPropertiesTable();
  MPT_EJ232->AddProperty("RINDEX"       , PhotonEnergyEJ232 , RefractiveIndexEJ232 , nEntriesEJ232);
  MPT_EJ232->AddProperty("SCINTILLATIONCOMPONENT1", PhotonEnergyEJ232 , FAST_EJ232_sorted           , nEntriesEJ232);
  MPT_EJ232->AddProperty("ABSLENGTH",     PhotonEnergyEJ232 , ABSL_EJ232_sorted           , nEntriesEJ232);
  MPT_EJ232->AddConstProperty("SCINTILLATIONYIELD", (8400.0/MeV * 5.51591522788642652e-01 * 0.5/1.4)/MeV);
  MPT_EJ232->AddConstProperty("RESOLUTIONSCALE", 1.0);
  MPT_EJ232->AddConstProperty("SCINTILLATIONTIMECONSTANT1",1.40*ns);
  //MPT_EJ232->AddConstProperty("SLOWTIMECONSTANT",1.40*ns);
  //MPT_EJ232->AddConstProperty("YIELDRATIO",1.0);
  if( fMaterialsListOpticalPhotonDisabled.find( "EJ232" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    EJ232->SetMaterialPropertiesTable(MPT_EJ232);
  }
  fMaterialsMap["EJ232"] = EJ232;
    
  //- EJ-280 Wave length shifter blue to green
  G4double sigalpha;
  G4Material* BC484 = new G4Material("BC484",density=1.023*g/cm3, nel=2,kStateSolid,293.15*kelvin,1.0*atmosphere);
  BC484->AddElement(H , natoms=11);
  BC484->AddElement(C , natoms=10);
  //BC484 optical
  const G4int nEntriesBC484 = 66;

  G4double AbsWLSfiberBC484[nEntriesBC484] =
    {
      1.0000e-05, 1.0000e-05, 1.0000e-05, 1.0000e-05, 1.0000e-05, 1.0000e-05, 1.0000e-05, 1.0000e-05, 1.0000e-05, 1.0000e-05,
      1.0000e-05, 1.0000e-05, 1.0000e-05, 1.8000e-01, 1.1912e-02, 7.4242e-03, 5.3838e-03, 4.1792e-03, 3.2238e-03, 2.5249e-03,
      1.9316e-03, 1.6394e-03, 1.3685e-03, 1.1648e-03, 9.7507e-04, 9.6058e-04, 9.5387e-04, 8.8492e-04, 8.0374e-04, 9.0456e-04,
      1.1453e-03, 1.3248e-03, 1.3103e-03, 1.3398e-03, 2.0575e-03, 5.3914e-03, 2.7426e-02, 8.8128e-02, 3.5300e+03, 3.8900e+03,
      4.0500e+03, 4.1400e+03, 4.2000e+03, 4.2200e+03, 4.2400e+03, 4.2600e+03, 4.2800e+03, 4.3000e+03, 4.3300e+03, 4.3500e+03,
      4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03,
      4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03, 4.3500e+03
    };
  G4double EmissionFibBC484[nEntriesBC484] =
    {
      0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
      0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
      0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00,
      0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 0.0000e+00, 9.9600e-03, 4.6108e-02, 1.4920e-01, 4.3929e-01,
      8.7405e-01, 9.8839e-01, 8.7609e-01, 7.5191e-01, 6.5761e-01, 6.1152e-01, 5.8572e-01, 5.6109e-01, 5.0687e-01, 4.6027e-01,
      4.1532e-01, 3.6851e-01, 3.2695e-01, 2.9074e-01, 2.5784e-01, 2.2847e-01, 1.9965e-01, 1.7665e-01, 1.5485e-01, 1.3599e-01,
      1.1646e-01, 1.0068e-01, 8.6331e-02, 7.1257e-02, 5.9146e-02, 4.7958e-02
    };

  G4double PhotonEnergyBC484[nEntriesBC484];
  G4double RefractiveIndexBC484[nEntriesBC484];
  G4double AbsWLSfiberBC484_sorted[nEntriesBC484];
  G4double EmissionFibBC484_sorted[nEntriesBC484];
  
  for(int ii = 0; ii < nEntriesBC484; ii++) {
    PhotonEnergyBC484[ii]    = 1240.0/(PhotonWaveLength_HC[nEntriesBC484-ii-1]) *eV;
    RefractiveIndexBC484[ii] = 1.58;
    AbsWLSfiberBC484_sorted[ii]     = AbsWLSfiberBC484[nEntriesBC484-ii-1]*cm;
    EmissionFibBC484_sorted[ii] = EmissionFibBC484[nEntriesBC484-ii-1];
  }
  G4MaterialPropertiesTable* MPT_BC484 = new G4MaterialPropertiesTable();
  MPT_BC484->AddProperty("RINDEX"              , PhotonEnergyBC484 , RefractiveIndexBC484 , nEntriesBC484);
  MPT_BC484->AddProperty("WLSABSLENGTH"        , PhotonEnergyBC484 , AbsWLSfiberBC484_sorted     , nEntriesBC484);
  MPT_BC484->AddProperty("WLSCOMPONENT"        , PhotonEnergyBC484 , EmissionFibBC484_sorted     , nEntriesBC484);
  MPT_BC484->AddConstProperty("WLSTIMECONSTANT", 3.0*ns);
  if( fMaterialsListOpticalPhotonDisabled.find( "BC484" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    BC484->SetMaterialPropertiesTable(MPT_BC484);
  }
  fMaterialsMap["BC484"] = BC484;

  //GSO
  G4Material* GSiO = new G4Material("GSiO", density=6.71*g/cm3, nel=3);
  GSiO->AddElement(Gd, natoms=2);
  GSiO->AddElement(Si, natoms=1);
  GSiO->AddElement(O , natoms=5);
  GSiO->GetIonisation()->SetBirksConstant(5.25); 
  fMaterialsMap["GSiO"] = GSiO;
  
  //Glass
  G4Material* Glass_HC = new G4Material("Glass_HC", density=1.032*g/cm3,2);
  Glass_HC->AddElement(C , 91.533*perCent);
  Glass_HC->AddElement(H ,  8.467*perCent);


  //no need to sort these since they're all the same
  G4double Glass_HC_RIND[nEntriesEJ232] =
    { 
      1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 
      1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 
      1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 
      1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 
      1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 
      1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 
      1.58 , 1.58 , 1.58 , 1.58 , 1.58 , 1.58 
    };

  G4double Glass_HC_AbsLength[nEntriesEJ232]=
    { 
      420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm ,
      420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm ,
      420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm ,
      420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm ,
      420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm ,
      420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm ,
      420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm , 420.0*cm 
    };

  G4MaterialPropertiesTable *Glass_HC_mt = new G4MaterialPropertiesTable();
  Glass_HC_mt->AddProperty("ABSLENGTH", PhotonEnergyBC484 , Glass_HC_AbsLength , nEntriesEJ232 );
  Glass_HC_mt->AddProperty("RINDEX"   , PhotonEnergyBC484 , Glass_HC_RIND      , nEntriesEJ232 );

  if( fMaterialsListOpticalPhotonDisabled.find( "Glass_HC" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    Glass_HC->SetMaterialPropertiesTable(Glass_HC_mt);
  }
  fMaterialsMap["Glass_HC"] = Glass_HC; 

  // -- Bialkali photocathode K2CsSb 
  G4Material* BialkaliK2CsSb = new G4Material ("BialkaliK2CsSb" , density= 2.9 *g/cm3 , 3 ); 
  BialkaliK2CsSb->AddElement( K  , 2 );
  BialkaliK2CsSb->AddElement( Cs , 1 );
  BialkaliK2CsSb->AddElement( Sb , 1 );
  fMaterialsMap["BialkaliK2CsSb"] = BialkaliK2CsSb;    

  G4Material* Paper = man->FindOrBuildMaterial("G4_TEFLON");
  
  G4double MilliPoreRefl[nEntriesEJ232] = {
    0.818, 0.828, 0.839, 0.850, 0.860, 0.868, 0.876, 0.883, 0.888, 0.894,
    0.898, 0.902, 0.906, 0.910, 0.913, 0.917, 0.920, 0.923, 0.926, 0.928,
    0.931, 0.933, 0.936, 0.938, 0.940, 0.941, 0.943, 0.945, 0.946, 0.947,
    0.948, 0.949, 0.949, 0.950, 0.950, 0.951, 0.951, 0.951, 0.951, 0.952,
    0.952, 0.952, 0.952, 0.951, 0.951, 0.951, 0.951, 0.951, 0.951, 0.951,
    0.951, 0.950, 0.950, 0.950, 0.950, 0.950, 0.949, 0.949, 0.949, 0.949,
    0.949, 0.948, 0.948, 0.948, 0.948, 0.947
  };
  G4double MilliPoreRefl_sorted[nEntriesEJ232];
  for( int ii = 0; ii<nEntriesEJ232; ii++ ){
    MilliPoreRefl_sorted[ii] = MilliPoreRefl[nEntriesEJ232-ii-1];
  }
  G4double MilliPoreRefrIndexl[nEntriesEJ232] = { 1.5 };
  G4double MilliPoreSS[nEntriesEJ232] = { 0.1 };
  G4double MilliPoreSL[nEntriesEJ232] = { 0.1 };
  G4double MilliPoreBK[nEntriesEJ232] = { 0.1 };

  G4MaterialPropertiesTable *Paper_MPT = new G4MaterialPropertiesTable();
  Paper_MPT->AddProperty("RINDEX"                , PhotonEnergyBC484 , MilliPoreRefrIndexl , nEntriesEJ232);
  Paper_MPT->AddProperty("REFLECTIVITY"          , PhotonEnergyBC484 , MilliPoreRefl_sorted       , nEntriesEJ232);
  Paper_MPT->AddProperty("SPECULARLOBECONSTANT"  , PhotonEnergyBC484 , MilliPoreSL , nEntriesEJ232);
  Paper_MPT->AddProperty("SPECULARSPIKECONSTANT" , PhotonEnergyBC484 , MilliPoreSS , nEntriesEJ232);
  Paper_MPT->AddProperty("BACKSCATTERCONSTANT"   , PhotonEnergyBC484 , MilliPoreBK , nEntriesEJ232);

  if( fMaterialsListOpticalPhotonDisabled.find( "Paper" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    Paper->SetMaterialPropertiesTable(Paper_MPT);
  }
  fMaterialsMap["Paper"] = Paper;

  G4double Foil_refl[nEntriesEJ232] = 
    {
      0.931, 0.932, 0.932, 0.932, 0.931, 0.930, 0.929, 0.928, 0.926, 0.924,
      0.922, 0.920, 0.918, 0.916, 0.913, 0.911, 0.909, 0.906, 0.904, 0.901,
      0.899, 0.896, 0.894, 0.892, 0.890, 0.887, 0.885, 0.883, 0.881, 0.880,
      0.878, 0.876, 0.874, 0.873, 0.872, 0.870, 0.869, 0.868, 0.867, 0.866,
      0.865, 0.864, 0.863, 0.863, 0.862, 0.862, 0.861, 0.861, 0.861, 0.860,
      0.860, 0.860, 0.860, 0.860, 0.860, 0.860, 0.860, 0.860, 0.860, 0.860,
      0.860, 0.860, 0.860, 0.860, 0.860, 0.861
    };

  G4double Foil_refl_sorted[nEntriesEJ232];
  for( int ii=0; ii<nEntriesEJ232; ii++ ){
    Foil_refl_sorted[ii] = Foil_refl[nEntriesEJ232-ii-1];
  }
  
  G4double Zero_refl[nEntriesEJ232];

  G4double Foil_efficiency[nEntriesEJ232];
  for(int ii = 0; ii < nEntriesBC484; ii++) 
    {
      Foil_efficiency[ii] = 0.2;
      Zero_refl[ii]       = 0.0;
    }
  G4OpticalSurface * OpFoilSurface = new G4OpticalSurface("FoilSurface", unified , polished , dielectric_metal);

  G4MaterialPropertiesTable *foil_mpt = new G4MaterialPropertiesTable();
  foil_mpt->AddProperty("REFLECTIVITY", PhotonEnergyBC484 , Foil_refl_sorted , nEntriesEJ232 );
  OpFoilSurface->SetMaterialPropertiesTable(foil_mpt);

  fOpticalSurfacesMap["Foil"] = OpFoilSurface;

  // -- Surface between far end of WLS and light absorber 
  G4double Absorber_reflectivity[nEntriesEJ232];
  G4double Absorber_efficiency[nEntriesEJ232];
  for(int ii = 0; ii < nEntriesBC484; ii++) {
    Absorber_efficiency[ii]   = 1.0;
    Absorber_reflectivity[ii] = 0.0;
  }

  G4OpticalSurface * AbsorberSurface = new G4OpticalSurface("AbsorberSurface", unified , polished , dielectric_metal);
  
  G4MaterialPropertiesTable *Absorber_mpt = new G4MaterialPropertiesTable();
  
  Absorber_mpt->AddProperty("REFLECTIVITY", PhotonEnergyBC484 , Absorber_reflectivity , nEntriesEJ232 );
  AbsorberSurface->SetMaterialPropertiesTable(Absorber_mpt);
  fOpticalSurfacesMap["Absorber"] = AbsorberSurface;

  // -- Surface between WLS and air
  G4double WLS_AirGapRefrIndex[nEntriesEJ232];
  G4double WLS_reflectivity[nEntriesEJ232];
  G4double WLS_SpecularLobe[nEntriesEJ232];
  G4double WLS_SpecularSpike[nEntriesEJ232];
  G4double WLS_Backscatter[nEntriesEJ232];
  
  for(int ii = 0; ii < nEntriesBC484; ii++) {
    WLS_AirGapRefrIndex[ii] = 1.0;
    WLS_reflectivity[ii]    = 0.0;
    WLS_SpecularLobe[ii]    = 0.1;
    WLS_SpecularSpike[ii]   = 0.9;
    WLS_Backscatter[ii]     = 0.0;
  }

  G4OpticalSurface *osWLSToAir = new G4OpticalSurface("osWLSToAir" , unified , ground , dielectric_dielectric , sigalpha=1.3*deg);
  
  G4MaterialPropertiesTable *osWLSToAir_mpt = new G4MaterialPropertiesTable();
  osWLSToAir_mpt->AddProperty("SPECULARLOBECONSTANT"   , PhotonEnergyBC484 , WLS_SpecularLobe    , nEntriesBC484);
  osWLSToAir_mpt->AddProperty("SPECULARSPIKECONSTANT"  , PhotonEnergyBC484 , WLS_SpecularSpike   , nEntriesBC484);
  osWLSToAir_mpt->AddProperty("BACKSCATTERCONSTANT"    , PhotonEnergyBC484 , WLS_Backscatter     , nEntriesBC484);
  
  osWLSToAir->SetMaterialPropertiesTable(osWLSToAir_mpt);
  fOpticalSurfacesMap["osWLSToAir"] = osWLSToAir;

  //  *****************************
  //  *  Large Angle Calorimeter  *
  //  *****************************
  G4double d_LAC_Scinti = 1.032*g/cm3;
  
  G4Material *LAC_Scinti = new G4Material( "LAC_Scinti", d_LAC_Scinti,  2 );
  LAC_Scinti->AddMaterial(PolyVinylToluene, fractionmass = 0.40);
  LAC_Scinti->AddMaterial(Anthracene, fractionmass = 0.60);
  fMaterialsMap["LAC_Scinti"] = LAC_Scinti;

  G4Material *Teflon = man->FindOrBuildMaterial("G4_TEFLON");
  fMaterialsMap["Teflon"] = Teflon;
  
  //   ************************
  //   *          C16         *
  //   ************************
  G4Material *Pyrex_Glass = man->FindOrBuildMaterial("G4_Pyrex_Glass");

  const G4int nentries_rindex_pyrex = 12;
  G4double Ephoton_rindex_pyrex[nentries_rindex_pyrex] =
    { 1.7551*eV, 1.8894*eV, 1.9261*eV, 1.9595*eV, 2.1042*eV, 2.2706*eV,
      2.5509*eV, 2.5833*eV, 2.8453*eV, 3.0640*eV, 3.3973*eV, 3.7115*eV };

  G4double Rindex_pyrex[nentries_rindex_pyrex] =
    { 1.5129, 1.5143, 1.5147, 1.5151, 1.5167, 1.5187,
      1.5224, 1.5228, 1.5267, 1.5302, 1.5363, 1.5427 };

  const G4int nentries_abslength_pyrex = 20;
  
  G4double Ephoton_abslength_pyrex[nentries_abslength_pyrex] =
    { 1.7714*eV, 1.8788*eV, 2.0000*eV, 2.1379*eV, 2.2711*eV,
      2.4800*eV, 2.6957*eV, 2.8440*eV, 2.9524*eV, 3.0617*eV,
      3.1000*eV, 3.1795*eV, 3.2632*eV, 3.3514*eV, 3.3973*eV,
      3.5429*eV, 3.7126*eV, 3.8750*eV, 4.0000*eV, 4.1333*eV };   

  G4double abslength_pyrex[nentries_abslength_pyrex] =
    { 561.6245*cm, 457.4576*cm, 457.4576*cm, 499.1244*cm, 561.6245*cm,
      457.4576*cm, 344.3622*cm, 322.0407*cm, 344.3622*cm, 344.3622*cm,
      322.0407*cm, 237.7600*cm, 144.0809*cm, 109.0256*cm, 83.8915*cm,
      29.8914*cm, 10.0400*cm, 3.8246*cm, 1.8024*cm, 0.8234*cm };

  MPT_temp = new G4MaterialPropertiesTable();

  MPT_temp->AddProperty("RINDEX", Ephoton_rindex_pyrex, Rindex_pyrex, nentries_rindex_pyrex );
  MPT_temp->AddProperty("ABSLENGTH", Ephoton_abslength_pyrex, abslength_pyrex, nentries_abslength_pyrex );

  if( fMaterialsListOpticalPhotonDisabled.find( "Pyrex_Glass" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    Pyrex_Glass->SetMaterialPropertiesTable( MPT_temp );
  }
  
  fMaterialsMap["Pyrex_Glass"] = Pyrex_Glass; 

  G4Material *SiO2_C16 = new G4Material("SiO2_C16", density = 0.12*g/cm3, nel = 2 );
  SiO2_C16->AddElement(elSi, 1);
  SiO2_C16->AddElement(elO, 2);
  fMaterialsMap["SiO2_C16"] = SiO2_C16;
  
  ///////////////////////////
  // DVCS calorimeter 
  ///////////////////////////
  const G4int nentries_SiPM = 31;
  G4double Ephoton_SiPM_data[nentries_SiPM] = 
    { 1.37760*eV, 1.40891*eV, 1.44168*eV, 1.47600*eV, 1.51200*eV, 
      1.54980*eV, 1.58954*eV, 1.63137*eV, 1.67546*eV, 1.72200*eV, 
      1.77120*eV, 1.82330*eV, 1.87855*eV, 1.93725*eV, 1.99974*eV, 
      2.06640*eV, 2.13766*eV, 2.21400*eV, 2.29600*eV, 2.38431*eV, 
      2.47968*eV, 2.58300*eV, 2.69531*eV, 2.81782*eV, 2.95200*eV, 
      3.09960*eV, 3.26274*eV, 3.44401*eV, 3.64659*eV, 3.75710*eV, 
      3.87451*eV };
  
  G4double SiPM_RefIndex[nentries_SiPM] = 
    { 1.41, 1.41, 1.41, 1.41, 1.41, 
      1.41, 1.41, 1.41, 1.41, 1.41, 
      1.41, 1.41, 1.41, 1.41, 1.41, 
      1.41, 1.41, 1.41, 1.41, 1.41, 
      1.41, 1.41, 1.41, 1.41, 1.41, 
      1.41, 1.41, 1.41, 1.41, 1.41, 
      1.41 };
  
  G4double SiPM_QE[nentries_SiPM] = 
    { 0.0419274, 0.0519399, 0.0619524, 0.0750939, 0.0869837, 
      0.100751, 0.115144, 0.131414, 0.149562, 0.168961, 
      0.188986, 0.210889, 0.231539, 0.262203, 0.299124, 
      0.329787, 0.359825, 0.394243, 0.428035, 0.459324, 
      0.482478, 0.491239, 0.500000, 0.488736, 0.469962, 
      0.441802, 0.379224, 0.300375, 0.169587, 0.0869837, 
      0.0319149 };
  
  G4Material* SiPM_Silicon = new G4Material(name="SiPM_Silicon", z=14., 28.086*g/mole, density=2.33*g/cm3);
  fMaterialsMap["SiPM_Silicon"] = SiPM_Silicon;
  
  MPT_temp = new G4MaterialPropertiesTable();
  
  MPT_temp->AddProperty("RINDEX", Ephoton_SiPM_data, SiPM_RefIndex, nentries_SiPM );
  MPT_temp->AddProperty("EFFICIENCY", Ephoton_SiPM_data, SiPM_QE, nentries_SiPM );

  if( fMaterialsListOpticalPhotonDisabled.find( "SiPM_Silicon" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    SiPM_Silicon->SetMaterialPropertiesTable( MPT_temp );
  }
  
  G4Material* PbF2 = new G4Material("PbF2", density=7.77*g/cm3, 2);
  PbF2->AddElement(elPb,1);
  PbF2->AddElement(elF, 2);
  fMaterialsMap["PbF2"] = PbF2;
  
  G4double PbF2refIndex[nentries_SiPM] = 
    { 1.74455, 1.74528, 1.74607, 1.74691, 1.74781, 
      1.74879, 1.74984, 1.75098, 1.75221, 1.75356, 
      1.75502, 1.75663, 1.75840, 1.76034, 1.76250, 
      1.76489, 1.76757, 1.77057, 1.77396, 1.77780, 
      1.78220, 1.78726, 1.79314, 1.80003, 1.80821, 
      1.81804, 1.83007, 1.84514, 1.86466, 1.87690, 
      1.89161 };
 
  MPT_temp = new G4MaterialPropertiesTable();

  MPT_temp->AddProperty("RINDEX", Ephoton_SiPM_data, PbF2refIndex, nentries_SiPM );
  
  if( fMaterialsListOpticalPhotonDisabled.find( "PbF2" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    PbF2->SetMaterialPropertiesTable( MPT_temp );
  }
    
  G4Material* PbWO4 = new G4Material("PbWO4", density=8.28*g/cm3, 3);
  PbWO4->AddElement(elPb,1);
  PbWO4->AddElement(elW, 1);
  PbWO4->AddElement(elO, 4);
  fMaterialsMap["PbWO4"] = PbWO4;
  
  G4double PbWO4refIndex[nentries_SiPM] = 
    { 2.16929, 2.17000, 2.17071, 2.17143, 2.17286, 
      2.17500, 2.17571, 2.17786, 2.18214, 2.18429, 
      2.18500, 2.18929, 2.19071, 2.19643, 2.19857, 
      2.20357, 2.20786, 2.21571, 2.22143, 2.22643, 
      2.23857, 2.24857, 2.26143, 2.27643, 2.29214, 
      2.31643, 2.34929, 2.41071, 2.51000, 2.52786, 
      2.53571 };

  G4double abslength_PbWO4[nentries_SiPM] = 
    { 175.0000*cm, 175.0000*cm, 87.5000*cm, 87.50000*cm, 58.3333*cm, 
       58.3333*cm,  43.7500*cm, 43.7500*cm, 35.00000*cm, 35.0000*cm, 
       29.1667*cm,  25.0000*cm, 25.0000*cm, 21.87500*cm, 21.8750*cm, 
       19.4444*cm,  17.5000*cm, 15.9091*cm, 14.58330*cm, 13.4615*cm, 
       12.5000*cm,  11.6667*cm, 10.2941*cm,  9.21053*cm,  7.29167*cm, 
        5.64516*cm,  3.01724*cm, 1.45833*cm, 0.562701*cm, 0.329567*cm, 
        0.201381*cm };
 
  MPT_temp = new G4MaterialPropertiesTable();

  MPT_temp->AddProperty("RINDEX", Ephoton_SiPM_data, PbWO4refIndex, nentries_SiPM );
  //MPT_temp->AddProperty("ABSLENGTH", Ephoton_SiPM_data, abslength_PbWO4, nentries_SiPM );
  
  //info taken at: http://scintillator.lbl.gov/
  MPT_temp->AddConstProperty("SCINTILLATIONYIELD", 200.0/MeV);
  MPT_temp->AddConstProperty("RESOLUTIONSCALE", 1.0);
  MPT_temp->AddConstProperty("FASTTIMECONSTANT",6.00*ns, true);
  MPT_temp->AddConstProperty("SLOWTIMECONSTANT",6.00*ns, true);
  
  if( fMaterialsListOpticalPhotonDisabled.find( "PbWO4" ) == fMaterialsListOpticalPhotonDisabled.end() ){
    PbWO4->SetMaterialPropertiesTable( MPT_temp );
  }

  ///////////////////////////////////////////
  //      CLAS Large-angle calorimeter:
  ///////////////////////////////////////////

  //Start with scintillator NE-110A (EJ208)

  G4Material *NE110A = new G4Material( "NE110A", density=1.032*g/cm3, 2 );
  NE110A->AddElement( elH, fractionmass = 8.474*perCent );
  NE110A->AddElement( elC, fractionmass = 91.526*perCent );

  fMaterialsMap["NE110A"] = NE110A;
  
  //Anything else?

  //Target collimators for helium-3 target: Tungsten powder epoxy mixture. 
  //For Epoxy let's use Epotek-301-1 properties from the particle data group:
  G4double epoxy_den = 1.190*g/cm3;
  G4double tungsten_den = 19.3*g/cm3;
  G4double collimator_den = 9.53*g/cm3; // from Bert Metzger's drawings (A09016-03-06-0211) // 10.0*g/cm3;

  // f*rho_epoxy + (1-f)*rho_W = rho_coll
  // f*(rho_epoxy - rho_W) = rho_coll - rho_W
  G4double massfrac_epoxy = (1.0 - collimator_den/tungsten_den)/(1.0 - epoxy_den/tungsten_den);
  
  G4Material *TargetCollimator_Material = new G4Material("TargetCollimator_Material", collimator_den, 5 );
  TargetCollimator_Material->AddElement( elH, fractionmass = 0.069894*massfrac_epoxy );
  TargetCollimator_Material->AddElement( elC, fractionmass = 0.689640*massfrac_epoxy );
  TargetCollimator_Material->AddElement( elN, fractionmass = 0.008936*massfrac_epoxy );
  TargetCollimator_Material->AddElement( elO, fractionmass = 0.231531*massfrac_epoxy );
  TargetCollimator_Material->AddElement( elW, fractionmass = 1.0-massfrac_epoxy );

  fMaterialsMap["TargetCollimator_Material"] = TargetCollimator_Material;

  // [for a test object] pure tungsten for a target collimator  
  G4Material *TargetBeamCollimator = new G4Material("TargetBeamCollimator_Material",tungsten_den,1); 
  TargetBeamCollimator->AddElement(elW,1); 
  fMaterialsMap["TargetBeamCollimator_Material"] = TargetBeamCollimator; 
  
}

G4Material *G4SBSDetectorConstruction::GetMaterial(G4String name){

  if( fMaterialsMap.empty() ) ConstructMaterials();

  map<G4String, G4Material*>::iterator it = fMaterialsMap.find( name );

  if( it != fMaterialsMap.end() ){ 
    return fMaterialsMap[name];
  } else {
    fprintf(stderr, "ERROR %s:%d - Material %s not found\n", __FILE__, __LINE__, name.data());
    exit(1);
    return NULL;
  }
}

G4OpticalSurface *G4SBSDetectorConstruction::GetOpticalSurface( G4String name ){
  if( fOpticalSurfacesMap.empty() ) ConstructMaterials();
  map<G4String, G4OpticalSurface*>::iterator it = fOpticalSurfacesMap.find( name );
  if( it != fOpticalSurfacesMap.end() ){
    return fOpticalSurfacesMap[name];
  } else {
    it = fOpticalSurfacesMap.begin();
    return it->second;
  }
}

//G4VPhysicalVolume* G4SBSDetectorConstruction::ConstructAll()
G4VPhysicalVolume* G4SBSDetectorConstruction::ConstructAll()
{
  G4cout << "\nG4SBSDetectorConstruction....\n" << G4endl;

  if( fMaterialsMap.empty() ) ConstructMaterials();

  //    TrackerIDnumber = 0; //Initialize TrackerIDnumber to zero. This gets incremented with each call to TrackerBuilder::BuildComponent!
  //TrackerArm.clear();  //Clear mapping of tracker modules to spectrometers (E arm or H arm)

  fSDman = G4SDManager::GetSDMpointer();
  //--------- Material definition was moved to ConstructMaterials()---------
  //--------- G4VSolid, G4LogicalVolume, G4VPhysicalVolume  ---------

  //--------------
  // World:
  //--------------
  //G4Box *WorldBox= new G4Box("WorldBox",20*m, 20*m, 28*m);
  G4Box *WorldBox= new G4Box("WorldBox",50*m, 50*m, 50*m);
  G4LogicalVolume *WorldLog=new G4LogicalVolume(WorldBox,GetMaterial("Air"),
						"WorldLogical", 0, 0, 0);
  G4PVPlacement *WorldPhys=new G4PVPlacement(0,G4ThreeVector(),
					     "WorldPhysical",
					     WorldLog,
					     0,false,0);

  // In the new version of ConstructAll(), now we call individual modular subsystem creation routines depending on fExpType:

  //All three types of experiments have a target:
  //Target builder is called first:
  fTargetBuilder->BuildComponent(WorldLog); 

  //Beamline builder is called second.
  //All three types of experiments have a beam line:
  fBeamlineBuilder->BuildComponent(WorldLog);

  fEArmBuilder->BuildComponent(WorldLog);
  fHArmBuilder->BuildComponent(WorldLog);
  // I'd suggest not to build RTPC (nor mTPC) in here but in TargetBuilder anyway (E. Fuchey 2018/02/12)
  // DetectorConstruction * RTPCWorld;
  //fRTPC->Construct(RTPCWorld);

  G4FieldManager *fm = new G4FieldManager(fGlobalField);

  // G4Mag_UsualEqRhs* fequation = new G4Mag_UsualEqRhs(fGlobalField);
  // G4MagIntegratorStepper *stepper = new G4ExplicitEuler(fequation, 8);
  // new G4ChordFinder(fGlobalField, 1.0e-2*mm, stepper);

  G4Mag_SpinEqRhs* fBMTequation = new G4Mag_SpinEqRhs(fGlobalField);
  //G4MagIntegratorStepper *pStepper = new G4ClassicalRK4(fBMTequation,12);
  G4MagIntegratorStepper *pStepper = new G4DormandPrince745(fBMTequation,12);
  G4ChordFinder *cftemp = new G4ChordFinder(fGlobalField, 1.0e-2*mm, pStepper);

  fm->SetChordFinder(cftemp);
  
  if( fUseGlobalField ){
    WorldLog->SetFieldManager(fm,true);
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
  WorldLog->SetVisAttributes(G4VisAttributes::GetInvisible());




  G4VisAttributes * dVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  dVisAtt->SetForceWireframe(true);
  //dLog->SetVisAttributes(dVisAtt);

  G4VisAttributes * xVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  xVisAtt->SetForceWireframe(true);
  //xLog->SetVisAttributes(xVisAtt);

  /*
    hcallog->SetVisAttributes(G4VisAttributes::GetInvisible());
    big48d48Log->SetVisAttributes(G4VisAttributes::GetInvisible());
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

void G4SBSDetectorConstruction::SetBigBiteField(int n, G4String fname){
  G4RotationMatrix rm;

  switch(n){
  case 1:
    rm.rotateY(-fEArmBuilder->fBBang); //rotation is negative about y for BB on beam left

    fbbfield = new G4SBSBigBiteField( G4ThreeVector(0.0, 0.0, fEArmBuilder->fBBdist),  rm, fname );

    fbbfield->fScaleFactor = fFieldScale_BB;
    fbbfield->fArm = G4SBS::kEarm;
    // Dimensions of the box
    fGlobalField->AddField(fbbfield);
    if( !fUseGlobalField ) fEArmBuilder->fUseLocalField = true;
    break;
  case 0:
    fGlobalField->DropField(fbbfield);
    //if( fUseGlobalField ) fEArmBuilder->fUseLocalField = false;
    fEArmBuilder->fUseLocalField = false; //this should be set regardless of whether the global field flag has been set, so that it disables creation of the local field in EArmBuilder!!!
    delete fbbfield;
    fbbfield = NULL;
    break;
  default:
    break;
  }
  return;
}
void G4SBSDetectorConstruction::Set48D48Field(int n){
  G4RotationMatrix rm;

  switch(n){
  case 1:
    if( f48d48field ){
      fGlobalField->DropField(f48d48field);
      delete f48d48field;
    }
    //rm.rotateY(-fHArmBuilder->f48D48ang);
    rm.rotateY(fHArmBuilder->f48D48ang); //rotation about y axis is positive for SBS on beam right
	    
    f48d48field = new G4SBSConstantField( 
					 G4ThreeVector(0.0, 0.0, fHArmBuilder->f48D48dist),  rm,
					 // Dimensions of the box
					 G4ThreeVector(469.9*mm/2+0.1*mm, 187.*cm/2.-263.7*mm,  1219.2*mm/2+0.1*mm), 
					 G4ThreeVector(f48D48_uniform_bfield, 0.0, 0.0)
					  );
    f48d48field->fScaleFactor = fFieldScale_SBS;
    f48d48field->fArm = G4SBS::kHarm;
    
    fGlobalField->AddField(f48d48field);
    if( !fUseGlobalField ) fHArmBuilder->fUseLocalField = true;
    break;
  case 0:
    fGlobalField->DropField(f48d48field);
    //if( fUseGlobalField ) fHArmBuilder->fUseLocalField = false;
    fHArmBuilder->fUseLocalField = false; //this should be set regardless of whether the global field flag has been set, so that it disables creation of the local field in HArmBuilder!!!
    delete f48d48field;
    f48d48field = NULL;
    break;
  default: //do nothing
    break;
  }
  return;
}

void G4SBSDetectorConstruction::SetBBDist(double a){ 
  fEArmBuilder->SetBBDist(a); 
  if( fbbfield ) fbbfield->SetOffset(G4ThreeVector(0.0, 0.0, a) ); 
}

void G4SBSDetectorConstruction::SetBBAng(double a){ 
  fEArmBuilder->SetBBAng(a); 
  G4RotationMatrix rm;
  rm.rotateY(-a); //for BB on beam left, rotation about Y is negative (CCW).
  if( fbbfield ) fbbfield->SetRM(rm); 
}

void G4SBSDetectorConstruction::Set48D48Dist(double a){ 
  fHArmBuilder->Set48D48Dist(a); 
  if( f48d48field )  f48d48field->SetOffset(G4ThreeVector(0.0, 0.0, a+ 48.0*2.54*cm/2.0 ) ); 
}

void G4SBSDetectorConstruction::Set48D48Ang(double a){ 
  fHArmBuilder->Set48D48Ang(a); 
  G4RotationMatrix rm;
  rm.rotateY(a); //for SBS on beam right, rotation about Y is positive (CW).
  if( f48d48field ) f48d48field->SetRM(rm); 
}

void G4SBSDetectorConstruction::SetUniformMagneticField48D48( double B ) { 
  f48D48_uniform_bfield = B;
  //Force scale factor to 1, in case something different is floating around?
  if( f48d48field ){
    G4SBSConstantField *f = dynamic_cast<G4SBSConstantField *>(f48d48field);
    if( f ){
      // Relative - is so protons bend up
      f->SetFieldVector(G4ThreeVector(f48D48_uniform_bfield, 0.0, 0.0));
    }
    //fHArmBuilder->fFieldStrength = f48D48_uniform_bfield;
  }
  fHArmBuilder->fFieldStrength = f48D48_uniform_bfield;
}

void G4SBSDetectorConstruction::AddToscaField( const char *fn, int flag ) { 

  if( f48d48field && flag != 1 ){ //this is either a global field definition or a "local" SBS field definition: delete local SBS uniform field if it exists
    fGlobalField->DropField( f48d48field );
    delete f48d48field;
    f48d48field = NULL;
  }

  //Check existence of local BB field. If it exists and flag == 1, then we definitely want to delete it.
  //Otherwise the situation is ambiguous, so we leave it alone. This isn't COMPLETELY idiot-proof, but
  //probably the best we can do without additional layers of complexity:
  if( fbbfield && flag == 1 ){
    fGlobalField->DropField( fbbfield );
    delete fbbfield;
    fbbfield = NULL;
  }
   
  fGlobalField->AddToscaField(fn, flag);

  //Depending on flag, optionally override the map header info for map coordinate transformation:
  if( flag == 1 ){
    G4cout << "Warning: overriding map header coordinate transformation info using BigBite angle and distance for TOSCA map " << fn << G4endl;

    fGlobalField->SetAngleAndDistance( fEArmBuilder->fBBang, fEArmBuilder->fBBdist, G4SBS::kEarm );

    //This will make sure that if the angle and distance get changed AFTER this command is invoked, the field map rotation matrix
    //and offset will also be changed:
    fGlobalField->SetOverride_Earm( true );
  }

  if( flag == 2 ){
    G4cout << "Warning: overriding map header coordinate transformation using SBS/48D48 angle and distance for TOSCA map " << fn << G4endl
	   << "This will have unintended consequences if SBS/48D48 angle and distance are not already properly set   " << G4endl;
    fGlobalField->SetAngleAndDistance( fHArmBuilder->f48D48ang, fHArmBuilder->f48D48dist, G4SBS::kHarm );

    //This will make sure that if the angle and distance get changed AFTER this command is invoked, the field map rotation matrix
    //and offset will also be changed:
    fGlobalField->SetOverride_Harm( true );
  }
    
  G4cout << "fGlobalField::AddToscaField done" << G4endl;
  
  //When creating for the first time, initialize overall scale factor based on fFieldScale_SBS (defaults to 1, is overridden by messenger command).
  //f48d48field->fScaleFactor = fFieldScale_SBS;
  
  if( !fUseGlobalField ){
    fUseGlobalField = true;
    fEArmBuilder->fUseLocalField = false;
    fHArmBuilder->fUseLocalField = false;
  }

  G4cout << "fUseGlobalField set" << G4endl;

}

void G4SBSDetectorConstruction::SetECALmapfilename( G4String sname ){
  fECALmapfilename = sname;
}

void G4SBSDetectorConstruction::SetCDetconfig( int cdetconfig ){
  fCDetOption = cdetconfig;
}

void G4SBSDetectorConstruction::SetC16Segmentation( int segmentC16 ){
  fSegmentC16 = segmentC16;
}

void G4SBSDetectorConstruction::SetSegmentThickC16( G4double thick ){
  fSegmentThickC16 = fabs(thick);
}

void G4SBSDetectorConstruction::SetDoseRateC16( G4double rate ){
  fDoseRateC16 = rate; 
}

void G4SBSDetectorConstruction::SetFieldScale_SBS( G4double v ){
  fFieldScale_SBS = v;
  //(fHArmBuilder->fFieldStrength) *= v;
  if( f48d48field ){
    f48d48field->fScaleFactor = v;
  }

  fGlobalField->ScaleFields( v, G4SBS::kHarm );
  
}

void G4SBSDetectorConstruction::SetFieldScale_BB( G4double v ){
  fFieldScale_BB = v;

  if( fbbfield ){
    fbbfield->fScaleFactor = v;
  }

  fGlobalField->ScaleFields( v, G4SBS::kEarm );
}

void G4SBSDetectorConstruction::SetFlipGEM( G4bool b ){
  fGEMflip = b;
}

void G4SBSDetectorConstruction::SetThresholdTimeWindowAndNTimeBins( G4String SDname, G4double Edefault, G4double Tdefault, G4int Ntbinsdefault ){
  G4SBSGEMSD *GEMSDptr;
  G4SBSCalSD *CalSDptr;
  G4SBSECalSD *ECalSDptr;

  G4double timewindow, threshold;
  G4int ntimebins;
  
  if( SDlist.find( SDname ) != SDlist.end() ){
    switch( SDtype[SDname] ){
    case G4SBS::kGEM: //For now, do nothing:
     
      break;
    case G4SBS::kCAL:
      CalSDptr = (G4SBSCalSD*) fSDman->FindSensitiveDetector( SDname, false );

      timewindow = SDgatewidth.find(SDname) != SDgatewidth.end() ? SDgatewidth[SDname] : Tdefault;
      threshold  = SDthreshold.find(SDname) != SDthreshold.end() ? SDthreshold[SDname] : Edefault;
      ntimebins = SDntimebins.find(SDname) != SDntimebins.end() ? SDntimebins[SDname] : Ntbinsdefault;

      CalSDptr->SetEnergyThreshold( threshold );
      CalSDptr->SetHitTimeWindow( timewindow );
      CalSDptr->SetNTimeBins( ntimebins );
     
      break;
      // *****
    case G4SBS::kECAL:
      ECalSDptr = (G4SBSECalSD*) fSDman->FindSensitiveDetector( SDname, false );

      timewindow = SDgatewidth.find(SDname) != SDgatewidth.end() ? SDgatewidth[SDname] : Tdefault;      
      threshold  = SDthreshold.find(SDname) != SDthreshold.end() ? SDthreshold[SDname] : Edefault;
      ntimebins = SDntimebins.find(SDname) != SDntimebins.end() ? SDntimebins[SDname] : Ntbinsdefault;

      ECalSDptr->SetPEThreshold( threshold );
      ECalSDptr->SetHitTimeWindow( timewindow );
      ECalSDptr->SetNTimeBins( ntimebins );
      
      break;
   default: //do nothing:

     break;
    }

  }

  // This part is the criminal
  // G4double wbin = timewindow/double(ntimebins);
  // for( int ibin=0; ibin<ntimebins; ibin++ ){
  //   // CalSDptr->hold_tbins.push_back( ibin*wbin + 0.5*wbin );
  //   CalSDptr->hold_tbins.push_back( 0.0 );
  // }
  
  return;
}

void G4SBSDetectorConstruction::SetTPCSolenoidField(){
  if( !fUseGlobalField ) fTargetBuilder->fUseLocalTPCSolenoid = true;
}

void G4SBSDetectorConstruction::InsertSDboundaryVolume( G4String bvname, G4String sdname ){

  if( SDboundaryVolumes.find( bvname ) == SDboundaryVolumes.end() ){ //first occurrence of this BV: clear sdlist before insertion:
    SDboundaryVolumes[bvname].clear();
  }
  SDboundaryVolumes[bvname].insert( sdname );
}
