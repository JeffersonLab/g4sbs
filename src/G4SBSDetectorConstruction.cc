#include "G4SBSDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4UserLimits.hh"

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

#include <vector>
#include <map>
//#include <pair>

using namespace std;

G4SBSDetectorConstruction::G4SBSDetectorConstruction()
{
    f48D48_uniform_bfield = 1.4*tesla;

    fbbfield = NULL;
    f48d48field = NULL;


    fTotalAbs = true;

    fExpType = kNeutronExp;

    ConstructMaterials(); //Now we want to construct all materials at the beginning, so that the physics tables can get built properly!!!

    fTargetBuilder   = new G4SBSTargetBuilder(this);
    fBeamlineBuilder = new G4SBSBeamlineBuilder(this);
    fEArmBuilder     = new G4SBSEArmBuilder(this);
    fHArmBuilder     = new G4SBSHArmBuilder(this);

    fGlobalField = new G4SBSGlobalField();
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

    fMaterialsMap["BlandAir"] = BlandAir;

    G4double a, iz, z;

    //Seamus's element definitions:
    a = 126.9*g/mole;
    G4Element* elI = new G4Element(name="Iodine",   symbol="I",  iz=53., a);
    a = 132.9*g/mole;
    G4Element* elCs= new G4Element(name="Cesium",   symbol="Cs", iz=55., a);

    G4Element *elH = new G4Element("Hydrogen", "H", 1, 1.007*g/mole );
    G4Element *elD = new G4Element("Deuterium", "D", 1, 2.014*g/mole );
    G4Element *el3He = new G4Element("Helium3", "3He", 2, 3.016*g/mole );
    G4Element *elC = new G4Element("Carbon", "C", 6, 12.011*g/mole ) ;
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

    fMaterialsMap["Vacuum"] = new G4Material(name="Vacuum", z=1., a=1.0*g/mole, density=1e-9*g/cm3);



    if( fMaterialsMap.find("Lead") == fMaterialsMap.end() ){ 
	fMaterialsMap["Lead"] = new G4Material(name="Lead", z=82., a=208.0*g/mole, density=11.34*g/cm3);
    }
    if( fMaterialsMap.find("Aluminum") == fMaterialsMap.end() ){ 
	fMaterialsMap["Aluminum"] = new G4Material(name="Aluminum", z=13., a=26.98*g/mole, density=2.7*g/cm3);
    }

    density = 4.51*g/cm3;
    G4Material* CsI = new G4Material(name="CsI", density, nel = 2);
    CsI->AddElement(elI, .5);
    CsI->AddElement(elCs,.5);

    fMaterialsMap["CsI"] = CsI;


    a = 55.85*g/mole;
    density = 7.87*g/cm3;
    fMaterialsMap["Fer"] = new G4Material(name="Fer", z=26., a, density);

    density = 1.29e-03*g/cm3;
    G4Material* Air = new G4Material(name="Air", density, nel=2);
    Air->AddElement(elN, .7);
    Air->AddElement(elO, .3);

    fMaterialsMap["Air"] = Air;

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
    fMaterialsMap["CO2"] = CO2;

    // 1.5 Atmosphere C4F8O for cerkenkov
    G4double density_C4F8O = 9.64*mg/cm3; // density at 1ATM
    G4Material* C4F8O = new G4Material("C4F8O", density_C4F8O*1.5, nel=3);
    C4F8O->AddElement(elC, 4);
    C4F8O->AddElement(elF, 8);
    C4F8O->AddElement(elO, 1);
    fMaterialsMap["C4F8O"] = C4F8O;

    G4double density_ArCO2 = .7*density_Ar + .3*density_CO2;
    // Use ArCO2
    G4Material *GEMgas= new G4Material("GEMgas", density_ArCO2, nel=2);
    GEMgas->AddMaterial(Argon, 0.7*density_Ar/density_ArCO2) ;
    GEMgas->AddMaterial(CO2, 0.3*density_CO2/density_ArCO2) ;

    fMaterialsMap["GEMgas"] = GEMgas;

    G4Material *Copper= new G4Material("Copper", z=29, a=   63.55*g/mole, density = 8.96*g/cm3);
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

    fMaterialsMap["CH"] = CH;

    //Target materials:
    double gasden = 10.5*atmosphere*(1.0079*2*g/Avogadro)/(300*kelvin*k_Boltzmann);
    G4Material *refH2 = new G4Material("refH2", gasden, 1 );
    refH2->AddElement(elH, 1);

    fMaterialsMap["refH2"] = refH2;

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
    G4Material* matBe = new G4Material("Berylium", 4., 9.012*g/mole, density=1.85*g/cm3);

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

    //Define AIR for RICH mother volume:
    G4Material *RICH_Air = man->FindOrBuildMaterial("G4_AIR");
    fMaterialsMap["RICH_Air"] = RICH_Air;

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

    UVT_Lucite->SetMaterialPropertiesTable( MPT_temp );
    // G4cout << "Material optical properties for material UVT_Lucite" << G4endl;
    // MPT_temp->DumpTable();

    fMaterialsMap["UVT_Lucite"] = UVT_Lucite;

    G4Material *Steel = man->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    fMaterialsMap["Steel"] = Steel;

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
    G4double Abslength_quartz[nentries_quartz] = {
	15.65792444*cm, 15.78550788*cm, 15.7794917*cm, 15.60910249*cm, 15.72664954*cm, 15.72488912*cm, 15.57290011*cm, 15.68021339*cm, 
	15.73546266*cm, 15.55685196*cm, 15.55490625*cm, 15.63907251*cm, 15.48113765*cm, 15.54074565*cm, 15.39638598*cm, 15.50169846*cm, 
	15.4950396*cm, 15.36125979*cm, 15.41113687*cm, 15.33874196*cm, 15.24165927*cm, 15.25602267*cm, 15.23330157*cm, 15.14071666*cm, 
	15.13642486*cm, 15.06590584*cm, 15.05023293*cm, 14.99006002*cm, 14.91826095*cm, 14.79500397*cm, 14.80590431*cm, 14.66396966*cm, 
	14.57959363*cm, 14.47194788*cm, 14.40952367*cm, 14.20967861*cm, 14.11981056*cm, 13.98888512*cm, 13.79714319*cm, 13.85187177*cm, 
	13.61931079*cm, 13.27721911*cm, 12.75893298*cm, 12.51276543*cm, 12.13753406*cm, 11.84114847*cm, 10.96932156*cm, 10.13778162*cm, 
	9.92333989*cm, 9.282597031*cm, 7.443241349*cm 
    };

    MPT_temp = new G4MaterialPropertiesTable();
    MPT_temp->AddProperty("RINDEX", Ephoton_quartz, Rindex_quartz, nentries_quartz);
    MPT_temp->AddProperty("ABSLENGTH", Ephoton_quartz, Abslength_quartz, nentries_quartz);

    QuartzWindow->SetMaterialPropertiesTable( MPT_temp );
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
    UVglass->SetMaterialPropertiesTable( MPT_temp );

    fMaterialsMap["UVglass"] = UVglass;

    G4double den_C4F10 = 10.62*mg/cm3; //Wow, that really is a heavy gas, heavier 
    				       //than HTCC mirror substrates made of Rohacell31!
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
    fMaterialsMap["C4F10_gas"] = C4F10_gas;

    //Quantum efficiency for PMT photocathode ("typical", from XP1911/UV data sheet):
    const G4int nentries_QE = 40;

    G4double Ephoton_QE[nentries_QE] = { 
	1.939182049*eV, 1.963579994*eV, 1.999929682*eV, 2.0346944*eV, 2.067643295*eV, 2.098539267*eV, 2.127141474*eV, 2.153230717*eV, 
	2.179967872*eV, 2.224858112*eV, 2.264313372*eV, 2.293899315*eV, 2.339752453*eV, 2.359406493*eV, 2.403825804*eV, 2.432973578*eV, 
	2.484622785*eV, 2.538512459*eV, 2.599590307*eV, 2.673824751*eV, 2.768699841*eV, 2.858865976*eV, 2.986530569*eV, 3.105387967*eV, 
	3.168432345*eV, 3.279402652*eV, 3.349799677*eV, 3.526419346*eV, 3.674006182*eV, 3.953218421*eV, 4.139635948*eV, 4.385211047*eV, 
	4.570741196*eV, 4.871970386*eV, 5.177284871*eV, 5.437857883*eV, 5.611675974*eV, 5.845211394*eV, 5.969437475*eV, 6.234402273*eV }; 

    G4double PMT_QuantumEfficiency[nentries_QE] = { 
	0.001154539, 0.001948441, 0.003714689, 0.006298763, 0.009646797, 0.013657231, 0.018010571, 0.022819293, 
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

    //G4double Ephot_Rcathode[2] = {1.77*eV, 6.20*eV};
    //G4double Rcathode[2] = {0.0, 0.0};

    MPT_temp = new G4MaterialPropertiesTable();
    MPT_temp->AddProperty("EFFICIENCY", Ephoton_QE, PMT_QuantumEfficiency, nentries_QE );
    MPT_temp->AddProperty("RINDEX", Ephoton_quartz, Rindex_quartz, nentries_quartz );
    MPT_temp->AddProperty("ABSLENGTH", Ephoton_quartz, Abslength_quartz, nentries_quartz );
    //MPT_temp->AddProperty("REFLECTIVITY", Ephot_Rcathode, Rcathode, 2 );

    Photocathode_material->SetMaterialPropertiesTable( MPT_temp );

    fMaterialsMap["Photocathode_material"] = Photocathode_material;

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
    fMaterialsMap["Aerogel"] = Aerogel;



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

    fOpticalSurfacesMap["MirrorDefault"] = DefaultOpticalSurface;

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

    if( fMaterialsMap.empty() ) ConstructMaterials();

    map<G4String, G4Material*>::iterator it = fMaterialsMap.find( name );

    if( it != fMaterialsMap.end() ){ 
	return fMaterialsMap[name];
    } else {
	it = fMaterialsMap.begin();
	return it->second;
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

G4VPhysicalVolume* G4SBSDetectorConstruction::ConstructAll()
{
    G4cout << "\nG4SBSDetectorConstruction....\n" << G4endl;

    if( fMaterialsMap.empty() ) ConstructMaterials();

    fSDman = G4SDManager::GetSDMpointer();
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
    fTargetBuilder->BuildComponent(WorldLog);

    //All three types of experiments have a beam line:
    fBeamlineBuilder->BuildComponent(WorldLog);

    fEArmBuilder->BuildComponent(WorldLog);
    fHArmBuilder->BuildComponent(WorldLog);

    G4FieldManager *fm = new G4FieldManager(fGlobalField);

    G4Mag_UsualEqRhs* fequation= new G4Mag_UsualEqRhs(fGlobalField);
    G4MagIntegratorStepper *stepper = new G4ExplicitEuler(fequation, 8);
    new G4ChordFinder(fGlobalField, 1.0*nm, stepper);
    WorldLog->SetFieldManager(fm,true);

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

void G4SBSDetectorConstruction::SetBigBiteField(int n){
    G4RotationMatrix rm;

    switch(n){
	case 1:
	    rm.rotateY(fEArmBuilder->fBBang);

	    fbbfield = new G4SBSBigBiteField( 
		    G4ThreeVector(0.0, 0.0, fEArmBuilder->fBBdist),  rm );
		    // Dimensions of the box
	    fGlobalField->AddField(fbbfield);
	    break;
	case 0:
	    fGlobalField->DropField(fbbfield);
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
	    rm.rotateY(fHArmBuilder->f48D48ang);

	    f48d48field = new G4SBSConstantField( 
		    G4ThreeVector(0.0, 0.0, fHArmBuilder->f48D48dist),  rm,
		    // Dimensions of the box
		    G4ThreeVector(469.9*mm/2+0.1*mm, 187.*cm/2.-263.7*mm,  1219.2*mm/2+0.1*mm), 
		    G4ThreeVector(f48D48_uniform_bfield, 0.0, 0.0)
		    );
	    fGlobalField->AddField(f48d48field);
	    break;
	case 0:
	    fGlobalField->DropField(f48d48field);
	    delete f48d48field;
	    f48d48field = NULL;
	    break;
	default:
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
    rm.rotateY(a);
    if( fbbfield ) fbbfield->SetRM(rm); 
}

void G4SBSDetectorConstruction::Set48D48Dist(double a){ 
    fHArmBuilder->Set48D48Dist(a); 
    if( f48d48field )  f48d48field->SetOffset(G4ThreeVector(0.0, 0.0, a+ 48.0*2.54*cm/2 ) ); 
}

void G4SBSDetectorConstruction::Set48D48Ang(double a){ 
    fHArmBuilder->Set48D48Ang(a); 
    G4RotationMatrix rm;
    rm.rotateY(-a);
    if( f48d48field ) f48d48field->SetRM(rm); 
}

void G4SBSDetectorConstruction::SetUniformMagneticField48D48( double B ) { 
    f48D48_uniform_bfield = B; 
    if( f48d48field ){
	G4SBSConstantField *f = dynamic_cast<G4SBSConstantField *>(f48d48field);
	if( f ){
	    f->SetFieldVector(G4ThreeVector(f48D48_uniform_bfield, 0.0, 0.0));
	}
    }
}

