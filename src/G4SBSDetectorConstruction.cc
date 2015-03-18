#include "G4SBSDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
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

#include <vector>
#include <map>
//#include <pair>

using namespace std;

G4SBSDetectorConstruction::G4SBSDetectorConstruction()
{
    f48D48_uniform_bfield = 1.4*tesla;

    fbbfield = NULL;
    f48d48field = NULL;
    //fIO = NULL;

    fTotalAbs = false;

    fExpType = kNeutronExp;

    ConstructMaterials(); //Now we want to construct all materials at the beginning, so that the physics tables can get built properly!!!

    fTargetBuilder   = new G4SBSTargetBuilder(this);
    fBeamlineBuilder = new G4SBSBeamlineBuilder(this);
    fEArmBuilder     = new G4SBSEArmBuilder(this);
    fHArmBuilder     = new G4SBSHArmBuilder(this);

    fHArmBuilder->fFieldStrength = f48D48_uniform_bfield;

    fGlobalField = new G4SBSGlobalField();

    fUseGlobalField = false;

    fLeadOption = 0;

    SDlist.clear();
    SDtype.clear();
    StepLimiterList.clear();
    
    //    TrackerIDnumber = 0;
    //TrackerArm.clear();
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

    G4Material *Vacuum =new G4Material(name="Vacuum", z=1., a=1.0*g/mole, density=1e-9*g/cm3);
    //Vacuum->SetMaterialPropertiesTable(Std_MPT);
    fMaterialsMap["Vacuum"] = Vacuum;

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
    fMaterialsMap["Iron"] = fMaterialsMap["Fer"];

    fMaterialsMap["Iron"] = man->FindOrBuildMaterial("G4_Fe");

    density = 1.29e-03*g/cm3;
    G4Material* Air = new G4Material(name="Air", density, nel=2);
    Air->AddElement(elN, .7);
    Air->AddElement(elO, .3);


    //    Air->SetMaterialPropertiesTable(Std_MPT);

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


    G4double CO2_RefractiveIndex[nEntries];
    for ( int i = 0; i < nEntries; ++i ) {
        CO2_RefractiveIndex[i]=1.00045;
        CO2_RefractiveIndex[i]=1+(CO2_RefractiveIndex[i]-1)*GC_Gas_Pressure/atmosphere;
    }
    G4MaterialPropertiesTable* MPCO2 = new G4MaterialPropertiesTable();
    MPCO2->AddProperty("RINDEX",PhotonEnergy, CO2_RefractiveIndex, nEntries);
    CO2->SetMaterialPropertiesTable(MPCO2);


    fMaterialsMap["CO2"] = CO2;

    // 1.5 Atmosphere C4F8O for cerkenkov
    G4double density_C4F8O = 9.64*mg/cm3; // density at 1ATM
    G4Material* C4F8O = new G4Material("C4F8O", density_C4F8O*1.5, nel=3);
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
    C4F8O->SetMaterialPropertiesTable(MPC4F8O);

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


    // Grinch Lucite

    G4Material* Lucite=man->FindOrBuildMaterial("G4_LUCITE");
    G4double Lucite_RefractiveIndex[nEntries];
    for ( int i = 0; i < nEntries; ++i ) {
        Lucite_RefractiveIndex[i]=1.5;
    }
    G4MaterialPropertiesTable* MPLucite = new G4MaterialPropertiesTable();
    MPLucite->AddProperty("RINDEX",PhotonEnergy, Lucite_RefractiveIndex, nEntries);
    Lucite->SetMaterialPropertiesTable(MPLucite);
    Lucite->SetName("Lucite");
    fMaterialsMap["Lucite"] = Lucite;


    G4Material* Glass=new G4Material("Glass",1.397*g/cm3,3);
    Glass->AddElement(H,8);
    Glass->AddElement(C,10);
    Glass->AddElement(O,4);
    fMaterialsMap["Glass"] = Glass;

    G4Material* matAl=man->FindOrBuildMaterial("G4_Al");
    matAl->SetName("Al");
    fMaterialsMap["Al"] = matAl;

    G4Element *Fe = man->FindOrBuildElement("Fe");
    G4Element *Cr = man->FindOrBuildElement("Cr");
    G4Element *Ni = man->FindOrBuildElement("Ni");
    G4Element *Mn = man->FindOrBuildElement("Mn");
    G4Element *P = man->FindOrBuildElement("P");
    G4Element *S = man->FindOrBuildElement("S");
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
    Quartz->SetMaterialPropertiesTable(mptQuartz);
    fMaterialsMap["Quartz"] = Quartz;



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

    fOpticalSurfacesMap["Mirrsurf"] = Mirrsurf;

    //G4cout << "Material properties for Mirrsurf:" << G4endl;

    
    //********************************************************************
    //************                 ECAL             **********************
    //********************************************************************  
    //Additional elements needed for the materials:
    G4Element* elK = new G4Element( "Potassium", "K", 19, 39.098*g/mole );
    G4Element* elAs = new G4Element( "Arsenic", "As", 33, 74.922*g/mole );
    G4Element* elPb = new G4Element( "Lead", "Pb", 82, 207.2*g/mole );
    G4Element* elBe = new G4Element( "Beryllium", "Be", 4, 9.012*g/mole );
    
    //Materials necessary to build TF1 aka lead-glass
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

    //Simulating annealing: http://hallaweb.jlab.org/12GeV/SuperBigBite/SBS-minutes/2014/Sergey_Abrahamyan_LGAnnealing_2014.pdf
    const G4int nentries_annealing_model=50;

    G4double Ephoton_annealing_model[nentries_annealing_model] = {
      1.513*eV, 1.531*eV, 1.548*eV, 1.569*eV, 1.590*eV,
      1.609*eV, 1.632*eV, 1.653*eV, 1.676*eV, 1.698*eV,
      1.724*eV, 1.746*eV, 1.773*eV, 1.797*eV, 1.825*eV,
      1.849*eV, 1.878*eV, 1.909*eV, 1.935*eV, 1.965*eV,
      2.001*eV, 2.033*eV, 2.069*eV, 2.100*eV, 2.135*eV,
      2.171*eV, 2.213*eV, 2.259*eV, 2.297*eV, 2.337*eV,
      2.386*eV, 2.431*eV, 2.481*eV, 2.530*eV, 2.587*eV,
      2.638*eV, 2.697*eV, 2.752*eV, 2.820*eV, 2.887*eV,
      2.954*eV, 3.024*eV, 3.102*eV, 3.179*eV, 3.261*eV,
      3.351*eV, 3.447*eV, 3.543*eV, 3.655*eV, 3.752*eV};

    G4double Absorption_avg[nentries_annealing_model] = {
      197.05*cm, 197.62*cm, 197.05*cm, 197.01*cm, 197.48*cm, 
      196.86*cm, 197.39*cm, 196.74*cm, 197.25*cm, 196.61*cm, 
      196.55*cm, 197.08*cm, 196.45*cm, 196.40*cm, 196.94*cm, 
      196.88*cm, 196.85*cm, 195.69*cm, 196.25*cm, 196.21*cm, 
      196.16*cm, 194.94*cm, 192.57*cm, 189.58*cm, 184.89*cm, 
      179.64*cm, 174.91*cm, 167.96*cm, 160.95*cm, 154.04*cm, 
      148.19*cm, 140.62*cm, 137.08*cm, 135.03*cm, 130.38*cm, 
      124.69*cm, 116.54*cm, 106.33*cm,  95.16*cm,  83.84*cm, 
       74.26*cm,  70.30*cm,  60.21*cm,  43.32*cm,  33.15*cm, 
       20.07*cm,   9.22*cm,   5.19*cm,   1.73*cm,   0.58*cm}; 

    //Values come from old GSTAR code written by K.Shestermanov 
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

    MPT_temp = new G4MaterialPropertiesTable();
    MPT_temp->AddProperty("RINDEX", Ephoton_ECAL_QE, Rindex_TF1, nentries_ecal_QE );
    MPT_temp->AddProperty("ABSLENGTH", Ephoton_ECAL_QE, Abslength_TF1, nentries_ecal_QE );

    TF1->SetMaterialPropertiesTable( MPT_temp );
    fMaterialsMap["TF1"] = TF1;

    //****  TF1 implementing annealing model  ****
    G4Material* TF1_anneal = new G4Material("TF1_anneal", 3.86*g/cm3, 4);
    TF1_anneal->AddMaterial(PbO, 0.512);
    TF1_anneal->AddMaterial(SiO2, 0.413);
    TF1_anneal->AddMaterial(K2O, 0.070);
    TF1_anneal->AddMaterial(As2O3, 0.005);

    MPT_temp = new G4MaterialPropertiesTable();
    MPT_temp->AddProperty("RINDEX", Ephoton_ECAL_QE, Rindex_TF1, nentries_ecal_QE );
    MPT_temp->AddProperty("ABSLENGTH", Ephoton_annealing_model, Absorption_avg, nentries_annealing_model );

    TF1_anneal->SetMaterialPropertiesTable( MPT_temp );
    fMaterialsMap["TF1_anneal"] = TF1_anneal;

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

    QuartzWindow_ECal->SetMaterialPropertiesTable( MPT_temp );
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

    Photocathode_material_ecal->SetMaterialPropertiesTable( MPT_temp );
    fMaterialsMap["Photocathode_material_ecal"] = Photocathode_material_ecal;

    //Define optical properties for AIR:
    G4Material *Special_Air = man->FindOrBuildMaterial("G4_AIR");
    MPT_temp = new G4MaterialPropertiesTable();

    G4double Rindex_air[nentries_ecal_QE] = {
      1.003, 1.003, 1.003, 1.003, 1.003,
      1.003, 1.003, 1.003, 1.003, 1.003,
      1.003, 1.003, 1.003, 1.003, 1.003,
      1.003, 1.003, 1.003, 1.003, 1.003,
      1.003, 1.003, 1.003, 1.003, 1.003,
      1.003, 1.003, 1.003, 1.003, 1.003,
      1.003, 1.003, 1.003, 1.003, 1.003,
      1.003, 1.003};  

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

    Special_Air->SetMaterialPropertiesTable( MPT_temp );
    fMaterialsMap["Special_Air"] = Special_Air;

    //CDet & Poly "filter"
    G4Material *Polyethylene = man->FindOrBuildMaterial("G4_POLYETHYLENE");
    fMaterialsMap["Polyethylene"] = Polyethylene;

    G4Material *PLASTIC_SC_VINYLTOLUENE = man->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    fMaterialsMap["PLASTIC_SC_VINYLTOLUENE"] = PLASTIC_SC_VINYLTOLUENE;

    //****************
    //***** HCAL *****
    //**************** 
    //Elements:
    G4Element *Gd = new G4Element("Gadolinium", "Gd" , z=64.0 , a=157.25*g/mole);

    //- EJ-232
    G4Material* EJ232 = new G4Material("EJ232",density=1.02*g/cm3, nel=2,kStateSolid,293.15*kelvin,1.0*atmosphere);
    EJ232->AddElement(H , natoms=11);
    EJ232->AddElement(C , natoms=10);

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

    G4double PhotonEnergyEJ232[nEntriesEJ232];
    G4double RefractiveIndexEJ232[nEntriesEJ232];
    for(int ii = 0; ii < nEntriesEJ232; ii++) {
      PhotonEnergyEJ232[ii] = 1240.0/(PhotonWaveLength_HC[ii]) *eV;
      RefractiveIndexEJ232[ii] = 1.58;
      ABSL_EJ232[ii]           = ABSL_EJ232[ii]*mm;
    }

    G4MaterialPropertiesTable* MPT_EJ232 = new G4MaterialPropertiesTable();
    MPT_EJ232->AddProperty("RINDEX"       , PhotonEnergyEJ232 , RefractiveIndexEJ232 , nEntriesEJ232);
    MPT_EJ232->AddProperty("FASTCOMPONENT", PhotonEnergyEJ232 , FAST_EJ232           , nEntriesEJ232);
    MPT_EJ232->AddProperty("ABSLENGTH",     PhotonEnergyEJ232 , ABSL_EJ232           , nEntriesEJ232);
    MPT_EJ232->AddConstProperty("SCINTILLATIONYIELD", (8400.0/MeV * 5.51591522788642652e-01 * 0.5/1.4)/MeV);
    MPT_EJ232->AddConstProperty("RESOLUTIONSCALE", 1.0);
    MPT_EJ232->AddConstProperty("FASTTIMECONSTANT",1.40*ns);
    MPT_EJ232->AddConstProperty("SLOWTIMECONSTANT",1.40*ns);
    MPT_EJ232->AddConstProperty("YIELDRATIO",1.0);
    EJ232->SetMaterialPropertiesTable(MPT_EJ232);
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

    for(int ii = 0; ii < nEntriesBC484; ii++) {
      PhotonEnergyBC484[ii]    = 1240.0/(PhotonWaveLength_HC[ii]) *eV;
      RefractiveIndexBC484[ii] = 1.58;
      AbsWLSfiberBC484[ii]     = AbsWLSfiberBC484[ii]*cm;        
    }
    G4MaterialPropertiesTable* MPT_BC484 = new G4MaterialPropertiesTable();
    MPT_BC484->AddProperty("RINDEX"              , PhotonEnergyBC484 , RefractiveIndexBC484 , nEntriesBC484);
    MPT_BC484->AddProperty("WLSABSLENGTH"        , PhotonEnergyBC484 , AbsWLSfiberBC484     , nEntriesBC484);
    MPT_BC484->AddProperty("WLSCOMPONENT"        , PhotonEnergyBC484 , EmissionFibBC484     , nEntriesBC484);
    MPT_BC484->AddConstProperty("WLSTIMECONSTANT", 3.0*ns);
    BC484->SetMaterialPropertiesTable(MPT_BC484);
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
    Glass_HC->SetMaterialPropertiesTable(Glass_HC_mt);
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
    G4double MilliPoreRefrIndexl[nEntriesEJ232] = { 1.5 };
    G4double MilliPoreSS[nEntriesEJ232] = { 0.1 };
    G4double MilliPoreSL[nEntriesEJ232] = { 0.1 };
    G4double MilliPoreBK[nEntriesEJ232] = { 0.1 };

    G4MaterialPropertiesTable *Paper_MPT = new G4MaterialPropertiesTable();
    Paper_MPT->AddProperty("RINDEX"                , PhotonEnergyBC484 , MilliPoreRefrIndexl , nEntriesEJ232);
    Paper_MPT->AddProperty("REFLECTIVITY"          , PhotonEnergyBC484 , MilliPoreRefl       , nEntriesEJ232);
    Paper_MPT->AddProperty("SPECULARLOBECONSTANT"  , PhotonEnergyBC484 , MilliPoreSL , nEntriesEJ232);
    Paper_MPT->AddProperty("SPECULARSPIKECONSTANT" , PhotonEnergyBC484 , MilliPoreSS , nEntriesEJ232);
    Paper_MPT->AddProperty("BACKSCATTERCONSTANT"   , PhotonEnergyBC484 , MilliPoreBK , nEntriesEJ232);

    Paper->SetMaterialPropertiesTable(Paper_MPT);
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
    G4double Zero_refl[nEntriesEJ232];

    G4double Foil_efficiency[nEntriesEJ232];
    for(int ii = 0; ii < nEntriesBC484; ii++) 
      {
	Foil_efficiency[ii] = 0.2;
	Zero_refl[ii]       = 0.0;
      }
    G4OpticalSurface * OpFoilSurface = new G4OpticalSurface("FoilSurface", unified , polished , dielectric_metal);

    G4MaterialPropertiesTable *foil_mpt = new G4MaterialPropertiesTable();
    foil_mpt->AddProperty("REFLECTIVITY", PhotonEnergyBC484 , Foil_refl , nEntriesEJ232 );
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
    G4Box *WorldBox= new G4Box("WorldBox",20*m, 20*m, 30*m);
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
	    if( !fUseGlobalField ) fEArmBuilder->fUseLocalField = true;
	    break;
	case 0:
	    fGlobalField->DropField(fbbfield);
	    if( fUseGlobalField ) fEArmBuilder->fUseLocalField = false;
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
	    if( !fUseGlobalField ) fHArmBuilder->fUseLocalField = true;
	    break;
	case 0:
	    fGlobalField->DropField(f48d48field);
	    if( fUseGlobalField ) fHArmBuilder->fUseLocalField = false;
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
	    // Relative - is so protons bend up
	    f->SetFieldVector(G4ThreeVector(f48D48_uniform_bfield, 0.0, 0.0));
	}
	fHArmBuilder->fFieldStrength = f48D48_uniform_bfield;
    }
}

void G4SBSDetectorConstruction::AddToscaField( const char *fn ) { 
    fGlobalField->AddToscaField(fn);

    if( !fUseGlobalField ){
	fUseGlobalField = true;
	fEArmBuilder->fUseLocalField = false;
	fHArmBuilder->fUseLocalField = false;
    }

}

void G4SBSDetectorConstruction::SetECALmapfilename( G4String s ){
  fECALmapfilename = s;
}

void G4SBSDetectorConstruction::SetHCALspecsfilename( G4String s){
  fHCALspecsfilename = s;
}
