// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class ESSNDetectorMessenger
// Online control of detector configuration via keyboard
// 20/05/13 JRMA adapted from SBS equivalent, under construction

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4NistManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "RTPC.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;
class G4VisAttributes;
//class DetectorArray;

class DetectorConstruction : public G4VUserDetectorConstruction
{
private:
  G4int fIsInteractive;   // batch(0) or interactive(1) mode
  G4Material*        fMaterialDefault;           
  G4Box*             fWorldSolid;        // pointer to the solid World 
  G4LogicalVolume*   fWorldLogic;        // pointer to the logical World
  G4VPhysicalVolume* fWorldPhysi;        // pointer to the physical World
  RTPC* fRTPC;                        // source container and shielding
  //DetectorArray* fDA;                    // array of detectors
  G4Material* fTargetMaterial;
  DetectorMessenger* fDetMessenger;  // pointer to the Messenger
  //
  //build detctor flags, should be set by DetectorSetup.mac 
  //or changed interactively eg, /SBS/det/useTAPS 0, then,  /SBS/det/update
  G4String fUseTarget; //Build the target, either "Cryo","Solid"
  G4String fRTPCFile;             // water tank setup configuration
  //G4String fDAFile;             // detector setup configuration
  G4double fTargetRadius;
  G4double fPlugLength;
  G4bool fIsOverlapVol;
  G4bool fIsSrcPb;
  G4bool fIsPlugApp[4];
  G4VisAttributes* fLgreen;     // non-standard colours
  G4VisAttributes* fLblue;
  G4VisAttributes* fGold;
  G4NistManager* fNistManager;
  G4Material* fPlexi;           // materials
  G4Material* fAir;
  G4Material* fWater;
  G4Material* fPolythene;
  G4Material* fGlass;
  G4Material* fPbGlass;
  G4Material* fPlaScint;
  G4Material* fLiqScint;
  G4Material* fYAP;
  G4Material* fPb;
  G4Material* fFe;
  G4Material* fCu;
  G4Material* fBe;
  G4Material* fAl;
  G4Material* fGe;
  G4Material* fW;
  G4Material* fKapton;
  G4Material* fVac;
  G4Material* fCoH2;
  G4Material* fCoD2;
  G4Material* fHe;
  G4Material* fAr;
  G4Material* fStainless;
  G4Material* fConcrete;
public:
  DetectorConstruction();
  ~DetectorConstruction(); 
  G4VPhysicalVolume* Construct();
  void UpdateGeometry();
  void DefineMaterials();
  void SetIsInteractive(G4int is){fIsInteractive=is;}
  //Set functions used by messenger class
  void SetPlugApp(G4int i){fIsPlugApp[i] = 1;}
  void SetTargetMaterial(G4String mat)
  {fTargetMaterial = fNistManager->FindOrBuildMaterial(mat);}
  void SetPlugLength(G4double len){fPlugLength = len;}
  void SetTargetRadius(G4double tr){fTargetRadius=tr*cm;}
  void SetWTFile(G4String file){fRTPCFile=file;}
  void SetIsOverlapVol( G4bool opt ){ fIsOverlapVol = opt; }
  void SetUseSrcPb( G4bool opt ){ fIsSrcPb = opt; }
  //TH3D* GetEloss(){ return fRTPC->GetEloss(); }
  //TH3D* GetElossG(){ return fRTPC->GetElossG(); }
  //TH3D* GetElossN(){ return fRTPC->GetElossN(); }
  G4double GetZext(){ return fRTPC->GetTz(); }
  G4double GetZst(){ return fRTPC->GetZst(); }
  //
  G4LogicalVolume* GetWorldLogic(){ return fWorldLogic; }
  G4LogicalVolume* GetWTLogic(){ return fRTPC->GetLogic(); }
  //G4int GetNArray(){ return fDA->GetNArray(); }
  const G4VPhysicalVolume* GetWorld() {return fWorldPhysi;}
  G4ThreeVector GetSrcPos(){ return fRTPC->GetSrcPos(); }
  G4ThreeVector GetSrcSize(){ return fRTPC->GetSrcSize(); }
  G4VisAttributes* GetLgreen(){ return fLgreen;}
  G4VisAttributes* GetLblue(){ return fLblue;}
  G4VisAttributes* GetGold(){ return fGold;}
  G4int* GetNhits(){ return fRTPC->GetNhits(); }
  G4int* GetHitID(){ return fRTPC->GetHitID(); }
  G4int* GetHits(){ return fRTPC->GetHits(); }

  // List of materials getters
  G4Material* GetPlexi(){
    if(!fPlexi)
      fPlexi = fNistManager->FindOrBuildMaterial("G4_PLEXIGLASS");
    return fPlexi;
  }
  G4Material* GetAir(){
    if(!fAir)
      fAir = fNistManager->FindOrBuildMaterial("G4_AIR");
    return fAir;
  }
  G4Material* GetWater(){
    if(!fWater)
      fWater = fNistManager->FindOrBuildMaterial("G4_WATER");
    return fWater;
  }
  G4Material* GetPolythene(){
    if(!fPolythene)
      fPolythene = fNistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
    return fPolythene;
  }
  G4Material* GetPlaScint(){
    if(!fPlaScint)
      fPlaScint=fNistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
    return fPlaScint;
  }
  G4Material* GetLiqScint(){
    if(!fLiqScint)
      fLiqScint=fNistManager->FindOrBuildMaterial("G4_XYLENE");
    return fLiqScint;
  }
  G4Material* GetGlass(){
    if(!fGlass)
      fGlass = fNistManager->FindOrBuildMaterial("G4_GLASS_PLATE");
    return fGlass;
  }
  G4Material* GetPbGlass(){
    if(!fGlass)
      fPbGlass = fNistManager->FindOrBuildMaterial("G4_GLASS_LEAD");
    return fPbGlass;
  }
  G4Material* GetYAP(){
    if( !fYAP ){
      fYAP = new G4Material("YAP",5.55*g/cm3, 3);
      fYAP->AddElement(fNistManager->FindOrBuildElement(39),1);
      fYAP->AddElement(fNistManager->FindOrBuildElement(13),1);
      fYAP->AddElement(fNistManager->FindOrBuildElement(8),3);
    }
    return fYAP;
  }
  G4Material* GetKapton(){
    if( !fKapton ){
      fKapton = new G4Material("Kapton",1.43*g/cm3, 4);
      fKapton->AddElement(fNistManager->FindOrBuildElement(6),22);
      fKapton->AddElement(fNistManager->FindOrBuildElement(7),2);
      fKapton->AddElement(fNistManager->FindOrBuildElement(8),5);
      fKapton->AddElement(fNistManager->FindOrBuildElement(1),10);
    }
    return fKapton;
  }
  G4Material* GetPb(){
    if(!fPb)
      fPb = fNistManager->FindOrBuildMaterial("G4_Pb");
    return fPb;
  }
  G4Material* GetFe(){
    if(!fFe)
      fFe = fNistManager->FindOrBuildMaterial("G4_Fe");
    return fFe;
  }
  G4Material* GetStainless(){
    if(!fStainless)
      fStainless = fNistManager->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    return fStainless;
  }
  G4Material* GetConcrete(){
    if(!fConcrete)
      fConcrete = fNistManager->FindOrBuildMaterial("G4_CONCRETE");
    return fConcrete;
  }
  G4Material* GetGe(){
    if(!fGe)
      fGe = fNistManager->FindOrBuildMaterial("G4_Ge");
    return fGe;
  }
  G4Material* GetBe(){
    if(!fBe)
      fBe = fNistManager->FindOrBuildMaterial("G4_Be");
    return fBe;
  }
  G4Material* GetAl(){
    if(!fAl)
      fAl = fNistManager->FindOrBuildMaterial("G4_Al");
    return fAl;
  }
  G4Material* GetW(){
    if(!fW)
      fW = fNistManager->FindOrBuildMaterial("G4_W");
    return fW;
  }
  G4Material* GetCu(){
    if(!fCu)
      fCu = fNistManager->FindOrBuildMaterial("G4_Cu");
    return fCu;
  }
  G4Material* GetAr(){
    if(!fAr)
      fAr = fNistManager->FindOrBuildMaterial("G4_Ar");
    return fAr;
  }
  G4Material* GetVac(){
    if(!fVac)
      fVac = fNistManager->FindOrBuildMaterial("G4_Galactic");
    return fVac;
  }
  // cold hydrogen gas
  G4Material* GetCoH2(G4double dens){
    if( !fCoH2 ){
      fCoH2 = new G4Material("CoH2",dens*g/cm3, 1);
      fCoH2->AddElement(fNistManager->FindOrBuildElement(1),1);
      //G4cout << fCoH2 << endl;
    }
    return fCoH2;
  }
  // cold deuterium gas
  G4Material* GetCoD2(G4double dens){
    if( !fCoD2 ){
      G4Isotope* d2 = new G4Isotope("2H",1,2,2.014*g/mole);
      G4Element* deut = new G4Element("Deut","2D",1);
      deut->AddIsotope(d2,100*perCent);
      fCoD2 = new G4Material("CoD2",dens*g/cm3, 1);
      fCoD2->AddElement(deut,1);
      //G4cout << fCoD2 << endl;
    }
    return fCoD2;
  }
  // He gas, specify density
  G4Material* GetHe(G4double dens){
    if(!fHe){
      fHe = new G4Material("GHe",dens*g/cm3, 1);
      fHe->AddElement(fNistManager->FindOrBuildElement(2),1);
      //G4cout << fHe << endl;
    }
    return fHe;
  }
};

#endif

