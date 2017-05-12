// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class ESSNDetectorMessenger
// Online control of detector configuration via keyboard
// 20/05/13 JRMA adapted from SBS equivalent, under construction

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4VVisManager.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Point3D.hh"
#include "G4Polyline.hh"
//
//#include "DetectorArray.hh"

DetectorConstruction::DetectorConstruction()
{
 
  // create commands for interactive definition of the calorimeter
  fDetMessenger = new DetectorMessenger(this);
  fIsInteractive=1;
  fRTPC = NULL;
  //fDA = NULL;
  fRTPCFile="data/RTPCparm.dat";
  //fDAFile="data/DAparameters.dat";
  fPlugLength = fTargetRadius = 0;
  fIsOverlapVol = false;
  for(G4int i=0; i<4; i++){ fIsPlugApp[i] = false; }
  fLgreen = new G4VisAttributes( G4Colour(0.0,0.75,0.0) );
  fLblue  = new G4VisAttributes( G4Colour(0.0,0.0,0.75) );
  fGold   = new G4VisAttributes( G4Colour(0.75,0.75,0.0));
  fNistManager=G4NistManager::Instance();
  fPlexi = fAir = fWater = fPolythene = fGlass = fPbGlass = fPlaScint =
    fLiqScint = fYAP = fPb = fFe = fCu = fBe = fAl = fGe = fW = fKapton =
    fVac = fCoH2 = fCoD2 = fHe = fAr = fStainless = fConcrete = NULL;
}

//----------------------------------------------------------------------------
DetectorConstruction::~DetectorConstruction(){
 delete fDetMessenger;
}

//----------------------------------------------------------------------------
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();
  //     
  // World
  //
  fWorldSolid = new G4Box("World", 5*m,5*m,30*m);            
  fWorldLogic = 
    new G4LogicalVolume(fWorldSolid,
			G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR")
			,"World");
                                   
  fWorldPhysi = new G4PVPlacement(0,			//no rotation
				  G4ThreeVector(),	//at (0,0,0)
				  fWorldLogic,		//its logical volume
				  "World",		//its name
				  0,			//its mother  volume
				  false,	       	//no boolean operation
				  0,			//copy number
				  fIsOverlapVol);
  fRTPC=new RTPC(this);
  fRTPC->ReadParameters(fRTPCFile);
  fRTPC->SetIsOverlapVol(fIsOverlapVol);
  fRTPC->SetIsSrcPb(fIsSrcPb);
  fRTPC->Construct(fWorldLogic);
  //fDA = new DetectorArray(this);
  //fDA->ReadParameters(fDAFile);
  //fDA->SetIsOverlapVol(fIsOverlapVol);
  //fDA->Construct(fWorldLogic);
  // Visualization attributes
  //
  fWorldLogic->SetVisAttributes (G4VisAttributes::Invisible);
  return fWorldPhysi;
}

//----------------------------------------------------------------------------
void DetectorConstruction::UpdateGeometry()
{
  G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
}
