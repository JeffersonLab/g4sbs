#include "DVCSCaloSD.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include "G4RunManager.hh"
#include "DVCSPrimaryGeneratorAction.hh"
#include "TCaloEvent.h"
#include "TCaloBlock.h"
#include "TDVCSGlobal2h"
#include "TDVCSEventMC.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

DVCSCaloSD::DVCSCaloSD(G4String name)
:G4VSensitiveDetector(name)
{
 gdvcs->SetWF(0);gdvcs->SetRun(8000);gdvcs->ForceUpdate();
 caloev=new TCaloEvent();

 //This is the tree to write
 tf=new TFile("output.root","RECREATE");
 t=new TTree("T","My tree");
 t->Branch("event_calo.","TCaloEvent",&caloev,16000,2);
 t->Branch("electron","TVector3",&electrowrite,16000,2);
 t->Branch("vertex","TVector3",&vertexwrite,16000,2);


 ev=new TDVCSEvent();
  
 G4String HCname;
 collectionName.insert(HCname="trackerCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DVCSCaloSD::~DVCSCaloSD(){ 

  // Write tree to root file

  tf->Write(); 
      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DVCSCaloSD::Initialize(G4HCofThisEvent* HCE)
{

  trackerCollection = new DVCSCaloHitsCollection
                          (SensitiveDetectorName,collectionName[0]); 
  static G4int HCID = -1;
  if(HCID<0)
    { 
      HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
  HCE->AddHitsCollection( HCID, trackerCollection ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool DVCSCaloSD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  
  if(edep==0.) return false;
  
  DVCSCaloHit* newHit = new DVCSCaloHit();
  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  newHit->SetBlockNb  (aStep->GetPreStepPoint()->GetTouchableHandle()
		       ->GetCopyNumber());
  newHit->SetEdep     (edep);
  newHit->SetPos      (aStep->GetPostStepPoint()->GetPosition());
  trackerCollection->insert( newHit );
  
  // newHit->Print();
  
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DVCSCaloSD::EndOfEvent(G4HCofThisEvent*)
{
   
  ev->Reset();
  
  G4double nblocks=208;
  
  
  //   TH1F *h=new TH1F("h","h",100,0,10);
  //  G4cout <<"Verbose level "<<verboseLevel<<G4endl;
  
  if (verboseLevel>0) 
    { 
      G4int NbHits = trackerCollection->entries();
      //   G4cout << "\n-------->Hits Collection: in this event they are " << NbHits 
      //   << " hits in the calorimeter block : " << G4endl;
      for (G4int i=0;i<NbHits;i++) (*trackerCollection)[i]->Print();
    }
  
  G4int NbHits = trackerCollection->entries();
  G4double eneb[208];
  
  for (G4int i=0;i<nblocks;i++) eneb[i]=0.; 
  for (G4int i=0;i<NbHits;i++)
    {
      eneb[(*trackerCollection)[i]->GetBlockNb()]+=(*trackerCollection)[i]->GetEdep();
      
    }
  G4double GeV=1000.*MeV;
  caloev->Reset();
  for (G4int i=0;i<nblocks;i++)
    {
      
      if(eneb[i]>0.001*GeV)	
	{	
	  //	G4cout<<"Block "<<i<<" - energy "<<G4BestUnit(eneb[i],"Energy")<<G4endl; 
	  G4BestUnit(eneb[i],"Energy");
	  //  h->Fill(eneb[i]);
	  TCaloBlock* block=caloev->AddBlock(i);
	  block->SetBlockEnergy(eneb[i]/GeV);
       	}
    }
  
  /***********************************Clusters***********************/ 
 
  caloev->DoClustering();

  for(G4int k=0;k<caloev->GetNbClusters();k++) caloev->GetCluster(k)->Analyze();
  
  ev->SetCaloEvent(caloev);
  // ev->SetVertex(0.,0.,3.);    // en cm

  if(caloev->GetNbClusters()!=0) t->Fill();
  
  G4RunManager *runman=G4RunManager::GetRunManager();
  DVCSPrimaryGeneratorAction *gen = (DVCSPrimaryGeneratorAction*)runman->GetUserPrimaryGeneratorAction ();
  
  electrowrite=&(gen->electron);
  
  vertexwrite=&(gen->vertex);
  
  caloev->Clear(); //We clear the write event
	  
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}

