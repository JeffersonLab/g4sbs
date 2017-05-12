// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class ArraySD
// Save the interesting parameter from the sensitive volumes
// ie the active regions of the detectors
// 03/05/14 JRMA

#include "G4RunManager.hh"
#include "ArraySD.hh"
#include "ArrayHit.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RootIO.hh"

#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4SDManager.hh"

#include "stdio.h"

ArraySD::ArraySD(G4String name, G4int Nelements, G4double* parm):
  G4VSensitiveDetector(name)
{
  collectionName.insert(G4String("ArraySDHits")+name);
  fCollection=NULL;

  fNelements=Nelements+1;//numbering starts from 1 not 0
  fHitID=new G4int[fNelements];
  for(G4int i=0;i<fNelements;i++)fHitID[i]=-1;
  fHits=new G4int[fNelements];
  for(G4int i=0;i<fNelements;i++)fHits[i]=0;
 
  fNhits=0;
  fHCID=-1;
  fEthr = 1.0;
  fTmax = 1000.0; // default max time to track
  pdgCode = trID = parID = stepID = 0;
  Tr1PDG = pTrID = nTrID = gTrID = -1;
  fZst = parm[0];
  fRHe1 = parm[1];
  fRHe2 = parm[2];
  fZPixG = (G4int)parm[3];
  fdZPixG = 2*fZst/fZPixG;
  fPhiPixG = (G4int)parm[4];
  //fdPhiPixG = 6.283185307/fPhiPixG;
  fdPhiPixG = twopi/fPhiPixG;
  fPix = new G4double[fZPixG*fPhiPixG];
  fRio = NULL;
  fIsStore = 1;
}

ArraySD::~ArraySD()
{

}


void ArraySD::Initialize(G4HCofThisEvent* HCE)
{
  //Hits collections are deleted in G4RunManager, DoEventLoop
  //This calls StackPreviousEvent which deletes the G4Event
  //This deletes the hit collection for the event and thus ArrayHit s
  //G4cout<<"void ArraySD::Initialize "<<SensitiveDetectorName<<G4endl;
  fCollection = new ArrayHitsCollection(SensitiveDetectorName,collectionName[0]);
  if(fHCID<0)
    fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  HCE->AddHitsCollection(fHCID,fCollection);
  for (G4int i=0;i<fNelements;i++){
    ArrayHit* hit = new ArrayHit;
    fCollection->insert(hit);
  }

  if( fIsStore ){
    if(!fRio){
      EventAction* eva = (EventAction*)
	(G4RunManager::GetRunManager()->GetUserEventAction());
      fRio = eva->GetRootIO();
      if(!fRio){
	fIsStore = 0;
	return;
      }
      Hparm* pixParm = new Hparm{fZPixG,-fZst,fZst,fPhiPixG,-pi,pi,0,0,0};
      fRio->SetH2(2,pixParm);
      fRio->SetH2(3,pixParm);
    }
  }
}


G4bool ArraySD::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{ 
  // 22 gamma, 11 electron, 2212 proton, 2112 neutron
  // Check 1st that track is still alive
  G4Track* aTrack = aStep->GetTrack();
  pdgCode = aTrack->GetDynamicParticle()->GetPDGcode();
  G4LorentzVector p4 = aTrack->GetDynamicParticle()->Get4Momentum();
  trID = aTrack->GetTrackID();
  if(trID == 1) Tr1PDG = pdgCode;
  parID = aTrack->GetParentID();
  stepID = aTrack->GetCurrentStepNumber();
  G4ThreeVector hpos = aStep->GetPreStepPoint()->GetPosition();
  G4int isgPar = 0;
  static G4double dedxT;
  G4double dedxM = 0;
  static G4int dedxN;
  G4int ispPar = 0; 
  if( parID == gTrID ){
    //printf("tra %d %d %d %d\n",pdgCode,trID,parID,stepID);
    //printf("ggg %d %d %d %d\n",pdgCode,trID,parID,stepID);
    isgPar = 1;
  }
  if( parID == pTrID ){
    //printf("tra %d %d %d %d\n",pdgCode,trID,parID,stepID);
    ispPar = 1;
  }
  if( pdgCode == 22 ){
    gTrID = trID;
    if(fRio){
      fRio->WriteH3(3,hpos.x(),hpos.y(),hpos.z(),p4.e());
      fRio->WriteH3(4,p4.x(),p4.y(),p4.z(),1);
      fRio->WriteH2(0,p4.z(),p4.perp());
      if(stepID == 1)fRio->WriteH1(0,p4.e());
    }
  }
  else if( pdgCode == 2212 ){
    if(pTrID == -1) pTrID = trID;
    if( fRio && (stepID == 1) ){
      fRio->WriteH3(5,hpos.x(),hpos.y(),hpos.z(),p4.e());
      fRio->WriteH2(1,p4.z(),p4.perp());
      fRio->WriteH1(1,p4.e());
      if((parID == 1) && (Tr1PDG == 11)) fRio->WriteH1(2,p4.e());
    }
  }
  else if( pdgCode != 11 ){
    nTrID = trID;
    //printf("nnn %d %d %d %d\n",pdgCode,trID,parID,stepID);
  }
  if( aTrack->GetTrackStatus() != fAlive )
    return false;
  // Check that there has been some energy deposit
  G4double edep = aStep->GetTotalEnergyDeposit();     // energy
  if ((edep/keV == 0.)){
    return false;
  }
  //G4double edepB = GetEffectiveEnergyDeposit(aStep);  // inc Birks atten
  G4double tHit = aStep->GetTrack()->GetGlobalTime(); // global time
  // This TouchableHistory is used to obtain the physical volume
  // of the hit
  G4TouchableHistory* theTouchable
    = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* volume=theTouchable->GetVolume();
  G4int id;
  // Get element copy number...if it exceeds EHCAL_CopyNo then the active
  // block is a segmented calorimeter (with 1 readout channel)
  // In that case get the element number from the mother logical volume
  id = volume->GetCopyNo();
  //  G4cout<<volume->GetName()<<" id "<<id <<" edep "<<edep/MeV<<G4endl;
  G4double x = hpos.x();
  G4double y = hpos.y();
  G4double z = hpos.z();
  G4double r = sqrt(x*x+y*y);
  if( id == 0 ){
    if( r > fRHe1 )
      id = 1;
  }
  if( (pdgCode == 2212) || (ispPar) ){
    id += 5;
  }
  G4double mb = 0;
  if( (id == 1)||(id == 6) ){
    if( fHitID[id] == -1){
      dedxT = 0;
      dedxN = 0;
    }
    G4ThreeVector hpa = hpos;
    G4ThreeVector hpb = aStep->GetPostStepPoint()->GetPosition();
    hpa.setRho(fRHe2);
    hpb.setRho(fRHe2);
    G4ThreeVector dhb = hpb - hpa;
    mb = dhb.mag();
    dedxT += edep/mb;
    dedxN++;
    dedxM = dedxT/dedxN;
  }
  if( fRio ){
    fRio->WriteH3(0,x,y,z,edep*MeV);
    if( isgPar ){
      fRio->WriteH3(1,x,y,z,edep*MeV);
    }
    if( pdgCode == 2212 ){
      fRio->WriteH3(2,x,y,z,edep*MeV);
    }
    //fRio->WriteH3(2,r,z,edep*MeV,1);
  }
  ArrayHit* xHit = (*fCollection)[id];
  if (fHitID[id]==-1){
    // Its a new hit...store start paramaters   
    xHit->SetPos(hpos);
    xHit->SetID(id);
    xHit->SetTime(tHit);
    //    Hit->SetPol(*PolNN);
    fHitID[id] = 0;
    fHits[fNhits++]=id;
  }
  // do following regardless of new or used hit buffer
  xHit->AddEnergy(edep);
  xHit->AddEnergyB(mb);
  xHit->SetDedx(dedxM);  // save mean dedx
  if( (id == 1)||(id == 6) ){
    G4ThreeVector hpb = aStep->GetPostStepPoint()->GetPosition();
    if(!xHit->AddStep(trID,parID,stepID,edep,&p4,&hpos,&hpb))
      //aTrack->SetTrackStatus(fStopAndKill);
      return false;
  }
  return true;
}

void ArraySD::EndOfEvent(G4HCofThisEvent*)
{
  /*
  ArrayHit* xHit;
  for(G4int i=0; i<fNhits; i++){
    xHit =  (*fCollection)[fHitID[fHits[i]]];
    xHit->SortSH(fEthr);;
  }
  // G4cout<<"EndOfEvent( "<<fHCID<<" "<<fCollection<<" "<<fNhits<<G4endl;
  if(fHCID<0) fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  if(fNhits>0)  HCE->AddHitsCollection(fHCID,fCollection);
  //G4cout<<"EndOfEvent( "<<G4endl;
 
  //reset hit arrays
  for (G4int i=0;i<fNhits;i++) 
    {
      fHitID[fHits[i]]=-1;
      fHits[i]=0;
    }
  fNhits=0;
  */
  pdgCode = trID = parID = stepID = 0;
  Tr1PDG = pTrID = nTrID = gTrID = -1;
  //printf("End\n");
}

void ArraySD::clear()
{} 

void ArraySD::DrawAll()
{} 

void ArraySD::PrintAll()
{} 

G4double ArraySD::GetEffectiveEnergyDeposit(const G4Step* aStep)
{
  // Weight energy deposit by scintillation-light quenching effects
  // Scintillation response is generally non-linear with energy
  // According to Birks it depends on de/dx
  //
  G4Material* material = aStep->GetTrack()->GetMaterial();
  G4double birks       = material->GetIonisation()->GetBirksConstant();
  G4double destep      = aStep->GetTotalEnergyDeposit();
  if( (!birks) || (!destep) ) return destep;
  G4double charge      = aStep->GetTrack()->GetDefinition()->GetPDGCharge();
  if( !charge ) return destep;
  G4double stepl       = aStep->GetStepLength();  
  if( !stepl ) return destep;
  //
  G4double response = destep/(1. + birks*destep/stepl);
  return response;
}









