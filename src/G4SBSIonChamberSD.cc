#include "G4SBSIonChamberSD.hh"

//______________________________________________________________________________
G4SBSIonChamberSD::G4SBSIonChamberSD(
                            const G4String& name,
                            const G4String& hitsCollectionName)
 : G4VSensitiveDetector(name),
   fHitsCollection(nullptr)
//   fNofCells(nofCells)
{
  collectionName.insert(hitsCollectionName);
}
//______________________________________________________________________________
G4SBSIonChamberSD::~G4SBSIonChamberSD()
{

}
//______________________________________________________________________________
void G4SBSIonChamberSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection = new G4SBSICHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  auto hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );

  // FIXME: This doesn't seem right (for us)
  // mimic g4sbs implementation for the GEMs -- will resize this vector dynamically  
  // Create hits
  // fNofCells for cells + one more for total sums 
  // for (G4int i=0; i<fNofCells+1; i++ ) {
  //   fHitsCollection->insert(new G4SBSICHit());
  // }

}
//______________________________________________________________________________
G4bool G4SBSIonChamberSD::ProcessHits(G4Step* step,G4TouchableHistory*)
{

  // energy deposit
  auto edep = step->GetTotalEnergyDeposit();

  // step length
  G4double stepLength = 0.;
  if ( step->GetTrack()->GetDefinition()->GetPDGCharge() != 0. ) {
    stepLength = step->GetStepLength();
  }

  if ( edep==0. && stepLength == 0. ) return false;

  // create a hit 
  G4SBSICHit *hit = new G4SBSICHit();

  auto touchable = (step->GetPreStepPoint()->GetTouchable());

  // Used to retrieve coordinate transformations relevant to spectrometer coordinate system:
  G4TouchableHistory* hist = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());

  // Get cell/layer ID 
  // auto layerNo = touchable->GetReplicaNumber(1);      // doesn't work 
  // auto layerNo = touchable->GetVolume()->GetCopyNo(); // works! 

  // G4int copyNo            = touchable->GetVolume()->GetCopyNo();
  // G4int motherCopyNo      = touchable->GetVolume(1)->GetCopyNo();
  // G4int grandMotherCopyNo = touchable->GetVolume(2)->GetCopyNo();
  // char msg[200];
  // sprintf(msg,"[G4SBSIonChamberSD]: layerNo = %d, copyNo = %d, mother copyNo = %d, gma copyNo = %d",
  //         layerNo,copyNo,motherCopyNo,grandMotherCopyNo); 
  // std::cout << msg << std::endl;

  // FIXME: Note that the layerNo != hit number!  
  // Get hit accounting data for this cell
  // auto hit = (*fHitsCollection)[layerNo];  // WARNING: This is a HitCollection! i.e., vector of hits 
  if ( ! hit ) {
    G4ExceptionDescription msg;
    // msg << "Cannot access hit " << layerNo; 
    msg << "Cannot access hit! ";
    G4Exception("G4SBSIonChamberSD::ProcessHits()",
      "MyCode0004", FatalException, msg);
  }

  // Get hit for total accounting
  // auto hitTotal = (*fHitsCollection)[fHitsCollection->entries()-1];

  // grab info from step 
  G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();     // in the lab coordinates
  G4ThreeVector mom = step->GetPreStepPoint()->GetMomentum();
  G4double E        = step->GetPreStepPoint()->GetTotalEnergy();
  G4double pMag     = step->GetPreStepPoint()->GetMomentum().mag(); // Momentum of particle that caused hit prior to the step

  G4int trackID     = step->GetTrack()->GetTrackID();
  G4int pid         = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  G4int mid         = step->GetTrack()->GetParentID();

  G4double beta     = step->GetPreStepPoint()->GetBeta();       // v/c of particle *prior* to step 
  G4double hitTime  = step->GetPreStepPoint()->GetGlobalTime(); // time right before the current step   

  // hit position 
  hit->SetPos(pos);
  // transform position into local coordinates of detector 
  G4AffineTransform aTrans = hist->GetHistory()->GetTopTransform();
  pos = aTrans.TransformPoint(pos);

  // set hit details
  hit->SetPID(pid);
  hit->SetMID(mid);
  hit->SetTrackID(trackID);
  hit->SetLabPos(pos);
  hit->SetMomentum(mom); 
  hit->SetMomentumMag(pMag); 
  hit->SetEdep(edep);
  hit->SetTrackLength(stepLength);
  hit->SetTotalEnergy(E);
  hit->SetBeta(beta);
  hit->SetHitTime(hitTime);
  // hit->SetLayer(layerNo);      
  // hit->Add(edep, stepLength);  

  // hitTotal->SetPos(pos);     
  // hitTotal->SetLabPos(pos);  
  // hitTotal->Add(edep, stepLength);
  // hitTotal->SetTotalEnergy(E);

  // now append to vector 
  fHitsCollection->insert(hit);

  return true;
}
//______________________________________________________________________________
void G4SBSIonChamberSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) {
     auto nofHits = fHitsCollection->entries();
     G4cout
       << G4endl
       << "-------->Hits Collection: in this event they are " << nofHits
       << " hits in the tracker chambers: " << G4endl;
     for ( std::size_t i=0; i<nofHits; ++i ) (*fHitsCollection)[i]->Print();
  }
}

