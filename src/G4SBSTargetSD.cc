#include "G4SBSTargetSD.hh"
//______________________________________________________________________________
G4SBSTargetSD::G4SBSTargetSD(
                            const G4String& name,
                            const G4String& hitsCollectionName)
 : G4VSensitiveDetector(name),
   fHitsCollection(nullptr)
{
  collectionName.insert(hitsCollectionName);
}
//______________________________________________________________________________
G4SBSTargetSD::~G4SBSTargetSD()
{

}
//______________________________________________________________________________
void G4SBSTargetSD::Initialize(G4HCofThisEvent* hce)
{
  // Create hits collection
  fHitsCollection = new G4SBSTargetHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection in hce
  auto hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection );

}
//______________________________________________________________________________
G4bool G4SBSTargetSD::ProcessHits(G4Step* step,G4TouchableHistory*)
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
  G4SBSTargetHit *hit = new G4SBSTargetHit();

  auto touchable = (step->GetPreStepPoint()->GetTouchable());

  // Used to retrieve coordinate transformations relevant to spectrometer coordinate system:
  G4TouchableHistory* hist = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());

  if(!hit){
    G4ExceptionDescription msg;
    msg << "Cannot access hit! ";
    G4Exception("G4SBSTargetSD::ProcessHits()","MyCode0004", FatalException, msg);
  }

  // grab info from step 
  G4ThreeVector mom = step->GetPreStepPoint()->GetMomentum();
  G4double E        = step->GetPreStepPoint()->GetTotalEnergy();
  G4double pMag     = step->GetPreStepPoint()->GetMomentum().mag(); // Momentum of particle that caused hit prior to the step

  G4int trackID     = step->GetTrack()->GetTrackID();
  G4int pid         = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  G4int mid         = step->GetTrack()->GetParentID();

  G4double beta     = step->GetPreStepPoint()->GetBeta();       // v/c of particle *prior* to step 
  G4double hitTime  = step->GetPreStepPoint()->GetGlobalTime(); // time right before the current step   

  // set hit details
  hit->SetPID(pid);
  hit->SetMID(mid);
  hit->SetTrackID(trackID);
  hit->SetMomentum(mom); 
  hit->SetMomentumMag(pMag); 
  hit->SetEdep(edep);
  hit->SetTrackLength(stepLength);
  hit->SetTotalEnergy(E);
  hit->SetBeta(beta);
  hit->SetHitTime(hitTime);

  // extra stuff 
  // G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();     // in the lab coordinates
  // hit->SetPos(pos);
  // // transform position into local coordinates of detector 
  // G4AffineTransform aTrans = hist->GetHistory()->GetTopTransform();
  // pos = aTrans.TransformPoint(pos);
  // hit->SetLabPos(pos);

  // now append to vector 
  fHitsCollection->insert(hit);

  return true;
}
//______________________________________________________________________________
void G4SBSTargetSD::EndOfEvent(G4HCofThisEvent*)
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

