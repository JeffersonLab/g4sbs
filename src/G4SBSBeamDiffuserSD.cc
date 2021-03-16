#include "G4SBSBeamDiffuserSD.hh"
//______________________________________________________________________________
G4SBSBeamDiffuserSD::G4SBSBeamDiffuserSD(
      const G4String& name,
      const G4String& hitsCollectionName)
   : G4VSensitiveDetector(name),
   fHitsCollection(nullptr)
{
   collectionName.insert(hitsCollectionName);
}
//______________________________________________________________________________
G4SBSBeamDiffuserSD::~G4SBSBeamDiffuserSD()
{

}
//______________________________________________________________________________
void G4SBSBeamDiffuserSD::Initialize(G4HCofThisEvent* hce)
{
   // Create hits collection
   fHitsCollection = new G4SBSBDHitsCollection(SensitiveDetectorName, collectionName[0]);

   // Add this collection in hce
   auto hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
   hce->AddHitsCollection( hcID, fHitsCollection );
}
//______________________________________________________________________________
G4bool G4SBSBeamDiffuserSD::ProcessHits(G4Step* step,G4TouchableHistory*)
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
   G4SBSBDHit *hit = new G4SBSBDHit();

   // Get plane ID 
   auto touchable = (step->GetPreStepPoint()->GetTouchable());
   auto planeNo   = touchable->GetVolume()->GetCopyNo();

   if ( ! hit ) {
      G4ExceptionDescription msg;
      msg << "Cannot access hit for plane " << planeNo; 
      G4Exception("G4SBSBeamDiffuserSD::ProcessHits()",
	    "MyCode0004", FatalException, msg);
   }

   // grab info from step 
   G4ThreeVector pos = step->GetPreStepPoint()->GetPosition();     // in the lab coordinates
   hit->SetPos(pos);                // position in detector coordinates 

   G4ThreeVector mom = step->GetPreStepPoint()->GetMomentum();
   G4double E        = step->GetPreStepPoint()->GetTotalEnergy();

   G4int trackID     = step->GetTrack()->GetTrackID();
   G4int pid         = step->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
   G4int mid         = step->GetTrack()->GetParentID();

   G4double beta     = step->GetPreStepPoint()->GetBeta();        
   G4double hitTime  = step->GetPreStepPoint()->GetGlobalTime(); // time right before the current step   

   // transform position into local coordinates of detector 
   G4TouchableHistory* hist = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());
   G4AffineTransform aTrans = hist->GetHistory()->GetTopTransform();
   pos = aTrans.TransformPoint(pos);
   hit->SetLabPos(pos);             // position in lab coordinates 

   // set all other hit details
   hit->SetPID(pid);                // particle ID 
   hit->SetMID(mid);                // material ID 
   hit->SetTrackID(trackID);        // track number 
   hit->SetPlane(planeNo);          // plane number  
   hit->SetEdep(edep);              // deposited energy  
   hit->SetTrackLength(stepLength); // length of the track 
   hit->SetTotalEnergy(E);          // set total energy
   hit->SetBeta(beta);              // v/c of paraticle
   hit->SetHitTime(hitTime);        // time of hit 

   // now append to vector 
   fHitsCollection->insert(hit);

   return true;
}
//______________________________________________________________________________
void G4SBSBeamDiffuserSD::EndOfEvent(G4HCofThisEvent*)
{
   if ( verboseLevel>1 ) {
      auto nofHits = fHitsCollection->entries();
      G4cout
	 << G4endl
	 << "-------->Beam Diffuser Hit Collection: in this event there are " << nofHits
	 << " hits: " << G4endl;
      for ( std::size_t i=0; i<nofHits; ++i ) (*fHitsCollection)[i]->Print();
   }
}  
