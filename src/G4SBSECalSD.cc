#include "G4SBSECalSD.hh"
#include "G4SBSECalHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"

G4SBSECalSD::G4SBSECalSD( G4String name, G4String collname ) : G4VSensitiveDetector(name) {
  collectionName.insert( collname );
}

G4SBSECalSD::~G4SBSECalSD(){;}

void G4SBSECalSD::Initialize( G4HCofThisEvent *HC ){
  static int HCID = -1;
  hitCollection = new G4SBSECalHitsCollection( SensitiveDetectorName, collectionName[0] );
  
  if( HCID < 0 ){
    HCID = GetCollectionID(0);
  }

  HC->AddHitsCollection( HCID, hitCollection );
}

G4bool G4SBSECalSD::ProcessHits( G4Step *aStep, G4TouchableHistory* ){
  G4double edep = aStep->GetTotalEnergyDeposit();

  //For the ECal, we only consider optical photons to be part of the hit
  if( aStep->GetTrack()->GetParticleDefinition() != G4OpticalPhoton::OpticalPhotonDefinition() )
    return false;
  
  G4SBSECalHit *newHit = new G4SBSECalHit();

  //Get pointer to track:
  G4Track *track = aStep->GetTrack();
  //Get pointers to pre-step point and post-step point:
  G4StepPoint *prestep = aStep->GetPreStepPoint();
//  G4StepPoint *poststep = aStep->GetPostStepPoint();

  //Where to get all the information needed for the hit? 
  //From get methods of G4Step, G4Track, and G4StepPoint

  newHit->SetTrackID( track->GetTrackID() );
  newHit->Setdx( aStep->GetStepLength() );
  newHit->SetTime( prestep->GetGlobalTime() );
  newHit->SetEdep( edep );
  newHit->Setenergy( prestep->GetTotalEnergy() );

  //Let's change this so that it refers to the local position and direction of the hit: 
  G4AffineTransform ECal_atrans = ( (G4TouchableHistory*) prestep->GetTouchable() )->GetHistory()->GetTopTransform();
  G4ThreeVector pos = prestep->GetPosition();
  pos = ECal_atrans.TransformPoint(pos); 
  newHit->SetPos( pos );

  //G4ThreeVector mom = prestep->GetMomentumDirection();
  //mom *= (RICH_atrans.NetRotation()).inverse();
  //newHit->SetDirection( mom );

  //Get PMT identifying information:

  int PMTno;
  newHit->SetPMTnumber( PMTno = prestep->GetPhysicalVolume()->GetCopyNo() );
  newHit->Setrownumber( newHit->calc_row( PMTno ) );
  newHit->Setcolnumber( newHit->calc_col( PMTno ) );

  newHit->SetLogicalVolume( prestep->GetPhysicalVolume()->GetLogicalVolume() );
  hitCollection->insert( newHit );

  return true;
}

void G4SBSECalSD::EndOfEvent( G4HCofThisEvent* ){
  ;
}

void G4SBSECalSD::clear()
{
} 

void G4SBSECalSD::DrawAll()
{
} 

void G4SBSECalSD::PrintAll()
{
} 

