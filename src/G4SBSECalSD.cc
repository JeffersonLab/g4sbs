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
  G4int HCID = -1;
  hitCollection = new G4SBSECalHitsCollection( SensitiveDetectorName, collectionName[0] );
  
  if( HCID < 0 ){
    HCID = GetCollectionID(0);
  }
  
  //  G4cout << "Adding hit collection " << collectionName[0] << " to HCE, HCID = " << HCID << G4endl;

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

  //Pos refers to global position:
  newHit->SetPos( pos );

  //Lpos refers to local position of the hit:
  pos = ECal_atrans.TransformPoint(pos); 
  newHit->SetLPos( pos );

  //G4ThreeVector mom = prestep->GetMomentumDirection();
  //mom *= (RICH_atrans.NetRotation()).inverse();
  //newHit->SetDirection( mom );

  //Get PMT identifying information:

  int PMTno;
  newHit->SetPMTnumber( PMTno = prestep->GetPhysicalVolume()->GetCopyNo() );
  newHit->Setrownumber( detmap.Row[PMTno] );
  newHit->Setcolnumber( detmap.Col[PMTno] );

  newHit->SetCellCoords( detmap.LocalCoord[PMTno] );
  //ECal_atrans is the transformation which, when applied to a position in the global coordinate system, gives the local coordinates of a point; 
  //To go the other way, is it sufficient to apply the inverse of the transformation to the (local) point (0,0,0)?
  G4AffineTransform ECal_atrans_inverse = ECal_atrans.Inverse();
  newHit->SetGlobalCellCoords( ECal_atrans_inverse.TransformPoint( G4ThreeVector(0,0,0) ) ); //We'll see if this works...

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

