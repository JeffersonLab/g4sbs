#include "G4SBSRICHSD.hh"
#include "G4SBSRICHHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"

G4SBSRICHSD::G4SBSRICHSD( G4String name, G4String collname ) : G4VSensitiveDetector(name) {
  collectionName.insert( collname );
}

G4SBSRICHSD::~G4SBSRICHSD(){;}

void G4SBSRICHSD::Initialize( G4HCofThisEvent *HC ){
  static int HCID = -1;
  hitCollection = new G4SBSRICHHitsCollection( SensitiveDetectorName, collectionName[0] );
  
  if( HCID < 0 ){
    HCID = GetCollectionID(0);
  }

  HC->AddHitsCollection( HCID, hitCollection );
}

G4bool G4SBSRICHSD::ProcessHits( G4Step *aStep, G4TouchableHistory* ){
  G4double edep = aStep->GetTotalEnergyDeposit();

  //For the RICH, we only consider optical photons to be part of the hit (this is what we are sensitive to)
  if( aStep->GetTrack()->GetParticleDefinition() != G4OpticalPhoton::OpticalPhotonDefinition() )
    return false;
  
  G4SBSRICHHit *newHit = new G4SBSRICHHit();

  //Get pointer to track:
  G4Track *track = aStep->GetTrack();
  //Get pointers to pre-step point and post-step point:
  G4StepPoint *prestep = aStep->GetPreStepPoint();
//  G4StepPoint *poststep = aStep->GetPostStepPoint();

  //Where to get all the information needed for the hit? 
  //From get methods of G4Step, G4Track, and G4StepPoint

  newHit->SetTrackID( track->GetTrackID() );
  newHit->SetTrackPID( track->GetParticleDefinition()->GetPDGEncoding() );
  newHit->SetMotherID( track->GetParentID() ); 
  //Getting the rest of the Mother particle information must wait until the end of the event we can access G4Event->G4TrajectoryContainer

  newHit->SetVertex( track->GetVertexPosition() );
  newHit->SetVertexDirection( track->GetVertexMomentumDirection() );

  //Record where was the optical photon produced:
  //Three cases are interesting: Aerogel tiles, RICH box containing C4F10 gas, or UVT-lucite of aerogel exit window
  G4String namevol_origin = track->GetLogicalVolumeAtVertex()->GetName();
  int origin_flag = 0; //default
  if( namevol_origin.contains("Aerogel_tile_log") ){  //Aerogel
    origin_flag = 1;
  } else if (namevol_origin.contains("RICHbox_log") ){ //C4F10
    origin_flag = 2;
  } else if (namevol_origin.contains("Aero_exitwindow") ){ //UVT lucite
    origin_flag = 3;
  }

  newHit->SetOriginVol( origin_flag );

  newHit->Setdx( aStep->GetStepLength() );
  newHit->SetTime( prestep->GetGlobalTime() );
  newHit->SetEdep( edep );
  newHit->Setenergy( prestep->GetTotalEnergy() );

  //Let's change this so that it refers to the local position and direction of the hit:
  
  G4AffineTransform RICH_atrans = ( (G4TouchableHistory*) prestep->GetTouchable() )->GetHistory()->GetTopTransform();

  G4ThreeVector pos = prestep->GetPosition();
  
  pos = RICH_atrans.TransformPoint(pos); 

  G4ThreeVector mom = prestep->GetMomentumDirection();

  mom *= (RICH_atrans.NetRotation()).inverse();

  newHit->SetPos( pos );
  newHit->SetDirection( mom );

  //Get PMT identifying information:

  int PMTno;

  newHit->SetPMTnumber( PMTno = prestep->GetPhysicalVolume()->GetCopyNo() );

  // G4cout << "Hit RICH PMT number " << newHit->GetPMTnumber() << G4endl;

  // G4cout << "Physical volume name = " << prestep->GetPhysicalVolume()->GetName() << G4endl;

  newHit->Setrownumber( newHit->calc_row( PMTno ) );
  newHit->Setcolnumber( newHit->calc_col( PMTno ) );

  newHit->SetLogicalVolume( prestep->GetPhysicalVolume()->GetLogicalVolume() );

  hitCollection->insert( newHit );

  return true;
}

void G4SBSRICHSD::EndOfEvent( G4HCofThisEvent* ){
  ;
}

void G4SBSRICHSD::clear()
{
} 

void G4SBSRICHSD::DrawAll()
{
} 

void G4SBSRICHSD::PrintAll()
{
} 

