#include "G4SBSmTPCSD.hh"
#include "G4SBSmTPCHit.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"

static const double PI=3.1415926535;

G4SBSmTPCSD::G4SBSmTPCSD( G4String name, G4String collname ) : G4VSensitiveDetector(name) {
  collectionName.insert( collname );
  detmap.SDname = name;
  detmap.clear();

  // temp from CAL SD
  fHitTimeWindow = 1000.0*ns; //"safe" default value for a calorimeter;
  fEnergyThreshold = 0.0*keV; //"safe" default value for a calorimeter;
  fNTimeBins = 500; 
}

G4SBSmTPCSD::~G4SBSmTPCSD(){;}

void G4SBSmTPCSD::Initialize( G4HCofThisEvent *HC ){
  G4int HCID = -1;
  hitCollection = new G4SBSmTPCHitsCollection( SensitiveDetectorName, collectionName[0] );
  
  if( HCID < 0 ){
    HCID = GetCollectionID(0);
  }

  HC->AddHitsCollection( HCID, hitCollection );
}

G4bool G4SBSmTPCSD::ProcessHits( G4Step *aStep, G4TouchableHistory* ){

  G4double edep = aStep->GetTotalEnergyDeposit();
  G4int pid = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();

  if( edep <= 0.0 || pid == 0 ) return false;
  // if( edep <= 0.0 ) return false;

  G4SBSmTPCHit *newHit = new G4SBSmTPCHit();

  //Get pointer to track:
  G4Track *track = aStep->GetTrack();

  // alternative option to stepping action,
  // to ensure protons are killed when hitting gem and readout
  // is to put in every SD:
  // if(theTrack->GetMaterial()->GetName() == "Kryptonite"){ theTrack->SetTrackStatus(fStopAndKill); }

  //Get pointers to pre-step point and post-step point:
  G4StepPoint *prestep = aStep->GetPreStepPoint();
  G4TouchableHistory* hist = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
//  G4StepPoint *poststep = aStep->GetPostStepPoint();

  newHit->SetTrackID( track->GetTrackID() );
  newHit->SetMotherID( track->GetParentID() );
  newHit->SetTrackPID( pid );
  newHit->SetTime( prestep->GetGlobalTime() );
  newHit->SetVertex( aStep->GetTrack()->GetVertexPosition() );
  newHit->SetEdep( edep );
  newHit->SetEnergy( prestep->GetTotalEnergy() );
  newHit->Setdx( aStep->GetStepLength() );

  G4ThreeVector pos = prestep->GetPosition();
  //Pos refers to local position:
  newHit->SetLPos( pos );
  // local momentum
  G4ThreeVector mom = prestep->GetMomentum();
  newHit->SetLMom( mom );
  // set direction
  newHit->SetLDirection( mom );

  //Let's change this so that it refers to the local position and direction of the hit: 
  G4AffineTransform mTPC_atrans = ( (G4TouchableHistory*) prestep->GetTouchable() )->GetHistory()->GetTopTransform();

  //Lpos refers to global position of the hit:
  pos = mTPC_atrans.TransformPoint(pos); 
  newHit->SetPos( pos );
  // and local momentum
  //mom *= (mTPC_atrans.NetRotation()).inverse();
  mom = track->GetVertexMomentumDirection()*(track->GetVertexKineticEnergy()+track->GetParticleDefinition()->GetPDGMass());
  newHit->SetMom( mom );
  //and local momentum direction
  newHit->SetLDirection( mom );


  //Get chamber identifying information:
  int cell;
  //newHit->SetPMTnumber( PMTno = prestep->GetPhysicalVolume()->GetCopyNo() );
  newHit->SetCell( cell = hist->GetVolume( detmap.depth )->GetCopyNo() );

  newHit->SetCellCoord( detmap.LocalCoord[cell] );
  //mTPC_atrans is the transformation which, when applied to a position in the global coordinate system, gives the local coordinates of a point; 
  //To go the other way, is it sufficient to apply the inverse of the transformation to the (local) point (0,0,0)?
  G4AffineTransform mTPC_atrans_inverse = mTPC_atrans.Inverse();
  newHit->SetGlobalCellCoord( mTPC_atrans_inverse.TransformPoint( G4ThreeVector(0,0,0) ) );
  // i think this is returning a strange value for z coord of hits



  //---------------- get strips and change in z for drift time calc
  // This is just a 1st order rough estimate
  // must implement properly
  G4StepPoint *poststep = aStep->GetPostStepPoint();
  G4ThreeVector postpos = poststep->GetPosition();
  // pre step is pos
  G4double Deltaz = abs(pos.getZ() - postpos.getZ());
  // G4double SegmentAngle = pos.angle(postpos);
  // G4double rXY = pos.perp();
  // G4double Segment = rXY * SegmentAngle;
  // G4double StripLength = sqrt(5.0*5.0 + 5.0*5.0);//5mm strip size, g4 positions in mm
  // G4int NStrip = (G4int) ceil(Segment/StripLength);
  // G4cout << "Pre Step is " << pos.getX() << " " << pos.getY() << " " << pos.getZ() <<
  //   " Post Step is " << postpos.getX() << " " << postpos.getY() << " " << postpos.getZ() << G4endl;
  // G4cout << "SegmentAngle " << SegmentAngle << G4endl;
  // G4cout << "Segment " << Segment << G4endl;
  // G4cout << "StripLength " << StripLength << G4endl;
  // G4cout << "NStrip " << Segment/StripLength << G4endl;
  // G4cout << "NStrip " << NStrip << G4endl;
  // G4cout << "Deltaz = " << Deltaz << G4endl;

  G4ThreeVector vDiff = postpos - pos;
  G4double rXY2 = vDiff.perp();
  G4int NStrip = (G4int) ceil(rXY2/5.0);
  // G4cout << "rXY2 " << rXY2 << G4endl;
  // G4cout << "NStrip2 " << NStrip2 << G4endl;


  newHit->SetZTravel( Deltaz );
  newHit->SetNStrips( NStrip );


  //---------------- get strips




  newHit->SetLogicalVolume( prestep->GetPhysicalVolume()->GetLogicalVolume() );

  newHit->SetOTrIdx( SDtracks.InsertOriginalTrackInformation( track ) );
  newHit->SetPTrIdx( SDtracks.InsertPrimaryTrackInformation( track ) );
  newHit->SetSDTrIdx( SDtracks.InsertSDTrackInformation( track ) );

  hitCollection->insert( newHit );

  return true;
}

void G4SBSmTPCSD::EndOfEvent( G4HCofThisEvent* ){
}

void G4SBSmTPCSD::clear()
{
} 

void G4SBSmTPCSD::DrawAll()
{
} 

void G4SBSmTPCSD::PrintAll()
{
} 
