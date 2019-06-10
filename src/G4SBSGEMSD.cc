#include "G4SBSGEMSD.hh"
#include "G4SBSGEMHit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//#include "G4Transportation.hh"

#include "G4Box.hh"

#include <iostream>

G4SBSGEMSD::G4SBSGEMSD( G4String name, G4String colname )
  : G4VSensitiveDetector(name)
{
    collectionName.insert(colname);
    //GEMTrackerIDs.clear();

    detmap.depth = 1;
    detmap.SDname = name;
    detmap.clear();

    SDtracks.Clear();
    SDtracks.SetSDname(name);
}

G4SBSGEMSD::~G4SBSGEMSD()
{
}

void G4SBSGEMSD::Initialize(G4HCofThisEvent*)
{
  hitCollection = new G4SBSGEMHitsCollection(SensitiveDetectorName,collectionName[0]);
  SDtracks.Clear();
}

G4bool G4SBSGEMSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{

  // G4cout << "Transporation magnetic moment enabled = "
  // 	 << G4Transportation::EnableUseMagneticMoment(true) << G4endl;
  
  double edep = aStep->GetTotalEnergyDeposit();
  const G4ParticleDefinition *particle = aStep->GetTrack()->GetParticleDefinition();
  
  // Charged geantinos always get tracked
  if( edep <= 0.0 && !(particle->GetParticleName()=="chargedgeantino") ) return false;

  //Used to retrieve coordinate transformations relevant to spectrometer coordinate system:
  G4TouchableHistory* hist = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

  //The sensitive region of the GEM is positioned inside a mother GEM chamber, this is why plane number is one level up the geometry hierarchy
  //from the sensitive region:
  G4int copyID = hist->GetReplicaNumber(1);
 
  //Each GEM chamber is placed inside a mother box, the entrance of which is defined as the Z position for tracking 
  double offset = ((G4Box *) hist->GetSolid(2))->GetZHalfLength();

  G4SBSGEMHit* hit = new G4SBSGEMHit();

  //We are interested in the coordinate of this hit relative to the mother box that contains all the individual GEM planes in a given tracker
  //This is two levels up the geometry hierarchy, meaning that we require the transformation two levels down from the top:
  G4AffineTransform aTrans = hist->GetHistory()->GetTransform(hist->GetHistoryDepth()-2);

  G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition(); //global position of prestep point
  G4ThreeVector outpos = aStep->GetPostStepPoint()->GetPosition(); //global position of poststep point
  G4ThreeVector gpos = pos;// variable for global position
  
  pos = aTrans.TransformPoint(pos); //local position in "tracker box"
  outpos = aTrans.TransformPoint(outpos); //local position in "tracker box"

  // Track Polarization in the same manner as pos
  G4ThreeVector polarization = aStep->GetPreStepPoint()->GetPolarization();
  polarization = aTrans.TransformAxis(polarization);

  if( aStep->GetTrack()->GetParentID() == 0 && aStep->GetTrack()->GetTrackID() == 2 ){
    //G4cout << "Pre-step point, polarization = " << polarization << G4endl;
  }
  //Momentum direction in global coordinates:
  G4ThreeVector mom = aStep->GetPreStepPoint()->GetMomentumDirection();
  
  G4ThreeVector thisdelta = aStep->GetPreStepPoint()->GetPosition() + aStep->GetDeltaPosition();
  thisdelta = aTrans.TransformPoint(thisdelta);
  thisdelta -= pos;

//  printf("%f %f %f\n", aStep->GetDeltaPosition().x()/mom.x(), aStep->GetDeltaPosition().y()/mom.y(), aStep->GetDeltaPosition().z()/mom.z());
//  printf("%f %f %f\n", mom.x(), mom.y(), mom.z());
  //mom *= aTrans.NetRotation();
  //mom *= aTrans.NetRotation().inverse();
  mom = aTrans.TransformAxis(mom); //Momentum direction in local coordinate system of "tracker box"

  pos.setZ( pos.getZ() + offset ); //Offset local z position by half of "tracker box" z width (so that z position is relative to "front" of tracker box)
  outpos.setZ( outpos.getZ() + offset ); //Offset local z position by half of "tracker box" z width (so that z position is relative to "front" of tracker box)

  // printf("pos x = %f y = %f z = %f\n", pos.x(), pos.y(), pos.z());
  // printf("outpos x = %f y = %f z = %f\n", outpos.x(), outpos.y(), outpos.z());
  
  hit->SetEdep(edep); //energy deposition
  hit->SetPos(pos); //position in tracker box coordinates
  hit->SetOutPos(outpos); //position in tracker box coordinates
  hit->SetGlobalPos(gpos); //position in lab frame
  hit->SetPolarization(polarization); //polarization
  hit->SetVertex(aStep->GetTrack()->GetVertexPosition()); //vertex location in global coordinates of particle that caused hit
  hit->SetMID(aStep->GetTrack()->GetParentID()); 
  hit->SetTrID(aStep->GetTrack()->GetTrackID());
  hit->SetPID(aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding());
  hit->SetMom(aStep->GetPreStepPoint()->GetMomentum().mag()); //Momentum of particle that caused hit prior to the step
  hit->SetDir(mom.getX()/mom.getZ(), mom.getY()/mom.getZ());
//  hit->SetDir(thisdelta.getX()/thisdelta.getZ(), thisdelta.getY()/thisdelta.getZ());
  hit->SetGEMID(copyID);
  hit->SetBeta( aStep->GetPreStepPoint()->GetBeta() ); //v/c of particle prior to the step
  hit->SetHittime( aStep->GetPreStepPoint()->GetGlobalTime() );

  G4Track *aTrack = aStep->GetTrack();

  
  hit->SetOTrIdx( SDtracks.InsertOriginalTrackInformation( aTrack ) );
  hit->SetPTrIdx( SDtracks.InsertPrimaryTrackInformation( aTrack ) ); 
  hit->SetSDTrIdx( SDtracks.InsertSDTrackInformation( aTrack ) );
 
  hitCollection->insert( hit );

  return true;
}

void G4SBSGEMSD::EndOfEvent(G4HCofThisEvent*HCE)
{
    G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    HCE->AddHitsCollection( HCID, hitCollection );

    //G4cout << "Adding hit collection " << collectionName[0] << " to HCE, HCID = " << HCID << G4endl;
}

void G4SBSGEMSD::clear()
{
} 

void G4SBSGEMSD::DrawAll()
{
} 

void G4SBSGEMSD::PrintAll()
{
} 
