#include "G4SBSGEMSD.hh"
#include "G4SBSGEMHit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "G4Box.hh"

#include <iostream>

G4SBSGEMSD::G4SBSGEMSD( G4String name, G4String colname )
  : G4VSensitiveDetector(name)
{
    collectionName.insert(colname);
    GEMTrackerIDs.clear();
}

G4SBSGEMSD::~G4SBSGEMSD()
{
}

void G4SBSGEMSD::Initialize(G4HCofThisEvent*)
{
  hitCollection = new G4SBSGEMHitsCollection(SensitiveDetectorName,collectionName[0]); 
}

G4bool G4SBSGEMSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  double edep = aStep->GetTotalEnergyDeposit();
  const G4ParticleDefinition *particle = aStep->GetTrack()->GetParticleDefinition();

  //Get logical volume name for Tracker ID assignment. Must do this before "MoveUpHistory" or it won't work properly
  G4String volname = (aStep->GetPreStepPoint()->GetPhysicalVolume() )->GetName();

  // Charged geantinos always get tracked
  if( edep <= 0.0 && !(particle->GetParticleName()=="chargedgeantino") ) return false;

  // Only return primary electron hits
//  int trid = aStep->GetTrack()->GetParentID();
//  if( trid != 0 ) return false;

//  G4ThreeVector mom = aStep->GetTrack()->GetMomentum();

  G4TouchableHistory* hist = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());

  // G4cout << G4endl << G4endl << "*******************************************************************" << G4endl;
  // G4cout << "GEM hit history depth = " << hist->GetHistory()->GetDepth() << endl;

  // G4cout << "GEMSD: Before MoveUpHistory, volume name = " << hist->GetVolume()->GetName() << 
  //   " prestep volume name = " << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() <<  G4endl;

  // for(G4int i=0; i<=hist->GetHistoryDepth(); i++){
  //   //G4ThreeVector trans = hist->GetTranslation(i);
  //   //G4RotationMatrix *rm = hist->GetRotation(i);
    
  //   G4cout << "at depth " << i << G4endl;
  //   G4cout << "Replica Number: " << hist->GetReplicaNumber(i) << G4endl;
  //   G4cout << "Physical volume = " << hist->GetVolume(i)->GetName() << G4endl;
  //   G4cout << "Logical volume = " << hist->GetVolume(i)->GetLogicalVolume()->GetName() << G4endl;
  //   // G4cout << "Translation (x,y,z)=(" << trans.x()/m << ", " << trans.y()/m << ", " << trans.z()/m << ")" << G4endl;
  //   // G4cout << "Rotation ( xx, xy, xz " << hist->GetRotation(i)->xx() << ", " << hist->GetRotation(i)->xy() << ", " << hist->GetRotation(i)->xz() << ", " << G4endl 
  //   // 	   << "          yx, yy, yz  " << hist->GetRotation(i)->yx() << ", " << hist->GetRotation(i)->yy() << ", " << hist->GetRotation(i)->yz() << ", " << G4endl 
  //   // 	   << "          zx, zy, zz )" << hist->GetRotation(i)->zx() << ", " << hist->GetRotation(i)->zy() << ", " << hist->GetRotation(i)->zz() << ", " << G4endl;

  //   G4AffineTransform transform = hist->GetHistory()->GetTransform(i);
  //   G4cout << "Net translation (x,y,z)=(" << transform.NetTranslation().x()/m << ", " << transform.NetTranslation().y()/m << ", " << transform.NetTranslation().z()/m << ")" << G4endl;
  //   G4cout << "Net rotation:" << G4endl
  // 	   << "| xx, xy, xz |   |" << transform.NetRotation().xx() << ", " << transform.NetRotation().xy() << ", " << transform.NetRotation().xz() << "|" << G4endl
  // 	   << "| yx, yy, yz | = |" << transform.NetRotation().yx() << ", " << transform.NetRotation().yy() << ", " << transform.NetRotation().yz() << "|" << G4endl
  // 	   << "| zx, zy, zz | = |" << transform.NetRotation().zx() << ", " << transform.NetRotation().zy() << ", " << transform.NetRotation().zz() << "|" << G4endl;
  // }      

  //hist->MoveUpHistory();
  
  // G4cout << "GEMSD: After MoveUpHistory, volume name = " << hist->GetVolume()->GetName() << 
  //   " prestep volume name = " << aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName() << G4endl;

  //  This is where the box position will be


  
  // G4AffineTransform bTrans = hist->GetHistory()->GetTransform(  hist->GetHistory()->GetDepth()-1 );
  // double boxz = bTrans.TransformPoint(hist->GetTranslation()).getZ();

  // G4cout << "boxz = " << boxz/m << G4endl;

  G4int copyID = hist->GetReplicaNumber(1);
  //  G4cout << "copyID = " << copyID << endl;

  //  Go to detector box
  //hist->MoveUpHistory();

  double offset = ((G4Box *) hist->GetSolid(2))->GetZHalfLength();
  //boxz += offset;

  G4SBSGEMHit* hit = new G4SBSGEMHit();

  //G4AffineTransform aTrans = hist->GetHistory()->GetTopTransform();
  G4AffineTransform aTrans = hist->GetHistory()->GetTransform(hist->GetHistoryDepth()-2);

  // G4cout << "Top transform:" << G4endl;
  // G4cout << "Net translation (x,y,z)=(" << aTrans.NetTranslation().x()/m << ", " << aTrans.NetTranslation().y()/m << ", " << aTrans.NetTranslation().z()/m << ")" << G4endl;
  // G4cout << "Net rotation:" << G4endl
  // 	 << "| xx, xy, xz |   |" << aTrans.NetRotation().xx() << ", " << aTrans.NetRotation().xy() << ", " << aTrans.NetRotation().xz() << "|" << G4endl
  // 	 << "| yx, yy, yz | = |" << aTrans.NetRotation().yx() << ", " << aTrans.NetRotation().yy() << ", " << aTrans.NetRotation().yz() << "|" << G4endl
  // 	 << "| zx, zy, zz | = |" << aTrans.NetRotation().zx() << ", " << aTrans.NetRotation().zy() << ", " << aTrans.NetRotation().zz() << "|" << G4endl;

  G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();

  pos = aTrans.TransformPoint(pos);

  // for(int i=0; i<=hist->GetHistoryDepth(); i++){
  //   G4AffineTransform boxTrans = hist->GetHistory()->GetTransform(i);
  //   G4ThreeVector trackerpos = boxTrans.TransformPoint(aStep->GetPreStepPoint()->GetPosition());

  //   G4cout << "at depth " << i << ", local trackerpos (x,y,z)=(" << trackerpos.x()/m << ", " << trackerpos.y()/m << ", " << trackerpos.z()/m << ")" << G4endl;
  // }
 // G4ThreeVector mom = aStep->GetPreStepPoint()->GetMomentum();
  //  G4ThreeVector mom = aStep->GetDeltaPosition();
  G4ThreeVector mom = aStep->GetPreStepPoint()->GetMomentumDirection();
  
  G4ThreeVector thisdelta = aStep->GetPreStepPoint()->GetPosition() + aStep->GetDeltaPosition();
  thisdelta = aTrans.TransformPoint(thisdelta);
  thisdelta -= pos;

//  printf("%f %f %f\n", aStep->GetDeltaPosition().x()/mom.x(), aStep->GetDeltaPosition().y()/mom.y(), aStep->GetDeltaPosition().z()/mom.z());
//  printf("%f %f %f\n", mom.x(), mom.y(), mom.z());
  //mom *= aTrans.NetRotation();
  //mom *= aTrans.NetRotation().inverse();
  mom = aTrans.TransformAxis(mom);

  //printf("mom:    %f %f %f\n", mom.unit().x(), mom.unit().y(), mom.unit().z());


  /*
  printf("mom:    %f %f %f\n", mom.unit().x(), mom.unit().y(), mom.unit().z());
  printf("thisdel %f %f %f\n", thisdelta.unit().x(), thisdelta.unit().y(), thisdelta.unit().z());
  */

  //std::cout << aTrans.NetRotation() << std::endl;

//  pos.setZ(  boxz );
  pos.setZ( pos.getZ() + offset );

  /*
  printf("pos = %f %f %f\n", pos.getX()/m, pos.getY()/m, pos.getZ()/m);
  printf("mom = %f %f %f (%f)\n", mom.getX()/GeV, mom.getY()/GeV, mom.getZ()/GeV, mom.mag()/GeV);
  */
  
  map<G4String,int>::iterator trackerID = GEMTrackerIDs.find( volname );
  
  if( trackerID != GEMTrackerIDs.end() ){
    hit->SetTrackerID( trackerID->second );

    //G4cout << "GEM Hit: Assigned tracker ID " << hit->GetTrackerID() << G4endl;
  } else {
    hit->SetTrackerID( -1 );

    //    G4cout << "GEM Hit: Assigned tracker ID " << hit->GetTrackerID() << G4endl;
  }

//  printf("Hit by %d at det %d,  (%f %f %f) %f?\n", trid, copyID, pos.getX()/cm, pos.getY()/cm, pos.getZ()/cm, boxz/cm );

  hit->SetEdep(edep);
  hit->SetPos(pos);
  hit->SetVertex(aStep->GetTrack()->GetVertexPosition());
  hit->SetMID(aStep->GetTrack()->GetParentID());
  hit->SetTrID(aStep->GetTrack()->GetTrackID());
  hit->SetPID(aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding());
  //hit->SetMom(aStep->GetTrack()->GetMomentum().mag());
  hit->SetMom(aStep->GetPreStepPoint()->GetMomentum().mag());
  hit->SetDir(mom.getX()/mom.getZ(), mom.getY()/mom.getZ());
//  hit->SetDir(thisdelta.getX()/thisdelta.getZ(), thisdelta.getY()/thisdelta.getZ());
  hit->SetGEMID(copyID);
  hit->SetBeta( aStep->GetPreStepPoint()->GetBeta() );
  hit->SetHittime( aStep->GetPreStepPoint()->GetGlobalTime() );

  hitCollection->insert( hit );

  return true;
}

void G4SBSGEMSD::EndOfEvent(G4HCofThisEvent*HCE)
{
    G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    HCE->AddHitsCollection( HCID, hitCollection );
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
