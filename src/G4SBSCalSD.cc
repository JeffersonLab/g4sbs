#include "G4SBSCalSD.hh"
#include "G4SBSCalHit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "G4Box.hh"
#include "G4SystemOfUnits.hh"

G4SBSCalSD::G4SBSCalSD( G4String name, G4String colname, G4int recordType)
  : G4VSensitiveDetector(name), fRecordType(recordType)
{
    collectionName.insert(colname);
    detmap.SDname = name;
    detmap.clear();

    fHitTimeWindow = 1000.0*ns; //"safe" default value for a calorimeter;
    fEnergyThreshold = 0.0*keV; //"safe" default value for a calorimeter;
    fNTimeBins = 500; 
}

G4SBSCalSD::~G4SBSCalSD()
{
}

void G4SBSCalSD::Initialize(G4HCofThisEvent*)
{
  hitCollection = new G4SBSCalHitsCollection(SensitiveDetectorName,collectionName[0]); 
}

G4bool G4SBSCalSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  //G4cout << "Processing CAL hits SDname = " << SensitiveDetectorName << G4endl;

  G4int pid = aStep->GetTrack()->GetParticleDefinition()->GetPDGEncoding();
  G4int mid = aStep->GetTrack()->GetParentID();
  double edep = aStep->GetTotalEnergyDeposit();

  // We can record different types of "hits" (0 is the default)
  // 0: All types that deposit non-zero energy (except special G4 particles)
  // 1: Same as 0, but records all primaries regardless if they interacted or not
  // 2: Only primaries
  bool checkBadTrack = (edep <= 0.0 || pid ==0);
  bool checkPrimary = (mid == 0);
  if( fRecordType == 0 ) {
    if( checkBadTrack ) return false;
  } else if ( fRecordType == 1 ) {
    if( checkBadTrack && !checkPrimary ) return false;
  } else if ( fRecordType == 2 ) {
    if( !checkPrimary ) return false;
    edep = 1; // Neutron's don't deposit much energy, but we still care about them
  }

  // We only keep things that deposited energy and are not special Geant4
  // particles, unless we are told to keep all "primary" particles.
  //if( (edep <= 0.0 || pid == 0) && !(fKeepAllPrimaries && mid == 0) ) return false;

//  if( edep <= 0.5*GeV ) return false;

  // Only return primary electron hits
  /*int trid = aStep->GetTrack()->GetParentID();
    if( trid != 0 ) return false; */

  G4SBSCalHit* hit = new G4SBSCalHit();

  G4ThreeVector pos = aStep->GetPreStepPoint()->GetPosition();
  G4ThreeVector mom = aStep->GetPreStepPoint()->GetMomentum();

  G4double E = aStep->GetPreStepPoint()->GetTotalEnergy();

  G4TouchableHistory* hist = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  //G4AffineTransform aTrans = hist->GetHistory()->GetTransform(hist->GetHistory()->GetDepth() );
  G4AffineTransform aTrans = hist->GetHistory()->GetTopTransform();

  hit->SetLabPos(pos); //global coordinates

  pos = aTrans.TransformPoint(pos);

  hit->SetPos(pos); //local coordinates
  
  hit->SetVertex(aStep->GetTrack()->GetVertexPosition());
  hit->SetTime(aStep->GetPreStepPoint()->GetGlobalTime());
  //hit->SetTime(aStep->GetPostStepPoint()->GetLocalTime());

  //  hit->SetTime(aStep->GetPostStepPoint()->GetLocalTime());

  //printf("Hit with E %g\n", edep/MeV);

  hit->SetEdep(edep);
  hit->SetEnergy(E);
  hit->SetLstep( aStep->GetStepLength() );
  hit->SetMomentum(mom);
  hit->SetPID(pid);
  hit->SetTrID(aStep->GetTrack()->GetTrackID());
  hit->SetMID(mid);

  //hit->SetCell( hist->GetVolume( DepthMap[this->GetName()][hit->GetCell()] )->GetCopyNo() );

  hit->SetCell( hist->GetVolume( detmap.depth )->GetCopyNo() );

  

  hit->SetRow( detmap.Row[hit->GetCell()] );
  hit->SetCol( detmap.Col[hit->GetCell()] );
  hit->SetPlane( detmap.Plane[hit->GetCell()] );
  hit->SetCellCoords( detmap.LocalCoord[hit->GetCell()] );
  hit->SetWire( detmap.Wire[hit->GetCell()] );

  // G4cout << "During CAL hit processing, SDname = " << SensitiveDetectorName << " physical volume name = " << hist->GetVolume( detmap.depth )->GetName() 
  // 	 << " copy number = " << hit->GetCell() << " (row,col)=(" << hit->GetRow() << ", " << hit->GetCol() << ")" << G4endl;

  //hit->SetGlobalCellCoords( detmap.GlobalCoord[hit->GetCell()] );
  G4AffineTransform aTrans_inverse = aTrans.Inverse();
  hit->SetGlobalCellCoords( aTrans_inverse.TransformPoint( G4ThreeVector(0,0,0) ) );

  /*
  printf("%f %f %f %f - dist %f beta %e momentum %f GeV\n", 
	  aStep->GetPreStepPoint()->GetLocalTime()/ns,
	  aStep->GetPostStepPoint()->GetLocalTime()/ns,
	  aStep->GetPreStepPoint()->GetGlobalTime()/ns,
	  aStep->GetPostStepPoint()->GetGlobalTime()/ns,
	  aStep->GetPreStepPoint()->GetPosition().mag()/m,
	  aStep->GetPreStepPoint()->GetPosition().mag()/aStep->GetPreStepPoint()->GetLocalTime()/(0.299792458*m/ns),
//aStep->GetTrack()->GetTotalEnergy()/GeV
	  edep/GeV
	  );
	  */

  hitCollection->insert( hit );

  return true;
}

void G4SBSCalSD::EndOfEvent(G4HCofThisEvent*HCE)
{
    G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    HCE->AddHitsCollection( HCID, hitCollection );

    //G4cout << "Adding hit collection " << collectionName[0] << " to HCE, HCID = " << HCID << G4endl;
}

void G4SBSCalSD::clear()
{
} 

void G4SBSCalSD::DrawAll()
{
} 

void G4SBSCalSD::PrintAll()
{
} 
