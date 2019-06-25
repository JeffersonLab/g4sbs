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
  detmap.SDname = name;
  detmap.clear();

  SDtracks.Clear();
  SDtracks.SetSDname(name);
}

G4SBSRICHSD::~G4SBSRICHSD(){;}

void G4SBSRICHSD::Initialize( G4HCofThisEvent *HC ){
  G4int HCID = -1;
  hitCollection = new G4SBSRICHHitsCollection( fullPathName.strip(G4String::leading,'/'), collectionName[0] );
  
  if( HCID < 0 ){
    HCID = GetCollectionID(0);
  }

  //  G4cout << "Adding hit collection " << collectionName[0] << " to HCE, HCID = " << HCID << ", SDname = " << SensitiveDetectorName << G4endl;

  HC->AddHitsCollection( HCID, hitCollection );

  SDtracks.Clear();
}

G4bool G4SBSRICHSD::ProcessHits( G4Step *aStep, G4TouchableHistory* ){
  //G4cout << "Processing RICH hits, SDname = " << SensitiveDetectorName << G4endl;

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

  G4TouchableHistory* hist = (G4TouchableHistory*) (prestep->GetTouchable());
  
  //Where to get all the information needed for the hit? 
  //From get methods of G4Step, G4Track, and G4StepPoint

  newHit->SetTrackID( track->GetTrackID() );
  newHit->SetTrackPID( track->GetParticleDefinition()->GetPDGEncoding() );
  newHit->SetMotherID( track->GetParentID() ); 
  //Getting the rest of the Mother particle information must wait until the end of the event we can access G4Event->G4TrajectoryContainer

  newHit->SetVertex( track->GetVertexPosition() );
  newHit->SetVertexDirection( track->GetVertexMomentumDirection() );

  //Record where was the optical photon produced:
  //Three cases are interesting: Aerogel tiles, RICH box containing C4F10 gas, or UVT-lucite of aerogel exit window (also PMT windows, PMT cathode volume, or PMT Quartz window):
  G4String namevol_origin = track->GetLogicalVolumeAtVertex()->GetName();
  int origin_flag = 0; //default
  if( namevol_origin.contains("Aerogel_tile_log") ){  //Aerogel
    origin_flag = 1;
  } else if (namevol_origin.contains("SBS_RICH_log") || namevol_origin.contains("GC_Tank_log") ){ //C4F10
    origin_flag = 2;
  } else if (namevol_origin.contains("Aero_exitwindow") ){ //UVT lucite
    origin_flag = 3;
  } else if (namevol_origin.contains("PMTwindow_log") || namevol_origin.contains("PMTcathode_log") || namevol_origin.contains("PMTQuartzWindow_log") || namevol_origin.contains("GC_PMT_Glass_log")  ){
    origin_flag = 4;
  }

  newHit->SetOriginVol( origin_flag );

  newHit->Setdx( aStep->GetStepLength() );
  newHit->SetTime( prestep->GetGlobalTime() );
  newHit->SetEdep( edep );
  newHit->SetEnergy( prestep->GetTotalEnergy() );

  //Let's change this so that it refers to the local position and direction of the hit:
  
  G4AffineTransform RICH_atrans = ( (G4TouchableHistory*) prestep->GetTouchable() )->GetHistory()->GetTopTransform();

  G4ThreeVector pos = prestep->GetPosition(); //Global position
  
  newHit->SetPos( pos ); //Global

  pos = RICH_atrans.TransformPoint(pos); 

  G4ThreeVector mom = prestep->GetMomentumDirection();

  newHit->SetDirection( mom ); //Global

  mom *= (RICH_atrans.NetRotation()).inverse();

  //local position and direction:
  newHit->SetLPos( pos );
  newHit->SetLDirection( mom );

  //Get PMT identifying information:

  //int PMTno;

  //G4cout << "RICH hit physical volume name, pre-step = " << prestep->GetPhysicalVolume()->GetName() << G4endl; 
  
  newHit->SetPMTnumber( hist->GetVolume( detmap.depth )->GetCopyNo() );

  
  //  newHit->Setrownumber( newHit->calc_row( PMTno ) );
  //  newHit->Setcolnumber( newHit->calc_col( PMTno ) );

  newHit->Setrownumber( detmap.Row[newHit->GetPMTnumber()] );
  newHit->Setcolnumber( detmap.Col[newHit->GetPMTnumber()] );

  // G4cout << "Hit RICH PMT number " << newHit->GetPMTnumber() << G4endl;
  // G4cout << "Physical volume name = " << prestep->GetPhysicalVolume()->GetName() << ", copy = " << newHit->GetPMTnumber()
  // 	 << " (row,col)=(" << newHit->Getrownumber() << ", " << newHit->Getcolnumber() << ")" << G4endl;

  
  newHit->SetCellCoord( detmap.LocalCoord[newHit->GetPMTnumber()] );
  
  G4AffineTransform RICH_atrans_inverse  = RICH_atrans.Inverse();
  
  newHit->SetGlobalCellCoord( RICH_atrans_inverse.TransformPoint( G4ThreeVector(0,0,0) ) );

  //newHit->SetLogicalVolume( prestep->GetPhysicalVolume()->GetLogicalVolume() );

  //Set quantum efficiency: default to 1 in case no material properties table has been defined:
  newHit->SetQuantumEfficiency( 1.0 );

  G4MaterialPropertiesTable *MPT = prestep->GetPhysicalVolume()->GetLogicalVolume()->GetMaterial()->GetMaterialPropertiesTable();

  if( MPT != NULL ){
    G4MaterialPropertyVector *QEvect = (G4MaterialPropertyVector*) MPT->GetProperty("EFFICIENCY");

    if( QEvect != NULL ){
      G4double Ephoton = newHit->GetEnergy();
      
      G4bool inrange = Ephoton >= QEvect->GetMinLowEdgeEnergy() && Ephoton <= QEvect->GetMaxLowEdgeEnergy();
      
      G4double QEtemp = QEvect->Value( Ephoton );

      if( !inrange ) QEtemp = 0.0;

      newHit->SetQuantumEfficiency(QEtemp);
    }
  }

  //G4cout << "G4SDname = " << GetName() << ", Ephoton = " << newHit->GetEnergy()/CLHEP::eV << ", QE = " << newHit->GetQuantumEfficiency() << G4endl;

  newHit->SetOTrIdx( SDtracks.InsertOriginalTrackInformation( track ) );
  newHit->SetPTrIdx( SDtracks.InsertPrimaryTrackInformation( track ) );
  newHit->SetSDTrIdx( SDtracks.InsertSDTrackInformation( track ) );
  
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

