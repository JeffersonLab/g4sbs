#define  _L_UNIT m
#define  _E_UNIT GeV
#define  _T_UNIT ns

#include "TMatrix.h"
#include "TVector.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "THashTable.h"

#include "G4SBSEventAction.hh"

#include "G4SBSEventGen.hh"
#include "G4SBSCalHit.hh"
#include "G4SBSGEMHit.hh"
#include "G4SBSGEMSD.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"


#include "G4SBSIO.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4MaterialPropertiesTable.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4TrajectoryContainer.hh"

#include <map>
//#include <unordered_map>
#include <vector>
#include <set>

using namespace std;

#define MAXHIT 2000

G4SBSEventAction::G4SBSEventAction()
{
    fGEMres = 70.0*um;

    // Load up resolution file if it exists

    int idx;

    for( idx = 0; idx < __MAXGEM; idx++ ){
	fGEMsigma[idx] = 1.0;
    }

    SDlist.clear();
}

G4SBSEventAction::~G4SBSEventAction()
{;}


void G4SBSEventAction::LoadSigmas(const char filename[] ){
    printf("Loading sigmas file %s\n", filename);
    FILE *sigmafile = fopen(filename, "r");

    int idx;
    if( !sigmafile ){ 
	printf("WARNING:  file %s not found.  Initializing to 1.0\n", filename);

	for( idx = 0; idx < __MAXGEM; idx++ ){
	    fGEMsigma[idx] = 1.0;
	}
	return;
    }

    int nread = 2;
    int gemidx;

    idx = 0;
    while( nread == 2 ){
	nread = fscanf(sigmafile, "%d%lf", &gemidx, &fGEMsigma[idx]);
	if( nread==2 && idx+1 == gemidx){ idx++; }
    }

    return;
}


void G4SBSEventAction::BeginOfEventAction(const G4Event*ev) {
   if( (ev->GetEventID()%1000)==0 ){
	printf("Event %8d\r", ev->GetEventID());
	fflush(stdout);
    }

    return;
}

void G4SBSEventAction::EndOfEventAction(const G4Event* evt )
{
  G4SDManager * SDman = G4SDManager::GetSDMpointer();

  //Let's add some more sophisticated checks so we don't print warnings on every event for sensitive detectors that don't exist:
  G4String GEMSDname = "G4SBS/GEM";
  G4String BBCalSDname = "G4SBS/BBCal";
  G4String HCalSDname = "G4SBS/HCAL";
  G4String RICHSDname = "G4SBS/RICH";

  bool GEMSD_exists = false;
  bool BBCalSD_exists = false;
  bool HCALSD_exists = false;
  bool RICHSD_exists = false; 
  
  G4bool warn = false;

  G4SBSGEMSD *GEMSDptr;

  //Now set flags for existence of sensitive detectors ***WITHOUT*** warnings:
  // if( (GEMSDptr = ( (G4SBSGEMSD*) SDman->FindSensitiveDetector( GEMSDname, warn ) ) ) )GEMSD_exists = true;
  // if( SDman->FindSensitiveDetector( BBCalSDname, warn ) ) BBCalSD_exists = true;
  // if( SDman->FindSensitiveDetector( HCalSDname, warn ) ) HCALSD_exists = true;
  // G4cout << "checking for RICH SD" << G4endl;
  // if( SDman->FindSensitiveDetector( RICHSDname, warn ) ) RICHSD_exists = true;
  // G4cout << "RICH detector exists = " << RICHSD_exists << G4endl;

  if( SDlist.find( GEMSDname ) != SDlist.end() ) GEMSD_exists = true;
  if( SDlist.find( BBCalSDname ) != SDlist.end() ) BBCalSD_exists = true;
  if( SDlist.find( HCalSDname ) != SDlist.end() ) HCALSD_exists = true;
  if( SDlist.find( RICHSDname ) != SDlist.end() ) RICHSD_exists = true;

  if( GEMSD_exists ) GEMSDptr = ( (G4SBSGEMSD*) SDman->FindSensitiveDetector( GEMSDname, false ) );

  bool anySD_exists = GEMSD_exists || BBCalSD_exists || HCALSD_exists || RICHSD_exists;

  if( !anySD_exists ) return;

  G4String colNam;

  //G4cout << "RICHCollID = " << RICHCollID << endl;
  //G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  
  MapTracks(evt);

  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  G4SBSCalHitsCollection* bbcalHC    = 0;
  G4SBSCalHitsCollection* hcalHC = 0;
  G4SBSGEMHitsCollection* gemHC = 0;
  G4SBSRICHHitsCollection *RICHHC = 0;
 
  bool hasbb   = false;
  bool hashcal = false;

  int i;

  cal_t caldata;

  caldata.bcndata = 0;
  caldata.hcndata = 0;

  if( BBCalSD_exists ) {
    
    bbcalCollID = SDman->GetCollectionID(colNam="BBCalcol");
    if( HCE && bbcalCollID >= 0 ){

      bbcalHC = (G4SBSCalHitsCollection*)(HCE->GetHC(bbcalCollID));
      
      if( bbcalHC->entries() > 0 ){
	hasbb = true;
	
	double xsum = 0.0;
	double ysum = 0.0;
	double esum = 0.0;
	
	caldata.bcndata = bbcalHC->entries();
	
	if( bbcalHC->entries() > MAXHITDATA ){
	  G4cerr << "WARNING:  Number of hits exceeds array length - truncating" << G4endl;
	  caldata.bcndata = MAXHITDATA;
	}
	
	
	for( i = 0; i < caldata.bcndata; i++ ){
	  xsum += (*bbcalHC)[i]->GetPos().x()*(*bbcalHC)[i]->GetEdep();
	  ysum += (*bbcalHC)[i]->GetPos().y()*(*bbcalHC)[i]->GetEdep();
	  esum += (*bbcalHC)[i]->GetEdep();
	  
	  caldata.bcx[i] = (*bbcalHC)[i]->GetPos().x()/_L_UNIT;
	  caldata.bcy[i] = (*bbcalHC)[i]->GetPos().y()/_L_UNIT;
	  caldata.bcz[i] = (*bbcalHC)[i]->GetPos().z()/_L_UNIT;
	  caldata.bce[i] = (*bbcalHC)[i]->GetEdep()/_E_UNIT;
	  caldata.bct[i] = (*bbcalHC)[i]->GetTime()/_T_UNIT;
	  
	  caldata.bcvx[i] = (*bbcalHC)[i]->GetVertex().x()/_L_UNIT;
	  caldata.bcvy[i] = (*bbcalHC)[i]->GetVertex().y()/_L_UNIT;
	  caldata.bcvz[i] = (*bbcalHC)[i]->GetVertex().z()/_L_UNIT;
	  
	  caldata.bcpid[i] = (*bbcalHC)[i]->GetPID();
	  caldata.bcmid[i] = (*bbcalHC)[i]->GetMID();
	  caldata.bctrid[i] = (*bbcalHC)[i]->GetTrID();
	}
      }
    }
  }

  if( HCALSD_exists ){

    hcalCollID  = SDman->GetCollectionID(colNam="HCALcol");
    
    if( HCE && hcalCollID >= 0 ){
      hcalHC  = (G4SBSCalHitsCollection*)(HCE->GetHC(hcalCollID));

      //if(hcalHC) {
      if( hcalHC->entries() > 0 ){
	hashcal = true;
	
	double xsum = 0.0;
	double ysum = 0.0;
	
	double xlsum = 0.0;
	double ylsum = 0.0;
	double zlsum = 0.0;
	
	double esum = 0.0;
	
	caldata.hcndata = hcalHC->entries();
	
	if( hcalHC->entries() > MAXHITDATA ){
	  G4cerr << "WARNING:  Number of hits exceeds array length - truncating" << G4endl;
	  caldata.hcndata = MAXHITDATA;
	}
	
	for( i = 0; i < caldata.hcndata; i++ ){
	  xsum += (*hcalHC)[i]->GetPos().x()*(*hcalHC)[i]->GetEdep();
	  ysum += (*hcalHC)[i]->GetPos().y()*(*hcalHC)[i]->GetEdep();
	  xlsum += (*hcalHC)[i]->GetLabPos().x()*(*hcalHC)[i]->GetEdep();
	  ylsum += (*hcalHC)[i]->GetLabPos().y()*(*hcalHC)[i]->GetEdep();
	  zlsum += (*hcalHC)[i]->GetLabPos().y()*(*hcalHC)[i]->GetEdep();
	  esum += (*hcalHC)[i]->GetEdep();
	  
	  // if( (*hcalHC)[i]->GetMID() == 0 ){
	  //   //trdata.hct = (*hcalHC)[i]->GetTime()/ns + CLHEP::RandGauss::shoot(0.0, fevgen->GetToFres());
	  // }
	  
	  caldata.hcx[i] = (*hcalHC)[i]->GetPos().x()/_L_UNIT;
	  caldata.hcy[i] = (*hcalHC)[i]->GetPos().y()/_L_UNIT;
	  caldata.hcz[i] = (*hcalHC)[i]->GetPos().z()/_L_UNIT;
	  caldata.hce[i] = (*hcalHC)[i]->GetEdep()/_E_UNIT;
	  caldata.hct[i] = (*hcalHC)[i]->GetTime()/_T_UNIT;
	  
	  caldata.hcvx[i] = (*hcalHC)[i]->GetVertex().x()/_L_UNIT;
	  caldata.hcvy[i] = (*hcalHC)[i]->GetVertex().y()/_L_UNIT;
	  caldata.hcvz[i] = (*hcalHC)[i]->GetVertex().z()/_L_UNIT;
	  
	  caldata.hcpid[i] = (*hcalHC)[i]->GetPID();
	  caldata.hcmid[i] = (*hcalHC)[i]->GetMID();
	  caldata.hctrid[i] = (*hcalHC)[i]->GetTrID();
	}
      
	
	// G4ThreeVector avglab = G4ThreeVector( trdata.hclx, trdata.hcly, trdata.hclz);
	
	
	// // Calculate expected time of flight
	// G4ThreeVector q3m = fevgen->GetBeamP()-fevgen->GetElectronP();
	// G4ThreeVector path = avglab-fevgen->GetV();
	// double hcd = path.mag();
	
	// trdata.hctex = hcd/(q3m.mag()*(0.3*m/ns)/sqrt(q3m.mag()*q3m.mag()+proton_mass_c2*proton_mass_c2))/_T_UNIT;
	// // Angular difference between q and reconstructed vector
	// double cosang = q3m.unit()*path.unit();
	// if( cosang > 1.0 ){ cosang = 1.0; } //  Apparent numerical problems in this dot product
	// trdata.hcdang = acos(cosang);
      }
    }
  }

  // If we don't have something in both arms end
  // and don't fill
  /*
  if( !hasbb && !hashcal ){
      return;
  }
  */


  // trdata.hcal = hashcal;
  // trdata.bb   = hasbb;

  // trdata.x = trdata.y = trdata.xp = trdata.yp = -1e9;
  // trdata.gemtr = 0;
  int idx, nhit, gid;

  hit_t hitdata;

  hitdata.ndata = 0;
  
  G4SBSTrackerOutput Toutput;

  int has_earm_track = 0;
  int has_harm_track = 0;

  if( GEMSD_exists ){

    gemCollID   = SDman->GetCollectionID(colNam="GEMcol");
    
    if( HCE && gemCollID >= 0 ){
      gemHC   = (G4SBSGEMHitsCollection*)(HCE->GetHC(gemCollID));
      // Need at least three hits to draw a line
      
      FillTrackData( evt, gemHC, Toutput ); 
     
      if( Toutput.ntracks > 0 ){
	for(int track=0; track<Toutput.ntracks; track++){
	  if( (GEMSDptr->GEMArmIDs)[Toutput.TrackerID[track]] == kEarm ) has_earm_track += 1;
	  if( (GEMSDptr->GEMArmIDs)[Toutput.TrackerID[track]] == kHarm ) has_harm_track += 1;
	}
      }

      nhit = 0;
      for( idx = 0; idx < gemHC->entries() && idx < MAXHIT; idx++ ){
	gid = (*gemHC)[idx]->GetGEMID();
	
	if( gid == 0 ) continue;
	
	// tx  =  -1e9;
	// ty  =  -1e9;
	// txp =  -1e9;
	// typ =  -1e9;
	
	// if( gid == 1 ){
	//   //  Project back to z = 0 plane
	//   tx  =  (*gemHC)[idx]->GetPos().getX() - (*gemHC)[idx]->GetXp()*(*gemHC)[idx]->GetPos().getZ();
	//   ty  =  (*gemHC)[idx]->GetPos().getY() - (*gemHC)[idx]->GetYp()*(*gemHC)[idx]->GetPos().getZ();
	//   txp =  (*gemHC)[idx]->GetXp();
	//   typ =  (*gemHC)[idx]->GetYp();
	// };
	
	//What is the purpose of the following line?
	//map |= (1 << (*gemHC)[idx]->GetGEMID()); 
	
	// Smear by resolution
	double lx = (*gemHC)[idx]->GetPos().getX() + CLHEP::RandGauss::shoot(0.0, fGEMres);
	double ly = (*gemHC)[idx]->GetPos().getY() + CLHEP::RandGauss::shoot(0.0, fGEMres);
	double lz = (*gemHC)[idx]->GetPos().getZ();
	
	hitdata.gid[nhit] = (*gemHC)[idx]->GetGEMID();
	hitdata.trid[nhit] = (*gemHC)[idx]->GetTrID();
	hitdata.trkrid[nhit] = (*gemHC)[idx]->GetTrackerID();
	hitdata.pid[nhit] = (*gemHC)[idx]->GetPID();
	hitdata.mid[nhit] = (*gemHC)[idx]->GetMID();
	hitdata.x[nhit] =  -ly/_L_UNIT; 
	hitdata.y[nhit] = lx/_L_UNIT;
	hitdata.z[nhit] = lz/_L_UNIT;
	hitdata.t[nhit] = (*gemHC)[idx]->GetHittime()/_T_UNIT;
	hitdata.p[nhit] = (*gemHC)[idx]->GetMom()/_E_UNIT;
	hitdata.edep[nhit] = (*gemHC)[idx]->GetEdep()/_E_UNIT;
	hitdata.pathl[nhit] = (*gemHC)[idx]->GetPathl()/_L_UNIT;
	
	hitdata.vx[nhit] =  (*gemHC)[idx]->GetVertex().getX()/_L_UNIT;
	hitdata.vy[nhit] =  (*gemHC)[idx]->GetVertex().getY()/_L_UNIT;
	hitdata.vz[nhit] =  (*gemHC)[idx]->GetVertex().getZ()/_L_UNIT;

	hitdata.vp[nhit] =  (*gemHC)[idx]->GetVertexMom().mag()/_E_UNIT;

	hitdata.vpx[nhit] =  (*gemHC)[idx]->GetVertexMom().getX()/_E_UNIT;
	hitdata.vpy[nhit] =  (*gemHC)[idx]->GetVertexMom().getY()/_E_UNIT;
	hitdata.vpz[nhit] =  (*gemHC)[idx]->GetVertexMom().getZ()/_E_UNIT;
	
	hitdata.tx[nhit] =  -(*gemHC)[idx]->GetPos().getY()/_L_UNIT;
	hitdata.ty[nhit] =   (*gemHC)[idx]->GetPos().getX()/_L_UNIT;
	hitdata.txp[nhit] =  -(*gemHC)[idx]->GetYp();
	hitdata.typ[nhit] = (*gemHC)[idx]->GetXp();
	
	//	  printf("GEM HIT %d (%f) %f %f\n", (*gemHC)[idx]->GetGEMID(), lz[nhit]/cm, lx[nhit]/cm, ly[nhit]/cm );
	nhit++;
	
	

	if( nhit == MAXHITDATA ){
	  G4cerr << "WARNING:  Number of hits exceeds array length - truncating" << G4endl;
	  break;
	}
      }
      hitdata.ndata = nhit;
    }
  }

  

  G4SBSRICHoutput richdata;

  richdata.timewindow = 10.0*ns; //group photons arriving at a PMT within timewindow into single hits
  richdata.threshold =  0.5; //Threshold on the number of photoelectrons

  if( RICHSD_exists ){ //RICH hit collection exists:
    RICHCollID  = SDman->GetCollectionID(colNam="RICHcoll");
    if( HCE && RICHCollID >= 0 ){
      RICHHC  = (G4SBSRICHHitsCollection*)(HCE->GetHC(RICHCollID));
      FillRICHData( evt, RICHHC, richdata );
    }
  }

  ev_t evdata = fIO->GetEventData();

  if( hasbb && has_earm_track != 0 ){
    evdata.earmaccept = 1;
  } else {
    evdata.earmaccept = 0;
  } 
  if( hashcal && has_harm_track != 0 ){
    evdata.harmaccept = 1;
  } else {
    evdata.harmaccept = 0;
  }

  fIO->SetEventData( evdata );
  //fIO->SetTrackData(trdata);
  fIO->SetCalData(caldata);
  fIO->SetHitData(hitdata);
  fIO->SetRICHData(richdata);
  fIO->SetTrackData( Toutput );

  fIO->FillTree();

  return;
}

void G4SBSEventAction::FillRICHData( const G4Event *evt, G4SBSRICHHitsCollection *hits, G4SBSRICHoutput &richoutput ){
  //Here is where we traverse the hit collection of the RICH and extract useful output data. 
  
  // set<int> TIDs; //list of all unique track IDs  
  //This is the total number of tracking steps in our sensitive volume!
  int nG4hits = hits->entries();

  //  cout << "Filling RICH data, nhits = " << nG4hits << endl;

  richoutput.Clear();

  set<int> TIDs_unique;  //set of all unique photon tracks involved in RICH hits.
  set<int> PMTs_unique;  //set of all unique PMTs with detected photons.
  set<int> mTIDs_unique;

  map<int,bool> Photon_used;
  map<int,bool> Photon_detected;
  map<int,double> Photon_energy;
  map<int,double> Photon_hittime;
  map<int,G4ThreeVector> Photon_position; // position at hit
  map<int,G4ThreeVector> Photon_direction; // direction at hit
  map<int,G4ThreeVector> Photon_vertex;  //emission vertex
  map<int,G4ThreeVector> Photon_vdirection;  //emission direction
  map<int,int> Photon_PMT;
  map<int,int> Photon_row;
  map<int,int> Photon_col;
  map<int,int> Photon_mTID;
  map<int,int> Photon_origvol;
  map<int,int> Photon_nsteps;

  map<int,int> PMT_Numphotoelectrons;
  map<int,int> PMT_row;
  map<int,int> PMT_col;
  map<int,double> PMT_hittime;
  map<int,double> PMT_hittime2; //hit time squared.
  map<int,double> PMT_rmstime;
  map<int,int> PMT_mTID;
  map<int,G4ThreeVector> PMT_pos;
  map<int,G4ThreeVector> PMT_dir;
  map<int,G4ThreeVector> PMT_vpos;
  map<int,G4ThreeVector> PMT_vdir;
  map<int,int> PMT_origvol;
  
  G4MaterialPropertiesTable *MPT;

  //G4SBSRICHoutput allhits; //"Unfiltered" data (later we will filter hits according to photoelectron threshold)
  //allhits.nhits_RICH = 0;
  //allhits.ntracks_RICH = 0;

  //int num_photons = 0;
  //Loop over all steps in the hits collection of the RICH:

  for( int step = 0; step < nG4hits; step++ ){
    //Retrieve all relevant information for this step:
    int pmt = (*hits)[step]->GetPMTnumber();
    int row = (*hits)[step]->Getrownumber();
    int col = (*hits)[step]->Getcolnumber();
    int tid = (*hits)[step]->GetTrackID();
    int mid = (*hits)[step]->GetMotherID();
    double Ephoton = (*hits)[step]->Getenergy();
    double Hittime = (*hits)[step]->GetTime();
    int origin_volume = (*hits)[step]->GetOriginVol();
    G4ThreeVector vertex = (*hits)[step]->GetVertex();
    G4ThreeVector vertexdirection = (*hits)[step]->GetVertexDirection();
    G4ThreeVector position = (*hits)[step]->GetPos();
    G4ThreeVector direction = (*hits)[step]->GetDirection();
    
    //First, we must ask: Is this a "new" photon track? It **should be** theoretically impossible for the same photon to 
    //be detected in two different PMTs, since the photocathodes don't have the necessary properties defined for the 
    //propagation of optical photons (RINDEX). On the other hand, it IS possible, and even likely, for two or more photons
    //to be detected by the same PMT. Therefore, we should consider photon tracks at the top level of the sorting logic:
      
    std::pair< set<int>::iterator, bool > testphotontrack = TIDs_unique.insert( tid );
    
    if( testphotontrack.second ){ //New photon track: determine whether this photon is detected:
      
      MPT = (*hits)[step]->GetLogicalVolume()->GetMaterial()->GetMaterialPropertiesTable();
      G4MaterialPropertyVector *QE = NULL;
      
      bool QE_defined = false;
      if( MPT != NULL ){
	QE = (G4MaterialPropertyVector*) MPT->GetProperty("EFFICIENCY");
	if( QE != NULL ) QE_defined = true;
      }
      
      bool photon_detected = true;
      bool inrange=false;
      //double qeff = QE->GetValue( Ephoton, inrange );
      if( QE_defined && G4UniformRand() > QE->GetValue( Ephoton, inrange ) ){ 
	//if quantum efficiency has been defined for the material in question, reject hit with probability 1 - QE:
	photon_detected = false;
      }
      
      Photon_used[ tid ] = !photon_detected; //If the photon is not detected, then we mark it as used. Otherwise, we mark it as unused, and it will be added to a PMT later.
      Photon_detected[ tid ] = photon_detected;
      Photon_energy[ tid ] = Ephoton;
      Photon_hittime[ tid ] = Hittime;
      Photon_position[ tid ] = position;
      Photon_direction[ tid ] = direction;
      Photon_vertex[ tid ] = vertex;
      Photon_vdirection[ tid ] = vertexdirection;
      Photon_PMT[ tid ] = pmt;
      Photon_row[ tid ] = row;
      Photon_col[ tid ] = col;
      Photon_mTID[ tid ] = mid; 
      Photon_origvol[ tid ] = origin_volume;
      Photon_nsteps[ tid ] = 1;
      
    } else { //existing photon, additional step. Increment averages of position, direction, time, etc for all steps of a detected photon. Don't bother for 
      //undetected photons...
      if( Photon_detected[ tid ] ){
	G4ThreeVector average_pos = (Photon_nsteps[ tid ] * Photon_position[ tid ] + position )/double( Photon_nsteps[tid] + 1 );
	Photon_position[tid] = average_pos;
	G4double average_time = (Photon_nsteps[ tid ] * Photon_hittime[ tid ] + Hittime )/double( Photon_nsteps[tid] + 1 );
	Photon_hittime[tid] = average_time; 
	Photon_nsteps[tid] += 1;
      }
    }
  }
  
  bool  remaining_hits = true;
  
  while( remaining_hits ) {
    
    remaining_hits = false;
    
    PMTs_unique.clear();

    for( set<int>::iterator it=TIDs_unique.begin(); it != TIDs_unique.end(); it++ ){
      int tid = *it;
      if( Photon_detected[ tid ] && !(Photon_used[ tid ] ) ){
	std::pair<set<int>::iterator,bool> testpmt = PMTs_unique.insert( Photon_PMT[tid] );

	int pmt = Photon_PMT[tid];

	if( testpmt.second ){ // new PMT;
	  
	  //Mark this photon track as used:
	  Photon_used[ tid ] = true;
	  
	  PMT_Numphotoelectrons[ pmt ] = 1;
	  PMT_row[ pmt ] = Photon_row[tid];
	  PMT_col[ pmt ] = Photon_col[tid];
	  PMT_hittime[ pmt ] = Photon_hittime[tid];
	  PMT_hittime2[ pmt ] = pow(Photon_hittime[tid],2);
	  PMT_mTID[ pmt ] = Photon_mTID[tid];
	  PMT_pos[ pmt ] = Photon_position[tid];
	  PMT_dir[ pmt ] = Photon_direction[tid];
	  PMT_vpos[ pmt ] = Photon_vertex[tid];
	  PMT_vdir[ pmt ] = Photon_vdirection[tid];
	  PMT_origvol[ pmt ] = Photon_origvol[tid];

	  mTIDs_unique.insert( Photon_mTID[tid] );

	} else if( fabs( Photon_hittime[tid] - PMT_hittime[pmt] ) <= richoutput.timewindow ){ //Existing pmt with multiple photon detections:
	  G4ThreeVector average_pos = (PMT_Numphotoelectrons[pmt] * PMT_pos[pmt] + Photon_position[tid])/double(PMT_Numphotoelectrons[pmt] + 1);
	  PMT_pos[pmt] = average_pos;
	  G4ThreeVector average_dir = (PMT_Numphotoelectrons[pmt] * PMT_dir[pmt] + Photon_direction[tid])/double(PMT_Numphotoelectrons[pmt] + 1);
	  PMT_dir[pmt] = average_dir;
	  G4ThreeVector average_vpos = (PMT_Numphotoelectrons[pmt] * PMT_vpos[pmt] + Photon_vertex[tid])/double(PMT_Numphotoelectrons[pmt] + 1);
	  PMT_vpos[pmt] = average_vpos;
	  G4ThreeVector average_vdir = (PMT_Numphotoelectrons[pmt] * PMT_vdir[pmt] + Photon_vdirection[tid])/double(PMT_Numphotoelectrons[pmt] + 1);
	  PMT_vdir[pmt] = average_vdir;
	  G4double average_hittime = (PMT_Numphotoelectrons[pmt] * PMT_hittime[pmt] + Photon_hittime[tid])/double(PMT_Numphotoelectrons[pmt] + 1 );
	  PMT_hittime[pmt] = average_hittime;
	  G4double average_hittime2 = (PMT_Numphotoelectrons[pmt] * PMT_hittime2[pmt] + pow(Photon_hittime[tid],2))/double(PMT_Numphotoelectrons[pmt] + 1 );
	  PMT_hittime2[pmt] = average_hittime2;

	  PMT_Numphotoelectrons[pmt] += 1;
	  
	  Photon_used[tid] = true;

	  mTIDs_unique.insert( Photon_mTID[tid] );
	}
	
	
	
	//If any photon is detected but not used, then remaining hits = true!
	if( !(Photon_used[tid] ) ) remaining_hits = true;
      }
    }
    
    //Now add hits to the output. PMTs are not required to be unique here, but the probability of multiple hits on the same PMT separated by more than richoutput.timewindow is very small. 
    for( set<int>::iterator it=PMTs_unique.begin(); it!=PMTs_unique.end(); it++ ){
      
      int pmt = *it;
      
      if( PMT_Numphotoelectrons[pmt] >= richoutput.threshold ){
	
	(richoutput.nhits_RICH)++;
	
	richoutput.PMTnumber.push_back( pmt );
	richoutput.row.push_back( PMT_row[pmt] );
	richoutput.col.push_back( PMT_col[pmt] );
	richoutput.NumPhotoelectrons.push_back( PMT_Numphotoelectrons[pmt] );
	richoutput.Time_avg.push_back( PMT_hittime[pmt] );
	richoutput.Time_rms.push_back( sqrt(fabs(PMT_hittime2[pmt] - pow(PMT_hittime[pmt],2) ) ) );
	richoutput.mTrackNo.push_back( PMT_mTID[pmt] ); 
	//For now, this is the mother track of the first photon detected by this PMT. This does not account for the possibility of the same PMT detecting photons produced by different tracks in the same event.
	//Hopefully, the probability of this occurrence is quite low?
	richoutput.xhit.push_back( PMT_pos[pmt].x()/_L_UNIT ); //Later we will go back to the hit definition and make sure this is the LOCAL position of the hit
	richoutput.yhit.push_back( PMT_pos[pmt].y()/_L_UNIT );
	richoutput.zhit.push_back( PMT_pos[pmt].z()/_L_UNIT );
	
	richoutput.pxhit.push_back( PMT_dir[pmt].x() ); //Later we will go back to the hit definition and make sure this is the LOCAL position of the hit
	richoutput.pyhit.push_back( PMT_dir[pmt].y() );
	richoutput.pzhit.push_back( PMT_dir[pmt].z() );
	
	richoutput.pvx.push_back( PMT_vpos[pmt].x()/_L_UNIT );
	richoutput.pvy.push_back( PMT_vpos[pmt].y()/_L_UNIT );
	richoutput.pvz.push_back( PMT_vpos[pmt].z()/_L_UNIT );
	
	richoutput.ppx.push_back( PMT_vdir[pmt].x() );
	richoutput.ppy.push_back( PMT_vdir[pmt].y() );
	richoutput.ppz.push_back( PMT_vdir[pmt].z() );
	
	richoutput.volume_flag.push_back( PMT_origvol[pmt] ); //Again, considering only the first photon detected by this PMT.
	
      }
    } 
    //PMTs_unique.clear();
  }
  
  G4TrajectoryContainer *tracklist = evt->GetTrajectoryContainer();

  if( !tracklist ) return;

  map<int,int> mtrackindex;

  for(set<int>::iterator it=mTIDs_unique.begin(); it!=mTIDs_unique.end(); it++ ){
    
    G4int mTrackID = *it;

    G4Trajectory *track = (G4Trajectory*) ( (*tracklist)[TrajectoryIndex[mTrackID]] );
    
    G4ThreeVector pinitial = track->GetInitialMomentum();
    
    richoutput.mPID.push_back( track->GetPDGEncoding() );
    richoutput.mTID.push_back( track->GetTrackID() );
    richoutput.mMID.push_back( track->GetParentID() );

    richoutput.mpx.push_back( pinitial.x()/_E_UNIT );
    richoutput.mpy.push_back( pinitial.y()/_E_UNIT);
    richoutput.mpz.push_back( pinitial.z()/_E_UNIT );
    
    G4ThreeVector vinitial = track->GetPoint(0)->GetPosition();
    
    richoutput.mvx.push_back( vinitial.x()/_L_UNIT );
    richoutput.mvy.push_back( vinitial.y()/_L_UNIT );
    richoutput.mvz.push_back( vinitial.z()/_L_UNIT );

    mtrackindex[mTrackID] = richoutput.ntracks_RICH;

    richoutput.ntracks_RICH++;
  }

  for(G4int ihit=0; ihit<richoutput.nhits_RICH; ihit++){
    set<int>::iterator pos = mTIDs_unique.find( richoutput.mTrackNo[ihit] );
    if( pos != mTIDs_unique.end() ){
      richoutput.mTrackNo[ihit] = mtrackindex[ *pos ];
    } else {
      richoutput.mTrackNo[ihit] = -1; //This should never happen!
    }
  }
}

void G4SBSEventAction::MapTracks( const G4Event *evt ){

  //This step is actually computationally expensive due to the large number of "trajectories" in any given event:

  TrajectoryIndex.clear();
  MotherTrackIDs.clear();

  G4TrajectoryContainer *tracklist = evt->GetTrajectoryContainer();
  
  if( tracklist ){

    for(unsigned int i=0; i<tracklist->size(); i++){
      G4Trajectory *track = (G4Trajectory*) (*tracklist)[i]; 
      // G4ThreeVector momentum = track->GetInitialMomentum();
      // G4ThreeVector pos = track->GetPoint(0)->GetPosition();
      
      // G4cout << "Trajectory " << i << ", ID = " << track->GetTrackID() <<", MID = " << track->GetParentID() << ", PID = " 
      // 	     << track->GetPDGEncoding() << ", (px,py,pz)=(" << momentum.x() << ", " << momentum.y() << ", " << momentum.z() 
      // 	     << "), (vx,vy,vz)=(" << pos.x() << ", " << pos.y() << ", " << pos.z() << ")" << G4endl; 
      G4int TrackID = track->GetTrackID();
      G4int MotherID = track->GetParentID();
      //TrackIDs.push_back( TrackID ); //This is the track ID number corresponding to trajectory number i
      TrajectoryIndex[ TrackID ] = i; 
      MotherTrackIDs[ TrackID ] = MotherID; //This is the mother track ID number corresponding to Track ID number 
      
    }
  }
}

void G4SBSEventAction::FillTrackData( const G4Event *evt, G4SBSGEMHitsCollection *HC, G4SBSTrackerOutput &Toutput){
  set<int> TrackerIDs_unique; //ID numbers of unique trackers with hits in this event:
  set<int> TrackTIDs_unique; //Track numbers of unique tracks causing GEM hits in this event:
  //set<int> GEMIDs; // ID numbers of unique GEM planes with hits in this event:

  map<int, map<int, vector<int> > > HitList;

  //Toutput.Clear();

  // map<int,int>     TrackerIDs; //mapping between hit indices and tracker IDs 
  // map<int,int>     GEMIDs;     //mapping between hit indices and plane IDs
  // map<int,int>     PIDs;       //mapping between hit indices and particle IDs
  // map<int,int>     MIDs;       //mapping between hit indices and Mother IDs ***Here we are only interested in hits with MID=0!***
  // map<int,int>     TIDs;

  // map<int,G4ThreeVector> HitPos;
  // map<int,G4ThreeVector> HitVertex;
  // map<int,double> HitXp;
  // map<int,double> HitYp;
  // map<int,double> HitP;
  // map<int,double> HitEdep;
  // map<int,double> HitTime;
  // map<int,double> HitBeta;

  G4int nhits = HC->entries();

  //First, map all the hits to Tracker IDs/GEM IDs:
  for(int i=0; i<nhits; i++){
    int tid =    (*HC)[i]->GetTrID();
    int trkrid = (*HC)[i]->GetTrackerID();
    //int pid =    (*HC)[i]->GetPID();
    int mid =    (*HC)[i]->GetMID();
    //int gid =    (*HC)[i]->GetGEMID();

    // G4ThreeVector pos = (*HC)[i]->GetPos();
    // G4ThreeVector vtx = (*HC)[i]->GetVertex();
    
    // double Xp = (*HC)[i]->GetXp();
    // double Yp = (*HC)[i]->GetYp();
    // double P  = (*HC)[i]->GetMom();
    double edep = (*HC)[i]->GetEdep();

    // double hittime = (*HC)[i]->GetHittime();
    // double beta    = (*HC)[i]->GetBeta();

    //Only consider hits with MID == 0 (primary particles):
    if( mid == 0 && edep > 0.0 ){
      TrackerIDs_unique.insert( trkrid );
      TrackTIDs_unique.insert( tid );
      
      HitList[trkrid][tid].push_back( i );

      // TrackerIDs[i] = trkrid;
      // GEMIDs[i]     = gid;
      // PIDs[i]       = pid;
      // MIDs[i]       = mid; 
      // TIDs[i]       = tid; 
      
      // HitPos[i]    = pos;
      // HitVertex[i] = vtx;
      // HitP[i]      = P;
      // HitXp[i]     = Xp;
      // HitYp[i]     = Yp;
      // HitEdep[i]   = edep;

      // HitTime[i]   = hittime;
      // HitBeta[i]   = beta;
    }
  }

  int nplanes_min = 3; //Minimum number of valid hits to define a track:

  //For the "true" track, we simply take the average of all points for x, y, xp, yp 
  //Fit procedure: minimize chi^2 defined as sum_i=1,N (xi - (x0 + xp*z))^2/sigmax^2 + (yi - (y0+yp*z))^2/sigmay^2:

  for(set<int>::iterator trkr=TrackerIDs_unique.begin(); trkr != TrackerIDs_unique.end(); trkr++ ){
    for(set<int>::iterator trk=TrackTIDs_unique.begin(); trk != TrackTIDs_unique.end(); trk++ ){
      int tracker = *trkr;
      int track   = *trk;
      //Define sums to keep track of: 
      int nhittrk = 0;
      int nplanetrk = 0;
      
      set<int> GEMIDs_unique; //List of unique GEM plane ids on this track:

      double x0avg = 0.0, y0avg = 0.0, xpavg = 0.0, ypavg = 0.0;
      double x0avg2 = 0.0, y0avg2 = 0.0, xpavg2 = 0.0, ypavg2 = 0.0;
      double tavg = 0.0, tavg2 = 0.0;
      double pavg = 0.0;

      TMatrixD M(4,4);
      TVectorD b(4);
      TVectorD btrue(4);

      //Initialize linear fitting matrices to zero:
      for(int i=0; i<4; i++){
	for(int j=0; j<4; j++){
	  M(i,i) = 0.0;
	}
	b(i) = 0.0;
	btrue(i) = 0.0;
      }
      
      nhittrk = HitList[tracker][track].size();

      if( nhittrk >= 3 ){
	
	vector<double> xsmear,ysmear,xtrue,ytrue,ztrue,hittime,beta;

	int PID = 0;

	for(int idx = 0; idx<nhittrk; idx++ ){ //First pass: determine number of unique GEM planes and compute average slopes of "true track". Also compute the sums for the linear fit:
	  int hit = HitList[tracker][track][idx];
	  
	  G4SBSGEMHit *Hit = (*HC)[hit];

	  if( idx == 0 ) PID = Hit->GetPID();

	  GEMIDs_unique.insert( Hit->GetGEMID() );

	  pavg += Hit->GetMom()/double(nhittrk);

	  xtrue.push_back( Hit->GetPos().x() );
	  ytrue.push_back( Hit->GetPos().y() );

	  xsmear.push_back( xtrue[idx] + CLHEP::RandGauss::shoot(0.0, fGEMres) );
	  ysmear.push_back( ytrue[idx] + CLHEP::RandGauss::shoot(0.0, fGEMres) );
	  
	  ztrue.push_back( Hit->GetPos().z() );
	 
	  hittime.push_back( Hit->GetHittime() );
	  beta.push_back( Hit->GetBeta() );
	  
	  //Indices of fit parameters in matrix are: 0 = x0, 1 = xp, 2 = y0, 3 = yp:
	  M(0,0) += pow(1.0/fGEMres,2); 
	  M(0,1) += pow(1.0/fGEMres,2) * ztrue[idx];
	  M(0,2) += 0.0; //No cross term between x0 and y0
	  M(0,3) += 0.0; //No cross term between x0 and yp
	  M(1,0) += pow(1.0/fGEMres,2) * ztrue[idx];
	  M(1,1) += pow(ztrue[idx]/fGEMres,2);
	  M(1,2) += 0.0;
	  M(1,3) += 0.0;
	  M(2,0) += 0.0;
	  M(2,1) += 0.0;
	  M(2,2) += pow(1.0/fGEMres,2);
	  M(2,3) += pow(1.0/fGEMres,2) * ztrue[idx];
	  M(3,0) += 0.0;
	  M(3,1) += 0.0;
	  M(3,2) += pow(1.0/fGEMres,2) * ztrue[idx];
	  M(3,3) += pow(ztrue[idx]/fGEMres,2);
	  
	  b(0) += xsmear[idx] / pow(fGEMres,2);
	  b(1) += ztrue[idx] * xsmear[idx] / pow(fGEMres,2);
	  b(2) += ysmear[idx] / pow(fGEMres,2);
	  b(3) += ztrue[idx] * ysmear[idx] / pow(fGEMres,2);
	  
	  btrue(0) += xtrue[idx] / pow(fGEMres,2);
	  btrue(1) += ztrue[idx] * xtrue[idx] / pow(fGEMres,2);
	  btrue(2) += ytrue[idx] / pow(fGEMres,2);
	  btrue(3) += ztrue[idx] * ytrue[idx] / pow(fGEMres,2);

	}

	TMatrixD Minv = M.Invert();
	TVectorD FitTrack = Minv * b;
	TVectorD TrueTrack = Minv * btrue;
	
	if( GEMIDs_unique.size() >= nplanes_min ){
	  Toutput.ntracks++;
	  Toutput.TrackerID.push_back(tracker);
	  Toutput.TrackTID.push_back(track);
	  Toutput.TrackPID.push_back( PID );
	  Toutput.NumHits.push_back( nhittrk );
	  Toutput.NumPlanes.push_back( GEMIDs_unique.size() );
	  
	  Toutput.TrackX.push_back( -TrueTrack(2)/_L_UNIT );
	  Toutput.TrackXp.push_back( -TrueTrack(3) );
	  Toutput.TrackY.push_back( TrueTrack(0)/_L_UNIT );
	  Toutput.TrackYp.push_back( TrueTrack(1) );

	  Toutput.TrackXfit.push_back( -FitTrack(2)/_L_UNIT );
	  Toutput.TrackXpfit.push_back( -FitTrack(3) );
	  Toutput.TrackYfit.push_back( FitTrack(0)/_L_UNIT );
	  Toutput.TrackYpfit.push_back( FitTrack(1) );

	  int ndf = 0;
	  double chi2 = 0.0, chi2_true = 0.0;
	  //Compute chi^2 and focal plane times projected to local origin:
	  for(int idx = 0; idx<nhittrk; idx++){

	    //tfp is hit time corrected for time of flight:
	    double tfp = hittime[idx] - ztrue[idx]*sqrt( 1.0 + pow(TrueTrack(1),2)+pow(TrueTrack(3),2) ) / (beta[idx]*c_light);
	    
	    chi2 += pow( (xsmear[idx] - (FitTrack(0) + FitTrack(1)*ztrue[idx] ) )/fGEMres, 2 );
	    chi2 += pow( (ysmear[idx] - (FitTrack(2) + FitTrack(3)*ztrue[idx] ) )/fGEMres, 2 );

	    chi2_true += pow( (xtrue[idx] - (TrueTrack(0) + TrueTrack(1)*ztrue[idx] ) )/fGEMres, 2 );
	    chi2_true += pow( (ytrue[idx] - (TrueTrack(2) + TrueTrack(3)*ztrue[idx] ) )/fGEMres, 2 );

	    tavg += tfp/double(nhittrk);

	    ndf += 2;
	  }

	  Toutput.Chi2fit.push_back( chi2 );
	  Toutput.Chi2true.push_back( chi2_true );
	  Toutput.NDF.push_back( ndf );
	  Toutput.TrackT.push_back( tavg/_T_UNIT );

	  Toutput.TrackP.push_back( pavg/_E_UNIT );
	}
      }
    }
  }
}
