#define  _L_UNIT CLHEP::m
#define  _E_UNIT CLHEP::GeV
#define  _T_UNIT CLHEP::ns

#include "TMatrix.h"
#include "TVector.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "THashTable.h"
#include "TH1F.h"

#include "G4SBSEventAction.hh" //RICHHit.hh, RICHoutput.hh, ECalHit.hh, ECaloutput.hh are here

#include "G4SBSEventGen.hh"
#include "G4SBSCalHit.hh"
#include "G4SBSGEMHit.hh"
#include "G4SBSGEMSD.hh"
#include "G4SBSCalSD.hh"
#include "G4SBSRICHSD.hh"
#include "G4SBSECalSD.hh"

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

G4SBSEventAction::G4SBSEventAction() : fEventStatusEvery(1000)
{
  // if recommended "in most cases" then fTreeFlag should be 1 *by default*, 
  // and one should ask *explicitely* to have it deactivated.
  // That would avoid to take 2TB (!!!!!!!!!) for a stupid beam background simulation !!!!
  fTreeFlag = 1;
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
   if( (ev->GetEventID()%fEventStatusEvery)==0 ){
     printf("Event %8d\r", ev->GetEventID());
     fflush(stdout);
   }

    return;
}

void G4SBSEventAction::EndOfEventAction(const G4Event* evt )
{

  G4SDManager * SDman = G4SDManager::GetSDMpointer();

  //SDman->ListTree();

  //Let's add some more sophisticated checks so we don't print warnings on every event for sensitive detectors that don't exist:
  // G4String GEMSDname = "G4SBS/GEM";
  // G4String BBCalSDname = "G4SBS/BBCal";
  // G4String HCalSDname = "G4SBS/HCAL";
  // G4String RICHSDname = "G4SBS/RICH";
  // G4String ECalSDname = "G4SBS/ECal";

  // bool GEMSD_exists = false;
  // bool BBCalSD_exists = false;
  // bool HCALSD_exists = false;
  // bool RICHSD_exists = false; 
  // bool ECalSD_exists = false;
  
  G4bool warn = false;

  G4SBSGEMSD *GEMSDptr;
  G4SBSCalSD *CalSDptr;
  G4SBSRICHSD *RICHSDptr;
  G4SBSECalSD *ECalSDptr;

  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  G4SBSCalHitsCollection* calHC = 0;
  G4SBSGEMHitsCollection* gemHC = 0;
  G4SBSRICHHitsCollection *RICHHC = 0;
  G4SBSECalHitsCollection *ECalHC = 0;

  G4SBSBeamDiffuserSD *BDSDptr; 
  G4SBSBDHitsCollection *bdHC = 0; 

  G4SBSIonChamberSD *ICSDptr; 
  G4SBSICHitsCollection *icHC = 0; 

  G4SBSTargetSD *genGCSDptr; 
  G4SBSTargetHitsCollection *gcHC = 0; 

  G4SBSTargetSD *genCUSDptr; 
  G4SBSTargetHitsCollection *cuHC = 0; 

  G4SBSTargetSD *genALSDptr; 
  G4SBSTargetHitsCollection *alHC = 0; 

  G4SBSTargetSD *gen3HESDptr; 
  G4SBSTargetHitsCollection *he3HC = 0; 

  MapTracks(evt);

  bool anyhits = false;
  bool has_earm_track=false;
  bool has_harm_track=false;
  bool has_earm_cal=false;
  bool has_harm_cal=false;

  G4SBSSDTrackOutput allsdtracks;

  allsdtracks.Clear();

  //G4cout << "End-of-event processing for event ID " << evt->GetEventID() << G4endl;
  
  //Loop over all sensitive detectors:
  for( set<G4String>::iterator d=SDlist.begin(); d!=SDlist.end(); ++d ){
    G4String colNam;

    G4SBS::SDet_t Det_type = SDtype[*d];
    // G4SBS::Arm_t Det_arm = SDarm[d->first];

    G4SBSGEMoutput gd;
    G4SBSTrackerOutput td;
    G4SBSCALoutput cd;
    G4SBSRICHoutput rd;
    G4SBSECaloutput ed;
    G4SBSSDTrackOutput sd;
    G4SBSSDTrackOutput *sdtemp;
    G4SBSBDoutput bd; 
    G4SBSICoutput ic; 
    G4SBSTargetoutput gc; 
    G4SBSTargetoutput cu; 
    G4SBSTargetoutput al; 
    G4SBSTargetoutput he3; 
    
    switch(Det_type){

    case G4SBS::kGEM:
      GEMSDptr = (G4SBSGEMSD*) SDman->FindSensitiveDetector( *d, false );

      if( GEMSDptr != NULL ){
	gemHC = (G4SBSGEMHitsCollection*) (HCE->GetHC(SDman->GetCollectionID(colNam=GEMSDptr->GetCollectionName(0))));
	
	if( gemHC != NULL ){

	  sd = GEMSDptr->SDtracks;

	  //This should only be called once, otherwise units will be wrong!
	  sd.ConvertToTreeUnits();

	  sdtemp = &sd;

	  map<G4String,G4bool>::iterator keep = fIO->GetKeepSDtracks().find( *d );

	  G4bool keepthis = keep != fIO->GetKeepSDtracks().end() && keep->second;
	  
	  if( fIO->GetKeepAllSDtracks() || keepthis ){
	    allsdtracks.Merge( sd );
	    sdtemp = &allsdtracks;
	  }
	      
	  fIO->SetSDtrackData( *d, sd );
	  
	  FillGEMData(evt, gemHC, gd, *sdtemp );
	  fIO->SetGEMData( *d, gd );
	  
	
	  anyhits = (anyhits || gd.nhits_GEM > 0);
	  
	  FillTrackData( gd, td );
	  
	  if( td.ntracks > 0 ){
	    if( (*d).contains("Earm") ) has_earm_track = true;
	    if( (*d).contains("Harm") ) has_harm_track = true;
	    
	  }
	  fIO->SetTrackData( *d, td );
	}
      }
      break;
    case G4SBS::kCAL:

      CalSDptr = (G4SBSCalSD*) SDman->FindSensitiveDetector( *d, false );

      fhistogram_index = fIO->histogram_index[*d];
      
      if( CalSDptr != NULL ){
	calHC = (G4SBSCalHitsCollection*) (HCE->GetHC(SDman->GetCollectionID(colNam=CalSDptr->GetCollectionName(0))));

	if( calHC != NULL ){
	  cd.timewindow = CalSDptr->GetTimeWindow();
	  cd.threshold =  CalSDptr->GetEnergyThreshold();
	  cd.ntimebins =  CalSDptr->GetNTimeBins();
	  cd.hold_tbins = CalSDptr->hold_tbins;

	  sd = CalSDptr->SDtracks;

	  //This should only be called once, otherwise units will be wrong!
	  sd.ConvertToTreeUnits();

	  sdtemp = &sd;

	  map<G4String,G4bool>::iterator keep = fIO->GetKeepSDtracks().find( *d );

	  G4bool keepthis = keep != fIO->GetKeepSDtracks().end() && keep->second;
	  
	  if( fIO->GetKeepAllSDtracks() || keepthis ){
	    allsdtracks.Merge( sd );
	    sdtemp = &allsdtracks;
	  }
	  
	  //allsdtracks.Merge( sd ); //This has to be called before FillCalData or the output won't make sense
	  
	  fIO->SetSDtrackData( *d, sd );

	  // G4cout << "SD name = " << *d << G4endl;
	  // G4cout << "Hits collection SD name = " << calHC->GetSDname() << G4endl << G4endl;
	  
	  FillCalData( evt, calHC, cd, *sdtemp );
	  
	  fIO->SetCalData( *d, cd );
	  
	  anyhits = (anyhits || cd.nhits_CAL > 0);
	  
	  if( cd.nhits_CAL > 0 ){
	    if( (*d).contains("Earm") ) has_earm_cal = true;
	    if( (*d).contains("Harm") ) has_harm_cal = true;
	  }
	}
      }
      break;
    case G4SBS::kRICH:
      
      RICHSDptr = (G4SBSRICHSD*) SDman->FindSensitiveDetector( *d, false );

      if( RICHSDptr != NULL ){
	RICHHC = (G4SBSRICHHitsCollection*) (HCE->GetHC(SDman->GetCollectionID(colNam=RICHSDptr->GetCollectionName(0))));

	sd = RICHSDptr->SDtracks;

	//This should only be called once, otherwise units will be wrong!
	sd.ConvertToTreeUnits();

	sdtemp = &sd;

	map<G4String,G4bool>::iterator keep = fIO->GetKeepSDtracks().find( *d );

	G4bool keepthis = keep != fIO->GetKeepSDtracks().end() && keep->second;
	  
	if( fIO->GetKeepAllSDtracks() || keepthis ){
	  allsdtracks.Merge( sd );
	  sdtemp = &allsdtracks;
	}
	
	//allsdtracks.Merge( sd );
	
	fIO->SetSDtrackData( *d, sd );
	
	FillRICHData( evt, RICHHC, rd, *sdtemp );
	
	fIO->SetRICHData( *d, rd );
	
	anyhits = (anyhits || rd.nhits_RICH > 0);
      }
	 

      break;
    case G4SBS::kECAL:

      ECalSDptr = (G4SBSECalSD*) SDman->FindSensitiveDetector( *d, false );
     

      if( ECalSDptr != NULL ){
	ECalHC = (G4SBSECalHitsCollection*) (HCE->GetHC(SDman->GetCollectionID(colNam=ECalSDptr->GetCollectionName(0))));

	// *****
	ed.timewindow = ECalSDptr->GetTimeWindow();
	ed.threshold =  ECalSDptr->GetPEThreshold();
	ed.ntimebins =  ECalSDptr->GetNTimeBins();
	// *****
	
	sd = ECalSDptr->SDtracks;

	//This should only be called once, otherwise units will be wrong!
	sd.ConvertToTreeUnits();

	sdtemp = &sd;

	map<G4String,G4bool>::iterator keep = fIO->GetKeepSDtracks().find( *d );

	G4bool keepthis = keep != fIO->GetKeepSDtracks().end() && keep->second;
	
	if( fIO->GetKeepAllSDtracks() || keepthis ){
	  allsdtracks.Merge( sd );
	  sdtemp = &allsdtracks;
	}
	
	//allsdtracks.Merge( sd );
	
	fIO->SetSDtrackData( *d, sd );
	
	FillECalData( ECalHC, ed, *sdtemp );
	
	fIO->SetECalData( *d, ed );
	
	anyhits = (anyhits || ed.nhits_ECal > 0);
      }
      break;
    case G4SBS::kBD:  
      // beam diffuser (BD)
      BDSDptr = (G4SBSBeamDiffuserSD*) SDman->FindSensitiveDetector(*d,false);
      if(BDSDptr!=NULL){
	 bdHC = (G4SBSBDHitsCollection*) (HCE->GetHC(SDman->GetCollectionID(colNam=BDSDptr->GetCollectionName(0))));
	 if(bdHC!=NULL){
	    FillBDData(evt,bdHC,bd); 
	    fIO->SetBDData(*d,bd); 
	    anyhits = (anyhits || bd.nhits_BD>0 );
	 }
      } 
      break;
    case G4SBS::kIC:  
      // ion chamber (IC)  
      ICSDptr = (G4SBSIonChamberSD*) SDman->FindSensitiveDetector(*d,false);
      if(ICSDptr!=NULL){
	 icHC = (G4SBSICHitsCollection*) (HCE->GetHC(SDman->GetCollectionID(colNam=ICSDptr->GetCollectionName(0))));
	 if(icHC!=NULL){
	    FillICData(evt,icHC,ic); 
	    fIO->SetICData(*d,ic); 
	    anyhits = (anyhits || ic.nhits_IC>0 );
	 }
      } 
      break;
    case G4SBS::kTarget_GEn_Glass:  
      // GEn target glass cell  
      genGCSDptr = (G4SBSTargetSD*) SDman->FindSensitiveDetector(*d,false);
      if(genGCSDptr!=NULL){
	 gcHC = (G4SBSTargetHitsCollection*) (HCE->GetHC(SDman->GetCollectionID(colNam=genGCSDptr->GetCollectionName(0))));
	 if(gcHC!=NULL){
	    FillGEnTargetData(evt,gcHC,gc); 
	    fIO->SetGEnTargetData_Glass(*d,gc); 
	    anyhits = (anyhits || gc.nhits_Target>0 );
	 }
      } 
      break;
    case G4SBS::kTarget_GEn_Cu:  
      // GEn target Cu  
      genCUSDptr = (G4SBSTargetSD*) SDman->FindSensitiveDetector(*d,false);
      if(genCUSDptr!=NULL){
	 cuHC = (G4SBSTargetHitsCollection*) (HCE->GetHC(SDman->GetCollectionID(colNam=genCUSDptr->GetCollectionName(0))));
	 if(cuHC!=NULL){
	    FillGEnTargetData(evt,cuHC,cu); 
	    fIO->SetGEnTargetData_Cu(*d,cu); 
	    anyhits = (anyhits || cu.nhits_Target>0 );
	 }
      } 
      break;
    case G4SBS::kTarget_GEn_Al:  
      // GEn target Al
      genALSDptr = (G4SBSTargetSD*) SDman->FindSensitiveDetector(*d,false);
      if(genALSDptr!=NULL){
	 alHC = (G4SBSTargetHitsCollection*) (HCE->GetHC(SDman->GetCollectionID(colNam=genALSDptr->GetCollectionName(0))));
	 if(alHC!=NULL){
	    FillGEnTargetData(evt,alHC,al); 
	    fIO->SetGEnTargetData_Al(*d,al); 
	    anyhits = (anyhits || al.nhits_Target>0 );
	 }
      } 
      break;
    case G4SBS::kTarget_GEn_3He:  
      // GEn target3He 
      gen3HESDptr = (G4SBSTargetSD*) SDman->FindSensitiveDetector(*d,false);
      if(gen3HESDptr!=NULL){
	 he3HC = (G4SBSTargetHitsCollection*) (HCE->GetHC(SDman->GetCollectionID(colNam=gen3HESDptr->GetCollectionName(0))));
	 if(he3HC!=NULL){
	    FillGEnTargetData(evt,he3HC,he3); 
	    fIO->SetGEnTargetData_3He(*d,he3); 
	    anyhits = (anyhits || he3.nhits_Target>0 );
	 }
      } 
      break;
    }
  }

  fIO->SetAllSDtrackData( allsdtracks );

  //This copy operation may be inefficient:
  ev_t evdata = fIO->GetEventData();
  evdata.earmaccept = 0;
  evdata.harmaccept = 0;

  if( has_earm_track && has_earm_cal )
    evdata.earmaccept = 1;
  if( has_harm_track && has_harm_cal )
    evdata.harmaccept = 1;


  fIO->SetEventData( evdata );

  if( fTreeFlag == 0 || anyhits ) fIO->FillTree();

  return;
}

void G4SBSEventAction::FillGEMData( const G4Event *evt, G4SBSGEMHitsCollection *hits, G4SBSGEMoutput &gemoutput, G4SBSSDTrackOutput &sdtracks ){
  gemoutput.Clear();
  gemoutput.timewindow = 1000.0*ns;
  gemoutput.threshold = 0.0*eV;

  //Need to map unique layers and then unique tracks within layers:
  map<int,set<int> > tracks_layers; //key is GEM layer ID, value is a set of unique tracks with energy deposition in the layer
  //For each unique particle track in each GEM layer, we want to tabulate sums/averages of coordinates, etc:
  map<int,map<int,int> > nsteps_track_layer; //number of steps by track/layer
  map<int,map<int,double> > x,y,z,t,t2,tmin,tmax;
  map<int,map<int,double> > tx,ty;
  map<int,map<int,double> > xin,yin,zin;
  map<int,map<int,double> > xout,yout,zout;
  map<int,map<int,double> > txp,typ;
  map<int,map<int,double> > xg,yg,zg;
  map<int,map<int,int> > mid,pid; //don't need one for plane, trid as these are already keys
  map<int,map<int,double> > vx,vy,vz;
  map<int,map<int,double> > p,edep,beta,pmin;//pmin: lowest p reached by the track. Does not make its way to output
  map<int,map<int,double> > polx,poly,polz;
  map<int,map<int,int> > OTrackIndices;
  map<int,map<int,int> > PTrackIndices;
  map<int,map<int,int> > SDTrackIndices;

  //G4String sdname = hits->GetSDname();
  //G4int nhit=0;
  //Loop over all "hits" (actually individual tracking steps):
  for( G4int i=0; i < hits->entries(); i++ ){

    int trid = (*hits)[i]->GetTrID();
    int gemID = (*hits)[i]->GetGEMID();
    
    std::pair<set<int>::iterator, bool> track = tracks_layers[gemID].insert( trid );

    if( track.second ){ //new track in this layer, first step:
      nsteps_track_layer[gemID][trid] = 1;
      //The coordinates below represent local GEM hit coordinates in the TRANSPORT system 
      //(+x = vertical down, +y = horizontal left when looking in the direction of particle motion)
      x[gemID][trid] =  (*hits)[i]->GetPos().x();
      y[gemID][trid] =  (*hits)[i]->GetPos().y();
      z[gemID][trid] =  (*hits)[i]->GetPos().z();
      
      polx[gemID][trid] = (*hits)[i]->GetPolarization().x();
      poly[gemID][trid] = (*hits)[i]->GetPolarization().y();
      polz[gemID][trid] = (*hits)[i]->GetPolarization().z();

      t[gemID][trid] =  (*hits)[i]->GetHittime();
      t2[gemID][trid] = pow( t[gemID][trid],2 );
      tmin[gemID][trid] = tmax[gemID][trid] = t[gemID][trid];
      
      //Track coordinates in TRANSPORT system:
      tx[gemID][trid] = x[gemID][trid];
      ty[gemID][trid] = y[gemID][trid];

      xin[gemID][trid] = (*hits)[i]->GetPos().x();
      yin[gemID][trid] = (*hits)[i]->GetPos().y();
      zin[gemID][trid] = (*hits)[i]->GetPos().z();

      xout[gemID][trid] = (*hits)[i]->GetOutPos().x();
      yout[gemID][trid] = (*hits)[i]->GetOutPos().y();
      zout[gemID][trid] = (*hits)[i]->GetOutPos().z();

      txp[gemID][trid] = (*hits)[i]->GetXp();
      typ[gemID][trid] = (*hits)[i]->GetYp();
      
      xg[gemID][trid] = (*hits)[i]->GetGlobalPos().x();
      yg[gemID][trid] = (*hits)[i]->GetGlobalPos().y();
      zg[gemID][trid] = (*hits)[i]->GetGlobalPos().z();
      
      mid[gemID][trid] = (*hits)[i]->GetMID();
      pid[gemID][trid] = (*hits)[i]->GetPID();
      
      //Vertex coordinates of track in GLOBAL coordinates:
      vx[gemID][trid] = (*hits)[i]->GetVertex().x();
      vy[gemID][trid] = (*hits)[i]->GetVertex().y();
      vz[gemID][trid] = (*hits)[i]->GetVertex().z();
      
      p[gemID][trid] = (*hits)[i]->GetMom();
      pmin[gemID][trid] = (*hits)[i]->GetMom();
      
      edep[gemID][trid] = (*hits)[i]->GetEdep();
   
      beta[gemID][trid] = (*hits)[i]->GetBeta();

      // OTrackIndices[gemID][trid] = (*hits)[i]->GetOTrIdx();
      // PTrackIndices[gemID][trid] = (*hits)[i]->GetPTrIdx();
      // SDTrackIndices[gemID][trid] = (*hits)[i]->GetSDTrIdx();

      //Changed SDtrackoutput class so that the hit "track indices" are now actually G4 track IDs
      //In principle there is no need to check for existence of the tracks in the list
      OTrackIndices[gemID][trid] = sdtracks.otracklist[(*hits)[i]->GetOTrIdx()];
      PTrackIndices[gemID][trid] = sdtracks.ptracklist[(*hits)[i]->GetPTrIdx()];
      SDTrackIndices[gemID][trid] = sdtracks.sdtracklist[(*hits)[i]->GetSDTrIdx()][hits->GetSDname()];
      
    } else { //existing track in this layer, additional step; increment sums and averages:
      int nstep = nsteps_track_layer[gemID][trid];
      double w = double(nstep)/(double(nstep+1) );

      //The coordinates below represent local GEM hit coordinates in the TRANSPORT system 
      //(+x = vertical down, +y = horizontal left when looking in the direction of particle motion)
      x[gemID][trid] = x[gemID][trid]*w +( (*hits)[i]->GetPos().x() )*(1.0-w);
      y[gemID][trid] = y[gemID][trid]*w +( (*hits)[i]->GetPos().y() )*(1.0-w);
      z[gemID][trid] = z[gemID][trid]*w +( (*hits)[i]->GetPos().z() )*(1.0-w);
      
      t[gemID][trid] = t[gemID][trid]*w +( (*hits)[i]->GetHittime() )*(1.0-w);
      t2[gemID][trid] += pow((*hits)[i]->GetHittime(),2);
      if( (*hits)[i]->GetHittime() < tmin[gemID][trid] ) tmin[gemID][trid] = (*hits)[i]->GetHittime();
      if( (*hits)[i]->GetHittime() > tmax[gemID][trid] ) tmax[gemID][trid] = (*hits)[i]->GetHittime();
      
      //Track coordinates in TRANSPORT system:
      tx[gemID][trid] = x[gemID][trid];
      ty[gemID][trid] = y[gemID][trid];
      //we now want tx and ty to be the entry point of the hit in the GEM gas layer...
      if((*hits)[i]->GetPos().z()<zin[gemID][trid]){
	xin[gemID][trid] = (*hits)[i]->GetPos().x();
	yin[gemID][trid] = (*hits)[i]->GetPos().y();
	zin[gemID][trid] = (*hits)[i]->GetPos().z();
      }
      //if((*hits)[i]->GetOutPos().z()>zout[gemID][trid]){// should I do that ???
      if((*hits)[i]->GetMom()<pmin[gemID][trid]){// or rather that ??? 
	// doesnt matter too much for signal, it may for background
	xout[gemID][trid] = (*hits)[i]->GetOutPos().x();
	yout[gemID][trid] = (*hits)[i]->GetOutPos().y();
	zout[gemID][trid] = (*hits)[i]->GetOutPos().z();
	pmin[gemID][trid] = (*hits)[i]->GetMom();
      }
      txp[gemID][trid] = txp[gemID][trid]*w +( (*hits)[i]->GetXp() )*(1.0-w);
      typ[gemID][trid] = typ[gemID][trid]*w +( (*hits)[i]->GetYp() )*(1.0-w);
      
      xg[gemID][trid] = xg[gemID][trid]*w +( (*hits)[i]->GetGlobalPos().x() )*(1.0-w);
      yg[gemID][trid] = yg[gemID][trid]*w +( (*hits)[i]->GetGlobalPos().y() )*(1.0-w);
      zg[gemID][trid] = zg[gemID][trid]*w +( (*hits)[i]->GetGlobalPos().z() )*(1.0-w);
      
      //For edep, we do the sum:
      edep[gemID][trid] += (*hits)[i]->GetEdep();
      
      nsteps_track_layer[gemID][trid]++;
    }
  }
  
  G4TrajectoryContainer *trajectorylist = evt->GetTrajectoryContainer(); //For particle history information:

  set<int> TIDs_unique; //all unique track IDs involved in GEM hits in this event (for filling particle history tree)

  for(map<int,set<int> >::iterator hit=tracks_layers.begin(); hit!=tracks_layers.end(); hit++ ){
    set<int> tracklist = hit->second;
    int gemID = hit->first;

    for(set<int>::iterator track=tracklist.begin(); track!=tracklist.end(); track++ ){
      int trackID = *track;
      if( edep[gemID][trackID] >= gemoutput.threshold ){
	gemoutput.plane.push_back( gemID );
	gemoutput.strip.push_back( 0 );
	//Difference between "x" and "tx" is that "x" is smeared by GEM coordinate resolution:
	gemoutput.x.push_back( (-y[gemID][trackID] + CLHEP::RandGauss::shoot(0.0,fGEMres) )/_L_UNIT );
	gemoutput.y.push_back( (x[gemID][trackID] + CLHEP::RandGauss::shoot(0.0,fGEMres) )/_L_UNIT );
	gemoutput.z.push_back( z[gemID][trackID]/_L_UNIT );
	gemoutput.polx.push_back( -poly[gemID][trackID] );
	gemoutput.poly.push_back(  polx[gemID][trackID] );
	gemoutput.polz.push_back(  polz[gemID][trackID] );
	gemoutput.t.push_back( t[gemID][trackID]/_T_UNIT );
	gemoutput.trms.push_back( sqrt(t2[gemID][trackID]/double(nsteps_track_layer[gemID][trackID]) - pow(t[gemID][trackID],2))/_T_UNIT );
	gemoutput.tmin.push_back( tmin[gemID][trackID]/_T_UNIT );
	gemoutput.tmax.push_back( tmax[gemID][trackID]/_T_UNIT );
	gemoutput.xin.push_back( -yin[gemID][trackID]/_L_UNIT );
	gemoutput.yin.push_back( xin[gemID][trackID]/_L_UNIT );
	gemoutput.zin.push_back( zin[gemID][trackID]/_L_UNIT );
	gemoutput.xout.push_back( -yout[gemID][trackID]/_L_UNIT );
	gemoutput.yout.push_back( xout[gemID][trackID]/_L_UNIT );
	gemoutput.zout.push_back( zout[gemID][trackID]/_L_UNIT );
	gemoutput.tx.push_back( -ty[gemID][trackID]/_L_UNIT );
	gemoutput.ty.push_back( tx[gemID][trackID]/_L_UNIT );
	gemoutput.txp.push_back( -typ[gemID][trackID] );
	gemoutput.typ.push_back( txp[gemID][trackID] );
	gemoutput.xg.push_back( xg[gemID][trackID]/_L_UNIT );
	gemoutput.yg.push_back( yg[gemID][trackID]/_L_UNIT );
	gemoutput.zg.push_back( zg[gemID][trackID]/_L_UNIT );
	gemoutput.trid.push_back( trackID );
	gemoutput.mid.push_back( mid[gemID][trackID] );
	gemoutput.pid.push_back( pid[gemID][trackID] );
	gemoutput.vx.push_back( vx[gemID][trackID]/_L_UNIT );
	gemoutput.vy.push_back( vy[gemID][trackID]/_L_UNIT );
	gemoutput.vz.push_back( vz[gemID][trackID]/_L_UNIT );
	gemoutput.p.push_back( p[gemID][trackID]/_E_UNIT );
	gemoutput.edep.push_back( edep[gemID][trackID]/_E_UNIT );
	gemoutput.beta.push_back( beta[gemID][trackID] );

	gemoutput.otridx.push_back( OTrackIndices[gemID][trackID] );
	gemoutput.ptridx.push_back( PTrackIndices[gemID][trackID] );
	gemoutput.sdtridx.push_back( SDTrackIndices[gemID][trackID] );
      
	if( trajectorylist ){ //Fill Particle History, starting with the particle itself and working all the way back to primary particles:
	  int MIDtemp = mid[gemID][trackID];
	  int TIDtemp = trackID;
	  int PIDtemp = pid[gemID][trackID];
	  int hitidx = gemoutput.nhits_GEM;
	  int nbouncetemp = 0;
	  do {
	    G4Trajectory *trajectory = (G4Trajectory*) (*trajectorylist)[TrajectoryIndex[TIDtemp]];
	    
	    PIDtemp = trajectory->GetPDGEncoding();
	    MIDtemp = MotherTrackIDs[TIDtemp];

	    std::pair<set<int>::iterator, bool > newtrajectory = TIDs_unique.insert( TIDtemp );
	    
	    if( newtrajectory.second ){ //This trajectory does not yet exist in the particle history of this detector for this event. Add it:
	      gemoutput.ParticleHistory.PID.push_back( PIDtemp );
	      gemoutput.ParticleHistory.MID.push_back( MIDtemp );
	      gemoutput.ParticleHistory.TID.push_back( TIDtemp );
	      gemoutput.ParticleHistory.hitindex.push_back( hitidx ); //Of course, this means that if a trajectory is involved in multiple hits in this detector, this variable will point to the first hit encountered only!
	      gemoutput.ParticleHistory.nbounce.push_back( nbouncetemp );
	      gemoutput.ParticleHistory.vx.push_back( (trajectory->GetPoint(0)->GetPosition() ).x()/_L_UNIT );
	      gemoutput.ParticleHistory.vy.push_back( (trajectory->GetPoint(0)->GetPosition() ).y()/_L_UNIT );
	      gemoutput.ParticleHistory.vz.push_back( (trajectory->GetPoint(0)->GetPosition() ).z()/_L_UNIT );
	      gemoutput.ParticleHistory.px.push_back( (trajectory->GetInitialMomentum() ).x()/_E_UNIT );
	      gemoutput.ParticleHistory.py.push_back( (trajectory->GetInitialMomentum() ).y()/_E_UNIT );
	      gemoutput.ParticleHistory.pz.push_back( (trajectory->GetInitialMomentum() ).z()/_E_UNIT );
	      gemoutput.ParticleHistory.npart++;
	    }
	    
	    TIDtemp = MIDtemp;
	    
	    nbouncetemp++;

	  } while( MIDtemp != 0 );
	}

	gemoutput.nhits_GEM++;
      }
    }
  }
}

void G4SBSEventAction::FillCalData( const G4Event *evt, G4SBSCalHitsCollection *hits, G4SBSCALoutput &caloutput, G4SBSSDTrackOutput &SDtracks ){
  //The "CAL" output class provides two kinds of information: 
  //1. Sum of energy deposition in a cell.
  //2. List of all particles in a cell with track ids, pids, mids, coordinates, vertices, etc. (also momentum and energy deposition)
  //3. Since the same particle can deposit energy in multiple cells, we have to consider this. 
  //4. Multiple particles can deposit energy in the same cell. 
  //5. We define a "hit" as the sum of all energy depositions in a cell per event, within a given time window that is the same for all physical placements of a given sensitive detector.
  //6. Hopefully, the geometry has been defined in such a way that each cell has a unique ID number!

  caloutput.Clear();
  // caloutput.threshold = 0.0*eV;

  //Is it reasonable to assume that all "hits" in a CalSD are chronologically ordered? Probably not, given the fact that large numbers of secondaries will be produced, and might be tracked before primary particles

  set<int> CellList; //List of all unique cells with hits
  map<int,vector<int> > steplist_cell; // list of all tracking steps in a cell
  //map<int,vector<int> > steptidx_cell; // time-ordered list of tracking steps in cell
  map<int,vector<int> > steplist_cell_timeordered; //key 1 = cell, key 2 = time ordering index, value = index in hit array
  map<int,int> nsteps_cell;
  map<int,vector<double> > xstep_cell, ystep_cell, zstep_cell, tstep_cell, edep_step_cell;
  map<int,double> tmin_cell; 
  map<int,int> Rows; //Mapping between cells and rows
  map<int,int> Cols; //Mapping between cells and columns
  map<int,int> Planes;
  map<int,int> Wires; 
  map<int,double> XCell; //Mapping between cells and x coordinates of cell centers (local)
  map<int,double> YCell; //Mapping between cells and y coordinates of cell centers (local)
  map<int,double> ZCell; //Mapping between cells and z coordinates of cell centers (local);
  map<int,double> XCellG,YCellG,ZCellG; //Mapping between cells and xyz coordinates of cell centers (global).
  map<int,int> nhits_cell; //Number of "hits" in a given cell:
  map<int,vector<int> > nsteps_hit_cell;
  map<int,vector<double> > xsum, ysum, zsum; //sum of local positions of tracks
  map<int,vector<double> > xsumg, ysumg, zsumg; //sum of global positions of tracks
  
  // ****

  caloutput.gatewidth = caloutput.timewindow;
  map<int,vector<vector<double> > > esum_tbin;

  // ****

  
  map<int,vector<double> > esum, t, t2, tmin, tmax;
  map<int,set<int> > OTrackIndices; //key = cell, value = list of all "OTracks" contributing to this hit in this cell
  map<int,set<int> > PTrackIndices; //key = cell, value = list of all "PTracks" contributing to this hit in this cell
  map<int,set<int> > SDTrackIndices; //key = cell, value = list of all "SDTracks" contributing to this hit in this cell

  map<int,set<int> > TrackIDs; //mapping between cells and a list of unique track IDs depositing energy in a cell
  map<int,map<int,int> > nsteps_track; //counting number of steps on a track
  map<int,map<int,double> > x,y,z,trt,E,trtmin,trtmax,L; //average coordinates, energy, path length for each unique track ID depositing energy in a cell:
  map<int,map<int,double> > vx,vy,vz; //production vertex coordinates of each unique track ID depositing energy in a cell
  map<int,map<int,int> > MID, PID; //mother ID and particle ID of unique tracks in each cell:
  map<int,map<int,double> > p, px, py, pz, edep; //initial momentum and total energy deposition of unique tracks in each cell:
  
  //Loop over all hits; in this loop, we want to sort tracking steps within individual cells chronologically:
  for( G4int hit=0; hit<hits->entries(); hit++ ){
    std::pair<set<int>::iterator, bool> newcell = CellList.insert( (*hits)[hit]->GetCell() );
    int cell = *(newcell.first);
    std::pair<set<int>::iterator, bool> newtrack = TrackIDs[cell].insert( (*hits)[hit]->GetTrID() );
    int track = *(newtrack.first);
    G4double Edep = (*hits)[hit]->GetEdep();
    G4int pid = (*hits)[hit]->GetPID();
    G4double steptime = (*hits)[hit]->GetTime(); //global, since start of event:
    
    if( pid != 0 && Edep > 0.0 ){ //exclude optical photons and other non-physical particles
      //Global hit information:
      if( newcell.second ){ //first hit in a new cell:
	nsteps_cell[cell] = 1;
	//nhits_cell[cell] = 1;
	Rows[cell] = (*hits)[hit]->GetRow();
	Cols[cell] = (*hits)[hit]->GetCol();
	Planes[cell] = (*hits)[hit]->GetPlane();
	Wires[cell] = (*hits)[hit]->GetWire();
	XCell[cell] = (*hits)[hit]->GetCellCoords().x();
	YCell[cell] = (*hits)[hit]->GetCellCoords().y();
	ZCell[cell] = (*hits)[hit]->GetCellCoords().z();
	XCellG[cell] = (*hits)[hit]->GetGlobalCellCoords().x();
	YCellG[cell] = (*hits)[hit]->GetGlobalCellCoords().y();
	ZCellG[cell] = (*hits)[hit]->GetGlobalCellCoords().z();
	
	steplist_cell_timeordered[cell].push_back(hit);

      // *****
	// G4cout << " **** hit **** " << hit << " ***** cell ***** " << cell << " ***** tid ***** " << track << " **** steptime **** " << steptime << endl;
      // *****
	
      } else { //additional step in existing cell:
	//double w = double(nsteps_cell[cell])/(double(nsteps_cell[cell]+1) );

	//loop over list of existing steps, find the lowest index for which steptime < tstep_j
	//G4int newstepindex = nsteps_cell[cell];
	steplist_cell_timeordered[cell].push_back(hit);
	G4int jidx = nsteps_cell[cell]-1;
	G4int jhit = steplist_cell_timeordered[cell][jidx];
	while( jidx >= 0 && steptime < (*hits)[jhit]->GetTime() ){
	  //a hit at an earlier position in the array came later than this hit:
	  steplist_cell_timeordered[cell][jidx] = hit;
	  steplist_cell_timeordered[cell][jidx+1] = jhit;
	  jidx--;
	  if( jidx < 0 ) break;
	  jhit =  steplist_cell_timeordered[cell][jidx];
	}
	tmin_cell[cell] = (*hits)[steplist_cell_timeordered[cell][0]]->GetTime();
	nsteps_cell[cell]++;

      // *****
	// G4cout << " **** hit **** " << hit << " ***** cell ***** " << cell << " ***** tid ***** " << track << " **** steptime **** " << steptime  tmin_cell[cell] << endl;
      // *****


      }
      //Track information:
      if( newtrack.second ){ //new track in this cell:
	nsteps_track[cell][track] = 1;
	x[cell][track] = Edep * (*hits)[hit]->GetPos().x(); //local
	y[cell][track] = Edep * (*hits)[hit]->GetPos().y(); //
	z[cell][track] = Edep * (*hits)[hit]->GetPos().z(); 
	trt[cell][track] = Edep * (*hits)[hit]->GetTime();
	E[cell][track] = (*hits)[hit]->GetEnergy();
	trtmin[cell][track] =  (*hits)[hit]->GetTime();
	trtmax[cell][track] =  (*hits)[hit]->GetTime();
	L[cell][track] = (*hits)[hit]->GetLstep(); //path length
	vx[cell][track] = (*hits)[hit]->GetVertex().x(); 
	vy[cell][track] = (*hits)[hit]->GetVertex().y(); 
	vz[cell][track] = (*hits)[hit]->GetVertex().z();
	MID[cell][track] = (*hits)[hit]->GetMID();
	PID[cell][track] = pid;
	p[cell][track] = (*hits)[hit]->GetMomentum().mag();
	px[cell][track] = (*hits)[hit]->GetMomentum().x();
	py[cell][track] = (*hits)[hit]->GetMomentum().y();
	pz[cell][track] = (*hits)[hit]->GetMomentum().z();
	edep[cell][track] = Edep;

	//This is the appropriate place to add this info, since each track in a CALSD
	//can only have exactly one set of otrack, ptrack, and sdtrack info, regardless of how many
	//steps
	OTrackIndices[cell].insert( SDtracks.otracklist[(*hits)[hit]->GetOTrIdx()] );
	PTrackIndices[cell].insert( SDtracks.ptracklist[(*hits)[hit]->GetPTrIdx()] );
	SDTrackIndices[cell].insert( SDtracks.sdtracklist[(*hits)[hit]->GetSDTrIdx()][hits->GetSDname()] );

	// G4cout << "Inserted SD track for SD " << hits->GetSDname() << ", cell ID, TID = " << cell << ", " << (*hits)[hit]->GetTrID()
	//        << ", SDTID = " << (*hits)[hit]->GetSDTrIdx()
	//        << ", SDTidx = " << SDtracks.sdtracklist[(*hits)[hit]->GetSDTrIdx()][hits->GetSDname()] << G4endl;
	
      } else { //additional step in this cell:
	//double w = double(nsteps_track[cell][track])/(double(nsteps_track[cell][track]+1) );
	// x[cell][track] = w * x[cell][track] + (1.0-w)*(*hits)[hit]->GetPos().x(); //local
	// y[cell][track] = w * y[cell][track] + (1.0-w)*(*hits)[hit]->GetPos().y();
	// z[cell][track] = w * z[cell][track] + (1.0-w)*(*hits)[hit]->GetPos().z();
	// trt[cell][track] = w * trt[cell][track] + (1.0-w)*(*hits)[hit]->GetTime();
	x[cell][track] += Edep * (*hits)[hit]->GetPos().x();
	y[cell][track] += Edep * (*hits)[hit]->GetPos().y();
	z[cell][track] += Edep * (*hits)[hit]->GetPos().z();
	trt[cell][track] += Edep * (*hits)[hit]->GetTime();
	if( (*hits)[hit]->GetTime() < trtmin[cell][track] ) trtmin[cell][track] = (*hits)[hit]->GetTime();
	if( (*hits)[hit]->GetTime() > trtmax[cell][track] ) trtmax[cell][track] = (*hits)[hit]->GetTime();
	L[cell][track] += (*hits)[hit]->GetLstep();
	edep[cell][track] += Edep;
      }
      
    }
  }
  
  set<int> TIDs_unique;
  G4TrajectoryContainer *trajectorylist = evt->GetTrajectoryContainer(); //For particle history information:

  //int hitindex=0;
  //map<int, set<int> > tracklist_hit; 

  TClonesArray *histpstemp = fIO->PulseShape_histograms;
  TClonesArray *histesumtemp = fIO->Esum_histograms;
  
  TH1F *hpstemp = ( (TH1F*) (*histpstemp)[fhistogram_index]);
  TH1F *hesumtemp = ( (TH1F*) (*histesumtemp)[fhistogram_index]);
  
  G4double esum_total = 0.0;
  
  //Now loop over all unique cells and tracks and fill CALoutput data structure:
  for( set<int>::iterator itcell=CellList.begin(); itcell != CellList.end(); itcell++ ){
    int cell = *itcell;

    vector<int> steplist = steplist_cell_timeordered[cell];

    //int hitindex=-1;

    nhits_cell[cell] = 0;

    map<int,int> firsthit_track;

    double tbin = 0.0;

    for( int istep=0; istep<steplist.size(); istep++ ){
      int jhit = steplist_cell_timeordered[cell][istep];
      G4double tstep = (*hits)[jhit]->GetTime();
      G4double estep = (*hits)[jhit]->GetEdep();
      G4double xstep = (*hits)[jhit]->GetPos().x();
      G4double ystep = (*hits)[jhit]->GetPos().y();
      G4double zstep = (*hits)[jhit]->GetPos().z();
      G4double xgstep = (*hits)[jhit]->GetLabPos().x();
      G4double ygstep = (*hits)[jhit]->GetLabPos().y();
      G4double zgstep = (*hits)[jhit]->GetLabPos().z();
      
      G4int pid = (*hits)[jhit]->GetPID();
      G4int hitindex = nhits_cell[cell] > 0 ? nhits_cell[cell]-1 : 0;

      G4int tidstep = (*hits)[jhit]->GetTrID();
      
      // G4bool newhit = false;
      
      if( istep == 0 || (nhits_cell[cell] > 0 && tstep > tmin[cell][hitindex] + caloutput.timewindow ) ){
	//This is either the first hit or a tracking step that fell outside the timing window (i.e., "gate") defined for this SD:
	nhits_cell[cell]++;
	nsteps_hit_cell[cell].push_back( 1 );
	//all quantities that are summed over the hit are energy-deposition-weighted:
	xsumg[cell].push_back( xgstep*estep );
	ysumg[cell].push_back( ygstep*estep );
	zsumg[cell].push_back( zgstep*estep );

	xsum[cell].push_back( xstep*estep );
	ysum[cell].push_back( ystep*estep );
	zsum[cell].push_back( zstep*estep );

	esum[cell].push_back( estep );
	t[cell].push_back( tstep*estep );
	t2[cell].push_back( pow(tstep,2)*estep );
	//tmin and tmax values are unweighted:
	tmin[cell].push_back( tstep );
	tmax[cell].push_back( tstep );

	// ******
	vector<double> esum_temp( caloutput.ntimebins );
	esum_tbin[cell].push_back( esum_temp );

	esum_tbin[cell][nhits_cell[cell] - 1][0] += estep ;

	// ******
	
      } else { //Add this step to the current hit:
	xsumg[cell][hitindex] += estep * xgstep;
	ysumg[cell][hitindex] += estep * ygstep;
	zsumg[cell][hitindex] += estep * zgstep;

	xsum[cell][hitindex] += estep * xstep;
	ysum[cell][hitindex] += estep * ystep;
	zsum[cell][hitindex] += estep * zstep;

	esum[cell][hitindex] += estep;
	t[cell][hitindex] += estep * tstep;
	t2[cell][hitindex] += estep * pow(tstep,2);
	tmin[cell][hitindex] = (tstep < tmin[cell][hitindex] ) ? tstep : tmin[cell][hitindex];
	tmax[cell][hitindex] = (tstep > tmax[cell][hitindex] ) ? tstep : tmax[cell][hitindex];

	// ******
	double wtbin = ( caloutput.timewindow - 0.0 )/double(caloutput.ntimebins);
      	int bin_tstep = int( (tstep - tmin[cell][hitindex])/wtbin );
	if ( bin_tstep >= 0 && bin_tstep < caloutput.ntimebins ) esum_tbin[cell][hitindex][bin_tstep] += estep;
	// ******

	nsteps_hit_cell[cell][hitindex]++;
      }

      hpstemp->Fill( tstep - tmin[cell][hitindex], estep );
      esum_total += estep;

      firsthit_track.insert( std::pair<int,int>( tidstep, nhits_cell[cell]-1 ) );
      
      //tracklist_hit[nhits_cell[cell]-1].insert(tidstep);
      
    }

    
    //G4int ngoodhits=0;

    map<int,int> goodhit_index;

    // G4double esum_total = 0.0;

    //G4cout << " Hello World " << endl;
    
    for( int ihit=0; ihit<nhits_cell[cell]; ihit++ ){
      if( esum[cell][ihit] >= caloutput.threshold ){
	//HitList[cell].insert( caloutput.nhits_CAL );
	caloutput.cell.push_back(cell);
	caloutput.row.push_back(Rows[cell]);
	caloutput.col.push_back(Cols[cell]);
	caloutput.plane.push_back(Planes[cell]);
	caloutput.wire.push_back(Wires[cell]);
	caloutput.xcell.push_back( XCell[cell]/_L_UNIT );
	caloutput.ycell.push_back( YCell[cell]/_L_UNIT );
	caloutput.zcell.push_back( ZCell[cell]/_L_UNIT );
	caloutput.xcellg.push_back( XCellG[cell]/_L_UNIT );
	caloutput.ycellg.push_back( YCellG[cell]/_L_UNIT );
	caloutput.zcellg.push_back( ZCellG[cell]/_L_UNIT );
	caloutput.sumedep.push_back( esum[cell][ihit]/_E_UNIT );

	caloutput.tavg.push_back( t[cell][ihit]/esum[cell][ihit]/_T_UNIT );
	caloutput.trms.push_back( sqrt( t2[cell][ihit]/esum[cell][ihit] - pow(t[cell][ihit]/esum[cell][ihit],2) )/_T_UNIT );
	caloutput.xhit.push_back( xsum[cell][ihit]/esum[cell][ihit]/_L_UNIT );
	caloutput.yhit.push_back( ysum[cell][ihit]/esum[cell][ihit]/_L_UNIT );
	caloutput.zhit.push_back( zsum[cell][ihit]/esum[cell][ihit]/_L_UNIT );

	caloutput.xhitg.push_back( xsumg[cell][ihit]/esum[cell][ihit]/_L_UNIT );
	caloutput.yhitg.push_back( ysumg[cell][ihit]/esum[cell][ihit]/_L_UNIT );
	caloutput.zhitg.push_back( zsumg[cell][ihit]/esum[cell][ihit]/_L_UNIT );
	
	// caloutput.tavg.push_back( t[cell]/_T_UNIT );
	// caloutput.trms.push_back( sqrt( t2[cell]/double(nsteps_cell[cell]) - pow(t[cell],2) )/_T_UNIT );
	caloutput.tmin.push_back( tmin[cell][ihit]/_T_UNIT );
	caloutput.tmax.push_back( tmax[cell][ihit]/_T_UNIT );

	// ************* ++++++ ************
	for ( int itbin=0; itbin<caloutput.ntimebins; itbin++ ){
	  esum_tbin[cell][ihit][itbin] /= _E_UNIT;
	}
	caloutput.edep_vs_time.push_back( esum_tbin[cell][ihit] );
	// *****
	
	//If there are multiple Otracks contributing to this hit, choose the one with the highest total energy:
	G4double maxE = 0.0;
	G4bool firsttrack=true;
	int otridx_final=-1;
	for( set<int>::iterator iotrk=OTrackIndices[cell].begin(); iotrk!=OTrackIndices[cell].end(); ++iotrk ){
	  G4double Eotrack = SDtracks.oenergy[*iotrk];
	  if( firsttrack || Eotrack > maxE ){
	    otridx_final = *iotrk;
	    maxE = Eotrack;
	    firsttrack = false;
	  }
	}

	//Primary tracks:
	maxE = 0.0;
	firsttrack = true;
	int ptridx_final=-1;
	for( set<int>::iterator iptrk=PTrackIndices[cell].begin(); iptrk!=PTrackIndices[cell].end(); ++iptrk ){
	  G4double Eptrack = SDtracks.penergy[*iptrk];
	  if( firsttrack || Eptrack > maxE ){
	    ptridx_final = *iptrk;
	    maxE = Eptrack;
	    firsttrack = false;
	  }
	}

	//SD boundary crossing tracks:
	maxE = 0.0;
	firsttrack = true;
	int sdtridx_final=-1;
	for( set<int>::iterator isdtrk=SDTrackIndices[cell].begin(); isdtrk!=SDTrackIndices[cell].end(); ++isdtrk ){
	  int sdidxtemp = *isdtrk;

	  if( sdidxtemp >= 0 && sdidxtemp <SDtracks.sdenergy.size() ){
	    G4double Esdtrack = SDtracks.sdenergy[sdidxtemp];
	    if( firsttrack || Esdtrack > maxE ){
	      sdtridx_final = *isdtrk;
	      maxE = Esdtrack;
	      firsttrack = false;
	    }
	  }
	}

	caloutput.otridx.push_back( otridx_final );
	caloutput.ptridx.push_back( ptridx_final );
	caloutput.sdtridx.push_back( sdtridx_final );
	
	goodhit_index[ihit] = caloutput.nhits_CAL;
	
	caloutput.nhits_CAL++;
      }

      // esum_total += esum[cell][ihit];
    }

    if( caloutput.nhits_CAL > 0 ){
      hesumtemp->Fill( esum_total );
    
      caloutput.Esum = esum_total/_E_UNIT;
    } else {
      caloutput.Esum = -100.0;
    }
    
    for( set<int>::iterator itrack = TrackIDs[cell].begin(); itrack != TrackIDs[cell].end(); itrack++ ){
      int track = *itrack;

      if( goodhit_index.find( firsthit_track[track] ) != goodhit_index.end() ){ //Track must be associated with at least one "GOOD" hit in this cell to be recorded:
	
	caloutput.ihit.push_back( goodhit_index[firsthit_track[track]] );
	caloutput.x.push_back( x[cell][track]/edep[cell][track]/_L_UNIT );
	caloutput.y.push_back( y[cell][track]/edep[cell][track]/_L_UNIT );
	caloutput.z.push_back( z[cell][track]/edep[cell][track]/_L_UNIT );
	caloutput.t.push_back( trt[cell][track]/edep[cell][track]/_T_UNIT );
	caloutput.dt.push_back( (trtmax[cell][track] - trtmin[cell][track])/_T_UNIT );
	caloutput.E.push_back( E[cell][track]/_E_UNIT );
	caloutput.L.push_back( L[cell][track]/_L_UNIT );
	caloutput.vx.push_back( vx[cell][track]/_L_UNIT );
	caloutput.vy.push_back( vy[cell][track]/_L_UNIT );
	caloutput.vz.push_back( vz[cell][track]/_L_UNIT );
	caloutput.mid.push_back( MID[cell][track] );
	caloutput.pid.push_back( PID[cell][track] );
	caloutput.trid.push_back( track );
	caloutput.p.push_back( p[cell][track]/_E_UNIT );
	caloutput.px.push_back( px[cell][track]/_E_UNIT );
	caloutput.py.push_back( py[cell][track]/_E_UNIT );
	caloutput.pz.push_back( pz[cell][track]/_E_UNIT );
	caloutput.edep.push_back( edep[cell][track]/_E_UNIT );
	caloutput.npart_CAL++;
	
	if( trajectorylist ){ //Fill Particle History, starting with the particle itself and working all the way back to primary particles:
	  int MIDtemp = MID[cell][track];
	  int TIDtemp = track;
	  int PIDtemp = PID[cell][track];
	  int hitidx = caloutput.nhits_CAL;
	  int nbouncetemp = 0;
	  do {
	    G4Trajectory *trajectory = (G4Trajectory*) (*trajectorylist)[TrajectoryIndex[TIDtemp]];
	    
	    PIDtemp = trajectory->GetPDGEncoding();
	    MIDtemp = MotherTrackIDs[TIDtemp];
	    
	    std::pair<set<int>::iterator, bool > newtrajectory = TIDs_unique.insert( TIDtemp );
	    
	    if( newtrajectory.second ){ //This trajectory does not yet exist in the particle history of this detector for this event. Add it:
	      caloutput.ParticleHistory.PID.push_back( PIDtemp );
	      caloutput.ParticleHistory.MID.push_back( MIDtemp );
	      caloutput.ParticleHistory.TID.push_back( TIDtemp );
	      caloutput.ParticleHistory.hitindex.push_back( hitidx ); //Of course, this means that if a trajectory is involved in multiple hits in this detector, this variable will point to the first hit encountered only!
	      caloutput.ParticleHistory.nbounce.push_back( nbouncetemp );
	      caloutput.ParticleHistory.vx.push_back( (trajectory->GetPoint(0)->GetPosition() ).x()/_L_UNIT );
	      caloutput.ParticleHistory.vy.push_back( (trajectory->GetPoint(0)->GetPosition() ).y()/_L_UNIT );
	      caloutput.ParticleHistory.vz.push_back( (trajectory->GetPoint(0)->GetPosition() ).z()/_L_UNIT );
	      caloutput.ParticleHistory.px.push_back( (trajectory->GetInitialMomentum() ).x()/_E_UNIT );
	      caloutput.ParticleHistory.py.push_back( (trajectory->GetInitialMomentum() ).y()/_E_UNIT );
	      caloutput.ParticleHistory.pz.push_back( (trajectory->GetInitialMomentum() ).z()/_E_UNIT );
	      caloutput.ParticleHistory.npart++;
	    }
	    
	    TIDtemp = MIDtemp;
	    
	    nbouncetemp++;
	    
	  } while( MIDtemp != 0 );
	}
      }
    }
  }
}

void G4SBSEventAction::FillECalData( G4SBSECalHitsCollection *hits, G4SBSECaloutput &ecaloutput, G4SBSSDTrackOutput &SDtracks )
{
  
  int nG4hits = hits->entries(); //Loop over nG4hits
  ecaloutput.Clear();

  set<int> TIDs_unique;
  set<int> PMTs_unique;

  map<int,bool> Photon_used;
  map<int,bool> Photon_detected;
  map<int,double> Photon_energy;
  map<int,double> Photon_hittime;
  map<int,int> Photon_PMT;
  map<int,int> Photon_row;
  map<int,int> Photon_col;
  map<int,int> Photon_plane;
  map<int,int> Photon_nsteps;
  map<int,G4ThreeVector> Photon_xpmt; //"local" xy and "global" (xyz) PMT coordinates
  map<int,G4ThreeVector> Photon_xgpmt;
  map<int,int> Photon_otridx;
  map<int,int> Photon_ptridx;
  map<int,int> Photon_sdtridx;

  map<int,int> PMT_Numphotoelectrons;
  map<int,int> PMT_row;
  map<int,int> PMT_col;
  map<int,int> PMT_plane;
  map<int,G4ThreeVector> PMT_x, PMT_xg; //"local" xy and "global" (xyz) PMT coordinates
  map<int,double> PMT_hittime;
  map<int,double> PMT_hittime2; //hit time squared.
  map<int,double> PMT_rmstime;
  map<int,double> PMT_tmin;
  map<int,double> PMT_tmax;
  //Let's handle optical photon detectors like ECAL (and RICH) similarly to CALSD: 
  map<int,set<int> > OTrackIndices;
  map<int,set<int> > PTrackIndices;
  map<int,set<int> > SDTrackIndices;

  
  // *****
  map<int,vector<int> > steplist_pmt_timeordered; //key 1 = pmt, key 2 = time ordering index, value = index in hit array
  set<int> PmtList;
  map<int,int> nsteps_pmt;
  map<int,int> nhits_pmt;
  map<int,vector<double> > tmin_pmt, tmin;
  map<int,vector<vector<double> > > npe_tbin;

  ecaloutput.gatewidth = ecaloutput.timewindow;
    
  //G4cout << " ******** timewindow ********* " << ecaloutput.timewindow << endl;
  //G4cout << " ******** threshold ********* " << ecaloutput.threshold << endl;
  //G4cout << " ******** ntimebins ********* " << ecaloutput.ntimebins << endl;
  // *****
  
  //G4MaterialPropertiesTable *MPT; 

  for( int step = 0; step < nG4hits; step++ ){
    //Retrieve all relevant information for this step:
    int pmt = (*hits)[step]->GetPMTnumber();
    int row = (*hits)[step]->Getrownumber();
    int col = (*hits)[step]->Getcolnumber();
    int plane = (*hits)[step]->Getplanenumber();
    int tid = (*hits)[step]->GetTrackID();
    double Ephoton = (*hits)[step]->Getenergy();
    double Hittime = (*hits)[step]->GetTime();

    G4ThreeVector xpmt = (*hits)[step]->GetCellCoords();
    G4ThreeVector xgpmt = (*hits)[step]->GetGlobalCellCoords();

    G4double QEphoton = (*hits)[step]->GetQuantumEfficiency();

    int otridx = SDtracks.otracklist[(*hits)[step]->GetOTrIdx()];
    int ptridx = SDtracks.ptracklist[(*hits)[step]->GetPTrIdx()];
    int sdtridx = SDtracks.sdtracklist[(*hits)[step]->GetSDTrIdx()][hits->GetSDname()];

    // *****

    // To get the time ordered list of tracks in a given PMT
    std::pair<set<int>::iterator, bool> newpmt = PmtList.insert( pmt );
    int pmt_un = *(newpmt.first);
    bool photon_det = G4UniformRand() <= QEphoton;
      
    
    if( newpmt.second ){ //first hit in a new pmt:


      if( photon_det ){ // I assume we only care about detected photons in this case
	nsteps_pmt[pmt_un] = 1;
	// PmtList.insert( pmt );
	steplist_pmt_timeordered[pmt_un].push_back(step);

      }

    }else { // additional hits in existing PMT

      if( photon_det ){

	steplist_pmt_timeordered[pmt_un].push_back(step);
	// it seems unnecessary to do the time ordering
	G4int jidx = nsteps_pmt[pmt_un]-1;
	G4int jhit = steplist_pmt_timeordered[pmt_un][jidx];
	while( jidx >= 0 && Hittime < (*hits)[jhit]->GetTime() ){
	  //a hit at an earlier position in the array came later than this hit:
	  steplist_pmt_timeordered[pmt_un][jidx] = step;
	  steplist_pmt_timeordered[pmt_un][jidx+1] = jhit;
	  jidx--;
	  if( jidx < 0 ) break;
	  jhit =  steplist_pmt_timeordered[pmt_un][jidx];
	}
	// tmin_pmt[pmt] = (*hits)[steplist_pmt_timeordered[pmt][0]]->GetTime();
	nsteps_pmt[pmt_un]++;
      }     
     
    }

    // *****



    
    //Following the method implemented during the RICH routine
    std::pair< set<int>::iterator, bool > photontrack = TIDs_unique.insert( tid );
   
    if( photontrack.second ){
     
      bool photon_detected = G4UniformRand() <= QEphoton;

      Photon_used[ tid ] = !photon_detected; //If the photon is not detected, then we mark it as used. Otherwise, we mark it as unused, and it will be added to a PMT later.
      Photon_detected[ tid ] = photon_detected;
      Photon_energy[ tid ] = Ephoton;
      Photon_hittime[ tid ] = Hittime;
     
      Photon_PMT[ tid ] = pmt;
      Photon_row[ tid ] = row;
      Photon_col[ tid ] = col;
      Photon_plane[ tid ] = plane;

      Photon_xpmt[tid] = xpmt;
      Photon_xgpmt[tid] = xgpmt;

      Photon_nsteps[ tid ] = 1;

      Photon_otridx[ tid ] = otridx;
      Photon_ptridx[ tid ] = ptridx;
      Photon_sdtridx[ tid ] = sdtridx;      

    } else { //existing photon, additional step. Increment averages of position, direction, time, etc for all steps of a detected photon. Don't bother for 
      //undetected photons...
      if( Photon_detected[ tid ] ){
	G4double average_time = (Photon_nsteps[ tid ] * Photon_hittime[ tid ] + Hittime )/double( Photon_nsteps[tid] + 1 );
	Photon_hittime[tid] = average_time; 
	Photon_nsteps[tid] += 1;

      }
    }
  }   
  
  bool  remaining_hits = true;
  // Used to make sure we only save the track info only once
  // when saving all tracks to the tree
  bool saved_track_info = false;
  
  while( remaining_hits ) {
    
    remaining_hits = false;
    
    PMTs_unique.clear();

    for( set<int>::iterator it=TIDs_unique.begin(); it != TIDs_unique.end(); it++ ){
      int tid = *it;

      
      // Save particle info for all tracks (doesn't matter if they were
      // ultimately not detected)
      if(!saved_track_info) {
        ecaloutput.npart_ECAL++;
        ecaloutput.part_PMT.push_back( Photon_PMT[tid] );
        ecaloutput.trid.push_back( tid );
        ecaloutput.E.push_back( Photon_energy[tid]/CLHEP::eV );
        ecaloutput.t.push_back( Photon_hittime[tid]/CLHEP::ns );
        ecaloutput.detected.push_back( Photon_detected[tid] );
      }

      if( Photon_detected[ tid ] && !(Photon_used[ tid ] ) ){
	
	std::pair<set<int>::iterator,bool> testpmt = PMTs_unique.insert( Photon_PMT[tid] );

	int pmt = Photon_PMT[tid];

	if( testpmt.second ){ // new PMT;
	  
	  //Mark this photon track as used:
	  Photon_used[ tid ] = true;
	  
	  PMT_Numphotoelectrons[ pmt ] = 1;
	  PMT_row[ pmt ] = Photon_row[tid];
	  PMT_col[ pmt ] = Photon_col[tid];
	  PMT_plane[ pmt ] = Photon_plane[tid];
	  PMT_hittime[ pmt ] = Photon_hittime[tid];
	  PMT_hittime2[ pmt ] = pow(Photon_hittime[tid],2);
	  PMT_tmin[ pmt ] = Photon_hittime[tid];
	  PMT_tmax[ pmt ] = Photon_hittime[tid];
	  PMT_x[ pmt ] = Photon_xpmt[tid];
	  PMT_xg[ pmt ] = Photon_xgpmt[tid];

	  OTrackIndices[ pmt ].insert( Photon_otridx[ tid ] );
	  PTrackIndices[ pmt ].insert( Photon_ptridx[ tid ] );
	  SDTrackIndices[ pmt ].insert( Photon_sdtridx[ tid ] );
				       
	  
	} else if( fabs( Photon_hittime[tid] - PMT_tmin[pmt] ) <= ecaloutput.timewindow ){ //Existing pmt with multiple photon detections:
	  G4double average_hittime = (PMT_Numphotoelectrons[pmt] * PMT_hittime[pmt] + Photon_hittime[tid])/double(PMT_Numphotoelectrons[pmt] + 1 );
	  PMT_hittime[pmt] = average_hittime;
	  G4double average_hittime2 = (PMT_Numphotoelectrons[pmt] * PMT_hittime2[pmt] + pow(Photon_hittime[tid],2))/double(PMT_Numphotoelectrons[pmt] + 1 );
	  PMT_hittime2[pmt] = average_hittime2;

	  PMT_Numphotoelectrons[pmt] += 1;

	  if( Photon_hittime[tid] > PMT_tmax[pmt] ) PMT_tmax[pmt] = Photon_hittime[tid];
	  if( Photon_hittime[tid] < PMT_tmin[pmt] ) PMT_tmin[pmt] = Photon_hittime[tid];
	  
	  Photon_used[tid] = true;

	  OTrackIndices[ pmt ].insert( Photon_otridx[ tid ] );
	  PTrackIndices[ pmt ].insert( Photon_ptridx[ tid ] );
	  SDTrackIndices[ pmt ].insert( Photon_sdtridx[ tid ] );
	}

	
	//If any photon is detected but not used, then remaining hits = true!
	if( !(Photon_used[tid] ) ) remaining_hits = true;
      }
    }


    //Now add hits to the output following the RICH example..
    for( set<int>::iterator it=PMTs_unique.begin(); it!=PMTs_unique.end(); it++ ){
           
      int pmt = *it;


      // *****
      // Storing pulse shape information
      vector<int> steplist = steplist_pmt_timeordered[pmt];	
      nhits_pmt[pmt] = 0;		
	
      for( int istep=0; istep<steplist.size(); istep++ ){
        int jhit = steplist_pmt_timeordered[pmt][istep];
        G4double tstep = (*hits)[jhit]->GetTime();
        G4int hitindex = nhits_pmt[pmt] > 0 ? nhits_pmt[pmt]-1 : 0;

        if( istep == 0 || (nhits_pmt[pmt] > 0 && tstep > tmin[pmt][hitindex] + ecaloutput.timewindow ) ){
          nhits_pmt[pmt]++;


          tmin[pmt].push_back( tstep );

          vector<double> npe_temp( ecaloutput.ntimebins );
          npe_tbin[pmt].push_back( npe_temp );

          npe_tbin[pmt][nhits_pmt[pmt] - 1][0] += 1. ;

        } else {

          double wtbin = ( ecaloutput.timewindow - 0.0 )/double(ecaloutput.ntimebins);
          int bin_tstep = int( (tstep - tmin[pmt][hitindex])/wtbin );
          if ( bin_tstep >= 0 && bin_tstep < ecaloutput.ntimebins ) npe_tbin[pmt][hitindex][bin_tstep] += 1.;
	  
        }
      } 
	  

      // *****

      
      if( PMT_Numphotoelectrons[pmt] >= ecaloutput.threshold ){
	
	(ecaloutput.nhits_ECal)++;
	
	ecaloutput.PMTnumber.push_back( pmt );
	ecaloutput.row.push_back( PMT_row[pmt] );
	ecaloutput.col.push_back( PMT_col[pmt] );
	ecaloutput.plane.push_back( PMT_plane[pmt] );
	ecaloutput.xcell.push_back( PMT_x[pmt].x()/_L_UNIT );
	ecaloutput.ycell.push_back( PMT_x[pmt].y()/_L_UNIT );
	ecaloutput.zcell.push_back( PMT_x[pmt].z()/_L_UNIT );
	ecaloutput.xgcell.push_back( PMT_xg[pmt].x()/_L_UNIT );
	ecaloutput.ygcell.push_back( PMT_xg[pmt].y()/_L_UNIT );
	ecaloutput.zgcell.push_back( PMT_xg[pmt].z()/_L_UNIT );
	ecaloutput.NumPhotoelectrons.push_back( PMT_Numphotoelectrons[pmt] ); //NumPhotoelectrons is vector<int> defined in ECaloutput.hh
	ecaloutput.Time_avg.push_back( PMT_hittime[pmt]/_T_UNIT );
	ecaloutput.Time_rms.push_back( sqrt(fabs(PMT_hittime2[pmt] - pow(PMT_hittime[pmt],2) ) )/_T_UNIT );	
	ecaloutput.Time_min.push_back( PMT_tmin[pmt]/_T_UNIT );
	ecaloutput.Time_max.push_back( PMT_tmax[pmt]/_T_UNIT );

	// *****
	for( int ihit=0; ihit<nhits_pmt[pmt]; ihit++ ){
	    ecaloutput.NPE_vs_time.push_back( npe_tbin[pmt][ihit] );
	}
	// *****

	//If there are multiple Otracks contributing to this hit, choose the one with the highest total energy:
	G4double maxE = 0.0;
	G4bool firsttrack=true;
	int otridx_final=-1;
	for( set<int>::iterator iotrk=OTrackIndices[pmt].begin(); iotrk!=OTrackIndices[pmt].end(); ++iotrk ){
	  G4double Eotrack = SDtracks.oenergy[*iotrk];
	  if( firsttrack || Eotrack > maxE ){
	    otridx_final = *iotrk;
	    maxE = Eotrack;
	    firsttrack = false;
	  }
	}

	//Primary tracks:
	maxE = 0.0;
	firsttrack = true;
	int ptridx_final=-1;
	for( set<int>::iterator iptrk=PTrackIndices[pmt].begin(); iptrk!=PTrackIndices[pmt].end(); ++iptrk ){
	  G4double Eptrack = SDtracks.penergy[*iptrk];
	  if( firsttrack || Eptrack > maxE ){
	    ptridx_final = *iptrk;
	    maxE = Eptrack;
	    firsttrack = false;
	  }
	}

	//SD boundary crossing tracks:
	maxE = 0.0;
	firsttrack = true;
	int sdtridx_final=-1;
	for( set<int>::iterator isdtrk=SDTrackIndices[pmt].begin(); isdtrk!=SDTrackIndices[pmt].end(); ++isdtrk ){
	  int sdidxtemp = *isdtrk;

	  if( sdidxtemp >= 0 && sdidxtemp <SDtracks.sdenergy.size() ){
	    G4double Esdtrack = SDtracks.sdenergy[sdidxtemp];
	    if( firsttrack || Esdtrack > maxE ){
	      sdtridx_final = *isdtrk;
	      maxE = Esdtrack;
	      firsttrack = false;
	    }
	  }
	}

	ecaloutput.otridx.push_back( otridx_final );
	ecaloutput.ptridx.push_back( ptridx_final );
	ecaloutput.sdtridx.push_back( sdtridx_final );
	
      }
    } //for    
  }//while
}//void

void G4SBSEventAction::FillRICHData( const G4Event *evt, G4SBSRICHHitsCollection *hits, G4SBSRICHoutput &richoutput, G4SBSSDTrackOutput &SDtracks ){
  //Here is where we traverse the hit collection of the RICH and extract useful output data. 
  
  // set<int> TIDs; //list of all unique track IDs  
  //This is the total number of tracking steps in our sensitive volume!
  int nG4hits = hits->entries();

  //cout << "Filling RICH data, nhits = " << nG4hits << endl;

  richoutput.Clear();

  set<int> TIDs_unique;  //set of all unique photon tracks involved in RICH hits.
  set<int> PMTs_unique;  //set of all unique PMTs with detected photons.
  set<int> mTIDs_unique; //set of all unique mother track IDs associated with RICH detected photons
  map<int,int> Nphe_mTID; // count the number of photoelectrons generated by a given mother track (compromise?)
  
  map<int,bool> Photon_used;
  map<int,bool> Photon_detected;
  map<int,double> Photon_energy;
  map<int,double> Photon_hittime;
  map<int,G4ThreeVector> Photon_position; // position at hit
  map<int,G4ThreeVector> Photon_direction; // direction at hit
  map<int,G4ThreeVector> Photon_vertex;  //emission vertex
  map<int,G4ThreeVector> Photon_vdirection;  //emission direction
  map<int,G4ThreeVector> Photon_xpmt; //Local coordinates
  map<int,G4ThreeVector> Photon_xgpmt; //Global coordinates
  map<int,int> Photon_PMT;
  map<int,int> Photon_row;
  map<int,int> Photon_col;
  map<int,int> Photon_mTID;
  map<int,int> Photon_origvol;
  map<int,int> Photon_nsteps;
  map<int,int> Photon_otridx;
  map<int,int> Photon_ptridx;
  map<int,int> Photon_sdtridx;

  map<int,int> PMT_Numphotoelectrons;
  map<int,int> PMT_row;
  map<int,int> PMT_col;
  map<int,double> PMT_hittime;
  map<int,double> PMT_hittime2; //hit time squared.
  map<int,double> PMT_rmstime;
  map<int,double> PMT_tmin; //earliest time
  map<int,double> PMT_tmax; //latest time:
  map<int,int> PMT_mTID;
  map<int,G4ThreeVector> PMT_pos;
  map<int,G4ThreeVector> PMT_dir;
  map<int,G4ThreeVector> PMT_vpos;
  map<int,G4ThreeVector> PMT_vdir;
  map<int,G4ThreeVector> PMT_x; //local coordinates
  map<int,G4ThreeVector> PMT_xg; //global coordinates
  map<int,int> PMT_origvol;

  //Let's handle optical photon detectors like ECAL (and RICH) similarly to CALSD: 
  map<int,set<int> > OTrackIndices;
  map<int,set<int> > PTrackIndices;
  map<int,set<int> > SDTrackIndices;
  
  //G4MaterialPropertiesTable *MPT;

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
    double Ephoton = (*hits)[step]->GetEnergy();
    double Hittime = (*hits)[step]->GetTime();
    int origin_volume = (*hits)[step]->GetOriginVol();
    G4ThreeVector vertex = (*hits)[step]->GetVertex();
    G4ThreeVector vertexdirection = (*hits)[step]->GetVertexDirection();
    G4ThreeVector position = (*hits)[step]->GetPos();
    G4ThreeVector direction = (*hits)[step]->GetDirection();
    G4ThreeVector Lposition = (*hits)[step]->GetLPos();
    G4ThreeVector Ldirection = (*hits)[step]->GetLDirection();
    G4ThreeVector xpmt = (*hits)[step]->GetCellCoord(); //"local" PMT coordinate (to e.g., detector mother volume)
    G4ThreeVector xgpmt = (*hits)[step]->GetGlobalCellCoord(); //global PMT coordinate

    G4double QEphoton = (*hits)[step]->GetQuantumEfficiency();
    
    //First, we must ask: Is this a "new" photon track? It **should be** theoretically impossible for the same photon to 
    //be detected in two different PMTs, since the photocathodes don't have the necessary properties defined for the 
    //propagation of optical photons (RINDEX). On the other hand, it IS possible, and even likely, for two or more photons
    //to be detected by the same PMT. Therefore, we should consider photon tracks at the top level of the sorting logic:
      
    std::pair< set<int>::iterator, bool > testphotontrack = TIDs_unique.insert( tid );
    
    if( testphotontrack.second ){ //New photon track: determine whether this photon is detected:
      
      bool photon_detected = G4UniformRand() <= QEphoton;
      
      Photon_used[ tid ] = !photon_detected; //If the photon is not detected, then we mark it as used. Otherwise, we mark it as unused, and it will be added to a PMT later.
      Photon_detected[ tid ] = photon_detected;
      Photon_energy[ tid ] = Ephoton;
      Photon_hittime[ tid ] = Hittime;
      Photon_position[ tid ] = position;
      Photon_direction[ tid ] = direction;
      Photon_vertex[ tid ] = vertex;
      Photon_vdirection[ tid ] = vertexdirection;
      Photon_xpmt[ tid ] = xpmt;
      Photon_xgpmt[ tid ] = xgpmt;
      Photon_PMT[ tid ] = pmt;
      Photon_row[ tid ] = row;
      Photon_col[ tid ] = col;
      Photon_mTID[ tid ] = mid; 
      Photon_origvol[ tid ] = origin_volume;
      Photon_nsteps[ tid ] = 1;

      Photon_otridx[ tid ] = SDtracks.otracklist[(*hits)[step]->GetOTrIdx()];
      Photon_ptridx[ tid ] = SDtracks.ptracklist[(*hits)[step]->GetPTrIdx()];
      Photon_sdtridx[ tid ] = SDtracks.sdtracklist[(*hits)[step]->GetSDTrIdx()][hits->GetSDname()];
      
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
	  PMT_tmin[ pmt ] = PMT_hittime[ pmt ];
	  PMT_tmax[ pmt ] = PMT_hittime[ pmt ];
	  PMT_mTID[ pmt ] = Photon_mTID[tid];
	  PMT_pos[ pmt ] = Photon_position[tid];
	  PMT_dir[ pmt ] = Photon_direction[tid];
	  PMT_vpos[ pmt ] = Photon_vertex[tid];
	  PMT_vdir[ pmt ] = Photon_vdirection[tid];
	  PMT_origvol[ pmt ] = Photon_origvol[tid];
	  PMT_x[ pmt ] = Photon_xpmt[tid];
	  PMT_xg[ pmt ] = Photon_xgpmt[tid];

	  OTrackIndices[ pmt ].insert( Photon_otridx[ tid ] );
	  PTrackIndices[ pmt ].insert( Photon_ptridx[ tid ] );
	  SDTrackIndices[ pmt ].insert( Photon_sdtridx[ tid ] );
	  
	  std::pair<set<int>::iterator,bool> newmother = mTIDs_unique.insert( Photon_mTID[tid] );

	  if( newmother.second ){
	    Nphe_mTID[ Photon_mTID[tid] ] = 1;
	  } else {
	    Nphe_mTID[ Photon_mTID[tid] ] += 1;
	  }
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

	  if( Photon_hittime[tid] < PMT_tmin[pmt] ) PMT_tmin[pmt] = Photon_hittime[tid];
	  if( Photon_hittime[tid] > PMT_tmax[pmt] ) PMT_tmax[pmt] = Photon_hittime[tid];

	  PMT_Numphotoelectrons[pmt] += 1;
	  
	  Photon_used[tid] = true;

	  OTrackIndices[ pmt ].insert( Photon_otridx[ tid ] );
	  PTrackIndices[ pmt ].insert( Photon_ptridx[ tid ] );
	  SDTrackIndices[ pmt ].insert( Photon_sdtridx[ tid ] );
	  
	  std::pair<set<int>::iterator,bool> newmother = mTIDs_unique.insert( Photon_mTID[tid] );

	  if( newmother.second ){
	    Nphe_mTID[ Photon_mTID[tid] ] = 1;
	  } else {
	    Nphe_mTID[ Photon_mTID[tid] ] += 1;
	  }
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
	richoutput.Time_avg.push_back( PMT_hittime[pmt]/_T_UNIT );
	richoutput.Time_rms.push_back( sqrt(fabs(PMT_hittime2[pmt] - pow(PMT_hittime[pmt],2) ) )/_T_UNIT );
	richoutput.Time_min.push_back( PMT_tmin[pmt]/_T_UNIT );
	richoutput.Time_max.push_back( PMT_tmax[pmt]/_T_UNIT );
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

	richoutput.xpmt.push_back( PMT_x[pmt].x()/_L_UNIT );
	richoutput.ypmt.push_back( PMT_x[pmt].y()/_L_UNIT );
	richoutput.zpmt.push_back( PMT_x[pmt].z()/_L_UNIT );
	richoutput.xgpmt.push_back( PMT_xg[pmt].x()/_L_UNIT );
	richoutput.ygpmt.push_back( PMT_xg[pmt].y()/_L_UNIT );
	richoutput.zgpmt.push_back( PMT_xg[pmt].z()/_L_UNIT );

	//If there are multiple Otracks contributing to this hit, choose the one with the highest total energy:
	G4double maxE = 0.0;
	G4bool firsttrack=true;
	int otridx_final=-1;
	for( set<int>::iterator iotrk=OTrackIndices[pmt].begin(); iotrk!=OTrackIndices[pmt].end(); ++iotrk ){
	  G4double Eotrack = SDtracks.oenergy[*iotrk];
	  if( firsttrack || Eotrack > maxE ){
	    otridx_final = *iotrk;
	    maxE = Eotrack;
	    firsttrack = false;
	  }
	}

	//Primary tracks:
	maxE = 0.0;
	firsttrack = true;
	int ptridx_final=-1;
	for( set<int>::iterator iptrk=PTrackIndices[pmt].begin(); iptrk!=PTrackIndices[pmt].end(); ++iptrk ){
	  G4double Eptrack = SDtracks.penergy[*iptrk];
	  if( firsttrack || Eptrack > maxE ){
	    ptridx_final = *iptrk;
	    maxE = Eptrack;
	    firsttrack = false;
	  }
	}

	//SD boundary crossing tracks:
	maxE = 0.0;
	firsttrack = true;
	int sdtridx_final=-1;
	for( set<int>::iterator isdtrk=SDTrackIndices[pmt].begin(); isdtrk!=SDTrackIndices[pmt].end(); ++isdtrk ){

	  int sdidxtemp = *isdtrk;

	  if( sdidxtemp >= 0 && sdidxtemp <SDtracks.sdenergy.size() ){
	    G4double Esdtrack = SDtracks.sdenergy[sdidxtemp];
	    if( firsttrack || Esdtrack > maxE ){
	      sdtridx_final = *isdtrk;
	      maxE = Esdtrack;
	      firsttrack = false;
	    }
	  }
	}

	richoutput.otridx.push_back( otridx_final );
	richoutput.ptridx.push_back( ptridx_final );
	richoutput.sdtridx.push_back( sdtridx_final );
	
      }
    }
    //PMTs_unique.clear();
  }
  
  G4TrajectoryContainer *tracklist = evt->GetTrajectoryContainer();

  if( !tracklist ) return;

  map<int,int> mtrackindex;

  set<int> MotherTrajectories;

  for(set<int>::iterator it=mTIDs_unique.begin(); it!=mTIDs_unique.end(); it++ ){ //this is a loop over all TIDs that produced photoelectrons:

    int PIDtemp = 0;
    int TIDtemp = *it;
    int MIDtemp = MotherTrackIDs[TIDtemp];
    int hitidx = -1;
    int nbouncetemp = 0;

    //G4cout << "Before history traversal, it = " << *it << G4endl;
    
    do {
      G4Trajectory *track = (G4Trajectory*) ( (*tracklist)[TrajectoryIndex[TIDtemp]] );
    
      G4ThreeVector pinitial = track->GetInitialMomentum();
      G4ThreeVector vinitial = track->GetPoint(0)->GetPosition();
    
      PIDtemp = track->GetPDGEncoding();
      MIDtemp = MotherTrackIDs[TIDtemp];
    
      std::pair<set<int>::iterator,bool> newtrajectory = MotherTrajectories.insert( TIDtemp ); //Whether this is a mother of a mother or not, this is the first time this TID is encountered!
      
      if( newtrajectory.second ){
	richoutput.ParticleHistory.npart++;
	richoutput.ParticleHistory.PID.push_back( PIDtemp );
	richoutput.ParticleHistory.MID.push_back( MIDtemp );
	richoutput.ParticleHistory.TID.push_back( TIDtemp );
	richoutput.ParticleHistory.nbounce.push_back( nbouncetemp );
	//determine whether this photon was detected as the first hit in any PMT:
	for(G4int ihit=0; ihit<richoutput.nhits_RICH; ihit++){
	  if( richoutput.mTrackNo[ihit] == *it ){ 
	    hitidx = ihit; 
	    break;
	  }
	}
	richoutput.ParticleHistory.hitindex.push_back( hitidx );
	richoutput.ParticleHistory.vx.push_back( vinitial.x()/_L_UNIT );
	richoutput.ParticleHistory.vy.push_back( vinitial.y()/_L_UNIT );
	richoutput.ParticleHistory.vz.push_back( vinitial.z()/_L_UNIT );
	richoutput.ParticleHistory.px.push_back( pinitial.x()/_E_UNIT );
	richoutput.ParticleHistory.py.push_back( pinitial.y()/_E_UNIT );
	richoutput.ParticleHistory.pz.push_back( pinitial.z()/_E_UNIT );
	if( mTIDs_unique.find( TIDtemp ) != mTIDs_unique.end() ){ //if this track is in the set of all tracks producing photodetections, record the number of photoelectrons associated:
	  richoutput.Nphe_part.push_back( Nphe_mTID[ *it ] );
	} else {
	  richoutput.Nphe_part.push_back( 0 );
	}
      }
      
      TIDtemp = MIDtemp;
      nbouncetemp++;
    } while ( MIDtemp != 0 );
    // G4cout << "After history traversal, it = " << *it << G4endl;
  }

  // for(G4int ihit=0; ihit<richoutput.nhits_RICH; ihit++){
  //   set<int>::iterator pos = mTIDs_unique.find( richoutput.mTrackNo[ihit] );
  //   if( pos != mTIDs_unique.end() ){
  //     richoutput.mTrackNo[ihit] = mtrackindex[ *pos ];
  //   } else {
  //     richoutput.mTrackNo[ihit] = -1; //This should never happen!
  //   }
  // }
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
      TrajectoryIndex[ TrackID ] = i; //Index of this track in the TrajectoryContainer array!
      MotherTrackIDs[ TrackID ] = MotherID; //This is the mother track ID number corresponding to Track ID number 
      
    }
  }
}

void G4SBSEventAction::FillTrackData( G4SBSGEMoutput gemdata, G4SBSTrackerOutput &Toutput ){
  //Note: gemdata have already been normalized to the correct units (meters, ns, GeV) and are already expressed in TRANSPORT coordinates:
  //Also note that gemdata.x and gemdata.y have already been smeared by coordinate resolution!

  set<int> TrackTIDs_unique; //Track numbers of unique tracks causing GEM hits in this event:
  
  map<int, vector<int> > HitList; //list of indices in the hit array of hits associated with a given track:

  G4int nhits = gemdata.nhits_GEM;

  //First, map all hits to track IDs:
  for(int i=0; i<nhits; i++){
    int tid = gemdata.trid[i];   
   
    double edep = gemdata.edep[i];

    if( edep > 0.0 ){
      TrackTIDs_unique.insert( tid );
      
      HitList[tid].push_back( i );
    }
  }

  int nplanes_min = 3; //Minimum number of valid hits to define a track:

  //For the "true" track, we simply take the average of all points for x, y, xp, yp 
  //Fit procedure: minimize chi^2 defined as sum_i=1,N (xi - (x0 + xp*z))^2/sigmax^2 + (yi - (y0+yp*z))^2/sigmay^2:

  for(set<int>::iterator trk=TrackTIDs_unique.begin(); trk != TrackTIDs_unique.end(); trk++ ){
    int track   = *trk;
    //Define sums to keep track of: 
    int nhittrk = 0;
    int nplanetrk = 0;
      
    set<int> GEMIDs_unique; //List of unique GEM plane ids on this track:

    double x0avg = 0.0, y0avg = 0.0, xpavg = 0.0, ypavg = 0.0;
    double x0avg2 = 0.0, y0avg2 = 0.0, xpavg2 = 0.0, ypavg2 = 0.0;
    double tavg = 0.0, tavg2 = 0.0;
    double pavg = 0.0;

    double polx_avg = 0.0, poly_avg = 0.0, polz_avg = 0.0; //compute average track polarization for each track.
    
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
      
    nhittrk = HitList[track].size();

    if( nhittrk >= nplanes_min ){
	
      vector<double> xsmear,ysmear,xtrue,ytrue,ztrue,hittime,beta;
      
      int PID = 0, MID = -1;
      
      double sigma = fGEMres/_L_UNIT;

      for(int idx = 0; idx<nhittrk; idx++ ){ //First pass: determine number of unique GEM planes and compute average slopes of "true track". Also compute the sums for the linear fit:
	int hit = HitList[track][idx];
	
	if( idx == 0 ) {
	  PID = gemdata.pid[hit];
	  MID = gemdata.mid[hit];
	}
	
	GEMIDs_unique.insert( gemdata.plane[hit] );
	
	pavg += gemdata.p[hit]/double(nhittrk);

	polx_avg += gemdata.polx[hit]/double(nhittrk);
	poly_avg += gemdata.poly[hit]/double(nhittrk);
	polz_avg += gemdata.polz[hit]/double(nhittrk);
	
	xtrue.push_back( gemdata.tx[hit] );
	ytrue.push_back( gemdata.ty[hit] );

	xsmear.push_back( gemdata.x[hit] );
	ysmear.push_back( gemdata.y[hit] );
	  
	ztrue.push_back( gemdata.z[hit] );
	 
	hittime.push_back( gemdata.t[hit] );
	beta.push_back( gemdata.beta[hit] );
	  
	//Indices of fit parameters in matrix are: 0 = x0, 1 = xp, 2 = y0, 3 = yp:
	M(0,0) += pow(1.0/sigma,2); 
	M(0,1) += pow(1.0/sigma,2) * ztrue[idx];
	M(0,2) += 0.0; //No cross term between x0 and y0
	M(0,3) += 0.0; //No cross term between x0 and yp
	M(1,0) += pow(1.0/sigma,2) * ztrue[idx];
	M(1,1) += pow(ztrue[idx]/sigma,2);
	M(1,2) += 0.0;
	M(1,3) += 0.0;
	M(2,0) += 0.0;
	M(2,1) += 0.0;
	M(2,2) += pow(1.0/sigma,2);
	M(2,3) += pow(1.0/sigma,2) * ztrue[idx];
	M(3,0) += 0.0;
	M(3,1) += 0.0;
	M(3,2) += pow(1.0/sigma,2) * ztrue[idx];
	M(3,3) += pow(ztrue[idx]/sigma,2);
	  
	b(0) += xsmear[idx] / pow(sigma,2);
	b(1) += ztrue[idx] * xsmear[idx] / pow(sigma,2);
	b(2) += ysmear[idx] / pow(sigma,2);
	b(3) += ztrue[idx] * ysmear[idx] / pow(sigma,2);
	  
	btrue(0) += xtrue[idx] / pow(sigma,2);
	btrue(1) += ztrue[idx] * xtrue[idx] / pow(sigma,2);
	btrue(2) += ytrue[idx] / pow(sigma,2);
	btrue(3) += ztrue[idx] * ytrue[idx] / pow(sigma,2);

      }
	
      if( GEMIDs_unique.size() >= nplanes_min ){

	TMatrixD Minv = M.Invert();
	TVectorD FitTrack = Minv * b;
	TVectorD TrueTrack = Minv * btrue;

	Toutput.ntracks++;
	
	Toutput.TrackTID.push_back(track);
	Toutput.TrackPID.push_back( PID );
	Toutput.TrackMID.push_back( MID );
	Toutput.NumHits.push_back( nhittrk );
	Toutput.NumPlanes.push_back( GEMIDs_unique.size() );
	  
	//Everything should already be expressed in the desired units in gemdata:
	Toutput.TrackX.push_back( TrueTrack(0) );
	Toutput.TrackXp.push_back( TrueTrack(1) );
	Toutput.TrackY.push_back( TrueTrack(2) );
	Toutput.TrackYp.push_back( TrueTrack(3) );

	Toutput.TrackXfit.push_back( FitTrack(0) );
	Toutput.TrackXpfit.push_back( FitTrack(1) );
	Toutput.TrackYfit.push_back( FitTrack(2) );
	Toutput.TrackYpfit.push_back( FitTrack(3) );

	int ndf = 0;
	double chi2 = 0.0, chi2_true = 0.0;
	//Compute chi^2 and focal plane times projected to local origin:
	for(int idx = 0; idx<nhittrk; idx++){

	  //tfp is hit time corrected for time of flight: ztrue has units of meters, while hittime has units of ns
	  //convert z to meters, then the tof correction term will have units of seconds, so we need to divide by _T_UNIT to get ns!
	  double tfp = hittime[idx] - ztrue[idx]*sqrt( 1.0 + pow(TrueTrack(1),2)+pow(TrueTrack(3),2) ) / (beta[idx]*c_light)*(_L_UNIT/_T_UNIT);
	    
	  chi2 += pow( (xsmear[idx] - (FitTrack(0) + FitTrack(1)*ztrue[idx] ) )/sigma, 2 );
	  chi2 += pow( (ysmear[idx] - (FitTrack(2) + FitTrack(3)*ztrue[idx] ) )/sigma, 2 );

	  chi2_true += pow( (xtrue[idx] - (TrueTrack(0) + TrueTrack(1)*ztrue[idx] ) )/sigma, 2 );
	  chi2_true += pow( (ytrue[idx] - (TrueTrack(2) + TrueTrack(3)*ztrue[idx] ) )/sigma, 2 );

	  tavg += tfp/double(nhittrk);

	  ndf += 2;
	}

	Toutput.TrackSx.push_back( polx_avg );
	Toutput.TrackSy.push_back( poly_avg );
	Toutput.TrackSz.push_back( polz_avg );
	
	Toutput.Chi2fit.push_back( chi2 );
	Toutput.Chi2true.push_back( chi2_true );
	Toutput.NDF.push_back( ndf - 4 );
	Toutput.TrackT.push_back( tavg ); //tavg already in ns!

	Toutput.TrackP.push_back( pavg ); //pavg already in GeV!

	Toutput.otridx.push_back( gemdata.otridx[HitList[track][0]] );
	Toutput.ptridx.push_back( gemdata.ptridx[HitList[track][0]] );
	Toutput.sdtridx.push_back( gemdata.sdtridx[HitList[track][0]] );
      }
    }
  }
}

void G4SBSEventAction::FillBDData(const G4Event *evt,G4SBSBDHitsCollection *hc,G4SBSBDoutput &out){
   // fill the G4SBSBDoutput class with hit data

   out.Clear();   // clear previous event data 

   int nstep=0,trackID=0,bdID=0;
   double w=0;
   std::map< int,std::set<int> > tracks_layers;           // key = BeamDiffuser layer ID, value = set of unique tracks with edep in layer 
   std::map< int,std::map<int,int> > nsteps_track_layer;  // number of steps by track/layer  
   std::map< int,std::map<int,int> > pid;                 // particle type
   std::map< int,std::map<int,int> > mid;                 // material/medium type (?)   
   std::map< int,std::map<int,double> > x,y,z,t,p,edep;
   std::map< int,std::map<int,double> > xg,yg,zg,beta;

   bool debug=false;
   char msg[200]; 

   int NHits = (int)hc->entries();
   if(debug){
      if(NHits!=0){
	 std::cout << "[G4SBSEventAction::FillBDData]: Found " << NHits << " hits!" << std::endl;
	 std::cout << "[G4SBSEventAction::FillBDData]: Printing first five hits: " << std::endl;
      }
   }

   // loop over all "hits" (i.e., individual tracking steps)
   for(int i=0;i<NHits;i++){
      // get track ID and BeamDiffuser plane ID 
      trackID = (*hc)[i]->GetTrackID();
      bdID    = (*hc)[i]->GetPlane();
      // now we examine the track
      std::pair<std::set<int>::iterator, bool> track = tracks_layers[bdID].insert(trackID);
      if( track.second ){
         // new track in this layer, first step
         nsteps_track_layer[bdID][trackID] = 1;
         // time of hit 
         t[bdID][trackID]    = (*hc)[i]->GetHitTime();
         // positional data (local coordinates of detector)  
         x[bdID][trackID]    = (*hc)[i]->GetPos().x();
         y[bdID][trackID]    = (*hc)[i]->GetPos().y();
         z[bdID][trackID]    = (*hc)[i]->GetPos().z();
         // positional data (global or lab coordinates) 
         xg[bdID][trackID]   = (*hc)[i]->GetLabPos().x();
         yg[bdID][trackID]   = (*hc)[i]->GetLabPos().y();
         zg[bdID][trackID]   = (*hc)[i]->GetLabPos().z();
         // energy and momentum 
         edep[bdID][trackID] = (*hc)[i]->GetEdep();
         p[bdID][trackID]    = (*hc)[i]->GetMom();  // momentum (magnitude) 
         beta[bdID][trackID] = (*hc)[i]->GetBeta();
         // Particle and material info  
         pid[bdID][trackID]  = (*hc)[i]->GetPID();
         mid[bdID][trackID]  = (*hc)[i]->GetMID();
      }else{
         // existing track in this layer, additional step; increment sums and averages
         nstep = nsteps_track_layer[bdID][trackID];
         w     = (double)nstep/( (double)(nstep+1) );
         // the coordinates below represent local BD hit coordinates 
         x[bdID][trackID] = x[bdID][trackID]*w +( (*hc)[i]->GetPos().x() )*(1.0-w);
         y[bdID][trackID] = y[bdID][trackID]*w +( (*hc)[i]->GetPos().y() )*(1.0-w);
         z[bdID][trackID] = z[bdID][trackID]*w +( (*hc)[i]->GetPos().z() )*(1.0-w);
         // the global coordinates 
         xg[bdID][trackID] = xg[bdID][trackID]*w +( (*hc)[i]->GetLabPos().x() )*(1.0-w);
         yg[bdID][trackID] = yg[bdID][trackID]*w +( (*hc)[i]->GetLabPos().y() )*(1.0-w);
         zg[bdID][trackID] = zg[bdID][trackID]*w +( (*hc)[i]->GetLabPos().z() )*(1.0-w);
         // for edep, we do the sum:
         edep[bdID][trackID] += (*hc)[i]->GetEdep();
         // increment 
         nsteps_track_layer[bdID][trackID]++;
      }
      if(debug){
	 sprintf(msg,"hit %04d, track %04d, plane %02d, edep = %.3lf keV",i+1,trackID,bdID,(*hc)[i]->GetEdep()/CLHEP::keV);
	 if((i+1)<5) std::cout << msg << std::endl;  // print first 5 hits 
      }
   }

   bdID = 0;
   trackID = 0;

   // G4TrajectoryContainer *trajectorylist = evt->GetTrajectoryContainer(); // for particle history information

   // for particle history details. mimics what is done for the GEMs 
   // int MIDtemp=0,TIDtemp=0,PIDtemp=0,hitidx=0,nbouncetemp=0; 
   // std::set<int> TIDs_unique; //all unique track IDs involved in BD hits in this event (for filling particle history tree)

   // now accumulate data into output class 
   for(std::map<int,std::set<int> >::iterator hit=tracks_layers.begin(); hit!=tracks_layers.end(); hit++ ){
      std::set<int> tracklist = hit->second;
      bdID = hit->first;
      for(std::set<int>::iterator track=tracklist.begin(); track!=tracklist.end(); track++ ){
         trackID = *track;
         out.plane.push_back( bdID );
         out.t.push_back( t[bdID][trackID]/_T_UNIT );
         // coordinates in detector system
         out.x.push_back( (-y[bdID][trackID])/_L_UNIT );
         out.y.push_back( (x[bdID][trackID])/_L_UNIT );
         out.z.push_back( z[bdID][trackID]/_L_UNIT );
         // coordinates in the hall 
         out.xg.push_back( xg[bdID][trackID]/_L_UNIT );
         out.yg.push_back( yg[bdID][trackID]/_L_UNIT );
         out.zg.push_back( zg[bdID][trackID]/_L_UNIT );
         out.trid.push_back( trackID );
         out.pid.push_back( pid[bdID][trackID] );
         out.mid.push_back( mid[bdID][trackID] );
         out.p.push_back( p[bdID][trackID]/_E_UNIT );
         out.beta.push_back( beta[bdID][trackID]/_E_UNIT );
         out.edep.push_back( edep[bdID][trackID]/_E_UNIT );
         // if( trajectorylist ){ 
         //    // fill Particle History, starting with the particle itself 
         //    // and working all the way back to primary particles:
         //    MIDtemp = mid[gemID][trackID];
         //    TIDtemp = trackID;
         //    PIDtemp = pid[gemID][trackID];
         //    hitidx = out.nhits;
         //    nbouncetemp = 0;
         //    do {
         //       G4Trajectory *trajectory = (G4Trajectory*) (*trajectorylist)[TrajectoryIndex[TIDtemp]];
         //       PIDtemp = trajectory->GetPDGEncoding();
         //       MIDtemp = MotherTrackIDs[TIDtemp];
         //       std::pair<set<int>::iterator, bool > newtrajectory = TIDs_unique.insert( TIDtemp );
         //       if( newtrajectory.second ){ 
         //          // this trajectory does not yet exist in the 
         //          // particle history of this detector for this event. Add it:
         //          out.ParticleHistory.PID.push_back( PIDtemp );
         //          out.ParticleHistory.MID.push_back( MIDtemp );
         //          out.ParticleHistory.TID.push_back( TIDtemp );
         //          // for hitindex: of course, this means that if a trajectory is 
         //          // involved in multiple hits in this detector, this variable 
         //          // will point to the first hit encountered only!
         //          out.ParticleHistory.hitindex.push_back( hitidx ); 
         //          out.ParticleHistory.nbounce.push_back( nbouncetemp );
         //          out.ParticleHistory.vx.push_back( (trajectory->GetPoint(0)->GetPosition() ).x()/_L_UNIT );
         //          out.ParticleHistory.vy.push_back( (trajectory->GetPoint(0)->GetPosition() ).y()/_L_UNIT );
         //          out.ParticleHistory.vz.push_back( (trajectory->GetPoint(0)->GetPosition() ).z()/_L_UNIT );
         //          out.ParticleHistory.px.push_back( (trajectory->GetInitialMomentum() ).x()/_E_UNIT );
         //          out.ParticleHistory.py.push_back( (trajectory->GetInitialMomentum() ).y()/_E_UNIT );
         //          out.ParticleHistory.pz.push_back( (trajectory->GetInitialMomentum() ).z()/_E_UNIT );
         //          out.ParticleHistory.npart++;
         //       }
         //       TIDtemp = MIDtemp;
         //       nbouncetemp++;
         //    } while( MIDtemp!=0 );
         // }
      }
      out.nhits_BD++;
   }
}

void G4SBSEventAction::FillICData(const G4Event *evt,G4SBSICHitsCollection *hc,G4SBSICoutput &out){
   // fill the G4SBSICoutput class with hit data

   out.Clear();   // clear previous event data 

   int nstep=0,trackID=0;
   double w=0;
   std::set<int>  tracks_layers;          // key = BeamDiffuser layer ID, value = set of unique tracks with edep in layer 
   std::map<int,int> nsteps_track_layer;  // number of steps by track/layer  
   std::map<int,int> nsteps_track;        // number of steps by track/layer  
   std::map<int,int> pid;                 // particle type
   std::map<int,int> mid;                 // material/medium type (?)   
   std::map<int,double> x,y,z,t,p,edep;
   std::map<int,double> xg,yg,zg,beta;

   bool debug=false;
   char msg[200]; 

   int NHits = (int)hc->entries();
   if(debug){
      if(NHits!=0){
	 std::cout << "[G4SBSEventAction::FillICData]: Found " << NHits << " hits!" << std::endl;
	 std::cout << "[G4SBSEventAction::FillICData]: Printing first five hits: " << std::endl;
      }
   }

   // loop over all "hits" (i.e., individual tracking steps)
   for(int i=0;i<NHits;i++){
      // get track ID  
      trackID = (*hc)[i]->GetTrackID();
      // now we examine the track
      std::pair<std::set<int>::iterator, bool> track = tracks_layers.insert(trackID);
      if( track.second ){
         // new track in this layer, first step
         nsteps_track_layer[trackID] = 1;
         // time of hit 
         t[trackID]    = (*hc)[i]->GetHitTime();
         // positional data (local coordinates of detector)  
         x[trackID]    = (*hc)[i]->GetPos().x();
         y[trackID]    = (*hc)[i]->GetPos().y();
         z[trackID]    = (*hc)[i]->GetPos().z();
         // positional data (global or lab coordinates) 
         xg[trackID]   = (*hc)[i]->GetLabPos().x();
         yg[trackID]   = (*hc)[i]->GetLabPos().y();
         zg[trackID]   = (*hc)[i]->GetLabPos().z();
         // energy and momentum 
         edep[trackID] = (*hc)[i]->GetEdep();
         p[trackID]    = (*hc)[i]->GetMomentumMag();  // momentum (magnitude) at pre-step 
         beta[trackID] = (*hc)[i]->GetBeta();
         // Particle and material info  
         pid[trackID]  = (*hc)[i]->GetPID();
         mid[trackID]  = (*hc)[i]->GetMID();
      }else{
         // existing track in this layer, additional step; increment sums and averages
         nstep = nsteps_track_layer[trackID];
         w     = (double)nstep/( (double)(nstep+1) );
         // the coordinates below represent local IC hit coordinates 
         x[trackID] = x[trackID]*w +( (*hc)[i]->GetPos().x() )*(1.0-w);
         y[trackID] = y[trackID]*w +( (*hc)[i]->GetPos().y() )*(1.0-w);
         z[trackID] = z[trackID]*w +( (*hc)[i]->GetPos().z() )*(1.0-w);
         // the global coordinates 
         xg[trackID] = xg[trackID]*w +( (*hc)[i]->GetLabPos().x() )*(1.0-w);
         yg[trackID] = yg[trackID]*w +( (*hc)[i]->GetLabPos().y() )*(1.0-w);
         zg[trackID] = zg[trackID]*w +( (*hc)[i]->GetLabPos().z() )*(1.0-w);
         // for edep, we do the sum:
         edep[trackID] += (*hc)[i]->GetEdep();
         // increment 
         nsteps_track_layer[trackID]++;
      }
      if(debug){
	 sprintf(msg,"hit %04d, track %04d, edep = %.3lf keV",i+1,trackID,(*hc)[i]->GetEdep()/CLHEP::keV);
	 if((i+1)<5) std::cout << msg << std::endl;  // print first 5 hits 
      }
   }

   trackID = 0;

   G4TrajectoryContainer *trajectorylist = evt->GetTrajectoryContainer(); // for particle history information

   // for particle history details. mimics what is done for the GEMs 
   int MIDtemp=0,TIDtemp=0,PIDtemp=0,hitidx=0,nbouncetemp=0; 
   std::set<int> TIDs_unique; // all unique track IDs involved in IC hits in this event (for filling particle history tree)

   // now accumulate data into output class 
   for(std::set<int>::iterator tid=tracks_layers.begin(); tid!=tracks_layers.end(); tid++ ){
      // std::set<int> tracklist = hit->second;
      trackID = *tid;
      // for(std::set<int>::iterator track=tracklist.begin(); track!=tracklist.end(); track++ ){
         // trackID = *track;
         out.t.push_back( t[trackID]/_T_UNIT );
         // coordinates in detector system
         out.x.push_back( (-y[trackID])/_L_UNIT );
         out.y.push_back( (x[trackID])/_L_UNIT );
         out.z.push_back( z[trackID]/_L_UNIT );
         // coordinates in the hall 
         out.xg.push_back( xg[trackID]/_L_UNIT );
         out.yg.push_back( yg[trackID]/_L_UNIT );
         out.zg.push_back( zg[trackID]/_L_UNIT );
         out.trid.push_back( trackID );
         out.pid.push_back( pid[trackID] );
         out.mid.push_back( mid[trackID] );
         out.p.push_back( p[trackID]/_E_UNIT );
         out.beta.push_back( beta[trackID]/_E_UNIT );
         out.edep.push_back( edep[trackID]/_E_UNIT );
         if( trajectorylist ){ 
            // fill Particle History, starting with the particle itself 
            // and working all the way back to primary particles:
            MIDtemp = mid[trackID];
            TIDtemp = trackID;
            PIDtemp = pid[trackID];
            hitidx  = out.nhits_IC;
            nbouncetemp = 0;
            do {
               G4Trajectory *trajectory = (G4Trajectory*) (*trajectorylist)[TrajectoryIndex[TIDtemp]];
               PIDtemp = trajectory->GetPDGEncoding();
               MIDtemp = MotherTrackIDs[TIDtemp];
               std::pair<set<int>::iterator, bool > newtrajectory = TIDs_unique.insert( TIDtemp );
               if( newtrajectory.second ){ 
                  // this trajectory does not yet exist in the 
                  // particle history of this detector for this event. Add it:
                  out.ParticleHistory.PID.push_back( PIDtemp );
                  out.ParticleHistory.MID.push_back( MIDtemp );
                  out.ParticleHistory.TID.push_back( TIDtemp );
                  // for hitindex: of course, this means that if a trajectory is 
                  // involved in multiple hits in this detector, this variable 
                  // will point to the first hit encountered only!
                  out.ParticleHistory.hitindex.push_back( hitidx ); 
                  out.ParticleHistory.nbounce.push_back( nbouncetemp );
                  out.ParticleHistory.vx.push_back( (trajectory->GetPoint(0)->GetPosition() ).x()/_L_UNIT );
                  out.ParticleHistory.vy.push_back( (trajectory->GetPoint(0)->GetPosition() ).y()/_L_UNIT );
                  out.ParticleHistory.vz.push_back( (trajectory->GetPoint(0)->GetPosition() ).z()/_L_UNIT );
                  out.ParticleHistory.px.push_back( (trajectory->GetInitialMomentum() ).x()/_E_UNIT );
                  out.ParticleHistory.py.push_back( (trajectory->GetInitialMomentum() ).y()/_E_UNIT );
                  out.ParticleHistory.pz.push_back( (trajectory->GetInitialMomentum() ).z()/_E_UNIT );
                  out.ParticleHistory.npart++;
               }
               TIDtemp = MIDtemp;
               nbouncetemp++;
            } while( MIDtemp!=0 );
         }
      // }
      out.nhits_IC++;
   }
}

void G4SBSEventAction::FillGEnTargetData(const G4Event *evt,G4SBSTargetHitsCollection *hc,G4SBSTargetoutput &out){
   // fill the G4SBSTargetoutput class with hit data

   out.Clear();   // clear previous event data 

   int nstep=0,trackID=0;
   double w=0;
   std::set<int>  tracks_layers;          // key = BeamDiffuser layer ID, value = set of unique tracks with edep in layer 
   std::map<int,int> nsteps_track_layer;  // number of steps by track/layer  
   std::map<int,int> nsteps_track;        // number of steps by track/layer  
   std::map<int,int> pid;                 // particle type
   std::map<int,int> mid;                 // material/medium type (?)   
   std::map<int,double> beta,edep,p,t,trackLen;
   // std::map<int,double> x,y,z,xg,yg,zg;

   bool debug=false;
   char msg[200]; 

   int NHits = (int)hc->entries();
   if(debug){
      if(NHits!=0){
	 std::cout << "[G4SBSEventAction::FillGEnTargetData]: Found " << NHits << " hits!" << std::endl;
	 std::cout << "[G4SBSEventAction::FillGEnTargetData]: Printing first five hits: " << std::endl;
      }
   }

   // loop over all "hits" (i.e., individual tracking steps)
   for(int i=0;i<NHits;i++){
      // get track ID  
      trackID = (*hc)[i]->GetTrackID();
      // now we examine the track
      std::pair<std::set<int>::iterator, bool> track = tracks_layers.insert(trackID);
      if( track.second ){
         // new track in this layer, first step
         nsteps_track_layer[trackID] = 1;
         // time of hit 
         t[trackID]    = (*hc)[i]->GetHitTime();
         // // positional data (local coordinates of detector)  
         // x[trackID]    = (*hc)[i]->GetPos().x();
         // y[trackID]    = (*hc)[i]->GetPos().y();
         // z[trackID]    = (*hc)[i]->GetPos().z();
         // // positional data (global or lab coordinates) 
         // xg[trackID]   = (*hc)[i]->GetLabPos().x();
         // yg[trackID]   = (*hc)[i]->GetLabPos().y();
         // zg[trackID]   = (*hc)[i]->GetLabPos().z();
         // energy and momentum 
         edep[trackID] = (*hc)[i]->GetEdep();
         p[trackID]    = (*hc)[i]->GetMomentumMag();  // momentum (magnitude) at pre-step 
         beta[trackID] = (*hc)[i]->GetBeta();
         trackLen[trackID] = (*hc)[i]->GetTrackLength(); 
         // Particle and material info  
         pid[trackID]  = (*hc)[i]->GetPID();
         mid[trackID]  = (*hc)[i]->GetMID();
      }else{
         // existing track in this layer, additional step; increment sums and averages
         nstep = nsteps_track_layer[trackID];
         w     = (double)nstep/( (double)(nstep+1) );
         // // the coordinates below represent local IC hit coordinates 
         // x[trackID] = x[trackID]*w +( (*hc)[i]->GetPos().x() )*(1.0-w);
         // y[trackID] = y[trackID]*w +( (*hc)[i]->GetPos().y() )*(1.0-w);
         // z[trackID] = z[trackID]*w +( (*hc)[i]->GetPos().z() )*(1.0-w);
         // // the global coordinates 
         // xg[trackID] = xg[trackID]*w +( (*hc)[i]->GetLabPos().x() )*(1.0-w);
         // yg[trackID] = yg[trackID]*w +( (*hc)[i]->GetLabPos().y() )*(1.0-w);
         // zg[trackID] = zg[trackID]*w +( (*hc)[i]->GetLabPos().z() )*(1.0-w);
         trackLen[trackID] += (*hc)[i]->GetTrackLength(); 
         // for edep, we do the sum:
         edep[trackID] += (*hc)[i]->GetEdep();
         // increment 
         nsteps_track_layer[trackID]++;
      }
      if(debug){
	 sprintf(msg,"hit %04d, track %04d, edep = %.3lf keV",i+1,trackID,(*hc)[i]->GetEdep()/CLHEP::keV);
	 if((i+1)<5) std::cout << msg << std::endl;  // print first 5 hits 
      }
   }

   trackID = 0;

   G4TrajectoryContainer *trajectorylist = evt->GetTrajectoryContainer(); // for particle history information

   // for particle history details. mimics what is done for the GEMs 
   int MIDtemp=0,TIDtemp=0,PIDtemp=0,hitidx=0,nbouncetemp=0; 
   std::set<int> TIDs_unique; // all unique track IDs involved in IC hits in this event (for filling particle history tree)

   // now accumulate data into output class 
   for(std::set<int>::iterator tid=tracks_layers.begin(); tid!=tracks_layers.end(); tid++ ){
      // std::set<int> tracklist = hit->second;
      trackID = *tid;
      // for(std::set<int>::iterator track=tracklist.begin(); track!=tracklist.end(); track++ ){
         // trackID = *track;
         out.t.push_back( t[trackID]/_T_UNIT );
         // // coordinates in detector system
         // out.x.push_back( (-y[trackID])/_L_UNIT );
         // out.y.push_back( (x[trackID])/_L_UNIT );
         // out.z.push_back( z[trackID]/_L_UNIT );
         // // coordinates in the hall 
         // out.xg.push_back( xg[trackID]/_L_UNIT );
         // out.yg.push_back( yg[trackID]/_L_UNIT );
         // out.zg.push_back( zg[trackID]/_L_UNIT );
         out.trid.push_back( trackID );
         out.pid.push_back( pid[trackID] );
         out.mid.push_back( mid[trackID] );
         out.p.push_back( p[trackID]/_E_UNIT );
         out.beta.push_back( beta[trackID]/_E_UNIT );
         out.edep.push_back( edep[trackID]/_E_UNIT );
         out.trackLength.push_back( trackLen[trackID]/_L_UNIT ); 
         if( trajectorylist ){ 
            // fill Particle History, starting with the particle itself 
            // and working all the way back to primary particles:
            MIDtemp = mid[trackID];
            TIDtemp = trackID;
            PIDtemp = pid[trackID];
            hitidx  = out.nhits_Target;
            nbouncetemp = 0;
            do {
               G4Trajectory *trajectory = (G4Trajectory*) (*trajectorylist)[TrajectoryIndex[TIDtemp]];
               PIDtemp = trajectory->GetPDGEncoding();
               MIDtemp = MotherTrackIDs[TIDtemp];
               std::pair<set<int>::iterator, bool > newtrajectory = TIDs_unique.insert( TIDtemp );
               if( newtrajectory.second ){ 
                  // this trajectory does not yet exist in the 
                  // particle history of this detector for this event. Add it:
                  out.ParticleHistory.PID.push_back( PIDtemp );
                  out.ParticleHistory.MID.push_back( MIDtemp );
                  out.ParticleHistory.TID.push_back( TIDtemp );
                  // for hitindex: of course, this means that if a trajectory is 
                  // involved in multiple hits in this detector, this variable 
                  // will point to the first hit encountered only!
                  out.ParticleHistory.hitindex.push_back( hitidx ); 
                  out.ParticleHistory.nbounce.push_back( nbouncetemp );
                  // out.ParticleHistory.vx.push_back( (trajectory->GetPoint(0)->GetPosition() ).x()/_L_UNIT );
                  // out.ParticleHistory.vy.push_back( (trajectory->GetPoint(0)->GetPosition() ).y()/_L_UNIT );
                  // out.ParticleHistory.vz.push_back( (trajectory->GetPoint(0)->GetPosition() ).z()/_L_UNIT );
                  out.ParticleHistory.px.push_back( (trajectory->GetInitialMomentum() ).x()/_E_UNIT );
                  out.ParticleHistory.py.push_back( (trajectory->GetInitialMomentum() ).y()/_E_UNIT );
                  out.ParticleHistory.pz.push_back( (trajectory->GetInitialMomentum() ).z()/_E_UNIT );
                  out.ParticleHistory.npart++;
               }
               TIDtemp = MIDtemp;
               nbouncetemp++;
            } while( MIDtemp!=0 );
         }
      // }
      out.nhits_Target++;
   }
}

