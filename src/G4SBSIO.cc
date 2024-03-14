#include "TObjString.h"
#include "TString.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TClonesArray.h>

#include "G4SBSGlobalField.hh"
#include "G4SBSRun.hh"
#include "G4SBSIO.hh"
#include "G4SBSCalSD.hh"
#include "G4SBSECalSD.hh"
#include "G4SBSEArmBuilder.hh"
#include "G4SBSHArmBuilder.hh"
//#include "G4SDManager.hh"
#include <assert.h>
#include "sbstypes.hh"

G4SBSIO::G4SBSIO(){
  fTree = NULL;
  //InitializeTree(); //We want experiment-dependent ROOT tree! Don't invoke until after fdetcon->ConstructAll() has been invoked!
  // Default filename
  strcpy(fFilename, "g4sbsout.root");
  fFile = NULL;

  //Sometimes these default parameters don't match what was actually used in simulation if certain commands
  //aren't invoked. Moreover, these default values don't match the default values in
  //G4SBSDetectorConstruction, G4SBSEarmBuilder, G4SBSHArmBuilder, etc.
  //Therefore, if commands aren't invoked, these defaults can be wrong/misleading.
  //In principle the best way to make sure there is one consistent set of default values is to grab this
  //information from the G4SBSDetectorConstruction whenever it changes!
  //Easiest (but not necessarily best?) way is to store a pointer to G4SBSIO as a data member of G4SBSDetectorConstruction so they can
  //directly talk to each other?
  gendata.Ebeam = 2.2;
  gendata.Ibeam = 1.0;
  gendata.thbb = 40.0*CLHEP::deg;
  gendata.dbb = 1.5;
  gendata.thsbs = 39.4*CLHEP::deg;
  gendata.dhcal = 17.0;
  gendata.voffhcal = 0.0;
  gendata.dsbs = 1.6;
  gendata.drich = 4.6;
  gendata.dsbstrkr = 4.3;

  KeepPartCALflags.clear();
  KeepHistoryflags.clear();
  //  KeepSDTracks.clear();
  //fKeepSDtracks = true; //by default.
  
  GEMdata.clear();
  CALdata.clear();
  richdata.clear();
  trackdata.clear();
  ecaldata.clear();
  sdtrackdata.clear();
  BDdata.clear(); 
  ICdata.clear(); 
  genTgtGCdata.clear(); 
  genTgtCUdata.clear(); 
  genTgtALdata.clear(); 
  genTgt3HEdata.clear(); 

  //Set SD track data recording to OFF by default:
  fKeepAllSDtracks = false;
  fKeepSDtracks.clear();

  //Set SD track data recording to OFF by default:
  fKeepAllPulseShape = false;
  fKeepPulseShape.clear();
  
  Esum_histograms = NULL;
  PulseShape_histograms = NULL;
    
  fUsePythia = false;
  
  fUseSIMC = false;

  fKineType = G4SBS::kElastic;
  
  fNhistograms = 0;

  fUsingScintillation = false;
  fUsingCerenkov = false;

  fWritePortableFieldMaps = false;
}

G4SBSIO::~G4SBSIO(){
  if( fTree ){delete fTree;}
  fTree = NULL;

  if( fFile ){delete fFile;}
  fFile = NULL;

  if( Esum_histograms ){ delete Esum_histograms; }
  if( PulseShape_histograms ){ delete PulseShape_histograms; }
}

void G4SBSIO::SetGEMData( G4String SDname, G4SBSGEMoutput gd ){
  GEMdata[SDname] = gd;
}

void G4SBSIO::SetTrackData( G4String SDname, G4SBSTrackerOutput td ){
  trackdata[SDname] = td;
}

void G4SBSIO::SetCalData( G4String SDname, G4SBSCALoutput cd ){
  CALdata[SDname] = cd;
}

void G4SBSIO::SetRICHData( G4String SDname, G4SBSRICHoutput rd ){
  richdata[SDname] = rd;
}

void G4SBSIO::SetECalData( G4String SDname, G4SBSECaloutput ed ){
  ecaldata[SDname] = ed;
}

void G4SBSIO::SetSDtrackData( G4String SDname, G4SBSSDTrackOutput td ){
  sdtrackdata[SDname] = td;
}

void G4SBSIO::SetBDData(G4String SDname,G4SBSBDoutput data){
   BDdata[SDname] = data;
}

void G4SBSIO::SetICData(G4String SDname,G4SBSICoutput data){
   ICdata[SDname] = data;
}

void G4SBSIO::SetGEnTargetData_Glass(G4String SDname,G4SBSTargetoutput data){
   genTgtGCdata[SDname] = data;
}

void G4SBSIO::SetGEnTargetData_Cu(G4String SDname,G4SBSTargetoutput data){
   genTgtCUdata[SDname] = data;
}

void G4SBSIO::SetGEnTargetData_Al(G4String SDname,G4SBSTargetoutput data){
   genTgtALdata[SDname] = data;
}

void G4SBSIO::SetGEnTargetData_3He(G4String SDname,G4SBSTargetoutput data){
   genTgt3HEdata[SDname] = data;
}

void G4SBSIO::InitializeTree(){
  if( fFile ){
    fFile->Close();
    delete fFile;
  }

  fFile = new TFile(fFilename, "RECREATE"); 

  if( fTree ){ delete fTree; }

  Esum_histograms = new TClonesArray("TH1F",10);
  PulseShape_histograms = new TClonesArray("TH1F",10);

  fNhistograms = 0;
    
  fTree = new TTree("T", "Geant4 SBS Simulation");

  // Let's stop changing the ev_t data structure, because it screws up reading of the tree in the future. If we want to store any other event-specific information,
  // then let's add dedicated tree branches to hold said information:
  fTree->Branch("ev", &evdata, "count/D:rate/D:solang/D:sigma/D:W2/D:xbj/D:Q2/D:th/D:ph/D:Aperp/D:Apar/D:Pt/D:Pl/D:vx/D:vy/D:vz/D:ep/D:np/D:epx/D:epy/D:epz/D:npx/D:npy/D:npz/D:nth/D:nph/D:pmperp/D:pmpar/D:pmparsm/D:z/D:phperp/D:phih/D:phiS/D:thetaS/D:MX2/D:Sx/D:Sy/D:Sz/D:s/D:t/D:u/D:costhetaCM/D:Egamma/D:nucl/I:fnucl/I:hadr/I:earmaccept/I:harmaccept/I");
  //fTree->Branch("tr", &trdata, "x/D:y/D:xp/D:yp/D:tx/D:ty/D:txp/D:typ/D:hcal/I:bb/I:gemtr/I:hcx/D:hcy/D:bcx/D:bcy/D:hct/D:hctex/D:hclx/D:hcly/D:hclz/D:hcdang/D");
  //fTree->Branch("gen", &gendata, "thbb/D:thsbs/D:dbb/D:dsbs/D:dhcal/D:voffhcal/D:drich/D:dsbstrkr/D:Ebeam/D");

  //Going forward, all new event-level variables and tree branches you want to add should go in their own branches. We start with the new SIDIS ones:
  //The beam and target polarization information will be potentially useful for all generators
  fTree->Branch( "TargPol", &fTargPol, "TargPol/D" );
  fTree->Branch( "TargThetaSpin", &fTargThetaSpin, "TargThetaSpin/D");
  fTree->Branch( "TargPhiSpin", &fTargPhiSpin, "TargPhiSpin/D");
  fTree->Branch( "BeamPol", &fBeamPol, "BeamPol/D");
  fTree->Branch( "BeamThetaSpin", &fBeamThetaSpin, "BeamThetaSpin/D" );
  fTree->Branch( "BeamPhiSpin", &fBeamPhiSpin, "BeamPhiSpin/D" );

  if( fKineType == G4SBS::kSIDIS ){ //Create branches for Collins and Sivers asymmetries (eventually others, like quark flavor contributions to cross sections/asymmetries/etc)
    fTree->Branch("AUT_Collins", &fAUT_Collins, "AUT_Collins/D");
    fTree->Branch("AUT_Sivers", &fAUT_Sivers, "AUT_Sivers/D" );
    fTree->Branch("AUT_Collins_min", &fAUT_Collins_min, "AUT_Collins_min/D");
    fTree->Branch("AUT_Collins_max", &fAUT_Collins_max, "AUT_Collins_max/D");
    fTree->Branch("AUT_Sivers_min", &fAUT_Sivers_min, "AUT_Sivers_min/D");
    fTree->Branch("AUT_Sivers_max", &fAUT_Sivers_max, "AUT_Sivers_max/D");
  }
  
  //Instead of having the same tree structure as before, we want to dynamically generate tree branches depending on what kinds of detectors are present: Since we already require the ROOT libraries, we might as well use TStrings:
    
  //For all tree branches representing data in sensitive detectors, we want to grab the information from fdetcon->SDlist
  //Later, we will add other kinds of sensitive detectors:

  bool keepanysdtracks = false;
  
  for( set<G4String>::iterator d = (fdetcon->SDlist).begin(); d != (fdetcon->SDlist).end(); d++ ){
    //for( G4int idet=0; idet<fdetcon->fSDman->G
    G4String SDname = *d;
    G4SBS::SDet_t SDtype = (fdetcon->SDtype)[SDname];

    G4cout << "Initializing tree branches for Sensitive Detector " << SDname.data() << G4endl;

    switch( SDtype ){
    case G4SBS::kGEM: //GEM: Add branches for the GEM AND "tracker" branches:
      //Create "GEM output" and "Tracker Output" data structures and associate them with this sensitive detector name:
      GEMdata[SDname] = G4SBSGEMoutput();
      trackdata[SDname] = G4SBSTrackerOutput();
	
      BranchGEM(SDname);
	
      break;
    case G4SBS::kCAL: //"CAL": Add appropriate branches:
      //Initialize "CAL output" data structure and associate with this sensitive detector:
      CALdata[SDname] = G4SBSCALoutput();
      BranchCAL(SDname);
      break;
    case G4SBS::kRICH: //"RICH"
      richdata[SDname] = G4SBSRICHoutput();
      // Note: if any optical photon production mechanisms are ON, create RICH detector branches.
      // Otherwise, we don't need them.
      if( fUsingCerenkov || fUsingScintillation ) BranchRICH(SDname);
      break;
    case G4SBS::kECAL: //"ECAL"
      // Note: if any optical photon production mechanisms are ON, create ECAL detector branches.
      // Otherwise, we don't need them.
      ecaldata[SDname] = G4SBSECaloutput();
      if( fUsingCerenkov || fUsingScintillation ) BranchECAL(SDname);
      break;
    case G4SBS::kBD: 
      // Beam Diffuser (BD) 
      BDdata[SDname] = G4SBSBDoutput(); 
      BranchBD(SDname); 
      break; 
    case G4SBS::kIC: 
      // Ion chamber (IC) 
      ICdata[SDname] = G4SBSICoutput(); 
      BranchIC(SDname); 
      break; 
    case G4SBS::kTarget_GEn_Glass: 
      // GEn target glass cell 
      genTgtGCdata[SDname] = G4SBSTargetoutput(); 
      BranchGEnTarget_Glass(SDname); 
      break; 
    case G4SBS::kTarget_GEn_Cu: 
      // GEn target Cu  
      genTgtCUdata[SDname] = G4SBSTargetoutput(); 
      BranchGEnTarget_Cu(SDname); 
      break; 
    case G4SBS::kTarget_GEn_Al: 
      // GEn target Al  
      genTgtALdata[SDname] = G4SBSTargetoutput(); 
      BranchGEnTarget_Al(SDname); 
      break;
    case G4SBS::kTarget_GEn_3He: 
      // GEn target 3He  
      genTgt3HEdata[SDname] = G4SBSTargetoutput(); 
      BranchGEnTarget_3He(SDname); 
      break;
    }

    map<G4String,G4bool>::iterator keepsdflag = fKeepSDtracks.find( SDname );
    
    if( fKeepAllSDtracks || (keepsdflag != fKeepSDtracks.end() && keepsdflag->second ) ){
      sdtrackdata[SDname] = G4SBSSDTrackOutput(SDname);

      allsdtrackdata = G4SBSSDTrackOutput("all");

      if( SDtype == G4SBS::kGEM || SDtype == G4SBS::kCAL ||
	  fUsingCerenkov || fUsingScintillation ){
	// NOTE: the above if condition is written this way because if the SDtype is ECAL or RICH, we only want to enable the SDtrack branches
	// if at least one optical photon-producing process is turned on!
	keepanysdtracks = true;
	//BranchSDTracks( SDname );
      }
    }
  }

  if( keepanysdtracks ){
    BranchSDTracks();
  }
  
  if( fUsePythia ){
    BranchPythia();
  }

  if( fUseSIMC ){
    BranchSIMC();
  }
  
  return;
}

void G4SBSIO::FillTree(){
  if( !fTree ){ 
    fprintf(stderr, "Error %s: %s line %d - Trying to fill non-existant tree\n", __PRETTY_FUNCTION__, __FILE__, __LINE__ );
    return; 
  }

  fTree->Fill();
}

void G4SBSIO::WriteTree(){
  assert( fFile );
  assert( fTree );
  if( !fFile->IsOpen() ){
    G4cerr << "ERROR: " << __FILE__ << " line " << __LINE__ << ": TFile not open" << G4endl;
    exit(1);
  }

  fFile->cd();
  fTree->Write("T", TObject::kOverwrite);

  Esum_histograms->Compress();
  Esum_histograms->Write();
  PulseShape_histograms->Compress();
  PulseShape_histograms->Write();
    
  G4SBSRun::GetRun()->GetData()->Write("run_data", TObject::kOverwrite);

  // Produce and write out field map graphics
  fGlobalField->DebugField( gendata.thbb, gendata.thsbs );

  //Before we delete all our field maps, let's write out some reusable "local ones" from the global:
  if( fdetcon->fUseGlobalField && fWritePortableFieldMaps ){ //We are using a global field definition with "TOSCA" maps, and the user wants to write portable ones based on a subset of the global field volume:
    bool writeSBS = false;
    bool writeBB = false;
    G4SBS::Exp_t exp = fdetcon->fExpType;

    //If this is an experiment with SBS, write a reusable "local" sbs field map:
    writeSBS = exp != G4SBS::kC16 && exp != G4SBS::kGEMHCtest; //SBS magnet is used in all expts. except C16 and GEM test
    writeBB  = exp != G4SBS::kGEp && exp != G4SBS::kC16 && exp != G4SBS::kTDIS && exp != G4SBS::kNDVCS
      && exp != G4SBS::kGEMHCtest && exp != G4SBS::kGEPpositron; //BigBite is used in all experiments except C16, gem test, TDIS, NDVCS, and GEP positron

    //Let's look at BigBite first: a 0.6 x 3 x 3 m grid with 2.5-cm spacing:
    //Should we make the grid volume and spacing user-configurable? I think not, in the name of idiot-proofing. We should provide a command to toggle the writing of these maps on and off
    //though.
    if( writeBB ){
      G4cout << "Writing portable BigBite field map... " << G4endl;
      // fGlobalField->WriteFieldMapSection( "database/BBfield_temp.table", G4SBS::kEarm, fdetcon->fEArmBuilder->fBBang,
      // 					  fdetcon->fEArmBuilder->fBBdist - 0.75*CLHEP::m, fdetcon->fEArmBuilder->fBBdist + 2.25*CLHEP::m,
      // 					  3.0*CLHEP::m, 0.6*CLHEP::m, 24, 120, 120 );

      fGlobalField->WriteFieldMapSection( "database/BBfield_temp.table", G4SBS::kEarm, fdetcon->fEArmBuilder->fBBang,
					  fdetcon->fEArmBuilder->fBBdist, -0.75*CLHEP::m, 2.25*CLHEP::m,
					  3.0*CLHEP::m, 0.6*CLHEP::m, 24, 120, 120 );
      
      G4cout << "done" << G4endl;
    }
    
    //Now SBS: here we use a 1 x 3 x 3.5 m grid with 2.5-cm spacing:
    //As in the case of BB, we make the grid NOT user-configurable. We may revisit this later.
    //We should still allow the user to turn off the writing of these maps, though:
    
    if( writeSBS ){
      G4cout << "Writing portable SBS field map..." << G4endl;
      
      // fGlobalField->WriteFieldMapSection( "database/SBSfield_temp.table", G4SBS::kHarm, fdetcon->fHArmBuilder->f48D48ang,
      // 					  fdetcon->fHArmBuilder->f48D48dist-1.0*CLHEP::m, fdetcon->fHArmBuilder->f48D48dist + 2.5*CLHEP::m,
      // 					  3.0*CLHEP::m, 1.0*CLHEP::m, 40, 120, 140 );

      fGlobalField->WriteFieldMapSection( "database/SBSfield_temp.table", G4SBS::kHarm, fdetcon->fHArmBuilder->f48D48ang,
					  fdetcon->fHArmBuilder->f48D48dist, -1.0*CLHEP::m, 2.5*CLHEP::m,
					  3.0*CLHEP::m, 1.0*CLHEP::m, 40, 120, 140 );
      
      G4cout << "done" << G4endl;
    }
  }
  
  for( vector<TH2F *>::iterator it = fGlobalField->fFieldPlots.begin(); it!= fGlobalField->fFieldPlots.end(); it++ ){
    (*it)->Write((*it)->GetName(), TObject::kOverwrite );
    delete (*it);
  }
  fGlobalField->fFieldPlots.clear();

  fTree->ResetBranchAddresses();
  delete fTree;
  fTree = NULL;

  fFile->Close();
  delete fFile;
  fFile = NULL;

  return;
}

void G4SBSIO::BranchGEM(G4String SDname="GEM"){
  TString branch_prefix = SDname.data();
  TString branch_name;
  
  branch_prefix.ReplaceAll("/",".");
 
  //Branches with raw GEM data:
  
  fTree->Branch( branch_name.Format( "%s.hit.nhits", branch_prefix.Data() ), &(GEMdata[SDname].nhits_GEM) );
  fTree->Branch( branch_name.Format( "%s.hit.plane", branch_prefix.Data() ), &(GEMdata[SDname].plane) );
  fTree->Branch( branch_name.Format( "%s.hit.strip", branch_prefix.Data() ), &(GEMdata[SDname].strip) );
  fTree->Branch( branch_name.Format( "%s.hit.x", branch_prefix.Data() ), &(GEMdata[SDname].x) );
  fTree->Branch( branch_name.Format( "%s.hit.y", branch_prefix.Data() ), &(GEMdata[SDname].y) );
  fTree->Branch( branch_name.Format( "%s.hit.z", branch_prefix.Data() ), &(GEMdata[SDname].z) );
  fTree->Branch( branch_name.Format( "%s.hit.polx", branch_prefix.Data() ), &(GEMdata[SDname].polx) );
  fTree->Branch( branch_name.Format( "%s.hit.poly", branch_prefix.Data() ), &(GEMdata[SDname].poly) );
  fTree->Branch( branch_name.Format( "%s.hit.polz", branch_prefix.Data() ), &(GEMdata[SDname].polz) );
  fTree->Branch( branch_name.Format( "%s.hit.t", branch_prefix.Data() ), &(GEMdata[SDname].t) );
  fTree->Branch( branch_name.Format( "%s.hit.trms", branch_prefix.Data() ), &(GEMdata[SDname].trms) );
  fTree->Branch( branch_name.Format( "%s.hit.tmin", branch_prefix.Data() ), &(GEMdata[SDname].tmin) );
  fTree->Branch( branch_name.Format( "%s.hit.tmax", branch_prefix.Data() ), &(GEMdata[SDname].tmax) );
  // fTree->Branch( branch_name.Format( "%s.hit.dx", branch_prefix.Data() ), &(GEMdata[SDname].dx) );
  // fTree->Branch( branch_name.Format( "%s.hit.dy", branch_prefix.Data() ), &(GEMdata[SDname].dy) );
  fTree->Branch( branch_name.Format( "%s.hit.tx", branch_prefix.Data() ), &(GEMdata[SDname].tx) );
  fTree->Branch( branch_name.Format( "%s.hit.ty", branch_prefix.Data() ), &(GEMdata[SDname].ty) );
  fTree->Branch( branch_name.Format( "%s.hit.xin", branch_prefix.Data() ), &(GEMdata[SDname].xin) );
  fTree->Branch( branch_name.Format( "%s.hit.yin", branch_prefix.Data() ), &(GEMdata[SDname].yin) );
  fTree->Branch( branch_name.Format( "%s.hit.zin", branch_prefix.Data() ), &(GEMdata[SDname].zin) );
  fTree->Branch( branch_name.Format( "%s.hit.xout", branch_prefix.Data() ), &(GEMdata[SDname].xout) );
  fTree->Branch( branch_name.Format( "%s.hit.yout", branch_prefix.Data() ), &(GEMdata[SDname].yout) );
  fTree->Branch( branch_name.Format( "%s.hit.zout", branch_prefix.Data() ), &(GEMdata[SDname].zout) );
  fTree->Branch( branch_name.Format( "%s.hit.txp", branch_prefix.Data() ), &(GEMdata[SDname].txp) );
  fTree->Branch( branch_name.Format( "%s.hit.typ", branch_prefix.Data() ), &(GEMdata[SDname].typ) );
  fTree->Branch( branch_name.Format( "%s.hit.xg", branch_prefix.Data() ), &(GEMdata[SDname].xg) );
  fTree->Branch( branch_name.Format( "%s.hit.yg", branch_prefix.Data() ), &(GEMdata[SDname].yg) );
  fTree->Branch( branch_name.Format( "%s.hit.zg", branch_prefix.Data() ), &(GEMdata[SDname].zg) );
  fTree->Branch( branch_name.Format( "%s.hit.trid", branch_prefix.Data() ), &(GEMdata[SDname].trid) );
  fTree->Branch( branch_name.Format( "%s.hit.mid", branch_prefix.Data() ), &(GEMdata[SDname].mid) );
  fTree->Branch( branch_name.Format( "%s.hit.pid", branch_prefix.Data() ), &(GEMdata[SDname].pid) );
  fTree->Branch( branch_name.Format( "%s.hit.vx", branch_prefix.Data() ), &(GEMdata[SDname].vx) );
  fTree->Branch( branch_name.Format( "%s.hit.vy", branch_prefix.Data() ), &(GEMdata[SDname].vy) );
  fTree->Branch( branch_name.Format( "%s.hit.vz", branch_prefix.Data() ), &(GEMdata[SDname].vz) );
  fTree->Branch( branch_name.Format( "%s.hit.p", branch_prefix.Data() ), &(GEMdata[SDname].p) );
  fTree->Branch( branch_name.Format( "%s.hit.edep", branch_prefix.Data() ), &(GEMdata[SDname].edep) );
  fTree->Branch( branch_name.Format( "%s.hit.beta", branch_prefix.Data() ), &(GEMdata[SDname].beta) );

  map<G4String,G4bool>::iterator keepsdflag = fKeepSDtracks.find( SDname );
    
  if( fKeepAllSDtracks || (keepsdflag != fKeepSDtracks.end() && keepsdflag->second ) ){
    //Add "SD track" indices:
    fTree->Branch( branch_name.Format( "%s.hit.otridx", branch_prefix.Data() ), &(GEMdata[SDname].otridx) );
    fTree->Branch( branch_name.Format( "%s.hit.ptridx", branch_prefix.Data() ), &(GEMdata[SDname].ptridx) );
    fTree->Branch( branch_name.Format( "%s.hit.sdtridx", branch_prefix.Data() ), &(GEMdata[SDname].sdtridx) );
  }
  
  //Branches with "Tracker output" data:
  fTree->Branch( branch_name.Format("%s.Track.ntracks",branch_prefix.Data() ), &(trackdata[SDname].ntracks) );
  fTree->Branch( branch_name.Format("%s.Track.TID",branch_prefix.Data() ), &(trackdata[SDname].TrackTID) );
  fTree->Branch( branch_name.Format("%s.Track.PID",branch_prefix.Data() ), &(trackdata[SDname].TrackPID) );
  fTree->Branch( branch_name.Format("%s.Track.MID",branch_prefix.Data() ), &(trackdata[SDname].TrackMID) );
  fTree->Branch( branch_name.Format("%s.Track.NumHits",branch_prefix.Data() ), &(trackdata[SDname].NumHits) );
  fTree->Branch( branch_name.Format("%s.Track.NumPlanes",branch_prefix.Data() ), &(trackdata[SDname].NumPlanes) );
  fTree->Branch( branch_name.Format("%s.Track.NDF",branch_prefix.Data() ), &(trackdata[SDname].NDF) );
  fTree->Branch( branch_name.Format("%s.Track.Chi2fit",branch_prefix.Data() ), &(trackdata[SDname].Chi2fit) );
  fTree->Branch( branch_name.Format("%s.Track.Chi2true",branch_prefix.Data() ), &(trackdata[SDname].Chi2true) );
  fTree->Branch( branch_name.Format("%s.Track.X",branch_prefix.Data() ), &(trackdata[SDname].TrackX) );
  fTree->Branch( branch_name.Format("%s.Track.Y",branch_prefix.Data() ), &(trackdata[SDname].TrackY) );
  fTree->Branch( branch_name.Format("%s.Track.Xp",branch_prefix.Data() ), &(trackdata[SDname].TrackXp) );
  fTree->Branch( branch_name.Format("%s.Track.Yp",branch_prefix.Data() ), &(trackdata[SDname].TrackYp) );
  fTree->Branch( branch_name.Format("%s.Track.T",branch_prefix.Data() ), &(trackdata[SDname].TrackT) );
  fTree->Branch( branch_name.Format("%s.Track.P",branch_prefix.Data() ), &(trackdata[SDname].TrackP) );
  fTree->Branch( branch_name.Format("%s.Track.Sx",branch_prefix.Data() ), &(trackdata[SDname].TrackSx) );
  fTree->Branch( branch_name.Format("%s.Track.Sy",branch_prefix.Data() ), &(trackdata[SDname].TrackSy) );
  fTree->Branch( branch_name.Format("%s.Track.Sz",branch_prefix.Data() ), &(trackdata[SDname].TrackSz) );
  fTree->Branch( branch_name.Format("%s.Track.Xfit",branch_prefix.Data() ), &(trackdata[SDname].TrackXfit) );
  fTree->Branch( branch_name.Format("%s.Track.Yfit",branch_prefix.Data() ), &(trackdata[SDname].TrackYfit) );
  fTree->Branch( branch_name.Format("%s.Track.Xpfit",branch_prefix.Data() ), &(trackdata[SDname].TrackXpfit) );
  fTree->Branch( branch_name.Format("%s.Track.Ypfit",branch_prefix.Data() ), &(trackdata[SDname].TrackYpfit) );

  //map<G4String,G4bool>::iterator keepsdflag = fKeepSDtracks.find( SDname );
    
  if( fKeepAllSDtracks || (keepsdflag != fKeepSDtracks.end() && keepsdflag->second ) ){
    //Add "SD track" indices:
    fTree->Branch( branch_name.Format( "%s.Track.otridx", branch_prefix.Data() ), &(trackdata[SDname].otridx) );
    fTree->Branch( branch_name.Format( "%s.Track.ptridx", branch_prefix.Data() ), &(trackdata[SDname].ptridx) );
    fTree->Branch( branch_name.Format( "%s.Track.sdtridx", branch_prefix.Data() ), &(trackdata[SDname].sdtridx) );
  }
  
  map<G4String,G4bool>::iterator it = KeepHistoryflags.find( SDname );

  if( it != KeepHistoryflags.end() && it->second ){
    //Branches with "Particle History" data:
    fTree->Branch( branch_name.Format("%s.part.npart", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.npart) );
    fTree->Branch( branch_name.Format("%s.part.PID", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.PID) );
    fTree->Branch( branch_name.Format("%s.part.MID", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.MID) );
    fTree->Branch( branch_name.Format("%s.part.TID", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.TID) );
    fTree->Branch( branch_name.Format("%s.part.nbounce", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.nbounce) );
    fTree->Branch( branch_name.Format("%s.part.hitindex", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.hitindex) );
    fTree->Branch( branch_name.Format("%s.part.vx", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.vx) );
    fTree->Branch( branch_name.Format("%s.part.vy", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.vy) );
    fTree->Branch( branch_name.Format("%s.part.vz", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.vz) );
    fTree->Branch( branch_name.Format("%s.part.px", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.px) );
    fTree->Branch( branch_name.Format("%s.part.py", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.py) );
    fTree->Branch( branch_name.Format("%s.part.pz", branch_prefix.Data() ), &(GEMdata[SDname].ParticleHistory.pz) );
  }

  return;
}

void G4SBSIO::BranchCAL( G4String SDname="CAL" ){
  TString branch_prefix = SDname.data();
  TString branch_name;
  
  branch_prefix.ReplaceAll("/",".");

  TString histname;

  G4SBSCalSD *SDtemp = (G4SBSCalSD*) fdetcon->fSDman->FindSensitiveDetector( SDname );

  G4double threshtemp = SDtemp->GetEnergyThreshold();
  G4double gatewidthtemp = SDtemp->GetTimeWindow();
  G4int ntimebinstemp = SDtemp->GetNTimeBins();

  TString hname_prefix = branch_prefix;
  hname_prefix.ReplaceAll(".","_");
  new( (*Esum_histograms)[fNhistograms] ) TH1F( histname.Format("%s_esum",hname_prefix.Data()), "",
						100, 0.0, 100.0*threshtemp );
  new( (*PulseShape_histograms)[fNhistograms] ) TH1F( histname.Format("%s_pulseshape",hname_prefix.Data()), "",
						      ntimebinstemp, 0.0,
						      gatewidthtemp );

  histogram_index[SDname] = fNhistograms;
  
  fNhistograms++;
			  
  //Define "hit" branches:
  fTree->Branch( branch_name.Format( "%s.det.esum", branch_prefix.Data() ), &(CALdata[SDname].Esum) );
  fTree->Branch( branch_name.Format( "%s.hit.nhits", branch_prefix.Data() ), &(CALdata[SDname].nhits_CAL) );
  fTree->Branch( branch_name.Format( "%s.hit.row", branch_prefix.Data() ), &(CALdata[SDname].row) );
  fTree->Branch( branch_name.Format( "%s.hit.col", branch_prefix.Data() ), &(CALdata[SDname].col) );
  fTree->Branch( branch_name.Format( "%s.hit.cell", branch_prefix.Data() ), &(CALdata[SDname].cell) );
  fTree->Branch( branch_name.Format( "%s.hit.plane", branch_prefix.Data() ), &(CALdata[SDname].plane) );
  fTree->Branch( branch_name.Format( "%s.hit.wire", branch_prefix.Data() ), &(CALdata[SDname].wire) );
  fTree->Branch( branch_name.Format( "%s.hit.xcell", branch_prefix.Data() ), &(CALdata[SDname].xcell) );
  fTree->Branch( branch_name.Format( "%s.hit.ycell", branch_prefix.Data() ), &(CALdata[SDname].ycell) );
  fTree->Branch( branch_name.Format( "%s.hit.zcell", branch_prefix.Data() ), &(CALdata[SDname].zcell) );
  fTree->Branch( branch_name.Format( "%s.hit.xcellg", branch_prefix.Data() ), &(CALdata[SDname].xcellg) );
  fTree->Branch( branch_name.Format( "%s.hit.ycellg", branch_prefix.Data() ), &(CALdata[SDname].ycellg) );
  fTree->Branch( branch_name.Format( "%s.hit.zcellg", branch_prefix.Data() ), &(CALdata[SDname].zcellg) );
  fTree->Branch( branch_name.Format( "%s.hit.xhit", branch_prefix.Data() ), &(CALdata[SDname].xhit) );
  fTree->Branch( branch_name.Format( "%s.hit.yhit", branch_prefix.Data() ), &(CALdata[SDname].yhit) );
  fTree->Branch( branch_name.Format( "%s.hit.zhit", branch_prefix.Data() ), &(CALdata[SDname].zhit) );
  fTree->Branch( branch_name.Format( "%s.hit.xhitg", branch_prefix.Data() ), &(CALdata[SDname].xhitg) );
  fTree->Branch( branch_name.Format( "%s.hit.yhitg", branch_prefix.Data() ), &(CALdata[SDname].yhitg) );
  fTree->Branch( branch_name.Format( "%s.hit.zhitg", branch_prefix.Data() ), &(CALdata[SDname].zhitg) );
  fTree->Branch( branch_name.Format( "%s.hit.sumedep", branch_prefix.Data() ), &(CALdata[SDname].sumedep) );
  fTree->Branch( branch_name.Format( "%s.hit.tavg", branch_prefix.Data() ), &(CALdata[SDname].tavg) );
  fTree->Branch( branch_name.Format( "%s.hit.trms", branch_prefix.Data() ), &(CALdata[SDname].trms) );
  fTree->Branch( branch_name.Format( "%s.hit.tmin", branch_prefix.Data() ), &(CALdata[SDname].tmin) );
  fTree->Branch( branch_name.Format( "%s.hit.tmax", branch_prefix.Data() ), &(CALdata[SDname].tmax) );

  // Fill in ROOT tree branch to hold Pulse Shape info 
  map<G4String,G4bool>::iterator keeppsflag = fKeepPulseShape.find( SDname );    
  if( fKeepAllPulseShape || (keeppsflag != fKeepPulseShape.end() && keeppsflag->second ) ){
    fTree->Branch( branch_name.Format( "%s.gatewidth", branch_prefix.Data() ), &(CALdata[SDname].gatewidth) );
    fTree->Branch( branch_name.Format( "%s.hit.edep_vs_time", branch_prefix.Data() ), &(CALdata[SDname].edep_vs_time) );
  }

  map<G4String,G4bool>::iterator keepsdflag = fKeepSDtracks.find( SDname );
    
  if( fKeepAllSDtracks || (keepsdflag != fKeepSDtracks.end() && keepsdflag->second ) ){
    //Add "SD track" indices:
    fTree->Branch( branch_name.Format( "%s.hit.otridx", branch_prefix.Data() ), &(CALdata[SDname].otridx) );
    fTree->Branch( branch_name.Format( "%s.hit.ptridx", branch_prefix.Data() ), &(CALdata[SDname].ptridx) );
    fTree->Branch( branch_name.Format( "%s.hit.sdtridx", branch_prefix.Data() ), &(CALdata[SDname].sdtridx) );
  }
  
  map<G4String,G4bool>::iterator it = KeepPartCALflags.find( SDname );

  if( it != KeepPartCALflags.end() && it->second ){
    //Define "particle" branches:
    fTree->Branch( branch_name.Format( "%s.npart_CAL", branch_prefix.Data() ), &(CALdata[SDname].npart_CAL) );
    fTree->Branch( branch_name.Format( "%s.ihit", branch_prefix.Data() ), &(CALdata[SDname].ihit) );
    fTree->Branch( branch_name.Format( "%s.x", branch_prefix.Data() ), &(CALdata[SDname].x) );
    fTree->Branch( branch_name.Format( "%s.y", branch_prefix.Data() ), &(CALdata[SDname].y) );
    fTree->Branch( branch_name.Format( "%s.z", branch_prefix.Data() ), &(CALdata[SDname].z) );
    fTree->Branch( branch_name.Format( "%s.t", branch_prefix.Data() ), &(CALdata[SDname].t) );
    fTree->Branch( branch_name.Format( "%s.E", branch_prefix.Data() ), &(CALdata[SDname].E) );
    fTree->Branch( branch_name.Format( "%s.dt", branch_prefix.Data() ), &(CALdata[SDname].dt) );
    fTree->Branch( branch_name.Format( "%s.L", branch_prefix.Data() ), &(CALdata[SDname].L) );
    fTree->Branch( branch_name.Format( "%s.vx", branch_prefix.Data() ), &(CALdata[SDname].vx) );
    fTree->Branch( branch_name.Format( "%s.vy", branch_prefix.Data() ), &(CALdata[SDname].vy) );
    fTree->Branch( branch_name.Format( "%s.vz", branch_prefix.Data() ), &(CALdata[SDname].vz) );
    fTree->Branch( branch_name.Format( "%s.trid", branch_prefix.Data() ),  &(CALdata[SDname].trid) );
    fTree->Branch( branch_name.Format( "%s.mid", branch_prefix.Data() ), &(CALdata[SDname].mid) );
    fTree->Branch( branch_name.Format( "%s.pid", branch_prefix.Data() ), &(CALdata[SDname].pid) );
    fTree->Branch( branch_name.Format( "%s.p", branch_prefix.Data() ), &(CALdata[SDname].p) );
    fTree->Branch( branch_name.Format( "%s.px", branch_prefix.Data() ), &(CALdata[SDname].px) );
    fTree->Branch( branch_name.Format( "%s.py", branch_prefix.Data() ), &(CALdata[SDname].py) );
    fTree->Branch( branch_name.Format( "%s.pz", branch_prefix.Data() ), &(CALdata[SDname].pz) );
    fTree->Branch( branch_name.Format( "%s.edep", branch_prefix.Data() ), &(CALdata[SDname].edep) );
  }

  it = KeepHistoryflags.find( SDname );

  if( it != KeepHistoryflags.end() && it->second ){
    //Branches with "Particle History" data:
    fTree->Branch( branch_name.Format("%s.part.npart", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.npart) );
    fTree->Branch( branch_name.Format("%s.part.PID", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.PID) );
    fTree->Branch( branch_name.Format("%s.part.MID", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.MID) );
    fTree->Branch( branch_name.Format("%s.part.TID", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.TID) );
    fTree->Branch( branch_name.Format("%s.part.nbounce", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.nbounce) );
    fTree->Branch( branch_name.Format("%s.part.hitindex", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.hitindex) );
    fTree->Branch( branch_name.Format("%s.part.vx", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.vx) );
    fTree->Branch( branch_name.Format("%s.part.vy", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.vy) );
    fTree->Branch( branch_name.Format("%s.part.vz", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.vz) );
    fTree->Branch( branch_name.Format("%s.part.px", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.px) );
    fTree->Branch( branch_name.Format("%s.part.py", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.py) );
    fTree->Branch( branch_name.Format("%s.part.pz", branch_prefix.Data() ), &(CALdata[SDname].ParticleHistory.pz) );
  }

  return;
}

void G4SBSIO::BranchRICH(G4String SDname="RICH"){
  TString branch_prefix = SDname.data();
  TString branch_name;
  branch_prefix.ReplaceAll("/",".");
  
  //Branches for "hits": 
  
  fTree->Branch( branch_name.Format("%s.hit.nhits", branch_prefix.Data() ), &(richdata[SDname].nhits_RICH) );
  fTree->Branch( branch_name.Format("%s.hit.PMT", branch_prefix.Data() ), &(richdata[SDname].PMTnumber) );
  fTree->Branch( branch_name.Format("%s.hit.row", branch_prefix.Data() ), &(richdata[SDname].row) );
  fTree->Branch( branch_name.Format("%s.hit.col", branch_prefix.Data() ), &(richdata[SDname].col) );
  fTree->Branch( branch_name.Format("%s.hit.xpmt", branch_prefix.Data() ), &(richdata[SDname].xpmt) );
  fTree->Branch( branch_name.Format("%s.hit.ypmt", branch_prefix.Data() ), &(richdata[SDname].ypmt) );
  fTree->Branch( branch_name.Format("%s.hit.zpmt", branch_prefix.Data() ), &(richdata[SDname].zpmt) );
  fTree->Branch( branch_name.Format( "%s.hit.xgpmt", branch_prefix.Data() ), &(richdata[SDname].xgpmt) );
  fTree->Branch( branch_name.Format( "%s.hit.ygpmt", branch_prefix.Data() ), &(richdata[SDname].ygpmt) );
  fTree->Branch( branch_name.Format( "%s.hit.zgpmt", branch_prefix.Data() ), &(richdata[SDname].zgpmt) );
  fTree->Branch( branch_name.Format("%s.hit.NumPhotoelectrons", branch_prefix.Data() ), &(richdata[SDname].NumPhotoelectrons) );
  fTree->Branch( branch_name.Format("%s.hit.Time_avg", branch_prefix.Data() ), &(richdata[SDname].Time_avg) );
  fTree->Branch( branch_name.Format("%s.hit.Time_rms", branch_prefix.Data() ), &(richdata[SDname].Time_rms) );
  fTree->Branch( branch_name.Format("%s.hit.Time_min", branch_prefix.Data() ), &(richdata[SDname].Time_min) );
  fTree->Branch( branch_name.Format("%s.hit.Time_max", branch_prefix.Data() ), &(richdata[SDname].Time_max) );
  fTree->Branch( branch_name.Format("%s.hit.mTrackNo", branch_prefix.Data() ), &(richdata[SDname].mTrackNo) );
  fTree->Branch( branch_name.Format("%s.hit.xhit", branch_prefix.Data() ), &(richdata[SDname].xhit) );
  fTree->Branch( branch_name.Format("%s.hit.yhit", branch_prefix.Data() ), &(richdata[SDname].yhit) );
  fTree->Branch( branch_name.Format("%s.hit.zhit", branch_prefix.Data() ), &(richdata[SDname].zhit) );
  fTree->Branch( branch_name.Format("%s.hit.pxhit", branch_prefix.Data() ), &(richdata[SDname].pxhit) );
  fTree->Branch( branch_name.Format("%s.hit.pyhit", branch_prefix.Data() ), &(richdata[SDname].pyhit) );
  fTree->Branch( branch_name.Format("%s.hit.pzhit", branch_prefix.Data() ), &(richdata[SDname].pzhit) );
  fTree->Branch( branch_name.Format("%s.hit.pvx", branch_prefix.Data() ), &(richdata[SDname].pvx) );
  fTree->Branch( branch_name.Format("%s.hit.pvy", branch_prefix.Data() ), &(richdata[SDname].pvy) );
  fTree->Branch( branch_name.Format("%s.hit.pvz", branch_prefix.Data() ), &(richdata[SDname].pvz) );
  fTree->Branch( branch_name.Format("%s.hit.ppx", branch_prefix.Data() ), &(richdata[SDname].ppx) );
  fTree->Branch( branch_name.Format("%s.hit.ppy", branch_prefix.Data() ), &(richdata[SDname].ppy) );
  fTree->Branch( branch_name.Format("%s.hit.ppz", branch_prefix.Data() ), &(richdata[SDname].ppz) );
  fTree->Branch( branch_name.Format("%s.hit.volume_flag", branch_prefix.Data() ), &(richdata[SDname].volume_flag) );

  map<G4String,G4bool>::iterator keepsdflag = fKeepSDtracks.find( SDname );
    
  if( fKeepAllSDtracks || (keepsdflag != fKeepSDtracks.end() && keepsdflag->second ) ){
    //Add "SD track" indices:
    fTree->Branch( branch_name.Format( "%s.hit.otridx", branch_prefix.Data() ), &(richdata[SDname].otridx) );
    fTree->Branch( branch_name.Format( "%s.hit.ptridx", branch_prefix.Data() ), &(richdata[SDname].ptridx) );
    fTree->Branch( branch_name.Format( "%s.hit.sdtridx", branch_prefix.Data() ), &(richdata[SDname].sdtridx) );
  }
  //Branches for "tracks": This might be reorganized later:
  // branch_name.Format( "%s.ntracks_RICH", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].ntracks_RICH), "ntracks_RICH/I" );
  // branch_name.Format( "%s.mPID", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mPID) );
  // branch_name.Format( "%s.mTID", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mTID) );
  // branch_name.Format( "%s.mMID", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mMID) );
  // branch_name.Format( "%s.mvx", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mvx) );
  // branch_name.Format( "%s.mvy", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mvy) );
  // branch_name.Format( "%s.mvz", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mvz) );
  // branch_name.Format( "%s.mpx", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mpx) );
  // branch_name.Format( "%s.mpy", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mpy) );
  // branch_name.Format( "%s.mpz", branch_prefix.Data() );
  // fTree->Branch( branch_name.Data(), &(richdata[SDname].mpz) );

  map<G4String,G4bool>::iterator it = KeepHistoryflags.find( SDname );

  if( it != KeepHistoryflags.end() && it->second ){
    //Branches with "Particle History" data:
    fTree->Branch( branch_name.Format("%s.part.npart", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.npart) );
    fTree->Branch( branch_name.Format("%s.part.PID", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.PID) );
    fTree->Branch( branch_name.Format("%s.part.MID", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.MID) );
    fTree->Branch( branch_name.Format("%s.part.TID", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.TID) );
    fTree->Branch( branch_name.Format("%s.part.nbounce", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.nbounce) );
    fTree->Branch( branch_name.Format("%s.part.hitindex", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.hitindex) );
    fTree->Branch( branch_name.Format("%s.part.vx", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.vx) );
    fTree->Branch( branch_name.Format("%s.part.vy", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.vy) );
    fTree->Branch( branch_name.Format("%s.part.vz", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.vz) );
    fTree->Branch( branch_name.Format("%s.part.px", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.px) );
    fTree->Branch( branch_name.Format("%s.part.py", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.py) );
    fTree->Branch( branch_name.Format("%s.part.pz", branch_prefix.Data() ), &(richdata[SDname].ParticleHistory.pz) );
    fTree->Branch( branch_name.Format("%s.part.Nphe_part", branch_prefix.Data() ), &(richdata[SDname].Nphe_part) );
  }
  return;
}

void G4SBSIO::BranchECAL(G4String SDname="ECAL"){
  TString branch_prefix = SDname.data();
  TString branch_name;
  branch_prefix.ReplaceAll("/",".");

  // *****
  G4SBSECalSD *SDtemp = (G4SBSECalSD*) fdetcon->fSDman->FindSensitiveDetector( SDname );

  G4double threshtemp = SDtemp->GetPEThreshold();
  G4double gatewidthtemp = SDtemp->GetTimeWindow();
  G4int ntimebinstemp = SDtemp->GetNTimeBins();
  // *****
  
  fTree->Branch( branch_name.Format("%s.hit.nhits", branch_prefix.Data() ), &(ecaldata[SDname].nhits_ECal) );
  fTree->Branch( branch_name.Format("%s.hit.PMT", branch_prefix.Data() ), &(ecaldata[SDname].PMTnumber) );
  fTree->Branch( branch_name.Format("%s.hit.row", branch_prefix.Data() ), &(ecaldata[SDname].row) );
  fTree->Branch( branch_name.Format("%s.hit.col", branch_prefix.Data() ), &(ecaldata[SDname].col) );
  fTree->Branch( branch_name.Format("%s.hit.plane", branch_prefix.Data() ), &(ecaldata[SDname].plane) );
  fTree->Branch( branch_name.Format("%s.hit.xcell", branch_prefix.Data() ), &(ecaldata[SDname].xcell) );
  fTree->Branch( branch_name.Format("%s.hit.ycell", branch_prefix.Data() ), &(ecaldata[SDname].ycell) );
  fTree->Branch( branch_name.Format("%s.hit.zcell", branch_prefix.Data() ), &(ecaldata[SDname].zcell) );
  fTree->Branch( branch_name.Format("%s.hit.xgcell", branch_prefix.Data() ), &(ecaldata[SDname].xgcell) );
  fTree->Branch( branch_name.Format("%s.hit.ygcell", branch_prefix.Data() ), &(ecaldata[SDname].ygcell) );
  fTree->Branch( branch_name.Format("%s.hit.zgcell", branch_prefix.Data() ), &(ecaldata[SDname].zgcell) );
  fTree->Branch( branch_name.Format("%s.hit.NumPhotoelectrons", branch_prefix.Data() ), &(ecaldata[SDname].NumPhotoelectrons) );
  fTree->Branch( branch_name.Format("%s.hit.Time_avg", branch_prefix.Data() ), &(ecaldata[SDname].Time_avg) );
  fTree->Branch( branch_name.Format("%s.hit.Time_rms", branch_prefix.Data() ), &(ecaldata[SDname].Time_rms) );
  fTree->Branch( branch_name.Format("%s.hit.Time_min", branch_prefix.Data() ), &(ecaldata[SDname].Time_min) );
  fTree->Branch( branch_name.Format("%s.hit.Time_max", branch_prefix.Data() ), &(ecaldata[SDname].Time_max) );

  // *****
  // Fill in ROOT tree branch to hold Pulse Shape info 
  map<G4String,G4bool>::iterator keeppsflag = fKeepPulseShape.find( SDname );    
  if( fKeepAllPulseShape || (keeppsflag != fKeepPulseShape.end() && keeppsflag->second ) ){
    fTree->Branch( branch_name.Format( "%s.gatewidth", branch_prefix.Data() ), &(ecaldata[SDname].gatewidth) );
    fTree->Branch( branch_name.Format( "%s.hit.NPE_vs_time", branch_prefix.Data() ), &(ecaldata[SDname].NPE_vs_time) );
  }
  // *****
  
  map<G4String,G4bool>::iterator keepsdflag = fKeepSDtracks.find( SDname );
    
  if( fKeepAllSDtracks || (keepsdflag != fKeepSDtracks.end() && keepsdflag->second ) ){
    //Add "SD track" indices:
    fTree->Branch( branch_name.Format( "%s.hit.otridx", branch_prefix.Data() ), &(ecaldata[SDname].otridx) );
    fTree->Branch( branch_name.Format( "%s.hit.ptridx", branch_prefix.Data() ), &(ecaldata[SDname].ptridx) );
    fTree->Branch( branch_name.Format( "%s.hit.sdtridx", branch_prefix.Data() ), &(ecaldata[SDname].sdtridx) );
  }
  
  map<G4String,G4bool>::iterator it = KeepPartCALflags.find( SDname );

  if( it != KeepPartCALflags.end() && it->second ){
    //Define "particle" branches:
    fTree->Branch( branch_name.Format( "%s.npart_ECAL", branch_prefix.Data() ), &(ecaldata[SDname].npart_ECAL) );
    fTree->Branch( branch_name.Format( "%s.part_PMT", branch_prefix.Data() ), &(ecaldata[SDname].part_PMT) );
    //fTree->Branch( branch_name.Format( "%s.ihit", branch_prefix.Data() ), &(ecaldata[SDname].ihit) );
    fTree->Branch( branch_name.Format( "%s.E", branch_prefix.Data() ), &(ecaldata[SDname].E) );
    fTree->Branch( branch_name.Format( "%s.t", branch_prefix.Data() ), &(ecaldata[SDname].t) );
    fTree->Branch( branch_name.Format( "%s.trid", branch_prefix.Data() ), &(ecaldata[SDname].trid) );
    fTree->Branch( branch_name.Format( "%s.detected", branch_prefix.Data() ), &(ecaldata[SDname].detected) );
  }
}

void G4SBSIO::BranchPythia(){
  //Per-event variables:
  fTree->Branch("primaries.Sigma",&(Primaries.Sigma),"primaries.Sigma/D");
  fTree->Branch("primaries.Ebeam",&(Primaries.Ebeam),"primaries.Ebeam/D");
  fTree->Branch("primaries.Eprime",&(Primaries.Eprime),"primaries.Eprime/D");
  fTree->Branch("primaries.Q2",&(Primaries.Q2),"primaries.Q2/D");
  fTree->Branch("primaries.xbj",&(Primaries.xbj),"primaries.xbj/D");
  fTree->Branch("primaries.y",&(Primaries.y),"primaries.y/D");
  fTree->Branch("primaries.W2",&(Primaries.W2),"primaries.W2/D");
  fTree->Branch("primaries.theta_e",&(Primaries.theta_e),"primaries.theta_e/D");
  fTree->Branch("primaries.phi_e",&(Primaries.phi_e),"primaries.phi_e/D");
  fTree->Branch("primaries.px_e",&(Primaries.px_e),"primaries.px_e/D");
  fTree->Branch("primaries.py_e",&(Primaries.py_e),"primaries.py_e/D");
  fTree->Branch("primaries.pz_e",&(Primaries.pz_e),"primaries.pz_e/D");
  fTree->Branch("primaries.vx_e",&(Primaries.vx_e),"primaries.vx_e/D");
  fTree->Branch("primaries.vy_e",&(Primaries.vy_e),"primaries.vy_e/D");
  fTree->Branch("primaries.vz_e",&(Primaries.vz_e),"primaries.vz_e/D");
  fTree->Branch("primaries.Egamma",&(Primaries.Egamma),"primaries.Egamma/D");
  fTree->Branch("primaries.theta_gamma",&(Primaries.theta_gamma),"primaries.theta_gamma/D");
  fTree->Branch("primaries.phi_gamma",&(Primaries.phi_gamma),"primaries.phi_gamma/D");
  fTree->Branch("primaries.px_gamma",&(Primaries.px_gamma),"primaries.px_gamma/D");
  fTree->Branch("primaries.py_gamma",&(Primaries.py_gamma),"primaries.py_gamma/D");
  fTree->Branch("primaries.pz_gamma",&(Primaries.pz_gamma),"primaries.pz_gamma/D");
  fTree->Branch("primaries.vx_gamma",&(Primaries.vx_gamma),"primaries.vx_gamma/D");
  fTree->Branch("primaries.vy_gamma",&(Primaries.vy_gamma),"primaries.vy_gamma/D");
  fTree->Branch("primaries.vz_gamma",&(Primaries.vz_gamma),"primaries.vz_gamma/D");
  //Primary particle arrays:
  fTree->Branch("Primaries.Nprimaries",&(Primaries.Nprimaries),"Nprimaries/I");
  fTree->Branch("Primaries.PID",&(Primaries.PID));
  fTree->Branch("Primaries.genflag",&(Primaries.genflag));
  fTree->Branch("Primaries.Px",&(Primaries.Px));
  fTree->Branch("Primaries.Py",&(Primaries.Py));
  fTree->Branch("Primaries.Pz",&(Primaries.Pz));
  fTree->Branch("Primaries.vx",&(Primaries.vx));
  fTree->Branch("Primaries.vy",&(Primaries.vy));
  fTree->Branch("Primaries.vz",&(Primaries.vz));
  fTree->Branch("Primaries.M",&(Primaries.M));
  fTree->Branch("Primaries.E",&(Primaries.E));
  fTree->Branch("Primaries.P",&(Primaries.P));
  fTree->Branch("Primaries.t",&(Primaries.t));
  fTree->Branch("Primaries.theta",&(Primaries.theta));
  fTree->Branch("Primaries.phi",&(Primaries.phi));
}

void G4SBSIO::BranchSIMC(){
  fTree->Branch("simc.sigma",&(SIMCprimaries.sigma),"simc.sigma/D");
  fTree->Branch("simc.Weight",&(SIMCprimaries.Weight),"simc.Weight/D");
  fTree->Branch("simc.Q2",&(SIMCprimaries.Q2),"simc.Q2/D");
  fTree->Branch("simc.xbj",&(SIMCprimaries.xbj),"simc.xbj/D");
  fTree->Branch("simc.nu",&(SIMCprimaries.nu),"simc.nu/D");
  fTree->Branch("simc.W",&(SIMCprimaries.W),"simc.W/D");
  fTree->Branch("simc.epsilon",&(SIMCprimaries.epsilon),"simc.epsilon/D");
  
  fTree->Branch("simc.Ebeam",&(SIMCprimaries.Ebeam),"simc.Ebeam/D");
  
  fTree->Branch("simc.fnucl",&(SIMCprimaries.fnucl),"simc.fnucl/I");
  fTree->Branch("simc.p_e",&(SIMCprimaries.p_e),"simc.p_e/D");
  fTree->Branch("simc.theta_e",&(SIMCprimaries.theta_e),"simc.theta_e/D");
  fTree->Branch("simc.phi_e",&(SIMCprimaries.phi_e),"simc.phi_e/D");
  fTree->Branch("simc.px_e",&(SIMCprimaries.px_e),"simc.px_e/D");
  fTree->Branch("simc.py_e",&(SIMCprimaries.py_e),"simc.py_e/D");
  fTree->Branch("simc.pz_e",&(SIMCprimaries.pz_e),"simc.pz_e/D");
  fTree->Branch("simc.p_n",&(SIMCprimaries.p_n),"simc.p_n/D");
  fTree->Branch("simc.theta_n",&(SIMCprimaries.theta_n),"simc.theta_n/D");
  fTree->Branch("simc.phi_n",&(SIMCprimaries.phi_n),"simc.phi_n/D");
  fTree->Branch("simc.px_n",&(SIMCprimaries.px_n),"simc.px_n/D");
  fTree->Branch("simc.py_n",&(SIMCprimaries.py_n),"simc.py_n/D");
  fTree->Branch("simc.pz_n",&(SIMCprimaries.pz_n),"simc.pz_n/D");
  
  fTree->Branch("simc.vx",&(SIMCprimaries.vx),"simc.vx/D");
  fTree->Branch("simc.vy",&(SIMCprimaries.vy),"simc.vy/D");
  fTree->Branch("simc.vz",&(SIMCprimaries.vz),"simc.vz/D");

  fTree->Branch("simc.veE",&(SIMCprimaries.veE),"simc.veE/D");
  fTree->Branch("simc.vetheta",&(SIMCprimaries.vetheta),"simc.vetheta/D");
}

void G4SBSIO::UpdateGenDataFromDetCon(){ //Go with whatever is in fdetcon as of run start for constant parameters of the run describing detector layout:
  gendata.thbb = fdetcon->fEArmBuilder->fBBang;
  gendata.dbb  = fdetcon->fEArmBuilder->fBBdist/CLHEP::m;
  gendata.thsbs = fdetcon->fHArmBuilder->f48D48ang;
  gendata.dsbs  = fdetcon->fHArmBuilder->f48D48dist/CLHEP::m;
  gendata.dhcal = fdetcon->fHArmBuilder->fHCALdist/CLHEP::m;
  gendata.voffhcal = fdetcon->fHArmBuilder->fHCALvertical_offset/CLHEP::m;
  gendata.hoffhcal = fdetcon->fHArmBuilder->fHCALhorizontal_offset/CLHEP::m;
  gendata.angoffhcal = fdetcon->fHArmBuilder->fHCALangular_offset;
  gendata.dlac = fdetcon->fHArmBuilder->fLACdist/CLHEP::m;
  gendata.vofflac = fdetcon->fHArmBuilder->fLACvertical_offset/CLHEP::m;
  gendata.hofflac = fdetcon->fHArmBuilder->fLAChorizontal_offset/CLHEP::m;
  gendata.drich = fdetcon->fHArmBuilder->fRICHdist/CLHEP::m;
  gendata.dsbstrkr = fdetcon->fHArmBuilder->fSBS_tracker_dist/CLHEP::m;
  gendata.sbstrkrpitch = fdetcon->fHArmBuilder->fSBS_tracker_pitch; //radians
}

void G4SBSIO::BranchSDTracks(){
  //TString branch_prefix = "AllSD";
  //TString branch_name;
  //branch_prefix.ReplaceAll("/",".");
  
  //  map<G4String,G4bool>::iterator k = KeepSDtracks.find( SDname );
  
  //if( sdtrackdata.find( SDname ) != sdtrackdata.end() && (fKeepSDtracks.find(SDname)->second || fKeepAllSDtracks) ){
  //"Original track" info:
  fTree->Branch(  "OTrack.ntracks", &(allsdtrackdata.notracks) );
  fTree->Branch(  "OTrack.TID", &(allsdtrackdata.otrid) );
  fTree->Branch(  "OTrack.MID", &(allsdtrackdata.omid) );
  fTree->Branch(  "OTrack.PID", &(allsdtrackdata.opid) );
  fTree->Branch(  "OTrack.MPID", &(allsdtrackdata.ompid) );
  fTree->Branch(  "OTrack.posx", &(allsdtrackdata.oposx) );
  fTree->Branch(  "OTrack.posy", &(allsdtrackdata.oposy) );
  fTree->Branch(  "OTrack.posz", &(allsdtrackdata.oposz) );
  fTree->Branch(  "OTrack.momx", &(allsdtrackdata.omomx) );
  fTree->Branch(  "OTrack.momy", &(allsdtrackdata.omomy) );
  fTree->Branch(  "OTrack.momz", &(allsdtrackdata.omomz) );
  fTree->Branch(  "OTrack.polx", &(allsdtrackdata.opolx) );
  fTree->Branch(  "OTrack.poly", &(allsdtrackdata.opoly) );
  fTree->Branch(  "OTrack.polz", &(allsdtrackdata.opolz) );
  fTree->Branch(  "OTrack.Etot", &(allsdtrackdata.oenergy) );
  fTree->Branch(  "OTrack.T", &(allsdtrackdata.otime) );

  //"Primary track" info:
  fTree->Branch(  "PTrack.ntracks", &(allsdtrackdata.nptracks) );
  fTree->Branch(  "PTrack.TID", &(allsdtrackdata.ptrid) );
  //fTree->Branch(  "PTrack.MID", &(allsdtrackdata.pmid) );
  fTree->Branch(  "PTrack.PID", &(allsdtrackdata.ppid) );
  fTree->Branch(  "PTrack.posx", &(allsdtrackdata.pposx) );
  fTree->Branch(  "PTrack.posy", &(allsdtrackdata.pposy) );
  fTree->Branch(  "PTrack.posz", &(allsdtrackdata.pposz) );
  fTree->Branch(  "PTrack.momx", &(allsdtrackdata.pmomx) );
  fTree->Branch(  "PTrack.momy", &(allsdtrackdata.pmomy) );
  fTree->Branch(  "PTrack.momz", &(allsdtrackdata.pmomz) );
  fTree->Branch(  "PTrack.polx", &(allsdtrackdata.ppolx) );
  fTree->Branch(  "PTrack.poly", &(allsdtrackdata.ppoly) );
  fTree->Branch(  "PTrack.polz", &(allsdtrackdata.ppolz) );
  fTree->Branch(  "PTrack.Etot", &(allsdtrackdata.penergy) );
  fTree->Branch(  "PTrack.T", &(allsdtrackdata.ptime) );
  //"SD boundary crossing track" info:
  fTree->Branch(  "SDTrack.ntracks", &(allsdtrackdata.nsdtracks) );
  fTree->Branch(  "SDTrack.TID", &(allsdtrackdata.sdtrid) );
  fTree->Branch(  "SDTrack.MID", &(allsdtrackdata.sdmid) );
  fTree->Branch(  "SDTrack.PID", &(allsdtrackdata.sdpid) );
  fTree->Branch(  "SDTrack.MPID", &(allsdtrackdata.sdmpid) );
  fTree->Branch(  "SDTrack.posx", &(allsdtrackdata.sdposx) );
  fTree->Branch(  "SDTrack.posy", &(allsdtrackdata.sdposy) );
  fTree->Branch(  "SDTrack.posz", &(allsdtrackdata.sdposz) );
  fTree->Branch(  "SDTrack.momx", &(allsdtrackdata.sdmomx) );
  fTree->Branch(  "SDTrack.momy", &(allsdtrackdata.sdmomy) );
  fTree->Branch(  "SDTrack.momz", &(allsdtrackdata.sdmomz) );
  fTree->Branch(  "SDTrack.polx", &(allsdtrackdata.sdpolx) );
  fTree->Branch(  "SDTrack.poly", &(allsdtrackdata.sdpoly) );
  fTree->Branch(  "SDTrack.polz", &(allsdtrackdata.sdpolz) );
  fTree->Branch(  "SDTrack.Etot", &(allsdtrackdata.sdenergy) );
  fTree->Branch(  "SDTrack.T", &(allsdtrackdata.sdtime) );
  //Add new vertex info:
  fTree->Branch(  "SDTrack.vx", &(allsdtrackdata.sdvx) );
  fTree->Branch(  "SDTrack.vy", &(allsdtrackdata.sdvy) );
  fTree->Branch(  "SDTrack.vz", &(allsdtrackdata.sdvz) );
  fTree->Branch(  "SDTrack.vnx", &(allsdtrackdata.sdvnx) );
  fTree->Branch(  "SDTrack.vny", &(allsdtrackdata.sdvny) );
  fTree->Branch(  "SDTrack.vnz", &(allsdtrackdata.sdvnz) );
  fTree->Branch(  "SDTrack.vEkin", &(allsdtrackdata.sdEkin) );
  //}
}

void G4SBSIO::BranchBD(G4String SDname){
   // create the branches for the Beam Diffuser (BD) 
   TString branch_name;
   TString branch_prefix = SDname.data();
   branch_prefix.ReplaceAll("/",".");
   // define branches
   fTree->Branch( branch_name.Format("%s.hit.nhits", branch_prefix.Data() ), &(BDdata[SDname].nhits_BD) );
   fTree->Branch( branch_name.Format("%s.hit.plane", branch_prefix.Data() ), &(BDdata[SDname].plane)    );
   fTree->Branch( branch_name.Format("%s.hit.trid" , branch_prefix.Data() ), &(BDdata[SDname].trid )    );
   fTree->Branch( branch_name.Format("%s.hit.mid"  , branch_prefix.Data() ), &(BDdata[SDname].mid  )    );
   fTree->Branch( branch_name.Format("%s.hit.pid"  , branch_prefix.Data() ), &(BDdata[SDname].pid  )    );
   fTree->Branch( branch_name.Format("%s.hit.x"    , branch_prefix.Data() ), &(BDdata[SDname].x    )    );
   fTree->Branch( branch_name.Format("%s.hit.y"    , branch_prefix.Data() ), &(BDdata[SDname].y    )    );
   fTree->Branch( branch_name.Format("%s.hit.z"    , branch_prefix.Data() ), &(BDdata[SDname].z    )    );
   fTree->Branch( branch_name.Format("%s.hit.t"    , branch_prefix.Data() ), &(BDdata[SDname].t    )    );
   fTree->Branch( branch_name.Format("%s.hit.xg"   , branch_prefix.Data() ), &(BDdata[SDname].xg   )    );
   fTree->Branch( branch_name.Format("%s.hit.yg"   , branch_prefix.Data() ), &(BDdata[SDname].yg   )    );
   fTree->Branch( branch_name.Format("%s.hit.zg"   , branch_prefix.Data() ), &(BDdata[SDname].zg   )    );
   fTree->Branch( branch_name.Format("%s.hit.p"    , branch_prefix.Data() ), &(BDdata[SDname].p    )    );
   fTree->Branch( branch_name.Format("%s.hit.edep" , branch_prefix.Data() ), &(BDdata[SDname].edep )    );
   fTree->Branch( branch_name.Format("%s.hit.beta" , branch_prefix.Data() ), &(BDdata[SDname].beta )    );
}

void G4SBSIO::BranchIC(G4String SDname){
   // create the branches for the Ion Chamber (IC)  
   TString branch_name;
   TString branch_prefix = SDname.data();
   branch_prefix.ReplaceAll("/",".");
   // define branches
   fTree->Branch( branch_name.Format("%s.hit.nhits", branch_prefix.Data() ), &(ICdata[SDname].nhits_IC) );
   fTree->Branch( branch_name.Format("%s.hit.trid" , branch_prefix.Data() ), &(ICdata[SDname].trid )    );
   fTree->Branch( branch_name.Format("%s.hit.mid"  , branch_prefix.Data() ), &(ICdata[SDname].mid  )    );
   fTree->Branch( branch_name.Format("%s.hit.pid"  , branch_prefix.Data() ), &(ICdata[SDname].pid  )    );
   fTree->Branch( branch_name.Format("%s.hit.x"    , branch_prefix.Data() ), &(ICdata[SDname].x    )    );
   fTree->Branch( branch_name.Format("%s.hit.y"    , branch_prefix.Data() ), &(ICdata[SDname].y    )    );
   fTree->Branch( branch_name.Format("%s.hit.z"    , branch_prefix.Data() ), &(ICdata[SDname].z    )    );
   fTree->Branch( branch_name.Format("%s.hit.t"    , branch_prefix.Data() ), &(ICdata[SDname].t    )    );
   fTree->Branch( branch_name.Format("%s.hit.xg"   , branch_prefix.Data() ), &(ICdata[SDname].xg   )    );
   fTree->Branch( branch_name.Format("%s.hit.yg"   , branch_prefix.Data() ), &(ICdata[SDname].yg   )    );
   fTree->Branch( branch_name.Format("%s.hit.zg"   , branch_prefix.Data() ), &(ICdata[SDname].zg   )    );
   fTree->Branch( branch_name.Format("%s.hit.p"    , branch_prefix.Data() ), &(ICdata[SDname].p    )    );
   fTree->Branch( branch_name.Format("%s.hit.edep" , branch_prefix.Data() ), &(ICdata[SDname].edep )    );
   fTree->Branch( branch_name.Format("%s.hit.beta" , branch_prefix.Data() ), &(ICdata[SDname].beta )    );

   map<G4String,G4bool>::iterator it = KeepHistoryflags.find( SDname );
   if( it != KeepHistoryflags.end() && it->second ){
      //Branches with "Particle History" data:
      fTree->Branch( branch_name.Format("%s.part.npart", branch_prefix.Data() ), &(ICdata[SDname].ParticleHistory.npart) );
      fTree->Branch( branch_name.Format("%s.part.PID"  , branch_prefix.Data() ), &(ICdata[SDname].ParticleHistory.PID) );
      fTree->Branch( branch_name.Format("%s.part.MID"  , branch_prefix.Data() ), &(ICdata[SDname].ParticleHistory.MID) );
      fTree->Branch( branch_name.Format("%s.part.TID"  , branch_prefix.Data() ), &(ICdata[SDname].ParticleHistory.TID) );
      fTree->Branch( branch_name.Format("%s.part.vx"   , branch_prefix.Data() ), &(ICdata[SDname].ParticleHistory.vx) );
      fTree->Branch( branch_name.Format("%s.part.vy"   , branch_prefix.Data() ), &(ICdata[SDname].ParticleHistory.vy) );
      fTree->Branch( branch_name.Format("%s.part.vz"   , branch_prefix.Data() ), &(ICdata[SDname].ParticleHistory.vz) );
      fTree->Branch( branch_name.Format("%s.part.px"   , branch_prefix.Data() ), &(ICdata[SDname].ParticleHistory.px) );
      fTree->Branch( branch_name.Format("%s.part.py"   , branch_prefix.Data() ), &(ICdata[SDname].ParticleHistory.py) );
      fTree->Branch( branch_name.Format("%s.part.pz"   , branch_prefix.Data() ), &(ICdata[SDname].ParticleHistory.pz) );
      // fTree->Branch( branch_name.Format("%s.part.nbounce", branch_prefix.Data() ), &(ICdata[SDname].ParticleHistory.nbounce) );
      // fTree->Branch( branch_name.Format("%s.part.hitindex", branch_prefix.Data() ), &(ICdata[SDname].ParticleHistory.hitindex) );
      // fTree->Branch( branch_name.Format("%s.part.Nphe_part", branch_prefix.Data() ), &(ICdata[SDname].Nphe_part) );
   }
}

void G4SBSIO::BranchGEnTarget_Glass(G4String SDname){
   // create the branches for the GEn target glass cell 
   TString branch_name;
   TString branch_prefix = SDname.data();
   branch_prefix.ReplaceAll("/",".");
   // define branches
   fTree->Branch( branch_name.Format("%s.hit.nhits", branch_prefix.Data() ), &(genTgtGCdata[SDname].nhits_Target) );
   fTree->Branch( branch_name.Format("%s.hit.trid" , branch_prefix.Data() ), &(genTgtGCdata[SDname].trid )    );
   fTree->Branch( branch_name.Format("%s.hit.mid"  , branch_prefix.Data() ), &(genTgtGCdata[SDname].mid  )    );
   fTree->Branch( branch_name.Format("%s.hit.pid"  , branch_prefix.Data() ), &(genTgtGCdata[SDname].pid  )    );
   // fTree->Branch( branch_name.Format("%s.hit.x"    , branch_prefix.Data() ), &(genTgtGCdata[SDname].x    )    );
   // fTree->Branch( branch_name.Format("%s.hit.y"    , branch_prefix.Data() ), &(genTgtGCdata[SDname].y    )    );
   // fTree->Branch( branch_name.Format("%s.hit.z"    , branch_prefix.Data() ), &(genTgtGCdata[SDname].z    )    );
   // fTree->Branch( branch_name.Format("%s.hit.t"    , branch_prefix.Data() ), &(genTgtGCdata[SDname].t    )    );
   // fTree->Branch( branch_name.Format("%s.hit.xg"   , branch_prefix.Data() ), &(genTgtGCdata[SDname].xg   )    );
   // fTree->Branch( branch_name.Format("%s.hit.yg"   , branch_prefix.Data() ), &(genTgtGCdata[SDname].yg   )    );
   // fTree->Branch( branch_name.Format("%s.hit.zg"   , branch_prefix.Data() ), &(genTgtGCdata[SDname].zg   )    );
   fTree->Branch( branch_name.Format("%s.hit.p"    , branch_prefix.Data() ), &(genTgtGCdata[SDname].p    )    );
   fTree->Branch( branch_name.Format("%s.hit.edep" , branch_prefix.Data() ), &(genTgtGCdata[SDname].edep )    );
   fTree->Branch( branch_name.Format("%s.hit.beta" , branch_prefix.Data() ), &(genTgtGCdata[SDname].beta )    );
   fTree->Branch( branch_name.Format("%s.hit.trackLength" , branch_prefix.Data() ), &(genTgtGCdata[SDname].beta )    );

   map<G4String,G4bool>::iterator it = KeepHistoryflags.find( SDname );
   if( it != KeepHistoryflags.end() && it->second ){
      //Branches with "Particle History" data:
      fTree->Branch( branch_name.Format("%s.part.npart", branch_prefix.Data() ), &(genTgtGCdata[SDname].ParticleHistory.npart) );
      fTree->Branch( branch_name.Format("%s.part.PID"  , branch_prefix.Data() ), &(genTgtGCdata[SDname].ParticleHistory.PID) );
      fTree->Branch( branch_name.Format("%s.part.MID"  , branch_prefix.Data() ), &(genTgtGCdata[SDname].ParticleHistory.MID) );
      fTree->Branch( branch_name.Format("%s.part.TID"  , branch_prefix.Data() ), &(genTgtGCdata[SDname].ParticleHistory.TID) );
      // fTree->Branch( branch_name.Format("%s.part.vx"   , branch_prefix.Data() ), &(genTgtGCdata[SDname].ParticleHistory.vx) );
      // fTree->Branch( branch_name.Format("%s.part.vy"   , branch_prefix.Data() ), &(genTgtGCdata[SDname].ParticleHistory.vy) );
      // fTree->Branch( branch_name.Format("%s.part.vz"   , branch_prefix.Data() ), &(genTgtGCdata[SDname].ParticleHistory.vz) );
      fTree->Branch( branch_name.Format("%s.part.px"   , branch_prefix.Data() ), &(genTgtGCdata[SDname].ParticleHistory.px) );
      fTree->Branch( branch_name.Format("%s.part.py"   , branch_prefix.Data() ), &(genTgtGCdata[SDname].ParticleHistory.py) );
      fTree->Branch( branch_name.Format("%s.part.pz"   , branch_prefix.Data() ), &(genTgtGCdata[SDname].ParticleHistory.pz) );
   }
}

void G4SBSIO::BranchGEnTarget_Al(G4String SDname){
   // create the branches for the GEn target glass cell, endcap (Al or Cu) 
   TString branch_name;
   TString branch_prefix = SDname.data();
   branch_prefix.ReplaceAll("/",".");
   // define branches
   fTree->Branch( branch_name.Format("%s.hit.nhits", branch_prefix.Data() ), &(genTgtALdata[SDname].nhits_Target) );
   fTree->Branch( branch_name.Format("%s.hit.trid" , branch_prefix.Data() ), &(genTgtALdata[SDname].trid )    );
   fTree->Branch( branch_name.Format("%s.hit.mid"  , branch_prefix.Data() ), &(genTgtALdata[SDname].mid  )    );
   fTree->Branch( branch_name.Format("%s.hit.pid"  , branch_prefix.Data() ), &(genTgtALdata[SDname].pid  )    );
   // fTree->Branch( branch_name.Format("%s.hit.x"    , branch_prefix.Data() ), &(genTgtALdata[SDname].x    )    );
   // fTree->Branch( branch_name.Format("%s.hit.y"    , branch_prefix.Data() ), &(genTgtALdata[SDname].y    )    );
   // fTree->Branch( branch_name.Format("%s.hit.z"    , branch_prefix.Data() ), &(genTgtALdata[SDname].z    )    );
   // fTree->Branch( branch_name.Format("%s.hit.t"    , branch_prefix.Data() ), &(genTgtALdata[SDname].t    )    );
   // fTree->Branch( branch_name.Format("%s.hit.xg"   , branch_prefix.Data() ), &(genTgtALdata[SDname].xg   )    );
   // fTree->Branch( branch_name.Format("%s.hit.yg"   , branch_prefix.Data() ), &(genTgtALdata[SDname].yg   )    );
   // fTree->Branch( branch_name.Format("%s.hit.zg"   , branch_prefix.Data() ), &(genTgtALdata[SDname].zg   )    );
   fTree->Branch( branch_name.Format("%s.hit.p"    , branch_prefix.Data() ), &(genTgtALdata[SDname].p    )    );
   fTree->Branch( branch_name.Format("%s.hit.edep" , branch_prefix.Data() ), &(genTgtALdata[SDname].edep )    );
   fTree->Branch( branch_name.Format("%s.hit.beta" , branch_prefix.Data() ), &(genTgtALdata[SDname].beta )    );
   fTree->Branch( branch_name.Format("%s.hit.trackLength" , branch_prefix.Data() ), &(genTgtALdata[SDname].trackLength )    );

   map<G4String,G4bool>::iterator it = KeepHistoryflags.find( SDname );
   if( it != KeepHistoryflags.end() && it->second ){
      //Branches with "Particle History" data:
      fTree->Branch( branch_name.Format("%s.part.npart", branch_prefix.Data() ), &(genTgtALdata[SDname].ParticleHistory.npart) );
      fTree->Branch( branch_name.Format("%s.part.PID"  , branch_prefix.Data() ), &(genTgtALdata[SDname].ParticleHistory.PID) );
      fTree->Branch( branch_name.Format("%s.part.MID"  , branch_prefix.Data() ), &(genTgtALdata[SDname].ParticleHistory.MID) );
      fTree->Branch( branch_name.Format("%s.part.TID"  , branch_prefix.Data() ), &(genTgtALdata[SDname].ParticleHistory.TID) );
      // fTree->Branch( branch_name.Format("%s.part.vx"   , branch_prefix.Data() ), &(genTgtALdata[SDname].ParticleHistory.vx) );
      // fTree->Branch( branch_name.Format("%s.part.vy"   , branch_prefix.Data() ), &(genTgtALdata[SDname].ParticleHistory.vy) );
      // fTree->Branch( branch_name.Format("%s.part.vz"   , branch_prefix.Data() ), &(genTgtALdata[SDname].ParticleHistory.vz) );
      fTree->Branch( branch_name.Format("%s.part.px"   , branch_prefix.Data() ), &(genTgtALdata[SDname].ParticleHistory.px) );
      fTree->Branch( branch_name.Format("%s.part.py"   , branch_prefix.Data() ), &(genTgtALdata[SDname].ParticleHistory.py) );
      fTree->Branch( branch_name.Format("%s.part.pz"   , branch_prefix.Data() ), &(genTgtALdata[SDname].ParticleHistory.pz) );
   }
}

void G4SBSIO::BranchGEnTarget_Cu(G4String SDname){
   // create the branches for the GEn target glass cell, endcap (Al or Cu) 
   TString branch_name;
   TString branch_prefix = SDname.data();
   branch_prefix.ReplaceAll("/",".");
   // define branches
   fTree->Branch( branch_name.Format("%s.hit.nhits", branch_prefix.Data() ), &(genTgtCUdata[SDname].nhits_Target) );
   fTree->Branch( branch_name.Format("%s.hit.trid" , branch_prefix.Data() ), &(genTgtCUdata[SDname].trid )    );
   fTree->Branch( branch_name.Format("%s.hit.mid"  , branch_prefix.Data() ), &(genTgtCUdata[SDname].mid  )    );
   fTree->Branch( branch_name.Format("%s.hit.pid"  , branch_prefix.Data() ), &(genTgtCUdata[SDname].pid  )    );
   // fTree->Branch( branch_name.Format("%s.hit.x"    , branch_prefix.Data() ), &(genTgtCUdata[SDname].x    )    );
   // fTree->Branch( branch_name.Format("%s.hit.y"    , branch_prefix.Data() ), &(genTgtCUdata[SDname].y    )    );
   // fTree->Branch( branch_name.Format("%s.hit.z"    , branch_prefix.Data() ), &(genTgtCUdata[SDname].z    )    );
   // fTree->Branch( branch_name.Format("%s.hit.t"    , branch_prefix.Data() ), &(genTgtCUdata[SDname].t    )    );
   // fTree->Branch( branch_name.Format("%s.hit.xg"   , branch_prefix.Data() ), &(genTgtCUdata[SDname].xg   )    );
   // fTree->Branch( branch_name.Format("%s.hit.yg"   , branch_prefix.Data() ), &(genTgtCUdata[SDname].yg   )    );
   // fTree->Branch( branch_name.Format("%s.hit.zg"   , branch_prefix.Data() ), &(genTgtCUdata[SDname].zg   )    );
   fTree->Branch( branch_name.Format("%s.hit.p"    , branch_prefix.Data() ), &(genTgtCUdata[SDname].p    )    );
   fTree->Branch( branch_name.Format("%s.hit.edep" , branch_prefix.Data() ), &(genTgtCUdata[SDname].edep )    );
   fTree->Branch( branch_name.Format("%s.hit.beta" , branch_prefix.Data() ), &(genTgtCUdata[SDname].beta )    );
   fTree->Branch( branch_name.Format("%s.hit.trackLength" , branch_prefix.Data() ), &(genTgtCUdata[SDname].trackLength )    );

   map<G4String,G4bool>::iterator it = KeepHistoryflags.find( SDname );
   if( it != KeepHistoryflags.end() && it->second ){
      //Branches with "Particle History" data:
      fTree->Branch( branch_name.Format("%s.part.npart", branch_prefix.Data() ), &(genTgtCUdata[SDname].ParticleHistory.npart) );
      fTree->Branch( branch_name.Format("%s.part.PID"  , branch_prefix.Data() ), &(genTgtCUdata[SDname].ParticleHistory.PID) );
      fTree->Branch( branch_name.Format("%s.part.MID"  , branch_prefix.Data() ), &(genTgtCUdata[SDname].ParticleHistory.MID) );
      fTree->Branch( branch_name.Format("%s.part.TID"  , branch_prefix.Data() ), &(genTgtCUdata[SDname].ParticleHistory.TID) );
      // fTree->Branch( branch_name.Format("%s.part.vx"   , branch_prefix.Data() ), &(genTgtCUdata[SDname].ParticleHistory.vx) );
      // fTree->Branch( branch_name.Format("%s.part.vy"   , branch_prefix.Data() ), &(genTgtCUdata[SDname].ParticleHistory.vy) );
      // fTree->Branch( branch_name.Format("%s.part.vz"   , branch_prefix.Data() ), &(genTgtCUdata[SDname].ParticleHistory.vz) );
      fTree->Branch( branch_name.Format("%s.part.px"   , branch_prefix.Data() ), &(genTgtCUdata[SDname].ParticleHistory.px) );
      fTree->Branch( branch_name.Format("%s.part.py"   , branch_prefix.Data() ), &(genTgtCUdata[SDname].ParticleHistory.py) );
      fTree->Branch( branch_name.Format("%s.part.pz"   , branch_prefix.Data() ), &(genTgtCUdata[SDname].ParticleHistory.pz) );
   }
}

void G4SBSIO::BranchGEnTarget_3He(G4String SDname){
   // create the branches for the GEn target glass cell, endcap (Al or Cu) 
   TString branch_name;
   TString branch_prefix = SDname.data();
   branch_prefix.ReplaceAll("/",".");
   // define branches
   fTree->Branch( branch_name.Format("%s.hit.nhits", branch_prefix.Data() ), &(genTgt3HEdata[SDname].nhits_Target) );
   fTree->Branch( branch_name.Format("%s.hit.trid" , branch_prefix.Data() ), &(genTgt3HEdata[SDname].trid )    );
   fTree->Branch( branch_name.Format("%s.hit.mid"  , branch_prefix.Data() ), &(genTgt3HEdata[SDname].mid  )    );
   fTree->Branch( branch_name.Format("%s.hit.pid"  , branch_prefix.Data() ), &(genTgt3HEdata[SDname].pid  )    );
   // fTree->Branch( branch_name.Format("%s.hit.x"    , branch_prefix.Data() ), &(genTgt3HEdata[SDname].x    )    );
   // fTree->Branch( branch_name.Format("%s.hit.y"    , branch_prefix.Data() ), &(genTgt3HEdata[SDname].y    )    );
   // fTree->Branch( branch_name.Format("%s.hit.z"    , branch_prefix.Data() ), &(genTgt3HEdata[SDname].z    )    );
   // fTree->Branch( branch_name.Format("%s.hit.t"    , branch_prefix.Data() ), &(genTgt3HEdata[SDname].t    )    );
   // fTree->Branch( branch_name.Format("%s.hit.xg"   , branch_prefix.Data() ), &(genTgt3HEdata[SDname].xg   )    );
   // fTree->Branch( branch_name.Format("%s.hit.yg"   , branch_prefix.Data() ), &(genTgt3HEdata[SDname].yg   )    );
   // fTree->Branch( branch_name.Format("%s.hit.zg"   , branch_prefix.Data() ), &(genTgt3HEdata[SDname].zg   )    );
   fTree->Branch( branch_name.Format("%s.hit.p"    , branch_prefix.Data() ), &(genTgt3HEdata[SDname].p    )    );
   fTree->Branch( branch_name.Format("%s.hit.edep" , branch_prefix.Data() ), &(genTgt3HEdata[SDname].edep )    );
   fTree->Branch( branch_name.Format("%s.hit.beta" , branch_prefix.Data() ), &(genTgt3HEdata[SDname].beta )    );
   fTree->Branch( branch_name.Format("%s.hit.trackLength" , branch_prefix.Data() ), &(genTgt3HEdata[SDname].trackLength )    );

   map<G4String,G4bool>::iterator it = KeepHistoryflags.find( SDname );
   if( it != KeepHistoryflags.end() && it->second ){
      //Branches with "Particle History" data:
      fTree->Branch( branch_name.Format("%s.part.npart", branch_prefix.Data() ), &(genTgt3HEdata[SDname].ParticleHistory.npart) );
      fTree->Branch( branch_name.Format("%s.part.PID"  , branch_prefix.Data() ), &(genTgt3HEdata[SDname].ParticleHistory.PID) );
      fTree->Branch( branch_name.Format("%s.part.MID"  , branch_prefix.Data() ), &(genTgt3HEdata[SDname].ParticleHistory.MID) );
      fTree->Branch( branch_name.Format("%s.part.TID"  , branch_prefix.Data() ), &(genTgt3HEdata[SDname].ParticleHistory.TID) );
      // fTree->Branch( branch_name.Format("%s.part.vx"   , branch_prefix.Data() ), &(genTgt3HEdata[SDname].ParticleHistory.vx) );
      // fTree->Branch( branch_name.Format("%s.part.vy"   , branch_prefix.Data() ), &(genTgt3HEdata[SDname].ParticleHistory.vy) );
      // fTree->Branch( branch_name.Format("%s.part.vz"   , branch_prefix.Data() ), &(genTgt3HEdata[SDname].ParticleHistory.vz) );
      fTree->Branch( branch_name.Format("%s.part.px"   , branch_prefix.Data() ), &(genTgt3HEdata[SDname].ParticleHistory.px) );
      fTree->Branch( branch_name.Format("%s.part.py"   , branch_prefix.Data() ), &(genTgt3HEdata[SDname].ParticleHistory.py) );
      fTree->Branch( branch_name.Format("%s.part.pz"   , branch_prefix.Data() ), &(genTgt3HEdata[SDname].ParticleHistory.pz) );
   }
}

