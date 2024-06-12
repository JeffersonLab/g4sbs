#include "TMatrixD.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TDecompSVD.h"
#include "TCut.h"
#include "TEventList.h"
//#include "g4sbs_tree.C"
//#include "gep_tree_with_spin.C"
//#include "genrp_tree.C"
//#include "gmn_tree.C"
#include "genrp_tree.C"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
//#include "a1n_tree.C"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TF2.h"
#include "TChainElement.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include "TObjArray.h"
#include "TObjString.h"
#include "G4SBSRunData.hh"
#include "G4SBSTextFile.hh"


const double PI = TMath::Pi();

//Optics fit code for GMN/GEN-RP

void Optics_GENRP( const char *inputfilename, const char *outputfilename, int NMAX=1000000, double BBscale=0.97, int downbend=0){
  
  
  ifstream inputfile(inputfilename);
  TFile *fout = new TFile( outputfilename, "RECREATE" );
  TChain *C = new TChain("T");

  set<TString> files;
  
  TString opticsfilename = outputfilename;
  opticsfilename.ReplaceAll(".root",".txt");  
  ofstream opticsfile(opticsfilename);

  TString TOFfilename = outputfilename;
  TOFfilename.ReplaceAll(".root",".dat");
  TOFfilename.Prepend("TOF_");
  ofstream TOFfile(TOFfilename);
  
  opticsfilename.Prepend("foptics_");
  ofstream fopticsfile(opticsfilename);
  
  TString currentline;
  while( currentline.ReadLine(inputfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());

      files.insert( currentline );
    }
  }

  TObjArray *FileList = C->GetListOfFiles();
  TIter next(FileList);

  TChainElement *chEl = 0;

  set<TString> bad_file_list;

  long ntries=0;
  long ngen=0;
  
  G4SBSRunData *rd;

  map<TString,double> BBdist_file;
  map<TString,double> BBang_file;
  map<TString,double> SBSang_file;
  map<TString,double> SBSdist_file;

  //map<TString,double> BBscale_file; 
  
  double A_pth1 = 0.28615*BBscale;
  double B_pth1 = 0.1976 + 0.4764 * 1.55;

  double BBdist_default = 1.63;
  
  while( (chEl=(TChainElement*)next() )){
    TFile newfile(chEl->GetTitle());
    newfile.GetObject("run_data",rd);
    if( rd ){
      //if( files.find( chEl->GetTitle() ) != files.end() ){
      ngen += rd->fNthrown;
      ntries += rd->fNtries;
      // ngen_file[chEl->GetTitle()] = rd->fNthrown;
      // ntries_file[chEl->GetTitle()] = rd->fNtries;
      BBang_file[chEl->GetTitle()] = rd->fBBtheta;
      SBSang_file[chEl->GetTitle()] = rd->fSBStheta;
      BBdist_file[chEl->GetTitle()] = rd->fBBdist;
      SBSdist_file[chEl->GetTitle()] = rd->fSBSdist;

      B_pth1 = 0.1976 + 0.4764 * rd->fBBdist;
      
      cout << "Bigbite theta = " << BBang_file[chEl->GetTitle()]*180.0/PI << " degrees" << endl;
      cout << "BigBite distance = " << BBdist_file[chEl->GetTitle()] << " meters" << endl;
      cout << "SBS theta = " << SBSang_file[chEl->GetTitle()]*180.0/PI << " degrees" << endl;
      cout << "SBS distance = " << SBSdist_file[chEl->GetTitle()] << " meters" << endl;

      BBdist_default = rd->fBBdist;
      
    } else {
      bad_file_list.insert( chEl->GetTitle());
    }
  }
  
  cout << "Total number of generated events = " << ngen << endl;
  cout << "Total number of attempts = " << ntries << endl;
  cout << "Generation efficiency = " << double(ngen)/double(ntries) << endl;
  //cout << "Total number of generated events, neutron = " << ngen_n << endl;
  
  TCut global_cut = "";
  while( currentline.ReadLine(inputfile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      global_cut += currentline.Data();
    }
  }

  if( downbend != 0 ){
    A_pth1 = 0.2885 * BBscale/0.9644;
    B_pth1 = 0.994; //to be refined later
  }
  
  TEventList *elist = new TEventList("elist");

  C->Draw(">>elist",global_cut);

  cout << "Number of events passing global cut = " << elist->GetN() << endl;

  // genrp_tree *T = new genrp_tree(C);
  genrp_tree *T = new genrp_tree(C);

  //Expansion is of the form: 
  // ( xptar yptar 1/p ytar ) = sum_{i+j+k+l+m<=N} C_{ijklm} xfp^i yfp^j xpfp^k ypfp^l xtar^m 
  // Fit chi^2 is sum_n=1^{Nevent} (xptar_n - sum_ijklm C_{ijklm} ...)^2, and similarly for yptar, ytar, 1/p, etc. 

  int order=4; 
  inputfile >> order; 
  int nbins = 200;
  inputfile >> nbins;
  int arm;  // 0 = BB, 1 = SBS
  inputfile >> arm;
  double chi2cut; 
  inputfile >> chi2cut; 
  double tracker_pitch_angle = 0.0;
  inputfile >> tracker_pitch_angle; //Positive value means that for up-bending particles, thetabend = tracker_pitch + thetatgt - thetafp  
  tracker_pitch_angle *= PI/180.0; //value is assumed to be given in degrees
  
  //double pcentral[2] = {2.5,5.0};

  double fpmin[4], fpmax[4];    //ranges for plotting of fp variables: order is: x, y, xp, yp
  double tgtmin[5], tgtmax[5];  //ranges for plotting of tgt variables: order is: xp, yp, y, p, x

  for(int i=0; i<4; i++){
    inputfile >> fpmin[i] >> fpmax[i];
  }
  for(int i=0; i<5; i++){
    inputfile >> tgtmin[i] >> tgtmax[i];
  }

  int xtar_flag=1; //include and fit xtar-dependent terms in the expansion or not? (also interpret flag as maximum order for xtar-dependent terms)
  inputfile >> xtar_flag;

  int pexpansion_flag=0; //default = p*thetabend expansion, alternative = 1/p expansion

  inputfile >> pexpansion_flag;

  double cutxptar_true = 0.2;
  double cutyptar_true = 0.1;
  double cutytar_true = 0.1;
  inputfile >> cutxptar_true;
  inputfile >> cutyptar_true;
  inputfile >> cutytar_true; 

  double dsieve=1.18; //m
  inputfile >> dsieve;

  
  
  fout->cd();
  
  //Histograms of "true" versions of fp track parameters:
  TH1D *hxfp = new TH1D("hxfp", "", nbins, fpmin[0], fpmax[0] );
  TH1D *hyfp = new TH1D("hyfp", "", nbins, fpmin[1], fpmax[1] );
  TH1D *hxpfp = new TH1D("hxpfp", "", nbins, fpmin[2], fpmax[2] );
  TH1D *hypfp = new TH1D("hypfp", "", nbins, fpmin[3], fpmax[3] );
  
  //Histograms of "recon" versions of fp track parameters:
  TH1D *hxfprecon = new TH1D("hxfprecon", "", nbins, fpmin[0], fpmax[0] );
  TH1D *hyfprecon = new TH1D("hyfprecon", "", nbins, fpmin[1], fpmax[1] );
  TH1D *hxpfprecon = new TH1D("hxpfprecon", "", nbins, fpmin[2], fpmax[2] );
  TH1D *hypfprecon = new TH1D("hypfprecon", "", nbins, fpmin[3], fpmax[3] );

  //Histograms of "recon" versions of fp track parameters:
  TH1D *hxfpdiff = new TH1D("hxfpdiff", "", nbins, -0.5e-3, 0.5e-3);
  TH1D *hyfpdiff = new TH1D("hyfpdiff", "", nbins, -0.5e-3, 0.5e-3);
  TH1D *hxpfpdiff = new TH1D("hxpfpdiff", "", nbins, -0.005, 0.005 );
  TH1D *hypfpdiff = new TH1D("hypfpdiff", "", nbins, -0.005, 0.005 );

  //Histograms of "true" versions of tgt track parameters:
  TH1D *hvz   = new TH1D("hvz","", nbins, -0.4, 0.4 );
  TH1D *hxtar = new TH1D("hxtar", "", nbins, tgtmin[4], tgtmax[4] );
  TH1D *hytar = new TH1D("hytar", "", nbins, tgtmin[2], tgtmax[2] );
  TH1D *hxptar = new TH1D("hxptar", "", nbins, tgtmin[0], tgtmax[0] );
  TH1D *hyptar = new TH1D("hyptar", "", nbins, tgtmin[1], tgtmax[1] );
  TH1D *hp     = new TH1D("hp", "", nbins, tgtmin[3], tgtmax[3] );

  //Histograms of "recon" versions of tgt track parameters:
  TH1D *hxtarrecon = new TH1D("hxtarrecon", "", nbins, tgtmin[4], tgtmax[4] );
  TH1D *hytarrecon = new TH1D("hytarrecon", "", nbins, tgtmin[2], tgtmax[2] );
  TH1D *hxptarrecon = new TH1D("hxptarrecon", "", nbins, tgtmin[0], tgtmax[0] );
  TH1D *hyptarrecon = new TH1D("hyptarrecon", "", nbins, tgtmin[1], tgtmax[1] );
  TH1D *hprecon     = new TH1D("hprecon", "", nbins, tgtmin[3], tgtmax[3] );
  TH1D *hvzrecon    = new TH1D("hvzrecon","",nbins,-0.4,0.4);

  TH1D *hxtardiff = new TH1D("hxtardiff", "", nbins, tgtmin[4], tgtmax[4] );
  TH1D *hytardiff = new TH1D("hytardiff", "", nbins, -0.05, 0.05 );
  TH1D *hxptardiff = new TH1D("hxptardiff", "", nbins, -0.02, 0.02 );
  TH1D *hyptardiff = new TH1D("hyptardiff", "", nbins, -0.04, 0.04 );
  TH1D *hpdiff     = new TH1D("hpdiff", "", nbins, -0.2,0.2 );

  TH2D *hxptardiff_xtar = new TH2D("hxptardiff_xtar", "", nbins, tgtmin[4], tgtmax[4], nbins, -0.02, 0.02 );
  TH2D *hxptardiff_xptar = new TH2D("hxptardiff_xptar", "", nbins, tgtmin[0], tgtmax[0], nbins, -0.02, 0.02 );
  TH2D *hxptardiff_yptar = new TH2D("hxptardiff_yptar", "", nbins, tgtmin[1], tgtmax[1], nbins, -0.02, 0.02 );
  TH2D *hxptardiff_ytar = new TH2D("hxptardiff_ytar", "", nbins, tgtmin[2], tgtmax[2], nbins, -0.02, 0.02 );
  TH2D *hxptardiff_p = new TH2D("hxptardiff_p", "", nbins, tgtmin[3], tgtmax[3], nbins, -0.02, 0.02 );

  TH2D *hyptardiff_xtar = new TH2D("hyptardiff_xtar", "", nbins, tgtmin[4], tgtmax[4], nbins, -0.04, 0.04 );
  TH2D *hyptardiff_xptar = new TH2D("hyptardiff_xptar", "", nbins, tgtmin[0], tgtmax[0], nbins, -0.04, 0.04 );
  TH2D *hyptardiff_yptar = new TH2D("hyptardiff_yptar", "", nbins, tgtmin[1], tgtmax[1], nbins, -0.04, 0.04 );
  TH2D *hyptardiff_ytar = new TH2D("hyptardiff_ytar", "", nbins, tgtmin[2], tgtmax[2], nbins, -0.04, 0.04 );
  TH2D *hyptardiff_p = new TH2D("hyptardiff_p", "", nbins, tgtmin[3], tgtmax[3], nbins, -0.04, 0.04 );

  TH2D *hytardiff_xtar = new TH2D("hytardiff_xtar", "", nbins, tgtmin[4], tgtmax[4], nbins, -0.05, 0.05 );
  TH2D *hytardiff_xptar = new TH2D("hytardiff_xptar", "", nbins, tgtmin[0], tgtmax[0], nbins, -0.05, 0.05 );
  TH2D *hytardiff_yptar = new TH2D("hytardiff_yptar", "", nbins, tgtmin[1], tgtmax[1], nbins, -0.05, 0.05 );
  TH2D *hytardiff_ytar = new TH2D("hytardiff_ytar", "", nbins, tgtmin[2], tgtmax[2], nbins, -0.05, 0.05 );
  TH2D *hytardiff_p = new TH2D("hytardiff_p", "", nbins, tgtmin[3], tgtmax[3], nbins, -0.05, 0.05 );
  
  TH2D *hpdiff_xtar = new TH2D("hpdiff_xtar", "", nbins, tgtmin[4], tgtmax[4], nbins, -0.2, 0.2 );
  TH2D *hpdiff_xptar = new TH2D("hpdiff_xptar", "", nbins, tgtmin[0], tgtmax[0], nbins, -0.2, 0.2 );
  TH2D *hpdiff_yptar = new TH2D("hpdiff_yptar", "", nbins, tgtmin[1], tgtmax[1], nbins, -0.2, 0.2 );
  TH2D *hpdiff_ytar = new TH2D("hpdiff_ytar", "", nbins, tgtmin[2], tgtmax[2], nbins, -0.2, 0.2 );
  TH2D *hpdiff_p = new TH2D("hpdiff_p", "", nbins, tgtmin[3], tgtmax[3], nbins, -0.2, 0.2 );
  
  TH2D *hpdiff_xfp = new TH2D("hpdiff_xfp", "", nbins, fpmin[0], fpmax[0], nbins, -0.2, 0.2 );
  TH2D *hpdiff_yfp = new TH2D("hpdiff_yfp", "", nbins, fpmin[1], fpmax[1], nbins, -0.2, 0.2 );
  TH2D *hpdiff_xpfp = new TH2D("hpdiff_xpfp", "", nbins, fpmin[2], fpmax[2], nbins, -0.2, 0.2 );
  TH2D *hpdiff_ypfp = new TH2D("hpdiff_ypfp", "", nbins, fpmin[3], fpmax[3], nbins, -0.2, 0.2 );

  TH1D *hpinvdiff = new TH1D("hpinvdiff", "", nbins, -0.05, 0.05 );

  TH1D *hvzdiff = new TH1D("hvzdiff", "", nbins, -0.3, 0.3 );

  TH2D *hvzdiff_xtar = new TH2D("hvzdiff_xtar", "", nbins, tgtmin[4],tgtmax[4], nbins, -0.3,0.3);
  TH2D *hvzdiff_ytar = new TH2D("hvzdiff_ytar", "", nbins, tgtmin[2],tgtmax[2], nbins, -0.3,0.3);
  TH2D *hvzdiff_xptar = new TH2D("hvzdiff_xptar", "", nbins, tgtmin[0],tgtmax[0], nbins, -0.3,0.3);
  TH2D *hvzdiff_yptar = new TH2D("hvzdiff_yptar", "", nbins, tgtmin[1],tgtmax[1], nbins, -0.3,0.3);
  TH2D *hvzdiff_p = new TH2D("hvzdiff_p", "", nbins, tgtmin[3],tgtmax[3], nbins, -0.3,0.3);

  TH2D *hxysieve = new TH2D("hxysieve","",250,-0.15,0.15,250,-0.35,0.35);

  // Diagnostic histos for forward optics
  TH1D *hxfpdiff_foptics = new TH1D("hxfpdiff_foptics","",250,-0.02,0.02);
  TH1D *hyfpdiff_foptics = new TH1D("hyfpdiff_foptics","",250,-0.02,0.02);
  TH1D *hxpfpdiff_foptics = new TH1D("hxpfpdiff_foptics","",250,-0.01,0.01);
  TH1D *hypfpdiff_foptics = new TH1D("hypfpdiff_foptics","",250,-0.01,0.01);

  TH2D *hxfp_recon_vs_true_foptics = new TH2D("hxfp_recon_vs_true_foptics","",250,-0.8,0.8,250,-0.8,0.8);
  TH2D *hyfp_recon_vs_true_foptics = new TH2D("hyfp_recon_vs_true_foptics","",250,-0.25,0.25,250,-0.25,0.25);
  TH2D *hxpfp_recon_vs_true_foptics = new TH2D("hxpfp_recon_vs_true_foptics","",250,-0.4,0.4,250,-0.4,0.4);
  TH2D *hypfp_recon_vs_true_foptics = new TH2D("hypfp_recon_vs_true_foptics","",250,-0.1,0.1,250,-0.1,0.1);

  
  
  int nparams = 0;

  vector<int> xtar_expon;
  vector<int> xfp_expon;
  vector<int> xpfp_expon;
  vector<int> yfp_expon;
  vector<int> ypfp_expon;

  for(int i=0; i<=order; i++){
    for(int j=0; j<=order-i; j++){
      for(int k=0; k<=order-i-j; k++){
	for(int l=0; l<=order-i-j-k; l++){
	  for(int m=0; m<=order-i-j-k-l; m++){
	    nparams++;
	    xtar_expon.push_back( i );
	    xfp_expon.push_back( m );
	    yfp_expon.push_back( l );
	    xpfp_expon.push_back( k );
	    ypfp_expon.push_back( j );
	  }
	}
      }
    }
  }

  //TMatrixD M_xptar(nparams,nparams), M_yptar(nparams,nparams), M_ytar(nparams,nparams), M_pinv(nparams,nparams);
  TMatrixD M(nparams,nparams), Mforward(nparams,nparams);
  TVectorD b_xptar(nparams), b_yptar(nparams), b_ytar(nparams), b_pinv(nparams);
  TVectorD b_xfp(nparams), b_yfp(nparams), b_xpfp(nparams), b_ypfp(nparams);

  //Set up equations for electron time-of-flight (really path length)
  TMatrixD M_etof(nparams,nparams); //technically we don't need another M for etof it is the same as for the optics (we're just adding another "b").
  TVectorD b_etof(nparams);
  //What is our expansion formalism for electron TOF in BigBite?
  //Should we expand in terms of target variables or focal plane variables?
  //We want etof = t0avg + \sum_{i,j,k,l,m=1}^i+j+k+l+m<=order 
  
  for(int i=0; i<nparams; i++){
    for(int j=0; j<nparams; j++){
      M(i,j) = 0.0;
      Mforward(i,j) = 0.0;
      M_etof(i,j) = 0.0;
    }
    b_xptar(i) = 0.0;
    b_yptar(i) = 0.0;
    b_ytar(i) = 0.0;
    b_pinv(i) = 0.0;

    b_xfp(i) = 0.0;
    b_yfp(i) = 0.0;
    b_xpfp(i) = 0.0;
    b_ypfp(i) = 0.0;

    b_etof(i) = 0.0;
  }

 
  
  long nevent = 0;
  while( T->GetEntry(elist->GetEntry(nevent++)) && nevent < NMAX ){
    if( nevent%1000 == 0 ){
      cout << nevent << endl;
    }

    TString fname = C->GetFile()->GetName();
    
    if( bad_file_list.find( fname ) == bad_file_list.end() ){
      
      double R0; 
      double theta0;
      if( arm == 0 ){ //BB:
	theta0 = BBang_file[fname]; //on beam left:
	R0 = BBdist_file[fname];
	B_pth1 = 0.1976 + 0.4764*R0;
	if( downbend != 0 ){
	  A_pth1 = 0.2885 * BBscale/0.9644;
	  B_pth1 = 0.994; //to be refined later
	}
	
      } else if( arm == 1 ){ //SBS:
	theta0 = SBSang_file[fname]; //on beam right
	R0 = SBSdist_file[fname];
      }

      //      cout << "fname, arm, theta0 = " << fname << ", " << arm << ", " << theta0 << endl;

      //Electron time-of-flight (to hodoscope):
      //In g4sbs, the hodo is about 3 meters downstream of the magnet:
      //Don't worry about SBS (yet).
      double etof_avg_BB = (3.0 + R0) / 3.0e8*1.e9; //convert to ns. 
      double Lpath_avg_BB = 3.0+R0;
      
      double vx, vy, vz, px, py, pz;
      double p, xptar, yptar, ytar, xtar;
      double xfp, yfp, xpfp, ypfp;

      double thetabend;
      
      vx = T->ev_vx;
      vy = T->ev_vy;
      vz = T->ev_vz;

      TVector3 vertex(vx,vy,vz);

      bool goodtrack = false;
      bool good_t_hodo = false;

      double thodo = 0.0;
    
      if( arm == 0 ){
	if( T->Earm_BBGEM_Track_ntracks == 1 && (*(T->Earm_BBGEM_Track_MID))[0] == 0
	    && (*(T->Earm_BBGEM_Track_P))[0]/T->ev_ep >= 0.95
	    && (*(T->Earm_BBGEM_Track_Chi2fit))[0]/(*(T->Earm_BBGEM_Track_NDF))[0] <= chi2cut ){ //BB

	  bool goodfirsthit=false;
	  for( int ihit=0; ihit<T->Earm_BBGEM_hit_nhits; ihit++ ){
	    if( (*(T->Earm_BBGEM_hit_mid))[ihit] == 0 && (*(T->Earm_BBGEM_hit_p))[ihit] >= 0.994*T->ev_ep ){
	      //&& (*(T->Earm_BBGEM_hit_plane))[ihit] == 1 ){
	      goodfirsthit = true;
	    }
	  }

	  //
	  if( T->Earm_BBHodoScint_det_esum > 0.01 ){

	    double Ehodomax = 0.0;
	    double ihitmax_hodo = -1;
	    for( int ihit=0; ihit<T->Earm_BBHodoScint_hit_nhits; ihit++ ){
	      double Ehit = T->Earm_BBHodoScint_hit_sumedep->at(ihit);
	      double thit = T->Earm_BBHodoScint_hit_tavg->at(ihit);

	      if( Ehit > Ehodomax ){
		Ehodomax = Ehit;
		ihitmax_hodo = ihit;
	      }
	    }
	    if( ihitmax_hodo >= 0 && Ehodomax > 0.01 ){
	      good_t_hodo = true;
	      thodo = T->Earm_BBHodoScint_hit_tavg->at(ihitmax_hodo);
	    }
	  }
	  
	  p = T->ev_ep; //Fit "true" momentum and "true" angles:

	  //px = T->ev_epx;
	  //py = T->ev_epy;
	  //pz = T->ev_epz;
      
	  TVector3 pvect( p*sin(T->ev_th)*cos(T->ev_ph), p*sin(T->ev_th)*sin(T->ev_ph), p*cos(T->ev_th) );
	  TVector3 BB_zaxis( sin(theta0), 0.0, cos(theta0) ); //BB is on beam left, global x axis points to beam left
	  TVector3 BB_xaxis(0,-1,0); //X axis of transport coordinates is vertically down:
	  //TVector3 BB_xaxis = (BB_yaxis.Cross(BB_zaxis)).Unit();
	  TVector3 BB_yaxis = (BB_zaxis.Cross(BB_xaxis)).Unit();

	  TVector3 BBGEM_zaxis = BB_zaxis;
	  TVector3 BBGEM_yaxis = BB_yaxis;
	  BBGEM_zaxis.Rotate( -tracker_pitch_angle, BBGEM_yaxis );
	  TVector3 BBGEM_xaxis = BBGEM_yaxis.Cross(BBGEM_zaxis).Unit();

	  // cout << "BBGEM axes = " << endl;
	  // BBGEM_xaxis.Print();
	  // BBGEM_yaxis.Print();
	  // BBGEM_zaxis.Print();
	  
	  //TVector3 pvect_BB = pvect.Dot(BB_xaxis)*BB_xaxis + pvect.Dot(BB_yaxis)*BB_yaxis + pvect.Dot(BB_zaxis)*BB_zaxis;
	  TVector3 pvect_BB( pvect.Dot(BB_xaxis), pvect.Dot(BB_yaxis), pvect.Dot(BB_zaxis) );

	  xptar = pvect_BB.X()/pvect_BB.Z();
	  yptar = pvect_BB.Y()/pvect_BB.Z();

	  //cout << "true yptar = " << yptar << ", theta0 = " << theta0 << endl;
      
	  TVector3 vertex_BB(vertex.Dot(BB_xaxis), vertex.Dot(BB_yaxis), vertex.Dot(BB_zaxis) );

	  ytar = vertex_BB.Y() - yptar * vertex_BB.Z();
	  xtar = vertex_BB.X() - xptar * vertex_BB.Z();

	  xfp = (*(T->Earm_BBGEM_Track_X))[0];
	  yfp = (*(T->Earm_BBGEM_Track_Y))[0];
	  xpfp = (*(T->Earm_BBGEM_Track_Xp))[0];
	  ypfp = (*(T->Earm_BBGEM_Track_Yp))[0];

	  hxfp->Fill(xfp);
	  hyfp->Fill(yfp);
	  hxpfp->Fill(xpfp);
	  hypfp->Fill(ypfp);
	  
	  hxfprecon->Fill( (*(T->Earm_BBGEM_Track_Xfit))[0] );
	  hyfprecon->Fill( (*(T->Earm_BBGEM_Track_Yfit))[0] );
	  hxpfprecon->Fill( (*(T->Earm_BBGEM_Track_Xpfit))[0] );
	  hypfprecon->Fill( (*(T->Earm_BBGEM_Track_Ypfit))[0] );
	
	  TVector3 pvect_fp_BB( xpfp, ypfp, 1.0 );
	  pvect_fp_BB = pvect_fp_BB.Unit();

	  TVector3 phat_fp_global = pvect_fp_BB.X() * BBGEM_xaxis + pvect_fp_BB.Y() * BBGEM_yaxis + pvect_fp_BB.Z() * BBGEM_zaxis;
	  TVector3 phat_targ_global = pvect.Unit();
	  
	  //TVector3 BendPlane_UnitNormal = phat_targ_global.Cross( phat_fp_global ).Unit();
	  thetabend = acos( phat_fp_global.Dot( phat_targ_global ) );

	  // cout << "BigBite thetabend exact = " << thetabend*57.3 << " degrees" << endl;
	  // cout << "BigBite thetabend approx = " << 57.3*(tracker_pitch_angle + atan(xptar)-atan(xpfp)) << endl;

	  double xface = xtar + xptar * R0;
	  double yface = ytar + yptar * R0;

	  double xsieve = xtar + xptar * dsieve;
	  double ysieve = ytar + yptar * dsieve;
	  
	  //shouldn't hardcode this but whatever:
	  if( fabs( xsieve ) < 0.24 &&
	      fabs( ysieve ) < 0.08 && fabs(xptar)<cutxptar_true && fabs(yptar)<cutyptar_true &&
	      fabs(ytar) < cutytar_true ){
	    goodtrack = true;
	  }
	  goodtrack = goodtrack && goodfirsthit;
	}
      } else if( arm == 1 ){
	if( T->Harm_CEPolFront_Track_ntracks==1 && (*(T->Harm_CEPolFront_Track_MID))[0] == 0 &&
	    (*(T->Harm_CEPolFront_Track_P))[0]/T->ev_ep >= 0.95 &&
	    (*(T->Harm_CEPolFront_Track_Chi2fit))[0]/(*(T->Harm_CEPolFront_Track_NDF))[0] <= chi2cut ){

	  bool goodfirsthit = false;
	  for( int ihit=0; ihit<T->Harm_CEPolFront_hit_nhits; ihit++ ){
	    if( (*(T->Harm_CEPolFront_hit_mid))[ihit]==0&&(*(T->Harm_CEPolFront_hit_p))[ihit] >= 0.997*T->ev_ep ){
	      goodfirsthit=true;
	    }
	  }

	  p = T->ev_ep;
	  
	  TVector3 pvect( p*sin(T->ev_th)*cos(T->ev_ph),
			  p*sin(T->ev_th)*sin(T->ev_ph),
			  p*cos(T->ev_th) );
	  
	  TVector3 SBS_zaxis( -sin(theta0), 0, cos(theta0) );
	  TVector3 SBS_xaxis(0,-1,0);
	  TVector3 SBS_yaxis = ( SBS_zaxis.Cross(SBS_xaxis)).Unit();

	  TVector3 SBSGEM_zaxis = SBS_zaxis;
	  TVector3 SBSGEM_yaxis = SBS_yaxis;
	  SBSGEM_zaxis.Rotate( -tracker_pitch_angle, SBSGEM_yaxis );
	  TVector3 SBSGEM_xaxis = SBSGEM_yaxis.Cross( SBSGEM_zaxis ).Unit();


	  TVector3 pvect_SBS( pvect.Dot(SBS_xaxis), pvect.Dot(SBS_yaxis), pvect.Dot(SBS_zaxis) );

	  xptar = pvect_SBS.X()/pvect_SBS.Z();
	  yptar = pvect_SBS.Y()/pvect_SBS.Z();

	  TVector3 vertex_SBS( vertex.Dot( SBS_xaxis ),
			       vertex.Dot( SBS_yaxis ),
			       vertex.Dot( SBS_zaxis ) );

	  ytar = vertex_SBS.Y() - yptar * vertex_SBS.Z();
	  xtar = vertex_SBS.X() - xptar * vertex_SBS.Z();

	  xfp = (*(T->Harm_CEPolFront_Track_X))[0];
	  yfp = (*(T->Harm_CEPolFront_Track_Y))[0];
	  xpfp = (*(T->Harm_CEPolFront_Track_Xp))[0];
	  ypfp = (*(T->Harm_CEPolFront_Track_Yp))[0];

	  hxfp->Fill(xfp);
	  hyfp->Fill(yfp);
	  hxpfp->Fill(xpfp);
	  hypfp->Fill(ypfp);
	  
	  hxfprecon->Fill( (*(T->Harm_CEPolFront_Track_Xfit))[0] );
	  hyfprecon->Fill( (*(T->Harm_CEPolFront_Track_Yfit))[0] );
	  hxpfprecon->Fill( (*(T->Harm_CEPolFront_Track_Xpfit))[0] );
	  hypfprecon->Fill( (*(T->Harm_CEPolFront_Track_Ypfit))[0] );

	  TVector3 pvect_fp_SBS( xpfp, ypfp, 1.0 );
	  pvect_fp_SBS = pvect_fp_SBS.Unit();

	  TVector3 phat_fp_global = pvect_fp_SBS.X() * SBSGEM_xaxis +
	    pvect_fp_SBS.Y() * SBSGEM_yaxis +
	    pvect_fp_SBS.Z() * SBSGEM_zaxis;

	  TVector3 phat_targ_global = pvect.Unit();

	  thetabend = acos( phat_fp_global.Dot( phat_targ_global ) );
	  
	  goodtrack = goodfirsthit && fabs(xptar)<cutxptar_true && fabs(yptar)<cutyptar_true &&
									       fabs(ytar)<cutytar_true;
									       
	  // }
	}
      }
      
      if( goodtrack ){ //Increment fit matrices:
      
	hxtar->Fill(xtar);
	hytar->Fill(ytar);
	hxptar->Fill(xptar);
	hyptar->Fill(yptar);
	hp->Fill(p);
	hvz->Fill(vz);
      
	//Are the equations below correct? 
	// chi^2 = sum_{i=1}^{N} (xtgt_i - sum_jklmn C_jklmn xfp^j yfp^k xpfp^l ypfp^m xtar^n)^2
	// Momentum reconstruction is a special case, we want to fit a 
	// functional form of 1/p = A(x,y,x',y') * thetabend 
	// in this case, chi^2 = sum_i (1/p - (A(x,y,x',y')*thetabend))^2
     
	//double thetabend = tracker_pitch_angle + atan(xptar) - atan(xpfp);
	
	vector<double> term(nparams);
	vector<double> fterm(nparams);
	//vector<double> tterm(nparams);
	//      vector<double> term_y(nparams);
	int ipar=0;
	for(int i=0; i<=order; i++){
	  for(int j=0; j<=order-i; j++){
	    for(int k=0; k<=order-i-j; k++){
	      for(int l=0; l<=order-i-j-k; l++){
		for(int m=0; m<=order-i-j-k-l; m++){
		  term[ipar] = pow(xfp,m)*pow(yfp,l)*pow(xpfp,k)*pow(ypfp,j)*pow(xtar,i);
		
		  b_xptar(ipar) += term[ipar]*xptar;
		
		  //if( i == 0 ){
		  b_yptar(ipar) += term[ipar]*yptar;
		  b_ytar(ipar) += term[ipar]*ytar;
		  //}
		  //b_ytar(ipar) += term[ipar]*T->ev_vz;
		  //} 
		  //b_pinv(ipar) += term[ipar]*(p/pcentral[arm]-1.0);

		  if( pexpansion_flag == 0 ){
		    b_pinv(ipar) += term[ipar]*p*thetabend;
		  } else if( arm == 0 && pexpansion_flag == 2 ){
		    double pth1 = A_pth1*( 1.0 + B_pth1 * xptar ); //GeV/c * rad
		    double delta = p*(thetabend)/pth1 - 1.0;
		    b_pinv(ipar) += term[ipar]*delta;
		  } else {
		    b_pinv(ipar) += term[ipar]/p;
		  }

		  //initially, let's try expanding etof in terms of focal plane variables:
		  //tterm[ipar] = term[ipar];

		  if( good_t_hodo && arm == 0 ){
		    b_etof[ipar] += term[ipar]*(thodo-etof_avg_BB);
		  } else {
		    
		  }
		  
		  fterm[ipar] = pow(xptar,m)*pow(yptar,l)*pow(ytar,k)*pow(1.0/p,j)*pow(xtar,i);
		  b_xfp[ipar] += fterm[ipar]*xfp;
		  b_yfp[ipar] += fterm[ipar]*yfp;
		  b_xpfp[ipar] += fterm[ipar]*xpfp;
		  b_ypfp[ipar] += fterm[ipar]*ypfp;
		  
		  ipar++;
		}
	      }
	    }
	  }
	}
      
	for(ipar=0; ipar<nparams; ipar++){
	  for(int jpar=0; jpar<nparams; jpar++){
	    //double term_ij = term[ipar]*term[jpar];
	    M(ipar,jpar) += term[ipar]*term[jpar];
	    Mforward(ipar,jpar) += fterm[ipar]*fterm[jpar];
	    //
	  }
	}
      
      }
    }
  }
    

  cout << "order, nparams = " << nparams << endl;

  if( xtar_flag != order ){
    for(int ipar=0; ipar<nparams; ipar++){
      double sumexpon = xtar_expon[ipar] + xfp_expon[ipar] + yfp_expon[ipar]
	+ xpfp_expon[ipar] + ypfp_expon[ipar];
      if( xtar_expon[ipar] > 0 ){
	M(ipar,ipar) = 1.0;
	b_xptar(ipar) = 0.0;
	b_yptar(ipar) = 0.0;
	b_ytar(ipar) = 0.0;
	b_pinv(ipar) = 0.0;
	b_etof(ipar) = 0.0;

	Mforward(ipar,ipar) = 1.0;
	b_xfp(ipar) = 0.0;
	b_yfp(ipar) = 0.0;
	b_xpfp(ipar) = 0.0;
	b_ypfp(ipar) = 0.0;
	for(int jpar=0; jpar<nparams; jpar++){
	  if( jpar != ipar ){
	    M(ipar,jpar) = 0.0;
	    Mforward(ipar,jpar) = 0.0;
	  }
	}
      }
    }
  }
  
  cout << "Setting up SVD for xptar:" << endl;
  TDecompSVD A_xptar(M);
  cout << "Setting up SVD for pinv:" << endl;
  TDecompSVD A_pinv(M);
  cout << "Setting up SVD for yptar:" << endl;
  TDecompSVD A_yptar(M);
  cout << "Setting up SVD for ytar:" << endl;
  TDecompSVD A_ytar(M);

  cout << "Setting up SVD for etof:" << endl;
  TDecompSVD A_etof(M);
  

  cout << "Solving for xptar coefficients:" << endl;
  bool good_xptar = A_xptar.Solve(b_xptar);
  cout << "xptar done, success = " << good_xptar << endl;
  cout << "Solving for yptar coefficients:" << endl;
  bool good_yptar = A_yptar.Solve(b_yptar);
  cout << "yptar done, success = " << good_yptar << endl;
  cout << "Solving for ytar coefficients:" << endl;
  bool good_ytar = A_ytar.Solve(b_ytar);
  cout << "ytar done, success = " << good_ytar << endl;
  cout << "Solving for 1/p coefficients:" << endl;
  bool good_pinv = A_pinv.Solve(b_pinv);
  cout << "1/p done, success = " << good_pinv << endl;

  if( arm == 0 ){
    bool good_etof = A_etof.Solve(b_etof);
    cout << "e^- TOF done, success = " << good_etof << endl;
  }
  // Solving for forward matrix elements
  TDecompSVD A_xfp(Mforward);
  TDecompSVD A_yfp(Mforward);
  TDecompSVD A_xpfp(Mforward);
  TDecompSVD A_ypfp(Mforward);

  A_xfp.Solve(b_xfp);
  A_yfp.Solve(b_yfp);
  A_xpfp.Solve(b_xpfp);
  A_ypfp.Solve(b_ypfp);

  if( arm == 0 && pexpansion_flag == 2 ){
    opticsfile << "bb.preconflag = 1" << endl << endl;
    TString stemp;
    opticsfile << stemp.Format("bb.A_pth1 = %15.9g",A_pth1) << endl;
    opticsfile << stemp.Format("bb.B_pth1 = %15.9g",B_pth1) << endl;
    opticsfile << "bb.C_pth1 = 0.0" << endl << endl;
  }

  if( arm == 0 ){
    TString stemp;
    double etof_avg = (BBdist_default + 3.0)/3.0e8*1.e9;
    //Write out the TOF expansion:
    TOFfile << "# Average TOF to hodoscope based on BigBite magnet distance" << endl;
    TOFfile << stemp.Format("bb.etof_avg = %15.9g", etof_avg) << endl << endl;
    TOFfile << "bb.etof_parameters = " << endl;
  }
  // Write out the backward optics file:
  // Order is: xptar, yptar, ytar, 1/p, exponents:
  opticsfile << nparams << endl;
  int ipar = 0;
  for(int i=0; i<=order; i++){
    for(int j=0; j<=order-i; j++){
      for(int k=0; k<=order-i-j; k++){
	for(int l=0; l<=order-i-j-k; l++){
	  for(int m=0; m<=order-i-j-k-l; m++){
	    char ccoeff[100];
	    sprintf( ccoeff, "%15.8g", b_xptar(ipar) );
	    opticsfile << ccoeff;
	    sprintf( ccoeff, "%15.8g", b_yptar(ipar) );
	    opticsfile << ccoeff;
	    sprintf( ccoeff, "%15.8g", b_ytar(ipar) );
	    opticsfile << ccoeff;
	    sprintf( ccoeff, "%15.8g", b_pinv(ipar) );
	    opticsfile << ccoeff;

	    sprintf( ccoeff, "%15.8g", b_etof(ipar) );
	    TOFfile << ccoeff;
	    
	    char cexpon[100];
	    sprintf( cexpon, "  %d %d %d %d %d", m, l, k, j, i );
	    opticsfile << cexpon << endl;

	    TOFfile << cexpon << endl;
	    
	    ipar++;
	  }
	}
      }
    }
  }
    
  // Write out the forward optics file:
  // Order is: xfp, yfp, xpfp, ypfp, exponents:
  fopticsfile << nparams << endl;
  ipar = 0;
  for(int i=0; i<=order; i++){
    for(int j=0; j<=order-i; j++){
      for(int k=0; k<=order-i-j; k++){
	for(int l=0; l<=order-i-j-k; l++){
	  for(int m=0; m<=order-i-j-k-l; m++){
	    char ccoeff[100];
	    sprintf( ccoeff, "%15.8g", b_xfp(ipar) );
	    fopticsfile << ccoeff;
	    sprintf( ccoeff, "%15.8g", b_yfp(ipar) );
	    fopticsfile << ccoeff;
	    sprintf( ccoeff, "%15.8g", b_xpfp(ipar) );
	    fopticsfile << ccoeff;
	    sprintf( ccoeff, "%15.8g", b_ypfp(ipar) );
	    fopticsfile << ccoeff;

	    char cexpon[100];
	    sprintf( cexpon, "  %d %d %d %d %d", m, l, k, j, i );
	    fopticsfile << cexpon << endl;
	    
	    ipar++;
	  }
	}
      }
    }
  }  

  double vx, vy, vz, px, py, pz;
  double p, xptar, yptar, ytar, xtar;
  double p_fit, xptar_fit, yptar_fit, ytar_fit; //Fit is reconstructed using fit coefficients, no smearing for detector resolution
  double p_recon, xptar_recon, yptar_recon, ytar_recon; //recon is reconstructed using fit coefficients, fp quantities smeared by det. resolution
  double pthetabend_true;
  double pthetabend_fit, pthetabend_recon;
  double pinv_fit, pinv_recon;
  double xfp, yfp, xpfp, ypfp;
  double xfp_fit, yfp_fit, xpfp_fit, ypfp_fit;
  double xfpforward, yfpforward, xpfpforward, ypfpforward;
  double vz_fit, vz_recon;
  double thetabend_true;
  double thetabend_fit;
  double thetabend_recon;
  double xtar_recon, xtar_fit;
  double etof, etof_recon; 
  
  TTree *tout = new TTree("tout","BigBite or SBS optics fit results, diagnostic ROOT tree");

  tout->Branch("etof",&etof);
  tout->Branch("etof_recon",&etof_recon);
  tout->Branch("vxtrue",&vx);
  tout->Branch("vytrue",&vy);
  tout->Branch("vztrue",&vz);
  tout->Branch("pxtrue",&px);
  tout->Branch("pytrue",&py);
  tout->Branch("pztrue",&pz);
  tout->Branch("ptrue",&p);
  tout->Branch("xptartrue",&xptar);
  tout->Branch("yptartrue",&yptar);
  tout->Branch("ytartrue",&ytar);
  tout->Branch("xtartrue",&xtar);
  tout->Branch("thetabend_true",&thetabend_true);
  tout->Branch("p_fit",&p_fit);
  tout->Branch("xptar_fit",&xptar_fit);
  tout->Branch("yptar_fit",&yptar_fit);
  tout->Branch("ytar_fit",&ytar_fit);
  tout->Branch("xtar_fit",&xtar_fit);
  tout->Branch("thetabend_fit",&thetabend_fit);
  tout->Branch("p_recon",&p_recon);
  tout->Branch("xptar_recon",&xptar_recon);
  tout->Branch("yptar_recon",&yptar_recon);
  tout->Branch("ytar_recon",&ytar_recon);
  tout->Branch("xtar_recon",&xtar_recon);
  tout->Branch("vz_fit",&vz_fit);
  tout->Branch("vz_recon",&vz_recon);
  tout->Branch("thetabend_recon",&thetabend_recon);
  tout->Branch("xfptrue",&xfp);
  tout->Branch("yfptrue",&yfp);
  tout->Branch("xpfptrue",&xpfp);
  tout->Branch("ypfptrue",&ypfp);
  tout->Branch("xfp_fit",&xfp_fit);
  tout->Branch("xpfp_fit",&xpfp_fit);
  tout->Branch("yfp_fit",&yfp_fit);
  tout->Branch("ypfp_fit",&ypfp_fit);
  tout->Branch("xfpforward",&xfpforward);
  tout->Branch("xpfpforward",&xpfpforward);
  tout->Branch("yfpforward",&yfpforward);
  tout->Branch("ypfpforward",&ypfpforward);
  
  cout << "arm = " << arm << endl;
  
  nevent = 0;
  while( T->GetEntry( elist->GetEntry(nevent++))){
    if( nevent%1000 == 0 ){
      cout << nevent << endl;
    }
    
    TString fname = C->GetFile()->GetName();
    
    if( bad_file_list.find(fname) == bad_file_list.end() ){
      
      double theta0;
      double R0;
      if( arm == 0 ){ //BB:
	theta0 = BBang_file[fname]; //on beam left:
	R0 = BBdist_file[fname];
      } else if( arm == 1 ){ //SBS:
	theta0 = SBSang_file[fname]; //on beam right
	R0 = SBSdist_file[fname];
      }

      double etof_avg_BB = (3.0 + R0) / 3.0e8*1.e9; //convert to ns. 
      
      pinv_fit = 0.0;
      pinv_recon = 0.0;

      vx = T->ev_vx;
      vy = T->ev_vy;
      vz = T->ev_vz;

      TVector3 vertex(vx,vy,vz);

      bool goodtrack = false;

      bool good_t_hodo = false;

      double thodo = 0.0;

      TVector3 spec_xaxis_fp,spec_yaxis_fp, spec_zaxis_fp;
      TVector3 spec_xaxis_tgt,spec_yaxis_tgt, spec_zaxis_tgt;
      
      if( arm == 0 ){
	if( T->Earm_BBGEM_Track_ntracks == 1 && (*(T->Earm_BBGEM_Track_MID))[0] == 0
	    && (*(T->Earm_BBGEM_Track_P))[0]/T->ev_ep >= 0.99 
	    && (*(T->Earm_BBGEM_Track_Chi2fit))[0]/(*(T->Earm_BBGEM_Track_NDF))[0] <= chi2cut ){ //BB

	  double Ehodomax = 0.0;
	  double ihitmax_hodo = -1;
	  for( int ihit=0; ihit<T->Earm_BBHodoScint_hit_nhits; ihit++ ){
	    double Ehit = T->Earm_BBHodoScint_hit_sumedep->at(ihit);
	    double thit = T->Earm_BBHodoScint_hit_tavg->at(ihit);
	    
	    if( Ehit > Ehodomax ){
	      Ehodomax = Ehit;
	      ihitmax_hodo = ihit;
	    }
	  }
	  if( ihitmax_hodo >= 0 && Ehodomax > 0.0 ){
	    good_t_hodo = true;
	    thodo = T->Earm_BBHodoScint_hit_tavg->at(ihitmax_hodo);
	  }
	  
	  p = T->ev_ep;
	  px = T->ev_epx;
	  py = T->ev_epy;
	  pz = T->ev_epz;
      
	  TVector3 pvect( p*sin(T->ev_th)*cos(T->ev_ph), p*sin(T->ev_th)*sin(T->ev_ph), p*cos(T->ev_th) );
	  TVector3 BB_zaxis( sin(theta0), 0.0, cos(theta0) ); //BB is on beam right, global x axis points to beam left
	  TVector3 BB_xaxis(0,-1,0); //X axis of transport coordinates is vertically down:
	  //TVector3 BB_xaxis = (BB_yaxis.Cross(BB_zaxis)).Unit();
	  TVector3 BB_yaxis = (BB_zaxis.Cross(BB_xaxis)).Unit();

	  spec_xaxis_tgt = BB_xaxis;
	  spec_yaxis_tgt = BB_yaxis;
	  spec_zaxis_tgt = BB_zaxis;
	  
	  spec_zaxis_fp = BB_zaxis;
	  spec_yaxis_fp = BB_yaxis;
	  spec_zaxis_fp.Rotate(-tracker_pitch_angle, spec_yaxis_fp);
	  spec_xaxis_fp = spec_yaxis_fp.Cross(spec_zaxis_fp).Unit();
	  
	  //TVector3 pvect_BB = pvect.Dot(BB_xaxis)*BB_xaxis + pvect.Dot(BB_yaxis)*BB_yaxis + pvect.Dot(BB_zaxis)*BB_zaxis;
	  TVector3 pvect_BB( pvect.Dot(BB_xaxis), pvect.Dot(BB_yaxis), pvect.Dot(BB_zaxis) );

	  xptar = pvect_BB.X()/pvect_BB.Z();
	  yptar = pvect_BB.Y()/pvect_BB.Z();
      
	  TVector3 vertex_BB(vertex.Dot(BB_xaxis), vertex.Dot(BB_yaxis), vertex.Dot(BB_zaxis) );

	  ytar = vertex_BB.Y() - yptar * vertex_BB.Z();
	  xtar = vertex_BB.X() - xptar * vertex_BB.Z();

	  xfp = (*(T->Earm_BBGEM_Track_X))[0];
	  yfp = (*(T->Earm_BBGEM_Track_Y))[0];
	  xpfp = (*(T->Earm_BBGEM_Track_Xp))[0];
	  ypfp = (*(T->Earm_BBGEM_Track_Yp))[0];

	  xfp_fit = (*(T->Earm_BBGEM_Track_Xfit))[0];
	  yfp_fit = (*(T->Earm_BBGEM_Track_Yfit))[0];
	  xpfp_fit = (*(T->Earm_BBGEM_Track_Xpfit))[0];
	  ypfp_fit = (*(T->Earm_BBGEM_Track_Ypfit))[0];

	  hxfpdiff->Fill(xfp_fit-xfp);
	  hyfpdiff->Fill(yfp_fit-yfp);
	  hxpfpdiff->Fill(xpfp_fit-xpfp);
	  hypfpdiff->Fill(ypfp_fit-ypfp);
	  
	  //	  TVector3 phat_targ = pvect_BB.Unit();
	  TVector3 phat_fp(xpfp,ypfp,1.0);
	  phat_fp = phat_fp.Unit();

	  TVector3 phat_fp_global = phat_fp.X() * spec_xaxis_fp + phat_fp.Y() * spec_yaxis_fp + phat_fp.Z() * spec_zaxis_fp;
	  
	  TVector3 phat_targ_global = pvect.Unit();

	  thetabend_true = acos( phat_fp_global.Dot( phat_targ_global ) );

	  
	  goodtrack = true;
	  
	}
      } else if( arm == 1 ){
	if( T->Harm_CEPolFront_Track_ntracks==1 && (*(T->Harm_CEPolFront_Track_MID))[0] == 0 &&
	    (*(T->Harm_CEPolFront_Track_P))[0]/T->ev_ep >= 0.95 &&
	    (*(T->Harm_CEPolFront_Track_Chi2fit))[0]/(*(T->Harm_CEPolFront_Track_NDF))[0] <= chi2cut ){

	  // bool goodfirsthit = false;
	  // for( int ihit=0; ihit<T->Harm_CEPolFront_hit_nhits; ihit++ ){
	  //   if( (*(T->Harm_CEPolFront_hit_mid))[ihit]==0&&(*(T->Harm_CEPolFront_hit_p))[ihit] >= 0.994*T->ev_ep ){
	  //     goodfirsthit=true;
	  //   }
	  // }

	  good_t_hodo = true;
	  //thodo = T->Harm_CEPolFront_Track_
	  
	  p = T->ev_ep;
	  px = T->ev_epx;
	  py = T->ev_epy;
	  pz = T->ev_epz;
	  
	  TVector3 pvect( p*sin(T->ev_th)*cos(T->ev_ph),
			  p*sin(T->ev_th)*sin(T->ev_ph),
			  p*cos(T->ev_th) );
	  
	  TVector3 SBS_zaxis( -sin(theta0), 0, cos(theta0) );
	  TVector3 SBS_xaxis(0,-1,0);
	  TVector3 SBS_yaxis = ( SBS_zaxis.Cross(SBS_xaxis)).Unit();

	  spec_xaxis_tgt = SBS_xaxis;
	  spec_yaxis_tgt = SBS_yaxis;
	  spec_zaxis_tgt = SBS_zaxis;
	  
	  spec_zaxis_fp = SBS_zaxis;
	  spec_yaxis_fp = SBS_yaxis;
	  spec_zaxis_fp.Rotate(-tracker_pitch_angle, spec_yaxis_fp);
	  spec_xaxis_fp = spec_yaxis_fp.Cross(spec_zaxis_fp).Unit();
	  
	  TVector3 SBSGEM_zaxis = SBS_zaxis;
	  TVector3 SBSGEM_yaxis = SBS_yaxis;
	  SBSGEM_zaxis.Rotate( -tracker_pitch_angle, SBSGEM_yaxis );
	  TVector3 SBSGEM_xaxis = SBSGEM_yaxis.Cross( SBSGEM_zaxis ).Unit();


	  TVector3 pvect_SBS( pvect.Dot(SBS_xaxis), pvect.Dot(SBS_yaxis), pvect.Dot(SBS_zaxis) );

	  xptar = pvect_SBS.X()/pvect_SBS.Z();
	  yptar = pvect_SBS.Y()/pvect_SBS.Z();

	  TVector3 vertex_SBS( vertex.Dot( SBS_xaxis ),
			       vertex.Dot( SBS_yaxis ),
			       vertex.Dot( SBS_zaxis ) );

	  ytar = vertex_SBS.Y() - yptar * vertex_SBS.Z();
	  xtar = vertex_SBS.X() - xptar * vertex_SBS.Z();

	  xfp = (*(T->Harm_CEPolFront_Track_X))[0];
	  yfp = (*(T->Harm_CEPolFront_Track_Y))[0];
	  xpfp = (*(T->Harm_CEPolFront_Track_Xp))[0];
	  ypfp = (*(T->Harm_CEPolFront_Track_Yp))[0];

	  xfp_fit = (*(T->Harm_CEPolFront_Track_Xfit))[0];
	  yfp_fit = (*(T->Harm_CEPolFront_Track_Yfit))[0];
	  xpfp_fit = (*(T->Harm_CEPolFront_Track_Xpfit))[0];
	  ypfp_fit = (*(T->Harm_CEPolFront_Track_Ypfit))[0];
	  
	  hxfpdiff->Fill(xfp_fit-xfp);
	  hyfpdiff->Fill(yfp_fit-yfp);
	  hxpfpdiff->Fill(xpfp_fit-xpfp);
	  hypfpdiff->Fill(ypfp_fit-ypfp);
	  
	  // hxfprecon->Fill( (*(T->Harm_CEPolFront_Track_Xfit))[0] );
	  // hyfprecon->Fill( (*(T->Harm_CEPolFront_Track_Yfit))[0] );
	  // hxpfprecon->Fill( (*(T->Harm_CEPolFront_Track_Xpfit))[0] );
	  // hypfprecon->Fill( (*(T->Harm_CEPolFront_Track_Ypfit))[0] );

	  TVector3 pvect_fp_SBS( xpfp, ypfp, 1.0 );
	  pvect_fp_SBS = pvect_fp_SBS.Unit();

	  TVector3 phat_fp_global = pvect_fp_SBS.X() * SBSGEM_xaxis +
	    pvect_fp_SBS.Y() * SBSGEM_yaxis +
	    pvect_fp_SBS.Z() * SBSGEM_zaxis;

	  TVector3 phat_targ_global = pvect.Unit();

	  thetabend_true = acos( phat_fp_global.Dot( phat_targ_global ) );
	  
	  //goodtrack = goodfirsthit;
	  goodtrack = true;
	  // }
	}
	
      }   
    
    
      if( goodtrack ){ //do reconstruction:
	int ipar = 0;

	int niter=1;
	
	xtar_recon = xtar;
	xtar_fit = xtar;

	if( xtar_flag != 0 ){ 
	  niter = 10;

	  xtar_recon = -vy;
	  xtar_fit   = -vy;
	} 

	etof = thodo;

	for( int iter=0; iter<niter; iter++){
	
	  ipar=0;

	  xptar_fit = 0.0;
	  yptar_fit = 0.0;
	  ytar_fit = 0.0;
	  pthetabend_fit = 0.0;
	  pinv_fit = 0.0;
	  
	  xptar_recon = 0.0;
	  yptar_recon = 0.0;
	  ytar_recon = 0.0;
	  pthetabend_recon = 0.0;
	  pinv_recon = 0.0;

	  xfpforward = 0.0;
	  yfpforward = 0.0;
	  xpfpforward = 0.0;
	  ypfpforward = 0.0;	  

	  etof_recon = etof_avg_BB; //start with average. The expansion is around this average
	  
	  ipar=0;
	  for(int i=0; i<=order; i++){
	    for(int j=0; j<=order-i; j++){
	      for(int k=0; k<=order-i-j; k++){
		for(int l=0; l<=order-i-j-k; l++){
		  for(int m=0; m<=order-i-j-k-l; m++){
		    double term = pow(xfp,m)*pow(yfp,l)*pow(xpfp,k)*pow(ypfp,j)*pow(xtar_fit,i);
		    xptar_fit += b_xptar(ipar)*term;
		    yptar_fit += b_yptar(ipar)*term;
		    ytar_fit += b_ytar(ipar)*term;
		    pthetabend_fit += b_pinv(ipar)*term;
		    pinv_fit += b_pinv(ipar)*term;
		    
		    term = pow(xfp_fit,m)*pow(yfp_fit,l)*pow(xpfp_fit,k)*pow(ypfp_fit,j)*pow(xtar_recon,i);
		    xptar_recon += b_xptar(ipar)*term;
		    yptar_recon += b_yptar(ipar)*term;
		    ytar_recon += b_ytar(ipar)*term;
		    pinv_recon += b_pinv(ipar)*term;
		    pthetabend_recon += b_pinv(ipar)*term;

		    double fterm = pow(xptar,m)*pow(yptar,l)*pow(ytar,k)*pow(1.0/p,j)*pow(xtar,i);
		    xfpforward += fterm * b_xfp(ipar);
		    yfpforward += fterm * b_yfp(ipar);
		    xpfpforward += fterm * b_xpfp(ipar);
		    ypfpforward += fterm * b_ypfp(ipar);

		    if( arm == 0 ) etof_recon += b_etof(ipar) * term;
		    
		    ipar++;
		  }
		}
	      }
	    }
	  }
            
	  //double thetabend_true = tracker_pitch_angle + atan(xptar) - atan(xpfp);

	  // for beam left, we have:
	  // ytar = -vz * sin theta - vz * costheta * yptar
	  // --> vz = -ytar/(sintheta + costheta * yptar)
	  // for beam right, we have:
	  // ytar = vz * sin theta - vz * costheta * yptar
	  
	  if( arm == 0 ){ //BB, beam left:
	    vz_fit = -ytar_fit / (sin(theta0) + cos(theta0)*yptar_fit);
	    vz_recon = -ytar_recon / (sin(theta0) + cos(theta0)*yptar_recon);
	  } else { //SBS, beam right:
	    vz_fit = ytar_fit / (sin(theta0) - cos(theta0)*yptar_fit);
	    vz_recon = ytar_recon / (sin(theta0) - cos(theta0)*yptar_recon);
	  }
	  
	  xtar_fit = -vy - vz_fit * cos(theta0) * xptar_fit;
	  xtar_recon = -vy - vz_recon * cos(theta0) * xptar_recon;
	  
	  //calculate "fit" and "recon" values of thetabend:
	  TVector3 phat_tgt_fit(xptar_fit, yptar_fit, 1.0 );
	  phat_tgt_fit = phat_tgt_fit.Unit();
	  
	  TVector3 phat_tgt_fit_global = phat_tgt_fit.X() * spec_xaxis_tgt +
	    phat_tgt_fit.Y() * spec_yaxis_tgt +
	    phat_tgt_fit.Z() * spec_zaxis_tgt;
	  
	  TVector3 phat_fp_fit(xpfp_fit, ypfp_fit, 1.0 );
	  phat_fp_fit = phat_fp_fit.Unit();
	  
	  TVector3 phat_fp_fit_global = phat_fp_fit.X() * spec_xaxis_fp +
	    phat_fp_fit.Y() * spec_yaxis_fp +
	    phat_fp_fit.Z() * spec_zaxis_fp;
	  
	  thetabend_fit = acos( phat_fp_fit_global.Dot( phat_tgt_fit_global ) );
	  
	  //cout << "Thetabend fit = " << 57.3 * thetabend_fit << endl;
	  
	  //calculate "fit" and "recon" values of thetabend:
	  TVector3 phat_tgt_recon(xptar_recon, yptar_recon, 1.0 );
	  phat_tgt_recon = phat_tgt_recon.Unit();
	  
	  TVector3 phat_tgt_recon_global = phat_tgt_recon.X() * spec_xaxis_tgt +
	    phat_tgt_recon.Y() * spec_yaxis_tgt +
	    phat_tgt_recon.Z() * spec_zaxis_tgt;
	  
	  TVector3 phat_fp_recon(xpfp_fit, ypfp_fit, 1.0 );
	  phat_fp_recon = phat_fp_recon.Unit();
	  
	  TVector3 phat_fp_recon_global = phat_fp_recon.X() * spec_xaxis_fp +
	    phat_fp_recon.Y() * spec_yaxis_fp +
	    phat_fp_recon.Z() * spec_zaxis_fp;
	  
	  thetabend_recon = acos( phat_fp_recon_global.Dot( phat_tgt_recon_global ) );
	  
	  //cout << "Thetabend recon = " << 57.3 * thetabend_recon << endl;
	  
	  if( pexpansion_flag == 0 ){
	    p_fit = pthetabend_fit/thetabend_fit;
	    p_recon = pthetabend_recon/thetabend_recon;
	  } else if ( arm == 0 && pexpansion_flag == 2 ) {
	    double delta_fit = pthetabend_fit;
	    double delta_recon = pthetabend_recon;
	    double pth1_fit = A_pth1*(1.0 + B_pth1 * xptar_fit );
	    double pth1_recon = A_pth1*(1.0 + B_pth1 * xptar_recon );

	    p_fit = pth1_fit * (1.0 + delta_fit)/thetabend_fit;
	    p_recon = pth1_recon * (1.0 + delta_recon)/thetabend_recon;
	  } else {
	    p_fit = 1.0/pthetabend_fit;
	    p_recon = 1.0/pthetabend_recon;
	  }
	  pinv_fit = 1.0/p_fit;
	  pinv_recon = 1.0/p_recon;

	  //p_fit = 1.0/pinv_fit;
	  //p_recon = 1.0/pinv_recon;
	  
	  // double pfirstorder_true, pfirstorder_recon;
	  // if(arm==0){
	  //   pfirstorder_true = (0.286378 + 0.141452*(xfp-0.8*xpfp))/thetabend_true;
	  //   pfirstorder_recon  = (0.286378 + 0.141452*(xfp_fit-0.8*xpfp_fit))/thetabend_recon;
	  // } else {
	  //   pfirstorder_true = -0.5161/thetabend_true;
	  //   pfirstorder_recon = -0.5161/thetabend_recon;
	  // }
	  //p_fit = pthetabend_fit/thetabend_fit;
	  //p_recon = pthetabend_recon/thetabend_recon;
	  
	  //p_fit = pthetabend_fit + pfirstorder_true;
	  //p_recon = pthetabend_recon + pfirstorder_recon;
	  //p_fit = 1.0/pinv_fit;
	  //p_recon = 1.0/pinv_recon;
	  //p_fit = pcentral[arm]*(1.0 + pinv_fit);
	  //p_recon = pcentral[arm]*(1.0 + pinv_recon);

	  // double vz_fit = ytar_fit / ( sin( theta;
	  // double vz_recon = ytar_recon;
	  //double vz_fit, vz_recon;
	  
	}

	hvzdiff->Fill( vz_recon - T->ev_vz );
	hvzdiff_xtar->Fill(xtar, vz_recon - T->ev_vz );
	hvzdiff_ytar->Fill(ytar, vz_recon - T->ev_vz );
	hvzdiff_xptar->Fill(xptar, vz_recon - T->ev_vz );
	hvzdiff_yptar->Fill(yptar, vz_recon - T->ev_vz );
	hvzdiff_p->Fill(p, vz_recon - T->ev_vz );

	// if( arm == 0 ){ //bigbite, beam right:
	// 	ytar_fit = vz_fit * (sin(theta0) - cos(theta0)*yptar_fit);
	// 	ytar_recon = vz_recon * (sin(theta0) - cos(theta0)*yptar_recon);
	// } else { //SBS, beam left:
	// 	ytar_fit = -vz_fit * (sin(theta0) + cos(theta0)*yptar_fit);
	// 	ytar_recon = -vz_recon * (sin(theta0) + cos(theta0)*yptar_recon);
	// }

	hpinvdiff->Fill( pinv_recon - 1.0/p );

	hytarrecon->Fill( ytar_recon );
	hyptarrecon->Fill( yptar_recon );
	hxptarrecon->Fill( xptar_recon );
	hprecon->Fill(p_recon );

	hytardiff->Fill( ytar_recon - ytar );
	hxptardiff->Fill( xptar_recon - xptar );
	hyptardiff->Fill( yptar_recon - yptar );
	hpdiff->Fill( (p_recon/p - 1.0) );
      
	hxptardiff_xtar->Fill( xtar, xptar_recon - xptar );
	hxptardiff_xptar->Fill( xptar, xptar_recon - xptar );
	hxptardiff_yptar->Fill( yptar, xptar_recon - xptar );
	hxptardiff_ytar->Fill( ytar, xptar_recon - xptar );
	hxptardiff_p->Fill( p, xptar_recon - xptar );

	hyptardiff_xtar->Fill(xtar, yptar_recon - yptar );
	hyptardiff_xptar->Fill( xptar, yptar_recon - yptar );
	hyptardiff_yptar->Fill( yptar, yptar_recon - yptar );
	hyptardiff_ytar->Fill( ytar, yptar_recon - yptar );
	hyptardiff_p->Fill( p, yptar_recon - yptar );
      
	hytardiff_xtar->Fill( xtar, ytar_recon - ytar );
	hytardiff_xptar->Fill( xptar, ytar_recon - ytar );
	hytardiff_yptar->Fill( yptar, ytar_recon - ytar );
	hytardiff_ytar->Fill( ytar, ytar_recon - ytar );
	hytardiff_p->Fill( p, ytar_recon - ytar );
      
	hpdiff_xtar->Fill( xtar, p_recon/p - 1.0 );
	hpdiff_xptar->Fill( xptar, p_recon/p - 1.0 );
	hpdiff_yptar->Fill( yptar, p_recon/p - 1.0 );
	hpdiff_ytar->Fill( ytar, p_recon/p - 1.0 );
	hpdiff_p->Fill( p, p_recon/p - 1.0 );

	hpdiff_xfp->Fill( xfp, p_recon/p - 1.0 );
	hpdiff_yfp->Fill( yfp, p_recon/p - 1.0 );
	hpdiff_xpfp->Fill( xpfp, p_recon/p - 1.0 );
	hpdiff_ypfp->Fill( ypfp, p_recon/p - 1.0 );

	hxtardiff->Fill( xtar_recon - xtar );
	hxtarrecon->Fill( xtar_recon );
	hvzrecon->Fill( vz_recon );

	hxysieve->Fill( ytar_recon + dsieve * yptar_recon, xtar_recon + dsieve * xptar_recon );

	// Diagnostic histos for forward optics fit
	hxfpdiff_foptics->Fill(xfpforward-xfp);
	hyfpdiff_foptics->Fill(yfpforward-yfp);
	hxpfpdiff_foptics->Fill(xpfpforward-xpfp);
	hypfpdiff_foptics->Fill(ypfpforward-ypfp);

	hxfp_recon_vs_true_foptics->Fill(xfp,xfpforward);
	hyfp_recon_vs_true_foptics->Fill(yfp,yfpforward);
	hxpfp_recon_vs_true_foptics->Fill(xpfp,xpfpforward);
	hypfp_recon_vs_true_foptics->Fill(ypfp,ypfpforward);
	
	tout->Fill();
      }
    }
  }

  elist->Delete();
  fout->Write();
}

