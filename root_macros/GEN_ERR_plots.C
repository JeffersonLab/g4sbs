#include "TH1D.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "G4SBSRunData.hh"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "gen_tree.C"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include "TCut.h"
#include "TEventList.h"

const double PI = TMath::Pi();
const double Mp = 0.938272;
const double Mn = 0.939565;
const double clight = 299792458.0; //m/s 

using namespace std;

void GEN_ERR_plots(const char *configfilename, const char *outputfilename="GEN_ERR_plots.root"){

  TFile *fout = new TFile(outputfilename,"RECREATE");
  
  TChain *C = new TChain("T");

  ifstream configfile(configfilename);

  TString currentline;

  //Things we need to configure:

  // 1) Optics file for reconstruction
  // 2) BigBite tracker pitch angle (we could be lazy and hardcode it)

  //Things we can extract from the run metadata:
  // Beam energy
  // SBS and BigBite angles
  // HCAL distance and vertical offset

  set<TString> list_of_files;
  map<TString,double> SBStheta_file;
  map<TString,double> BBtheta_file;
  map<TString,double> HCALdist_file;
  map<TString,double> HCALvoff_file;
  map<TString,double> HCALhoff_file; //for GEN we can ignore this safely:
  map<TString,double> Ebeam_file;
  //For normalization:
  map<TString,int> Ngen_file;
  map<TString,double> Lumi_file;
  map<TString,double> Genvol_file;

  double BBtheta = 34.0 * PI/180.0;
  double SBStheta = 17.5 * PI/180.0;
  double HCALdist = 17.0; //m
  double HCALvoff = 0.0;
  double HCALhoff = 0.0;
  double Ebeam = 8.8; //GeV
  double Lumi = 1.8e37; //electrons/s * nucleons/cm2
  double GenVol = 0.17; //solid angle generated:
  
  double ngen_total = 0.0;
  //Parse list of files:
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      //Check each file for the existence of the run_data object:
      TFile *ftemp = new TFile(currentline,"READ");

      G4SBSRunData *rd;

      ftemp->GetObject("run_data",rd);

      if( rd ){
	C->Add(currentline);
	list_of_files.insert(currentline);
	SBStheta_file[currentline] = rd->fSBStheta;
	BBtheta_file[currentline] = rd->fBBtheta;
	HCALdist_file[currentline] = rd->fHCALdist;
	HCALvoff_file[currentline] = rd->fHCALvoff;
	HCALhoff_file[currentline] = rd->fHCALhoff;
	Ebeam_file[currentline] = rd->fBeamE;
	Ngen_file[currentline] = rd->fNtries;
	//Increment total number of events:
	ngen_total += rd->fNtries;
	
	Lumi_file[currentline] = rd->fLuminosity;
	Genvol_file[currentline] = rd->fGenVol;

	//Initialize defaults:
	BBtheta = rd->fBBtheta;
	SBStheta = rd->fSBStheta;
	HCALdist = rd->fHCALdist;
	HCALvoff = rd->fHCALvoff;
	HCALhoff = rd->fHCALhoff;
	Ebeam = rd->fBeamE;
	Lumi = rd->fLuminosity;
	GenVol = rd->fGenVol;
      }

      ftemp->Close();
      ftemp->Delete();
    }
  }

  TCut globalcut = "";
  
  while( currentline.ReadLine(configfile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }
  }
  
  double BBtrackerpitch = 10.0*PI/180.0;
  double cut_W2 = 1.5; //GeV^2
  double cut_pmissperp = 0.1; //GeV
  double cut_pmisspar  = 1.5; //GeV
  double cut_thetapq   = 10.0*PI/180.0; //10 degrees;
  double cut_beta      = 0.05; //beta(TOF)-betaexpected (electron)
  double cut_En        = 0.5; //En_recon/En_expect - 1
  
  double gain_cor = 0.895; //constant of proportionality between BBcal sum and reconstructed P.
  
  vector<double> xptar_coeff, yptar_coeff, ytar_coeff, ptheta_coeff;
  vector<int> xfpexpon,yfpexpon,xpfpexpon,ypfpexpon,xtarexpon,sumexpon;

  //for first-order optics calculations, we are interested in computing
  // x'(x_5, E), y'(y_5,E)
  // x_5 = x + x' z5
  // xptar ~= -.00128 + 0.532*x -0.406 * x'
  // yptar ~= -1.3e-4 -0.05 y + 1.023 y'
  // ytar ~= -3.9e-4 + 1.12 y -2.69 y' ~= 0 --> 
  // pthetabend ~= 0.286 + 0.157 x -0.161 x'
  // thetabend = pthetabend/E ~= (10 deg. + xptar - x')
  // --> 0.286 + 0.157 x - 0.161 x' = E*(10 deg. + 
  double Cxptar_x,Cxptar_xp;
  double Cyptar_y,Cyptar_yp;
  double Cytar_y, Cytar_yp;
  double Cpth_x, Cpth_xp;
  double xptar_0, yptar_0, ytar_0;
  double pth_0; //constant term;

  fout->cd();
  TH1D *hxpdiff_firstorder = new TH1D("hxpdiff_firstorder","",250,-.05,.05);
  TH1D *hypdiff_firstorder = new TH1D("hypdiff_firstorder","",250,-.05,.05);
  
  TH2D *hdx1storder_vs_plane = new TH2D("hdx1storder_vs_plane","",4,0.5,4.5,250,-0.075,0.075);
  TH2D *hdy1storder_vs_plane = new TH2D("hdy1storder_vs_plane","",4,0.5,4.5,250,-0.075,0.075);
  
  TH1D *hEsum_HCAL = new TH1D("hEsum_HCAL","",250,0.0,1.0);
  TH1D *hEoverP_BBCAL = new TH1D("hEoverP_BBCAL","",250,0.0,2.0);

  TH1D *hW2_all_proton = new TH1D("hW2_all_proton","",250,-2.0,7.0);
  TH1D *hW2_cut_proton = new TH1D("hW2_cut_proton","",250,-2.0,7.0);
  TH1D *hW2_anticut_proton = new TH1D("hW2_anticut_proton","",250,-2.0,7.0);
  
  TH1D *hW2_all_neutron = new TH1D("hW2_all_neutron","",250,-2.0,7.0);
  TH1D *hW2_cut_neutron = new TH1D("hW2_cut_neutron","",250,-2.0,7.0);
  TH1D *hW2_anticut_neutron = new TH1D("hW2_anticut_neutron","",250,-2.0,7.0);

  

  TH1D *hpmisspar_all_proton= new TH1D("hpmisspar_all_proton","",250,-2.,2.);
  TH1D *hpmissperp_all_proton = new TH1D("hpmissperp_all_proton","",250,0.0,2.0);

  TH1D *hpmisspar_cut_proton= new TH1D("hpmisspar_cut_proton","",250,-2.,2.);
  TH1D *hpmissperp_cut_proton = new TH1D("hpmissperp_cut_proton","",250,0.0,2.0);
  
  TH1D *hpmisspar_all_neutron= new TH1D("hpmisspar_all_neutron","",250,-2.,2.);
  TH1D *hpmissperp_all_neutron = new TH1D("hpmissperp_all_neutron","",250,0.0,2.0);

  TH1D *hpmisspar_cut_neutron= new TH1D("hpmisspar_cut_neutron","",250,-2.,2.);
  TH1D *hpmissperp_cut_neutron = new TH1D("hpmissperp_cut_neutron","",250,0.0,2.0);

  TH2D *hpmisspar_vs_W2_proton = new TH2D("hpmisspar_vs_W2_proton","",200,-2.0,7.0,200,-2.0,2.0);
  TH2D *hpmisspar_vs_W2_neutron = new TH2D("hpmisspar_vs_W2_neutron","",200,-2.0,7.0,200,-2.0,2.0);

  TH2D *hpmissperp_vs_W2_proton = new TH2D("hpmissperp_vs_W2_proton","",200,-2.0,7.0,200,0.0,2.0);
  TH2D *hpmissperp_vs_W2_neutron = new TH2D("hpmissperp_vs_W2_neutron","",200,-2.0,7.0,200,0.0,2.0);

  //missing momenta using smeared neutron momentum via TOF
  
  TH1D *hpmissparsm_neutron = new TH1D("hpmissparsm_neutron","",250,-4.,4.);
  TH1D *hpmissparsm_proton = new TH1D("hpmissparsm_proton","",250,-4.,4.);

  TH1D *hpmissperpsm_neutron = new TH1D("hpmissperpsm_neutron","",250,-4.,4.);
  TH1D *hpmissperpsm_proton = new TH1D("hpmissperpsm_proton","",250,-4.,4.);

  TH1D *hnhitclust_HCAL = new TH1D("hnhitclust_HCAL","",50,0.5,50.5);
  TH1D *hEclust_HCAL_p = new TH1D("hEclust_HCAL_max_p","",250,0.0,1.0);
  TH1D *hEclust_HCAL_n = new TH1D("hEclust_HCAL_max_n","",250,0.0,1.0);
  TH1D *hnclust_HCAL = new TH1D("hnclust_HCAL","",21,-0.5,20.5);

  //
  TH1D *hEclust_HCAL_all_p = new TH1D("hEclust_HCAL_all_p","",250,0.0,1.0);
  TH1D *hEclust_HCAL_all_n = new TH1D("hEclust_HCAL_all_n","",250,0.0,1.0);
  
  
  TH1D *hdt_HCAL = new TH1D("hdt_HCAL","",250,-10.0,10.0);
  TH1D *hdx_HCAL = new TH1D("hdx_HCAL","",250,-0.3,0.3);
  TH1D *hdy_HCAL = new TH1D("hdy_HCAL","",250,-0.3,0.3);

  TH1D *hbeta_HCAL_p = new TH1D("hbeta_HCAL_p","",250,0.0,2.0);

  TH1D *hdbeta_HCAL_p = new TH1D("hdbeta_HCAL_p","",250,-0.1,0.1);

  TH1D *hbeta_HCAL_n = new TH1D("hbeta_HCAL_n","",250,0.0,2.0);

  TH1D *hdbeta_HCAL_n = new TH1D("hdbeta_HCAL_n","",250,-0.1,0.1);

  TH1D *hsinthetapq_p = new TH1D("hsinthetapq_p","",250,0.0,0.5);
  TH1D *hsinthetapq_n = new TH1D("hsinthetapq_n","",250,0.0,0.5);

  TH1D *hEoverP_HCAL = new TH1D("hEoverP_HCAL","",250,0.0,0.5);
  TH1D *hErecon_HCAL = new TH1D("hErecon_HCAL","",250,0.0,12.0);
  TH1D *hdErecon_HCAL = new TH1D("hdErecon_HCAL","",250,-1.0,1.0);
  
  //hpmisspar_vs_W2_proton->SetBit(TH1::kCanRebin, false);
  
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
    TObjArray *tokens = currentline.Tokenize(" ");

    int ntokens = tokens->GetEntries();
    if( ntokens > 1 ){
      TString skey = ( (TObjString*) (*tokens)[0] )->GetString();

      if( skey == "cut_beta" ){
	TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	cut_beta = sval.Atof();
      }
      
      if( skey == "cut_W2" ){
	TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	cut_W2 = sval.Atof();
      }

      if( skey == "cut_pmissperp" ){
	TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	cut_pmissperp = sval.Atof();
      }

      if( skey == "cut_pmisspar" ){
	TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	cut_pmisspar = sval.Atof();
      }

      if( skey == "cut_thetapq" ){
	TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	cut_thetapq = sval.Atof() * PI/180.0;
      }
      
      if( skey == "BBtrackerpitch" ){ ///BigBite tracker pitch angle, assumed to be given in degrees:
	TString sval = ( (TObjString*) (*tokens)[1] )->GetString();
	BBtrackerpitch = sval.Atof()*PI/180.0;
      }

      if( skey == "BBopticsfile" ){ //BigBite optics model
	TString sfname = ( (TObjString*) (*tokens)[1] )->GetString();
	ifstream opticsfile(sfname.Data() );

	if( opticsfile ) cout << "Opened file " << sfname << endl;

	int nterms;
	opticsfile >> nterms;
						   
	cout << "nterms = " << nterms << endl;
						  
						   
	xptar_coeff.resize(nterms);
	yptar_coeff.resize(nterms);
	ytar_coeff.resize(nterms);
	ptheta_coeff.resize(nterms);

	xfpexpon.resize(nterms);
	yfpexpon.resize(nterms);
	xpfpexpon.resize(nterms);
	ypfpexpon.resize(nterms);
	xtarexpon.resize(nterms);
	sumexpon.resize(nterms);
	
	for( int iterm=0; iterm<nterms; iterm++ ){
	  opticsfile >> xptar_coeff[iterm] >> yptar_coeff[iterm] >> ytar_coeff[iterm] >> ptheta_coeff[iterm]
		     >> xfpexpon[iterm] >> yfpexpon[iterm] >> xpfpexpon[iterm] >> ypfpexpon[iterm] >> xtarexpon[iterm];

	  //This will help in choosing the first-order terms:
	  sumexpon[iterm] = xfpexpon[iterm]+yfpexpon[iterm]+xpfpexpon[iterm]+ypfpexpon[iterm]+xtarexpon[iterm];

	  
	  
	  //Collect constant and first-order terms:
	  if( sumexpon[iterm] == 0 ){
	    xptar_0 = xptar_coeff[iterm];
	    yptar_0 = yptar_coeff[iterm];
	    ytar_0 = ytar_coeff[iterm];
	    pth_0 = ptheta_coeff[iterm];
	  }
	  if( sumexpon[iterm] == 1 ){
	    if( xfpexpon[iterm] == 1 ){ //x terms:
	      Cxptar_x = xptar_coeff[iterm];
	      Cpth_x = ptheta_coeff[iterm];
	    }
	    if( xpfpexpon[iterm] == 1 ){
	      Cxptar_xp = xptar_coeff[iterm];
	      Cpth_xp = ptheta_coeff[iterm];
	    }
	    if( yfpexpon[iterm] == 1 ){
	      Cyptar_y = yptar_coeff[iterm];
	      Cytar_y = ytar_coeff[iterm]; 
	    }

	    if( ypfpexpon[iterm] == 1 ){
	      Cyptar_yp = yptar_coeff[iterm];
	      Cytar_yp = ytar_coeff[iterm]; 
	    }

	  }
	}
	
      }
      
    }
    
    delete tokens;
  }

  TEventList *elist = new TEventList("elist");

  C->Draw(">>elist",globalcut);
  
  gen_tree *T = new gen_tree(C);

  long nevent=0;

  int oldtreenum = -1;

  TRandom3 num(0);
  
  while( C->GetEntry(elist->GetEntry(nevent++)) ){
    if( nevent % 1000 == 0 ) cout << "Event " << nevent << endl;

    int treenum = C->GetTreeNumber();
    
    if( treenum != oldtreenum ){ //new tree or first event: update constant parameters:
      TString fname = C->GetFile()->GetName();
      BBtheta = BBtheta_file[fname];
      SBStheta = SBStheta_file[fname];
      HCALdist = HCALdist_file[fname];
      HCALvoff = HCALvoff_file[fname];
      HCALhoff = HCALhoff_file[fname];
      Ebeam = Ebeam_file[fname];
      Lumi = Lumi_file[fname];
      GenVol = Genvol_file[fname];

      oldtreenum = treenum;
      
      cout << "New file " << fname << endl;
    }

    double weight = T->ev_sigma * Lumi * GenVol / ngen_total;

    TVector3 BB_zaxis( sin(BBtheta), 0, cos(BBtheta) );
    TVector3 BB_xaxis( 0, -1, 0 );
    TVector3 BB_yaxis = (BB_zaxis.Cross(BB_xaxis)).Unit();

    TVector3 BBGEM_zaxis = BB_zaxis;
    TVector3 BBGEM_yaxis = BB_yaxis;
    BBGEM_zaxis.Rotate( -BBtrackerpitch, BBGEM_yaxis );
    TVector3 BBGEM_xaxis = BBGEM_yaxis.Cross(BBGEM_zaxis).Unit();

    TVector3 HCAL_zaxis( -sin(SBStheta),0,cos(SBStheta));
    TVector3 HCAL_yaxis( 0,1,0 );
    TVector3 HCAL_xaxis = HCAL_yaxis.Cross(HCAL_zaxis).Unit();
    
    double EPS = T->Earm_BBPSTF1_det_esum;
    double ESH = T->Earm_BBSHTF1_det_esum;

    if( EPS + ESH > 0.0 && T->Earm_BBGEM_Track_ntracks==1&&(*(T->Earm_BBGEM_Track_MID))[0]==0 ){
    
      double npemean_PS=200.*EPS;
      double npemean_SH=300.*ESH;
      
      double EPS_smear = num.Gaus( npemean_PS, sqrt(npemean_PS) )/200.;
      double ESH_smear = num.Gaus( npemean_SH, sqrt(npemean_SH) )/300.;
      
      double Esum_BBCAL = EPS_smear + ESH_smear;

      double xfp = (*(T->Earm_BBGEM_Track_Xfit))[0];
      double yfp = (*(T->Earm_BBGEM_Track_Yfit))[0];
      double xpfp = (*(T->Earm_BBGEM_Track_Xpfit))[0];
      double ypfp = (*(T->Earm_BBGEM_Track_Ypfit))[0];

      double xptar_recon = 0.0, yptar_recon = 0.0, ytar_recon = 0.0, pthetabend_recon = 0.0;
      double xtar_recon = 0.0;

      //target reconstruction from smeared FP quantities:
      for( int iterm = 0; iterm<sumexpon.size(); iterm++ ){
	double term =
	  pow(xfp,xfpexpon[iterm])*
	  pow(yfp,yfpexpon[iterm])*
	  pow(xpfp,xpfpexpon[iterm])*
	  pow(ypfp,ypfpexpon[iterm])*
	  pow(xtar_recon,xtarexpon[iterm]);

	xptar_recon += xptar_coeff[iterm]*term;
	yptar_recon += yptar_coeff[iterm]*term;
	ytar_recon += ytar_coeff[iterm]*term;
	pthetabend_recon += ptheta_coeff[iterm]*term;
      }

      // cout << "pthetabend_recon = " << pthetabend_recon << endl;
      // cout << "xptar recon, yptar recon, ytar recon = " << xptar_recon << ", " << yptar_recon << ", " << ytar_recon << endl;

      //used later to calculate nucleon TOF:
      double vz_recon = -ytar_recon/(sin(BBtheta)+atan(yptar_recon)*cos(BBtheta) );
      
      TVector3 phat_tgt_local(xptar_recon,yptar_recon,1.0);
      phat_tgt_local = phat_tgt_local.Unit();
      TVector3 phat_tgt =
	phat_tgt_local.X() * BB_xaxis +
	phat_tgt_local.Y() * BB_yaxis +
	phat_tgt_local.Z() * BB_zaxis; 
      
      TVector3 phat_fp_local(xpfp,ypfp,1.0);
      phat_fp_local = phat_fp_local.Unit();

      TVector3 phat_fp =
	phat_fp_local.X() * BBGEM_xaxis +
	phat_fp_local.Y() * BBGEM_yaxis +
	phat_fp_local.Z() * BBGEM_zaxis;

      double thetabend_recon = acos( phat_fp.Dot( phat_tgt ) );

      double p_recon = pthetabend_recon/thetabend_recon;
      
      //cout << "p_recon, ptrue = " << p_recon << ", " << (*(T->Earm_BBGEM_Track_P))[0] << endl;
      
      hEoverP_BBCAL->Fill( Esum_BBCAL/p_recon, weight );
      if( T->Harm_HCalScint_det_esum > 0 ) hEsum_HCAL->Fill( T->Harm_HCalScint_det_esum, weight );

      TVector3 ep = p_recon*phat_tgt;
      TLorentzVector eP(ep, p_recon);

      TLorentzVector Pbeam(0,0,Ebeam,Ebeam);

      double MN = Mp;
      if( T->ev_fnucl == 0 ) MN = Mn;
      
      TLorentzVector Pinucleon(0,0,0,MN);

      TLorentzVector q = Pbeam - eP;

      TLorentzVector Pfnucleon = Pinucleon + q;

      //Need to do nucleon reconstruction:
      //Find hit with max. energy deposit:
      set<int> unused_hits; //list of hits not yet assigned to a cluster by position in good hit array
      vector<int> hitlist_hcal; //list of hits above threshold by position in tree hit array
      //map<int,int> hitidx_hcal; //key = position in hit list, val = position in root tree hit array
      vector<double> xhit_hcal,yhit_hcal,Ehit_hcal,thit_hcal;
      vector<int> row_hcal,col_hcal;

      //For island algorithm, we need to be able to look up cells by row and column. So we need
      //a map indexed by unique cell number:
      map<int,int> hitcell_hcal; //key = cell, val = hit index in good hit array:
      
      double thresh_hcal = 0.007;

      int ngood=0;
      
      int imax_hcal=-1;
      double emax_hcal = 0.0;

      //If we want here we can "cheat" and only include hits attributable to the primary nucleon:

      TVector3 SDTrack_pos; //global position of SD track boundary crossing:
      
      for( int ihit=0; ihit<T->Harm_HCalScint_hit_nhits; ihit++ ){
	double edep = (*(T->Harm_HCalScint_hit_sumedep))[ihit];
	int MID = (*(T->SDTrack_MID))[(*(T->Harm_HCalScint_hit_sdtridx))[ihit]];

	int PrTID = (*(T->PTrack_TID))[(*(T->Harm_HCalScint_hit_ptridx))[ihit]];

	int PID = (*(T->SDTrack_PID))[(*(T->Harm_HCalScint_hit_sdtridx))[ihit]];

	int PrPID = (*(T->PTrack_PID))[(*(T->Harm_HCalScint_hit_ptridx))[ihit]];

	
	//if( edep >= thresh_hcal && MID == PrTID && PID == 22 && PrPID == 111 ){
	if( edep >= thresh_hcal ){
	  hitlist_hcal.push_back( ihit );
	  xhit_hcal.push_back( (*(T->Harm_HCalScint_hit_xcell))[ihit] );
	  yhit_hcal.push_back( (*(T->Harm_HCalScint_hit_ycell))[ihit] );
	  Ehit_hcal.push_back( edep );
	  thit_hcal.push_back( (*(T->Harm_HCalScint_hit_tavg))[ihit] );
	  //hitused_hcal.push_back( false );
	  row_hcal.push_back( (*(T->Harm_HCalScint_hit_row))[ihit] );
	  col_hcal.push_back( (*(T->Harm_HCalScint_hit_col))[ihit] );
	  hitcell_hcal[(*(T->Harm_HCalScint_hit_cell))[ihit]] = ngood;
	  if( edep > emax_hcal ){
	    imax_hcal = ngood;
	    emax_hcal = edep;
	    SDTrack_pos.SetXYZ( (*(T->SDTrack_posx))[(*(T->Harm_HCalScint_hit_sdtridx))[ihit] ],
				(*(T->SDTrack_posy))[(*(T->Harm_HCalScint_hit_sdtridx))[ihit] ],
				(*(T->SDTrack_posz))[(*(T->Harm_HCalScint_hit_sdtridx))[ihit] ] );
						     
	  }

	  // hitrow[(*(T->Harm_HCalScint_hit_row))[ihit]] = ngood;
	  // hitcol[(*(T->Harm_HCalScint_hit_col))[ihit]] = ngood;

	  unused_hits.insert(ngood);
	  
	  ngood++;
	}
      }

      int nclust=0;
      vector<int> nhitclust;
      vector<double> xclust, yclust, Eclust,tclust;
      
      if( imax_hcal >= 0 ){
	bool foundclust = true;
	while( foundclust ){
	  foundclust=false;

	  if( imax_hcal >= 0 ){

	    vector<int> hitsincluster;
	    double xsum = 0.0, ysum = 0.0, esum = 0.0, tsum = 0.0;
	    int nhittemp = 0;

	    tsum += Ehit_hcal[imax_hcal]*thit_hcal[imax_hcal];
	    xsum += Ehit_hcal[imax_hcal]*xhit_hcal[imax_hcal];
	    ysum += Ehit_hcal[imax_hcal]*yhit_hcal[imax_hcal];
	    esum += Ehit_hcal[imax_hcal];

	    nhittemp++;

	    unused_hits.erase( imax_hcal );

	    int rowmax = row_hcal[imax_hcal], colmax=col_hcal[imax_hcal];

	    int cellmax = colmax + 12*rowmax;

	    hitsincluster.push_back(imax_hcal);
	  
	    int ihittemp=0;
	    while( ihittemp < nhittemp ){ //island algorithm for clustering hits:
	      //search left, right, up, down until we run out of unused hits
	      //contiguous with any hit already added to the cluster:
	      int rowtemp = row_hcal[hitsincluster[ihittemp]];
	      int coltemp = col_hcal[hitsincluster[ihittemp]];

	      int rowtest=rowtemp;
	      int coltest=coltemp-1;
	      int celltest = coltest + 12*rowtest;

	      std::map<int,int>::iterator foundhit = hitcell_hcal.find(celltest);

	      if( foundhit != hitcell_hcal.end() && coltest >= 0 && coltest <12 && rowtest >= 0 && rowtest < 24){ //new hit
	      
		int hitidx = foundhit->second;
		if( unused_hits.find(hitidx) != unused_hits.end() ){
	      
		  hitsincluster.push_back( hitidx );
		  nhittemp++;
		  xsum += Ehit_hcal[hitidx]*xhit_hcal[hitidx];
		  ysum += Ehit_hcal[hitidx]*yhit_hcal[hitidx];
		  tsum += Ehit_hcal[hitidx]*thit_hcal[hitidx];
		  esum += Ehit_hcal[hitidx];

		  unused_hits.erase(hitidx);
		}
	      }

	      coltest=coltemp+1;
	      celltest=coltest + 12*rowtest;

	      foundhit = hitcell_hcal.find(celltest);

	      if( foundhit != hitcell_hcal.end() && coltest >= 0 && coltest <12 && rowtest >= 0 && rowtest < 24){ //new hit
	      
		int hitidx = foundhit->second;
		if( unused_hits.find(hitidx) != unused_hits.end() ){
	      
		  hitsincluster.push_back( hitidx );
		  nhittemp++;
		  xsum += Ehit_hcal[hitidx]*xhit_hcal[hitidx];
		  ysum += Ehit_hcal[hitidx]*yhit_hcal[hitidx];
		  tsum += Ehit_hcal[hitidx]*thit_hcal[hitidx];
		  esum += Ehit_hcal[hitidx];

		  unused_hits.erase(hitidx);
		}
	      }

	      coltest=coltemp;
	      rowtest=rowtemp-1;
	      celltest=coltest + 12*rowtest;

	      foundhit = hitcell_hcal.find(celltest);

	      if( foundhit != hitcell_hcal.end() && coltest >= 0 && coltest <12 && rowtest >= 0 && rowtest < 24){ //new hit
	      
		int hitidx = foundhit->second;
		if( unused_hits.find(hitidx) != unused_hits.end() ){
	      
		  hitsincluster.push_back( hitidx );
		  nhittemp++;
		  xsum += Ehit_hcal[hitidx]*xhit_hcal[hitidx];
		  ysum += Ehit_hcal[hitidx]*yhit_hcal[hitidx];
		  tsum += Ehit_hcal[hitidx]*thit_hcal[hitidx];
		  esum += Ehit_hcal[hitidx];

		  unused_hits.erase(hitidx);
		}
	      }

	      coltest=coltemp;
	      rowtest=rowtemp+1;
	      celltest=coltest + 12*rowtest;

	      foundhit = hitcell_hcal.find(celltest);

	      if( foundhit != hitcell_hcal.end() && coltest >= 0 && coltest <12 && rowtest >= 0 && rowtest < 24 ){ //new hit
	      
		int hitidx = foundhit->second;
		if( unused_hits.find(hitidx) != unused_hits.end() ){
	      
		  hitsincluster.push_back( hitidx );
		  nhittemp++;
		  xsum += Ehit_hcal[hitidx]*xhit_hcal[hitidx];
		  ysum += Ehit_hcal[hitidx]*yhit_hcal[hitidx];
		  tsum += Ehit_hcal[hitidx]*thit_hcal[hitidx];
		  esum += Ehit_hcal[hitidx];

		  unused_hits.erase(hitidx);
		}
	      }
	      ihittemp++; //go to next available hit.
	    } //while (ihittemp < nhittemp)

	    foundclust = true;

	    nhitclust.push_back( nhittemp );
	    xclust.push_back( xsum/esum );
	    yclust.push_back( ysum/esum );
	    tclust.push_back( tsum/esum );
	    Eclust.push_back( esum );
	    nclust++;
	    
	  } //if (imax_hcal >= 0 ) subsequent maxima

	  imax_hcal = -1;
	  emax_hcal = 0.0;
	  //Loop over unused hits and see if we can find another maximum:
	  for(set<int>::iterator ihit=unused_hits.begin(); ihit!=unused_hits.end(); ++ihit ){
	    int hitidx = *ihit;
	    if( Ehit_hcal[hitidx] > emax_hcal ){
	      emax_hcal = Ehit_hcal[hitidx];
	      imax_hcal = hitidx;
	    }
	  }
      	  
	} //while (foundclust)
      } //if(imax_hcal >= 0) first maximum

      if( nclust > 0 ){ //choose cluster with maximum energy:
	int imax_clust=0;
	double Emax_clust = 0.0;
	for( int iclust=0; iclust<nclust; iclust++ ){
	  if( Eclust[iclust] > Emax_clust ){
	    imax_clust = iclust;
	    Emax_clust = Eclust[iclust];
	  }
	}

	hnhitclust_HCAL->Fill( nhitclust[imax_clust], weight );
	hnclust_HCAL->Fill( nclust, weight );
	
	//double sampling_frac = 0.08;

	//Now, since the energy measurement is almost useless, let's calculate the nucleon
	//momentum via TOF, ASSUMING a neutron:

	TVector3 hitpos(xclust[imax_clust],yclust[imax_clust],HCALdist);

	TVector3 hitpos_global = hitpos.X() * HCAL_xaxis + hitpos.Y() * HCAL_yaxis + hitpos.Z() * HCAL_zaxis;

	TVector3 HCAL_center_pos = HCALdist * HCAL_zaxis;

	TVector3 SDpos_rel = SDTrack_pos - HCAL_center_pos;

	double SDx = SDpos_rel.Dot(HCAL_xaxis);
	double SDy = SDpos_rel.Dot(HCAL_yaxis);

	hdx_HCAL->Fill( xclust[imax_clust] - SDx, weight );
	hdy_HCAL->Fill( yclust[imax_clust] - SDy, weight );

	//Include resolution of reconstructed z vertex:
	TVector3 vertex_global( T->ev_vx, T->ev_vy, vz_recon ); 

	TVector3 FlightPath = hitpos_global - vertex_global;

	double PathLength_Neutron = FlightPath.Mag(); //in meters

	double TOF = tclust[imax_clust]-2.83; //in ns

	double v_neutron = PathLength_Neutron/TOF*1.e9;
	double beta_neutron = v_neutron/clight;

	//Identify HCAL energy with the kinetic energy:
	double En_recon = Eclust[imax_clust]/0.0685 + Mn;
	double pn_recon = sqrt(pow(En_recon,2)-pow(Mn,2));

	double Tneutron_true = sqrt(pow(T->ev_np,2)+pow(Mn,2))-Mn;
	
	if( T->ev_fnucl == 0 ){
	  hbeta_HCAL_n->Fill( beta_neutron, weight );

	  hEclust_HCAL_n->Fill( Eclust[imax_clust], weight );

	  hEoverP_HCAL->Fill( Eclust[imax_clust]/Tneutron_true, weight );
	  hErecon_HCAL->Fill( En_recon, weight );
	  //hdErecon_HCAL->Fill( En_recon/(Tneutron_true+Mn)-1.0, weight );
	  for( int iclust=0; iclust<nclust; iclust++ ){
	    if( iclust != imax_clust ){
	      hEclust_HCAL_all_n->Fill( Eclust[iclust], weight );
	    }
	  }
	} else {
	  hbeta_HCAL_p->Fill( beta_neutron, weight );
	  hEclust_HCAL_p->Fill( Eclust[imax_clust], weight );
	  
	  hEoverP_HCAL->Fill( Eclust[imax_clust]/Tneutron_true, weight );
	  hErecon_HCAL->Fill( En_recon, weight );
	  //hdErecon_HCAL->Fill( En_recon/(Tneutron_true+Mn)-1.0, weight );
	  
	  for( int iclust=0; iclust<nclust; iclust++ ){
	    if( iclust != imax_clust ){
	      hEclust_HCAL_all_p->Fill( Eclust[iclust], weight );
	    }
	  }
	}
	
	//E = gamma m = m/sqrt(1-beta^2)
	//p = beta E
	// double gamma_neutron = 1.0/sqrt(1.0-pow(beta_neutron,2));
	// double p_neutron = gamma_neutron * beta_neutron * Mn;

	// //TVector3 Pvect_neutron = p_neutron * FlightPath.Unit();

	// double E_neutron = gamma_neutron * Mn;

	// TVector3 Pn_TOF = p_neutron * FlightPath.Unit();
	// TLorentzVector Pneutron_TOF(Pn_TOF, E_neutron);
	
	//	TLorentzVector Pneutron( Pvect_neutron, E_neutron );

	//Since the TOF resolution isn't really going to be sufficient to get the neutron energy/momentum right, let's try something different
	//Use the electron kinematics to predict the neutron momentum assuming elastic scattering on free neutron at rest:
	double Q2_recon = -q.M2();

	//	cout << "Q2 recon, Q2 = " << Q2_recon << ", " << T->ev_Q2 << endl;
	
	double pnexpect = sqrt(Q2_recon*(1.+Q2_recon/(4.*pow(Mn,2))));
	
	//	cout << "pnexpect, np = " << pnexpect << ", " << T->ev_np << endl;

	TVector3 Pvect_n = pnexpect * FlightPath.Unit();

	// cout << "thHCAL, nth = " << Pvect_n.Theta()*180.0/PI << ", " << T->ev_nth * 180.0/PI << endl;
	// cout << "phHCAL, nph = " << Pvect_n.Phi()*180.0/PI << ", " << T->ev_nph * 180.0/PI << endl;
	
	TLorentzVector Pneutron( Pvect_n, sqrt(pow(pnexpect,2)+pow(Mn,2)));

	double beta_expect = Pneutron.Beta();

	//	double beta_true = 
	
	if( T->ev_fnucl == 0 ){
	  hdbeta_HCAL_n->Fill( beta_neutron - beta_expect, weight );
	  hdErecon_HCAL->Fill( En_recon/Pneutron.E()-1.0, weight );
	} else {
	  hdbeta_HCAL_p->Fill( beta_neutron - beta_expect, weight );
	  hdErecon_HCAL->Fill( En_recon/Pneutron.E()-1.0, weight );
	}
	  
	double TOFexpect = PathLength_Neutron/beta_expect/clight*1.e9;

	hdt_HCAL->Fill( TOF - TOFexpect, weight );
	
	TLorentzVector Pmiss_4vect = Pneutron - q;

	TVector3 Pmiss = Pmiss_4vect.Vect();

	double pmiss_par = Pmiss.Dot( (q.Vect().Unit() ) );

	double pmiss_perp = (Pmiss - pmiss_par*q.Vect().Unit()).Mag();

	double W2_recon = Pfnucleon.M2();

	double thetapq = acos( FlightPath.Unit().Dot( q.Vect().Unit() ) );

	double sinthetapq = sin(thetapq);
	
	if( T->ev_fnucl == 0 ){ //neutron
	  hW2_all_neutron->Fill( W2_recon, weight );
	  hsinthetapq_n->Fill(sinthetapq, weight );
	  if( fabs( pmiss_par ) <= cut_pmisspar && pmiss_perp <= cut_pmissperp && fabs( beta_neutron - beta_expect) <= cut_beta && thetapq <= cut_thetapq && fabs( En_recon/Pneutron.E() -1.0 ) <= cut_En ) {
	    hW2_cut_neutron->Fill( W2_recon, weight );
	  } else {
	    hW2_anticut_neutron->Fill( W2_recon, weight );
	  }

	  hpmisspar_all_neutron->Fill( pmiss_par, weight );
	  hpmissperp_all_neutron->Fill( pmiss_perp, weight );

	  hpmisspar_vs_W2_neutron->Fill( W2_recon, pmiss_par, weight );
	  hpmissperp_vs_W2_neutron->Fill( W2_recon, pmiss_perp, weight );
	  
	  if( W2_recon < cut_W2 ){
	    if( pmiss_perp <= cut_pmissperp ){
	      hpmisspar_cut_neutron->Fill( pmiss_par, weight );
	    }
	    if( fabs( pmiss_par ) <= cut_pmisspar ){
	      hpmissperp_cut_neutron->Fill( pmiss_perp, weight );
	    }
	  }
	  
	  
	} else { //proton
	  hW2_all_proton->Fill( W2_recon, weight );
	  hsinthetapq_p->Fill(sinthetapq, weight );
	  if( fabs( pmiss_par ) <= cut_pmisspar && pmiss_perp <= cut_pmissperp && fabs( beta_neutron - beta_expect) <= cut_beta && thetapq <= cut_thetapq && fabs( En_recon/Pneutron.E() -1.0 ) <= cut_En ){
	    hW2_cut_proton->Fill( W2_recon, weight );
	  } else {
	    hW2_anticut_proton->Fill( W2_recon, weight );
	  }

	  hpmisspar_all_proton->Fill( pmiss_par, weight );
	  hpmissperp_all_proton->Fill( pmiss_perp, weight );

	  hpmisspar_vs_W2_proton->Fill( W2_recon, pmiss_par, weight );
	  hpmissperp_vs_W2_proton->Fill( W2_recon, pmiss_perp, weight );
	  
	  if( W2_recon < cut_W2 ){
	    if( pmiss_perp <= cut_pmissperp ){
	      hpmisspar_cut_proton->Fill( pmiss_par, weight );
	    }
	    if( fabs( pmiss_par ) <= cut_pmisspar ){
	      hpmissperp_cut_proton->Fill( pmiss_perp, weight );
	    }
	  }
	  
	}
	
      }
      
      //Calculate missing parallel and perp momenta:
      
      
      
      
      bool hitlayer5 = false;

      double x5,y5,z5;
      //Now, about computing thost first-order xp and yp diffs:
      double xfGEM[4],yfGEM[4],zfGEM[4];
      
      for( int ihit=0; ihit<T->Earm_BBGEM_hit_nhits; ihit++ ){
	//check for hits in the fifth layer caused by primary electron:
	if( (*(T->Earm_BBGEM_hit_mid))[ihit] == 0 ){
	  if( (*(T->Earm_BBGEM_hit_plane))[ihit] == 5 ){
	    hitlayer5 = true;
	    x5 = (*(T->Earm_BBGEM_hit_x))[ihit];
	    y5 = (*(T->Earm_BBGEM_hit_y))[ihit];
	    z5 = (*(T->Earm_BBGEM_hit_z))[ihit];
	    //Let's assume 100 um smearing = 0.1 mm = 1e-4 m 
	    x5 = num.Gaus(x5,1.e-4);
	    y5 = num.Gaus(y5,1.e-4);
	  } else {
	    int plane = (*(T->Earm_BBGEM_hit_plane))[ihit];
	    xfGEM[plane-1] = num.Gaus( (*(T->Earm_BBGEM_hit_x))[ihit], 1.e-4);
	    yfGEM[plane-1] = num.Gaus( (*(T->Earm_BBGEM_hit_y))[ihit], 1.e-4);
	    zfGEM[plane-1] = (*(T->Earm_BBGEM_hit_z))[ihit];
	  }
	}
      }

      if( hitlayer5 ){ //then attempt to calculate x' and y' first first order optics:
	double ECAL = Esum_BBCAL/gain_cor;
	// xptar = xptar_0 + Cxptar_x * (x5 - z5 x') + Cxptar_xp * x'
	// thetabend = pthetabend/E = pth_0/E + Cpth_x * (x5 - z5 x') + Cpth_xp * x'
	// thetabend = 10 deg. + xptar - x'
	// 10 deg. + xptar - x' = Cpth_x/E * (x5-z5 x') + Cpth_xp/E x'
	// xptar = xptar_0 + Cxptar_x * (x5 - z5 x') + Cxptar_xp * x'
	// 10 deg. + [xptar_0 + Cxpx * (x5 - z5 x')  + Cxpxp * x'] - x' = Cpth_x/E * (x5 - z5 x') + Cpth_xp/E*x'
	double xp_firstorder = (BBtrackerpitch + xptar_0 + Cxptar_x * x5 - Cpth_x/ECAL * x5 - pth_0/ECAL )/
	  (1.0 + Cpth_xp/ECAL - Cpth_x/ECAL*z5 + Cxptar_x * z5 - Cxptar_xp );

	//cout << "x' first order, x' true = " << xp_firstorder << ", " << xpfp << endl;

	//y5 = y + y' z5
	//For ytar = 0 assumption, we have y' Cyy' +  (y5- y'z5) Cyy + y0 ~= 0
	// y'(Cyy' - Cyy z5) = -(y0 + Cyy y5)
	
	double yp_firstorder = -(ytar_0 + Cytar_y * y5)/(Cytar_yp - Cytar_y * z5);
	
	hxpdiff_firstorder->Fill( xp_firstorder - xpfp, weight );
	hypdiff_firstorder->Fill( yp_firstorder - ypfp, weight );

	for( int iplane=1; iplane<=4; iplane++ ){
	  hdx1storder_vs_plane->Fill( iplane, x5 + xp_firstorder * (zfGEM[iplane-1]-z5) - xfGEM[iplane-1], weight );
	  hdy1storder_vs_plane->Fill( iplane, y5 + yp_firstorder * (zfGEM[iplane-1]-z5) - yfGEM[iplane-1], weight );
	}
      }
    }   
  }

  elist->Delete();

  fout->Write();
  
}
