#include "g4sbs_a1n_tree.C"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TChainElement.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "G4SBSRunData.hh"

const double Mp = 0.938272;

void DIS_rates_DIS(const char *configfilename, const char *outfilename){
  
  
  TChain *C = new TChain("T");

  TString currentline;

  ifstream infile(configfilename);

  set<TString> files;

  map<TString, long> ngen_file;
  map<TString, long> ntries_file; 
  
  while( currentline.ReadLine( infile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      
      C->Add( currentline );
      
      files.insert(currentline);
    }
  }

  TObjArray *FileList = C->GetListOfFiles();
  TIter next(FileList);

  TChainElement *chEl = 0;

  set<TString> bad_file_list;

  long ntries=0;
  long ngen=0;
  
  G4SBSRunData *rd;
  
  while( (chEl=(TChainElement*)next() )){
    TFile newfile(chEl->GetTitle());
    newfile.GetObject("run_data",rd);
    if( rd ){
      if( files.find( chEl->GetTitle() ) != files.end() ){
	ngen += rd->fNthrown;
	ntries += rd->fNtries;
	ngen_file[chEl->GetTitle()] = rd->fNthrown;
	ntries_file[chEl->GetTitle()] = rd->fNtries;
      }
    } else {
      bad_file_list.insert( chEl->GetTitle());
    }
  }

  cout << "Total number of generated events = " << ngen << endl;
  cout << "Total number of attempts = " << ntries << endl;
  cout << "Generation efficiency = " << double(ngen)/double(ntries) << endl;
  //cout << "Total number of generated events, neutron = " << ngen_n << endl;
  
  TH1D::SetDefaultSumw2();
  TH2D::SetDefaultSumw2();
  
  double Emin_SBS = 1.5; //GeV
  double Emin_BB = 0.8;

  int binflag = 0, nbins;
  vector<double> bins;
  infile >> binflag >> nbins;
  bins.resize(nbins+1);
  if( binflag == 0 ){ //fixed-width bins in xbj:
    double min, max;
    infile >> min >> max;
    for( int bin=0; bin<=nbins; bin++ ){
      bins[bin] = min + bin*(max-min)/double(nbins);
    }
  } else { //variable-width bins in xbj:
    bins.resize(nbins+1);
    for( int bin=0; bin<=nbins; bin++ ){
      infile >> bins[bin];
    }
  }

  double Ibeam=50e-6; //Amps
  double rho_tgt=1.3e-3, Ltgt=55.0; //g/cm^3 and cm
  double N_A = 6.022e23; //avogadro's number
  double e = 1.602e-19; //electron charge
  double Mmol_3He = 3.016; //g/mol
  double Mmol_H2 = 1.008; //g/mol

  //For H2, assume pressure of 10.5 atm:
  double targpressure = 10.5 * 101325.0; // J / m^3

  // P = rho R T / M --> rho = M * P / (RT)
  double RT = 8.314 * 300.0; // J / mol 
  // 1 atm = 101325 N/m^2 = 101325 J/m^3
  
  int tgt_flag = 3; 
  
  infile >> Ibeam >> Ltgt >> tgt_flag;
  Ibeam *= 1.e-6; //assumed to be given in muA.
  //             (e-/s) *     (g/cm^2)   / (g/mol)  * (atoms/mol) = e- * atoms /cm^2/s 
  double Lumi = Ibeam/e * rho_tgt * Ltgt / Mmol_3He * N_A; //default scenario:

  switch( tgt_flag ){
  case 1: //H2:
    rho_tgt = targpressure/RT * Mmol_H2 * 2.0 / 1.e6; //g/cm^3: 2 atoms/molecule.
    Lumi = Ibeam/e * rho_tgt * Ltgt / Mmol_H2 * N_A;
    break;
  case 3: //"Helium-3":
  default:
    rho_tgt = targpressure/RT * Mmol_3He / 1.e6;
    Lumi = Ibeam/e * rho_tgt * Ltgt / Mmol_3He * N_A;
    break;
  }

  
  cout << "For an assumed beam current of " << Ibeam*1e6 << " muA and Ltgt = " << Ltgt << " cm, electron-nucleus luminosity = " << Lumi
       << endl;

  double Pbeam = 0.85;
  double Ptgt = 0.6;
  infile >> Pbeam >> Ptgt;
  infile >> Emin_SBS >> Emin_BB;
  
  TFile *fout = new TFile(outfilename,"RECREATE");
  
  TH1D *hrate_p_xbj_BB = new TH1D("hrate_p_xbj_BB","",nbins,&(bins[0]));
  TH1D *hrate_n_xbj_BB = new TH1D("hrate_n_xbj_BB","",nbins,&(bins[0]));
  TH1D *hrate_3He_xbj_BB = new TH1D("hrate_3He_xbj_BB","",nbins,&(bins[0]));

  TH1D *hrate_p_xbj_BB_fixed_width = new TH1D("hrate_p_xbj_BB_fixed_width","",50,0.0,1.0);
  TH1D *hrate_n_xbj_BB_fixed_width = new TH1D("hrate_n_xbj_BB_fixed_width","",50,0.0,1.0);
  TH1D *hrate_3He_xbj_BB_fixed_width = new TH1D("hrate_3He_xbj_BB_fixed_width","",50,0.0,1.0);
  
  TProfile *hQ2_xbj_BB_prof = new TProfile("hQ2_xbj_BB_prof","",nbins,&(bins[0]));
  TProfile *hEprime_xbj_BB_prof = new TProfile("hEprime_xbj_BB_prof","",nbins,&(bins[0]));
  TProfile *hxmean_xbj_BB_prof = new TProfile("hxmean_xbj_BB_prof","",nbins,&(bins[0]));
  TProfile *hW_xbj_BB_prof = new TProfile("hW_xbj_BB_prof","",nbins,&(bins[0]));
  TProfile *hy_xbj_BB_prof = new TProfile("hy_xbj_BB_prof","",nbins,&(bins[0]));
  TProfile *heps_xbj_BB_prof = new TProfile("heps_xbj_BB_prof","",nbins,&(bins[0]));
  
  TH2D *hQ2_xbj_BB = new TH2D("hQ2_xbj_BB","",100,0,1,100,1,12);
  TH2D *hEprime_xbj_BB = new TH2D("hEprime_xbj_BB","",100,0,1,100,0,11);

  TH1D *hrate_p_xbj_SBS = new TH1D("hrate_p_xbj_SBS","",nbins,&(bins[0]));
  TH1D *hrate_n_xbj_SBS = new TH1D("hrate_n_xbj_SBS","",nbins,&(bins[0]));
  TH1D *hrate_3He_xbj_SBS = new TH1D("hrate_3He_xbj_SBS","",nbins,&(bins[0]));
  
  TH1D *hrate_p_xbj_SBS_fixed_width = new TH1D("hrate_p_xbj_SBS_fixed_width","",50,0.0,1.0);
  TH1D *hrate_n_xbj_SBS_fixed_width = new TH1D("hrate_n_xbj_SBS_fixed_width","",50,0.0,1.0);
  TH1D *hrate_3He_xbj_SBS_fixed_width = new TH1D("hrate_3He_xbj_SBS_fixed_width","",50,0.0,1.0);

  // TH1D *hdA1n_500h_BB = new TH1D("hdA1n_500h_BB", "", nbins, &(bins[0]) );
  // TH1D *hdA1n_500h_SBS = new TH1D("hdA1n_500h_SBS", "", nbins, &(bins[0]) );
  
  TProfile *hQ2_xbj_SBS_prof = new TProfile("hQ2_xbj_SBS_prof","",nbins,&(bins[0]));
  TProfile *hEprime_xbj_SBS_prof = new TProfile("hEprime_xbj_SBS_prof","",nbins,&(bins[0]));
  TProfile *hxmean_xbj_SBS_prof = new TProfile("hxmean_xbj_SBS_prof","",nbins,&(bins[0]));
  TProfile *hW_xbj_SBS_prof = new TProfile("hW_xbj_SBS_prof","",nbins,&(bins[0]));
  TProfile *hy_xbj_SBS_prof = new TProfile("hy_xbj_SBS_prof","",nbins,&(bins[0]));
  TProfile *heps_xbj_SBS_prof = new TProfile("heps_xbj_SBS_prof","",nbins,&(bins[0]));
  TH2D *hQ2_xbj_SBS = new TH2D("hQ2_xbj_SBS","",100,0,1,100,1,12);
  TH2D *hEprime_xbj_SBS = new TH2D("hEprime_xbj_SBS","",100,0,1,100,0,11);
  
  
  g4sbs_a1n_tree *T = new g4sbs_a1n_tree(C);

  long nevent = 0;

  while( T->GetEntry(nevent++) ){
    if( nevent%1000 == 0 ) cout << nevent << endl;

    //sigma is in mb, we want to express in cm^2. 1 barn = 1e-24 cm^2, 1 mb = 1e-27 cm^2
    //double sigma = T->primaries_Sigma*1e-27;
    //figure out whether this is a proton file or a neutron file.

    //double sigma = T->ev_sigma; 
    double sigma = T->ev_sigma; //In the g4sbs DIS generator, sigma is given in cm^2/sr (per-nucleon)
    double genvol = T->ev_solang; //In the g4sbs DIS generator, ev.solang = genvol/fNevt
    //double rate = T->ev_rate; //In the g4sbs DIS generator, ev.rate = 
    
    double weight_p = 0.0; //This will be normalized to an event rate:
    double weight_n = 0.0;
    double weight_total = 0.0;

    int nucl = T->ev_nucl; //1 = proton, 0 = neutron
    
    TString fname = C->GetFile()->GetName();

    //cout << "Current filename = " << fname << endl;
    
    // if( proton_files.find( fname ) != proton_files.end() ){
    //   weight_p = 2.0 * sigma * Lumi / double(ngen_p);
    // } else if( neutron_files.find( fname ) != neutron_files.end() ){
    //   weight_n = sigma * Lumi / double(ngen_n);  
    // }

    
    
    if ( files.find( fname ) != files.end() ){
      // sigma is the per-nucleon cross section
      // genvol is the phase space volume divided by the number of "generated" events
      // Needs to be corrected for the simulation efficiency.
      // Lumi is the luminosity in e*atoms/cm^2/s  
      //weight_total = sigma * genvol * double(ngen_file[fname]) * Lumi / double( ntries );
      //In the SBS DIS generator, a neutron (proton) is chosen with probability 1/3 (2/3) for Helium-3.
      //Therefore, we just need to multiply the luminosity by the number of nucleons/atom to achieve sigma_3He = 2sigma_p + sigma_n
      if( nucl == 1 ){ //proton
	weight_p = sigma * genvol * double(ngen_file[fname]) * 3. * Lumi / double( ntries );
      } else if ( nucl == 0 ){ //neutron
	weight_n = sigma * genvol * double(ngen_file[fname]) * 3. * Lumi / double( ntries );
      }
    }
    
    //cout << "event " << nevent << ", weight_p, weight_n = " << weight_p << ", " << weight_n << endl;
    
    //Select events in which the primary electron made a track in SBS:
    if( T->Harm_SBSGEM_Track_ntracks > 0 ){
      for( int track=0; track<T->Harm_SBSGEM_Track_ntracks; track++ ){
	if( (*(T->Harm_SBSGEM_Track_PID))[track] == 11 &&
	    (*(T->Harm_SBSGEM_Track_MID))[track] == 0 &&
	    T->ev_ep >= Emin_SBS && T->ev_Q2 > 1.0 &&
	    T->ev_W2 > 4.0 ){ //Then this is the primary electron track:
	  double xbj = T->ev_xbj;
	  double Q2 = T->ev_Q2;
	  //double y = T->evy;
	  double Eprime = T->ev_ep;
	  double W = sqrt(T->ev_W2);
	  double etheta = T->ev_th;
	  double y = 1.0 - Eprime/T->gen_Ebeam;

	  //What if, instead of the "true" xbj, we use the "reconstructed" value?
	  //xbj = Q2/(2.*Mp*(T->gen_Ebeam - Eprime));
	  
	  double gamma2 = pow(2.*Mp*xbj,2)/Q2; // Q^2 = 2M nu x --> 2Mx = Q^2/nu --> gamma^2 = Q^2 / nu^2 

	  double epsilon = 1.0/( 1. + 2*(1+1./gamma2)*pow(tan(etheta/2.),2) );
	  
	  hrate_p_xbj_SBS->Fill( xbj, weight_p );
	  hrate_n_xbj_SBS->Fill( xbj, weight_n );
	  hrate_3He_xbj_SBS->Fill( xbj, weight_p + weight_n ); //=2p + n;

	  hrate_p_xbj_SBS_fixed_width->Fill( xbj, weight_p );
	  hrate_n_xbj_SBS_fixed_width->Fill( xbj, weight_n );
	  hrate_3He_xbj_SBS_fixed_width->Fill( xbj, weight_p + weight_n ); //=2p + n;
	  
	  hQ2_xbj_SBS->Fill( xbj, Q2, weight_p + weight_n );
	  hEprime_xbj_SBS->Fill( xbj, Eprime, weight_p + weight_n );

	  hQ2_xbj_SBS_prof->Fill( xbj, Q2, weight_p + weight_n );
	  hEprime_xbj_SBS_prof->Fill( xbj, Eprime, weight_p + weight_n );
	  hxmean_xbj_SBS_prof->Fill( xbj, xbj, weight_p + weight_n );
	  hW_xbj_SBS_prof->Fill( xbj, W, weight_p + weight_n );
	  hy_xbj_SBS_prof->Fill( xbj, y, weight_p + weight_n );
	  heps_xbj_SBS_prof->Fill( xbj, epsilon, weight_p + weight_n );
	}
      }
    }

    //Select events in which the primary electron made a track in SBS:
    if( T->Earm_BBGEM_Track_ntracks > 0 ){
      for( int track=0; track<T->Earm_BBGEM_Track_ntracks; track++ ){
	if( (*(T->Earm_BBGEM_Track_PID))[track] == 11 &&
	    (*(T->Earm_BBGEM_Track_MID))[track] == 0 &&
	    T->ev_ep >= Emin_BB && T->ev_Q2 > 1.0 &&
	    T->ev_W2 > 4.0 ){ //Then this is with a high degree of certainty the primary scattered electron track:
	  double xbj = T->ev_xbj;
	  double Q2 = T->ev_Q2;
	  //double y = T->evy;
	  double Eprime = T->ev_ep;
	  double W = sqrt(T->ev_W2);
	  double etheta = T->ev_th;
	  double y = 1.0 - Eprime/T->gen_Ebeam;

	  //xbj = Q2/(2.*Mp*(T->gen_Ebeam - Eprime));
	  
	  double gamma2 = pow(2.*Mp*xbj,2)/Q2; // Q^2 = 2M nu x --> 2Mx = Q^2/nu --> gamma^2 = Q^2 / nu^2 

	  double epsilon = 1.0/( 1. + 2*(1+1./gamma2)*pow(tan(etheta/2.),2) );
	  
	  hrate_p_xbj_BB->Fill( xbj, weight_p );
	  hrate_n_xbj_BB->Fill( xbj, weight_n );
	  hrate_3He_xbj_BB->Fill( xbj, weight_p + weight_n ); //=2p + n;

	  hrate_p_xbj_BB_fixed_width->Fill( xbj, weight_p );
	  hrate_n_xbj_BB_fixed_width->Fill( xbj, weight_n );
	  hrate_3He_xbj_BB_fixed_width->Fill( xbj, weight_p + weight_n ); //=2p + n;
	  
	  hQ2_xbj_BB->Fill( xbj, Q2, weight_p + weight_n );
	  hEprime_xbj_BB->Fill( xbj, Eprime, weight_p + weight_n );

	  hQ2_xbj_BB_prof->Fill( xbj, Q2, weight_p + weight_n );
	  hEprime_xbj_BB_prof->Fill( xbj, Eprime, weight_p + weight_n );
	  hxmean_xbj_BB_prof->Fill( xbj, xbj, weight_p + weight_n );
	  hW_xbj_BB_prof->Fill( xbj, W, weight_p + weight_n );
	  hy_xbj_BB_prof->Fill( xbj, y, weight_p + weight_n );
	  heps_xbj_BB_prof->Fill( xbj, epsilon, weight_p + weight_n );
	}
      }
    }
  }

  fout->Write();
  
}
