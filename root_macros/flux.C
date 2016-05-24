#include "TChain.h"
#include "TTree.h"
#include "G4SBSRunData.hh"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObjArray.h"
#include "TChainElement.h"
#include "TMath.h"
#include "C16_tree_with_flux.C"
#include <set>
#include <iostream>


void flux(const char *rootfilename, const char *outputfilename){

  gStyle->SetMarkerStyle(20);
  
  TH1::SetDefaultSumw2();
  
  TChain *C = new TChain("T");
  C->Add(rootfilename);

  G4SBSRunData *rd;

  long ngen = 0;
  int nfiles = 0;
  
  TObjArray *FileList = C->GetListOfFiles();
  TIter next(FileList);

  TChainElement *chEl = 0;

  set<TString> bad_file_list;
  
  while( (chEl=(TChainElement*)next() )){
    TFile newfile(chEl->GetTitle());
    newfile.GetObject("run_data",rd);
    if( rd ){
      ngen += rd->fNtries;
      nfiles++;
    } else {
      bad_file_list.insert( chEl->GetTitle());
    }
    //cout << chEl->GetTitle() << endl;
  }

  cout << "number of generated events = " << ngen << endl;

  C16_tree_with_flux *T = new C16_tree_with_flux(C);
  
  TFile *fout = new TFile(outputfilename,"RECREATE");

  long nevent=0;

  double cutgamma = 0.05*1.e-3; //energy cutoff in GeV for photons
  double cutelectron = 5.e-3; //range cutoff in GeV for electrons:
  
  TH1D *hphoton_rate_vs_theta = new TH1D("hphoton_rate_vs_theta","",80,0.0,40.0);
  TH1D *helectron_rate_vs_theta = new TH1D("helectron_rate_vs_theta","",80,0.0,40.0);
  TH1D *hphoton_rate_vs_theta_Edep = new TH1D("hphoton_rate_vs_theta_Edep","",80,0.0,40.0);
  TH1D *helectron_rate_vs_theta_Edep = new TH1D("helectron_rate_vs_theta_Edep","",80,0.0,40.0);
  TH1D *hother_rate_vs_theta = new TH1D("hother_rate_vs_theta","",80,0.0,40.0);
  TH1D *hother_rate_vs_theta_Edep = new TH1D("hother_rate_vs_theta_Edep","",80,0.0,40.0);

  TH1D *htotal_rate_vs_theta = new TH1D("htotal_rate_vs_theta","",80,0.0,40.0);
  TH1D *htotal_rate_vs_theta_Edep = new TH1D("htotal_rate_vs_theta_Edep","",80,0.0,40.0);

  double PI = TMath::Pi();

  double binwidth_theta = 40.0/80.0*PI/180.0;

  double Ibeam = 20e-6;
  double e = 1.602e-19;
  
  while( T->GetEntry( nevent++ ) ){
    TFile *f = ( (TChain*) (T->fChain) )->GetFile();

    TString fname = f->GetName();
    
    if( bad_file_list.find( fname ) == bad_file_list.end() ){ 
    
      if( nevent%1000 == 0 ) {
	cout << nevent << endl;
	cout << "Current file = " << f->GetName() << endl;
      }
      
      for( int part=0; part<T->FLUX_npart_CAL; part++ ){
	double E = (*(T->FLUX_E))[part];
	double px = (*(T->FLUX_px))[part];
	double py = (*(T->FLUX_py))[part];
	double pz = (*(T->FLUX_pz))[part];
	double p = (*(T->FLUX_p))[part];
	double theta = acos( pz/p );
	double phi = atan2( py, px );

	int bintheta = hphoton_rate_vs_theta->FindBin( theta * 180.0/PI );

	double dcostheta = cos( PI/180.0*hphoton_rate_vs_theta->GetBinLowEdge(bintheta) ) -
	  cos( PI/180.0*hphoton_rate_vs_theta->GetBinLowEdge(bintheta+1) );
	double dOmega = 2.0*PI *dcostheta;

	double weight = Ibeam/e/double(ngen)/dOmega;
	
	int pid = (*(T->FLUX_pid))[part];

	htotal_rate_vs_theta->Fill( theta*180.0/PI, weight );
	htotal_rate_vs_theta_Edep->Fill( theta*180.0/PI, weight*E );
	
	if( pid == 22 && E >= cutgamma ){//photon
	  hphoton_rate_vs_theta->Fill( theta*180.0/PI, weight );
	  hphoton_rate_vs_theta_Edep->Fill( theta*180.0/PI, weight*E );
	} else if( fabs(pid) == 11 && E >= cutelectron ){ //e+/e-
	  helectron_rate_vs_theta->Fill( theta*180.0/PI, weight );
	  helectron_rate_vs_theta_Edep->Fill( theta*180.0/PI, weight*E );
	} else { //other
	  hother_rate_vs_theta->Fill( theta*180.0/PI, weight );
	  hother_rate_vs_theta_Edep->Fill( theta*180.0/PI, weight*E );
	}
      }

    }
  }

  hphoton_rate_vs_theta->GetXaxis()->SetTitle("#theta (#circ)");
  hphoton_rate_vs_theta->GetYaxis()->SetTitle("#frac{dN_{#gamma}}{dt d#Omega} (Photons/s/sr)");

  helectron_rate_vs_theta->GetXaxis()->SetTitle("#theta (#circ)");
  helectron_rate_vs_theta->GetYaxis()->SetTitle("#frac{dN_{e}}{dt d#Omega} (Electrons/s/sr)");

  hother_rate_vs_theta->GetXaxis()->SetTitle("#theta (#circ)");
  hother_rate_vs_theta->GetYaxis()->SetTitle("#frac{dN_{other}}{dt d#Omega} (Particles/s/sr)");

  htotal_rate_vs_theta->GetXaxis()->SetTitle("#theta (#circ)");
  htotal_rate_vs_theta->GetYaxis()->SetTitle("#frac{dN_{total}}{dt d#Omega} (Particles/s/sr)");

  hphoton_rate_vs_theta_Edep->GetXaxis()->SetTitle("#theta (#circ)");
  hphoton_rate_vs_theta_Edep->GetYaxis()->SetTitle("#frac{dE}{dt d#Omega} (GeV/s/sr)");

  helectron_rate_vs_theta_Edep->GetXaxis()->SetTitle("#theta (#circ)");
  helectron_rate_vs_theta_Edep->GetYaxis()->SetTitle("#frac{dE}{dt d#Omega} (GeV/s/sr)");

  hother_rate_vs_theta_Edep->GetXaxis()->SetTitle("#theta (#circ)");
  hother_rate_vs_theta_Edep->GetYaxis()->SetTitle("#frac{dE}{dt d#Omega} (GeV/s/sr)");

  htotal_rate_vs_theta_Edep->GetXaxis()->SetTitle("#theta (#circ)");
  htotal_rate_vs_theta_Edep->GetYaxis()->SetTitle("#frac{dE}{dt d#Omega} (GeV/s/sr)");
  
  fout->Write();
}
