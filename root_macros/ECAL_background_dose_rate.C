#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include <vector>
#include <set>
#include <map>
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "G4SBSRunData.hh"
#include "TObjArray.h"
#include "TChainElement.h"
#include "TClonesArray.h"
#include "gep_tree_July2015.C"

void ECAL_background_dose_rate(const char *rootfilename, const char *outfilename){

  G4SBSRunData *rd;
  
  TChain *C = new TChain("T");
  C->Add(rootfilename);

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
  
  gep_tree_July2015 *T = new gep_tree_July2015( C );

  TFile *fout = new TFile(outfilename,"RECREATE");
  
  //For ECAL, we want to compute the total dose rate vs x/y at the face of the calorimeter,
  // Averaged over the face of ECAL as a function of depth;
  TH1D *hdose_rate_vs_depth_ECAL = new TH1D("hdose_rate_vs_depth_ECAL","",100,0.0,50.0);
  //TH1D *hdose_rate_vs_depth_ECAL_local = new TH1D("hdose_rate_vs_depth_ECAL_local","",100,0.0,50);
  TH2D *hdose_rate_vs_xy_ECAL = new TH2D("hdose_rate_vs_xy_ECAL","",100,-80.0,80.0,100,-200.0,200.0);

  TH2D *hdose_rate_vs_row_col_ECAL = new TH2D("hdose_rate_vs_row_col_ECAL","",32,0.5,32.5,85,0.5,85.5);
  
  double Ibeam = 75.0e-6; //Amps
  double e = 1.602e-19; //C
  double rho = 3.86; //g/cm^3;
  double thick = 40.0; //cm

  //We want to compute the total volume/mass of ECAL to normalize to a dose rate:
  double area = 999.0 * 4.2 * 4.2 + 697 * 4.0 * 4.0 + 81 * 3.8 * 3.8; //cm^2
  double volume = area * thick; //cm^3
  double totalmass = rho * volume / 1000.0; //kg
  double avgblockmass = totalmass / (999.0 + 697.0 + 81.0);

  double weight = Ibeam/e/double(ngen);

  double binwidth_Z = 0.5; //cm
  double volume_zbin = area * binwidth_Z;
  double mass_zbin = volume_zbin * rho / 1000.0; //kg
  
  long nevent=0;

  while( T->GetEntry( nevent++ ) ){
    TFile *f = ( (TChain*) (T->fChain) )->GetFile();

    double Dbb = T->gen_dbb*100.0 - thick;
    double thbb = T->gen_thbb;
    
    if( bad_file_list.find( f->GetName() ) == bad_file_list.end() ){
      if( nevent%1000 == 0 ) cout << nevent << endl;

      for( int hit=0; hit<T->Earm_ECalTF1_hit_nhits; hit++ ){
	int row = (*(T->Earm_ECalTF1_hit_row))[hit];
	int col = (*(T->Earm_ECalTF1_hit_col))[hit];
	double edep = (*(T->Earm_ECalTF1_hit_sumedep))[hit]; //Is in GeV
	double xhit = (*(T->Earm_ECalTF1_hit_xhit))[hit]*100.0;
	double yhit = (*(T->Earm_ECalTF1_hit_yhit))[hit]*100.0;
	double zhit = (*(T->Earm_ECalTF1_hit_zhit))[hit]*100.0;

	double Zdepth = xhit * sin(thbb) + zhit * cos(thbb) - Dbb;
	double xcalo = xhit * cos(thbb) - zhit * sin(thbb);
	double ycalo = yhit;

	double doserate = edep * e * 1.e9 / totalmass * weight * 3600.0 * 100.0/1000.0; //should be in krad/h

	double doserate_zbin = edep * e * 1.e9 / mass_zbin * weight * 3600.0 * 0.1;
	
	hdose_rate_vs_depth_ECAL->Fill( Zdepth, doserate_zbin );
	//hdose_rate_vs_depth_ECAL_local->Fill( Zdepth, doserate * binwidth_Z/thick );
	hdose_rate_vs_xy_ECAL->Fill( xcalo, ycalo, doserate );

	hdose_rate_vs_row_col_ECAL->Fill( col, row, doserate * totalmass/avgblockmass );
	
      }
      
    }
  }

  fout->Write();
  fout->Close();
}
