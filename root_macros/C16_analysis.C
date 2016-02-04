#include "C16_tree.C"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
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
#include "TProfile.h"

//#include "TIter.h"

void C16_analysis( const char *infilename, const char *outputfilename ){

  TFile *fout = new TFile(outputfilename,"RECREATE");
  
  TChain *C = new TChain("T");
  C->Add(infilename);

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

  C16_tree *T = new C16_tree(C);
  
  fout->cd();
  
  long nevent=0;

  TH1D *hNphesum_4x4_all = new TH1D("hNphesum_4x4_all","",250,0.0,2500.0);
  TH1D *hNphesum_4x4_cut = new TH1D("hNphesum_4x4_cut","",250,0.0,2500.0);

  TH2D *hNphe_vs_sumedep = new TH2D("hNphe_vs_sumedep","",50,0.0,2.0,50,0.0,2000.0);

  TH2D *Dose_rate_vs_row_col = new TH2D("Dose_rate_vs_row_col","",4,0.5,4.5,4,0.5,4.5);
  TH1D *Dose_rate_vs_Zdepth = new TH1D("Dose_rate_vs_Zdepth","",25,0.0,50.0);
  TClonesArray *histos_Dose_rate = new TClonesArray("TH1D",16);
  TClonesArray *histos_Dose_rate_vs_plane = new TClonesArray("TH1D",16);
  // double total_dose_rate_row_col_plane[4][4][10];
  // double total_dose_rate_row_col[4][4];

  TClonesArray *histos_nphe_vs_edep = new TClonesArray("TProfile",16);
  
  int cell=0;
  for( int row=1; row<=4; row++ ){
    for( int col=1; col<=4; col++ ){
      TString histname;
      histname.Form( "Dose_rate_vs_Z_row%d_col%d", row, col );
      new( (*histos_Dose_rate)[cell] ) TH1D( histname.Data(), "", 25,0.0,50.0 );
      for( int plane=1; plane<=10; plane++ ){
	// total_dose_rate_row_col_plane[row-1][col-1][plane-1] = 0.0;
	// total_dose_rate_row_col[row-1][col-1] = 0.0;
      }
      histname.Form( "Dose_rate_vs_plane_row%d_col%d", row, col );
      new( (*histos_Dose_rate_vs_plane)[cell] ) TH1D( histname.Data(), "", 10, 0.5,10.5 );

      histname.Form( "nphe_vs_edep_row%d_col%d", row, col );
      new( (*histos_nphe_vs_edep)[cell] ) TProfile( histname.Data(), "", 50, 0.0, 1.5);
      
      cell++;
    }
  }
  
  //   double zoffset = -.2479;  

  double Ibeam = 20.0e-6; //In amperes.
  double e = 1.602e-19; //electron charge
  double rho = 3.86; //g/cm^3
  double segthick = 3.4; //cm
  double segvolume = 4.2*4.2*segthick;
  double blockvolume = 4.2*4.2*34.3; 
  double lastsegvolume = 4.2*4.2*3.7;
  //We want dose rate in krad/h; 1 rad = 0.01 J/kg
  double weight = Ibeam/e/double(ngen); // Beam current in electrons/second divided by N generated events: gives rate in Hz!

  double mass_block = blockvolume*rho / 1000.0; //kg
  double mass_seg = segvolume*rho / 1000.0; //kg
  double mass_lastseg = lastsegvolume*rho / 1000.0; //kg

  double mass_total = mass_block * 16.0;
  double mass_zbin = 4.2*4.2*2.0*rho/1000.0; //kg
  
  while( T->GetEntry( nevent++ ) ){
    TFile *f = ( (TChain*) (T->fChain) )->GetFile();

    TString fname = f->GetName();
    
    if( bad_file_list.find( fname ) == bad_file_list.end() ){ 
    
      if( nevent%1000 == 0 ) {
	cout << nevent << endl;
	cout << "Current file = " << f->GetName() << endl;
      }
    
      double sum_nphe = 0.0;
    
      double nphe_max = 0.0;
      int rowmax = -1, colmax = -1;
      for( int hit=0; hit<T->Earm_C16_hit_nhits; hit++ ){
	sum_nphe += (*(T->Earm_C16_hit_NumPhotoelectrons))[hit];
	double sumedep_same_cell = 0.0;
	for( int jhit=0; jhit<T->Earm_C16TF1_hit_nhits; jhit++ ){
	  if( (*(T->Earm_C16TF1_hit_row))[jhit] == (*(T->Earm_C16_hit_row))[hit] &&
	      (*(T->Earm_C16TF1_hit_col))[jhit] == (*(T->Earm_C16_hit_col))[hit] ){
	    sumedep_same_cell += (*(T->Earm_C16TF1_hit_sumedep))[jhit];	    
	  }
	}
	hNphe_vs_sumedep->Fill( sumedep_same_cell, (*(T->Earm_C16_hit_NumPhotoelectrons))[hit] );
	double nphe = double( (*(T->Earm_C16_hit_NumPhotoelectrons))[hit] );
	if( nphe > nphe_max ){
	  nphe_max = nphe;
	  rowmax = (*(T->Earm_C16_hit_row))[hit] + 1;
	  colmax = (*(T->Earm_C16_hit_col))[hit] + 1;
	}

	int row = (*(T->Earm_C16_hit_row))[hit];
	int col = (*(T->Earm_C16_hit_col))[hit];
	int plane = (*(T->Earm_C16_hit_plane))[hit];
	int cell = col + 4*row;

	( (TH2D*) (*histos_nphe_vs_edep)[cell] )->Fill( sumedep_same_cell, (*(T->Earm_C16_hit_NumPhotoelectrons))[hit] );
	
      }
    
      hNphesum_4x4_all->Fill( sum_nphe );
      if( rowmax >= 2 && rowmax <= 3 &&
	  colmax >= 2 && colmax <= 3 ){ //ignore clusters with max at edge
	hNphesum_4x4_cut->Fill( sum_nphe );
      }
      double Dbb = T->gen_dbb*100.0; //convert to cm
      double thbb = T->gen_thbb;
      //Loop over TF1 hits:
      for( int hit=0; hit<T->Earm_C16TF1_hit_nhits; hit++ ){
	int row = (*(T->Earm_C16TF1_hit_row))[hit];
	int col = (*(T->Earm_C16TF1_hit_col))[hit];
	int plane = (*(T->Earm_C16TF1_hit_plane))[hit];
	int cell = col + 4*row;
	double edep = (*(T->Earm_C16TF1_hit_sumedep))[hit];
	double xhit = (*(T->Earm_C16TF1_hit_xhit))[hit]*100.0; //convert to cm
	double yhit = (*(T->Earm_C16TF1_hit_yhit))[hit]*100.0; //convert to cm
	double zhit = (*(T->Earm_C16TF1_hit_zhit))[hit]*100.0; //convert to cm
	
	double Zdepth = xhit * sin(thbb) + zhit * cos(thbb) - Dbb;
	
	//Compute the contribution of this hit to the dose rate:

	//edep is given in GeV:
	//weight is in Hz
	//mass_block is given in kg
	// 1 rad = 0.01 J/kg; 1 J/kg = 100 rad
	//Need to convert edep to J: 1 eV = 1.602e-19 J; 1 GeV = 1.602e-10 J
	double doserate = edep * e * 1.e9 / mass_block  * weight * 3600.0 * 100.0/1000.0; //now this is in J/(kg*s); now convert to krad/h 	
	double doserate_zbin = edep * e * 1.e9 / mass_zbin * weight * 3600.0 * 100.0/1000.0;
	
	//total_dose_rate_row_col_plane[row+1][col+1][plane+1] += doserate;
	
	Dose_rate_vs_Zdepth->Fill( Zdepth, doserate_zbin / 16.0 );
	Dose_rate_vs_row_col->Fill( col+1, row+1, doserate );
	( (TH1D*) (*histos_Dose_rate)[cell] )->Fill( Zdepth, doserate_zbin );
	if( plane+1 < 10 ){
	  ( (TH1D*) (*histos_Dose_rate_vs_plane)[cell] )->Fill( plane+1, doserate * mass_block / mass_seg );
	} else {
	  ( (TH1D*) (*histos_Dose_rate_vs_plane)[cell] )->Fill( plane+1, doserate * mass_block / mass_lastseg );
	}
      }
    }
  }

  cout << "Finished event loop " << endl;

  // histos_Dose_rate->Compress();
  // histos_Dose_rate_vs_plane->Compress();

  
  
  //fout->cd();
  fout->Write();
  //fout->Close();

}
