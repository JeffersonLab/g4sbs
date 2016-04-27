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

void C16_analysis( const char *infilename, const char *outputfilename ){
  
  TChain *C = new TChain("T");
  C->Add(infilename);

  TFile *fout = new TFile(outputfilename,"RECREATE");

  C16_tree *T = new C16_tree(C);

  long nevent=0;

  TH1D *hNphesum_4x4_all = new TH1D("hNphesum_4x4_all","",250,0.0,2500.0);
  TH1D *hNphesum_4x4_cut = new TH1D("hNphesum_4x4_cut","",250,0.0,2500.0);

  TH2D *hNphe_vs_sumedep = new TH2D("hNphe_vs_sumedep","",250,0.0,2.0,250,0.0,2500.0);
  
  while( T->GetEntry( nevent++ ) ){
    if( nevent%1000 == 0 ) cout << nevent << endl;

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
    }
    
    hNphesum_4x4_all->Fill( sum_nphe );
    if( rowmax >= 2 && rowmax <= 3 &&
	colmax >= 2 && colmax <= 3 ){ //ignore clusters with max at edge
      hNphesum_4x4_cut->Fill( sum_nphe );
    }
  }

  fout->cd();
  fout->Write();
  fout->Close();

}
