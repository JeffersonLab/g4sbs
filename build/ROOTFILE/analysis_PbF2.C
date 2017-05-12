#include "TMatrixD.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <pair>
#include "TString.h"
#include "TDecompSVD.h"
#include "TCut.h"
#include "TEventList.h"
#include "g4sbs_PbF2_tree.C"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TF2.h"
 

void analysis_PbF2(const char* rootfile) //000
{
  TFile* f = new TFile(rootfile);//000
  
  TChain *C1 = (TChain*)f->Get("T");
  
  g4sbs_PbF2_tree *T1 = new g4sbs_PbF2_tree(C1);
  
  Long64_t NEvts1 = C1->GetEntries();
    
  cout << NEvts1 << " entries" << endl;
  
  for(Long64_t nevent = 0; nevent<NEvts1; nevent++){
    if( nevent%10 == 0 ){//000
	cout << nevent << endl;
    }
    
    T1->GetEntry(nevent);
    
    cout << "======== number of hits in BBSH: " << T1->Earm_BBSH_hit_nhits << endl;
    // loop on BBSHTF1 hits
    for(Int_t i = 0; i<T1->Earm_BBSHTF1_hit_nhits; i++){
      cout << "BBSHTF1 row: " << T1->Earm_BBSHTF1_hit_row->at(i) << ", column: " << T1->Earm_BBSHTF1_hit_col->at(i) << endl;
      cout << "energy deposited (GeV) in hit " << i << ": "<< T1->Earm_BBSHTF1_hit_sumedep->at(i) << endl;
    }

//00000000000000000000000000000000000000000000000000000000000000000000000000000000000
    cout << "======= number of hits in PbF2: " << T1->Harm_BBSHPbF2_hit_nhits << endl;
    // loop on BBSHTF1 hits
    for(Int_t i = 0; i<T1->Harm_BBSHTF1_hit_nhits; i++){
      cout << "PbF2 row: " << T1->Harm_BBSHPbF2_hit_row->at(i) << ", column: " << T1->Harm_BBSHPbF2_hit_col->at(i) << endl;
      cout << "energy deposited (GeV) in hit_PbF2 " << i << ": "<< T1->Harm_BBSHPbF2_hit_sumedep->at(i) << endl;
    }
//00000000000000000000000000000000000000000000000000000000000000000000000000000000000

   
    // loop on BBSH hits
    for(Int_t i = 0; i<T1->Earm_BBSH_hit_nhits; i++){
      cout << "BBSH row: " << T1->Earm_BBSH_hit_row->at(i) << ", column: " << T1->Earm_BBSH_hit_col->at(i) << endl;
      cout << "number of photoelectrons in hit " << i << ": " << T1->Earm_BBSH_hit_NumPhotoelectrons->at(i) << endl;
      //T1->Earm_BBSH_hit_NumPhotoelectrons->at(i): 
      //procedure to access element i of your std::vector object
    }
    
  }
  
}

