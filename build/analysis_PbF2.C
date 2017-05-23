#include "TMatrixD.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
//#include <pair>
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
#include "TTree.h"

void analysis_PbF2(const char* rootfile)
{
  TFile* f = new TFile(rootfile);
  
  TChain *C1 = (TChain*)f->Get("T");
  
  g4sbs_PbF2_tree *T1 = new g4sbs_PbF2_tree(C1);
  
  Long64_t NEvts1 = C1->GetEntries();
    
  cout << NEvts1 << " entries" << endl;

 
  TFile *fout=new TFile("kin_variables_PbF2.root","recreate");////////////////
  TTree *treeout=new TTree("T","Tree wf-analyzed ");/////////////////////


Double_t ener[208];
Double_t Edep;
Double_t Edep_Max;
Double_t ebeam=6.0;
Double_t M=0.93827;
Double_t Mpi=0.1349;
Double_t Me=0.000511;
//TLorentzVector p(0,0,0,M);
//TLorentzVector k(0,0,ebeam,ebeam);
Int_t row_edepmaxcal, col_edepmaxcal;
Int_t row_edepcal, col_edepcal;

Int_t evt; treeout->Branch("evt", &evt, "evt/I");
Double_t thbb; treeout->Branch("thbb", &thbb, "thbb/D");//Central angle of electron arm (BigBite or ECAL) in radians
Double_t xbj; treeout->Branch("xbj", &xbj, "xbj/D");
Double_t W2; treeout->Branch("W2", &W2, "W2/D");
Double_t Q2; treeout->Branch("Q2", &Q2, "Q2/D");
//Polar angle of scattered electron (radians) 
Double_t th; treeout->Branch("th", &th, "th/D");
//azimuthal angle of scattered electron (radians)
Double_t phi; treeout->Branch("phi", &phi, "phi/D");//angle between hadron production plane and lepton scattering plane in radians.
Double_t vv; treeout->Branch("vv", &vv, "vv/D");
//component of final-state electron momentum in GeV. 
Double_t k_primx; treeout->Branch("k_primx", &k_primx, "k_primx/D");
Double_t k_primy; treeout->Branch("k_primy", &k_primy, "k_primy/D");
Double_t k_primz; treeout->Branch("k_primz", &k_primz, "k_primz/D");
Double_t p_primx; treeout->Branch("p_primx", &p_primx, "p_primx/D");
Double_t p_primy; treeout->Branch("p_primy", &p_primy, "p_primy/D");
Double_t p_primz; treeout->Branch("p_primz", &p_primz, "p_primz/D");
//component of final-state nucleon momentum in GeV.

treeout->Branch("Edep", &Edep, "Edep/D");
treeout->Branch("row_edepcal", &row_edepcal, "row_edepcal/I");
treeout->Branch("col_edepcal", &col_edepcal, "col_edepcal/I");
treeout->Branch("ener", &ener, "ener[208]/D");

TH2D* h1_thetaVsphi_e =  new TH2D("h1_thetaVsphi_e", "", 360, -180, 180, 180, 0, 180);


  for(Long64_t nevent = 0; nevent<NEvts1; nevent++){
    if( nevent%10 == 0 ){cout << nevent << endl;}
    
    T1->GetEntry(nevent);

 h1_thetaVsphi_e->Fill(T1->ev_ph, T1->ev_th);
//////////////////////////////////////////////////////////////////////////////
   
  /*  cout << "======== number of hits in HCal: " << T1->Earm_BBSH_hit_nhits << endl;
    // loop on BBSHTF1 hits
    for(Int_t i = 0; i<T1->Earm_BBSHTF1_hit_nhits; i++){
      cout << "BBSHTF1 row: " << T1->Earm_BBSHTF1_hit_row->at(i) << ", column: " << T1->Earm_BBSHTF1_hit_col->at(i) << endl;
      cout << "energy deposited (GeV) in hit " << i << ": "<< T1->Earm_BBSHTF1_hit_sumedep->at(i) << endl;




    }*/

//00000000000000000000000000000000000000000000000000000000000000000000000000000000000
   // cout << "======= number of hits in PbF2: " << T1->BBSHPbF2_hit_nhits << endl;

 // loop on BBSHTF1 hits
 Edep=0.;
 for(Int_t j=0;j<208;j++){ener[j]=0.;}

    for(Int_t i = 0; i<T1->BBSHPbF2_hit_nhits; i++){
 //cout << "PbF2 row: " << T1->BBSHPbF2_hit_row->at(i) << ",column: " << T1->BBSHPbF2_hit_col->at(i) << endl;
 //cout << "energy deposited (GeV) in hit_PbF2 " << i << ": "<< T1->BBSHPbF2_hit_sumedep->at(i) << endl;


	row_edepcal = T1->BBSHPbF2_hit_row->at(i);
	col_edepcal = T1->BBSHPbF2_hit_col->at(i);

ener[T1->BBSHPbF2_hit_cell->at(i)]=T1->BBSHPbF2_hit_sumedep->at(i);

Edep +=T1->BBSHPbF2_hit_sumedep->at(i);

}

//00000000000000000000000000000000000000000000000000000000000000000000000000000000000
 for(Int_t i = 0; i<T1->BBcal_hit_nhits; i++){


 //Edep = T1->BBcal_hit_NumPhotoelectrons->at(i);

 /*Edep_Max=-10;
      if(Edep>Edep_Max){
	Edep_Max = Edep;
	row_edepmaxcal = T1->BBcal_hit_row->at(i);
	col_edepmaxcal = T1->BBcal_hit_col->at(i);
      }
*/

    }
   
    // loop on BBSH hits
    for(Int_t i = 0; i<T1->Earm_BBSH_hit_nhits; i++){
 //cout << "BBSH row: " << T1->Earm_BBSH_hit_row->at(i) << ", column: " << T1->Earm_BBSH_hit_col->at(i) << endl;
  //    cout << "number of photoelectrons in hit " << i << ": " << T1->Earm_BBSH_hit_NumPhotoelectrons->at(i) << endl;
      //T1->Earm_BBSH_hit_NumPhotoelectrons->at(i): 
      //procedure to access element i of your std::vector object
    }

//////////////////////////////////variables communes////////////////////////////////////:

            evt=T1->ev_count;
            thbb=T1->gen_thbb;
            xbj=T1->ev_xbj;
            W2=T1->ev_W2;
            Q2=T1->ev_Q2;
            th=T1->ev_th*180./TMath::Pi();
             phi=T1->ev_ph*180./TMath::Pi();
          //  phi=T1->ev_phih;
            vv=T1->ev_vz;
            k_primx=T1->ev_epx;
            k_primy=T1->ev_epy;
            k_primz=T1->ev_epz;

            p_primx=T1->ev_npx;
            p_primx=T1->ev_npy;
            p_primx=T1->ev_npz;

//cout << "epx " << k_primx << "  "<<T1->ev_epx<<endl;
    treeout->Fill();////////////////

  }//end evts
   TCanvas *c1 = new TCanvas("c1","c1",800,1000);
   c1->cd();

   h1_thetaVsphi_e->Draw();
 
 h1_thetaVsphi_e->Write();

 treeout->Write();////////////////



  
}//end macro

