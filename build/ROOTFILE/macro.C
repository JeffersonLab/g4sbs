#include "TCanvas.h"
 #include "TROOT.h"
 #include "TGraphErrors.h"
 #include "TF1.h"
 #include "TH1F.h"
 #include "TLegend.h"
 #include "TMath.h"
 #include "TLorentzVector.h"
 #include "TChain.h"
 #include "TTree.h"
 #include "TStyle.h"
 #include "TSystem.h"
 #include "TRandom.h"
 #include "TFile.h"
 #include "TMatrixD.h"
// #include "fonction.h"
 #include <iostream>
 #include <fstream>

void macro(){
TChain *tree=new TChain("T");

 tree->Add("/home/meriem/g4sbs/build/test.root");

Int_t nentries = tree->GetEntries(); 
  cout<<"nentries = "<<nentries<<endl;
//nentries = 10000;

Int_t Earm_BBGEM_hit_nhits;tree->SetBranchAddress("Earm_BBGEM_hit_nhits",&Earm_BBGEM_hit_nhits);
//Double_t xbj;tree->SetBranchAddress("xbj",&xbj);
/*
Int_t run;tree->SetBranchAddress("run",&run);
Double_t kp[3];tree->SetBranchAddress("kp",&kp);
Double_t q1[3];tree->SetBranchAddress("q1",&q1);
Double_t ene1;tree->SetBranchAddress("ene1",&ene1);
Double_t xc1;tree->SetBranchAddress("xc1",&xc1);
Double_t yc1;tree->SetBranchAddress("yc1",&yc1);
Double_t size1;tree->SetBranchAddress("size1",&size1);
Double_t mm2;tree->SetBranchAddress("mm2",&mm2);
Double_t dp;tree->SetBranchAddress("dp",&dp); // the va
Double_t the_s;tree->SetBranchAddress("the_s",&the_s);
Double_t phi_s;tree->SetBranchAddress("phi_s",&phi_s);
Double_t vv;tree->SetBranchAddress("vv",&vv); //vertex
*/
////////////////////////////////////////////
Double_t ebeam;
ebeam=2.2;
TFile *fout=new TFile("variables.root","recreate");
TTree *treeout=new TTree("T","Tree wf-analyzed ");


//treeout->Branch("run",&run,"run/I");
//treeout->Branch("t",&t,"t/D");
//treeout->Branch("Earm_BBGEM_hit_nhits",&Earm_BBGEM_hit_nhits,"Earm_BBGEM_hit_nhits/I");
//treeout->Branch("xbj",&xbj,"xbj/D");


 for(Int_t evt=0;evt<nentries;evt++){
//if(evt%10==0) cout<<evt<<"/"<<nentries<<" "<<run<<endl;
tree->GetEntry(evt);
//if(xc1>-21.5&&xc1<12.2&&yc1<21.5&&yc1>-21.4&&ene1>1.&&mm2>0.&&mm2<2.){





treeout->Fill();

//}//if



}//ev loop










}
