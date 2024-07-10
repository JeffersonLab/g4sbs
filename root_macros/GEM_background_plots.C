#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "G4SBSRunData.hh"
//#include "TIter.h"
#include "TChainElement.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TVector3.h"

#include <iostream>
#include <vector> 
#include <set>
#include <map>
#include "gep_ntuple.C"

using namespace std;

//New version for final polarimeter layout:

void GEM_background_plots( const char *rootfilename, const char *outputfilename ){

  //double nplanes[2] = {8,8};
  //double wplanes[3] = {40.0,60.0,60.0};
  //double hplanes[3] = {150.0,200.0,200.0};

  int nplanes_FT = 8;
  int nplanes_FPP = 8;

  double thresh_CDET = 0.005; //5 MeV CDET threshold
  
  vector<double> wplanes_FT(nplanes_FT,40.0);
  vector<double> hplanes_FT(nplanes_FT,150.0);

  wplanes_FT[6] = 60.0;
  wplanes_FT[7] = 60.0;
  hplanes_FT[6] = 200.0;
  hplanes_FT[7] = 200.0;

  vector<double> wplanes_FPP(nplanes_FPP,60.0);
  vector<double> hplanes_FPP(nplanes_FPP,200.0);
  
  
  TH1D::SetDefaultSumw2();

  TChain *C = new TChain("T");
  C->Add( rootfilename );
  
  TFile *fout = new TFile( outputfilename, "RECREATE" );

  G4SBSRunData *rundata;

  TObjArray *filelist = C->GetListOfFiles();

  TIter next(filelist);

  TChainElement *chEl = 0;

  gep_ntuple *T = new gep_ntuple(C);

  long nevent = 0;

  set<TString> bad_file_list;

  long nevents_generated = 0;
  while( (chEl=(TChainElement*)next() )){
    TFile newfile(chEl->GetTitle());
    newfile.GetObject("run_data",rundata);
    if( rundata ){
      nevents_generated += rundata->fNtries;
      cout << "Added file " << chEl->GetTitle() << endl;
    } else {
      bad_file_list.insert(chEl->GetTitle());
      cout << "File " << chEl->GetTitle() << " added to bad files list" << endl;
    }
  }

  cout << "Number of generated events = " << nevents_generated << endl;

  fout->cd();

  TH1D *hrate_layer_FT = new TH1D("hrate_layer_FT","SBS front tracker; GEM layer; Hit rate (Hz/cm^{2})",8,0.5,8.5);
  TH1D *hrate_layer_FPP1 = new TH1D("hrate_layer_FPP1","SBS back tracker; GEM layer; Hit rate (Hz/cm^{2})",8,0.5,8.5);

  TH1D *hrate_edep_layer_FT = new TH1D("hrate_edep_layer_FT","SBS front tracker primary ionization region; GEM layer; Energy deposit rate (MeV/cm^{2}/s)",8,0.5,8.5);
  TH1D *hrate_edep_layer_FPP1 = new TH1D("hrate_edep_layer_FPP1","SBS back tracker primary ionization region; GEM layer; Energy deposit rate (MeV/cm^{2}/s)",8,0.5,8.5);

  TH1D *hrate_versus_channel_CDET = new TH1D("hrate_versus_channel_CDET", "CDET (energy threshold = 5 MeV); channel ; Hit rate (Hz)",2401,-0.5,2400.5);
  
  //TH1D *hrate_layer_FPP2 = new TH1D("hrate_layer_FPP2","",5,0.5,5.5);
  TH2D *hvx_vz_layer1_FT = new TH2D("hvx_vz_layer1_FT","FT layer 1; vertex z (m); vertex x (m)",200,-1,7,200,-2,1);
  TH1D *hrate_vs_x_layer1_FT = new TH1D("hrate_vs_x_layer1_FT","FT layer 1; hit x (m); hit rate (Hz/cm^{2})",200,-.75,.75);
  TH1D *hrate_vs_y_layer1_FT = new TH1D("hrate_vs_y_layer1_FT","FT layer 1; hit y (m); hit rate (Hz/cm^{2})",200,-0.2,0.2);
  TH1D *hrate_vs_p_layer1_FT_wide = new TH1D("hrate_vs_p_layer1_FT_wide","FT layer 1; p (GeV); Hit rate (Hz/cm^{2})",200,0.0,8.0);
  TH1D *hrate_vs_p_layer1_FT_zoom = new TH1D("hrate_vs_p_layer1_FT_zoom","FT layer 1; p (GeV); Hit rate (Hz/cm^{2})",200,0.0,0.01);
  //TH1D *hrate_vs_p_electron_layer1_FT_wide = new TH1D("hrate_vs_p_electron_layer1_FT","",200,0.0,8.0);

  TH1D *hrate_vs_p_gamma_FT = new TH1D("hrate_vs_p_gamma_FT","FT #gamma-induced hits; p (GeV); rate (Hz/cm^{2}/layer)",200,0.0,0.0002);
  TH1D *hrate_vs_p_electron_FT = new TH1D("hrate_vs_p_electron_FT","FT e^{-}-induced hits; p (GeV); rate (Hz/cm^{2}/layer)",200,0.0,.1);
  TH1D *hrate_vs_p_pion_FT = new TH1D("hrate_vs_p_pion_FT","FT charged pion-induced hits; p (GeV); rate (Hz/cm^{2}/layer)",200,0.0,8.0);
  TH1D *hrate_vs_p_proton_FT = new TH1D("hrate_vs_p_proton_FT","FT proton-induced hits; p (GeV); rate (Hz/cm^{2}/layer)",200,0.0,8.0);
  TH1D *hrate_vs_p_other_FT = new TH1D("hrate_vs_p_other_FT","FT other-induced hits; p (GeV); rate (Hz/cm^{2}/layer)",200,0.0,8.0);
  
  TH1D *hrate_vs_p_gamma_FPP1 = new TH1D("hrate_vs_p_gamma_FPP1","FPP #gamma-induced hits; p (GeV); rate (Hz/cm^{2}/layer)",200,0.0,.0002);
  TH1D *hrate_vs_p_electron_FPP1 = new TH1D("hrate_vs_p_electron_FPP1","FPP e^{-}-induced hits; p (GeV); rate (Hz/cm^{2}/layer)",200,0.0,.1);
  TH1D *hrate_vs_p_pion_FPP1 = new TH1D("hrate_vs_p_pion_FPP1","FPP pion-induced hits; p (GeV); rate (Hz/cm^{2}/layer)",200,0.0,8.0);
  TH1D *hrate_vs_p_proton_FPP1 = new TH1D("hrate_vs_p_proton_FPP1","FPP proton-induced hits; p (GeV); rate (Hz/cm^{2}/layer)",200,0.0,8.0);
  TH1D *hrate_vs_p_other_FPP1 = new TH1D("hrate_vs_p_other_FPP1","FPP other-induced hits; p (GeV); rate (Hz/cm^{2}/layer)",200,0.0,8.0);

  
  
  // TH1D *hrate_vs_p_gamma_FPP2 = new TH1D("hrate_vs_p_gamma_FPP2","",200,0.0,.0002);
  // TH1D *hrate_vs_p_electron_FPP2 = new TH1D("hrate_vs_p_electron_FPP2","",200,0.0,.1);
  // TH1D *hrate_vs_p_pion_FPP2 = new TH1D("hrate_vs_p_pion_FPP2","",200,0.0,8.0);
  // TH1D *hrate_vs_p_proton_FPP2 = new TH1D("hrate_vs_p_proton_FPP2","",200,0.0,8.0);
  // TH1D *hrate_vs_p_other_FPP2 = new TH1D("hrate_vs_p_other_FPP2","",200,0.0,8.0);
  
  
  TH1D *hedep_FT = new TH1D("hedep_FT","Front Tracker; Energy deposit (GeV); Hit rate (Hz/cm^{2}/layer)",200,0.0,2e-5);
  TH1D *hedep_FPP1 = new TH1D("hedep_FPP1","Back Tracker; Energy deposit (GeV); Hit rate (Hz/cm^{2}/layer)",200,0.0,2e-5);
  //TH1D *hedep_FPP2 = new TH1D("hedep_FPP2","",200,0.0,2e-5);

  TH2D *hedep_vs_p_FT = new TH2D("hedep_vs_p_FT","",200,0.0,0.1,200,0.0,2e-5);
  TH2D *hedep_vs_p_FPP1 = new TH2D("hedep_vs_p_FPP1","",200,0.0,0.1,200,0.0,2e-5);
  //TH2D *hedep_vs_p_FPP2 = new TH2D("hedep_vs_p_FPP2","",200,0.0,0.1,200,0.0,2e-5);
  
  
  TH1D *hedep_FT_electron = new TH1D("hedep_FT_electron","",200,0.0,2e-5);
  TH1D *hedep_FT_other = new TH1D("hedep_FT_other","",200,0.0,2e-5);

  TH1D *hedep_FPP1_electron = new TH1D("hedep_FPP1_electron","",200,0.0,2e-5);
  TH1D *hedep_FPP1_other = new TH1D("hedep_FPP1_other","",200,0.0,2e-5);

  // TH1D *hedep_FPP2_electron = new TH1D("hedep_FPP2_electron","",200,0.0,2e-5);
  // TH1D *hedep_FPP2_other = new TH1D("hedep_FPP2_other","",200,0.0,2e-5);
  
  TH2D *hvx_vz_FT = new TH2D("hvx_vz_FT","SBS front tracker; Hit vertex z (m); Hit vertex x (m)",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FT_onebounce = new TH2D("hvx_vz_FT_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FT_gamma_onebounce = new TH2D("hvx_vz_FT_gamma_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FT_electron_onebounce = new TH2D("hvx_vz_FT_electron_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FT_pion_onebounce = new TH2D("hvx_vz_FT_pion_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FT_proton_onebounce = new TH2D("hvx_vz_FT_proton_onebounce","",200,-1,7,200,-2,1);
  // //TH2D *hvx_vz_FT_neutron_onebounce = new TH2D("hvx_vz_FT_neutron_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FT_other_onebounce = new TH2D("hvx_vz_FT_other_onebounce","",200,-1,7,200,-2,1);
  // TH1D *hp_onebounce_FT = new TH1D("hp_onebounce_FT","",200,0.0,8.0);

  // TH2D *hvx_vz_bounce1_FT = new TH2D("hvx_vz_bounce1_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm3_FT = new TH2D("hvx_vz_bounceNm3_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm2_FT = new TH2D("hvx_vz_bounceNm2_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm1_FT = new TH2D("hvx_vz_bounceNm1_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceN_FT = new TH2D("hvx_vz_bounceN_FT", "", 200, -1, 7, 200, -2, 1);
  
  // TH2D *hvx_vz_bounce1_gamma_FT = new TH2D("hvx_vz_bounce1_gamma_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm3_gamma_FT = new TH2D("hvx_vz_bounceNm3_gamma_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm2_gamma_FT = new TH2D("hvx_vz_bounceNm2_gamma_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm1_gamma_FT = new TH2D("hvx_vz_bounceNm1_gamma_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceN_gamma_FT = new TH2D("hvx_vz_bounceN_gamma_FT", "", 200, -1, 7, 200, -2, 1);
  
  // TH2D *hvx_vz_bounce1_electron_FT = new TH2D("hvx_vz_bounce1_electron_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm3_electron_FT = new TH2D("hvx_vz_bounceNm3_electron_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm2_electron_FT = new TH2D("hvx_vz_bounceNm2_electron_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm1_electron_FT = new TH2D("hvx_vz_bounceNm1_electron_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceN_electron_FT = new TH2D("hvx_vz_bounceN_electron_FT", "", 200, -1, 7, 200, -2, 1);

  // TH2D *hvx_vz_bounce1_other_FT = new TH2D("hvx_vz_bounce1_other_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm3_other_FT = new TH2D("hvx_vz_bounceNm3_other_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm2_other_FT = new TH2D("hvx_vz_bounceNm2_other_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm1_other_FT = new TH2D("hvx_vz_bounceNm1_other_FT", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceN_other_FT = new TH2D("hvx_vz_bounceN_other_FT", "", 200, -1, 7, 200, -2, 1);

  TH2D *hvx_vz_FPP1 = new TH2D("hvx_vz_FPP1","SBS back tracker; Hit vertex z (m); Hit vertex x(m)",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FPP1_onebounce = new TH2D("hvx_vz_FPP1_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FPP1_gamma_onebounce = new TH2D("hvx_vz_FPP1_gamma_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FPP1_electron_onebounce = new TH2D("hvx_vz_FPP1_electron_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FPP1_pion_onebounce = new TH2D("hvx_vz_FPP1_pion_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FPP1_proton_onebounce = new TH2D("hvx_vz_FPP1_proton_onebounce","",200,-1,7,200,-2,1);
  // //TH2D *hvx_vz_FPP1_neutron_onebounce = new TH2D("hvx_vz_FPP1_neutron_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FPP1_other_onebounce = new TH2D("hvx_vz_FPP1_other_onebounce","",200,-1,7,200,-2,1);
  // TH1D *hp_onebounce_FPP1 = new TH1D("hp_onebounce_FPP1","",200,0.0,8.0);

  // TH2D *hvx_vz_bounce1_FPP1 = new TH2D("hvx_vz_bounce1_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm3_FPP1 = new TH2D("hvx_vz_bounceNm3_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm2_FPP1 = new TH2D("hvx_vz_bounceNm2_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm1_FPP1 = new TH2D("hvx_vz_bounceNm1_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceN_FPP1 = new TH2D("hvx_vz_bounceN_FPP1", "", 200, -1, 7, 200, -2, 1);
  
  // TH2D *hvx_vz_bounce1_gamma_FPP1 = new TH2D("hvx_vz_bounce1_gamma_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm3_gamma_FPP1 = new TH2D("hvx_vz_bounceNm3_gamma_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm2_gamma_FPP1 = new TH2D("hvx_vz_bounceNm2_gamma_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm1_gamma_FPP1 = new TH2D("hvx_vz_bounceNm1_gamma_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceN_gamma_FPP1 = new TH2D("hvx_vz_bounceN_gamma_FPP1", "", 200, -1, 7, 200, -2, 1);
  
  // TH2D *hvx_vz_bounce1_electron_FPP1 = new TH2D("hvx_vz_bounce1_electron_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm3_electron_FPP1 = new TH2D("hvx_vz_bounceNm3_electron_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm2_electron_FPP1 = new TH2D("hvx_vz_bounceNm2_electron_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm1_electron_FPP1 = new TH2D("hvx_vz_bounceNm1_electron_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceN_electron_FPP1 = new TH2D("hvx_vz_bounceN_electron_FPP1", "", 200, -1, 7, 200, -2, 1);

  // TH2D *hvx_vz_bounce1_other_FPP1 = new TH2D("hvx_vz_bounce1_other_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm3_other_FPP1 = new TH2D("hvx_vz_bounceNm3_other_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm2_other_FPP1 = new TH2D("hvx_vz_bounceNm2_other_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm1_other_FPP1 = new TH2D("hvx_vz_bounceNm1_other_FPP1", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceN_other_FPP1 = new TH2D("hvx_vz_bounceN_other_FPP1", "", 200, -1, 7, 200, -2, 1);
  
  //TH2D *hvx_vz_FPP2 = new TH2D("hvx_vz_FPP2","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FPP2_onebounce = new TH2D("hvx_vz_FPP2_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FPP2_gamma_onebounce = new TH2D("hvx_vz_FPP2_gamma_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FPP2_electron_onebounce = new TH2D("hvx_vz_FPP2_electron_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FPP2_pion_onebounce = new TH2D("hvx_vz_FPP2_pion_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FPP2_proton_onebounce = new TH2D("hvx_vz_FPP2_proton_onebounce","",200,-1,7,200,-2,1);
  // //TH2D *hvx_vz_FPP2_neutron_onebounce = new TH2D("hvx_vz_FPP2_neutron_onebounce","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_FPP2_other_onebounce = new TH2D("hvx_vz_FPP2_other_onebounce","",200,-1,7,200,-2,1);
  // TH1D *hp_onebounce_FPP2 = new TH1D("hp_onebounce_FPP2","",200,0.0,8.0);

  // TH2D *hvx_vz_bounce1_FPP2 = new TH2D("hvx_vz_bounce1_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm3_FPP2 = new TH2D("hvx_vz_bounceNm3_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm2_FPP2 = new TH2D("hvx_vz_bounceNm2_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm1_FPP2 = new TH2D("hvx_vz_bounceNm1_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceN_FPP2 = new TH2D("hvx_vz_bounceN_FPP2", "", 200, -1, 7, 200, -2, 1);
  
  // TH2D *hvx_vz_bounce1_gamma_FPP2 = new TH2D("hvx_vz_bounce1_gamma_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm3_gamma_FPP2 = new TH2D("hvx_vz_bounceNm3_gamma_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm2_gamma_FPP2 = new TH2D("hvx_vz_bounceNm2_gamma_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm1_gamma_FPP2 = new TH2D("hvx_vz_bounceNm1_gamma_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceN_gamma_FPP2 = new TH2D("hvx_vz_bounceN_gamma_FPP2", "", 200, -1, 7, 200, -2, 1);
  
  // TH2D *hvx_vz_bounce1_electron_FPP2 = new TH2D("hvx_vz_bounce1_electron_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm3_electron_FPP2 = new TH2D("hvx_vz_bounceNm3_electron_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm2_electron_FPP2 = new TH2D("hvx_vz_bounceNm2_electron_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm1_electron_FPP2 = new TH2D("hvx_vz_bounceNm1_electron_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceN_electron_FPP2 = new TH2D("hvx_vz_bounceN_electron_FPP2", "", 200, -1, 7, 200, -2, 1);

  // TH2D *hvx_vz_bounce1_other_FPP2 = new TH2D("hvx_vz_bounce1_other_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm3_other_FPP2 = new TH2D("hvx_vz_bounceNm3_other_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm2_other_FPP2 = new TH2D("hvx_vz_bounceNm2_other_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceNm1_other_FPP2 = new TH2D("hvx_vz_bounceNm1_other_FPP2", "", 200, -1, 7, 200, -2, 1);
  // TH2D *hvx_vz_bounceN_other_FPP2 = new TH2D("hvx_vz_bounceN_other_FPP2", "", 200, -1, 7, 200, -2, 1);
  
  // TH1D *hnbounce_FT = new TH1D("hnbounce_FT","",16,-0.5,15.5);
  // TH1D *hnbounce_FPP1 = new TH1D("hnbounce_FPP1","",16,-0.5,15.5);
  // TH1D *hnbounce_FPP2 = new TH1D("hnbounce_FPP2","",16,-0.5,15.5);

  // TH2D *hvx_vz_mother_e_gamma_FT = new TH2D("hvx_vz_mother_e_gamma_FT","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_mother_e_e_FT = new TH2D("hvx_vz_mother_e_e_FT","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_mother_e_other_FT = new TH2D("hvx_vz_mother_e_other_FT","",200,-1,7,200,-2,1);
  // TH1D *hp_mother_e_gamma_FT = new TH1D("hp_mother_e_gamma_FT","",200,0.0,1.0);
  // TH1D *hp_mother_e_e_FT = new TH1D("hp_mother_e_e_FT","",200,0.0,2.0);
  // TH1D *hp_mother_e_other_FT = new TH1D("hp_mother_e_other_FT","",200,0.0,8.0);

  // TH2D *hvx_vz_mother_e_gamma_FPP1 = new TH2D("hvx_vz_mother_e_gamma_FPP1","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_mother_e_e_FPP1 = new TH2D("hvx_vz_mother_e_e_FPP1","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_mother_e_other_FPP1 = new TH2D("hvx_vz_mother_e_other_FPP1","",200,-1,7,200,-2,1);
  // TH1D *hp_mother_e_gamma_FPP1 = new TH1D("hp_mother_e_gamma_FPP1","",200,0.0,1.0);
  // TH1D *hp_mother_e_e_FPP1 = new TH1D("hp_mother_e_e_FPP1","",200,0.0,2.0);
  // TH1D *hp_mother_e_other_FPP1 = new TH1D("hp_mother_e_other_FPP1","",200,0.0,8.0);

  // TH2D *hvx_vz_mother_e_gamma_FPP2 = new TH2D("hvx_vz_mother_e_gamma_FPP2","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_mother_e_e_FPP2 = new TH2D("hvx_vz_mother_e_e_FPP2","",200,-1,7,200,-2,1);
  // TH2D *hvx_vz_mother_e_other_FPP2 = new TH2D("hvx_vz_mother_e_other_FPP2","",200,-1,7,200,-2,1);
  // TH1D *hp_mother_e_gamma_FPP2 = new TH1D("hp_mother_e_gamma_FPP2","",200,0.0,1.0);
  // TH1D *hp_mother_e_e_FPP2 = new TH1D("hp_mother_e_e_FPP2","",200,0.0,2.0);
  // TH1D *hp_mother_e_other_FPP2 = new TH1D("hp_mother_e_other_FPP2","",200,0.0,8.0);

  TH2D *hvy_vx_FT = new TH2D("hvy_vx_FT","Front tracker; hit vertex x (m); hit vertex y (m)",200,-3,3,200,-3,3);
  TH2D *hvy_vx_FPP1 = new TH2D("hvy_vx_FPP1","Back tracker; hit vertex x (m); hit vertex y (m)",200,-3,3,200,-3,3);
  //TH2D *hvy_vx_FPP2 = new TH2D("hvy_vx_FPP2","",200,-3,3,200,-3,3);

  TH2D *hvy_vz_FT = new TH2D("hvy_vz_FT","Front tracker; hit vertex z (m); hit vertex y (m)",200,-1,7,200,-3,3);
  TH2D *hvy_vz_FPP1 = new TH2D("hvy_vz_FPP1","Back tracker; hit vertex z (m); hit vertex y (m)",200,-1,7,200,-3,3);
  //TH2D *hvy_vz_FPP2 = new TH2D("hvy_vz_FPP2","",200,-1,7,200,-3,3);

  //SD track origin variables:
  TH2D *hSD_vx_vz_FT = new TH2D("hSD_vx_vz_FT","Front tracker; SD origin z (m); SD origin x (m)",200,-1,7,200,-2,1);
  TH2D *hSD_vy_vx_FT = new TH2D("hSD_vy_vx_FT","Front tracker; SD origin x (m); SD origin y (m)",200,-3,3,200,-3,3);
  TH2D *hSD_vy_vz_FT = new TH2D("hSD_vy_vz_FT","Front tracker; SD origin z (m); SD origin y (m)",200,-1,7,200,-3,3);
  
  //SD track origin variables:
  TH2D *hSD_vx_vz_FPP = new TH2D("hSD_vx_vz_FPP","Back tracker; SD origin z (m); SD origin x (m)",200,-1,7,200,-2,1);
  TH2D *hSD_vy_vx_FPP = new TH2D("hSD_vy_vx_FPP","Back tracker; SD origin x (m); SD origin y (m)",200,-3,3,200,-3,3);
  TH2D *hSD_vy_vz_FPP = new TH2D("hSD_vy_vz_FPP","Back tracker; SD origin z (m); SD origin y (m)",200,-1,7,200,-3,3);

  TH1D *hSD_vz_FT_all = new TH1D("hSD_vz_FT_all", "Front tracker all hits; vertex z (m); hit rate (Hz/cm^{2}/layer)", 200, -1, 7 );
  TH1D *hSD_vz_FT_gamma = new TH1D("hSD_vz_FT_gamma", "Front tracker #gamma-induced; vertex z (m); hit rate (Hz/cm^{2}/layer)", 200, -1, 7 );
  TH1D *hSD_vz_FT_electron = new TH1D("hSD_vz_FT_electron", "Front tracker e^{+}/e^{-}-induced; vertex z (m); hit rate (Hz/cm^{2}/layer)", 200, -1, 7 );
  TH1D *hSD_vz_FT_other = new TH1D("hSD_vz_FT_other", "Front-tracker hadron-induced; vertex z (m); hit rate (Hz/cm^{2}/layer)", 200, -1, 7 );

  TH1D *hSD_vz_FPP_all = new TH1D("hSD_vz_FPP_all", "Back tracker all hits; vertex z (m); hit rate (Hz/cm^{2}/layer)", 200, -1, 7 );
  TH1D *hSD_vz_FPP_gamma = new TH1D("hSD_vz_FPP_gamma", "Back tracker #gamma-induced; vertex z (m); hit rate (Hz/cm^{2}/layer)", 200, -1, 7 );
  TH1D *hSD_vz_FPP_electron = new TH1D("hSD_vz_FPP_electron", "Back tracker e^{+}/e^{-}-induced; vertex z (m); hit rate (Hz/cm^{2}/layer)", 200, -1, 7 );
  TH1D *hSD_vz_FPP_other = new TH1D("hSD_vz_FPP_other", "Back tracker hadron-induced; vertex z (m); hit rate (Hz/cm^{2}/layer)", 200, -1, 7 );
  
  double weight = 50e-6/1.602e-19/double(nevents_generated);
  // double weight_FT = weight / wplanes[0] / hplanes[0];
  // double weight_FPP1 = weight / wplanes[1] / hplanes[1];
  //double weight_FPP2 = weight / wplanes[2] / hplanes[2];

  while( T->GetEntry(nevent++) ){
    if( nevent%10000 == 0 ) cout << nevent << endl;
    
    int treenum = T->fChain->GetTreeNumber();
    
    TFile *ftemp = C->GetFile();
    TString fname = ftemp->GetName();
    if( bad_file_list.find(fname) == bad_file_list.end() ){
      
      for( int hit=0; hit<T->Harm_FT_hit_nhits; hit++ ){
	int plane = (*(T->Harm_FT_hit_plane))[hit];

	double edep_MeV = (*(T->Harm_FT_hit_edep))[hit] * 1000.0 ;
	
	double hitweight = weight/wplanes_FT[plane-1]/hplanes_FT[plane-1];
	if( (*(T->Harm_FT_hit_plane))[hit]  == 1 ){
	  hvx_vz_layer1_FT->Fill( (*(T->Harm_FT_hit_vz))[hit], (*(T->Harm_FT_hit_vx))[hit], hitweight );
	  //To get the rates versus position, we need to multiply the rate per unit area by the number of bins along
	  //that dimension:
	  
	  hrate_vs_x_layer1_FT->Fill( (*(T->Harm_FT_hit_x))[hit], hitweight*200.0);
	  hrate_vs_y_layer1_FT->Fill( (*(T->Harm_FT_hit_y))[hit], hitweight*200.0);
	  hrate_vs_p_layer1_FT_wide->Fill( (*(T->Harm_FT_hit_p))[hit], hitweight );
	  hrate_vs_p_layer1_FT_zoom->Fill( (*(T->Harm_FT_hit_p))[hit], hitweight );
	}
	hrate_layer_FT->Fill( (*(T->Harm_FT_hit_plane))[hit], hitweight );
	hrate_edep_layer_FT->Fill( plane, hitweight * edep_MeV );
	if( (*(T->Harm_FT_hit_pid))[hit] == 22 ){ //gammas:
	  hrate_vs_p_gamma_FT->Fill( (*(T->Harm_FT_hit_p))[hit], hitweight/double(nplanes_FT) );
	} else if( abs( (*(T->Harm_FT_hit_pid))[hit] ) == 11 ){ //e+/e-
	  hrate_vs_p_electron_FT->Fill( (*(T->Harm_FT_hit_p))[hit], hitweight/double(nplanes_FT) );
	} else if( fabs( (*(T->Harm_FT_hit_pid))[hit] ) == 211 ){ //pi+/pi-
	  hrate_vs_p_pion_FT->Fill( (*(T->Harm_FT_hit_p))[hit], hitweight/double(nplanes_FT) );
	} else if( (*(T->Harm_FT_hit_pid))[hit] == 2212 ){ //protons
	  hrate_vs_p_proton_FT->Fill( (*(T->Harm_FT_hit_p))[hit], hitweight/double(nplanes_FT) );
	} else { //all others
	  hrate_vs_p_other_FT->Fill( (*(T->Harm_FT_hit_p))[hit], hitweight/double(nplanes_FT) );
	}
	
	hedep_FT->Fill( (*(T->Harm_FT_hit_edep))[hit], hitweight/double(nplanes_FT) );
	if( abs( (*(T->Harm_FT_hit_pid))[hit] ) == 11 ){
	  hedep_FT_electron->Fill( (*(T->Harm_FT_hit_edep))[hit], hitweight/double(nplanes_FT) );
	} else {
	  hedep_FT_other->Fill( (*(T->Harm_FT_hit_edep))[hit], hitweight/double(nplanes_FT) );
	}
	//hedep_layer_FT->Fill( (*(T->Harm_FT_hit_plane))[hit], (*(T->Harm_FT_hit_edep))[hit], weight_FT );
	hedep_vs_p_FT->Fill( (*(T->Harm_FT_hit_p))[hit], (*(T->Harm_FT_hit_edep))[hit], hitweight/double(nplanes_FT) );
	
	hvx_vz_FT->Fill( (*(T->Harm_FT_hit_vz))[hit], (*(T->Harm_FT_hit_vx))[hit], hitweight/double(nplanes_FT) );
	hvy_vx_FT->Fill( (*(T->Harm_FT_hit_vx))[hit], (*(T->Harm_FT_hit_vy))[hit], hitweight/double(nplanes_FT) );
	hvy_vz_FT->Fill( (*(T->Harm_FT_hit_vz))[hit], (*(T->Harm_FT_hit_vy))[hit], hitweight/double(nplanes_FT) );

	int sdtridx = (*(T->Harm_FT_hit_sdtridx))[hit];
	
	hSD_vx_vz_FT->Fill( (*(T->SDTrack_vz))[sdtridx], (*(T->SDTrack_vx))[sdtridx], hitweight/double(nplanes_FT) );
	hSD_vy_vx_FT->Fill( (*(T->SDTrack_vx))[sdtridx], (*(T->SDTrack_vy))[sdtridx], hitweight/double(nplanes_FT) );
	hSD_vy_vz_FT->Fill( (*(T->SDTrack_vz))[sdtridx], (*(T->SDTrack_vy))[sdtridx], hitweight/double(nplanes_FT) );

	bool isgamma = ( (*(T->SDTrack_PID))[sdtridx] == 22 );
	bool iselectron = abs( (*(T->SDTrack_PID))[sdtridx] ) == 11;
	bool isother = !(isgamma || iselectron);

	hSD_vz_FT_all->Fill( (*(T->SDTrack_vz))[sdtridx], hitweight/double(nplanes_FT) );
	if( isgamma ){
	  hSD_vz_FT_gamma->Fill( (*(T->SDTrack_vz))[sdtridx], hitweight/double(nplanes_FT) );
	} else if( iselectron ){
	  hSD_vz_FT_electron->Fill( (*(T->SDTrack_vz))[sdtridx], hitweight/double(nplanes_FT) );
	} else {
	  hSD_vz_FT_other->Fill( (*(T->SDTrack_vz))[sdtridx], hitweight/double(nplanes_FT) );
	}
	  
      }

      for( int hit=0; hit<T->Harm_FPP1_hit_nhits; hit++ ){

	int plane = (*(T->Harm_FPP1_hit_plane))[hit];

	double edep_MeV = (*(T->Harm_FPP1_hit_edep))[hit] * 1000.0;
	
	double hitweight = weight/wplanes_FPP[plane-1]/hplanes_FPP[plane-1];
	
	hrate_layer_FPP1->Fill( (*(T->Harm_FPP1_hit_plane))[hit], hitweight );
	hrate_edep_layer_FPP1->Fill( plane, hitweight * edep_MeV );
	
	if( (*(T->Harm_FPP1_hit_pid))[hit] == 22 ){ 
	  hrate_vs_p_gamma_FPP1->Fill( (*(T->Harm_FPP1_hit_p))[hit], hitweight/double(nplanes_FPP) );
	} else if( abs( (*(T->Harm_FPP1_hit_pid))[hit] ) == 11 ){
	  hrate_vs_p_electron_FPP1->Fill( (*(T->Harm_FPP1_hit_p))[hit], hitweight/double(nplanes_FPP) );
	} else if( abs( (*(T->Harm_FPP1_hit_pid))[hit] ) == 211 ){
	  hrate_vs_p_pion_FPP1->Fill( (*(T->Harm_FPP1_hit_p))[hit], hitweight/double(nplanes_FPP) );
	} else if( (*(T->Harm_FPP1_hit_pid))[hit] == 2212 ){
	  hrate_vs_p_proton_FPP1->Fill( (*(T->Harm_FPP1_hit_p))[hit], hitweight/double(nplanes_FPP) );
	} else {
	  hrate_vs_p_other_FPP1->Fill( (*(T->Harm_FPP1_hit_p))[hit], hitweight/double(nplanes_FPP) );
	}

	hedep_FPP1->Fill( (*(T->Harm_FPP1_hit_edep))[hit], hitweight/double(nplanes_FPP) );
	//hedep_layer_FPP1->Fill( (*(T->Harm_FPP1_hit_plane))[hit], (*(T->Harm_FPP1_hit_edep))[hit], weight_FPP1 );
	if( abs( (*(T->Harm_FPP1_hit_pid))[hit] ) == 11 ){
	  hedep_FPP1_electron->Fill( (*(T->Harm_FPP1_hit_edep))[hit], hitweight/double(nplanes_FPP) );
	} else {
	  hedep_FPP1_other->Fill( (*(T->Harm_FPP1_hit_edep))[hit], hitweight/double(nplanes_FPP) );
	}
	hedep_vs_p_FPP1->Fill( (*(T->Harm_FPP1_hit_p))[hit], (*(T->Harm_FPP1_hit_edep))[hit], hitweight/double(nplanes_FPP) );

	hvx_vz_FPP1->Fill( (*(T->Harm_FPP1_hit_vz))[hit], (*(T->Harm_FPP1_hit_vx))[hit], hitweight/double(nplanes_FPP) );
	hvy_vx_FPP1->Fill( (*(T->Harm_FPP1_hit_vx))[hit], (*(T->Harm_FPP1_hit_vy))[hit], hitweight/double(nplanes_FPP) );
	hvy_vz_FPP1->Fill( (*(T->Harm_FPP1_hit_vz))[hit], (*(T->Harm_FPP1_hit_vy))[hit], hitweight/double(nplanes_FPP) );

	int sdtridx = (*(T->Harm_FPP1_hit_sdtridx))[hit];
	
	hSD_vx_vz_FPP->Fill( (*(T->SDTrack_vz))[sdtridx], (*(T->SDTrack_vx))[sdtridx], hitweight/double(nplanes_FPP) );
	hSD_vy_vx_FPP->Fill( (*(T->SDTrack_vx))[sdtridx], (*(T->SDTrack_vy))[sdtridx], hitweight/double(nplanes_FPP) );
	hSD_vy_vz_FPP->Fill( (*(T->SDTrack_vz))[sdtridx], (*(T->SDTrack_vy))[sdtridx], hitweight/double(nplanes_FPP) );

	bool isgamma = ( (*(T->SDTrack_PID))[sdtridx] == 22 );
	bool iselectron = abs( (*(T->SDTrack_PID))[sdtridx] ) == 11;
	bool isother = !(isgamma || iselectron);

	hSD_vz_FPP_all->Fill( (*(T->SDTrack_vz))[sdtridx], hitweight/double(nplanes_FPP) );
	if( isgamma ){
	  hSD_vz_FPP_gamma->Fill( (*(T->SDTrack_vz))[sdtridx], hitweight/double(nplanes_FPP) );
	} else if( iselectron ){
	  hSD_vz_FPP_electron->Fill( (*(T->SDTrack_vz))[sdtridx], hitweight/double(nplanes_FPP) );
	} else {
	  hSD_vz_FPP_other->Fill( (*(T->SDTrack_vz))[sdtridx], hitweight/double(nplanes_FPP) );
	}
	
      }

      for( int ihit=0; ihit<T->Earm_CDET_Scint_hit_nhits; ihit++ ){
	double edep = (*(T->Earm_CDET_Scint_hit_sumedep))[ihit];

	double hitweight = weight;
	
	if( edep >= thresh_CDET ){
	  hrate_versus_channel_CDET->Fill( (*(T->Earm_CDET_Scint_hit_cell))[ihit], hitweight );
	}
      }
      
    }
  }

  //C->Delete();
  fout->cd();
  fout->Write();
}
