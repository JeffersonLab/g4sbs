#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "sidis_tree.C"
#include "G4SBSRunData.hh"
#include "TObjArray.h"
#include "TChainElement.h"
#include "TRandom3.h"
#include <iostream>

void SIDIS_beam_background_plots(const char *rootfilename, const char *outfilename, double Ibeam=60e-6 ){
  TChain *C = new TChain("T");

  C->Add(rootfilename);

  TRandom3 num(0);

  double ngen_total = 0.0;

  TObjArray *filelist = C->GetListOfFiles();

  TIter next(filelist);

  //The only metadata we need to grab for this purpose is the number of generated events:
  TChainElement *chEl=0;
  while( (chEl = (TChainElement*) next())) {
    TFile F(chEl->GetTitle(),"READ");

    G4SBSRunData *rd;

    F.GetObject("run_data",rd);
    
    ngen_total += double(rd->fNtries);

    F.Close();
  }

  cout << "Total number of generated events = " << ngen_total << endl;

  sidis_tree *T = new sidis_tree(C);

  TFile *fout = new TFile(outfilename,"RECREATE");
  //Create the basic histograms. We will also want to look at the sources of the background as well:
  TH1D *hBBGEM_rate_vs_layer = new TH1D("hBBGEM_rate_vs_layer","BigBite GEM rates vs layer", 5, 0.5, 5.5);
  TH1D *hSBSGEM_rate_vs_layer = new TH1D("hSBSGEM_rate_vs_layer","SBS GEM rates vs layer", 5, 0.5, 5.5);

  TH1D *hBBGEM_occupancy_vs_layer = new TH1D("hBBGEM_occupancy_vs_layer","BigBite GEM occupancies vs layer", 5, 0.5, 5.5);
  TH1D *hSBSGEM_occupancy_vs_layer = new TH1D("hSBSGEM_occupancy_vs_layer","SBS GEM occupancies vs layer", 5, 0.5, 5.5);
  
  TH2D *hBBGEM_rate_vs_x_layer = new TH2D("hBBGEM_rate_vs_x_layer","BigBite GEM rates vs x by layer (Hz/cm^{2})", 5,0.5,5.5, 100,-110.,110.);
  TH2D *hSBSGEM_rate_vs_x_layer = new TH2D("hSBSGEM_rate_vs_x_layer","SBS GEM rates vs x by layer (Hz/cm^{2})", 5,0.5,5.5,100,-110.,110.);

  TH2D *hBBGEM_occupancy_vs_x_layer = new TH2D("hBBGEM_occupancy_vs_x_layer","BigBite GEM occupancy vs x by layer", 5,0.5,5.5, 100,-110.,110.);
  TH2D *hSBSGEM_occupancy_vs_x_layer = new TH2D("hSBSGEM_occupancy_vs_x_layer","SBS GEM occupancy vs x by layer", 5,0.5,5.5,100,-110.,110.);

  TH2D *hBBGEM_rate_vs_y_layer = new TH2D("hBBGEM_rate_vs_y_layer","BigBite GEM rates vs y by layer (Hz/cm^{2})", 5,0.5,5.5, 100,-35.,35.);
  TH2D *hSBSGEM_rate_vs_y_layer = new TH2D("hSBSGEM_rate_vs_y_layer","SBS GEM rates vs y by layer (Hz/cm^{2})", 5,0.5,5.5,100,-35.,35.);

  TH2D *hBBGEM_occupancy_vs_y_layer = new TH2D("hBBGEM_occupancy_vs_y_layer","BigBite GEM occupancy vs x by layer", 5,0.5,5.5, 100,-35.,35.);
  TH2D *hSBSGEM_occupancy_vs_y_layer = new TH2D("hSBSGEM_occupancy_vs_y_layer","SBS GEM occupancy vs x by layer", 5,0.5,5.5,100,-35.,35.);
  
  TH2D *hBBGEM_vx_vz_background = new TH2D("hBBGEM_vx_vz_background", "origin x vs z of GEM background hits", 200, -1.,6., 200, -2.5, 2.5);
  TH2D *hBBGEM_vy_vz_background = new TH2D("hBBGEM_vy_vz_background", "origin y vs z of GEM background hits", 200, -1.,6., 200, -2.5, 2.5);
  TH2D *hBBGEM_vy_vx_background = new TH2D("hBBGEM_vy_vx_background", "origin y vs x of GEM background hits", 200, -2,2,200,-2,2);
  
  TH1D *hrate_vs_PMT_RICH = new TH1D("hrate_vs_PMT_RICH", "RICH PMT counting rates", 1934, 0.5,1934.5 );
  TH1D *hrate_vs_PMT_GRINCH = new TH1D("hrate_vs_PMT_GRINCH", "GRINCH PMT counting rates", 550, 0.5, 550.5 );
  TH1D *hRICH_PMT_occupancy = new TH1D("hRICH_PMT_occupancy", "RICH PMT occupancy (10 ns window)", 1934, 0.5, 1934.5 );
  TH1D *hGRINCH_PMT_occupancy = new TH1D("hGRINCH_PMT_occupancy", "GRINCH PMT occupancy (10 ns window)", 550, 0.5, 550.5 );

  TH2D *hRICH_vx_vz_background = new TH2D("hRICH_vx_vz_background", "origin x vs z of RICH background hits", 200, -1, 6, 200, -2.5, 2.5 );
  TH2D *hRICH_vy_vz_background = new TH2D("hRICH_vy_vz_background", "origin y vs z of RICH background hits", 200, -1, 6, 200, -2.5, 2.5 );
  TH2D *hRICH_vy_vx_background = new TH2D("hRICH_vy_vx_background", "origin y vs x of RICH background hits", 200, -2, 2, 200, -2, 2 );

  TH1D *hrate_vs_edep_HCAL = new TH1D("hrate_vs_edep_HCAL", "Rate vs energy deposit (unintegrated, Hz/(10 MeV))", 200, 0.0, 2.0 );
  TH1D *hrate_vs_edep_BBCAL = new TH1D("hrate_vs_edep_BBCAL", "Rate vs energy deposit (unintegrated, Hz/(10 MeV))", 500, 0.0, 5.0 );

  TH2D *hHCAL_vx_vz_background = new TH2D("hHCAL_vx_vz_background", "origin x vs z of HCAL background hits", 200, -1, 9, 200, -3, 3 );
  TH2D *hBBPS_vx_vz_background = new TH2D("hBBPS_vx_vz_background", "origin x vs z of BB pre-shower background hits", 200, -1, 6, 200, -2.5, 2.5 );
  TH2D *hBBSH_vx_vz_background = new TH2D("hBBSH_vx_vz_background", "origin x vs z of BB shower background hits", 200, -1, 6, 200, -2.5, 2.5 );

  double bbgem_area_cm2[5] = {40.*150.,40.*150.,40.*150.,40.*150.,60.*200.};
  double sbsgem_area_cm2[5] = {60.*200.,60.*200.,60.*200.,60.*200.,60.*200.};

  double strip_pitch = 0.04; //0.4 mm

  int bbgem_nmodules[5] = {3,3,3,3,4};
  int sbsgem_nmodules[5] = {4,4,4,4,4};
  
  double bbgem_Lx[5] = {150., 150., 150., 150., 200.};
  double bbgem_Ly[5] = {40., 40., 40., 40., 60.};
  double sbsgem_Lx[5] = {200., 200., 200., 200., 200.};
  double sbsgem_Ly[5] = {60., 60., 60., 60., 60. };
  
  double Xbinwidth = 220.0/100.0;
  double Ybinwidth = 70.0/100.0;
  
  long nevent=0; 

  double weight_rate = Ibeam/1.602e-19/ngen_total;

  //For occupancy, we calculate the total hit rate times the average strip multiplicity per hit divided by the total number of strips 
  
  while( C->GetEntry( nevent++ ) ){
    if ( nevent % 1000 == 0 ) cout << "event " << nevent << endl;
    for( int ihit=0; ihit<T->Earm_BBGEM_hit_nhits; ihit++ ){
      int layer = (*(T->Earm_BBGEM_hit_plane))[ihit];
      double xhit = (*(T->Earm_BBGEM_hit_x))[ihit]*100.;
      double yhit = (*(T->Earm_BBGEM_hit_y))[ihit]*100.;
      hBBGEM_rate_vs_layer->Fill( layer, weight_rate/bbgem_area_cm2[layer-1] );  
      hBBGEM_rate_vs_x_layer->Fill( layer, xhit, weight_rate/Xbinwidth/bbgem_Ly[layer-1] );
      hBBGEM_rate_vs_y_layer->Fill( layer, yhit, weight_rate/Ybinwidth/bbgem_Lx[layer-1] );
      hBBGEM_occupancy_vs_layer->Fill( layer, weight_rate * 4. * 325e-9 * strip_pitch/bbgem_Lx[layer-1] );
      hBBGEM_occupancy_vs_x_layer->Fill( layer, xhit, weight_rate * 4. * 325e-9 * strip_pitch/Xbinwidth );
      hBBGEM_occupancy_vs_y_layer->Fill( layer, yhit, weight_rate * 4. * 325e-9 * strip_pitch/Ybinwidth/double(bbgem_nmodules[layer-1]) );

      hBBGEM_vx_vz_background->Fill( (*(T->SDTrack_vz))[(*(T->Earm_BBGEM_hit_sdtridx))[ihit]], (*(T->SDTrack_vx))[(*(T->Earm_BBGEM_hit_sdtridx))[ihit]], weight_rate );  
      hBBGEM_vy_vz_background->Fill( (*(T->SDTrack_vz))[(*(T->Earm_BBGEM_hit_sdtridx))[ihit]], (*(T->SDTrack_vy))[(*(T->Earm_BBGEM_hit_sdtridx))[ihit]], weight_rate );
      hBBGEM_vy_vx_background->Fill( (*(T->SDTrack_vx))[(*(T->Earm_BBGEM_hit_sdtridx))[ihit]], (*(T->SDTrack_vy))[(*(T->Earm_BBGEM_hit_sdtridx))[ihit]], weight_rate );
    }

    for( int ihit=0; ihit<T->Harm_SBSGEM_hit_nhits; ihit++ ){
      int layer = (*(T->Harm_SBSGEM_hit_plane))[ihit];
      double xhit = (*(T->Harm_SBSGEM_hit_x))[ihit]*100.;
      double yhit = (*(T->Harm_SBSGEM_hit_y))[ihit]*100.;
      hSBSGEM_rate_vs_layer->Fill( layer, weight_rate/sbsgem_area_cm2[layer-1] );  
      hSBSGEM_rate_vs_x_layer->Fill( layer, xhit, weight_rate/Xbinwidth/sbsgem_Ly[layer-1] );
      hSBSGEM_rate_vs_y_layer->Fill( layer, yhit, weight_rate/Ybinwidth/sbsgem_Lx[layer-1] );
      hSBSGEM_occupancy_vs_layer->Fill( layer, weight_rate * 4. * 325e-9 * strip_pitch/sbsgem_Lx[layer-1] );
      hSBSGEM_occupancy_vs_x_layer->Fill( layer, xhit, weight_rate * 4. * 325e-9 * strip_pitch/Xbinwidth );
      hSBSGEM_occupancy_vs_y_layer->Fill( layer, yhit, weight_rate * 4. * 325e-9 * strip_pitch/Ybinwidth/double(sbsgem_nmodules[layer-1]) );
    }
    
    for( int ihit=0; ihit<T->Earm_GRINCH_hit_nhits; ihit++ ){
      int PMT = (*(T->Earm_GRINCH_hit_col))[ihit] + 9*(*(T->Earm_GRINCH_hit_row))[ihit];
      hrate_vs_PMT_GRINCH->Fill( PMT+1, weight_rate );
      hGRINCH_PMT_occupancy->Fill( PMT+1, weight_rate * 10e-9 );
    }

    for( int ihit=0; ihit<T->Harm_RICH_hit_nhits; ihit++ ){
      int PMT = (*(T->Harm_RICH_hit_PMT))[ihit];
      hrate_vs_PMT_RICH->Fill( PMT+1, weight_rate );
      hRICH_PMT_occupancy->Fill( PMT+1, weight_rate*10e-9 );

      hRICH_vx_vz_background->Fill( (*(T->SDTrack_vz))[(*(T->Harm_RICH_hit_sdtridx))[ihit]], (*(T->SDTrack_vx))[(*(T->Harm_RICH_hit_sdtridx))[ihit]], weight_rate );
      hRICH_vy_vz_background->Fill( (*(T->SDTrack_vz))[(*(T->Harm_RICH_hit_sdtridx))[ihit]], (*(T->SDTrack_vy))[(*(T->Harm_RICH_hit_sdtridx))[ihit]], weight_rate );
      hRICH_vy_vx_background->Fill( (*(T->SDTrack_vx))[(*(T->Harm_RICH_hit_sdtridx))[ihit]], (*(T->SDTrack_vy))[(*(T->Harm_RICH_hit_sdtridx))[ihit]], weight_rate );
    }

    //For HCAL we can just use the sum of energy deposit in the detector:
    hrate_vs_edep_HCAL->Fill( T->Harm_HCalScint_det_esum, weight_rate );

    //For BigBite calorimeter, we also want to smear for photoelectron statistics:
    //we'll use a canonical assumption of 300 pe/GeV

    double edep_PS = T->Earm_BBPSTF1_det_esum;
    double edep_SH = T->Earm_BBSHTF1_det_esum;

    double mean_npe_PS = 300.*edep_PS;
    double mean_npe_SH = 300.*edep_SH;

    double npesmear_PS = num.Gaus( mean_npe_PS, sqrt(mean_npe_PS) );
    double npesmear_SH = num.Gaus( mean_npe_SH, sqrt(mean_npe_SH) );

    double edep_smear_PS = npesmear_PS/300.;
    double edep_smear_SH = npesmear_SH/300.;

    double etot_BBCAL = edep_smear_PS + edep_smear_SH;

    hrate_vs_edep_BBCAL->Fill( etot_BBCAL, weight_rate );

    if( T->Harm_HCalScint_det_esum >= 0.1 ){ //plot x vs z distribution of the background:
      for( int ihit=0; ihit<T->Harm_HCalScint_hit_nhits; ihit++ ){
	hHCAL_vx_vz_background->Fill( (*(T->SDTrack_vz))[(*(T->Harm_HCalScint_hit_sdtridx))[ihit]], (*(T->SDTrack_vx))[(*(T->Harm_HCalScint_hit_sdtridx))[ihit]], 
				      weight_rate );
      }
    }

    if( etot_BBCAL >= 0.8 ){ //plot x vs z distribution of the background: 
      for( int ihit=0; ihit<T->Earm_BBPSTF1_hit_nhits; ihit++ ){
	hBBPS_vx_vz_background->Fill( (*(T->SDTrack_vz))[(*(T->Earm_BBPSTF1_hit_sdtridx))[ihit]],
				      (*(T->SDTrack_vx))[(*(T->Earm_BBPSTF1_hit_sdtridx))[ihit]], 
				      weight_rate );
      }

      for( int ihit=0; ihit<T->Earm_BBSHTF1_hit_nhits; ihit++ ){
	hBBSH_vx_vz_background->Fill( (*(T->SDTrack_vz))[(*(T->Earm_BBSHTF1_hit_sdtridx))[ihit]],
				      (*(T->SDTrack_vx))[(*(T->Earm_BBSHTF1_hit_sdtridx))[ihit]], 
				      weight_rate );
      }
      
    }
  }
  
  fout->Write();
  
}
