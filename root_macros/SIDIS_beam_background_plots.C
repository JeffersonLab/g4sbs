#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "sidis_tree.C"
#include "G4SBSRunData.hh"
#include "TObjArray.h"
#include "TChainElement.h"

void SIDIS_beam_background_plots(const char *rootfilename, const char *outfilename, double Ibeam=60e-6 ){
  TChain *C = new TChain("T");

  C->Add(rootfilename);

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
  }
  
  fout->Write();
  
}
