#include "TTree.h"
#include "TChain.h"
#include "G4SBSRunData.hh"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "WAPP_tree.C"
#include "TObjArray.h"
#include "TChainElement.h"

#include <iostream>
#include <vector>

using namespace std;

void plot_GENRP_GEMbackgrounds(const char *infilename, const char *outfilename, double Ibeam = 30e-6){
  TH1::SetDefaultSumw2();

  TFile *fout = new TFile(outfilename,"RECREATE");
  
  TChain *C = new TChain("T");

  C->Add(infilename);

  //int nfiles = C->GetNtrees();

  TObjArray *FileList = C->GetListOfFiles();

  long ngen_total = 0;
  
  for( int ifile=0; ifile<FileList->GetEntries(); ifile++ ){

    TString fname = ( (TChainElement*) (*FileList)[ifile] )->GetTitle();

    TFile *ftemp = new TFile(fname, "READ");

    G4SBSRunData *rd;

    ftemp->GetObject("run_data",rd);

    ngen_total += rd->fNtries;

    cout << "Added file " << fname << " ngen total = " << ngen_total << endl;

    ftemp->Close();
    ftemp->Delete();
  }

  cout << "Finished adding files, ngen = " << ngen_total << endl;

  //Now we have BB GEMs, and GEN-RP GEMs to worry about initially:

  double BBGEM_area_cm2[5] = {40.0*150.0, 40.0*150.0, 40.0*150.0, 40.0*150.0, 60.0*200.0 };
  double CEPolFrontGEM_area_cm2[4] = {40.0*150.0, 40.0*150.0, 60.0*200.0, 60.0*200.0 };
  double CEPolRearGEM_area_cm2[4] = {60.0*200.0, 60.0*200.0, 60.0*200.0, 60.0*200.0 };
  double PRPolBSGEM_area_cm2[2] = {60.0*200.0, 60.0*200.0};
  double PRPolFSGEM_area_cm2[2] = {60.0*200.0, 60.0*200.0};

  double BBGEM_LX[5] = {150.,150.,150.,150.,200.};
  double CEPolFrontGEM_LX[4] = {150.,150.,200.,200.};
  double CEPolRearGEM_LX[4] = {200.,200.,200.,200.};
  double PRPolBSGEM_LX[2] = {200.,200.};
  double PRPolFSGEM_LX[2] = {200.,200.};

  double BBGEM_LY[5] = {40.,40.,40.,40.,60.};
  double CEPolFrontGEM_LY[4] = {40.,40.,60.,60.};
  double CEPolRearGEM_LY[4] = {60.,60.,60.,60.};
  double PRPolBSGEM_LY[2] = {60.,60.};
  double PRPolFSGEM_LY[2] = {60.,60.};
  
  fout->cd();
  
  TH1D *hitrate_vs_layer_BBGEM = new TH1D("hitrate_vs_layer_BBGEM","",5,0.5,5.5);
  TH1D *hitrate_vs_layer_CEPolFrontGEM = new TH1D("hitrate_vs_layer_CEPolFrontGEM","",4,0.5,4.5);
  TH1D *hitrate_vs_layer_CEPolRearGEM = new TH1D("hitrate_vs_layer_CEPolRearGEM","",4,0.5,4.5);
  TH1D *hitrate_vs_layer_PRPolBeamSideGEM = new TH1D("hitrate_vs_layer_PRPolBeamSideGEM","",2,0.5,2.5);
  TH1D *hitrate_vs_layer_PRPolFarSideGEM = new TH1D("hitrate_vs_layer_PRPolFarSideGEM","",2,0.5,2.5);

  TH2D *hitrate_vs_X_BBGEM = new TH2D("hitrate_vs_X_BBGEM","Hit rate (Hz/cm^{2})",5,0.5,5.5,100,-1.05,1.05);
  TH2D *hitrate_vs_X_CEPolFrontGEM = new TH2D("hitrate_vs_X_CEPolFrontGEM","Hit rate (Hz/cm^{2})",4,0.5,4.5,100,-1.05,1.05);
  TH2D *hitrate_vs_X_CEPolRearGEM = new TH2D("hitrate_vs_X_CEPolRearGEM","Hit rate (Hz/cm^{2})",4,0.5,4.5,100,-1.05,1.05);
  TH2D *hitrate_vs_X_PRPolBeamSideGEM = new TH2D("hitrate_vs_X_PRPolBeamSideGEM","Hit rate (Hz/cm^{2})",4,0.5,4.5,100,-1.05,1.05);
  TH2D *hitrate_vs_X_PRPolFarSideGEM = new TH2D("hitrate_vs_X_PRPolFarSideGEM","Hit rate (Hz/cm^{2})",4,0.5,4.5,100,-1.05,1.05);

  TH2D *hitrate_vs_Y_BBGEM = new TH2D("hitrate_vs_Y_BBGEM","Hit rate (Hz/cm^{2})",5,0.5,5.5,100,-0.31,0.31);
  TH2D *hitrate_vs_Y_CEPolFrontGEM = new TH2D("hitrate_vs_Y_CEPolFrontGEM","Hit rate (Hz/cm^{2})",4,0.5,4.5,100,-0.31,0.31);
  TH2D *hitrate_vs_Y_CEPolRearGEM = new TH2D("hitrate_vs_Y_CEPolRearGEM","Hit rate (Hz/cm^{2})",4,0.5,4.5,100,-0.31,0.31);
  TH2D *hitrate_vs_Y_PRPolBeamSideGEM = new TH2D("hitrate_vs_Y_PRPolBeamSideGEM","Hit rate (Hz/cm^{2})",4,0.5,4.5,100,-0.31+1.36,0.31+1.36);
  TH2D *hitrate_vs_Y_PRPolFarSideGEM = new TH2D("hitrate_vs_Y_PRPolFarSideGEM","Hit rate (Hz/cm^{2})",4,0.5,4.5,100,-0.31+1.36,0.31+1.36);

  hitrate_vs_X_BBGEM->SetXTitle("BigBite GEM layer");
  hitrate_vs_X_BBGEM->SetYTitle("Hit X (m)");

  hitrate_vs_Y_BBGEM->SetXTitle("BigBite GEM layer");
  hitrate_vs_Y_BBGEM->SetYTitle("Hit Y (m)");

  hitrate_vs_X_CEPolFrontGEM->SetXTitle("C.E. Pol. Front GEM layer");
  hitrate_vs_X_CEPolFrontGEM->SetYTitle("Hit X (m)");

  hitrate_vs_Y_CEPolFrontGEM->SetXTitle("C.E. Pol. Front GEM layer");
  hitrate_vs_Y_CEPolFrontGEM->SetYTitle("Hit Y (m)");

  hitrate_vs_X_CEPolRearGEM->SetXTitle("C.E. Pol. Rear GEM layer");
  hitrate_vs_X_CEPolRearGEM->SetYTitle("Hit X (m)");

  hitrate_vs_Y_CEPolRearGEM->SetXTitle("C.E. Pol. Rear GEM layer");
  hitrate_vs_Y_CEPolRearGEM->SetYTitle("Hit Y (m)");

  hitrate_vs_X_PRPolBeamSideGEM->SetXTitle("P.R. Pol. Beam Side GEM layer");
  hitrate_vs_X_PRPolBeamSideGEM->SetYTitle("Hit X (m)");

  hitrate_vs_Y_PRPolBeamSideGEM->SetXTitle("P.R. Pol. Beam Side GEM layer");
  hitrate_vs_Y_PRPolBeamSideGEM->SetYTitle("Hit Y (m)");

  hitrate_vs_X_PRPolFarSideGEM->SetXTitle("P.R. Pol. Far Side GEM layer");
  hitrate_vs_X_PRPolFarSideGEM->SetYTitle("Hit X (m)");

  hitrate_vs_Y_PRPolFarSideGEM->SetXTitle("P.R. Pol. Far Side GEM layer");
  hitrate_vs_Y_PRPolFarSideGEM->SetYTitle("Hit Y (m)");
  
  double Xbinwidth = 210.0/100.0;
  double Ybinwidth = 62.0/100.0;
  
  WAPP_tree *T = new WAPP_tree(C);

  //Now we're ready to loop over all the relevant detectors:

  long nevent=0; 
  
  while( C->GetEntry( nevent ++ ) ){

    //    cout << "Starting event loop, nevent = " << nevent << endl;
    
    if( nevent % 1000 == 0 ) cout << "event " << nevent << endl;
    
    double weight = Ibeam/double(ngen_total)/1.602e-19;

    for( int ihit=0; ihit<T->Earm_BBGEM_hit_nhits; ihit++ ){
      int plane = (*(T->Earm_BBGEM_hit_plane))[ihit];
      hitrate_vs_layer_BBGEM->Fill( plane, weight/BBGEM_area_cm2[plane-1] );
      hitrate_vs_X_BBGEM->Fill( plane, (*(T->Earm_BBGEM_hit_x))[ihit], weight/Xbinwidth/BBGEM_LY[plane-1] );
      hitrate_vs_Y_BBGEM->Fill( plane, (*(T->Earm_BBGEM_hit_y))[ihit], weight/Ybinwidth/BBGEM_LX[plane-1] );
    }

    for( int ihit=0; ihit<T->Harm_CEPolFront_hit_nhits; ihit++ ){
      int plane = (*(T->Harm_CEPolFront_hit_plane))[ihit];
      hitrate_vs_layer_CEPolFrontGEM->Fill( plane, weight/CEPolFrontGEM_area_cm2[plane-1] );
      hitrate_vs_X_CEPolFrontGEM->Fill( plane, (*(T->Harm_CEPolFront_hit_x))[ihit], weight/Xbinwidth/CEPolFrontGEM_LY[plane-1] );
      hitrate_vs_Y_CEPolFrontGEM->Fill( plane, (*(T->Harm_CEPolFront_hit_y))[ihit], weight/Ybinwidth/CEPolFrontGEM_LX[plane-1] );
    }

    for( int ihit=0; ihit<T->Harm_CEPolRear_hit_nhits; ihit++ ){
      int plane = (*(T->Harm_CEPolRear_hit_plane))[ihit];
      hitrate_vs_layer_CEPolRearGEM->Fill( plane, weight/CEPolRearGEM_area_cm2[plane-1] );
      hitrate_vs_X_CEPolRearGEM->Fill( plane, (*(T->Harm_CEPolRear_hit_x))[ihit], weight/Xbinwidth/CEPolRearGEM_LY[plane-1] );
      hitrate_vs_Y_CEPolRearGEM->Fill( plane, (*(T->Harm_CEPolRear_hit_y))[ihit], weight/Ybinwidth/CEPolRearGEM_LX[plane-1] );
    }

    for( int ihit=0; ihit<T->Harm_PRPolGEMBeamSide_hit_nhits; ihit++ ){
      int plane = (*(T->Harm_PRPolGEMBeamSide_hit_plane))[ihit];
      hitrate_vs_layer_PRPolBeamSideGEM->Fill( plane, weight/PRPolBSGEM_area_cm2[plane-1] );
      hitrate_vs_X_PRPolBeamSideGEM->Fill( plane, (*(T->Harm_PRPolGEMBeamSide_hit_x))[ihit], weight/Xbinwidth/PRPolBSGEM_LY[plane-1] );
      hitrate_vs_Y_PRPolBeamSideGEM->Fill( plane, (*(T->Harm_PRPolGEMBeamSide_hit_z))[ihit], weight/Ybinwidth/PRPolBSGEM_LX[plane-1] );
    }

    for( int ihit=0; ihit<T->Harm_PRPolGEMFarSide_hit_nhits; ihit++ ){
      int plane = (*(T->Harm_PRPolGEMFarSide_hit_plane))[ihit];
      hitrate_vs_layer_PRPolFarSideGEM->Fill( plane, weight/PRPolFSGEM_area_cm2[plane-1] );
      hitrate_vs_X_PRPolFarSideGEM->Fill( plane, (*(T->Harm_PRPolGEMFarSide_hit_x))[ihit], weight/Xbinwidth/PRPolFSGEM_LY[plane-1] );
      hitrate_vs_Y_PRPolFarSideGEM->Fill( plane, (*(T->Harm_PRPolGEMFarSide_hit_z))[ihit], weight/Ybinwidth/PRPolFSGEM_LX[plane-1] );
    }
    
  }
  
  fout->Write();

}
