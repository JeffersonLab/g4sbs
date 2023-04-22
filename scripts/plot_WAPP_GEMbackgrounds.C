#include "TTree.h"
#include "TChain.h"
#include "G4SBSRunData.hh"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "WAPP_tree.C"
#include "TObjArray.h"
#include "TChainElement.h"
#include "TRandom3.h"

#include <iostream>
#include <vector>

using namespace std;

void plot_GENRP_GEMbackgrounds(const char *infilename, const char *outfilename, double Ibeam = 30e-6, double thresh_PS=0.1, double thresh_SH=0.5, double thresh_HCAL=0.08){
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

  TH1D *hedep_HCAL = new TH1D("hedep_HCAL","",250,0.0,0.5);
  TH1D *hedep_BBPS = new TH1D("hedep_BBPS","",250,0.0,1.0);
  TH1D *hedep_BBSH = new TH1D("hedep_BBSH","",250,0.0,2.5);

  TH1D *hedep_BBSH_cutPS = new TH1D("hedep_BBSH_cutPS","",250,0.0,2.5);
  TH1D *hedep_BBSH_anticutPS = new TH1D("hedep_BBSH_anticutPS","",250,0.0,2.5);
  TH1D *hedep_BBPS_cutSH = new TH1D("hedep_BBPS_cutSH","",250,0.0,2.5);
  TH1D *hedep_BBPS_anticutSH = new TH1D("hedep_BBPS_anticutSH","",250,0.0,2.5);

  
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

  double rate_BB_pion = 0.0;
  double rate_BB_electron = 0.0;
  double rate_HCAL = 0.0;

  TRandom3 num(0);
  
  while( C->GetEntry( nevent++ ) ){

    //    cout << "Starting event loop, nevent = " << nevent << endl;
    
    if( nevent % 1000 == 0 ) cout << "event " << nevent << endl;
    
    double weight = Ibeam/double(ngen_total)/1.602e-19;

    hedep_HCAL->Fill( T->Harm_HCalScint_det_esum, weight );
    // hedep_BBPS->Fill( T->Earm_BBPSTF1_det_esum, weight );
    // hedep_BBSH->Fill( T->Earm_BBSHTF1_det_esum, weight );

    double BBSH_sums[26];
    double BBPS_sums[25];
    for( int i=0; i<26; i++ ){
      BBSH_sums[i] = 0.0;
      if( i<25 ) BBPS_sums[i] = 0.0;
    }

    for( int ihit=0; ihit<T->Earm_BBPSTF1_hit_nhits; ihit++ ){
      int row = (*(T->Earm_BBPSTF1_hit_row))[ihit];
      double edep = (*(T->Earm_BBPSTF1_hit_sumedep))[ihit];

      double npe_mean = 300.0*edep;
      double npe_smear = num.Gaus( npe_mean, sqrt(npe_mean) );

      double esmear = npe_smear/300.0;

      if( row < 25 ){
	BBPS_sums[row] += esmear;
      }
      if( row > 0 ) {
	BBPS_sums[row-1] += esmear; 
      }
    }

    for( int ihit=0; ihit<T->Earm_BBSHTF1_hit_nhits; ihit++ ){
      int row = (*(T->Earm_BBSHTF1_hit_row))[ihit];
      double edep = (*(T->Earm_BBSHTF1_hit_sumedep))[ihit];

      double npe_mean = 300.0*edep;
      double npe_smear = num.Gaus( npe_mean, sqrt(npe_mean) );

      double esmear = npe_smear/300.0;
      
      if( row < 26 ) BBSH_sums[row] += esmear;
      if( row > 0 ) BBSH_sums[row-1] += esmear; 
    }

    double EmaxSH=0.0;
    int imax_SH=-1;
    for( int isum_SH=0; isum_SH<26; isum_SH++ ){
      if( BBSH_sums[isum_SH] > EmaxSH ){
	EmaxSH = BBSH_sums[isum_SH];
	imax_SH = isum_SH;
      }
    }

    if( imax_SH >= 0 ){
      hedep_BBSH->Fill( EmaxSH, weight );
    }
    
    bool trigger = false;
    
    if( imax_SH >= 0 && EmaxSH >= thresh_SH ){ //potentially trigger, require that none of the nearest-neighbor PS groups is above threshold:
      trigger = true;

      int ipsmin,ipsmax;
      if( imax_SH < 19 ){
	ipsmin=imax_SH-1;
	ipsmax=imax_SH+1;
      } else {
	ipsmin=imax_SH-2;
	ipsmax=imax_SH;
      }

      //crude approximation to using the PS as a veto:
      for( int iPS=ipsmin; iPS<=ipsmax; iPS++ ){
	if( iPS>=0&&iPS<25 && BBPS_sums[iPS] >= thresh_PS ){
	  trigger = false;
	}
      }
    }

    if( trigger ) {
      rate_BB_pion += weight;
      hedep_BBSH_anticutPS->Fill( EmaxSH, weight );
    } else if( imax_SH >= 0 ){
      hedep_BBSH_cutPS->Fill( EmaxSH, weight );
    }

    double EmaxPS=0.0;
    int imax_PS=-1;
    for( int isum_PS=0; isum_PS<25; isum_PS++ ){
      if( BBPS_sums[isum_PS] > EmaxPS ) {
	EmaxPS = BBPS_sums[isum_PS];
	imax_PS=isum_PS;
      }
    }

    if( imax_PS>=0 ){
      hedep_BBPS->Fill( EmaxPS, weight );
      if( EmaxSH >= thresh_SH ){
	hedep_BBPS_cutSH->Fill( EmaxPS, weight );
      } else {
	hedep_BBPS_anticutSH->Fill( EmaxPS, weight );
      }
    }
    
      
      

    // if( EPS >= thresh_PS ){
    //   hedep_BBSH_cutPS->Fill( ESH, weight );
    // } else {
    //   hedep_BBSH_anticutPS->Fill( ESH, weight );
    // }

    // if( ESH >= thresh_SH ){
    //   hedep_BBPS_cutSH->Fill( EPS, weight );
    // } else {
    //   hedep_BBPS_anticutSH->Fill( EPS, weight );
    // }

    // if( ESH >= thresh_SH && EPS <= thresh_PS ){ //pion logic:
    //   rate_BB_pion += weight;
    // }

    // if( ESH + EPS >= thresh_SH ){ //electron logic
    //   rate_BB_electron += weight;
    // }

    if( T->Harm_HCalScint_det_esum >= thresh_HCAL ){
      rate_HCAL += weight;
    }
    
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

  cout << "BigBite singles rate with pion logic = " << rate_BB_pion << endl;
  cout << "BigBite singles rate with electron logic = " << rate_BB_electron << endl;
  cout << "HCAL singles rate = " << rate_HCAL << endl;
  
  fout->Write();

}
