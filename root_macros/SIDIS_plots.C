#include "sidis_tree.C"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TCut.h"
#include "G4SBSRunData.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TEventList.h"
#include "TMath.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TVector3.h"
#include <set>

void SIDIS_plots( const char *infilename, double ndays=40.0, const char *outfilename="SIDIS_plots_temp.root", int pi0flag=0 ){

  TFile *fout = new TFile(outfilename,"RECREATE" );
  
  TChain *C = new TChain("T");
  
  ifstream infile(infilename);

  //first read list of files:

  long ntries_total = 0;

  vector<G4SBSRunData *> RunData; //grab the run metadata object for each file:
  
  TString currentline;
  while( currentline.ReadLine( infile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());

      G4SBSRunData *rdtemp;
      TFile *ftemp = new TFile( currentline.Data(), "READ" );
      ftemp->GetObject( "run_data", rdtemp );
      
      RunData.push_back( rdtemp );
      
      ntries_total += rdtemp->fNtries;
    }
  }

  TCut cut = "";
  while( currentline.ReadLine( infile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      cut += currentline.Data();
    }
  }

  TEventList *elist = new TEventList("elist","list of events");

  C->Draw(">>elist", cut );
  
  sidis_tree *T = new sidis_tree( C );

  long nevent=0;

  double PI = TMath::Pi();

  fout->cd();
  //Start making histograms:
  TH1D *hxbj = new TH1D("hxbj","Bjorken x",100,0,1);
  TH1D *hQ2 = new TH1D("hQ2","Q^{2} (GeV^{2})",100,0,15);
  TH1D *hy = new TH1D("hy", "y = #nu/E", 100,0,1);
  TH1D *hz  = new TH1D("hz","z hadron", 100,0,1);
  TH1D *hpT = new TH1D("hpT","p_{T} hadron (GeV)",100,0,2);
  TH1D *hW = new TH1D("hW","W (GeV)",100,0,7);
  TH1D *hMX = new TH1D("hMX", "M_{X} (GeV)", 100, 0, 5);

  TH1D *hphih = new TH1D("hphih", "#phi_{h} (rad)", 90, -PI, PI );
  TH1D *hphiS = new TH1D("hphiS", "#phi_{S} (rad)", 90, -PI, PI );
  TH1D *hphiSiv = new TH1D("hphiSiv","#phi_{h}-#phi_{S}", 90, -PI, PI );
  TH1D *hphiColl = new TH1D("hphiColl", "#phi_{h}+#phi_{S}", 90, -PI, PI );

  
  TH1D *hsinphiSiv = new TH1D("hsinphiSiv","sin(#phi_{h}-#phi_{S})",100,-1,1);
  TH1D *hsinphiColl = new TH1D("hsinphiColl", "sin(#phi_{h}+#phi_{S})",100,-1,1);
  TH1D *hsinphiPretz = new TH1D("hsinphiPretz", "sin(3#phi_{h}-#phi_{S})",100,-1,1);
  
  TH2D *hQ2_vs_x = new TH2D("hQ2_vs_x", "Q^{2} vs x", 100,0,1,100,0,15);
  TH2D *hz_vs_x = new TH2D("hz_vs_x", "z vs x", 100, 0, 1, 100, 0, 1 );
  TH2D *hpT_vs_x = new TH2D("hpT_vs_x", "p_{T} vs x", 100, 0, 1, 100, 0, 2);
  TH2D *hsinphiSiv_vs_x = new TH2D("hsinphiSiv_vs_x","sin(#phi_{h}-#phi_{S}) vs. x", 100, 0, 1, 100, -1, 1 );
  TH2D *hsinphiColl_vs_x = new TH2D("hsinphiColl_vs_x","sin(#phi_{h}+#phi_{S}) vs. x", 100, 0, 1, 100, -1, 1 );

  TH2D *hphiSiv_vs_x = new TH2D("hphiSiv_vs_x","#phi_{h}-#phi_{S} vs x", 100, 0, 1, 100, -PI, PI );
  TH2D *hphiColl_vs_x = new TH2D("hphiColl_vs_x", "#phi_{h}+#phi_{S} vs x", 100, 0, 1, 100, -PI, PI );
  
  TH2D *hphih_vs_x = new TH2D("hphih_vs_x", "#phi_{h} vs x", 100, 0, 1, 100, -PI, PI );

  TH1D *hxbj_p = new TH1D("hxbj_p", "Bjorken x (struck proton)", 100, 0, 1 );
  TH1D *hxbj_n = new TH1D("hxbj_n", "Bjorken x (struck neutron)", 100, 0, 1 );
  
  while( C->GetEntry( elist->GetEntry( nevent++ ) ) ){
    double Ebeam, SBStheta, BBtheta, Lumi;
    int ifile = C->GetTreeNumber();

    Ebeam = RunData[ifile]->fBeamE;
    SBStheta = RunData[ifile]->fSBStheta;
    BBtheta = RunData[ifile]->fBBtheta;
    Lumi = RunData[ifile]->fLuminosity;

    //This will give a total yield:
    double weight = T->ev_sigma * Lumi/double(ntries_total);

    if( nevent % 1000 == 0 ){
      cout << "event, Ebeam, SBStheta, BBtheta, Lumi = " << nevent << ", " << Ebeam << ", " << SBStheta*57.3 << ", " << BBtheta*57.3 << ", " << Lumi << endl;
    }

    //Now implement detector cuts for filling histograms:
    //Hadron arm:
    // HCAL threshold > 0.1 GeV.
    // Good pion track in SBS GEMs
    // RICH hits > 3 (for pions)
    
    // Shower + preshower threshold > 0.8 GeV
    // Good electron track in BB
    // GRINCH hits > 4

    double HCALsum = T->Harm_HCalScint_det_esum;
    double BBCALsum = T->Earm_BBPSTF1_det_esum+T->Earm_BBSHTF1_det_esum;
    int GRINCH_nhits = T->Earm_GRINCH_hit_nhits;
    int RICH_nhits = T->Harm_RICH_hit_nhits;
    
    bool goodHarm = T->Harm_SBSGEM_Track_ntracks==1 && (*(T->Harm_SBSGEM_Track_MID))[0]==0 && (*(T->Harm_SBSGEM_Track_P))[0]/T->ev_np >= 0.95 && HCALsum >= 0.1 && RICH_nhits >= 3;

    if( pi0flag != 0 ){ //then we need some different analysis:
      //we need to loop over all HCAL hits and find the ones caused by the pi0 decay photons:
      //Idea is to require both photons to hit HCAL:
      goodHarm = false;

      vector<int> TIDs;
      set<int> unique_photon_TIDs;
      map<int,double> photon_xhit,photon_yhit,photon_zhit,photon_Esum;
      for( int ihit=0; ihit<T->Harm_HCalScint_hit_nhits; ihit++ ){
	int otridx = (*(T->Harm_HCalScint_hit_otridx))[ihit];
	int sdtridx = (*(T->Harm_HCalScint_hit_sdtridx))[ihit];
	int ptridx = (*(T->Harm_HCalScint_hit_ptridx))[ihit];
	int tid = (*(T->SDTrack_TID))[sdtridx];
	//This is one of the pi0 decay photons:
	if( (*(T->SDTrack_PID))[sdtridx] == 22 && (*(T->SDTrack_MID))[sdtridx]==2 ){
	  auto newtrack = unique_photon_TIDs.insert( tid );
	  if( newtrack.second ){
	    photon_xhit[tid] = (*(T->SDTrack_posx))[sdtridx];
	    photon_yhit[tid] = (*(T->SDTrack_posy))[sdtridx];
	    photon_zhit[tid] = (*(T->SDTrack_posz))[sdtridx];
	    photon_Esum[tid] = (*(T->Harm_HCalScint_hit_sumedep))[ihit];
	    TIDs.push_back( tid );
	  } else {
	    // photon_xhit[tid] += (*(T->SDTrack_posx))[sdtridx];
	    // photon_yhit[tid] += (*(T->SDTrack_posy))[sdtridx];
	    // photon_zhit[tid] += (*(T->SDTrack_posz))[sdtridx];
	    photon_Esum[tid] += (*(T->Harm_HCalScint_hit_sumedep))[ihit]; 
	  }
	}
      }

      if( unique_photon_TIDs.size() == 2 ){
	double E1 = photon_Esum[TIDs[0]];
	TVector3 pos1(photon_xhit[TIDs[0]],
		      photon_yhit[TIDs[0]],
		      photon_zhit[TIDs[0]] );
	double E2 = photon_Esum[TIDs[1]];
	TVector3 pos2(photon_xhit[TIDs[1]],
		      photon_yhit[TIDs[1]],
		      photon_zhit[TIDs[1]] );

	if( (E1 > 0.1 || E2 > 0.1) && (pos1-pos2).Mag() >= 0.3 &&
	    T->Harm_SBSGEM_Track_ntracks == 0 ){
	  goodHarm = true;
	}
      }
    }
    
    bool goodEarm = BBCALsum >= 0.8 && T->Earm_BBGEM_Track_ntracks == 1 && (*(T->Earm_BBGEM_Track_MID))[0]==0 && GRINCH_nhits > 4;
    
    if( goodEarm && goodHarm ){ //passed trigger threshold, track, and Cherenkov cuts:
      hxbj->Fill( T->ev_xbj, weight );
      hQ2->Fill( T->ev_Q2, weight );
      hy->Fill( 1.0 - T->ev_ep/Ebeam, weight );
      hz->Fill( T->ev_z, weight );
      hpT->Fill( T->ev_phperp, weight );
      hW->Fill( sqrt(T->ev_W2), weight );
      hMX->Fill( sqrt(T->ev_MX2), weight );
      hphih->Fill( T->ev_phih, weight );
      hphiS->Fill( T->ev_phiS, weight );

      double phiSiv = T->ev_phih - T->ev_phiS;
      if( phiSiv > PI ) phiSiv -= 2.0*PI;
      if( phiSiv < -PI ) phiSiv += 2.0*PI;

      double phiColl = T->ev_phih + T->ev_phiS;
      if( phiColl > PI ) phiColl -= 2.0*PI;
      if( phiColl < -PI ) phiColl += 2.0*PI;

      double phiPretz = 3.*T->ev_phih - T->ev_phiS;
      
      hphiSiv->Fill( phiSiv, weight );
      hphiColl->Fill( phiColl, weight );
      hsinphiSiv->Fill( sin(phiSiv), weight );
      hsinphiColl->Fill( sin(phiColl), weight );
      hsinphiPretz->Fill( sin(phiPretz), weight );
      
      hQ2_vs_x->Fill( T->ev_xbj, T->ev_Q2, weight );
      hz_vs_x->Fill( T->ev_xbj, T->ev_z, weight );
      hpT_vs_x->Fill( T->ev_xbj, T->ev_phperp, weight );
      hsinphiSiv_vs_x->Fill( T->ev_xbj, sin(phiSiv), weight );
      hsinphiColl_vs_x->Fill( T->ev_xbj, sin(phiColl), weight );
      hphih_vs_x->Fill( T->ev_xbj, T->ev_phih, weight );

      hphiSiv_vs_x->Fill( T->ev_xbj, phiSiv, weight );
      hphiColl_vs_x->Fill( T->ev_xbj, phiColl, weight );
      
      if( T->ev_nucl == 0 ){ //neutron
	hxbj_n->Fill( T->ev_xbj, weight );
      } else { //proton
	hxbj_p->Fill( T->ev_xbj, weight );
      }
    }
    
  }
  
  elist->Delete();
  fout->Write();
}
