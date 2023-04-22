#include "gen_tree.C"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "G4SBSRunData.hh"
#include "TChainElement.h"
#include "TObjArray.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include <iostream>
#include <set>

const double PI = TMath::Pi();

const double mu_n = -1.913;

void gen_Aphys_plots( const char *inputfilename, const char *outputfilename, double thetapoldeg = 118.0, double phipoldeg=180.0, double Wmin=0.8, double Wmax=1.15 ){
  double thetapol = thetapoldeg*PI/180.0;
  double phipol   = phipoldeg*PI/180.0;

  double PBPT = 0.85*0.47;
  
  TVector3 TargetPol_Unit( sin(thetapol)*cos(phipol),
			   sin(thetapol)*sin(phipol),
			   cos(thetapol) );
  
  TChain *C = new TChain("T");

  C->Add(inputfilename);

  TFile *fout = new TFile(outputfilename,"RECREATE");

  TObjArray *FileList = C->GetListOfFiles();
  TIter next(FileList);

  TChainElement *chEl = 0;

  set<TString> bad_file_list;

  long ntries=0;
  long ngen=0;
  
  G4SBSRunData *rd;

  map<TString,double> BBang_file;
  map<TString,double> SBSang_file;
  map<TString,double> Ebeam_file;

  double BBtheta,SBStheta,Ebeam;

  bool onegoodrun = false;
  
  while( (chEl=(TChainElement*)next() )){
    TFile newfile(chEl->GetTitle());
    newfile.GetObject("run_data",rd);
    if( rd ){

      onegoodrun = true;
      
      ngen += rd->fNthrown;
      ntries += rd->fNtries;
      // ngen_file[chEl->GetTitle()] = rd->fNthrown;
      // ntries_file[chEl->GetTitle()] = rd->fNtries;
      BBang_file[chEl->GetTitle()] = rd->fBBtheta;
      SBSang_file[chEl->GetTitle()] = rd->fSBStheta;
      Ebeam_file[chEl->GetTitle()] = rd->fBeamE;

      BBtheta = BBang_file[chEl->GetTitle()];
      SBStheta = SBSang_file[chEl->GetTitle()];
      Ebeam = Ebeam_file[chEl->GetTitle()];
    } else {
      bad_file_list.insert( chEl->GetTitle());
    }
  }

  if( !onegoodrun ){
    cout << "failed to extract run metadata for any runs, exiting..." << endl;
  }

  fout->cd();
  
  TH1D *hAphys = new TH1D("hAphys","",200,-1.05,1.05);
  TProfile *hAphys_Q2 = new TProfile("hAphys_Q2","",200,0.0,4.0);
  TProfile *hAphys_W = new TProfile("hAphys_W","",200,0.0,1.9);
  TProfile *hAphys_phie = new TProfile("hAphys_phie","",200,-90.0,90.0);
  TH1D *hQ2 = new TH1D("hQ2","",200,0.0,4.0);
  TH1D *hepsilon = new TH1D("hepsilon","",200,-0.05,1.05);
  TH1D *hW  = new TH1D("hW","",200,0.0,1.9);
  TH1D *hAperp = new TProfile("hAperp","",200,-1.05,1.05);
  TH1D *hApar  = new TProfile("hApar","",200,-1.05,1.05);
  TH1D *hRatio = new TProfile("hRatio","",200,-1.05,1.05);
  TH1D *hPx = new TH1D("hPx","",200,-1.05,1.05);
  TH1D *hPz = new TH1D("hPz","",200,-1.05,1.05);
  TProfile *hAperp_Q2 = new TProfile("hAperp_Q2","",200,0.0,4.0);
  TProfile *hApar_Q2  = new TProfile("hApar_Q2","",200,0.0,4.0);
  TProfile *hRatio_Q2 = new TProfile("hRatio_Q2","",200,0.0,4.0);

  TProfile *hAperp_W = new TProfile("hAperp_W","",200,0.0,1.9);
  TProfile *hApar_W  = new TProfile("hApar_W","",200,0.0,1.9);
  TProfile *hRatio_W = new TProfile("hRatio_W","",200,0.0,1.9);

  TProfile *hAperp_phie = new TProfile("hAperp_phie","",200,-90.0,90.0);
  TProfile *hApar_phie  = new TProfile("hApar_phie","",200,-90.0,90.0);
  TProfile *hRatio_phie = new TProfile("hRatio_phie","",200,-90.0,90.0);

  TH1D *hFFratio_plus = new TH1D("hFFratio_plus","",200,-2.0,2.0);
  TH1D *hFFratio_minus = new TH1D("hFFratio_minus","",200,-2.0,2.0);
  
  TProfile *hFFratio_plus_Q2 = new TProfile("hFFratio_plus_Q2","",200,0.0,4.0);
  TProfile *hFFratio_plus_W  = new TProfile("hFFratio_plus_W","",200,0.0,1.9);
  TProfile *hFFratio_plus_phie = new TProfile("hFFratio_plus_phie","",200,-90.0,90.0);

  TProfile *hFFratio_minus_Q2 = new TProfile("hFFratio_minus_Q2","",200,0.0,4.0);
  TProfile *hFFratio_minus_W  = new TProfile("hFFratio_minus_W","",200,0.0,1.9);
  TProfile *hFFratio_minus_phie = new TProfile("hFFratio_minus_phie","",200,-90.0,90.0);
  
  TH2D *hW_Q2 = new TH2D("hW_Q2","",200,0.0,1.9,200,0.0,4.0);
  TH2D *hepsilon_Q2 = new TH2D("hepsilon_Q2","",200,0.0,4.0,200,-0.05,1.05);
  
  gen_tree *T = new gen_tree(C);

  long nevent=0;

  int treenum = 0;
  int oldtreenum = -1;
  
  while( C->GetEntry( nevent++ ) ){
    if( nevent%1000 == 0 ) cout << nevent << endl;

    treenum = C->GetTreeNumber();
    
    if( treenum != oldtreenum ){

      oldtreenum = treenum;
      
      TString fname = C->GetFile()->GetName();
    
      if( bad_file_list.find( fname ) == bad_file_list.end() ){
      
        Ebeam = Ebeam_file[fname];
	BBtheta = BBang_file[fname];
	SBStheta = SBSang_file[fname];

      } 
    }

    TVector3 SBS_zaxis( -sin(SBStheta),0,cos(SBStheta) );
    TVector3 SBS_xaxis( 0,-1,0 );
    TVector3 SBS_yaxis = SBS_zaxis.Cross(SBS_xaxis).Unit();

    double Q2 = T->ev_Q2;
    double W  = sqrt(T->ev_W2);
    TVector3 pe( T->ev_epx, T->ev_epy, T->ev_epz );
    TVector3 pN( T->ev_npx, T->ev_npy, T->ev_npz );

    double M = 0.938272;
    if( T->ev_fnucl == 0 ) M = 0.939565;
    
    TLorentzVector Pe( pe, T->ev_ep );
    TLorentzVector PN( pN, sqrt(pow(T->ev_np,2)+pow(M,2)) );

    TLorentzVector Pbeam( 0, 0, Ebeam, Ebeam );
    TLorentzVector q = Pbeam - Pe;
    
    double tau = Q2/(4.*pow(M,2));
    double etheta = T->ev_th;

    double epsilon = pow(1.+2.0*(1.0+tau)*pow(tan(etheta/2.),2),-1);
    //Compute theta* and phi*:

    TVector3 qvect = q.Vect();
    
    TVector3 qhat = qvect.Unit();

    TVector3 yhat = qvect.Cross( Pbeam.Vect() ).Unit();
    TVector3 xhat = yhat.Cross(qhat).Unit();

    double Px = TargetPol_Unit.Dot( xhat );
    double Pz = TargetPol_Unit.Dot( qhat );

    double Aphys = ( T->ev_Aperp * Px + T->ev_Apar * Pz );

    double Ratio = T->ev_Aperp / T->ev_Apar;

    //Aphys = -1/(1+epsilon/tau*r^2)*(sqrt(2eps*(1-eps)/tau)*r*Px + sqrt(1-eps^2) Pz);
    //-Aphys*(1+epsilon/tau*r^2) = sqrt(2eps*(1-eps)/tau)*Px*r + sqrt(1-eps^2)*Pz)
    // 0 = [Aphys+sqrt(1-eps^2)*Pz] + [sqrt(2eps*(1-eps)/tau)*Px]*r + Aphys*epsilon/tau*r^2
    double C = Aphys+sqrt(1.0-pow(epsilon,2))*Pz; //A in quadratic formula
    double B = sqrt(2.0*epsilon*(1.0-epsilon)/tau)*Px; //B in quadratic formula
    double A = epsilon/tau*Aphys;
    double FFratio_plus = (-B + sqrt(pow(B,2)-4.0*A*C))/(2.0*A);
    double FFratio_minus = (-B - sqrt(pow(B,2)-4.0*A*C))/(2.0*A);
    
    if( T->Earm_BBGEM_Track_ntracks == 1 && (*(T->Earm_BBGEM_Track_MID))[0] == 0 && T->ev_fnucl == 0 && W > Wmin && W < Wmax ){
      hAphys->Fill( Aphys, T->ev_rate );
      hAphys_Q2->Fill( Q2, Aphys, T->ev_rate );
      hAphys_W->Fill( W, Aphys, T->ev_rate );
      hAphys_phie->Fill( T->ev_ph * 180.0/PI , Aphys, T->ev_rate );
      hQ2->Fill( Q2, T->ev_rate );
      hW->Fill( W, T->ev_rate );
      hepsilon->Fill( epsilon, T->ev_rate );
      hAperp->Fill( T->ev_Aperp, T->ev_rate );
      hApar->Fill( T->ev_Apar, T->ev_rate );
      hRatio->Fill( Ratio, T->ev_rate );
      hAperp_Q2->Fill(Q2, T->ev_Aperp, T->ev_rate );
      hAperp_W->Fill(W,  T->ev_Aperp, T->ev_rate );
      hAperp_phie->Fill( T->ev_ph * 180.0/PI, T->ev_Aperp, T->ev_rate );

      hApar_Q2->Fill(Q2, T->ev_Apar, T->ev_rate );
      hApar_W->Fill(W,  T->ev_Apar, T->ev_rate );
      hApar_phie->Fill( T->ev_ph * 180.0/PI, T->ev_Apar, T->ev_rate );

      hRatio_Q2->Fill(Q2, Ratio, T->ev_rate );
      hRatio_W->Fill(W,  Ratio, T->ev_rate );
      hRatio_phie->Fill( T->ev_ph * 180.0/PI, Ratio, T->ev_rate );

      hW_Q2->Fill( Q2, W, T->ev_rate );
      hepsilon_Q2->Fill( Q2, epsilon, T->ev_rate );

      hFFratio_plus_Q2->Fill( Q2, FFratio_plus, T->ev_rate );
      hFFratio_plus_W->Fill( W, FFratio_plus, T->ev_rate );
      hFFratio_plus_phie->Fill( T->ev_ph*180.0/PI, FFratio_plus, T->ev_rate );

      hFFratio_minus_Q2->Fill( Q2, FFratio_minus, T->ev_rate );
      hFFratio_minus_W->Fill( W, FFratio_minus, T->ev_rate );
      hFFratio_minus_phie->Fill( T->ev_ph*180.0/PI, FFratio_minus, T->ev_rate );

      hFFratio_plus->Fill( FFratio_plus, T->ev_rate );
      hFFratio_minus->Fill( FFratio_minus, T->ev_rate );
      
      hPx->Fill(Px,T->ev_rate);
      hPz->Fill(Pz,T->ev_rate);
    }     
  }

  double Pxavg = hPx->GetMean();
  double Pzavg = hPz->GetMean();
  double Q2avg = hQ2->GetMean();
  double epsavg = hepsilon->GetMean();
  double Aphysavg = hAphys->GetMean();
  double tauavg = Q2avg/(4.*pow(0.939565,2));
  
  double Aavg = epsavg/tauavg * Aphysavg;
  double Bavg = sqrt(2.0*epsavg*(1.0-epsavg)/tauavg)*Pxavg;
  double Cavg = Aphysavg+sqrt(1.0-pow(epsavg,2))*Pzavg;

  double FFratio_minus_avg = (-Bavg - sqrt(pow(Bavg,2)-4.0*Aavg*Cavg))/(2.0*Aavg);

  cout << "FFratio minus solution from acceptance-average Aphys, Q2, epsilon, Px, Pz = " << FFratio_minus_avg << endl;
  cout << "Acceptance-average of event-by-event FFratio minus solution = " << hFFratio_minus->GetMean() << endl;

  cout << "FFratio(acceptance-average kinematics/Aphys)/acceptance-average FF ratio - 1 = "
       << FFratio_minus_avg/hFFratio_minus->GetMean() - 1.0 << endl;
  
  fout->Write();
}
