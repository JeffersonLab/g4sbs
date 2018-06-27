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

void gen_Aphys_plots( const char *inputfilename, const char *outputfilename, double thetapoldeg = 118.0, double phipoldeg=180.0 ){
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
  TProfile *hAphys_W2 = new TProfile("hAphys_W2","",200,0.0,1.9);
  TProfile *hAphys_phie = new TProfile("hAphys_phie","",200,-90.0,90.0);
  
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

    //    double epsilon = pow(1.+2.0*(1.0+tau)*pow(tan(etheta/2.),2),-1);
    //Compute theta* and phi*:

    TVector3 qvect = q.Vect();
    
    TVector3 qhat = qvect.Unit();

    TVector3 yhat = qvect.Cross( Pbeam.Vect() ).Unit();
    TVector3 xhat = yhat.Cross(qhat).Unit();

    double Px = TargetPol_Unit.Dot( xhat );
    double Pz = TargetPol_Unit.Dot( qhat );

    double Aphys = ( T->ev_Aperp * Px + T->ev_Apar * Pz );

    if( T->Earm_BBGEM_Track_ntracks == 1 && (*(T->Earm_BBGEM_Track_MID))[0] == 0 && T->ev_fnucl == 0){
      hAphys->Fill( Aphys, T->ev_rate );
      hAphys_Q2->Fill( T->ev_Q2, Aphys, T->ev_rate );
      hAphys_W2->Fill( T->ev_W2, Aphys, T->ev_rate );
      hAphys_phie->Fill( T->ev_ph * 180.0/PI , Aphys, T->ev_rate );
    }     
  }

  fout->Write();
}
