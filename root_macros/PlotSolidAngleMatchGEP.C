#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include "TCanvas.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "G4SBSRunData.hh"
#include "TMath.h"

#include "TVector3.h"
#include "gep_tree_with_spin.C"

#include <iostream>
#include <fstream>

void PlotSolidAngleMatchGEP( const char *rootfilename, const char *outfilename="temp.root", double thresh_ECAL=1.5, double thresh_HCAL=0.07 ){
  TFile *fout = new TFile(outfilename,"RECREATE");

  //We should be able to do this without a config file (read meta-data from root file)

  TChain *C = new TChain("T");
  C->Add(rootfilename);

  gep_tree_with_spin *T = new gep_tree_with_spin(C);

  //Histograms!

  TH1D *hE_ECAL = new TH1D("hE_ECAL","All events; ECAL energy dep (GeV);",300,0,6.0);
  TH1D *hE_ECAL_pcut = new TH1D("hE_ECAL_pcut", "Proton arm cut; ECAL energy dep (GeV);",300,0.0,6.0);
  TH1D *hE_ECAL_panticut = new TH1D("hE_ECAL_panticut", "Proton arm anticut; ECAL energy dep (GeV)", 300,0,6);

  TH1D *hQ2 = new TH1D("hQ2", "All events; Q^{2} (GeV^{2});",300,0,15.0);
  TH1D *hQ2cutp = new TH1D("hQ2cutp", "Proton arm cuts; Q^{2} (GeV^{2});",300,0,15.0);
  TH1D *hQ2cute = new TH1D("hQ2cute", "Electron arm cuts; Q^{2} (GeV^{2});",300,0,15.0);
  TH1D *hQ2cutep = new TH1D("hQ2cutep", "E+P cuts; Q^{2} (GeV^{2});",300,0,15.0);

  TH2D *hECALxy = new TH2D("hECALxy", "All events; ECAL y projected (m); ECAL x projected (m)", 150,-0.75,0.75, 300, -1.5,1.5);
  TH2D *hECALxy_cutp = new TH2D("hECALxy_cutp", "Proton arm cuts; ECAL y projected (m); ECAL x projected (m)", 150,-0.75,0.75, 300, -1.5,1.5);
  TH2D *hECALxy_acutp = new TH2D("hECALxy_acutp", "Proton arm anticuts; ECAL y projected (m); ECAL x projected (m)", 150,-0.75,0.75, 300, -1.5,1.5);
  TH2D *hECALxy_cute = new TH2D("hECALxy_cute", "Electron arm cuts; ECAL y projected (m); ECAL x projected (m)", 150,-0.75,0.75, 300, -1.5,1.5);
  TH2D *hECALxy_acute = new TH2D("hECALxy_acute", "Electron arm anticuts; ECAL y projected (m); ECAL x projected (m)", 150,-0.75,0.75, 300, -1.5,1.5);
  TH2D *hECALxy_cutep = new TH2D("hECALxy_cutep", "E+P arm cuts; ECAL y projected (m); ECAL x projected (m)", 150,-0.75,0.75, 300, -1.5,1.5);
  TH2D *hECALxy_cute_acutp = new TH2D("hECALxy_cute_acutp", "E+!P; ECAL y projected (m); ECAL x projected (m)", 150,-0.75,0.75, 300, -1.5,1.5);
  TH2D *hECALxy_cutp_acute = new TH2D("hECALxy_cutp_acute", "P+!E; ECAL y projected (m); ECAL x projected (m)", 150,-0.75,0.75, 300, -1.5,1.5);

  TH1D *hECALx = new TH1D("hECALx", "All events; ECAL x projected (m);", 300,-1.5,1.5);
  TH1D *hECALx_cute = new TH1D("hECALx_cute", "E arm cut; ECAL x projected (m);", 300,-1.5,1.5);
  TH1D *hECALx_acute = new TH1D("hECALx_acute", "E arm anticut; ECAL x projected (m);", 300,-1.5,1.5);
  TH1D *hECALx_cutp = new TH1D("hECALx_cutp", "P arm cut; ECAL x projected (m);", 300,-1.5,1.5);
  TH1D *hECALx_acutp = new TH1D("hECALx_acutp", "P arm anticut; ECAL x projected (m);", 300,-1.5,1.5);
  TH1D *hECALx_cutep = new TH1D("hECALx_cutep", "E+P arm cuts; ECAL x projected (m);", 300,-1.5,1.5);
 
  TH1D *hECALy = new TH1D("hECALy", "All events; ECAL y projected (m);", 150,-0.75,0.75);
  TH1D *hECALy_cute = new TH1D("hECALy_cute", "E arm cut; ECAL y projected (m);", 150,-0.75,0.75);
  TH1D *hECALy_acute = new TH1D("hECALy_acute", "E arm anticut; ECAL y projected (m);", 150,-0.75,0.75);
  TH1D *hECALy_cutp = new TH1D("hECALy_cutp", "P arm cut; ECAL y projected (m);", 150,-0.75,0.75);
  TH1D *hECALy_acutp = new TH1D("hECALy_acutp", "P arm anticut; ECAL y projected (m);", 150,-0.75,0.75);
  TH1D *hECALy_cutep = new TH1D("hECALy_cutep", "E+P arm cuts; ECAL y projected (m);", 150,-0.75,0.75);
  
  //Now we also need 1D and 2D histograms for theta, phi, th vs phi etc for various cuts to set event generation limits

  TH1D *hetheta = new TH1D("hetheta", "All events; #theta_{e} (deg)", 500, 0, 50);
  TH1D *hetheta_cute = new TH1D("hetheta_cute", "ECAL cut; #theta_{e} (deg)", 500, 0, 50);
  TH1D *hetheta_anticute = new TH1D("hetheta_anticute", "ECAL anticut; #theta_{e} (deg)", 500, 0, 50);
  TH1D *hetheta_cutp = new TH1D("hetheta_cutp", "P arm cut; #theta_{e} (deg)", 500, 0, 50);
  TH1D *hetheta_anticutp = new TH1D("hetheta_anticutp", "P arm anticut; #theta_{e} (deg)", 500, 0, 50);
  TH1D *hetheta_cutep = new TH1D("hetheta_cutep", "E+P arm cuts; #theta_{e} (deg)", 500, 0, 50);
  
  TH1D *hephi = new TH1D("hephi", "All events; #phi_{e} (deg)", 450, -45, 45);
  TH1D *hephi_cute = new TH1D("hephi_cute", "ECAL cut; #phi_{e} (deg)", 450, -45, 45);
  TH1D *hephi_anticute = new TH1D("hephi_anticute", "ECAL anticut; #phi_{e} (deg)", 450, -45, 45);
  TH1D *hephi_cutp = new TH1D("hephi_cutp", "P arm cut; #phi_{e} (deg)", 450, -45, 45);
  TH1D *hephi_anticutp = new TH1D("hephi_anticutp", "P arm anticut; #phi_{e} (deg)", 450, -45, 45);
  TH1D *hephi_cutep = new TH1D("hephi_cutep", "E+P arm cuts; #phi_{e} (deg)", 450, -45, 45);
  
  TH2D *hethph = new TH2D("hethph", "All events; #theta_{e} (deg); #phi_{e} (deg)", 250, 0, 50, 225, -45, 45 );
  TH2D *hethph_cute = new TH2D("hethph_cute", "ECAL cut; #theta_{e} (deg); #phi_{e} (deg)", 250, 0, 50, 225, -45, 45 );
  TH2D *hethph_cutp = new TH2D("hethph_cutp", "P arm cuts; #theta_{e} (deg); #phi_{e} (deg)", 250, 0, 50, 225, -45, 45 );
  TH2D *hethph_cutep = new TH2D("hethph_cutep", "E+P arm cuts; #theta_{e} (deg); #phi_{e} (deg)", 250, 0, 50, 225, -45, 45 );
  TH2D *hethph_anticute = new TH2D("hethph_anticute", "ECAL anticut; #theta_{e} (deg); #phi_{e} (deg)", 250, 0, 50, 225, -45, 45 );
  TH2D *hethph_anticutp = new TH2D("hethph_anticutp", "P arm anticuts; #theta_{e} (deg); #phi_{e} (deg)", 250, 0, 50, 225, -45, 45 );

  TH2D *hethph_cute_acutp = new TH2D("hethph_cute_acutp", "ECAL cut; #theta_{e} (deg); #phi_{e} (deg)", 250, 0, 50, 225, -45, 45 );
  TH2D *hethph_cutp_acute = new TH2D("hethph_cutp_acute", "P arm cuts; #theta_{e} (deg); #phi_{e} (deg)", 250, 0, 50, 225, -45, 45 );
  
  double ECALdist=5.0, ECALtheta=36.1*TMath::DegToRad();
  double SBSdist=1.6, SBStheta=28.5;
  double HCALdist=10.0;
  double Ebeam=4.359;

  long nevent=0;

  int treenum=-1, oldtreenum=-1;
  
  while(C->GetEntry(nevent++) ){
    if( nevent%1000 == 0 ) cout << nevent << endl;
    
    treenum = C->GetTreeNumber();
    if( treenum != oldtreenum ){
      oldtreenum = treenum;
      G4SBSRunData *rdtemp;
      TFile *ftemp = C->GetFile();

      ftemp->GetObject( "run_data", rdtemp );
      if( rdtemp ){ //read relevant parameters from file:
	ECALdist = rdtemp->fBBdist;
	ECALtheta = rdtemp->fBBtheta;
	SBSdist = rdtemp->fSBSdist;
	SBStheta = rdtemp->fSBStheta;
	HCALdist = rdtemp->fHCALdist;
	Ebeam = rdtemp->fBeamE;
      }
    }

    bool FTtrack = false, FPPtrack = false, HCALcut = false, ECALcut = false;
    HCALcut = T->Harm_HCalScint_det_esum>thresh_HCAL;
    ECALcut = T->Earm_ECalTF1_det_esum>thresh_ECAL;

    for( int itr=0; itr<T->Harm_FT_Track_ntracks; itr++ ){
      if( T->Harm_FT_Track_MID->at(itr) == 0 ) FTtrack = true;
    }
    for( int itr=0; itr<T->Harm_FPP1_Track_ntracks; itr++ ){
      if( T->Harm_FPP1_Track_MID->at(itr) == 0 ) FPPtrack = true;
    }

    hE_ECAL->Fill( T->Earm_ECalTF1_det_esum );
    bool pcut = FTtrack && FPPtrack && HCALcut;
    if( pcut ) {
      hE_ECAL_pcut->Fill( T->Earm_ECalTF1_det_esum );
    } else {
      hE_ECAL_panticut->Fill( T->Earm_ECalTF1_det_esum );
    }

    hetheta->Fill( T->ev_th*TMath::RadToDeg() );
    hephi->Fill( T->ev_ph*TMath::RadToDeg() );
    hethph->Fill( T->ev_th*TMath::RadToDeg(), T->ev_ph*TMath::RadToDeg() );

    if( ECALcut ){
      hetheta_cute->Fill( T->ev_th*TMath::RadToDeg() );
      hephi_cute->Fill( T->ev_ph*TMath::RadToDeg() );
      hethph_cute->Fill( T->ev_th*TMath::RadToDeg(), T->ev_ph*TMath::RadToDeg() );
      if( pcut ){
	hetheta_cutep->Fill( T->ev_th*TMath::RadToDeg() );
	hephi_cutep->Fill( T->ev_ph*TMath::RadToDeg() );
	hethph_cutep->Fill( T->ev_th*TMath::RadToDeg(), T->ev_ph*TMath::RadToDeg() );
      } else {
	hethph_cute_acutp->Fill( T->ev_th*TMath::RadToDeg(), T->ev_ph*TMath::RadToDeg() );
      }
      
    } else {
      hetheta_anticute->Fill( T->ev_th*TMath::RadToDeg() );
      hephi_anticute->Fill( T->ev_ph*TMath::RadToDeg() );
      hethph_anticute->Fill( T->ev_th*TMath::RadToDeg(), T->ev_ph*TMath::RadToDeg() );
      if( pcut ){
	hethph_cutp_acute->Fill( T->ev_th*TMath::RadToDeg(), T->ev_ph*TMath::RadToDeg() );
      }
    }

    if( pcut ){
      hetheta_cutp->Fill( T->ev_th*TMath::RadToDeg() );
      hephi_cutp->Fill( T->ev_ph*TMath::RadToDeg() );
      hethph_cutp->Fill( T->ev_th*TMath::RadToDeg(), T->ev_ph*TMath::RadToDeg() );
    } else {
      hetheta_anticutp->Fill( T->ev_th*TMath::RadToDeg() );
      hephi_anticutp->Fill( T->ev_ph*TMath::RadToDeg() );
      hethph_anticutp->Fill( T->ev_th*TMath::RadToDeg(), T->ev_ph*TMath::RadToDeg() );
    }
    hQ2->Fill( T->ev_Q2 );
    if( pcut ){
      hQ2cutp->Fill( T->ev_Q2 );
    }
    if( ECALcut ) hQ2cute->Fill( T->ev_Q2 );
    if( ECALcut && pcut ) hQ2cutep->Fill( T->ev_Q2 );

    TVector3 ehat( sin( T->ev_th ) * cos( T->ev_ph ),
		   sin( T->ev_th ) * sin( T->ev_ph ),
		   cos( T->ev_th ) );

    TVector3 ECALzaxis(sin(ECALtheta), 0, cos(ECALtheta) );
    TVector3 ECALxaxis(cos(ECALtheta), 0, -sin(ECALtheta) );
    TVector3 ECALyaxis(0,1,0);

    TVector3 ECAL_origin = ECALdist * ECALzaxis;

    TVector3 vertex(T->ev_vx, T->ev_vy, T->ev_vz );

    double sint = (ECAL_origin - vertex).Dot( ECALzaxis ) / (ehat.Dot(ECALzaxis) );

    TVector3 intersect = vertex + sint * ehat;
    
    double xECAL = (intersect - ECAL_origin).Dot( ECALxaxis );
    double yECAL = (intersect - ECAL_origin).Dot( ECALyaxis );

    hECALx->Fill( yECAL );
    hECALy->Fill( xECAL );
    hECALxy->Fill( xECAL, yECAL );
    if( pcut ){
      hECALxy_cutp->Fill( xECAL, yECAL );
      hECALx_cutp->Fill( yECAL );
      hECALy_cutp->Fill( xECAL );
    } else {
      hECALxy_acutp->Fill( xECAL, yECAL );
      hECALx_acutp->Fill( yECAL );
      hECALy_acutp->Fill( xECAL );
    }

    if( ECALcut ){
      hECALxy_cute->Fill( xECAL, yECAL );
      hECALx_cute->Fill( yECAL );
      hECALy_cute->Fill( xECAL );
    } else {
      hECALxy_acute->Fill( xECAL, yECAL );
      hECALx_acute->Fill( yECAL );
      hECALy_acute->Fill( xECAL );
    };

    if( ECALcut && pcut ){
      hECALxy_cutep->Fill( xECAL, yECAL );
      hECALx_cutep->Fill( yECAL );
      hECALy_cutep->Fill( xECAL );
    }

    if( ECALcut && !pcut ){
      hECALxy_cute_acutp->Fill( xECAL, yECAL );
    }

    if( pcut && !ECALcut ){
      hECALxy_cutp_acute->Fill( xECAL, yECAL );
    }
      
    
  }

  fout->Write();
}
