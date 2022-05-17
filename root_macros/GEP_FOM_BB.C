#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "G4SBSRunData.hh"
#include "gep_BB_tree.C"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TEventList.h"
#include "TF2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set> 

const double Mp = 0.9382720813;
const double mu_p = 2.79284734462;

double KellyFunc(double *x, double *par){
  double Q2 = x[0];
  double tau = Q2/(4.*pow(Mp,2));

  return (1.0 + par[0]*tau)/(1.0 + par[1]*tau + par[2]*pow(tau,2) + par[3]*pow(tau,3)); //returns GEp or GMp/mu
}

//Analyzing power parametrization: x = p_T, y = p
TF2 *Ayfunc = new TF2("Ayfunc","([0]+[1]/y)*x*exp(-[2]*pow(x,2))",0.0,2.0, 1.0,15.0);
TF1 *GEPfunc = new TF1("GEPfunc",KellyFunc, 0.0,40.0,4);
TF1 *GMPfunc = new TF1("GMPfunc",KellyFunc, 0.0,40.0,4);

void GEP_FOM_quick_and_dirty(const char *configfilename, const char *outfilename="GEP_FOM.root"){
  ifstream configfile(configfilename);

  double geppar[4] = {-0.01,12.16,9.7,37.0};
  double gmppar[4] = {0.093,11.07,19.1,5.6};

  GEPfunc->SetParameters(geppar);
  GMPfunc->SetParameters(gmppar);
  
  double pTmin = 0.06;
  double pTmax = 1.2;

  double Aymax_pslope = 0.348;
  double Aymax_pintercept = 0.0342;
  double Ay_pT_b = (2.376*pow(0.189,-2) + 2.776*pow(0.212,-2) + 2.292*pow(0.145,-2))/(pow(0.189,-2)+pow(0.212,-2)+pow(0.145,-2));

  double beampol = 0.85;

  Ayfunc->SetParameter( 0, Aymax_pintercept * sqrt(2.0*exp(1.0)*Ay_pT_b) );
  Ayfunc->SetParameter( 1, Aymax_pslope * sqrt(2.0*exp(1.0)*Ay_pT_b) );
  Ayfunc->SetParameter( 2, Ay_pT_b );
  
  Ayfunc->Draw("cont4z");

  if( !configfile ) return;

  TString currentline;

  TChain *C = new TChain("T");

  TFile *ftemp;
  G4SBSRunData *rdtemp;

  double Ngen_total = 0.0;
  
  map<TString, double> Ebeam_file;
  map<TString, double> SBStheta_file;
  map<TString, double> SBStracker_pitch_file;
  map<TString, double> Ngen_file;
  map<TString, double> Norm_file; //normalization is defined such that rate = dsigma * Normalization
  //map<TString, double> EvtWeight_file; //should be norm_file * Ngen_file/Ngen_total, weight for individual events
  map<TString,double> Lumi_file;
  
  double Ebeam_default = 11.0;
  double Lumi = 6e38;
  
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){ //Each line is assumed to be one ROOT file:
    if( !currentline.BeginsWith("#") ){
      ftemp = new TFile(currentline,"READ");
      if( !ftemp->IsZombie() ){
	ftemp->GetObject("run_data",rdtemp);
	if( rdtemp != NULL ){
	  Ebeam_file[currentline] = rdtemp->fBeamE;
	  Ngen_file[currentline] = rdtemp->fNtries;
	  Norm_file[currentline] = rdtemp->fNormalization;
	  Ngen_total += Ngen_file[currentline];
	  Lumi_file[currentline] = rdtemp->fLuminosity;
	  
	  Ebeam_default = Ebeam_file[currentline];
	  SBStheta_file[currentline] = rdtemp->fSBStheta;
	  SBStracker_pitch_file[currentline] = rdtemp->fSBSTrackerPitch;
	  
	  C->Add(currentline);

	  cout << "added file " << currentline << " ngen = " << Ngen_file[currentline]
	       << " total ngen = " << Ngen_total << endl;
	}
      }
    }
  }

  //for( map<TString,double>::iterator it=Ngen_fi
  
  if( C->GetNtrees() == 0 ) return;

  TCut globalcut = "";
  
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }
  }

  TEventList *elist_temp = new TEventList("elist_temp");

  C->Draw(">>elist_temp",globalcut);

  gep_BB_tree *T = new gep_BB_tree(C);

  long nevent = 0;

  double NAy2_sum = 0.0;
  double Ngoodevent_sum = 0.0;

  double PT_sum = 0.0;
  double PL_sum = 0.0;

  double sinchi_sum = 0.0;
  double Q2_sum = 0.0;
  double epsilon_sum = 0.0;
  double kinfact_sum = 0.0;
  double FFratio_sum = 0.0;
  
  //read configuration parameters:
  double Ndays = 30.0;
  double effrecon = 0.7;
  configfile >> Ndays;
  configfile >> beampol;
  configfile >> effrecon;

  // int FPPoption=2;
  // configfile >> FPPoption;

  

  TFile *fout = new TFile(outfilename,"RECREATE");
  
  TH1D *hQ2 = new TH1D("hQ2","",150,0.0,15.0);
  TH1D *hepsilon = new TH1D("hepsilon","",250,0.0,1.0);
  TH1D *hetheta = new TH1D("hethetadeg","",250,0.0,90.0);
  TH1D *hEprime = new TH1D("hEprime","",250,0.0,11.0);
  TH1D *hptheta = new TH1D("hpthetadeg","",250,0.0,90.0);
  TH1D *hpp     = new TH1D("hpp","",250,0.0,11.0);
  TH1D *hpT_FPP1     = new TH1D("hpT_FPP1","",250,0.0,2.0);
  TH1D *hpT_FPP2     = new TH1D("hpT_FPP2","",250,0.0,2.0);
  
  int treenumber = -1, oldtreenumber = -1;

  double Ebeam_current = Ebeam_default;
  double SBStheta = 16.9*TMath::Pi()/180.0;
  double SBStracker_pitch = 5.0*TMath::Pi()/180.0;
  //double EventWeight_current = 1.0;

  double weight = 1.0;
  
  while( T->GetEntry( elist_temp->GetEntry( nevent++ ) ) ){
    if( (nevent-1) % 1000 == 0 ) cout << "Event " << nevent << endl;

    //double weight = T->ev_rate;
    
    int treenum = C->GetTreeNumber();
    if( treenum != oldtreenumber ){
      cout << "New tree found, new tree number = " << treenum
	   << " old tree number = " << oldtreenumber
	   << " event = " << nevent << endl;
      
      TString fname_temp = C->GetFile()->GetName();
      Ebeam_current = Ebeam_file[fname_temp];
      Lumi = Lumi_file[fname_temp];
      SBStheta = SBStheta_file[fname_temp];
      SBStracker_pitch = SBStracker_pitch_file[fname_temp];
      oldtreenumber = treenum;
    }

    weight = T->ev_sigma * T->ev_solang * Lumi / Ngen_total;
    
    TVector3 SBS_zaxis( -sin(SBStheta), 0, cos(SBStheta) );
    TVector3 SBS_xaxis(0, -1, 0 );
    TVector3 SBS_yaxis = SBS_zaxis.Cross(SBS_xaxis).Unit();

    TVector3 ppvect_global(T->ev_npx, T->ev_npy, T->ev_npz );
    TVector3 ppunit_global = ppvect_global.Unit();

    TVector3 ppunit_sbs( ppunit_global.Dot( SBS_xaxis ),
			 ppunit_global.Dot( SBS_yaxis ),
			 ppunit_global.Dot( SBS_zaxis ) );

    TVector3 SBS_FT_zaxis( -sin(SBStracker_pitch) ,0,cos(SBStracker_pitch) );
    TVector3 SBS_FT_yaxis( 0,1,0 );
    TVector3 SBS_FT_xaxis = SBS_FT_yaxis.Cross(SBS_FT_zaxis).Unit();

    // SBS_FT_zaxis.Print();
    // SBS_FT_xaxis.Print();
    
    //Get FT Track:
    int idx_FT_track = -1;
    for( int itrack=0; itrack<T->Harm_FT_Track_ntracks; itrack++ ){
      int MID = (*(T->Harm_FT_Track_MID))[itrack];
      if( MID == 0 ) {
	idx_FT_track = itrack;
	break;
      }
    }

    //Get FPP1 Track:
    int idx_FPP1_track = -1;

    // for( int itrack=0; itrack<T->Harm_FPP1_Track_ntracks; itrack++ ){
    //   int MID = (*(T->Harm_FPP1_Track_MID))[itrack];
    //   if( MID == 0 || T->Harm_FPP1_Track_ntracks == 1 ) {
    // 	idx_FPP1_track = itrack;
    // 	break;
    //   }
    // }
    if( T->Harm_FPP1_Track_ntracks == 1 ) idx_FPP1_track = 0;

    //Get FPP2 Track:
    int idx_FPP2_track = -1;

    // for( int itrack=0; itrack<T->Harm_FPP2_Track_ntracks; itrack++ ){
    //   int MID = (*(T->Harm_FPP2_Track_MID))[itrack];
    //   if( MID == 0 || T->Harm_FPP1_Track_ntracks == 1 ) {
    // 	idx_FPP2_track = itrack;
    // 	break;
    //   }
    // }
    //if( T->Harm_FPP2_Track_ntracks == 1 ) idx_FPP2_track = 0;

    if( idx_FT_track >= 0 ){
      TVector3 FT_track( (*(T->Harm_FT_Track_Xp))[idx_FT_track],
			 (*(T->Harm_FT_Track_Yp))[idx_FT_track],
			 1.0 );

      FT_track = FT_track.Unit();

      TVector3 FT_track_SBS = FT_track.X() * SBS_FT_xaxis +
	FT_track.Y() * SBS_FT_yaxis +
	FT_track.Z() * SBS_FT_zaxis;

      double thetabend = acos( FT_track_SBS.Dot( ppunit_sbs ) );
      
      
      double pp_FT = (*(T->Harm_FT_Track_P))[idx_FT_track];
      double etheta = T->ev_th;
      double Eprime = T->ev_ep;
      double ptheta = T->ev_nth;
      double pp = T->ev_np;

      // cout << "SBS theta = " << SBStheta * 57.3 << endl;
      // cout << "FT track in SBS target coords = " << endl;
      //FT_track_SBS.Print();
      // cout << "Target track in SBS target coords = " << endl;
      //ppunit_sbs.Print();
      
      //      cout << "pp, thetabend, pp*thetabend, BdL = " << pp << ", " << thetabend * 57.3 << ", " << pp*thetabend
      //	   << ", " << pp*thetabend/0.3 << endl;

      double gamma = sqrt(1.+pow(pp/Mp,2));
      double chi = gamma*(mu_p - 1.0)*thetabend;
      double sinchi = sin(chi);

      double PT = T->ev_Pt;
      double PL = T->ev_Pl;
      
      double Q2 = T->ev_Q2;

      double tau = Q2/(4.*pow(Mp,2));
      double epsilon = pow( 1. + 2.*(1.+tau)*pow(tan(etheta/2.),2), -1 );

      double kinfact = mu_p * sqrt(tau*(1.+epsilon)/(2.*epsilon));
      double R = -kinfact * PT/PL;
      
      hQ2->Fill( Q2, weight );
      hepsilon->Fill( epsilon, weight );
      hetheta->Fill( etheta * 180.0/TMath::Pi(), weight );
      hEprime->Fill( Eprime, weight );
      hptheta->Fill( ptheta * 180.0/TMath::Pi(), weight );
      hpp->Fill( pp, weight );
      
      if( idx_FPP1_track >= 0 ){ 
	TVector3 FPP1_track( (*(T->Harm_FPP1_Track_Xp))[idx_FPP1_track],
			     (*(T->Harm_FPP1_Track_Yp))[idx_FPP1_track],
			     1.0 );
	FPP1_track = FPP1_track.Unit();

	double thetaFPP1 = acos( FPP1_track.Dot( FT_track ) );
	double pT1 = pp_FT*sin(thetaFPP1);
	double Ay1 = Ayfunc->Eval( pT1, pp_FT );

	bool goodFPP1 = false;
	
	if( pT1 >= pTmin && pT1 <= pTmax ){
	  NAy2_sum += weight * pow( beampol * Ay1, 2 );
	  goodFPP1 = true;
	  Ngoodevent_sum += weight;
	  PT_sum += PT*weight * pow( beampol * Ay1, 2 );
	  PL_sum += PL*weight * pow( beampol * Ay1, 2 );
	  sinchi_sum += sinchi * weight * pow( beampol * Ay1, 2 );
	  Q2_sum += Q2 * weight * pow( beampol * Ay1, 2 );
	  epsilon_sum += epsilon * weight * pow( beampol * Ay1, 2 );
	  kinfact_sum += mu_p * sqrt(tau*(1.+epsilon)/(2.*epsilon)) * weight * pow( beampol * Ay1, 2 );
	  FFratio_sum += R * weight * pow( beampol * Ay1, 2 );
	}
	
	hpT_FPP1->Fill( pT1, weight );
	
	// if( idx_FPP2_track >= 0 && pT1 <= pTmax ){
	//   TVector3 FPP2_track( (*(T->Harm_FPP2_Track_Xp))[idx_FPP2_track],
	// 		       (*(T->Harm_FPP2_Track_Yp))[idx_FPP2_track],
	// 		       1.0 );

	//   FPP2_track = FPP2_track.Unit();

	//   double thetaFPP2 = acos( FPP2_track.Dot( FPP1_track ) );
	//   double pT2 = pp_FT*sin(thetaFPP2);
	//   double Ay2 = Ayfunc->Eval( pT2, pp_FT );

	//   hpT_FPP2->Fill( pT2, weight );

	//   if( pT2 >= pTmin && pT2 <= pTmax && (!goodFPP1||pT1<0.1) ){
	//     NAy2_sum += weight * pow( beampol * Ay2, 2 );
	//     Ngoodevent_sum += weight;
	//     PT_sum += PT*weight * pow( beampol * Ay2, 2 );
	//     PL_sum += PL*weight * pow( beampol * Ay2, 2 );
	//     sinchi_sum += sinchi * weight * pow( beampol * Ay2, 2 );
	//     Q2_sum += Q2 * weight * pow( beampol * Ay2, 2 );
	//     epsilon_sum += epsilon * weight * pow( beampol * Ay2, 2 );
	//     kinfact_sum += mu_p * sqrt(tau*(1.+epsilon)/(2.*epsilon)) * weight * pow( beampol * Ay2, 2 );
	//     FFratio_sum += R * weight * pow( beampol * Ay2, 2 );
	//   }
	  
	// }
      }
    
    }
  }

  Ngoodevent_sum *= Ndays * 24. * 3600.;
  NAy2_sum *= Ndays * 24. * 3600.;
  PT_sum *= Ndays * 24. * 3600.;
  PL_sum *= Ndays * 24. * 3600.;
  sinchi_sum *= Ndays * 24. * 3600.;
  Q2_sum *= Ndays * 24. * 3600.;
  epsilon_sum *= Ndays * 24. * 3600.;
  kinfact_sum *= Ndays * 24. * 3600.;
  FFratio_sum *= Ndays * 24. * 3600.;

  PT_sum /= NAy2_sum;
  PL_sum /= NAy2_sum;
  sinchi_sum /= NAy2_sum;
  Q2_sum /= NAy2_sum;
  epsilon_sum /= NAy2_sum;
  kinfact_sum /= NAy2_sum;
  FFratio_sum /= NAy2_sum;

  cout << "N good events total for " << Ndays << " days of running = " << Ngoodevent_sum << endl;
  cout << "FOM = sum_{i=1}^N (PeAy)^2 = " << NAy2_sum << endl;
  cout << "Weighted-average analyzing power using GEp-III parametrization = "
       << sqrt(NAy2_sum/Ngoodevent_sum)/beampol << endl;

  cout << "Average Q^2 = " << hQ2->GetMean() << endl;
  cout << "Average epsilon = " << hepsilon->GetMean() << endl;
  cout << "Total coincidence event rate passing trigger cuts = " << hQ2->Integral() << endl;
  cout << "Combined FPP efficiency = " << Ngoodevent_sum / ( hQ2->Integral() * Ndays * 24. * 3600. ) << endl;

  cout << "Assuming " << 100.*effrecon << "% reconstruction efficiency:" << endl;
  cout << "Projected Uncertainty in cos(phi), sin(phi) asymmetries = " << sqrt(2.0/NAy2_sum/effrecon) << endl;
  cout << "weighted average PT = " << PT_sum << endl;
  cout << "weighted average PL = " << PL_sum << endl;
  cout << "weighted average sin(chi) = " << sinchi_sum << endl;
  cout << "weighted average Q^2 (polarimeter events only) = " << Q2_sum << endl;
  cout << "weighted average epsilon (FPP events only) = " << epsilon_sum << endl;
  cout << "weighted average kinematic factor = " << kinfact_sum << endl;

  cout << "weighted average FF ratio (Puckett fit, PRC 2017) = " << FFratio_sum << endl; 

  double dPT = sqrt(2.0/NAy2_sum/effrecon);
  double dPL = dPT/sinchi_sum;

  double dR = fabs(FFratio_sum)*sqrt(pow( dPT/PT_sum, 2 ) + pow( dPL/PL_sum, 2 ) );

  cout << "Assuming " << 100.*effrecon << "% reconstruction efficiency and " << 100.*beampol << "% beam polarization:" << endl;
  cout << "Projected FF ratio uncertainty (absolute Delta (mu GE/GM) ) = " << dR << endl;
  
  
  elist_temp->Delete();
  fout->Write();
}
