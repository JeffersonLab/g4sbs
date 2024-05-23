#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "G4SBSRunData.hh"
#include "gep_tree_new.C"
//#include "gep_tree_singleFPP.C"
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
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set> 

const double Mp = 0.9382720813;
const double mu_p = 2.79284734462;
const double PI = TMath::Pi();

double KellyFunc(double *x, double *par){
  double Q2 = x[0];
  double tau = Q2/(4.*pow(Mp,2));

  return (1.0 + par[0]*tau)/(1.0 + par[1]*tau + par[2]*pow(tau,2) + par[3]*pow(tau,3)); //returns GEp or GMp/mu
}
//sclose and zclose calculation:
void calc_sclose_zclose( TVector3 Track1_Coord, TVector3 Track2_Coord, TVector3 Track1_Slope, TVector3 Track2_Slope, double &sclose, double &zclose ){
  double x1 = Track1_Coord.X() - Track1_Coord.Z()*Track1_Slope.X();
  double y1 = Track1_Coord.Y() - Track1_Coord.Z()*Track1_Slope.Y();
  double x2 = Track2_Coord.X() - Track2_Coord.Z()*Track2_Slope.X();
  double y2 = Track2_Coord.Y() - Track2_Coord.Z()*Track2_Slope.Y();

  double xp1 = Track1_Slope.X();
  double yp1 = Track1_Slope.Y();
  double xp2 = Track2_Slope.X();
  double yp2 = Track2_Slope.Y();

  TMatrixD Mclose(2,2);
  TVectorD bclose(2);

  Mclose(0,0) = 1.0 + pow(xp1,2) + pow(yp1,2);
  Mclose(0,1) = -(1.0 + xp1*xp2 + yp1*yp2);
  Mclose(1,0) = Mclose(0,1);
  Mclose(1,1) = 1.0 + pow(xp2,2) + pow(yp2,2);

  bclose(0) = xp1*(x2-x1) + yp1*(y2-y1);
  bclose(1) = xp2*(x1-x2) + yp2*(y1-y2);

  TVectorD zClose = Mclose.Invert() * bclose;

  double z1 = zClose(0);
  double z2 = zClose(1);

  double sClose2 = pow( x1 + xp1*z1 - (x2 + xp2*z2), 2 ) + pow( y1 + yp1*z1 - (y2 + yp2*z2), 2 ) + pow( z1-z2, 2 );

  sclose = sqrt(sClose2);
  zclose = 0.5*(z1 + z2 );
}

bool conetest( TVector3 Track1_Coord, TVector3 Track1_Slope, double theta, double zclose, double zback, double Lx=2.0, double Ly=0.6, double xcenter=0.0, double ycenter=0.0 ){
  double xfp, yfp, xpfp, ypfp;

  xfp = Track1_Coord.X() - Track1_Coord.Z() * Track1_Slope.X();
  yfp = Track1_Coord.Y() - Track1_Coord.Z() * Track1_Slope.Y();
  xpfp = Track1_Slope.X();
  ypfp = Track1_Slope.Y();

  double xclose = xfp + xpfp*zclose;
  double yclose = yfp + ypfp*zclose;

  double xpplus = (xpfp + tan(theta))/(1.-xpfp*tan(theta));
  double xpminus = (xpfp - tan(theta))/(1.+xpfp*tan(theta));
  double ypplus = (ypfp + tan(theta))/(1.-ypfp*tan(theta));
  double ypminus = (ypfp - tan(theta))/(1.+ypfp*tan(theta));

  double xmax = xclose + xpplus * (zback - zclose);
  double xmin = xclose + xpminus * (zback - zclose);
  double ymax = yclose + ypplus * (zback - zclose);
  double ymin = yclose + ypminus * (zback - zclose);

  return ( fabs( xmax - xcenter ) <= Lx/2.0 && fabs( xmin - xcenter ) <= Lx/2.0 && fabs( ymax - ycenter ) <= Ly/2.0 && fabs( ymin - ycenter ) <= Ly/2.0 );
  
}

//Analyzing power parametrization: x = p_T, y = p
TF2 *Ayfunc = new TF2("Ayfunc","([0]+[1]/y)*x*exp(-[2]*pow(x,2))",0.0,2.0, 1.0,15.0);
TF1 *GEPfunc = new TF1("GEPfunc",KellyFunc, 0.0,40.0,4);
TF1 *GMPfunc = new TF1("GMPfunc",KellyFunc, 0.0,40.0,4);

void GEP_FOM_quick_and_dirty(const char *configfilename, const char *outfilename="GEP_FOM.root", int rsampling=0, double maxweight=3.215e-34){
  ifstream configfile(configfilename);

  double geppar[4] = {-0.01,12.16,9.7,37.0};
  double gmppar[4] = {0.093,11.07,19.1,5.6};

  GEPfunc->SetParameters(geppar);
  GMPfunc->SetParameters(gmppar);
  
  double pTmin = 0.07;
  double pTmax = 1.0;

  double smax_FPP1 = 0.005; //meters
  double smax_FPP2 = 0.01; //meters

  double zmin_FPP1 = 0.6;
  double zmax_FPP1 = 1.25;

  double zmin_FPP2 = 1.6;
  double zmax_FPP2 = 2.4;

  //Later we will want to make these parameters user-configurable
  
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

  
  //gep_tree *T = new gep_tree(C);
  //don't initialize the "gep_tree" yet. We need some kind of flag to select the single-FPP or double-FPP options

  int NFPP = 2; 
  configfile >> NFPP;

  NFPP = std::max(1,std::min(NFPP,2));

  TEventList *elist_temp = new TEventList("elist_temp");

  C->Draw(">>elist_temp",globalcut);

  auto T =  new gep_tree_new( C );

  //if( treeflag == 1 ) T = T1; 

  long nevent = 0;

  double NAy2_sum = 0.0;
  double Ngoodevent_sum = 0.0;
  double Ngoodevent_sum1 = 0.0;
  double Ngoodevent_sum2 = 0.0;

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

  configfile >> smax_FPP1 >> smax_FPP2;
  configfile >> zmin_FPP1 >> zmax_FPP1;
  configfile >> zmin_FPP2 >> zmax_FPP2;

  int conetflag = 0; //default = don't require cone test
  configfile >> conetflag;
  
  double zback_FPP1=1.622, zback_FPP2=2.761;
  configfile >> zback_FPP1 >> zback_FPP2;

  
  
  TFile *fout = new TFile(outfilename,"RECREATE");
  
  TH1D *hECAL_sum = new TH1D("hECAL_sum","ECAL;Energy deposit in Pb-glass (GeV);",200,0.0,7.0);
  TH1D *hHCAL_sum_all = new TH1D("hHCAL_sum_all","HCAL; Energy deposit in Scintillator (GeV)",200,0.0,1.0);
  TH1D *hHCAL_sum_goodFPP1 = new TH1D("hHCAL_sum_goodFPP1","HCAL (good FPP events); Energy deposit in scintillator (GeV);",200,0.0,1.0);
  TH1D *hHCAL_sum_goodFPP2 = new TH1D("hHCAL_sum_goodFPP2","",200,0.0,1.0);

  TH2D *hPvstheta_FPP1 = new TH2D("hPvstheta_FPP1","",200,0.0,15.0,200,0.0,1.0);
  TH2D *hPvstheta_FPP2 = new TH2D("hPvstheta_FPP2","",200,0.0,15.0,200,0.0,1.0);

  TH1D *hQ2 = new TH1D("hQ2","Good coincidences; Q^{2} (GeV^{2})",150,0.0,15.0);
  TH1D *hepsilon = new TH1D("hepsilon",";#varepsilon;",250,0.0,1.0);
  TH1D *hetheta = new TH1D("hethetadeg",";#theta_{e} (deg)",250,0.0,90.0);
  TH1D *hEprime = new TH1D("hEprime",";E'_{e} (GeV);",250,0.0,11.0);
  TH1D *hptheta = new TH1D("hpthetadeg",";#theta_{p} (deg)",250,0.0,90.0);
  TH1D *hpp     = new TH1D("hpp",";p_{p} (GeV);",250,0.0,11.0);
  TH1D *hpT_FPP1     = new TH1D("hpT_FPP1",";p_{T} #equiv p_p sin(#theta_{FPP});",250,0.0,2.0);
  TH1D *hpT_FPP2     = new TH1D("hpT_FPP2","",250,0.0,2.0);
  TH1D *htheta_FPP1  = new TH1D("htheta_FPP1",";#theta_{FPP} (deg);",200,0.0,20.0); //plot in degrees
  TH1D *hphi_FPP1    = new TH1D("hphi_FPP1", ";#varphi_{FPP} (deg);", 36,-180.0,180.0);
  TH1D *htheta_FPP2  = new TH1D("htheta_FPP2","",250,0.0,15.0); //plot in degrees
  TH1D *hphi_FPP2    = new TH1D("hphi_FPP2", "", 36,-180.0,180.0);
  TH1D *htheta_FPP21 = new TH1D("htheta_FPP21","",250,0.0,15.0);
  TH1D *hphi_FPP21    = new TH1D("hphi_FPP21", "", 36,-180.0,180.0);
  
  TH1D *hsclose_FPP1 = new TH1D("hsclose_FPP1",";s_{close} (mm);",250,0.0,10.0); //millimeters
  TH1D *hsclose_FPP2 = new TH1D("hsclose_FPP2","",250,0.0,20.0); //millimeters
  TH1D *hsclose_FPP21 = new TH1D("hsclose_FPP21","",250,0.0,10.0);

  TH1D *hzclose_FPP1 = new TH1D("hzclose_FPP1",";z_{close} (m);",250,0.0,3.0); //meters
  TH1D *hzclose_FPP2 = new TH1D("hzclose_FPP2","",250,0.0,3.0); //meters
  TH1D *hzclose_FPP21 = new TH1D("hzclose_FPP21","",250,0.0,3.0);

  TH2D *hzclose_theta_FPP1 = new TH2D("hzclose_theta_FPP1",";z_{close} (m);#theta_{FPP} (deg)",150,0.0,3.0, 200,0.0,20.0);
  TH2D *hzclose_theta_FPP2 = new TH2D("hzclose_theta_FPP2","",250,0.0,3.0, 250,0.0,15.0);
  TH2D *hzclose_theta_FPP21 = new TH2D("hzclose_theta_FPP21","",250,0.0,3.0, 250,0.0,15.0);
  //Next step; add HCAL-based constraints

  //Add some diagnostic plots to check the conetest calculation:
  TH1D *htheta_FPP1_conetfail = new TH1D("htheta_FPP1_conetfail","",250,0.0,15.0);
  TH1D *htheta_FPP2_conetfail = new TH1D("htheta_FPP2_conetfail","",250,0.0,15.0);

  //Add some diagnostic plots to check the conetest calculation:
  TH1D *htheta_FPP1_conetpass = new TH1D("htheta_FPP1_conetpass","",250,0.0,15.0);
  TH1D *htheta_FPP2_conetpass = new TH1D("htheta_FPP2_conetpass","",250,0.0,15.0);

  //Add some diagnostic plots to check the conetest calculation:
  TH2D *hzclose_theta_FPP1_conetfail = new TH2D("hzclose_theta_FPP1_conetfail","",250,0.0,3.0,250,0.0,15.0);
  TH2D *hzclose_theta_FPP2_conetfail = new TH2D("hzclose_theta_FPP2_conetfail","",250,0.0,3.0,250,0.0,15.0);

  TH1D *hphiplus_FPP1 = new TH1D("hphiplus_FPP1","FPP1;#varphi_{FPP} (deg);",36,-180.,180.);
  TH1D *hphiminus_FPP1 = new TH1D("hphiminus_FPP1","FPP1;#varphi_{FPP} (deg);",36,-180.,180.);

  TH1D *hchideg = new TH1D("hchideg","Precession angle #chi (deg)",160,25.,41.0);
  TH2D *hchidegvsp = new TH2D("hchidegvsp","Precession versus proton momentum;p_{p} (GeV);#chi (deg)",250,2.0,3.5,160,25.,41.);
  
  //Let's add zclose and zclose calculations, and also 
  
  int treenumber = -1, oldtreenumber = -1;

  double Ebeam_current = Ebeam_default;
  double SBStheta = 16.9*TMath::Pi()/180.0;
  double SBStracker_pitch = 0.0*TMath::Pi()/180.0;
  //double EventWeight_current = 1.0;

  TRandom3 num(0);
  
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

    if( rsampling != 0 ){
      weight = T->ev_solang * maxweight * Lumi / Ngen_total;
    }
    
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

    if( T->Harm_FPP2_Track_ntracks == 1 && NFPP == 2 ) idx_FPP2_track = 0;

    if( idx_FT_track >= 0 ){
      TVector3 FT_track( (*(T->Harm_FT_Track_Xpfit))[idx_FT_track],
			 (*(T->Harm_FT_Track_Ypfit))[idx_FT_track],
			 1.0 );

      FT_track = FT_track.Unit();

      TVector3 FT_coord( (*(T->Harm_FT_Track_Xfit))[idx_FT_track],
			 (*(T->Harm_FT_Track_Yfit))[idx_FT_track],
			 0.0 );
      
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

      hchideg->Fill( chi * 180.0/TMath::Pi() );
      hchidegvsp->Fill( pp, chi * 180.0/TMath::Pi() );
      
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

      hECAL_sum->Fill( T->Earm_ECalTF1_det_esum, weight );
      hHCAL_sum_all->Fill( T->Harm_HCalScint_det_esum, weight );

      bool goodFPP1 = false;
      
      if( idx_FPP1_track >= 0 ){ 
	TVector3 FPP1_track( (*(T->Harm_FPP1_Track_Xpfit))[idx_FPP1_track],
			     (*(T->Harm_FPP1_Track_Ypfit))[idx_FPP1_track],
			     1.0 );
	FPP1_track = FPP1_track.Unit();

	TVector3 FPP1_coord( (*(T->Harm_FPP1_Track_Xfit))[idx_FPP1_track],
			     (*(T->Harm_FPP1_Track_Yfit))[idx_FPP1_track],
			     0.0 );
	
	double thetaFPP1 = acos( FPP1_track.Dot( FT_track ) );
	double pT1 = pp_FT*sin(thetaFPP1);
	double Ay1 = Ayfunc->Eval( pT1, pp_FT );

	//Calculate phi angles:
	TVector3 yaxistemp(0,1,0);
	TVector3 xaxistemp = yaxistemp.Cross(FT_track).Unit();
	yaxistemp = FT_track.Cross(xaxistemp).Unit();

	double phiFPP1 = TMath::ATan2( FPP1_track.Dot(yaxistemp), FPP1_track.Dot(xaxistemp) );
	
	double scloseFPP1,zcloseFPP1;

	calc_sclose_zclose( FT_coord, FPP1_coord, FT_track, FPP1_track, scloseFPP1, zcloseFPP1 );

	hsclose_FPP1->Fill( scloseFPP1*1000.0, weight );
	htheta_FPP1->Fill( thetaFPP1*180./PI, weight );
	
	hPvstheta_FPP1->Fill( thetaFPP1*180.0/PI, (*(T->Harm_FPP1_Track_P))[idx_FPP1_track]/pp_FT, weight );

	bool conetest1 = conetest( FT_coord, FT_track, thetaFPP1, zcloseFPP1, zback_FPP1 );

	if( !conetest1 ) {
	  htheta_FPP1_conetfail->Fill( thetaFPP1*180.0/PI, weight );
	  hzclose_theta_FPP1_conetfail->Fill( zcloseFPP1, thetaFPP1*180./PI, weight );
	} else {
	  htheta_FPP1_conetpass->Fill( thetaFPP1*180.0/PI, weight );
	}
	
	//Only fill zclose histograms if sclose < smax:
	if( scloseFPP1 <= smax_FPP1 ){
	  if( pT1 >= pTmin ) hzclose_FPP1->Fill( zcloseFPP1, weight );
	  if( conetest1 ) hzclose_theta_FPP1->Fill( zcloseFPP1, thetaFPP1*180./PI, weight );
	}

	
	if( pT1 >= pTmin && pT1 <= pTmax && scloseFPP1 <= smax_FPP1 && zcloseFPP1 >= zmin_FPP1
	    && zcloseFPP1 <= zmax_FPP1 && (conetest1||conetflag == 0 ) ){
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

	  hphi_FPP1->Fill( phiFPP1*180.0/PI, weight );

	  //Sample helicity; asymmetry goes like:
	  // Asym = Ay * Pe * (PT cos(phi) +PL sin(chi) sin(phi));
	  // Asym ~= PyFPP cos(phi) - PxFPP sin(phi)
	  // PyFPP ~= PT
	  // PxFPP ~= -PL sin chi
	  double Asym = Ay1 * beampol * (PT*cos(phiFPP1)+PL*sinchi*sin(phiFPP1));
	  //Asym = (f+ - f-)/(f+ + f-)
	  double prob_hplus = 0.5*(1.0+Asym);

	  if( num.Uniform() < prob_hplus ){
	    hphiplus_FPP1->Fill( phiFPP1*180.0/PI, weight );
	  } else {
	    hphiminus_FPP1->Fill( phiFPP1*180.0/PI, weight );
	  }
	  

	  hHCAL_sum_goodFPP1->Fill( T->Harm_HCalScint_det_esum, weight );

	  Ngoodevent_sum1 += weight;
	}
	
	hpT_FPP1->Fill( pT1, weight );
      }
      if( idx_FPP2_track >= 0 ){
	TVector3 FPP2_track( (*(T->Harm_FPP2_Track_Xpfit))[idx_FPP2_track],
			     (*(T->Harm_FPP2_Track_Ypfit))[idx_FPP2_track],
			     1.0 );
	
	FPP2_track = FPP2_track.Unit();
	
	TVector3 FPP2_coord( (*(T->Harm_FPP2_Track_Xfit))[idx_FPP2_track],
			     (*(T->Harm_FPP2_Track_Yfit))[idx_FPP2_track],
			     0.0 );
	
	double thetaFPP2 = acos( FPP2_track.Dot( FT_track ) );
	double pT2 = pp_FT*sin(thetaFPP2);
	double Ay2 = Ayfunc->Eval( pT2, pp_FT );
	
	hpT_FPP2->Fill( pT2, weight );
	
	//Calculate phi angles:
	TVector3 yaxistemp(0,1,0);
	TVector3 xaxistemp = yaxistemp.Cross(FT_track).Unit();
	yaxistemp = FT_track.Cross(xaxistemp).Unit();
	  
	double phiFPP2 = TMath::ATan2( FPP2_track.Dot(yaxistemp), FPP2_track.Dot(xaxistemp) );
	
	double scloseFPP2,zcloseFPP2;
	
	calc_sclose_zclose( FT_coord, FPP2_coord, FT_track, FPP2_track, scloseFPP2, zcloseFPP2 );
	
	hsclose_FPP2->Fill( scloseFPP2*1000.0, weight );
	htheta_FPP2->Fill( thetaFPP2*180./PI, weight );

	hPvstheta_FPP2->Fill( thetaFPP2*180.0/PI, (*(T->Harm_FPP2_Track_P))[idx_FPP2_track]/pp_FT, weight );

	bool conetest2 = conetest( FT_coord, FT_track, thetaFPP2, zcloseFPP2, zback_FPP2 );

	if( !conetest2 ) {
	  htheta_FPP2_conetfail->Fill( thetaFPP2*180.0/PI, weight );
	  hzclose_theta_FPP2_conetfail->Fill( zcloseFPP2, thetaFPP2*180./PI, weight );
	} else {
	  htheta_FPP2_conetpass->Fill( thetaFPP2*180.0/PI, weight );
	}
	
	if( scloseFPP2 <= smax_FPP2 ){
	  if( pT2 >= pTmin ) hzclose_FPP2->Fill( zcloseFPP2, weight );
	  if( conetest2 ) hzclose_theta_FPP2->Fill( zcloseFPP2, thetaFPP2*180./PI, weight );
	}
	  
	if( pT2 >= pTmin && pT2 <= pTmax && scloseFPP2 <= smax_FPP2 &&
	    zcloseFPP2 >= zmin_FPP2 && zcloseFPP2 <= zmax_FPP2 && (conetest2 || conetflag==0) && !goodFPP1 ){
	  NAy2_sum += weight * pow( beampol * Ay2, 2 );
	  Ngoodevent_sum += weight;
	  PT_sum += PT*weight * pow( beampol * Ay2, 2 );
	  PL_sum += PL*weight * pow( beampol * Ay2, 2 );
	  sinchi_sum += sinchi * weight * pow( beampol * Ay2, 2 );
	  Q2_sum += Q2 * weight * pow( beampol * Ay2, 2 );
	  epsilon_sum += epsilon * weight * pow( beampol * Ay2, 2 );
	  kinfact_sum += mu_p * sqrt(tau*(1.+epsilon)/(2.*epsilon)) * weight * pow( beampol * Ay2, 2 );
	  FFratio_sum += R * weight * pow( beampol * Ay2, 2 );

	  hHCAL_sum_goodFPP2->Fill( T->Harm_HCalScint_det_esum, weight );
	  hphi_FPP2->Fill( phiFPP2*180.0/PI, weight );

	  Ngoodevent_sum2 += weight;
	}
	
      }
    }
    
  }
  

  Ngoodevent_sum *= Ndays * 24. * 3600.;
  Ngoodevent_sum1 *= Ndays * 24. * 3600.;
  Ngoodevent_sum2 *= Ndays * 24. * 3600.; 
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

  cout << "FPP1 efficiency = " << Ngoodevent_sum1 / ( hQ2->Integral() * Ndays * 24. * 3600. ) << endl;
  cout << "FPP2 efficiency = " << Ngoodevent_sum2 / ( hQ2->Integral() * Ndays * 24. * 3600. ) << endl;
  
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

  cout << "P_T +/- dPT = " << PT_sum << " +/- " << dPT << endl;
  cout << "P_L +/- dPL = " << PL_sum << " +/- " << dPL << endl;
  
  double dR = fabs(FFratio_sum)*sqrt(pow( dPT/PT_sum, 2 ) + pow( dPL/PL_sum, 2 ) );

  cout << "Assuming " << 100.*effrecon << "% reconstruction efficiency:" << endl;
  cout << "Projected FF ratio uncertainty (absolute Delta (mu GE/GM) ) = " << dR << endl;

  TH1D *hphidiff = new TH1D(*hphiplus_FPP1);
  hphidiff->SetName("hAsym_FPP1");
  hphidiff->Add(hphiminus_FPP1,-1.);
  hphidiff->Divide(hphi_FPP1);
  
  
  elist_temp->Delete();
  fout->Write();
}
