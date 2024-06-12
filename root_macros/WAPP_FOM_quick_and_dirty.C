#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "G4SBSRunData.hh"
#include "WAPP_tree.C"
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
#include "TRandom3.h"

const double Mp = 0.9382720813;
const double mu_p = 2.79284734462;
const double Mpi = 0.13957;
const double Mn = 0.939565;
//Not sure if we need this for WAPP
// double KellyFunc(double *x, double *par){
//   double Q2 = x[0];
//   double tau = Q2/(4.*pow(Mp,2));

//   return (1.0 + par[0]*tau)/(1.0 + par[1]*tau + par[2]*pow(tau,2) + par[3]*pow(tau,3)); //returns GEp or GMp/mu
// }

//Analyzing power parametrization: x = p_T, y = p
TF2 *Ayfunc = new TF2("Ayfunc","([0]+[1]/y)*x*exp(-[2]*pow(x,2))",0.0,2.0, 1.0,15.0);
//Not sure we need this for WAPP:
// TF1 *GEPfunc = new TF1("GEPfunc",KellyFunc, 0.0,40.0,4);
// TF1 *GMPfunc = new TF1("GMPfunc",KellyFunc, 0.0,40.0,4);

void WAPP_FOM_quick_and_dirty(const char *configfilename, const char *outfilename="WAPP_FOM.root"){
  TRandom3 num(0);
  ifstream configfile(configfilename);

  // double geppar[4] = {-0.01,12.16,9.7,37.0};
  // double gmppar[4] = {0.093,11.07,19.1,5.6};

  // GEPfunc->SetParameters(geppar);
  // GMPfunc->SetParameters(gmppar);
  
  double pTmin = 0.15;
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
  map<TString, double> BBtheta_file;
  map<TString, double> SBStracker_pitch_file;
  map<TString, double> Ngen_file;
  map<TString, double> Norm_file; //normalization is defined such that rate = dsigma * Normalization
  //map<TString, double> EvtWeight_file; //should be norm_file * Ngen_file/Ngen_total, weight for individual events
  map<TString,double> Lumi_file;
  
  double Ebeam_default = 6.6;
  double Lumi = 5e37;
  
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

	  BBtheta_file[currentline] = rdtemp->fBBtheta;
	  
	  
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

  WAPP_tree *T = new WAPP_tree(C);

  long nevent = 0;

  double NAy2_sum = 0.0;
  double Ngoodevent_sum = 0.0;

  //double PT_sum = 0.0;
  //double PL_sum = 0.0;
  
  double sinchi_sum = 0.0;
  //  double Q2_sum = 0.0;
  // double epsilon_sum = 0.0;
  // double kinfact_sum = 0.0;
  // double FFratio_sum = 0.0;
  
  //read configuration parameters:
  double Ndays = 5.0;
  double effrecon = 0.9;
  configfile >> Ndays;
  configfile >> beampol;
  configfile >> effrecon;

  double thresh_SH = 0.5;
  double thresh_PS = 0.1;
  configfile >> thresh_SH;
  configfile >> thresh_PS;
  
  double gamma_pol = 0.9; //rough average polarization of photon beam as percentage of electron
  //beam polarization
  
  beampol *= gamma_pol;

  double K_LL = 0.8; //For pi-

  double K_LL_sum = 0.0;
  
  TFile *fout = new TFile(outfilename,"RECREATE");
  
  // TH1D *hQ2 = new TH1D("hQ2","",150,0.0,15.0);
  // TH1D *hepsilon = new TH1D("hepsilon","",250,0.0,1.0);
  // TH1D *hetheta = new TH1D("hethetadeg","",250,0.0,90.0);
  // TH1D *hEprime = new TH1D("hEprime","",250,0.0,11.0);
  // TH1D *hptheta = new TH1D("hpthetadeg","",250,0.0,90.0);

  TH1D *hs = new TH1D("hs","s (GeV^{2})",150,6.0,18.0);
  TH1D *ht = new TH1D("ht","-t (GeV^{2})",150,3.0,10.0);
  TH1D *hu = new TH1D("hu","-u (GeV^{2})", 150,0.0,10.0);
  TH1D *hsrecon = new TH1D("hsrecon","",150,6,18);
  TH1D *htrecon = new TH1D("htrecon","-t (GeV^{2})",150,3.0,10.0);
  TH1D *hurecon = new TH1D("hurecon","-u (GeV^{2})", 150,0.0,10.0);
  
  TH1D *hEgamma = new TH1D("hEgamma","E_{#gamma} (GeV)",150,3.9,6.7);
  TH1D *hcosthetaCM = new TH1D("hcosthetaCM","cos(#theta_{CM})",150,-1.0,1.0);
  //TH1D *hcosthetaCM_recon = new TH1D("hcosthetaCM_
  
  TH1D *hEgammadiff_pion = new TH1D("hEgammadiff_pi","E_{#gamma} (recon, #pi^{-} - true)",150,-1.0,1.0);
  TH2D *hEgamma_pion_vs_true = new TH2D("hEgamma_pion_vs_true","",150,4.0,6.6,150,2.0,8.6);
  TH1D *hEgammadiff_p = new TH1D("hEgammadiff_p","E_{#gamma} (recon, p - true)",150,-1.0,1.0);
  TH1D *hEgammadiff_ppi = new TH1D("hEgammadiff_ppi","E_{#gamma}^{recon} (p - #pi)",150,-1.0,1.0);
  TH1D *hEgammadiff_combined = new TH1D("hEgammadiff_combined","",150,-1,1);
  
  TH1D *hEmiss = new TH1D( "hEmiss", "E_{#gamma} + M_{n} - (E_{p} + E_{#pi})",150,-1.0,1.0);
  TH1D *hPmiss_perp = new TH1D("hPmiss_perp", "", 150,0.0,0.3);
  TH1D *hPmiss_par = new TH1D("hPmiss_par", "", 150,-0.5,0.5);
  TH1D *hPmiss = new TH1D("hPmiss","",150,0.0,0.5);
  TH1D *hMX2 = new TH1D("hMX2", "Missing mass squared", 150.0, 0.0, 0.4);
  TH1D *hpT2 = new TH1D("hpT2", "Squared transverse momentum of #pi^{-} + p", 150.0,0.0,0.1);

  TH1D *hdeltaphi = new TH1D("hdeltaphi","#phi_{p} - #phi_{#pi} - #pi",150,-0.2,0.2);
  TH1D *hacoplanarity = new TH1D("hacoplanarity","",300,0.0,0.3);
  TH2D *hphipi_vs_phip = new TH2D("hphipi_vs_phip","#phi_{#pi} vs. #phi_{p} - #pi",150,-30,30,150,-30,30 );
  
  TH1D *hpionP = new TH1D("hpionP","p_{#pi}",150,1.0,6.0);
  TH1D *hnucleonP = new TH1D("hnucleonP","p_{p}",150,2.0,6.6);

  TH1D *hsinchi = new TH1D("hsinchi","sin(#chi)",150,0.6,1.0);

  TH2D *hsinchi_vs_p = new TH2D("hsinchi_vs_p","sin(#chi) vs. p_{p}",100,2.5,5.0, 100, 0.6,1.0);

  TH1D *hsinchidiff = new TH1D( "hsinchidiff","",150,-0.2,0.2);
  TH2D *hsinchi_recon_vs_true = new TH2D("hsinchi_recon_vs_true","",150,0.6,1,150,0.6,1);
  
  //  TH1D *hpp     = new TH1D("hpp","",250,0.0,11.0);
  //We don't have double polarimeter for WAPP
  TH1D *hpT_FPP1     = new TH1D("hpT_FPP1","",100,0.0,1.2);
  //TH1D *hpT_FPP2     = new TH1D("hpT_FPP2","",250,0.0,2.0);

  TH1D *htheta_FPP1 = new TH1D("htheta_FPP1","",100,0.0,30.0);

  TH1D *hphiplus_FPP1 = new TH1D("hphiplus_FPP1","",36,-TMath::Pi(),TMath::Pi());
  TH1D *hphiminus_FPP1 = new TH1D("hphiminus_FPP1","",36,-TMath::Pi(),TMath::Pi());

  TH1D *hEBBCAL_smear = new TH1D("hEBBCAL_smear",";Total shower energy (GeV);",300,0,3.0);
  TH1D *hEPS_smear = new TH1D("hEPS_smear",";Preshower energy (GeV);",300,0,1.5);
  TH1D *hESH_smear = new TH1D("hESH_smear",";Shower energy (GeV);",300,0,3.0);
  TH2D *hESH_EPS_smear = new TH2D("hESH_EPS_smear",";Preshower Energy (GeV); Shower energy (GeV)", 300,0,1.5,300,0,3);

  
  int treenumber = -1, oldtreenumber = -1;

  double Ebeam_current = Ebeam_default;
  double SBStheta = 24.7*TMath::Pi()/180.0;
  double SBStracker_pitch = 0.0*TMath::Pi()/180.0;
  //double EventWeight_current = 1.0;

  double BBtheta = 42.5*TMath::Pi()/180.0;
  
  double weight = 1.0;
  
  while( T->GetEntry( elist_temp->GetEntry( nevent++ ) ) ){
    if( (nevent-1) % 1000 == 0 ) cout << "Event " << nevent << endl;


    TVector3 ppivect(T->ev_epx, T->ev_epy, T->ev_epz);
    double Epi = sqrt(pow(T->ev_ep,2)+pow(Mpi,2) );

    TLorentzVector Ppi_lab_true(ppivect,Epi);
    double ppi = Ppi_lab_true.P();
    double pitheta = Ppi_lab_true.Theta();
    double piphi = Ppi_lab_true.Phi();
    
    double Egamma_ppi_true = (2.0*Mn*Epi + pow(Mp,2)-pow(Mn,2) - pow(Mpi,2) )/(2.*Mn + 2.*ppi*cos(pitheta) - 2.*Epi );

    TVector3 ppvect(T->ev_npx, T->ev_npy, T->ev_npz );
    double Ep = sqrt(pow(T->ev_np,2)+pow(Mp,2));
    TLorentzVector Pp_lab_true(ppvect,Ep);

    double pp = Pp_lab_true.P();
    double ptheta = Pp_lab_true.Theta();
    double pphi = Pp_lab_true.Phi();
    
    double Egamma_pp_true = (2.0*Mn*Ep + pow(Mpi,2)-pow(Mp,2)-pow(Mn,2) )/(2.*Mn + 2.*pp*cos(ptheta) - 2.*Ep );

    

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
      BBtheta = BBtheta_file[fname_temp];
      oldtreenumber = treenum;
    }

    //compute "smeared" four-vectors:
    double ppi_smear = num.Gaus( ppi, 0.01*ppi ); //1% momentum resolution for BB
    double pitheta_smear = num.Gaus(pitheta, 0.002); //2 mrad polar angle
    double piphi_smear = num.Gaus(piphi,0.002/sin(BBtheta)); //2mrad / sin(BBtheta) for phi resolution

    double pp_smear = num.Gaus( pp, 0.007*pp ); //0.7% fractional momentum resolution for SBS
    double ptheta_smear = num.Gaus(ptheta,0.002); //~1 mrad polar angle
    double pphi_smear = num.Gaus(pphi,0.002/sin(SBStheta));

    if( pphi_smear < 0.0 ) pphi_smear += 2.0*TMath::Pi();
    
    //vertex correlation:
    double vzBB_smear = num.Gaus( T->ev_vz, 0.002/sin(BBtheta) );
    double vzSBS_smear = num.Gaus( T->ev_vz, 0.002/sin(SBStheta) );

    
    
    TVector3 pivect_smear( ppi_smear*sin(pitheta_smear)*cos(piphi_smear), ppi_smear*sin(pitheta_smear)*sin(piphi_smear), ppi_smear*cos(pitheta_smear) );
    double Epi_smear = sqrt(pow(ppi_smear,2) + pow(Mpi,2) );
    TLorentzVector Ppi_lab_smear(pivect_smear, Epi_smear);

    TVector3 pvect_smear( pp_smear*sin(ptheta_smear)*cos(pphi_smear), pp_smear*sin(ptheta_smear)*sin(pphi_smear), pp_smear*cos(ptheta_smear) );
    double Ep_smear = sqrt(pow(pp_smear,2) + pow(Mp,2) );
    TLorentzVector Pp_lab_smear( pvect_smear, Ep_smear);

    double Egamma_ppi_smear = (2.0*Mn*Epi_smear + pow(Mp,2)-pow(Mn,2)-pow(Mpi,2) )/(2.*Mn+2.*ppi_smear*cos(pitheta_smear) - 2.*Epi_smear);
    double Egamma_pp_smear = (2.0*Mn*Ep_smear + pow(Mpi,2)-pow(Mn,2)-pow(Mp,2) )/(2.*Mn+2.*pp_smear*cos(ptheta_smear)-2.*Ep_smear);

    TLorentzVector Pgamma_ppi_smear( 0.0, 0.0, Egamma_ppi_smear, Egamma_ppi_smear );
    TLorentzVector Pgamma_pp_smear( 0.0, 0.0, Egamma_pp_smear, Egamma_pp_smear );

    TLorentzVector PN( 0.0, 0.0, 0.0, Mn );

    TLorentzVector Pmiss = Pgamma_pp_smear + PN - Ppi_lab_smear - Pp_lab_smear;

    TLorentzVector Q4vect = Pgamma_pp_smear - Ppi_lab_smear; 

    TLorentzVector Psum = Ppi_lab_smear + Pp_lab_smear;
    
    double s_recon = Psum.M2();
    
    //s_Nrest = M_N^2 + 2M_N Egamma

    //Egamma = (s_Nrest - M_N^2)/(2M_N)
    //double Egamma_recon_combined = (s_recon - pow(Mn,2))/(2.*Mn);
    double Egamma_recon_combined = (s_recon - pow(Mn,2))/(2.0*(Psum.E() - Psum.Pz()));

    TLorentzVector Pgamma_combined(0,0,Egamma_recon_combined, Egamma_recon_combined);

    TVector3 norm_piplane = (Ppi_lab_smear.Vect().Unit().Cross(Pgamma_combined.Vect().Unit())).Unit();
    TVector3 norm_pplane = -(Pp_lab_smear.Vect().Unit().Cross(Pgamma_combined.Vect().Unit())).Unit();

    double acoplanarity = acos(norm_piplane.Dot(norm_pplane));
    
    //PN = Psum - Pgamma_combined; doing this makes PN zero by definition, not what we actually want

    Q4vect = Pgamma_combined - Ppi_lab_smear;

    double t_recon = (Pp_lab_smear-PN).M2();

    //s + t + u = MN^2 + Mp^2 + Mpi^2
    double u_recon = pow(Mn,2)+pow(Mp,2)+pow(Mpi,2) - s_recon - t_recon;

    Pmiss = Pgamma_combined + PN - Psum;
    
    
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
    for( int itrack=0; itrack<T->Harm_CEPolFront_Track_ntracks; itrack++ ){
      int MID = (*(T->Harm_CEPolFront_Track_MID))[itrack];
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
    if( T->Harm_CEPolRear_Track_ntracks == 1 ) idx_FPP1_track = 0;

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

    double ESH = T->Earm_BBSHTF1_det_esum;
    double EPS = T->Earm_BBPSTF1_det_esum;

    double ESHsm = num.Gaus( 200.0*ESH, sqrt(fabs(200.*ESH)) )/200.0;
    double EPSsm = num.Gaus( 200.0*EPS, sqrt(fabs(200.*EPS)) )/200.0;

    if( T->ev_Egamma >= 4.0 ){
    
      hEBBCAL_smear->Fill( ESHsm + EPSsm, weight );
      hEPS_smear->Fill( EPSsm, weight );
      hESH_smear->Fill( ESHsm, weight );
      hESH_EPS_smear->Fill( EPSsm, ESHsm, weight );

    }
    bool BBtrig = (EPSsm <= thresh_PS && ESHsm >= thresh_SH);
    if( thresh_PS < 0. ) BBtrig = ((EPSsm + ESHsm) >= thresh_SH);
    
    if( idx_FT_track >= 0 && BBtrig ){
      TVector3 FT_track( (*(T->Harm_CEPolFront_Track_Xp))[idx_FT_track],
			 (*(T->Harm_CEPolFront_Track_Yp))[idx_FT_track],
			 1.0 );

      FT_track = FT_track.Unit();

      TVector3 FT_track_SBS = FT_track.X() * SBS_FT_xaxis +
	FT_track.Y() * SBS_FT_yaxis +
	FT_track.Z() * SBS_FT_zaxis;

      double thetabend = acos( FT_track_SBS.Dot( ppunit_sbs ) );

      TVector3 FT_track_pol( (*(T->Harm_CEPolFront_Track_Sx))[idx_FT_track],
			     (*(T->Harm_CEPolFront_Track_Sy))[idx_FT_track],
			     (*(T->Harm_CEPolFront_Track_Sz))[idx_FT_track] );

      FT_track_pol = FT_track_pol.Unit();

      double chi_true = acos( FT_track_pol.Dot(FT_track) );
			     
      double pp_FT = (*(T->Harm_CEPolFront_Track_P))[idx_FT_track];
      double pitheta = T->ev_th;
      double Eprime = T->ev_ep;
      double ptheta = T->ev_nth;
      double pp = T->ev_np;

      hs->Fill( T->ev_s, weight );
      ht->Fill( -T->ev_t, weight );
      hu->Fill( -T->ev_u, weight );
      hEgamma->Fill( T->ev_Egamma, weight );
      hcosthetaCM->Fill( T->ev_costhetaCM, weight );

      hsrecon->Fill(s_recon,weight);
      htrecon->Fill(-t_recon,weight);
      hurecon->Fill(-u_recon,weight);
      
      hEgammadiff_combined->Fill( Egamma_recon_combined - T->ev_Egamma, weight );
      hEgammadiff_pion->Fill( Egamma_ppi_smear - T->ev_Egamma, weight );
      hEgamma_pion_vs_true->Fill( T->ev_Egamma, Egamma_ppi_true, weight );
      hEgammadiff_p->Fill( Egamma_pp_smear - T->ev_Egamma, weight );
      hEgammadiff_ppi->Fill( Egamma_pp_smear - Egamma_ppi_smear, weight );
      
      hpionP->Fill( T->ev_ep, weight );
      hnucleonP->Fill( T->ev_np, weight );

      //hdeltaphi->Fill( pphi_smear 
      hphipi_vs_phip->Fill( 180.0/TMath::Pi()*(pphi_smear-TMath::Pi()), 180.0/TMath::Pi()*piphi_smear, weight );
      hacoplanarity->Fill( acoplanarity, weight );
      hdeltaphi->Fill( pphi_smear - piphi_smear - TMath::Pi(), weight );
      
      hEmiss->Fill( Egamma_pp_smear + Mn - Psum.E(), weight );
      TVector3 q3vect = Q4vect.Vect();

      double Pmiss_par = Pmiss.Vect().Dot( q3vect.Unit() );
      double Pmiss_perp = (Pmiss.Vect() - q3vect.Unit() * (Pmiss.Vect().Dot(q3vect.Unit()))).Mag();

      hPmiss_par->Fill( Pmiss_par, weight );
      hPmiss_perp->Fill( Pmiss_perp, weight );

      hPmiss->Fill( Pmiss.Vect().Mag(), weight );
      
      hMX2->Fill( -Pmiss.M2(), weight );

      hpT2->Fill( Psum.Perp2(), weight );
      
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

      hsinchi->Fill( sinchi, weight );
      hsinchi_vs_p->Fill( T->ev_np, sinchi, weight );

      hsinchidiff->Fill( sinchi - sin(chi_true), weight );
      hsinchi_recon_vs_true->Fill( sinchi, sin(chi_true), weight );
      
      // double PT = T->ev_Pt;
      // double PL = T->ev_Pl;
      
      //     double Q2 = T->ev_Q2;

      //      double tau = Q2/(4.*pow(Mp,2));
      // double epsilon = pow( 1. + 2.*(1.+tau)*pow(tan(etheta/2.),2), -1 );

      //      double kinfact = mu_p * sqrt(tau*(1.+epsilon)/(2.*epsilon));
      //double R = -kinfact * PT/PL;
      
      // hQ2->Fill( Q2, weight );
      // hepsilon->Fill( epsilon, weight );
      // hetheta->Fill( etheta * 180.0/TMath::Pi(), weight );
      // hEprime->Fill( Eprime, weight );
      // hptheta->Fill( ptheta * 180.0/TMath::Pi(), weight );
      // hpp->Fill( pp, weight );
      
      if( idx_FPP1_track >= 0 ){ 
	TVector3 FPP1_track( (*(T->Harm_CEPolRear_Track_Xp))[idx_FPP1_track],
			     (*(T->Harm_CEPolRear_Track_Yp))[idx_FPP1_track],
			     1.0 );
	FPP1_track = FPP1_track.Unit();

	TVector3 yaxis(0,1,0);
	TVector3 xaxis = yaxis.Cross(FT_track).Unit();
	yaxis = FT_track.Cross(xaxis).Unit();
	
	double thetaFPP1 = acos( FPP1_track.Dot( FT_track ) );
	double pT1 = pp_FT*sin(thetaFPP1);
	double Ay1 = Ayfunc->Eval( pT1, pp_FT );

	double phiFPP1 = atan2(FPP1_track.Dot(yaxis),FPP1_track.Dot(xaxis));
	
	bool goodFPP1 = false;
	
	if( pT1 >= pTmin && pT1 <= pTmax ){
	  NAy2_sum += weight * pow( beampol * Ay1, 2 );
	  goodFPP1 = true;
	  Ngoodevent_sum += weight;
	  //PT_sum += PT*weight * pow( beampol * Ay1, 2 );
	  //PL_sum += PL*weight * pow( beampol * Ay1, 2 );
	  K_LL_sum += K_LL*weight * pow( beampol*Ay1,2 );
	  sinchi_sum += sinchi * weight * pow( beampol * Ay1, 2 );
	  //Q2_sum += Q2 * weight * pow( beampol * Ay1, 2 );
	  //epsilon_sum += epsilon * weight * pow( beampol * Ay1, 2 );
	  //kinfact_sum += mu_p * sqrt(tau*(1.+epsilon)/(2.*epsilon)) * weight * pow( beampol * Ay1, 2 );
	  //FFratio_sum += R * weight * pow( beampol * Ay1, 2 );
	  double asym_phi = beampol*Ay1*K_LL*(FT_track_pol.Dot(yaxis)*cos(phiFPP1) - FT_track_pol.Dot(xaxis)*sin(phiFPP1));

	  //How to fill asymmetry plot: A = (f+-f-)/(f+ + f-)
	  hphiplus_FPP1->Fill( phiFPP1, 0.5*weight*(1.+asym_phi) );
	  hphiminus_FPP1->Fill( phiFPP1, 0.5*weight*(1.-asym_phi) );
	  
	}
	
	hpT_FPP1->Fill( pT1, weight );

	htheta_FPP1->Fill( thetaFPP1*180.0/TMath::Pi(), weight );
	
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
  // PT_sum *= Ndays * 24. * 3600.;
  // PL_sum *= Ndays * 24. * 3600.;
  sinchi_sum *= Ndays * 24. * 3600.;
  // Q2_sum *= Ndays * 24. * 3600.;
  // epsilon_sum *= Ndays * 24. * 3600.;
  // kinfact_sum *= Ndays * 24. * 3600.;
  // FFratio_sum *= Ndays * 24. * 3600.;

  K_LL_sum *= Ndays * 24. * 3600.;

  // PT_sum /= NAy2_sum;
  // PL_sum /= NAy2_sum;
  sinchi_sum /= NAy2_sum;
  //Q2_sum /= NAy2_sum;
  // epsilon_sum /= NAy2_sum;
  // kinfact_sum /= NAy2_sum;
  // FFratio_sum /= NAy2_sum;

  K_LL_sum /= NAy2_sum;

  cout << "N good events total for " << Ndays << " days of running = " << Ngoodevent_sum << endl;
  cout << "FOM = sum_{i=1}^N (PeAy)^2 = " << NAy2_sum << endl;
  cout << "Weighted-average analyzing power using GEp-III parametrization = "
       << sqrt(NAy2_sum/Ngoodevent_sum)/beampol << endl;

  //  cout << "Average Q^2 = " << hQ2->GetMean() << endl;
  // cout << "Average epsilon = " << hepsilon->GetMean() << endl;
  cout << "Total coincidence event rate passing trigger cuts = " << ht->Integral() << endl;
  cout << "Combined FPP efficiency = " << Ngoodevent_sum / (ht->Integral() * Ndays * 24. * 3600. ) << endl;

  cout << "Assuming " << 100.*effrecon << "% reconstruction efficiency:" << endl;
  cout << "Projected Uncertainty in cos(phi), sin(phi) asymmetries = " << sqrt(2.0/NAy2_sum/effrecon) << endl;
  // cout << "weighted average PT = " << PT_sum << endl;
  // cout << "weighted average PL = " << PL_sum << endl;
  cout << "weighted average sin(chi) = " << sinchi_sum << endl;
  //cout << "weighted average Q^2 (polarimeter events only) = " << Q2_sum << endl;
  //cout << "weighted average epsilon (FPP events only) = " << epsilon_sum << endl;
  //cout << "weighted average kinematic factor = " << kinfact_sum << endl;

  //  cout << "weighted average FF ratio (Puckett fit, PRC 2017) = " << FFratio_sum << endl; 

  double dPT = sqrt(2.0/NAy2_sum/effrecon);
  double dPL = dPT/sinchi_sum;

  //double dR = fabs(FFratio_sum)*sqrt(pow( dPT/PT_sum, 2 ) + pow( dPL/PL_sum, 2 ) );

  cout << "Assuming " << 100.*effrecon << "% reconstruction efficiency:" << endl;
  //cout << "Projected FF ratio uncertainty (absolute Delta (mu GE/GM) ) = " << dR << endl;

  cout << "Projected KLL absolute statistical uncertainty for gamma n --> pi- p = " << dPL << endl;
  cout << "Relative statistical uncertainty = " << dPL/K_LL << endl;
  

  TH1D *hphisum_FPP1 = new TH1D( *hphiplus_FPP1 );
  hphisum_FPP1->SetName("hphisum_FPP1");
  hphisum_FPP1->Reset();
  
  TH1D *hphidiff_FPP1 = new TH1D( *hphiplus_FPP1 );
  hphidiff_FPP1->SetName("hphidiff_FPP1");
  hphidiff_FPP1->Reset();

  hphisum_FPP1->Add( hphiplus_FPP1, hphiminus_FPP1 );
  hphidiff_FPP1->Add( hphiplus_FPP1, hphiminus_FPP1, 1, -1 );

  TH1D *hasym_FPP1 = new TH1D( *hphiplus_FPP1 );
  hasym_FPP1->SetName("hasym_FPP1");
  hasym_FPP1->Reset();

  hasym_FPP1->Divide( hphidiff_FPP1, hphisum_FPP1, 1, 1 );

  //now compute appropriate errors on the bin contents:

  for( int ibin=1; ibin<=36; ibin++ ){
    double A = hasym_FPP1->GetBinContent( ibin );
    double N = hphisum_FPP1->GetBinContent( ibin ) * Ndays * 24. * 3600.;
    double dA = sqrt( (1.0-pow(A,2))/N );

    //hasym_FPP1->SetBinContent
    hasym_FPP1->SetBinError( ibin, dA );
  }

  TH1D *heffpi_vs_thresh_BBCAL = new TH1D("heffpi_vs_thresh_BBCAL",";threshold (GeV);Charged pion effficiency",300,0,3.0);

  for( int i=1; i<=300; i++ ){
    heffpi_vs_thresh_BBCAL->SetBinContent( i, hEBBCAL_smear->Integral(i,300)/hEBBCAL_smear->GetEntries() );
  }
  
  elist_temp->Delete();
  fout->Write();
}
