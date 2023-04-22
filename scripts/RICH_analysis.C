#include "g4sbs_tree.C"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TString.h"
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include "TCut.h"
#include "TEventList.h"
#include "TMath.h"
#include "TVector3.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TDirectory.h"

double IRT_func(double *x, double *par){
  double R = par[0];
  double alpha = par[1];
  double e = par[2];
  double d = par[3];

  double beta = x[0];

  return e*d*sin( alpha - 2.0*beta ) + R*( e*sin(beta) - d*sin(alpha - beta ) );
}

void RICH_analysis( const char *inputfilename, const char *outputfilename){
  double mpi = 0.13957018;
  double mK  = 0.493677;
  double mp  = 0.938272046;

  //  double n_aero = 1.0304;
  double n_aero = 1.03128; //effective average from fit to thetaC (true) vs. p
  double n_gas  = 1.00137;

  double aero_yieldcoeff = 145.5; // average hits in aerogel = this number * sin2 thetaC
  double gas_yieldcoeff  = 4904.6; // average hits in gas = this number * sin2 thetaC

  double PI = TMath::Pi();

  TF1 *IRT_thetaCfunc = new TF1( "IRT_thetaCfunc", IRT_func, 0.0, PI, 4 ); 

  ifstream inputfile(inputfilename);
  TFile *fout = new TFile(outputfilename,"RECREATE");

  TChain *C = new TChain("T");

  //Read list of files:
  TString currentline;
  while( currentline.ReadLine(inputfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add( currentline.Data() );
    }
  }

  TCut globalcut = "";
  
  while( currentline.ReadLine(inputfile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline.Data();
    }
  }

  currentline.ReadLine(inputfile);

  ifstream sbsopticsfile(currentline.Data());
  //ifstream bbopticsfile("optics_BB.txt");

  vector<vector<double> > SBScoeff(4);
  vector<vector<int> > SBSexpon(5);

  //vector<vector<double> > BBcoeff(4);
  //vector<vector<int> > BBexpon(5);

  //int nparams_BB=0;
  int nparams_SBS=0;
  sbsopticsfile >> nparams_SBS;
  for(int i=0; i<5; i++){
    if( i < 4) SBScoeff[i].resize(nparams_SBS);
    SBSexpon[i].resize(nparams_SBS);
  }

  for(int ipar=0; ipar<nparams_SBS; ipar++){
    for(int icoeff=0; icoeff<4; icoeff++){
      sbsopticsfile >> SBScoeff[icoeff][ipar];
    }
    for(int iexpon=0; iexpon<5; iexpon++){
      sbsopticsfile >> SBSexpon[iexpon][ipar];
    }
  }

  int true_particle_flag=0;
  inputfile >> true_particle_flag;
  if( true_particle_flag < 0 ) true_particle_flag = 0;
  if( true_particle_flag > 2 ) true_particle_flag = 2;

  int nbins=200;
  inputfile >> nbins;

  double ebeam, sbstheta, sbsdtheta, sbsdphi;
  double pmin, pmax;
  inputfile >> ebeam >> sbstheta >> sbsdtheta >> sbsdphi;
  inputfile >> pmin >> pmax;
  
  double dRICH, dTRKR;
  inputfile >> dRICH >> dTRKR;

  double sigma_thetaC, sigma_pp, nsigma_cut; //cut on sigma for average angle resolution:
  inputfile >> sigma_thetaC >> sigma_pp >> nsigma_cut;

  double Lsth = 1.0; //"likelihood" for sub-threshold particles:
  inputfile >> Lsth;

  double Lratio_cut=1.0;
  inputfile >> Lratio_cut;

  double t0, tcut;
  inputfile >> t0 >> tcut;

  //There are three relevant points in space for RICH ray-tracing/reconstruction:
  //1. Position of tracker FP relative to RICH entry window (dRICH - dTRKR)
  //2. Position of spherical mirror center relative to RICH entry window
  //3. Position of center of PMT plane
  //Let's work in spectrometer coordinates, where +X = vertically down, +Y = beam left, and +Z = direction of particle motion (spectrometer axis)
  double zRICH = dRICH - dTRKR; //This is the relative z spacing between tracker focal plane and RICH entry window
  
  //Need to compute RICH center coordinates (in m)
  double Lx_aero = 5.0*0.114 + 4.0*2.54e-5;
  double xoff_aero = 0.07191;

  double MirrorRadius = 2.2;
  TVector3 MirrorCenterCoords( 0.0, 1.36403-(xoff_aero+Lx_aero/2.0), -1.03112 + zRICH ); //relative to center of entry window
  TVector3 FocalPointPosition( 0.0, 1.19350-(xoff_aero+Lx_aero/2.0), 0.42521 + zRICH ); //relative to center of entry window
  TVector3 HalfAeroPos( 0.0, 0.0, zRICH + 0.001 + 3.63*0.0113 );
  TVector3 GasStartPos( 0.0, 0.0, zRICH + 0.001 + 5.0*0.0113 + 0.0032 );

  TVector3 FocalPlanePosGlobal( dTRKR * sin( sbstheta*PI/180.0 ), 0.0, dTRKR * cos( sbstheta*PI/180.0 ) );
  FocalPlanePosGlobal.Print();

  TVector3 PMT_zaxis( 0.0, -sin(50.0*PI/180.0), cos(50.0*PI/180.0) );
  TVector3 PMT_xaxis( 1, 0, 0 );
  TVector3 PMT_yaxis = ( PMT_zaxis.Cross(PMT_xaxis) ).Unit();

  //Parameters needed to compute PMT coordinates:
  double ymin_PMT = -72.5376/100.0, ymax_PMT = 72.5376/100.0;
  double xmin_PMT[2] = { -29.083/100.0, -30.24632/100.0 };
  double xmax_PMT[2] = { 29.083/100.0, 30.24632/100.0 };
  int nrows_PMT[2] = {26, 27};

  //sbstheta *= PI/180.0;
  //sbsdtheta *= PI/180.0;
  //sbsdphi *= PI/180.0;

  TEventList *elist= new TEventList("elist");

  C->Draw(">>elist",globalcut);
  
  g4sbs_tree *T = new g4sbs_tree(C);

  TH1D::SetDefaultSumw2();

  //Initially, we want to determine the number of aerogel and gas photons produced by primary particle tracks as a function of theta, phi and momentum within SBS acceptance:
  TH1D *hp = new TH1D("hp","",nbins,pmin,pmax);
  TH1D *hth = new TH1D("hth","",nbins,sbstheta-sbsdtheta,sbstheta+sbsdtheta);
  TH1D *hph = new TH1D("hph","",nbins,-sbsdphi,sbsdphi);
  TH2D *hthph = new TH2D("hthph","",nbins,sbstheta-sbsdtheta,sbstheta+sbsdtheta, nbins, -sbsdphi, sbsdphi );
  TH2D *hpth = new TH2D("hpth","",nbins,sbstheta-sbsdtheta,sbstheta+sbsdtheta, nbins, pmin, pmax );
  TH2D *hpph = new TH2D("hpph","",nbins,-sbsdphi,sbsdphi,nbins,pmin,pmax);

  TH1D *hxfp = new TH1D("hxfp", "", nbins, -1, 1);
  TH1D *hxpfp = new TH1D("hxpfp", "", nbins, -.7,.7);
  TH1D *hyfp = new TH1D("hyfp", "", nbins, -.3, .3 );
  TH1D *hypfp = new TH1D("hypfp", "", nbins, -.1, .1 );
  
  TH1D *Nhits_tot  = new TH1D("Nhits_tot","",51,-0.5,50.5);
  TH1D *Nhits_aero = new TH1D("Nhits_aero","",51,-0.5,50.5);
  TH1D *Nhits_gas  = new TH1D("Nhits_gas","",51,-0.5,50.5);
  TH1D *Nphe_tot   = new TH1D("Nphe_tot","",16,-0.5,15.5);

  TH1D *Nphe_aero  = new TH1D("Nphe_aero","",11,-0.5,10.5);
  TH1D *Nphe_gas   = new TH1D("Nphe_gas","",16,-0.5,15.5);

  TH1D *htrackx    = new TH1D("htrackx","",nbins,-1.0,1.0);

  TProfile *Nhits_tot_p = new TProfile("Nhits_tot_p","",nbins,pmin,pmax);
  TProfile *Nhits_tot_th = new TProfile("Nhits_tot_th","",nbins,sbstheta-sbsdtheta,sbstheta+sbsdtheta);
  TProfile *Nhits_tot_ph = new TProfile("Nhits_tot_ph","",nbins,-sbsdphi,sbsdphi);
  
  TProfile *Nhits_aero_p = new TProfile("Nhits_aero_p","",nbins,pmin,pmax);
  TProfile *Nhits_aero_th = new TProfile("Nhits_aero_th","",nbins,sbstheta-sbsdtheta,sbstheta+sbsdtheta);
  TProfile *Nhits_aero_ph = new TProfile("Nhits_aero_ph","",nbins,-sbsdphi,sbsdphi);

  TProfile *Nhits_gas_p = new TProfile("Nhits_gas_p","",nbins,pmin,pmax);
  TProfile *Nhits_gas_th = new TProfile("Nhits_gas_th","",nbins,sbstheta-sbsdtheta,sbstheta+sbsdtheta);
  TProfile *Nhits_gas_ph = new TProfile("Nhits_gas_ph","",nbins,-sbsdphi,sbsdphi);

  TH2D *hNhits_tot_p = new TH2D("hNhits_tot_p","",nbins,pmin,pmax,51,-0.5,50.5);
  TH2D *hNhits_tot_th = new TH2D("hNhits_tot_th","",nbins,sbstheta-sbsdtheta,sbstheta+sbsdtheta,51,-0.5,50.5);
  TH2D *hNhits_tot_ph = new TH2D("hNhits_tot_ph","",nbins,-sbsdphi,sbsdphi,51,-0.5,50.5);

  TH2D *hNhits_aero_p = new TH2D("hNhits_aero_p","",nbins,pmin,pmax,51,-0.5,50.5);
  TH2D *hNhits_aero_th = new TH2D("hNhits_aero_th","",nbins,sbstheta-sbsdtheta,sbstheta+sbsdtheta,51,-0.5,50.5);
  TH2D *hNhits_aero_ph = new TH2D("hNhits_aero_ph","",nbins,-sbsdphi,sbsdphi,51,-0.5,50.5);

  TH2D *hNhits_gas_p = new TH2D("hNhits_gas_p","",nbins,pmin,pmax,51,-0.5,50.5);
  TH2D *hNhits_gas_th = new TH2D("hNhits_gas_th","",nbins,sbstheta-sbsdtheta,sbstheta+sbsdtheta,51,-0.5,50.5);
  TH2D *hNhits_gas_ph = new TH2D("hNhits_gas_ph","",nbins,-sbsdphi,sbsdphi,51,-0.5,50.5);

  TH1D *hpdiff = new TH1D("hpdiff","",nbins,-0.2, 0.2);
  TH2D *hpdiff_xptar = new TH2D("hpdiff_xptar", "", nbins, -.3, .3, nbins, -.2, .2 );
  TH2D *hpdiff_yptar = new TH2D("hpdiff_yptar", "", nbins, -.2, .2, nbins, -.2, .2 );
  TH2D *hpdiff_ytar = new TH2D("hpdiff_ytar", "", nbins, -.3, .3, nbins, -.2, .2 );
  TH2D *hpdiff_xtar = new TH2D("hpdiff_xtar", "", nbins, -.2, .2, nbins, -.2, .2 );
  TH2D *hpdiff_p    = new TH2D("hpdiff_p", "", nbins, 0.5, 11.0, nbins, -.2, .2 );
 
  TH1D *hvzdiff = new TH1D("hvzdiff", "", nbins, -0.1, 0.1 );
  TH2D *hvzdiff_xptar = new TH2D("hvzdiff_xptar", "", nbins, -0.3, 0.3, nbins, -0.1, 0.1 );
  TH2D *hvzdiff_yptar = new TH2D("hvzdiff_yptar", "", nbins, -0.15, 0.15, nbins, -0.1, 0.1 );
  TH2D *hvzdiff_ytar = new TH2D("hvzdiff_ytar", "", nbins, -0.2, 0.2, nbins, -0.1, 0.1 );
  TH2D *hvzdiff_xtar = new TH2D("hvzdiff_xtar", "", nbins, -0.15, 0.15, nbins, -0.1, 0.1 );
  TH2D *hvzdiff_p    = new TH2D("hvzdiff_p", "", nbins, 0.5, 11.0, nbins, -0.1, 0.1 );

  TH1D *hxtardiff = new TH1D("hxtardiff", "", nbins, -0.03, 0.03 );

  //Histograms for Cherenkov angle: true, reconstructed and "theoretical":
  TH1D *hthetaC = new TH1D("hthetaC", "", nbins, 0.0, 20.0 );
  TH2D *hthetaC_p = new TH2D("hthetaC_p", "", nbins, 0.5, 11.0, nbins, 0.0, 20.0 );

  TH1D *hthetaC_recon = new TH1D("hthetaC_recon", "", nbins, 0.0, 20.0 );
  TH2D *hthetaC_recon_p = new TH2D("hthetaC_recon_p", "", nbins, 0.5, 11.0, nbins, 0.0, 20.0 );
  
  TH1D *hthetaC_diff = new TH1D("hthetaC_diff", "", nbins, -3.0, 3.0 );
  TH1D *hthetaC_diffexpect = new TH1D("hthetaC_diffexpect", "", nbins, -3.0, 3.0 ); //This one includes momentum resolution
  TH1D *hthetaC_diffexpect_aero = new TH1D("hthetaC_diffexpect_aero", "", nbins, -3.0, 3.0 );
  TH1D *hthetaC_diffexpect_gas = new TH1D("hthetaC_diffexpect_gas", "", nbins, -3.0, 3.0 );

  TH1D *hthetaC_diffexpect_true = new TH1D("hthetaC_diffexpect_true", "", nbins, -3.0, 3.0 );
  TH1D *hthetaC_diffexpect_true_aero = new TH1D("hthetaC_diffexpect_true_aero", "", nbins, -3.0, 3.0 );
  TH1D *hthetaC_diffexpect_true_gas = new TH1D("hthetaC_diffexpect_true_gas", "", nbins, -3.0, 3.0 );

  TH1D *hthetaC_diff_aero = new TH1D("hthetaC_diff_aero", "", nbins, -3.0, 3.0 );
  TH1D *hthetaC_diff_gas  = new TH1D("hthetaC_diff_gas", "", nbins, -3.0, 3.0 );
  TH2D *hthetaC_diff_p = new TH2D("hthetaC_diff_p", "", nbins, 0.5, 11.0, nbins, -3.0, 3.0 );
  TH2D *hthetaC_diff_xfp = new TH2D("hthetaC_diff_xfp", "", nbins, -1, 1, nbins, -3.0, 3.0 );
  TH2D *hthetaC_diff_yfp = new TH2D("hthetaC_diff_yfp", "", nbins, -.3, .3, nbins, -3.0, 3.0 );
  TH2D *hthetaC_diff_xpfp = new TH2D("hthetaC_diff_xpfp", "", nbins, -.7, .7, nbins, -3.0, 3.0 );
  TH2D *hthetaC_diff_ypfp = new TH2D("hthetaC_diff_ypfp", "", nbins, -.1, .1, nbins, -3.0, 3.0 );
  TH2D *hthetaC_recon_vs_true = new TH2D("hthetaC_recon_vs_true", "", nbins, 0.0, 20.0, nbins, 0.0, 20.0 );

  TH1D *hemission_point = new TH1D("hemission_point", "", nbins, 0.0, 1.2);
  TH1D *hemission_point_aero = new TH1D("hemission_point_aero", "", nbins, 0.0, 1.2);
  TH1D *hemission_point_gas = new TH1D("hemission_point_gas", "", nbins, 0.0, 1.2);

  TH2D *hemission_point_xfp = new TH2D("hemission_point_xfp", "", nbins, -1, 1, nbins, 0.0, 1.2 );
  TH2D *hemission_point_xfp_aero = new TH2D("hemission_point_xfp_aero", "", nbins, -1, 1, nbins, 0.0, 1.2 );
  TH2D *hemission_point_xfp_gas = new TH2D("hemission_point_xfp_gas", "", nbins, -1, 1, nbins, 0.0, 1.2 );

  TH2D *hemission_point_xpfp = new TH2D("hemission_point_xpfp", "", nbins, -0.7, 0.7, nbins, 0.0, 1.2 );
  TH2D *hemission_point_xpfp_aero = new TH2D("hemission_point_xpfp_aero", "", nbins, -0.7, 0.7, nbins, 0.0, 1.2 );
  TH2D *hemission_point_xpfp_gas = new TH2D("hemission_point_xpfp_gas", "", nbins, -0.7, 0.7, nbins, 0.0, 1.2 );

  TH2D *hemission_point_yfp = new TH2D("hemission_point_yfp", "", nbins, -.3,.3, nbins, 0.0, 1.2 );
  TH2D *hemission_point_yfp_aero = new TH2D("hemission_point_yfp_aero", "", nbins, -.3,.3, nbins, 0.0, 1.2 );
  TH2D *hemission_point_yfp_gas = new TH2D("hemission_point_yfp_gas", "", nbins, -.3,.3, nbins, 0.0, 1.2 );

  TH2D *hemission_point_ypfp = new TH2D("hemission_point_ypfp", "", nbins, -.1, .1, nbins, 0.0, 1.2 );
  TH2D *hemission_point_ypfp_aero = new TH2D("hemission_point_ypfp_aero", "", nbins, -.1, .1, nbins, 0.0, 1.2 );
  TH2D *hemission_point_ypfp_gas = new TH2D("hemission_point_ypfp_gas", "", nbins, -.1, .1, nbins, 0.0, 1.2 );

  // TH1D *hemission_point_diff = new TH1D("hemission_point_diff", "", nbins, -0.6, 0.6 );
  // TH1D *hemission_point_diff_aero = new TH1D("hemission_point_diff_aero", "", nbins, -0.6, 0.6 );
  // TH1D *hemission_point_diff_gas = new TH1D("hemission_point_diff_gas", "", nbins, -0.6, 0.6 );

  TProfile *htest_point_xfp_aero = new TProfile("htest_point_xfp_aero", "", nbins, -1, 1);
  TProfile *htest_point_yfp_aero = new TProfile("htest_point_yfp_aero", "", nbins, -.3, .3);
  TProfile *htest_point_xpfp_aero = new TProfile("htest_point_xpfp_aero", "", nbins, -.7, .7);
  TProfile *htest_point_ypfp_aero = new TProfile("htest_point_ypfp_aero", "", nbins, -.1, .1);

  TProfile *htest_point_xfp_gas = new TProfile("htest_point_xfp_gas", "", nbins, -1, 1);
  TProfile *htest_point_yfp_gas = new TProfile("htest_point_yfp_gas", "", nbins, -.3, .3);
  TProfile *htest_point_xpfp_gas = new TProfile("htest_point_xpfp_gas", "", nbins, -.7, .7);
  TProfile *htest_point_ypfp_gas = new TProfile("htest_point_ypfp_gas", "", nbins, -.1, .1);

  TH1D *hPID = new TH1D("hPID", "", 4, -0.5, 3.5 );
  TH1D *hLogLratio = new TH1D("hLogLratio","",nbins,0.0,20.0);
  TH2D *hPID_p = new TH2D("hPID_p","", nbins, pmin, pmax, 4, -0.5, 3.5 );
  TH2D *hPID_xfp = new TH2D("hPID_xfp","", nbins, -1, 1, 4, -0.5, 3.5 );
  TH2D *hPID_yfp = new TH2D("hPID_yfp","", nbins, -.3, .3, 4, -0.5, 3.5 );
  TH2D *hPID_xpfp = new TH2D("hPID_xpfp","", nbins, -.7, .7, 4, -0.5, 3.5 );
  TH2D *hPID_ypfp = new TH2D("hPID_ypfp","", nbins, -.1, .1, 4, -0.5, 3.5 );

  //TH1D *hPnoID = new TH1D("hPnoID", "", nbins, 
  

  sbstheta *= PI/180.0;
  
  //TH1D *hxptardiff = new TH1D("hxptardiff", "", nbins, -0.3, 0.3 );

  TVector3 sbs_xaxis(0,-1,0);
  TVector3 sbs_zaxis( sin(sbstheta), 0.0, cos(sbstheta) );
  TVector3 sbs_yaxis = (sbs_zaxis.Cross(sbs_xaxis)).Unit();

  // beta_t n = 1, beta_t = 1/n
  // beta_t^2 = 1/n^2 = p_t^2 / p_t^2 + m^2
  // n^2 p_t^2 = p_t^2 + m^2 --> p_t^2 * ( n^2 - 1 ) = m^2 --> p_t = m/sqrt(n^2 - 1)
  double gas_threshold[3] = { mpi/sqrt( pow(n_gas,2) - 1.0 ), mK/sqrt(pow(n_gas,2)-1.0 ), mp/sqrt(pow(n_gas,2)-1.0 ) };
  double aero_threshold[3] = { mpi/sqrt( pow(n_aero,2) - 1.0 ), mK/sqrt(pow(n_aero,2)-1.0 ), mp/sqrt(pow(n_aero,2)-1.0 ) };

  long nevent = 0;
  while( T->GetEntry( elist->GetEntry(nevent++) ) ){
    if( nevent%1000 == 0 ) cout << nevent << endl;
    if( T->ev_harmaccept != 0 ){
      //hp->Fill( T->ev_np );
      hth->Fill( T->ev_nth * 180.0/PI );
      hph->Fill( T->ev_nph * 180.0/PI );
      
      hthph->Fill( T->ev_nth * 180.0/PI, T->ev_nph * 180.0/PI );
      hpth->Fill( T->ev_nth * 180.0/PI, T->ev_np );
      hpph->Fill( T->ev_nph * 180.0/PI, T->ev_np );
      
      //     Nhits_tot->Fill( T->RICH_nhits );
     
      //Let's retrieve track information:
      bool goodtrack = false;
      double xfp, yfp, xpfp, ypfp, ptrack; 
      double xfpfit, yfpfit, xpfpfit, ypfpfit, precon;

      //Reconstruct momentum: two iterations, correct for xtar on second iteration:
      double xtar = -T->ev_vy;
      double xptar, yptar, ytar, xptar_recon, yptar_recon, ytar_recon, vz_recon, xtar_recon;

      TVector3 pvect( T->ev_npx, T->ev_npy, T->ev_npz );
      TVector3 vertex( T->ev_vx, T->ev_vy, T->ev_vz );
      
      TVector3 pspec( pvect.Dot(sbs_xaxis), pvect.Dot(sbs_yaxis), pvect.Dot(sbs_zaxis) );
      TVector3 vspec( vertex.Dot(sbs_xaxis), vertex.Dot(sbs_yaxis), vertex.Dot(sbs_zaxis) );
      
      xptar = pspec.X()/pspec.Z();
      yptar = pspec.Y()/pspec.Z();
      xtar = vspec.X() - xptar * vspec.Z();
      ytar = vspec.Y() - yptar * vspec.Z();
      
      for(int itrack=0; itrack<T->ntracks; itrack++){
	if( (*(T->trackerid))[itrack] == 1 ){
	  goodtrack = true;
	  xfp = (*(T->trackx))[itrack];
	  yfp = (*(T->tracky))[itrack];
	  xpfp = (*(T->trackxp))[itrack];
	  ypfp = (*(T->trackyp))[itrack];
	  ptrack = (*(T->trackp))[itrack];
	  
	  xfpfit = (*(T->trackxfit))[itrack];
	  yfpfit = (*(T->trackyfit))[itrack];
	  xpfpfit = (*(T->trackxpfit))[itrack];
	  ypfpfit = (*(T->trackypfit))[itrack];
	  
	  double recon_sum[4] = {0,0,0,0};

	  double pinv;

	  xtar_recon = vspec.X(); //This is the vertical beam position, assumed known from BPM/raster

	  for(int iter=0; iter<3; iter++){
	    for( int coeff=0; coeff<4; coeff++){
	      recon_sum[coeff] = 0.0;
	    }
	    
	    for(int par=0; par<nparams_SBS; par++){
	      for(int coeff=0; coeff<4; coeff++){
		recon_sum[coeff] += SBScoeff[coeff][par] * pow( xfpfit, SBSexpon[0][par] ) * pow( yfpfit, SBSexpon[1][par] ) * pow( xpfpfit, SBSexpon[2][par] ) * pow(ypfpfit, SBSexpon[3][par] ) * pow( xtar_recon, SBSexpon[4][par] );
	      }
	    }
	    
	    xptar_recon = recon_sum[0]; 
	    yptar_recon = recon_sum[1];
	    ytar_recon = recon_sum[2];
	    precon = 1.0/recon_sum[3];
	    // ytar = -vz*sin(sbstheta) - yptar*zspec, zspec = vz * cos(sbstheta)
	    // ytar = vz*(-sin(theta0) - yptar*cos(theta0)) --> vz = -ytar/(sin(theta0)+yptar*cos(theta0))
	    // ytar = -vz*(sin(theta0) + yptar*cos(theta0))
	    vz_recon = -ytar_recon/( sin( sbstheta) + yptar_recon*cos(sbstheta) );

	    xtar_recon = vspec.X() - vz_recon * cos( sbstheta ) * xptar;
	  }

	  hvzdiff->Fill( vz_recon - T->ev_vz );
	  hvzdiff_xptar->Fill( xptar, vz_recon - T->ev_vz );
	  hvzdiff_yptar->Fill( yptar, vz_recon - T->ev_vz );
	  hvzdiff_ytar->Fill( ytar, vz_recon - T->ev_vz );
	  hvzdiff_xtar->Fill( xtar, vz_recon - T->ev_vz );
	  hvzdiff_p->Fill( ptrack, vz_recon - T->ev_vz );
	  //	  cout << "(vz, vz_recon)=(" << T->ev_vz << ", " << vz_recon << ")" << endl;
	  
	  hxtardiff->Fill( xtar_recon - xtar );
	  
	  hpdiff->Fill( precon/ptrack - 1.0 );
	  hpdiff_xptar->Fill( xptar, precon/ptrack - 1.0 );
	  hpdiff_yptar->Fill( yptar, precon/ptrack - 1.0 );
	  hpdiff_ytar->Fill( ytar, precon/ptrack - 1.0 );
	  hpdiff_xtar->Fill( xtar, precon/ptrack - 1.0 );
	  hpdiff_p->Fill( ptrack, precon/ptrack - 1.0 );

	  break;
	}
      }

      hp->Fill(ptrack);
      hxfp->Fill(xfp);
      hyfp->Fill(yfp);
      hxpfp->Fill(xpfp);
      hypfp->Fill(ypfp);

      if( goodtrack ){
	//if( goodtrack ) cout << "ptrack, precon, precon/ptrack - 1 = " << ptrack << ", " << precon << ", " << precon/ptrack - 1.0 << endl;
	//hpdiff->Fill( precon/ptrack-1.0 );
	
	TVector3 trackpos_spec( xfp, yfp, 0.0 );
	TVector3 phat_spec( xpfp, ypfp, 1.0 );
	TVector3 ptrack_spec = ptrack * phat_spec.Unit();

	TVector3 phat_spec_recon( xpfpfit, ypfpfit, 1.0 );
	TVector3 ptrack_spec_recon = precon * phat_spec_recon.Unit();
	TVector3 trackpos_spec_recon(xfpfit, yfpfit, 0.0 );

	//Calculate two trial emission points: 
	// 1. at half the thickness of aerogel.
	// 2. at half the thickness of gas. 
	// For this we need the intersection point of the reconstructed track with the spherical mirror:
	// Sphere equation: R^2 = ( rtrack - rcenter )^2
	// We also need to know the intersection of the track with the plane at half aerogel thickness:
	// R^2 = ( r0track + s * ntrack - rcenter )^2
	// R^2 = (r0track - rcenter)^2 + s^2 + 2*s*ntrack dot (r0track - rcenter)
	
	double A = 1.0;
	double B = 2.0*ptrack_spec_recon.Unit().Dot( trackpos_spec_recon - MirrorCenterCoords );
	double C = (trackpos_spec_recon - MirrorCenterCoords).Mag2() - pow(MirrorRadius,2);

	//A and C are both positive, B is also (usually) positive, therefore the + solution should usually be the unambiguously correct one

	//s is mirror intersection point:
	double s = (-B + sqrt(pow(B,2) - 4.0*A*C) )/(2.0*A);

	//The distance along the track to half the aerogel is simple:
	//0 = nhat dot (r - r0) = nhat dot ( r0track + s * ntrack - r0plane )
	double s_halfaero = (HalfAeroPos - trackpos_spec_recon).Z() / ptrack_spec_recon.Unit().Z();
	double s_startgas = (GasStartPos - trackpos_spec_recon).Z() / ptrack_spec_recon.Unit().Z();
	double s_halfgas  = 0.5*(s + s_startgas);

	double s_startaero_true = (zRICH + 0.001 - trackpos_spec.Z() ) / ptrack_spec.Unit().Z();

	TVector3 emission_point_aero = trackpos_spec_recon + s_halfaero * ptrack_spec_recon.Unit();
	TVector3 emission_point_gas  = trackpos_spec_recon + s_halfgas * ptrack_spec_recon.Unit();

	//For three particle-type hypotheses (pion, kaon, proton), calculate the "likelihood" that 
	//RICH hit pattern corresponds to that particle based on mass and momentum of particle

	double betarecon[3] = { precon/sqrt(pow(precon,2)+pow(mpi,2) ),
				precon/sqrt(pow(precon,2)+pow(mK,2) ),
				precon/sqrt(pow(precon,2)+pow(mp,2) ) };

	double thetaC_gas_expect[3] = { acos( TMath::Min(1.0,1.0/(n_gas*betarecon[0])) ),
					acos( TMath::Min(1.0,1.0/(n_gas*betarecon[1])) ),
					acos( TMath::Min(1.0,1.0/(n_gas*betarecon[2])) ) };
	double thetaC_aero_expect[3] = { acos( TMath::Min(1.0,1.0/(n_aero*betarecon[0])) ),
					acos( TMath::Min(1.0,1.0/(n_aero*betarecon[1])) ),
					acos( TMath::Min(1.0,1.0/(n_aero*betarecon[2])) ) };

	int nhits_gas_ptype[3] = {0,0,0};
	int nhits_aero_ptype[3] = {0,0,0};

	double thetaC_sum_gas_ptype[3] = {0.0, 0.0, 0.0};
	double thetaC_sum_aero_ptype[3] = {0.0, 0.0, 0.0};

	int nhits_tot  = 0;
	int nhits_aero = 0;
	int nhits_gas  = 0;
	//int nphe_aero = 0;
	//int nphe_gas  = 0;
	//int nphe_tot  = 0;
	for(int hit = 0; hit<T->RICH_nhits; hit++){
	  if( (*(T->RICH_mID))[hit] == 2 && fabs( (*(T->RICH_tavg))[hit] - t0 )< tcut ){
	    //  nphe_tot += (*(T->RICH_nphe))[hit];
	    nhits_tot += 1;
	    Nphe_tot->Fill( (*(T->RICH_nphe))[hit] );
	    if( (*(T->RICH_vol))[hit] == 1 ){ // aerogel:
	      nhits_aero += 1;
	      //    nphe_aero  += (*(T->RICH_nphe))[hit];
	      Nphe_aero->Fill( (*(T->RICH_nphe))[hit] );
	    }
	    if( (*(T->RICH_vol))[hit] == 2 ){ // gas:
	      nhits_gas += 1;
	      //nphe_gas  += (*(T->RICH_nphe))[hit];
	      Nphe_gas->Fill( (*(T->RICH_nphe))[hit] );
	    }
	    
	    TVector3 RICH_vertexhit( (*(T->RICH_vxhit))[hit], (*(T->RICH_vyhit))[hit], (*(T->RICH_vzhit))[hit] );
	    RICH_vertexhit -= FocalPlanePosGlobal;
	    TVector3 RICH_vertexhit_spec( RICH_vertexhit.Dot(sbs_xaxis), RICH_vertexhit.Dot( sbs_yaxis), RICH_vertexhit.Dot(sbs_zaxis) );

	    //Calculate the "true" point of closest approach between the emission vertex and the track: 
	    //D^2 = (xpoint - x0track - s*ntrack)^2; 2D dD = 2*(xpoint - x0track - s*ntrack) dot ntrack * ds = 0
	    // 2DdD = d/ds [( xpoint - x0track )^2 + s^2 * ntrack^2 - 2*s*(xpoint - x0track)dot ntrack] = 2s ds - 2*(xpoint - x0track) dot ntrack
	    // --> sclose = (xpoint - x0track) dot ntrack
	    double sclose = (RICH_vertexhit_spec - trackpos_spec).Dot( ptrack_spec.Unit() );

	   
	    //cout << "sclose, s_startaero_true = " << sclose << ", " << s_startaero_true << endl;

	    hemission_point->Fill( sclose - s_startaero_true );
	    hemission_point_xfp->Fill( xfp, sclose - s_startaero_true );
	    hemission_point_yfp->Fill( yfp, sclose - s_startaero_true );
	    hemission_point_xpfp->Fill( xpfp, sclose - s_startaero_true );
	    hemission_point_ypfp->Fill( ypfp, sclose - s_startaero_true );
	    
	    if( (*(T->RICH_vol))[hit] == 1 ){
	      hemission_point_aero->Fill( sclose - s_startaero_true );
	      hemission_point_xfp_aero->Fill( xfp, sclose - s_startaero_true );
	      hemission_point_yfp_aero->Fill( yfp, sclose - s_startaero_true );
	      hemission_point_xpfp_aero->Fill( xpfp, sclose - s_startaero_true );
	      hemission_point_ypfp_aero->Fill( ypfp, sclose - s_startaero_true );

	      htest_point_xfp_aero->Fill( xfp, s_halfaero - s_startaero_true );
	      htest_point_yfp_aero->Fill( yfp, s_halfaero - s_startaero_true );
	      htest_point_xpfp_aero->Fill( xpfp, s_halfaero - s_startaero_true );
	      htest_point_ypfp_aero->Fill( ypfp, s_halfaero - s_startaero_true );
	    }
	    if( (*(T->RICH_vol))[hit] == 2 ){
	      hemission_point_gas->Fill( sclose - s_startaero_true );
	      hemission_point_xfp_gas->Fill( xfp, sclose - s_startaero_true );
	      hemission_point_yfp_gas->Fill( yfp, sclose - s_startaero_true );
	      hemission_point_xpfp_gas->Fill( xpfp, sclose - s_startaero_true );
	      hemission_point_ypfp_gas->Fill( ypfp, sclose - s_startaero_true );

	      htest_point_xfp_gas->Fill( xfp, s_halfgas - s_startaero_true );
	      htest_point_yfp_gas->Fill( yfp, s_halfgas - s_startaero_true );
	      htest_point_xpfp_gas->Fill( xpfp, s_halfgas - s_startaero_true );
	      htest_point_ypfp_gas->Fill( ypfp, s_halfgas - s_startaero_true );
	    }

	    TVector3 RICH_vphit( (*(T->RICH_vpxhit))[hit], (*(T->RICH_vpyhit))[hit], (*(T->RICH_vpzhit))[hit] );
	    TVector3 RICH_vphit_spec( RICH_vphit.Dot(sbs_xaxis), RICH_vphit.Dot(sbs_yaxis), RICH_vphit.Dot(sbs_zaxis) );
	    
	    double thetaC_hit = acos( RICH_vphit_spec.Unit().Dot( phat_spec.Unit() ) );

	    hthetaC->Fill( thetaC_hit * 180.0/PI );
	    hthetaC_p->Fill( ptrack, thetaC_hit * 180.0/PI );

	    //now we need to calculate the coordinates of the detection point for each hit:
	    int rowhit = (*(T->RICH_row))[hit];
	    int colhit = (*(T->RICH_col))[hit];

	    //In going from RICH coordinates to transport coordinates, ysbs = xrich, xsbs = -yrich
	    double ypmt = xmin_PMT[colhit%2] + rowhit * ( xmax_PMT[colhit%2] - xmin_PMT[colhit%2] )/(double(nrows_PMT[colhit%2]-1));
	    double xpmt = -(ymin_PMT + colhit*(ymax_PMT - ymin_PMT)/72.0);

	    TVector3 hitpos = FocalPointPosition + PMT_xaxis * xpmt + PMT_yaxis * ypmt;
	    //Inverse ray-tracing uses the following geometrical facts about the optics:
	    //Plane of incidence contains the emission point, the reflection point, the detection point, and the sphere center
	    //Angle of incidence equals angle of reflection, such that the angle between CE and CR equals the angle between CD and CR
	    TVector3 CE_aero = emission_point_aero - MirrorCenterCoords;
	    TVector3 CD = hitpos - MirrorCenterCoords;
	    TVector3 CE_gas = emission_point_gas - MirrorCenterCoords;

	    //Reconstruct the emission angle for two hypotheses: 1 = aerogel, 2 = gas.
	    TVector3 nplaneinc_aero = CE_aero.Cross( CD ).Unit();
	    TVector3 nplaneinc_gas  = CE_gas.Cross( CD ).Unit();

	    //angle between the rays from center to emission point and center to detection point for aerogel and gas assumptions:
	    double alpha_aero = acos( CE_aero.Unit().Dot( CD.Unit() ) );
	    double alpha_gas  = acos( CE_gas.Unit().Dot( CD.Unit() ) );

	    IRT_thetaCfunc->SetParameter( 0, MirrorRadius );
	    IRT_thetaCfunc->SetParameter( 1, alpha_aero );
	    IRT_thetaCfunc->SetParameter( 2, CE_aero.Mag() );
	    IRT_thetaCfunc->SetParameter( 3, CD.Mag() );

	    double beta_aero = IRT_thetaCfunc->GetX( 0.0, alpha_aero );
	    //double theta_aero = atan2( CE_aero.Mag()*sin(beta_aero), MirrorRadius - CE_aero.Mag()*cos( beta_aero ) ); //Angle of incidence of ray with mirror
	    TVector3 vperp_aero = nplaneinc_aero.Cross( CE_aero ).Unit();
	    
	    TVector3 ReflectPoint_aero = MirrorCenterCoords + MirrorRadius * ( cos(beta_aero) * CE_aero.Unit() + sin(beta_aero) * vperp_aero );

	    IRT_thetaCfunc->SetParameter( 0, MirrorRadius );
	    IRT_thetaCfunc->SetParameter( 1, alpha_gas );
	    IRT_thetaCfunc->SetParameter( 2, CE_gas.Mag() );
	    IRT_thetaCfunc->SetParameter( 3, CD.Mag() );
	    
	    double beta_gas = IRT_thetaCfunc->GetX( 0.0, alpha_gas );
	    TVector3 vperp_gas = nplaneinc_gas.Cross( CE_gas ).Unit();
	    
	    TVector3 ReflectPoint_gas = MirrorCenterCoords + MirrorRadius * ( cos(beta_gas) * CE_gas.Unit() + sin(beta_gas) * vperp_gas );
	    
	    double thetaC_aero = acos( (ReflectPoint_aero - emission_point_aero).Unit().Dot( ptrack_spec_recon.Unit() ) );
	    double thetaC_gas = acos( (ReflectPoint_gas - emission_point_gas).Unit().Dot( ptrack_spec_recon.Unit() ) );
	    
	    //Correct reconstructed aerogel angle for refraction at aerogel-gas boundary:
	    
	    TVector3 emitted_ray_aero = (ReflectPoint_aero - emission_point_aero).Unit();

	    double thetaZ_aero = acos( emitted_ray_aero.Z() );
	    //from gas into aerogel: ngas sin thetagas = naero sin thetaaero
	    //assuming that thetaX and thetaY stay the same, thetaZ becomes smaller
	    double thetaZ_refracted = asin( n_gas/n_aero * sin(thetaZ_aero) );

	    TVector3 zaxis_temp(0,0,1);
	    TVector3 yaxis_temp = zaxis_temp.Cross( emitted_ray_aero ).Unit(); //perpendicular to the plane of incidence
	    TVector3 xaxis_temp = yaxis_temp.Cross( zaxis_temp ).Unit();
	    //we want to compute the ray that is parallel to the plane of incidence and makes an angle thetaZ_refracted with the Z axis:
	    TVector3 refracted_ray = cos(thetaZ_refracted) * zaxis_temp + sin(thetaZ_refracted) * xaxis_temp;
	   
	    thetaC_aero = acos(refracted_ray.Unit().Dot( ptrack_spec_recon.Unit() ) );

	    for( int ptype=0; ptype<3; ptype++ ){
	      //gas:
	      if( precon >= gas_threshold[ptype] - nsigma_cut*sigma_pp*precon && fabs( thetaC_gas - thetaC_gas_expect[ptype] ) <= nsigma_cut * sigma_thetaC ){
		nhits_gas_ptype[ptype] += 1;
		thetaC_sum_gas_ptype[ptype] += thetaC_gas;
	      }

	      //aerogel:
	      if( precon >= aero_threshold[ptype] - nsigma_cut*sigma_pp*precon && fabs( thetaC_aero - thetaC_aero_expect[ptype] ) <= nsigma_cut * sigma_thetaC ){
		nhits_aero_ptype[ptype] += 1;
		thetaC_sum_aero_ptype[ptype] += thetaC_aero;
	      }
	      
	    }

	    

	    if( (*(T->RICH_vol))[hit] == 1 ){ //photon was produced in aerogel, fill histograms for aerogel:
	      hthetaC_diffexpect->Fill( (thetaC_aero-thetaC_aero_expect[true_particle_flag])*180.0/PI );
	      hthetaC_diffexpect_aero->Fill( (thetaC_aero-thetaC_aero_expect[true_particle_flag])*180.0/PI );

	      hthetaC_diffexpect_true->Fill( (thetaC_hit-thetaC_aero_expect[true_particle_flag])*180.0/PI );
	      hthetaC_diffexpect_true_aero->Fill( (thetaC_hit-thetaC_aero_expect[true_particle_flag])*180.0/PI );

	      hthetaC_recon->Fill( thetaC_aero*180.0/PI );
	      hthetaC_recon_p->Fill( precon, thetaC_aero*180.0/PI );
	      hthetaC_recon_vs_true->Fill( thetaC_hit*180.0/PI, thetaC_aero*180.0/PI );
	      
	      hthetaC_diff_aero->Fill( (thetaC_aero-thetaC_hit)*180.0/PI );

	      hthetaC_diff->Fill( (thetaC_aero-thetaC_hit)*180.0/PI );
	      hthetaC_diff_p->Fill( precon, (thetaC_aero-thetaC_hit)*180.0/PI );
	      hthetaC_diff_xfp->Fill( xfp, (thetaC_aero-thetaC_hit)*180.0/PI );
	      hthetaC_diff_yfp->Fill( yfp, (thetaC_aero-thetaC_hit)*180.0/PI );
	      hthetaC_diff_xpfp->Fill( xpfp, (thetaC_aero-thetaC_hit)*180.0/PI );
	      hthetaC_diff_ypfp->Fill( ypfp, (thetaC_aero-thetaC_hit)*180.0/PI );
	    } else if( (*(T->RICH_vol))[hit] == 2){ //photon was produced in gas, fill gas histograms:

	      hthetaC_diffexpect_true->Fill( (thetaC_hit-thetaC_gas_expect[true_particle_flag])*180.0/PI );
	      hthetaC_diffexpect_true_gas->Fill( (thetaC_hit-thetaC_gas_expect[true_particle_flag])*180.0/PI );

	      hthetaC_diffexpect->Fill( (thetaC_gas-thetaC_gas_expect[true_particle_flag])*180.0/PI );
	      hthetaC_diffexpect_gas->Fill( (thetaC_gas-thetaC_gas_expect[true_particle_flag])*180.0/PI );

	      hthetaC_recon->Fill( thetaC_gas*180.0/PI );
	      hthetaC_recon_p->Fill( ptrack, thetaC_gas*180.0/PI );
	      hthetaC_recon_vs_true->Fill( thetaC_hit*180.0/PI, thetaC_gas*180.0/PI );

	      hthetaC_diff_gas->Fill( (thetaC_gas-thetaC_hit)*180.0/PI );

	      hthetaC_diff->Fill( (thetaC_gas-thetaC_hit)*180.0/PI );
	      hthetaC_diff_p->Fill( ptrack, (thetaC_gas-thetaC_hit)*180.0/PI );
	      hthetaC_diff_xfp->Fill( xfp, (thetaC_gas-thetaC_hit)*180.0/PI );
	      hthetaC_diff_yfp->Fill( yfp, (thetaC_gas-thetaC_hit)*180.0/PI );
	      hthetaC_diff_xpfp->Fill( xpfp, (thetaC_gas-thetaC_hit)*180.0/PI );
	      hthetaC_diff_ypfp->Fill( ypfp, (thetaC_gas-thetaC_hit)*180.0/PI );
	    }

	  }
	}

	// double Lgas[3] = {1,1,1};
	// double Laero[3] = {1,1,1};
	//double Lcombined[3] = {0,0,0};
	
	double Lgas[3] = {1,1,1};
	double Laero[3] = {1,1,1};
	double Lcombined[3] = {0,0,0};

	int bestPID_aero = -1, bestPID_gas = -1, bestPID_combined=-1;
	int bestPID_combined_hits_required = -1;

	double Lmax_combined_hits_required = 0.0;

	double Lmax_aero = 0.0, Lmax_gas = 0.0, Lmax_combined = 0.0;

	double Lmin = exp( -0.5*pow( nsigma_cut, 2 ) );

	//Lmin = 0.0;

	double LogLratio_min = 0.0;

	double LogLratio = LogLratio_min;
	double LogLratio_hits_required = LogLratio_min;

	int num_candidates=0;

	//if( precon >= aero_threshold[1] && precon <= gas_threshold[0] ) cout << "RICH summary:" << endl;
	//Now let's try to assign PID: We have thresholds and reconstructed emission angles:
	for( int ptype=0; ptype<3; ptype++){

	  // if( precon >= aero_threshold[1] && precon <= gas_threshold[0] ){
	  //   cout << "(ptype,precon, nhits_aero, nhits_gas, thetaC_aero, thetaC_gas)=(" << ptype+1 << ", " << precon << ", " 
	  // 	 << nhits_aero_ptype[ptype] << ", " << nhits_gas_ptype[ptype] << ", " 
	  // 	 << thetaC_sum_aero_ptype[ptype] << ", " << thetaC_sum_gas_ptype[ptype] << ")" << endl;
	  // }

	  //If a particle is above threshold AND has hits, then assign likelihood for that radiator:
	  //compute the average number of hits expected for each particle-type/radiator hypothesis:
	  double nexpect_aero = aero_yieldcoeff * pow(sin(thetaC_aero_expect[ptype]),2);
	  double nexpect_gas  = gas_yieldcoeff * pow(sin(thetaC_gas_expect[ptype]),2);

	  double Prob_aero = TMath::PoissonI( nhits_aero_ptype[ptype], nexpect_aero );
	  double Prob_gas  = TMath::PoissonI( nhits_gas_ptype[ptype], nexpect_gas );

	  if( nhits_gas_ptype[ptype] > 0 ){
	    //The way this is set up, nhits will only be non-zero if precon exceeds the threshold - nsigma * sigma_p, therefore, if the condition above is not satisfied, it guarantees that no hits are assigned to this track for this radiator. 
	    double thetaC_avg = thetaC_sum_gas_ptype[ptype]/double(nhits_gas_ptype[ptype]);
	    double sigma = sigma_thetaC / sqrt( double(nhits_gas_ptype[ptype]) );
	    Lgas[ptype] = exp( -0.5*pow( (thetaC_avg - thetaC_gas_expect[ptype])/sigma, 2 ) );
	  } else if( thetaC_gas_expect[ptype] > 0.0 ){ //Probability that no hits would occur based on expected angle:
	    double nexpect = gas_yieldcoeff * pow(sin(thetaC_gas_expect[ptype]),2);
	    //Poisson statistics: P(n) = a^n e^-a/n!
	    //Lgas[ptype] = TMath::Min(TMath::PoissonI(0,nexpect), Lmin);
	    Lgas[ptype] = TMath::PoissonI(0,nexpect);
	  } else { //particle was below threshold: 
	    Lgas[ptype] = Lsth;
	  }

	  if( nhits_aero_ptype[ptype] > 0 ){
	    double thetaC_avg = thetaC_sum_aero_ptype[ptype]/double(nhits_aero_ptype[ptype]);
	    double sigma = sigma_thetaC / sqrt(double(nhits_aero_ptype[ptype]));
	    Laero[ptype] = exp( -0.5*pow( (thetaC_avg - thetaC_aero_expect[ptype])/sigma, 2) );
	  } else if( thetaC_aero_expect[ptype] > 0.0 ){ //Particle was above threshold, but no hits observed. 
	    double nexpect = aero_yieldcoeff * pow(sin(thetaC_aero_expect[ptype]),2);
	    //Laero[ptype] = TMath::Min(TMath::PoissonI(0,nexpect), Lmin);
	    Laero[ptype] = TMath::PoissonI(0,nexpect);
	  } else { //particle was below threshold:
	    Laero[ptype] = Lsth;
	  }

	  //Now let us add a premium likelihood bonus to particles with hits in both radiators:
	  //int hasboth = 0;
	  int hasgas = 0, hasaero=0;
	  if( nhits_aero_ptype[ptype] > 0 ) hasaero = 1;
	  if( nhits_gas_ptype[ptype] > 0 ) hasgas = 1;

	  Lcombined[ptype] = Lgas[ptype]*Laero[ptype]*Prob_aero*Prob_gas + hasaero + hasgas;

	  if( Lgas[ptype] > Lmax_gas ){
	    Lmax_gas = Lgas[ptype];
	    bestPID_gas = ptype;
	  }
	  if( Laero[ptype] > Lmax_aero ){
	    Lmax_aero = Laero[ptype];
	    bestPID_aero = ptype;
	  }
	  if( Lcombined[ptype] > Lmax_combined ){
	    if( Lmax_combined > 0.0 ){ 
	      LogLratio = log( Lcombined[ptype] / Lmax_combined );
	    }

	    Lmax_combined = Lcombined[ptype];
	    bestPID_combined = ptype;
	  }

	  if( Lcombined[ptype] > Lmax_combined_hits_required && ( nhits_gas_ptype[ptype] > 0 || nhits_aero_ptype[ptype] > 0 ) ){
	    if( Lmax_combined_hits_required > 0.0 ){
	      LogLratio_hits_required = log( Lcombined[ptype] / Lmax_combined_hits_required );
	    }
	    Lmax_combined_hits_required = Lcombined[ptype];
	    bestPID_combined_hits_required = ptype;
	  }
	  
	}
	double Lc_max = 0.0;
	int bestPID = -1;

	if( LogLratio_hits_required > LogLratio_min ){
	  hLogLratio->Fill( LogLratio_hits_required );
	}

	// We reject events with no RICH signals; Safest way to go: high-purity, but relatively inefficient PID
	if( Lmax_combined_hits_required > 0 && bestPID_combined_hits_required >= 0 ){
	  if( (LogLratio_hits_required > LogLratio_min && LogLratio_hits_required > Lratio_cut) || LogLratio_hits_required == LogLratio_min ){
	    bestPID = bestPID_combined_hits_required;
	  } 
	} // else if( Lmax_combined > 0.0 && bestPID_combined >= 0 ){ //Then highest possible confidence, at least one particle type with gas and 
	//   //aerogel rings!
	bestPID = bestPID_combined;
	// } 
	
	//	hp->Fill( ptrack );

	hPID->Fill( bestPID + 1 );
	hPID_p->Fill( ptrack, bestPID + 1 );
	hPID_xfp->Fill( xfp, bestPID + 1 );
	hPID_yfp->Fill( yfp, bestPID + 1 );
	hPID_xpfp->Fill( xpfp, bestPID + 1 );
	hPID_ypfp->Fill( ypfp, bestPID + 1 );

	//Nphe_tot->Fill( nphe_tot );
	//Nphe_aero->Fill( nphe_aero );
	//Nphe_gas->Fill( nphe_gas );
	Nhits_aero->Fill( nhits_aero );
	Nhits_gas->Fill( nhits_gas );
	Nhits_tot->Fill( nhits_tot );

	Nhits_tot_p->Fill( T->ev_np, nhits_tot );
	Nhits_aero_p->Fill( T->ev_np, nhits_aero );
	Nhits_gas_p->Fill( T->ev_np, nhits_gas );

	hNhits_tot_p->Fill( T->ev_np, nhits_tot );
	hNhits_aero_p->Fill( T->ev_np, nhits_aero );
	hNhits_gas_p->Fill( T->ev_np, nhits_gas );
      
	Nhits_tot_th->Fill( T->ev_nth * 180.0/PI, nhits_tot );
	Nhits_aero_th->Fill( T->ev_nth * 180.0/PI, nhits_aero );
	Nhits_gas_th->Fill( T->ev_nth * 180.0/PI, nhits_gas );
      
	Nhits_tot_ph->Fill( T->ev_nph * 180.0/PI, nhits_tot );
	Nhits_aero_ph->Fill( T->ev_nph * 180.0/PI, nhits_aero );
	Nhits_gas_ph->Fill( T->ev_nph * 180.0/PI, nhits_gas );

	hNhits_tot_th->Fill( T->ev_nth * 180.0/PI, nhits_tot );
	hNhits_aero_th->Fill( T->ev_nth * 180.0/PI, nhits_aero );
	hNhits_gas_th->Fill( T->ev_nth * 180.0/PI, nhits_gas );
      
	hNhits_tot_ph->Fill( T->ev_nph * 180.0/PI, nhits_tot );
	hNhits_aero_ph->Fill( T->ev_nph * 180.0/PI, nhits_aero );
	hNhits_gas_ph->Fill( T->ev_nph * 180.0/PI, nhits_gas );
      }
    }
  }

  TH1D *htemp;

  vector<double> xtemp, ytemp_noID, ytemp_pi, ytemp_K, ytemp_pro;
  vector<double> extemp, eytemp_noID, eytemp_pi, eytemp_K, eytemp_pro;

  for(int i=1; i<=hPID_p->GetNbinsX(); i++){
    hPID_p->ProjectionY( "htemp", i, i );
    gDirectory->GetObject( "htemp", htemp );
    
    double tot = htemp->Integral();
    double noID = htemp->GetBinContent(1);
    double pi = htemp->GetBinContent(2);
    double K = htemp->GetBinContent(3);
    double pro = htemp->GetBinContent(4);

    double dtot = sqrt(tot);
    double dnoID = sqrt(noID);
    double dpi   = sqrt(pi);
    double dK    = sqrt(K);
    double dpro   = sqrt(pro);

    // if( htemp->Integral() > 0 ){
    //   htemp->Scale(1.0/htemp->Integral());
    // }
    if( tot > 0 ){
      
      xtemp.push_back( hp->GetBinCenter(i) );
      ytemp_noID.push_back( noID/tot );
      ytemp_pi.push_back( pi/tot );
      ytemp_K.push_back( K/tot );
      ytemp_pro.push_back( pro/tot );
      
      extemp.push_back(0.0);
      
      eytemp_noID.push_back( TMath::Max(1.0/tot,sqrt( noID/tot*(1.0-noID/tot)/tot ) ) );
      eytemp_pi.push_back( TMath::Max(1.0/tot,sqrt( pi/tot*(1.0-pi/tot)/tot) ) );
      eytemp_K.push_back( TMath::Max(1.0/tot,sqrt( K/tot*(1.0-K/tot)/tot) ) );
      eytemp_pro.push_back( TMath::Max(1.0/tot,sqrt( pro/tot*(1.0-pro/tot)/tot) ) );

    }
  }
  TGraphErrors *gPnoID_p = new TGraphErrors(xtemp.size(), &(xtemp[0]), &(ytemp_noID[0]), &(extemp[0]), &(eytemp_noID[0]) );
  TGraphErrors *gPpi_p = new TGraphErrors(xtemp.size(), &(xtemp[0]), &(ytemp_pi[0]), &(extemp[0]), &(eytemp_pi[0]) );
  TGraphErrors *gPK_p = new TGraphErrors(xtemp.size(), &(xtemp[0]), &(ytemp_K[0]), &(extemp[0]), &(eytemp_K[0]) );
  TGraphErrors *gPpro_p = new TGraphErrors(xtemp.size(), &(xtemp[0]), &(ytemp_pro[0]), &(extemp[0]), &(eytemp_pro[0]) );

  gPnoID_p->SetMarkerStyle(20);
  gPnoID_p->SetMarkerColor(1);
  gPnoID_p->SetLineColor(1);

  gPpi_p->SetMarkerStyle(21);
  gPpi_p->SetMarkerColor(2);
  gPpi_p->SetLineColor(2);

  gPK_p->SetMarkerStyle(22);
  gPK_p->SetMarkerColor(4);
  gPK_p->SetLineColor(4);

  gPpro_p->SetMarkerStyle(23);
  gPpro_p->SetMarkerColor(6);
  gPpro_p->SetLineColor(6);

  gPnoID_p->SetTitle("");
  gPpi_p->SetTitle("");
  gPK_p->SetTitle("");
  gPpro_p->SetTitle("");

  gPnoID_p->Write("gPnoID_p");
  gPpi_p->Write("gPpi_p");
  gPK_p->Write("gPK_p");
  gPpro_p->Write("gPpro_p");

  elist->Delete();
  fout->Write();
  fout->Close();
}
