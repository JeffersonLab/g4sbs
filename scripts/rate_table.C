#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "g4sbs_tree.C"
#include <iostream>
#include <fstream>
#include <vector>
#include "TSystem.h"
#include "TString.h"
#include "TEventList.h"
#include "TCut.h"
#include "TObject.h"
#include "TROOT.h"
#include "G4SBSRunData.hh"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TRandom.h"
//class G4SBSRunData;
//#include "G4SBSDict.h"
//include "G4SBSDict.cxx"

double get_R(double x, double Q2, double &dR){
  double R;
  //double Ra, Rb, Rc;
  double theta, q2thr;
  theta = 1.0 + 12.0*Q2/(Q2+1.0)*pow(0.125,2)/(pow(0.125,2)+pow(x,2));

  double a[6] = {0.0485, 0.5470, 2.0621, -0.3804, 0.5090, -0.0285};
  double b[6] = {0.0481, 0.6114, -0.3509, -0.4611, 0.7172, -0.0317};
  double c[6] = {0.0577, 0.4644, 1.8288, 12.3708, -43.1043, 41.7415};

  double Theta = 1.0 + 12.0*Q2/(Q2+1.0)*pow(0.125,2)/(pow(0.125,2)+pow(x,2));
  double Q2thr = c[3]*x + c[4]*pow(x,2) + c[5]*pow(x,3);

  double Ra = a[0]/log(Q2/0.04) * Theta + a[1]/pow( pow(Q2,4) + pow(a[2],4), 0.25 )*( 1.0 + a[3]*x + a[4]*pow(x,2) )*pow(x,a[5]);
  double Rb = b[0]/log(Q2/0.04) * Theta + ( b[1]/Q2 + b[2]/(pow(Q2,2)+pow(0.3,2)) )*(1.0 + b[3]*x + b[4]*pow(x,2) )*pow(x,b[5] );
  double Rc = c[0]/log(Q2/0.04) * Theta + c[1]*pow( pow(Q2-Q2thr,2) + pow(c[2],2), -0.5 ); 
  
  dR = 0.0078 - 0.013*x + (0.070 - 0.39*x + 0.70*pow(x,2) )/(1.7 + Q2);

  return (Ra + Rb + Rc)/3.0;
}

void rate_table( const char *inputfilename, const char *outputfilename, const char *rootfilename ){

  // TString rootfilename(outputfilename);
  // rootfilename.ReplaceAll(".txt",".root");

  TFile *fout = new TFile(rootfilename, "RECREATE");

  TRandom3 num(0);

  double PI = TMath::Pi();

  double pixelsize = 0.15;
  double Hcal_height = 3.3;
  double Hcal_width = 1.65;
  //double mpi0 = 0.134977;
  
  //gSystem->Load("libg4sbsroot.so");

  int nbins_x, nbins_z;
  double xmin, xmax, zmin, zmax;
  
  double PB,PT;

  double chargesep_m, chargesep_b;

  ifstream inputfile(inputfilename);
  ofstream outputfile(outputfilename);
  
  //We already know rate, the only thing we need is duration of the expt.
  // to get the number of events:
  double ndays;

  int pi0flag = 0;
 
  TChain *C = new TChain("T");

  int nfiles=0;
  vector<int> ngen;
  vector<int> ntries;
  vector<double> weights; //correct for ratio of generated events to number of tries:

  double sum_weights = 0.0;

  G4SBSRunData *rtemp;

  //TObject *rtemp;

  TString currentline;
  while(currentline.ReadLine(inputfile)&&!currentline.BeginsWith("endlist")){
    if( !currentline.BeginsWith("#") ){
      TFile Ftemp(currentline.Data(),"READ");

      if( !Ftemp.IsZombie() ){
	Ftemp.GetObject("run_data",rtemp);
	
	if( rtemp ){

	  ngen.push_back( rtemp->fNthrown );
	  ntries.push_back( rtemp->fNtries );
	  weights.push_back( double( rtemp->fNthrown )/double( rtemp->fNtries ) );
	  cout << "weight = " << weights[nfiles] << endl;
	  C->Add(currentline.Data());
	
	  sum_weights += weights[nfiles];
	  nfiles++;
	}
      }
    }
  }

  

  TCut cut = "";

  while(currentline.ReadLine(inputfile)&&!currentline.BeginsWith("endcut")){
    if( !currentline.BeginsWith("#") ){
      cut += currentline.Data();
    }
  }

  TEventList *elist = new TEventList("elist");

  C->Draw(">>elist",cut);
  
  g4sbs_tree *T = new g4sbs_tree(C);

  inputfile >> nbins_x >> xmin >> xmax;
  inputfile >> nbins_z >> zmin >> zmax;
  inputfile >> PB >> PT;
  inputfile >> ndays;
 
  inputfile >> pi0flag;
  
  chargesep_m = 0.25;
  chargesep_b = 0.0;

  inputfile >> chargesep_m >> chargesep_b;

  currentline.ReadLine( inputfile );
  ifstream sbsopticsfile_up(currentline.Data());
  currentline.ReadLine( inputfile );
  ifstream sbsopticsfile_down(currentline.Data());
  currentline.ReadLine( inputfile );
  ifstream bbopticsfile(currentline.Data());

  vector<vector<double> > SBScoeff_up(4);
  vector<vector<int> > SBSexpon_up(5);

  vector<vector<double> > SBScoeff_down(4);
  vector<vector<int> > SBSexpon_down(5);

  vector<vector<double> > BBcoeff(4);
  vector<vector<int> > BBexpon(5);

  int nparams_SBS_down=0;
  sbsopticsfile_down >> nparams_SBS_down;
  for(int i=0; i<5; i++){
    if( i < 4) SBScoeff_down[i].resize(nparams_SBS_down);
    SBSexpon_down[i].resize(nparams_SBS_down);
  }

  for(int ipar=0; ipar<nparams_SBS_down; ipar++){
    for(int icoeff=0; icoeff<4; icoeff++){
      sbsopticsfile_down >> SBScoeff_down[icoeff][ipar];
    }
    for(int iexpon=0; iexpon<5; iexpon++){
      sbsopticsfile_down >> SBSexpon_down[iexpon][ipar];
    }
  }

  int nparams_SBS_up=0;
  sbsopticsfile_up >> nparams_SBS_up;
  for(int i=0; i<5; i++){
    if( i < 4) SBScoeff_up[i].resize(nparams_SBS_up);
    SBSexpon_up[i].resize(nparams_SBS_up);
  }

  for(int ipar=0; ipar<nparams_SBS_up; ipar++){
    for(int icoeff=0; icoeff<4; icoeff++){
      sbsopticsfile_up >> SBScoeff_up[icoeff][ipar];
    }
    for(int iexpon=0; iexpon<5; iexpon++){
      sbsopticsfile_up >> SBSexpon_up[iexpon][ipar];
    }
  }

  int nparams_BB=0;
  bbopticsfile >> nparams_BB;
  for(int i=0; i<5; i++){
    if( i < 4) BBcoeff[i].resize(nparams_BB);
    BBexpon[i].resize(nparams_BB);
  }

  for(int ipar=0; ipar<nparams_BB; ipar++){
    for(int icoeff=0; icoeff<4; icoeff++){
      bbopticsfile >> BBcoeff[icoeff][ipar];
    }
    for(int iexpon=0; iexpon<5; iexpon++){
      bbopticsfile >> BBexpon[iexpon][ipar];
    }
  }
  
  for(int i=0; i<nfiles; i++){
    
    weights[i] /= nfiles;
    weights[i] *= ndays * 24.0 * 3600.0;
  }
  vector<vector<int> >    Nsim_xz(nbins_x);
  vector<vector<double> > nevents_xz(nbins_x);
  vector<vector<double> > nevents_neutron_xz(nbins_x);
  vector<vector<double> > nevents_proton_xz(nbins_x);
  vector<vector<double> > avg_Q2(nbins_x);
  vector<vector<double> > avg_xbj(nbins_x);
  vector<vector<double> > avg_z(nbins_x);
  vector<vector<double> > avg_pT(nbins_x);
  vector<vector<double> > avg_phi(nbins_x);
  vector<vector<double> > dA_LL_He(nbins_x);
  vector<vector<double> > dA1_He(nbins_x);
  vector<vector<double> > dA_LL_n(nbins_x);
  vector<vector<double> > dA1_n(nbins_x);
  vector<vector<double> > avg_Pkin(nbins_x); //kinematic factor to go from A_LL to A1

  for(int i=0; i<nbins_x; i++){
    Nsim_xz[i].resize(nbins_z);
    nevents_xz[i].resize(nbins_z);
    nevents_neutron_xz[i].resize(nbins_z);
    nevents_proton_xz[i].resize(nbins_z);
    avg_Q2[i].resize(nbins_z);
    avg_xbj[i].resize(nbins_z);
    avg_z[i].resize(nbins_z);
    avg_pT[i].resize(nbins_z);
    avg_phi[i].resize(nbins_z);
    avg_Pkin[i].resize(nbins_z);
    dA_LL_He[i].resize(nbins_z);
    dA1_He[i].resize(nbins_z);
    dA_LL_n[i].resize(nbins_z);
    dA1_n[i].resize(nbins_z);
    for(int j=0; j<nbins_z; j++){
      Nsim_xz[i][j] = 0.0;
      nevents_xz[i][j] = 0.0;
      nevents_neutron_xz[i][j] = 0.0;
      nevents_proton_xz[i][j] = 0.0;
      avg_Q2[i][j] = 0.0;
      avg_xbj[i][j] = 0.0;
      avg_z[i][j] = 0.0;
      avg_pT[i][j] = 0.0;
      avg_phi[i][j] = 0.0;
      avg_Pkin[i][j] = 0.0;
    }
  }
  
  fout->cd();

  TH1D *hQ2 = new TH1D("hQ2", "", 200, 0.0, 12.0 );
  TH1D *hxbj = new TH1D("hxbj", "", 200, 0.0, 1.0 );
  TH1D *hz   = new TH1D("hz", "", 200, 0.0, 1.0 );
  TH1D *hPhperp = new TH1D("hPhperp", "", 200, 0.0, 2.0 );
  TH1D *hphih   = new TH1D("hphih", "", 200, -PI, PI );

  TH1D *hmpi0 = new TH1D("hmpi0","",200,0.0, 0.3);
  
  TH2D *hQ2_xbj = new TH2D("hQ2_xbj","",100,0.0, 1.0, 100, 0.0, 12.0 );
  TH2D *hz_xbj  = new TH2D("hz_xbj","",100,0.0,1.0,100,0.0,1.0);
  TH2D *hPhperp_xbj = new TH2D("hPhperp_xbj","",100,0.0,1.0,100,0.0,2.0);
  TH2D *hPhperp_z   = new TH2D("hPhperp_z","",100,0.0,1.0,100,0.0,2.0);
  TH2D *hphih_xbj   = new TH2D("hphih_xbj","",100,0.0,1.0,100,-PI,PI);
  TH2D *hphih_z     = new TH2D("hphih_z","",100,0.0,1.0,100,-PI,PI);
  TH2D *hphih_Phperp = new TH2D("hphih_Phperp","",100,0.0,1.6,100,-PI,PI);
  
  TH1D *hevzdiff = new TH1D("hevzdiff","",200,-0.2,0.2);
  TH1D *hhvzdiff = new TH1D("hhvzdiff","",200,-0.2,0.2);

  TH1D* hQ2diff = new TH1D("hQ2diff","",200,-1.0,1.0);
  TH1D* hxbjdiff = new TH1D("hxbjdiff","",200,-0.3,0.3);
  TH1D* hzdiff   = new TH1D("hzdiff","",200,-0.1,0.1);
  TH1D* hPhperpdiff = new TH1D("hPhperpdiff","",200,-0.2,0.2);
  TH1D* hphihdiff = new TH1D("hphihdiff","",200,-0.5,0.5);

  TH2D *hphihcorr = new TH2D("hphihcorr","",200,-PI,PI,200,-PI,PI);

  TH2D *hphihdiff_pT = new TH2D("hphihdiff_pT","",200,0.0,1.6,200,-0.5,0.5);

  double exfp, eyfp, expfp, eypfp, eptrack;
  double hxfp, hyfp, hxpfp, hypfp, hptrack;
  
  double extar, eytar, exptar, eyptar, eprecon;
  double hxtar, hytar, hxptar, hyptar, hprecon;

  TVector3 xaxis(1,0,0);
  TVector3 yaxis(0,1,0);
  TVector3 zaxis(0,0,1);

  long nevent = 0;
  
  while( T->GetEntry( elist->GetEntry(nevent++)) ){
    if( nevent%1000 == 0 ){
      cout << nevent << ", tree number = " << T->fChain->GetTreeNumber() << endl;
    }
    
    double Ebeam = T->gen_Ebeam;
    
    double bbtheta = T->gen_thbb;
    double sbstheta = T->gen_thhcal;

    TVector3 BB_zaxis( -sin(bbtheta), 0, cos(bbtheta) );
    TVector3 BB_xaxis( 0, -1, 0 );
    TVector3 BB_yaxis = (BB_zaxis.Cross(BB_xaxis)).Unit();
    
    TVector3 SBS_zaxis( sin(sbstheta), 0, cos(sbstheta) );
    TVector3 SBS_xaxis( 0, -1, 0 );
    TVector3 SBS_yaxis = (SBS_zaxis.Cross(SBS_xaxis)).Unit();
    
    TVector3 vertex( T->ev_vx, T->ev_vy, T->ev_vz );
    TVector3 vspec_BB(vertex.Dot( BB_xaxis ), vertex.Dot( BB_yaxis ), vertex.Dot( BB_zaxis ) );
    TVector3 vspec_SBS(vertex.Dot( SBS_xaxis ), vertex.Dot( SBS_yaxis ), vertex.Dot( SBS_zaxis ) );
    //Let's do reconstructed kinematics: 
    //Two iterations for electron and hadron arms:
    extar = -T->ev_vy;
    hxtar = -T->ev_vy;

    double mpi = 0.13957018; 
    double mpi0 = 0.1349766;
    double mK   = 0.493677;
    double mp   = 0.938272046;
    double mn   = 0.939565379;

    double Mh = mpi; //default to mpi;

    bool good_earm=false, good_harm=false;
    for(int itrack=0; itrack<T->ntracks; itrack++){
      if( (*(T->trackerid))[itrack] == 1 ){
	good_harm = true;
	hxfp = (*(T->trackxfit))[itrack];
	hyfp = (*(T->trackyfit))[itrack];
	hxpfp = (*(T->trackxpfit))[itrack];
	hypfp = (*(T->trackypfit))[itrack];
	hptrack = (*(T->trackp))[itrack];

	//determine up-bending or down-bending track:
	
	
	int pid = (*(T->trackpid))[itrack];

	if( abs(pid) == 211 ){ //charged pions:
	  Mh = mpi;
	} else if( abs(pid) == 111 ){ //pi0:
	  Mh = mpi0;
	} else if( abs(pid) == 321 ){ //K+/-
	  Mh = mK;
	} else if( abs(pid) == 2212 ){ //proton:
	  Mh = mp;
	} else if( abs(pid) == 2112 ){ //neutron:
	  Mh = mn;
	}
	double recon_sum[4] = {0,0,0,0};
	double pinv;
	
	//hxtar = 0.0;

	for(int iter=0; iter<2; iter++){
	  for(int coeff=0; coeff<4; coeff++){
	    recon_sum[coeff] = 0.0;
	  }

	  if( hxpfp < chargesep_b + chargesep_m * hxfp ){ //up-bending:

	    for(int par=0; par<nparams_SBS_up; par++){
	      for(int coeff=0; coeff<4; coeff++){
		//is this an up-bending or downbending track?
	      
		recon_sum[coeff] += SBScoeff_up[coeff][par] * pow( hxfp, SBSexpon_up[0][par] ) * pow( hyfp, SBSexpon_up[1][par] ) * pow( hxpfp, SBSexpon_up[2][par] ) *
		  pow( hypfp, SBSexpon_up[3][par] ) * pow( hxtar, SBSexpon_up[4][par] );
	      }
	    }
	  } else { //down-bending:
	    for(int par=0; par<nparams_SBS_down; par++){
	      for(int coeff=0; coeff<4; coeff++){
		//is this an up-bending or downbending track?
	      
		recon_sum[coeff] += SBScoeff_down[coeff][par] * pow( hxfp, SBSexpon_down[0][par] ) * pow( hyfp, SBSexpon_down[1][par] ) * pow( hxpfp, SBSexpon_down[2][par] ) *
		  pow( hypfp, SBSexpon_down[3][par] ) * pow( hxtar, SBSexpon_down[4][par] );
	      }
	    }
	  }
	  
	  hxptar = recon_sum[0];
	  hyptar = recon_sum[1];
	  hytar = recon_sum[2];
	  hprecon = 1.0/recon_sum[3];
	  
	  double vztemp = -hytar / ( sin(sbstheta) + hyptar*cos(sbstheta) );

	  //	  cout << "iteration, vzrecon, vztrue = " << iter << ", " << vztemp << ", " << T->ev_vz << endl;

	  if( iter == 1 ){
	    hhvzdiff->Fill( vztemp - T->ev_vz );
	  }

	  hxtar = vspec_SBS.X() - vztemp * cos( sbstheta ) * hxptar;
	}
      }
      if( (*(T->trackerid))[itrack] == 0 ){
	good_earm = true;
	exfp = (*(T->trackxfit))[itrack];
	eyfp = (*(T->trackyfit))[itrack];
	expfp = (*(T->trackxpfit))[itrack];
	eypfp = (*(T->trackypfit))[itrack];
	eptrack = (*(T->trackp))[itrack];

	double recon_sum[4] = {0,0,0,0};
	double pinv;
	
	extar = 0.0;

	for(int iter=0; iter<2; iter++){
	  for(int coeff=0; coeff<4; coeff++){
	    recon_sum[coeff] = 0.0;
	  }
	  
	  for(int par=0; par<nparams_BB; par++){
	    for(int coeff=0; coeff<4; coeff++){
	      recon_sum[coeff] += BBcoeff[coeff][par] * pow( exfp, BBexpon[0][par] ) * pow( eyfp, BBexpon[1][par] ) * pow( expfp, BBexpon[2][par] ) *
		pow( eypfp, BBexpon[3][par] ) * pow( extar, BBexpon[4][par] );
	    }
	  }
	  
	  exptar = recon_sum[0];
	  eyptar = recon_sum[1];
	  eytar = recon_sum[2];
	  eprecon = 1.0/recon_sum[3];
	  
	  double vztemp = eytar / ( sin(bbtheta) - eyptar*cos(bbtheta) );

	  // cout << "iteration, vzrecon, vztrue = " << iter << ", " << vztemp << ", " << T->ev_vz << endl;

	  if( iter == 1 ){
	    hevzdiff->Fill( vztemp - T->ev_vz );
	  }

	  extar = vspec_BB.X() - vztemp * cos( bbtheta ) * exptar;
	}
	
      }
    }
    
    TLorentzVector TwoPhotons;

    int pi0accept = 1;
    
    if( pi0flag != 0 ){ //Check if we have acceptance for pi0s: require two photons from pi0 decay to hit the calorimeter; AND require a minimum of 1-pixel 
      // separation between the two hits:
      // We also want reconstructed quantities:
      TVector3 vertex_earm( T->ev_vx, T->ev_vy, eytar/( sin(bbtheta) - eyptar*cos(bbtheta) ) ); //This includes resolution!

      pi0accept = 0;
      
      int ngoodhits = 0;
      
      vector<int> goodhits;
      vector<int> ixgood, iygood;
      vector<double> xgood, ygood, egood, tgood;
      vector<double> xcgood, ycgood;
      
      for(int hit=0; hit<T->hc_ndata; hit++){
	if( (T->hc_mid)[hit] == 2 && (T->hc_e)[hit] >= 0.5 ){
	  
	  goodhits.push_back(hit);
	  
	  ixgood.push_back( int(  ( (T->hc_x)[hit] + Hcal_width/2.0)/pixelsize  ) );
	  iygood.push_back( int(  ( (T->hc_y)[hit] + Hcal_height/2.0)/pixelsize  ) );
	  
	  xgood.push_back( (T->hc_x)[hit] );
	  ygood.push_back( (T->hc_y)[hit] );
	  xcgood.push_back( -Hcal_width/2.0 + (ixgood[ngoodhits]+0.5)*pixelsize );
	  ycgood.push_back( -Hcal_height/2.0 + (iygood[ngoodhits]+0.5)*pixelsize );
	  egood.push_back( (T->hc_e)[hit] );
	  tgood.push_back( (T->hc_t)[hit] );
	  
	  ngoodhits++;
	}
      }
      
      if( ngoodhits >= 2 ){
	
	double mindiff; //locate hit pair with minimum invariant mass difference:
	
	int best1=-1, best2=-1;
	
	for(int i=0; i<ngoodhits; i++){
	  for(int j=i+1; j<ngoodhits; j++){
	    TVector3 ray1( xgood[i]*cos(T->gen_thhcal) + T->gen_dhcal*sin(T->gen_thhcal) - T->ev_vx, 
			   ygood[i] - T->ev_vy,
			   -xgood[i]*sin(T->gen_thhcal) + T->gen_dhcal*cos(T->gen_thhcal) - T->ev_vz );
	    TVector3 ray2( xgood[j]*cos(T->gen_thhcal) + T->gen_dhcal*sin(T->gen_thhcal) - T->ev_vx, 
			   ygood[j] - T->ev_vy,
			   -xgood[j]*sin(T->gen_thhcal) + T->gen_dhcal*cos(T->gen_thhcal) - T->ev_vz );
	    
	    TLorentzVector p1( ray1.Unit()*egood[i], egood[i] );
	    TLorentzVector p2( ray2.Unit()*egood[j], egood[j] );
	    
	    double mass = (p1+p2).M();
	    
	    //cout << "i,j,mass = " << i << ", " << j << ", " << mass << endl;
	    
	    if( (i == 0 && j == 1 ) || fabs(mass - mpi0) < mindiff ){
	      mindiff = fabs(mass - mpi0);
	      best1 = i;
	      best2 = j;
	    }
	  }
	}
	
	if( best1 >= 0 && best2 > best1 ){
	  if( abs(ixgood[best1] - ixgood[best2]) > 1 ||  
	      abs(iygood[best1] - iygood[best2]) > 1 ){

	    TVector3 ray1recon( xcgood[best1]*cos(sbstheta) + T->gen_dhcal*sin(sbstheta) - vertex_earm.X(), 
				ycgood[best1] - vertex_earm.Y(), 
				-xcgood[best1]*sin(sbstheta) + T->gen_dhcal*cos(sbstheta) - vertex_earm.Z() );
	    
	    TVector3 ray2recon( xcgood[best2]*cos(sbstheta) + T->gen_dhcal*sin(sbstheta) - vertex_earm.X(), 
				ycgood[best2] - vertex_earm.Y(), 
				-xcgood[best2]*sin(sbstheta) + T->gen_dhcal*cos(sbstheta) - vertex_earm.Z() );

	    double E1recon = egood[best1] * (1.0 + num.Gaus(0.0,0.14/sqrt(egood[best1])));
	    double E2recon = egood[best2] * (1.0 + num.Gaus(0.0,0.14/sqrt(egood[best2])));

	    TLorentzVector photon1( ray1recon.Unit()*E1recon, E1recon );
	    TLorentzVector photon2( ray2recon.Unit()*E2recon, E2recon );

	    TwoPhotons = photon1 + photon2;

	    // cout << "True pi0:" << endl;
	    // TLorentzVector TruePi0( T->ev_npx, T->ev_npy, T->ev_npz, sqrt(pow(T->ev_np,2)+pow(mpi0,2)) );
	    // TruePi0.Print();
	    // cout << "Reconstructed pi0:" << endl;
	    // TwoPhotons.Print();

	    pi0accept = 1;
	  }
	}
      }
    }
    
    if( good_earm && (good_harm || pi0accept) ){

      TVector3 ephat_spec( exptar, eyptar, 1.0 );
      ephat_spec = ephat_spec.Unit();
      TVector3 hphat_spec( hxptar, hyptar, 1.0 );
      hphat_spec = hphat_spec.Unit();
      TVector3 ephat = ephat_spec.X() * BB_xaxis + ephat_spec.Y() * BB_yaxis + ephat_spec.Z() * BB_zaxis;
      TVector3 hphat = hphat_spec.X() * SBS_xaxis + hphat_spec.Y() * SBS_yaxis + hphat_spec.Z() * SBS_zaxis;
      
      TVector3 epvect = eprecon * ephat;
      TVector3 hpvect = hprecon * hphat;
      
      double Eh = sqrt( pow( Mh,2 ) + pow( hprecon, 2 ) );

      TLorentzVector Kprime( epvect, eprecon );
      TLorentzVector Ph( hpvect, Eh );
      
      if( pi0flag != 0 && pi0accept == 1 ) Ph = TwoPhotons;

      TLorentzVector K( 0, 0, Ebeam, Ebeam );
      TLorentzVector qrecon = K - Kprime;
      TLorentzVector P( 0, 0, 0, mp ); //use proton mass even though scattering on a neutron:
      
      double Q2_recon = -qrecon.M2();
      double xbj_recon = Q2_recon/(2.0*P.Dot(qrecon));
      double z_recon = P.Dot(Ph)/P.Dot(qrecon);

      TVector3 Ph_perp_vect = Ph.Vect() - (Ph.Vect() ).Dot( qrecon.Vect() ) / ( (qrecon.Vect()).Mag2() )  * qrecon.Vect();

      double Ph_perp_recon = Ph_perp_vect.Mag();
      
      TVector3 reaction_plane_yaxis = (K.Vect().Cross(Kprime.Vect())).Unit();
      TVector3 reaction_plane_xaxis = (reaction_plane_yaxis.Cross(qrecon.Vect())).Unit();
      
      double phih_recon = atan2( Ph_perp_vect.Dot( reaction_plane_yaxis ), Ph_perp_vect.Dot( reaction_plane_xaxis ) );

      double rate = T->ev_rate;
      double etheta = T->ev_th;
      double xbj = T->ev_xbj;
      double Q2  = T->ev_Q2;
      double W2  = T->ev_W2;
      double z   = T->ev_z;
      double pT  = T->ev_phperp;
      double phih = T->ev_phih;

      double Mp = 0.938272;
      double nu = Q2/(2.0*Mp*xbj);
    
      double y = nu/Ebeam;
      etheta = acos(TMath::Max(-1.0,TMath::Min(1.0,1.0 - Q2/(2.0*pow(Ebeam,2)*(1.0-y)))));

      double gamma = 2.0*Mp*xbj/sqrt(Q2);
      double epsilon = pow(1.0 + 2.0*(1.0 + pow(nu,2)/Q2)*pow(tan(etheta/2.0),2), -1.0 );

      double dR;
    
      double R = get_R(xbj, Q2, dR );
    
      double D = (1.0 - (1.0-y)*epsilon)/(1.0 + epsilon*R);
      double eta = 2.0*gamma*(1.0-y)/(2.0-y);
      double Pkin = D*(1.0+gamma*eta)*(1.0+R)/(1.0+pow(gamma,2));
    
      // x = Q2/2P dot q --> P dot q = Q2/2x = M nu 
      // W2 = (P+q)^2 = M^2 -Q2 + 2Pdotq
      // y = nu/Ebeam = P dot q / P dot k 
   
      int bin_x = int( (xbj-xmin)*double(nbins_x)/(xmax-xmin) );
      int bin_z = int( (z - zmin)*double(nbins_z)/(zmax-zmin) );

      

      if( bin_x >= 0 && bin_x < nbins_x && 
	  bin_z >= 0 && bin_z < nbins_z && 
	  pi0accept != 0 ){
	Nsim_xz[bin_x][bin_z] += 1;
	nevents_xz[bin_x][bin_z] += rate*weights[T->fChain->GetTreeNumber()];
	if( T->ev_nucl == 0 ){ 
	  nevents_neutron_xz[bin_x][bin_z] += rate*weights[T->fChain->GetTreeNumber()];
	} else {
	  nevents_proton_xz[bin_x][bin_z] += rate*weights[T->fChain->GetTreeNumber()];
	}
	avg_Q2[bin_x][bin_z] += Q2*rate*weights[T->fChain->GetTreeNumber()];
	avg_xbj[bin_x][bin_z] += xbj*rate*weights[T->fChain->GetTreeNumber()];
	avg_z[bin_x][bin_z] += z*rate*weights[T->fChain->GetTreeNumber()];
	avg_pT[bin_x][bin_z] += pT*rate*weights[T->fChain->GetTreeNumber()];
	avg_phi[bin_x][bin_z] += phih*rate*weights[T->fChain->GetTreeNumber()];
	avg_Pkin[bin_x][bin_z] += Pkin*rate*weights[T->fChain->GetTreeNumber()];
      }

      double Weight = rate*weights[T->fChain->GetTreeNumber()];

      hQ2->Fill( Q2, Weight );
      hxbj->Fill( xbj, Weight );
      hz->Fill( z, Weight );
      hPhperp->Fill(pT, Weight );
      hphih->Fill( phih, Weight );

      if( pi0flag != 0 ){
	hmpi0->Fill( TwoPhotons.M(), Weight );
      }

      hQ2_xbj->Fill( xbj, Q2, Weight );
      hz_xbj->Fill( xbj, z, Weight );
      hPhperp_xbj->Fill( xbj, pT, Weight );
      hPhperp_z->Fill( z, pT, Weight );
      hphih_xbj->Fill( xbj, phih, Weight );
      hphih_z->Fill( z, phih, Weight );
      hphih_Phperp->Fill( pT, phih, Weight );

      hQ2diff->Fill( Q2_recon - Q2 );
      hxbjdiff->Fill( xbj_recon - xbj );
      hzdiff->Fill( z_recon - z );
      hPhperpdiff->Fill( Ph_perp_recon - pT );
      hphihdiff->Fill( phih_recon - phih );

      hphihcorr->Fill( phih, phih_recon );

      hphihdiff_pT->Fill( pT, phih_recon - phih );
      
    }
  }
  
  char header[255];
  sprintf( header, "%5s %5s %5s %5s %8s %8s %8s %8s %12s %8s %12s %12s %12s %12s", "#xmin", "xmax", "zmin", "zmax", "<xbj>", "<Q2>", "<z>", 
	   "<pT>", "Nevent", "fD", "dALL_He3", "dALL_n", "dA1_He3", "dA1_n" );

  outputfile << header << endl;

  for(int i=0; i<nbins_x; i++){
    for(int j=0; j<nbins_z; j++){
      char cline[255];
      double xlow_bin = xmin + i*(xmax-xmin)/double(nbins_x);
      double xhigh_bin = xlow_bin + (xmax-xmin)/double(nbins_x);
      double zlow_bin = zmin + j*(zmax-zmin)/double(nbins_z);
      double zhigh_bin = zlow_bin + (zmax-zmin)/double(nbins_z);
      double Q2avg_bin = avg_Q2[i][j]/nevents_xz[i][j];
      double xavg_bin = avg_xbj[i][j]/nevents_xz[i][j];
      double zavg_bin = avg_z[i][j]/nevents_xz[i][j];
      double pTavg_bin = avg_pT[i][j]/nevents_xz[i][j];
      double Pkinavg_bin = avg_Pkin[i][j]/nevents_xz[i][j];

      double dilution = nevents_neutron_xz[i][j]/nevents_xz[i][j];

      dA_LL_He[i][j] = 1.0/(PB*PT*sqrt(nevents_xz[i][j]));
      dA_LL_n[i][j] = dA_LL_He[i][j]/(0.86*dilution);
      dA1_He[i][j] = dA_LL_He[i][j]/Pkinavg_bin;
      dA1_n[i][j] = dA_LL_n[i][j]/Pkinavg_bin;
      
      sprintf(cline, "%5.3g %5.3g %5.3g %5.3g %8.5g %8.5g %8.5g %8.5g %12.5g %8.5g %12.5g %12.5g %12.5g %12.5g", xlow_bin, xhigh_bin, zlow_bin, zhigh_bin, 
	      xavg_bin, Q2avg_bin, zavg_bin, pTavg_bin, nevents_xz[i][j], dilution, dA_LL_He[i][j], dA_LL_n[i][j], dA1_He[i][j], dA1_n[i][j] );
      if( Nsim_xz[i][j] >= 10 && dA1_n[i][j] <= 0.2 ) outputfile << cline << endl;
    }
  }

  fout->Write();
  fout->Close();
  
}
