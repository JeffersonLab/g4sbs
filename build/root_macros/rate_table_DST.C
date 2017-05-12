#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
//#include "g4sbs_tree.C"
#include "SIDIS_DST.C"
#include <iostream>
#include <fstream>
#include <vector>
#include "TSystem.h"
#include "TString.h"
#include "TEventList.h"
#include "TCut.h"
#include "TObject.h"
#include "TROOT.h"
//#include "G4SBSRunData.hh"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TRandom.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TClonesArray.h"
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

  int nbins_x, nbins_z, nbins_pT;
  double xmin, xmax, zmin, zmax, pTmin, pTmax;
  
  double PB,PT;

  double chargesep_m, chargesep_b;

  ifstream inputfile(inputfilename);
  ofstream outputfile(outputfilename);
  
  //We already know rate, the only thing we need is duration of the expt.
  // to get the number of events:
  double ndays;

  int pi0flag = 0;
 
  TChain *C = new TChain("Tout");

  int nfiles=0;
  vector<int> ngen;
  vector<int> ntries;
  vector<double> weights; //correct for ratio of generated events to number of tries:

  double sum_weights = 0.0;

  //  G4SBSRunData *rtemp;

  //TObject *rtemp;

  TString currentline;
  while(currentline.ReadLine(inputfile)&&!currentline.BeginsWith("endlist")){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      TString fname = ( (TObjString*) (*tokens)[0] )->GetString();
 
      C->Add(fname.Data());

      TString sndays = ( (TObjString*) (*tokens)[1] )->GetString();
      for(int ifile=nfiles; ifile<C->GetNtrees(); ifile++){
	weights.push_back( sndays.Atof() * 24.0 * 3600.0 );
      }
      nfiles = C->GetNtrees();
    }
  }

  TCut cut = "";

  while(currentline.ReadLine(inputfile)&&!currentline.BeginsWith("endcut")){
    if( !currentline.BeginsWith("#") ){
      cut += currentline.Data();
    }
  }

  inputfile >> nbins_x >> xmin >> xmax;
  inputfile >> nbins_z >> zmin >> zmax;
  inputfile >> nbins_pT >> pTmin >> pTmax;
  inputfile >> PB >> PT;

  TEventList *elist = new TEventList("elist");

  C->Draw(">>elist",cut);
  
  //g4sbs_tree *T = new g4sbs_tree(C);
  SIDIS_DST *T = new SIDIS_DST(C);

  
  //inputfile >> ndays;
 
  //inputfile >> pi0flag;
  
  // vector<vector<int> >    Nsim_xz(nbins_x);
  // vector<vector<double> > nevents_xz(nbins_x);
  // vector<vector<double> > nevents_neutron_xz(nbins_x);
  // vector<vector<double> > nevents_proton_xz(nbins_x);
  // vector<vector<double> > avg_Q2(nbins_x);
  // vector<vector<double> > avg_xbj(nbins_x);
  // vector<vector<double> > avg_z(nbins_x);
  // vector<vector<double> > avg_pT(nbins_x);
  // vector<vector<double> > avg_phi(nbins_x);
  // vector<vector<double> > dA_LL_He(nbins_x);
  // vector<vector<double> > dA1_He(nbins_x);
  // vector<vector<double> > dA_LL_n(nbins_x);
  // vector<vector<double> > dA1_n(nbins_x);
  // vector<vector<double> > avg_Pkin(nbins_x); //kinematic factor to go from A_LL to A1

  int Nsim[nbins_x][nbins_z][nbins_pT];
  double nevents[nbins_x][nbins_z][nbins_pT];
  double nevents_neutron[nbins_x][nbins_z][nbins_pT];
  double nevents_proton[nbins_x][nbins_z][nbins_pT];
  double avg_Q2[nbins_x][nbins_z][nbins_pT];
  double avg_xbj[nbins_x][nbins_z][nbins_pT];
  double avg_z[nbins_x][nbins_z][nbins_pT];
  double avg_pT[nbins_x][nbins_z][nbins_pT];
  double avg_phi[nbins_x][nbins_z][nbins_pT];
  double dA_LL_He[nbins_x][nbins_z][nbins_pT];
  double dA1_He[nbins_x][nbins_z][nbins_pT];
  double dA_LL_n[nbins_x][nbins_z][nbins_pT];
  double dA1_n[nbins_x][nbins_z][nbins_pT];
  double avg_Pkin[nbins_x][nbins_z][nbins_pT];

  for(int i=0; i<nbins_x; i++){
    for(int j=0; j<nbins_z; j++){
      for(int k=0; k<nbins_pT; k++){
	Nsim[i][j][k] = 0;
	nevents[i][j][k] = 0.0;
	nevents_neutron[i][j][k] = 0.0;
	nevents_proton[i][j][k] = 0.0;
	avg_Q2[i][j][k] = 0.0;
	avg_xbj[i][j][k] = 0.0;
	avg_z[i][j][k] = 0.0;
	avg_pT[i][j][k] = 0.0;
	avg_phi[i][j][k] = 0.0;
	dA_LL_He[i][j][k] = 0.0;
	dA1_He[i][j][k] = 0.0;
	dA_LL_n[i][j][k] = 0.0;
	dA1_n[i][j][k] = 0.0;
	avg_Pkin[i][j][k] = 0.0;
      }
    }
  }

  // for(int i=0; i<nbins_x; i++){
  //   Nsim_xz[i].resize(nbins_z);
  //   nevents_xz[i].resize(nbins_z);
  //   nevents_neutron_xz[i].resize(nbins_z);
  //   nevents_proton_xz[i].resize(nbins_z);
  //   avg_Q2[i].resize(nbins_z);
  //   avg_xbj[i].resize(nbins_z);
  //   avg_z[i].resize(nbins_z);
  //   avg_pT[i].resize(nbins_z);
  //   avg_phi[i].resize(nbins_z);
  //   avg_Pkin[i].resize(nbins_z);
  //   dA_LL_He[i].resize(nbins_z);
  //   dA1_He[i].resize(nbins_z);
  //   dA_LL_n[i].resize(nbins_z);
  //   dA1_n[i].resize(nbins_z);
  //   for(int j=0; j<nbins_z; j++){
  //     Nsim_xz[i][j] = 0.0;
  //     nevents_xz[i][j] = 0.0;
  //     nevents_neutron_xz[i][j] = 0.0;
  //     nevents_proton_xz[i][j] = 0.0;
  //     avg_Q2[i][j] = 0.0;
  //     avg_xbj[i][j] = 0.0;
  //     avg_z[i][j] = 0.0;
  //     avg_pT[i][j] = 0.0;
  //     avg_phi[i][j] = 0.0;
  //     avg_Pkin[i][j] = 0.0;
  //   }
  // }
  
  fout->cd();

  TClonesArray *histos_pT_phih_xz = new TClonesArray("TH2D",nbins_x*nbins_z);

  for(int i=0; i<nbins_x; i++){
    for(int j=0; j<nbins_z; j++){
      TString histname;
      histname.Form("hpT_phih_x%d_z%d",i,j);

      new( (*histos_pT_phih_xz)[j+i*nbins_z] ) TH2D(histname.Data(),"",200,-PI,PI,200,0.0,2.0);
    }
  }

  TH1D *hQ2 = new TH1D("hQ2", "", 200, 0.0, 15.0 );
  TH1D *hxbj = new TH1D("hxbj", "", 200, 0.0, 1.0 );
  TH1D *hz   = new TH1D("hz", "", 200, 0.0, 1.0 );
  TH1D *hPhperp = new TH1D("hPhperp", "", 200, 0.0, 2.0 );
  TH1D *hphih   = new TH1D("hphih", "", 200, -PI, PI );

  TH1D *hmpi0 = new TH1D("hmpi0","",200,0.0, 0.3);
  
  TH2D *hQ2_xbj = new TH2D("hQ2_xbj","",200,0.0, 1.0, 200, 0.0, 15.0 );
  TH2D *hz_xbj  = new TH2D("hz_xbj","",200,0.0,1.0,200,0.0,1.0);
  TH2D *hPhperp_xbj = new TH2D("hPhperp_xbj","",200,0.0,1.0,200,0.0,2.0);
  TH2D *hPhperp_z   = new TH2D("hPhperp_z","",200,0.0,1.0,200,0.0,2.0);
  TH2D *hphih_xbj   = new TH2D("hphih_xbj","",200,0.0,1.0,200,-PI,PI);
  TH2D *hphih_z     = new TH2D("hphih_z","",200,0.0,1.0,200,-PI,PI);
  TH2D *hphih_Phperp = new TH2D("hphih_Phperp","",200,0.0,1.6,200,-PI,PI);
  
  TH2D *hQ2_xbj_nosmear = new TH2D("hQ2_xbj_nosmear","",200,0.0, 1.0, 200, 0.0, 15.0 );
  TH2D *hz_xbj_nosmear  = new TH2D("hz_xbj_nosmear","",200,0.0,1.0,200,0.0,1.0);
  TH2D *hPhperp_xbj_nosmear = new TH2D("hPhperp_xbj_nosmear","",200,0.0,1.0,200,0.0,2.0);
  TH2D *hPhperp_z_nosmear   = new TH2D("hPhperp_z_nosmear","",200,0.0,1.0,200,0.0,2.0);
  TH2D *hphih_xbj_nosmear   = new TH2D("hphih_xbj_nosmear","",200,0.0,1.0,200,-PI,PI);
  TH2D *hphih_z_nosmear     = new TH2D("hphih_z_nosmear","",200,0.0,1.0,200,-PI,PI);
  TH2D *hphih_Phperp_nosmear = new TH2D("hphih_Phperp_nosmear","",200,0.0,1.6,200,-PI,PI);

  TH1D *hevzdiff = new TH1D("hevzdiff","",200,-0.2,0.2);
  TH1D *hhvzdiff = new TH1D("hhvzdiff","",200,-0.2,0.2);
  TH1D *hehvzdiff = new TH1D("hehvzdiff","",200,-0.2,0.2);

  TH1D *hepdiff = new TH1D("hepdiff","",200,-0.15,0.15);
  TH1D *hhpdiff = new TH1D("hhpdiff","",200,-0.1,0.1);

  TH1D* hQ2diff = new TH1D("hQ2diff","",200,-1.0,1.0);
  TH1D* hxbjdiff = new TH1D("hxbjdiff","",200,-0.2,0.2);
  TH1D* hzdiff   = new TH1D("hzdiff","",200,-0.1,0.1);
  TH1D* hPhperpdiff = new TH1D("hPhperpdiff","",200,-0.2,0.2);
  TH1D* hphihdiff = new TH1D("hphihdiff","",200,-0.5,0.5);

  TH1D *hQ2diff_nofermi = new TH1D("hQ2diff_nofermi","",200,-1.0,1.0);
  TH1D *hxbjdiff_nofermi = new TH1D("hxbjdiff_nofermi","",200,-0.2,0.2);
  TH1D *hzdiff_nofermi = new TH1D("hzdiff_nofermi","",200,-0.1,0.1);
  TH1D *hPhperpdiff_nofermi = new TH1D("hPhperpdiff_nofermi","",200,-0.2,0.2);
  TH1D *hphihdiff_nofermi = new TH1D("hphihdiff_nofermi","",200,-0.5,0.5);

  TH2D *hphihcorr = new TH2D("hphihcorr","",200,-PI,PI,200,-PI,PI);
  TH2D *hphihdiff_pT = new TH2D("hphihdiff_pT","",200,0.0,1.6,200,-0.5,0.5);


  long nevent = 0;
  
  while( T->GetEntry( elist->GetEntry(nevent++)) ){
    if( nevent%1000 == 0 ){
      cout << nevent << ", tree number = " << T->fChain->GetTreeNumber() << endl;
    }
    
    double Ebeam = T->Ebeam;
    
    double bbtheta = T->ThetaBB;
    double sbstheta = T->ThetaSBS;

    double mpi = 0.13957018; 
    double mpi0 = 0.1349766;
    double mK   = 0.493677;
    double mp   = 0.938272046;
    double mn   = 0.939565379;

    double Mh = mpi; //default to mpi;
	
    int pid = T->Hadron;

    if( abs(pid) == 1 ){ //charged pions:
      Mh = mpi;
    } else if( abs(pid) == 0 ){ //pi0:
      Mh = mpi0;
    } else if( abs(pid) == 2 ){ //K+/-
      Mh = mK;
    } else if( abs(pid) == 3 ){ //proton:
      Mh = mp;
    } else if( pid == 4 ){ //neutron: not currently implemented in G4SBS SIDIS generator
      Mh = mn;
    }

    double rate = T->Weight;
    double etheta = T->eth;
    double xbj = T->x;
    double Q2  = T->Q2;
    double W2  = T->W2;
    double z   = T->z;
    double pT  = T->pT;
    double phih = T->phih;

    double Mp = 0.938272;
    double nu = Q2/(2.0*Mp*xbj);
    
    double y = T->y;
    //etheta = acos(TMath::Max(-1.0,TMath::Min(1.0,1.0 - Q2/(2.0*pow(Ebeam,2)*(1.0-y)))));

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
    int bin_pT = int( (pT - pTmin)*double(nbins_pT)/(pTmax-pTmin) );

    
      
    if( bin_x >= 0 && bin_x < nbins_x && 
	bin_z >= 0 && bin_z < nbins_z && 
	bin_pT >= 0 && bin_pT < nbins_pT ){
      

      Nsim[bin_x][bin_z][bin_pT] += 1;
      nevents[bin_x][bin_z][bin_pT] += rate*weights[T->fChain->GetTreeNumber()];
      if( T->Nucleon == 0 ){ 
	nevents_neutron[bin_x][bin_z][bin_pT] += rate*weights[T->fChain->GetTreeNumber()];
      } else {
	nevents_proton[bin_x][bin_z][bin_pT] += rate*weights[T->fChain->GetTreeNumber()];
      }
      avg_Q2[bin_x][bin_z][bin_pT] += Q2*rate*weights[T->fChain->GetTreeNumber()];
      avg_xbj[bin_x][bin_z][bin_pT] += xbj*rate*weights[T->fChain->GetTreeNumber()];
      avg_z[bin_x][bin_z][bin_pT] += z*rate*weights[T->fChain->GetTreeNumber()];
      avg_pT[bin_x][bin_z][bin_pT] += pT*rate*weights[T->fChain->GetTreeNumber()];
      avg_phi[bin_x][bin_z][bin_pT] += phih*rate*weights[T->fChain->GetTreeNumber()];
      avg_Pkin[bin_x][bin_z][bin_pT] += Pkin*rate*weights[T->fChain->GetTreeNumber()];
    }
    
    double Weight = rate*weights[T->fChain->GetTreeNumber()];

    if( bin_x >= 0 && bin_x < nbins_x && 
	bin_z >= 0 && bin_z < nbins_z ){
      ( (TH2D*) (*histos_pT_phih_xz)[bin_z+nbins_z*bin_x] )->Fill(phih,pT,Weight);
    }

    hQ2->Fill( Q2, Weight );
    hxbj->Fill( xbj, Weight );
    hz->Fill( z, Weight );
    hPhperp->Fill(pT, Weight );
    hphih->Fill( phih, Weight );
    
    hmpi0->Fill( T->mpi0recon, Weight );

    hQ2_xbj->Fill( xbj, Q2, Weight );
    hz_xbj->Fill( xbj, z, Weight );
    hPhperp_xbj->Fill( xbj, pT, Weight );
    hPhperp_z->Fill( z, pT, Weight );
    hphih_xbj->Fill( xbj, phih, Weight );
    hphih_z->Fill( z, phih, Weight );
    hphih_Phperp->Fill( pT, phih, Weight );
    
    hevzdiff->Fill( T->evzrecon - T->vz, Weight );
    hhvzdiff->Fill( T->hvzrecon - T->vz, Weight );
    hehvzdiff->Fill( T->evzrecon - T->hvzrecon, Weight );

    hepdiff->Fill( T->eprecon/T->ep - 1.0, Weight );
    hhpdiff->Fill( T->hprecon/T->hp - 1.0, Weight );

    hQ2diff->Fill( T->Q2recon - Q2, Weight );
    hxbjdiff->Fill( T->xrecon - xbj, Weight );
    hzdiff->Fill( T->zrecon - z, Weight );
    hPhperpdiff->Fill( T->pTrecon - pT, Weight );
    hphihdiff->Fill( T->phihrecon - phih, Weight );
    
    hphihcorr->Fill( phih, T->phihrecon, Weight );
    hphihdiff_pT->Fill( pT, T->phihrecon - phih, Weight );
    
    hQ2diff_nofermi->Fill( T->Q2recon - 2.0*T->Ebeam*T->ep*(1.0-cos(T->eth)), Weight );
    hxbjdiff_nofermi->Fill( T->xrecon - 2.0*T->Ebeam*T->ep*(1.0-cos(T->eth))/(2.0*mp*(Ebeam - T->ep) ), Weight );
    hzdiff_nofermi->Fill( T->zrecon - sqrt(pow(T->hp,2)+pow(Mh,2))/(Ebeam-T->ep), Weight );
    
    double Q2_nosmear = 2.0*T->Ebeam*T->ep*(1.0-cos(T->eth));
    double x_nosmear = Q2_nosmear/(2.0*mp*(Ebeam-T->ep));
    double z_nosmear = sqrt(pow(T->hp,2)+pow(Mh,2))/(Ebeam-T->ep); 

    double W2_nosmear = pow(mp,2) + Q2_nosmear*(1.0-x_nosmear)/x_nosmear;


    TLorentzVector ElectronPrime( TVector3( T->ep*sin(T->eth)*cos(T->eph), T->ep*sin(T->eth)*sin(T->eph), T->ep*cos(T->eth)), T->ep );
    TLorentzVector HadronPrime( TVector3( T->hp*sin(T->hth)*cos(T->hph), T->hp*sin(T->hth)*sin(T->hph), T->hp*cos(T->hth)), sqrt(pow(T->hp,2)+pow(Mh,2)) );
    TLorentzVector Beam( 0.0, 0.0, Ebeam, Ebeam );
    TLorentzVector NucleonRest( 0.0, 0.0, 0.0, mp );
    TLorentzVector q = Beam - ElectronPrime;
    TLorentzVector Xprime = NucleonRest + q - HadronPrime;
    
    double MX2_nosmear = Xprime.M2();

    //How to compute phih and Phperp?
    TVector3 lepton_plane_yaxis = (q.Vect().Cross(Beam.Vect())).Unit();
    TVector3 lepton_plane_zaxis = (q.Vect()).Unit();
    TVector3 lepton_plane_xaxis = (lepton_plane_yaxis.Cross(lepton_plane_zaxis)).Unit();
    
    TVector3 HadronPerp = HadronPrime.Vect() - HadronPrime.Vect().Dot(lepton_plane_zaxis) * lepton_plane_zaxis;
    hPhperpdiff_nofermi->Fill( T->pTrecon - HadronPerp.Mag(), Weight );
    hphihdiff_nofermi->Fill( T->phihrecon - atan2( HadronPerp.Dot(lepton_plane_yaxis), HadronPerp.Dot(lepton_plane_xaxis) ), Weight );

    double phih_nosmear = atan2( HadronPerp.Dot(lepton_plane_yaxis), HadronPerp.Dot(lepton_plane_xaxis) );
    
    if( W2_nosmear >= 4.0 && MX2_nosmear >= 2.25 && Q2_nosmear >= 1.0 ){

      hQ2_xbj_nosmear->Fill( x_nosmear, Q2_nosmear, Weight );
      hz_xbj_nosmear->Fill( x_nosmear, z_nosmear, Weight );

      hPhperp_xbj_nosmear->Fill( x_nosmear, HadronPerp.Mag(), Weight );
      hPhperp_z_nosmear->Fill( z_nosmear, HadronPerp.Mag(), Weight );
      hphih_xbj_nosmear->Fill( x_nosmear, phih_nosmear, Weight );
      hphih_z_nosmear->Fill( z_nosmear, phih_nosmear, Weight );
      hphih_Phperp_nosmear->Fill( HadronPerp.Mag(), phih_nosmear, Weight );

    }
  }

  char header[255];
  sprintf( header, "%5s %5s %5s %5s %5s %5s %8s %8s %8s %8s %12s %8s %12s %12s %12s %12s", "#xmin", "xmax", "zmin", "zmax", "pTmin", "pTmax", "<xbj>", "<Q2>", "<z>", 
	   "<pT>", "Nevent", "fD", "dALL_He3", "dALL_n", "dA1_He3", "dA1_n" );
  
  outputfile << header << endl;

  for(int i=0; i<nbins_x; i++){
    for(int j=0; j<nbins_z; j++){
      for(int k=0; k<nbins_pT; k++){
	char cline[255];
	double xlow_bin = xmin + i*(xmax-xmin)/double(nbins_x);
	double xhigh_bin = xlow_bin + (xmax-xmin)/double(nbins_x);
	double zlow_bin = zmin + j*(zmax-zmin)/double(nbins_z);
	double zhigh_bin = zlow_bin + (zmax-zmin)/double(nbins_z);
	double pTlow_bin = pTmin + k*(pTmax-pTmin)/double(nbins_pT);
	double pThigh_bin = pTlow_bin + (pTmax-pTmin)/double(nbins_pT);
	double Q2avg_bin = avg_Q2[i][j][k]/nevents[i][j][k];
	double xavg_bin = avg_xbj[i][j][k]/nevents[i][j][k];
	double zavg_bin = avg_z[i][j][k]/nevents[i][j][k];
	double pTavg_bin = avg_pT[i][j][k]/nevents[i][j][k];
	double Pkinavg_bin = avg_Pkin[i][j][k]/nevents[i][j][k];
	
	double dilution = nevents_neutron[i][j][k]/nevents[i][j][k];

	dA_LL_He[i][j][k] = 1.0/(PB*PT*sqrt(nevents[i][j][k]));
	dA_LL_n[i][j][k] = dA_LL_He[i][j][k]/(0.86*dilution);
	dA1_He[i][j][k] = dA_LL_He[i][j][k]/Pkinavg_bin;
	dA1_n[i][j][k] = dA_LL_n[i][j][k]/Pkinavg_bin;
	
	sprintf(cline, "%5.3g %5.3g %5.3g %5.3g %5.3g %5.3g %8.5g %8.5g %8.5g %8.5g %12.5g %8.5g %12.5g %12.5g %12.5g %12.5g", xlow_bin, xhigh_bin, zlow_bin, zhigh_bin, pTlow_bin, pThigh_bin,
		xavg_bin, Q2avg_bin, zavg_bin, pTavg_bin, nevents[i][j][k], dilution, dA_LL_He[i][j][k], dA_LL_n[i][j][k], dA1_He[i][j][k], dA1_n[i][j][k] );
	if( Nsim[i][j][k] >= 10 && dA1_n[i][j][k] <= 0.2 ) outputfile << cline << endl;
      }
    }
  }
  
  elist->Delete();

  //histos_pT_phih_xz->Write();

  fout->Write();
  fout->Close();
  
}
