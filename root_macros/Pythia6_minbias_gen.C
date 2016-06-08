#include "TPythia6.h"
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMCParticle.h"
#include "TGraph.h"

using namespace std;

const double Mp = 0.938272046;
const double me = 0.511e-3;

void Pythia6_minbias_gen( const char *outputfilename, double Ebeam=11.0, long ngen=125000 ){
  TFile *Fout = new TFile(outputfilename,"RECREATE");
  
  //ifstream infile(configfilename);

  // double Ebeam = 11.0;
  // infile >> Ebeam;

  // long ngen = 10000;
  // infile >> ngen;
  
  double pbeam = sqrt(pow(Ebeam,2)-pow(me,2));
  
  TVector3 k(0,0,pbeam);
  
  TLorentzVector kbeam( k, Ebeam );
  TVector3 Ptarget(0,0,0);
  TLorentzVector Pproton(Ptarget,Mp);

  TLorentzVector Ptotal = kbeam + Pproton;
  double Ecm = Ptotal.M();

  if( Ebeam < 1.0 ) return;
  
  TPythia6 Generator;
  Generator.SetPARP(2,1.0); //Set minimum CM energy to 1 GeV instead of the default 10 GeV:

  Generator.Initialize( "FIXT", "gamma/e-", "p", 11.0 );

  int Nparticles;
  float Q2,xbj,y,W2;
  vector<int> status;
  vector<int> pid;
  vector<int> parent;
  vector<int> fchild;
  vector<int> lchild;
  vector<double> px,py,pz,E,M,theta,phi,vx,vy,vz,t,tau;

  TTree *Tout = new TTree("Tout","Pythia6 min-bias events");
  Tout->Branch("Nparticles",&Nparticles,"Nparticles/I");
  Tout->Branch("Q2",&Q2,"Q2/F");
  Tout->Branch("xbj",&xbj,"xbj/F");
  Tout->Branch("y", &y, "y/F");
  Tout->Branch("W2", &W2, "W2/F");
  Tout->Branch("status",&status);
  Tout->Branch("pid",&pid);
  Tout->Branch("parent",&parent);
  Tout->Branch("fchild",&fchild);
  Tout->Branch("lchild",&lchild);
  Tout->Branch("px",&px);
  Tout->Branch("py",&py);
  Tout->Branch("pz",&pz);
  Tout->Branch("vx",&vx);
  Tout->Branch("vy",&vy);
  Tout->Branch("vz",&vz);
  Tout->Branch("E",&E);
  Tout->Branch("M",&M);
  Tout->Branch("theta",&theta);
  Tout->Branch("phi",&phi);
  Tout->Branch("t",&t);
  Tout->Branch("tau",&tau);

  vector<double> event_num;
  vector<double> sigma;

  TClonesArray *Particles = ( (TClonesArray*) Generator.GetListOfParticles() );
  
  for(long ievent=1; ievent<=ngen; ievent++ ){
    status.clear();
    pid.clear();
    parent.clear();
    fchild.clear();
    lchild.clear();
    px.clear();
    py.clear();
    pz.clear();
    E.clear();
    M.clear();
    theta.clear();
    phi.clear();
    vx.clear();
    vy.clear();
    vz.clear();
    t.clear();
    tau.clear();
    Generator.GenerateEvent();
    if( ievent%100 == 0 || ievent == ngen ){
      cout << "Event number " << ievent << endl;
      sigma.push_back( Generator.GetPARI(1) ); //total cross section in mb
      event_num.push_back( double(ievent) );
    }

    Nparticles = Particles->GetEntries();

    //0 is the beam electron, 1 is the target proton
    //2 is the scattered electron, and 3 is the virtual photon
    
    //compute global kinematic variables from header information:
    Q2 = pow( ( (TMCParticle*) (*Particles)[3] )->GetMass(), 2 );
    //xbj = Q2/2Pdot q = Q2/(2*Mp*Egamma);
    xbj = Q2/(2.0*( (TMCParticle*) (*Particles)[1] )->GetMass()*( (TMCParticle*) (*Particles)[3] )->GetEnergy() );
    //y = nu/Ebeam:
    y = ( (TMCParticle*) (*Particles)[3] )->GetEnergy()/( (TMCParticle*) (*Particles)[0] )->GetEnergy();
    //W2 = (P+q)^2 = Mp^2 - Q2 + 2Pdot q = Mp^2 -Q2 + 2Mp*Egamma
    W2 = pow(( (TMCParticle*) (*Particles)[1] )->GetMass(),2) - Q2 + 2.0*( (TMCParticle*) (*Particles)[1] )->GetMass()*( (TMCParticle*) (*Particles)[3] )->GetEnergy();
    
    for(int ipart=0; ipart<Nparticles; ipart++ ){
      status.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetKS() );
      pid.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetKF() );
      parent.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetParent() );
      fchild.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetFirstChild() );
      lchild.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetLastChild() );
      px.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetPx() );
      py.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetPy() );
      pz.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetPz() );
      vx.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetVx() );
      vy.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetVy() );
      vz.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetVz() );
      M.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetMass() );
      E.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetEnergy() );
      t.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetTime() );
      tau.push_back( ( (TMCParticle*) (*Particles)[ipart] )->GetLifetime() );
      theta.push_back( acos( pz[ipart]/sqrt(pow(px[ipart],2)+pow(py[ipart],2)+pow(pz[ipart],2))) );
      phi.push_back( TMath::ATan2( py[ipart], px[ipart] ) );
    }

    Tout->Fill();
  }

  TGraph *graph_sigma = new TGraph( event_num.size(), &(event_num[0]), &(sigma[0] ) );

  graph_sigma->SetMarkerStyle(20);
  graph_sigma->SetTitle("Total cross section (mb) vs events generated");
  graph_sigma->Write("graph_sigma");
  
  Fout->cd();
  Tout->Write();
  
}
