#include "gep_tree_with_spin.C"
#include "TChain.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TVector3.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TRotation.h"
#include "TString.h"
#include "TMath.h"

const double Mp = 0.938272046; //GeV/c^2
const double mu_p = 2.792847356; //mu_N
const double kappa_p = mu_p - 1.0; //anomalous magnetic moment
const double PI = TMath::Pi();
const double SBS_thetabend = 5.0*PI/180.0; //central bend angle

void G4SBS_spin_transport( const char *inputfilename, const char *outputfilename ){
  ifstream infile(inputfilename);
  
  TFile *fout = new TFile(outputfilename,"RECREATE");
  TTree *Tout = new TTree("Tout","SBS spin transport");

  double beta, gamma; //relativistic factors for the proton.
  double xfp, yfp, xpfp, ypfp, t; //include time-of-flight!
  double xtar, ytar, xptar, yptar, p;
  double chi, chiphi; //usual dipole precession angles
  //Euler angles for the trajectory bend:
  double phitrack, thetatrack, psitrack;
  double Pxtg, Pytg, Pztg; //In TRANSPORT coordinates!
  double Pxfp, Pyfp, Pzfp; //In TRANSPORT coordinates!

  Tout->Branch("beta",&beta);
  Tout->Branch("gamma",&gamma);
  Tout->Branch("xfp",&xfp);
  Tout->Branch("yfp",&yfp);
  Tout->Branch("xpfp",&xpfp);
  Tout->Branch("ypfp",&ypfp);
  Tout->Branch("t",&t);
  Tout->Branch("xtar",&xtar);
  Tout->Branch("ytar",&ytar);
  Tout->Branch("xptar",&xptar);
  Tout->Branch("yptar",&yptar);
  Tout->Branch("p",&p);
  Tout->Branch("chi",&chi);
  Tout->Branch("chiphi",&chiphi);
  Tout->Branch("phitrack",&phitrack);
  Tout->Branch("thetatrack",&thetatrack);
  Tout->Branch("psitrack",&psitrack);
  Tout->Branch("Pxtg",&Pxtg);
  Tout->Branch("Pytg",&Pytg);
  Tout->Branch("Pztg",&Pztg);
  Tout->Branch("Pxfp",&Pxfp);
  Tout->Branch("Pyfp",&Pyfp);
  Tout->Branch("Pzfp",&Pzfp);

  TChain *C = new TChain("T");
  
  TString currentline;
  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endlist") ){
    C->Add(currentline.Data());
  }

  gep_tree_with_spin *T = new gep_tree_with_spin(C);

  long nevent=0;
  
  while( T->GetEntry(nevent++) ){
    if( nevent%10000 == 0 ) cout << nevent << endl;

    if( T->Harm_FT_Track_ntracks == 1 &&
	(*(T->Harm_FT_Track_MID))[0] == 0 &&
	(*(T->Harm_FT_Track_P))[0]/T->ev_ep >= 0.995 ){ //good track in SBS:
      xfp = (*(T->Harm_FT_Track_X))[0];
      yfp = (*(T->Harm_FT_Track_Y))[0];
      xpfp = (*(T->Harm_FT_Track_Xp))[0];
      ypfp = (*(T->Harm_FT_Track_Yp))[0];
      t = (*(T->Harm_FT_Track_T))[0];
      p = T->ev_ep; //use initial momentum rather than "track" momentum
      beta = p/sqrt(pow(p,2)+pow(Mp,2));
      gamma = sqrt(1.0 + pow(p/Mp,2));

      //Now need target quantities!
      double px = T->ev_epx;
      double py = T->ev_epy;
      double pz = T->ev_epz;

      double SBStheta = T->gen_thsbs; //on beam right!

      TVector3 phall(px,py,pz);
      TVector3 phall_unit = phall.Unit();
      TVector3 SBSzaxis( -sin(SBStheta),0,cos(SBStheta) );
      TVector3 SBSxaxis(0,-1,0);
      TVector3 SBSyaxis = (SBSzaxis.Cross(SBSxaxis)).Unit();

      TVector3 punit_TRANSPORT( phall_unit.Dot( SBSxaxis ),
				phall_unit.Dot( SBSyaxis ),
				phall_unit.Dot( SBSzaxis ) );

      TVector3 vertex_hall( T->ev_vx, T->ev_vy, T->ev_vz ); //in m
      
      xptar = punit_TRANSPORT.X()/punit_TRANSPORT.Z();
      yptar = punit_TRANSPORT.Y()/punit_TRANSPORT.Z();
      
      TVector3 vertex_TRANSPORT( vertex_hall.Dot( SBSxaxis ),
				 vertex_hall.Dot( SBSyaxis ),
				 vertex_hall.Dot( SBSzaxis ) );

      xtar = vertex_TRANSPORT.X() - xptar * vertex_TRANSPORT.Z();
      ytar = vertex_TRANSPORT.Y() - yptar * vertex_TRANSPORT.Z();

      double thetabend = SBS_thetabend + atan(xptar) - atan(xpfp);
      double phibend = atan(ypfp) - atan(yptar);

      chi = gamma * kappa_p * thetabend;
      chiphi = gamma * kappa_p * phibend;

      //Compute Euler angles: for this we need to compute the axes of the tgt and fp comoving coordinate systems:
      TVector3 pfp_TRANSPORT( xpfp, ypfp, 1.0 );
      pfp_TRANSPORT = pfp_TRANSPORT.Unit();

      TVector3 SBSzaxis_fp( -sin(SBS_thetabend), 0, cos(SBS_thetabend) );
      TVector3 SBSxaxis_fp( cos(SBS_thetabend), 0, sin(SBS_thetabend) );
      TVector3 SBSyaxis_fp(0,1,0);

      TVector3 pfp_bent =
	pfp_TRANSPORT.X() * SBSxaxis_fp +
	pfp_TRANSPORT.Y() * SBSyaxis_fp +
	pfp_TRANSPORT.Z() * SBSzaxis_fp;

      TVector3 xaxis_comoving_fp = SBSxaxis_fp;
      TVector3 zaxis_comoving_fp = pfp_bent;
      TVector3 yaxis_comoving_fp = (zaxis_comoving_fp.Cross( xaxis_comoving_fp )).Unit();
      xaxis_comoving_fp = (yaxis_comoving_fp.Cross(zaxis_comoving_fp)).Unit();

      TVector3 xaxis_comoving_tgt(1,0,0);
      TVector3 zaxis_comoving_tgt = punit_TRANSPORT;
      TVector3 yaxis_comoving_tgt = (zaxis_comoving_tgt.Cross(xaxis_comoving_tgt)).Unit();
      xaxis_comoving_tgt = (yaxis_comoving_tgt.Cross(zaxis_comoving_tgt)).Unit();

      phitrack = atan2( yaxis_comoving_fp.Dot( xaxis_comoving_tgt ), -zaxis_comoving_fp.Dot( xaxis_comoving_tgt ) );
      psitrack = atan2( xaxis_comoving_fp.Dot( yaxis_comoving_tgt ), xaxis_comoving_fp.Dot( zaxis_comoving_tgt ) );
      thetatrack = acos( xaxis_comoving_fp.Dot( xaxis_comoving_tgt ) );

      Pxtg = T->ev_Sx;
      Pytg = T->ev_Sy;
      Pztg = T->ev_Sz;

      Pxfp = (*(T->Harm_FT_Track_Sx))[0];
      Pyfp = (*(T->Harm_FT_Track_Sy))[0];
      Pzfp = (*(T->Harm_FT_Track_Sz))[0];

      Tout->Fill();
    }
  }
   
  fout->Write();
  
}
