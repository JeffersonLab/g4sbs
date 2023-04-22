#include "TMatrixD.h"
#include "TVectorD.h"
#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include "TDecompSVD.h"
#include "TCut.h"
#include "TEventList.h"
//#include "g4sbs_tree.C"
#include "gep_tree_with_spin.C"
#include "TH1D.h"
#include "TH2D.h"
#include "TProfile.h"
#include "TVector3.h"
#include "TMath.h"
#include "TF2.h"

const double PI = TMath::Pi();

void G4SBS_optics_fit( const char *inputfilename, const char *outputfilename, int NMAX=1000000){
  double L0[2] = {0.8, 1.2192}; //meters, length of effective field boundary for first-order analytic model for momentum reconstruction:
  double alpha0[2] = {20.0*PI/180.0, 0.0}; //tilt of back of effective field boundary for first-order analytic model:
  
  ifstream inputfile(inputfilename);
  TFile *fout = new TFile( outputfilename, "RECREATE" );
  TChain *C = new TChain("T");

  TString opticsfilename = outputfilename;
  opticsfilename.ReplaceAll(".root",".txt");

  ofstream opticsfile(opticsfilename);
  
  TString currentline;
  while( currentline.ReadLine(inputfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      C->Add(currentline.Data());
    }
  }

  TCut global_cut = "";
  while( currentline.ReadLine(inputfile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      global_cut += currentline.Data();
    }
  }
  
  TEventList *elist = new TEventList("elist");

  C->Draw(">>elist",global_cut);

  gep_tree_with_spin *T = new gep_tree_with_spin(C);

  //Expansion is of the form: 
  // ( xptar yptar 1/p ytar ) = sum_{i+j+k+l+m<=N} C_{ijklm} xfp^i yfp^j xpfp^k ypfp^l xtar^m 
  // Fit chi^2 is sum_n=1^{Nevent} (xptar_n - sum_ijklm C_{ijklm} ...)^2, and similarly for yptar, ytar, 1/p, etc. 

  int order=4; 
  inputfile >> order; 
  int nbins = 200;
  inputfile >> nbins;
  int arm;  // 0 = BB, 1 = SBS
  inputfile >> arm;
  double chi2cut; 
  inputfile >> chi2cut; 
  double tracker_pitch_angle = 0.0;
  inputfile >> tracker_pitch_angle; //Positive value means that for up-bending particles, thetabend = tracker_pitch + thetatgt - thetafp  
  tracker_pitch_angle *= PI/180.0; //value is assumed to be given in degrees
  
  double pcentral[2] = {2.5,5.0};

  double fpmin[4], fpmax[4];    //ranges for plotting of fp variables: order is: x, y, xp, yp
  double tgtmin[5], tgtmax[5];  //ranges for plotting of tgt variables: order is: xp, yp, y, p, x

  for(int i=0; i<4; i++){
    inputfile >> fpmin[i] >> fpmax[i];
  }
  for(int i=0; i<5; i++){
    inputfile >> tgtmin[i] >> tgtmax[i];
  }

  //Histograms of "true" versions of fp track parameters:
  TH1D *hxfp = new TH1D("hxfp", "", nbins, fpmin[0], fpmax[0] );
  TH1D *hyfp = new TH1D("hyfp", "", nbins, fpmin[1], fpmax[1] );
  TH1D *hxpfp = new TH1D("hxpfp", "", nbins, fpmin[2], fpmax[2] );
  TH1D *hypfp = new TH1D("hypfp", "", nbins, fpmin[3], fpmax[3] );
  
  //Histograms of "recon" versions of fp track parameters:
  TH1D *hxfprecon = new TH1D("hxfprecon", "", nbins, fpmin[0], fpmax[0] );
  TH1D *hyfprecon = new TH1D("hyfprecon", "", nbins, fpmin[1], fpmax[1] );
  TH1D *hxpfprecon = new TH1D("hxpfprecon", "", nbins, fpmin[2], fpmax[2] );
  TH1D *hypfprecon = new TH1D("hypfprecon", "", nbins, fpmin[3], fpmax[3] );

  //Histograms of "recon" versions of fp track parameters:
  TH1D *hxfpdiff = new TH1D("hxfpdiff", "", nbins, -0.5e-3, 0.5e-3);
  TH1D *hyfpdiff = new TH1D("hyfpdiff", "", nbins, -0.5e-3, 0.5e-3);
  TH1D *hxpfpdiff = new TH1D("hxpfpdiff", "", nbins, -0.005, 0.005 );
  TH1D *hypfpdiff = new TH1D("hypfpdiff", "", nbins, -0.005, 0.005 );

  //Histograms of "true" versions of tgt track parameters:
  TH1D *hvz   = new TH1D("hvz","", nbins, -0.4, 0.4 );
  TH1D *hxtar = new TH1D("hxtar", "", nbins, tgtmin[4], tgtmax[4] );
  TH1D *hytar = new TH1D("hytar", "", nbins, tgtmin[2], tgtmax[2] );
  TH1D *hxptar = new TH1D("hxptar", "", nbins, tgtmin[0], tgtmax[0] );
  TH1D *hyptar = new TH1D("hyptar", "", nbins, tgtmin[1], tgtmax[1] );
  TH1D *hp     = new TH1D("hp", "", nbins, tgtmin[3], tgtmax[3] );

  //Histograms of "recon" versions of tgt track parameters:
  //TH1D *hxtarrecon = new TH1D("hxtarrecon", "", nbins, tgtmin[4], tgtmax[4] );
  TH1D *hytarrecon = new TH1D("hytarrecon", "", nbins, tgtmin[2], tgtmax[2] );
  TH1D *hxptarrecon = new TH1D("hxptarrecon", "", nbins, tgtmin[0], tgtmax[0] );
  TH1D *hyptarrecon = new TH1D("hyptarrecon", "", nbins, tgtmin[1], tgtmax[1] );
  TH1D *hprecon     = new TH1D("hprecon", "", nbins, tgtmin[3], tgtmax[3] );
  
  // TH1D *hxtardiff = new TH1D("hxtardiff", "", nbins, tgtmin[4], tgtmax[4] );
  TH1D *hytardiff = new TH1D("hytardiff", "", nbins, -0.025, 0.025 );
  TH1D *hxptardiff = new TH1D("hxptardiff", "", nbins, -0.02, 0.02 );
  TH1D *hyptardiff = new TH1D("hyptardiff", "", nbins, -0.02, 0.02 );
  TH1D *hpdiff     = new TH1D("hpdiff", "", nbins, -0.2,0.2 );

  TH2D *hxptardiff_xtar = new TH2D("hxptardiff_xtar", "", nbins, tgtmin[4], tgtmax[4], nbins, -0.02, 0.02 );
  TH2D *hxptardiff_xptar = new TH2D("hxptardiff_xptar", "", nbins, tgtmin[0], tgtmax[0], nbins, -0.02, 0.02 );
  TH2D *hxptardiff_yptar = new TH2D("hxptardiff_yptar", "", nbins, tgtmin[1], tgtmax[1], nbins, -0.02, 0.02 );
  TH2D *hxptardiff_ytar = new TH2D("hxptardiff_ytar", "", nbins, tgtmin[2], tgtmax[2], nbins, -0.02, 0.02 );
  TH2D *hxptardiff_p = new TH2D("hxptardiff_p", "", nbins, tgtmin[3], tgtmax[3], nbins, -0.02, 0.02 );

  TH2D *hyptardiff_xtar = new TH2D("hyptardiff_xtar", "", nbins, tgtmin[4], tgtmax[4], nbins, -0.02, 0.02 );
  TH2D *hyptardiff_xptar = new TH2D("hyptardiff_xptar", "", nbins, tgtmin[0], tgtmax[0], nbins, -0.02, 0.02 );
  TH2D *hyptardiff_yptar = new TH2D("hyptardiff_yptar", "", nbins, tgtmin[1], tgtmax[1], nbins, -0.02, 0.02 );
  TH2D *hyptardiff_ytar = new TH2D("hyptardiff_ytar", "", nbins, tgtmin[2], tgtmax[2], nbins, -0.02, 0.02 );
  TH2D *hyptardiff_p = new TH2D("hyptardiff_p", "", nbins, tgtmin[3], tgtmax[3], nbins, -0.02, 0.02 );

  TH2D *hytardiff_xtar = new TH2D("hytardiff_xtar", "", nbins, tgtmin[4], tgtmax[4], nbins, -0.025, 0.025 );
  TH2D *hytardiff_xptar = new TH2D("hytardiff_xptar", "", nbins, tgtmin[0], tgtmax[0], nbins, -0.025, 0.025 );
  TH2D *hytardiff_yptar = new TH2D("hytardiff_yptar", "", nbins, tgtmin[1], tgtmax[1], nbins, -0.025, 0.025 );
  TH2D *hytardiff_ytar = new TH2D("hytardiff_ytar", "", nbins, tgtmin[2], tgtmax[2], nbins, -0.025, 0.025 );
  TH2D *hytardiff_p = new TH2D("hytardiff_p", "", nbins, tgtmin[3], tgtmax[3], nbins, -0.025, 0.025 );
  
  TH2D *hpdiff_xtar = new TH2D("hpdiff_xtar", "", nbins, tgtmin[4], tgtmax[4], nbins, -0.2, 0.2 );
  TH2D *hpdiff_xptar = new TH2D("hpdiff_xptar", "", nbins, tgtmin[0], tgtmax[0], nbins, -0.2, 0.2 );
  TH2D *hpdiff_yptar = new TH2D("hpdiff_yptar", "", nbins, tgtmin[1], tgtmax[1], nbins, -0.2, 0.2 );
  TH2D *hpdiff_ytar = new TH2D("hpdiff_ytar", "", nbins, tgtmin[2], tgtmax[2], nbins, -0.2, 0.2 );
  TH2D *hpdiff_p = new TH2D("hpdiff_p", "", nbins, tgtmin[3], tgtmax[3], nbins, -0.2, 0.2 );
  
  TH2D *hpdiff_xfp = new TH2D("hpdiff_xfp", "", nbins, fpmin[0], fpmax[0], nbins, -0.2, 0.2 );
  TH2D *hpdiff_yfp = new TH2D("hpdiff_yfp", "", nbins, fpmin[1], fpmax[1], nbins, -0.2, 0.2 );
  TH2D *hpdiff_xpfp = new TH2D("hpdiff_xpfp", "", nbins, fpmin[2], fpmax[2], nbins, -0.2, 0.2 );
  TH2D *hpdiff_ypfp = new TH2D("hpdiff_ypfp", "", nbins, fpmin[3], fpmax[3], nbins, -0.2, 0.2 );

  TH1D *hpinvdiff = new TH1D("hpinvdiff", "", nbins, -0.05, 0.05 );

  TH1D *hvzdiff = new TH1D("hvzdiff", "", nbins, -0.3, 0.3 );

  TH2D *hvzdiff_xtar = new TH2D("hvzdiff_xtar", "", nbins, tgtmin[4],tgtmax[4], nbins, -0.3,0.3);
  TH2D *hvzdiff_ytar = new TH2D("hvzdiff_ytar", "", nbins, tgtmin[2],tgtmax[2], nbins, -0.3,0.3);
  TH2D *hvzdiff_xptar = new TH2D("hvzdiff_xptar", "", nbins, tgtmin[0],tgtmax[0], nbins, -0.3,0.3);
  TH2D *hvzdiff_yptar = new TH2D("hvzdiff_yptar", "", nbins, tgtmin[1],tgtmax[1], nbins, -0.3,0.3);
  TH2D *hvzdiff_p = new TH2D("hvzdiff_p", "", nbins, tgtmin[3],tgtmax[3], nbins, -0.3,0.3);
  
  int nparams = 0;

  vector<int> xtar_expon;
  vector<int> xfp_expon;
  vector<int> xpfp_expon;

  for(int i=0; i<=order; i++){
    for(int j=0; j<=order-i; j++){
      for(int k=0; k<=order-i-j; k++){
	for(int l=0; l<=order-i-j-k; l++){
	  for(int m=0; m<=order-i-j-k-l; m++){
	    nparams++;
	    xtar_expon.push_back( i );
	    xfp_expon.push_back( m );
	    xpfp_expon.push_back( k );
	  }
	}
      }
    }
  }

  //TMatrixD M_xptar(nparams,nparams), M_yptar(nparams,nparams), M_ytar(nparams,nparams), M_pinv(nparams,nparams);
  TMatrixD M(nparams,nparams);
  TVectorD b_xptar(nparams), b_yptar(nparams), b_ytar(nparams), b_pinv(nparams);
  
  for(int i=0; i<nparams; i++){
    for(int j=0; j<nparams; j++){
      M(i,j) = 0.0;
    }
    b_xptar(i) = 0.0;
    b_yptar(i) = 0.0;
    b_ytar(i) = 0.0;
    b_pinv(i) = 0.0;
  }

  long nevent = 0;
  while( T->GetEntry(elist->GetEntry(nevent++)) && nevent < NMAX ){
    if( nevent%1000 == 0 ){
      cout << nevent << endl;
    }
    
    double theta0;
    if( arm == 0 ){ //BB:
      theta0 = T->gen_thbb; //on beam left:
    } else if( arm == 1 ){ //SBS:
      theta0 = T->gen_thsbs; //on beam right
    }

    double vx, vy, vz, px, py, pz;
    double p, xptar, yptar, ytar, xtar;
    double xfp, yfp, xpfp, ypfp;

    vx = T->ev_vx;
    vy = T->ev_vy;
    vz = T->ev_vz;

    TVector3 vertex(vx,vy,vz);

    
    
    if( arm == 0 ){ //BB
      p = T->ev_ep;

      //px = T->ev_epx;
      //py = T->ev_epy;
      //pz = T->ev_epz;
      
      TVector3 pvect( p*sin(T->ev_th)*cos(T->ev_ph), p*sin(T->ev_th)*sin(T->ev_ph), p*cos(T->ev_th) );
      TVector3 BB_zaxis( sin(theta0), 0.0, cos(theta0) ); //BB is on beam left, global x axis points to beam left
      TVector3 BB_xaxis(0,-1,0); //X axis of transport coordinates is vertically down:
      //TVector3 BB_xaxis = (BB_yaxis.Cross(BB_zaxis)).Unit();
      TVector3 BB_yaxis = (BB_zaxis.Cross(BB_xaxis)).Unit();
      
      //TVector3 pvect_BB = pvect.Dot(BB_xaxis)*BB_xaxis + pvect.Dot(BB_yaxis)*BB_yaxis + pvect.Dot(BB_zaxis)*BB_zaxis;
      TVector3 pvect_BB( pvect.Dot(BB_xaxis), pvect.Dot(BB_yaxis), pvect.Dot(BB_zaxis) );

      xptar = pvect_BB.X()/pvect_BB.Z();
      yptar = pvect_BB.Y()/pvect_BB.Z();
      
      TVector3 vertex_BB(vertex.Dot(BB_xaxis), vertex.Dot(BB_yaxis), vertex.Dot(BB_zaxis) );

      ytar = vertex_BB.Y() - yptar * vertex_BB.Z();
      xtar = vertex_BB.X() - xptar * vertex_BB.Z();
      
    } else {
      p = T->ev_np;
      px = T->ev_npx;
      py = T->ev_npy;
      pz = T->ev_npz;
      
      TVector3 pvect(p*sin(T->ev_nth)*cos(T->ev_nph), p*sin(T->ev_nth)*sin(T->ev_nph), p*cos(T->ev_nth) );
      TVector3 SBS_zaxis( -sin(theta0), 0, cos(theta0) );
      TVector3 SBS_xaxis(0,-1,0);
      TVector3 SBS_yaxis = (SBS_zaxis.Cross(SBS_xaxis)).Unit();

      TVector3 pvect_SBS( pvect.Dot(SBS_xaxis), pvect.Dot(SBS_yaxis), pvect.Dot(SBS_zaxis) );
      
      xptar = pvect_SBS.X()/pvect_SBS.Z();
      yptar = pvect_SBS.Y()/pvect_SBS.Z();

      TVector3 vertex_SBS( vertex.Dot(SBS_xaxis), vertex.Dot(SBS_yaxis), vertex.Dot(SBS_zaxis) );
      ytar = vertex_SBS.Y() - yptar * vertex_SBS.Z();
      xtar = vertex_SBS.X() - xptar * vertex_SBS.Z();
      
    }    

   

    //Fill focal plane track variables:
    for(int itrack=0; itrack<T->ntracks; itrack++){
      if( (*(T->trackerid))[itrack] == arm && (*(T->trackchi2true))[itrack]/(*(T->trackndf))[itrack] < chi2cut ){
	//
	p = (*(T->trackp))[itrack];
	
	xfp = (*(T->trackx))[itrack];
	yfp = (*(T->tracky))[itrack];
	xpfp = (*(T->trackxp))[itrack];
	ypfp = (*(T->trackyp))[itrack];
	
	hxfp->Fill(xfp);
	hyfp->Fill(yfp);
	hxpfp->Fill(xpfp);
	hypfp->Fill(ypfp);

	hxfprecon->Fill( (*(T->trackxfit))[itrack] );
	hyfprecon->Fill( (*(T->trackyfit))[itrack] );
	hxpfprecon->Fill( (*(T->trackxpfit))[itrack] );
	hypfprecon->Fill( (*(T->trackypfit))[itrack] );

	goodtrack = true;
	break;
      }
    }

    if( goodtrack ){ //Increment fit matrices:
      
      hxtar->Fill(xtar);
      hytar->Fill(ytar);
      hxptar->Fill(xptar);
      hyptar->Fill(yptar);
      hp->Fill(p);
      hvz->Fill(vz);
      
      //Are the equations below correct? 
      // chi^2 = sum_{i=1}^{N} (xtgt_i - sum_jklmn C_jklmn xfp^j yfp^k xpfp^l ypfp^m xtar^n)^2
      // Momentum reconstruction is a special case, we want to fit a 
      // functional form of 1/p = A * thetabend + sum C_jklmn ...
      // in this case, chi^2 = sum_i (1/p - (A*thetabend + sum_jklmn C_jklmn * term))^2
      
      //In the case of the momentum, what we want to fit is the "small" corrections as a function of kinematics, to a first order analytical model in 
      //which we treat the dipole field in the dispersive plane as a constant field with an effective boundary. In both cases, the front and back edges of the field boundary are planes, (lines in one dimension). In SBS case, both front and back edges of field boundary are vertical. In BB case, back edge of field bdry is 
      // tilted at an angle of 20 degrees:
      double thetabend = tracker_pitch_angle + atan(xptar) - atan(xpfp);

      double p_firstorder;
      if( arm == 0 ){
	p_firstorder = (0.286378 + 0.141452*(xfp-0.8*xpfp))/thetabend;
      } else {
	p_firstorder = -0.5161/thetabend;
      }

      vector<double> term(nparams);
      //      vector<double> term_y(nparams);
      int ipar=0;
      for(int i=0; i<=order; i++){
	for(int j=0; j<=order-i; j++){
	  for(int k=0; k<=order-i-j; k++){
	    for(int l=0; l<=order-i-j-k; l++){
	      for(int m=0; m<=order-i-j-k-l; m++){
		term[ipar] = pow(xfp,m)*pow(yfp,l)*pow(xpfp,k)*pow(ypfp,j)*pow(xtar,i);
		
		b_xptar(ipar) += term[ipar]*xptar;
		
		if( i == 0 ){
		  b_yptar(ipar) += term[ipar]*yptar;
		  b_ytar(ipar) += term[ipar]*ytar;
		}
		  //b_ytar(ipar) += term[ipar]*T->ev_vz;
		  //} 
		//b_pinv(ipar) += term[ipar]*(p/pcentral[arm]-1.0);
		b_pinv(ipar) += term[ipar]/p;
		ipar++;
	      }
	    }
	  }
	}
      }
      
      for(ipar=0; ipar<nparams; ipar++){
	for(int jpar=0; jpar<nparams; jpar++){
	  //double term_ij = term[ipar]*term[jpar];
	  M(ipar,jpar) += term[ipar]*term[jpar];
	}
      }
      
    }
  }

  cout << "order, nparams = " << nparams << endl;

  cout << "Setting up SVD for xptar:" << endl;
  TDecompSVD A_xptar(M);
  cout << "Setting up SVD for pinv:" << endl;
  TDecompSVD A_pinv(M);
  
  for(int ipar=0; ipar<nparams; ipar++){
    if( xtar_expon[ipar] > 0 ){
      M(ipar,ipar) = 1.0;
      for(int jpar=0; jpar<nparams; jpar++){
  	if( jpar != ipar ) M(ipar,jpar) = 0.0;
      }
    }
  }
  
  cout << "Setting up SVD for yptar:" << endl;
  TDecompSVD A_yptar(M);
  cout << "Setting up SVD for ytar:" << endl;
  TDecompSVD A_ytar(M);
  

  cout << "Solving for xptar coefficients:" << endl;
  bool good_xptar = A_xptar.Solve(b_xptar);
  cout << "xptar done, success = " << good_xptar << endl;
  cout << "Solving for yptar coefficients:" << endl;
  bool good_yptar = A_yptar.Solve(b_yptar);
  cout << "yptar done, success = " << good_yptar << endl;
  cout << "Solving for ytar coefficients:" << endl;
  bool good_ytar = A_ytar.Solve(b_ytar);
  cout << "ytar done, success = " << good_ytar << endl;
  cout << "Solving for 1/p coefficients:" << endl;
  bool good_pinv = A_pinv.Solve(b_pinv);
  cout << "1/p done, success = " << good_pinv << endl;

  //Write out the optics file:
  // Order is: xptar, yptar, ytar, 1/p, exponents:
  opticsfile << nparams << endl;
  int ipar = 0;
  for(int i=0; i<=order; i++){
    for(int j=0; j<=order-i; j++){
      for(int k=0; k<=order-i-j; k++){
	for(int l=0; l<=order-i-j-k; l++){
	  for(int m=0; m<=order-i-j-k-l; m++){
	    char ccoeff[100];
	    sprintf( ccoeff, "%15.8g", b_xptar(ipar) );
	    opticsfile << ccoeff;
	    sprintf( ccoeff, "%15.8g", b_yptar(ipar) );
	    opticsfile << ccoeff;
	    sprintf( ccoeff, "%15.8g", b_ytar(ipar) );
	    opticsfile << ccoeff;
	    sprintf( ccoeff, "%15.8g", b_pinv(ipar) );
	    opticsfile << ccoeff;

	    char cexpon[100];
	    sprintf( cexpon, "  %d %d %d %d %d", m, l, k, j, i );
	    opticsfile << cexpon << endl;
	    
	    ipar++;
	  }
	}
      }
    }
  }	   
 

  cout << "arm = " << arm << endl;

  nevent = 0;
  while( T->GetEntry( elist->GetEntry(nevent++))){
    if( nevent%1000 == 0 ){
      cout << nevent << endl;
    }

    double theta0;
    if( arm == 0 ){ //BB:
      theta0 = T->gen_thbb; //on beam right
    } else if( arm == 1 ){ //SBS:
      theta0 = T->gen_thhcal; //on beam left
    }
    
    double vx, vy, vz, px, py, pz;
    double p, xptar, yptar, ytar, xtar;
    double p_fit, xptar_fit, yptar_fit, ytar_fit; //Fit is reconstructed using fit coefficients, no smearing for detector resolution
    double p_recon, xptar_recon, yptar_recon, ytar_recon; //recon is reconstructed using fit coefficients, fp quantities smeared by det. resolution
    //double pthetabend_fit, pthetabend_recon;
    double pinv_fit, pinv_recon;
    double xfp, yfp, xpfp, ypfp;
    double xfp_fit, yfp_fit, xpfp_fit, ypfp_fit;

    pinv_fit = 0.0;
    pinv_recon = 0.0;

    vx = T->ev_vx;
    vy = T->ev_vy;
    vz = T->ev_vz;

    TVector3 vertex(vx,vy,vz);

    if( arm == 0 ){
      p = T->ev_ep;
      //px = T->ev_epx;
      //py = T->ev_epy;
      //pz = T->ev_epz;
      
      TVector3 pvect( p*sin(T->ev_th)*cos(T->ev_ph), p*sin(T->ev_th)*sin(T->ev_ph), p*cos(T->ev_th) );
      TVector3 BB_zaxis( -sin(theta0), 0.0, cos(theta0) ); //BB is on beam right, global x axis points to beam left
      TVector3 BB_xaxis(0,-1,0); //X axis of transport coordinates is vertically down:
      //TVector3 BB_xaxis = (BB_yaxis.Cross(BB_zaxis)).Unit();
      TVector3 BB_yaxis = (BB_zaxis.Cross(BB_xaxis)).Unit();
      
      //TVector3 pvect_BB = pvect.Dot(BB_xaxis)*BB_xaxis + pvect.Dot(BB_yaxis)*BB_yaxis + pvect.Dot(BB_zaxis)*BB_zaxis;
      TVector3 pvect_BB( pvect.Dot(BB_xaxis), pvect.Dot(BB_yaxis), pvect.Dot(BB_zaxis) );

      xptar = pvect_BB.X()/pvect_BB.Z();
      yptar = pvect_BB.Y()/pvect_BB.Z();
      
      TVector3 vertex_BB(vertex.Dot(BB_xaxis), vertex.Dot(BB_yaxis), vertex.Dot(BB_zaxis) );

      ytar = vertex_BB.Y() - yptar * vertex_BB.Z();
      xtar = vertex_BB.X() - xptar * vertex_BB.Z();
      
    } else {
      p = T->ev_np;
      px = T->ev_npx;
      py = T->ev_npy;
      pz = T->ev_npz;
      
      //TVector3 pvect(px,py,pz);
      TVector3 pvect( p*sin(T->ev_nth)*cos(T->ev_nph), p*sin(T->ev_nth)*sin(T->ev_nph), p*cos(T->ev_nth) );
      TVector3 SBS_zaxis( sin(theta0), 0, cos(theta0) );
      TVector3 SBS_xaxis(0,-1,0);
      TVector3 SBS_yaxis = (SBS_zaxis.Cross(SBS_xaxis)).Unit();

      TVector3 pvect_SBS( pvect.Dot(SBS_xaxis), pvect.Dot(SBS_yaxis), pvect.Dot(SBS_zaxis) );
      
      xptar = pvect_SBS.X()/pvect_SBS.Z();
      yptar = pvect_SBS.Y()/pvect_SBS.Z();

      TVector3 vertex_SBS( vertex.Dot(SBS_xaxis), vertex.Dot(SBS_yaxis), vertex.Dot(SBS_zaxis) );
      ytar = vertex_SBS.Y() - yptar * vertex_SBS.Z();
      xtar = vertex_SBS.X() - xptar * vertex_SBS.Z();
      
    }
    
    bool goodtrack = false;
    
    //Fill focal plane track variables:
    for(int itrack=0; itrack<T->ntracks; itrack++){
      if( (*(T->trackerid))[itrack] == arm && (*(T->trackchi2true))[itrack]/(*(T->trackndf))[itrack] < chi2cut ){
	p = (*(T->trackp))[itrack];
	//
	xfp = (*(T->trackx))[itrack];
	yfp = (*(T->tracky))[itrack];
	xpfp = (*(T->trackxp))[itrack];
	ypfp = (*(T->trackyp))[itrack];
	
	xfp_fit = (*(T->trackxfit))[itrack];
	xpfp_fit = (*(T->trackxpfit))[itrack];
	yfp_fit = (*(T->trackyfit))[itrack];
	ypfp_fit = (*(T->trackypfit))[itrack];

	hxfpdiff->Fill(xfp_fit-xfp);
	hyfpdiff->Fill(yfp_fit-yfp);
	hxpfpdiff->Fill(xpfp_fit-xpfp);
	hypfpdiff->Fill(ypfp_fit-ypfp);

	goodtrack = true;
	break;
      }
    }
    
    if( goodtrack ){ //Increment fit matrices:
      int ipar = 0;
      xptar_fit = 0.0;
      yptar_fit = 0.0;
      ytar_fit = 0.0;
      //pthetabend_fit = 0.0;

      xptar_recon = 0.0;
      yptar_recon = 0.0;
      ytar_recon = 0.0;
      //      pthetabend_recon = 0.0;
      
      ipar=0;
      for(int i=0; i<=order; i++){
	for(int j=0; j<=order-i; j++){
	  for(int k=0; k<=order-i-j; k++){
	    for(int l=0; l<=order-i-j-k; l++){
	      for(int m=0; m<=order-i-j-k-l; m++){
		double term = pow(xfp,m)*pow(yfp,l)*pow(xpfp,k)*pow(ypfp,j)*pow(xtar,i);
		xptar_fit += b_xptar(ipar)*term;
		yptar_fit += b_yptar(ipar)*term;
		ytar_fit += b_ytar(ipar)*term;
		//pthetabend_fit += b_pinv(ipar)*term;
		pinv_fit += b_pinv(ipar)*term;

		term = pow(xfp_fit,m)*pow(yfp_fit,l)*pow(xpfp_fit,k)*pow(ypfp_fit,j)*pow(xtar,i);
		xptar_recon += b_xptar(ipar)*term;
		yptar_recon += b_yptar(ipar)*term;
		ytar_recon += b_ytar(ipar)*term;
		pinv_recon += b_pinv(ipar)*term;
		//pthetabend_recon += b_pinv(ipar)*term;
		ipar++;
	      }
	    }
	  }
	}
      }
            
      double thetabend_true = tracker_pitch_angle + atan(xptar) - atan(xpfp);
      double thetabend_fit = tracker_pitch_angle + atan(xptar_fit) - atan(xpfp);
      double thetabend_recon = tracker_pitch_angle + atan(xptar_recon) - atan(xpfp_fit);

      double pfirstorder_true, pfirstorder_recon;
      if(arm==0){
	pfirstorder_true = (0.286378 + 0.141452*(xfp-0.8*xpfp))/thetabend_true;
	pfirstorder_recon  = (0.286378 + 0.141452*(xfp_fit-0.8*xpfp_fit))/thetabend_recon;
      } else {
	pfirstorder_true = -0.5161/thetabend_true;
	pfirstorder_recon = -0.5161/thetabend_recon;
      }
      //p_fit = pthetabend_fit/thetabend_fit;
      //p_recon = pthetabend_recon/thetabend_recon;

      //p_fit = pthetabend_fit + pfirstorder_true;
      //p_recon = pthetabend_recon + pfirstorder_recon;
      p_fit = 1.0/pinv_fit;
      p_recon = 1.0/pinv_recon;
      //p_fit = pcentral[arm]*(1.0 + pinv_fit);
      //p_recon = pcentral[arm]*(1.0 + pinv_recon);

      // double vz_fit = ytar_fit / ( sin( theta;
      // double vz_recon = ytar_recon;
      double vz_fit, vz_recon;

      if( arm == 0 ){ //BB, beam right:
	vz_fit = ytar_fit / (sin(theta0) - cos(theta0)*yptar_fit);
	vz_recon = ytar_recon / (sin(theta0) - cos(theta0)*yptar_recon);
      } else { //SBS, beam left:
	vz_fit = -ytar_fit / (sin(theta0) + cos(theta0)*yptar_fit);
	vz_recon = -ytar_recon / (sin(theta0) + cos(theta0)*yptar_recon);
      }

      hvzdiff->Fill( vz_fit - T->ev_vz );
      hvzdiff_xtar->Fill(xtar, vz_fit - T->ev_vz );
      hvzdiff_ytar->Fill(ytar, vz_fit - T->ev_vz );
      hvzdiff_xptar->Fill(xptar, vz_fit - T->ev_vz );
      hvzdiff_yptar->Fill(yptar, vz_fit - T->ev_vz );
      hvzdiff_p->Fill(p, vz_fit - T->ev_vz );

      // if( arm == 0 ){ //bigbite, beam right:
      // 	ytar_fit = vz_fit * (sin(theta0) - cos(theta0)*yptar_fit);
      // 	ytar_recon = vz_recon * (sin(theta0) - cos(theta0)*yptar_recon);
      // } else { //SBS, beam left:
      // 	ytar_fit = -vz_fit * (sin(theta0) + cos(theta0)*yptar_fit);
      // 	ytar_recon = -vz_recon * (sin(theta0) + cos(theta0)*yptar_recon);
      // }

      hpinvdiff->Fill( pinv_fit - 1.0/p );

      hytarrecon->Fill( ytar_fit );
      hyptarrecon->Fill( yptar_fit );
      hxptarrecon->Fill( xptar_fit );
      hprecon->Fill(p_fit );

      hytardiff->Fill( ytar_fit - ytar );
      hxptardiff->Fill( xptar_fit - xptar );
      hyptardiff->Fill( yptar_fit - yptar );
      hpdiff->Fill( (p_fit/p - 1.0) );
      
      hxptardiff_xtar->Fill( xtar, xptar_fit - xptar );
      hxptardiff_xptar->Fill( xptar, xptar_fit - xptar );
      hxptardiff_yptar->Fill( yptar, xptar_fit - xptar );
      hxptardiff_ytar->Fill( ytar, xptar_fit - xptar );
      hxptardiff_p->Fill( p, xptar_fit - xptar );

      hyptardiff_xtar->Fill(xtar, yptar_fit - yptar );
      hyptardiff_xptar->Fill( xptar, yptar_fit - yptar );
      hyptardiff_yptar->Fill( yptar, yptar_fit - yptar );
      hyptardiff_ytar->Fill( ytar, yptar_fit - yptar );
      hyptardiff_p->Fill( p, yptar_fit - yptar );
      
      hytardiff_xtar->Fill( xtar, ytar_fit - ytar );
      hytardiff_xptar->Fill( xptar, ytar_fit - ytar );
      hytardiff_yptar->Fill( yptar, ytar_fit - ytar );
      hytardiff_ytar->Fill( ytar, ytar_fit - ytar );
      hytardiff_p->Fill( p, ytar_fit - ytar );
      
      hpdiff_xtar->Fill( xtar, p_fit/p - 1.0 );
      hpdiff_xptar->Fill( xptar, p_fit/p - 1.0 );
      hpdiff_yptar->Fill( yptar, p_fit/p - 1.0 );
      hpdiff_ytar->Fill( ytar, p_fit/p - 1.0 );
      hpdiff_p->Fill( p, p_fit/p - 1.0 );

      hpdiff_xfp->Fill( xfp, p_fit/p - 1.0 );
      hpdiff_yfp->Fill( yfp, p_fit/p - 1.0 );
      hpdiff_xpfp->Fill( xpfp, p_fit/p - 1.0 );
      hpdiff_ypfp->Fill( ypfp, p_fit/p - 1.0 );
      
    }
  }

  elist->Delete();
  fout->Write();
}
