#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "gep_optics_tree.C"
#include "TDecompSVD.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include <vector>
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRotation.h"
#include "TH1D.h"
#include "TH2D.h"

//double PI = TMath::Pi();
const double Mp = 0.938272046; //GeV/c^2
const double mu_p = 2.792847356; //mu_N
const double kappa_p = mu_p - 1.0; //anomalous magnetic moment
const double PI = TMath::Pi();
const double SBS_thetabend = 5.0*PI/180.0; //central bend angle

//const double SBS_thetabend = 5.0*PI/180.0;

void Fit_SBS_optics_and_spin( const char *rootfilename, int order_optics=4, int order_spin=4, const char *outfilename="SBS_optics_fitresult.root", int use_xtar_flag=0 ){
  int ncoeff_optics = 0;

  int ncoeff_optics_nonzero = 0;
  
  vector<int> xtar_expon_optics, xfp_expon_optics, xpfp_expon_optics;
  
  for( int i=0; i<=order_optics; i++ ){
    for( int j=0; j<=order_optics-i; j++ ){
      for( int k=0; k<=order_optics-i-j; k++ ){
	for( int l=0; l<=order_optics-i-j-k; l++ ){
	  for( int m=0; m<=order_optics-i-j-k-l; m++ ){
	    ncoeff_optics++;
	    xtar_expon_optics.push_back( i );
	    xfp_expon_optics.push_back( m );
	    xpfp_expon_optics.push_back( k );
	    if( i == 0 ) ncoeff_optics_nonzero++;
	  }
	}
      }
    }
  }
  
  cout << "nterms optics fit = " << ncoeff_optics << endl;
  
  
  int ncoeff_spin = 0;
  int ncoeff_spin_nonzero = 0;
  
  vector<int> xtar_expon_spin, xfp_expon_spin, xpfp_expon_spin;

  for( int i=0; i<=order_spin; i++ ){
    for( int j=0; j<=order_spin-i; j++ ){
      for( int k=0; k<=order_spin-i-j; k++ ){
	for( int l=0; l<=order_spin-i-j-k; l++ ){
	  for( int m=0; m<=order_spin-i-j-k-l; m++ ){
	    ncoeff_spin++;
	    xtar_expon_spin.push_back( i );
	    //    xfp_expon_spin.push_back( m );
	    //xpfp_expon_spin.push_back( k );
	    if( i == 0 ) ncoeff_spin_nonzero++;
	  }
	}
      }
    }
  }

  TMatrixD Mforward(ncoeff_optics,ncoeff_optics);
  TVectorD b_xfp(ncoeff_optics), b_yfp(ncoeff_optics), b_xpfp(ncoeff_optics), b_ypfp(ncoeff_optics);
  
  TMatrixD Moptics(ncoeff_optics,ncoeff_optics);
  TVectorD b_xptar(ncoeff_optics), b_yptar(ncoeff_optics), b_ytar(ncoeff_optics), b_pinv(ncoeff_optics);

  for( int ipar=0; ipar<ncoeff_optics; ipar++ ){
    b_xptar(ipar) = 0.0;
    b_yptar(ipar) = 0.0;
    b_ytar(ipar) = 0.0;
    b_pinv(ipar) = 0.0;

    b_xfp(ipar) = 0.0;
    b_yfp(ipar) = 0.0;
    b_xpfp(ipar) = 0.0;
    b_ypfp(ipar) = 0.0;
    for( int jpar=0; jpar<ncoeff_optics; jpar++ ){
      Moptics(ipar,jpar) = 0.0;
      Mforward(ipar,jpar) = 0.0;
    }
  }
  
  TMatrixD Mspinx(ncoeff_spin,ncoeff_spin);
  TMatrixD Mspiny(ncoeff_spin,ncoeff_spin);
  TMatrixD Mspinz(ncoeff_spin,ncoeff_spin);
  TVectorD b_Sxx(ncoeff_spin), b_Sxy(ncoeff_spin), b_Sxz(ncoeff_spin);
  TVectorD b_Syx(ncoeff_spin), b_Syy(ncoeff_spin), b_Syz(ncoeff_spin);
  TVectorD b_Szx(ncoeff_spin), b_Szy(ncoeff_spin), b_Szz(ncoeff_spin);

  TMatrixD b_Sjx(ncoeff_spin,3);
  TMatrixD b_Sjy(ncoeff_spin,3);
  TMatrixD b_Sjz(ncoeff_spin,3);
  
  for( int ipar=0; ipar<ncoeff_spin; ipar++ ){
    b_Sxx(ipar) = 0.0;
    b_Sxy(ipar) = 0.0;
    b_Sxz(ipar) = 0.0;

    b_Syx(ipar) = 0.0;
    b_Syy(ipar) = 0.0;
    b_Syz(ipar) = 0.0;

    b_Szx(ipar) = 0.0;
    b_Szy(ipar) = 0.0;
    b_Szz(ipar) = 0.0;
    
    for( int jpar=0; jpar<ncoeff_spin; jpar++ ){
      Mspinx(ipar,jpar) = 0.0;
      Mspiny(ipar,jpar) = 0.0;
      Mspinz(ipar,jpar) = 0.0;
    }
  }
  
  TChain *C = new TChain("Tout");
  C->Add(rootfilename);

  gep_optics_tree *T = new gep_optics_tree(C);
  
  long nevent=0;
  while( T->GetEntry(nevent++) ){
    if( nevent%100 == 0 ) cout << nevent << endl;

    //optics fitting:
    //

    //First compute bend angle:

    TVector3 zaxis_fp(-sin(SBS_thetabend),0,cos(SBS_thetabend));
    TVector3 yaxis_fp(0,1,0);
    TVector3 xaxis_fp = yaxis_fp.Cross(zaxis_fp).Unit();

    TVector3 nhat_fp(T->xpfp, T->ypfp, 1.0 );
    nhat_fp = nhat_fp.Unit();
    TVector3 nhat_fp_global = nhat_fp.X()*xaxis_fp +
      nhat_fp.Y()*yaxis_fp +
      nhat_fp.Z()*zaxis_fp;

    TVector3 nhat_tgt(T->xptar, T->yptar, 1.0);
    nhat_tgt = nhat_tgt.Unit();

    double thetabend = acos(nhat_fp_global.Dot(nhat_tgt));

    TVector3 bend_axis = nhat_tgt.Cross( nhat_fp_global ).Unit();
    
    //    cout << "(p (GeV), thetabend (deg) ) = (" << T->p << ", " << thetabend*180.0/PI << ")" << endl;
    
    vector<double> term(ncoeff_optics);
    vector<double> fterm(ncoeff_optics);
    
    int icoeff=0;
    for( int i=0; i<=order_optics; i++ ){
      for( int j=0; j<=order_optics-i; j++ ){
	for( int k=0; k<=order_optics-i-j; k++ ){
	  for( int l=0; l<=order_optics-i-j-k; l++ ){
	    for( int m=0; m<=order_optics-i-j-k-l; m++ ){
	      term[icoeff] = pow(T->xfp,m)*pow(T->yfp,l)*pow(T->xpfp,k)*pow(T->ypfp,j)*pow(T->xtar,i);
	      fterm[icoeff] = pow(T->xptar,m)*pow(T->yptar,l)*pow(T->ytar,k)*pow(1.0/T->p,j)*pow(T->xtar,i);
	      
	      b_xptar[icoeff] += term[icoeff]*T->xptar;
	      //if( i == 0 ){ //don't include xtar-dependent terms in yptar, ytar fits:
	      b_yptar[icoeff] += term[icoeff]*T->yptar;
	      b_ytar[icoeff] += term[icoeff]*T->ytar;
		//}
	      //	      b_pinv[icoeff] += term[icoeff]/T->p;
	      b_pinv[icoeff] += term[icoeff]*T->p*thetabend;

	      b_xfp[icoeff] += fterm[icoeff]*T->xfp;
	      b_yfp[icoeff] += fterm[icoeff]*T->yfp;
	      b_xpfp[icoeff] += fterm[icoeff]*T->xpfp;
	      b_ypfp[icoeff] += fterm[icoeff]*T->ypfp;
	      icoeff++;
	    }
	  }
	}
      }
    }

    for( icoeff=0; icoeff<ncoeff_optics; icoeff++ ){
      for( int jcoeff=0; jcoeff<ncoeff_optics; jcoeff++ ){
	Moptics(icoeff,jcoeff) += term[icoeff]*term[jcoeff];
	Mforward(icoeff,jcoeff) += fterm[icoeff]*fterm[jcoeff];
      }
    }

    //gamma = E/m = sqrt(p^2+m^2)/m = sqrt(1+p^2/m^2)
    double gamma = sqrt(1.0 + pow(T->p/Mp,2));
    double chi = gamma*kappa_p*thetabend; //
    double thetaspin = chi + thetabend;
    
    //compute the spin transport matrix in dipole approximation:
    //double chi_transport = T->chi + atan(T->xptar) - atan(T->xpfp);
    //T->chi is defined as gamma * kappa_p * (SBS_thetabend + atan(xptar)-atan(xpfp))
    // We are thus expanding the deviations of the actual spin transport relative to the dipole approximation
    // double Sxx_dipole = cos( chi_transport );
    // double Sxy_dipole = 0.0;
    // double Sxz_dipole = -sin( chi_transport );
    // double Syx_dipole = 0.0;
    // double Syy_dipole = 1.0;
    // double Syz_dipole =

    //    bend_axis.Print();
    
    TRotation Rdipole;
    //Rdipole.RotateY( -chi_transpo;
    Rdipole.Rotate( thetaspin, bend_axis );

    // cout << Rdipole.XX() << ", " << Rdipole.XY() << ", " << Rdipole.XZ() << endl
    // 	 << Rdipole.YX() << ", " << Rdipole.YY() << ", " << Rdipole.YZ() << endl
    // 	 << Rdipole.ZX() << ", " << Rdipole.ZY() << ", " << Rdipole.ZZ() << endl;
    
    Rdipole.RotateY( SBS_thetabend );
    
    // cout << Rdipole.XX() << ", " << Rdipole.XY() << ", " << Rdipole.XZ() << endl
    // 	 << Rdipole.YX() << ", " << Rdipole.YY() << ", " << Rdipole.YZ() << endl
    // 	 << Rdipole.ZX() << ", " << Rdipole.ZY() << ", " << Rdipole.ZZ() << endl;
      
    if( T->Pxtg == 1. ){ //what should the expansion for Spin matrix elements look like?
      // Here we are expanding the deviations from the dipole approximation in powers of
      // xptar, yptar, ytar, and 1/p. Is 1/p really the best expansion variable for the spin transport matrices?
      // Naively it seems like yes, it is.
      vector<double> term(ncoeff_spin);
      int icoeff=0;
      for( int i=0; i<=order_spin; i++ ){
	for( int j=0; j<=order_spin-i; j++ ){
	  for( int k=0; k<=order_spin-i-j; k++ ){
	    for( int l=0; l<=order_spin-i-j-k; l++ ){
	      for( int m=0; m<=order_spin-i-j-k-l; m++ ){
		term[icoeff] = pow(T->xptar,m)*pow(T->yptar,l)*pow(T->xtar,i)*pow(T->ytar,k)*pow(1.0/T->p,j);
		b_Sxx[icoeff] += (T->Pxfp - Rdipole.XX() ) * term[icoeff];
		b_Syx[icoeff] += (T->Pyfp - Rdipole.YX() ) * term[icoeff];
		b_Szx[icoeff] += (T->Pzfp - Rdipole.ZX() ) * term[icoeff];
		icoeff++;
	      }
	    }
	  }
	}
      }

      for( icoeff=0; icoeff<ncoeff_spin; icoeff++ ){
	for( int jcoeff=0; jcoeff<ncoeff_spin; jcoeff++ ){
	  Mspinx(icoeff,jcoeff) += term[icoeff]*term[jcoeff];
	}
      }
    } else if( T->Pytg == 1. ){
      vector<double> term(ncoeff_spin);
      int icoeff=0;
      for( int i=0; i<=order_spin; i++ ){
	for( int j=0; j<=order_spin-i; j++ ){
	  for( int k=0; k<=order_spin-i-j; k++ ){
	    for( int l=0; l<=order_spin-i-j-k; l++ ){
	      for( int m=0; m<=order_spin-i-j-k-l; m++ ){
		term[icoeff] = pow(T->xptar,m)*pow(T->yptar,l)*pow(T->xtar,i)*pow(T->ytar,k)*pow(1.0/T->p,j);
		b_Sxy[icoeff] += (T->Pxfp - Rdipole.XY() ) * term[icoeff];
		b_Syy[icoeff] += (T->Pyfp - Rdipole.YY() ) * term[icoeff];
		b_Szy[icoeff] += (T->Pzfp - Rdipole.ZY() ) * term[icoeff];
		icoeff++;
	      }
	    }
	  }
	}
      }

      for( icoeff=0; icoeff<ncoeff_spin; icoeff++ ){
	for( int jcoeff=0; jcoeff<ncoeff_spin; jcoeff++ ){
	  Mspiny(icoeff,jcoeff) += term[icoeff]*term[jcoeff];
	}
      }
    } else if( T->Pztg == 1. ){
      vector<double> term(ncoeff_spin);
      int icoeff=0;
      for( int i=0; i<=order_spin; i++ ){
	for( int j=0; j<=order_spin-i; j++ ){
	  for( int k=0; k<=order_spin-i-j; k++ ){
	    for( int l=0; l<=order_spin-i-j-k; l++ ){
	      for( int m=0; m<=order_spin-i-j-k-l; m++ ){
		term[icoeff] = pow(T->xptar,m)*pow(T->yptar,l)*pow(T->xtar,i)*pow(T->ytar,k)*pow(1.0/T->p,j);
		b_Sxz[icoeff] += (T->Pxfp - Rdipole.XZ() ) * term[icoeff];
		b_Syz[icoeff] += (T->Pyfp - Rdipole.YZ() ) * term[icoeff];
		b_Szz[icoeff] += (T->Pzfp - Rdipole.ZZ() ) * term[icoeff];
		icoeff++;
	      }
	    }
	  }
	}
      }

      for( icoeff=0; icoeff<ncoeff_spin; icoeff++ ){
	for( int jcoeff=0; jcoeff<ncoeff_spin; jcoeff++ ){
	  Mspinz(icoeff,jcoeff) += term[icoeff]*term[jcoeff];
	}
      }
    }
  }

  if( use_xtar_flag == 0 ){ //Set all diagonal coefficients to 1 and all off-diagonal coefficients to zero for Moptics, AND set all b_XX elements to zero: 
    for( int icoeff=0; icoeff<ncoeff_optics; icoeff++){
      if( xtar_expon_optics[icoeff] > 0 ){
	Moptics(icoeff,icoeff)=1.0;
	Mforward(icoeff,icoeff)=1.0;
	for(int jcoeff=0; jcoeff<ncoeff_optics; jcoeff++){
	  if( jcoeff != icoeff ) {
	    Moptics(icoeff,jcoeff) = 0.0;
	    Mforward(icoeff,jcoeff) = 0.0;
	  }
	}
	b_xptar[icoeff] = 0.0;
	b_yptar[icoeff] = 0.0;
	b_ytar[icoeff] = 0.0;
	b_pinv[icoeff] = 0.0;

	b_xfp[icoeff] = 0.0;
	b_yfp[icoeff] = 0.0;
	b_xpfp[icoeff] = 0.0;
	b_ypfp[icoeff] = 0.0;
      }
    }

    for( int icoeff=0; icoeff<ncoeff_spin; icoeff++ ){
      if( xtar_expon_spin[icoeff] > 0 ){
	Mspinx(icoeff,icoeff) = 1.0;
	Mspiny(icoeff,icoeff) = 1.0;
	Mspinz(icoeff,icoeff) = 1.0;
	for( int jcoeff=0; jcoeff<ncoeff_spin; jcoeff++ ){
	  if( jcoeff != icoeff ){
	    Mspinx(icoeff,jcoeff) = 0.0;
	    Mspiny(icoeff,jcoeff) = 0.0;
	    Mspinz(icoeff,jcoeff) = 0.0;
	  }
	}
	b_Sxx(icoeff) = 0.0;
	b_Syx(icoeff) = 0.0;
	b_Szx(icoeff) = 0.0;
	
	b_Sxy(icoeff) = 0.0;
	b_Syy(icoeff) = 0.0;
	b_Szy(icoeff) = 0.0;
	
	b_Sxz(icoeff) = 0.0;
	b_Syz(icoeff) = 0.0;
	b_Szz(icoeff) = 0.0;
	
      }
    }
  }
  
  TDecompSVD A_xptar(Moptics);
  TDecompSVD A_pinv(Moptics);
  
  // for( int icoeff=0; icoeff<ncoeff_optics; icoeff++){
  //   if( xtar_expon_optics[icoeff] > 0 ){
  //     Moptics(icoeff,icoeff) = 1.0;
  //     for( int jcoeff=0; jcoeff<ncoeff_optics; jcoeff++ ){
  // 	if( jcoeff != icoeff ) Moptics(icoeff,jcoeff) = 0.0;
  //     }
  //     b_ytar(icoeff) = 0.0;
  //     b_yptar(icoeff) = 0.0;
  //   }
  // }

  TDecompSVD A_yptar(Moptics);
  TDecompSVD A_ytar(Moptics);

  cout << "solving xptar" << endl;
  bool good_xptar = A_xptar.Solve(b_xptar);
  cout << "xptar done, success = " << good_xptar << endl;
  cout << "solving yptar" << endl;
  bool good_yptar = A_yptar.Solve(b_yptar);
  cout << "yptar done, success = " << good_yptar << endl;
  cout << "solving ytar" << endl;
  bool good_ytar = A_ytar.Solve(b_ytar);
  cout << "ytar done, success = " << good_ytar << endl;
  cout << "solving pinv" << endl;
  bool good_pinv = A_pinv.Solve(b_pinv);
  cout << "pinv done, success = " << good_pinv << endl;

  TDecompSVD A_xfp(Mforward);
  TDecompSVD A_yfp(Mforward);
  TDecompSVD A_xpfp(Mforward);
  TDecompSVD A_ypfp(Mforward);

  A_xfp.Solve(b_xfp);
  A_yfp.Solve(b_yfp);
  A_xpfp.Solve(b_xpfp);
  A_ypfp.Solve(b_ypfp);
  
  TDecompSVD A_Sxx(Mspinx);
  TDecompSVD A_Syx(Mspinx);
  TDecompSVD A_Szx(Mspinx);

  A_Sxx.Solve( b_Sxx );
  cout << "Sxx done" << endl;
  A_Syx.Solve( b_Syx );
  cout << "Syx done" << endl;
  A_Szx.Solve( b_Szx );
  cout << "Szx done" << endl;
  
  TDecompSVD A_Sxy(Mspiny);
  TDecompSVD A_Syy(Mspiny);
  TDecompSVD A_Szy(Mspiny);

  A_Sxy.Solve( b_Sxy );
  cout << "Sxy done" << endl;
  A_Syy.Solve( b_Syy );
  cout << "Syy done" << endl;
  A_Szy.Solve( b_Szy );
  cout << "Szy done" << endl;
  
  TDecompSVD A_Sxz(Mspinz);
  TDecompSVD A_Syz(Mspinz);
  TDecompSVD A_Szz(Mspinz);

  A_Sxz.Solve( b_Sxz );
  cout << "Sxz done" << endl;
  A_Syz.Solve( b_Syz );
  cout << "Syz done" << endl;
  A_Szz.Solve( b_Szz );
  cout << "Szz done" << endl;
  
  TFile *fout = new TFile(outfilename,"RECREATE");
  TTree *Tout = new TTree("Tout","SBS optics fit result");

  double xfp,yfp,xpfp,ypfp,xfprecon,yfprecon,xpfprecon,ypfprecon;
  double xptar,yptar,ytar,xtar,p,thetabend;
  double xptarrecon,yptarrecon,ytarrecon,xtarrecon,precon,thetabendrecon;
  double xfpforward,yfpforward,xpfpforward,ypfpforward;

  Tout->Branch("xfp",&xfp);
  Tout->Branch("yfp",&yfp);
  Tout->Branch("xpfp",&xpfp);
  Tout->Branch("ypfp",&ypfp);
  Tout->Branch("xfprecon",&xfprecon);
  Tout->Branch("yfprecon",&yfprecon);
  Tout->Branch("xpfprecon",&xpfprecon);
  Tout->Branch("ypfprecon",&ypfprecon);
  Tout->Branch("xptar",&xptar);
  Tout->Branch("yptar",&yptar);
  Tout->Branch("ytar",&ytar);
  Tout->Branch("xtar",&xtar);
  Tout->Branch("p",&p);
  Tout->Branch("thetabend",&thetabend);
  Tout->Branch("xptarrecon",&xptarrecon);
  Tout->Branch("yptarrecon",&yptarrecon);
  Tout->Branch("ytarrecon",&ytarrecon);
  Tout->Branch("xtarrecon",&xtarrecon);
  Tout->Branch("precon",&precon);
  Tout->Branch("thetabendrecon",&thetabendrecon);
  Tout->Branch("xfpforward",&xfpforward);
  Tout->Branch("yfpforward",&yfpforward);
  Tout->Branch("xpfpforward",&xpfpforward);
  Tout->Branch("ypfpforward",&ypfpforward);

  
  double Pxtg, Pytg, Pztg;
  double Pxfp, Pyfp, Pzfp;
  double Pxfpdipole,Pyfpdipole,Pzfpdipole;
  double Pxfprecon,Pyfprecon,Pzfprecon;
  Tout->Branch("Pxtg",&Pxtg);
  Tout->Branch("Pytg",&Pytg);
  Tout->Branch("Pztg",&Pztg);

  Tout->Branch("Pxfp",&Pxfp);
  Tout->Branch("Pyfp",&Pyfp);
  Tout->Branch("Pzfp",&Pzfp);

  Tout->Branch("Pxfprecon",&Pxfprecon);
  Tout->Branch("Pyfprecon",&Pyfprecon);
  Tout->Branch("Pzfprecon",&Pzfprecon);

  Tout->Branch("Pxfpdipole",&Pxfpdipole);
  Tout->Branch("Pyfpdipole",&Pyfpdipole);
  Tout->Branch("Pzfpdipole",&Pzfpdipole);
  double det;
  Tout->Branch("Determinant",&det);

  TH1D *hxptardiff = new TH1D("hxptardiff","",250,-0.01,0.01);
  TH1D *hyptardiff = new TH1D("hyptardiff","",250,-0.01,0.01);
  TH1D *hytardiff = new TH1D("hytardiff","",250,-0.02,0.02);
  TH1D *hpdiff = new TH1D("hpdiff","",250,-0.1,0.1);
  TH1D *hthetabend_diff = new TH1D("hthetabend_diff","",250,-.01,.01);

  TH1D *hxfpdiff = new TH1D("hxfpdiff","",250,-0.02,0.02);
  TH1D *hyfpdiff = new TH1D("hyfpdiff","",250,-0.02,0.02);
  TH1D *hxpfpdiff = new TH1D("hxpfpdiff","",250,-0.01,0.01);
  TH1D *hypfpdiff = new TH1D("hypfpdiff","",250,-0.01,0.01);
  
  
  TH2D *hxptar_recon_vs_true = new TH2D("hxptar_recon_vs_true","",250,-0.3,0.3,250,-0.3,0.3);
  TH2D *hyptar_recon_vs_true = new TH2D("hyptar_recon_vs_true","",250,-0.1,0.1,250,-0.1,0.1);
  TH2D *hytar_recon_vs_true = new TH2D("hytar_recon_vs_true","",250,-0.1,0.1,250,-0.1,0.1);
  TH2D *hp_recon_vs_true = new TH2D("hp_recon_vs_true","",250,0.0,12.5,250,0.0,12.5);
  TH2D *hthetabend_recon_vs_true = new TH2D("hthetabend_recon_vs_true","",250,0.0,0.5,250,0.0,0.5);
  TH2D *hpdiff_p = new TH2D("hpdiff_p","",250,0.0,12.5,250,-0.1,0.1);

  
  TH2D *hxfp_recon_vs_true = new TH2D("hxfp_recon_vs_true","",250,-0.8,0.8,250,-0.8,0.8);
  TH2D *hyfp_recon_vs_true = new TH2D("hyfp_recon_vs_true","",250,-0.25,0.25,250,-0.25,0.25);
  TH2D *hxpfp_recon_vs_true = new TH2D("hxpfp_recon_vs_true","",250,-0.4,0.4,250,-0.4,0.4);
  TH2D *hypfp_recon_vs_true = new TH2D("hypfp_recon_vs_true","",250,-0.1,0.1,250,-0.1,0.1);
  
  
  TH1D *hPxfpdiff = new TH1D("hPxfpdiff","",250,-0.05,0.05);
  TH1D *hPyfpdiff = new TH1D("hPyfpdiff","",250,-0.05,0.05);
  TH1D *hPzfpdiff = new TH1D("hPzfpdiff","",250,-0.05,0.05);

  TH2D *hPxfp_recon_vs_true = new TH2D("hPxfp_recon_vs_true","",250,-1.05,1.05,250,-1.05,1.05);
  TH2D *hPyfp_recon_vs_true = new TH2D("hPyfp_recon_vs_true","",250,-1.05,1.05,250,-1.05,1.05);
  TH2D *hPzfp_recon_vs_true = new TH2D("hPzfp_recon_vs_true","",250,-1.05,1.05,250,-1.05,1.05);

  TH1D *hSxx = new TH1D("hSxx","",250,-1.1,1.1);
  TH1D *hSxy = new TH1D("hSxy","",250,-1.1,1.1);
  TH1D *hSxz = new TH1D("hSxz","",250,-1.1,1.1);

  TH1D *hSyx = new TH1D("hSyx","",250,-1.1,1.1);
  TH1D *hSyy = new TH1D("hSyy","",250,-1.1,1.1);
  TH1D *hSyz = new TH1D("hSyz","",250,-1.1,1.1);

  TH1D *hSzx = new TH1D("hSzx","",250,-1.1,1.1);
  TH1D *hSzy = new TH1D("hSzy","",250,-1.1,1.1);
  TH1D *hSzz = new TH1D("hSzz","",250,-1.1,1.1);

  TH1D *hDeterminant = new TH1D("hDeterminant","",250,0.95,1.05);
  
  TString opticsfilename = outfilename;
  opticsfilename.ReplaceAll(".root","_optics_coeff.txt");

  TString spinfilename = outfilename;
  spinfilename.ReplaceAll(".root","_spin_coeff.txt");

  //Write "self-documenting" optics database file:
  ofstream opticsfile(opticsfilename.Data());

  TString optics_header;

  opticsfile << "# number of nonzero terms = " << endl
	     << ncoeff_optics_nonzero << endl;
  opticsfile << "# Expansion is: u = sum_{i,j,k,l,m}=0^{i+j+k+l+m<=order} Cu_{ijklm} xfp^i * yfp^j * xpfp^k * ypfp^l * xtar^m, " << endl
	     << "# u = xptar, yptar, ytar, p*thetabend. Columns are:" << endl; 
  optics_header.Form( "# %15s %15s %15s %15s   i j k l m",
		      "Cxptar", "Cyptar", "Cytar", "Cpthetabend" );
  
  opticsfile << optics_header << endl;

  int ordertemp = 0;

  
  int ipar=0;
  while( ordertemp <= order_optics ){
    ipar=0;
    for( int i=0; i<=order_optics; i++ ){
      for( int j=0; j<=order_optics-i; j++ ){
	for( int k=0; k<=order_optics-i-j; k++ ){
	  for( int l=0; l<=order_optics-i-j-k; l++ ){
	    for( int m=0; m<=order_optics-i-j-k-l; m++ ){
	      TString opticsline;
	      opticsline.Form("  %15.8g %15.8g %15.8g %15.8g   %d %d %d %d %d",
			    b_xptar(ipar), b_yptar(ipar),
			    b_ytar(ipar), b_pinv(ipar),
			    m, l, k, j, i );
	      if( i+j+k+l+m == ordertemp && (i==0 || use_xtar_flag > 0 )){
		opticsfile << opticsline << endl;
	      }
	      ipar++;
	    }
	  }
	}
      }
    }
    ordertemp++;
  }

  TString fopticsfilename = outfilename;
  fopticsfilename.ReplaceAll(".root","_forward_optics_coeff.txt");
  //Write "self-documenting" optics database file:
  ofstream fopticsfile(fopticsfilename.Data());

  TString foptics_header;

  fopticsfile << "# number of nonzero terms = " << endl
	      << ncoeff_optics_nonzero << endl;
  fopticsfile << "# Expansion is: u = sum_{i,j,k,l,m}=0^{i+j+k+l+m<=order} Cu_{ijklm} xptar^i * yptar^j * ytar^k * (1/p)^l * xtar^m, " << endl
	     << "# u = xfp, yfp, xpfp, ypfp. Columns are:" << endl; 
  optics_header.Form( "# %15s %15s %15s %15s   i j k l m",
		      "Cxfp", "Cyfp", "Cxpfp", "Cypfp" );
  
  fopticsfile << optics_header << endl;

  ipar=0;

  ordertemp=0;

  while( ordertemp <= order_optics ){
    ipar=0;
    for( int i=0; i<=order_optics; i++ ){
      for( int j=0; j<=order_optics-i; j++ ){
	for( int k=0; k<=order_optics-i-j; k++ ){
	  for( int l=0; l<=order_optics-i-j-k; l++ ){
	    for( int m=0; m<=order_optics-i-j-k-l; m++ ){
	      TString opticsline;
	      opticsline.Form("  %15.8g %15.8g %15.8g %15.8g   %d %d %d %d %d",
			      b_xfp(ipar), b_yfp(ipar),
			      b_xpfp(ipar), b_ypfp(ipar),
			      m, l, k, j, i );

	      if( i+j+k+l+m == ordertemp && (i==0 || use_xtar_flag > 0 )){
		
		fopticsfile << opticsline << endl;
	      }
	      ipar++;
	    }
	  }
	}
      }
    }
    ordertemp++;
  }

  
  
  
  ofstream spinfile(spinfilename.Data());

  spinfile << "# number of nonzero terms = " << endl
	   << ncoeff_spin_nonzero << endl;
  
  spinfile << "# Expansion is: S_{ab} = sum_{i,j,k,l,m}=0^{i+j+k+l+m<=order} Cab_{ijklm} xptar^i yptar^j ytar^k (1/p)^l xtar^m, " << endl
	   << "# ab = xx,xy,xz,yx,yy,yz,zx,zy,zz: Columns are:" << endl;

  TString spin_header;
  spin_header.Form("# %15s %15s %15s %15s %15s %15s %15s %15s %15s   i j k l m",
		   "Cxx", "Cxy", "Cxz",
		   "Cyx", "Cyy", "Cyz",
		   "Czx", "Czy", "Czz" );

  spinfile << spin_header << endl;

  ordertemp = 0;
  while(ordertemp <= order_spin ){
    ipar=0;
    for( int i=0; i<=order_spin; i++ ){
      for( int j=0; j<=order_spin-i; j++ ){
	for( int k=0; k<=order_spin-i-j; k++ ){
	  for( int l=0; l<=order_spin-i-j-k; l++ ){
	    for( int m=0; m<=order_spin-i-j-k-l; m++ ){
	      TString spinline;
	      spinline.Form("  %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g %15.8g   %d %d %d %d %d",
			    b_Sxx(ipar), b_Sxy(ipar), b_Sxz(ipar),
			    b_Syx(ipar), b_Syy(ipar), b_Syz(ipar),
			    b_Szx(ipar), b_Szy(ipar), b_Szz(ipar),
			    m,l,k,j,i );
	      if( i+j+k+l+m == ordertemp && (i == 0 || use_xtar_flag > 0 ) ){
		spinfile << spinline << endl;
	      }
	      
	      ipar++;
	    }
	  }
	}
      }
    }
    ordertemp++;
  }

  //Do reconstruction again and histogram results of fit:
  
  nevent=0;
  while( T->GetEntry(nevent++) ){
    if( nevent%100 == 0 ) cout << nevent << endl;
    xfp = T->xfp;
    yfp = T->yfp;
    xpfp = T->xpfp;
    ypfp = T->ypfp;
    xfprecon = T->xfprecon;
    yfprecon = T->yfprecon;
    xpfprecon = T->xpfprecon;
    ypfprecon = T->ypfprecon;
    
    xptar = T->xptar;
    yptar = T->yptar;
    ytar = T->ytar;
    xtar = T->xtar;
    p = T->p;

    Pxtg = T->Pxtg;
    Pytg = T->Pytg;
    Pztg = T->Pztg;
    Pxfp = T->Pxfp;
    Pyfp = T->Pyfp;
    Pzfp = T->Pzfp;
    
    TVector3 zaxis_FT(-sin(SBS_thetabend),0,cos(SBS_thetabend));
    TVector3 yaxis_FT(0,1,0);
    TVector3 xaxis_FT = yaxis_FT.Cross(zaxis_FT).Unit();

    TVector3 nhat_fp(xpfp,ypfp,1.0);
    nhat_fp = nhat_fp.Unit();
    TVector3 nhat_fp_global = nhat_fp.X() * xaxis_FT +
      nhat_fp.Y() * yaxis_FT +
      nhat_fp.Z() * zaxis_FT;

    TVector3 nhat_tgt(xptar,yptar,1.0);

    nhat_tgt = nhat_tgt.Unit();

    thetabend = acos(nhat_fp_global.Dot(nhat_tgt) );

    TVector3 bend_axis_track = nhat_tgt.Cross(nhat_fp_global).Unit();
    TRotation Rtrack;

    Rtrack.Rotate( thetabend, bend_axis_track );
    
    //compute the spin transport matrix in dipole approximation:
    //double chi_transport = T->chi + atan(T->xptar) - atan(T->xpfp);
    // double Sxx_dipole = cos( chi_transport );
    // double Sxy_dipole = 0.0;
    // double Sxz_dipole = -sin( chi_transport );
    // double Syx_dipole = 0.0;
    // double Syy_dipole = 1.0;
    // double Syz_dipole =
    
    double pinv_sum = 0.0;
    double xptar_sum = 0.0;
    double yptar_sum = 0.0;
    double ytar_sum = 0.0;

    double xfp_sum = 0.0;
    double yfp_sum = 0.0;
    double xpfp_sum = 0.0;
    double ypfp_sum = 0.0;
    
    int icoeff=0;
    for( int i=0; i<=order_optics; i++ ){
      for( int j=0; j<=order_optics-i; j++ ){
	for( int k=0; k<=order_optics-i-j; k++ ){
	  for( int l=0; l<=order_optics-i-j-k; l++ ){
	    for( int m=0; m<=order_optics-i-j-k-l; m++ ){
	      //double term = pow(T->xfp,m)*pow(T->yfp,l)*pow(T->xpfp,k)*pow(T->ypfp,j)*pow(T->xtar,i);
	      //use smeared fp quantities to account for effects of tracking resolution
	      double term = pow(xfprecon,m)*pow(yfprecon,l)*pow(xpfprecon,k)*pow(ypfprecon,j)*pow(T->xtar,i);    
	      xptar_sum += term * b_xptar[icoeff];
	      yptar_sum += term * b_yptar[icoeff];
	      ytar_sum += term * b_ytar[icoeff];
	      pinv_sum += term * b_pinv[icoeff];

	      double fterm = pow(T->xptar,m)*pow(T->yptar,l)*pow(T->ytar,k)*pow(1.0/T->p,j)*pow(T->xtar,i);
	      xfp_sum += fterm * b_xfp[icoeff];
	      yfp_sum += fterm * b_yfp[icoeff];
	      xpfp_sum += fterm * b_xpfp[icoeff];
	      ypfp_sum += fterm * b_ypfp[icoeff];
	      
	      icoeff++;
	    }
	  }
	}
      }
    }

    xfpforward = xfp_sum;
    yfpforward = yfp_sum;
    xpfpforward = xpfp_sum;
    ypfpforward = ypfp_sum;

    hxfpdiff->Fill(xfpforward-xfp);
    hyfpdiff->Fill(yfpforward-yfp);
    hxpfpdiff->Fill(xpfpforward-xpfp);
    hypfpdiff->Fill(ypfpforward-ypfp);

    hxfp_recon_vs_true->Fill(xfp,xfpforward);
    hyfp_recon_vs_true->Fill(yfp,yfpforward);
    hxpfp_recon_vs_true->Fill(xpfp,xpfpforward);
    hypfp_recon_vs_true->Fill(ypfp,ypfpforward);
 
    xptarrecon = xptar_sum;
    yptarrecon = yptar_sum;
    ytarrecon = ytar_sum;
    //precon = 1.0/pinv_sum;
    //now pinv_sum = p*thetabend, so we need to compute the reconstructed and true thetabend!

    hxptardiff->Fill(xptarrecon-xptar);
    hyptardiff->Fill(yptarrecon-yptar);
    hytardiff->Fill(ytarrecon-ytar);

    hxptar_recon_vs_true->Fill(xptar,xptarrecon);
    hyptar_recon_vs_true->Fill(yptar,yptarrecon);
    hytar_recon_vs_true->Fill(ytar,ytarrecon);
    
    TVector3 nhat_fp_recon(xpfprecon,ypfprecon,1.);
    nhat_fp_recon = nhat_fp_recon.Unit();

    TVector3 nhat_fp_recon_global = nhat_fp_recon.X() * xaxis_FT +
      nhat_fp_recon.Y() * yaxis_FT +
      nhat_fp_recon.Z() * zaxis_FT;

    TVector3 nhat_tgt_recon(xptarrecon,yptarrecon,1.);
    nhat_tgt_recon = nhat_tgt_recon.Unit();

    thetabendrecon = acos( nhat_tgt_recon.Dot(nhat_fp_recon_global) );
    
    precon = pinv_sum/thetabendrecon;

    hpdiff->Fill(precon/T->p-1.0);
    hp_recon_vs_true->Fill(T->p,precon);
    hpdiff_p->Fill( T->p, precon/T->p-1.0);
    
    hthetabend_diff->Fill(thetabendrecon-thetabend);
    hthetabend_recon_vs_true->Fill(thetabend,thetabendrecon);
    
    TVector3 bend_axis_recon = nhat_tgt_recon.Cross( nhat_fp_recon_global ).Unit();

    double gammarecon = sqrt(1.0+pow(precon/Mp,2));
    double chirecon = gammarecon*kappa_p*thetabendrecon;

    double thetaspinrecon = chirecon + thetabendrecon;

    TRotation Rdipole;
    Rdipole.Rotate( thetaspinrecon, bend_axis_recon );
    Rdipole.RotateY( SBS_thetabend );
    
    double Sxx_sum = Rdipole.XX(), Sxy_sum = Rdipole.XY(), Sxz_sum = Rdipole.XZ();
    double Syx_sum = Rdipole.YX(), Syy_sum = Rdipole.YY(), Syz_sum = Rdipole.YZ();
    double Szx_sum = Rdipole.ZX(), Szy_sum = Rdipole.ZY(), Szz_sum = Rdipole.ZZ();

    icoeff=0;
    for( int i=0; i<=order_spin; i++ ){
      for( int j=0; j<=order_spin-i; j++ ){
	for( int k=0; k<=order_spin-i-j; k++ ){
	  for( int l=0; l<=order_spin-i-j-k; l++ ){
	    for( int m=0; m<=order_spin-i-j-k-l; m++ ){
	      //double term = pow(T->xptar,m)*pow(T->yptar,l)*pow(T->xtar,k)*pow(T->ytar,j)*pow(1.0/T->p,i);
	      double term = pow(xptarrecon,m)*pow(yptarrecon,l)*pow(T->xtar,i)*pow(T->ytar,k)*pow(1.0/precon,j);
	      Sxx_sum += b_Sxx[icoeff] * term;
	      Syx_sum += b_Syx[icoeff] * term;
	      Szx_sum += b_Szx[icoeff] * term;

	      Sxy_sum += b_Sxy[icoeff] * term;
	      Syy_sum += b_Syy[icoeff] * term;
	      Szy_sum += b_Szy[icoeff] * term;

	      Sxz_sum += b_Sxz[icoeff] * term;
	      Syz_sum += b_Syz[icoeff] * term;
	      Szz_sum += b_Szz[icoeff] * term;

	      icoeff++;
	    }
	  }
	}
      }
    }

    TMatrixD Rot(3,3);
    Rot(0,0) = Sxx_sum; Rot(0,1) = Sxy_sum; Rot(0,2) = Sxz_sum;
    Rot(1,0) = Syx_sum; Rot(1,1) = Syy_sum; Rot(1,2) = Syz_sum;
    Rot(2,0) = Szx_sum; Rot(2,1) = Szy_sum; Rot(2,2) = Szz_sum;

    TVectorD Ptg(3);
    Ptg(0) = Pxtg;
    Ptg(1) = Pytg;
    Ptg(2) = Pztg;

    TVectorD Pfp(3);
    Pfp = Rot * Ptg;
    Pxfprecon = Pfp(0);
    Pyfprecon = Pfp(1);
    Pzfprecon = Pfp(2);

    TVector3 Pfpdipole(Pxtg,Pytg,Pztg);
    Pfpdipole *= Rdipole;

    // Pxfpdipole = Pfpdipole.Dot( xaxis_FT );
    // Pyfpdipole = Pfpdipole.Dot( yaxis_FT );
    // Pzfpdipole = Pfpdipole.Dot( zaxis_FT );

    Pxfpdipole = Pfpdipole.X();
    Pyfpdipole = Pfpdipole.Y();
    Pzfpdipole = Pfpdipole.Z();
    
    //This is in the rotated coordinate system: Transform back to FP coordinates:
    

    hPxfpdiff->Fill(Pxfprecon-Pxfp);
    hPyfpdiff->Fill(Pyfprecon-Pyfp);
    hPzfpdiff->Fill(Pzfprecon-Pzfp);

    hPxfp_recon_vs_true->Fill(Pxfp,Pxfprecon);
    hPyfp_recon_vs_true->Fill(Pyfp,Pyfprecon);
    hPzfp_recon_vs_true->Fill(Pzfp,Pzfprecon);

    hSxx->Fill( Sxx_sum );
    hSxy->Fill( Sxy_sum );
    hSxz->Fill( Sxz_sum );

    hSyx->Fill( Syx_sum );
    hSyy->Fill( Syy_sum );
    hSyz->Fill( Syz_sum );

    hSzx->Fill( Szx_sum );
    hSzy->Fill( Szy_sum );
    hSzz->Fill( Szz_sum );

    //cout << "Determinant = " << Rot.Determinant() << endl;

    det = Rot.Determinant();

    hDeterminant->Fill(det);
    
    Tout->Fill();
  }

  hxptardiff->GetXaxis()->SetTitle("#Deltax'_{tar}");
  hyptardiff->GetXaxis()->SetTitle("#Deltay'_{tar}");
  hytardiff->GetXaxis()->SetTitle("#Deltay_{tar} (m)");
  hpdiff->GetXaxis()->SetTitle("p_{recon}/p_{true}-1");

  hthetabend_diff->GetXaxis()->SetTitle("#Delta#theta_{bend} (rad)");

  hxfpdiff->GetXaxis()->SetTitle("#Deltax_{fp} (m)");
  hyfpdiff->GetXaxis()->SetTitle("#Deltay_{fp} (m)");
  hxpfpdiff->GetXaxis()->SetTitle("#Deltax'_{fp}" );
  hypfpdiff->GetXaxis()->SetTitle("#Deltay'_{fp}" );

  hxptar_recon_vs_true->GetXaxis()->SetTitle("x'_{tar} (true)");
  hxptar_recon_vs_true->GetYaxis()->SetTitle("x'_{tar} (reconstructed)");
  
  hyptar_recon_vs_true->GetXaxis()->SetTitle("y'_{tar} (true)");
  hyptar_recon_vs_true->GetYaxis()->SetTitle("y'_{tar} (reconstructed)");

  hytar_recon_vs_true->GetXaxis()->SetTitle("y_{tar}^{true} (m)");
  hytar_recon_vs_true->GetYaxis()->SetTitle("y_{tar}^{recon} (m)");

  hp_recon_vs_true->GetXaxis()->SetTitle("p_{true} (GeV)");
  hp_recon_vs_true->GetYaxis()->SetTitle("p_{recon} (GeV)");

  hthetabend_recon_vs_true->GetXaxis()->SetTitle("#theta_{bend}^{true} (rad)");
  hthetabend_recon_vs_true->GetYaxis()->SetTitle("#theta_{bend}^{recon} (rad)");

  hpdiff_p->GetXaxis()->SetTitle("p_{true} (GeV)");
  hpdiff_p->GetYaxis()->SetTitle("p_{recon}/p_{true}-1 (%)");

  hxfp_recon_vs_true->GetXaxis()->SetTitle("x_{fp}^{true} (m)");
  hxfp_recon_vs_true->GetYaxis()->SetTitle("x_{fp}^{recon} (m)");

  hyfp_recon_vs_true->GetXaxis()->SetTitle("y_{fp}^{true} (m)");
  hyfp_recon_vs_true->GetYaxis()->SetTitle("y_{fp}^{recon} (m)");

  hxpfp_recon_vs_true->GetXaxis()->SetTitle("(x'_{fp})_{true} (m)");
  hxpfp_recon_vs_true->GetYaxis()->SetTitle("(x'_{fp})_{recon} (m)");

  hypfp_recon_vs_true->GetXaxis()->SetTitle("(y'_{fp})_{true} (m)");
  hypfp_recon_vs_true->GetYaxis()->SetTitle("(y'_{fp})_{recon} (m)");

  hPxfpdiff->GetXaxis()->SetTitle("#DeltaP_{x}");
  hPyfpdiff->GetXaxis()->SetTitle("#DeltaP_{y}");
  hPzfpdiff->GetXaxis()->SetTitle("#DeltaP_{z}");

  hPxfp_recon_vs_true->GetXaxis()->SetTitle("P_{x}^{true}");
  hPyfp_recon_vs_true->GetXaxis()->SetTitle("P_{y}^{true}");
  hPzfp_recon_vs_true->GetXaxis()->SetTitle("P_{z}^{true}");

  hPxfp_recon_vs_true->GetYaxis()->SetTitle("P_{x}^{recon}");
  hPyfp_recon_vs_true->GetYaxis()->SetTitle("P_{y}^{recon}");
  hPzfp_recon_vs_true->GetYaxis()->SetTitle("P_{z}^{recon}");

  hSxx->GetXaxis()->SetTitle("S_{xx}");
  hSxy->GetXaxis()->SetTitle("S_{xy}");
  hSxz->GetXaxis()->SetTitle("S_{xz}");

  hSyx->GetXaxis()->SetTitle("S_{yx}");
  hSyy->GetXaxis()->SetTitle("S_{yy}");
  hSyz->GetXaxis()->SetTitle("S_{yz}");

  hSzx->GetXaxis()->SetTitle("S_{zx}");
  hSzy->GetXaxis()->SetTitle("S_{zy}");
  hSzz->GetXaxis()->SetTitle("S_{zz}");
  
  fout->Write();

  hDeterminant->GetXaxis()->SetTitle("Determinant");
  
}
