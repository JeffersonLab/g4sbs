#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "gep_optics_tree.C"
#include "TDecompSVD.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include <vector>
#include <iostream>
#include "TString.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRotation.h"

//double PI = TMath::Pi();
const double Mp = 0.938272046; //GeV/c^2
const double mu_p = 2.792847356; //mu_N
const double kappa_p = mu_p - 1.0; //anomalous magnetic moment
const double PI = TMath::Pi();
const double SBS_thetabend = 5.0*PI/180.0; //central bend angle

//const double SBS_thetabend = 5.0*PI/180.0;

void Fit_SBS_optics_and_spin( const char *rootfilename, int order_optics=4, int order_spin=4, const char *outfilename="SBS_optics_fitresult.root" ){
  int ncoeff_optics = 0;

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
	  }
	}
      }
    }
  }

  cout << "nterms optics fit = " << ncoeff_optics << endl;
  
  int ncoeff_spin = 0;

  vector<int> xtar_expon_spin, xfp_expon_spin, xpfp_expon_spin;

  for( int i=0; i<=order_spin; i++ ){
    for( int j=0; j<=order_spin-i; j++ ){
      for( int k=0; k<=order_spin-i-j; k++ ){
	for( int l=0; l<=order_spin-i-j-k; l++ ){
	  for( int m=0; m<=order_spin-i-j-k-l; m++ ){
	    ncoeff_spin++;
	    xtar_expon_spin.push_back( i );
	    xfp_expon_spin.push_back( m );
	    xpfp_expon_spin.push_back( k );
	  }
	}
      }
    }
  }

  TMatrixD Moptics(ncoeff_optics,ncoeff_optics);
  TVectorD b_xptar(ncoeff_optics), b_yptar(ncoeff_optics), b_ytar(ncoeff_optics), b_pinv(ncoeff_optics);

  for( int ipar=0; ipar<ncoeff_optics; ipar++ ){
    b_xptar(ipar) = 0.0;
    b_yptar(ipar) = 0.0;
    b_ytar(ipar) = 0.0;
    b_pinv(ipar) = 0.0;
    for( int jpar=0; jpar<ncoeff_optics; jpar++ ){
      Moptics(ipar,jpar) = 0.0;
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
    vector<double> term(ncoeff_optics);
    int icoeff=0;
    for( int i=0; i<=order_optics; i++ ){
      for( int j=0; j<=order_optics-i; j++ ){
	for( int k=0; k<=order_optics-i-j; k++ ){
	  for( int l=0; l<=order_optics-i-j-k; l++ ){
	    for( int m=0; m<=order_optics-i-j-k-l; m++ ){
	      term[icoeff] = pow(T->xfp,m)*pow(T->yfp,l)*pow(T->xpfp,k)*pow(T->ypfp,j)*pow(T->xtar,i);

	      b_xptar[icoeff] += term[icoeff]*T->xptar;
	      if( i == 0 ){ //don't include xtar-dependent terms in yptar, ytar fits:
		b_yptar[icoeff] += term[icoeff]*T->yptar;
		b_ytar[icoeff] += term[icoeff]*T->ytar;
	      }
	      b_pinv[icoeff] += term[icoeff]/T->p;
	      icoeff++;
	    }
	  }
	}
      }
    }

    for( icoeff=0; icoeff<ncoeff_optics; icoeff++ ){
      for( int jcoeff=0; jcoeff<ncoeff_optics; jcoeff++ ){
	Moptics(icoeff,jcoeff) += term[icoeff]*term[jcoeff];
      }
    }

    //compute the spin transport matrix in dipole approximation:
    double chi_transport = T->chi + atan(T->xptar) - atan(T->xpfp);
    // double Sxx_dipole = cos( chi_transport );
    // double Sxy_dipole = 0.0;
    // double Sxz_dipole = -sin( chi_transport );
    // double Syx_dipole = 0.0;
    // double Syy_dipole = 1.0;
    // double Syz_dipole =
    TRotation Rdipole;
    Rdipole.RotateY( -chi_transport );

    // cout << Rdipole.XX() << ", " << Rdipole.XY() << ", " << Rdipole.XZ() << endl
    // 	 << Rdipole.YX() << ", " << Rdipole.YY() << ", " << Rdipole.YZ() << endl
    // 	 << Rdipole.ZX() << ", " << Rdipole.ZY() << ", " << Rdipole.ZZ() << endl;
      
    if( T->Pxtg == 1. ){ //what should the expansion for Spin matrix elements look like? 
      vector<double> term(ncoeff_spin);
      int icoeff=0;
      for( int i=0; i<=order_spin; i++ ){
	for( int j=0; j<=order_spin-i; j++ ){
	  for( int k=0; k<=order_spin-i-j; k++ ){
	    for( int l=0; l<=order_spin-i-j-k; l++ ){
	      for( int m=0; m<=order_spin-i-j-k-l; m++ ){
		term[icoeff] = pow(T->xptar,m)*pow(T->yptar,l)*pow(T->xtar,k)*pow(T->ytar,j)*pow(1.0/T->p,i);
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
		term[icoeff] = pow(T->xptar,m)*pow(T->yptar,l)*pow(T->xtar,k)*pow(T->ytar,j)*pow(1.0/T->p,i);
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
		term[icoeff] = pow(T->xptar,m)*pow(T->yptar,l)*pow(T->xtar,k)*pow(T->ytar,j)*pow(1.0/T->p,i);
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

  TDecompSVD A_xptar(Moptics);
  TDecompSVD A_pinv(Moptics);

  for( int icoeff=0; icoeff<ncoeff_optics; icoeff++){
    if( xtar_expon_optics[icoeff] > 0 ){
      Moptics(icoeff,icoeff) = 1.0;
      for( int jcoeff=0; jcoeff<ncoeff_optics; jcoeff++ ){
	if( jcoeff != icoeff ) Moptics(icoeff,jcoeff) = 0.0;
      }
    }
  }

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
  double xfp,yfp,xpfp,ypfp;
  double xptar,yptar,ytar,xtar,p;
  double xptarrecon,yptarrecon,ytarrecon,xtarrecon,precon;

  Tout->Branch("xfp",&xfp);
  Tout->Branch("yfp",&yfp);
  Tout->Branch("xpfp",&xpfp);
  Tout->Branch("ypfp",&ypfp);
  Tout->Branch("xptar",&xptar);
  Tout->Branch("yptar",&yptar);
  Tout->Branch("ytar",&ytar);
  Tout->Branch("xtar",&xtar);
  Tout->Branch("p",&p);
  Tout->Branch("xptarrecon",&xptarrecon);
  Tout->Branch("yptarrecon",&yptarrecon);
  Tout->Branch("ytarrecon",&ytarrecon);
  Tout->Branch("xtarrecon",&xtarrecon);
  Tout->Branch("precon",&precon);

  double Pxtg, Pytg, Pztg;
  double Pxfp, Pyfp, Pzfp;
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
  double det;
  Tout->Branch("Determinant",&det);
  
  nevent=0;
  while( T->GetEntry(nevent++) ){
    if( nevent%100 == 0 ) cout << nevent << endl;
    xfp = T->xfp;
    yfp = T->yfp;
    xpfp = T->xpfp;
    ypfp = T->ypfp;
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

    //compute the spin transport matrix in dipole approximation:
    double chi_transport = T->chi + atan(T->xptar) - atan(T->xpfp);
    // double Sxx_dipole = cos( chi_transport );
    // double Sxy_dipole = 0.0;
    // double Sxz_dipole = -sin( chi_transport );
    // double Syx_dipole = 0.0;
    // double Syy_dipole = 1.0;
    // double Syz_dipole =
    TRotation Rdipole;
    Rdipole.RotateY( -chi_transport );
    
    double pinv_sum = 0.0;
    double xptar_sum = 0.0;
    double yptar_sum = 0.0;
    double ytar_sum = 0.0;
    int icoeff=0;
    for( int i=0; i<=order_optics; i++ ){
      for( int j=0; j<=order_optics-i; j++ ){
	for( int k=0; k<=order_optics-i-j; k++ ){
	  for( int l=0; l<=order_optics-i-j-k; l++ ){
	    for( int m=0; m<=order_optics-i-j-k-l; m++ ){
	      double term = pow(T->xfp,m)*pow(T->yfp,l)*pow(T->xpfp,k)*pow(T->ypfp,j)*pow(T->xtar,i);
	      xptar_sum += term * b_xptar[icoeff];
	      yptar_sum += term * b_yptar[icoeff];
	      ytar_sum += term * b_ytar[icoeff];
	      pinv_sum += term * b_pinv[icoeff];
	      icoeff++;
	    }
	  }
	}
      }
    }
    xptarrecon = xptar_sum;
    yptarrecon = yptar_sum;
    ytarrecon = ytar_sum;
    precon = 1.0/pinv_sum;

    double Sxx_sum = Rdipole.XX(), Sxy_sum = Rdipole.XY(), Sxz_sum = Rdipole.XZ();
    double Syx_sum = Rdipole.YX(), Syy_sum = Rdipole.YY(), Syz_sum = Rdipole.YZ();
    double Szx_sum = Rdipole.ZX(), Szy_sum = Rdipole.ZY(), Szz_sum = Rdipole.ZZ();

    icoeff=0;
    for( int i=0; i<=order_spin; i++ ){
      for( int j=0; j<=order_spin-i; j++ ){
	for( int k=0; k<=order_spin-i-j; k++ ){
	  for( int l=0; l<=order_spin-i-j-k; l++ ){
	    for( int m=0; m<=order_spin-i-j-k-l; m++ ){
	      double term = pow(T->xptar,m)*pow(T->yptar,l)*pow(T->xtar,k)*pow(T->ytar,j)*pow(1.0/T->p,i);
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
    

    //cout << "Determinant = " << Rot.Determinant() << endl;

    det = Rot.Determinant();
    
    Tout->Fill();
  }
  fout->Write();
  
  
}
