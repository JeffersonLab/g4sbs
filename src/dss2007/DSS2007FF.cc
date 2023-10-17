#include "DSS2007FF.hh"
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>

using namespace std;

DSS2007FF::DSS2007FF(){
  //Initialize Q2 and x values of the grid:
  double Q2[NQ] = { 1.0, 1.25, 1.5, 2.5, 4.0, 
		    6.4, 10.0, 15.0, 25.0, 40.0, 
		    64.0, 100.0, 180.0, 320.0, 580.0, 
		    1000.0, 1800.0, 3200.0, 5800.0, 
		    10000.0, 18000.0, 32000.0, 58000.0, 100000.0 };
  double x[Nx] = { 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.095, 0.1, 0.125,
		   0.15, 0.175, 0.2, 0.225, 0.25, 0.275, 0.3, 0.325, 0.35, 0.375, 0.4, 0.45, 0.5, 0.55, 
		   0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.93, 1.0 };

  double exp1[Nparton] = {4,4,4,7,7,4,4,4,4};
  double exp2[Nparton] = {0.5, 0.5, 0.5, 0.3, 0.3, 0.3, 0.5, 0.5, 0.5 };

  for(int i=0; i<Nparton; i++){
    exponent[0][i] = exp1[i]; //power of 1-x in xD
    exponent[1][i] = exp2[i]; //power of x in xD
  }

  for( int i=0; i<NQ; i++){ Qtable[i] = Q2[i]; }
  for( int i=0; i<Nx; i++){ xtable[i] = x[i]; }

  LoadInterpolationGrids();
}

void DSS2007FF::SetGridPath(string fn){ //Re-load interpolation grids when this method is invoked.
  fgridpath = fn;
  LoadInterpolationGrids();
}

void DSS2007FF::LoadInterpolationGrids(){ //Only leading order for now:

  string fnamepi = fgridpath;
  fnamepi += "/PILO.GRID";
  string fnamek = fgridpath;
  fnamek += "/KALO.GRID";
  string fnamepro = fgridpath;
  fnamepro += "/PROLO.GRID";
  
  ifstream pilo(fnamepi.c_str());
  ifstream kalo(fnamek.c_str());
  ifstream prolo(fnamepro.c_str());
  
  for(int ix=0; ix<Nx-1; ix++){
    for(int iq=0; iq<NQ; iq++){
      for(int ip=0; ip<9; ip++){
	pilo >> Parton[0][ip][iq][ix];
	kalo >> Parton[1][ip][iq][ix];
	prolo >> Parton[2][ip][iq][ix];
	// cout << "ix,iq,ip,pionq=" << ix << ", " << iq << ", " << ip << ", " << Parton[0][ip][iq][ix] << endl;
	// cout << "ix,iq,ip,kaonq=" << ix << ", " << iq << ", " << ip << ", " << Parton[0][ip][iq][ix] << endl;
	// cout << "ix,iq,ip,protonq=" << ix << ", " << iq << ", " << ip << ", " << Parton[0][ip][iq][ix] << endl;
      }
    }
  }

  //Initialize arrays for the interpolation grids:
  for(int ih=0; ih<3; ih++){
    for(int ip=0; ip<Nparton; ip++){
      for(int iq=0; iq<NQ; iq++){
	for(int ix=0; ix<Nx-1; ix++){
	  double x = xtable[ix]; 
	  double F = pow( 1.0 - x, exponent[0][ip] ) * pow( x, exponent[1][ip] );
	  xD[ih][ip][ix][iq] = Parton[ih][ip][iq][ix] / F;
	}
	xD[ih][ip][Nx-1][iq] = 0.0;
      }
    }
  }
}

void DSS2007FF::Interpolate(int ihadron, double x, double Q2, vector<double> &xDout ){
  
  if( ihadron < 0 ){ 
    ihadron = 0;
    cout << "ihadron < 0 not supported, using pi" << endl;
  } else if( ihadron > 2 ){
    ihadron = 2;
    cout << "ihadron > 2 not supported, using proton" << endl;
  }

  xDout.clear();
  xDout.resize(Nparton);

  //Force x and Q2 inside the range of the grid:
  double logx = log(min(max(x,xtable[0]),xtable[Nx-1]) );
  double logQ2 = log(min(max(Q2,Qtable[0]),Qtable[NQ-1]) );
 
  //Perform two-dimensional linear interpolation of the table in log(x) and log(Q2)
  double logxtable[Nx];
  double logQ2table[NQ];
  for(int i=0; i<Nx; i++){ logxtable[i] = log(xtable[i]); }
  for(int i=0; i<NQ; i++){ logQ2table[i] = log(Qtable[i]); }

  //Initial guesses assume fixed bin width:
  int ix = int( (logx - logxtable[0])*double(Nx-1)/(logxtable[Nx-1]-logxtable[0]) );
  int iq = int( (logQ2 - logQ2table[0])*double(NQ-1)/(logQ2table[NQ-1]-logQ2table[0]) );

  //Increment or decrement bin numbers until we are in the correct bin:
  while( logx > logxtable[ix+1] ){ ix++; } //log(x) is above the bin upper edge, increment bin by 1
  while( logx < logxtable[ix] ){ ix--; }   //log(x) is below the bin lower edge, decrement bin by 1
  while( logQ2 > logQ2table[iq+1] ){ iq++; } 
  while( logQ2 < logQ2table[iq] ){ iq--; }

  //Get the corners of the grid point:
  
  double f11, f12, f21, f22;
 
  for(int iparton=0; iparton<Nparton; iparton++){
    f11 = xD[ihadron][iparton][ix][iq]; //low-x, low-Q
    f12 = xD[ihadron][iparton][ix][iq+1]; //low-x, high-Q
    f21 = xD[ihadron][iparton][ix+1][iq]; //high-x, low-Q
    f22 = xD[ihadron][iparton][ix+1][iq+1]; //high-x, high-Q

    double xfrac = (logx - logxtable[ix])/( logxtable[ix+1]-logxtable[ix] );
    double Qfrac = (logQ2 - logQ2table[iq])/( logQ2table[iq+1]-logQ2table[iq] );

    double fxint_Qlow = xfrac * f21 + (1.0-xfrac) * f11; //Interpolate in x at low Q2 point
    double fxint_Qhigh = xfrac * f22 + (1.0-xfrac) * f12; //Interpolate in x at high Q2 point
    double fxyint      = Qfrac * fxint_Qhigh + (1.0-Qfrac) * fxint_Qlow; //Interpolate in Q2 between low and high interpolated x points:

    // cout << "log(z), log(Q2), zfrac, Q2frac, f11, f12, f21, f22, f = " << logx << ", " 
    // 	 << logQ2 << ", " << xfrac << ", " << Qfrac << ", " << f11 << ", " << f12 << ", " << f21 << ", " << f22 
    // 	 << ", " << fxyint << endl;

    xDout[iparton] = fxyint;
  }
}

void DSS2007FF::GetFFs( int ihadron, int icharge, double z, double Q2, vector<double> &Dqh ){

  vector<double> xDtemp;

  Interpolate( ihadron, z, Q2, xDtemp );

  //Order for parton grids is:
  // 0: utot
  // 1: dtot
  // 2: stot
  // 3: ctot
  // 4: btot
  // 5: gluon
  // 6: uvalence
  // 7: dvalence
  // 8: svalence

  for(int i=0; i<Nparton; i++){
    xDtemp[i] *= pow(1.0-z,exponent[0][i])*pow(z, exponent[1][i] );
    //cout << "iparton, xD(iparton) = " << i << ", " << xDtemp[i] << endl;
  }

  double Up = (xDtemp[0]+xDtemp[6])/2.0;
  double UBp = (xDtemp[0]-xDtemp[6])/2.0;
  double Dp = (xDtemp[1]+xDtemp[7])/2.0;
  double DBp = (xDtemp[1]-xDtemp[7])/2.0;
  double Sp = (xDtemp[2]+xDtemp[8])/2.0;
  double SBp = (xDtemp[2]-xDtemp[8])/2.0;

  double Cp = xDtemp[3]/2.0;
  double Bp = xDtemp[4]/2.0;

  Dqh.resize(Nparton);
  
  if( icharge == 1 ){ //+charge
    Dqh[0] = Up; //Duh
    Dqh[1] = UBp; //Dubarh
    Dqh[2] = Dp; //Ddh
    Dqh[3] = DBp; //Ddbarh
    Dqh[4] = Sp;
    Dqh[5] = SBp;
  } else if( icharge == -1 ){ //-charge: exchange quarks and antiquarks
    Dqh[0] = UBp; //Duh
    Dqh[1] = Up; //Dubarh
    Dqh[2] = DBp; //Ddh
    Dqh[3] = Dp; //Ddbarh
    Dqh[4] = SBp;
    Dqh[5] = Sp;
  } else { //neutral: take average of plus and minus (isospin symmetry assumption):
    Dqh[0] = (Up+UBp)/2.0; //Duh
    Dqh[1] = (Up+UBp)/2.0; //Dubarh
    Dqh[2] = (Dp+DBp)/2.0; //Ddh
    Dqh[3] = (Dp+DBp)/2.0; //Ddbarh
    Dqh[4] = (Sp+SBp)/2.0; 
    Dqh[5] = (Sp+SBp)/2.0;
  }
  Dqh[6] = Cp;
  Dqh[7] = Bp;
  Dqh[8] = xDtemp[5]/2.0;

  for(int i=0; i<Nparton; i++){
    Dqh[i] /= z;
  }
  
}
