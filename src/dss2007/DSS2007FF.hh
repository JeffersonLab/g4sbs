#ifndef DSS2007FF_hh 
#define DSS2007FF_hh 1

#include <vector>

using namespace std;

class DSS2007FF {
public:
  DSS2007FF();
  
  void Interpolate(int ihadron, double x, double Q2, vector<double> &Dout);
  void LoadInterpolationGrids();
  //void GetFFs( int ihadron, int icharge, int order, double x, double Q2, vector<double> &Dqh );
  
  void GetFFs( int ihadron, int icharge, double x, double Q2, vector<double> &Dqh );
private:
  static const int Nparton = 9;
  static const int Nx = 35;
  static const int NQ = 24; 

  double exponent[2][Nparton];
  double xD[3][Nparton][Nx][NQ];
  double Parton[3][Nparton][NQ][Nx-1]; //Three FFs for pion, kaon, proton:
  double Qtable[NQ], xtable[Nx];
  
};

#endif
