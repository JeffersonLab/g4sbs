#ifndef G4SBSSIMCPi0Output_h
#define G4SBSSIMCPi0Output_h 1

#include <vector>
using namespace std;

class G4SBSSIMCPi0Output {
public:
  G4SBSSIMCPi0Output();
  ~G4SBSSIMCPi0Output();
  void Clear();
  void ConvertToTreeUnits();
  
  int fnucl; // Final-state nucleon type: 1 = proton, 0 = neutron

  double sigma, Weight;
  double Q2; //negative of virtual photon invariant mass:
  double xbj; //Usual Bjorken x variable;
  double nu; //electron energy loss in lab frame
  double W; //photon-nucleon invariant mass.
  double epsilon; //virtual photon longitudinal polarization
  
  double Ebeam;
  double p_g1, theta_g1, phi_g1, px_g1, py_g1, pz_g1;
  double p_g2, theta_g2, phi_g2, px_g2, py_g2, pz_g2;
  double vx, vy, vz;
  double veE, vetheta; //scattered e- kinematics at vertex
};

#endif
