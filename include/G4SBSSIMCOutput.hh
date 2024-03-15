#ifndef G4SBSSIMCOutput_h
#define G4SBSSIMCOutput_h 1

#include <vector>
using namespace std;

class G4SBSSIMCOutput {
public:
  G4SBSSIMCOutput();
  ~G4SBSSIMCOutput();
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
  double p_e, theta_e, phi_e, px_e, py_e, pz_e;
  double p_n, theta_n, phi_n, px_n, py_n, pz_n;
  double vx, vy, vz;
  double veE, vetheta; //scattered e- kinematics at vertex
};

#endif
