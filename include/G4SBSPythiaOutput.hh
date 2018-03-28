#ifndef G4SBSPythiaOutput_h
#define G4SBSPythiaOutput_h 1

#include <vector>
using namespace std;

class G4SBSPythiaOutput {
public:
  G4SBSPythiaOutput();
  ~G4SBSPythiaOutput();
  void Clear();
  void ConvertToTreeUnits();

  double Sigma; //PYTHIA6 calculated total cross section
  double SigmaDiff; //Exclusive cross section difference (+)-(-)
  double Ebeam; //Incident beam energy
  //The outgoing electron:
  double Eprime; //Outgoing electron energy
  double theta_e, phi_e, px_e, py_e, pz_e, vx_e, vy_e, vz_e; //
  //The virtual photon:
  double Egamma; //Virtual photon energy;
  double theta_gamma, phi_gamma, px_gamma, py_gamma, pz_gamma, vx_gamma, vy_gamma, vz_gamma;
  //The kinematic invariants:
  double Q2; //negative of virtual photon invariant mass:
  double xbj; //Usual Bjorken x variable;
  double y; //fractional electron energy loss in lab frame
  double W2; //photon-nucleon invariant mass.
  //Exclusive variable
  double Delta2; // quadrimomentum transfer (usually called t)
  double phigg; // angle between virtual photon and whatever is produced (e.g. photon for DVCS)

  int Nprimaries; //Number of final state particles:
  vector<int> PID, genflag; //PID code and flag (1 if particle was actually generated, 0 otherwise)
  vector<double> Px, Py, Pz, M, E, P, t, vx, vy, vz, theta, phi;
  
};

#endif
