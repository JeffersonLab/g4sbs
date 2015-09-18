#include "G4SBSPythiaOutput.hh"

G4SBSPythiaOutput::G4SBSPythiaOutput(){
  Clear();
}

G4SBSPythiaOutput::~G4SBSPythiaOutput(){
  ;
}

void G4SBSPythiaOutput::Clear(){
  Ebeam = 0.0;
  Eprime = 0.0;
  theta_e = 0.0;
  phi_e = 0.0;
  px_e = 0.0;
  py_e = 0.0;
  pz_e = 0.0;
  vx_e = 0.0;
  vy_e = 0.0;
  vz_e = 0.0;
  Egamma = 0.0;
  theta_gamma = 0.0;
  phi_gamma = 0.0;
  px_gamma = 0.0;
  py_gamma = 0.0;
  pz_gamma = 0.0;
  vx_gamma = 0.0;
  vy_gamma = 0.0;
  vz_gamma = 0.0;
  Q2 = 0.0;
  xbj = 0.0;
  y = 0.0;
  W2 = 0.0;

  Nprimaries = 0;
  PID.clear();
  genflag.clear();
  Px.clear();
  Py.clear();
  Pz.clear();
  M.clear();
  E.clear();
  P.clear();
  t.clear();
  vx.clear();
  vy.clear();
  vz.clear();
  theta.clear();
  phi.clear();
}

