#include "G4SBSPythiaOutput.hh"
#include "G4SystemOfUnits.hh"

G4SBSPythiaOutput::G4SBSPythiaOutput(){
  Sigma = 1.0; //default to 1.0 mb
  SigmaDiff = 0.0; //default to 0.0
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
  Delta2 = 0.0;
  phigg = 0.0;

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

void G4SBSPythiaOutput::ConvertToTreeUnits(){ //This is called once per event after primary vertices are generated. When this is called, all quantities should be in GEANT4 standard units of MeV, ns, cm;
  Ebeam /= GeV;
  Eprime /= GeV;
  px_e /= GeV;
  py_e /= GeV;
  pz_e /= GeV;
  vx_e /= m;
  vy_e /= m;
  vz_e /= m;
  Egamma /= GeV;
  px_gamma /= GeV;
  py_gamma /= GeV;
  pz_gamma /= GeV;
  vx_gamma /= m;
  vy_gamma /= m;
  vz_gamma /= m;
  Q2 /= (GeV*GeV);
  W2 /= (GeV*GeV);
  Delta2 /= (GeV*GeV);

  for( int i=0; i<Nprimaries; i++ ){
    Px[i] /= GeV;
    Py[i] /= GeV;
    Pz[i] /= GeV;
    M[i] /= GeV;
    E[i] /= GeV;
    P[i] /= GeV;
    t[i] /= ns;
    vx[i] /= m;
    vy[i] /= m;
    vz[i] /= m;
  }

}
