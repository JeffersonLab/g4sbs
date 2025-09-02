#include "G4SBSSIMCPi0Output.hh"
#include "G4SystemOfUnits.hh"

G4SBSSIMCPi0Output::G4SBSSIMCPi0Output(){
  sigma = 1.0; //default to 1.0 mb
  Clear();
}

G4SBSSIMCPi0Output::~G4SBSSIMCPi0Output(){
  ;
}

void G4SBSSIMCPi0Output::Clear(){
  fnucl = 0;

  sigma = 0.0;
  Weight = 0.0;
  Q2 = 0.0;
  xbj = 0.0;
  nu = 0.0;
  W = 0.0;
  epsilon = 0.0;
  
  Ebeam = 0.0;
  p_g1 = 0.0;
  theta_g1 = 0.0;
  phi_g1 = 0.0;
  px_g1 = 0.0;
  py_g1 = 0.0;
  pz_g1 = 0.0;
  p_g2 = 0.0;
  theta_g2 = 0.0;
  phi_g2 = 0.0;
  px_g2 = 0.0;
  py_g2 = 0.0;
  pz_g2 = 0.0;
  
  vx = 0.0;
  vy = 0.0;
  vz = 0.0;
}

void G4SBSSIMCPi0Output::ConvertToTreeUnits(){ //This is called once per event after primary vertices are generated. When this is called, all quantities should be in GEANT4 standard units of MeV, ns, cm;
  //sigma /= (1.0/(cm*cm));
  //Q2 /= (GeV*GeV);
  //W /= (GeV);
  //nu /= (GeV);
  
  Ebeam /= GeV;
  
  p_g1 /= GeV;
  px_g1 /= GeV;
  py_g1 /= GeV;
  pz_g1 /= GeV;
  p_g2 /= GeV;
  px_g2 /= GeV;
  py_g2 /= GeV;
  pz_g2 /= GeV;

  vx /= m;
  vy /= m;
  vz /= m;
}
