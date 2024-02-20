#include "G4SBSSIMCOutput.hh"
#include "G4SystemOfUnits.hh"

G4SBSSIMCOutput::G4SBSSIMCOutput(){
  sigma = 1.0; //default to 1.0 mb
  Clear();
}

G4SBSSIMCOutput::~G4SBSSIMCOutput(){
  ;
}

void G4SBSSIMCOutput::Clear(){
  fnucl = 0;

  sigma = 0.0;
  Weight = 0.0;
  Q2 = 0.0;
  xbj = 0.0;
  nu = 0.0;
  W = 0.0;
  epsilon = 0.0;
  
  Ebeam = 0.0;
  p_e = 0.0;
  theta_e = 0.0;
  phi_e = 0.0;
  px_e = 0.0;
  py_e = 0.0;
  pz_e = 0.0;
  p_n = 0.0;
  theta_n = 0.0;
  phi_n = 0.0;
  px_n = 0.0;
  py_n = 0.0;
  pz_n = 0.0;
  
  vx = 0.0;
  vy = 0.0;
  vz = 0.0;
}

void G4SBSSIMCOutput::ConvertToTreeUnits(){ //This is called once per event after primary vertices are generated. When this is called, all quantities should be in GEANT4 standard units of MeV, ns, cm;
  //sigma /= (1.0/(cm*cm));
  //Q2 /= (GeV*GeV);
  //W /= (GeV);
  //nu /= (GeV);
  
  Ebeam /= GeV;
  
  p_e /= GeV;
  px_e /= GeV;
  py_e /= GeV;
  pz_e /= GeV;
  p_n /= GeV;
  px_n /= GeV;
  py_n /= GeV;
  pz_n /= GeV;

  vx /= m;
  vy /= m;
  vz /= m;
}
