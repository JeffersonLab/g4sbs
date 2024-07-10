#include "G4SBSAcquMCOutput.hh"
#include "G4SystemOfUnits.hh"

G4SBSAcquMCOutput::G4SBSAcquMCOutput(){
  Clear();
}

G4SBSAcquMCOutput::~G4SBSAcquMCOutput(){
  ;
}

void G4SBSAcquMCOutput::Clear(){
  Vx = 0.0;
  Vy = 0.0;
  Vz = 0.0;
  Px = 0.0;
  Py = 0.0;
  Pz = 0.0;
  Pt = 0.0;
  E = 0.0;
  pid = 2212; // hard coded input protons for now to check mtpc
  // mass = 938.2720813;//MeV

}

void G4SBSAcquMCOutput::ConvertToTreeUnits(){ //This is called once per event after primary vertices are generated. When this is called, all quantities should be in GEANT4 standard units of MeV, ns, cm;
  // Ebeam /= GeV;
  // Eprime /= GeV;
  // px_e /= GeV;
  // py_e /= GeV;
  // pz_e /= GeV;
  // vx_e /= m;
  // vy_e /= m;
  // vz_e /= m;
  // Egamma /= GeV;
  // px_gamma /= GeV;
  // py_gamma /= GeV;
  // pz_gamma /= GeV;
  // vx_gamma /= m;
  // vy_gamma /= m;
  // vz_gamma /= m;
  // Q2 /= (GeV*GeV);
  // W2 /= (GeV*GeV);
  // Delta2 /= (GeV*GeV);

  Vx /= m;
  Vy /= m;
  Vz /= m;
  Px /= GeV;
  Py /= GeV;
  Pz /= GeV;
  Pt /= GeV;
  E /= GeV;

  // for( int i=0; i<Nprimaries; i++ ){
  //   Px[i] /= GeV;
  //   Py[i] /= GeV;
  //   Pz[i] /= GeV;
  //   M[i] /= GeV;
  //   E[i] /= GeV;
  //   P[i] /= GeV;
  //   t[i] /= ns;
  //   vx[i] /= m;
  //   vy[i] /= m;
  //   vz[i] /= m;
  // }
  //Need to check
}
