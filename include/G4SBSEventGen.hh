#ifndef G4SBSEVENTGEN_HH
#define G4SBSEVENTGEN_HH

#include "globals.hh"
#include "sbstypes.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandGauss.h"
#include "G4SBSIO.hh"
#include "DSS2007FF.hh"
#include "G4SBSPythiaOutput.hh"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "Pythia6_tree.h"
#include "G4SBSESEPP.hh"

#define MAXMOMPT 1000 // N points for targets momentum distribution interpolations

class G4SBSEventGen {
public:
  G4SBSEventGen();
  ~G4SBSEventGen();
  
  double GetBeamE(){ return fBeamE; }
  G4ThreeVector GetBeamP(){ return fBeamP; }
  G4ThreeVector GetBeamPol(){ return fBeamPol; }
  
  G4ThreeVector GetV(){ return fVert; }
  
  double GetElectronE(){ return fElectronE; }
  double GetNucleonE(){ return fNucleonE; }
  double GetHadronE(){ return fHadronE; }
 
  G4ThreeVector GetElectronP(){ return fElectronP; }
  G4ThreeVector GetNucleonP(){ return fNucleonP; }
  G4ThreeVector GetHadronP(){ return fHadronP; }
 
  Nucl_t GetNucleonType(){ return fNuclType; }
  Nucl_t GetFinalNucleon(){ return fFinalNucl; }

  Hadron_t GetHadronType(){ return fHadronType; }

  double GetPt(){ return fPt; }
  double GetPl(){ return fPl; }
  
  bool GenerateEvent();
  
  ev_t GetEventData();
  
  void SetNevents(int n){fNevt = n;}
  void SetBeamCur(double c){fBeamCur = c;}
  void SetBeamE(double c){fBeamE= c; fBeamP = G4ThreeVector(0.0, 0.0, c); }
  void SetRunTime(double t){fRunTime = t;}
  
  void SetKine(Kine_t t ){fKineType = t;}
  Kine_t GetKine(){return fKineType;}
  
  void SetTarget(Targ_t t ){fTargType = t;}
  void SetTargLen(double len){fTargLen = len;}
  void SetTargDen(double den){fTargDen = den;}
  
  void SetRasterX(double v){fRasterX = v;}
  void SetRasterY(double v){fRasterY = v;}
  
  void SetThMin(double v){fThMin = v;}
  void SetThMax(double v){fThMax = v;}
  void SetPhMin(double v){fPhMin = v;}
  void SetPhMax(double v){fPhMax = v;}
  
  void SetEeMin(double v){fEeMin = v;}
  void SetEeMax(double v){fEeMax = v;}

  void SetEhadMin(double v){fEhadMin = v;}
  void SetEhadMax(double v){fEhadMax = v;}
  void SetThMin_had(double v){fThMin_had = v; }
  void SetThMax_had(double v){fThMax_had = v; }
  void SetPhMin_had(double v){fPhMin_had = v; }
  void SetPhMax_had(double v){fPhMax_had = v; }

  void SetHadronType( Hadron_t h ){fHadronType = h; }

  void SetHCALDist(double v){ fHCALdist = v;}
  
  double GetHcalDist(){ return fHCALdist; }
  double GetToFres(){ return fToFres; }

  void SetPythiaEvent( G4SBSPythiaOutput ev ){ fPythiaEvent = ev; }
  G4SBSPythiaOutput GetPythiaEvent(){ return fPythiaEvent; }

  Pythia6_tree *GetPythiaTree(){ return fPythiaTree; }
  TChain *GetPythiaChain(){ return fPythiaChain; }
  
  void LoadPythiaChain(G4String fname);
private:
  double fElectronE, fNucleonE, fHadronE, fBeamE;
  G4ThreeVector fElectronP, fNucleonP, fBeamP, fVert;
  G4ThreeVector fHadronP;
  G4ThreeVector fBeamPol;

  double fWeight, fQ2, fW2, fxbj, fSigma, fAperp, fApar;
  double fPt, fPl;  // born-approx polarization componenets
  int fhel;         // electron beam helicity
  
  //Define additional kinematic quantities for SIDIS:
  double fz, fPh_perp, fphi_h, fphi_S, fMx;

  double fBeamCur;
  double fRunTime;
  int    fNevt;
  double Wfact;
  
  Nucl_t fNuclType, fFinalNucl;
  Targ_t fTargType;
  Kine_t fKineType;
  
  // Which hadron species are we considering for pi/K SIDIS?
  Hadron_t fHadronType; //Currently available: pi+/-/0, K+/-, p/pbar

  double fThMin, fThMax, fPhMin, fPhMax; //Angular generation limits for electron arm
  double fEeMin, fEeMax; //Electron energy generation limits
  double fThMin_had, fThMax_had, fPhMin_had, fPhMax_had; //Angular generation limits for hadron arm 
  double fEhadMin, fEhadMax; //Hadron (total) energy generation limits (for SIDIS case)--Later we will want to add exclusive hadron production.
  double fTargLen, fRasterX, fRasterY, fTargDen;
  double fPmisspar, fPmissperp, fPmissparSm;
  double fHCALdist, fToFres;

  G4LorentzVector GetInitialNucl( Targ_t, Nucl_t );
  
  bool GenerateElastic( Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateInelastic( Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateDIS( Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateFlat( Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateBeam( Nucl_t, G4LorentzVector, G4LorentzVector );
  
  bool GenerateSIDIS( Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateWiser( Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateGun(); //The "GenerateGun" routine generates generic particles of any type, flat in costheta, phi and p within user-specified limits.
  bool GeneratePythia(); //Generates primaries from a ROOT Tree containing PYTHIA6 events.
  bool GenerateESEPP();

  double deutpdist( double );
  double he3pdist( Nucl_t, double );
  
  DSS2007FF fFragFunc; //Class to calculate fragmentation functions using DSS2007

  //void LoadTargetData(); //why is this here? Not implemented...

  //TFile *fPythiaFile;
  //TTree *fPythiaTree;
  long fchainentry;
  TChain *fPythiaChain;
  Pythia6_tree *fPythiaTree;
  
  G4SBSPythiaOutput fPythiaEvent;
  G4SBSESEPP *fESEPP;
};

#endif//G4SBSEVENTGEN_HH
