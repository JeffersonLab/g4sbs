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
// TDIS Acqu MC
#include "G4SBSAcquMCOutput.hh"
#include "AcquMCTree.h"

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

  // TDIS addition
  G4ThreeVector PiMake();
 
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
  //void SetTargRadLen
  
  void SetRasterX(double v){fRasterX = v;}
  void SetRasterY(double v){fRasterY = v;}
  
  void SetThMin(double v){fThMin = v;}
  void SetThMax(double v){fThMax = v;}
  // void SetQ2min(double v){fQ2min = v;}
  // void SetQ2max(double v){fQ2max = v;}
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
  void SetCosmicsPointer( G4ThreeVector point ){fCosmPointer = point;}
  void SetCosmicsPointerRadius( G4double radius );
  void UpdateCosmicsCeilingRadius();
  void SetCosmicsMaxAngle( G4double maxangle ){fCosmicsMaxAngle = maxangle;};

  //Initialize constant quantities so we aren't doing these calculations every event:
  //void SetConstantsInitialized( G4bool b ){ fConstantsInitialized = b; }
  //G4bool ConstantsAreInitialized(){ return fConstantsInitialized; } 
  void InitializeConstants();
  
  void SetHadronType( Hadron_t h ){fHadronType = h; }

  void SetHCALDist(double v){ fHCALdist = v;}
  
  double GetHcalDist(){ return fHCALdist; }
  double GetToFres(){ return fToFres; }

  void SetPythiaEvent( G4SBSPythiaOutput ev ){ fPythiaEvent = ev; }
  G4SBSPythiaOutput GetPythiaEvent(){ return fPythiaEvent; }

  Pythia6_tree *GetPythiaTree(){ return fPythiaTree; }
  TChain *GetPythiaChain(){ return fPythiaChain; }
  
  void LoadPythiaChain(G4String fname);
  void SetExclPythiaXSOption(G4int XSOption){fExclPyXSoption = XSOption;}

  //TDIS AcquMC
  void SetAcquMCEvent( G4SBSAcquMCOutput ev ){ fAcquMCEvent = ev; }
  G4SBSAcquMCOutput GetAcquMCEvent(){ return fAcquMCEvent; }

  AcquMCTree *GetAcquMCTree(){ return fAcquMCTree; }
  TChain *GetAcquMCChain(){ return fAcquMCChain; }
  
  void LoadAcquMCChain(G4String fname);

  void Initialize();

  G4bool GetRejectionSamplingFlag(){ return fRejectionSamplingFlag; }
  //G4bool GetRejectionSamplingInitialized(){ return fRejectionSamplingInitialized; }
  G4bool GetInitialized(){ return fInitialized; }
  void SetInitialized( G4bool b ){ fInitialized = b; }
  void SetRejectionSamplingFlag( G4bool b ){ fRejectionSamplingFlag = b; }
  void SetMaxWeight( G4double w ){ fMaxWeight = w; }
  void SetNeventsWeightCheck( G4int n ){ fNeventsWeightCheck = n; } //Number of "pre-events" used to initialize rejection sampling
  //void SetRejectionSamplingInitialized( G4bool b ){ fRejectionSamplingInitialized = b; }

  double GetGenVol(){ return fGenVol; }
  double GetLumi(){ return fLumi; }
  double GetMaxWeight(){ return fMaxWeight; }

  void InitializePythia6_Tree();

  //TDIS AcquMC
  void InitializeAcquMC_Tree();

  int GetNfoils() const { return fNfoils; }
  // std::vector<double> GetZfoil() const { return fZfoil; }
  // std::vector<double> GetThickFoil() const { return fThickFoil; }

  void SetNfoils( G4int nfoil );
  void SetFoilZandThick( const std::vector<G4double> foilz, const std::vector<G4double> foilthick );

  G4double GetPionPhoto_tmin() const { return fPionPhoto_tmin; }
  G4double GetPionPhoto_tmax() const { return fPionPhoto_tmax; }
  
  void SetPionPhoto_tmin( G4double tmin ){ fPionPhoto_tmin = tmin; }
  void SetPionPhoto_tmax( G4double tmax ){ fPionPhoto_tmax = tmax; }
  void SetUseRadiator( G4bool b ){fUseRadiator = b; }
  void SetRadthickX0( G4double thick ){ fRadiatorThick_X0 = thick; }
  
private:

  void InitializeRejectionSampling(); //Make private so it can only be called by G4SBSEventGen::Initialize()

  // double fElectronE, fNucleonE, fHadronE, fBeamE;
  // TDIS addition
  double fElectronE, fNucleonE, fHadronE, fBeamE, fNeutronE, fProton1E, fProton2E;
  // G4ThreeVector fElectronP, fNucleonP, fBeamP, fVert;
  // TDIS addition
  G4ThreeVector fElectronP, fNucleonP, fBeamP, fVert, fNeutronP, fProton1P, fProton2P;
  G4ThreeVector fHadronP;
  G4ThreeVector fBeamPol;
  
  //Define parameters for cosmics generator
  G4ThreeVector fCosmPointer;
  G4double fPointerZoneRadiusMax;
  G4double fCosmicsMaxAngle;
  G4double fCosmicsCeilingRadius;
  
  double fWeight, fQ2, fW2, fxbj, fSigma, fAperp, fApar;
  double fPt, fPl;  // born-approx polarization componenets
  int fhel;         // electron beam helicity

  //Adding these for pion photoproduction, but they can also be defined and calculated for the
  //elastic generators in principle:
  double fs, ft, fu, fcosthetaCM, fEgamma_lab;
  
  //Define additional kinematic quantities for SIDIS:
  double fz, fPh_perp, fphi_h, fphi_S, fMx;

  //Define additional kinematic quantities for TDIS
  double fxpi, ftpi, fxd, fnu, fya, fy, ff2p, ff2pi, fxa, fPtTDIS, fypi, fSigmaDIS, fSigmaTDIS;
  
  double fBeamCur;
  double fRunTime;
  long    fNevt;   //number of primary events to be generated
  //long    fNtries; //number of "tries" to generate an event (to keep track of efficiency of MC generation).
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
  double fTargLen, fRasterX, fRasterY, fTargDen; //Targ density is given in atoms or molecules/unit volume
  //double fTargRadLen; //Radiation length of target material, regardless of thickness
  double fTargRadLen; //Radiation length of target material
  double fTargUpstreamWindowRadLen;
  double fTargZatomic; //atomic number of target for purposes of any bremsstrahlung calculations:
  // set<G4String> G4TargetMaterialNames; 
  
  // G4ThreeVector fTargOffset;
  // G4ThreeVector fBeamOffset;
  // G4ThreeVector fBeamDirection;
  //For multi-foil optics target event generation, I guess we need to add (redundant) target foil thickness information here:
  G4int fNfoils;
  std::vector<std::pair<G4double, G4double> > fFoilZandThick;
  G4double fTotalThickFoil;
  std::vector<G4double> fFoilZfraction;
  
  double fPmisspar, fPmissperp, fPmissparSm;
  double fHCALdist, fToFres;
  //TDIS
  double fphi;


  double fGenVol; //Phase space generation volume
  double fLumi;   //Luminosity

  // Use radiator for photoproduction? 
  G4bool fUseRadiator;
  G4double fRadiatorThick_X0;
  
  G4double fPionPhoto_tmin, fPionPhoto_tmax; //Convert polar angle generation limits to generation limits in -t for pion photoproduction

  //void InitializePionPhotoLimits(Nucl_t);
  
  //New parameters for Bremsstrahlung generation and photoproduction:
  bool GeneratePionPhotoproduction( Nucl_t, G4LorentzVector, G4LorentzVector ); //exclusive pion photoproduction:
  //G4bool fConstantsInitialized;
  
  G4LorentzVector GetInitialNucl( Targ_t, Nucl_t );
  
  bool GenerateElastic( Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateInelastic( Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateDIS( Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateFlat( Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateBeam( Nucl_t, G4LorentzVector, G4LorentzVector );
  
  bool GenerateSIDIS( Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateTDIS( Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateWiser( Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateGun(); //The "GenerateGun" routine generates generic particles of any type, flat in costheta, phi and p within user-specified limits.
  bool GeneratePythia(); //Generates primaries from a ROOT Tree containing PYTHIA6 events.
  // TDIS AcquMC
  bool GenerateAcquMC(); //Generates primaries from a ROOT Tree containing AcquMC events.
  bool GenerateCosmics(); //Generates muons from the top of the world geometry, directed towards a point in space

  G4bool fRejectionSamplingFlag; //Flag to turn on rejection sampling;
  G4double fMaxWeight; //Maximum event weight within generation limits
  G4int fNeventsWeightCheck; //Number of "pre-events" to generate in order to check weights
  //G4bool fRejectionSamplingInitialized; //Flag to indicate whether rejection sampling has been initialized.

  G4bool fInitialized; //consolidate initialization of constant event generator parameters:
  
  double deutpdist( double );
  double he3pdist( Nucl_t, double );

  double f2p (double);
  double f2pi (double, double, double);
  double c12pdist( double );
  
  DSS2007FF fFragFunc; //Class to calculate fragmentation functions using DSS2007

  //void LoadTargetData(); //why is this here? Not implemented...

  //TFile *fPythiaFile;
  //TTree *fPythiaTree;
  long fchainentry;
  TChain *fPythiaChain;
  Pythia6_tree *fPythiaTree;

  map<G4String, G4double> fPythiaSigma;
  
  G4SBSPythiaOutput fPythiaEvent;

  G4int fExclPyXSoption; //Flag to choose "Exclusive pythia" event-by-event cross section events

  // TDIS AcquMC
  long fAcquMCchainentry;
  TChain *fAcquMCChain;
  AcquMCTree *fAcquMCTree;
  G4SBSAcquMCOutput fAcquMCEvent;

  G4double TriangleFunc(G4double a, G4double b, G4double c );
};

#endif//G4SBSEVENTGEN_HH
