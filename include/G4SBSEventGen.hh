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
#include "G4SBSSIMCOutput.hh"
#include "G4SBSUtil.hh"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "simc_tree.h"
#include "Pythia6_tree.h"
// TDIS Acqu MC
#include "G4SBSAcquMCOutput.hh"
#include "AcquMCTree.h"

#define MAXMOMPT 1000 // N points for targets momentum distribution interpolations

#include "G4SBSTDISGen.hh"
class G4SBSTDISGen;


class G4SBSEventGen {
public:
  G4SBSEventGen();
  ~G4SBSEventGen();
  
  double GetBeamE(){ return fBeamE; }
  G4ThreeVector GetBeamP(){ return fBeamP; }
  G4ThreeVector GetBeamPol(){ return fBeamPolMagnitude * fBeamPolDirection; }
  G4ThreeVector GetBeamPolDirection(){ return fBeamPolDirection; }
  G4double GetBeamPolMagnitude(){ return fBeamPolMagnitude; }
  
  G4ThreeVector GetV(){ return fVert; }
  
  double GetElectronE(){ return fElectronE; }
  double GetNucleonE(){ return fNucleonE; }
  double GetHadronE(){ return fHadronE; }
  double GetProtonSpecE(){ return fProtonSpecE; }
 
  G4ThreeVector GetElectronP(){ return fElectronP; }
  G4ThreeVector GetNucleonP(){ return fNucleonP; }
  G4ThreeVector GetHadronP(){ return fHadronP; }
  // // // // HEAD
  G4ThreeVector GetProtonSpecP(){ return fProtonSpecP; }

  // TDIS addition
  G4ThreeVector PiMake();
 
  // // // // 

  G4ThreeVector GetTargPol(){ return fTargPolMagnitude * fTargPolDirection; }
  G4ThreeVector GetTargPolDirection(){ return fTargPolDirection; }
  G4double GetTargPolMagnitude(){ return fTargPolMagnitude; }

  G4bool GetRandomizeTargetSpin() { return fRandomizeTargetSpin; }
  G4int GetNumTargetSpinDirections() { return fNumTargetSpinDirections; }
  vector<G4double> GetTargetThetaSpin() { return fTargetThetaSpin; }
  vector<G4double> GetTargetPhiSpin() { return fTargetPhiSpin; }
  G4double GetTargetThetaSpin(G4int ispin);
  G4double GetTargetPhiSpin(G4int ispin);
  
  // // // // 11a33984f47772444ffb08222f8a978d2bee837e
  G4SBS::Nucl_t GetNucleonType(){ return fNuclType; }
  G4SBS::Nucl_t GetFinalNucleon(){ return fFinalNucl; }

  G4SBS::Hadron_t GetHadronType(){ return fHadronType; }

  G4double GetAUT_Collins(){ return fAUT_Collins; }
  G4double GetAUT_Sivers(){ return fAUT_Sivers; }

  G4double GetAUT_Collins_min(){ return fAUT_Collins_min; }
  G4double GetAUT_Collins_max(){ return fAUT_Collins_max; }

  G4double GetAUT_Sivers_min(){ return fAUT_Sivers_min; }
  G4double GetAUT_Sivers_max(){ return fAUT_Sivers_max; }
  
  double GetPt(){ return fPt; }
  double GetPl(){ return fPl; }
  
  bool GenerateEvent();
  
  ev_t GetEventData();
  ev_tdis_t GetTDISEventData();
  
  void SetNevents(int n){fNevt = n;}
  void SetFirstEvent(long n);

  long GetFirstEvent() const { return fFirstEvent; }
  
  void SetBeamCur(double c){fBeamCur = c;}
  void SetBeamE(double c){fBeamE= c; fBeamP = G4ThreeVector(0.0, 0.0, c); }
  void SetRunTime(double t){fRunTime = t;}
  
  void SetKine(G4SBS::Kine_t t ){fKineType = t;}
  G4SBS::Kine_t GetKine(){return fKineType;}
  
  G4SBS::Targ_t GetTarget(){return fTargType;}
  
  void SetTarget(G4SBS::Targ_t t ){fTargType = t;}
  void SetTargLen(double len){fTargLen = len;}
  void SetTargDen(double den){fTargDen = den;}
  //void SetTargRadLen
  
  void SetRasterX(double v){fRasterX = v;}
  void SetRasterY(double v){fRasterY = v;}
  
  void SetRasterRadius(double v){fCircularRasterRadius = v;}
  void SetBeamSpotSize(double v){fBeamSpotSize = v;}

  // D Flay (Aug 2020) 
  void SetBeamOffsetX(double v) { fBeamOffsetX = v;} 
  void SetBeamOffsetY(double v) { fBeamOffsetY = v;}

  // D Flay (Oct 2020) 
  void SetBeamAngleX(double v)  { fBeamAngleX = v; }
  void SetBeamAngleY(double v)  { fBeamAngleY = v; }
  void SetBeamAngleZ(double v)  { fBeamAngleZ = v; }
  
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
  
  void SetHadronType( G4SBS::Hadron_t h ){fHadronType = h; }

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

  void SetSIMCEvent( G4SBSSIMCOutput ev ){ fSIMCEvent = ev; }
  G4SBSSIMCOutput GetSIMCEvent(){ return fSIMCEvent; }

  simc_tree *GetSIMCTree(){ return fSIMCTree; }
  TChain *GetSIMCChain(){ return fSIMCChain; }
  
  void LoadSIMCChain(G4String fname);


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
  void InitializeSIMC_Tree();

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
  // // // // HEAD
  // // // // 
  //void SetTargPol( G4ThreeVector pol ){ fTargPolDirection = pol.unit(); }
  void SetTargPol( G4ThreeVector pol ){ //version with single three-vector argument; takes magnitude of vector as polarization, direction:
    fTargPolDirection = pol.unit();
    fTargPolMagnitude = std::max(0.0, std::min(1.0, pol.mag() ) );
  }
  void SetTargPol( G4ThreeVector dir, G4double mag ) //version with two arguments, explicitly set direction and magnitude
  {
    fTargPolDirection = dir.unit();
    if( mag < 0.0 ){
      fTargPolDirection *= -1.0;
      mag = fabs(mag);
    }
    fTargPolMagnitude = std::max(0.0,std::min(1.0,mag));
  };
  void SetTargPolDir( G4ThreeVector dir ){ fTargPolDirection = dir.unit(); } 
  void SetTargPolMag( G4double mag ){
    fTargPolMagnitude = std::min(1.0,fabs(mag));
    if( mag < 0.0 ){
      fTargPolDirection *= -1.0;
    }
  }

  void SetBeamPol( G4ThreeVector pol ){
    fBeamPolDirection = pol.unit();
    fBeamPolMagnitude = std::max(0.0, std::min(1.0,pol.mag()));
  }
  void SetBeamPol( G4ThreeVector dir, G4double mag )
  {
    fBeamPolDirection = dir.unit();
    if( mag < 0.0 ){
      fBeamPolDirection *= -1.0;
      mag = fabs(mag);
    }
    fBeamPolMagnitude = std::max(0.0, std::min(1.0,mag));
  };
  void SetBeamPolDir( G4ThreeVector dir ){ fBeamPolDirection = dir.unit(); }
  void SetBeamPolMag( G4double mag ){
    fBeamPolMagnitude = std::min(1.0,fabs(mag));
    if( mag < 0.0 ){
      fBeamPolDirection *= -1.0;
    }
  }

  void SetTargZoffset( double z ){ fTargZoffset = z; }
  double GetTargZoffset() const { return fTargZoffset; }
  
  void SetRandomizeTargetSpin( G4bool flag ){ fRandomizeTargetSpin = flag; }
  void SetNumTargetSpinDirections( G4int nspin );
  void SetTargetThetaSpinVector( vector<double> &thspin ){ fTargetThetaSpin = thspin; }
  void SetTargetPhiSpinVector( vector<double> &phspin ){ fTargetPhiSpin = phspin; }

  //Setters for single element of array:
  void SetTargetThetaSpin( G4int ispin, G4double theta );
  void SetTargetPhiSpin( G4int ispin, G4double phi );
  //void SetTargetThetaSpin( G4int ispin, G4double value )

  void SetAUT_Collins_Sivers( G4double Acoll, G4double Asiv ){
    fAUT_Collins = Acoll;
    fAUT_Sivers = Asiv;
  }
  // // // // 11a33984f47772444ffb08222f8a978d2bee837e
  
private:

  //TDIS (I think I used these variables to cross-check I am storing the correct values)
  G4LorentzVector tElectron_f; //scattered 4-vector electron (CA)
  G4LorentzVector tNucleon_f; //final 4-vector nucleon (CA)
  
  void InitializeRejectionSampling(); //Make private so it can only be called by G4SBSEventGen::Initialize()

  // double fElectronE, fNucleonE, fHadronE, fBeamE;
  // TDIS addition
  double fElectronE, fNucleonE, fHadronE, fBeamE, fProtonSpecE, fNeutronE, fProton1E, fProton2E;
  // G4ThreeVector fElectronP, fNucleonP, fBeamP, fVert;
  // TDIS addition
  G4ThreeVector fElectronP, fNucleonP, fBeamP, fVert, fProtonSpecP, fNeutronP, fProton1P, fProton2P;
  G4ThreeVector fHadronP;
  G4ThreeVector fBeamPolDirection;
  G4ThreeVector fTargPolDirection;

  G4double fBeamPolMagnitude;
  G4double fTargPolMagnitude;

  //Define parameters for randomized target spin generation (with a discrete number of directions):
  G4bool fRandomizeTargetSpin;
  G4int  fNumTargetSpinDirections;
  // vectors of target spin angles:
  vector<G4double> fTargetThetaSpin;
  vector<G4double> fTargetPhiSpin; 
  
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

  //Define quantities to hold "true" values of Collins/Sivers asymmetries for target SSA simulations:
  double fAUT_Collins, fAUT_Sivers; //These will be the undiluted asymmetry moments for the struck nucleon
  double fAUT_Collins_min, fAUT_Collins_max; //"min" and "max" values from Alexei's parameter sets
  double fAUT_Sivers_min, fAUT_Sivers_max; //"min" and "max" values from Alexei's parameter sets
  //double fAUT_Collins_nucleus, fAUT_Sivers_nucleus; //These will be the effective asymmetries for the target nucleus, 
  //Now, since we know on an event-by-event basis which nucleon was the struck nucleon, these can be the (undiluted) asymmetries for the struck nucleon
  // The ACTUAL asymmetry that we generate for any given event (in terms of cross section)
  // will be diluted by the "effective polarization" of the proton and neutron in Helium-3.
  // If we have a vector-polarized deuterium target, we basically assume that deuteron vector
  // polarization = proton + neutron polarization; i.e., both proton and neutron are polarized
  // to the same degree along the same direction (this is an approximation)
  
  //Define additional kinematic quantities for SIDIS:
  double fz, fPh_perp, fphi_h, fphi_S, fTheta_S, fMx;

  //Define additional kinematic quantities for TDIS
  double fxpi, ftpi, fxd, fnu, fya, fy, ff2p, ff2pi, fxa, fPtTDIS, fypi, fSigmaDIS, fSigmaTDIS;
  
  double fBeamCur;
  double fRunTime;
  long    fNevt;   //number of primary events to be generated
  //long    fNtries; //number of "tries" to generate an event (to keep track of efficiency of MC generation).
  double Wfact;
  
  G4SBS::Nucl_t fNuclType, fFinalNucl;
  G4SBS::Targ_t fTargType;
  G4SBS::Kine_t fKineType;
  
  // Which hadron species are we considering for pi/K SIDIS?
  G4SBS::Hadron_t fHadronType; //Currently available: pi+/-/0, K+/-, p/pbar

  double fThMin, fThMax, fPhMin, fPhMax; //Angular generation limits for electron arm

  double fEeMin, fEeMax; //Electron energy generation limits
  double fThMin_had, fThMax_had, fPhMin_had, fPhMax_had; //Angular generation limits for hadron arm 
  double fEhadMin, fEhadMax; //Hadron (total) energy generation limits (for SIDIS case)--Later we will want to add exclusive hadron production.
  double fTargLen, fRasterX, fRasterY, fTargDen; //Targ density is given in atoms or molecules/unit volume
  double fTargZoffset;
  
  double fBeamSpotSize;
  double fCircularRasterRadius;
  //double fTargRadLen; //Radiation length of target material, regardless of thickness
  double fTargRadLen; //Radiation length of target material
  double fTargUpstreamWindowRadLen;
  double fTargZatomic; //atomic number of target for purposes of any bremsstrahlung calculations:
  // set<G4String> G4TargetMaterialNames; 
 
  // D. Flay (8/25/20).  beam pointing 
  double fBeamOffsetX,fBeamOffsetY;

  // D. Flay (10/15/20). beam angle 
  double fBeamAngleX,fBeamAngleY,fBeamAngleZ;  
 
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

  //void InitializePionPhotoLimits(G4SBS::Nucl_t);
  
  //New parameters for Bremsstrahlung generation and photoproduction:
  bool GeneratePionPhotoproduction( G4SBS::Nucl_t, G4LorentzVector, G4LorentzVector ); //exclusive pion photoproduction:
  //G4bool fConstantsInitialized;
  
  G4LorentzVector GetInitialNucl( G4SBS::Targ_t, G4SBS::Nucl_t );
  
  bool GenerateElastic( G4SBS::Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateInelastic( G4SBS::Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateDIS( G4SBS::Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateFlat( G4SBS::Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateBeam( G4SBS::Nucl_t, G4LorentzVector, G4LorentzVector );
  
  bool GenerateSIDIS( G4SBS::Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateTDIS( G4SBS::Nucl_t, G4LorentzVector, G4LorentzVector );
  bool GenerateWiser( G4SBS::Nucl_t, G4LorentzVector, G4LorentzVector );

  bool GenerateGun(); //The "GenerateGun" routine generates generic particles of any type, flat in costheta, phi and p within user-specified limits.
  bool GeneratePythia(); //Generates primaries from a ROOT Tree containing PYTHIA6 events.
  // TDIS AcquMC
  bool GenerateAcquMC(); //Generates primaries from a ROOT Tree containing AcquMC events.
  bool GenerateCosmics(); //Generates muons from the top of the world geometry, directed towards a point in space
  bool GenerateSIMC(); //Generates primaries from a ROOT Tree containing PYTHIA6 events.
  
  //AJRP: June 5, 2021: calculate soffer bounds for transversity calculations:

  G4bool fSofferGridInitialized;
  vector<double> fSofferGrid;

  G4bool fTransversityInitialized;
  vector<double> fTran_a, fTran_b, fTran_n, fTran_m2; //Transversity parameters: size of these is 201*6
  G4bool fSiversInitialized;
  vector<double> fSiv_a, fSiv_b, fSiv_n, fSiv_m2; //Sivers parameters: also 201*6
  
  G4bool fCollinsInitialized;
  vector<double> fColl_a, fColl_b, fColl_n, fColl_m2; //Collins parameters: also 201*6
  
  void SofferBound( G4double x, G4double Q2, vector<double> &partons );
  void Transversity( G4double x, G4double Q2, vector<double> &partons, int iset=0 );
  void Sivers( G4double x, vector<double> &partons, int iset=0 );
  void Collins( G4double z, vector<double> &partons, int iset=0 );

  double AUT_Collins( G4double x, G4double y, G4double Q2, G4double z, G4double PT, vector<double> pdf_unpol, vector<double> fragfunc_unpol, G4SBS::Nucl_t nucl, G4SBS::Hadron_t had, int iset=0 );
  double AUT_Sivers( G4double x, G4double y, G4double Q2, G4double z, G4double PT, vector<double> pdf_unpol, vector<double> fragfunc_unpol, G4SBS::Nucl_t nucl, int iset=0 ); 
  
  G4double fSIDISkperp2_avg; //default 0.25 GeV^2
  G4double fSIDISpperp2_avg; //default 0.20 GeV^2
  
  
  // D Flay (10/15/20).  Generate random beam angle based on non-zero file input.  works for beam generator only
  void CalculateBeamAnglesAndPositions(G4double bd_L,std::vector<G4double> &R,std::vector<G4double> &P);  

  G4bool fRejectionSamplingFlag; //Flag to turn on rejection sampling;
  G4double fMaxWeight; //Maximum event weight within generation limits
  G4int fNeventsWeightCheck; //Number of "pre-events" to generate in order to check weights
  //G4bool fRejectionSamplingInitialized; //Flag to indicate whether rejection sampling has been initialized.

  G4bool fInitialized; //consolidate initialization of constant event generator parameters:
  
  double deutpdist( double );
  double he3pdist( G4SBS::Nucl_t, double );

  double f2p (double);
  double f2pi (double, double, double);
  double c12pdist( double );
  double o16pdist( double );
  
  DSS2007FF fFragFunc; //Class to calculate fragmentation functions using DSS2007

  //void LoadTargetData(); //why is this here? Not implemented...

  //TFile *fPythiaFile;
  //TTree *fPythiaTree;
  long fFirstEvent; //option to start at some event index other than zero
  long fchainentry;
  TChain *fPythiaChain;
  Pythia6_tree *fPythiaTree;

  map<G4String, G4double> fPythiaSigma;
  
  G4SBSPythiaOutput fPythiaEvent;

  // // // // HEAD
  G4int fExclPyXSoption; //Flag to choose "Exclusive pythia" event-by-event cross section events

  G4int counter;//(CA)
  
  // TDIS AcquMC
  long fAcquMCchainentry;
  TChain *fAcquMCChain;
  AcquMCTree *fAcquMCTree;
  G4SBSAcquMCOutput fAcquMCEvent;
  // // // // 
  TChain *fSIMCChain;
  simc_tree *fSIMCTree;
  
  G4SBSSIMCOutput fSIMCEvent;
  // // // // 11a33984f47772444ffb08222f8a978d2bee837e

  G4double TriangleFunc(G4double a, G4double b, G4double c );
};

#endif//G4SBSEVENTGEN_HH
