#include "TBuffer.h"
#include "TString.h"
#include "THashTable.h"
#include "TGraph.h"
#include "TChain.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TChainElement.h"

#include "G4SBSEventGen.hh"
#include "G4RotationMatrix.hh"
#include "G4SBSInelastic.hh"
#include "G4SBSDIS.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4PionZero.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4SystemOfUnits.hh"


#include "wiser_pion.h"

#include <errno.h>

using namespace CLHEP;

G4SBSEventGen::G4SBSEventGen(){
  //As default values, these don't make sense:
  fThMin = 0.0*deg;
  fThMax = 180.00*deg;
  fPhMin = -180.0*deg;
  fPhMax = 180.0*deg;

  //////////////////////////////////

  fKineType = kElastic;
  fTargType = kH2;
  fTargLen  = 60.0*cm;
  fTargDen  = 10.5*atmosphere/(296.0*kelvin*k_Boltzmann); // This is actually in molecules/unit volume = number density

  //Default SIDIS hadron type to pi+:
  fHadronType = kPiPlus;

  fRasterX  = 0.0*mm;
  fRasterY  = 0.0*mm;

  fBeamE = 2.2*GeV;
  fBeamP = G4ThreeVector( 0.0, 0.0, fBeamE );

  fBeamPol = G4ThreeVector( 0.0, 0.0, 1.0 );
  fhel = 1;

  fCosmPointer = G4ThreeVector(0.0*m, 0.0*m, 0.0*m);
  fPointerZoneRadiusMax = 1.0*m;
  fCosmicsCeilingRadius = 49.0*m;
  fCosmicsMaxAngle = 90.0*deg;
  
  fVert = G4ThreeVector();

  Wfact = 1.0;

  fNevt = 0.0;

  fBeamCur  = 20.0e-6*ampere;     //This is actually in electrons/second
  fRunTime  = 10.0*24.0*3600.0*s; //Ten days 

  fLumi = fTargDen*Wfact * fTargLen * fBeamCur/(e_SI*ampere*second);
  //e_SI  = 1.6e-19
  //ampere * second = coulomb = 1/e_SI --> e_SI * ampere * second = 1
  //fBeamCur is already expressed in electrons/second. So we're good to go
  //fLumi here is now expressed in nucleons * electrons per unit area per second. 

  //Phase space event generation volume. Since we default to elastic, this is simply the solid angle:
  fGenVol = (fPhMax-fPhMin)*(cos(fThMin)-cos(fThMax));
  
  fHCALdist = 17.0*m;

  fToFres = 0.5*ns;

  // init DIS cteq pdf
  initcteqpdf();

  fFragFunc = DSS2007FF();

  fFragFunc.SetGridPath( "." );
  char *G4SBS_env_var = std::getenv("G4SBS");
  if( G4SBS_env_var != NULL ){
    string gridpathname = G4SBS_env_var;
    gridpathname += "/share/DSS2007_GRIDS";
    fFragFunc.SetGridPath(gridpathname);
  }
  
  fEeMin = 0.5*GeV;
  fEeMax = 11.0*GeV;
  fEhadMin = 0.5*GeV;
  fEhadMax = 11.0*GeV;

  // Selecting a broad range of these so they're more inclusive
  fThMin_had = 5.0*deg;
  fThMax_had = 60.0*deg;
  fPhMin_had = 0.0*deg;
  fPhMax_had = 360.0*deg;

  fPythiaChain = NULL;
  fPythiaTree = NULL;
  fchainentry = 0;

  fInitialized = false;
  
  //fRejectionSamplingInitialized = false;
  fRejectionSamplingFlag = false;
  //fMaxWeight = 1.0;
  fMaxWeight = cm2; 
  
  fNeventsWeightCheck = 0;
}


G4SBSEventGen::~G4SBSEventGen(){
  delete fPythiaChain;
  delete fPythiaTree;
}

void G4SBSEventGen::LoadPythiaChain( G4String fname ){
  if( fPythiaChain != NULL ){
    fPythiaChain->Add( fname );
  } else { //First file:
    fPythiaChain = new TChain("Tout");
    fPythiaChain->Add(fname);
    fchainentry = 0;
  } 

  // TFile *ftemp = new TFile( fname, "READ" );
  // TGraph *gtemp;

  // ftemp->GetObject( "graph_sigma", gtemp );

  // if( gtemp ){
  //   //fPythiaEvent.Sigma = gtemp->GetY()[gtemp->GetN()-1];
  //   fPythiaSigma[fname] = gtemp->GetY()[gtemp->GetN()-1]*millibarn;
  // } else {
  //   //fPythiaEvent.Sigma = cm2;
  //   fPythiaSigma[fname] = cm2;
  // }

  // ftemp->Close();
  // delete ftemp;
}

void G4SBSEventGen::Initialize(){
  //Initialize weight factor to convert molecules or atoms number density to number density of nucleons in luminosity calculation:
  switch(fTargType){
  case kH2:
    Wfact = 2.0;
    break;
  case kD2:
    Wfact = 4.0;
    break;
  case kNeutTarg:
    Wfact = 1.0;
    break;
  case kLH2:
    Wfact = 1.0;
    break;
  case kLD2:
    Wfact = 2.0;
    break;
  case k3He:
    Wfact = 3.0;
    break;
  case kCfoil:
    Wfact = 12.0;
    break;
  default:
    Wfact = 1.0;
    break;
  }

  fLumi = fBeamCur / (e_SI*ampere*second) * fTargDen * Wfact * fTargLen;

  //G4cout << "Luminosity = " << fLumi*cm2*s << " cm^{-2} s^{-1}" << G4endl;
  
  fGenVol = (fPhMax - fPhMin)*(cos(fThMin)-cos(fThMax));
  //This expression works for elastic and inelastic, and any other generator that is differential in solid angle only.
  //The inelastic generator returns (Emax-Emin)*dsig/(dE'dOmega_e), so it would be double-counting to multiply by (Emax-Emin)
  //in the calculation of GenVol

  fMaxWeight = cm2; //Maxweight is only relevant when using rejection sampling to produce events distributed according to the cross section.

  //The following generators return cross sections differential in both solid angle AND energy:
  if( fKineType == kDIS || fKineType == kSIDIS || fKineType == kGun ||
      fKineType == kFlat ){
    //All of these generators throw flat in "electron arm" solid angle and energy:
    fGenVol *= (fEeMax - fEeMin);
    fMaxWeight /= GeV;
  }

  if( fKineType == kWiser ){ //Wiser generates flat in "hadron" solid angle and energy; i.e., it uses
    //the "hadron" event generation limits:
    fGenVol = (fPhMax_had - fPhMin_had)*(cos(fThMin_had)-cos(fThMax_had))*(fEhadMax-fEhadMin);
    fMaxWeight /= GeV;
  }

  if( fKineType == kSIDIS ){ //SIDIS generates flat in both electron and hadron solid angle and energy:
    fGenVol *= (fPhMax_had - fPhMin_had)*(cos(fThMin_had)-cos(fThMax_had))*(fEhadMax-fEhadMin);
    fMaxWeight /= GeV;
    G4cout << "Generation volume = " << fGenVol/pow(GeV,2) << "sr*GeV^2" << G4endl;
  }

  if( fRejectionSamplingFlag ){
    InitializeRejectionSampling();
  }
  
  fInitialized = true;
}

bool G4SBSEventGen::GenerateEvent(){
  // Generate initial electron
  // Insert radiative effects: Where are the radiative effects?

  if( !fInitialized ) Initialize();
  
  double Mp = proton_mass_c2;

  G4LorentzVector ei( fBeamP, fBeamE );
  G4LorentzVector ni; 

  // Generate initial nucleon - target dependent
  
  Nucl_t thisnucl;
  //Wfact = 0.0;

  bool success = false; 

  //AJRP: Wfact is now initialized in G4SBSEventGen::InitializeConstants(), invoked at start of run
  switch( fTargType ) {
  case kH2:
    thisnucl = kProton;
    ni = G4LorentzVector(Mp);
    //    Wfact = 1.0;
    // 2 Here because we do molecules/cm3 for density AJRP: is this comment really correct?
    // It appears to be: for a gaseous H2 target, you have two atoms/molecule, so if TargDen is given in molecules/volume, then
    // you have twice as many protons per unit volume. So let's actually change the code to reflect this.
    //Wfact = 2.0;
    break;
  case kD2:
    if( CLHEP::RandFlat::shootInt(2) == 0 ){
      thisnucl = kNeutron;
    } else {
      thisnucl = kProton;
    }

    ni = GetInitialNucl( fTargType, thisnucl );
    //   Wfact = 2.0;
    // AJRP: Based on same considerations discussed above, this should be changed to 4:
    //Wfact = 4.0;
    break;
  case kNeutTarg:
    thisnucl = kNeutron;
    ni = G4LorentzVector(Mp);
    //Wfact = 1.0;
    break;
  case kLH2:
    thisnucl = kProton;
    ni = G4LorentzVector(Mp);
    //Wfact = 1.0;
    //AJRP: for liquid hydrogen, we compute the number density using the mass density, Avogadro's number, and the molar mass, so
    //Wfact = 1 is appropriate here
    break;
  case kLD2:
    if( CLHEP::RandFlat::shootInt(2) == 0 ){
      thisnucl = kNeutron;
    } else {
      thisnucl = kProton;
    }

    ni = GetInitialNucl( fTargType, thisnucl );
    //Wfact = 2.0;
    //AJRP: for liquid deuterium, we compute the number density using the mass density, Avogadro's number, and the molar mass, so
    //Wfact = 2 is appropriate here
    break;
  case k3He:
    if( CLHEP::RandFlat::shootInt(3) == 0 ){
      thisnucl = kNeutron;
    } else {
      thisnucl = kProton;
    }
    ni = GetInitialNucl( fTargType, thisnucl );
    //Wfact = 3.0;
    //AJRP: 3He gas is monatomic, so Wfact = 3 is appropriate here
    break;
  case kCfoil:
    if( CLHEP::RandFlat::shootInt(2) == 0 ){
      thisnucl = kNeutron;
    } else {
      thisnucl = kProton;
    }
    ni = GetInitialNucl( fTargType, thisnucl );
    //Wfact = 3.0;
    //AJRP: 3He gas is monatomic, so Wfact = 3 is appropriate here
    break;
  case kOptics:
    if( CLHEP::RandFlat::shootInt(2) == 0 ){
      thisnucl = kNeutron;
    } else {
      thisnucl = kProton;
    }
    ni = GetInitialNucl( fTargType, thisnucl );
    //Wfact = 3.0;
    //AJRP: 3He gas is monatomic, so Wfact = 3 is appropriate here
    break;
  default:
    thisnucl = kProton;
    ni = G4LorentzVector(Mp);
    //Wfact = 1.0;
  }

  if( fTargType != kOptics ){
    fVert = G4ThreeVector(CLHEP::RandFlat::shoot(-fRasterX/2.0, fRasterX/2.0),
			  CLHEP::RandFlat::shoot(-fRasterY/2.0, fRasterY/2.0),
			  CLHEP::RandFlat::shoot(-fTargLen/2.0, fTargLen/2.0));
  } else { //vertex generation for multi-foil optics target:
    G4double beamx = CLHEP::RandFlat::shoot(-fRasterX/2.0, fRasterX/2.0 );
    G4double beamy = CLHEP::RandFlat::shoot(-fRasterY/2.0, fRasterY/2.0 );

    G4double zfrac = CLHEP::RandFlat::shoot();

    G4double beamz = 0.0;
    
    for( int ifoil=0; ifoil<fNfoils; ifoil++ ){
      if( fFoilZfraction[ifoil] <= zfrac && zfrac < fFoilZfraction[ifoil+1] ){
	//linearly interpolate within this zfrac bin:

	G4double zfoiltemp = fFoilZandThick[ifoil].first;
	G4double foilthicktemp = fFoilZandThick[ifoil].second;

	beamz = zfoiltemp - foilthicktemp/2.0 + foilthicktemp*(zfrac - fFoilZfraction[ifoil])/(fFoilZfraction[ifoil+1]-fFoilZfraction[ifoil]);
	
	break;
      }
    }

    fVert = G4ThreeVector(beamx, beamy, beamz );
    
    //    std::vector<double> zfoil_un
    
  }
    
  fNuclType = thisnucl;

  switch(fKineType){
  case kElastic:
    success = GenerateElastic( thisnucl, ei, ni );
    break;
  case kInelastic:
    success = GenerateInelastic( thisnucl, ei, ni );
    break;
  case kDIS:
    success = GenerateDIS( thisnucl, ei, ni );
    break;
  case kSIDIS:
    success = GenerateSIDIS( thisnucl, ei, ni );
    break;
  case kFlat:
    success = GenerateFlat( thisnucl, ei, ni );
    break;
  case kBeam:
    fVert.setZ( -5.0*m ); // Set at something upstream if just simple beam
    success = GenerateBeam( thisnucl, ei, ni );
    break;
  case kGun:
    success = GenerateGun();
    break;
  case kWiser:
    success = GenerateWiser( thisnucl, ei, ni );
    break;
  case kPYTHIA6:
    success = GeneratePythia();
    break;
  case kCosmics:
    success = GenerateCosmics();
    break;
  default:
    success = GenerateElastic( thisnucl, ei, ni );
    break;
  }

  if( fRejectionSamplingFlag && fInitialized && fSigma > fMaxWeight ) {
    G4cout << "Warning: fSigma > MaxWeight, fSigma/Maxweight = "
	   << fSigma/fMaxWeight << G4endl;
    //fMaxWeight = fSigma;
  }

  // How to normalize? events are thrown flat in phase space, and then accepted or rejected with probability
  // fSigma/fMaxWeight.
  // Overall normalization should be proportional to:
  // xsec * luminosity * phase space volume. The appropriate cross section to use for the normalization is
  // fMaxWeight
  
  if( fRejectionSamplingFlag && fInitialized ){
    success = success && fSigma/fMaxWeight >= CLHEP::RandFlat::shoot();
  }
  
  return success;
}

// bool G4SBSEventGen::GenerateElastic( Nucl_t nucl, G4LorentzVector ei, G4LorentzVector ni ){
//   double Mp = proton_mass_c2;

//   G4ThreeVector pboost = -1.0*(ni.boostVector());

//   G4LorentzVector eip = ei.boost(pboost);
//   ei.boost(-pboost);

//   // Rotation that puts z down eip
//   // Orthogonal vector with z
//   G4ThreeVector rotax = (eip.vect().cross(G4ThreeVector(0.0, 0.0, 1.0))).unit();
//   G4RotationMatrix prot;

//   prot.rotate(-eip.vect().theta(), rotax);

//   eip = G4LorentzVector(eip.e(), G4ThreeVector(0.0, 0.0, eip.e()));

//   G4LorentzVector nip = G4LorentzVector( Mp );
//   // Now we have our boost and way to get back, calculate elastic scattering

//   G4ThreeVector efp3, nfp3, qfp3;
//   G4LorentzVector efp, nfp, q, qf;

//   //  Now do real physics

//   double th = acos( CLHEP::RandFlat::shoot(cos(fThMax), cos(fThMin)) );
//   double ph = CLHEP::RandFlat::shoot(fPhMin, fPhMax );

//   double eprime = (Mp*eip.e())/(Mp + eip.e()*(1.0-cos(th)));

//   /*
//     printf("nucleon p = %f, mass = %f\n", ni.vect().mag()/GeV, ni.m()/GeV);
//     printf("beam e= %f, eprime = %f\n", ei.e()/GeV, eprime/GeV);

//     printf("th = %f, phi = %f\n", th/deg, ph/deg);
//   */

//   efp3.setRThetaPhi(eprime, th, ph );
//   efp = G4LorentzVector( efp3, efp3.mag() );


//   q = eip-efp;

//   nfp3 = q.vect();


//   nfp = G4LorentzVector( nfp3, sqrt(Mp*Mp + nfp3.mag2()));

//   //printf("nucleon f p = %f, ang = %f deg, phi = %f deg, mass = %f\n", nfp.vect().mag()/GeV, nfp.theta()/deg, nfp.phi()/deg,  nfp.m()/GeV);

//   fQ2 = -q.mag2();
    
//   //  Do cross sections and asymmetries

//   double GE, GM, GD;

//   double tau = fQ2/(4.0*Mp*Mp);
//   double alpha = fine_structure_const;

//   GD = pow(1.0 + fQ2/(0.71*GeV*GeV), -2.0);

//   switch( nucl ){
//   case kNeutron:
//     // Our fit
//     GE = (1.520*tau + 2.629*tau*tau + 3.055*tau*tau*tau)*GD/(1.0+5.222*tau+0.040*tau*tau+11.438*tau*tau*tau);
//     // Kelly
//     GM = -1.913*(1.0+2.33*tau)/(1.0 + 14.72*tau + 24.20*tau*tau + 84.1*tau*tau*tau );
//     break;
//   default:
//     // Kelly
//     GE = (1.0-0.24*tau)/(1.0 + 10.98*tau + 12.82*tau*tau + 21.97*tau*tau*tau );
//     // Kelly
//     GM = 2.79*(1.0+0.12*tau)/(1.0 + 10.97*tau + 18.86*tau*tau + 6.55*tau*tau*tau );
//     break;
//   }

//   double dsdx_Mott = pow( cos(th/2.0)*alpha/(2.0*eip.e()*sin(th/2.0)*sin(th/2.0)), 2.0)*hbarc*hbarc;
//   fSigma    = dsdx_Mott*(efp.e()/eip.e())*( (GE*GE+tau*GM*GM)/(1.0+tau) + 2.0*tau*GM*GM*tan(th/2.0)*tan(th/2.0) ); // Dimensions of area


//   fApar  = -(2.0*tau*sqrt(1.0+tau+pow((1.0+tau)*tan(th/2.0),2.0)  )*tan(th/2.0))/
//     (pow(GE/GM,2.0) + (tau + 2.0*tau*(1.0+tau)*pow(tan(th/2.0),2.0)  ));
//   fAperp = -(GE/GM)*2.0*sqrt(tau*(tau+1.0))*tan(th/2.0)/
//     (pow(GE/GM,2.0) + (tau + 2.0*tau*(1.0+tau)*pow(tan(th/2.0),2.0)  ));

//   // Calculate longitudinal / transverse polarization components 
//   double r = GE / GM;
//   double epsilon = pow(1.0 + 2.0*(1.0+tau)*tan(th/2.0)*tan(th/2.0), -1);
//   fPt = ( -fhel*fBeamPol.z()*sqrt( (2.0*epsilon*(1.0-epsilon))/tau) ) * ( r / (1.0+epsilon*r*r/tau) );
//   fPl = ( fhel*fBeamPol.z()*sqrt(1.0-epsilon*epsilon) ) / ( 1.0+epsilon*r*r/tau );

//   // Boost back
    
//   efp3 = prot*efp3;
//   G4LorentzVector ef(efp3, efp3.mag());
//   ef = ef.boost(-pboost);

//   qf = ei - ef;
//   G4ThreeVector qf3 = qf.vect();

//   nfp3 = prot*nfp3;
//   G4LorentzVector nf(nfp3, sqrt(Mp*Mp + nfp3.mag2()) );
//   nf = nf.boost(-pboost);
//   G4ThreeVector nf3 = nf.vect();

//   fPmisspar  = (qf3-nf3)*qf3/qf3.mag();

//   double beta = nf3.mag()/sqrt(nf3.mag2()+Mp*Mp);
//   double tofsm  = beta*fHCALdist/(0.3*m/ns) + CLHEP::RandGauss::shoot(0.0, fToFres);
//   double betasm = fHCALdist/tofsm/(0.3*m/ns);
//   double psm    = Mp*betasm/sqrt(1.0-betasm*betasm);

//   G4ThreeVector nf3sm = (psm/nf3.mag())*nf3;
//   fPmissparSm  = (qfp3-nf3sm)*qf3/qf3.mag();

//   fPmissperp = ((qf3-nf3) - fPmisspar*qf3/qf3.mag()).mag();

//   fW2 = (qf+nip).mag2();
//   fxbj = 1.0;

//   fElectronP = ef.vect();
//   fElectronE = ef.e();

//   fNucleonP = nf.vect();
//   fNucleonE = nf.e();
//   //    printf("nfp_e = %f GeV\n", nfp.e()/GeV);

//   fFinalNucl = nucl;
//   return true;
// }

bool G4SBSEventGen::GenerateElastic( Nucl_t nucl, G4LorentzVector ei, G4LorentzVector ni ){
  G4double Mp = proton_mass_c2;

  G4ThreeVector boost_Nrest = ni.boostVector();

  G4LorentzVector Pisum_lab = ei + ni;

  G4LorentzVector ei_Nrest = ei;
  G4LorentzVector ni_Nrest = ni;

  ei_Nrest.boost( -boost_Nrest );
  ni_Nrest.boost( -boost_Nrest );

  G4double Ebeam_Nrest = ei_Nrest.e();

  G4ThreeVector efp3, nfp3, qfp3;
  G4LorentzVector efp, nfp, q, qf;

  G4double Ebeam_lab = ei.e();

  //Q2 = -(k-k')^2 = 2k dot k' = -2EE'(1-cos theta)
  
  //Generate both theta and phi angles in the LAB frame:
  G4double th = acos( CLHEP::RandFlat::shoot(cos(fThMax), cos(fThMin)) );
  G4double ph = CLHEP::RandFlat::shoot(fPhMin, fPhMax); 

  //unit vector in the direction of the scattered electron in the LAB frame:
  G4ThreeVector kfhat_lab( sin(th)*cos(ph),sin(th)*sin(ph), cos(th) );

  //Outgoing energy of the scattered electron in the LAB frame accounting for the initial nucleon motion (no off-shell or binding energy corrections, just Fermi momentum)
  G4double Eprime_lab = (ei.e()*(ni.e()-ni.pz()))/(ei.e()*(1.-cos(th))+ni.e()-ni.vect().dot(kfhat_lab));
  G4double Pprime_lab = sqrt(pow(Eprime_lab,2)-ei.m2());
  
  G4ThreeVector kf_lab = Pprime_lab*kfhat_lab;

  //Four-momentum of scattered electron in the LAB frame:
  G4LorentzVector ef_lab( kf_lab, Eprime_lab );

  //q vector in the LAB frame:
  G4LorentzVector q_lab = ei - ef_lab;
  G4double Q2 = -q_lab.m2();

  //Calculate four-momentum of scattered electron boosted to the nucleon REST frame for cross section calculation:
  G4LorentzVector ef_Nrest = ef_lab;
  ef_Nrest.boost( -boost_Nrest );

  G4LorentzVector nf_lab_test = ni + ei - ef_lab;
  
  //Calculate dsigma/dOmega_e in the nucleon rest frame:
  

  //G4LorentzVector q_Nrest = ei_Nrest - ef_Nrest;
  
  fQ2 = Q2;
    
  //  Do cross sections and asymmetries

  G4double GE, GM, GD;

  G4double tau = fQ2/(4.0*Mp*Mp);
  G4double alpha = fine_structure_const;

  G4double th_Nrest = acos( ei_Nrest.vect().unit().dot( ef_Nrest.vect().unit()) );
  
  GD = pow(1.0 + fQ2/(0.71*GeV*GeV), -2.0);

  switch( nucl ){
  case kNeutron:
    // Our fit
    GE = (1.520*tau + 2.629*tau*tau + 3.055*tau*tau*tau)*GD/(1.0+5.222*tau+0.040*tau*tau+11.438*tau*tau*tau);
    // Kelly
    GM = -1.913*(1.0+2.33*tau)/(1.0 + 14.72*tau + 24.20*tau*tau + 84.1*tau*tau*tau );
    break;
  default:
    // Kelly
    GE = (1.0-0.24*tau)/(1.0 + 10.98*tau + 12.82*tau*tau + 21.97*tau*tau*tau );
    // Kelly
    GM = 2.79*(1.0+0.12*tau)/(1.0 + 10.97*tau + 18.86*tau*tau + 6.55*tau*tau*tau );
    break;
  }

  //Differential cross section dsigma/dOmega_e in the nucleon rest frame:
  double dsdx_Mott = pow( cos(th_Nrest/2.0)*alpha/(2.0*ei_Nrest.e()*sin(th_Nrest/2.0)*sin(th_Nrest/2.0)), 2.0)*hbarc*hbarc;
  fSigma    = dsdx_Mott*(ef_Nrest.e()/ei_Nrest.e())*( (GE*GE+tau*GM*GM)/(1.0+tau) + 2.0*tau*GM*GM*tan(th_Nrest/2.0)*tan(th_Nrest/2.0) ); // Dimensions of area


  fApar  = -(2.0*tau*sqrt(1.0+tau+pow((1.0+tau)*tan(th_Nrest/2.0),2.0)  )*tan(th_Nrest/2.0))/
    (pow(GE/GM,2.0) + (tau + 2.0*tau*(1.0+tau)*pow(tan(th_Nrest/2.0),2.0)  ));
  fAperp = -(GE/GM)*2.0*sqrt(tau*(tau+1.0))*tan(th_Nrest/2.0)/
    (pow(GE/GM,2.0) + (tau + 2.0*tau*(1.0+tau)*pow(tan(th_Nrest/2.0),2.0)  ));

  // Calculate longitudinal / transverse polarization components 
  double r = GE / GM;
  double epsilon = pow(1.0 + 2.0*(1.0+tau)*tan(th_Nrest/2.0)*tan(th_Nrest/2.0), -1);
  fPt = ( -fhel*fBeamPol.z()*sqrt( (2.0*epsilon*(1.0-epsilon))/tau) ) * ( r / (1.0+epsilon*r*r/tau) );
  fPl = ( fhel*fBeamPol.z()*sqrt(1.0-epsilon*epsilon) ) / ( 1.0+epsilon*r*r/tau );

  // Boost back

  G4LorentzVector q_Nrest = ei_Nrest - ef_Nrest;
  //G4cout << "Q2 (lab) = " << fQ2/pow(GeV,2) << " GeV^2, Q2 (Nrest) = " << -q_Nrest.m2()/pow(GeV,2) << " GeV^2" << G4endl;
  
  G4LorentzVector nf_Nrest = ni_Nrest + q_Nrest;

  G4LorentzVector nf_lab = nf_Nrest;
  nf_lab.boost( boost_Nrest );

  // G4cout << "Lab frame nucleon momentum boosted from rest frame = " << nf_lab << G4endl;

  // G4cout << "Lab frame nucleon momentum from lab frame two-body kinematics = " << nf_lab_test << G4endl;
  //fPmisspar  = (qf3-nf3)*qf3/qf3.mag();
  fPmisspar = (q_lab.vect() - nf_lab.vect()).dot( q_lab.vect().unit() );
  
  G4double beta = nf_lab.vect().mag()/nf_lab.e(); //beta = p/e
  double tofsm  = beta*fHCALdist/(0.3*m/ns) + CLHEP::RandGauss::shoot(0.0, fToFres);
  double betasm = fHCALdist/tofsm/(0.3*m/ns);
  double psm    = Mp*betasm/sqrt(1.0-betasm*betasm);

  G4ThreeVector nf3sm = (psm/nf_lab.vect().mag())*nf_lab.vect();
  fPmissparSm  = (q_lab.vect()-nf3sm).dot(q_lab.vect().unit() );

  fPmissperp = ((q_lab.vect()-nf_lab.vect()) - fPmisspar*q_lab.vect().unit() ).mag();

  double costheta_eN_lab = (ei.vect().unit() ).dot( ni.vect().unit() );
  double betaN_lab = ni.beta();
  double gammaN_lab = ni.gamma();
  
  double flux_Nrest = 4.0*ni.m()*ei_Nrest.e();
  double flux_lab = 4.0*ei.e()*ni.e()*sqrt( 2.0*(1.0-betaN_lab*costheta_eN_lab) - pow(gammaN_lab,-2) );

  // G4cout << "Initial nucleon four-momentum = " << ni << G4endl;
  // G4cout << "Q2, GeV^2 = " << fQ2/pow(GeV,2) << G4endl;
  // G4cout << "nucleon rest frame differential cross section = " << fSigma/nanobarn << " nb/sr" << G4endl;

  // G4cout << "Flux factor, nucleon rest frame = " << flux_Nrest << G4endl;
  // G4cout << "Flux factor, lab frame = " << flux_lab << G4endl;
  
  fSigma *= flux_Nrest/flux_lab;

  //G4cout << "Lab frame differential cross section = " << fSigma/nanobarn << " nb/sr" << G4endl;
  
  //fW2 = (q_lab+ni).mag2();

  fW2 = (q_lab+ni_Nrest).m2();

  // G4cout << "fW2 = " << fW2/pow(GeV,2) << G4endl;
  // G4cout << "fW2 (true, lab) = " << (q_lab+ni).m2()/pow(GeV,2) << G4endl;
  fxbj = 1.0;

  fElectronP = ef_lab.vect();
  fElectronE = ef_lab.e();

  fNucleonP = nf_lab.vect();
  fNucleonE = nf_lab.e();
  //    printf("nfp_e = %f GeV\n", nfp.e()/GeV);

  fFinalNucl = nucl;
  return true;
}

// bool G4SBSEventGen::GenerateInelastic( Nucl_t nucl, G4LorentzVector ei, G4LorentzVector ni ){
//   //double minE = 0.1*GeV;
//   //This generator needs clean-up, in particular to correct the cross section calculation for the non-collinear boost to the nucleon rest frame.
  
//   G4double minE = fEeMin;
//   G4double maxE = fEeMax;
  
//   double Mp = proton_mass_c2;
//   double mpi = 0.140*GeV;

//   G4ThreeVector pboost = -1.0*(ni.boostVector());

//   G4LorentzVector eip = ei.boost(pboost);
//   ei.boost(-pboost);

//   // Rotation that puts z down eip
//   // Orthogonal vector with z
//   G4ThreeVector rotax = (eip.vect().cross(G4ThreeVector(0.0, 0.0, 1.0))).unit();
//   G4RotationMatrix prot;

//   prot.rotate(-eip.vect().theta(), rotax);

//   eip = G4LorentzVector(eip.e(), G4ThreeVector(0.0, 0.0, eip.e()));

//   G4LorentzVector nip = G4LorentzVector( Mp );
//   // Now we have our boost and way to get back, calculate elastic scattering

//   G4ThreeVector efp3, nfp3, qfp3;
//   G4LorentzVector efp, nfp, q, qf;

//   //  Now do real physics

//   double th = acos( CLHEP::RandFlat::shoot(cos(fThMax), cos(fThMin)) );
//   double ph = CLHEP::RandFlat::shoot(fPhMin, fPhMax );

//   //double eprime = CLHEP::RandFlat::shoot(minE, eip.e()-mpi);
//   G4double eprime = CLHEP::RandFlat::shoot(minE, maxE );
  
//   /*
//     printf("nucleon p = %f, mass = %f\n", ni.vect().mag()/GeV, ni.m()/GeV);
//     printf("beam e= %f, eprime = %f\n", ei.e()/GeV, eprime/GeV);

//     printf("th = %f, phi = %f\n", th/deg, ph/deg);
//   */

//   efp3.setRThetaPhi(eprime, th, ph );
//   efp = G4LorentzVector( efp3, efp3.mag() );

//   q = eip-efp;

//   G4ThreeVector hrest3;
//   G4LorentzVector hrest;

//   // Invariant mass system
//   hrest = G4LorentzVector( q.vect(), Mp+q.e() );

//   // This is the invariant mass of the system
//   // Let's assume single pion decay from that

//   double W2 = hrest.mag2();

//   if( W2 < pow(Mp + mpi,2.0) ){
//     // Kinematically not so good - abort
//     fSigma = 0.0;
//     fApar  = 0.0;
//     fAperp = 0.0;
//     fFinalNucl = fNuclType;

//     fPmisspar  = -1e9;
//     fPmissparSm  = -1e9;

//     fPmissperp = -1e9;

//     fW2 = W2;
//     fxbj = -1.0;

//     fElectronP = G4ThreeVector();
//     fElectronE = 0.0;

//     fNucleonP = G4ThreeVector();
//     fNucleonE = 0.0;

//     return false;
//   }

//   double W  = sqrt(W2);

//   double thpi = acos( CLHEP::RandFlat::shoot(-1,1) );
//   double phpi = CLHEP::RandFlat::shoot(0.0, 2.0*3.14159);

//   // Working in the hadronic system rest frame, we isotropically
//   // decay the pion - 1/3 of the time we change charge
//   // from charged pion decay (simple Clebsch Gordon coefficients from delta (I=3/2))

//   if( CLHEP::RandFlat::shoot() < 2.0/3.0 ){
//     fFinalNucl = nucl;
//   } else {
//     fFinalNucl = nucl==kProton?kNeutron:kProton;
//   }


//   double ppi = sqrt(pow(W2 - Mp*Mp - mpi*mpi,2.0) - 4.0*mpi*mpi*Mp*Mp)/(2.0*W);

//   double Ecm = sqrt(ppi*ppi+Mp*Mp);

//   G4ThreeVector ncm3;
//   ncm3.setRThetaPhi(sqrt(Ecm*Ecm-Mp*Mp), thpi, phpi);
//   G4LorentzVector ncm(ncm3, Ecm);

//   /*
//     printf("Ecm = %f (pcm %f)\n", Ecm/GeV, ncm3.mag()/GeV);
//     printf("hrest boost = %f %f %f\n",hrest.boostVector().x(), hrest.boostVector().y(), hrest.boostVector().z());
//     printf("hrest boost mag = %f\n",hrest.boostVector().mag());

//     printf("ncm before = %f %f %f\n", ncm.vect().x()/GeV, ncm.vect().y()/GeV, ncm.vect().z()/GeV);
//   */
//   ncm.boost(hrest.boostVector());
//   //printf("ncm after = %f %f %f\n", ncm.vect().x()/GeV, ncm.vect().y()/GeV, ncm.vect().z()/GeV);
//   nfp  = ncm;
//   nfp3 = ncm.vect();

//   //printf("nucleon f p = %f, ang = %f deg, phi = %f deg, mass = %f\n", nfp.vect().mag()/GeV, nfp.theta()/deg, nfp.phi()/deg,  nfp.m()/GeV);

//   fQ2 = -q.mag2();
    
//   //  Do cross sections and asymmetries

//   fSigma = 0.0;
//   if( nucl == kProton ){
//     //	printf("sigma p! %f %f %f\n", eip.e()/GeV, th/deg, eprime/GeV);
//     //fSigma    = sigma_p(eip.e()/GeV, th/rad, eprime/GeV)*((eip.e()-minE-mpi)/GeV)*nanobarn; // Dimensions of area
//     fSigma    = sigma_p(eip.e()/GeV, th/rad, eprime/GeV)*((maxE-minE)/GeV)*nanobarn; // Dimensions of area
//   }
//   if( nucl == kNeutron ){
//     //fSigma    = sigma_n(eip.e()/GeV, th/rad, eprime/GeV)*((eip.e()-minE-mpi)/GeV)*nanobarn; // Dimensions of area
//     fSigma    = sigma_n(eip.e()/GeV, th/rad, eprime/GeV)*((maxE-minE-mpi)/GeV)*nanobarn; // Dimensions of area
//   }
//   //    printf("fSigma = %e\n", fSigma);

//   if( fSigma != fSigma ) fSigma = 0.0;

//   fApar  = 0.0;
//   fAperp = 0.0;

//   // Boost back
    
//   efp3 = prot*efp3;
//   G4LorentzVector ef(efp3, efp3.mag());
//   ef = ef.boost(-pboost);

//   qf = ei - ef;
//   G4ThreeVector qf3 = qf.vect();

//   nfp3 = prot*nfp3;
//   G4LorentzVector nf(nfp3, sqrt(Mp*Mp + nfp3.mag2()) );
//   nf = nf.boost(-pboost);
//   G4ThreeVector nf3 = nf.vect();

//   fPmisspar  = (qf3-nf3)*qf3/qf3.mag();

//   double beta = nf3.mag()/sqrt(nf3.mag2()+Mp*Mp);
//   double tofsm  = beta*fHCALdist/(0.3*m/ns) + CLHEP::RandGauss::shoot(0.0, fToFres);
//   double betasm = fHCALdist/tofsm/(0.3*m/ns);
//   double psm    = Mp*betasm/sqrt(1.0-betasm*betasm);

//   G4ThreeVector nf3sm = (psm/nf3.mag())*nf3;
//   fPmissparSm  = (qf3-nf3sm)*qf3/qf3.mag();

//   fPmissperp = ((qf3-nf3) - fPmisspar*qf3/qf3.mag()).mag();

//   fW2 = (qf+nip).mag2();
//   fxbj = fQ2/(2.0*Mp*qf.e());

//   fElectronP = ef.vect();
//   fElectronE = ef.e();

//   fNucleonP = nf.vect();
//   fNucleonE = nf.e();

//   return true;
// }

bool G4SBSEventGen::GenerateInelastic( Nucl_t nucl, G4LorentzVector ei, G4LorentzVector ni ){
  //double minE = 0.1*GeV;
  //This generator needs clean-up, in particular to correct the cross section calculation for the non-collinear boost to the nucleon rest frame.
  
  G4double minE = fEeMin;
  G4double maxE = fEeMax;
  
  double Mp = proton_mass_c2;
  double mpi = 0.140*GeV;

  G4LorentzVector Pisum_lab = ei + ni;

  G4ThreeVector boost_Nrest = ni.boostVector();

  G4LorentzVector ei_Nrest = ei;
  G4LorentzVector ni_Nrest = ni;

  ei_Nrest.boost( -boost_Nrest );
  ni_Nrest.boost( -boost_Nrest );

  G4double Ebeam_Nrest = ei_Nrest.e();

  //  Now do real physics

  //Generate electron  angles and energy in the LAB frame:
  //These will then be boosted to the nucleon rest frame to compute the differential cross section.
  
  G4double eth = acos( CLHEP::RandFlat::shoot(cos(fThMax), cos(fThMin)) );
  G4double eph = CLHEP::RandFlat::shoot(fPhMin, fPhMax );

  //double eprime = CLHEP::RandFlat::shoot(minE, eip.e()-mpi);
  G4double Eeprime = CLHEP::RandFlat::shoot(minE, maxE );
  G4double Peprime = sqrt(pow(Eeprime,2) - ei.m2() );
  /*
    printf("nucleon p = %f, mass = %f\n", ni.vect().mag()/GeV, ni.m()/GeV);
    printf("beam e= %f, eprime = %f\n", ei.e()/GeV, eprime/GeV);

    printf("th = %f, phi = %f\n", th/deg, ph/deg);
  */

  G4LorentzVector ef_lab( Eeprime, G4ThreeVector( Peprime*sin(eth)*cos(eph),
						  Peprime*sin(eth)*sin(eph),
						  Peprime*cos(eth) ) );

  G4LorentzVector q_lab = ei - ef_lab;

  G4double Q2 = -q_lab.m2();

  // Boost generated outgoing electron four-momentum to the nucleon rest frame:
  G4LorentzVector ef_Nrest = ef_lab;
  ef_Nrest.boost( -boost_Nrest );

  G4LorentzVector q_Nrest = ei_Nrest - ef_Nrest; 

  //Calculate the boosted value of Bjorken x:
  G4double x = -q_Nrest.m2()/(2.0*ni_Nrest.dot( q_Nrest ) );

  //This is P + q evaluated in the nucleon rest frame:
  //G4LorentzVector P_GammaN_Nrest = ni_Nrest + q_Nrest;
  
  //G4ThreeVector boost_GammaN_Nrest = P_GammaN_Nrest.boostVector();

  //Actually, we can boost directly from the virtual photon-nucleon rest frame to the lab frame:
  G4LorentzVector P_GammaN_lab = q_lab + ni; 
  
  G4double W2 = P_GammaN_lab.mag2();

  if( W2 < pow(Mp + mpi,2.0) ){
    // Kinematically not so good - abort
    fSigma = 0.0;
    fApar  = 0.0;
    fAperp = 0.0;
    fFinalNucl = fNuclType;

    fPmisspar  = -1e9;
    fPmissparSm  = -1e9;

    fPmissperp = -1e9;

    fW2 = W2;
    fxbj = -1.0;

    fElectronP = G4ThreeVector();
    fElectronE = 0.0;

    fNucleonP = G4ThreeVector();
    fNucleonE = 0.0;

    return false;
  }

  G4ThreeVector boost_GammaN_lab = P_GammaN_lab.boostVector();
  
  G4double W  = sqrt(W2);

  //Now we isotropically decay the photon-nucleon system assuming an Npi final state:
  
  G4double thpi = acos( CLHEP::RandFlat::shoot(-1,1) );
  G4double phpi = CLHEP::RandFlat::shoot(0.0, CLHEP::twopi);

  // Working in the hadronic system rest frame, we isotropically
  // decay the pion - 1/3 of the time we change charge
  // from charged pion decay (simple Clebsch Gordon coefficients from delta (I=3/2))

  if( CLHEP::RandFlat::shoot() < 2.0/3.0 ){
    fFinalNucl = nucl;
  } else {
    fFinalNucl = nucl==kProton?kNeutron:kProton;
  }

  G4double Epi_GammaNrest = (W2 - pow(mpi,2) - pow(Mp,2))/(2.0*W);
  G4double EN_GammaNrest = W - Epi_GammaNrest;
  G4double PN_GammaNrest = sqrt(pow(EN_GammaNrest,2)-pow(Mp,2));

  //This is the final nucleon momentum in the virtual photon-nucleon rest frame
  //It can be boosted directly to the lab frame: 
  G4LorentzVector Pfnucleon_GammaNrest( EN_GammaNrest,
					PN_GammaNrest *
					G4ThreeVector( sin(thpi)*cos(phpi),
						       sin(thpi)*sin(phpi),
						       cos(thpi) ) );

  G4LorentzVector Pfnucleon_lab = Pfnucleon_GammaNrest; 
  
  Pfnucleon_lab.boost( boost_GammaN_lab );
  
  
  
  fQ2 = Q2;
  fW2 = W2;
  fxbj = x;
  
  //  Do cross sections and asymmetries

  //This gives the cross section in the nucleon rest frame in units of area/unit solid angle:
  fSigma = 0.0;
  if( nucl == kProton ){
    //	printf("sigma p! %f %f %f\n", eip.e()/GeV, th/deg, eprime/GeV);
    //fSigma    = sigma_p(eip.e()/GeV, th/rad, eprime/GeV)*((eip.e()-minE-mpi)/GeV)*nanobarn; // Dimensions of area
    fSigma    = sigma_p(ei_Nrest.e()/GeV, eth/rad, ef_Nrest.e()/GeV)*((maxE-minE)/GeV)*nanobarn; // Dimensions of area
  }
  if( nucl == kNeutron ){
    //fSigma    = sigma_n(eip.e()/GeV, th/rad, eprime/GeV)*((eip.e()-minE-mpi)/GeV)*nanobarn; // Dimensions of area
    fSigma    = sigma_n(ei_Nrest.e()/GeV, eth/rad, ef_Nrest.e()/GeV)*((maxE-minE)/GeV)*nanobarn; // Dimensions of area
  }
  //    printf("fSigma = %e\n", fSigma);

  if( fSigma != fSigma ) fSigma = 0.0;

  fApar  = 0.0;
  fAperp = 0.0;

  // Boost back:
  //Now that we have the DIS cross section in the nucleon rest frame, we have
  //to account for the modification of the flux factor due to the fact that the collision in the lab frame is no longer collinear
  //(this is the only part of the cross section that is not Lorentz-invariant--it transforms like a cross-sectional area!):
  // Ratio F(lab)/F(Nrest) = 4Ee_lab*En_lab*| v_e - v_n |_lab/4M E_e_Nrest
  // |v_e - v_n| = sqrt( v_e^2 + v_n^2 - 2v_e v_n cos( theta_en ) )
  // v_e^2 = 1, v_n^2 = p_n^2/E_n^2, v_e v_n = p_n/E_n:
  // |v_e - v_n| = sqrt( 1 + p_n^2/E_n^2 - 2p_n/E_n * cos( theta_en ) );
  // p_n^2/E_n^2 = 1 - m_n^2/E_n^2
  // |v_e - v_n| = sqrt( 2*(1 - p_n/E_n*cos(theta_en)) - m_n^2/E_n^2 ) = sqrt( 2*(1 - beta_n * cos(thetaen)) - 1/gamman^2);
  // After boosting to the nucleon rest frame, p_n/E_n = 0 and m_n^2/E_n^2 = 1, so |v_e - v_n| is sqrt( 2(1-beta cos( theta_en )) - 1 ) = 1 (or c when working in regular units).
  //Since dsigma ~ 1/F, we need to multiply the ratio of dsigma lab = dsigma Nrest * Frest/Flab

  double costheta_eN_lab = (ei.vect().unit() ).dot( ni.vect().unit() );
  double betaN_lab = ni.beta();
  double gammaN_lab = ni.gamma();

  double flux_Nrest = 4.0*ni.m()*ei_Nrest.e();
  double flux_lab = 4.0*ei.e()*ni.e()*sqrt( 2.0*(1.0-betaN_lab*costheta_eN_lab) - pow(gammaN_lab,-2) );

  fSigma *= flux_Nrest/flux_lab; //The lines above already converted the cross section to GEANT4 units. Now this has dimensions of area, is expressed in the lab frame, and is differential in solid angle only!

  //The last step is to compute the "true" and "smeared" missing energy and momentum quantities in the lab frame by boosting the
  //final nucleon momentum to the lab frame:
  
  // efp3 = prot*efp3;
  // G4LorentzVector ef(efp3, efp3.mag());
  //ef = ef.boost(-pboost);

  //qf = ei - ef;
  //G4ThreeVector qf3 = qf.vect();

  // nfp3 = prot*nfp3;
  // G4LorentzVector nf(nfp3, sqrt(Mp*Mp + nfp3.mag2()) );
  // nf = nf.boost(-pboost);
  // G4ThreeVector nf3 = nf.vect();

  G4ThreeVector Pmiss = q_lab.vect() - Pfnucleon_lab.vect();
  
  fPmisspar  = ( Pmiss ).dot( q_lab.vect().unit() );
  //fPmisspar  = (qf3-nf3)*qf3/qf3.mag();

  // double beta = nf3.mag()/sqrt(nf3.mag2()+Mp*Mp);
  // // double tofsm  = beta*fHCALdist/(0.3*m/ns) + CLHEP::RandGauss::shoot(0.0, fToFres);
  // double betasm = fHCALdist/tofsm/(0.3*m/ns);
  // double psm    = Mp*betasm/sqrt(1.0-betasm*betasm);

  G4double beta = Pfnucleon_lab.beta();
  G4double tofsm = fHCALdist/(beta*0.3*m/ns) + CLHEP::RandGauss::shoot( 0.0, fToFres );
  G4double betasm = fHCALdist/tofsm/(0.3*m/ns);
  G4double psm = Mp*betasm/sqrt(1.0-pow(betasm,2));
  
  //G4ThreeVector nf3sm = (psm/nf3.mag())*nf3;
  // fPmissparSm  = (qf3-nf3sm)*qf3/qf3.mag();
  fPmissparSm = (q_lab.vect() - psm * Pfnucleon_lab.vect().unit() ).dot( q_lab.vect().unit() );
  
  //fPmissperp = ((qf3-nf3) - fPmisspar*qf3/qf3.mag()).mag();
  fPmissperp = (Pmiss - Pmiss.dot( q_lab.vect().unit() )*q_lab.vect().unit() ).mag();
  
  // fW2 = (qf+nip).mag2();
  // fxbj = fQ2/(2.0*Mp*qf.e());

  fElectronP = ef_lab.vect();
  fElectronE = ef_lab.e();

  fNucleonP = Pfnucleon_lab.vect();
  fNucleonE = Pfnucleon_lab.e();

  return true;
}

bool G4SBSEventGen::GenerateDIS( Nucl_t nucl, G4LorentzVector ei, G4LorentzVector ni ){
  //double minE = 0.*GeV;
  G4double minE = fEeMin;
  G4double maxE = fEeMax;
  double Mp = proton_mass_c2;

  // G4ThreeVector pboost = -1.0*(ni.boostVector());
  // Here we are basically assuming DIS collisions on an on-shell nucleon with simplistic Fermi smearing by sampling the nucleon
  // momentum distribution in 3He or 2H!
  // G4LorentzVector eip = ei.boost(pboost); //eip is incident electron 4-mom boosted to nucleon rest frame

  G4LorentzVector Pisum_lab = ei + ni;

  G4ThreeVector boost_Nrest = ni.boostVector();

  G4LorentzVector ei_Nrest = ei;
  G4LorentzVector ni_Nrest = ni;

  ei_Nrest.boost( -boost_Nrest );
  ni_Nrest.boost( -boost_Nrest );

  double Ebeam_Nrest = ei_Nrest.e();

  //Generate electron kinematics, checking whether an event is kinematically allowed:
  
  //Generate electron  angles and energy in the LAB frame:
  //These will then be boosted to the nucleon rest frame to compute the differential cross section.
  
  //Throw flat in costheta and phi:
  double etheta = acos( CLHEP::RandFlat::shoot( cos( fThMax ), cos( fThMin ) ) ); //same as DIS case.
  double ephi = CLHEP::RandFlat::shoot( fPhMin, fPhMax );
  
  //G4cout << "Generated (etheta, ephi) = (" << etheta/deg << ", " << ephi/deg << ")" << G4endl;

  double Eeprime = CLHEP::RandFlat::shoot( fEeMin, fEeMax );
  double Peprime = sqrt(pow(Eeprime,2) - ei.m2() );

  //G4cout << "Generated Eeprime, Peprime = " << Eeprime/GeV << ", " << Peprime/GeV << G4endl;
   
  G4LorentzVector ef_lab( Eeprime, G4ThreeVector( Peprime*sin(etheta)*cos(ephi), Peprime*sin(etheta)*sin(ephi), Peprime*cos(etheta) ) );

  G4LorentzVector q_lab = ei - ef_lab;

  double Q2 = -q_lab.m2();
  
  G4LorentzVector ef_Nrest = ef_lab;
  ef_Nrest.boost( -boost_Nrest );

  G4LorentzVector q_Nrest = ei_Nrest - ef_Nrest;

  double x = -q_Nrest.m2()/(2.0*ni_Nrest.dot( q_Nrest ) );

  G4LorentzVector X_lab = ni + q_lab; //Four-momentum of unobserved final state, to check whether event is kinematically allowed:

  G4LorentzVector Pfsum_lab = X_lab + ef_lab;
  if( Pfsum_lab.m2() > Pisum_lab.m2() || Pfsum_lab.e() > Pisum_lab.e() ||
      x > 1.0 || x < 0.0 ){ //Check energy and momentum conservation for the generated kinematics:
    fSigma = 0.0;
    fHadronE = 0.0;
    fHadronP = G4ThreeVector();
    fElectronE = 0.0;
    fElectronP = G4ThreeVector();
    fWeight = 0.0;
    fxbj = -1.0;
    fz = -1.0;
    //fW2 = (ni_Nrest + q_Nrest).m2();
    return false;
  }
  
  //What is the purpose of this line? There appears to be no purpose
  // if( CLHEP::RandFlat::shoot() < 2.0/3.0 ){
  //   fFinalNucl = nucl;
  // }

  G4double eth_Nrest = acos( ei_Nrest.vect().unit().dot(ef_Nrest.vect().unit()) );
  G4double eprime_Nrest = ef_Nrest.e();
  
  fQ2 = Q2;
    
  //  Do cross sections and asymmetries

  fSigma = 0.0;
  if( nucl == kProton ){
    //	printf("sigma p! %f %f %f\n", eip.e()/GeV, th/deg, eprime/GeV);
    //fSigma    = dissigma_p(ei_Nrest.e()/GeV, th/rad, eprime/GeV)*((eip.e()-minE)/GeV)*nanobarn; // Dimensions of area
    fSigma = dissigma_p( ei_Nrest.e()/GeV, eth_Nrest/rad, eprime_Nrest/GeV ); //This is in units of nb/GeV/sr in the nucleon rest frame!
  }
  if( nucl == kNeutron ){
    //fSigma    = dissigma_n(eip.e()/GeV, th/rad, eprime/GeV)*((eip.e()-minE)/GeV)*nanobarn; // Dimensions of area
    fSigma = dissigma_n( ei_Nrest.e()/GeV, eth_Nrest/rad, eprime_Nrest/GeV ); //This is in units of nb/GeV/sr in the nucleon rest frame!
  }
  //    printf("fSigma = %e\n", fSigma);

  if( fSigma != fSigma ) fSigma = 0.0; //Why?

  //Now that we have the DIS cross section in the nucleon rest frame, we have
  //to account for the modification of the flux factor due to the fact that the collision in the lab frame is no longer collinear
  //(this is the only part of the cross section that is not Lorentz-invariant--it transforms like a cross-sectional area!):
  // Ratio F(lab)/F(Nrest) = 4Ee_lab*En_lab*| v_e - v_n |_lab/4M E_e_Nrest
  // |v_e - v_n| = sqrt( v_e^2 + v_n^2 - 2v_e v_n cos( theta_en ) )
  // v_e^2 = 1, v_n^2 = p_n^2/E_n^2, v_e v_n = p_n/E_n:
  // |v_e - v_n| = sqrt( 1 + p_n^2/E_n^2 - 2p_n/E_n * cos( theta_en ) );
  // p_n^2/E_n^2 = 1 - m_n^2/E_n^2
  // |v_e - v_n| = sqrt( 2*(1 - p_n/E_n*cos(theta_en)) - m_n^2/E_n^2 ) = sqrt( 2*(1 - beta_n * cos(thetaen)) - 1/gamman^2);
  // After boosting to the nucleon rest frame, p_n/E_n = 0 and m_n^2/E_n^2 = 1, so |v_e - v_n| is sqrt( 2(1-beta cos( theta_en )) - 1 ) = 1 (or c when working in regular units).
  //Since dsigma ~ 1/F, we need to multiply the ratio of dsigma lab = dsigma Nrest * Frest/Flab
  double costheta_eN_lab = (ei.vect().unit() ).dot( ni.vect().unit() );
  double betaN_lab = ni.beta();
  double gammaN_lab = ni.gamma();
  
  double flux_Nrest = 4.0*ni.m()*ei_Nrest.e();
  double flux_lab = 4.0*ei.e()*ni.e()*sqrt( 2.0*(1.0-betaN_lab*costheta_eN_lab) - pow(gammaN_lab,-2) );

  //fSigma *= flux_Nrest/flux_lab * nanobarn * ( fEeMax - fEeMin ) / GeV; //Now this is expressed in nb/sr in the LAB frame, and is a per-nucleon cross section!
  //AJRP: energy generation limit is part of phase space generation volume; don't
  // include it here.
  fSigma *= flux_Nrest/flux_lab * nanobarn/GeV; //Now this is expressed in GEANT4 units (MeV/mm) in the LAB frame, and is a per-nucleon cross section!
  
  fApar  = 0.0;
  fAperp = 0.0;

  // Boost back
    
  // efp3 = prot*efp3;
  // G4LorentzVector ef(efp3, efp3.mag());
  // ef = ef.boost(-pboost);

  // qf = ei - ef;
  // G4ThreeVector qf3 = qf.vect();

  fPmisspar  = 1e-3;

  fPmissparSm  = -1e9;

  fPmissperp = 1e-3;

  // fW2 = (qf+nip).mag2();
  // fxbj = fQ2/(2.0*Mp*qf.e());

  fW2 = X_lab.m2(); //Lorentz-invariant
  fxbj = x;         //Lorentz-invariant
  
  /*
    printf("qf.e = %f (%f)\n", qf.e()/GeV, 6.6-ef.e()/GeV);
    printf("ef = %f (%f)\n", ef.e()/GeV, ef.vect().mag()/GeV);
    printf("nip.e = %f, nip.p.mag = %f\n", nip.e()/GeV, nip.vect().mag()/GeV);

    printf("fQ2 = %f and should be %f and %f\n", fQ2/GeV/GeV, -qf.mag2()/GeV/GeV, 2.0*6.6*GeV*ef.e()*(1.0-cos(ef.vect().theta()))/GeV/GeV );
  */

  fElectronP = ef_lab.vect();
  fElectronE = ef_lab.e();

  fNucleonP = G4ThreeVector();
  fNucleonE = -1e9;  // This ensures we won't generate a nucleon event

  return true;
}

bool G4SBSEventGen::GenerateSIDIS( Nucl_t nucl, G4LorentzVector ei, G4LorentzVector ni ){
  //Get hadron mass:
  
  double Mh;
  int icharge = 1;
  int ihadron = 0;
  switch( fHadronType ){
  case kPiPlus:
    Mh = G4PionPlus::PionPlusDefinition()->GetPDGMass();
    icharge = 1;
    ihadron = 0;
    break;
  case kPiMinus:
    Mh = G4PionMinus::PionMinusDefinition()->GetPDGMass();
    icharge = -1;
    ihadron = 0;
    break;
  case kPi0:
    Mh = G4PionZero::PionZeroDefinition()->GetPDGMass();
    icharge = 0;
    ihadron = 0;
    break;
  case kKPlus:
    Mh = G4KaonPlus::KaonPlusDefinition()->GetPDGMass();
    icharge = 1;
    ihadron = 1;
    break;
  case kKMinus:
    Mh = G4KaonMinus::KaonMinusDefinition()->GetPDGMass();
    icharge = -1;
    ihadron = 1;
    break;
  case kP:
    Mh = G4Proton::ProtonDefinition()->GetPDGMass();
    icharge = 1;
    ihadron = 2;
    break;
  case kPbar:
    Mh = G4AntiProton::AntiProtonDefinition()->GetPDGMass();
    icharge = -1;
    ihadron = 2;
    break;
  default:
    Mh = G4PionPlus::PionPlusDefinition()->GetPDGMass();
    icharge = 1;
    ihadron = 0;
    break;
  }

  G4LorentzVector Pisum_lab = ei + ni;

  //The nucleon could be a proton or a neutron. It has initial 4-momentum ni:
  //Boost to the nucleon rest frame:
  G4ThreeVector boost_Nrest = ni.boostVector();

  //G4LorentzVector ei_Nrest = ei.boost( -boost_Nrest ); 

  //G4LorentzVector ni_Nrest = ni.boost( -boost_Nrest ); //This should equal (M, 0, 0, 0);
  
  //Just in case, copy ei and ni to ei_Nrest and ni_Nrest before boosting:
  G4LorentzVector ei_Nrest = ei;
  G4LorentzVector ni_Nrest = ni;

  ei_Nrest.boost( -boost_Nrest );
  ni_Nrest.boost( -boost_Nrest );

  double Ebeam_Nrest = ei_Nrest.e();

  //Generate electron and hadron kinematics, checking whether an event is kinematically allowed:
  
  //Generate electron and hadron angles and energies in the LAB frame:
  //These will then be boosted to the nucleon rest frame to compute the differential cross section.
  
  //Throw flat in costheta and phi:
  double etheta = acos( CLHEP::RandFlat::shoot( cos( fThMax ), cos( fThMin ) ) ); //same as DIS case.
  double ephi = CLHEP::RandFlat::shoot( fPhMin, fPhMax );
  
  //G4cout << "Generated (etheta, ephi) = (" << etheta/deg << ", " << ephi/deg << ")" << G4endl;

  double Eeprime = CLHEP::RandFlat::shoot( fEeMin, fEeMax );
  double Peprime = sqrt(pow(Eeprime,2) - ei.m2() );

  //G4cout << "Generated Eeprime, Peprime = " << Eeprime/GeV << ", " << Peprime/GeV << G4endl;
   
  G4LorentzVector ef_lab( Eeprime, G4ThreeVector( Peprime*sin(etheta)*cos(ephi), Peprime*sin(etheta)*sin(ephi), Peprime*cos(etheta) ) );

  double htheta = acos( CLHEP::RandFlat::shoot( cos( fThMax_had ), cos( fThMin_had ) ) );
  double hphi = CLHEP::RandFlat::shoot( fPhMin_had, fPhMax_had );

  double Eh = CLHEP::RandFlat::shoot( fEhadMin, fEhadMax );

  //G4cout << "Generated (Eh, htheta, hphi)=(" << Eh/GeV << ", " << htheta/deg << ", " << hphi/deg << ")" << G4endl;

  //For now we assume that Eh > Mh:
  double Ph = sqrt( pow(Eh,2)-pow(Mh,2) );
  

  //G4cout << "Generated Ph = " << Ph/GeV << G4endl;

  G4LorentzVector Phad_lab( Eh, G4ThreeVector( Ph*sin(htheta)*cos(hphi), Ph*sin(htheta)*sin(hphi), Ph*cos(htheta) ) );

  //Check energy and momentum conservation: 
  //1. the sum of outgoing electron and hadron energies cannot exceed the incoming electron and nucleon energies (assuming the collision takes place on a single nucleon):
  //2. 
  G4LorentzVector Pfsum_lab = Phad_lab + ef_lab;
  //G4LorentzVector Pisum_lab = ei + ni;

  if( Pfsum_lab.m2() > Pisum_lab.m2() || Pfsum_lab.e() > Pisum_lab.e() ){
    fSigma = 0.0;
    fHadronE = 0.0;
    fHadronP = G4ThreeVector();
    fElectronE = 0.0;
    fElectronP = G4ThreeVector();
    fWeight = 0.0;
    fxbj = -1.0;
    fz = -1.0;
    //fW2 = (ni_Nrest + q_Nrest).m2();
    return false;
  }

  
  

  //To compute cross section, boost outgoing electron and hadron momenta to nucleon rest frame:
  
  // G4LorentzVector ef_Nrest = ef_lab.boost( -boost_Nrest );
  // G4LorentzVector Phad_Nrest = Phad_lab.boost( -boost_Nrest );
  
  //Just in case, copy before boosting so we don't modify lab-frame quantities:
  G4LorentzVector ef_Nrest = ef_lab;
  G4LorentzVector Phad_Nrest = Phad_lab;

  ef_Nrest.boost( -boost_Nrest );
  Phad_Nrest.boost( -boost_Nrest );

  //Q2 is Lorentz-invariant and depends only on initial and outgoing electron momenta, which are measured in the lab:
  double Q2 = -(ei - ef_lab).m2();

  // At this moment, we have the boosted four-momenta of the initial electron and the outgoing electron and hadron in the nucleon rest frame
  // Now, let us compute the five-fold differential cross section for SIDIS in this frame. The cross section in the lab frame will be modified 
  // in several ways because the collision is non-collinear:
  //   1. kinematics are modified
  //   2. Flux factor is modified by the non-collinear boost: F = 1/( (2E_A) (2E_B) | v_B - v_A | ) (relative to 1/4ME_e) 
  // Therefore, to obtain the cross section, we evaluate the SIDIS structure functions at the modified values of x (and z and Ph_perp), and then 
  // and we compute the modified flux factor 
  // Ingredients are: 
  // 1. PDFs (get from cteq6 routines)
  // 2. FFs (get from DSS2007 routines)


  G4LorentzVector q_Nrest = ei_Nrest - ef_Nrest; //Four-momentum transfer evaluated in the nucleon rest frame
  double x = -q_Nrest.m2() / (2.0*ni_Nrest.dot( q_Nrest ) );

  //double W2 

  //G4cout << "(x, Q2)=(" << x << ", " << Q2/pow(GeV,2) << ")" << endl;

  if( x >= 1.0 || x <= 0.0 ){ //Kinematically forbidden --> abort:
    fSigma = 0.0;
    fHadronE = 0.0;
    fHadronP = G4ThreeVector();
    fElectronE = 0.0;
    fElectronP = G4ThreeVector();
    fWeight = 0.0;
    fxbj = -1.0;
    fz = -1.0;
    fW2 = (ni_Nrest + q_Nrest).m2();
    return false;
  }
  
  //Get PDFs: sqrt(Q2) has units of energy, we should divide by GeV as argument to CTEQ: 
  double u = cteq_pdf_evolvepdf(__dis_pdf, 1, x, sqrt(Q2)/GeV );
  double d = cteq_pdf_evolvepdf(__dis_pdf, 2, x, sqrt(Q2)/GeV );
  double ubar = cteq_pdf_evolvepdf(__dis_pdf, -1, x, sqrt(Q2)/GeV );
  double dbar = cteq_pdf_evolvepdf(__dis_pdf, -2, x, sqrt(Q2)/GeV );
  double st = cteq_pdf_evolvepdf(__dis_pdf, 3, x, sqrt(Q2)/GeV );
  double sbar = st;
  
  //Gaussian model for transverse momentum widths of quark distribution (kperp) and fragmentation (pperp):
  double kperp2_avg = 0.25 * GeV * GeV;
  double pperp2_avg = 0.20 * GeV * GeV;
  
  //Compute SIDIS kinematic quantities:
  // z = P dot Ph / P dot q:
  double z = ni_Nrest.dot( Phad_Nrest ) / ni_Nrest.dot( q_Nrest ); //This quantity is also Lorentz-invariant

  //G4cout << "z = " << z << endl;

  if( z > 1.0 ){ //Kinematically forbidden --> abort:
    fSigma = 0.0;
    fHadronE = 0.0;
    fHadronP = G4ThreeVector();
    fElectronE = 0.0;
    fElectronP = G4ThreeVector();
    fWeight = 0.0;
    fxbj = -1.0;
    fz = -1.0;
    fW2 = (ni_Nrest + q_Nrest).m2();
    return false ;
  }

  fxbj = x;
  fQ2  = Q2;

  vector<double> Dqh;
  fFragFunc.GetFFs( ihadron, icharge, z, sqrt(Q2)/GeV, Dqh );

  // G4cout << "Got fragmentation functions..." << G4endl;
  
  // for( int iparton=0; iparton<6; iparton++ ){
  //   G4cout << "iparton, z, Q2, Dhq = " << iparton << ", " << z << ", " << Q2/pow(GeV,2) << ", " << Dqh[iparton] << endl;
  // }
  
  //Pperp = ph - (ph dot q) * q/q^2 
  G4ThreeVector phad_Nrest_vect = Phad_Nrest.vect();
  G4ThreeVector q_Nrest_vect = q_Nrest.vect();

  //double qvect2_Nrest = q_Nrest_vect.mag2();

  G4ThreeVector Phad_perp = phad_Nrest_vect - ( phad_Nrest_vect.dot(q_Nrest_vect)/q_Nrest_vect.mag2() ) * q_Nrest_vect ;

  double Ph_perp = sqrt( Phad_perp.mag2() );

  double b = 1.0/( pow(z,2)*kperp2_avg + pperp2_avg );

  double e_u = 2.0/3.0;
  double e_d = -1.0/3.0;
  double e_s = -1.0/3.0;
  
  G4ThreeVector zaxis = q_Nrest_vect.unit();
  G4ThreeVector yaxis = (zaxis.cross( ei_Nrest.vect().unit() ) ).unit();
  G4ThreeVector xaxis = (yaxis.cross(zaxis) ).unit();

  //This is phi hadron, the azimuthal angle between the hadron production plane and the lepton scattering plane:
  double phi_hadron = atan2( phad_Nrest_vect.dot( yaxis ), phad_Nrest_vect.dot( xaxis ) );
  fphi_h = phi_hadron;

  fW2 = (ni_Nrest + q_Nrest).mag2();
  fMx = (ni_Nrest + q_Nrest - Phad_Nrest ).mag2();
  
  //Compute SIDIS structure function for a proton:
  double H2 = x * b/twopi*exp(-b*pow(Ph_perp,2)) * ( pow(e_u,2) * (u * Dqh[0] + ubar * Dqh[1]) + 
						     pow(e_d,2) * (d * Dqh[2] + dbar * Dqh[3]) + 
						     pow(e_s,2) * (st * Dqh[4] + sbar * Dqh[5]) );
  
  if( nucl == kNeutron ){ //Interchange u and d quarks: the d quark density in a neutron = u quark density in a proton etc.
    H2 = x * b/twopi*exp(-b*pow(Ph_perp,2)) * ( pow(e_u,2) * (d * Dqh[0] + dbar * Dqh[1]) + 
						pow(e_d,2) * (u * Dqh[2] + ubar * Dqh[3]) + 
						pow(e_s,2) * (st * Dqh[4] + sbar * Dqh[5]) );
  }

  double H1 = H2/(2.0*x); //Callan-Gross relation

  double nu_Nrest = q_Nrest.e();
  //Compute the e- scattering angle in the nucleon rest frame:
  G4ThreeVector ki_Nrest = ei_Nrest.vect();
  G4ThreeVector kf_Nrest = ef_Nrest.vect();
  
  double etheta_Nrest = acos( ki_Nrest.unit().dot( kf_Nrest.unit() ) );

  double theta_pq_Nrest = acos( phad_Nrest_vect.unit().dot( q_Nrest_vect.unit() ) );

  double sigsemi = 4.0*pow(fine_structure_const,2)*hbarc_squared * pow(ef_Nrest.e(),2)/pow(Q2,2) * ( 2.0*H1/proton_mass_c2 * pow(sin(etheta_Nrest/2.0),2) + H2/nu_Nrest * pow(cos(etheta_Nrest/2.0),2) );
  double jacobian = 2.0*phad_Nrest_vect.mag2() * cos( theta_pq_Nrest ) / nu_Nrest;
  sigsemi *= jacobian; 
  //This jacobian factor converts the cross section from d5sig/dE'dOmega_e dz dPh_perp^2 dphi_h  to 
  // d5sig/dE'dOmega_e dE_h dOmega_h. 
  
  //Finally, we have the modification of the flux factor (this is the only part of the cross section that is not Lorentz-invariant--it transforms like a cross-sectional area!):
  // Ratio F(lab)/F(Nrest) = 4Ee_lab*En_lab*| v_e - v_n |_lab/4M E_e_Nrest
  // |v_e - v_n| = sqrt( v_e^2 + v_n^2 - 2v_e v_n cos( theta_en ) )
  // v_e^2 = 1, v_n^2 = p_n^2/E_n^2, v_e v_n = p_n/E_n:
  // |v_e - v_n| = sqrt( 1 + p_n^2/E_n^2 - 2p_n/E_n * cos( theta_en ) );
  // p_n^2/E_n^2 = 1 - m_n^2/E_n^2
  // |v_e - v_n| = sqrt( 2*(1 - p_n/E_n*cos(theta_en)) - m_n^2/E_n^2 );
  // After boosting to the nucleon rest frame, p_n/E_n = 0 and m_n^2/E_n^2 = 1, so |v_e - v_n| is sqrt( 2(1-beta cos( theta_en )) - 1 )
  //Since dsigma ~ 1/F, we need to multiply the ratio of dsigma lab = dsigma Nrest * Frest/Flab
  double costheta_eN_lab = (ei.vect().unit() ).dot( ni.vect().unit() );
  double betaN_lab = ni.beta();
  double gammaN_lab = ni.gamma();
  
  double flux_Nrest = 4.0*ni.m()*ei_Nrest.e();
  double flux_lab = 4.0*ei.e()*ni.e()*sqrt( 2.0*(1.0-betaN_lab*costheta_eN_lab) - pow(gammaN_lab,-2) );
  
  sigsemi *= flux_Nrest/flux_lab; //This is the cross section dsig/dEe' dOmega_e dE_h dOmega_h in units of area/energy^2
  
  fSigma = sigsemi;

  // G4cout << "(x, Q2, z, phperp, MX2)=(" << x << ", " << Q2/pow(GeV,2) << ", " << z << ", " 
  // 	 << Ph_perp/GeV << ", " << fMx/pow(GeV,2) << ")" << ", fSigma = "<< fSigma * pow(GeV,2) / pow(cm,2) << " cm^2/GeV^2/sr^2" << G4endl;

  //These are the four-momenta of outgoing hadron and electron needed for generation of primary particles in GEANT4:
  fHadronE = Phad_lab.e();
  fHadronP = Phad_lab.vect();

  fElectronE = ef_lab.e();
  fElectronP = ef_lab.vect();

  //Record additional SIDIS kinematic variables:
  fz = z; 
  fPh_perp = Ph_perp;
  
  return true;
}


bool G4SBSEventGen::GenerateWiser( Nucl_t nucl, G4LorentzVector ei, G4LorentzVector ni ){

  
  double htheta = acos( CLHEP::RandFlat::shoot( cos( fThMax_had ), cos( fThMin_had ) ) );
  double hphi = CLHEP::RandFlat::shoot( fPhMin_had, fPhMax_had );
  double Eh = CLHEP::RandFlat::shoot( fEhadMin, fEhadMax );

  double Mh;
  int icharge, ihadron;

  switch( fHadronType ){
  case kPiPlus:
    Mh = G4PionPlus::PionPlusDefinition()->GetPDGMass();
    icharge = 1;
    ihadron = 0;
    break;
  case kPiMinus:
    Mh = G4PionMinus::PionMinusDefinition()->GetPDGMass();
    icharge = -1;
    ihadron = 0;
    break;
  case kPi0:
    Mh = G4PionZero::PionZeroDefinition()->GetPDGMass();
    icharge = 0;
    ihadron = 0;
    break;
  case kKPlus:
    Mh = G4KaonPlus::KaonPlusDefinition()->GetPDGMass();
    icharge = 1;
    ihadron = 1;
    break;
  case kKMinus:
    Mh = G4KaonMinus::KaonMinusDefinition()->GetPDGMass();
    icharge = -1;
    ihadron = 1;
    break;
  case kP:
    Mh = G4Proton::ProtonDefinition()->GetPDGMass();
    icharge = 1;
    ihadron = 2;
    break;
  case kPbar:
    Mh = G4AntiProton::AntiProtonDefinition()->GetPDGMass();
    icharge = -1;
    ihadron = 2;
    break;
  default:
    Mh = G4PionPlus::PionPlusDefinition()->GetPDGMass();
    icharge = 1;
    ihadron = 0;
    break;
  }

  double Ph = sqrt(Eh*Eh - Mh*Mh);

  G4LorentzVector Phad_lab( Eh, G4ThreeVector( Ph*sin(htheta)*cos(hphi), Ph*sin(htheta)*sin(hphi), Ph*cos(htheta) ) );

  double intrad = 0.05;

  // This is from 130um wall thickness and 7 cm rad len for GE180
  double glasswallradlen = 0.002;

  double targlen_passed = fVert.z() + fTargLen/2.0;

  double rad_len;
  switch( fTargType ){
  case kH2:
    rad_len = glasswallradlen + targlen_passed/(52.*g/cm2)/(fTargDen*g/6.022e23);
  case kNeutTarg:
    rad_len = glasswallradlen;
    break;
  case k3He:
    rad_len = glasswallradlen + targlen_passed/(94.32*3.*g/cm2/4.)/(fTargDen*3*g/6.022e23);
    break;
  case kLH2:
    rad_len = targlen_passed/(63.04*cm/0.071);
    break;
  case kLD2:
    rad_len = targlen_passed/(125.97*cm/0.169);
    break;
  }


  double sigpip = wiser_sigma(ei.e()/GeV, Phad_lab.vect().mag()/GeV, htheta, rad_len*4.0/3.0 + intrad, 0)*nanobarn/GeV;
  double sigpim = wiser_sigma(ei.e()/GeV, Phad_lab.vect().mag()/GeV, htheta, rad_len*4.0/3.0 + intrad, 1)*nanobarn/GeV;


  if( fNuclType == kProton ){
    switch( fHadronType ){
    case kPiPlus:
      fSigma = sigpip;
      break;
    case kPiMinus:
      fSigma = sigpim;
      break;
    case kPi0:
      fSigma = (sigpip +sigpim)*2.0;
      break;
    default:
      fSigma = 0;
      break;
    }
  } else {
    switch( fHadronType ){
    case kPiPlus:
      fSigma = sigpim;
      break;
    case kPiMinus:
      fSigma = sigpip;
      break;
    case kPi0:
      fSigma = (sigpip +sigpim)*2.0;
      break;
    default:
      fSigma = 0;
      break;
    }
  }

  //AJRP: Why is this line here? To be consistent with the other generators, we keep fSigma as the differential cross section,
  // and multiply in the phase space volume at the stage of a rate calculation/normalization; commented out:
  //  fSigma *= (fEhadMax-fEhadMin)*(cos( fThMax_had) - cos( fThMin_had ))*(fPhMax_had-fPhMin_had)/(cos(fThMax)-cos(fThMin) )/(fPhMax-fPhMin);

  fApar  = 0.0;
  fAperp = 0.0;
  fFinalNucl = fNuclType;

  fPmisspar  = -1e9;
  fPmissparSm  = -1e9;

  fPmissperp = -1e9;

  fW2 = 0;
  fxbj = -1.0;

  fElectronP =Phad_lab.vect();
  fElectronE =Phad_lab.e();

  fHadronE = Phad_lab.e();
  fHadronP = Phad_lab.vect();

  return true;
}


bool G4SBSEventGen::GenerateFlat( Nucl_t nucl, G4LorentzVector ei, G4LorentzVector ni){
  // Ignore initial nucleon
  double Mp = proton_mass_c2;

  // Initial angles
  double th      = acos( CLHEP::RandFlat::shoot(cos(fThMax), cos(fThMin)) );
  double ph      = CLHEP::RandFlat::shoot(fPhMin, fPhMax );
  //double eprime  = CLHEP::RandFlat::shoot(0.0, ei.e());
  double eprime = CLHEP::RandFlat::shoot( fEeMin, fEeMax );
  
  G4ThreeVector efp3, nfp3, qfp3;
  efp3.setRThetaPhi(eprime, th, ph );
    
  G4LorentzVector efp = G4LorentzVector( efp3, efp3.mag() );
  G4LorentzVector q = ei-efp;
  G4LorentzVector nfp = ni+q;

  qfp3 = q.vect();
  nfp3 = nfp.vect();

  fQ2 = -q.mag2();
  fPmisspar  = (qfp3-nfp3)*qfp3/qfp3.mag();

  double beta = nfp3.mag()/sqrt(nfp3.mag2()+Mp*Mp);
  double tofsm  = beta*fHCALdist/(0.3*m/ns) + CLHEP::RandGauss::shoot(0.0, fToFres);
  double betasm = fHCALdist/tofsm/(0.3*m/ns);
  double psm    = Mp*betasm/sqrt(1.0-betasm*betasm);

  G4ThreeVector nfp3sm = (psm/nfp3.mag())*nfp3;
  fPmissparSm  = (qfp3-nfp3sm)*qfp3/qfp3.mag();

  fPmissperp = ((qfp3-nfp3) - fPmisspar*qfp3/qfp3.mag()).mag();

  fW2 = (q+ni).mag2();
  fxbj = fQ2/(2.0*Mp*(ei.e()-efp.e()));

  fElectronP = efp.vect();
  fElectronE = efp.e();

  fNucleonP = nfp.vect();
  fNucleonE = sqrt(Mp*Mp + nfp.mag2());

  fSigma    = cm2/GeV; //1.0 cm2/GeV What else to use?
  fApar     = 0.0;
  fAperp    = 0.0;

  fFinalNucl = nucl;

  return true;
}


bool G4SBSEventGen::GenerateBeam( Nucl_t nucl, G4LorentzVector ei, G4LorentzVector ){

  fQ2 = 0.0;
  fPmisspar  = 0.0;

  fPmissparSm  = 0.0;

  fPmissperp = 0.0;

  fW2 = 0.0;
  fxbj = 0.0;

  fElectronP = ei.vect();
  fElectronE = ei.e();

  fNucleonP = G4ThreeVector();
  fNucleonE = proton_mass_c2;

  fSigma    = 1.0;
  fApar     = 0.0;
  fAperp    = 0.0;

  fFinalNucl = nucl;

  return true;
}

bool G4SBSEventGen::GenerateGun(){

  // // G4cout << "Gun generator: (Emin,Emax,thmin,thmax,phmin,phmax)=("
  // // 	 << fEeMin/GeV << ", " << fEeMax/GeV << ", " 
  // // 	 << fThMin/deg << ", " << fThMax/deg << ", " 
  // // 	 << fPhMin/deg << ", " << fPhMax/deg << ")" << G4endl;

  G4double ep = CLHEP::RandFlat::shoot( fEeMin, fEeMax );
  G4double etheta = acos( CLHEP::RandFlat::shoot( cos(fThMax), cos(fThMin) ) );
  G4double ephi   = CLHEP::RandFlat::shoot( fPhMin, fPhMax );

  //Generate a random momentum and vertex and store the results in fElectronP and fVert:
  fElectronP.set( ep*sin(etheta)*cos(ephi), ep*sin(etheta)*sin(ephi), ep*cos(etheta) );
  
  // G4cout << "Gun generator: Actual p, theta, phi = " << ep/GeV << ", " << etheta/deg << ", " << ephi/deg << G4endl;
  // G4cout << "Gun generator: (px, py, pz)=(" << fElectronP.x()/GeV << ", " << fElectronP.y()/GeV << ", " << fElectronP.z()/GeV << ")" << G4endl;

  // fVert.set( CLHEP::RandFlat::shoot( -fRasterX/2.0, fRasterX/2.0 ),
  // 	     CLHEP::RandFlat::shoot( -fRasterY/2.0, fRasterY/2.0 ),
  // 	     CLHEP::RandFlat::shoot( -fTargLen/2.0, fTargLen/2.0 ) );

  return true;
}

double G4SBSEventGen::deutpdist( double p ){
  // Fit to Bernheim data 
  double thisp = p/GeV;
  if( p < 0.0 ) return 0.0;


  double a0 = 1.199e-6/7.8e-10;
  double b0 = -6.0522;
  double c0 = 7.202e2;

  double a1 = 1.6e-9/7.8e-10;
  double b1 = 17.448;

  if( p < 0.048*GeV ){
    return a0*thisp*thisp*exp(-thisp*b0-thisp*thisp*c0);
  } else {
    return a1*exp(-thisp*b1);
  }
}

double G4SBSEventGen::he3pdist( Nucl_t nucl, double p ){
  // Fits to AV18 supplied by Misak

  double thisp = p/GeV; // Work in units of GeV

  if( p < 0.0 ) return 0.0;

  double ap0, bp0, cp0, ap1, bp1, ap2, bp2, cp1, cp2, cp3, cp4, cp5;
  double an0, bn0, cn0, an1, bn1, an2, bn2, cn1, cn2, cn3, cn4, cn5;
    
  // Proton coeffs

  ap0 = 4.17471e+02;
  bp0 = 1.86153e+01;
  cp0 = 5.83664e+01;

  ap1 = 2.24341e+00;
  bp1 = 1.70263e+01;

  ap2 = 2.59093e-03;
  bp2 = -6.88859e-01;

  cp1 = -5.42591e+00;
  cp2 = -1.91745e+01;
  cp3 = 1.25873e+02;
  cp4 = 3.31128e+02;
  cp5 = -1.60996e+03;

  double pp0 = 0.6;
  double pnorm = 0.43;

  // Neutron coeffs
  an0 = 1.31412e+02;
  bn0 = 1.59806e+01;
  cn0 = 4.14437e+01;

  an1 = 9.93423e-01;
  bn1 = 1.51888e+01;

  an2 = 9.32978e-01;
  bn2 = 1.00013e+01;

  cn1 = 7.34485e+00;
  cn2 = -5.69287e+00;
  cn3 = -1.68647e+02;
  cn4 = 5.79143e+01;
  cn5 = 1.59614e+03;

  double pn0 = 0.55;
  double nnorm = 0.185;
    
  // Proton case

  if( nucl == kProton ){
    if( thisp < 0.14 ){
      return ap0*thisp*thisp*exp(-thisp*bp0-thisp*thisp*cp0)/pnorm;
    } else if (thisp < 0.3 ){
      return ap1*exp(-thisp*bp1)/pnorm;
    } else if (thisp < 0.85){
      return ap2*exp(-thisp*bp2)*(1.0 + cp1*(thisp-pp0) + cp2*pow(thisp-pp0,2.0)
				  + cp3*pow(thisp-pp0,3.0) + cp4*pow(thisp-pp0,4.0) + cp5*pow(thisp-pp0,5.0) )/pnorm;
    } else {
      return 0.0;
    }
  }

  // Neutron case
    
  if( nucl == kNeutron ){
    if( thisp < 0.18 ){
      return an0*thisp*thisp*exp(-thisp*bn0-thisp*thisp*cn0)/nnorm;
    } else if (thisp < 0.32 ){
      return an1*exp(-thisp*bn1)/nnorm;
    } else if (thisp < 0.85){
      return an2*exp(-thisp*bn2)*(1.0 + cn1*(thisp-pn0) + cn2*pow(thisp-pn0,2.0)
				  + cn3*pow(thisp-pn0,3.0) + cn4*pow(thisp-pn0,4.0) + cn5*pow(thisp-pn0,5.0) )/nnorm;
    } else {
      return 0.0;
    }
  }


  return 0.0;
}

double G4SBSEventGen::c12pdist( double p ){
  //C12 being symmetric, we don't need to make a distinction between p and n
  //Add code here:
  double thisp = p/GeV;
  if( p < 0.0 ) return 0.0;


  double a0 = 1.199e-6/7.8e-10;
  double b0 = -6.0522;
  double c0 = 7.202e2;

  double a1 = 1.6e-9/7.8e-10;
  double b1 = 17.448;

  if( p < 0.048*GeV ){
    return a0*thisp*thisp*exp(-thisp*b0-thisp*thisp*c0);
  } else {
    return a1*exp(-thisp*b1);
  }
}


G4LorentzVector G4SBSEventGen::GetInitialNucl( Targ_t targ, Nucl_t nucl ){

  double PMAX;
   
  switch( targ ){
  case kLD2:
    PMAX = 0.35*GeV;
    break;
  case k3He:
    PMAX = 0.85*GeV;
    break;
  case kCfoil:
    PMAX = 1.00*GeV;// PMAX to be adjusted for carbon 12
    break;
  default:
    PMAX = 0.0;
    break;
  }
   
  G4ThreeVector p;
  double theta, phi, psample;

  theta = acos( CLHEP::RandFlat::shoot(-1.0,1.0) );
  phi   = CLHEP::RandFlat::shoot(2.0*pi);

  psample = CLHEP::RandFlat::shoot(PMAX);

  if( targ == k3He ){
    while( CLHEP::RandFlat::shoot() > he3pdist( nucl, psample) ){
      psample = CLHEP::RandFlat::shoot(PMAX);
    }
  }
  if( targ == kLD2 ){
    while( CLHEP::RandFlat::shoot() > deutpdist( psample) ){
      psample = CLHEP::RandFlat::shoot(PMAX);
    }
  }
  if( targ == kCfoil ){
    while( CLHEP::RandFlat::shoot() > c12pdist( psample) ){
      psample = CLHEP::RandFlat::shoot(PMAX);
    }
  }

  p.setRThetaPhi( psample, theta, phi );

  return G4LorentzVector( p, sqrt(p.mag2() + pow(proton_mass_c2,2.0) ) );
}


/*
  double G4SBSEventGen::he3pdist( double p, double x0, double s ){
  return p*p*exp( -1.0*pow( p - x0,2.0)/(2.0*s*s));
  }

  G4LorentzVector G4SBSEventGen::GetInitial3He( Nucl_t nucl ){
  double p_WIDTH = 0.0572372*GeV;
  double p_center = -0.014848*GeV;

  double n_WIDTH = 0.0633468*GeV;
  double n_center = -0.0079127*GeV;

  double fNeutronMax = he3pdist( sqrt( 2.0*n_WIDTH*n_WIDTH + n_center*n_center/4.0) + n_center/2.0, n_center, n_WIDTH );
  double fProtonMax  = he3pdist( sqrt( 2.0*p_WIDTH*p_WIDTH + p_center*p_center/4.0) + p_center/2.0, p_center, p_WIDTH );

  double PMAX = 0.30*GeV;

  G4ThreeVector p;
  double theta, phi, psample;

  theta = acos( CLHEP::RandFlat::shoot(-1.0,1.0) );
  phi   = CLHEP::RandFlat::shoot(2.0*pi);

  psample = -1e9;
  if( nucl == kProton ){
  psample = CLHEP::RandFlat::shoot(PMAX);
  while( CLHEP::RandFlat::shoot() > he3pdist( psample, p_center, p_WIDTH )/fProtonMax ){
  psample = CLHEP::RandFlat::shoot(PMAX);
  }
  }

  if( nucl == kNeutron ){
  psample = CLHEP::RandFlat::shoot(PMAX);
  while( CLHEP::RandFlat::shoot() > he3pdist( psample, n_center, n_WIDTH )/fNeutronMax ){
  psample = CLHEP::RandFlat::shoot(PMAX);
  }
  }

  p.setRThetaPhi( psample, theta, phi );

  return G4LorentzVector( p, sqrt(p.mag2() + pow(proton_mass_c2,2.0) ) );

  }
*/
bool G4SBSEventGen::GeneratePythia(){
  
  fPythiaTree->GetEntry(fchainentry++);

  G4String fnametemp = ( (TChain*) fPythiaTree->fChain )->GetFile()->GetName();

  G4double sigmatemp = fPythiaSigma[fnametemp];
  
  if( fchainentry % 1000 == 0 ) G4cout << "Passed event " << fchainentry << " in PYTHIA6 tree" << G4endl;
  
  //Populate the pythiaoutput data structure:
  fPythiaEvent.Clear();
  //fPythiaEvent.Nprimaries = fPythiaTree->Nparticles;
  fPythiaEvent.Sigma = sigmatemp/cm2;
  fPythiaEvent.Ebeam = (*(fPythiaTree->E))[0]*GeV;
  fPythiaEvent.Eprime = (*(fPythiaTree->E))[2]*GeV;
  fPythiaEvent.theta_e = (*(fPythiaTree->theta))[2]*radian;
  fPythiaEvent.phi_e = (*(fPythiaTree->phi))[2]*radian;
  fPythiaEvent.px_e = (*(fPythiaTree->px))[2]*GeV;
  fPythiaEvent.py_e = (*(fPythiaTree->py))[2]*GeV;
  fPythiaEvent.pz_e = (*(fPythiaTree->pz))[2]*GeV;
  fPythiaEvent.vx_e = (*(fPythiaTree->vx))[2]*mm + fVert.x();
  fPythiaEvent.vy_e = (*(fPythiaTree->vy))[2]*mm + fVert.y();
  fPythiaEvent.vz_e = (*(fPythiaTree->vz))[2]*mm + fVert.z();
  fPythiaEvent.Egamma = (*(fPythiaTree->E))[3]*GeV;
  fPythiaEvent.theta_gamma = (*(fPythiaTree->theta))[3]*radian;
  fPythiaEvent.phi_gamma = (*(fPythiaTree->phi))[3]*radian;
  fPythiaEvent.px_gamma = (*(fPythiaTree->px))[3]*GeV;
  fPythiaEvent.py_gamma = (*(fPythiaTree->py))[3]*GeV;
  fPythiaEvent.pz_gamma = (*(fPythiaTree->pz))[3]*GeV;
  fPythiaEvent.vx_gamma = (*(fPythiaTree->vx))[3]*mm + fVert.x();
  fPythiaEvent.vy_gamma = (*(fPythiaTree->vy))[3]*mm + fVert.y();
  fPythiaEvent.vz_gamma = (*(fPythiaTree->vz))[3]*mm + fVert.z();
  fPythiaEvent.Q2 = fPythiaTree->Q2*GeV*GeV;
  fPythiaEvent.xbj = fPythiaTree->xbj;
  fPythiaEvent.y   = fPythiaTree->y;
  fPythiaEvent.W2  = fPythiaTree->W2*GeV*GeV;

  int ngood = 0;

  int ngen = 0;
  
  for( int i=0; i<fPythiaTree->Nparticles; i++ ){
    //Only fill the first four particles (event header info) and final-state particles (primaries to be generated):
    if( i<4 || (*(fPythiaTree->status))[i] == 1 ){
      fPythiaEvent.PID.push_back( (*(fPythiaTree->pid))[i] );
      //fPythiaEvent.genflag.push_back( 0 ); //Later we will generate in GEANT4 according to user-defined cuts!
      fPythiaEvent.Px.push_back( (*(fPythiaTree->px))[i]*GeV );
      fPythiaEvent.Py.push_back( (*(fPythiaTree->py))[i]*GeV );
      fPythiaEvent.Pz.push_back( (*(fPythiaTree->pz))[i]*GeV );
      fPythiaEvent.M.push_back( (*(fPythiaTree->M))[i]*GeV );
      fPythiaEvent.E.push_back( (*(fPythiaTree->E))[i]*GeV );
      fPythiaEvent.P.push_back( sqrt(pow(fPythiaEvent.Px[ngood],2)+pow(fPythiaEvent.Py[ngood],2)+pow(fPythiaEvent.Pz[ngood],2)) ); //units were already set for this
      fPythiaEvent.t.push_back( (*(fPythiaTree->t))[i]*mm/c_light );
      fPythiaEvent.vx.push_back( (*(fPythiaTree->vx))[i]*mm + fVert.x() );
      fPythiaEvent.vy.push_back( (*(fPythiaTree->vy))[i]*mm + fVert.y() );
      fPythiaEvent.vz.push_back( (*(fPythiaTree->vz))[i]*mm + fVert.z() );
      fPythiaEvent.theta.push_back( (*(fPythiaTree->theta))[i]*radian );
      fPythiaEvent.phi.push_back( (*(fPythiaTree->phi))[i]*radian );
      if( (*(fPythiaTree->status))[i] == 1 &&
	  ( (fPythiaEvent.theta[ngood] >= fThMin && fPythiaEvent.theta[ngood] <= fThMax &&
	     fPythiaEvent.phi[ngood] >= fPhMin && fPythiaEvent.phi[ngood] <= fPhMax &&
	     fPythiaEvent.E[ngood] >= fEeMin && fPythiaEvent.E[ngood] <= fEeMax) ||
	    (fPythiaEvent.theta[ngood] >= fThMin_had && fPythiaEvent.theta[ngood] <= fThMax_had &&
	     fPythiaEvent.phi[ngood] >= fPhMin_had && fPythiaEvent.phi[ngood] <= fPhMax_had &&
	     fPythiaEvent.E[ngood] >= fEhadMin && fPythiaEvent.E[ngood] <= fEhadMax) ) ){ 
	fPythiaEvent.genflag.push_back( 1 );
	ngen++;
	//G4cout << "located good event with one or more primary particles within generation limits" << G4endl;
      } else {
	fPythiaEvent.genflag.push_back( 0 );
      }
      ngood++;
    }
  }
  fPythiaEvent.Nprimaries = fPythiaEvent.PID.size();

  if( ngen == 0 ) return false;
  
  return true;
}

ev_t G4SBSEventGen::GetEventData(){
  ev_t data;

  //Why are we calculating these constant quantities every event?
  //Let's do it once, we can do it via the SBS messenger when we invoke
  // /g4sbs/run command
  // double lumin    = fTargDen*Wfact // Nucleons/Volume
  //   *fTargLen       // Nuclei/area
  //   *fBeamCur/(e_SI*ampere*second);
  //AJRP: moved luminosity calculation to Initialize()
    
  // printf("density = %e N/m3\n", fTargDen*m3);
  // printf("density = %e N/cm3\n", fTargDen*cm3);
  // printf("targlen = %f m\n", fTargLen/m);
  // printf("%e e-/s (I = %f uA)\n", fBeamCur/(e_SI*ampere), fBeamCur/(1e-6*ampere) );
  //printf("luminosity = %e Hz/cm2\n", lumin*second*cm2);
  // printf("e_SI = %e, ampere = %f, \n", e_SI);

  //double genvol   = (fPhMax-fPhMin)*(cos(fThMin)-cos(fThMax));
  //AJRP: moved genvol calculation to Initialize()
  double thisrate = fSigma*fLumi*fGenVol/fNevt;

  //Again: moved genvol calculation to Initialize()
  // if( fKineType == kSIDIS ){ //Then fSigma is dsig/dOmega_e dE'_e dOmega_h dE'_h
  //   genvol *= (fPhMax_had - fPhMin_had)*( cos(fThMin_had) - cos(fThMax_had) );
  //   genvol *= (fEeMax - fEeMin);
  //   genvol *= (fEhadMax - fEhadMin);

  //   thisrate = fSigma*lumin*genvol/fNevt;
      
  // }

  data.count  = thisrate*fRunTime;
  data.rate   = thisrate*second;
  //data.solang = genvol/fNevt; 
  data.solang = fGenVol; //Makes no sense to normalize by number of events here.
  // if( fKineType == kSIDIS ){ //convert genvol to units of GeV^2 in SIDIS case
  //   data.solang /= pow(GeV,2); AJRP: moved to InitializeConstants()
  // }
  data.sigma = fSigma/cm2;

  if( fKineType == kDIS || fKineType == kWiser){
    data.sigma = fSigma/cm2*GeV;
    data.solang = fGenVol/GeV;
  }
  
  if( fKineType == kSIDIS ){ //The SIDIS cross section is also differential in e- energy and hadron energy and has units of area/energy^2/sr^2, so we also need to express it in the correct energy units:
    data.sigma = fSigma/cm2*pow(GeV,2);
    data.solang = fGenVol/pow(GeV,2); //The phase space generation volume has units of energy^2 for SIDIS
  }
  data.Aperp  = fAperp;
  data.Apar   = fApar;
  data.Pt     = fPt;
  data.Pl     = fPl;
  data.W2     = fW2/(GeV*GeV);
  data.xbj    = fxbj;
  data.Q2     = fQ2/(GeV*GeV);
  data.th     = fElectronP.theta()/rad;
  data.ph     = fElectronP.phi()/rad;
  data.vx     = fVert.x()/m;
  data.vy     = fVert.y()/m;
  data.vz     = fVert.z()/m;
  data.ep     = fElectronP.mag()/GeV;
  data.np     = fNucleonP.mag()/GeV;
  data.epx    = fElectronP.x()/GeV;
  data.epy    = fElectronP.y()/GeV;
  data.epz    = fElectronP.z()/GeV;
  data.npx    = fNucleonP.x()/GeV;
  data.npy    = fNucleonP.y()/GeV;
  data.npz    = fNucleonP.z()/GeV;
  data.nth    = fNucleonP.theta()/rad;
  data.nph    = fNucleonP.phi()/rad;

  data.z      = fz;
  data.phperp = fPh_perp/GeV;
  data.phih   = fphi_h;
  data.MX     = fMx/pow(GeV,2);

  if( fKineType == kSIDIS ){ //Then replace final nucleon variables with final hadron variables:
    data.np = fHadronP.mag()/GeV;
    data.npx = fHadronP.x()/GeV;
    data.npy = fHadronP.y()/GeV;
    data.npz = fHadronP.z()/GeV;
    data.nth = fHadronP.theta()/rad;
    data.nph = fHadronP.phi()/rad;
      
    switch( fHadronType ){
    case kPiPlus:
      data.hadr = 1;
      break;
    case kPiMinus:
      data.hadr = -1;
      break;
    case kPi0:
      data.hadr = 0;
      break;
    case kKPlus:
      data.hadr = 2;
      break;
    case kKMinus:
      data.hadr = -2;
      break;
    case kP:
      data.hadr = 3;
      break;
    case kPbar:
      data.hadr = -3;
      break;
    default:
      data.hadr = 1;
      break;
    }
  }
  data.pmpar  = fPmisspar/GeV;
  data.pmparsm= fPmissparSm/GeV;
  data.pmperp = fPmissperp/GeV;
    
  switch( fNuclType ){
  case( kProton ):
    data.nucl   = 1;
    break;
  case( kNeutron):
    data.nucl   = 0;
    break;
  default:
    data.nucl   = -1;
    break;
  }
    
  switch( fFinalNucl ){
  case( kProton ):
    data.fnucl   = 1;
    break;
  case( kNeutron):
    data.fnucl   = 0;
    break;
  default:
    data.fnucl   = -1;
    break;
  }
    
  data.earmaccept = 0;
  data.harmaccept = 0;

  return data;
}

void G4SBSEventGen::InitializeRejectionSampling(){

  fRejectionSamplingFlag = fRejectionSamplingFlag &&
    (fKineType == kElastic || fKineType == kInelastic || fKineType == kDIS ||
     fKineType == kSIDIS || fKineType == kWiser );
  
  if( fRejectionSamplingFlag ){
  
    fInitialized = true;
  
    G4cout << "Initializing rejection sampling..." << G4endl;
  
    fMaxWeight = 0.0;

    if( fNeventsWeightCheck<100000 ){
      fNeventsWeightCheck = 100000; 
    }

    fRejectionSamplingFlag = false;
    for( G4int i=0; i<fNeventsWeightCheck; ++i ){
      if( i % 1000 == 0 ) G4cout << "Estimating max. event weight within generation limits, pre-event = " << i
				 << ", max weight = " << fMaxWeight 
				 <<  G4endl;

      //Generate 1 event:
      while( !GenerateEvent() ){}
    
      fMaxWeight = ( fSigma > fMaxWeight ) ? fSigma : fMaxWeight; 
    }  

    if( fKineType == kSIDIS ){
      G4cout << "Initialized Rejection sampling, max. weight = " << fMaxWeight/(nanobarn/steradian/(GeV*GeV))
	     << G4endl;
    } else {
      G4cout << "Initialized Rejection sampling, max. weight = " << fMaxWeight/(nanobarn/steradian) << " nb/sr" << G4endl;
    }
    fRejectionSamplingFlag = true;
  }
}

void G4SBSEventGen::SetCosmicsPointerRadius( G4double radius ){
  fPointerZoneRadiusMax = radius;
  if(fPointerZoneRadiusMax>min(50.0*m-fabs(fCosmPointer.x()),50.0*m-fabs(fCosmPointer.z()))){
    fPointerZoneRadiusMax = min(50.0*m-fabs(fCosmPointer.x()),50.0*m-fabs(fCosmPointer.z()));
    G4cout << "Warning: you can potentially shoot a cosmics to a point out of 'world' boundaries. Cosmics pointer radius value set to (m): " << fPointerZoneRadiusMax/m << G4endl;
  }
  if(fPointerZoneRadiusMax<0){
    fPointerZoneRadiusMax = 0.0*m;
    G4cout << "Warning: your pointer radius is nonsense (<0). New cosmic pointer radius value set to (m): " << fPointerZoneRadiusMax << G4endl;
  }
}

void G4SBSEventGen::UpdateCosmicsCeilingRadius(){
  fCosmicsCeilingRadius = min(50.0*m-fabs(fCosmPointer.x())-fPointerZoneRadiusMax,50.0*m-fabs(fCosmPointer.z())-fPointerZoneRadiusMax);
}

bool G4SBSEventGen::GenerateCosmics(){
  //G4cout << "Cosmics generated !" << endl;
  
  G4double ep = CLHEP::RandFlat::shoot( fEeMin, fEeMax );
  
  G4double radius2 = CLHEP::RandFlat::shoot( 0.0, fPointerZoneRadiusMax*fPointerZoneRadiusMax); //Add param to configure the pointing area
    G4double phi2 = CLHEP::RandFlat::shoot( -180.0*deg, +180*deg);
  
  //cout << "fCosmPointer: x y z: " << fCosmPointer.x() << " " << fCosmPointer.y() << " " << fCosmPointer.z() << endl;
  
  G4double xptr = fCosmPointer.x()+sin(phi2)*sqrt(radius2);
  G4double zptr = fCosmPointer.z()+cos(phi2)*sqrt(radius2);
  
  //cout << " x,z ptr " << xptr << " " << zptr << endl;
  
  G4double costheta2 = CLHEP::RandFlat::shoot( pow(cos(fCosmicsMaxAngle), 2), 1.0); //Add param to configure the pointing area
  //G4double radius = CLHEP::RandFlat::shoot( 0.0, fCosmicsCeilingRadius*fCosmicsCeilingRadius);
  G4double phi = CLHEP::RandFlat::shoot( -180.0*deg, +180*deg);
  
  /*
  G4double xvtx = xptr+fCosmicsCeilingRadius*sin(phi)*sqrt(radius);
  G4double yvtx = fCosmicsCeiling;
  G4double zvtx = zptr+fCosmicsCeilingRadius*cos(phi)*sqrt(radius);
  */
  
  G4double yvtx = fCosmPointer.y()+fCosmicsCeilingRadius*sqrt(costheta2);
  G4double xvtx = xptr+fCosmicsCeilingRadius*sin(phi)*sqrt(1-costheta2);
  G4double zvtx = zptr+fCosmicsCeilingRadius*cos(phi)*sqrt(1-costheta2);
  
  //cout << " x,z vtx " << xvtx << " " << zvtx << endl;
  
  fVert.set(xvtx, yvtx, zvtx); //Add param to configure the ceiling ?
  //fVert.set(-4.526*m, +5.0*m, +17.008*m);//for test
  
  //cout << " fVert x,z  " << fVert.x() << " " << fVert.z() << endl;
  
  double norm = sqrt(pow(xptr-fVert.x(), 2) + pow(fCosmPointer.y()-fVert.y(), 2) + pow(zptr-fVert.z(), 2)); 
  
  fElectronP.set( ep*(xptr-fVert.x())/norm, ep*(fCosmPointer.y()-fVert.y())/norm, ep*(zptr-fVert.z())/norm );
    
  return true;
}

void G4SBSEventGen::InitializePythia6_Tree(){

  TObjArray *FileList = fPythiaChain->GetListOfFiles();
  TIter next(FileList);

  TChainElement *chEl = 0;

  TGraph *gtemp;
  
  while( (chEl = (TChainElement*) next()) ){
    TFile newfile(chEl->GetTitle(),"READ");
    newfile.GetObject("graph_sigma",gtemp);

    if( gtemp ){
      fPythiaSigma[chEl->GetTitle()] = (gtemp->GetY()[gtemp->GetN()-1])*millibarn;
    } else {
      fPythiaSigma[chEl->GetTitle()] = 1.0*cm2; //
    }
    newfile.Close();
  }
  
  fPythiaTree = new Pythia6_tree( fPythiaChain );
  
  if( !fPythiaTree ){
    G4cout << "Failed to initialize PYTHIA6 tree, aborting... " << G4endl;
    exit(-1);
  }
}

void G4SBSEventGen::SetNfoils( int nfoil ){
  fNfoils = nfoil;
  //fZfoil.clear();
  //fThickFoil.clear();
  //fZfoil.resize(nfoil);
  //ThickFoil.resize(nfoil);
  fFoilZandThick.clear();
}

void G4SBSEventGen::SetFoilZandThick( const std::vector<double> foilz, const std::vector<double> foilthick ){
  if( foilz.size() != fNfoils || foilthick.size() != fNfoils ){
    G4cout << "Foil Z and thickness values given don't match number of foils, exiting..." << G4endl;
    exit(-1);
  }
  fFoilZandThick.clear();

  fTotalThickFoil = 0.0;
  
  for( int ifoil=0; ifoil<fNfoils; ifoil++ ){
    fFoilZandThick.push_back( std::make_pair(foilz[ifoil], foilthick[ifoil] ) );
    fTotalThickFoil += foilthick[ifoil];
  }

  //Now sort by Z foil in ascending order:

  std::sort( fFoilZandThick.begin(), fFoilZandThick.end() );

  //Now loop over foils and compute total thickness and fraction:
  fFoilZfraction.clear();
  // fFoilZfraction.resize(fNfoils);

  G4double thicksum = 0.0;
  
  for( int ifoil=0; ifoil<fNfoils; ifoil++ ){
    fFoilZfraction.push_back( thicksum );
    thicksum += fFoilZandThick[ifoil].second / fTotalThickFoil; 
  }
  fFoilZfraction.push_back( 1.0 );
}
