#include "G4SBSPhysicsList.hh"

#include "G4LossTableManager.hh"
//#include "G4PhysListFactory.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
//#include "HadronPhysicsQGSP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4OpticalParameters.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4ProcessManager.hh"
//#include "G4DataQuestionaire.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4Transportation.hh"
#include "G4PhysicsListHelper.hh"

G4SBSPhysicsList::G4SBSPhysicsList() : G4VModularPhysicsList() {
  
  G4LossTableManager::Instance();

  G4int verb = 0;
  SetVerboseLevel(verb);

  RegisterPhysics( new G4DecayPhysics(verb) );
  RegisterPhysics( new G4EmStandardPhysics(verb) );
  RegisterPhysics( new G4EmExtraPhysics(verb) );
  RegisterPhysics( new G4HadronElasticPhysics(verb) );
  RegisterPhysics( new G4IonPhysics(verb) );
  RegisterPhysics( new G4NeutronTrackingCut(verb) );
  RegisterPhysics( new G4HadronPhysicsFTFP_BERT(verb) );
  RegisterPhysics( new G4StoppingPhysics(verb) );
  RegisterPhysics( G4SBSOpticalPhysics = new G4OpticalPhysics(verb) );
  RegisterPhysics( new G4StepLimiterPhysics );

  // //G4SBSParticleList = new G4DecayPhysics("decays");
  // G4SBSParticleList = new G4DecayPhysics(verb);

  // G4SBSEMPhysics = new G4EmStandardPhysics();
  
  // G4SBSHadronicPhysics.clear();
  // G4SBSHadronicPhysics.push_back( new G4EmExtraPhysics(verb) );
  // G4SBSHadronicPhysics.push_back( new G4HadronElasticPhysics(verb) );
  // G4SBSHadronicPhysics.push_back( new G4IonPhysics(verb) );
  // G4SBSHadronicPhysics.push_back( new G4NeutronTrackingCut(verb) );
  // G4SBSHadronicPhysics.push_back( new HadronPhysicsFTFP_BERT(verb) );
  // G4SBSHadronicPhysics.push_back( new G4StoppingPhysics(verb) );
  // G4SBSHadronicPhysics.push_back( new G4StepLimiterBuilder(verb) );

  //G4SBSOpticalPhysics = new G4OpticalPhysics(verb);

  G4OpticalParameters *optParam = G4OpticalParameters::Instance();

  optParam->SetDefaults();
  
  // ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->SetScintillationYieldFactor(1.0);
  // ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->SetScintillationExcitationRatio(0.0);
  // ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->SetMaxNumPhotonsPerStep(300);
  // ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->SetMaxBetaChangePerStep(10.0);
  // ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->SetTrackSecondariesFirst( kCerenkov, false );
  // ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->SetTrackSecondariesFirst( kScintillation, false );

  //G4SBSHadronicPhysics.push_back( opt_phys );

  defaultCutValue = 0.7*mm;
  // cutGamma     = defaultCutValue;
  // cutElectron  = defaultCutValue;
  // cutPositron  = defaultCutValue;
  // cutProton    = defaultCutValue;
  
  G4Transportation::EnableUseMagneticMoment(true);

  //Set optical photon production to be OFF by default:
  ToggleCerenkov(false);
  ToggleScintillation(false);
  
}

//Destructor:
G4SBSPhysicsList::~G4SBSPhysicsList(){
  //delete CerenkovProcess;
  //delete AbsorptionProcess;
  //delete RayleighProcess;
  //delete BoundaryProcess;
  ;
}

void G4SBSPhysicsList::ToggleCerenkov(G4bool useckov=true){
  //( (G4OpticalPhysics*) G4SBSOpticalPhysics )->Configure( kCerenkov, useckov );
  G4OpticalParameters *optParam = G4OpticalParameters::Instance();
  optParam->SetProcessActivation( "Cerenkov", useckov );
}

void G4SBSPhysicsList::ToggleScintillation( G4bool usescint=true ){
  //( (G4OpticalPhysics*) G4SBSOpticalPhysics )->Configure( kScintillation, usescint );
  G4OpticalParameters *optParam = G4OpticalParameters::Instance();
  optParam->SetProcessActivation( "Scintillation", usescint );
}

void G4SBSPhysicsList::SetCuts()
{
  // SetParticleCuts( cutGamma, G4Gamma::Gamma() );
  // SetParticleCuts( cutElectron, G4Electron::Electron() );
  // SetParticleCuts( cutPositron, G4Positron::Positron() );
  // SetParticleCuts( cutProton, G4Proton::Proton() );

  SetCutsWithDefault();
}

// void G4SBSPhysicsList::ConstructProcess(){
//   G4VModularPhysicsList::ConstructProcess();

//   AddTransportation();
  
//   G4Transportation::EnableUseMagneticMoment(true);
  
  
// }
