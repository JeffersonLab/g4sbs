#include "G4SBSPhysicsList.hh"

#include "G4LossTableManager.hh"
//#include "G4PhysListFactory.hh"

#include "G4DecayPhysics.hh"
#include "G4EMStandardPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "HadronPhysicsFTFP_BERT.hh"
#include "G4OpticalPhysics.hh"
#include "G4StepLimiterBuilder.hh"
#include "G4ProcessManager.hh"

G4SBSPhysicsList::G4SBSPhysicsList() : G4VModularPhysicsList() {
  
  G4LossTableManager::Instance();

  G4int verb = 0;
  SetVerboseLevel(verb);

  G4SBSParticleList = new G4DecayPhysics("decays");

  G4SBSEMPhysics = new G4EmStandardPhysics();
  
  G4SBSHadronicPhysics.clear();
  G4SBSHadronicPhysics.push_back( new G4EmExtraPhysics(verb) );
  G4SBSHadronicPhysics.push_back( new G4HadronElasticPhysics(verb) );
  G4SBSHadronicPhysics.push_back( new G4IonPhysics(verb) );
  G4SBSHadronicPhysics.push_back( new G4NeutronTrackingCut(verb) );
  G4SBSHadronicPhysics.push_back( new HadronPhysicsFTFP_BERT(verb) );
  G4SBSHadronicPhysics.push_back( new G4StepLimiterBuilder(verb) );

  G4SBSOpticalPhysics = new G4OpticalPhysics(verb);
  
  ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->SetScintillationYieldFactor(1.0);
  ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->SetScintillationExcitationRatio(0.0);
  ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->SetMaxNumPhotonsPerStep(300);
  ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->SetMaxBetaChangePerStep(10.0);
  ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->SetTrackSecondariesFirst( kCerenkov, true );
  ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->SetTrackSecondariesFirst( kScintillation, true );

  //G4SBSHadronicPhysics.push_back( opt_phys );

  defaultCutValue = 0.7*mm;
  cutGamma     = defaultCutValue;
  cutElectron  = defaultCutValue;
  cutPositron  = defaultCutValue;
  cutProton    = defaultCutValue;

}

//Destructor:
G4SBSPhysicsList::~G4SBSPhysicsList(){
  //delete CerenkovProcess;
  //delete AbsorptionProcess;
  //delete RayleighProcess;
  //delete BoundaryProcess;
  ;
}

void G4SBSPhysicsList::SetCuts()
{
  SetParticleCuts( cutGamma, G4Gamma::Gamma() );
  SetParticleCuts( cutElectron, G4Electron::Electron() );
  SetParticleCuts( cutPositron, G4Positron::Positron() );
  SetParticleCuts( cutProton, G4Proton::Proton() );
}

void G4SBSPhysicsList::ConstructParticle(){
  G4SBSParticleList->ConstructParticle();
}

void G4SBSPhysicsList::ConstructProcess(){
  AddTransportation();
  G4SBSEMPhysics->ConstructProcess();
  G4SBSParticleList->ConstructProcess();
  for( size_t i=0; i<G4SBSHadronicPhysics.size(); i++){
    (G4SBSHadronicPhysics[i])->ConstructProcess();
  }
  G4SBSOpticalPhysics->ConstructProcess();
  
}

void G4SBSPhysicsList::SetOpticalPhysicsProcessActive( G4int pindex, G4bool isactive ){
  
  //Get the name of the process in question:
  //G4VProcess *process;
  G4String processname = "";
  
  bool goodindex = true;
  if( pindex < 0 || pindex >= kNoProcess ) {
    goodindex = false;
    return;
  }

  switch(pindex){
  case kCerenkov:
    //process = ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->GetCerenkovProcess();
    processname = "Cerenkov";
    break;
  case kScintillation:
    //process = ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->GetScintillationProcess();
    processname = "Scintillation";
    break;
  case kWLS:
    //process = ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->GetOpWLSProcess();
    processname = "OpWLS";
    break;
  case kAbsorption:
    //process = ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->GetOpAbsorptionProcess();
    processname = "OpAbsorption";
    break;
  case kRayleigh:
    //process = ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->GetOpRayleighProcess();
    processname = "OpRayleigh";
    break;
  case kBoundary:
    //process = ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->GetOpBoundaryProcess();
    processname = "OpBoundary";
    break;
  case kMieHG:
    //process = ( (G4OpticalPhysics*) G4SBSOpticalPhysics )->GetOpMieHGProcess();
    processname = "OpMie";
    break;
  case kNoProcess: //SHOULDN'T be possible to get here:
  default:
    break;
  }

  if( processname != "" ){

    theParticleIterator->reset();  
    
    while( (*theParticleIterator)() ){
      G4ParticleDefinition *particle = theParticleIterator->value();
      G4ProcessManager *pman = particle->GetProcessManager();
      
      G4ProcessVector *processlist = pman->GetProcessList();

      for( int i=0; i<processlist->size(); i++ ){
	if( ( (G4VProcess*) (*processlist)[i] )->GetProcessName() == processname && 
	    ( (G4VProcess*) (*processlist)[i] )->IsApplicable(*particle) ){
	  G4cout << "Setting process " << ( (G4VProcess*) (*processlist)[i] )->GetProcessName() 
		 << " activation to " << isactive << " for particle " << particle->GetParticleName() << G4endl;
	  pman->SetProcessActivation( ( (G4VProcess*) (*processlist)[i] ), isactive );
	}
      }

      
    }
  }
}
