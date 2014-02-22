#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

#include "G4SBSOpticalPhysics.hh"

//#include "G4VPhysicsConstructor.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTableIterator.hh"

#include "G4ProcessManager.hh"

G4SBSOpticalPhysics::G4SBSOpticalPhysics() : G4VPhysicsConstructor("Optical") 
{
 
  //theWLSProcess                = NULL;
  theScintProcess              = NULL;
  theCerenkovProcess           = NULL;
  theBoundaryProcess           = NULL;
  theAbsorptionProcess         = NULL;
  theRayleighScattering        = NULL;
  //theMieHGScatteringProcess    = NULL;

}

G4SBSOpticalPhysics::~G4SBSOpticalPhysics() { }

#include "G4OpticalPhoton.hh"

void G4SBSOpticalPhysics::ConstructParticle()
{
  G4OpticalPhoton::OpticalPhotonDefinition();
}

#include "G4ProcessManager.hh"
//#include "G4VPhysicsConstructor.hh"
void G4SBSOpticalPhysics::ConstructProcess()
{
  G4cout << "G4SBSOpticalPhysics: Add Optical Physics Processes" << G4endl;

  theCerenkovProcess = new G4Cerenkov();
  theCerenkovProcess->SetMaxNumPhotonsPerStep(300);
  theCerenkovProcess->SetTrackSecondariesFirst(true);

  theScintProcess = new G4Scintillation();
  theScintProcess->SetScintillationYieldFactor(1.);
  theScintProcess->SetScintillationExcitationRatio(0.0);
  theScintProcess->SetTrackSecondariesFirst(true);

  theAbsorptionProcess = new G4OpAbsorption();
  theRayleighScattering = new G4OpRayleigh();
  theBoundaryProcess = new G4OpBoundaryProcess();

  G4ProcessManager *pManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();
  if (!pManager) {
    std::ostringstream o;
    o << "Optical photon without a Process Manager";
    G4Exception("WLSOpticalPhysics::ConstructProcess()","",
		FatalException,o.str().c_str());
  }
  
  pManager->AddDiscreteProcess( theAbsorptionProcess );
  pManager->AddDiscreteProcess( theBoundaryProcess );
  
  G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
  theScintProcess->AddSaturation(emSaturation);
  
  //aParticleIterator->reset();
  G4ParticleTable::G4PTblDicIterator *aParticleIterator = theParticleTable->GetIterator();

  aParticleIterator->reset();

  while ( (*aParticleIterator)() ){
    
    G4ParticleDefinition* particle = aParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    pManager = particle->GetProcessManager();
    if (!pManager) {
      std::ostringstream o;
      o << "Particle " << particleName << "without a Process Manager";
      G4Exception("WLSOpticalPhysics::ConstructProcess()","",
		  FatalException,o.str().c_str());
    }
    
    if(theCerenkovProcess->IsApplicable(*particle)){
      pManager->AddProcess(theCerenkovProcess);
      pManager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
    }
    if(theScintProcess->IsApplicable(*particle)){
      pManager->AddProcess(theScintProcess);
      pManager->SetProcessOrderingToLast(theScintProcess,idxAtRest);
      pManager->SetProcessOrderingToLast(theScintProcess,idxPostStep);
    }
    
  }
}
