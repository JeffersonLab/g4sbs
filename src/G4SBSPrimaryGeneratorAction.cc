#include "TBuffer.h"
#include "TMatrixTBase.h"
#include "TString.h"
#include "THashTable.h"

#include "G4SBSPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SBSEventGen.hh"
#include "G4SBSIO.hh"
#include "G4SBSRunAction.hh"
#include "sbstypes.hh"
#include "globals.hh"

G4SBSPrimaryGeneratorAction::G4SBSPrimaryGeneratorAction()
{
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle;

  particle = particleTable->FindParticle(particleName="e-");

  particleGun->SetParticleDefinition(particle);

  particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(-40.0*deg),0.0,cos(-40.0*deg)));
  particleGun->SetParticleEnergy(1.0*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));

  sbsgen = new G4SBSEventGen();

  fUseGeantino = false;
}

G4SBSPrimaryGeneratorAction::~G4SBSPrimaryGeneratorAction()
{
  delete particleGun;
}

void G4SBSPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  G4ParticleDefinition* particle;

  ev_t evdata;

  // Several different types of scattering
  // Let's start with e'N elastic

  //  Roll up random values
  int ntries = 1;
  while( !sbsgen->GenerateEvent() ){ ntries++; }

  G4cout << "Got event, ntries = " << ntries << G4endl;

  int ntries_run = RunAction->GetNtries();
  RunAction->SetNtries( ntries_run + ntries );

  evdata = sbsgen->GetEventData();
  fIO->SetEventData(evdata);

  if( !fUseGeantino ){
      particle = particleTable->FindParticle(particleName="e-");
  } else {
      particle = particleTable->FindParticle(particleName="chargedgeantino");
  }

  particleGun->SetParticleDefinition(particle);

  particleGun->SetParticleMomentumDirection(sbsgen->GetElectronP().unit() );
  particleGun->SetParticleEnergy(sbsgen->GetElectronE());
  particleGun->SetParticlePosition(sbsgen->GetV());
	  
  /*
  particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(-40.0*deg),0.0,cos(-40.0*deg)));
  particleGun->SetParticleEnergy(1.0*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
  */

  // Not necessarily kinematically allowed
  particleGun->GeneratePrimaryVertex(anEvent);

  if( sbsgen->GetKine() != kSIDIS ){ //Then we are generating a final nucleon
    switch( sbsgen->GetFinalNucleon() ){
    case kProton:
      particle = particleTable->FindParticle(particleName="proton");
      break;
    case kNeutron:
      particle = particleTable->FindParticle(particleName="neutron");
      break;
    default:
      particle = particleTable->FindParticle(particleName="geantino");
      break;
    } 
    particleGun->SetParticleDefinition(particle);

    // Ensure we're doing something sensible for Geant4
    if( sbsgen->GetNucleonE()-particle->GetPDGMass() > 0.0 ) {
      particleGun->SetParticleMomentumDirection(sbsgen->GetNucleonP().unit() );
      // This is KINETIC energy
      particleGun->SetParticleEnergy(sbsgen->GetNucleonE()-particle->GetPDGMass());
      particleGun->SetParticlePosition(sbsgen->GetV());
    }
  } else { //SIDIS case: generate a final hadron:
    switch( sbsgen->GetHadronType() ){
    case kPiPlus:
      particle = particleTable->FindParticle(particleName="pi+");
      break;
    case kPiMinus:
      particle = particleTable->FindParticle(particleName="pi-");
      break;
    case kPi0:
      particle = particleTable->FindParticle(particleName="pi0");
      break;
    case kKPlus:
      particle = particleTable->FindParticle(particleName="kaon+");
      break;
    case kKMinus:
      particle = particleTable->FindParticle(particleName="kaon-");
      break;
    default:
      particle = particleTable->FindParticle(particleName="pi+");
      break;
    }

    particleGun->SetParticleDefinition( particle );
    if( sbsgen->GetHadronE()-particle->GetPDGMass() > 0.0 ) {
      particleGun->SetParticleMomentumDirection( sbsgen->GetHadronP().unit() );
      particleGun->SetParticleEnergy( sbsgen->GetHadronE()-particle->GetPDGMass() );
      particleGun->SetParticlePosition( sbsgen->GetV() );
    }
  }

  /*
  particleGun->SetParticleMomentumDirection(G4ThreeVector(sin(39.4*deg),0.0,cos(39.4*deg)));
  particleGun->SetParticleEnergy(1.0*GeV);
  particleGun->SetParticlePosition(G4ThreeVector(0.*cm,0.*cm,0.*cm));
  */

  // Only do final nucleon for generators other than
  // the generic beam generator
  if( sbsgen->GetKine() != kBeam ){
      particleGun->GeneratePrimaryVertex(anEvent);
  }


}

G4ParticleGun* G4SBSPrimaryGeneratorAction::GetParticleGun()
{
  return particleGun;
} 

