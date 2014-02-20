#ifndef G4SBSOpticalPhysics_h
#define G4SBSOpticalPhysics_h 1

#include "globals.hh"

//#include "G4OpWLS.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"

#include "G4OpMieHG.hh"
#include "G4OpRayleigh.hh"
#include "G4OpAbsorption.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4VPhysicsConstructor.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTableIterator.hh"

class G4SBSOpticalPhysics : public G4VPhysicsConstructor 
{
public:
  
  G4SBSOpticalPhysics();
  virtual ~G4SBSOpticalPhysics();
  
  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:
  
  //G4ParticleTable::G4PTblDicIterator *theParticleIterator;

  //G4OpWLS*             theWLSProcess;
  G4Cerenkov*          theCerenkovProcess;
  G4Scintillation*     theScintProcess;
  G4OpAbsorption*      theAbsorptionProcess;
  G4OpRayleigh*        theRayleighScattering;
  //G4OpMieHG*           theMieHGScatteringProcess;
  G4OpBoundaryProcess* theBoundaryProcess;

};

#endif
