
#ifndef G4SBSPrimaryGeneratorAction_h
#define G4SBSPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include "sbstypes.hh"//CA

class G4ParticleGun;
class G4ParticleDefinition;
class G4Event;
class G4SBSEventGen;
class G4SBSIO;
class G4SBSRunAction;
//class G4ThreeVector;

class G4SBSPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  G4SBSPrimaryGeneratorAction();
  ~G4SBSPrimaryGeneratorAction();

public:
  void GeneratePrimaries(G4Event* anEvent);
  G4ParticleGun* GetParticleGun();
  void SetIO( G4SBSIO *io ){ fIO = io; }

  G4SBSEventGen *GetEvGen(){ return sbsgen; }

  G4ThreeVector GetGunPolarization(){ return GunPolarization; }
  
  void SetUseGeantino(bool b){ fUseGeantino = b; }
  void SetRunAction( G4SBSRunAction *r ){ RunAction = r; }
  
  void SetParticleType( G4ParticleDefinition *ptype ){ GunParticleType = ptype; }
  void SetParticleName( G4String pname ){ GunParticleName = pname; }
  void SetGunPolarization( G4ThreeVector S ){ GunPolarization = S; }
  
 // I don't like this solution. I would like that the information came from the TDIS generator
  void tSetHadronType(G4SBS::Hadron_t h ){tHadronType = h; }//CA
private:

  G4SBS::Hadron_t tHadronType; //CA
  
  G4ParticleGun* particleGun;
  G4String GunParticleName;
  G4ParticleDefinition *GunParticleType; //Which particle type will the gun use?
  G4SBSEventGen* sbsgen;
  G4SBSRunAction *RunAction;
  G4SBSIO *fIO;
 
  G4ThreeVector GunPolarization;

  bool fUseGeantino;
};

#endif


