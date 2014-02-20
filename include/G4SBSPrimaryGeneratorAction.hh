
#ifndef G4SBSPrimaryGeneratorAction_h
#define G4SBSPrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;
class G4SBSEventGen;
class G4SBSIO;
class G4SBSRunAction;

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

    void SetUseGeantino(bool b){ fUseGeantino = b; }
  void SetRunAction( G4SBSRunAction *r ){ RunAction = r; }

  private:
    G4ParticleGun* particleGun;
    G4SBSEventGen* sbsgen;
  G4SBSRunAction *RunAction;
    G4SBSIO *fIO;

  //store a pointer to G4SBSRunAction:
  
  

    bool fUseGeantino;
};

#endif


