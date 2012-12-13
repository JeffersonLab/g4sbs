
#ifndef G4SBSRunAction_h
#define G4SBSRunAction_h 1

#include "globals.hh"
#include "G4UserRunAction.hh"

class G4Timer;
class G4Run;
class G4SBSIO;

class G4SBSRunAction : public G4UserRunAction
{
  public:
    G4SBSRunAction();
    ~G4SBSRunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);

    void SetIO( G4SBSIO *io ){ fIO = io; }

  private:
    G4Timer* timer;

    G4SBSIO *fIO;
};

#endif

