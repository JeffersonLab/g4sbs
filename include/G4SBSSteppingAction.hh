
#ifndef G4SBSSteppingAction_h
#define G4SBSSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"

class G4SBSSteppingAction : public G4UserSteppingAction
{
  public:
    G4SBSSteppingAction();
    virtual ~G4SBSSteppingAction(){};

    virtual void UserSteppingAction(const G4Step*);

  private:
    G4bool drawFlag;

  public:
    inline void SetDrawFlag(G4bool val)
    { drawFlag = val; };
};

#endif
