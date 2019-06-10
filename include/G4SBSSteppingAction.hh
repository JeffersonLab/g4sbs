
#ifndef G4SBSSteppingAction_h
#define G4SBSSteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "G4SBSDetectorConstruction.hh"
#include "globals.hh"

class G4SBSSteppingAction : public G4UserSteppingAction
{
public:
  G4SBSSteppingAction();
  virtual ~G4SBSSteppingAction(){};

  virtual void UserSteppingAction(const G4Step*);

  void Initialize( G4SBSDetectorConstruction *fdc );
private:
  G4bool drawFlag;

  map<G4String, set<G4String> > fSDboundaryVolumes; //mapping between sensitive detector names and entrance boundary volumes. In this class, the key value is the boundary volume name, the mapped value is a set of unique SD names associated with the boundary volume.
  
public:
  inline void SetDrawFlag(G4bool val)
  { drawFlag = val; };
};

#endif
