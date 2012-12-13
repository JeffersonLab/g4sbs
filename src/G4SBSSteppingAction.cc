#include "G4SBSSteppingAction.hh"
//#include "G4SBSSteppingActionMessenger.hh"

#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SteppingManager.hh"

G4SBSSteppingAction::G4SBSSteppingAction()
:drawFlag(false)
{
///  new G4SBSSteppingActionMessenger(this);
}

void G4SBSSteppingAction::UserSteppingAction(const G4Step*)
{
}


