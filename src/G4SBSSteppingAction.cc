#include "G4SBSSteppingAction.hh"
//#include "G4SBSSteppingActionMessenger.hh"

#include "G4VVisManager.hh"
#include "G4Polyline.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SteppingManager.hh"

#include "G4StepPoint.hh"
#include "G4Material.hh"
#include "G4VProcess.hh"

G4SBSSteppingAction::G4SBSSteppingAction()
:drawFlag(false)
{
///  new G4SBSSteppingActionMessenger(this);
}

void G4SBSSteppingAction::UserSteppingAction(const G4Step* s)
{
  // G4StepPoint *prestep = s->GetPreStepPoint();

  // G4Material *mat = prestep->GetMaterial();
  
  // G4cout << "Before step, material = " << mat->GetName() << G4endl;
  
  // G4StepPoint *poststep = s->GetPostStepPoint();
  
  // mat = poststep->GetMaterial();
  
  // G4cout << "After step, material = " << mat->GetName() << G4endl;

  // G4MaterialPropertiesTable *mpt = mat->GetMaterialPropertiesTable();

  // if( mpt != NULL ){
  //   //mpt->DumpTable();
  // }
  
  // G4cout << poststep->GetProcessDefinedStep()->GetProcessName() << G4endl;

  //  G4cout << process->GetProcessName() << G4endl;

}


