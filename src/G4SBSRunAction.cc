
// Make this appear first!
#include "G4Timer.hh"

#include "G4SBSRunAction.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4SBSIO.hh"

G4SBSRunAction::G4SBSRunAction()
{
  timer = new G4Timer;
}

G4SBSRunAction::~G4SBSRunAction()
{
  delete timer;
}

void G4SBSRunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  //  timer->Start();
  fIO->InitializeTree();
}

void G4SBSRunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;
  //       << " " << *timer << G4endl;

  fIO->WriteTree();
}

