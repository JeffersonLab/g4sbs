// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class RunAction
// MC run control
// 20/05/13 JRMA adapted from SBS equivalent, under construction
//
#include "RunAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

RunAction::RunAction()
{
  fEventAction=NULL;
}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{ 
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  fEventAction= (EventAction*)(G4RunManager::GetRunManager()->GetUserEventAction());
  fEventAction->PrepareOutput();
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  fEventAction->CloseOutput();

}

