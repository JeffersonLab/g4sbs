#include "G4SBSRun.hh"
#include "G4SBSIO.hh"
// Make this appear first!
#include "G4Timer.hh"

#include "G4SBSRunAction.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

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
  Ntries = 0; //Keep track of total number of tries to throw Nevt events:
  fIO->InitializeTree();

  G4SBSRunData *rmrundata = G4SBSRun::GetRun()->GetData();

//  rmrundata->SetBeamE( G4SBSBeamTarget::GetBeamTarget()->fBeamE/GeV );
  rmrundata->SetNthrown( aRun->GetNumberOfEventToBeProcessed() );

  gen_t gendata = fIO->GetGenData();
  rmrundata->SetBeamE( gendata.Ebeam );
  rmrundata->SetBBtheta( gendata.thbb );
  rmrundata->SetBBdist( gendata.dbb );
  rmrundata->SetSBStheta( gendata.thsbs );
  rmrundata->SetSBSdist( gendata.dsbs );
  rmrundata->SetHCALdist( gendata.dhcal );
  rmrundata->SetHCALvoff( gendata.voffhcal );
  rmrundata->SetRICHdist( gendata.drich );
  rmrundata->SetSBSTrackerdist( gendata.dsbstrkr );
  
  rmrundata->Print();

}

void G4SBSRunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;
  //       << " " << *timer << G4endl;
  G4cout << "number of tries = " << Ntries << G4endl;

  G4SBSRun::GetRun()->GetData()->SetNtries( Ntries );

  G4SBSRunData *rmrundata = G4SBSRun::GetRun()->GetData();
  
  rmrundata->CalcNormalization();
  
  fIO->WriteTree();
}

