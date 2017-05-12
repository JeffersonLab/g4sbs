// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class EventAction
// Data saving procedures event by event
// 05/05/14 JRMA
//
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RootIO.hh"
#include "ArrayHit.hh"
//#include "SBSVisHit.hh"
#include "RunAction.hh"
#include "EventActionMessenger.hh"
#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"
#include "G4RunManager.hh"
#include <iomanip>

Hparm H3p[] = {
  {200, -200, 200, 200, -200, 200, 200, -200, 200},
  {200, -200, 200, 200, -200, 200, 200, -200, 200},
  {200, -200, 200, 200, -200, 200, 200, -200, 200},
  {200, -200, 200, 200, -200, 200, 200, -200, 200},
  {200, -200, 200, 200, -200, 200, 200,    0,8000},
  {200, -200, 200, 200,    0,8000, 200,    0, 200},
  {200, -200, 200, 200, -200, 200, 200,    0, 200},
  {200, -200, 200, 200, -200, 200, 200,    0, 200}
};
Hparm H2p[] = {
  {400, 0, 400, 400, 0, 10, 0,0,0},
  {400, 0, 400, 400, 0, 100, 0,0,0},
  {200, -200, 200, 200, -200, 200, 0,0,0},
  {200, -200, 200, 200, -200, 200, 0,0,0},
  {200, -200, 200, 200, -200, 200, 0,0,0},
  {200, -200, 200, 200, -200, 200, 0,0,0},
  {200, -200, 200, 200, -200, 200, 0,0,0},
  {200, -200, 200, 200, -200, 200, 0,0,0}
};
Hparm H1p[] = {
  {11000, 0000, 11000, 0,0,0,0,0,0},
  {11000, 0000, 11000, 0,0,0,0,0,0},
  {11000, 0000, 11000, 0,0,0,0,0,0},
  {11000, 0000, 11000, 0,0,0,0,0,0},
  {11000, 0000, 11000, 0,0,0,0,0,0},
  {11000, 0000, 11000, 0,0,0,0,0,0},
  {11000, 0000, 11000, 0,0,0,0,0,0},
  {11000, 0000, 11000, 0,0,0,0,0,0}
};


//
//------------------------------------------------------------------------------
EventAction::EventAction(RunAction* run)
:frunAct(run),fdrawFlag("all"),fprintModulo(1),feventMessenger(0)
{
  feventMessenger = new EventActionMessenger(this);
  fIsInteractive=1;
  // hits collections
  fCollID = -1;
  fHitDrawOpt="edep";
  fOutFileName=G4String("NULL");
  fprintModulo=1000;
  fRTPCOut=NULL;
  //f4Vector = NULL;
  fPGA = (PrimaryGeneratorAction*)
    (G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
}
//
//------------------------------------------------------------------------------
EventAction::~EventAction()
{
  delete feventMessenger;
}

//
//------------------------------------------------------------------------------
void EventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
  if (evtNb%fprintModulo ==0) { 
    G4cout <<evtNb%fprintModulo<< "---> Begin of event: " << evtNb << G4endl;
  }
}
//
//------------------------------------------------------------------------------
void EventAction::EndOfEventAction(const G4Event* evt)
{
  //write to the output ntuple if it exists
  //if not need to set file via /ESSN/event/setOutputFile XXX.root
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  if(fRTPCOut){
    fRTPCOut->WriteHit(HCE);
  }
}  
//
//------------------------------------------------------------------------------
G4int EventAction::PrepareOutput(){
  if(fOutFileName == "NULL") {
    G4cout<<"G4int EventAction::PrepareOutput() no output file"
	  <<G4endl;
    G4cout<<"Please supply filename if output required..."<<G4endl;
    G4cout<<"/RTPC/event/SetOutputFile XXX.root"<<G4endl;
    return 0;
  }
  printf("EventAction Output will be written to %s\n",fOutFileName.data());
  //Create output tree
  //This is curently made in the same format as the cbsim output
  if( fPGA->GetRootIO() ){
    fRTPCOut = fPGA->GetRootIO();
    fRTPCOut->SetOutput(fOutFileName.data());
    for(G4int i=0; i<6; i++){ fRTPCOut->SetH3(i,H3p+i); }
    for(G4int i=0; i<3; i++){
      fRTPCOut->SetH1(i,H1p+i);
      fRTPCOut->SetH2(i,H2p+i);
    }
    fRTPCOut->SetHits(fDet->GetNhits(), fDet->GetHitID(), fDet->GetHits());
  }
  else{
    fRTPCOut = new RootIO(fOutFileName.data(),EOutputFile);
    for(G4int i=0; i<3; i++){ fRTPCOut->SetH3(i,H3p+i); }
  }
  return 1;
}
//
//------------------------------------------------------------------------------
void  EventAction::CloseOutput(){
  if(!fRTPCOut) return;
  fRTPCOut->Close();
}
