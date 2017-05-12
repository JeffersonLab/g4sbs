// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class SteppingMessenger
// Control of max time for tracking
// 20/05/13 JRMA adapted from SBS equivalent, under construction
//
#include "SteppingMessenger.hh"
#include "SteppingAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "Randomize.hh"


SteppingMessenger::SteppingMessenger(SteppingAction* SBSGun)
:SBSAction(SBSGun)
{
  fStepDir = new G4UIdirectory("/RTPC/Step/");
  fStepDir->SetGuidance("Stepping control");

  SetMaxTimeCmd =
    new G4UIcmdWithADoubleAndUnit("/RTPC/Step/SetMaxTime",this);
  SetMaxTimeCmd->SetGuidance("Set the max.time for tracking");
  SetMaxTimeCmd->SetParameterName("MaxTime",false);
  SetMaxTimeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  SetMaxTimeCmd->SetUnitCategory("Time");
}

SteppingMessenger::~SteppingMessenger()
{
  delete fStepDir;
  delete SetMaxTimeCmd;
}

void SteppingMessenger::SetNewValue( G4UIcommand* command,
					G4String newValue)
{ 
  if( command == SetMaxTimeCmd )
    { SBSAction->SetMaxTime(SetMaxTimeCmd->GetNewDoubleValue(newValue));}
}



