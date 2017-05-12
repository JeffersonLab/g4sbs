// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class ESSNDetectorMessenger
// Online control of detector configuration via keyboard
// 20/05/13 JRMA adapted from SBS equivalent, under construction

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"



DetectorMessenger::DetectorMessenger(
                                           DetectorConstruction* SBSDet)
:fDetector(SBSDet)
{ 
  fSBSDir = new G4UIdirectory("/RTPC/");
  fSBSDir->SetGuidance("UI commands to configure Hall-A detectors");
  
  fDataDir = new G4UIdirectory("/RTPC/data/");
  fDataDir->SetGuidance("detector setup parameters directory");

  fUpdateCmd = new G4UIcmdWithoutParameter("/RTPC/det/update",this);
  fUpdateCmd->SetGuidance("Update calorimeter geometry.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);

  fOverlapVolCmd = new G4UIcmdWithoutParameter("/RTPC/det/CheckOverlap",this);
  fOverlapVolCmd->SetGuidance("Turn on overlapping volume checking");
  fOverlapVolCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fUseSrcPbCmd = new G4UIcmdWithoutParameter("/RTPC/det/UseSrcPb",this);
  fUseSrcPbCmd->SetGuidance("Insert Pb shield round source");
  fUseSrcPbCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fPlugAppCmd=new G4UIcmdWithAnInteger("/RTPC/det/PlugApperture",this);
  fPlugAppCmd->SetGuidance("Plug aperture with polythene block");
  fPlugAppCmd->SetParameterName("PlugApperture",false);
  fPlugAppCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fTargetMatCmd=new G4UIcmdWithAString("/RTPC/det/targetMaterial",this);
  fTargetMatCmd->SetGuidance("Select the target material");
  fTargetMatCmd->SetParameterName("TargetMaterial",false);
  fTargetMatCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fPlugLengthCmd = new G4UIcmdWithADouble("/RTPC/det/PlugLength",this);
  fPlugLengthCmd->SetGuidance("Set the length of polythene plug in mm");
  fPlugLengthCmd->SetParameterName("PlugLength",false);
  fPlugLengthCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fTargetRadiusCmd = new G4UIcmdWithADouble("/RTPC/det/targetRadius",this);
  fTargetRadiusCmd->SetGuidance("Set the radius of the target in mm");
  fTargetRadiusCmd->SetParameterName("TargetRadius",false);
  fTargetRadiusCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}



DetectorMessenger::~DetectorMessenger()
{
  delete fPlugAppCmd;
  delete fTargetMatCmd;
  delete fUpdateCmd;
  delete fOverlapVolCmd;
  delete fUseSrcPbCmd;
}



void DetectorMessenger::SetNewValue(G4UIcommand* command,G4String newValue)
{ 
  if( command == fUpdateCmd )
    { fDetector->UpdateGeometry(); }
  
  if( command == fUseSrcPbCmd )
    { fDetector->SetUseSrcPb(true); }

  if( command == fOverlapVolCmd )
    { fDetector->SetIsOverlapVol(true); }
  
  if( command == fPlugAppCmd )
    { fDetector->SetPlugApp(fPlugAppCmd->GetNewIntValue(newValue));}
  
  if( command == fTargetMatCmd )
    { fDetector->SetTargetMaterial(newValue);}

  if( command == fPlugLengthCmd )
    { fDetector->SetPlugLength(fPlugLengthCmd->
				    GetNewDoubleValue(newValue));}

  if( command == fTargetRadiusCmd )
    { fDetector->SetTargetRadius(fTargetRadiusCmd->
				    GetNewDoubleValue(newValue));}

}


