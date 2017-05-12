// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class PrimaryGeneratorMessenger
// Online control of particle gun
// 20/05/13 JRMA adapted from SBS equivalent, under construction

#include "PrimaryGeneratorMessenger.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3Vector.hh"
#include "Randomize.hh"


PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(
                                          PrimaryGeneratorAction* Gun)
:Action(Gun)
{
  gunDir = new G4UIdirectory("/RTPC/Gen/");
  gunDir->SetGuidance("PrimaryGenerator control");

  SetRootInCmd = new G4UIcmdWithAString("/RTPC/Gen/Input",this);
  SetRootInCmd->SetGuidance("Set input ROOT file");
  SetRootInCmd->SetParameterName("inputfile",false);
  SetRootInCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetRoot2DCmd = new G4UIcmdWithAString("/RTPC/Gen/Hist2D",this);
  SetRoot2DCmd->SetGuidance("Set 2D histogram name for sampling");
  SetRoot2DCmd->SetParameterName("TH2D-name",false);
  SetRoot2DCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetNTrackCmd = new G4UIcmdWithAnInteger("/RTPC/Gen/NToBeTracked",this);
  SetNTrackCmd->SetGuidance("Set the number of generated particles to be tracked");
  SetNTrackCmd->SetParameterName("Ntrack",false);
  SetNTrackCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetTrackCmd = new G4UIcmdWithAnInteger("/RTPC/Gen/Track",this);
  SetTrackCmd->SetGuidance("Set the index of the particles to be tracked, this comes from the ntuple branch name.");
  SetTrackCmd->SetParameterName("Track",false);
  SetTrackCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetModeCmd = new G4UIcmdWithAnInteger("/RTPC/Gen/Mode",this);
  SetModeCmd->SetGuidance("Set the mode of the generator, command line, phase space or ROOT");
  SetModeCmd->SetParameterName("Mode",false);
  SetModeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetWindowCmd = new G4UIcmdWithAnInteger("/RTPC/Gen/Window",this);
  SetWindowCmd->SetGuidance("Set 1 to generate target window events, else 0");
  SetWindowCmd->SetParameterName("Window",false);
  SetWindowCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetParticleCmd = new G4UIcmdWith3Vector("/RTPC/Gen/SetParticle",this);
  SetParticleCmd->SetGuidance("Input PDG code, min and max kinetic energy");
  SetParticleCmd->SetParameterName("PDGcode","TMin","TMax",false,false);
  SetParticleCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  SetPDirCmd = new G4UIcmdWith3Vector("/RTPC/Gen/SetPDir",this);
  SetPDirCmd->SetGuidance("Input particle direction (3-vector)");
  SetPDirCmd->SetParameterName("Px","Py","Pz",false,false);
  SetPDirCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete gunDir;
  delete SetRootInCmd;
  delete SetRoot2DCmd;
  delete SetNTrackCmd;
  delete SetTrackCmd;
  delete SetModeCmd;
  delete SetWindowCmd;
  delete SetParticleCmd;
}

void PrimaryGeneratorMessenger::SetNewValue(
                                        G4UIcommand* command, G4String newValue)
{ 
  if( command == SetRootInCmd )
    { Action->SetUpROOTInput(newValue);}
  if( command == SetRoot2DCmd )
    { Action->Set2Dname(newValue);}
  if( command == SetNTrackCmd )
    { Action->SetNTracked(SetNTrackCmd->GetNewIntValue(newValue));}
  if( command == SetTrackCmd )
    { Action->SetTracked(SetTrackCmd->GetNewIntValue(newValue));}
  if( command == SetModeCmd )
    { Action->SetMode(SetModeCmd->GetNewIntValue(newValue));}
  if( command == SetWindowCmd )
    { Action->SetIsWindow(SetWindowCmd->GetNewIntValue(newValue));}
  if( command == SetParticleCmd )
    { Action->SetParticle(SetParticleCmd->GetNew3VectorValue(newValue)); }
  if( command == SetPDirCmd )
    { Action->SetPDir(SetPDirCmd->GetNew3VectorValue(newValue)); }
}



