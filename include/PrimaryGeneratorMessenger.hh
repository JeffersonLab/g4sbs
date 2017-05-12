// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class PrimaryGeneratorMessenger
// Online control of particle gun
// 20/05/13 JRMA adapted from SBS equivalent, under construction

#ifndef PrimaryGeneratorMessenger_h
#define PrimaryGeneratorMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWith3Vector;
class PrimaryGeneratorMessenger: public G4UImessenger
{
public:
  PrimaryGeneratorMessenger(PrimaryGeneratorAction*);
  ~PrimaryGeneratorMessenger();
  void SetNewValue(G4UIcommand*, G4String);  
private:
  PrimaryGeneratorAction* Action;
  G4UIdirectory*  gunDir;
  G4UIcmdWithAString* SetRootInCmd;
  G4UIcmdWithAString* SetRoot2DCmd;
  G4UIcmdWithAnInteger* SetNTrackCmd;
  G4UIcmdWithAnInteger* SetTrackCmd;
  G4UIcmdWithAnInteger* SetModeCmd;
  G4UIcmdWithAnInteger* SetWindowCmd;
  G4UIcmdWith3Vector* SetParticleCmd;
  G4UIcmdWith3Vector* SetPDirCmd;
};

#endif

