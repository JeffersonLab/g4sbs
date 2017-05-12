// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class PhysicsListMessenger
// Online control of physics used to describe particle interactions in matter
// 20/05/13 JRMA adapted from SBS equivalent, under construction

#ifndef PhysicsListMessenger_h
#define PhysicsListMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class PhysicsList;
class G4UIcmdWithADouble;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3Vector;
class G4UIdirectory;
class PolHadronicProcess;
class G4HadronicProcess;
class PolNucleonRotate;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsListMessenger: public G4UImessenger
{
public:
  
  PhysicsListMessenger(PhysicsList* );
  virtual ~PhysicsListMessenger();
    
  void SetNewValue(G4UIcommand*, G4String);
    
private:
  
  PhysicsList* fpPhysicsList;
  // physics cuts and listing
  G4UIcmdWithADoubleAndUnit* fgammaCutCmd;
  G4UIcmdWithADoubleAndUnit* felectCutCmd;
  G4UIcmdWithADoubleAndUnit* fposCutCmd;
  G4UIcmdWithADoubleAndUnit* fallCutCmd;
  G4UIcmdWithAString*        fpListCmd;
  G4UIcmdWithoutParameter*   flistCmd;  
  G4UIdirectory* fphysDir;
  // polarised nucleon stuff
  G4UIcmdWith3Vector* fSpinPrecN;  
  PolHadronicProcess* fNproc;
  PolNucleonRotate* fPolRot;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

