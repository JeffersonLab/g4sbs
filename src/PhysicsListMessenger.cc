// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class ESSNPhysicsListMessenger
// Online control of physics used to describe particle interactions in matter
// 20/05/13 JRMA adapted from SBS equivalent, under construction

#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
//#include "PolHadronElasticPhysicsN.hh"
//#include "PolHadronInelasticPhysics.hh"
//#include "PolHadronicProcess.hh"
//#include "PolWHadronElasticProcess.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIdirectory.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
:fpPhysicsList(pPhys)
{   
  fphysDir = new G4UIdirectory("/RTPC/physics/");
  fphysDir->SetGuidance("physics control");
  
  fgammaCutCmd = new G4UIcmdWithADoubleAndUnit("/RTPC/physics/CutGamma",this);  
  fgammaCutCmd->SetGuidance("Set gamma cut.");
  fgammaCutCmd->SetParameterName("Gcut",false);
  fgammaCutCmd->SetUnitCategory("Length");
  fgammaCutCmd->SetRange("Gcut>0.0");
  fgammaCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  felectCutCmd = new G4UIcmdWithADoubleAndUnit("/RTPC/physics/CutEl",this);  
  felectCutCmd->SetGuidance("Set electron cut.");
  felectCutCmd->SetParameterName("Ecut",false);
  felectCutCmd->SetUnitCategory("Length");
  felectCutCmd->SetRange("Ecut>0.0");
  felectCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fposCutCmd = new G4UIcmdWithADoubleAndUnit("/RTPC/physics/CutPos",this);
  fposCutCmd->SetGuidance("Set positron cut.");
  fposCutCmd->SetParameterName("Pcut",false);
  fposCutCmd->SetUnitCategory("Length");
  fposCutCmd->SetRange("Pcut>0.0");
  fposCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fallCutCmd = new G4UIcmdWithADoubleAndUnit("/RTPC/physics/CutsAll",this);
  fallCutCmd->SetGuidance("Set cut for all.");
  fallCutCmd->SetParameterName("cut",false);
  fallCutCmd->SetUnitCategory("Length");
  fallCutCmd->SetRange("cut>0.0");
  fallCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fpListCmd = new G4UIcmdWithAString("/RTPC/physics/Physics",this);
  fpListCmd->SetGuidance("Add modula physics list.");
  fpListCmd->SetParameterName("PList",false);
  fpListCmd->AvailableForStates(G4State_PreInit);

  flistCmd = new G4UIcmdWithoutParameter("/RTPC/physics/ListPhysics",this);
  flistCmd->SetGuidance("Available Physics Lists");
  flistCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fSpinPrecN = new G4UIcmdWith3Vector("/RTPC/pol/N-spin-precess",this);
  fSpinPrecN->AvailableForStates(G4State_PreInit,G4State_Idle);
  fSpinPrecN->SetGuidance("Enable nucleon spin precession in a magnetic field");
  fSpinPrecN->SetGuidance("Set proton and neutron anomalous magnetic moment");
  fSpinPrecN->SetParameterName("muP","muN","spare",true,true); 

  fNproc = NULL;
  fPolRot = NULL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete fgammaCutCmd;
  delete felectCutCmd;
  delete fposCutCmd;
  delete fallCutCmd;
  delete fpListCmd;
  delete flistCmd;
  delete fphysDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  //  PolHadronicProcess* proc;
  if( command == fgammaCutCmd )
    fpPhysicsList->SetCutForGamma(fgammaCutCmd->GetNewDoubleValue(newValue));

  else if( command == felectCutCmd )
    fpPhysicsList->SetCutForElectron(felectCutCmd->GetNewDoubleValue(newValue));

  else if( command == fposCutCmd )
   fpPhysicsList->SetCutForPositron(fposCutCmd->GetNewDoubleValue(newValue));

  else if( command == fallCutCmd )
    {
      G4double cut = fallCutCmd->GetNewDoubleValue(newValue);
      fpPhysicsList->SetCutForGamma(cut);
      fpPhysicsList->SetCutForElectron(cut);
      fpPhysicsList->SetCutForPositron(cut);
    }
  // Add physics process
  else if( command == fpListCmd ) {
    G4String name = newValue;
    if(name == "PHYSLIST") {
      char* path = getenv(name);
      if (path) name = G4String(path);
      else {
        G4cout << "### PhysicsListMessenger WARNING: "
	       << " environment variable PHYSLIST is not defined"
	       << G4endl;
	return; 
      }
    }
    fpPhysicsList->AddPhysicsList(name);
  }
  // List physics processes
  else if( command == flistCmd )
    fpPhysicsList->List();
  // Set Hadronic process for any subsequent polarisation control cmds
  // Default analysing-power value
  // Analysing power model
  else if( command == fSpinPrecN ) {
    G4ThreeVector mu = fSpinPrecN->GetNew3VectorValue(newValue);
    fpPhysicsList->SetSpinPrecN(mu.x(),mu.y());
  }    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
