#include "TString.h"  //  These need to come first or the
#include "TBuffer.h"  //  compiler bitches about shadowed variables
#include "TMatrixTBase.h"  
#include "THashTable.h"  
#include "CLHEP/Random/Random.h"

#include "G4SBSRunAction.hh"
#include "G4SBSPrimaryGeneratorAction.hh"
#include "G4SBSEventAction.hh"
#include "G4SBSSteppingAction.hh"
#include "G4SBSRun.hh"
#include "G4SBSRunData.hh"


//------------
// Geometries:
//------------
#include "G4SBSDetectorConstruction.hh"

#include "G4SBSIO.hh"
#include "G4SBSMessenger.hh"

//-----------------------------------
// G4SBSPhysicsList (makes use of the
// G4ParameterisationManagerProcess):
//-----------------------------------
#include "G4SBSPhysicsList.hh"

#include "G4UImanager.hh"
#include "G4RunManager.hh"
#include "G4RunManagerKernel.hh"

#include "G4UnitsTable.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#ifdef __APPLE__
#include "unistd.h"
#endif

bool parseArgument(std::string arg, std::string &name, std::string &value)
{
  size_t length = arg.length();

  // There are two possible types of parameters.
  // 1) Parameters with a value (--name=value)
  // 2) Parameters with no value (--name), also known as flags.

  // Minimum size for the argument is 3, two dashes and at least one
  // other character
  if(length<3)
    return false; // Not long enough to be a valid parameter

  // Now, look for the double dash at the beginning that will signal
  // that a parameter is being specified
  size_t pos1 = arg.find_first_of("--");
  if(pos1 != std::string::npos && pos1 == 0) { // Double dash found

    // First, identify if it is parameter type 1 by looking for the equal sign
    size_t pos2 = arg.find_first_of("=");
    if(pos2 != std::string::npos && pos2 < length-1 && pos2>0) { // Param type 1
      // Split name from value using on both sides of the equal sign
      name = arg.substr(2,pos2-2);
      value = arg.substr(pos2+1,length);
      return true;
    } else { // No equal sign means parameter type 2
      name = arg.substr(2,length);
      return true;
    }
  }
  return false;
}

void executeMacro(G4String macro, G4UImanager *UImanager)
{
  if(macro.length()>0) {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  }
}

bool getBool(std::string value, bool default_value)
{
  if(value.compare("true") == 0) {
    return true;
  } else if(value.compare("false") == 0) {
    return false;
  } else {
    return default_value;
  }
}

int main(int argc, char** argv)
{
  //Allow user to apply some (perhaps all?) commands in the pre-initialization phase:

  //-------------------------------
  // Initialize the CLHEP random engine used by
  // "shoot" type functions
  //-------------------------------

  G4bool batch_mode = false; //Run in interactive mode by default
  
  G4String preinit_macro = "";
  G4String postinit_macro = "";
  bool flag_gui = false; // Display GUI flag

  //-------------------------------
  // Parse command line arguments.
  //
  // In order to maintain backwards compatibility with the old method
  // let's keep track of a flag. If any of the new paramethers are found
  // we will require the new way. If none are found, we'll proceed with
  // the old method.
  //-------------------------------
  bool paramsFound = false;
  for(int i = 1; i < argc; i++) {
    std::string paramName;
    std::string paramValue;
    bool validParam = parseArgument(argv[i],paramName,paramValue);
    if(validParam) {
      paramsFound = true;
      if(paramName.compare("pre") == 0) {
        preinit_macro = paramValue;
      } else if (paramName.compare("post") == 0) {
        postinit_macro = paramValue;
      } else if (paramName.compare("gui") == 0) {
        flag_gui = getBool(paramValue,true);
      }
    }
  }

  // If no valid parameters were found, revert to the old method, in which
  // there are two possible parameters to specify the pre-init macro and
  // post init-macro.
  if(!paramsFound) {
    if( argc == 2 ){ //if only one command line argument is given, assume all commands are to be applied post-initialization:
      postinit_macro = argv[1];
    } else if( argc > 2 ){ //assume first argument is pre-init commands, second argument is post-init commands:
      preinit_macro = argv[1];
      postinit_macro = argv[2];
    }
  }

  bool use_gui = false;
  // If no postinit macro specified, turn on the GUI
  if(postinit_macro == "") {
    use_gui = true;
  }
  // If the GUI flag/parameter is passed, always display GUI
  if(flag_gui) {
    use_gui = true;
  }

  CLHEP::HepRandom::createInstance();

  unsigned int seed = time(0) + (int) getpid();
  unsigned int devrandseed = 0;
  FILE *fdrand = fopen("/dev/urandom", "r");
  if( fdrand ){
    fread(&devrandseed, sizeof(int), 1, fdrand);
    seed += devrandseed;
    fclose(fdrand);
  }

  CLHEP::HepRandom::setTheSeed(seed);

  G4SBSRunData *rundata = G4SBSRun::GetRun()->GetData();

  rundata->SetSeed(seed);

  G4SBSIO *io = new G4SBSIO();



  //-------------------------------
  // Initialization of Run manager
  //-------------------------------
  G4cout << "RunManager construction starting...." << G4endl;
  G4RunManager * runManager = new G4RunManager;

  G4SBSMessenger *sbsmess = new G4SBSMessenger();
  sbsmess->SetIO(io);

  G4VModularPhysicsList *physicslist = new G4SBSPhysicsList; 
  runManager->SetUserInitialization(physicslist);

  sbsmess->SetPhysList( (G4SBSPhysicsList*) physicslist );

  // Detector/mass geometry and parallel geometry(ies):
  G4SBSDetectorConstruction* detector = new G4SBSDetectorConstruction();
  sbsmess->SetDetCon((G4SBSDetectorConstruction *) detector);
  // This lets us output graphical field maps
  io->SetGlobalField( ((G4SBSDetectorConstruction *) detector)->GetGlobalField() );

  runManager->SetUserInitialization(detector);
  
  

  //-------------------------------
  // UserAction classes
  //-------------------------------
  G4UserRunAction* run_action = new G4SBSRunAction;
  ((G4SBSRunAction *) run_action)->SetIO(io);
  runManager->SetUserAction(run_action);
  //
  G4VUserPrimaryGeneratorAction* gen_action = new G4SBSPrimaryGeneratorAction;
  ((G4SBSPrimaryGeneratorAction *) gen_action)->SetIO(io);
  sbsmess->SetEvGen(((G4SBSPrimaryGeneratorAction *) gen_action)->GetEvGen());
  sbsmess->SetPriGen((G4SBSPrimaryGeneratorAction *) gen_action);

  ( (G4SBSPrimaryGeneratorAction*) gen_action )->SetRunAction( (G4SBSRunAction*) run_action );
  
  runManager->SetUserAction(gen_action);
  //
  G4UserEventAction* event_action = new G4SBSEventAction;
  ((G4SBSEventAction *) event_action)->SetIO(io);
  ((G4SBSEventAction *) event_action)->SetEvGen(((G4SBSPrimaryGeneratorAction *) gen_action)->GetEvGen());
  runManager->SetUserAction(event_action);
  sbsmess->SetEvAct((G4SBSEventAction *) event_action);
  //
  G4UserSteppingAction* stepping_action = new G4SBSSteppingAction;
  runManager->SetUserAction(stepping_action);

  G4UImanager * UImanager = G4UImanager::GetUIpointer();

  if( preinit_macro != "" ){
    rundata->SetPreInitMacroFile(preinit_macro);
    executeMacro(preinit_macro,UImanager);
  }

  if( postinit_macro != "") {
    rundata->SetMacroFile(postinit_macro);
  }

  // Initialize Run manager
  runManager->Initialize();
 
  //----------------
  // Visualization:
  //----------------
#ifdef G4VIS_USE
  G4cout << "Instantiating Visualization Manager......." << G4endl;
  G4VisManager* visManager = new G4VisExecutive;
  visManager -> Initialize ();
#endif

  // Setup commands
  //
 
  /*
  UImanager->ApplyCommand("/Step/Verbose 0");
  UImanager->ApplyCommand("/tracking/Verbose 1");
  UImanager->ApplyCommand("/gun/particle e-");
  UImanager->ApplyCommand("/gun/energy 100 MeV");
  UImanager->ApplyCommand("/gun/direction 0 0 1");
  UImanager->ApplyCommand("/gun/position 0 0 0");
  UImanager->ApplyCommand("/gun/direction 0 .3 1.");
  */

  if( use_gui ) {
    //--------------------------
    // Define (G)UI
    //--------------------------

#ifdef G4UI_USE
    G4UIExecutive * ui = new G4UIExecutive(argc,argv);
    UImanager->SetSession( ui->GetSession() ); 
#endif

    executeMacro(postinit_macro, UImanager);

#ifdef G4UI_USE
    ui->SessionStart();

    delete ui;
#endif
  } else {
    // Run the postinit macro if one is specified
    executeMacro(postinit_macro, UImanager);
  }

  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !

#ifdef G4VIS_USE
//  delete visManager;
#endif
//  delete runManager;

  return 0;
}
