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

#include "TString.h"

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
  // if( argc == 2 ){ //if only one command line argument is given, assume all commands are to be applied post-initialization:
  //   postinit_macro = argv[1];
  // } else if( argc > 2 ){ //assume first argument is pre-init commands, second argument is post-init commands:
  //   preinit_macro = argv[1];
  //   postinit_macro = argv[2];
  // }

  if( argc > 1 ){ //command-line arguments were given: preserve previous behavior while adding new dials to turn:
    TString argument;
    G4int nmacros = 0; //ignore more than two macro definitions:
    G4bool force_batch = false;

    TString macrodefs[2];
    
    for( int iarg=1; iarg<argc; iarg++ ){
      argument = argv[iarg];

      if( argument.BeginsWith("batch=") ){
	argument.ReplaceAll("batch=","");
	std::istringstream is( argument.Data() );
	force_batch = true;
	if( argument.Contains("true")||argument.Contains("false") ){ //true/false
	  is >> std::boolalpha >> batch_mode;
	} else { //1/0
	  is >> batch_mode;
	}
      }
      
      if( argument.EndsWith(".mac") && nmacros < 2 ){ //macro file:
	//Parse command line arguments:
	macrodefs[nmacros] = argument;
	nmacros++;
      }
    }

    if( nmacros == 0 ){
      batch_mode = false; //no macro --> meaningless to run in batch!
    } else if( nmacros == 1 ){
      if( !macrodefs[0].BeginsWith("preinit=") ){
      //only one macro definition was given and it was not forced to pre-init:
      //default to post-init execution:
	macrodefs[0].ReplaceAll("postinit=","");
	postinit_macro = macrodefs[0].Data();
	if( !force_batch ) batch_mode = true; //default to batch mode unless explicitly overridden by user
      } else {
	macrodefs[0].ReplaceAll("preinit=","");
	preinit_macro = macrodefs[0].Data();
	batch_mode = false; 
      }
    } 
    if( nmacros == 2 ){ //assume first macrodef is preinit and second is post-init unless
      //explicitly overridden by user:
      TString m1 = macrodefs[0];
      TString m2 = macrodefs[1];

      TString mac1( m1(m1.Index("=")+1,m1.Length()-m1.Index("=")-1) );
      TString mac2( m2(m2.Index("=")+1,m2.Length()-m2.Index("=")-1) );
      
      //under what conditions is m1 the postinit and m2 the pre-init?
      // 1) m1 IS labeled as postinit and m2 is NOT explicitly labeled as postinit
      // 2) m2 is explicitly labeled as preinit and m1 is NOT explicitly labeled as preinit
      // These are the only cases that can override the default behavior of m1 being pre-init and
      // m2 being post-init:
      if( (m1.BeginsWith("postinit=") && !m2.BeginsWith("postinit=") ) ||
	  (m2.BeginsWith("preinit=") && !m1.BeginsWith("preinit=") ) ){
	postinit_macro = mac1.Data();
	preinit_macro = mac2.Data();
      } else {
	preinit_macro = mac1.Data();
	postinit_macro = mac2.Data();
      }
      //Run in batch mode unless explicitly overridden by user:
      if( !force_batch ) batch_mode = true;
    }
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
    G4String command = "/control/execute ";
    command += preinit_macro;
    UImanager->ApplyCommand(command);
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

  G4cout << "batch mode = " << batch_mode << G4endl;
  
  
  if( !batch_mode )
  {
    //--------------------------
    // Define (G)UI
    //--------------------------

    
    
#ifdef G4UI_USE
    G4UIExecutive * ui = new G4UIExecutive(argc,argv);
    
    UImanager->SetSession( ui->GetSession() ); 
    
    if( postinit_macro != "" ){
      rundata->SetMacroFile(postinit_macro);
      G4String command = "/control/execute ";
      command += postinit_macro;
      UImanager->ApplyCommand(command);
    }

    ui->SessionStart();
			  
    delete ui;
#endif
  } else { //batch mode:
      
    rundata->SetMacroFile(postinit_macro);
    G4String command = "/control/execute ";
    G4String fileName = postinit_macro;
    UImanager->ApplyCommand(command+fileName);
    
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
