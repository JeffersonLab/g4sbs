#include "TBuffer.h"
#include "TString.h"
#include "TMatrixTBase.h"
#include "THashTable.h"
#include "G4SBSMessenger.hh"
#include "G4SBSRun.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWith3Vector.hh"

#include "G4SBSDetectorConstruction.hh"
#include "G4SBSIO.hh"
#include "G4SBSEventGen.hh"
#include "G4SBSEventAction.hh"

#include "G4MagneticField.hh"
#include "G4SBSGlobalField.hh"
#include "G4SBSMagneticField.hh"
#include "G4SBSToscaField.hh"
#include "G4SBSBigBiteField.hh"
#include "G4SBSPrimaryGeneratorAction.hh"
#include "G4SBSPhysicsList.hh"
#include "G4OpticalPhysics.hh"
#include "G4ParticleGun.hh"

#include "G4SBSBeamlineBuilder.hh"
#include "G4SBSTargetBuilder.hh"
#include "G4SBSEArmBuilder.hh"
#include "G4SBSHArmBuilder.hh"

#include "G4SBSSteppingAction.hh"
#include "G4SBSTrackingAction.hh"

#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4UImanager.hh"
#include "G4RunManager.hh"
#ifdef G4SBS_USE_GDML
#include "G4GDMLParser.hh"
#endif
#include "G4VPhysicalVolume.hh"

#ifdef __APPLE__
#include <unistd.h>
#endif

using namespace CLHEP;

G4SBSMessenger::G4SBSMessenger(){
  fExpType = G4SBS::kGMN; //default to GMN

  runCmd = new G4UIcmdWithAnInteger("/g4sbs/run",this);
  runCmd->SetGuidance("Run simulation with x events");
  runCmd->SetParameterName("nevt", false);

  printCmd = new G4UIcmdWithAnInteger("/g4sbs/print",this); 
  printCmd->SetGuidance("Print the line number (arg = number)"); 
  printCmd->SetParameterName("print",false); 

  gemconfigCmd = new G4UIcmdWithAnInteger("/g4sbs/gemconfig",this);
  gemconfigCmd->SetGuidance("BigBite GEM layout: option 1 (default), 2 or 3");
  gemconfigCmd->SetParameterName("gemconfig", false);

  shieldconfigCmd = new G4UIcmdWithAnInteger("/g4sbs/bbshieldconfig",this);
  shieldconfigCmd->SetGuidance("BB Ecal shielding layout: option 0 (none), 1 (default), 2 (+10cm Al + 3cm steel on side), 3 (+3cm steel + 3cm steel on side)");
  shieldconfigCmd->SetParameterName("bbshieldconfig", false);

  bbpsconfigCmd = new G4UIcmdWithAnInteger("/g4sbs/bbpsconfig",this);
  bbpsconfigCmd->SetGuidance("BB PS config: option 0 (old geometry), 1 (new modules, 25 blocks), 2 (new modules, 26 blocks)");
  bbpsconfigCmd->SetParameterName("bbpsconfig", false);
  
  CDetconfigCmd = new G4UIcmdWithAnInteger("/g4sbs/CDetconfig",this);
  CDetconfigCmd->SetGuidance("CDet Geometry Options(integer 1,2,3, or 4): Option 1 (default) = Simplest option, only material budget with no SD assigned, Option 2 = Flat, 2 planes with sensitive regions, Option 3 = Top/Bottom 'modules' are angled relative to central 'module', and Option 4 = Each bar is shimmed in order to optimize normal incidence");
  CDetconfigCmd->SetParameterName("CDetconfig",false);

  flipGEMCmd = new G4UIcmdWithABool("/g4sbs/flipGEM",this);
  flipGEMCmd->SetGuidance("Reverse GEM orientation front-to-back (bool). Applies to ALL GEMs!");
  flipGEMCmd->SetParameterName("flipGEM", false );
  
  ECALmapfileCmd = new G4UIcmdWithAString("/g4sbs/ECALmap",this);
  ECALmapfileCmd->SetGuidance("Name of text file listing active ECAL cells (assumed to be located in database/)");
  ECALmapfileCmd->SetParameterName("ECALmapfile",false);

  fileCmd = new G4UIcmdWithAString("/g4sbs/filename",this);
  fileCmd->SetGuidance("Output ROOT filename");
  fileCmd->SetParameterName("filename", false);

  sigfileCmd = new G4UIcmdWithAString("/g4sbs/sigmafile",this);
  sigfileCmd->SetGuidance("File containing GEM coordinate resolutions by chamber ID number");
  sigfileCmd->SetParameterName("sigmafile", false);

  tgtCmd = new G4UIcmdWithAString("/g4sbs/target",this);
  tgtCmd->SetGuidance("Target type from LH2, LD2, H2, D2, 3He, (fictional) neutron target");
  tgtCmd->SetParameterName("targtype", false);

  // D. Flay (7/28/20) 
  // for GEn 3He target Helmholtz coil configuration 
  GENTargetHelmholtzCmd = new G4UIcmdWithAnInteger("/g4sbs/targgenhhconfig",this);
  GENTargetHelmholtzCmd->SetGuidance("GEn 3He target Helmholts coil configuration based on central Q2 value"); 
  GENTargetHelmholtzCmd->SetGuidance("146 => Q2 = 1.46 (GeV/c)^2, 368 => Q2 = 3.68 (GeV/c)^2, 677 => Q2 = 6.77 (GeV/c)^2, 1018 => Q2 = 10.18 (GeV/c)^2, "); 
  GENTargetHelmholtzCmd->SetParameterName("targgenhhconfig",false); // user must provide an integer value, non-argument not allowed 
  GENTargetHelmholtzCmd->SetDefaultValue(G4SBS::kGEN_146);          // probably not utilized since we require an input value 

  // D. Flay (9/29/20) 
  // for GEn 3He target rotational misalignment 
  GENTargetRXCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/targgenDRX",this); 
  GENTargetRXCmd->SetGuidance("GEn 3He target rotational misalignment relative to x axis"); 
  GENTargetRXCmd->SetParameterName("targgenDRX",true); // second argument = omittable?  
  GENTargetRYCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/targgenDRY",this); 
  GENTargetRYCmd->SetGuidance("GEn 3He target rotational misalignment relative to y axis"); 
  GENTargetRYCmd->SetParameterName("targgenDRY",true); 
  GENTargetRZCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/targgenDRZ",this); 
  GENTargetRZCmd->SetGuidance("GEn 3He target rotational misalignment relative to z axis"); 
  GENTargetRZCmd->SetParameterName("targgenDRZ",true);

  // D. Flay (10/9/20) 
  // for GEn 3He target collimators
  // all: if enabled, *allows* the collimators to be built.  can toggle on/off individual ones, see below  
  GENTargetColCmd = new G4UIcmdWithABool("/g4sbs/targgenColEnable",this); 
  GENTargetColCmd->SetGuidance("GEn 3He target collimator enable.  If enabled, allows collimators to be built"); 
  GENTargetColCmd->SetParameterName("targgenColEnable",true);
  // A (upstream) 
  GENTargetColACmd = new G4UIcmdWithABool("/g4sbs/targgenColEnableA",this); 
  GENTargetColACmd->SetGuidance("GEn 3He target collimator A enable"); 
  GENTargetColACmd->SetParameterName("targgenColEnableA",true); 
  // B (downstream, 1st) 
  GENTargetColBCmd = new G4UIcmdWithABool("/g4sbs/targgenColEnableB",this); 
  GENTargetColBCmd->SetGuidance("GEn 3He target collimator B enable"); 
  GENTargetColBCmd->SetParameterName("targgenColEnableB",true); 
  // C (downstream, 2nd, furthest from target) 
  GENTargetColCCmd = new G4UIcmdWithABool("/g4sbs/targgenColEnableC",this); 
  GENTargetColCCmd->SetGuidance("GEn 3He target collimator C enable"); 
  GENTargetColCCmd->SetParameterName("targgenColEnableC",true);

  // D. Flay (12/9/20) 
  // for enabling the GEn target as a sensitive detector  
  GENTargetSDEnableCmd = new G4UIcmdWithABool("/g4sbs/targgenSDEnable",this); 
  GENTargetSDEnableCmd->SetGuidance("GEn 3He target SD enable"); 
  GENTargetSDEnableCmd->SetParameterName("targgenSDEnable",true);

  // D. Flay (10/15/20) 
  // beam angular misalignment 
  // - horizontal (x)  
  beamAngleXcmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamAngleX",this);
  beamAngleXcmd->SetGuidance("Rotate beam momentum about the x axis");
  beamAngleXcmd->SetParameterName("beamAngleX",true);  // this is omittable 
  // - vertical (y)
  beamAngleYcmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamAngleY",this);
  beamAngleYcmd->SetGuidance("Rotate beam momentum about the y axis");
  beamAngleYcmd->SetParameterName("beamAngleY",true);   
  // - axial (z)
  beamAngleZcmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamAngleZ",this);
  beamAngleZcmd->SetGuidance("Rotate beam momentum about the z axis");
  beamAngleZcmd->SetParameterName("beamAngleZ",true);  // must provide input

  // D. Flay (10/15/20) 
  // ruidmentary ion chamber 
  ionChamberEnableCmd = new G4UIcmdWithABool("/g4sbs/ionChamberEnable",this);
  ionChamberEnableCmd->SetGuidance("Enable an ion chamber, including sensitive detector capability");
  ionChamberEnableCmd->SetParameterName("ionChamberEnable",true); // this is omittable 
  // coordinates
  // -x 
  ionChamberXCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/ionChamberX",this);
  ionChamberXCmd->SetGuidance("Ion chamber x coordinate");
  ionChamberXCmd->SetParameterName("ionChamberX",true);  
  // - y 
  ionChamberYCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/ionChamberY",this);
  ionChamberYCmd->SetGuidance("Ion chamber y coordinate");
  ionChamberYCmd->SetParameterName("ionChamberY",true);   
  // - z 
  ionChamberZCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/ionChamberZ",this);
  ionChamberZCmd->SetGuidance("Ion chamber z coordinate");
  ionChamberZCmd->SetParameterName("ionChamberZ",true); 
  // rotation 
  // -x 
  ionChamberRXCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/ionChamberRX",this);
  ionChamberRXCmd->SetGuidance("Ion chamber angle about x");
  ionChamberRXCmd->SetParameterName("ionChamberRX",true);  
  // - y 
  ionChamberRYCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/ionChamberRY",this);
  ionChamberRYCmd->SetGuidance("Ion chamber angle about y");
  ionChamberRYCmd->SetParameterName("ionChamberRY",true);   
  // - z 
  ionChamberRZCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/ionChamberRZ",this);
  ionChamberRZCmd->SetGuidance("Ion chamber angle about z");
  ionChamberRZCmd->SetParameterName("ionChamberRZ",true);   

  // D. Flay (11/5/20) 
  // ruidmentary beam collimator for the GEn target 
  // downstream
  // enable  
  beamCollimatorEnableDnCmd = new G4UIcmdWithABool("/g4sbs/beamCollimatorEnable_dnstr",this);
  beamCollimatorEnableDnCmd->SetGuidance("Enable a beam collimator for the GEn target");
  beamCollimatorEnableDnCmd->SetParameterName("beamCollimatorEnable_dnstr",true); // must provide input 
  // geometry  
  // -length 
  beamCollimatorLDnCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamCollimatorL_dnstr",this);
  beamCollimatorLDnCmd->SetGuidance("Beam collimator length");
  beamCollimatorLDnCmd->SetParameterName("beamCollimatorL_dnstr",true);  
  // - y 
  beamCollimatorDminDnCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamCollimatorDmin_dnstr",this);
  beamCollimatorDminDnCmd->SetGuidance("Beam collimator diameter (inner)");
  beamCollimatorDminDnCmd->SetParameterName("beamCollimatorDmin_dnstr",true);   
  // - z 
  beamCollimatorDmaxDnCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamCollimatorDmax_dnstr",this);
  beamCollimatorDmaxDnCmd->SetGuidance("Beam collimator diameter (outer)");
  beamCollimatorDmaxDnCmd->SetParameterName("beamCollimatorDmax_dnstr",true);   
  // coordinates
  // -x 
  beamCollimatorXDnCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamCollimatorX_dnstr",this);
  beamCollimatorXDnCmd->SetGuidance("Beam collimator coordinate");
  beamCollimatorXDnCmd->SetParameterName("beamCollimatorX_dnstr",true);  
  // - y 
  beamCollimatorYDnCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamCollimatorY_dnstr",this);
  beamCollimatorYDnCmd->SetGuidance("Beam collimator y coordinate");
  beamCollimatorYDnCmd->SetParameterName("beamCollimatorY_dnstr",true);   
  // - z 
  beamCollimatorZDnCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamCollimatorZ_dnstr",this);
  beamCollimatorZDnCmd->SetGuidance("Beam collimator z coordinate");
  beamCollimatorZDnCmd->SetParameterName("beamCollimatorZ_dnstr",true);   
  // upstream
  // enable  
  beamCollimatorEnableUpCmd = new G4UIcmdWithABool("/g4sbs/beamCollimatorEnable_upstr",this);
  beamCollimatorEnableUpCmd->SetGuidance("Enable a beam collimator for the GEn target");
  beamCollimatorEnableUpCmd->SetParameterName("beamCollimatorEnable_upstr",true); 
  // geometry  
  // -length 
  beamCollimatorLUpCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamCollimatorL_upstr",this);
  beamCollimatorLUpCmd->SetGuidance("Beam collimator length");
  beamCollimatorLUpCmd->SetParameterName("beamCollimatorL_upstr",true);  
  // - y 
  beamCollimatorDminUpCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamCollimatorDmin_upstr",this);
  beamCollimatorDminUpCmd->SetGuidance("Beam collimator diameter (inner)");
  beamCollimatorDminUpCmd->SetParameterName("beamCollimatorDmin_upstr",true);   
  // - z 
  beamCollimatorDmaxUpCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamCollimatorDmax_upstr",this);
  beamCollimatorDmaxUpCmd->SetGuidance("Beam collimator diameter (outer)");
  beamCollimatorDmaxUpCmd->SetParameterName("beamCollimatorDmax_upstr",true);   
  // coordinates
  // -x 
  beamCollimatorXUpCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamCollimatorX_upstr",this);
  beamCollimatorXUpCmd->SetGuidance("Beam collimator coordinate");
  beamCollimatorXUpCmd->SetParameterName("beamCollimatorX_upstr",true);  
  // - y 
  beamCollimatorYUpCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamCollimatorY_upstr",this);
  beamCollimatorYUpCmd->SetGuidance("Beam collimator y coordinate");
  beamCollimatorYUpCmd->SetParameterName("beamCollimatorY_upstr",true);   
  // - z 
  beamCollimatorZUpCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamCollimatorZ_upstr",this);
  beamCollimatorZUpCmd->SetGuidance("Beam collimator z coordinate");
  beamCollimatorZUpCmd->SetParameterName("beamCollimatorZ_upstr",true);   

  kineCmd = new G4UIcmdWithAString("/g4sbs/kine",this);
  kineCmd->SetGuidance("Kinematics from elastic, inelastic, flat, dis, beam, sidis, wiser, gun, pythia6, simc, wapp");
  kineCmd->SetParameterName("kinetype", false);

  PYTHIAfileCmd = new G4UIcmdWithAString("/g4sbs/pythia6file",this);
  PYTHIAfileCmd->SetGuidance("Name of ROOT file containing PYTHIA6 events as a ROOT tree");
  PYTHIAfileCmd->SetParameterName("fname",false);
  
  SIMCfileCmd = new G4UIcmdWithAString("/g4sbs/simcfile",this);
  SIMCfileCmd->SetGuidance("Name of ROOT file containing SIMC events as a ROOT tree");
  SIMCfileCmd->SetParameterName("fname",false);
  
  expCmd = new G4UIcmdWithAString("/g4sbs/exp",this);
  expCmd->SetGuidance("Experiment type from gep, gmn, gen, a1n, sidis, C16, tdis, ndvcs, genrp");
  expCmd->SetParameterName("exptype", false);

  GunParticleCmd = new G4UIcmdWithAString("/g4sbs/particle",this);
  GunParticleCmd->SetGuidance("Particle type for gun generator (valid GEANT4 particle names)");
  GunParticleCmd->SetParameterName("ptype", false );

  HadrCmd = new G4UIcmdWithAString("/g4sbs/hadron",this);
  HadrCmd->SetGuidance("Hadron type h for SIDIS N(e,e'h)X generator: pi+/pi-/K+/K-/p/pbar possible. Also, Nucleon type N for SIMC A(e,e'N) generator: p/n are possible.");
  HadrCmd->SetParameterName("hadrontype", false );

  RejectionSamplingCmd = new G4UIcommand("/g4sbs/rejectionsampling",this);
  RejectionSamplingCmd->SetGuidance("Turn on rejection sampling in event generator");
  RejectionSamplingCmd->SetGuidance("Produces events distributed according to xsec");
  RejectionSamplingCmd->SetGuidance("Applies only to generators: elastic, inelastic, dis, sidis, wiser");
  RejectionSamplingCmd->SetGuidance("Usage: /g4sbs/rejectionsampling flag nevent");
  RejectionSamplingCmd->SetGuidance("flag (bool) = toggle on/off");
  RejectionSamplingCmd->SetGuidance("nevent = Number of \"pre-events\" (not generated) to estimate maximum event weight within generation limits");
  RejectionSamplingCmd->SetParameter( new G4UIparameter("flag",'b',false) );
  RejectionSamplingCmd->SetParameter( new G4UIparameter("N", 'i', true) );
  RejectionSamplingCmd->GetParameter(1)->SetDefaultValue(100000);
				     
  
  bigfieldCmd = new G4UIcmdWithAnInteger("/g4sbs/48d48field", this);
  bigfieldCmd->SetGuidance("0 = turn off SBS constant magnetic field, 1 = turn on SBS constant magnetic field");
  bigfieldCmd->SetParameterName("48d48field", false);

  bbfieldCmd = new G4UIcommand("/g4sbs/bbfield", this);
  bbfieldCmd->SetGuidance("Turn on Bigbite field (requires field map)");
  bbfieldCmd->SetGuidance("usage: /g4sbs/bbfield flag fname");
  bbfieldCmd->SetGuidance("flag = 1/0 for on/off");
  bbfieldCmd->SetGuidance("fname = field map file name (D = map_696A.dat)");
  bbfieldCmd->SetParameter(new G4UIparameter("bbfield",'i', false) );
  bbfieldCmd->SetParameter(new G4UIparameter("fname",'s', true) );
  bbfieldCmd->GetParameter(1)->SetDefaultValue("map_696A.dat");

  // bbfield_fnameCmd = new G4UIcmdWithAString("/g4sbs/bbfieldmapfname",this);
  // bbfield_fnameCmd->SetGuidance("BigBite field map file name (if non-standard name/location)");
  // bbfield_fnameCmd->SetParameterName("bbfieldfname",false);
  
  //tosfieldCmd = new G4UIcmdWithAString("/g4sbs/tosfield", this);
  tosfieldCmd = new G4UIcommand("/g4sbs/tosfield", this );
  tosfieldCmd->SetGuidance("Add TOSCA field map to global field definition");
  tosfieldCmd->SetGuidance("Usage: /g4sbs/tosfield fname flag");
  tosfieldCmd->SetGuidance("fname = field map file name");
  tosfieldCmd->SetGuidance("flag = 0 (default): transformation according to header information in map file");
  tosfieldCmd->SetGuidance("flag = 1 (BigBite): transformation according to BigBite angle/magnet distance" );
  tosfieldCmd->SetGuidance("flag = 2 (SBS): transformation according to SBS angle/magnet distance" );
  tosfieldCmd->SetParameter(new G4UIparameter("tosfield", 's', false) ); //map file name
  tosfieldCmd->SetParameter(new G4UIparameter("fieldflag", 'i', true) );
  tosfieldCmd->GetParameter(1)->SetDefaultValue(0);

  eventStatusEveryCmd = new G4UIcmdWithAnInteger("/g4sbs/eventstatusevery", this);
  eventStatusEveryCmd->SetGuidance("Print event status at every N entries");
  eventStatusEveryCmd->SetParameterName("eventstatusevery", false);

  geantinoCmd = new G4UIcmdWithABool("/g4sbs/shootgeantino", this);
  geantinoCmd->SetGuidance("Shoot a geantino instead of e-");
  geantinoCmd->SetParameterName("shootgeantino", false);

  invertCmd = new G4UIcmdWithABool("/g4sbs/invertfield", this);
  invertCmd->SetGuidance("invert global field polarity (inverts all global fields)");
  invertCmd->SetParameterName("invert", false);

  totalabsCmd = new G4UIcmdWithABool("/g4sbs/totalabs", this);
  totalabsCmd->SetGuidance("Magnet materials are total absorbers (default=false)");
  totalabsCmd->SetParameterName("totalabs", false);

  checkOverlapCmd = new G4UIcmdWithABool("/g4sbs/checkOverlap", this);
  checkOverlapCmd->SetGuidance("Check if volumes overlap (default=false)");
  checkOverlapCmd->SetParameterName("checkOverlap", false);

  tgtLenCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/targlen",this);
  tgtLenCmd->SetGuidance("Target length along beam direction");
  tgtLenCmd->SetParameterName("targlen", false);

  tgtDenCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/targden",this);
  tgtDenCmd->SetGuidance("Target density");
  tgtDenCmd->SetParameterName("targden", false);

  tgtPresCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/targpres",this);
  tgtPresCmd->SetGuidance("Gas Target pressure (applies to H2, 3He)");
  tgtPresCmd->SetParameterName("targpres", false);

  tgtDiamCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/targdiam",this);
  tgtDiamCmd->SetGuidance("Cryotarget or gas target cell diameter (assumed cylindrical), applies to LH2/LD2/H2/3He");
  tgtDiamCmd->SetParameterName("targdiam", false);
  
  // SchamGasTgtCmd = new G4UIcmdWithAnInteger("/g4sbs/schbrflag",this);
  // SchamGasTgtCmd->SetGuidance("Build evacuated scattering chamber for gas target? (1=yes, 0=no)");
  // SchamGasTgtCmd->SetParameterName("schbrflag",false);

  beamcurCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamcur",this);
  beamcurCmd->SetGuidance("Beam current");
  beamcurCmd->SetParameterName("beamcur", false);

  runtimeCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/runtime",this);
  runtimeCmd->SetGuidance("Run time (used to convert event rate to counts, units must be defined in G4SystemOfUnits)");
  runtimeCmd->SetParameterName("runtime", false);

  rasterxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/rasterx",this);
  rasterxCmd->SetGuidance("Raster x size");
  rasterxCmd->SetParameterName("rasterx",true); // is omittable

  rasteryCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/rastery",this);
  rasteryCmd->SetGuidance("Raster y size");
  rasteryCmd->SetParameterName("rastery",true);
  
  rasterrCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/rasterR",this);
  rasterrCmd->SetGuidance("Raster radius size");
  rasterrCmd->SetParameterName("size",false);
  
  beamspotsizeCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamspotsize",this);
  beamspotsizeCmd->SetGuidance("beam spot size");
  beamspotsizeCmd->SetParameterName("size", false);
  
  // D. Flay 8/25/20.  Beam pointing and beam diffuser   
  // - horizontal (x)  
  beamOffsetXcmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamoffsetx",this);
  beamOffsetXcmd->SetGuidance("Set beam offset along the horizontal (x) direction");
  beamOffsetXcmd->SetParameterName("beamoffsetx",true);  // is omittable 
  // - vertical (y)
  beamOffsetYcmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamoffsety",this);
  beamOffsetYcmd->SetGuidance("Set beam offset along the vertical (y) direction");
  beamOffsetYcmd->SetParameterName("beamoffsety",true);   
  // beam dump 
  beamDumpCmd = new G4UIcmdWithABool("/g4sbs/beamDumpEnable",this);
  beamDumpCmd->SetGuidance("Enable the Beam Dump");
  beamDumpCmd->SetParameterName("beamDumpEnable",true);
  // beam diffuser 
  beamDiffuserCmd = new G4UIcmdWithABool("/g4sbs/beamDiffuserEnable",this);
  beamDiffuserCmd->SetGuidance("Enable the Beam Diffuser device");
  beamDiffuserCmd->SetParameterName("beamDiffuserEnable",true);

  tgtNfoilCmd = new G4UIcmdWithAnInteger("/g4sbs/Nfoil",this);
  tgtNfoilCmd->SetGuidance("Number of foils for optics target");
  tgtNfoilCmd->SetGuidance("Only has any effect if target = optics");
  tgtNfoilCmd->SetGuidance("Default = 1");
  tgtNfoilCmd->SetParameterName("nfoil",true);
  tgtNfoilCmd->SetDefaultValue(1);

  tgtFoilThickCmd = new G4UIcommand("/g4sbs/foilthick",this);
  tgtFoilThickCmd->SetGuidance("Foil thickness for multi-foil target");
  tgtFoilThickCmd->SetGuidance("Usage: /g4sbs/foilthick thick1 thick2 ... thickN unit");
  tgtFoilThickCmd->SetGuidance("Separate by whitespace, units at the end");
  tgtFoilThickCmd->SetGuidance("number of entries must match number of foils (see /g4sbs/Nfoil)");
  tgtFoilThickCmd->SetParameter( new G4UIparameter("foilthicklist",'s',false) );

  tgtFoilZCmd = new G4UIcommand("/g4sbs/foilz",this);
  tgtFoilZCmd->SetGuidance("Foil z positions along beamline for multi-foil target");
  tgtFoilZCmd->SetGuidance("Usage: /g4sbs/foilz z1 z2 ... zN unit");
  tgtFoilZCmd->SetGuidance("separate entries by whitespace, units at the end");
  tgtFoilZCmd->SetGuidance("number of entries must match number of foils (see /g4sbs/Nfoil)");
  tgtFoilZCmd->SetParameter( new G4UIparameter("foilthicklist",'s',false) );
  
  beamECmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamE",this);
  beamECmd->SetGuidance("Beam Energy");
  beamECmd->SetParameterName("energy", false);

  bbangCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/bbang",this);
  bbangCmd->SetGuidance("BigBite angle");
  bbangCmd->SetParameterName("angle", false);

  bbdistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/bbdist",this);
  bbdistCmd->SetGuidance("BigBite distance, target to front face of magnet yoke");
  bbdistCmd->SetParameterName("dist", false);

  hcalangCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/sbsang",this);
  hcalangCmd->SetGuidance("SBS angle");
  hcalangCmd->SetParameterName("angle", false);

  sbstrkrpitchCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/sbstrkrpitch",this);
  sbstrkrpitchCmd->SetGuidance("SBS tracker pitch angle (tilt toward up-bending particles)");
  sbstrkrpitchCmd->SetParameterName("angle", false);

  sbstrkrdistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/sbstrkrdist",this);
  sbstrkrdistCmd->SetGuidance( "SBS tracker distance from target to projection of center of first tracker plane onto horizontal plane");
  sbstrkrdistCmd->SetParameterName("sbstrkrdist",false );
  
  dvcsecalmatCmd = new G4UIcmdWithAString("/g4sbs/dvcsecalmat",this);
  dvcsecalmatCmd->SetGuidance("DVCS ECal material: 'PbF2' or 'PbWO4'");
  dvcsecalmatCmd->SetParameterName("dvcsecalmatname", false);

  GRINCH_gas_Cmd = new G4UIcmdWithAString("/g4sbs/grinchgas",this);
  GRINCH_gas_Cmd->SetGuidance("Gas for GRINCH detector: choose from C4F10, C4F8O, CF4, CO2, SF6 (C4F8, C3F8 coming soon)");
  GRINCH_gas_Cmd->SetParameterName("grinchgasname", false );

  RICH_gas_Cmd = new G4UIcmdWithAString("/g4sbs/richgas",this);
  RICH_gas_Cmd->SetGuidance("Gas for RICH detector: choose from C4F10, C4F8O, CF4, CO2, SF6 (C4F8, C3F8 coming soon)");
  RICH_gas_Cmd->SetParameterName("richgasname", false );
  
  hcaldistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hcaldist",this);
  hcaldistCmd->SetGuidance("HCAL distance");
  hcaldistCmd->SetParameterName("dist", false);

  hcalvoffsetCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hcalvoffset",this);
  hcalvoffsetCmd->SetGuidance("HCAL vertical offset");
  hcalvoffsetCmd->SetParameterName("dist", false);

  hcalhoffsetCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hcalhoffset",this);
  hcalhoffsetCmd->SetGuidance("HCAL horizontal offset relative to SBS center line (+ = TOWARD beam line)");
  hcalhoffsetCmd->SetParameterName("dist", false);

  hcalangoffsetCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hcalangoffset",this);
  hcalangoffsetCmd->SetGuidance("HCAL angular offset relative to exit beamline (+ = away from beamline)");
  hcalangoffsetCmd->SetParameterName("angle", false);

  CDetReadyCmd = new G4UIcmdWithABool("/g4sbs/cdetready",this);
  CDetReadyCmd->SetGuidance("Will CDet be ready or not for the experiment");
  CDetReadyCmd->SetParameterName("dist", false);

  lacdistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/lacdist",this);
  lacdistCmd->SetGuidance("LAC distance");
  lacdistCmd->SetParameterName("dist", false);

  lacvoffsetCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/lacvoffset",this);
  lacvoffsetCmd->SetGuidance("LAC vertical offset");
  lacvoffsetCmd->SetParameterName("dist", false);

  lachoffsetCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/lachoffset",this);
  lachoffsetCmd->SetGuidance("LAC horizontal offset relative to SBS center line (+ = TOWARD beam line)");
  lachoffsetCmd->SetParameterName("dist", false);
  
  
  hmagdistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/48D48dist",this);
  hmagdistCmd->SetGuidance("48D48 distance");
  hmagdistCmd->SetParameterName("dist", false);

  gemresCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/gemres",this);
  gemresCmd->SetGuidance("GEM coordinate resolution");
  gemresCmd->SetParameterName("dist", false);

  // Detector position commands

  cerDisCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/cerdist",this);
  cerDisCmd->SetGuidance("Cerenkov distance from front GEM");
  cerDisCmd->SetParameterName("dist", false);

  cerDepCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/cerdepth",this);
  cerDepCmd->SetGuidance("Cerenkov gas depth");
  cerDepCmd->SetParameterName("dist", false);

  gemSepCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/gemsep",this);
  gemSepCmd->SetGuidance("GEM separation from front to back set (BigBite GEMs)");
  gemSepCmd->SetParameterName("dist", false);

  bbCalDistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/bbcaldist",this);
  bbCalDistCmd->SetGuidance("BigBite calorimeter distance from front GEM");
  bbCalDistCmd->SetParameterName("dist", false);

  thminCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/thmin",this);
  thminCmd->SetGuidance("Minimum electron generation polar angle");
  thminCmd->SetParameterName("angle", false);

  thmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/thmax",this);
  thmaxCmd->SetGuidance("Maximum electron generation polar angle");
  thmaxCmd->SetParameterName("angle", false);

  phminCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/phmin",this);
  phminCmd->SetGuidance("Minimum electron generation azimuthal angle");
  phminCmd->SetParameterName("angle", false);

  phmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/phmax",this);
  phmaxCmd->SetGuidance("Maximum electron generation azimuthal angle");
  phmaxCmd->SetParameterName("angle", false);

  HthminCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hthmin",this);
  HthminCmd->SetGuidance("Minimum hadron generation polar angle (SIDIS generator)");
  HthminCmd->SetParameterName("htheta",false);

  HthmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hthmax",this);
  HthmaxCmd->SetGuidance("Maximum hadron generation polar angle (SIDIS generator)");
  HthmaxCmd->SetParameterName("htheta",false);

  HphminCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hphmin",this);
  HphminCmd->SetGuidance("Minimum hadron generation azimuthal angle (SIDIS generator)");
  HphminCmd->SetParameterName("htheta",false);

  HphmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hphmax",this);
  HphmaxCmd->SetGuidance("Maximum hadron generation azimuthal angle (SIDIS generator)");
  HphmaxCmd->SetParameterName("htheta",false);
    
  EhminCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/ehmin",this);
  EhminCmd->SetGuidance("Minimum hadron generation energy (SIDIS generator)");
  EhminCmd->SetParameterName("ehmin",false);
    
  EhmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/ehmax",this);
  EhmaxCmd->SetGuidance("Maximum hadron generation energy (SIDIS generator)");
  EhmaxCmd->SetParameterName("ehmax",false);

  EeminCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/eemin",this);
  EeminCmd->SetGuidance("Minimum electron generation energy (SIDIS generator)");
  EeminCmd->SetParameterName("eemin",false);
    
  EemaxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/eemax",this);
  EemaxCmd->SetGuidance("Maximum electron generation energy (SIDIS generator)");
  EemaxCmd->SetParameterName("eemax",false);

  PionPhoto_tminCmd = new G4UIcmdWithADouble("/g4sbs/tmin",this);
  PionPhoto_tminCmd->SetGuidance("Minimum -t value for pion photoproduction event generation");
  PionPhoto_tminCmd->SetGuidance("Assumed to be positive and given in units of GeV^2");
  PionPhoto_tminCmd->SetParameterName("tmin",false);

  PionPhoto_tmaxCmd = new G4UIcmdWithADouble("/g4sbs/tmax",this);
  PionPhoto_tmaxCmd->SetGuidance("Minimum -t value for pion photoproduction event generation");
  PionPhoto_tmaxCmd->SetGuidance("Assumed to be given in units of GeV^2");
  PionPhoto_tmaxCmd->SetParameterName("tmax",false);

  PionPhoto_useradCmd = new G4UIcmdWithABool("/g4sbs/userad", this );
  PionPhoto_useradCmd->SetGuidance("Use radiator for pion photoproduction event generation (Default = true)");
  PionPhoto_useradCmd->SetParameterName("userad",true);
  PionPhoto_useradCmd->SetDefaultValue(true);

  PionPhoto_radthickCmd = new G4UIcmdWithADouble("/g4sbs/radthick", this );
  PionPhoto_radthickCmd->SetGuidance("Radiator thickness, in units of X_0 (default = 0.06)");
  PionPhoto_radthickCmd->SetParameterName("radthickX0",true);
  PionPhoto_radthickCmd->SetDefaultValue(0.06);

  PionPhoto_radzCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/radz",this);
  PionPhoto_radzCmd->SetGuidance("How far upstream of the target is the radiator?");
  PionPhoto_radzCmd->SetGuidance("+Z = distance upstream (no minus sign needed)");
  PionPhoto_radzCmd->SetParameterName("radzoff",true);
  PionPhoto_radzCmd->SetDefaultValue( 10.0*cm );
  
  RICHdistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/richdist",this);
  RICHdistCmd->SetGuidance("SBS RICH distance from target");
  RICHdistCmd->SetParameterName("dist",false);

  RICHhoffsetCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/richhoffset",this);
  RICHhoffsetCmd->SetGuidance("SBS RICH horizontal offset (wrt SBS axis, + = toward beamline)");
  RICHhoffsetCmd->SetParameterName("hoffset",false);

  RICHvoffsetCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/richvoffset",this);
  RICHvoffsetCmd->SetGuidance("SBS RICH vertical offset (wrt SBS axis, + = up)");
  RICHvoffsetCmd->SetParameterName("voffset",false);
  
  RICHaeroCmd = new G4UIcmdWithABool("/g4sbs/userichaero",this);
  RICHaeroCmd->SetGuidance("Toggle use of RICH aerogel (default = true)" );
  RICHaeroCmd->SetParameterName("useaero",true);
  RICHaeroCmd->SetDefaultValue(true);

  RICHSnoutExtensionCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/richsnoutext",this);
  RICHSnoutExtensionCmd->SetGuidance("extend the RICH snout for electron mode (default 0");
  RICHSnoutExtensionCmd->SetParameterName("snoutextension",false);

  SBSMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/sbsmagfield",this);
  SBSMagFieldCmd->SetGuidance("SBS uniform magnetic field value");
  SBSMagFieldCmd->SetParameterName("sbsbfield",false);

  EARM_ScaleFieldCmd = new G4UIcmdWithADouble("/g4sbs/scalebbfield",this);
  EARM_ScaleFieldCmd->SetGuidance("Scale factor applied to BigBite magnetic field");
  EARM_ScaleFieldCmd->SetParameterName("bbfieldscale",false);

  HARM_ScaleFieldCmd = new G4UIcmdWithADouble("/g4sbs/scalesbsfield",this);
  HARM_ScaleFieldCmd->SetGuidance("Scale factor applied to SBS magnetic field");
  HARM_ScaleFieldCmd->SetParameterName("sbsfieldscale",false);
  
  SBSFieldClampOptionCmd = new G4UIcmdWithAnInteger("/g4sbs/sbsclampopt",this);
  SBSFieldClampOptionCmd->SetGuidance("SBS field clamp configuration: 0=no clamp, 3=Front clamp only (GMN, GEN)), 2=Front and rear clamps (GEP, SIDIS)");
  SBSFieldClampOptionCmd->SetParameterName("sbsclampoption",false);

  SBSBeamlineConfCmd = new G4UIcmdWithAnInteger("/g4sbs/beamlineconfig",this);
  SBSBeamlineConfCmd->SetGuidance("SBS beamline configuration: 1: GEp; 2: GEn; 3 (def): GMn 2-4 pass; 4; GMn 5 pass");
  SBSBeamlineConfCmd->SetParameterName("beamlineconf",false);
  
  SBSLeadOptionCmd = new G4UIcmdWithAnInteger("/g4sbs/uselead",this);
  SBSLeadOptionCmd->SetGuidance("SBS beamline lead shielding configuration: 0= nope 1=yes");
  SBSLeadOptionCmd->SetParameterName("uselead",false);

  GENRPAnalyzerOptionCmd = new G4UIcmdWithAnInteger("/g4sbs/genrpAnalyzer",this);
  GENRPAnalyzerOptionCmd->SetGuidance("GEnRP Analyzer configuration: 0=none+no beamline PR; 1=none, 2=Cu+Gla(para), 3=Cu+Gla(perp), 4=Cu+CGEN");
  GENRPAnalyzerOptionCmd->SetParameterName("genrpAnalyzer",false);

  GEPFPP1_CH2thickCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/FPP1CH2thick",this);
  GEPFPP1_CH2thickCmd->SetGuidance("CH2 thickness for first analyzer (GEP only)");
  GEPFPP1_CH2thickCmd->SetGuidance("0 < FPP1 CH2 thick < 150 cm");
  GEPFPP1_CH2thickCmd->SetParameterName("CH2thick1",false);
  //GEPFPP1_CH2thickCmd->SetRange("0.0 <= CH2thick1 && CH2thick1 <= 60.0*cm"

  GEPFPP2_CH2thickCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/FPP2CH2thick",this);
  GEPFPP2_CH2thickCmd->SetGuidance("CH2 thickness for first analyzer (GEP only)");
  GEPFPP2_CH2thickCmd->SetGuidance("0 < FPP2 CH2 thick < 150 cm");
  GEPFPP2_CH2thickCmd->SetParameterName("CH2thick2",false);

  GEPFPPoptionCmd = new G4UIcmdWithAnInteger("/g4sbs/gepfppoption",this);
  GEPFPPoptionCmd->SetGuidance("GEP FPP option:");
  GEPFPPoptionCmd->SetGuidance("1 = One analyzer, 8 (FT) + 8 (FPP) GEM trackers");
  GEPFPPoptionCmd->SetGuidance("2 = Two analyzers, 6 (FT) + 5 (FPP1) + 5 (FPP2) GEM trackers (default)");
  GEPFPPoptionCmd->SetGuidance("3 = Same as 2, but with second analyzer replaced by 3.5\" steel from GEN-RP");
  GEPFPPoptionCmd->SetParameterName("gepfppoption",true);
  GEPFPPoptionCmd->SetDefaultValue(2);

  HadronFilterCmd = new G4UIcmdWithABool("/g4sbs/usehadronfilter",this);
  HadronFilterCmd->SetGuidance("Use target shield wall for GEP/GEN/GEN-RP/SIDIS (yes/no)");

  HadronFilterThickCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hadronfilterthick",this);
  HadronFilterThickCmd->SetGuidance("target shield wall thickness");
  HadronFilterThickCmd->SetParameterName("hadronfilterthick",false);

  HadronFilterMaterialCmd = new G4UIcmdWithAString("/g4sbs/hadronfiltermaterial",this);
  HadronFilterMaterialCmd->SetGuidance("Material for target shield wall");
  HadronFilterMaterialCmd->SetGuidance("any valid material defined in G4SBSDetectorConstruction::ConstructMaterials()");
  HadronFilterMaterialCmd->SetParameterName("hadronfiltermaterial",false);
  
  BLneutronDetsCmd = new G4UIcmdWithABool("/g4sbs/BLneutronDets",this);
  BLneutronDetsCmd->SetGuidance("Setup neutron detectors along the beamline");
  BLneutronDetsCmd->SetParameterName("switch", false);
  
  GEMfrontendCmd = new G4UIcmdWithABool("/g4sbs/buildGEMfrontend",this);
  GEMfrontendCmd->SetGuidance("build GEM front end for GMn or GEp");
  GEMfrontendCmd->SetParameterName("switch", false);
  
  GEMfrontendDistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/GEMfrontendDist",this);
  GEMfrontendDistCmd->SetGuidance("GEM front end distance");
  GEMfrontendDistCmd->SetParameterName("GEMfrontendDist", false);
  
  GEMfrontendPosAngleCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/GEMfrontendPosang",this);
  GEMfrontendPosAngleCmd->SetGuidance("GEM front end position angle");
  GEMfrontendPosAngleCmd->SetParameterName("GEMfrontendPosAngle", false);
  
  GEMfrontendRotAngleCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/GEMfrontendRotang",this);
  GEMfrontendRotAngleCmd->SetGuidance("GEM front end position angle");
  GEMfrontendRotAngleCmd->SetParameterName("GEMfrontendRotAngle", false);
  
  SetGrinchPMTglassHitsCmd = new G4UIcmdWithABool("/g4sbs/GrinchPMTglassHits",this);
  SetGrinchPMTglassHitsCmd->SetGuidance("build GEM front end for GMn or GEp");
  SetGrinchPMTglassHitsCmd->SetParameterName("switch", false);  
  
  buildSBSsieveCmd = new G4UIcmdWithABool("/g4sbs/buildSBSsieve",this);
  buildSBSsieveCmd->SetGuidance("Use SBS sieve (true or false, false by default)");
  buildSBSsieveCmd->SetParameterName("buildSBSsieve",false);

  //SSeeds - converted to integer input for multiple sieve plate option 10.4.20
  buildBBsieveCmd = new G4UIcmdWithAnInteger("/g4sbs/buildBBsieve",this);
  buildBBsieveCmd->SetGuidance("BB Ecal shielding layout: option 0 (none), 1 (straight holes and slots), 2 (Holes at dispersive angle in x and y)");
  buildBBsieveCmd->SetParameterName("buildBBsieve", false);
  
  TreeFlagCmd = new G4UIcmdWithAnInteger("/g4sbs/treeflag",this);
  TreeFlagCmd->SetGuidance("G4SBS ROOT tree filling: 0=keep all, 1=keep only evts w/hits in sensitive volumes");
  TreeFlagCmd->SetParameterName("treeflag",false);

  SBS_FT_absorberCmd = new G4UIcmdWithABool("/g4sbs/FTabsorberflag",this);
  SBS_FT_absorberCmd->SetGuidance("Turn on sheet of absorber material in front of SBS FT (only applicable to GEP)");
  SBS_FT_absorberCmd->SetParameterName("FTabsflag",false);

  SBS_FT_absorberMaterialCmd = new G4UIcmdWithAString("/g4sbs/FTabsorbermaterial",this);
  SBS_FT_absorberMaterialCmd->SetGuidance("Set material of FT absorber (default = aluminum)" );
  SBS_FT_absorberMaterialCmd->SetGuidance("Valid material names defined in G4SBSDetectorConstruction::ConstructMaterials()");
  SBS_FT_absorberMaterialCmd->SetParameterName("FTabsmat",true);
  SBS_FT_absorberMaterialCmd->SetDefaultValue("Aluminum");

  SBS_FT_absorberThickCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/FTabsorberthick",this);
  SBS_FT_absorberThickCmd->SetGuidance("Set thickness of absorber in front of FT (default = 1 inch)");
  SBS_FT_absorberThickCmd->SetParameterName("FTabsthick",true);
  SBS_FT_absorberThickCmd->SetDefaultValue(2.54*cm);
  
  // Earm_CAL_part_cmd = new G4UIcmdWithABool("/g4sbs/keep_part_earm_cal",this);
  // Earm_CAL_part_cmd->SetGuidance("Keep particle info in root tree for electron arm calorimeter (default = false)");
  // Earm_CAL_part_cmd->SetParameterName("keeppart",true);
  // Earm_CAL_part_cmd->SetDefaultValue( true );

  // Harm_CAL_part_cmd = new G4UIcmdWithABool("/g4sbs/keep_part_harm_cal",this);
  // Harm_CAL_part_cmd->SetGuidance("Keep particle info in root tree for hadron arm calorimeter (default = false)" );
  // Harm_CAL_part_cmd->SetParameterName("keeppart",true);
  // Harm_CAL_part_cmd->SetDefaultValue( true );

  KeepPartCALcmd = new G4UIcommand("/g4sbs/keepcalparticles",this);
  KeepPartCALcmd->SetGuidance("Keep particle information in ROOT tree for sensitive volumes defined as calorimeters");
  KeepPartCALcmd->SetGuidance("Usage: /g4sbs/keepcalparticles name flag");
  KeepPartCALcmd->SetGuidance("name = sensitive detector name");
  KeepPartCALcmd->SetGuidance("flag = true/false (default = false)");
  KeepPartCALcmd->SetParameter( new G4UIparameter("sdname",'s',false) );
  KeepPartCALcmd->SetParameter( new G4UIparameter("flag",'b',true) );
  KeepPartCALcmd->GetParameter(1)->SetDefaultValue(false);

  KeepHistorycmd = new G4UIcommand("/g4sbs/keephistory",this);
  KeepHistorycmd->SetGuidance("Keep entire particle history in the root tree for any sensitive detector");
  KeepHistorycmd->SetGuidance("Usage: /g4sbs/keephistory name flag");
  KeepHistorycmd->SetGuidance("name = sensitive detector name");
  KeepHistorycmd->SetGuidance("flag = true/false (default = false)" );
  KeepHistorycmd->SetGuidance("note: /tracking/storeTrajectory 1 must be invoked for availability of particle history at end of event");
  KeepHistorycmd->SetParameter( new G4UIparameter("sdname", 's', false ) );
  KeepHistorycmd->SetParameter( new G4UIparameter("flag",'b',true) );
  KeepHistorycmd->GetParameter(1)->SetDefaultValue(false);

  LimitStepCALcmd = new G4UIcommand("/g4sbs/steplimit",this);
  LimitStepCALcmd->SetGuidance("Set max step length to zero for particles entering sensitive volume (stops and kills particles)");
  LimitStepCALcmd->SetGuidance("Usage: /g4sbs/steplimit name flag" );
  LimitStepCALcmd->SetGuidance("name = sensitive detector name");
  LimitStepCALcmd->SetGuidance("flag = true/false (default = true)");
  //LimitStepCALcmd->SetGuidance("flag = true/false (default = false)" );
  LimitStepCALcmd->SetParameter( new G4UIparameter("sdname", 's', false ) );
  LimitStepCALcmd->SetParameter( new G4UIparameter("flag",'b',true) );
  LimitStepCALcmd->GetParameter(1)->SetDefaultValue(true);

  // **********
  SD_EnergyThresholdCmd = new G4UIcommand("/g4sbs/threshold",this );

  SD_EnergyThresholdCmd->SetGuidance("Set hit energy threshold by sensitive detector name (only valid for kCAL, kGEM)");
  SD_EnergyThresholdCmd->SetGuidance("Usage: /g4sbs/threshold SDname ethresh unit");
  SD_EnergyThresholdCmd->SetParameter( new G4UIparameter("sdname", 's', false) );
  SD_EnergyThresholdCmd->SetParameter( new G4UIparameter("threshold", 'd', false) );
  SD_EnergyThresholdCmd->SetParameter( new G4UIparameter("unit", 's', false) );

  SD_TimeWindowCmd = new G4UIcommand("/g4sbs/timewindow",this);
  SD_TimeWindowCmd->SetGuidance( "Set hit timing window by sensitive detector name (only valid for kCAL, kGEM, kECAL)" );
  SD_TimeWindowCmd->SetGuidance( "Usage: /g4sbs/timewindow SDname twindow unit" );
  SD_TimeWindowCmd->SetParameter( new G4UIparameter("sdname", 's', false ) );
  SD_TimeWindowCmd->SetParameter( new G4UIparameter("timewindow", 'd', false ) );
  SD_TimeWindowCmd->SetParameter( new G4UIparameter("unit", 's', false) );

  SD_NTimeBinsCmd = new G4UIcommand("/g4sbs/ntimebins",this);
  SD_NTimeBinsCmd->SetGuidance( "Set number of time bins by sensitive detector name (only valid for kCAL, kGEM, kECAL)" );
  SD_NTimeBinsCmd->SetGuidance( "Usage: /g4sbs/ntimebins" );
  SD_NTimeBinsCmd->SetParameter( new G4UIparameter("sdname", 's', false ) );
  SD_NTimeBinsCmd->SetParameter( new G4UIparameter("ntimebins", 'i', false ) );

  KeepPulseShapeCmd = new G4UIcommand("/g4sbs/keeppulseshapeinfo",this);
  KeepPulseShapeCmd->SetGuidance("Toggle recording of Pulse Shape info in the tree by SD name");
  KeepPulseShapeCmd->SetGuidance("Usage: /g4sbs/keeppulseshapeinfo SDname flag");
  KeepPulseShapeCmd->SetGuidance("SDname = sensitive detector name");
  KeepPulseShapeCmd->SetGuidance("flag = true/false or 0/1 (default = false/0)");
  KeepPulseShapeCmd->SetParameter( new G4UIparameter("sdname", 's', false ) );
  KeepPulseShapeCmd->SetParameter( new G4UIparameter("flag", 'b', false) );
  KeepPulseShapeCmd->GetParameter(1)->SetDefaultValue(false);    
  // **********

  KeepSDtrackcmd = new G4UIcommand("/g4sbs/keepsdtrackinfo",this);
  KeepSDtrackcmd->SetGuidance("Toggle recording of SD track info in the tree by SD name");
  KeepSDtrackcmd->SetGuidance("Usage: /g4sbs/keepsdtrackinfo SDname flag");
  KeepSDtrackcmd->SetGuidance("SDname = sensitive detector name");
  KeepSDtrackcmd->SetGuidance("flag = true/false or 0/1 (default = true/1)");
  KeepSDtrackcmd->SetParameter( new G4UIparameter("sdname", 's', false ) );
  KeepSDtrackcmd->SetParameter( new G4UIparameter("flag", 'b', true) );
  KeepSDtrackcmd->GetParameter(1)->SetDefaultValue(true);
  
  // DisableOpticalPhysicsCmd = new G4UIcmdWithABool("/g4sbs/useopticalphysics", this );
  // DisableOpticalPhysicsCmd->SetGuidance("toggle optical physics on/off");
  // DisableOpticalPhysicsCmd->SetGuidance("default = true (ON)");
  // DisableOpticalPhysicsCmd->SetParameterName("optphys",true);
  
  UseCerenkovCmd = new G4UIcmdWithABool( "/g4sbs/useckov",this );
  UseCerenkovCmd->SetGuidance( "Toggle Cerenkov process on/off (default = ON)" );
  UseCerenkovCmd->SetParameterName("useckov",true);
  UseCerenkovCmd->AvailableForStates(G4State_PreInit);
  
  UseScintCmd = new G4UIcmdWithABool( "/g4sbs/usescint",this );
  UseScintCmd->SetGuidance( "Toggle Scintillation process on/off (default = ON)" );
  UseScintCmd->SetParameterName("usescint",true);
  UseScintCmd->AvailableForStates(G4State_PreInit);

  DisableOpticalPhotonProductionByMaterialCmd = new G4UIcmdWithAString( "/g4sbs/DisableOpticalPhotonByMaterial",this );
  DisableOpticalPhotonProductionByMaterialCmd->SetGuidance( "Disable optical photon production by material name" );
  DisableOpticalPhotonProductionByMaterialCmd->SetGuidance( "Works by preventing definition of optical properties for a given material" );
  DisableOpticalPhotonProductionByMaterialCmd->SetGuidance( "Use with caution" );
  DisableOpticalPhotonProductionByMaterialCmd->SetParameterName("material",false);
  DisableOpticalPhotonProductionByMaterialCmd->AvailableForStates(G4State_PreInit);
  
  
  FluxCmd = new G4UIcmdWithABool("/g4sbs/fluxcalc",this);
  FluxCmd->SetGuidance( "Compute particle flux as a function of angles, energy");
  FluxCmd->SetParameterName( "fluxcalc", false);  

  TargPolDirectionCmd = new G4UIcmdWith3Vector("/g4sbs/targpoldirection", this );
  TargPolDirectionCmd->SetGuidance("Set target polarization direction");
  TargPolDirectionCmd->SetGuidance("Three-vector arguments are x,y,z components of polarization");
  TargPolDirectionCmd->SetGuidance("Automatically converted to unit vector internally");
  TargPolDirectionCmd->SetGuidance("Assumed to be given in global coordinate system");
  TargPolDirectionCmd->SetParameterName("Px","Py","Pz",false);
  
  TargPolMagnitudeCmd = new G4UIcmdWithADouble("/g4sbs/targpolmag", this );
  TargPolMagnitudeCmd->SetGuidance("Set target polarization magnitude");
  TargPolMagnitudeCmd->SetGuidance("0 <= Ptarg <= 1");
  TargPolMagnitudeCmd->SetParameterName("Ptgt",true);
  TargPolMagnitudeCmd->SetDefaultValue(1.0);

  BeamPolDirectionCmd = new G4UIcmdWith3Vector("/g4sbs/beampoldirection", this );
  BeamPolDirectionCmd->SetGuidance("Set beam polarization direction");
  BeamPolDirectionCmd->SetGuidance("Three-vector arguments are x,y,z components of polarization");
  BeamPolDirectionCmd->SetGuidance("Automatically converted to unit vector internally");
  BeamPolDirectionCmd->SetGuidance("Assumed to be given in global coordinate system");
  BeamPolDirectionCmd->SetParameterName("Px","Py","Pz",false);
  
  BeamPolMagnitudeCmd = new G4UIcmdWithADouble("/g4sbs/beampolmag", this );
  BeamPolMagnitudeCmd->SetGuidance("Set beam polarization magnitude");
  BeamPolMagnitudeCmd->SetGuidance("Not yet used by anything, but anticipated for use in polarized cross section calculation");
  BeamPolMagnitudeCmd->SetGuidance("0 <= Pbeam <= 1");
  BeamPolMagnitudeCmd->SetParameterName("Pbeam",true);
  BeamPolMagnitudeCmd->SetDefaultValue(1.0);

  RandomizeTargetSpinCmd = new G4UIcmdWithABool("/g4sbs/randomizetargetspin",this);
  RandomizeTargetSpinCmd->SetGuidance("Turn on randomization of target spin direction in event generator");
  RandomizeTargetSpinCmd->SetParameterName("randomizespin", true );
  RandomizeTargetSpinCmd->SetDefaultValue(true);

  NumSpinStatesTargCmd = new G4UIcmdWithAnInteger("/g4sbs/numtargspinstates",this);
  NumSpinStatesTargCmd->SetGuidance("Specify number of target spin directions to randomly generate");
  //NumSpinStatesTargCmd->SetGuidance("");
  NumSpinStatesTargCmd->SetParameterName("Ntargspin", false);

  TargThetaSpinCmd = new G4UIcommand("/g4sbs/targthetaspin",this);
  TargThetaSpinCmd->SetGuidance("Specify polar angle(s) of target spin for randomized target spin direction in event generator");
  TargThetaSpinCmd->SetGuidance("Number of angles must match number of target spin states (see /g4sbs/numtargspinstates)");
  TargThetaSpinCmd->SetGuidance("Entries separated by whitespace, units at the end of the string (degrees or radians)"); 
  TargThetaSpinCmd->SetParameter( new G4UIparameter("thetaspinlist",'s',false) ); //parameter is of string type

  TargPhiSpinCmd = new G4UIcommand("/g4sbs/targphispin",this);
  TargPhiSpinCmd->SetGuidance("Specify azimuthal angle(s) of target spin for randomized target spin direction in event generator");
  TargPhiSpinCmd->SetGuidance("Number of angles must match number of target spin states (see /g4sbs/numtargspinstates)");
  TargPhiSpinCmd->SetGuidance("Entries separated by whitespace, units at the end of the string (degrees or radians)"); 
  TargPhiSpinCmd->SetParameter( new G4UIparameter("phispinlist",'s',false) ); //parameter is of string type
  
  GunPolarizationCommand = new G4UIcmdWith3Vector( "/g4sbs/gunpol", this );
  GunPolarizationCommand->SetGuidance( "Set particle polarization for gun generator:" );
  GunPolarizationCommand->SetGuidance( "Three-vector arguments are x,y,z components of polarization" );
  GunPolarizationCommand->SetGuidance( "Automatically converted to unit vector internally" );
  GunPolarizationCommand->SetGuidance( "Assumed to be given in TRANSPORT coordinates" );
  GunPolarizationCommand->SetParameterName("Sx","Sy","Sz",false);
  
  SegmentC16Cmd = new G4UIcmdWithAnInteger( "/g4sbs/segmentTF1", this );
  SegmentC16Cmd->SetGuidance( "Longitudinally segment the TF1 lead glass for ECAL/C16 into N segments for thermal annealing model" );
  SegmentC16Cmd->SetGuidance( "0 = OFF (one segment, optical properties based on no rad damage, no temperature increase)" );
  SegmentC16Cmd->SetGuidance( "1..N: Use N segments; linear temperature profile assumed" );
  SegmentC16Cmd->SetParameterName("N",false);

  SegmentThickC16Cmd = new G4UIcmdWithADoubleAndUnit( "/g4sbs/segthickTF1", this );
  SegmentThickC16Cmd->SetGuidance( "Longitudinal segment thickness for TF1 lead glass of ECAL/C16 for thermal annealing model (default = 4 cm)" );
  SegmentThickC16Cmd->SetParameterName("thick",false);

  DoseRateCmd = new G4UIcmdWithADouble("/g4sbs/doserate", this );
  DoseRateCmd->SetGuidance( "Overall scale factor for dose rate in lead-glass for ECAL/C16 (depth profile is hard-coded!)");
  //DoseRateCmd->SetGuidance( "Assumed to be given in units of krad/hour" ); //Note 1 rad = 0.01 J/kg
  DoseRateCmd->SetParameterName("rate",false);
  
  CosmicsPointerCommand = new G4UIcmdWith3VectorAndUnit( "/g4sbs/cosmicpointer", this );
  CosmicsPointerCommand->SetGuidance( "Set pointer for cosmics:" );
  CosmicsPointerCommand->SetGuidance( "Three-vector arguments are x,y,z;" );
  CosmicsPointerCommand->SetGuidance( "please provide unit" );
  CosmicsPointerCommand->SetParameterName("x_ptr","y_ptr","z_ptr", false);

  CosmicsPointerRadiusCommand = new G4UIcmdWithADoubleAndUnit( "/g4sbs/cosmicpointerradius", this );
  CosmicsPointerRadiusCommand->SetGuidance( "Set pointer radius for cosmics (please provide unit)" );
  CosmicsPointerRadiusCommand->SetParameterName("radius",false);

  CosmicsMaxAngleCommand = new G4UIcmdWithADoubleAndUnit( "/g4sbs/cosmicmaxangle", this );
  CosmicsMaxAngleCommand->SetGuidance( "Set max zenithal angle for cosmics (please provide unit)" );
  CosmicsMaxAngleCommand->SetParameterName("maxangle",false);

  WriteFieldMapCmd = new G4UIcmdWithABool( "/g4sbs/writefieldmaps", this );
  WriteFieldMapCmd->SetGuidance( "Toggle writing of \"portable\" field maps for BB+SBS" );
  WriteFieldMapCmd->SetParameterName("writemapsflag", false );

  UseGEMshieldCmd = new G4UIcmdWithABool( "/g4sbs/usegemshielding", this );
  UseGEMshieldCmd->SetGuidance( "Include thin aluminum GEM shielding (for noise reduction)" );
  UseGEMshieldCmd->SetParameterName("GEMshieldflag", false );

  GEMshieldThickCmd = new G4UIcmdWithADoubleAndUnit( "/g4sbs/gemshieldthick", this );
  GEMshieldThickCmd->SetGuidance( "Thickness of GEM aluminum shielding (give value and unit, default = 50 um)");
  GEMshieldThickCmd->SetParameterName( "gemshieldthick", true );
  GEMshieldThickCmd->SetDefaultValue( 50.0*CLHEP::um );

  GEMshieldAirGapThickCmd = new G4UIcmdWithADoubleAndUnit( "/g4sbs/gemshieldairgapthick", this );
  GEMshieldAirGapThickCmd->SetGuidance( "Thickness of air gap between aluminum shielding (in front) and GEM (give value and unit, default = 0)" );
  GEMshieldAirGapThickCmd->SetParameterName( "gemshieldairgapthick", true );
  GEMshieldAirGapThickCmd->SetDefaultValue( 0.0*CLHEP::um );

  EnableBigBitePlateCmd = new G4UIcmdWithABool( "/g4sbs/setbigbiteplate", this );
  EnableBigBitePlateCmd->SetGuidance( "setup a plate in front of BigBite, default false" );
  EnableBigBitePlateCmd->SetParameterName( "setbigbiteplate", false );
  EnableBigBitePlateCmd->SetDefaultValue( false );
  
  SetBigBitePlateThicknessCmd = new G4UIcmdWithADoubleAndUnit( "/g4sbs/bigbiteplatethick", this );
  SetBigBitePlateThicknessCmd->SetGuidance( "set bigbite plate thickness, default 2.54cm" );
  SetBigBitePlateThicknessCmd->SetParameterName( "bigbiteplatethick", false );
  SetBigBitePlateThicknessCmd->SetDefaultValue( 2.54*CLHEP::cm );
  
  SetBigBitePlateMaterialCmd = new G4UIcmdWithAString( "/g4sbs/bigbiteplatematerial", this );
  SetBigBitePlateMaterialCmd->SetGuidance( "set bigbite plate material, default CH2" );
  SetBigBitePlateMaterialCmd->SetParameterName( "bigbiteplatematerial", false );
  SetBigBitePlateMaterialCmd->SetDefaultValue( "CH2" );
  
}

G4SBSMessenger::~G4SBSMessenger(){
}


void G4SBSMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue){
  char cmdstr[255];

  if(cmd==printCmd){
     G4int lineNo = printCmd->GetNewIntValue(newValue); 
     std::cout << "*************************** The line number is " << lineNo << std::endl;
  }

  if( cmd == runCmd ){
	
    G4VPhysicalVolume* pWorld;

    G4int nevt = runCmd->GetNewIntValue(newValue);

    //If the generator is PYTHIA, don't try to generate more events than we have available:
    if( fevgen->GetKine() == G4SBS::kPYTHIA6 ){
      if( fevgen->GetPythiaChain()->GetEntries() < nevt ){
	nevt = fevgen->GetPythiaChain()->GetEntries();
      }
      fevgen->InitializePythia6_Tree();
    }
    if( fevgen->GetKine() == G4SBS::kSIMC ){
      if( fevgen->GetSIMCChain()->GetEntries() < nevt ){
	nevt = fevgen->GetSIMCChain()->GetEntries();
      }
      fevgen->InitializeSIMC_Tree();
    }

    //    G4double TargMassDensity;
    G4double TargNumberDensity; 
    
    switch(fdetcon->fTargType){ //Initialize fTargDen correctly and consistently with Material definition in fdetcon->ConstructMaterials:
    case G4SBS::kH2:
      TargNumberDensity = fdetcon->GetMaterial("refH2")->GetTotNbOfAtomsPerVolume();
      break;
    case G4SBS::kD2:
      TargNumberDensity = fdetcon->GetMaterial("refD2")->GetTotNbOfAtomsPerVolume();
      break;
    case G4SBS::kNeutTarg:
      TargNumberDensity = fdetcon->GetMaterial("refN2")->GetTotNbOfAtomsPerVolume();
      break;
    case G4SBS::kLH2:
      TargNumberDensity = fdetcon->GetMaterial("LH2")->GetTotNbOfAtomsPerVolume();
      break;
    case G4SBS::kLD2:
      TargNumberDensity = fdetcon->GetMaterial("LD2")->GetTotNbOfAtomsPerVolume();
      break;
    case G4SBS::k3He:
      TargNumberDensity = fdetcon->GetMaterial("pol3He")->GetTotNbOfAtomsPerVolume();
      break;
    case G4SBS::kCfoil:
      TargNumberDensity = fdetcon->GetMaterial("Carbon")->GetTotNbOfAtomsPerVolume();
      break;
    default: //assume H2 gas:
      TargNumberDensity = fdetcon->GetMaterial("refH2")->GetTotNbOfAtomsPerVolume();
      break;
    }
    fevgen->SetTargDen(TargNumberDensity);
    
    fevgen->SetNevents(nevt);
    fevgen->Initialize();

    //For optics target, copy target foil information from targetbuilder to evgen:
    if( fdetcon->fTargType == G4SBS::kOptics ){
      fevgen->SetNfoils( fdetcon->fTargetBuilder->GetNtargetFoils() );
      fevgen->SetFoilZandThick( fdetcon->fTargetBuilder->GetFoilZpos(), fdetcon->fTargetBuilder->GetFoilThick() );
    }
    // if( fevgen->GetRejectionSamplingFlag() && !fevgen->GetRejectionSamplingInitialized() ){
    //   fevgen->InitializeRejectionSampling();
    // }

    G4SBS::Kine_t kinetype = fevgen->GetKine(); 
    if( kinetype == G4SBS::kDIS || kinetype == G4SBS::kWiser ){ //Processes with xsec in units of area/energy/solid angle; i.e., nb/GeV/sr
      G4SBSRun::GetRun()->GetData()->SetGenVol( fevgen->GetGenVol()/GeV );
      //if( fevgen->GetRejectionSamplingFlag() ){
      G4SBSRun::GetRun()->GetData()->SetMaxWeight( fevgen->GetMaxWeight()/cm2 * GeV );
	//}
    } else if ( kinetype == G4SBS::kSIDIS ){ //Processes with xsec differential in area/energy^2/solid angle^2:
      G4SBSRun::GetRun()->GetData()->SetGenVol( fevgen->GetGenVol()/(GeV*GeV) );
      // if( fevgen->GetRejectionSamplingFlag() ){
      G4SBSRun::GetRun()->GetData()->SetMaxWeight( fevgen->GetMaxWeight()/cm2 * pow(GeV,2) );
      //}
    } else { //Processes with xsec differential in solid angle only:
      G4SBSRun::GetRun()->GetData()->SetGenVol( fevgen->GetGenVol() );
      //if( fevgen->GetRejectionSamplingFlag() ){
      G4SBSRun::GetRun()->GetData()->SetMaxWeight( fevgen->GetMaxWeight()/cm2 );
      //}
    }
    G4SBSRun::GetRun()->GetData()->SetLuminosity( fevgen->GetLumi()*cm2*s );
    //G4SBSRun::GetRun()->GetData()->SetMaxWeight( fevgen->GetMaxWeight() );
    
    //Clean out and rebuild the detector geometry from scratch: 

    G4SolidStore::GetInstance()->Clean();
    G4LogicalVolumeStore::GetInstance()->Clean();
    G4PhysicalVolumeStore::GetInstance()->Clean();
	
    G4RunManager::GetRunManager()->DefineWorldVolume(pWorld = fdetcon->ConstructAll());
    G4RunManager::GetRunManager()->GeometryHasBeenModified();

    //Copy sensitive detector list and type codes to event action:
    fevact->SDlist = fdetcon->SDlist;
    fevact->SDtype = fdetcon->SDtype;
    //fevact->SDarm = fdetcon->SDarm;

    fIO->SetDetCon( fdetcon );
    //The following was moved to BeginOfRunAction.
    //fIO->InitializeTree();
    //ftrkact->SetDetCon( fdetcon );
    //Copy list of "analyzer" and "target" volumes and mapping between SD names and "boundary volumes" to
    //the following was moved to BeginOfRunAction
    // ftrkact->Initialize( fdetcon );
    // fstepact->Initialize( fdetcon );
    
    // G4cout << "InitializeTree() successful" << G4endl;

    // Clobber old gdml if it exists and write out the
    // present geometry
    // Save geometry to GDML file
#ifdef G4SBS_USE_GDML
    G4GDMLParser parser;
    unlink("g4sbs.gdml");
    parser.Write("g4sbs.gdml", pWorld);
#endif
    // Run the simulation
    G4UImanager * UImanager = G4UImanager::GetUIpointer();
    sprintf(cmdstr, "/run/beamOn %d", nevt);
    UImanager->ApplyCommand(cmdstr);
  }

  if( cmd == fileCmd ){
    fIO->SetFilename(newValue.data());
    G4SBSRun::GetRun()->GetData()->SetFileName(newValue);
  }

  if( cmd == sigfileCmd ){
    fevact->LoadSigmas(newValue.data());
  }

  if( cmd == gemconfigCmd ){
    G4int gemconfval = gemconfigCmd->GetNewIntValue(newValue);
    fdetcon->fEArmBuilder->SetGEMConfig(gemconfval);
  }

  if( cmd == shieldconfigCmd ){
    G4int shieldconfval = shieldconfigCmd->GetNewIntValue(newValue);
    fdetcon->fEArmBuilder->SetShieldConfig(shieldconfval);
  }
  
  if( cmd == bbpsconfigCmd ){
    G4int bbpsconfval = shieldconfigCmd->GetNewIntValue(newValue);
    fdetcon->fEArmBuilder->SetBBPSOption(bbpsconfval);
  }
  
  if( cmd == CDetconfigCmd ){
    G4int cdetconf = CDetconfigCmd->GetNewIntValue(newValue);
    fdetcon->SetCDetconfig(cdetconf);
  }

  if( cmd == flipGEMCmd ){
    G4bool val = flipGEMCmd->GetNewBoolValue(newValue);
    fdetcon->SetFlipGEM(val);
  }
  
  if( cmd == SegmentC16Cmd ){
    G4int segmentC16 = SegmentC16Cmd->GetNewIntValue(newValue);

    //G4cout << "segmentC16 = " << segmentC16 << G4endl;
    fdetcon->SetC16Segmentation( segmentC16 );
  }

  if( cmd == SegmentThickC16Cmd ){
    G4double thick = SegmentThickC16Cmd->GetNewDoubleValue(newValue);
    fdetcon->SetSegmentThickC16( thick );
  }

  if( cmd == DoseRateCmd ){
    G4double rate = DoseRateCmd->GetNewDoubleValue(newValue);
    fdetcon->SetDoseRateC16( rate );
  } 
  
  if( cmd == ECALmapfileCmd ){
    fdetcon->SetECALmapfilename( newValue );
  }

  if( cmd == kineCmd ){
    bool validcmd = false;

    G4SBS::Kine_t kinetemp = G4SBS::kElastic;
    
    if( newValue.compareTo("elastic") == 0 ){
      kinetemp = G4SBS::kElastic;
      // fevgen->SetKine(G4SBS::kElastic);
      // fIO->SetKine(G4SBS::kElastic);
      validcmd = true;
    }
    if( newValue.compareTo("inelastic") == 0 ){
      kinetemp = G4SBS::kInelastic;
      // fevgen->SetKine(G4SBS::kInelastic);
      //fevgen->SetRejectionSamplingFlag(false);
      validcmd = true;
    }
    if( newValue.compareTo("flat") == 0 ){
      kinetemp = G4SBS::kFlat;
      //fevgen->SetKine(G4SBS::kFlat);
      //fevgen->SetMaxWeight( cm2 );
      fevgen->SetRejectionSamplingFlag(false);
      validcmd = true;
    }
    if( newValue.compareTo("dis") == 0 ){
      fevgen->SetKine(G4SBS::kDIS);
      //fevgen->SetMaxWeight( cm2/GeV );
      validcmd = true;
    }
    if( newValue.compareTo("beam") == 0 ){
      kinetemp = G4SBS::kBeam;
      //fevgen->SetKine(G4SBS::kBeam);
      fevgen->SetRejectionSamplingFlag(false);
      validcmd = true;
    }
    if( newValue.compareTo("sidis") == 0 ){
      kinetemp = G4SBS::kSIDIS;
      //fevgen->SetKine( G4SBS::kSIDIS );
      //fevgen->SetMaxWeight( cm2/pow(GeV,2) );
      validcmd = true;
    }
    if( newValue.compareTo("wiser") == 0 ){
      kinetemp = G4SBS::kWiser;
      //      fevgen->SetKine( G4SBS::kWiser);
      //fevgen->SetMaxWeight( cm2/GeV );
      validcmd = true;
    }
    if( newValue.compareTo("gun") == 0 ){
      kinetemp = G4SBS::kGun;
      //   fevgen->SetKine( G4SBS::kGun );
      fevgen->SetRejectionSamplingFlag(false);
      validcmd = true;
    }
    if( newValue.compareTo("cosmics") == 0 ){
      kinetemp = G4SBS::kCosmics;
      //     fevgen->SetKine( G4SBS::kCosmics );
      validcmd = true;
    }

    if( newValue.compareTo("pythia6") == 0 ){
      kinetemp = G4SBS::kPYTHIA6;
      //fevgen->SetKine( G4SBS::kPYTHIA6 );
      fIO->SetUsePythia6( true );
      fevgen->SetRejectionSamplingFlag(false);
      validcmd = true;
    }
     if( newValue.compareTo("simc") == 0 ){
      kinetemp = G4SBS::kSIMC;
      fIO->SetUseSIMC( true );
      fevgen->SetRejectionSamplingFlag(false);
      validcmd = true;
    }
   if (newValue.compareTo("gmnelasticcheck") == 0 ){
      kinetemp = G4SBS::kGMnElasticCheck;
      //fevgen->SetKine(G4SBS::kGMnElasticCheck);
      fevgen->SetRejectionSamplingFlag(false);
      //fevgen->SetMaxWeight( cm2 );
      validcmd = true;
    }

    if( newValue.compareTo("wapp") == 0 ){ //wide angle pion photoproduction
      kinetemp = G4SBS::kPionPhoto;
      //fevgen->SetKine(G4SBS::kPionPhoto);
      validcmd = true;
    }

    if( !validcmd ){
      fprintf(stderr, "%s: %s line %d - Error: kinematic type %s not valid\n", __PRETTY_FUNCTION__, __FILE__, __LINE__, newValue.data());
      exit(1);
    } else { //valid kinematics given: 
      fevgen->SetKine(kinetemp);
      fIO->SetKine(kinetemp); //This is necessary because G4SBSIO and G4SBSEventGen cannot directly talk to each other
      G4SBSRun::GetRun()->GetData()->SetGenName(newValue.data());
    }

    //After any change of kinematics, set this flag to false so that rejection sampling will be re-initialized:
    fevgen->SetInitialized(false);
  }

  if( cmd == PYTHIAfileCmd ){
    fevgen->LoadPythiaChain( newValue );
  }

  if( cmd == SIMCfileCmd ){
    fevgen->LoadSIMCChain( newValue );
  }

  if( cmd == expCmd ){
    bool validcmd = false;
    if( newValue.compareTo("gep") == 0 ){
      fExpType = G4SBS::kGEp;
      validcmd = true;
    }

    if( newValue.compareTo("gep_bb") == 0 ){
      fExpType = G4SBS::kGEp_BB;
      validcmd = true;
    }
    
    if( newValue.compareTo("gepeplus") == 0 ){
      fExpType = G4SBS::kGEPpositron;
      validcmd = true;
    }
    if( newValue.compareTo("gmn") == 0 ){
      fExpType = G4SBS::kGMN;
      validcmd = true;
    }
    if( newValue.compareTo("gen") == 0 ){
      fExpType = G4SBS::kGEN;
      validcmd = true;
    }
    if( newValue.compareTo("genrp") == 0 ){
      fExpType = G4SBS::kGEnRP;
      validcmd = true;
    }
    if( newValue.compareTo("a1n") == 0 ){
      fExpType = G4SBS::kA1n; //"A1n" experiment type for new proposal with both SBS and BigBite in electron mode to detect DIS electrons at high-x: requires some geometry modifications on SBS side, including RICH w/CO2 instead of C4F10 and no aerogel, AND with a non-zero pitch angle for the SBS tracker. Also: HCAL + LAC.
      validcmd = true;
    }
    //AJP: Add SIDIS as a valid experiment type:
    if( newValue.compareTo("sidis") == 0 ){
      fExpType = G4SBS::kSIDISExp;
      validcmd = true;
    }
    if( newValue.compareTo("C16") == 0 ){
      fExpType = G4SBS::kC16;
      validcmd = true;
    }
    if( newValue.compareTo("tdis") == 0 ){
      fExpType = G4SBS::kTDIS;
      validcmd = true;
    }
    if( newValue.compareTo("ndvcs") == 0 ){
      fExpType = G4SBS::kNDVCS;
      validcmd = true;
    }
    if( newValue.compareTo("hcgem") == 0 ){
      fExpType = G4SBS::kGEMHCtest;
      validcmd = true;
    }
    if( newValue.compareTo("all") == 0 ){
      fExpType = G4SBS::kALL;
      validcmd = true;
    }

    if( validcmd ){
      fdetcon->SetExpType( fExpType );
      G4SBSRun::GetRun()->GetData()->SetExpType(newValue.data());
    }

    if( !validcmd ){
      fprintf(stderr, "%s: %s line %d - Error: kinematic type %s not valid\n", __PRETTY_FUNCTION__, __FILE__, __LINE__, newValue.data());
      exit(1);
    }
  }

  if( cmd == GunParticleCmd ){
    bool validcmd = false;
    fprigen->SetParticleName( newValue );
  }

  if( cmd == HadrCmd ){
    bool validcmd = false;
    if( newValue.compareTo("pi+") == 0 ){
      fevgen->SetHadronType( G4SBS::kPiPlus );
      validcmd = true;
    }
    if( newValue.compareTo("pi-") == 0 ){
      fevgen->SetHadronType( G4SBS::kPiMinus );
      validcmd = true;
    }
    if( newValue.compareTo("pi0") == 0 ){
      fevgen->SetHadronType( G4SBS::kPi0 );
      validcmd = true;
    }
    if( newValue.compareTo("K+") == 0 ){
      fevgen->SetHadronType( G4SBS::kKPlus );
      validcmd = true;
    }
    if( newValue.compareTo("K-") == 0 ){
      fevgen->SetHadronType( G4SBS::kKMinus );
      validcmd = true; 
    }
    if( newValue.compareTo("p") == 0 ){
      fevgen->SetHadronType( G4SBS::kP );
      validcmd = true;
    } 
    if( newValue.compareTo("pbar") == 0 ){
      fevgen->SetHadronType( G4SBS::kPbar );
      validcmd = true;
    }
    if( newValue.compareTo("n") == 0 ){
      fevgen->SetHadronType( G4SBS::kN );
      validcmd = true;
    }

    if( !validcmd ){
      fprintf(stderr, "%s: %s line %d - Error: Hadron type %s not valid\n", __PRETTY_FUNCTION__, __FILE__, __LINE__, newValue.data());
      exit(1);
    }

    fevgen->SetInitialized(false);

  }

  if( cmd == RejectionSamplingCmd ){ //Turn on rejection sampling for built-in event generators (has no effect for Pythia/beam/etc)
    std::istringstream is(newValue);
    G4bool flag;
    G4int N;
    if( newValue.contains("true") || newValue.contains("false") ){
      is >> std::boolalpha >> flag >> N;
    } else {
      is >> flag >> N;
    }
    
    fevgen->SetRejectionSamplingFlag(flag);
    fevgen->SetNeventsWeightCheck( N );
    fevgen->SetInitialized( false );
    //    if( flag ) fevgen->InitializeRejectionSampling();
  }
  
  if( cmd == tgtCmd ){
    bool validcmd = false;
    if( newValue.compareTo("LH2") == 0 ){
      fevgen->SetTarget(G4SBS::kLH2);
      fdetcon->SetTarget(G4SBS::kLH2);

      G4double den = (0.071*g/cm3)*Avogadro/(1.008*g/mole);
      fevgen->SetTargDen(den);
      fdetcon->fTargetBuilder->SetTargDen(den);
      validcmd = true;
    }
    if( newValue.compareTo("H2") == 0 ){
      fevgen->SetTarget(G4SBS::kH2);
      fdetcon->SetTarget(G4SBS::kH2);

      G4double den = 10.0*atmosphere/(296.0*kelvin*k_Boltzmann); //Should this be hard-coded? I think not. On the other hand, this provides a sensible default value, soooo....
      fevgen->SetTargDen(den);
      fdetcon->fTargetBuilder->SetTargDen(den);
      validcmd = true;
    }
    if( newValue.compareTo("D2") == 0 ){
      fevgen->SetTarget(G4SBS::kD2);
      fdetcon->SetTarget(G4SBS::kD2);

      G4double den = 1.0*atmosphere/(77.0*kelvin*k_Boltzmann);
      fevgen->SetTargDen(den);
      fdetcon->fTargetBuilder->SetTargDen(den);
      validcmd = true;
    }
    if( newValue.compareTo("LD2") == 0 ){
      fevgen->SetTarget(G4SBS::kLD2);
      fdetcon->SetTarget(G4SBS::kLD2);

      G4double den = (162.4*kg/m3)*Avogadro/(2.014*g/mole);
      fevgen->SetTargDen(den);
      fdetcon->fTargetBuilder->SetTargDen(den);
      validcmd = true;
    }
    if( newValue.compareTo("3He") == 0 ){
      fevgen->SetTarget(G4SBS::k3He);
      fdetcon->SetTarget(G4SBS::k3He);

      G4double den = 10.0*atmosphere/(296.0*kelvin*k_Boltzmann);
      fevgen->SetTargDen(den);
      fdetcon->fTargetBuilder->SetTargDen(den);
      validcmd = true;

    }
    if( newValue.compareTo("Neutron") == 0 ){
      fevgen->SetTarget(G4SBS::kNeutTarg);
      fdetcon->SetTarget(G4SBS::kNeutTarg);

      G4double den = 10.0*atmosphere/(296.0*kelvin*k_Boltzmann);
      fevgen->SetTargDen(den);
      fdetcon->fTargetBuilder->SetTargDen(den);
      validcmd = true;
    }
    if( newValue.compareTo("Cfoil") == 0 ){
      fevgen->SetTarget(G4SBS::kCfoil);
      fdetcon->SetTarget(G4SBS::kCfoil);
      validcmd = true;
    }

    if( newValue.compareTo("optics") == 0 ){
      fevgen->SetTarget(G4SBS::kOptics);
      fdetcon->SetTarget(G4SBS::kOptics);
      validcmd = true;
      //fdetcon->fTargetBuilder->SetNtargetFoils(1); //default to one carbon foil at Z = 0;
    }
    
    if( !validcmd ){
      fprintf(stderr, "%s: %s line %d - Error: target type %s not valid\n", __PRETTY_FUNCTION__, __FILE__, __LINE__, newValue.data());
      exit(1);
    }

    fevgen->SetInitialized(false);
    
  }

  // D. Flay (7/28/20) 
  // GEn 3He target Helmholtz coil configuration 
  if( cmd == GENTargetHelmholtzCmd ){
     G4int genTgtHHconf = GENTargetHelmholtzCmd->GetNewIntValue(newValue);  
     fdetcon->SetGEnTargetHelmholtzConfig(genTgtHHconf);
  }

  // D. Flay (9/29/20) 
  // GEn 3He target angular misalignment 
  if( cmd == GENTargetRXCmd ){
     G4double genTgtDRX = GENTargetRXCmd->GetNewDoubleValue(newValue);
     fdetcon->SetGEnTargetDRX(genTgtDRX); 
  }
  if( cmd == GENTargetRYCmd ){
     G4double genTgtDRY = GENTargetRYCmd->GetNewDoubleValue(newValue); 
     fdetcon->SetGEnTargetDRY(genTgtDRY); 
  }
  if( cmd == GENTargetRZCmd ){
     G4double genTgtDRZ = GENTargetRZCmd->GetNewDoubleValue(newValue); 
     fdetcon->SetGEnTargetDRZ(genTgtDRZ); 
  }

  // D. Flay (10/9/20) 
  // GEn 3He collimators 
  if( cmd == GENTargetColCmd ){
     G4bool tcEnable = GENTargetColCmd->GetNewBoolValue(newValue); 
     fdetcon->SetGEnTargetCollimatorEnable(tcEnable);
  }
  if( cmd == GENTargetColACmd ){
     G4bool tcaEnable = GENTargetColACmd->GetNewBoolValue(newValue); 
     fdetcon->SetGEnTargetCollimatorAEnable(tcaEnable);
  }
  if( cmd == GENTargetColBCmd ){
     G4bool tcbEnable = GENTargetColBCmd->GetNewBoolValue(newValue); 
     fdetcon->SetGEnTargetCollimatorBEnable(tcbEnable);
  }
  if( cmd == GENTargetColCCmd ){
     G4bool tccEnable = GENTargetColCCmd->GetNewBoolValue(newValue); 
     fdetcon->SetGEnTargetCollimatorCEnable(tccEnable);
  }

  if( cmd == GENTargetSDEnableCmd ){
     G4bool genSDEnable = GENTargetSDEnableCmd->GetNewBoolValue(newValue); 
     fdetcon->SetGEnTargetSDEnable(genSDEnable);  
  }

  // D. Flay (8/25/20) 
  // beam offset 
  if(cmd==beamOffsetXcmd){
     G4double bpx = beamOffsetXcmd->GetNewDoubleValue(newValue);
     fevgen->SetBeamOffsetX(bpx);
  }
  if(cmd==beamOffsetYcmd){
     G4double bpy = beamOffsetYcmd->GetNewDoubleValue(newValue);
     fevgen->SetBeamOffsetY(bpy);
  }
  // beam dump
  if(cmd==beamDumpCmd){
     G4bool bdEnable = beamDumpCmd->GetNewBoolValue(newValue); 
     fdetcon->SetBeamDumpEnable(bdEnable); 
  } 
  // beam diffuser
  if(cmd==beamDiffuserCmd){
     G4bool bdEnable = beamDiffuserCmd->GetNewBoolValue(newValue); 
     fdetcon->SetBeamDiffuserEnable(bdEnable); 
  } 

  if( cmd == bigfieldCmd ){
    G4int n = bigfieldCmd->GetNewIntValue(newValue);
    fdetcon->Set48D48Field(n);
  }

  // D. Flay (10/15/20) 
  // beam angular alignment  
  if(cmd==beamAngleXcmd){
     G4double bax = beamAngleXcmd->GetNewDoubleValue(newValue);
     fevgen->SetBeamAngleX(bax);
  }
  if(cmd==beamAngleYcmd){
     G4double bay = beamAngleYcmd->GetNewDoubleValue(newValue);
     fevgen->SetBeamAngleY(bay);
  }
  if(cmd==beamAngleZcmd){
     G4double baz = beamAngleZcmd->GetNewDoubleValue(newValue);
     fevgen->SetBeamAngleZ(baz);
  }

  // D. Flay (10/15/20) 
  // ion chamber enable  
  if( cmd == ionChamberEnableCmd ){ 
     G4bool icEnable = ionChamberEnableCmd->GetNewBoolValue(newValue);
     fdetcon->SetIonChamberEnable(icEnable);  
  }
  if( cmd == ionChamberXCmd ){ 
     G4double icx = ionChamberXCmd->GetNewDoubleValue(newValue);
     fdetcon->SetIonChamberX(icx);  
  }
  if( cmd == ionChamberYCmd ){ 
     G4double icy = ionChamberYCmd->GetNewDoubleValue(newValue);
     fdetcon->SetIonChamberY(icy);  
  }
  if( cmd == ionChamberZCmd ){ 
     G4double icz = ionChamberZCmd->GetNewDoubleValue(newValue);
     fdetcon->SetIonChamberZ(icz);  
  }
  if( cmd == ionChamberRXCmd ){ 
     G4double icrx = ionChamberRXCmd->GetNewDoubleValue(newValue);
     fdetcon->SetIonChamberRX(icrx);  
  }
  if( cmd == ionChamberRYCmd ){ 
     G4double icry = ionChamberRYCmd->GetNewDoubleValue(newValue);
     fdetcon->SetIonChamberRY(icry);  
  }
  if( cmd == ionChamberRZCmd ){ 
     G4double icrz = ionChamberRZCmd->GetNewDoubleValue(newValue);
     fdetcon->SetIonChamberRZ(icrz);  
  }

  // D. Flay (11/5/20) 
  // [GEn target] beam collimator enable  
  if( cmd == beamCollimatorEnableDnCmd ){ 
     G4bool bcEnable_dn = beamCollimatorEnableDnCmd->GetNewBoolValue(newValue);
     fdetcon->SetBeamCollimatorEnable_dnstr(bcEnable_dn);  
  }
  if( cmd == beamCollimatorLDnCmd ){ 
     G4double bcl_dn = beamCollimatorLDnCmd->GetNewDoubleValue(newValue);
     fdetcon->SetBeamCollimatorL_dnstr(bcl_dn);  
  }
  if( cmd == beamCollimatorDminDnCmd ){ 
     G4double bcdmin_dn = beamCollimatorDminDnCmd->GetNewDoubleValue(newValue);
     fdetcon->SetBeamCollimatorDmin_dnstr(bcdmin_dn);  
  }
  if( cmd == beamCollimatorDmaxDnCmd ){ 
     G4double bcdmax_dn = beamCollimatorDmaxDnCmd->GetNewDoubleValue(newValue);
     fdetcon->SetBeamCollimatorDmax_dnstr(bcdmax_dn);  
  }
  if( cmd == beamCollimatorXDnCmd ){ 
     G4double bcx_dn = beamCollimatorXDnCmd->GetNewDoubleValue(newValue);
     fdetcon->SetBeamCollimatorX_dnstr(bcx_dn);  
  }
  if( cmd == beamCollimatorYDnCmd ){ 
     G4double bcy_dn = beamCollimatorYDnCmd->GetNewDoubleValue(newValue);
     fdetcon->SetBeamCollimatorY_dnstr(bcy_dn);  
  }
  if( cmd == beamCollimatorZDnCmd ){ 
     G4double bcz_dn = beamCollimatorZDnCmd->GetNewDoubleValue(newValue);
     fdetcon->SetBeamCollimatorZ_dnstr(bcz_dn);  
  }
  // upstream
  if( cmd == beamCollimatorEnableUpCmd ){ 
     G4bool bcEnable_up = beamCollimatorEnableUpCmd->GetNewBoolValue(newValue);
     fdetcon->SetBeamCollimatorEnable_upstr(bcEnable_up);  
  }
  if( cmd == beamCollimatorLUpCmd ){ 
     G4double bcl_up = beamCollimatorLUpCmd->GetNewDoubleValue(newValue);
     fdetcon->SetBeamCollimatorL_upstr(bcl_up);  
  }
  if( cmd == beamCollimatorDminUpCmd ){ 
     G4double bcdmin_up = beamCollimatorDminUpCmd->GetNewDoubleValue(newValue);
     fdetcon->SetBeamCollimatorDmin_upstr(bcdmin_up);  
  }
  if( cmd == beamCollimatorDmaxUpCmd ){ 
     G4double bcdmax_up = beamCollimatorDmaxUpCmd->GetNewDoubleValue(newValue);
     fdetcon->SetBeamCollimatorDmax_upstr(bcdmax_up);  
  }
  if( cmd == beamCollimatorXUpCmd ){ 
     G4double bcx_up = beamCollimatorXUpCmd->GetNewDoubleValue(newValue);
     fdetcon->SetBeamCollimatorX_upstr(bcx_up);  
  }
  if( cmd == beamCollimatorYDnCmd ){ 
     G4double bcy_up = beamCollimatorYUpCmd->GetNewDoubleValue(newValue);
     fdetcon->SetBeamCollimatorY_upstr(bcy_up);  
  }
  if( cmd == beamCollimatorZUpCmd ){ 
     G4double bcz_up = beamCollimatorZUpCmd->GetNewDoubleValue(newValue);
     fdetcon->SetBeamCollimatorZ_upstr(bcz_up);  
  }

  if( cmd == bbfieldCmd ){
    std::istringstream is(newValue);

    //G4int n = bbfieldCmd->GetNewIntValue(newValue);
    G4int n;
    G4String fname;

    is >> n >> fname;
    
    fdetcon->SetBigBiteField(n, fname);
  }
  
  if( cmd == tosfieldCmd ){
    std::istringstream is(newValue);

    G4String fname;
    G4int flag;
    is >> fname >> flag;
    fdetcon->AddToscaField( fname.data(), flag );
    //    fdetcon->AddToscaField(newValue.data());
  }

  if( cmd == eventStatusEveryCmd ){
    fevact->SetEventStatusEvery(eventStatusEveryCmd->GetNewIntValue(newValue));
  }

  if( cmd == geantinoCmd ){
    G4bool b = geantinoCmd->GetNewBoolValue(newValue);
    fprigen->SetUseGeantino(b);
    fdetcon->GetGlobalField()->SetInvertField(b);
  }

  if( cmd == invertCmd ){
    G4bool b = invertCmd->GetNewBoolValue(newValue);
    fdetcon->GetGlobalField()->SetInvertField(b);
  }

  if( cmd == totalabsCmd ){
    G4bool b = totalabsCmd->GetNewBoolValue(newValue);
    fdetcon->SetTotalAbs(b);
  }

  if( cmd == checkOverlapCmd ){
    G4bool b = checkOverlapCmd->GetNewBoolValue(newValue);
    fdetcon->SetCheckOverlap(b);
  }

  if( cmd == tgtLenCmd ){
    G4double len = tgtLenCmd->GetNewDoubleValue(newValue);
    fevgen->SetTargLen(len);
    fdetcon->fTargetBuilder->SetTargLen(len);
    fevgen->SetInitialized(false);
  }

  if( cmd == tgtDenCmd ){
    G4double den = tgtDenCmd->GetNewDoubleValue(newValue);
    fevgen->SetTargDen(den);
    fdetcon->fTargetBuilder->SetTargDen(den);
    fevgen->SetInitialized(false);
  }
  if( cmd == tgtPresCmd ){
    G4double pre = tgtPresCmd->GetNewDoubleValue(newValue);
    G4double den = pre/(296.0*kelvin*k_Boltzmann); //molecules/unit volume
    fevgen->SetTargDen(den);
    fdetcon->fTargetBuilder->SetTargDen(den);
    fevgen->SetInitialized(false);
  }

  if( cmd == tgtDiamCmd ){
    G4double Dcell = tgtDiamCmd->GetNewDoubleValue(newValue);
    fdetcon->fTargetBuilder->SetTargDiameter(Dcell);
  }

  // if( cmd == SchamGasTgtCmd ){
  //   G4int flag = SchamGasTgtCmd->GetNewIntValue( newValue );
  //   fdetcon->fTargetBuilder->SetSchamFlag( flag );
  // }

  if( cmd == beamcurCmd ){
    G4double v = beamcurCmd->GetNewDoubleValue(newValue);
    printf("Setting beam current to %f uA\n", v/microampere);
    fevgen->SetBeamCur(v);
    fIO->SetBeamCur(v/microampere);
    fevgen->SetInitialized(false);
  }
  if( cmd == runtimeCmd ){
    G4double v = runtimeCmd->GetNewDoubleValue(newValue);
    fevgen->SetRunTime(v);
    fevgen->SetInitialized(false);
  }

  if( cmd == rasterxCmd ){
    G4double v = rasterxCmd->GetNewDoubleValue(newValue);
    fevgen->SetRasterX(v);
  }

  if( cmd == rasteryCmd ){
    G4double v = rasteryCmd->GetNewDoubleValue(newValue); 
    fevgen->SetRasterY(v);
  }

  if( cmd == rasterrCmd ){
    G4double v = rasterrCmd->GetNewDoubleValue(newValue);
    fevgen->SetRasterRadius(v);
  }
  
  if( cmd == beamspotsizeCmd ){
    G4double v = beamspotsizeCmd->GetNewDoubleValue(newValue);
    fevgen->SetBeamSpotSize(v);
  }
  
  if( cmd == tgtNfoilCmd ){
    G4int n = tgtNfoilCmd->GetNewIntValue(newValue);
    fdetcon->fTargetBuilder->SetNtargetFoils(n);
  }

  if( cmd == tgtFoilThickCmd ){
    std::istringstream is(newValue);
    G4int nfoil = fdetcon->fTargetBuilder->GetNtargetFoils();

    std::vector<double> thicktemp(nfoil);

    G4String unit;

    bool success = true;
    
    for( G4int ifoil=0; ifoil<nfoil; ifoil++ ){
      is >> thicktemp[ifoil];
      if( is.eof() || is.fail() || is.bad() ) {
	success = false;
	break;
      }
    }
    is >> unit;
    if( is.fail() || is.bad() || !is.eof() ) {
      success = false;
      exit(-1);
    }
    
    for( G4int ifoil=0; ifoil<nfoil; ifoil++ ){
      fdetcon->fTargetBuilder->SetFoilThick( ifoil, thicktemp[ifoil]*cmd->ValueOf(unit) );
    }
  }

  if( cmd == tgtFoilZCmd ){
    std::istringstream is(newValue);
    G4int nfoil = fdetcon->fTargetBuilder->GetNtargetFoils();

    std::vector<double> Ztemp(nfoil);

    G4String unit;

    bool success = true;
    
    for( G4int ifoil=0; ifoil<nfoil; ifoil++ ){
      is >> Ztemp[ifoil];
      if( is.eof() || is.fail() || is.bad() ) {
	success = false;
	break;
      }
    }
    is >> unit;
    if( is.fail() || is.bad() || !is.eof() ) {
      success = false;
      exit(-1);
    }
    
    for( G4int ifoil=0; ifoil<nfoil; ifoil++ ){
      fdetcon->fTargetBuilder->SetFoilZpos( ifoil, Ztemp[ifoil]*cmd->ValueOf(unit) );
    }
  }
    
  
  if( cmd == beamECmd ){
    G4double v = beamECmd->GetNewDoubleValue(newValue);
    fevgen->SetBeamE(v);
    fIO->SetBeamE(v);

    //    G4SBSRun::GetRun()->GetData()->SetBeamE(v/GeV); //redundant with fIO
    //after any command affecting the kinematics or cross section of the built-in event generators, re-initialize rejection sampling:
    fevgen->SetInitialized(false);
  }

  if( cmd == bbangCmd ){
    G4double v = bbangCmd->GetNewDoubleValue(newValue);
    printf("Setting BB ang to %f deg\n", v/deg);
    fdetcon->SetBBAng(v);
    fIO->SetBigBiteTheta(v);

    //If TOSCA map override flag is set, update the angle and distance for the field map:
    if( fdetcon->fGlobalField->GetOverride_Earm() ){
      fdetcon->fGlobalField->SetAngleAndDistance( fdetcon->fEArmBuilder->fBBang, fdetcon->fEArmBuilder->fBBdist, G4SBS::kEarm );
    }
  }

  if( cmd == bbdistCmd ){
    G4double v = bbdistCmd->GetNewDoubleValue(newValue);
    fdetcon->SetBBDist(v);
    fIO->SetBigBiteDist(v);

    if( fdetcon->fGlobalField->GetOverride_Earm() ){
      fdetcon->fGlobalField->SetAngleAndDistance( fdetcon->fEArmBuilder->fBBang, fdetcon->fEArmBuilder->fBBdist, G4SBS::kEarm );
    }
    
  }

  if( cmd == hcalangCmd ){
    G4double v = hcalangCmd->GetNewDoubleValue(newValue);
    fdetcon->Set48D48Ang(v);
    fIO->SetSBSTheta(v);

    if( fdetcon->fGlobalField->GetOverride_Harm() ){
      fdetcon->fGlobalField->SetAngleAndDistance( fdetcon->fHArmBuilder->f48D48ang, fdetcon->fHArmBuilder->f48D48dist, G4SBS::kHarm );
    }
    
  }

  if( cmd == sbstrkrpitchCmd ){
    G4double v = sbstrkrpitchCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetTrackerPitch(v);
    fIO->SetSBStrkrPitch( v );
    //G4SBSRun::GetRun()->GetData()->SetSBSTrackerPitch( v );
  }

  if( cmd == sbstrkrdistCmd ){
    G4double d = sbstrkrdistCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetTrackerDist(d);
    fIO->SetSBStrkrDist( d );
    //G4SBSRun::GetRun()->GetData()->SetSBSTrackerDist( d );
  }
  
  if( cmd == dvcsecalmatCmd ){
    fdetcon->fEArmBuilder->SetDVCSECalMaterial(newValue);
  }

  if( cmd == GRINCH_gas_Cmd ){
    G4String gasname = newValue;
    
    gasname.toUpper();

    if( gasname.contains( "C4F10" ) ){
      gasname = "C4F10_gas";
    } else if( gasname.contains( "C4F8O" ) ){
      gasname = "C4F8O";
    } else if( gasname.contains( "CF4" ) ){
      gasname = "CF4_gas";
    } else if( gasname.contains( "SF6" ) ){
      gasname = "SF6_gas";
    } else if( gasname.contains( "CO2" ) ){
      gasname = "CO2";
    } else if( gasname.contains( "C4F8" ) ){
      gasname = "C4F8_gas";
    } else { //default to C4F10 if no valid name given:
      gasname = "C4F10_gas";
      G4cout << "WARNING: invalid GRINCH gas option, defaulting to C4F10" << G4endl;
    }
    
    G4cout << "GRINCH gas name = " << gasname << G4endl;
    
    fdetcon->fEArmBuilder->SetGRINCHgas( gasname );
  }

  if( cmd == RICH_gas_Cmd ){
    G4String gasname = newValue;
    
    gasname.toUpper();

    //G4cout << "gasname = " << gasname << G4endl;

    //G4cout << gasname.index( "C4F10" ) << G4endl;
    
    if( gasname.contains( "C4F10" ) ){
      gasname = "C4F10_gas";
    } else if( gasname.contains( "C4F8O" ) ){
      gasname = "C4F8O";
    } else if( gasname.contains( "CF4" ) ){
      gasname = "CF4_gas";
    } else if( gasname.contains( "SF6" ) ){
      gasname = "SF6_gas";
    } else if( gasname.contains( "CO2" ) ){
      gasname = "CO2";
    } else if( gasname.contains( "C4F8" ) ){
      gasname = "C4F8_gas";
    } else { //default to C4F10 if no valid name given:
      gasname = "C4F10_gas";
      G4cout << "WARNING: invalid RICH gas option, defaulting to C4F10" << G4endl;
    }

    G4cout << "/g4sbs/richgas invoked, setting RICH gas to " << gasname << G4endl;
    
    fdetcon->fHArmBuilder->SetRICHgas( gasname );
  }
  
  if( cmd == hcaldistCmd ){
    G4double v = hcaldistCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetHCALDist(v);
    fevgen->SetHCALDist(v);
    fIO->SetHcalDist(v);
  }

  if( cmd == hcalvoffsetCmd ){
    G4double v = hcalvoffsetCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetHCALVOffset(v);
    //fevgen->SetHCALDist(v);
    fIO->SetHcalVOffset(v);
  }

  if( cmd == hcalhoffsetCmd ){
    G4double v = hcalhoffsetCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetHCALHOffset(v);
    //fevgen->SetHCALDist(v);
    fIO->SetHcalHOffset(v);
  }

  if( cmd == hcalangoffsetCmd ){
    G4double v = hcalangoffsetCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetHCALAngOffset(v);
    //fevgen->SetHCALDist(v);
    fIO->SetHcalAngOffset(v);
  }

  if( cmd == CDetReadyCmd ){
    G4bool v = CDetReadyCmd->GetNewBoolValue(newValue);
    fdetcon->fHArmBuilder->SetCDetReady(v);
    //fIO->SetCDetReady(v);
  }

  if( cmd == lacdistCmd ){
    G4double v = lacdistCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetLACDist(v);
    fIO->SetLACDist( v );
  }

  if( cmd == lacvoffsetCmd ){
    G4double v = lacvoffsetCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetLACVOffset(v);
    fIO->SetLACVOffset( v );
  }

  if( cmd == lachoffsetCmd ){
    G4double v = lachoffsetCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetLACHOffset(v);
    fIO->SetLACHOffset( v );
  }

  if( cmd == hmagdistCmd ){
    G4double v = hmagdistCmd->GetNewDoubleValue(newValue);
    fdetcon->Set48D48Dist(v);
    fIO->SetSBSDist( v );

    if( fdetcon->fGlobalField->GetOverride_Harm() ){
      fdetcon->fGlobalField->SetAngleAndDistance( fdetcon->fHArmBuilder->f48D48ang, fdetcon->fHArmBuilder->f48D48dist, G4SBS::kHarm );
    }
  } 

  if( cmd == cerDepCmd ){
    G4double v = cerDepCmd->GetNewDoubleValue(newValue);
    fdetcon->fEArmBuilder->SetCerDepth(v);
  }

  if( cmd == cerDisCmd ){
    G4double v = cerDisCmd->GetNewDoubleValue(newValue);
    fdetcon->fEArmBuilder->SetCerDist(v);
  }

  if( cmd == gemSepCmd ){
    G4double v = gemSepCmd->GetNewDoubleValue(newValue);
    fdetcon->fEArmBuilder->SetGEMSep(v);
  }

  if( cmd == bbCalDistCmd ){
    G4double v = bbCalDistCmd->GetNewDoubleValue(newValue);
    fdetcon->fEArmBuilder->SetBBCalDist(v);
  }

  if( cmd == thminCmd ){
    G4double v = thminCmd->GetNewDoubleValue(newValue);
    fevgen->SetThMin(v);
    fevgen->SetInitialized(false);
  }
  if( cmd == thmaxCmd ){
    G4double v = thmaxCmd->GetNewDoubleValue(newValue);
    fevgen->SetThMax(v);
    fevgen->SetInitialized(false);
  }
  if( cmd == phminCmd ){
    G4double v = phminCmd->GetNewDoubleValue(newValue);
    fevgen->SetPhMin(v);
    fevgen->SetInitialized(false);
  }
  if( cmd == phmaxCmd ){
    G4double v = phmaxCmd->GetNewDoubleValue(newValue);
    fevgen->SetPhMax(v);
    fevgen->SetInitialized(false);
  }
  if( cmd == HthminCmd ){
    G4double v = HthminCmd->GetNewDoubleValue(newValue);
    fevgen->SetThMin_had(v);
    fevgen->SetInitialized(false);
  }
  if( cmd == HthmaxCmd ){
    G4double v = HthmaxCmd->GetNewDoubleValue(newValue);
    fevgen->SetThMax_had(v);
    fevgen->SetInitialized(false);
  }

  if( cmd == HphminCmd ){
    G4double v = HphminCmd->GetNewDoubleValue(newValue);
    fevgen->SetPhMin_had(v);
    fevgen->SetInitialized(false);
  }
  
  if( cmd == HphmaxCmd ){
    G4double v = HphmaxCmd->GetNewDoubleValue(newValue);
    fevgen->SetPhMax_had(v);
    fevgen->SetInitialized(false);
  }

  if( cmd == EhminCmd ){
    G4double v = EhminCmd->GetNewDoubleValue(newValue);
    fevgen->SetEhadMin(v);
    fevgen->SetInitialized(false);
  }
  if( cmd == EhmaxCmd ){
    G4double v = EhmaxCmd->GetNewDoubleValue(newValue);
    fevgen->SetEhadMax(v);
    fevgen->SetInitialized(false);
  }
  if( cmd == EeminCmd ){
    G4double v = EeminCmd->GetNewDoubleValue(newValue);
    fevgen->SetEeMin(v);
    fevgen->SetInitialized(false);
  }
  if( cmd == EemaxCmd ){
    G4double v = EemaxCmd->GetNewDoubleValue(newValue);
    fevgen->SetEeMax(v);
    fevgen->SetInitialized(false);
  }

  if( cmd == PionPhoto_tminCmd ){
    G4double v = PionPhoto_tminCmd->GetNewDoubleValue(newValue);
    fevgen->SetPionPhoto_tmin( v );
  }

  if( cmd == PionPhoto_tmaxCmd ){
    G4double v = PionPhoto_tmaxCmd->GetNewDoubleValue(newValue);
    fevgen->SetPionPhoto_tmax( v );
  }

  if( cmd == PionPhoto_useradCmd ){
    G4bool b = PionPhoto_useradCmd->GetNewBoolValue(newValue);
    fevgen->SetUseRadiator( b );
    fdetcon->fTargetBuilder->SetUseRad( b );
  }

  if( cmd == PionPhoto_radthickCmd ){
    G4double v = PionPhoto_radthickCmd->GetNewDoubleValue(newValue);
    fevgen->SetRadthickX0( v );
    fdetcon->fTargetBuilder->SetRadThick( v );
  }

  if ( cmd == PionPhoto_radzCmd ){
    G4double v = PionPhoto_radzCmd->GetNewDoubleValue(newValue);
    fdetcon->fTargetBuilder->SetRadZoffset( v );
  }
  
  if( cmd == gemresCmd ){
    G4double v = gemresCmd->GetNewDoubleValue(newValue);
    fevact->SetGEMRes(v);
  }

  if( cmd == RICHdistCmd ){
    G4double v = RICHdistCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetRICHdist(v);
    fIO->SetRICHDist( v );
  }

  if( cmd == RICHhoffsetCmd ){
    G4double v = RICHhoffsetCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetRICHHoffset( v );
  }

  if( cmd == RICHvoffsetCmd ){
    G4double v = RICHvoffsetCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetRICHVoffset( v );
  }

  if( cmd == RICHaeroCmd ){
    G4bool b = RICHaeroCmd->GetNewBoolValue(newValue);
    fdetcon->fHArmBuilder->SetRICH_use_aerogel( b );
  }
  
  if( cmd == RICHSnoutExtensionCmd ){
    G4double v = RICHSnoutExtensionCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetRICHSnoutExtension(v);
  }

  if( cmd == SBSMagFieldCmd ){
    G4double v = SBSMagFieldCmd->GetNewDoubleValue(newValue);
    fdetcon->SetUniformMagneticField48D48( v );
  }

  if( cmd == EARM_ScaleFieldCmd ){
    G4double v = EARM_ScaleFieldCmd->GetNewDoubleValue( newValue );
    fdetcon->SetFieldScale_BB( v );
    // if( fdetcon->GetBBField() != NULL ){
    //   fdetcon->GetBBField()->fScaleFactor = s;
    // }
  }

  if( cmd == HARM_ScaleFieldCmd ){
    G4double v = HARM_ScaleFieldCmd->GetNewDoubleValue( newValue );
    fdetcon->SetFieldScale_SBS(v);
    // if( fdetcon->Get48D48Field() != NULL ){
    //   fdetcon->Get48D48Field()->fScaleFactor = s;
    //   G4cout << "Setting SBS magnetic field scale factor to " << s << G4endl;
    // }

    // if( fdetcon->fUseGlobalField ){
      
    // }
  }
  
  if( cmd == SBSFieldClampOptionCmd ){
    G4int i = SBSFieldClampOptionCmd->GetNewIntValue(newValue);
    fdetcon->fHArmBuilder->SetFieldClampConfig48D48( i );
  }

  if( cmd == SBSBeamlineConfCmd ){
    G4int i = SBSBeamlineConfCmd->GetNewIntValue(newValue);
    fdetcon->fBeamlineConf = i;
  }
  
  if( cmd == SBSLeadOptionCmd ){
    G4int i = SBSLeadOptionCmd->GetNewIntValue(newValue);
    fdetcon->fLeadOption = i;
  }

  if( cmd == GENRPAnalyzerOptionCmd ){
    G4int i = GENRPAnalyzerOptionCmd->GetNewIntValue(newValue);
    fdetcon->fHArmBuilder->SetGENRPAnalyzerOption(i);
  }

  if( cmd == GEPFPP1_CH2thickCmd ){
    G4double ch2thick = GEPFPP1_CH2thickCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetFPP_CH2thick(1,ch2thick);
  }

  if( cmd == GEPFPP2_CH2thickCmd ){
    G4double ch2thick = GEPFPP2_CH2thickCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetFPP_CH2thick(2,ch2thick);
  }

  if( cmd == GEPFPPoptionCmd ){
    G4int i = GEPFPPoptionCmd->GetNewIntValue(newValue);
    fdetcon->fHArmBuilder->SetGEPFPPoption( i );
  }

  if( cmd == HadronFilterCmd ){
    G4bool flag = HadronFilterCmd->GetNewBoolValue(newValue);
    fdetcon->fTargetBuilder->EnableHadronFilter(flag);
  }

  if( cmd == HadronFilterThickCmd ){
    G4double shieldthick = HadronFilterThickCmd->GetNewDoubleValue(newValue);
    fdetcon->fTargetBuilder->SetHadronFilterThick( shieldthick );
  }

  if( cmd == HadronFilterMaterialCmd ){
    G4String matname = newValue;
    fdetcon->fTargetBuilder->SetHadronFilterMaterial( matname );
  }
  
  if( cmd == BLneutronDetsCmd ){
    G4bool v = BLneutronDetsCmd->GetNewBoolValue(newValue);
    fdetcon->fBLneutronDet = v;
  }
  
  if( cmd == GEMfrontendCmd ){
    G4bool v = GEMfrontendCmd->GetNewBoolValue(newValue);
    fdetcon->fEArmBuilder->SetGEMfrontend(v);
  }
  
  if( cmd == GEMfrontendDistCmd ){
    G4double v = GEMfrontendDistCmd->GetNewDoubleValue(newValue);
    fdetcon->fEArmBuilder->fGEMfrontendDist = v;
  }
  
  if( cmd == GEMfrontendPosAngleCmd ){
    G4double v = GEMfrontendPosAngleCmd->GetNewDoubleValue(newValue);
    fdetcon->fEArmBuilder->fGEMfrontendPosAngle = v;
  }
  
 if( cmd == GEMfrontendRotAngleCmd ){
    G4double v = GEMfrontendRotAngleCmd->GetNewDoubleValue(newValue);
    fdetcon->fEArmBuilder->fGEMfrontendRotAngle = v;
  }
  
  if( cmd == SetGrinchPMTglassHitsCmd ){
    G4bool v = SetGrinchPMTglassHitsCmd->GetNewBoolValue(newValue);
    fdetcon->fEArmBuilder->SetGrinchPMTglassHits(v);
  }
  
  if( cmd == buildSBSsieveCmd ){
    G4bool b = buildSBSsieveCmd->GetNewBoolValue(newValue);
    fdetcon->fHArmBuilder->SetSBSSieve(b);
  }
  //SSeeds - Multiple sieve plate option 10.4.20
  if( cmd == buildBBsieveCmd ){
    G4int bbsieveconfval = buildBBsieveCmd->GetNewIntValue(newValue);
    fdetcon->fEArmBuilder->SetBBSieve(bbsieveconfval);
  }
  
  if( cmd == TreeFlagCmd ){
    G4int flag = TreeFlagCmd->GetNewIntValue(newValue);
    fevact->SetTreeFlag( flag );
  }

  if( cmd == SBS_FT_absorberCmd ){
    G4bool flag = SBS_FT_absorberCmd->GetNewBoolValue(newValue);

    fdetcon->fHArmBuilder->SetFTuseabsorber( flag );
  }

  if( cmd == SBS_FT_absorberMaterialCmd ){
    G4String matname = newValue;

    fdetcon->fHArmBuilder->SetFTabsmaterial( matname );
  }

  if( cmd == SBS_FT_absorberThickCmd ){
    G4double absthick = SBS_FT_absorberThickCmd->GetNewDoubleValue(newValue);

    fdetcon->fHArmBuilder->SetFTabsthick( absthick );
  }
  
  // if( cmd == Earm_CAL_part_cmd ){ 
  //   G4bool flag = Earm_CAL_part_cmd->GetNewBoolValue( newValue );
  //   fIO->SetEarmCALpart_flag( flag );
  // }
    
  // if( cmd == Harm_CAL_part_cmd ){
  //   G4bool flag = Harm_CAL_part_cmd->GetNewBoolValue( newValue );
  //   fIO->SetHarmCALpart_flag( flag );
  // }

  if( cmd == KeepPartCALcmd ){
    std::istringstream is(newValue);

    G4bool flag;
    G4String SDname;

    //Let's do (somewhat) intelligent parsing of the string here:
    if( newValue.contains("true") || newValue.contains("false") ){ //parse with the "boolalpha" flag:
      is >> SDname >> std::boolalpha >> flag;
    } else { //assume that the boolean parameter is given as 1 or 0:
      is >> SDname >> flag;
    }

    fIO->KeepPartCALflags[SDname] = flag;
  }

  if( cmd == KeepHistorycmd ){
    std::istringstream is(newValue);
    G4bool flag;
    G4String SDname;

    //Let's do (somewhat) intelligent parsing of the string here:
    if( newValue.contains("true") || newValue.contains("false") ){ //parse with the "boolalpha" flag:
      is >> SDname >> std::boolalpha >> flag;
    } else { //assume that the boolean parameter is given as 1 or 0:
      is >> SDname >> flag;
    }
    //is >> SDname >> std::boolalpha >> flag;

    G4cout << "/g4sbs/keephistory invoked, (SDname, flag)=(" << SDname << ", " << flag << ")" << G4endl;
    
    fIO->KeepHistoryflags[SDname] = flag;
  }

  if( cmd == LimitStepCALcmd ){
    std::istringstream is(newValue);
    G4bool flag;
    G4String SDname;

    //Let's do (somewhat) intelligent parsing of the string here:
    if( newValue.contains("true") || newValue.contains("false") ){ //parse with the "boolalpha" flag:
      is >> SDname >> std::boolalpha >> flag;
    } else { //assume that the boolean parameter is given as 1 or 0:
      is >> SDname >> flag;
    }
    
    //is >> SDname >> std::boolalpha >> flag;

    G4cout << "/g4sbs/steplimit invoked" << G4endl;
    G4cout << "newValue = " << newValue << G4endl;
    G4cout << "SDname = " << SDname << G4endl;
    G4cout << "flag = " << flag << G4endl;
    
    if( flag ){
      (fdetcon->StepLimiterList).insert(SDname);
    } else {
      (fdetcon->StepLimiterList).erase(SDname);
    }
  }

  // ******
  if( cmd == SD_EnergyThresholdCmd ){ //store the SDname and dimensioned threshold value in a map<G4String,G4double> assigned to fdetcon?
    std::istringstream is(newValue);

    G4String SDname;
    G4double ethresh;
    G4String unit;

    is >> SDname >> ethresh >> unit;
    
    fdetcon->SDthreshold[SDname] = ethresh*cmd->ValueOf(unit);

    G4cout << "Set Energy threshold for SD name = " << SDname << " to " << fdetcon->SDthreshold[SDname]/MeV << " MeV" << G4endl;
    
  }

  if( cmd == SD_TimeWindowCmd ){ 
    std::istringstream is(newValue);

    G4String SDname;
    G4double timewindow;
    G4String unit;

    is >> SDname >> timewindow >> unit;
    
    fdetcon->SDgatewidth[SDname] = timewindow*cmd->ValueOf(unit);

    G4cout << "Set time window for SD name = " << SDname << " to " << fdetcon->SDgatewidth[SDname]/ns << " ns" << G4endl;
  }

  if( cmd == SD_NTimeBinsCmd ){ 
    std::istringstream is(newValue);

    G4String SDname;
    G4int ntimebins;
 
    is >> SDname >> ntimebins;
    
    fdetcon->SDntimebins[SDname] = ntimebins;

    G4cout << "Set number of time bins for SD name = " << SDname << " to " << fdetcon->SDntimebins[SDname] << G4endl;
  }

  if( cmd == KeepPulseShapeCmd ){ //

    std::istringstream is(newValue);

    G4String SDname;
    G4bool flag;

    is >> SDname;
    
    //Let's do (somewhat) intelligent parsing of the string here:
    if( newValue.contains("true") || newValue.contains("false") ){ //parse with the "boolalpha" flag:
      is >> std::boolalpha >> flag;
    } else { //assume that the boolean parameter is given as 1 or 0:
      is >> flag;
    }
    
    //is >> SDname >> flag; 
    fIO->SetKeepPulseShape( SDname, flag );
    if( SDname == "all" ) fIO->SetKeepAllPulseShape(flag);
   
  }
  // ******

  if( cmd == KeepSDtrackcmd ){ //
    //newValue.toLower();
    
    std::istringstream is(newValue);

    G4String SDname;
    G4bool flag;

    is >> SDname;
    
    //Let's do (somewhat) intelligent parsing of the string here:
    if( newValue.contains("true") || newValue.contains("false") ){ //parse with the "boolalpha" flag:
      is >> std::boolalpha >> flag;
    } else { //assume that the boolean parameter is given as 1 or 0:
      is >> flag;
    }
    
    //is >> SDname >> flag; 
    fIO->SetKeepSDtracks( SDname, flag );
    if( SDname == "all" ) fIO->SetKeepAllSDtracks(flag);
    
  }
  // if( cmd == DisableOpticalPhysicsCmd ){
  //   G4bool b = DisableOpticalPhysicsCmd->GetNewBoolValue(newValue);
  //   if( b ){ 
  //     //if( fphyslist->GetOpticalPhysics() != NULL ) fphyslist->RemovePhysics( fphyslist->GetOpticalPhysics() );
  //     //fphyslist->SetOpticalPhysics( new G4OpticalPhysics(0) );
  //     //G4VPhysicsConstructor *ctemp;
  //     fphyslist->ReplacePhysics( new G4OpticalPhysics(0) );
  //     //fphyslist->SetOpticalPhysics(ctemp);
  //   } else {
  //     G4VPhysicsConstructor *ctemp = new G4OpticalPhysics(0);
  //     fphyslist->RemovePhysics( ctemp );
  //     delete ctemp;
  //   }
  //   G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  // }

  if( cmd == FluxCmd ){
    G4bool b = FluxCmd->GetNewBoolValue(newValue);
    fdetcon->fTargetBuilder->SetFlux(b);
    fIO->KeepPartCALflags["FLUX"] = b;
  }

  if( cmd == TargPolDirectionCmd ){
    G4ThreeVector v = TargPolDirectionCmd->GetNew3VectorValue(newValue);
    fdetcon->fTargetBuilder->SetTargPolDir( v.unit() );
    fevgen->SetTargPolDir( v.unit() );
  }

  if( cmd == TargPolMagnitudeCmd ){
    G4double pol = TargPolMagnitudeCmd->GetNewDoubleValue( newValue );
    fdetcon->fTargetBuilder->SetTargPolMag( pol );
    fevgen->SetTargPolMag( pol );
  }

  if( cmd == BeamPolDirectionCmd ){
    G4ThreeVector v = BeamPolDirectionCmd->GetNew3VectorValue(newValue);
    //fdetcon->fTargetBuilder->SetTargPolDir( v.unit() );
    fevgen->SetBeamPolDir( v.unit() );
  }

  if( cmd == BeamPolMagnitudeCmd ){
    G4double pol = BeamPolMagnitudeCmd->GetNewDoubleValue( newValue );
    //fdetcon->fTargetBuilder->SetTargPolMag( pol );
    fevgen->SetBeamPolMag( pol );
  }

  if( cmd == RandomizeTargetSpinCmd ){
    G4bool flag = RandomizeTargetSpinCmd->GetNewBoolValue( newValue );
    fevgen->SetRandomizeTargetSpin( flag );
  }

  if( cmd == NumSpinStatesTargCmd ){
    G4int nspin = NumSpinStatesTargCmd->GetNewIntValue( newValue );
    fevgen->SetNumTargetSpinDirections( nspin );
  }

  if( cmd == TargThetaSpinCmd ){ //there must be at least nspin entries or this will fail:
    std::istringstream is(newValue);
    G4int nspin = fevgen->GetNumTargetSpinDirections();

    std::vector<G4double> thspintemp(nspin);

    G4String unit;

    bool success = true;
    for( G4int ispin=0; ispin<nspin; ispin++ ){
      is >> thspintemp[ispin];
      if( is.eof() || is.fail() || is.bad() ){
	success = false;
	exit(-1);
      }
    }

    
    is >> unit; 
    if( is.fail() || is.bad() || !is.eof() ) {
      success = false;
      exit(-1);
    }

    for( G4int ispin=0; ispin<nspin; ispin++ ){
      fevgen->SetTargetThetaSpin( ispin, thspintemp[ispin]*cmd->ValueOf(unit) );
    }
  }

  if( cmd == TargPhiSpinCmd ){ //there must be at least nspin entries or this will fail:
    std::istringstream is(newValue);
    G4int nspin = fevgen->GetNumTargetSpinDirections();

    std::vector<G4double> phspintemp(nspin);

    G4String unit;

    bool success = true;
    for( G4int ispin=0; ispin<nspin; ispin++ ){
      is >> phspintemp[ispin];
      if( is.eof() || is.fail() || is.bad() ){
	success = false;
	exit(-1);
      }
    }

    
    is >> unit; 
    if( is.fail() || is.bad() || !is.eof() ) {
      success = false;
      exit(-1);
    }

    for( G4int ispin=0; ispin<nspin; ispin++ ){
      fevgen->SetTargetPhiSpin( ispin, phspintemp[ispin]*cmd->ValueOf(unit) );
    }
  }
  
  if( cmd == UseCerenkovCmd ){
    G4bool b = UseCerenkovCmd->GetNewBoolValue(newValue);
    fphyslist->ToggleCerenkov(b);
    fIO->SetUsingCerenkov(b);
  }

  if( cmd == UseScintCmd ){
    G4bool b = UseScintCmd->GetNewBoolValue(newValue);
    fphyslist->ToggleScintillation(b);
    fIO->SetUsingScintillation(b);
  }

  if( cmd == DisableOpticalPhotonProductionByMaterialCmd ){
    G4String materialname = newValue;
    fdetcon->SetOpticalPhotonDisabled( materialname );
  }

  if( cmd == GunPolarizationCommand ){
    G4ThreeVector pol = GunPolarizationCommand->GetNew3VectorValue(newValue);
    fprigen->SetGunPolarization( pol.unit() );
  }

  if( cmd == CosmicsPointerCommand ){
    G4ThreeVector point = CosmicsPointerCommand->GetNew3VectorValue(newValue);
    fevgen->SetCosmicsPointer( point );
    fevgen->UpdateCosmicsCeilingRadius();
  }

  if( cmd == CosmicsPointerRadiusCommand ){
    G4double radius = CosmicsPointerRadiusCommand->GetNewDoubleValue(newValue);
    fevgen->SetCosmicsPointerRadius( radius );
    fevgen->UpdateCosmicsCeilingRadius();
  }
  
  if( cmd == CosmicsMaxAngleCommand ){
    G4double maxangle = CosmicsMaxAngleCommand->GetNewDoubleValue(newValue);
    fevgen->SetCosmicsMaxAngle( maxangle );
  }

  if( cmd == WriteFieldMapCmd ){
    G4bool flag = WriteFieldMapCmd->GetNewBoolValue(newValue);
    fIO->SetWriteFieldMaps(flag);
  }

  if( cmd == UseGEMshieldCmd ){
    G4bool flag = UseGEMshieldCmd->GetNewBoolValue(newValue);
    fdetcon->SetGEMuseAlshield(flag);
  }

  if( cmd == GEMshieldThickCmd ){
    G4double shieldthick = GEMshieldThickCmd->GetNewDoubleValue(newValue);
    fdetcon->SetGEMAlShieldThick( shieldthick );
  }

  if( cmd == GEMshieldAirGapThickCmd ){
    G4double airgapthick = GEMshieldAirGapThickCmd->GetNewDoubleValue(newValue);
    fdetcon->SetGEMAirGapThick( airgapthick );
  }
  
  if( cmd == EnableBigBitePlateCmd ){
    G4bool flag = EnableBigBitePlateCmd->GetNewBoolValue(newValue);
    fdetcon->fTargetBuilder->EnableBigBitePlate(flag);
  }
  
  if( cmd == SetBigBitePlateThicknessCmd ){
    G4double platethick = SetBigBitePlateThicknessCmd->GetNewDoubleValue(newValue);
    fdetcon->fTargetBuilder->SetBigBitePlateThickness(platethick);
  }
  
  if( cmd == SetBigBitePlateMaterialCmd ){
    G4String platematerial = newValue;
    fdetcon->fTargetBuilder->SetBigBitePlateMaterial(platematerial);
  }
  
  
}
