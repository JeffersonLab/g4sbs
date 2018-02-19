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
  fExpType = kNeutronExp;

  runCmd = new G4UIcmdWithAnInteger("/g4sbs/run",this);
  runCmd->SetGuidance("Run simulation with x events");
  runCmd->SetParameterName("nevt", false);

  gemconfigCmd = new G4UIcmdWithAnInteger("/g4sbs/gemconfig",this);
  gemconfigCmd->SetGuidance("BigBite GEM layout: option 1 (default), 2 or 3");
  gemconfigCmd->SetParameterName("gemconfig", false);

  shieldconfigCmd = new G4UIcmdWithAnInteger("/g4sbs/bbshieldconfig",this);
  shieldconfigCmd->SetGuidance("BB Ecal shielding layout: option 0 (none), 1 (default), 2 (+10cm Al + 3cm steel on side), 3 (+3cm steel + 3cm steel on side)");
  shieldconfigCmd->SetParameterName("bbshieldconfig", false);

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

  kineCmd = new G4UIcmdWithAString("/g4sbs/kine",this);
  kineCmd->SetGuidance("Kinematics from elastic, inelastic, flat, dis, beam, sidis, wiser, gun, pythia6");
  kineCmd->SetParameterName("kinetype", false);

  PYTHIAfileCmd = new G4UIcmdWithAString("/g4sbs/pythia6file",this);
  PYTHIAfileCmd->SetGuidance("Name of ROOT file containing PYTHIA6 events as a ROOT tree");
  PYTHIAfileCmd->SetParameterName("fname",false);
  
  expCmd = new G4UIcmdWithAString("/g4sbs/exp",this);
  expCmd->SetGuidance("Experiment type from gep, gmn, gen, a1n, sidis, C16, tdis, ndvcs");
  expCmd->SetParameterName("exptype", false);

  GunParticleCmd = new G4UIcmdWithAString("/g4sbs/particle",this);
  GunParticleCmd->SetGuidance("Particle type for gun generator (valid GEANT4 particle names)");
  GunParticleCmd->SetParameterName("ptype", false );

  HadrCmd = new G4UIcmdWithAString("/g4sbs/hadron",this);
  HadrCmd->SetGuidance("Hadron type h for SIDIS N(e,e'h)X generator: pi+/pi-/K+/K-/p/pbar possible");
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

  bbfieldCmd = new G4UIcmdWithAnInteger("/g4sbs/bbfield", this);
  bbfieldCmd->SetGuidance("Turn on Bigbite field (requires field map)");
  bbfieldCmd->SetParameterName("bbfield", false);

  tosfieldCmd = new G4UIcmdWithAString("/g4sbs/tosfield", this);
  tosfieldCmd->SetGuidance("Use SBS TOSCA field map from file");
  tosfieldCmd->SetParameterName("tosfield", false);

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
  rasterxCmd->SetParameterName("size", false);

  rasteryCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/rastery",this);
  rasteryCmd->SetGuidance("Raster y size");
  rasteryCmd->SetParameterName("size", false);

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
  
  dvcsecalmatCmd = new G4UIcmdWithAString("/g4sbs/dvcsecalmat",this);
  dvcsecalmatCmd->SetGuidance("DVCS ECal material: 'PbF2' or 'PbWO4'");
  dvcsecalmatCmd->SetParameterName("dvcsecalmatname", false);
  
  hcaldistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hcaldist",this);
  hcaldistCmd->SetGuidance("HCAL distance");
  hcaldistCmd->SetParameterName("dist", false);

  hcalvoffsetCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hcalvoffset",this);
  hcalvoffsetCmd->SetGuidance("HCAL vertical offset");
  hcalvoffsetCmd->SetParameterName("dist", false);

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

  RICHdistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/richdist",this);
  RICHdistCmd->SetGuidance("SBS RICH distance from target");
  RICHdistCmd->SetParameterName("dist",false);

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
  SBSFieldClampOptionCmd->SetGuidance("SBS field clamp configuration: 0=no clamp, 1=BigBite(default), 2=GEp");
  SBSFieldClampOptionCmd->SetParameterName("sbsclampoption",false);

  SBSBeamlineConfCmd = new G4UIcmdWithAnInteger("/g4sbs/beamlineconfig",this);
  SBSBeamlineConfCmd->SetGuidance("SBS beamline configuration: 1: GEp; 2: GEn; 3 (def): GMn 2-4 pass; 4; GMn 5 pass");
  SBSBeamlineConfCmd->SetParameterName("beamlineconf",false);
  
  SBSLeadOptionCmd = new G4UIcmdWithAnInteger("/g4sbs/uselead",this);
  SBSLeadOptionCmd->SetGuidance("SBS beamline lead shielding configuration: 0= nope 1=yes");
  SBSLeadOptionCmd->SetParameterName("uselead",false);

  buildSBSsieveCmd = new G4UIcmdWithABool("/g4sbs/buildSBSsieve",this);
  buildSBSsieveCmd->SetGuidance("Use SBS sieve (true or false, false by default)");
  buildSBSsieveCmd->SetParameterName("buildSBSsieve",false);

  buildBBsieveCmd = new G4UIcmdWithABool("/g4sbs/buildBBsieve",this);
  buildBBsieveCmd->SetGuidance("Use BB sieve (true or false, false by default)");
  buildBBsieveCmd->SetParameterName("buildBBsieve",false);
  
  TreeFlagCmd = new G4UIcmdWithAnInteger("/g4sbs/treeflag",this);
  TreeFlagCmd->SetGuidance("G4SBS ROOT tree filling: 0=keep all, 1=keep only evts w/hits in sensitive volumes");
  TreeFlagCmd->SetParameterName("treeflag",false);

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

  FluxCmd = new G4UIcmdWithABool("/g4sbs/fluxcalc",this);
  FluxCmd->SetGuidance( "Compute particle flux as a function of angles, energy");
  FluxCmd->SetParameterName( "fluxcalc", false);  
  
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
}

G4SBSMessenger::~G4SBSMessenger(){
}


void G4SBSMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue){
  char cmdstr[255];

  if( cmd == runCmd ){
	
    G4VPhysicalVolume* pWorld;

    G4int nevt = runCmd->GetNewIntValue(newValue);

    //If the generator is PYTHIA, don't try to generate more events than we have available:
    if( fevgen->GetKine() == kPYTHIA6 && fevgen->GetPythiaChain()->GetEntries() < nevt ){
      nevt = fevgen->GetPythiaChain()->GetEntries();
    }
    
    fevgen->SetNevents(nevt);
    fevgen->Initialize();
    
    // if( fevgen->GetRejectionSamplingFlag() && !fevgen->GetRejectionSamplingInitialized() ){
    //   fevgen->InitializeRejectionSampling();
    // }

    Kine_t kinetype = fevgen->GetKine(); 
    if( kinetype == kDIS || kinetype == kWiser ){ //Processes with xsec in units of area/energy/solid angle; i.e., nb/GeV/sr
      G4SBSRun::GetRun()->GetData()->SetGenVol( fevgen->GetGenVol()/GeV );
      //if( fevgen->GetRejectionSamplingFlag() ){
      G4SBSRun::GetRun()->GetData()->SetMaxWeight( fevgen->GetMaxWeight()/cm2 * GeV );
	//}
    } else if ( kinetype == kSIDIS ){ //Processes with xsec differential in area/energy^2/solid angle^2:
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
    //fIO->InitializeTree();

    G4cout << "InitializeTree() successful" << G4endl;

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
    if( newValue.compareTo("elastic") == 0 ){
      fevgen->SetKine(kElastic);
      validcmd = true;
    }
    if( newValue.compareTo("inelastic") == 0 ){
      fevgen->SetKine(kInelastic);
      //fevgen->SetRejectionSamplingFlag(false);
      validcmd = true;
    }
    if( newValue.compareTo("flat") == 0 ){
      fevgen->SetKine(kFlat);
      //fevgen->SetMaxWeight( cm2 );
      fevgen->SetRejectionSamplingFlag(false);
      validcmd = true;
    }
    if( newValue.compareTo("dis") == 0 ){
      fevgen->SetKine(kDIS);
      //fevgen->SetMaxWeight( cm2/GeV );
      validcmd = true;
    }
    if( newValue.compareTo("beam") == 0 ){
      fevgen->SetKine(kBeam);
      fevgen->SetRejectionSamplingFlag(false);
      validcmd = true;
    }
    if( newValue.compareTo("sidis") == 0 ){
      fevgen->SetKine( kSIDIS );
      //fevgen->SetMaxWeight( cm2/pow(GeV,2) );
      validcmd = true;
    }
    if( newValue.compareTo("wiser") == 0 ){
      fevgen->SetKine( kWiser);
      //fevgen->SetMaxWeight( cm2/GeV );
      validcmd = true;
    }
    if( newValue.compareTo("gun") == 0 ){
      fevgen->SetKine( kGun );
      fevgen->SetRejectionSamplingFlag(false);
      validcmd = true;
    }

    if( newValue.compareTo("pythia6") == 0 ){
      fevgen->SetKine( kPYTHIA6 );
      fIO->SetUsePythia6( true );
      fevgen->SetRejectionSamplingFlag(false);
      validcmd = true;
    }
    if (newValue.compareTo("gmnelasticcheck") == 0 ){
      fevgen->SetKine(kGMnElasticCheck);
      fevgen->SetRejectionSamplingFlag(false);
      //fevgen->SetMaxWeight( cm2 );
      validcmd = true;
    }

    if( !validcmd ){
      fprintf(stderr, "%s: %s line %d - Error: kinematic type %s not valid\n", __PRETTY_FUNCTION__, __FILE__, __LINE__, newValue.data());
      exit(1);
    } else {
      G4SBSRun::GetRun()->GetData()->SetGenName(newValue.data());
    }

    //After any change of kinematics, set this flag to false so that rejection sampling will be re-initialized:
    fevgen->SetInitialized(false);
  }

  if( cmd == PYTHIAfileCmd ){
    fevgen->LoadPythiaChain( newValue );
  }

  if( cmd == expCmd ){
    bool validcmd = false;
    if( newValue.compareTo("gep") == 0 ){
      fExpType = kGEp;
      validcmd = true;
    }
    if( newValue.compareTo("gmn") == 0 ){
      fExpType = kNeutronExp;
      validcmd = true;
    }
    if( newValue.compareTo("gen") == 0 ){
      fExpType = kNeutronExp;
      validcmd = true;
    }
    if( newValue.compareTo("a1n") == 0 ){
      fExpType = kA1n; //"A1n" experiment type for new proposal with both SBS and BigBite in electron mode to detect DIS electrons at high-x: requires some geometry modifications on SBS side, including RICH w/CO2 instead of C4F10 and no aerogel, AND with a non-zero pitch angle for the SBS tracker. Later: HCAL replaced by CLAS LAC?
      validcmd = true;
    }
    //AJP: Add SIDIS as a valid experiment type:
    if( newValue.compareTo("sidis") == 0 ){
      fExpType = kSIDISExp;
      validcmd = true;
    }
    if( newValue.compareTo("C16") == 0 ){
      fExpType = kC16;
      validcmd = true;
    }
    if( newValue.compareTo("tdis") == 0 ){
      fExpType = kTDIS;
      validcmd = true;
    }
    if( newValue.compareTo("ndvcs") == 0 ){
      fExpType = kNDVCS;
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
      fevgen->SetHadronType( kPiPlus );
      validcmd = true;
    }
    if( newValue.compareTo("pi-") == 0 ){
      fevgen->SetHadronType( kPiMinus );
      validcmd = true;
    }
    if( newValue.compareTo("pi0") == 0 ){
      fevgen->SetHadronType( kPi0 );
      validcmd = true;
    }
    if( newValue.compareTo("K+") == 0 ){
      fevgen->SetHadronType( kKPlus );
      validcmd = true;
    }
    if( newValue.compareTo("K-") == 0 ){
      fevgen->SetHadronType( kKMinus );
      validcmd = true; 
    }
    if( newValue.compareTo("p") == 0 ){
      fevgen->SetHadronType( kP );
      validcmd = true;
    } 
    if( newValue.compareTo("pbar") == 0 ){
      fevgen->SetHadronType( kPbar );
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
      fevgen->SetTarget(kLH2);
      fdetcon->SetTarget(kLH2);

      G4double den = (0.071*g/cm3)*Avogadro/(1.008*g/mole);
      fevgen->SetTargDen(den);
      fdetcon->fTargetBuilder->SetTargDen(den);
      validcmd = true;
    }
    if( newValue.compareTo("H2") == 0 ){
      fevgen->SetTarget(kH2);
      fdetcon->SetTarget(kH2);

      G4double den = 10.0*atmosphere/(296.0*kelvin*k_Boltzmann);
      fevgen->SetTargDen(den);
      fdetcon->fTargetBuilder->SetTargDen(den);
      validcmd = true;
    }
    if( newValue.compareTo("D2") == 0 ){
      fevgen->SetTarget(kD2);
      fdetcon->SetTarget(kD2);

      G4double den = 1.0*atmosphere/(77.0*kelvin*k_Boltzmann);
      fevgen->SetTargDen(den);
      fdetcon->fTargetBuilder->SetTargDen(den);
      validcmd = true;
    }
    if( newValue.compareTo("LD2") == 0 ){
      fevgen->SetTarget(kLD2);
      fdetcon->SetTarget(kLD2);

      G4double den = (162.4*kg/m3)*Avogadro/(2.014*g/mole);
      fevgen->SetTargDen(den);
      fdetcon->fTargetBuilder->SetTargDen(den);
      validcmd = true;
    }
    if( newValue.compareTo("3He") == 0 ){
      fevgen->SetTarget(k3He);
      fdetcon->SetTarget(k3He);

      G4double den = 10.0*atmosphere/(296.0*kelvin*k_Boltzmann);
      fevgen->SetTargDen(den);
      fdetcon->fTargetBuilder->SetTargDen(den);
      validcmd = true;

    }
    if( newValue.compareTo("Neutron") == 0 ){
      fevgen->SetTarget(kNeutTarg);
      fdetcon->SetTarget(kNeutTarg);

      G4double den = 10.0*atmosphere/(296.0*kelvin*k_Boltzmann);
      fevgen->SetTargDen(den);
      fdetcon->fTargetBuilder->SetTargDen(den);
      validcmd = true;
    }

    if( !validcmd ){
      fprintf(stderr, "%s: %s line %d - Error: target type %s not valid\n", __PRETTY_FUNCTION__, __FILE__, __LINE__, newValue.data());
      exit(1);
    }

    fevgen->SetInitialized(false);
    
  }

  if( cmd == bigfieldCmd ){
    G4int n = bigfieldCmd->GetNewIntValue(newValue);
    fdetcon->Set48D48Field(n);
  }

  if( cmd == bbfieldCmd ){
    G4int n = bbfieldCmd->GetNewIntValue(newValue);
    fdetcon->SetBigBiteField(n);
  }

  if( cmd == tosfieldCmd ){
    fdetcon->AddToscaField(newValue.data());
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
    G4double den = pre/(296.0*kelvin*k_Boltzmann);
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

  if( cmd == beamECmd ){
    G4double v = beamECmd->GetNewDoubleValue(newValue);
    fevgen->SetBeamE(v);
    fIO->SetBeamE(v);

    G4SBSRun::GetRun()->GetData()->SetBeamE(v/GeV);
    //after any command affecting the kinematics or cross section of the built-in event generators, re-initialize rejection sampling:
    fevgen->SetInitialized(false);
  }

  if( cmd == bbangCmd ){
    G4double v = bbangCmd->GetNewDoubleValue(newValue);
    printf("Setting BB ang to %f deg\n", v/deg);
    fdetcon->SetBBAng(v);
    fIO->SetBigBiteTheta(v);
  }

  if( cmd == bbdistCmd ){
    G4double v = bbdistCmd->GetNewDoubleValue(newValue);
    fdetcon->SetBBDist(v);
    fIO->SetBigBiteDist(v);
  }

  if( cmd == hcalangCmd ){
    G4double v = hcalangCmd->GetNewDoubleValue(newValue);
    fdetcon->Set48D48Ang(v);
    fIO->SetSBSTheta(v);
  }

  if( cmd == sbstrkrpitchCmd ){
    G4double v = sbstrkrpitchCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetTrackerPitch(v);
    G4SBSRun::GetRun()->GetData()->SetSBSTrackerPitch( v );
  }
  
  if( cmd == dvcsecalmatCmd ){
    fdetcon->fEArmBuilder->SetDVCSECalMaterial(newValue);
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
    fevgen->SetHCALDist(v);
    fIO->SetHcalVOffset(v);
  }


  if( cmd == hmagdistCmd ){
    G4double v = hmagdistCmd->GetNewDoubleValue(newValue);
    fdetcon->Set48D48Dist(v);
    fIO->SetSBSDist( v );
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

  if( cmd == gemresCmd ){
    G4double v = gemresCmd->GetNewDoubleValue(newValue);
    fevact->SetGEMRes(v);
  }

  if( cmd == RICHdistCmd ){
    G4double v = RICHdistCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetRICHdist(v);
    fIO->SetRICHDist( v );
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

  if( cmd == buildSBSsieveCmd ){
    fdetcon->fHArmBuilder->SetSBSSieve(newValue);
  }

  if( cmd == buildBBsieveCmd ){
    fdetcon->fEArmBuilder->SetBBSieve(newValue);
  }

  if( cmd == TreeFlagCmd ){
    G4int flag = TreeFlagCmd->GetNewIntValue(newValue);
    fevact->SetTreeFlag( flag );
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
  
  if( cmd == UseCerenkovCmd ){
    G4bool b = UseCerenkovCmd->GetNewBoolValue(newValue);
    fphyslist->ToggleCerenkov(b);
  }

  if( cmd == UseScintCmd ){
    G4bool b = UseScintCmd->GetNewBoolValue(newValue);
    fphyslist->ToggleScintillation(b);
  }

  if( cmd == GunPolarizationCommand ){
    G4ThreeVector pol = GunPolarizationCommand->GetNew3VectorValue(newValue);
    fprigen->SetGunPolarization( pol.unit() );
  }
}
