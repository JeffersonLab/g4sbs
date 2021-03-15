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
#include "G4SBSECal.hh"

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
  fExpType = kGMN; //default to GMN

  runCmd = new G4UIcmdWithAnInteger("/g4sbs/run",this);
  runCmd->SetGuidance("Run simulation with x events");
  runCmd->SetParameterName("nevt", false);

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
  GENTargetHelmholtzCmd->SetGuidance("1 => Q2 = 1.46 (GeV/c)^2, 2 => Q2 = 3.68 (GeV/c)^2, 3 => Q2 = 6.77 (GeV/c)^2, 4 => Q2 = 10.18 (GeV/c)^2, "); 
  GENTargetHelmholtzCmd->SetParameterName("targgenhhconfig",false); // user must provide an integer value, non-argument not allowed 
  GENTargetHelmholtzCmd->SetDefaultValue(kSBS_GEN_146);             // probably not utilized since we require an input value 

  kineCmd = new G4UIcmdWithAString("/g4sbs/kine",this);
  kineCmd->SetGuidance("Kinematics from elastic, inelastic, flat, dis, beam, sidis, wiser, gun, pythia6, wapp, tdiskin, AcquMC");
  kineCmd->SetParameterName("kinetype", false);

  // TDIS
  AcquMCfileCmd = new G4UIcmdWithAString("/g4sbs/acqumcfile",this);
  AcquMCfileCmd->SetGuidance("Name of ROOT file containing AcquMC events as a ROOT tree");
  AcquMCfileCmd->SetParameterName("fname",false);

  PYTHIAfileCmd = new G4UIcmdWithAString("/g4sbs/pythia6file",this);
  PYTHIAfileCmd->SetGuidance("Name of ROOT file containing PYTHIA6 events as a ROOT tree");
  PYTHIAfileCmd->SetParameterName("fname",false);

  exclPythiaXSoptCmd = new G4UIcmdWithAnInteger("/g4sbs/exclpythiaXSoption",this);
  exclPythiaXSoptCmd->SetGuidance("Indicate which cross section should be used by pythia");
  exclPythiaXSoptCmd->SetGuidance("1: total unp. XS; 2: BH2; 3:DVCS2; 4:I; 5: XS(+); 6: XS(-);");
  exclPythiaXSoptCmd->SetGuidance("Must be called *before* /g4sbs/pythia6file to have an effect ");
  exclPythiaXSoptCmd->SetParameterName("exclpyXSopt", false);

  expCmd = new G4UIcmdWithAString("/g4sbs/exp",this);
  expCmd->SetGuidance("Experiment type from gep, gmn, gen, a1n, sidis, C16, tdis, ndvcs, genrp");
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
  
  dvcsecalhoffsetCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/dvcsecalhoffset",this);
  dvcsecalhoffsetCmd->SetGuidance("DVCS ECal horizontal offset");
  dvcsecalhoffsetCmd->SetParameterName("dvcsecalhoffset", false);

  dvcsecalvoffsetCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/dvcsecalvoffset",this);
  dvcsecalvoffsetCmd->SetGuidance("DVCS ECal vertical offset");
  dvcsecalvoffsetCmd->SetParameterName("dvcsecalvoffset", false);

  dvcsecalnrowsCmd = new G4UIcmdWithAnInteger("/g4sbs/dvcsecalnrows",this);
  dvcsecalnrowsCmd->SetGuidance("number of rows for (PbF2) DVCS ECal");
  dvcsecalnrowsCmd->SetParameterName("nrows", false);
  
  dvcsecalncolsCmd = new G4UIcmdWithAnInteger("/g4sbs/dvcsecalncols",this);
  dvcsecalncolsCmd->SetGuidance("number of cols for (PbF2) DVCS ECal");
  dvcsecalncolsCmd->SetParameterName("ncols", false);
  
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
  hcalvoffsetCmd->SetParameterName("hcalvoffset", false);

  hcalhoffsetCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hcalhoffset",this);
  hcalhoffsetCmd->SetGuidance("HCAL horizontal offset");
  hcalhoffsetCmd->SetParameterName("hcalhoffset", false);

  hcalhoffsetCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hcalhoffset",this);
  hcalhoffsetCmd->SetGuidance("HCAL horizontal offset relative to SBS center line (+ = TOWARD beam line)");
  hcalhoffsetCmd->SetParameterName("dist", false);

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
  SBSFieldClampOptionCmd->SetGuidance("SBS field clamp configuration: 0=no clamp, 1=BigBite(default), 2=GEp");
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
  GEPFPP1_CH2thickCmd->SetGuidance("0 < FPP1 CH2 thick < 60 cm");
  GEPFPP1_CH2thickCmd->SetParameterName("CH2thick1",false);
  //GEPFPP1_CH2thickCmd->SetRange("0.0 <= CH2thick1 && CH2thick1 <= 60.0*cm"

  GEPFPP2_CH2thickCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/FPP2CH2thick",this);
  GEPFPP2_CH2thickCmd->SetGuidance("CH2 thickness for first analyzer (GEP only)");
  GEPFPP2_CH2thickCmd->SetGuidance("0 < FPP2 CH2 thick < 60 cm");
  GEPFPP2_CH2thickCmd->SetParameterName("CH2thick2",false);

  GEPFPPoptionCmd = new G4UIcmdWithAnInteger("/g4sbs/gepfppoption",this);
  GEPFPPoptionCmd->SetGuidance("GEP FPP option:");
  GEPFPPoptionCmd->SetGuidance("1 = One analyzer, 8 (FT) + 8 (FPP) GEM trackers");
  GEPFPPoptionCmd->SetGuidance("2 = Two analyzers, 6 (FT) + 5 (FPP1) + 5 (FPP2) GEM trackers (default)");
  GEPFPPoptionCmd->SetParameterName("gepfppoption",true);
  GEPFPPoptionCmd->SetDefaultValue(2);
  
  BLneutronDetsCmd = new G4UIcmdWithABool("/g4sbs/BLneutronDets",this);
  BLneutronDetsCmd->SetGuidance("Setup neutron detectors along the beamline");
  BLneutronDetsCmd->SetParameterName("switch", false);
  
  GEMfrontendCmd = new G4UIcmdWithABool("/g4sbs/buildGEMfrontend",this);
  GEMfrontendCmd->SetGuidance("build GEM front end for GMn or GEp");
  GEMfrontendCmd->SetParameterName("switch", false);

  SetGrinchPMTglassHitsCmd = new G4UIcmdWithABool("/g4sbs/GrinchPMTglassHits",this);
  SetGrinchPMTglassHitsCmd->SetGuidance("build GEM front end for GMn or GEp");
  SetGrinchPMTglassHitsCmd->SetParameterName("switch", false);  
  
  buildSBSsieveCmd = new G4UIcmdWithABool("/g4sbs/buildSBSsieve",this);
  buildSBSsieveCmd->SetGuidance("Use SBS sieve (true or false, false by default)");
  buildSBSsieveCmd->SetParameterName("buildSBSsieve",false);

  buildBBsieveCmd = new G4UIcmdWithABool("/g4sbs/buildBBsieve",this);
  buildBBsieveCmd->SetGuidance("Use BB sieve (true or false, false by default)");
  buildBBsieveCmd->SetParameterName("buildBBsieve",false);
  
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

  // Commands related to solenoid of TPC
  SolUniFieldCmd = new G4UIcmdWithABool("/g4sbs/solunifield", this );
  SolUniFieldCmd->SetGuidance("Switch on uniform field for solenoid" );
  SolUniFieldCmd->SetParameterName("UniformSolField" , false);

  // SolUniFieldMagCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/solunimag", this );
  SolUniFieldMagCmd = new G4UIcmdWithADouble("/g4sbs/solunimag", this );
  SolUniFieldMagCmd->SetGuidance("Set magnitude of uniform field for solenoid in tesla" );
  SolUniFieldMagCmd->SetParameterName("UniformSolFieldMag", false);

  SolTosFieldCmd = new G4UIcmdWithABool("/g4sbs/soltoscafield", this );
  SolTosFieldCmd->SetGuidance("Switch on tosca field for solenoid" );
  SolTosFieldCmd->SetParameterName("ToscaSolField", false);

  SolTosFieldScaleCmd = new G4UIcmdWithADouble("/g4sbs/soltoscascale", this );
  SolTosFieldScaleCmd->SetGuidance("Scale the tosca field for solenoid" );
  SolTosFieldScaleCmd->SetParameterName("ToscaSolFieldScale", false);

  SolTosFieldOffsetCmd= new G4UIcmdWithADoubleAndUnit("/g4sbs/soltoscaoffset", this );
  SolTosFieldOffsetCmd->SetGuidance("Set offset of tosca field for solenoid in mm" );
  SolTosFieldOffsetCmd->SetGuidance("Requires hard coded sol_map_03.dat map" );
  SolTosFieldOffsetCmd->SetParameterName("ToscaSolFieldOffset", false);

  mTPCHeGasRatioCmd = new G4UIcmdWithADouble("/g4sbs/mtpchegasratio", this );
  mTPCHeGasRatioCmd->SetGuidance("Set the fraction of He in the mTPC gas mix. Default 0.9" );
  mTPCHeGasRatioCmd->SetParameterName("mTPCHeGasRatio", false);

  mTPCCH4GasRatioCmd = new G4UIcmdWithADouble("/g4sbs/mtpcch4gasratio", this );
  mTPCCH4GasRatioCmd->SetGuidance("Set the fraction of CH4 in the mTPC gas mix. Default 0.1" );
  mTPCCH4GasRatioCmd->SetParameterName("mTPCCH4GasRatio", false);

  mTPCGasTempCmd = new G4UIcmdWithADouble("/g4sbs/mtpcgastemp", this );
  mTPCGasTempCmd->SetGuidance("Set the temp of the mTPC gas mix in K. Default 296.15K" );
  mTPCGasTempCmd->SetParameterName("mTPCGasTemp", false);

  mTPCGasPressureCmd = new G4UIcmdWithADouble("/g4sbs/mtpcgaspressure", this );
  mTPCGasPressureCmd->SetGuidance("Set the pressure of the mTPC gas mix in atm. Default 0.1atm" );
  mTPCGasPressureCmd->SetParameterName("mTPCGasPressure", false);

  mTPCTgtThickCmd = new G4UIcmdWithADouble("/g4sbs/mtpctargetthick", this );
  mTPCTgtThickCmd->SetGuidance("Set the thickness of the mTPC gas target in mm. Default 0.1atm" );
  mTPCTgtThickCmd->SetParameterName("mTPCTargetThickness", false);

  mTPCkryptoCmd = new G4UIcmdWithABool("/g4sbs/mtpckrypto", this );
  mTPCkryptoCmd->SetGuidance("Set mTPC RO and GEMs as full absorbers" );
  mTPCkryptoCmd->SetParameterName("mTPCkrypto", false);

  mTPCRoomTempCmd = new G4UIcmdWithABool("/g4sbs/setmtpcroomtemp", this );
  mTPCRoomTempCmd->SetGuidance("Set mTPC materials at room temp parameters" );
  mTPCRoomTempCmd->SetParameterName("SetmTPCroomTemp" , false);

  TDIStgtWallThickCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/tdistgtwallthick", this );
  TDIStgtWallThickCmd->SetGuidance("Set TDIS wall target thickness" );
  TDIStgtWallThickCmd->SetParameterName("TDIStgtWallThickness" , false);
 
}

G4SBSMessenger::~G4SBSMessenger(){
}


void G4SBSMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue){
  char cmdstr[255];

  if( cmd == runCmd ){
	
    G4VPhysicalVolume* pWorld;

    G4int nevt = runCmd->GetNewIntValue(newValue);

    //If the generator is PYTHIA, don't try to generate more events than we have available:
    if( fevgen->GetKine() == kPYTHIA6 ){
      if( fevgen->GetPythiaChain()->GetEntries() < nevt ){
	nevt = fevgen->GetPythiaChain()->GetEntries();
      }
      fevgen->InitializePythia6_Tree();
    }

    // TDIS
    if( fevgen->GetKine() == kAcquMC ){
      if( fevgen->GetAcquMCChain()->GetEntries() < nevt ){
	nevt = fevgen->GetAcquMCChain()->GetEntries();
      }
      fevgen->InitializeAcquMC_Tree();
    }
    //    G4double TargMassDensity;
    G4double TargNumberDensity; 
    
    switch(fdetcon->fTargType){ //Initialize fTargDen correctly and consistently with Material definition in fdetcon->ConstructMaterials:
    case kH2:
      TargNumberDensity = fdetcon->GetMaterial("refH2")->GetTotNbOfAtomsPerVolume();
      break;
    case kD2:
      TargNumberDensity = fdetcon->GetMaterial("refD2")->GetTotNbOfAtomsPerVolume();
      break;
    case kNeutTarg:
      TargNumberDensity = fdetcon->GetMaterial("refN2")->GetTotNbOfAtomsPerVolume();
      break;
    case kLH2:
      TargNumberDensity = fdetcon->GetMaterial("LH2")->GetTotNbOfAtomsPerVolume();
      break;
    case kLD2:
      TargNumberDensity = fdetcon->GetMaterial("LD2")->GetTotNbOfAtomsPerVolume();
      break;
    case k3He:
      TargNumberDensity = fdetcon->GetMaterial("pol3He")->GetTotNbOfAtomsPerVolume();
      break;
    case kCfoil:
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
    if( fdetcon->fTargType == kOptics ){
      fevgen->SetNfoils( fdetcon->fTargetBuilder->GetNtargetFoils() );
      fevgen->SetFoilZandThick( fdetcon->fTargetBuilder->GetFoilZpos(), fdetcon->fTargetBuilder->GetFoilThick() );
    }
    // if( fevgen->GetRejectionSamplingFlag() && !fevgen->GetRejectionSamplingInitialized() ){
    //   fevgen->InitializeRejectionSampling();
    // }

    Kine_t kinetype = fevgen->GetKine(); 
    if( kinetype == kDIS || kinetype == kWiser){ //Processes with xsec in units of area/energy/solid angle; i.e., nb/GeV/sr
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
      // TDIS, think it is in solid angle only so in here is ok, but must check
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
    if( newValue.compareTo("cosmics") == 0 ){
      fevgen->SetKine( kCosmics );
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
    // TDIS addition
    if (newValue.compareTo("tdiskin") == 0 ){
      fevgen->SetKine(kTDISKin);
      // fevgen->SetRejectionSamplingFlag(false);
      //fevgen->SetMaxWeight( cm2 );
      validcmd = true;
    }
    // TDIS AcquMC
    if( newValue.compareTo("AcquMC") == 0 ){
      fevgen->SetKine( kAcquMC );
      fIO->SetUseAcquMC( true );
      fevgen->SetRejectionSamplingFlag(false);
      validcmd = true;
    }

    if( newValue.compareTo("wapp") == 0 ){ //wide angle pion photoproduction
      fevgen->SetKine(kPionPhoto);
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
  //TDIS AcquMC
  if( cmd == AcquMCfileCmd ){
    fevgen->LoadAcquMCChain( newValue );
  }
  
  if( cmd == exclPythiaXSoptCmd){
    G4int xsopt = exclPythiaXSoptCmd->GetNewIntValue(newValue);
    fevgen->SetExclPythiaXSOption(xsopt);
    if(xsopt>0)fIO->SetExclPythia6( true );
  }

  if( cmd == expCmd ){
    bool validcmd = false;
    if( newValue.compareTo("gep") == 0 ){
      fExpType = kGEp;
      validcmd = true;
    }
    if( newValue.compareTo("gepeplus") == 0 ){
      fExpType = kGEPpositron;
      validcmd = true;
    }
    if( newValue.compareTo("gmn") == 0 ){
      fExpType = kGMN;
      validcmd = true;
    }
    if( newValue.compareTo("gen") == 0 ){
      fExpType = kGEN;
      validcmd = true;
    }
    if( newValue.compareTo("genrp") == 0 ){
      fExpType = kGEnRP;
      validcmd = true;
    }
    if( newValue.compareTo("a1n") == 0 ){
      fExpType = kA1n; //"A1n" experiment type for new proposal with both SBS and BigBite in electron mode to detect DIS electrons at high-x: requires some geometry modifications on SBS side, including RICH w/CO2 instead of C4F10 and no aerogel, AND with a non-zero pitch angle for the SBS tracker. Also: HCAL + LAC.
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
    if( newValue.compareTo("hcgem") == 0 ){
      fExpType = kGEMHCtest;
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

      // G4double den = 10.0*atmosphere/(296.0*kelvin*k_Boltzmann); //Should this be hard-coded? I think not. On the other hand, this provides a sensible default value, soooo....
      // TDIS temp check
      G4double den = 7.5*atmosphere/(300.0*kelvin*k_Boltzmann);
      fevgen->SetTargDen(den);
      fdetcon->fTargetBuilder->SetTargDen(den);
      validcmd = true;
    }
    if( newValue.compareTo("D2") == 0 ){
      fevgen->SetTarget(kD2);
      fdetcon->SetTarget(kD2);

      // G4double den = 1.0*atmosphere/(77.0*kelvin*k_Boltzmann);
      // TDIS temp check
      G4double den = 7.5*atmosphere/(300.0*kelvin*k_Boltzmann);
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
    if( newValue.compareTo("Cfoil") == 0 ){
      fevgen->SetTarget(kCfoil);
      fdetcon->SetTarget(kCfoil);
      validcmd = true;
    }

    if( newValue.compareTo("optics") == 0 ){
      fevgen->SetTarget(kOptics);
      fdetcon->SetTarget(kOptics);
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

  if( cmd == bigfieldCmd ){
    G4int n = bigfieldCmd->GetNewIntValue(newValue);
    fdetcon->Set48D48Field(n);
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
    fIO->SetSBStrkrPitch( v );
    //G4SBSRun::GetRun()->GetData()->SetSBSTrackerPitch( v );
  }

  if( cmd == sbstrkrdistCmd ){
    G4double d = sbstrkrdistCmd->GetNewDoubleValue(newValue);
    fdetcon->fHArmBuilder->SetTrackerDist(d);
    fIO->SetSBStrkrDist( d );
    //G4SBSRun::GetRun()->GetData()->SetSBSTrackerDist( d );
  }
  
  if( cmd == dvcsecalhoffsetCmd ){
    G4double v = dvcsecalhoffsetCmd->GetNewDoubleValue(newValue);
    fdetcon->fECal->SetDVCSECalHOffset(v);
  }
  
  if( cmd == dvcsecalvoffsetCmd ){
    G4double v = dvcsecalvoffsetCmd->GetNewDoubleValue(newValue);
    fdetcon->fECal->SetDVCSECalVOffset(v);
  }
  
  if( cmd == dvcsecalnrowsCmd ){
    G4double v = dvcsecalnrowsCmd->GetNewIntValue(newValue);
    fdetcon->fECal->SetDVCSECalNRows(v);
  }
  
  if( cmd == dvcsecalncolsCmd ){
    G4double v = dvcsecalncolsCmd->GetNewIntValue(newValue);
    fdetcon->fECal->SetDVCSECalNCols(v);
  }
  
  if( cmd == dvcsecalmatCmd ){
    fdetcon->fECal->SetDVCSECalMaterial(newValue);
    if(newValue == "PbWO4"){
     fdetcon->fECal->SetDVCSECalNRows(31);
     fdetcon->fECal->SetDVCSECalNCols(36);
    }
  }

  if( cmd == GRINCH_gas_Cmd ){
    G4String gasname = newValue;
    
    gasname.toUpper();

    if( gasname.index( "C4F10" ) != gasname.npos ){
      gasname = "C4F10_gas";
    } else if( gasname.index( "C4F8O" ) != gasname.npos ){
      gasname = "C4F8O";
    } else if( gasname.index( "CF4" ) != gasname.npos ){
      gasname = "CF4_gas";
    } else if( gasname.index( "SF6" ) != gasname.npos ){
      gasname = "SF6_gas";
    } else if( gasname.index( "CO2" ) != gasname.npos ){
      gasname = "CO2";
    } else if( gasname.index( "C4F8" ) != gasname.npos ){
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

    G4cout << "gasname = " << gasname << G4endl;

    G4cout << gasname.index( "C4F10" ) << G4endl;
    
    if( gasname.index( "C4F10" ) != gasname.npos ){
      gasname = "C4F10_gas";
    } else if( gasname.index( "C4F8O" ) != gasname.npos ){
      gasname = "C4F8O";
    } else if( gasname.index( "CF4" ) != gasname.npos ){
      gasname = "CF4_gas";
    } else if( gasname.index( "SF6" ) != gasname.npos ){
      gasname = "SF6_gas";
    } else if( gasname.index( "CO2" ) != gasname.npos ){
      gasname = "CO2";
    } else if( gasname.index( "C4F8" ) != gasname.npos ){
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
  
  if( cmd == BLneutronDetsCmd ){
    G4bool v = BLneutronDetsCmd->GetNewBoolValue(newValue);
    fdetcon->fBLneutronDet = v;
  }
  
  if( cmd == GEMfrontendCmd ){
    G4bool v = GEMfrontendCmd->GetNewBoolValue(newValue);
    fdetcon->fEArmBuilder->SetGEMfrontend(v);
  }
  
  if( cmd == SetGrinchPMTglassHitsCmd ){
    G4bool v = SetGrinchPMTglassHitsCmd->GetNewBoolValue(newValue);
    fdetcon->fEArmBuilder->SetGrinchPMTglassHits(v);
  }
  
  if( cmd == buildSBSsieveCmd ){
    G4bool b = buildSBSsieveCmd->GetNewBoolValue(newValue);
    fdetcon->fHArmBuilder->SetSBSSieve(b);
  }

  if( cmd == buildBBsieveCmd ){
    G4bool b = buildBBsieveCmd->GetNewBoolValue(newValue);
    fdetcon->fEArmBuilder->SetBBSieve(b);
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

  if( cmd == SD_EnergyThresholdCmd ){ //store the SDname and dimensioned threshold value in a map<G4String,G4double> assigned to fdetcon?
    std::istringstream is(newValue);

    G4String SDname;
    G4double ethresh;
    G4String unit;

    is >> SDname >> ethresh >> unit;
    
    fdetcon->SDthreshold[SDname] = ethresh*cmd->ValueOf(unit);

    G4cout << "Set Energy threshold for SD name = " << SDname << " to " << fdetcon->SDthreshold[SDname]/MeV << " MeV" << G4endl;
    
  }

  if( cmd == SD_TimeWindowCmd ){ //store the SDname and dimensioned threshold value in a map<G4String,G4double> assigned to fdetcon?
    std::istringstream is(newValue);

    G4String SDname;
    G4double timewindow;
    G4String unit;

    is >> SDname >> timewindow >> unit;
    
    fdetcon->SDgatewidth[SDname] = timewindow*cmd->ValueOf(unit);

    G4cout << "Set time window for SD name = " << SDname << " to " << fdetcon->SDgatewidth[SDname]/ns << " ns" << G4endl;
  }

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
  }

  if( cmd == TargPolMagnitudeCmd ){
    G4double pol = TargPolMagnitudeCmd->GetNewDoubleValue( newValue );
    fdetcon->fTargetBuilder->SetTargPolMag( pol );
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

  // Commands related to solenoid of mTPC, 
  if( cmd == SolUniFieldCmd ){
    G4bool soluniflag = SolUniFieldCmd->GetNewBoolValue(newValue);
    fdetcon->SetTPCSolenoidField();
    fdetcon->fTargetBuilder->SetSolUni(soluniflag);
  }
  if( cmd == SolUniFieldMagCmd ){
    G4double solunimag = SolUniFieldMagCmd->GetNewDoubleValue(newValue);
    fdetcon->fTargetBuilder->SetSolUniMag(solunimag);
  }
  if( cmd == SolTosFieldCmd ){
    G4bool soltosflag = SolTosFieldCmd->GetNewBoolValue(newValue);
    fdetcon->SetTPCSolenoidField();
    fdetcon->fTargetBuilder->SetSolTosca(soltosflag);
  }
  if( cmd == SolTosFieldScaleCmd ){
    G4double soltosscale = SolTosFieldScaleCmd->GetNewDoubleValue(newValue);
    fdetcon->fTargetBuilder->SetSolToscaScale(soltosscale);
  }
  if( cmd == SolTosFieldOffsetCmd ){
    G4double soltosoffset = SolTosFieldOffsetCmd->GetNewDoubleValue(newValue);
    fdetcon->fTargetBuilder->SetSolToscaOffset(soltosoffset);
  }

  //change mtpc sensitive gas settings
  if( cmd == mTPCHeGasRatioCmd ){
    G4double mTPCHeGasRatio = mTPCHeGasRatioCmd->GetNewDoubleValue(newValue);
    fdetcon->SetmTPCHeGasRatio(mTPCHeGasRatio);
  }
  if( cmd == mTPCCH4GasRatioCmd ){
    G4double mTPCCH4GasRatio = mTPCCH4GasRatioCmd->GetNewDoubleValue(newValue);
    fdetcon->SetmTPCCH4GasRatio(mTPCCH4GasRatio);
  }
  if( cmd == mTPCGasTempCmd ){
    G4double mTPCGasTemp = mTPCGasTempCmd->GetNewDoubleValue(newValue);
    fdetcon->SetmTPCGasTemp(mTPCGasTemp);
  }
  if( cmd == mTPCGasPressureCmd ){
    G4double mTPCGasPressure = mTPCGasPressureCmd->GetNewDoubleValue(newValue);
    fdetcon->SetmTPCGasPressure(mTPCGasPressure);
  }
  // mtpc implementation target wall thickness
  if( cmd == mTPCTgtThickCmd ){
    G4double mTPCTgtThick = mTPCTgtThickCmd->GetNewDoubleValue(newValue);
    fdetcon->fTargetBuilder->SetmTPCTgtWallThick(mTPCTgtThick);
  }

  if( cmd == mTPCkryptoCmd ){
    G4bool setmtpckrypto = mTPCkryptoCmd->GetNewBoolValue(newValue);
    fdetcon->fTargetBuilder->SetmTPCkrypto(setmtpckrypto);
  }

  if( cmd == mTPCRoomTempCmd ){
    G4bool setroomtemp = mTPCRoomTempCmd->GetNewBoolValue(newValue);
    //fdetcon->fTargetBuilder->SetmTPCmatAtRoomTemp(setroomtemp);
  }
  
  if( cmd == TDIStgtWallThickCmd ){
    G4double tdistgtwallthick = TDIStgtWallThickCmd->GetNewDoubleValue(newValue);
    fdetcon->fTargetBuilder->SetTDIStgtWallThick(tdistgtwallthick);
  }
  
  

}
