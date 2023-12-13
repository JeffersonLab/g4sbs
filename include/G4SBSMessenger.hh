#ifndef G4SBSMessenger_HH
#define G4SBSMessenger_HH

#include "globals.hh"
#include "sbstypes.hh"
#include "G4UImessenger.hh"
#include "G4UIcommand.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

class G4SBSIO;
class G4SBSEventGen;
class G4SBSDetectorConstruction;
class G4SBSEventAction;
class G4SBSPrimaryGeneratorAction;
class G4SBSPhysicsList;
class G4SBSSteppingAction;
class G4SBSTrackingAction;

class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWith3Vector;

class G4SBSMessenger : public G4UImessenger {
public:
  G4SBSMessenger();
  ~G4SBSMessenger();
  
  void SetIO( G4SBSIO *io ){ fIO = io; }
  void SetEvGen( G4SBSEventGen *eg ){ fevgen = eg; }
  void SetPriGen( G4SBSPrimaryGeneratorAction *pg ){ fprigen = pg; }
  void SetDetCon( G4SBSDetectorConstruction *dc ){ fdetcon= dc; }
  void SetEvAct( G4SBSEventAction *ev ){ fevact = ev; }
  void SetPhysList( G4SBSPhysicsList *pl ){ fphyslist = pl; }
  void SetTrackingAction( G4SBSTrackingAction *tr ){ ftrkact = tr; }
  void SetSteppingAction( G4SBSSteppingAction *st ){ fstepact = st; }
  
  void SetNewValue(G4UIcommand* cmd, G4String newValue);

private:
  G4SBSIO *fIO;
  G4SBSEventGen *fevgen;
  G4SBSDetectorConstruction *fdetcon;
  G4SBSEventAction *fevact;
  G4SBSPrimaryGeneratorAction *fprigen;
  G4SBSPhysicsList *fphyslist;
  G4SBSTrackingAction *ftrkact;
  G4SBSSteppingAction *fstepact;
  

  G4SBS::Exp_t fExpType;
  
  G4UIcmdWithAnInteger *printCmd;  
  G4UIcmdWithAnInteger *runCmd;
  G4UIcmdWithAString   *fileCmd;
  G4UIcmdWithAString   *tgtCmd;
  
  G4UIcmdWithAString   *sigfileCmd;
  
  G4UIcmdWithAString   *kineCmd;
  G4UIcmdWithAString   *PYTHIAfileCmd; 
  G4UIcmdWithAString   *SIMCfileCmd; 
  G4UIcmdWithAString   *expCmd;
  
  G4UIcmdWithAString   *GunParticleCmd;
  G4UIcmdWithAString   *HadrCmd;
  G4UIcommand    *RejectionSamplingCmd;

  
  G4UIcmdWithAnInteger *bigfieldCmd;
  G4UIcommand *bbfieldCmd;
  //G4UIcmdWithAString *bbfield_fnameCmd; //Set filename for BB magnetic field map
  //  G4UIcmdWithAString *tosfieldCmd;
  G4UIcommand *tosfieldCmd;
  
  G4UIcmdWithAnInteger *eventStatusEveryCmd;

  G4UIcmdWithABool *geantinoCmd;
  G4UIcmdWithABool *invertCmd;
  G4UIcmdWithABool *totalabsCmd;
  G4UIcmdWithABool *checkOverlapCmd;
  
  G4UIcmdWithAnInteger *gemconfigCmd;
  G4UIcmdWithAnInteger *shieldconfigCmd;
  G4UIcmdWithAnInteger *bbpsconfigCmd;
  G4UIcmdWithAnInteger *CDetconfigCmd;

  G4UIcmdWithABool *flipGEMCmd;
  
  G4UIcmdWithAString *ECALmapfileCmd; //Set name of text file with list of active rows and columns for ECAL

  //Flag to build evacuated scattering chamber for gas target:
  //G4UIcmdWithAnInteger *SchamGasTgtCmd;

  G4UIcmdWithADoubleAndUnit *tgtLenCmd;
  G4UIcmdWithADoubleAndUnit *tgtDenCmd;
  G4UIcmdWithADoubleAndUnit *tgtPresCmd;
  G4UIcmdWithADoubleAndUnit *tgtDiamCmd;
  G4UIcmdWithADoubleAndUnit *beamcurCmd;
  G4UIcmdWithADoubleAndUnit *runtimeCmd;
  G4UIcmdWithADoubleAndUnit *rasterxCmd;
  G4UIcmdWithADoubleAndUnit *rasteryCmd;
  G4UIcmdWithADoubleAndUnit *rasterrCmd;
  G4UIcmdWithADoubleAndUnit *beamspotsizeCmd;
  
  //commands controlling pion photoproduction event generation:
  G4UIcmdWithADouble *PionPhoto_tminCmd;
  G4UIcmdWithADouble *PionPhoto_tmaxCmd;
  G4UIcmdWithABool *PionPhoto_useradCmd;
  G4UIcmdWithADouble *PionPhoto_radthickCmd;
  G4UIcmdWithADoubleAndUnit *PionPhoto_radzCmd;
  
  //Optics targets: presently all assumed to be Carbon
  G4UIcmdWithAnInteger *tgtNfoilCmd;
  G4UIcommand *tgtFoilThickCmd; //Foil thickness
  G4UIcommand *tgtFoilZCmd;  //Foil position along z axis:
  
  G4UIcmdWithADoubleAndUnit *beamECmd;
  
  G4UIcmdWithADoubleAndUnit *bbangCmd;
  G4UIcmdWithADoubleAndUnit *bbdistCmd;
  
  G4UIcmdWithADoubleAndUnit *hcaldistCmd;
  G4UIcmdWithADoubleAndUnit *hcalvoffsetCmd;
  G4UIcmdWithADoubleAndUnit *hcalhoffsetCmd;
  G4UIcmdWithADoubleAndUnit *hcalangoffsetCmd;

  G4UIcmdWithABool *CDetReadyCmd;   //Cerenkov

  G4UIcmdWithADoubleAndUnit *lacdistCmd;
  G4UIcmdWithADoubleAndUnit *lacvoffsetCmd;
  G4UIcmdWithADoubleAndUnit *lachoffsetCmd;
  
  
  G4UIcmdWithADoubleAndUnit *hmagdistCmd;
  G4UIcmdWithADoubleAndUnit *hcalangCmd;
  //Add command to set pitch angle for SBS tracker + RICH in "electron mode" 
  G4UIcmdWithADoubleAndUnit *sbstrkrpitchCmd;
  G4UIcmdWithADoubleAndUnit *sbstrkrdistCmd;
  
  G4UIcmdWithAString   *dvcsecalmatCmd;

  G4UIcmdWithAString   *GRINCH_gas_Cmd;
  G4UIcmdWithAString   *RICH_gas_Cmd; 
  
  //These commands set angle generation limits for the electron:
  G4UIcmdWithADoubleAndUnit *thminCmd;
  G4UIcmdWithADoubleAndUnit *thmaxCmd;
  G4UIcmdWithADoubleAndUnit *phminCmd;
  G4UIcmdWithADoubleAndUnit *phmaxCmd;
  //But for inclusive and semi-inclusive reactions, we also need to define energy generation limits for electron and hadron 
  // AND angle generation limits for the hadron:
  G4UIcmdWithADoubleAndUnit *HthminCmd;
  G4UIcmdWithADoubleAndUnit *HthmaxCmd;
  G4UIcmdWithADoubleAndUnit *HphminCmd;
  G4UIcmdWithADoubleAndUnit *HphmaxCmd;
  
  G4UIcmdWithADoubleAndUnit *EhminCmd;
  G4UIcmdWithADoubleAndUnit *EhmaxCmd;
  G4UIcmdWithADoubleAndUnit *EeminCmd;
  G4UIcmdWithADoubleAndUnit *EemaxCmd;

  G4UIcmdWithADoubleAndUnit *cerDepCmd;
  G4UIcmdWithADoubleAndUnit *cerDisCmd;
  G4UIcmdWithADoubleAndUnit *gemSepCmd;
  G4UIcmdWithADoubleAndUnit *bbCalDistCmd;
  
  G4UIcmdWithADoubleAndUnit *gemresCmd;
  
  // Commands needed to specify RICH positioning:
  G4UIcmdWithADoubleAndUnit *RICHdistCmd; //Set RICH distance
  G4UIcmdWithADoubleAndUnit *RICHhoffsetCmd; //Set RICH horizontal offset
  G4UIcmdWithADoubleAndUnit *RICHvoffsetCmd; //Set RICH vertical offset
  G4UIcmdWithABool          *RICHaeroCmd; //Toggle use of RICH aerogel
  G4UIcmdWithADoubleAndUnit *RICHSnoutExtensionCmd; //Set RICH snout extension
  
  // Commands to set configurable properties of SBS:
  G4UIcmdWithADoubleAndUnit *SBSMagFieldCmd;

  //Set overall scale factors for magnetic fields:
  G4UIcmdWithADouble *EARM_ScaleFieldCmd;
  G4UIcmdWithADouble *HARM_ScaleFieldCmd; 
  
  G4UIcmdWithAnInteger      *SBSFieldClampOptionCmd;
  G4UIcmdWithAnInteger      *SBSBeamlineConfCmd;
  G4UIcmdWithAnInteger      *SBSLeadOptionCmd;
  G4UIcmdWithAnInteger      *GENRPAnalyzerOptionCmd;

  G4UIcmdWithADoubleAndUnit  *GEPFPP1_CH2thickCmd;
  G4UIcmdWithADoubleAndUnit  *GEPFPP2_CH2thickCmd;

  G4UIcmdWithAnInteger       *GEPFPPoptionCmd;

  G4UIcmdWithABool *HadronFilterCmd;
  G4UIcmdWithADoubleAndUnit *HadronFilterThickCmd;
  G4UIcmdWithAString *HadronFilterMaterialCmd;
  
  // D. Flay 7/28/20: command to set GEn 3He target Helmholtz coils and magnetic shield orientations 
  G4UIcmdWithAnInteger       *GENTargetHelmholtzCmd; 

  // D. Flay 9/29/20: commands to set the angular misalignment of the GEn 3He target 
  G4UIcmdWithADoubleAndUnit *GENTargetRXCmd; 
  G4UIcmdWithADoubleAndUnit *GENTargetRYCmd; 
  G4UIcmdWithADoubleAndUnit *GENTargetRZCmd;

  // D. Flay 10/9/20: commands to enable/disable the GEn 3He target collimators 
  G4UIcmdWithABool *GENTargetColCmd;
  G4UIcmdWithABool *GENTargetColACmd;
  G4UIcmdWithABool *GENTargetColBCmd;
  G4UIcmdWithABool *GENTargetColCCmd;
  
  // D. Flay 12/9/20: command to enable the GEn target as an SD 
  G4UIcmdWithABool *GENTargetSDEnableCmd; 
 
  // D. Flay 8/25/20 
  // command to set the beam offset 
  G4UIcmdWithADoubleAndUnit *beamOffsetXcmd; 
  G4UIcmdWithADoubleAndUnit *beamOffsetYcmd;
  // command to enable the beam diffuser 
  G4UIcmdWithABool *beamDumpCmd;  
  G4UIcmdWithABool *beamDiffuserCmd;

  // D. Flay 10/15/20 
  // beam angular alignment 
  G4UIcmdWithADoubleAndUnit *beamAngleXcmd;
  G4UIcmdWithADoubleAndUnit *beamAngleYcmd;
  G4UIcmdWithADoubleAndUnit *beamAngleZcmd; 
  // command to enable the ion chamber  
  G4UIcmdWithABool *ionChamberEnableCmd; 
  G4UIcmdWithADoubleAndUnit *ionChamberXCmd; 
  G4UIcmdWithADoubleAndUnit *ionChamberYCmd; 
  G4UIcmdWithADoubleAndUnit *ionChamberZCmd; 
  G4UIcmdWithADoubleAndUnit *ionChamberRXCmd; 
  G4UIcmdWithADoubleAndUnit *ionChamberRYCmd; 
  G4UIcmdWithADoubleAndUnit *ionChamberRZCmd; 

  // D. Flay 11/5/20 
  // beam collimator (for GEn)
  // downstream 
  G4UIcmdWithABool *beamCollimatorEnableDnCmd; 
  G4UIcmdWithADoubleAndUnit *beamCollimatorLDnCmd; 
  G4UIcmdWithADoubleAndUnit *beamCollimatorDminDnCmd; 
  G4UIcmdWithADoubleAndUnit *beamCollimatorDmaxDnCmd; 
  G4UIcmdWithADoubleAndUnit *beamCollimatorXDnCmd; 
  G4UIcmdWithADoubleAndUnit *beamCollimatorYDnCmd; 
  G4UIcmdWithADoubleAndUnit *beamCollimatorZDnCmd; 
  // upstream
  G4UIcmdWithABool *beamCollimatorEnableUpCmd; 
  G4UIcmdWithADoubleAndUnit *beamCollimatorLUpCmd; 
  G4UIcmdWithADoubleAndUnit *beamCollimatorDminUpCmd; 
  G4UIcmdWithADoubleAndUnit *beamCollimatorDmaxUpCmd; 
  G4UIcmdWithADoubleAndUnit *beamCollimatorXUpCmd; 
  G4UIcmdWithADoubleAndUnit *beamCollimatorYUpCmd; 
  G4UIcmdWithADoubleAndUnit *beamCollimatorZUpCmd; 
 
  G4UIcmdWithABool *BLneutronDetsCmd;
  G4UIcmdWithABool *GEMfrontendCmd;

  G4UIcmdWithADoubleAndUnit *GEMfrontendDistCmd;
  G4UIcmdWithADoubleAndUnit *GEMfrontendPosAngleCmd;
  G4UIcmdWithADoubleAndUnit *GEMfrontendRotAngleCmd;
    
  G4UIcmdWithABool *SetGrinchPMTglassHitsCmd;

  G4UIcmdWithABool *buildSBSsieveCmd; //Build the SBS Sieve slit
  //SSeeds 10.4.20
  //Command to set sieve plate in BigBite 0: nothing; 1: old plate; 2: new plate
  G4UIcmdWithAnInteger *buildBBsieveCmd;
  
  G4UIcmdWithAnInteger      *TreeFlagCmd; //Set criteria for filling output root tree

  G4UIcmdWithABool *SBS_FT_absorberCmd; //Command to turn on absorber material in front of SBS FT.
  G4UIcmdWithAString *SBS_FT_absorberMaterialCmd; //Command to set material of SBS FT absorber material (default is aluminum)

  G4UIcmdWithADoubleAndUnit *SBS_FT_absorberThickCmd;
  
  // G4UIcmdWithABool *Earm_CAL_part_cmd;
  // G4UIcmdWithABool *Harm_CAL_part_cmd;

  G4UIcommand *KeepPartCALcmd; //Command to keep extra particle trajectory information in the ROOT tree by sensitive detector name
  G4UIcommand *KeepHistorycmd; //Command to store particle history information in the ROOT tree by sensitive detector name
  G4UIcommand *LimitStepCALcmd; //Command to turn on step limiter physics for sensitive volumes defined as calorimeters, by detector name.

  G4UIcommand *SD_EnergyThresholdCmd;   //Set hit energy threshold for "Calorimeter" type sensitive detectors
  G4UIcommand *SD_TimeWindowCmd;
  G4UIcommand *SD_NTimeBinsCmd;

  G4UIcommand *KeepPulseShapeCmd; //Flag to turn on recording of Pulse Shape info  
  G4UIcommand *KeepSDtrackcmd; //Flag to turn on recording of "sensitive detector" track info
  
  //Commands to activate/de-activate parts of the optical physics list (which are CPU intensive!!!)
  G4UIcmdWithABool *UseCerenkovCmd;   //Cerenkov
  G4UIcmdWithABool *UseScintCmd;      //Scintillation
  G4UIcmdWithAString *DisableOpticalPhotonProductionByMaterialCmd;
  // G4UIcmdWithABool *UseOpRayleighCmd; //Rayleigh for optical photons
  // G4UIcmdWithABool *UseOpAbsorbCmd;   //optical absorption
  // G4UIcmdWithABool *UseOpBdryCmd;     //optical boundary process (reflection/refraction/absorption)
  // G4UIcmdWithABool *UseOpWLSCmd;      //Wavelength shifting of optical photons
  // G4UIcmdWithABool *UseOpMieHGCmd;    //Mie scattering;
  // G4UIcmdWithABool *DisableOpticalPhysicsCmd; //disable CPU-intensive optical photon physics

  G4UIcmdWithABool *FluxCmd; //Make sphere around target and use to compute flux of particles

  //Command to define (fixed) target spin orientation
  G4UIcmdWith3Vector *TargPolDirectionCmd;
  G4UIcmdWithADouble *TargPolMagnitudeCmd;
  G4UIcmdWithADouble *BeamPolMagnitudeCmd; 
  G4UIcmdWith3Vector *BeamPolDirectionCmd;

  //For the SIDIS generator, we probably want the ability to randomize the target polarization direction: 
  G4UIcmdWithABool        *RandomizeTargetSpinCmd; //randomize the target spin in the event generator
  G4UIcmdWithAnInteger    *NumSpinStatesTargCmd;   // Number of target spin states: positive integer = discrete number of states, 0 = randomize in all three dimensions, -1 = randomize in the plane perpendicular to the beam direction
  G4UIcommand             *TargThetaSpinCmd;       // Command to read in vector of target theta spin angles
  G4UIcommand             *TargPhiSpinCmd;         // Command to read in vector of target phi spin angles
  
  // Command to set particle polarization for spin transport calculations:
  // ONLY relevant for particle gun generator!
  G4UIcmdWith3Vector *GunPolarizationCommand;
  G4UIcmdWithAnInteger *SegmentC16Cmd;
  G4UIcmdWithADoubleAndUnit *SegmentThickC16Cmd;
  G4UIcmdWithADouble *DoseRateCmd;
  G4UIcmdWith3VectorAndUnit *CosmicsPointerCommand;
  G4UIcmdWithADoubleAndUnit *CosmicsPointerRadiusCommand;
  G4UIcmdWithADoubleAndUnit *CosmicsMaxAngleCommand;

  G4UIcmdWithABool *WriteFieldMapCmd;

  G4UIcmdWithABool *UseGEMshieldCmd;
  G4UIcmdWithADoubleAndUnit *GEMshieldThickCmd;
  G4UIcmdWithADoubleAndUnit *GEMshieldAirGapThickCmd;
  
  G4UIcmdWithABool *EnableBigBitePlateCmd;
  G4UIcmdWithADoubleAndUnit *SetBigBitePlateThicknessCmd;
  G4UIcmdWithAString *SetBigBitePlateMaterialCmd;
};

#endif//G4SBSMessenger_HH

