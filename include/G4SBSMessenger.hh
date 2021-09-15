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

class G4UIcmdWithAnInteger;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;
class G4UIcmdWithABool;
class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;

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

  void SetNewValue(G4UIcommand* cmd, G4String newValue);

private:
  G4SBSIO *fIO;
  G4SBSEventGen *fevgen;
  G4SBSDetectorConstruction *fdetcon;
  G4SBSEventAction *fevact;
  G4SBSPrimaryGeneratorAction *fprigen;
  G4SBSPhysicsList *fphyslist;
  

  Exp_t fExpType;
  
  G4UIcmdWithAnInteger *runCmd;
  G4UIcmdWithAString   *fileCmd;
  G4UIcmdWithAString   *tgtCmd;
  
  G4UIcmdWithAString   *sigfileCmd;
  
  G4UIcmdWithAString   *kineCmd;
  G4UIcmdWithAString   *PYTHIAfileCmd; 
  G4UIcmdWithAString   *expCmd;
  
  G4UIcmdWithAString   *GunParticleCmd;
  G4UIcmdWithAString   *HadrCmd;
  G4UIcommand    *RejectionSamplingCmd;

  
  G4UIcmdWithAnInteger *bigfieldCmd;
  G4UIcmdWithAnInteger *bbfieldCmd;
  G4UIcmdWithAString *tosfieldCmd;

  G4UIcmdWithAnInteger *eventStatusEveryCmd;

  G4UIcmdWithABool *geantinoCmd;
  G4UIcmdWithABool *invertCmd;
  G4UIcmdWithABool *totalabsCmd;
  G4UIcmdWithABool *checkOverlapCmd;
  
  G4UIcmdWithAnInteger *gemconfigCmd;
  G4UIcmdWithAnInteger *shieldconfigCmd;
  G4UIcmdWithAnInteger *CDetconfigCmd;
  G4UIcmdWithAnInteger *hcalcosmicconfigCmd;

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
  
  G4UIcmdWithADoubleAndUnit *beamECmd;
  
  G4UIcmdWithADoubleAndUnit *bbangCmd;
  G4UIcmdWithADoubleAndUnit *bbdistCmd;
  
  G4UIcmdWithADoubleAndUnit *hcaldistCmd;
  G4UIcmdWithADoubleAndUnit *hcalvoffsetCmd;
  G4UIcmdWithADoubleAndUnit *hcalhoffsetCmd;

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
  
  // Commands to set configurable properties of SBS:
  G4UIcmdWithADoubleAndUnit *SBSMagFieldCmd;

  //Set overall scale factors for magnetic fields:
  G4UIcmdWithADouble *EARM_ScaleFieldCmd;
  G4UIcmdWithADouble *HARM_ScaleFieldCmd; 
  
  G4UIcmdWithAnInteger      *SBSFieldClampOptionCmd;
  G4UIcmdWithAnInteger      *SBSBeamlineConfCmd;
  G4UIcmdWithAnInteger      *SBSLeadOptionCmd;

  G4UIcmdWithADoubleAndUnit  *GEPFPP1_CH2thickCmd;
  G4UIcmdWithADoubleAndUnit  *GEPFPP2_CH2thickCmd;
  
  G4UIcmdWithABool *BLneutronDetsCmd;
  
  G4UIcmdWithABool *buildSBSsieveCmd; //Build the SBS Sieve slit
  G4UIcmdWithABool *buildBBsieveCmd; //Build the BB Sieve slit
  
  G4UIcmdWithAnInteger      *TreeFlagCmd; //Set criteria for filling output root tree

  
  
  // G4UIcmdWithABool *Earm_CAL_part_cmd;
  // G4UIcmdWithABool *Harm_CAL_part_cmd;

  G4UIcommand *KeepPartCALcmd; //Command to keep extra particle trajectory information in the ROOT tree by sensitive detector name
  G4UIcommand *KeepHistorycmd; //Command to store particle history information in the ROOT tree by sensitive detector name
  G4UIcommand *LimitStepCALcmd; //Command to turn on step limiter physics for sensitive volumes defined as calorimeters, by detector name.

  G4UIcommand *SD_EnergyThresholdCmd;   //Set hit energy threshold for "Calorimeter" type sensitive detectors
  G4UIcommand *SD_TimeWindowCmd; 
  
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
  
  // Command to set particle polarization for spin transport calculations:
  // ONLY relevant for particle gun generator!
  G4UIcmdWith3Vector *GunPolarizationCommand;
  G4UIcmdWith3VectorAndUnit *GunVertexCmd;
  G4UIcmdWithAnInteger *SegmentC16Cmd;
  G4UIcmdWithADoubleAndUnit *SegmentThickC16Cmd;
  G4UIcmdWithADouble *DoseRateCmd;
  G4UIcmdWith3VectorAndUnit *CosmicsPointerCommand;
  G4UIcmdWithADoubleAndUnit *CosmicsPointerRadiusCommand;
  G4UIcmdWithADoubleAndUnit *CosmicsMaxAngleCommand;

  // Command to add components that should not be built
  G4UIcmdWithAString *disableComponentCmd;
};

#endif//G4SBSMessenger_HH























