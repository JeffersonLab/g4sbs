#ifndef __G4SBSTargetBuilder_hh
#define __G4SBSTargetBuilder_hh

#include "G4SBSComponent.hh"
#include "sbstypes.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
// #include "G4SBSGEMSD.hh"
// #include "G4SBSmTPCSD.hh"
// #include "G4SBSCalSD.hh"

#include "G4SBSmTPC.hh" //CA
class G4SBSmTPC;//CA


#include "G4SBSPartParameters.hh"

#include "G4SBSIonChamberSD.hh"
#include "G4SBSTargetSD.hh"

class G4DetectorConstruction;

class G4SBSTargetBuilder: public G4SBSComponent {
public:
  G4SBSTargetBuilder(G4SBSDetectorConstruction *);
  ~G4SBSTargetBuilder();

  void BuildComponent(G4LogicalVolume *);

  void BuildCryoTarget(G4LogicalVolume *);
  // void BuildC16CryoTarget(G4LogicalVolume *);// EFuchey: 2017/02/10: Now defunct \-> 
  // Replaced by : BuildC16ScatCham
  void BuildTDISTarget(G4LogicalVolume *);
  void BuildGasTarget(G4LogicalVolume *);
  // void BuildTPC(G4LogicalVolume *, G4double);
  // Montgomery July 2018, mtpc implement
  // void BuildmTPCWalls(G4LogicalVolume *, G4double, G4double, G4double, G4double);
  // void BuildmTPCReadouts(G4LogicalVolume *, G4double, G4double, G4double,  G4double);//, G4SBSmTPCSD *);
  // void BuildmTPCGEMs(G4LogicalVolume *, G4double, G4double, G4double, G4double);
  // void BuildmTPCGasCells(G4LogicalVolume *, G4double, G4double, G4double, G4double);//, G4SBSmTPCSD *);//, G4SBSmTPCSD *);
  // void BuildmTPCHVPlanes(G4LogicalVolume *, G4double, G4double, G4double,  G4double, G4SBSmTPCSD *);


  // EFuchey: 2017/02/10:  This function is now meant to build the cryotarget and target cell only.
  // This function takes as input the mother logical volume, a rotation matrix, and a 3-vector offset.
  void BuildStandardCryoTarget(G4LogicalVolume *, G4RotationMatrix *, G4ThreeVector);
  void BuildCfoil(G4LogicalVolume *, G4RotationMatrix *, G4ThreeVector);

  void BuildOpticsTarget(G4LogicalVolume *, G4RotationMatrix *, G4ThreeVector);
  
  // EFuchey: 2017/02/10: Added those functions to build scattering chamber separately from target,
  // and avoid, if possible, duplicates of the code actually building the target.
  void BuildStandardScatCham(G4LogicalVolume *);
  void BuildGEpScatCham(G4LogicalVolume *);
  void BuildC16ScatCham(G4LogicalVolume *);
  //Add "Toy" Scattering chamber for new proposal development:
  void BuildToyScatCham(G4LogicalVolume *);

  // GEn Polarized 3He target (D Flay, July 2020) 
  void BuildGEnTarget(G4LogicalVolume *motherLog); 
  void BuildGEnTarget_GlassCell(G4LogicalVolume *motherLog); 
  void BuildGEnTarget_EndWindows(G4LogicalVolume *motherLog); 
  void BuildGEnTarget_EndWindows_solidCu(G4LogicalVolume *motherLog); 
  void BuildGEnTarget_PolarizedHe3(G4LogicalVolume *motherLog); 
  void BuildGEnTarget_HelmholtzCoils(const int config,const std::string type,G4LogicalVolume *motherLog); 
  void BuildGEnTarget_Shield(const int config,G4LogicalVolume *motherLog); 
  void BuildGEnTarget_LadderPlate(G4LogicalVolume *motherLog); 
  void BuildGEnTarget_PickupCoils(G4LogicalVolume *motherLog);
  void BuildGEnTarget_Collimators(G4LogicalVolume *motherLog,G4double z0=0);  
  void BuildGEnTarget_Collimator_A(G4LogicalVolume *motherLog,G4double z0=0);  
  void BuildGEnTarget_Collimator_B(G4LogicalVolume *motherLog,G4double z0=0);  
  void BuildGEnTarget_Collimator_C(G4LogicalVolume *motherLog,G4double z0=0); 
  void BuildGEnTarget_Collimator_Table(G4LogicalVolume *motherLog,G4double z0=0);
  // test items 
  void BuildGEnTarget_IonChamber(G4LogicalVolume *motherLog);
  void BuildGEnTarget_BeamCollimator(G4LogicalVolume *motherLog,int type=0); // default is downstream  

  void CheckZPos(G4LogicalVolume *logicMother,G4double z0); // dummy function to check z positoning  
  
  void SetTarget(G4SBS::Targ_t t){fTargType = t;}
  void SetTargLen(G4double len){ fTargLen = len;}
  void SetTargDen(G4double den){ fTargDen = den;} //Currently, fTargDen has NO effect!
  void SetTargDiameter(G4double D){ fTargDiameter = D; }
  void SetSchamFlag(int flag){ fSchamFlag = flag; }

  int GetSchamFlag() const { return fSchamFlag; }
  G4double GetTargLen() const { return fTargLen; }
  G4double GetTargDiameter() const { return fTargDiameter; }
  //G4LogicalVolume *BuildSnoutWindows(G4Box *, G4double, G4double, G4double, G4double, G4double);

  G4bool GetFlux() const { return fFlux; }
  void SetFlux(G4bool b){fFlux = b;}

  G4bool fUseLocalTPCSolenoid;
  void SetSolUni(G4bool soluniflag){fSolUni = soluniflag;}
  void SetSolUniMag(G4double solunimag){fSolUniMag = solunimag;}
  void SetSolTosca(G4bool soltosflag){fSolTosca = soltosflag;}
  void SetSolToscaScale(G4double soltosscale){fSolToscaScale = soltosscale;}
  void SetSolToscaOffset(G4double soltosoffset){ fSolToscaOffset= soltosoffset;}

  // //mTPC functions
  // void SetmTPCTgtWallThick(double thickness){ftdis_tgt_wallthick = thickness;};
  // void SetmTPCkrypto(bool iskrypto){fmTPCkrypto = iskrypto;};
  // //void SetmTPCmatAtRoomTemp(bool set){};
  
  void SetTDIStgtWallThick(double th){fTDIStgtWallThick = th;};
  
  G4ThreeVector GetTargPolDir() const { return fTargPolDir; }
  G4double GetTargPolMag() const { return fTargPolMag; }
  void SetTargPolDir( G4ThreeVector pdir ){ fTargPolDir = pdir.unit(); }
  void SetTargPolMag( G4double pmag ){ fTargPolMag = pmag; }

  G4int GetNtargetFoils() const { return fNtargetFoils; }
  void SetNtargetFoils(G4int n);

  std::vector<G4double> GetFoilThick() const { return fFoilThick; }
  std::vector<G4double> GetFoilZpos() const { return fFoilZpos; }
  void SetFoilThick(G4int, G4double);
  void SetFoilZpos(G4int, G4double);

  void SetUseRad( G4bool b ){fUseRad = b; }
  void SetRadThick( G4double v ){ fRadThick = v; }
  void SetRadZoffset( G4double v ){ fRadZoffset = v; }

  void BuildRadiator(G4LogicalVolume *, G4RotationMatrix *, G4ThreeVector );
  
private:
  //Multi-foil solid targets (only Carbon available for now):
  G4int fNtargetFoils;
  std::vector<G4double> fFoilThick; //foil thickness
  std::vector<G4double> fFoilZpos; //foil Z position along beamline

  G4double fTargLen;
  G4double fTargDen;
  G4double fTargDiameter; //diameter of cryotarget
  G4ThreeVector fTargPos; //Note: fTargPos currently has no effect!
  G4ThreeVector fTargDir; //Note: fTargDir currently has no effect!

  G4ThreeVector fTargPolDir; //direction of target polarization
  G4double      fTargPolMag; //magnitude of target polarization.
  
  int fSchamFlag;

  G4bool fFlux;

  G4bool fUseRad; //use radiator?
  G4double fRadThick; //Thickness in units of X0;
  G4double fRadZoffset;  //Distance upstream of target

  G4double fGEn_GLASS_TUBE_LENGTH; // length of GEn 3He glass tube (target cell)  
  
  // Montgomery 2018, tdis solenoid implement
  G4bool fSolUni;
  G4double fSolUniMag;
  G4bool fSolTosca;
  G4double fSolToscaScale;
  G4double fSolToscaOffset;
  
  G4double fTDIStgtWallThick;

  // // Montgomery July 2018, mtpc implement
  G4double ftdis_tgt_diam;
  G4double ftdis_tgt_wallthick;
  G4double ftdis_tgt_len;
  
  // // variables for mtpc construction
  // // taken from M. carmignotto gemc mtpc implementation
  // // inner electrode at r=5cm
  // G4double fmTPC_inelectrode_r;
  // G4double fmTPC_inelectrode_kaptonthick;
  // G4double fmTPC_inelectrode_authick;
  // // outer electrode at r=15cm
  // G4double fmTPC_outelectrode_r;
  // G4double fmTPC_outelectrode_kaptonthick;
  // G4double fmTPC_outelectrode_authick;
  // // mtpc chambers
  // G4double fmTPC_cell_len;
  // G4int fmTPC_Ncells;
  // // readout discs
  // G4double fmTPC_readout_thick;
  // // GEMs
  // G4int fmTPC_Ngems;
  // G4double fmTPC_gem_surf1thick;
  // G4double fmTPC_gem_dielecthick;
  // G4double fmTPC_gem_surf2thick;
  // G4double fmTPC_gap_readoutGEM;
  // G4double fmTPC_gap_GEMGEM;
  // // HV
  // G4double fmTPC_HV_thick;
  // //GasSettings
  // G4bool fmTPCkrypto;
  G4bool fChkOvLaps;
  
  //G4SBSmTPC *mTPC;//CA

  G4SBS::Targ_t fTargType;
};

#endif//__G4SBSTargetBuilder_hh
