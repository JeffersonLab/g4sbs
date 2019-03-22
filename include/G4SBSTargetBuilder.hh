#ifndef __G4SBSTargetBuilder_hh
#define __G4SBSTargetBuilder_hh

#include "G4SBSComponent.hh"
#include "sbstypes.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4SBSGEMSD.hh"
#include "G4SBSmTPCSD.hh"


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
  void BuildTPC(G4LogicalVolume *, G4double);
  // Montgomery July 2018, mtpc implement
  void BuildmTPCWalls(G4LogicalVolume *, G4double, G4double, G4double, G4double);
  void BuildmTPCReadouts(G4LogicalVolume *, G4double, G4double, G4double,  G4double);
  void BuildmTPCGEMs(G4LogicalVolume *, G4double, G4double, G4double, G4double);
  void BuildmTPCGasCells(G4LogicalVolume *, G4double, G4double, G4double, G4double, G4SBSmTPCSD *);

  // EFuchey: 2017/02/10:  This function is now meant to build the cryotarget and target cell only.
  // This function takes as input the mother logical volume, a rotation matrix, and a 3-vector offset.
  void BuildStandardCryoTarget(G4LogicalVolume *, G4RotationMatrix *, G4ThreeVector);

  // EFuchey: 2017/02/10: Added those functions to build scattering chamber separately from target,
  // and avoid, if possible, duplicates of the code actually building the target.
  void BuildStandardScatCham(G4LogicalVolume *);
  void BuildGEpScatCham(G4LogicalVolume *);
  void BuildC16ScatCham(G4LogicalVolume *);
  
  void SetTarget(Targ_t t){fTargType = t;}
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

  //mTPC functions
  void SetmTPCTgtWallThick(double thickness){ftdis_tgt_wallthick = thickness;};
  void SetmTPCkrypto(bool iskrypto){fmTPCkrypto = iskrypto;};


  //void SetmTPCmatAtRoomTemp(bool set){};
  void SetTDIStgtWallThick(double th){fTDIStgtWallThick = th;};
  
private:
  G4double fTargLen;
  G4double fTargDen;
  G4double fTargDiameter; //diameter of cryotarget
  G4ThreeVector fTargPos; //Note: fTargPos currently has no effect!
  G4ThreeVector fTargDir; //Note: fTargDir currently has no effect!
  int fSchamFlag;

  G4bool fFlux;
  
  Targ_t fTargType;

  // Montgomery 2018, tdis solenoid implement
  G4bool fSolUni;
  G4double fSolUniMag;
  G4bool fSolTosca;
  G4double fSolToscaScale;
  G4double fSolToscaOffset;
  
  G4double fTDIStgtWallThick;

  // Montgomery July 2018, mtpc implement
  G4double ftdis_tgt_diam;
  G4double ftdis_tgt_wallthick;
  G4double ftdis_tgt_len;
  // variables for mtpc construction
  // taken from M. carmignotto gemc mtpc implementation
  // inner electrode at r=5cm
  G4double fmTPC_inelectrode_r;
  G4double fmTPC_inelectrode_kaptonthick;
  G4double fmTPC_inelectrode_authick;
  // outer electrode at r=15cm
  G4double fmTPC_outelectrode_r;
  G4double fmTPC_outelectrode_kaptonthick;
  G4double fmTPC_outelectrode_authick;
  // mtpc chambers
  G4double fmTPC_cell_len;
  G4int fmTPC_Ncells;
  // readout discs
  G4double fmTPC_readout_thick;
  // GEMs
  G4int fmTPC_Ngems;
  G4double fmTPC_gem_surf1thick;
  G4double fmTPC_gem_dielecthick;
  G4double fmTPC_gem_surf2thick;
  G4double fmTPC_gap_readoutGEM;
  G4double fmTPC_gap_GEMGEM;
  //GasSettings
  G4bool fmTPCkrypto;


};

#endif//__G4SBSTargetBuilder_hh
