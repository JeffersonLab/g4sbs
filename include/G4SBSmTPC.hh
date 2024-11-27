#ifndef __G4SBSmTPC_hh
#define __G4SBSmTPC_hh

#include "sbstypes.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4SBSGEMSD.hh"
#include "G4SBSmTPCSD.hh"
#include "G4SBSCalSD.hh"
#include "G4SBSComponent.hh"
#include "G4SBSTPCTOSCAField2D.hh"
#include "G4Material.hh"

class G4SBSDetectorConstruction;

class G4SBSmTPC: public G4SBSComponent {
public:
  G4SBSmTPC(G4SBSDetectorConstruction *);
  ~G4SBSmTPC();
  
  void BuildComponent(G4LogicalVolume *);
  
  // EPAF 2024/11/18: depreacated
  // void BuildTDISTarget(G4LogicalVolume *);
  // void BuildGasTarget(G4LogicalVolume *);
  //void BuildTPC(G4LogicalVolume *, G4double);
  // Montgomery July 2018, mtpc implement
  void BuildmTPCWalls(G4LogicalVolume *, G4double, G4double, G4double, G4double);
  void BuildmTPCReadouts(G4LogicalVolume *, G4double, G4double, G4double,  G4double);//, G4SBSmTPCSD *);
  void BuildmTPCGEMs(G4LogicalVolume *, G4double, G4double, G4double, G4double);
  void BuildmTPCGasCells(G4LogicalVolume *, G4double, G4double, G4double, G4double);//, G4SBSmTPCSD *);//, G4SBSmTPCSD *);

  // EPAF 2024/11/18: Methods for MLPC design
  void BuildMLPCCathode(G4LogicalVolume *mother, G4double z);
  void BuildMLPCG10Support(G4LogicalVolume *mother, G4double z);
  void BuildMLPCWires(G4LogicalVolume *mother, G4double z, G4double angle, G4int WirePlaneNum);

  // EPAF 2024/11/18: depreacated
  // // EFuchey: 2017/02/10:  This function is now meant to build the cryotarget and target cell only.
  // // This function takes as input the mother logical volume, a rotation matrix, and a 3-vector offset.
  // void BuildStandardCryoTarget(G4LogicalVolume *, G4RotationMatrix *, G4ThreeVector);
  
  // // EFuchey: 2017/02/10: Added those functions to build scattering chamber separately from target,
  // // and avoid, if possible, duplicates of the code actually building the target.
  // void BuildStandardScatCham(G4LogicalVolume *);
  // void BuildGEpScatCham(G4LogicalVolume *);
  // void BuildC16ScatCham(G4LogicalVolume *);
  
  // //  void SetTarget(Targ_t t){fTargType = t;} //defined in sbstype.hh
  // void SetTargLen(G4double len){ fTargLen = len;}
  // void SetTargDen(G4double den){ fTargDen = den;} //Currently, fTargDen has NO effect!
  // void SetTargDiameter(G4double D){ fTargDiameter = D; }
  // void SetSchamFlag(int flag){ fSchamFlag = flag; }
  
  // int GetSchamFlag() const { return fSchamFlag; }
  // G4double GetTargLen() const { return fTargLen; }
  // G4double GetTargDiameter() const { return fTargDiameter; }
  // //G4LogicalVolume *BuildSnoutWindows(G4Box *, G4double, G4double, G4double, G4double, G4double);
  
  // G4bool GetFlux() const { return fFlux; }
  // void SetFlux(G4bool b){fFlux = b;}
  
  // G4bool fUseLocalTPCSolenoid;
  // void SetSolUni(G4bool soluniflag){fSolUni = soluniflag;}
  // void SetSolUniMag(G4double solunimag){fSolUniMag = solunimag;}
  // void SetSolTosca(G4bool soltosflag){fSolTosca = soltosflag;}
  // void SetSolToscaScale(G4double soltosscale){fSolToscaScale = soltosscale;}
  // void SetSolToscaOffset(G4double soltosoffset){ fSolToscaOffset= soltosoffset;}
  
  //mTPC functions
  //void SetmTPCTgtWallThick(double thickness){ftdis_tgt_wallthick = thickness;};
  void SetTPCdesign(int design){fTPCdesign = design;};
  void SetmTPCkrypto(bool iskrypto){fmTPCkrypto = iskrypto;};
  
  
  //void SetmTPCmatAtRoomTemp(bool set){};
  //void SetTDIStgtWallThick(double th){fTDIStgtWallThick = th;};
  
private:
  // G4double fTargLen;
  // G4double fTargDen;
  // G4double fTargDiameter; //diameter of cryotarget
  // G4ThreeVector fTargPos; //Note: fTargPos currently has no effect!
  // G4ThreeVector fTargDir; //Note: fTargDir currently has no effect!
  // int fSchamFlag;

  G4int fTPCdesign;
  
  G4double fZpos;
  G4bool fFlux;

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
  G4int   fmTPC_Ncells;
  // readout discs
  G4double fmTPC_readout_thick;
  // GEMs
  G4int fmTPC_Ngems;
  G4double fmTPC_gem_surf1thick;
  G4double fmTPC_gem_dielecthick;
  G4double fmTPC_gem_surf2thick;
  G4double fmTPC_gap_readoutGEM;
  G4double fmTPC_gap_GEMGEM;
  // HV
  G4double fmTPC_HV_thick;

  // Parameters for MLPC design
  G4int fMLPC_NplaneWires;
  G4int fMLPC_NWiresPerHalfPlane;
  G4double fMLPC_InterWireGap;
  G4double fMLPC_InnerWireDistance;
  G4double fMLPC_SupportThickness;
  G4double fMLPC_RinSupport;
  G4double fMLPC_RoutSupport;
  G4double fMLPC_CathodeRin;
  G4double fMLPC_CathodeMylarThickness;
  G4double fMLPC_CathodeAlThickness;
  G4double fMLPC_WirePlaneThickness;

  //GasSettings
  G4bool fmTPCkrypto;
  G4bool fChkOvLaps;

  //User step limit of particle in the drift gas (CA)
  G4double fmTPCstep;


};

//extern G4SBSmTPC *mTPC;


#endif//__G4SBSmTPC_hh
