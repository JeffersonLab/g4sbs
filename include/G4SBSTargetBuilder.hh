#ifndef __G4SBSTargetBuilder_hh
#define __G4SBSTargetBuilder_hh

#include "G4SBSComponent.hh"
#include "sbstypes.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"

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
  
  // EFuchey: 2017/02/10:  This function is now meant to build the cryotarget and target cell only.
  // This function takes as input the mother logical volume, a rotation matrix, and a 3-vector offset.
  void BuildStandardCryoTarget(G4LogicalVolume *, G4RotationMatrix *, G4ThreeVector);
  void BuildCfoil(G4LogicalVolume *, G4RotationMatrix *, G4ThreeVector);

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

  G4ThreeVector GetTargPolDir() const { return fTargPolDir; }
  G4double GetTargPolMag() const { return fTargPolMag; }
  void SetTargPolDir( G4ThreeVector pdir ){ fTargPolDir = pdir.unit(); }
  void SetTargPolMag( G4double pmag ){ fTargPolMag = pmag; }
  
private:
  G4double fTargLen;
  G4double fTargDen;
  G4double fTargDiameter; //diameter of cryotarget
  G4ThreeVector fTargPos; //Note: fTargPos currently has no effect!
  G4ThreeVector fTargDir; //Note: fTargDir currently has no effect!

  G4ThreeVector fTargPolDir; //direction of target polarization
  G4double      fTargPolMag; //magnitude of target polarization.
  
  int fSchamFlag;

  G4bool fFlux;
  
  Targ_t fTargType;
};

#endif//__G4SBSTargetBuilder_hh
