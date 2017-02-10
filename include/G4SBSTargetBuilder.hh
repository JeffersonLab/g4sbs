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
  void BuildGasTarget(G4LogicalVolume *);
  
  // EFuchey: 2017/02/10:  This function is now meant to build the cryotarget and target cell only.
  // This function takes as input the mother logical volume, a rotation angle, and a Y offset.
  void BuildStandardCryoTarget(G4LogicalVolume *, G4double, G4double);

  // EFuchey: 2017/02/10: Added those functions to build scattering chamber separately from target,
  // and avoid, if possible, duplicates of the code actually building the target.
  void BuildStandardScatCham(G4LogicalVolume *);
  void BuildGEpScatCham(G4LogicalVolume *);
  void BuildC16ScatCham(G4LogicalVolume *);
  
  void SetTarget(Targ_t t){fTargType = t;}
  void SetTargLen(double len){ fTargLen = len;}
  void SetTargDen(double den){ fTargDen = den;} //Currently, fTargDen has NO effect!
  void SetSchamFlag(int flag){ fSchamFlag = flag; }

  int GetSchamFlag() const { return fSchamFlag; }
  double GetTargLen() const { return fTargLen; }
  //G4LogicalVolume *BuildSnoutWindows(G4Box *, G4double, G4double, G4double, G4double, G4double);

  G4bool GetFlux() const { return fFlux; }
  void SetFlux(G4bool b){fFlux = b;}
  
private:
  double fTargLen;
  double fTargDen;
  G4ThreeVector fTargPos;
  G4ThreeVector fTargDir;
  int fSchamFlag;

  G4bool fFlux;
  
  Targ_t fTargType;
};

#endif//__G4SBSTargetBuilder_hh
