#ifndef __G4SBSTargetBuilder_hh
#define __G4SBSTargetBuilder_hh

#include "G4SBSComponent.hh"
#include "sbstypes.hh"
#include "G4ThreeVector.hh"

class G4DetectorConstruction;

class G4SBSTargetBuilder: public G4SBSComponent {
public:
  G4SBSTargetBuilder(G4SBSDetectorConstruction *);
  ~G4SBSTargetBuilder();

  void BuildComponent(G4LogicalVolume *);

  G4LogicalVolume* BuildScatteringChamber(G4LogicalVolume *);
  void BuildCryoTarget(G4LogicalVolume *);
  void BuildGasTarget(G4LogicalVolume *);

  void SetTarget(Targ_t t){fTargType = t;}
  void SetTargLen(double len){ fTargLen = len;}
  void SetTargDen(double den){ fTargDen = den;}
  void SetSchamFlag(int flag){ fSchamFlag = flag; }

  int GetSchamFlag() const { return fSchamFlag; }
  double GetTargLen() const { return fTargLen; }

private:
  double fTargLen;
  double fTargDen;
  G4ThreeVector fTargPos;
  G4ThreeVector fTargDir;
  int fSchamFlag;

  Targ_t fTargType;
};

#endif//__G4SBSTargetBuilder_hh
