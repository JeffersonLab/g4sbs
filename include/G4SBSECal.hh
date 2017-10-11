#ifndef __G4SBSECal_hh
#define __G4SBSECal_hh

#include "G4SBSComponent.hh"

class G4LogicalVolume;

class G4SBSECal: public G4SBSComponent {
public:
  G4SBSECal(G4SBSDetectorConstruction *);
  ~G4SBSECal();

  void BuildComponent(G4LogicalVolume *);
  
  void SetAng(double a){ fAng = a; }
  void SetDist(double a){ fDist= a; }

  G4LogicalVolume* MakeSuperModule(G4double, G4double, G4double);
  
  void MakeECal_new(G4LogicalVolume *);
  void MakeBigCal(G4LogicalVolume *);
  void MakeC16(G4LogicalVolume *);
  
  int fnzsegments_leadglass_ECAL;
  int fnzsegments_leadglass_C16;

  double fAng;
  double fDist;
  
private:
};

#endif//__G4SBSECal_hh
