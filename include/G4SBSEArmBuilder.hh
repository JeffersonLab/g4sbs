#ifndef __G4SBSEArmBuilder_hh
#define __G4SBSEArmBuilder_hh

#include "G4SBSComponent.hh"

class G4LogicalVolume;
class G4SBSBigBiteField;

class G4SBSEArmBuilder: public G4SBSComponent {
public:
  G4SBSEArmBuilder(G4SBSDetectorConstruction *);
  ~G4SBSEArmBuilder();

  void BuildComponent(G4LogicalVolume *);

  void SetBBAng(double a){ fBBang = a; }
  void SetBBDist(double a){ fBBdist= a; }

  void SetCerDepth(double a){ fCerDepth = a; }
  void SetCerDist(double a){fCerDist = a;}

  void SetGEMSep(double a){fGEMDist = a;}

  void SetBBCalDist(double a){ fBBCaldist= a; }
  void SetGEMConfig(int gc ){ fGEMOption = gc; }
  //void SetCDetconfig(int cdetc){ fCDetOption = cdetc; }

  void MakeBigBite(G4LogicalVolume *);
  void MakeBigCal(G4LogicalVolume *);
  void MakeC16(G4LogicalVolume *);
  
  void MakeCDET(G4double, G4double, G4LogicalVolume *); //

  double fBBang;
  double fBBdist;
  double fBBCaldist;

  G4SBSBigBiteField *fbbfield;

  double fRICHdist; //distance from target of RICH detector

  double fCerDepth;
  double fCerDist;
  double fGEMDist;

  int  fGEMOption;
  //int  fCDetOption;

  int fnzsegments_leadglass_ECAL;
  int fnzsegments_leadglass_C16;
  
  bool fUseLocalField;

private:
};

#endif//__G4SBSEArmBuilder_hh
