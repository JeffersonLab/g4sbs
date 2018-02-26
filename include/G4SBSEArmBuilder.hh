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
  void SetBBSieve(bool a){fBuildBBSieve = a;};
  
  void SetBBCalDist(double a){ fBBCaldist= a; }
  void SetGEMConfig(int gc ){ fGEMOption = gc; }
  //void SetCDetconfig(int cdetc){ fCDetOption = cdetc; }
  void SetShieldConfig(int sc ){ fShieldOption = sc; }
  
  void MakeBigBite(G4LogicalVolume *);
  void MakeBigCal(G4LogicalVolume *);
  void MakeC16(G4LogicalVolume *);
  void MakeDVCSECal(G4LogicalVolume *);
  
  void MakeCDET(G4double, G4double, G4LogicalVolume *); //
  void MakeGMnGEMShielding(G4LogicalVolume *);

  void SetDVCSECalMaterial(G4String str){ fDVCSECalMaterial = str; }
  void SetDVCSECalHOffset(double a){ fDVCSECALhorizontal_offset = a; }
  
  void MakeBBSieveSlit(G4LogicalVolume *);

  void SetGRINCHgas( G4String str ){ fGRINCHgas = str; }

  double fBBang;
  double fBBdist;
  double fBBCaldist;
  
  G4SBSBigBiteField *fbbfield;

  double fRICHdist; //distance from target of RICH detector

  double fCerDepth;
  double fCerDist;
  double fGEMDist;
  
  // EFuchey: 2017/03/02: flag for BBECal shielding option: 
  // 0: nothing; 1: default 1/4 in SS+0.5mm; 2: 10cm Al + 3cm SS on the side; 3: 10cm Al + 3cm SS on the side; 
  int fShieldOption;
  
  int  fGEMOption;
  //int  fCDetOption;

  int fnzsegments_leadglass_ECAL;
  int fnzsegments_leadglass_C16;
  
  bool fUseLocalField;

  G4String fDVCSECalMaterial;
  double fDVCSECALhorizontal_offset;  // Horizontal offset (from center) of DVCS ECal
  
  bool fBuildBBSieve;

  G4String fGRINCHgas;

private:
};

#endif//__G4SBSEArmBuilder_hh
