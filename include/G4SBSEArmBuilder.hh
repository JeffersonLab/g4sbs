#ifndef __G4SBSEArmBuilder_hh
#define __G4SBSEArmBuilder_hh

#include "G4SBSComponent.hh"
#include "G4ThreeVector.hh"

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
  void SetBBPSOption(int psopt ){ fBBPSOption = psopt; }
  
  void MakeBigBite(G4LogicalVolume *);
  
  void MakeGMnGEMShielding(G4LogicalVolume *);

  void SetDVCSECalMaterial(G4String str){ fDVCSECalMaterial = str; }
  void MakeBBSieveSlit(G4LogicalVolume *, G4ThreeVector);
  
  void SetGRINCHgas( G4String str ){ fGRINCHgas = str; }
  void SetGrinchPMTglassHits(bool b ){ fTurnOnGrinchPMTglassHits = b; }
  
  void SetGEMfrontend(bool b ){ fBuildGEMfrontend = b; }

  void MakeHallCGEM(G4LogicalVolume *);
  
  double fBBang;
  double fBBdist;
  double fBBCaldist;

  //G4SBSBigBiteField *fbbfield; //Why do we need this in both EArmBuilder and DetectorConstruction?
  
  double fRICHdist; //distance from target of RICH detector

  double fCerDepth;
  double fCerDist;
  double fGEMDist;
  
  // EFuchey: 2017/03/02: flag for BBECal shielding option: 
  // 0: nothing; 1: default 1/4 in SS+0.5mm; 2: 10cm Al + 3cm SS on the side; 3: 10cm Al + 3cm SS on the side; 
  int fShieldOption;
  //New geometry for BB PS: 0: default (old blocks, 27 rows); 1: new blocks, 25 rows; 2: new blocks; 26 rows
  int fBBPSOption;
  
  int  fGEMOption;
  //int  fCDetOption;

  int fnzsegments_leadglass_ECAL;
  int fnzsegments_leadglass_C16;
  
  bool fUseLocalField;

  bool fBuildBBSieve;

  G4String fGRINCHgas;
  G4bool fTurnOnGrinchPMTglassHits;
  G4bool fBuildGEMfrontend;
  
  G4String fDVCSECalMaterial;
  
private:
};

#endif//__G4SBSEArmBuilder_hh
