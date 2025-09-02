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
  void SetVOff(double a){ fVOff = a; }
  void SetHOff(double a){ fHOff = a; }

  void SetDVCSECalMaterial(G4String str){ fDVCSECalMaterial = str; }
  void SetDVCSECalNRows(G4int nrows){ fDVCSNrows = nrows; }
  void SetDVCSECalNCols(G4int ncols){ fDVCSNcols = ncols; }
  void SetDVCSECalHOffset(double a){ fDVCSECALhorizontal_offset = a; }
  void SetDVCSECalVOffset(double a){ fDVCSECALvertical_offset = a; }
  
  G4LogicalVolume* MakeSuperModule(G4double, G4double, G4double);
  
  void MakeECal_new(G4LogicalVolume *);
  void MakeBigCal(G4LogicalVolume *);
  void MakeC16(G4LogicalVolume *);
  void MakeDVCSECal(G4LogicalVolume *);
  
  int fnzsegments_leadglass_ECAL;
  int fnzsegments_leadglass_C16;
  
  double fAng;
  double fDist;

  G4String fDVCSECalMaterial;
  G4int fDVCSNrows;
  G4int fDVCSNcols;
  G4double fDVCSECALhorizontal_offset;  // Horizontal offset (from center) of DVCS ECal
  G4double fDVCSECALvertical_offset;  // Horizontal offset (from center) of DVCS ECal

  double fVOff;
  double fHOff;
  
private:
};

#endif//__G4SBSECal_hh
