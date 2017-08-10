#ifndef __G4SBSCDet_hh
#define __G4SBSCDet_hh

#include "G4SBSComponent.hh"

class G4LogicalVolume;
class G4SBSBigBiteField;

class G4SBSCDet: public G4SBSComponent {
public:
  G4SBSCDet(G4SBSDetectorConstruction *);
  ~G4SBSCDet();

  void BuildComponent(G4LogicalVolume *);
  
  void SetR0(double a){fR0 = a;}
  void SetZ0(double a){fZ0 = a;}
  void SetPlanesHOffset(double a){fPlanesHoffset = a;}
  
  void MakeCDET(G4LogicalVolume *);
  
  double fR0;
  double fZ0;
  double fPlanesHoffset;

private:
};

#endif//__G4SBSCDet_hh
