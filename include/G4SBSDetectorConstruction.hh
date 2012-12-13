#ifndef G4SBSDetectorConstruction_h
#define G4SBSDetectorConstruction_h 1

#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "sbstypes.hh"
#include "globals.hh"

#include "G4SBSBigBiteField.hh"

class G4SBSDetectorMessenger;

class G4SBSDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  G4SBSDetectorConstruction();
  ~G4SBSDetectorConstruction();
  
public:
  G4VPhysicalVolume* Construct();
  G4VPhysicalVolume* ConstructAll();
  G4VPhysicalVolume* ConstructAllGEp();

  void SetTarget(Targ_t t){fTargType = t;}
  void SetTargLen(double len){ fTargLen = len;}
  void SetTargDen(double den){ fTargDen = den;}

  void SetBBAng(double a){ fBBang = a; }
  void SetBBDist(double a){ fBBdist= a; fbbfield->SetOffset(a); }

  void SetHCALAng(double a){ f48D48ang = a; }
  void SetHCALDist(double a){ fHCALdist= a; }
  void Set48D48Dist(double a){ f48D48dist= a; }

  void SetCerDepth(double a){ fCerDepth = a; }
  void SetCerDist(double a){fCerDist = a;}

  void SetGEMSep(double a){fGEMDist = a;}

  void SetBBCalDist(double a){ fBBCaldist= a; }

  void SetGEMConfig(int gc ){ fGEMOption = gc; }

  G4SBSBigBiteField *GetBBField(){ return fbbfield; }

private:
  // messeneger
  G4SBSDetectorMessenger* theMessenger;
  G4SBSBigBiteField *fbbfield;
  
  double fBBang;
  double fBBdist;
  double fBBCaldist;

  double f48D48ang;
  double f48D48dist;
  double fHCALdist;

  double fCerDepth;
  double fCerDist;
  double fGEMDist;

  double fTargLen;
  double fTargDen;

  Targ_t fTargType;

  int  fGEMOption;

};


#endif


