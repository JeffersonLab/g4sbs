#ifndef G4SBSDetectorConstruction_h
#define G4SBSDetectorConstruction_h 1

#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "sbstypes.hh"
#include "globals.hh"


class G4SBSGlobalField;
class G4SBSMagneticField;
class G4SBSBigBiteField;
class G4SBS48D48Field;


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

  void ConstructBeamline(G4LogicalVolume*);
  void ConstructTarget(G4LogicalVolume*);
  void Make48D48(G4LogicalVolume*, double);

  void SetTarget(Targ_t t){fTargType = t;}
  void SetTargLen(double len){ fTargLen = len;}
  void SetTargDen(double den){ fTargDen = den;}

  void SetBBAng(double a);
  void SetBBDist(double a);

  void SetHCALAng(double a);
  void SetHCALDist(double a){ fHCALdist= a;   }
  void Set48D48Dist(double a);

  void SetCerDepth(double a){ fCerDepth = a; }
  void SetCerDist(double a){fCerDist = a;}

  void SetGEMSep(double a){fGEMDist = a;}

  void SetBBCalDist(double a){ fBBCaldist= a; }

  void SetGEMConfig(int gc ){ fGEMOption = gc; }

  G4SBSGlobalField *GetGlobalField(){ return fGlobalField; }
  G4SBSBigBiteField *GetBBField(){ return fbbfield; }
  G4SBS48D48Field *Get48D48Field(){ return f48d48field; }

  void Set48D48Field(int n);

  void SetTotalAbs(bool b){ fTotalAbs= b; }


private:
  // messeneger
  G4SBSDetectorMessenger* theMessenger;

  G4SBSGlobalField *fGlobalField;

  G4SBSBigBiteField *fbbfield;
  G4SBS48D48Field *f48d48field;
  
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

  bool fTotalAbs;

  int  fGEMOption;

};


#endif


