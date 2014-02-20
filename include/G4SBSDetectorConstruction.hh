#ifndef G4SBSDetectorConstruction_h
#define G4SBSDetectorConstruction_h 1

#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "sbstypes.hh"
#include "globals.hh"

#include "G4SBSBigBiteField.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4SDManager.hh"
#include <map>

using namespace std;
//class G4SBSDetectorMessenger;

class G4SBSDetectorConstruction : public G4VUserDetectorConstruction
{
public:
  G4SBSDetectorConstruction();
  ~G4SBSDetectorConstruction();
  
public:
  G4VPhysicalVolume* Construct();
  G4VPhysicalVolume* ConstructAll();
  //G4VPhysicalVolume* ConstructAllGEp();
  //  G4VPhysicalVolume* ConstructAllSIDIS();

  void ConstructMaterials();  //Construct materials and optical surfaces
  G4Material *GetMaterial( G4String );
  G4OpticalSurface *GetOpticalSurface( G4String );

  //Build the various subsystems (BeamLine, Target, SBS Magnet, HCAL, BigCal, RICH, GEMs, etc.):
  void ConstructBeamline(G4LogicalVolume*);
  void ConstructTarget(G4LogicalVolume*);
  void Make48D48(G4LogicalVolume*, double);
  void MakeBigBite(G4LogicalVolume*);
  void MakeHCAL(G4LogicalVolume*, G4double);
  void MakeBigCal(G4LogicalVolume*);
  void MakeTracker(G4LogicalVolume*, G4RotationMatrix*, G4ThreeVector, G4int, vector<double>, vector<double>, vector<double> );
  void MakeFPP(G4LogicalVolume*, G4RotationMatrix*, G4ThreeVector );
  void MakeRICH(G4LogicalVolume *);
  

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

  void SetTotalAbs(bool b){ fTotalAbs= b; }

  void SetRICHdist( double d ){ fRICHdist = d; } //Set RICH detector distance

  void SetExpType( Exp_t et ){ fExpType = et; }

private:
  // messeneger
  //G4SBSDetectorMessenger* theMessenger;

  G4SDManager *SDman; //Due to splitting things into different routines, we need to store a pointer to the G4SDManager.

  map<G4String, G4Material*> MaterialsMap;
  map<G4String, G4OpticalSurface*> OpticalSurfacesMap;

  G4SBSBigBiteField *fbbfield;
  
  double fBBang;
  double fBBdist;
  double fBBCaldist;

  double f48D48ang;
  double f48D48dist;

  //Let's define some additional configurable properties of 48D48:
  double f48D48_uniform_bfield; //set magnitude (and polarity) of SBS magnetic field.
  int f48D48fieldclamp_config; //Configuration of field clamp. There could be several of these.

  double fHCALdist;

  double fRICHdist; //distance from target of RICH detector

  double fCerDepth;
  double fCerDist;
  double fGEMDist;

  double fTargLen;
  double fTargDen;

  Targ_t fTargType;

  bool fTotalAbs;

  int  fGEMOption;

  Exp_t fExpType;

};


#endif


