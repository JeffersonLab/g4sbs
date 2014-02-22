#ifndef G4SBSDetectorConstruction_h
#define G4SBSDetectorConstruction_h 1

#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
#include "sbstypes.hh"
#include "globals.hh"

#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4SDManager.hh"
#include <map>

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

  void ConstructMaterials();  //Construct materials and optical surfaces
  G4Material *GetMaterial( G4String );
  G4OpticalSurface *GetOpticalSurface( G4String );

  //Build the various subsystems (BeamLine, Target, SBS Magnet, HCAL, BigCal, RICH, GEMs, etc.):
  void ConstructBeamline(G4LogicalVolume*);
  void ConstructTarget(G4LogicalVolume*);
  
  G4SBSBigBiteField *GetBBField(){ return fbbfield; }
  G4SBS48D48Field *Get48D48Field(){ return f48d48field; }

  void SetBigBiteField(int n);
  void Set48D48Field(int n);

  void SetTotalAbs(bool b){ fTotalAbs= b; }

  void SetExpType( Exp_t et ){ fExpType = et; }
  
  void SetUniformMagneticField48D48( double B ){ f48D48_uniform_bfield = B; }

private:
  // messeneger

  G4SDManager *fSDman; 

  map<G4String, G4Material*> fMaterialsMap;
  map<G4String, G4OpticalSurface*> fOpticalSurfacesMap;

  G4SBSBigBiteField *fbbfield;
  G4SBSMagneticField *f48d48field;
  


  //Let's define some additional configurable properties of 48D48:
  double f48D48_uniform_bfield; //set magnitude (and polarity) of SBS magnetic field (direction is fixed)




  bool fTotalAbs;

  Exp_t fExpType;

};


gicalVolume*, double);
  void MakeSBSFieldClamps(G4LogicalVolume*);
    void MakeHCAL(G4LogicalVolume*, G4double);

#endif


