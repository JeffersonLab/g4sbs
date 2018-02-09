#ifndef G4SBSDetectorConstruction_h
#define G4SBSDetectorConstruction_h 1

#include "G4LogicalVolume.hh"
#include "G4VUserDetectorConstruction.hh"
//#include "G4SBSIO.hh"
#include "sbstypes.hh"
#include "globals.hh"

#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4SDManager.hh"
#include <map>
#include <set>

using namespace std;

class G4SBSGlobalField;
class G4SBSMagneticField;
class G4SBSBigBiteField;
class G4SBS48D48Field;

class G4SBSBeamlineBuilder;
class G4SBSTargetBuilder;
class G4SBSEArmBuilder;
class G4SBSHArmBuilder;
class RTPC;

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

  G4SBSMagneticField *GetBBField(){ return fbbfield; }
  G4SBSMagneticField *Get48D48Field(){ return f48d48field; }
  G4SBSGlobalField *GetGlobalField(){ return fGlobalField; }

  void SetBigBiteField(int n);
  void Set48D48Field(int n);

  void SetTotalAbs(bool b){ fTotalAbs= b; }
  void SetCheckOverlap(bool b){ fCheckOverlap = b; }

  void SetExpType( Exp_t et ){ fExpType = et; }
  void SetTarget( Targ_t tg ){ fTargType = tg; }

  void SetUniformMagneticField48D48( double B );

  int  fBeamlineConf;
  int  fLeadOption;

  Exp_t fExpType;
  Targ_t fTargType;

  //map<G4String, G4VSensitiveDetector*> SDlist; //List of all sensitive detectors in the run. This is redundant with G4SDManager; just use that instead.
  set<G4String> SDlist;
  map<G4String, SDet_t> SDtype; //Mapping of sensitive detector names to sensitive detector types
  //map<G4String, Arm_t> SDarm; //Mapping of sensitive detector names to spectrometer arms
  set<G4String> StepLimiterList; //List of sensitive detectors for which G4UserLimits are defined to stop all particles entering (only allowed for calorimeters!)
  // map<G4String, G4double> SDgatewidth; //Time window for accumulating hit signal
  // map<G4String, G4double> SDthreshold; //threshold (energy deposition or photoelectrons) for recording a hit
  // map<G4String, G4int>    SDntimebins; //Time bins for "pulse shape" histogram
  
  
  G4SDManager *fSDman; 
  bool fTotalAbs;
  bool fCheckOverlap; //< Check if volumes overlap

  G4SBSBeamlineBuilder *fBeamlineBuilder;
  G4SBSTargetBuilder   *fTargetBuilder;
  G4SBSEArmBuilder     *fEArmBuilder;
  G4SBSHArmBuilder     *fHArmBuilder;
  RTPC *fRTPC;
  void SetBBDist( double);
  void SetBBAng( double);
  void Set48D48Dist( double);
  void Set48D48Ang( double);

  void AddToscaField(const char *);

  bool fUseGlobalField;

  G4SBSGlobalField *fGlobalField;

  void SetECALmapfilename( G4String );
  G4String GetECALmapfilename(){ return fECALmapfilename; }

  // int TrackerIDnumber;
  // map<int,Arm_t> TrackerArm; //Is tracker in E arm or P arm?

  void SetCDetconfig( int );
  int GetCDetConfigOption() { return fCDetOption; }

  void SetC16Segmentation( int );
  int GetC16Segmentation() { return fSegmentC16; }

  void SetFieldScale_SBS( G4double );
  G4double GetFieldScale_SBS() { return fFieldScale_SBS; }

  void SetFieldScale_BB( G4double );
  G4double GetFieldScale_BB() { return fFieldScale_BB; }

  void SetSegmentThickC16( G4double );
  G4double GetSegmentThickC16(){ return fSegmentThickC16; }

  void SetDoseRateC16( G4double );
  G4double GetDoseRateC16(){ return fDoseRateC16; }

  void SetFlipGEM( G4bool );
  G4bool GetFlipGEM(){ return fGEMflip; }
  
private:

  map<G4String, G4Material*> fMaterialsMap;
  map<G4String, G4OpticalSurface*> fOpticalSurfacesMap;

  G4SBSMagneticField *fbbfield;
  G4SBSMagneticField *f48d48field;
  
  G4String fECALmapfilename;

  //Let's define some additional configurable properties of 48D48:
  G4double f48D48_uniform_bfield; //set magnitude (and polarity) of SBS magnetic field (direction is fixed)
  G4double fFieldScale_SBS;
  G4double fFieldScale_BB;
  
  G4int fCDetOption;
  G4int fSegmentC16;

  G4double fSegmentThickC16;
  G4double fDoseRateC16; //Dose rate at z = 0 of lead-glass.

  G4bool fGEMflip;

};


#endif


