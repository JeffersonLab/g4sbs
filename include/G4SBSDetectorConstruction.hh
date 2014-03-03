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

using namespace std;

class G4SBSGlobalField;
class G4SBSMagneticField;
class G4SBSBigBiteField;
class G4SBS48D48Field;

class G4SBSBeamlineBuilder;
class G4SBSTargetBuilder;
class G4SBSEArmBuilder;
class G4SBSHArmBuilder;


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

	void SetExpType( Exp_t et ){ fExpType = et; }
	void SetTarget( Targ_t tg ){ fTargType = tg; }

	void SetUniformMagneticField48D48( double B );

	int  fLeadOption;

	Exp_t fExpType;
	Targ_t fTargType;

  map<G4String, G4VSensitiveDetector*> SDlist; //List of all sensitive detectors in the run.
	G4SDManager *fSDman; 
	bool fTotalAbs;

	G4SBSBeamlineBuilder *fBeamlineBuilder;
	G4SBSTargetBuilder   *fTargetBuilder;
	G4SBSEArmBuilder     *fEArmBuilder;
	G4SBSHArmBuilder     *fHArmBuilder;

	void SetBBDist( double);
	void SetBBAng( double);
	void Set48D48Dist( double);
	void Set48D48Ang( double);

	void AddToscaField(const char *);

	bool fUseGlobalField;

	G4SBSGlobalField *fGlobalField;

  int TrackerIDnumber;
  map<int,Arm_t> TrackerArm; //Is tracker in E arm or P arm?

    private:

	map<G4String, G4Material*> fMaterialsMap;
	map<G4String, G4OpticalSurface*> fOpticalSurfacesMap;

	G4SBSMagneticField *fbbfield;
	G4SBSMagneticField *f48d48field;



	//Let's define some additional configurable properties of 48D48:
	double f48D48_uniform_bfield; //set magnitude (and polarity) of SBS magnetic field (direction is fixed)



};


#endif


