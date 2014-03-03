#ifndef __G4SBSHArmBuilder_hh
#define __G4SBSHArmBuilder_hh

#include "G4SBSComponent.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class G4SBSHArmBuilder: public G4SBSComponent {
    public:
	G4SBSHArmBuilder(G4SBSDetectorConstruction *);
	~G4SBSHArmBuilder();

	void BuildComponent(G4LogicalVolume *);

	void Set48D48Ang(double a){ f48D48ang = a; }
	void SetHCALDist(double a){ fHCALdist= a;   }
	void Set48D48Dist(double a){ f48D48dist = a; }
	void SetRICHdist( double d ){ fRICHdist = d; } //Set RICH detector distance
	void SetFieldClampConfig48D48( int option ){ f48D48_fieldclamp_config = option; }



	void Make48D48(G4LogicalVolume*, double);
	void MakeSBSFieldClamps(G4LogicalVolume*);
	void MakeHCAL(G4LogicalVolume*, G4double);
	void MakeFPP(G4LogicalVolume*, G4RotationMatrix*, G4ThreeVector );
	void MakeRICH(G4LogicalVolume *);

	double f48D48ang;
	double f48D48dist;
	int f48D48_fieldclamp_config; //Configuration of field clamp. There could be several of these.
	double fHCALdist;
	double fRICHdist;

	// Useful constants
	double f48D48depth;
	double f48D48height;
	double f48D48width;

	bool fUseLocalField;

	double fFieldStrength;
    private:

};

#endif//__G4SBSHArmBuilder_hh
