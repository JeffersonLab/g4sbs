#ifndef __G4SBSEArmBuilder_hh
#define __G4SBSEArmBuilder_hh

#include "G4SBSComponent.hh"

class G4SBSEArmBuilder: public G4SBSComponent {
    public:
	G4SBSEArmBuilder(G4SBSDetectorConstruction *);
	~G4SBSEArmBuilder();

	void BuildComponent(G4LogicalVolume *);

	void SetBBAng(double a);
	void SetBBDist(double a);

	void SetCerDepth(double a){ fCerDepth = a; }
	void SetCerDist(double a){fCerDist = a;}

	void SetGEMSep(double a){fGEMDist = a;}

	void SetBBCalDist(double a){ fBBCaldist= a; }
	void SetGEMConfig(int gc ){ fGEMOption = gc; }


    private:

	double fBBang;
	double fBBdist;
	double fBBCaldist;


	double fRICHdist; //distance from target of RICH detector

	double fCerDepth;
	double fCerDist;
	double fGEMDist;

	int  fGEMOption;

}

#endif//__G4SBSEArmBuilder_hh
