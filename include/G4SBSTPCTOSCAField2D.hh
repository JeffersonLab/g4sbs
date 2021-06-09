// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class G4SBSTPCTOSCAField2D
// TOSCA generated field. Phi-symmetric field
// 13/05/14 JRMA
#ifndef G4SBSTPCTOSCAField2D_hh
#define G4SBSTPCTOSCAField2D_hh

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4String.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

class G4SBSTPCTOSCAField2D : public G4MagneticField
{
private:
  // Storage space for the table
  vector< vector< double > > BrField;
  vector< vector< double > > BzField;
  // The dimensions of the table
  unsigned int nx,ny,nz; 
  // The physical limits of the defined region
  double Rmin,Rmax,Zmin,Zmax;
  // The physical extent of the defined region
  double dR,dZ;
  double fXoffset, fYoffset, fZoffset;
public:
  // G4SBSTPCTOSCAField2D(const char*, double, double, double, double=1.0 );
  G4SBSTPCTOSCAField2D(G4double xoff, G4double yoff, G4double zoff, G4double scale=1.0 );
 ~G4SBSTPCTOSCAField2D();
  void  GetFieldValue( const  double Point[4],
		       double *Bfield          ) const;
};

#endif//G4SBSTPCTOSCAField2D_hh
