// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class TOSCAField2D
// TOSCA generated field. Phi-symmetric field
// 13/05/14 JRMA

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

class TOSCAField2D : public G4MagneticField
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
  TOSCAField2D(const char*, double, double, double, double=1.0 );
  void  GetFieldValue( const  double Point[4],
		       double *Bfield          ) const;
};

