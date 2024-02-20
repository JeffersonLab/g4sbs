// SBS a Geant-4 Based Model of Hall-A Experiments with 11 GeV
// J.R.M Annand, University of Glasgow
// Class TOSCAField3D
// TOSCA generated field map
// Original by D.J.Hamilton
// 23/02/12 JRMA magnet "aperture" dimension changeable,
//               e.g. where field extends beyond aperture
// 27/02/12 JRMA implement x,y,z-offset to coordinates..covert to cm
//               change to c-style printf, scanf i/o
//

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

class TOSCAField3D : public G4MagneticField
{
  
  // Storage space for the table
  vector< vector< vector< double > > > xField;
  vector< vector< vector< double > > > yField;
  vector< vector< vector< double > > > zField;
  // The dimensions of the table
  int nx,ny,nz; 
  // The physical limits of the defined region
  double minx, maxx, miny, maxy, minz, maxz;
  // The physical extent of the defined region
  double dx, dy, dz;
  double fXoffset, fYoffset, fZoffset;
  bool invertX, invertY, invertZ;

public:
  TOSCAField3D(const char*, double, double, double, double=1.0 );
  void  GetFieldValue( const  double Point[4],
		       double *Bfield          ) const;
};

