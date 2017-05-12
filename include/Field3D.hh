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

using namespace std;

class Field3D : public G4MagneticField
{
protected:
  G4double fZmin;
  G4double fBmin, fBmax;
public:
  Field3D(double, double, double );
  void  GetFieldValue( const  G4double Point[4],
		       G4double *Bfield          ) const;
};

