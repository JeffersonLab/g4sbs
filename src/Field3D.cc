// SBS a Geant-4 Based Model of Hall-A Experiments with 11 GeV
// J.R.M Annand, University of Glasgow
// Class Field3D
// TOSCA generated field map
// Original by D.J.Hamilton
// 23/02/12 JRMA magnet "aperture" dimension changeable,
//               e.g. where field extends beyond aperture
// 27/02/12 JRMA implement x,y,z-offset to coordinates..covert to cm
//               change to c-style printf, scanf i/o
//

#include "Field3D.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

Field3D::Field3D( G4double z, G4double bmin, G4double bmax) 
{
  fZmin = z;
  fBmin = bmin;
  fBmax = bmax;
}

void Field3D::GetFieldValue(const double point[4],
				      double *Bfield ) const
{
  // 2T/200mm
  //double x = point[0];
  //double y = point[1];
  double z = point[2];
  Bfield[0] = 0.0;
  Bfield[1] = 0.0;
  if( z < fZmin ) Bfield[2] = fBmax*tesla;
  else Bfield[2] = (fBmax + (z - fZmin)*0.5*(fBmax-fBmin)/fZmin)*tesla;
}

