// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class BonusBField
// Read Field map for Bonus solenoid (origin unknown). Phi-symmetric field
// 13/05/14 JRMA

#if !defined(BonusBField_H)
#define BonusBField_H

#include  <math.h>
#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

class BonusBField : public G4MagneticField
{
public:
  BonusBField(char*, double, double, double, double);
  virtual ~BonusBField();
  int jk(int j, int k) const { return j + k*Nr; } 
  virtual void GetFieldValue( const double Point[4],double *Bfield) const;
private:
  double* BrField;
  double* BzField;
  double Rmin,Rmax;
  double Zmin,Zmax;
  double Xoff,Yoff,Zoff;
  double ScaleFactor;
  double dR,dZ;
  int Nr,Nz;
};

#endif // !defined(BonusBField_H)
