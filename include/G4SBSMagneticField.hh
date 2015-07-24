#ifndef G4SBSMagneticField_hh
#define G4SBSMagneticField_hh

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4RotationMatrix.hh"
#include "sbstypes.hh"
#include <vector>

using namespace std;

class G4SBSMagneticField : public G4MagneticField {
public:
  G4SBSMagneticField(G4ThreeVector, G4RotationMatrix);
  virtual ~G4SBSMagneticField();

  virtual void GetFieldValue( const  G4double Point[3], G4double *Bfield ) const = 0;

  void InvertField( G4bool b ){ fInverted = b; }

  void SetOffset( G4ThreeVector off ){ fOffset = off; }
  void SetRM( G4RotationMatrix rm ){ frm = rm; }

  G4bool fInverted;
  G4double fScaleFactor; //Set overall scale factor for any magnetic field

  Arm_t fArm; //kEarm or kHarm: defaults to Earm.
  
protected:
  
  G4ThreeVector fOffset;
  
  G4RotationMatrix frm;
};


////////////////////////////////////////////////////////////////////////////

class G4SBSMappedField: public G4SBSMagneticField {
public:
  G4SBSMappedField(G4ThreeVector, G4RotationMatrix,  const char *);
  virtual ~G4SBSMappedField();

protected:
  virtual void ReadField() = 0;

  char fFilename[255];

  G4int fN[3];
  G4int fNx, fNy, fNz; //Total points in x,y,z
  G4double fMin[3], fMax[3];
  //G4double fxmin,fxmax,fymin,fymax,fzmin,fzmax; //minimum, maximum x,y,z values
  //assumed to be a uniform, rectangular grid!
  //double ****fFieldVal;
  //how to define the grid:
  //vector<double> Bx, By, Bz; //array of field values; later convert to grid
  //map<G4int,G4ThreeVector> fBfield;
  vector<G4ThreeVector> fBfield;//store in a single vector, 1D array.
};

#endif//G4SBSMagneticField_hh
