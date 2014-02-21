#ifndef G4SBSConstantField_hh
#define G4SBSConstantField_hh

#include "globals.hh"
#include "G4SBSMagneticField.hh"
#include "G4RotationMatrix.hh"

class G4SBSConstantField : public G4SBSMagneticField {
    public:
	G4SBSConstantField(G4ThreeVector , G4RotationMatrix *, G4ThreeVector, G4ThreeVector);
	~G4SBSConstantField();

	void GetFieldValue( const  double Point[3], double *Bfield ) const;
    private:

	G4ThreeVector fFieldVal;
	G4ThreeVector fBoxDim;
};

#endif//G4SBSConstantField_hh
