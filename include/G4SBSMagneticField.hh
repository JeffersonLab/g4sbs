#ifndef G4SBSMagneticField_hh
#define G4SBSMagneticField_hh

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4RotationMatrix.hh"

class G4SBSMagneticField : public G4MagneticField {
    public:
	G4SBSMagneticField(G4ThreeVector, G4RotationMatrix *);
	virtual ~G4SBSMagneticField();

	virtual void GetFieldValue( const  double Point[3], double *Bfield ) const = 0;

	void InvertField( G4bool b ){ fInverted = b; }

	void SetOffset( G4ThreeVector off ){ fOffset = off; }
	void SetRM( G4RotationMatrix *rm ){ frm = rm; }

    protected:
	G4bool fInverted;

	G4ThreeVector fOffset;
	G4RotationMatrix *frm;
};


////////////////////////////////////////////////////////////////////////////

class G4SBSMappedField: public G4SBSMagneticField {
    public:
	G4SBSMappedField(G4ThreeVector, G4RotationMatrix *, const char *);
	virtual ~G4SBSMappedField();

    protected:
	virtual void ReadField() = 0;

	char fFilename[255];

	int fN[3];
	double fMin[3], fMax[3];

	double ****fFieldVal;
};

#endif//G4SBSMagneticField_hh
