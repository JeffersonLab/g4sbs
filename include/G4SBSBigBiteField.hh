#ifndef G4SBSBigBiteField_hh
#define G4SBSBigBiteField_hh

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4RotationMatrix.hh"

#define MAXPT 120


class G4SBSBigBiteField : public G4SBSMappedField {
    public:
	G4SBSBigBiteField(G4ThreeVector , G4RotationMatrix *);
	~G4SBSBigBiteField();

	void GetFieldValue( const  double Point[3], double *Bfield ) const;
    private:
	void ReadField();
};

#endif//G4SBSBigBiteField_hh
