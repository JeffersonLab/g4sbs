#ifndef G4SBS48D48Field_hh
#define G4SBS48D48Field_hh

#include "globals.hh"
#include "G4SBSMagneticField.hh"
#include "G4RotationMatrix.hh"

class G4SBS48D48Field : public G4SBSMappedField {
    public:
	G4SBS48D48Field(G4ThreeVector , G4RotationMatrix *);
	~G4SBS48D48Field();

	void GetFieldValue( const  double Point[3], double *Bfield ) const;

	bool GoodParams();

    private:
	void ReadField();
};

#endif//G4SBS48D48Field_hh
