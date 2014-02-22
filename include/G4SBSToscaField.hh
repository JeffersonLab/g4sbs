#ifndef G4SBSToscaField_hh
#define G4SBSToscaField_hh

#include "globals.hh"
#include "G4SBSMagneticField.hh"
#include "G4RotationMatrix.hh"

class G4SBSToscaField : public G4SBSMappedField {
    public:
	G4SBSToscaField(const char *);
	~G4SBSToscaField();

	void GetFieldValue( const  double Point[3], double *Bfield ) const;

    private:
	void ReadField();
};

#endif//G4SBSToscaField_hh
