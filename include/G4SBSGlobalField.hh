#ifndef G4SBSGlobalField_hh
#define G4SBSGlobalField_hh

#include "globals.hh"
#include <vector>

#include "G4SBSMagneticField.hh"


class G4SBSGlobalField : public G4MagneticField {
    public:
	G4SBSGlobalField();
	~G4SBSGlobalField();

	void GetFieldValue( const  double Point[3], double *Bfield ) const;

	void SetInvertField( G4bool b );

	void AddField( G4SBSMagneticField *f );
    private:
	std::vector<G4SBSMagneticField *> fFields;

	bool fInverted;
};

#endif//G4SBSGlobalField_hh
