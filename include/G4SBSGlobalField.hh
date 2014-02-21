#ifndef G4SBSGlobalField_hh
#define G4SBSGlobalField_hh

#include "globals.hh"
#include "G4MagneticField.hh"
#include "G4SBSMagneticField.hh"

#define MAXPT 120


class G4SBSGlobalField : public G4MagneticField {
    public:
	G4SBSGlobalField();
	~G4SBSGlobalField();

	void GetFieldValue( const  double Point[3], double *Bfield ) const;

	void SetInvertField( G4bool b );

	void AddField( G4SBSMagneticField *f );
    private:
	vector<G4SBSMagneticField *> fFields;

	bool fInverted;
};

#endif//G4SBSGlobalField_hh
