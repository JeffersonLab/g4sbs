#ifndef G4SBSToscaField_hh
#define G4SBSToscaField_hh

#include "globals.hh"
#include "G4SBSMagneticField.hh"
#include "G4RotationMatrix.hh"
#include "G4String.hh"

class G4SBSToscaField : public G4SBSMappedField {
public:
  G4SBSToscaField(G4String);
  ~G4SBSToscaField();

  void GetFieldValue( const  double Point[3], double *Bfield ) const;

  int GetIndex( int, int, int ) const;
  
private:
  void ReadField();
};

#endif//G4SBSToscaField_hh
