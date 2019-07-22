#ifndef G4SBSBigBiteField_hh
#define G4SBSBigBiteField_hh

#include "globals.hh"
#include "G4SBSMagneticField.hh"
#include "G4MagneticField.hh"
#include "G4RotationMatrix.hh"


class G4SBSBigBiteField : public G4SBSMappedField {
public:
  G4SBSBigBiteField(G4ThreeVector, G4RotationMatrix, G4String filename="map_696A.dat");
  ~G4SBSBigBiteField();

  void GetFieldValue( const G4double Point[3], G4double *Bfield ) const;
  
  G4int GetIndex( G4int, G4int, G4int ) const;  
  
private:
  void ReadField();

};

#endif//G4SBSBigBiteField_hh
