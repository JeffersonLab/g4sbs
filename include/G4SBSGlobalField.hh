#ifndef G4SBSGlobalField_hh
#define G4SBSGlobalField_hh

#include "globals.hh"
#include "sbstypes.hh"
#include <vector>

#include "G4SBSMagneticField.hh"
#include "G4SystemOfUnits.hh"

class TH2F;

class G4SBSGlobalField : public G4MagneticField {
public:
  G4SBSGlobalField();
  ~G4SBSGlobalField();

  void GetFieldValue( const  double Point[3], double *Bfield ) const;

  void SetInvertField( G4bool b );

  void AddToscaField(const char *); //return a pointer to the added field so that we can store it in G4SBSDetectorConstruction

  void AddField( G4SBSMagneticField *f );
  void DropField( G4SBSMagneticField *f );

  void DebugField(G4double thEarm=33.0*deg, G4double thHarm=14.8*deg);

  void ScaleFields( G4double, G4SBS::Arm_t );
  
  std::vector<TH2F *> fFieldPlots;

  bool fInverted;

  //Utility method to write out a "local" field map for SBS, BigBite, or both, that can be rotated, translated and/or scaled with the magnet(s)
  //We will still have to use a global field definition with these maps, and we should write them out in the TOSCA format:
  //Parameters needed:
  // z min and z max (along spectrometer axis)
  // map height
  // map width
  // number of grid points along x, y, z (this will define grid spacing)
  // Grid spacing
  // number of 
  void WriteFieldMapSection( const char *fname, G4SBS::Arm_t arm, G4double theta, G4double zmin, G4double zmax, G4double h, G4double w, G4int nx, G4int ny, G4int nz );
  
  private:
  std::vector<G4SBSMagneticField *> fFields;
  //std::vector<G4double> fScaleFactor;
};

#endif//G4SBSGlobalField_hh
