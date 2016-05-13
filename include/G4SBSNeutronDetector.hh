#ifndef G4SBSNEUTRONDETECTOR_HH
#define G4SBSNEUTRONDETECTOR_HH

#include "G4SBSComponent.hh"
#include <vector>

#define NBARTYPES 6
#define NCASSETTE_TYPES 6
#define NLAYERS 9
#define NPLATES 6
#define MAX_CASS 10

class G4LogicalVolume;

class G4SBSNeutronDetector : public G4SBSComponent {

public :
  G4SBSNeutronDetector(G4SBSDetectorConstruction *);
  ~G4SBSNeutronDetector(){;}

  enum BarType_t { kBarCMU = 0, kBarUVA, kBarJLab, kBarGlasgow, kBarVetoLong, kBarVetoShort };
  enum CassType_t { kCassCMU = 0, kCassUVA, kCassJLab, kCassGlasgow, kCassVeto1, kCassVeto2, kCassNull };

  void  BuildComponent(G4LogicalVolume *); 
   
  G4int GetNRows(int i) { return NRows[i]; } 
  void SetNDDist( double d ){ fNDdist = d; }
  void SetNDAngle( double a ){ fNDang = a; }
  void SetNDTarg( int a ) {fTarget = a; }

private:
    
  G4LogicalVolume* ConstructND(G4LogicalVolume*);     

  G4int NRows[NLAYERS];
  G4LogicalVolume* logicPMT[NLAYERS];

  G4LogicalVolume *logBlock[NBARTYPES][NLAYERS];
  G4LogicalVolume *logBar[NBARTYPES][NLAYERS];

  G4LogicalVolume *logCassette[NCASSETTE_TYPES][NLAYERS];

  // Geometry descriptors
  double X[NBARTYPES], Y[NBARTYPES], Z[NBARTYPES], attlen[NBARTYPES];
  double lgSize[NBARTYPES], lgNearDepth[NBARTYPES], lgFarDepth[NBARTYPES], lgLength[NBARTYPES], lgcylLen[NBARTYPES];

  double casX[NCASSETTE_TYPES], casY[NCASSETTE_TYPES], casZ[NCASSETTE_TYPES];
  int Nbars[NCASSETTE_TYPES];
  double barspacing[NCASSETTE_TYPES], casbarZ[NCASSETTE_TYPES], casAlthick[NCASSETTE_TYPES];
  double  casFe1thick[NCASSETTE_TYPES], casFe2thick[NCASSETTE_TYPES];
  double  casFe1length[NCASSETTE_TYPES], casFe2length[NCASSETTE_TYPES], casFe2height[NCASSETTE_TYPES];

  std::vector<G4Material*> plateMat;
  int Ncass[NLAYERS];
  double layerBottom[NLAYERS], layerZ[NLAYERS], pmtRad[NLAYERS];
  double plateZ[NPLATES], plateThick[NPLATES];

  double plateHeight, plateWidth;

  double casbarY[NCASSETTE_TYPES];
  CassType_t casType[NLAYERS][MAX_CASS];

  double fNDdist;
  double fNDang;
  double NDdistance;
  double NDangle;

  double fThreshold[NBARTYPES];
  int fTarget;
};
#endif
