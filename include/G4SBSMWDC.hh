#ifndef G4SBSMWDC_HH
#define G4SBSMWDC_HH

#include "G4SBSComponent.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4VisAttributes.hh"

#include <vector>
#include <map>

class G4SBSMWDC : public G4SBSComponent {
public:
  G4SBSMWDC(G4SBSDetectorConstruction *);
  ~G4SBSMWDC();

  void BuildComponent(G4LogicalVolume *);
  void BuildComponent(G4LogicalVolume *, G4RotationMatrix *, 
		      G4ThreeVector, G4String );  

  G4LogicalVolume* BuildX(double,double,int,int);
  G4LogicalVolume* BuildUorV(double,double,G4String,int,int);

private:

  // GEn 2008 Data from Seamus' Thesis
  std::vector<int> fChamberNumber;
  std::vector<int> fNplanes;
  std::vector<int> fNwires;
  std::vector<double> fNwirespacing;
  std::vector<double> fNheight;
  std::vector<double> fNwidth;
  std::vector<double> fDist_z0;

  // Key is chamber #, the vector corresponds to the
  // plane type (X,U,V)
  std::vector<G4String> fChamber0, fChamber1, fChamber2;
  std::map<int, std::vector<G4String> > fGEn_Setup; 
  std::map<int, G4LogicalVolume*> fCathodes;
  
  double fFieldD;        // diameter of field wire 
  double fSignalD;       // diameter of signal wire
  double fWireSep;       // distance b/t adjacent field wires
  double fPlaneThick;    // thickness of one plane (X,U,V)
  double fCathodeThick;  // thickness of mylar + 2*copper cathode
  double fMylarThick;    // thickness of mylar wrt cathode
  double fCuThick;       // thickness of cu wtt cathode

  double fCath2WireDist;
  double fGlassThick;

  // VISUALS:
  G4VisAttributes* mylarVisAtt;
  G4VisAttributes* cuVisAtt;
  G4VisAttributes* glassVisAtt;
  G4VisAttributes* gasVisAtt;
};
#endif
