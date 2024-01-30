// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class RTPC
// Geometry and materials of the RTPC
// 03/05/14 JRMA
// 11/01/16 JRMA Update more "realistic" geometry

#ifndef RTPC_h
#define RTPC_h 1
#include "globals.hh"
#include "G4ThreeVector.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorConstruction;
class ArraySD;
//
class RTPC
{
public:
  RTPC(DetectorConstruction*);
  ~RTPC();
  virtual G4VPhysicalVolume* Construct(G4LogicalVolume * );
  virtual void ConEndCap();
  virtual void ConGEM();
  virtual void ConTarget();
  virtual void ConWire();
  virtual void ConGas();
  virtual void ConBeamLine();
  virtual void ConShell();
  virtual G4int ReadParameters( G4String );
  virtual void DefaultInit(){}
  virtual void SetMagField();

  //G4VPhysicalVolume* GetPhysi(){return fPCol;};
  G4LogicalVolume* GetLogic(){return 0;}
  G4LogicalVolume* GetMotherLogic(){return fMaw;}
  void BuildSrc();
  void SetIsInteractive(G4int is){fIsInteractive=is;}
  void SetIsOverlapVol( G4bool overlap ){ fIsOverlapVol = overlap; }
  void SetIsSrcPb( G4bool opt ){ fIsSrcPb = opt; }
  G4ThreeVector GetSrcPos(){ return G4ThreeVector(0,0,0); }
  G4ThreeVector GetSrcSize(){ return G4ThreeVector(0,0,0); }
  G4double GetTz(){ return fTz; }
  G4double GetZst(){ return fZst; }
  //TH3D* GetEloss();
  //TH3D* GetElossG();
  //TH3D* GetElossN();
  G4int* GetNhits();
  G4int* GetHitID();
  G4int* GetHits();

protected:
  DetectorConstruction* fRtag;
  ArraySD* fArrSD;
  G4int fNtarget;                       // 1 or 2....hydrogen or deuterium
  G4int fNgas;                          // 1 or 2....hydrogen or He
  G4double fTx, fTy, fTz;               // tank outer dimensions
  G4double fRst, fTst, fZst, fMst;      // target straw dimensions & material
  G4double fRst1;                       // outer radius straw ext. to beamline
  G4double fBz;                         // Z component of magnetic field
  G4double fZHe;                        // length of He in RTPC
  G4double fZrtpc;
  G4double fRHe1,fRHe1a,fRHe2,fRw,fOvHe;// He radii and wire radius
  //G4double fRblI,fRblO;               // inner outer beam line radius
  G4double fRbl,fTbl,fSbl,fZbl;         // inner outer beam line radius
  G4double fShI,fShO,fShZI,fShZO,fShTh; // outer shell dimensions
  G4double fWTh;                        // target window thickness
  G4int fWMat;                          // window material (Al or Be)
  G4double fBmClen1,fBmCr1,fBmClen2,fBmCr2; // collimators
  G4double fRbaf, fTbaf;                // radius & thickness of baffle
  G4double fTsh, fZsh;                  // downstream shield thickness & offset
  G4int fNwI,fNwO;                      //# field wires, inner and outer rings
  G4int fNbl;                           //# extra beamline segments after target
  G4int fIsSep;
  G4double fRG,fTG;                     // GEM space and outer can thickness
  G4double fRGtot;                      // total radius RTPC
  G4int fZPixG, fPhiPixG;               // GEM readout pads # pixels in Z & Phi
  G4double fTend1,fTend2;               // end cap foil thicknesses
  G4double fXmin, fBmin, fBmax;         // Non uniform field
  G4double fTXoff, fTYoff, fTZoff;      // Magnetic fieldmap offsets
  G4double fBScaleFac;                  // Scaling factor for magnetic field
  char* fFieldMap;                      // Field map file
  G4int fBFieldType;                    // model of B field
  G4double fTDens;                      // H2 target density
  G4double fHeDens;                     // He density
  G4int fVerbose;                       // verbose level
  G4int fIsInteractive;                 // batch(0) or interactive(1) mode
  G4LogicalVolume* fMaw;                // Logical volume of the mother
  G4LogicalVolume* fLrtpc;              // Logical volume of the RTPC
  G4LogicalVolume* fLHe2;               // Logical volume of RTPC gas  
  G4LogicalVolume* fLBfield;            // Logical volume of magnetic field
  G4VPhysicalVolume* fPWT;              // Physical volume for this detector
  G4Material* fMgas;                    // RTPC gas
  G4bool fIsOverlapVol;                 // if true check for overlaps
  G4bool fIsSrcPb;                      // is there a Pb shield around the src
};

#endif
