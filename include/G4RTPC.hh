//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4RTPC.hh 75215 2013-10-29 16:07:06Z gcosmo $
// 
/// \file G4RTPC.hh
/// \brief Definition of the G4RTPC class

#ifndef G4RTPC_h
#define G4RTPC_h 1
#include "G4RotationMatrix.hh"
#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

//class G4VPhysicalVolume;
//class G4GlobalMagFieldMessenger;
class G4LogicalVolume;

/// Detector defconstruction class to define materials and geometry.
/// The calorimeter is a box made of a given number of layers. A layer consists
/// of an absorber plate and of a detection gap. The layer is replicated.
///<
/// Four parameters define the geometry of the calorimeter :
///
/// - the thickness of an absorber plate,
/// - the thickness of a gap,
/// - the number of layers,
/// - the transverse size of the calorimeter (the input face is a square).
///
/// In addition a transverse uniform magnetic field is defined 
/// via G4GlobalMagFieldMessenger class.

#include "G4SBSComponent.hh"

class G4RTPC : public G4SBSComponent
{
  public:
G4RTPC(G4SBSDetectorConstruction *);
G4RTPC();
  ~G4RTPC();


    void BuildComponent(G4LogicalVolume *);
   // void MakeRTPC(G4LogicalVolume *);
   void ConEndCap();
   void ConShell();
   void ConTarget();
   void ConGas();
   void ConWire();
   void ConGEM();
   void ConBeamLine();
   G4int ReadParameters( G4String );//***************
   void DefaultInit(){}
   void Beamline(G4LogicalVolume *);
   void SetMagField();
   // G4int ReadParameters();//***************
   // void DefaultInit(){};
    //virtual G4VPhysicalVolume* Construct();
    //virtual void ConstructSDandField();

    //void SetCalAngle(G4double val);//**************
   // G4double GetCalAngle() { return fCalAngle; }//**************
    
    
  private:
    // methods
    //
   //void DefineCommands();//*****************

    // data members
    //
   // static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger; 
                                      // magnetic field messenger
    G4double fCalAngle;//***********
    G4RotationMatrix* fCalRotation;//*************

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
    //G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps
};
   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

