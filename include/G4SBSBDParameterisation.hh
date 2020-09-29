// Beam Diffuser
// - Parameterization of the Hall A or C diffuser device 
// - Repeated volumes of varying size, separated by some distance
// - Geometry and materials based on JLab-TN-16-024 
// - Author: D. Flay (JLab) 

#ifndef G4SBS_BEAM_DIFFUSER_PARAMETERISATION_HH
#define G4SBS_BEAM_DIFFUSER_PARAMETERISATION_HH

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4PVParameterised.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

namespace diffuser {
   enum colorIndex {
      kWhite   = 0,
      kGray    = 1,
      kGrey    = 2,
      kBlack   = 3,
      kBrown   = 4,
      kRed     = 5,
      kGreen   = 6,
      kBlue    = 7,
      kCyan    = 8,
      kMagenta = 9,
      kYellow  = 10
   };
}

class G4SBSBDParameterisation: public G4VPVParameterisation {

   public:
      G4SBSBDParameterisation(char Hall='A',G4ThreeVector r0=G4ThreeVector(0,0,0) );
      virtual ~G4SBSBDParameterisation();

      // position and rotation 
      void ComputeTransformation(const G4int copyNo,G4VPhysicalVolume *physVol) const;

      // size of the object 
      virtual void ComputeDimensions(G4Tubs &plate,const G4int copyNo,const G4VPhysicalVolume *physVol) const;

      // shape of the object 
      // virtual G4VSolid *ComputeSolid(const G4int copyNo,G4VPhysicalVolume *physVol);

      // NOTE: Don't need this since the shape does not change 
      // material, sensitivity, visual attributes 
      // G4VTouchable should not be used for ordinary parameterisation
      // virtual G4Material *ComputeMaterial(const G4int copyNo,G4VPhysicalVolume *physVol,
      //                                     const G4VTouchable *parentTouch=0);

      // initialize private member arrays based on scalar members
      void InitParameters();

      // setters
      void SetHall(char h)        { fHall       = h; }
      void SetGap(double v)       { fGap        = v; }
      void SetWidth(double v)     { fWidth      = v; }
      void SetRadiusMin(double v) { fRadius_min = v; }
      void SetNPlanes(int v)      { fNPlanes    = v; }

      // getters 
      char     GetHall()        const { return fHall;       }
      G4int    GetNPlanes()     const { return fNPlanes;    }
      G4double GetGap()         const { return fGap;        }
      G4double GetWidth()       const { return fWidth;      }
      G4double GetRadiusMin()   const { return fRadius_min; }
      G4double GetRadiusMax()   const { return fRadius_max; }

   private:
      G4ThreeVector fR0;      // origin of device relative to mother volume 
      double fTotalThickness; // total thickness of assembly 
      double fGap;            // separation of the plates
      double fWidth;          // width of a plate 
      double fRadius_min;     // inner radius of a plate 
      double fRadius_max;     // outer radius of a plate (derived from width and inner radius)  
      double *fThickness;     // plate thicknesses
      double *fStartPhi;      // start angles 
      double *fDeltaPhi;      // step angles
      int    *fColor;         // plate colors  
      int    fNPlanes;        // number of planes 
      char   fHall;           // A or C 

};

#endif

