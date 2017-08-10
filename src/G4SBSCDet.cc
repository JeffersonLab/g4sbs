#include "G4SBSCDet.hh"

#include "G4SBSHArmBuilder.hh"
#include "G4SBSDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4TwoVector.hh"
#include "G4GenericTrap.hh"
#include "G4Polycone.hh"

#include "G4SBSCalSD.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4SBSRICHSD.hh"
#include "G4SBSECalSD.hh"
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>

#include "sbstypes.hh"

#include "G4PhysicalConstants.hh"

#include "TString.h"

// To supress errors with TString, system of units should be included last
#include "G4SystemOfUnits.hh"

using namespace std;

G4SBSCDet::G4SBSCDet(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
  fR0 = 4050.0*cm;
  fZ0 = -426.5*cm;
  fPlanesHOffset = 0.0;
}

G4SBSCDet::~G4SBSCDet(){;}

void G4SBSCDet::BuildComponent(G4LogicalVolume *motherlog){
  MakeCDET( motherlog);
}

void G4SBSCDet::MakeCDET( G4LogicalVolume *mother ){
  //fR0 is the nominal distance from target to the start of CDET
  //fZ0 is the z position of the start of CDET relative to mother
  
  G4double Lx_scint = 51.0*cm;
  G4double Ly_scint = 0.5*cm;
  G4double Lz_scint = 4.0*cm;

  G4double HoleDiameter = 0.3*cm;
  G4double WLSdiameter      = 0.2*cm;
  G4double WLScladding_thick = 0.03*WLSdiameter/2.0;

  G4double mylar_thick = 0.25*0.001*2.54*cm; //.25 mil thickness of mylar
  //G4double mylar_thick = 0.1*mm;
  
  // G4int NColumns = 2;
  // G4int NRows    = 196;
  
  G4Box *Scint_strip = new G4Box("Scint_strip", (Lx_scint + mylar_thick)/2.0, Ly_scint/2.0 + mylar_thick, Lz_scint/2.0 + mylar_thick );
  G4Box *Scint_strip_nowrap = new G4Box("Scint_strip_nowrap", Lx_scint/2.0, Ly_scint/2.0, Lz_scint/2.0 );
  
  //Need to make a subtraction solid to define the mylar wrapping:
  // Need one (x) end to be open:
  G4Box *Scint_wrap_cutout = new G4Box("Scint_wrap_cutout", (Lx_scint + mylar_thick)/2.0 + 1.0*cm, Ly_scint/2.0, Lz_scint/2.0 );

  //Want to leave mylar_thick at +x end: cutout is 1 cm longer in x than module:
  // Lx_scint/2 + mylar_thick/2 - (x0 + Lx_scint/2 + mylar_thick/2 + 1.0 cm) = mylar_thick -->
  // x0 = mylar_thick - 1.0 cm
  G4SubtractionSolid *scint_wrap = new G4SubtractionSolid( "scint_wrap", Scint_strip, Scint_wrap_cutout, 0, G4ThreeVector( -(mylar_thick + 1.0*cm), 0, 0 ) );

  G4RotationMatrix *rot_fiber = new G4RotationMatrix;

  rot_fiber->rotateY( 90.0*deg );
  
  //This is to be used 
  G4Tubs *ScintHole = new G4Tubs( "ScintHole", 0.0, HoleDiameter/2.0, 1.01 * Lx_scint/2.0, 0.0, twopi );
  G4SubtractionSolid *Scint_strip_with_hole = new G4SubtractionSolid( "Scint_strip_with_hole", Scint_strip_nowrap, ScintHole, rot_fiber, G4ThreeVector(0,0,0) );
  
  G4Tubs *WLSfiber = new G4Tubs( "WLSfiber",                        0.0, 0.97*WLSdiameter/2.0, Lx_scint/2.0, 0.0, twopi );
  G4Tubs *WLScladding = new G4Tubs( "WLScladding", 0.97*WLSdiameter/2.0,      WLSdiameter/2.0, Lx_scint/2.0, 0.0, twopi );
  
  G4LogicalVolume *Scint_module = new G4LogicalVolume( Scint_strip, GetMaterial("Special_Air"), "Scint_module" );

  G4LogicalVolume *ScintWrapLog = new G4LogicalVolume( scint_wrap, GetMaterial("Mylar"), "ScintWrapLog" );
  //Make Mylar reflective:
  new G4LogicalSkinSurface( "CDET_mylar_wrap_osurf", ScintWrapLog, GetOpticalSurface("Mirrsurf") ); 
  
  G4LogicalVolume *ScintStripLog = new G4LogicalVolume( Scint_strip_with_hole, GetMaterial("CDET_BC408"), "ScintStripLog" );
  G4LogicalVolume *WLSFiberLog = new G4LogicalVolume( WLSfiber, GetMaterial("BCF_92"), "WLSFiberLog" );
  G4LogicalVolume *WLSCladdingLog = new G4LogicalVolume( WLScladding, GetMaterial("CDET_Acrylic"), "WLSCladdingLog" );

  new G4PVPlacement( 0, G4ThreeVector( 0,0,0 ), ScintWrapLog, "ScintWrapPhys", Scint_module, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( -mylar_thick/2.0, 0, 0 ), ScintStripLog, "ScintStripPhys", Scint_module, false, 0 );
  new G4PVPlacement( rot_fiber, G4ThreeVector( -mylar_thick/2.0, 0, 0 ), WLSFiberLog, "WLSFiberPhys", Scint_module, false, 0 );
  new G4PVPlacement( rot_fiber, G4ThreeVector( -mylar_thick/2.0, 0, 0 ), WLSCladdingLog, "WLSCladdingPhys", Scint_module, false, 0 );
  
  G4Tubs *CDET_pmt_cathode = new G4Tubs( "CDET_pmt_cathode", 0.0, WLSdiameter/2.0, 0.1*cm, 0.0, twopi );
  G4LogicalVolume *CDET_pmt_cathode_log = new G4LogicalVolume( CDET_pmt_cathode, GetMaterial("Photocathode_CDet"), "CDET_pmt_cathode_log" );

  G4SBSECalSD *cdet_sd = NULL;
  G4SDManager *sdman = fDetCon->fSDman;

  G4String sdname = "Earm/CDET";
  G4String collname = "CDETHitsCollection";
  
  if( !( cdet_sd = (G4SBSECalSD*) sdman->FindSensitiveDetector(sdname) ) ){
    G4cout << "Adding CDET sensitive detector to sdman..." << G4endl;
    cdet_sd = new G4SBSECalSD( sdname, collname );
    fDetCon->fSDman->AddNewDetector( cdet_sd );
    (fDetCon->SDlist).insert( sdname );
    fDetCon->SDtype[sdname] = kECAL;
    (cdet_sd->detmap).depth = 0;
    CDET_pmt_cathode_log->SetSensitiveDetector( cdet_sd );
  }

  sdname = "Earm/CDET_Scint";
  collname = "CDET_ScintHitsCollection";

  G4SBSCalSD *cdet_scint_sd = NULL;
  
  if( !( cdet_scint_sd = (G4SBSCalSD*) sdman->FindSensitiveDetector( sdname ) ) ){
    G4cout << "Adding CDET Scint sensitive detector to sdman..." << G4endl;
    cdet_scint_sd = new G4SBSCalSD( sdname, collname );
    fDetCon->fSDman->AddNewDetector( cdet_scint_sd );
    (fDetCon->SDlist).insert( sdname );
    fDetCon->SDtype[sdname] = kCAL;
    (cdet_scint_sd->detmap).depth = 1;
    ScintStripLog->SetSensitiveDetector( cdet_scint_sd );
  }
  
  //Now we need to define the coordinates of the "modules":
  //horizontal position within mother:
  G4double x0_modules[3] = { 0.5*(-70.986+31.014)*cm,
			     0.5*(-63.493+38.507)*cm,
			     0.5*(-56.0+46.0)*cm };
			     // -0.5*(-63.493+38.507)*cm,
			     // -0.5*(-70.986+31.014)*cm };

  //Number of rows per module:
  //G4int Nrow_module[3] = { 98, 98, 98, 98, 98, 98 };
  G4int Nrow_total = 588;

  //G4double y0[2][Nrow_total], pitch_angle[2][Nrow_total];

  //Nominal distance to CDET used to define the projective geometry:
  //G4double R0_CDET = 405.0*cm;
  //Nominal distance to planes:
  G4double R0_planes[2] = { fR0 + Lz_scint/2.0 + 1.0*cm,
			    fR0 + 3.0*Lz_scint/2.0 + 2.0*cm }; //allow for some small (1 cm) gaps between CH2 and start of 1st plane and between 1st and second planes...

  G4int istrip=0;
  
  for( int plane=0; plane<2; plane++ ){
    //step size in vertical angle:
    G4double dalpha = 2.0 * atan( (Ly_scint/2.0 + mylar_thick)/(R0_planes[plane] - Lz_scint/2.0 - mylar_thick) );
    for( int col=0; col<2; col++ ){
      //G4double ysum = Ly_scint/2.0;
      for( int row=0; row<294; row++ ){
      //for(int row=0; row<1; row++ ){
	G4double alpha = (row + 0.5 ) * dalpha;
	//for the first two strips, make them horizontal:
	G4int imod = 2 - row/98;	

	//hopefully, this expression won't lead to overlaps of strips?
	G4ThreeVector pos_strip( x0_modules[imod] + fPlanesHOffset/2.0*pow(-1, plane) + ( col - 0.5 )*(Lx_scint+mylar_thick), R0_planes[plane] * tan( alpha ), fZ0 + R0_planes[plane] - fR0 );

	G4RotationMatrix *rot_strip = new G4RotationMatrix;
	rot_strip->rotateY( col*pi );
	rot_strip->rotateX( alpha*pow(-1,col) );
	//rot_strip->rotateX( alpha );
	char physname[255];
	sprintf( physname, "Scint_module_phys_plane%d_row%d_col%d", plane+1, row+295, col+1 ); //In this construction, row varies from 295-588 for the top half
	G4String pvolname = physname;
	new G4PVPlacement( rot_strip, pos_strip, Scint_module, pvolname, mother, false, istrip );

	G4ThreeVector pos_pmt( x0_modules[imod] + fPlanesHOffset/2.0*pow(-1, plane) + pow(-1,col+1)*(Lx_scint+mylar_thick+0.1*cm), R0_planes[plane] * tan( alpha ), fZ0 + R0_planes[plane] - fR0 );
	sprintf( physname, "CDET_pmt_phys_plane%d_row%d_col%d", plane+1, row+295, col+1 );
	pvolname = physname;
	new G4PVPlacement( rot_fiber, pos_pmt, CDET_pmt_cathode_log, pvolname, mother, false, istrip );

	(cdet_sd->detmap).Row[istrip] = row+295;
	(cdet_sd->detmap).Col[istrip] = col+1;
	(cdet_sd->detmap).Plane[istrip] = plane+1;
	(cdet_sd->detmap).LocalCoord[istrip] = pos_strip;

	(cdet_scint_sd->detmap).Row[istrip] = row+295;
	(cdet_scint_sd->detmap).Col[istrip] = col+1;
	(cdet_scint_sd->detmap).Plane[istrip] = plane+1;
	(cdet_scint_sd->detmap).LocalCoord[istrip] = pos_pmt;
	
	istrip++;

	//Next: make bottom half:

	alpha *= -1.0;

	pos_strip.setY( R0_planes[plane] * tan(alpha) );

	G4RotationMatrix *rot_strip2 = new G4RotationMatrix;
	rot_strip2->rotateY( col * pi );
	rot_strip2->rotateX( alpha*pow(-1,col) );

	sprintf( physname, "Scint_module_phys_plane%d_row%d_col%d", plane+1, 294-row, col+1 ); //In this construction, row varies from 294 down to 1 for the bottom half:

	pvolname = physname;

	new G4PVPlacement( rot_strip2, pos_strip, Scint_module, pvolname, mother, false, istrip );

	pos_pmt.setY( R0_planes[plane] * tan(alpha) );
	sprintf( physname, "CDET_pmt_phys_plane%d_row%_col%d", plane+1, 294-row, col+1 );
	pvolname = physname;
	new G4PVPlacement( rot_fiber, pos_pmt, CDET_pmt_cathode_log, pvolname, mother, false, istrip );

	(cdet_sd->detmap).Row[istrip] = 294-row;
	(cdet_sd->detmap).Col[istrip] = col+1;
	(cdet_sd->detmap).Plane[istrip] = plane+1;
	(cdet_sd->detmap).LocalCoord[istrip] = pos_strip;

	(cdet_scint_sd->detmap).Row[istrip] = 294-row;
	(cdet_scint_sd->detmap).Col[istrip] = col+1;
	(cdet_scint_sd->detmap).Plane[istrip] = plane+1;
	(cdet_scint_sd->detmap).LocalCoord[istrip] = pos_pmt;
	
	istrip++;
	
      }
    }
  }

  Scint_module->SetVisAttributes( G4VisAttributes::Invisible );
  
  G4VisAttributes *scintstrip_visatt = new G4VisAttributes( G4Colour( 0.8, 0, 0.8 ) );
  ScintStripLog->SetVisAttributes( scintstrip_visatt );

  G4VisAttributes *WLSfiber_visatt = new G4VisAttributes( G4Colour( 0, 0.8, 0.8 ) );
  WLSFiberLog->SetVisAttributes( WLSfiber_visatt );

  G4VisAttributes *scintwrap_visatt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0 ) );
  ScintWrapLog->SetVisAttributes( scintwrap_visatt );
  
}
