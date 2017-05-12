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
// $Id: G4CalDetectorConstruction.cc 101905 2016-12-07 11:34:39Z gunter $
// 
/// \file G4CalDetectorConstruction.cc
/// \brief Implementation of the G4CalDetectorConstruction class

#include "G4CalDetectorConstruction.hh"

#include "G4SBSHArmBuilder.hh"
#include "G4SBSDetectorConstruction.hh"
#include "G4SBSRICHSD.hh"
#include "G4SBSGlobalField.hh"
#include "G4SBSTrackerBuilder.hh"
#include "G4SBSCalSD.hh"
#include "G4SBSECalSD.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"

#include "G4Trd.hh"
#include "G4TwoVector.hh"
#include "G4GenericTrap.hh"
#include "G4Polycone.hh"

#include "G4FieldManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <string>
#include <sstream>
#include "G4Element.hh"
#include "G4PhysicalConstants.hh"

#include "sbstypes.hh"

#include "G4PhysicalConstants.hh"

#include "TString.h"

#include "G4SystemOfUnits.hh"
using namespace std;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*G4ThreadLocal 
G4GlobalMagFieldMessenger* G4CalDetectorConstruction::fMagFieldMessenger = nullptr; */

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4CalDetectorConstruction::G4CalDetectorConstruction(G4SBSDetectorConstruction *dc):G4SBSComponent(dc)
   //fCheckOverlaps(true)
  //fCalAngle(30.*deg), fCalRotation(nullptr)
{
  fCalAngle=20*deg;
  fCalRotation = new G4RotationMatrix();
  fCalRotation->rotateY(fCalAngle);
 
  // define commands for this class
  //DefineCommands();//*******************

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4CalDetectorConstruction::~G4CalDetectorConstruction()
{ 
//delete fCalRotation;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4CalDetectorConstruction::BuildComponent(G4LogicalVolume *worldlog){

///oooooooooooooooje fais ci dessous si je contruit un detecteur s'apelle MAkeBigBite et pour des cas
  Exp_t exptype = fDetCon->fExpType;

  //  The neutron experiments and the SIDIS experiment use BigBite:
  //------------ BigBite: -----------------------------------------------------
  //if( exptype == kNeutronExp || exptype == kSIDISExp ) 
    //	{
  MakeCalo( worldlog);
   // }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//void G4SBSEArmBuilder::MakeBigBite(G4LogicalVolume *worldlog){

void G4CalDetectorConstruction::MakeCalo(G4LogicalVolume *worldlog)
{
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
// Liquid argon material(Si o veux rajouter un materiel avec des carateristiques)
 /* G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
 new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                kStateGas, 2.73*kelvin, 3.e-18*pascal);*/
  
  G4double a;  // atomic mass
  G4double zz;  // z=mean number of protons
  G4double density;
  G4int ncomponents;


  /*G4Element* elPb = new G4Element("Lead","Pb",zz=82.,a=207.2*g/mole);
  G4Element* elF = new G4Element("Fluor","F",zz=9.,a=18.998*g/mole);

  G4Material* PbF2 = new G4Material("PbF2",density=7.77*g/cm3,ncomponents=2);
  PbF2->AddElement(elPb,1);
  PbF2->AddElement(elF, 2);*/

  // Geometry parameters
  //G4int nofLayers = 1;
  //G4double absoThickness = 10.*mm;
  //G4double gapThickness =  5.*mm;

  G4double calorThickness =18.6*cm;
  G4double calorSizeXY  = 3.5*cm;

 /////ooooooooooooooooMother volume de calorimetreoooooooooooooooooooooooooo

  G4double worldSizeXY = 100. * calorSizeXY;// 
  G4double worldSizeZ  = 200. * calorThickness; 
  
  // Get materials
 auto defaultMaterial = G4Material::GetMaterial("Air");
  /*auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
  auto gapMaterial = G4Material::GetMaterial("liquidArgon");
  */
  if ( ! defaultMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("G4CalDetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
   
  //     
  // World
  //
  G4Box *worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2.0, worldSizeXY/2.0, worldSizeZ/2.0); // its size
                         
  G4LogicalVolume *worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
  G4double x = -2.0*m* std::sin(fCalAngle);
  G4double z = 2.0*m* std::cos(fCalAngle); 
                            
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(x,0.,z),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 worldlog,         // its mother  volume*************
                 false,            // no boolean operation
                 0,                // copy number
                 0);  // checking overlaps 
   
  ////////////////////////////////////////////////////////                               
  // 13 columns, 16 rows 
  double mylarthickness = 0.0020*cm, airthickness = 0.0040*cm;
  double mylar_air_sum = mylarthickness + airthickness; 
  double bbmodule_x = 3.5*cm, bbmodule_y = 3.5*cm;  
  double bbpmtz = 0.20*cm;

  double calheight = 16*3.5*cm;
  double calwidth  = 13*3.5*cm;
  double caldepth  = 18.6*cm + 2*bbpmtz + mylar_air_sum;

  G4Box *bbshowerbox = new G4Box("bbshowerbox", calwidth/2.0, calheight/2.0, caldepth/2.0);
  G4LogicalVolume *bbshowerlog = new G4LogicalVolume(bbshowerbox, defaultMaterial, "bbshowerlog");
  //new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth+bbhododepth+caldepth/2.0), bbshowerlog, "bbshowerphys", bbdetLog, false, 0);
  new G4PVPlacement( 0, G4ThreeVector(0,0,-30.0), bbshowerlog, "bbshowerphys", worldLV, false, 0 );
  

  // Calo module: 
  double bbPbF2_x = bbmodule_x - 2*mylar_air_sum;
  double bbPbF2_y = bbmodule_y - 2*mylar_air_sum;
  double bbPbF2_z = caldepth - 2*bbpmtz - mylar_air_sum;

  G4Box *calmodbox = new G4Box("calmodbox", bbmodule_x/2.0, bbmodule_y/2.0, caldepth/2.0);
  G4LogicalVolume *calmodlog = new G4LogicalVolume(calmodbox, GetMaterial("Special_Air"), "calmodlog");

  G4Box *tempbox = new G4Box("tempbox", bbmodule_x/2.0, bbmodule_y/2.0, (caldepth-2*bbpmtz)/2.0);


  // calorimeter box Subtraction
    G4Box *calmodbox_sub = new G4Box( "calmodbox_sub", (bbmodule_x-2*mylarthickness)/2.0, (bbmodule_y-2*mylarthickness)/2.0, (caldepth-2*bbpmtz)/2.0 );

  G4SubtractionSolid *bbmylarwrap = new G4SubtractionSolid( "bbmylarwrap", tempbox, calmodbox_sub, 0, G4ThreeVector(0.0, 0.0, mylarthickness) );
  G4LogicalVolume *bbmylarwraplog = new G4LogicalVolume( bbmylarwrap, defaultMaterial, "bbmylarwraplog" ); 
  
 // new G4LogicalSkinSurface( "BB Mylar Skin", bbmylarwraplog, GetOpticalSurface("Mirrsurf") );
  // Make Lead Glass 
  G4Box *bbPbF2box = new G4Box( "bbPbF2box", bbPbF2_x/2.0, bbPbF2_y/2.0, bbPbF2_z/2.0 );
  G4LogicalVolume *bbPbF2log = new G4LogicalVolume( bbPbF2box, GetMaterial("PbF2"), "bbPbF2log" );

  // Shower PbF2 SD of type CAL
  G4SDManager *sdman = fDetCon->fSDman;

  G4String BBSHPbF2SDname = "BBSHPbF2";
  G4String BBSHPbF2collname = "BBSHPbF2HitsCollection";
  G4SBSCalSD *BBSHPbF2SD = NULL;

  if( !((G4SBSCalSD*) sdman->FindSensitiveDetector(BBSHPbF2SDname)) ) {
    G4cout << "Adding G4Calo Shower PbF2 Sensitive Detector to SDman..." << G4endl;
    BBSHPbF2SD = new G4SBSCalSD( BBSHPbF2SDname, BBSHPbF2collname );
    sdman->AddNewDetector( BBSHPbF2SD );
    (fDetCon->SDlist).insert( BBSHPbF2SDname );
    fDetCon->SDtype[BBSHPbF2SDname] = kCAL;
    (BBSHPbF2SD->detmap).depth = 1;
  }
  bbPbF2log->SetSensitiveDetector( BBSHPbF2SD ); 

//////////////////

  if( (fDetCon->StepLimiterList).find( BBSHPbF2SDname ) != (fDetCon->StepLimiterList).end() ){
    bbPbF2log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }
  // Make PMT/Window
  double pmtrad = 1.5*cm;
  G4Tubs *bbPMT = new G4Tubs( "bbPMT", 0.0*cm, pmtrad, bbpmtz/2.0, 0.0, twopi );//geometri pmt
  G4LogicalVolume *bbpmtwindowlog = new G4LogicalVolume( bbPMT, GetMaterial("QuartzWindow_ECal"), "bbpmtwindowlog" );
  G4LogicalVolume *bbpmtcathodecallog = new G4LogicalVolume( bbPMT, GetMaterial("Photocathode_BB"), "bbpmtcathodecallog" );

  // Shower PMT SD of type ECAL
  G4String BBcalSDname = "BBcal";
  G4String BBcalcollname = "BBcalHitsCollection";
  G4SBSECalSD *BBcalSD = NULL;

  if( !((G4SBSECalSD*) sdman->FindSensitiveDetector(BBcalSDname)) ) {
    G4cout << "Adding BB Shower PMT Sensitive Detector to SDman..." << G4endl;
    BBcalSD = new G4SBSECalSD( BBcalSDname, BBcalcollname );
    sdman->AddNewDetector( BBcalSD );
    (fDetCon->SDlist).insert(BBcalSDname);
    fDetCon->SDtype[BBcalSDname] = kECAL;
    (BBcalSD->detmap).depth = 1;
  }
  bbpmtcathodecallog->SetSensitiveDetector( BBcalSD );

  // Put everything in a calo Module
  int shower_copy_number = 0;

  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-bbpmtz)/2.0), bbpmtcathodecallog,"bbcathodephys", calmodlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-3*bbpmtz)/2.0), bbpmtwindowlog, "bbwindowphys", calmodlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-4*bbpmtz-bbPbF2_z)/2.0), bbPbF2log, "bbPbF2phys", calmodlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, -bbpmtz), bbmylarwraplog, "bbmylarphys", calmodlog, false, 0 );

  int bbscol = 13;
  int bbsrow = 16;
  for( int l=0; l<bbscol; l++ ) {
    for( int j=0; j<bbsrow; j++ ) {

      (BBcalSD->detmap).Col[shower_copy_number] = l;
      (BBcalSD->detmap).Row[shower_copy_number] = j;
      (BBSHPbF2SD->detmap).Col[shower_copy_number] = l;
      (BBSHPbF2SD->detmap).Row[shower_copy_number] = j;
      double xtemp = (calwidth - bbmodule_x)/2.0 - l*bbmodule_x;
      double ytemp = (calheight - bbmodule_y)/2.0 - j*bbmodule_y;

 new G4PVPlacement(0, G4ThreeVector(xtemp,ytemp,0.0), calmodlog, "calphys", bbshowerlog, false, shower_copy_number);
      
      (BBcalSD->detmap).LocalCoord[shower_copy_number] = G4ThreeVector( xtemp,ytemp,(caldepth-bbpmtz)/2.0  );
      (BBSHPbF2SD->detmap).LocalCoord[shower_copy_number] = G4ThreeVector( xtemp, ytemp, (caldepth-4*bbpmtz-bbPbF2_z)/2.0 );

      shower_copy_number++;
    }
  }

///coloration
  //G4VisAttributes *mycalmodbox_colour = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0) );
  calmodlog->SetVisAttributes( G4VisAttributes::Invisible );

  //TF1
  G4VisAttributes *PbF2_colour = new G4VisAttributes(G4Colour( 1.0, 0.0, 0.0 ) );
  bbPbF2log->SetVisAttributes(PbF2_colour);

  //PMTcathode
  G4VisAttributes *PMT_colour = new G4VisAttributes(G4Colour( 0.0, 1.0, 0.0 ));
  bbpmtcathodecallog->SetVisAttributes(PMT_colour);

  
  //                                        
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  bbshowerlog->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  //return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

/*void G4CalDetectorConstruction::ConstructSDandField()
{ 
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);
  
  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
