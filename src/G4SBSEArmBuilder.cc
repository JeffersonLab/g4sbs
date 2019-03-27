#include "G4SBSEArmBuilder.hh"

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

#include "G4SBSECal.hh"
#include "G4SBSGrinch.hh"
#include "G4SBSBigBiteField.hh"
#include "G4SBSTrackerBuilder.hh"
#include "G4SBSCalSD.hh"
#include "G4SBSGlobalField.hh"

#include "G4FieldManager.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"
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

G4SBSEArmBuilder::G4SBSEArmBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
  fBBang  = 40.0*deg;
  fBBdist = 1.5*m;

  /*
  G4double frontGEM_depth = 20.*cm;
  G4double backGEM_depth = 10.*cm;
  
  fCerDepth = 92.0*cm;
  //fCerDist  =  7.0*cm;
  //fCerDist = 22.0*cm;
  fCerDist = frontGEM_depth + 2.*cm;
  
  // fBBCaldist = 20*cm + fCerDepth;
  // fGEMDist   = 10*cm + fCerDepth;
  // fGEMOption = 1;

  fBBCaldist = fCerDist + fCerDepth + backGEM_depth + 5.*cm;
  fGEMDist   = fCerDist + fCerDepth + 0.5*backGEM_depth;
  fGEMOption = 1;
  */
  G4double frontGEM_depth = 60.*cm;
  G4double backGEM_depth = 11.59*cm;
  
  fCerDepth = 88.9*cm;
  fCerDist = frontGEM_depth - 8.571*cm + 1.811*cm;
  
  //NB: fBBCalDist now designates the distance to the shielding
  fBBCaldist = fCerDist + fCerDepth + backGEM_depth;
  fGEMDist   = fCerDist + fCerDepth + 0.5*backGEM_depth;
  fGEMOption = 2;
  /**/
  fShieldOption = 1;
  
  fUseLocalField = false;

  fnzsegments_leadglass_ECAL = 1;
  fnzsegments_leadglass_C16 = 1;

  fBuildBBSieve = false;
  
  assert(fDetCon);
  
  fbbfield =  NULL;

  fGRINCHgas = "C4F10_gas"; //default to C4F10;
  
}

G4SBSEArmBuilder::~G4SBSEArmBuilder(){;}

void G4SBSEArmBuilder::BuildComponent(G4LogicalVolume *worldlog){
  Exp_t exptype = fDetCon->fExpType;
  G4SBSECal* ECal = fDetCon->fECal;
  
  //  The neutron experiments and the SIDIS experiment use BigBite:
  //------------ BigBite: -----------------------------------------------------
  if( exptype == kNeutronExp || exptype == kSIDISExp || exptype == kA1n || exptype == kGEnRP ) 
    {
      MakeBigBite( worldlog );
      if(fBuildBBSieve)
	MakeBBSieveSlit(worldlog);
      ECal->~G4SBSECal();
    }
  if( exptype == kGEp ) //Subsystems unique to the GEp experiment include FPP and BigCal:
    {
      ECal->SetAng(fBBang);
      ECal->SetDist(fBBdist);
      ECal->BuildComponent(worldlog);
      //MakeBigCal( worldlog );
    }
  if( exptype == kC16 ) 
    {
      ECal->SetAng(fBBang);
      ECal->SetDist(fBBdist);
      ECal->BuildComponent(worldlog);
      //MakeC16( worldlog );
    }
  if( exptype == kNeutronExp || exptype == kGEnRP )  MakeGMnGEMShielding( worldlog );
  
  if( exptype == kNDVCS ){
    ECal->SetAng(fBBang);
    ECal->SetDist(fBBdist);
    ECal->BuildComponent(worldlog);
  }
}

void G4SBSEArmBuilder::MakeBigBite(G4LogicalVolume *worldlog){
  //Lines of code used to build BigBite moved to their own method:
  printf("BigBite at %f deg\n", fBBang/deg);

  G4RotationMatrix *bbrm = new G4RotationMatrix;
  bbrm->rotateY(-fBBang);

  G4RotationMatrix *bbykrm = new G4RotationMatrix;
  bbykrm->rotateX(90.0*deg);

  double motherdepth = 600.0*cm;

  double bbmagwidth  = 1670.0*mm; 
  double bbmagheight = 486.0*mm;

  double eps = 0.1*mm;

  // Mother box
  G4Box *bbmotherBox= new G4Box("bbmotherBox", bbmagwidth/2.0 + eps, 250*cm, motherdepth/2.0);

  // We need to account for offsets so we can fit BigBite and detectors in without running into
  // the target
  // Need 70 cm clearance from front of the spectrometer to front of the mother volume

  double clear = 70.0*cm;

  double motherr = fBBdist + motherdepth/2.0 - clear;

  std::vector<G4TwoVector> bbyokepts;
  bbyokepts.push_back( G4TwoVector( bbmagheight, 0.0*mm));
  bbyokepts.push_back( G4TwoVector( bbmagheight, 464.0*mm));
  bbyokepts.push_back( G4TwoVector( -bbmagheight, 805.0*mm));
  bbyokepts.push_back( G4TwoVector( -bbmagheight, 0.0*mm));
  bbyokepts.push_back( G4TwoVector( bbmagheight, 0.0*mm));
  bbyokepts.push_back( G4TwoVector( bbmagheight, 464.0*mm));
  bbyokepts.push_back( G4TwoVector( -bbmagheight, 805.0*mm));
  bbyokepts.push_back( G4TwoVector( -bbmagheight, 0.0*mm));

  G4GenericTrap *bbyokeTrap = new G4GenericTrap("bbyokeTrap",
						bbmagwidth/2.0, bbyokepts );


  std::vector<G4TwoVector> bbairpts;
  bbairpts.push_back( G4TwoVector(  bbmagheight-133.1*mm, 0.0*mm - eps));
  bbairpts.push_back( G4TwoVector(  bbmagheight, 464.1*mm));
  bbairpts.push_back( G4TwoVector( -bbmagheight, 805.1*mm));
  bbairpts.push_back( G4TwoVector( -bbmagheight, 0.0*mm - eps));

  bbairpts.push_back( G4TwoVector(  bbmagheight-133.1*mm, 0.0*mm - eps));
  bbairpts.push_back( G4TwoVector(  bbmagheight, 464.1*mm));
  bbairpts.push_back( G4TwoVector( -bbmagheight, 805.1*mm));
  bbairpts.push_back( G4TwoVector( -bbmagheight, 0.0*mm - eps));

  double gapsize = 250.0*mm;

  G4GenericTrap *bbairTrap = new G4GenericTrap("bbairTrap",
					       gapsize/2.0, bbairpts );

  double coilsize = 320.0*mm;
  double coilwidth = 90.0*mm;

  std::vector<G4TwoVector> bbcoilpts;
  bbcoilpts.push_back( G4TwoVector(  bbmagheight-133.0*mm+coilsize, 0.0*mm-coilsize));
  bbcoilpts.push_back( G4TwoVector(  bbmagheight+coilsize, 464.0*mm+coilsize*(1.0-tan(20.0*deg))));
  bbcoilpts.push_back( G4TwoVector( -bbmagheight-coilsize, 805.0*mm+coilsize*(1.0+tan(20.0*deg))));
  bbcoilpts.push_back( G4TwoVector( -bbmagheight-coilsize, 0.0*mm-coilsize));

  bbcoilpts.push_back( G4TwoVector(  bbmagheight-133.0*mm+coilsize, 0.0*mm-coilsize));
  bbcoilpts.push_back( G4TwoVector(  bbmagheight+coilsize, 464.0*mm+coilsize*(1.0-tan(20.0*deg))));
  bbcoilpts.push_back( G4TwoVector( -bbmagheight-coilsize, 805.0*mm+coilsize*(1.0+tan(20.0*deg))));
  bbcoilpts.push_back( G4TwoVector( -bbmagheight-coilsize, 0.0*mm-coilsize));


  G4GenericTrap *bbcoilTrap = new G4GenericTrap("bbyokeTrap",
						coilwidth, bbcoilpts );

  double coilrot = 40.0*mrad;

  G4RotationMatrix *coilrot1 = new G4RotationMatrix;
  coilrot1->rotateX(-coilrot );
  G4RotationMatrix *coilrot2 = new G4RotationMatrix;
  coilrot2->rotateX( coilrot );

  G4ThreeVector zcoil1 = G4ThreeVector( 0.0, 0.0, coilwidth+gapsize/2.0+(805.0*mm/2.0+coilsize)*tan(coilrot));
  G4ThreeVector zcoil2 = G4ThreeVector( 0.0, 0.0, -coilwidth-gapsize/2.0-(805.0*mm/2.0+coilsize)*tan(coilrot));

  // "Guitar" cuts
  std::vector<G4TwoVector> bbleftcutpts;
  bbleftcutpts.push_back( G4TwoVector(  0.0, ((bbmagwidth+gapsize)/2.0+coilwidth)/2.0 ));
  bbleftcutpts.push_back( G4TwoVector(  bbmagheight,  bbmagwidth/2.0));
  bbleftcutpts.push_back( G4TwoVector(  0.0,  bbmagwidth/2.0 + 1.0*cm));
  bbleftcutpts.push_back( G4TwoVector(  -bbmagheight*2.0/3.0,  bbmagwidth/2.0));

  bbleftcutpts.push_back( G4TwoVector(  0.0, ((bbmagwidth+gapsize)/2.0+coilwidth)/2.0 ));
  bbleftcutpts.push_back( G4TwoVector(  bbmagheight,  bbmagwidth/2.0));
  bbleftcutpts.push_back( G4TwoVector(  0.0,  bbmagwidth/2.0 + 1.0*cm));
  bbleftcutpts.push_back( G4TwoVector(  -bbmagheight*2.0/3.0,  bbmagwidth/2.0));


  G4GenericTrap *bbleftcutTrap = new G4GenericTrap("bbleftcutTrap",
						   2010.1*mm, bbleftcutpts );

  G4RotationMatrix *leftcutrot = new G4RotationMatrix;
  leftcutrot->rotateX( -90*deg );

  G4RotationMatrix *leftcutrot2 = new G4RotationMatrix;
  leftcutrot2->rotateZ( -90*deg );

  std::vector<G4TwoVector> bbrightcutpts;
  bbrightcutpts.push_back( G4TwoVector(  0.0, -((bbmagwidth+gapsize)/2.0+coilwidth)/2.0 ));
  bbrightcutpts.push_back( G4TwoVector(  bbmagheight,  -bbmagwidth/2.0));
  bbrightcutpts.push_back( G4TwoVector(  0.0, -( bbmagwidth/2.0 + 1.0*cm)));
  bbrightcutpts.push_back( G4TwoVector(  -bbmagheight*2.0/3.0,  -bbmagwidth/2.0));

  bbrightcutpts.push_back( G4TwoVector(  0.0, -((bbmagwidth+gapsize)/2.0+coilwidth)/2.0 ));
  bbrightcutpts.push_back( G4TwoVector(  bbmagheight,  -bbmagwidth/2.0));
  bbrightcutpts.push_back( G4TwoVector(  0.0, -( bbmagwidth/2.0 + 1.0*cm)));
  bbrightcutpts.push_back( G4TwoVector(  -bbmagheight*2.0/3.0,  -bbmagwidth/2.0));

  G4GenericTrap *bbrightcutTrap = new G4GenericTrap("bbrightcutTrap",
						    2010.1*mm, bbrightcutpts );

  G4RotationMatrix *rightcutrot = new G4RotationMatrix;
  rightcutrot->rotateX( 90*deg );

  G4RotationMatrix *rightcutrot2 = new G4RotationMatrix;
  rightcutrot2->rotateZ( -90*deg );

  G4UnionSolid *fullyoke;

  fullyoke = new G4UnionSolid("yoke_coil1", bbyokeTrap, bbcoilTrap, coilrot1, zcoil1 );
  fullyoke = new G4UnionSolid("yoke_coils", fullyoke, bbcoilTrap, coilrot2, zcoil2 );

  double topheight = 66.0*cm;
  G4Box *yoketopbox = new G4Box("yoketopbox", topheight/2.0, 464.0*mm/2.0, bbmagwidth/2.0);
  G4ThreeVector topboxpos = G4ThreeVector( topheight/2.0 + bbmagheight, 464.0*mm/2.0, 0.0);

  double bottomheight = 82.0*cm;
  G4Box *yokebottombox = new G4Box("yokebotbox", bottomheight/2.0, 805.0*mm/2.0, bbmagwidth/2.0);
  G4ThreeVector bottomboxpos = G4ThreeVector( -bottomheight/2.0 - bbmagheight, 805.0*mm/2.0, 0.0);

  fullyoke = new G4UnionSolid("yoke_top", fullyoke, yoketopbox, 0, topboxpos );
  fullyoke = new G4UnionSolid("yokefull", fullyoke, yokebottombox, 0, bottomboxpos );

  G4SubtractionSolid* yokewgap = new G4SubtractionSolid("yoke_with_gap", fullyoke, bbairTrap);
  yokewgap = new G4SubtractionSolid("yoke_with_gap_lcut", yokewgap, bbleftcutTrap, leftcutrot, G4ThreeVector());
  yokewgap = new G4SubtractionSolid("yoke_with_gap_lrcut", yokewgap, bbrightcutTrap, leftcutrot, G4ThreeVector());

  G4LogicalVolume *bbyokewgapLog=new G4LogicalVolume(yokewgap, GetMaterial("Fer"),
						     "bbyokewgapLog", 0, 0, 0);

  if( fDetCon->fTotalAbs ){
    bbyokewgapLog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }

  G4RotationMatrix *yokerm = new G4RotationMatrix;
  yokerm->rotateY(90.0*deg);
  yokerm->rotateZ(-90.0*deg);
  yokerm->rotateX(180.0*deg);

  // Cut mother box
  G4SubtractionSolid* bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxL", bbmotherBox, bbleftcutTrap, leftcutrot2, G4ThreeVector(-10*eps, 0.0, -motherdepth/2.0+clear));
  //bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR", bbmothercutBox, bbrightcutTrap, rightcutrot2, G4ThreeVector(10*eps, 0.0, -motherdepth/2.0+clear));
  //EFuchey: 2017/04/17: commented this line to avoid to have some of the GRINCH PMTs outside of the mother volume.
  //Besides, it is not necessary, as it removes some volume from the mother box which would not interfere with anything anyhow
  
  G4Box *frontboxcut = new G4Box("frontboxcut",(bbmagwidth-gapsize)/4.0-coilwidth, 250*cm + eps, clear/2);
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR_fL", bbmothercutBox, frontboxcut, 0, G4ThreeVector( ((bbmagwidth+gapsize)/4.0+coilwidth)+ 10*eps +5*cm, 0, -motherdepth/2.0+clear/2 - 10*eps ));
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR_fLR", bbmothercutBox, frontboxcut, 0, G4ThreeVector( -((bbmagwidth+gapsize)/4.0+coilwidth)-10*eps -5*cm, 0, -motherdepth/2.0+clear/2 - 10*eps ));

  G4Box *bottomboxcut = new G4Box("bottomboxcut",bbmagwidth+eps, 250*cm/2, motherdepth+eps);
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR_fLR_floor", bbmothercutBox, bottomboxcut, 0, G4ThreeVector( 0.0, -1.25*m*2-20*cm, 0.0));
  
  // EFuchey: 2017/05/15: Added this cut to leave some room for the lead shielding for GMn.
  // I checked it was not interfering with the BB magent volume nor with the BB field volume...
  G4Box *beamshieldcut = new G4Box("beamshieldcut", 30.0*cm/2.0 + eps, 300*cm + eps, 10.0*cm/2.0 + eps);
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR_fLR_shieldcut", bbmothercutBox, beamshieldcut, 0, G4ThreeVector( -(bbmagwidth+gapsize)/4.0+coilwidth+15.0*cm, 0.0, -motherdepth/2.0+2.5*cm));
					  
  

  //   Make logical volume for mother vol and place
  G4LogicalVolume *bbmotherLog=new G4LogicalVolume(bbmothercutBox,GetMaterial("Air"),
						   "bbmotherLog", 0, 0, 0);
  
  new G4PVPlacement(bbrm, G4ThreeVector(motherr*sin(fBBang), 0.0, motherr*cos(fBBang)),
		    bbmotherLog, "bbmotherPhys", worldlog, 0,false,0);


  new G4PVPlacement(yokerm,G4ThreeVector(0.0, 0.0, -motherdepth/2.0+clear),
		    bbyokewgapLog, "bbyokewgapPhysical", bbmotherLog, 0,false,0);

  //  Bigbite field log volume
  G4LogicalVolume *bbfieldLog=new G4LogicalVolume(bbairTrap, GetMaterial("Air"),
						  "bbfieldLog", 0, 0, 0);


  fbbfield = new G4SBSBigBiteField( G4ThreeVector(0.0, 0.0, fBBdist), *bbrm );

  fbbfield->fInverted = fDetCon->fGlobalField->fInverted;
  fbbfield->fScaleFactor = fDetCon->GetFieldScale_BB();
  

  G4FieldManager *bbfm = new G4FieldManager(fbbfield);
  //new G4ChordFinder(fbbfield);
  bbfm->SetDetectorField(fbbfield);
  bbfm->CreateChordFinder(fbbfield);
  
  if( fUseLocalField ){
    bbmotherLog->SetFieldManager(bbfm,true);
  }

  //does this volume serve any purpose? Apparently not
  new G4PVPlacement(0, G4ThreeVector(), bbfieldLog, "bbfieldPhysical", bbyokewgapLog, 0,false,0);

  //--------- BigBite Detector Volumes ------------
  //
  // Mother volume is at 10 degrees with BigBite

  double detboxang    = 10.0*deg;
  double detboxheight = 2.5*m;
  double detboxdepth  = 4.0*m;
  double detboxplace  = 0.8*m; // From midplane pivot

  G4RotationMatrix *bbdetrot = new G4RotationMatrix();
  // Y is "down"
  //bbdetrot->rotateZ(180.0*deg);
  //bbdetrot->rotateX(-detboxang);

  bbdetrot->rotateX( detboxang );

  // 325mm is about half the depth of the field volume
  // this is in mother box coords
  double midplanez    = -motherdepth/2.0+clear+325.0*mm;

  G4Box *bbdetbox = new G4Box("bbdetbox", bbmagwidth/2.0, detboxheight/2.0, detboxdepth/2.0);
  G4LogicalVolume *bbdetLog=new G4LogicalVolume(bbdetbox, GetMaterial("Air"),
						"bbdetLog", 0, 0, 0);
  new G4PVPlacement(bbdetrot, G4ThreeVector(0.0, 
					    (detboxplace+detboxdepth/2.0)*sin(detboxang),
					    (detboxplace+detboxdepth/2.0)*cos(detboxang)+midplanez),
		    bbdetLog, "bbdetPhysical", bbmotherLog, 0,false,0);

  //  Just interested in the GEMs for now:

  double detoffset = 0.05*m -detboxdepth/2.0; //z offset of GEM plane positions within BB detector volume: "global" GEM plane z = detoffset + gemz[plane]

  int i;
  int ngem = 0;
  double gemdsep;
  // double *gemz;
  // double *gemw;
  // double *gemh;
  vector<double> gemz, gemw, gemh;

  switch( fGEMOption ){
  case 1:
    ngem = 6;
    gemdsep = 0.05*m;
    break;
  case 2:
    ngem = 5;
    gemdsep = 0.15*m;
    break;
  case 3:
    ngem = 3;
    gemdsep = 0.35*m;
    break;
  default:
    ngem = 4;
    gemdsep = 0.05*m;
    break;
  }

  gemz.resize(ngem);
  gemw.resize(ngem);
  gemh.resize(ngem);
  //
  // GEM option 1
  double gemz_opt1[] = { 0.0*cm, gemdsep, 2.*gemdsep, 3.*gemdsep, fGEMDist, fGEMDist+gemdsep};
  double gemw_opt1[] = { 40.0*cm, 40.0*cm, 40.0*cm, 40.0*cm, 60.0*cm, 60.0*cm };
  double gemh_opt1[] = { 150.0*cm, 150.0*cm, 150.0*cm, 150.0*cm, 200.0*cm, 200.0*cm };

  // GEM option 2
  double gemz_opt2[] = { 0.0*cm, gemdsep, 2.0*gemdsep, 3.0*gemdsep, fGEMDist};
  double gemw_opt2[] = { 40.0*cm, 40.0*cm, 40.0*cm, 40.0*cm, 60.0*cm };
  double gemh_opt2[] = { 150.0*cm, 150.0*cm, 150.0*cm, 150.0*cm, 200.0*cm };

  // GEM option 3
  double gemz_opt3[] = { 0.0*cm, gemdsep, gemdsep*2.0};
  double gemw_opt3[] = { 40.0*cm, 50.0*cm, 50.0*cm};
  double gemh_opt3[] = { 150.0*cm, 200.0*cm, 200.0*cm};

  for( i = 0; i < ngem; i++ ){
    if( fGEMOption == 1 ){
      for( i = 0; i < ngem; i++ ){
	gemz[i] = gemz_opt1[i];
	gemw[i] = gemw_opt1[i];
	gemh[i] = gemh_opt1[i];
      }
    }
    if( fGEMOption == 2 ){
      for( i = 0; i < ngem; i++ ){
	gemz[i] = gemz_opt2[i];
	gemw[i] = gemw_opt2[i];
	gemh[i] = gemh_opt2[i];
      }
    }
    if( fGEMOption == 3 ){
      for( i = 0; i < ngem; i++ ){
	gemz[i] = gemz_opt3[i];
	gemw[i] = gemw_opt3[i];
	gemh[i] = gemh_opt3[i];
      }
    }
  }

  G4RotationMatrix *rot_identity = new G4RotationMatrix;

  G4SBSTrackerBuilder trackerbuilder(fDetCon);

    

  //This routine creates and positions GEM planes in bbdetLog:

  //------------------------------------------- BigBite GEMs: ----------------------------------------//

  //  (fDetCon->TrackerArm)[fDetCon->TrackerIDnumber] = kEarm;
    
  trackerbuilder.BuildComponent(bbdetLog, rot_identity, G4ThreeVector( 0.0, 0.0, detoffset ), ngem, gemz, gemw, gemh, "Earm/BBGEM" );
  //----- Note: Lines of code that are common to the construction of all individual GEM planes/modules were moved to MakeTracker() -----// 
  //----- All we do here in MakeBigBite() is define the number of planes, their z positions, and their transverse dimensions ------//

  double mylarthickness = 0.0020*cm, airthickness = 0.0040*cm;
  double mylar_air_sum = mylarthickness + airthickness;
  double bbpmtz = 0.20*cm;
  
  // **** BIGBITE CALORIMETER MOTHER VOLUME ****:
  G4double bbcal_box_height = 27*8.5*cm;
  G4double bbcal_box_width  = 2.0*37.0*cm;
  //G4double bbcal_box_depth  = (8.5+2.5+37.0)*cm;
  G4double bbcal_box_depth  = (8.5+8.89+37.0)*cm;//8.89 cm (3.5") is the size of the gap between the PS and the SH
  
  // Big Bite Calorimeter shielding.
  // 
  // EFuchey: 2017/03/02: flag for BBECal shielding option: 
  // 0: nothing; 1: default 1/4 in SS+0.5mm; 2: 10cm Al + 3cm SS on the side; 3: 10cm Al + 3cm SS on the side; 
  
  // Default front plate: 0.25" steel + 0.5mm mu metal
  G4double bbcal_shield_thick = 6.85*mm + 9.525*cm;
  G4double Al_thick = 10.0*cm;
  G4double SS_thick = 2.0*cm;
  if(fShieldOption==2)bbcal_shield_thick+=max(0.0, Al_thick-9.0*cm);
  if(fShieldOption==4){
    Al_thick = Al_thick/2.0;
    SS_thick = SS_thick/2.0;
  }
  
  G4Box *bbcalshieldbox = new G4Box( "bbcalshieldbox", bbmagwidth/2.0-2.0*cm, bbcal_box_height/2.0, bbcal_shield_thick/2.0 );
  G4LogicalVolume *bbcal_shield_log = new G4LogicalVolume(bbcalshieldbox, GetMaterial("Air"), "bbcal_shield_log");
  bbcal_shield_log->SetVisAttributes( G4VisAttributes::Invisible );
  
  G4Box *bbcalfrontmufoil = new G4Box( "bbcalfrontmufoil", bbcal_box_width/2.0, bbcal_box_height/2.0, 0.5*mm/2.0 );
  G4LogicalVolume *bbcal_front_mufoil_log = new G4LogicalVolume(bbcalfrontmufoil, GetMaterial("mu-metal"), "bbcal_front_mufoil_log");
  //bbcal_front_mufoil_log->SetVisAttributes( G4Colour(0.,1.0, 0.0) );
  
  G4Box *bbcalfrontsteelplate = new G4Box( "bbcalfrontsteelplate", bbcal_box_width/2.0, bbcal_box_height/2.0, 6.35*mm/2.0 );
  G4LogicalVolume *bbcal_front_steelplate_log = new G4LogicalVolume(bbcalfrontsteelplate, GetMaterial("Steel"), "bbcal_front_steelplate_log");
  bbcal_front_steelplate_log->SetVisAttributes( G4Colour(0.0, 1.0, 1.0) );
  
  // Additional shielding:
  // Attempt 1: option 2: Al 
  G4Box *bbcalshield_al = new G4Box( "bbcalshield_al", bbcal_box_width/2.0+5.0*cm, bbcal_box_height/2.0, Al_thick/2.0 );
  G4LogicalVolume *bbcal_shield_al_log = new G4LogicalVolume(bbcalshield_al, GetMaterial("Aluminum"), "bbcal_shield_al_log");
  bbcal_shield_al_log->SetVisAttributes( G4Colour(0.8, 0.8, 0.8) );
  
  // Attempt 2: option 3: steel 
  G4Box *bbcalshield_ss = new G4Box( "bbcalshield_ss", bbcal_box_width/2.0+5.0*cm, bbcal_box_height/2.0, SS_thick/2.0 );
  G4LogicalVolume *bbcal_shield_ss_log = new G4LogicalVolume(bbcalshield_ss, GetMaterial("Steel"), "bbcal_shield_ss_log");
  bbcal_shield_ss_log->SetVisAttributes( G4Colour(0.4, 0.4, 0.4) );
  
  //shielding on the side option > 2
  G4Box *bbcalshield_side_ss = new G4Box( "bbcalshield_side_ss", 3.0*cm/2.0, detboxheight/2.0, detboxdepth/3.0 );
  G4LogicalVolume *bbcal_shield_side_ss_log = new G4LogicalVolume(bbcalshield_side_ss, GetMaterial("Steel"), "bbcal_shield_side_ss_log");
  bbcal_shield_side_ss_log->SetVisAttributes( G4Colour(1.0, 1.0, 0.9) );
  
  if(fShieldOption){
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, detoffset + fBBCaldist + bbcal_shield_thick/2.0 ), bbcal_shield_log, "bbcal_shield_phys", bbdetLog, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, +bbcal_shield_thick/2.0-0.5*mm/2.0 ), bbcal_front_mufoil_log, "bbcal_front_mufoil_phys", bbcal_shield_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, +bbcal_shield_thick/2.0-0.5*mm-6.35*mm/2.0 ), bbcal_front_steelplate_log, "bbcal_front_steelplate_phys", bbcal_shield_log, false, 0 );
    
    switch(fShieldOption){
    case(2):
      //new G4PVPlacement( 0, G4ThreeVector( 0, 0, -bbcal_shield_thick/2.0+Al_thick/2.0), bbcal_shield_al_log, "bbcal_shield_al_phys", bbcal_shield_log, false, 0 );
      new G4PVPlacement( 0, G4ThreeVector( 0, 0, +bbcal_shield_thick/2.0-0.5*mm-6.35*mm-0.525*cm-Al_thick/2.0), bbcal_shield_al_log, "bbcal_shield_al_phys", bbcal_shield_log, false, 0 );
      //new G4PVPlacement( 0, G4ThreeVector( (-bbmagwidth+3.0*cm)/2.0, 0, -detboxdepth/4.0), bbcal_shield_side_ss_log, "bbcal_shield_side_ss_phys", bbdetLog, false, 0 );
      break;
    case(3):
      new G4PVPlacement( 0, G4ThreeVector( 0, 0, +bbcal_shield_thick/2.0-0.5*mm-6.35*mm-0.525*cm-SS_thick/2.0), bbcal_shield_ss_log, "bbcal_shield_ss_phys", bbcal_shield_log, false, 0 );
      //new G4PVPlacement( 0, G4ThreeVector( (-bbmagwidth+3.0*cm)/2.0, 0, -detboxdepth/4.0), bbcal_shield_side_ss_log, "bbcal_shield_side_ss_phys", bbdetLog, false, 0 );
      break;
    case(4):
      new G4PVPlacement( 0, G4ThreeVector( 0, 0, +bbcal_shield_thick/2.0-0.5*mm-6.35*mm-0.525*cm-SS_thick/2.0), bbcal_shield_ss_log, "bbcal_shield_ss_phys", bbcal_shield_log, false, 0 );
      new G4PVPlacement( 0, G4ThreeVector( 0, 0, +bbcal_shield_thick/2.0-0.5*mm-6.35*mm-0.525*cm-SS_thick-Al_thick/2.0), bbcal_shield_al_log, "bbcal_shield_al_phys", bbcal_shield_log, false, 0 );
      //new G4PVPlacement( 0, G4ThreeVector( (-bbmagwidth+3.0*cm)/2.0, 0, -detboxdepth/4.0), bbcal_shield_side_ss_log, "bbcal_shield_side_ss_phys", bbdetLog, false, 0 );
    default:
      break;
    }
  }
  
  // BB Ecal
  G4Box *bbcalbox = new G4Box( "bbcalbox", bbcal_box_width/2.0, bbcal_box_height/2.0, bbcal_box_depth/2.0+mm );
  G4LogicalVolume *bbcal_mother_log = new G4LogicalVolume(bbcalbox, GetMaterial("Air"), "bbcal_mother_log");
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, detoffset + fBBCaldist + bbcal_shield_thick + bbcal_box_depth/2.0 ), bbcal_mother_log, "bbcal_mother_phys", bbdetLog, false, 0 ); 

  bbcal_mother_log->SetVisAttributes( G4VisAttributes::Invisible );
  
  //option to "turn off" BBCAL (make total absorber)
  if( (fDetCon->StepLimiterList).find( "bbcal_mother_log" ) != (fDetCon->StepLimiterList).end() ){
    bbcal_mother_log->SetUserLimits( new G4UserLimits(0,0,0,DBL_MAX,DBL_MAX) );
    G4String SDname = "Earm/BBCal";
    G4String collname = "BBCalHitsCollection";
    G4SBSCalSD *BBCalSD = NULL;
    if( !( (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(SDname)) ){ //add sensitivity:
      G4cout << "Adding BBCal sensitive detector to SDman..." << G4endl;
      BBCalSD = new G4SBSCalSD( SDname, collname );
      fDetCon->fSDman->AddNewDetector( BBCalSD );
      (fDetCon->SDlist).insert( SDname );
      fDetCon->SDtype[SDname] = kCAL;
      (BBCalSD->detmap).depth = 0;
      bbcal_mother_log->SetSensitiveDetector( BBCalSD );
      //(BBCalSD->detmap).Row[0] = 0;
      //(BBCalSD->detmap).Col[0] = 0;
      
    }
  }

  // **** BIGBITE PRESHOWER **** 
  // 2 columns, 27 rows

  double psheight = 27*8.5*cm;
  double pswidth  = 2.0*37.0*cm;
  double psdepth  = 8.5*cm;

  G4Box *bbpsbox = new G4Box("bbpsbox", pswidth/2.0, psheight/2.0, psdepth/2.0 );
  G4LogicalVolume *bbpslog = new G4LogicalVolume(bbpsbox, GetMaterial("Air"), "bbpslog");
  //new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth/2.0), bbpslog, "bbpsphys", bbdetLog, false, 0);
  new G4PVPlacement(0, G4ThreeVector( 0, 0, -bbcal_box_depth/2.0 + psdepth/2.0 ), bbpslog, "bbpsphys", bbcal_mother_log, false, 0 );
  
  //placement of second mu-metal foil behind the PS
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, -bbcal_box_depth/2.0 + psdepth + 0.5*mm/2.0), bbcal_front_mufoil_log, "bbcal_back_mufoil_phys", bbcal_mother_log, false, 0 );
  // Preshower module - geometry will be assigned after Shower

  // **** BIGBITE HODOSCOPE **** 
  // Scintillator box - same dimensions as preshower
  double bbhododepth = 8.8*cm;// logic volume...
  double bbslat_length = 60.0*cm;
  double bbslat_section = 2.5*cm;
  G4Box *bbhodobox = new G4Box("bbhodobox", pswidth/2.0, psheight/2.0, bbhododepth/2.0 );
  G4LogicalVolume *bbhodolog = new G4LogicalVolume( bbhodobox, GetMaterial("Air"), "bbhodolog" );
  //new G4PVPlacement(0, G4ThreeVector(0.0,0.0, detoffset+fBBCaldist+psdepth+bbhododepth/2.0), bbhodolog, "bbhodophys", bbdetLog, false, 0);
  //new G4PVPlacement( 0, G4ThreeVector(0,0, -bbcal_box_depth/2.0 + psdepth + bbhododepth/2.0 ), bbhodolog, "bbhodophys", bbcal_mother_log, false, 0 );
  //new G4PVPlacement( 0, G4ThreeVector(0,0, -bbcal_box_depth/2.0 + psdepth + 0.217*2.54 + bbhododepth/2.0 ), bbhodolog, "bbhodophys", bbcal_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0, -bbcal_box_depth/2.0 + psdepth + bbhododepth/2.0+0.5*mm ), bbhodolog, "bbhodophys", bbcal_mother_log, false, 0 );
  bbhodolog->SetVisAttributes(G4VisAttributes::Invisible);
  
  //
  G4Box *bbhodoslatbox = new G4Box("bbhodoslat", bbslat_length/2.0, bbslat_section/2.0, bbslat_section/2.0);
  G4LogicalVolume *bbhodoslatlog = new G4LogicalVolume( bbhodoslatbox, GetMaterial("BBHodo_Scinti"), "bbhodoslatlog" );
  bbhodoslatlog->SetVisAttributes(G4Colour(0.0, 1.0, 0.0));
  
  G4SDManager *sdman = fDetCon->fSDman;

  G4String BBHodoScintSDname = "Earm/BBHodoScint";
  G4String BBHodoScintcollname = "BBHodoScintHitsCollection";
  G4SBSCalSD *BBHodoScintSD = NULL;
  
  if( !((G4SBSCalSD*) sdman->FindSensitiveDetector(BBHodoScintSDname)) ) {
    G4cout << "Adding BB Hodoscope Scint Sensitive Detector to SDman..." << G4endl;
    BBHodoScintSD = new G4SBSCalSD( BBHodoScintSDname, BBHodoScintcollname );
    
    //    BBHodoScintSD->SetEnergyThreshold( 5.0*MeV ); //Based on simulations of quasi-elastic scattering
    //BBHodoScintSD->SetNTimeBins( 60 ); //0.5 ns/bin
    sdman->AddNewDetector( BBHodoScintSD );
    (fDetCon->SDlist).insert( BBHodoScintSDname );
    fDetCon->SDtype[BBHodoScintSDname] = kCAL;
    (BBHodoScintSD->detmap).depth = 0;

    G4double ethresh_default = 5.0*MeV;
    G4double timewindow_default = 30.0*ns;
    
    fDetCon->SetTimeWindowAndThreshold( BBHodoScintSDname, ethresh_default, timewindow_default );
  }
  bbhodoslatlog->SetSensitiveDetector( BBHodoScintSD ); 

  ofstream mapfile("database/BBhodo_map.txt");

  TString currentline;

  currentline.Form("# %15s, %15s, %15s, %18s, %18s",
		   "Cell", "Row", "Column", "Xcenter", "Ycenter" );
  
  G4int n_bbhodoslats = 90;
  for(int i_bbhslat = 0; i_bbhslat<n_bbhodoslats; i_bbhslat++){
    G4double y_slat = n_bbhodoslats*bbslat_section/2.0-(G4double(i_bbhslat)+0.5)*bbslat_section;
    G4double z_slat = -bbhododepth/2.0-0.5*mm+0.217*2.54*cm+bbslat_section/2.0;
    new G4PVPlacement( 0, G4ThreeVector(0, y_slat, z_slat), bbhodoslatlog, "bbhodoslatphys", bbhodolog, false, i_bbhslat );
    (BBHodoScintSD->detmap).Col[i_bbhslat] = 0;
    (BBHodoScintSD->detmap).Row[i_bbhslat] = i_bbhslat;
    (BBHodoScintSD->detmap).LocalCoord[i_bbhslat] = G4ThreeVector( 0, y_slat,  z_slat);

    currentline.Form( "  %15d, %15d, %15d, %18.3f, %18.3f",
		      i_bbhslat, i_bbhslat, 0, 0.0/cm, y_slat/cm );

    mapfile << currentline << endl;
  }

  mapfile.close();
  //0.217" is the gap between the PS and the hodoscope

  // **** BIGBITE SHOWER ****
  // 7 columns, 27 rows
  double calheight = 27*8.5*cm;
  double calwidth  = 7*8.5*cm;
  double caldepth  = 37.0*cm;
  G4Box *bbshowerbox = new G4Box("bbshowerbox", calwidth/2.0, calheight/2.0, caldepth/2.0);
  G4LogicalVolume *bbshowerlog = new G4LogicalVolume(bbshowerbox, GetMaterial("Air"), "bbshowerlog");
  //new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth+bbhododepth+caldepth/2.0), bbshowerlog, "bbshowerphys", bbdetLog, false, 0);
  //new G4PVPlacement( 0, G4ThreeVector( 0, 0, -bbcal_box_depth/2.0 + psdepth + bbhododepth + caldepth/2.0), bbshowerlog, "bbshowerphys", bbcal_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, +bbcal_box_depth/2.0 - caldepth/2.0), bbshowerlog, "bbshowerphys", bbcal_mother_log, false, 0 );
  

  // Shower module:
  double bbmodule_x = 8.5*cm, bbmodule_y = 8.5*cm;  
  double bbTF1_x = bbmodule_x - 2*mylar_air_sum;
  double bbTF1_y = bbmodule_y - 2*mylar_air_sum;
  double bbTF1_z = caldepth - 2*bbpmtz - mylar_air_sum;

  G4Box *showermodbox = new G4Box("showermodbox", bbmodule_x/2.0, bbmodule_y/2.0, caldepth/2.0);
  G4LogicalVolume *showermodlog = new G4LogicalVolume(showermodbox, GetMaterial("Special_Air"), "showermodlog");

  G4Box *tempbox = new G4Box("tempbox", bbmodule_x/2.0, bbmodule_y/2.0, (caldepth-2*bbpmtz)/2.0);

  // Subtraction
  G4Box *showermodbox_sub = new G4Box( "showermodbox_sub", (bbmodule_x-2*mylarthickness)/2.0, (bbmodule_y-2*mylarthickness)/2.0, (caldepth-2*bbpmtz)/2.0 );
  G4SubtractionSolid *bbmylarwrap = new G4SubtractionSolid( "bbmylarwrap", tempbox, showermodbox_sub, 0, G4ThreeVector(0.0, 0.0, mylarthickness) );
  G4LogicalVolume *bbmylarwraplog = new G4LogicalVolume( bbmylarwrap, GetMaterial("Mylar"), "bbmylarwraplog" ); 
  new G4LogicalSkinSurface( "BB Mylar Skin", bbmylarwraplog, GetOpticalSurface("Mirrsurf") );

  // Make Lead Glass 
  G4Box *bbTF1box = new G4Box( "bbTF1box", bbTF1_x/2.0, bbTF1_y/2.0, bbTF1_z/2.0 );
  G4LogicalVolume *bbTF1log = new G4LogicalVolume( bbTF1box, GetMaterial("TF5"), "bbTF1log" );
  
  // Shower TF1 SD of type CAL
  //G4SDManager *sdman = fDetCon->fSDman;

  G4String BBSHTF1SDname = "Earm/BBSHTF1";
  G4String BBSHTF1collname = "BBSHTF1HitsCollection";
  G4SBSCalSD *BBSHTF1SD = NULL;

  if( !((G4SBSCalSD*) sdman->FindSensitiveDetector(BBSHTF1SDname)) ) {
    G4cout << "Adding BB Shower TF1 Sensitive Detector to SDman..." << G4endl;
    BBSHTF1SD = new G4SBSCalSD( BBSHTF1SDname, BBSHTF1collname );
      
    sdman->AddNewDetector( BBSHTF1SD );
    (fDetCon->SDlist).insert( BBSHTF1SDname );
    fDetCon->SDtype[BBSHTF1SDname] = kCAL;
    (BBSHTF1SD->detmap).depth = 1;

    G4double threshold_default = 10.0*MeV; //1% of 1 GeV
    G4double timewindow_default = 50.0*ns; //We could use 10 ns here if we wanted, but also have to consider pulse shape. 
    
    fDetCon->SetTimeWindowAndThreshold( BBSHTF1SDname, threshold_default, timewindow_default );
  }
  bbTF1log->SetSensitiveDetector( BBSHTF1SD ); 

  if( (fDetCon->StepLimiterList).find( BBSHTF1SDname ) != (fDetCon->StepLimiterList).end() ){
    bbTF1log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }

  // Make PMT/Window
  double pmtrad = 7.62/2.0*cm;
  G4Tubs *bbPMT = new G4Tubs( "bbPMT", 0.0*cm, pmtrad, bbpmtz/2.0, 0.0, twopi );
  G4LogicalVolume *bbpmtwindowlog = new G4LogicalVolume( bbPMT, GetMaterial("QuartzWindow_ECal"), "bbpmtwindowlog" );
  G4LogicalVolume *bbpmtcathodelog = new G4LogicalVolume( bbPMT, GetMaterial("Photocathode_BB"), "bbpmtcathodelog" );

  // Shower PMT SD of type ECAL
  G4String BBSHSDname = "Earm/BBSH";
  G4String BBSHcollname = "BBSHHitsCollection";
  G4SBSECalSD *BBSHSD = NULL;

  if( !((G4SBSECalSD*) sdman->FindSensitiveDetector(BBSHSDname)) ) {
    G4cout << "Adding BB Shower PMT Sensitive Detector to SDman..." << G4endl;
    BBSHSD = new G4SBSECalSD( BBSHSDname, BBSHcollname );
    sdman->AddNewDetector( BBSHSD );
    (fDetCon->SDlist).insert(BBSHSDname);
    fDetCon->SDtype[BBSHSDname] = kECAL;
    (BBSHSD->detmap).depth = 1;
  }
  bbpmtcathodelog->SetSensitiveDetector( BBSHSD );

  // Put everything in a BB Shower Module
  int shower_copy_number = 0;
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-bbpmtz)/2.0), bbpmtcathodelog,"bbcathodephys", showermodlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-3*bbpmtz)/2.0), bbpmtwindowlog, "bbwindowphys", showermodlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-4*bbpmtz-bbTF1_z)/2.0), bbTF1log, "bbTF1phys", showermodlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, -bbpmtz), bbmylarwraplog, "bbmylarphys", showermodlog, false, 0 );

  mapfile.open("database/BBSH_blockmap.txt");

  currentline.Form("# %15s, %15s, %15s, %18s, %18s",
		   "Cell", "Row", "Column", "Xcenter", "Ycenter" );

  mapfile << currentline << endl;
  
  int bbscol = 7;
  int bbsrow = 27;
  for( int l=0; l<bbscol; l++ ) {
    for( int j=0; j<bbsrow; j++ ) {
      (BBSHSD->detmap).Col[shower_copy_number] = l;
      (BBSHSD->detmap).Row[shower_copy_number] = j;
      (BBSHTF1SD->detmap).Col[shower_copy_number] = l;
      (BBSHTF1SD->detmap).Row[shower_copy_number] = j;
      double xtemp = (calwidth - bbmodule_x)/2.0 - l*bbmodule_x;
      double ytemp = (calheight - bbmodule_y)/2.0 - j*bbmodule_y;

      currentline.Form( "  %15d, %15d, %15d, %18.3f, %18.3f",
			shower_copy_number, j, l, xtemp/cm, ytemp/cm );

      mapfile << currentline << endl;
      
      new G4PVPlacement(0, G4ThreeVector(xtemp,ytemp,0.0), showermodlog, "showermodphys", bbshowerlog, false, shower_copy_number);
      
      (BBSHSD->detmap).LocalCoord[shower_copy_number] = G4ThreeVector( xtemp,ytemp,(caldepth-bbpmtz)/2.0  );
      (BBSHTF1SD->detmap).LocalCoord[shower_copy_number] = G4ThreeVector( xtemp, ytemp, (caldepth-4*bbpmtz-bbTF1_z)/2.0 );
      shower_copy_number++;
    }
  }

  mapfile.close();

  mapfile.open("database/BBPS_blockmap.txt");

  currentline.Form("# %15s, %15s, %15s, %18s, %18s",
		   "Cell", "Row", "Column", "Xcenter", "Ycenter" );

  mapfile << currentline << endl;

  // ****Preshower Continued****
  // Reusing modules from Shower (same variables), rotated by either +/- 90 deg depending on column #
  G4Box *preshowermodbox = new G4Box( "preshowermodbox", bbmodule_x/2.0, bbmodule_y/2.0, caldepth/2.0 );
  G4LogicalVolume *preshowermodlog = new G4LogicalVolume( preshowermodbox, GetMaterial("Special_Air"), "preshowermodlog" );
 
  // Preshower TF1 SD of type CAL
  G4LogicalVolume *bbpsTF1log = new G4LogicalVolume( bbTF1box, GetMaterial("TF5"), "bbpsTF1log" );

  G4String BBPSTF1SDname = "Earm/BBPSTF1";
  G4String BBPSTF1collname = "BBPSTF1HitsCollection";
  G4SBSCalSD *BBPSTF1SD = NULL;

  if( !((G4SBSCalSD*) sdman->FindSensitiveDetector(BBPSTF1SDname)) ) {
    G4cout << "Adding BB Preshower TF1 Sensitive Detector to SDman..." << G4endl;
    BBPSTF1SD = new G4SBSCalSD( BBPSTF1SDname, BBPSTF1collname );
    
    sdman->AddNewDetector( BBPSTF1SD );
    (fDetCon->SDlist).insert( BBPSTF1SDname );
    fDetCon->SDtype[BBPSTF1SDname] = kCAL;
    (BBPSTF1SD->detmap).depth = 1;

    //Photoelectron yield is approximately 500/GeV (or so)
    G4double threshold_default = 10.0*MeV; //1% of 1 GeV
    G4double timewindow_default = 50.0*ns; //We could use 10 ns here if we wanted, but also have to consider pulse shape.

    fDetCon->SetTimeWindowAndThreshold( BBPSTF1SDname, threshold_default, timewindow_default );
  }
  bbpsTF1log->SetSensitiveDetector( BBPSTF1SD ); 

  if( (fDetCon->StepLimiterList).find( BBPSTF1SDname ) != (fDetCon->StepLimiterList).end() ){
    bbpsTF1log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }

  // Preshower PMT SD of type ECAL
  G4LogicalVolume *bbpspmtcathodelog = new G4LogicalVolume( bbPMT, GetMaterial("Photocathode_material_ecal"), "bbpspmtcathodelog" );

  G4String BBPSSDname = "Earm/BBPS";
  G4String BBPScollname = "BBPSHitsCollection";
  G4SBSECalSD *BBPSSD = NULL;

  if( !((G4SBSECalSD*) sdman->FindSensitiveDetector(BBPSSDname)) ) {
    G4cout << "Adding BB Preshower PMT Sensitive Detector to SDman..." << G4endl;
    BBPSSD = new G4SBSECalSD( BBPSSDname, BBPScollname );
    sdman->AddNewDetector( BBPSSD );
    (fDetCon->SDlist).insert(BBPSSDname);
    fDetCon->SDtype[BBPSSDname] = kECAL;
    (BBPSSD->detmap).depth = 1;
  }
  bbpspmtcathodelog->SetSensitiveDetector( BBPSSD );

  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-bbpmtz)/2.0), bbpspmtcathodelog,"bbpscathodephys", preshowermodlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-3*bbpmtz)/2.0), bbpmtwindowlog, "bbpswindowphys", preshowermodlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-4*bbpmtz-bbTF1_z)/2.0), bbpsTF1log, "bbpsTF1phys", preshowermodlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, -bbpmtz), bbmylarwraplog, "bbpsmylarphys", preshowermodlog, false, 0 );
  
  G4RotationMatrix *bbpsrm_col1 = new G4RotationMatrix;
  bbpsrm_col1->rotateY(-90.0*deg);
  G4RotationMatrix *bbpsrm_col2 = new G4RotationMatrix;
  bbpsrm_col2->rotateY(90.0*deg);
  
  int bbpscol = 2;
  int bbpsrow = 27;
  int ps_copy_number = 0;
  for(int l=0; l<bbpscol; l++) {
    for(int j=0; j<bbpsrow; j++) {
      double xtemp = (pswidth-caldepth)/2.0 - l*caldepth;
      double ytemp = (psheight-bbmodule_y)/2.0 - j*bbmodule_y;
      (BBPSSD->detmap).Col[ps_copy_number] = l;
      (BBPSSD->detmap).Row[ps_copy_number] = j;
      (BBPSTF1SD->detmap).Col[ps_copy_number] = l;
      (BBPSTF1SD->detmap).Row[ps_copy_number] = j;

      currentline.Form( "  %15d, %15d, %15d, %18.3f, %18.3f",
			  ps_copy_number, j, l, xtemp/cm, ytemp/cm );

      mapfile << currentline << endl;
      
      if(l==0) { 
	new G4PVPlacement( bbpsrm_col1, G4ThreeVector(xtemp,ytemp,0.0), preshowermodlog, "preshowermodphys", bbpslog, false, ps_copy_number );
	(BBPSSD->detmap).LocalCoord[ps_copy_number] = G4ThreeVector(xtemp+caldepth/2.0-bbpmtz/2.0, ytemp, 0.0);
	(BBPSTF1SD->detmap).LocalCoord[ps_copy_number] = G4ThreeVector(xtemp,ytemp,0.0);
	ps_copy_number++;
      }
      if(l==1) {
	new G4PVPlacement( bbpsrm_col2, G4ThreeVector(xtemp,ytemp,0.0), preshowermodlog, "preshowermodphys", bbpslog, false, ps_copy_number );
	(BBPSSD->detmap).LocalCoord[ps_copy_number] = G4ThreeVector(xtemp-caldepth/2.0+bbpmtz/2.0, ytemp, 0.0);
	(BBPSTF1SD->detmap).LocalCoord[ps_copy_number] = G4ThreeVector(xtemp,ytemp,0.0);
	ps_copy_number++;
      }
    }
  }

  mapfile.close();
  
  //--------- Visualization attributes -------------------------------
  //Mother volumes
  bbdetLog->SetVisAttributes( G4VisAttributes::Invisible );
  bbfieldLog->SetVisAttributes( G4VisAttributes::Invisible );
  bbmotherLog->SetVisAttributes( G4VisAttributes::Invisible );
  
  bbpslog->SetVisAttributes( G4VisAttributes::Invisible ); 
  bbshowerlog->SetVisAttributes( G4VisAttributes::Invisible );

  //Mylar
  G4VisAttributes *mylar_colour = new G4VisAttributes(G4Colour( 0.5, 0.5, 0.5 ) );
  bbmylarwraplog->SetVisAttributes(mylar_colour);

  //Air
  showermodlog->SetVisAttributes( G4VisAttributes::Invisible );
  preshowermodlog->SetVisAttributes( G4VisAttributes::Invisible );

  //TF1
  G4VisAttributes *TF1_colour = new G4VisAttributes(G4Colour( 0.8, 0.8, 0.0 ) );
  bbTF1log->SetVisAttributes(TF1_colour);
  bbpsTF1log->SetVisAttributes(TF1_colour);

  //PMTcathode
  G4VisAttributes *PMT_colour = new G4VisAttributes(G4Colour( G4Colour::Blue() ));
  PMT_colour->SetForceLineSegmentsPerCircle( 12 );
  bbpmtcathodelog->SetVisAttributes(PMT_colour);
  bbpspmtcathodelog->SetVisAttributes(PMT_colour);

  //Yoke
  G4VisAttributes * yokeVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  bbyokewgapLog->SetVisAttributes(yokeVisAtt);

//G4LogicalVolume* bbcallog = new G4LogicalVolume(bbcalbox, GetMaterial("Lead"), "bbcallog");
  
  //new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth+caldepth/2.0+5.0*cm), bbcallog,
  //"bbcalphys", bbdetLog, false, 0, false);

  // G4String BBCalSDname = "Earm/BBCal";
  // G4String BBCalcolname = "BBCalHitsCollection";
  // G4SBSCalSD* BBCalSD;

  // if( !(BBCalSD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(BBCalSDname)) ){
  //   BBCalSD = new G4SBSCalSD( BBCalSDname, BBCalcolname );
  //   fDetCon->fSDman->AddNewDetector(BBCalSD);
  //   (fDetCon->SDlist).insert( BBCalSDname );
  //   fDetCon->SDtype[BBCalSDname] = kCAL;
  //   //fDetCon->SDarm[BBCalSDname] = kEarm;
  // }

  // bbcallog->SetSensitiveDetector(BBCalSD);
  // bbcallog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );

  // (BBCalSD->detmap).depth=0;
  // (BBCalSD->detmap).Row[0] = 0;
  // (BBCalSD->detmap).Col[0] = 0;
  // (BBCalSD->detmap).LocalCoord[0] = G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth+caldepth/2.0+5.0*cm);

  // G4cout << "fDetCon->StepLimiterList.size() == " << (fDetCon->StepLimiterList).size() << G4endl;
  // for( set<G4String>::iterator itlist=(fDetCon->StepLimiterList).begin(); itlist != (fDetCon->StepLimiterList).end(); itlist++){
  //   G4cout << "step limiter list element = " << *itlist << G4endl;
  // }

  //(BBCalSD->detmap).GlobalCoord[0] = G4ThreeVector(0.0,0.0,0.0); //To be set later 

  //--------- BigBite Cerenkov ------------------------------

  /*
  //  double cer_mirrorthick = 0.635*mm;
  double cer_mirrorthick = 3.00*mm;
  //double cer_winthick_in   = 1.0*mm;
  double cer_winthick_in   = 0.1*mm;
  double cer_winthick_out  = 0.2*mm;

  double cer_width  =  50.0*cm;
  double cer_height = 200.0*cm;

  //  G4Box *cer_winbox = new G4Box("cer_winbox", cer_width/2.0, cer_height/2.0, cer_winthick/2.0 );
  G4Box *cer_winbox_in = new G4Box("cer_winbox_in", cer_width/2.0, cer_height/2.0, cer_winthick_in/2.0 );
  G4Box *cer_winbox_out = new G4Box("cer_winbox_out", cer_width/2.0, cer_height/2.0, cer_winthick_out/2.0 );
  G4Box *cer_mirbox = new G4Box("cer_mirbox", cer_width/2.0, cer_height/2.0, cer_mirrorthick/2.0 );
  G4Box *cer_gasbox = new G4Box("cer_gasbox", cer_width/2.0, cer_height/2.0, fCerDepth/2.0 );

  G4LogicalVolume* cer_winlog_in = new G4LogicalVolume(cer_winbox_in, GetMaterial("Aluminum"), "cer_winlog_in");
  G4LogicalVolume* cer_winlog_out = new G4LogicalVolume(cer_winbox_out, GetMaterial("Aluminum"), "cer_winlog_out");
  //  G4LogicalVolume* cer_mirlog = new G4LogicalVolume(cer_mirbox, SiO2, "cer_mirlog");
  G4LogicalVolume* cer_mirlog = new G4LogicalVolume(cer_mirbox, GetMaterial("Acrylic"), "cer_mirlog");
  G4LogicalVolume* cer_gaslog = new G4LogicalVolume(cer_gasbox, GetMaterial("C4F8O"), "cer_gaslog");

  double thisz = detoffset+fCerDist;

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, thisz + cer_winthick_in/2.0 ), cer_winlog_in, "cerwin1", bbdetLog, false, 0, false);
  thisz += cer_winthick_in;
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, thisz + fCerDepth/2.0 ), cer_gaslog, "cergas", bbdetLog, false, 0, false);
  thisz += fCerDepth;
  // new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, thisz + cer_mirrorthick/2.0 ), cer_mirlog, "cermir", bbdetLog, false, 0, false);
  if( fCerDepth > 20.0*cm ){
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, fCerDepth/2.0-20.0*cm ), cer_mirlog, "cermir", cer_gaslog, false, 0, false);
  } else {
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, fCerDepth/2.0 ), cer_mirlog, "cermir", cer_gaslog, false, 0, false);
  }
  //  thisz += cer_mirrorthick;
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, thisz + cer_winthick_out/2.0 ), cer_winlog_out, "cerwin2", bbdetLog, false, 0, false);
  */

  G4SBSGrinch *grinch = new G4SBSGrinch(fDetCon);
  grinch->SetZOffset( detoffset + fCerDist );
  grinch->SetCerDepth( fCerDepth);
  grinch->BuildComponent(bbdetLog);
  grinch->SetGrinchGas( fGRINCHgas );

}

void G4SBSEArmBuilder::MakeGMnGEMShielding( G4LogicalVolume *motherlog ){
  /////////////////////////////////////////////////////////////////
  //
  // Working with very little information here.. apparently it is 
  // "similar" to the hut that dasuni built but rotated..

  // Got coordinates from Alan which apparently correspond to the center of
  // the hut, I am assuming that the hut is also rotated. Need confirmation on this

  G4SDManager *sdman = fDetCon->fSDman;

  G4double GboxX = 52.0*2.54*cm;
  G4double GboxY = 26.0*2.54*cm;
  G4double GboxZ = 52.0*2.54*cm;
   
  G4double SPlateX = 120.0*2.54*cm;
  G4double SPlateY =   5.0*2.54*cm;
  G4double SPlateZ =  43.0*2.54*cm;

  // Vertical plate
  G4double GPlateX1 = 96.0*2.54*cm;
  G4double GPlateY1 = GboxY - 2.5*2.54*cm;
  G4double GPlateZ1 = 7.5*2.54*cm;

  // Horizontal plate to match the heights
  G4double GPlateX2 = 96.0*2.54*cm;
  G4double GPlateY2 = 2.5*2.54*cm;
  G4double GPlateZ2 = GboxX;

  // Make all the parts:
  G4Box *GreenBox = new G4Box( "GreenBox",GboxX/2.0, GboxY/2.0, GboxZ/2.0);
  G4LogicalVolume *GreenBox_log = new G4LogicalVolume( GreenBox, GetMaterial("Steel"), 
						       "GreenBox_log" );

  G4Box *SteelPlate = new G4Box( "SteelPlate", SPlateX/2.0, SPlateY/2.0, SPlateZ/2.0);
  G4LogicalVolume *SteelPlate_log = new G4LogicalVolume( SteelPlate, GetMaterial("Steel"), 
							 "SteelPlate_log" );

  G4Box *GreenPlate1 = new G4Box( "GreenPlate1", GPlateX1/2.0, GPlateY1/2.0, GPlateZ1/2.0);
  G4LogicalVolume *GreenPlate1_log = new G4LogicalVolume( GreenPlate1, GetMaterial("Steel"), 
							  "GreenPlate1_log" );

  G4Box *GreenPlate2 = new G4Box( "GreenPlate2", GPlateX2/2.0, GPlateY2/2.0, GPlateZ2/2.0);
  G4LogicalVolume *GreenPlate2_log = new G4LogicalVolume( GreenPlate2, GetMaterial("Steel"), 
							  "GreenPlate2_log" );

  // Make a Mother Volume to house everything:
  G4double ShieldMotherX = 2.0*GboxX + GPlateX1;
  G4double ShieldMotherY = GboxY + SPlateY;
  G4double ShieldMotherZ = GboxZ;

  G4Box *ShieldBox = new G4Box( "ShieldBox", ShieldMotherX/2.0, ShieldMotherY/2.0, 
				ShieldMotherZ/2.0 );
  G4LogicalVolume *ShieldLog = new G4LogicalVolume( ShieldBox, GetMaterial("Air"), "ShieldLog");
  G4VisAttributes *temp = new G4VisAttributes(G4Colour(0.0,0.6,0.0));
  //temp->SetForceWireframe(true);
  ShieldLog->SetVisAttributes(G4VisAttributes::Invisible);
 
  // And place everything within the Mother:
  new G4PVPlacement( 0, G4ThreeVector(-ShieldMotherX/2.0 + GboxX/2.0, -SPlateY/2.0, 0.0 ), 
		     GreenBox_log, "LeftBox",  ShieldLog, false, 0 );

  new G4PVPlacement( 0, G4ThreeVector( ShieldMotherX/2.0 - GboxX/2.0, -SPlateY/2.0, 0.0 ), 
		     GreenBox_log, "RightBox", ShieldLog, false, 1 );

  new G4PVPlacement( 0, G4ThreeVector( 0.0, ShieldMotherY/2.0 - SPlateY/2.0, ShieldMotherZ/2.0 - SPlateZ/2.0 ), 
		     SteelPlate_log, "TopPlate", ShieldLog, false, 0 );

  new G4PVPlacement( 0, G4ThreeVector( 0.0, -ShieldMotherY/2.0 + GPlateY1/2.0 + GPlateY2, -ShieldMotherZ/2.0 + GPlateZ1/2.0 ), 
		     GreenPlate1_log, "VerticalPlate", ShieldLog, false, 0 );

  new G4PVPlacement( 0, G4ThreeVector( 0.0, -ShieldMotherY/2.0 + GPlateY2/2.0, 0.0), 
		     GreenPlate2_log, "BottomPlate", ShieldLog, false, 0 );

  // In order to calculate the dose, we need a SD of type CAL:
  G4double ElecX = 150.0*cm;
  G4double ElecY = 40.0*cm;
  G4double ElecZ = 0.5*cm;
 
  G4Box *Electronics = new G4Box( "Electronics" , ElecX/2.0, ElecY/2.0, ElecZ/2.0);
  G4LogicalVolume *Electronics_log = new G4LogicalVolume( Electronics , GetMaterial("Silicon"), "Electronics_log" );
  
  G4String GEMElectronicsname = "Earm/GEMElectronics";
  G4String GEMElectronicscollname = "GEMElectronicsHitsCollection";
  G4SBSCalSD *GEMElecSD = NULL;

  GEMElectronicsname += "GMn";
  GEMElectronicscollname += "GMn";

  if( !( (G4SBSCalSD*) sdman->FindSensitiveDetector(GEMElectronicsname) )){
    G4cout << "Adding GEM electronics Sensitive Detector to SDman..." << G4endl;
    GEMElecSD = new G4SBSCalSD( GEMElectronicsname, GEMElectronicscollname );
    sdman->AddNewDetector(GEMElecSD);
    (fDetCon->SDlist).insert(GEMElectronicsname);
    fDetCon->SDtype[GEMElectronicsname] = kCAL;
    (GEMElecSD->detmap).depth = 1;

    fDetCon->SetTimeWindowAndThreshold( GEMElectronicsname );
  }
  Electronics_log->SetSensitiveDetector( GEMElecSD );
  
  if( (fDetCon->StepLimiterList).find( GEMElectronicsname ) != (fDetCon->StepLimiterList).end() ){
    Electronics_log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }

  // Place the electronics in our hut:
  // new G4PVPlacement( 0, G4ThreeVector(0.0, -ShieldMotherY/2.0 + GPlateY2 + ElecY/2.0, ShieldMotherZ/2.0 - GPlateZ1 - ElecZ/2.0),
  // 		     Electronics_log, "Electronics", ShieldLog, false, 0);
  new G4PVPlacement( 0, G4ThreeVector(0.0, -1.25*cm, 0.0),
		     Electronics_log, "GMn_Electronics", ShieldLog, false, 0);

  // Numbers come from email exchange with Alan Gavalya - he says the coordinate
  // system is such that z points upstream, but did not elaborate on x/y. I made the
  // assumption that y is "up" and x is beam-left

  G4double inch = 2.54*cm;
  double x =  190.0795 * inch;
  double y = -105.6100 * inch; // + ShieldMotherY/2.0;
  double z =  187.0807 * inch;
  G4ThreeVector pos_mom(x,y,z);

  G4ThreeVector pos_temp(x,0,z);
  G4ThreeVector punit = pos_temp.unit();
  double theta = acos(punit.z());

  G4RotationMatrix *hutrm = new G4RotationMatrix;
  hutrm->rotateY(-theta);
  
  // **** Estimation ****
  // Bogdan literally told me to look at a printed engineering document 
  // and find the r / theta / hut dimensions by using a ruler...I wasted 5 minutes
  // of my time and this is the result (pretty close to Alan's #s though)
 
  // double r = 273.7 * inch;
  // double th = 45.8*3.141592/180.0;
  // G4ThreeVector estimate(r*sin(th),y,r*cos(th));
  // G4RotationMatrix *hutrm_est = new G4RotationMatrix;
  // hutrm_est->rotateY(-th);

  //new G4PVPlacement( hutrm_est, estimate, ShieldLog, "ShieldMother", motherlog, false, 0 );

  new G4PVPlacement( hutrm, pos_mom, ShieldLog, "ShieldMother", motherlog, false, 0 );

  // VISUALS:
  G4VisAttributes *BlockAtt = new G4VisAttributes(G4Colour(0.6,0.6,0.6));

  G4VisAttributes *RoofAtt = new G4VisAttributes(G4Colour(0.3,0.3,0.3));


  G4VisAttributes *GreenBoxAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  GreenBox_log->SetVisAttributes(BlockAtt);

  G4VisAttributes *SteelPlateAtt = new G4VisAttributes(G4Colour(0.0,0.5,1.0));
  SteelPlate_log->SetVisAttributes(RoofAtt);

  G4VisAttributes *GreenPlateAtt = new G4VisAttributes(G4Colour(0.7,0.9,0.3));
  GreenPlate1_log->SetVisAttributes(RoofAtt);
  GreenPlate2_log->SetVisAttributes(RoofAtt);

  G4VisAttributes *ElecAtt = new G4VisAttributes(G4Colour(0.8,0.0,0.0));
  ElecAtt->SetForceWireframe(true);
  Electronics_log->SetVisAttributes(ElecAtt);
}

//Sieve slit
void G4SBSEArmBuilder::MakeBBSieveSlit(G4LogicalVolume *motherlog)
{
  printf("Building BB sieve slit...\n");
}
 
