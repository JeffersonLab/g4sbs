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
#include "G4TwoVector.hh"
#include "G4GenericTrap.hh"
#include "G4Polycone.hh"

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

#include "sbstypes.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

using namespace std;

G4SBSEArmBuilder::G4SBSEArmBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
  fBBang  = 40.0*deg;
  fBBdist = 1.5*m;


  fCerDepth = 92.0*cm;
  fCerDist  =  7.0*cm;

  fBBCaldist = 20*cm + fCerDepth;
  fGEMDist   = 10*cm + fCerDepth;
  fGEMOption = 1;

  fUseLocalField = false;

  assert(fDetCon);

  fbbfield =  NULL;  
}

G4SBSEArmBuilder::~G4SBSEArmBuilder(){;}

void G4SBSEArmBuilder::BuildComponent(G4LogicalVolume *worldlog){
  Exp_t exptype = fDetCon->fExpType;

  //  The neutron experiments and the SIDIS experiment use BigBite:
  //------------ BigBite: -----------------------------------------------------
  if( exptype == kNeutronExp || exptype == kSIDISExp ) 
    {
      MakeBigBite( worldlog );
    }
  if( exptype == kGEp ) //Subsystems unique to the GEp experiment include FPP and BigCal:
    {
      MakeBigCal( worldlog );
    }
}

void G4SBSEArmBuilder::MakeBigBite(G4LogicalVolume *worldlog){
  //Lines of code used to build BigBite moved to their own method:
  printf("BigBite at %f deg\n", fBBang/deg);

  G4RotationMatrix *bbrm = new G4RotationMatrix;
  bbrm->rotateY(fBBang);

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
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR", bbmothercutBox, bbrightcutTrap, rightcutrot2, G4ThreeVector(10*eps, 0.0, -motherdepth/2.0+clear));

  G4Box *frontboxcut = new G4Box("frontboxcut",(bbmagwidth-gapsize)/4.0-coilwidth, 250*cm + eps, clear/2);
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR_fL", bbmothercutBox, frontboxcut, 0, G4ThreeVector( ((bbmagwidth+gapsize)/4.0+coilwidth)+ 10*eps +5*cm, 0, -motherdepth/2.0+clear/2 - 10*eps ));
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR_fLR", bbmothercutBox, frontboxcut, 0, G4ThreeVector( -((bbmagwidth+gapsize)/4.0+coilwidth)-10*eps -5*cm, 0, -motherdepth/2.0+clear/2 - 10*eps ));

  G4Box *bottomboxcut = new G4Box("bottomboxcut",bbmagwidth+eps, 250*cm/2, motherdepth+eps);
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR_fLR_floor", bbmothercutBox, bottomboxcut, 0, G4ThreeVector( 0.0, -1.25*m*2-20*cm, 0.0));


  //   Make logical volume for mother vol and place
  G4LogicalVolume *bbmotherLog=new G4LogicalVolume(bbmothercutBox,GetMaterial("Air"),
						   "bbmotherLog", 0, 0, 0);

  new G4PVPlacement(bbrm, G4ThreeVector(-motherr*sin(fBBang), 0.0, motherr*cos(fBBang)),
		    bbmotherLog, "bbmotherPhys", worldlog, 0,false,0);


  new G4PVPlacement(yokerm,G4ThreeVector(0.0, 0.0, -motherdepth/2.0+clear),
		    bbyokewgapLog, "bbyokewgapPhysical", bbmotherLog, 0,false,0);

  //  Bigbite field log volume
  G4LogicalVolume *bbfieldLog=new G4LogicalVolume(bbairTrap, GetMaterial("Air"),
						  "bbfieldLog", 0, 0, 0);


  fbbfield = new G4SBSBigBiteField( G4ThreeVector(0.0, 0.0, fBBdist), *bbrm );

  fbbfield->fInverted = fDetCon->fGlobalField->fInverted;


  G4FieldManager *bbfm = new G4FieldManager(fbbfield);
  new G4ChordFinder(fbbfield);

  if( fUseLocalField ){
    bbmotherLog->SetFieldManager(bbfm,true);
  }

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
    ngem = 4;
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
  double gemz_opt1[] = { 0.0*cm, gemdsep, fGEMDist, fGEMDist+gemdsep};
  double gemw_opt1[] = { 40.0*cm, 40.0*cm, 50.0*cm, 50.0*cm };
  double gemh_opt1[] = { 150.0*cm, 150.0*cm, 200.0*cm, 200.0*cm };

  // GEM option 2
  double gemz_opt2[] = { 0.0*cm, gemdsep, 2.0*gemdsep, 3.0*gemdsep, fGEMDist};
  double gemw_opt2[] = { 40.0*cm, 40.0*cm, 40.0*cm, 40.0*cm, 50.0*cm };
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

  // **** BIGBITE PRESHOWER **** 
  // 2 columns, 27 rows

  double psheight = 27*8.5*cm;
  double pswidth  = 2.0*37.0*cm;
  double psdepth  = 8.5*cm;

  G4Box *bbpsbox = new G4Box("bbpsbox", pswidth/2.0, psheight/2.0, psdepth/2.0 );
  G4LogicalVolume *bbpslog = new G4LogicalVolume(bbpsbox, GetMaterial("Air"), "bbpslog");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth/2.0), bbpslog, "bbpsphys", bbdetLog, false, 0);

  // Preshower module - geometry will be assigned after Shower

  // **** BIGBITE HODOSCOPE **** 
  // Scintillator box - same dimensions as preshower
  double bbhododepth = 2.5*cm;
  G4Box *bbhodobox = new G4Box("bbhodobox", pswidth/2.0, psheight/2.0, bbhododepth/2.0 );
  G4LogicalVolume *bbhodolog = new G4LogicalVolume( bbhodobox, GetMaterial("PLASTIC_SC_VINYLTOLUENE"), "bbhodolog" );
  new G4PVPlacement(0, G4ThreeVector(0.0,0.0, detoffset+fBBCaldist+psdepth+bbhododepth/2.0), bbhodolog, "bbhodophys", bbdetLog, false, 0);

  // **** BIGBITE SHOWER ****
  // 7 columns, 27 rows
  double calheight = 27*8.5*cm;
  double calwidth  = 7*8.5*cm;
  double caldepth  = 37.0*cm;
  G4Box *bbshowerbox = new G4Box("bbshowerbox", calwidth/2.0, calheight/2.0, caldepth/2.0);
  G4LogicalVolume *bbshowerlog = new G4LogicalVolume(bbshowerbox, GetMaterial("Air"), "bbshowerlog");
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth+bbhododepth+caldepth/2.0), bbshowerlog, "bbshowerphys", bbdetLog, false, 0);

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
  G4LogicalVolume *bbTF1log = new G4LogicalVolume( bbTF1box, GetMaterial("TF1"), "bbTF1log" );
  
  // Shower TF1 SD of type CAL
  G4SDManager *sdman = fDetCon->fSDman;

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

      new G4PVPlacement(0, G4ThreeVector(xtemp,ytemp,0.0), showermodlog, "showermodphys", bbshowerlog, false, shower_copy_number);
      
      (BBSHSD->detmap).LocalCoord[shower_copy_number] = G4ThreeVector( xtemp,ytemp,(caldepth-bbpmtz)/2.0  );
      (BBSHTF1SD->detmap).LocalCoord[shower_copy_number] = G4ThreeVector( xtemp, ytemp, (caldepth-4*bbpmtz-bbTF1_z)/2.0 );
      shower_copy_number++;
    }
  }

  // ****Preshower Continued****
  // Reusing modules from Shower (same variables), rotated by either +/- 90 deg depending on column #
  G4Box *preshowermodbox = new G4Box( "preshowermodbox", bbmodule_x/2.0, bbmodule_y/2.0, caldepth/2.0 );
  G4LogicalVolume *preshowermodlog = new G4LogicalVolume( preshowermodbox, GetMaterial("Special_Air"), "preshowermodlog" );
 
  // Preshower TF1 SD of type CAL
  G4LogicalVolume *bbpsTF1log = new G4LogicalVolume( bbTF1box, GetMaterial("TF1"), "bbpsTF1log" );

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

}

void G4SBSEArmBuilder::MakeBigCal(G4LogicalVolume *worldlog){

  printf("BigCal at %f deg\n", fBBang/deg);

  // Pointer to SDmanager, used frequently in this routine
  G4SDManager *sdman = fDetCon->fSDman;

  // Electron Arm - order of materials - all "right next" to each other in z-hat_spectrometer coordinates
  //   1) 20cm Polyethylene 
  //   2) 2 4cm CDet planes
  //   3) 2.4cm Al
  //   4) ECal

  // Rotations, offsets, and distances from Target:
  G4RotationMatrix *bbrm = new G4RotationMatrix;
  bbrm->rotateY(fBBang);

  //   -some of these variables are used in ECAL Section exclusively.
  G4double x_ecal = 246.402*cm, y_ecal = 370.656*cm, z_ecal = 45.406*cm; //45.006cm(module) + 2*0.2(pmt depth)

  G4double z_Al = 2.40*cm;         //Al depth
  G4double polydepth = 20.00*cm;   //Polyethylene
  G4double CD_depth = 4.00*cm;     //Used in Option 1
  G4double x_earm = 200.070*cm, y_earm = 315.900*cm, z_earm = 30.40*cm; 

  double bbr = fBBdist - z_ecal/2.0; //backside of ECal should be located at fBBdist
  double offset = 15*cm;             //Motivation - match SBS acceptance, used to declare distances
                                     //from target and frequently in ECAL section


  /*************************************************************************************
   ********************************        CDet      ***********************************
   *************************************************************************************/

  // 4 Options/Configurations available:
  //    /g4sbs/CDetconfig 1   -   Just boxes containing nothing but one material
  //    /g4sbs/CDetconfig 2   -   392 bars make up one module, 3 modules per plane. 
  //                              PMTs/Scint are SD, planes are FLAT
  //    /g4sbs/CDetconfig 3   -   Same as 2, except Top/Bottom modules are rotated 
  //                              by some angle CDetPhi such that a line
  //                              running from the origin to the center of a module 
  //                              is perpendular to the face of a module
  //    /g4sbs/CDetconfig 4   -   All bars are rotated such that a line running from 
  //                              the origin to the  center of a "long" bar
  //                              is perpendicular to the face of a CDet bar. A "long" 
  //                              bar is defined at two small bars 
  //                              connected long ways, or one CDet row.
  G4double x_bar_CDet = 0.5*cm;
  G4double y_bar_CDet = 4.0*cm;
  G4double z_bar_CDet = 51*cm;

  G4double mylar_CDet       = 0.025*cm;
  G4double pmt_CDet_depth   = 0.010*cm;

  G4double BoreHoleDiameter = 0.30*cm;
  G4double WLSDiameter      = 0.20*cm; 

  G4int CDetColumns = 2;
  G4int CDetRows    = 196;

  // Setting up the Scintillating Bars by defining a mother volume
  G4Box *CDetBox = new G4Box("CDetBox",(x_bar_CDet+2*mylar_CDet)/2.0,
			     (y_bar_CDet+2*mylar_CDet)/2.0,
			     (z_bar_CDet+mylar_CDet)/2.0 );
  G4LogicalVolume *CDetLog = new G4LogicalVolume( CDetBox, GetMaterial("Special_Air"), "CDetLog" );
  
  G4Box *CDetScintBox_temp = new G4Box("CDetScintBox_temp", x_bar_CDet/2.0, y_bar_CDet/2.0, z_bar_CDet/2.0);

  // Setting up Mylar Wrap
  G4SubtractionSolid *CDetMylarBox = new G4SubtractionSolid( "CDetMylar", CDetBox, CDetScintBox_temp, 0, G4ThreeVector(0.0,0.0,mylar_CDet/2.0) );
  G4LogicalVolume *CDetMylarLog = new G4LogicalVolume( CDetMylarBox, GetMaterial("Mylar"), "CDetMylarLog" );
  new G4LogicalSkinSurface( "CDet_Mylar_Skin", CDetMylarLog, GetOpticalSurface("Mirrsurf") );
    
  // Bore out cylindrical section of our scintillating material, CDetScintBox.
  G4Tubs *BoreHole = new G4Tubs( "BoreHole", 0.0, BoreHoleDiameter/2.0, 1.001*z_bar_CDet/2.0, 0.0, twopi );
  G4SubtractionSolid *CDScintBoreBox = new G4SubtractionSolid( "CDScintBoreBox", CDetScintBox_temp, BoreHole, 0, G4ThreeVector() );
  G4LogicalVolume *CDScintBoreLog = new G4LogicalVolume( CDScintBoreBox, GetMaterial("EJ232"), "CDScintBoreLog" );

  // Make a Lightguide
  G4Tubs *WLSguide = new G4Tubs( "WLSguide", 0.0, WLSDiameter/2.0, z_bar_CDet/2.0, 0.0, twopi );
  G4LogicalVolume *WLSguideLog = new G4LogicalVolume( WLSguide, GetMaterial("BC484"), "WLSguideLog" );

  // Place these materials inside CDetBox.
  new G4PVPlacement(0, G4ThreeVector(), CDetMylarLog, "CDetMylarLogPhys", CDetLog, false, 0); 
  new G4PVPlacement(0, G4ThreeVector(0.0,0.0,mylar_CDet/2.0), CDScintBoreLog, "CDScintBorePhys", CDetLog, false, 0);
  new G4PVPlacement(0, G4ThreeVector(0.0,0.0,mylar_CDet/2.0), WLSguideLog, "WLSguidePhys", CDetLog, false, 0);

  // Make some PMT's
  G4Tubs *CDetPMT = new G4Tubs( "CDetPMT", 0.0, WLSDiameter/2.0, pmt_CDet_depth/2.0, 0.0, twopi );
  G4LogicalVolume *CDetPMTLog = new G4LogicalVolume( CDetPMT, GetMaterial("Photocathode_CDet"), "CDetPMTLog" );

  // CDetPMTLog is the SD, assigned to ECalSD which detects optical photons
  G4String CDetSDname = "Earm/CDet";
  G4String CDetcollname = "CDetHitsCollection";
  G4SBSECalSD *CDetSD = NULL;

  if( !( (G4SBSECalSD*) sdman->FindSensitiveDetector(CDetSDname)) ) {
    G4cout << "Adding CDet PMT Sensitive Detector to SDman... " << G4endl;
    CDetSD = new G4SBSECalSD( CDetSDname, CDetcollname );
    sdman->AddNewDetector( CDetSD );
    (fDetCon->SDlist).insert( CDetSDname );
    fDetCon->SDtype[CDetSDname] = kECAL;
    //fDetCon->SDarm[CDetSDname] = kEarm;
    //(CDetSD->detmap).depth = 1;            // needs to be updated!!!!! and changed for each option most likely
  }
  CDetPMTLog->SetSensitiveDetector( CDetSD );
    
  // Make the Final CDet Bar, housing CDetLog and CDetPMTLog!
  G4double CDetBarX = x_bar_CDet + 2*mylar_CDet;
  G4double CDetBarY = y_bar_CDet + 2*mylar_CDet;
  G4double CDetBarZ = z_bar_CDet + mylar_CDet + pmt_CDet_depth;
  G4Box *CDetBarBox = new G4Box( "CDetBarBox", CDetBarX/2.0, CDetBarY/2.0, CDetBarZ/2.0 );
  G4LogicalVolume *CDetBarLog = new G4LogicalVolume( CDetBarBox, GetMaterial("Air"), "CDetBarLog" );

  // Placing..
  new G4PVPlacement(0, G4ThreeVector(0.0,0.0,(CDetBarZ-pmt_CDet_depth)/2.0), CDetPMTLog,"CDetPMTPhys", CDetBarLog, false, 0);
  new G4PVPlacement(0, G4ThreeVector(0.0,0.0,-pmt_CDet_depth/2.0), CDetLog, "CDetPhys", CDetBarLog, false, 0);

  G4RotationMatrix *CDetBarRot = new G4RotationMatrix;
  CDetBarRot->rotateZ(90*deg);
  //CDetBarRot->rotateX(90*deg);
  //TEST to see if module looks like I expect
  //new G4PVPlacement(CDetBarRot,G4ThreeVector(-5*m, 0*m, 3.0*m), CDetBarLog, "Ex", worldlog, false, 0,true );

  // Lets define a CDet MODULE - 2 cols, 196 rows => 392 bars in 1 module:
  G4double CDetModuleX = CDetColumns*CDetBarZ; 
  G4double CDetModuleY = CDetRows*CDetBarX;
  G4double CDetModuleZ = CDetBarY;
  G4Box *CDetModuleBox = new G4Box( "CDetModuleBox", CDetModuleX/2.0, CDetModuleY/2.0, CDetModuleZ/2.0 );
  G4LogicalVolume *CDetModuleLog = new G4LogicalVolume( CDetModuleBox, GetMaterial("Air"), "CDetModuleLog" );

  // Place Bars in a Module, rotating based on Column. 
  G4RotationMatrix *CDetBar_Col1 = new G4RotationMatrix;
  CDetBar_Col1->rotateZ(90.0*deg);
  CDetBar_Col1->rotateX(90.0*deg);
  G4RotationMatrix *CDetBar_Col2 = new G4RotationMatrix;
  CDetBar_Col2->rotateZ(-90.0*deg);
  CDetBar_Col2->rotateX(90.0*deg);
  G4int cdet_copy_number = 0;
  
  // Make modules for Options 2 & 3
  if( fDetCon->GetCDetConfigOption() == 2 || fDetCon->GetCDetConfigOption() == 3 ) {
    (CDetSD->detmap).depth = 1; 
    for(int l=0; l<CDetColumns; l++) {
      for(int j=0; j<CDetRows; j++) {
	double xtemp = (CDetModuleX - CDetBarZ)/2.0 - l*CDetBarZ;
	double ytemp = (CDetModuleY - CDetBarX)/2.0 - j*CDetBarX;

	(CDetSD->detmap).Col[cdet_copy_number] = l;
	(CDetSD->detmap).Row[cdet_copy_number] = j;
	if(l==0) {
	  new G4PVPlacement( CDetBar_Col1, G4ThreeVector(xtemp,ytemp,0.0), CDetBarLog, "CDetBarPhys", CDetModuleLog, false, cdet_copy_number );
	  (CDetSD->detmap).LocalCoord[cdet_copy_number] = G4ThreeVector(xtemp+CDetBarZ/2.0,ytemp,0.0);
	  cdet_copy_number++;
	}
	if(l==1) {
	  new G4PVPlacement( CDetBar_Col2, G4ThreeVector(xtemp,ytemp,0.0), CDetBarLog, "CDetBarPhys", CDetModuleLog, false, cdet_copy_number );
	  (CDetSD->detmap).LocalCoord[cdet_copy_number] = G4ThreeVector(xtemp-CDetBarZ/2.0,ytemp,0.0);
	  cdet_copy_number++;
	}
      }
    }
  }
  // 3 Modules make up a "Plane" and CDet consists of two planes:
  G4Box *AlBox = new G4Box( "AlBox", x_earm/2.0, y_earm/2.0, z_Al/2.0 );
  G4LogicalVolume *Allog = new G4LogicalVolume ( AlBox, GetMaterial("RICHAluminum"), "Allog" );

  G4double PlaneX = x_earm;
  G4double PlaneY = 3*CDetModuleY; //no 2
  G4double PlaneZ = CDetModuleZ + CDetModuleY*sin(40*deg); //sin25

  // Mother Volume - houses 2 planes, Al, Polyethylene
  G4Box *PlaneBox = new G4Box( "PlaneBox", PlaneX/2.0, PlaneY/2.0, PlaneZ/2.0 );
  G4LogicalVolume *PlaneLog = new G4LogicalVolume( PlaneBox, GetMaterial("Air"), "PlaneLog" );

  // Initialize some Geometry Variables:
  double Layeroffset = 0.0;                // Space inbetween planes, needed for Options 3 & 4
  double TopBottomModuleOffset = 0.0;      // Offset top/bottom module
  double MiddleModuleOffset    = 0.0;      // Offset middle module
  double CDetPhi = 0.0;                    // Angle that rotates bars or modules
  double CDetFrontRad = 0.0;               // Distance between origin & front of CDet - used for rotations of modules/bars                          
  double cdr = 0.0;                        // Distance between origin and center of CDet - used to place mother

  // Make the modules for Option 4, redefine some variables due to rotations
  double CDetBarX4 = 2*CDetBarZ; 
  double CDetBarY4 = CDetBarX;
  double CDetBarZ4 = CDetBarY;
  G4Box *CDetModuleBox4 = new G4Box( "CDetModuleBox4", CDetBarX4/2.0, CDetBarY4/2.0, CDetBarZ4/2.0 );
  G4LogicalVolume *CDetModuleLog4 = new G4LogicalVolume( CDetModuleBox4, GetMaterial("Air"), "CDetModuleLog4");
  // Place two bars in one of these modules
  new G4PVPlacement(CDetBar_Col1, G4ThreeVector((CDetBarX4-CDetBarZ)/2.0,0.0,0.0), CDetBarLog, "CDetBarPhys4", CDetModuleLog4, false, 0, true);
  new G4PVPlacement(CDetBar_Col2, G4ThreeVector((-CDetBarX4+CDetBarZ)/2.0,0.0,0.0), CDetBarLog, "CDetBarPhys4", CDetModuleLog4, false, 0, true);

  int CDetOpt4_copyn = 0;

  if( fDetCon->GetCDetConfigOption() == 4 ) {
    (CDetSD->detmap).depth = 2;  
    cdr = fBBdist - z_ecal - PlaneZ/2.0;  
    CDetFrontRad = cdr + PlaneZ/2.0 - CDetBarZ4/2.0;
    CDetPhi = 2*atan( CDetBarY4/(2*CDetFrontRad) );

    // start at center, work up/down:
    for(int j=0; j<(1.5*CDetRows); j++) { 

      G4RotationMatrix *CDetBar4Rottop = new G4RotationMatrix;
      double angle = j*CDetPhi;
      CDetBar4Rottop->rotateX(angle);

      double dispy = ((CDetBarZ4/2.0)*sin( angle ) - (CDetBarY4/2.0)*(1-cos( angle ))); 
      double dispz = ((CDetBarZ4/2.0)*(1-cos( angle )) + (CDetBarY4/2.0)*sin( angle ));
      double ytemp = j*CDetBarX;
      double ztemp = PlaneZ/2.0 - z_Al - CDetBarZ4/2.0;

      // Additional placement correction since modules are placed relative to top left (TL) corner of previous module:
      // currently this is just an approximation
      double TLZ = 0.5*j*CDetBarY4*sin(angle)*cos(angle);
      double TLY = 0.1*j*CDetBarY4*sin(angle)*sin(angle);
      new G4PVPlacement(CDetBar4Rottop, G4ThreeVector(0.0,ytemp+dispy-TLY,-dispz-TLZ+ztemp), CDetModuleLog4, "CDetBarPhys4",PlaneLog, false, CDetOpt4_copyn,true);
      // odd rows including 0
      (CDetSD->detmap).LocalCoord[CDetOpt4_copyn] = G4ThreeVector(0.0,ytemp+dispy-TLY,-dispz-TLZ+ztemp);  
      (CDetSD->detmap).Row[CDetOpt4_copyn] = j;                                                            
      CDetOpt4_copyn++;

      G4RotationMatrix *CDetBar4Rotbot = new G4RotationMatrix;
      angle = -j*CDetPhi;
      CDetBar4Rotbot->rotateX(angle);
      if(j>0) {
	new G4PVPlacement(CDetBar4Rotbot, G4ThreeVector(0.0,-(ytemp+dispy-TLY),-dispz-TLZ+ztemp), CDetModuleLog4, "CDetBarPhys4",PlaneLog, false, CDetOpt4_copyn,true);
	// even rows excluding 0
	(CDetSD->detmap).LocalCoord[CDetOpt4_copyn] = G4ThreeVector(0.0,-(ytemp+dispy-TLY),-dispz-TLZ+ztemp);  
	(CDetSD->detmap).Row[CDetOpt4_copyn] = j;
	CDetOpt4_copyn++;
      }
    }

    // Place Option 4
    G4ThreeVector CDetPos( cdr*sin(-fBBang)+(offset+ (4.212/2.0)*cm )*cos(fBBang), 
			   -(4.212/2.0)*cm, 
			   cdr*cos(-fBBang)+(offset+(4.212/2.0)*cm)*sin(fBBang) );
    
    new G4PVPlacement( bbrm, CDetPos, PlaneLog, "CDet", worldlog, false, 0); \
    new G4PVPlacement( 0, G4ThreeVector(0.0,0.0,(PlaneZ-z_Al)/2.0),Allog,"AlOpt4",PlaneLog,false,0);
  }

  // Place Options 2 & 3
  if( fDetCon->GetCDetConfigOption() == 2 || fDetCon->GetCDetConfigOption() == 3 ) {
    cdr = fBBdist - z_ecal - PlaneZ/2.0;                   //distance from origin to center of CDet Plane
    CDetFrontRad = cdr + PlaneZ/2.0 - CDetModuleZ;         //origin to CDet Front
  
    switch( fDetCon->GetCDetConfigOption() ) {
    case 2:
      CDetPhi = 0.0;
      Layeroffset = 0.0;
      break;
    case 3:
      CDetPhi = 2.0*atan( CDetModuleY/(2.0*CDetFrontRad) );
      Layeroffset = 0.30*cm;
      break;
    default:
      CDetPhi = 0.0;
      Layeroffset = 0.0;
      break;
    }
    G4RotationMatrix *CDetModuleRotTop = new G4RotationMatrix;
    CDetModuleRotTop->rotateX(CDetPhi);

    G4RotationMatrix *CDetModuleRotBottom = new G4RotationMatrix;
    CDetModuleRotBottom->rotateX(-CDetPhi);

    G4ThreeVector CDetPos( cdr*sin(-fBBang)+(offset+ (4.212/2.0)*cm )*cos(fBBang), 
    			   -(4.212/2.0)*cm, 
    			   cdr*cos(-fBBang)+(offset+(4.212/2.0)*cm)*sin(fBBang) );

    new G4PVPlacement( bbrm, CDetPos, PlaneLog, "CDet", worldlog, false, 0);

    double DisplaceZ = (CDetModuleY/2.0)*sin(CDetPhi) + (CDetModuleZ/2.0)*( 1-cos(CDetPhi) ); 
    double DisplaceY = (CDetModuleY/2.0)*( 1-cos(CDetPhi) ) - (CDetModuleZ/2.0)*sin(CDetPhi);

    G4ThreeVector postop(TopBottomModuleOffset, CDetModuleY-DisplaceY, (PlaneZ-CDetModuleZ)/2.0-DisplaceZ-z_Al);
    G4ThreeVector posbot(TopBottomModuleOffset, -CDetModuleY+DisplaceY, (PlaneZ-CDetModuleZ)/2.0-DisplaceZ-z_Al);

    //Layer 2 (closest to ECal)
    new G4PVPlacement( 0,G4ThreeVector(0.0,0.0,(PlaneZ-z_Al)/2.0), Allog, "AlShield", PlaneLog, false, 0 );

    new G4PVPlacement( CDetModuleRotTop, postop, CDetModuleLog, "CDetModulePhys", PlaneLog, false, 1, true );
    new G4PVPlacement( 0, G4ThreeVector(-MiddleModuleOffset, 0.0, PlaneZ/2.0-CDetModuleZ/2.0-z_Al),
		       CDetModuleLog, "CDetModulePhys", PlaneLog, false, 2, true); 
    new G4PVPlacement( CDetModuleRotBottom, posbot, CDetModuleLog, "CDetModulePhys", PlaneLog, false, 3, true );

    //Layer 1(closest to Target)
    CDetFrontRad = cdr + PlaneZ/2.0 - 2*CDetModuleZ; 

    switch( fDetCon->GetCDetConfigOption() ) {
    case 2:
      CDetPhi = 0.0;
      Layeroffset = 0.0;
      break;
    case 3:
      CDetPhi = 2.0*atan( CDetModuleY/(2.0*CDetFrontRad) );
      Layeroffset = 0.30*cm;
      break;
    default:
      CDetPhi = 0.0;
      Layeroffset = 0.0;
      break;
    }
  
    G4RotationMatrix *CDetModuleRotTopLayer1 = new G4RotationMatrix;
    CDetModuleRotTopLayer1->rotateX(CDetPhi);
    G4RotationMatrix *CDetModuleRotBotLayer1 = new G4RotationMatrix;
    CDetModuleRotBotLayer1->rotateX(-CDetPhi);

    DisplaceZ = (CDetModuleY/2.0)*sin(CDetPhi) + (CDetModuleZ/2.0)*( 1-cos(CDetPhi) ); 
    DisplaceY = (CDetModuleY/2.0)*( 1-cos(CDetPhi) ) - (CDetModuleZ/2.0)*sin(CDetPhi);

    postop.set(TopBottomModuleOffset, CDetModuleY-DisplaceY, (PlaneZ-3*CDetModuleZ)/2.0-DisplaceZ-Layeroffset-z_Al);
    posbot.set(TopBottomModuleOffset, -CDetModuleY+DisplaceY, (PlaneZ-3*CDetModuleZ)/2.0-DisplaceZ-Layeroffset-z_Al);

    new G4PVPlacement( CDetModuleRotTopLayer1, postop, CDetModuleLog, "CDetModulePhys", PlaneLog, false, 4, true );
    new G4PVPlacement( 0, G4ThreeVector(-MiddleModuleOffset, 0.0, (PlaneZ-3*CDetModuleZ)/2.0 - Layeroffset-z_Al),
		       CDetModuleLog, "CDetModulePhys", PlaneLog, false, 5, true); 
    new G4PVPlacement( CDetModuleRotBotLayer1, posbot, CDetModuleLog, "CDetModulePhys", PlaneLog, false, 6, true );
  }

  // VISUALIZATIONS
  CDetLog->SetVisAttributes( G4VisAttributes::Invisible );
  CDetBarLog->SetVisAttributes( G4VisAttributes::Invisible );
  CDetModuleLog->SetVisAttributes( G4VisAttributes::Invisible );
  //PlaneLog->SetVisAttributes( G4VisAttributes::Invisible );
  CDetModuleLog4->SetVisAttributes( G4VisAttributes::Invisible );

  G4VisAttributes *CDetPlaneVis = new G4VisAttributes( G4Colour::Green() );
  CDetPlaneVis->SetForceWireframe(true);
  PlaneLog->SetVisAttributes(CDetPlaneVis);

  G4VisAttributes *CDetMod4Vis = new G4VisAttributes( G4Colour::Green() );
  CDetMod4Vis->SetForceWireframe(true);
  //CDetModuleLog4->SetVisAttributes(CDetMod4Vis);
 
  // Mylar
  G4VisAttributes *CDetMylarVis = new G4VisAttributes( G4Colour(0.5,0.5,0.5) );
  CDetMylarLog->SetVisAttributes( CDetMylarVis );
  //CDetMylarLog->SetVisAttributes( G4VisAttributes::Invisible );

  // Waveguide
  G4VisAttributes *CDetWLSVis = new G4VisAttributes( G4Colour(0.0, 1.0, 1.0) );
  WLSguideLog->SetVisAttributes( CDetWLSVis );

  // Scintillating Material
  G4VisAttributes *CDetScintVis = new G4VisAttributes( G4Colour( 0.5, 0.0, 0.8) );
  CDScintBoreLog->SetVisAttributes( CDetScintVis );
  
  // PMT
  G4VisAttributes *CDetPMTVis = new G4VisAttributes( G4Colour(0.3,0.3,0.3) );
  CDetPMTLog->SetVisAttributes( CDetPMTVis );
	
  // Al
  G4VisAttributes *AlVis = new G4VisAttributes(G4Colour(0.9, 0.9, 0.9));
  AlVis->SetForceWireframe(true);
  Allog->SetVisAttributes(AlVis);	    

  /*************************************************************************************
   ********************************        ECAL      ***********************************
   *************************************************************************************/
 
  //Dimensions are different than ecal_log because ECal was shrunk significantly when the crescent shape was developed
  //x_earm = 47.5*4.212 (the .5 comes from staggering effects)
  //y_earm = 75*4.212
  //z_earm = 20 + 2*4 + 2.40

  // Make CDet Option 1
  G4Box *earm_mother_box = new G4Box("earm_mother_box", x_earm/2.0, y_earm/2.0, z_earm/2.0);
  G4LogicalVolume *earm_mother_log = new G4LogicalVolume(earm_mother_box, GetMaterial("Air"), "earm_mother_log");

  // Polyethylene Box
  G4Box *polybox = new G4Box("polybox", x_earm/2.0, y_earm/2.0, polydepth/2.0 );
  G4LogicalVolume *polybox_log = new G4LogicalVolume( polybox, GetMaterial("Polyethylene"), "polybox_log");

  // Coordinate Detector - 2 planes
  G4Box *CD_box = new G4Box("CD_box", x_earm/2.0, y_earm/2.0, CD_depth/2.0 );
  G4LogicalVolume *CD_log = new G4LogicalVolume(CD_box, GetMaterial("PLASTIC_SC_VINYLTOLUENE"), "CD_log");

  // Aluminum Shielding in front of ECAL - make it big enough to enclose crescent
  G4Box *Al_box = new G4Box( "Al_box", x_earm/2.0, y_earm/2.0, z_Al/2.0 );
  G4LogicalVolume *Al_log = new G4LogicalVolume ( Al_box, GetMaterial("RICHAluminum"), "Al_log" );

  // Place everything inside the Earm Mother Volume
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, z_earm/2.0 - z_Al/2.0 ), Al_log, "Al", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, z_earm/2.0 - z_Al - CD_depth/2.0 ), CD_log, "plane1", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, z_earm/2.0 - z_Al - CD_depth - CD_depth/2.0 ), CD_log, "plane2", earm_mother_log, false, 1 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, z_earm/2.0 - z_Al - 2*CD_depth - polydepth/2.0 ), polybox_log, "Polyethylene", earm_mother_log, false, 0 );

  if( fDetCon->GetCDetConfigOption() == 1 ) {
    cdr = fBBdist - z_ecal - z_earm/2.0;
    //Note the offsets in x,y,z of 4.212/2.0
    new G4PVPlacement( bbrm, G4ThreeVector( cdr*sin(-fBBang)+(offset + (4.212/2.0)*cm )*cos(fBBang), 
					    -(4.212/2.0)*cm, 
					    cdr*cos(-fBBang) + (offset+(4.212/2.0)*cm)*sin(fBBang) ), 
		       earm_mother_log, "CDet & Filter", worldlog, false, 0);
  }
  if( fDetCon->GetCDetConfigOption() == 2 ) {
    cdr = fBBdist - z_ecal - PlaneZ/2.0;  
    new G4PVPlacement( 0, G4ThreeVector(0.0,0.0, PlaneZ/2.0 - z_Al - 2*CDetModuleZ - polydepth/2.0), 
		       polybox_log, "Polyphys", PlaneLog, false, 0 );
  }
  
  //Define a Mother Volume to place the ECAL modules
  G4Box *ecal_box = new G4Box( "ecal_box", x_ecal/2.0, y_ecal/2.0, z_ecal/2.0 );
  G4LogicalVolume *ecal_log = new G4LogicalVolume( ecal_box, GetMaterial("Air"), "ecal_log" );
  new G4PVPlacement( bbrm, G4ThreeVector( bbr*sin(-fBBang)+offset*cos(fBBang), 0.0, bbr*cos(-fBBang)+offset*sin(fBBang) ), ecal_log, "ECal Mother", worldlog, false, 0 );

  //Dimensions of mylar/air wrapping:
  G4double mylar_wrapping_size = 0.00150*cm;
  G4double air_wrapping_size = 0.00450*cm;
  G4double mylar_plus_air = mylar_wrapping_size + air_wrapping_size;

  //TYPE1 MODULE - 4.200 x 4.200 x 45.000 cm^3 + layer of air (0.0045cm) + layer of mylar (0.0015cm)
  G4double x_module_type1 = 4.2000*cm + (2*mylar_plus_air); //4.212
  G4double y_module_type1 = 4.2000*cm + (2*mylar_plus_air); //4.212
  G4double z_module_type1 = 45.0000*cm + mylar_plus_air;    //45.006 - only one side in z has mylar+air
  G4Box *module_type1 = new G4Box( "module_type1", x_module_type1/2.0, y_module_type1/2.0, z_module_type1/2.0 );

  //Define a new mother volume which will house the contents of a module(i.e. TF1 & Mylar)
  G4LogicalVolume *module_log_type1 = new G4LogicalVolume( module_type1, GetMaterial("Special_Air"), "module_log_type1" );

  //MYLAR
  //Subtraction solid to leave us with a 0.0015cm "mylar wrapping" with one end open
  G4double x_sub = x_module_type1 - (2*mylar_wrapping_size);
  G4double y_sub = y_module_type1 - (2*mylar_wrapping_size);
  G4double z_sub = z_module_type1; //going to be shifted in G4SubtractionSolid in order to get 0.0015cm of mylar on one end
  G4Box *mylar_module1_subtract = new G4Box( "mylar_module1_subtract", x_sub/2.0, y_sub/2.0, z_sub/2.0 );
  G4SubtractionSolid *mylarwrap = new G4SubtractionSolid( "mylarwrap", module_type1, mylar_module1_subtract, 0, G4ThreeVector(0.0, 0.0, mylar_wrapping_size));
  G4LogicalVolume *mylar_wrap_log = new G4LogicalVolume( mylarwrap, GetMaterial("Mylar"), "mylar_wrap_log" );

  //TF1
  G4double x_TF1 = 4.200*cm, y_TF1 = 4.200*cm, z_TF1 = 45.000*cm;
  G4Box *TF1_box = new G4Box( "TF1_box", x_TF1/2.0, y_TF1/2.0, z_TF1/2.0 );
  G4LogicalVolume *TF1_log = new G4LogicalVolume ( TF1_box, GetMaterial("TF1_anneal"), "TF1_log" );

  G4String ECalTF1SDname = "Earm/ECalTF1";
  G4String ECalTF1collname = "ECalTF1HitsCollection";
  G4SBSCalSD *ECalTF1SD = NULL;
    
  if( !( (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(ECalTF1SDname) ) ){
    G4cout << "Adding ECal TF1 Sensitive Detector to SDman..." << G4endl;
    ECalTF1SD = new G4SBSCalSD( ECalTF1SDname, ECalTF1collname );
    fDetCon->fSDman->AddNewDetector( ECalTF1SD );
    (fDetCon->SDlist).insert(ECalTF1SDname);
    fDetCon->SDtype[ECalTF1SDname] = kCAL;
    //fDetCon->SDarm[ECalTF1SDname] = kEarm;

    (ECalTF1SD->detmap).depth = 1;
  }
  TF1_log->SetSensitiveDetector( ECalTF1SD );

  if( (fDetCon->StepLimiterList).find( ECalTF1SDname ) != (fDetCon->StepLimiterList).end() ){
    TF1_log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }

  //Place TF1 mother & Mylar inside module_log_type1 which is already full of Special_Air
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, 0.0), mylar_wrap_log, "Mylar_Wrap", module_log_type1, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, mylar_plus_air/2.0), TF1_log, "TF1", module_log_type1, false, 0 );

  //Want blocks to be optically independent => assign an optical skin to the Mylar
  new G4LogicalSkinSurface( "Mylar_skin", mylar_wrap_log, GetOpticalSurface( "Mirrsurf" ));
    
  //Build PMT detector which will be placed at the end of the module
  //Starting with the Quartz window
  G4double PMT_window_radius = 1.25*cm;
  G4double PMT_window_depth = 0.20*cm;
  G4Tubs *PMT_window = new G4Tubs( "PMT_window", 0.0*cm, PMT_window_radius, PMT_window_depth/2.0, 0.0, twopi );
  G4LogicalVolume *PMT_window_log = new G4LogicalVolume( PMT_window, GetMaterial("QuartzWindow_ECal"), "PMT_window_log" );

  //PMT
  G4double PMT_radius = 1.25*cm;
  G4double PMT_depth = 0.20*cm;
  G4Tubs *PMTcathode_ecal = new G4Tubs( "PMTcathode_ecal", 0.0*cm, PMT_radius, PMT_depth/2.0, 0.0, twopi );
  G4LogicalVolume *PMTcathode_ecal_log = new G4LogicalVolume( PMTcathode_ecal, GetMaterial("Photocathode_material_ecal"), "PMTcathode_ecal_log" );

  //PMTcathode_ecal_log is the sensitive detector, assigned to ECalSD which detects optical photons
  G4String ECalSDname = "Earm/ECal";
  G4String ECalcollname = "ECalHitsCollection";
  G4SBSECalSD *ECalSD = NULL;

  if( !( (G4SBSECalSD*) sdman->FindSensitiveDetector(ECalSDname) ) ){
    G4cout << "Adding ECal PMT Sensitive Detector to SDman..." << G4endl;
    ECalSD = new G4SBSECalSD( ECalSDname, ECalcollname );
    sdman->AddNewDetector( ECalSD );
    (fDetCon->SDlist).insert(ECalSDname);
    fDetCon->SDtype[ECalSDname] = kECAL;
    //fDetCon->SDarm[ECalSDname] = kEarm;
    (ECalSD->detmap).depth = 0;
  }
  PMTcathode_ecal_log->SetSensitiveDetector( ECalSD );

  G4int x_number_ecal = 58;
  G4int y_number_ecal = 88;   

  int maxrow = -1;
  int maxcol = -1;
  int minrow = y_number_ecal;
  int mincol = x_number_ecal;

  G4String filename = "database/";
  filename += fDetCon->GetECALmapfilename();
  ifstream ecal_map_file;
  
  if( fDetCon->GetECALmapfilename() != "" ){
    ecal_map_file.open(filename.data());
  }
 
  set<pair<int,int> > active_cells;

  G4String Line;

  if( ecal_map_file.is_open() ){
    while ( Line.readLine( ecal_map_file ) ){
      //G4cout << "Read Line from ecal map file: \"" <<  Line.data() << "\"" << G4endl;
      if( !(Line[0] == '#') ){
	int row, col;
	sscanf( Line.data(), "%d %d", &row, &col );
	if( active_cells.size() == 0 || row > maxrow ) maxrow = row;
	if( active_cells.size() == 0 || row < minrow ) minrow = row;
	if( active_cells.size() == 0 || col > maxcol ) maxcol = col;
	if( active_cells.size() == 0 || col < mincol ) mincol = col;
	active_cells.insert(make_pair(row,col));
      } 
    }
  } else {
    for( G4int i=0; i<x_number_ecal; i++){
      for( G4int j=0; j<y_number_ecal; j++){
	active_cells.insert(make_pair(j,i));
	int row=j, col=i; 
	if( active_cells.size() == 0 || row > maxrow ) maxrow = row;
	if( active_cells.size() == 0 || row < minrow ) minrow = row;
	if( active_cells.size() == 0 || col > maxcol ) maxcol = col;
	if( active_cells.size() == 0 || col < mincol ) mincol = col;
      }
    }
  }

  G4double x_position = x_ecal/2.0-x_module_type1/2.0 , y_position = y_ecal/2.0-y_module_type1/2.0;
  
  G4int copy_number_PMT = 0;  //label modules

  //Need a Steel module to fill voids
  G4int x_steel = x_module_type1, y_steel = y_module_type1, z_steel = z_module_type1;
  G4Box *steel_box = new G4Box("steel_box", x_steel/2.0, y_steel/2.0, z_steel/2.0);
  G4LogicalVolume *steel_log = new G4LogicalVolume( steel_box, GetMaterial("Steel"), "steel_log");

  //And a half Steel module to fill in voids caused by staggering
  G4Box *steel_box_half = new G4Box("steel_box_half", 0.5*(x_steel/2.0), y_steel/2.0, z_steel/2.0);
  G4LogicalVolume *steel_log_half = new G4LogicalVolume( steel_box_half, GetMaterial("Steel"), "steel_log");
  //int module_number; //counts modules if needed

  //Iterating to build ECal 
 	
  for( G4int i=0; i<y_number_ecal; i++ ){
    G4double ytemp = y_position - i*(y_module_type1);

    for(G4int j=0; j<x_number_ecal; j++){	
      G4double xtemp_even = x_position - j*(x_module_type1);
      G4double xtemp_odd = x_position - x_module_type1/2.0 - j*(x_module_type1);
	  
      pair<int,int> rowcol( i, j );
	  
      if( active_cells.find( rowcol ) != active_cells.end() ){ //This row and column are active:
	(ECalSD->detmap).Row[copy_number_PMT] = i;
	(ECalSD->detmap).Col[copy_number_PMT] = j;
	(ECalTF1SD->detmap).Row[copy_number_PMT] = i;
	(ECalTF1SD->detmap).Col[copy_number_PMT] = j;
	//EVEN ROWS
	if(i%2==0){ 
	  //Type1 Module
	  new G4PVPlacement(0, G4ThreeVector( xtemp_even, ytemp, -PMT_depth), module_log_type1, "Type1Module", ecal_log, false, copy_number_PMT );
	  //PMT_window
	  new G4PVPlacement( 0, G4ThreeVector(xtemp_even, ytemp, z_ecal/2.0-PMT_depth-PMT_window_depth/2.0), PMT_window_log,"PMT_window_pv", ecal_log, false, copy_number_PMT );
	  //PMT_cathode
	  new G4PVPlacement( 0, G4ThreeVector(xtemp_even, ytemp, z_ecal/2.0-PMT_depth/2.0), PMTcathode_ecal_log, "PMTcathode_pv", ecal_log, false, copy_number_PMT );
	      
	  (ECalSD->detmap).LocalCoord[copy_number_PMT] = G4ThreeVector(xtemp_even,ytemp,z_ecal/2.0-PMT_depth/2.0);
	  //(ECalSD->detmap).GlobalCoord[copy_number_PMT] = G4ThreeVector(0.0,0.0,0.0);
	  (ECalTF1SD->detmap).LocalCoord[copy_number_PMT] = G4ThreeVector(xtemp_even,ytemp,-PMT_depth);
	  //(ECalTF1SD->detmap).GlobalCoord[copy_number_PMT] = G4ThreeVector(0.0,0.0,0.0);
	}

	//ODD ROWS
	else{
	  //Type1 Module
	  new G4PVPlacement( 0, G4ThreeVector( xtemp_odd, ytemp, -PMT_depth), module_log_type1, "Type1Module", ecal_log, false, copy_number_PMT );
	  //PMT_window
	  new G4PVPlacement( 0, G4ThreeVector(xtemp_odd, ytemp, z_ecal/2.0-PMT_depth-PMT_window_depth/2.0), PMT_window_log,"PMT_window_pv", ecal_log, false, copy_number_PMT );
	  //PMT_cathode
	  new G4PVPlacement( 0, G4ThreeVector(xtemp_odd, ytemp, z_ecal/2.0-PMT_depth/2.0), PMTcathode_ecal_log, "PMTcathode_pv", ecal_log, false, copy_number_PMT );
	      
	  (ECalSD->detmap).LocalCoord[copy_number_PMT] = G4ThreeVector(xtemp_odd,ytemp,z_ecal/2.0-PMT_depth/2.0);
	  //(ECalSD->detmap).GlobalCoord[copy_number_PMT] = G4ThreeVector(0.0,0.0,0.0);
	  (ECalTF1SD->detmap).LocalCoord[copy_number_PMT] = G4ThreeVector(xtemp_odd,ytemp,-PMT_depth);
	  //(ECalTF1SD->detmap).GlobalCoord[copy_number_PMT] = G4ThreeVector(0.0,0.0,0.0);
	}
	//	module_number++;
	copy_number_PMT++;
      } else { //Steel Filler:
	if( i>=minrow-2 && i<=maxrow+2 && j>=mincol-2 && j<=maxcol+2 ){
	  if(i%2==0){
	    new G4PVPlacement(0, G4ThreeVector(xtemp_even, ytemp, -PMT_depth ), steel_log, "steel_pv", ecal_log, false, 0);
	    
	  } else {
	    new G4PVPlacement(0, G4ThreeVector(xtemp_odd, ytemp, -PMT_depth ), steel_log, "steel_pv", ecal_log, false, 0);   
	  }
	}
      }
    }

    if( i >= minrow-2 && i <= maxrow+2 ){

      int colleft = max(0,mincol-2);
      int colright = min(maxcol+2,x_number_ecal-1);

      double x_left = x_position - (colleft)*x_module_type1 + x_steel/4.0;
      double x_right = x_position - (colright+0.5)*x_module_type1 - x_steel/4.0;

      if( i%2 == 0 ){
	// new G4PVPlacement(0, G4ThreeVector( (x_earm/2.0)-(x_steel/4.0), ytemp, -PMT_depth ), steel_log_half, "steel_pv", ecal_log, false, 0);
	new G4PVPlacement(0, G4ThreeVector( x_right, ytemp, -PMT_depth ), steel_log_half, "steel_pv", ecal_log, false, 0);
      } else {
	// new G4PVPlacement(0, G4ThreeVector(-(x_earm/2.0)+(x_steel/4.0), ytemp, -PMT_depth ), steel_log_half, "steel_pv", ecal_log, false, 0);
	new G4PVPlacement(0, G4ThreeVector( x_left, ytemp, -PMT_depth ), steel_log_half, "steel_pv", ecal_log, false, 0);
      }
    }
  }
  //END MODULE PLACEMENT

  //VISUALIZATIONS
  //Make mother volume invisible
  ecal_log->SetVisAttributes( G4VisAttributes::Invisible );
	
  //Mylar
  G4VisAttributes *mylar_colour = new G4VisAttributes(G4Colour( 0.5, 0.5, 0.5 ) );
  mylar_wrap_log->SetVisAttributes(G4VisAttributes::Invisible);
  
  //Air
  G4VisAttributes *air_colour = new G4VisAttributes(G4Colour( G4Colour::Blue() ));
  module_log_type1->SetVisAttributes(G4VisAttributes::Invisible);

  //TF1
  G4VisAttributes *TF1_colour = new G4VisAttributes(G4Colour( 0.8, 0.8, 0.0 ) );
  TF1_log->SetVisAttributes(TF1_colour);

  //Steel
  G4VisAttributes *steel_colour = new G4VisAttributes(G4Colour( 0.4, 0.4, 0.4 ) );
  steel_log->SetVisAttributes(steel_colour);
  steel_log_half->SetVisAttributes(steel_colour);

  //PMTcathode
  G4VisAttributes *PMT_colour = new G4VisAttributes(G4Colour( G4Colour::Blue() ));
  PMT_colour->SetForceLineSegmentsPerCircle( 12 );
  PMTcathode_ecal_log->SetVisAttributes(PMT_colour);

  //Electron Arm Package - Houses Polyethylene, CDet, Al
  earm_mother_log->SetVisAttributes(G4VisAttributes::Invisible);

  //Polyethylene
  G4VisAttributes *poly_colour = new G4VisAttributes(G4Colour( 0.2,0.3,0.4 ));
  poly_colour->SetForceWireframe(true);
  polybox_log->SetVisAttributes(poly_colour);

  //CDet
  G4VisAttributes *CD_colour = new G4VisAttributes(G4Colour(0.5, 0.4, 0.1));
  CD_colour->SetForceWireframe(true);
  CD_log->SetVisAttributes(CD_colour);

  //Al
  G4VisAttributes *Al_colour = new G4VisAttributes(G4Colour(0.9, 0.9, 0.9));
  Al_colour->SetForceWireframe(true);
  Al_log->SetVisAttributes(Al_colour);
		    
}
