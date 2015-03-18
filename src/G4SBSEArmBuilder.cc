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

  // BigBite Preshower 

  double psheight = 27*8.5*cm;
  double pswidth  = 2.0*37.0*cm;
  double psdepth  = 8.5*cm;

  G4Box *bbpsbox = new G4Box("bbpsbox", pswidth/2.0, psheight/2.0, psdepth/2.0 );

  G4LogicalVolume* bbpslog = new G4LogicalVolume(bbpsbox, GetMaterial("Air"), "bbpslog");

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth/2.0), bbpslog,
		    "bbpsphys", bbdetLog, false, 0, false);

  // BigBite Shower

  double calheight = 27*8.5*cm;
  double calwidth  = 7*8.5*cm;
  double caldepth  = 37.0*cm;

  G4Box *bbcalbox = new G4Box("bbcalbox", calwidth/2.0, calheight/2.0, caldepth/2.0 );
  G4LogicalVolume* bbcallog = new G4LogicalVolume(bbcalbox, GetMaterial("Lead"), "bbcallog");

  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth+caldepth/2.0+5.0*cm), bbcallog,
		    "bbcalphys", bbdetLog, false, 0, false);

  G4String BBCalSDname = "Earm/BBCal";
  G4String BBCalcolname = "BBCalHitsCollection";
  G4SBSCalSD* BBCalSD;

  if( !(BBCalSD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(BBCalSDname)) ){
    BBCalSD = new G4SBSCalSD( BBCalSDname, BBCalcolname );
    fDetCon->fSDman->AddNewDetector(BBCalSD);
    (fDetCon->SDlist).insert( BBCalSDname );
    fDetCon->SDtype[BBCalSDname] = kCAL;
    //fDetCon->SDarm[BBCalSDname] = kEarm;
  }

  bbcallog->SetSensitiveDetector(BBCalSD);

  // G4cout << "fDetCon->StepLimiterList.size() == " << (fDetCon->StepLimiterList).size() << G4endl;
  // for( set<G4String>::iterator itlist=(fDetCon->StepLimiterList).begin(); itlist != (fDetCon->StepLimiterList).end(); itlist++){
  //   G4cout << "step limiter list element = " << *itlist << G4endl;
  // }
  
  if( (fDetCon->StepLimiterList).find(BBCalSDname) != (fDetCon->StepLimiterList).end() ){
    G4cout << "Creating user limits for " << BBCalSDname << G4endl;
    bbcallog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }
  
  (BBCalSD->detmap).depth=0;
  (BBCalSD->detmap).Row[0] = 0;
  (BBCalSD->detmap).Col[0] = 0;
  (BBCalSD->detmap).LocalCoord[0] = G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth+caldepth/2.0+5.0*cm);
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

  //--------- Visualization attributes -------------------------------
  bbdetLog->SetVisAttributes(G4VisAttributes::Invisible);
  bbfieldLog->SetVisAttributes(G4VisAttributes::Invisible);
  bbmotherLog->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes * yokeVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  //  yokeVisAtt->SetForceWireframe(true);
  bbyokewgapLog->SetVisAttributes(yokeVisAtt);

  /*
    G4VisAttributes * alVisAtt
    = new G4VisAttributes(G4Colour(0.1,0.1,0.1));
    cer_winlog_in->SetVisAttributes(alVisAtt);
    cer_winlog_out->SetVisAttributes(alVisAtt);

    G4VisAttributes * gasVisAtt
    = new G4VisAttributes(G4Colour(0.6,0.6,1.0));
    gasVisAtt->SetForceWireframe(true);
    cer_gaslog->SetVisAttributes(gasVisAtt);
  */

  G4VisAttributes * psVisAtt
    = new G4VisAttributes(G4Colour(0.3,0.9,0.3));
  psVisAtt->SetForceWireframe(true);
  bbpslog->SetVisAttributes(psVisAtt);

  G4VisAttributes * bbcalVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.6,0.0));
  bbcallog->SetVisAttributes(bbcalVisAtt);

}

void G4SBSEArmBuilder::MakeBigCal(G4LogicalVolume *worldlog){

  //-----------------------------
  //  BigCal: Currently just a box with dimensions and sensitivity

  printf("BigCal at %f deg\n", fBBang/deg);

  // Ecal will act as BBcal detector

  double bigcalheight = (24*4.5+32*4.0)*cm;
  double bigcalwidth  = 44.10*2.54*cm;
  double bigcaldepth  = 15.75*2.54*cm;
  //double bbr = fBBdist+bigcaldepth/2.0;
  /*
    printf("BigCal at %f deg\n", fBBang/deg);

    
    double CH2depth = 15.0*cm;
    double CHdepth  = 6.0*cm;

    G4Box *CH2box = new G4Box("ch2box", bigcalwidth/2.0, bigcalheight/2.0, CH2depth/2.0 );
    G4LogicalVolume* ch2boxlog = new G4LogicalVolume(CH2box, GetMaterial("CH2"), "ch2log");
    G4Box *CHbox = new G4Box("chbox", bigcalwidth/2.0, bigcalheight/2.0, CHdepth/2.0 );
    G4LogicalVolume* chboxlog = new G4LogicalVolume(CHbox, GetMaterial("CH"), "chlog");

    double ch2r = bbr - bigcaldepth/2.0 - CHdepth - CH2depth/2.0;
    double chr = bbr - bigcaldepth/2.0 - CHdepth/2.0;


    new G4PVPlacement(bbrm, G4ThreeVector(ch2r*sin(-fBBang), 0.0, ch2r*cos(-fBBang) ), ch2boxlog,
    "ch2boxphys", worldlog, false, 0, false);
    new G4PVPlacement(bbrm, G4ThreeVector(chr*sin(-fBBang), 0.0, chr*cos(-fBBang) ), chboxlog,
    "chboxphys", worldlog, false, 0, false);


    G4VisAttributes * chVisAtt
    = new G4VisAttributes(G4Colour(1.0,0.6,0.0));
    chboxlog->SetVisAttributes(chVisAtt);

    G4VisAttributes * ch2VisAtt
    = new G4VisAttributes(G4Colour(1.0,0.8,0.0));
    ch2boxlog->SetVisAttributes(ch2VisAtt);
  */

  /*************************************************************************************
   ********************************        ECAL      ***********************************
   *************************************************************************************/
   
  //Define  coordinates, rotations, and offsets:
  G4RotationMatrix *bbrm = new G4RotationMatrix;
  bbrm->rotateY(fBBang);

  G4double x_ecal = 246.402*cm, y_ecal = 370.656*cm, z_ecal = 45.406*cm; //45.006cm(module) + 2*0.2(pmt depth)
  double bbr = fBBdist - z_ecal/2.0; //backside of ECal should be located at fBBdist
  double offset = 15*cm; //Motivation - match SBS acceptance

  //Electron Arm - order of materials - all "right next" to each other in z-hat_spectrometer coordinates
  //1)20cm Polyethylene 
  //2)2 4cm CDet planes
  //3)2.4cm Al
  //4)ECal - geometry specified below

  //Dimensions are different than ecal_log because ECal was shrunk significantly when the crescent shape was developed
  //x_earm = 47.5*4.212 (the .5 comes from staggering effects)
  //y_earm = 75*4.212
  //z_earm = 20 + 2*4 + 2.40

  G4double z_Al = 2.40*cm;
  G4double polydepth = 20.00*cm;
  G4double CD_depth = 4.00*cm;
  G4double x_earm = 200.070*cm, y_earm = 315.900*cm, z_earm = 30.40*cm; 

  double cdr = fBBdist - z_ecal - z_earm/2.0;
  G4Box *earm_mother_box = new G4Box("earm_mother_box", x_earm/2.0, y_earm/2.0, z_earm/2.0);
  G4LogicalVolume *earm_mother_log = new G4LogicalVolume(earm_mother_box, GetMaterial("Air"), "earm_mother_log");

  //Polyethylene Box
  G4Box *polybox = new G4Box("polybox", x_earm/2.0, y_earm/2.0, polydepth/2.0 );
  G4LogicalVolume *polybox_log = new G4LogicalVolume( polybox, GetMaterial("Polyethylene"), "polybox_log");

  //Coordinate Detector - 2 planes
  G4Box *CD_box = new G4Box("CD_box", x_earm/2.0, y_earm/2.0, CD_depth/2.0 );
  G4LogicalVolume *CD_log = new G4LogicalVolume(CD_box, GetMaterial("PLASTIC_SC_VINYLTOLUENE"), "CD_log");

  //Aluminum Shielding in front of ECAL - make it big enough to enclose crescent
  G4Box *Al_box = new G4Box( "Al_box", x_earm/2.0, y_earm/2.0, z_Al/2.0 );
  G4LogicalVolume *Al_log = new G4LogicalVolume ( Al_box, GetMaterial("RICHAluminum"), "Al_log" );

  //Place everything inside the Earm Mother Volume
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, z_earm/2.0 - z_Al/2.0 ), Al_log, "Al", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, z_earm/2.0 - z_Al - CD_depth/2.0 ), CD_log, "plane1", earm_mother_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, z_earm/2.0 - z_Al - CD_depth - CD_depth/2.0 ), CD_log, "plane2", earm_mother_log, false, 1 );
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, z_earm/2.0 - z_Al - 2*CD_depth - polydepth/2.0 ), polybox_log, "Polyethylene", earm_mother_log, false, 0 );
  //Note the offsets in x,y,z of 4.212/2.0
  new G4PVPlacement( bbrm, G4ThreeVector( cdr*sin(-fBBang)+(offset+ (4.212/2.0)*cm )*cos(fBBang), -(4.212/2.0)*cm, cdr*cos(-fBBang)+(offset+(4.212/2.0)*cm)*sin(fBBang) ), earm_mother_log, "CDet & Filter", worldlog, false, 0);

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
  G4LogicalVolume *module_log_type1 = new G4LogicalVolume( module_type1, GetMaterial("ECal_Air"), "module_log_type1" );

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

  //Place TF1 mother & Mylar inside module_log_type1 which is already full of ECal_Air
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
  G4SDManager *sdman = fDetCon->fSDman;

  G4String ECalSDname = "Earm/ECal";
  G4String ECalcollname = "ECalHitsCollection";
  G4SBSECalSD *ECalSD = NULL;

  if( !( (G4SBSECalSD*) sdman->FindSensitiveDetector(ECalSDname) ) ){
    G4cout << "Adding ECal sensitive detector to SDman..." << G4endl;
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

