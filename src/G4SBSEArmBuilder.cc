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
  fCerDist = frontGEM_depth - 8.571*cm + 1.811*cm;//this shall be about right
  
  //NB: fBBCalDist now designates the distance to the shielding
  //fix: add an extra 1.136 inch between the back of the GRINCH and the "GEM frame"
  fBBCaldist = fCerDist + fCerDepth + 1.136*2.54*cm + backGEM_depth;
  fGEMDist   = fCerDist + fCerDepth + 1.136*2.54*cm + 0.5*backGEM_depth;
  fGEMOption = 2;
  /**/
  fShieldOption = 1;
  
  fUseLocalField = false;

  fnzsegments_leadglass_ECAL = 1;
  fnzsegments_leadglass_C16 = 1;

  fBuildBBSieve = 0;

  assert(fDetCon);
  
  fDVCSECalMaterial = G4String("PbF2");
  
  //  fbbfield =  NULL;

  fGRINCHgas = "C4F8_gas"; //default to C4F8;
  fTurnOnGrinchPMTglassHits = false;// turn it off by default
  
  fBuildGEMfrontend = false; // do not build it by default  
  fBBPSOption = 2;
}

G4SBSEArmBuilder::~G4SBSEArmBuilder(){;}

void G4SBSEArmBuilder::BuildComponent(G4LogicalVolume *worldlog){
  G4SBS::Exp_t exptype = fDetCon->fExpType;

  //  The neutron experiments and the SIDIS experiment use BigBite:
  //------------ BigBite: -----------------------------------------------------
  if( exptype == G4SBS::kGMN || exptype == G4SBS::kGEN || exptype == G4SBS::kSIDISExp || exptype == G4SBS::kA1n  || exptype == G4SBS::kTDIS || exptype == G4SBS::kGEnRP ) 
    {
      MakeBigBite( worldlog );
      //Move sieve slit construction to MakeBigBite subroutine:
      //      if(fBuildBBSieve)
      //	MakeBBSieveSlit(worldlog);
    }
  if( exptype == G4SBS::kGEp || exptype == G4SBS::kGEPpositron ) //Subsystems unique to the GEp experiment include FPP and BigCal:
    {
      G4SBSECal* ECal = new G4SBSECal(fDetCon);
      ECal->SetAng(fBBang);
      ECal->SetDist(fBBdist);
      ECal->BuildComponent(worldlog);
      //MakeBigCal( worldlog );
    }
  if( exptype == G4SBS::kC16 ) 
    {
      G4SBSECal* ECal = new G4SBSECal(fDetCon);
      ECal->SetAng(fBBang);
      ECal->SetDist(fBBdist);
      ECal->BuildComponent(worldlog);
      //MakeC16( worldlog );
    }
  if( (exptype == G4SBS::kGMN || exptype == G4SBS::kGEnRP) && fBuildGEMfrontend )  MakeGMnGEMShielding_update( worldlog );
  
  if( exptype == G4SBS::kNDVCS ){
    MakeDVCSECal(worldlog);
  }
  
  if( exptype ==  G4SBS::kGEMHCtest){
    MakeHallCGEM(worldlog);
  }
}

void G4SBSEArmBuilder::MakeBigBite(G4LogicalVolume *worldlog){
  bool chkoverlap = false;
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
  G4SubtractionSolid* bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxL", bbmotherBox, bbleftcutTrap, leftcutrot2, G4ThreeVector(-10*eps, 0.0, -motherdepth/2.0-1.0*m+clear));
  //bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR", bbmothercutBox, bbrightcutTrap, rightcutrot2, G4ThreeVector(10*eps, 0.0, -motherdepth/2.0+clear));
  //EPAF: 2017/04/17: commented this line to avoid to have some of the GRINCH PMTs outside of the mother volume.
  //Besides, it is not necessary, as it removes some volume from the mother box which would not interfere with anything anyhow
  
  G4Box *frontboxcut = new G4Box("frontboxcut",(bbmagwidth-gapsize)/4.0-coilwidth, 250*cm + eps, clear/2);
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR_fL", bbmothercutBox, frontboxcut, 0, G4ThreeVector( ((bbmagwidth+gapsize)/4.0+coilwidth)+ 10*eps +5*cm, 0, -motherdepth/2.0+clear/2 - 10*eps ));
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR_fLR", bbmothercutBox, frontboxcut, 0, G4ThreeVector( -((bbmagwidth+gapsize)/4.0+coilwidth)-10*eps -5*cm, 0, -motherdepth/2.0+clear/2 - 10*eps ));

  G4Box *bottomboxcut = new G4Box("bottomboxcut",bbmagwidth+eps, 250*cm/2, motherdepth+eps);
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR_fLR_floor", bbmothercutBox, bottomboxcut, 0, G4ThreeVector( 0.0, -1.25*m*2-20*cm, 0.0));
  
  // EPAF: 2017/05/15: Added this cut to leave some room for the lead shielding for GMn.
  // I checked it was not interfering with the BB magent volume nor with the BB field volume...
  G4Box *beamshieldcut = new G4Box("beamshieldcut", 30.0*cm/2.0 + eps, 300*cm + eps, 10.0*cm/2.0 + eps);
  bbmothercutBox = new G4SubtractionSolid("bbmothercutBoxLR_fLR_shieldcut", bbmothercutBox, beamshieldcut, 0, G4ThreeVector( -(bbmagwidth+gapsize)/4.0+coilwidth+15.0*cm, 0.0, -motherdepth/2.0+2.5*cm));
					  

  //   Make logical volume for mother vol and place
  G4LogicalVolume *bbmotherLog=new G4LogicalVolume(bbmothercutBox,GetMaterial("Air"),
						   "bbmotherLog", 0, 0, 0);

  
  new G4PVPlacement(bbrm, G4ThreeVector(motherr*sin(fBBang), 0.0, motherr*cos(fBBang)),
		    bbmotherLog, "bbmotherPhys", worldlog, false,0, chkoverlap);


  new G4PVPlacement(yokerm,G4ThreeVector(0.0, 0.0, -motherdepth/2.0+clear),
		    bbyokewgapLog, "bbyokewgapPhysical", bbmotherLog, false,0, chkoverlap);


  
  //Sieve plate position is 13.37 inches ~= 34 cm upstream of front of magnet yoke:
  //Thickness of sieve plate is 1.5 inches
  //SSeeds - added second sieve option to accomodate new plate 10.4.20
  //sseeds - added third and fourth for optics studies 11.5.20. Will need to find more efficient way to implement
  
 if( fBuildBBSieve == 1 ){ 
    G4ThreeVector BBSievePos(0,0,-motherdepth/2.0+36.0*cm-0.75*2.54*cm);
    MakeBBSieveSlit( bbmotherLog, BBSievePos );
 }
 else if( fBuildBBSieve == 2 ){  
   G4ThreeVector BBSievePos(0,0,-motherdepth/2.0+36.0*cm-0.75*2.54*cm); //Not sure where 0.75" comes from - sseeds
    MakeNewBBSieveSlit( bbmotherLog, BBSievePos );
 }
 else if( fBuildBBSieve == 3 ){  
   G4ThreeVector BBSievePos(0,0,-motherdepth/2.0+36.0*cm-0.75*2.54*cm); //Not sure where 0.75" comes from - sseeds
   MakeThirdBBSieveSlit( bbmotherLog, BBSievePos );
 }
 else if( fBuildBBSieve == 4 ){  
   G4ThreeVector BBSievePos(0,0,-motherdepth/2.0+36.0*cm-0.75*2.54*cm); //Not sure where 0.75" comes from - sseeds
   MakeFourthBBSieveSlit( bbmotherLog, BBSievePos );
 }
 else {
   //cout << "Invalid sieve entry. Please enter 1 (old - straight holes and slots) or 2 (new - holes with angles in dispersive direction) or 3 (new - holes without angles). No sieve constructed.\n";
 }
  
  //  Bigbite field log volume
  G4LogicalVolume *bbfieldLog=new G4LogicalVolume(bbairTrap, GetMaterial("Air"),
						  "bbfieldLog", 0, 0, 0);
  
  //NOTE that the invocation of the commented out command below can potentially cause
  // G4SBSBigBiteField::ReadField() to
  // be invoked twice in the same run. Not a big deal, but also not desirable behavior
  // NOW BB field is only created ONCE when G4SBSDetectorConstruction::SetBigBiteField() is invoked.
  // fbbfield = new G4SBSBigBiteField( G4ThreeVector(0.0, 0.0, fBBdist), *bbrm, bbfieldmap_fname );

  if( fUseLocalField ){
    G4SBSBigBiteField *bbfield = (G4SBSBigBiteField*) fDetCon->GetBBField();

    if( bbfield ){
  
      bbfield->fInverted = fDetCon->fGlobalField->fInverted;
      bbfield->fScaleFactor = fDetCon->GetFieldScale_BB();
  

      G4FieldManager *bbfm = new G4FieldManager(bbfield);
      //new G4ChordFinder(fbbfield);
      bbfm->SetDetectorField(bbfield);
      bbfm->CreateChordFinder(bbfield);
  
  
      bbmotherLog->SetFieldManager(bbfm,true);
    }
  }

  //does this volume serve any purpose? Apparently not
  new G4PVPlacement(0, G4ThreeVector(), bbfieldLog, "bbfieldPhysical", bbyokewgapLog, false,0, chkoverlap);

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
		    bbdetLog, "bbdetPhysical", bbmotherLog, false,0, chkoverlap);

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
  double gemw_opt2[] = { 40.96*cm, 40.96*cm, 40.96*cm, 40.96*cm, 61.44*cm };
  double gemh_opt2[] = { 153.6*cm, 153.6*cm, 153.6*cm, 153.6*cm, 204.8*cm };

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
	G4cout << "GEM " << i << " (mm) " <<detoffset+gemz_opt2[i]+detboxdepth/2.0 << endl;
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

  //how are these used?
  double mylarthickness = 0.0020*cm, airthickness = 0.0040*cm;
  double mylar_air_sum = mylarthickness + airthickness;
  double bbpmtz = 0.20*cm;
  
  G4cout << " BBPS option (0: old geometry; 1: new modules, 25 blocks; 2: new modules, 26 blocks)" << fBBPSOption << endl;
  
  // **** BIGBITE CALORIMETER MOTHER VOLUME ****:
  // EPAF: 2020/02/28: updateing the detector package geometry to account for the last developments:
  // - 1/4" steel plate;
  // - PS;
  // - 1/4" steel plate;
  // - 2mm Al plate;
  // - hodoscope;
  // - honeycomb;
  // - shower;
  
  //AJRP 05/19/19: re-working BBCAL geometry so user tracking action and stepping action classes
  //behave properly. Put everything (including detectors/shielding) inside a single mother volume
  G4double bbcal_box_height = 27*8.5*cm;
  if(fBBPSOption==2)bbcal_box_height = 26*9.0*cm;
  G4double bbcal_box_width  = 2.0*37.0*cm;
  //G4double bbcal_box_depth  = (8.5+2.5+37.0)*cm;
  // space for optional shielding + steel plate + mu-metal + PS block + space for hodoscope + Shower blocks
  G4double shielding_space = 3.75*2.54*cm;
  G4double bbcal_box_depth  = shielding_space+(0.25*2.54 + 0.05 + 8.5 + 3.5*2.54 +37.0)*cm;//8.89 cm (3.5") is the size of the gap between the PS and the SH

  // Big Bite Calorimeter shielding.
  // 
  // EPAF: 2017/03/02: flag for BBECal shielding option: 
  // 0: nothing; 1: default 1/4 in SS+0.5mm; 2: 10cm Al + 3cm SS on the side; 3: 10cm Al + 3cm SS on the side; 
  
  //AJRP: cleaning up shielding: total thickness of cal mother box will be the total thickness of all BB calorimeter
  //layers (including shielding)
  //NOTE: "SIDE" shielding never actually gets placed as coded; thus, ignore:
  G4double bbcal_total_thick=bbcal_box_depth;
  G4double Al_thick = 9.5*cm;
  G4double SS_thick = 2.0*cm;
  // Default front plate: 0.25" steel + 0.5mm mu metal
  //NO: the space for the shielding is available no matter what!!!
  /*
  if( fShieldOption > 0 ){
    //bbcal_total_thick += 0.25*2.54*cm + 0.5*mm; //Steel + mu-metal
    switch( fShieldOption ){
    case 2: //"10 cm Al in front "
      //bbcal_total_thick += Al_thick;
      break;
    case 3: //"2 cm SS in front "
      //bbcal_total_thick += SS_thick;
      break;
    case 4: //"1 cm SS in front + 5cm Al in front "
      Al_thick /= 2.0;
      SS_thick /= 2.0;
      //bbcal_total_thick += Al_thick + SS_thick; 
      break;
    default:
      break;
    }
  }
  */
  if( fShieldOption==4 ){
    Al_thick /= 2.0;
    SS_thick /= 2.0;
  }
  
  //  G4double bbcal_shield_thick = 6.85*mm + 9.525*cm;
  // G4double Al_thick = 10.0*cm;
  // G4double SS_thick = 2.0*cm;
  // if(fShieldOption==2)bbcal_shield_thick+=max(0.0, Al_thick-9.0*cm);
  // if(fShieldOption==4){
  //   Al_thick = Al_thick/2.0;
  //   SS_thick = SS_thick/2.0;
  // }

  // BB Ecal
  G4Box *bbcalbox = new G4Box( "bbcalbox", bbcal_box_width/2.0, bbcal_box_height/2.0, bbcal_total_thick/2.0+0.1*mm );
  G4LogicalVolume *bbcal_mother_log = new G4LogicalVolume(bbcalbox, GetMaterial("Air"), "bbcal_mother_log");
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, detoffset + fBBCaldist + bbcal_total_thick/2.0 ), bbcal_mother_log, "bbcal_mother_phys", bbdetLog, false, 0 , chkoverlap); 

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
      fDetCon->SDtype[SDname] = G4SBS::kCAL;
      (BBCalSD->detmap).depth = 0;
      bbcal_mother_log->SetSensitiveDetector( BBCalSD );
      //(BBCalSD->detmap).Row[0] = 0;
      //(BBCalSD->detmap).Col[0] = 0;
      
    }
  }

  //Get rid of dedicated shielding box, it's not needed
  // G4Box *bbcalshieldbox = new G4Box( "bbcalshieldbox", bbmagwidth/2.0-2.0*cm, bbcal_box_height/2.0, bbcal_shield_thick/2.0 );
  // G4LogicalVolume *bbcal_shield_log = new G4LogicalVolume(bbcalshieldbox, GetMaterial("Air"), "bbcal_shield_log");
  // bbcal_shield_log->SetVisAttributes( G4VisAttributes::Invisible );
  
  // G4Box *bbcalfrontmufoil = new G4Box( "bbcalfrontmufoil", bbcal_box_width/2.0, bbcal_box_height/2.0, 0.5*mm/2.0 );
  // G4LogicalVolume *bbcal_front_mufoil_log = new G4LogicalVolume(bbcalfrontmufoil, GetMaterial("mu-metal"), "bbcal_front_mufoil_log");
  // //bbcal_front_mufoil_log->SetVisAttributes( G4Colour(0.,1.0, 0.0) );
  
  G4Box *hodosupportAlplate = new G4Box( "hodosupportAlplate", bbcal_box_width/2.0, bbcal_box_height/2.0, 2.0*mm/2.0 );
  G4LogicalVolume *hodosupportAlplate_log = new G4LogicalVolume(hodosupportAlplate, GetMaterial("Aluminum"), "hodosupportAlplate_log");
  hodosupportAlplate_log->SetVisAttributes( G4Colour(0.8,0.8, 0.8) );

  G4Box *honeycombAlplate = new G4Box( "honeycombAlplate", bbcal_box_width/2.0, bbcal_box_height/2.0, 0.5*mm/2.0 );
  G4LogicalVolume *honeycombAlplate_log = new G4LogicalVolume(honeycombAlplate, GetMaterial("Aluminum"), "honeycombAlplate_log");
  honeycombAlplate_log->SetVisAttributes( G4Colour(0.8,0.8, 0.8) );
  
  G4Box *bbcalfrontsteelplate = new G4Box( "bbcalfrontsteelplate", bbcal_box_width/2.0, bbcal_box_height/2.0, 6.35*mm/2.0 );
  G4LogicalVolume *bbcal_front_steelplate_log = new G4LogicalVolume(bbcalfrontsteelplate, GetMaterial("Steel"), "bbcal_front_steelplate_log");
  bbcal_front_steelplate_log->SetVisAttributes( G4Colour(0.0, 0.0, 1.0) );
  
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

  //start at the front of the BB cal mother volume and work our way back:
  G4double zpos_temp = -0.5*bbcal_total_thick;
  if(fShieldOption > 0){
    // new G4PVPlacement( 0, G4ThreeVector( 0, 0, detoffset + fBBCaldist + bbcal_shield_thick/2.0 ), bbcal_shield_log, "bbcal_shield_phys", bbdetLog, false, 0 );
    
    switch(fShieldOption){
    case 2: //aluminum shielding in front ONLY (10 cm)
      zpos_temp += bbcalshield_al->GetZHalfLength();
      new G4PVPlacement( 0, G4ThreeVector( 0, 0, zpos_temp), bbcal_shield_al_log, "bbcal_shield_al_phys", bbcal_mother_log, false, 0, chkoverlap );
      zpos_temp += bbcalshield_al->GetZHalfLength();
      //new G4PVPlacement( 0, G4ThreeVector( (-bbmagwidth+3.0*cm)/2.0, 0, -detboxdepth/4.0), bbcal_shield_side_ss_log, "bbcal_shield_side_ss_phys", bbdetLog, false, 0 );
      break;
    case 3: //mutually exclusive with case 2: SS shielding in front ONLY (2 cm)
      zpos_temp += bbcalshield_ss->GetZHalfLength();
      new G4PVPlacement( 0, G4ThreeVector( 0, 0, zpos_temp), bbcal_shield_ss_log, "bbcal_shield_ss_phys", bbcal_mother_log, false, 0, chkoverlap );
      zpos_temp += bbcalshield_ss->GetZHalfLength();
      //new G4PVPlacement( 0, G4ThreeVector( (-bbmagwidth+3.0*cm)/2.0, 0, -detboxdepth/4.0), bbcal_shield_side_ss_log, "bbcal_shield_side_ss_phys", bbdetLog, false, 0 );2.-
      break;
    case 4: //SS + Al:
      zpos_temp += bbcalshield_ss->GetZHalfLength();
      new G4PVPlacement( 0, G4ThreeVector( 0, 0, zpos_temp), bbcal_shield_ss_log, "bbcal_shield_ss_phys", bbcal_mother_log, false, 0, chkoverlap );
      zpos_temp += bbcalshield_ss->GetZHalfLength() + bbcalshield_al->GetZHalfLength();
      new G4PVPlacement( 0, G4ThreeVector( 0, 0, zpos_temp), bbcal_shield_al_log, "bbcal_shield_al_phys", bbcal_mother_log, false, 0, chkoverlap );
      //new G4PVPlacement( 0, G4ThreeVector( (-bbmagwidth+3.0*cm)/2.0, 0, -detboxdepth/4.0), bbcal_shield_side_ss_log, "bbcal_shield_side_ss_phys", bbdetLog, false, 0 );
      zpos_temp += bbcalshield_al->GetZHalfLength();
    default:  //do nothing here; no "extra" shielding:
      break;
    }
    zpos_temp = -0.5*bbcal_total_thick+shielding_space;
    //Place front mu-metal foil and front steel plate?
    //EPAF:2020/02/28: remove all mu metal foils
    //zpos_temp += bbcalfrontmufoil->GetZHalfLength();
    //new G4PVPlacement( 0, G4ThreeVector( 0, 0, zpos_temp), bbcal_front_mufoil_log, "bbcal_front_mufoil_phys", bbcal_mother_log, false, 0, chkoverlap );
    zpos_temp += //bbcalfrontmufoil->GetZHalfLength() + 
      bbcalfrontsteelplate->GetZHalfLength();
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, zpos_temp ), bbcal_front_steelplate_log, "bbcal_front_steelplate_phys", bbcal_mother_log, false, 0, chkoverlap );
    zpos_temp += bbcalfrontsteelplate->GetZHalfLength();
  }
  
  //At this stage, zpos_temp = -bbcal_total_thick/2 + total shielding thickness! 
  // NO: zpos_temp = -bbcal_total_thick/2+shielding_space+ front mu-metal foil and front steel plate thickness

  // **** BIGBITE PRESHOWER **** 
  // 2 columns, 27 rows

  double psheight = 27*8.5*cm;
  double pswidth  = 2.0*37.0*cm;
  double psdepth  = 8.5*cm;

  if(fBBPSOption>=1){
    psdepth = 9.0*cm;
    psheight = psdepth*25;
    if(fBBPSOption==2)psheight+=psdepth;
  }

  G4Box *bbpsbox = new G4Box("bbpsbox", pswidth/2.0, psheight/2.0, psdepth/2.0 );
  G4LogicalVolume *bbpslog = new G4LogicalVolume(bbpsbox, GetMaterial("Air"), "bbpslog");
  //new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fBBCaldist+psdepth/2.0), bbpslog, "bbpsphys", bbdetLog, false, 0);
  G4cout << "bb ps front (mm) "<< detoffset + fBBCaldist+bbcal_total_thick/2.+zpos_temp+detboxdepth/2.0 << endl;
  zpos_temp += psdepth/2.0;
  new G4PVPlacement(0, G4ThreeVector( 0, 0, zpos_temp ), bbpslog, "bbpsphys", bbcal_mother_log, false, 0, chkoverlap );
  zpos_temp += psdepth/2.0 + bbcalfrontsteelplate->GetZHalfLength();
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, zpos_temp ), bbcal_front_steelplate_log, "bbcal_front_steelplate_phys", bbcal_mother_log, false, 0, chkoverlap );
  zpos_temp += bbcalfrontsteelplate->GetZHalfLength();
  
  
  //placement of second mu-metal foil behind the PS
  //zpos_temp += psdepth/2.0 + bbcalfrontmufoil->GetZHalfLength();
  //new G4PVPlacement( 0, G4ThreeVector( 0, 0, zpos_temp), bbcal_front_mufoil_log, "bbcal_back_mufoil_phys", bbcal_mother_log, false, 0, chkoverlap );
  // Preshower module - geometry will be assigned after Shower
  
  //Aluminum foil for Al "honeycomb"
  zpos_temp += hodosupportAlplate->GetZHalfLength()/2.0+2.54*cm;
  new G4PVPlacement( 0, G4ThreeVector(0,0, zpos_temp ), hodosupportAlplate_log, "hodosupportAlplate_phys", bbcal_mother_log, false, 0, chkoverlap );
  zpos_temp += hodosupportAlplate->GetZHalfLength()/2.0;
  //zpos_temp += 1.905*cm-2*bbcalfrontmufoil->GetZHalfLength();
  //new G4PVPlacement( 0, G4ThreeVector(0,0, zpos_temp ), honeycombAlplate_log, "honeycombplatefront2_phys", bbcal_mother_log, false, 0, chkoverlap );
  
  // **** BIGBITE HODOSCOPE **** 
  // Scintillator box - same dimensions as preshower
  G4int n_bbhodoslats = 90;
  double bbslat_length = 60.0*cm;
  double bbslat_section = 2.5*cm;
  // logic volume dimensions...
  double bbhododepth = bbslat_section;
  double bbhodoheight= bbslat_section*n_bbhodoslats;
  double bbhodowidth = bbslat_length;
  
  G4Box *bbhodobox = new G4Box("bbhodobox", pswidth/2.0, psheight/2.0, bbhododepth/2.0 );
  G4LogicalVolume *bbhodolog = new G4LogicalVolume( bbhodobox, GetMaterial("Air"), "bbhodolog" );
  //new G4PVPlacement(0, G4ThreeVector(0.0,0.0, detoffset+fBBCaldist+psdepth+bbhododepth/2.0), bbhodolog, "bbhodophys", bbdetLog, false, 0);
  //new G4PVPlacement( 0, G4ThreeVector(0,0, -bbcal_box_depth/2.0 + psdepth + bbhododepth/2.0 ), bbhodolog, "bbhodophys", bbcal_mother_log, false, 0 );
  //new G4PVPlacement( 0, G4ThreeVector(0,0, -bbcal_box_depth/2.0 + psdepth + 0.217*2.54 + bbhododepth/2.0 ), bbhodolog, "bbhodophys", bbcal_mother_log, false, 0 );
  G4cout << "BB hodo front (mm) " << detoffset + fBBCaldist+bbcal_total_thick/2.+zpos_temp+detboxdepth/2.0+1.0*mm << endl;
  zpos_temp += bbhododepth/2.0+1.0*mm;
  new G4PVPlacement( 0, G4ThreeVector(0,0, zpos_temp ), bbhodolog, "bbhodophys", bbcal_mother_log, false, 0, chkoverlap );
  bbhodolog->SetVisAttributes(G4VisAttributes::Invisible);
  // zpos_temp += bbhododepth/2.0+honeycombAlplate->GetZHalfLength();
  // new G4PVPlacement( 0, G4ThreeVector(0,0, zpos_temp ), honeycombAlplate_log, "honeycombplateback1_phys", bbcal_mother_log, false, 0, chkoverlap );
  // zpos_temp += 1.905*cm-2*bbcalfrontmufoil->GetZHalfLength();
  // new G4PVPlacement( 0, G4ThreeVector(0,0, zpos_temp ), honeycombAlplate_log, "honeycombplateback2_phys", bbcal_mother_log, false, 0, chkoverlap );
 
  //zpos_temp +=3.175*cm;
  
  //
  G4Box *bbhodoslatbox = new G4Box("bbhodoslatbox", bbslat_length/2.0, bbslat_section/2.0, bbslat_section/2.0);
  /*
  G4Box* mylarboxhollow = new G4Box("mylarboxhollow", bbslat_length/2.0, bbslat_section/2.0-mylarthickness, bbslat_section/2.0-mylarthickness);
  G4SubtractionSolid* mylarwrap = new G4SubtractionSolid("mylarwrap", VetoElemBox, MylarBoxHollow, 0, G4ThreeVector(0.0, 0.0, 0.0));
  
  G4LogicalVolume *mylarwrap_log = new G4LogicalVolume( MylarWrap, GetMaterial("Mylar"), "mylarwrap_log" );
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), MylarWrapLog, "mylarwrap_phys", VetoElemLog, false, 0, checkOL);
   
  G4Box *bbhodoscintbox = new G4Box("bbhodoscint", bbslat_length/2.0, bbslat_section/2.0-mylar_air_sum, bbslat_section/2.0-mylar_air_sum);
  */
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
    fDetCon->SDtype[BBHodoScintSDname] = G4SBS::kCAL;
    (BBHodoScintSD->detmap).depth = 0;

    G4double ethresh_default = 0.0*MeV;
    G4double timewindow_default = 30.0*ns;
    
    fDetCon->SetTimeWindowAndThreshold( BBHodoScintSDname, ethresh_default, timewindow_default );
  }
  bbhodoslatlog->SetSensitiveDetector( BBHodoScintSD ); 

  //Record track info at both the calorimeter mother volume boundary AND the hodoscope mother volume boundary:
  fDetCon->InsertSDboundaryVolume( bbcal_mother_log->GetName(), BBHodoScintSDname );
  //fDetCon->InsertSDboundaryVolume( bbhodolog->GetName(),       BBHodoScintSDname );
  ofstream mapfile("database/BBhodo_map.txt");

  TString currentline;

  currentline.Form("# %15s, %15s, %15s, %18s, %18s",
		   "Cell", "Row", "Column", "Xcenter", "Ycenter" );
  
  for(int i_bbhslat = 0; i_bbhslat<n_bbhodoslats; i_bbhslat++){
    G4double y_slat = n_bbhodoslats*bbslat_section/2.0-(G4double(i_bbhslat)+0.5)*bbslat_section;
    //G4double z_slat = -bbhododepth/2.0-0.5*mm+0.217*2.54*cm+bbslat_section/2.0;
    G4double z_slat = 0.0;// in the middle -sound about right provided it is sandwiched between 2 honeycombs - which will add in later
    new G4PVPlacement( 0, G4ThreeVector(0, y_slat, z_slat), bbhodoslatlog, "bbhodoslatphys", bbhodolog, false, i_bbhslat, chkoverlap );
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
  
  zpos_temp = bbcal_total_thick/2.0 - caldepth - 1.905*cm+honeycombAlplate->GetZHalfLength()/2.0;
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, zpos_temp), honeycombAlplate_log, "honeycombplate1_phys", bbcal_mother_log, false, 0, chkoverlap );
  zpos_temp+=1.905*cm-honeycombAlplate->GetZHalfLength();
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, +bbcal_total_thick/2.0 - caldepth + honeycombAlplate->GetZHalfLength()/2.0-1.905*cm), honeycombAlplate_log, "honeycombplate2_phys", bbcal_mother_log, false, 0, chkoverlap );
  zpos_temp+=honeycombAlplate->GetZHalfLength()/2.0;

  G4cout << "BB SH front (m) " << detoffset + fBBCaldist+bbcal_total_thick-caldepth+detboxdepth/2.0 << endl;
  
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, +bbcal_total_thick/2.0 - caldepth/2.0), bbshowerlog, "bbshowerphys", bbcal_mother_log, false, 0, chkoverlap );
  

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
    fDetCon->SDtype[BBSHTF1SDname] = G4SBS::kCAL;
    (BBSHTF1SD->detmap).depth = 1;

    G4double threshold_default = 0.0*MeV; //1% of 1 GeV
    G4double timewindow_default = 50.0*ns; //We could use 10 ns here if we wanted, but also have to consider pulse shape. 
    
    fDetCon->SetTimeWindowAndThreshold( BBSHTF1SDname, threshold_default, timewindow_default );
  }
  bbTF1log->SetSensitiveDetector( BBSHTF1SD );

  //fDetCon->InsertSDboundaryVolume( bbshowerlog->GetName(), BBSHTF1SDname );
  fDetCon->InsertSDboundaryVolume( bbcal_mother_log->GetName(), BBSHTF1SDname );
  
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
    fDetCon->SDtype[BBSHSDname] = G4SBS::kECAL;
    (BBSHSD->detmap).depth = 1;
  }
  bbpmtcathodelog->SetSensitiveDetector( BBSHSD );

  fDetCon->InsertSDboundaryVolume( bbcal_mother_log->GetName(), BBSHSDname );
  
  // Put everything in a BB Shower Module
  int shower_copy_number = 0;
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-bbpmtz)/2.0), bbpmtcathodelog,"bbcathodephys", showermodlog, false, 0, chkoverlap );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-3*bbpmtz)/2.0), bbpmtwindowlog, "bbwindowphys", showermodlog, false, 0, chkoverlap );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-4*bbpmtz-bbTF1_z)/2.0), bbTF1log, "bbTF1phys", showermodlog, false, 0, chkoverlap );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, -bbpmtz), bbmylarwraplog, "bbmylarphys", showermodlog, false, 0, chkoverlap );

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
      
      new G4PVPlacement(0, G4ThreeVector(xtemp,ytemp,0.0), showermodlog, "showermodphys", bbshowerlog, false, shower_copy_number, chkoverlap);
      
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
  
  if(fBBPSOption>0){
    bbmodule_x = bbmodule_y = 9.0*cm;
  }
  // ****Preshower Continued****
  // Reusing modules from Shower (same variables), rotated by either +/- 90 deg depending on column #
  G4Box *preshowermodbox = new G4Box( "preshowermodbox", bbmodule_x/2.0, bbmodule_y/2.0, caldepth/2.0 );
  G4LogicalVolume *preshowermodlog = new G4LogicalVolume( preshowermodbox, GetMaterial("Special_Air"), "preshowermodlog" );
 
  // Preshower TF1 SD of type CAL
  G4LogicalVolume *bbpsTF1log;
  if(fBBPSOption==0){
    bbpsTF1log = new G4LogicalVolume( bbTF1box, GetMaterial("TF5"), "bbpsTF1log" );
  }else{
    bbpsTF1log = new G4LogicalVolume( bbTF1box, GetMaterial("TF1"), "bbpsTF1log" );
  }
  G4cout << "BBPS blocks material is " << bbpsTF1log->GetMaterial()->GetName() << G4endl;

  G4String BBPSTF1SDname = "Earm/BBPSTF1";
  G4String BBPSTF1collname = "BBPSTF1HitsCollection";
  G4SBSCalSD *BBPSTF1SD = NULL;

  if( !((G4SBSCalSD*) sdman->FindSensitiveDetector(BBPSTF1SDname)) ) {
    G4cout << "Adding BB Preshower TF1 Sensitive Detector to SDman..." << G4endl;
    BBPSTF1SD = new G4SBSCalSD( BBPSTF1SDname, BBPSTF1collname );
    
    sdman->AddNewDetector( BBPSTF1SD );
    (fDetCon->SDlist).insert( BBPSTF1SDname );
    fDetCon->SDtype[BBPSTF1SDname] = G4SBS::kCAL;
    (BBPSTF1SD->detmap).depth = 1;

    //Photoelectron yield is approximately 500/GeV (or so)
    G4double threshold_default = 0.0*MeV; //1% of 1 GeV
    G4double timewindow_default = 50.0*ns; //We could use 10 ns here if we wanted, but also have to consider pulse shape.

    fDetCon->SetTimeWindowAndThreshold( BBPSTF1SDname, threshold_default, timewindow_default );
  }
  bbpsTF1log->SetSensitiveDetector( BBPSTF1SD );

  fDetCon->InsertSDboundaryVolume( bbcal_mother_log->GetName(), BBPSTF1SDname );

  if( (fDetCon->StepLimiterList).find( BBPSTF1SDname ) != (fDetCon->StepLimiterList).end() ){
    bbpsTF1log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }

  // Preshower PMT SD of type ECAL
  //G4LogicalVolume *bbpspmtcathodelog = new G4LogicalVolume( bbPMT, GetMaterial("Photocathode_material_ecal"), "bbpspmtcathodelog" );

  G4LogicalVolume *bbpspmtcathodelog = new G4LogicalVolume( bbPMT, GetMaterial("Photocathode_BB"), "bbpspmtcathodelog" );
  
  G4String BBPSSDname = "Earm/BBPS";
  G4String BBPScollname = "BBPSHitsCollection";
  G4SBSECalSD *BBPSSD = NULL;

  if( !((G4SBSECalSD*) sdman->FindSensitiveDetector(BBPSSDname)) ) {
    G4cout << "Adding BB Preshower PMT Sensitive Detector to SDman..." << G4endl;
    BBPSSD = new G4SBSECalSD( BBPSSDname, BBPScollname );
    sdman->AddNewDetector( BBPSSD );
    (fDetCon->SDlist).insert(BBPSSDname);
    fDetCon->SDtype[BBPSSDname] = G4SBS::kECAL;
    (BBPSSD->detmap).depth = 1;
  }
  bbpspmtcathodelog->SetSensitiveDetector( BBPSSD );

  fDetCon->InsertSDboundaryVolume( bbcal_mother_log->GetName(), BBPSSDname );
  
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-bbpmtz)/2.0), bbpspmtcathodelog,"bbpscathodephys", preshowermodlog, false, 0, chkoverlap );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-3*bbpmtz)/2.0), bbpmtwindowlog, "bbpswindowphys", preshowermodlog, false, 0, chkoverlap );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-4*bbpmtz-bbTF1_z)/2.0), bbpsTF1log, "bbpsTF1phys", preshowermodlog, false, 0, chkoverlap );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, -bbpmtz), bbmylarwraplog, "bbpsmylarphys", preshowermodlog, false, 0, chkoverlap );
  
  G4RotationMatrix *bbpsrm_col1 = new G4RotationMatrix;
  bbpsrm_col1->rotateY(-90.0*deg);
  G4RotationMatrix *bbpsrm_col2 = new G4RotationMatrix;
  bbpsrm_col2->rotateY(90.0*deg);
  
  int bbpscol = 2;
  int bbpsrow = 27;
  if(fBBPSOption==1)bbpsrow = 25;
  if(fBBPSOption==2)bbpsrow = 26;
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
	new G4PVPlacement( bbpsrm_col1, G4ThreeVector(xtemp,ytemp,0.0), preshowermodlog, "preshowermodphys", bbpslog, false, ps_copy_number, chkoverlap );
	(BBPSSD->detmap).LocalCoord[ps_copy_number] = G4ThreeVector(xtemp+caldepth/2.0-bbpmtz/2.0, ytemp, 0.0);
	(BBPSTF1SD->detmap).LocalCoord[ps_copy_number] = G4ThreeVector(xtemp,ytemp,0.0);
	ps_copy_number++;
      }
      if(l==1) {
	new G4PVPlacement( bbpsrm_col2, G4ThreeVector(xtemp,ytemp,0.0), preshowermodlog, "preshowermodphys", bbpslog, false, ps_copy_number, chkoverlap );
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
  //   fDetCon->SDtype[BBCalSDname] = G4SBS::kCAL;
  //   //fDetCon->SDarm[BBCalSDname] = G4SBS::kEarm;
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
  grinch->SetCerDepth( fCerDepth );
  grinch->SetGrinchGas( fGRINCHgas );
  grinch->SetTurnOnPMTglassHits( fTurnOnGrinchPMTglassHits );
  grinch->BuildComponent( bbdetLog );
  
  //Shielding for UVA GEM
  G4Box *Shield_backgem_box = new G4Box("Shield_backgem_box", 0.5*65*cm, 0.5*210*cm, 0.5*2.54*cm);
  G4LogicalVolume *Shield_backgem_log = new G4LogicalVolume(Shield_backgem_box, GetMaterial("CH2"), "Shield_backgem_log");//GetMaterial("CDET_Acrylic") ???
  //new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fCerDist+fCerDepth+0.51*2.54*cm ), Shield_backgem_log, "", bbdetLog, false, 0, true);
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, detoffset+fGEMDist - 0.51*2.54*cm -1.6*cm ), Shield_backgem_log, "", bbdetLog, false, 0, true);
}

void G4SBSEArmBuilder::MakeDVCSECal(G4LogicalVolume *motherlog){
  G4cout << "Building DVCS ECal with following material: " << fDVCSECalMaterial << endl;
  G4bool defined_mat = false;
  G4int Nrows, Ncols;
  G4double dvcsblkmodule_x, dvcsblkmodule_y;
  G4double caldepth;
  G4double EcalHoffset, EcalVoffset;
  if(fDVCSECalMaterial=="PbF2"){
    Nrows = 16;
    Ncols = 13;
    dvcsblkmodule_x = 3.00*cm;
    dvcsblkmodule_y = 3.00*cm;
    caldepth = 18.6*cm+2*2.0*cm;
    EcalHoffset = +1.5*3.0*cm;
    EcalVoffset = 0;
    defined_mat = true;
  }
  if(fDVCSECalMaterial=="PbWO4"){
    Nrows = 31;
    Ncols = 36;
    dvcsblkmodule_x = 2.05*cm;
    dvcsblkmodule_y = 2.05*cm; 
    caldepth = 18.0*cm+2*2.0*cm;
    EcalHoffset = 0;
    EcalVoffset = 0;
    defined_mat = true;
  }
  if(!defined_mat){
    G4cout << "Warning: Invalid DVCS ECal material: " << fDVCSECalMaterial 
	   << "; Use 'PbF2' or 'PbWO4' " << G4endl;
    return;
  }
  
  ////////////////////////////////////////////////////////                               
  // 36 columns, 31 rows 
  G4double mylarthickness = 0.0020*cm, airthickness = 0.0040*cm;
  G4double mylar_air_sum = mylarthickness + airthickness; 
  G4double dvcsblkpmtz = 0.20*cm;

  G4double calheight = Nrows*dvcsblkmodule_x+2*2.0*cm;
  G4double calwidth  = Ncols*dvcsblkmodule_y+2*2.0*cm;

  G4Box *dvcsblkecalbox = new G4Box("dvcsblkecalbox", calwidth/2.0, calheight/2.0, caldepth/2.0);
  G4LogicalVolume *dvcsblkecallog = new G4LogicalVolume(dvcsblkecalbox, GetMaterial("Air"), "dvcsblkecallog");
  G4ThreeVector dvcsblkecal_pos(EcalHoffset, EcalVoffset, fBBdist+caldepth/2.0);
  dvcsblkecal_pos.rotateY(fBBang);
  G4RotationMatrix* dvcsblkecal_rm = new G4RotationMatrix();
  dvcsblkecal_rm->rotateY(-fBBang);
  new G4PVPlacement( dvcsblkecal_rm, dvcsblkecal_pos, dvcsblkecallog, "dvcsblkecalphys", motherlog, false, 0 );
  
  // Calo module: 
  double DVCSblk_x = dvcsblkmodule_x - 2*mylar_air_sum;
  double DVCSblk_y = dvcsblkmodule_y - 2*mylar_air_sum;
  double DVCSblk_z = caldepth -2*2.0*cm;
  
  G4Box *dvcsblkmodbox = new G4Box("dvcsblkmodbox", dvcsblkmodule_x/2.0, dvcsblkmodule_y/2.0, caldepth/2.0);
  G4LogicalVolume *dvcsblkmodlog = new G4LogicalVolume(dvcsblkmodbox, GetMaterial("Special_Air"), "dvcsblkmodlog");

  G4Box *tempbox = new G4Box("tempbox", dvcsblkmodule_x/2.0, dvcsblkmodule_y/2.0, (caldepth-2*dvcsblkpmtz)/2.0);


  // calorimeter box Subtraction
    G4Box *dvcsblkmodbox_sub = new G4Box( "dvcsblkmodbox_sub", (dvcsblkmodule_x-2*mylarthickness)/2.0, (dvcsblkmodule_y-2*mylarthickness)/2.0, (caldepth-2*dvcsblkpmtz)/2.0 );

  G4SubtractionSolid *dvcsblkmylarwrap = new G4SubtractionSolid( "dvcsblkmylarwrap", tempbox, dvcsblkmodbox_sub, 0, G4ThreeVector(0.0, 0.0, mylarthickness) );
  G4LogicalVolume *dvcsblkmylarwraplog = new G4LogicalVolume( dvcsblkmylarwrap, GetMaterial("Mylar"), "dvcsblkmylarwraplog" ); 
  
 // new G4LogicalSkinSurface( "DVCSBLK Mylar Skin", dvcsblkmylarwraplog, GetOpticalSurface("Mirrsurf") );
  // Make Lead Glass 
  G4Box *DVCSblkbox = new G4Box( "DVCSblkbox", DVCSblk_x/2.0, DVCSblk_y/2.0, DVCSblk_z/2.0 );
  G4LogicalVolume *DVCSblklog = new G4LogicalVolume( DVCSblkbox, GetMaterial(fDVCSECalMaterial.data()), "DVCSblklog" );

  // Shower DVCSblk SD of type CAL
  G4SDManager *sdman = fDetCon->fSDman;

  G4String DVCSblkSDname = "DVCSblk";
  G4String DVCSblkcollname = "DVCSblkHitsCollection";
  G4SBSCalSD *DVCSblkSD = NULL;

  if( !((G4SBSCalSD*) sdman->FindSensitiveDetector(DVCSblkSDname)) ) {
    G4cout << "Adding DVCSblk Sensitive Detector to SDman..." << G4endl;
    DVCSblkSD = new G4SBSCalSD( DVCSblkSDname, DVCSblkcollname );
    
    sdman->AddNewDetector( DVCSblkSD );
    (fDetCon->SDlist).insert( DVCSblkSDname );
    fDetCon->SDtype[DVCSblkSDname] = G4SBS::kCAL;
    (DVCSblkSD->detmap).depth = 1;

    G4double threshold_default = 0.0*MeV; //1% of 1 GeV
    G4double timewindow_default = 100.0*ns; //We could use 10 ns here if we wanted, but also have to consider pulse shape. 

    fDetCon->SetTimeWindowAndThreshold( DVCSblkSDname, threshold_default, timewindow_default );
  }
  DVCSblklog->SetSensitiveDetector( DVCSblkSD ); 

  fDetCon->InsertSDboundaryVolume( dvcsblkecallog->GetName(), DVCSblkSDname );
//////////////////

  if( (fDetCon->StepLimiterList).find( DVCSblkSDname ) != (fDetCon->StepLimiterList).end() ){
    DVCSblklog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }
  // Make PMT/Window
  double pmtsize = 2.0*cm;
  G4Box *dvcsblkpmt = new G4Box( "dvcsblkpmt", pmtsize/2.0, pmtsize/2.0, dvcsblkpmtz/2.0 );
  G4LogicalVolume *dvcsblkpmtwindowlog = new G4LogicalVolume( dvcsblkpmt, GetMaterial("QuartzWindow_ECal"), "dvcsblkpmtwindowlog" );
  G4LogicalVolume *dvcsblkpmtcathodecallog = new G4LogicalVolume( dvcsblkpmt, GetMaterial("Photocathode_material_ecal"), "dvcsblkpmtcathodecallog" );

  // Shower PMT SD of type ECAL
  G4String DVCSblkecalSDname = "DVCSblkEcal";
  G4String DVCSblkecalcollname = "DVCSblkEcalHitsCollection";
  G4SBSECalSD *DVCSblkecalSD = NULL;

  if( !((G4SBSECalSD*) sdman->FindSensitiveDetector(DVCSblkecalSDname)) ) {
    G4cout << "Adding DVCSblkEcal Sensitive Detector to SDman..." << G4endl;
    DVCSblkecalSD = new G4SBSECalSD( DVCSblkecalSDname, DVCSblkecalcollname );
    sdman->AddNewDetector( DVCSblkecalSD );
    (fDetCon->SDlist).insert(DVCSblkecalSDname);
    fDetCon->SDtype[DVCSblkecalSDname] = G4SBS::kECAL;
    (DVCSblkecalSD->detmap).depth = 1;
  }
  dvcsblkpmtcathodecallog->SetSensitiveDetector( DVCSblkecalSD );

  fDetCon->InsertSDboundaryVolume( dvcsblkecallog->GetName(), DVCSblkecalSDname );

  // Put everything in a calo Module
  int mod_copy_number = 0;

  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-dvcsblkpmtz)/2.0-2.0*cm), dvcsblkpmtcathodecallog,"bbcathodephys", dvcsblkmodlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-3*dvcsblkpmtz)/2.0-2.0*cm), dvcsblkpmtwindowlog, "bbwindowphys", dvcsblkmodlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, (caldepth-4*dvcsblkpmtz-DVCSblk_z)/2.0-2.0*cm), DVCSblklog, "DVCSblkphys", dvcsblkmodlog, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, -dvcsblkpmtz-2.0*cm), dvcsblkmylarwraplog, "dvcsblkmylarphys", dvcsblkmodlog, false, 0 );

  for( int l=0; l<Ncols; l++ ) {
    for( int j=0; j<Nrows; j++ ) {

      (DVCSblkSD->detmap).Col[mod_copy_number] = l;
      (DVCSblkSD->detmap).Row[mod_copy_number] = j;
      (DVCSblkecalSD->detmap).Col[mod_copy_number] = l;
      (DVCSblkecalSD->detmap).Row[mod_copy_number] = j;
      double xtemp = (calwidth - dvcsblkmodule_x)/2.0 - 2.0*cm - l*dvcsblkmodule_x;
      double ytemp = (calheight - dvcsblkmodule_y)/2.0 - 2.0*cm - j*dvcsblkmodule_y;

      new G4PVPlacement(0, G4ThreeVector(xtemp,ytemp,0.0), dvcsblkmodlog, "calphys", dvcsblkecallog, false, mod_copy_number);
      
      (DVCSblkSD->detmap).LocalCoord[mod_copy_number] = G4ThreeVector( xtemp,ytemp,(caldepth-dvcsblkpmtz)/2.0  );
      (DVCSblkecalSD->detmap).LocalCoord[mod_copy_number] = G4ThreeVector( xtemp, ytemp, (caldepth-4*dvcsblkpmtz-DVCSblk_z)/2.0 );

      mod_copy_number++;
    }
  }

  // Visualization attributes
  G4VisAttributes *DVCSblkecalbox_visatt = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7) );
  DVCSblkecalbox_visatt->SetForceWireframe(true);
  dvcsblkecallog->SetVisAttributes( DVCSblkecalbox_visatt );
    
  //G4VisAttributes *mydvcsblkmodbox_visatt = new G4VisAttributes(G4Colour(1.0, 0.0, 0.0) );
  dvcsblkmodlog->SetVisAttributes( DVCSblkecalbox_visatt );// G4VisAttributes::Invisible );//
  
  dvcsblkmylarwraplog->SetVisAttributes( G4VisAttributes::Invisible );
  
  //TF1
  G4VisAttributes *DVCSblk_visatt = new G4VisAttributes(G4Colour( 1.0, 1.0, 0.0 ) );
  DVCSblklog->SetVisAttributes( DVCSblk_visatt);

  //PMTcathode
  G4VisAttributes *PMT_visatt = new G4VisAttributes(G4Colour( 0.0, 0.0, 1.0 ));
  dvcsblkpmtcathodecallog->SetVisAttributes( PMT_visatt);
}

void G4SBSEArmBuilder::MakeC16( G4LogicalVolume *motherlog ){
  printf("C16 at %f deg\n", fBBang/deg);

  G4SDManager *sdman = fDetCon->fSDman;

  // Module specs
  G4double width_42 = 4.2*cm;
  G4double depth_leadglass = 34.3*cm;
  G4int nrows = 4;
  G4int ncols = 4;

  // Wave-guide specs
  G4double depth_WG  = 15.0*cm;
  G4double radius_WG = 1.25*cm;
  
  // PMT specs
  G4double depth_ecal_pmt = 0.3*cm;
  G4double radius_ecal_pmt = 1.25*cm;

  // We should be using aluminized foil instead of mylar.
  G4double alum_thick = 0.001*2.54*cm;
  G4double air_thick = alum_thick;

  // Make a C16 Mother Volume to house all the modules+WG+PMTs:
  G4double C16_width  = nrows*width_42;
  G4double C16_depth  = depth_leadglass + depth_WG + depth_ecal_pmt;

  // Geometric information & C16 placing
  G4double r_C16 = fBBdist + C16_depth/2.0;      // Distance from target to front of C16
  G4double angle_C16 = fBBang;                    // C16 BB angle
  G4ThreeVector R_C16( r_C16*sin( angle_C16 ), 0.0, r_C16*cos( angle_C16 ) );

  G4RotationMatrix *bbrm_C16 = new G4RotationMatrix;
  bbrm_C16->rotateY( -angle_C16 );

  G4Box *C16_Box = new G4Box( "C16_Box", C16_width/2.0, C16_width/2.0, C16_depth/2.0 );
  G4LogicalVolume *C16_Log = new G4LogicalVolume( C16_Box, GetMaterial("Special_Air"), "C16_Log" );
  new G4PVPlacement( bbrm_C16, R_C16, C16_Log, "C16_Phys", motherlog, false, 0 );

  // Getting the PMTs set up - use same PMT as ECal, same sensitivity
  G4Tubs *ecal_PMT = new G4Tubs( "ecal_PMT", 0.0, radius_ecal_pmt, depth_ecal_pmt/2.0, 0.0, twopi );
  G4LogicalVolume *ecal_PMT_log = new G4LogicalVolume( ecal_PMT, GetMaterial("Photocathode_material_ecal"), "ecal_PMT_log" );

  G4String C16SDname = "Earm/C16";
  G4String collname = "C16HitsCollection";
  G4SBSECalSD *C16SD = NULL;

  if( !((G4SBSECalSD*) sdman->FindSensitiveDetector( C16SDname )) ){
    G4cout << "Adding C16 PMT sensitive detector to SDman..." << G4endl;
    C16SD = new G4SBSECalSD( C16SDname, collname );
    sdman->AddNewDetector( C16SD );
    (fDetCon->SDlist).insert( C16SDname );
    fDetCon->SDtype[C16SDname] = G4SBS::kECAL;
    (C16SD->detmap).depth = 0;
  }

  ecal_PMT_log->SetSensitiveDetector( C16SD );

  fDetCon->InsertSDboundaryVolume( C16_Log->GetName(), C16SDname );
  
  // Getting the Wave Guides set up, *****need to UPDATE material properties*****
  G4Tubs *C16_WG = new G4Tubs( "C16_WG", 0.0, radius_WG, depth_WG/2.0, 0.0, twopi );
  G4LogicalVolume *C16_WG_Log = new G4LogicalVolume( C16_WG, GetMaterial("Pyrex_Glass"), "C16_WG_Log" );
  //new G4LogicalSkinSurface( "C16_WG_Skin", C16_WG_Log, GetOpticalSurface("osWLSToAir") );

  // Make a trapzoidal aluminum piece that starts at the end of the TF1 and ends
  // half way down the WG
  G4double  dx1 = width_42 / 2.0;
  G4double  dx2 = (radius_WG + 2*alum_thick + 2*air_thick);
  G4double  dy1 = width_42 / 2.0;
  G4double  dy2 = (radius_WG + 2*alum_thick + 2*air_thick);
  G4double  dz  = depth_WG / 4.0;
  G4Trd *Al_endpiece = new G4Trd( "Al_endpiece", dx1, dx2, dy1, dy2, dz );

  // Make a subtraction in order to get Trapezoidal Al foil
  dx1 = (width_42 - 2*alum_thick - 2*air_thick) / 2.0;
  dx2 = radius_WG;
  dy1 = (width_42 - 2*alum_thick - 2*air_thick) / 2.0;
  dy2 = radius_WG;
  dz = (depth_WG + 0.1*cm) / 4.0;
  G4Trd *Al_wrap_endpiece_sub = new G4Trd( "Al_wrap_endpiece_sub", dx1, dx2, dy1, dy2, dz );
  G4SubtractionSolid *Al_wrap_endpiece = new G4SubtractionSolid( "Al_wrap_endpiece", Al_endpiece, Al_wrap_endpiece_sub, 0, G4ThreeVector(0.0,0.0,0.0) );
  G4LogicalVolume *Al_wrap_endpiece_log = new G4LogicalVolume( Al_wrap_endpiece, GetMaterial("Aluminum"), "Al_wrap_endpiece_log" );
  new G4LogicalSkinSurface( "Al_endpiece_skin", Al_wrap_endpiece_log, GetOpticalSurface("Mirrsurf") );

  // There is a 10cm Foam Glass wrap around everything except the PMT area
  G4double foam_thick = 10.0*cm;
  G4double foam_XY = C16_width + 2*foam_thick;
  G4double foam_Z  = C16_depth + foam_thick - depth_ecal_pmt;
  G4Box *foam_box = new G4Box( "foam_box", foam_XY/2.0, foam_XY/2.0, foam_Z/2.0 );

  // Cut out the inside, but keep 10cm on one longitudinal face
  G4Box *foam_sub = new G4Box( "foam_sub", C16_width/2.0, C16_width/2.0, foam_Z/2.0 );
  G4SubtractionSolid *foam_wrap = new G4SubtractionSolid( "foam_wrap", foam_box, foam_sub, 0, G4ThreeVector(0.0,0.0,foam_thick) );
  G4LogicalVolume *foam_wrap_log = new G4LogicalVolume( foam_wrap, GetMaterial("SiO2_C16"), "foam_wrap_log" );
  G4ThreeVector R_TF1_Foam( (depth_ecal_pmt/2.0+foam_thick/2.0)*sin( angle_C16 ), 0.0, (depth_ecal_pmt/2.0+foam_thick/2.0)*cos( angle_C16 ) );
  G4ThreeVector R_Foam = R_C16 - R_TF1_Foam;
  new G4PVPlacement(bbrm_C16, R_Foam, foam_wrap_log, "Foam_phys", motherlog, false, 0 ,true);

  // There is a small Aluminium plate 5" upstream of the lead-glass face
  G4double Al_depth = 0.25*2.54*cm;
  G4double Al_width = 6.0*2.54*cm;
  G4double Distfrom_TF1 = 5.0*2.54*cm + C16_depth/2.0;
  G4ThreeVector R_TF1_Al( Distfrom_TF1*sin( angle_C16 ), 0.0, Distfrom_TF1*cos( angle_C16 ) );
  G4ThreeVector R_Al = R_C16 - R_TF1_Al;

  G4Box *Al_Plate_Box = new G4Box( "Al_Plate_Box", Al_width/2.0, Al_width/2.0, Al_depth/2.0 );
  G4LogicalVolume *Al_Plate_Log = new G4LogicalVolume( Al_Plate_Box, GetMaterial("Aluminum"), "Al_Plate_Log" );
  new G4PVPlacement( bbrm_C16, R_Al, Al_Plate_Log, "Al_Plate", motherlog, false, 0 );

  // VISUALS
  G4VisAttributes *TF1visatt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0 ) );
  G4VisAttributes *Alvisatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
  G4VisAttributes *ECALpmtvisatt = new G4VisAttributes( G4Colour( 0, 0, 1 ) );
  G4VisAttributes *C16WG_visatt = new G4VisAttributes( G4Colour(0.54, 0.53, 0.79) );
  G4VisAttributes *Foam_visatt = new G4VisAttributes( G4Colour( 0.0, 0.6, 0.6) );
  // C16 Mother:
  C16_Log->SetVisAttributes( G4VisAttributes::Invisible );
  // Al Plate & Foil:
  Alvisatt->SetForceWireframe(true);
  Al_Plate_Log->SetVisAttributes( Alvisatt );
  Al_wrap_endpiece_log->SetVisAttributes( Alvisatt );
  // PMTs:
  ecal_PMT_log->SetVisAttributes( ECALpmtvisatt );
  // WaveGuides:
  C16_WG_Log->SetVisAttributes( C16WG_visatt );
  // Foam Wrap
  Foam_visatt->SetForceWireframe(true);
  foam_wrap_log->SetVisAttributes( Foam_visatt );

  /////////////////////////////////////////////////////////////////////////////////
  //   There are two options to build the TF1 ( /g4sbs/segmentC16 int )          //
  //     /g4sbs/segmentC16 N segments the TF1 into N "equal" sections, used       //
  //     to analyze dose rate as a function of longitudinal dimension of module  //
  //                                                                             //
  //     /g4sbs/segmentC16 0 builds the normal ECal modules                      //
  /////////////////////////////////////////////////////////////////////////////////

  if( fDetCon->GetC16Segmentation() <= 0 ) { //Default case: room temperature, no radiation damage:
    // Make a C16 Module  
    G4Box *Module_42 = new G4Box( "Module_42", width_42/2.0, width_42/2.0, depth_leadglass/2.0 );
    G4LogicalVolume *Module_42_log = new G4LogicalVolume( Module_42, GetMaterial("Special_Air"), "Module_42_log" );
    // Next, we want to make a subtraction solid for the Al:
    G4Box *Al_42 = new G4Box( "Al_42", (width_42 - alum_thick)/2.0, (width_42 - alum_thick)/2.0, depth_leadglass/2.0 + 1.0*cm );
    //
    G4SubtractionSolid *Al_wrap_42 = new G4SubtractionSolid( "Al_wrap_42", Module_42, Al_42, 0, G4ThreeVector( 0, 0, alum_thick + 1.0*cm ) );
    G4LogicalVolume *Al_wrap_42_log = new G4LogicalVolume( Al_wrap_42, GetMaterial("Aluminum"), "Al_wrap_42_log" );
    // Make lead-glass
    G4Box *LeadGlass_42 = new G4Box("LeadGlass_42", (width_42 - alum_thick - air_thick)/2.0, (width_42 - alum_thick - air_thick)/2.0, (depth_leadglass - alum_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_42_log = new G4LogicalVolume( LeadGlass_42, GetMaterial("TF1"), "LeadGlass_42_log" );

    // Define Sensitive Detector for lead-glass of type kCAL
    G4String C16TF1SDname = "Earm/C16TF1";
    G4String C16TF1collname = "C16TF1HitsCollection";
    G4SBSCalSD *C16TF1SD = NULL;
    
    if( !( C16TF1SD = (G4SBSCalSD*) sdman->FindSensitiveDetector(C16TF1SDname) ) ){
      G4cout << "Adding C16 TF1 Sensitive Detector to SDman..." << G4endl;
      C16TF1SD = new G4SBSCalSD( C16TF1SDname, C16TF1collname );

      fDetCon->fSDman->AddNewDetector( C16TF1SD );
      (fDetCon->SDlist).insert( C16TF1SDname );
      fDetCon->SDtype[C16TF1SDname] = G4SBS::kCAL;
      (C16TF1SD->detmap).depth = 1;

      G4double default_threshold = 0.0*MeV;
      G4double default_timewindow = 100.0*ns;

      fDetCon->SetTimeWindowAndThreshold( C16TF1SDname, default_threshold, default_timewindow );
    }
    // Assign "kCAL" sensitivity to the lead-glass:
    LeadGlass_42_log->SetSensitiveDetector( C16TF1SD );

    fDetCon->InsertSDboundaryVolume( C16_Log->GetName(), C16TF1SDname );
    
    // Place lead-glass and Al wrap inside module:
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, (alum_thick + air_thick)/2.0 ), LeadGlass_42_log, "LeadGlass_42_phys", Module_42_log, false, 0 );
    // Al:
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), Al_wrap_42_log, "Al_wrap_42_phys", Module_42_log, false, 0 );
    new G4LogicalSkinSurface( "Al_skin_42", Al_wrap_42_log, GetOpticalSurface("Mirrsurf") );

    // Construct C16
    G4int copyID = 0;
    for( G4int i=0; i<nrows; i++ ) {
      for( G4int j=0; j<ncols; j++ ) {
	G4double tempx = C16_width/2.0 - width_42/2.0 - j*width_42;
	G4double tempy = C16_width/2.0 - width_42/2.0 - i*width_42;
	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy, -depth_WG/2.0 - depth_ecal_pmt/2.0), Module_42_log, "C16_Module", C16_Log, false, copyID );
	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy, -C16_depth/2.0 + depth_leadglass + depth_WG/4.0), Al_wrap_endpiece_log, "Aluminum_wrap_endpiece_phys", C16_Log, false, copyID );
	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy, -depth_ecal_pmt/2.0 + depth_leadglass/2.0), C16_WG_Log, "C16_WG", C16_Log, false, copyID );
	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy,  depth_leadglass/2.0 + depth_WG/2.0), ecal_PMT_log, "C16_PMT", C16_Log, false, copyID );

	(C16TF1SD->detmap).Row[copyID] = i;
	(C16TF1SD->detmap).Col[copyID] = j;
	(C16TF1SD->detmap).LocalCoord[copyID] = G4ThreeVector( tempx, tempy, 0.0 );

	(C16SD->detmap).Row[copyID] = i;
	(C16SD->detmap).Col[copyID] = j;
	(C16SD->detmap).LocalCoord[copyID] = G4ThreeVector( tempx, tempy, 0.0 );
	
	copyID++;
      }
    }
    // Set Visuals
    Module_42_log->SetVisAttributes( G4VisAttributes::Invisible );
    LeadGlass_42_log->SetVisAttributes( TF1visatt );
    Al_wrap_42_log->SetVisAttributes( Alvisatt );  
  } else {
    // The strategy is to place the aluminum wrap within C16_Log (the mother volume), then iteratively
    // place TF1 modules in the longitudinal direction. Therefore, we can define different material properties
    // to each TF1 "segment" within a C16 Module.

    // Make a C16 Module  
    G4Box *Module_42 = new G4Box( "Module_42", width_42/2.0, width_42/2.0, depth_leadglass/2.0 );
    G4LogicalVolume *Module_42_log = new G4LogicalVolume( Module_42, GetMaterial("Special_Air"), "Module_42_log" );
    
    // Next, we want to make a subtraction solid for the Al:
    G4Box *Al_42 = new G4Box( "Al_42", (width_42 - alum_thick)/2.0, (width_42 - alum_thick)/2.0, depth_leadglass/2.0 + 1.0*cm );
    G4SubtractionSolid *Al_wrap_42 = new G4SubtractionSolid( "Al_wrap_42", Module_42, Al_42, 0, G4ThreeVector( 0, 0, alum_thick + 1.0*cm ) );
    G4LogicalVolume *Al_wrap_42_log = new G4LogicalVolume( Al_wrap_42, GetMaterial("Aluminum"), "Al_wrap_42_log" );
    new G4LogicalSkinSurface( "Al_skin_42", Al_wrap_42_log, GetOpticalSurface("Mirrsurf") );

    // Make TF1 - Logical Volume will be defined iteratively below in order to change material 
    // properties based on segmentation 
    //G4double Nsegments = 10.0;
    //    G4double Nsegments = G4double(fDetCon->GetC16Segmentation());
    //G4double segment_depth = (depth_leadglass - alum_thick - air_thick) / Nsegments;
    //Note: Because segment thickness is always the same for the detector and the material definition,
    // as long as nsegments * segmentthick > depth of lead glass, the material definition will always exist!
    G4double segment_depth = fDetCon->GetSegmentThickC16();

    G4int Nsegments = G4int( (depth_leadglass - alum_thick - air_thick ) / segment_depth ); //truncate the remainder.
    
    G4double remainder = (depth_leadglass - alum_thick - air_thick ) - Nsegments * segment_depth; //remainder is always non-negative!

    if( Nsegments > fDetCon->GetC16Segmentation() ){
      Nsegments = fDetCon->GetC16Segmentation();
      remainder = (depth_leadglass - alum_thick - air_thick ) - Nsegments * segment_depth;
    }
    
    if( Nsegments == 0 ){
      Nsegments = 1;
      segment_depth = (depth_leadglass - alum_thick - air_thick );
      remainder = 0.0;
    }
    
    G4Box *Segments_TF1 = new G4Box( "Segments_TF1", (width_42 - 2.0*alum_thick - 2.0*air_thick)/2.0, 
					   (width_42 - 2.0*alum_thick - 2.0*air_thick)/2.0, segment_depth/2.0 );

    //Always add the remainder to the last segment, since optical properties will be varying less with z deeper in the glass anyway!
    G4Box *LastSegment_TF1 = new G4Box( "LastSegment_TF1",
					(width_42 - 2.0*alum_thick - 2.0*air_thick)/2.0, 
					(width_42 - 2.0*alum_thick - 2.0*air_thick)/2.0,
					(segment_depth + remainder)/2.0 );  
    // TF1 is a Sensitive Detector of type kCAL. Sensitivity will be assigned to a LV
    // within the loop below:
    G4String C16TF1SDname = "Earm/C16TF1";
    G4String C16TF1collname = "C16TF1HitsCollection";
    G4SBSCalSD *C16TF1SD = NULL;
    
    if( !( C16TF1SD = (G4SBSCalSD*) sdman->FindSensitiveDetector(C16TF1SDname) ) ){
      G4cout << "Adding C16 TF1 Segmented Sensitive Detector to SDman..." << G4endl;
      C16TF1SD = new G4SBSCalSD( C16TF1SDname, C16TF1collname );
      
      fDetCon->fSDman->AddNewDetector( C16TF1SD );
      (fDetCon->SDlist).insert( C16TF1SDname );
      fDetCon->SDtype[C16TF1SDname] = G4SBS::kCAL;
      (C16TF1SD->detmap).depth = 0;

      G4double default_threshold = 0.0*MeV;
      G4double default_timewindow = 100.0*ns;

      fDetCon->SetTimeWindowAndThreshold( C16TF1SDname, default_threshold, default_timewindow );
    }

    G4int cell_number = 0 ;    // cell #
    G4int TF1_number = 0 ;     // TF1 identifier
    for( G4int i = 0; i < nrows; i++ ) {
      for( G4int j = 0; j < ncols; j++ ) {
	G4double tempx = C16_width/2.0 - width_42/2.0 - j*width_42;
	G4double tempy = C16_width/2.0 - width_42/2.0 - i*width_42;
	
	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy, -C16_depth/2.0 + depth_leadglass/2.0), 
			   Al_wrap_42_log, "Aluminum_wrap_phys", C16_Log, false, cell_number );

	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy, -C16_depth/2.0 + depth_leadglass + depth_WG/4.0), 
			   Al_wrap_endpiece_log, "Aluminum_wrap_endpiece_phys", C16_Log, false, cell_number );

	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy, C16_depth/2.0 - depth_ecal_pmt - depth_WG/2.0), 
			   C16_WG_Log, "C16_WG_phys", C16_Log, false, cell_number );

	new G4PVPlacement( 0, G4ThreeVector(tempx, tempy, C16_depth/2.0 - depth_ecal_pmt/2.0), 
			   ecal_PMT_log, "C16_PMT_phys", C16_Log, false, cell_number );

	(C16SD->detmap).Row[cell_number] = i;
	(C16SD->detmap).Col[cell_number] = j;
	(C16SD->detmap).LocalCoord[cell_number] = G4ThreeVector( tempx, tempy, 0.0 );
	
	cell_number++;

	for( G4int planeN = 0; planeN < Nsegments; planeN++ ) {
	  // ostringstream temp, temp1, temp2;
	  // temp << planeN;// segment #
	  // temp1 << i;    // row #
	  // temp2 << j;    // col #
	  // G4String tempstring  = temp.str();
	  // G4String tempstring1 = temp1.str();
	  // G4String tempstring2 = temp2.str();
	  // G4String seg_TF1_name_log  = "TF1_log_seg_"  + tempstring + "_row_" + tempstring1 + "_col_" + tempstring2;
	  // G4String seg_TF1_name_phys = "TF1_phys_seg_" + tempstring + "_row_" + tempstring1 + "_col_" + tempstring2;
	  // G4String seg_TF1_material = "TF1_anneal_" + tempstring;

	  TString material_name, lv_name, pv_name;
	  material_name.Form( "TF1_anneal_C16_row%d_col%d_z%d", i+1, j+1, planeN );

	  lv_name.Form( "TF1_log_z%d", planeN );
	  pv_name.Form( "TF1_pv_z%d", planeN );
	  
	  // Assign each TF1 segment a kCAL Sensitivity
	  G4LogicalVolume *Segments_TF1_log;
	  if( planeN + 1 < Nsegments || Nsegments == 1 ){
	    Segments_TF1_log = new G4LogicalVolume( Segments_TF1, GetMaterial(material_name.Data()), lv_name.Data() );
	  } else {
	    Segments_TF1_log = new G4LogicalVolume( LastSegment_TF1, GetMaterial(material_name.Data()), lv_name.Data() );
	  }
	    
	  Segments_TF1_log->SetSensitiveDetector( C16TF1SD );

	  G4VisAttributes *Segment_VisAtt = new G4VisAttributes( G4Colour( 0.8*(planeN/15.0)+0.20, 0.8*(planeN/15.0)+0.20, 0.0 ) );
	  Segments_TF1_log->SetVisAttributes( Segment_VisAtt );

	  // Place the TF1 segments longitudinally down the module
	  // Therefore, seg_0 corresponds to the face of C16
	  G4double tempz;
	  if( planeN + 1 < Nsegments ){
	    tempz = -C16_depth/2.0 + air_thick + alum_thick + (planeN + 0.5)*segment_depth;
	  } else {
	    tempz = -C16_depth/2.0 + air_thick + alum_thick + (Nsegments-1)*segment_depth +
	      0.5* ( segment_depth + remainder );
	  }
	  //G4double tempz = -C16_depth/2.0 + segment_depth/2.0 + air_thick + alum_thick + planeN*segment_depth;
	    
	  new G4PVPlacement(0, G4ThreeVector(tempx, tempy, tempz), Segments_TF1_log, pv_name.Data(), C16_Log, false, TF1_number );
	  
	  // Record useful information using Detmap:
	  (C16TF1SD->detmap).Row[TF1_number] = i;
	  (C16TF1SD->detmap).Col[TF1_number] = j;
	  (C16TF1SD->detmap).LocalCoord[TF1_number] = G4ThreeVector(tempx,tempy,tempz);
	  (C16TF1SD->detmap).Plane[TF1_number] = planeN;
	  TF1_number++;
	}
      }
    }
    // Set Visuals
    Module_42_log->SetVisAttributes( G4VisAttributes::Invisible );
    Al_wrap_42_log->SetVisAttributes( Alvisatt );
  }
}

void G4SBSEArmBuilder::MakeBigCal(G4LogicalVolume *motherlog){
  printf("BigCal at %f deg\n", fBBang/deg);

  // Pointer to SDmanager, used frequently in this routine
  G4SDManager *sdman = fDetCon->fSDman;

  //Parameters of ECAL for GEP:
  G4double xfpmin_cal = -160.0*cm;
  G4double xfpmax_cal = 160.0*cm;

  //number of available blocks:
  G4int nblocks_42 = 1000;
  G4int nblocks_40 = 700;
  G4int nblocks_38 = 300;

  //assume for now that all lead-glass blocks are 40 cm deep:
  G4double width_42 = 4.2*cm;
  G4double width_40 = 4.0*cm;
  G4double width_38 = 3.8*cm;

  G4double depth_42 = 34.3*cm;
  G4double depth_40 = 40.0*cm;
  G4double depth_38 = 45.0*cm;
  
  G4double depth_leadglass = 45.0*cm;
  G4double depth_ecal_pmt = 0.3*cm;
  G4double depth_lightguide_short = 15.0*cm;
  G4double radius_ecal_pmt = 1.25*cm;
  G4double depth_ecal_frontplate = 2.54*cm;
  G4double depth_CH2 = 20.0*cm; //This goes directly in front of CDET:
  G4double depth_CDET = 45.0*cm; // CDET starts 45 cm in front of ECAL:
  
  //Define "Earm" box a bit wider than total ECAL area:
  G4double width_earm = 150.0*cm;
  G4double height_earm = 340.0*cm;
  G4double depth_earm = depth_CH2 + depth_CDET + depth_ecal_pmt + depth_leadglass + depth_lightguide_short; // 
  
  G4int nSuperRows = 20;
  //These are the "yfp" coordinates of the box lower and upper edges. +yfp
  //Note: 0-9 and 10-19 are mirror images of each other!
  G4double yfpmin_super_rows[20] = { -56.0*cm,
				   -64.0*cm, -64.0*cm, -64.0*cm, -64.0*cm, -64.0*cm, -64.0*cm,
				   -56.0*cm,
				   -48.0*cm, -48.0*cm, -48.0*cm, -48.0*cm,
				   -56.0*cm,
				   -64.0*cm, -64.0*cm, -64.0*cm, -64.0*cm, -64.0*cm, -64.0*cm,
				   -56.0*cm };

  G4double yfpmax_super_rows[20] = { -24.0*cm, -8.0*cm, 8.0*cm, 24.0*cm, 32.0*cm,
				     48.0*cm, 48.0*cm, 56.0*cm, 56.0*cm, 56.0*cm,
				     56.0*cm, 56.0*cm, 56.0*cm, 48.0*cm, 48.0*cm,
				     32.0*cm, 24.0*cm, 8.0*cm, -8.0*cm, -24.0*cm };

  //G4double ycalo_min = -160.0*cm; //This can be adjusted if necessary!
  G4double ycalo_min = -164.0*cm;
  
  G4double xcalomin_super_rows[20], xcalomax_super_rows[20], width_super_rows[20];
  G4int ncol42_row[20], ncol40_row[20], ncol38_row[20];
  G4int nblock42_super_row[20], nblock40_super_row[20], nblock38_super_row[20];
  for(G4int row=0; row<20; row++){//
    xcalomin_super_rows[row] = yfpmin_super_rows[row];
    xcalomax_super_rows[row] = yfpmax_super_rows[row];
    width_super_rows[row] = xcalomax_super_rows[row] - xcalomin_super_rows[row];
    ncol42_row[row] = G4int( width_super_rows[row]/width_42 ) + 1; //truncate and add 1:
    ncol40_row[row] = G4int( width_super_rows[row]/width_40 ) + 1;
    ncol38_row[row] = G4int( width_super_rows[row]/width_38 ) + 1;

    //compute number of blocks required to fill a "super row":
    //Note that this is always a multiple of 4
    nblock42_super_row[row] = 4*ncol42_row[row];
    nblock40_super_row[row] = 4*ncol40_row[row];
    nblock38_super_row[row] = 4*ncol38_row[row];
  }
				     
  // Rotations, offsets, and distances from Target:
  G4RotationMatrix *bbrm = new G4RotationMatrix;
  bbrm->rotateY(-fBBang);

  G4Box *earm_mother_box = new G4Box( "earm_mother_box", width_earm/2.0, height_earm/2.0, (depth_earm+1.0*mm)/2.0 );
  G4LogicalVolume *earm_mother_log = new G4LogicalVolume( earm_mother_box, GetMaterial("Air"), "earm_mother_log");

  //If "earm_mother_log" is in the step limiter list, make ECAL a total-absorber with "kCAL" sensitivity:
  if( (fDetCon->StepLimiterList).find( "earm_mother_log" ) != (fDetCon->StepLimiterList).end() ){
    earm_mother_log->SetUserLimits( new G4UserLimits( 0.0, 0.0, 0.0, DBL_MAX, DBL_MAX ) );

    G4String sdname = "Earm/ECAL_box";
    G4String collname = "ECAL_boxHitsCollection";
    G4SBSCalSD *earm_mother_SD = NULL;
    if( !( (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(sdname) ) ){
      G4cout << "Adding ECAL_box sensitive detector to SDman..." << G4endl;
      earm_mother_SD = new G4SBSCalSD( sdname, collname );
      fDetCon->fSDman->AddNewDetector( earm_mother_SD );
      (fDetCon->SDlist).insert( sdname );
      fDetCon->SDtype[sdname] = G4SBS::kCAL;
      (earm_mother_SD->detmap).depth = 0;
   
      earm_mother_log->SetSensitiveDetector( earm_mother_SD );
    }
  }

  //fBBdist should now be interpreted to mean the distance from the origin to the ****FRONT**** of the ECAL lead-glass:
  //G4double zback_ECAL = depth_earm/2.0 - depth_ecal_pmt;
  G4double zfront_ECAL = depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass;
  G4double R_Earm = fBBdist - zfront_ECAL;

  G4ThreeVector pos_ECAL( R_Earm*sin(fBBang), 0.0, R_Earm*cos(fBBang) );
  
  new G4PVPlacement( bbrm, pos_ECAL, earm_mother_log, "earm_mother_phys", motherlog, false, 0 );

  //assume lead-glass is surrounded by 1-mil (.001") thickness of mylar:
  //assume lead-glass is also surrounded by another mil of air:
  G4double mylar_thick = 0.001*2.54*cm;
  G4double air_thick = mylar_thick;
  
  G4double copper_thick = 0.25*mm;
  G4double al_thick = 0.50*mm;
  G4VisAttributes* hc_visAtt = new G4VisAttributes();

  // G4double hcf_thick = copper_thick;
  // const char* hcf_mat_name = "Copper";
  // hc_visAtt->SetColour(1.0, 0.5, 0.0);
  // G4double hcf_thick = al_thick;
  // const char* hcf_mat_name = "Aluminum";
  // hc_visAtt->SetColour(0.7, 0.7, 0.7);
  G4double hcf_thick = 0.0;
  const char* hcf_mat_name = "Special_Air";
  hc_visAtt->SetVisibility(0);
  
  //EPAF 2017-01-11: Declaring sensitive detector for light guide 
  // shall be temporary, and not end in the repo...
  // G4String ECalLGSDname = "Earm/ECalLG";
  // G4String ECalLGcollname = "ECalLGHitsCollection";
  // G4SBSCalSD *ECalLGSD = NULL;
  // if( !( ECalLGSD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(ECalLGSDname) ) ){
  //   G4cout << "Adding ECal light guide Sensitive Detector to SDman..." << G4endl;
  //   ECalLGSD = new G4SBSCalSD( ECalLGSDname, ECalLGcollname );
  //   fDetCon->fSDman->AddNewDetector( ECalLGSD );
  //   (fDetCon->SDlist).insert(ECalLGSDname);
  //   fDetCon->SDtype[ECalLGSDname] = kCAL;
  //   (ECalLGSD->detmap).depth = 1;//?????
  // }
  
  //Now place things in ECAL:
  //Start with the lead-glass modules, PMTs and light guides:
  G4Tubs *LightGuide_42 = new G4Tubs("LightGuide_42", 0.0, 2.5*cm/2.0, (depth_lightguide_short+depth_38-depth_42)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LightGuide_42_log = new G4LogicalVolume( LightGuide_42, GetMaterial("Pyrex_Glass"), "LightGuide_42_log" );
  
  G4Tubs *LightGuide_40 = new G4Tubs("LightGuide_40", 0.0, 2.5*cm/2.0, (depth_lightguide_short+depth_38-depth_40)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LightGuide_40_log = new G4LogicalVolume( LightGuide_40, GetMaterial("Pyrex_Glass"), "LightGuide_40_log" );

  G4Tubs *LightGuide_38 = new G4Tubs("LightGuide_38", 0.0, 2.5*cm/2.0, (depth_lightguide_short+depth_38-depth_38)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LightGuide_38_log = new G4LogicalVolume( LightGuide_38, GetMaterial("Pyrex_Glass"), "LightGuide_38_log" );

  //EPAF 2017-01-11: Need to make sensitive the three volumes above, to measure their dose.
  // shall be temporary, and not end in the repo...
  // LightGuide_42_log->SetSensitiveDetector( ECalLGSD );
  // LightGuide_40_log->SetSensitiveDetector( ECalLGSD );
  // LightGuide_38_log->SetSensitiveDetector( ECalLGSD );
  
  G4Tubs *LGWrap_42 = new G4Tubs( "LGWrap_42", 2.5*cm/2.0+air_thick, 2.5*cm/2.0 + air_thick + mylar_thick, (depth_lightguide_short+depth_38-depth_42)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LGWrap_42_log = new G4LogicalVolume( LGWrap_42, GetMaterial("Mylar"), "LGWrap_42_log" );

  G4Tubs *LGWrap_40 = new G4Tubs( "LGWrap_40", 2.5*cm/2.0+air_thick, 2.5*cm/2.0 + air_thick + mylar_thick, (depth_lightguide_short+depth_38-depth_40)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LGWrap_40_log = new G4LogicalVolume( LGWrap_40, GetMaterial("Mylar"), "LGWrap_40_log" );

  G4Tubs *LGWrap_38 = new G4Tubs( "LGWrap_38", 2.5*cm/2.0+air_thick, 2.5*cm/2.0 + air_thick + mylar_thick, (depth_lightguide_short+depth_38-depth_38)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LGWrap_38_log = new G4LogicalVolume( LGWrap_38, GetMaterial("Mylar"), "LGWrap_38_log" );
  
  G4Tubs *LG42 = new G4Tubs("LG42", 0.0, 2.5*cm/2.0 + air_thick + mylar_thick+0.1*mm, (depth_lightguide_short+depth_38-depth_42)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LG42_log = new G4LogicalVolume( LG42, GetMaterial("Special_Air"), "LG42_log" );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LightGuide_42_log, "LightGuide_42_phys", LG42_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LGWrap_42_log, "LGWrap_42_phys", LG42_log, false, 0 );

  G4Tubs *LG40 = new G4Tubs("LG40", 0.0, 2.5*cm/2.0 + air_thick + mylar_thick+0.1*mm, (depth_lightguide_short+depth_38-depth_40)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LG40_log = new G4LogicalVolume( LG40, GetMaterial("Special_Air"), "LG40_log" );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LightGuide_40_log, "LightGuide_40_phys", LG40_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LGWrap_40_log, "LGWrap_40_phys", LG40_log, false, 0 );
  
  G4Tubs *LG38 = new G4Tubs("LG38", 0.0, 2.5*cm/2.0 + air_thick + mylar_thick+0.1*mm, (depth_lightguide_short+depth_38-depth_38)/2.0, 0.0*deg, 360.0*deg );
  G4LogicalVolume *LG38_log = new G4LogicalVolume( LG38, GetMaterial("Special_Air"), "LG38_log" );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LightGuide_38_log, "LightGuide_38_phys", LG38_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector(0,0,0), LGWrap_38_log, "LGWrap_38_phys", LG38_log, false, 0 );

  G4Box *Module_42 = new G4Box( "Module_42", width_42/2.0, width_42/2.0, depth_42/2.0 );
  G4LogicalVolume *Module_42_log = new G4LogicalVolume( Module_42, GetMaterial("Special_Air"), "Module_42_log" );

  G4Box *Module_40 = new G4Box( "Module_40", width_40/2.0, width_40/2.0, depth_40/2.0 );
  G4LogicalVolume *Module_40_log = new G4LogicalVolume( Module_40, GetMaterial("Special_Air"), "Module_40_log" );

  G4Box *Module_38 = new G4Box( "Module_38", width_38/2.0, width_38/2.0, depth_38/2.0 );
  G4LogicalVolume *Module_38_log = new G4LogicalVolume( Module_38, GetMaterial("Special_Air"), "Module_38_log" );
  
  /*
  //Next, we want to make a subtraction solid for the mylar:
  G4Box *Mylar_42 = new G4Box( "Mylar_42", width_42/2.0 - mylar_thick, width_42/2.0 - mylar_thick, depth_42/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_42 = new G4SubtractionSolid( "Mylar_wrap_42", Module_42, Mylar_42, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm ) );
  G4LogicalVolume *Mylar_wrap_42_log = new G4LogicalVolume( Mylar_wrap_42, GetMaterial("Mylar"), "Mylar_wrap_42_log" );
  
  G4Box *Mylar_40 = new G4Box( "Mylar_40", width_40/2.0 - mylar_thick, width_40/2.0 - mylar_thick, depth_40/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_40 = new G4SubtractionSolid( "Mylar_wrap_40", Module_40, Mylar_40, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm ) );
  G4LogicalVolume *Mylar_wrap_40_log = new G4LogicalVolume( Mylar_wrap_40, GetMaterial("Mylar"), "Mylar_wrap_40_log" );

  G4Box *Mylar_38 = new G4Box( "Mylar_38", width_38/2.0 - mylar_thick, width_38/2.0 - mylar_thick, depth_38/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_38 = new G4SubtractionSolid( "Mylar_wrap_38", Module_38, Mylar_38, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm) );
  G4LogicalVolume *Mylar_wrap_38_log = new G4LogicalVolume( Mylar_wrap_38, GetMaterial("Mylar"), "Mylar_wrap_38_log" );
  */
  G4Box *Module_42_k = new G4Box( "Module_42_k", width_42/2.0, width_42/2.0, depth_42/2.0-0.25*mm );
  G4Box *Module_40_k = new G4Box( "Module_42_k", width_40/2.0, width_40/2.0, depth_40/2.0-0.25*mm );
  G4Box *Module_38_k = new G4Box( "Module_42_k", width_38/2.0, width_38/2.0, depth_38/2.0-0.25*mm );
  
  // Solid for heat conducting foils: ouside of the wrapping, supposedly
  G4Box *hc_42 = new G4Box( "hc_42", width_42/2.0 - hcf_thick, width_42/2.0 - hcf_thick, depth_42/2.0 );
  G4SubtractionSolid *hc_foil_42 = new G4SubtractionSolid( "hc_foil_42", Module_42_k, hc_42, 0, G4ThreeVector( 0, 0, 0 ) );
  G4LogicalVolume *hc_foil_42_log = new G4LogicalVolume( hc_foil_42, GetMaterial(hcf_mat_name), "hc_foil_42_log" );
  hc_foil_42_log->SetVisAttributes(hc_visAtt);

  G4Box *hc_40 = new G4Box( "hc_40", width_40/2.0 - hcf_thick, width_40/2.0 - hcf_thick, depth_40/2.0 );
  G4SubtractionSolid *hc_foil_40 = new G4SubtractionSolid( "hc_foil_40", Module_40_k, hc_40, 0, G4ThreeVector( 0, 0, 0 ) );
  G4LogicalVolume *hc_foil_40_log = new G4LogicalVolume( hc_foil_40, GetMaterial(hcf_mat_name), "hc_foil_40_log" );
  hc_foil_40_log->SetVisAttributes(hc_visAtt);
  
  G4Box *hc_38 = new G4Box( "hc_38", width_38/2.0 - hcf_thick, width_38/2.0 - hcf_thick, depth_38/2.0 );
  G4SubtractionSolid *hc_foil_38 = new G4SubtractionSolid( "hc_foil_38", Module_38_k, hc_38, 0, G4ThreeVector( 0, 0, 0 ) );
  G4LogicalVolume *hc_foil_38_log = new G4LogicalVolume( hc_foil_38, GetMaterial(hcf_mat_name), "hc_foil_38_log" );
  hc_foil_38_log->SetVisAttributes(hc_visAtt);
  
  //Next, we want to make a subtraction solid for the mylar:
  G4Box *Mylar_42 = new G4Box( "Mylar_42", width_42/2.0 - hcf_thick - mylar_thick, width_42/2.0 - hcf_thick - mylar_thick, depth_42/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_42 = new G4SubtractionSolid( "Mylar_wrap_42", hc_42, Mylar_42, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm ) );
  G4LogicalVolume *Mylar_wrap_42_log = new G4LogicalVolume( Mylar_wrap_42, GetMaterial("Mylar"), "Mylar_wrap_42_log" );
  
  G4Box *Mylar_40 = new G4Box( "Mylar_40", width_40/2.0 - hcf_thick - mylar_thick, width_40/2.0 - hcf_thick - mylar_thick, depth_40/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_40 = new G4SubtractionSolid( "Mylar_wrap_40", hc_40, Mylar_40, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm ) );
  G4LogicalVolume *Mylar_wrap_40_log = new G4LogicalVolume( Mylar_wrap_40, GetMaterial("Mylar"), "Mylar_wrap_40_log" );

  G4Box *Mylar_38 = new G4Box( "Mylar_38", width_38/2.0 - hcf_thick - mylar_thick, width_38/2.0 - hcf_thick - mylar_thick, depth_38/2.0 + 1.0*cm );
  G4SubtractionSolid *Mylar_wrap_38 = new G4SubtractionSolid( "Mylar_wrap_38", hc_38, Mylar_38, 0, G4ThreeVector( 0, 0, mylar_thick + 1.0*cm) );
  G4LogicalVolume *Mylar_wrap_38_log = new G4LogicalVolume( Mylar_wrap_38, GetMaterial("Mylar"), "Mylar_wrap_38_log" );
  
  
  ////// Define Sensitive Detector for lead-glass:
  G4String ECalTF1SDname = "Earm/ECalTF1";
  G4String ECalTF1collname = "ECalTF1HitsCollection";
  G4SBSCalSD *ECalTF1SD = NULL;
    
  if( !( ECalTF1SD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(ECalTF1SDname) ) ){
    G4cout << "Adding ECal TF1 Sensitive Detector to SDman..." << G4endl;
    ECalTF1SD = new G4SBSCalSD( ECalTF1SDname, ECalTF1collname );
    
    fDetCon->fSDman->AddNewDetector( ECalTF1SD );
    (fDetCon->SDlist).insert(ECalTF1SDname);
    fDetCon->SDtype[ECalTF1SDname] = G4SBS::kCAL;
    //fDetCon->SDarm[ECalTF1SDname] = G4SBS::kEarm;

    (ECalTF1SD->detmap).depth = 1;

    G4double default_timewindow = 100.0*ns;
    G4double default_threshold  = 0.0*MeV;

    fDetCon->SetTimeWindowAndThreshold( ECalTF1SDname, default_threshold, default_timewindow );
  }

  fDetCon->InsertSDboundaryVolume( earm_mother_log->GetName(), ECalTF1SDname );
  
  //Make lead-glass and place in modules:
  
  if( fDetCon->GetC16Segmentation() <= 0 ){
    /*
    // EPAF:2017/03/03 Was that correct anyhow ??? Don't think so...
    // block section shall be: module_section - 2*mylar_thick - 2*airthick ( - 2*hcf_thick, but I added that)
    G4Box *LeadGlass_42 = new G4Box("LeadGlass_42", width_42/2.0 - mylar_thick - air_thick, width_42/2.0 - mylar_thick - air_thick, (depth_42 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_42_log = new G4LogicalVolume( LeadGlass_42, GetMaterial("TF1"), "LeadGlass_42_log" );

    G4Box *LeadGlass_40 = new G4Box("LeadGlass_40", width_40/2.0 - mylar_thick - air_thick, width_40/2.0 - mylar_thick - air_thick, (depth_40 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_40_log = new G4LogicalVolume( LeadGlass_40, GetMaterial("TF1"), "LeadGlass_40_log" );

    G4Box *LeadGlass_38 = new G4Box("LeadGlass_38", width_38/2.0 - mylar_thick - air_thick, width_38/2.0 - mylar_thick - air_thick, (depth_38 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_38_log = new G4LogicalVolume( LeadGlass_38, GetMaterial("TF1"), "LeadGlass_38_log" );
    */
    G4Box *LeadGlass_42 = new G4Box("LeadGlass_42", width_42/2.0 - hcf_thick - mylar_thick - air_thick, width_42/2.0 - hcf_thick - mylar_thick - air_thick, (depth_42 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_42_log = new G4LogicalVolume( LeadGlass_42, GetMaterial("TF1"), "LeadGlass_42_log" );

    G4Box *LeadGlass_40 = new G4Box("LeadGlass_40", width_40/2.0 - hcf_thick - mylar_thick - air_thick, width_40/2.0 - hcf_thick - mylar_thick - air_thick, (depth_40 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_40_log = new G4LogicalVolume( LeadGlass_40, GetMaterial("TF1"), "LeadGlass_40_log" );

    G4Box *LeadGlass_38 = new G4Box("LeadGlass_38", width_38/2.0 - hcf_thick - mylar_thick - air_thick, width_38/2.0 - hcf_thick - mylar_thick - air_thick, (depth_38 - mylar_thick - air_thick)/2.0 );
    G4LogicalVolume *LeadGlass_38_log = new G4LogicalVolume( LeadGlass_38, GetMaterial("TF1"), "LeadGlass_38_log" );
    
    
    //Assign "kCAL" sensitivity to the lead-glass:
    LeadGlass_42_log->SetSensitiveDetector( ECalTF1SD );
    LeadGlass_40_log->SetSensitiveDetector( ECalTF1SD );
    LeadGlass_38_log->SetSensitiveDetector( ECalTF1SD );

    G4VisAttributes *TF1visatt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0 ) );
    LeadGlass_42_log->SetVisAttributes( TF1visatt );
    LeadGlass_40_log->SetVisAttributes( TF1visatt );
    LeadGlass_38_log->SetVisAttributes( TF1visatt );
    
    //Positioning of lead-glass in module:
    // z + Lz/2 - m/2 - a/2 = Lz/2 --> z = m/2 + a/2
    //lead-glass:
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, (mylar_thick + air_thick)/2.0 ), LeadGlass_42_log, "LeadGlass_42_phys", Module_42_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, (mylar_thick + air_thick)/2.0 ), LeadGlass_40_log, "LeadGlass_40_phys", Module_40_log, false, 0 );
    new G4PVPlacement( 0, G4ThreeVector( 0, 0, (mylar_thick + air_thick)/2.0 ), LeadGlass_38_log, "LeadGlass_38_phys", Module_38_log, false, 0 );
  } else {
    G4int nsegments = fDetCon->GetC16Segmentation();
    G4double segthick = fDetCon->GetSegmentThickC16();

    //Segment the lead-glass longitudinally so that different optical properties can be defined:
    G4int nseg_42 = G4int( (depth_42 - mylar_thick - air_thick)/segthick );
    if( nseg_42 > nsegments ){
      nseg_42 = nsegments;
    }
    G4double remainder_42 = (depth_42 - mylar_thick - air_thick) - nseg_42 * segthick; //remainder is always non-negative!

    if( nseg_42 == 0 ){
      nseg_42 = 1;
      segthick = (depth_42 - mylar_thick - air_thick);
      remainder_42 = 0.0;
    }

    G4Box *segment_42 = new G4Box( "segment_42", (width_42 - 2.*(mylar_thick+air_thick) )/2.0, (width_42 - 2.*(mylar_thick+air_thick) )/2.0, segthick/2.0 );
    G4Box *lastsegment_42 = new G4Box( "lastsegment_42", (width_42 - 2.*(mylar_thick+air_thick) )/2.0, (width_42 - 2.*(mylar_thick+air_thick) )/2.0, (segthick+remainder_42)/2.0 );

    G4double ztemp = -depth_42/2.0 + mylar_thick + air_thick;
    
    for( int seg42=0; seg42<nseg_42; seg42++ ){
      TString material_name, lv_name, pv_name;
      material_name.Form( "TF1_anneal_ECAL_z%d", seg42 );
      lv_name.Form("segment_42_log_%d", seg42 );
      pv_name.Form("segment_42_phys_%d", seg42 );
      G4LogicalVolume *Seg42_log;
      if( seg42 + 1 < nseg_42 || nseg_42 == 1 ){
	Seg42_log = new G4LogicalVolume( segment_42, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += segthick/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg42_log, pv_name.Data(), Module_42_log, false, 0 );
	ztemp += segthick/2.0;
      } else {
	Seg42_log = new G4LogicalVolume( lastsegment_42, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += (segthick+remainder_42)/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg42_log, pv_name.Data(), Module_42_log, false, 0 );
	ztemp += (segthick+remainder_42)/2.0;
      }
      Seg42_log->SetSensitiveDetector( ECalTF1SD );

      G4VisAttributes *Segment_VisAtt = new G4VisAttributes( G4Colour( 0.8*(seg42/15.0)+0.20, 0.8*(seg42/15.0)+0.20, 0.0 ) );
      Seg42_log->SetVisAttributes( Segment_VisAtt );
    }


    //Segment the lead-glass longitudinally so that different optical properties can be defined:
    G4int nseg_40 = G4int( (depth_40 - mylar_thick - air_thick)/segthick );
    if( nseg_40 > nsegments ){
      nseg_40 = nsegments;
    }
    G4double remainder_40 = (depth_40 - mylar_thick - air_thick) - nseg_40 * segthick; //remainder is always non-negative!

    if( nseg_40 == 0 ){
      nseg_40 = 1;
      segthick = (depth_40 - mylar_thick - air_thick);
      remainder_40 = 0.0;
    }

    G4Box *segment_40 = new G4Box( "segment_40", (width_40 - 2.*(mylar_thick+air_thick) )/2.0, (width_40 - 2.*(mylar_thick+air_thick) )/2.0, segthick/2.0 );
    G4Box *lastsegment_40 = new G4Box( "lastsegment_40", (width_40 - 2.*(mylar_thick+air_thick) )/2.0, (width_40 - 2.*(mylar_thick+air_thick) )/2.0, (segthick+remainder_40)/2.0 );

    ztemp = -depth_40/2.0 + mylar_thick + air_thick;
    
    for( int seg40=0; seg40<nseg_40; seg40++ ){
      TString material_name, lv_name, pv_name;
      material_name.Form( "TF1_anneal_ECAL_z%d", seg40 );
      lv_name.Form("segment_40_log_%d", seg40 );
      pv_name.Form("segment_40_phys_%d", seg40 );
      G4LogicalVolume *Seg40_log;
      if( seg40 + 1 < nseg_40 || nseg_40 == 1 ){
	Seg40_log = new G4LogicalVolume( segment_40, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += segthick/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg40_log, pv_name.Data(), Module_40_log, false, 0 );
	ztemp += segthick/2.0;
      } else {
	Seg40_log = new G4LogicalVolume( lastsegment_40, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += (segthick+remainder_40)/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg40_log, pv_name.Data(), Module_40_log, false, 0 );
	ztemp += (segthick+remainder_40)/2.0;
      }
      Seg40_log->SetSensitiveDetector( ECalTF1SD );

      G4VisAttributes *Segment_VisAtt = new G4VisAttributes( G4Colour( 0.8*(seg40/15.0)+0.20, 0.8*(seg40/15.0)+0.20, 0.0 ) );
      Seg40_log->SetVisAttributes( Segment_VisAtt );
    }

    //Segment the lead-glass longitudinally so that different optical properties can be defined:
    G4int nseg_38 = G4int( (depth_38 - mylar_thick - air_thick)/segthick );
    if( nseg_38 > nsegments ){
      nseg_38 = nsegments;
    }
    G4double remainder_38 = (depth_38 - mylar_thick - air_thick) - nseg_38 * segthick; //remainder is always non-negative!

    if( nseg_38 == 0 ){
      nseg_38 = 1;
      segthick = (depth_38 - mylar_thick - air_thick);
      remainder_38 = 0.0;
    }

    G4Box *segment_38 = new G4Box( "segment_38", (width_38 - 2.*(mylar_thick+air_thick) )/2.0, (width_38 - 2.*(mylar_thick+air_thick) )/2.0, segthick/2.0 );
    G4Box *lastsegment_38 = new G4Box( "lastsegment_38", (width_38 - 2.*(mylar_thick+air_thick) )/2.0, (width_38 - 2.*(mylar_thick+air_thick) )/2.0, (segthick+remainder_38)/2.0 );

    ztemp = -depth_38/2.0 + mylar_thick + air_thick;
    
    for( int seg38=0; seg38<nseg_38; seg38++ ){
      TString material_name, lv_name, pv_name;
      material_name.Form( "TF1_anneal_ECAL_z%d", seg38 );
      lv_name.Form("segment_38_log_%d", seg38 );
      pv_name.Form("segment_38_phys_%d", seg38 );
      G4LogicalVolume *Seg38_log;
      if( seg38 + 1 < nseg_38 || nseg_38 == 1 ){
	Seg38_log = new G4LogicalVolume( segment_38, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += segthick/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg38_log, pv_name.Data(), Module_38_log, false, 0 );
	ztemp += segthick/2.0;
      } else {
	Seg38_log = new G4LogicalVolume( lastsegment_38, GetMaterial( material_name.Data() ), lv_name.Data() );
	ztemp += (segthick+remainder_38)/2.0;
	new G4PVPlacement( 0, G4ThreeVector( 0, 0, ztemp ), Seg38_log, pv_name.Data(), Module_38_log, false, 0 );
	ztemp += (segthick+remainder_38)/2.0;
      }
      Seg38_log->SetSensitiveDetector( ECalTF1SD );

      G4VisAttributes *Segment_VisAtt = new G4VisAttributes( G4Colour( 0.8*(seg38/15.0)+0.20, 0.8*(seg38/15.0)+0.20, 0.0 ) );
      Seg38_log->SetVisAttributes( Segment_VisAtt );
    }
    
  }
    
  //Place lead-glass and mylar wrap inside module:
  
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), hc_foil_42_log, "hc_foil_42_phys", Module_42_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), hc_foil_40_log, "hc_foil_40_phys", Module_40_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), hc_foil_38_log, "hc_foil_38_phys", Module_38_log, false, 0 );
  
  //mylar:
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), Mylar_wrap_42_log, "Mylar_wrap_42_phys", Module_42_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), Mylar_wrap_40_log, "Mylar_wrap_40_phys", Module_40_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), Mylar_wrap_38_log, "Mylar_wrap_38_phys", Module_38_log, false, 0 );

  new G4LogicalSkinSurface( "Mylar_skin_42", Mylar_wrap_42_log, GetOpticalSurface("Mirrsurf") );
  new G4LogicalSkinSurface( "Mylar_skin_40", Mylar_wrap_40_log, GetOpticalSurface("Mirrsurf") );
  new G4LogicalSkinSurface( "Mylar_skin_38", Mylar_wrap_38_log, GetOpticalSurface("Mirrsurf") );

  new G4LogicalSkinSurface( "lgwrap_42", LGWrap_42_log, GetOpticalSurface("Mirrsurf") );
  new G4LogicalSkinSurface( "lgwrap_40", LGWrap_40_log, GetOpticalSurface("Mirrsurf") );
  new G4LogicalSkinSurface( "lgwrap_38", LGWrap_38_log, GetOpticalSurface("Mirrsurf") );

  

  //Filling order:
  //G4int index_super_rows[20] = {0,19,1,18,2,17,3,16,4,15,5,14,6,13,7,12,8,11,9,10};
  //Now, start placing the modules:
  //Approach: Put larger blocks on bottom, smaller blocks on top? 

  //Next: PMTs:
  G4Tubs *ecal_PMT = new G4Tubs( "ecal_PMT", 0.0, radius_ecal_pmt, depth_ecal_pmt/2.0, 0.0, twopi );
  G4LogicalVolume *ecal_PMT_log = new G4LogicalVolume( ecal_PMT, GetMaterial("Photocathode_material_ecal"), "ecal_PMT_log" );

  G4String sdname = "Earm/ECAL";
  G4String collname = "ECALHitsCollection";

  G4SBSECalSD *ECalSD = NULL;

  if( !( (G4SBSECalSD*) sdman->FindSensitiveDetector( sdname ) ) ){
    G4cout << "Adding ECAL PMT sensitive detector to SDman..." << G4endl;
    ECalSD = new G4SBSECalSD( sdname, collname );
    sdman->AddNewDetector( ECalSD );
    (fDetCon->SDlist).insert( sdname );
    fDetCon->SDtype[sdname] = G4SBS::kECAL;
    (ECalSD->detmap).depth = 0;
  }

  ecal_PMT_log->SetSensitiveDetector( ECalSD );

  fDetCon->InsertSDboundaryVolume( earm_mother_log->GetName(), sdname );
  
  int lastrow42 = 0;
  int nused42 = 0, nused40=0, nused38=0;
  int icell=0, global_row=0, global_col=0;

  int nrows42 = 0, nrows40 = 0, nrows38 = 0;

  G4double ysum = ycalo_min;

  G4VisAttributes *Alvisatt = new G4VisAttributes( G4Colour( 0.4, 0.4, 0.4 ) );

  //Fill out bottom with Al:
  G4Box *bottom_Al = new G4Box( "bottom_Al", width_earm/2.0, (ycalo_min + height_earm/2.0)/2.0, depth_leadglass/2.0 );
  G4LogicalVolume *bottom_Al_log = new G4LogicalVolume( bottom_Al, GetMaterial("Al"), "bottom_Al_log" );
  new G4PVPlacement( 0, G4ThreeVector( 0, (ycalo_min - height_earm/2.0)/2.0, depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 ), bottom_Al_log, "bottom_Al_phys", earm_mother_log, false, 0 ); 

  // bottom_Al_log->SetVisAttributes( Alvisatt );
  bottom_Al_log->SetVisAttributes( G4VisAttributes::Invisible );
  for( int super_row=0; super_row<20; super_row++ ){
    for( int sub_row=0; sub_row<4; sub_row++ ){
      global_row = sub_row + 4*super_row;
      int ncol;

      G4double xlow_row, xhigh_row;
      //check whether there are enough 4.2 cm blocks to fill this row:
      if( nused42 + ncol42_row[super_row] <= nblocks_42 ){
	//Fill this row with 4.2-cm blocks:
	ncol = ncol42_row[super_row];
	for( int col=0; col<ncol; col++ ){
	  G4ThreeVector modpos( xcalomin_super_rows[super_row] + (col + 0.5*(1+sub_row%2))*width_42,
				ysum + 0.5*width_42,
				zfront_ECAL + depth_42/2.0 );
				// depth_earm/2.0 - depth_ecal_pmt - depth_leadglass/2.0 );
	  new G4PVPlacement( 0, modpos, Module_42_log, "Module_42_phys", earm_mother_log, false, icell );

	  (ECalTF1SD->detmap).Row[icell] = global_row;
	  (ECalTF1SD->detmap).Col[icell] = col;
	  (ECalTF1SD->detmap).LocalCoord[icell] = modpos;

	  G4ThreeVector pmtpos( modpos.x(), modpos.y(), depth_earm/2.0 - depth_ecal_pmt/2.0 );
	  new G4PVPlacement( 0, pmtpos, ecal_PMT_log, "ecal_PMT_phys", earm_mother_log, false, icell );

	  (ECalSD->detmap).Row[icell] = global_row;
	  (ECalSD->detmap).Col[icell] = col;
	  (ECalSD->detmap).LocalCoord[icell] = pmtpos;

	  //Add light-guide with mylar wrap:
	  G4ThreeVector LGpos( modpos.x(), modpos.y(), modpos.z() + depth_42/2.0 + LG42->GetZHalfLength() );
	  new G4PVPlacement( 0, LGpos, LG42_log, "LG42_phys", earm_mother_log, false, icell );
	  // new G4PVPlacement( 0, LGpos, LightGuide_42_log, "LightGuide_42_phys", earm_mother_log, false, icell );
	  // new G4PVPlacement( 0, LGpos, LGWrap_42_log, "LGWrap_42_phys", earm_mother_log, false, icell );

	  // //EPAF 2017-01-12: Need to make sensitive the three volumes above, to measure their dose.
	  // // shall be temporary, and not end in the repo...
	  // (ECalLGSD->detmap).Row[icell] = global_row;
	  // (ECalLGSD->detmap).Col[icell] = col;
	  // (ECalLGSD->detmap).LocalCoord[icell] = modpos;
	  
	  if( col == 0 ) xlow_row = modpos.x() - 0.5*width_42;
	  if( col+1 == ncol ) xhigh_row = modpos.x() + 0.5*width_42;
	  
	  icell++;
	}

	ysum += width_42;
	
	nused42 += ncol42_row[super_row];
	nrows42++;
	
	char prefix[255];
	
	//Fill out row with Al:
	if( xlow_row > -width_earm/2.0 ){
	  sprintf( prefix, "Albox1_row%d", super_row );
	  G4String boxname = prefix;
	  
	  G4Box *Al_filler1 = new G4Box( boxname, (xlow_row + width_earm/2.0)/2.0, width_42/2.0, depth_leadglass/2.0 );
	  
	  G4String logname = boxname + "log";
	  
	  G4LogicalVolume *Al_filler1_log = new G4LogicalVolume( Al_filler1, GetMaterial("Al"), logname );
	  
	  //Al_filler1_log->SetVisAttributes( Alvisatt );
	  Al_filler1_log->SetVisAttributes( G4VisAttributes::Invisible );

	  G4ThreeVector pos( (xlow_row - width_earm/2.0)/2.0,
			     ysum - 0.5*width_42,
			     depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 );

	  G4String physname = boxname + "phys";
	  
	  new G4PVPlacement( 0, pos, Al_filler1_log, physname, earm_mother_log, false, 0 );
	}
	if( xhigh_row < width_earm/2.0 ){
	  sprintf( prefix, "Albox2_row%d", super_row );
	  G4String boxname = prefix;
	  
	  G4Box *Al_filler2 = new G4Box( boxname, (-xhigh_row + width_earm/2.0)/2.0, width_42/2.0, depth_leadglass/2.0 );
	  
	  G4String logname = boxname + "log";
	  
	  G4LogicalVolume *Al_filler2_log = new G4LogicalVolume( Al_filler2, GetMaterial("Al"), logname );
	  //Al_filler2_log->SetVisAttributes( Alvisatt );

	  Al_filler2_log->SetVisAttributes( G4VisAttributes::Invisible );
	  
	  //G4ThreeVector pos( (xhigh_row + width_earm/2.0)/2.0, modpos.y(), modpos.z() );
	  G4ThreeVector pos( (xhigh_row + width_earm/2.0)/2.0,
			     ysum - 0.5*width_42,
			     depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 );
	  
	  G4String physname = boxname + "phys";
	  
	  new G4PVPlacement( 0, pos, Al_filler2_log, physname, earm_mother_log, false, 0 );
	}
	
      } else if( nused40 + ncol40_row[super_row] <= nblocks_40 ){
	ncol = ncol40_row[super_row];
	for( int col=0; col<ncol; col++ ){
	  G4ThreeVector modpos( xcalomin_super_rows[super_row] + (col + 0.5*(1+sub_row%2))*width_40,
				ysum + 0.5*width_40,
				zfront_ECAL + depth_40/2.0 );
				//depth_earm/2.0 - depth_ecal_pmt - depth_leadglass/2.0 );
	  new G4PVPlacement( 0, modpos, Module_40_log, "Module_40_phys", earm_mother_log, false, icell );

	  (ECalTF1SD->detmap).Row[icell] = global_row;
	  (ECalTF1SD->detmap).Col[icell] = col;
	  (ECalTF1SD->detmap).LocalCoord[icell] = modpos;

	  G4ThreeVector pmtpos( modpos.x(), modpos.y(), depth_earm/2.0 - depth_ecal_pmt/2.0 );
	  new G4PVPlacement( 0, pmtpos, ecal_PMT_log, "ecal_PMT_phys", earm_mother_log, false, icell );

	  (ECalSD->detmap).Row[icell] = global_row;
	  (ECalSD->detmap).Col[icell] = col;
	  (ECalSD->detmap).LocalCoord[icell] = pmtpos;

	  //Add light-guide with mylar wrap:
	  G4ThreeVector LGpos( modpos.x(), modpos.y(), modpos.z() + depth_40/2.0 + LightGuide_40->GetZHalfLength() );
	  new G4PVPlacement( 0, LGpos, LG40_log, "LG40_phys", earm_mother_log, false, icell );
	  // new G4PVPlacement( 0, LGpos, LightGuide_40_log, "LightGuide_40_phys", earm_mother_log, false, icell );
	  // new G4PVPlacement( 0, LGpos, LGWrap_40_log, "LGWrap_40_phys", earm_mother_log, false, icell );
	  
	  //EPAF 2017-01-12: Need to make sensitive the three volumes above, to measure their dose.
	  // shall be temporary, and not end in the repo...
	  // (ECalLGSD->detmap).Row[icell] = global_row;
	  // (ECalLGSD->detmap).Col[icell] = col;
	  // (ECalLGSD->detmap).LocalCoord[icell] = modpos;
	  
	  if( col == 0 ) xlow_row = modpos.x() - 0.5*width_40;
	  if( col+1 == ncol ) xhigh_row = modpos.x() + 0.5*width_40;
	  
	  icell++;
	}
	ysum += width_40;
	nused40 += ncol40_row[super_row];
	nrows40++;

	char prefix[255];
	//Fill out row with Al:
	if( xlow_row > -width_earm/2.0 ){
	  sprintf( prefix, "Albox1_row%d", super_row );
	  G4String boxname = prefix;
	  
	  G4Box *Al_filler1 = new G4Box( boxname, (xlow_row + width_earm/2.0)/2.0, width_40/2.0, depth_leadglass/2.0 );
	  
	  G4String logname = boxname + "log";
	  
	  G4LogicalVolume *Al_filler1_log = new G4LogicalVolume( Al_filler1, GetMaterial("Al"), logname );
	  //Al_filler1_log->SetVisAttributes( Alvisatt );

	  Al_filler1_log->SetVisAttributes( G4VisAttributes::Invisible );
	  
	  G4ThreeVector pos( (xlow_row - width_earm/2.0)/2.0,
			     ysum - 0.5*width_40,
			     depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 );

	  G4String physname = boxname + "phys";
	  
	  new G4PVPlacement( 0, pos, Al_filler1_log, physname, earm_mother_log, false, 0 );
	}
	if( xhigh_row < width_earm/2.0 ){
	  sprintf( prefix, "Albox2_row%d", super_row );
	  G4String boxname = prefix;
	  
	  G4Box *Al_filler2 = new G4Box( boxname, (-xhigh_row + width_earm/2.0)/2.0, width_40/2.0, depth_leadglass/2.0 );
	  
	  G4String logname = boxname + "log";
	  
	  G4LogicalVolume *Al_filler2_log = new G4LogicalVolume( Al_filler2, GetMaterial("Al"), logname );
	  //Al_filler2_log->SetVisAttributes( Alvisatt );

	  Al_filler2_log->SetVisAttributes( G4VisAttributes::Invisible );
	  
	  //G4ThreeVector pos( (xhigh_row + width_earm/2.0)/2.0, modpos.y(), modpos.z() );
	  G4ThreeVector pos( (xhigh_row + width_earm/2.0)/2.0,
			     ysum - 0.5*width_40,
			     depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 );
	  
	  G4String physname = boxname + "phys";
	  
	  new G4PVPlacement( 0, pos, Al_filler2_log, physname, earm_mother_log, false, 0 );
	}
	
      } else if( nused38 + ncol38_row[super_row] <= nblocks_38 ){
	ncol = ncol38_row[super_row];
	for( int col=0; col<ncol; col++ ){
	  G4ThreeVector modpos( xcalomin_super_rows[super_row] + (col + 0.5*(1+sub_row%2))*width_38,
				ysum + 0.5*width_38,
				zfront_ECAL + depth_38/2.0 );
				//depth_earm/2.0 - depth_ecal_pmt - depth_leadglass/2.0 );
	  new G4PVPlacement( 0, modpos, Module_38_log, "Module_38_phys", earm_mother_log, false, icell );

	  (ECalTF1SD->detmap).Row[icell] = global_row;
	  (ECalTF1SD->detmap).Col[icell] = col;
	  (ECalTF1SD->detmap).LocalCoord[icell] = modpos;

	  G4ThreeVector pmtpos( modpos.x(), modpos.y(), depth_earm/2.0 - depth_ecal_pmt/2.0 );
	  new G4PVPlacement( 0, pmtpos, ecal_PMT_log, "ecal_PMT_phys", earm_mother_log, false, icell );

	  (ECalSD->detmap).Row[icell] = global_row;
	  (ECalSD->detmap).Col[icell] = col;
	  (ECalSD->detmap).LocalCoord[icell] = pmtpos;

	   //Add light-guide with mylar wrap:
	  G4ThreeVector LGpos( modpos.x(), modpos.y(), modpos.z() + depth_38/2.0 + LightGuide_38->GetZHalfLength() );
	  new G4PVPlacement( 0, LGpos, LG38_log, "LG38_phys", earm_mother_log, false, icell );
	  // new G4PVPlacement( 0, LGpos, LightGuide_38_log, "LightGuide_38_phys", earm_mother_log, false, icell );
	  // new G4PVPlacement( 0, LGpos, LGWrap_38_log, "LGWrap_38_phys", earm_mother_log, false, icell );

	  //EPAF 2017-01-12: Need to make sensitive the three volumes above, to measure their dose.
	  // shall be temporary, and not end in the repo...
	  // (ECalLGSD->detmap).Row[icell] = global_row;
	  // (ECalLGSD->detmap).Col[icell] = col;
	  // (ECalLGSD->detmap).LocalCoord[icell] = modpos;
	  
	  if( col == 0 ) xlow_row = modpos.x() - 0.5*width_38;
	  if( col+1 == ncol ) xhigh_row = modpos.x() + 0.5*width_38;
	  
	  icell++;
	}
	ysum += width_38;
	nused38 += ncol38_row[super_row];
	nrows38++;

	char prefix[255];
	//Fill out row with Al:
	if( xlow_row > -width_earm/2.0 ){
	  sprintf( prefix, "Albox1_row%d", super_row );
	  G4String boxname = prefix;
	  
	  G4Box *Al_filler1 = new G4Box( boxname, (xlow_row + width_earm/2.0)/2.0, width_38/2.0, depth_leadglass/2.0 );
	  
	  G4String logname = boxname + "log";
	  
	  G4LogicalVolume *Al_filler1_log = new G4LogicalVolume( Al_filler1, GetMaterial("Al"), logname );
	  //Al_filler1_log->SetVisAttributes( Alvisatt );

	  Al_filler1_log->SetVisAttributes( G4VisAttributes::Invisible );
	  
	  G4ThreeVector pos( (xlow_row - width_earm/2.0)/2.0,
			     ysum - 0.5*width_38,
			     depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 );

	  G4String physname = boxname + "phys";
	  
	  new G4PVPlacement( 0, pos, Al_filler1_log, physname, earm_mother_log, false, 0 );
	}
	if( xhigh_row < width_earm/2.0 ){
	  sprintf( prefix, "Albox2_row%d", super_row );
	  G4String boxname = prefix;
	  
	  G4Box *Al_filler2 = new G4Box( boxname, (-xhigh_row + width_earm/2.0)/2.0, width_38/2.0, depth_leadglass/2.0 );
	  
	  G4String logname = boxname + "log";
	  
	  G4LogicalVolume *Al_filler2_log = new G4LogicalVolume( Al_filler2, GetMaterial("Al"), logname );
	  //Al_filler2_log->SetVisAttributes( Alvisatt );

	  Al_filler2_log->SetVisAttributes( G4VisAttributes::Invisible );
	  
	  //G4ThreeVector pos( (xhigh_row + width_earm/2.0)/2.0, modpos.y(), modpos.z() );
	  G4ThreeVector pos( (xhigh_row + width_earm/2.0)/2.0,
			     ysum - 0.5*width_38,
			     depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 );
	  
	  G4String physname = boxname + "phys";
	  
	  new G4PVPlacement( 0, pos, Al_filler2_log, physname, earm_mother_log, false, 0 );
	}
	
      }
    }
  }

  //Fill out top with Al:
  G4Box *top_Al = new G4Box( "top_Al", width_earm/2.0, (-ysum + height_earm/2.0)/2.0, depth_leadglass/2.0 );
  G4LogicalVolume *top_Al_log = new G4LogicalVolume( top_Al, GetMaterial("Al"), "top_Al_log" );
  new G4PVPlacement( 0, G4ThreeVector( 0, (ysum + height_earm/2.0)/2.0, depth_earm/2.0 - depth_ecal_pmt - depth_lightguide_short - depth_leadglass/2.0 ), top_Al_log, "top_Al_phys", earm_mother_log, false, 0 ); 

  //  top_Al_log->SetVisAttributes( Alvisatt );
  top_Al_log->SetVisAttributes( G4VisAttributes::Invisible );
  
  //Next: Put front Aluminum plate in front of ECAL (make wireframe):
  G4Box *ECAL_FrontPlate = new G4Box( "ECAL_FrontPlate", width_earm/2.0, height_earm/2.0, depth_ecal_frontplate/2.0 );
  G4LogicalVolume *ECAL_FrontPlate_log = new G4LogicalVolume( ECAL_FrontPlate, GetMaterial("Al"), "ECAL_FrontPlate_log" );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, zfront_ECAL - depth_ecal_frontplate/2.0 ), ECAL_FrontPlate_log, "ECAL_FrontPlate_phys", earm_mother_log, false, 0 );

  //Next: CH2 filter:
  G4Box *CH2_filter = new G4Box( "CH2_filter", width_earm/2.0, height_earm/2.0, depth_CH2/2.0 );
  G4LogicalVolume *CH2_filter_log = new G4LogicalVolume( CH2_filter, GetMaterial("Polyethylene"), "CH2_filter_log" );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, zfront_ECAL - depth_CDET - depth_CH2/2.0 ), CH2_filter_log, "CH2_filter_phys", earm_mother_log, false, 0 );

  G4double z0_CDET = -depth_earm/2.0 + depth_CH2;
  //G4double R0_CDET = R_Earm - depth_leadglass - depth_CDET;
  //G4double R0_CDET = fBBdist - depth_leadglass - depth_CDET;
  G4double R0_CDET = fBBdist - depth_CDET;
  
  MakeCDET( R0_CDET, z0_CDET, earm_mother_log );

  //Visualization:
  
  G4VisAttributes *FrontPlate_visatt = new G4VisAttributes( G4Colour( 0.7, 0.7, 0.7 ) );
  FrontPlate_visatt->SetForceWireframe(true);
  ECAL_FrontPlate_log->SetVisAttributes( FrontPlate_visatt );

  G4VisAttributes *CH2_visatt = new G4VisAttributes( G4Colour( 0, 0.6, 0.6 ) );
  CH2_visatt->SetForceWireframe(true);

  CH2_filter_log->SetVisAttributes(CH2_visatt);
  
  earm_mother_log->SetVisAttributes( G4VisAttributes::Invisible );
  Module_42_log->SetVisAttributes( G4VisAttributes::Invisible );
  Module_40_log->SetVisAttributes( G4VisAttributes::Invisible );
  Module_38_log->SetVisAttributes( G4VisAttributes::Invisible );

  G4VisAttributes *Mylarvisatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
  Mylarvisatt->SetForceWireframe(true);
  Mylar_wrap_42_log->SetVisAttributes( Mylarvisatt );
  Mylar_wrap_40_log->SetVisAttributes( Mylarvisatt );
  Mylar_wrap_38_log->SetVisAttributes( Mylarvisatt );

  G4VisAttributes *ECALpmtvisatt = new G4VisAttributes( G4Colour( 0, 0, 1 ) );
  ecal_PMT_log->SetVisAttributes( ECALpmtvisatt );
  
}

void G4SBSEArmBuilder::MakeCDET( G4double R0, G4double z0, G4LogicalVolume *mother ){
  //R0 is the nominal distance from target to the start of CDET
  //z0 is the z position of the start of CDET relative to mother
  
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
    fDetCon->SDtype[sdname] = G4SBS::kECAL;
    (cdet_sd->detmap).depth = 0;
    CDET_pmt_cathode_log->SetSensitiveDetector( cdet_sd );
  }

  fDetCon->InsertSDboundaryVolume( mother->GetName(), sdname );
  
  sdname = "Earm/CDET_Scint";
  collname = "CDET_ScintHitsCollection";

  G4SBSCalSD *cdet_scint_sd = NULL;
  
  if( !( cdet_scint_sd = (G4SBSCalSD*) sdman->FindSensitiveDetector( sdname ) ) ){
    G4cout << "Adding CDET Scint sensitive detector to sdman..." << G4endl;
    cdet_scint_sd = new G4SBSCalSD( sdname, collname );
    
    fDetCon->fSDman->AddNewDetector( cdet_scint_sd );
    (fDetCon->SDlist).insert( sdname );
    fDetCon->SDtype[sdname] = G4SBS::kCAL;
    (cdet_scint_sd->detmap).depth = 1;
    ScintStripLog->SetSensitiveDetector( cdet_scint_sd );

    G4double default_timewindow = 50.0*ns;
    G4double default_threshold  = 4.0*MeV;

    fDetCon->SetTimeWindowAndThreshold( sdname, default_threshold, default_timewindow );
  }

  fDetCon->InsertSDboundaryVolume( mother->GetName(), sdname );
  
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
  G4double R0_planes[2] = { R0 + Lz_scint/2.0 + 1.0*cm,
			    R0 + 3.0*Lz_scint/2.0 + 2.0*cm }; //allow for some small (1 cm) gaps between CH2 and start of 1st plane and between 1st and second planes...

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
	G4ThreeVector pos_strip( x0_modules[imod] + ( col - 0.5 )*(Lx_scint+mylar_thick), R0_planes[plane] * tan( alpha ), z0 + R0_planes[plane] - R0 );

	G4RotationMatrix *rot_strip = new G4RotationMatrix;
	rot_strip->rotateY( col*pi );
	rot_strip->rotateX( alpha*pow(-1,col) );
	//rot_strip->rotateX( alpha );
	char physname[255];
	sprintf( physname, "Scint_module_phys_plane%d_row%d_col%d", plane+1, row+295, col+1 ); //In this construction, row varies from 295-588 for the top half
	G4String pvolname = physname;
	new G4PVPlacement( rot_strip, pos_strip, Scint_module, pvolname, mother, false, istrip );

	G4ThreeVector pos_pmt( x0_modules[imod] + pow(-1,col+1)*(Lx_scint+mylar_thick+0.1*cm), R0_planes[plane] * tan( alpha ), z0 + R0_planes[plane] - R0 );
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

// void G4SBSEArmBuilder::MakeBigCal(G4LogicalVolume *worldlog){

//   printf("BigCal at %f deg\n", fBBang/deg);

//   // Pointer to SDmanager, used frequently in this routine
//   G4SDManager *sdman = fDetCon->fSDman;

//   // Electron Arm - order of materials - all "right next" to each other in z-hat_spectrometer coordinates
//   //   1) 20cm Polyethylene 
//   //   2) 2 4cm CDet planes
//   //   3) 2.4cm Al
//   //   4) ECal

//   // Rotations, offsets, and distances from Target:
//   G4RotationMatrix *bbrm = new G4RotationMatrix;
//   bbrm->rotateY(fBBang);

//   //   -some of these variables are used in ECAL Section exclusively.
//   G4double x_ecal = 246.402*cm, y_ecal = 370.656*cm, z_ecal = 45.406*cm; //45.006cm(module) + 2*0.2(pmt depth)

//   G4double z_Al = 2.40*cm;         //Al depth
//   G4double polydepth = 20.00*cm;   //Polyethylene
//   G4double CD_depth = 4.00*cm;     //Used in Option 1
//   G4double x_earm = 200.070*cm, y_earm = 315.900*cm, z_earm = 30.40*cm; 

//   double bbr = fBBdist - z_ecal/2.0; //backside of ECal should be located at fBBdist
//   double offset = 15*cm;             //Motivation - match SBS acceptance, used to declare distances
//                                      //from target and frequently in ECAL section


//   /*************************************************************************************
//    ********************************        CDet      ***********************************
//    *************************************************************************************/

//   // 4 Options/Configurations available:
//   //    /g4sbs/CDetconfig 1   -   Just boxes containing nothing but one material
//   //    /g4sbs/CDetconfig 2   -   392 bars make up one module, 3 modules per plane. 
//   //                              PMTs/Scint are SD, planes are FLAT
//   //    /g4sbs/CDetconfig 3   -   Same as 2, except Top/Bottom modules are rotated 
//   //                              by some angle CDetPhi such that a line
//   //                              running from the origin to the center of a module 
//   //                              is perpendular to the face of a module
//   //    /g4sbs/CDetconfig 4   -   All bars are rotated such that a line running from 
//   //                              the origin to the  center of a "long" bar
//   //                              is perpendicular to the face of a CDet bar. A "long" 
//   //                              bar is defined as two small bars 
//   //                              connected long ways, or one CDet row.

//   G4double x_bar_CDet = 0.5*cm;
//   G4double y_bar_CDet = 4.0*cm;
//   G4double z_bar_CDet = 51*cm;

//   G4double mylar_CDet       = 0.025*cm;
//   G4double pmt_CDet_depth   = 0.010*cm;

//   G4double BoreHoleDiameter = 0.30*cm;
//   G4double WLSDiameter      = 0.20*cm; 

//   G4int CDetColumns = 2;
//   G4int CDetRows    = 196;

//   // Setting up the Scintillating Bars by defining a mother volume
//   G4Box *CDetBox = new G4Box("CDetBox",(x_bar_CDet+2*mylar_CDet)/2.0,
// 			     (y_bar_CDet+2*mylar_CDet)/2.0,
// 			     (z_bar_CDet+mylar_CDet)/2.0 );
//   G4LogicalVolume *CDetLog = new G4LogicalVolume( CDetBox, GetMaterial("Special_Air"), "CDetLog" );
  
//   G4Box *CDetScintBox_temp = new G4Box("CDetScintBox_temp", x_bar_CDet/2.0, y_bar_CDet/2.0, z_bar_CDet/2.0);

//   // Setting up Mylar Wrap
//   G4SubtractionSolid *CDetMylarBox = new G4SubtractionSolid( "CDetMylar", CDetBox, CDetScintBox_temp, 0, G4ThreeVector(0.0,0.0,mylar_CDet/2.0) );
//   G4LogicalVolume *CDetMylarLog = new G4LogicalVolume( CDetMylarBox, GetMaterial("Mylar"), "CDetMylarLog" );
//   new G4LogicalSkinSurface( "CDet_Mylar_Skin", CDetMylarLog, GetOpticalSurface("Mirrsurf") );
    
//   // Bore out cylindrical section of our scintillating material, CDetScintBox.
//   G4Tubs *BoreHole = new G4Tubs( "BoreHole", 0.0, BoreHoleDiameter/2.0, 1.001*z_bar_CDet/2.0, 0.0, twopi );
//   G4SubtractionSolid *CDScintBoreBox = new G4SubtractionSolid( "CDScintBoreBox", CDetScintBox_temp, BoreHole, 0, G4ThreeVector() );
//   G4LogicalVolume *CDScintBoreLog = new G4LogicalVolume( CDScintBoreBox, GetMaterial("EJ232"), "CDScintBoreLog" );

//   // Make a Lightguide
//   G4Tubs *WLSguide = new G4Tubs( "WLSguide", 0.0, WLSDiameter/2.0, z_bar_CDet/2.0, 0.0, twopi );
//   G4LogicalVolume *WLSguideLog = new G4LogicalVolume( WLSguide, GetMaterial("BC484"), "WLSguideLog" );

//   // Place these materials inside CDetBox.
//   new G4PVPlacement(0, G4ThreeVector(), CDetMylarLog, "CDetMylarLogPhys", CDetLog, false, 0); 
//   new G4PVPlacement(0, G4ThreeVector(0.0,0.0,mylar_CDet/2.0), CDScintBoreLog, "CDScintBorePhys", CDetLog, false, 0);
//   new G4PVPlacement(0, G4ThreeVector(0.0,0.0,mylar_CDet/2.0), WLSguideLog, "WLSguidePhys", CDetLog, false, 0);

//   // Make some PMTs
//   G4Tubs *CDetPMT = new G4Tubs( "CDetPMT", 0.0, WLSDiameter/2.0, pmt_CDet_depth/2.0, 0.0, twopi );
//   G4LogicalVolume *CDetPMTLog = new G4LogicalVolume( CDetPMT, GetMaterial("Photocathode_CDet"), "CDetPMTLog" );

//   // CDetPMTLog is the SD, assigned to ECalSD which detects optical photons
//   G4String CDetSDname = "Earm/CDet";
//   G4String CDetcollname = "CDetHitsCollection";
//   G4SBSECalSD *CDetSD = NULL;

//   if( !( (G4SBSECalSD*) sdman->FindSensitiveDetector(CDetSDname)) ) {
//     G4cout << "Adding CDet PMT Sensitive Detector to SDman... " << G4endl;
//     CDetSD = new G4SBSECalSD( CDetSDname, CDetcollname );
//     sdman->AddNewDetector( CDetSD );
//     (fDetCon->SDlist).insert( CDetSDname );
//     fDetCon->SDtype[CDetSDname] = G4SBS::kECAL;
//     //fDetCon->SDarm[CDetSDname] = G4SBS::kEarm;
//     //(CDetSD->detmap).depth = 1;            // needs to be updated!!!!! and changed for each option most likely
//   }
//   CDetPMTLog->SetSensitiveDetector( CDetSD );
    
//   // Make the Final CDet Bar, housing CDetLog and CDetPMTLog!
//   G4double CDetBarX = x_bar_CDet + 2*mylar_CDet;
//   G4double CDetBarY = y_bar_CDet + 2*mylar_CDet;
//   G4double CDetBarZ = z_bar_CDet + mylar_CDet + pmt_CDet_depth;
//   G4Box *CDetBarBox = new G4Box( "CDetBarBox", CDetBarX/2.0, CDetBarY/2.0, CDetBarZ/2.0 );
//   G4LogicalVolume *CDetBarLog = new G4LogicalVolume( CDetBarBox, GetMaterial("Air"), "CDetBarLog" );

//   // Placing..
//   new G4PVPlacement(0, G4ThreeVector(0.0,0.0,(CDetBarZ-pmt_CDet_depth)/2.0), CDetPMTLog,"CDetPMTPhys", CDetBarLog, false, 0);
//   new G4PVPlacement(0, G4ThreeVector(0.0,0.0,-pmt_CDet_depth/2.0), CDetLog, "CDetPhys", CDetBarLog, false, 0);

//   G4RotationMatrix *CDetBarRot = new G4RotationMatrix;
//   CDetBarRot->rotateZ(90*deg);
//   //CDetBarRot->rotateX(90*deg);
//   //TEST to see if module looks like I expect
//   //new G4PVPlacement(CDetBarRot,G4ThreeVector(-5*m, 0*m, 3.0*m), CDetBarLog, "Ex", worldlog, false, 0,true );

//   // Lets define a CDet MODULE - 2 cols, 196 rows => 392 bars in 1 module:
//   G4double CDetModuleX = CDetColumns*CDetBarZ; 
//   G4double CDetModuleY = CDetRows*CDetBarX;
//   G4double CDetModuleZ = CDetBarY;
//   G4Box *CDetModuleBox = new G4Box( "CDetModuleBox", CDetModuleX/2.0, CDetModuleY/2.0, CDetModuleZ/2.0 );
//   G4LogicalVolume *CDetModuleLog = new G4LogicalVolume( CDetModuleBox, GetMaterial("Air"), "CDetModuleLog" );

//   // Place Bars in a Module, rotating based on Column. 
//   G4RotationMatrix *CDetBar_Col1 = new G4RotationMatrix;
//   CDetBar_Col1->rotateZ(90.0*deg);
//   CDetBar_Col1->rotateX(90.0*deg);
//   G4RotationMatrix *CDetBar_Col2 = new G4RotationMatrix;
//   CDetBar_Col2->rotateZ(-90.0*deg);
//   CDetBar_Col2->rotateX(90.0*deg);
//   G4int cdet_copy_number = 0;
  
//   // Make modules for Options 2 & 3
//   if( fDetCon->GetCDetConfigOption() == 2 || fDetCon->GetCDetConfigOption() == 3 ) {
//     (CDetSD->detmap).depth = 1; 
//     for(int l=0; l<CDetColumns; l++) {
//       for(int j=0; j<CDetRows; j++) {
// 	double xtemp = (CDetModuleX - CDetBarZ)/2.0 - l*CDetBarZ;
// 	double ytemp = (CDetModuleY - CDetBarX)/2.0 - j*CDetBarX;

// 	(CDetSD->detmap).Col[cdet_copy_number] = l;
// 	(CDetSD->detmap).Row[cdet_copy_number] = j;
// 	if(l==0) {
// 	  new G4PVPlacement( CDetBar_Col1, G4ThreeVector(xtemp,ytemp,0.0), CDetBarLog, "CDetBarPhys", CDetModuleLog, false, cdet_copy_number );
// 	  (CDetSD->detmap).LocalCoord[cdet_copy_number] = G4ThreeVector(xtemp+CDetBarZ/2.0,ytemp,0.0);
// 	  cdet_copy_number++;
// 	}
// 	if(l==1) {
// 	  new G4PVPlacement( CDetBar_Col2, G4ThreeVector(xtemp,ytemp,0.0), CDetBarLog, "CDetBarPhys", CDetModuleLog, false, cdet_copy_number );
// 	  (CDetSD->detmap).LocalCoord[cdet_copy_number] = G4ThreeVector(xtemp-CDetBarZ/2.0,ytemp,0.0);
// 	  cdet_copy_number++;
// 	}
//       }
//     }
//   }
//   // 3 Modules make up a "Plane" and CDet consists of two planes:
//   G4Box *AlBox = new G4Box( "AlBox", x_earm/2.0, y_earm/2.0, z_Al/2.0 );
//   G4LogicalVolume *Allog = new G4LogicalVolume ( AlBox, GetMaterial("RICHAluminum"), "Allog" );

//   G4double PlaneX = x_earm;
//   G4double PlaneY = 3*CDetModuleY; //no 2
//   G4double PlaneZ = CDetModuleZ + CDetModuleY*sin(40*deg); //sin25

//   // Mother Volume - houses 2 planes, Al, Polyethylene
//   G4Box *PlaneBox = new G4Box( "PlaneBox", PlaneX/2.0, PlaneY/2.0, PlaneZ/2.0 );
//   G4LogicalVolume *PlaneLog = new G4LogicalVolume( PlaneBox, GetMaterial("Air"), "PlaneLog" );

//   // Initialize some Geometry Variables:
//   double Layeroffset = 0.0;                // Space inbetween planes, needed for Options 3 & 4
//   double TopBottomModuleOffset = 0.0;      // Offset top/bottom module
//   double MiddleModuleOffset    = 0.0;      // Offset middle module
//   double CDetPhi = 0.0;                    // Angle that rotates bars or modules
//   double CDetFrontRad = 0.0;               // Distance between origin & front of CDet - used for rotations of modules/bars                          
//   double cdr = 0.0;                        // Distance between origin and center of CDet - used to place mother

//   // Make the modules for Option 4, redefine some variables due to rotations
//   double CDetBarX4 = 2*CDetBarZ; 
//   double CDetBarY4 = CDetBarX;
//   double CDetBarZ4 = CDetBarY;
//   G4Box *CDetModuleBox4 = new G4Box( "CDetModuleBox4", CDetBarX4/2.0, CDetBarY4/2.0, CDetBarZ4/2.0 );
//   G4LogicalVolume *CDetModuleLog4 = new G4LogicalVolume( CDetModuleBox4, GetMaterial("Air"), "CDetModuleLog4");
//   // Place two bars in one of these modules
//   new G4PVPlacement(CDetBar_Col1, G4ThreeVector((CDetBarX4-CDetBarZ)/2.0,0.0,0.0), CDetBarLog, "CDetBarPhys4", CDetModuleLog4, false, 0);
//   new G4PVPlacement(CDetBar_Col2, G4ThreeVector((-CDetBarX4+CDetBarZ)/2.0,0.0,0.0), CDetBarLog, "CDetBarPhys4", CDetModuleLog4, false, 0);

//   int CDetOpt4_copyn = 0;

//   if( fDetCon->GetCDetConfigOption() == 4 ) {
//     (CDetSD->detmap).depth = 2;  
//     cdr = fBBdist - z_ecal - PlaneZ/2.0;  
//     CDetFrontRad = cdr + PlaneZ/2.0 - CDetBarZ4/2.0;
//     CDetPhi = 2*atan( CDetBarY4/(2*CDetFrontRad) );

//     // start at center, work up/down:
//     for(int j=0; j<(1.5*CDetRows); j++) { 

//       G4RotationMatrix *CDetBar4Rottop = new G4RotationMatrix;
//       double angle = j*CDetPhi;
//       CDetBar4Rottop->rotateX(angle);

//       double dispy = ((CDetBarZ4/2.0)*sin( angle ) - (CDetBarY4/2.0)*(1-cos( angle ))); 
//       double dispz = ((CDetBarZ4/2.0)*(1-cos( angle )) + (CDetBarY4/2.0)*sin( angle ));
//       double ytemp = j*CDetBarX;
//       double ztemp = PlaneZ/2.0 - z_Al - CDetBarZ4/2.0;

//       // Additional placement correction since modules are placed relative to top left (TL) corner of previous module:
//       // currently this is just an approximation
//       double TLZ = 0.565*j*CDetBarY4*sin(angle)*cos(angle);
//       double TLY = 0.1059*j*CDetBarY4*sin(angle)*sin(angle);
//       new G4PVPlacement(CDetBar4Rottop, G4ThreeVector(0.0,ytemp+dispy-TLY,-dispz-TLZ+ztemp), CDetModuleLog4, "CDetBarPhys4",PlaneLog, false, CDetOpt4_copyn,true);
//       // odd rows including 0
//       (CDetSD->detmap).LocalCoord[CDetOpt4_copyn] = G4ThreeVector(0.0,ytemp+dispy-TLY,-dispz-TLZ+ztemp);  
//       (CDetSD->detmap).Row[CDetOpt4_copyn] = j;                                                            
//       CDetOpt4_copyn++;

//       G4RotationMatrix *CDetBar4Rotbot = new G4RotationMatrix;
//       angle = -j*CDetPhi;
//       CDetBar4Rotbot->rotateX(angle);
//       if(j>0) {
// 	new G4PVPlacement(CDetBar4Rotbot, G4ThreeVector(0.0,-(ytemp+dispy-TLY),-dispz-TLZ+ztemp), CDetModuleLog4, "CDetBarPhys4",PlaneLog, false, CDetOpt4_copyn,true);
// 	// even rows excluding 0
// 	(CDetSD->detmap).LocalCoord[CDetOpt4_copyn] = G4ThreeVector(0.0,-(ytemp+dispy-TLY),-dispz-TLZ+ztemp);  
// 	(CDetSD->detmap).Row[CDetOpt4_copyn] = j;
// 	CDetOpt4_copyn++;
//       }
//     }

//     // Place Option 4
//     G4ThreeVector CDetPos( cdr*sin(-fBBang)+(offset+ (4.212/2.0)*cm )*cos(fBBang), 
// 			   -(4.212/2.0)*cm, 
// 			   cdr*cos(-fBBang)+(offset+(4.212/2.0)*cm)*sin(fBBang) );
    
//     new G4PVPlacement( bbrm, CDetPos, PlaneLog, "CDet", worldlog, false, 0); \
//     new G4PVPlacement( 0, G4ThreeVector(0.0,0.0,(PlaneZ-z_Al)/2.0),Allog,"AlOpt4",PlaneLog,false,0);
//   }

//   // Place Options 2 & 3
//   if( fDetCon->GetCDetConfigOption() == 2 || fDetCon->GetCDetConfigOption() == 3 ) {
//     cdr = fBBdist - z_ecal - PlaneZ/2.0;                   //distance from origin to center of CDet Plane
//     CDetFrontRad = cdr + PlaneZ/2.0 - CDetModuleZ;         //origin to CDet Front
  
//     switch( fDetCon->GetCDetConfigOption() ) {
//     case 2:
//       CDetPhi = 0.0;
//       Layeroffset = 0.0;
//       break;
//     case 3:
//       CDetPhi = 2.0*atan( CDetModuleY/(2.0*CDetFrontRad) );
//       Layeroffset = 0.30*cm;
//       break;
//     default:
//       CDetPhi = 0.0;
//       Layeroffset = 0.0;
//       break;
//     }
//     G4RotationMatrix *CDetModuleRotTop = new G4RotationMatrix;
//     CDetModuleRotTop->rotateX(CDetPhi);

//     G4RotationMatrix *CDetModuleRotBottom = new G4RotationMatrix;
//     CDetModuleRotBottom->rotateX(-CDetPhi);

//     G4ThreeVector CDetPos( cdr*sin(-fBBang)+(offset+ (4.212/2.0)*cm )*cos(fBBang), 
//     			   -(4.212/2.0)*cm, 
//     			   cdr*cos(-fBBang)+(offset+(4.212/2.0)*cm)*sin(fBBang) );

//     new G4PVPlacement( bbrm, CDetPos, PlaneLog, "CDet", worldlog, false, 0);

//     double DisplaceZ = (CDetModuleY/2.0)*sin(CDetPhi) + (CDetModuleZ/2.0)*( 1-cos(CDetPhi) ); 
//     double DisplaceY = (CDetModuleY/2.0)*( 1-cos(CDetPhi) ) - (CDetModuleZ/2.0)*sin(CDetPhi);

//     G4ThreeVector postop(TopBottomModuleOffset, CDetModuleY-DisplaceY, (PlaneZ-CDetModuleZ)/2.0-DisplaceZ-z_Al);
//     G4ThreeVector posbot(TopBottomModuleOffset, -CDetModuleY+DisplaceY, (PlaneZ-CDetModuleZ)/2.0-DisplaceZ-z_Al);

//     //Layer 2 (closest to ECal)
//     new G4PVPlacement( 0,G4ThreeVector(0.0,0.0,(PlaneZ-z_Al)/2.0), Allog, "AlShield", PlaneLog, false, 0 );

//     new G4PVPlacement( CDetModuleRotTop, postop, CDetModuleLog, "CDetModulePhys", PlaneLog, false, 1, true );
//     new G4PVPlacement( 0, G4ThreeVector(-MiddleModuleOffset, 0.0, PlaneZ/2.0-CDetModuleZ/2.0-z_Al),
// 		       CDetModuleLog, "CDetModulePhys", PlaneLog, false, 2, true); 
//     new G4PVPlacement( CDetModuleRotBottom, posbot, CDetModuleLog, "CDetModulePhys", PlaneLog, false, 3, true );

//     //Layer 1(closest to Target)
//     CDetFrontRad = cdr + PlaneZ/2.0 - 2*CDetModuleZ; 

//     switch( fDetCon->GetCDetConfigOption() ) {
//     case 2:
//       CDetPhi = 0.0;
//       Layeroffset = 0.0;
//       break;
//     case 3:
//       CDetPhi = 2.0*atan( CDetModuleY/(2.0*CDetFrontRad) );
//       Layeroffset = 0.30*cm;
//       break;
//     default:
//       CDetPhi = 0.0;
//       Layeroffset = 0.0;
//       break;
//     }
  
//     G4RotationMatrix *CDetModuleRotTopLayer1 = new G4RotationMatrix;
//     CDetModuleRotTopLayer1->rotateX(CDetPhi);
//     G4RotationMatrix *CDetModuleRotBotLayer1 = new G4RotationMatrix;
//     CDetModuleRotBotLayer1->rotateX(-CDetPhi);

//     DisplaceZ = (CDetModuleY/2.0)*sin(CDetPhi) + (CDetModuleZ/2.0)*( 1-cos(CDetPhi) ); 
//     DisplaceY = (CDetModuleY/2.0)*( 1-cos(CDetPhi) ) - (CDetModuleZ/2.0)*sin(CDetPhi);

//     postop.set(TopBottomModuleOffset, CDetModuleY-DisplaceY, (PlaneZ-3*CDetModuleZ)/2.0-DisplaceZ-Layeroffset-z_Al);
//     posbot.set(TopBottomModuleOffset, -CDetModuleY+DisplaceY, (PlaneZ-3*CDetModuleZ)/2.0-DisplaceZ-Layeroffset-z_Al);

//     new G4PVPlacement( CDetModuleRotTopLayer1, postop, CDetModuleLog, "CDetModulePhys", PlaneLog, false, 4, true );
//     new G4PVPlacement( 0, G4ThreeVector(-MiddleModuleOffset, 0.0, (PlaneZ-3*CDetModuleZ)/2.0 - Layeroffset-z_Al),
// 		       CDetModuleLog, "CDetModulePhys", PlaneLog, false, 5, true); 
//     new G4PVPlacement( CDetModuleRotBotLayer1, posbot, CDetModuleLog, "CDetModulePhys", PlaneLog, false, 6, true );
//   }

//   // VISUALIZATIONS
//   CDetLog->SetVisAttributes( G4VisAttributes::Invisible );
//   CDetBarLog->SetVisAttributes( G4VisAttributes::Invisible );
//   CDetModuleLog->SetVisAttributes( G4VisAttributes::Invisible );
//   PlaneLog->SetVisAttributes( G4VisAttributes::Invisible );
//   CDetModuleLog4->SetVisAttributes( G4VisAttributes::Invisible );

//   //G4VisAttributes *CDetPlaneVis = new G4VisAttributes( G4Colour::Green() );
//   //CDetPlaneVis->SetForceWireframe(true);
//   //PlaneLog->SetVisAttributes(CDetPlaneVis);

//   G4VisAttributes *CDetMod4Vis = new G4VisAttributes( G4Colour::Green() );
//   CDetMod4Vis->SetForceWireframe(true);
//   //CDetModuleLog4->SetVisAttributes(CDetMod4Vis);
 
//   // Mylar
//   G4VisAttributes *CDetMylarVis = new G4VisAttributes( G4Colour(0.5,0.5,0.5) );
//   CDetMylarLog->SetVisAttributes( CDetMylarVis );
//   //CDetMylarLog->SetVisAttributes( G4VisAttributes::Invisible );

//   // Waveguide
//   G4VisAttributes *CDetWLSVis = new G4VisAttributes( G4Colour(0.0, 1.0, 1.0) );
//   WLSguideLog->SetVisAttributes( CDetWLSVis );

//   // Scintillating Material
//   G4VisAttributes *CDetScintVis = new G4VisAttributes( G4Colour( 0.5, 0.0, 0.8) );
//   CDScintBoreLog->SetVisAttributes( CDetScintVis );
  
//   // PMT
//   G4VisAttributes *CDetPMTVis = new G4VisAttributes( G4Colour(0.3,0.3,0.3) );
//   CDetPMTLog->SetVisAttributes( CDetPMTVis );
	
//   // Al
//   G4VisAttributes *AlVis = new G4VisAttributes(G4Colour(0.9, 0.9, 0.9));
//   AlVis->SetForceWireframe(true);
//   Allog->SetVisAttributes(AlVis);

//   /*************************************************************************************
//    ********************************        ECAL      ***********************************
//    *************************************************************************************/
 
//   //Dimensions are different than ecal_log because ECal was shrunk significantly when the crescent shape was developed
//   //x_earm = 47.5*4.212 (the .5 comes from staggering effects)
//   //y_earm = 75*4.212
//   //z_earm = 20 + 2*4 + 2.40

//   // Make CDet Option 1
//   G4Box *earm_mother_box = new G4Box("earm_mother_box", x_earm/2.0, y_earm/2.0, z_earm/2.0);
//   G4LogicalVolume *earm_mother_log = new G4LogicalVolume(earm_mother_box, GetMaterial("Air"), "earm_mother_log");

//   if( (fDetCon->StepLimiterList).find( "earm_mother_log" ) != (fDetCon->StepLimiterList).end() ){
//     earm_mother_log->SetUserLimits( new G4UserLimits( 0.0, 0.0, 0.0, DBL_MAX, DBL_MAX ) );

//     G4String sdname = "Earm/ECAL_box";
//     G4String collname = "ECAL_boxHitsCollection";
//     G4SBSCalSD *earm_mother_SD = NULL;
//     if( !( (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(sdname) ) ){
//       G4cout << "Adding ECAL_box sensitive detector to SDman..." << G4endl;
//       earm_mother_SD = new G4SBSCalSD( sdname, collname );
//       fDetCon->fSDman->AddNewDetector( earm_mother_SD );
//       (fDetCon->SDlist).insert( sdname );
//       fDetCon->SDtype[sdname] = G4SBS::kCAL;
//       (earm_mother_SD->detmap).depth = 0;
//     }
//     earm_mother_log->SetSensitiveDetector( earm_mother_SD );
//   }

//   // Polyethylene Box
//   G4Box *polybox = new G4Box("polybox", x_earm/2.0, y_earm/2.0, polydepth/2.0 );
//   G4LogicalVolume *polybox_log = new G4LogicalVolume( polybox, GetMaterial("Polyethylene"), "polybox_log");

//   // Coordinate Detector - 2 planes
//   G4Box *CD_box = new G4Box("CD_box", x_earm/2.0, y_earm/2.0, CD_depth/2.0 );
//   G4LogicalVolume *CD_log = new G4LogicalVolume(CD_box, GetMaterial("PLASTIC_SC_VINYLTOLUENE"), "CD_log");

//   // Aluminum Shielding in front of ECAL - make it big enough to enclose crescent
//   G4Box *Al_box = new G4Box( "Al_box", x_earm/2.0, y_earm/2.0, z_Al/2.0 );
//   G4LogicalVolume *Al_log = new G4LogicalVolume ( Al_box, GetMaterial("RICHAluminum"), "Al_log" );

//   // Place everything inside the Earm Mother Volume
//   new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, z_earm/2.0 - z_Al/2.0 ), Al_log, "Al", earm_mother_log, false, 0 );
//   new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, z_earm/2.0 - z_Al - CD_depth/2.0 ), CD_log, "plane1", earm_mother_log, false, 0 );
//   new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, z_earm/2.0 - z_Al - CD_depth - CD_depth/2.0 ), CD_log, "plane2", earm_mother_log, false, 1 );
//   new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, z_earm/2.0 - z_Al - 2*CD_depth - polydepth/2.0 ), polybox_log, "Polyethylene", earm_mother_log, false, 0 );

//   if( fDetCon->GetCDetConfigOption() == 1 ) {
//     cdr = fBBdist - z_ecal - z_earm/2.0;
//     //Note the offsets in x,y,z of 4.212/2.0
//     new G4PVPlacement( bbrm, G4ThreeVector( cdr*sin(-fBBang)+(offset + (4.212/2.0)*cm )*cos(fBBang), 
// 					    -(4.212/2.0)*cm, 
// 					    cdr*cos(-fBBang) + (offset+(4.212/2.0)*cm)*sin(fBBang) ), 
// 		       earm_mother_log, "CDet & Filter", worldlog, false, 0);
//   }
//   if( fDetCon->GetCDetConfigOption() == 2 ) {
//     cdr = fBBdist - z_ecal - PlaneZ/2.0;  
//     new G4PVPlacement( 0, G4ThreeVector(0.0,0.0, PlaneZ/2.0 - z_Al - 2*CDetModuleZ - polydepth/2.0), 
// 		       polybox_log, "Polyphys", PlaneLog, false, 0 );
//   }
  
//   // ECAL CONSTRUCTION - START

//   //Define a Mother Volume to place the ECAL modules
//   G4Box *ecal_box = new G4Box( "ecal_box", x_ecal/2.0, y_ecal/2.0, z_ecal/2.0 );
//   G4LogicalVolume *ecal_log = new G4LogicalVolume( ecal_box, GetMaterial("Air"), "ecal_log" );
//   new G4PVPlacement( bbrm, G4ThreeVector( bbr*sin(-fBBang)+offset*cos(fBBang), 0.0, bbr*cos(-fBBang)+offset*sin(fBBang) ), ecal_log, "ECal Mother", worldlog, false, 0 );

//   //Dimensions of mylar/air wrapping:
//   G4double mylar_wrapping_size = 0.00150*cm;
//   G4double air_wrapping_size = 0.00450*cm;
//   G4double mylar_plus_air = mylar_wrapping_size + air_wrapping_size;

//   //TYPE1 MODULE - 4.200 x 4.200 x 45.000 cm^3 + layer of air (0.0045cm) + layer of mylar (0.0015cm)
//   G4double x_module_type1 = 4.2000*cm + (2*mylar_plus_air); //4.212
//   G4double y_module_type1 = 4.2000*cm + (2*mylar_plus_air); //4.212
//   G4double z_module_type1 = 45.0000*cm + mylar_plus_air;    //45.006 - only one side in z has mylar+air
//   G4Box *module_type1 = new G4Box( "module_type1", x_module_type1/2.0, y_module_type1/2.0, z_module_type1/2.0 );

//   //Define a new mother volume which will house the contents of a module(i.e. TF1 & Mylar)
//   G4LogicalVolume *module_log_type1 = new G4LogicalVolume( module_type1, GetMaterial("Special_Air"), "module_log_type1" );

//   //MYLAR
//   //Subtraction solid to leave us with a 0.0015cm "mylar wrapping" with one end open
//   G4double x_sub = x_module_type1 - (2*mylar_wrapping_size);
//   G4double y_sub = y_module_type1 - (2*mylar_wrapping_size);
//   G4double z_sub = z_module_type1; //going to be shifted in G4SubtractionSolid in order to get 0.0015cm of mylar on one end
//   G4Box *mylar_module1_subtract = new G4Box( "mylar_module1_subtract", x_sub/2.0, y_sub/2.0, z_sub/2.0 );
//   G4SubtractionSolid *mylarwrap = new G4SubtractionSolid( "mylarwrap", module_type1, mylar_module1_subtract, 0, G4ThreeVector(0.0, 0.0, mylar_wrapping_size));
//   G4LogicalVolume *mylar_wrap_log = new G4LogicalVolume( mylarwrap, GetMaterial("Mylar"), "mylar_wrap_log" );

//   //TF1
//   G4double x_TF1 = 4.200*cm, y_TF1 = 4.200*cm, z_TF1 = 45.000*cm;
//   G4Box *TF1_box = new G4Box( "TF1_box", x_TF1/2.0, y_TF1/2.0, z_TF1/2.0 );
//   G4LogicalVolume *TF1_log = new G4LogicalVolume ( TF1_box, GetMaterial("TF1"), "TF1_log" );

//   G4String ECalTF1SDname = "Earm/ECalTF1";
//   G4String ECalTF1collname = "ECalTF1HitsCollection";
//   G4SBSCalSD *ECalTF1SD = NULL;
    
//   if( !( (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(ECalTF1SDname) ) ){
//     G4cout << "Adding ECal TF1 Sensitive Detector to SDman..." << G4endl;
//     ECalTF1SD = new G4SBSCalSD( ECalTF1SDname, ECalTF1collname );
//     fDetCon->fSDman->AddNewDetector( ECalTF1SD );
//     (fDetCon->SDlist).insert(ECalTF1SDname);
//     fDetCon->SDtype[ECalTF1SDname] = G4SBS::kCAL;
//     //fDetCon->SDarm[ECalTF1SDname] = G4SBS::kEarm;

//     (ECalTF1SD->detmap).depth = 1;
//   }
//   TF1_log->SetSensitiveDetector( ECalTF1SD );

//   if( (fDetCon->StepLimiterList).find( ECalTF1SDname ) != (fDetCon->StepLimiterList).end() ){
//     TF1_log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
//   }

//   //Place TF1 mother & Mylar inside module_log_type1 which is already full of Special_Air
//   new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, 0.0), mylar_wrap_log, "Mylar_Wrap", module_log_type1, false, 0 );
//   new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, mylar_plus_air/2.0), TF1_log, "TF1", module_log_type1, false, 0 );

//   //Want blocks to be optically independent => assign an optical skin to the Mylar
//   new G4LogicalSkinSurface( "Mylar_skin", mylar_wrap_log, GetOpticalSurface( "Mirrsurf" ));
    
//   //Build PMT detector which will be placed at the end of the module
//   //Starting with the Quartz window
//   G4double PMT_window_radius = 1.25*cm;
//   G4double PMT_window_depth = 0.20*cm;
//   G4Tubs *PMT_window = new G4Tubs( "PMT_window", 0.0*cm, PMT_window_radius, PMT_window_depth/2.0, 0.0, twopi );
//   G4LogicalVolume *PMT_window_log = new G4LogicalVolume( PMT_window, GetMaterial("QuartzWindow_ECal"), "PMT_window_log" );

//   //PMT
//   G4double PMT_radius = 1.25*cm;
//   G4double PMT_depth = 0.20*cm;
//   G4Tubs *PMTcathode_ecal = new G4Tubs( "PMTcathode_ecal", 0.0*cm, PMT_radius, PMT_depth/2.0, 0.0, twopi );
//   G4LogicalVolume *PMTcathode_ecal_log = new G4LogicalVolume( PMTcathode_ecal, GetMaterial("Photocathode_material_ecal"), "PMTcathode_ecal_log" );

//   //PMTcathode_ecal_log is the sensitive detector, assigned to ECalSD which detects optical photons
//   G4String ECalSDname = "Earm/ECal";
//   G4String ECalcollname = "ECalHitsCollection";
//   G4SBSECalSD *ECalSD = NULL;

//   if( !( (G4SBSECalSD*) sdman->FindSensitiveDetector(ECalSDname) ) ){
//     G4cout << "Adding ECal PMT Sensitive Detector to SDman..." << G4endl;
//     ECalSD = new G4SBSECalSD( ECalSDname, ECalcollname );
//     sdman->AddNewDetector( ECalSD );
//     (fDetCon->SDlist).insert(ECalSDname);
//     fDetCon->SDtype[ECalSDname] = G4SBS::kECAL;
//     //fDetCon->SDarm[ECalSDname] = G4SBS::kEarm;
//     (ECalSD->detmap).depth = 0;
//   }
//   PMTcathode_ecal_log->SetSensitiveDetector( ECalSD );

//   G4int x_number_ecal = 58;
//   G4int y_number_ecal = 88;   

//   int maxrow = -1;
//   int maxcol = -1;
//   int minrow = y_number_ecal;
//   int mincol = x_number_ecal;

//   G4String filename = "database/";
//   filename += fDetCon->GetECALmapfilename();
//   ifstream ecal_map_file;
  
//   if( fDetCon->GetECALmapfilename() != "" ){
//     ecal_map_file.open(filename.data());
//   }
 
//   set<pair<int,int> > active_cells;

//   G4String Line;

//   if( ecal_map_file.is_open() ){
//     while ( Line.readLine( ecal_map_file ) ){
//       //G4cout << "Read Line from ecal map file: \"" <<  Line.data() << "\"" << G4endl;
//       if( !(Line[0] == '#') ){
// 	int row, col;
// 	sscanf( Line.data(), "%d %d", &row, &col );
// 	if( active_cells.size() == 0 || row > maxrow ) maxrow = row;
// 	if( active_cells.size() == 0 || row < minrow ) minrow = row;
// 	if( active_cells.size() == 0 || col > maxcol ) maxcol = col;
// 	if( active_cells.size() == 0 || col < mincol ) mincol = col;
// 	active_cells.insert(make_pair(row,col));
//       } 
//     }
//   } else {
//     for( G4int i=0; i<x_number_ecal; i++){
//       for( G4int j=0; j<y_number_ecal; j++){
// 	active_cells.insert(make_pair(j,i));
// 	int row=j, col=i; 
// 	if( active_cells.size() == 0 || row > maxrow ) maxrow = row;
// 	if( active_cells.size() == 0 || row < minrow ) minrow = row;
// 	if( active_cells.size() == 0 || col > maxcol ) maxcol = col;
// 	if( active_cells.size() == 0 || col < mincol ) mincol = col;
//       }
//     }
//   }

//   G4double x_position = x_ecal/2.0-x_module_type1/2.0 , y_position = y_ecal/2.0-y_module_type1/2.0;
  
//   G4int copy_number_PMT = 0;  //label modules

//   //Need a Steel module to fill voids
//   G4int x_steel = x_module_type1, y_steel = y_module_type1, z_steel = z_module_type1;
//   G4Box *steel_box = new G4Box("steel_box", x_steel/2.0, y_steel/2.0, z_steel/2.0);
//   G4LogicalVolume *steel_log = new G4LogicalVolume( steel_box, GetMaterial("Steel"), "steel_log");

//   //And a half Steel module to fill in voids caused by staggering
//   G4Box *steel_box_half = new G4Box("steel_box_half", 0.5*(x_steel/2.0), y_steel/2.0, z_steel/2.0);
//   G4LogicalVolume *steel_log_half = new G4LogicalVolume( steel_box_half, GetMaterial("Steel"), "steel_log");
//   //int module_number; //counts modules if needed

//   //Iterating to build ECal 
 	
//   for( G4int i=0; i<y_number_ecal; i++ ){
//     G4double ytemp = y_position - i*(y_module_type1);

//     for(G4int j=0; j<x_number_ecal; j++){	
//       G4double xtemp_even = x_position - j*(x_module_type1);
//       G4double xtemp_odd = x_position - x_module_type1/2.0 - j*(x_module_type1);
	  
//       pair<int,int> rowcol( i, j );
	  
//       if( active_cells.find( rowcol ) != active_cells.end() ){ //This row and column are active:
// 	(ECalSD->detmap).Row[copy_number_PMT] = i;
// 	(ECalSD->detmap).Col[copy_number_PMT] = j;
// 	(ECalTF1SD->detmap).Row[copy_number_PMT] = i;
// 	(ECalTF1SD->detmap).Col[copy_number_PMT] = j;
// 	//EVEN ROWS
// 	if(i%2==0){ 
// 	  //Type1 Module
// 	  new G4PVPlacement(0, G4ThreeVector( xtemp_even, ytemp, -PMT_depth), module_log_type1, "Type1Module", ecal_log, false, copy_number_PMT );
// 	  //PMT_window
// 	  new G4PVPlacement( 0, G4ThreeVector(xtemp_even, ytemp, z_ecal/2.0-PMT_depth-PMT_window_depth/2.0), PMT_window_log,"PMT_window_pv", ecal_log, false, copy_number_PMT );
// 	  //PMT_cathode
// 	  new G4PVPlacement( 0, G4ThreeVector(xtemp_even, ytemp, z_ecal/2.0-PMT_depth/2.0), PMTcathode_ecal_log, "PMTcathode_pv", ecal_log, false, copy_number_PMT );
	      
// 	  (ECalSD->detmap).LocalCoord[copy_number_PMT] = G4ThreeVector(xtemp_even,ytemp,z_ecal/2.0-PMT_depth/2.0);
// 	  //(ECalSD->detmap).GlobalCoord[copy_number_PMT] = G4ThreeVector(0.0,0.0,0.0);
// 	  (ECalTF1SD->detmap).LocalCoord[copy_number_PMT] = G4ThreeVector(xtemp_even,ytemp,-PMT_depth);
// 	  //(ECalTF1SD->detmap).GlobalCoord[copy_number_PMT] = G4ThreeVector(0.0,0.0,0.0);
// 	}

// 	//ODD ROWS
// 	else{
// 	  //Type1 Module
// 	  new G4PVPlacement( 0, G4ThreeVector( xtemp_odd, ytemp, -PMT_depth), module_log_type1, "Type1Module", ecal_log, false, copy_number_PMT );
// 	  //PMT_window
// 	  new G4PVPlacement( 0, G4ThreeVector(xtemp_odd, ytemp, z_ecal/2.0-PMT_depth-PMT_window_depth/2.0), PMT_window_log,"PMT_window_pv", ecal_log, false, copy_number_PMT );
// 	  //PMT_cathode
// 	  new G4PVPlacement( 0, G4ThreeVector(xtemp_odd, ytemp, z_ecal/2.0-PMT_depth/2.0), PMTcathode_ecal_log, "PMTcathode_pv", ecal_log, false, copy_number_PMT );
	      
// 	  (ECalSD->detmap).LocalCoord[copy_number_PMT] = G4ThreeVector(xtemp_odd,ytemp,z_ecal/2.0-PMT_depth/2.0);
// 	  //(ECalSD->detmap).GlobalCoord[copy_number_PMT] = G4ThreeVector(0.0,0.0,0.0);
// 	  (ECalTF1SD->detmap).LocalCoord[copy_number_PMT] = G4ThreeVector(xtemp_odd,ytemp,-PMT_depth);
// 	  //(ECalTF1SD->detmap).GlobalCoord[copy_number_PMT] = G4ThreeVector(0.0,0.0,0.0);
// 	}
// 	//	module_number++;
// 	copy_number_PMT++;
//       } else { //Steel Filler:
// 	if( i>=minrow-2 && i<=maxrow+2 && j>=mincol-2 && j<=maxcol+2 ){
// 	  if(i%2==0){
// 	    new G4PVPlacement(0, G4ThreeVector(xtemp_even, ytemp, -PMT_depth ), steel_log, "steel_pv", ecal_log, false, 0);
	    
// 	  } else {
// 	    new G4PVPlacement(0, G4ThreeVector(xtemp_odd, ytemp, -PMT_depth ), steel_log, "steel_pv", ecal_log, false, 0);   
// 	  }
// 	}
//       }
//     }

//     if( i >= minrow-2 && i <= maxrow+2 ){

//       int colleft = max(0,mincol-2);
//       int colright = min(maxcol+2,x_number_ecal-1);

//       double x_left = x_position - (colleft)*x_module_type1 + x_steel/4.0;
//       double x_right = x_position - (colright+0.5)*x_module_type1 - x_steel/4.0;

//       if( i%2 == 0 ){
// 	// new G4PVPlacement(0, G4ThreeVector( (x_earm/2.0)-(x_steel/4.0), ytemp, -PMT_depth ), steel_log_half, "steel_pv", ecal_log, false, 0);
// 	new G4PVPlacement(0, G4ThreeVector( x_right, ytemp, -PMT_depth ), steel_log_half, "steel_pv", ecal_log, false, 0);
//       } else {
// 	// new G4PVPlacement(0, G4ThreeVector(-(x_earm/2.0)+(x_steel/4.0), ytemp, -PMT_depth ), steel_log_half, "steel_pv", ecal_log, false, 0);
// 	new G4PVPlacement(0, G4ThreeVector( x_left, ytemp, -PMT_depth ), steel_log_half, "steel_pv", ecal_log, false, 0);
//       }
//     }
//   }
//   //END MODULE PLACEMENT

//   //VISUALIZATIONS
//   //Make mother volume invisible
//   ecal_log->SetVisAttributes( G4VisAttributes::Invisible );
	
//   //Mylar
//   G4VisAttributes *mylar_colour = new G4VisAttributes(G4Colour( 0.5, 0.5, 0.5 ) );
//   mylar_wrap_log->SetVisAttributes(G4VisAttributes::Invisible);
  
//   //Air
//   G4VisAttributes *air_colour = new G4VisAttributes(G4Colour( G4Colour::Blue() ));
//   module_log_type1->SetVisAttributes(G4VisAttributes::Invisible);

//   //TF1
//   G4VisAttributes *TF1_colour = new G4VisAttributes(G4Colour( 0.8, 0.8, 0.0 ) );
//   TF1_log->SetVisAttributes(TF1_colour);

//   //Steel
//   G4VisAttributes *steel_colour = new G4VisAttributes(G4Colour( 0.4, 0.4, 0.4 ) );
//   steel_log->SetVisAttributes(steel_colour);
//   steel_log_half->SetVisAttributes(steel_colour);

//   //PMTcathode
//   G4VisAttributes *PMT_colour = new G4VisAttributes(G4Colour( G4Colour::Blue() ));
//   PMT_colour->SetForceLineSegmentsPerCircle( 12 );
//   PMTcathode_ecal_log->SetVisAttributes(PMT_colour);

//   //Electron Arm Package - Houses Polyethylene, CDet, Al
//   earm_mother_log->SetVisAttributes(G4VisAttributes::Invisible);

//   //Polyethylene
//   G4VisAttributes *poly_colour = new G4VisAttributes(G4Colour( 0.2,0.3,0.4 ));
//   poly_colour->SetForceWireframe(true);
//   polybox_log->SetVisAttributes(poly_colour);

//   //CDet
//   G4VisAttributes *CD_colour = new G4VisAttributes(G4Colour(0.5, 0.4, 0.1));
//   CD_colour->SetForceWireframe(true);
//   CD_log->SetVisAttributes(CD_colour);

//   //Al
//   G4VisAttributes *Al_colour = new G4VisAttributes(G4Colour(0.9, 0.9, 0.9));
//   Al_colour->SetForceWireframe(true);
//   Al_log->SetVisAttributes(Al_colour);
		    
// } 

void G4SBSEArmBuilder::MakeGMnGEMShielding_update( G4LogicalVolume *motherLog ){
  // GEM electronics hut 
  // date: 7/16/21 
  // Drawing: A00000-06-00-0200-P001 
  // coordinate system: +x = beam left, +y = up, +z = downstream 

  G4double inch  = 2.54*cm;
  G4double foot  = 12.*inch; 

  // outer dimensions of the shielding blocks
  G4double GboxX = 104.0*inch;
  G4double GboxY = 52.0*inch;
  G4double GboxZ = 130*inch;
  G4Box *shield_tmp = new G4Box("shield_tmp",GboxX/2.,GboxY/2.,GboxZ/2.); 

  // define a cut to cut out the interior 
  G4double GboxX_cut = 52.0*inch;
  G4double GboxY_cut = 70.0*inch;
  G4double GboxZ_cut = 130.0*inch;
  G4Box *shield_tmp_cut = new G4Box("shield_tmp_cut",GboxX_cut/2.,GboxY_cut/2.,GboxZ_cut/2.);

  G4ThreeVector P_cut = G4ThreeVector(0,-9*inch,-26.*inch);
  G4SubtractionSolid *bunker = new G4SubtractionSolid("gemBunker",shield_tmp,shield_tmp_cut,0,P_cut);

  G4LogicalVolume *bunker_log = new G4LogicalVolume(bunker,GetMaterial("Steel"),"gemBunker_log"); 

  bool checkOverlaps = true;

  // placement of the electronics hut (approximate)
  G4double x0 = 3.*m;      // FIXME: Should be 3.2 m
  G4double y0 = -10.*foot; // about 10 feet below the beamline is the floor
  G4double z0 = 7.0*m;     // FIXME: Should be 7.7 m  
  G4ThreeVector Pb = G4ThreeVector(x0,y0+GboxY/2.,z0);  
  G4RotationMatrix *rmb = new G4RotationMatrix();
  rmb->rotateY(90.*deg);

  new G4PVPlacement(rmb,
                    Pb, 
		    bunker_log,         // logical volume
                    "gemBunker_phys",   // name 
                    motherLog,          // logical mother 
                    true,               // boolean?  
                    0,                  // copy number 
                    checkOverlaps);     // check overlaps  

  // shielding plates (goes on top of bunker)
  // small 
  G4double x_sm = 96.0*inch; 
  G4double y_sm = 2.5*inch; 
  G4double z_sm = 21.5*inch;
  G4Box *smallPlate = new G4Box("smallPlate",x_sm/2.,y_sm/2.,z_sm/2.);  
  // medium 
  G4double x_med = 96.0*inch; 
  G4double y_med = 2.5*inch; 
  G4double z_med = 23.5*inch;
  G4Box *mediumPlate = new G4Box("mediumPlate",x_med/2.,y_med/2.,z_med/2.);  
  // large 
  G4double x_lg = 120.0*inch; 
  G4double y_lg = 2.5*inch; 
  G4double z_lg = 21.5*inch;
  G4Box *largePlate = new G4Box("largePlate",x_lg/2.,y_lg/2.,z_lg/2.);  

  // union the plates
  G4double dx = 0.5*(GboxX - z_lg); 
  G4double dy = 0.1*inch; 
  G4double dz = 13.0*inch; 
  // large plates 
  // pair 1
  G4double zp1 = 0;
  G4ThreeVector P1b  = G4ThreeVector(0,-y_lg,zp1);
  G4UnionSolid *bunkerTop = new G4UnionSolid("bt_0",largePlate,largePlate,0,P1b);
  // G4ThreeVector P1t  = G4ThreeVector(0,y_lg,zp1);
  // bunkerTop = new G4UnionSolid("bt_1",bunkerTop,largePlate,0,P1t);
  // pair 2 
  G4double zp2 = z_lg;
  G4ThreeVector P2b  = G4ThreeVector(0,-y_lg,zp2);
  bunkerTop = new G4UnionSolid("bt_2",bunkerTop,largePlate,0,P2b);
  G4ThreeVector P2t  = G4ThreeVector(0,0,zp2);
  bunkerTop = new G4UnionSolid("bt_3",bunkerTop,largePlate,0,P2t);
  // pair 3: med on top, sm on bottom 
  G4double zp3 = zp2 + 0.5*(z_sm + z_lg);
  G4ThreeVector P3b = G4ThreeVector(0,-y_sm,zp3);
  bunkerTop = new G4UnionSolid("bt_4",bunkerTop,smallPlate,0,P3b);
  G4ThreeVector P3t = G4ThreeVector(0,0,zp3);
  bunkerTop = new G4UnionSolid("bt_5",bunkerTop,mediumPlate,0,P3t);
  // pair 4: med on top, sm on bottom 
  G4double zp4 = zp3 + 0.5*(z_sm + z_med) + dz;
  G4ThreeVector P4b = G4ThreeVector(0,-y_sm,zp4);
  bunkerTop = new G4UnionSolid("bt_6",bunkerTop,smallPlate,0,P4b);
  G4ThreeVector P4t = G4ThreeVector(0,0,zp4);
  bunkerTop = new G4UnionSolid("gemBunkerTop",bunkerTop,mediumPlate,0,P4t);

  G4VisAttributes *visTop = new G4VisAttributes();
  visTop->SetColour( G4Colour::Magenta() );

  G4LogicalVolume *bunkerTop_log = new G4LogicalVolume(bunkerTop,GetMaterial("Steel"),"bunkerTop_log");
  bunkerTop_log->SetVisAttributes(visTop);  

  // placement
  G4double xt = x0 + GboxZ/2. - z_lg/2.; 
  G4double yt = y0 + GboxY + 1.5*y_lg; 
  G4double zt = z0 + (1./8.)*(x_lg/2.-GboxX/2.); 
  G4ThreeVector Pbt = G4ThreeVector(xt,yt,zt);

  new G4PVPlacement(rmb,                   // rotation
                    Pbt,                   // position  
		    bunkerTop_log,         // logical volume
                    "gemBunkerTop_phys",   // name 
                    motherLog,             // logical mother 
                    true,                  // boolean?  
                    0,                     // copy number 
                    checkOverlaps);        // check overlaps  

   // In order to calculate the dose, we need a SD of type CAL:
   // total volume = 101.6 x 101.6 x 2.54 cm^3 = 26219.302 cm^3 = 2638.961 in^3 
   // modified dimensions to match the same volume but fit inside the bunker 
   G4double ElecX = 38*inch; 
   G4double ElecY = 50*inch; 
   G4double ElecZ = 1.38*inch; 
  
   G4Box *Electronics = new G4Box( "Electronics" , ElecX/2.0, ElecY/2.0, ElecZ/2.0);

   G4VisAttributes *visElec = new G4VisAttributes();
   visElec->SetColour( G4Colour::Red() ); 

   G4LogicalVolume *Electronics_log = new G4LogicalVolume( Electronics , GetMaterial("Silicon"), "Electronics_log" );
   Electronics_log->SetVisAttributes(visElec);  

   G4String GEMElectronicsname = "Earm/GEMElectronics";
   G4String GEMElectronicscollname = "GEMElectronicsHitsCollection";
   G4SBSCalSD *GEMElecSD = NULL;
 
   GEMElectronicsname += "GMn";
   GEMElectronicscollname += "GMn"; 

   G4SDManager *sdman = fDetCon->fSDman;
 
   if( !( (G4SBSCalSD*) sdman->FindSensitiveDetector(GEMElectronicsname) )){
     G4cout << "Adding GEM electronics Sensitive Detector to SDman..." << G4endl;
     GEMElecSD = new G4SBSCalSD( GEMElectronicsname, GEMElectronicscollname );
     sdman->AddNewDetector(GEMElecSD);
     (fDetCon->SDlist).insert(GEMElectronicsname);
     fDetCon->SDtype[GEMElectronicsname] = G4SBS::kCAL;
     (GEMElecSD->detmap).depth = 1;
 
     fDetCon->SetTimeWindowAndThreshold( GEMElectronicsname );
   }
   Electronics_log->SetSensitiveDetector( GEMElecSD );
   
   if( (fDetCon->StepLimiterList).find( GEMElectronicsname ) != (fDetCon->StepLimiterList).end() ){
     Electronics_log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
   }
 
   // Place the electronics in our hut:
   G4double xg = x0; 
   G4double yg = y0 + ElecY/2.; 
   G4double zg = z0; 
   G4ThreeVector Pgem = G4ThreeVector(xg,yg,zg);

   new G4PVPlacement(rmb,                 // rotation relative to mother
                     Pgem,                // position relative to mother
 		     Electronics_log,     // logical volume
                     "GMn_Electronics",   // physical name
                     motherLog,           // logical mother
                     false,               // boolean?
                     0,                   // copy number
                     checkOverlaps);      // check overlaps 

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
  // G4double ElecX = 150.0*cm;
  // G4double ElecY = 40.0*cm;
  // G4double ElecZ = 0.5*cm;

  // D. Flay, 7/14/21 
  // Updated geometry of electronics 
  G4double ElecX = 101.6*cm; 
  G4double ElecY = 101.6*cm; 
  G4double ElecZ = 2.54*cm; 
 
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
    fDetCon->SDtype[GEMElectronicsname] = G4SBS::kCAL;
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
void G4SBSEArmBuilder::MakeBBSieveSlit(G4LogicalVolume *motherlog, G4ThreeVector pos)
{
  G4double bbsievew = 14.75*2.54*cm;
  G4double bbsieveh = 27.50*2.54*cm;
  G4double bbsieved = 1.50*2.54*cm;
  
  //BigBite Sieve Slit is a box with holes:
  G4Box *BBsieveplate_box = new G4Box("BBsieveplate_box", bbsievew/2.0, bbsieveh/2.0, bbsieved/2.0 );

  //Next, need to make cuts:
  G4double bbsieve_holediameter = 0.750*2.54*cm;
  G4double bbsieve_holeradius = bbsieve_holediameter/2.0;

  G4double bbsieve_holeangle1 = 18.0*deg;
  
  G4Tubs *bbsieveholecut = new G4Tubs("bbsieveholecut", 0.0, bbsieve_holeradius, bbsieved/cos(18.0*deg),0.0,360.0*deg);

  const G4int bbsieve_nrows = 13;
  const G4int bbsieve_ncols = 7;

  //Now we are going to attempt to cut holes in the sieve plate. Then we will deal with the slots:
  G4double bbsieve_holespaceY = 1.50*2.54*cm;
  G4double bbsieve_holespaceX = 49.2*mm; 

  //There are four areas where we need horizontal slots; these are in row 4, columns 2-3 and 4-5
  //row 10, columns 3-5, and row 12, column 5-6. These can be accomplished by the union of a box with two cylinders
 
  
  //We will need five rotation matrices to make the cuts (actually four since one is the identity):
  //G4RotationMatrix *bbsrot1 = new G4RotationMatrix;
  //bbsrot1->rotateX( -18.0*deg );
  //G4RotationMatrix *bbsrot2 = new G4RotationMatrix;
  //bbsrot2->rotateX( -9.0*deg );

  //G4RotationMatrix *bbsrot3 = new G4RotationMatrix;
  //bbsrot3->rotateX( +9.0*deg );
  //G4RotationMatrix *bbsrot4 = new G4RotationMatrix;
  //bbsrot4->rotateX( +18.0*deg );

  G4Box *Slot1Box = new G4Box("Slot1Box", bbsieve_holespaceY/2.0, bbsieve_holeradius, bbsieved/cos(18.0*deg) );
  G4Box *Slot2Box = new G4Box("Slot2Box", bbsieve_holespaceY, bbsieve_holeradius, bbsieved/cos(18.0*deg) );
			      
  
  G4double alpha1 = 9.0*deg;
  G4double alpha2 = 18.0*deg;
  
  G4ThreeVector bbsieveTopFrontCornerPos( 0, bbsieveh/2.0, -bbsieved/2.0 );
  G4ThreeVector bbsieveBottomFrontCornerPos( 0, -bbsieveh/2.0, -bbsieved/2.0 );

  G4SubtractionSolid *BBSievePlateCut;

  bool first=true;
  
  for( int irow=-6; irow<=6; irow++ ){

    for( int icol=-3; icol<=3; icol++ ){
      G4ThreeVector xhat_temp, yhat_temp, zhat_temp, r0_hole; 
      
      G4double stemp;

      G4RotationMatrix *rot_temp = new G4RotationMatrix;
      
      if( irow < -4 ){ //bottom two rows:
	xhat_temp = G4ThreeVector(0.0, cos(alpha2), sin(alpha2) );
	yhat_temp = G4ThreeVector(1.0, 0.0, 0.0 );
	zhat_temp = xhat_temp.cross(yhat_temp).unit();
	
	r0_hole = bbsieveBottomFrontCornerPos + xhat_temp*(1.942*2.54*cm + (irow+6)*1.938*2.54*cm);
	//Now compute the intersection of this with the z=0 plane:
	// (r0hole + s*nhat) dot uhat = 0;
	//s nhat dot uhat = -r0hole.uhat;
	stemp = -r0_hole.z() / zhat_temp.z();

	rot_temp->rotateX( -alpha2 );
	
      } else if( irow < -1 ){ //next three rows: -4, -3, -2 
	xhat_temp = G4ThreeVector(0.0, cos(alpha1), sin(alpha1) );
	yhat_temp = G4ThreeVector(1.0, 0.0, 0.0 );
	zhat_temp = xhat_temp.cross(yhat_temp).unit();

	r0_hole = bbsieveBottomFrontCornerPos + xhat_temp*(149.3*mm + (irow+4)*1.938*2.54*cm);

	stemp = -r0_hole.z() / zhat_temp.z();

	rot_temp->rotateX( -alpha1 );
      } else if( irow < 2 ){ //middle three rows: -1, 0, 1
	xhat_temp = G4ThreeVector(0,1,0);
	yhat_temp = G4ThreeVector(1,0,0);
	zhat_temp = xhat_temp.cross(yhat_temp).unit();

	r0_hole = xhat_temp*irow*1.938*2.54*cm;

	stemp = 0.0;
      } else if( irow < 5 ){ //next three rows:
	xhat_temp = G4ThreeVector(0.0, -cos(alpha1), sin(alpha1) );
	yhat_temp = G4ThreeVector(1,0,0);
	zhat_temp = xhat_temp.cross(yhat_temp).unit();

	r0_hole = bbsieveTopFrontCornerPos + xhat_temp*(149.3*mm + (4-irow)*1.938*2.54*cm);

	stemp = -r0_hole.z()/zhat_temp.z();

	rot_temp->rotateX( alpha1 );
      } else { //top two rows:
	xhat_temp = G4ThreeVector(0.0, -cos(alpha2), sin(alpha2) );
	yhat_temp = G4ThreeVector(1,0,0);
	zhat_temp = xhat_temp.cross(yhat_temp).unit();

	r0_hole = bbsieveTopFrontCornerPos + xhat_temp*(1.942*2.54*cm + (6-irow)*1.938*2.54*cm);

	rot_temp->rotateX( alpha2 );
      }
      
      G4ThreeVector holepos = r0_hole + stemp * zhat_temp + icol*bbsieve_holespaceY*yhat_temp;
    
      G4SubtractionSolid *NextCut;
      if( first ){
	first = false;
	NextCut = new G4SubtractionSolid( "NextCut", BBsieveplate_box, bbsieveholecut, rot_temp, holepos );
      } else {
	NextCut = new G4SubtractionSolid( "NextCut", BBSievePlateCut, bbsieveholecut, rot_temp, holepos );
      }

      //Copy "NextCut" to BBSievePlateCut and then delete NextCut
      //BBSievePlateCut = new G4SubtractionSolid( *NextCut );

      BBSievePlateCut = NextCut;
      BBSievePlateCut->SetName("BBSievePlateCut");

      if( irow == -3 && (icol == -2 || icol == 1 ) ){ //Cut two rectangular slots between columns -2--1 and columns 1-2
	holepos += 0.5*bbsieve_holespaceY*yhat_temp;
	NextCut = new G4SubtractionSolid( "NextCut", BBSievePlateCut, Slot1Box, rot_temp, holepos );

	BBSievePlateCut = NextCut;
	BBSievePlateCut->SetName("BBSievePlateCut");
	// delete NextCut;      
      }

      if( irow == 3 && icol == 0 ){
	NextCut = new G4SubtractionSolid( "NextCut", BBSievePlateCut, Slot2Box, rot_temp, holepos );

	BBSievePlateCut = NextCut;
	BBSievePlateCut->SetName("BBSievePlateCut");
	
      }

      if( irow == 5 && icol == -2 ){
	holepos += 0.5*bbsieve_holespaceY*yhat_temp;
	NextCut = new G4SubtractionSolid( "NextCut", BBSievePlateCut, Slot1Box, rot_temp, holepos );

	BBSievePlateCut = NextCut;
	BBSievePlateCut->SetName("BBSievePlateCut");	
      }
    }
    
    
    
  }

  G4LogicalVolume *BBSievePlate_log = new G4LogicalVolume( BBSievePlateCut, GetMaterial("Lead"), "BBSievePlate_log" );

  //Next we have to figure out where this sits in relation to the BB magnet: Naively we want it just in front of the coils:
  new G4PVPlacement( 0, pos, BBSievePlate_log, "BBSievePlate_phys", motherlog, 0, false, 0 );
  
  printf("Building BB sieve slit...\n");
}
 

//Sieve slit designed by Holly S. for optics studies 10.5.20
void G4SBSEArmBuilder::MakeNewBBSieveSlit(G4LogicalVolume *motherlog, G4ThreeVector pos)
{
  printf("Building New BB sieve slit...\n");
  
  //Plate
  //Plate dims - Assuming dimensions of box from extremes of holes and previous depth.
  G4double inch = 2.54*cm;

  G4double bbsievew = 14.75*inch;
  G4double bbsieveh = 27.50*inch;
  G4double bbsieved = 1.50*inch;
  
  //BigBite Sieve Slit is a box with holes:
  G4Box *Fullplate = new G4Box("Fullplate", bbsievew/2.0, bbsieveh/2.0, bbsieved/2.0 );

  //Hole Dimensions
  G4double HoleCenter_r = 0.125*inch;
  G4double Hole_r = 0.25*inch;
  G4double Hole_z = 3.0*inch; //Large enough to leave no solid volume at extreme displacements from center
  
  //Hole solids
  G4Tubs *Hole = new G4Tubs("Hole", 0, Hole_r, Hole_z, 0.0*deg, 360.0*deg);
  G4Tubs *HoleCenter = new G4Tubs("HoleCenter", 0, HoleCenter_r, Hole_z, 0.0*deg, 360.0*deg);

  //angular spacing of holes at the nominal distance of 51.825 inches
  G4double angspace_y = 1.65776*deg;
  G4double angspace_x = 1.06114*deg; 

  //G4SubtractionSolid *NextCut;
  G4SubtractionSolid *SievePlateCut;
  
  G4bool first = true;

  //Cut holes in the plate:

  cout << "cutting holes" << endl;
  
  
  for(G4int iy=-8; iy<=8; iy++ ){
    for(G4int ix=-5; ix<=5; ix++ ){
      G4double xangle = ix*angspace_x;
      G4double yangle = iy*angspace_y;

      //G4SubtractionSolid 
      
      //if( ix != 0 || iy != 0 ){ //compute rotation matrix for all holes but center hole:
      G4ThreeVector holeaxis( tan(xangle), tan(yangle), 1.0 );
      holeaxis = holeaxis.unit();
      
      G4ThreeVector zaxis(0.,0.,1.0);
      
      G4ThreeVector rotationaxis = (zaxis.cross(holeaxis)).unit();
      G4double rotationangle = acos( zaxis.dot(holeaxis) );

      G4RotationMatrix *holerot = new G4RotationMatrix;
      if( !(ix == 0 && iy == 0 ) ) {
	holerot->rotate( -rotationangle, rotationaxis );
      }
      
      //G4double platecenter_z = 68.8976*inch + bbsieved/2.0; //Assuming old placement
      G4double platecenter_z = 54.0157*inch + bbsieved/2.0; //Accounting for BB shift to 1.75m

      
      G4ThreeVector origin(0,0,-platecenter_z );
      G4double holedist = platecenter_z * sqrt(1.0 + pow(tan(xangle),2)+pow(tan(yangle),2) ); 
      
      G4ThreeVector holecenterpos = origin + holedist * holeaxis;
      
      cout << "(row,col,holeposx,holeposy,holeposz,rotationangle)=("
	   << iy << ", " << ix << ", " << holecenterpos.x()/inch << ", "
	   << holecenterpos.y()/inch << ", "
	   << holecenterpos.z()/inch << ", "
	   << rotationangle/deg << ")" << endl;

      cout << "rotation matrix = " << endl;
      holerot->print(cout);
      
      G4SubtractionSolid *NextCut;
      if( first ){
	first = false;
	if( !(ix==0&&iy==0) ){
	  NextCut = new G4SubtractionSolid( "NextCut", Fullplate, Hole, holerot, holecenterpos );
	} else { //it should be impossible to end up here:
	  NextCut = new G4SubtractionSolid( "NextCut", Fullplate, HoleCenter, holerot, holecenterpos );
	}
      } else {
	if( (iy==-3&&ix==-3)||(iy==0&&ix==0)||(iy==2&&ix==2) ){ //Cutting smaller holes
	  NextCut = new G4SubtractionSolid( "NextCut", SievePlateCut, HoleCenter, holerot, holecenterpos );
	}else{
	  if( iy==-3&&ix==1){
	    cout << "Skipping hole at y = -3 and x = 1" << endl; //Leaving one hole out
	  }else{
	    NextCut = new G4SubtractionSolid( "NextCut", SievePlateCut, Hole, holerot, holecenterpos );
	  }
	}
      }

      SievePlateCut = NextCut;
      SievePlateCut->SetName("NewBBSievePlateCut");
    }
  }

  G4cout << "finished holes..." << endl;
  
  G4LogicalVolume *NewBBSievePlate_log = new G4LogicalVolume(SievePlateCut, GetMaterial("Lead"), "New_BB_SievePlate_log");

  //For BB, plate is oriented along z axis and offset is passed into function.
  //G4double SBSsieve_dist = offset_z - ThickPlate_z/2.0;
  //Placement:
  //G4ThreeVector sievepos = SBSsieve_dist * G4ThreeVector( -sin(f48D48ang), 0, cos(f48D48ang) );
  //G4RotationMatrix *sieverot = new G4RotationMatrix;
  //sieverot->rotateY(f48D48ang);

  new G4PVPlacement(0, pos, NewBBSievePlate_log, "NewBBSievePlate_phys", motherlog, false, 0);
  
  G4cout << "Sieve plate finished" << G4endl;
  
}
 //Sieve slit designed by Holly S. - remove angled holes. For optics studies 10.29.20
void G4SBSEArmBuilder::MakeThirdBBSieveSlit(G4LogicalVolume *motherlog, G4ThreeVector pos)
{
  printf("Building new BB sieve slit without angled holes...\n");
  
  //Plate
  //Plate dims - Assuming dimensions of box from extremes of holes and previous depth.
  G4double inch = 2.54*cm;

  G4double bbsievew = 14.75*inch;
  G4double bbsieveh = 27.50*inch;
  G4double bbsieved = 1.50*inch;
  
  //BigBite Sieve Slit is a box with holes:
  G4Box *Fullplate = new G4Box("Fullplate", bbsievew/2.0, bbsieveh/2.0, bbsieved/2.0 );

  //Hole Dimensions
  G4double HoleCenter_r = 0.125*inch;
  G4double Hole_r = 0.25*inch;
  G4double Hole_z = 3.0*inch; //Large enough to leave no solid volume at extreme displacements from center
  
  //Hole solids
  G4Tubs *Hole = new G4Tubs("Hole", 0, Hole_r, Hole_z, 0.0*deg, 360.0*deg);
  G4Tubs *HoleCenter = new G4Tubs("HoleCenter", 0, HoleCenter_r, Hole_z, 0.0*deg, 360.0*deg);

  //angular spacing of holes at the nominal distance of 51.825 inches
  G4double angspace_y = 1.65776*deg;
  G4double angspace_x = 1.06114*deg; 

  //G4SubtractionSolid *NextCut;
  G4SubtractionSolid *SievePlateCut;
  
  G4bool first = true;

  //Cut holes in the plate:

  cout << "cutting holes" << endl;
  
  
  for(G4int iy=-8; iy<=8; iy++ ){
    for(G4int ix=-5; ix<=5; ix++ ){
      G4double xangle = ix*angspace_x;
      G4double yangle = iy*angspace_y;

      //G4SubtractionSolid 
      
      //if( ix != 0 || iy != 0 ){ //compute rotation matrix for all holes but center hole:
      G4ThreeVector holeaxis( tan(xangle), tan(yangle), 1.0 );
      holeaxis = holeaxis.unit();
      
      G4ThreeVector zaxis(0.,0.,1.0);
      
      G4ThreeVector rotationaxis = (zaxis.cross(holeaxis)).unit();
      G4double rotationangle = acos( zaxis.dot(holeaxis) );

      G4RotationMatrix *holerot = new G4RotationMatrix;  
      //if( !(ix == 0 && iy == 0 ) ) {  //removing the rotation of the holes in the plate. 10.29.20
      //holerot->rotate( -rotationangle, rotationaxis );
      //}
      
      //G4double platecenter_z = 68.8976*inch + bbsieved/2.0; //Assuming old placement
      G4double platecenter_z = 54.0157*inch + bbsieved/2.0; //Accounting for BB shift to 1.75m

      
      G4ThreeVector origin(0,0,-platecenter_z );
      G4double holedist = platecenter_z * sqrt(1.0 + pow(tan(xangle),2)+pow(tan(yangle),2) ); 
      
      G4ThreeVector holecenterpos = origin + holedist * holeaxis;
      
      cout << "(row,col,holeposx,holeposy,holeposz,rotationangle)=("
	   << iy << ", " << ix << ", " << holecenterpos.x()/inch << ", "
	   << holecenterpos.y()/inch << ", "
	   << holecenterpos.z()/inch << ", "
	   << rotationangle/deg << ")" << endl;

      cout << "rotation matrix = " << endl;
      holerot->print(cout);
      
      G4SubtractionSolid *NextCut;
      if( first ){
	first = false;
	if( !(ix==0&&iy==0) ){
	  NextCut = new G4SubtractionSolid( "NextCut", Fullplate, Hole, holerot, holecenterpos );
	} else { //it should be impossible to end up here:
	  NextCut = new G4SubtractionSolid( "NextCut", Fullplate, HoleCenter, holerot, holecenterpos );
	}
      } else {
	if( (iy==-3&&ix==-3)||(iy==0&&ix==0)||(iy==2&&ix==2) ){ //Cutting smaller holes
	  NextCut = new G4SubtractionSolid( "NextCut", SievePlateCut, HoleCenter, holerot, holecenterpos );
	}else{
	  if( iy==-3&&ix==1){
	    cout << "Skipping hole at y = -3 and x = 1" << endl; //Leaving one hole out
	  }else{
	    NextCut = new G4SubtractionSolid( "NextCut", SievePlateCut, Hole, holerot, holecenterpos );
	  }
	}
      }

      SievePlateCut = NextCut;
      SievePlateCut->SetName("NewBBSievePlateCut");
    }
  }

  G4cout << "finished holes..." << endl;
  
  G4LogicalVolume *NewBBSievePlate_log = new G4LogicalVolume(SievePlateCut, GetMaterial("Lead"), "New_BB_SievePlate_log");

  //For BB, plate is oriented along z axis and offset is passed into function.
  //G4double SBSsieve_dist = offset_z - ThickPlate_z/2.0;
  //Placement:
  //G4ThreeVector sievepos = SBSsieve_dist * G4ThreeVector( -sin(f48D48ang), 0, cos(f48D48ang) );
  //G4RotationMatrix *sieverot = new G4RotationMatrix;
  //sieverot->rotateY(f48D48ang);

  new G4PVPlacement(0, pos, NewBBSievePlate_log, "NewBBSievePlate_phys", motherlog, false, 0);
  
  G4cout << "Sieve plate finished" << G4endl;
}

 //Sieve slit designed by Holly S. - angled holes y, modified spacing. For optics studies 11.5.20
void G4SBSEArmBuilder::MakeFourthBBSieveSlit(G4LogicalVolume *motherlog, G4ThreeVector pos)
{
  printf("Building new BB sieve slit without angled holes...\n");
  
  //Plate
  //Plate dims - Assuming dimensions of box from extremes of holes and previous depth.
  G4double inch = 2.54*cm;

  G4double bbsievew = 14.75*inch;
  G4double bbsieveh = 27.50*inch;
  G4double bbsieved = 1.50*inch;
  
  //BigBite Sieve Slit is a box with holes:
  G4Box *Fullplate = new G4Box("Fullplate", bbsievew/2.0, bbsieveh/2.0, bbsieved/2.0 );

  //Hole Dimensions
  G4double HoleCenter_r = 0.125*inch;
  G4double Hole_r = 0.25*inch;
  G4double Hole_z = 3.0*inch; //Large enough to leave no solid volume at extreme displacements from center
  
  //Hole solids
  G4Tubs *Hole = new G4Tubs("Hole", 0, Hole_r, Hole_z, 0.0*deg, 360.0*deg);
  G4Tubs *HoleCenter = new G4Tubs("HoleCenter", 0, HoleCenter_r, Hole_z, 0.0*deg, 360.0*deg);

  G4double y_space = 1.25*inch;
  G4double x_space = 7.0/8.0*inch;

  G4double xangle = 0.0*deg;

  //G4SubtractionSolid *NextCut;
  G4SubtractionSolid *SievePlateCut;
  
  G4bool first = true;

  //Cut holes in the plate:

  cout << "cutting holes" << endl;
  
  
  for(G4int iy=-9; iy<=9; iy++ ){
    for(G4int ix=-6; ix<=6; ix++ ){

      G4double yangle = atan((iy*y_space)/(1.172*m)); //Computed from sieve distance at 1.55m

      //G4SubtractionSolid 
      G4RotationMatrix *holerot = new G4RotationMatrix;   
      holerot->rotateX(yangle); //May safely simplify since xangle is zero
      holerot->rotateY(xangle); 
      
      G4ThreeVector holecenterpos(ix*x_space, iy*y_space, 0);
      
      cout << "(row,col,holeposy,holeposx,holeposz,rotationy,rotationx)=("
	   << iy << ", " << ix << ", " << holecenterpos.x()/inch << ", "
	   << holecenterpos.y()/inch << ", "
	   << holecenterpos.z()/inch << ", "
	   << xangle/deg << ", "
	   << yangle/deg << ")" << endl;
      
      G4SubtractionSolid *NextCut;
      if( first ){
	first = false;
	if( !(ix==0&&iy==0) ){
	  NextCut = new G4SubtractionSolid( "NextCut", Fullplate, Hole, holerot, holecenterpos );
	} else { //it should be impossible to end up here:
	  NextCut = new G4SubtractionSolid( "NextCut", Fullplate, HoleCenter, holerot, holecenterpos );
	}
      } else {
	if( (iy==-3&&ix==-3)||(iy==0&&ix==0)||(iy==2&&ix==2) ){ //Cutting smaller holes
	  NextCut = new G4SubtractionSolid( "NextCut", SievePlateCut, HoleCenter, holerot, holecenterpos );
	}else{
	  if( iy==-3&&ix==1){
	    cout << "Skipping hole at y = -3 and x = 1" << endl; //Leaving one hole out
	  }else{
	    NextCut = new G4SubtractionSolid( "NextCut", SievePlateCut, Hole, holerot, holecenterpos );
	  }
	}
      }

      SievePlateCut = NextCut;
      SievePlateCut->SetName("NewBBSievePlateCut");
    }
  }

  cout << "finished holes..." <<endl;
  G4cout << "finished holes..." << endl;
  
  G4LogicalVolume *NewBBSievePlate_log = new G4LogicalVolume(SievePlateCut, GetMaterial("Lead"), "New_BB_SievePlate_log");

  new G4PVPlacement(0, pos, NewBBSievePlate_log, "NewBBSievePlate_phys", motherlog, false, 0);
  
  cout << "Sieve plate finished" << G4endl;
  G4cout << "Sieve plate finished" << G4endl;
}

void G4SBSEArmBuilder::MakeHallCGEM(G4LogicalVolume *motherlog){
  G4SBSTrackerBuilder trackerbuilder(fDetCon);
  
  //This routine creates and positions GEM plane in Hall
  
  //---Hall C GEM-------//
  G4double z0 = fBBdist;//m
  
  G4Box* HCGEMBox =  new G4Box("HCGEMBox", 10.0*cm, 10.0*cm, 2.0*cm);
  
  G4LogicalVolume* HCGEMlog = new G4LogicalVolume(HCGEMBox, GetMaterial("Air"),
						  "HCGEMLog", 0, 0, 0);
  HCGEMlog->SetVisAttributes(G4VisAttributes::Invisible);
  G4RotationMatrix *HCGEMrot = new G4RotationMatrix();
  
  HCGEMrot->rotateY(-fBBang);
  new G4PVPlacement(HCGEMrot, G4ThreeVector(z0*sin(fBBang), 0.0, z0*cos(fBBang)),
		    HCGEMlog, "HCGEMPhysical", motherlog, 0,false,0);
  
  G4int ngem = 1;
  vector<double> gemz, gemw, gemh;
  gemz.resize(ngem);
  gemw.resize(ngem);
  gemh.resize(ngem);
  
  gemz[0] = 0.0;;
  gemw[0] = 15.36*cm;
  gemh[0] = 15.36*cm;
  
  G4RotationMatrix *rot_identity = new G4RotationMatrix;
  
  trackerbuilder.BuildComponent(HCGEMlog, rot_identity, G4ThreeVector( 0.0, 0.0, 0.0 ), 1, gemz, gemw, gemh, "Earm/HCGEM" );
  
}
 
