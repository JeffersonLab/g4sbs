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
#include "G4GenericTrap.hh"
#include "G4SBSRICHSD.hh"
#include "G4SBSGlobalField.hh"

#include "G4OpticalSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4TwoVector.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SBSTrackerBuilder.hh"
#include "G4Sphere.hh"

#include "G4SBSCalSD.hh"
#include "G4Box.hh"
#include "sbstypes.hh"

#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"

#include <vector>
#include <map>

using namespace std;

G4SBSHArmBuilder::G4SBSHArmBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
    f48D48ang  = 39.4*deg;
    f48D48dist = 2.8*m;
    f48D48_fieldclamp_config = 1; //0 = No field clamps. 1 = GEp (default). 2 = BigBite experiments:

    fHCALdist  = 17.0*m;

    fRICHdist  = 15.0*m;

    f48D48depth = 1219.2*mm;
    f48D48width = 2324.1*mm;
    f48D48height = 3721.1*mm;

    fUseLocalField = false;
    fFieldStrength = 1.4*tesla;

    assert(fDetCon);
}

G4SBSHArmBuilder::~G4SBSHArmBuilder(){
}

void G4SBSHArmBuilder::BuildComponent(G4LogicalVolume *worldlog){
  Exp_t exptype = fDetCon->fExpType;


    // All three types of experiments have a 48D48 magnet:
   Make48D48(worldlog, f48D48dist + f48D48depth/2 );

    //--------------- HCAL --------------------------
    //All the experiments use HCAL:

    G4double HCAL_vertical_offset = 0.0*cm; //Neutron/SIDIS experiments have no vertical offset for HCAL (in Neutron case because it is detecting neutrons, which don't bend in a magnetic field, and in SIDIS case because we are detecting +/- charged hadrons simultaneously, want to have symmetric acceptance).
    if( exptype == kGEp ) HCAL_vertical_offset = 49.7*cm; //A number like this, which represents a positioning offset, shouldn't be hard-coded!


    MakeHCAL( worldlog, HCAL_vertical_offset );

    //The SIDIS experiment uses a  RICH in SBS:
    //--------- RICH (experimental): -------------------------
    if( exptype == kSIDISExp ) //SIDIS experiment requires a RICH detector and a tracker for SBS: 
    {
      //Let's make a simple tracker: 5 planes of GEMs, equally spaced in z, separation in z between planes of 10 cm. Then total length of tracker is ~50 cm + about 1.6 cm
      G4double SBStracker_dist = fRICHdist - 0.3*m;
      G4ThreeVector SBStracker_pos( SBStracker_dist * sin( f48D48ang ), 0.0, SBStracker_dist * cos( f48D48ang ) );

	G4RotationMatrix *SBStracker_rot_I = new G4RotationMatrix(G4RotationMatrix::IDENTITY);

	//Just a test:
	//SBStracker_rot_I->rotateY( 14.0*deg );

	G4RotationMatrix *SBStracker_rot = new G4RotationMatrix;
	SBStracker_rot->rotateY( -f48D48ang );

      G4Box *SBStracker_box = new G4Box("SBStracker_box", 32.0*cm, 102.0*cm, 22.0*cm );

      G4LogicalVolume *SBStracker_log = new G4LogicalVolume( SBStracker_box, GetMaterial("Air"), "SBStracker_log" );
	
      //For consistency with BigBite, place "Tracker" volume before GEMs are placed in it:
      new G4PVPlacement( SBStracker_rot, SBStracker_pos, SBStracker_log, "SBStracker_phys", worldlog, false, 0 );
	
      int ngems_SBStracker = 5;
      vector<double> zplanes_SBStracker, wplanes_SBStracker, hplanes_SBStracker;

      G4double zspacing_SBStracker = 10.0*cm;
      G4double zoffset_SBStracker = -20.0*cm;

      for(int i=0; i<ngems_SBStracker; i++ ){
	zplanes_SBStracker.push_back( zoffset_SBStracker + i*zspacing_SBStracker );
	wplanes_SBStracker.push_back( 60.0*cm );
	hplanes_SBStracker.push_back( 200.0*cm );
      }

      G4SBSTrackerBuilder trackerbuilder(fDetCon);

      (fDetCon->TrackerArm)[fDetCon->TrackerIDnumber] = kHarm; //H arm is "1"

      trackerbuilder.BuildComponent( SBStracker_log, SBStracker_rot_I, G4ThreeVector(0,0,0), 
				     ngems_SBStracker, zplanes_SBStracker, wplanes_SBStracker, hplanes_SBStracker, (fDetCon->TrackerIDnumber)++ );

      MakeRICH( worldlog );

      SBStracker_log->SetVisAttributes(G4VisAttributes::Invisible);
    }
    //---------------------------------------------------------
    if( exptype == kGEp ) //Subsystems unique to the GEp experiment include FPP and BigCal:
    {
	//Let's make a box and then put the FPP in it:
	//Define the rotation matrix for the FPP (pitch angle of 5 deg relative to vertical): 
	G4double sbsboxpitch = 5.0*deg;
	G4RotationMatrix *SBS_FPP_rm = new G4RotationMatrix;
	SBS_FPP_rm->rotateY( -f48D48ang );
	SBS_FPP_rm->rotateX( sbsboxpitch );

	//FPP box: 
	double sbsdepth  = 3.0*m;
	double sbswidth  = 2.0*m;
	double sbsheight = 2.1*m;

	double sbsr = fHCALdist-4.106*m + sbsheight*sin(sbsboxpitch)/2+sbsdepth/2;

	G4Box *sbsbox = new G4Box("sbsbox", sbswidth/2.0, sbsheight/2.0, sbsdepth/2.0 );
	G4LogicalVolume* sbslog = new G4LogicalVolume(sbsbox, GetMaterial("Air"), "sbslog");

	sbslog->SetVisAttributes( G4VisAttributes::Invisible );
	//Now position and orient the FPP "box":
	new G4PVPlacement(SBS_FPP_rm, G4ThreeVector(sbsr*sin(f48D48ang), (sbsr-f48D48dist)*sin(sbsboxpitch), sbsr*cos(f48D48ang) ), sbslog,
		"sbsphys", worldlog, false, 0, false);

	G4RotationMatrix *rot_I = new G4RotationMatrix;

	double detoffset = 0.05*m - sbsdepth/2.0;

	MakeFPP( sbslog, rot_I, G4ThreeVector( 0.0, 0.0, detoffset) );
    }


}


void G4SBSHArmBuilder::Make48D48( G4LogicalVolume *worldlog, double r48d48 ){
    G4String name;

    double bigcoilwidth = 214.5*mm;
    double bigcoilheight = 263.7*mm;

    double notchdepth = 25*cm;

    G4Box *biggap  = new G4Box("biggap",  469.9*mm/2+0.1*mm, 187.*cm/2.-bigcoilheight,  f48D48depth/2+0.1*mm);

    std::vector<G4TwoVector> bigpoly;
    bigpoly.push_back( G4TwoVector(-f48D48width/2.0,  f48D48depth/2.0 ));
    bigpoly.push_back( G4TwoVector(-f48D48width/2.0, -f48D48depth/2.0 ));
    bigpoly.push_back( G4TwoVector( f48D48width/2.0, -f48D48depth/2.0 ));
    bigpoly.push_back( G4TwoVector( f48D48width/2.0, f48D48depth/2.0 - notchdepth*sqrt(2.0) ));
    bigpoly.push_back( G4TwoVector( f48D48width/2.0 - notchdepth*sqrt(2.0), f48D48depth/2.0  ));

    G4ExtrudedSolid *bigbox_ext = new G4ExtrudedSolid("bigbox_ext", bigpoly, f48D48height/2.0, G4TwoVector(), 1.0, G4TwoVector(), 1.0);
    G4RotationMatrix *bigboxrm = new G4RotationMatrix;
    bigboxrm->rotateY(-f48D48ang);
    bigboxrm->rotateX( -90.*deg);
    bigboxrm->rotateZ( 180.*deg);

    G4RotationMatrix *bigboxaddrm = new G4RotationMatrix;
    bigboxaddrm->rotateZ( -180.*deg);
    bigboxaddrm->rotateX( 90.*deg);

    //moved definition of clamp gaps to field clamp method:
    // G4Box *bclampgap  = new G4Box("bclampgap",  23.*cm, 65.*cm,  12.*cm/2.);
    // G4Box *fclampgap  = new G4Box("fclampgap",  11.*cm, 35.*cm,  12.*cm/2.);

    G4SubtractionSolid* bigbase = new G4SubtractionSolid("bigbase", bigbox_ext, biggap, bigboxaddrm, G4ThreeVector());

    double coilgapwidth = (60.*cm - bigcoilheight)*2;
    double coilgapheight = 160*cm-bigcoilheight;

    G4Box *bigcoilbase = new G4Box("bigcoilbase", (bigcoilheight+coilgapwidth/2.0)/2.0, bigcoilheight+coilgapheight/2, bigcoilwidth/2.0);
    G4Box *bigcoilgap = new G4Box("bigcoilgap", coilgapwidth/4.0+1.0*mm, coilgapheight/2, bigcoilwidth/2.0+0.1*mm);

    //  double coilspace = 6.63*mm;
    double coilspace = 20.63*mm;



    std::vector<G4TwoVector> woundpoly;
    woundpoly.push_back( G4TwoVector(0.0,  -f48D48depth/2.0 -coilspace ));
    woundpoly.push_back( G4TwoVector(0.0, -f48D48depth/2.0 -coilspace-bigcoilwidth));
    woundpoly.push_back( G4TwoVector( f48D48width/2.0+coilspace+bigcoilwidth, -f48D48depth/2.0 -coilspace-bigcoilwidth));
    woundpoly.push_back( G4TwoVector( f48D48width/2.0+coilspace+bigcoilwidth, f48D48depth/2.0 - notchdepth*sqrt(2.0) ));

    woundpoly.push_back( G4TwoVector( f48D48width/2.0+coilspace+bigcoilwidth  - 2.0*(coilspace+bigcoilwidth)*sin(pi/8.)*sin(pi/8.) , 
		f48D48depth/2.0 - notchdepth*sqrt(2.0) + 2.0*(coilspace+bigcoilwidth)*sin(pi/8)*cos(pi/8) ));

    woundpoly.push_back( G4TwoVector( f48D48width/2.0- notchdepth*sqrt(2.0) + 2.0*(coilspace+bigcoilwidth)*sin(pi/8)*cos(pi/8) , 
		f48D48depth/2.0 + coilspace+bigcoilwidth  - 2.0*(coilspace+bigcoilwidth)*sin(pi/8)*sin(pi/8) ));

    woundpoly.push_back( G4TwoVector( f48D48width/2.0 - notchdepth*sqrt(2.0), f48D48depth/2.0 +coilspace+bigcoilwidth ));

    ////
    woundpoly.push_back( G4TwoVector(0.0,  f48D48depth/2.0 +coilspace+bigcoilwidth));
    woundpoly.push_back( G4TwoVector(0.0,  f48D48depth/2.0 +coilspace ));

    woundpoly.push_back( G4TwoVector( f48D48width/2.0 - notchdepth*sqrt(2.0), f48D48depth/2.0 +coilspace ));

    // arc here
    woundpoly.push_back( G4TwoVector( f48D48width/2.0+coilspace  - 2.0*(coilspace)*sin(pi/8.)*sin(pi/8.) , 
		f48D48depth/2.0 - notchdepth*sqrt(2.0) + 2.0*(coilspace)*sin(pi/8)*cos(pi/8) ));

    woundpoly.push_back( G4TwoVector( f48D48width/2.0- notchdepth*sqrt(2.0) + 2.0*(coilspace)*sin(pi/8)*cos(pi/8) , 
		f48D48depth/2.0 + coilspace - 2.0*(coilspace)*sin(pi/8)*sin(pi/8) ));

    woundpoly.push_back( G4TwoVector( f48D48width/2.0+coilspace, f48D48depth/2.0 - notchdepth*sqrt(2.0) ));

    woundpoly.push_back( G4TwoVector( f48D48width/2.0+coilspace, f48D48depth/2.0 +coilspace - notchdepth*sqrt(2.0) ));
    woundpoly.push_back( G4TwoVector( f48D48width/2.0+coilspace, -f48D48depth/2.0 -coilspace));

    G4ExtrudedSolid *woundcoil_ext = new G4ExtrudedSolid("woundcoil_ext", woundpoly, bigcoilheight/2.0, G4TwoVector(), 1.0, G4TwoVector(), 1.0);

    // Pull out left side of gap
    G4SubtractionSolid *bigcoil = new G4SubtractionSolid("bigcoil", bigcoilbase, bigcoilgap, 0, G4ThreeVector(-coilgapwidth/4.0, 0.0, 0.0) );

    //  double coilfrontback = 150*cm;
    double coilfrontback = f48D48depth+2.0*coilspace+bigcoilwidth;

    G4Box *bigcoilthr = new G4Box("bigcoilthr", bigcoilwidth,  bigcoilheight/2,  coilfrontback/2.0 );

    // Sum together coils


    // Sum together base iron plus coils

    G4UnionSolid* big48d48;

    G4Box *bigbeamslot = new G4Box("bigbeamslot",  f48D48width/2, 15.5*cm, 2.0*m ); 

    big48d48 = new G4UnionSolid("big48d48_1", bigbase, bigcoilthr, bigboxaddrm, 
	    G4ThreeVector(0.0, 0.0, (coilgapheight+bigcoilheight)/2.0));
    big48d48 = new G4UnionSolid("big48d48_2", big48d48, bigcoilthr, bigboxaddrm, 
	    G4ThreeVector(0.0, 0.0, -(coilgapheight+bigcoilheight)/2.0));

    big48d48 = new G4UnionSolid("big48d48_3", big48d48, bigcoil, bigboxaddrm, 
	    G4ThreeVector(-(bigcoilheight+coilgapwidth/2.0)/2.0-1.0*mm, coilfrontback/2.0, 0.0));
    big48d48 = new G4UnionSolid("big48d48_4", big48d48, bigcoil, bigboxaddrm, 
	    G4ThreeVector(-(bigcoilheight+coilgapwidth/2.0)/2.0-1.0*mm, -coilfrontback/2.0, 0.0));

    big48d48 = new G4UnionSolid("big48d48_5", big48d48, woundcoil_ext, 0, 
	    G4ThreeVector( 1.0*mm, 0.0,  coilgapheight/2.+bigcoilheight/2.0));
    big48d48 = new G4UnionSolid("big48d48_6", big48d48, woundcoil_ext, 0, 
	    G4ThreeVector( 1.0*mm, 0.0,  -coilgapheight/2.-bigcoilheight/2.0));


    //  Cut out slot - from magnet center to inside of cut is ~35cm
    G4SubtractionSolid *big48d48_wslot = new G4SubtractionSolid("big48d48_5", big48d48, bigbeamslot, bigboxaddrm, 
	    G4ThreeVector( f48D48width/2+35*cm, 0.0, 0.0) );

    G4LogicalVolume *big48d48Log=new G4LogicalVolume(big48d48_wslot, GetMaterial("Fer"),
	    "b48d48Log", 0, 0, 0);

    if( fDetCon->fTotalAbs ){
	big48d48Log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
    }

    G4RotationMatrix *bigrm = new G4RotationMatrix;
    bigrm->rotateY(-f48D48ang);


    new G4PVPlacement(bigboxrm, 
	    G4ThreeVector(-r48d48*sin(-f48D48ang), 0.0, r48d48*cos(-f48D48ang)),
	    big48d48Log, "big48d48Physical", worldlog, 0,false,0);

    G4LogicalVolume *bigfieldLog=new G4LogicalVolume(biggap, GetMaterial("Air"),
	    "bigfieldLog", 0, 0, 0);

    // Associate magnetic field with gap

    double sign = 1.0;
    if( fDetCon->fGlobalField->fInverted ) sign = -1.0;
    G4UniformMagField* magField
	= new G4UniformMagField(G4ThreeVector(sign*fFieldStrength*cos(f48D48ang), 0.0, -sign*fFieldStrength*sin(f48D48ang)));

    G4FieldManager *bigfm = new G4FieldManager(magField);
    bigfm->SetDetectorField(magField);
    bigfm->CreateChordFinder(magField);

    if( fUseLocalField ){
	bigfieldLog->SetFieldManager(bigfm,true);
    }


    new G4PVPlacement(bigrm, 
	    G4ThreeVector(-r48d48*sin(-f48D48ang), 0.0, r48d48*cos(-f48D48ang)),
	    bigfieldLog, "bigfieldPhysical", worldlog, 0,false,0);


    if( fDetCon->fExpType == kGEp ){
	// Addtional iron inside the field region

	std::vector<G4TwoVector> leftverts;

	leftverts.push_back( G4TwoVector( -12*cm, -45*cm ) );
	leftverts.push_back( G4TwoVector( -12*cm,  45*cm ) );
	leftverts.push_back( G4TwoVector( -23.5*cm, 45*cm ) );
	leftverts.push_back( G4TwoVector( -23.5*cm, -45*cm ) );

	leftverts.push_back( G4TwoVector( -17.75*cm, -61*cm ) );
	leftverts.push_back( G4TwoVector( -17.75*cm, 61*cm ) );
	leftverts.push_back( G4TwoVector( -23.5*cm, 61*cm ) );
	leftverts.push_back( G4TwoVector( -23.5*cm, -61*cm ) );

	std::vector<G4TwoVector> rightverts;

	rightverts.push_back( G4TwoVector( 23.5*cm, -45*cm ) );
	rightverts.push_back( G4TwoVector( 23.5*cm, 45*cm ) );
	rightverts.push_back( G4TwoVector( 12*cm,  45*cm ) );
	rightverts.push_back( G4TwoVector( 12*cm, -45*cm ) );

	rightverts.push_back( G4TwoVector( 23.5*cm, -61*cm ) );
	rightverts.push_back( G4TwoVector( 23.5*cm, 61*cm ) );
	rightverts.push_back( G4TwoVector( 17.75*cm, 61*cm ) );
	rightverts.push_back( G4TwoVector( 17.75*cm, -61*cm ) );


	double slabdepth= 61*cm;
	G4GenericTrap *leftslab = new G4GenericTrap("leftslab", slabdepth, leftverts );
	G4GenericTrap *rightslab = new G4GenericTrap("rightslab", slabdepth, rightverts );

	G4LogicalVolume *leftslabLog=new G4LogicalVolume(leftslab, GetMaterial("Fer"), "leftslabLog", 0, 0, 0);
	G4LogicalVolume *rightslabLog=new G4LogicalVolume(rightslab, GetMaterial("Fer"), "rightslabLog", 0, 0, 0);

	if( fDetCon->fTotalAbs ){
		leftslabLog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
		rightslabLog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
	}

	new G4PVPlacement(bigrm, 
	    G4ThreeVector(-r48d48*sin(-f48D48ang), 0.0, r48d48*cos(-f48D48ang)),
	    leftslabLog, "leftslabPhysical", worldlog, 0,false,0);
	new G4PVPlacement(bigrm, 
	    G4ThreeVector(-r48d48*sin(-f48D48ang), 0.0, r48d48*cos(-f48D48ang)),
	    rightslabLog, "rightslabPhysical", worldlog, 0,false,0);

	G4VisAttributes * slabVisAtt
	    = new G4VisAttributes(G4Colour(1.0,0.1,0.0));

	leftslabLog->SetVisAttributes(slabVisAtt);
	rightslabLog->SetVisAttributes(slabVisAtt);
    }

    // Clamps

    // The positioning and acceptance gaps were taken directly from CAD
    // The position widths for the beam pipe holes are fudged around so
    // that it doesn't interfere with the beam pipe
    //
    if( f48D48_fieldclamp_config != 0 ){

	MakeSBSFieldClamps(worldlog);

    }

    bigfieldLog->SetVisAttributes(G4VisAttributes::Invisible);
    //  backclampLog->SetVisAttributes(G4VisAttributes::Invisible);
    //  frontclampLog->SetVisAttributes(G4VisAttributes::Invisible);
}

void G4SBSHArmBuilder::MakeSBSFieldClamps( G4LogicalVolume *motherlog ){
  
  if( f48D48_fieldclamp_config == 1 ){ //BigBite:

    double clampdepth = 10.*cm;
    double clampoffset = 45*cm/2;

    G4Box *bclampgap  = new G4Box("bclampgap",  23.*cm, 65.*cm,  12.*cm/2.);

    double fclampgapwidth = (11+80)*cm;
    G4Box *fclampgap  = new G4Box("fclampgap",  fclampgapwidth/2, 35.*cm,  12.*cm/2.);

    double frontclampwidth = 160*cm;
    double frontclampheight = 3*m;

    double frontclampz = -100*cm + clampdepth/2.0;
    double backclampz  =  100*cm - clampdepth/2.0;

    G4Box *frontclampbase  = new G4Box("frontclampbase", frontclampwidth/2, frontclampheight/2,  clampdepth/2);
    G4SubtractionSolid *frontclamp = new G4SubtractionSolid("frontclamp1", frontclampbase, fclampgap, 0,
	    G4ThreeVector( -(frontclampwidth-fclampgapwidth)/2 -0.1*mm, 0,0 ) );

    // G4Box *frontclampbeamhole  = new G4Box("frontclampbeamhole", 20.*cm/2, 20.*cm/2,  clampdepth/2+2*cm);
    //  frontclamp = new G4SubtractionSolid("frontclamp2", frontclamp, frontclampbeamhole, 0, G4ThreeVector(-55*cm, 0, 0) );
    //  frontclamp = new G4SubtractionSolid("frontclamp3", frontclamp, frontclampbeamhole, 0, G4ThreeVector(-28*cm, 0, 0) );

    G4Box *frontclampbeamhole  = new G4Box("frontclampbeamhole", 50.*cm/2, 14.*cm/2,  clampdepth/2+2*cm);
    frontclamp = new G4SubtractionSolid("frontclamp2", frontclamp, frontclampbeamhole, 0, G4ThreeVector(-55*cm+clampoffset, 0, 0) );

//    G4Box *frontclampecalhole  = new G4Box("frontclampecalhole", 100.*cm/2, 236.*cm/4.,  clampdepth/2+2*cm);
//    frontclamp = new G4SubtractionSolid("frontclamp3", frontclamp, frontclampecalhole, 0, G4ThreeVector(-120*cm+clampoffset, 0, 0) );

    // This awful extrusion //

    double extang= 16*deg;
    std::vector<G4TwoVector> faceverts;
    faceverts.push_back( G4TwoVector( -80.0*cm, -clampdepth/2.0-25.1*cm ) );
    faceverts.push_back( G4TwoVector( -80.0*cm, -clampdepth/2.0-25.1*cm + clampdepth*cos(extang)) );
    faceverts.push_back( G4TwoVector( -11.0*cm, -clampdepth/2 + clampdepth*cos(extang) ) );
    faceverts.push_back( G4TwoVector( -11.0*cm, -clampdepth/2 ) );

    faceverts.push_back( G4TwoVector( -80.0*cm, -clampdepth/2.0-25.1*cm ) );
    faceverts.push_back( G4TwoVector( -80.0*cm, -clampdepth/2.0-25.1*cm + clampdepth*cos(extang)) );
    faceverts.push_back( G4TwoVector( -11.0*cm, -clampdepth/2 + clampdepth*cos(extang) ) );
    faceverts.push_back( G4TwoVector( -11.0*cm, -clampdepth/2 ) );

    std::vector<G4TwoVector> topverts;
    topverts.push_back( G4TwoVector( -80.0*cm, -clampdepth/2.0-25.1*cm ) );
    topverts.push_back( G4TwoVector( -80.0*cm, -clampdepth/2.0 + clampdepth*cos(extang)) );
    topverts.push_back( G4TwoVector( -11.0*cm, -clampdepth/2.0 + clampdepth*cos(extang)) );
    topverts.push_back( G4TwoVector( -11.0*cm, -clampdepth/2 ) );

    topverts.push_back( G4TwoVector( -80.0*cm, -clampdepth/2.0-25.1*cm ) );
    topverts.push_back( G4TwoVector( -80.0*cm, -clampdepth/2.0 + clampdepth*cos(extang)) );
    topverts.push_back( G4TwoVector( -11.0*cm, -clampdepth/2.0 + clampdepth*cos(extang)) );
    topverts.push_back( G4TwoVector( -11.0*cm, -clampdepth/2 ) );

    double extheight = 90*cm;
    G4GenericTrap *extface = new G4GenericTrap("extface", extheight/2-clampdepth-0.5*mm, faceverts );
    // Cut out hole for face at 16deg
    G4RotationMatrix *faceholerm = new G4RotationMatrix();
    faceholerm->rotateZ(-f48D48ang);
    faceholerm->rotateX(90*deg);

    G4Tubs *facehole;
   
    if( fDetCon->fTargType == kLH2 || fDetCon->fTargType == kLD2 ){
	facehole = new G4Tubs("facehole", 0.0, 5.5*cm, 40.*cm, 0, 360*deg);
    } else {
	facehole = new G4Tubs("facehole", 0.0, 8*cm, 40.*cm, 0, 360*deg);
    }
    G4SubtractionSolid *extface_whole = NULL;
   
    extface_whole = new G4SubtractionSolid("extface_whole", extface, facehole, faceholerm, G4ThreeVector(  tan(-f48D48ang)*(f48D48dist+frontclampz + f48D48depth/2.0), 0.0, 0.0));

    G4GenericTrap *exttop  = new G4GenericTrap("exttop", clampdepth/2, topverts );

    G4RotationMatrix *extrot = new G4RotationMatrix();
    extrot->rotateX(-90*deg);


    G4UnionSolid *frontclampun = new G4UnionSolid("frontclamp3", frontclamp, exttop, extrot, G4ThreeVector(0,  (extheight-clampdepth+2*mm)/2, 0 ));
    frontclampun = new G4UnionSolid("frontclamp4", frontclampun, exttop, extrot, G4ThreeVector(0, -(extheight-clampdepth+2*mm)/2, 0 ));


    G4LogicalVolume *frontclampLog=new G4LogicalVolume(frontclampun, GetMaterial("Fer"), "frontclampLog", 0, 0, 0);

    G4LogicalVolume *frontextfaceLog= NULL;
    if( fDetCon->fExpType == kGEp || fDetCon->fExpType == kNeutronExp){
	frontextfaceLog = new G4LogicalVolume(extface_whole, GetMaterial("Fer"), "frontextfaceLog", 0, 0, 0);
    } else {
	frontextfaceLog = new G4LogicalVolume(extface, GetMaterial("Fer"), "frontextfaceLog", 0, 0, 0);
    }

    if( fDetCon->fTotalAbs ){
	frontclampLog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
	frontextfaceLog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
    }

    double backclampheight = 292.0*cm;

    double backclampwidth = (150+117)*cm;
    G4Box *backclampbase  = new G4Box("backclampbase", backclampwidth/2, backclampheight/2,  clampdepth/2);

    double backclampaddheight = 55.0*cm;
    G4Box *backclampadd  = new G4Box("backclampadd", backclampwidth/2, backclampaddheight/2,  clampdepth);
    G4UnionSolid *backclampfull = new G4UnionSolid("backclampfull1", backclampbase, backclampadd, 0,
	    G4ThreeVector( 0.0, backclampheight/2.0 - backclampaddheight/2.0, clampdepth/2.0-0.1*mm ));
    backclampfull = new G4UnionSolid("backclampfull1", backclampfull, backclampadd, 0,
	    G4ThreeVector( 0.0, -backclampheight/2.0 + backclampaddheight/2.0, clampdepth/2.0 ));

    G4SubtractionSolid *backclamp = new G4SubtractionSolid("backclamp1", backclampfull, bclampgap, 0, 
	    G4ThreeVector(clampoffset, 0, 0) );

    G4Box *backclampbeamhole  = new G4Box("backclampbeamhole", (151-55)*cm/2., 51.*cm,  clampdepth/2+2*cm);
    backclamp = new G4SubtractionSolid("backclamp2", backclamp, backclampbeamhole, 0, G4ThreeVector(-backclampwidth/2+(151-55)*cm/2 - 0.1*mm, 0, 0) );

    G4LogicalVolume *backclampLog=new G4LogicalVolume(backclamp, GetMaterial("Fer"), "backclampLog", 0, 0, 0);

    if( fDetCon->fTotalAbs ){
	backclampLog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
    }


    G4RotationMatrix *rot = new G4RotationMatrix;
    rot->rotateY( -f48D48ang );

    double r48d48 = f48D48dist + 1219.2*mm/2.0;

    new G4PVPlacement(rot, 
	    G4ThreeVector(-(r48d48+frontclampz)*sin(-f48D48ang), 0.0, (r48d48+frontclampz)*cos(-f48D48ang)),
	    frontclampLog, "frontclampPhysical", motherlog, 0,false,0);

    G4RotationMatrix *rotextface = new G4RotationMatrix();
    rotextface->rotateY(-f48D48ang);
    rotextface->rotateX(-90*deg);

    // Face extrusion
    new G4PVPlacement(rotextface, 
	    G4ThreeVector(-(r48d48+frontclampz)*sin(-f48D48ang), 0.0, (r48d48+frontclampz)*cos(-f48D48ang)),
	    frontextfaceLog, "extfacePhysical", motherlog, 0,false,0);

    // Back clamp is GEp only?
    if( fDetCon->fExpType == kGEp ){
	new G4PVPlacement(rot, 
		G4ThreeVector(-(r48d48+backclampz)*sin(-f48D48ang)-cos(-f48D48ang)*clampoffset, 0.0, (r48d48+backclampz)*cos(-f48D48ang)-sin(-f48D48ang)*clampoffset),
		backclampLog, "backclampPhysical", motherlog, 0,false,0);
    }

    G4VisAttributes * clampVisAtt
	= new G4VisAttributes(G4Colour(0.8,1.0,0.4));

    frontclampLog->SetVisAttributes(clampVisAtt);
    frontextfaceLog->SetVisAttributes(clampVisAtt);
    backclampLog->SetVisAttributes(clampVisAtt);
  } else if( f48D48_fieldclamp_config == 2 ){ //GEp
   
    G4double FrontClamp_width = 62.99*2.54*cm;
    G4double FrontClamp_height = 118.11*2.54*cm;
    G4double FrontClamp_depth = 3.94*2.54*cm;

    G4double FrontClamp_notch_width = 35.83*2.54*cm;
    G4double FrontClamp_notch_height = 27.56*2.54*cm;

    G4Box *FrontClamp_Box = new G4Box( "FrontClamp_Box", FrontClamp_width/2.0, FrontClamp_height/2.0, FrontClamp_depth/2.0 );
    G4Box *FrontClamp_Notch = new G4Box( "FrontClamp_Notch", FrontClamp_notch_width/2.0 + 1.0*cm, FrontClamp_notch_height/2.0, FrontClamp_depth/2.0 + 1.0*cm );

    //Position notch so that the distance from the right edge of the clamp to the left edge of the notch equals notch width:
    //xnotch + wnotch/2 + 1 cm + wclamp/2 = wnotch --> xnotch = wnotch/2- 1 cm - wclamp/2
    G4double xnotch = FrontClamp_notch_width/2.0 - 1.0*cm - FrontClamp_width/2.0;

    G4SubtractionSolid *FrontClamp = new G4SubtractionSolid( "FrontClamp", FrontClamp_Box, FrontClamp_Notch, 0, G4ThreeVector( xnotch, 0.0, 0.0 ) );

    G4LogicalVolume *FrontClamp_log = new G4LogicalVolume( FrontClamp, GetMaterial("Fer"), "FrontClamp_log" );

    G4double FrontClamp_zoffset = 13.40*2.54*cm + FrontClamp_depth/2.0;

    G4double FrontClamp_r = f48D48dist - FrontClamp_zoffset;
    
    G4double FrontClamp_xshift = 14.22*2.54*cm; //x offset of left edge of front clamp relative to left edge of magnet

    G4double FrontClamp_xoffset = f48D48width/2.0 - FrontClamp_width/2.0 - FrontClamp_xshift; //net offset in position needed to align left edge of front clamp at correct position.

    G4RotationMatrix *clamp_rot = new G4RotationMatrix;

    clamp_rot->rotateY( -f48D48ang );

    new G4PVPlacement( clamp_rot, 
		       G4ThreeVector( FrontClamp_r*sin( f48D48ang ) + FrontClamp_xoffset*cos(f48D48ang), 0.0, FrontClamp_r*cos(f48D48ang) - FrontClamp_xoffset*sin(f48D48ang) ), 
		       FrontClamp_log, "FrontClamp_phys", motherlog, false, 0, false );
 
    ////REAR CLAMP:
    G4double RearClamp_width = 105.12*2.54*cm;
    G4double RearClamp_height = 114.96*2.54*cm;
    G4double RearClamp_depth = 5.91*2.54*cm;

    G4Box *RearClamp_Box = new G4Box("RearClamp_Box", RearClamp_width/2.0, RearClamp_height/2.0, RearClamp_depth/2.0 );
    
    G4double RearClamp_GapWidth = 18.11*2.54*cm;
    G4double RearClamp_GapHeight = 51.18*2.54*cm;
    G4double RearClamp_NotchWidth = 37.84*2.54*cm;
    G4double RearClamp_NotchHeight = 39.37*2.54*cm;

    G4double RearClamp_GapX = 36.61*2.54*cm; //This is the distance from the left edge of the rear clamp to the left edge of the opening.

    G4Box *RearClamp_Gap = new G4Box("RearClamp_Gap", RearClamp_GapWidth/2.0, RearClamp_GapHeight/2.0, RearClamp_depth/2.0 + 1.0*cm );

    G4SubtractionSolid *RearClamp_GapCutout = new G4SubtractionSolid( "RearClamp_GapCutout", RearClamp_Box, RearClamp_Gap, 0, 
								      G4ThreeVector( RearClamp_width/2.0 - RearClamp_GapX - RearClamp_GapWidth/2.0, 0.0, 0.0 ) );
    
    G4Box *RearClamp_Notch = new G4Box("RearClamp_Notch", RearClamp_NotchWidth/2.0+1.0*cm, RearClamp_NotchHeight/2.0, RearClamp_depth/2.0 + 1.0*cm );
    
    xnotch = RearClamp_NotchWidth/2.0 - 1.0*cm - RearClamp_width/2.0;

    G4SubtractionSolid *RearClamp = new G4SubtractionSolid( "RearClamp", RearClamp_GapCutout, RearClamp_Notch, 0, 
							    G4ThreeVector( xnotch, 0.0, 0.0 ) );

    G4LogicalVolume *RearClamp_log = new G4LogicalVolume( RearClamp, GetMaterial("Fer"), "RearClamp_log" );

    G4double RearClamp_zoffset = 11.43*2.54*cm + RearClamp_depth/2.0; 
    G4double RearClamp_xoffset = f48D48width/2.0 - RearClamp_width/2.0; 
    G4double RearClamp_r = f48D48dist + f48D48depth + RearClamp_zoffset;

    G4ThreeVector RearClamp_pos( RearClamp_r*sin(f48D48ang) + RearClamp_xoffset*cos(f48D48ang), 0.0, RearClamp_r*cos(f48D48ang) - RearClamp_xoffset * sin(f48D48ang) );

    new G4PVPlacement( clamp_rot, RearClamp_pos, RearClamp_log, "RearClamp_phys", motherlog, false, 0, false );

    G4VisAttributes * clampVisAtt
      = new G4VisAttributes(G4Colour(0.8,1.0,0.4));

    FrontClamp_log->SetVisAttributes(clampVisAtt);
    RearClamp_log->SetVisAttributes(clampVisAtt);
  }
}


void G4SBSHArmBuilder::MakeHCAL( G4LogicalVolume *motherlog, G4double VerticalOffset=0.0*cm ){
    //G4SDManager* SDman = G4SDManager::GetSDMpointer();

    // HCAL 
    double hcalheight = 330.0*cm;
    double hcalwidth  = 165.0*cm;
    double hcaldepth  = 101.0*cm;
    double hcalr = fHCALdist+hcaldepth/2.0;

    G4RotationMatrix *hcalrm = new G4RotationMatrix;
    hcalrm->rotateY(-f48D48ang);

    G4Box *hcalbox = new G4Box("hcalbox", hcalwidth/2.0, hcalheight/2.0, hcaldepth/2.0 );
    G4LogicalVolume* hcallog = new G4LogicalVolume(hcalbox, GetMaterial("Lead"), "hcallog");

    new G4PVPlacement(hcalrm, G4ThreeVector(hcalr*sin(f48D48ang), VerticalOffset, hcalr*cos(f48D48ang) ), hcallog,
	    "hcalphys", motherlog, false, 0, false);

    G4String HCALSDname = "G4SBS/HCAL";
    G4String HCALcolname = "HCALcol";
    G4SBSCalSD* HCalSD;

    G4SDManager *sdman = fDetCon->fSDman;

  if( !(HCalSD = (G4SBSCalSD*) sdman->FindSensitiveDetector(HCALSDname)) ){
    HCalSD = new G4SBSCalSD( HCALSDname, HCALcolname );
    sdman->AddNewDetector(HCalSD);
    fDetCon->SDlist[HCALSDname] = HCalSD;
  }

  sdman->AddNewDetector(HCalSD);
  hcallog->SetSensitiveDetector(HCalSD);
  hcallog->SetUserLimits(  new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );

  G4VisAttributes * hcalVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.6,0.0));

  hcallog->SetVisAttributes(hcalVisAtt);
}

void G4SBSHArmBuilder::MakeRICH( G4LogicalVolume *motherlog ){

    //*********************************************************************************************************************************//
    //                                  BEGIN GEOMETRY DEFINITION FOR SBS RICH COUNTER                                                 //
    //*********************************************************************************************************************************//


    //First, define a RICH box that will serve as the mother volume for the positioning of the RICH components relative to each other and 
    //as the containment volume for the C4F10 radiator gas:
    G4double RICHbox_dx=165.0*cm, RICHbox_dy=283.0*cm, RICHbox_dz=127.0*cm;

    G4Box *RICHbox = new G4Box( "RICHbox", RICHbox_dx/2.0, RICHbox_dy/2.0, RICHbox_dz/2.0 );
    G4LogicalVolume *RICHbox_log = new G4LogicalVolume( RICHbox, GetMaterial("C4F10_gas"), "RICHbox_log" );

    //We also want to define exterior walls of the RICH box as 1-inch thick aluminum: we will do a subtraction solid of RICHbox and RICHwalls 
    //as well as the entry windows...
    G4double RICHwall_dx = RICHbox_dx - 2.0*2.54*cm;
    G4double RICHwall_dy = RICHbox_dy - 2.0*2.54*cm;
    G4double RICHwall_dz = RICHbox_dz - 2.0*2.54*cm;
    G4Box *RICHwall = new G4Box( "RICHwall", RICHwall_dx/2.0, RICHwall_dy/2.0, RICHwall_dz/2.0 );

    G4SubtractionSolid *RICHbox_wall = new G4SubtractionSolid("RICHbox_wall", RICHbox, RICHwall );

    // The stacking of the tiles is defined by three parameters, nx, ny, and nz, the number of tiles along each dimension. x is assumed to be the horizontal direction, y the 
    // vertical direction, and z the nominal direction of particle motion. The design constraint is that nx*ny*nz <= 850:
    G4int nx_aero=5, ny_aero=17, nz_aero=5;

    //define tile dimensions: 
    G4double aero_dx = (11.4/2.0)*cm, aero_dy=(11.4/2.0)*cm, aero_dz=(1.13/2.0)*cm;

    G4Box *Aerogel_tile_solid = new G4Box("Aerogel_tile", aero_dx, aero_dy, aero_dz );  
    G4LogicalVolume *Aerogel_tile_log = new G4LogicalVolume( Aerogel_tile_solid, GetMaterial("Aerogel"), "Aerogel_tile_log" );

    //Assume 1 mil gap between tiles, filled with tedlar

    //Decide where we want center of aerogel coordinate system to be located. 
    //Assume that the geometric center of the combined aerogel box is at the origin of RICHbox

    G4double tilegap = 0.00254*cm; //gap between tiles in x and y:

    G4double Lx_aero = 2.0*nx_aero*aero_dx + tilegap*(nx_aero-1);
    G4double Ly_aero = 2.0*ny_aero*aero_dy + tilegap*(ny_aero-1);
    G4double Lz_aero = 2.0*nz_aero*aero_dz;

    G4double RICH_entrywindow_thick = 0.1*cm;
    G4double RICH_entrywindow_dz = RICH_entrywindow_thick/2.0;
    G4double RICH_entrywindow_dx = Lx_aero/2.0;
    G4double RICH_entrywindow_dy = Ly_aero/2.0;

    //1 mm-thick aluminum entry window for RICH, assumed to have same dimensions as aerogel tiles. 
    G4Box *RICH_entrywindow = new G4Box("RICH_entrywindow", RICH_entrywindow_dx, RICH_entrywindow_dy, RICH_entrywindow_dz );
    G4LogicalVolume *RICH_entrywindow_log = new G4LogicalVolume( RICH_entrywindow, GetMaterial("RICHAluminum"), "RICH_entrywindow_log" );

    G4double Aero_exit_dx = Lx_aero/2.0;
    G4double Aero_exit_dy = Ly_aero/2.0;
    G4double Aero_exit_dz = (0.32/2.0)*cm;
    //3.2 mm-thick UVT-lucite exit window for aerogel:
    G4Box *Aero_exit = new G4Box( "Aero_exit", Aero_exit_dx, Aero_exit_dy, Aero_exit_dz );
    G4LogicalVolume *Aero_exitwindow = new G4LogicalVolume( Aero_exit, GetMaterial("UVT_Lucite"), "Aero_exitwindow" );

    //1 mm-thick aluminum exit window for RICH, dimensions are up to us, but let's start with HERMES case:
    G4double RICH_exitwindow_thick = 0.1*cm;
    G4double RICH_exitwindow_dx = (59.0/2.0)*cm;
    G4double RICH_exitwindow_dy = (257.0/2.0)*cm;
    G4double RICH_exitwindow_dz = RICH_exitwindow_thick/2.0;

    G4Box *RICH_exit = new G4Box( "RICH_exit", RICH_exitwindow_dx, RICH_exitwindow_dy, RICH_exitwindow_dz );
    G4LogicalVolume *RICH_exitwindow = new G4LogicalVolume( RICH_exit, GetMaterial("RICHAluminum"), "RICH_exitwindow" );

    G4double aero_xoffset = 7.191*cm;

    G4double z0_entrywindow = -RICHbox_dz/2.0 + RICH_entrywindow_dz;
    G4double y0_entrywindow = 0.0;
    G4double x0_entrywindow = -RICHbox_dx/2.0 + aero_xoffset + Lx_aero/2.0;

    //Position entry and exit windows inside RICHbox:
    new G4PVPlacement( 0, 
	    G4ThreeVector( x0_entrywindow, y0_entrywindow, z0_entrywindow), 
	    RICH_entrywindow_log, 
	    "RICH_entrywindow_pv", 
	    RICHbox_log, 
	    false, 
	    0 );


    G4double x0_aeroexit = x0_entrywindow;
    G4double y0_aeroexit = 0.0;
    G4double z0_aeroexit = z0_entrywindow + RICH_entrywindow_dz + Lz_aero + Aero_exit_dz;

    //This is aerogel exit window.
    new G4PVPlacement( 0, 
	    G4ThreeVector( x0_aeroexit, y0_aeroexit, z0_aeroexit ),
	    Aero_exitwindow, 
	    "Aero_exitwindow_pv",
	    RICHbox_log, 
	    false,
	    0 );

    G4double x0_RICHexit = -RICHbox_dx/2.0 + aero_xoffset + RICH_exitwindow_dx;
    G4double y0_RICHexit = 0.0;
    G4double z0_RICHexit = RICHbox_dz/2.0 - RICH_exitwindow_dz;

    new G4PVPlacement( 0, 
	    G4ThreeVector( x0_RICHexit, y0_RICHexit, z0_RICHexit ),
	    RICH_exitwindow, 
	    "RICH_exitwindow_pv",
	    RICHbox_log, 
	    false, 
	    0 );

    //We need to define cutouts from the 1"-thick Aluminum box for the entry and exit windows:
    G4Box *RICHentry_cutout = new G4Box( "RICHentry_cutout", RICH_entrywindow_dx, RICH_entrywindow_dy, 10.0*cm );
    G4SubtractionSolid *RICHbox_wall_entrycut = new G4SubtractionSolid( "RICHbox_wall_entrycut", RICHbox_wall, RICHentry_cutout, 0, G4ThreeVector( x0_entrywindow, y0_entrywindow, -RICHwall_dz/2.0 ) );

    G4Box *RICHexit_cutout = new G4Box( "RICHexit_cutout", RICH_exitwindow_dx, RICH_exitwindow_dy, 10.0*cm );
    G4SubtractionSolid *RICHbox_wall_entryexitcut = new G4SubtractionSolid( "RICHbox_wall_entryexitcut", RICHbox_wall_entrycut, RICHexit_cutout, 0, G4ThreeVector( x0_RICHexit, y0_RICHexit, RICHwall_dz/2.0 ) );

    G4LogicalVolume *RICH_container_walls = new G4LogicalVolume( RICHbox_wall_entryexitcut, GetMaterial("RICHAluminum"), "RICH_container_walls" );

    new G4PVPlacement( 0, G4ThreeVector(0,0,0), RICH_container_walls, "RICH_container_walls_placement", RICHbox_log, false, 0 );

    //We also need to define the Tedlar spacers: Let us define horizontal and vertical spacers. Our convention will be that the 
    //vertical spacer fills the corner region:
    G4double horizontal_spacer_dx = aero_dx, horizontal_spacer_dy = tilegap/2.0, horizontal_spacer_dz=Lz_aero/2.0;
    G4double vertical_spacer_dx = tilegap/2.0, vertical_spacer_dy = Ly_aero/2.0, vertical_spacer_dz=Lz_aero/2.0;

    G4Box *Horizontal_spacer = new G4Box( "Horizontal_spacer", horizontal_spacer_dx, horizontal_spacer_dy, horizontal_spacer_dz );
    G4Box *Vertical_spacer = new G4Box( "Vertical_spacer", vertical_spacer_dx, vertical_spacer_dy, vertical_spacer_dz );

    G4LogicalVolume *Horizontal_spacer_log = new G4LogicalVolume( Horizontal_spacer, GetMaterial("Tedlar"), "Horizontal_spacer_log" );
    G4LogicalVolume *Vertical_spacer_log = new G4LogicalVolume( Vertical_spacer, GetMaterial("Tedlar"), "Vertical_spacer_log" );

    G4int icopy = 1;
    G4int icopy_vertical_spacer = 1;
    G4int icopy_horizontal_spacer = 1;

    G4String tilename;

    //Next: position aerogel tiles: 
    for(G4int ix=0; ix<nx_aero; ix++ ){

	G4double xtemp = -RICHbox_dx/2.0 + aero_xoffset + (ix+0.5)*2.0*aero_dx + ix*tilegap;

	for(G4int iy=0; iy<ny_aero; iy++ ){

	    G4double ytemp = -Ly_aero/2.0 + (iy+0.5)*2.0*aero_dy + iy*tilegap;

	    for(G4int iz=0; iz<nz_aero; iz++ ){
		//compute center coordinates of aerogel tiles and position in RICHbox:

		G4double ztemp = -RICHbox_dz/2.0 + RICH_entrywindow_thick + (iz+0.5)*2.0*aero_dz;

		tilename = "Aerogel_tile_pv_";

		char ccopy[20];

		sprintf( ccopy, "%d", icopy );

		tilename += ccopy;

		new G4PVPlacement( 0, 
			G4ThreeVector( xtemp, ytemp, ztemp ), 
			Aerogel_tile_log, 
			tilename,
			RICHbox_log,
			false,
			icopy++ );

	    }
	    if( iy > 0 && iy+1 < ny_aero ){
		G4double xspacer = xtemp;
		G4double yspacer = ytemp + aero_dy + tilegap/2.0;
		G4double zspacer = -RICHbox_dz/2.0 + RICH_entrywindow_thick + Lz_aero/2.0;

		new G4PVPlacement( 0, 
			G4ThreeVector( xspacer, yspacer, zspacer ), 
			Horizontal_spacer_log, 
			"Horizontal_spacer_pv",
			RICHbox_log, 
			false, 
			icopy_horizontal_spacer++ );

	    }
	} 
	//position vertical tedlar spacers:
	if( ix>0 && ix+1 < nx_aero ){
	    //vertical spacer position is equal to tile position + half tile width + half gap width:
	    G4double xspacer = xtemp + aero_dx + tilegap/2.0;
	    G4double yspacer = 0.0;
	    G4double zspacer = -RICHbox_dz/2.0 + RICH_entrywindow_thick + Lz_aero/2.0;

	    new G4PVPlacement( 0, 
		    G4ThreeVector( xspacer, yspacer, zspacer ),
		    Vertical_spacer_log, 
		    "Vertical_spacer_pv",
		    RICHbox_log, 
		    false,
		    icopy_vertical_spacer++ );

	}
    }

    //Next, let's try to define the mirror. For this, we can probably use a "spherical shell section" without resorting to 
    //solid operations. Alternatively we could use polycone.

    G4double MirrorCenter_x = -RICHbox_dx/2.0 + 136.403*cm;
    G4double MirrorCenter_y = 0.0*cm;
    G4double MirrorCenter_z = -RICHbox_dz/2.0 - 103.112*cm;

    G4double MirrorShell_thick = 0.1932*cm; //This corresponds to 1% X0 of graphite.

    G4double MirrorRadius = 220.0*cm;
    //G4double MirrorRadius_XZproject = 180.190*cm;

    //This is relative to RICHbox
    G4double Mirror_xmin = -RICHbox_dx/2.0 + 5.540*cm;
    G4double Mirror_xmax = 79.505*cm + Mirror_xmin;

    G4double Mirror_zmin = MirrorCenter_z + 123.890*cm;
    G4double Mirror_zmax = RICHbox_dz/2.0;

    G4double Mirror_ymin = -126.2*cm;
    G4double Mirror_ymax = 126.2*cm;

    //Let's make a spherical shell and a "cut box"

    //The spherical shell should go from 0 to 90 degrees in theta and +/- 90 degrees in phi, and then we form the intersection 
    // with the cut box defining the minimum and maximum z planes:
    //Include the entire forward hemisphere for simplicity:
    G4Sphere *RICH_mirror_shell = new G4Sphere( "RICH_mirror_shell", MirrorRadius, MirrorRadius + MirrorShell_thick, 0.0, twopi, 0.0, halfpi );
    G4Box *RICH_mirror_cutbox = new G4Box( "RICH_mirror_cutbox", (Mirror_xmax-Mirror_xmin)/2.0, (Mirror_ymax-Mirror_ymin)/2.0, (Mirror_zmax-Mirror_zmin)/2.0 );

    //Now, we want to make the intersection of the two solids, so we need to express the coordinates of the center of the box in the mirror
    // coordinate system. 
    G4ThreeVector MirrorCenterCoords( MirrorCenter_x, MirrorCenter_y, MirrorCenter_z );
    G4ThreeVector BoxCenterCoords( 0.5*(Mirror_xmin+Mirror_xmax), 0.0, 0.5*(Mirror_zmin+Mirror_zmax) );
    G4ThreeVector RelativeCoords = BoxCenterCoords - MirrorCenterCoords;

    G4IntersectionSolid *Mirror_solid = new G4IntersectionSolid( "Mirror_solid", RICH_mirror_shell, RICH_mirror_cutbox, 0, RelativeCoords );
    G4LogicalVolume *Mirror_log = new G4LogicalVolume( Mirror_solid, GetMaterial("MirrorComposite"), "Mirror_log" );

    new G4PVPlacement( 0, 
	    MirrorCenterCoords, 
	    Mirror_log, 
	    "Mirror_pv",
	    RICHbox_log, 
	    false,
	    0 );

    new G4LogicalSkinSurface( "Mirrskin", Mirror_log, GetOpticalSurface("Mirrsurf") );

    //What is left? We've done the mirror, the aerogel, the gas box. All that remains is the PMTs and the structure of the containment vessel. Let's start with the PMTs: 

    ////////////////////////////////////////////////////////////////////////
    //                         !!!PMTS!!!                                 //
    ////////////////////////////////////////////////////////////////////////


    //Define the PMT windows as 1 mm-thick discs of "UVglass":
    G4Tubs *PMTwindow = new G4Tubs( "PMTwindow", 0.0*cm, (1.86/2.0)*cm, 0.05*cm, 0.0, twopi ); 
    //Define the PMT photocathode as a thin disc of 15 mm
    G4Tubs *PMTcathode = new G4Tubs( "PMTcathode", 0.0*cm, (1.50/2.0)*cm, 0.025*cm, 0.0, twopi );
    G4Tubs *PMTtube    = new G4Tubs( "PMTtube", (1.66/2.0)*cm, (1.86/2.0)*cm, ((8.7+0.3)/2.0)*cm, 0.0, twopi );
    G4Tubs *PMTendcap  = new G4Tubs( "PMTendcap", 0.0*cm, (1.66/2.0)*cm, 0.15*cm, 0.0, twopi );

    //"Quartz window" is a different, sealed window that separates the PMT from the C4F10 environment.
    G4Tubs *PMTQuartzWindow = new G4Tubs( "PMTQuartzWindow", 0.0*cm, (2.13/2.0)*cm, 0.15*cm, 0.0, twopi );
    //CollectionCone is a light-collecting cone that increases the effective collection efficiency:
    G4Cons *CollectionCone = new G4Cons( "CollectionCone", 0.75*cm, (2.13/2.0)*cm, (2.13/2.0)*cm, (2.13/2.0)*cm, 0.75*cm, 0.0, twopi );

    G4double PMT_total_length = 10.9*cm;
    //    G4double PMT_max_radius = 1.065*cm;

    G4LogicalVolume *PMTwindow_log  = new G4LogicalVolume( PMTwindow, GetMaterial("UVglass"), "PMTwindow_log" );
    G4LogicalVolume *PMTcathode_log = new G4LogicalVolume( PMTcathode, GetMaterial("Photocathode_material"), "PMTcathode_log" );

    //PMTcathode_log is the sensitive detector for the RICH:

    //  G4SDManager *fSDman = G4SDManager::GetSDMpointer();
    G4SDManager *sdman = fDetCon->fSDman;

    G4String RICHSDname = "G4SBS/RICH";
    G4String RICHcollname = "RICHcoll";
    G4SBSRICHSD *RICHSD = NULL;

    if( !( RICHSD = (G4SBSRICHSD*) sdman->FindSensitiveDetector(RICHSDname) ) ){
	G4cout << "Adding RICH sensitive detector to SDman..." << G4endl;
	RICHSD = new G4SBSRICHSD( RICHSDname, RICHcollname );
	sdman->AddNewDetector( RICHSD );
	fDetCon->SDlist[RICHSDname] = RICHSD;
    }

    PMTcathode_log->SetSensitiveDetector( RICHSD ); //This assigns the sensitive detector type "RICHSD" to the logical volume PMTcathode!
    //We make this a hollow cylinder with length and radius approximately equal to that of the PMT housing, made of steel 
    //to approximate the material shielding the PMT.
    G4LogicalVolume *PMTtube_log    = new G4LogicalVolume( PMTtube, GetMaterial("Steel"), "PMTtube_log" ); 
    G4LogicalVolume *PMTendcap_log  = new G4LogicalVolume( PMTendcap, GetMaterial("Steel"), "PMTendcap_log" );
    G4LogicalVolume *PMTquartzwindow_log = new G4LogicalVolume( PMTQuartzWindow, GetMaterial("QuartzWindow"), "PMTQuartzWindow_log" );
    G4LogicalVolume *CollectionCone_log = new G4LogicalVolume( CollectionCone, GetMaterial("Steel"), "CollectionCone_log" );
    //Define a logical skin surface for the collection cone and assign it the same reflectivity as the mirror:
    new G4LogicalSkinSurface( "Coneskin", CollectionCone_log, GetOpticalSurface("Mirrsurf") );

    //Within the RICHbox, each PMT assembly unit is rotated so that its z-axis makes an angle of 50 degrees with respect to the 
    //local z axis of the RICHbox. Therefore, we rotate by an angle of 
    G4double PMT_rotation_angle = 50.0*degree;
    G4RotationMatrix *rot_PMT = new G4RotationMatrix;
    rot_PMT->rotateY( PMT_rotation_angle );

    G4int icopy_PMT_assembly = 0;

    G4double xfp = 119.350*cm - RICHbox_dx/2.0;
    G4double yfp = 0.0;
    G4double zfp = 42.521*cm - RICHbox_dz/2.0;

    G4ThreeVector focalpoint_position( xfp, yfp, zfp );

    G4ThreeVector PMT_zaxis( -sin(PMT_rotation_angle), 0.0, cos(PMT_rotation_angle) );
    G4ThreeVector PMT_yaxis( 0, 1, 0 );
    G4ThreeVector PMT_xaxis( (PMT_yaxis.cross( PMT_zaxis ) ).unit() );

    G4double ymin_PMT = -72.5376*cm, ymax_PMT = 72.5376*cm;
    G4double xmin_PMT[2] = { -29.083*cm, -30.24632*cm };
    G4double xmax_PMT[2] = { 29.083*cm, 30.24632*cm };
    G4int nrows_PMT[2] = {26, 27};

    for( G4int icol=0; icol<=72; icol++){
	G4int evenoddcol = icol%2;
	for( G4int irow=0; irow<nrows_PMT[evenoddcol]; irow++ ){
	    G4double xtemp = xmin_PMT[evenoddcol] + irow * ( xmax_PMT[evenoddcol] - xmin_PMT[evenoddcol] )/( G4double(nrows_PMT[evenoddcol]-1) );
	    G4double ytemp = ymin_PMT + icol*(ymax_PMT-ymin_PMT)/( 72.0 );

	    G4ThreeVector PMT_position = focalpoint_position - PMT_zaxis * PMT_total_length/2.0 + xtemp * PMT_xaxis + ytemp * PMT_yaxis;

	    //Place PMT components inside RICHbox.
	    G4ThreeVector Pos_temp;
	    //Steel tube (mainly for visualization and shielding
	    G4double ztube = -PMT_total_length/2.0 + 4.5*cm;
	    Pos_temp = PMT_position + ztube * PMT_zaxis;
	    new G4PVPlacement( rot_PMT, Pos_temp, PMTtube_log, "PMTtube_pv", RICHbox_log, false, icopy_PMT_assembly );
	    //Endcap of steel tube (keep optical photons originating from behind PMTs from hitting the cathode):
	    G4double zendcap = -PMT_total_length/2.0 + 0.15*cm;
	    Pos_temp = PMT_position + zendcap * PMT_zaxis;
	    new G4PVPlacement( rot_PMT, Pos_temp, PMTendcap_log, "PMTendcap_pv", RICHbox_log, false, icopy_PMT_assembly );
	    //Photocathode (this is the sensitive part!!):
	    G4double zcathode = ztube + 4.475*cm;
	    Pos_temp = PMT_position + zcathode * PMT_zaxis;
	    new G4PVPlacement( rot_PMT, Pos_temp, PMTcathode_log, "PMTcathode_pv", RICHbox_log, false, icopy_PMT_assembly );
	    //UV-glass PMT window:
	    G4double zwindow = zcathode + 0.075*cm;
	    Pos_temp = PMT_position + zwindow * PMT_zaxis;
	    new G4PVPlacement( rot_PMT, Pos_temp, PMTwindow_log, "PMTwindow_pv", RICHbox_log, false, icopy_PMT_assembly );
	    //Quartz window between PMT and gas:
	    G4double zquartz = zwindow + 0.2*cm;
	    Pos_temp = PMT_position + zquartz * PMT_zaxis;
	    new G4PVPlacement( rot_PMT, Pos_temp, PMTquartzwindow_log, "PMTquartzwindow_pv", RICHbox_log, false, icopy_PMT_assembly );
	    //Light collection cone:
	    G4double zcone = zquartz + 0.9*cm;
	    Pos_temp = PMT_position + zcone * PMT_zaxis;
	    new G4PVPlacement( rot_PMT, Pos_temp, CollectionCone_log, "CollectionCone_pv", RICHbox_log, false, icopy_PMT_assembly );

	    // G4VPhysicalVolume *PMT_placement = new G4PVPlacement( rot_PMT, 
	    // 							   PMT_position, 
	    // 							   PMT_assembly, 
	    // 							   "PMT_placement", 
	    // 							   RICHbox_log, 
	    // 							   false, 
	    // 							   icopy_PMT_assembly++ );
	    icopy_PMT_assembly++;


	}
    }


    ////////////////////////////////////////////////////////////////////////
    //                         !!!END OF PMTS!!!                                 //
    ////////////////////////////////////////////////////////////////////////

    //At the end we have to define the translation and rotation to apply for correct positioning of RICHbox:
    //Rotation is easy, same as HCAL, we rotate about the Y axis by -f48D48ang:

    G4RotationMatrix *rot_RICH = new G4RotationMatrix;
    rot_RICH->rotateY( -f48D48ang );

    //We want the center of the RICH entry window to be located at a distance equal to fRICHdist along the line at angle f48D48ang from the origin. For this condition to be satisfied, the center of the RICH box must be offset from this line:
    G4ThreeVector RICHcoord_global( fRICHdist*sin( f48D48ang ), 0.0, fRICHdist*cos( f48D48ang ) );

    G4ThreeVector RICH_zaxis( RICHcoord_global.unit() );
    G4ThreeVector RICH_yaxis( 0.0, 1.0, 0.0 );
    G4ThreeVector RICH_xaxis( (RICH_yaxis.cross( RICH_zaxis )).unit() );

    //RICH center coordinates

    G4ThreeVector RICH_centercoord_global = RICHcoord_global - x0_entrywindow * RICH_xaxis - y0_entrywindow * RICH_yaxis - z0_entrywindow * RICH_zaxis;

    new G4PVPlacement( rot_RICH, 
	    RICH_centercoord_global,
	    RICHbox_log,
	    "RICHbox_pv",
	    motherlog, 
	    false, 
	    0 );


    G4VisAttributes *RICHbox_vis = new G4VisAttributes( G4Colour( 1.0, 1.0, 1.0 ) );
    RICHbox_vis->SetForceWireframe(true);

    RICHbox_log->SetVisAttributes( RICHbox_vis ); 

    //Set color and transparency for RICH windows (Aluminum)
    G4VisAttributes *RICHwindow_visatt = new G4VisAttributes( G4Colour( 0.75,0.75,0.75) );

    //RICH entry and exit windows are not inherently interesting, so we force them to wireframe:
    RICHwindow_visatt->SetForceWireframe(true);

    RICH_exitwindow->SetVisAttributes( RICHwindow_visatt );
    RICH_entrywindow_log->SetVisAttributes( RICHwindow_visatt );

    //Set aerogel exit window to a magenta color (equal parts red and blue) and also wireframe:
    G4VisAttributes *Lucitewindow_visatt = new G4VisAttributes( G4Colour( 1.0,0.0,1.0 ) );
    Lucitewindow_visatt->SetForceWireframe(true); 

    Aero_exitwindow->SetVisAttributes( Lucitewindow_visatt );

    G4VisAttributes *aero_tile_visatt = new G4VisAttributes( G4Colour( 0.0, 0.8, 0.8 ) );
    Aerogel_tile_log->SetVisAttributes( aero_tile_visatt );

    G4VisAttributes *tedlar_vis = new G4VisAttributes( G4Colour(0.3,0.3,0.3) );
    Vertical_spacer_log->SetVisAttributes( tedlar_vis );
    Horizontal_spacer_log->SetVisAttributes( tedlar_vis );

    G4VisAttributes *mirror_vis = new G4VisAttributes( G4Colour(0.0, 0.5, 1.0) );
    Mirror_log->SetVisAttributes( mirror_vis );

    //PMT_assembly->SetVisAttributes( G4VisAttributes::GetInvisible() );

    //  G4VisAttributes for PMT assemblies:

    G4VisAttributes *PMTtube_vis = new G4VisAttributes( G4Colour( 0.4, 0.4, 0.4 ) );
    PMTtube_vis->SetForceLineSegmentsPerCircle( 12 );
    PMTtube_log->SetVisAttributes( PMTtube_vis );
    PMTendcap_log->SetVisAttributes( PMTtube_vis );

    G4VisAttributes *PMTwindow_vis = new G4VisAttributes( G4Colour::Cyan() );
    PMTwindow_vis->SetForceLineSegmentsPerCircle( 12 );
    PMTwindow_vis->SetForceWireframe( true );
    PMTwindow_log->SetVisAttributes( PMTwindow_vis );

    //  PMTcathode_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
    G4VisAttributes *PMTcathode_vis = new G4VisAttributes( G4Colour::Blue() );
    PMTcathode_vis->SetForceLineSegmentsPerCircle( 12 );
    PMTcathode_log->SetVisAttributes( PMTcathode_vis );

    G4VisAttributes *PMTquartzwindow_vis = new G4VisAttributes( G4Colour::Green() );
    PMTquartzwindow_vis->SetForceLineSegmentsPerCircle( 12 );
    PMTquartzwindow_vis->SetForceWireframe( true );
    PMTquartzwindow_log->SetVisAttributes( PMTquartzwindow_vis );

    G4VisAttributes *CollectionCone_vis = new G4VisAttributes( G4Colour::Red() );
    CollectionCone_vis->SetForceLineSegmentsPerCircle( 12 );
    CollectionCone_log->SetVisAttributes( CollectionCone_vis );

    G4VisAttributes *RICHwalls_vis = new G4VisAttributes( G4Colour::Gray() );
    RICHwalls_vis->SetForceWireframe( true );
    RICH_container_walls->SetVisAttributes( RICHwalls_vis );

}



void G4SBSHArmBuilder::MakeFPP( G4LogicalVolume *Mother, G4RotationMatrix *rot, G4ThreeVector pos ){
    //FPP consists of GEM tracker interspersed w/CH2 analyzers:

  int i, j;
  //int ngem = 0;
  
  vector<double> gemz, gemw, gemh;
  
  int ntracker = 3; //FT, FPP1, FPP2
  int ngem[3] = {6,4,4};
  
  //int TrackerID = 2;
  G4SBSTrackerBuilder trackerbuilder(fDetCon);

  for( i = 0; i<ntracker; i++){
    gemz.resize( ngem[i] );
    gemw.resize( ngem[i] );
    gemh.resize( ngem[i] );
    for( j = 0; j < ngem[i]; j++ ){
      if( i == 0 ){
	gemz[j] = ((double) j)*9*cm;
	gemw[j] = 40.0*cm;
	gemh[j] = 150.0*cm;
      } else if( i == 1 ){
	gemz[j] = ((double) j)*16*cm + 1.2*m;
	//	  gemz[i] = pairspac*((i-6)/2) + (i%2)*gemdsep + 1.2*m;
	gemw[j] = 50.0*cm;
	gemh[j] = 200.0*cm;
      } else {
	gemz[j] = ((double) j)*16.*cm + 2.316*m;
	//	  gemz[i] = pairspac*((i-10)/2) + (i%2)*gemdsep + 2.316*m;
	gemw[j] = 50.0*cm;
	gemh[j] = 200.0*cm;
      }
    }
    (fDetCon->TrackerArm)[fDetCon->TrackerIDnumber] = kHarm; //1 is H arm.  
    trackerbuilder.BuildComponent( Mother, rot, pos, ngem[i], gemz, gemw, gemh, (fDetCon->TrackerIDnumber)++ );
  }
  //CH2 analyzers:
  
  double anaheight = 200.0*cm;
  double anawidth  = 44.0*2.54*cm;
  double anadepth  = 22.0*2.54*cm;
  
  G4Box *anabox = new G4Box("anabox", anawidth/2.0, anaheight/2.0, anadepth/2.0 );
  G4LogicalVolume* analog = new G4LogicalVolume(anabox, GetMaterial("CH2"), "analog");
  
  G4ThreeVector Ana1_pos = pos + G4ThreeVector( 0.0, 0.0, 58.53*cm + anadepth/2.0 );
  G4ThreeVector Ana2_pos = pos + G4ThreeVector( 0.0, 0.0, 170.3*cm + anadepth/2.0 );
  
  new G4PVPlacement(0, Ana1_pos, analog,
		    "anaphys1", Mother, false, 0, false);
  new G4PVPlacement(0, Ana2_pos, analog,
		    "anaphys1", Mother, false, 0, false);

  G4VisAttributes *CH2anavisatt = new G4VisAttributes( G4Colour(0.0, 0.0, 1.0) );
  CH2anavisatt->SetForceWireframe(true);

  analog->SetVisAttributes( CH2anavisatt );
  
}
