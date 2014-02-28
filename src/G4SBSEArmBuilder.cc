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

#include "G4SBSTrackerBuilder.hh"
#include "G4SBSCalSD.hh"

G4SBSEArmBuilder::G4SBSEArmBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
    fBBang  = 40.0*deg;
    fBBdist = 1.5*m;

    fBBCaldist = 0.8*m;
    fCerDepth = 60.0*cm;
    fCerDist  =  7.0*cm;

    fGEMDist  = 70.0*cm;
    fGEMOption = 1;

    assert(fDetCon);
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

    (fDetCon->TrackerArm)[fDetCon->TrackerIDnumber] = kEarm;
    
    trackerbuilder.BuildComponent(bbdetLog, rot_identity, G4ThreeVector( 0.0, 0.0, detoffset ), ngem, gemz, gemw, gemh, (fDetCon->TrackerIDnumber)++ );
    //----- Note: Lines of code that are common to the construction of all individual GEM planes/modules were moved to MakeTracker() -----// 
    //----- All we do here in MakeBigBite() is define the number of planes, their z positions, and their transverse dimensions ------//

    // BigBite Preshower 
    // AJP: Why make a preshower box if it's just going to be full of air and not a sensitive detector?

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

    G4String BBCalSDname = "G4SBS/BBCal";
    G4String BBCalcolname = "BBCalcol";
    G4SBSCalSD* BBCalSD;

    if( !(BBCalSD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(BBCalSDname)) ){
	BBCalSD = new G4SBSCalSD( BBCalSDname, BBCalcolname );
	fDetCon->fSDman->AddNewDetector(BBCalSD);
	fDetCon->SDlist[BBCalSDname] = BBCalSD;
    }

    bbcallog->SetSensitiveDetector(BBCalSD);
    bbcallog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );


    //--------- BigBite Cerenkov ------------------------------

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

    /*
     *  FIXME  - NO CERENKOV */
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

    //--------- Visualization attributes -------------------------------
    bbdetLog->SetVisAttributes(G4VisAttributes::Invisible);
    bbfieldLog->SetVisAttributes(G4VisAttributes::Invisible);
    bbmotherLog->SetVisAttributes(G4VisAttributes::Invisible);

    G4VisAttributes * yokeVisAtt
	= new G4VisAttributes(G4Colour(0.0,0.0,1.0));
    //  yokeVisAtt->SetForceWireframe(true);
    bbyokewgapLog->SetVisAttributes(yokeVisAtt);

    G4VisAttributes * alVisAtt
	= new G4VisAttributes(G4Colour(0.1,0.1,0.1));
    cer_winlog_in->SetVisAttributes(alVisAtt);
    cer_winlog_out->SetVisAttributes(alVisAtt);

    G4VisAttributes * gasVisAtt
	= new G4VisAttributes(G4Colour(0.6,0.6,1.0));
    gasVisAtt->SetForceWireframe(true);
    cer_gaslog->SetVisAttributes(gasVisAtt);

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

    G4RotationMatrix *bbrm = new G4RotationMatrix;
    bbrm->rotateY(fBBang);

    // Ecal will act as BBcal detector

    double bigcalheight = (24*4.5+32*4.0)*cm;
    double bigcalwidth  = 44.10*2.54*cm;
    double bigcaldepth  = 15.75*2.54*cm;
    double bbr = fBBdist+bigcaldepth/2.0;

    double CH2depth = 15.0*cm;
    double CHdepth  = 6.0*cm;

    G4Box *bigcalbox = new G4Box("bigcalbox", bigcalwidth/2.0, bigcalheight/2.0, bigcaldepth/2.0 );
    G4LogicalVolume* bigcallog = new G4LogicalVolume(bigcalbox, GetMaterial("Lead"), "bigcallog");

    G4Box *CH2box = new G4Box("ch2box", bigcalwidth/2.0, bigcalheight/2.0, CH2depth/2.0 );
    G4LogicalVolume* ch2boxlog = new G4LogicalVolume(CH2box, GetMaterial("CH2"), "ch2log");
    G4Box *CHbox = new G4Box("chbox", bigcalwidth/2.0, bigcalheight/2.0, CHdepth/2.0 );
    G4LogicalVolume* chboxlog = new G4LogicalVolume(CHbox, GetMaterial("CH"), "chlog");

    double ch2r = bbr - bigcaldepth/2.0 - CHdepth - CH2depth/2.0;
    double chr = bbr - bigcaldepth/2.0 - CHdepth/2.0;

    new G4PVPlacement(bbrm, G4ThreeVector(bbr*sin(-fBBang), 0.0, bbr*cos(-fBBang) ), bigcallog,
	    "bigcalphys", worldlog, false, 0, false);
    new G4PVPlacement(bbrm, G4ThreeVector(ch2r*sin(-fBBang), 0.0, ch2r*cos(-fBBang) ), ch2boxlog,
	    "ch2boxphys", worldlog, false, 0, false);
    new G4PVPlacement(bbrm, G4ThreeVector(chr*sin(-fBBang), 0.0, chr*cos(-fBBang) ), chboxlog,
	    "chboxphys", worldlog, false, 0, false);

    G4String BBCALSDname = "G4SBS/BBCal";
    G4String BBCALcolname = "BBCalcol";
    G4SBSCalSD* BBCalSD;

    if( !(BBCalSD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(BBCALSDname)) ){
	BBCalSD = new G4SBSCalSD( BBCALSDname, BBCALcolname );
	fDetCon->fSDman->AddNewDetector(BBCalSD);
	fDetCon->SDlist[BBCALSDname] = BBCalSD;
    }

    fDetCon->fSDman->AddNewDetector(BBCalSD);
    bigcallog->SetSensitiveDetector(BBCalSD);
    bigcallog->SetUserLimits(  new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );

    G4VisAttributes * bcVisAtt
	= new G4VisAttributes(G4Colour(0.0,1.0,1.0));
    bigcallog->SetVisAttributes(bcVisAtt);

    G4VisAttributes * chVisAtt
	= new G4VisAttributes(G4Colour(1.0,0.6,0.0));
    chboxlog->SetVisAttributes(chVisAtt);

    G4VisAttributes * ch2VisAtt
	= new G4VisAttributes(G4Colour(1.0,0.8,0.0));
    ch2boxlog->SetVisAttributes(ch2VisAtt);

}


