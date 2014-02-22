#include "G4SBSTargetBuilder.hh" 

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
#include "G4Sphere.hh"
#include "G4TwoVector.hh"
#include "G4GenericTrap.hh"
#include "G4Polycone.hh"

G4SBSTargetBuilder::G4SBSTargetBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
    fTargLen = 60.0*cm;
    fTargType = kH2;
    fTargDen = 10.5*atmosphere/(300*kelvin*k_Boltzmann);
}

G4SBSTargetBuilder::~G4SBSTargetBuilder(){;}

void G4SBSTargetBuilder::BuildComponent(G4LogicalVolume *worldlog){

    //Material definition was moved to ConstructMaterials();
    //--------- Glass target cell -------------------------------

    double wallthick = 1.61*mm;
    double capthick  = 0.126*mm;
    double cellradius    = 0.75*2.54*cm/2.0;

    G4Tubs *targ_tube = new G4Tubs("targ_tube", cellradius-wallthick, cellradius, fTargLen/2.0, 0.*deg, 360.*deg );
    G4Tubs *targ_cap = new G4Tubs("targ_cap", 0.0, cellradius, capthick/2.0, 0.*deg, 360.*deg );

    G4LogicalVolume* targ_tube_log = new G4LogicalVolume(targ_tube, GetMaterial("GE180"),"targ_tube_log");
    G4LogicalVolume* targ_cap_log = new G4LogicalVolume(targ_cap, GetMaterial("GE180"),"targ_cap_log");

    /* FIXME
     * */
    if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
	new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), targ_tube_log,
		"targ_tube_phys", worldlog, false, 0);

	new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, fTargLen/2.0+capthick/2.0), targ_cap_log,
		"targ_cap_phys1", worldlog, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -fTargLen/2.0-capthick/2.0), targ_cap_log,
		"targ_cap_phys2", worldlog, false, 0);
    }


    /**/


    // gas
    G4Tubs *gas_tube = new G4Tubs("gas_tube", 0.0, cellradius-wallthick,fTargLen/2.0, 0.*deg, 360.*deg );
    G4LogicalVolume* gas_tube_log;


    if( fTargType == kH2 || fTargType == kNeutTarg ){
	gas_tube_log = new G4LogicalVolume(gas_tube, GetMaterial("refH2"), "gas_tube_log");
    }
    if( fTargType == k3He ){
	gas_tube_log = new G4LogicalVolume(gas_tube, GetMaterial("pol3He"), "gas_tube_log");
    }

    /*
     * FIXME*/
    if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
	new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), gas_tube_log,
		"gas_tube_phys", worldlog, false, 0);
    }
    /**/

    //--------- Cryo target cell -------------------------------

    wallthick   = 635*um;
    double upcapthick    = 100*um;
    double downcapthick  = 125*um;
    double targconeang = 15.*deg;
    double cellupradius = 2.0*cm;

    double celldownradius = fTargLen*sin(targconeang);
    double cellconelen = fTargLen*cos(targconeang);
    double cellconeang = atan((celldownradius-cellupradius)/cellconelen);

    // Aluminum shell sphere
    G4Sphere *shellsph = new G4Sphere("shellsph", 0, fTargLen, 0, 360.0*deg, 0, targconeang);
    // Aluminum shell cone
    G4Cons *shellcon = new G4Cons("shellcon", 0.0, cellupradius, 0.0, celldownradius,  cellconelen/2, 0.0, 360.0*deg);
    // Union 
    G4UnionSolid *cryoshell = new G4UnionSolid("cryoshell", shellcon, shellsph, 0, G4ThreeVector(0,0,-cellconelen/2));

    double cryoupradius =  cellupradius - wallthick/cos(cellconeang);
    double cryoconelen = cellconelen - upcapthick - downcapthick; 

    double cryodownradius = cryoupradius + cryoconelen*tan(cellconeang);

    double cryooffset = cryoconelen/2+upcapthick- cellconelen/2;

    double cryoang = 14.95*deg;

    // Cryo sphere
    G4Sphere *cryosph = new G4Sphere("cryosph", 0, fTargLen-downcapthick, 0*deg, 360.0*deg, 0, cryoang);
    // Cryo cone
    G4Cons *cryocon = new G4Cons("shellcon", 0.0, cryoupradius, 0.0, cryodownradius,  cryoconelen/2,  0.0*deg, 360.0*deg);
    // Union 
    G4UnionSolid *cryovol1 = new G4UnionSolid("cryovol1", cryocon, cryosph, 0, G4ThreeVector(0.0,0.0,-cryoconelen/2.-upcapthick ));

    double trimboxsize = 50.0*cm;
    G4Box *cryotrimbox = new G4Box("cryotrimbox", trimboxsize, trimboxsize, trimboxsize);

    G4SubtractionSolid *cryovol = new G4SubtractionSolid("cryovol", cryovol1, cryotrimbox, 0, G4ThreeVector(0.0,0.0,-trimboxsize-cryoconelen/2 ));

    /*
       targ_tube = new G4Tubs("targ_tube", cellradius-wallthick, cellradius, fTargLen/2.0, 0.*deg, 360.*deg );
       G4Tubs *targ_ucap = new G4Tubs("targ_ucap", 0.0, cellradius, upcapthick/2.0, 0.*deg, 360.*deg );
       G4Tubs *targ_dcap = new G4Tubs("targ_dcap", 0.0, cellradius, downcapthick/2.0, 0.*deg, 360.*deg );
       targ_tube_log = new G4LogicalVolume(targ_tube, Aluminum,"targ_tube_log");
       */

    targ_tube_log = new G4LogicalVolume(cryoshell, GetMaterial("Aluminum"),"targ_tube_log");

    //////////////////////////////////////////////////////////////////

    G4double entpipe_rin = 31.75*mm;
    G4double extpipe_rin = 41.28*mm;

    G4double extpipestart = 2.06*m;
    // 2.06m is where the main exit pipe starts
    G4double extpipe_len;

    double sheight = 1.2*m;

    double swallthick   = 0.38*mm;
    double swallrad     = 1.143*m/2;
    double swallrad_in  = 1.041*m/2;

    double hcal_ang_min = -55*deg;
    double hcal_ang_max = -7*deg;
    double hcal_win_h = 0.4*m;

    double bb_ang_min = 18*deg;
    double bb_ang_max = 80*deg;
    double bb_win_h = 0.5*m;

    if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
	// Gas target
	extpipe_len = extpipestart - 1.0*m;
    } else {
	// Cryotarget
	extpipe_len = extpipestart -  swallrad;
    }

    G4Tubs *swall = new G4Tubs("scham_wall", swallrad_in, swallrad, sheight/2, 0.*deg, 360.*deg );

    // Cut out for windows

    G4Tubs *swall_hcalcut = new G4Tubs("scham_wall_hcalcut", swallrad_in-2*cm, swallrad+2*cm, hcal_win_h, hcal_ang_min, hcal_ang_max-hcal_ang_min);
    G4Tubs *swall_bbcut = new G4Tubs("scham_wall_bbcut", swallrad_in-2*cm, swallrad+2*cm, bb_win_h, bb_ang_min, bb_ang_max-bb_ang_min);

    G4SubtractionSolid *swallcut = new G4SubtractionSolid("swallcut1", swall, swall_hcalcut);
    swallcut = new G4SubtractionSolid("swallcut2", swallcut, swall_bbcut);

    G4Tubs *swall_hcalwin = new G4Tubs("scham_wall_hcalwin", swallrad_in, swallrad_in+swallthick, hcal_win_h, hcal_ang_min, hcal_ang_max-hcal_ang_min);
    G4Tubs *swall_bbwin = new G4Tubs("scham_wall_bbwin", swallrad_in, swallrad_in+swallthick, bb_win_h, bb_ang_min, bb_ang_max-bb_ang_min);

    G4Tubs *exttube = new G4Tubs("exitpipetube", extpipe_rin, extpipe_rin+0.120*cm, extpipe_len/2, 0.*deg, 360.*deg );
    G4Tubs *extvactube = new G4Tubs("exitpipetube_vac", 0.0, extpipe_rin, extpipe_len, 0.*deg, 360.*deg );
    G4LogicalVolume *extpipe_log = new G4LogicalVolume(exttube, GetMaterial("Aluminum"),"extpipe_log");
    G4LogicalVolume *extvac_log = new G4LogicalVolume(extvactube, GetMaterial("Vacuum"),"extvac_log");

    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, extpipestart-extpipe_len/2), extpipe_log, "extpipe_phys", worldlog, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, extpipestart-extpipe_len/2), extvac_log, "extvacpipe_phys", worldlog, false, 0);

    if( fTargType == kH2 || fTargType == k3He || fTargType == kNeutTarg ){
	// Add in exit Al window

	double extwin_thick = 5.0e-4*cm;

	G4Tubs *extwin = new G4Tubs("ent_win", 0.0, extpipe_rin, extwin_thick/2, 0.*deg, 360.*deg );
	G4LogicalVolume *ext_winlog = new G4LogicalVolume(extwin, GetMaterial("Aluminum"), "entwin_log", 0, 0, 0);
	new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -extpipe_len/2 + extwin_thick/2), ext_winlog, "extwin_phys", extvac_log,false,0);

	ext_winlog->SetVisAttributes(new G4VisAttributes(G4Colour(0.6,0.6,0.6)));
    }


    //  Place exit pipe tube

    G4Tubs *swall_enthole = new G4Tubs("scham_wall_enthole", 0.0, entpipe_rin, 20.0*cm, 0.*deg, 360.*deg );
    G4Tubs *swall_exthole = new G4Tubs("scham_wall_exthole", 0.0, extpipe_rin, 20.0*cm, 0.*deg, 360.*deg );

    G4RotationMatrix *chamholerot = new G4RotationMatrix;
    chamholerot->rotateY(90.0*deg);

    //  Cut holes in the scattering chamber
    G4SubtractionSolid* swall_holes = new G4SubtractionSolid("swall_enthole", swallcut, swall_enthole, chamholerot, G4ThreeVector(-(swallrad+swallrad_in)/2, 0.0, 0.0) );
    swall_holes = new G4SubtractionSolid("swall_holes", swall_holes, swall_exthole, chamholerot, G4ThreeVector((swallrad+swallrad_in)/2, 0.0, 0.0) );

    G4LogicalVolume *swall_log = new G4LogicalVolume(swall_holes, GetMaterial("Aluminum"),"scham_wall_log");

    G4LogicalVolume *sc_hcalwin_log = new G4LogicalVolume(swall_hcalwin, GetMaterial("Aluminum"),"sc_hcalwin_log");
    G4LogicalVolume *sc_bbwin_log = new G4LogicalVolume(swall_bbwin, GetMaterial("Aluminum"),"sc_bbwin_log");

    G4RotationMatrix *schamrot = new G4RotationMatrix;
    schamrot->rotateX(-90.0*deg);
    schamrot->rotateZ(-90.0*deg);

    G4RotationMatrix *targrot = new G4RotationMatrix;
    targrot->rotateY(-90.0*deg);

    G4Tubs *chamber_inner = new G4Tubs("chamber_inner", 0.0, swallrad_in,  sheight/2, 0*deg, 360*deg );
    G4LogicalVolume* chamber_inner_log = new G4LogicalVolume(chamber_inner, GetMaterial("Vacuum"), "cham_inner_log");

    // Top and bottom
    G4Tubs *sc_topbottom = new G4Tubs("scham_topbottom", 0.0, swallrad, (swallrad-swallrad_in)/2, 0.*deg, 360.*deg );
    G4LogicalVolume* sc_topbottom_log = new G4LogicalVolume(sc_topbottom, GetMaterial("Aluminum"), "scham_topbottom_log");

    /*
     * FIXME*/
    if( fTargType == kLH2 || fTargType == kLD2 ){
	new G4PVPlacement(targrot, G4ThreeVector(0.0, 0.0, 0.0), targ_tube_log,
		"targ_tube_phys", chamber_inner_log, false, 0);

	/*
	   new G4PVPlacement(targrot, G4ThreeVector(fTargLen/2.0+downcapthick/2.0, 0.0, 0.0), targ_dcap_log,
	   "targ_dcap_phys", chamber_inner_log, false, 0);
	   new G4PVPlacement(targrot, G4ThreeVector(-fTargLen/2.0-upcapthick/2.0, 0.0, 0.0), targ_ucap_log,
	   "targ_ucap_phys", chamber_inner_log, false, 0);
	   */

	// Scattering chamber
	new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), swall_log,
		"scham_wall_phys", worldlog, false, 0);
	new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), chamber_inner_log,
		"chamber_inner_phys", worldlog, false, 0);

	new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), sc_hcalwin_log,
		"sc_hcalwin_phys", worldlog, false, 0);
	new G4PVPlacement(schamrot, G4ThreeVector(0.0, 0.0, 0.0), sc_bbwin_log,
		"sc_bbwin_phys", worldlog, false, 0);

	new G4PVPlacement(schamrot, G4ThreeVector(0.0, sheight/2.0 + (swallrad-swallrad_in)/2, 0.0), sc_topbottom_log,
		"scham_top_phys", worldlog, false, 0);
	new G4PVPlacement(schamrot, G4ThreeVector(0.0, -sheight/2.0 - (swallrad-swallrad_in)/2, 0.0), sc_topbottom_log,
		"scham_bot_phys", worldlog, false, 0);
    }
    /**/

    //  G4Tubs *cryo_tube = new G4Tubs("cryo_tube", 0.0, cellradius-wallthick,fTargLen/2.0, 0.*deg, 360.*deg );
    G4LogicalVolume* cryo_tube_log = NULL;


    if( fTargType == kLH2 ){
	cryo_tube_log = new G4LogicalVolume(cryovol, GetMaterial("LH2mat"), "cryo_tube_log");
    }
    if( fTargType == kLD2 ){
	cryo_tube_log = new G4LogicalVolume(cryovol, GetMaterial("LD2mat"), "cryo_tube_log");
    }

    /*
     * FIXME */
    if( fTargType == kLD2 || fTargType == kLH2 ){
	new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, cryooffset), cryo_tube_log,
		"cryo_tube_phys", targ_tube_log, false, 0);
    }
    /*  */


    //  Vis attributes
    chamber_inner_log->SetVisAttributes(G4VisAttributes::Invisible);
    G4VisAttributes * schamVisAtt
	= new G4VisAttributes(G4Colour(0.7,0.7,1.0));
    //      = new G4VisAttributes(G4VisAttributes::Invisible);
    swall_log->SetVisAttributes(schamVisAtt);
    sc_topbottom_log->SetVisAttributes(schamVisAtt);

    G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.6,0.6,0.6));
    extpipe_log->SetVisAttributes(pipeVisAtt);
    extvac_log->SetVisAttributes( G4VisAttributes::Invisible );

    G4VisAttributes *cryoVisAtt= new G4VisAttributes(G4Colour(1.0,0.0,0.0));
    if( cryo_tube_log ){
	cryo_tube_log->SetVisAttributes(cryoVisAtt);
    }

    G4VisAttributes *winVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
    sc_hcalwin_log->SetVisAttributes(winVisAtt);
    sc_bbwin_log->SetVisAttributes(winVisAtt);

    return;

}
