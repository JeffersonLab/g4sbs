#include "BeamLine.hh"
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
#include "G4Polycone.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "DetectorConstruction.hh"



Beamline::Beamline(DetectorConstruction* rectagg)
{
  fRtag = rectagg;
}

Beamline::~Beamline(){;}

void Beamline::BuildComponent(G4LogicalVolume *worldlog){
  //  Entry iron tube to shield beam from stray field
  //
  G4double tRmin = 5.*cm;
  G4double tRmax = 7.*cm;
  G4double tDzz = 20.*cm/2;
  G4double tSPhi = 0.*deg;
  G4double tDphi = 360.*deg;
  //G4RotationMatrix *rm_snout = new G4RotationMatrix();
  
  //rm_snout->rotateY(90*deg);
  //G4LogicalVolume* snoutvaclog = fDetCon->GetSnoutVacLog();

 G4Tubs *iron_ent_tube =
   new G4Tubs("iron_ent_tube", tRmin, tRmax, tDzz, tSPhi, tDphi);

  G4LogicalVolume *iron_ent_tube_log =
    new G4LogicalVolume(iron_ent_tube, fRtag->GetFe(), "iron_ent_tube",
			0, 0, 0 );

  //new G4PVPlacement(rm_snout,G4ThreeVector(110.936*cm+tDzz,0.0,0.0),
  //		    iron_ent_tube_log, "iron_ent_tube_phys", snoutvaclog ,
  //		    false,0,fOvLap);
  //  next entrance correction magnet
  G4Box* box_1 = new G4Box("EnMag_1",50.*cm/2., 50.*cm/2., 16.*cm/2.);
  G4Box* box_2 = new G4Box("EnMag_2",20.*cm/2., 30.*cm/2., 17.*cm/2.); // aperture to pass beam
  G4SubtractionSolid* EnMag = new G4SubtractionSolid("EnMag", box_1, box_2);   
  G4LogicalVolume * EnMag_log =
    new G4LogicalVolume(EnMag , fRtag->GetFe(), "EnMag", 0, 0, 0); 
  //new G4PVPlacement(rm_snout,G4ThreeVector(142.935*cm, 0.0, 0.0), EnMag_log,
  //		    "EnMag_phys", snoutvaclog, false, 0, fOvLap);

  // exit correction magnet 

  G4Box* box_3 = new G4Box("ExtMag_1",50.*cm/2., 54.*cm/2., 40.*cm/2.);
  G4Box* box_4 = new G4Box("ExtMag_2",25.5*cm/2., 40.*cm/2., 41.*cm/2.); // aperture to pass beam
  G4SubtractionSolid* ExtMag = new G4SubtractionSolid("ExtMag", box_3, box_4);   
  G4LogicalVolume * ExtMag_log =
    new G4LogicalVolume(ExtMag ,fRtag->GetFe(), "ExtMag", 0, 0, 0); 
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, 421.435*cm+20.*cm), ExtMag_log,
		    "ExtMag_phys", worldlog, false, 0, fOvLap);

  // iron conical tube on the beamline inside SBS magnet split

  G4double tcRmin1 = 7.4*cm;
  G4double tcRmax1 = 8.68*cm;
  G4double tcRmin2 = 10.55*cm;
  G4double tcRmax2 = 11.82*cm;
  G4double tcDzz = 120.*cm/2;
  G4double tcSPhi = 0.*deg;
  G4double tcDphi = 360.*deg;
  
 G4Cons *iron_con_tube =
   new G4Cons("iron_con_tube", tcRmin1, tcRmax1,tcRmin2, tcRmax2, tcDzz,
	      tcSPhi, tcDphi);

  G4LogicalVolume *iron_con_tube_log =
    new G4LogicalVolume(iron_con_tube,fRtag->GetFe(), "iron_con_tube",
			0, 0, 0 );

  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, 170.944*cm + tcDzz),
		    iron_con_tube_log, "iron_con_tube_phys", worldlog,
		    false,0,fOvLap);

  G4VisAttributes *McorrVisAtt= new G4VisAttributes(G4Colour(0.9,0.9,0.9));
  iron_con_tube_log->SetVisAttributes(McorrVisAtt);
  ExtMag_log->SetVisAttributes(McorrVisAtt);
  EnMag_log->SetVisAttributes(McorrVisAtt);
  iron_ent_tube_log->SetVisAttributes(McorrVisAtt);

  double swallrad = 1.143*m/2;
  double swallrad_inner = 1.041/2.0*m; 
  double beamheight = 10.0*12*2.54*cm; // 10 feet off the ground

  // Stainless
  G4double ent_len = 10*m;
  G4double ent_rin = 31.75*mm;
  G4double ent_rou = ent_rin+0.120*mm;

  G4Tubs *ent_tube = new G4Tubs("ent_tube", ent_rin, ent_rou, ent_len/2, 0.*deg, 360.*deg );
  G4Tubs *ent_vac  = new G4Tubs("ent_vac", 0.0, ent_rin, ent_len/2, 0.*deg, 360.*deg );

  //We want to subtract this cylinder from the entry tube/pipe:
  G4Tubs *cut_cylinder = new G4Tubs("cut_cylinder", 0.0, swallrad, 1.0*m, 0.0*deg, 360.0*deg );

  G4RotationMatrix *cut_cylinder_rot = new G4RotationMatrix;
  cut_cylinder_rot->rotateX( -90.0*deg );

  G4SubtractionSolid *ent_tube_cut = new G4SubtractionSolid( "ent_tube_cut", ent_tube, cut_cylinder, cut_cylinder_rot, 
							     G4ThreeVector( 0.0, 0.0, ent_len/2.0 + swallrad_inner ) );
  G4SubtractionSolid *ent_vac_cut = new G4SubtractionSolid( "ent_vac_cut", ent_vac, cut_cylinder, cut_cylinder_rot, 
							    G4ThreeVector( 0.0, 0.0, ent_len/2.0 + swallrad_inner ) );

  G4LogicalVolume *entLog =
    new G4LogicalVolume(ent_tube, fRtag->GetStainless(),"ent_log", 0, 0, 0);
  G4LogicalVolume *entvacLog =
    new G4LogicalVolume(ent_vac, fRtag->GetVac(), "entvac_log", 0, 0, 0);
  G4LogicalVolume *entLog_cut =
    new G4LogicalVolume(ent_tube_cut, fRtag->GetStainless(), "ent_log_cut", 0, 0, 0);
  G4LogicalVolume *entvacLog_cut =
    new G4LogicalVolume(ent_vac_cut, fRtag->GetVac(), "entvac_log_cut", 0, 0, 0);

    // Cryotarget - up against the chamber wall
  /*
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner),
		    entLog_cut, "ent_phys", worldlog, false,0,fOvLap);
  new G4PVPlacement(0,G4ThreeVector(0.0, 0.0, -ent_len/2-swallrad_inner),
		    entvacLog_cut, "entvac_phys", worldlog,false,0,fOvLap);
  */

  int nsec = 7;
  //  Definition taken from GEN_10M.opc by Bogdan to z = 5.92.  2mm thickness assumed
  G4double exit_z[]   = { 162.2*cm, 592.2*cm, 609.84*cm,609.85*cm, 1161.02*cm, 1161.03*cm,2725.66*cm };
  G4double exit_z_vac[] = { 162.5*cm, 592.5*cm, 610.24*cm,610.35*cm, 1161.52*cm, 1161.53*cm,2726.46*cm };

  G4double exit_zero[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  G4double exit_rin[] = { 4.8*cm, 14.8*cm,15.24*cm, 30.48*cm,  30.48*cm,45.72*cm, 45.72*cm };
  G4double exit_rou[] = { 5.0*cm, 15.0*cm,15.558*cm,30.798*cm,30.798*cm, 46.038*cm, 46.038*cm  };


  G4Polycone *ext_cone = new G4Polycone("ext_cone", 0.0*deg, 360.0*deg, nsec, exit_z, exit_rin, exit_rou);
  G4Polycone *ext_vac  = new G4Polycone("ext_vac ", 0.0*deg, 360.0*deg, nsec, exit_z_vac, exit_zero, exit_rin);

  G4LogicalVolume *extLog =
    new G4LogicalVolume(ext_cone, fRtag->GetAl(), "ext_log", 0, 0, 0);
  G4LogicalVolume *extvacLog =
    new G4LogicalVolume(ext_vac, fRtag->GetVac(), "extvac_log", 0, 0, 0);

  new G4PVPlacement(0,G4ThreeVector(), extLog, "ext_phys", worldlog, false,0,fOvLap);
  new G4PVPlacement(0,G4ThreeVector(), extvacLog, "extvac_phys", worldlog,false,0,fOvLap);


  G4VisAttributes *extVisAtt= new G4VisAttributes(G4Colour(0.9,0.1,0.9));
  extLog->SetVisAttributes(extVisAtt);
  /*
  double floorthick = 1.0*m;
  G4Tubs *floor_tube = new G4Tubs("floor_tube", 0.0, 20*m, floorthick/2, 0.*deg, 360.*deg );

  G4RotationMatrix *floorrm = new G4RotationMatrix;
  floorrm->rotateX(90*deg);

  G4LogicalVolume *floorLog =
    new G4LogicalVolume(floor_tube, fRtag->GetConcrete(), "floor_log", 0, 0, 0);
  new G4PVPlacement(floorrm, G4ThreeVector(0.0, -floorthick/2 - beamheight, 0.0), floorLog, "floor_phys", worldlog, false, 0, fOvLap);
  */

  extvacLog->SetVisAttributes(G4VisAttributes::Invisible);
  entvacLog->SetVisAttributes(G4VisAttributes::Invisible);

  entvacLog_cut->SetVisAttributes(G4VisAttributes::Invisible);

  G4VisAttributes *pipeVisAtt= new G4VisAttributes(G4Colour(0.2,0.6,0.2));

  extLog->SetVisAttributes(pipeVisAtt);
  entLog->SetVisAttributes(pipeVisAtt);
  entLog_cut->SetVisAttributes(pipeVisAtt);

  /*    G4VisAttributes *floorVisAtt= new G4VisAttributes(G4Colour(0.9,0.9,0.9));
	floorLog->SetVisAttributes(floorVisAtt); */
  //floorLog->SetVisAttributes(G4VisAttributes::Invisible);
  
  //MakeGEnClamp(worldlog);
  //MakeGEnLead(worldlog);
  return;
}

