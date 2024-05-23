#include "G4SBSHArmBuilder.hh"
#include "G4SBSDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4LogicalBorderSurface.hh"
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

#include "G4SBSCDet.hh"
#include "G4SBSCalSD.hh"
#include "G4SBSECalSD.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "sbstypes.hh"

#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"
#include "G4Mag_SpinEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4DormandPrince745.hh"
#include "TString.h"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <vector>
#include <map>
#include <stdlib.h>
#include <cmath>

using namespace std;

G4SBSHArmBuilder::G4SBSHArmBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
  f48D48ang  = 39.4*deg;
  f48D48dist = 2.8*m;
  f48D48_fieldclamp_config = 2; //0 = No field clamps. 2 = GEp (default). 1 = BigBite experiments:

  fHCALdist  = 17.0*m;
  fLACdist   = 14.0*m;
  fHCALvertical_offset = 0.0*cm;
  fHCALhorizontal_offset = 0.0*cm;
  fHCALangular_offset = 0.0*deg;

  fLACvertical_offset = 0.0*cm;
  fLAChorizontal_offset = 0.0*cm;

  fRICHdist  = 5.5*m;
  fRICHhorizontal_offset = 0.0*cm;
  fRICHvertical_offset = 0.0*cm;
  fRICHgas   = "C4F8_gas"; //default to C4F8;
  fRICH_use_aerogel = true;
  fRichSnoutExtension = 0;
  

  fGEnRP_analyzer_option = 1; //0=none+no beamline PR; 1=none(default), 2=Cu+Gla(para), 3=Cu+Gla(perp), 4=Cu+CGEN

  f48D48depth = 1219.2*mm;
  f48D48width = 2324.1*mm;
  f48D48height = 3721.1*mm;

  fUseLocalField = false;
  fFieldStrength = 1.4*tesla;

  fSBS_tracker_dist = 5.1*m;
  fSBS_tracker_pitch = 0.0*deg;
  
  fBuildSBSSieve = false;

  fCH2thickFPP[0] = 22.0*2.54*cm;
  fCH2thickFPP[1] = 22.0*2.54*cm;

  //Default GEP FPP configuration:
  fGEPFPPoption = 2; // default configuration is two 22-inch CH2 blocks with trackers in between:
                     //alternate configuration to be tested is one analyzer with the FPP1 GEMs re-distributed to FT and FPP2

  //default to two analyzers, with same transverse dimensions as the polarimeter GEMs, each 22 inches thick:
  fGEP_CH2width[0] = fGEP_CH2width[1] = 28.88*2.54*cm;
  fGEP_CH2height[0] = fGEP_CH2height[1] = 84.0*2.54*cm;

  fGEP_CH2yoff[0] = fGEP_CH2yoff[1] = 0.0*cm;
  //fGEP_trkryoff[0] = fGEP_trkryoff[1] = fGEP_trkryoff[2] = 0.0*cm;
  
  fGEP_CH2zpos[0] = 58.53*cm;
  fGEP_CH2zpos[1] = 170.3*cm;

  fGEP_CH2yoff[0] = fGEP_CH2yoff[1] = 0.0*cm;
  
  
  fFTuseabsorber = false;
  fFTabsthick = 2.54*cm;
  fFTabsmaterial = "Aluminum";
  
  //fCDetReady = true;
  fCDetReady = false; //default to false for both of these
  
  fUseNeutronVeto=false;
  
  assert(fDetCon);
}

G4SBSHArmBuilder::~G4SBSHArmBuilder(){
}

void G4SBSHArmBuilder::BuildComponent(G4LogicalVolume *worldlog){
  G4SBS::Exp_t exptype = fDetCon->fExpType;
  G4SBS::Targ_t tgttype = fDetCon->fTargType;

  // Build the 48D48 magnet and HCAL
  // All three types of experiments have a 48D48 magnet:
  if( exptype != G4SBS::kC16 && exptype != G4SBS::kGEMHCtest ) {
    Make48D48(worldlog, f48D48dist + f48D48depth/2. );
    if(fBuildSBSSieve)
      MakeSBSSieveSlit(worldlog);
    
    //MakeHCAL( worldlog, fHCALvertical_offset );
    if( exptype != G4SBS::kA1n  && exptype != G4SBS::kTDIS && exptype != G4SBS::kNDVCS )
      MakeHCALV2( worldlog, fHCALvertical_offset );
  }

  // Now build special components for experiments
  if( exptype == G4SBS::kSIDISExp) {
    //SIDIS experiment requires a RICH detector and a tracker for SBS: 
    MakeTracker(worldlog);
    //MakeRICH( worldlog );
    MakeRICH_new( worldlog );
  } else if ( exptype == G4SBS::kGEp || exptype == G4SBS::kGEPpositron || exptype == G4SBS::kGEp_BB ) {
    //Subsystems unique to the GEp experiment include FPP and BigCal:
    MakeGEpFPP(worldlog);
  } 
  
  if( exptype == G4SBS::kALL) {
    
    MakeTracker(worldlog, 8);
    
  }
  

  if ( exptype == G4SBS::kA1n  || exptype == G4SBS::kTDIS || exptype == G4SBS::kNDVCS ) {
    //A1n case is similar to SIDIS, except now we want to use the SBS in
    //"electron mode"; meaning we want to remove the aerogel from the RICH,
    //and replace the RICH gas with CO2, and we also want to have a non-zero
    //pitch angle for the SBS tracker. We assume (for NOW) that the RICH can
    //be supported at some non-zero "pitch" angle:
    
    //TDIS and nDVCS will also use SBS in electron mode
    MakeElectronModeSBS(worldlog);
    //MakeTracker_A1n(worldlog);
    //MakeRICH_new( worldlog );
  }

  // Build Test detector ****************************************
  if( exptype == G4SBS::kC16 ) {
    MakeTest( worldlog );
  }

  // Build CDET (as needed)
  //if( (exptype == G4SBS::kGMN || exptype == G4SBS::kGEnRP ) && (tgttype==kLH2 || tgttype==kLD2)){
  if( exptype == G4SBS::kGMN || exptype == G4SBS::kGEnRP ){
    //plugging in CDET for GMn  
    G4double depth_HCal_shield = 7.62*cm; //3 inches
    G4double depth_CH2 = 20.0*cm; //This goes directly in front of CDET:
    G4double depth_CDET = 40.0*cm;
    G4double dist_from_hcal = 35.0*cm;
    
    G4Box *CH2_filter = new G4Box( "CH2_filter", 150.0*cm/2.0, 340.0*cm/2.0, depth_CH2/2.0 );
    G4LogicalVolume *CH2_filter_log = new G4LogicalVolume( CH2_filter, GetMaterial("Polyethylene"), "CH2_filter_log" );

    G4RotationMatrix *HArmRot = new G4RotationMatrix;
    HArmRot->rotateY(f48D48ang);

    G4ThreeVector CH2_pos( ( fHCALdist-dist_from_hcal-(depth_CDET+depth_CH2)/2.0 ) * sin( -f48D48ang ),  fHCALvertical_offset, ( fHCALdist-dist_from_hcal-(depth_CDET+depth_CH2)/2.0 ) * cos( -f48D48ang ) );

    //new G4PVPlacement( HArmRot, CH2_pos, CH2_filter_log, "CH2_filter_phys", worldlog, false, 0 );
    
    if(fCDetReady){

      new G4PVPlacement( HArmRot, CH2_pos, CH2_filter_log, "CH2_filter_phys", worldlog, false, 0 );
      
      G4Box* CDetmother = new G4Box("CDetmother", 2.5*m/2.0, 3.0*m/2, depth_CDET/2.0);
      G4LogicalVolume *CDetmother_log = new G4LogicalVolume( CDetmother, GetMaterial("Air"), "CDetmother_log" );
      
      HArmRot->rotateY(180.0*deg);
      
      //G4ThreeVector CDetmother_pos( ( fHCALdist-0.30*m ) * sin( -f48D48ang ),  fHCALvertical_offset, ( fHCALdist-0.30*m ) * cos( -f48D48ang ) );
      G4ThreeVector CDetmother_pos( -15.0*cm,  fHCALvertical_offset, fHCALdist-dist_from_hcal+5.0*cm );
      CDetmother_pos.rotateY(-f48D48ang);
      new G4PVPlacement(HArmRot, CDetmother_pos, CDetmother_log, "CDetmother_phys", worldlog, false, 0);

      G4VisAttributes *CH2_visatt = new G4VisAttributes( G4Colour( 0, 0.6, 0.6 ) );
      CH2_visatt->SetForceWireframe(true);
      CH2_filter_log->SetVisAttributes(CH2_visatt);
      
      //G4VisAttributes* VisAtt = new G4VisAttributes( G4Colour(1, 1, 1) );
      //VisAtt->SetForceWireframe(true);
      //CDetmother_log->SetVisAttributes( VisAtt );
      CDetmother_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
      
      G4double z0_CDET = -0.15*m;
      G4double planes_hoffset = 0.84*m;
      G4double planes_interdist = 1.0*cm;//20.0*cm;
      
      G4SBSCDet* CDet = new G4SBSCDet(fDetCon);
      CDet->SetArmName("Harm");
      CDet->SetR0(fHCALdist + z0_CDET);
      CDet->SetZ0(z0_CDET);
      CDet->SetPlanesHOffset(planes_hoffset);
      CDet->SetPlanesInterDistance(planes_interdist);
      CDet->BuildComponent( CDetmother_log );
      //MakeCDET( CDetmother_log, z0_CDET, planes_hoffset );
    }else if( fUseNeutronVeto ){
      new G4PVPlacement( HArmRot, CH2_pos, CH2_filter_log, "CH2_filter_phys", worldlog, false, 0 );
      
      MakeNeutronVeto(worldlog, dist_from_hcal);
    }
    //Add a plate on the side of HCal
    
    double SideShield_Width = 1574.75*mm;
    double SideShield_Height = 3716.02*mm;
    double SideShield_Thick = 38.1*mm;
    double SideShield_XOffset = +1898.65*mm/2.0+SideShield_Thick/2.0;
    
    G4Box *HCal_sideshield = new G4Box( "HCal_sideshield", SideShield_Thick/2.0, SideShield_Height/2.0, SideShield_Width/2.0 );
    G4LogicalVolume *HCal_sideshield_log = new G4LogicalVolume( HCal_sideshield, GetMaterial("Steel"), "HCal_sideshield_log" );
    //HCal_sideshield_log
    G4ThreeVector HCal_sideshield_pos( SideShield_XOffset,  fHCALvertical_offset, fHCALdist+SideShield_Width/2.0 );
    HCal_sideshield_pos.rotateY(-f48D48ang);
    
    G4bool checkOverlap = fDetCon->fCheckOverlap;

    //new G4PVPlacement( HArmRot, HCal_sideshield_pos, HCal_sideshield_log, "HCal_sideshield_phys", worldlog, false, 0, checkOverlap );
    HCal_sideshield_log->SetVisAttributes( G4Colour(0.7, 0.7, 0.7) );
  }

  // Build GEn-RP polarimeter
  if( exptype == G4SBS::kGEnRP || exptype == G4SBS::kGEN) 
    MakePolarimeterGEnRP( worldlog );
}

void G4SBSHArmBuilder::MakeGEpFPP(G4LogicalVolume *worldlog)
{

  G4SBS::Exp_t exptype = fDetCon->fExpType;
  
  //If exptype is GEp_BB option, we build a slightly different version of the FPP:
  
  //Let's make a box and then put the FPP in it:
  //Define the rotation matrix for the FPP (pitch angle of 5 deg relative to vertical): 
  G4double sbsboxpitch = 5.0*deg;

  if( exptype == G4SBS::kGEp_BB ){
    sbsboxpitch = 3.0*deg;
  }

  G4double sbsvoff = 0.0*cm;
  
  if( fGEPFPPoption != 2 ){
    sbsboxpitch = 0.0*deg;
    //sbsvoff = 10.0*cm;
  }
  
  SetTrackerPitch( sbsboxpitch );
  
  G4RotationMatrix *SBS_FPP_rm = new G4RotationMatrix;
  SBS_FPP_rm->rotateY( f48D48ang );
  SBS_FPP_rm->rotateX( sbsboxpitch );

  //FPP box: 
  double sbsdepth  = 3.5*m;
  double sbswidth  = 2.0*m;
  double sbsheight = 3.5*m;

  //double sbsr = fHCALdist - 4.106*m + sbsheight*sin(sbsboxpitch)/2.0 + sbsdepth/2.0;
  //double sbsr = f48D48dist + 1.694*m + sbsheight*sin(sbsboxpitch)/2.0 + sbsdepth/2.0;
  double sbsr = f48D48dist + 48.0*2.54*cm + 0.82*m + sbsdepth/2.0;
  
  SetTrackerDist( sbsr );
  
  G4Box *sbsbox = new G4Box("sbsbox", sbswidth/2.0, sbsheight/2.0, sbsdepth/2.0 );
  G4LogicalVolume* sbslog = new G4LogicalVolume(sbsbox, GetMaterial("Air"), "sbslog");

  sbslog->SetVisAttributes( G4VisAttributes::GetInvisible() );
  //Now position and orient the FPP "box":

  G4ThreeVector sbsbox_pos( -sbsr*sin(f48D48ang), (sbsr-f48D48dist)*sin(sbsboxpitch) + sbsvoff, sbsr*cos(f48D48ang) );
  
  new G4PVPlacement(SBS_FPP_rm, sbsbox_pos, sbslog,
		    "sbsphys", worldlog, false, 0, false);

  G4RotationMatrix *rot_I = new G4RotationMatrix;

  //First GEM is 5 cm downstream of edge of SBS box. 
  double detoffset = 0.05*m - sbsdepth/2.0;

  MakeFPP( sbslog, rot_I, G4ThreeVector( 0.0, 0.0, detoffset) );

  //optional absorber in front of FT:
  if( fFTuseabsorber ){
    G4Box *FTabs_box = new G4Box( "FTabs_box", 20.0*cm, 75.0*cm, fFTabsthick/2.0 );

    G4LogicalVolume *FTabs_log = new G4LogicalVolume( FTabs_box,  GetMaterial(fFTabsmaterial), "FTabs_log" );
    
    G4ThreeVector SBS_FPP_zaxis(0,0,1);

    SBS_FPP_zaxis *= SBS_FPP_rm->inverse();

    G4cout << "SBS FPP z axis = ("
	   << SBS_FPP_zaxis.x() << ", "
	   << SBS_FPP_zaxis.y() << ", "
	   << SBS_FPP_zaxis.z() << ")" << G4endl;
    
    G4ThreeVector FTabs_pos = sbsbox_pos - (sbsdepth/2.0 + 1.0*mm + fFTabsthick/2.0)*SBS_FPP_zaxis;

    new G4PVPlacement( SBS_FPP_rm, FTabs_pos, FTabs_log, "FTabs_phys", worldlog, false, 0, false );
  }
}

void G4SBSHArmBuilder::Make48D48( G4LogicalVolume *worldlog, double r48d48 ){

  G4RotationMatrix *bigrm = new G4RotationMatrix;
  //bigrm->rotateY(-f48D48ang);
  bigrm->rotateY(f48D48ang);
  G4ThreeVector SBS_xaxis( cos(f48D48ang), 0, sin(f48D48ang) );
  G4ThreeVector SBS_yaxis(0,1,0);
  G4ThreeVector SBS_zaxis( -sin(f48D48ang), 0, cos(f48D48ang) );
  
  G4String name;

  double bigcoilwidth = 214.5*mm;
  double bigcoilheight = 263.7*mm;

  double notchdepth = 25*cm;

  G4Box *biggap  = new G4Box("biggap",  469.9*mm/2+0.1*mm, 187.*cm/2.-bigcoilheight,  f48D48depth/2+0.1*mm);

  std::vector<G4TwoVector> bigpoly;
  // bigpoly.push_back( G4TwoVector(-f48D48width/2.0,  f48D48depth/2.0 ));
  // bigpoly.push_back( G4TwoVector(-f48D48width/2.0, -f48D48depth/2.0 ));
  // bigpoly.push_back( G4TwoVector( f48D48width/2.0, -f48D48depth/2.0 ));
  // bigpoly.push_back( G4TwoVector( f48D48width/2.0, f48D48depth/2.0 - notchdepth*sqrt(2.0) ));
  // bigpoly.push_back( G4TwoVector( f48D48width/2.0 - notchdepth*sqrt(2.0), f48D48depth/2.0  ));

  //How does this change when moving things to beam right?
  //polygon should be in clockwise order:
  //let the z axis of the extruded solid point vertically upward. Start with top right point:
  bigpoly.push_back( G4TwoVector( -f48D48width/2.0, f48D48depth/2.0 ) ); //back right
  bigpoly.push_back( G4TwoVector( -f48D48width/2.0, -f48D48depth/2.0 ) ); //front right
  bigpoly.push_back( G4TwoVector( f48D48width/2.0 - notchdepth*sqrt(2.0), -f48D48depth/2.0 ) ); //front left:
  bigpoly.push_back( G4TwoVector( f48D48width/2.0, -f48D48depth/2.0 + notchdepth*sqrt(2.0) ) ); //front left
  bigpoly.push_back( G4TwoVector( f48D48width/2.0, f48D48depth/2.0 ) ); //back left:

  G4ExtrudedSolid *bigbox_ext = new G4ExtrudedSolid("bigbox_ext", bigpoly, f48D48height/2.0, G4TwoVector(), 1.0, G4TwoVector(), 1.0);
  G4RotationMatrix *bigboxrm = new G4RotationMatrix;

  //We want the z axis to point vertically up!
  bigboxrm->rotateY(f48D48ang);
  bigboxrm->rotateX( -90.*deg);
  //bigboxrm->rotateZ( 180.*deg);

  G4RotationMatrix *bigboxaddrm = new G4RotationMatrix;
  //bigboxaddrm->rotateZ( -180.*deg);
  bigboxaddrm->rotateX( 90.*deg);

  //moved definition of clamp gaps to field clamp method:
  // G4Box *bclampgap  = new G4Box("bclampgap",  23.*cm, 65.*cm,  12.*cm/2.);
  // G4Box *fclampgap  = new G4Box("fclampgap",  11.*cm, 35.*cm,  12.*cm/2.);

  G4SubtractionSolid* bigbase = new G4SubtractionSolid("bigbase", bigbox_ext, biggap, bigboxaddrm, G4ThreeVector());

  double coilgapwidth = (60.*cm - bigcoilheight)*2; //= 120 cm - 2*coilheight
  double coilgapheight = 160*cm-bigcoilheight; //=160 cm - 1*coilheight

  //bigcoilbase has width of (height+ 60 cm - height ) = 60cm,
  //height/2 = bheight + 80cm - bheight/2 = 80 cm + bheight/2 --> height = 160 cm + bheight
  //depth = bigcoil width
  //bigcoil gap has width of coilgapwidth/2 + 2 mm = 60 cm - bheight + 2 mm

  //G4Box *bigcoilbase = new G4Box( "bigcoilbase", 60.0*cm/2.0, 
  
  G4Box *bigcoilbase = new G4Box("bigcoilbase", (bigcoilheight+coilgapwidth/2.0)/2.0, bigcoilheight+coilgapheight/2, bigcoilwidth/2.0);
  G4Box *bigcoilgap = new G4Box("bigcoilgap", coilgapwidth/4.0, coilgapheight/2, bigcoilwidth/2.0+0.1*mm);

  //  double coilspace = 6.63*mm;
  double coilspace = 20.63*mm;

  std::vector<G4TwoVector> woundpoly_outer;
  std::vector<G4TwoVector> woundpoly_inner;
  
  woundpoly_outer.push_back( G4TwoVector( 0.0, f48D48depth/2.0 + coilspace + bigcoilwidth ) );
  woundpoly_outer.push_back( G4TwoVector( f48D48width/2.0 + coilspace + bigcoilwidth, f48D48depth/2.0 + coilspace + bigcoilwidth ) );
  woundpoly_outer.push_back( G4TwoVector( f48D48width/2.0 + coilspace + bigcoilwidth, -f48D48depth/2.0 + notchdepth*sqrt(2.0) ) );
  woundpoly_outer.push_back( G4TwoVector( f48D48width/2.0 + coilspace + bigcoilwidth - 2.0*(coilspace+bigcoilwidth)*sin(pi/8.)*sin(pi/8.),
					  -f48D48depth/2.0 + notchdepth*sqrt(2.0) - 2.0*(coilspace+bigcoilwidth)*sin(pi/8.)*cos(pi/8.) ) );
  woundpoly_outer.push_back( G4TwoVector( f48D48width/2.0 - notchdepth*sqrt(2.0) + 2.0*(coilspace+bigcoilwidth)*sin(pi/8.)*cos(pi/8.),
					  -f48D48depth/2.0 - coilspace - bigcoilwidth + 2.0*(coilspace+bigcoilwidth)*sin(pi/8.)*sin(pi/8.) ) );
  woundpoly_outer.push_back( G4TwoVector( f48D48width/2.0 - notchdepth*sqrt(2.0), -f48D48depth/2.0 - coilspace - bigcoilwidth ) );
  woundpoly_outer.push_back( G4TwoVector( 0.0, -f48D48depth/2.0 - coilspace - bigcoilwidth ) );

  //inner polygon: coil will be subtraction outer - inner:
  woundpoly_inner.push_back( G4TwoVector( bigcoilwidth, f48D48depth/2.0 + coilspace ) );
  woundpoly_inner.push_back( G4TwoVector( f48D48width/2.0 + coilspace, f48D48depth/2.0 + coilspace ) );
  woundpoly_inner.push_back( G4TwoVector( f48D48width/2.0 + coilspace, -f48D48depth/2.0 + notchdepth*sqrt(2.0) ) );
  woundpoly_inner.push_back( G4TwoVector( f48D48width/2.0 + coilspace - 2.0*(coilspace)*sin(pi/8.)*sin(pi/8.),
					  -f48D48depth/2.0 + notchdepth*sqrt(2.0) - 2.0*(coilspace)*sin(pi/8.)*cos(pi/8.) ) );
  woundpoly_inner.push_back( G4TwoVector( f48D48width/2.0 - notchdepth*sqrt(2.0) + 2.0*(coilspace)*sin(pi/8.)*cos(pi/8.),
					  -f48D48depth/2.0 - coilspace + 2.0*(coilspace)*sin(pi/8.)*sin(pi/8.) ) );
  woundpoly_inner.push_back( G4TwoVector( f48D48width/2.0 - notchdepth*sqrt(2.0), -f48D48depth/2.0 - coilspace ) );
  woundpoly_inner.push_back( G4TwoVector( bigcoilwidth, -f48D48depth/2.0 - coilspace ) );
  
  G4ExtrudedSolid *woundcoil_outer = new G4ExtrudedSolid( "woundcoil_outer", woundpoly_outer, bigcoilheight/2.0, G4TwoVector(), 1.0, G4TwoVector(), 1.0);
  G4ExtrudedSolid *woundcoil_inner = new G4ExtrudedSolid( "woundcoil_inner", woundpoly_inner, bigcoilheight/2.0+mm, G4TwoVector(), 1.0, G4TwoVector(), 1.0);
  G4SubtractionSolid *woundcoil = new G4SubtractionSolid( "woundcoil", woundcoil_outer, woundcoil_inner, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *woundcoil_log = new G4LogicalVolume( woundcoil, GetMaterial("Copper"), "woundcoil_log" );
  G4VisAttributes *coil_color = new G4VisAttributes( G4Colour( 0.5, 0.15, 0.75 ) );
  woundcoil_log->SetVisAttributes( coil_color );

  G4RotationMatrix *rot_coil = new G4RotationMatrix;
  rot_coil->rotateY(f48D48ang);
  rot_coil->rotateX(-90.0*deg);

  if( fDetCon->fTotalAbs ){
    woundcoil_log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }
  
  G4ThreeVector woundcoil_pos1 = (coilgapheight+bigcoilheight)/2.0 * SBS_yaxis + r48d48 * SBS_zaxis;
  new G4PVPlacement( rot_coil, woundcoil_pos1, woundcoil_log, "woundcoil_top_phys", worldlog, false, 0 );
  G4ThreeVector woundcoil_pos2 = -(coilgapheight+bigcoilheight)/2.0 * SBS_yaxis + r48d48 * SBS_zaxis;
  new G4PVPlacement( rot_coil, woundcoil_pos2, woundcoil_log, "woundcoil_bottom_phys", worldlog, false, 1 );
  //////////////////////		       
  
  // woundpoly.push_back( G4TwoVector(0.0,  -f48D48depth/2.0 -coilspace )); 
  // woundpoly.push_back( G4TwoVector(0.0, -f48D48depth/2.0 -coilspace-bigcoilwidth));
  // woundpoly.push_back( G4TwoVector( f48D48width/2.0+coilspace+bigcoilwidth, -f48D48depth/2.0 -coilspace-bigcoilwidth));
  // woundpoly.push_back( G4TwoVector( f48D48width/2.0+coilspace+bigcoilwidth, f48D48depth/2.0 - notchdepth*sqrt(2.0) ));

  // woundpoly.push_back( G4TwoVector( f48D48width/2.0+coilspace+bigcoilwidth  - 2.0*(coilspace+bigcoilwidth)*sin(pi/8.)*sin(pi/8.) , 
  // 				    f48D48depth/2.0 - notchdepth*sqrt(2.0) + 2.0*(coilspace+bigcoilwidth)*sin(pi/8)*cos(pi/8) ));

  // woundpoly.push_back( G4TwoVector( f48D48width/2.0- notchdepth*sqrt(2.0) + 2.0*(coilspace+bigcoilwidth)*sin(pi/8)*cos(pi/8) , 
  // 				    f48D48depth/2.0 + coilspace+bigcoilwidth  - 2.0*(coilspace+bigcoilwidth)*sin(pi/8)*sin(pi/8) ));

  // woundpoly.push_back( G4TwoVector( f48D48width/2.0 - notchdepth*sqrt(2.0), f48D48depth/2.0 +coilspace+bigcoilwidth ));

  // ////
  // woundpoly.push_back( G4TwoVector(0.0,  f48D48depth/2.0 +coilspace+bigcoilwidth));
  // woundpoly.push_back( G4TwoVector(0.0,  f48D48depth/2.0 +coilspace ));

  // woundpoly.push_back( G4TwoVector( f48D48width/2.0 - notchdepth*sqrt(2.0), f48D48depth/2.0 +coilspace ));

  // // arc here
  // woundpoly.push_back( G4TwoVector( f48D48width/2.0+coilspace  - 2.0*(coilspace)*sin(pi/8.)*sin(pi/8.) , 
  // 				    f48D48depth/2.0 - notchdepth*sqrt(2.0) + 2.0*(coilspace)*sin(pi/8)*cos(pi/8) ));

  // woundpoly.push_back( G4TwoVector( f48D48width/2.0- notchdepth*sqrt(2.0) + 2.0*(coilspace)*sin(pi/8)*cos(pi/8) , 
  // 				    f48D48depth/2.0 + coilspace - 2.0*(coilspace)*sin(pi/8)*sin(pi/8) ));

  // woundpoly.push_back( G4TwoVector( f48D48width/2.0+coilspace, f48D48depth/2.0 - notchdepth*sqrt(2.0) ));

  // woundpoly.push_back( G4TwoVector( f48D48width/2.0+coilspace, f48D48depth/2.0 +coilspace - notchdepth*sqrt(2.0) ));
  // woundpoly.push_back( G4TwoVector( f48D48width/2.0+coilspace, -f48D48depth/2.0 -coilspace));

  //G4ExtrudedSolid *woundcoil_ext = new G4ExtrudedSolid("woundcoil_ext", woundpoly, bigcoilheight/2.0, G4TwoVector(), 1.0, G4TwoVector(), 1.0);

  // Pull out left side of gap
  G4SubtractionSolid *bigcoil = new G4SubtractionSolid("bigcoil", bigcoilbase, bigcoilgap, 0, G4ThreeVector(coilgapwidth/4.0, 0.0, 0.0) );
  G4LogicalVolume *bigcoil_log = new G4LogicalVolume(bigcoil, GetMaterial("Copper"), "bigcoil_log" );
  bigcoil_log->SetVisAttributes( coil_color );

  if( fDetCon->fTotalAbs ){
    bigcoil_log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }
  
  G4ThreeVector bigcoil_pos1 = -(bigcoilheight+coilgapwidth/2.0)/2.0*SBS_xaxis + (f48D48dist - coilspace - bigcoilwidth/2.0)*SBS_zaxis;
  G4ThreeVector bigcoil_pos2 = -(bigcoilheight+coilgapwidth/2.0)/2.0*SBS_xaxis + (f48D48dist + f48D48depth + coilspace + bigcoilwidth/2.0)*SBS_zaxis;
  
  new G4PVPlacement( bigrm, bigcoil_pos1, bigcoil_log, "bigcoil_front_phys", worldlog, false, 0 );
  new G4PVPlacement( bigrm, bigcoil_pos2, bigcoil_log, "bigcoil_back_phys", worldlog, false, 1 );
  //  double coilfrontback = 150*cm;
  double coilfrontback = f48D48depth+2.0*coilspace;

  G4Box *bigcoilthr = new G4Box("bigcoilthr", bigcoilwidth/2.0,  bigcoilheight/2,  coilfrontback/2.0 );
  G4LogicalVolume *bigcoilthr_log = new G4LogicalVolume(bigcoilthr, GetMaterial("Copper"), "bigcoilthr_log" );

  if( fDetCon->fTotalAbs ){
    bigcoilthr_log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }
  
  bigcoilthr_log->SetVisAttributes( coil_color );

  G4ThreeVector bigcoilthr_pos1 = -bigcoilwidth/2.0 * SBS_xaxis + (coilgapheight+bigcoilheight)/2.0 * SBS_yaxis + r48d48 * SBS_zaxis;
  new G4PVPlacement( bigrm, bigcoilthr_pos1, bigcoilthr_log, "bigcoilthr_phys_top", worldlog, false, 0 );
  G4ThreeVector bigcoilthr_pos2 = -bigcoilwidth/2.0 * SBS_xaxis - (coilgapheight+bigcoilheight)/2.0 * SBS_yaxis + r48d48 * SBS_zaxis;
  new G4PVPlacement( bigrm, bigcoilthr_pos2, bigcoilthr_log, "bigcoilthr_phys_bottom", worldlog, false, 1 );
  
  // Sum together coils


  // Sum together base iron plus coils

  //G4UnionSolid* big48d48;

  //G4Box *bigbeamslot = new G4Box("bigbeamslot",  f48D48width/2, 15.5*cm, 2.0*m ); 
  G4Box *bigbeamslot = new G4Box("bigbeamslot", f48D48width/2.0, 2.0*m, 15.5*cm );
  
  G4double slot_angle = 6.71*deg;
  G4double slot_depth_front = 58.85*2.54*cm;

  G4ThreeVector beamslot_xaxis( cos(slot_angle), -sin(slot_angle), 0 );
  G4ThreeVector beamslot_yaxis( sin(slot_angle), cos(slot_angle), 0 );
  G4ThreeVector beamslot_zaxis = (beamslot_xaxis.cross(beamslot_yaxis)).unit();

  G4ThreeVector beamslot_posrel( -f48D48width/2.0 + slot_depth_front + f48D48depth/2.0 * tan(slot_angle), 0, 0 );
  beamslot_posrel += f48D48width/2.0 * beamslot_xaxis;
  
  G4RotationMatrix *rot_temp = new G4RotationMatrix;
  rot_temp->rotateZ( slot_angle );
  

  // big48d48 = new G4UnionSolid("big48d48_1", bigbase, bigcoilthr, bigboxaddrm, 
  // 			      G4ThreeVector(0.0, 0.0, (coilgapheight+bigcoilheight)/2.0));
  // big48d48 = new G4UnionSolid("big48d48_2", big48d48, bigcoilthr, bigboxaddrm, 
  // 			      G4ThreeVector(0.0, 0.0, -(coilgapheight+bigcoilheight)/2.0));

  // big48d48 = new G4UnionSolid("big48d48_3", big48d48, bigcoil, bigboxaddrm, 
  // 			      G4ThreeVector(-(bigcoilheight+coilgapwidth/2.0)/2.0-1.0*mm, coilfrontback/2.0, 0.0));
  // big48d48 = new G4UnionSolid("big48d48_4", big48d48, bigcoil, bigboxaddrm, 
  // 			      G4ThreeVector(-(bigcoilheight+coilgapwidth/2.0)/2.0-1.0*mm, -coilfrontback/2.0, 0.0));

  // big48d48 = new G4UnionSolid("big48d48_5", big48d48, woundcoil_ext, 0, 
  // 			      G4ThreeVector( 1.0*mm, 0.0,  coilgapheight/2.+bigcoilheight/2.0));
  // big48d48 = new G4UnionSolid("big48d48_6", big48d48, woundcoil_ext, 0, 
  // 			      G4ThreeVector( 1.0*mm, 0.0,  -coilgapheight/2.-bigcoilheight/2.0));


  //  Cut out slot - from magnet center to inside of cut is ~35cm
  G4SubtractionSolid *big48d48_wslot = new G4SubtractionSolid("big48d48_wslot", bigbase, bigbeamslot, rot_temp, 
							      beamslot_posrel );

  G4LogicalVolume *big48d48Log=new G4LogicalVolume(big48d48_wslot, GetMaterial("Fer"),
						   "b48d48Log", 0, 0, 0);

  G4VisAttributes *magnet_visatt = new G4VisAttributes( G4Colour( 0.75, 0.75, 0.75 ) );
  big48d48Log->SetVisAttributes( magnet_visatt );
  
  if( fDetCon->fTotalAbs ){
    big48d48Log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }


  new G4PVPlacement(bigboxrm, 
  		    G4ThreeVector(-r48d48*sin(f48D48ang), 0.0, r48d48*cos(f48D48ang)),
  		    big48d48Log, "big48d48Physical", worldlog, 0,false,0);

  G4LogicalVolume *bigfieldLog=new G4LogicalVolume(biggap, GetMaterial("Air"),
						   "bigfieldLog", 0, 0, 0);

  //AJRP 9/29/2020: moved SBS sieve plate definition to its own method:



  // Dipole gap shielding for GEnRP

  if( fDetCon->fExpType == G4SBS::kGEnRP ){
    // Defining parametermes for rectangular portion
    G4double h_gapshield = 30.*cm;  // guesstimate
    G4double w_gapshield = 68.5*cm;  // guesstimate
    G4double d_gapshield = 60.*cm;
  
    G4double d1 = r48d48 + f48D48depth/2.0 - d_gapshield/2.0;
    G4double d2 = f48D48width/2.0 - w_gapshield/2.0;
    G4double a_req = atan(d1/d2)*(180./M_PI)*deg + f48D48ang - 90.*deg;
    G4double d3 = (tan(f48D48ang)*d1 - d2)*cos(f48D48ang);
    G4double r_req = d3/sin(a_req);

    G4Box *gapshield = new G4Box("gapshield", w_gapshield/2.0, d_gapshield/2.0, h_gapshield/2.0);

    // Defining parameters for wedge shaped portion
    G4double w_wedge = h_gapshield;
    G4double d_wedge = d_gapshield;
    G4double h_wedge = d_wedge*tan(slot_angle);

    G4Trap *wedge = new G4Trap("wedge", w_wedge, d_wedge, h_wedge, 0.001*mm);

    G4RotationMatrix *wedgerm = new G4RotationMatrix;
    wedgerm->rotateZ(180.*deg);
    wedgerm->rotateX(180.*deg);
  
    G4UnionSolid *dgapshld = new G4UnionSolid("dgapshld", gapshield, wedge, wedgerm,
    					    G4ThreeVector(- w_gapshield/2.0 - h_wedge/4.0, 0.0, 0.0) );
    G4LogicalVolume *dgapshld_log = new G4LogicalVolume(dgapshld, GetMaterial("Lead"), "dgapshld_log" );
    new G4PVPlacement(bigboxrm, G4ThreeVector(-r_req*sin(a_req), 0.0, r_req*cos(a_req)),
		      dgapshld_log, "big48d48Physical", worldlog, 0,false,0);

    G4VisAttributes *Leadcolor = new G4VisAttributes(G4Colour(0.4,0.4,0.4));
    dgapshld_log->SetVisAttributes( Leadcolor );
  }


  // Associate magnetic field with gap

  double sign = 1.0;
  if( fDetCon->fGlobalField->fInverted ) sign = -1.0;

  //Since G4UniformMagField does not inherit from G4SBSMagneticField, we need to apply the scale factor here:
  G4double FieldMag = fFieldStrength * fDetCon->GetFieldScale_SBS();
  
  G4UniformMagField* magField
    = new G4UniformMagField(G4ThreeVector(sign*FieldMag*cos(f48D48ang), 0.0, sign*FieldMag*sin(f48D48ang)));

  G4FieldManager *bigfm = new G4FieldManager(magField);

  //Use BMT equation for mag. field tracking here too, so we can do spin tracking in SBS uniform field
  //in the event that global field is not used
  G4Mag_SpinEqRhs* fBMTequation = new G4Mag_SpinEqRhs(magField);
  G4MagIntegratorStepper *pStepper = new G4ClassicalRK4(fBMTequation,12);
  //G4MagIntegratorStepper *pStepper = new G4DormandPrince745(fBMTequation,12);
  G4ChordFinder *cftemp = new G4ChordFinder(magField, 1.0e-2*mm, pStepper);
  
  bigfm->SetDetectorField(magField);
  bigfm->SetChordFinder(cftemp);

  if( fUseLocalField ){
    bigfieldLog->SetFieldManager(bigfm,true);
  }


  new G4PVPlacement(bigrm, 
  		    G4ThreeVector(-r48d48*sin(f48D48ang), 0.0, r48d48*cos(f48D48ang)),
  		    bigfieldLog, "bigfieldPhysical", worldlog, 0,false,0);


  if( fDetCon->fExpType == G4SBS::kGEp || fDetCon->fExpType == G4SBS::kGEPpositron ){
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

    //These are the "pole shims":

    double slabdepth= 61*cm; //actually slap half-depth: total depth = 122 cm
    G4GenericTrap *leftslab = new G4GenericTrap("leftslab", slabdepth, leftverts );
    G4GenericTrap *rightslab = new G4GenericTrap("rightslab", slabdepth, rightverts );

    G4LogicalVolume *leftslabLog=new G4LogicalVolume(leftslab, GetMaterial("Fer"), "leftslabLog", 0, 0, 0);
    G4LogicalVolume *rightslabLog=new G4LogicalVolume(rightslab, GetMaterial("Fer"), "rightslabLog", 0, 0, 0);

    if( fDetCon->fTotalAbs ){
      leftslabLog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
      rightslabLog->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
    }

    new G4PVPlacement(bigrm, 
    		      G4ThreeVector(-r48d48*sin(f48D48ang), 0.0, r48d48*cos(f48D48ang)),
    		      leftslabLog, "leftslabPhysical", worldlog, 0,false,0);
    new G4PVPlacement(bigrm, 
    		      G4ThreeVector(-r48d48*sin(f48D48ang), 0.0, r48d48*cos(f48D48ang)),
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

  bigfieldLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  backclampLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  //  frontclampLog->SetVisAttributes(G4VisAttributes::GetInvisible());
}

void G4SBSHArmBuilder::MakeSBSFieldClamps( G4LogicalVolume *motherlog ){
  
  if( false ){ //BigBite: Field clamp version 1 (AFAIK) is obsolete. Never build this...

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
   
    if( fDetCon->fTargType == G4SBS::kLH2 || fDetCon->fTargType == G4SBS::kLD2 ){
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
    if( fDetCon->fExpType == G4SBS::kGEp || fDetCon->fExpType == G4SBS::kGMN || fDetCon->fExpType == G4SBS::kGEN || fDetCon->fExpType == G4SBS::kGEnRP ){
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
    if( fDetCon->fExpType == G4SBS::kGEp ){
      new G4PVPlacement(rot, 
			G4ThreeVector(-(r48d48+backclampz)*sin(-f48D48ang)-cos(-f48D48ang)*clampoffset, 0.0, (r48d48+backclampz)*cos(-f48D48ang)-sin(-f48D48ang)*clampoffset),
			backclampLog, "backclampPhysical", motherlog, 0,false,0);
    }

    G4VisAttributes * clampVisAtt
      = new G4VisAttributes(G4Colour(0.8,1.0,0.4));

    frontclampLog->SetVisAttributes(clampVisAtt);
    frontextfaceLog->SetVisAttributes(clampVisAtt);
    backclampLog->SetVisAttributes(clampVisAtt);
  } else if( f48D48_fieldclamp_config != 0 ) { // 2 == GEp, 3 == GMn

    //Any non-zero value causes front clamp to be built: "2" causes both front and rear clamps to be built:
    G4double FrontClamp_width = 62.99*2.54*cm;
    G4double FrontClamp_height = 118.11*2.54*cm;
    G4double FrontClamp_depth = 3.94*2.54*cm;

    G4double FrontClamp_notch_width = 35.83*2.54*cm;
    G4double FrontClamp_notch_height = 27.56*2.54*cm;

    G4Box *FrontClamp_Box = new G4Box( "FrontClamp_Box", FrontClamp_width/2.0, FrontClamp_height/2.0, FrontClamp_depth/2.0 );
    G4Box *FrontClamp_Notch = new G4Box( "FrontClamp_Notch", FrontClamp_notch_width/2.0 + 1.0*cm, FrontClamp_notch_height/2.0, FrontClamp_depth/2.0 + 1.0*cm );

    //Position notch so that the distance from the right edge of the clamp to the left edge of the notch equals notch width:
    //xnotch + wnotch/2 + 1 cm + wclamp/2 = wnotch --> xnotch = wnotch/2- 1 cm - wclamp/2

    //xnotch - wnotch/2 - 1 cm = wclamp/2 - wnotch --> xnotch = wclamp/2 - wnotch/2 + 1 cm
    G4double xnotch = -FrontClamp_notch_width/2.0 + 1.0*cm + FrontClamp_width/2.0;

    G4SubtractionSolid *FrontClamp = new G4SubtractionSolid( "FrontClamp", FrontClamp_Box, FrontClamp_Notch, 0, G4ThreeVector( xnotch, 0.0, 0.0 ) );

    G4LogicalVolume *FrontClamp_log = new G4LogicalVolume( FrontClamp, GetMaterial("Fer"), "FrontClamp_log" );
    fDetCon->InsertAnalyzerVolume( FrontClamp_log-> GetName() );
    if(fDetCon->fTotalAbs) {
      FrontClamp_log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
    }

    G4double FrontClamp_zoffset = 13.40*2.54*cm + FrontClamp_depth/2.0;

    G4double FrontClamp_r = f48D48dist - FrontClamp_zoffset;
    
    G4double FrontClamp_xshift = 14.22*2.54*cm; //x offset of left edge of front clamp relative to left edge of magnet

    G4double FrontClamp_xoffset = -f48D48width/2.0 + FrontClamp_width/2.0 + FrontClamp_xshift; //net offset in position needed to align left edge of front clamp at correct position.

    G4RotationMatrix *clamp_rot = new G4RotationMatrix;

    clamp_rot->rotateY( f48D48ang );

    new G4PVPlacement( clamp_rot, 
		       G4ThreeVector( -FrontClamp_r*sin( f48D48ang ) + FrontClamp_xoffset*cos(f48D48ang), 0.0, FrontClamp_r*cos(f48D48ang) + FrontClamp_xoffset*sin(f48D48ang) ), 
		       FrontClamp_log, "FrontClamp_phys", motherlog, false, 0, false );

    G4VisAttributes * clampVisAtt
      = new G4VisAttributes(G4Colour(0.8,1.0,0.4));

    FrontClamp_log->SetVisAttributes(clampVisAtt);
 
    // (Note jc2): Stickng this if statement here so we can draw the
    // back fieldclamp only for GEp, but not GMn. However, I'm leaving
    // the indentation as is because I don't want git to make me the
    // author of all these changes when all I did is indent things :/
    // Sorry in advance to whoemever is angry at the lack of indentation,
    // blame Linus for making git :D
    if( f48D48_fieldclamp_config == 2 ) { // GEp uses a back clamp too
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
									G4ThreeVector( -RearClamp_width/2.0 + RearClamp_GapX + RearClamp_GapWidth/2.0, 0.0, 0.0 ) );
    
      G4Box *RearClamp_Notch = new G4Box("RearClamp_Notch", RearClamp_NotchWidth/2.0+1.0*cm, RearClamp_NotchHeight/2.0, RearClamp_depth/2.0 + 1.0*cm );
    
      xnotch = -RearClamp_NotchWidth/2.0 + 1.0*cm + RearClamp_width/2.0;

      G4SubtractionSolid *RearClamp = new G4SubtractionSolid( "RearClamp", RearClamp_GapCutout, RearClamp_Notch, 0, 
							      G4ThreeVector( xnotch, 0.0, 0.0 ) );

      G4LogicalVolume *RearClamp_log = new G4LogicalVolume( RearClamp, GetMaterial("Fer"), "RearClamp_log" );

      fDetCon->InsertAnalyzerVolume( RearClamp_log->GetName() );
      
      if(fDetCon->fTotalAbs) {
	RearClamp_log->SetUserLimits( new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
      }

      G4double RearClamp_zoffset = 11.43*2.54*cm + RearClamp_depth/2.0; 
      G4double RearClamp_xoffset = -f48D48width/2.0 + RearClamp_width/2.0; 
      G4double RearClamp_r = f48D48dist + f48D48depth + RearClamp_zoffset;

      G4ThreeVector RearClamp_pos( -RearClamp_r*sin(f48D48ang) + RearClamp_xoffset*cos(f48D48ang), 0.0, RearClamp_r*cos(f48D48ang) + RearClamp_xoffset * sin(f48D48ang) );

      new G4PVPlacement( clamp_rot, RearClamp_pos, RearClamp_log, "RearClamp_phys", motherlog, false, 0, false );

      RearClamp_log->SetVisAttributes(clampVisAtt);
    } // Rear clamp log (GEp)

    //Make lead shielding in clamp:
    // G4double angtrap = 10.0*deg;
    // G4double Trap_DZ = 70.0*cm; //length in z
    // G4double Trap_theta = 0.0; //polar angle between face at -dz/2 and +dz/2
    // G4double Trap_phi = 0.0; //azimuthal angle between face at -dz/2 and +dz/2
    // G4double Trap_H1 = 13.6*cm; //length in y at -dz
    // G4double Trap_BL1 = 16.0*cm; //width at -H1/2 and -dz/2
    // G4double Trap_TL1 = 20.0*cm; //width at +H1/2 and -dz/2
    // G4double Trap_alpha1 = angtrap; //
    // G4double Trap_H2 = 13.6*cm;
    // G4double Trap_BL2 = 16.0*cm;
    // G4double Trap_TL2 = 20.0*cm;
    // G4double Trap_alpha2 = angtrap;

    if( fDetCon->fLeadOption == 1 &&  fDetCon->fExpType == G4SBS::kGEp ){
      //Let us redefine this guy so that the sides make proper angles:
      G4double Trap_DZ = 13.6*cm;
      G4double Trap_Width1 = 15*cm;
      G4double poleshim_angle = atan( (17.75-12.0)/122.0 ); //relative to SBS central axis
      G4double ang_sbs_edge = 16.9*deg - poleshim_angle;
      G4double dist_poleshim = 1.60*m;
      G4double x_poleshim = 12.0*cm;

      G4ThreeVector zaxis_temp( -sin(16.9*deg), 0, cos(16.9*deg) );
      G4ThreeVector yaxis_temp(0,1,0);
      G4ThreeVector xaxis_temp = (yaxis_temp.cross(zaxis_temp)).unit(); 

      G4double zstart_leadinsert = -102.9*cm + dist_poleshim + f48D48depth/2.0 - Trap_DZ/2.0;
      G4double xstart1_leadinsert = x_poleshim + tan(poleshim_angle)*(zstart_leadinsert - dist_poleshim );

      G4double xstart2_leadinsert = xstart1_leadinsert + tan(poleshim_angle)*Trap_DZ;
      G4double xstop1_leadinsert = xstart1_leadinsert + Trap_Width1;
      G4double xstop2_leadinsert = xstop1_leadinsert + tan(16.9*deg)*Trap_DZ;

      G4ThreeVector posrel_leadinsert( 0.25*(xstart1_leadinsert+xstart2_leadinsert+xstop1_leadinsert+xstop2_leadinsert), 0, zstart_leadinsert + 0.5*Trap_DZ );
      G4double Trap_Width2 = xstop2_leadinsert - xstart2_leadinsert;

      G4cout << "Trap_Width2 = " << Trap_Width2/cm << " cm" << G4endl;
      
      G4double Theta_leadinsert = atan( 0.5*(xstart2_leadinsert+xstop2_leadinsert - xstart1_leadinsert - xstop1_leadinsert)/Trap_DZ );
      G4double Phi_leadinsert = 0.0;

      G4double Trap_Dy = 70.0*cm;

      // G4Trap *FrontClampLeadInsert = new G4Trap( "FrontClampLeadInsert", Trap_DZ/2.0, Trap_theta, Trap_phi, Trap_H1/2.0, Trap_BL1/2.0, Trap_TL1/2.0, Trap_alpha1,
      //					       Trap_H2/2.0, Trap_BL2/2.0, Trap_TL2/2.0, Trap_alpha2 );

      G4Trap *FrontClampLeadInsert =  new G4Trap( "FrontClampLeadInsert", Trap_DZ/2.0, Theta_leadinsert, Phi_leadinsert, Trap_Dy/2.0, Trap_Width1/2.0, Trap_Width1/2.0, 0.0, Trap_Dy/2.0, Trap_Width2/2.0, Trap_Width2/2.0, 0.0 );


      G4LogicalVolume *FrontClampLeadInsert_log = new G4LogicalVolume( FrontClampLeadInsert, GetMaterial("Lead"), "FrontClampLeadInsert_log" );
      G4VisAttributes *lead_visatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
      FrontClampLeadInsert_log->SetVisAttributes( lead_visatt );
    
      //G4ThreeVector FrontClampLeadInsert_posrel( +17.5*cm, 0.0, -102.9*cm ); //Position relative to SBS magnet
      //rotation matrix has z axis along +y and x axis along +x, y axis along -z. This is a rotation about the x axis by 90 deg:
      G4RotationMatrix *rot_lead = new G4RotationMatrix;
      rot_lead->rotateY( f48D48ang );
      //rot_lead->rotateZ( 180.0*deg );
      //rot_lead->rotateX( -90.0*deg );
      //rot_lead->rotateX( 90.0*deg );
   
      G4ThreeVector SBS_xaxis( cos( f48D48ang ), 0, sin(f48D48ang ) );
      G4ThreeVector SBS_yaxis(0,1,0);
      G4ThreeVector SBS_zaxis( -sin( f48D48ang ), 0, cos(f48D48ang ) );
    
      G4ThreeVector FrontClampLeadInsert_pos = posrel_leadinsert.x() * SBS_xaxis + posrel_leadinsert.y() * SBS_yaxis + posrel_leadinsert.z() * SBS_zaxis;

      // Add lead inserts IFF lead option is turned on:
      //if( fDetCon->fLeadOption == 1 ){
      new G4PVPlacement( rot_lead, FrontClampLeadInsert_pos, FrontClampLeadInsert_log, "FrontClampLeadInsert_phys", motherlog, false, 0, false );
    
      
      //Add lead bar:
      G4Box *FrontClampLeadBar = new G4Box( "FrontClampLeadBar", 25.0*cm, 10.0*cm, 2.3*cm );
      G4LogicalVolume *FrontClampLeadBar_log = new G4LogicalVolume( FrontClampLeadBar, GetMaterial("Lead"), "FrontClampLeadBar_log" );
      FrontClampLeadBar_log->SetVisAttributes( lead_visatt );
      G4ThreeVector FrontClampLeadBar_posrel( 0.0, FrontClamp_notch_height/2.0 + 10.0*cm, -107.4*cm );
    
      G4ThreeVector FrontClampLeadBar_pos = FrontClampLeadBar_posrel.x()*SBS_xaxis + FrontClampLeadBar_posrel.y() * SBS_yaxis +
	( FrontClampLeadBar_posrel.z() + f48D48dist + 24.0*2.54*cm)*SBS_zaxis;

      //if( fDetCon->fLeadOption == 1 ){
      new G4PVPlacement( clamp_rot, FrontClampLeadBar_pos, FrontClampLeadBar_log, "FrontClampLeadBar_phys", motherlog, false, 0, false );
    }
    
  }
}

void G4SBSHArmBuilder::MakeHCALV2( G4LogicalVolume *motherlog,
    G4double VerticalOffset)
{
  //******************************************************
  //****************         HCAL         ****************
  //****************************************************** 
  // Set to true to enable checking of overlapping volumes at construction
  G4bool checkOverlap = fDetCon->fCheckOverlap;

  // Determine how "realistic" we want to make HCal. Include the various small
  // holes and plates?
  const G4bool const_FrontPinHole = false; // Mounting pin hole on front plate

  // Modify Vahe's implementation so that the scintillators are now staggered
  // Also updated with latest measurements based on what CMU is constructing

  // Dimensions as defined by CAD drawing HCALJ-0002
  // This contains only the "sensitive" parts of the module:
  //   Scintillators + Fe Absorbers + Wavelenght Shifter + Front/Back Plate
  // Excludes:
  //  - Extending cylindrical piece of light guide
  //  - PMT and associated enclosure
  G4double dim_ModuleX    =  152.40*CLHEP::mm;
  G4double dim_ModuleY    =  152.40*CLHEP::mm;
  G4double dim_ModuleZ    = 1555.70*CLHEP::mm;
  G4double dim_ModCanZ    = 1238.00*CLHEP::mm;//948.00*CLHEP::mm;

  // Also define the "CAN" dimensions, since there is still ~1.5 mm of steel
  // in the top and bottom  about ~3 mm of steel at the sides.
  G4double dim_ModCanThickX = 2*1.4478*CLHEP::mm; // Twice as thick in sides
  G4double dim_ModCanThickY =   1.4478*CLHEP::mm; // Top and bottom thickness
  G4double dim_ModCanX0    = dim_ModuleX-dim_ModCanThickX;
  G4double dim_ModCanY0    = dim_ModuleY-dim_ModCanThickY;

  // Dimensions as defined by CAD drawing HCALJ-00011
  G4double dim_ScintX     =   69.35*CLHEP::mm;
  G4double dim_ScintY     =  147.00*CLHEP::mm;
  G4double dim_ScintZ     =    9.90*CLHEP::mm;

  // Dimensions as defined by CAD drawing HCALJ-00010-V2
  // (simplified, does not include the small triangular-cuts on outer corners)
  G4double dim_ThinAbsorbX    =   69.60*CLHEP::mm;
  G4double dim_ThinAbsorbY    =  148.80*CLHEP::mm;
  G4double dim_ThinAbsorbZ    =    6.35*CLHEP::mm;

  // Dimensions as defined by CAD drawing HCALJ-009
  // This is the 1/2 inch thick absorbers (just twice the width of
  // the normal absorbers defined above)
  // For simplicity, these will also serve as the ribs
  G4double dim_AbsorbX =   69.60*CLHEP::mm;
  G4double dim_AbsorbY =  148.80*CLHEP::mm;
  G4double dim_AbsorbZ =   12.70*CLHEP::mm;

  // Dimensions of wavelength shifter (deduced from light guide CAD sheet)
  G4double dim_WaveShiftX     =   5.00*CLHEP::mm;
  G4double dim_WaveShiftY     = 135.02*CLHEP::mm;
  G4double dim_WaveShiftZ     = 948.00*CLHEP::mm;
  G4double dim_WaveShiftRodX  =   6.00*CLHEP::mm;
  G4double dim_WaveShiftRodY  =   5.00*CLHEP::mm;
  G4double dim_WaveShiftRodZ  =   dim_WaveShiftZ;

  // Dimensions based on the real thickness of the plastic spacers used.
  // (This is a spacer that goes in between the last scintillator of a 
  // substack and the rib)
  G4double dim_ShimGapSpacerX =  69.60*CLHEP::mm;
  G4double dim_ShimGapSpacerY = 148.80*CLHEP::mm;
  G4double dim_ShimGapSpacerZ =   4.00*CLHEP::mm;

  // Dimensions of mylar wrapping (from Vahe's HCalo code)
  //G4double dim_MylarThickness = 0.2*CLHEP::mm;
  // Got these from a sample of Mylar Eric had around (it was ~2-2.4 mil)
  G4double dim_MylarThickness = 0.0508*CLHEP::mm;

  // Specify the number of rows and columns for standard HCAL configuration
  G4double num_rows = 24;
  G4double num_cols = 12;

  // Specify center-to-center distance between modules
  // (These numbers are from an email from Gregg)
  G4double dist_ModuleCToCX     =  154.94*CLHEP::mm;
  G4double dist_ModuleCToCY     =  158.75*CLHEP::mm;

  // Specify a minimum distance between objects (an "air" gap)
  G4double dim_AirGap = 0*0.05*CLHEP::mm;

  // Specify number of Scint-Fe pairs per substack
  const int numSubstacks = 3;
  int numPairs[numSubstacks] = {13, 14, 13};

  // Dimensions for Front and Back Plate as defined by CAD drawing HCALJ-00012
  // and HCALJ-00016
  G4double dim_EndPl0X          =  dim_ModuleX;
  G4double dim_EndPl0Y          =  dim_ModuleY;
  G4double dim_EndPl0Z          =    6.35*CLHEP::mm;
  G4double dim_EndPl1X          =  149.35*CLHEP::mm;
  G4double dim_EndPl1Y          =  149.35*CLHEP::mm;
  G4double dim_EndPl1Z          =   12.70*CLHEP::mm;
  G4double dim_EndPlZ           =  dim_EndPl0Z+dim_EndPl1Z;
  G4double dim_FrontPlHoleR     =    9.525/2.*CLHEP::mm;
  G4double dim_FrontPlHoleZ     =   15.875*CLHEP::mm;
  G4double dim_BackPlHoleR      =   53.950/2.*CLHEP::mm;
  G4double dim_BackPlHoleZ      =  dim_EndPlZ;
  G4double dim_BackPlCutX       =   89.00*CLHEP::mm;
  G4double dim_BackPlCutY       =   89.00*CLHEP::mm;
  G4double dim_BackPlCutZ       =   11.20*CLHEP::mm;

  // Dimensions for the light guide. The first part that transforms rectangle
  // to cylindrical extension is represented as a trapezoid for simplicity.
  // (This is similar to how Vahe defined it in his code).
  G4double dim_LightGuide0X0    = dim_WaveShiftX;
  G4double dim_LightGuide0Y0    = dim_WaveShiftY;
  G4double dim_LightGuide0X1    =   20.05*CLHEP::mm;
  G4double dim_LightGuide0Y1    =   43.20*CLHEP::mm;
  G4double dim_LightGuide0Z     =  262.00*CLHEP::mm;
  G4double dim_LightGuide1R     =   47.63/2*CLHEP::mm;
  G4double dim_LightGuide1Z     =  102.22*CLHEP::mm;

  // These numbers come from Vahe for the XP2262 PMT
  G4double dim_PMTCathodeR      =   27.0*CLHEP::mm;
  G4double dim_PMTCathodeZ      =    1.0*CLHEP::mm;

  // Dimensions for the Module Mylar Wrap
  G4double dim_ModuleMylar0X     = dim_EndPl1X;
  G4double dim_ModuleMylar0Y     = dim_EndPl1Y;
  G4double dim_ModuleMylar0Z     = dim_ModCanZ-2*dim_EndPl1Z;
  G4double dim_ModuleMylar1X     = dim_ModuleMylar0X-2*dim_MylarThickness;
  G4double dim_ModuleMylar1Y     = dim_ModuleMylar0Y-2*dim_MylarThickness;
  G4double dim_ModuleMylar1Z     = dim_ModuleMylar0Z;

  // A small imperfection to use when cutting solids. It seems that if the
  // solid and its associated cut share the same surface, they are not
  // drawn properly.
  G4double dim_CutFudgeFactor  = 0.005*CLHEP::mm;

  // Specify the dimensions of each "external" horizontal shim
  G4double dim_ExternalShimX   = 2*927.1*CLHEP::mm;
  G4double dim_ExternalShimY   = 3.51*CLHEP::mm;
  G4double dim_ExternalShimZ   = 882.65*CLHEP::mm;
  G4double dim_ExternalShimTabZ = 44.45*CLHEP::mm;

  // Specify the thickness of the front "mounting" plate that is just in
  // front of each HCal sub sub-assembly (a steel plate)
  // According to Eric Day, this should be about 0.75 inches
  // (and for simplicity, the x,y dimensions will be the whole of HCal)
  G4double dim_HCALFrontPlateZ = 19.05*CLHEP::mm;

  // Specify the dimensions of the HCAL mother box/volume
  G4double dim_HCALX = (num_cols-1)*dist_ModuleCToCX + dim_ModuleX;
  G4double dim_HCALY = (num_rows-1)*dist_ModuleCToCY + dim_ModuleY;
  G4double dim_HCALZ = dim_ModuleZ+dim_HCALFrontPlateZ;

  // Specify the distance from entrance to HCAL to center of target
  G4double dist_HCalRadius = fHCALdist + dim_HCALZ/2.0;
  G4double dist_HCALX = -dist_HCalRadius*sin(f48D48ang+fHCALangular_offset);
  G4double dist_HCALY = VerticalOffset;
  G4double dist_HCALZ = dist_HCalRadius*cos(f48D48ang+fHCALangular_offset);

  // Specify the rotation matrix for this HCAL
  G4RotationMatrix *rot_HCAL= new G4RotationMatrix;
  rot_HCAL->rotateY(f48D48ang+fHCALangular_offset);

  // Define the solid for the module container and "CAN"
  G4Box *sol_Module = new G4Box("sol_Module",
      dim_ModuleX/2., dim_ModuleY/2., dim_ModuleZ/2.);
  G4Box *sol_ModCan0 = new G4Box("sol_ModCan0",
      dim_ModuleX/2., dim_ModuleY/2., dim_ModCanZ/2.);
  G4Box *sol_ModCan1 = new G4Box("sol_ModCan0",
      dim_ModCanX0/2., dim_ModCanY0/2., dim_ModCanZ/2.);
  G4SubtractionSolid *sol_ModCan = new G4SubtractionSolid("sol_ModCan",
      sol_ModCan0, sol_ModCan1,0, G4ThreeVector(0,0,0));

  // Define solids for scintillator, absorber and shim gap spacer
  G4Box *sol_Scint  = new G4Box("sol_Scint",
      dim_ScintX/2., dim_ScintY/2., dim_ScintZ/2.);
  G4Box *sol_ThinAbsorb = new G4Box("sol_ThinAbsorb",
      dim_ThinAbsorbX/2., dim_ThinAbsorbY/2., dim_ThinAbsorbZ/2.);
  G4Box *sol_Absorb = new G4Box("sol_Absorb",
      dim_AbsorbX/2., dim_AbsorbY/2., dim_AbsorbZ/2.);
  G4Box *sol_ShimGapSpacer = new G4Box("sol_ShimGapSpacer",
      dim_ShimGapSpacerX/2., dim_ShimGapSpacerY/2., dim_ShimGapSpacerZ/2.);

  // Define solid for WaveShift and the rods at the edge
  G4Box *sol_WaveShift = new G4Box("sol_WaveShift",
      dim_WaveShiftX/2., dim_WaveShiftY/2., dim_WaveShiftZ/2.);
  G4Box *sol_WaveShiftRod = new G4Box("sol_WaveShiftRod",
      dim_WaveShiftRodX/2., dim_WaveShiftRodY/2., dim_WaveShiftRodZ/2.);

  // The Front plate is a "union" of two plates and a hole cut out
  // (and the back plate has an additional square cut)
  G4Box *sol_EndPl0 = new G4Box("sol_EndPl0",
      dim_EndPl0X/2.,dim_EndPl0Y/2.,dim_EndPl0Z/2.);
  G4Box *sol_EndPl1 = new G4Box("sol_EndPl1",
      dim_EndPl1X/2.,dim_EndPl1Y/2.,dim_EndPl1Z/2.);
  // The front and end plates with "trimmed" region
  G4UnionSolid *sol_FrontPl0 = new G4UnionSolid(
      (const_FrontPinHole?"sol_FrontPl0":"sol_FrontPl"),
      sol_EndPl0,sol_EndPl1,0,G4ThreeVector(0.,0.,
         (dim_EndPl0Z+dim_EndPl1Z)/2.));
  G4UnionSolid *sol_BackPl0 = new G4UnionSolid("sol_BackPl0",
      sol_EndPl0,sol_EndPl1,0,G4ThreeVector(0.,0.,
        -(dim_EndPl0Z+dim_EndPl1Z)/2.));
  // The front hole (cuts 10um bigger to avoid vizualization errors
  // with shared surfaces)
  G4Tubs *sol_FrontPlHole = new G4Tubs("sol_FrontPlHole",0.0,dim_FrontPlHoleR,
      (dim_FrontPlHoleZ+dim_CutFudgeFactor)/2.,0.0*CLHEP::deg,360.0*CLHEP::deg);
  // The back hole and cut
  G4Tubs *sol_BackPlHole = new G4Tubs("sol_BackPlHole",0.0,dim_BackPlHoleR,
      (dim_BackPlHoleZ+dim_CutFudgeFactor)/2.,0.0*CLHEP::deg,360.0*CLHEP::deg);
  G4Box *sol_BackPlCut = new G4Box("sol_BackPlCut",
      dim_BackPlCutX/2.,dim_BackPlCutY/2.,(dim_BackPlCutZ+dim_CutFudgeFactor)/2.);
  G4UnionSolid *sol_BackPlHoleCut  = new G4UnionSolid("sol_BackPlHoleCut",
      sol_BackPlHole,sol_BackPlCut,0,
      G4ThreeVector(0.0,0.0, (dim_BackPlHoleZ-dim_BackPlCutZ)/2.));
  // Cut the pieces on front and back plates
  G4VSolid *sol_FrontPl = 0;
  if(const_FrontPinHole ) {
    sol_FrontPl = new G4SubtractionSolid("sol_FrontPl",
        sol_FrontPl0, sol_FrontPlHole,0,
        G4ThreeVector(0.,0., -(dim_EndPl0Z-dim_FrontPlHoleZ)/2.));
  } else {
    sol_FrontPl = sol_FrontPl0;
  }
  G4SubtractionSolid *sol_BackPl = new G4SubtractionSolid("sol_BackPl",
      sol_BackPl0,sol_BackPlHoleCut,0,
      G4ThreeVector(0.,0., (dim_EndPl0Z-dim_BackPlHoleZ)/2.));

  // Light Guide
  G4Trd *sol_LightGuide0 = new G4Trd("sol_LightGuide0",
      dim_LightGuide0X0/2.,dim_LightGuide0X1/2.,
      dim_LightGuide0Y0/2.,dim_LightGuide0Y1/2.,
      dim_LightGuide0Z/2.);
  G4Tubs *sol_LightGuide1 = new G4Tubs("sol_LightGuide1",0.0,dim_LightGuide1R,
      dim_LightGuide1Z/2.,0.0*CLHEP::deg,360.0*CLHEP::deg);
  G4UnionSolid *sol_LightGuide = new G4UnionSolid("sol_LightGuide",
      sol_LightGuide0,sol_LightGuide1,0,G4ThreeVector(0.,0.,
        (dim_LightGuide0Z+dim_LightGuide1Z)/2.));

  // PMT (cathode)
  G4Tubs *sol_PMTCathode = new G4Tubs("sol_PMTCathode",0.0,dim_PMTCathodeR,
      dim_PMTCathodeZ/2.,0.0*CLHEP::deg,360.0*CLHEP::deg);

  // Mylar wrapping around the module. The "ribbon" style wrapping is
  // not implemented for simplicity. Apparently implementing it really
  // takes a toll on Geant4. Instead, the proper reflection of the various
  // components will be set later, in liue of adding so much mylar.
  G4Box *sol_ModuleMylar0 = new G4Box("sol_ModuleMylar0",
      dim_ModuleMylar0X/2.,dim_ModuleMylar0Y/2.,dim_ModuleMylar0Z/2.);
  G4Box *sol_ModuleMylar1 = new G4Box("sol_ModuleMylar1",
      dim_ModuleMylar1X/2.,dim_ModuleMylar1Y/2.,
      (dim_ModuleMylar1Z+dim_CutFudgeFactor)/2.);
  G4SubtractionSolid *sol_ModuleMylar = new G4SubtractionSolid("sol_ModuleMylar",
      sol_ModuleMylar0,sol_ModuleMylar1,0,G4ThreeVector(0.,0.,0.));

  // Define the external shims (horizontal direction) solid and volume
  G4Box *sol_ExternalShim = new G4Box("sol_ExternalShim",
      dim_ExternalShimX/2.,dim_ExternalShimY/2.,dim_ExternalShimZ/2.);
  G4LogicalVolume *log_ExternalShim = new G4LogicalVolume( sol_ExternalShim,
      GetMaterial("Steel"), "log_ExternalShim");

  // Define logical volumes for Sensitive-Modules, Scintillators, Absorbers,
  // and Front/Back plates
  G4LogicalVolume *log_Module = new G4LogicalVolume( sol_Module,
      GetMaterial("Special_Air"), "log_Module" );

  // Module Can
  G4LogicalVolume *log_ModCan = new G4LogicalVolume( sol_ModCan,
      GetMaterial("Steel"), "log_ModCan");
  // Scintillator
  G4LogicalVolume *log_Scint = new G4LogicalVolume( sol_Scint,
      GetMaterial("EJ232"), "log_Scint");

  // Absorbers and wrap
  G4LogicalVolume *log_ThinAbsorb = new G4LogicalVolume( sol_ThinAbsorb,
      GetMaterial("Iron"), "log_ThinAbsorb");
  G4LogicalVolume *log_Absorb = new G4LogicalVolume( sol_Absorb,
      GetMaterial("Iron"), "log_Absorb");

  // The wrapped spacer
  G4LogicalVolume *log_ShimGapSpacer = new G4LogicalVolume( sol_ShimGapSpacer,
      GetMaterial("Air"), "log_ShimGapSpacer");

  // The other volumes
  G4LogicalVolume *log_WaveShift = new G4LogicalVolume( sol_WaveShift,
      GetMaterial("BC484"), "log_WaveShift");
  G4LogicalVolume *log_WaveShiftRod = new G4LogicalVolume( sol_WaveShiftRod,
      GetMaterial("Steel"), "log_WaveShiftRod");
  G4LogicalVolume *log_FrontPl = new G4LogicalVolume( sol_FrontPl,
      GetMaterial("Steel"), "log_FrontPl");
  G4LogicalVolume *log_BackPl = new G4LogicalVolume( sol_BackPl,
      GetMaterial("Steel"), "log_BackPl");
  G4LogicalVolume *log_LightGuide = new G4LogicalVolume( sol_LightGuide,
      //GetMaterial("Acrylic"), "log_LightGuide");
      GetMaterial("BC484"), "log_LightGuide");
  G4LogicalVolume *log_PMTCathode = new G4LogicalVolume( sol_PMTCathode,
      GetMaterial("Glass_HC"), "log_PMTCathode");
  G4LogicalVolume *log_ModuleMylar = new G4LogicalVolume( sol_ModuleMylar,
      GetMaterial("Mylar"), "log_ModuleMylar");

  // Define the solid and logical volume for the FrontPlate steel mount
  G4Box *sol_HCALFrontPlate = new G4Box("sol_HCALFrontPlate",
      dim_HCALX/2., dim_HCALY/2., dim_HCALFrontPlateZ/2.);
  G4LogicalVolume *log_HCALFrontPlate = new G4LogicalVolume( sol_HCALFrontPlate,
      GetMaterial("Steel"), "log_HCALFrontPlate" );

  // Define solid and logical volume for mother HCAL box/volume
  // (make it slightly larger than needed so that setting the stepping limits
  // works)
  G4Box *sol_HCAL = new G4Box("sol_HCAL",
      dim_HCALX/2.+0.01*mm, dim_HCALY/2.+0.01*mm, dim_HCALZ/2.+0.01*mm);
  G4LogicalVolume *log_HCAL = new G4LogicalVolume( sol_HCAL,
      GetMaterial("Air"), "log_HCAL" );

  // Position information for each element in the module
  G4double posZ = -dim_ModuleZ/2.;

  // Position the CAN inside the module
  posZ += dim_EndPl0Z + dim_ModCanZ/2.;
  new G4PVPlacement(0,G4ThreeVector(0,0,posZ),log_ModCan,"pModCan",
      log_Module,false,0,checkOverlap);
  posZ -= dim_EndPl0Z + dim_ModCanZ/2.;

  // Then position the front plate and back plate
  posZ += dim_EndPl0Z/2.;
  G4double pos_EndOfFrontPlZ = posZ + dim_EndPl0Z/2. + dim_EndPl1Z;
  new G4PVPlacement(0,G4ThreeVector(0,0,posZ),log_FrontPl,"pFrontPl",
    log_Module,false,0,checkOverlap);
  posZ += dim_EndPl0Z/2.+dim_ModCanZ+dim_EndPl0Z/2.;
  new G4PVPlacement(0,G4ThreeVector(0,0,posZ),log_BackPl,"pBackPl",
    log_Module,false,0,checkOverlap);
  posZ = pos_EndOfFrontPlZ;

  // Position the mylar wrapping around the module
  posZ += dim_ModuleMylar0Z/2.;
  new G4PVPlacement(0,G4ThreeVector(0,0,posZ),
      log_ModuleMylar,"pModuleMylar",log_Module,false,0.,checkOverlap);
  posZ = pos_EndOfFrontPlZ;

  // Position the wavelength shifter in the middle
  posZ += dim_WaveShiftZ/2;
  new G4PVPlacement(0, G4ThreeVector(0.,0.,posZ),
      log_WaveShift,"pWaveShift",log_Module, false, 0, checkOverlap);
  new G4PVPlacement(0, G4ThreeVector(0.,(dim_WaveShiftY+dim_WaveShiftRodY)/2.,
        posZ),log_WaveShiftRod,"pWaveShiftRod1",
      log_Module, false, 0, checkOverlap);
  new G4PVPlacement(0, G4ThreeVector(0.,-(dim_WaveShiftY+dim_WaveShiftRodY)/2.,
        posZ),log_WaveShiftRod,"pWaveShiftRod1",
      log_Module, false, 0, checkOverlap);
  posZ += dim_WaveShiftZ/2;

  // Position the Light Guides at the end of the WaveShift
  posZ += dim_LightGuide0Z/2;
  new G4PVPlacement(0,G4ThreeVector(0.,0.,posZ), log_LightGuide,"pLightGuide",
      log_Module,false,0, checkOverlap);
  posZ += dim_LightGuide0Z/2 + dim_LightGuide1Z;

  // Record the relative Position the Cathode at the end of the LightGuide
  posZ += dim_PMTCathodeZ/2;
  G4double posPMTCathodeZ = posZ;
  new G4PVPlacement(0,G4ThreeVector(0.,0.,posPMTCathodeZ),
      log_PMTCathode,"pPMTCathode", log_Module, false, 0, checkOverlap);

  // Reset the initial Z position back to the front of the first plate
  posZ = pos_EndOfFrontPlZ;

  // A copy number for each element placed in the module
  int copyNo = 0;

  // Position information for each "side" of the stack
  G4double posLeftZ = posZ;
  G4double posRightZ = posZ;

  G4double posScintX  = dim_ScintX/2.0  + dim_WaveShiftRodX/2. + dim_AirGap;
  G4double posScintY  = 0.0;
  G4double posAbsorbX = dim_AbsorbX/2.0 + dim_WaveShiftRodX/2. + dim_AirGap +
    dim_MylarThickness;
  G4double posAbsorbY = 0.0;

  G4ThreeVector posLeft(-posAbsorbX, posAbsorbY,posLeftZ);
  G4ThreeVector posRight(posAbsorbX, posAbsorbY,posRightZ);

  // Rotation matrix for the wrapped scintillators on both sides
  G4RotationMatrix rotLeft;
  G4RotationMatrix rotRight;
  rotLeft.rotateZ(180*CLHEP::deg);
  
  // Place stack elements in the module
  for(int sub = 0; sub < numSubstacks; sub++) {

    // Place all Scintillator-Fe pairs
    for(int pair = 0; pair < numPairs[sub]; pair++) {
      // The first pair usually doesn't start with an absorber
      // except for the first substack
      if(sub==0 && pair==0) {
        posLeftZ  += dim_ThinAbsorbZ/2.;
        posRightZ += dim_AbsorbZ/2.;
        posLeft.set(-posAbsorbX,posAbsorbY,posLeftZ);
        posRight.set(posAbsorbX,posAbsorbY,posRightZ);
        new G4PVPlacement(0, posLeft, log_ThinAbsorb, "pAbsorb",
            log_Module, false, copyNo++, checkOverlap);
        new G4PVPlacement(0, posRight, log_Absorb, "pAbsorb",
            log_Module, false, copyNo++, checkOverlap);
        posRightZ += dim_AbsorbZ/2. + dim_ThinAbsorbZ/2.;
        posRight.set(posAbsorbX,posAbsorbY,posRightZ);
        new G4PVPlacement(0, posRight, log_ThinAbsorb, "pAbsorb",
            log_Module, false, copyNo++, checkOverlap);
        posLeftZ  += dim_ThinAbsorbZ/2.;
        posRightZ += dim_ThinAbsorbZ/2.;
      } else if(pair>0) {
        // Position an absorber in front of every other pair after the first
        posLeftZ  += dim_AbsorbZ/2. + dim_AirGap;
        posRightZ += dim_AbsorbZ/2. + dim_AirGap;
        posLeft.set(-posAbsorbX,posAbsorbY,posLeftZ);
        posRight.set(posAbsorbX,posAbsorbY,posRightZ);
        new G4PVPlacement(0, posLeft, log_Absorb, "pAbsorb",
            log_Module,false, copyNo++, checkOverlap);
        new G4PVPlacement(0, posRight, log_Absorb, "pAbsorb",
            log_Module,false, copyNo++, checkOverlap);
        posLeftZ  += dim_AbsorbZ/2.;
        posRightZ += dim_AbsorbZ/2.;
      }

      // Now set the scintillator
      posLeftZ  += dim_ScintZ/2. + dim_MylarThickness/2. + dim_AirGap;
      posRightZ += dim_ScintZ/2. + dim_MylarThickness/2. + dim_AirGap;
      posLeft.set(-posScintX,posScintY,posLeftZ);
      posRight.set(posScintX,posScintY,posRightZ);
      new G4PVPlacement(G4Transform3D(rotLeft, posLeft), log_Scint,
          "pScint", log_Module,false, copyNo++, checkOverlap);
      new G4PVPlacement(G4Transform3D(rotRight, posRight), log_Scint,
          "pScint", log_Module,false, copyNo++, checkOverlap);
      posLeftZ  += dim_ScintZ/2.+dim_MylarThickness/2.;
      posRightZ += dim_ScintZ/2.+dim_MylarThickness/2.;

    }
    // A spacer (or shim gap) is at the end of each substack but
    // before the last rib.
    posLeftZ  += dim_ShimGapSpacerZ/2. + dim_AirGap;
    posRightZ += dim_ShimGapSpacerZ/2. + dim_AirGap;
    posLeft.set(-posAbsorbX,posAbsorbY,posLeftZ);
    posRight.set(posAbsorbX,posAbsorbY,posRightZ);
    new G4PVPlacement(0, posLeft, log_ShimGapSpacer, "pShimGapSpacer",
        log_Module, false, copyNo++, checkOverlap);
    new G4PVPlacement(0, posRight, log_ShimGapSpacer, "pShimGapSaper",
        log_Module, false, copyNo++, checkOverlap);
    posLeftZ  += dim_ShimGapSpacerZ/2.;
    posRightZ += dim_ShimGapSpacerZ/2.;

    // Place a "rib" at the end of each substack
    posLeftZ  += dim_AbsorbZ/2. + dim_AirGap;
    posRightZ += dim_AbsorbZ/2. + dim_AirGap;
    posLeft.set(-posAbsorbX,posAbsorbY,posLeftZ);
    posRight.set(posAbsorbX,posAbsorbY,posRightZ);
    new G4PVPlacement(0, posLeft, log_Absorb, "pAbsorb", log_Module,
        false, copyNo++, checkOverlap);
    new G4PVPlacement(0, posRight, log_Absorb, "pAbsorb", log_Module,
        false, copyNo++, checkOverlap);
    posLeftZ  += dim_AbsorbZ/2.;
    posRightZ += dim_AbsorbZ/2.;
  }

  // Now register the Sensitive Detectors
  G4SDManager *sdman = fDetCon->fSDman;

  // Scintillators in HCAL
  G4String HCalScintSDName = "Harm/HCalScint";
  G4String HCalScintCollName = "HCalScintHitsCollection";
  G4SBSCalSD *HCalScintSD = NULL;

  if( !((G4SBSCalSD*) sdman->FindSensitiveDetector(HCalScintSDName)) ){
    G4cout << "Adding HCal Scintillator Sensitive Detector to SDman..." << G4endl;
    HCalScintSD = new G4SBSCalSD( HCalScintSDName, HCalScintCollName );
    sdman->AddNewDetector(HCalScintSD);
    (fDetCon->SDlist).insert(HCalScintSDName);
    fDetCon->SDtype[HCalScintSDName] = G4SBS::kCAL;
    (HCalScintSD->detmap).depth = 1;

    G4int default_ntbins = 25;

    //This will be overridden if the command /g4sbs/threshold has been invoked for this detector:
    fDetCon->SetThresholdTimeWindowAndNTimeBins( HCalScintSDName, 0.0*MeV, 500.0*ns, 25 );
  }
  log_Scint->SetSensitiveDetector(HCalScintSD);
  fDetCon->InsertSDboundaryVolume( log_HCAL->GetName(), HCalScintSDName );

  // PMTs in HCAL (assigned to ECalSD which detects optical photons)
  G4String HCalSDName = "Harm/HCal";
  G4String HCalCollName = "HCalHitsCollection";
  G4SBSECalSD *HCalSD = NULL;

  if( !((G4SBSECalSD*) sdman->FindSensitiveDetector(HCalSDName)) ){
    G4cout << "Adding HCal PMT Sensitive Detector to SDman..." << G4endl;
    HCalSD = new G4SBSECalSD( HCalSDName, HCalCollName );
    sdman->AddNewDetector(HCalSD);
    (fDetCon->SDlist).insert(HCalSDName);
    fDetCon->SDtype[HCalSDName] = G4SBS::kECAL;
    (HCalSD->detmap).depth = 1;

    // *****
    fDetCon->SetThresholdTimeWindowAndNTimeBins( HCalSDName, 0.0*MeV, 250.0*ns, 25 );
    // *****
  }
  log_PMTCathode->SetSensitiveDetector(HCalSD);

  fDetCon->InsertSDboundaryVolume( log_HCAL->GetName(), HCalSDName );
  
  // Set step limit to all of HCAL?
  if( (fDetCon->StepLimiterList).find( log_HCAL->GetName() ) != (fDetCon->StepLimiterList).end() ){
    G4UserLimits *limits = new G4UserLimits(0.0,0.0,0.0,DBL_MAX, DBL_MAX);
    //log_HCAL->SetUserLimits(new G4UserLimits(0.0,0.0,0.0,DBL_MAX, DBL_MAX) );
    log_HCAL->SetUserLimits(limits);
    //log_Scint->SetUserLimits(limits);
    //log_HCALFrontPlate->SetUserLimits(limits);
    G4cout << "Making " << log_HCAL->GetName() << " a total absorber" << G4endl;
    G4String sdname = "Harm/HCALbox";
    G4String collname = "HCAL_boxHitsCollection";
    G4SBSCalSD *HCALboxSD = NULL;
    if( !((G4SBSCalSD*) sdman->FindSensitiveDetector(sdname)) ){
      G4cout << "Adding HCALbox sensitive detector to SDman..." << G4endl;
      HCALboxSD = new G4SBSCalSD( sdname, collname );
      sdman->AddNewDetector( HCALboxSD );
      (fDetCon->SDlist).insert( sdname );
      fDetCon->SDtype[sdname] = G4SBS::kCAL;
      (HCALboxSD->detmap).depth = 0;
      log_HCAL->SetSensitiveDetector( HCALboxSD );
    }
  }

  // Place the HCALFrontPlate steel mount
  new G4PVPlacement(0, G4ThreeVector(0,0,-(dim_HCALZ-dim_HCALFrontPlateZ)/2.),
      log_HCALFrontPlate,"log_HCALFrontPlate",log_HCAL,false,0,checkOverlap);

  ofstream mapfile("database/HCAL_map.txt");

  TString currentline;
  currentline.Form("# %15s, %15s, %15s, %18s, %18s",
		   "Cell", "Row", "Column", "Xcenter (cm)", "Ycenter (cm)" );

  mapfile << currentline << endl;
  
  // Set the initial vertical position for a module as the top of HCAL
  // Initial horizontal position would be on the left
  G4double posModX;
  G4double posModY = dim_HCALY/2. - dim_ModuleY/2.;
  G4double posModZ = dim_HCALFrontPlateZ/2.;
  G4double posExternalShimZ = -dim_HCALZ/2. + dim_HCALFrontPlateZ + dim_ExternalShimZ/2.;
  copyNo = 0;
  // Construct physical volumes for each of the sensitive modules
  for(int row = 0; row <  num_rows; row++ ) {
    posModX = dim_HCALX/2. - dim_ModuleX/2.;
    for(int col = 0; col <  num_cols; col++ ) {
      // Place module inside HCAL
      new G4PVPlacement(0, G4ThreeVector(posModX,posModY,posModZ),
          log_Module, "log_Module", log_HCAL, false, copyNo, checkOverlap);

      // Configure Sensitive Detector Map
      (HCalSD->detmap).Row[copyNo] = row;
      (HCalSD->detmap).Col[copyNo] = col;
      (HCalSD->detmap).LocalCoord[copyNo] =
        G4ThreeVector(posModX,posModY,posPMTCathodeZ);

      (HCalScintSD->detmap).Row[copyNo] = row;
      (HCalScintSD->detmap).Col[copyNo] = col;
      (HCalScintSD->detmap).LocalCoord[copyNo] =
        G4ThreeVector(posModX,posModY,posModZ);

      currentline.Form("  %15d, %15d, %15d, %18.3f, %18.3f",
		       copyNo, row, col, posModX/cm, posModY/cm );
      mapfile << currentline << endl;
      
      // Increment horizontal position for next module
      posModX -= dist_ModuleCToCX;

      // Increment the copyNo
      copyNo++;
    }
    // Increment vertical position for next module. For most instances
    // we'll need a shim spacer before doing so
    posModY -= dist_ModuleCToCY/2.;
    if(row < num_rows -1)
      new G4PVPlacement(0, G4ThreeVector(0,posModY,posExternalShimZ),
          log_ExternalShim, "log_ExternalShim", log_HCAL, false, copyNo,
          checkOverlap);

    posModY -= dist_ModuleCToCY/2.;
  }

  mapfile.close();

  G4ThreeVector HCAL_zaxis( dist_HCALX, 0.0, dist_HCALZ );
  G4ThreeVector HCAL_yaxis( 0.0,        1.0, 0.0        );

  HCAL_zaxis = HCAL_zaxis.unit();
  G4ThreeVector HCAL_xaxis = HCAL_yaxis.cross( HCAL_zaxis ).unit();

  G4ThreeVector HCAL_pos =
    HCAL_zaxis * dist_HCalRadius +
    HCAL_xaxis * fHCALhorizontal_offset +
    HCAL_yaxis * fHCALvertical_offset; 
  
  // Lastly, place the HCAL volume
  new G4PVPlacement(rot_HCAL, HCAL_pos,
      log_HCAL, "HCal Mother", motherlog, false, 0, checkOverlap);

  // Apply the step limiter, if requested
  if( (fDetCon->StepLimiterList).find( HCalScintSDName ) != (fDetCon->StepLimiterList).end() ){
    log_Scint->SetUserLimits(new G4UserLimits(0.0,0.0,0.0,DBL_MAX, DBL_MAX) );
  }

  // Register optical properties of materials. For this, I think it will be
  // sufficient to specify optical surfaces between scintillators-iron and
  // scintillator-waveshift, waveshift-lightguide, lightguide-pmtcathode

  // Define optical skin of light guide
  new G4LogicalSkinSurface("skin_LightGuide",log_LightGuide,
      GetOpticalSurface("osWLSToAir"));
  // Make the mylar pieces reflective
  new G4LogicalSkinSurface("skin_ThinAbsorb",log_ThinAbsorb,
      GetOpticalSurface("Foil") ); 
  new G4LogicalSkinSurface("skin_Absorb",log_Absorb,
      GetOpticalSurface("Foil") ); 
  new G4LogicalSkinSurface("skin_ModuleMylar",log_ModuleMylar,
      GetOpticalSurface("Foil") ); 
  new G4LogicalSkinSurface("skin_ShimGapSpacer",log_ShimGapSpacer,
      GetOpticalSurface("Foil") ); 
  // Define optical surface of a reflector
  /*
  G4OpticalSurface *sur_Reflector = new G4OpticalSurface("sur_Reflector");
  sur_Reflector->SetType(dielectric_metal);
  sur_Reflector->SetFinish(polished);
  */

  // Scint-Absorber border
  /*
  G4LogicalBorderSurface *sur_ScintAbsorb = new G4LogicalBorderSurface(
      "surf_ScintAbsorb",log_Scint,log_Absorb,sur_Reflector);
  G4LogicalBorderSurface *sur_ScintThinAbsorb = new G4LogicalBorderSurface(
      "surf_ScintThinAbsorb",log_Scint,log_ThinAbsorb,sur_Reflector);
      */

  //--------------------------------------------------------------------------
  // Set visualization properties
  double trans = 0.35; // Translucency of light guide + scintillators

  // Module
  //G4VisAttributes *vis_Module  = new G4VisAttributes(G4Colour::Blue());
  //vis_Module ->SetForceWireframe(true);
  //log_Module->SetVisAttributes(vis_Module);
  log_Module->SetVisAttributes(G4VisAttributes::GetInvisible());

  // Hide the Can
  //G4VisAttributes *vis_ModCan = new G4VisAttributes(G4Colour(0.78,0.92,0.27));
  //log_ModCan->SetVisAttributes(vis_ModCan);
  log_ModCan->SetVisAttributes(G4VisAttributes::GetInvisible());

  // Module Mylar (invisible)
  log_ModuleMylar->SetVisAttributes(G4VisAttributes::GetInvisible());

  // Absorber
  G4VisAttributes *vis_Absorb = new G4VisAttributes(G4Colour(0.83,0.84,0.85));
  log_Absorb->SetVisAttributes(vis_Absorb);
  log_ThinAbsorb->SetVisAttributes(vis_Absorb);

  // Scintillator and Mylar
  G4VisAttributes *vis_Scint = new G4VisAttributes(G4Colour(0.0,1.0,1.0,trans));
  log_Scint->SetVisAttributes(vis_Scint);

  // Wavelength Shifter & Rod
  G4VisAttributes *vis_WaveShift = new G4VisAttributes(G4Colour(0.54,0.53,0.79,trans));
  log_WaveShift->SetVisAttributes(vis_WaveShift);
  G4VisAttributes *vis_WaveShiftRod = new G4VisAttributes(G4Colour(0.3,0.3,0.3));
  log_WaveShiftRod->SetVisAttributes(vis_WaveShiftRod);

  // LightGuide
  G4VisAttributes *vis_LightGuide = new G4VisAttributes(G4Colour(1.0,1.0,1.0,trans));
  log_LightGuide->SetVisAttributes(vis_LightGuide);

  // PMTCathode
  G4VisAttributes *vis_PMTCathode = new G4VisAttributes(G4Colour(1.0,.41,71));
  log_PMTCathode->SetVisAttributes(vis_PMTCathode);

  // FrontPl & BackPl
  G4VisAttributes *vis_FrontPl = new G4VisAttributes(G4Colour(0.3,0.3,0.3));
  log_FrontPl->SetVisAttributes(vis_FrontPl);
  G4VisAttributes *vis_BackPl = new G4VisAttributes(G4Colour(0.3,0.3,0.3));
  log_BackPl->SetVisAttributes(vis_BackPl);

  // Shim Gap Spacer
  G4VisAttributes *vis_ShimGapSpacer = new G4VisAttributes(
      G4Colour(1.0,1.0,0.0));
  log_ShimGapSpacer->SetVisAttributes(vis_ShimGapSpacer);

  // FrontPlate (steel mounting plate)
  G4VisAttributes *vis_HCALFrontPlate = new G4VisAttributes(
      G4Colour(1.0,0.0,0.0));
  log_HCALFrontPlate->SetVisAttributes(vis_HCALFrontPlate);

  // HCAL enclosure should be invisible
  log_HCAL->SetVisAttributes(G4VisAttributes::GetInvisible());
}

void G4SBSHArmBuilder::MakeHCAL( G4LogicalVolume *motherlog, G4double VerticalOffset=0.0*cm ){


  //******************************************************
  //****************         HCAL         ****************
  //****************************************************** 

  //Code adopted from Vahe, specifically HCalo.cc && HCaloMaterials.cc

  /////////////// NOTE: THIS CODE IS NO LONGER USED/OBSOLETE AND WON'T BE MAINTAINED (SEE MAKEHCALV2) //////////////

  G4RotationMatrix *mRotateZ = new G4RotationMatrix;
  mRotateZ->rotateZ( 90 *degree );
  
  double AlFoilThick    = 0.02*cm;
  double IronPlThick    = 1.27*cm;
  double ScinPlThick    = 1.0*cm;
  double PlateX         = 14.6*cm; //difference of 1 mm: presumably to make room for WLS? 14.6+14.5 = 29.1 cm. 
  double PlateY         = 14.5*cm;
  double TotalPlatesL   = 92.0*cm;    
  double ModuleL        = 111.0*cm; //total length of module:
  double ModuleX        = 15.3*cm; //transverse dimensions of module:
  double ModuleY        = 15.3*cm; //transverse dimensions of module:     
  double LightGuideX    = 14.2*cm;
  double LightGuideY    = 0.5*cm;
  double LightGuideZ    = 92.0*cm;
  double ScinToLgGap    = 0.1*cm;
  double ContainerThick = 0.3*cm;
  int NRows             = 24;
  int NColumns          = 12;
  int NumberOfLayers    = 40;
  G4double PlateGaps    = (TotalPlatesL-(NumberOfLayers*(IronPlThick+ScinPlThick)))/(2*NumberOfLayers - 1); //Gap between each Fe/Scint plate

  //  G4double CaloX = ModuleX * NRows * 1.001;
  //G4double CaloY = ModuleY * NColumns * 1.001;

  //Interchange Y <--> X to avoid need for confusing 90-degree rotation!
  G4double CaloX = ModuleX * NColumns * 1.001; 
  G4double CaloY = ModuleY * NRows * 1.001;
  
  //G4double CaloL = ModuleL;
  //G4double CaloL = (LightGuideZ + 2.0 * LightGuideX + 0.5*cm) * 1.001; //this isn't actually used!
  G4double ModuleLtotal = LightGuideZ + 2.0 * LightGuideX + ContainerThick;
  G4double CaloL = (ModuleLtotal + 0.5*cm)* 1.001; //Add PMT photocathode thickness:

  //We want hcalr to be the distance from the origin to the surface of HCAL:
  //double hcaldepth  = 101.0*cm;
  G4double hcaldepth = ModuleLtotal + 0.5*cm;
  double hcalr = fHCALdist + hcaldepth/2.0;
  
  G4double PlateXHalf = PlateX/2.0 - ScinToLgGap - LightGuideY/2.0; // = 7.3 cm - 0.1 cm - 0.25 cm = 6.95 cm
  
  //G4Box *solModule = new G4Box( "solModule", ModuleX/2.0, ModuleY/2.0, ModuleL/2.0 );
  //increase module box length so that it contains lightguide and PMT photocathode. We will need to change some of the positioning arguments of sub-volumes accordingly.
  //This is to prevent geometry overlaps!
  G4Box *solModule = new G4Box( "solModule", ModuleX/2.0, ModuleY/2.0, ModuleLtotal/2.0 );
  G4LogicalVolume *logModule = new G4LogicalVolume( solModule, GetMaterial("Special_Air"), "logModule" );

  G4Box *solIronPl = new G4Box( "solIronPl", PlateXHalf/2.0, PlateY/2.0, IronPlThick/2.0 );
  G4LogicalVolume *logIronPl = new G4LogicalVolume( solIronPl, GetMaterial("Iron"), "logIronPl" );

  // ****Scintillator**** 
  // is a Sensitive Detector of type CAL
  G4Box *solScinPl = new G4Box( "solScinPl" , PlateXHalf/2.0, PlateY/2.0, ScinPlThick/2.0 );
  G4LogicalVolume *logScinPl = new G4LogicalVolume( solScinPl, GetMaterial("EJ232"), "logScinPl" );
  
  G4SDManager *sdman = fDetCon->fSDman;

  G4String HCalScintSDname = "Harm/HCalScint";
  G4String HCalScintcollname = "HCalScintHitsCollection";
  G4SBSCalSD *HCalScintSD = NULL;

  if( !((G4SBSCalSD*) sdman->FindSensitiveDetector(HCalScintSDname)) ){
    G4cout << "Adding HCal Scintillator Sensitive Detector to SDman..." << G4endl;
    HCalScintSD = new G4SBSCalSD( HCalScintSDname, HCalScintcollname );
    sdman->AddNewDetector(HCalScintSD);
    (fDetCon->SDlist).insert(HCalScintSDname);
    fDetCon->SDtype[HCalScintSDname] = G4SBS::kCAL;
    //fDetCon->SDarm[HCalScintSDname] = kHarm;

    (HCalScintSD->detmap).depth = 1;

    G4int default_ntbins = 25;
    
    fDetCon->SetThresholdTimeWindowAndNTimeBins( HCalScintSDname, 0.0*MeV, 100.0*ns, 25 );
  }
  logScinPl->SetSensitiveDetector(HCalScintSD);

  
  
  if( (fDetCon->StepLimiterList).find( HCalScintSDname ) != (fDetCon->StepLimiterList).end() ){
    logScinPl->SetUserLimits(  new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
  }

  // Scintillator Wrap
  G4Box *sBox1 = new G4Box("sBox1", (PlateX + AlFoilThick)/2.0,
			   (PlateY + 0.5*AlFoilThick)/2.0,     
			   (ScinPlThick + AlFoilThick)/2.0 ); //Why only half the foil thickness in y? sbox1 dimensions are (7.31 cm, 7.255 cm, 0.51 cm)

  G4double FoilThickness = AlFoilThick;
  G4double DeltaX = PlateX/2.0 + FoilThickness;
  G4double DeltaY = PlateY/2.0 + FoilThickness/2.0; //Why only half thickness in y?
  G4double DeltaZ = ScinPlThick/2.0 + FoilThickness;

  G4Box* sBox2 = new G4Box("sBoxSc2", DeltaX, DeltaY, DeltaZ); //( 7.32 cm, 7.26 cm, 0.52 cm ); 

  //Total thickness:
  // sBox1 = (14.62 cm, 14.51 cm, 1.02 cm);
  // sBox2 = (14.64 cm, 14.52 cm, 1.04 cm);
  
  G4ThreeVector pos(0.0, FoilThickness/2.0+AlFoilThick/4.0, 0.0); //( 0, 0.015 cm, 0 ); 

  G4SubtractionSolid* solScinPlWrap = new G4SubtractionSolid("sScinPlWr", sBox2, sBox1, 0, pos);
  //The subtraction solid shifts the y position of sBox1 up by .015 cm, so that on the +Y side we
  //The thickness of the remaining solid is
  // (0.01 at +X edge, 0.01 at -X edge, 0.02 at -Y edge, 0 at +Y edge, and 0.01 at -Z edge and 0.01 at +Z edge): Is that what we want?
  // Seems incorrect, but probably doesn't matter much.

  //This never gets placed, so ignore geometry!
  G4LogicalVolume *logScinPlWrap = new G4LogicalVolume( solScinPlWrap, GetMaterial("Aluminum"), "lScinPlWr" );
    
  G4double Xpos, Ypos, Zpos;
  pos.set(0.0,0.0,0.0);

  //add offset to correct for new module thickness:
  //We want z0 to be z0 = -(LightGuideZ + 2*LightGuideX + 0.5 cm)/2;
  // instead everything here is -ModuleL/2
  // z0desired = -ModuleL/2 + zoffset --> zoffset = z0desired + ModuleL/2
  //G4double zoffset = ModuleL/2.0 - (TotalPlatesL + 2.0*LightGuideX + 0.5*cm)/2.0;
  
  for( int ii=0; ii<NumberOfLayers; ii++ ) {
    G4double iron_gap = 2.0 * ii* PlateGaps; 

    Xpos = ( ScinToLgGap +  LightGuideY/2.0 + PlateXHalf/2.0 ); // = 0.1 cm + 0.25 cm + 6.95/2 cm
    Zpos = (-ModuleLtotal/2 + ContainerThick + IronPlThick/2 + ii * (IronPlThick + ScinPlThick)) + iron_gap;
    pos.set(Xpos, 0.0, Zpos);
    new G4PVPlacement( 0, pos, logIronPl, "FePl", logModule, false, ii );

    pos.set(-Xpos, 0.0, Zpos); 
    new G4PVPlacement( 0, pos, logIronPl, "FePl", logModule, false, NumberOfLayers + ii );
      
    G4double scin_gap = ( 1 + 2 * (ii)) * PlateGaps;
    Zpos = (-ModuleLtotal/2 + ContainerThick + IronPlThick + ScinPlThick/2.0 
	    + ii * (IronPlThick + ScinPlThick)) + scin_gap;
    //pos.setZ(Zpos + zoffset);
    pos.set(Xpos,0.0,Zpos);      
    new G4PVPlacement( 0, pos, logScinPl, "ScPlL", logModule, false, ii );
     
    pos.set(-Xpos,0.0,Zpos); 
    new G4PVPlacement( 0, pos, logScinPl, "ScPlR", logModule, false, ii );
    //pos.setY(-AlFoilThick); //this command appears to have no effect!     
  }

  //Are these statements relevant?
  pos.setZ(0.0);
  Ypos = (PlateY/2 + ScinToLgGap + LightGuideY/2)*cm;
  Zpos = (-(ModuleL - TotalPlatesL)/2 + ContainerThick )*cm;

  // ****Lightguide****
  G4Box* sBox = new G4Box( "sBox" , LightGuideX/2.0, LightGuideY/2.0, LightGuideZ/2.0 );

  G4double pDx1   =  LightGuideX;
  G4double pDx2   =  LightGuideX;
  G4double pDy1   =  LightGuideY; //5 mm = full length at -dz/2
  G4double pDx3   =  2.7*cm;
  G4double pDx4   =  2.7*cm;
  G4double pDy2   =  2.7*cm;
  G4double pDz    =  2.0*LightGuideX; //28.4 cm = full length along z.
  G4double pTheta =  0*degree; 
  G4double pPhi   =  90*degree;
  G4double pAlp1  =  0*degree;
  G4double pAlp2  =  pAlp1;

  Zpos = LightGuideZ/2.0 + pDz/2.0; // = 46 cm + 14.2 cm = 60.2 cm
  pos.set(0.0, 0.0, Zpos);

  G4Trap *solTrap = new G4Trap( "sTrap1",
				pDz/2,   pTheta,
				pPhi,    pDy1/2,
				pDx1/2,  pDx2/2,
				pAlp1,   pDy2/2,
				pDx3/2,  pDx4/2,
				pAlp2);

  //*****TEST*****
  //G4LogicalVolume *test = new G4LogicalVolume(solTrap,GetMaterial("Air"),"test");
  //new G4PVPlacement(lightg,G4ThreeVector(3*m,3*m,3*m),test,"testphys",motherlog,false,0);

  G4UnionSolid *sol = new G4UnionSolid( "usol", sBox, solTrap, 0, pos );
  //G4LogicalVolume *test1 = new G4LogicalVolume(sol,GetMaterial("Air"),"test");
  //new G4PVPlacement(0,G4ThreeVector(3.2*m,3.2*m,3.2*m),test1,"test1phys",motherlog,false,0);
  //*****END*****

  Zpos = (-(ModuleLtotal - TotalPlatesL)/2 + ContainerThick );
  pos.set( 0.0, 0.0, Zpos );

  G4ThreeVector pos_lg = pos;
  
  //NOTE: "Ligd" OVERLAPS WITH MOTHERVOLUME
  G4LogicalVolume *logLightG = new G4LogicalVolume( sol, GetMaterial("BC484"), "lLiGd");  
  new G4PVPlacement(mRotateZ , pos , logLightG , "Ligd" , logModule , false , 0 , true );
  new G4LogicalSkinSurface( "Lightguide Skin", logLightG, GetOpticalSurface("osWLSToAir") );   
  //G4LogicalBorderSurface* WLSToAir = new G4LogicalBorderSurface("WLSToAir", phyLightG , phyModule , OpWLSToAir);
  
  G4Box *solWLSPaper = new G4Box( "sLiGd" , LightGuideX/2.0, 0.01 *cm, LightGuideZ/2.0);
  G4LogicalVolume *logWLSPaper = new G4LogicalVolume(solWLSPaper, GetMaterial("Paper"), "lPaper");
 
  Ypos = PlateY/2 + ScinToLgGap + LightGuideY + 2.0*0.01;
  pos.setY(Ypos);
  //  phyWLSPaper = new G4PVPlacement(0 , pos , logWLSPaper , "Paper" ,logModule  , false , 0);
  
  // ****PMT****
  double radiuscath   = 2.7*cm;
  pDz    = 2.0*LightGuideX;
  double LightGYpos = PlateY/2 + ScinToLgGap + LightGuideY/2;
  G4double PMT_L = 0.5*cm;

  G4Tubs *solCathod = new G4Tubs("scath", 0.0*cm, radiuscath/2.0 , PMT_L/2.0 , 0.0*deg , 360.0*deg);
  G4LogicalVolume *logCathod = new G4LogicalVolume(solCathod,GetMaterial("Glass_HC"), "lcath");

  // logCathod is the SD, assigned to ECalSD which detects optical photons
  G4String HCalSDname = "Harm/HCal";
  G4String HCalcollname = "HCalHitsCollection";
  G4SBSECalSD* HCalSD = NULL;

  if( !((G4SBSECalSD*) sdman->FindSensitiveDetector(HCalSDname)) ) {
    G4cout << "Adding HCal PMT Sensitive Detector to SDman..." << G4endl;
    HCalSD = new G4SBSECalSD( HCalSDname, HCalcollname );
    sdman->AddNewDetector(HCalSD);
    (fDetCon->SDlist).insert(HCalSDname);
    fDetCon->SDtype[HCalSDname] = G4SBS::kECAL;
    //fDetCon->SDarm[HCalSDname] = G4SBS::kHarm;
    (HCalSD->detmap).depth = 0;
  }
  logCathod->SetSensitiveDetector(HCalSD);

  pos.set( 0.0, 0.0, 0.0);
  G4Box* abox0 = new G4Box("abox0", ModuleX/2.0, ModuleY/2.0, ModuleL/2.0 );
  G4Box* abox1 = new G4Box("abox1", ModuleX/2.0 - ContainerThick, ModuleY/2.0 - ContainerThick, 
  			   ModuleL/2.0 - ContainerThick );

  G4SubtractionSolid* solContar = new G4SubtractionSolid("hollow-box", abox0 , abox1, 0, pos);

  //G4ThreeVector lg_container_offset = pos_lg - pos;
  G4ThreeVector pos_container( 0.0, 0.0, -ModuleLtotal/2.0 + ModuleL/2.0 );
  G4ThreeVector lg_container_offset = pos_lg - pos_container; 
  
  G4SubtractionSolid *container_minus_lightguide = new G4SubtractionSolid( "container-lightguide", solContar, sol, 0, lg_container_offset ); //prevent overlap between steel container and light-guide!
  
  G4LogicalVolume *logContar = new G4LogicalVolume( container_minus_lightguide, GetMaterial("Steel"), "lCont" );  
  new G4PVPlacement(mRotateZ, pos_container, logContar, "Cont", logModule, false, 0);

  //test to see if a module builds with PMT
  //new G4PVPlacement(0 , G4ThreeVector(3*m,3*m,3*m+LightGuideZ/2.0+(3*LightGuideX)/2.0-2*ContainerThick-radiuscath/2.0), logCathod , "physcathode" , motherlog , false , 0);
  //new G4PVPlacement(0,G4ThreeVector(3*m,3*m,3*m),logModule,"phys",motherlog,false,0);

  const G4int    NSupportsPerModule  = 8;
  const G4double SupportPlateDy      = 0.0*cm; 
  const G4double SupportPlateDyExtra = ( NColumns / NSupportsPerModule - 1.0 ) * SupportPlateDy * cm; //-- 3 extra plates

  //Make Mother Volume which will house all modules

  G4RotationMatrix *hcalrm = new G4RotationMatrix;
  hcalrm->rotateY(f48D48ang);
  //hcalrm->rotateZ(90*degree); //Do we really need this? Is this a needless complication?
  G4RotationMatrix *modrot = new G4RotationMatrix;
  modrot->rotateZ(-90*degree );

  //G4Box *solCalo = new G4Box( "sHCalo", CaloX/2.0, (CaloY+SupportPlateDyExtra)/2.0, LightGuideZ/2.0+(3*LightGuideX)/2.0+2*ContainerThick+radiuscath);
  //Since SupportPlateDyExtra is zero, don't needlessly confuse by including it!
  G4Box *solCalo = new G4Box( "sHCalo", CaloX/2.0, CaloY/2.0, CaloL/2.0 );
  G4LogicalVolume *logCalo = new G4LogicalVolume( solCalo, GetMaterial("Air"), "lHCalo" );
  new G4PVPlacement( hcalrm, G4ThreeVector( -hcalr*sin(f48D48ang), VerticalOffset, hcalr*cos(f48D48ang) ), logCalo, "HCal Mother", motherlog, false, 0, false);

  if( (fDetCon->StepLimiterList).find( "lHCalo" ) != (fDetCon->StepLimiterList).end() ){
    logCalo->SetUserLimits( new G4UserLimits( 0.0, 0.0, 0.0, DBL_MAX, DBL_MAX ) );
    G4String sdname = "Harm/HCAL_box";
    G4String collname = "HCAL_boxHitsCollection";
    G4SBSCalSD *HCALboxSD = NULL;
    if( !((G4SBSCalSD*) sdman->FindSensitiveDetector(sdname)) ){
      G4cout << "Adding HCALbox sensitive detector to SDman..." << G4endl;
      HCALboxSD = new G4SBSCalSD( sdname, collname );
      sdman->AddNewDetector( HCALboxSD );
      (fDetCon->SDlist).insert( sdname );
      fDetCon->SDtype[sdname] = G4SBS::kCAL;
      (HCALboxSD->detmap).depth = 0;
      logCalo->SetSensitiveDetector( HCALboxSD );
    }
  }

  fDetCon->InsertSDboundaryVolume( logCalo->GetName(), HCalScintSDname );
  fDetCon->InsertSDboundaryVolume( logCalo->GetName(), HCalSDname );
  
  G4int copyid = 0;
  for(int ii = 0; ii < NColumns; ii++) {
    for(int jj = 0; jj < NRows; jj++) {

      G4double xtemp = -CaloX/2.0/1.001 + (ModuleX/2.0 + ModuleX*ii);
      G4double ytemp = -CaloY/2.0/1.001 + (ModuleY/2.0 + ModuleY*jj);
      G4double zmodule = -CaloL/2.0/1.001 + ModuleLtotal/2.0;
      //pos.set( xtemp, ytemp, zmodule);
      G4double zcathode = zmodule + ModuleLtotal/2.0 + 0.25*cm;
      
      //G4ThreeVector pos_cathode( xtemp, ytemp, 
      
      //     new G4PVPlacement( 0, G4ThreeVector(xtemp, ytemp, LightGuideZ/2.0+(3*LightGuideX)/2.0-2*ContainerThick-radiuscath/2.0),
      //		 logCathod, "physcathode", logCalo, false, copyid );
      new G4PVPlacement( 0, G4ThreeVector(xtemp, ytemp, zcathode),
			 logCathod, "physcathode", logCalo, false, copyid ); 
      //new G4PVPlacement(modrot, pos, logModule, "module", logCalo, false, copyid);
      new G4PVPlacement(0, G4ThreeVector(xtemp,ytemp,zmodule), logModule, "module", logCalo, false, copyid);
      
      (HCalSD->detmap).Row[copyid] = jj;
      (HCalSD->detmap).Col[copyid] = ii;
      (HCalSD->detmap).LocalCoord[copyid] = G4ThreeVector(xtemp, ytemp, 0.0);

      (HCalScintSD->detmap).Row[copyid] = jj;
      (HCalScintSD->detmap).Col[copyid] = ii;
      (HCalScintSD->detmap).LocalCoord[copyid] = G4ThreeVector(xtemp,ytemp,0.0);

      copyid++;
    }
  }

  //--- Front plate
  G4double FrontPlatedZ = 2.0*2.54*cm;
  pos.set(0.0,0.0,-(CaloL+FrontPlatedZ*1.1)/2.0 );
  G4Box *solFrontPlate = new G4Box( "sFrontPlate" , CaloX/2.0, (CaloY+SupportPlateDyExtra)/2.0, FrontPlatedZ/2.0 );
  G4LogicalVolume *logFrontPlate = new G4LogicalVolume( solFrontPlate , GetMaterial("Iron") , "lFrontPlate" );
  //didn't place it yet


  // Visualization
  
  // Iron
  G4VisAttributes * IronPlVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  logIronPl->SetVisAttributes(IronPlVisAtt);
 
  // Scintillator
  G4VisAttributes * ScinPlVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,1.0));
  logScinPl->SetVisAttributes(ScinPlVisAtt);

  // LightGuide
  G4VisAttributes * LightGVisAtt = new G4VisAttributes(G4Colour(0.54, 0.53, 0.79));
  LightGVisAtt->SetForceSolid(true);
  logLightG->SetVisAttributes(LightGVisAtt);

  // PMT
  logCathod->SetVisAttributes(G4Colour::Blue());
  
  // Module container vis
  G4VisAttributes *logContarVis  = new G4VisAttributes(G4Colour::Blue());
  logContarVis->SetForceWireframe(true);
  logContarVis->SetColor(G4Colour::Grey());
  logContar->SetVisAttributes(logContarVis);  

  G4VisAttributes *logModuleVis  = new G4VisAttributes(G4Colour::Blue());
  logModuleVis ->SetForceWireframe(true);
  logModule->SetVisAttributes(logModuleVis);

  // Scint Plate
  logScinPlWrap->SetVisAttributes(G4VisAttributes::GetInvisible());

  // Mother Volume
  //G4VisAttributes * caloVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  //caloVisAtt->SetForceWireframe(true);
  logCalo->SetVisAttributes(G4VisAttributes::GetInvisible());
}

// void G4SBSHArmBuilder::MakeTracker( G4LogicalVolume *worldlog)
// {
//   // 2016/06/23 (jc2): Moved SIDIS tracker to its own function
//   //Let's make a simple tracker: 5 planes of GEMs, equally spaced in z, separation in z between planes of 10 cm. Then total length of tracker is ~50 cm + about 1.6 cm
//   //G4double SBStracker_dist = fRICHdist - 0.3*m;
//   G4double SBStracker_dist = fSBS_tracker_dist;
//   G4ThreeVector SBStracker_pos( -SBStracker_dist * sin( f48D48ang ), 0.0, SBStracker_dist * cos( f48D48ang ) );

//   G4RotationMatrix *SBStracker_rot_I = new G4RotationMatrix(G4RotationMatrix::IDENTITY);

//   //Just a test:
//   //SBStracker_rot_I->rotateY( 14.0*deg );

//   G4RotationMatrix *SBStracker_rot = new G4RotationMatrix;
//   SBStracker_rot->rotateY( f48D48ang );

//   G4Box *SBStracker_box = new G4Box("SBStracker_box", 32.0*cm, 102.0*cm, 22.0*cm );

//   G4LogicalVolume *SBStracker_log = new G4LogicalVolume( SBStracker_box, GetMaterial("Air"), "SBStracker_log" );

//   //For consistency with BigBite, place "Tracker" volume before GEMs are placed in it:
//   new G4PVPlacement( SBStracker_rot, SBStracker_pos, SBStracker_log, "SBStracker_phys", worldlog, false, 0 );

//   int ngems_SBStracker = 5;
//   vector<double> zplanes_SBStracker, wplanes_SBStracker, hplanes_SBStracker;

//   G4double zspacing_SBStracker = 10.0*cm;
//   G4double zoffset_SBStracker = -20.0*cm;

//   for(int i=0; i<ngems_SBStracker; i++ ){
//     zplanes_SBStracker.push_back( zoffset_SBStracker + i*zspacing_SBStracker );
//     wplanes_SBStracker.push_back( 60.0*cm );
//     hplanes_SBStracker.push_back( 200.0*cm );
//   }

//   G4SBSTrackerBuilder trackerbuilder(fDetCon);

//   //(fDetCon->TrackerArm)[fDetCon->TrackerIDnumber] = kHarm; //H arm is "1"

//   trackerbuilder.BuildComponent( SBStracker_log, SBStracker_rot_I, G4ThreeVector(0,0,0), 
//       ngems_SBStracker, zplanes_SBStracker, wplanes_SBStracker, hplanes_SBStracker, G4String("Harm/SBSGEM") );

//   SBStracker_log->SetVisAttributes(G4VisAttributes::GetInvisible());
// }

void G4SBSHArmBuilder::MakeTracker(G4LogicalVolume *motherlog, G4int nplanes){
  //      G4double SBStracker_dist = fRICHdist - 0.3*m; //distance to the front of the SBS tracker
  
  G4ThreeVector SBS_midplane_pos( -(f48D48dist + 0.5*f48D48depth)*sin(f48D48ang), 0.0, (f48D48dist+0.5*f48D48depth)*cos(f48D48ang) );

  G4RotationMatrix *SBStracker_rot_I = new G4RotationMatrix(G4RotationMatrix::IDENTITY);

  //Just a test:
  //SBStracker_rot_I->rotateY( 14.0*deg );

  G4RotationMatrix *SBStracker_rot = new G4RotationMatrix;
  SBStracker_rot->rotateY( f48D48ang );
  SBStracker_rot->rotateX( fSBS_tracker_pitch );

  int ngems_SBStracker = nplanes;
  
  G4double zspacing_SBStracker = 10.0*cm;
  G4double zoffset_SBStracker = -(ngems_SBStracker-1)*zspacing_SBStracker/2.0;
  G4double LtrackerBox = (ngems_SBStracker-1)*zspacing_SBStracker+4.0*cm;
  printf("%f %f\n", zoffset_SBStracker, LtrackerBox/2.0);
  
  G4Box *SBStracker_box = new G4Box("SBStracker_box", 32.0*cm, 102.0*cm, LtrackerBox/2.0 );

  G4LogicalVolume *SBStracker_log = new G4LogicalVolume( SBStracker_box, GetMaterial("Air"), "SBStracker_log" );

  //G4double RICH_yoffset = (fRICHdist - (f48D48dist + 0.5*f48D48depth) )*sin( fSBS_tracker_pitch );
  //We want the center of the first tracker plane to be directly above the point at a distance D along the SBS center line
  //from the target:
  // xmidplane + L*sin

  G4double Tracker_yoffset = (fSBS_tracker_dist - (f48D48dist + 0.5*f48D48depth) )*tan(fSBS_tracker_pitch);

  G4ThreeVector FirstPlane_pos( -fSBS_tracker_dist*sin(f48D48ang),
				Tracker_yoffset,
				fSBS_tracker_dist*cos(f48D48ang) );
  
  G4ThreeVector SBS_tracker_axis = (FirstPlane_pos - SBS_midplane_pos).unit();
			     
  
  // G4ThreeVector RICH_pos( -fRICHdist*sin(f48D48ang), RICH_yoffset, fRICHdist*cos(f48D48ang) );

  // G4ThreeVector SBS_tracker_axis = (RICH_pos - SBS_midplane_pos).unit();
  // G4ThreeVector SBS_tracker_pos = RICH_pos - 0.3*m*SBS_tracker_axis;

  G4ThreeVector SBS_tracker_pos = FirstPlane_pos - zoffset_SBStracker*SBS_tracker_axis;
  
  printf("sbs_tracker_pos: %f, %f, %f\n", SBS_tracker_pos.x(), SBS_tracker_pos.y(), SBS_tracker_pos.z() );
  
  new G4PVPlacement( SBStracker_rot, SBS_tracker_pos, SBStracker_log, "SBStracker_phys", motherlog, false, 0 );

  vector<double> zplanes_SBStracker, wplanes_SBStracker, hplanes_SBStracker;

  for(int i=0; i<ngems_SBStracker; i++ ){
    zplanes_SBStracker.push_back( zoffset_SBStracker + i*zspacing_SBStracker );
    wplanes_SBStracker.push_back( 60.0*cm );
    hplanes_SBStracker.push_back( 200.0*cm );
  }

  G4SBSTrackerBuilder trackerbuilder(fDetCon);

  //(fDetCon->TrackerArm)[fDetCon->TrackerIDnumber] = kHarm; //H arm is "1"

  trackerbuilder.BuildComponent( SBStracker_log, SBStracker_rot_I, G4ThreeVector(0,0,0), 
      ngems_SBStracker, zplanes_SBStracker, wplanes_SBStracker, hplanes_SBStracker, G4String("Harm/SBSGEM") );

  SBStracker_log->SetVisAttributes(G4VisAttributes::GetInvisible());
}

void G4SBSHArmBuilder::MakeElectronModeSBS(G4LogicalVolume *motherlog){
  MakeTracker( motherlog );
  MakeRICH_new( motherlog , true);

  G4SBS::Exp_t exptype = fDetCon->fExpType;
  if( exptype == G4SBS::kA1n ){ //HCAL AND LAC:
    MakeHCALV2( motherlog, fHCALvertical_offset );
    MakeLAC( motherlog );
  } else { //TDIS or NDVCS: LAC only:
    MakeLAC( motherlog );
  }
}

void G4SBSHArmBuilder::MakeRICH_new( G4LogicalVolume *motherlog, bool extended_snout){
  if(!extended_snout)fRichSnoutExtension = 0;

  G4RotationMatrix *rot_RICH = new G4RotationMatrix;
  rot_RICH->rotateY( f48D48ang );
  rot_RICH->rotateX( fSBS_tracker_pitch );
  rot_RICH->rotateZ( 180.0*deg );

  G4double RICH_yoffset = (fRICHdist - (f48D48dist + 0.5*f48D48depth) )*sin( fSBS_tracker_pitch );
  G4ThreeVector RICHcoord_global( -fRICHdist*sin( f48D48ang ), RICH_yoffset, fRICHdist*cos( f48D48ang ) );

  G4ThreeVector SBS_midplane_pos( -(f48D48dist + 0.5*f48D48depth)*sin(f48D48ang), 0.0, (f48D48dist+0.5*f48D48depth)*cos(f48D48ang) );
  
  G4ThreeVector RICH_zaxis = (RICHcoord_global - SBS_midplane_pos).unit();

  G4ThreeVector RICH_xaxis( cos(f48D48ang), 0.0, sin(f48D48ang) );
  G4ThreeVector RICH_yaxis = (RICH_zaxis.cross(RICH_xaxis)).unit();
  
  //  G4ThreeVector RICH_zaxis( RICHcoord_global.unit() );
  //G4ThreeVector RICH_yaxis( 0.0, 1.0, 0.0 );
  //G4ThreeVector RICH_xaxis( (RICH_yaxis.cross( RICH_zaxis ) ).unit() );

  
  
  G4double RICHbox_w = 164.0*cm, RICHbox_h = 284.0*cm, RICHbox_thick = 126.0*cm;
  if(extended_snout)RICHbox_thick+=fRichSnoutExtension;
  //RICHbox_h = 10.0*m;
  
  G4Box *RICHbox = new G4Box("RICHbox", RICHbox_w/2.0, RICHbox_h/2.0, RICHbox_thick/2.0  );

  G4String RadiatorGas_Name = fRICHgas;
  // if( electronmode ){
  //   printf("SBS in electron mode: using CO2 radiator for RICH\n");
  //   RadiatorGas_Name = "CO2";
  // }
  
  G4LogicalVolume *RICHbox_log = new G4LogicalVolume( RICHbox, GetMaterial(RadiatorGas_Name), "SBS_RICH_log" );

  G4cout << "Building RICHbox with gas " << RadiatorGas_Name << G4endl;
  //At the end, we will rotate it by 180 degrees about z:
  
  //Define the origin with the x and z coordinates at "bottom left" corner of the box (everything will be centered in y)
  G4ThreeVector origin( -RICHbox_w/2.0, 0.0, -RICHbox_thick/2.0 + 0.75*mm); //Add 0.75 mm to account for thickness of entry window!
  if(extended_snout)origin.setZ(origin.z()+fRichSnoutExtension);
  G4cout << "origin.z() (mm) " << origin.z()/mm << G4endl;
  
  G4double inch = 2.54*cm;
  
  //Mounting plate for PMT array:
  G4double MountPlate_width = 111.0*inch;
  G4double MountPlate_height = 42.28*inch;
  G4double MountPlate_thick = 1.50*inch;

  G4Box *MountPlate = new G4Box( "MountPlate", MountPlate_width/2.0, MountPlate_height/2.0, MountPlate_thick/2.0 );

  G4double MountPlate_opening_width = 58.17*inch;
  G4double MountPlate_opening_height = 25.25*inch;

  G4Box *MountPlate_window = new G4Box( "MountPlate_window", MountPlate_opening_width/2.0, MountPlate_opening_height/2.0, MountPlate_thick/2.0 + 1.0*cm );

  G4double MountPlate_window_center_x = 0.0;
  G4double MountPlate_window_center_y = 11.2699*inch + MountPlate_opening_height/2.0 - MountPlate_height/2.0;
 
  G4ThreeVector MountPlate_window_center( 0.0, MountPlate_window_center_y, 0.0 );

  G4SubtractionSolid *MountPlate_cut = new G4SubtractionSolid( "MountPlate_cut", MountPlate, MountPlate_window, 0, MountPlate_window_center );

  G4LogicalVolume *MountPlate_log = new G4LogicalVolume( MountPlate_cut, GetMaterial("Al"), "MountPlate_log" );

  
  
  G4RotationMatrix *MountPlate_rot = new G4RotationMatrix;
  MountPlate_rot->rotateZ( 90.0*deg );
  MountPlate_rot->rotateX( -50.0*deg );

  G4double mountplate_frontedge_z = 1.137*cm;
  G4double mountplate_frontedge_x = 80.361*cm;

  G4ThreeVector Det_zaxis( -sin( 50.0*deg ), 0.0, cos( 50.0*deg ) );
  G4ThreeVector Det_yaxis( 0, 1, 0 );
  G4ThreeVector Det_xaxis( (Det_yaxis.cross( Det_zaxis ) ).unit() );

  G4ThreeVector mountplate_frontedge_pos( mountplate_frontedge_x, 0.0, mountplate_frontedge_z );

  G4ThreeVector mountplate_center_pos = mountplate_frontedge_pos + origin + MountPlate_height/2.0 * Det_xaxis + MountPlate_thick/2.0 * Det_zaxis;

  new G4PVPlacement( MountPlate_rot, mountplate_center_pos, MountPlate_log, "SBSRICHMountPlate_pv", RICHbox_log, false, 0 );
  
  //PMT array plate #1:
  G4double ArrayPlate1_width = 66.0*inch;
  G4double ArrayPlate1_height = 34.625*inch;
  G4double ArrayPlate1_thick = 3.625*inch;
  G4Box *ArrayPlate1 = new G4Box("ArrayPlate1", ArrayPlate1_width/2.0, ArrayPlate1_height/2.0, ArrayPlate1_thick/2.0 );

  G4double ArrayPlate1_opening_width = 58.0*inch;
  G4double ArrayPlate1_opening_height = 24.687*inch;

  G4Box *ArrayPlate1_window = new G4Box("ArrayPlate1_window", ArrayPlate1_opening_width/2.0, ArrayPlate1_opening_height/2.0, ArrayPlate1_thick/2.0 + 1.0*cm );

  G4SubtractionSolid *ArrayPlate1_cut = new G4SubtractionSolid( "ArrayPlate1_cut", ArrayPlate1, ArrayPlate1_window, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *ArrayPlate1_log = new G4LogicalVolume( ArrayPlate1_cut, GetMaterial("Steel"), "ArrayPlate1_log" );
  
  //PMT array plate #2:
  G4double ArrayPlate2_width = 63.0*inch;
  G4double ArrayPlate2_height = 30.625*inch;
  G4double ArrayPlate2_thick = 3.625*inch;

  G4Box *ArrayPlate2 = new G4Box( "ArrayPlate2", ArrayPlate2_width/2.0, ArrayPlate2_height/2.0, ArrayPlate2_thick/2.0 );

  G4double ArrayPlate2_opening_width = 58.0*inch;
  G4double ArrayPlate2_opening_height = 24.687*inch;

  G4Box *ArrayPlate2_window = new G4Box("ArrayPlate2_window", ArrayPlate2_opening_width/2.0, ArrayPlate2_opening_height/2.0, ArrayPlate2_thick/2.0 + 1.0*cm );

  G4SubtractionSolid *ArrayPlate2_cut = new G4SubtractionSolid( "ArrayPlate2_cut", ArrayPlate2, ArrayPlate2_window, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *ArrayPlate2_log = new G4LogicalVolume( ArrayPlate2_cut, GetMaterial("Steel"), "ArrayPlate2_log" );

  G4ThreeVector arrayplate1_pos = mountplate_center_pos + Det_xaxis * MountPlate_window_center_y - Det_zaxis * 0.5 * ArrayPlate1_thick ;

  G4ThreeVector arrayplate2_pos = arrayplate1_pos - Det_zaxis * 0.5 * ( ArrayPlate1_thick + ArrayPlate2_thick );

  new G4PVPlacement( MountPlate_rot, arrayplate1_pos, ArrayPlate1_log, "SBSRICHarrayplate1_pv", RICHbox_log, false, 0 );
  new G4PVPlacement( MountPlate_rot, arrayplate2_pos, ArrayPlate2_log, "SBSRICHarrayplate2_pv", RICHbox_log, false, 0 );
  
  //O-ring plate:
  G4double Oring_plate_width = 111.02*inch;
  G4double Oring_plate_height = 43.5*inch;
  G4double Oring_plate_thick = 0.75*inch;
  
  G4Box *Oring_plate = new G4Box("Oring_plate", Oring_plate_width/2.0, Oring_plate_height/2.0, Oring_plate_thick/2.0 );
  
  G4double Oring_plate_opening_height = 35.12*inch;
  G4double Oring_plate_opening_width = 104.400*inch;

  G4Box *Oring_plate_opening = new G4Box("Oring_plate_opening", Oring_plate_opening_width/2.0, Oring_plate_opening_height/2.0, Oring_plate_thick/2.0+1.0*cm );

  G4double Oring_plate_opening_x0 = 3.31*inch + Oring_plate_opening_width/2.0;
  G4double Oring_plate_opening_y0 = 5.00*inch + Oring_plate_opening_height/2.0 - Oring_plate_height/2.0;

  G4ThreeVector Oring_plate_opening_center( 0.0, Oring_plate_opening_y0, 0.0 );

  G4SubtractionSolid *Oring_plate_cut = new G4SubtractionSolid( "Oring_plate_cut", Oring_plate, Oring_plate_opening, 0,
								Oring_plate_opening_center );

  G4LogicalVolume *Oring_plate_log = new G4LogicalVolume( Oring_plate_cut, GetMaterial("Al"), "Oring_plate_log" );

  G4double oring_frontedge_z = 12.60*mm;
  G4double oring_frontedge_x = 75.5*cm;

  G4ThreeVector oring_plate_frontedge_pos( oring_frontedge_x, 0.0, oring_frontedge_z );

  G4ThreeVector oring_plate_center_pos = oring_plate_frontedge_pos + origin + Oring_plate_height/2.0 * Det_xaxis + Oring_plate_thick/2.0 * Det_zaxis;

  new G4PVPlacement( MountPlate_rot, oring_plate_center_pos, Oring_plate_log, "SBSRICHOringPlate_pv", RICHbox_log, false, 0 );  
  
  //Next: Bottom plate:
  G4double bottomplate_width = 111.02*inch;
  G4double bottomplate_height = 46.21*inch;
  G4double bottomplate_thick = 1.0*inch;

  G4Box *BottomPlate = new G4Box( "BottomPlate", bottomplate_width/2.0, bottomplate_height/2.0, bottomplate_thick/2.0 );
  G4LogicalVolume *BottomPlate_log = new G4LogicalVolume( BottomPlate, GetMaterial("Al"), "BottomPlate_log" );

  G4double bottomplate_xcenter = 3.00*cm + bottomplate_thick/2.0;
  G4double bottomplate_zcenter = 1.50*inch + bottomplate_height/2.0;

  G4ThreeVector bottomplate_pos = origin + G4ThreeVector( bottomplate_xcenter, 0.0, bottomplate_zcenter );
  
  G4RotationMatrix *bottomplate_rot = new G4RotationMatrix;
  bottomplate_rot->rotateZ( 90.0*deg );
  bottomplate_rot->rotateX( 90.0*deg );

  new G4PVPlacement( bottomplate_rot, bottomplate_pos, BottomPlate_log, "SBSRICHbottomplate_pv", RICHbox_log, false, 0 );
  
  //Next: Front plate:
  G4double frontplate_width = 111.02*inch;
  G4double frontplate_height = 27.75*inch;
  G4double frontplate_thick = 1.00*inch;

  G4Box *FrontPlate = new G4Box("FrontPlate", frontplate_width/2.0, frontplate_height/2.0, frontplate_thick/2.0 );

  G4double frontplate_window_width = 187.7*cm;
  G4double frontplate_window_height = 46.4*cm;

  G4Box *FrontPlate_window = new G4Box("FrontPlate_window", frontplate_window_width/2.0, frontplate_window_height/2.0, frontplate_thick/2.0+1.0*cm );

  G4double frontplate_window_y0 = 6.60*cm + frontplate_window_height/2.0 - frontplate_height/2.0;

  G4SubtractionSolid *FrontPlate_cut = new G4SubtractionSolid( "FrontPlate_cut", FrontPlate, FrontPlate_window, 0, G4ThreeVector( 0, frontplate_window_y0, 0 ) );

  G4LogicalVolume *FrontPlate_log = new G4LogicalVolume( FrontPlate_cut, GetMaterial("Al"), "FrontPlate_log" );

  G4double frontplate_z0 = 0.50*inch + frontplate_thick/2.0;
  G4double frontplate_x0 = frontplate_height/2.0;
  G4ThreeVector frontplate_pos = origin + G4ThreeVector( frontplate_x0, 0.0, frontplate_z0 );

  G4RotationMatrix *frontplate_rot = new G4RotationMatrix;
  frontplate_rot->rotateZ( 90.0 * deg );

  new G4PVPlacement( frontplate_rot, frontplate_pos, FrontPlate_log, "SBSRICHFrontPlate_pv", RICHbox_log, false, 0 );
  
  //Next: front angle wedge:

  G4double FrontWedge_width = 4.0*inch;
  G4double FrontWedge_angle = 40.0*deg;
  G4double FrontWedge_height = 4.0*inch;

  G4double FrontWedge_h0 = 1.0*inch;
  G4double FrontWedge_h1 = FrontWedge_h0 + FrontWedge_width * tan( FrontWedge_angle );
  
  G4double FrontWedge_length = 111.02*inch;

  G4double Frontwedge_pgon_x[5] = { 0.0, 0.0, (FrontWedge_height - FrontWedge_h0)/tan( FrontWedge_angle ), FrontWedge_width, FrontWedge_width };
  G4double Frontwedge_pgon_y[5] = { 0.0, FrontWedge_h0, FrontWedge_height, FrontWedge_height, 0.0 };
  vector<G4TwoVector> frontwedge_pgon;
  for( G4int i=0; i<5; i++ ){
    frontwedge_pgon.push_back( G4TwoVector( Frontwedge_pgon_x[i], Frontwedge_pgon_y[i] ) );
  }
  
  //G4Trd *FrontWedge = new G4Trd( "FrontWedge", FrontWedge_length/2.0, FrontWedge_length/2.0, FrontWedge_h0/2.0, FrontWedge_h1/2.0, FrontWedge_width/2.0 );
  G4ExtrudedSolid *FrontWedge = new G4ExtrudedSolid( "FrontWedge", frontwedge_pgon, FrontWedge_length/2.0, G4TwoVector(0,0), 1.0, G4TwoVector(0,0), 1.0 );
  
  G4LogicalVolume *FrontWedge_log = new G4LogicalVolume( FrontWedge, GetMaterial("Al"), "FrontWedge_log" );

  G4double frontwedge_x = frontplate_height;
  G4double frontwedge_z = 0.50*inch;
  G4double frontwedge_y = 0.0;

  G4ThreeVector frontwedge_pos = origin + G4ThreeVector( frontwedge_x, frontwedge_y, frontwedge_z );

  G4RotationMatrix *frontwedge_rot = new G4RotationMatrix;
  //frontwedge_rot->rotateZ( 90.0*deg );
  frontwedge_rot->rotateY( 90.0*deg );
  frontwedge_rot->rotateX( 90.0*deg );
  
  new G4PVPlacement( frontwedge_rot, frontwedge_pos, FrontWedge_log, "SBSRICHfrontwedge", RICHbox_log, false, 0 );
  
  //Next: Back window plate:
  G4double BackPlate_width = 111.02*inch;
  G4double BackPlate_height = 34.95*inch;
  G4double BackPlate_thick = 1.00*inch;

  G4Box *BackPlate = new G4Box( "BackPlate", BackPlate_width/2.0, BackPlate_height/2.0, BackPlate_thick/2.0 );

  G4double BackPlate_window_height = 59.0*cm;
  G4double BackPlate_window_width = 257.0*cm;

  G4Box *BackPlate_Window = new G4Box( "BackPlate_Window", BackPlate_window_width/2.0, BackPlate_window_height/2.0, BackPlate_thick/2.0+1.0*cm );
  
  G4double backwindow_y0 = 11.6*cm + BackPlate_window_height/2.0 - BackPlate_height/2.0;

  G4SubtractionSolid *BackPlate_cut = new G4SubtractionSolid( "BackPlate_cut", BackPlate, BackPlate_Window, 0, G4ThreeVector( 0, backwindow_y0, 0 ) );

  G4LogicalVolume *BackPlate_log = new G4LogicalVolume( BackPlate_cut, GetMaterial("Al"), "BackPlate_log" );

  G4double BackPlate_z = 121.18*cm + BackPlate_thick/2.0;
  G4double BackPlate_y = 0.0;
  G4double BackPlate_x = 0.0 + BackPlate_height/2.0;

  G4ThreeVector backplate_pos = origin + G4ThreeVector( BackPlate_x, BackPlate_y, BackPlate_z );
  G4RotationMatrix *backplate_rot = new G4RotationMatrix;
  backplate_rot->rotateZ( 90.0*deg );

  new G4PVPlacement( backplate_rot, backplate_pos, BackPlate_log, "SBSRICH_BackPlate_pv", RICHbox_log, false, 0 );

  G4LogicalVolume* snout_log;
  if(extended_snout && fRichSnoutExtension>0){
    G4double snout_width = 198.8*cm;
    G4double snout_height = 55.5*cm;
    G4double snout_thick = inch/16.0;
    
    G4Box *snoutbox = new G4Box("snoutbox", snout_width/2.0, snout_height/2.0, fRichSnoutExtension/2.0 );
    G4Box *snouthollow = new G4Box("snouthollow", snout_width/2.0-snout_thick, snout_height/2.0-snout_thick, fRichSnoutExtension/2.0+snout_thick );
    
    G4SubtractionSolid* snout_sol = new G4SubtractionSolid( "snout_sol", snoutbox, snouthollow, 0, G4ThreeVector(0,0,0) );
    snout_log = new G4LogicalVolume(snout_sol, GetMaterial("Al"), "snout_log");
    
    G4double snout_x0 = 6.60*cm + snout_height/2.0-5.55*cm;
    G4double snout_z0 = -0.75*mm - fRichSnoutExtension/2.0;
    G4ThreeVector snout_pos = origin + G4ThreeVector( snout_x0, 0.0, snout_z0 );
    
    G4RotationMatrix *snout_rot = new G4RotationMatrix;

    snout_rot->rotateZ( 90.0*deg );

    new G4PVPlacement( snout_rot, snout_pos, snout_log, "SBSRICH_snout_extension_pv", RICHbox_log, false, 0 );
    
    //G4Box *snoutgas = new G4Box("snoutgas", snout_width/2.0-snout_thick, snout_height/2.0-snout_thick, fRichSnoutExtension/2.0);
    //G4LogicalVolume* snoutgas_log = new G4LogicalVolume(snout_sol, GetMaterial("Al"), "snout_log");
    
  }
  
  //Next: Front Window frame:
  G4double frontframe_width = 198.8*cm;
  G4double frontframe_height = 55.5*cm;
  G4double frontframe_thick = 0.5*inch;
  
  G4Box *FrontFrame = new G4Box("FrontFrame", frontframe_width/2.0, frontframe_height/2.0, frontframe_thick/2.0 );

  G4double frontframe_window_width = 187.7*cm;
  G4double frontframe_window_height = 46.4*cm;

  G4Box *FrontFrame_window = new G4Box( "FrontFrame_window", frontframe_window_width/2.0, frontframe_window_height/2.0, frontframe_thick/2.0+1.0*cm );
  
  G4SubtractionSolid *FrontFrame_cut = new G4SubtractionSolid( "FrontFrame_cut", FrontFrame, FrontFrame_window, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *FrontFrame_log = new G4LogicalVolume( FrontFrame_cut, GetMaterial("Al"), "FrontFrame_log" );

  G4double frontframe_x0 = 6.60*cm + frontframe_window_height/2.0;
  G4double frontframe_z0 = -0.75*mm + frontframe_thick/2.0;
  G4ThreeVector frontframe_pos = origin + G4ThreeVector( frontframe_x0, 0.0, frontframe_z0 );
  if(extended_snout)frontframe_pos.setZ(frontframe_pos.z()-fRichSnoutExtension);
  
  G4RotationMatrix *frontframe_rot = new G4RotationMatrix;

  frontframe_rot->rotateZ( 90.0*deg );

  new G4PVPlacement( frontframe_rot, frontframe_pos, FrontFrame_log, "SBSRICH_FrontFrame_pv", RICHbox_log, false, 0 );
  
  //Next: Front window itself:

  G4Box *FrontWindow = new G4Box( "FrontWindow", 198.8*cm/2.0, 55.5*cm/2.0, 0.75*mm/2.0 );
  G4LogicalVolume *FrontWindow_log = new G4LogicalVolume( FrontWindow, GetMaterial("Al"), "FrontWindow_log" );

  G4ThreeVector frontwindow_pos = origin + G4ThreeVector( frontframe_x0, 0.0, frontframe_z0 + (0.75*mm + frontframe_thick)/2.0 );
  if(extended_snout)frontwindow_pos.setZ(frontwindow_pos.z()-fRichSnoutExtension);
  
  new G4PVPlacement( frontframe_rot, frontwindow_pos, FrontWindow_log, "SBSRICH_FrontWindow_pv", RICHbox_log, false, 0 );
  
  //Back window frame:
  G4double backframe_width = 268.1*cm;
  G4double backframe_height = 68.1*cm;
  G4double backframe_thick = 0.5*inch;

  G4Box *BackFrame = new G4Box("BackFrame", backframe_width/2.0, backframe_height/2.0, backframe_thick/2.0 );

  G4double backwindow_width = 257.0*cm;
  G4double backwindow_height = 59.0*cm;

  G4Box *BackFrame_window = new G4Box("BackFrame_window", backwindow_width/2.0, backwindow_height/2.0, backframe_thick/2.0+mm );

  G4SubtractionSolid *BackFrame_cut = new G4SubtractionSolid( "BackFrame_cut", BackFrame, BackFrame_window, 0, G4ThreeVector(0,0,0) );

  G4LogicalVolume *BackFrame_log = new G4LogicalVolume( BackFrame_cut, GetMaterial("Al"), "BackFrame_log" );

  G4double backframe_z0 = 121.18*cm + BackPlate_thick + 0.75*mm + 0.5*backframe_thick;
  G4double backframe_x0 = 11.6*cm + BackPlate_window_height/2.0;

  G4ThreeVector backframe_pos = origin + G4ThreeVector( backframe_x0, 0.0, backframe_z0 );
  new G4PVPlacement( frontframe_rot, backframe_pos, BackFrame_log, "SBSRICH_backframe_pv", RICHbox_log, false, 0 );
  
  //Back window:
  G4Box *BackWindow = new G4Box( "BackWindow", 268.1*cm/2.0, 68.1*cm/2.0, 0.75*mm/2.0 );
  G4LogicalVolume *BackWindow_log = new G4LogicalVolume( BackWindow, GetMaterial("Al"), "BackWindow_log" );

  G4ThreeVector backwindow_pos = backframe_pos - G4ThreeVector( 0.0, 0.0, 0.5*backframe_thick + 0.75*mm/2.0 );

  new G4PVPlacement( frontframe_rot, backwindow_pos, BackWindow_log, "SBSRICH_backwindow_pv", RICHbox_log, false, 0 );
  
  //TOPWELDMENT: neglect channels, use average thickness of
  G4double topweld_average_thick = 0.351*inch;
  G4double topweld_width = 111.02*inch;
  G4double topweld_height = 51.0*cm;

  G4Box *topweld = new G4Box("topweld", topweld_width/2.0, topweld_height/2.0, topweld_average_thick/2.0 );
  G4LogicalVolume *topweld_log = new G4LogicalVolume( topweld, GetMaterial("Al"), "topweld_log" );

  G4RotationMatrix *topweld_rot = new G4RotationMatrix;
  topweld_rot->rotateZ( 90.0*deg );
  topweld_rot->rotateX( 40.0*deg );

  G4double topweld_edge_z0 = 81.795*cm;
  G4double topweld_edge_x0 = 140.581*cm;

  G4ThreeVector topweld_zaxis( sin( 40.0*deg ), 0.0, cos( 40.0*deg ) );
  G4ThreeVector topweld_yaxis( 0, 1, 0 );
  G4ThreeVector topweld_xaxis( (topweld_yaxis.cross( topweld_zaxis ) ).unit() );

  G4ThreeVector topweld_edge_pos = origin + G4ThreeVector( topweld_edge_x0, 0.0, topweld_edge_z0 );
  G4ThreeVector topweld_center_pos = topweld_edge_pos + 0.5*topweld_average_thick*topweld_zaxis - 0.5*topweld_height * topweld_xaxis;

  new G4PVPlacement( topweld_rot, topweld_center_pos, topweld_log, "SBSRICH_topweld_pv", RICHbox_log, false, 0 );
  
  //Square tube:
  G4double sqtube_width = 111.02*inch;
  G4double sqtube_height = 5.0*inch;
  G4double sqtube_thick = 5.0*inch;

  G4Box *SqTube = new G4Box("SqTube", sqtube_width/2.0, sqtube_height/2.0, sqtube_thick/2.0 );

  G4double sqtube_wallthick = 0.188*inch;
  
  G4Box *SqTube_hole = new G4Box( "SqTube_hole", sqtube_width/2.0+mm, sqtube_height/2.0-sqtube_wallthick, sqtube_thick/2.0-sqtube_wallthick );

  G4SubtractionSolid *SqTube_cut = new G4SubtractionSolid( "SqTube_cut", SqTube, SqTube_hole, 0, G4ThreeVector(0,0,0) );
  G4LogicalVolume *SqTube_log = new G4LogicalVolume( SqTube_cut, GetMaterial("Al"), "SqTube_log");

  G4double sqtube_z0 = 121.18*cm + 1.00*inch - sqtube_thick/2.0;
  G4double sqtube_x0 = 88.763*cm + sqtube_height/2.0;
  G4ThreeVector sqtube_pos = origin + G4ThreeVector( sqtube_x0, 0.0, sqtube_z0 );

  G4RotationMatrix *rot_sqtube = new G4RotationMatrix;
  rot_sqtube->rotateZ( 90.0*deg );

  new G4PVPlacement( rot_sqtube, sqtube_pos, SqTube_log, "SBSRICH_sqtube_pv", RICHbox_log, false, 0 );
  
  //Now side walls:
  vector<G4TwoVector> sidewall_pgon;

  G4double x_pgon[7] = {0.0, 0.0, 78.0*cm, 117.4*cm - 4.0*inch, 117.4*cm - 4.0*inch, 117.4*cm, 117.4*cm};
  G4double y_pgon[7] = {0.0, 72.15*cm, 137.6*cm, 102.725*cm, 85.765*cm, 85.765*cm, 0.0 };

  for(G4int i=0; i<7; i++){
    sidewall_pgon.push_back( G4TwoVector( x_pgon[i], y_pgon[i] ) );
  }

  G4double sidewall_thick = 0.25*inch;

  G4ExtrudedSolid *SideWall = new G4ExtrudedSolid( "SideWall", sidewall_pgon, sidewall_thick/2.0,
						   G4TwoVector( 0, 0 ), 1.0, G4TwoVector( 0, 0 ), 1.0 );
  G4LogicalVolume *SideWall_log = new G4LogicalVolume( SideWall, GetMaterial("Al"), "SideWall_log" );

  G4ThreeVector sidewall1_pos = origin + G4ThreeVector( 3.0*cm, (111.02*inch + sidewall_thick)/2.0, 3.8*cm );
  G4ThreeVector sidewall2_pos = origin + G4ThreeVector( 3.0*cm, -(111.02*inch + sidewall_thick)/2.0, 3.8*cm );

  G4RotationMatrix *sidewall1_rot = new G4RotationMatrix;
  sidewall1_rot->rotateY( 90.0*deg );
  sidewall1_rot->rotateX( 90.0*deg );

  // G4RotationMatrix *sidewall2_rot = new G4RotationMatrix;
  // sidewall2_rot->rotateY( -90.0*deg );
  
  new G4PVPlacement( sidewall1_rot, sidewall1_pos, SideWall_log, "SBSRICH_sidewall1", RICHbox_log, false, 0 );
  new G4PVPlacement( sidewall1_rot, sidewall2_pos, SideWall_log, "SBSRICH_sidewall2", RICHbox_log, false, 0 );
  //Now start building aerogel wall:

  G4int nx_aero=5, ny_aero=17, nz_aero=5;
  
  //define tile dimensions:
  G4double aero_tile_w = 11.4*cm, aero_tile_thick = 1.13*cm;
  
  G4double tile_gap = (57.912*cm - nx_aero * aero_tile_w)/G4double( nx_aero+1 ); //six spacers:

  G4Box *aerogel_tile = new G4Box( "Aerogel_tile", aero_tile_w/2.0, aero_tile_w/2.0, aero_tile_thick/2.0 );
  G4LogicalVolume *Aerogel_tile_log = new G4LogicalVolume( aerogel_tile, GetMaterial("Aerogel"), "Aerogel_tile_log" );

  G4double Width_aerogel_wall = nx_aero * (aero_tile_w + tile_gap) + tile_gap; 
  G4double Height_aerogel_wall = ny_aero * (aero_tile_w + tile_gap) + tile_gap;
  G4double Thick_aerogel_wall = nz_aero * aero_tile_thick;
  
  G4Box *Aerogel_wall_container = new G4Box("Aerogel_wall_container", Width_aerogel_wall/2.0, Height_aerogel_wall/2.0, Thick_aerogel_wall/2.0 );
  G4LogicalVolume *Aerogel_wall_container_log = new G4LogicalVolume( Aerogel_wall_container, GetMaterial("Air"), "Aerogel_wall_container_log" );

  G4Box *horizontal_spacer = new G4Box( "horizontal_spacer", aero_tile_w/2.0, tile_gap/2.0, Thick_aerogel_wall/2.0 );
  G4Box *vertical_spacer = new G4Box( "vertical_spacer", tile_gap/2.0, Height_aerogel_wall/2.0, Thick_aerogel_wall/2.0 );

  G4LogicalVolume *h_spacer_log = new G4LogicalVolume( horizontal_spacer, GetMaterial("Tedlar"), "h_spacer_log" );
  G4LogicalVolume *v_spacer_log = new G4LogicalVolume( vertical_spacer, GetMaterial("Tedlar"), "v_spacer_log" );

  TString pv_name;

  G4int tile_copy = 0;
  
  for( G4int col = 0; col <= nx_aero; col++ ){
    //position vertical spacer (one instance per column)
    G4double xspacer = -Width_aerogel_wall/2.0 + 0.5*tile_gap + col * ( aero_tile_w + tile_gap );
    G4double yspacer = 0.0;
    G4double zspacer = 0.0;

    new G4PVPlacement( 0, G4ThreeVector(xspacer,yspacer,zspacer), v_spacer_log, pv_name.Format( "v_spacer_pv_col%d", col ).Data(),
		       Aerogel_wall_container_log, false, col );
    
    for( G4int row = 0; row <= ny_aero; row++ ){
      //position horizontal spacer:
      xspacer = -Width_aerogel_wall/2.0 + tile_gap + 0.5*aero_tile_w + col * ( aero_tile_w + tile_gap );
      yspacer = -Height_aerogel_wall/2.0 + 0.5*tile_gap + row * (aero_tile_w + tile_gap );
      zspacer = 0.0;

      if( col < nx_aero ){
	new G4PVPlacement( 0, G4ThreeVector(xspacer,yspacer,zspacer), h_spacer_log, pv_name.Format( "h_spacer_pv_col%d_row%d", col, row ).Data(),
			   Aerogel_wall_container_log, false, row + col * ny_aero );
      }
      
      for( G4int iz = 0; iz < nz_aero; iz++ ){
	//position aerogel tile:
	G4double x_tile = -Width_aerogel_wall/2.0 + tile_gap + 0.5*aero_tile_w + col * ( aero_tile_w + tile_gap );
	G4double y_tile = -Height_aerogel_wall/2.0 + tile_gap + 0.5*aero_tile_w + row * ( aero_tile_w + tile_gap );
	G4double z_tile = -Thick_aerogel_wall/2.0 + (iz + 0.5)*aero_tile_thick;

	if( col < nx_aero && row < ny_aero ){
	  new G4PVPlacement( 0, G4ThreeVector(x_tile, y_tile, z_tile), Aerogel_tile_log, pv_name.Format( "Aerogel_tile_pv_%d", tile_copy ).Data(),
			     Aerogel_wall_container_log, false, tile_copy );
	  tile_copy++;
	}
      }
    }
  }

  G4double x0_aerogel_wall = 7.191*cm + Width_aerogel_wall/2.0;
  G4double y0_aerogel_wall = 0.0;
  G4double z0_aerogel_wall = 2.625*inch + Thick_aerogel_wall/2.0;

  G4ThreeVector pos_aerogel_wall( x0_aerogel_wall, y0_aerogel_wall, z0_aerogel_wall );
  pos_aerogel_wall += origin;

  
  

  //Also need 1 mm Al aerogel entry window and 3.2 mm UVT-lucite exit window:
  G4Box *aero_entry_window = new G4Box("aerogel_entry_window", Width_aerogel_wall/2.0, Height_aerogel_wall/2.0, 1.0*mm/2.0 );
  G4LogicalVolume *aero_entry_log = new G4LogicalVolume( aero_entry_window, GetMaterial("Al"), "aero_entry_log" );
  G4ThreeVector aero_entry_window_pos = pos_aerogel_wall + G4ThreeVector( 0.0, 0.0, -0.5*(Thick_aerogel_wall + 1.0*mm ) );
  
  
  
  //3.2 mm UVT-lucite aerogel exit window:
  G4Box *aero_exit_window = new G4Box( "aerogel_exit_window", Width_aerogel_wall/2.0, Height_aerogel_wall/2.0, 3.2*mm/2.0 );
  G4LogicalVolume *aero_exit_window_log = new G4LogicalVolume( aero_exit_window, GetMaterial("UVT_Lucite"), "Aero_exitwindow" );
  
  G4ThreeVector aero_exit_window_pos = pos_aerogel_wall + G4ThreeVector( 0.0, 0.0, 0.5*(Thick_aerogel_wall + 3.2*mm) );
  
  //For A1n ("Electron mode"), do not create/place aerogel wall components:

  if( fRICH_use_aerogel ){
    new G4PVPlacement( 0, pos_aerogel_wall, Aerogel_wall_container_log, "Aerogel_wall_container_pv", RICHbox_log, false, 0 );
    new G4PVPlacement( 0, aero_entry_window_pos, aero_entry_log, "SBSRICH_aero_entry_pv", RICHbox_log, false, 0 );
    new G4PVPlacement( 0, aero_exit_window_pos, aero_exit_window_log, "SBSRICH_aero_exit_pv", RICHbox_log, false, 0 );

  }
    
  G4double x0_mirror_center = 136.403*cm;
  G4double z0_mirror_center = -98.032*cm;
  
  G4double MirrorRadius = 220.0*cm;

  G4double Mirror_xmin = 5.542*cm;
  G4double Mirror_xmax = 85.045*cm;

  G4double Mirror_height = Mirror_xmax - Mirror_xmin; //along x axis
  G4double Mirror_width = 252.4*cm; //along y axis
  G4double Mirror_shellthick = 0.1932*cm; //~1% X0 of graphite:

  G4Sphere *RICH_mirror_shell = new G4Sphere( "RICH_mirror_shell", MirrorRadius, MirrorRadius + Mirror_shellthick, 0.0, twopi, 0.0, halfpi );
  G4Box *RICH_mirror_cutbox = new G4Box( "RICH_mirror_cutbox", Mirror_height/2.0, Mirror_width/2.0, RICHbox_thick/2.0 );

  G4ThreeVector box_centercoords( 0.5*(Mirror_xmin + Mirror_xmax) + origin.x(), 0.0, 0.0 );
  G4ThreeVector mirror_centercoords = origin + G4ThreeVector( x0_mirror_center, 0.0, z0_mirror_center );

  G4ThreeVector relative_coords = box_centercoords - mirror_centercoords;

  G4IntersectionSolid *Mirror_solid = new G4IntersectionSolid( "Mirror_solid", RICH_mirror_shell, RICH_mirror_cutbox, 0, relative_coords );
  G4LogicalVolume *Mirror_log = new G4LogicalVolume( Mirror_solid, GetMaterial("MirrorComposite"), "SBS_RICH_Mirror_log" );

  new G4PVPlacement( 0, mirror_centercoords, Mirror_log, "SBS_RICH_Mirror_pv", RICHbox_log, false, 0 );

  //Define the optical property (reflectivity) of the surface:
  new G4LogicalSkinSurface( "SBS_RICH_Mirrskin", Mirror_log, GetOpticalSurface("Mirrsurf") );

  //Now for the PMTs:
  ////////////////////////////////////////////////////////////////////////
  //                         !!!PMTS!!!                                 //
  ////////////////////////////////////////////////////////////////////////

  //cylinder to house PMTs:
  G4Tubs *PMTcylinder = new G4Tubs( "PMTcylinder", 0.0*cm, (1.86/2.0)*cm, 4.5*cm, 0.0, twopi ); //cylinder full of vacuum, mother volume for PMT
  G4LogicalVolume *PMTcylinder_log = new G4LogicalVolume( PMTcylinder, GetMaterial("BlandAir"), "PMTcylinder_log" );

  PMTcylinder_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  
  //Define the PMT windows as 1 mm-thick discs of "UVglass"; how thick are the windows really? 1 mm is a guess; let's go with 2 mm just to be conservative
  G4Tubs *PMTwindow = new G4Tubs( "PMTwindow", 0.0*cm, (1.66/2.0)*cm, 0.1*cm, 0.0, twopi ); 
  //Define the PMT photocathode as a thin disc of 10 micron thickness:
  G4Tubs *PMTcathode = new G4Tubs( "PMTcathode", 0.0*cm, (1.50/2.0)*cm, 0.005*mm, 0.0, twopi );
  //Define PMTtube as a stainless-steel tube that should butt up against collection cone to optically isolate PMTs from each other:
  G4Tubs *PMTtube    = new G4Tubs( "PMTtube", (1.66/2.0)*cm, (1.86/2.0)*cm, 4.5*cm, 0.0, twopi );
  G4Tubs *PMTendcap  = new G4Tubs( "PMTendcap", 0.0*cm, (1.66/2.0)*cm, 0.15*cm, 0.0, twopi ); //"end cap" for PMT

  //"Quartz window" is a different, sealed window that separates the PMT from the C4F10 environment. 2 mm thick
  G4Tubs *PMTQuartzWindow = new G4Tubs( "PMTQuartzWindow", 0.0*cm, (1.66/2.0)*cm, 0.1*cm, 0.0, twopi );
  G4Tubs *PMTWindowAirGap = new G4Tubs( "PMTWindowAirGap", 0.0*cm, (1.66/2.0)*cm, (0.5/2.0)*mm, 0.0, twopi ); //Air gap between quartz window and PMT entry window
  //CollectionCone is a light-collecting cone that increases the effective collection efficiency:
  G4Cons *CollectionCone = new G4Cons( "CollectionCone", 0.75*cm, (2.13/2.0)*cm, (2.13/2.0)*cm, (2.13/2.0)*cm, 0.75*cm, 0.0, twopi );

  //Total length = 9.0 cm (tube) + 1.5 cm (collection cone)
  G4double PMT_total_length = 10.5*cm;
  //    G4double PMT_max_radius = 1.065*cm;

  G4LogicalVolume *PMTwindow_log  = new G4LogicalVolume( PMTwindow, GetMaterial("UVglass"), "PMTwindow_log" );
  G4LogicalVolume *PMTcathode_log = new G4LogicalVolume( PMTcathode, GetMaterial("Photocathode_material_RICH"), "PMTcathode_log" );
  G4LogicalVolume *PMTWindowAirGap_log = new G4LogicalVolume( PMTWindowAirGap, GetMaterial("RICH_air"), "PMTWindowAirGap_log" ); //RICH_air is just air with a refractive index defined in the range of wavelengths of interest.
  
  //PMTcathode_log is the sensitive detector for the RICH:

  //  G4SDManager *fSDman = G4SDManager::GetSDMpointer();
  G4SDManager *sdman = fDetCon->fSDman;

  G4String RICHSDname = "Harm/RICH";
  G4String RICHcollname = "RICHHitsCollection";
  G4SBSRICHSD *RICHSD = NULL;

  if( !( RICHSD = (G4SBSRICHSD*) sdman->FindSensitiveDetector(RICHSDname) ) ){
    G4cout << "Adding RICH sensitive detector to SDman..." << G4endl;
    RICHSD = new G4SBSRICHSD( RICHSDname, RICHcollname );
    sdman->AddNewDetector( RICHSD );
    (fDetCon->SDlist).insert(RICHSDname);
    fDetCon->SDtype[RICHSDname] = G4SBS::kRICH;
    //fDetCon->SDarm[RICHSDname] = G4SBS::kHarm;

    PMTcathode_log->SetSensitiveDetector( RICHSD ); //This assigns the sensitive detector type "RICHSD" to the logical volume PMTcathode!
    (RICHSD->detmap).depth = 1;
  }

  fDetCon->InsertSDboundaryVolume( RICHbox_log->GetName(), RICHSDname );
  //We make this a hollow cylinder with length and radius approximately equal to that of the PMT housing, made of steel 
  //to approximate the material shielding the PMT.
  G4LogicalVolume *PMTtube_log    = new G4LogicalVolume( PMTtube, GetMaterial("Steel"), "PMTtube_log" ); 
  G4LogicalVolume *PMTendcap_log  = new G4LogicalVolume( PMTendcap, GetMaterial("Steel"), "PMTendcap_log" );
  G4LogicalVolume *PMTquartzwindow_log = new G4LogicalVolume( PMTQuartzWindow, GetMaterial("QuartzWindow"), "PMTQuartzWindow_log" );

  //Now we position PMT components inside PMT cylinder:
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), PMTtube_log, "PMTtube_pv", PMTcylinder_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, (-4.5+0.15)*cm ), PMTendcap_log, "PMTendcap_pv", PMTcylinder_log, false, 0 );
  //PMT photocathode is located at +4.5 cm - 2 mm - 0.5 mm - 2 mm - 0.005 mm 
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, (4.5-0.45-5e-4)*cm ), PMTcathode_log, "PMTcathode_pv", PMTcylinder_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, (4.5-0.35)*cm ), PMTwindow_log, "PMTwindow_pv", PMTcylinder_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, (4.5-0.2-0.025)*cm ), PMTWindowAirGap_log, "PMTWindowAirGap_pv", PMTcylinder_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, (4.5-0.1)*cm ), PMTquartzwindow_log, "PMTquartzwindow_pv", PMTcylinder_log, false, 0 );
  
  G4LogicalVolume *CollectionCone_log = new G4LogicalVolume( CollectionCone, GetMaterial("Steel"), "CollectionCone_log" );
  //Define a logical skin surface for the collection cone and assign it the same reflectivity as the mirror:
  new G4LogicalSkinSurface( "Coneskin", CollectionCone_log, GetOpticalSurface("Mirrsurf") );

  //Within the RICHbox, each PMT assembly unit is rotated so that its z-axis makes an angle of 50 degrees with respect to the 
  //local z axis of the RICHbox. Therefore, we rotate by an angle of 
  G4double PMT_rotation_angle = 50.0*degree;
  G4RotationMatrix *rot_PMT = new G4RotationMatrix;
  rot_PMT->rotateY( PMT_rotation_angle );

  G4int icopy_PMT_assembly = 0;

  // G4double xfp = 119.350*cm - RICHbox_dx/2.0;
  // G4double yfp = 0.0;
  // G4double zfp = 42.521*cm - RICHbox_dz/2.0;

  G4double xfp = 119.350*cm;
  G4double yfp = 0.0;
  G4double zfp = 47.601*cm;

  G4ThreeVector focalpoint_position = origin + G4ThreeVector( xfp, yfp, zfp );

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
      new G4PVPlacement( rot_PMT, Pos_temp, PMTcylinder_log, "SBS_RICH_PMT_assembly", RICHbox_log, false, icopy_PMT_assembly );
      G4double zcone = ztube + 5.25*cm;
      Pos_temp = PMT_position + zcone * PMT_zaxis;
      new G4PVPlacement( rot_PMT, Pos_temp, CollectionCone_log, "CollectionCone_pv", RICHbox_log, false, icopy_PMT_assembly );
      
      // new G4PVPlacement( rot_PMT, Pos_temp, PMTtube_log, "PMTtube_pv", RICHbox_log, false, icopy_PMT_assembly );
      // //Endcap of steel tube (keep optical photons originating from behind PMTs from hitting the cathode):
      // G4double zendcap = -PMT_total_length/2.0 + 0.15*cm;
      // Pos_temp = PMT_position + zendcap * PMT_zaxis;
      // new G4PVPlacement( rot_PMT, Pos_temp, PMTendcap_log, "PMTendcap_pv", RICHbox_log, false, icopy_PMT_assembly );
      // //Photocathode (this is the sensitive part!!):
      // G4double zcathode = ztube + 4.475*cm;
      // Pos_temp = PMT_position + zcathode * PMT_zaxis;
      // new G4PVPlacement( rot_PMT, Pos_temp, PMTcathode_log, "PMTcathode_pv", RICHbox_log, false, icopy_PMT_assembly );
      // //UV-glass PMT window:
      // G4double zwindow = zcathode + 0.075*cm;
      // Pos_temp = PMT_position + zwindow * PMT_zaxis;
      // new G4PVPlacement( rot_PMT, Pos_temp, PMTwindow_log, "PMTwindow_pv", RICHbox_log, false, icopy_PMT_assembly );
      // //Quartz window between PMT and gas:
      // G4double zquartz = zwindow + 0.2*cm;
      // Pos_temp = PMT_position + zquartz * PMT_zaxis;
      // new G4PVPlacement( rot_PMT, Pos_temp, PMTquartzwindow_log, "PMTquartzwindow_pv", RICHbox_log, false, icopy_PMT_assembly );
      // //Light collection cone:
      // G4double zcone = zquartz + 0.9*cm;
      // Pos_temp = PMT_position + zcone * PMT_zaxis;
      // new G4PVPlacement( rot_PMT, Pos_temp, CollectionCone_log, "CollectionCone_pv", RICHbox_log, false, icopy_PMT_assembly );

      (RICHSD->detmap).depth = 1;
      (RICHSD->detmap).Row[icopy_PMT_assembly] = irow;
      (RICHSD->detmap).Col[icopy_PMT_assembly] = icol;
      (RICHSD->detmap).LocalCoord[icopy_PMT_assembly] = G4ThreeVector(xtemp,ytemp,0.0);
      //G4ThreeVector pos_cathode_local = PMT_position + zcathode * PMT_zaxis;
      // (RICHSD->detmap).GlobalCoord[icopy_PMT_assembly] = RICH_centercoord_global + 
      // 	pos_cathode_local.X() * RICH_xaxis + 
      // 	pos_cathode_local.Y() * RICH_yaxis + 
      // 	pos_cathode_local.Z() * RICH_zaxis;
	      
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
  //                         !!!END OF PMTS!!!                          //
  ////////////////////////////////////////////////////////////////////////
  

  
  G4double x0_RICH = 6.6*cm + frontframe_window_height/2.0;
  G4double y0_RICH = 0.0;
  G4double z0_RICH = frontframe_thick;

  if(extended_snout)z0_RICH = z0_RICH-fRichSnoutExtension; 
    
  G4ThreeVector RICH_offset( x0_RICH, y0_RICH, z0_RICH ); //coordinates of center of entry window relative to origin.
  
  RICH_offset += origin;
  

  G4ThreeVector RICH_centercoord_global = RICHcoord_global - (-RICH_offset.x() * RICH_xaxis + RICH_offset.y() * RICH_yaxis + RICH_offset.z() * RICH_zaxis) + fRICHhorizontal_offset * RICH_xaxis + fRICHvertical_offset * RICH_yaxis;
  
  //We want to position the RICH box so that the center of the entry window is aligned with the SBS axis:
  new G4PVPlacement( rot_RICH, RICH_centercoord_global, RICHbox_log, "SBS_RICH_pv", motherlog, false, 0 );

  //Visualization attributes:
  G4VisAttributes *aero_tile_visatt = new G4VisAttributes( G4Colour( 0.0, 0.8, 0.8 ) );
  Aerogel_tile_log->SetVisAttributes( aero_tile_visatt );

  G4VisAttributes *tedlar_vis = new G4VisAttributes( G4Colour(0.3,0.3,0.3) );
  tedlar_vis->SetForceWireframe(true);
  v_spacer_log->SetVisAttributes( tedlar_vis );
  h_spacer_log->SetVisAttributes( tedlar_vis );

  Aerogel_wall_container_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  RICHbox_log->SetVisAttributes( G4VisAttributes::GetInvisible() );

  G4VisAttributes *RICHbox_visatt = new G4VisAttributes( G4Colour( 0.5, 0.5, 0.5 ) );
  RICHbox_visatt->SetForceWireframe(true);
  MountPlate_log->SetVisAttributes( RICHbox_visatt );
  ArrayPlate1_log->SetVisAttributes( RICHbox_visatt );
  ArrayPlate2_log->SetVisAttributes( RICHbox_visatt );
  Oring_plate_log->SetVisAttributes( RICHbox_visatt );
  BottomPlate_log->SetVisAttributes( RICHbox_visatt );
  FrontPlate_log->SetVisAttributes( RICHbox_visatt );
  FrontWedge_log->SetVisAttributes( RICHbox_visatt );
  BackPlate_log->SetVisAttributes( RICHbox_visatt );
  FrontFrame_log->SetVisAttributes( RICHbox_visatt );
  BackFrame_log->SetVisAttributes( RICHbox_visatt );
  topweld_log->SetVisAttributes( RICHbox_visatt );
  SqTube_log->SetVisAttributes( RICHbox_visatt );
  SideWall_log->SetVisAttributes( RICHbox_visatt );
  aero_entry_log->SetVisAttributes( RICHbox_visatt );
  G4VisAttributes *Window_visatt = new G4VisAttributes( G4Colour( 0.1, 0.8, 0.2 ) );
  Window_visatt->SetForceWireframe(true);
  FrontWindow_log->SetVisAttributes( Window_visatt );
  BackWindow_log->SetVisAttributes( Window_visatt );

  G4VisAttributes *lucite_visatt = new G4VisAttributes( G4Colour( 0.5, 0.05, 0.5 ) );
  lucite_visatt->SetForceWireframe(true);
  aero_exit_window_log->SetVisAttributes( lucite_visatt );

  G4VisAttributes *mirror_visatt = new G4VisAttributes( G4Colour( 0.0, 0.2, 0.8 ) );
  Mirror_log->SetVisAttributes( mirror_visatt );

  G4VisAttributes *PMT_window_visatt = new G4VisAttributes( G4Colour( 0.0, 1.0, 1.0 ) );
  PMT_window_visatt->SetForceWireframe(true);

  G4VisAttributes *PMT_collectionCone_visatt = new G4VisAttributes( G4Colour( 0.6, 0.3, 0.6 ) );
  G4VisAttributes *PMT_cathode_visatt = new G4VisAttributes( G4Colour( 0.5, 0.0, 1.0 ) );
  G4VisAttributes *PMTtube_visatt = new G4VisAttributes( G4Colour( 0.7, 0.7, 0.7 ) );

  PMTwindow_log->SetVisAttributes( PMT_window_visatt );
  PMTcathode_log->SetVisAttributes( PMT_cathode_visatt );
  PMTtube_log->SetVisAttributes( PMTtube_visatt );
  PMTendcap_log->SetVisAttributes( PMTtube_visatt );
  PMTquartzwindow_log->SetVisAttributes( PMT_window_visatt );
  
  if(extended_snout && fRichSnoutExtension>0){
    snout_log->SetVisAttributes( RICHbox_visatt );
  }
  
}


//////////// NOTE: THIS VERSION OF CDET CONSTRUCTION IS ACTUALLY NO LONGER USED IN ORDER TO AVOID CODE DUPLICATION/MAINTENANCE ISSUES /////////////////
void G4SBSHArmBuilder::MakeCDET( G4LogicalVolume *mother, G4double z0, G4double PlanesHOffset ){
  //z0 is the z position of the start of CDET relative to the HCal surface
  
  //R0 is the nominal distance from target to the start of CDET
  G4double R0 = fHCALdist + z0;
    
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

  G4String sdname = "Harm/CDET";
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

  sdname = "Harm/CDET_Scint";
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

    fDetCon->SetThresholdTimeWindowAndNTimeBins( sdname, 0.0*MeV, 50.0*ns, 25 );
  }
  
  //Now we need to define the coordinates of the "modules":
  //horizontal position within mother:
  G4double x0_modules[3] = { 0.5*(-70.986+31.014)*cm,
			     0.5*(-63.493+38.507)*cm,
			     0.5*(-56.0+46.0)*cm };
			     // -0.5*(-63.493+38.507)*cm,
			     // -0.5*(-70.986+31.014)*cm };
  //G4double x0_modules[3] = { 0.0, 0.0, 0.0 };
  
  //Number of rows per module:
  //G4int Nrow_module[3] = { 98, 98, 98, 98, 98, 98 };
  G4int Nrow_total = 588;

  //G4double y0[2][Nrow_total], pitch_angle[2][Nrow_total];

  //Nominal distance to CDET used to define the projective geometry:
  //G4double R0_CDET = 405.0*cm;
  //Nominal distance to planes:
  G4double R0_planes[2] = { R0 + Lz_scint/2.0 + 1.0*cm,
			    R0 + 3.0*Lz_scint/2.0 + 20.0*cm }; //allow for some small (1 cm) gaps between CH2 and start of 1st plane and between 1st and second planes...

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
	G4ThreeVector pos_strip( x0_modules[imod] + PlanesHOffset/2.0*pow(-1, plane) + ( col - 0.5 )*(Lx_scint+mylar_thick), R0_planes[plane] * tan( alpha ), z0 + R0_planes[plane] - R0 );

	G4RotationMatrix *rot_strip = new G4RotationMatrix;
	rot_strip->rotateY( col*pi );
	rot_strip->rotateX( alpha*pow(-1,col) );
	//rot_strip->rotateX( alpha );
	char physname[255];
	sprintf( physname, "Scint_module_phys_plane%d_row%d_col%d", plane+1, row+295, col+1 ); //In this construction, row varies from 295-588 for the top half
	G4String pvolname = physname;
	new G4PVPlacement( rot_strip, pos_strip, Scint_module, pvolname, mother, false, istrip );

	G4ThreeVector pos_pmt( x0_modules[imod] + PlanesHOffset/2.0*pow(-1, plane) + pow(-1,col+1)*(Lx_scint+mylar_thick+0.1*cm), R0_planes[plane] * tan( alpha ), z0 + R0_planes[plane] - R0 );
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

  Scint_module->SetVisAttributes( G4VisAttributes::GetInvisible() );
  
  G4VisAttributes *scintstrip_visatt = new G4VisAttributes( G4Colour( 0.8, 0, 0.8 ) );
  ScintStripLog->SetVisAttributes( scintstrip_visatt );

  G4VisAttributes *WLSfiber_visatt = new G4VisAttributes( G4Colour( 0, 0.8, 0.8 ) );
  WLSFiberLog->SetVisAttributes( WLSfiber_visatt );

  G4VisAttributes *scintwrap_visatt = new G4VisAttributes( G4Colour( 0.8, 0.8, 0 ) );
  ScintWrapLog->SetVisAttributes( scintwrap_visatt );
  
}


/////////////////////// NOTE: THIS METHOD IS ACTUALLY NO LONGER USED //////////////////////////////
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

  //////////////////////////////// BEGIN definition of RICH global coordinate transformations //////////////////////////////
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
  //////////////////////////////// END definition of RICH global coordinate transformations /////////////////////////////////

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

  //cylinder to house PMTs:
  G4Tubs *PMTcylinder = new G4Tubs( "PMTcylinder", 0.0*cm, (1.86/2.0)*cm, 4.5*cm, 0.0, twopi );
  G4LogicalVolume *PMTcylinder_log = new G4LogicalVolume( PMTcylinder, GetMaterial("BlandAir"), "PMTcylinder_log" );
  
  //Define the PMT windows as 1 mm-thick discs of "UVglass":
  G4Tubs *PMTwindow = new G4Tubs( "PMTwindow", 0.0*cm, (1.66/2.0)*cm, 0.05*cm, 0.0, twopi ); 
  //Define the PMT photocathode as a thin disc of 0.5 mm-thickness
  G4Tubs *PMTcathode = new G4Tubs( "PMTcathode", 0.0*cm, (1.50/2.0)*cm, 0.025*cm, 0.0, twopi );
  //Define PMTtube as a stainless-steel tube that should butt up against collection cone to optically isolate PMTs from each other:
  G4Tubs *PMTtube    = new G4Tubs( "PMTtube", (1.66/2.0)*cm, (1.86/2.0)*cm, 4.5*cm, 0.0, twopi );
  G4Tubs *PMTendcap  = new G4Tubs( "PMTendcap", 0.0*cm, (1.66/2.0)*cm, 0.15*cm, 0.0, twopi ); //end cap for PMT

  //"Quartz window" is a different, sealed window that separates the PMT from the C4F10 environment.
  G4Tubs *PMTQuartzWindow = new G4Tubs( "PMTQuartzWindow", 0.0*cm, (1.66/2.0)*cm, 0.15*cm, 0.0, twopi );
  //CollectionCone is a light-collecting cone that increases the effective collection efficiency:
  G4Cons *CollectionCone = new G4Cons( "CollectionCone", 0.75*cm, (2.13/2.0)*cm, (2.13/2.0)*cm, (2.13/2.0)*cm, 0.75*cm, 0.0, twopi );

  //Total length = 9.0 cm (tube) + 1.5 cm (collection cone)
  G4double PMT_total_length = 10.5*cm;
  //    G4double PMT_max_radius = 1.065*cm;

  G4LogicalVolume *PMTwindow_log  = new G4LogicalVolume( PMTwindow, GetMaterial("UVglass"), "PMTwindow_log" );
  G4LogicalVolume *PMTcathode_log = new G4LogicalVolume( PMTcathode, GetMaterial("Photocathode_material_RICH"), "PMTcathode_log" );

  //PMTcathode_log is the sensitive detector for the RICH:

  //  G4SDManager *fSDman = G4SDManager::GetSDMpointer();
  G4SDManager *sdman = fDetCon->fSDman;

  G4String RICHSDname = "Harm/RICH";
  G4String RICHcollname = "RICHHitsCollection";
  G4SBSRICHSD *RICHSD = NULL;

  if( !( RICHSD = (G4SBSRICHSD*) sdman->FindSensitiveDetector(RICHSDname) ) ){
    G4cout << "Adding RICH sensitive detector to SDman..." << G4endl;
    RICHSD = new G4SBSRICHSD( RICHSDname, RICHcollname );
    sdman->AddNewDetector( RICHSD );
    (fDetCon->SDlist).insert(RICHSDname);
    fDetCon->SDtype[RICHSDname] = G4SBS::kRICH;
    //fDetCon->SDarm[RICHSDname] = G4SBS::kHarm;

    PMTcathode_log->SetSensitiveDetector( RICHSD ); //This assigns the sensitive detector type "RICHSD" to the logical volume PMTcathode!
    (RICHSD->detmap).depth = 1;
  }
  //We make this a hollow cylinder with length and radius approximately equal to that of the PMT housing, made of steel 
  //to approximate the material shielding the PMT.
  G4LogicalVolume *PMTtube_log    = new G4LogicalVolume( PMTtube, GetMaterial("Steel"), "PMTtube_log" ); 
  G4LogicalVolume *PMTendcap_log  = new G4LogicalVolume( PMTendcap, GetMaterial("Steel"), "PMTendcap_log" );
  G4LogicalVolume *PMTquartzwindow_log = new G4LogicalVolume( PMTQuartzWindow, GetMaterial("QuartzWindow"), "PMTQuartzWindow_log" );

  //Now we position PMT components inside PMT cylinder:
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, 0 ), PMTtube_log, "PMTtube_pv", PMTcylinder_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, (-4.5+0.15)*cm ), PMTendcap_log, "PMTendcap_pv", PMTcylinder_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, (4.5-0.4-0.025)*cm ), PMTcathode_log, "PMTcathode_pv", PMTcylinder_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, (4.5-0.3-0.05)*cm ), PMTwindow_log, "PMTwindow_pv", PMTcylinder_log, false, 0 );
  new G4PVPlacement( 0, G4ThreeVector( 0, 0, (4.5-0.15)*cm ), PMTquartzwindow_log, "PMTquartzwindow_pv", PMTcylinder_log, false, 0 );
  
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
      new G4PVPlacement( rot_PMT, Pos_temp, PMTcylinder_log, "SBS_RICH_PMT_assembly", RICHbox_log, false, icopy_PMT_assembly );
      G4double zcone = ztube + 5.25*cm;
      Pos_temp = PMT_position + zcone * PMT_zaxis;
      new G4PVPlacement( rot_PMT, Pos_temp, CollectionCone_log, "CollectionCone_pv", RICHbox_log, false, icopy_PMT_assembly );
      
      // new G4PVPlacement( rot_PMT, Pos_temp, PMTtube_log, "PMTtube_pv", RICHbox_log, false, icopy_PMT_assembly );
      // //Endcap of steel tube (keep optical photons originating from behind PMTs from hitting the cathode):
      // G4double zendcap = -PMT_total_length/2.0 + 0.15*cm;
      // Pos_temp = PMT_position + zendcap * PMT_zaxis;
      // new G4PVPlacement( rot_PMT, Pos_temp, PMTendcap_log, "PMTendcap_pv", RICHbox_log, false, icopy_PMT_assembly );
      // //Photocathode (this is the sensitive part!!):
      // G4double zcathode = ztube + 4.475*cm;
      // Pos_temp = PMT_position + zcathode * PMT_zaxis;
      // new G4PVPlacement( rot_PMT, Pos_temp, PMTcathode_log, "PMTcathode_pv", RICHbox_log, false, icopy_PMT_assembly );
      // //UV-glass PMT window:
      // G4double zwindow = zcathode + 0.075*cm;
      // Pos_temp = PMT_position + zwindow * PMT_zaxis;
      // new G4PVPlacement( rot_PMT, Pos_temp, PMTwindow_log, "PMTwindow_pv", RICHbox_log, false, icopy_PMT_assembly );
      // //Quartz window between PMT and gas:
      // G4double zquartz = zwindow + 0.2*cm;
      // Pos_temp = PMT_position + zquartz * PMT_zaxis;
      // new G4PVPlacement( rot_PMT, Pos_temp, PMTquartzwindow_log, "PMTquartzwindow_pv", RICHbox_log, false, icopy_PMT_assembly );
      // //Light collection cone:
      // G4double zcone = zquartz + 0.9*cm;
      // Pos_temp = PMT_position + zcone * PMT_zaxis;
      // new G4PVPlacement( rot_PMT, Pos_temp, CollectionCone_log, "CollectionCone_pv", RICHbox_log, false, icopy_PMT_assembly );

      (RICHSD->detmap).depth = 1;
      (RICHSD->detmap).Row[icopy_PMT_assembly] = irow;
      (RICHSD->detmap).Col[icopy_PMT_assembly] = icol;
      (RICHSD->detmap).LocalCoord[icopy_PMT_assembly] = G4ThreeVector(xtemp,ytemp,0.0);
      //G4ThreeVector pos_cathode_local = PMT_position + zcathode * PMT_zaxis;
      // (RICHSD->detmap).GlobalCoord[icopy_PMT_assembly] = RICH_centercoord_global + 
      // 	pos_cathode_local.X() * RICH_xaxis + 
      // 	pos_cathode_local.Y() * RICH_yaxis + 
      // 	pos_cathode_local.Z() * RICH_zaxis;
	      
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
  //                         !!!END OF PMTS!!!                          //
  ////////////////////////////////////////////////////////////////////////

    
  //Place completed RICH geometry in the "world" volume:
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

  PMTcylinder_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  
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

  G4SBS::Exp_t exptype = fDetCon->fExpType;

  if( exptype == G4SBS::kGEp_BB ) {
    fGEPFPPoption = 4;
  }
  
  //Make analyzers first:
  // double anaheight = fGEP_CH2height;
  // double anawidth  = fGEP_CH2width;
  // //double anawidth  = 44.0*2.54*cm;
  
  // double anadepth  = 22.0*2.54*cm;
  // double ana1depth = fCH2thickFPP[0];
  // double ana2depth = fCH2thickFPP[1];

  //Don't assume these are the same:

  int nana = 2;
  const int ntrackermax = 3;
  int ntracker = 3;
  //int ngem[ntracker] = {6,5,5};
  int ngem[ntrackermax] = {6,5,5}; //defaults
  double GEM_z_spacing[ntrackermax];
  double trkr_zpos[ntrackermax];
  double trkr_yoff[ntrackermax] = {0.0, 0.0, 0.0};

  
  
  vector<G4String> SDnames;
  vector<G4String> AnalyzerMaterials;

  double zavg,zspace;
  
  switch( fGEPFPPoption ){ //1 = 8-plane front tracker of 6x(40x150 cm^2)+2x(60x200 cm^2) and 8-plane back tracker of 6x(60x200 cm^2)
  case 1:
  default:
    nana = 1;
    ntracker = 2;
    ngem[0] = 8;
    ngem[1] = 8;
    GEM_z_spacing[0] = 13.2*cm;
    GEM_z_spacing[1] = 11.3*cm;
    trkr_zpos[0] = 0.0;
    //trkr_zpos[1] = trkr_zpos[0] + ngem[0]*GEM_z_spacing[0] + fCH2thickFPP[0] + GEM_z_spacing[1];
    trkr_zpos[1] = 2.074*m;
    //trkr_yoff[0] = 
      //fGEP_CH2zpos[0] = trkr_zpos[0] + ngem[0]*GEM_z_spacing[0]; //upstream edge of CH2
    fGEP_CH2zpos[0] = 1.518*m - 0.5*fCH2thickFPP[0];
    SDnames.push_back( "Harm/FT" );
    SDnames.push_back( "Harm/FPP1" );
    fGEP_CH2width[0] = 28.88*2.54*cm; //about 73 cm
    fGEP_CH2height[0] = 84.0*2.54*cm; //about 213.36 cm
    fGEP_CH2yoff[0] = 25.78*cm;
    AnalyzerMaterials.push_back( G4String("CH2") );
    trkr_yoff[0] = 16.5*cm;
    trkr_yoff[1] = 34.1*cm;
    break;
  case 3: //6-plane FT + 22" CH2 + 5-plane FPP1 + let's use 3.5" steel + 5-plane FPP2:
    nana = 2;
    ntracker = 3;
    GEM_z_spacing[0] = 9.0*cm;
    GEM_z_spacing[1] = 10.0*cm;
    GEM_z_spacing[2] = 10.0*cm;

    fGEP_CH2zpos[0] = 58.53*cm;
    fGEP_CH2zpos[1] = 170.3*cm;
    fGEP_CH2width[0] = 28.88*2.54*cm; //about 73 cm
    fGEP_CH2height[0] = 84.0*2.54*cm; //about 213.36 cm
    fGEP_CH2width[1] = 28.88*2.54*cm; //about 73 cm
    fGEP_CH2height[1] = 84.0*2.54*cm; //about 213.36 cm

    //fCH2thickFPP[1] = 3.5*2.54*cm; //GEN-RP steel analyzer thickness: let's not actually
    //hard-code this8.
    
    SDnames.push_back("Harm/FT");
    SDnames.push_back("Harm/FPP1");
    SDnames.push_back("Harm/FPP2");

    zavg = 0.5*(fGEP_CH2zpos[0]+fCH2thickFPP[0] + fGEP_CH2zpos[1]); // mid-point between two analyzers 
    zspace = fGEP_CH2zpos[1] - (fGEP_CH2zpos[0]+fCH2thickFPP[0]); //spacing between first and second analyzers:

    trkr_zpos[0] = 0.0*cm;
    trkr_zpos[1] = zavg - 2.0*GEM_z_spacing[1];
    trkr_zpos[2] = fGEP_CH2zpos[1] + fCH2thickFPP[1] + GEM_z_spacing[2];

    AnalyzerMaterials.push_back( "CH2" );
    AnalyzerMaterials.push_back( "Steel" );
    
    break;
  case 4:
    //4 = GEP BigBite option: since 4 U/V layers and 1 X/Y layer will be in BigBite, what is left over will be 9 X/Y layers 60x200. Let's try the following layout:
    // 9 X/Y layers 50x120 in the front, plus 6 X/Y 60x150 in the back.
    nana = 1;
    ntracker = 2;

    ngem[0] = 6;
    ngem[1] = 6;
    GEM_z_spacing[0] = 9.0*cm;
    GEM_z_spacing[1] = 10.0*cm;

    trkr_zpos[0] = 0.0;
    trkr_zpos[1] = trkr_zpos[0] + ngem[0]*GEM_z_spacing[0] + fCH2thickFPP[0] + GEM_z_spacing[1];
    fGEP_CH2zpos[0] = trkr_zpos[0] + ngem[0]*GEM_z_spacing[0]; //upstream edge of CH2;
    SDnames.push_back( "Harm/FT" );
    SDnames.push_back( "Harm/FPP1" );
    fGEP_CH2width[0] = 28.88*2.54*cm; //about 73 cm
    fGEP_CH2height[0] = 84.0*2.54*cm; //about 213.36 cm

    AnalyzerMaterials.push_back( G4String("CH2") );
    
    break;
  case 2:
    //default: //2 = original layout: 6-plane FT of (40x150) cm^2 plus FPP1 and FPP2 trackers:
    nana = 2;
    ntracker = 3;
    GEM_z_spacing[0] = 9.0*cm;
    GEM_z_spacing[1] = 10.0*cm;
    GEM_z_spacing[2] = 10.0*cm;
    fGEP_CH2zpos[0] = 58.53*cm;
    fGEP_CH2zpos[1] = 170.3*cm;
    fGEP_CH2width[0] = 28.88*2.54*cm; //about 73 cm
    fGEP_CH2height[0] = 84.0*2.54*cm; //about 213.36 cm
    fGEP_CH2width[1] = 28.88*2.54*cm; //about 73 cm
    fGEP_CH2height[1] = 84.0*2.54*cm; //about 213.36 cm
    
    SDnames.push_back("Harm/FT");
    SDnames.push_back("Harm/FPP1");
    SDnames.push_back("Harm/FPP2");
    
    zavg = 0.5*(fGEP_CH2zpos[0]+fCH2thickFPP[0] + fGEP_CH2zpos[1]); // mid-point between two analyzers 
    zspace = fGEP_CH2zpos[1] - (fGEP_CH2zpos[0]+fCH2thickFPP[0]); //spacing between first and second analyzers:
    trkr_zpos[0] = 0.0*cm;
    trkr_zpos[1] = zavg - 2.0*GEM_z_spacing[1];
    trkr_zpos[2] = fGEP_CH2zpos[1] + fCH2thickFPP[1] + GEM_z_spacing[2];

    AnalyzerMaterials.push_back( "CH2" );
    AnalyzerMaterials.push_back( "CH2" );
    
    break;
  }

  G4VisAttributes *CH2anavisatt = new G4VisAttributes( G4Colour(0.0, 0.0, 1.0) );
  CH2anavisatt->SetForceWireframe(true);
  
  char ananame[50];
  char analogname[50];
  
  for( G4int ana=0; ana<nana; ana++ ){
    sprintf( ananame, "ana%dbox", ana );

    G4String anaboxname(ananame);
    G4String analogname = ananame;
    G4String anaphysname = ananame;

    analogname.append("_log");
    anaphysname.append("_phys");
    
    G4Box *anabox_temp = new G4Box( anaboxname, fGEP_CH2width[ana]/2.0, fGEP_CH2height[ana]/2.0, fCH2thickFPP[ana]/2.0 );
    G4LogicalVolume *analog_temp = new G4LogicalVolume( anabox_temp, GetMaterial(AnalyzerMaterials[ana]), analogname );
    fDetCon->InsertAnalyzerVolume( analog_temp->GetName() );

    analog_temp->SetVisAttributes( CH2anavisatt );
    
    G4ThreeVector anapos_temp = pos + G4ThreeVector( 0.0, fGEP_CH2yoff[ana], fGEP_CH2zpos[ana] + fCH2thickFPP[ana]/2.0 );
    new G4PVPlacement( 0, anapos_temp, analog_temp, anaphysname, Mother, false, 0, false );
  }

  // double zavg = 0.5*(170.3*cm + 58.53*cm+anadepth); //midpoint between first and second analyzers
  // double zspace = 170.3*cm - (58.53*cm+anadepth); //available space between first and second analyzers
  
  int i, j;
  //int ngem = 0;
  
  vector<double> gemz, gemw, gemh;
  
  // int ntracker = 3; //FT, FPP1, FPP2
  // int ngem[3] = {6,5,5};

  //we want equal air gaps between gem planes in FPP1. one of the GEMs will be placed at the midpoint.
  // ngap = ngem + 1 
  
  //double GEM_z_spacing[3] = {9.0*cm, zspace/double(ngem[1]+1), zspace/double(ngem[1]+1) };
  
  // vector<G4String> SDnames; 
  // SDnames.push_back("Harm/FT");
  // SDnames.push_back("Harm/FPP1");
  // SDnames.push_back("Harm/FPP2");

  //int TrackerID = 2;
  G4SBSTrackerBuilder trackerbuilder(fDetCon);

  for( i = 0; i<ntracker; i++){
    gemz.resize( ngem[i] );
    gemw.resize( ngem[i] );
    gemh.resize( ngem[i] );
    for( j = 0; j < ngem[i]; j++ ){
      gemz[j] = trkr_zpos[i] + ((double) j)*GEM_z_spacing[i];

      if( fGEPFPPoption != 4 ){
	if( i == 0 && j < 6 ){ //FT 40x150:
	  //	gemz[j] = ((double) j)*GEM_z_spacing[i];
	  gemw[j] = 40.0*cm;
	  gemh[j] = 150.0*cm;
	} else { //200x60
	  //gemz[j] = ((double) j)*10.0*cm + 1.2*m;
	  //gemz[j] = zavg + double(j-2)*GEM_z_spacing[i];
	  //	  gemz[i] = pairspac*((i-6)/2) + (i%2)*gemdsep + 1.2*m;
	  gemw[j] = 60.0*cm;
	  gemh[j] = 200.0*cm;
	}
      } else {
	if( i == 0 && j<2 ){
	  gemw[j] = 40.0*cm;
	  gemh[j] = 150.0*cm;
	} else {
	  gemw[j] = 60.0*cm;
	  gemh[j] = 200.0*cm;
	}
      }
    }
    //(fDetCon->TrackerArm)[fDetCon->TrackerIDnumber] = kHarm; //1 is H arm.
    G4bool ispolarimeter = false;
    if( i > 0 ) ispolarimeter = true;
    
    trackerbuilder.BuildComponent( Mother, rot, pos, ngem[i], gemz, gemw, gemh, SDnames[i], ispolarimeter, trkr_yoff[i] );
  }

  //CH2 analyzers:
  // ana1log->SetVisAttributes( CH2anavisatt );
  // ana2log->SetVisAttributes( CH2anavisatt );
  
}

//Sieve slit
void G4SBSHArmBuilder::MakeSBSSieveSlit(G4LogicalVolume *motherlog)
{
  printf("Building SBS sieve slit...\n");
  
  ////SIEVE PLATE OPTION - SSEEDS. LAST UPDATE - 9.5.20

  ////Dims per C. Soova JLab Drawing "HALL A - TEMPLATE - A00000 MAGNET - 48D48 DIPOLE SBS SIEVE PLATE ASSY SBS SIEVE PLATE"
  ////Strategy is to make a union solid of a thick and a thin plate to make base plate without holes. With this base plate, punch holes by creating a series of 77 subtraction solids where each has the previous subtraction solid as it's input. 
  ////Untested code complete 9.1.20 - waiting to test with visualization
  ////Determined placement by centering on bigbite arm (30.0 deg). Needs confirmation. 9.4.20
  ////9.5.20 - Debugging - y-direction holes appear to be rotated with incorrect reversal of rotation. Attempts to correct with sign flip over array of sub_yAng does not appear to change the output in visualization.
  ////9.5.20 - Debugging complete. Reversed order of x angles as fed by row in subtraction solids and rotated plate to -30.0 deg

  //Plate
  //Plate dims - Will union thick and thin plates to create final solid
  G4double inch = 2.54*cm;
  G4double ThinPlate_z = 1.0*inch;  //per print 0.995 - 0.985 inch 
  G4double ThinPlate_y = 30.0*inch; 
  G4double ThinPlate_x = 10.375*inch; 
  G4double ThickPlate_z = 2.0*inch; 
  G4double ThickPlate_y = 29.0*inch; 
  G4double ThickPlate_x = 9.375*inch;   
  
  //Plate positioning and orientation dims
  //G4double offset_z = 50.825*inch; //distance to center of thick plate face from center of target: 51.825", thickness of plate: 2.0"; offset_z = 51.825 - 2.0/2
  //G4double PlateAngDisp_theta = 30.0*deg;
  // G4double PlateAngDisp_theta = -f48D48ang;
  //G4double PlateAngDisp_phi = 180.0*deg;
  // G4RotationMatrix *SievePlateRot = new G4RotationMatrix; //rotates plate across each of its own axes
  // SievePlateRot->rotateZ(0*deg);
  // SievePlateRot->rotateX(0*deg);
  // SievePlateRot->rotateY(f48D48ang);

  // //Hole dims - x by y: 7 holes by 11 holes. Center hole: row 6, column 4. 
  // G4double holeSpace_y = 2.716*inch;
  // G4double holeSpace_x = 1.176*inch;
  //G4double offset_z = r48d48-30*cm-f48D48depth/2;
  G4double offset_z = f48D48dist - 1.6*m + 51.825*inch; //distance to BACK of plate.
  //G4double Ang_y = 3.0*deg; //Displacement angle in y direction - 3.0 per print
  //G4double Ang_x = 1.3*deg; //Displacement angle in x direction - 1.3 per print
  // G4double Ang_y = atan(holeSpace_y/offset_z);
  // G4double Ang_x = atan(holeSpace_x/offset_z);

  G4double HoleCenter_r = 0.125*inch;
  G4double Hole_r = 0.25*inch;
  G4double Hole_z = 3.0*inch; //Large enough to leave no solid volume at extreme displacements from center
  
  //Plate solids
  G4Box *ThinPlate = new G4Box("ThinPlate", ThinPlate_x/2.0, ThinPlate_y/2.0, ThinPlate_z/2.0);  
  G4Box *ThickPlate = new G4Box("ThickPlate", ThickPlate_x/2.0, ThickPlate_y/2.0, ThickPlate_z/2.0);

  //Union of plates
  G4VSolid *FullPlate = new G4UnionSolid("FullPlate", ThinPlate, ThickPlate);

  //Hole solids
  G4Tubs *Hole = new G4Tubs("Hole", 0, Hole_r, Hole_z, 0.0*deg, 360.0*deg);
  G4Tubs *HoleCenter = new G4Tubs("HoleCenter", 0, HoleCenter_r, Hole_z, 0.0*deg, 360.0*deg);

  //angular spacing of holes at the nominal distance of 51.825 inches
  G4double angspace_y = 3.0*deg;
  G4double angspace_x = 1.3*deg; 

  //G4SubtractionSolid *NextCut;
  G4SubtractionSolid *SievePlateCut;
  
  G4bool first = true;

  //Cut holes in the plate:

  cout << "cutting holes" << endl;
  
  
  for(G4int iy=-5; iy<=5; iy++ ){
    for(G4int ix=-3; ix<=3; ix++ ){
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
      
      G4double platecenter_z = 51.825*inch - ThickPlate_z/2.0;
      
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
	  NextCut = new G4SubtractionSolid( "NextCut", FullPlate, Hole, holerot, holecenterpos );
	} else { //it should be impossible to end up here:
	  NextCut = new G4SubtractionSolid( "NextCut", FullPlate, HoleCenter, holerot, holecenterpos );
	}
      } else {
	if( !(ix==0&&iy==0) ){
	  NextCut = new G4SubtractionSolid( "NextCut", SievePlateCut, Hole, holerot, holecenterpos );
	} else { //it should be impossible to end up here:
	  NextCut = new G4SubtractionSolid( "NextCut", SievePlateCut, HoleCenter, holerot, holecenterpos );
	}
      }

      SievePlateCut = NextCut;
      SievePlateCut->SetName("SBSSievePlateCut");
    }
  }

  G4cout << "finished holes..." << endl;
  
  G4LogicalVolume *SBSSievePlate_log = new G4LogicalVolume(SievePlateCut, GetMaterial("Lead"), "SBS_SievePlate_log");

  G4double SBSsieve_dist = offset_z - ThickPlate_z/2.0;
  //Placement:
  G4ThreeVector sievepos = SBSsieve_dist * G4ThreeVector( -sin(f48D48ang), 0, cos(f48D48ang) );
  G4RotationMatrix *sieverot = new G4RotationMatrix;
  sieverot->rotateY(f48D48ang);

  new G4PVPlacement(sieverot, sievepos, SBSSievePlate_log, "SBSSievePlate_phys", motherlog, false, 0);
  
  G4cout << "Sieve plate finished" << G4endl;
  
}
 
void G4SBSHArmBuilder::MakeLAC( G4LogicalVolume *motherlog ){

  //Define geometry parameters needed to build LAC:
  G4double W_LAC = 241.0*cm; //outer width of LAC module:
  G4double H_LAC = 446.0*cm; //outer height of LAC module
  //G4double D_LAC = 56.1*cm; //depth of LAC module
  G4double Winner_LAC = 217.0*cm; //inner width of LAC module
  G4double Hinner_LAC = 400.0*cm; //outer width of LAC module
  G4double PbStrip_thick = 0.20*cm; //thickness of lead strips
  G4double ScintStrip_thick = 1.5*cm; //thickness of NE110A strips

  G4int Nlayers_total = 33;
  G4int Nlayers_inner = 17;
  G4int Nlayers_outer = 16;

  G4int NstripsX = 24;
  G4int NstripsY = 40;
  
  G4double D_LAC = 55.9*cm; 
  
  G4Box *LAC_mother_box = new G4Box( "LAC_mother_box", 0.5*W_LAC+0.1*mm, 0.5*H_LAC+0.1*mm, 0.5*D_LAC+0.1*mm );

  G4LogicalVolume *log_LAC = new G4LogicalVolume( LAC_mother_box, GetMaterial("Air"), "log_LAC" );

  //G4VisAttributes *LAC_visatt = new G4VisAttributes( G4VisAttributes::GetInvisible() );

  log_LAC->SetVisAttributes( G4VisAttributes::GetInvisible() );
  
  G4double dLACbox = fLACdist + 0.5*D_LAC;

  G4ThreeVector LAC_pos( -dLACbox*sin(f48D48ang), fLACvertical_offset, dLACbox*cos(f48D48ang ) );

  G4ThreeVector LAC_zaxis( -sin(f48D48ang), 0.0, cos(f48D48ang) );
  G4ThreeVector LAC_yaxis(0,1,0);
  G4ThreeVector LAC_xaxis = LAC_yaxis.cross(LAC_zaxis).unit();

  LAC_pos = LAC_pos + fLACvertical_offset * LAC_yaxis + fLAChorizontal_offset * LAC_xaxis; 
  
  G4RotationMatrix *rot_LAC = new G4RotationMatrix;
  rot_LAC->rotateY(f48D48ang);

  new G4PVPlacement( rot_LAC, LAC_pos, log_LAC, "phys_LAC", motherlog, false, 0 );

  
  
  //Define sensitivity:
  // grab SD manager:
  G4SDManager *sdman = fDetCon->fSDman;
  
  G4String LACScintSDname  = "Harm/LACScint";
  G4String LACScintCollName = "LACScintHitsCollection";
  G4SBSCalSD *LACScintSD = NULL;

  if( !( (G4SBSCalSD*) sdman->FindSensitiveDetector(LACScintSDname)) ){
    G4cout << "Adding LAC Scint. Sensitive Detector to SDman..." << G4endl;
    LACScintSD = new G4SBSCalSD( LACScintSDname, LACScintCollName );
    sdman->AddNewDetector( LACScintSD );
    (fDetCon->SDlist).insert(LACScintSDname);
    fDetCon->SDtype[LACScintSDname] = G4SBS::kCAL;
    (LACScintSD->detmap).depth = 0;

    fDetCon->SetThresholdTimeWindowAndNTimeBins( LACScintSDname, 0.0*MeV, 100.0*ns, 25 );
  }
  fDetCon->InsertSDboundaryVolume( log_LAC->GetName(), LACScintSDname );
  
  //Now start populating the layers. In the absence of a better "guess", I will assume that there is one more layer's worth of "short" strips than "long" strips:

  TString ScintStrip_boxname,ScintStrip_logname,ScintStrip_physname;
  TString LeadSheet_boxname, LeadSheet_logname, LeadSheet_physname;

  G4int icopy_LACscint = 0; //Global copy number of scintillator strips
  G4int istack         = 0; //What we are after here is the sum of all energy depositions in individual "stacks" (each stack is coupled to a PMT at both ends)
  
  G4VisAttributes *LACscint_visatt = new G4VisAttributes( G4Colour(0.05, 0.9, 0.7) );
  G4VisAttributes *PbSheet_visatt = new G4VisAttributes( G4Colour( 0.3, 0.3, 0.3 ) );
  PbSheet_visatt->SetForceWireframe(true);

  ofstream mapfile("database/LAC_map.txt");

  TString currentline;
  currentline.Form( "# %10s, %10s, %10s, %10s, %18s, %18s, %18s",
		    "Stack", "Row", "Column", "Plane", "Xcenter (cm)", "Ycenter (cm)", "Zcenter (cm)" );

  mapfile << currentline << endl;
  
  for( G4int ilayer=0; ilayer<Nlayers_total; ilayer++ ){

    G4double zscint_layer = -D_LAC/2.0 + (ilayer+0.5)*ScintStrip_thick + ilayer*PbStrip_thick; //z of center of scint. strip
    G4double zlead_layer = zscint_layer + 0.5*(ScintStrip_thick + PbStrip_thick );
    G4double Wlayer = Winner_LAC + (zscint_layer+D_LAC/2.0)*(W_LAC - Winner_LAC)/D_LAC; //projective geometry
    G4double Hlayer = Hinner_LAC + (zscint_layer+D_LAC/2.0)*(H_LAC - Hinner_LAC)/D_LAC; //projective geometry

    LeadSheet_boxname.Form("LACleadbox_layer%d", ilayer);
    LeadSheet_logname.Form("LACleadlog_layer%d", ilayer);
    LeadSheet_physname.Form("LACleadphys_layer%d", ilayer );
    
    G4Box *leadbox_temp = new G4Box(LeadSheet_boxname.Data(), Wlayer/2.0, Hlayer/2.0, PbStrip_thick/2.0 );
    G4LogicalVolume *leadlog_temp = new G4LogicalVolume( leadbox_temp, GetMaterial("Lead"), LeadSheet_logname.Data() );

    if( ilayer+1 < Nlayers_total ){
      new G4PVPlacement( 0, G4ThreeVector(0,0,zlead_layer), leadlog_temp, LeadSheet_physname.Data(), log_LAC, false, 0 ); //no copy number since not sensitive
    }
    
    G4double stripwidth,striplength;
    
    if( ilayer % 2 == 0 ){ //"short" (horizontal) strips
      G4double stripW = Wlayer;
      G4double stripH = Hlayer/G4double(NstripsY);
      ScintStrip_boxname.Form( "LACscintbox_layer%d", ilayer );
      //All strips in a given layer have the same dimensions:
      G4Box *strip_temp = new G4Box( ScintStrip_boxname.Data(),
				     stripW/2.0, stripH/2.0, ScintStrip_thick/2.0 );

      ScintStrip_logname.Form( "LACscintlog_layer%d", ilayer );
      G4LogicalVolume *log_strip_temp = new G4LogicalVolume( strip_temp, GetMaterial("NE110A"), ScintStrip_logname.Data() );  

      log_strip_temp->SetVisAttributes( LACscint_visatt );      
      log_strip_temp->SetSensitiveDetector( LACScintSD );
      
      for( G4int istrip=0; istrip<NstripsY; istrip++ ){
	ScintStrip_physname.Form( "LACscintphys_layer%d_strip%d", ilayer, istrip );
	G4ThreeVector strip_pos( 0.0, -Hlayer/2.0 + (istrip+0.5)*stripH, zscint_layer );

	G4int stack_layer = (ilayer < Nlayers_inner) ? 0 : 1;
	
	istack = istrip + (NstripsX+NstripsY)*stack_layer;
	
	new G4PVPlacement( 0, strip_pos, log_strip_temp, ScintStrip_physname.Data(), log_LAC, false, istack );
	
	//Now initialize the detector map quantities:
	(LACScintSD->detmap).Plane[istack] = stack_layer;
	(LACScintSD->detmap).Wire[istack]  = ilayer;
	(LACScintSD->detmap).Row[istack]   = istrip;
	(LACScintSD->detmap).Col[istack]   = -1;
	(LACScintSD->detmap).LocalCoord[istack] = strip_pos;

	//icopy_LACscint++;
      }
    } else { //"long" (vertical) strips
      G4double stripW = Wlayer/G4double(NstripsX);
      G4double stripH = Hlayer;
      ScintStrip_boxname.Form( "LACscintbox_layer%d", ilayer );
      G4Box *strip_temp = new G4Box( ScintStrip_boxname.Data(),
				     stripW/2.0, stripH/2.0, ScintStrip_thick/2.0 );

      ScintStrip_logname.Form( "LACscintlog_layer%d", ilayer );
      G4LogicalVolume *log_strip_temp = new G4LogicalVolume( strip_temp, GetMaterial("NE110A"), ScintStrip_logname.Data() );

      log_strip_temp->SetVisAttributes( LACscint_visatt );
      log_strip_temp->SetSensitiveDetector( LACScintSD );

      for( G4int istrip = 0; istrip<NstripsX; istrip++ ){
	ScintStrip_physname.Form( "LACscintphys_layer%d_strip%d", ilayer, istrip );
	G4ThreeVector strip_pos( -Wlayer/2.0 + (istrip+0.5)*stripW, 0.0, zscint_layer );

	G4int stack_layer = (ilayer < Nlayers_inner) ? 0 : 1;
	
	istack = istrip + NstripsY + (NstripsX+NstripsY)*stack_layer;
	
	new G4PVPlacement( 0, strip_pos, log_strip_temp, ScintStrip_physname.Data(), log_LAC, false, istack );

	//Now initialize the detector map quantities:
	(LACScintSD->detmap).Plane[istack] = stack_layer;
	(LACScintSD->detmap).Wire[istack]  = ilayer;
	(LACScintSD->detmap).Row[istack]   = -1;
	(LACScintSD->detmap).Col[istack]   = istrip;
	(LACScintSD->detmap).LocalCoord[istack] = strip_pos;

	//icopy_LACscint++;
      } 
    }
  }

  G4SBSDetMap dtemp = (LACScintSD->detmap);

  for( map<G4int,G4int>::iterator it=(dtemp.Plane).begin(); it != (dtemp.Plane).end(); ++it ){
    G4int istack = it->first;

    G4int irow = dtemp.Row[istack];
    G4int icol = dtemp.Col[istack];
    G4int iplane = dtemp.Plane[istack];

    G4ThreeVector Rtemp = dtemp.LocalCoord[istack];

    currentline.Form( "  %10d, %10d, %10d, %10d, %18.3f, %18.3f, %18.3f",
		      istack, irow, icol, iplane, Rtemp.x()/cm, Rtemp.y()/cm, Rtemp.z()/cm );
    mapfile << currentline << endl;
    
  }

  mapfile.close();
  
}

void G4SBSHArmBuilder::SetFPP_CH2thick( int ifpp, double CH2thick ){
  double fthickmin = 0.0*cm;
  double fthickmax = 150.0*cm;

  ifpp = ifpp >= 1 ? ( ifpp <= 2 ? ifpp : 2 ) : 1;
  
  fCH2thickFPP[ifpp-1] = CH2thick > fthickmin ? ( CH2thick < fthickmax ? CH2thick : fthickmax ) : fthickmin;
  
}

void G4SBSHArmBuilder::MakePolarimeterGEnRP(G4LogicalVolume *worldlog)
{
  // ----------------------------------------------------------------------------------------------------------------------
  // Make GEn-RP Polarimeter Box
  // Local coordinate system -- x, y relative to central ray and z positions relative to back face of rear field clamp
  // ----------------------------------------------------------------------------------------------------------------------
  
  G4RotationMatrix *SBS_FPP_rm = new G4RotationMatrix;
  SBS_FPP_rm->rotateY( f48D48ang );
  
  double sbsdepth  = 3.0*m; 
  double sbswidth  = 3.0*m; 
  double sbsheight = 2.1*m; 
  double sbsr      = f48D48dist + f48D48depth + (11.43+5.91)*2.54 *cm + sbsdepth/2.0;    // (11.43+5.91)*2.54 *cm from rear clamp code above 
  
  G4Box*           sbsbox = new G4Box("sbsbox", sbswidth/2.0, sbsheight/2.0, sbsdepth/2.0 );
  G4LogicalVolume* sbslog = new G4LogicalVolume(sbsbox, GetMaterial("Air"), "sbslog");
  sbslog->SetVisAttributes( G4VisAttributes::GetInvisible() );
  
  new G4PVPlacement(SBS_FPP_rm, G4ThreeVector(-sbsr*sin(f48D48ang), 0.0, sbsr*cos(f48D48ang) ), 
		    sbslog, "sbsphys", worldlog, false, 0, false);

  // -----------------------------------------------------------------------

  G4RotationMatrix* rot_I = new G4RotationMatrix;
  double        detoffset = -sbsdepth/2.0; 
  G4ThreeVector       pos = G4ThreeVector( 0.0, 0.0, detoffset);

  G4SBSTrackerBuilder trackerbuilder(fDetCon);
  vector<double> gemz, gemw, gemh;  
  
  G4SDManager *sdman = fDetCon->fSDman;
  
  G4double ethresh_default    = 0.05*MeV; // default threshold for scintillator detectors
  G4double timewindow_default = 30.0*ns;  // default time window for scintillator detectors

  // ----------------------------------------------------------------------------------------------------------------------
  // Make Cu analyzer
  // ----------------------------------------------------------------------------------------------------------------------

  double cuanadist   = 60.0*cm; // distance between back face of rear field clamp and front face of Cu analyzer
  
  double cuanaheight = 198.12*cm; 
  double cuanawidth  = 60.96*cm;
  //double cuanadepth  = 4.0*cm; 
  double cuanadepth = 8.89*cm;
  
  G4ThreeVector cuana_pos = pos + G4ThreeVector( 0.0, 0.0, (cuanadist + cuanadepth/2.0) ); 
  
  if( fGEnRP_analyzer_option >= 2 && fDetCon->fExpType != G4SBS::kGEN) {
    G4Box*           cuanabox  = new G4Box("cuanabox", cuanawidth/2.0, cuanaheight/2.0, cuanadepth/2.0 );
    //G4LogicalVolume* cuanalog  = new G4LogicalVolume(cuanabox, GetMaterial("CH2"), "cuanalog");
    //G4LogicalVolume* cuanalog  = new G4LogicalVolume(cuanabox, GetMaterial("Copper"), "cuanalog");

    //Use generic stainless steel for now:
    G4LogicalVolume* cuanalog  = new G4LogicalVolume(cuanabox, GetMaterial("Steel"), "cuanalog");
    //AJRP: Did we intend for the "Copper" analyzer to have CH2 as the material?
    
    new G4PVPlacement(0, cuana_pos, cuanalog,"cuanaphys", sbslog, false, 0, false);
    
    G4VisAttributes *cuanavisatt = new G4VisAttributes( G4Colour(0.0, 0.0, 1.0) );
    cuanavisatt->SetForceWireframe(true);
    cuanalog->SetVisAttributes( cuanavisatt );

    fDetCon->InsertAnalyzerVolume( cuanalog->GetName() );
  }
  
  // ----------------------------------------------------------------------------------------------------------------------
  // Make CE Polarimeter Trackers (2 INFN + 2 UVa before Cu analyzer, 4 UVa after)
  // ---------------------------------------------------------------------------------------------------------------------- 
  
  // int ntracker_ce      = 3; 
  // int ngem_ce[3]       = {2,2,4};
  // double cegem_spacing = 10.0 *cm;

  //SBS tracking geometry notes:
  // 1. SBS detector box starts just downstream of rear clamp; i.e. front edge of SBS box starts at the back edge of rear clamp, which according to "sbsr" above is 11.43+5.91 inches downstream of the back edge of the magnet.
  // 2. according to the GEP detector stack drawings (I know, I know, not the same as GEN), the first SBS tracking layer (INFN) starts 55.91 inches downstream of the CENTER of the SBS magnet, which is 37 cm downstream of the start of the SBS box:
  
  if( fDetCon->fExpType == G4SBS::kGEN ){ //build ONE tracker, not two. 
    int ntracker = 1;
    int ngem = 8;

    G4String SDname = "Harm/SBSGEM";

    vector<double> zplanes_SBStracker(ngem), wplanes_SBStracker(ngem), hplanes_SBStracker(ngem);

    for( int igem=0; igem<ngem; igem++ ){
      double ztemp;
      if( igem < 4 ){
	ztemp = igem*10.0*cm;
      } else {
	ztemp = (igem+1)*10.0*cm;
      }
	
      zplanes_SBStracker[igem] = detoffset + 37.0*cm + ztemp;
      if( igem<2 ){
	wplanes_SBStracker[igem] = 40.0*cm;
	hplanes_SBStracker[igem] = 150.0*cm;
      } else {
	wplanes_SBStracker[igem] = 60.0*cm;
	hplanes_SBStracker[igem] = 200.0*cm;
      }
    }

    trackerbuilder.BuildComponent( sbslog, rot_I, G4ThreeVector(0,0,0), ngem, zplanes_SBStracker, wplanes_SBStracker, hplanes_SBStracker, SDname );
    
  } else {
  
    int ntracker_ce      = 2; 
    int ngem_ce[2]       = {4,4};
    double cegem_spacing = 10.0 *cm;
  
    vector<G4String> SDnames_ce; 
    //SDnames_ce.push_back("Harm/CEPolFront1");
    //SDnames_ce.push_back("Harm/CEPolFront2");
    //SDnames_ce.push_back("Harm/CEPolRear");

    // SDnames_ce.push_back("Harm/CEPolFront1");
    SDnames_ce.push_back("Harm/CEPolFront");
    SDnames_ce.push_back("Harm/CEPolRear");
  
    for( int i = 0; i<ntracker_ce; i++){ 
      gemz.resize( ngem_ce[i] ); 
      gemw.resize( ngem_ce[i] );
      gemh.resize( ngem_ce[i] );
      for( int j = 0; j < ngem_ce[i]; j++ ){
	if( i == 0 ){
	  gemz[j] = -cuanadepth/2.0 + (double)(j-4)*cegem_spacing;
	  if( j < 2 ){
	    gemw[j] = 40.0*cm;
	    gemh[j] = 150.0*cm;
	  } else {
	    gemw[j] = 60.0*cm;
	    gemh[j] = 200.0*cm;
	  }
	} else {// else if( i == 1 ){
	  // 	gemz[j] = -cuanadepth/2.0 + (double)(j-2)*cegem_spacing;
	  // 	gemw[j] = 60.0*cm;
	  // 	gemh[j] = 200.0*cm;
	  // }
	  //else {
	  gemz[j] = cuanadepth/2.0 + (double)(j+1)*cegem_spacing;
	  gemw[j] = 60.0*cm;
	  gemh[j] = 200.0*cm;
	}
      }

      bool ispol = i > 0;
      
      trackerbuilder.BuildComponent( sbslog, rot_I, cuana_pos, ngem_ce[i], gemz, gemw, gemh, SDnames_ce[i], ispol );
    }
  }
  
  // ----------------------------------------------------------------------------------------------------------------------
  // Make Active Analyzer
  // ----------------------------------------------------------------------------------------------------------------------

  double actanadist = 120.0 *cm; // distance between back face of rear field clamp and front face of active analyzer

  int    nactanabarsx   = 4;
  int    nactanabarsz   = 8;
  double actanabarwidth = 4.0 *cm;
  double actanabardepth = 4.0 *cm;
  double actanaheight   = 25.0 *cm;

  double actanawidth = actanabarwidth * nactanabarsx;
  double actanadepth = actanabardepth * nactanabarsz;

  if( fGEnRP_analyzer_option == 4 ) {  // CGEN configuration
    nactanabarsx   = 2;
    nactanabarsz   = 1;
    actanabarwidth = 10.0 *cm;
    actanabardepth = 10.0 *cm;
    actanaheight   = 100.0 *cm;

    actanawidth = 60.0 *cm;
    actanadepth = actanabardepth * nactanabarsz;
  }    
  if( fDetCon->fExpType != G4SBS::kGEN) {  
    G4ThreeVector actana_pos = pos + G4ThreeVector( 0.0, 0.0, (actanadist + actanadepth/2.0) );
    
    G4Box*           actanabox  = new G4Box("actanabox", actanawidth/2.0, actanaheight/2.0, actanadepth/2.0 );
    G4LogicalVolume* actanalog  = new G4LogicalVolume(actanabox, GetMaterial("Air"), "actanalog");
    
    if( fGEnRP_analyzer_option == 2 ) {   // long axis parallel to neutron direction
      G4RotationMatrix* rot_ana  = new G4RotationMatrix;
      rot_ana->rotateX( 90.0 *deg );
      new G4PVPlacement(rot_ana, actana_pos, actanalog,"actanaphys", sbslog, false, 0, false);  
    }
    else if( fGEnRP_analyzer_option == 3 || fGEnRP_analyzer_option == 4 )  // long axis perpendicular to neutron direction or CGEN configuration
      new G4PVPlacement(0, actana_pos, actanalog,"actanaphys", sbslog, false, 0, false);  
    
    actanalog ->SetVisAttributes( G4VisAttributes::GetInvisible() );
    
    G4String ActAnScintSDname   = "Harm/ActAnScint";
    G4String ActAnScintcollname = "ActAnScintHitsCollection";
    G4SBSCalSD *ActAnScintSD    = NULL;
    
    G4Box *actanabarbox = new G4Box("actanabarbox", actanabarwidth/2.0, actanaheight/2.0, actanabardepth/2.0 );
    G4LogicalVolume *actanabarlog = new G4LogicalVolume( actanabarbox, GetMaterial("BBHodo_Scinti"), "actanabarlog" );
    
    
    if( !((G4SBSCalSD*) sdman->FindSensitiveDetector(ActAnScintSDname)) ) {
      G4cout << "Adding SBS PR Polarimeter Scint BeamSide Sensitive Detector to SDman..." << G4endl;
      ActAnScintSD = new G4SBSCalSD( ActAnScintSDname, ActAnScintcollname );
      
      sdman->AddNewDetector( ActAnScintSD );
      (fDetCon->SDlist).insert( ActAnScintSDname );
      fDetCon->SDtype[ActAnScintSDname] = G4SBS::kCAL;
      (ActAnScintSD->detmap).depth = 0;
    
      fDetCon->SetThresholdTimeWindowAndNTimeBins( ActAnScintSDname, ethresh_default, timewindow_default, 25 );
      fDetCon->InsertSDboundaryVolume( actanalog->GetName(), ActAnScintSDname );
    }
    G4VisAttributes *actanavisatt = new G4VisAttributes( G4Colour(0.0, 1.0, 0.0));
    actanavisatt->SetForceWireframe(true);
    actanabarlog->SetVisAttributes( actanavisatt );
    actanabarlog->SetSensitiveDetector( ActAnScintSD ); 
    
    fDetCon->InsertAnalyzerVolume( actanabarlog->GetName() );
  
    int barindex = 0; // SD detectors indexing (left to right, bottom to top)
    for(int ix = 0; ix < nactanabarsx; ix++) {
      double xbar = (((double)(ix-nactanabarsx/2.0) * actanabarwidth) + actanabarwidth/2.0);
      for(int iz = 0; iz < nactanabarsz; iz++) {
	double zbar = (((double)(iz-nactanabarsz/2.0) * actanabardepth) + actanabardepth/2.0);
	if( fGEnRP_analyzer_option == 2 || fGEnRP_analyzer_option == 3 ) // Glasgow
	  new G4PVPlacement( 0, G4ThreeVector(xbar, 0, zbar), actanabarlog, "actanabarphys", actanalog, false, barindex );
	else if( fGEnRP_analyzer_option == 4 ) { // CGEN 
	  xbar = (2*ix-1)*(actanawidth/2.0) + (1-2*ix)*actanabarwidth/2.0;
	  new G4PVPlacement( 0, G4ThreeVector(xbar, 0, 0), actanabarlog, "actanabarphys", actanalog, false, barindex );
	}      
	barindex++;
      }
    }
  }
  
  if( fGEnRP_analyzer_option == 2 ) // long axis parallel to neutron direction (needed to maintain correct centering)
    actanadepth = actanaheight;
  
  // ----------------------------------------------------------------------------------------------------------------------
  // Make PR Polarimeter Trackers (2 UVa each)
  // ----------------------------------------------------------------------------------------------------------------------
  
  double prgem_perpdist = 50.0 *cm;  // pependicular distance to first GEM from analyzer centre
  double prgem_spacing  = 10.0 *cm;

  gemz.resize( 2 );
  gemw.resize( 2 );
  gemh.resize( 2 );
  
  // -----------------------------------------------------------------------
  // beamline side gem detectors

  G4RotationMatrix* rot_bs  = new G4RotationMatrix;
  rot_bs->rotateY( 90.0 *deg );
  
  G4ThreeVector prgem_posbs = pos + G4ThreeVector( prgem_perpdist, 0.0, (actanadist + actanadepth/2.) );

  if( fGEnRP_analyzer_option != 0 && fDetCon->fExpType != G4SBS::kGEN) {  
    for( int j = 0; j < 2; j++ ){
      gemz[j] = (double)(j-1)*prgem_spacing; 
      gemw[j] = 60.0*cm;
      gemh[j] = 200.0*cm;
    }
    trackerbuilder.BuildComponent( sbslog, rot_bs, prgem_posbs, 2, gemz, gemw, gemh, G4String("Harm/PRPolGEMBeamSide"));
  }

  // -----------------------------------------------------------------------
  // far side gem detectors

  G4RotationMatrix* rot_fs  = new G4RotationMatrix;
  rot_fs->rotateY( -90.0 *deg );
  
  G4ThreeVector prgem_posfs = pos + G4ThreeVector( -prgem_perpdist, 0.0, (actanadist + actanadepth/2.0) );
  
  if( fDetCon->fExpType != G4SBS::kGEN) {  
    for( int j = 0; j < 2; j++ ){
      gemz[j] = (double)(j-1)*prgem_spacing; 
      gemw[j] = 60.0*cm;
      gemh[j] = 200.0*cm;
    }
    trackerbuilder.BuildComponent( sbslog, rot_fs, prgem_posfs, 2, gemz, gemw, gemh, G4String("Harm/PRPolGEMFarSide"));
  }
  
  // ----------------------------------------------------------------------------------------------------------------------
  // Make PR Polarimeter Hodoscopes (1 each)
  // ----------------------------------------------------------------------------------------------------------------------

  double prscint_perpdist = 75 *cm; // pependicular distance to front face of scintillator plane from analyzer centre 

  int    nprbars          = 24;
  double prbarheight      = 8.6 *cm;
  double prscintheight    = nprbars*prbarheight;
  double prscintwidth     = 50.0 *cm; 
  double prscintbsdepth   = 3.0 *mm;  // put thin dE detectors on beamline side
  double prscintfsdepth   = 30.0 *mm; // put thicker E detectors on far side 

  // -----------------------------------------------------------------------
  // beamline side scintillator detectors

  G4ThreeVector prscintbs_pos = pos + G4ThreeVector( (prscint_perpdist + prscintbsdepth/2.0), 0.0, (actanadist + actanadepth/2.0) );

  if(fDetCon->fExpType != G4SBS::kGEN) {  
    if( fGEnRP_analyzer_option != 0){ 
      
      G4Box *prscintbsbox = new G4Box("prscintbsbox", prscintwidth/2.0, prscintheight/2.0, prscintbsdepth/2.0 );
      
      G4LogicalVolume *prscintbslog = new G4LogicalVolume( prscintbsbox, GetMaterial("Air"), "prscintbslog" );
      
      new G4PVPlacement( rot_bs, prscintbs_pos, prscintbslog,"prscintbsphys", sbslog, false, 0, false);
      prscintbslog->SetVisAttributes( G4VisAttributes::GetInvisible() );
      
      G4Box *prbarbsbox = new G4Box("prbarbsbox", prscintwidth/2.0, prbarheight/2.0, prscintbsdepth/2.0 );
      G4LogicalVolume *prbarbslog = new G4LogicalVolume( prbarbsbox, GetMaterial("BBHodo_Scinti"), "prbarbslog" );
      
      G4String PRPolScintBSSDname   = "Harm/PRPolScintBeamSide";
      G4String PRPolScintBScollname = "PRPolScintBSHitsCollection";
      G4SBSCalSD *PRPolScintBSSD    = NULL;
      
      if( !((G4SBSCalSD*) sdman->FindSensitiveDetector(PRPolScintBSSDname)) ) {
	G4cout << "Adding SBS PR Polarimeter Scint BeamSide Sensitive Detector to SDman..." << G4endl;
	PRPolScintBSSD = new G4SBSCalSD( PRPolScintBSSDname, PRPolScintBScollname );
	
	sdman->AddNewDetector( PRPolScintBSSD );
	(fDetCon->SDlist).insert( PRPolScintBSSDname );
	fDetCon->SDtype[PRPolScintBSSDname] = G4SBS::kCAL;
	(PRPolScintBSSD->detmap).depth = 0;
	
	fDetCon->SetThresholdTimeWindowAndNTimeBins( PRPolScintBSSDname, ethresh_default, timewindow_default, 25 );
	fDetCon->InsertSDboundaryVolume( prscintbslog->GetName(), PRPolScintBSSDname );
      }
      prbarbslog->SetVisAttributes(G4Colour(0.0, 1.0, 0.0));
      prbarbslog->SetSensitiveDetector( PRPolScintBSSD ); 
      
      // SD detectors indexed from bottom to top
      for(int i = 0; i < nprbars; i++) {
	double ybar = (((double)(i-nprbars/2.0) * prbarheight) + prbarheight/2.0);
	new G4PVPlacement( 0, G4ThreeVector(0, ybar, 0), prbarbslog, "prbarbsphys", prscintbslog, false, i );
      }
    }
    
    // -----------------------------------------------------------------------
    // far side scintillator detectors
    
    G4ThreeVector prscintfs_pos = pos + G4ThreeVector( -(prscint_perpdist + prscintfsdepth/2.0), 0.0, (actanadist + actanadepth/2.0) );
    
    G4Box *prscintfsbox = new G4Box("prscintfsbox", prscintwidth/2.0, prscintheight/2.0, prscintfsdepth/2.0 );
    
    G4LogicalVolume *prscintfslog = new G4LogicalVolume( prscintfsbox, GetMaterial("Air"), "prscintfslog" );
  
    G4Box *prbarfsbox = new G4Box("prbarfsbox", prscintwidth/2.0, prbarheight/2.0, prscintfsdepth/2.0 );
    G4LogicalVolume *prbarfslog = new G4LogicalVolume( prbarfsbox, GetMaterial("BBHodo_Scinti"), "prbarfslog" );
    
    new G4PVPlacement( rot_fs, prscintfs_pos, prscintfslog,"prscintfsphys", sbslog, false, 0, false);
    prscintfslog->SetVisAttributes( G4VisAttributes::GetInvisible() );
    
    G4String PRPolScintFSSDname   = "Harm/PRPolScintFarSide";
    G4String PRPolScintFScollname = "PRPolScintFSHitsCollection";
    G4SBSCalSD *PRPolScintFSSD    = NULL;
    
    if( !((G4SBSCalSD*) sdman->FindSensitiveDetector(PRPolScintFSSDname)) ) {
      G4cout << "Adding SBS PR Polarimeter Scint Far Side Sensitive Detector to SDman..." << G4endl;
      PRPolScintFSSD = new G4SBSCalSD( PRPolScintFSSDname, PRPolScintFScollname );
      
      sdman->AddNewDetector( PRPolScintFSSD );
      (fDetCon->SDlist).insert( PRPolScintFSSDname );
      fDetCon->SDtype[PRPolScintFSSDname] = G4SBS::kCAL;
      (PRPolScintFSSD->detmap).depth = 0;
      
      fDetCon->SetThresholdTimeWindowAndNTimeBins( PRPolScintFSSDname, ethresh_default, timewindow_default, 25 );
      fDetCon->InsertSDboundaryVolume( prscintfslog->GetName(), PRPolScintFSSDname );
    }
    prbarfslog->SetVisAttributes(G4Colour(0.0, 1.0, 0.0));
    prbarfslog->SetSensitiveDetector( PRPolScintFSSD ); 
    
    // SD detectors indexed from bottom to top 
    for(int i = 0; i < nprbars; i++) {
      double ybar = (((double)(i-nprbars/2.0) * prbarheight) + prbarheight/2.0);
      new G4PVPlacement( 0, G4ThreeVector(0, ybar, 0), prbarfslog, "prbarfsphys", prscintfslog, false, i );
    }
  }

}

void G4SBSHArmBuilder::MakeNeutronVeto(G4LogicalVolume* worldlog, G4double dist_from_hcal)
{
  bool checkOL = false;
  //build a "simple" veto with the old BIG HAND veto elements: 
  //dimensions x, y, z = 180*15*5 cm^3 (or is it actually 72*6*2 ci ?) 
  // 24 needed to cover HCAL (maybe 25 for extra coverage ?)
  G4double VetoElemWidth = 72.0*2.54*cm;
  G4double VetoElemHeight = 6.0*2.54*cm;
  G4double VetoElemDepth = 2.0*2.54*cm;
  G4double planes_interdist = 5.4*cm;//20.0*cm;
  G4int Nplanes = 1;
  G4int Nelem = 24;
 
  G4RotationMatrix *HArmRot = new G4RotationMatrix;
  HArmRot->rotateY(f48D48ang);
  
  G4Box* VetoMotherBox = new G4Box("VetoMotherBox", VetoElemWidth/2.0, Nelem*VetoElemHeight/2.0, (Nplanes*VetoElemDepth+(Nplanes-1)*planes_interdist)/2.0);
  G4LogicalVolume* VetoMotherLog = new G4LogicalVolume(VetoMotherBox, GetMaterial("Air"), "VetoMotherLog");

  G4ThreeVector Veto_pos( ( fHCALdist-dist_from_hcal+VetoMotherBox->GetZHalfLength() ) * sin( -f48D48ang ),  fHCALvertical_offset, ( fHCALdist-dist_from_hcal+VetoMotherBox->GetZHalfLength() ) * cos( -f48D48ang ) );
      
  new G4PVPlacement(HArmRot, Veto_pos, VetoMotherLog, "VetoMother_phys", worldlog, false, 0, checkOL);
      
  G4double mylarthickness = 0.0020*cm, airthickness = 0.0040*cm;
  G4double mylar_air_sum = mylarthickness + airthickness;
      
  G4Box* VetoElemBox = new G4Box("VetoElemBox", VetoElemWidth/2.0, VetoElemHeight/2.0, VetoElemDepth/2.0);
  //Let's do the wrap... assuming similar to BB calo blocks
  G4LogicalVolume* VetoElemLog = new G4LogicalVolume(VetoElemBox, GetMaterial("Air"), "VetoElemLog");
      
  G4Box* MylarBoxHollow = new G4Box("MylarBoxHollow", VetoElemWidth/2.0, VetoElemHeight/2.0-mylarthickness, VetoElemDepth/2.0-mylarthickness);
  G4SubtractionSolid* MylarWrap = new G4SubtractionSolid("MylarWrap", VetoElemBox, MylarBoxHollow, 0, G4ThreeVector(0.0, 0.0, 0.0));
      
  G4LogicalVolume *MylarWrapLog = new G4LogicalVolume( MylarWrap, GetMaterial("Mylar"), "MylarWrapLog" );
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), MylarWrapLog, "MylarWrapPhys", VetoElemLog, false, 0, checkOL);
      
  G4Box* VetoScint = new G4Box("VetoScint", VetoElemWidth/2.0, VetoElemHeight/2.0-mylar_air_sum, VetoElemDepth/2.0-mylar_air_sum);
  G4LogicalVolume *VetoScintLog = new G4LogicalVolume( VetoScint, GetMaterial("Acrylic"), "VetoScintLog" );
  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), VetoScintLog, "VetoScintPhys", VetoElemLog, false, 0, checkOL);
      
  G4String sdname = "Harm/Veto_Scint";
  G4String collname = "Veto_ScintHitsCollection";
      
  G4SBSCalSD *veto_scint_sd = NULL;
      
  if( !( veto_scint_sd = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector( sdname ) ) ){
    G4cout << "Adding Veto Scint sensitive detector to sdman..." << G4endl;
    veto_scint_sd = new G4SBSCalSD( sdname, collname );
    fDetCon->fSDman->AddNewDetector( veto_scint_sd );
    (fDetCon->SDlist).insert( sdname );
    fDetCon->SDtype[sdname] = G4SBS::kCAL;
    (veto_scint_sd->detmap).depth = 1;
    VetoScintLog->SetSensitiveDetector( veto_scint_sd );
	
    fDetCon->SetThresholdTimeWindowAndNTimeBins( sdname, 0.0*MeV, 50.0*ns, 25 );
  }
      
  fDetCon->InsertSDboundaryVolume( VetoMotherLog->GetName(), sdname );
      
  for(int j = 0; j<Nplanes; j++){
    for(int i = 0; i<Nelem; i++){
      //G4ThreeVector ElemPos(0.0, ((-G4double(Nelem)+1)/2.0+i)*VetoElemHeight, (j-0.5)*(planes_interdist+VetoElemDepth) );//only works for NPlanes = 2....
      G4ThreeVector ElemPos(0.0, ((-G4double(Nelem)+1)/2.0+i)*VetoElemHeight, 0.0);//only works for NPlanes = 1....
      new G4PVPlacement(0, ElemPos, VetoElemLog, "VetoElemPhys", VetoMotherLog, false, Nelem*j+i, checkOL);
    }
  }
      	
  VetoMotherLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  VetoElemLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  MylarWrapLog->SetVisAttributes(G4VisAttributes::GetInvisible());
  VetoScintLog->SetVisAttributes(G4Colour(0.0, 1.0, 1.0, 1.0));
}
  
//*****************************************

void G4SBSHArmBuilder::MakeTest( G4LogicalVolume *worldlog)
{
  //Scintillator specs
  G4double scin_vert = 50.0*cm;
  G4double scin_width = 60.0*cm;
  G4double scin_depth = 1.5*2.54*cm; 
  G4double space_in_front = 1.0*cm; // space between motherlog boundary and scintillator 

  //Test Mother Volume to house all Scintilator+GEMs+Calorimeter
  G4double Test_vert = scin_vert + 1.0*cm;
  G4double Test_width = scin_width + 1.0*cm;
  G4double Test_depth = 1.0*cm + scin_depth + 78.0*cm + 50.0*cm + 1.0*cm; // 78cm for GEMs & 50cm for Calorimeter (roughly)

  //Geometric info & placement of Test detector
  G4double r_Test = fHCALdist + Test_depth/2.0;
  G4double angle_Test = 360.0*deg - f48D48ang;
  G4ThreeVector R_Test( r_Test*sin( angle_Test ), 0.0, r_Test*cos( angle_Test ));

  G4RotationMatrix *TestRot = new G4RotationMatrix;
  TestRot->rotateY( -angle_Test );

  G4Box *Test_Box = new G4Box("Test_Box", Test_width/2.0, Test_vert/2.0, Test_depth/2.0);
  G4LogicalVolume *Test_Log = new G4LogicalVolume(Test_Box, GetMaterial("Air"), "Test_Log");
  new G4PVPlacement( TestRot, R_Test, Test_Log, "Test_Phys", worldlog, 0, false, 0);

  //Geometric info and placement of Scintillator
  G4Box *scin_Box = new G4Box("scin_Box", scin_width/2.0, scin_vert/2.0, scin_depth/2.0);
  G4LogicalVolume *scin_Log = new G4LogicalVolume(scin_Box, GetMaterial("POLYSTYRENE"), "scin_Log");
  new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, -Test_depth/2.0 + scin_depth/2.0 + space_in_front ), scin_Log, "scin_phys", Test_Log, 0, false, 0);

  //Making GEMs
  G4SBSTrackerBuilder trackerbuilder(fDetCon);
  G4RotationMatrix *GEMRot = new G4RotationMatrix;
  GEMRot->rotateY( 0.0*deg );
  
  //GEM specs
  int ngem = 5; 		// 5 GEMs
  G4double gemz_spacing = 13.0*cm; // spacing between GEMs is 13cm
  G4ThreeVector gem_pos( 0.0, 0.0, -Test_depth/2.0 + scin_depth + 5.0*gemz_spacing);
  vector<double> gemz, gemw, gemh;
  gemz.resize( ngem );
  gemw.resize( ngem );
  gemh.resize( ngem );
  for(int i=0; i<ngem; i++){
    gemz[i] = -Test_depth/2.0 + space_in_front + scin_depth + gemz_spacing + ((double)i)*gemz_spacing;
    gemw[i] = 60.0*cm;
    gemh[i] = 50.0*cm;
}  
  trackerbuilder.BuildComponent( Test_Log, GEMRot, gem_pos, ngem, gemz, gemw, gemh, "Test/HC" );

  //Visual
  Test_Log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  // G4VisAttributes *Testvisatt = new G4VisAttributes( G4Colour( 0.0, 0.0, 1.0 ) );
  // Testvisatt->SetForceWireframe(true);
  // Test_Log->SetVisAttributes( Testvisatt );
  G4VisAttributes *Scinvisatt = new G4VisAttributes( G4Colour( 0.0, 1.0, 0.0 ) );
  scin_Log->SetVisAttributes( Scinvisatt );
}
