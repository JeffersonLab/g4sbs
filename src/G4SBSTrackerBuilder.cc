#include "G4SBSTrackerBuilder.hh"

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
#include "G4RotationMatrix.hh"

#include "G4SBSGEMSD.hh"
#include "sbstypes.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

G4SBSTrackerBuilder::G4SBSTrackerBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){

}

G4SBSTrackerBuilder::~G4SBSTrackerBuilder(){;}

void G4SBSTrackerBuilder::BuildComponent(G4LogicalVolume *){
  // This shouldn't get called
  ;
}

//This routine allows us to flexibly position GEM modules without code duplication:
void G4SBSTrackerBuilder::BuildComponent(G4LogicalVolume *Mother, G4RotationMatrix *rot, G4ThreeVector pos, unsigned int nplanes, vector<double> zplanes, vector<double> wplanes, vector<double> hplanes, G4String SDname ) 
{
  //This routine will create and position a GEM tracker consisting of nplanes planes centered at position pos oriented with rotation rot wrt logical volume Mother. 
  //The list of z coordinates, widths and heights of the planes are passed as arguments:

  G4bool GEMflip = fDetCon->GetFlipGEM();
  
  // How should we interpret the rotation Matrix rot? It is the rotation that orients the z axis of the tracker with respect to the mother volume. 
  // Since pos is the nominal position of the tracker with respect to the mother volume, 
  // the positioning of the centers of the tracker planes should be pos + zplane * tracker_zaxis
  G4ThreeVector zaxis(0,0,1);
  zaxis *= rot->inverse();
  
  G4String TrackerPrefix = SDname; 

  //Create sensitive detector for this tracker:
  G4String GEMSDname = SDname;
  G4String GEMSDname_nopath = SDname;
  GEMSDname_nopath.remove(0,SDname.last('/')+1);
  G4String GEMcolname = GEMSDname_nopath; //We have to remove all the directory structure from the Hits collection name or else GEANT4 SDmanager routines will not handle correctly.
  GEMcolname += "HitsCollection";

  G4SBSGEMSD* GEMSD;

  if( !(GEMSD = (G4SBSGEMSD*) fDetCon->fSDman->FindSensitiveDetector(GEMSDname)) ){ //Make sure SD with this name doesn't already exist
    GEMSD = new G4SBSGEMSD( GEMSDname, GEMcolname );
    fDetCon->fSDman->AddNewDetector(GEMSD);
    (fDetCon->SDlist).insert(GEMSDname);
    fDetCon->SDtype[GEMSDname] = G4SBS::kGEM;
  }

  fDetCon->InsertSDboundaryVolume( Mother->GetName(), GEMSDname );
  
  if( !( nplanes > 0 && zplanes.size() == nplanes && wplanes.size() == nplanes && hplanes.size() == nplanes ) ){
    G4cout << "MakeTracker called with invalid arguments, doing nothing..." << G4endl;
    return;
  }

  //Best way to handle "flipping" of GEMs? Boxes are symmetric, just flip ordering of z positions!
  
  //Define z extent of various GEM layer components:

  const unsigned int nlayers = 24;

  double gempz[nlayers] = {
    120.0*um, 3.0*mm, 120.0*um, 5.0*um, 50.0*um, 3.0*mm, // cover + cathode
    5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM0
    5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM1
    5.0*um, 50.0*um, 5.0*um, 2.0*mm, // GEM2
    10.0*um, 50.0*um, 180.0*um, 120.0*um, 3.0*mm, 120.0*um // readout + honeycomb
  };

  //Define GEM materials:

  G4Material *gempm[nlayers] = {
    GetMaterial("NEMAG10"), GetMaterial("NOMEX"), GetMaterial("NEMAG10"), GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("GEMgas"),
    GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("Copper"), GetMaterial("GEMgas"),
    GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("Copper"), GetMaterial("GEMgas"),
    GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("Copper"), GetMaterial("GEMgas"),
    GetMaterial("Copper"), GetMaterial("Kapton"), GetMaterial("NEMAG10"), GetMaterial("NEMAG10"), GetMaterial("NOMEX"), GetMaterial("NEMAG10")
  };

  unsigned int gidx, gpidx;

  //Determine the largest transverse dimensions of a plane in this tracker:
  double gemmaxw = 0.0;
  double gemmaxh = 0.0;

  for( gidx = 0; gidx < nplanes; gidx++ ){
    if( wplanes[gidx] > gemmaxw ) gemmaxw = wplanes[gidx];
    if( hplanes[gidx] > gemmaxh ) gemmaxh = hplanes[gidx];
  }

  //Determine the total z extent of a plane in this tracker:
  double gempzsum = 0.0;
  for( gpidx = 0; gpidx < nlayers; gpidx++ ){
    gempzsum += gempz[gpidx];
  }

  //implement aluminum shielding: 
  G4bool useAlshield = fDetCon->GetGEMuseAlshield();

  // G4Box *alshield_box;
  // G4LogicalVolume *alshield_log;
  
  G4double Al_thick=0.0, Airgap_thick=0.0;
  if( useAlshield ){
    //G4cout << "building aluminum shielding for GEMs..."; 
    
    Al_thick = fDetCon->GetGEMAlShieldThick();
    Airgap_thick = fDetCon->GetGEMAirGapThick();
    

    //G4cout << " complete" << G4endl;
  }

  // char gpname[20][50][3][255];
  // char gemname[10][255], gemboxname[10][255], gemlogname[10][255];
  char cgidx[255];
  char cgpidx[255];

  G4Box *gpbox;
  G4LogicalVolume *gplog;

  G4VisAttributes *gemvisatt = new G4VisAttributes( G4Colour( 1.0, 0.0, 0.0 ) );

  gemvisatt->SetForceWireframe(true);

  G4VisAttributes *gemsdvisatt = new G4VisAttributes( G4Colour( 1.0, 0.0, 0.0 ) );

  
  
  for( gidx = 0; gidx < nplanes; gidx++ ){
    sprintf( cgidx, "%02d", gidx );

    // sprintf(gemboxname[gidx], "gembox_%02d", gidx);
    // sprintf(gemlogname[gidx], "gemlog_%02d", gidx);
    G4String gemboxname =  TrackerPrefix + G4String("_gembox_");
    gemboxname += cgidx;
    G4String gemlogname = TrackerPrefix + G4String("_gemlog_");
    gemlogname += cgidx;

    G4Box *gembox = new G4Box( gemboxname, wplanes[gidx]/2.0, hplanes[gidx]/2.0, gempzsum/2.0 );
    G4LogicalVolume *gemlog = new G4LogicalVolume( gembox, GetMaterial("Air"), gemlogname, 0, 0, 0 );

    gemlog->SetVisAttributes( gemvisatt );

    double ztemp = 0.0;
    double sign = 1.0;
    if( GEMflip ) { //this reverses ordering of GEM materials front-to-back!
      ztemp = gempzsum;
      sign = -1.0;
    }

    if( useAlshield ){ //put Aluminum shielding upstream of GEM:
      G4String shieldboxname = TrackerPrefix + G4String("alshield_gem") + cgidx + G4String("box");
      G4String shieldlogname = TrackerPrefix + G4String("alshield_gem") + cgidx + G4String("log");
      
      G4String frontshieldphysname = TrackerPrefix + G4String("alshieldfront_gem") + cgidx + G4String("phys");
      G4String backshieldphysname = TrackerPrefix + G4String("alshieldback_gem") + cgidx + G4String("phys");

      G4cout << "Building and placing GEM aluminum shielding, plane = " << gidx << G4endl;

      G4Box *alshield_box = new G4Box( shieldboxname, wplanes[gidx]/2.0, hplanes[gidx]/2.0, Al_thick/2.0 );
      G4LogicalVolume *alshield_log = new G4LogicalVolume( alshield_box, GetMaterial("Aluminum"), shieldlogname, 0, 0, 0 );

      G4VisAttributes *alshield_visatt = new G4VisAttributes( G4Colour( 0.6, 0.6, 0.6 ) );
      //alshield_visatt->SetForceWireframe( true );
      
      alshield_log->SetVisAttributes( alshield_visatt );
      
      //position the front shielding upstream of the GEM with possible air gap:
      G4ThreeVector frontshieldpos = pos + (zplanes[gidx] - gempzsum/2.0 - Airgap_thick - Al_thick/2.0) * zaxis;
      
      new G4PVPlacement( 0, frontshieldpos, alshield_log, frontshieldphysname, Mother, true, gidx+1, false );

      //position the back shielding immediately downstream of the GEM with no air gap:
      G4ThreeVector backshieldpos = frontshieldpos + (gempzsum + Al_thick + Airgap_thick ) * zaxis;

      new G4PVPlacement( 0, backshieldpos, alshield_log, backshieldphysname, Mother, true, gidx+1, false );
    }

    double zdrift; //z position of center of drift region relative to GEM plane box:
    
    for( gpidx = 0; gpidx < nlayers; gpidx++ ){
      sprintf( cgpidx, "_%02d_%03d_", gidx, gpidx );
      
      G4String gempboxname = TrackerPrefix + G4String("_gemplane") + cgpidx + G4String("box");
      G4String gemplogname = TrackerPrefix + G4String("_gemplane") + cgpidx + G4String("log");
      G4String gempphysname = TrackerPrefix + G4String("_gemplane") + cgpidx + G4String("phy");

      ztemp += sign*gempz[gpidx]/2.0;

      gpbox = new G4Box( gempboxname, wplanes[gidx]/2.0, hplanes[gidx]/2.0, gempz[gpidx]/2.0 );
      gplog = new G4LogicalVolume( gpbox, gempm[gpidx], gemplogname, 0, 0, 0 );

      new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, ztemp - gempzsum/2.0 ), gplog, gempphysname, gemlog, false, 0, false ); 

      

      //Assign sensitive volume: why 5?  // SPR: This is the gas drift region
      if( gpidx == 5 ){
	gplog->SetSensitiveDetector(GEMSD);
	gplog->SetVisAttributes( gemsdvisatt );
	// Until we implement actual strips/wires in the GEM construction, the detmap is irrelevant for the GEMs
	zdrift = ztemp - gempzsum/2.0;
      } else {
	gplog->SetVisAttributes( G4VisAttributes::Invisible );
      }

      ztemp += sign*gempz[gpidx]/2.0;
      
    }

    G4String gemname = TrackerPrefix + G4String("_gemphys_") + cgidx;
    //Now place the fully constructed GEM plane in the mother logical volume:
    //G4ThreeVector plane_pos = pos + G4ThreeVector( 0.0, 0.0, zplanes[gidx] );
    G4ThreeVector plane_pos = pos + zplanes[gidx] * zaxis;
    //Now we are positioning the GEM AFTER positioning all of its components inside its logical volume:
    if( gidx == 0 ){ //set Z offset of GEMSD to equal the center of the drift region of first GEM
      GEMSD->SetZoffset( plane_pos.getZ() + zdrift );
    }
    new G4PVPlacement( rot, plane_pos, gemlog, gemname, Mother, true, gidx+1, false );
  }
}


