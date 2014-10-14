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

G4SBSTrackerBuilder::G4SBSTrackerBuilder(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){

}

G4SBSTrackerBuilder::~G4SBSTrackerBuilder(){;}

void G4SBSTrackerBuilder::BuildComponent(G4LogicalVolume *){
  // This shouldn't get called
  ;
}

//This routine allows us to flexibly position GEM modules without code duplication:
void G4SBSTrackerBuilder::BuildComponent(G4LogicalVolume *Mother, G4RotationMatrix *rot, G4ThreeVector pos, unsigned int nplanes, vector<double> zplanes, vector<double> wplanes, vector<double> hplanes, G4int TrackerID=0) 
{
  //This routine will create and position a GEM tracker consisting of nplanes planes centered at position pos oriented with rotation rot wrt logical volume Mother. 
  //The list of z coordinates, widths and heights of the planes are passed as arguments:

  //    G4String MotherName = Mother->GetName();

  //  G4SDManager* fSDman = G4SDManager::GetSDMpointer();

  // How should we interpret the rotation Matrix rot? It is the rotation that orients the z axis of the tracker with respect to the mother volume. Since pos is the nominal position of the tracker with respect to the mother volume, 
  // the positioning of the centers of the tracker planes should be pos + zplane * tracker_zaxis
  G4ThreeVector zaxis(0,0,1);
  zaxis *= rot->inverse();

  char ctrid[20];
  sprintf( ctrid, "%02d", TrackerID );
  
  G4String TrackerPrefix = G4String("Tracker_") + ctrid; 

  //Define sensitive detector, if not already done:
  G4String GEMSDname = "G4SBS/GEM";
  G4String GEMcolname = "GEMcol";
  G4SBSGEMSD* GEMSD;

  if( !(GEMSD = (G4SBSGEMSD*) fDetCon->fSDman->FindSensitiveDetector(GEMSDname)) ){
    GEMSD = new G4SBSGEMSD( GEMSDname, GEMcolname );
    fDetCon->fSDman->AddNewDetector(GEMSD);
    fDetCon->SDlist[GEMSDname] = GEMSD;
  }

  if( !( nplanes > 0 && zplanes.size() == nplanes && wplanes.size() == nplanes && hplanes.size() == nplanes ) ){
    G4cout << "MakeTracker called with invalid arguments, doing nothing..." << G4endl;
    return;
  }

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
    for( gpidx = 0; gpidx < nlayers; gpidx++ ){
      sprintf( cgpidx, "_%02d_%03d_trkr%02d_", gidx, gpidx, TrackerID );
      // sprintf(gpname[gidx][gpidx][0], "gemplane_%02d_%03d_box", gidx, gpidx );
      // sprintf(gpname[gidx][gpidx][1], "gemplane_%02d_%03d_log", gidx, gpidx );
      // sprintf(gpname[gidx][gpidx][2], "gemplane_%02d_%03d_phy", gidx, gpidx );
      G4String gempboxname = TrackerPrefix + G4String("_gemplane") + cgpidx + G4String("box");
      G4String gemplogname = TrackerPrefix + G4String("_gemplane") + cgpidx + G4String("log");
      G4String gempphysname = TrackerPrefix + G4String("_gemplane") + cgpidx + G4String("phy");

      ztemp += gempz[gpidx]/2.0;

      gpbox = new G4Box( gempboxname, wplanes[gidx]/2.0, hplanes[gidx]/2.0, gempz[gpidx]/2.0 );
      gplog = new G4LogicalVolume( gpbox, gempm[gpidx], gemplogname, 0, 0, 0 );

      //gplog->SetVisAttributes( 

      new G4PVPlacement( 0, G4ThreeVector( 0.0, 0.0, ztemp - gempzsum/2.0 ), gplog, gempphysname, gemlog, false, 0, false ); //??

      ztemp += gempz[gpidx]/2.0;

      //Assign sensitive volume: why 5?  // SPR: This is the gas drift region
      if( gpidx == 5 ){
	gplog->SetSensitiveDetector(GEMSD);
	gplog->SetVisAttributes( gemsdvisatt );

	//		G4cout << "Assigning sensitive volume " << gempphysname << " to tracker ID " << TrackerID << G4endl;

	(GEMSD->GEMTrackerIDs)[ gempphysname ] = TrackerID;
	(GEMSD->GEMArmIDs)[TrackerID]          = (fDetCon->TrackerArm)[TrackerID];
      } else {
	gplog->SetVisAttributes( G4VisAttributes::Invisible );
      }
    }

    G4String gemname = TrackerPrefix + G4String("_gemphys_") + cgidx;
    //Now place the fully constructed GEM plane in the mother logical volume:
    //G4ThreeVector plane_pos = pos + G4ThreeVector( 0.0, 0.0, zplanes[gidx] );
    G4ThreeVector plane_pos = pos + zplanes[gidx] * zaxis;
    //Now we are positioning the GEM AFTER positioning all of its components inside its logical volume:
    new G4PVPlacement( rot, plane_pos, gemlog, gemname, Mother, true, gidx+1, false );
  }
}


