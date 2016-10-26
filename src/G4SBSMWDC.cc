#include "G4SBSMWDC.hh"

#include "G4SBSDetectorConstruction.hh"

#include "G4UserLimits.hh"
#include "G4SDManager.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4PVPlacement.hh"
#include "TMath.h"

#include "sbstypes.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "TVector3.h"

#include <cmath>

G4SBSMWDC::G4SBSMWDC(G4SBSDetectorConstruction *dc):G4SBSComponent(dc){
  fFieldD = 90.0*um;
  fSignalD = 20.0*um;
  fWireSep = 1.0*cm;
  fCath2WireDist = 3.0*mm;

  fCathodeThick = 12.0*um;
  fCuThick = 120.0*nm;
  fMylarThick = 12.0*um - 2.0*fCuThick;
  fPlaneThick = 2.0*fCathodeThick + 2.0*fCath2WireDist;

  fGasWinThick = 12.0*um;

  // **NOTE** -- abs val of angles are assumed to be equal!
  fUtheta = -30.0*deg;
  fVtheta =  30.0*deg;

  int chamberN[3] = {  1,   2,   3};
  int nplanes[3]  = {  6,   3,   6};
  int nwires[3]   = {140, 200, 200};

  double wirespace[3] = {1.0*cm, 1.0*cm, 1.0*cm};     // cm
  double height[3]    = {1.40*m, 2.00*m, 2.0*m };     // m
  double width[3]     = {0.35*m, 0.50*m, 0.50*m};     // m
  double z0_dis[3]    = {0.00*m, 0.3598*m, 0.705*m};     // m - numbers are taken from pg 89 Seamus Thesis

  G4String pattern[15] = { "U","U","X","X","V","V",
			   "U","X","V",
			   "U","U","X","X","V","V" };

  // *** Distanes between MWDC Planes need some special care ***
  // number calculated by comparing pg 89 of Seamus' Thesis to g4sbs, 
  // see spreadsheet for more information on MWDC positions:
  fSpacer = 0.000376*m; 

  // These are spacing differences between g4sbs and page 89 of Seamus' thesis where
  // a survey recorded z0 distances for each plane. I think these distances are important
  // as they are on the order of wire spacing => let's incorporate it into the simulation:
  // chamber 2, plane 2 has negative spacing so I just made it zero to avoid overlaps
  double offsets[15] = {0.00*m, 0.00*m, 0.00*m, 0.0032*m, 0.0032*m, 0.0032*m, 
			0.00*m, 0.00*m, 0.0099*m,
			0.00*m, 0.00*m, 0.00*m, 0.0074*m, 0.0079*m, 0.0074*m };
  
  double diff_survey_g4sbs[15] = {0.00*m, fSpacer, fSpacer, fSpacer, fSpacer, fSpacer, 
				  0.00*m, fSpacer, fSpacer,
				  0.00*m, fSpacer, 0.00*m, fSpacer, fSpacer, fSpacer};

  std::vector<double> temp0,temp1,temp2,tempp0,tempp1,tempp2;
  // Fill our GEn map:
  for(int i=0; i<nplanes[0]; i++){
    fChamber0.push_back( pattern[i] );
    temp0.push_back( diff_survey_g4sbs[i] );
    tempp0.push_back( offsets[i] );
  }
  for(int i=nplanes[0]; i<(nplanes[0]+nplanes[1]); i++){
    fChamber1.push_back( pattern[i] );
    temp1.push_back( diff_survey_g4sbs[i] );
    tempp1.push_back( offsets[i] );
  }
  for(int i=(nplanes[0]+nplanes[1]); i<(nplanes[0]+nplanes[1]+nplanes[2]); i++){
    fChamber2.push_back( pattern[i] );
    temp2.push_back( diff_survey_g4sbs[i] );
    tempp2.push_back( offsets[i] );
  }

  fGEn_Setup[0] = fChamber0;
  fGEn_Setup[1] = fChamber1;
  fGEn_Setup[2] = fChamber2;
  fSpace_Survey[0] = temp0;
  fSpace_Survey[1] = temp1;
  fSpace_Survey[2] = temp2;
  fSpace_Offset[0] = tempp0;
  fSpace_Offset[1] = tempp1;
  fSpace_Offset[2] = tempp2;

  // Fill all the vectors
  for(int i=0; i<3; i++){
    fChamberNumber.push_back( chamberN[i] );
    fNplanes.push_back( nplanes[i] );
    fNwires.push_back( nwires[i] );
    fNwirespacing.push_back( wirespace[i] );
    fNheight.push_back( height[i] );
    fNwidth.push_back( width[i] );
    fDist_z0.push_back( z0_dis[i] );
  }

  // Needed to build chamber first, then fill it with planes => Get depth of chamber needed to avoid seg fault:
  fSpace_Total_Chamber[0] = (0.000000 + 0.000000 + 0.000000 + 0.003200 + 0.003200 + 0.003200)*m + fSpacer*6.0;
  fSpace_Total_Chamber[1] = (0.000000 + 0.000000 + 0.009900)*m + fSpacer*3.0;
  fSpace_Total_Chamber[2] = (0.000000 + 0.000000 + 0.000000 + 0.007400 + 0.007900 + 0.007400)*m + fSpacer*6.0;

  mylarVisAtt = new G4VisAttributes(G4Colour( 0.5, 0.5, 0.5 ) );
  mylarVisAtt->SetForceWireframe(true);
  
  cuVisAtt = new G4VisAttributes(G4Colour(0.76,0.47,0.043));
  cuVisAtt->SetForceWireframe(true);

  winVisAtt = new G4VisAttributes(G4Colour(0.25,0.86,0.94));
  winVisAtt->SetForceWireframe(true);

  gasVisAtt = new G4VisAttributes(G4Colour(1.0,0.0,0.0,0.08));
  //gasVisAtt->SetForceWireframe(true);

  mothVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  mothVisAtt->SetForceWireframe(true);

  chamVisAtt = new G4VisAttributes(G4Colour(G4Colour::Magenta()));
  chamVisAtt->SetForceWireframe(true);

  sigwireVisAtt = new G4VisAttributes( G4Colour(0.0,0.8,0.0) );
  fieldwireVisAtt = new G4VisAttributes(G4Colour(G4Colour::Red()));
 
  // Define rotations for the wires
  fWireRotX = new G4RotationMatrix;
  fWireRotX->rotateY(90.0*deg);

  fWireRotU = new G4RotationMatrix;
  fWireRotU->rotateY(90.0*deg);
  fWireRotU->rotateX(fUtheta);

  fWireRotV = new G4RotationMatrix;
  fWireRotV->rotateY(90.0*deg);
  fWireRotV->rotateX(fVtheta);
}

G4SBSMWDC::~G4SBSMWDC(){;}

void G4SBSMWDC::BuildComponent(G4LogicalVolume *){
  ;
}

// Note: realworld was used for testing, nothing more. It was difficult to see if my building algorithm 
//       was doing what I thought it was doing, so I needed a way to test it.
void G4SBSMWDC::BuildComponent( G4LogicalVolume* realworld, G4LogicalVolume* world, G4RotationMatrix* rot, 
				G4ThreeVector pos, G4String SDname ) {
  double mX = 1.0*m;
  double mY = 2.1*m;
  //double mZ = 1.0*m;
  double mZ = 0.80*m;
  G4Box* mother = new G4Box("mother", mX/2.0, mY/2.0, mZ/2.0);
  G4LogicalVolume* mother_log = new G4LogicalVolume(mother,GetMaterial("Air"),"mother_log");
  mother_log->SetVisAttributes( G4VisAttributes::Invisible );
  mother_log->SetVisAttributes( chamVisAtt );
  G4ThreeVector origin(0.0,0.0,0.0);

  ///////////////////////////////////////////////////////////////////////////////
  // Define MWDC Sensitivity
  ///////////////////////////////////////////////////////////////////////////////

  G4String MWDCSDname = SDname;
  G4String MWDCSDname_nopath = SDname;
  MWDCSDname_nopath.remove(0,SDname.last('/')+1);
  G4String MWDCcolname = MWDCSDname_nopath;
  MWDCcolname += "HitsCollection";

  if( !(fMWDCSD = (G4SBSMWDCSD*) fDetCon->fSDman->FindSensitiveDetector(MWDCSDname)) ){
    fMWDCSD = new G4SBSMWDCSD( MWDCSDname, MWDCcolname );
    fDetCon->fSDman->AddNewDetector(fMWDCSD);
    (fDetCon->SDlist).insert(MWDCSDname);
    fDetCon->SDtype[MWDCSDname] = kMWDC;
  }


  ///////////////////////////////////////////////////////////////////////////////
  // Cathodes:
  ///////////////////////////////////////////////////////////////////////////////

  // Generate the Cathodes for chambers 0-2, and
  // put into a vector

  for(unsigned int i=0; i<fChamberNumber.size(); i++){
    G4Box* cathode_moth = new G4Box("temp", fNwidth[i]/2.0, fNheight[i]/2.0, fCathodeThick/2.0);
    G4Box* copper_box = new G4Box("copper_box", fNwidth[i]/2.0, fNheight[i]/2.0, fCuThick/2.0);
    G4Box* mylar_box = new G4Box("mylar_box", fNwidth[i]/2.0, fNheight[i]/2.0, fMylarThick/2.0);

    G4LogicalVolume* cathode_log = new G4LogicalVolume(cathode_moth, GetMaterial("Air"), "cathode_log");
    G4LogicalVolume* copper_log = new G4LogicalVolume(copper_box, GetMaterial("Copper"), "copper_log");
    G4LogicalVolume* mylar_log = new G4LogicalVolume(mylar_box, GetMaterial("Mylar"), "mylar_log");

    // Colors:
    cathode_log->SetVisAttributes( G4VisAttributes::Invisible );
    copper_log->SetVisAttributes( cuVisAtt );
    mylar_log->SetVisAttributes( mylarVisAtt );

    double z = fCathodeThick/2.0 - fCuThick/2.0;
    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, z), copper_log, "copper_phys_back", cathode_log, false, 0);
    z = -fCathodeThick/2.0 + fCuThick/2.0;
    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, z), copper_log, "copper_phys_front", cathode_log, false, 1);
    new G4PVPlacement(0, origin, mylar_log, "mylar_phys", cathode_log, false, 0);
    
    fCathodes[i] = cathode_log;
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////
  // CHAMBERS
  /////////////////////////////////////////////////////////////////////////////////////////////////////

  std::map<int,std::vector<G4String> >::iterator mit;
  std::vector<G4String>::iterator vit;
  char temp_name[255];
  int copyID = 0;
  G4LogicalVolume* test_log;

  for( mit = fGEn_Setup.begin(); mit != fGEn_Setup.end(); mit++ ) {
    // Make a Chamber:
    int chamber_number = mit->first;
    int num_planes = (mit->second).size();
    double chamber_thick = fPlaneThick*num_planes + 2.0*fGasWinThick + fSpace_Total_Chamber[chamber_number];

    // I added an "arbitrary" 2.0*cm to the chamber height, this is a result of offsetting
    // the signal / field wires.
    sprintf(temp_name, "chamber_%1d_box", chamber_number);
    G4Box *chamber_temp = new G4Box(temp_name, fNwidth[chamber_number]/2.0, 
    				    fNheight[chamber_number]/2.0 + 2.0*cm/2.0, chamber_thick/2.0);
  
    sprintf(temp_name, "chamber_%1d_log", chamber_number);
    G4LogicalVolume* chamber_log = new G4LogicalVolume(chamber_temp, GetMaterial("Air"), temp_name);
    chamber_log->SetVisAttributes( chamVisAtt );
    
    // Make the gas windows, and place it at front / back of a chamber:    
    G4Box* gas_winbox = new G4Box("gas_winbox", fNwidth[chamber_number]/2.0, 
				  fNheight[chamber_number]/2.0, fGasWinThick/2.0);
    G4LogicalVolume* gas_winlog = new G4LogicalVolume(gas_winbox, GetMaterial("Mylar"), 
						      "gas_winlog" );
    gas_winlog->SetVisAttributes( winVisAtt );

    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, -chamber_thick/2.0 + fGasWinThick/2.0),
    		      gas_winlog, "gas_phys_front", chamber_log, false, 0);
    new G4PVPlacement(0, G4ThreeVector(0.0, 0.0,  chamber_thick/2.0 - fGasWinThick/2.0),
    		      gas_winlog, "gas_phys_back", chamber_log, false, 0);
    
    int planeN = 0;

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // PLANES
    ///////////////////////////////////////////////////////////////////////////////////////////////

    std::vector<double> plane_positions_from_survey = fSpace_Survey[chamber_number];
    std::vector<double> plane_offset = fSpace_Offset[chamber_number];
    int weight = 1;
     
    for( vit = mit->second.begin(); vit != mit->second.end(); vit++ ) {

      // This annoying routine was the only way I could get the planes to match Seamus' Thesis. 
      // When I did a simple difference between G4SBS numbers and Seamus' numbers, that gave me the 
      // "offsets" described in the Constructor. However, in order to get the placement to match,
      // I needed some stupid averaging - hence the introduction of "weight"
      double plane_placement = 0.0*m;
      for( int pln=0; pln<=planeN; pln++ ){
	if(pln<3){
	  plane_placement += plane_positions_from_survey[pln] + plane_offset[pln];
	} else {
	  plane_placement += plane_positions_from_survey[pln] + plane_offset[pln] / double(weight);
	}
      }

      if( planeN >=3 ) weight++;

      G4String plane_type = *vit;
      G4LogicalVolume* plane_log;
   
      if( plane_type == "X" ) {
      	plane_log = BuildX( fNwidth[chamber_number], fNheight[chamber_number], 
			    chamber_number, planeN, copyID );
      } else {
	plane_log = BuildUorV( fNwidth[chamber_number], fNheight[chamber_number], plane_type,
			       chamber_number, planeN, copyID );
	// Tests:
	// if( chamber_number==1 ){
	//   if( planeN == 2 ) {
	//     test_log = plane_log;
	//   }
	// }
      }
 
      // Place the built plane inside a chamber:
      double z = -chamber_thick/2.0 + fGasWinThick + (planeN+0.5)*fPlaneThick + plane_placement;
      sprintf(temp_name, "chamber%1d_plane%1d_log", chamber_number, planeN);
      new G4PVPlacement(0, G4ThreeVector(0.0,0.0,z), plane_log, temp_name, chamber_log, false, copyID);

      planeN++;
      copyID++;
    }
    // Place the Chamber inside our mother:
    sprintf(temp_name, "chamber%1d_phys",chamber_number);    
    new G4PVPlacement(0,G4ThreeVector(0.0,0.0, -mZ/2.0 + chamber_thick/2.0 + fDist_z0[chamber_number]), 
    		      chamber_log, temp_name, mother_log, false, chamber_number);
  }
  // Place the MWDC (Chambers 0-2, all planes) inside the world:
  new G4PVPlacement(rot, pos, mother_log, "MWDC_mother_phys", world, false, 0);

  // Test to see what a single plane looks like:
  //new G4PVPlacement( rot, G4ThreeVector(0.0, 2.0*m, 2.0*m), test_log, "test", realworld, false, 0 ); 
}


// Arguments: width = plane width
//            height = plane height
//            chamber = chamber #
//            planeN = plane number in chamber
//            planeN is also used as a counting device to get signal/field wire offset. 
//            Successive planes of the same type have offset = 0.5*cm in order to increase
//            the resolution
//            copyID = global plane number - needed for DetMap

G4LogicalVolume* G4SBSMWDC::BuildX(double width, double height, int chamber, int planeN, int global_planeN) {
  char temp_name[255];
  sprintf(temp_name, "X_chamber%1d_plane%1d_box", chamber, planeN);
 
  G4Box* Xbox = new G4Box(temp_name, width/2.0, height/2.0, fPlaneThick/2.0);

  sprintf(temp_name, "X_chamber%1d_plane%1d_log", chamber, planeN);
  G4LogicalVolume* Xlog = new G4LogicalVolume( Xbox, GetMaterial("MWDC_gas"), temp_name );
  Xlog->SetVisAttributes( gasVisAtt );
  Xlog->SetSensitiveDetector( fMWDCSD );

  // Put in the cathodes
  sprintf(temp_name, "X_chamber%1d_plane%1d_cathodefront", chamber, planeN);
  new G4PVPlacement(0, G4ThreeVector(0.0,0.0,-fPlaneThick/2.0 + fCathodeThick/2.0), fCathodes[chamber], 
  		    temp_name, Xlog, false, 0);
  sprintf(temp_name, "X_chamber%1d_plane%1d_cathodeback", chamber, planeN);
  new G4PVPlacement(0, G4ThreeVector(0.0,0.0,fPlaneThick/2.0 - fCathodeThick/2.0), fCathodes[chamber], 
  		    temp_name, Xlog, false, 0);

  // Now Let's make the signal / field wires:
  G4Tubs *signaltub = new G4Tubs( "signal_tub", 0.0*cm, fSignalD/2.0, width/2.0, 0.0, twopi );
  G4Tubs *fieldtub = new G4Tubs( "field_tub", 0.0*cm, fFieldD/2.0, width/2.0, 0.0, twopi );

  // Signal wire is actually gold plated-tungsten, not sure on the dimensions so I simply chose Tungsten
  G4LogicalVolume* signal_log = new G4LogicalVolume( signaltub, GetMaterial("Tungsten"), "signal_log" );
  signal_log->SetVisAttributes( sigwireVisAtt );

  // Field wire is actually copper-beryllium, but beryllium (according to google) is < 3% so I ignored it
  G4LogicalVolume* field_log = new G4LogicalVolume( fieldtub, GetMaterial("Copper"), "field_log" );
  field_log->SetVisAttributes( fieldwireVisAtt );

  // Offset for two adjacent planes of same type
  double offset = 0.0*cm;
  if( (planeN+1) % 2 == 0 && fNplanes[chamber] > 3) {
    offset = fWireSep / 2.0; // should be 0.5*cm for GEn
  }

  ///////////////////////////////////////////////////////////////////////////////
  // Generate wire mesh:
  ///////////////////////////////////////////////////////////////////////////////

  // Used to give unique ID number to field wires:
  // signal wires are 0->fNWires[chamber]-1
  // field wires are fNwires[chamber]->2*fNwires[chamber]-1
  int field_count = fNwires[chamber];
  G4ThreeVector w0_pos(0.0,0.0,0.0);

  for( int i=0; i<fNwires[chamber]; i++ ){
    double y_signal = -height/2.0 + i*fWireSep + offset + fSignalD/2.0;
    double y_field  = y_signal + fWireSep/2.0;

    sprintf(temp_name, "X_chamber%1d_plane%1d_signalwire%1d", chamber, planeN, i);
    new G4PVPlacement( fWireRotX, G4ThreeVector(0.0, y_signal, 0.0), signal_log, temp_name,
		       Xlog, false, i );

    sprintf(temp_name, "X_chamber%1d_plane%1d_fieldwire%1d", chamber, planeN, field_count); 
    new G4PVPlacement( fWireRotX, G4ThreeVector(0.0, y_field, 0.0), field_log, temp_name,
		       Xlog, false, field_count );

    field_count++;

    // Grab the position of the first signal wire for DetMap:
    if( i == 0 ) {
      w0_pos = G4ThreeVector(0.0, y_signal, 0.0);
    }
  }

  // Fill the DetMap:
  (fMWDCSD->detmap).w0[global_planeN] = w0_pos;
  (fMWDCSD->detmap).WireSpacing[global_planeN] = fNwirespacing[chamber];
  // These coordinates will be rotated in EventAction, so it makes sense that Px is the only component
  (fMWDCSD->detmap).Px[global_planeN] = 1.0; 
  (fMWDCSD->detmap).Py[global_planeN] = 0.0; 
  (fMWDCSD->detmap).Plane_nhat[global_planeN] = G4ThreeVector(0.0, 1.0, 0.0);

 return Xlog;
}

G4LogicalVolume* G4SBSMWDC::BuildUorV(double width, double height, G4String type, int chamber, int planeN, int global_planeN) {
  char temp_name[255];
  char temp_namelog[255];

  if( type == "U" ) {
    sprintf(temp_name, "U_chamber%1d_plane%1d_box", chamber, planeN);
    sprintf(temp_namelog, "U_chamber%1d_plane%1d_log", chamber, planeN);
  } else {
    sprintf(temp_name, "V_chamber%1d_plane%1d_box", chamber, planeN);
    sprintf(temp_namelog, "V_chamber%1d_plane%1d_log", chamber, planeN);
  }
 
  G4Box* UorVbox = new G4Box(temp_name, width/2.0, height/2.0, fPlaneThick/2.0);

  G4LogicalVolume* UorVlog = new G4LogicalVolume( UorVbox, GetMaterial("MWDC_gas"), temp_namelog );
  UorVlog->SetVisAttributes( gasVisAtt );
  UorVlog->SetSensitiveDetector( fMWDCSD );

  // Put in the cathodes
  if( type == "U" ) {
    sprintf(temp_name, "U_chamber%1d_plane%1d_cathodefront", chamber, planeN);
    sprintf(temp_name, "U_chamber%1d_plane%1d_cathodeback", chamber, planeN);
  } else {
    sprintf(temp_name, "V_chamber%1d_plane%1d_cathodefront", chamber, planeN);
    sprintf(temp_name, "V_chamber%1d_plane%1d_cathodeback", chamber, planeN);
  }

  new G4PVPlacement(0, G4ThreeVector(0.0,0.0,-fPlaneThick/2.0 + fCathodeThick/2.0), fCathodes[chamber], 
  		    temp_name, UorVlog, false, 0 );

  new G4PVPlacement(0, G4ThreeVector(0.0,0.0,fPlaneThick/2.0 - fCathodeThick/2.0), fCathodes[chamber], 
  		    temp_name, UorVlog, false, 0 );

  double offset = 0.0*cm;
  if( (planeN+1) % 2 == 0 && fNplanes[chamber] > 3) {
    offset = fWireSep / 2.0;
  }

  // Vertical offsets should be the same (assuming the angles are the same),
  // so its easier to work with a positive angle:
  double utheta = fabs(fUtheta);
  double vtheta = fabs(fVtheta);
  TVector3 nhat(0.0,0.0,0.0);
  TVector3 nhat_top(0.0,0.0,0.0);

  if( type == "U" ) {
    TVector3 uhat( cos(utheta), sin(utheta), 0.0 );
    TVector3 uhat_top( -cos(utheta), -sin(utheta), 0.0 );
    nhat = uhat;
    nhat_top = uhat_top;
  } else {
    TVector3 vhat( -cos(vtheta), sin(vtheta), 0.0 );
    TVector3 vhat_top( cos(vtheta), -sin(vtheta), 0.0 );
    nhat = vhat;
    nhat_top = vhat_top;
  }
  double initial_offset = 0.25*cm;

  double length = width / TMath::Cos( utheta );
  double wiresep = fWireSep / TMath::Cos( utheta );

  int field_count = fNwires[chamber];

  double y_start = -0.5*height - 0.5*length*sin( utheta );      // start iteration from here  
  double y_normal = -0.5*height + 0.5*length*sin( utheta );     // where lengths are normal
  double y_lowest = y_start;
  double y_highest = 0.5*height + 0.5*length*sin( utheta ); 
  double y_normal_top = 0.5*height - 0.5*length*sin( utheta );

  G4ThreeVector w0_pos(0.0,0.0,0.0);

  for( int i=0; i<fNwires[chamber]-2; i++ ){
      double y_signal = y_start + i*wiresep + offset + fSignalD/2.0 + initial_offset;   // increment through the detector
      double y_field  = y_signal + wiresep/2.0;

      // need to find the y difference b/t y_lowest and y_signal:
      double diff = 0.0;
      double diff_field = 0.0;

      if( y_signal < y_normal ) {
	diff = fabs(y_lowest - y_signal);
      }
      if( y_signal > y_normal_top ) {
	diff = fabs(y_highest - y_signal);
      }
      
      if( y_field < y_normal ) {
	diff_field = fabs(y_lowest - y_field);
      }
      if( y_field > y_normal_top ) {
	diff_field = fabs(y_highest - y_field);
      }
      
      double newlength = 0.0;
      double newlength_field = 0.0;
      if( type == "U" ) {
	newlength = diff / sin( utheta ); 
	newlength_field = diff_field / sin( utheta );
      } else {
	newlength = diff / sin( vtheta ); 
	newlength_field = diff_field / sin( vtheta );
      }

      if( fabs(y_signal) <= fabs(y_normal) ) {
	newlength = length;
      }
      if( fabs(y_field) <= fabs(y_normal) ){
	newlength_field = length;
      }

      G4Tubs *signaltub = new G4Tubs( "wire_tub", 0.0*cm, fSignalD/2.0, newlength/2.0, 0.0, twopi );
      G4LogicalVolume* signal_log = new G4LogicalVolume( signaltub, GetMaterial("Tungsten"), "signal_log" );
      signal_log->SetVisAttributes( sigwireVisAtt );

      G4Tubs *fieldtub = new G4Tubs( "field_tub", 0.0*cm, fFieldD/2.0, newlength_field/2.0, 0.0, twopi );
      G4LogicalVolume* field_log = new G4LogicalVolume( fieldtub, GetMaterial("Copper"), "field_log" );
      field_log->SetVisAttributes( fieldwireVisAtt );

      // Place the Wires:
      TVector3 zero(0.0,0.0,0.0);
      double offset_from_ysig = length/2.0 - newlength;
      TVector3 start(0.0, y_signal, 0.0);
      TVector3 move = zero;
      TVector3 wireC = zero;

      double offset_from_ysig_field = length/2.0 - newlength_field;
      TVector3 start_field(0.0, y_field, 0.0);
      TVector3 move_field = zero;
      TVector3 wireC_field = zero;

      if( y_signal < y_normal ) {
	move = start + offset_from_ysig * nhat;
	wireC = move + newlength/2.0 * nhat;
      }
      if( y_signal > y_normal_top ) {
	move = start + offset_from_ysig * nhat_top;
	wireC = move + newlength/2.0 * nhat_top;
      }
      if( fabs(y_signal) <= fabs(y_normal) ) {
	wireC = start;
      }

      if( y_field < y_normal ) {
	move_field = start_field + offset_from_ysig_field * nhat;
	wireC_field = move_field + newlength_field/2.0 * nhat;
      }
      if( y_field > y_normal_top ) {
	move_field = start_field + offset_from_ysig_field * nhat_top;
	wireC_field = move_field + newlength_field/2.0 * nhat_top;
      }
      if( fabs(y_field) <= fabs(y_normal) ){
	wireC_field = start_field;
      }

      if( type == "U" ) {
	sprintf(temp_name, "U_chamber%1d_plane%1d_signalwire%1d", chamber, planeN, i);
	new G4PVPlacement( fWireRotU, G4ThreeVector(wireC.X(), wireC.Y(), wireC.Z()), signal_log, 
			   temp_name, UorVlog, false, i );

	sprintf(temp_name, "U_chamber%1d_plane%1d_fieldwire%1d", chamber, planeN, field_count);
	new G4PVPlacement( fWireRotU, G4ThreeVector(wireC_field.X(), wireC_field.Y(), wireC_field.Z()), field_log, 
			   temp_name, UorVlog, false, field_count );

      } else {
	sprintf(temp_name, "V_chamber%1d_plane%1d_signalwire%1d", chamber, planeN, i);
	new G4PVPlacement( fWireRotV, G4ThreeVector(wireC.X(), wireC.Y(), wireC.Z()), signal_log, 
			   temp_name, UorVlog, false, i );

	sprintf(temp_name, "V_chamber%1d_plane%1d_fieldwire%1d", chamber, planeN, field_count);
	new G4PVPlacement( fWireRotV, G4ThreeVector(wireC_field.X(), wireC_field.Y(), wireC_field.Z()), field_log, 
			   temp_name, UorVlog, false, field_count );
      }

      if( i == 0 ){
	w0_pos = G4ThreeVector(wireC.X(), wireC.Y(), wireC.Z());
      }
     
      field_count++;
  }

  // Fill the DetMap:
  (fMWDCSD->detmap).w0[global_planeN] = w0_pos;
  (fMWDCSD->detmap).WireSpacing[global_planeN] = fNwirespacing[chamber];
  if( type == "U" ) {
    (fMWDCSD->detmap).Px[global_planeN] = cos( fabs(fUtheta) );
    (fMWDCSD->detmap).Py[global_planeN] = sin( fabs(fUtheta) );
    (fMWDCSD->detmap).Plane_nhat[global_planeN] = 
      G4ThreeVector( -cos(fabs(fUtheta)), sin(fabs(fUtheta)), 0.0 );
  }
  if( type == "V" ) {
    (fMWDCSD->detmap).Px[global_planeN] = cos( fabs(fVtheta) );
    (fMWDCSD->detmap).Py[global_planeN] = -sin( fabs(fVtheta) );
    (fMWDCSD->detmap).Plane_nhat[global_planeN]
      = G4ThreeVector( cos(fabs(fUtheta)), sin(fabs(fUtheta)), 0.0 );
  }     
  return UorVlog;
}
