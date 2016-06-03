#include "G4SBSNeutronDetector.hh"
#include "G4SBSDetectorConstruction.hh"

#include "G4SBSECalSD.hh"
#include "G4SBSCalSD.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SDManager.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "TString.h"
#include <map>
#include <iostream>
#include <vector>

#define LAYERTOSHOW 2


G4SBSNeutronDetector::G4SBSNeutronDetector(G4SBSDetectorConstruction *dc) : G4SBSComponent(dc) {
  int i, j;
  G4double inch = 2.54*cm;

  fNDdist = fDetCon->GetNDdist() * m;
  fNDang  = fDetCon->GetNDang() * deg;
  fTarget = fDetCon->GetGEnTarget();

  // Scintillator block dimensions
  //  Veto long length is 100 cm (despite what Tim's document says)
  //                         CMU       UVA    JLab       Glas      Vlong     Vshort
  double x[NBARTYPES]    = {  15.0*cm,  10.0*cm,  10.0*cm,  20.0*cm,  11.0*cm,  11.0*cm};
  double y[NBARTYPES]    = { 180.0*cm, 160.0*cm, 160.0*cm, 190.0*cm, 100.0*cm,  70.0*cm};
  double z[NBARTYPES]    = {   5.0*cm,  10.0*cm,  10.0*cm,  10.0*cm,   2.0*cm,   2.0*cm};
  double atln[NBARTYPES] = {   1.0*m ,   1.0*m ,   1.0*m ,   1.0*m ,   1.0*m ,   1.0*m };

  // Lightguides
  double lgx[NBARTYPES]={      2.0*inch,   2.0*inch,   2.0*inch,   2.0*inch,   1.1*inch,   1.1*inch};
  double lgfdep[NBARTYPES]={   5.0*cm,      10.0*cm,   10.0*cm,    10.0*cm,   2.0*cm,     2.0*cm};
  double lgndep[NBARTYPES]={   2.0*inch,   2.0*inch,   2.0*inch,   2.0*inch,   2.0*cm,     2.0*cm};
  double lgang[NBARTYPES] ={  16.0*deg,    20.0*deg,  20.0*deg,   30.0*deg,    16.0*deg,  16.0*deg};
  double lglen[NBARTYPES];

  for( i = 0; i < NBARTYPES; i++ ){
    lglen[i] = (x[i]-lgx[i])/(2.0*tan(lgang[i]));
  }

  // Cassette dimensions
  //                                   CMU         UVA      JLab         Glas         V1         V2
  double casx[NCASSETTE_TYPES]  = {  78.26*cm, 104.46*cm, 104.46*cm, 101.90*cm,  88.55*cm,  88.55*cm};
  double casy[NCASSETTE_TYPES]  = { 365.76*cm, 365.75*cm, 365.75*cm, 365.75*cm, 365.75*cm, 365.75*cm};
  double casz[NCASSETTE_TYPES]  = {  10.00*cm,  14.00*cm,  14.00*cm, 15.90*cm,   7.80*cm,   7.80*cm};

  int nbars[NCASSETTE_TYPES] = { 5, 10, 10, 5, 8, 8 };

  // Fe1 is really styrofoam for vetos
  double casfe1th[NCASSETTE_TYPES] = {   1.27*cm,   1.27*cm, 1.27*cm,  1.27*cm,   1.27*cm,   1.27*cm};
  //	double casfe1ln[NCASSETTE_TYPES]={ 182.88*cm, 162.56*cm,162.56*cm,193.04*cm,365.75*cm, 365.75*cm};
  double casfe1ln[NCASSETTE_TYPES] = { 182.88*cm, 162.56*cm,162.56*cm,193.04*cm,182.88*cm, 182.88*cm};

  // This does not exist for vetos
  double casfe2th[NCASSETTE_TYPES]={   1.27*cm,   1.27*cm, 1.27*cm,  1.27*cm,   1.27*cm,   1.27*cm};
  double casfe2ln[NCASSETTE_TYPES]={ 170.82*cm, 151.84*cm,151.84*cm,180.31*cm,  170.82*cm, 170.82*cm};
  double casfe2ht[NCASSETTE_TYPES]={  48.26*cm,  64.42*cm, 64.42*cm, 62.84*cm,  60.00*cm,  60.00*cm};

  // 	This is the back plate
  double casalth[NCASSETTE_TYPES]={   0.64*cm,   0.64*cm, 0.64*cm,   0.64*cm,   1.27*cm,   1.27*cm}; 

  // Everything relative to 164.5cm from Tim's document
  //                         CMU       UVA    JLab       Glas      V1     V2
  double relbary[NCASSETTE_TYPES] = { 0.0*cm, 0.0*cm, 0.0*cm, 17.4*cm, -2.5*cm, 2.6*cm };
  double bary[NCASSETTE_TYPES];

  // position of ND center relative to 164.5 measurement
  double yoffset = -20.0*cm;
  for( i = 0; i < NCASSETTE_TYPES; i++ ){
    bary[i] = relbary[i] + yoffset;
  }


  // Put 1mm in front of al plate for normal bars
  //
  // For veto 1, put in center of al plate flush with back and plate
  // 1mm behind front iron plate, backwards for veto 2
  double casbarz[NCASSETTE_TYPES];
  for( i = 0; i < NCASSETTE_TYPES; i++ ){
    switch(i){
    case kCassVeto1:
      casbarz[i] = ((casz[i]/2.0 - casalth[i]) 
		    + (-casz[i]/2.0 + casfe1th[i] + casalth[i] + 1.0*mm) )/2.0;
      break;
    case kCassVeto2:
      casbarz[i] = ((-casz[i]/2.0 + casalth[i]) 
		    + (casz[i]/2.0 - casfe1th[i] - casalth[i]) - 1.0*mm )/2.0;
      break;
    default:
      casbarz[i] = casz[i]/2.0 - casalth[i] - 10.0*mm - z[i]/2.0;
      break;
    }
  }

  // 0 is between 23rd and 24th bar on v2
  double offset = 259.37*cm;
  // V1    V2      N1         N2      N3       N4       N5     N6     N7
  double relbot[NLAYERS]  = {5.1*cm, 0.0*cm, 17.0*cm, 17.1*cm, 17.1*cm, 16.8*cm, 0.0*cm, -0.64*cm, 0.0*cm};
  double laybot[NLAYERS];
  for( i = 0; i < NLAYERS; i++ ){
    laybot[i] = offset-relbot[i];
  }
  double pmtrad[NLAYERS]  = {0.55*inch, 0.55*inch, 1.0*inch, 1.0*inch, 1.0*inch, 1.0*inch, 1.0*inch, 1.0*inch, 1.0*inch};

  // Relative to front of iron plate between vetos and N1
  // V1      V2       N1     N2     N3       N4     N5       N6     N7
  double layz[NLAYERS]={-21.7*cm, -13.5*cm, 8.3*cm,19.7*cm,31.1*cm,42.5*cm,53.9*cm,71.1*cm,88.3*cm};


  int       ncass[NLAYERS] = {6, 6, 6, 5, 6, 5, 5, 5, 5 };
  CassType_t castypes[NLAYERS][MAX_CASS] = 
    {       { kCassVeto1,kCassVeto1,kCassVeto1,kCassVeto1,kCassVeto1,kCassVeto1 },    // V1
	    { kCassVeto2,kCassVeto2,kCassVeto2,kCassVeto2,kCassVeto2,kCassVeto2 },    // V2
	    { kCassCMU, kCassCMU, kCassCMU, kCassCMU, kCassCMU, kCassGlasgow }, // N1
	    { kCassCMU, kCassCMU, kCassCMU, kCassCMU, kCassCMU, kCassNull},             // N2
	    { kCassCMU, kCassCMU, kCassCMU, kCassCMU, kCassCMU, kCassGlasgow }, // N3
	    { kCassCMU, kCassCMU, kCassCMU, kCassCMU, kCassCMU,    kCassNull},          // N4
	    { kCassUVA, kCassUVA, kCassUVA, kCassUVA, kCassGlasgow,kCassNull},          // N5
	    { kCassUVA, kCassUVA, kCassUVA, kCassUVA, kCassGlasgow,kCassNull},          // N6
	    { kCassUVA, kCassUVA, kCassUVA, kCassUVA, kCassGlasgow,kCassNull}           // N7
    };

  double platez[NPLATES]        = {-33.62*cm, -28.57*cm, 0.0*cm, 1.28*cm, 106.1*cm, 107.38*cm };
  double platethk[NPLATES]      = {5.04*cm, 1.27*cm, 1.27*cm, 2.54*cm, 1.27*cm, 2.54*cm};
  plateWidth  = 370.0*cm;
  plateHeight = 583.9*cm;


  // Set Values
  for( i = 0; i < NBARTYPES; i++ ){
    attlen[i] = atln[i];
    X[i]  = x[i];
    Y[i]  = y[i];
    Z[i]  = z[i];

    lgSize[i]= lgx[i];
    lgFarDepth[i]   = lgfdep[i];
    lgNearDepth[i]   = lgndep[i];
    lgLength[i]   = lglen[i];
    lgcylLen[i] = 3.0*inch;
  }

  for( i = 0; i < NCASSETTE_TYPES; i++ ){
    Nbars[i]= nbars[i];
    casX[i] = casx[i];
    casY[i] = casy[i];
    casZ[i] = casz[i];
    casbarY[i] = bary[i];
    barspacing[i] = casX[i]/Nbars[i];

    casbarZ[i] = casbarz[i];

    casAlthick[i] = casalth[i];
    casFe1length[i] = casfe1ln[i];
    casFe1thick[i] = casfe1th[i];

    casFe2length[i] = casfe2ln[i];
    casFe2thick[i] = casfe2th[i];
    casFe2height[i] = casfe2ht[i];
  }

  for( i = 0; i < NLAYERS; i++ ){
    Ncass[i] = ncass[i];
    layerBottom[i] = laybot[i];
    layerZ[i] = layz[i];

    for( j = 0; j < MAX_CASS; j++ ){
      casType[i][j] = castypes[i][j];
    }

    pmtRad[i] = pmtrad[i];
  }

  for( i = 0; i < NPLATES; i++ ){
    plateZ[i] = platez[i];
    plateThick[i] = platethk[i];
  }

  for( i = 0; i < NBARTYPES; i++ ){
    fThreshold[i]  = 10000.0;
  }
}

void G4SBSNeutronDetector::BuildComponent(G4LogicalVolume* world) {
  ConstructND( world );
}

void G4SBSNeutronDetector::ConstructND( G4LogicalVolume* world) {
  G4double inch = 2.54*cm;
  
  char barname[NBARTYPES][255] = { "CMU", "UVA", "JLAB", "GLAS", "VETOL", "VETOS" };
  char casname[NCASSETTE_TYPES][255] = { "CMU", "UVA", "JLAB", "GLAS", "VETO1", "VETO2" };

  G4Material *platemat[NPLATES] = {GetMaterial("Lead"), GetMaterial("Iron"), 
				   GetMaterial("Iron"), GetMaterial("Iron"), 
				   GetMaterial("Iron"), GetMaterial("Iron")};

  // VISUALS:
  G4VisAttributes* leadVisAtt = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  leadVisAtt->SetForceWireframe(true);
  G4VisAttributes* alVisAtt = new G4VisAttributes(G4Colour(0.5,0.5,0.5));
  alVisAtt->SetForceWireframe(true);

  G4VisAttributes* lgVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,1.0));
  G4VisAttributes* pmtVisAtt = new G4VisAttributes(G4Colour(0.1,0.1,0.1));
  G4VisAttributes* ironVisAtt = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  ironVisAtt->SetForceWireframe(true);
  G4VisAttributes* styroVisAtt = new G4VisAttributes(G4Colour(0.8,0.8,0.8));

  G4VisAttributes* invisVisAtt = new G4VisAttributes(G4Colour(0.0,0.0,0.0));
  invisVisAtt->SetVisibility(false);
  G4VisAttributes *ND_MotherVis = new G4VisAttributes(G4Colour::Blue());
  ND_MotherVis ->SetForceWireframe(true);
  G4VisAttributes *vetoVis  = new G4VisAttributes(G4Colour::Cyan());
  G4VisAttributes *ndVis = new G4VisAttributes(G4Colour(0.49,0.0,1.0));

  G4VisAttributes* barVisAtt     = new G4VisAttributes(G4Colour(1.0,1.0,0.0));
  G4VisAttributes* vetoVisAtt    = new G4VisAttributes(G4Colour(1.0,0.0,0.0));
  G4VisAttributes* veto2VisAtt   = new G4VisAttributes(G4Colour(1.0,0.4,0.2));
  G4VisAttributes* GlasVisAtt = new G4VisAttributes( G4Colour(0.78,0.78,0.0) );
  G4VisAttributes* CMUVisAtt = new G4VisAttributes( G4Colour(0.78,0.0,0.0) );
  G4VisAttributes* UVAVisAtt = new G4VisAttributes( G4Colour(0.0,0.784,0.0) );

  
  for( int i = 0; i < NPLATES; i++ ){
    plateMat.push_back( platemat[i] );
  }

  // BUILD NEUTRON DETECTOR
  NDangle = fNDang;
  NDdistance = fNDdist;

  // ************************************************************
  // this needs to be updated, need to find this information    *
  const G4double theta_norm = 30.17*deg;      
                
  //const G4double theta_norm = fNDang;                           
  const G4double bob_z0_dist = 64.0*cm;                         
  const G4double bob_targ_dist = NDdistance - bob_z0_dist;   
  //const G4double bob_targ_dist = NDdistance;
  //*************************************************************

  G4double ND_x = 3.0*m;
  G4double ND_y = 3.0*m;
  G4double ND_z = 2.0*m;

  G4RotationMatrix* ndRot = new G4RotationMatrix();
  ndRot->rotateY( -theta_norm );
  ndRot->rotateZ( 90.0*deg );
  G4Box* neutronarm_box = new G4Box("neutronarm_box",ND_x,ND_y,ND_z);
  G4LogicalVolume *neutronarm_log = new G4LogicalVolume(neutronarm_box,
							GetMaterial("Air"),"neutronarm_log");
  
  neutronarm_log->SetVisAttributes( ND_MotherVis );
  // Place Neutron Detector Mother Volume in the World:
  double x3 = bob_targ_dist*sin(NDangle)+bob_z0_dist*sin(theta_norm);
  double y3 = 0.0;
  double z3 = bob_targ_dist*cos(NDangle)+bob_z0_dist*cos(theta_norm);
  new G4PVPlacement(ndRot, G4ThreeVector(x3,y3,z3), neutronarm_log, "neutronarm", world, false, 0);


  ///////////////////////////////////////////////////////////////////////////////////////////////////
  // Make all the details.....                                                                     //
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  
  G4LogicalVolume *lightguidelog, *lightguidecyllog;
  G4Box  *solidBlock, *solidBar;
  G4Trd  *solidLG;
  G4Tubs *solidLGcyl, *solidPMTglass, *solidPMTcath, *solidPMT;
  G4LogicalVolume *logicPMTglass, *logicPMTcath;

  double lgoffset, barl;
  char blockName[255];

  G4RotationMatrix *lgleftrot = new G4RotationMatrix();
  lgleftrot->rotateY(90.0*deg );
  lgleftrot->rotateX(90.0*deg );
  G4RotationMatrix *lgrightrot = new G4RotationMatrix();
  lgrightrot->rotateY(90.0*deg );
  lgrightrot->rotateX(-90.0*deg );
  int layer;
  int type;

  double cathlen, glasslen;
  cathlen = 6.0*inch;
  glasslen= 1.0*inch;

  // Make two types of PMT Cathodes which will be a SD
  G4SDManager *sdman = fDetCon->fSDman;

  G4Tubs* veto_tub = new G4Tubs("veto_tub", 0.0, 0.55*inch*1.3, cathlen/2.0, 0.0, twopi );
  G4Tubs* others_tub = new G4Tubs("others_tub", 0.0, 1.0*inch*1.3, cathlen/2.0, 0.0, twopi );
  G4LogicalVolume* veto_log = new G4LogicalVolume(veto_tub, GetMaterial("Glass_HC"), "veto_log");  //this was Aluminum
  G4LogicalVolume* others_log = new G4LogicalVolume(others_tub, GetMaterial("Glass_HC"), "others_log"); //this was Aluminum
  
  // Set up Sensitive Detector:
  G4String NDSDname = "Harm/ND";
  G4String NDcollname = "NDHitsCollection";
  G4SBSECalSD* NDSD = NULL;

  // The depths are irrelevant for ND:

  if( !((G4SBSECalSD*) sdman->FindSensitiveDetector(NDSDname)) ) {
    G4cout << "Adding ND PMT Sensitive Detector to SDman..." << G4endl;
    NDSD = new G4SBSECalSD( NDSDname, NDcollname );
    sdman->AddNewDetector(NDSD);
    (fDetCon->SDlist).insert(NDSDname);
    fDetCon->SDtype[NDSDname] = kECAL;
    (NDSD->detmap).depth = 3; 
  }
  veto_log->SetSensitiveDetector(NDSD);
  others_log->SetSensitiveDetector(NDSD);

  // Need to make a bunch of scintillating bars, then assign them 
  // Cal sensitivity:
  G4String NDScintSDname = "Harm/NDScint";
  G4String NDScintcollname = "NDScintHitsCollection";
  G4SBSCalSD *NDScintSD = NULL;

  if( !((G4SBSCalSD*) sdman->FindSensitiveDetector(NDScintSDname)) ){
    G4cout << "Adding ND Scintillator Sensitive Detector to SDman..." << G4endl;
    NDScintSD = new G4SBSCalSD( NDScintSDname, NDScintcollname );
    sdman->AddNewDetector(NDScintSD);
    (fDetCon->SDlist).insert(NDScintSDname);
    fDetCon->SDtype[NDScintSDname] = kCAL;
    (NDScintSD->detmap).depth = 3;
  }

  // Necessary to make the scintillators a SD using Seamus' setup:
  map<int,G4LogicalVolume*> scint_map;
  G4Box* temp;
  G4LogicalVolume* temp_log;
  G4double mylar_thick = 0.001*2.54*cm; 

  for( type = kBarCMU; type < NBARTYPES; type++ ){
    
    // Make a Mother that will house mylar wrapping + scintillator:
    temp = new G4Box("temp", X[type]/2.0, Y[type]/2.0, Z[type]/2.0 );
    temp_log = new G4LogicalVolume(temp, GetMaterial("Air"), blockName); 

    // Make Scintillator, assign Sensitivity:
    sprintf(blockName, "%s_scint_temp", barname[type] );
    G4Box *scint_temp = new G4Box(blockName, X[type]/2.0 - mylar_thick,
				  Y[type]/2.0, Z[type]/2.0 - mylar_thick);

    sprintf(blockName, "%s_scint_log_temp", barname[type] );
    G4LogicalVolume *scint_log_temp = new G4LogicalVolume(scint_temp, GetMaterial("MWDC_Scint"), blockName); //changed material from Scintillator to MWDC_Scint
    scint_log_temp->SetSensitiveDetector(NDScintSD);

    if( (fDetCon->StepLimiterList).find( NDScintSDname ) != (fDetCon->StepLimiterList).end() ){
      temp_log->SetUserLimits(  new G4UserLimits(0.0, 0.0, 0.0, DBL_MAX, DBL_MAX) );
    }

    // Need to subtract scint_temp from temp in order to get a wrapping,
    // make it bigger in y-hat to avoid "undefined" subtraction behaviour:
    G4Box *subbox = new G4Box("subtemp", X[type]/2.0 - mylar_thick,
			      Y[type]/2.0 + 1.0*cm, Z[type]/2.0 - mylar_thick);

    sprintf(blockName, "%s_mylar_temp", barname[type] );
    G4SubtractionSolid *mylar_wrap = new G4SubtractionSolid( blockName, temp, subbox, 0, G4ThreeVector(0.0,0.0,0.0) );
 
    sprintf(blockName, "%s_mylar_log", barname[type] );
    G4LogicalVolume *mylar = new G4LogicalVolume(mylar_wrap, GetMaterial("Mylar"), blockName);

    sprintf(blockName, "%s_mylar_skin", barname[type] );
    new G4LogicalSkinSurface( blockName, mylar, GetOpticalSurface("Mirrsurf") );

    mylar->SetVisAttributes( G4VisAttributes::Invisible );
    temp_log->SetVisAttributes( G4VisAttributes::Invisible );
  
    // Place in the mother
    sprintf(blockName, "%s_mylar_phys", barname[type] );
    new G4PVPlacement( 0, G4ThreeVector(0.0,0.0,0.0), mylar, blockName, temp_log, false, 0 );
    sprintf(blockName, "%s_scint_phys", barname[type] );
    new G4PVPlacement( 0, G4ThreeVector(0.0,0.0,0.0), scint_log_temp, blockName, temp_log, false, 0 );

    switch(type){
    case kBarVetoShort:
      scint_log_temp->SetVisAttributes( vetoVisAtt  );
      break;
    case kBarVetoLong:
      scint_log_temp->SetVisAttributes( veto2VisAtt  );
      break;
    case kBarUVA:
      scint_log_temp->SetVisAttributes( UVAVisAtt  );
      break;
    case kBarCMU:
      scint_log_temp->SetVisAttributes( CMUVisAtt  );
      break;
    case kBarGlasgow:
      scint_log_temp->SetVisAttributes( GlasVisAtt  );
      break;
    default:
      scint_log_temp->SetVisAttributes( barVisAtt  );
    break;
    }
    scint_map[type] = temp_log;
  }
  
  // Get Scint Position within each bar type
  // Tim's document states that there are 340 bars. This code generates
  // 341 - the discrepancy is ***i think*** Plane 3, top row.

  // key is type of bar, value is the position. pmt needs two values for each bar, 
  // with the exception of vetos. They only have one PMT
  map<int,G4ThreeVector> scint_type_pos_logBar;
  map<int,vector<G4ThreeVector> > pmt_type_pos_logBar;

  // Seamus makes a mother volume called solidPMT, then inserts glass + cathode.
  for( layer = 0; layer < NLAYERS; layer++ ){
    
    sprintf(blockName, "solpmt_%1d",layer+1 );
    solidPMT= new G4Tubs(blockName, 0.0, pmtRad[layer]*1.3, (cathlen+glasslen)/2.0, 0.0, twopi  );

    sprintf(blockName, "solpmtgl_%1d",layer+1 );
    solidPMTglass= new G4Tubs(blockName, 0.0, pmtRad[layer], glasslen/2.0, 0.0, twopi  );
    sprintf(blockName, "logpmtgl_%1d",layer+1 );
    
    logicPMTglass = new G4LogicalVolume( solidPMTglass, GetMaterial("Glass_HC"), blockName ); // changed glass material
    logicPMTglass->SetVisAttributes( lgVisAtt );

    sprintf(blockName, "solpmtca_%1d",layer+1 );
    solidPMTcath = new G4Tubs(blockName, 0.0, pmtRad[layer]*1.3, cathlen/2.0, 0.0, twopi  );
    sprintf(blockName, "logpmtca_%1d",layer+1 );
    //logicPMTcath = new G4LogicalVolume( solidPMTcath, GetMaterial("Aluminum"), blockName );

    // Necessary as vetos are different from the "others" :
    if( layer==0 || layer==1 ) {
      logicPMTcath = veto_log;
      logicPMTcath->SetVisAttributes( vetoVis );
    } else {
      logicPMTcath = others_log;
      logicPMTcath->SetVisAttributes( ndVis );
    }

    // Here is an array of mother volumes which will be used during the cassette routine:
    sprintf(blockName, "logpmt_%1d",layer+1 );
    logicPMT[layer] = new G4LogicalVolume( solidPMT, GetMaterial("Air"), blockName );
    logicPMT[layer]->SetVisAttributes( invisVisAtt );

    ///////////////////////////////////////////////////////////
    // Place stuff in "PMT Mother", then put in an array 
    // called logicPMT where placement = layer.
    ///////////////////////////////////////////////////////////
    sprintf(blockName, "physpmtgl_%1d",layer+1 );
    new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, -cathlen/2.0), logicPMTglass, blockName, logicPMT[layer], false, 0 );
    sprintf(blockName, "physpmtca_%1d",layer+1 );
    new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, glasslen/2.0), logicPMTcath, blockName, logicPMT[layer], false, 0 );


    // Make the Scintillating Bars:
    for( type = kBarCMU; type < NBARTYPES; type++ ){

      sprintf(blockName, "%ssolblock_%1d", barname[type], layer+1 );
      solidBlock = new G4Box(blockName, X[type]/2.0, Y[type]/2.0, Z[type]/2.0 );

      sprintf(blockName, "%slogblock_%1d", barname[type], layer+1 );
  
      //  logBlock[type][layer] = new G4LogicalVolume(solidBlock, GetMaterial("Scintillator"), blockName);
      // make these before all the iterating, so we can assign a sensitivity.
      logBlock[type][layer] = scint_map[type]; 

      // Lightguide
      sprintf(blockName, "%ssollg_%1d", barname[type], layer+1 );
      solidLG = new G4Trd( blockName,  lgFarDepth[type]/2.0, lgNearDepth[type]/2.0, X[type]/2.0, lgSize[type]/2.0, lgLength[type]/2.0 );
      sprintf(blockName, "%sloglg_%1d", barname[type], layer+1 );
      lightguidelog = new G4LogicalVolume(solidLG, GetMaterial("pglass"), blockName);
      lightguidelog->SetVisAttributes(  lgVisAtt );

      sprintf(blockName, "%ssollgcyl_%1d", barname[type], layer+1 );
      solidLGcyl = new G4Tubs(blockName, 0.0, lgSize[type]/2.0, lgcylLen[type]/2.0, 0.0, twopi );
      sprintf(blockName, "%sloglgcyl_%1d", barname[type], layer+1 );
      lightguidecyllog = new G4LogicalVolume(solidLGcyl, GetMaterial("pglass"), blockName);
      lightguidecyllog->SetVisAttributes(  lgVisAtt );


      // Bar

      sprintf(blockName, "%ssolbar_%1d", barname[type], layer+1 );

      if( type != kBarVetoShort && type != kBarVetoLong ){
	barl = Y[type]+2.0*lgLength[type]+2.0*lgcylLen[type]+2.0*glasslen+2.0*cathlen;
	solidBar = new G4Box(blockName, X[type]/2.0, barl/2.0, pmtRad[layer]*2.0 > Z[type]?  pmtRad[layer] : Z[type]/2.0 );
      } else {
	barl = Y[type]+lgLength[type]+lgcylLen[type]+glasslen+cathlen;
	solidBar = new G4Box(blockName, X[type]/2.0, barl/2.0, pmtRad[layer]*2.0 > Z[type]?  pmtRad[layer] : Z[type]/2.0 );
      }

      sprintf(blockName, "%slogbar_%1d", barname[type], layer+1 );
      // Changed the material from Air to Special_Air: 
      logBar[type][layer] = new G4LogicalVolume(solidBar, GetMaterial("Special_Air"), blockName);
      logBar[type][layer]->SetVisAttributes(  G4VisAttributes::Invisible   );

      //////////////////////////////////////////////////////////////////////////////////////////////
      // GRAB THE SCINTILLATOR INFORMATION --- THESE ARE DUMMY CONTAINERS, WILL BE USED
      // TO FILL THE DETMAP VARIABLES IN A ROUTINE BELOW!!!                            
      //////////////////////////////////////////////////////////////////////////////////////////////
      // Make a "logBar" mother volume which will house the Scintillator SD, PMT SD, & LightGuides
      // Placing the bars (logBlock) inside logBar:
      sprintf(blockName, "%sphysblock_%1d", barname[type], layer+1 );
	    
      switch( type ){
      case kBarVetoLong:
	lgoffset = -(barl-Y[type])/2.0;
	new G4PVPlacement(0, G4ThreeVector(0.0, lgoffset,  0.0), logBlock[type][layer],
			  blockName, logBar[type][layer], false, type );
	scint_type_pos_logBar[type] =  G4ThreeVector(0.0, lgoffset,  0.0);
	break;

      case kBarVetoShort:
	lgoffset = (barl-Y[type])/2.0;
	new G4PVPlacement(0, G4ThreeVector(0.0, lgoffset, 0.0), logBlock[type][layer],
			  blockName, logBar[type][layer], false, type );
	scint_type_pos_logBar[type] =  G4ThreeVector(0.0, lgoffset,  0.0);
	break;

      default:
	lgoffset = 0.0;
	new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0), logBlock[type][layer],
			  blockName, logBar[type][layer], false, type );
	scint_type_pos_logBar[type] =  G4ThreeVector(0.0, lgoffset,  0.0);
	break;
      }

      // Placing the Waveguides and PMTs inside logBar:
      // Get the appropriate rotations, make two copies of PMT + Waveguides (either side)

      //////////////////////////////////////////////////////////////////////////////////////////////
      // GRAB THE PMT INFORMATION --- THESE ARE DUMMY CONTAINERS, WILL BE USED
      // TO FILL THE DETMAP VARIABLES IN A ROUTINE BELOW!!!                            
      //////////////////////////////////////////////////////////////////////////////////////////////
      G4ThreeVector pmt_pos(0.0,0.0,0.0);

      if( type != kBarVetoShort ){
	sprintf(blockName, "%sphyslg_l_%1d", barname[type], layer+1 );
	new G4PVPlacement( lgleftrot, G4ThreeVector(0.0, lgoffset+Y[type]/2.0+lgLength[type]/2.0, 0.0), 
	lightguidelog, blockName, logBar[type][layer], false, 0 );

	sprintf(blockName, "%sphyslgcyl_l_%1d", barname[type], layer+1 );
	new G4PVPlacement( lgleftrot, G4ThreeVector(0.0, lgoffset+Y[type]/2.0+lgLength[type]+lgcylLen[type]/2.0, 0.0), 
	lightguidecyllog, blockName, logBar[type][layer], false, 0 );

	sprintf(blockName, "%sphyspmt_l_%1d", barname[type], layer+1 );
	new G4PVPlacement( lgleftrot, G4ThreeVector(0.0, lgoffset+Y[type]/2.0+lgLength[type]+lgcylLen[type]+(cathlen+glasslen)/2.0, 0.0), 
	logicPMT[layer], blockName, logBar[type][layer], false, 0 );

	pmt_pos = G4ThreeVector(0.0, lgoffset+Y[type]/2.0+lgLength[type]+lgcylLen[type]+(cathlen+glasslen)/2.0, 0.0 );
      }

      pmt_type_pos_logBar[type].push_back(pmt_pos);

      if( type != kBarVetoLong ){
	sprintf(blockName, "%sphyslg_r_%1d", barname[type], layer+1 );
	new G4PVPlacement( lgrightrot, G4ThreeVector(0.0, lgoffset-Y[type]/2.0-lgLength[type]/2.0, 0.0), 
	lightguidelog, blockName, logBar[type][layer], false, 1 );

	sprintf(blockName, "%sphyslgcyl_r_%1d", barname[type], layer+1 );
	new G4PVPlacement( lgrightrot, G4ThreeVector(0.0, lgoffset-Y[type]/2.0-lgLength[type]-lgcylLen[type]/2.0, 0.0), 
	lightguidecyllog, blockName, logBar[type][layer], false, 1 );

	sprintf(blockName, "%sphyspmt_r_%1d", barname[type], layer+1 );
	new G4PVPlacement( lgrightrot, G4ThreeVector(0.0, lgoffset-Y[type]/2.0-lgLength[type]-lgcylLen[type]-(cathlen+glasslen)/2.0, 0.0), 
	logicPMT[layer], blockName, logBar[type][layer], false, 1 );

	pmt_pos = G4ThreeVector(0.0, lgoffset-Y[type]/2.0-lgLength[type]-lgcylLen[type]-(cathlen+glasslen)/2.0, 0.0);
      }

      pmt_type_pos_logBar[type].push_back(pmt_pos);	    
    }
  }
  // Cut up the PMT map to kill duplicates, insert into pmt_types:
  // There should be 6 types, all of which having a vector size of 2
  // except for the VETO short / long layers.
  map<int,vector<G4ThreeVector> > pmt_types;

  map<int,vector<G4ThreeVector> >::iterator tit;
  vector<G4ThreeVector>::iterator it1;

  for(tit=pmt_type_pos_logBar.begin(); tit!=pmt_type_pos_logBar.end(); tit++){
    int count = 1;
    for(it1=tit->second.begin(); it1!=tit->second.end(); it1++){
      if(count <= 2 ) {
	if( tit->first != 4 && tit->first != 5 ) {
	  pmt_types[tit->first].push_back( *it1 );
	  count++;
	} else {
	  if( tit->first == 4 ){
	    count = 3;
	    pmt_types[tit->first].push_back( *it1 );
	  } else {
	    if( count==2 ){
	      pmt_types[tit->first].push_back( *it1 );
	    }
	    count++;
	    
	  }
	}

      }
    }
  }
  // Check to make sure that pmt_types is correct - here is the output:
  // type = 0 size = 2
  //   (0,1238.08,0)
  //   (0,-1238.08,0)
  // type = 1 size = 2
  //   (0,1032.69,0)
  //   (0,-1032.69,0)
  // type = 2 size = 2
  //   (0,1032.69,0)
  //   (0,-1032.69,0)
  // type = 3 size = 2
  //   (0,1244.31,0)
  //   (0,-1244.31,0)
  // type = 4 size = 1
  //   (0,609.644,0)
  // type = 5 size = 1
  //   (0,-459.644,0)


  //////////////////////////////////////////////////////////////////////////////
  //  Cassettes
  // We will be placing iron/lead plates + logBars inside mother volumes
  // referred to as logCassette
  //////////////////////////////////////////////////////////////////////////////
  G4LogicalVolume *logplate;
  double platez;
  int ctype;

  // Get Scint Position within each Cassette:
  map<int,vector<G4ThreeVector> > scint_casstype_pos_logCass;
  map<int,vector<G4ThreeVector> > pmt_casstype_pos_logCass;
  
  for( layer = 0; layer < NLAYERS; layer++ ){
    for( ctype = 0; ctype < NCASSETTE_TYPES; ctype++ ){

      sprintf(blockName, "%ssolcassette_%1d", casname[ctype], layer+1 );
      solidBlock = new G4Box(blockName, casX[ctype]/2.0, casY[ctype]/2.0, casZ[ctype]/2.0 );
      sprintf(blockName, "%slogcassette_%1d", casname[ctype], layer+1 );
      logCassette[ctype][layer] = new G4LogicalVolume(solidBlock, GetMaterial("Air"), blockName);

      // Place plates
      // Iron 1 this is the plate closest to front of scintillator
      sprintf(blockName, "%ssolFe1casplate_%1d", casname[ctype], layer+1 );
      solidBlock = new G4Box(blockName, casX[ctype]/2.0, casFe1length[ctype]/2.0, casFe1thick[ctype]/2.0 );
      sprintf(blockName, "%slogFe1casplate_%1d", casname[ctype], layer+1 );

      // All the Veto plates are actually aluminum, despite
      // other documentation (Bogdan, Jan 8, 2009)
      logplate = 0;
      if( ctype == kCassVeto1 || ctype == kCassVeto2 ){
	logplate   = new G4LogicalVolume(solidBlock, GetMaterial("Styro"), blockName );
	logplate->SetVisAttributes( styroVisAtt );
      } else {
	logplate   = new G4LogicalVolume(solidBlock, GetMaterial("Iron"), blockName );
	logplate->SetVisAttributes( ironVisAtt );
      }

      sprintf(blockName, "%sphysFe1casplate_%1d", casname[ctype], layer+1 );

      if( ctype == kCassVeto1 || ctype == kCassVeto2 ){

	// Vetos only have this plate as padding
	platez = -casZ[ctype]/2.0 + casAlthick[ctype] + casFe1thick[ctype]/2.0 + 1.0*mm;
	// Veto 2 is backwards - styrofoam in back
	if( ctype == kCassVeto2 ){
	  platez *= -1.0;
	}
	new G4PVPlacement( 0, G4ThreeVector(0.0, casbarY[ctype], platez ), logplate,
			   blockName, logCassette[ctype][layer], false, 0 );
      } else {
	// Regular bars have plate 1 in front of bars
	platez = -casZ[ctype]/2.0 + casFe1thick[ctype]/2.0 + casFe2thick[ctype];
	new G4PVPlacement( 0, G4ThreeVector(0.0, casbarY[ctype], platez ), logplate,
			   blockName, logCassette[ctype][layer], false, 0 );
      }


      // Iron 2
      // This one has holes in it - this is for full layers ONLY
      //

      if( ctype != kCassVeto1 && ctype != kCassVeto2 ){
	double edge = 9.75*inch;
	double holewidth1, holewidth2;
	double holepos1, holepos2;

	G4Box *fullplate, *hole1, *hole2;
	// Need a full solid and then 2 representing the holes

	sprintf(blockName, "%ssolFe2casfullplate_%1d", casname[ctype], layer+1 );
	fullplate = new G4Box(blockName, casX[ctype]/2.0, casY[ctype]/2.0, casFe2thick[ctype]/2.0 );

	holewidth1 = casY[ctype]/2.0 - casbarY[ctype] - casFe2length[ctype]/2.0 - edge;
	holewidth2 = casY[ctype]/2.0 + casbarY[ctype] - casFe2length[ctype]/2.0 - edge;

	holepos1 = casY[ctype]/2.0 - edge - holewidth1/2.0;;
	holepos2 = -casY[ctype]/2.0 + edge + holewidth2/2.0;;

	sprintf(blockName, "%ssolFe2cashole1plate_%1d", casname[ctype], layer+1 );
	hole1 = new G4Box(blockName,casFe2height[ctype]/2.0,  holewidth1/2.0, casFe2thick[ctype]/2.0 );
	sprintf(blockName, "%ssolFe2cashole2plate_%1d", casname[ctype], layer+1 );
	hole2 = new G4Box(blockName,casFe2height[ctype]/2.0,  holewidth2/2.0, casFe2thick[ctype]/2.0 );
	G4ThreeVector hole1trans( holepos1, 0.0, 0.0 );
	G4ThreeVector hole2trans( holepos2, 0.0, 0.0 );

	sprintf(blockName, "%ssolFe2cas1holeplate_%1d", casname[ctype], layer+1 );
	G4SubtractionSolid *ironplate2_1hole = new G4SubtractionSolid( blockName, fullplate, hole1, 0, hole1trans );
	sprintf(blockName, "%ssolFe2casplate_%1d", casname[ctype], layer+1 );
	G4SubtractionSolid *ironplate2 = new G4SubtractionSolid( blockName, ironplate2_1hole, hole2, 0, hole2trans );
	// Now we have the solid - place
	sprintf(blockName, "%slogFe2casplate_%1d", casname[ctype], layer+1 );

	logplate   = new G4LogicalVolume(ironplate2, GetMaterial("Iron"), blockName );
	logplate->SetVisAttributes(ironVisAtt );

	platez = -casZ[ctype]/2.0 + casFe2thick[ctype]/2.0;

	new G4PVPlacement(0, G4ThreeVector( 0.0, 0.0, platez ), logplate,
			  blockName, logCassette[ctype][layer], false, 0 );

	// Aluminum - this is the back plate for normal layers

	sprintf(blockName, "%ssolAlcasplate_%1d", casname[ctype], layer+1 );
	solidBlock = new G4Box(blockName, casX[ctype]/2.0, casY[ctype]/2.0, casAlthick[ctype]/2.0 );
	sprintf(blockName, "%slogAlcasplate_%1d", casname[ctype], layer+1 );
	logplate   = new G4LogicalVolume(solidBlock, GetMaterial("Aluminum"), blockName );
	logplate->SetVisAttributes(  alVisAtt );

	// Al goes in back
	platez = casZ[ctype]/2.0 - casAlthick[ctype]/2.0;

	sprintf(blockName, "%sphysAlcasplate_%1d", casname[ctype], layer+1 );
	new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, platez), logplate,
			   blockName, logCassette[ctype][layer], false, 0 );
      }

      // This the main veto Al plate and supporting strips
      if( ctype == kCassVeto1 || ctype == kCassVeto2 ){

	sprintf(blockName, "%ssolAlcasplate2_%1d", casname[ctype], layer+1 );
	solidBlock = new G4Box(blockName, casX[ctype]/2.0, casY[ctype]/2.0, casAlthick[ctype]/2.0 );
	sprintf(blockName, "%slogAlcasplate2_%1d", casname[ctype], layer+1 );
	logplate   = new G4LogicalVolume(solidBlock, GetMaterial("Aluminum"), blockName );
	logplate->SetVisAttributes(  alVisAtt );

	// This plate is flush with plate 1 (+~mm)
	platez =  -casZ[ctype]/2.0 + casAlthick[ctype]/2.0; 
	if( ctype == kCassVeto2 ){
	  platez *= -1.0;
	}

	sprintf(blockName, "%sphysAlcasplate2_%1d", casname[ctype], layer+1 );
	new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, platez), logplate,
			   blockName, logCassette[ctype][layer], false, 0 );


	// Aluminum strips
	sprintf(blockName, "%ssolAlstrip_%1d", casname[ctype], layer+1 );
	solidBlock = new G4Box(blockName, casX[ctype]/2.0, 3.0*2.54*cm/2.0, 0.25*2.54*cm/2.0 );
	logplate   = new G4LogicalVolume(solidBlock, GetMaterial("Aluminum"), blockName );
	logplate->SetVisAttributes(  alVisAtt );

	if( ctype == kCassVeto1 ){
	  platez = casbarZ[ctype] + Z[ctype]/2.0 + 2.54*cm/8.0;
	} else {
	  //					platez *= -1.0;
	  platez = casbarZ[ctype] - Z[ctype]/2.0 - 2.54*cm/8.0;
	}

	barl = Y[kBarVetoShort] + Y[kBarVetoLong];
	double bardiff = (Y[kBarVetoLong] - Y[kBarVetoShort])/2.0;

	// Four strips
	sprintf(blockName, "%sphysAlstrip1_%1d", casname[ctype], layer+1 );
	new G4PVPlacement( 0, G4ThreeVector(0.0, barl/2.0-2.0*2.54*cm + casbarY[ctype] , platez), logplate,
			   blockName, logCassette[ctype][layer], false, 0 );

	sprintf(blockName, "%sphysAlstrip2_%1d", casname[ctype], layer+1 );
	if( ctype == kCassVeto1 ){
	  new G4PVPlacement( 0, G4ThreeVector( 0.0, casbarY[ctype] + bardiff, platez), logplate,
			     blockName, logCassette[ctype][layer], false, 0 );
	  new G4PVPlacement( 0, G4ThreeVector(0.0, casbarY[ctype] - bardiff, platez), logplate,
			     blockName, logCassette[ctype][layer], false, 0 );
	} else {
	  new G4PVPlacement( 0, G4ThreeVector(0.0, casbarY[ctype] + bardiff, platez), logplate,
			     blockName, logCassette[ctype][layer], false, 0 );
	  new G4PVPlacement( 0, G4ThreeVector(0.0, casbarY[ctype] - bardiff, platez), logplate,
			     blockName, logCassette[ctype][layer], false, 0 );
	}

	sprintf(blockName, "%sphysAlstrip3_%1d", casname[ctype], layer+1 );
	new G4PVPlacement( 0, G4ThreeVector(0.0,-barl/2.0+2.0*2.54*cm + casbarY[ctype], platez), logplate,
			   blockName, logCassette[ctype][layer], false, 0 );
      }

      double x;
      double vetolen, vetooff;
 
      /////////////////////////////////////////////////////////////////////////////////////////////////////////
      // Place bars in the CASSETTES:
      /////////////////////////////////////////////////////////////////////////////////////////////////////////

      // Grab the necessary information Sensitive Detector Information:
      G4ThreeVector scint_pos_in_bar = scint_type_pos_logBar[ctype];
      vector<G4ThreeVector> pmt_pos_in_bar = pmt_types[ctype];

      if( ctype == kCassVeto1 || ctype == kCassVeto2 ){
	for( int bar = 0; bar < Nbars[ctype]; bar++ ){

	  x = barspacing[ctype]*bar + X[kBarVetoLong]/2.0 - casX[ctype]/2.0;
	  x*=-1;
	  // Full veto length
	  vetolen = Y[kBarVetoLong] + Y[kBarVetoShort] +lgLength[kBarVetoShort] + lgLength[kBarVetoLong] + lgcylLen[kBarVetoShort] + lgcylLen[kBarVetoLong] +2.0*glasslen+2.0*cathlen;

	  barl = Y[kBarVetoShort] + lgLength[kBarVetoShort] + lgcylLen[kBarVetoShort] +  glasslen + cathlen;

	  vetooff = -(vetolen - barl)/2.0 + casbarY[ctype];

	  sprintf(blockName, "%sphysbar_%1d_%1d_%1d", barname[kBarVetoShort], layer+1, ctype, bar );
	  new G4PVPlacement( 0, G4ThreeVector(x, vetooff, casbarZ[ctype]), logBar[kBarVetoShort][layer],
	  		     blockName, logCassette[ctype][layer], false, bar );
	  
	  G4ThreeVector vetoshortpos(x, vetooff, casbarZ[ctype]);
	  scint_casstype_pos_logCass[kBarVetoShort].push_back( vetoshortpos + scint_pos_in_bar );
	  pmt_casstype_pos_logCass[kBarVetoShort].push_back( vetoshortpos + pmt_pos_in_bar[0] );
	  

	  barl = Y[kBarVetoLong] + lgLength[kBarVetoLong] + lgcylLen[kBarVetoLong] +  glasslen + cathlen;
	  vetooff = (vetolen - barl)/2.0 + casbarY[ctype];

	  sprintf(blockName, "%sphysbar_%1d_%1d_%1d", barname[kBarVetoLong], layer+1, ctype, bar );
	  new G4PVPlacement( 0, G4ThreeVector(x, vetooff, casbarZ[ctype]), logBar[kBarVetoLong][layer],
			     blockName, logCassette[ctype][layer], false, bar );

	  G4ThreeVector vetolongpos(x, vetooff, casbarZ[ctype]);
	  scint_casstype_pos_logCass[kBarVetoLong].push_back( vetolongpos + scint_pos_in_bar );
	  pmt_casstype_pos_logCass[kBarVetoLong].push_back( vetolongpos + pmt_pos_in_bar[0] );
 
	}
      } else {
	for( int bar = 0; bar < Nbars[ctype]; bar++ ){
	  x = barspacing[ctype]*bar + X[ctype]/2.0 - casX[ctype]/2.0;
	  x*=-1;
	  sprintf(blockName, "%sphysbar_%1d_%1d_%1d", barname[ctype], layer+1, ctype, bar );
	  new G4PVPlacement( 0, G4ThreeVector(x, casbarY[ctype], casbarZ[ctype]), logBar[ctype][layer],
			     blockName, logCassette[ctype][layer], false, bar );
	  scint_casstype_pos_logCass[ctype].push_back(  G4ThreeVector(x, casbarY[ctype], casbarZ[ctype]) );
	  
	  // These should have two PMTs, so we need to iterate through the vector:
	  for( it1 = pmt_pos_in_bar.begin(); it1 != pmt_pos_in_bar.end(); it1++ ) {
	    pmt_casstype_pos_logCass[ctype].push_back( G4ThreeVector(x, casbarY[ctype], casbarZ[ctype]) + *it1   );
	  }

	}
      }
 
    }
  }
  
  // There is over counting in my vectors, therefore lets remove everything and
  // use scint_logCass / pmt_logCass instead.
  // key = cassette type, value = vector of scintillator positions within the cassette
  map<int,vector<G4ThreeVector> > scint_logCass, pmt_logCass;  
  map<int,vector<G4ThreeVector> >::iterator vit;
  vector<G4ThreeVector>::iterator vit2;

  for(vit = scint_casstype_pos_logCass.begin(); vit!=scint_casstype_pos_logCass.end(); vit++){
    int counter = 1;
    for( vit2=vit->second.begin(); vit2!=vit->second.end(); vit2++ ){
      if( counter <= Nbars[vit->first]) { 
	scint_logCass[vit->first].push_back( *vit2 );
	counter++;
      }
    }   
  }
  for(vit = pmt_casstype_pos_logCass.begin(); vit!=pmt_casstype_pos_logCass.end(); vit++){
    int counter = 1;
    for( vit2=vit->second.begin(); vit2!=vit->second.end(); vit2++ ){
      if( counter <= 2.0*Nbars[vit->first]) { 
	pmt_logCass[vit->first].push_back( *vit2 );
	counter++;
      }
    }   
  }

  // Testing the output.. Looks pretty good.
  // for(vit = pmt_logCass.begin(); vit!=pmt_logCass.end(); vit++){
  //   cout << vit->first << "    total = " << (vit->second).size() << endl;
  //   for( vit2=vit->second.begin(); vit2!=vit->second.end(); vit2++ ){
  //     cout << *vit2 << endl;
  //   }
  // }
  // sample output:
  // 3    total = 10
  //   (-409.5,1218.31,13.1)
  //   (-409.5,-1270.31,13.1)
  //   (-205.7,1218.31,13.1)
  //   (-205.7,-1270.31,13.1)
  //   (-1.9,1218.31,13.1)
  //   (-1.9,-1270.31,13.1)
  //   (201.9,1218.31,13.1)
  //   (201.9,-1270.31,13.1)
  //   (405.7,1218.31,13.1)
  //   (405.7,-1270.31,13.1)

    
  // Spacer configuration
  double spacer[NLAYERS][6] =
  // Tim's document suggests there's no spacer for the bottom
  // V1 cassette
  // These numbers also include the deviations from expressed
  // cassette height
    { { 1.30*cm, 1.30*cm, 1.46*cm, 0.99*cm, 0.99*cm, 0.99*cm },
      { 0.19*cm, 1.30*cm, 1.35*cm, 1.30*cm, 1.46*cm, 1.30*cm },

      { 0.00*cm, 0.12*cm, 1.27*cm, 1.27*cm, 1.43*cm, 1.27*cm },
      { 0.48*cm, 1.75*cm, 1.43*cm, 1.27*cm, 1.27*cm, 1.27*cm },
      { 1.12*cm, 1.27*cm, 1.12*cm, 1.27*cm, 1.64*cm, 1.27*cm },
      { 0.00*cm, 0.32*cm, 1.27*cm, 1.59*cm, 0.52*cm, 1.27*cm },

      { 1.27*cm, 0.95*cm, 1.27*cm, 1.27*cm, 1.27*cm, 1.27*cm },
      { 0.00*cm, 0.95*cm, 0.32*cm, 1.43*cm, 1.27*cm, 1.27*cm },
      { 0.95*cm, 0.16*cm, 1.13*cm, 1.59*cm, 1.27*cm, 1.27*cm },
    };
		  
  double x;
  G4VPhysicalVolume *phys;

  int scint_dummy = 1;
  int pmt_dummy = 1;

  int barNo = 1;


  for( layer = 0; layer < NLAYERS; layer++ ){
    // Build layers from cassettes

    //		printf("Layer %d", layer );
    NRows[layer] = 0;
    //barno = 0;
    // This was done backwards...
    x = layerBottom[layer];
    int number = 1;
    int rowNo_pmt = 1;
    int rowNo = 1;
 

    for( int cass = 0; cass < Ncass[layer]; cass++ ){

      // We need to grab to correct cassette vector, then add
      // to our current positions
      vector<G4ThreeVector> scint_pos_cass = scint_logCass[casType[layer][cass]];
      vector<G4ThreeVector> pmt_pos_cass = pmt_logCass[casType[layer][cass]];
     
      x -= casX[casType[layer][cass]]/2.0;

      // Fronts of cassettes are flush with layer z

      sprintf(blockName, "physcas_%1d_%1d", cass, layer+1 );
      G4ThreeVector cassette_placement( x, 0.0, layerZ[layer]+casZ[casType[layer][cass]]/2.0 );
      phys = new G4PVPlacement( 0, cassette_placement,
				logCassette[casType[layer][cass]][layer], blockName, neutronarm_log, false, barNo );

      barNo += Nbars[casType[layer][cass]];

      x -= (casX[casType[layer][cass]]/2.0 + spacer[layer][cass]);
      // set up NRows here too
      NRows[layer] += Nbars[casType[layer][cass]];

      ///////////////////////////////////////////////////////////////////////
      // DETMAP 
      ///////////////////////////////////////////////////////////////////////
      // Finish up gathering information for our DetMap:
      for( vit2=scint_pos_cass.begin(); vit2!=scint_pos_cass.end(); vit2++ ){	  
	G4ThreeVector block_position = *vit2;
	// Fill the detmap for the scintillators ( Col = Scintillator type since it has no real meaning ):
	(NDScintSD->detmap).LocalCoord[scint_dummy] = block_position + cassette_placement;
	(NDScintSD->detmap).Plane[scint_dummy] = layer;
	// This is not really column, but scintillator type
	(NDScintSD->detmap).Col[scint_dummy] = casType[layer][cass];
	(NDScintSD->detmap).Row[scint_dummy] = rowNo;
       
	rowNo++;
	scint_dummy++;
      }
      for( vit2=pmt_pos_cass.begin(); vit2!=pmt_pos_cass.end(); vit2++ ){	
	G4ThreeVector block_position = *vit2;
	(NDSD->detmap).LocalCoord[pmt_dummy] = block_position + cassette_placement;
	(NDSD->detmap).Plane[pmt_dummy] = layer;
	// This is not really column, but block type
	(NDSD->detmap).Col[pmt_dummy] = casType[layer][cass];
	
	(NDSD->detmap).Row[pmt_dummy] = rowNo_pmt;

	if( number%2==0 ) {
	  rowNo_pmt++;
	}
	number++;
    
	pmt_dummy++;
      }
      
    }
  }

  // Test:
  //map<int,G4ThreeVector>::iterator vit3;
  //for(vit3 = final_scint.begin(); vit3!=final_scint.end(); vit3++){
    //cout << vit3->first << "    " << vit3->second << endl;
  //}
  // map<int,int>::iterator tet;
  // for( tet=((NDSD->detmap).Plane).begin(); tet!=((NDSD->detmap).Plane).end(); tet++ ){
  //   cout << tet->first << "     " << tet->second << endl;
  // }

  int plate;
			  // printf("Generating plates\n");
  for( plate = 0; plate < NPLATES; plate++ ){
    sprintf(blockName, "solplate_%1d",plate+1 );
    solidBlock = new G4Box(blockName,plateHeight/2.0,  plateWidth/2.0, plateThick[plate]/2.0 );

    sprintf(blockName, "logplate_%1d",plate+1 );
    logplate   = new G4LogicalVolume(solidBlock, plateMat[plate], blockName );
    
    if( plateMat[plate]==GetMaterial("Lead") ) {
      logplate->SetVisAttributes( leadVisAtt );
    } else {
      logplate->SetVisAttributes( ironVisAtt );
    }    

    phys = new G4PVPlacement( 0, G4ThreeVector(0.0, 0.0, plateZ[plate] + plateThick[plate]/2.0), logplate, blockName, neutronarm_log, false,0 );
  }

  // Cassettes are invisible
  for( layer = 0; layer < NLAYERS; layer++ ){
    for( type = kBarCMU; type < NCASSETTE_TYPES; type++ ){
      logCassette[type][layer]->SetVisAttributes( G4VisAttributes::Invisible );
    }
  }

  // Bars are yellow, vetos are red
  // for( layer = 0; layer < NLAYERS; layer++ ){
  //   for( type = kBarCMU; type < NBARTYPES; type++ ){
  //     logBar[type][layer]->SetVisAttributes(  G4VisAttributes::Invisible   );
  //     switch(type){
  //     case kBarVetoShort:
  // 	logBlock[type][layer]->SetVisAttributes( vetoVisAtt  );
  // 	break;
  //     case kBarVetoLong:
  // 	logBlock[type][layer]->SetVisAttributes( veto2VisAtt  );
  // 	break;
  //     case kBarUVA:
  // 	logBlock[type][layer]->SetVisAttributes( UVAVisAtt  );
  // 	break;
  //     case kBarCMU:
  // 	logBlock[type][layer]->SetVisAttributes( CMUVisAtt  );
  // 	break;
  //     case kBarGlasgow:
  // 	logBlock[type][layer]->SetVisAttributes( GlasVisAtt  );
  // 	break;
  //     default:
  // 	logBlock[type][layer]->SetVisAttributes( barVisAtt  );
  // 	break;
  //     }
  //   }
  // }

}
