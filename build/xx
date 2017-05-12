  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");
// Liquid argon material(Si o veux rajouter un materiel avec des carateristiques)
 /* G4double a;  // mass of a mole;
  G4double z;  // z=mean number of protons;  
  G4double density; 
 new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                kStateGas, 2.73*kelvin, 3.e-18*pascal);*/
  
  G4double a;  // atomic mass
  G4double zz;  // z=mean number of protons
  G4double density;
  G4int ncomponents;


  /*G4Element* elPb = new G4Element("Lead","Pb",zz=82.,a=207.2*g/mole);
  G4Element* elF = new G4Element("Fluor","F",zz=9.,a=18.998*g/mole);

  G4Material* PbF2 = new G4Material("PbF2",density=7.77*g/cm3,ncomponents=2);
  PbF2->AddElement(elPb,1);
  PbF2->AddElement(elF, 2);*/

  // Geometry parameters
  //G4int nofLayers = 1;
  //G4double absoThickness = 10.*mm;
  //G4double gapThickness =  5.*mm;

  G4double calorThickness =18.6*cm;
  G4double calorSizeXY  = 3.5*cm;

 /////ooooooooooooooooMother volume de calorimetreoooooooooooooooooooooooooo

  G4double worldSizeXY = 100. * calorSizeXY;// 
  G4double worldSizeZ  = 200. * calorThickness; 
  
  // Get materials
 auto defaultMaterial = G4Material::GetMaterial("Air");
  /*auto absorberMaterial = G4Material::GetMaterial("G4_Pb");
  auto gapMaterial = G4Material::GetMaterial("liquidArgon");
  */
  if ( ! defaultMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined."; 
    G4Exception("G4RTPC::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }  
   
  //     
  // World
  //
  /*G4Box *worldS 
    = new G4Box("World",           // its name
                 worldSizeXY/2.0, worldSizeXY/2.0, worldSizeZ/2.0); // its size
                         
  G4LogicalVolume *worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name
  G4double x = -5.0*m;//* std::sin(fCalAngle);
  G4double z = -5.0*m;//* std::cos(fCalAngle); 
                            
  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(x,0.,z),  // at (0,0,0)
                 worldLV,          // its logical volume                         
                 "World",          // its name
                 worldlog,         // its mother  volume*************
                 false,            // no boolean operation
                 0,                // copy number
                 0);  // checking overlaps 
   */
  ////////////////////////////////////////////////////////        



  //fMaw = NULL;//*********
  //  fPWT=NULL;
 // fIsInteractive=1;//*********
  //fIsSrcPb = false;//*********
  //fIsOverlapVol = false;//*********
  //fNtarget = 1;//*********
  //fNgas = 2;//*********
  /*double fXmin = fBmin = fBmax = fBz = 0.0;
  double fTx = fTy = fTz = fRst = fTst = fZst = fRHe1 = fRHe1a = fRHe2 = fRw = fRbl =
    fTbl = fSbl = fZbl =
    double fRbaf = fTbaf = fRG = fTG = fTend1 = fTend2 = fTsh = fZsh = 0;
   double fNwI = fNwO = fNbl = 0;
 double fShI = fShO = fShZI = fShZO = fShTh = 0.0;
 //fFieldMap = NULL;
  //fBFieldType = EUniField;
  double fBScaleFac = 1.0;
 // fMst = 0; // 0 = kapton, 1 = Al
  double fOvHe = 0.0;
 //fWMat = 0;*/
////////////////////////////:parameters moher GEM (RTPCContainer)///////////////////
        double mD2GasL = 380*mm;
        double mPCBReadOutR = 87.7936*mm;
        double mBedPlateHighEdge=105.0*mm;
	double RTPCContainerL = mD2GasL+120*mm;       
	double RTPCContainerRin = 0.0*m; 
	double RTPCContainerRout = max(mPCBReadOutR+50*mm,mBedPlateHighEdge+12*mm); 

  // Cylindrical triple GEM chamber surrounding He gas volume of RTPC
  //
 G4double  fRHe2=0;
  G4double  fZHe = RTPCContainerL; //**********************
  G4double rg = fRHe2;       // start at outer radius of He volume
 /* G4double tgkap = 0.05*mm;  // 50um kapton GEM foils
  G4double tgkap1 = 0.2*mm;  // 200um kapton readout foil
  G4double tgcu = 0.005*mm;  // 5 um Cu cladding of foils
  G4double tgap = 2.0*mm;    // radial spacing between foils
  G4double rginner = rg;     // inner GEM radius
  G4double rgouter;
 
  // Inner GEM                                                // Rout:
  G4Tubs* G1a = new G4Tubs("G1a",rg,rg+tgcu,fZHe,0.0,360*deg);//tgcu=0.005
  G4LogicalVolume* LG1a =
    new G4LogicalVolume(G1a,GetMaterial("Copper"),"LG1a",0,0,0);
  rg += tgcu;
  G4Tubs* G1b = new G4Tubs("G1b",rg,rg+tgkap,fZHe,0.0,360*deg);//0.005+0.05
  G4LogicalVolume* LG1b =
    new G4LogicalVolume(G1b,GetMaterial("Kapton"),"LG1b",0,0,0);
  rg += tgkap;
  G4Tubs* G1c = new G4Tubs("G1c",rg,rg+tgcu,fZHe,0.0,360*deg);//0.055+0.005 = 0.06
  G4LogicalVolume* LG1c =
    new G4LogicalVolume(G1c,GetMaterial("Copper"),"LG1c",0,0,0);
  LG1a->SetVisAttributes(G4VisAttributes::Invisible);
  LG1b->SetVisAttributes(G4VisAttributes::Invisible);
  LG1c->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
  // Middle GEM                                               // Rout:
  rg = fRHe2 + tgap;  // rg = 2.mm
  G4Tubs* G2a = new G4Tubs("G2a",rg,rg+tgcu,fZHe,0.0,360*deg);//2.005 
  G4LogicalVolume* LG2a =
    new G4LogicalVolume(G2a,GetMaterial("Copper"),"LG2a",0,0,0);
  rg += tgcu;
  G4Tubs* G2b = new G4Tubs("G2b",rg,rg+tgkap,fZHe,0.0,360*deg);//2.005+0.05  =  2.055
  G4LogicalVolume* LG2b =
    new G4LogicalVolume(G2b,GetMaterial("Kapton"),"LG2b",0,0,0);
  rg += tgkap;
  G4Tubs* G2c = new G4Tubs("G2c",rg,rg+tgcu,fZHe,0.0,360*deg);// 2.06
  G4LogicalVolume* LG2c =
    new G4LogicalVolume(G2c,GetMaterial("Copper"),"LG2c",0,0,0);
  LG2a->SetVisAttributes(G4VisAttributes::Invisible);
  LG2b->SetVisAttributes(G4VisAttributes::Invisible);
  LG2c->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
  // Outer GEM
  rg = fRHe2 + 2*tgap;    //rg = 4.0 mm
  G4Tubs* G3a = new G4Tubs("G3a",rg,rg+tgcu,fZHe,0.0,360*deg); // 4.005
  G4LogicalVolume* LG3a =
    new G4LogicalVolume(G3a,GetMaterial("Copper"),"LG3a",0,0,0);
  rg += tgcu;
  G4Tubs* G3b = new G4Tubs("G3b",rg,rg+tgkap,fZHe,0.0,360*deg);// 4.055
  G4LogicalVolume* LG3b =
    new G4LogicalVolume(G3b,GetMaterial("Kapton"),"LG3b",0,0,0);
  rg += tgkap;
  G4Tubs* G3c = new G4Tubs("G3c",rg,rg+tgcu,fZHe,0.0,360*deg);// 4.06
  G4LogicalVolume* LG3c =
    new G4LogicalVolume(G3c,GetMaterial("Copper"),"LG3c",0,0,0);
  LG3a->SetVisAttributes(G4VisAttributes::Invisible);
  LG3b->SetVisAttributes(G4VisAttributes::Invisible);
  LG3c->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));
  // Redaout Foil
  rg = fRHe2 + 3*tgap;   // rg = 6.0 mm
  G4Tubs* G4a = new G4Tubs("G4a",rg,rg+tgcu,fZHe,0.0,360*deg);// 6.005
  G4LogicalVolume* LG4a =
    new G4LogicalVolume(G4a,GetMaterial("Copper"),"LG4a",0,0,0);
  rg += tgcu;
  G4Tubs* G4b = new G4Tubs("G4b",rg,rg+tgkap1,fZHe,0.0,360*deg);//6.205
  G4LogicalVolume* LG4b =
    new G4LogicalVolume(G4b,GetMaterial("Kapton"),"LG4b",0,0,0);
  rg += tgkap1;
  G4Tubs* G4c = new G4Tubs("G4c",rg,rg+tgcu,fZHe,0.0,360*deg);//6.21
  rgouter = rg + tgcu;//6.21
  //
  G4LogicalVolume* LG4c =
    new G4LogicalVolume(G4c,GetMaterial("Copper"),"LG4c",0,0,0);
  LG4a->SetVisAttributes(G4VisAttributes::Invisible);
  LG4b->SetVisAttributes(G4VisAttributes::Invisible);
  LG4c->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));*/
  //////////////////////////////////////////////////////////////////////:                       
             //////////////////////////
	     //The mother box of GEM components
	      //////////////////////////
  double rgouter = 6.25*mm;
 double rginner = 0.*mm;
	//In this way we can place everything into the box as if placed them in the hall 
	// Size of this tub, large enough but not too big

	/*if(mSetupSolenoid==1)  {RTPCContainerRout=0.50*m; RTPCContainerL=1.2*m;}
	else if(mSetupSolenoid>=2)  {RTPCContainerRout=1.00*m; RTPCContainerL=1.6*m;} */

	G4Tubs* RTPCContainer = new G4Tubs("RTPCContainer",rginner,
		rgouter,fZHe/2,0.0*deg,360.*deg);

	G4LogicalVolume* RTPCContainerLogical = new G4LogicalVolume(RTPCContainer, 
		defaultMaterial, "logRTPCContainer", 0, 0, 0);
	//RTPCContainerLogical->SetVisAttributes(blueVisAtt);

	//the position at the hall
	double mTargetXOffset=-5.0*m;;// en principe 0 (mais juste un essai pour voir leRTPC
	double mTargetYOffset=0.0*m;
	double mTargetZOffset=5.0*m;;
	double RTPCContainerPosX=mTargetXOffset;
	double RTPCContainerPosY=mTargetYOffset;
	double RTPCContainerPosZ=mTargetZOffset;
	G4PVPlacement* phyRTPCContainer= new G4PVPlacement(0,
		G4ThreeVector(RTPCContainerPosX,RTPCContainerPosY,RTPCContainerPosZ), RTPCContainerLogical,"RTPCContainner",worldlog,0,0);
///coloration
 
///////////////////////////////////////////////////////////////////////////////////////////:
  // Place all layers of GEM foil and readout in GEM volume
 /* new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG1a,"PG1a",RTPCContainerLogical,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG1b,"PG1b",RTPCContainerLogical,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG1c,"PG1c",RTPCContainerLogical,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG2a,"PG2a",RTPCContainerLogical,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG2b,"PG2b",RTPCContainerLogical,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG2c,"PG2b",RTPCContainerLogical,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG3a,"PG3a",RTPCContainerLogical,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG3b,"PG3b",RTPCContainerLogical,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG3c,"PG3c",RTPCContainerLogical,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG4a,"PG4a",RTPCContainerLogical,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG4b,"PG4b",RTPCContainerLogical,0,0);
  new G4PVPlacement(0,G4ThreeVector(0,0,0),
		    LG4c,"PG4c",RTPCContainerLogical,0,0);
*/

