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
#include "G4Trd.hh"
#include "G4Cons.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4TwoVector.hh"
#include "G4GenericTrap.hh"
#include "G4Polycone.hh"
#include "G4ExtrudedSolid.hh"
#include "G4SBSCalSD.hh"
#include "G4SBSGEMSD.hh"


#include "G4SBSTPCTOSCAField2D.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4SBSmTPC.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4MaterialTable.hh"
#include "G4MaterialPropertiesTable.hh" 
#include "sbstypes.hh"

// *dc is the detector constructor method which carries the information
// of all detectors (check (CA) )

G4SBSmTPC::G4SBSmTPC(G4SBSDetectorConstruction *dc):G4SBSComponent(dc)
//G4SBSmTPC::G4SBSmTPC(G4LogicalVolume *physiParent)
//G4SBSmTPC::G4SBSmTPC()
{
  assert(fDetCon);
  G4cout << "<G4SBSmTPC::G4SBSmTPC>: creating mTPC" << G4endl;
  
  fTPCdesign = 0;
  
  // variables for mtpc construction
  // taken from M. carmignotto gemc mtpc implementation
  // inner electrode at r=5cm
  fmTPC_inelectrode_r = 50.0*mm; //5cm of inner electrode
  fmTPC_inelectrode_kaptonthick = 0.012*mm; //12um kapton
  fmTPC_inelectrode_authick = 0.0001*mm; //0.1um Au
  // outer electrode at r=15cm
  fmTPC_outelectrode_r = 150.0*mm; //5cm of inner electrode
  fmTPC_outelectrode_kaptonthick = 0.012*mm; //12um kapton
  fmTPC_outelectrode_authick = 0.0001*mm; //0.1um Au
  // mtpc chambers
  fmTPC_cell_len = 50.0*mm; //5cm length cells
  fmTPC_Ncells = 10; //10 cells
  // readout discs
  fmTPC_readout_thick = 0.130*mm; //130um thick "kryptonite" for readout disc
  // GEMs
  fmTPC_Ngems = 2;
  fmTPC_gem_surf1thick = 0.005*mm; // 5um copper surface
  fmTPC_gem_dielecthick = 0.05*mm; // 50um dielectric
  fmTPC_gem_surf2thick = 0.005*mm; //5um copper surface
  fmTPC_gap_readoutGEM = 0.001*mm;
  fmTPC_gap_GEMGEM = 0.001*mm;
  // HV
  fmTPC_HV_thick = 0.05*mm; // 50um say gold?
  //
  fmTPCkrypto = true;//false;//by default
  fChkOvLaps = false;//true;//

  fMLPC_NplaneWires = 150;
  fMLPC_NWiresPerHalfPlane = 130;
  fMLPC_InterWireGap = 1.0*mm;
  fMLPC_InnerWireDistance = 2.0*cm;
  fMLPC_SupportThickness = 1.6*mm;
  fMLPC_RinSupport = 15.*cm;
  fMLPC_RoutSupport = 16.*cm;
  fMLPC_CathodeRin = 1.5*cm;
  fMLPC_CathodeMylarThickness = 0.001*mm;//1.0 um Mylar
  fMLPC_CathodeAlThickness =    0.00001*mm; //0.01um Al;
  fMLPC_WirePlaneThickness =    0.008*mm;//8 um wires;
  

  //step user limit of the particle in the gas detector 
  // -1*mm -> default G4 step (no user limit)
  // user step size should be greater than 0*mm!
  fmTPCstep = -1*mm;

    if(fmTPCstep != -1*mm)
      {
	G4cout<<"****************************************************************************************"<<G4endl;
	G4cout<<"            WARNING!! STEP LIMIT IN DRIFT GAS: "<< fmTPCstep<< " mm"<<G4endl;
	G4cout<<"****************************************************************************************"<<G4endl;
      }

}

void G4SBSmTPC::BuildComponent(G4LogicalVolume *motherlog){
  // Check if "Kryptonite" is activated
  
  if(fTPCdesign==1){//MLPC design
    G4double gap = 0.0*mm;// small gap...
    G4double MLPC_z_total = (fMLPC_SupportThickness*2+fMLPC_CathodeMylarThickness+fMLPC_CathodeAlThickness*2+fMLPC_WirePlaneThickness+gap*4)*fMLPC_NplaneWires+fMLPC_CathodeMylarThickness+fMLPC_CathodeAlThickness*2+fMLPC_SupportThickness*2+gap*4;

    G4cout << "MLPC total length = " << MLPC_z_total/cm << " cm " << G4endl; 
    
    // make a mother shell for mtpc filled with the drift gas
    G4Tubs* MLPCmother_solid = 
      new G4Tubs("MLPCmother_solid", 
		 fDetCon->fTargetBuilder->GetTargDiameter()/2,
		 fMLPC_RoutSupport,
		 MLPC_z_total*0.5, 
		 0.*deg, 
		 360.*deg );
  
    G4LogicalVolume* MLPCmother_log = 
      new G4LogicalVolume(MLPCmother_solid, 
			  GetMaterial("mTPCgas"),
			  "MLPCmother_log");
    
    new G4PVPlacement(0, 
		      G4ThreeVector(0.0, 0.0, fZpos), 
		      MLPCmother_log,
		      "MLPCmother_phys", 
		      motherlog, 
		      false, 
		      0, 
		      fChkOvLaps);

    G4double zpos = -MLPC_z_total*0.5 + fMLPC_SupportThickness*0.5+gap;
    G4double theta = 0.0*deg;
    G4cout << " Z position = " << zpos/mm << " mm " << endl;
    //starting with one G10 support
    BuildMLPCG10Support(MLPCmother_log, zpos);
    //loop on number of wire planes:
    // 1 layer of cathode;
    // 1 support;
    // 1 layer of wires, rotated by an angle theta;
    // 1 support;
    for(int k = 0; k<fMLPC_NplaneWires; k++){
      zpos+= (fMLPC_SupportThickness+fMLPC_CathodeMylarThickness)*0.5+fMLPC_CathodeAlThickness+gap;
      G4cout << " Z position = " << zpos/mm << " mm " << endl;
      BuildMLPCCathode(MLPCmother_log, zpos);
      zpos+= (fMLPC_CathodeMylarThickness+fMLPC_SupportThickness)*0.5+fMLPC_CathodeAlThickness+gap;
      G4cout << " Z position = " << zpos/mm << " mm " << endl;
      BuildMLPCG10Support(MLPCmother_log, zpos);
      zpos+= (fMLPC_SupportThickness+fMLPC_WirePlaneThickness)*0.5+gap;
      G4cout << " Z position = " << zpos/mm << " mm " << endl;
      BuildMLPCWires(MLPCmother_log, zpos, theta, k);
      theta+= 120.0*deg;
      zpos+= (fMLPC_WirePlaneThickness+fMLPC_SupportThickness)*0.5+gap;
      G4cout << " Z position = " << zpos/mm << " mm " << endl;
      BuildMLPCG10Support(MLPCmother_log, zpos);
    }
    //after building all wires, capping with a cathode plane and a G10 support plane
    zpos+= (fMLPC_SupportThickness+fMLPC_CathodeMylarThickness)*0.5+fMLPC_CathodeAlThickness+gap;
    G4cout << " Z position = " << zpos/mm << " mm " << endl;
    BuildMLPCCathode(MLPCmother_log, zpos);
    zpos+= (fMLPC_CathodeMylarThickness+fMLPC_SupportThickness)*0.5+fMLPC_CathodeAlThickness+gap;
    G4cout << " Z position = " << zpos/mm << " mm " << endl;
    BuildMLPCG10Support(MLPCmother_log, zpos);
    
    MLPCmother_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
  }else{
    if(fmTPCkrypto)G4cout << "Superman is dead!" << G4endl;

    // Montgomery July 2018
    // implementing geometry atm as a exact copy of implementation in gemc by park/carmignotto
    // There is no end cap on mTPC beyond readout discs and also no cathode planes
    // variables for building detector
    // total length
    double mTPC_z_total =  fmTPC_cell_len * fmTPC_Ncells;
    // centre of 1st cell
    double mTPC_centre_cell1 = -1.0 * fmTPC_cell_len * (fmTPC_Ncells-1) / 2.0;
    // radii of disks
    // inner is inner radius plus the inner electrode material and equivalent for outer
    double mTPC_rIN = fmTPC_inelectrode_r + fmTPC_inelectrode_kaptonthick + fmTPC_inelectrode_authick;
    double mTPC_rOUT = fmTPC_outelectrode_r - fmTPC_outelectrode_kaptonthick - fmTPC_outelectrode_authick;
   
    // make a mother shell for mtpc filled with the drift gas
    G4Tubs* mTPCmother_solid = 
      new G4Tubs("mTPCmother_solid", 
		 fDetCon->fTargetBuilder->GetTargDiameter()/2,
		 fmTPC_outelectrode_r, 
		 mTPC_z_total/2.0, 
		 0.*deg, 
		 360.*deg );
  
    G4LogicalVolume* mTPCmother_log = 
      new G4LogicalVolume(mTPCmother_solid, 
			  GetMaterial("ref4He"),
			  "mTPCmother_log");
  
    new G4PVPlacement(0, 
		      G4ThreeVector(0.0, 0.0, fZpos), 
		      mTPCmother_log,
		      "mTPCmother_phys", 
		      motherlog, 
		      false, 
		      0, 
		      fChkOvLaps);
   
  
    // G4String mInnerGasSDname = "SBS/InnerGas";
    // G4String mInnerGasSDcolname = "mInnerGasHitsCollection";

    // G4SBSmTPCSD *mInnerGasSD;
    // if( !(mInnerGasSD = (G4SBSmTPCSD*) fDetCon->fSDman->FindSensitiveDetector(mInnerGasSDname)) ){ //Make sure SD with this name doesn't already exist
    //   mInnerGasSD = new G4SBSmTPCSD( mInnerGasSDname, mInnerGasSDcolname );
    //   fDetCon->fSDman->AddNewDetector(mInnerGasSD);
    //   (fDetCon->SDlist).insert(mInnerGasSDname);
    //   fDetCon->SDtype[mInnerGasSDname] = G4SBS::kmTPC;
    // }
  
    G4Tubs* mInnerGas_solid = 
      new G4Tubs("mInnerGas_solid", 
		 fDetCon->fTargetBuilder->GetTargDiameter()/2,
		 fmTPC_inelectrode_r, 
		 mTPC_z_total/2.0, 
		 0.*deg, 
		 360.*deg );
  
    G4LogicalVolume* mInnerGas_log = 
      new G4LogicalVolume(mInnerGas_solid, 
			  GetMaterial("ref4He"),
			  "mInnerGas_log");
    //mInnerGas_log->SetSensitiveDetector(mInnerGasSD);  
  
    new G4PVPlacement(0, 
		      G4ThreeVector(0.0, 0.0, fZpos), 
		      mInnerGas_log,
		      "mInnerGas_phys", 
		      mTPCmother_log, 
		      false, 
		      0, 
		      fChkOvLaps);
  
   
    //Call for the construction of the different mTPC parts

    G4cout << "<G4SBSmTPC::G4SBSmTPC>: BuildmTPCWalls -> Rin = " << mTPC_rIN/mm << " , Rout " << mTPC_rOUT/mm << G4endl;  
    // make the field electrodes and boundary walls at the inner and outer radii
    BuildmTPCWalls(mTPCmother_log, mTPC_z_total, fZpos, mTPC_rIN, mTPC_rOUT);

    G4cout << "<G4SBSmTPC::G4SBSmTPC>: BuildmTPCReadouts" << G4endl;
    // build the readout discs and the gap between readout disc and gems (1 per cell)
    BuildmTPCReadouts(mTPCmother_log, mTPC_centre_cell1, fmTPC_cell_len, mTPC_rIN,  mTPC_rOUT);//, mTPCReadoutSD);

    G4cout << "<G4SBSmTPC::G4SBSmTPC>: BuildmTPCGEMs" << G4endl;
    // build the gem detectors
    BuildmTPCGEMs(mTPCmother_log, mTPC_centre_cell1, fmTPC_cell_len, mTPC_rIN,  mTPC_rOUT);

    G4cout << "<G4SBSmTPC::G4SBSmTPC>: BuildmTPCGasCells" << G4endl;
    // build the sensitive gas cells
    BuildmTPCGasCells(mTPCmother_log, mTPC_centre_cell1, fmTPC_cell_len, mTPC_rIN,  mTPC_rOUT);//, mTPCSD);//, mTPCHVSD);
   
 
 
    // Visualization attributes
    G4VisAttributes *tgt_mTPCmother_visatt = new G4VisAttributes( G4Colour( 0.0, 1.0, 1.0 ) );
    tgt_mTPCmother_visatt->SetForceWireframe(true);
    //mTPCmother_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
    mTPCmother_log->SetVisAttributes( tgt_mTPCmother_visatt );  
  
    // TPCinnergas_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
   
    // G4VisAttributes *tpcwalls_visatt = new G4VisAttributes( G4Colour( 1.0, 1.0, 1.0 ) );
    // tpcwalls_visatt->SetForceWireframe(true);
    // TPCinnerwall_log->SetVisAttributes( tpcwalls_visatt );
    // TPCouterwall_log->SetVisAttributes( tpcwalls_visatt );
   
    // // G4VisAttributes *tpcgas_visatt = new G4VisAttributes( G4Colour( 1.0, 1.0, 0.0, 0.1 ) );
    // // TPCgas_log->SetVisAttributes( tpcgas_visatt );
  }
  
}
 
void G4SBSmTPC::BuildmTPCWalls(G4LogicalVolume *motherlog, G4double mtpctotallength, G4double mtpczpos, G4double mtpcinnerR, G4double mtpcouterR){
  // make inner and outer boundary layers
  // these comprise outer kapton walls and electrodes on inner and outer raddii to set up field
   
  // inner wall, has two layers: kapton and electrode
  // layer closest to 40cm target, is kapton

  G4Tubs* mTPCinnerwall1_solid = 
    new G4Tubs("mTPCinnerwall1_solid", 
	       fmTPC_inelectrode_r, 
	       fmTPC_inelectrode_r+fmTPC_inelectrode_kaptonthick, 
	       mtpctotallength/2.0, 
	       0.*deg, 
	       360.*deg );

  G4LogicalVolume* mTPCinnerwall1_log = 
    new G4LogicalVolume(mTPCinnerwall1_solid, 
			 GetMaterial("Kapton"),
			"mTPCinnerwall1_log");

  new G4PVPlacement(0, 
		    G4ThreeVector(0.0, 0.0, mtpczpos), 
		    mTPCinnerwall1_log,
		    "mTPCinnerwall1_phys", 
		    motherlog, 
		    false, 
		    0, 
		    fChkOvLaps);
   
  // next layer before tpc gas is electrode

  G4Tubs* mTPCinnerwall2_solid = 
    new G4Tubs("mTPCinnerwall2_solid", 
	       fmTPC_inelectrode_r+fmTPC_inelectrode_kaptonthick, 
	       mtpcinnerR, mtpctotallength/2.0,
	       0.*deg, 
	       360.*deg );

  G4LogicalVolume* mTPCinnerwall2_log = 
    new G4LogicalVolume(mTPCinnerwall2_solid, 
			 GetMaterial("Au"),
			"mTPCinnerwall2_log");

  new G4PVPlacement(0, 
		    G4ThreeVector(0.0, 0.0, mtpczpos), 
		    mTPCinnerwall2_log,
		    "mTPCinnerwall2_phys", 
		    motherlog, 
		    false, 
		    0, 
		    fChkOvLaps);
   
  // outer wall, has two layers: electrode and kapton
  // layer closest to inner gas is electrode

  G4Tubs* mTPCouterwall1_solid = 
    new G4Tubs("mTPCouterwall1_solid", 
	       mtpcouterR, mtpcouterR+fmTPC_outelectrode_authick, 
	       mtpctotallength/2.0, 
	       0.*deg, 
	       360.*deg );

  G4LogicalVolume* mTPCouterwall1_log = 
    new G4LogicalVolume(mTPCouterwall1_solid, 
			 GetMaterial("Au"),
			"mTPCouterwall1_log");
  
  new G4PVPlacement(0, 
		    G4ThreeVector(0.0, 0.0, mtpczpos), 
		    mTPCouterwall1_log,
		    "mTPCouterwall1_phys", 
		    motherlog, 
		    false, 
		    0, 
		    fChkOvLaps);
   
  //outer most layer is kapton

  G4Tubs* mTPCouterwall2_solid = 
    new G4Tubs("mTPCouterwall2_solid", 
	       mtpcouterR+fmTPC_outelectrode_authick, 
	       fmTPC_outelectrode_r, 
	       mtpctotallength/2.0, 
	       0.*deg, 
	       360.*deg );

  G4LogicalVolume* mTPCouterwall2_log = 
    new G4LogicalVolume(mTPCouterwall2_solid, 
			 GetMaterial("Kapton"),
			"mTPCouterwall2_log");
  
  new G4PVPlacement(0, 
		    G4ThreeVector(0.0, 0.0, mtpczpos), 
		    mTPCouterwall2_log,
		    "mTPCouterwall2_phys", 
		    motherlog, 
		    false, 
		    0, 
		    fChkOvLaps);

  if(fmTPCkrypto)mTPCouterwall2_log->SetUserLimits( new G4UserLimits( 0.0, 0.0, 0.0, DBL_MAX, DBL_MAX ) );
  
  //Visualization attributes:
  G4VisAttributes *mtpc_kaptonboudary_visatt = new G4VisAttributes( G4Colour( 0.0, 1.0, 1.0) );
  mtpc_kaptonboudary_visatt->SetForceWireframe(true);
  mTPCinnerwall1_log->SetVisAttributes( mtpc_kaptonboudary_visatt );
  mTPCouterwall2_log->SetVisAttributes( mtpc_kaptonboudary_visatt );
   
   
  G4VisAttributes *mtpc_electrodeboudary_visatt = new G4VisAttributes( G4Colour( 1.0, 0.0, 0.0) );
  mtpc_kaptonboudary_visatt->SetForceWireframe(true);
  mTPCinnerwall2_log->SetVisAttributes( mtpc_electrodeboudary_visatt );
  mTPCouterwall1_log->SetVisAttributes( mtpc_electrodeboudary_visatt );
   
   
}
 

void G4SBSmTPC::BuildmTPCReadouts(G4LogicalVolume *motherlog, G4double centrecell1, G4double celllength, G4double innerR,  G4double outerR)
{//, G4SBSmTPCSD* mtpcreadoutSD){
  //build readout discs, one per cell, even numbered cells have it on the "LHS", odd ones on "RHS"

  printf( "G4SBSmTPC::BuildmTPCReadouts: Entering...\n");

  G4String mTPCReadoutSDname = "SBS/mTPCReadout";
  G4String mTPCReadoutcolname = "mTPCReadoutHitsCollection";
  
  // G4SBSCalSD *mTPCReadoutSD;
  // if( !(mTPCReadoutSD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(mTPCReadoutSDname)) ){ //Make sure SD with this name doesn't already exist
  //   mTPCReadoutSD = new G4SBSCalSD( mTPCReadoutSDname, mTPCReadoutcolname );
  //   fDetCon->fSDman->AddNewDetector(mTPCReadoutSD);
  //   (fDetCon->SDlist).insert(mTPCReadoutSDname);
  //   fDetCon->SDtype[mTPCReadoutSDname] = kCAL;
  // }

  
  G4Tubs* mTPCReadoutDisc_solid; 
  G4LogicalVolume* mTPCReadoutDisc_log;
  G4Tubs* mTPCReadoutGEMGap_solid; 
  G4LogicalVolume* mTPCReadoutGEMGap_log;

  G4VisAttributes *mtpc_readout_visatt = new G4VisAttributes( G4Colour( 1.0, 0.0, 1.0) );
  G4VisAttributes *mtpc_readoutgemgap_visatt = new G4VisAttributes( G4Colour( 1.0, 0.0, 1.0) );
  mtpc_readoutgemgap_visatt->SetForceWireframe(true);
  mtpc_readout_visatt->SetForceWireframe(true);

  // loop over each cell/chamber of mTPC
  for(int incCell=0; incCell<fmTPC_Ncells; incCell++)
    {
      double mTPC_CentreCell = centrecell1 + incCell*celllength;
      // make the readout discs
      double mTPC_zpos = 0.0;
      
      if(incCell % 2 == 0)
	{
	  mTPC_zpos = mTPC_CentreCell - (celllength/2.0) + (fmTPC_readout_thick/2.0);
	}
      else
	{
	  mTPC_zpos = mTPC_CentreCell + (celllength/2.0) - (fmTPC_readout_thick/2.0);
	}
      mTPCReadoutDisc_solid = 
	new G4Tubs("mTPCReadoutDisc_solid", 
		   innerR, outerR, 
		   fmTPC_readout_thick/2.0, 
		   0.*deg, 
		   360.*deg);

      mTPCReadoutDisc_log = 
	new G4LogicalVolume(mTPCReadoutDisc_solid, 
			     GetMaterial("BonusPCB"),
			    "mTPCReadoutDisc_log");    
      
      // FOR MOMENT PUT AS BONUS PCB MATERIAL, TOOK FROM MATERIALS IN GEMC
      if(fmTPCkrypto)mTPCReadoutDisc_log->SetUserLimits( new G4UserLimits( 0.0, 0.0, 0.0, DBL_MAX, DBL_MAX ) );
     
      new G4PVPlacement(0, 
			G4ThreeVector(0.0, 0.0, mTPC_zpos), 
			mTPCReadoutDisc_log, 
			"mTPCReadoutDisc_phys", 
			motherlog, 
			false, 
			incCell, 
			fChkOvLaps);

      mTPCReadoutDisc_log->SetVisAttributes( mtpc_readout_visatt );
      // set readout disc as sensitive
      //   mTPCReadoutDisc_log->SetSensitiveDetector(mTPCReadoutSD);
      
      // now we want to make a gap between readout disc and where gem will go, material should be same as mTPC gas
      double mTPC_zposgap = 0.0;
      double mTPC_edgecell = 0.0;
      
      
      if(incCell % 2 == 0)
	{
	  mTPC_edgecell = mTPC_CentreCell - celllength/2.0;
	  mTPC_zposgap = mTPC_edgecell + fmTPC_readout_thick + fmTPC_gap_readoutGEM/2.0;
	}
      else
	{
	  mTPC_edgecell = mTPC_CentreCell + celllength/2.0;
	  mTPC_zposgap = mTPC_edgecell - fmTPC_readout_thick - fmTPC_gap_readoutGEM/2.0;
	}
      
      
      //    G4cout<<" innerR: "<< innerR<< " outerR: "<<outerR<< " fmTPC_gap_readoutGEM/2.0: "<<fmTPC_gap_readoutGEM/2.0<<G4endl;
      
      mTPCReadoutGEMGap_solid = 
	new G4Tubs("mTPCReadoutGEMGap_solid", 
		   innerR, 
		   outerR, 
		   fmTPC_gap_readoutGEM/2.0, 
		   0.*deg, 
		   360.*deg);
      
      mTPCReadoutGEMGap_log = 
	new G4LogicalVolume(mTPCReadoutGEMGap_solid, 
			     GetMaterial("ref4He"),
			    "mTPCReadoutGEMGap_log");    
      
      new G4PVPlacement(0, 
			G4ThreeVector(0.0, 0.0, mTPC_zposgap), 
			mTPCReadoutGEMGap_log, 
			"mTPCReadoutGEMGap_phys", 
			motherlog, 
			false, 
			incCell, 
			fChkOvLaps);
      
      mTPCReadoutGEMGap_log->SetVisAttributes( mtpc_readoutgemgap_visatt );
      
    }//loop over mTPC cells/chambers
  
  printf("G4SBSmTPC::BuildmTPCReadouts: Leaving...\n");
 
}

void G4SBSmTPC::BuildmTPCGEMs(G4LogicalVolume *motherlog, G4double centrecell1, G4double celllength, G4double mtpcinnerR, G4double mtpcouterR)
{
  printf("G4SBSmTPC::BuildmTPCGEMs: Entering...\n");
  
  G4Tubs* mTPCGEMfoil_solid;
  G4LogicalVolume* mTPCGEMfoil_log;
  
  G4Tubs* mTPCGEMSurf1_solid; 
  G4LogicalVolume* mTPCGEMSurf1_log;
  G4Tubs* mTPCGEMDielec_solid; 
  G4LogicalVolume* mTPCGEMDielec_log;
  G4Tubs* mTPCGEMSurf2_solid; 
  G4LogicalVolume* mTPCGEMSurf2_log;
  G4Tubs* mTPCGEMGap_solid; 
  G4LogicalVolume* mTPCGEMGap_log;

  G4VisAttributes *mtpc_gem_visatt = new G4VisAttributes( G4Colour( 0.0, 1.0, 0.0) );
  mtpc_gem_visatt->SetForceWireframe(true);
  G4VisAttributes *mtpc_gemgap_visatt = new G4VisAttributes( G4Colour( 1.0, 1.0, 1.0) );
  mtpc_gemgap_visatt->SetForceWireframe(true);

  int counter = -1;

  G4String mTPCGEMfoil_solidname = "mTPCGEMfoil_solid";
  G4String mTPCGEMfoil_logname = "mTPCGEMfoil_log";
  G4String mTPCGEMfoil_physname = "mTPCGEMfoil_phys";

  mTPCGEMfoil_solid = 
    new G4Tubs(mTPCGEMfoil_solidname, 
	       mtpcinnerR, 
	       mtpcouterR, 
	       (fmTPC_gem_surf1thick+fmTPC_gem_dielecthick+fmTPC_gem_surf2thick)/2.0, 
	       0.*deg, 
	       360.*deg);

  mTPCGEMfoil_log = 
    new G4LogicalVolume(mTPCGEMfoil_solid, 
			 GetMaterial("Air"),
			mTPCGEMfoil_logname);    

  double zpossurf1 = -(fmTPC_gem_surf1thick+fmTPC_gem_dielecthick)/2.;

  G4String mTPCGEMSurf1_solidname = "mTPCGEMSurf1_solid";
  G4String mTPCGEMSurf1_logname = "mTPCGEMSurf1_log";
  G4String mTPCGEMSurf1_physname = "mTPCGEMSurf1_phys";

  mTPCGEMSurf1_solid = 
    new G4Tubs(mTPCGEMSurf1_solidname, 
	       mtpcinnerR, 
	       mtpcouterR, 
	       fmTPC_gem_surf1thick/2.0, 
	       0.*deg, 
	       360.*deg);

  mTPCGEMSurf1_log = 
    new G4LogicalVolume(mTPCGEMSurf1_solid, 
			 GetMaterial("Copper"),
			mTPCGEMSurf1_logname);    

  //place "surf1" into "foil"

  new G4PVPlacement(0, 
		    G4ThreeVector(0.0, 0.0, zpossurf1), 
		    mTPCGEMSurf1_log, 
		    mTPCGEMSurf1_physname, 
		    mTPCGEMfoil_log, 
		    false,
		    0, 
		    fChkOvLaps); 
   
  mTPCGEMSurf1_log->SetVisAttributes( mtpc_gem_visatt );
  
  double zposdielec = 0.0;
  G4String mTPCGEMDielec_solidname = "mTPCGEMDielec_solid";
  G4String mTPCGEMDielec_logname = "mTPCGEMDielec_log";
  G4String mTPCGEMDielec_physname = "mTPCGEMDielec_phys";

  mTPCGEMDielec_solid = 
    new G4Tubs(mTPCGEMDielec_solidname, 
	       mtpcinnerR, 
	       mtpcouterR, 
	       fmTPC_gem_dielecthick/2.0, 
	       0.*deg, 
	       360.*deg);

  mTPCGEMDielec_log = 
    new G4LogicalVolume(mTPCGEMDielec_solid, 
			 GetMaterial("Kapton"),
			mTPCGEMDielec_logname);    
  //place "dielec" into "foil"

  new G4PVPlacement(0, 
		    G4ThreeVector(0.0, 0.0, zposdielec), 
		    mTPCGEMDielec_log, 
		    mTPCGEMDielec_physname, 
		    mTPCGEMfoil_log, 
		    false,
		    0, 
		    fChkOvLaps);

  mTPCGEMDielec_log->SetVisAttributes( mtpc_gem_visatt );
  
  double zpossurf2 = +(fmTPC_gem_surf2thick+fmTPC_gem_dielecthick)/2.;
  G4String mTPCGEMSurf2_solidname = "mTPCGEMSurf2_solid";
  G4String mTPCGEMSurf2_logname = "mTPCGEMSurf2_log";
  G4String mTPCGEMSurf2_physname = "mTPCGEMSurf2_phys";

  mTPCGEMSurf2_solid = 
    new G4Tubs(mTPCGEMSurf2_solidname, 
	       mtpcinnerR, 
	       mtpcouterR, 
	       fmTPC_gem_surf2thick/2.0, 
	       0.*deg, 
	       360.*deg);

  mTPCGEMSurf2_log = 
    new G4LogicalVolume(mTPCGEMSurf2_solid, 
			 GetMaterial("Copper"),
			mTPCGEMSurf2_logname);
    
  //place "surf2" into "foil"
  new G4PVPlacement(0, 
		    G4ThreeVector(0.0, 0.0, zpossurf2), 
		    mTPCGEMSurf2_log, 
		    mTPCGEMSurf2_physname, 
		    mTPCGEMfoil_log, 
		    false,
		    0, 
		    fChkOvLaps);    

  mTPCGEMSurf2_log->SetVisAttributes( mtpc_gem_visatt );

  if(fmTPCkrypto)mTPCGEMfoil_log->SetUserLimits( new G4UserLimits( 0.0, 0.0, 0.0, DBL_MAX, DBL_MAX ) );
    
  // loop over each cell/chamber of mTPC
  for(int incCell=0; incCell<fmTPC_Ncells; incCell++){

    double mTPC_CentreCell = centrecell1 + incCell*celllength;
    double mTPC_edgecell = 0.0;
    if(incCell % 2 == 0){
      mTPC_edgecell = mTPC_CentreCell - celllength/2.0;
    }
    else
    {
      mTPC_edgecell = mTPC_CentreCell + celllength/2.0;
    }
      
  // now loop over how many GEMs per cell
  for(int incGEM=0; incGEM<fmTPC_Ngems; incGEM++)
    {
      counter++;
      //first, place the "GEM foil" - the master volume which will contain all the others - easier to handle stuff such as steplimiter
      double zposfoil = 0.0;
      if(incCell % 2 == 0)
	{
	  zposfoil = mTPC_edgecell + fmTPC_readout_thick + fmTPC_gap_readoutGEM + fmTPC_gem_surf1thick + fmTPC_gem_dielecthick/2.0
	    + incGEM*(fmTPC_gap_GEMGEM + fmTPC_gem_surf1thick + fmTPC_gem_dielecthick + fmTPC_gem_surf2thick);
	}
      else
	{
	  zposfoil = mTPC_edgecell - fmTPC_readout_thick - fmTPC_gap_readoutGEM - fmTPC_gem_surf1thick - fmTPC_gem_dielecthick/2.0
	    - incGEM*(fmTPC_gap_GEMGEM + fmTPC_gem_surf1thick + fmTPC_gem_dielecthick + fmTPC_gem_surf2thick);
	}
	  
      new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, zposfoil), mTPCGEMfoil_log, mTPCGEMfoil_physname, motherlog, false,counter, fChkOvLaps);
	  
      // gaps between gems, but not last one which is flush with end of cell
      if(incGEM != (fmTPC_Ngems-1))
	{
	  double zposgap = 0.0;
	  if(incCell % 2 == 0)
	    {
	      zposgap =  mTPC_edgecell + fmTPC_readout_thick + fmTPC_gap_readoutGEM + fmTPC_gem_surf1thick + fmTPC_gem_dielecthick + fmTPC_gem_surf2thick + fmTPC_gap_GEMGEM/2.0
		+ incGEM*(fmTPC_gap_GEMGEM + fmTPC_gem_surf1thick + fmTPC_gem_dielecthick + fmTPC_gem_surf2thick);
	    }
	  else
	    {
	      zposgap =  mTPC_edgecell - fmTPC_readout_thick - fmTPC_gap_readoutGEM - fmTPC_gem_surf1thick - fmTPC_gem_dielecthick - fmTPC_gem_surf2thick - fmTPC_gap_GEMGEM/2.0
		- incGEM*(fmTPC_gap_GEMGEM + fmTPC_gem_surf1thick + fmTPC_gem_dielecthick + fmTPC_gem_surf2thick);
	    }
	  G4String mTPCGEMGap_solidname = "mTPCGEMGap_solid";
	  G4String mTPCGEMGap_logname = "mTPCGEMGap_log";
	  G4String mTPCGEMGap_physname = "mTPCGEMGap_phys";
	  mTPCGEMGap_solid = 
	    new G4Tubs(mTPCGEMGap_solidname, 
		       mtpcinnerR, mtpcouterR, 
		       fmTPC_gap_GEMGEM/2.0, 
		       0.*deg, 
		       360.*deg);

	  mTPCGEMGap_log = 
	    new G4LogicalVolume(mTPCGEMGap_solid, 
				 GetMaterial("ref4He"),
				mTPCGEMGap_logname);
	  
	  new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, zposgap), 
			    mTPCGEMGap_log, 
			    mTPCGEMGap_physname, 
			    motherlog, 
			    false, 
			    counter, 
			    fChkOvLaps);

	  mTPCGEMGap_log->SetVisAttributes( mtpc_gemgap_visatt );

	}//if not the last gem which has no gap
    }//loop over gems per cell
}//loop over mTPC cells/chambers
  
printf("G4SBSmTPC::BuildmTPCGEMs: Leaving...\n");

}



void G4SBSmTPC::BuildmTPCGasCells(G4LogicalVolume *motherlog, G4double centrecell1, G4double celllength, G4double mtpcinnerR, G4double mtpcouterR)
{//, G4SBSmTPCSD* mtpcSD){//, G4SBSmTPCSD* mtpchvSD){

  printf("G4SBSmTPC::BuildmTPCGasCells: Entering...\n");

  // set up temp sd for readoutHV discs
  // set up SD, for moment only make gas cells sensitive as do not want to record info in gems and readout right now
  G4String mTPCSDname = "SBS/mTPC";
  G4String mTPCcolname = "mTPCHitsCollection";
   
  G4SBSmTPCSD* mTPCSD;
  if( !(mTPCSD = (G4SBSmTPCSD*) fDetCon->fSDman->FindSensitiveDetector(mTPCSDname)) ){ //Make sure SD with this name doesn't already exist
    mTPCSD = new G4SBSmTPCSD( mTPCSDname, mTPCcolname );
    fDetCon->fSDman->AddNewDetector(mTPCSD);
    (fDetCon->SDlist).insert(mTPCSDname);
    fDetCon->SDtype[mTPCSDname] = G4SBS::kmTPC;
  }
   
   
  G4String mTPCHVSDname = "SBS/mTPCHV";
  G4String mTPCHVcolname = "mTPCHVHitsCollection";

  G4SBSCalSD *mTPCHVSD;
  if( !(mTPCHVSD = (G4SBSCalSD*) fDetCon->fSDman->FindSensitiveDetector(mTPCHVSDname)) ){ //Make sure SD with this name doesn't already exist
    mTPCHVSD = new G4SBSCalSD( mTPCHVSDname, mTPCHVcolname );
    fDetCon->fSDman->AddNewDetector(mTPCHVSD);
    (fDetCon->SDlist).insert(mTPCHVSDname);
    fDetCon->SDtype[mTPCHVSDname] = G4SBS::kCAL;
  }
   
  double HVThickness = fmTPC_HV_thick;
   
  // double CellGasLength = celllength - (fmTPC_readout_thick + fmTPC_gap_readoutGEM + (fmTPC_Ngems-1)*fmTPC_gap_GEMGEM
  // 				       + fmTPC_Ngems*(fmTPC_gem_surf1thick + fmTPC_gem_dielecthick + fmTPC_gem_surf2thick));
  double CellGasLength = celllength - (fmTPC_readout_thick + 
				       fmTPC_gap_readoutGEM + 
				       (fmTPC_Ngems-1)*fmTPC_gap_GEMGEM + 
				       fmTPC_Ngems *
				       (fmTPC_gem_surf1thick + fmTPC_gem_dielecthick + fmTPC_gem_surf2thick)
				       ) - HVThickness/2.0;

   
  double zposHVDisc = 0.0;
  //  double mTPCHVDisc_Centre = 0.0;
  G4Tubs* mTPCHVDisc_solid; 
  G4LogicalVolume* mTPCHVDisc_log;
  G4VisAttributes *mtpc_HVDisc_visatt = new G4VisAttributes( G4Colour( 1.0, 0.0, 1.0) );
  mtpc_HVDisc_visatt->SetForceWireframe(true);

  printf(" G4SBSmTPC::BuildmTPCGasCells: test \n");

  G4double GasLayerThick = 5.0*mm; // pad size

  //  cout << "mTPC inner " << mtpcinnerR  << " mTPC outer " << mtpcouterR << endl;

  //G4int NGasLayers = ceil( (mtpcouterR-mtpcinnerR)/GasLayerThick );
  //GasLayerThick = (mtpcouterR-mtpcinnerR)/NGasLayers;
 
  // cout << "NGasLayers " << NGasLayers << " GasLayerThick " << GasLayerThick << endl;
  
  double zposGasCell = 0.0;
  double mTPC_CentreCell = 0.0;
  G4Tubs* mTPCGasCell_solid; 
  G4LogicalVolume* mTPCGasCell_log;
  G4VisAttributes *mtpc_gas_visatt = new G4VisAttributes( G4Colour( 0.0, 0.0, 1.0) );
  mtpc_gas_visatt->SetForceWireframe(true);

  // loop over each cell/chamber of mTPC
  for(int incCell=0; incCell<fmTPC_Ncells; incCell++)
    {
      mTPC_CentreCell = centrecell1 + incCell*celllength;
      if(incCell % 2 == 0)
	{
	  zposGasCell = mTPC_CentreCell + celllength/2.0 - CellGasLength/2.0 - HVThickness/2.0;
	  // HV disc, only have 5
	  zposHVDisc = zposGasCell + CellGasLength/2.0 + HVThickness/2.0;

	  mTPCHVDisc_solid = 
	    new G4Tubs("mTPCHVDisc_solid", 
		       mtpcinnerR, 
		       mtpcouterR, 
		       HVThickness/2.0, 
		       0.*deg, 
		       360.*deg);

	  mTPCHVDisc_log = 
	    new G4LogicalVolume(mTPCHVDisc_solid, 
				 GetMaterial("Au"),
				"mTPCHVDisc_log"); 
	  
	  // if(fmTPCkrypto)mTPCHVDisc_log->SetUserLimits( new G4UserLimits( 0.0, 0.0, 0.0, DBL_MAX, DBL_MAX ) );
	  
	  new G4PVPlacement(0, 
			    G4ThreeVector(0.0, 0.0, zposHVDisc), 
			    mTPCHVDisc_log, 
			    "mTPCHVDisc_phys", 
			    motherlog, 
			    false, 
			    incCell, 
			    fChkOvLaps);
	  mTPCHVDisc_log->SetVisAttributes( mtpc_HVDisc_visatt );
	  // set cell as a sensitive detector
	  mTPCHVDisc_log->SetSensitiveDetector(mTPCHVSD);
	}
      else
	{
	  zposGasCell = mTPC_CentreCell - celllength/2.0 + CellGasLength/2.0 + HVThickness/2.0;
	}
      mTPCGasCell_solid = 
	new G4Tubs("mTPCGasCell_solid", 
		   mtpcinnerR,
		   mtpcouterR, 
		   CellGasLength/2.0, 
		   0.*deg, 
		   360.*deg);
      mTPCGasCell_log = 
	    new G4LogicalVolume(mTPCGasCell_solid, 
				 GetMaterial("mTPCgas"),
				"mTPCGasCell_log");
      mTPCGasCell_log->SetSensitiveDetector(mTPCSD);
      
      new G4PVPlacement(0, 
			G4ThreeVector(0.0, 0.0, zposGasCell), 
			mTPCGasCell_log, 
			"mTPCGasCell_phy",
			motherlog, 
			false, 
			incCell, 
			fChkOvLaps);

      //Shorten the step length (CA)
      if(fmTPCstep != -1*mm)
	{
	  mTPCGasCell_log->SetUserLimits( new G4UserLimits(fmTPCstep, DBL_MAX , DBL_MAX, 0, 0) );
	}

      /*
      for(G4int i_gl = 0; i_gl<NGasLayers; i_gl++)
	{
	  mTPCGasCell_solid = 
	    new G4Tubs("mTPCGasCell_solid", 
		       mtpcinnerR+i_gl*GasLayerThick, 
		       min(mtpcinnerR+(i_gl+1)*GasLayerThick, mtpcouterR), 
		       CellGasLength/2.0, 
		       0.*deg, 
		       360.*deg);

	  mTPCGasCell_log = 
	    new G4LogicalVolume(mTPCGasCell_solid, 
				 GetMaterial("ref4He"),
				"mTPCGasCell_log");    

	  mTPCGasCell_log->SetSensitiveDetector(mTPCSD);
	  G4int layernum = incCell*20+i_gl+1;
	  G4String layername = G4String("mTPCGasCell_phy") + G4String(layernum);

	  new G4PVPlacement(0, 
			    G4ThreeVector(0.0, 0.0, zposGasCell), 
			    mTPCGasCell_log, 
			    layername, 
			    motherlog, 
			    false, 
			    layernum, 
			    fChkOvLaps);

	  //cout << " incCell " << incCell << " i_gl " << i_gl  << " logical mTPC volume number " << layernum << endl;
	  mTPCGasCell_log->SetVisAttributes( mtpc_gas_visatt );
	  // set cell as a sensitive detector

	}
      */
    }

  printf("G4SBSmTPC::BuildmTPCGasCells: Leaving...\n");
}


// EPAF 2024/11/18
// Methods to build MLPC design components.
void G4SBSmTPC::BuildMLPCWires(G4LogicalVolume *motherlog, G4double zPos, G4double angle, G4int WirePlaneNum)
{
  G4Tubs* MLPCWiresGas_solid = 
    new G4Tubs("MLPCWiresGas_solid", 
	       fDetCon->fTargetBuilder->GetTargDiameter()/2,
	       fMLPC_RinSupport, 
	       fMLPC_SupportThickness+fMLPC_WirePlaneThickness*0.5, 
	       0.*deg, 
	       360.*deg );
  
  //mother volume for the wires...
  G4LogicalVolume* MLPCWiresGas_log = 
    new G4LogicalVolume(MLPCWiresGas_solid,
			GetMaterial("mTPCgas"),
			"MLPCWiresGas_log");
  
  new G4PVPlacement(0, 
		    G4ThreeVector(0.0, 0.0, zPos), 
		    MLPCWiresGas_log,
		    "MLPCWiresGas_phys", 
		    motherlog, 
		    false, 
		    WirePlaneNum, 
		    fChkOvLaps);

  G4VisAttributes *tpc_gas_visatt = new G4VisAttributes( G4Colour( 0.0, 0.75, 1.0) );
  tpc_gas_visatt->SetForceWireframe(true);

  G4VisAttributes *wire_visatt = new G4VisAttributes( G4Colour( 0.8, 0.75, 0.0) );
  
  MLPCWiresGas_log->SetVisAttributes( tpc_gas_visatt );
  
  //Split gas in cells for wires... will be tedious
  G4Box* WireCell_box = 
    new G4Box("WireCell_box",
	      fMLPC_InterWireGap*0.5,
	      fMLPC_RinSupport,
	      fMLPC_SupportThickness+fMLPC_WirePlaneThickness*0.5);
   
  G4String mTPCSDname = "SBS/MLPC";
  G4String mTPCcolname = "mTPCHitsCollection";
  
  G4SBSmTPCSD* mTPCSD;
  if( !(mTPCSD = (G4SBSmTPCSD*) fDetCon->fSDman->FindSensitiveDetector(mTPCSDname)) ){ //Make sure SD with this name doesn't already exist
    mTPCSD = new G4SBSmTPCSD( mTPCSDname, mTPCcolname );
    fDetCon->fSDman->AddNewDetector(mTPCSD);
    (fDetCon->SDlist).insert(mTPCSDname);
    fDetCon->SDtype[mTPCSDname] = G4SBS::kmTPC;
  }
  
  for(int i = 0; i<fMLPC_NWiresPerHalfPlane; i++){
    G4double WireDistanceToCenter = fMLPC_InnerWireDistance+fMLPC_InterWireGap*i;

    G4Tubs* Wire_solid =
      new G4Tubs("Wire_solid",
		 0,
		 fMLPC_WirePlaneThickness*0.5,
		 fMLPC_RinSupport*cos( asin( (WireDistanceToCenter+fMLPC_WirePlaneThickness*0.5)/fMLPC_RinSupport ) ),
		 0.*deg,
		 360.*deg);
    
    G4IntersectionSolid* WireCell_solid_0 =
      new G4IntersectionSolid("WireCell_solid_0",
			      MLPCWiresGas_solid,
			      WireCell_box,
			      0,
			      G4ThreeVector(WireDistanceToCenter, 0, 0));
    
    G4RotationMatrix* RotationWire = new G4RotationMatrix();
    RotationWire->rotateX(90.0*deg);
    
    G4SubtractionSolid* WireCell_solid =
      new G4SubtractionSolid("WireCell_solid",
			     WireCell_solid_0,
			     Wire_solid,
			     RotationWire,
			     G4ThreeVector(WireDistanceToCenter, 0, 0));
    
    G4LogicalVolume* WireCell_log = 
      new G4LogicalVolume(WireCell_solid,
			  GetMaterial("mTPCgas"),
			  "WireCell_log");
    
    
    G4LogicalVolume* Wire_log = 
      new G4LogicalVolume(Wire_solid,
			  GetMaterial("Tungsten"),
			  "Wire_log");
    
    WireCell_log->SetVisAttributes( tpc_gas_visatt );
    Wire_log->SetVisAttributes( wire_visatt );

    //G4cout << " Wire distance to center: " << WireDistanceToCenter/mm << " mm " << G4endl; 

    for(int j = 0; j<2; j++){
      G4int WireUniqueID = WirePlaneNum*fMLPC_NWiresPerHalfPlane*2+i*2+j;
      //G4cout << " j "<<  j << " WireUniqueID " << WireUniqueID << " cos angle??? " << cos(angle+j*180.0*deg) << " sin angle??? " << sin(angle+j*180.0*deg) << G4endl; 

      G4RotationMatrix* RotationPlane = new G4RotationMatrix();
      RotationPlane->rotateZ(-angle+j*180.0*deg);

      new G4PVPlacement(RotationPlane,//0
			G4ThreeVector(0, 0, 0),
			WireCell_log,
			"WireCell", 
			MLPCWiresGas_log, 
			false, 
			WireUniqueID, 
			fChkOvLaps
			);
      
      WireCell_log->SetSensitiveDetector(mTPCSD);

      RotationWire = new G4RotationMatrix();
      RotationWire->rotateX(90.0*deg);
      RotationWire->rotateY(angle);
      
      new G4PVPlacement(RotationWire,
			G4ThreeVector(WireDistanceToCenter*cos(angle+j*180.0*deg), WireDistanceToCenter*sin(angle+j*180.0*deg), 0), 
			Wire_log,
			"Wire", 
			MLPCWiresGas_log, 
			false, 
			WireUniqueID, 
			fChkOvLaps
			);
      
    }
  }
}

void G4SBSmTPC::BuildMLPCCathode(G4LogicalVolume *motherlog, G4double zPos)
{
  G4Tubs* MLPCCathode_solid_base = 
    new G4Tubs("MLPCCathode_solid_base", 
	       fMLPC_CathodeRin, 
	       fMLPC_RoutSupport, 
	       fMLPC_CathodeMylarThickness*0.5+fMLPC_CathodeAlThickness, 
	       0.*deg, 
	       360.*deg );
  
  G4LogicalVolume* MLPCCathode_log = 
    new G4LogicalVolume(MLPCCathode_solid_base,
			GetMaterial("Mylar"),
			"MLPCCathode_log");

  G4Tubs* MLPCCathodeAlCoating_solid_base = 
    new G4Tubs("MLPCCathodeAlCoating_solid_base", 
	       fMLPC_CathodeRin, 
	       fMLPC_RoutSupport, 
	       fMLPC_CathodeAlThickness*0.5, 
	       0.*deg, 
	       360.*deg );
  
  G4LogicalVolume* MLPCCathodeAlCoating_log = 
    new G4LogicalVolume(MLPCCathodeAlCoating_solid_base,
			GetMaterial("Aluminum"),
			"MLPCCathodeAlCoating_log");
  
  for(int k = 0; k<2; k++){
    new G4PVPlacement(0, 
		      G4ThreeVector(0.0, 0.0, pow(-1, k) * (fMLPC_CathodeMylarThickness*0.5-fMLPC_CathodeAlThickness*0.5) ), 
		      MLPCCathodeAlCoating_log,
		      "MLPCCathodeAlCoating_phys", 
		      MLPCCathode_log, 
		      false, 
		      0, 
		      fChkOvLaps);
  }
  
  new G4PVPlacement(0, 
		    G4ThreeVector(0.0, 0.0, zPos), 
		    MLPCCathode_log,
		    "MLPCCathode_phys", 
		    motherlog, 
		    false, 
		    0, 
		    fChkOvLaps);

  G4VisAttributes *cathode_visatt = new G4VisAttributes( G4Colour( 0.8, 0.75, 0.0) );
  cathode_visatt->SetForceWireframe(true);
  MLPCCathode_log->SetVisAttributes( cathode_visatt );

  MLPCCathodeAlCoating_log->SetVisAttributes( G4VisAttributes::GetInvisible() );
}

void G4SBSmTPC::BuildMLPCG10Support(G4LogicalVolume *motherlog, G4double zPos)
{
  G4Tubs* MLPCSupport_solid = 
    new G4Tubs("MLPCSupport_solid", 
	       fMLPC_RinSupport, 
	       fMLPC_RoutSupport, 
	       fMLPC_SupportThickness*0.5, 
	       0.*deg, 
	       360.*deg );
  
  G4LogicalVolume* MLPCSupport_log = 
    new G4LogicalVolume(MLPCSupport_solid,
			GetMaterial("NEMAG10"),
			"MLPCSupport_log");
  
  new G4PVPlacement(0, 
		    G4ThreeVector(0.0, 0.0, zPos), 
		    MLPCSupport_log,
		    "MLPCSupport_phys", 
		    motherlog, 
		    false, 
		    0, 
		    fChkOvLaps);

  G4VisAttributes *support_visatt = new G4VisAttributes( G4Colour( 1.0, 1.0, 1.0) );
  support_visatt->SetForceWireframe(true);
  MLPCSupport_log->SetVisAttributes( support_visatt );
 
}


G4SBSmTPC::~G4SBSmTPC()
{
  printf("G4SBSmTPC::~G4SBSmTPC: Finished");
}
