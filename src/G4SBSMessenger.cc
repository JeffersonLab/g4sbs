#include "TBuffer.h"
#include "TString.h"
#include "TMatrixTBase.h"
#include "THashTable.h"
#include "G4SBSMessenger.hh"
#include "G4SBSRun.hh"

#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"

#include "G4SBSDetectorConstruction.hh"
#include "G4SBSIO.hh"
#include "G4SBSEventGen.hh"
#include "G4SBSEventAction.hh"

#include "G4MagneticField.hh"
#include "G4SBSGlobalField.hh"
#include "G4SBSMagneticField.hh"
#include "G4SBSToscaField.hh"
#include "G4SBSBigBiteField.hh"
#include "G4SBSPrimaryGeneratorAction.hh"
#include "G4SBSPhysicsList.hh"
#include "G4OpticalPhysics.hh"

#include "G4SBSBeamlineBuilder.hh"
#include "G4SBSTargetBuilder.hh"
#include "G4SBSEArmBuilder.hh"
#include "G4SBSHArmBuilder.hh"

#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"

#include "G4UImanager.hh"
#include "G4RunManager.hh"
#ifdef G4SBS_USE_GDML
#include "G4GDMLParser.hh"
#endif
#include "G4VPhysicalVolume.hh"

using namespace CLHEP;

G4SBSMessenger::G4SBSMessenger(){
    fExpType = kNeutronExp;

    runCmd = new G4UIcmdWithAnInteger("/g4sbs/run",this);
    runCmd->SetGuidance("Run simulation with x events");
    runCmd->SetParameterName("nevt", false);

    gemconfigCmd = new G4UIcmdWithAnInteger("/g4sbs/gemconfig",this);
    gemconfigCmd->SetGuidance("Change between GEM configurations");
    gemconfigCmd->SetParameterName("gemconfig", false);

    fileCmd = new G4UIcmdWithAString("/g4sbs/filename",this);
    fileCmd->SetGuidance("Output filename");
    fileCmd->SetParameterName("filename", false);

    sigfileCmd = new G4UIcmdWithAString("/g4sbs/sigmafile",this);
    sigfileCmd->SetGuidance("GEM Sigma filename");
    sigfileCmd->SetParameterName("sigmafile", false);

    tgtCmd = new G4UIcmdWithAString("/g4sbs/target",this);
    tgtCmd->SetGuidance("Target type from LH2, LD2, H2, 3He");
    tgtCmd->SetParameterName("targtype", false);

    kineCmd = new G4UIcmdWithAString("/g4sbs/kine",this);
    kineCmd->SetGuidance("Kinematic type");
    kineCmd->SetParameterName("kinetype", false);

    expCmd = new G4UIcmdWithAString("/g4sbs/exp",this);
    expCmd->SetGuidance("Experiment type");
    expCmd->SetParameterName("exptype", false);

    HadrCmd = new G4UIcmdWithAString("/g4sbs/hadron",this);
    HadrCmd->SetGuidance("Hadron type h for SIDIS N(e,e'h)X reaction");
    HadrCmd->SetParameterName("hadrontype", false );

    bigfieldCmd = new G4UIcmdWithAnInteger("/g4sbs/48d48field", this);
    bigfieldCmd->SetGuidance("48d48 magnet field");
    bigfieldCmd->SetParameterName("48d48field", false);

    bbfieldCmd = new G4UIcmdWithAnInteger("/g4sbs/bbfield", this);
    bbfieldCmd->SetGuidance("Bigbite field");
    bbfieldCmd->SetParameterName("bbfield", false);

    tosfieldCmd = new G4UIcmdWithAString("/g4sbs/tosfield", this);
    tosfieldCmd->SetGuidance("Tosca field");
    tosfieldCmd->SetParameterName("tosfield", false);

    geantinoCmd = new G4UIcmdWithABool("/g4sbs/shootgeantino", this);
    geantinoCmd->SetGuidance("Shoot a geantino instead of e-");
    geantinoCmd->SetParameterName("shootgeantino", false);

    invertCmd = new G4UIcmdWithABool("/g4sbs/invertfield", this);
    invertCmd->SetGuidance("invert field polarity");
    invertCmd->SetParameterName("invert", false);

    totalabsCmd = new G4UIcmdWithABool("/g4sbs/totalabs", this);
    totalabsCmd->SetGuidance("Magnet materials are total absorbers");
    totalabsCmd->SetParameterName("totalabs", false);

    tgtLenCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/targlen",this);
    tgtLenCmd->SetGuidance("Target length");
    tgtLenCmd->SetParameterName("targlen", false);

    tgtDenCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/targden",this);
    tgtDenCmd->SetGuidance("Target density");
    tgtDenCmd->SetParameterName("targden", false);

    tgtPresCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/targpres",this);
    tgtPresCmd->SetGuidance("Gaseous Target pressure");
    tgtPresCmd->SetParameterName("targpres", false);

    beamcurCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamcur",this);
    beamcurCmd->SetGuidance("Beam current");
    beamcurCmd->SetParameterName("beamcur", false);

    runtimeCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/runtime",this);
    runtimeCmd->SetGuidance("Run time");
    runtimeCmd->SetParameterName("runtime", false);

    rasterxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/rasterx",this);
    rasterxCmd->SetGuidance("Raster x size");
    rasterxCmd->SetParameterName("size", false);

    rasteryCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/rastery",this);
    rasteryCmd->SetGuidance("Raster y size");
    rasteryCmd->SetParameterName("size", false);

    beamECmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/beamE",this);
    beamECmd->SetGuidance("Beam Energy");
    beamECmd->SetParameterName("energy", false);

    bbangCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/bbang",this);
    bbangCmd->SetGuidance("BigBite angle");
    bbangCmd->SetParameterName("angle", false);

    bbdistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/bbdist",this);
    bbdistCmd->SetGuidance("BigBite distance");
    bbdistCmd->SetParameterName("dist", false);

    hcalangCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hcalang",this);
    hcalangCmd->SetGuidance("HCAL angle");
    hcalangCmd->SetParameterName("angle", false);

    hcaldistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hcaldist",this);
    hcaldistCmd->SetGuidance("HCAL distance");
    hcaldistCmd->SetParameterName("dist", false);

    hmagdistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/48D48dist",this);
    hmagdistCmd->SetGuidance("48D48 distance");
    hmagdistCmd->SetParameterName("dist", false);

    gemresCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/gemres",this);
    gemresCmd->SetGuidance("GEM resolution");
    gemresCmd->SetParameterName("dist", false);

    // Detector position commands

    cerDisCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/cerdist",this);
    cerDisCmd->SetGuidance("Cerenkov distance from front GEM");
    cerDisCmd->SetParameterName("dist", false);

    cerDepCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/cerdepth",this);
    cerDepCmd->SetGuidance("Cerenkov gas depth");
    cerDepCmd->SetParameterName("dist", false);

    gemSepCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/gemsep",this);
    gemSepCmd->SetGuidance("GEM separation from front to back set");
    gemSepCmd->SetParameterName("dist", false);

    bbCalDistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/bbcaldist",this);
    bbCalDistCmd->SetGuidance("BigBite caloriter distance from front GEM");
    bbCalDistCmd->SetParameterName("dist", false);

    thminCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/thmin",this);
    thminCmd->SetGuidance("Minimum electron generation polar angle");
    thminCmd->SetParameterName("angle", false);

    thmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/thmax",this);
    thmaxCmd->SetGuidance("Maximum electron generation polar angle");
    thmaxCmd->SetParameterName("angle", false);

    phminCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/phmin",this);
    phminCmd->SetGuidance("Minimum electron generation azimuthal angle");
    phminCmd->SetParameterName("angle", false);

    phmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/phmax",this);
    phmaxCmd->SetGuidance("Maximum electron generation azimuthal angle");
    phmaxCmd->SetParameterName("angle", false);

    HthminCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hthmin",this);
    HthminCmd->SetGuidance("Minimum hadron generation polar angle");
    HthminCmd->SetParameterName("htheta",false);

    HthmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hthmax",this);
    HthmaxCmd->SetGuidance("Maximum hadron generation polar angle");
    HthmaxCmd->SetParameterName("htheta",false);

    HphminCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hphmin",this);
    HphminCmd->SetGuidance("Minimum hadron generation azimuthal angle");
    HphminCmd->SetParameterName("htheta",false);

    HphmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/hphmax",this);
    HphmaxCmd->SetGuidance("Maximum hadron generation azimuthal angle");
    HphmaxCmd->SetParameterName("htheta",false);
    
    EhminCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/ehmin",this);
    EhminCmd->SetGuidance("Minimum hadron generation energy");
    EhminCmd->SetParameterName("ehmin",false);
    
    EhmaxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/ehmax",this);
    EhmaxCmd->SetGuidance("Maximum hadron generation energy");
    EhmaxCmd->SetParameterName("ehmax",false);

    EeminCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/eemin",this);
    EeminCmd->SetGuidance("Minimum electron generation energy");
    EeminCmd->SetParameterName("eemin",false);
    
    EemaxCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/eemax",this);
    EemaxCmd->SetGuidance("Maximum electron generation energy");
    EemaxCmd->SetParameterName("eemax",false);

    RICHdistCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/richdist",this);
    RICHdistCmd->SetGuidance("RICH distance from target");
    RICHdistCmd->SetParameterName("dist",false);

    SBSMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/g4sbs/sbsmagfield",this);
    SBSMagFieldCmd->SetGuidance("SBS magnetic field setting");
    SBSMagFieldCmd->SetParameterName("sbsbfield",false);

    SBSFieldClampOptionCmd = new G4UIcmdWithAnInteger("/g4sbs/sbsclampopt",this);
    SBSFieldClampOptionCmd->SetGuidance("SBS field clamp configuration: 1=BigBite(default), 2=GEp");
    SBSFieldClampOptionCmd->SetParameterName("sbsclampoption",false);

    SBSLeadOptionCmd = new G4UIcmdWithAnInteger("/g4sbs/uselead",this);
    SBSLeadOptionCmd->SetGuidance("SBS lead configuration: 0= nope 1=yaes");
    SBSLeadOptionCmd->SetParameterName("uselead",false);

    //Optical physics toggle commands:
    UseCerenkovCmd = new G4UIcmdWithABool( "/g4sbs/usecerenkov",this );
    UseCerenkovCmd->SetGuidance("Activate Cherenkov radiation (default true)");
    UseCerenkovCmd->SetParameterName("useckov",true);
    UseCerenkovCmd->SetDefaultValue(true);

    UseScintCmd = new G4UIcmdWithABool( "/g4sbs/usescint",this );
    UseScintCmd->SetGuidance("Activate Scintillation (default false)");
    UseScintCmd->SetParameterName("usescint",true);
    UseScintCmd->SetDefaultValue(false);

    UseOpRayleighCmd = new G4UIcmdWithABool( "/g4sbs/userayleigh",this );
    UseOpRayleighCmd->SetGuidance("Activate Rayleigh scattering of optical photons (default true)");
    UseOpRayleighCmd->SetParameterName("useoprayleigh",true);
    UseOpRayleighCmd->SetDefaultValue(true);
    
    UseOpAbsorbCmd = new G4UIcmdWithABool( "/g4sbs/useabsorb",this );
    UseOpAbsorbCmd->SetGuidance("Activate absorption of optical photons (default true)");
    UseOpAbsorbCmd->SetParameterName("useabsorb",true);
    UseOpAbsorbCmd->SetDefaultValue(true);

    UseOpBdryCmd = new G4UIcmdWithABool( "/g4sbs/useboundary",this );
    UseOpBdryCmd->SetGuidance("Activate optical boundary processes (default true)");
    UseOpBdryCmd->SetParameterName("usebdry",true);
    UseOpBdryCmd->SetDefaultValue(true);

    UseOpWLSCmd = new G4UIcmdWithABool( "/g4sbs/usewls",this );
    UseOpWLSCmd->SetGuidance("Activate optical wavelength shifting process (default false)");
    UseOpWLSCmd->SetParameterName("usewls",true);
    UseOpWLSCmd->SetDefaultValue(false);

    UseOpMieHGCmd = new G4UIcmdWithABool( "/g4sbs/usemie",this );
    UseOpMieHGCmd->SetGuidance("Activate Mie scattering (default false)");
    UseOpMieHGCmd->SetParameterName("usemie",true);
    UseOpMieHGCmd->SetDefaultValue(false);
}

G4SBSMessenger::~G4SBSMessenger(){
}


void G4SBSMessenger::SetNewValue(G4UIcommand* cmd, G4String newValue){
    char cmdstr[255];

    if( cmd == runCmd ){
	
	G4VPhysicalVolume* pWorld;

	G4int nevt = runCmd->GetNewIntValue(newValue);
	fevgen->SetNevents(nevt);
	
	//Clean out and rebuild the detector geometry from scratch: 

	G4SolidStore::GetInstance()->Clean();
	G4LogicalVolumeStore::GetInstance()->Clean();
	G4PhysicalVolumeStore::GetInstance()->Clean();
	
	G4RunManager::GetRunManager()->DefineWorldVolume(pWorld = fdetcon->ConstructAll());
	G4RunManager::GetRunManager()->GeometryHasBeenModified();

	fevact->SDlist = fdetcon->SDlist;

	// Clobber old gdml if it exists and write out the
	// present geometry
	// Save geometry to GDML file
#ifdef G4SBS_USE_GDML
	G4GDMLParser parser;
	unlink("g4sbs.gdml");
	parser.Write("g4sbs.gdml", pWorld);
#endif
	// Run the simulation
	G4UImanager * UImanager = G4UImanager::GetUIpointer();
	sprintf(cmdstr, "/run/beamOn %d", nevt);
	UImanager->ApplyCommand(cmdstr);
    }

    if( cmd == fileCmd ){
	fIO->SetFilename(newValue.data());
    }

    if( cmd == sigfileCmd ){
	fevact->LoadSigmas(newValue.data());
    }

    if( cmd == gemconfigCmd ){
	int gemconfval = gemconfigCmd->GetNewIntValue(newValue);
	fdetcon->fEArmBuilder->SetGEMConfig(gemconfval);
    }

    if( cmd == kineCmd ){
	bool validcmd = false;
	if( newValue.compareTo("elastic") == 0 ){
	    fevgen->SetKine(kElastic);
	    validcmd = true;
	}
	if( newValue.compareTo("inelastic") == 0 ){
	    fevgen->SetKine(kInelastic);
	    validcmd = true;
	}
	if( newValue.compareTo("flat") == 0 ){
	    fevgen->SetKine(kFlat);
	    validcmd = true;
	}
	if( newValue.compareTo("dis") == 0 ){
	    fevgen->SetKine(kDIS);
	    validcmd = true;
	}
	if( newValue.compareTo("beam") == 0 ){
	    fevgen->SetKine(kBeam);
	    validcmd = true;
	}
	if( newValue.compareTo("sidis") == 0 ){
	  fevgen->SetKine( kSIDIS );
	  validcmd = true;
	}
	if( !validcmd ){
	    fprintf(stderr, "%s: %s line %d - Error: kinematic type %s not valid\n", __PRETTY_FUNCTION__, __FILE__, __LINE__, newValue.data());
	    exit(1);
	} else {
	    G4SBSRun::GetRun()->GetData()->SetGenName(newValue.data());
	}
    }

    if( cmd == expCmd ){
	bool validcmd = false;
	if( newValue.compareTo("gep") == 0 ){
	    fExpType = kGEp;
	    validcmd = true;
	}
	if( newValue.compareTo("gmn") == 0 ){
	    fExpType = kNeutronExp;
	    validcmd = true;
	}
	if( newValue.compareTo("gen") == 0 ){
	    fExpType = kNeutronExp;
	    validcmd = true;
	}
	if( newValue.compareTo("a1n") == 0 ){
	    fExpType = kNeutronExp;
	    validcmd = true;
	}
	//AJP: Add SIDIS as a valid experiment type:
	if( newValue.compareTo("sidis") == 0 ){
	  fExpType = kSIDISExp;
	  validcmd = true;
	}

	if( validcmd ){
	    fdetcon->SetExpType( fExpType );
	    G4SBSRun::GetRun()->GetData()->SetExpType(newValue.data());
	}

	if( !validcmd ){
	    fprintf(stderr, "%s: %s line %d - Error: kinematic type %s not valid\n", __PRETTY_FUNCTION__, __FILE__, __LINE__, newValue.data());
	    exit(1);
	}
    }


    if( cmd == HadrCmd ){
      bool validcmd = false;
      if( newValue.compareTo("pi+") == 0 ){
	fevgen->SetHadronType( kPiPlus );
	validcmd = true;
      }
      if( newValue.compareTo("pi-") == 0 ){
	fevgen->SetHadronType( kPiMinus );
	validcmd = true;
      }
      if( newValue.compareTo("pi0") == 0 ){
	fevgen->SetHadronType( kPi0 );
	validcmd = true;
      }
      if( newValue.compareTo("K+") == 0 ){
	fevgen->SetHadronType( kKPlus );
	validcmd = true;
      }
      if( newValue.compareTo("K-") == 0 ){
	fevgen->SetHadronType( kKMinus );
	validcmd = true; 
      }
      if( newValue.compareTo("p") == 0 ){
	fevgen->SetHadronType( kP );
	validcmd = true;
      } 
      if( newValue.compareTo("pbar") == 0 ){
	fevgen->SetHadronType( kPbar );
	validcmd = true;
      }

      if( !validcmd ){
	  fprintf(stderr, "%s: %s line %d - Error: Hadron type %s not valid\n", __PRETTY_FUNCTION__, __FILE__, __LINE__, newValue.data());
	  exit(1);
      }


    }

    if( cmd == tgtCmd ){
	bool validcmd = false;
	if( newValue.compareTo("LH2") == 0 ){
	    fevgen->SetTarget(kLH2);
	    fdetcon->SetTarget(kLH2);

	    G4double den = (0.071*g/cm3)*Avogadro/(1.008*g/mole);
	    fevgen->SetTargDen(den);
	    fdetcon->fTargetBuilder->SetTargDen(den);
	    validcmd = true;
	}
	if( newValue.compareTo("H2") == 0 ){
	    fevgen->SetTarget(kH2);
	    fdetcon->SetTarget(kH2);

	    G4double den = 10.0*atmosphere/(296.0*kelvin*k_Boltzmann);
	    fevgen->SetTargDen(den);
	    fdetcon->fTargetBuilder->SetTargDen(den);
	    validcmd = true;
	}
	if( newValue.compareTo("LD2") == 0 ){
	    fevgen->SetTarget(kLD2);
	    fdetcon->SetTarget(kLD2);

	    G4double den = (162.4*kg/m3)*Avogadro/(2.014*g/mole);
	    fevgen->SetTargDen(den);
	    fdetcon->fTargetBuilder->SetTargDen(den);
	    validcmd = true;
	}
	if( newValue.compareTo("3He") == 0 ){
	    fevgen->SetTarget(k3He);
	    fdetcon->SetTarget(k3He);

	    G4double den = 10.0*atmosphere/(296.0*kelvin*k_Boltzmann);
	    fevgen->SetTargDen(den);
	    fdetcon->fTargetBuilder->SetTargDen(den);
	    validcmd = true;

	}
	if( newValue.compareTo("Neutron") == 0 ){
	    fevgen->SetTarget(kNeutTarg);
	    fdetcon->SetTarget(kNeutTarg);

	    G4double den = 10.0*atmosphere/(296.0*kelvin*k_Boltzmann);
	    fevgen->SetTargDen(den);
	    fdetcon->fTargetBuilder->SetTargDen(den);
	    validcmd = true;
	}

	if( !validcmd ){
	    fprintf(stderr, "%s: %s line %d - Error: target type %s not valid\n", __PRETTY_FUNCTION__, __FILE__, __LINE__, newValue.data());
	    exit(1);
	}

    }

    if( cmd == bigfieldCmd ){
	G4int n = bigfieldCmd->GetNewIntValue(newValue);
	fdetcon->Set48D48Field(n);
    }

    if( cmd == bbfieldCmd ){
	G4int n = bbfieldCmd->GetNewIntValue(newValue);
	fdetcon->SetBigBiteField(n);
    }

    if( cmd == tosfieldCmd ){
	fdetcon->AddToscaField(newValue.data());
    }

    if( cmd == geantinoCmd ){
	G4bool b = geantinoCmd->GetNewBoolValue(newValue);
	fprigen->SetUseGeantino(b);
	fdetcon->GetGlobalField()->SetInvertField(b);
    }

    if( cmd == invertCmd ){
	G4bool b = invertCmd->GetNewBoolValue(newValue);
	fdetcon->GetGlobalField()->SetInvertField(b);
    }

    if( cmd == totalabsCmd ){
	G4bool b = totalabsCmd->GetNewBoolValue(newValue);
	fdetcon->SetTotalAbs(b);
    }

    if( cmd == tgtLenCmd ){
	G4double len = tgtLenCmd->GetNewDoubleValue(newValue);
	fevgen->SetTargLen(len);
	fdetcon->fTargetBuilder->SetTargLen(len);
    }

    if( cmd == tgtDenCmd ){
	G4double den = tgtDenCmd->GetNewDoubleValue(newValue);
	fevgen->SetTargDen(den);
	fdetcon->fTargetBuilder->SetTargDen(den);
    }
    if( cmd == tgtPresCmd ){
	G4double pre = tgtPresCmd->GetNewDoubleValue(newValue);
	G4double den = pre/(296.0*kelvin*k_Boltzmann);
	fevgen->SetTargDen(den);
	fdetcon->fTargetBuilder->SetTargDen(den);
    }

    if( cmd == beamcurCmd ){
	G4double v = beamcurCmd->GetNewDoubleValue(newValue);
	printf("Setting beam current to %f uA\n", v/microampere);
	fevgen->SetBeamCur(v);
    }
    if( cmd == runtimeCmd ){
	G4double v = runtimeCmd->GetNewDoubleValue(newValue);
	fevgen->SetRunTime(v);
    }

    if( cmd == rasterxCmd ){
	G4double v = rasterxCmd->GetNewDoubleValue(newValue);
	fevgen->SetRasterX(v);
    }

    if( cmd == rasteryCmd ){
	G4double v = rasteryCmd->GetNewDoubleValue(newValue);
	fevgen->SetRasterY(v);
    }

    if( cmd == beamECmd ){
	G4double v = beamECmd->GetNewDoubleValue(newValue);
	fevgen->SetBeamE(v);
	fIO->SetBeamE(v);

	G4SBSRun::GetRun()->GetData()->SetBeamE(v/GeV);
    }

    if( cmd == bbangCmd ){
	G4double v = bbangCmd->GetNewDoubleValue(newValue);
	printf("Setting BB ang to %f deg\n", v/deg);
	fdetcon->SetBBAng(v);
	fIO->SetBigBiteTheta(v);
    }

    if( cmd == bbdistCmd ){
	G4double v = bbdistCmd->GetNewDoubleValue(newValue);
	fdetcon->SetBBDist(v);
	fIO->SetBigBiteDist(v);
    }

    if( cmd == hcalangCmd ){
	G4double v = hcalangCmd->GetNewDoubleValue(newValue);
	fdetcon->Set48D48Ang(v);
	fIO->SetHcalTheta(v);
    }

    if( cmd == hcaldistCmd ){
	G4double v = hcaldistCmd->GetNewDoubleValue(newValue);
	fdetcon->fHArmBuilder->SetHCALDist(v);
	fevgen->SetHCALDist(v);
	fIO->SetHcalDist(v);
    }

    if( cmd == hmagdistCmd ){
	G4double v = hmagdistCmd->GetNewDoubleValue(newValue);
	fdetcon->Set48D48Dist(v);
    }

    if( cmd == cerDepCmd ){
	G4double v = cerDepCmd->GetNewDoubleValue(newValue);
	fdetcon->fEArmBuilder->SetCerDepth(v);
    }

    if( cmd == cerDisCmd ){
	G4double v = cerDisCmd->GetNewDoubleValue(newValue);
	fdetcon->fEArmBuilder->SetCerDist(v);
    }

    if( cmd == gemSepCmd ){
	G4double v = gemSepCmd->GetNewDoubleValue(newValue);
	fdetcon->fEArmBuilder->SetGEMSep(v);
    }

    if( cmd == bbCalDistCmd ){
	G4double v = bbCalDistCmd->GetNewDoubleValue(newValue);
	fdetcon->fEArmBuilder->SetBBCalDist(v);
    }

    if( cmd == thminCmd ){
	G4double v = thminCmd->GetNewDoubleValue(newValue);
	fevgen->SetThMin(v);
    }
    if( cmd == thmaxCmd ){
	G4double v = thmaxCmd->GetNewDoubleValue(newValue);
	fevgen->SetThMax(v);
    }
    if( cmd == phminCmd ){
	G4double v = phminCmd->GetNewDoubleValue(newValue);
	fevgen->SetPhMin(v);
    }
    if( cmd == phmaxCmd ){
	G4double v = phmaxCmd->GetNewDoubleValue(newValue);
	fevgen->SetPhMax(v);
    }
    if( cmd == HthminCmd ){
      G4double v = HthminCmd->GetNewDoubleValue(newValue);
      fevgen->SetThMin_had(v);
    }
    if( cmd == HthmaxCmd ){
      G4double v = HthmaxCmd->GetNewDoubleValue(newValue);
      fevgen->SetThMax_had(v);
    }

    if( cmd == HphminCmd ){
      G4double v = HphminCmd->GetNewDoubleValue(newValue);
      fevgen->SetPhMin_had(v);
    }
    if( cmd == HphmaxCmd ){
      G4double v = HphmaxCmd->GetNewDoubleValue(newValue);
      fevgen->SetPhMax_had(v);
    }

    if( cmd == EhminCmd ){
      G4double v = EhminCmd->GetNewDoubleValue(newValue);
      fevgen->SetEhadMin(v);
    }
    if( cmd == EhmaxCmd ){
      G4double v = EhmaxCmd->GetNewDoubleValue(newValue);
      fevgen->SetEhadMax(v);
    }
    if( cmd == EeminCmd ){
      G4double v = EeminCmd->GetNewDoubleValue(newValue);
      fevgen->SetEeMin(v);
    }
    if( cmd == EemaxCmd ){
      G4double v = EemaxCmd->GetNewDoubleValue(newValue);
      fevgen->SetEeMax(v);
    }

    if( cmd == gemresCmd ){
	G4double v = gemresCmd->GetNewDoubleValue(newValue);
	fevact->SetGEMRes(v);
    }

    if( cmd == RICHdistCmd ){
      G4double v = RICHdistCmd->GetNewDoubleValue(newValue);
      fdetcon->fHArmBuilder->SetRICHdist(v);
    }

    if( cmd == SBSMagFieldCmd ){
      G4double v = SBSMagFieldCmd->GetNewDoubleValue(newValue);
      fdetcon->SetUniformMagneticField48D48( v );
    }

    if( cmd == SBSFieldClampOptionCmd ){
      G4int i = SBSFieldClampOptionCmd->GetNewIntValue(newValue);
      fdetcon->fHArmBuilder->SetFieldClampConfig48D48( i );
    }

    if( cmd == SBSLeadOptionCmd ){
      G4int i = SBSLeadOptionCmd->GetNewIntValue(newValue);
      fdetcon->fLeadOption = i;
    }

    if( cmd == UseCerenkovCmd ){
      G4bool isactive = UseCerenkovCmd->GetNewBoolValue(newValue);
      fphyslist->SetOpticalPhysicsProcessActive( kCerenkov, isactive );
    }

    if( cmd == UseScintCmd ){
      G4bool isactive = UseScintCmd->GetNewBoolValue(newValue);
      fphyslist->SetOpticalPhysicsProcessActive( kScintillation, isactive );
    }
    
    if( cmd == UseOpRayleighCmd ){
      G4bool isactive = UseOpRayleighCmd->GetNewBoolValue(newValue);
      fphyslist->SetOpticalPhysicsProcessActive( kRayleigh, isactive );
    }

    if( cmd == UseOpAbsorbCmd ){
      G4bool isactive = UseOpAbsorbCmd->GetNewBoolValue(newValue);
      fphyslist->SetOpticalPhysicsProcessActive( kAbsorption, isactive );
    }

    if( cmd == UseOpBdryCmd ){
      G4bool isactive = UseOpBdryCmd->GetNewBoolValue(newValue);
      fphyslist->SetOpticalPhysicsProcessActive( kBoundary, isactive );
    }

    if( cmd == UseOpWLSCmd ){
      G4bool isactive = UseOpWLSCmd->GetNewBoolValue(newValue);
      fphyslist->SetOpticalPhysicsProcessActive( kWLS, isactive );
    }

    if( cmd == UseOpMieHGCmd ){
      G4bool isactive = UseOpMieHGCmd->GetNewBoolValue(newValue);
      fphyslist->SetOpticalPhysicsProcessActive( kMieHG, isactive );
    }
}
