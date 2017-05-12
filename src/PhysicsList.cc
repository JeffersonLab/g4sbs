// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class PhysicsList
// Control of physics models used to describe particle interactions in matter
// 20/05/13 JRMA adapted from SBS equivalent, under construction

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "G4DecayPhysics.hh"
#include "G4EmStandardPhysics_option2.hh"
#include "G4EmStandardPhysics_option1.hh"
#include "G4EmStandardPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "EmPhysics.hh"
#include "G4HadronDElasticPhysics.hh"
//#include "G4HadronQElasticPhysics.hh"
#include "G4HadronHElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"
//#include "G4QStoppingPhysics.hh"
//#include "G4LHEPStoppingPhysics.hh"
#include "G4IonBinaryCascadePhysics.hh"
#include "G4IonPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmProcessOptions.hh"
#include "G4BertiniElectroNuclearBuilder.hh"

//#include "HadronPhysicsLHEP.hh"
//#include "HadronPhysicsLHEP_EMV.hh"
#include "G4HadronInelasticQBBC.hh"
//#include "G4HadronPhysicsQGSC_BERT.hh"
//#include "G4HadronPhysicsQGSP.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "G4HadronPhysicsQGSP_BIC.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"

#include "G4IonPhysics.hh"
#include "G4LossTableManager.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4SystemOfUnits.hh"

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  G4LossTableManager::Instance();
  defaultCutValue = 1.*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;
  dump            = false;
  fIsSpinPrecN    = false;
  pMessenger = new PhysicsListMessenger(this);
  // Particles
  particleList = new G4DecayPhysics("decays");
  // EM physics
  emPhysicsList = new G4EmStandardPhysics();
}


PhysicsList::~PhysicsList()
{
  delete pMessenger;
  delete particleList;
  delete emPhysicsList;
  for(size_t i=0; i<hadronPhys.size(); i++) {
    delete hadronPhys[i];
  }
}

void PhysicsList::ConstructParticle()
{
  particleList->ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  //  thePLHelper->theTransportationProcess->EnableUseMagneticMoment()
  emPhysicsList->ConstructProcess();
  particleList->ConstructProcess();
  for(size_t i=0; i<hadronPhys.size(); i++) {
    hadronPhys[i]->ConstructProcess();
  }
}

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>0)
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  if (name == "emstandard") {

    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics();
    
  } else if (name == "emTest") {
    
    delete emPhysicsList;
    emPhysicsList = new EmPhysics();
    
  } else if (name == "emTest1") {

    //delete emPhysicsList;
    //emPhysicsList = new G4EmExtraPhysics();
    hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
    //G4BertiniElectroNuclearBuilder* xxx = new G4BertiniElectroNuclearBuilder();
    //emPhysicsList = new G4BertiniElectroNuclearBuilder();
    //xxx->Build();

  } else if (name == "emTest2") {

    delete emPhysicsList;
    emPhysicsList = new G4EmExtraPhysics();
    //emPhysicsList = new G4BertiniElectroNuclearBuilder();

  } else  if (name == "emstandard_opt2") {

    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option2();

  } else if (name == "emstandard_opt1") {

    delete emPhysicsList;
    emPhysicsList = new G4EmStandardPhysics_option1();
    
  } else if (name == "FTFP_EMV") {

    AddPhysicsList("emstandard_opt1");
    AddPhysicsList("FTFP");

  } else if (name == "LHEP") {

    SetBuilderList2();
    //hadronPhys.push_back( new HadronPhysicsLHEP("hadron"));
    dump = true;

  } else if (name == "LHEP_EMV") {

    AddPhysicsList("emstandard_opt1");
    SetBuilderList3();
    //hadronPhys.push_back( new HadronPhysicsLHEP_EMV("hadron"));
    dump = true;

  } else if (name == "QBBC") {

    SetBuilderList0();
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,false,false,false,true));

  } else if (name == "QBBC_DEL") {

    SetBuilderList5();
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,false,false,false,true));

  } else if (name == "QBBC_HEL") {

    SetBuilderList6();
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,false,false,false,true));

  } else if (name == "QBBC_HP") {

    SetBuilderList0();
    hadronPhys.push_back( new G4HadronInelasticQBBC("QBBC",verboseLevel,
						    false,false,false,true,true));
  } else if (name == "QGSP_BERT_HP") {

    SetBuilderList1();
    hadronPhys.push_back( new G4HadronPhysicsQGSP_BERT_HP("hadron",true));

  } else if (name == "QGSP_BIC") {

    SetBuilderList0();
    hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC("hadron",true));
    dump = true;

  } else if (name == "NoHadronic") {

    SetNoHadList();
    dump = true;

  } else if (name == "QGSP_BIC_HP") {

    SetBuilderList0();
    hadronPhys.push_back( new G4HadronPhysicsQGSP_BIC_HP("hadron",true));
    dump = true;

  } else {

    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
  }
}

void PhysicsList::SetBuilderList0()
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  hadronPhys.push_back( new G4HadronElasticPhysics(verboseLevel));
  //hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4IonBinaryCascadePhysics("ionBIC"));
  hadronPhys.push_back( new G4NeutronTrackingCut("Neutron tracking cut"));
}

void PhysicsList::SetNoHadList()
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  //hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4NeutronTrackingCut("Neutron tracking cut"));
}

void PhysicsList::SetBuilderList1()
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  hadronPhys.push_back( new G4HadronElasticPhysics(verboseLevel));
  //hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4IonPhysics("ion"));
  hadronPhys.push_back( new G4NeutronTrackingCut("Neutron tracking cut"));
}

void PhysicsList::SetBuilderList2()
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  hadronPhys.push_back( new G4HadronElasticPhysics(verboseLevel));
  hadronPhys.push_back( new G4IonPhysics("ion"));
}

void PhysicsList::SetBuilderList3()
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  hadronPhys.push_back( new G4HadronElasticPhysics(verboseLevel));
  //hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4IonPhysics("ion"));
}

void PhysicsList::SetBuilderList4()
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  //hadronPhys.push_back( new G4HadronQElasticPhysics("elastic",verboseLevel));
  //hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4IonPhysics("ion"));
  hadronPhys.push_back( new G4NeutronTrackingCut("Neutron tracking cut"));
}

void PhysicsList::SetBuilderList5()
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  //hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4IonBinaryCascadePhysics("ionBIC"));
  hadronPhys.push_back( new G4NeutronTrackingCut("Neutron tracking cut"));
}

void PhysicsList::SetBuilderList6()
{
  hadronPhys.push_back( new G4EmExtraPhysics("extra EM"));
  hadronPhys.push_back( new G4HadronHElasticPhysics(verboseLevel));
  //hadronPhys.push_back( new G4QStoppingPhysics("stopping",verboseLevel));
  hadronPhys.push_back( new G4IonBinaryCascadePhysics("ionBIC"));
  hadronPhys.push_back( new G4NeutronTrackingCut("Neutron tracking cut"));
}

void PhysicsList::SetCuts()
{

  if (verboseLevel >0){
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");

  if (verboseLevel>0) DumpCutValuesTable();
}

void PhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

void PhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

void PhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

//----------------------------------------------------------------------------
void PhysicsList::List()
{
  G4cout << "### PhysicsLists available: FTFC FTFP FTFP_BERT FTFP_EMV LHEP LHEP_BERT LHEP_EMV "
	 << G4endl;
  G4cout << "                            LHEP_PRECO_HP QBBC QBBC_DEL QBBC_HEL QBBC_HP QGSC "
	 << G4endl; 
  G4cout << "                            QGSC_BERT QGSC_EFLOW QGSC_EMV QGSP QGSP_BERT QGSP_BER_EMV "
	 << G4endl; 
  G4cout << "                            QGSP_BERT_HP QGSP_BIC QGSP_BIC_HP QGSP_DIF " 
	 << G4endl; 
  G4cout << "                            QGSP_EMV QGSP_EMX QGSP_NQE QGSP_QEL  "
	 << G4endl; 
  G4cout << "                            QGSP_BIC_POL QGSP_BIC_POLN NoHadronic"
	 << G4endl; 
}

