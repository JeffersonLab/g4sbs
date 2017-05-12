// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class PhysicsList
// Control of physics models used to describe particle interactions in matter
// 20/05/13 JRMA adapted from SBS equivalent, under construction

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class PhysicsListMessenger;

class PhysicsList: public G4VModularPhysicsList
{
public:

  PhysicsList();
  virtual ~PhysicsList();
  void ConstructParticle();
  void SetCuts();
  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);
  void AddPhysicsList(const G4String& name);
  void ConstructProcess();
  void List();
private:
  void SetNoHadList();
  void SetBuilderList0();
  void SetBuilderList1();
  void SetBuilderList2();
  void SetBuilderList3();
  void SetBuilderList4();
  void SetBuilderList5();
  void SetBuilderList6();
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  G4VPhysicsConstructor*  emPhysicsList;
  G4VPhysicsConstructor*  particleList;
  std::vector<G4VPhysicsConstructor*>  hadronPhys;
    
  PhysicsListMessenger* pMessenger;
  G4bool dump;
  G4bool fIsSpinPrecN;
  G4double fmuP;
  G4double fmuN;
public:
  void SetSpinPrecN( G4double muP, G4double muN )
  { fIsSpinPrecN = true; fmuP = muP; fmuN = muN; }
};

#endif

