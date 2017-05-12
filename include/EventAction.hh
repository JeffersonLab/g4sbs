// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class EventAction
// Data saving procedures event by event
// 05/05/14 JRMA
//
#ifndef EventAction_h
#define EventAction_h 1
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "DetectorConstruction.hh"

class RunAction;
class EventActionMessenger;
class PrimaryGeneratorAction;
class RootIO;

class EventAction : public G4UserEventAction
{
public:
  EventAction(RunAction*);
  ~EventAction();
public:
  void BeginOfEventAction(const G4Event*);
  void EndOfEventAction(const G4Event*);
  void SetDrawFlag   (G4String val)  {fdrawFlag = val;}
  void SetPrintModulo(G4int    val)  {fprintModulo = val;}
  void SetCollID(G4int val){fCollID=val;}
  void SetIsInteractive(G4int is){fIsInteractive=is;}
  void SetHitDrawOpt(G4String val){fHitDrawOpt=val;}
  void SetOutFileName(G4String name){fOutFileName=name;}
  G4int PrepareOutput();
  void CloseOutput();
  void SetDet(DetectorConstruction* det){ fDet = det; }
  RootIO* GetRootIO(){ return fRTPCOut; }
private:
  RunAction*  frunAct;
  RootIO* fRTPCOut;
  DetectorConstruction* fDet;
  G4int fIsInteractive;        
  G4String  fdrawFlag;
  G4int     fprintModulo;
  G4int     fDrawMode;
  G4String fHitDrawOpt;
  EventActionMessenger*  feventMessenger;
  G4int fCollID;
  G4String fOutFileName;
  PrimaryGeneratorAction* fPGA;
};

#endif

    
