// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class ArraySD
// Save the interesting parameter from the sensitive volumes
// ie the active regions of the detectors
// 03/05/14 JRMA

#ifndef ArraySD_h
#define ArraySD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class G4HCofThisEvent;
class G4Step;
class RootIO;

#include "ArrayHit.hh"

class ArraySD : public G4VSensitiveDetector
{
public:
  
  ArraySD(G4String,G4int,G4double*);
  ~ArraySD();
  
  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step* astep,G4TouchableHistory* ROHist);
  void EndOfEvent(G4HCofThisEvent*);
  void clear();
  void DrawAll();
  void PrintAll();
  G4double GetEffectiveEnergyDeposit(const G4Step*);
  void SetEthr(G4double ethr){ fEthr = ethr; }
  G4double GetEThr(){ return fEthr; }
  void SetTmax(G4double tm){ if(tm > 0.0) fTmax = tm; }
  G4double GetTmax(){ return fTmax; }
  G4int* GetNhits(){ return &fNhits; }
  G4int* GetHitID(){ return fHitID; }
  G4int* GetHits(){ return fHits; }
private:
  ArrayHitsCollection*  fCollection;
  RootIO* fRio;
  G4int fHCID;
  G4int fNelements;
  G4int * fHitID;
  G4int * fHits;
  G4int fNhits;
  G4double fEthr;
  G4double fTmax;
  G4int pdgCode;
  G4int trID;
  G4int parID;
  G4int stepID;
  G4int gTrID;
  G4int pTrID;
  G4int nTrID;
  G4int Tr1PDG;
  G4double fZst;
  G4double fRHe1;
  G4double fRHe2;
  G4int fZPixG;
  G4int fPhiPixG;
  G4double fdZPixG;
  G4double fdPhiPixG;
  G4double* fPix;
  G4int fIsStore;
};

#endif






