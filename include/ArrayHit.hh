// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class ArrayHit
// Define the interesting parameter from the sensitive volumes
// ie the active regions of the detectors
// 06/05/14 JRMA

#ifndef ArrayHit_h
#define ArrayHit_h 1
#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "TMath.h"

enum { ENSt=1000 };            // max number of subhits

struct HitStep{
  G4int trID;
  G4int parID;
  G4int stepID;
  G4double dE;
  G4LorentzVector p4; // start momentum
  G4ThreeVector pre;  // prestep pos
  G4ThreeVector post; // poststep pos
};

class ArrayHit : public G4VHit
{
public:
  
  ArrayHit();
  ~ArrayHit();
  ArrayHit(const ArrayHit&);
  const ArrayHit& operator=(const ArrayHit&);
  //  int operator==(const ArrayHit&) const;
  inline void* operator new(size_t);
  inline void  operator delete(void*);
  virtual void Draw();
  void Print();

protected:
  
  G4double fEdep;     // Energy deposited in detector
  G4double fEdepB;    // Birks weighted energy deposited in detector
  G4ThreeVector fPos; // Position of the hit (in what frame?)
  G4int fID;          // ID of detector hit
  //G4int fPDG;       // particle type
  G4double fTime;     // global time of hit
  G4double fTimeE;    // hit time over threshold
  G4ThreeVector fPol; // spin alignment
  G4int fNst;         // no. times element hit in event (sub-hits)
  HitStep fStepList[ENSt];   // List of steps
  //G4double fTimeS[ENSH];    // stored times of sub-hits
  //G4double fEdepBS[ENSH];   // stored energy deposits in sub-hits
  //G4int fISH[ENSH];         // for sorting times

public:
  
  inline void SetEnergy(G4double de) {fEdep = de;};
  inline void SetEnergyB(G4double de) {fEdepB = de;};
  inline void AddEnergy(G4double de) {fEdep += de;};
  inline void AddEnergyB(G4double de) {fEdepB += de;};
  inline void SetPos(G4ThreeVector pos) {fPos=pos;};
  inline void SetID(G4int i) {fID = i;};
  //inline void SetPDG(G4int i) {fPDG = i;};
  inline void SetTime(G4double t) {fTime = t;};
  inline void SetTimeE(G4double t) {if( !fTimeE )fTimeE = t;};
  inline void SetDedx(G4double dedx){fTimeE = dedx;};
  inline G4bool AddStep(G4int trid, G4int parid, G4int stepid, G4double de,
		      G4LorentzVector* p4, G4ThreeVector* pre, 
		      G4ThreeVector* post)
  {
    if(fNst >= ENSt-1) return false;
    HitStep* st = fStepList + fNst;
    st->trID = trid;
    st->parID = parid;
    st->stepID = stepid;
    st->dE = de;
    st->p4.setX(p4->x());
    st->p4.setY(p4->y());
    st->p4.setZ(p4->z());
    st->p4.setE(p4->e());
    st->pre.setX(pre->x());
    st->pre.setY(pre->y());
    st->pre.setZ(pre->z());
    st->post.setX(post->x());
    st->post.setY(post->y());
    st->post.setZ(post->z());
    fNst++;
    return true;
  } 
  //inline void SetPol(G4ThreeVector p) { fPol = p; p.set(0,0,0);}
  inline void ResetSH(){ fNst = 0; }
  inline G4double GetEdep()     { return fEdep; };
  inline G4double GetEdepB()     { return fEdepB; };
  inline G4ThreeVector GetPos() { return fPos; };
  inline G4int    GetID(){ return fID; };
  //inline G4int    GetPDG(){ return fPDG; };
  inline G4double    GetTime()   { return fTime; };
  inline G4double    GetTimeE()   { return fTimeE; };
  inline HitStep* GetHitStep(){ return fStepList; }
  inline HitStep* GetHitStep(G4int i){ return fStepList+i; }
  inline G4int GetNSt(){ return fNst; }
  //inline G4ThreeVector GetPol() { return fPol;};
  //inline void AddSH( G4double, G4double );
  //inline void SortSH( G4double );
  inline void Clear(){
    fEdep = fEdepB = fTime = fTimeE = 0.0;
    fID = -1;
    fNst = 0;
    fPos.set(0,0,0);
  }
};


typedef G4THitsCollection<ArrayHit> ArrayHitsCollection;

extern G4Allocator<ArrayHit> ArrayHitAllocator;


inline void* ArrayHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) ArrayHitAllocator.MallocSingle();
  return aHit;
}


inline void ArrayHit::operator delete(void* aHit)
{
  ArrayHitAllocator.FreeSingle((ArrayHit*) aHit);
}

/*
inline void ArrayHit::AddSH(G4double e, G4double t)
{
  if( fNst >= ENSH ) return;
  fTimeS[fNst] = t;
  fEdepBS[fNst] = e;
  fNst++;
}
inline void ArrayHit::SortSH(G4double thresh)
{
  TMath::Sort(fNst,fTimeS,fISH,kFALSE);
  G4double esum = 0.0;
  G4int j;
  for(G4int i=0; i<fNst; i++){
    j = fISH[i];
    esum += fEdepBS[j];
    if(esum >= thresh) break;
  }
  fTimeE = fTimeS[j];
}
*/
#endif










