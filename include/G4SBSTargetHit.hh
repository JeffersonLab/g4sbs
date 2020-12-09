#ifndef G4SBS_TARGET_HIT_HH
#define G4SBS_TARGET_HIT_HH

// a class for keeping track of hits in a target system 

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4Threading.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>
#include <fstream> 

class G4SBSTargetHit : public G4VHit
{
  public:
    G4SBSTargetHit();
    G4SBSTargetHit(const G4SBSTargetHit&);
    virtual ~G4SBSTargetHit();

    // operators
    const G4SBSTargetHit& operator=(const G4SBSTargetHit&);
    G4bool operator==(const G4SBSTargetHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();
    
    void PrintToFile(bool v=true)      { fPrintToCSV = v; }
    void PrintToCSV();  

    // setter methods 
    void Add(G4double de, G4double dl);
    void SetTotalEnergy(G4double E);
    void SetEdep(G4double edep);
    void SetBeta(G4double beta);
    void SetHitTime(G4double time);
    void SetTrackLength(G4double len);
    void SetMomentumMag(G4double pmag); 

    void SetTrackID(G4int trackID);
    void SetPID(G4int pid);
    void SetMID(G4int mid);

    void SetPos(G4ThreeVector v);
    void SetLabPos(G4ThreeVector v);
    void SetMomentum(G4ThreeVector m);

    // getter methods
    G4double GetEdep()        const;
    G4double GetTrackLength() const;
    G4double GetTotalEnergy() const;
    G4double GetMomentumMag() const;
    G4double GetHitTime()     const;
    G4double GetBeta()        const;

    G4int GetTrackID() const;
    G4int GetPID()     const;
    G4int GetMID()     const;

    G4ThreeVector GetPos()      const;
    G4ThreeVector GetLabPos()   const;
    G4ThreeVector GetMomentum() const;

  private:
    G4double fEdep;        // Energy deposit in the sensitive volume
    G4double fTrackLength; // Track length in the sensitive volume
    G4double fEtot;        // Total energy (at pre-step)
    G4double fBeta;        // Particle speed 
    G4double fHitTime;     // Time of hit 
    G4double fPmag;        // momentum magnitude  

    G4int fTrackID;        // Track number 
    G4int fPID;            // Particle type 
    G4int fMID;            // Material type 

    G4ThreeVector fPos;    // Local hit coordinate 
    G4ThreeVector fLabPos; // Global hit coordinate 
    G4ThreeVector fMom;    // Momentum

    bool fPrintToCSV;      // print data to csv? 
    unsigned long fCntr;   // a counter  

};

using G4SBSTargetHitsCollection = G4THitsCollection<G4SBSTargetHit>;
extern G4ThreadLocal G4Allocator<G4SBSTargetHit>* G4SBSTargetHitAllocator;

//______________________________________________________________________________
inline void* G4SBSTargetHit::operator new(size_t)
{
  if (!G4SBSTargetHitAllocator) {
    G4SBSTargetHitAllocator = new G4Allocator<G4SBSTargetHit>;
  }
  void *hit;
  hit = (void *) G4SBSTargetHitAllocator->MallocSingle();
  return hit;
}
//______________________________________________________________________________
inline void G4SBSTargetHit::operator delete(void *hit)
{
  if (!G4SBSTargetHitAllocator) {
    G4SBSTargetHitAllocator = new G4Allocator<G4SBSTargetHit>;
  }
  G4SBSTargetHitAllocator->FreeSingle((G4SBSTargetHit*) hit);
}
//______________________________________________________________________________
inline void G4SBSTargetHit::SetTrackID(G4int trackID){
   fTrackID = trackID;
}
//______________________________________________________________________________
inline void G4SBSTargetHit::Add(G4double de, G4double dl) {
  fEdep        += de;
  fTrackLength += dl;
}
//______________________________________________________________________________
inline void G4SBSTargetHit::SetTrackLength(G4double len){
   fTrackLength = len;
}
//______________________________________________________________________________
inline void G4SBSTargetHit::SetEdep(G4double edep){
   fEdep = edep;
}
//______________________________________________________________________________
inline void G4SBSTargetHit::SetTotalEnergy(G4double E){
   fEtot = E;
}
//______________________________________________________________________________
inline void G4SBSTargetHit::SetBeta(G4double beta){
   fBeta = beta;
}
//______________________________________________________________________________
inline void G4SBSTargetHit::SetHitTime(G4double time){
   fHitTime = time;
}
//______________________________________________________________________________
inline void G4SBSTargetHit::SetMomentumMag(G4double pmag){
   fPmag = pmag; 
}
//______________________________________________________________________________
inline void G4SBSTargetHit::SetPos(G4ThreeVector v){
   fPos = v;
}
//______________________________________________________________________________
inline void G4SBSTargetHit::SetLabPos(G4ThreeVector v){
   fLabPos = v;
}
//______________________________________________________________________________
inline void G4SBSTargetHit::SetMomentum(G4ThreeVector m){
   fMom = m;
}
//______________________________________________________________________________
inline void G4SBSTargetHit::SetPID(G4int pid){
   fPID = pid;
}
//______________________________________________________________________________
inline void G4SBSTargetHit::SetMID(G4int mid){
   fMID = mid;
}
//______________________________________________________________________________
inline G4double G4SBSTargetHit::GetEdep() const {
  return fEdep;
}
//______________________________________________________________________________
inline G4double G4SBSTargetHit::GetTrackLength() const {
  return fTrackLength;
}
//______________________________________________________________________________
inline G4double G4SBSTargetHit::GetTotalEnergy() const{
   return fEtot;
}
//______________________________________________________________________________
inline G4double G4SBSTargetHit::GetBeta() const{
   return fBeta;
}
//______________________________________________________________________________
inline G4double G4SBSTargetHit::GetHitTime() const{
   return fHitTime;
}
//______________________________________________________________________________
inline G4ThreeVector G4SBSTargetHit::GetPos() const{
   return fPos;
}
//______________________________________________________________________________
inline G4ThreeVector G4SBSTargetHit::GetLabPos() const{
   return fLabPos;
}
//______________________________________________________________________________
inline G4ThreeVector G4SBSTargetHit::GetMomentum() const{
   return fMom;
}
//______________________________________________________________________________
inline G4double G4SBSTargetHit::GetMomentumMag() const{
   return fPmag;
}
//______________________________________________________________________________
inline G4int G4SBSTargetHit::GetPID() const{
   return fPID;
}
//______________________________________________________________________________
inline G4int G4SBSTargetHit::GetMID() const{
   return fMID;
}
//______________________________________________________________________________
inline G4int G4SBSTargetHit::GetTrackID() const{
   return fTrackID;
}

#endif 
