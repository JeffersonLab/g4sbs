// BeamDiffuser (BD) hit class 
// - Defines data members to store the the energy deposit and track lengths
//   of charged particles in a selected layer
// - Author: D. Flay (JLab) 

#ifndef G4SBS_BEAM_DIFFUSER_HIT_HH
#define G4SBS_BEAM_DIFFUSER_HIT_HH

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

class G4SBSBDHit : public G4VHit
{
  public:
    G4SBSBDHit();
    G4SBSBDHit(const G4SBSBDHit&);
    virtual ~G4SBSBDHit();

    // operators
    const G4SBSBDHit& operator=(const G4SBSBDHit&);
    G4bool operator==(const G4SBSBDHit&) const;

    inline void* operator new(size_t);
    inline void  operator delete(void*);

    // methods from base class
    virtual void Draw() {}
    virtual void Print();

    // setter methods 
    void SetTotalEnergy(G4double E);
    void SetEdep(G4double edep);
    void SetBeta(G4double beta);
    void SetHitTime(G4double time);
    void SetTrackLength(G4double len); 

    void SetTrackID(G4int trackID);
    void SetPlane(G4int i);
    void SetPID(G4int pid);
    void SetMID(G4int mid);
    void SetVerbosity(G4int v); 

    void SetPos(G4ThreeVector v);
    void SetLabPos(G4ThreeVector v);
    void SetMomentum(G4ThreeVector m);

    // getter methods
    G4double GetEdep()        const;
    G4double GetTrackLength() const;
    G4double GetTotalEnergy() const;
    G4double GetMom()         const;
    G4double GetHitTime()     const;
    G4double GetBeta()        const;

    G4int GetTrackID()   const;
    G4int GetPlane()     const;
    G4int GetPID()       const;
    G4int GetMID()       const;
    G4int GetVerbosity() const;

    G4ThreeVector GetPos()      const;
    G4ThreeVector GetLabPos()   const;
    G4ThreeVector GetMomentum() const;

  private:
    G4double fEdep;        // Energy deposit in the sensitive volume
    G4double fTrackLength; // Track length in the sensitive volume
    G4double fEtot;        // Total energy (at pre-step)
    G4double fBeta;        // Particle speed 
    G4double fHitTime;     // Time of hit  

    G4int fTrackID;        // Track number 
    G4int fPlane;          // Plane number
    G4int fPID;            // Particle type 
    G4int fMID;            // Material type
    G4int fVerbosity;  

    G4ThreeVector fPos;    // Local hit coordinate 
    G4ThreeVector fLabPos; // Global hit coordinate 
    G4ThreeVector fMom;    // Momentum 

};

using G4SBSBDHitsCollection = G4THitsCollection<G4SBSBDHit>;
extern G4ThreadLocal G4Allocator<G4SBSBDHit>* G4SBSBDHitAllocator;

//______________________________________________________________________________
inline void* G4SBSBDHit::operator new(size_t)
{
  if (!G4SBSBDHitAllocator) {
    G4SBSBDHitAllocator = new G4Allocator<G4SBSBDHit>;
  }
  void *hit;
  hit = (void *) G4SBSBDHitAllocator->MallocSingle();
  return hit;
}
//______________________________________________________________________________
inline void G4SBSBDHit::operator delete(void *hit)
{
  if (!G4SBSBDHitAllocator) {
    G4SBSBDHitAllocator = new G4Allocator<G4SBSBDHit>;
  }
  G4SBSBDHitAllocator->FreeSingle((G4SBSBDHit*) hit);
}
//______________________________________________________________________________
inline void G4SBSBDHit::SetVerbosity(G4int v){
   fVerbosity = v;
}
//______________________________________________________________________________
inline void G4SBSBDHit::SetTrackID(G4int trackID){
   fTrackID = trackID;
}
//______________________________________________________________________________
inline void G4SBSBDHit::SetEdep(G4double edep){
   fEdep = edep;
}
//______________________________________________________________________________
inline void G4SBSBDHit::SetTotalEnergy(G4double E){
   fEtot = E;
}
//______________________________________________________________________________
inline void G4SBSBDHit::SetBeta(G4double beta){
   fBeta = beta;
}
//______________________________________________________________________________
inline void G4SBSBDHit::SetHitTime(G4double time){
   fHitTime = time;
}
//______________________________________________________________________________
inline void G4SBSBDHit::SetTrackLength(G4double len){
   fTrackLength = len;
}
//______________________________________________________________________________
inline void G4SBSBDHit::SetPos(G4ThreeVector v){
   fPos = v;
}
//______________________________________________________________________________
inline void G4SBSBDHit::SetLabPos(G4ThreeVector v){
   fLabPos = v;
}
//______________________________________________________________________________
inline void G4SBSBDHit::SetMomentum(G4ThreeVector m){
   fMom = m;
}
//______________________________________________________________________________
inline void G4SBSBDHit::SetPlane(G4int i){
   fPlane = i;
}
//______________________________________________________________________________
inline void G4SBSBDHit::SetPID(G4int pid){
   fPID = pid;
}
//______________________________________________________________________________
inline void G4SBSBDHit::SetMID(G4int mid){
   fMID = mid;
}
//______________________________________________________________________________
inline G4double G4SBSBDHit::GetEdep() const {
  return fEdep;
}
//______________________________________________________________________________
inline G4double G4SBSBDHit::GetTrackLength() const {
  return fTrackLength;
}
//______________________________________________________________________________
inline G4double G4SBSBDHit::GetTotalEnergy() const{
   return fEtot;
}
//______________________________________________________________________________
inline G4double G4SBSBDHit::GetBeta() const{
   return fBeta;
}
//______________________________________________________________________________
inline G4double G4SBSBDHit::GetHitTime() const{
   return fHitTime;
}
//______________________________________________________________________________
inline G4ThreeVector G4SBSBDHit::GetPos() const{
   return fPos;
}
//______________________________________________________________________________
inline G4ThreeVector G4SBSBDHit::GetLabPos() const{
   return fLabPos;
}
//______________________________________________________________________________
inline G4ThreeVector G4SBSBDHit::GetMomentum() const{
   return fMom;
}
//______________________________________________________________________________
inline G4double G4SBSBDHit::GetMom() const{
   double x      = fMom.x();
   double y      = fMom.y();
   double z      = fMom.z();
   double sum_sq = x*x + y*y + z*z;
   return sqrt(sum_sq);
}
//______________________________________________________________________________
inline G4int G4SBSBDHit::GetVerbosity() const{
   return fVerbosity;
}
//______________________________________________________________________________
inline G4int G4SBSBDHit::GetPlane() const{
   return fPlane;
}
//______________________________________________________________________________
inline G4int G4SBSBDHit::GetPID() const{
   return fPID;
}
//______________________________________________________________________________
inline G4int G4SBSBDHit::GetMID() const{
   return fMID;
}
//______________________________________________________________________________
inline G4int G4SBSBDHit::GetTrackID() const{
   return fTrackID;
}

#endif
