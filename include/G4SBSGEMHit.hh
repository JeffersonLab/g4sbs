#ifndef G4SBSGEMHit_h
#define G4SBSGEMHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4SBSGEMHit : public G4VHit
{
  public:

      G4SBSGEMHit();
      ~G4SBSGEMHit();
      G4SBSGEMHit(const G4SBSGEMHit &right);
      const G4SBSGEMHit& operator=(const G4SBSGEMHit &right);
      G4int operator==(const G4SBSGEMHit &right) const;


      inline void *operator new(size_t);
      inline void operator delete(void *aHit);
      void *operator new(size_t,void*p){return p;}
#ifndef G4NOT_ISO_DELETES
      void operator delete(void *,void*){}
#endif

      void Draw();
      void Print();

  private:
      G4ThreeVector pos;
      G4ThreeVector vert;
      G4int	GEMID;
  G4int TrackerID;
      G4int	trid;
      G4int	mid;
      G4int	pid;

      G4double xp;
      G4double yp;
      G4double p, edep;
  G4double hittime;
  G4double beta; //v/c, for timing:
  public:
      inline void SetVertex(G4ThreeVector v)
      { vert= v;};
      inline void SetPos(G4ThreeVector v)
      { pos = v;};
      inline void SetDir(G4double x, G4double y)
      { xp = x; yp =y;};
      inline void SetGEMID(G4int i)
      { GEMID = i;};
  inline void SetTrackerID(G4int i) 
  { TrackerID= i; }
      inline void SetTrID(G4int i)
      { trid= i;};
      inline void SetMID(G4int i)
      { mid= i;};
      inline void SetPID(G4int i)
      { pid= i;};
      inline void SetMom(G4double x)
      { p=x;};
      inline void SetEdep(G4double x)
      {edep=x;};
  inline void SetHittime( G4double t )
  { hittime = t; }
  inline void SetBeta( G4double b )
  { beta = b; }

      inline G4ThreeVector GetPos()
      { return pos;};
      inline G4ThreeVector GetVertex()
      { return vert;};
      inline G4int GetGEMID()
      { return GEMID;};
  inline G4int GetTrackerID()
  { return TrackerID; };
      inline G4int GetPID()
      { return pid;};
      inline G4int GetMID()
      { return mid;};
      inline G4int GetTrID()
      { return trid;};
      inline G4double GetXp(){ return xp; }
      inline G4double GetYp(){ return yp; }
      inline G4double GetMom(){ return p; }
      inline G4double GetEdep(){ return edep; }
  inline G4double GetHittime(){ return hittime; }
  inline G4double GetBeta(){ return beta; }
};

typedef G4THitsCollection<G4SBSGEMHit> G4SBSGEMHitsCollection;

extern G4Allocator<G4SBSGEMHit> G4SBSGEMHitAllocator;

inline void* G4SBSGEMHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) G4SBSGEMHitAllocator.MallocSingle();
  return aHit;
}

inline void G4SBSGEMHit::operator delete(void *aHit)
{
  G4SBSGEMHitAllocator.FreeSingle((G4SBSGEMHit*) aHit);
}

#endif
