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
      G4int	GEMID;

      G4double xp;
      G4double yp;

  public:
      inline void SetPos(G4ThreeVector v)
      { pos = v;};
      inline void SetDir(G4double x, G4double y)
      { xp = x; yp =y;};
      inline void SetGEMID(G4int i)
      { GEMID = i;};

      inline G4ThreeVector GetPos()
      { return pos;};
      inline G4int GetGEMID()
      { return GEMID;};
      inline G4double GetXp(){ return xp; }
      inline G4double GetYp(){ return yp; }

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
