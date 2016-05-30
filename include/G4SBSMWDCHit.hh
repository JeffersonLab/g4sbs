#ifndef G4SBSMWDCHit_h
#define G4SBSMWDCHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4SBSMWDCHit : public G4VHit
{
public:

  G4SBSMWDCHit();
  ~G4SBSMWDCHit();
  G4SBSMWDCHit(const G4SBSMWDCHit &right);
  const G4SBSMWDCHit& operator=(const G4SBSMWDCHit &right);
  G4int operator==(const G4SBSMWDCHit &right) const;


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
  G4ThreeVector globalpos; //global position, for drawing
  G4ThreeVector vert;
  G4ThreeVector polarization;
  G4int	MWDCID;
  //G4int TrackerID;

  G4int	trid;
  G4int	mid;
  G4int	pid;

  G4double xp;
  G4double yp;
  G4double p, edep;
  G4double hittime;
  G4double beta; //v/c, for timing:

public:
  inline void SetPos(G4ThreeVector v)
  { pos = v;};
  inline void SetGlobalPos( G4ThreeVector v )
  { globalpos = v; }
  inline void SetVertex(G4ThreeVector v)
  { vert = v;};
  inline void SetPolarization(G4ThreeVector v)
  { polarization = v;}
  inline void SetDir(G4double x, G4double y)
  { xp = x; yp =y;};
  inline void SetMWDCID(G4int i)
  { MWDCID = i;};
  // inline void SetTrackerID(G4int i) 
  // { TrackerID= i; }
  inline void SetTrID(G4int i)
  { trid = i;};
  inline void SetMID(G4int i)
  { mid = i;};
  inline void SetPID(G4int i)
  { pid = i;};
  inline void SetMom(G4double x)
  { p = x;};
  inline void SetEdep(G4double x)
  { edep = x;};
  inline void SetHittime( G4double t )
  { hittime = t; }
  inline void SetBeta( G4double b )
  { beta = b; }

  inline G4ThreeVector GetPos()
  { return pos;};
  inline G4ThreeVector GetGlobalPos()
  { return globalpos; }
  inline G4ThreeVector GetVertex()
  { return vert;};
  inline G4ThreeVector GetPolarization()
  { return polarization;}
  inline G4int GetMWDCID()
  { return MWDCID;};
  // inline G4int GetTrackerID()
  // { return TrackerID; };
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

typedef G4THitsCollection<G4SBSMWDCHit> G4SBSMWDCHitsCollection;

extern G4Allocator<G4SBSMWDCHit> G4SBSMWDCHitAllocator;

inline void* G4SBSMWDCHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) G4SBSMWDCHitAllocator.MallocSingle();
  return aHit;
}

inline void G4SBSMWDCHit::operator delete(void *aHit)
{
  G4SBSMWDCHitAllocator.FreeSingle((G4SBSMWDCHit*) aHit);
}

#endif
