#ifndef G4SBSCalHit_h
#define G4SBSCalHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"

class G4SBSCalHit : public G4VHit
{
public:

  G4SBSCalHit();
  ~G4SBSCalHit();
  G4SBSCalHit(const G4SBSCalHit &right);
  const G4SBSCalHit& operator=(const G4SBSCalHit &right);
  G4int operator==(const G4SBSCalHit &right) const;


  inline void *operator new(size_t);
  inline void operator delete(void *aHit);
  void *operator new(size_t,void*p){return p;}
#ifndef G4NOT_ISO_DELETES
  void operator delete(void *,void*){}
#endif

  void Draw();
  void Print();

private:
  G4ThreeVector pos; //Local coordinate of the hit
  G4ThreeVector vertex; //vertex coordinate of particle producing hit.
  G4ThreeVector labpos; //Global coordinate of the hit.
  G4ThreeVector mom; //momentum of the particle before the step

  double hittime; //time of the hit (since start of event)
  double edep; //energy deposition of the hit
  double energy; //Initial energy of the particle prior to the hit.
  double Lstep; //spatial length of the step (magnitude)

  G4int cell, row, col, plane, wire; //Channel information 
  G4ThreeVector CellCoords; //"local" coordinates of center of cell in which hit occurred
  G4ThreeVector GlobalCellCoords; //"global" coordinates of center of cell in which hit occurred

  G4int pid, mid, trid; //pid of particle causing hit, mother track id and track id.

public:
  inline void SetPos(G4ThreeVector v)
  { pos = v;};

  inline void SetVertex(G4ThreeVector v)
  { vertex = v;};

  inline void SetLabPos(G4ThreeVector v)
  { labpos = v;};

  inline void SetMomentum(G4ThreeVector v)
  { mom = v; };

  inline void SetTime(G4double t)
  { hittime = t; };

  inline void SetEdep(G4double e)
  { edep = e;};

  inline void SetEnergy(G4double E)
  { energy = E; }

  inline void SetLstep(G4double L)
  { Lstep = L; }

  inline void SetCell(G4int c)
  { cell = c; }
  
  inline void SetRow(G4int r)
  { row = r; }
  
  inline void SetCol(G4int c)
  { col = c; }

  inline void SetPlane(G4int p)
  { plane = p; }

  inline void SetWire(G4int w)
  { wire = w; }
  
  inline void SetCellCoords( G4ThreeVector x )
  { CellCoords = x; }

  inline void SetGlobalCellCoords( G4ThreeVector x )
  { GlobalCellCoords = x; }
  

  inline void SetPID(G4int p)
  { pid = p;};

  inline void SetMID(G4int mother)
  { mid = mother;};

  inline void SetTrID(G4int t)
  { trid = t;};

  
  ////////GET methods://////////
  inline G4ThreeVector GetPos()
  { return pos;};

  inline G4ThreeVector GetVertex()
  { return vertex;};

  inline G4ThreeVector GetLabPos()
  { return labpos;};

  inline G4ThreeVector GetMomentum()
  { return mom; };

  inline double GetTime()
  { return hittime;};

  inline G4double GetEdep()
  { return edep;};

  inline G4double GetEnergy()
  { return energy; };

  inline G4double GetLstep()
  { return Lstep; }
  
  inline G4int GetCell() { return cell; }
  inline G4int GetRow() { return row; }
  inline G4int GetCol() { return col; }
  inline G4int GetPlane() { return plane; }
  inline G4int GetWire() { return wire; }
  
  inline G4ThreeVector GetCellCoords() { return CellCoords; }
  inline G4ThreeVector GetGlobalCellCoords() { return GlobalCellCoords; }

  inline G4int GetPID()
  { return pid;};
  
  inline G4int GetMID()
  { return mid;};

  inline G4int GetTrID()
  { return trid;};
  
};

typedef G4THitsCollection<G4SBSCalHit> G4SBSCalHitsCollection;

extern G4Allocator<G4SBSCalHit> G4SBSCalHitAllocator;

inline void* G4SBSCalHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *) G4SBSCalHitAllocator.MallocSingle();
  return aHit;
}

inline void G4SBSCalHit::operator delete(void *aHit)
{
  G4SBSCalHitAllocator.FreeSingle((G4SBSCalHit*) aHit);
}

#endif
