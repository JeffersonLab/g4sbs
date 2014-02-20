#include "G4SBSGEMHit.hh"

#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4AttValue.hh"
#include "G4AttDef.hh"
#include "G4AttCheck.hh"

G4Allocator<G4SBSGEMHit> G4SBSGEMHitAllocator;

G4SBSGEMHit::G4SBSGEMHit()
{pos = G4ThreeVector(); GEMID = -1; xp = -1e9; yp = -1e9;}

G4SBSGEMHit::~G4SBSGEMHit()
{;}

G4SBSGEMHit::G4SBSGEMHit(const G4SBSGEMHit &right)
  : G4VHit()
{
  pos = right.pos;
  xp = right.xp;
  yp = right.yp;
  GEMID = right.GEMID;

  vert = right.vert;
  trid = right.trid;
  mid = right.mid;
  pid = right.pid;
  p = right.p;
  edep = right.edep;
}

const G4SBSGEMHit& G4SBSGEMHit::operator=(const G4SBSGEMHit &right)
{
  pos = right.pos;
  GEMID= right.GEMID;
  xp = right.xp;
  yp = right.yp;
  p = right.p;
  edep = right.edep;

  vert = right.vert;
  trid = right.trid;
  mid = right.mid;
  pid = right.pid;

  return *this;
}

G4int G4SBSGEMHit::operator==(const G4SBSGEMHit &right) const
{
  return (this==&right) ? 1 : 0;
}

void G4SBSGEMHit::Draw()
{
}

void G4SBSGEMHit::Print()
{
}


