#include "G4SBSCalHit.hh"

#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"

G4Allocator<G4SBSCalHit> G4SBSCalHitAllocator;

G4SBSCalHit::G4SBSCalHit()
{pos = G4ThreeVector();
    energy = 0;
    mid = -1;
    pid = -1e9;
}

G4SBSCalHit::~G4SBSCalHit()
{;}

G4SBSCalHit::G4SBSCalHit(const G4SBSCalHit &right)
  : G4VHit()
{
  pos = right.pos;
  vertex = right.vertex;
  energy = right.energy;
  pid = right.pid;
  mid = right.mid;
  trid = right.trid;
  cell = right.cell;
  row = right.row;
  col = right.col;
}

const G4SBSCalHit& G4SBSCalHit::operator=(const G4SBSCalHit &right)
{
  pos = right.pos;
  vertex = right.vertex;
  energy = right.energy;
  pid = right.pid;
  mid = right.mid;
  trid = right.trid;
  cell = right.cell;
  row = right.row;
  col = right.col;
  return *this;
}

G4int G4SBSCalHit::operator==(const G4SBSCalHit &right) const
{
  return (this==&right) ? 1 : 0;
}

void G4SBSCalHit::Draw()
{
}

void G4SBSCalHit::Print()
{
}


