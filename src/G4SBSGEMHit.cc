#include "G4SBSGEMHit.hh"

#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4Circle.hh"
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
  vert = right.vert;
  polarization = right.polarization;
  GEMID = right.GEMID;
  //TrackerID = right.TrackerID;
  trid = right.trid;
  mid = right.mid;
  pid = right.pid;
  xp = right.xp;
  yp = right.yp;
  p = right.p;
  edep = right.edep;
  hittime = right.hittime;
  beta = right.beta;
}

const G4SBSGEMHit& G4SBSGEMHit::operator=(const G4SBSGEMHit &right)
{
  pos = right.pos;
  vert = right.vert;
  polarization = right.polarization;
  GEMID = right.GEMID;
  //TrackerID = right.TrackerID;
  trid = right.trid;
  mid = right.mid;
  pid = right.pid;
  xp = right.xp;
  yp = right.yp;
  p = right.p;
  edep = right.edep;
  hittime = right.hittime;
  beta = right.beta;

  return *this;
}

G4int G4SBSGEMHit::operator==(const G4SBSGEMHit &right) const
{
  return (this==&right) ? 1 : 0;
}

void G4SBSGEMHit::Draw()
{
  G4VVisManager *pVVisManager = G4VVisManager::GetConcreteInstance();
  if( pVVisManager ){
    G4Circle circle( pos );
    circle.SetScreenSize(5);
    circle.SetFillStyle( G4Circle::filled );
    G4VisAttributes attribs(G4Colour(0.0,1.0,0.0) );
    circle.SetVisAttributes( attribs );
    pVVisManager->Draw(circle);
  }
}

void G4SBSGEMHit::Print()
{
}


