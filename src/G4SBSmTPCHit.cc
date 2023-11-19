#include "G4SBSmTPCHit.hh"

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

G4Allocator<G4SBSmTPCHit> *G4SBSmTPCHitAllocator;

G4SBSmTPCHit::G4SBSmTPCHit() : G4VHit() //Is invocation of the G4VHit() constructor needed?
{;}

G4SBSmTPCHit::~G4SBSmTPCHit()
{;}

G4SBSmTPCHit::G4SBSmTPCHit(const G4SBSmTPCHit &hit ) : G4VHit() //copy constructor:
{
  //copy all private data members:
  fTrackID = hit.GetTrackID();
  fMotherID = hit.GetMotherID();
  fTrackPID = hit.GetTrackPID();
  // fMotherPID = hit.GetMotherPID();
  fpos = hit.GetPos();
  fLpos = hit.GetLPos();
  fdirection = hit.GetDirection();
  fLdirection = hit.GetLDirection();
  fvertex = hit.GetVertex();
  // fMothervertex = hit.GetMotherVertex();
  fedep = hit.GetEdep();
  fdx = hit.Getdx();
  fenergy = hit.GetEnergy();
  ftime = hit.GetTime();
  fmom = hit.GetMom();
  fLmom = hit.GetLMom();
  fCell = hit.GetCell();
  // fPMTnumber = hit.GetPMTnumber();
  // frownumber = hit.Getrownumber();
  // fcolnumber = hit.Getcolnumber();

  // fdirection = hit.GetDirection();
  // fvertexdirection = hit.GetVertexDirection();
  // fMothervertexdirection = hit.GetMotherVertexDirection();

  // foriginvol = hit.GetOriginVol();
  fCellCoord = hit.GetCellCoord();
  fGlobalCellCoord = hit.GetGlobalCellCoord();
  fvolume_log = hit.GetLogicalVolume();
  fZTravel = hit.GetZTravel();
  fNStrips = hit.GetNStrips();

  otridx = hit.GetOTrIdx();
  ptridx = hit.GetPTrIdx();
  sdtridx = hit.GetSDTrIdx();
}

const G4SBSmTPCHit& G4SBSmTPCHit::operator=(const G4SBSmTPCHit &hit) //assignment operator
{
  //copy all private data members and return the object pointed by this:
  fTrackID = hit.GetTrackID();
  fMotherID = hit.GetMotherID();
  fTrackPID = hit.GetTrackPID();
  fpos = hit.GetPos();
  fLpos = hit.GetLPos();
  fdirection = hit.GetDirection();
  fLdirection = hit.GetLDirection();
  fvertex = hit.GetVertex();
  fedep = hit.GetEdep();
  fdx = hit.Getdx();
  fenergy = hit.GetEnergy();
  ftime = hit.GetTime();
  fmom = hit.GetMom();
  fLmom = hit.GetLMom();
  fCell = hit.GetCell();
  fCellCoord = hit.GetCellCoord();
  fGlobalCellCoord = hit.GetGlobalCellCoord();
  fvolume_log = hit.GetLogicalVolume();
  fZTravel = hit.GetZTravel();
  fNStrips = hit.GetNStrips();
  
  otridx = hit.GetOTrIdx();
  ptridx = hit.GetPTrIdx();
  sdtridx = hit.GetSDTrIdx();

  return *this;
}

G4int G4SBSmTPCHit::operator==(const G4SBSmTPCHit &hit ) const
{
  return (this==&hit) ? 1 : 0;
}

//std::map<G4String, G4AttDef> G4SBSmTPCHit::fAttDefs;

void G4SBSmTPCHit::Draw() 
{
  G4VVisManager *pVVisManager = G4VVisManager::GetConcreteInstance();
  if( pVVisManager ){
    G4Circle circle( fpos );
    circle.SetScreenSize(5);
    circle.SetFillStyle( G4Circle::filled );
    G4VisAttributes attribs(G4Colour(0.0,1.0,0.0) );
    circle.SetVisAttributes( attribs );
    pVVisManager->Draw(circle);
  }
}

void G4SBSmTPCHit::Print() {;}

