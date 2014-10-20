#include "G4SBSECalHit.hh"

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

G4Allocator<G4SBSECalHit> *G4SBSECalHitAllocator;

G4SBSECalHit::G4SBSECalHit() : G4VHit()
{;}

G4SBSECalHit::~G4SBSECalHit()
{;}

G4SBSECalHit::G4SBSECalHit(const G4SBSECalHit &hit ) : G4VHit()
{
  //copy all private data members:
  fTrackID = hit.GetTrackID();
  //fMotherID = hit.GetMotherID();
  //fTrackPID = hit.GetTrackPID();
  //fMotherPID = hit.GetMotherPID();
  fpos = hit.GetPos();
  //fLpos = hit.GetLPos();
  //fvertex = hit.GetVertex();
  //fMothervertex = hit.GetMotherVertex();
  fedep = hit.GetEdep();
  fdx = hit.Getdx();
  fenergy = hit.Getenergy();
  ftime = hit.GetTime();
  fPMTnumber = hit.GetPMTnumber();
  frownumber = hit.Getrownumber();
  fcolnumber = hit.Getcolnumber();

  //fdirection = hit.GetDirection();
  //fvertexdirection = hit.GetVertexDirection();
  //fMothervertexdirection = hit.GetMotherVertexDirection();

  //foriginvol = hit.GetOriginVol();
 
}

const G4SBSECalHit& G4SBSECalHit::operator=(const G4SBSECalHit &hit)
{
  //copy all private data members and return the object pointed by this:
  fTrackID = hit.GetTrackID();
  //fMotherID = hit.GetMotherID();
  //fTrackPID = hit.GetTrackPID();
  //fMotherPID = hit.GetMotherPID();
  fpos = hit.GetPos();
  //fLpos = hit.GetLPos();
  //fvertex = hit.GetVertex();
  //fMothervertex = hit.GetMotherVertex();
  fedep = hit.GetEdep();
  fdx = hit.Getdx();
  fenergy = hit.Getenergy();
  ftime = hit.GetTime();
  fPMTnumber = hit.GetPMTnumber();
  frownumber = hit.Getrownumber();
  fcolnumber = hit.Getcolnumber();

  //fdirection = hit.GetDirection();
  //fvertexdirection = hit.GetVertexDirection();
  //fMothervertexdirection = hit.GetMotherVertexDirection();

  //foriginvol = hit.GetOriginVol();

  return *this;
}

G4int G4SBSECalHit::operator==(const G4SBSECalHit &hit ) const
{
  return (this==&hit) ? 1 : 0;
}

void G4SBSECalHit::Draw() 
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

void G4SBSECalHit::Print() {;}

G4int G4SBSECalHit::calc_row( G4int PMT ){

  G4int ecal_row = PMT/58; //0-75
  return ecal_row;

}

G4int G4SBSECalHit::calc_col( G4int PMT ){

  G4int ecal_col = PMT%58; //ranges from 0-19 as PMT goes to 20*76=1520 PMTs
  return ecal_col;

}
