#include "G4SBSRICHHit.hh"

#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4ios.hh"

G4Allocator<G4SBSRICHHit> *G4SBSRICHHitAllocator;

G4SBSRICHHit::G4SBSRICHHit() : G4VHit() //Is invocation of the G4VHit() constructor needed?
{;}

G4SBSRICHHit::~G4SBSRICHHit()
{;}

G4SBSRICHHit::G4SBSRICHHit(const G4SBSRICHHit &hit ) : G4VHit() //copy constructor:
{
  //copy all private data members:
  fTrackID = hit.GetTrackID();
  fMotherID = hit.GetMotherID();
  fTrackPID = hit.GetTrackPID();
  fMotherPID = hit.GetMotherPID();
  fpos = hit.GetPos();
  fLpos = hit.GetLPos();
  fvertex = hit.GetVertex();
  fMothervertex = hit.GetMotherVertex();
  fedep = hit.GetEdep();
  fdx = hit.Getdx();
  fenergy = hit.GetEnergy();
  ftime = hit.GetTime();
  fPMTnumber = hit.GetPMTnumber();
  frownumber = hit.Getrownumber();
  fcolnumber = hit.Getcolnumber();

  fdirection = hit.GetDirection();
  fvertexdirection = hit.GetVertexDirection();
  fMothervertexdirection = hit.GetMotherVertexDirection();

  //foriginvol = hit.GetOriginVol();
  fQuantumEfficiency = hit.GetQuantumEfficiency();

  otridx = hit.GetOTrIdx();
  ptridx = hit.GetPTrIdx();
  sdtridx = hit.GetSDTrIdx();
}

const G4SBSRICHHit& G4SBSRICHHit::operator=(const G4SBSRICHHit &hit) //assignment operator
{
  //copy all private data members and return the object pointed by this:
  fTrackID = hit.GetTrackID();
  fMotherID = hit.GetMotherID();
  fTrackPID = hit.GetTrackPID();
  fMotherPID = hit.GetMotherPID();
  fpos = hit.GetPos();
  fLpos = hit.GetLPos();
  fvertex = hit.GetVertex();
  fMothervertex = hit.GetMotherVertex();
  fedep = hit.GetEdep();
  fdx = hit.Getdx();
  fenergy = hit.GetEnergy();
  ftime = hit.GetTime();
  fPMTnumber = hit.GetPMTnumber();
  frownumber = hit.Getrownumber();
  fcolnumber = hit.Getcolnumber();

  fdirection = hit.GetDirection();
  fvertexdirection = hit.GetVertexDirection();
  fMothervertexdirection = hit.GetMotherVertexDirection();

  //foriginvol = hit.GetOriginVol();
  fQuantumEfficiency = hit.GetQuantumEfficiency();

  otridx = hit.GetOTrIdx();
  ptridx = hit.GetPTrIdx();
  sdtridx = hit.GetSDTrIdx();
  
  return *this;
}

G4int G4SBSRICHHit::operator==(const G4SBSRICHHit &hit ) const
{
  return (this==&hit) ? 1 : 0;
}

//std::map<G4String, G4AttDef> G4SBSRICHHit::fAttDefs;

void G4SBSRICHHit::Draw() 
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

void G4SBSRICHHit::Print() {;}


// G4int G4SBSRICHHit::calc_row( G4int PMT ){
//   //Bad practice to have these hard-coded, since we want to study the effects of changes later, but for now we go ahead:
//   //G4int nrow[2] = {26, 27};
  
//   //total of 73 columns, 

//   //G4int super_col = PMT/53;
//   G4int super_row = PMT%53; //ranges from 0..52
//   //if 
//   G4int sub_col = super_row/26;
//   G4int sub_row = super_row%26; 
  
//   return sub_row + 26*(sub_col/2); //sub_col/2 returns 0 unless super_row==52, in which case it returns 1, so for this one result row number == 26 
// }

// G4int G4SBSRICHHit::calc_col( G4int PMT ){
//   G4int super_col = PMT/53; 
//   G4int super_row = PMT%53; 
//   G4int sub_col = (super_row-1)/26;
//   //G4int sub_row = super_row%26;

//   return 2*super_col + sub_col; //should have value from 0..72
// }
