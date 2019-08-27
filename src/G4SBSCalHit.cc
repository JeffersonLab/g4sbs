#include "G4SBSCalHit.hh"

#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4ios.hh"
#include "G4Circle.hh"

G4Allocator<G4SBSCalHit> G4SBSCalHitAllocator;

G4SBSCalHit::G4SBSCalHit()
{pos = G4ThreeVector();
  energy = 0;
  edep = 0.0;
  mid = -1;
  pid = -1e9;
}

G4SBSCalHit::~G4SBSCalHit()
{;}
//copy constructor:
G4SBSCalHit::G4SBSCalHit(const G4SBSCalHit &right)
  : G4VHit()
{
  pos = right.pos;
  vertex = right.vertex;
  labpos = right.labpos;
  mom = right.mom;
  
  hittime = right.hittime;
  edep = right.edep;
  energy = right.energy;

  cell = right.cell;
  row = right.row;
  col = right.col;
  plane = right.plane;
  wire = right.wire;
  CellCoords = right.CellCoords;
  GlobalCellCoords = right.GlobalCellCoords;

  pid = right.pid;
  mid = right.mid;
  trid = right.trid;

  //Original track info:
  otridx = right.otridx;
  
  //Primary track info:
  ptridx = right.ptridx;
  
  //SD boundary crossing track info:

  sdtridx = right.sdtridx;
  
}
//assignment operator:
const G4SBSCalHit& G4SBSCalHit::operator=(const G4SBSCalHit &right)
{
  pos = right.pos;
  vertex = right.vertex;
  labpos = right.labpos;
  mom = right.mom;
  
  hittime = right.hittime;
  edep = right.edep;
  energy = right.energy;

  cell = right.cell;
  row = right.row;
  col = right.col;
  plane = right.plane;
  wire = right.wire;
  CellCoords = right.CellCoords;
  GlobalCellCoords = right.GlobalCellCoords;

  pid = right.pid;
  mid = right.mid;
  trid = right.trid;

  //Original track info:
  otridx = right.otridx;
  
  //Primary track info:
  ptridx = right.ptridx;
  
  //SD boundary crossing track info:

  sdtridx = right.sdtridx;
 
  
  return *this;
}

G4int G4SBSCalHit::operator==(const G4SBSCalHit &right) const
{
  return (this==&right) ? 1 : 0;
}

void G4SBSCalHit::Draw()
{
  G4VVisManager *pVVisManager = G4VVisManager::GetConcreteInstance();
  if( pVVisManager ){
    G4Circle circle( labpos );
    circle.SetScreenSize(5);
    circle.SetFillStyle( G4Circle::filled );
    G4VisAttributes attribs(G4Colour(0.0,1.0,0.0) );
    circle.SetVisAttributes( attribs );
    pVVisManager->Draw(circle);
  }
}

void G4SBSCalHit::Print()
{
}

// void G4SBSCalHit::SetOriginalTrackInformation( G4Track *oTrack ){
//   G4SBSTrackInformation *oTrackInfo = (G4SBSTrackInformation*) oTrack->GetUserInformation();

//   otrid = oTrackInfo->GetOriginalTrackID();
//   opid  = oTrackInfo->GetOriginalDefinition()->GetPDGEncoding();
//   omid  = oTrackInfo->GetOriginalParentID();

//   opos  = oTrackInfo->GetOriginalPosition();
//   omom  = oTrackInfo->GetOriginalMomentum();
//   opol  = oTrackInfo->GetOriginalPolarization();

//   oenergy = oTrackInfo->GetOriginalEnergy();
//   otime   = oTrackInfo->GetOriginalTime();

// }

// void G4SBSCalHit::SetPrimaryTrackInformation( G4Track *pTrack ){
//   G4SBSTrackInformation *pTrackInfo = (G4SBSTrackInformation*) pTrack->GetUserInformation();
  
//   ptrid = oTrackInfo->GetPrimaryTrackID();
//   ppid  = oTrackInfo->GetPrimaryDefinition()->GetPDGEncoding();
//   //pmid  = oTrackInfo->GetPrimaryParentID();

//   ppos  = oTrackInfo->GetPrimaryPosition();
//   pmom  = oTrackInfo->GetPrimaryMomentum();
//   ppol  = oTrackInfo->GetPrimaryPolarization();

//   penergy = oTrackInfo->GetPrimaryEnergy();
//   ptime   = oTrackInfo->GetPrimaryTime();

// }

// void G4SBSCalHit::SetSDTrackInformation( G4Track *sTrack, G4String SDname ){
//   G4SBSTrackInformation *sTrackInfo = (G4SBSTrackInformation*) sTrack->GetUserInformation();

//   set<G4String> SDlist = sTrackInfo->fSDlist; 

//   if( SDlist.find( SDname ) != SDlist.end() ){
//     sdtrid = (sTrackInfo->fSDTrackID)[SDname];
//     sdmid  = (sTrackInfo->fSDParentID)[SDname];
//     sdpid = (sTrackInfo->fSDDefinition)[SDname]->GetPDGEncoding();
//     sdpos = (sTrackInfo->fSDPosition)[SDname];
//     sdmom = (sTrackInfo->fSDMomentum)[SDname];
//     sdpol = (sTrackInfo->fSDPolarization)[SDname];
//     sdenergy = (sTrackInfo->fSDEnergy)[SDname];
//     sdtime = (sTrackInfo->fSDTime)[SDname];
//   }
  
// }
