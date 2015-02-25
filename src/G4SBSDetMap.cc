#include "G4SBSDetMap.hh"

G4SBSDetMap::G4SBSDetMap(){
  SDname = "";
  clear();
}

G4SBSDetMap::G4SBSDetMap( G4String name ){
  SDname = name;
  clear();
}

G4SBSDetMap::~G4SBSDetMap()
{;}

void G4SBSDetMap::clear(){
  depth = 0;
  Plane.clear();
  Wire.clear();
  Row.clear();
  Col.clear();
  LocalCoord.clear();
  //GlobalCoord.clear();
}
