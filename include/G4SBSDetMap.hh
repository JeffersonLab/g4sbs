#ifndef G4SBSDetMap_h
#define G4SBSDetMap_h 1

#include <map>
#include <set>
#include <vector>
#include "G4String.hh"
#include "G4ThreeVector.hh"

using namespace std;

class G4SBSDetMap 
{ //Class to store mapping between "channel numbers" and information related to physical location for sensitive detectors; e.g., rows, columns, x, y, z, depth, etc.
public:
  G4SBSDetMap();
  G4SBSDetMap( G4String SDname );
  ~G4SBSDetMap();

  void clear();

  // Each of the following maps is a map<G4String, map<int,X> > where G4String is the sensitive detector name, 
  // the key is "copy number" for the sensitive detector physical volume, and the additional information can be used to identify the physical location of hits.

  G4String SDname;
  G4int depth;
  
  map<G4int,G4int>  Plane; //"Plane" information (for wire-chamber/tracking type detectors)
  map<G4int,G4int>  Wire;  //"Wire" information 
  map<G4int,G4int>  Row; //"Row" information (for calorimeter type detectors)
  map<G4int,G4int>  Col; //"Column" information (for calorimeters)
  map<G4int,G4ThreeVector>  LocalCoord; //local coordinates of the center of the sensitive volume (relative to some "mother" logical volume that contains the 
                                        //sensitive detector
  //map<G4int,G4ThreeVector>  GlobalCoord; //Global coordinates of the center of the sensitive volume 

};


#endif
