#ifndef G4SBS_UTILITIES_HH
#define G4SBS_UTILITIES_HH

// a namespace of useful functions that don't necessarily fit into 
// a specific class's scope

#include <cstdlib> 
#include <vector>
#include "G4ThreeVector.hh" 

namespace G4SBS {
   namespace Util {
      void RotateVector(std::vector<G4double> R,G4ThreeVector P,G4ThreeVector &W); 
   } //::Util
} //::G4SBS

#endif 
