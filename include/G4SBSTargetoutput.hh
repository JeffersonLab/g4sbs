// C++ std output formatting that gets attached 
// to an output ROOT tree branch 
#ifndef G4SBS_TARGET_OUTPUT_HH
#define G4SBS_TARGET_OUTPUT_HH

#include <cstdlib>
#include <vector>

#include "G4SBSParticleOutput.hh"

class G4SBSTargetoutput {

   public:
      G4SBSTargetoutput();
      ~G4SBSTargetoutput();

      void Clear();

      int nhits_Target;                // number of hits 
      std::vector<int> trid;           // track ID  
      std::vector<int> pid;            // particle type 
      std::vector<int> mid;            // material type 
      std::vector<double> t,x,y,z;     // time, local coordinates (x,y,z) 
      std::vector<double> xg,yg,zg;    // lab coordinates  
      std::vector<double> p,edep,beta; // momentum, deposited energy, beta = v/c  

      G4SBSParticleOutput ParticleHistory;

};

#endif
