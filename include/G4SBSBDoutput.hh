#ifndef G4SBS_BEAM_DIFFUSER_OUTPUT_HH
#define G4SBS_BEAM_DIFFUSER_OUTPUT_HH

// C++ std output formatting that gets attached 
// to an output ROOT tree branch 

#include <cstdlib>
#include <vector>

class G4SBSBDoutput {

   public:
      G4SBSBDoutput();
      ~G4SBSBDoutput();

      void Clear();

      int nhits;                       // number of hits 
      std::vector<int> plane;          // which plane 
      std::vector<int> trid;           // track ID  
      std::vector<int> pid;            // particle type 
      std::vector<int> mid;            // material type 
      std::vector<double> t,x,y,z;     // time, local coordinates (x,y,z) 
      std::vector<double> xg,yg,zg;    // lab coordinates  
      std::vector<double> p,edep,beta; // momentum, deposited energy, beta = v/c  

      // G4SBSParticleOutput ParticleHistory;

};

#endif
