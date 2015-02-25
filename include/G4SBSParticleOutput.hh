#ifndef G4SBSParticleOutput_h
#define G4SBSParticleOutput_h 1

#include <vector>

using namespace std;

class G4SBSParticleOutput {
public:
  G4SBSParticleOutput();
  ~G4SBSParticleOutput();
  void Clear();

  //The goal of this class is to store vertex, momentum and PID information for all particles directly or indirectly involved in producing hits in 
  //sensitive detectors:
  //This is most useful for background simulations:
  int npart;
  vector<int> PID, MID, TID; //Particle ID, mother track ID, Track ID
  vector<int> nbounce; //Number of interactions ("bounces") from primary mother particle to this particle's role in causing a hit in a sensitive detector
  vector<int> hitindex;
  vector<double> vx,vy,vz;
  vector<double> px, py, pz;
};

#endif
