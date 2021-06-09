#ifndef G4SBSAcquMCOutput_h
#define G4SBSAcquMCOutput_h 1

#include <vector>
using namespace std;

class G4SBSAcquMCOutput {
public:
  G4SBSAcquMCOutput();
  ~G4SBSAcquMCOutput();
  void Clear();
  void ConvertToTreeUnits();

  // vertex, momentum, total energymomentum and energy of primaries
  double Vx, Vy, Vz, Px, Py, Pz, Pt, E;
  int pid;
  // double mass;
};

#endif
