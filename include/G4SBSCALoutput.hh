#ifndef G4SBSCALoutput_h
#define G4SBSCALoutput_h 1

#include <vector>
#include "G4SBSParticleOutput.hh"
#include "G4SBSSDTrackOutput.hh"

using namespace std;

class G4SBSCALoutput {
public:
  G4SBSCALoutput();
  ~G4SBSCALoutput();
  
  void Clear();
  
  double timewindow, threshold;
  int ntimebins;

  vector< vector<double> > edep_vs_time;  //time dependence of energy deposition within timewindow
  double gatewidth;

  double Esum; //energy deposition sum over entire sensitive detector
  vector<double> timebins; //divide timewindow into ntimebins bins and histogram the energy in each bin (summed over whole detector)
  
  //"Hit" is defined as the sum of all energy deposition in a given cell during timewindow, provided the total energy deposition is above "threshold":
  int nhits_CAL; //Number of cells with above-threshold energy deposition:
  vector<int> row, col, plane, wire, cell; //"row" and "column" of cells
  vector<double> xcell, ycell, zcell; //"local" x and y coordinates of center of cell (assumes a rectangular, planar geometry)
  vector<double> xcellg,ycellg,zcellg; //"global" xyz coordinates of center of cell (no assumptions on geometry)
  vector<double> xhit, yhit, zhit; //weighted "average" local position of energy deposition
  vector<double> xhitg, yhitg, zhitg; //weighted "average" global position of energy deposition
  vector<double> sumedep, tavg, trms, tmin, tmax; //Sum of energy deposition, average, rms, min and max global times of energy depositions in this cell
  
  //"Part" keeps track of all unique particles depositing energy in a "calorimeter" sensitive volume:
  int npart_CAL; //Number of particles depositing energy in a given cell
  vector<int> ihit; //hit index associated with this particle 
  vector<double> x, y, z, t, E, dt, L; //average global xyz coordinates of steps of this particle in this cell, average time, total time, path length, and initial energy of particles
  vector<double> vx, vy, vz; //global vertex coordinates of particles in this hit
  vector<int> trid, mid, pid; //track ID, mother track ID and particle ID info
  vector<double> p, edep; //initial momentum and total energy deposition of particles.
  vector<double> px,py,pz; //initial momentum components

  //Add bookkeeping indices for "original", "primary", and "SD boundary" tracks:
  vector<int> otridx,ptridx,sdtridx;
  
  bool keeppart;
  
  G4SBSParticleOutput ParticleHistory;

};

#endif
