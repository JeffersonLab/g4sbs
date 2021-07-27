#ifndef G4SBSmTPCoutput_h
#define G4SBSmTPCoutput_h 1

#include <vector>
#include "G4SBSParticleOutput.hh"
using namespace std;

class G4SBSmTPCoutput {
public:
  G4SBSmTPCoutput();
  ~G4SBSmTPCoutput();
  
  void Clear();
  
  double timewindow, threshold;
  int ntimebins;

  double Esum; //energy deposition sum over entire sensitive detector
  vector<double> timebins; //divide timewindow into ntimebins bins and histogram the energy in each bin (summed over whole detector)
  vector<double> Edep_vs_time; //time dependence of energy deposition within timewindow 
  
  //"Hit" is defined as the sum of all energy deposition in a given cell during timewindow, provided the total energy deposition is above "threshold":
  int nhits_mTPC; //Number of cells with above-threshold energy deposition:
  // vector<int> row, col, plane, wire, cell; //"row" and "column" of cells
  vector<int> cell; //"row" and "column" of cells
  vector<double> xcell, ycell, zcell; //"local" x and y coordinates of center of cell (assumes a rectangular, planar geometry)
  vector<double> xcellg,ycellg,zcellg; //"global" xyz coordinates of center of cell (no assumptions on geometry)
  vector<double> xhit, yhit, zhit; //weighted "average" local position of energy deposition
  vector<double> xhitg, yhitg, zhitg; //weighted "average" global position of energy deposition
  vector<double> sumedep, tavg, trms, tmin, tmax; //Sum of energy deposition, average, rms, min and max global times of energy depositions in this cell
  // vector<double> Ehit;
  // vector<double> pxhit, pyhit, pzhit;
  vector<int> trid, mid, pid; //track ID, mother track ID and particle ID info
  // here we want to add step length of hit
  vector<double> hitL;
  
  //"Part" keeps track of all unique particles depositing energy in a "calorimeter" sensitive volume:
  int npart_mTPC; //Number of particles depositing energy in a given cell
  vector<int> ihit; //hit index associated with this particle 
  vector<double> x, y, z, t, E, dt, L; //average global xyz coordinates of steps of this particle in this cell, average time, total time, path length, and initial energy of particles
  vector<double> vx, vy, vz; //global vertex coordinates of particles in this hit
  vector<int> trid_, mid_, pid_; //track ID, mother track ID and particle ID info
  vector<double> p, edep; //initial momentum and total energy deposition of particles.
  vector<double> px,py,pz; //momentum components at this hit
  vector<double> px_v,py_v,pz_v; //momentum components at vertex
  vector<double> ztravel; //z travel in hit for drift time calc
  vector<int> nstrips; //estimated number of readout pads hit

  // bool keeppart;
  
  G4SBSParticleOutput ParticleHistory;

};

#endif
