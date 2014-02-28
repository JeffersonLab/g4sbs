#ifndef G4SBSTrackerOutput_hh 
#define G4SBSTrackerOutput_hh 1

#include <vector>

using namespace std;

class G4SBSTrackerOutput {
public:
  G4SBSTrackerOutput();
  ~G4SBSTrackerOutput();

  void Clear();

  // What variables should be stored here? 
  // The tracks created by primary particles, and ONLY by primary particles, in all defined tracker modules:
  // An implicit assumption is made that tracks are straight lines
  int ntracks; //Number of tracks made by primary particles in tracker modules:
  vector<int> TrackerID; //Tracker ID number of the track
  vector<int> TrackTID; //Track ID of the particle that caused this track
  vector<int> TrackPID; //PDG encoding of particle type which caused this track 
  //No need to store mother ID info because we only output tracks made by primary particles
  vector<int> NumHits; //Total number of "hits" (actually, tracking steps that result in energy deposition in a sensitive volume) in this tracker
  vector<int> NumPlanes; //Total number of unique planes with at least one "hit" in this track.
  vector<int> NDF;
  vector<double> Chi2fit; 
  vector<double> Chi2true;

  //"True" track info:
  vector<double> TrackX; //Track x at "focal plane" in local tracker coordinates
  vector<double> TrackY; //Track y at "focal plane" in local tracker coordinates
  vector<double> TrackXp; //Track dx/dz in local tracker coordinates
  vector<double> TrackYp; //Track dy/dz in local tracker coordinates
  vector<double> TrackT; //Track time at "focal plane" in spectrometer coordinates
  
  //"Reconstructed" track info:
  vector<double> TrackXfit; //reconstructed (straight-line fit to hits smeared by position resolution)
  vector<double> TrackYfit; //reconstructed 
  vector<double> TrackXpfit; //reconstructed 
  vector<double> TrackYpfit; //reconstructed
  //vector<double> TrackTfit; //reconstructed

};

#endif
