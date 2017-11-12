#ifndef G4SBSSEAMUSMC_HH
#define G4SBSSEAMUSMC_HH

#include <vector>
#include <iostream>
#include "TString.h"

// This class is very similar to ESEPP, but there are some subtleties. 
// I'd rather reuse code then attempt to get clever and miss something important

class G4SBSSeamusMC{
private:

  // List of outgoing particle energies (GeV):
  // fEp = Outgoing Electron
  // fPp = Outgoing Nucleon
  // fGp = Outgoing Pion (if applicable)

  std::vector<double> fEp, fPp, fGp;
  
  // List of outgoing angles (Radians):
  std::vector<double> fEth, fEphi, fPth, fPphi, fGth, fGphi;

  // Initial Nucleon changes between neutron & proton
  // Final pion can be pi+, pi-, pi0;
  // EventType means elastic or inelastic (no pion generated vs. pion generated)
  std::vector<int> fInitialNucleon, fFinalPion, fEventType;
  std::vector<double> fdsdx;
 
  TString fFileName;      // Name of File to Open
  int fUserEvents;        // Nevents defined in .mac file
  int fMin, fMax;         // User may specify a range of events within file, useful for the Farm or multiprocess hack  
  bool fRange;            // An event range has been specified => changes the way the input file gets handled 
  int fNevents;           // N events actually found in .dat files

  void OpenFile(TString);
  void Clear();   
  
public:
  G4SBSSeamusMC();
  ~G4SBSSeamusMC(){};

  // Global index to access vectors on an event basis:
  static unsigned int fEvent;

  // Keep track of global index within EventGen:
  void SetUpNextEvent(){ fEvent++; }
  int CurrentEvent(){ return fEvent; }
  int GetEventN() { return fNevents; }

  void LoadFile(int, TString, int, int);
  
  // Access containers:
  double GetEprime(unsigned int index);
  double GetEtheta(unsigned int index);
  double GetEphi(unsigned int index);
  double GetPprime(unsigned int index);
  double GetPtheta(unsigned int index);
  double GetPphi(unsigned int index);
  double GetGprime(unsigned int index);
  double GetGtheta(unsigned int index);
  double GetGphi(unsigned int index);
  double Getdsdx(unsigned int index);

  int GetInitialNucleon(unsigned int index);
  int GetFinalPion(unsigned int index);
  int GetEventType(unsigned int index);
};
#endif
