#ifndef G4SBSESEPP_hh
#define G4SBSESEPP_hh

#include <vector>
#include <iostream>
#include "TString.h"

class G4SBSESEPP{
private:

  // List of outgoing particle energies (MeV):
  // fEp = Outgoing Electron
  // fPp = Outgoing Proton
  // fGp = Outgoing Brem Photon
  std::vector<double> fEp, fPp, fGp; // should be const..

  // List of outgoing angles (Radians):
  std::vector<double> fEth, fEphi, fPth, fPphi, fGth, fGphi;

  TString fRadFileName;   // File name for Radiative Events
  TString fRosenFileName; // File name for Rosenbluth Events
  bool fRad;              // Radiative Flag
  bool fRosenbluth;       // Rosenbluth Flag
  int fUserEvents;        // Nevents defined in .mac file
  
  void LoadRadEvents();        // Easier to read
  void LoadRosenbluthEvents(); // " "
  void Clear();                // 
  void OpenFile(TString);      //

  int fRadEvents,fRosenEvents; // N events actually found in .dat files

public:
  G4SBSESEPP();
  ~G4SBSESEPP(){};

  // Keep track of global index within EventGen:
  void SetUpNextEvent(){ fEvent++; }
  int CurrentEvent(){ return fEvent; }

  // Load files within EventGen:
  void LoadFiles(int,TString,bool,bool);

  // Global index to access vectors on an event basis:
  static unsigned int fEvent;

  // Get # events found in files:
  int GetEventN(){ return fRadEvents + fRosenEvents; }

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
};
#endif
