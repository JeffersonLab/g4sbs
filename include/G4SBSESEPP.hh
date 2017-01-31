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
  bool fRosenbluth;       // Do we want Rosenbluth events as well?
  bool fRad;              // Do we want Radiative Events?
  int fUserEvents;        // Nevents defined in .mac file
  
  void LoadRadEvents();
  void LoadRosenbluthEvents();
  void Clear();
  void OpenFile(TString);

public:
  G4SBSESEPP(bool,bool); // # of events, name of file, rad, rosenbluth
  ~G4SBSESEPP(){};

  void SetUpNextEvent(){ fEvent++; }
  int CurrentEvent(){ return fEvent; }
  void LoadFiles(int);
  void SetName(TString);

  // Global index to access vectors on an event basis
  static unsigned int fEvent;

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
