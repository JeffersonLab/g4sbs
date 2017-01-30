#ifndef G4SBSESEPP_hh
#define G4SBSESEPP_hh

#include <vector>
#include <iostream>
#include "TString.h"

class G4SBSESEPP{
private:

public:
  G4SBSESEPP();
  ~G4SBSESEPP(){};
  void Clear();
  void ConvertToTreeUnits();
  void LoadESEPPOutput(TString);

  void SetUpNextEvent(){ fEvent++; }
  int CurrentEvent(){ return fEvent; }

  // List of outgoing particle energies (MeV):
  std::vector<double> fEp, fPp, fGp;
  // List of outgoing angles (Radians):
  std::vector<double> fEth, fEphi, fPth, fPphi, fGth, fGphi;

  static unsigned int fEvent; // Global index to access event level info
  int GetNEvents() { return fEp.size(); }

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
