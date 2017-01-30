#include "G4SBSESEPP.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>

unsigned int G4SBSESEPP::fEvent = 0;

G4SBSESEPP::G4SBSESEPP(){
  Clear();
  LoadESEPPOutput("test_v1_e-.dat");
}

void G4SBSESEPP::Clear(){
  fEp.clear();
  fPp.clear();
  fGp.clear();
  fEth.clear();
  fEphi.clear();
  fPth.clear();
  fPphi.clear();
  fGth.clear();
  fGphi.clear();
  return;
}

void G4SBSESEPP::LoadESEPPOutput(TString esepp){

  TString temp = "database/" + esepp;

  const char* filename = temp;

  std::ifstream input(filename);
  if( input.is_open() ){
    // Energies of outgoing electron, hadron, and brem photon:
    double Eprime, Pprime, Gprime; 
  
    // The outgoing angles for the three particles:
    double Eth, Ephi;
    double Pth, Pphi;
    double Gth, Gphi;

    Eprime = 0.0;
    Pprime = 0.0;
    Gprime = 0.0;
    Eth = 0.0;
    Ephi = 0.0;
    Pth = 0.0;
    Pphi = 0.0;
    Gth = 0.0;
    Gphi = 0.0;

    while(input 
	  >> Eprime >> Eth >> Ephi 
	  >> Pprime >> Pth >> Pphi 
	  >> Gprime >> Gth >> Gphi ){
      fEp.push_back( Eprime );
      fPp.push_back( Pprime );
      fGp.push_back( Gprime );
      fEth.push_back( Eth );
      fEphi.push_back( Ephi );
      fPth.push_back( Pth );
      fPphi.push_back( Pphi );
      fGth.push_back( Gth );
      fGphi.push_back( Gphi );
    }
  } else {
    std::cerr << "Error: " << esepp << " within G4SBSESEPP.cc did not "
	      << "open properly. Check for typos..." << std::endl;
  }
  input.close();

  std::cout << "Loaded " << fEp.size() << " ESEPP events from " 
	    << temp << std::endl;

  return;
}

double G4SBSESEPP::GetEprime(unsigned int index){
  if( index >= fEp.size() ) {
    std::cerr << "Error accessing EPrime due to bad index!" << std::endl;
    return 0.0;
  } else {
    return fEp[index];
  }
}

double G4SBSESEPP::GetEtheta(unsigned int index){
  if( index >= fEth.size() ) {
    std::cerr << "Error accessing Etheta due to bad index!" << std::endl;
    return 0.0;
  } else {
    return fEth[index];
  }
}

double G4SBSESEPP::GetEphi(unsigned int index){
  if( index >= fEphi.size() ) {
    std::cerr << "Error accessing EPrime due to bad index!" << std::endl;
    return 0.0;
  } else {
    return fEphi[index];
  }
}
