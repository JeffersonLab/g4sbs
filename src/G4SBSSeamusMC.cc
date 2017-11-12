#include "G4SBSSeamusMC.hh"
#include <fstream>

G4SBSSeamusMC::G4SBSSeamusMC(){
  Clear();
  fFileName = "";
  fUserEvents = 0;
  fNevents = 0;
  fRange = false;
}

void G4SBSSeamusMC::Clear(){
  fEp.clear();
  fPp.clear();
  fGp.clear();
  fEth.clear();
  fEphi.clear();
  fPth.clear();
  fPphi.clear();
  fGth.clear();
  fGphi.clear();
  fdsdx.clear();
  fInitialNucleon.clear();
  fFinalPion.clear();
  fEventType.clear();
}

////////////////////////////////////////////////////////////////////////////////
// The all important global index - accesses all containers within this
// class and gets incremented once per event. Doesn't have to be static...
//
unsigned int G4SBSSeamusMC::fEvent = 0;

////////////////////////////////////////////////////////////////////////////////
// Allow G4SBSEventGen access
//
void G4SBSSeamusMC::LoadFile(int n, TString name, int min, int max){
  fUserEvents = n;
  fFileName = "seamus_mc/" + name;
  fMin = min;
  fMax = max;

  // If uninitialized, range is not specified and deactivated. The signal
  // for deactivation is fMax=0, fMin=0 => fMax-fMin==0 for default behavior
  if( fMax-fMin == 0 ){
    fRange = false;
  } else {
    fRange = true;
  }

  if( fRange ){
    // If fMax is larger than fMin, swap
    if( fMax - fMin < 0 ) {
      int temp = fMax;
      fMax = fMin;
      fMin = temp;
    }
    fUserEvents = fMax - fMin;
  }

  if( fRange ){
    std::cout << "Loading " << fUserEvents 
	      << " Seamus Inelastic MC events in the range [" 
	      << fMin << ", " << fMax << ")..." << std::endl;
  } else {
    std::cout << "Loading " << fUserEvents << " Seamus Inelastic MC events..." << std::endl;
  }

  OpenFile( fFileName );
}

////////////////////////////////////////////////////////////////////////////////
// Open .dat files - handle containers - change N events to be generated
// if necessary
//
void  G4SBSSeamusMC::OpenFile(TString name){
  const char* filename = name;
  int temp = 0;  // line # for accepted data
  int dummy = 0; // line # in text file
  std::ifstream input(filename);

  double GeV2MeV = 1000.0;

  if( input.is_open() ){ 
    double Eprime, Pprime, Gprime; 
    double Eth, Ephi, Pth, Pphi, Gth, Gphi;
    double cross_section;
    int Nucl, Pion, Type;

    Eprime = Pprime = Gprime = 0.0;
    Eth = Pth = Gth = 0.0;
    Ephi = Pphi = Gphi = 0.0;
    Nucl = Pion = Type = 0;
    cross_section = 0.0;

    while(input >> Eprime >> Eth >> Ephi 
	  >> Pprime >> Pth >> Pphi 
	  >> Gprime >> Gth >> Gphi 
	  >> Nucl >> Pion >> Type >> cross_section ){

      // I am breaking it up this way so no sleight-of-hand is performed and 
      // the logic is crystal clear.
      // Only accept lines of an ESEPP text file within a particular range
      // if fRange == true, otherwise, take everything until end of file or
      // code has reached user's input 

      // If a range is specified
      if( fRange ){
	if( dummy >= fMin && dummy < fMax ){
	  // Only take # of requested events
	  if( temp >= fUserEvents ) break;

	  fEp.push_back( Eprime*GeV2MeV );
	  fPp.push_back( Pprime*GeV2MeV );
	  fGp.push_back( Gprime*GeV2MeV );

	  fEth.push_back( Eth );
	  fPth.push_back( Pth );
	  fGth.push_back( Gth );

	  fEphi.push_back( Ephi );
	  fPphi.push_back( Pphi );	 
	  fGphi.push_back( Gphi );

	  fInitialNucleon.push_back( Nucl );
	  fFinalPion.push_back( Pion );
	  fEventType.push_back( Type );

	  fdsdx.push_back( cross_section );
	  printf("%d \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %1.2f \t %d \t %d \t %d \n", 
	  	 fMin+temp, Eprime, Eth, Ephi, Pprime, Pth, Pphi, Gprime, Gth, Gphi, Nucl, Pion, Type);
	  temp++;
	}
      } else {
	// Only take # of requested events
	if( temp >= fUserEvents ) break;

	fEp.push_back( Eprime );
	fPp.push_back( Pprime );
	fGp.push_back( Gprime );
	fEth.push_back( Eth );
	fEphi.push_back( Ephi );
	fPth.push_back( Pth );
	fPphi.push_back( Pphi );
	fGth.push_back( Gth );
	fGphi.push_back( Gphi );
	fInitialNucleon.push_back( Nucl );
	fFinalPion.push_back( Pion );
	fEventType.push_back( Type );
	fdsdx.push_back( cross_section );
	temp++;
      }
      // Count all lines regardless for a range of events:
      dummy++;
    }
  } else {
    std::cerr << "Error: " << name << " did not "
	      << "open properly. Check for typos in .mac file or " 
	      << "add the /g4sbs/SeamusMCfile <filename> command. Put file "
	      << "in the seamus_mc/ directory." << std::endl;
  }
  input.close();

  if( temp != fUserEvents ) {
    std::cout << "Notice: Seamus MC file " << name << " does not have "
	      << fUserEvents << " events, and will only be generating " 
	      << temp << " events." << std::endl;
  }
  
  std::cout << "Loaded " << temp << " Seamus MC inelastic events from " 
	    << name << "." << std::endl;
  
  fNevents = temp;
}

////////////////////////////////////////////////////////////////////////////////
// Some get-methods to access the containers
//
double G4SBSSeamusMC::GetEprime(unsigned int index){
  if( index >= fEp.size() ) {
    std::cerr << "Error accessing Eprime due to bad index!" << std::endl;
    return 0.0;
  } else {
    return fEp[index];
  }
}

double G4SBSSeamusMC::GetEtheta(unsigned int index){
  if( index >= fEth.size() ) {
    std::cerr << "Error accessing Etheta due to bad index!" << std::endl;
    return 0.0;
  } else {
    return fEth[index];
  }
}

double G4SBSSeamusMC::GetEphi(unsigned int index){
  if( index >= fEphi.size() ) {
    std::cerr << "Error accessing Ephi due to bad index!" << std::endl;
    return 0.0;
  } else {
    return fEphi[index];
  }
}

double G4SBSSeamusMC::GetPprime(unsigned int index){
  if( index >= fPp.size() ) {
    std::cerr << "Error accessing Pprime due to bad index!" << std::endl;
    return 0.0;
  } else {
    return fPp[index];
  }
}

double G4SBSSeamusMC::GetPtheta(unsigned int index){
  if( index >= fPth.size() ) {
    std::cerr << "Error accessing Ptheta due to bad index!" << std::endl;
    return 0.0;
  } else {
    return fPth[index];
  }
}

double G4SBSSeamusMC::GetPphi(unsigned int index){
  if( index >= fPphi.size() ) {
    std::cerr << "Error accessing Pphi due to bad index!" << std::endl;
    return 0.0;
  } else {
    return fPphi[index];
  }
}

double G4SBSSeamusMC::GetGprime(unsigned int index){
  if( index >= fGp.size() ) {
    std::cerr << "Error accessing Gprime due to bad index!" << std::endl;
    return 0.0;
  } else {
    return fGp[index];
  }
}

double G4SBSSeamusMC::GetGtheta(unsigned int index){
  if( index >= fGth.size() ) {
    std::cerr << "Error accessing Gtheta due to bad index!" << std::endl;
    return 0.0;
  } else {
    return fGth[index];
  }
}

double G4SBSSeamusMC::GetGphi(unsigned int index){
  if( index >= fGphi.size() ) {
    std::cerr << "Error accessing Gphi due to bad index!" << std::endl;
    return 0.0;
  } else {
    return fGphi[index];
  }
}

// These methods require different protocols b/c a value of 0 may be normal:
int G4SBSSeamusMC::GetInitialNucleon(unsigned int index){
  if( index >= fInitialNucleon.size() ) {
    std::cerr << "Error accessing fInitialNucleon due to bad index!" << std::endl;
    return -100;
  } else {
    return fInitialNucleon[index];
  }
}

int G4SBSSeamusMC::GetFinalPion(unsigned int index){
  if( index >= fFinalPion.size() ) {
    std::cerr << "Error accessing fFinalPion due to bad index!" << std::endl;
    return -100;
  } else {
    return fFinalPion[index];
  }
}

int G4SBSSeamusMC::GetEventType(unsigned int index){
  if( index >= fEventType.size() ) {
    std::cerr << "Error accessing fFinalPion due to bad index!" << std::endl;
    return -100;
  } else {
    return fEventType[index];
  }
}

double G4SBSSeamusMC::Getdsdx(unsigned int index){
  if( index >= fdsdx.size() ) {
    std::cerr << "Error accessing fdsdx due to bad index!" << std::endl;
    return 0.0;
  } else {
    return fdsdx[index];
  }
}
