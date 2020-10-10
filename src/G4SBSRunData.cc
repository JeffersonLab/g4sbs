#include "G4SBSRunData.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
//#include "G4UImanager.hh"

#include <string.h>
#include <sstream>
#include <errno.h>
#include <sys/stat.h>

#include "TObjArray.h"
#include "TObjString.h"

#ifdef __APPLE__
#include "unistd.h"
#endif

G4SBSRunData::G4SBSRunData(){
    fNthrown = -1;
    fNtries = -1;
    fBeamE   = -1e9;
    fExpType[0]  = '\0';
    fGenName[0]  = '\0';
    fGitInfo[0]  = '\0';
    fHostName[0] = '\0';
    fMagData.clear();
}

G4SBSRunData::~G4SBSRunData(){
}

void G4SBSRunData::Init(){
    fNthrown = 0;
    fNtries = 0;
    fBeamE   = 0;
    fNormalization = 1.0;
    fGenVol = 4.0*pi;
    fLuminosity = 0.0;
    fBBtheta = 30.0*degree;
    fBBdist = 1.5*meter;
    fSBStheta = 15.0*degree;
    fSBSdist  = 3.0*meter;
    fHCALdist = 10.0*meter;
    fHCALvoff = 0.0*centimeter;
    fRICHdist = 5.5*meter;
    fSBSTrackerdist = 4.5*meter;
    fSBSTrackerPitch = 0.0*degree;
    
    strcpy(fExpType, "default");
    strcpy(fGenName, "default");
    strcpy(fGitInfo, gGitInfoStr);
    if(gethostname(fHostName,__RUNSTR_LEN) == -1){
	fprintf(stderr, "%s line %d: ERROR could not get hostname\n", __PRETTY_FUNCTION__ ,  __LINE__ );
	fprintf(stderr, "%s\n",strerror(errno));
    }
    if(getcwd(fRunPath,__RUNSTR_LEN) == NULL){
	fprintf(stderr, "%s line %d: ERROR could not get current working directory\n", __PRETTY_FUNCTION__ ,  __LINE__ );
	fprintf(stderr, "%s\n",strerror(errno));
    }
}

void G4SBSRunData::Print() const { Print(NULL); }

void G4SBSRunData::Print(Option_t *) const {
    printf("git repository info\n-------------------------------------------------\n%s-------------------------------------------------\n\n", fGitInfo);
    printf("Run at %s on %s\n", fRunTime.AsString("ls"), fHostName);
    printf("Run Path %s\n", fRunPath);
    printf("N generated = %ld\n", fNthrown);
    printf("N tries     = %ld\n", fNtries);
    printf("Beam Energy = %f GeV\n", fBeamE);
    printf("Experiment  = %s\n", fExpType);
    printf("Generator   = %s\n", fGenName);
    printf("Normalization = %8.5g\n", fNormalization);
    printf("Rejection sampling OFF: rate = dsigma * Normalization = dsigma x (phsp. vol) x luminosity/Ntries\n");
    printf("Rejection sampling ON: rate = (max xsec * Normalization) = maxweightx(phsp.vol)xluminosity/Ntries\n");
    printf("Generation Volume = %8.5g\n", fGenVol);
    printf("Luminosity = %8.5g Hz/cm^2\n", fLuminosity );
    printf("BigBite central angle = %8.5g degrees\n", fBBtheta/degree );
    printf("BigBite magnet distance = %8.5g meters\n", fBBdist );
    printf("SBS central angle = %8.5g degrees\n", fSBStheta/degree );
    printf("SBS magnet distance = %8.5g meters\n", fSBSdist );
    printf("HCAL distance = %8.5g meters\n", fHCALdist );
    printf("HCAL vertical offset = %8.5g cm\n", fHCALvoff*m/cm );
    printf("RICH distance = %8.5g meters\n", fRICHdist );
    printf("SBS tracker distance = %8.5g meters\n", fSBSTrackerdist );
    printf("SBS tracker pitch angle = %8.5g degrees\n", fSBSTrackerPitch/degree );
    printf("Max weight for rejection sampling = %8.5g (diff. xsec units) \n", fMaxWeight );
    
    printf("Field maps:\n");
    for( unsigned int i = 0; i < fMagData.size(); i++ ){
	printf("\t%s\n", fMagData[i].filename);
	printf("\t%s\n", fMagData[i].hashsum);
	printf("\t%s\n\n", fMagData[i].timestamp.AsString("ls"));
    }

    printf("Pre-Init Macro run:\n---------------------------------------\n");

    fPreInitMacro.Print();

    printf("Macro run:\n-------------------------------------------------\n");

    fMacro.Print();

    printf("External macros called (%3d total):\n------------------------\n",
        int(fExternalMacros.size()));
    for(size_t i = 0; i < fExternalMacros.size(); i++) {
      printf("External Macro %3d out of %3d:\n",
          int(i+1),int(fExternalMacros.size()));
      fExternalMacros[i].Print();
    }

}

TString G4SBSRunData::FindMacro( const char *fn ){

  struct stat fdata_temp;
  
  bool found_file = (stat( fn, &fdata_temp ) == 0 );

  if( found_file ){
    return TString(fn);
  } else {
    TObjArray *tokens = fMacroPath.Tokenize(":");
    for( int n=0; n<tokens->GetEntries(); n++ ){
      TString spathtemp = ( (TObjString*) (*tokens)[n] )->GetString();

      TString fnametemp = spathtemp + "/" + fn;

      found_file = (stat( fnametemp.Data(), &fdata_temp ) == 0 );

      if( found_file ) {
	return fnametemp;
      }
    }
    return TString(fn);
  }
}

void G4SBSRunData::SetMacroFile( const char *fn ){

  TString fullpathname = FindMacro( fn );
  
  fMacro = G4SBSTextFile( fullpathname.Data() );
  FindExternalMacros( fMacro );
}

void G4SBSRunData::SetPreInitMacroFile( const char *fn ){
  
  TString fullpathname = FindMacro( fn );

  fPreInitMacro = G4SBSTextFile( fullpathname.Data() );
  FindExternalMacros( fPreInitMacro );
}

void G4SBSRunData::FindExternalMacros(G4SBSTextFile macro){
  // Turn the file into an input stream for easier parsing
  std::istringstream macro_string(macro.GetBuffer());

  //  G4UImanager * UImanager = G4UImanager::GetUIpointer();
  
  // Command which calls external macros (this is what we are trying to find)
  std::string cmd("/control/execute ");
  std::string line;
  size_t p; // position of found string
  size_t c; // position of a comment character (#)
  while(std::getline(macro_string,line) && !macro_string.eof()) {
    p = line.find(cmd);
    if(p != std::string::npos) { // cmd found, now see if it is commented out
      c = line.find("#");
      if(c > p ) { // A comment before the cmd was not found, so store macro

	//TString macrofilename = line.substr(p+cmd.length(),line.length()).c_str();
	
	TString fullpathname = FindMacro( line.substr(p+cmd.length(),line.length()).c_str() );
	
	G4SBSTextFile ext_macro( fullpathname.Data() );
	
        // G4SBSTextFile ext_macro(
        //     line.substr(p+cmd.length(),line.length()).c_str());
        fExternalMacros.push_back(ext_macro);
        // Now check if this macro calls another
        FindExternalMacros(ext_macro);
      }
    }
  }
}

void G4SBSRunData::CalcNormalization(){
  SetNormalization( fMaxWeight * fGenVol * fLuminosity / double(fNtries) );
}

ClassImp(G4SBSRunData)











