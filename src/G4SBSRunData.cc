#include "G4SBSRunData.hh"

#include <string.h>
#include <sstream>
#include <errno.h>

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
}

G4SBSRunData::~G4SBSRunData(){
}

void G4SBSRunData::Init(){
    fNthrown = 0;
    fNtries = 0;
    fBeamE   = 0;
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


void G4SBSRunData::FindExternalMacros(G4SBSTextFile macro){
  // Turn the file into an input stream for easier parsing
  std::istringstream macro_string(macro.GetBuffer());

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
        G4SBSTextFile ext_macro(
            line.substr(p+cmd.length(),line.length()).c_str());
        fExternalMacros.push_back(ext_macro);
        // Now check if this macro calls another
        FindExternalMacros(ext_macro);
      }
    }
  }
}

ClassImp(G4SBSRunData)











