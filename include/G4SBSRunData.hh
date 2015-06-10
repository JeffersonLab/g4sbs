#ifndef __G4SBSRUNDATA_HH
#define __G4SBSRUNDATA_HH

#include "TObject.h"

#include <vector>
#include <string>
#include <sbstypes.hh>

#include "gitinfo.hh"
#include "G4SBSTextFile.hh"

/*!
 * All the information on the run
 * This will get put into the output
 * stream
 */

class TGeoManager;

class G4SBSRunData : public TObject {

public:
  G4SBSRunData();
  ~G4SBSRunData();

  unsigned long long int GetNthrown(){ return fNthrown; }
  void SetNthrown(unsigned long long int n){ fNthrown = n; }
  void SetNtries(unsigned long long int n){ fNtries = n; } //Number of tries to generate Nthrown events.

  void Init();

  void SetGenName(const char *n){ strcpy(fGenName, n); }
  const char *GetGenName(){ return fGenName; }

  void SetExpType(const char *n){ strcpy(fExpType, n); }
  const char *GetExpType(){ return fExpType; }

  void SetBeamE(double E){ fBeamE = E; }
  void SetSeed(unsigned int seed){ fSeed = seed; }

  void AddMagData(filedata_t d){fMagData.push_back(d);}
  void SetMacroFile(const char *fn){ fMacro = G4SBSTextFile(fn); }
  void SetPreInitMacroFile(const char *fn){ fPreInitMacro = G4SBSTextFile(fn); }

  void Print(Option_t *) const;
  void Print() const;

  TTimeStamp fRunTime;

  long int  fNthrown;
  long int fNtries;
  unsigned int  fSeed;
  double fBeamE;
  char fExpType[__RUNSTR_LEN];
  char fGenName[__RUNSTR_LEN];
  char fGitInfo[__GITMAXINFO_SIZE];

  char fHostName[__RUNSTR_LEN];
  char fRunPath[__RUNSTR_LEN];

  G4SBSTextFile              fMacro;
  G4SBSTextFile   fPreInitMacro; //PreInit Commands

  std::vector<filedata_t> fMagData;

  ClassDef(G4SBSRunData, 1);
};

#endif//__G4SBSRUNDATA_HH
