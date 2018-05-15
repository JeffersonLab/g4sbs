#ifndef __G4SBSRUNDATA_HH
#define __G4SBSRUNDATA_HH

#include "TObject.h"

#include <vector>
#include <string>
#include <sbstypes.hh>

#include "gitinfo.hh"
#include "G4SBSTextFile.hh"


#define __LINE_STRLEN 4096
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

  void SetNormalization( double N ){ fNormalization = N; }
  void SetGenVol( double V ){ fGenVol = V; }
  void SetMaxWeight( double w ){ fMaxWeight = w; }
  void SetLuminosity( double Lumi ){ fLuminosity = Lumi; }
  void SetBBtheta( double theta ){ fBBtheta = theta; }
  void SetBBdist( double d ){ fBBdist = d; }
  void SetSBStheta( double theta ){ fSBStheta = theta; }
  void SetSBSdist( double d ){ fSBSdist = d; }
  void SetHCALdist( double d ){ fHCALdist = d; }
  void SetHCALvoff( double y ){ fHCALvoff = y; }
  void SetRICHdist( double d ){ fRICHdist = d; }
  void SetSBSTrackerDist( double d ){ fSBSTrackerdist = d; }
  void SetSBSTrackerPitch( double a ){ fSBSTrackerPitch = a; }

  void CalcNormalization();
  
  void AddMagData(filedata_t d){fMagData.push_back(d);}

  /**
   * \fn FindExternalMacros
   * \brief Finds and reads in external macros called from current macro.
   */
  void FindExternalMacros(G4SBSTextFile macro);

  void SetMacroFile(const char *fn){
    fMacro = G4SBSTextFile(fn);
    FindExternalMacros(fMacro);
  }
  void SetPreInitMacroFile(const char *fn){
    fPreInitMacro = G4SBSTextFile(fn);
    FindExternalMacros(fPreInitMacro);
  }

  void Print(Option_t *) const;
  void Print() const;

  TTimeStamp fRunTime;

  long int  fNthrown;
  long int fNtries;
  unsigned int  fSeed;
  double fBeamE; //GeV
  double fNormalization; //Normalization constant to convert observed counts to a rate. This accounts for efficiency of Monte Carlo generation, phase space volume, luminosity, etc
  double fGenVol; //Phase space generation volume; process-dependent;
  double fLuminosity; //Hz/cm^2
  double fMaxWeight; //Maximum calculated event weight (cross section), in appropriate units.
  double fBBtheta; //rad
  double fBBdist;  //m
  double fSBStheta; //rad
  double fSBSdist; //m
  double fHCALdist; //m
  double fHCALvoff; //m
  double fRICHdist; //m
  double fSBSTrackerdist; //m
  double fSBSTrackerPitch; //rad
  
  char fExpType[__RUNSTR_LEN];
  char fGenName[__RUNSTR_LEN];
  char fGitInfo[__GITMAXINFO_SIZE];

  char fHostName[__RUNSTR_LEN];
  char fRunPath[__RUNSTR_LEN];

  G4SBSTextFile              fMacro;
  G4SBSTextFile   fPreInitMacro; //PreInit Commands
  std::vector<G4SBSTextFile> fExternalMacros; ///< External macros called by
                                              ///< the pre-init or post macro

  std::vector<filedata_t> fMagData;

  ClassDef(G4SBSRunData, 1);
};

#endif//__G4SBSRUNDATA_HH
