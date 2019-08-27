#ifndef G4SBSIO_H
#define G4SBSIO_H

#include "TROOT.h"
#include "TObject.h"
#include "THashTable.h"
#include "TClonesArray.h"
#include "G4Run.hh"
#include "G4SBSRICHoutput.hh"
#include "G4SBSECaloutput.hh"
#include "G4SBSTrackerOutput.hh"
#include "G4SBSCALoutput.hh"
#include "G4SBSGEMoutput.hh"
#include "G4SBSmTPCoutput.hh"
#include "G4SBSSDTrackOutput.hh"
#include "G4SBSDetectorConstruction.hh"
#include "G4SBSPythiaOutput.hh"
#include "G4SBSAcquMCOutput.hh"

class TFile;
class TTree;
class G4SBSGlobalField;


#define MAXHITDATA 2000

//These aren't really "event"-level quantities, as they are constants describing the setup, and should be stored in the "rundata" object.
typedef struct {
  Double_t thbb, thsbs, dbb, dsbs, dhcal,voffhcal, hoffhcal, drich, dsbstrkr, sbstrkrpitch, dlac, vofflac, hofflac, Ebeam;
} gen_t;


//"count", "rate", "sigma" are all redundant, should really only store one to the tree. 
//Also "solang" is a "run" level quantity, shouldn't really be written to the tree every event, but whatever.
typedef struct {
  // Double_t count, rate, solang, sigma, W2, xbj, Q2, th, ph;
  Double_t count, rate, solang, sigma, W2, xbj, Q2, th, ph, sigmaDIS, sigmaTDIS; // TDIS
  Double_t Aperp, Apar;
  Double_t Pt, Pl;
  Double_t vx, vy, vz;
  Double_t ep, np;
  Double_t p1p, p2p, pip; // TDIS
  Double_t epx, epy, epz;
  Double_t npx, npy, npz;
  Double_t p1px, p1py, p1pz; // TDIS
  Double_t p2px, p2py, p2pz; // TDIS
  Double_t pipx, pipy, pipz; // TDIS
  Double_t nth, nph;
  Double_t p1th, p1ph, p2th, p2ph, pith, piph; // TDIS
  Double_t pmperp, pmpar, pmparsm;
  Double_t z, phperp, phih, MX;
  Double_t Sx, Sy, Sz; //polarization: only meaningful for gun generator!
  Double_t xpi, tpi, xa,pt, nu, ya, y, f2p, f2pi, ypi; // TDIS
  Int_t nucl, fnucl;
  Int_t hadr;
  Int_t earmaccept, harmaccept;
} ev_t;

typedef struct {
  Double_t x, y, xp, yp;
  Double_t tx, ty, txp, typ;
  Int_t hcal, bb, gemtr;
  Double_t hcx, hcy, bcx, bcy;
  Double_t hct, hctex;
  Double_t hclx, hcly, hclz, hcdang;
} tr_t;

typedef struct {
  Int_t ndata;
  Int_t gid[MAXHITDATA];
  Int_t trkrid[MAXHITDATA];
  Double_t x[MAXHITDATA], y[MAXHITDATA], z[MAXHITDATA], t[MAXHITDATA];
  Double_t dx[MAXHITDATA], dy[MAXHITDATA];
  Double_t tx[MAXHITDATA], ty[MAXHITDATA];
  Double_t txp[MAXHITDATA], typ[MAXHITDATA];
  Int_t trid[MAXHITDATA], mid[MAXHITDATA], pid[MAXHITDATA];
  Double_t vx[MAXHITDATA], vy[MAXHITDATA], vz[MAXHITDATA];
  Double_t p[MAXHITDATA], edep[MAXHITDATA];
} hit_t;

typedef struct {
  Int_t hcndata, bcndata;
  Double_t bcx[MAXHITDATA], bcy[MAXHITDATA], bcz[MAXHITDATA], bce[MAXHITDATA], bct[MAXHITDATA];
  Double_t hcx[MAXHITDATA], hcy[MAXHITDATA], hcz[MAXHITDATA], hce[MAXHITDATA], hct[MAXHITDATA];

  Double_t bcvx[MAXHITDATA], bcvy[MAXHITDATA], bcvz[MAXHITDATA];
  Double_t hcvx[MAXHITDATA], hcvy[MAXHITDATA], hcvz[MAXHITDATA];
  
  Int_t hctrid[MAXHITDATA], hcmid[MAXHITDATA], hcpid[MAXHITDATA];
  Int_t bctrid[MAXHITDATA], bcmid[MAXHITDATA], bcpid[MAXHITDATA];

  Int_t row[MAXHITDATA], col[MAXHITDATA], xcell[MAXHITDATA], ycell[MAXHITDATA];

} cal_t;

class G4SBSIO {
public:
  G4SBSIO();
  ~G4SBSIO();
  
  void SetFilename(const char *fn){strcpy(fFilename, fn);}
  //void SetTrackData(tr_t td){ trdata = td; }
  //void SetCalData(cal_t cd){ caldata = cd; }
  void SetEventData(ev_t ed){ evdata = ed; }
  //void SetHitData(hit_t ht){ hitdata = ht; }
  //void SetRICHData( G4SBSRICHoutput rd ) { richdata = rd; }
  //void SetTrackData( G4SBSTrackerOutput td ){ trackdata = td; }
  //void SetGEMData( G4SBSGEMoutput gd ){ GEMdata = gd; }
  //void 

  void SetGEMData( G4String, G4SBSGEMoutput );
  void SetTrackData( G4String, G4SBSTrackerOutput );
  void SetCalData( G4String, G4SBSCALoutput );
  void SetRICHData( G4String, G4SBSRICHoutput );
  void SetECalData( G4String, G4SBSECaloutput );
  void SetmTPCData( G4String, G4SBSmTPCoutput );
  void SetSDtrackData( G4String, G4SBSSDTrackOutput );
  inline void SetAllSDtrackData( G4SBSSDTrackOutput sd ){ allsdtrackdata = sd; } 

  //inline G4SBSSDTrackOutput GetSDtrackData( G4String sdname ){ return sdtrackdata[sdname]; }
  
  //void SetECalData( G4SBSECaloutput ed ){ ecaldata = ed; }

  
  void FillTree();
  void WriteTree();
  
  void SetBeamE(double E){ gendata.Ebeam = E/CLHEP::GeV; }
  void SetBigBiteTheta(double th){ gendata.thbb = th; }
  void SetBigBiteDist(double d){ gendata.dbb = d/CLHEP::m; }
  void SetSBSTheta(double th){ gendata.thsbs = th; }
  void SetHcalDist(double d){ gendata.dhcal = d/CLHEP::m; }
  void SetHcalVOffset(double d){ gendata.voffhcal = d/CLHEP::m; }
  void SetHcalHOffset(double d){ gendata.hoffhcal = d/CLHEP::m; }
  void SetLACDist( double d){ gendata.dlac = d/CLHEP::m; }
  void SetLACVOffset( double d ){ gendata.vofflac = d/CLHEP::m; }
  void SetLACHOffset( double d ){ gendata.hofflac = d/CLHEP::m; }
  void SetSBSDist(double d){ gendata.dsbs = d/CLHEP::m; }
  void SetRICHDist(double d){ gendata.drich = d/CLHEP::m; }
  void SetSBStrkrDist(double d){ gendata.dsbstrkr = d/CLHEP::m; }
  void SetSBStrkrPitch(double a){ gendata.sbstrkrpitch = a; } //radians
  
  void SetGlobalField(G4SBSGlobalField *gf){fGlobalField = gf; }
  
  ev_t GetEventData(){ return evdata; }
  gen_t GetGenData(){ return gendata; }
  
  void InitializeTree();
  void BranchGEM(G4String s);
  void BranchCAL(G4String s);
  void BranchRICH(G4String s);
  //void BranchTracker(G4String s);
  void BranchECAL(G4String s);
  void BranchmTPC(G4String s);
  void BranchPythia();
  // TDIS
  void BranchAcquMC();
  //void BranchSDTracks(G4String s);
  void BranchSDTracks();
  
  void SetDetCon(G4SBSDetectorConstruction *dc ){ fdetcon = dc; }

  // void SetEarmCALpart_flag( G4bool b ){ EarmCALpart_flag = b; }
  // void SetHarmCALpart_flag( G4bool b ){ HarmCALpart_flag = b; }
  map<G4String,G4bool> KeepPartCALflags;
  map<G4String,G4bool> KeepHistoryflags;
  //map<G4String,G4bool> KeepSDtracks;
  
  
  void SetPythiaOutput( G4SBSPythiaOutput p ){ Primaries = p; }
  void SetUsePythia6( G4bool b ){ fUsePythia = b; }
  void SetExclPythia6( G4bool b ){ fExclPythia = b; }

  // TDIS
  void SetAcquMCOutput( G4SBSAcquMCOutput p ){ AcquMCPrimaries = p; }
  void SetUseAcquMC( G4bool b ){ fUseAcquMC = b; }

  map<G4String,G4int> histogram_index; //map with key = SDname, val = histogram index in TClonesArray

  TClonesArray *Esum_histograms;
  TClonesArray *PulseShape_histograms;

  void UpdateGenDataFromDetCon(); //Check and correct any mismatch between constant parameters defined during geometry construction and default values

  G4SBSDetectorConstruction *GetDetCon() { return fdetcon; }

  //Let the user enable this by a UI command either for individual detectors or all detectors, but
  //turn off by default:
  void SetKeepSDtracks( G4String sdname, G4bool b=true ){ fKeepSDtracks[sdname] = b; }
  void SetKeepAllSDtracks( G4bool b ){ fKeepAllSDtracks=b; }

  void SetUsingCerenkov( G4bool b ){ fUsingCerenkov = b; }
  void SetUsingScintillation( G4bool b ){ fUsingScintillation = b; }

  G4bool GetKeepAllSDtracks() const { return fKeepAllSDtracks; }
  map<G4String,G4bool> GetKeepSDtracks() const { return fKeepSDtracks; }
  
private:
  TFile *fFile;
  TTree *fTree;
 
  G4SBSDetectorConstruction *fdetcon;

  G4bool fKeepAllSDtracks;
  map<G4String,G4bool> fKeepSDtracks;
  
  ev_t evdata;
  gen_t gendata;
  //tr_t trdata;
  // cal_t caldata;
  // hit_t hitdata;

  map<G4String,G4SBSGEMoutput> GEMdata;
  map<G4String,G4SBSCALoutput> CALdata;
  map<G4String,G4SBSRICHoutput> richdata;
  map<G4String,G4SBSTrackerOutput> trackdata;
  map<G4String,G4SBSECaloutput> ecaldata;
  map<G4String,G4SBSmTPCoutput> mTPCdata;
  map<G4String,G4SBSSDTrackOutput> sdtrackdata;

  G4SBSSDTrackOutput allsdtrackdata;
  
  G4bool fUsePythia;
  G4bool fExclPythia;
  G4SBSPythiaOutput Primaries;

  // TDIS
  G4bool fUseAcquMC;
  G4SBSAcquMCOutput AcquMCPrimaries;
  
  G4SBSGlobalField *fGlobalField;
  
  char fFilename[255];

  G4int fNhistograms;
  
  // G4bool EarmCALpart_flag;
  // G4bool HarmCALpart_flag;

  //Flags indicating usage of optical photon simulation: if false, don't create tree branches for optical photon
  //sensitive detectors since they will be empty anyway:
  G4bool fUsingCerenkov;
  G4bool fUsingScintillation;
  
};

#endif//G4SBSIO_H
