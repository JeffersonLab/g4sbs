#ifndef G4SBSIO_H
#define G4SBSIO_H

#include "TROOT.h"
#include "TObject.h"
#include "THashTable.h"
#include "G4Run.hh"
#include "G4SystemOfUnits.hh"
#include "G4SBSRICHoutput.hh"
#include "G4SBSTrackerOutput.hh"
#include "G4SBSECaloutput.hh"

class TFile;
class TTree;
class G4SBSGlobalField;


#define MAXHITDATA 2000

typedef struct {
    Double_t thbb, thhcal, dbb, dhcal, Ebeam;
} gen_t;


typedef struct {
    Double_t count, rate, solang, sigma, W2, xbj, Q2, th, ph;
    Double_t Aperp, Apar;
    Double_t vx, vy, vz;
    Double_t ep, np;
    Double_t epx, epy, epz;
    Double_t npx, npy, npz;
    Double_t nth, nph;
    Double_t pmperp, pmpar, pmparsm;
  Double_t z, phperp, phih, MX;
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
} cal_t;

class G4SBSIO {
public:
  G4SBSIO();
  ~G4SBSIO();
  
  void SetFilename(const char *fn){strcpy(fFilename, fn);}
  //void SetTrackData(tr_t td){ trdata = td; }
  void SetCalData(cal_t cd){ caldata = cd; }
  void SetEventData(ev_t ed){ evdata = ed; }
  void SetHitData(hit_t ht){ hitdata = ht; }
  void SetRICHData( G4SBSRICHoutput rd ) { richdata = rd; }
  void SetTrackData( G4SBSTrackerOutput td ){ trackdata = td; }
  void SetECalData( G4SBSECaloutput ed ){ ecaldata = ed; }
  
  void FillTree();
  void WriteTree();
  
  void SetBeamE(double E){ gendata.Ebeam = E/GeV; }
  void SetBigBiteTheta(double th){ gendata.thbb = th; }
  void SetBigBiteDist(double d){ gendata.dbb = d/m; }
  void SetHcalTheta(double th){ gendata.thhcal = th; }
  void SetHcalDist(double d){ gendata.dhcal = d/m; }
  
  void SetGlobalField(G4SBSGlobalField *gf){fGlobalField = gf; }
  
  ev_t GetEventData(){ return evdata; }
  
  void InitializeTree();
private:
  TFile *fFile;
  TTree *fTree;
  
  ev_t evdata;
  gen_t gendata;
  //tr_t trdata;
  cal_t caldata;
  hit_t hitdata;
  
  G4SBSRICHoutput richdata;
  G4SBSTrackerOutput trackdata;
  G4SBSECaloutput ecaldata;

  G4SBSGlobalField *fGlobalField;
  
  char fFilename[255];
  
};

#endif//G4SBSIO_H
