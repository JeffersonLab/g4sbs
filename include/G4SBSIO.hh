#ifndef G4SBSIO_H
#define G4SBSIO_H

#include "TROOT.h"
#include "TObject.h"
#include "G4Run.hh"

class TFile;
class TTree;

#define MAXHITDATA 100

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
    Int_t nucl, fnucl;
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
	Double_t x[MAXHITDATA], y[MAXHITDATA], z[MAXHITDATA];
	Double_t dx[MAXHITDATA], dy[MAXHITDATA];
	Double_t tx[MAXHITDATA], ty[MAXHITDATA];
	Double_t txp[MAXHITDATA], typ[MAXHITDATA];
} hit_t;


class G4SBSIO {
    public:
	G4SBSIO();
	~G4SBSIO();

	void SetFilename(const char *fn){strcpy(fFilename, fn);}
	void SetTrackData(tr_t td){ trdata = td; }
	void SetEventData(ev_t ed){ evdata = ed; }
	void SetHitData(hit_t ht){ hitdata = ht; }
	void FillTree();
	void WriteTree();

	void SetBeamE(double E){ gendata.Ebeam = E/GeV; }
	void SetBigBiteTheta(double th){ gendata.thbb = th; }
	void SetBigBiteDist(double d){ gendata.dbb = d/m; }
	void SetHcalTheta(double th){ gendata.thhcal = th; }
	void SetHcalDist(double d){ gendata.dhcal = d/m; }

	ev_t GetEventData(){ return evdata; }

	void InitializeTree();
    private:
	TFile *fFile;
	TTree *fTree;

	ev_t evdata;
	gen_t gendata;
	tr_t trdata;
	hit_t hitdata;

	char fFilename[255];

};

#endif//G4SBSIO_H
