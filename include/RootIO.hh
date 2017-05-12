// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class RootIO
// Write energy, time, position and ID collected in a sebsitive detector
// element hit to a root TTree
// 05/05/14 JRMA

#ifndef RootIO_h
#define RootIO_h 1

#include "ArrayHit.hh"
#include "G4HCofThisEvent.hh"

enum{ EMaxNPart = 256 };
enum {EInputFile, E2DInputFile, EOutputFile, EInputNoFile};

class TFile;
class TTree;
class TLorentzVector;
class TH3D;
class TH2D;
class TH1D;

// parameters for a ROOT histogram
struct Hparm{
  G4int nX;
  G4double Xmin;
  G4double Xmax;
  G4int nY;
  G4double Ymin;
  G4double Ymax;
  G4int nZ;
  G4double Zmin;
  G4double Zmax;
};

class RootIO
{
public:
  RootIO(const char*, G4int, const char* = NULL);
  RootIO(G4int, G4int*);
  ~RootIO();
protected:
  TFile* fFileI;    //Root input file
  TFile* fFileO;    //Root output file
  TTree* fTreeI;    //ROOT input tree
  TTree* fTreeO;    //ROOT output tree
  G4float fBeam[5]; //beam branch Px,Py,Pz(all unit),Pt,E
  G4float fdircos[EMaxNPart][3]; //direction cosines of generated particles
  G4float *felab;    //Energy of initial generatd particles
  G4float feleak;    //Energy leaking out of system (NOT CURRENTLY IMPLEMENTED)
  G4float fetot;     //Total energy deposited in all detectors
  G4float fWgt[8];
  G4int *fidpart;    //g3 id number of initial generated particle
  G4int fNhits;      // # of hits
  G4int fNpart;      // # generated particles (not necessarily same # tracked)
  G4int fNbranch;
  G4int fNevent;
  G4int fIevent;
  G4float *fplab;   // momentum of original generated particles
  G4float *fvertex; // Vertex position
  G4float **f4Vector;
  G4int *fPartType;
  G4int fIsInput;
  // Hits in detector elements
  G4int fArrayTot; //total number of constructed tof bars
  G4int fNArr; //Number of hits in Array
  G4int *fArr_i; //hit bar indexes
  G4float *fArr_e;  // hit bar energy deposits
  G4float *fArr_ew; // hit bar weighted energy deposits
  G4float *fArr_t;  // hit bar time
  G4float *fArr_te; // hit pickoff time above energy threshold
  G4float *fArr_x;  // x hit position
  G4float *fArr_y;  // y hit position
  G4float *fArr_z;  // z hit position
  G4int fStepTot;
  G4int fNSt;
  G4int *fStID;
  G4int *fStTrID;
  G4int *fStParID;
  G4int *fStStepID;
  G4float *fStdE;
  G4float *fStPx;
  G4float *fStPy;
  G4float *fStPz;
  G4float *fStPE;
  G4float *fStPreX;
  G4float *fStPreY;
  G4float *fStPreZ;
  G4float *fStPostX;
  G4float *fStPostY;
  G4float *fStPostZ;
  TH3D* fH3[8];
  TH2D* fH2[8];
  TH1D* fH1[8];
  TH2D* fH2s;
  G4int* fNHits;
  G4int* fHitID;
  G4int* fHits;
  //
public:
  void SetFile(TFile*);
  TFile* GetFile(){return fFileO;}
  void SetTree(TTree* t){fTreeO=t;}
  TTree* GetTree(){return fTreeO;}
  void SetOutput(const char*);
  void SetInput(const char*);
  void Set2DInput(const char*, const char*);
  //void WriteTree();
  void WriteHit(G4HCofThisEvent* );
  void WriteGenInput();
  void WriteH3(G4int,G4double,G4double,G4double,G4double=1);
  void WriteH2(G4int,G4double,G4double,G4double=1);
  void WriteH1(G4int,G4double,G4double=1);
  void SetH3(G4int, Hparm*);
  void SetH2(G4int, Hparm*);
  void SetH1(G4int, Hparm*);
  void Close();
  void GetEvent();
  G4int GetNevent(){ return fNevent; }
  G4int GetNpart(){ return fNpart; }
  G4int* GetPartType(){ return fPartType; }
  G4float* GetPos(){ return fvertex; }
  G4float** Get4Vector(){ return f4Vector; }
  TH3D* GetH3(G4int i){ return fH3[i]; }
  TH2D* GetH2(G4int i){ return fH2[i]; }
  TH1D* GetH1(G4int i){ return fH1[i]; }
  void SetWgt(G4double wgt, G4int j){ fWgt[j] = wgt; }
  void Sample2D(G4double&, G4double& );
  void SetHits(G4int* nh, G4int* hid, G4int* hits){
    fNHits = nh; fHitID = hid; fHits = hits;
  }

};



#endif










