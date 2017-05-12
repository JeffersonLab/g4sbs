// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class PrimaryGeneratorAction
// Generation of particles
// 20/05/13 JRMA adapted from SBS equivalent, under construction
// 14/02/15 JRMA add event generation fron ROOT histogram

#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "globals.hh"
#include "Randomize.hh"

class G4ParticleGun;
class G4Event;
class PrimaryGeneratorMessenger;
class DetectorConstruction;
class RootIO;

const G4int nfermi = 201;
const G4double PI = 3.141592654;

//Event generator mode
enum { EPGA_g4, EPGA_multi, EPGA_ROOT, EPGA_Ebeam, EPGA_2DSample,
       EPGA_TDISp, EPGA_TDISn};

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction(DetectorConstruction*);    
  ~PrimaryGeneratorAction();
  void GeneratePrimaries(G4Event*);
  void SetUpROOTInput(G4String filename);
  void SetUpEbeam();
  //  void GenEbeam(G4Event*);
  void GenVertex();
  //TH1D* GenDsig(G4double);
  void GenBrem(G4Event*);
  void Gen2DSample(G4Event*);
  void InitTDIS();
  void TrackTDIS(G4Event*);
  void GenTDIS();
  void InitDfermi();
  void GetDfermi();
  G4double F2pi(G4double,G4double,G4double);
  G4double URand(G4double x1, G4double x2){
    G4double r = G4UniformRand();
    return x1 + r*(x2 - x1);
  }
private:
  G4ParticleGun*  fParticleGun;     //pointer to particle gun
  G4ParticleTable* fParticleTable;  //pointer to particle table
  G4ParticleDefinition** fPDef;     //pointer to particle definitions
  PrimaryGeneratorMessenger* fPGMessenger; //messenger of this class 
  DetectorConstruction* fDC;        // detector construction class
  G4int fNpart;                     // Number of particles in ntuple
  G4float * fPos;        // vertex position from ntuple, can't be double!
  G4float ** f4Vector;   // 4 vector components from the ntuple branches
  G4float *fMass;              // Masses of the generated particles 
  G4int *fPartType;            // Array of G3 particle types 
  G4int *fTrackThis;           // Array carrying the index of particles
  G4int fNTracked;             // No. of tracked particles in input file
  G4int fNToBeTcount;          // Counter for setting fTrackThis array
  G4int fNevent;               // event number for the ROOT tree
  G4int fMode;                 // event mode selector
  G4ThreeVector fSrcPos;
  G4ThreeVector fSrcSize;
  G4ThreeVector* fPDir;
  G4ThreeVector fPfermi;
  G4double* fPTmin;
  G4double* fPTmax;
  G4double fMe,fEe,fEg;
  G4double fRe,fThe;
  G4int fIpart;
  G4int fIsWindow;
  RootIO* fRootInput;
  char* fName;
  char* fHname;
  //
  G4double fMn, fMp, fMd;
  G4int fTtype;                          // 0 = p, 1 = n
  G4double fSigDIS, fSigTDIS;
  G4double fPfdis[nfermi];               // 2H fermi momentum
  G4double fPpfermi[nfermi];
  G4double fnu, fQ2, fxa, fya;
  G4double ftpi, fypi, fFpi, fF2N, fxpi;
  G4double fxbj, fMx2, fy, fz;
  G4LorentzVector fPei,fPef,fPp1f,fPp2f,fPnf,fPpif;
  G4LorentzVector* fPtdis[6];
public:
  void SetMode(G4int mode){fMode=mode;}
  void SetParticle(G4ThreeVector part){
    fPDef[fIpart] = fParticleTable->FindParticle((G4int)part.x());
    fPTmin[fIpart] = part.y();
    fPTmax[fIpart] = part.z();
  } 
  void SetPDir(G4ThreeVector dir){
    if( fIpart == 0 ){
      fRe = dir.x();
      fThe = dir.y();
    }
    fPDir[fIpart] = dir.unit(); fIpart++;
  }
  void SetNTracked(G4int n)
  { 
    fNTracked=n;fTrackThis=new G4int[n];
    fPTmin = new G4double[n]; fPTmax = new G4double[n];
    fPDef=new G4ParticleDefinition*[n];
    fPDir = new G4ThreeVector[n];
  }
  void SetTracked(G4int ipart)
  {
    if( ipart >= fNpart ) return;
    if(!fTrackThis) return;
    if(fNToBeTcount >= fNpart) return;
    fTrackThis[fNToBeTcount++]=ipart;
    if(fNToBeTcount > fNTracked)fNTracked++;
  }
  void Set2Dname(G4String name){ 
    char* fn = (char*)name.data();
    fHname = new char(strlen(fn)+1);
    strcpy(fHname,fn);
  }
  void SaveTrack(G4int j){
    G4LorentzVector* p = fPtdis[j];
    f4Vector[j][0] = p->x();
    f4Vector[j][1] = p->y();
    f4Vector[j][2] = p->z();
    f4Vector[j][3] = p->t();
  }
  void SetIsWindow(G4int i){ fIsWindow = i; }
  // void SetSeed(G4int seed){ fRandom->SetSeed(seed); }
  G4int GetMode(){return fMode;}
  G4int *GetTracked(){ return fTrackThis; }
  G4int GetNEvents();
  G4int GetNpart(){return fNpart;}
  G4float* GetVertex(){return fPos;}
  G4int* GetPartType(){return fPartType;}
  RootIO* GetRootIO(){ return fRootInput; }
};


#endif


