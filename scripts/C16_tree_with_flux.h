//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Feb 10 11:58:39 2016 by ROOT version 5.34/32
// from TTree T/Geant4 SBS Simulation
// found on file: g4sbsout.root
//////////////////////////////////////////////////////////

#ifndef C16_tree_with_flux_h
#define C16_tree_with_flux_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class C16_tree_with_flux {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        ev_count;
   Double_t        ev_rate;
   Double_t        ev_solang;
   Double_t        ev_sigma;
   Double_t        ev_W2;
   Double_t        ev_xbj;
   Double_t        ev_Q2;
   Double_t        ev_th;
   Double_t        ev_ph;
   Double_t        ev_Aperp;
   Double_t        ev_Apar;
   Double_t        ev_Pt;
   Double_t        ev_Pl;
   Double_t        ev_vx;
   Double_t        ev_vy;
   Double_t        ev_vz;
   Double_t        ev_ep;
   Double_t        ev_np;
   Double_t        ev_epx;
   Double_t        ev_epy;
   Double_t        ev_epz;
   Double_t        ev_npx;
   Double_t        ev_npy;
   Double_t        ev_npz;
   Double_t        ev_nth;
   Double_t        ev_nph;
   Double_t        ev_pmperp;
   Double_t        ev_pmpar;
   Double_t        ev_pmparsm;
   Double_t        ev_z;
   Double_t        ev_phperp;
   Double_t        ev_phih;
   Double_t        ev_MX2;
   Double_t        ev_Sx;
   Double_t        ev_Sy;
   Double_t        ev_Sz;
   Int_t           ev_nucl;
   Int_t           ev_fnucl;
   Int_t           ev_hadr;
   Int_t           ev_earmaccept;
   Int_t           ev_harmaccept;
   Double_t        gen_thbb;
   Double_t        gen_thsbs;
   Double_t        gen_dbb;
   Double_t        gen_dsbs;
   Double_t        gen_dhcal;
   Double_t        gen_drich;
   Double_t        gen_dsbstrkr;
   Double_t        gen_Ebeam;
   Int_t           Earm_C16_hit_nhits;
   vector<int>     *Earm_C16_hit_PMT;
   vector<int>     *Earm_C16_hit_row;
   vector<int>     *Earm_C16_hit_col;
   vector<int>     *Earm_C16_hit_plane;
   vector<double>  *Earm_C16_hit_xcell;
   vector<double>  *Earm_C16_hit_ycell;
   vector<double>  *Earm_C16_hit_zcell;
   vector<double>  *Earm_C16_hit_xgcell;
   vector<double>  *Earm_C16_hit_ygcell;
   vector<double>  *Earm_C16_hit_zgcell;
   vector<int>     *Earm_C16_hit_NumPhotoelectrons;
   vector<double>  *Earm_C16_hit_Time_avg;
   vector<double>  *Earm_C16_hit_Time_rms;
   vector<double>  *Earm_C16_hit_Time_min;
   vector<double>  *Earm_C16_hit_Time_max;
   Int_t           Earm_C16TF1_hit_nhits;
   vector<int>     *Earm_C16TF1_hit_row;
   vector<int>     *Earm_C16TF1_hit_col;
   vector<int>     *Earm_C16TF1_hit_cell;
   vector<int>     *Earm_C16TF1_hit_plane;
   vector<double>  *Earm_C16TF1_hit_xcell;
   vector<double>  *Earm_C16TF1_hit_ycell;
   vector<double>  *Earm_C16TF1_hit_zcell;
   vector<double>  *Earm_C16TF1_hit_xcellg;
   vector<double>  *Earm_C16TF1_hit_ycellg;
   vector<double>  *Earm_C16TF1_hit_zcellg;
   vector<double>  *Earm_C16TF1_hit_xhit;
   vector<double>  *Earm_C16TF1_hit_yhit;
   vector<double>  *Earm_C16TF1_hit_zhit;
   vector<double>  *Earm_C16TF1_hit_sumedep;
   vector<double>  *Earm_C16TF1_hit_tavg;
   vector<double>  *Earm_C16TF1_hit_trms;
   vector<double>  *Earm_C16TF1_hit_tmin;
   vector<double>  *Earm_C16TF1_hit_tmax;
   Int_t           FLUX_hit_nhits;
   vector<int>     *FLUX_hit_row;
   vector<int>     *FLUX_hit_col;
   vector<int>     *FLUX_hit_cell;
   vector<int>     *FLUX_hit_plane;
   vector<double>  *FLUX_hit_xcell;
   vector<double>  *FLUX_hit_ycell;
   vector<double>  *FLUX_hit_zcell;
   vector<double>  *FLUX_hit_xcellg;
   vector<double>  *FLUX_hit_ycellg;
   vector<double>  *FLUX_hit_zcellg;
   vector<double>  *FLUX_hit_xhit;
   vector<double>  *FLUX_hit_yhit;
   vector<double>  *FLUX_hit_zhit;
   vector<double>  *FLUX_hit_sumedep;
   vector<double>  *FLUX_hit_tavg;
   vector<double>  *FLUX_hit_trms;
   vector<double>  *FLUX_hit_tmin;
   vector<double>  *FLUX_hit_tmax;
   Int_t           FLUX_npart_CAL;
   vector<int>     *FLUX_ihit;
   vector<double>  *FLUX_x;
   vector<double>  *FLUX_y;
   vector<double>  *FLUX_z;
   vector<double>  *FLUX_t;
   vector<double>  *FLUX_E;
   vector<double>  *FLUX_dt;
   vector<double>  *FLUX_L;
   vector<double>  *FLUX_vx;
   vector<double>  *FLUX_vy;
   vector<double>  *FLUX_vz;
   vector<int>     *FLUX_trid;
   vector<int>     *FLUX_mid;
   vector<int>     *FLUX_pid;
   vector<double>  *FLUX_p;
   vector<double>  *FLUX_px;
   vector<double>  *FLUX_py;
   vector<double>  *FLUX_pz;
   vector<double>  *FLUX_edep;

   // List of branches
   TBranch        *b_ev;   //!
   TBranch        *b_gen;   //!
   TBranch        *b_Earm_C16_hit_nhits;   //!
   TBranch        *b_Earm_C16_hit_PMT;   //!
   TBranch        *b_Earm_C16_hit_row;   //!
   TBranch        *b_Earm_C16_hit_col;   //!
   TBranch        *b_Earm_C16_hit_plane;   //!
   TBranch        *b_Earm_C16_hit_xcell;   //!
   TBranch        *b_Earm_C16_hit_ycell;   //!
   TBranch        *b_Earm_C16_hit_zcell;   //!
   TBranch        *b_Earm_C16_hit_xgcell;   //!
   TBranch        *b_Earm_C16_hit_ygcell;   //!
   TBranch        *b_Earm_C16_hit_zgcell;   //!
   TBranch        *b_Earm_C16_hit_NumPhotoelectrons;   //!
   TBranch        *b_Earm_C16_hit_Time_avg;   //!
   TBranch        *b_Earm_C16_hit_Time_rms;   //!
   TBranch        *b_Earm_C16_hit_Time_min;   //!
   TBranch        *b_Earm_C16_hit_Time_max;   //!
   TBranch        *b_Earm_C16TF1_hit_nhits;   //!
   TBranch        *b_Earm_C16TF1_hit_row;   //!
   TBranch        *b_Earm_C16TF1_hit_col;   //!
   TBranch        *b_Earm_C16TF1_hit_cell;   //!
   TBranch        *b_Earm_C16TF1_hit_plane;   //!
   TBranch        *b_Earm_C16TF1_hit_xcell;   //!
   TBranch        *b_Earm_C16TF1_hit_ycell;   //!
   TBranch        *b_Earm_C16TF1_hit_zcell;   //!
   TBranch        *b_Earm_C16TF1_hit_xcellg;   //!
   TBranch        *b_Earm_C16TF1_hit_ycellg;   //!
   TBranch        *b_Earm_C16TF1_hit_zcellg;   //!
   TBranch        *b_Earm_C16TF1_hit_xhit;   //!
   TBranch        *b_Earm_C16TF1_hit_yhit;   //!
   TBranch        *b_Earm_C16TF1_hit_zhit;   //!
   TBranch        *b_Earm_C16TF1_hit_sumedep;   //!
   TBranch        *b_Earm_C16TF1_hit_tavg;   //!
   TBranch        *b_Earm_C16TF1_hit_trms;   //!
   TBranch        *b_Earm_C16TF1_hit_tmin;   //!
   TBranch        *b_Earm_C16TF1_hit_tmax;   //!
   TBranch        *b_FLUX_hit_nhits;   //!
   TBranch        *b_FLUX_hit_row;   //!
   TBranch        *b_FLUX_hit_col;   //!
   TBranch        *b_FLUX_hit_cell;   //!
   TBranch        *b_FLUX_hit_plane;   //!
   TBranch        *b_FLUX_hit_xcell;   //!
   TBranch        *b_FLUX_hit_ycell;   //!
   TBranch        *b_FLUX_hit_zcell;   //!
   TBranch        *b_FLUX_hit_xcellg;   //!
   TBranch        *b_FLUX_hit_ycellg;   //!
   TBranch        *b_FLUX_hit_zcellg;   //!
   TBranch        *b_FLUX_hit_xhit;   //!
   TBranch        *b_FLUX_hit_yhit;   //!
   TBranch        *b_FLUX_hit_zhit;   //!
   TBranch        *b_FLUX_hit_sumedep;   //!
   TBranch        *b_FLUX_hit_tavg;   //!
   TBranch        *b_FLUX_hit_trms;   //!
   TBranch        *b_FLUX_hit_tmin;   //!
   TBranch        *b_FLUX_hit_tmax;   //!
   TBranch        *b_FLUX_npart_CAL;   //!
   TBranch        *b_FLUX_ihit;   //!
   TBranch        *b_FLUX_x;   //!
   TBranch        *b_FLUX_y;   //!
   TBranch        *b_FLUX_z;   //!
   TBranch        *b_FLUX_t;   //!
   TBranch        *b_FLUX_E;   //!
   TBranch        *b_FLUX_dt;   //!
   TBranch        *b_FLUX_L;   //!
   TBranch        *b_FLUX_vx;   //!
   TBranch        *b_FLUX_vy;   //!
   TBranch        *b_FLUX_vz;   //!
   TBranch        *b_FLUX_trid;   //!
   TBranch        *b_FLUX_mid;   //!
   TBranch        *b_FLUX_pid;   //!
   TBranch        *b_FLUX_p;   //!
   TBranch        *b_FLUX_px;   //!
   TBranch        *b_FLUX_py;   //!
   TBranch        *b_FLUX_pz;   //!
   TBranch        *b_FLUX_edep;   //!

   C16_tree_with_flux(TTree *tree=0);
   virtual ~C16_tree_with_flux();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef C16_tree_with_flux_cxx
C16_tree_with_flux::C16_tree_with_flux(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("g4sbsout.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("g4sbsout.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

C16_tree_with_flux::~C16_tree_with_flux()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t C16_tree_with_flux::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t C16_tree_with_flux::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void C16_tree_with_flux::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Earm_C16_hit_PMT = 0;
   Earm_C16_hit_row = 0;
   Earm_C16_hit_col = 0;
   Earm_C16_hit_plane = 0;
   Earm_C16_hit_xcell = 0;
   Earm_C16_hit_ycell = 0;
   Earm_C16_hit_zcell = 0;
   Earm_C16_hit_xgcell = 0;
   Earm_C16_hit_ygcell = 0;
   Earm_C16_hit_zgcell = 0;
   Earm_C16_hit_NumPhotoelectrons = 0;
   Earm_C16_hit_Time_avg = 0;
   Earm_C16_hit_Time_rms = 0;
   Earm_C16_hit_Time_min = 0;
   Earm_C16_hit_Time_max = 0;
   Earm_C16TF1_hit_row = 0;
   Earm_C16TF1_hit_col = 0;
   Earm_C16TF1_hit_cell = 0;
   Earm_C16TF1_hit_plane = 0;
   Earm_C16TF1_hit_xcell = 0;
   Earm_C16TF1_hit_ycell = 0;
   Earm_C16TF1_hit_zcell = 0;
   Earm_C16TF1_hit_xcellg = 0;
   Earm_C16TF1_hit_ycellg = 0;
   Earm_C16TF1_hit_zcellg = 0;
   Earm_C16TF1_hit_xhit = 0;
   Earm_C16TF1_hit_yhit = 0;
   Earm_C16TF1_hit_zhit = 0;
   Earm_C16TF1_hit_sumedep = 0;
   Earm_C16TF1_hit_tavg = 0;
   Earm_C16TF1_hit_trms = 0;
   Earm_C16TF1_hit_tmin = 0;
   Earm_C16TF1_hit_tmax = 0;
   FLUX_hit_row = 0;
   FLUX_hit_col = 0;
   FLUX_hit_cell = 0;
   FLUX_hit_plane = 0;
   FLUX_hit_xcell = 0;
   FLUX_hit_ycell = 0;
   FLUX_hit_zcell = 0;
   FLUX_hit_xcellg = 0;
   FLUX_hit_ycellg = 0;
   FLUX_hit_zcellg = 0;
   FLUX_hit_xhit = 0;
   FLUX_hit_yhit = 0;
   FLUX_hit_zhit = 0;
   FLUX_hit_sumedep = 0;
   FLUX_hit_tavg = 0;
   FLUX_hit_trms = 0;
   FLUX_hit_tmin = 0;
   FLUX_hit_tmax = 0;
   FLUX_ihit = 0;
   FLUX_x = 0;
   FLUX_y = 0;
   FLUX_z = 0;
   FLUX_t = 0;
   FLUX_E = 0;
   FLUX_dt = 0;
   FLUX_L = 0;
   FLUX_vx = 0;
   FLUX_vy = 0;
   FLUX_vz = 0;
   FLUX_trid = 0;
   FLUX_mid = 0;
   FLUX_pid = 0;
   FLUX_p = 0;
   FLUX_px = 0;
   FLUX_py = 0;
   FLUX_pz = 0;
   FLUX_edep = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev", &ev_count, &b_ev);
   fChain->SetBranchAddress("gen", &gen_thbb, &b_gen);
   fChain->SetBranchAddress("Earm.C16.hit.nhits", &Earm_C16_hit_nhits, &b_Earm_C16_hit_nhits);
   fChain->SetBranchAddress("Earm.C16.hit.PMT", &Earm_C16_hit_PMT, &b_Earm_C16_hit_PMT);
   fChain->SetBranchAddress("Earm.C16.hit.row", &Earm_C16_hit_row, &b_Earm_C16_hit_row);
   fChain->SetBranchAddress("Earm.C16.hit.col", &Earm_C16_hit_col, &b_Earm_C16_hit_col);
   fChain->SetBranchAddress("Earm.C16.hit.plane", &Earm_C16_hit_plane, &b_Earm_C16_hit_plane);
   fChain->SetBranchAddress("Earm.C16.hit.xcell", &Earm_C16_hit_xcell, &b_Earm_C16_hit_xcell);
   fChain->SetBranchAddress("Earm.C16.hit.ycell", &Earm_C16_hit_ycell, &b_Earm_C16_hit_ycell);
   fChain->SetBranchAddress("Earm.C16.hit.zcell", &Earm_C16_hit_zcell, &b_Earm_C16_hit_zcell);
   fChain->SetBranchAddress("Earm.C16.hit.xgcell", &Earm_C16_hit_xgcell, &b_Earm_C16_hit_xgcell);
   fChain->SetBranchAddress("Earm.C16.hit.ygcell", &Earm_C16_hit_ygcell, &b_Earm_C16_hit_ygcell);
   fChain->SetBranchAddress("Earm.C16.hit.zgcell", &Earm_C16_hit_zgcell, &b_Earm_C16_hit_zgcell);
   fChain->SetBranchAddress("Earm.C16.hit.NumPhotoelectrons", &Earm_C16_hit_NumPhotoelectrons, &b_Earm_C16_hit_NumPhotoelectrons);
   fChain->SetBranchAddress("Earm.C16.hit.Time_avg", &Earm_C16_hit_Time_avg, &b_Earm_C16_hit_Time_avg);
   fChain->SetBranchAddress("Earm.C16.hit.Time_rms", &Earm_C16_hit_Time_rms, &b_Earm_C16_hit_Time_rms);
   fChain->SetBranchAddress("Earm.C16.hit.Time_min", &Earm_C16_hit_Time_min, &b_Earm_C16_hit_Time_min);
   fChain->SetBranchAddress("Earm.C16.hit.Time_max", &Earm_C16_hit_Time_max, &b_Earm_C16_hit_Time_max);
   fChain->SetBranchAddress("Earm.C16TF1.hit.nhits", &Earm_C16TF1_hit_nhits, &b_Earm_C16TF1_hit_nhits);
   fChain->SetBranchAddress("Earm.C16TF1.hit.row", &Earm_C16TF1_hit_row, &b_Earm_C16TF1_hit_row);
   fChain->SetBranchAddress("Earm.C16TF1.hit.col", &Earm_C16TF1_hit_col, &b_Earm_C16TF1_hit_col);
   fChain->SetBranchAddress("Earm.C16TF1.hit.cell", &Earm_C16TF1_hit_cell, &b_Earm_C16TF1_hit_cell);
   fChain->SetBranchAddress("Earm.C16TF1.hit.plane", &Earm_C16TF1_hit_plane, &b_Earm_C16TF1_hit_plane);
   fChain->SetBranchAddress("Earm.C16TF1.hit.xcell", &Earm_C16TF1_hit_xcell, &b_Earm_C16TF1_hit_xcell);
   fChain->SetBranchAddress("Earm.C16TF1.hit.ycell", &Earm_C16TF1_hit_ycell, &b_Earm_C16TF1_hit_ycell);
   fChain->SetBranchAddress("Earm.C16TF1.hit.zcell", &Earm_C16TF1_hit_zcell, &b_Earm_C16TF1_hit_zcell);
   fChain->SetBranchAddress("Earm.C16TF1.hit.xcellg", &Earm_C16TF1_hit_xcellg, &b_Earm_C16TF1_hit_xcellg);
   fChain->SetBranchAddress("Earm.C16TF1.hit.ycellg", &Earm_C16TF1_hit_ycellg, &b_Earm_C16TF1_hit_ycellg);
   fChain->SetBranchAddress("Earm.C16TF1.hit.zcellg", &Earm_C16TF1_hit_zcellg, &b_Earm_C16TF1_hit_zcellg);
   fChain->SetBranchAddress("Earm.C16TF1.hit.xhit", &Earm_C16TF1_hit_xhit, &b_Earm_C16TF1_hit_xhit);
   fChain->SetBranchAddress("Earm.C16TF1.hit.yhit", &Earm_C16TF1_hit_yhit, &b_Earm_C16TF1_hit_yhit);
   fChain->SetBranchAddress("Earm.C16TF1.hit.zhit", &Earm_C16TF1_hit_zhit, &b_Earm_C16TF1_hit_zhit);
   fChain->SetBranchAddress("Earm.C16TF1.hit.sumedep", &Earm_C16TF1_hit_sumedep, &b_Earm_C16TF1_hit_sumedep);
   fChain->SetBranchAddress("Earm.C16TF1.hit.tavg", &Earm_C16TF1_hit_tavg, &b_Earm_C16TF1_hit_tavg);
   fChain->SetBranchAddress("Earm.C16TF1.hit.trms", &Earm_C16TF1_hit_trms, &b_Earm_C16TF1_hit_trms);
   fChain->SetBranchAddress("Earm.C16TF1.hit.tmin", &Earm_C16TF1_hit_tmin, &b_Earm_C16TF1_hit_tmin);
   fChain->SetBranchAddress("Earm.C16TF1.hit.tmax", &Earm_C16TF1_hit_tmax, &b_Earm_C16TF1_hit_tmax);
   fChain->SetBranchAddress("FLUX.hit.nhits", &FLUX_hit_nhits, &b_FLUX_hit_nhits);
   fChain->SetBranchAddress("FLUX.hit.row", &FLUX_hit_row, &b_FLUX_hit_row);
   fChain->SetBranchAddress("FLUX.hit.col", &FLUX_hit_col, &b_FLUX_hit_col);
   fChain->SetBranchAddress("FLUX.hit.cell", &FLUX_hit_cell, &b_FLUX_hit_cell);
   fChain->SetBranchAddress("FLUX.hit.plane", &FLUX_hit_plane, &b_FLUX_hit_plane);
   fChain->SetBranchAddress("FLUX.hit.xcell", &FLUX_hit_xcell, &b_FLUX_hit_xcell);
   fChain->SetBranchAddress("FLUX.hit.ycell", &FLUX_hit_ycell, &b_FLUX_hit_ycell);
   fChain->SetBranchAddress("FLUX.hit.zcell", &FLUX_hit_zcell, &b_FLUX_hit_zcell);
   fChain->SetBranchAddress("FLUX.hit.xcellg", &FLUX_hit_xcellg, &b_FLUX_hit_xcellg);
   fChain->SetBranchAddress("FLUX.hit.ycellg", &FLUX_hit_ycellg, &b_FLUX_hit_ycellg);
   fChain->SetBranchAddress("FLUX.hit.zcellg", &FLUX_hit_zcellg, &b_FLUX_hit_zcellg);
   fChain->SetBranchAddress("FLUX.hit.xhit", &FLUX_hit_xhit, &b_FLUX_hit_xhit);
   fChain->SetBranchAddress("FLUX.hit.yhit", &FLUX_hit_yhit, &b_FLUX_hit_yhit);
   fChain->SetBranchAddress("FLUX.hit.zhit", &FLUX_hit_zhit, &b_FLUX_hit_zhit);
   fChain->SetBranchAddress("FLUX.hit.sumedep", &FLUX_hit_sumedep, &b_FLUX_hit_sumedep);
   fChain->SetBranchAddress("FLUX.hit.tavg", &FLUX_hit_tavg, &b_FLUX_hit_tavg);
   fChain->SetBranchAddress("FLUX.hit.trms", &FLUX_hit_trms, &b_FLUX_hit_trms);
   fChain->SetBranchAddress("FLUX.hit.tmin", &FLUX_hit_tmin, &b_FLUX_hit_tmin);
   fChain->SetBranchAddress("FLUX.hit.tmax", &FLUX_hit_tmax, &b_FLUX_hit_tmax);
   fChain->SetBranchAddress("FLUX.npart_CAL", &FLUX_npart_CAL, &b_FLUX_npart_CAL);
   fChain->SetBranchAddress("FLUX.ihit", &FLUX_ihit, &b_FLUX_ihit);
   fChain->SetBranchAddress("FLUX.x", &FLUX_x, &b_FLUX_x);
   fChain->SetBranchAddress("FLUX.y", &FLUX_y, &b_FLUX_y);
   fChain->SetBranchAddress("FLUX.z", &FLUX_z, &b_FLUX_z);
   fChain->SetBranchAddress("FLUX.t", &FLUX_t, &b_FLUX_t);
   fChain->SetBranchAddress("FLUX.E", &FLUX_E, &b_FLUX_E);
   fChain->SetBranchAddress("FLUX.dt", &FLUX_dt, &b_FLUX_dt);
   fChain->SetBranchAddress("FLUX.L", &FLUX_L, &b_FLUX_L);
   fChain->SetBranchAddress("FLUX.vx", &FLUX_vx, &b_FLUX_vx);
   fChain->SetBranchAddress("FLUX.vy", &FLUX_vy, &b_FLUX_vy);
   fChain->SetBranchAddress("FLUX.vz", &FLUX_vz, &b_FLUX_vz);
   fChain->SetBranchAddress("FLUX.trid", &FLUX_trid, &b_FLUX_trid);
   fChain->SetBranchAddress("FLUX.mid", &FLUX_mid, &b_FLUX_mid);
   fChain->SetBranchAddress("FLUX.pid", &FLUX_pid, &b_FLUX_pid);
   fChain->SetBranchAddress("FLUX.p", &FLUX_p, &b_FLUX_p);
   fChain->SetBranchAddress("FLUX.px", &FLUX_px, &b_FLUX_px);
   fChain->SetBranchAddress("FLUX.py", &FLUX_py, &b_FLUX_py);
   fChain->SetBranchAddress("FLUX.pz", &FLUX_pz, &b_FLUX_pz);
   fChain->SetBranchAddress("FLUX.edep", &FLUX_edep, &b_FLUX_edep);
   Notify();
}

Bool_t C16_tree_with_flux::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void C16_tree_with_flux::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t C16_tree_with_flux::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef C16_tree_with_flux_cxx
