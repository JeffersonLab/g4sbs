#ifndef g4sbs_PbF2_tree_h
#define g4sbs_PbF2_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// ----------------------------
// This class is useful to activate variables from the tree 
// stored in g4sbs output root files. 
// 
// It is more particularly dedicated to unfold files obtained with GEp setup.
// It includes information of CDET, GEp ECal, FT, FPP1&2, HCal.
// For more info check the following link: 
// https://hallaweb.jlab.org/wiki/index.php/Documentation_of_g4sbs#ROOT_Tree_Structure 
// 

// Header file for the classes stored in the TTree if any.
#include <vector>
// Fixed size dimensions of array or collections stored in the TTree if any.

class g4sbs_PbF2_tree {
 public :
  TTree          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  // Declaration of leaf types
   
  // Event variables
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
   
  // (0) Big Bite Shower variables: can be a starting point for PbF2:
  // Just copy-paste this block from /* 0 */ to /* *0* */
  // and change the variables name as you want
  // (ex: Earm_BBSHPbF2_hit_nhits, Earm_BBSH_hit_PMT, etc.)
  // You will also need to copy-paste after line marked with a (1)
  /* 0 */
  Int_t                 Earm_BBSH_hit_nhits;
  std::vector<int>     *Earm_BBSH_hit_PMT;
  std::vector<int>     *Earm_BBSH_hit_row;
  std::vector<int>     *Earm_BBSH_hit_col;
  std::vector<int>     *Earm_BBSH_hit_plane;
  std::vector<double>  *Earm_BBSH_hit_xcell;
  std::vector<double>  *Earm_BBSH_hit_ycell;
  std::vector<double>  *Earm_BBSH_hit_zcell;
  std::vector<double>  *Earm_BBSH_hit_xgcell;
  std::vector<double>  *Earm_BBSH_hit_ygcell;
  std::vector<double>  *Earm_BBSH_hit_zgcell;
  std::vector<int>     *Earm_BBSH_hit_NumPhotoelectrons;
  std::vector<double>  *Earm_BBSH_hit_Time_avg;
  std::vector<double>  *Earm_BBSH_hit_Time_rms;
  std::vector<double>  *Earm_BBSH_hit_Time_min;
  std::vector<double>  *Earm_BBSH_hit_Time_max;
   
  Int_t                 Earm_BBSHTF1_hit_nhits;
  std::vector<int>     *Earm_BBSHTF1_hit_row;
  std::vector<int>     *Earm_BBSHTF1_hit_col;
  std::vector<int>     *Earm_BBSHTF1_hit_cell;
  std::vector<int>     *Earm_BBSHTF1_hit_plane;
  std::vector<double>  *Earm_BBSHTF1_hit_xcell;
  std::vector<double>  *Earm_BBSHTF1_hit_ycell;
  std::vector<double>  *Earm_BBSHTF1_hit_zcell;
  std::vector<double>  *Earm_BBSHTF1_hit_xcellg;
  std::vector<double>  *Earm_BBSHTF1_hit_ycellg;
  std::vector<double>  *Earm_BBSHTF1_hit_zcellg;
  std::vector<double>  *Earm_BBSHTF1_hit_xhit;
  std::vector<double>  *Earm_BBSHTF1_hit_yhit;
  std::vector<double>  *Earm_BBSHTF1_hit_zhit;
  std::vector<double>  *Earm_BBSHTF1_hit_sumedep;
  std::vector<double>  *Earm_BBSHTF1_hit_tavg;
  std::vector<double>  *Earm_BBSHTF1_hit_trms;
  std::vector<double>  *Earm_BBSHTF1_hit_tmin;
  std::vector<double>  *Earm_BBSHTF1_hit_tmax;
  /* *0* */
 //0000000000000000000000000000000000000000000000000000000
  Int_t                 Harm_BBSHPbF2_hit_nhits;
  std::vector<int>     *Harm_BBSHPbF2_hit_row;
  std::vector<int>     *Harm_BBSHPbF2_hit_col;
  std::vector<int>     *Harm_BBSHPbF2_hit_cell;
  std::vector<int>     *Harm_BBSHPbF2_hit_plane;
  std::vector<double>  *Harm_BBSHPbF2_hit_xcell;
  std::vector<double>  *Harm_BBSHPbF2_hit_ycell;
  std::vector<double>  *Harm_BBSHPbF2_hit_zcell;
  std::vector<double>  *Harm_BBSHPbF2_hit_xcellg;
  std::vector<double>  *Harm_BBSHPbF2_hit_ycellg;
  std::vector<double>  *Harm_BBSHPbF2_hit_zcellg;
  std::vector<double>  *Harm_BBSHPbF2_hit_xhit;
  std::vector<double>  *Harm_BBSHPbF2_hit_yhit;
  std::vector<double>  *Harm_BBSHPbF2_hit_zhit;
  std::vector<double>  *Harm_BBSHPbF2_hit_sumedep;
  std::vector<double>  *Harm_BBSHPbF2_hit_tavg;
  std::vector<double>  *Harm_BBSHPbF2_hit_trms;
  std::vector<double>  *Harm_BBSHPbF2_hit_tmin;
  std::vector<double>  *Harm_BBSHPbF2_hit_tmax; 
//0000000000000000000000000000000000000000000000000000000

  // (1) List of branches: branches are waht store the values in the TTree Object.
  TBranch        *b_ev;   //!
  TBranch        *b_gen;   //!
  
  // Copy past this block from  /* 1 */ to /* *1* */ and change also the branch names. 
  // Then go to line (2)
  /* 1 */ 
  TBranch        *b_Earm_BBSH_hit_nhits;   //!
  TBranch        *b_Earm_BBSH_hit_PMT;   //!
  TBranch        *b_Earm_BBSH_hit_row;   //!
  TBranch        *b_Earm_BBSH_hit_col;   //!
  TBranch        *b_Earm_BBSH_hit_plane;   //!
  TBranch        *b_Earm_BBSH_hit_xcell;   //!
  TBranch        *b_Earm_BBSH_hit_ycell;   //!
  TBranch        *b_Earm_BBSH_hit_zcell;   //!
  TBranch        *b_Earm_BBSH_hit_xgcell;   //!
  TBranch        *b_Earm_BBSH_hit_ygcell;   //!
  TBranch        *b_Earm_BBSH_hit_zgcell;   //!
  TBranch        *b_Earm_BBSH_hit_NumPhotoelectrons;   //!
  TBranch        *b_Earm_BBSH_hit_Time_avg;   //!
  TBranch        *b_Earm_BBSH_hit_Time_rms;   //!
  TBranch        *b_Earm_BBSH_hit_Time_min;   //!
  TBranch        *b_Earm_BBSH_hit_Time_max;   //!
   
  TBranch        *b_Earm_BBSHTF1_hit_nhits;   //!
  TBranch        *b_Earm_BBSHTF1_hit_row;   //!
  TBranch        *b_Earm_BBSHTF1_hit_col;   //!
  TBranch        *b_Earm_BBSHTF1_hit_cell;   //!
  TBranch        *b_Earm_BBSHTF1_hit_plane;   //!
  TBranch        *b_Earm_BBSHTF1_hit_xcell;   //!
  TBranch        *b_Earm_BBSHTF1_hit_ycell;   //!
  TBranch        *b_Earm_BBSHTF1_hit_zcell;   //!
  TBranch        *b_Earm_BBSHTF1_hit_xcellg;   //!
  TBranch        *b_Earm_BBSHTF1_hit_ycellg;   //!
  TBranch        *b_Earm_BBSHTF1_hit_zcellg;   //!
  TBranch        *b_Earm_BBSHTF1_hit_xhit;   //!
  TBranch        *b_Earm_BBSHTF1_hit_yhit;   //!
  TBranch        *b_Earm_BBSHTF1_hit_zhit;   //!
  TBranch        *b_Earm_BBSHTF1_hit_sumedep;   //!
  TBranch        *b_Earm_BBSHTF1_hit_tavg;   //!
  TBranch        *b_Earm_BBSHTF1_hit_trms;   //!
  TBranch        *b_Earm_BBSHTF1_hit_tmin;   //!
  TBranch        *b_Earm_BBSHTF1_hit_tmax;   //!
  /* *1* */
   /////branch names of CalPbF2 //00000000000000000000000000000000000000000

  TBranch        *b_Harm_BBSHPbF2_hit_nhits;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_row;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_col;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_cell;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_plane;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_xcell;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_ycell;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_zcell;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_xcellg;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_ycellg;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_zcellg;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_xhit;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_yhit;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_zhit;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_sumedep;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_tavg;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_trms;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_tmin;   //!
  TBranch        *b_Harm_BBSHPbF2_hit_tmax;   //!
//0000000000000000000000000000000000000000000000000000000000000000000000000000000000000


  g4sbs_PbF2_tree(TTree *tree=0);
  
  virtual ~g4sbs_PbF2_tree();
  virtual Int_t    Cut(Long64_t entry);
  virtual Int_t    GetEntry(Long64_t entry);
  virtual Long64_t LoadTree(Long64_t entry);
  virtual void     Init(TTree *tree);
  virtual void     Loop();
  virtual Bool_t   Notify();
  virtual void     Show(Long64_t entry = -1);
};
#endif

#ifdef g4sbs_PbF2_tree_cxx
// g4sbs_PbF2_tree constructor: the tree will be the 
// the boolean is a flag to consider(true) or ignore(false) the ECal_box and HCal_box data
g4sbs_PbF2_tree::g4sbs_PbF2_tree(TTree *tree) : fChain(0) 
		       {
			 // if parameter tree is not specified (or zero), connect the file
			 // used to generate this class and read the Tree.
			 if (tree == 0) {
			   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("gep_spin_transport_Sx.root");
			   if (!f || !f->IsOpen()) {
			     f = new TFile("gep_spin_transport_Sx.root");
			   }
			   f->GetObject("T",tree);
			 }
			 Init(tree);
		       }

//default destructor
g4sbs_PbF2_tree::~g4sbs_PbF2_tree()
{
  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

//overload of the TTree::GetEntry(Long64_t) function
Int_t g4sbs_PbF2_tree::GetEntry(Long64_t entry)
{
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}
Long64_t g4sbs_PbF2_tree::LoadTree(Long64_t entry)
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

void g4sbs_PbF2_tree::Init(TTree *tree)
{
  // The Init() function is called when the selector needs to initialize
  // a new tree or chain. Typically here the branch addresses and branch
  // pointers of the tree will be set.
  // It is normally not necessary to make changes to the generated
  // code, but the routine can be extended by the user if needed.
  // Init() will be called many times when running on PROOF
  // (once per file to be processed).
  
  // (2) here you initialize the "std::vector" variables:
  // not doing so will cause the program to crash.
  // copy paste lines from /* 2 */  to /* *2* */
  // and obviously give them the same name as the variables you gave in the block /* 0 */
  // Then go to line (3)
  // Set object pointer
  /* 2 */ 
  Earm_BBSH_hit_PMT = 0;
  Earm_BBSH_hit_row = 0;
  Earm_BBSH_hit_col = 0;
  Earm_BBSH_hit_plane = 0;
  Earm_BBSH_hit_xcell = 0;
  Earm_BBSH_hit_ycell = 0;
  Earm_BBSH_hit_zcell = 0;
  Earm_BBSH_hit_xgcell = 0;
  Earm_BBSH_hit_ygcell = 0;
  Earm_BBSH_hit_zgcell = 0;
  Earm_BBSH_hit_NumPhotoelectrons = 0;
  Earm_BBSH_hit_Time_avg = 0;
  Earm_BBSH_hit_Time_rms = 0;
  Earm_BBSH_hit_Time_min = 0;
  Earm_BBSH_hit_Time_max = 0;
   
  Earm_BBSHTF1_hit_row = 0;
  Earm_BBSHTF1_hit_col = 0;
  Earm_BBSHTF1_hit_cell = 0;
  Earm_BBSHTF1_hit_plane = 0;
  Earm_BBSHTF1_hit_xcell = 0;
  Earm_BBSHTF1_hit_ycell = 0;
  Earm_BBSHTF1_hit_zcell = 0;
  Earm_BBSHTF1_hit_xcellg = 0;
  Earm_BBSHTF1_hit_ycellg = 0;
  Earm_BBSHTF1_hit_zcellg = 0;
  Earm_BBSHTF1_hit_xhit = 0;
  Earm_BBSHTF1_hit_yhit = 0;
  Earm_BBSHTF1_hit_zhit = 0;
  Earm_BBSHTF1_hit_sumedep = 0;
  Earm_BBSHTF1_hit_tavg = 0;
  Earm_BBSHTF1_hit_trms = 0;
  Earm_BBSHTF1_hit_tmin = 0;
  Earm_BBSHTF1_hit_tmax = 0;
  /* *2* */
//00000000000000000000000000000000000000000000000000000
  Harm_BBSHPbF2_hit_row = 0;
  Harm_BBSHPbF2_hit_col = 0;
  Harm_BBSHPbF2_hit_cell = 0;
  Harm_BBSHPbF2_hit_plane = 0;
  Harm_BBSHPbF2_hit_xcell = 0;
  Harm_BBSHPbF2_hit_ycell = 0;
  Harm_BBSHPbF2_hit_zcell = 0;
  Harm_BBSHPbF2_hit_xcellg = 0;
  Harm_BBSHPbF2_hit_ycellg = 0;
  Harm_BBSHPbF2_hit_zcellg = 0;
  Harm_BBSHPbF2_hit_xhit = 0;
  Harm_BBSHPbF2_hit_yhit = 0;
  Harm_BBSHPbF2_hit_zhit = 0;
  Harm_BBSHPbF2_hit_sumedep = 0;
  Harm_BBSHPbF2_hit_tavg = 0;
  Harm_BBSHPbF2_hit_trms = 0;
  Harm_BBSHPbF2_hit_tmin = 0;
  Harm_BBSHPbF2_hit_tmax = 0;
//00000000000000000000000000000000000000000000000000000
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);

  // (3) here you affect the variables to their adresses, 
  // and to the information which is stored in the TTree: 
  // the character chains need to be spelt exactly as they are spelt in the TTree when you open it with TBrowser
  // Again, copy-paste the block from /* 3 */ to /* *3* */,
  // change the character chains (1st argument) with the variable names you want to read in your output, 
  // the variables (2nd argument) with what you have declared in block /* 0 */ to /* *0* */,
  // and the branches (3nd argument) with what you have declared in block /* 1 */ to /* *1* */,
  fChain->SetBranchAddress("ev", &ev_count, &b_ev);
  fChain->SetBranchAddress("gen", &gen_thbb, &b_gen);
  
  /* 3 */
  fChain->SetBranchAddress("Earm.BBSH.hit.nhits", &Earm_BBSH_hit_nhits, &b_Earm_BBSH_hit_nhits);
  fChain->SetBranchAddress("Earm.BBSH.hit.PMT", &Earm_BBSH_hit_PMT, &b_Earm_BBSH_hit_PMT);
  fChain->SetBranchAddress("Earm.BBSH.hit.row", &Earm_BBSH_hit_row, &b_Earm_BBSH_hit_row);
  fChain->SetBranchAddress("Earm.BBSH.hit.col", &Earm_BBSH_hit_col, &b_Earm_BBSH_hit_col);
  fChain->SetBranchAddress("Earm.BBSH.hit.plane", &Earm_BBSH_hit_plane, &b_Earm_BBSH_hit_plane);
  fChain->SetBranchAddress("Earm.BBSH.hit.xcell", &Earm_BBSH_hit_xcell, &b_Earm_BBSH_hit_xcell);
  fChain->SetBranchAddress("Earm.BBSH.hit.ycell", &Earm_BBSH_hit_ycell, &b_Earm_BBSH_hit_ycell);
  fChain->SetBranchAddress("Earm.BBSH.hit.zcell", &Earm_BBSH_hit_zcell, &b_Earm_BBSH_hit_zcell);
  fChain->SetBranchAddress("Earm.BBSH.hit.xgcell", &Earm_BBSH_hit_xgcell, &b_Earm_BBSH_hit_xgcell);
  fChain->SetBranchAddress("Earm.BBSH.hit.ygcell", &Earm_BBSH_hit_ygcell, &b_Earm_BBSH_hit_ygcell);
  fChain->SetBranchAddress("Earm.BBSH.hit.zgcell", &Earm_BBSH_hit_zgcell, &b_Earm_BBSH_hit_zgcell);
  fChain->SetBranchAddress("Earm.BBSH.hit.NumPhotoelectrons", &Earm_BBSH_hit_NumPhotoelectrons, &b_Earm_BBSH_hit_NumPhotoelectrons);
  fChain->SetBranchAddress("Earm.BBSH.hit.Time_avg", &Earm_BBSH_hit_Time_avg, &b_Earm_BBSH_hit_Time_avg);
  fChain->SetBranchAddress("Earm.BBSH.hit.Time_rms", &Earm_BBSH_hit_Time_rms, &b_Earm_BBSH_hit_Time_rms);
  fChain->SetBranchAddress("Earm.BBSH.hit.Time_min", &Earm_BBSH_hit_Time_min, &b_Earm_BBSH_hit_Time_min);
  fChain->SetBranchAddress("Earm.BBSH.hit.Time_max", &Earm_BBSH_hit_Time_max, &b_Earm_BBSH_hit_Time_max);
  
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.nhits", &Earm_BBSHTF1_hit_nhits, &b_Earm_BBSHTF1_hit_nhits);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.row", &Earm_BBSHTF1_hit_row, &b_Earm_BBSHTF1_hit_row);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.col", &Earm_BBSHTF1_hit_col, &b_Earm_BBSHTF1_hit_col);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.cell", &Earm_BBSHTF1_hit_cell, &b_Earm_BBSHTF1_hit_cell);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.plane", &Earm_BBSHTF1_hit_plane, &b_Earm_BBSHTF1_hit_plane);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.xcell", &Earm_BBSHTF1_hit_xcell, &b_Earm_BBSHTF1_hit_xcell);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.ycell", &Earm_BBSHTF1_hit_ycell, &b_Earm_BBSHTF1_hit_ycell);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.zcell", &Earm_BBSHTF1_hit_zcell, &b_Earm_BBSHTF1_hit_zcell);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.xcellg", &Earm_BBSHTF1_hit_xcellg, &b_Earm_BBSHTF1_hit_xcellg);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.ycellg", &Earm_BBSHTF1_hit_ycellg, &b_Earm_BBSHTF1_hit_ycellg);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.zcellg", &Earm_BBSHTF1_hit_zcellg, &b_Earm_BBSHTF1_hit_zcellg);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.xhit", &Earm_BBSHTF1_hit_xhit, &b_Earm_BBSHTF1_hit_xhit);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.yhit", &Earm_BBSHTF1_hit_yhit, &b_Earm_BBSHTF1_hit_yhit);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.zhit", &Earm_BBSHTF1_hit_zhit, &b_Earm_BBSHTF1_hit_zhit);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.sumedep", &Earm_BBSHTF1_hit_sumedep, &b_Earm_BBSHTF1_hit_sumedep);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.tavg", &Earm_BBSHTF1_hit_tavg, &b_Earm_BBSHTF1_hit_tavg);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.trms", &Earm_BBSHTF1_hit_trms, &b_Earm_BBSHTF1_hit_trms);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.tmin", &Earm_BBSHTF1_hit_tmin, &b_Earm_BBSHTF1_hit_tmin);
  fChain->SetBranchAddress("Earm.BBSHTF1.hit.tmax", &Earm_BBSHTF1_hit_tmax, &b_Earm_BBSHTF1_hit_tmax);
   /* *3* */
//00000000000000000000000000000000000000000000000000000
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.nhits", &Harm_BBSHPbF2_hit_nhits, &b_Harm_BBSHPbF2_hit_nhits);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.row", &Harm_BBSHPbF2_hit_row, &b_Harm_BBSHPbF2_hit_row);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.col", &Harm_BBSHPbF2_hit_col, &b_Harm_BBSHPbF2_hit_col);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.cell", &Harm_BBSHPbF2_hit_cell, &b_Harm_BBSHPbF2_hit_cell);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.plane", &Harm_BBSHPbF2_hit_plane, &b_Harm_BBSHPbF2_hit_plane);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.xcell", &Harm_BBSHPbF2_hit_xcell, &b_Harm_BBSHPbF2_hit_xcell);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.ycell", &Harm_BBSHPbF2_hit_ycell, &b_Harm_BBSHPbF2_hit_ycell);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.zcell", &Harm_BBSHPbF2_hit_zcell, &b_Harm_BBSHPbF2_hit_zcell);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.xcellg", &Harm_BBSHPbF2_hit_xcellg, &b_Harm_BBSHPbF2_hit_xcellg);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.ycellg", &Harm_BBSHPbF2_hit_ycellg, &b_Harm_BBSHPbF2_hit_ycellg);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.zcellg", &Harm_BBSHPbF2_hit_zcellg, &b_Harm_BBSHPbF2_hit_zcellg);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.xhit", &Harm_BBSHPbF2_hit_xhit, &b_Harm_BBSHPbF2_hit_xhit);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.yhit", &Harm_BBSHPbF2_hit_yhit, &b_Harm_BBSHPbF2_hit_yhit);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.zhit", &Harm_BBSHPbF2_hit_zhit, &b_Harm_BBSHPbF2_hit_zhit);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.sumedep", &Harm_BBSHPbF2_hit_sumedep, &b_Harm_BBSHPbF2_hit_sumedep);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.tavg", &Harm_BBSHPbF2_hit_tavg, &b_Harm_BBSHPbF2_hit_tavg);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.trms", &Harm_BBSHPbF2_hit_trms, &b_Harm_BBSHPbF2_hit_trms);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.tmin", &Harm_BBSHPbF2_hit_tmin, &b_Harm_BBSHPbF2_hit_tmin);
  fChain->SetBranchAddress("Harm.BBSHPbF2.hit.tmax", &Harm_BBSHPbF2_hit_tmax, &b_Harm_BBSHPbF2_hit_tmax); 
//00000000000000000000000000000000000000000000000000000

 Notify();
}
 

Bool_t g4sbs_PbF2_tree::Notify()
{
  // The Notify() function is called when a new file is opened. This
  // can be either for a new TTree in a TChain or when when a new TTree
  // is started when using PROOF. It is normally not necessary to make changes
  // to the generated code, but the routine can be extended by the
  // user if needed. The return value is currently not used.

  return kTRUE;
}

void g4sbs_PbF2_tree::Show(Long64_t entry)
{
  // Print contents of entry.
  // If entry is not specified, print current entry
  if (!fChain) return;
  fChain->Show(entry);
}

Int_t g4sbs_PbF2_tree::Cut(Long64_t entry)
{
  // This function may be called from Loop.
  // returns  1 if entry is accepted.
  // returns -1 otherwise.
  return 1;
}

#endif // #ifdef g4sbs_PbF2_tree_cxx


