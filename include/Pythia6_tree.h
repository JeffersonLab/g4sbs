//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Sep 18 12:44:19 2015 by ROOT version 5.34/32
// from TTree Tout/Pythia6.4 min-bias events
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef Pythia6_tree_h
#define Pythia6_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

using namespace std;
// Fixed size dimensions of array or collections stored in the TTree if any.

class Pythia6_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           Nparticles;
   Float_t         Q2;
   Float_t         xbj;
   Float_t         y;
   Float_t         W2;
   vector<int>     *status;
   vector<int>     *pid;
   vector<int>     *parent;
   vector<int>     *fchild;
   vector<int>     *lchild;
   vector<double>  *px;
   vector<double>  *py;
   vector<double>  *pz;
   vector<double>  *E;
   vector<double>  *M;
   vector<double>  *theta;
   vector<double>  *phi;
   vector<double>  *vx;
   vector<double>  *vy;
   vector<double>  *vz;
   vector<double>  *t;
   vector<double>  *tau;

   // List of branches
   TBranch        *b_Nparticles;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_xbj;   //!
   TBranch        *b_y;   //!
   TBranch        *b_W2;   //!
   TBranch        *b_status;   //!
   TBranch        *b_pid;   //!
   TBranch        *b_parent;   //!
   TBranch        *b_fchild;   //!
   TBranch        *b_lchild;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_E;   //!
   TBranch        *b_M;   //!
   TBranch        *b_theta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_t;   //!
   TBranch        *b_tau;   //!

   Pythia6_tree(TTree *tree=0);
   virtual ~Pythia6_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef Pythia6_tree_cxx
Pythia6_tree::Pythia6_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("test.root");
      }
      f->GetObject("Tout",tree);

   }
   Init(tree);
}

Pythia6_tree::~Pythia6_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t Pythia6_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t Pythia6_tree::LoadTree(Long64_t entry)
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

void Pythia6_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   status = 0;
   pid = 0;
   parent = 0;
   fchild = 0;
   lchild = 0;
   px = 0;
   py = 0;
   pz = 0;
   E = 0;
   M = 0;
   theta = 0;
   phi = 0;
   vx = 0;
   vy = 0;
   vz = 0;
   t = 0;
   tau = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("Nparticles", &Nparticles, &b_Nparticles);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("xbj", &xbj, &b_xbj);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("W2", &W2, &b_W2);
   fChain->SetBranchAddress("status", &status, &b_status);
   fChain->SetBranchAddress("pid", &pid, &b_pid);
   fChain->SetBranchAddress("parent", &parent, &b_parent);
   fChain->SetBranchAddress("fchild", &fchild, &b_fchild);
   fChain->SetBranchAddress("lchild", &lchild, &b_lchild);
   fChain->SetBranchAddress("px", &px, &b_px);
   fChain->SetBranchAddress("py", &py, &b_py);
   fChain->SetBranchAddress("pz", &pz, &b_pz);
   fChain->SetBranchAddress("E", &E, &b_E);
   fChain->SetBranchAddress("M", &M, &b_M);
   fChain->SetBranchAddress("theta", &theta, &b_theta);
   fChain->SetBranchAddress("phi", &phi, &b_phi);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("tau", &tau, &b_tau);
   Notify();
}

Bool_t Pythia6_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void Pythia6_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t Pythia6_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef Pythia6_tree_cxx
