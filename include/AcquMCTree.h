//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Dec 13 14:39:48 2018 by ROOT version 6.10/02
// from TTree h1/TMCUserGenerator
// found on file: /scratch/HallA/SBS/RTPC/acqu/data/pLow_100kEv_100MeVMom_40cmTar.root
//////////////////////////////////////////////////////////

#ifndef AcquMCTree_h
#define AcquMCTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class AcquMCTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         X_vtx;
   Float_t         Y_vtx;
   Float_t         Z_vtx;
   Float_t         Px_bm;
   Float_t         Py_bm;
   Float_t         Pz_bm;
   Float_t         Pt_bm;
   Float_t         En_bm;
   Float_t         Px_l0114;
   Float_t         Py_l0114;
   Float_t         Pz_l0114;
   Float_t         Pt_l0114;
   Float_t         En_l0114;

   // List of branches
   TBranch        *b_X_vtx;   //!
   TBranch        *b_Y_vtx;   //!
   TBranch        *b_Z_vtx;   //!
   TBranch        *b_Px_bm;   //!
   TBranch        *b_Py_bm;   //!
   TBranch        *b_Pz_bm;   //!
   TBranch        *b_Pt_bm;   //!
   TBranch        *b_En_bm;   //!
   TBranch        *b_Px_l0114;   //!
   TBranch        *b_Py_l0114;   //!
   TBranch        *b_Pz_l0114;   //!
   TBranch        *b_Pt_l0114;   //!
   TBranch        *b_En_l0114;   //!

   AcquMCTree(TTree *tree=0);
   virtual ~AcquMCTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef AcquMCTree_cxx
AcquMCTree::AcquMCTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/scratch/HallA/SBS/RTPC/acqu/data/pLow_100kEv_100MeVMom_40cmTar.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/scratch/HallA/SBS/RTPC/acqu/data/pLow_100kEv_100MeVMom_40cmTar.root");
      }
      f->GetObject("h1",tree);

   }
   Init(tree);
}

AcquMCTree::~AcquMCTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t AcquMCTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t AcquMCTree::LoadTree(Long64_t entry)
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

void AcquMCTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("X_vtx", &X_vtx, &b_X_vtx);
   fChain->SetBranchAddress("Y_vtx", &Y_vtx, &b_Y_vtx);
   fChain->SetBranchAddress("Z_vtx", &Z_vtx, &b_Z_vtx);
   fChain->SetBranchAddress("Px_bm", &Px_bm, &b_Px_bm);
   fChain->SetBranchAddress("Py_bm", &Py_bm, &b_Py_bm);
   fChain->SetBranchAddress("Pz_bm", &Pz_bm, &b_Pz_bm);
   fChain->SetBranchAddress("Pt_bm", &Pt_bm, &b_Pt_bm);
   fChain->SetBranchAddress("En_bm", &En_bm, &b_En_bm);
   fChain->SetBranchAddress("Px_l0114", &Px_l0114, &b_Px_l0114);
   fChain->SetBranchAddress("Py_l0114", &Py_l0114, &b_Py_l0114);
   fChain->SetBranchAddress("Pz_l0114", &Pz_l0114, &b_Pz_l0114);
   fChain->SetBranchAddress("Pt_l0114", &Pt_l0114, &b_Pt_l0114);
   fChain->SetBranchAddress("En_l0114", &En_l0114, &b_En_l0114);
   Notify();
}

Bool_t AcquMCTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void AcquMCTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t AcquMCTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef AcquMCTree_cxx
