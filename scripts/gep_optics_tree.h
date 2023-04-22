//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 11 18:31:29 2019 by ROOT version 6.14/06
// from TTree Tout/SBS optics and spin transport for GEP
// found on file: SBS_optics_and_spin_inputfile.root
//////////////////////////////////////////////////////////

#ifndef gep_optics_tree_h
#define gep_optics_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class gep_optics_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        beta;
   Double_t        gamma;
   Double_t        xfp;
   Double_t        yfp;
   Double_t        xpfp;
   Double_t        ypfp;
   Double_t        xfprecon;
   Double_t        yfprecon;
   Double_t        xpfprecon;
   Double_t        ypfprecon;
   Double_t        t;
   Double_t        xtar;
   Double_t        ytar;
   Double_t        xptar;
   Double_t        yptar;
   Double_t        p;
   Double_t        chi;
   Double_t        chiphi;
   Double_t        phitrack;
   Double_t        thetatrack;
   Double_t        psitrack;
   Double_t        Pxtg;
   Double_t        Pytg;
   Double_t        Pztg;
   Double_t        Pxfp;
   Double_t        Pyfp;
   Double_t        Pzfp;
   Double_t        Pxfpgeom;
   Double_t        Pyfpgeom;
   Double_t        Pzfpgeom;
   Double_t        Pxfpdipole;
   Double_t        Pyfpdipole;
   Double_t        Pzfpdipole;

   // List of branches
   TBranch        *b_beta;   //!
   TBranch        *b_gamma;   //!
   TBranch        *b_xfp;   //!
   TBranch        *b_yfp;   //!
   TBranch        *b_xpfp;   //!
   TBranch        *b_ypfp;   //!
   TBranch        *b_xfprecon;   //!
   TBranch        *b_yfprecon;   //!
   TBranch        *b_xpfprecon;   //!
   TBranch        *b_ypfprecon;   //!
   TBranch        *b_t;   //!
   TBranch        *b_xtar;   //!
   TBranch        *b_ytar;   //!
   TBranch        *b_xptar;   //!
   TBranch        *b_yptar;   //!
   TBranch        *b_p;   //!
   TBranch        *b_chi;   //!
   TBranch        *b_chiphi;   //!
   TBranch        *b_phitrack;   //!
   TBranch        *b_thetatrack;   //!
   TBranch        *b_psitrack;   //!
   TBranch        *b_Pxtg;   //!
   TBranch        *b_Pytg;   //!
   TBranch        *b_Pztg;   //!
   TBranch        *b_Pxfp;   //!
   TBranch        *b_Pyfp;   //!
   TBranch        *b_Pzfp;   //!
   TBranch        *b_Pxfpgeom;   //!
   TBranch        *b_Pyfpgeom;   //!
   TBranch        *b_Pzfpgeom;   //!
   TBranch        *b_Pxfpdipole;   //!
   TBranch        *b_Pyfpdipole;   //!
   TBranch        *b_Pzfpdipole;   //!

   gep_optics_tree(TTree *tree=0);
   virtual ~gep_optics_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gep_optics_tree_cxx
gep_optics_tree::gep_optics_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("SBS_optics_and_spin_inputfile.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("SBS_optics_and_spin_inputfile.root");
      }
      f->GetObject("Tout",tree);

   }
   Init(tree);
}

gep_optics_tree::~gep_optics_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gep_optics_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gep_optics_tree::LoadTree(Long64_t entry)
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

void gep_optics_tree::Init(TTree *tree)
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

   fChain->SetBranchAddress("beta", &beta, &b_beta);
   fChain->SetBranchAddress("gamma", &gamma, &b_gamma);
   fChain->SetBranchAddress("xfp", &xfp, &b_xfp);
   fChain->SetBranchAddress("yfp", &yfp, &b_yfp);
   fChain->SetBranchAddress("xpfp", &xpfp, &b_xpfp);
   fChain->SetBranchAddress("ypfp", &ypfp, &b_ypfp);
   fChain->SetBranchAddress("xfprecon", &xfprecon, &b_xfprecon);
   fChain->SetBranchAddress("yfprecon", &yfprecon, &b_yfprecon);
   fChain->SetBranchAddress("xpfprecon", &xpfprecon, &b_xpfprecon);
   fChain->SetBranchAddress("ypfprecon", &ypfprecon, &b_ypfprecon);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("xtar", &xtar, &b_xtar);
   fChain->SetBranchAddress("ytar", &ytar, &b_ytar);
   fChain->SetBranchAddress("xptar", &xptar, &b_xptar);
   fChain->SetBranchAddress("yptar", &yptar, &b_yptar);
   fChain->SetBranchAddress("p", &p, &b_p);
   fChain->SetBranchAddress("chi", &chi, &b_chi);
   fChain->SetBranchAddress("chiphi", &chiphi, &b_chiphi);
   fChain->SetBranchAddress("phitrack", &phitrack, &b_phitrack);
   fChain->SetBranchAddress("thetatrack", &thetatrack, &b_thetatrack);
   fChain->SetBranchAddress("psitrack", &psitrack, &b_psitrack);
   fChain->SetBranchAddress("Pxtg", &Pxtg, &b_Pxtg);
   fChain->SetBranchAddress("Pytg", &Pytg, &b_Pytg);
   fChain->SetBranchAddress("Pztg", &Pztg, &b_Pztg);
   fChain->SetBranchAddress("Pxfp", &Pxfp, &b_Pxfp);
   fChain->SetBranchAddress("Pyfp", &Pyfp, &b_Pyfp);
   fChain->SetBranchAddress("Pzfp", &Pzfp, &b_Pzfp);
   fChain->SetBranchAddress("Pxfpgeom", &Pxfpgeom, &b_Pxfpgeom);
   fChain->SetBranchAddress("Pyfpgeom", &Pyfpgeom, &b_Pyfpgeom);
   fChain->SetBranchAddress("Pzfpgeom", &Pzfpgeom, &b_Pzfpgeom);
   fChain->SetBranchAddress("Pxfpdipole", &Pxfpdipole, &b_Pxfpdipole);
   fChain->SetBranchAddress("Pyfpdipole", &Pyfpdipole, &b_Pyfpdipole);
   fChain->SetBranchAddress("Pzfpdipole", &Pzfpdipole, &b_Pzfpdipole);
   Notify();
}

Bool_t gep_optics_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gep_optics_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gep_optics_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gep_optics_tree_cxx
