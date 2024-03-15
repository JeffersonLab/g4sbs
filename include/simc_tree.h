//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Mar 14 14:52:38 2024 by ROOT version 6.26/10
// from TTree h10/h10
// found on file: 0p733sf_sbs7_sbs85p_simc_deep.root
//////////////////////////////////////////////////////////

#ifndef simc_tree_h
#define simc_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class simc_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         hsdelta;
   Float_t         hsyptar;
   Float_t         hsxptar;
   Float_t         hsytar;
   Float_t         hsxfp;
   Float_t         hsxpfp;
   Float_t         hsyfp;
   Float_t         hsypfp;
   Float_t         hsdeltai;
   Float_t         hsyptari;
   Float_t         hsxptari;
   Float_t         hsytari;
   Float_t         ssdelta;
   Float_t         ssyptar;
   Float_t         ssxptar;
   Float_t         ssytar;
   Float_t         ssxfp;
   Float_t         ssxpfp;
   Float_t         ssyfp;
   Float_t         ssypfp;
   Float_t         ssdeltai;
   Float_t         ssyptari;
   Float_t         ssxptari;
   Float_t         ssytari;
   Float_t         q;
   Float_t         nu;
   Float_t         Q2;
   Float_t         W;
   Float_t         epsilon;
   Float_t         Em;
   Float_t         Pm;
   Float_t         thetapq;
   Float_t         phipq;
   Float_t         corrsing;
   Float_t         Pmx;
   Float_t         Pmy;
   Float_t         Pmz;
   Float_t         PmPar;
   Float_t         PmPer;
   Float_t         PmOop;
   Float_t         fry;
   Float_t         radphot;
   Float_t         sigcc;
   Float_t         Weight;
   Float_t         p_e;
   Float_t         ux_e;
   Float_t         uy_e;
   Float_t         uz_e;
   Float_t         p_p;
   Float_t         ux_p;
   Float_t         uy_p;
   Float_t         uz_p;
   Float_t         th_e;
   Float_t         ph_e;
   Float_t         th_p;
   Float_t         ph_p;
   Float_t         vxi;
   Float_t         vyi;
   Float_t         vzi;
   Float_t         ebeam;
   Float_t         veE;
   Float_t         vetheta;
   Float_t         vQ2;
   Float_t         vnu;

   // List of branches
   TBranch        *b_hsdelta;   //!
   TBranch        *b_hsyptar;   //!
   TBranch        *b_hsxptar;   //!
   TBranch        *b_hsytar;   //!
   TBranch        *b_hsxfp;   //!
   TBranch        *b_hsxpfp;   //!
   TBranch        *b_hsyfp;   //!
   TBranch        *b_hsypfp;   //!
   TBranch        *b_hsdeltai;   //!
   TBranch        *b_hsyptari;   //!
   TBranch        *b_hsxptari;   //!
   TBranch        *b_hsytari;   //!
   TBranch        *b_ssdelta;   //!
   TBranch        *b_ssyptar;   //!
   TBranch        *b_ssxptar;   //!
   TBranch        *b_ssytar;   //!
   TBranch        *b_ssxfp;   //!
   TBranch        *b_ssxpfp;   //!
   TBranch        *b_ssyfp;   //!
   TBranch        *b_ssypfp;   //!
   TBranch        *b_ssdeltai;   //!
   TBranch        *b_ssyptari;   //!
   TBranch        *b_ssxptari;   //!
   TBranch        *b_ssytari;   //!
   TBranch        *b_q;   //!
   TBranch        *b_nu;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_W;   //!
   TBranch        *b_epsilon;   //!
   TBranch        *b_Em;   //!
   TBranch        *b_Pm;   //!
   TBranch        *b_thetapq;   //!
   TBranch        *b_phipq;   //!
   TBranch        *b_corrsing;   //!
   TBranch        *b_Pmx;   //!
   TBranch        *b_Pmy;   //!
   TBranch        *b_Pmz;   //!
   TBranch        *b_PmPar;   //!
   TBranch        *b_PmPer;   //!
   TBranch        *b_PmOop;   //!
   TBranch        *b_fry;   //!
   TBranch        *b_radphot;   //!
   TBranch        *b_sigcc;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_p_e;   //!
   TBranch        *b_ux_e;   //!
   TBranch        *b_uy_e;   //!
   TBranch        *b_uz_e;   //!
   TBranch        *b_p_p;   //!
   TBranch        *b_ux_p;   //!
   TBranch        *b_uy_p;   //!
   TBranch        *b_uz_p;   //!
   TBranch        *b_th_e;   //!
   TBranch        *b_ph_e;   //!
   TBranch        *b_th_p;   //!
   TBranch        *b_ph_p;   //!
   TBranch        *b_vxi;   //!
   TBranch        *b_vyi;   //!
   TBranch        *b_vzi;   //!
   TBranch        *b_ebeam;   //!
   TBranch        *b_veE;   //!
   TBranch        *b_vetheta;   //!
   TBranch        *b_vQ2;   //!
   TBranch        *b_vnu;   //!

   simc_tree(TTree *tree=0);
   virtual ~simc_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef simc_tree_cxx
simc_tree::simc_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("0p733sf_sbs7_sbs85p_simc_deep.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("0p733sf_sbs7_sbs85p_simc_deep.root");
      }
      f->GetObject("h10",tree);

   }
   Init(tree);
}

simc_tree::~simc_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t simc_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t simc_tree::LoadTree(Long64_t entry)
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

void simc_tree::Init(TTree *tree)
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

   fChain->SetBranchAddress("hsdelta", &hsdelta, &b_hsdelta);
   fChain->SetBranchAddress("hsyptar", &hsyptar, &b_hsyptar);
   fChain->SetBranchAddress("hsxptar", &hsxptar, &b_hsxptar);
   fChain->SetBranchAddress("hsytar", &hsytar, &b_hsytar);
   fChain->SetBranchAddress("hsxfp", &hsxfp, &b_hsxfp);
   fChain->SetBranchAddress("hsxpfp", &hsxpfp, &b_hsxpfp);
   fChain->SetBranchAddress("hsyfp", &hsyfp, &b_hsyfp);
   fChain->SetBranchAddress("hsypfp", &hsypfp, &b_hsypfp);
   fChain->SetBranchAddress("hsdeltai", &hsdeltai, &b_hsdeltai);
   fChain->SetBranchAddress("hsyptari", &hsyptari, &b_hsyptari);
   fChain->SetBranchAddress("hsxptari", &hsxptari, &b_hsxptari);
   fChain->SetBranchAddress("hsytari", &hsytari, &b_hsytari);
   fChain->SetBranchAddress("ssdelta", &ssdelta, &b_ssdelta);
   fChain->SetBranchAddress("ssyptar", &ssyptar, &b_ssyptar);
   fChain->SetBranchAddress("ssxptar", &ssxptar, &b_ssxptar);
   fChain->SetBranchAddress("ssytar", &ssytar, &b_ssytar);
   fChain->SetBranchAddress("ssxfp", &ssxfp, &b_ssxfp);
   fChain->SetBranchAddress("ssxpfp", &ssxpfp, &b_ssxpfp);
   fChain->SetBranchAddress("ssyfp", &ssyfp, &b_ssyfp);
   fChain->SetBranchAddress("ssypfp", &ssypfp, &b_ssypfp);
   fChain->SetBranchAddress("ssdeltai", &ssdeltai, &b_ssdeltai);
   fChain->SetBranchAddress("ssyptari", &ssyptari, &b_ssyptari);
   fChain->SetBranchAddress("ssxptari", &ssxptari, &b_ssxptari);
   fChain->SetBranchAddress("ssytari", &ssytari, &b_ssytari);
   fChain->SetBranchAddress("q", &q, &b_q);
   fChain->SetBranchAddress("nu", &nu, &b_nu);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("W", &W, &b_W);
   fChain->SetBranchAddress("epsilon", &epsilon, &b_epsilon);
   fChain->SetBranchAddress("Em", &Em, &b_Em);
   fChain->SetBranchAddress("Pm", &Pm, &b_Pm);
   fChain->SetBranchAddress("thetapq", &thetapq, &b_thetapq);
   fChain->SetBranchAddress("phipq", &phipq, &b_phipq);
   fChain->SetBranchAddress("corrsing", &corrsing, &b_corrsing);
   fChain->SetBranchAddress("Pmx", &Pmx, &b_Pmx);
   fChain->SetBranchAddress("Pmy", &Pmy, &b_Pmy);
   fChain->SetBranchAddress("Pmz", &Pmz, &b_Pmz);
   fChain->SetBranchAddress("PmPar", &PmPar, &b_PmPar);
   fChain->SetBranchAddress("PmPer", &PmPer, &b_PmPer);
   fChain->SetBranchAddress("PmOop", &PmOop, &b_PmOop);
   fChain->SetBranchAddress("fry", &fry, &b_fry);
   fChain->SetBranchAddress("radphot", &radphot, &b_radphot);
   fChain->SetBranchAddress("sigcc", &sigcc, &b_sigcc);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("p_e", &p_e, &b_p_e);
   fChain->SetBranchAddress("ux_e", &ux_e, &b_ux_e);
   fChain->SetBranchAddress("uy_e", &uy_e, &b_uy_e);
   fChain->SetBranchAddress("uz_e", &uz_e, &b_uz_e);
   fChain->SetBranchAddress("p_p", &p_p, &b_p_p);
   fChain->SetBranchAddress("ux_p", &ux_p, &b_ux_p);
   fChain->SetBranchAddress("uy_p", &uy_p, &b_uy_p);
   fChain->SetBranchAddress("uz_p", &uz_p, &b_uz_p);
   fChain->SetBranchAddress("th_e", &th_e, &b_th_e);
   fChain->SetBranchAddress("ph_e", &ph_e, &b_ph_e);
   fChain->SetBranchAddress("th_p", &th_p, &b_th_p);
   fChain->SetBranchAddress("ph_p", &ph_p, &b_ph_p);
   fChain->SetBranchAddress("vxi", &vxi, &b_vxi);
   fChain->SetBranchAddress("vyi", &vyi, &b_vyi);
   fChain->SetBranchAddress("vzi", &vzi, &b_vzi);
   fChain->SetBranchAddress("ebeam", &ebeam, &b_ebeam);
   fChain->SetBranchAddress("veE", &veE, &b_veE);
   fChain->SetBranchAddress("vetheta", &vetheta, &b_vetheta);
   fChain->SetBranchAddress("vQ2", &vQ2, &b_vQ2);
   fChain->SetBranchAddress("vnu", &vnu, &b_vnu);
   Notify();
}

Bool_t simc_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void simc_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t simc_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef simc_tree_cxx
