//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 26 15:48:20 2025 by ROOT version 6.30/04
// from TTree h10/h10
// found on file: /w/halla-scshelf2102/tdis/rmontgom/nDVCS/simc_out/nps_sbs_sidis_pi0_test_job_1.root
//////////////////////////////////////////////////////////

#ifndef simc_pi0_tree_h
#define simc_pi0_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class simc_pi0_tree {
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
   Float_t         missmass;
   Float_t         ppi;
   Float_t         t;
   Float_t         fry;
   Float_t         radphot;
   Float_t         siglab;
   Float_t         sigcent;
   Float_t         Weight;
   Float_t         decdist;
   Float_t         Mhadron;
   Float_t         z;
   Float_t         zi;
   Float_t         pt2;
   Float_t         pt2i;
   Float_t         xbj;
   Float_t         xbji;
   Float_t         thqi;
   Float_t         sighad;
   Float_t         jacobian;
   Float_t         centjac;
   Float_t         pfermi;
   Float_t         xfermi;
   Float_t         phipqi;
   Float_t         xcal_gamma1;
   Float_t         ycal_gamma1;
   Float_t         Egamma1;
   Float_t         Pgamma1x;
   Float_t         Pgamma1y;
   Float_t         Pgamma1z;
   Float_t         xcal_gamma2;
   Float_t         ycal_gamma2;
   Float_t         Egamma2;
   Float_t         Pgamma2x;
   Float_t         Pgamma2y;
   Float_t         Pgamma2z;
   Float_t         p_e;
   Float_t         ux_e;
   Float_t         uy_e;
   Float_t         uz_e;

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
   TBranch        *b_missmass;   //!
   TBranch        *b_ppi;   //!
   TBranch        *b_t;   //!
   TBranch        *b_fry;   //!
   TBranch        *b_radphot;   //!
   TBranch        *b_siglab;   //!
   TBranch        *b_sigcent;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_decdist;   //!
   TBranch        *b_Mhadron;   //!
   TBranch        *b_z;   //!
   TBranch        *b_zi;   //!
   TBranch        *b_pt2;   //!
   TBranch        *b_pt2i;   //!
   TBranch        *b_xbj;   //!
   TBranch        *b_xbji;   //!
   TBranch        *b_thqi;   //!
   TBranch        *b_sighad;   //!
   TBranch        *b_jacobian;   //!
   TBranch        *b_centjac;   //!
   TBranch        *b_pfermi;   //!
   TBranch        *b_xfermi;   //!
   TBranch        *b_phipqi;   //!
   TBranch        *b_xcal_gamma1;   //!
   TBranch        *b_ycal_gamma1;   //!
   TBranch        *b_Egamma1;   //!
   TBranch        *b_Pgamma1x;   //!
   TBranch        *b_Pgamma1y;   //!
   TBranch        *b_Pgamma1z;   //!
   TBranch        *b_xcal_gamma2;   //!
   TBranch        *b_ycal_gamma2;   //!
   TBranch        *b_Egamma2;   //!
   TBranch        *b_Pgamma2x;   //!
   TBranch        *b_Pgamma2y;   //!
   TBranch        *b_Pgamma2z;   //!
   TBranch        *b_p_e;   //!
   TBranch        *b_ux_e;   //!
   TBranch        *b_uy_e;   //!
   TBranch        *b_uz_e;   //!

   simc_pi0_tree(TTree *tree=0);
   virtual ~simc_pi0_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef simc_pi0_tree_cxx
simc_pi0_tree::simc_pi0_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/w/halla-scshelf2102/tdis/rmontgom/nDVCS/simc_out/nps_sbs_sidis_pi0_test_job_1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/w/halla-scshelf2102/tdis/rmontgom/nDVCS/simc_out/nps_sbs_sidis_pi0_test_job_1.root");
      }
      f->GetObject("h10",tree);

   }
   Init(tree);
}

simc_pi0_tree::~simc_pi0_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t simc_pi0_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t simc_pi0_tree::LoadTree(Long64_t entry)
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

void simc_pi0_tree::Init(TTree *tree)
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
   fChain->SetBranchAddress("missmass", &missmass, &b_missmass);
   fChain->SetBranchAddress("ppi", &ppi, &b_ppi);
   fChain->SetBranchAddress("t", &t, &b_t);
   fChain->SetBranchAddress("fry", &fry, &b_fry);
   fChain->SetBranchAddress("radphot", &radphot, &b_radphot);
   fChain->SetBranchAddress("siglab", &siglab, &b_siglab);
   fChain->SetBranchAddress("sigcent", &sigcent, &b_sigcent);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("decdist", &decdist, &b_decdist);
   fChain->SetBranchAddress("Mhadron", &Mhadron, &b_Mhadron);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("zi", &zi, &b_zi);
   fChain->SetBranchAddress("pt2", &pt2, &b_pt2);
   fChain->SetBranchAddress("pt2i", &pt2i, &b_pt2i);
   fChain->SetBranchAddress("xbj", &xbj, &b_xbj);
   fChain->SetBranchAddress("xbji", &xbji, &b_xbji);
   fChain->SetBranchAddress("thqi", &thqi, &b_thqi);
   fChain->SetBranchAddress("sighad", &sighad, &b_sighad);
   fChain->SetBranchAddress("jacobian", &jacobian, &b_jacobian);
   fChain->SetBranchAddress("centjac", &centjac, &b_centjac);
   fChain->SetBranchAddress("pfermi", &pfermi, &b_pfermi);
   fChain->SetBranchAddress("xfermi", &xfermi, &b_xfermi);
   fChain->SetBranchAddress("phipqi", &phipqi, &b_phipqi);
   fChain->SetBranchAddress("xcal_gamma1", &xcal_gamma1, &b_xcal_gamma1);
   fChain->SetBranchAddress("ycal_gamma1", &ycal_gamma1, &b_ycal_gamma1);
   fChain->SetBranchAddress("Egamma1", &Egamma1, &b_Egamma1);
   fChain->SetBranchAddress("Pgamma1x", &Pgamma1x, &b_Pgamma1x);
   fChain->SetBranchAddress("Pgamma1y", &Pgamma1y, &b_Pgamma1y);
   fChain->SetBranchAddress("Pgamma1z", &Pgamma1z, &b_Pgamma1z);
   fChain->SetBranchAddress("xcal_gamma2", &xcal_gamma2, &b_xcal_gamma2);
   fChain->SetBranchAddress("ycal_gamma2", &ycal_gamma2, &b_ycal_gamma2);
   fChain->SetBranchAddress("Egamma2", &Egamma2, &b_Egamma2);
   fChain->SetBranchAddress("Pgamma2x", &Pgamma2x, &b_Pgamma2x);
   fChain->SetBranchAddress("Pgamma2y", &Pgamma2y, &b_Pgamma2y);
   fChain->SetBranchAddress("Pgamma2z", &Pgamma2z, &b_Pgamma2z);
   fChain->SetBranchAddress("p_e", &p_e, &b_p_e);
   fChain->SetBranchAddress("ux_e", &ux_e, &b_ux_e);
   fChain->SetBranchAddress("uy_e", &uy_e, &b_uy_e);
   fChain->SetBranchAddress("uz_e", &uz_e, &b_uz_e);
   Notify();
}

Bool_t simc_pi0_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void simc_pi0_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t simc_pi0_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef simc_pi0_tree_cxx
