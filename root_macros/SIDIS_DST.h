//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 30 15:55:18 2014 by ROOT version 5.34/18
// from TTree Tout/SIDIS reduced data file
// found on file: DST_pip_ebeam11gev_sbs10deg_0.root
//////////////////////////////////////////////////////////

#ifndef SIDIS_DST_h
#define SIDIS_DST_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class SIDIS_DST {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        Weight;
   Int_t           Hadron;
   Int_t           Nucleon;
   Double_t        ThetaBB;
   Double_t        ThetaSBS;
   Double_t        Ebeam;
   Double_t        vx;
   Double_t        vy;
   Double_t        vz;
   Double_t        ep;
   Double_t        eth;
   Double_t        eph;
   Double_t        Eh;
   Double_t        hp;
   Double_t        hth;
   Double_t        hph;
   Double_t        x;
   Double_t        y;
   Double_t        z;
   Double_t        Q2;
   Double_t        pT;
   Double_t        phih;
   Double_t        W2;
   Double_t        MX2;
   Double_t        exfp;
   Double_t        eyfp;
   Double_t        expfp;
   Double_t        eypfp;
   Double_t        hxfp;
   Double_t        hyfp;
   Double_t        hxpfp;
   Double_t        hypfp;
   Double_t        extar;
   Double_t        eytar;
   Double_t        exptar;
   Double_t        eyptar;
   Double_t        evzrecon;
   Double_t        eprecon;
   Double_t        ethrecon;
   Double_t        ephrecon;
   Double_t        hxtar;
   Double_t        hytar;
   Double_t        hxptar;
   Double_t        hyptar;
   Double_t        hvzrecon;
   Double_t        hprecon;
   Double_t        hthrecon;
   Double_t        hphrecon;
   Double_t        xrecon;
   Double_t        yrecon;
   Double_t        zrecon;
   Double_t        Q2recon;
   Double_t        pTrecon;
   Double_t        phihrecon;
   Double_t        W2recon;
   Double_t        MX2recon;
   Double_t        mpi0recon;

   // List of branches
   TBranch        *b_Weight;   //!
   TBranch        *b_Hadron;   //!
   TBranch        *b_Nucleon;   //!
   TBranch        *b_ThetaBB;   //!
   TBranch        *b_ThetaSBS;   //!
   TBranch        *b_Ebeam;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_ep;   //!
   TBranch        *b_eth;   //!
   TBranch        *b_eph;   //!
   TBranch        *b_Eh;   //!
   TBranch        *b_hp;   //!
   TBranch        *b_hth;   //!
   TBranch        *b_hph;   //!
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_pT;   //!
   TBranch        *b_phih;   //!
   TBranch        *b_W2;   //!
   TBranch        *b_MX2;   //!
   TBranch        *b_exfp;   //!
   TBranch        *b_eyfp;   //!
   TBranch        *b_expfp;   //!
   TBranch        *b_eypfp;   //!
   TBranch        *b_hxfp;   //!
   TBranch        *b_hyfp;   //!
   TBranch        *b_hxpfp;   //!
   TBranch        *b_hypfp;   //!
   TBranch        *b_extar;   //!
   TBranch        *b_eytar;   //!
   TBranch        *b_exptar;   //!
   TBranch        *b_eyptar;   //!
   TBranch        *b_evzrecon;   //!
   TBranch        *b_eprecon;   //!
   TBranch        *b_ethrecon;   //!
   TBranch        *b_ephrecon;   //!
   TBranch        *b_hxtar;   //!
   TBranch        *b_hytar;   //!
   TBranch        *b_hxptar;   //!
   TBranch        *b_hyptar;   //!
   TBranch        *b_hvzrecon;   //!
   TBranch        *b_hprecon;   //!
   TBranch        *b_hthrecon;   //!
   TBranch        *b_hphrecon;   //!
   TBranch        *b_xrecon;   //!
   TBranch        *b_yrecon;   //!
   TBranch        *b_zrecon;   //!
   TBranch        *b_Q2recon;   //!
   TBranch        *b_pTrecon;   //!
   TBranch        *b_phihrecon;   //!
   TBranch        *b_W2recon;   //!
   TBranch        *b_MX2recon;   //!
   TBranch        *b_mpi0recon;   //!

   SIDIS_DST(TTree *tree=0);
   virtual ~SIDIS_DST();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef SIDIS_DST_cxx
SIDIS_DST::SIDIS_DST(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("DST_pip_ebeam11gev_sbs10deg_0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("DST_pip_ebeam11gev_sbs10deg_0.root");
      }
      f->GetObject("Tout",tree);

   }
   Init(tree);
}

SIDIS_DST::~SIDIS_DST()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t SIDIS_DST::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t SIDIS_DST::LoadTree(Long64_t entry)
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

void SIDIS_DST::Init(TTree *tree)
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

   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("Hadron", &Hadron, &b_Hadron);
   fChain->SetBranchAddress("Nucleon", &Nucleon, &b_Nucleon);
   fChain->SetBranchAddress("ThetaBB", &ThetaBB, &b_ThetaBB);
   fChain->SetBranchAddress("ThetaSBS", &ThetaSBS, &b_ThetaSBS);
   fChain->SetBranchAddress("Ebeam", &Ebeam, &b_Ebeam);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("ep", &ep, &b_ep);
   fChain->SetBranchAddress("eth", &eth, &b_eth);
   fChain->SetBranchAddress("eph", &eph, &b_eph);
   fChain->SetBranchAddress("Eh", &Eh, &b_Eh);
   fChain->SetBranchAddress("hp", &hp, &b_hp);
   fChain->SetBranchAddress("hth", &hth, &b_hth);
   fChain->SetBranchAddress("hph", &hph, &b_hph);
   fChain->SetBranchAddress("x", &x, &b_x);
   fChain->SetBranchAddress("y", &y, &b_y);
   fChain->SetBranchAddress("z", &z, &b_z);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
   fChain->SetBranchAddress("pT", &pT, &b_pT);
   fChain->SetBranchAddress("phih", &phih, &b_phih);
   fChain->SetBranchAddress("W2", &W2, &b_W2);
   fChain->SetBranchAddress("MX2", &MX2, &b_MX2);
   fChain->SetBranchAddress("exfp", &exfp, &b_exfp);
   fChain->SetBranchAddress("eyfp", &eyfp, &b_eyfp);
   fChain->SetBranchAddress("expfp", &expfp, &b_expfp);
   fChain->SetBranchAddress("eypfp", &eypfp, &b_eypfp);
   fChain->SetBranchAddress("hxfp", &hxfp, &b_hxfp);
   fChain->SetBranchAddress("hyfp", &hyfp, &b_hyfp);
   fChain->SetBranchAddress("hxpfp", &hxpfp, &b_hxpfp);
   fChain->SetBranchAddress("hypfp", &hypfp, &b_hypfp);
   fChain->SetBranchAddress("extar", &extar, &b_extar);
   fChain->SetBranchAddress("eytar", &eytar, &b_eytar);
   fChain->SetBranchAddress("exptar", &exptar, &b_exptar);
   fChain->SetBranchAddress("eyptar", &eyptar, &b_eyptar);
   fChain->SetBranchAddress("evzrecon", &evzrecon, &b_evzrecon);
   fChain->SetBranchAddress("eprecon", &eprecon, &b_eprecon);
   fChain->SetBranchAddress("ethrecon", &ethrecon, &b_ethrecon);
   fChain->SetBranchAddress("ephrecon", &ephrecon, &b_ephrecon);
   fChain->SetBranchAddress("hxtar", &hxtar, &b_hxtar);
   fChain->SetBranchAddress("hytar", &hytar, &b_hytar);
   fChain->SetBranchAddress("hxptar", &hxptar, &b_hxptar);
   fChain->SetBranchAddress("hyptar", &hyptar, &b_hyptar);
   fChain->SetBranchAddress("hvzrecon", &hvzrecon, &b_hvzrecon);
   fChain->SetBranchAddress("hprecon", &hprecon, &b_hprecon);
   fChain->SetBranchAddress("hthrecon", &hthrecon, &b_hthrecon);
   fChain->SetBranchAddress("hphrecon", &hphrecon, &b_hphrecon);
   fChain->SetBranchAddress("xrecon", &xrecon, &b_xrecon);
   fChain->SetBranchAddress("yrecon", &yrecon, &b_yrecon);
   fChain->SetBranchAddress("zrecon", &zrecon, &b_zrecon);
   fChain->SetBranchAddress("Q2recon", &Q2recon, &b_Q2recon);
   fChain->SetBranchAddress("pTrecon", &pTrecon, &b_pTrecon);
   fChain->SetBranchAddress("phihrecon", &phihrecon, &b_phihrecon);
   fChain->SetBranchAddress("W2recon", &W2recon, &b_W2recon);
   fChain->SetBranchAddress("MX2recon", &MX2recon, &b_MX2recon);
   fChain->SetBranchAddress("mpi0recon", &mpi0recon, &b_mpi0recon);
   Notify();
}

Bool_t SIDIS_DST::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void SIDIS_DST::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t SIDIS_DST::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef SIDIS_DST_cxx
