//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Mar 12 10:18:27 2014 by ROOT version 5.34/12
// from TTree T/Geant4 SBS Simulation
// found on file: sidis_pip_11gev_upbend2.root
//////////////////////////////////////////////////////////

#ifndef g4sbs_tree_h
#define g4sbs_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class g4sbs_tree {
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
   Int_t           ev_nucl;
   Int_t           ev_fnucl;
   Int_t           ev_hadr;
   Int_t           ev_earmaccept;
   Int_t           ev_harmaccept;
   Double_t        gen_thbb;
   Double_t        gen_thhcal;
   Double_t        gen_dbb;
   Double_t        gen_dhcal;
   Double_t        gen_Ebeam;
   Int_t           ht_ndata;
   Int_t           ht_gid[146];   //[ht.ndata]
   Int_t           ht_trkrid[146];   //[ht.ndata]
   Double_t        ht_x[146];   //[ht.ndata]
   Double_t        ht_y[146];   //[ht.ndata]
   Double_t        ht_z[146];   //[ht.ndata]
   Double_t        ht_t[146];   //[ht.ndata]
   Double_t        ht_vx[146];   //[ht.ndata]
   Double_t        ht_vy[146];   //[ht.ndata]
   Double_t        ht_vz[146];   //[ht.ndata]
   Double_t        ht_dx[146];   //[ht.ndata]
   Double_t        ht_dy[146];   //[ht.ndata]
   Double_t        ht_p[146];   //[ht.ndata]
   Int_t           ht_trid[146];   //[ht.ndata]
   Int_t           ht_pid[146];   //[ht.ndata]
   Int_t           ht_mid[146];   //[ht.ndata]
   Double_t        ht_edep[146];   //[ht.ndata]
   Double_t        ht_tx[146];   //[ht.ndata]
   Double_t        ht_ty[146];   //[ht.ndata]
   Double_t        ht_txp[146];   //[ht.ndata]
   Double_t        ht_typ[146];   //[ht.ndata]
   Int_t           hc_ndata;
   Double_t        hc_x[298];   //[hc.ndata]
   Double_t        hc_y[298];   //[hc.ndata]
   Double_t        hc_z[298];   //[hc.ndata]
   Double_t        hc_e[298];   //[hc.ndata]
   Double_t        hc_t[298];   //[hc.ndata]
   Double_t        hc_vx[298];   //[hc.ndata]
   Double_t        hc_vy[298];   //[hc.ndata]
   Double_t        hc_vz[298];   //[hc.ndata]
   Int_t           hc_trid[298];   //[hc.ndata]
   Int_t           hc_mid[298];   //[hc.ndata]
   Int_t           hc_pid[298];   //[hc.ndata]
   Int_t           bc_ndata;
   Double_t        bc_x[12];   //[bc.ndata]
   Double_t        bc_y[12];   //[bc.ndata]
   Double_t        bc_z[12];   //[bc.ndata]
   Double_t        bc_e[12];   //[bc.ndata]
   Double_t        bc_t[12];   //[bc.ndata]
   Double_t        bc_vx[12];   //[bc.ndata]
   Double_t        bc_vy[12];   //[bc.ndata]
   Double_t        bc_vz[12];   //[bc.ndata]
   Int_t           bc_trid[12];   //[bc.ndata]
   Int_t           bc_mid[12];   //[bc.ndata]
   Int_t           bc_pid[12];   //[bc.ndata]
   Int_t           ntracks;
   vector<int>     *trackerid;
   vector<int>     *trackid;
   vector<int>     *trackpid;
   vector<int>     *tracknhits;
   vector<int>     *tracknplanes;
   vector<int>     *trackndf;
   vector<double>  *trackchi2;
   vector<double>  *trackchi2true;
   vector<double>  *trackx;
   vector<double>  *tracky;
   vector<double>  *trackxp;
   vector<double>  *trackyp;
   vector<double>  *trackt;
   vector<double>  *trackp;
   vector<double>  *trackxfit;
   vector<double>  *trackyfit;
   vector<double>  *trackxpfit;
   vector<double>  *trackypfit;
   Int_t           RICH_nhits;
   vector<int>     *RICH_pmt;
   vector<int>     *RICH_row;
   vector<int>     *RICH_col;
   vector<int>     *RICH_nphe;
   vector<double>  *RICH_tavg;
   vector<double>  *RICH_trms;
   vector<int>     *RICH_mID;
   vector<int>     *RICH_vol;
   vector<double>  *RICH_xhit;
   vector<double>  *RICH_yhit;
   vector<double>  *RICH_zhit;
   vector<double>  *RICH_pxhit;
   vector<double>  *RICH_pyhit;
   vector<double>  *RICH_pzhit;
   vector<double>  *RICH_vxhit;
   vector<double>  *RICH_vyhit;
   vector<double>  *RICH_vzhit;
   vector<double>  *RICH_vpxhit;
   vector<double>  *RICH_vpyhit;
   vector<double>  *RICH_vpzhit;
   Int_t           RICH_ntracks;
   vector<int>     *RICH_mPID;
   vector<int>     *RICH_mTID;
   vector<int>     *RICH_mMID;
   vector<double>  *RICH_mvx;
   vector<double>  *RICH_mvy;
   vector<double>  *RICH_mvz;
   vector<double>  *RICH_mpx;
   vector<double>  *RICH_mpy;
   vector<double>  *RICH_mpz;

   // List of branches
   TBranch        *b_ev;   //!
   TBranch        *b_gen;   //!
   TBranch        *b_ht_ndata;   //!
   TBranch        *b_ht_gid;   //!
   TBranch        *b_ht_trkrid;   //!
   TBranch        *b_ht_x;   //!
   TBranch        *b_ht_y;   //!
   TBranch        *b_ht_z;   //!
   TBranch        *b_ht_t;   //!
   TBranch        *b_ht_vx;   //!
   TBranch        *b_ht_vy;   //!
   TBranch        *b_ht_vz;   //!
   TBranch        *b_ht_dx;   //!
   TBranch        *b_ht_dy;   //!
   TBranch        *b_ht_p;   //!
   TBranch        *b_ht_trid;   //!
   TBranch        *b_ht_pid;   //!
   TBranch        *b_ht_mid;   //!
   TBranch        *b_ht_edep;   //!
   TBranch        *b_ht_tx;   //!
   TBranch        *b_ht_ty;   //!
   TBranch        *b_ht_txp;   //!
   TBranch        *b_ht_typ;   //!
   TBranch        *b_hc_ndata;   //!
   TBranch        *b_hc_x;   //!
   TBranch        *b_hc_y;   //!
   TBranch        *b_hc_z;   //!
   TBranch        *b_hc_e;   //!
   TBranch        *b_hc_t;   //!
   TBranch        *b_hc_vx;   //!
   TBranch        *b_hc_vy;   //!
   TBranch        *b_hc_vz;   //!
   TBranch        *b_hc_trid;   //!
   TBranch        *b_hc_mid;   //!
   TBranch        *b_hc_pid;   //!
   TBranch        *b_bc_ndata;   //!
   TBranch        *b_bc_x;   //!
   TBranch        *b_bc_y;   //!
   TBranch        *b_bc_z;   //!
   TBranch        *b_bc_e;   //!
   TBranch        *b_bc_t;   //!
   TBranch        *b_bc_vx;   //!
   TBranch        *b_bc_vy;   //!
   TBranch        *b_bc_vz;   //!
   TBranch        *b_bc_trid;   //!
   TBranch        *b_bc_mid;   //!
   TBranch        *b_bc_pid;   //!
   TBranch        *b_ntracks;   //!
   TBranch        *b_trackerid;   //!
   TBranch        *b_trackid;   //!
   TBranch        *b_trackpid;   //!
   TBranch        *b_tracknhits;   //!
   TBranch        *b_tracknplanes;   //!
   TBranch        *b_trackndf;   //!
   TBranch        *b_trackchi2;   //!
   TBranch        *b_trackchi2true;   //!
   TBranch        *b_trackx;   //!
   TBranch        *b_tracky;   //!
   TBranch        *b_trackxp;   //!
   TBranch        *b_trackyp;   //!
   TBranch        *b_trackt;   //!
   TBranch        *b_trackp;   //!
   TBranch        *b_trackxfit;   //!
   TBranch        *b_trackyfit;   //!
   TBranch        *b_trackxpfit;   //!
   TBranch        *b_trackypfit;   //!
   TBranch        *b_nhits_RICH;   //!
   TBranch        *b_RICH_pmt;   //!
   TBranch        *b_RICH_row;   //!
   TBranch        *b_RICH_col;   //!
   TBranch        *b_RICH_nphe;   //!
   TBranch        *b_RICH_tavg;   //!
   TBranch        *b_RICH_trms;   //!
   TBranch        *b_RICH_mID;   //!
   TBranch        *b_RICH_vol;   //!
   TBranch        *b_RICH_xhit;   //!
   TBranch        *b_RICH_yhit;   //!
   TBranch        *b_RICH_zhit;   //!
   TBranch        *b_RICH_pxhit;   //!
   TBranch        *b_RICH_pyhit;   //!
   TBranch        *b_RICH_pzhit;   //!
   TBranch        *b_RICH_vxhit;   //!
   TBranch        *b_RICH_vyhit;   //!
   TBranch        *b_RICH_vzhit;   //!
   TBranch        *b_RICH_vpxhit;   //!
   TBranch        *b_RICH_vpyhit;   //!
   TBranch        *b_RICH_vpzhit;   //!
   TBranch        *b_ntracks_RICH;   //!
   TBranch        *b_RICH_mPID;   //!
   TBranch        *b_RICH_mTID;   //!
   TBranch        *b_RICH_mMID;   //!
   TBranch        *b_RICH_mvx;   //!
   TBranch        *b_RICH_mvy;   //!
   TBranch        *b_RICH_mvz;   //!
   TBranch        *b_RICH_mpx;   //!
   TBranch        *b_RICH_mpy;   //!
   TBranch        *b_RICH_mpz;   //!

   g4sbs_tree(TTree *tree=0);
   virtual ~g4sbs_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef g4sbs_tree_cxx
g4sbs_tree::g4sbs_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("sidis_pip_11gev_upbend2.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("sidis_pip_11gev_upbend2.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

g4sbs_tree::~g4sbs_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t g4sbs_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t g4sbs_tree::LoadTree(Long64_t entry)
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

void g4sbs_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   trackerid = 0;
   trackid = 0;
   trackpid = 0;
   tracknhits = 0;
   tracknplanes = 0;
   trackndf = 0;
   trackchi2 = 0;
   trackchi2true = 0;
   trackx = 0;
   tracky = 0;
   trackxp = 0;
   trackyp = 0;
   trackt = 0;
   trackp = 0;
   trackxfit = 0;
   trackyfit = 0;
   trackxpfit = 0;
   trackypfit = 0;
   RICH_pmt = 0;
   RICH_row = 0;
   RICH_col = 0;
   RICH_nphe = 0;
   RICH_tavg = 0;
   RICH_trms = 0;
   RICH_mID = 0;
   RICH_vol = 0;
   RICH_xhit = 0;
   RICH_yhit = 0;
   RICH_zhit = 0;
   RICH_pxhit = 0;
   RICH_pyhit = 0;
   RICH_pzhit = 0;
   RICH_vxhit = 0;
   RICH_vyhit = 0;
   RICH_vzhit = 0;
   RICH_vpxhit = 0;
   RICH_vpyhit = 0;
   RICH_vpzhit = 0;
   RICH_mPID = 0;
   RICH_mTID = 0;
   RICH_mMID = 0;
   RICH_mvx = 0;
   RICH_mvy = 0;
   RICH_mvz = 0;
   RICH_mpx = 0;
   RICH_mpy = 0;
   RICH_mpz = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev", &ev_count, &b_ev);
   fChain->SetBranchAddress("gen", &gen_thbb, &b_gen);
   fChain->SetBranchAddress("ht.ndata", &ht_ndata, &b_ht_ndata);
   fChain->SetBranchAddress("ht.gid", ht_gid, &b_ht_gid);
   fChain->SetBranchAddress("ht.trkrid", ht_trkrid, &b_ht_trkrid);
   fChain->SetBranchAddress("ht.x", ht_x, &b_ht_x);
   fChain->SetBranchAddress("ht.y", ht_y, &b_ht_y);
   fChain->SetBranchAddress("ht.z", ht_z, &b_ht_z);
   fChain->SetBranchAddress("ht.t", ht_t, &b_ht_t);
   fChain->SetBranchAddress("ht.vx", ht_vx, &b_ht_vx);
   fChain->SetBranchAddress("ht.vy", ht_vy, &b_ht_vy);
   fChain->SetBranchAddress("ht.vz", ht_vz, &b_ht_vz);
   fChain->SetBranchAddress("ht.dx", ht_dx, &b_ht_dx);
   fChain->SetBranchAddress("ht.dy", ht_dy, &b_ht_dy);
   fChain->SetBranchAddress("ht.p", ht_p, &b_ht_p);
   fChain->SetBranchAddress("ht.trid", ht_trid, &b_ht_trid);
   fChain->SetBranchAddress("ht.pid", ht_pid, &b_ht_pid);
   fChain->SetBranchAddress("ht.mid", ht_mid, &b_ht_mid);
   fChain->SetBranchAddress("ht.edep", ht_edep, &b_ht_edep);
   fChain->SetBranchAddress("ht.tx", ht_tx, &b_ht_tx);
   fChain->SetBranchAddress("ht.ty", ht_ty, &b_ht_ty);
   fChain->SetBranchAddress("ht.txp", ht_txp, &b_ht_txp);
   fChain->SetBranchAddress("ht.typ", ht_typ, &b_ht_typ);
   fChain->SetBranchAddress("hc.ndata", &hc_ndata, &b_hc_ndata);
   fChain->SetBranchAddress("hc.x", hc_x, &b_hc_x);
   fChain->SetBranchAddress("hc.y", hc_y, &b_hc_y);
   fChain->SetBranchAddress("hc.z", hc_z, &b_hc_z);
   fChain->SetBranchAddress("hc.e", hc_e, &b_hc_e);
   fChain->SetBranchAddress("hc.t", hc_t, &b_hc_t);
   fChain->SetBranchAddress("hc.vx", hc_vx, &b_hc_vx);
   fChain->SetBranchAddress("hc.vy", hc_vy, &b_hc_vy);
   fChain->SetBranchAddress("hc.vz", hc_vz, &b_hc_vz);
   fChain->SetBranchAddress("hc.trid", hc_trid, &b_hc_trid);
   fChain->SetBranchAddress("hc.mid", hc_mid, &b_hc_mid);
   fChain->SetBranchAddress("hc.pid", hc_pid, &b_hc_pid);
   fChain->SetBranchAddress("bc.ndata", &bc_ndata, &b_bc_ndata);
   fChain->SetBranchAddress("bc.x", bc_x, &b_bc_x);
   fChain->SetBranchAddress("bc.y", bc_y, &b_bc_y);
   fChain->SetBranchAddress("bc.z", bc_z, &b_bc_z);
   fChain->SetBranchAddress("bc.e", bc_e, &b_bc_e);
   fChain->SetBranchAddress("bc.t", bc_t, &b_bc_t);
   fChain->SetBranchAddress("bc.vx", bc_vx, &b_bc_vx);
   fChain->SetBranchAddress("bc.vy", bc_vy, &b_bc_vy);
   fChain->SetBranchAddress("bc.vz", bc_vz, &b_bc_vz);
   fChain->SetBranchAddress("bc.trid", bc_trid, &b_bc_trid);
   fChain->SetBranchAddress("bc.mid", bc_mid, &b_bc_mid);
   fChain->SetBranchAddress("bc.pid", bc_pid, &b_bc_pid);
   fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
   fChain->SetBranchAddress("trackerid", &trackerid, &b_trackerid);
   fChain->SetBranchAddress("trackid", &trackid, &b_trackid);
   fChain->SetBranchAddress("trackpid", &trackpid, &b_trackpid);
   fChain->SetBranchAddress("tracknhits", &tracknhits, &b_tracknhits);
   fChain->SetBranchAddress("tracknplanes", &tracknplanes, &b_tracknplanes);
   fChain->SetBranchAddress("trackndf", &trackndf, &b_trackndf);
   fChain->SetBranchAddress("trackchi2", &trackchi2, &b_trackchi2);
   fChain->SetBranchAddress("trackchi2true", &trackchi2true, &b_trackchi2true);
   fChain->SetBranchAddress("trackx", &trackx, &b_trackx);
   fChain->SetBranchAddress("tracky", &tracky, &b_tracky);
   fChain->SetBranchAddress("trackxp", &trackxp, &b_trackxp);
   fChain->SetBranchAddress("trackyp", &trackyp, &b_trackyp);
   fChain->SetBranchAddress("trackt", &trackt, &b_trackt);
   fChain->SetBranchAddress("trackp", &trackp, &b_trackp);
   fChain->SetBranchAddress("trackxfit", &trackxfit, &b_trackxfit);
   fChain->SetBranchAddress("trackyfit", &trackyfit, &b_trackyfit);
   fChain->SetBranchAddress("trackxpfit", &trackxpfit, &b_trackxpfit);
   fChain->SetBranchAddress("trackypfit", &trackypfit, &b_trackypfit);
   fChain->SetBranchAddress("RICH_nhits", &RICH_nhits, &b_nhits_RICH);
   fChain->SetBranchAddress("RICH_pmt", &RICH_pmt, &b_RICH_pmt);
   fChain->SetBranchAddress("RICH_row", &RICH_row, &b_RICH_row);
   fChain->SetBranchAddress("RICH_col", &RICH_col, &b_RICH_col);
   fChain->SetBranchAddress("RICH_nphe", &RICH_nphe, &b_RICH_nphe);
   fChain->SetBranchAddress("RICH_tavg", &RICH_tavg, &b_RICH_tavg);
   fChain->SetBranchAddress("RICH_trms", &RICH_trms, &b_RICH_trms);
   fChain->SetBranchAddress("RICH_mID", &RICH_mID, &b_RICH_mID);
   fChain->SetBranchAddress("RICH_vol", &RICH_vol, &b_RICH_vol);
   fChain->SetBranchAddress("RICH_xhit", &RICH_xhit, &b_RICH_xhit);
   fChain->SetBranchAddress("RICH_yhit", &RICH_yhit, &b_RICH_yhit);
   fChain->SetBranchAddress("RICH_zhit", &RICH_zhit, &b_RICH_zhit);
   fChain->SetBranchAddress("RICH_pxhit", &RICH_pxhit, &b_RICH_pxhit);
   fChain->SetBranchAddress("RICH_pyhit", &RICH_pyhit, &b_RICH_pyhit);
   fChain->SetBranchAddress("RICH_pzhit", &RICH_pzhit, &b_RICH_pzhit);
   fChain->SetBranchAddress("RICH_vxhit", &RICH_vxhit, &b_RICH_vxhit);
   fChain->SetBranchAddress("RICH_vyhit", &RICH_vyhit, &b_RICH_vyhit);
   fChain->SetBranchAddress("RICH_vzhit", &RICH_vzhit, &b_RICH_vzhit);
   fChain->SetBranchAddress("RICH_vpxhit", &RICH_vpxhit, &b_RICH_vpxhit);
   fChain->SetBranchAddress("RICH_vpyhit", &RICH_vpyhit, &b_RICH_vpyhit);
   fChain->SetBranchAddress("RICH_vpzhit", &RICH_vpzhit, &b_RICH_vpzhit);
   fChain->SetBranchAddress("RICH_ntracks", &RICH_ntracks, &b_ntracks_RICH);
   fChain->SetBranchAddress("RICH_mPID", &RICH_mPID, &b_RICH_mPID);
   fChain->SetBranchAddress("RICH_mTID", &RICH_mTID, &b_RICH_mTID);
   fChain->SetBranchAddress("RICH_mMID", &RICH_mMID, &b_RICH_mMID);
   fChain->SetBranchAddress("RICH_mvx", &RICH_mvx, &b_RICH_mvx);
   fChain->SetBranchAddress("RICH_mvy", &RICH_mvy, &b_RICH_mvy);
   fChain->SetBranchAddress("RICH_mvz", &RICH_mvz, &b_RICH_mvz);
   fChain->SetBranchAddress("RICH_mpx", &RICH_mpx, &b_RICH_mpx);
   fChain->SetBranchAddress("RICH_mpy", &RICH_mpy, &b_RICH_mpy);
   fChain->SetBranchAddress("RICH_mpz", &RICH_mpz, &b_RICH_mpz);
   Notify();
}

Bool_t g4sbs_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void g4sbs_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t g4sbs_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef g4sbs_tree_cxx
