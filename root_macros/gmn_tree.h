//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Jul 10 22:51:55 2017 by ROOT version 6.08/02
// from TTree T/Geant4 SBS Simulation
// found on file: beam_bkgd_1695.root
//////////////////////////////////////////////////////////

#ifndef gmn_tree_h
#define gmn_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class gmn_tree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

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
   Double_t        gen_voffhcal;
   Double_t        gen_drich;
   Double_t        gen_dsbstrkr;
   Double_t        gen_Ebeam;
   Int_t           Earm_BBGEM_hit_nhits;
   vector<int>     *Earm_BBGEM_hit_plane;
   vector<int>     *Earm_BBGEM_hit_strip;
   vector<double>  *Earm_BBGEM_hit_x;
   vector<double>  *Earm_BBGEM_hit_y;
   vector<double>  *Earm_BBGEM_hit_z;
   vector<double>  *Earm_BBGEM_hit_polx;
   vector<double>  *Earm_BBGEM_hit_poly;
   vector<double>  *Earm_BBGEM_hit_polz;
   vector<double>  *Earm_BBGEM_hit_t;
   vector<double>  *Earm_BBGEM_hit_trms;
   vector<double>  *Earm_BBGEM_hit_tmin;
   vector<double>  *Earm_BBGEM_hit_tmax;
   vector<double>  *Earm_BBGEM_hit_tx;
   vector<double>  *Earm_BBGEM_hit_ty;
   vector<double>  *Earm_BBGEM_hit_txp;
   vector<double>  *Earm_BBGEM_hit_typ;
   vector<double>  *Earm_BBGEM_hit_xg;
   vector<double>  *Earm_BBGEM_hit_yg;
   vector<double>  *Earm_BBGEM_hit_zg;
   vector<int>     *Earm_BBGEM_hit_trid;
   vector<int>     *Earm_BBGEM_hit_mid;
   vector<int>     *Earm_BBGEM_hit_pid;
   vector<double>  *Earm_BBGEM_hit_vx;
   vector<double>  *Earm_BBGEM_hit_vy;
   vector<double>  *Earm_BBGEM_hit_vz;
   vector<double>  *Earm_BBGEM_hit_p;
   vector<double>  *Earm_BBGEM_hit_edep;
   vector<double>  *Earm_BBGEM_hit_beta;
   Int_t           Earm_BBGEM_Track_ntracks;
   vector<int>     *Earm_BBGEM_Track_TID;
   vector<int>     *Earm_BBGEM_Track_PID;
   vector<int>     *Earm_BBGEM_Track_MID;
   vector<int>     *Earm_BBGEM_Track_NumHits;
   vector<int>     *Earm_BBGEM_Track_NumPlanes;
   vector<int>     *Earm_BBGEM_Track_NDF;
   vector<double>  *Earm_BBGEM_Track_Chi2fit;
   vector<double>  *Earm_BBGEM_Track_Chi2true;
   vector<double>  *Earm_BBGEM_Track_X;
   vector<double>  *Earm_BBGEM_Track_Y;
   vector<double>  *Earm_BBGEM_Track_Xp;
   vector<double>  *Earm_BBGEM_Track_Yp;
   vector<double>  *Earm_BBGEM_Track_T;
   vector<double>  *Earm_BBGEM_Track_P;
   vector<double>  *Earm_BBGEM_Track_Sx;
   vector<double>  *Earm_BBGEM_Track_Sy;
   vector<double>  *Earm_BBGEM_Track_Sz;
   vector<double>  *Earm_BBGEM_Track_Xfit;
   vector<double>  *Earm_BBGEM_Track_Yfit;
   vector<double>  *Earm_BBGEM_Track_Xpfit;
   vector<double>  *Earm_BBGEM_Track_Ypfit;
   Int_t           Earm_BBHodoScint_hit_nhits;
   vector<int>     *Earm_BBHodoScint_hit_row;
   vector<int>     *Earm_BBHodoScint_hit_col;
   vector<int>     *Earm_BBHodoScint_hit_cell;
   vector<int>     *Earm_BBHodoScint_hit_plane;
   vector<double>  *Earm_BBHodoScint_hit_xcell;
   vector<double>  *Earm_BBHodoScint_hit_ycell;
   vector<double>  *Earm_BBHodoScint_hit_zcell;
   vector<double>  *Earm_BBHodoScint_hit_xcellg;
   vector<double>  *Earm_BBHodoScint_hit_ycellg;
   vector<double>  *Earm_BBHodoScint_hit_zcellg;
   vector<double>  *Earm_BBHodoScint_hit_xhit;
   vector<double>  *Earm_BBHodoScint_hit_yhit;
   vector<double>  *Earm_BBHodoScint_hit_zhit;
   vector<double>  *Earm_BBHodoScint_hit_sumedep;
   vector<double>  *Earm_BBHodoScint_hit_tavg;
   vector<double>  *Earm_BBHodoScint_hit_trms;
   vector<double>  *Earm_BBHodoScint_hit_tmin;
   vector<double>  *Earm_BBHodoScint_hit_tmax;
   Int_t           Earm_BBPS_hit_nhits;
   vector<int>     *Earm_BBPS_hit_PMT;
   vector<int>     *Earm_BBPS_hit_row;
   vector<int>     *Earm_BBPS_hit_col;
   vector<int>     *Earm_BBPS_hit_plane;
   vector<double>  *Earm_BBPS_hit_xcell;
   vector<double>  *Earm_BBPS_hit_ycell;
   vector<double>  *Earm_BBPS_hit_zcell;
   vector<double>  *Earm_BBPS_hit_xgcell;
   vector<double>  *Earm_BBPS_hit_ygcell;
   vector<double>  *Earm_BBPS_hit_zgcell;
   vector<int>     *Earm_BBPS_hit_NumPhotoelectrons;
   vector<double>  *Earm_BBPS_hit_Time_avg;
   vector<double>  *Earm_BBPS_hit_Time_rms;
   vector<double>  *Earm_BBPS_hit_Time_min;
   vector<double>  *Earm_BBPS_hit_Time_max;
   Int_t           Earm_BBPSTF1_hit_nhits;
   vector<int>     *Earm_BBPSTF1_hit_row;
   vector<int>     *Earm_BBPSTF1_hit_col;
   vector<int>     *Earm_BBPSTF1_hit_cell;
   vector<int>     *Earm_BBPSTF1_hit_plane;
   vector<double>  *Earm_BBPSTF1_hit_xcell;
   vector<double>  *Earm_BBPSTF1_hit_ycell;
   vector<double>  *Earm_BBPSTF1_hit_zcell;
   vector<double>  *Earm_BBPSTF1_hit_xcellg;
   vector<double>  *Earm_BBPSTF1_hit_ycellg;
   vector<double>  *Earm_BBPSTF1_hit_zcellg;
   vector<double>  *Earm_BBPSTF1_hit_xhit;
   vector<double>  *Earm_BBPSTF1_hit_yhit;
   vector<double>  *Earm_BBPSTF1_hit_zhit;
   vector<double>  *Earm_BBPSTF1_hit_sumedep;
   vector<double>  *Earm_BBPSTF1_hit_tavg;
   vector<double>  *Earm_BBPSTF1_hit_trms;
   vector<double>  *Earm_BBPSTF1_hit_tmin;
   vector<double>  *Earm_BBPSTF1_hit_tmax;
   Int_t           Earm_BBSH_hit_nhits;
   vector<int>     *Earm_BBSH_hit_PMT;
   vector<int>     *Earm_BBSH_hit_row;
   vector<int>     *Earm_BBSH_hit_col;
   vector<int>     *Earm_BBSH_hit_plane;
   vector<double>  *Earm_BBSH_hit_xcell;
   vector<double>  *Earm_BBSH_hit_ycell;
   vector<double>  *Earm_BBSH_hit_zcell;
   vector<double>  *Earm_BBSH_hit_xgcell;
   vector<double>  *Earm_BBSH_hit_ygcell;
   vector<double>  *Earm_BBSH_hit_zgcell;
   vector<int>     *Earm_BBSH_hit_NumPhotoelectrons;
   vector<double>  *Earm_BBSH_hit_Time_avg;
   vector<double>  *Earm_BBSH_hit_Time_rms;
   vector<double>  *Earm_BBSH_hit_Time_min;
   vector<double>  *Earm_BBSH_hit_Time_max;
   Int_t           Earm_BBSHTF1_hit_nhits;
   vector<int>     *Earm_BBSHTF1_hit_row;
   vector<int>     *Earm_BBSHTF1_hit_col;
   vector<int>     *Earm_BBSHTF1_hit_cell;
   vector<int>     *Earm_BBSHTF1_hit_plane;
   vector<double>  *Earm_BBSHTF1_hit_xcell;
   vector<double>  *Earm_BBSHTF1_hit_ycell;
   vector<double>  *Earm_BBSHTF1_hit_zcell;
   vector<double>  *Earm_BBSHTF1_hit_xcellg;
   vector<double>  *Earm_BBSHTF1_hit_ycellg;
   vector<double>  *Earm_BBSHTF1_hit_zcellg;
   vector<double>  *Earm_BBSHTF1_hit_xhit;
   vector<double>  *Earm_BBSHTF1_hit_yhit;
   vector<double>  *Earm_BBSHTF1_hit_zhit;
   vector<double>  *Earm_BBSHTF1_hit_sumedep;
   vector<double>  *Earm_BBSHTF1_hit_tavg;
   vector<double>  *Earm_BBSHTF1_hit_trms;
   vector<double>  *Earm_BBSHTF1_hit_tmin;
   vector<double>  *Earm_BBSHTF1_hit_tmax;
   Int_t           Earm_GC_PMT_Glass_hit_nhits;
   vector<int>     *Earm_GC_PMT_Glass_hit_row;
   vector<int>     *Earm_GC_PMT_Glass_hit_col;
   vector<int>     *Earm_GC_PMT_Glass_hit_cell;
   vector<int>     *Earm_GC_PMT_Glass_hit_plane;
   vector<double>  *Earm_GC_PMT_Glass_hit_xcell;
   vector<double>  *Earm_GC_PMT_Glass_hit_ycell;
   vector<double>  *Earm_GC_PMT_Glass_hit_zcell;
   vector<double>  *Earm_GC_PMT_Glass_hit_xcellg;
   vector<double>  *Earm_GC_PMT_Glass_hit_ycellg;
   vector<double>  *Earm_GC_PMT_Glass_hit_zcellg;
   vector<double>  *Earm_GC_PMT_Glass_hit_xhit;
   vector<double>  *Earm_GC_PMT_Glass_hit_yhit;
   vector<double>  *Earm_GC_PMT_Glass_hit_zhit;
   vector<double>  *Earm_GC_PMT_Glass_hit_sumedep;
   vector<double>  *Earm_GC_PMT_Glass_hit_tavg;
   vector<double>  *Earm_GC_PMT_Glass_hit_trms;
   vector<double>  *Earm_GC_PMT_Glass_hit_tmin;
   vector<double>  *Earm_GC_PMT_Glass_hit_tmax;
   Int_t           Earm_GEMElectronicsGMn_hit_nhits;
   vector<int>     *Earm_GEMElectronicsGMn_hit_row;
   vector<int>     *Earm_GEMElectronicsGMn_hit_col;
   vector<int>     *Earm_GEMElectronicsGMn_hit_cell;
   vector<int>     *Earm_GEMElectronicsGMn_hit_plane;
   vector<double>  *Earm_GEMElectronicsGMn_hit_xcell;
   vector<double>  *Earm_GEMElectronicsGMn_hit_ycell;
   vector<double>  *Earm_GEMElectronicsGMn_hit_zcell;
   vector<double>  *Earm_GEMElectronicsGMn_hit_xcellg;
   vector<double>  *Earm_GEMElectronicsGMn_hit_ycellg;
   vector<double>  *Earm_GEMElectronicsGMn_hit_zcellg;
   vector<double>  *Earm_GEMElectronicsGMn_hit_xhit;
   vector<double>  *Earm_GEMElectronicsGMn_hit_yhit;
   vector<double>  *Earm_GEMElectronicsGMn_hit_zhit;
   vector<double>  *Earm_GEMElectronicsGMn_hit_sumedep;
   vector<double>  *Earm_GEMElectronicsGMn_hit_tavg;
   vector<double>  *Earm_GEMElectronicsGMn_hit_trms;
   vector<double>  *Earm_GEMElectronicsGMn_hit_tmin;
   vector<double>  *Earm_GEMElectronicsGMn_hit_tmax;
   Int_t           Earm_GRINCH_hit_nhits;
   vector<int>     *Earm_GRINCH_hit_PMT;
   vector<int>     *Earm_GRINCH_hit_row;
   vector<int>     *Earm_GRINCH_hit_col;
   vector<double>  *Earm_GRINCH_hit_xpmt;
   vector<double>  *Earm_GRINCH_hit_ypmt;
   vector<double>  *Earm_GRINCH_hit_zpmt;
   vector<double>  *Earm_GRINCH_hit_xgpmt;
   vector<double>  *Earm_GRINCH_hit_ygpmt;
   vector<double>  *Earm_GRINCH_hit_zgpmt;
   vector<int>     *Earm_GRINCH_hit_NumPhotoelectrons;
   vector<double>  *Earm_GRINCH_hit_Time_avg;
   vector<double>  *Earm_GRINCH_hit_Time_rms;
   vector<double>  *Earm_GRINCH_hit_Time_min;
   vector<double>  *Earm_GRINCH_hit_Time_max;
   vector<int>     *Earm_GRINCH_hit_mTrackNo;
   vector<double>  *Earm_GRINCH_hit_xhit;
   vector<double>  *Earm_GRINCH_hit_yhit;
   vector<double>  *Earm_GRINCH_hit_zhit;
   vector<double>  *Earm_GRINCH_hit_pxhit;
   vector<double>  *Earm_GRINCH_hit_pyhit;
   vector<double>  *Earm_GRINCH_hit_pzhit;
   vector<double>  *Earm_GRINCH_hit_pvx;
   vector<double>  *Earm_GRINCH_hit_pvy;
   vector<double>  *Earm_GRINCH_hit_pvz;
   vector<double>  *Earm_GRINCH_hit_ppx;
   vector<double>  *Earm_GRINCH_hit_ppy;
   vector<double>  *Earm_GRINCH_hit_ppz;
   vector<int>     *Earm_GRINCH_hit_volume_flag;
   Int_t           Harm_CDET_hit_nhits;
   vector<int>     *Harm_CDET_hit_PMT;
   vector<int>     *Harm_CDET_hit_row;
   vector<int>     *Harm_CDET_hit_col;
   vector<int>     *Harm_CDET_hit_plane;
   vector<double>  *Harm_CDET_hit_xcell;
   vector<double>  *Harm_CDET_hit_ycell;
   vector<double>  *Harm_CDET_hit_zcell;
   vector<double>  *Harm_CDET_hit_xgcell;
   vector<double>  *Harm_CDET_hit_ygcell;
   vector<double>  *Harm_CDET_hit_zgcell;
   vector<int>     *Harm_CDET_hit_NumPhotoelectrons;
   vector<double>  *Harm_CDET_hit_Time_avg;
   vector<double>  *Harm_CDET_hit_Time_rms;
   vector<double>  *Harm_CDET_hit_Time_min;
   vector<double>  *Harm_CDET_hit_Time_max;
   Int_t           Harm_CDET_Scint_hit_nhits;
   vector<int>     *Harm_CDET_Scint_hit_row;
   vector<int>     *Harm_CDET_Scint_hit_col;
   vector<int>     *Harm_CDET_Scint_hit_cell;
   vector<int>     *Harm_CDET_Scint_hit_plane;
   vector<double>  *Harm_CDET_Scint_hit_xcell;
   vector<double>  *Harm_CDET_Scint_hit_ycell;
   vector<double>  *Harm_CDET_Scint_hit_zcell;
   vector<double>  *Harm_CDET_Scint_hit_xcellg;
   vector<double>  *Harm_CDET_Scint_hit_ycellg;
   vector<double>  *Harm_CDET_Scint_hit_zcellg;
   vector<double>  *Harm_CDET_Scint_hit_xhit;
   vector<double>  *Harm_CDET_Scint_hit_yhit;
   vector<double>  *Harm_CDET_Scint_hit_zhit;
   vector<double>  *Harm_CDET_Scint_hit_sumedep;
   vector<double>  *Harm_CDET_Scint_hit_tavg;
   vector<double>  *Harm_CDET_Scint_hit_trms;
   vector<double>  *Harm_CDET_Scint_hit_tmin;
   vector<double>  *Harm_CDET_Scint_hit_tmax;
   Int_t           Harm_HCal_hit_nhits;
   vector<int>     *Harm_HCal_hit_PMT;
   vector<int>     *Harm_HCal_hit_row;
   vector<int>     *Harm_HCal_hit_col;
   vector<int>     *Harm_HCal_hit_plane;
   vector<double>  *Harm_HCal_hit_xcell;
   vector<double>  *Harm_HCal_hit_ycell;
   vector<double>  *Harm_HCal_hit_zcell;
   vector<double>  *Harm_HCal_hit_xgcell;
   vector<double>  *Harm_HCal_hit_ygcell;
   vector<double>  *Harm_HCal_hit_zgcell;
   vector<int>     *Harm_HCal_hit_NumPhotoelectrons;
   vector<double>  *Harm_HCal_hit_Time_avg;
   vector<double>  *Harm_HCal_hit_Time_rms;
   vector<double>  *Harm_HCal_hit_Time_min;
   vector<double>  *Harm_HCal_hit_Time_max;
   Int_t           Harm_HCalScint_hit_nhits;
   vector<int>     *Harm_HCalScint_hit_row;
   vector<int>     *Harm_HCalScint_hit_col;
   vector<int>     *Harm_HCalScint_hit_cell;
   vector<int>     *Harm_HCalScint_hit_plane;
   vector<double>  *Harm_HCalScint_hit_xcell;
   vector<double>  *Harm_HCalScint_hit_ycell;
   vector<double>  *Harm_HCalScint_hit_zcell;
   vector<double>  *Harm_HCalScint_hit_xcellg;
   vector<double>  *Harm_HCalScint_hit_ycellg;
   vector<double>  *Harm_HCalScint_hit_zcellg;
   vector<double>  *Harm_HCalScint_hit_xhit;
   vector<double>  *Harm_HCalScint_hit_yhit;
   vector<double>  *Harm_HCalScint_hit_zhit;
   vector<double>  *Harm_HCalScint_hit_sumedep;
   vector<double>  *Harm_HCalScint_hit_tavg;
   vector<double>  *Harm_HCalScint_hit_trms;
   vector<double>  *Harm_HCalScint_hit_tmin;
   vector<double>  *Harm_HCalScint_hit_tmax;

   // List of branches
   TBranch        *b_ev;   //!
   TBranch        *b_gen;   //!
   TBranch        *b_Earm_BBGEM_hit_nhits;   //!
   TBranch        *b_Earm_BBGEM_hit_plane;   //!
   TBranch        *b_Earm_BBGEM_hit_strip;   //!
   TBranch        *b_Earm_BBGEM_hit_x;   //!
   TBranch        *b_Earm_BBGEM_hit_y;   //!
   TBranch        *b_Earm_BBGEM_hit_z;   //!
   TBranch        *b_Earm_BBGEM_hit_polx;   //!
   TBranch        *b_Earm_BBGEM_hit_poly;   //!
   TBranch        *b_Earm_BBGEM_hit_polz;   //!
   TBranch        *b_Earm_BBGEM_hit_t;   //!
   TBranch        *b_Earm_BBGEM_hit_trms;   //!
   TBranch        *b_Earm_BBGEM_hit_tmin;   //!
   TBranch        *b_Earm_BBGEM_hit_tmax;   //!
   TBranch        *b_Earm_BBGEM_hit_tx;   //!
   TBranch        *b_Earm_BBGEM_hit_ty;   //!
   TBranch        *b_Earm_BBGEM_hit_txp;   //!
   TBranch        *b_Earm_BBGEM_hit_typ;   //!
   TBranch        *b_Earm_BBGEM_hit_xg;   //!
   TBranch        *b_Earm_BBGEM_hit_yg;   //!
   TBranch        *b_Earm_BBGEM_hit_zg;   //!
   TBranch        *b_Earm_BBGEM_hit_trid;   //!
   TBranch        *b_Earm_BBGEM_hit_mid;   //!
   TBranch        *b_Earm_BBGEM_hit_pid;   //!
   TBranch        *b_Earm_BBGEM_hit_vx;   //!
   TBranch        *b_Earm_BBGEM_hit_vy;   //!
   TBranch        *b_Earm_BBGEM_hit_vz;   //!
   TBranch        *b_Earm_BBGEM_hit_p;   //!
   TBranch        *b_Earm_BBGEM_hit_edep;   //!
   TBranch        *b_Earm_BBGEM_hit_beta;   //!
   TBranch        *b_Earm_BBGEM_Track_ntracks;   //!
   TBranch        *b_Earm_BBGEM_Track_TID;   //!
   TBranch        *b_Earm_BBGEM_Track_PID;   //!
   TBranch        *b_Earm_BBGEM_Track_MID;   //!
   TBranch        *b_Earm_BBGEM_Track_NumHits;   //!
   TBranch        *b_Earm_BBGEM_Track_NumPlanes;   //!
   TBranch        *b_Earm_BBGEM_Track_NDF;   //!
   TBranch        *b_Earm_BBGEM_Track_Chi2fit;   //!
   TBranch        *b_Earm_BBGEM_Track_Chi2true;   //!
   TBranch        *b_Earm_BBGEM_Track_X;   //!
   TBranch        *b_Earm_BBGEM_Track_Y;   //!
   TBranch        *b_Earm_BBGEM_Track_Xp;   //!
   TBranch        *b_Earm_BBGEM_Track_Yp;   //!
   TBranch        *b_Earm_BBGEM_Track_T;   //!
   TBranch        *b_Earm_BBGEM_Track_P;   //!
   TBranch        *b_Earm_BBGEM_Track_Sx;   //!
   TBranch        *b_Earm_BBGEM_Track_Sy;   //!
   TBranch        *b_Earm_BBGEM_Track_Sz;   //!
   TBranch        *b_Earm_BBGEM_Track_Xfit;   //!
   TBranch        *b_Earm_BBGEM_Track_Yfit;   //!
   TBranch        *b_Earm_BBGEM_Track_Xpfit;   //!
   TBranch        *b_Earm_BBGEM_Track_Ypfit;   //!
   TBranch        *b_Earm_BBHodoScint_hit_nhits;   //!
   TBranch        *b_Earm_BBHodoScint_hit_row;   //!
   TBranch        *b_Earm_BBHodoScint_hit_col;   //!
   TBranch        *b_Earm_BBHodoScint_hit_cell;   //!
   TBranch        *b_Earm_BBHodoScint_hit_plane;   //!
   TBranch        *b_Earm_BBHodoScint_hit_xcell;   //!
   TBranch        *b_Earm_BBHodoScint_hit_ycell;   //!
   TBranch        *b_Earm_BBHodoScint_hit_zcell;   //!
   TBranch        *b_Earm_BBHodoScint_hit_xcellg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_ycellg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_zcellg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_xhit;   //!
   TBranch        *b_Earm_BBHodoScint_hit_yhit;   //!
   TBranch        *b_Earm_BBHodoScint_hit_zhit;   //!
   TBranch        *b_Earm_BBHodoScint_hit_sumedep;   //!
   TBranch        *b_Earm_BBHodoScint_hit_tavg;   //!
   TBranch        *b_Earm_BBHodoScint_hit_trms;   //!
   TBranch        *b_Earm_BBHodoScint_hit_tmin;   //!
   TBranch        *b_Earm_BBHodoScint_hit_tmax;   //!
   TBranch        *b_Earm_BBPS_hit_nhits;   //!
   TBranch        *b_Earm_BBPS_hit_PMT;   //!
   TBranch        *b_Earm_BBPS_hit_row;   //!
   TBranch        *b_Earm_BBPS_hit_col;   //!
   TBranch        *b_Earm_BBPS_hit_plane;   //!
   TBranch        *b_Earm_BBPS_hit_xcell;   //!
   TBranch        *b_Earm_BBPS_hit_ycell;   //!
   TBranch        *b_Earm_BBPS_hit_zcell;   //!
   TBranch        *b_Earm_BBPS_hit_xgcell;   //!
   TBranch        *b_Earm_BBPS_hit_ygcell;   //!
   TBranch        *b_Earm_BBPS_hit_zgcell;   //!
   TBranch        *b_Earm_BBPS_hit_NumPhotoelectrons;   //!
   TBranch        *b_Earm_BBPS_hit_Time_avg;   //!
   TBranch        *b_Earm_BBPS_hit_Time_rms;   //!
   TBranch        *b_Earm_BBPS_hit_Time_min;   //!
   TBranch        *b_Earm_BBPS_hit_Time_max;   //!
   TBranch        *b_Earm_BBPSTF1_hit_nhits;   //!
   TBranch        *b_Earm_BBPSTF1_hit_row;   //!
   TBranch        *b_Earm_BBPSTF1_hit_col;   //!
   TBranch        *b_Earm_BBPSTF1_hit_cell;   //!
   TBranch        *b_Earm_BBPSTF1_hit_plane;   //!
   TBranch        *b_Earm_BBPSTF1_hit_xcell;   //!
   TBranch        *b_Earm_BBPSTF1_hit_ycell;   //!
   TBranch        *b_Earm_BBPSTF1_hit_zcell;   //!
   TBranch        *b_Earm_BBPSTF1_hit_xcellg;   //!
   TBranch        *b_Earm_BBPSTF1_hit_ycellg;   //!
   TBranch        *b_Earm_BBPSTF1_hit_zcellg;   //!
   TBranch        *b_Earm_BBPSTF1_hit_xhit;   //!
   TBranch        *b_Earm_BBPSTF1_hit_yhit;   //!
   TBranch        *b_Earm_BBPSTF1_hit_zhit;   //!
   TBranch        *b_Earm_BBPSTF1_hit_sumedep;   //!
   TBranch        *b_Earm_BBPSTF1_hit_tavg;   //!
   TBranch        *b_Earm_BBPSTF1_hit_trms;   //!
   TBranch        *b_Earm_BBPSTF1_hit_tmin;   //!
   TBranch        *b_Earm_BBPSTF1_hit_tmax;   //!
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
   TBranch        *b_Earm_GC_PMT_Glass_hit_nhits;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_row;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_col;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_cell;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_plane;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_xcell;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_ycell;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_zcell;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_xcellg;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_ycellg;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_zcellg;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_xhit;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_yhit;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_zhit;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_sumedep;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_tavg;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_trms;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_tmin;   //!
   TBranch        *b_Earm_GC_PMT_Glass_hit_tmax;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_nhits;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_row;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_col;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_cell;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_plane;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_xcell;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_ycell;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_zcell;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_xcellg;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_ycellg;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_zcellg;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_xhit;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_yhit;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_zhit;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_sumedep;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_tavg;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_trms;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_tmin;   //!
   TBranch        *b_Earm_GEMElectronicsGMn_hit_tmax;   //!
   TBranch        *b_Earm_GRINCH_hit_nhits;   //!
   TBranch        *b_Earm_GRINCH_hit_PMT;   //!
   TBranch        *b_Earm_GRINCH_hit_row;   //!
   TBranch        *b_Earm_GRINCH_hit_col;   //!
   TBranch        *b_Earm_GRINCH_hit_xpmt;   //!
   TBranch        *b_Earm_GRINCH_hit_ypmt;   //!
   TBranch        *b_Earm_GRINCH_hit_zpmt;   //!
   TBranch        *b_Earm_GRINCH_hit_xgpmt;   //!
   TBranch        *b_Earm_GRINCH_hit_ygpmt;   //!
   TBranch        *b_Earm_GRINCH_hit_zgpmt;   //!
   TBranch        *b_Earm_GRINCH_hit_NumPhotoelectrons;   //!
   TBranch        *b_Earm_GRINCH_hit_Time_avg;   //!
   TBranch        *b_Earm_GRINCH_hit_Time_rms;   //!
   TBranch        *b_Earm_GRINCH_hit_Time_min;   //!
   TBranch        *b_Earm_GRINCH_hit_Time_max;   //!
   TBranch        *b_Earm_GRINCH_hit_mTrackNo;   //!
   TBranch        *b_Earm_GRINCH_hit_xhit;   //!
   TBranch        *b_Earm_GRINCH_hit_yhit;   //!
   TBranch        *b_Earm_GRINCH_hit_zhit;   //!
   TBranch        *b_Earm_GRINCH_hit_pxhit;   //!
   TBranch        *b_Earm_GRINCH_hit_pyhit;   //!
   TBranch        *b_Earm_GRINCH_hit_pzhit;   //!
   TBranch        *b_Earm_GRINCH_hit_pvx;   //!
   TBranch        *b_Earm_GRINCH_hit_pvy;   //!
   TBranch        *b_Earm_GRINCH_hit_pvz;   //!
   TBranch        *b_Earm_GRINCH_hit_ppx;   //!
   TBranch        *b_Earm_GRINCH_hit_ppy;   //!
   TBranch        *b_Earm_GRINCH_hit_ppz;   //!
   TBranch        *b_Earm_GRINCH_hit_volume_flag;   //!
   TBranch        *b_Harm_CDET_hit_nhits;   //!
   TBranch        *b_Harm_CDET_hit_PMT;   //!
   TBranch        *b_Harm_CDET_hit_row;   //!
   TBranch        *b_Harm_CDET_hit_col;   //!
   TBranch        *b_Harm_CDET_hit_plane;   //!
   TBranch        *b_Harm_CDET_hit_xcell;   //!
   TBranch        *b_Harm_CDET_hit_ycell;   //!
   TBranch        *b_Harm_CDET_hit_zcell;   //!
   TBranch        *b_Harm_CDET_hit_xgcell;   //!
   TBranch        *b_Harm_CDET_hit_ygcell;   //!
   TBranch        *b_Harm_CDET_hit_zgcell;   //!
   TBranch        *b_Harm_CDET_hit_NumPhotoelectrons;   //!
   TBranch        *b_Harm_CDET_hit_Time_avg;   //!
   TBranch        *b_Harm_CDET_hit_Time_rms;   //!
   TBranch        *b_Harm_CDET_hit_Time_min;   //!
   TBranch        *b_Harm_CDET_hit_Time_max;   //!
   TBranch        *b_Harm_CDET_Scint_hit_nhits;   //!
   TBranch        *b_Harm_CDET_Scint_hit_row;   //!
   TBranch        *b_Harm_CDET_Scint_hit_col;   //!
   TBranch        *b_Harm_CDET_Scint_hit_cell;   //!
   TBranch        *b_Harm_CDET_Scint_hit_plane;   //!
   TBranch        *b_Harm_CDET_Scint_hit_xcell;   //!
   TBranch        *b_Harm_CDET_Scint_hit_ycell;   //!
   TBranch        *b_Harm_CDET_Scint_hit_zcell;   //!
   TBranch        *b_Harm_CDET_Scint_hit_xcellg;   //!
   TBranch        *b_Harm_CDET_Scint_hit_ycellg;   //!
   TBranch        *b_Harm_CDET_Scint_hit_zcellg;   //!
   TBranch        *b_Harm_CDET_Scint_hit_xhit;   //!
   TBranch        *b_Harm_CDET_Scint_hit_yhit;   //!
   TBranch        *b_Harm_CDET_Scint_hit_zhit;   //!
   TBranch        *b_Harm_CDET_Scint_hit_sumedep;   //!
   TBranch        *b_Harm_CDET_Scint_hit_tavg;   //!
   TBranch        *b_Harm_CDET_Scint_hit_trms;   //!
   TBranch        *b_Harm_CDET_Scint_hit_tmin;   //!
   TBranch        *b_Harm_CDET_Scint_hit_tmax;   //!
   TBranch        *b_Harm_HCal_hit_nhits;   //!
   TBranch        *b_Harm_HCal_hit_PMT;   //!
   TBranch        *b_Harm_HCal_hit_row;   //!
   TBranch        *b_Harm_HCal_hit_col;   //!
   TBranch        *b_Harm_HCal_hit_plane;   //!
   TBranch        *b_Harm_HCal_hit_xcell;   //!
   TBranch        *b_Harm_HCal_hit_ycell;   //!
   TBranch        *b_Harm_HCal_hit_zcell;   //!
   TBranch        *b_Harm_HCal_hit_xgcell;   //!
   TBranch        *b_Harm_HCal_hit_ygcell;   //!
   TBranch        *b_Harm_HCal_hit_zgcell;   //!
   TBranch        *b_Harm_HCal_hit_NumPhotoelectrons;   //!
   TBranch        *b_Harm_HCal_hit_Time_avg;   //!
   TBranch        *b_Harm_HCal_hit_Time_rms;   //!
   TBranch        *b_Harm_HCal_hit_Time_min;   //!
   TBranch        *b_Harm_HCal_hit_Time_max;   //!
   TBranch        *b_Harm_HCalScint_hit_nhits;   //!
   TBranch        *b_Harm_HCalScint_hit_row;   //!
   TBranch        *b_Harm_HCalScint_hit_col;   //!
   TBranch        *b_Harm_HCalScint_hit_cell;   //!
   TBranch        *b_Harm_HCalScint_hit_plane;   //!
   TBranch        *b_Harm_HCalScint_hit_xcell;   //!
   TBranch        *b_Harm_HCalScint_hit_ycell;   //!
   TBranch        *b_Harm_HCalScint_hit_zcell;   //!
   TBranch        *b_Harm_HCalScint_hit_xcellg;   //!
   TBranch        *b_Harm_HCalScint_hit_ycellg;   //!
   TBranch        *b_Harm_HCalScint_hit_zcellg;   //!
   TBranch        *b_Harm_HCalScint_hit_xhit;   //!
   TBranch        *b_Harm_HCalScint_hit_yhit;   //!
   TBranch        *b_Harm_HCalScint_hit_zhit;   //!
   TBranch        *b_Harm_HCalScint_hit_sumedep;   //!
   TBranch        *b_Harm_HCalScint_hit_tavg;   //!
   TBranch        *b_Harm_HCalScint_hit_trms;   //!
   TBranch        *b_Harm_HCalScint_hit_tmin;   //!
   TBranch        *b_Harm_HCalScint_hit_tmax;   //!

   gmn_tree(TTree *tree=0);
   virtual ~gmn_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gmn_tree_cxx
gmn_tree::gmn_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("beam_bkgd_1695.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("beam_bkgd_1695.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

gmn_tree::~gmn_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gmn_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gmn_tree::LoadTree(Long64_t entry)
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

void gmn_tree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Earm_BBGEM_hit_plane = 0;
   Earm_BBGEM_hit_strip = 0;
   Earm_BBGEM_hit_x = 0;
   Earm_BBGEM_hit_y = 0;
   Earm_BBGEM_hit_z = 0;
   Earm_BBGEM_hit_polx = 0;
   Earm_BBGEM_hit_poly = 0;
   Earm_BBGEM_hit_polz = 0;
   Earm_BBGEM_hit_t = 0;
   Earm_BBGEM_hit_trms = 0;
   Earm_BBGEM_hit_tmin = 0;
   Earm_BBGEM_hit_tmax = 0;
   Earm_BBGEM_hit_tx = 0;
   Earm_BBGEM_hit_ty = 0;
   Earm_BBGEM_hit_txp = 0;
   Earm_BBGEM_hit_typ = 0;
   Earm_BBGEM_hit_xg = 0;
   Earm_BBGEM_hit_yg = 0;
   Earm_BBGEM_hit_zg = 0;
   Earm_BBGEM_hit_trid = 0;
   Earm_BBGEM_hit_mid = 0;
   Earm_BBGEM_hit_pid = 0;
   Earm_BBGEM_hit_vx = 0;
   Earm_BBGEM_hit_vy = 0;
   Earm_BBGEM_hit_vz = 0;
   Earm_BBGEM_hit_p = 0;
   Earm_BBGEM_hit_edep = 0;
   Earm_BBGEM_hit_beta = 0;
   Earm_BBGEM_Track_TID = 0;
   Earm_BBGEM_Track_PID = 0;
   Earm_BBGEM_Track_MID = 0;
   Earm_BBGEM_Track_NumHits = 0;
   Earm_BBGEM_Track_NumPlanes = 0;
   Earm_BBGEM_Track_NDF = 0;
   Earm_BBGEM_Track_Chi2fit = 0;
   Earm_BBGEM_Track_Chi2true = 0;
   Earm_BBGEM_Track_X = 0;
   Earm_BBGEM_Track_Y = 0;
   Earm_BBGEM_Track_Xp = 0;
   Earm_BBGEM_Track_Yp = 0;
   Earm_BBGEM_Track_T = 0;
   Earm_BBGEM_Track_P = 0;
   Earm_BBGEM_Track_Sx = 0;
   Earm_BBGEM_Track_Sy = 0;
   Earm_BBGEM_Track_Sz = 0;
   Earm_BBGEM_Track_Xfit = 0;
   Earm_BBGEM_Track_Yfit = 0;
   Earm_BBGEM_Track_Xpfit = 0;
   Earm_BBGEM_Track_Ypfit = 0;
   Earm_BBHodoScint_hit_row = 0;
   Earm_BBHodoScint_hit_col = 0;
   Earm_BBHodoScint_hit_cell = 0;
   Earm_BBHodoScint_hit_plane = 0;
   Earm_BBHodoScint_hit_xcell = 0;
   Earm_BBHodoScint_hit_ycell = 0;
   Earm_BBHodoScint_hit_zcell = 0;
   Earm_BBHodoScint_hit_xcellg = 0;
   Earm_BBHodoScint_hit_ycellg = 0;
   Earm_BBHodoScint_hit_zcellg = 0;
   Earm_BBHodoScint_hit_xhit = 0;
   Earm_BBHodoScint_hit_yhit = 0;
   Earm_BBHodoScint_hit_zhit = 0;
   Earm_BBHodoScint_hit_sumedep = 0;
   Earm_BBHodoScint_hit_tavg = 0;
   Earm_BBHodoScint_hit_trms = 0;
   Earm_BBHodoScint_hit_tmin = 0;
   Earm_BBHodoScint_hit_tmax = 0;
   Earm_BBPS_hit_PMT = 0;
   Earm_BBPS_hit_row = 0;
   Earm_BBPS_hit_col = 0;
   Earm_BBPS_hit_plane = 0;
   Earm_BBPS_hit_xcell = 0;
   Earm_BBPS_hit_ycell = 0;
   Earm_BBPS_hit_zcell = 0;
   Earm_BBPS_hit_xgcell = 0;
   Earm_BBPS_hit_ygcell = 0;
   Earm_BBPS_hit_zgcell = 0;
   Earm_BBPS_hit_NumPhotoelectrons = 0;
   Earm_BBPS_hit_Time_avg = 0;
   Earm_BBPS_hit_Time_rms = 0;
   Earm_BBPS_hit_Time_min = 0;
   Earm_BBPS_hit_Time_max = 0;
   Earm_BBPSTF1_hit_row = 0;
   Earm_BBPSTF1_hit_col = 0;
   Earm_BBPSTF1_hit_cell = 0;
   Earm_BBPSTF1_hit_plane = 0;
   Earm_BBPSTF1_hit_xcell = 0;
   Earm_BBPSTF1_hit_ycell = 0;
   Earm_BBPSTF1_hit_zcell = 0;
   Earm_BBPSTF1_hit_xcellg = 0;
   Earm_BBPSTF1_hit_ycellg = 0;
   Earm_BBPSTF1_hit_zcellg = 0;
   Earm_BBPSTF1_hit_xhit = 0;
   Earm_BBPSTF1_hit_yhit = 0;
   Earm_BBPSTF1_hit_zhit = 0;
   Earm_BBPSTF1_hit_sumedep = 0;
   Earm_BBPSTF1_hit_tavg = 0;
   Earm_BBPSTF1_hit_trms = 0;
   Earm_BBPSTF1_hit_tmin = 0;
   Earm_BBPSTF1_hit_tmax = 0;
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
   Earm_GC_PMT_Glass_hit_row = 0;
   Earm_GC_PMT_Glass_hit_col = 0;
   Earm_GC_PMT_Glass_hit_cell = 0;
   Earm_GC_PMT_Glass_hit_plane = 0;
   Earm_GC_PMT_Glass_hit_xcell = 0;
   Earm_GC_PMT_Glass_hit_ycell = 0;
   Earm_GC_PMT_Glass_hit_zcell = 0;
   Earm_GC_PMT_Glass_hit_xcellg = 0;
   Earm_GC_PMT_Glass_hit_ycellg = 0;
   Earm_GC_PMT_Glass_hit_zcellg = 0;
   Earm_GC_PMT_Glass_hit_xhit = 0;
   Earm_GC_PMT_Glass_hit_yhit = 0;
   Earm_GC_PMT_Glass_hit_zhit = 0;
   Earm_GC_PMT_Glass_hit_sumedep = 0;
   Earm_GC_PMT_Glass_hit_tavg = 0;
   Earm_GC_PMT_Glass_hit_trms = 0;
   Earm_GC_PMT_Glass_hit_tmin = 0;
   Earm_GC_PMT_Glass_hit_tmax = 0;
   Earm_GEMElectronicsGMn_hit_row = 0;
   Earm_GEMElectronicsGMn_hit_col = 0;
   Earm_GEMElectronicsGMn_hit_cell = 0;
   Earm_GEMElectronicsGMn_hit_plane = 0;
   Earm_GEMElectronicsGMn_hit_xcell = 0;
   Earm_GEMElectronicsGMn_hit_ycell = 0;
   Earm_GEMElectronicsGMn_hit_zcell = 0;
   Earm_GEMElectronicsGMn_hit_xcellg = 0;
   Earm_GEMElectronicsGMn_hit_ycellg = 0;
   Earm_GEMElectronicsGMn_hit_zcellg = 0;
   Earm_GEMElectronicsGMn_hit_xhit = 0;
   Earm_GEMElectronicsGMn_hit_yhit = 0;
   Earm_GEMElectronicsGMn_hit_zhit = 0;
   Earm_GEMElectronicsGMn_hit_sumedep = 0;
   Earm_GEMElectronicsGMn_hit_tavg = 0;
   Earm_GEMElectronicsGMn_hit_trms = 0;
   Earm_GEMElectronicsGMn_hit_tmin = 0;
   Earm_GEMElectronicsGMn_hit_tmax = 0;
   Earm_GRINCH_hit_PMT = 0;
   Earm_GRINCH_hit_row = 0;
   Earm_GRINCH_hit_col = 0;
   Earm_GRINCH_hit_xpmt = 0;
   Earm_GRINCH_hit_ypmt = 0;
   Earm_GRINCH_hit_zpmt = 0;
   Earm_GRINCH_hit_xgpmt = 0;
   Earm_GRINCH_hit_ygpmt = 0;
   Earm_GRINCH_hit_zgpmt = 0;
   Earm_GRINCH_hit_NumPhotoelectrons = 0;
   Earm_GRINCH_hit_Time_avg = 0;
   Earm_GRINCH_hit_Time_rms = 0;
   Earm_GRINCH_hit_Time_min = 0;
   Earm_GRINCH_hit_Time_max = 0;
   Earm_GRINCH_hit_mTrackNo = 0;
   Earm_GRINCH_hit_xhit = 0;
   Earm_GRINCH_hit_yhit = 0;
   Earm_GRINCH_hit_zhit = 0;
   Earm_GRINCH_hit_pxhit = 0;
   Earm_GRINCH_hit_pyhit = 0;
   Earm_GRINCH_hit_pzhit = 0;
   Earm_GRINCH_hit_pvx = 0;
   Earm_GRINCH_hit_pvy = 0;
   Earm_GRINCH_hit_pvz = 0;
   Earm_GRINCH_hit_ppx = 0;
   Earm_GRINCH_hit_ppy = 0;
   Earm_GRINCH_hit_ppz = 0;
   Earm_GRINCH_hit_volume_flag = 0;
   Harm_CDET_hit_PMT = 0;
   Harm_CDET_hit_row = 0;
   Harm_CDET_hit_col = 0;
   Harm_CDET_hit_plane = 0;
   Harm_CDET_hit_xcell = 0;
   Harm_CDET_hit_ycell = 0;
   Harm_CDET_hit_zcell = 0;
   Harm_CDET_hit_xgcell = 0;
   Harm_CDET_hit_ygcell = 0;
   Harm_CDET_hit_zgcell = 0;
   Harm_CDET_hit_NumPhotoelectrons = 0;
   Harm_CDET_hit_Time_avg = 0;
   Harm_CDET_hit_Time_rms = 0;
   Harm_CDET_hit_Time_min = 0;
   Harm_CDET_hit_Time_max = 0;
   Harm_CDET_Scint_hit_row = 0;
   Harm_CDET_Scint_hit_col = 0;
   Harm_CDET_Scint_hit_cell = 0;
   Harm_CDET_Scint_hit_plane = 0;
   Harm_CDET_Scint_hit_xcell = 0;
   Harm_CDET_Scint_hit_ycell = 0;
   Harm_CDET_Scint_hit_zcell = 0;
   Harm_CDET_Scint_hit_xcellg = 0;
   Harm_CDET_Scint_hit_ycellg = 0;
   Harm_CDET_Scint_hit_zcellg = 0;
   Harm_CDET_Scint_hit_xhit = 0;
   Harm_CDET_Scint_hit_yhit = 0;
   Harm_CDET_Scint_hit_zhit = 0;
   Harm_CDET_Scint_hit_sumedep = 0;
   Harm_CDET_Scint_hit_tavg = 0;
   Harm_CDET_Scint_hit_trms = 0;
   Harm_CDET_Scint_hit_tmin = 0;
   Harm_CDET_Scint_hit_tmax = 0;
   Harm_HCal_hit_PMT = 0;
   Harm_HCal_hit_row = 0;
   Harm_HCal_hit_col = 0;
   Harm_HCal_hit_plane = 0;
   Harm_HCal_hit_xcell = 0;
   Harm_HCal_hit_ycell = 0;
   Harm_HCal_hit_zcell = 0;
   Harm_HCal_hit_xgcell = 0;
   Harm_HCal_hit_ygcell = 0;
   Harm_HCal_hit_zgcell = 0;
   Harm_HCal_hit_NumPhotoelectrons = 0;
   Harm_HCal_hit_Time_avg = 0;
   Harm_HCal_hit_Time_rms = 0;
   Harm_HCal_hit_Time_min = 0;
   Harm_HCal_hit_Time_max = 0;
   Harm_HCalScint_hit_row = 0;
   Harm_HCalScint_hit_col = 0;
   Harm_HCalScint_hit_cell = 0;
   Harm_HCalScint_hit_plane = 0;
   Harm_HCalScint_hit_xcell = 0;
   Harm_HCalScint_hit_ycell = 0;
   Harm_HCalScint_hit_zcell = 0;
   Harm_HCalScint_hit_xcellg = 0;
   Harm_HCalScint_hit_ycellg = 0;
   Harm_HCalScint_hit_zcellg = 0;
   Harm_HCalScint_hit_xhit = 0;
   Harm_HCalScint_hit_yhit = 0;
   Harm_HCalScint_hit_zhit = 0;
   Harm_HCalScint_hit_sumedep = 0;
   Harm_HCalScint_hit_tavg = 0;
   Harm_HCalScint_hit_trms = 0;
   Harm_HCalScint_hit_tmin = 0;
   Harm_HCalScint_hit_tmax = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev", &ev_count, &b_ev);
   fChain->SetBranchAddress("gen", &gen_thbb, &b_gen);
   fChain->SetBranchAddress("Earm.BBGEM.hit.nhits", &Earm_BBGEM_hit_nhits, &b_Earm_BBGEM_hit_nhits);
   fChain->SetBranchAddress("Earm.BBGEM.hit.plane", &Earm_BBGEM_hit_plane, &b_Earm_BBGEM_hit_plane);
   fChain->SetBranchAddress("Earm.BBGEM.hit.strip", &Earm_BBGEM_hit_strip, &b_Earm_BBGEM_hit_strip);
   fChain->SetBranchAddress("Earm.BBGEM.hit.x", &Earm_BBGEM_hit_x, &b_Earm_BBGEM_hit_x);
   fChain->SetBranchAddress("Earm.BBGEM.hit.y", &Earm_BBGEM_hit_y, &b_Earm_BBGEM_hit_y);
   fChain->SetBranchAddress("Earm.BBGEM.hit.z", &Earm_BBGEM_hit_z, &b_Earm_BBGEM_hit_z);
   fChain->SetBranchAddress("Earm.BBGEM.hit.polx", &Earm_BBGEM_hit_polx, &b_Earm_BBGEM_hit_polx);
   fChain->SetBranchAddress("Earm.BBGEM.hit.poly", &Earm_BBGEM_hit_poly, &b_Earm_BBGEM_hit_poly);
   fChain->SetBranchAddress("Earm.BBGEM.hit.polz", &Earm_BBGEM_hit_polz, &b_Earm_BBGEM_hit_polz);
   fChain->SetBranchAddress("Earm.BBGEM.hit.t", &Earm_BBGEM_hit_t, &b_Earm_BBGEM_hit_t);
   fChain->SetBranchAddress("Earm.BBGEM.hit.trms", &Earm_BBGEM_hit_trms, &b_Earm_BBGEM_hit_trms);
   fChain->SetBranchAddress("Earm.BBGEM.hit.tmin", &Earm_BBGEM_hit_tmin, &b_Earm_BBGEM_hit_tmin);
   fChain->SetBranchAddress("Earm.BBGEM.hit.tmax", &Earm_BBGEM_hit_tmax, &b_Earm_BBGEM_hit_tmax);
   fChain->SetBranchAddress("Earm.BBGEM.hit.tx", &Earm_BBGEM_hit_tx, &b_Earm_BBGEM_hit_tx);
   fChain->SetBranchAddress("Earm.BBGEM.hit.ty", &Earm_BBGEM_hit_ty, &b_Earm_BBGEM_hit_ty);
   fChain->SetBranchAddress("Earm.BBGEM.hit.txp", &Earm_BBGEM_hit_txp, &b_Earm_BBGEM_hit_txp);
   fChain->SetBranchAddress("Earm.BBGEM.hit.typ", &Earm_BBGEM_hit_typ, &b_Earm_BBGEM_hit_typ);
   fChain->SetBranchAddress("Earm.BBGEM.hit.xg", &Earm_BBGEM_hit_xg, &b_Earm_BBGEM_hit_xg);
   fChain->SetBranchAddress("Earm.BBGEM.hit.yg", &Earm_BBGEM_hit_yg, &b_Earm_BBGEM_hit_yg);
   fChain->SetBranchAddress("Earm.BBGEM.hit.zg", &Earm_BBGEM_hit_zg, &b_Earm_BBGEM_hit_zg);
   fChain->SetBranchAddress("Earm.BBGEM.hit.trid", &Earm_BBGEM_hit_trid, &b_Earm_BBGEM_hit_trid);
   fChain->SetBranchAddress("Earm.BBGEM.hit.mid", &Earm_BBGEM_hit_mid, &b_Earm_BBGEM_hit_mid);
   fChain->SetBranchAddress("Earm.BBGEM.hit.pid", &Earm_BBGEM_hit_pid, &b_Earm_BBGEM_hit_pid);
   fChain->SetBranchAddress("Earm.BBGEM.hit.vx", &Earm_BBGEM_hit_vx, &b_Earm_BBGEM_hit_vx);
   fChain->SetBranchAddress("Earm.BBGEM.hit.vy", &Earm_BBGEM_hit_vy, &b_Earm_BBGEM_hit_vy);
   fChain->SetBranchAddress("Earm.BBGEM.hit.vz", &Earm_BBGEM_hit_vz, &b_Earm_BBGEM_hit_vz);
   fChain->SetBranchAddress("Earm.BBGEM.hit.p", &Earm_BBGEM_hit_p, &b_Earm_BBGEM_hit_p);
   fChain->SetBranchAddress("Earm.BBGEM.hit.edep", &Earm_BBGEM_hit_edep, &b_Earm_BBGEM_hit_edep);
   fChain->SetBranchAddress("Earm.BBGEM.hit.beta", &Earm_BBGEM_hit_beta, &b_Earm_BBGEM_hit_beta);
   fChain->SetBranchAddress("Earm.BBGEM.Track.ntracks", &Earm_BBGEM_Track_ntracks, &b_Earm_BBGEM_Track_ntracks);
   fChain->SetBranchAddress("Earm.BBGEM.Track.TID", &Earm_BBGEM_Track_TID, &b_Earm_BBGEM_Track_TID);
   fChain->SetBranchAddress("Earm.BBGEM.Track.PID", &Earm_BBGEM_Track_PID, &b_Earm_BBGEM_Track_PID);
   fChain->SetBranchAddress("Earm.BBGEM.Track.MID", &Earm_BBGEM_Track_MID, &b_Earm_BBGEM_Track_MID);
   fChain->SetBranchAddress("Earm.BBGEM.Track.NumHits", &Earm_BBGEM_Track_NumHits, &b_Earm_BBGEM_Track_NumHits);
   fChain->SetBranchAddress("Earm.BBGEM.Track.NumPlanes", &Earm_BBGEM_Track_NumPlanes, &b_Earm_BBGEM_Track_NumPlanes);
   fChain->SetBranchAddress("Earm.BBGEM.Track.NDF", &Earm_BBGEM_Track_NDF, &b_Earm_BBGEM_Track_NDF);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Chi2fit", &Earm_BBGEM_Track_Chi2fit, &b_Earm_BBGEM_Track_Chi2fit);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Chi2true", &Earm_BBGEM_Track_Chi2true, &b_Earm_BBGEM_Track_Chi2true);
   fChain->SetBranchAddress("Earm.BBGEM.Track.X", &Earm_BBGEM_Track_X, &b_Earm_BBGEM_Track_X);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Y", &Earm_BBGEM_Track_Y, &b_Earm_BBGEM_Track_Y);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Xp", &Earm_BBGEM_Track_Xp, &b_Earm_BBGEM_Track_Xp);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Yp", &Earm_BBGEM_Track_Yp, &b_Earm_BBGEM_Track_Yp);
   fChain->SetBranchAddress("Earm.BBGEM.Track.T", &Earm_BBGEM_Track_T, &b_Earm_BBGEM_Track_T);
   fChain->SetBranchAddress("Earm.BBGEM.Track.P", &Earm_BBGEM_Track_P, &b_Earm_BBGEM_Track_P);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Sx", &Earm_BBGEM_Track_Sx, &b_Earm_BBGEM_Track_Sx);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Sy", &Earm_BBGEM_Track_Sy, &b_Earm_BBGEM_Track_Sy);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Sz", &Earm_BBGEM_Track_Sz, &b_Earm_BBGEM_Track_Sz);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Xfit", &Earm_BBGEM_Track_Xfit, &b_Earm_BBGEM_Track_Xfit);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Yfit", &Earm_BBGEM_Track_Yfit, &b_Earm_BBGEM_Track_Yfit);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Xpfit", &Earm_BBGEM_Track_Xpfit, &b_Earm_BBGEM_Track_Xpfit);
   fChain->SetBranchAddress("Earm.BBGEM.Track.Ypfit", &Earm_BBGEM_Track_Ypfit, &b_Earm_BBGEM_Track_Ypfit);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.nhits", &Earm_BBHodoScint_hit_nhits, &b_Earm_BBHodoScint_hit_nhits);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.row", &Earm_BBHodoScint_hit_row, &b_Earm_BBHodoScint_hit_row);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.col", &Earm_BBHodoScint_hit_col, &b_Earm_BBHodoScint_hit_col);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.cell", &Earm_BBHodoScint_hit_cell, &b_Earm_BBHodoScint_hit_cell);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.plane", &Earm_BBHodoScint_hit_plane, &b_Earm_BBHodoScint_hit_plane);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.xcell", &Earm_BBHodoScint_hit_xcell, &b_Earm_BBHodoScint_hit_xcell);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.ycell", &Earm_BBHodoScint_hit_ycell, &b_Earm_BBHodoScint_hit_ycell);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.zcell", &Earm_BBHodoScint_hit_zcell, &b_Earm_BBHodoScint_hit_zcell);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.xcellg", &Earm_BBHodoScint_hit_xcellg, &b_Earm_BBHodoScint_hit_xcellg);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.ycellg", &Earm_BBHodoScint_hit_ycellg, &b_Earm_BBHodoScint_hit_ycellg);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.zcellg", &Earm_BBHodoScint_hit_zcellg, &b_Earm_BBHodoScint_hit_zcellg);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.xhit", &Earm_BBHodoScint_hit_xhit, &b_Earm_BBHodoScint_hit_xhit);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.yhit", &Earm_BBHodoScint_hit_yhit, &b_Earm_BBHodoScint_hit_yhit);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.zhit", &Earm_BBHodoScint_hit_zhit, &b_Earm_BBHodoScint_hit_zhit);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.sumedep", &Earm_BBHodoScint_hit_sumedep, &b_Earm_BBHodoScint_hit_sumedep);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.tavg", &Earm_BBHodoScint_hit_tavg, &b_Earm_BBHodoScint_hit_tavg);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.trms", &Earm_BBHodoScint_hit_trms, &b_Earm_BBHodoScint_hit_trms);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.tmin", &Earm_BBHodoScint_hit_tmin, &b_Earm_BBHodoScint_hit_tmin);
   fChain->SetBranchAddress("Earm.BBHodoScint.hit.tmax", &Earm_BBHodoScint_hit_tmax, &b_Earm_BBHodoScint_hit_tmax);
   fChain->SetBranchAddress("Earm.BBPS.hit.nhits", &Earm_BBPS_hit_nhits, &b_Earm_BBPS_hit_nhits);
   fChain->SetBranchAddress("Earm.BBPS.hit.PMT", &Earm_BBPS_hit_PMT, &b_Earm_BBPS_hit_PMT);
   fChain->SetBranchAddress("Earm.BBPS.hit.row", &Earm_BBPS_hit_row, &b_Earm_BBPS_hit_row);
   fChain->SetBranchAddress("Earm.BBPS.hit.col", &Earm_BBPS_hit_col, &b_Earm_BBPS_hit_col);
   fChain->SetBranchAddress("Earm.BBPS.hit.plane", &Earm_BBPS_hit_plane, &b_Earm_BBPS_hit_plane);
   fChain->SetBranchAddress("Earm.BBPS.hit.xcell", &Earm_BBPS_hit_xcell, &b_Earm_BBPS_hit_xcell);
   fChain->SetBranchAddress("Earm.BBPS.hit.ycell", &Earm_BBPS_hit_ycell, &b_Earm_BBPS_hit_ycell);
   fChain->SetBranchAddress("Earm.BBPS.hit.zcell", &Earm_BBPS_hit_zcell, &b_Earm_BBPS_hit_zcell);
   fChain->SetBranchAddress("Earm.BBPS.hit.xgcell", &Earm_BBPS_hit_xgcell, &b_Earm_BBPS_hit_xgcell);
   fChain->SetBranchAddress("Earm.BBPS.hit.ygcell", &Earm_BBPS_hit_ygcell, &b_Earm_BBPS_hit_ygcell);
   fChain->SetBranchAddress("Earm.BBPS.hit.zgcell", &Earm_BBPS_hit_zgcell, &b_Earm_BBPS_hit_zgcell);
   fChain->SetBranchAddress("Earm.BBPS.hit.NumPhotoelectrons", &Earm_BBPS_hit_NumPhotoelectrons, &b_Earm_BBPS_hit_NumPhotoelectrons);
   fChain->SetBranchAddress("Earm.BBPS.hit.Time_avg", &Earm_BBPS_hit_Time_avg, &b_Earm_BBPS_hit_Time_avg);
   fChain->SetBranchAddress("Earm.BBPS.hit.Time_rms", &Earm_BBPS_hit_Time_rms, &b_Earm_BBPS_hit_Time_rms);
   fChain->SetBranchAddress("Earm.BBPS.hit.Time_min", &Earm_BBPS_hit_Time_min, &b_Earm_BBPS_hit_Time_min);
   fChain->SetBranchAddress("Earm.BBPS.hit.Time_max", &Earm_BBPS_hit_Time_max, &b_Earm_BBPS_hit_Time_max);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.nhits", &Earm_BBPSTF1_hit_nhits, &b_Earm_BBPSTF1_hit_nhits);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.row", &Earm_BBPSTF1_hit_row, &b_Earm_BBPSTF1_hit_row);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.col", &Earm_BBPSTF1_hit_col, &b_Earm_BBPSTF1_hit_col);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.cell", &Earm_BBPSTF1_hit_cell, &b_Earm_BBPSTF1_hit_cell);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.plane", &Earm_BBPSTF1_hit_plane, &b_Earm_BBPSTF1_hit_plane);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.xcell", &Earm_BBPSTF1_hit_xcell, &b_Earm_BBPSTF1_hit_xcell);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.ycell", &Earm_BBPSTF1_hit_ycell, &b_Earm_BBPSTF1_hit_ycell);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.zcell", &Earm_BBPSTF1_hit_zcell, &b_Earm_BBPSTF1_hit_zcell);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.xcellg", &Earm_BBPSTF1_hit_xcellg, &b_Earm_BBPSTF1_hit_xcellg);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.ycellg", &Earm_BBPSTF1_hit_ycellg, &b_Earm_BBPSTF1_hit_ycellg);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.zcellg", &Earm_BBPSTF1_hit_zcellg, &b_Earm_BBPSTF1_hit_zcellg);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.xhit", &Earm_BBPSTF1_hit_xhit, &b_Earm_BBPSTF1_hit_xhit);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.yhit", &Earm_BBPSTF1_hit_yhit, &b_Earm_BBPSTF1_hit_yhit);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.zhit", &Earm_BBPSTF1_hit_zhit, &b_Earm_BBPSTF1_hit_zhit);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.sumedep", &Earm_BBPSTF1_hit_sumedep, &b_Earm_BBPSTF1_hit_sumedep);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.tavg", &Earm_BBPSTF1_hit_tavg, &b_Earm_BBPSTF1_hit_tavg);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.trms", &Earm_BBPSTF1_hit_trms, &b_Earm_BBPSTF1_hit_trms);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.tmin", &Earm_BBPSTF1_hit_tmin, &b_Earm_BBPSTF1_hit_tmin);
   fChain->SetBranchAddress("Earm.BBPSTF1.hit.tmax", &Earm_BBPSTF1_hit_tmax, &b_Earm_BBPSTF1_hit_tmax);
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
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.nhits", &Earm_GC_PMT_Glass_hit_nhits, &b_Earm_GC_PMT_Glass_hit_nhits);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.row", &Earm_GC_PMT_Glass_hit_row, &b_Earm_GC_PMT_Glass_hit_row);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.col", &Earm_GC_PMT_Glass_hit_col, &b_Earm_GC_PMT_Glass_hit_col);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.cell", &Earm_GC_PMT_Glass_hit_cell, &b_Earm_GC_PMT_Glass_hit_cell);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.plane", &Earm_GC_PMT_Glass_hit_plane, &b_Earm_GC_PMT_Glass_hit_plane);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.xcell", &Earm_GC_PMT_Glass_hit_xcell, &b_Earm_GC_PMT_Glass_hit_xcell);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.ycell", &Earm_GC_PMT_Glass_hit_ycell, &b_Earm_GC_PMT_Glass_hit_ycell);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.zcell", &Earm_GC_PMT_Glass_hit_zcell, &b_Earm_GC_PMT_Glass_hit_zcell);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.xcellg", &Earm_GC_PMT_Glass_hit_xcellg, &b_Earm_GC_PMT_Glass_hit_xcellg);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.ycellg", &Earm_GC_PMT_Glass_hit_ycellg, &b_Earm_GC_PMT_Glass_hit_ycellg);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.zcellg", &Earm_GC_PMT_Glass_hit_zcellg, &b_Earm_GC_PMT_Glass_hit_zcellg);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.xhit", &Earm_GC_PMT_Glass_hit_xhit, &b_Earm_GC_PMT_Glass_hit_xhit);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.yhit", &Earm_GC_PMT_Glass_hit_yhit, &b_Earm_GC_PMT_Glass_hit_yhit);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.zhit", &Earm_GC_PMT_Glass_hit_zhit, &b_Earm_GC_PMT_Glass_hit_zhit);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.sumedep", &Earm_GC_PMT_Glass_hit_sumedep, &b_Earm_GC_PMT_Glass_hit_sumedep);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.tavg", &Earm_GC_PMT_Glass_hit_tavg, &b_Earm_GC_PMT_Glass_hit_tavg);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.trms", &Earm_GC_PMT_Glass_hit_trms, &b_Earm_GC_PMT_Glass_hit_trms);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.tmin", &Earm_GC_PMT_Glass_hit_tmin, &b_Earm_GC_PMT_Glass_hit_tmin);
   fChain->SetBranchAddress("Earm.GC_PMT_Glass.hit.tmax", &Earm_GC_PMT_Glass_hit_tmax, &b_Earm_GC_PMT_Glass_hit_tmax);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.nhits", &Earm_GEMElectronicsGMn_hit_nhits, &b_Earm_GEMElectronicsGMn_hit_nhits);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.row", &Earm_GEMElectronicsGMn_hit_row, &b_Earm_GEMElectronicsGMn_hit_row);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.col", &Earm_GEMElectronicsGMn_hit_col, &b_Earm_GEMElectronicsGMn_hit_col);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.cell", &Earm_GEMElectronicsGMn_hit_cell, &b_Earm_GEMElectronicsGMn_hit_cell);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.plane", &Earm_GEMElectronicsGMn_hit_plane, &b_Earm_GEMElectronicsGMn_hit_plane);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.xcell", &Earm_GEMElectronicsGMn_hit_xcell, &b_Earm_GEMElectronicsGMn_hit_xcell);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.ycell", &Earm_GEMElectronicsGMn_hit_ycell, &b_Earm_GEMElectronicsGMn_hit_ycell);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.zcell", &Earm_GEMElectronicsGMn_hit_zcell, &b_Earm_GEMElectronicsGMn_hit_zcell);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.xcellg", &Earm_GEMElectronicsGMn_hit_xcellg, &b_Earm_GEMElectronicsGMn_hit_xcellg);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.ycellg", &Earm_GEMElectronicsGMn_hit_ycellg, &b_Earm_GEMElectronicsGMn_hit_ycellg);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.zcellg", &Earm_GEMElectronicsGMn_hit_zcellg, &b_Earm_GEMElectronicsGMn_hit_zcellg);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.xhit", &Earm_GEMElectronicsGMn_hit_xhit, &b_Earm_GEMElectronicsGMn_hit_xhit);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.yhit", &Earm_GEMElectronicsGMn_hit_yhit, &b_Earm_GEMElectronicsGMn_hit_yhit);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.zhit", &Earm_GEMElectronicsGMn_hit_zhit, &b_Earm_GEMElectronicsGMn_hit_zhit);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.sumedep", &Earm_GEMElectronicsGMn_hit_sumedep, &b_Earm_GEMElectronicsGMn_hit_sumedep);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.tavg", &Earm_GEMElectronicsGMn_hit_tavg, &b_Earm_GEMElectronicsGMn_hit_tavg);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.trms", &Earm_GEMElectronicsGMn_hit_trms, &b_Earm_GEMElectronicsGMn_hit_trms);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.tmin", &Earm_GEMElectronicsGMn_hit_tmin, &b_Earm_GEMElectronicsGMn_hit_tmin);
   fChain->SetBranchAddress("Earm.GEMElectronicsGMn.hit.tmax", &Earm_GEMElectronicsGMn_hit_tmax, &b_Earm_GEMElectronicsGMn_hit_tmax);
   fChain->SetBranchAddress("Earm.GRINCH.hit.nhits", &Earm_GRINCH_hit_nhits, &b_Earm_GRINCH_hit_nhits);
   fChain->SetBranchAddress("Earm.GRINCH.hit.PMT", &Earm_GRINCH_hit_PMT, &b_Earm_GRINCH_hit_PMT);
   fChain->SetBranchAddress("Earm.GRINCH.hit.row", &Earm_GRINCH_hit_row, &b_Earm_GRINCH_hit_row);
   fChain->SetBranchAddress("Earm.GRINCH.hit.col", &Earm_GRINCH_hit_col, &b_Earm_GRINCH_hit_col);
   fChain->SetBranchAddress("Earm.GRINCH.hit.xpmt", &Earm_GRINCH_hit_xpmt, &b_Earm_GRINCH_hit_xpmt);
   fChain->SetBranchAddress("Earm.GRINCH.hit.ypmt", &Earm_GRINCH_hit_ypmt, &b_Earm_GRINCH_hit_ypmt);
   fChain->SetBranchAddress("Earm.GRINCH.hit.zpmt", &Earm_GRINCH_hit_zpmt, &b_Earm_GRINCH_hit_zpmt);
   fChain->SetBranchAddress("Earm.GRINCH.hit.xgpmt", &Earm_GRINCH_hit_xgpmt, &b_Earm_GRINCH_hit_xgpmt);
   fChain->SetBranchAddress("Earm.GRINCH.hit.ygpmt", &Earm_GRINCH_hit_ygpmt, &b_Earm_GRINCH_hit_ygpmt);
   fChain->SetBranchAddress("Earm.GRINCH.hit.zgpmt", &Earm_GRINCH_hit_zgpmt, &b_Earm_GRINCH_hit_zgpmt);
   fChain->SetBranchAddress("Earm.GRINCH.hit.NumPhotoelectrons", &Earm_GRINCH_hit_NumPhotoelectrons, &b_Earm_GRINCH_hit_NumPhotoelectrons);
   fChain->SetBranchAddress("Earm.GRINCH.hit.Time_avg", &Earm_GRINCH_hit_Time_avg, &b_Earm_GRINCH_hit_Time_avg);
   fChain->SetBranchAddress("Earm.GRINCH.hit.Time_rms", &Earm_GRINCH_hit_Time_rms, &b_Earm_GRINCH_hit_Time_rms);
   fChain->SetBranchAddress("Earm.GRINCH.hit.Time_min", &Earm_GRINCH_hit_Time_min, &b_Earm_GRINCH_hit_Time_min);
   fChain->SetBranchAddress("Earm.GRINCH.hit.Time_max", &Earm_GRINCH_hit_Time_max, &b_Earm_GRINCH_hit_Time_max);
   fChain->SetBranchAddress("Earm.GRINCH.hit.mTrackNo", &Earm_GRINCH_hit_mTrackNo, &b_Earm_GRINCH_hit_mTrackNo);
   fChain->SetBranchAddress("Earm.GRINCH.hit.xhit", &Earm_GRINCH_hit_xhit, &b_Earm_GRINCH_hit_xhit);
   fChain->SetBranchAddress("Earm.GRINCH.hit.yhit", &Earm_GRINCH_hit_yhit, &b_Earm_GRINCH_hit_yhit);
   fChain->SetBranchAddress("Earm.GRINCH.hit.zhit", &Earm_GRINCH_hit_zhit, &b_Earm_GRINCH_hit_zhit);
   fChain->SetBranchAddress("Earm.GRINCH.hit.pxhit", &Earm_GRINCH_hit_pxhit, &b_Earm_GRINCH_hit_pxhit);
   fChain->SetBranchAddress("Earm.GRINCH.hit.pyhit", &Earm_GRINCH_hit_pyhit, &b_Earm_GRINCH_hit_pyhit);
   fChain->SetBranchAddress("Earm.GRINCH.hit.pzhit", &Earm_GRINCH_hit_pzhit, &b_Earm_GRINCH_hit_pzhit);
   fChain->SetBranchAddress("Earm.GRINCH.hit.pvx", &Earm_GRINCH_hit_pvx, &b_Earm_GRINCH_hit_pvx);
   fChain->SetBranchAddress("Earm.GRINCH.hit.pvy", &Earm_GRINCH_hit_pvy, &b_Earm_GRINCH_hit_pvy);
   fChain->SetBranchAddress("Earm.GRINCH.hit.pvz", &Earm_GRINCH_hit_pvz, &b_Earm_GRINCH_hit_pvz);
   fChain->SetBranchAddress("Earm.GRINCH.hit.ppx", &Earm_GRINCH_hit_ppx, &b_Earm_GRINCH_hit_ppx);
   fChain->SetBranchAddress("Earm.GRINCH.hit.ppy", &Earm_GRINCH_hit_ppy, &b_Earm_GRINCH_hit_ppy);
   fChain->SetBranchAddress("Earm.GRINCH.hit.ppz", &Earm_GRINCH_hit_ppz, &b_Earm_GRINCH_hit_ppz);
   fChain->SetBranchAddress("Earm.GRINCH.hit.volume_flag", &Earm_GRINCH_hit_volume_flag, &b_Earm_GRINCH_hit_volume_flag);
   fChain->SetBranchAddress("Harm.CDET.hit.nhits", &Harm_CDET_hit_nhits, &b_Harm_CDET_hit_nhits);
   fChain->SetBranchAddress("Harm.CDET.hit.PMT", &Harm_CDET_hit_PMT, &b_Harm_CDET_hit_PMT);
   fChain->SetBranchAddress("Harm.CDET.hit.row", &Harm_CDET_hit_row, &b_Harm_CDET_hit_row);
   fChain->SetBranchAddress("Harm.CDET.hit.col", &Harm_CDET_hit_col, &b_Harm_CDET_hit_col);
   fChain->SetBranchAddress("Harm.CDET.hit.plane", &Harm_CDET_hit_plane, &b_Harm_CDET_hit_plane);
   fChain->SetBranchAddress("Harm.CDET.hit.xcell", &Harm_CDET_hit_xcell, &b_Harm_CDET_hit_xcell);
   fChain->SetBranchAddress("Harm.CDET.hit.ycell", &Harm_CDET_hit_ycell, &b_Harm_CDET_hit_ycell);
   fChain->SetBranchAddress("Harm.CDET.hit.zcell", &Harm_CDET_hit_zcell, &b_Harm_CDET_hit_zcell);
   fChain->SetBranchAddress("Harm.CDET.hit.xgcell", &Harm_CDET_hit_xgcell, &b_Harm_CDET_hit_xgcell);
   fChain->SetBranchAddress("Harm.CDET.hit.ygcell", &Harm_CDET_hit_ygcell, &b_Harm_CDET_hit_ygcell);
   fChain->SetBranchAddress("Harm.CDET.hit.zgcell", &Harm_CDET_hit_zgcell, &b_Harm_CDET_hit_zgcell);
   fChain->SetBranchAddress("Harm.CDET.hit.NumPhotoelectrons", &Harm_CDET_hit_NumPhotoelectrons, &b_Harm_CDET_hit_NumPhotoelectrons);
   fChain->SetBranchAddress("Harm.CDET.hit.Time_avg", &Harm_CDET_hit_Time_avg, &b_Harm_CDET_hit_Time_avg);
   fChain->SetBranchAddress("Harm.CDET.hit.Time_rms", &Harm_CDET_hit_Time_rms, &b_Harm_CDET_hit_Time_rms);
   fChain->SetBranchAddress("Harm.CDET.hit.Time_min", &Harm_CDET_hit_Time_min, &b_Harm_CDET_hit_Time_min);
   fChain->SetBranchAddress("Harm.CDET.hit.Time_max", &Harm_CDET_hit_Time_max, &b_Harm_CDET_hit_Time_max);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.nhits", &Harm_CDET_Scint_hit_nhits, &b_Harm_CDET_Scint_hit_nhits);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.row", &Harm_CDET_Scint_hit_row, &b_Harm_CDET_Scint_hit_row);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.col", &Harm_CDET_Scint_hit_col, &b_Harm_CDET_Scint_hit_col);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.cell", &Harm_CDET_Scint_hit_cell, &b_Harm_CDET_Scint_hit_cell);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.plane", &Harm_CDET_Scint_hit_plane, &b_Harm_CDET_Scint_hit_plane);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.xcell", &Harm_CDET_Scint_hit_xcell, &b_Harm_CDET_Scint_hit_xcell);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.ycell", &Harm_CDET_Scint_hit_ycell, &b_Harm_CDET_Scint_hit_ycell);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.zcell", &Harm_CDET_Scint_hit_zcell, &b_Harm_CDET_Scint_hit_zcell);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.xcellg", &Harm_CDET_Scint_hit_xcellg, &b_Harm_CDET_Scint_hit_xcellg);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.ycellg", &Harm_CDET_Scint_hit_ycellg, &b_Harm_CDET_Scint_hit_ycellg);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.zcellg", &Harm_CDET_Scint_hit_zcellg, &b_Harm_CDET_Scint_hit_zcellg);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.xhit", &Harm_CDET_Scint_hit_xhit, &b_Harm_CDET_Scint_hit_xhit);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.yhit", &Harm_CDET_Scint_hit_yhit, &b_Harm_CDET_Scint_hit_yhit);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.zhit", &Harm_CDET_Scint_hit_zhit, &b_Harm_CDET_Scint_hit_zhit);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.sumedep", &Harm_CDET_Scint_hit_sumedep, &b_Harm_CDET_Scint_hit_sumedep);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.tavg", &Harm_CDET_Scint_hit_tavg, &b_Harm_CDET_Scint_hit_tavg);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.trms", &Harm_CDET_Scint_hit_trms, &b_Harm_CDET_Scint_hit_trms);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.tmin", &Harm_CDET_Scint_hit_tmin, &b_Harm_CDET_Scint_hit_tmin);
   fChain->SetBranchAddress("Harm.CDET_Scint.hit.tmax", &Harm_CDET_Scint_hit_tmax, &b_Harm_CDET_Scint_hit_tmax);
   fChain->SetBranchAddress("Harm.HCal.hit.nhits", &Harm_HCal_hit_nhits, &b_Harm_HCal_hit_nhits);
   fChain->SetBranchAddress("Harm.HCal.hit.PMT", &Harm_HCal_hit_PMT, &b_Harm_HCal_hit_PMT);
   fChain->SetBranchAddress("Harm.HCal.hit.row", &Harm_HCal_hit_row, &b_Harm_HCal_hit_row);
   fChain->SetBranchAddress("Harm.HCal.hit.col", &Harm_HCal_hit_col, &b_Harm_HCal_hit_col);
   fChain->SetBranchAddress("Harm.HCal.hit.plane", &Harm_HCal_hit_plane, &b_Harm_HCal_hit_plane);
   fChain->SetBranchAddress("Harm.HCal.hit.xcell", &Harm_HCal_hit_xcell, &b_Harm_HCal_hit_xcell);
   fChain->SetBranchAddress("Harm.HCal.hit.ycell", &Harm_HCal_hit_ycell, &b_Harm_HCal_hit_ycell);
   fChain->SetBranchAddress("Harm.HCal.hit.zcell", &Harm_HCal_hit_zcell, &b_Harm_HCal_hit_zcell);
   fChain->SetBranchAddress("Harm.HCal.hit.xgcell", &Harm_HCal_hit_xgcell, &b_Harm_HCal_hit_xgcell);
   fChain->SetBranchAddress("Harm.HCal.hit.ygcell", &Harm_HCal_hit_ygcell, &b_Harm_HCal_hit_ygcell);
   fChain->SetBranchAddress("Harm.HCal.hit.zgcell", &Harm_HCal_hit_zgcell, &b_Harm_HCal_hit_zgcell);
   fChain->SetBranchAddress("Harm.HCal.hit.NumPhotoelectrons", &Harm_HCal_hit_NumPhotoelectrons, &b_Harm_HCal_hit_NumPhotoelectrons);
   fChain->SetBranchAddress("Harm.HCal.hit.Time_avg", &Harm_HCal_hit_Time_avg, &b_Harm_HCal_hit_Time_avg);
   fChain->SetBranchAddress("Harm.HCal.hit.Time_rms", &Harm_HCal_hit_Time_rms, &b_Harm_HCal_hit_Time_rms);
   fChain->SetBranchAddress("Harm.HCal.hit.Time_min", &Harm_HCal_hit_Time_min, &b_Harm_HCal_hit_Time_min);
   fChain->SetBranchAddress("Harm.HCal.hit.Time_max", &Harm_HCal_hit_Time_max, &b_Harm_HCal_hit_Time_max);
   fChain->SetBranchAddress("Harm.HCalScint.hit.nhits", &Harm_HCalScint_hit_nhits, &b_Harm_HCalScint_hit_nhits);
   fChain->SetBranchAddress("Harm.HCalScint.hit.row", &Harm_HCalScint_hit_row, &b_Harm_HCalScint_hit_row);
   fChain->SetBranchAddress("Harm.HCalScint.hit.col", &Harm_HCalScint_hit_col, &b_Harm_HCalScint_hit_col);
   fChain->SetBranchAddress("Harm.HCalScint.hit.cell", &Harm_HCalScint_hit_cell, &b_Harm_HCalScint_hit_cell);
   fChain->SetBranchAddress("Harm.HCalScint.hit.plane", &Harm_HCalScint_hit_plane, &b_Harm_HCalScint_hit_plane);
   fChain->SetBranchAddress("Harm.HCalScint.hit.xcell", &Harm_HCalScint_hit_xcell, &b_Harm_HCalScint_hit_xcell);
   fChain->SetBranchAddress("Harm.HCalScint.hit.ycell", &Harm_HCalScint_hit_ycell, &b_Harm_HCalScint_hit_ycell);
   fChain->SetBranchAddress("Harm.HCalScint.hit.zcell", &Harm_HCalScint_hit_zcell, &b_Harm_HCalScint_hit_zcell);
   fChain->SetBranchAddress("Harm.HCalScint.hit.xcellg", &Harm_HCalScint_hit_xcellg, &b_Harm_HCalScint_hit_xcellg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.ycellg", &Harm_HCalScint_hit_ycellg, &b_Harm_HCalScint_hit_ycellg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.zcellg", &Harm_HCalScint_hit_zcellg, &b_Harm_HCalScint_hit_zcellg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.xhit", &Harm_HCalScint_hit_xhit, &b_Harm_HCalScint_hit_xhit);
   fChain->SetBranchAddress("Harm.HCalScint.hit.yhit", &Harm_HCalScint_hit_yhit, &b_Harm_HCalScint_hit_yhit);
   fChain->SetBranchAddress("Harm.HCalScint.hit.zhit", &Harm_HCalScint_hit_zhit, &b_Harm_HCalScint_hit_zhit);
   fChain->SetBranchAddress("Harm.HCalScint.hit.sumedep", &Harm_HCalScint_hit_sumedep, &b_Harm_HCalScint_hit_sumedep);
   fChain->SetBranchAddress("Harm.HCalScint.hit.tavg", &Harm_HCalScint_hit_tavg, &b_Harm_HCalScint_hit_tavg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.trms", &Harm_HCalScint_hit_trms, &b_Harm_HCalScint_hit_trms);
   fChain->SetBranchAddress("Harm.HCalScint.hit.tmin", &Harm_HCalScint_hit_tmin, &b_Harm_HCalScint_hit_tmin);
   fChain->SetBranchAddress("Harm.HCalScint.hit.tmax", &Harm_HCalScint_hit_tmax, &b_Harm_HCalScint_hit_tmax);
   Notify();
}

Bool_t gmn_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gmn_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gmn_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gmn_tree_cxx
