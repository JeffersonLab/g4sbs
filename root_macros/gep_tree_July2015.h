//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jan 29 13:00:48 2016 by ROOT version 5.34/32
// from TTree T/Geant4 SBS Simulation
// found on file: gep_background_job1.root
//////////////////////////////////////////////////////////

#ifndef gep_tree_July2015_h
#define gep_tree_July2015_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.

class gep_tree_July2015 {
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
   Double_t        gen_thsbs;
   Double_t        gen_dbb;
   Double_t        gen_dsbs;
   Double_t        gen_dhcal;
   Double_t        gen_drich;
   Double_t        gen_dsbstrkr;
   Double_t        gen_Ebeam;
   Int_t           Earm_CDET_hit_nhits;
   vector<int>     *Earm_CDET_hit_PMT;
   vector<int>     *Earm_CDET_hit_row;
   vector<int>     *Earm_CDET_hit_col;
   vector<int>     *Earm_CDET_hit_plane;
   vector<double>  *Earm_CDET_hit_xcell;
   vector<double>  *Earm_CDET_hit_ycell;
   vector<double>  *Earm_CDET_hit_zcell;
   vector<double>  *Earm_CDET_hit_xgcell;
   vector<double>  *Earm_CDET_hit_ygcell;
   vector<double>  *Earm_CDET_hit_zgcell;
   vector<int>     *Earm_CDET_hit_NumPhotoelectrons;
   vector<double>  *Earm_CDET_hit_Time_avg;
   vector<double>  *Earm_CDET_hit_Time_rms;
   vector<double>  *Earm_CDET_hit_Time_min;
   vector<double>  *Earm_CDET_hit_Time_max;
   Int_t           Earm_CDET_Scint_hit_nhits;
   vector<int>     *Earm_CDET_Scint_hit_row;
   vector<int>     *Earm_CDET_Scint_hit_col;
   vector<int>     *Earm_CDET_Scint_hit_plane;
   vector<double>  *Earm_CDET_Scint_hit_xcell;
   vector<double>  *Earm_CDET_Scint_hit_ycell;
   vector<double>  *Earm_CDET_Scint_hit_zcell;
   vector<double>  *Earm_CDET_Scint_hit_xcellg;
   vector<double>  *Earm_CDET_Scint_hit_ycellg;
   vector<double>  *Earm_CDET_Scint_hit_zcellg;
   vector<double>  *Earm_CDET_Scint_hit_xhit;
   vector<double>  *Earm_CDET_Scint_hit_yhit;
   vector<double>  *Earm_CDET_Scint_hit_zhit;
   vector<double>  *Earm_CDET_Scint_hit_sumedep;
   vector<double>  *Earm_CDET_Scint_hit_tavg;
   vector<double>  *Earm_CDET_Scint_hit_trms;
   vector<double>  *Earm_CDET_Scint_hit_tmin;
   vector<double>  *Earm_CDET_Scint_hit_tmax;
   Int_t           Earm_CDET_Scint_part_npart;
   vector<int>     *Earm_CDET_Scint_part_PID;
   vector<int>     *Earm_CDET_Scint_part_MID;
   vector<int>     *Earm_CDET_Scint_part_TID;
   vector<int>     *Earm_CDET_Scint_part_nbounce;
   vector<int>     *Earm_CDET_Scint_part_hitindex;
   vector<double>  *Earm_CDET_Scint_part_vx;
   vector<double>  *Earm_CDET_Scint_part_vy;
   vector<double>  *Earm_CDET_Scint_part_vz;
   vector<double>  *Earm_CDET_Scint_part_px;
   vector<double>  *Earm_CDET_Scint_part_py;
   vector<double>  *Earm_CDET_Scint_part_pz;
   Int_t           Earm_ECAL_hit_nhits;
   vector<int>     *Earm_ECAL_hit_PMT;
   vector<int>     *Earm_ECAL_hit_row;
   vector<int>     *Earm_ECAL_hit_col;
   vector<int>     *Earm_ECAL_hit_plane;
   vector<double>  *Earm_ECAL_hit_xcell;
   vector<double>  *Earm_ECAL_hit_ycell;
   vector<double>  *Earm_ECAL_hit_zcell;
   vector<double>  *Earm_ECAL_hit_xgcell;
   vector<double>  *Earm_ECAL_hit_ygcell;
   vector<double>  *Earm_ECAL_hit_zgcell;
   vector<int>     *Earm_ECAL_hit_NumPhotoelectrons;
   vector<double>  *Earm_ECAL_hit_Time_avg;
   vector<double>  *Earm_ECAL_hit_Time_rms;
   vector<double>  *Earm_ECAL_hit_Time_min;
   vector<double>  *Earm_ECAL_hit_Time_max;
   Int_t           Earm_ECalTF1_hit_nhits;
   vector<int>     *Earm_ECalTF1_hit_row;
   vector<int>     *Earm_ECalTF1_hit_col;
   vector<int>     *Earm_ECalTF1_hit_plane;
   vector<double>  *Earm_ECalTF1_hit_xcell;
   vector<double>  *Earm_ECalTF1_hit_ycell;
   vector<double>  *Earm_ECalTF1_hit_zcell;
   vector<double>  *Earm_ECalTF1_hit_xcellg;
   vector<double>  *Earm_ECalTF1_hit_ycellg;
   vector<double>  *Earm_ECalTF1_hit_zcellg;
   vector<double>  *Earm_ECalTF1_hit_xhit;
   vector<double>  *Earm_ECalTF1_hit_yhit;
   vector<double>  *Earm_ECalTF1_hit_zhit;
   vector<double>  *Earm_ECalTF1_hit_sumedep;
   vector<double>  *Earm_ECalTF1_hit_tavg;
   vector<double>  *Earm_ECalTF1_hit_trms;
   vector<double>  *Earm_ECalTF1_hit_tmin;
   vector<double>  *Earm_ECalTF1_hit_tmax;
   Int_t           Earm_ECalTF1_part_npart;
   vector<int>     *Earm_ECalTF1_part_PID;
   vector<int>     *Earm_ECalTF1_part_MID;
   vector<int>     *Earm_ECalTF1_part_TID;
   vector<int>     *Earm_ECalTF1_part_nbounce;
   vector<int>     *Earm_ECalTF1_part_hitindex;
   vector<double>  *Earm_ECalTF1_part_vx;
   vector<double>  *Earm_ECalTF1_part_vy;
   vector<double>  *Earm_ECalTF1_part_vz;
   vector<double>  *Earm_ECalTF1_part_px;
   vector<double>  *Earm_ECalTF1_part_py;
   vector<double>  *Earm_ECalTF1_part_pz;
   Int_t           Harm_FPP1_hit_nhits;
   vector<int>     *Harm_FPP1_hit_plane;
   vector<int>     *Harm_FPP1_hit_strip;
   vector<double>  *Harm_FPP1_hit_x;
   vector<double>  *Harm_FPP1_hit_y;
   vector<double>  *Harm_FPP1_hit_z;
   vector<double>  *Harm_FPP1_hit_t;
   vector<double>  *Harm_FPP1_hit_trms;
   vector<double>  *Harm_FPP1_hit_tmin;
   vector<double>  *Harm_FPP1_hit_tmax;
   vector<double>  *Harm_FPP1_hit_tx;
   vector<double>  *Harm_FPP1_hit_ty;
   vector<double>  *Harm_FPP1_hit_txp;
   vector<double>  *Harm_FPP1_hit_typ;
   vector<int>     *Harm_FPP1_hit_trid;
   vector<int>     *Harm_FPP1_hit_mid;
   vector<int>     *Harm_FPP1_hit_pid;
   vector<double>  *Harm_FPP1_hit_vx;
   vector<double>  *Harm_FPP1_hit_vy;
   vector<double>  *Harm_FPP1_hit_vz;
   vector<double>  *Harm_FPP1_hit_p;
   vector<double>  *Harm_FPP1_hit_edep;
   vector<double>  *Harm_FPP1_hit_beta;
   Int_t           Harm_FPP1_Track_ntracks;
   vector<int>     *Harm_FPP1_Track_TID;
   vector<int>     *Harm_FPP1_Track_PID;
   vector<int>     *Harm_FPP1_Track_MID;
   vector<int>     *Harm_FPP1_Track_NumHits;
   vector<int>     *Harm_FPP1_Track_NumPlanes;
   vector<int>     *Harm_FPP1_Track_NDF;
   vector<double>  *Harm_FPP1_Track_Chi2fit;
   vector<double>  *Harm_FPP1_Track_Chi2true;
   vector<double>  *Harm_FPP1_Track_X;
   vector<double>  *Harm_FPP1_Track_Y;
   vector<double>  *Harm_FPP1_Track_Xp;
   vector<double>  *Harm_FPP1_Track_Yp;
   vector<double>  *Harm_FPP1_Track_T;
   vector<double>  *Harm_FPP1_Track_P;
   vector<double>  *Harm_FPP1_Track_Xfit;
   vector<double>  *Harm_FPP1_Track_Yfit;
   vector<double>  *Harm_FPP1_Track_Xpfit;
   vector<double>  *Harm_FPP1_Track_Ypfit;
   Int_t           Harm_FPP1_part_npart;
   vector<int>     *Harm_FPP1_part_PID;
   vector<int>     *Harm_FPP1_part_MID;
   vector<int>     *Harm_FPP1_part_TID;
   vector<int>     *Harm_FPP1_part_nbounce;
   vector<int>     *Harm_FPP1_part_hitindex;
   vector<double>  *Harm_FPP1_part_vx;
   vector<double>  *Harm_FPP1_part_vy;
   vector<double>  *Harm_FPP1_part_vz;
   vector<double>  *Harm_FPP1_part_px;
   vector<double>  *Harm_FPP1_part_py;
   vector<double>  *Harm_FPP1_part_pz;
   Int_t           Harm_FPP2_hit_nhits;
   vector<int>     *Harm_FPP2_hit_plane;
   vector<int>     *Harm_FPP2_hit_strip;
   vector<double>  *Harm_FPP2_hit_x;
   vector<double>  *Harm_FPP2_hit_y;
   vector<double>  *Harm_FPP2_hit_z;
   vector<double>  *Harm_FPP2_hit_t;
   vector<double>  *Harm_FPP2_hit_trms;
   vector<double>  *Harm_FPP2_hit_tmin;
   vector<double>  *Harm_FPP2_hit_tmax;
   vector<double>  *Harm_FPP2_hit_tx;
   vector<double>  *Harm_FPP2_hit_ty;
   vector<double>  *Harm_FPP2_hit_txp;
   vector<double>  *Harm_FPP2_hit_typ;
   vector<int>     *Harm_FPP2_hit_trid;
   vector<int>     *Harm_FPP2_hit_mid;
   vector<int>     *Harm_FPP2_hit_pid;
   vector<double>  *Harm_FPP2_hit_vx;
   vector<double>  *Harm_FPP2_hit_vy;
   vector<double>  *Harm_FPP2_hit_vz;
   vector<double>  *Harm_FPP2_hit_p;
   vector<double>  *Harm_FPP2_hit_edep;
   vector<double>  *Harm_FPP2_hit_beta;
   Int_t           Harm_FPP2_Track_ntracks;
   vector<int>     *Harm_FPP2_Track_TID;
   vector<int>     *Harm_FPP2_Track_PID;
   vector<int>     *Harm_FPP2_Track_MID;
   vector<int>     *Harm_FPP2_Track_NumHits;
   vector<int>     *Harm_FPP2_Track_NumPlanes;
   vector<int>     *Harm_FPP2_Track_NDF;
   vector<double>  *Harm_FPP2_Track_Chi2fit;
   vector<double>  *Harm_FPP2_Track_Chi2true;
   vector<double>  *Harm_FPP2_Track_X;
   vector<double>  *Harm_FPP2_Track_Y;
   vector<double>  *Harm_FPP2_Track_Xp;
   vector<double>  *Harm_FPP2_Track_Yp;
   vector<double>  *Harm_FPP2_Track_T;
   vector<double>  *Harm_FPP2_Track_P;
   vector<double>  *Harm_FPP2_Track_Xfit;
   vector<double>  *Harm_FPP2_Track_Yfit;
   vector<double>  *Harm_FPP2_Track_Xpfit;
   vector<double>  *Harm_FPP2_Track_Ypfit;
   Int_t           Harm_FPP2_part_npart;
   vector<int>     *Harm_FPP2_part_PID;
   vector<int>     *Harm_FPP2_part_MID;
   vector<int>     *Harm_FPP2_part_TID;
   vector<int>     *Harm_FPP2_part_nbounce;
   vector<int>     *Harm_FPP2_part_hitindex;
   vector<double>  *Harm_FPP2_part_vx;
   vector<double>  *Harm_FPP2_part_vy;
   vector<double>  *Harm_FPP2_part_vz;
   vector<double>  *Harm_FPP2_part_px;
   vector<double>  *Harm_FPP2_part_py;
   vector<double>  *Harm_FPP2_part_pz;
   Int_t           Harm_FT_hit_nhits;
   vector<int>     *Harm_FT_hit_plane;
   vector<int>     *Harm_FT_hit_strip;
   vector<double>  *Harm_FT_hit_x;
   vector<double>  *Harm_FT_hit_y;
   vector<double>  *Harm_FT_hit_z;
   vector<double>  *Harm_FT_hit_t;
   vector<double>  *Harm_FT_hit_trms;
   vector<double>  *Harm_FT_hit_tmin;
   vector<double>  *Harm_FT_hit_tmax;
   vector<double>  *Harm_FT_hit_tx;
   vector<double>  *Harm_FT_hit_ty;
   vector<double>  *Harm_FT_hit_txp;
   vector<double>  *Harm_FT_hit_typ;
   vector<int>     *Harm_FT_hit_trid;
   vector<int>     *Harm_FT_hit_mid;
   vector<int>     *Harm_FT_hit_pid;
   vector<double>  *Harm_FT_hit_vx;
   vector<double>  *Harm_FT_hit_vy;
   vector<double>  *Harm_FT_hit_vz;
   vector<double>  *Harm_FT_hit_p;
   vector<double>  *Harm_FT_hit_edep;
   vector<double>  *Harm_FT_hit_beta;
   Int_t           Harm_FT_Track_ntracks;
   vector<int>     *Harm_FT_Track_TID;
   vector<int>     *Harm_FT_Track_PID;
   vector<int>     *Harm_FT_Track_MID;
   vector<int>     *Harm_FT_Track_NumHits;
   vector<int>     *Harm_FT_Track_NumPlanes;
   vector<int>     *Harm_FT_Track_NDF;
   vector<double>  *Harm_FT_Track_Chi2fit;
   vector<double>  *Harm_FT_Track_Chi2true;
   vector<double>  *Harm_FT_Track_X;
   vector<double>  *Harm_FT_Track_Y;
   vector<double>  *Harm_FT_Track_Xp;
   vector<double>  *Harm_FT_Track_Yp;
   vector<double>  *Harm_FT_Track_T;
   vector<double>  *Harm_FT_Track_P;
   vector<double>  *Harm_FT_Track_Xfit;
   vector<double>  *Harm_FT_Track_Yfit;
   vector<double>  *Harm_FT_Track_Xpfit;
   vector<double>  *Harm_FT_Track_Ypfit;
   Int_t           Harm_FT_part_npart;
   vector<int>     *Harm_FT_part_PID;
   vector<int>     *Harm_FT_part_MID;
   vector<int>     *Harm_FT_part_TID;
   vector<int>     *Harm_FT_part_nbounce;
   vector<int>     *Harm_FT_part_hitindex;
   vector<double>  *Harm_FT_part_vx;
   vector<double>  *Harm_FT_part_vy;
   vector<double>  *Harm_FT_part_vz;
   vector<double>  *Harm_FT_part_px;
   vector<double>  *Harm_FT_part_py;
   vector<double>  *Harm_FT_part_pz;
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
   Int_t           Harm_HCalScint_part_npart;
   vector<int>     *Harm_HCalScint_part_PID;
   vector<int>     *Harm_HCalScint_part_MID;
   vector<int>     *Harm_HCalScint_part_TID;
   vector<int>     *Harm_HCalScint_part_nbounce;
   vector<int>     *Harm_HCalScint_part_hitindex;
   vector<double>  *Harm_HCalScint_part_vx;
   vector<double>  *Harm_HCalScint_part_vy;
   vector<double>  *Harm_HCalScint_part_vz;
   vector<double>  *Harm_HCalScint_part_px;
   vector<double>  *Harm_HCalScint_part_py;
   vector<double>  *Harm_HCalScint_part_pz;

   // List of branches
   TBranch        *b_ev;   //!
   TBranch        *b_gen;   //!
   TBranch        *b_Earm_CDET_hit_nhits;   //!
   TBranch        *b_Earm_CDET_hit_PMT;   //!
   TBranch        *b_Earm_CDET_hit_row;   //!
   TBranch        *b_Earm_CDET_hit_col;   //!
   TBranch        *b_Earm_CDET_hit_plane;   //!
   TBranch        *b_Earm_CDET_hit_xcell;   //!
   TBranch        *b_Earm_CDET_hit_ycell;   //!
   TBranch        *b_Earm_CDET_hit_zcell;   //!
   TBranch        *b_Earm_CDET_hit_xgcell;   //!
   TBranch        *b_Earm_CDET_hit_ygcell;   //!
   TBranch        *b_Earm_CDET_hit_zgcell;   //!
   TBranch        *b_Earm_CDET_hit_NumPhotoelectrons;   //!
   TBranch        *b_Earm_CDET_hit_Time_avg;   //!
   TBranch        *b_Earm_CDET_hit_Time_rms;   //!
   TBranch        *b_Earm_CDET_hit_Time_min;   //!
   TBranch        *b_Earm_CDET_hit_Time_max;   //!
   TBranch        *b_Earm_CDET_Scint_hit_nhits;   //!
   TBranch        *b_Earm_CDET_Scint_hit_row;   //!
   TBranch        *b_Earm_CDET_Scint_hit_col;   //!
   TBranch        *b_Earm_CDET_Scint_hit_plane;   //!
   TBranch        *b_Earm_CDET_Scint_hit_xcell;   //!
   TBranch        *b_Earm_CDET_Scint_hit_ycell;   //!
   TBranch        *b_Earm_CDET_Scint_hit_zcell;   //!
   TBranch        *b_Earm_CDET_Scint_hit_xcellg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_ycellg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_zcellg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_xhit;   //!
   TBranch        *b_Earm_CDET_Scint_hit_yhit;   //!
   TBranch        *b_Earm_CDET_Scint_hit_zhit;   //!
   TBranch        *b_Earm_CDET_Scint_hit_sumedep;   //!
   TBranch        *b_Earm_CDET_Scint_hit_tavg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_trms;   //!
   TBranch        *b_Earm_CDET_Scint_hit_tmin;   //!
   TBranch        *b_Earm_CDET_Scint_hit_tmax;   //!
   TBranch        *b_Earm_CDET_Scint_part_npart;   //!
   TBranch        *b_Earm_CDET_Scint_part_PID;   //!
   TBranch        *b_Earm_CDET_Scint_part_MID;   //!
   TBranch        *b_Earm_CDET_Scint_part_TID;   //!
   TBranch        *b_Earm_CDET_Scint_part_nbounce;   //!
   TBranch        *b_Earm_CDET_Scint_part_hitindex;   //!
   TBranch        *b_Earm_CDET_Scint_part_vx;   //!
   TBranch        *b_Earm_CDET_Scint_part_vy;   //!
   TBranch        *b_Earm_CDET_Scint_part_vz;   //!
   TBranch        *b_Earm_CDET_Scint_part_px;   //!
   TBranch        *b_Earm_CDET_Scint_part_py;   //!
   TBranch        *b_Earm_CDET_Scint_part_pz;   //!
   TBranch        *b_Earm_ECAL_hit_nhits;   //!
   TBranch        *b_Earm_ECAL_hit_PMT;   //!
   TBranch        *b_Earm_ECAL_hit_row;   //!
   TBranch        *b_Earm_ECAL_hit_col;   //!
   TBranch        *b_Earm_ECAL_hit_plane;   //!
   TBranch        *b_Earm_ECAL_hit_xcell;   //!
   TBranch        *b_Earm_ECAL_hit_ycell;   //!
   TBranch        *b_Earm_ECAL_hit_zcell;   //!
   TBranch        *b_Earm_ECAL_hit_xgcell;   //!
   TBranch        *b_Earm_ECAL_hit_ygcell;   //!
   TBranch        *b_Earm_ECAL_hit_zgcell;   //!
   TBranch        *b_Earm_ECAL_hit_NumPhotoelectrons;   //!
   TBranch        *b_Earm_ECAL_hit_Time_avg;   //!
   TBranch        *b_Earm_ECAL_hit_Time_rms;   //!
   TBranch        *b_Earm_ECAL_hit_Time_min;   //!
   TBranch        *b_Earm_ECAL_hit_Time_max;   //!
   TBranch        *b_Earm_ECalTF1_hit_nhits;   //!
   TBranch        *b_Earm_ECalTF1_hit_row;   //!
   TBranch        *b_Earm_ECalTF1_hit_col;   //!
   TBranch        *b_Earm_ECalTF1_hit_plane;   //!
   TBranch        *b_Earm_ECalTF1_hit_xcell;   //!
   TBranch        *b_Earm_ECalTF1_hit_ycell;   //!
   TBranch        *b_Earm_ECalTF1_hit_zcell;   //!
   TBranch        *b_Earm_ECalTF1_hit_xcellg;   //!
   TBranch        *b_Earm_ECalTF1_hit_ycellg;   //!
   TBranch        *b_Earm_ECalTF1_hit_zcellg;   //!
   TBranch        *b_Earm_ECalTF1_hit_xhit;   //!
   TBranch        *b_Earm_ECalTF1_hit_yhit;   //!
   TBranch        *b_Earm_ECalTF1_hit_zhit;   //!
   TBranch        *b_Earm_ECalTF1_hit_sumedep;   //!
   TBranch        *b_Earm_ECalTF1_hit_tavg;   //!
   TBranch        *b_Earm_ECalTF1_hit_trms;   //!
   TBranch        *b_Earm_ECalTF1_hit_tmin;   //!
   TBranch        *b_Earm_ECalTF1_hit_tmax;   //!
   TBranch        *b_Earm_ECalTF1_part_npart;   //!
   TBranch        *b_Earm_ECalTF1_part_PID;   //!
   TBranch        *b_Earm_ECalTF1_part_MID;   //!
   TBranch        *b_Earm_ECalTF1_part_TID;   //!
   TBranch        *b_Earm_ECalTF1_part_nbounce;   //!
   TBranch        *b_Earm_ECalTF1_part_hitindex;   //!
   TBranch        *b_Earm_ECalTF1_part_vx;   //!
   TBranch        *b_Earm_ECalTF1_part_vy;   //!
   TBranch        *b_Earm_ECalTF1_part_vz;   //!
   TBranch        *b_Earm_ECalTF1_part_px;   //!
   TBranch        *b_Earm_ECalTF1_part_py;   //!
   TBranch        *b_Earm_ECalTF1_part_pz;   //!
   TBranch        *b_Harm_FPP1_hit_nhits;   //!
   TBranch        *b_Harm_FPP1_hit_plane;   //!
   TBranch        *b_Harm_FPP1_hit_strip;   //!
   TBranch        *b_Harm_FPP1_hit_x;   //!
   TBranch        *b_Harm_FPP1_hit_y;   //!
   TBranch        *b_Harm_FPP1_hit_z;   //!
   TBranch        *b_Harm_FPP1_hit_t;   //!
   TBranch        *b_Harm_FPP1_hit_trms;   //!
   TBranch        *b_Harm_FPP1_hit_tmin;   //!
   TBranch        *b_Harm_FPP1_hit_tmax;   //!
   TBranch        *b_Harm_FPP1_hit_tx;   //!
   TBranch        *b_Harm_FPP1_hit_ty;   //!
   TBranch        *b_Harm_FPP1_hit_txp;   //!
   TBranch        *b_Harm_FPP1_hit_typ;   //!
   TBranch        *b_Harm_FPP1_hit_trid;   //!
   TBranch        *b_Harm_FPP1_hit_mid;   //!
   TBranch        *b_Harm_FPP1_hit_pid;   //!
   TBranch        *b_Harm_FPP1_hit_vx;   //!
   TBranch        *b_Harm_FPP1_hit_vy;   //!
   TBranch        *b_Harm_FPP1_hit_vz;   //!
   TBranch        *b_Harm_FPP1_hit_p;   //!
   TBranch        *b_Harm_FPP1_hit_edep;   //!
   TBranch        *b_Harm_FPP1_hit_beta;   //!
   TBranch        *b_Harm_FPP1_Track_ntracks;   //!
   TBranch        *b_Harm_FPP1_Track_TID;   //!
   TBranch        *b_Harm_FPP1_Track_PID;   //!
   TBranch        *b_Harm_FPP1_Track_MID;   //!
   TBranch        *b_Harm_FPP1_Track_NumHits;   //!
   TBranch        *b_Harm_FPP1_Track_NumPlanes;   //!
   TBranch        *b_Harm_FPP1_Track_NDF;   //!
   TBranch        *b_Harm_FPP1_Track_Chi2fit;   //!
   TBranch        *b_Harm_FPP1_Track_Chi2true;   //!
   TBranch        *b_Harm_FPP1_Track_X;   //!
   TBranch        *b_Harm_FPP1_Track_Y;   //!
   TBranch        *b_Harm_FPP1_Track_Xp;   //!
   TBranch        *b_Harm_FPP1_Track_Yp;   //!
   TBranch        *b_Harm_FPP1_Track_T;   //!
   TBranch        *b_Harm_FPP1_Track_P;   //!
   TBranch        *b_Harm_FPP1_Track_Xfit;   //!
   TBranch        *b_Harm_FPP1_Track_Yfit;   //!
   TBranch        *b_Harm_FPP1_Track_Xpfit;   //!
   TBranch        *b_Harm_FPP1_Track_Ypfit;   //!
   TBranch        *b_Harm_FPP1_part_npart;   //!
   TBranch        *b_Harm_FPP1_part_PID;   //!
   TBranch        *b_Harm_FPP1_part_MID;   //!
   TBranch        *b_Harm_FPP1_part_TID;   //!
   TBranch        *b_Harm_FPP1_part_nbounce;   //!
   TBranch        *b_Harm_FPP1_part_hitindex;   //!
   TBranch        *b_Harm_FPP1_part_vx;   //!
   TBranch        *b_Harm_FPP1_part_vy;   //!
   TBranch        *b_Harm_FPP1_part_vz;   //!
   TBranch        *b_Harm_FPP1_part_px;   //!
   TBranch        *b_Harm_FPP1_part_py;   //!
   TBranch        *b_Harm_FPP1_part_pz;   //!
   TBranch        *b_Harm_FPP2_hit_nhits;   //!
   TBranch        *b_Harm_FPP2_hit_plane;   //!
   TBranch        *b_Harm_FPP2_hit_strip;   //!
   TBranch        *b_Harm_FPP2_hit_x;   //!
   TBranch        *b_Harm_FPP2_hit_y;   //!
   TBranch        *b_Harm_FPP2_hit_z;   //!
   TBranch        *b_Harm_FPP2_hit_t;   //!
   TBranch        *b_Harm_FPP2_hit_trms;   //!
   TBranch        *b_Harm_FPP2_hit_tmin;   //!
   TBranch        *b_Harm_FPP2_hit_tmax;   //!
   TBranch        *b_Harm_FPP2_hit_tx;   //!
   TBranch        *b_Harm_FPP2_hit_ty;   //!
   TBranch        *b_Harm_FPP2_hit_txp;   //!
   TBranch        *b_Harm_FPP2_hit_typ;   //!
   TBranch        *b_Harm_FPP2_hit_trid;   //!
   TBranch        *b_Harm_FPP2_hit_mid;   //!
   TBranch        *b_Harm_FPP2_hit_pid;   //!
   TBranch        *b_Harm_FPP2_hit_vx;   //!
   TBranch        *b_Harm_FPP2_hit_vy;   //!
   TBranch        *b_Harm_FPP2_hit_vz;   //!
   TBranch        *b_Harm_FPP2_hit_p;   //!
   TBranch        *b_Harm_FPP2_hit_edep;   //!
   TBranch        *b_Harm_FPP2_hit_beta;   //!
   TBranch        *b_Harm_FPP2_Track_ntracks;   //!
   TBranch        *b_Harm_FPP2_Track_TID;   //!
   TBranch        *b_Harm_FPP2_Track_PID;   //!
   TBranch        *b_Harm_FPP2_Track_MID;   //!
   TBranch        *b_Harm_FPP2_Track_NumHits;   //!
   TBranch        *b_Harm_FPP2_Track_NumPlanes;   //!
   TBranch        *b_Harm_FPP2_Track_NDF;   //!
   TBranch        *b_Harm_FPP2_Track_Chi2fit;   //!
   TBranch        *b_Harm_FPP2_Track_Chi2true;   //!
   TBranch        *b_Harm_FPP2_Track_X;   //!
   TBranch        *b_Harm_FPP2_Track_Y;   //!
   TBranch        *b_Harm_FPP2_Track_Xp;   //!
   TBranch        *b_Harm_FPP2_Track_Yp;   //!
   TBranch        *b_Harm_FPP2_Track_T;   //!
   TBranch        *b_Harm_FPP2_Track_P;   //!
   TBranch        *b_Harm_FPP2_Track_Xfit;   //!
   TBranch        *b_Harm_FPP2_Track_Yfit;   //!
   TBranch        *b_Harm_FPP2_Track_Xpfit;   //!
   TBranch        *b_Harm_FPP2_Track_Ypfit;   //!
   TBranch        *b_Harm_FPP2_part_npart;   //!
   TBranch        *b_Harm_FPP2_part_PID;   //!
   TBranch        *b_Harm_FPP2_part_MID;   //!
   TBranch        *b_Harm_FPP2_part_TID;   //!
   TBranch        *b_Harm_FPP2_part_nbounce;   //!
   TBranch        *b_Harm_FPP2_part_hitindex;   //!
   TBranch        *b_Harm_FPP2_part_vx;   //!
   TBranch        *b_Harm_FPP2_part_vy;   //!
   TBranch        *b_Harm_FPP2_part_vz;   //!
   TBranch        *b_Harm_FPP2_part_px;   //!
   TBranch        *b_Harm_FPP2_part_py;   //!
   TBranch        *b_Harm_FPP2_part_pz;   //!
   TBranch        *b_Harm_FT_hit_nhits;   //!
   TBranch        *b_Harm_FT_hit_plane;   //!
   TBranch        *b_Harm_FT_hit_strip;   //!
   TBranch        *b_Harm_FT_hit_x;   //!
   TBranch        *b_Harm_FT_hit_y;   //!
   TBranch        *b_Harm_FT_hit_z;   //!
   TBranch        *b_Harm_FT_hit_t;   //!
   TBranch        *b_Harm_FT_hit_trms;   //!
   TBranch        *b_Harm_FT_hit_tmin;   //!
   TBranch        *b_Harm_FT_hit_tmax;   //!
   TBranch        *b_Harm_FT_hit_tx;   //!
   TBranch        *b_Harm_FT_hit_ty;   //!
   TBranch        *b_Harm_FT_hit_txp;   //!
   TBranch        *b_Harm_FT_hit_typ;   //!
   TBranch        *b_Harm_FT_hit_trid;   //!
   TBranch        *b_Harm_FT_hit_mid;   //!
   TBranch        *b_Harm_FT_hit_pid;   //!
   TBranch        *b_Harm_FT_hit_vx;   //!
   TBranch        *b_Harm_FT_hit_vy;   //!
   TBranch        *b_Harm_FT_hit_vz;   //!
   TBranch        *b_Harm_FT_hit_p;   //!
   TBranch        *b_Harm_FT_hit_edep;   //!
   TBranch        *b_Harm_FT_hit_beta;   //!
   TBranch        *b_Harm_FT_Track_ntracks;   //!
   TBranch        *b_Harm_FT_Track_TID;   //!
   TBranch        *b_Harm_FT_Track_PID;   //!
   TBranch        *b_Harm_FT_Track_MID;   //!
   TBranch        *b_Harm_FT_Track_NumHits;   //!
   TBranch        *b_Harm_FT_Track_NumPlanes;   //!
   TBranch        *b_Harm_FT_Track_NDF;   //!
   TBranch        *b_Harm_FT_Track_Chi2fit;   //!
   TBranch        *b_Harm_FT_Track_Chi2true;   //!
   TBranch        *b_Harm_FT_Track_X;   //!
   TBranch        *b_Harm_FT_Track_Y;   //!
   TBranch        *b_Harm_FT_Track_Xp;   //!
   TBranch        *b_Harm_FT_Track_Yp;   //!
   TBranch        *b_Harm_FT_Track_T;   //!
   TBranch        *b_Harm_FT_Track_P;   //!
   TBranch        *b_Harm_FT_Track_Xfit;   //!
   TBranch        *b_Harm_FT_Track_Yfit;   //!
   TBranch        *b_Harm_FT_Track_Xpfit;   //!
   TBranch        *b_Harm_FT_Track_Ypfit;   //!
   TBranch        *b_Harm_FT_part_npart;   //!
   TBranch        *b_Harm_FT_part_PID;   //!
   TBranch        *b_Harm_FT_part_MID;   //!
   TBranch        *b_Harm_FT_part_TID;   //!
   TBranch        *b_Harm_FT_part_nbounce;   //!
   TBranch        *b_Harm_FT_part_hitindex;   //!
   TBranch        *b_Harm_FT_part_vx;   //!
   TBranch        *b_Harm_FT_part_vy;   //!
   TBranch        *b_Harm_FT_part_vz;   //!
   TBranch        *b_Harm_FT_part_px;   //!
   TBranch        *b_Harm_FT_part_py;   //!
   TBranch        *b_Harm_FT_part_pz;   //!
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
   TBranch        *b_Harm_HCalScint_part_npart;   //!
   TBranch        *b_Harm_HCalScint_part_PID;   //!
   TBranch        *b_Harm_HCalScint_part_MID;   //!
   TBranch        *b_Harm_HCalScint_part_TID;   //!
   TBranch        *b_Harm_HCalScint_part_nbounce;   //!
   TBranch        *b_Harm_HCalScint_part_hitindex;   //!
   TBranch        *b_Harm_HCalScint_part_vx;   //!
   TBranch        *b_Harm_HCalScint_part_vy;   //!
   TBranch        *b_Harm_HCalScint_part_vz;   //!
   TBranch        *b_Harm_HCalScint_part_px;   //!
   TBranch        *b_Harm_HCalScint_part_py;   //!
   TBranch        *b_Harm_HCalScint_part_pz;   //!

   gep_tree_July2015(TTree *tree=0);
   virtual ~gep_tree_July2015();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gep_tree_July2015_cxx
gep_tree_July2015::gep_tree_July2015(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("gep_background_job1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("gep_background_job1.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

gep_tree_July2015::~gep_tree_July2015()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gep_tree_July2015::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gep_tree_July2015::LoadTree(Long64_t entry)
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

void gep_tree_July2015::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Earm_CDET_hit_PMT = 0;
   Earm_CDET_hit_row = 0;
   Earm_CDET_hit_col = 0;
   Earm_CDET_hit_plane = 0;
   Earm_CDET_hit_xcell = 0;
   Earm_CDET_hit_ycell = 0;
   Earm_CDET_hit_zcell = 0;
   Earm_CDET_hit_xgcell = 0;
   Earm_CDET_hit_ygcell = 0;
   Earm_CDET_hit_zgcell = 0;
   Earm_CDET_hit_NumPhotoelectrons = 0;
   Earm_CDET_hit_Time_avg = 0;
   Earm_CDET_hit_Time_rms = 0;
   Earm_CDET_hit_Time_min = 0;
   Earm_CDET_hit_Time_max = 0;
   Earm_CDET_Scint_hit_row = 0;
   Earm_CDET_Scint_hit_col = 0;
   Earm_CDET_Scint_hit_plane = 0;
   Earm_CDET_Scint_hit_xcell = 0;
   Earm_CDET_Scint_hit_ycell = 0;
   Earm_CDET_Scint_hit_zcell = 0;
   Earm_CDET_Scint_hit_xcellg = 0;
   Earm_CDET_Scint_hit_ycellg = 0;
   Earm_CDET_Scint_hit_zcellg = 0;
   Earm_CDET_Scint_hit_xhit = 0;
   Earm_CDET_Scint_hit_yhit = 0;
   Earm_CDET_Scint_hit_zhit = 0;
   Earm_CDET_Scint_hit_sumedep = 0;
   Earm_CDET_Scint_hit_tavg = 0;
   Earm_CDET_Scint_hit_trms = 0;
   Earm_CDET_Scint_hit_tmin = 0;
   Earm_CDET_Scint_hit_tmax = 0;
   Earm_CDET_Scint_part_PID = 0;
   Earm_CDET_Scint_part_MID = 0;
   Earm_CDET_Scint_part_TID = 0;
   Earm_CDET_Scint_part_nbounce = 0;
   Earm_CDET_Scint_part_hitindex = 0;
   Earm_CDET_Scint_part_vx = 0;
   Earm_CDET_Scint_part_vy = 0;
   Earm_CDET_Scint_part_vz = 0;
   Earm_CDET_Scint_part_px = 0;
   Earm_CDET_Scint_part_py = 0;
   Earm_CDET_Scint_part_pz = 0;
   Earm_ECAL_hit_PMT = 0;
   Earm_ECAL_hit_row = 0;
   Earm_ECAL_hit_col = 0;
   Earm_ECAL_hit_plane = 0;
   Earm_ECAL_hit_xcell = 0;
   Earm_ECAL_hit_ycell = 0;
   Earm_ECAL_hit_zcell = 0;
   Earm_ECAL_hit_xgcell = 0;
   Earm_ECAL_hit_ygcell = 0;
   Earm_ECAL_hit_zgcell = 0;
   Earm_ECAL_hit_NumPhotoelectrons = 0;
   Earm_ECAL_hit_Time_avg = 0;
   Earm_ECAL_hit_Time_rms = 0;
   Earm_ECAL_hit_Time_min = 0;
   Earm_ECAL_hit_Time_max = 0;
   Earm_ECalTF1_hit_row = 0;
   Earm_ECalTF1_hit_col = 0;
   Earm_ECalTF1_hit_plane = 0;
   Earm_ECalTF1_hit_xcell = 0;
   Earm_ECalTF1_hit_ycell = 0;
   Earm_ECalTF1_hit_zcell = 0;
   Earm_ECalTF1_hit_xcellg = 0;
   Earm_ECalTF1_hit_ycellg = 0;
   Earm_ECalTF1_hit_zcellg = 0;
   Earm_ECalTF1_hit_xhit = 0;
   Earm_ECalTF1_hit_yhit = 0;
   Earm_ECalTF1_hit_zhit = 0;
   Earm_ECalTF1_hit_sumedep = 0;
   Earm_ECalTF1_hit_tavg = 0;
   Earm_ECalTF1_hit_trms = 0;
   Earm_ECalTF1_hit_tmin = 0;
   Earm_ECalTF1_hit_tmax = 0;
   Earm_ECalTF1_part_PID = 0;
   Earm_ECalTF1_part_MID = 0;
   Earm_ECalTF1_part_TID = 0;
   Earm_ECalTF1_part_nbounce = 0;
   Earm_ECalTF1_part_hitindex = 0;
   Earm_ECalTF1_part_vx = 0;
   Earm_ECalTF1_part_vy = 0;
   Earm_ECalTF1_part_vz = 0;
   Earm_ECalTF1_part_px = 0;
   Earm_ECalTF1_part_py = 0;
   Earm_ECalTF1_part_pz = 0;
   Harm_FPP1_hit_plane = 0;
   Harm_FPP1_hit_strip = 0;
   Harm_FPP1_hit_x = 0;
   Harm_FPP1_hit_y = 0;
   Harm_FPP1_hit_z = 0;
   Harm_FPP1_hit_t = 0;
   Harm_FPP1_hit_trms = 0;
   Harm_FPP1_hit_tmin = 0;
   Harm_FPP1_hit_tmax = 0;
   Harm_FPP1_hit_tx = 0;
   Harm_FPP1_hit_ty = 0;
   Harm_FPP1_hit_txp = 0;
   Harm_FPP1_hit_typ = 0;
   Harm_FPP1_hit_trid = 0;
   Harm_FPP1_hit_mid = 0;
   Harm_FPP1_hit_pid = 0;
   Harm_FPP1_hit_vx = 0;
   Harm_FPP1_hit_vy = 0;
   Harm_FPP1_hit_vz = 0;
   Harm_FPP1_hit_p = 0;
   Harm_FPP1_hit_edep = 0;
   Harm_FPP1_hit_beta = 0;
   Harm_FPP1_Track_TID = 0;
   Harm_FPP1_Track_PID = 0;
   Harm_FPP1_Track_MID = 0;
   Harm_FPP1_Track_NumHits = 0;
   Harm_FPP1_Track_NumPlanes = 0;
   Harm_FPP1_Track_NDF = 0;
   Harm_FPP1_Track_Chi2fit = 0;
   Harm_FPP1_Track_Chi2true = 0;
   Harm_FPP1_Track_X = 0;
   Harm_FPP1_Track_Y = 0;
   Harm_FPP1_Track_Xp = 0;
   Harm_FPP1_Track_Yp = 0;
   Harm_FPP1_Track_T = 0;
   Harm_FPP1_Track_P = 0;
   Harm_FPP1_Track_Xfit = 0;
   Harm_FPP1_Track_Yfit = 0;
   Harm_FPP1_Track_Xpfit = 0;
   Harm_FPP1_Track_Ypfit = 0;
   Harm_FPP1_part_PID = 0;
   Harm_FPP1_part_MID = 0;
   Harm_FPP1_part_TID = 0;
   Harm_FPP1_part_nbounce = 0;
   Harm_FPP1_part_hitindex = 0;
   Harm_FPP1_part_vx = 0;
   Harm_FPP1_part_vy = 0;
   Harm_FPP1_part_vz = 0;
   Harm_FPP1_part_px = 0;
   Harm_FPP1_part_py = 0;
   Harm_FPP1_part_pz = 0;
   Harm_FPP2_hit_plane = 0;
   Harm_FPP2_hit_strip = 0;
   Harm_FPP2_hit_x = 0;
   Harm_FPP2_hit_y = 0;
   Harm_FPP2_hit_z = 0;
   Harm_FPP2_hit_t = 0;
   Harm_FPP2_hit_trms = 0;
   Harm_FPP2_hit_tmin = 0;
   Harm_FPP2_hit_tmax = 0;
   Harm_FPP2_hit_tx = 0;
   Harm_FPP2_hit_ty = 0;
   Harm_FPP2_hit_txp = 0;
   Harm_FPP2_hit_typ = 0;
   Harm_FPP2_hit_trid = 0;
   Harm_FPP2_hit_mid = 0;
   Harm_FPP2_hit_pid = 0;
   Harm_FPP2_hit_vx = 0;
   Harm_FPP2_hit_vy = 0;
   Harm_FPP2_hit_vz = 0;
   Harm_FPP2_hit_p = 0;
   Harm_FPP2_hit_edep = 0;
   Harm_FPP2_hit_beta = 0;
   Harm_FPP2_Track_TID = 0;
   Harm_FPP2_Track_PID = 0;
   Harm_FPP2_Track_MID = 0;
   Harm_FPP2_Track_NumHits = 0;
   Harm_FPP2_Track_NumPlanes = 0;
   Harm_FPP2_Track_NDF = 0;
   Harm_FPP2_Track_Chi2fit = 0;
   Harm_FPP2_Track_Chi2true = 0;
   Harm_FPP2_Track_X = 0;
   Harm_FPP2_Track_Y = 0;
   Harm_FPP2_Track_Xp = 0;
   Harm_FPP2_Track_Yp = 0;
   Harm_FPP2_Track_T = 0;
   Harm_FPP2_Track_P = 0;
   Harm_FPP2_Track_Xfit = 0;
   Harm_FPP2_Track_Yfit = 0;
   Harm_FPP2_Track_Xpfit = 0;
   Harm_FPP2_Track_Ypfit = 0;
   Harm_FPP2_part_PID = 0;
   Harm_FPP2_part_MID = 0;
   Harm_FPP2_part_TID = 0;
   Harm_FPP2_part_nbounce = 0;
   Harm_FPP2_part_hitindex = 0;
   Harm_FPP2_part_vx = 0;
   Harm_FPP2_part_vy = 0;
   Harm_FPP2_part_vz = 0;
   Harm_FPP2_part_px = 0;
   Harm_FPP2_part_py = 0;
   Harm_FPP2_part_pz = 0;
   Harm_FT_hit_plane = 0;
   Harm_FT_hit_strip = 0;
   Harm_FT_hit_x = 0;
   Harm_FT_hit_y = 0;
   Harm_FT_hit_z = 0;
   Harm_FT_hit_t = 0;
   Harm_FT_hit_trms = 0;
   Harm_FT_hit_tmin = 0;
   Harm_FT_hit_tmax = 0;
   Harm_FT_hit_tx = 0;
   Harm_FT_hit_ty = 0;
   Harm_FT_hit_txp = 0;
   Harm_FT_hit_typ = 0;
   Harm_FT_hit_trid = 0;
   Harm_FT_hit_mid = 0;
   Harm_FT_hit_pid = 0;
   Harm_FT_hit_vx = 0;
   Harm_FT_hit_vy = 0;
   Harm_FT_hit_vz = 0;
   Harm_FT_hit_p = 0;
   Harm_FT_hit_edep = 0;
   Harm_FT_hit_beta = 0;
   Harm_FT_Track_TID = 0;
   Harm_FT_Track_PID = 0;
   Harm_FT_Track_MID = 0;
   Harm_FT_Track_NumHits = 0;
   Harm_FT_Track_NumPlanes = 0;
   Harm_FT_Track_NDF = 0;
   Harm_FT_Track_Chi2fit = 0;
   Harm_FT_Track_Chi2true = 0;
   Harm_FT_Track_X = 0;
   Harm_FT_Track_Y = 0;
   Harm_FT_Track_Xp = 0;
   Harm_FT_Track_Yp = 0;
   Harm_FT_Track_T = 0;
   Harm_FT_Track_P = 0;
   Harm_FT_Track_Xfit = 0;
   Harm_FT_Track_Yfit = 0;
   Harm_FT_Track_Xpfit = 0;
   Harm_FT_Track_Ypfit = 0;
   Harm_FT_part_PID = 0;
   Harm_FT_part_MID = 0;
   Harm_FT_part_TID = 0;
   Harm_FT_part_nbounce = 0;
   Harm_FT_part_hitindex = 0;
   Harm_FT_part_vx = 0;
   Harm_FT_part_vy = 0;
   Harm_FT_part_vz = 0;
   Harm_FT_part_px = 0;
   Harm_FT_part_py = 0;
   Harm_FT_part_pz = 0;
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
   Harm_HCalScint_part_PID = 0;
   Harm_HCalScint_part_MID = 0;
   Harm_HCalScint_part_TID = 0;
   Harm_HCalScint_part_nbounce = 0;
   Harm_HCalScint_part_hitindex = 0;
   Harm_HCalScint_part_vx = 0;
   Harm_HCalScint_part_vy = 0;
   Harm_HCalScint_part_vz = 0;
   Harm_HCalScint_part_px = 0;
   Harm_HCalScint_part_py = 0;
   Harm_HCalScint_part_pz = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev", &ev_count, &b_ev);
   fChain->SetBranchAddress("gen", &gen_thbb, &b_gen);
   fChain->SetBranchAddress("Earm.CDET.hit.nhits", &Earm_CDET_hit_nhits, &b_Earm_CDET_hit_nhits);
   fChain->SetBranchAddress("Earm.CDET.hit.PMT", &Earm_CDET_hit_PMT, &b_Earm_CDET_hit_PMT);
   fChain->SetBranchAddress("Earm.CDET.hit.row", &Earm_CDET_hit_row, &b_Earm_CDET_hit_row);
   fChain->SetBranchAddress("Earm.CDET.hit.col", &Earm_CDET_hit_col, &b_Earm_CDET_hit_col);
   fChain->SetBranchAddress("Earm.CDET.hit.plane", &Earm_CDET_hit_plane, &b_Earm_CDET_hit_plane);
   fChain->SetBranchAddress("Earm.CDET.hit.xcell", &Earm_CDET_hit_xcell, &b_Earm_CDET_hit_xcell);
   fChain->SetBranchAddress("Earm.CDET.hit.ycell", &Earm_CDET_hit_ycell, &b_Earm_CDET_hit_ycell);
   fChain->SetBranchAddress("Earm.CDET.hit.zcell", &Earm_CDET_hit_zcell, &b_Earm_CDET_hit_zcell);
   fChain->SetBranchAddress("Earm.CDET.hit.xgcell", &Earm_CDET_hit_xgcell, &b_Earm_CDET_hit_xgcell);
   fChain->SetBranchAddress("Earm.CDET.hit.ygcell", &Earm_CDET_hit_ygcell, &b_Earm_CDET_hit_ygcell);
   fChain->SetBranchAddress("Earm.CDET.hit.zgcell", &Earm_CDET_hit_zgcell, &b_Earm_CDET_hit_zgcell);
   fChain->SetBranchAddress("Earm.CDET.hit.NumPhotoelectrons", &Earm_CDET_hit_NumPhotoelectrons, &b_Earm_CDET_hit_NumPhotoelectrons);
   fChain->SetBranchAddress("Earm.CDET.hit.Time_avg", &Earm_CDET_hit_Time_avg, &b_Earm_CDET_hit_Time_avg);
   fChain->SetBranchAddress("Earm.CDET.hit.Time_rms", &Earm_CDET_hit_Time_rms, &b_Earm_CDET_hit_Time_rms);
   fChain->SetBranchAddress("Earm.CDET.hit.Time_min", &Earm_CDET_hit_Time_min, &b_Earm_CDET_hit_Time_min);
   fChain->SetBranchAddress("Earm.CDET.hit.Time_max", &Earm_CDET_hit_Time_max, &b_Earm_CDET_hit_Time_max);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.nhits", &Earm_CDET_Scint_hit_nhits, &b_Earm_CDET_Scint_hit_nhits);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.row", &Earm_CDET_Scint_hit_row, &b_Earm_CDET_Scint_hit_row);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.col", &Earm_CDET_Scint_hit_col, &b_Earm_CDET_Scint_hit_col);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.plane", &Earm_CDET_Scint_hit_plane, &b_Earm_CDET_Scint_hit_plane);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.xcell", &Earm_CDET_Scint_hit_xcell, &b_Earm_CDET_Scint_hit_xcell);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.ycell", &Earm_CDET_Scint_hit_ycell, &b_Earm_CDET_Scint_hit_ycell);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.zcell", &Earm_CDET_Scint_hit_zcell, &b_Earm_CDET_Scint_hit_zcell);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.xcellg", &Earm_CDET_Scint_hit_xcellg, &b_Earm_CDET_Scint_hit_xcellg);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.ycellg", &Earm_CDET_Scint_hit_ycellg, &b_Earm_CDET_Scint_hit_ycellg);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.zcellg", &Earm_CDET_Scint_hit_zcellg, &b_Earm_CDET_Scint_hit_zcellg);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.xhit", &Earm_CDET_Scint_hit_xhit, &b_Earm_CDET_Scint_hit_xhit);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.yhit", &Earm_CDET_Scint_hit_yhit, &b_Earm_CDET_Scint_hit_yhit);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.zhit", &Earm_CDET_Scint_hit_zhit, &b_Earm_CDET_Scint_hit_zhit);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.sumedep", &Earm_CDET_Scint_hit_sumedep, &b_Earm_CDET_Scint_hit_sumedep);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.tavg", &Earm_CDET_Scint_hit_tavg, &b_Earm_CDET_Scint_hit_tavg);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.trms", &Earm_CDET_Scint_hit_trms, &b_Earm_CDET_Scint_hit_trms);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.tmin", &Earm_CDET_Scint_hit_tmin, &b_Earm_CDET_Scint_hit_tmin);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.tmax", &Earm_CDET_Scint_hit_tmax, &b_Earm_CDET_Scint_hit_tmax);
   fChain->SetBranchAddress("Earm.CDET_Scint.part.npart", &Earm_CDET_Scint_part_npart, &b_Earm_CDET_Scint_part_npart);
   fChain->SetBranchAddress("Earm.CDET_Scint.part.PID", &Earm_CDET_Scint_part_PID, &b_Earm_CDET_Scint_part_PID);
   fChain->SetBranchAddress("Earm.CDET_Scint.part.MID", &Earm_CDET_Scint_part_MID, &b_Earm_CDET_Scint_part_MID);
   fChain->SetBranchAddress("Earm.CDET_Scint.part.TID", &Earm_CDET_Scint_part_TID, &b_Earm_CDET_Scint_part_TID);
   fChain->SetBranchAddress("Earm.CDET_Scint.part.nbounce", &Earm_CDET_Scint_part_nbounce, &b_Earm_CDET_Scint_part_nbounce);
   fChain->SetBranchAddress("Earm.CDET_Scint.part.hitindex", &Earm_CDET_Scint_part_hitindex, &b_Earm_CDET_Scint_part_hitindex);
   fChain->SetBranchAddress("Earm.CDET_Scint.part.vx", &Earm_CDET_Scint_part_vx, &b_Earm_CDET_Scint_part_vx);
   fChain->SetBranchAddress("Earm.CDET_Scint.part.vy", &Earm_CDET_Scint_part_vy, &b_Earm_CDET_Scint_part_vy);
   fChain->SetBranchAddress("Earm.CDET_Scint.part.vz", &Earm_CDET_Scint_part_vz, &b_Earm_CDET_Scint_part_vz);
   fChain->SetBranchAddress("Earm.CDET_Scint.part.px", &Earm_CDET_Scint_part_px, &b_Earm_CDET_Scint_part_px);
   fChain->SetBranchAddress("Earm.CDET_Scint.part.py", &Earm_CDET_Scint_part_py, &b_Earm_CDET_Scint_part_py);
   fChain->SetBranchAddress("Earm.CDET_Scint.part.pz", &Earm_CDET_Scint_part_pz, &b_Earm_CDET_Scint_part_pz);
   fChain->SetBranchAddress("Earm.ECAL.hit.nhits", &Earm_ECAL_hit_nhits, &b_Earm_ECAL_hit_nhits);
   fChain->SetBranchAddress("Earm.ECAL.hit.PMT", &Earm_ECAL_hit_PMT, &b_Earm_ECAL_hit_PMT);
   fChain->SetBranchAddress("Earm.ECAL.hit.row", &Earm_ECAL_hit_row, &b_Earm_ECAL_hit_row);
   fChain->SetBranchAddress("Earm.ECAL.hit.col", &Earm_ECAL_hit_col, &b_Earm_ECAL_hit_col);
   fChain->SetBranchAddress("Earm.ECAL.hit.plane", &Earm_ECAL_hit_plane, &b_Earm_ECAL_hit_plane);
   fChain->SetBranchAddress("Earm.ECAL.hit.xcell", &Earm_ECAL_hit_xcell, &b_Earm_ECAL_hit_xcell);
   fChain->SetBranchAddress("Earm.ECAL.hit.ycell", &Earm_ECAL_hit_ycell, &b_Earm_ECAL_hit_ycell);
   fChain->SetBranchAddress("Earm.ECAL.hit.zcell", &Earm_ECAL_hit_zcell, &b_Earm_ECAL_hit_zcell);
   fChain->SetBranchAddress("Earm.ECAL.hit.xgcell", &Earm_ECAL_hit_xgcell, &b_Earm_ECAL_hit_xgcell);
   fChain->SetBranchAddress("Earm.ECAL.hit.ygcell", &Earm_ECAL_hit_ygcell, &b_Earm_ECAL_hit_ygcell);
   fChain->SetBranchAddress("Earm.ECAL.hit.zgcell", &Earm_ECAL_hit_zgcell, &b_Earm_ECAL_hit_zgcell);
   fChain->SetBranchAddress("Earm.ECAL.hit.NumPhotoelectrons", &Earm_ECAL_hit_NumPhotoelectrons, &b_Earm_ECAL_hit_NumPhotoelectrons);
   fChain->SetBranchAddress("Earm.ECAL.hit.Time_avg", &Earm_ECAL_hit_Time_avg, &b_Earm_ECAL_hit_Time_avg);
   fChain->SetBranchAddress("Earm.ECAL.hit.Time_rms", &Earm_ECAL_hit_Time_rms, &b_Earm_ECAL_hit_Time_rms);
   fChain->SetBranchAddress("Earm.ECAL.hit.Time_min", &Earm_ECAL_hit_Time_min, &b_Earm_ECAL_hit_Time_min);
   fChain->SetBranchAddress("Earm.ECAL.hit.Time_max", &Earm_ECAL_hit_Time_max, &b_Earm_ECAL_hit_Time_max);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.nhits", &Earm_ECalTF1_hit_nhits, &b_Earm_ECalTF1_hit_nhits);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.row", &Earm_ECalTF1_hit_row, &b_Earm_ECalTF1_hit_row);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.col", &Earm_ECalTF1_hit_col, &b_Earm_ECalTF1_hit_col);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.plane", &Earm_ECalTF1_hit_plane, &b_Earm_ECalTF1_hit_plane);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.xcell", &Earm_ECalTF1_hit_xcell, &b_Earm_ECalTF1_hit_xcell);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.ycell", &Earm_ECalTF1_hit_ycell, &b_Earm_ECalTF1_hit_ycell);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.zcell", &Earm_ECalTF1_hit_zcell, &b_Earm_ECalTF1_hit_zcell);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.xcellg", &Earm_ECalTF1_hit_xcellg, &b_Earm_ECalTF1_hit_xcellg);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.ycellg", &Earm_ECalTF1_hit_ycellg, &b_Earm_ECalTF1_hit_ycellg);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.zcellg", &Earm_ECalTF1_hit_zcellg, &b_Earm_ECalTF1_hit_zcellg);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.xhit", &Earm_ECalTF1_hit_xhit, &b_Earm_ECalTF1_hit_xhit);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.yhit", &Earm_ECalTF1_hit_yhit, &b_Earm_ECalTF1_hit_yhit);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.zhit", &Earm_ECalTF1_hit_zhit, &b_Earm_ECalTF1_hit_zhit);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.sumedep", &Earm_ECalTF1_hit_sumedep, &b_Earm_ECalTF1_hit_sumedep);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.tavg", &Earm_ECalTF1_hit_tavg, &b_Earm_ECalTF1_hit_tavg);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.trms", &Earm_ECalTF1_hit_trms, &b_Earm_ECalTF1_hit_trms);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.tmin", &Earm_ECalTF1_hit_tmin, &b_Earm_ECalTF1_hit_tmin);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.tmax", &Earm_ECalTF1_hit_tmax, &b_Earm_ECalTF1_hit_tmax);
   fChain->SetBranchAddress("Earm.ECalTF1.part.npart", &Earm_ECalTF1_part_npart, &b_Earm_ECalTF1_part_npart);
   fChain->SetBranchAddress("Earm.ECalTF1.part.PID", &Earm_ECalTF1_part_PID, &b_Earm_ECalTF1_part_PID);
   fChain->SetBranchAddress("Earm.ECalTF1.part.MID", &Earm_ECalTF1_part_MID, &b_Earm_ECalTF1_part_MID);
   fChain->SetBranchAddress("Earm.ECalTF1.part.TID", &Earm_ECalTF1_part_TID, &b_Earm_ECalTF1_part_TID);
   fChain->SetBranchAddress("Earm.ECalTF1.part.nbounce", &Earm_ECalTF1_part_nbounce, &b_Earm_ECalTF1_part_nbounce);
   fChain->SetBranchAddress("Earm.ECalTF1.part.hitindex", &Earm_ECalTF1_part_hitindex, &b_Earm_ECalTF1_part_hitindex);
   fChain->SetBranchAddress("Earm.ECalTF1.part.vx", &Earm_ECalTF1_part_vx, &b_Earm_ECalTF1_part_vx);
   fChain->SetBranchAddress("Earm.ECalTF1.part.vy", &Earm_ECalTF1_part_vy, &b_Earm_ECalTF1_part_vy);
   fChain->SetBranchAddress("Earm.ECalTF1.part.vz", &Earm_ECalTF1_part_vz, &b_Earm_ECalTF1_part_vz);
   fChain->SetBranchAddress("Earm.ECalTF1.part.px", &Earm_ECalTF1_part_px, &b_Earm_ECalTF1_part_px);
   fChain->SetBranchAddress("Earm.ECalTF1.part.py", &Earm_ECalTF1_part_py, &b_Earm_ECalTF1_part_py);
   fChain->SetBranchAddress("Earm.ECalTF1.part.pz", &Earm_ECalTF1_part_pz, &b_Earm_ECalTF1_part_pz);
   fChain->SetBranchAddress("Harm.FPP1.hit.nhits", &Harm_FPP1_hit_nhits, &b_Harm_FPP1_hit_nhits);
   fChain->SetBranchAddress("Harm.FPP1.hit.plane", &Harm_FPP1_hit_plane, &b_Harm_FPP1_hit_plane);
   fChain->SetBranchAddress("Harm.FPP1.hit.strip", &Harm_FPP1_hit_strip, &b_Harm_FPP1_hit_strip);
   fChain->SetBranchAddress("Harm.FPP1.hit.x", &Harm_FPP1_hit_x, &b_Harm_FPP1_hit_x);
   fChain->SetBranchAddress("Harm.FPP1.hit.y", &Harm_FPP1_hit_y, &b_Harm_FPP1_hit_y);
   fChain->SetBranchAddress("Harm.FPP1.hit.z", &Harm_FPP1_hit_z, &b_Harm_FPP1_hit_z);
   fChain->SetBranchAddress("Harm.FPP1.hit.t", &Harm_FPP1_hit_t, &b_Harm_FPP1_hit_t);
   fChain->SetBranchAddress("Harm.FPP1.hit.trms", &Harm_FPP1_hit_trms, &b_Harm_FPP1_hit_trms);
   fChain->SetBranchAddress("Harm.FPP1.hit.tmin", &Harm_FPP1_hit_tmin, &b_Harm_FPP1_hit_tmin);
   fChain->SetBranchAddress("Harm.FPP1.hit.tmax", &Harm_FPP1_hit_tmax, &b_Harm_FPP1_hit_tmax);
   fChain->SetBranchAddress("Harm.FPP1.hit.tx", &Harm_FPP1_hit_tx, &b_Harm_FPP1_hit_tx);
   fChain->SetBranchAddress("Harm.FPP1.hit.ty", &Harm_FPP1_hit_ty, &b_Harm_FPP1_hit_ty);
   fChain->SetBranchAddress("Harm.FPP1.hit.txp", &Harm_FPP1_hit_txp, &b_Harm_FPP1_hit_txp);
   fChain->SetBranchAddress("Harm.FPP1.hit.typ", &Harm_FPP1_hit_typ, &b_Harm_FPP1_hit_typ);
   fChain->SetBranchAddress("Harm.FPP1.hit.trid", &Harm_FPP1_hit_trid, &b_Harm_FPP1_hit_trid);
   fChain->SetBranchAddress("Harm.FPP1.hit.mid", &Harm_FPP1_hit_mid, &b_Harm_FPP1_hit_mid);
   fChain->SetBranchAddress("Harm.FPP1.hit.pid", &Harm_FPP1_hit_pid, &b_Harm_FPP1_hit_pid);
   fChain->SetBranchAddress("Harm.FPP1.hit.vx", &Harm_FPP1_hit_vx, &b_Harm_FPP1_hit_vx);
   fChain->SetBranchAddress("Harm.FPP1.hit.vy", &Harm_FPP1_hit_vy, &b_Harm_FPP1_hit_vy);
   fChain->SetBranchAddress("Harm.FPP1.hit.vz", &Harm_FPP1_hit_vz, &b_Harm_FPP1_hit_vz);
   fChain->SetBranchAddress("Harm.FPP1.hit.p", &Harm_FPP1_hit_p, &b_Harm_FPP1_hit_p);
   fChain->SetBranchAddress("Harm.FPP1.hit.edep", &Harm_FPP1_hit_edep, &b_Harm_FPP1_hit_edep);
   fChain->SetBranchAddress("Harm.FPP1.hit.beta", &Harm_FPP1_hit_beta, &b_Harm_FPP1_hit_beta);
   fChain->SetBranchAddress("Harm.FPP1.Track.ntracks", &Harm_FPP1_Track_ntracks, &b_Harm_FPP1_Track_ntracks);
   fChain->SetBranchAddress("Harm.FPP1.Track.TID", &Harm_FPP1_Track_TID, &b_Harm_FPP1_Track_TID);
   fChain->SetBranchAddress("Harm.FPP1.Track.PID", &Harm_FPP1_Track_PID, &b_Harm_FPP1_Track_PID);
   fChain->SetBranchAddress("Harm.FPP1.Track.MID", &Harm_FPP1_Track_MID, &b_Harm_FPP1_Track_MID);
   fChain->SetBranchAddress("Harm.FPP1.Track.NumHits", &Harm_FPP1_Track_NumHits, &b_Harm_FPP1_Track_NumHits);
   fChain->SetBranchAddress("Harm.FPP1.Track.NumPlanes", &Harm_FPP1_Track_NumPlanes, &b_Harm_FPP1_Track_NumPlanes);
   fChain->SetBranchAddress("Harm.FPP1.Track.NDF", &Harm_FPP1_Track_NDF, &b_Harm_FPP1_Track_NDF);
   fChain->SetBranchAddress("Harm.FPP1.Track.Chi2fit", &Harm_FPP1_Track_Chi2fit, &b_Harm_FPP1_Track_Chi2fit);
   fChain->SetBranchAddress("Harm.FPP1.Track.Chi2true", &Harm_FPP1_Track_Chi2true, &b_Harm_FPP1_Track_Chi2true);
   fChain->SetBranchAddress("Harm.FPP1.Track.X", &Harm_FPP1_Track_X, &b_Harm_FPP1_Track_X);
   fChain->SetBranchAddress("Harm.FPP1.Track.Y", &Harm_FPP1_Track_Y, &b_Harm_FPP1_Track_Y);
   fChain->SetBranchAddress("Harm.FPP1.Track.Xp", &Harm_FPP1_Track_Xp, &b_Harm_FPP1_Track_Xp);
   fChain->SetBranchAddress("Harm.FPP1.Track.Yp", &Harm_FPP1_Track_Yp, &b_Harm_FPP1_Track_Yp);
   fChain->SetBranchAddress("Harm.FPP1.Track.T", &Harm_FPP1_Track_T, &b_Harm_FPP1_Track_T);
   fChain->SetBranchAddress("Harm.FPP1.Track.P", &Harm_FPP1_Track_P, &b_Harm_FPP1_Track_P);
   fChain->SetBranchAddress("Harm.FPP1.Track.Xfit", &Harm_FPP1_Track_Xfit, &b_Harm_FPP1_Track_Xfit);
   fChain->SetBranchAddress("Harm.FPP1.Track.Yfit", &Harm_FPP1_Track_Yfit, &b_Harm_FPP1_Track_Yfit);
   fChain->SetBranchAddress("Harm.FPP1.Track.Xpfit", &Harm_FPP1_Track_Xpfit, &b_Harm_FPP1_Track_Xpfit);
   fChain->SetBranchAddress("Harm.FPP1.Track.Ypfit", &Harm_FPP1_Track_Ypfit, &b_Harm_FPP1_Track_Ypfit);
   fChain->SetBranchAddress("Harm.FPP1.part.npart", &Harm_FPP1_part_npart, &b_Harm_FPP1_part_npart);
   fChain->SetBranchAddress("Harm.FPP1.part.PID", &Harm_FPP1_part_PID, &b_Harm_FPP1_part_PID);
   fChain->SetBranchAddress("Harm.FPP1.part.MID", &Harm_FPP1_part_MID, &b_Harm_FPP1_part_MID);
   fChain->SetBranchAddress("Harm.FPP1.part.TID", &Harm_FPP1_part_TID, &b_Harm_FPP1_part_TID);
   fChain->SetBranchAddress("Harm.FPP1.part.nbounce", &Harm_FPP1_part_nbounce, &b_Harm_FPP1_part_nbounce);
   fChain->SetBranchAddress("Harm.FPP1.part.hitindex", &Harm_FPP1_part_hitindex, &b_Harm_FPP1_part_hitindex);
   fChain->SetBranchAddress("Harm.FPP1.part.vx", &Harm_FPP1_part_vx, &b_Harm_FPP1_part_vx);
   fChain->SetBranchAddress("Harm.FPP1.part.vy", &Harm_FPP1_part_vy, &b_Harm_FPP1_part_vy);
   fChain->SetBranchAddress("Harm.FPP1.part.vz", &Harm_FPP1_part_vz, &b_Harm_FPP1_part_vz);
   fChain->SetBranchAddress("Harm.FPP1.part.px", &Harm_FPP1_part_px, &b_Harm_FPP1_part_px);
   fChain->SetBranchAddress("Harm.FPP1.part.py", &Harm_FPP1_part_py, &b_Harm_FPP1_part_py);
   fChain->SetBranchAddress("Harm.FPP1.part.pz", &Harm_FPP1_part_pz, &b_Harm_FPP1_part_pz);
   fChain->SetBranchAddress("Harm.FPP2.hit.nhits", &Harm_FPP2_hit_nhits, &b_Harm_FPP2_hit_nhits);
   fChain->SetBranchAddress("Harm.FPP2.hit.plane", &Harm_FPP2_hit_plane, &b_Harm_FPP2_hit_plane);
   fChain->SetBranchAddress("Harm.FPP2.hit.strip", &Harm_FPP2_hit_strip, &b_Harm_FPP2_hit_strip);
   fChain->SetBranchAddress("Harm.FPP2.hit.x", &Harm_FPP2_hit_x, &b_Harm_FPP2_hit_x);
   fChain->SetBranchAddress("Harm.FPP2.hit.y", &Harm_FPP2_hit_y, &b_Harm_FPP2_hit_y);
   fChain->SetBranchAddress("Harm.FPP2.hit.z", &Harm_FPP2_hit_z, &b_Harm_FPP2_hit_z);
   fChain->SetBranchAddress("Harm.FPP2.hit.t", &Harm_FPP2_hit_t, &b_Harm_FPP2_hit_t);
   fChain->SetBranchAddress("Harm.FPP2.hit.trms", &Harm_FPP2_hit_trms, &b_Harm_FPP2_hit_trms);
   fChain->SetBranchAddress("Harm.FPP2.hit.tmin", &Harm_FPP2_hit_tmin, &b_Harm_FPP2_hit_tmin);
   fChain->SetBranchAddress("Harm.FPP2.hit.tmax", &Harm_FPP2_hit_tmax, &b_Harm_FPP2_hit_tmax);
   fChain->SetBranchAddress("Harm.FPP2.hit.tx", &Harm_FPP2_hit_tx, &b_Harm_FPP2_hit_tx);
   fChain->SetBranchAddress("Harm.FPP2.hit.ty", &Harm_FPP2_hit_ty, &b_Harm_FPP2_hit_ty);
   fChain->SetBranchAddress("Harm.FPP2.hit.txp", &Harm_FPP2_hit_txp, &b_Harm_FPP2_hit_txp);
   fChain->SetBranchAddress("Harm.FPP2.hit.typ", &Harm_FPP2_hit_typ, &b_Harm_FPP2_hit_typ);
   fChain->SetBranchAddress("Harm.FPP2.hit.trid", &Harm_FPP2_hit_trid, &b_Harm_FPP2_hit_trid);
   fChain->SetBranchAddress("Harm.FPP2.hit.mid", &Harm_FPP2_hit_mid, &b_Harm_FPP2_hit_mid);
   fChain->SetBranchAddress("Harm.FPP2.hit.pid", &Harm_FPP2_hit_pid, &b_Harm_FPP2_hit_pid);
   fChain->SetBranchAddress("Harm.FPP2.hit.vx", &Harm_FPP2_hit_vx, &b_Harm_FPP2_hit_vx);
   fChain->SetBranchAddress("Harm.FPP2.hit.vy", &Harm_FPP2_hit_vy, &b_Harm_FPP2_hit_vy);
   fChain->SetBranchAddress("Harm.FPP2.hit.vz", &Harm_FPP2_hit_vz, &b_Harm_FPP2_hit_vz);
   fChain->SetBranchAddress("Harm.FPP2.hit.p", &Harm_FPP2_hit_p, &b_Harm_FPP2_hit_p);
   fChain->SetBranchAddress("Harm.FPP2.hit.edep", &Harm_FPP2_hit_edep, &b_Harm_FPP2_hit_edep);
   fChain->SetBranchAddress("Harm.FPP2.hit.beta", &Harm_FPP2_hit_beta, &b_Harm_FPP2_hit_beta);
   fChain->SetBranchAddress("Harm.FPP2.Track.ntracks", &Harm_FPP2_Track_ntracks, &b_Harm_FPP2_Track_ntracks);
   fChain->SetBranchAddress("Harm.FPP2.Track.TID", &Harm_FPP2_Track_TID, &b_Harm_FPP2_Track_TID);
   fChain->SetBranchAddress("Harm.FPP2.Track.PID", &Harm_FPP2_Track_PID, &b_Harm_FPP2_Track_PID);
   fChain->SetBranchAddress("Harm.FPP2.Track.MID", &Harm_FPP2_Track_MID, &b_Harm_FPP2_Track_MID);
   fChain->SetBranchAddress("Harm.FPP2.Track.NumHits", &Harm_FPP2_Track_NumHits, &b_Harm_FPP2_Track_NumHits);
   fChain->SetBranchAddress("Harm.FPP2.Track.NumPlanes", &Harm_FPP2_Track_NumPlanes, &b_Harm_FPP2_Track_NumPlanes);
   fChain->SetBranchAddress("Harm.FPP2.Track.NDF", &Harm_FPP2_Track_NDF, &b_Harm_FPP2_Track_NDF);
   fChain->SetBranchAddress("Harm.FPP2.Track.Chi2fit", &Harm_FPP2_Track_Chi2fit, &b_Harm_FPP2_Track_Chi2fit);
   fChain->SetBranchAddress("Harm.FPP2.Track.Chi2true", &Harm_FPP2_Track_Chi2true, &b_Harm_FPP2_Track_Chi2true);
   fChain->SetBranchAddress("Harm.FPP2.Track.X", &Harm_FPP2_Track_X, &b_Harm_FPP2_Track_X);
   fChain->SetBranchAddress("Harm.FPP2.Track.Y", &Harm_FPP2_Track_Y, &b_Harm_FPP2_Track_Y);
   fChain->SetBranchAddress("Harm.FPP2.Track.Xp", &Harm_FPP2_Track_Xp, &b_Harm_FPP2_Track_Xp);
   fChain->SetBranchAddress("Harm.FPP2.Track.Yp", &Harm_FPP2_Track_Yp, &b_Harm_FPP2_Track_Yp);
   fChain->SetBranchAddress("Harm.FPP2.Track.T", &Harm_FPP2_Track_T, &b_Harm_FPP2_Track_T);
   fChain->SetBranchAddress("Harm.FPP2.Track.P", &Harm_FPP2_Track_P, &b_Harm_FPP2_Track_P);
   fChain->SetBranchAddress("Harm.FPP2.Track.Xfit", &Harm_FPP2_Track_Xfit, &b_Harm_FPP2_Track_Xfit);
   fChain->SetBranchAddress("Harm.FPP2.Track.Yfit", &Harm_FPP2_Track_Yfit, &b_Harm_FPP2_Track_Yfit);
   fChain->SetBranchAddress("Harm.FPP2.Track.Xpfit", &Harm_FPP2_Track_Xpfit, &b_Harm_FPP2_Track_Xpfit);
   fChain->SetBranchAddress("Harm.FPP2.Track.Ypfit", &Harm_FPP2_Track_Ypfit, &b_Harm_FPP2_Track_Ypfit);
   fChain->SetBranchAddress("Harm.FPP2.part.npart", &Harm_FPP2_part_npart, &b_Harm_FPP2_part_npart);
   fChain->SetBranchAddress("Harm.FPP2.part.PID", &Harm_FPP2_part_PID, &b_Harm_FPP2_part_PID);
   fChain->SetBranchAddress("Harm.FPP2.part.MID", &Harm_FPP2_part_MID, &b_Harm_FPP2_part_MID);
   fChain->SetBranchAddress("Harm.FPP2.part.TID", &Harm_FPP2_part_TID, &b_Harm_FPP2_part_TID);
   fChain->SetBranchAddress("Harm.FPP2.part.nbounce", &Harm_FPP2_part_nbounce, &b_Harm_FPP2_part_nbounce);
   fChain->SetBranchAddress("Harm.FPP2.part.hitindex", &Harm_FPP2_part_hitindex, &b_Harm_FPP2_part_hitindex);
   fChain->SetBranchAddress("Harm.FPP2.part.vx", &Harm_FPP2_part_vx, &b_Harm_FPP2_part_vx);
   fChain->SetBranchAddress("Harm.FPP2.part.vy", &Harm_FPP2_part_vy, &b_Harm_FPP2_part_vy);
   fChain->SetBranchAddress("Harm.FPP2.part.vz", &Harm_FPP2_part_vz, &b_Harm_FPP2_part_vz);
   fChain->SetBranchAddress("Harm.FPP2.part.px", &Harm_FPP2_part_px, &b_Harm_FPP2_part_px);
   fChain->SetBranchAddress("Harm.FPP2.part.py", &Harm_FPP2_part_py, &b_Harm_FPP2_part_py);
   fChain->SetBranchAddress("Harm.FPP2.part.pz", &Harm_FPP2_part_pz, &b_Harm_FPP2_part_pz);
   fChain->SetBranchAddress("Harm.FT.hit.nhits", &Harm_FT_hit_nhits, &b_Harm_FT_hit_nhits);
   fChain->SetBranchAddress("Harm.FT.hit.plane", &Harm_FT_hit_plane, &b_Harm_FT_hit_plane);
   fChain->SetBranchAddress("Harm.FT.hit.strip", &Harm_FT_hit_strip, &b_Harm_FT_hit_strip);
   fChain->SetBranchAddress("Harm.FT.hit.x", &Harm_FT_hit_x, &b_Harm_FT_hit_x);
   fChain->SetBranchAddress("Harm.FT.hit.y", &Harm_FT_hit_y, &b_Harm_FT_hit_y);
   fChain->SetBranchAddress("Harm.FT.hit.z", &Harm_FT_hit_z, &b_Harm_FT_hit_z);
   fChain->SetBranchAddress("Harm.FT.hit.t", &Harm_FT_hit_t, &b_Harm_FT_hit_t);
   fChain->SetBranchAddress("Harm.FT.hit.trms", &Harm_FT_hit_trms, &b_Harm_FT_hit_trms);
   fChain->SetBranchAddress("Harm.FT.hit.tmin", &Harm_FT_hit_tmin, &b_Harm_FT_hit_tmin);
   fChain->SetBranchAddress("Harm.FT.hit.tmax", &Harm_FT_hit_tmax, &b_Harm_FT_hit_tmax);
   fChain->SetBranchAddress("Harm.FT.hit.tx", &Harm_FT_hit_tx, &b_Harm_FT_hit_tx);
   fChain->SetBranchAddress("Harm.FT.hit.ty", &Harm_FT_hit_ty, &b_Harm_FT_hit_ty);
   fChain->SetBranchAddress("Harm.FT.hit.txp", &Harm_FT_hit_txp, &b_Harm_FT_hit_txp);
   fChain->SetBranchAddress("Harm.FT.hit.typ", &Harm_FT_hit_typ, &b_Harm_FT_hit_typ);
   fChain->SetBranchAddress("Harm.FT.hit.trid", &Harm_FT_hit_trid, &b_Harm_FT_hit_trid);
   fChain->SetBranchAddress("Harm.FT.hit.mid", &Harm_FT_hit_mid, &b_Harm_FT_hit_mid);
   fChain->SetBranchAddress("Harm.FT.hit.pid", &Harm_FT_hit_pid, &b_Harm_FT_hit_pid);
   fChain->SetBranchAddress("Harm.FT.hit.vx", &Harm_FT_hit_vx, &b_Harm_FT_hit_vx);
   fChain->SetBranchAddress("Harm.FT.hit.vy", &Harm_FT_hit_vy, &b_Harm_FT_hit_vy);
   fChain->SetBranchAddress("Harm.FT.hit.vz", &Harm_FT_hit_vz, &b_Harm_FT_hit_vz);
   fChain->SetBranchAddress("Harm.FT.hit.p", &Harm_FT_hit_p, &b_Harm_FT_hit_p);
   fChain->SetBranchAddress("Harm.FT.hit.edep", &Harm_FT_hit_edep, &b_Harm_FT_hit_edep);
   fChain->SetBranchAddress("Harm.FT.hit.beta", &Harm_FT_hit_beta, &b_Harm_FT_hit_beta);
   fChain->SetBranchAddress("Harm.FT.Track.ntracks", &Harm_FT_Track_ntracks, &b_Harm_FT_Track_ntracks);
   fChain->SetBranchAddress("Harm.FT.Track.TID", &Harm_FT_Track_TID, &b_Harm_FT_Track_TID);
   fChain->SetBranchAddress("Harm.FT.Track.PID", &Harm_FT_Track_PID, &b_Harm_FT_Track_PID);
   fChain->SetBranchAddress("Harm.FT.Track.MID", &Harm_FT_Track_MID, &b_Harm_FT_Track_MID);
   fChain->SetBranchAddress("Harm.FT.Track.NumHits", &Harm_FT_Track_NumHits, &b_Harm_FT_Track_NumHits);
   fChain->SetBranchAddress("Harm.FT.Track.NumPlanes", &Harm_FT_Track_NumPlanes, &b_Harm_FT_Track_NumPlanes);
   fChain->SetBranchAddress("Harm.FT.Track.NDF", &Harm_FT_Track_NDF, &b_Harm_FT_Track_NDF);
   fChain->SetBranchAddress("Harm.FT.Track.Chi2fit", &Harm_FT_Track_Chi2fit, &b_Harm_FT_Track_Chi2fit);
   fChain->SetBranchAddress("Harm.FT.Track.Chi2true", &Harm_FT_Track_Chi2true, &b_Harm_FT_Track_Chi2true);
   fChain->SetBranchAddress("Harm.FT.Track.X", &Harm_FT_Track_X, &b_Harm_FT_Track_X);
   fChain->SetBranchAddress("Harm.FT.Track.Y", &Harm_FT_Track_Y, &b_Harm_FT_Track_Y);
   fChain->SetBranchAddress("Harm.FT.Track.Xp", &Harm_FT_Track_Xp, &b_Harm_FT_Track_Xp);
   fChain->SetBranchAddress("Harm.FT.Track.Yp", &Harm_FT_Track_Yp, &b_Harm_FT_Track_Yp);
   fChain->SetBranchAddress("Harm.FT.Track.T", &Harm_FT_Track_T, &b_Harm_FT_Track_T);
   fChain->SetBranchAddress("Harm.FT.Track.P", &Harm_FT_Track_P, &b_Harm_FT_Track_P);
   fChain->SetBranchAddress("Harm.FT.Track.Xfit", &Harm_FT_Track_Xfit, &b_Harm_FT_Track_Xfit);
   fChain->SetBranchAddress("Harm.FT.Track.Yfit", &Harm_FT_Track_Yfit, &b_Harm_FT_Track_Yfit);
   fChain->SetBranchAddress("Harm.FT.Track.Xpfit", &Harm_FT_Track_Xpfit, &b_Harm_FT_Track_Xpfit);
   fChain->SetBranchAddress("Harm.FT.Track.Ypfit", &Harm_FT_Track_Ypfit, &b_Harm_FT_Track_Ypfit);
   fChain->SetBranchAddress("Harm.FT.part.npart", &Harm_FT_part_npart, &b_Harm_FT_part_npart);
   fChain->SetBranchAddress("Harm.FT.part.PID", &Harm_FT_part_PID, &b_Harm_FT_part_PID);
   fChain->SetBranchAddress("Harm.FT.part.MID", &Harm_FT_part_MID, &b_Harm_FT_part_MID);
   fChain->SetBranchAddress("Harm.FT.part.TID", &Harm_FT_part_TID, &b_Harm_FT_part_TID);
   fChain->SetBranchAddress("Harm.FT.part.nbounce", &Harm_FT_part_nbounce, &b_Harm_FT_part_nbounce);
   fChain->SetBranchAddress("Harm.FT.part.hitindex", &Harm_FT_part_hitindex, &b_Harm_FT_part_hitindex);
   fChain->SetBranchAddress("Harm.FT.part.vx", &Harm_FT_part_vx, &b_Harm_FT_part_vx);
   fChain->SetBranchAddress("Harm.FT.part.vy", &Harm_FT_part_vy, &b_Harm_FT_part_vy);
   fChain->SetBranchAddress("Harm.FT.part.vz", &Harm_FT_part_vz, &b_Harm_FT_part_vz);
   fChain->SetBranchAddress("Harm.FT.part.px", &Harm_FT_part_px, &b_Harm_FT_part_px);
   fChain->SetBranchAddress("Harm.FT.part.py", &Harm_FT_part_py, &b_Harm_FT_part_py);
   fChain->SetBranchAddress("Harm.FT.part.pz", &Harm_FT_part_pz, &b_Harm_FT_part_pz);
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
   fChain->SetBranchAddress("Harm.HCalScint.part.npart", &Harm_HCalScint_part_npart, &b_Harm_HCalScint_part_npart);
   fChain->SetBranchAddress("Harm.HCalScint.part.PID", &Harm_HCalScint_part_PID, &b_Harm_HCalScint_part_PID);
   fChain->SetBranchAddress("Harm.HCalScint.part.MID", &Harm_HCalScint_part_MID, &b_Harm_HCalScint_part_MID);
   fChain->SetBranchAddress("Harm.HCalScint.part.TID", &Harm_HCalScint_part_TID, &b_Harm_HCalScint_part_TID);
   fChain->SetBranchAddress("Harm.HCalScint.part.nbounce", &Harm_HCalScint_part_nbounce, &b_Harm_HCalScint_part_nbounce);
   fChain->SetBranchAddress("Harm.HCalScint.part.hitindex", &Harm_HCalScint_part_hitindex, &b_Harm_HCalScint_part_hitindex);
   fChain->SetBranchAddress("Harm.HCalScint.part.vx", &Harm_HCalScint_part_vx, &b_Harm_HCalScint_part_vx);
   fChain->SetBranchAddress("Harm.HCalScint.part.vy", &Harm_HCalScint_part_vy, &b_Harm_HCalScint_part_vy);
   fChain->SetBranchAddress("Harm.HCalScint.part.vz", &Harm_HCalScint_part_vz, &b_Harm_HCalScint_part_vz);
   fChain->SetBranchAddress("Harm.HCalScint.part.px", &Harm_HCalScint_part_px, &b_Harm_HCalScint_part_px);
   fChain->SetBranchAddress("Harm.HCalScint.part.py", &Harm_HCalScint_part_py, &b_Harm_HCalScint_part_py);
   fChain->SetBranchAddress("Harm.HCalScint.part.pz", &Harm_HCalScint_part_pz, &b_Harm_HCalScint_part_pz);
   Notify();
}

Bool_t gep_tree_July2015::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gep_tree_July2015::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gep_tree_July2015::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gep_tree_July2015_cxx
