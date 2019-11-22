//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 24 12:15:03 2019 by ROOT version 6.16/00
// from TTree T/Geant4 SBS Simulation
// found on file: gep_12GeV2_pythia6_opticalphotons_1kevents.root
//////////////////////////////////////////////////////////

#ifndef gep_pythia6_tree_h
#define gep_pythia6_tree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class gep_pythia6_tree {
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
   Double_t        Earm_CDET_Scint_det_esum;
   Int_t           Earm_CDET_Scint_hit_nhits;
   vector<int>     *Earm_CDET_Scint_hit_row;
   vector<int>     *Earm_CDET_Scint_hit_col;
   vector<int>     *Earm_CDET_Scint_hit_cell;
   vector<int>     *Earm_CDET_Scint_hit_plane;
   vector<int>     *Earm_CDET_Scint_hit_wire;
   vector<double>  *Earm_CDET_Scint_hit_xcell;
   vector<double>  *Earm_CDET_Scint_hit_ycell;
   vector<double>  *Earm_CDET_Scint_hit_zcell;
   vector<double>  *Earm_CDET_Scint_hit_xcellg;
   vector<double>  *Earm_CDET_Scint_hit_ycellg;
   vector<double>  *Earm_CDET_Scint_hit_zcellg;
   vector<double>  *Earm_CDET_Scint_hit_xhit;
   vector<double>  *Earm_CDET_Scint_hit_yhit;
   vector<double>  *Earm_CDET_Scint_hit_zhit;
   vector<double>  *Earm_CDET_Scint_hit_xhitg;
   vector<double>  *Earm_CDET_Scint_hit_yhitg;
   vector<double>  *Earm_CDET_Scint_hit_zhitg;
   vector<double>  *Earm_CDET_Scint_hit_sumedep;
   vector<double>  *Earm_CDET_Scint_hit_tavg;
   vector<double>  *Earm_CDET_Scint_hit_trms;
   vector<double>  *Earm_CDET_Scint_hit_tmin;
   vector<double>  *Earm_CDET_Scint_hit_tmax;
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
   Double_t        Earm_ECalTF1_det_esum;
   Int_t           Earm_ECalTF1_hit_nhits;
   vector<int>     *Earm_ECalTF1_hit_row;
   vector<int>     *Earm_ECalTF1_hit_col;
   vector<int>     *Earm_ECalTF1_hit_cell;
   vector<int>     *Earm_ECalTF1_hit_plane;
   vector<int>     *Earm_ECalTF1_hit_wire;
   vector<double>  *Earm_ECalTF1_hit_xcell;
   vector<double>  *Earm_ECalTF1_hit_ycell;
   vector<double>  *Earm_ECalTF1_hit_zcell;
   vector<double>  *Earm_ECalTF1_hit_xcellg;
   vector<double>  *Earm_ECalTF1_hit_ycellg;
   vector<double>  *Earm_ECalTF1_hit_zcellg;
   vector<double>  *Earm_ECalTF1_hit_xhit;
   vector<double>  *Earm_ECalTF1_hit_yhit;
   vector<double>  *Earm_ECalTF1_hit_zhit;
   vector<double>  *Earm_ECalTF1_hit_xhitg;
   vector<double>  *Earm_ECalTF1_hit_yhitg;
   vector<double>  *Earm_ECalTF1_hit_zhitg;
   vector<double>  *Earm_ECalTF1_hit_sumedep;
   vector<double>  *Earm_ECalTF1_hit_tavg;
   vector<double>  *Earm_ECalTF1_hit_trms;
   vector<double>  *Earm_ECalTF1_hit_tmin;
   vector<double>  *Earm_ECalTF1_hit_tmax;
   Int_t           Harm_FPP1_hit_nhits;
   vector<int>     *Harm_FPP1_hit_plane;
   vector<int>     *Harm_FPP1_hit_strip;
   vector<double>  *Harm_FPP1_hit_x;
   vector<double>  *Harm_FPP1_hit_y;
   vector<double>  *Harm_FPP1_hit_z;
   vector<double>  *Harm_FPP1_hit_polx;
   vector<double>  *Harm_FPP1_hit_poly;
   vector<double>  *Harm_FPP1_hit_polz;
   vector<double>  *Harm_FPP1_hit_t;
   vector<double>  *Harm_FPP1_hit_trms;
   vector<double>  *Harm_FPP1_hit_tmin;
   vector<double>  *Harm_FPP1_hit_tmax;
   vector<double>  *Harm_FPP1_hit_tx;
   vector<double>  *Harm_FPP1_hit_ty;
   vector<double>  *Harm_FPP1_hit_xin;
   vector<double>  *Harm_FPP1_hit_yin;
   vector<double>  *Harm_FPP1_hit_zin;
   vector<double>  *Harm_FPP1_hit_xout;
   vector<double>  *Harm_FPP1_hit_yout;
   vector<double>  *Harm_FPP1_hit_zout;
   vector<double>  *Harm_FPP1_hit_txp;
   vector<double>  *Harm_FPP1_hit_typ;
   vector<double>  *Harm_FPP1_hit_xg;
   vector<double>  *Harm_FPP1_hit_yg;
   vector<double>  *Harm_FPP1_hit_zg;
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
   vector<double>  *Harm_FPP1_Track_Sx;
   vector<double>  *Harm_FPP1_Track_Sy;
   vector<double>  *Harm_FPP1_Track_Sz;
   vector<double>  *Harm_FPP1_Track_Xfit;
   vector<double>  *Harm_FPP1_Track_Yfit;
   vector<double>  *Harm_FPP1_Track_Xpfit;
   vector<double>  *Harm_FPP1_Track_Ypfit;
   Int_t           Harm_FPP2_hit_nhits;
   vector<int>     *Harm_FPP2_hit_plane;
   vector<int>     *Harm_FPP2_hit_strip;
   vector<double>  *Harm_FPP2_hit_x;
   vector<double>  *Harm_FPP2_hit_y;
   vector<double>  *Harm_FPP2_hit_z;
   vector<double>  *Harm_FPP2_hit_polx;
   vector<double>  *Harm_FPP2_hit_poly;
   vector<double>  *Harm_FPP2_hit_polz;
   vector<double>  *Harm_FPP2_hit_t;
   vector<double>  *Harm_FPP2_hit_trms;
   vector<double>  *Harm_FPP2_hit_tmin;
   vector<double>  *Harm_FPP2_hit_tmax;
   vector<double>  *Harm_FPP2_hit_tx;
   vector<double>  *Harm_FPP2_hit_ty;
   vector<double>  *Harm_FPP2_hit_xin;
   vector<double>  *Harm_FPP2_hit_yin;
   vector<double>  *Harm_FPP2_hit_zin;
   vector<double>  *Harm_FPP2_hit_xout;
   vector<double>  *Harm_FPP2_hit_yout;
   vector<double>  *Harm_FPP2_hit_zout;
   vector<double>  *Harm_FPP2_hit_txp;
   vector<double>  *Harm_FPP2_hit_typ;
   vector<double>  *Harm_FPP2_hit_xg;
   vector<double>  *Harm_FPP2_hit_yg;
   vector<double>  *Harm_FPP2_hit_zg;
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
   vector<double>  *Harm_FPP2_Track_Sx;
   vector<double>  *Harm_FPP2_Track_Sy;
   vector<double>  *Harm_FPP2_Track_Sz;
   vector<double>  *Harm_FPP2_Track_Xfit;
   vector<double>  *Harm_FPP2_Track_Yfit;
   vector<double>  *Harm_FPP2_Track_Xpfit;
   vector<double>  *Harm_FPP2_Track_Ypfit;
   Int_t           Harm_FT_hit_nhits;
   vector<int>     *Harm_FT_hit_plane;
   vector<int>     *Harm_FT_hit_strip;
   vector<double>  *Harm_FT_hit_x;
   vector<double>  *Harm_FT_hit_y;
   vector<double>  *Harm_FT_hit_z;
   vector<double>  *Harm_FT_hit_polx;
   vector<double>  *Harm_FT_hit_poly;
   vector<double>  *Harm_FT_hit_polz;
   vector<double>  *Harm_FT_hit_t;
   vector<double>  *Harm_FT_hit_trms;
   vector<double>  *Harm_FT_hit_tmin;
   vector<double>  *Harm_FT_hit_tmax;
   vector<double>  *Harm_FT_hit_tx;
   vector<double>  *Harm_FT_hit_ty;
   vector<double>  *Harm_FT_hit_xin;
   vector<double>  *Harm_FT_hit_yin;
   vector<double>  *Harm_FT_hit_zin;
   vector<double>  *Harm_FT_hit_xout;
   vector<double>  *Harm_FT_hit_yout;
   vector<double>  *Harm_FT_hit_zout;
   vector<double>  *Harm_FT_hit_txp;
   vector<double>  *Harm_FT_hit_typ;
   vector<double>  *Harm_FT_hit_xg;
   vector<double>  *Harm_FT_hit_yg;
   vector<double>  *Harm_FT_hit_zg;
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
   vector<double>  *Harm_FT_Track_Sx;
   vector<double>  *Harm_FT_Track_Sy;
   vector<double>  *Harm_FT_Track_Sz;
   vector<double>  *Harm_FT_Track_Xfit;
   vector<double>  *Harm_FT_Track_Yfit;
   vector<double>  *Harm_FT_Track_Xpfit;
   vector<double>  *Harm_FT_Track_Ypfit;
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
   Double_t        Harm_HCalScint_det_esum;
   Int_t           Harm_HCalScint_hit_nhits;
   vector<int>     *Harm_HCalScint_hit_row;
   vector<int>     *Harm_HCalScint_hit_col;
   vector<int>     *Harm_HCalScint_hit_cell;
   vector<int>     *Harm_HCalScint_hit_plane;
   vector<int>     *Harm_HCalScint_hit_wire;
   vector<double>  *Harm_HCalScint_hit_xcell;
   vector<double>  *Harm_HCalScint_hit_ycell;
   vector<double>  *Harm_HCalScint_hit_zcell;
   vector<double>  *Harm_HCalScint_hit_xcellg;
   vector<double>  *Harm_HCalScint_hit_ycellg;
   vector<double>  *Harm_HCalScint_hit_zcellg;
   vector<double>  *Harm_HCalScint_hit_xhit;
   vector<double>  *Harm_HCalScint_hit_yhit;
   vector<double>  *Harm_HCalScint_hit_zhit;
   vector<double>  *Harm_HCalScint_hit_xhitg;
   vector<double>  *Harm_HCalScint_hit_yhitg;
   vector<double>  *Harm_HCalScint_hit_zhitg;
   vector<double>  *Harm_HCalScint_hit_sumedep;
   vector<double>  *Harm_HCalScint_hit_tavg;
   vector<double>  *Harm_HCalScint_hit_trms;
   vector<double>  *Harm_HCalScint_hit_tmin;
   vector<double>  *Harm_HCalScint_hit_tmax;
   Double_t        primaries_Sigma;
   Double_t        primaries_Ebeam;
   Double_t        primaries_Eprime;
   Double_t        primaries_Q2;
   Double_t        primaries_xbj;
   Double_t        primaries_y;
   Double_t        primaries_W2;
   Double_t        primaries_theta_e;
   Double_t        primaries_phi_e;
   Double_t        primaries_px_e;
   Double_t        primaries_py_e;
   Double_t        primaries_pz_e;
   Double_t        primaries_vx_e;
   Double_t        primaries_vy_e;
   Double_t        primaries_vz_e;
   Double_t        primaries_Egamma;
   Double_t        primaries_theta_gamma;
   Double_t        primaries_phi_gamma;
   Double_t        primaries_px_gamma;
   Double_t        primaries_py_gamma;
   Double_t        primaries_pz_gamma;
   Double_t        primaries_vx_gamma;
   Double_t        primaries_vy_gamma;
   Double_t        primaries_vz_gamma;
   Int_t           Primaries_Nprimaries;
   vector<int>     *Primaries_PID;
   vector<int>     *Primaries_genflag;
   vector<double>  *Primaries_Px;
   vector<double>  *Primaries_Py;
   vector<double>  *Primaries_Pz;
   vector<double>  *Primaries_vx;
   vector<double>  *Primaries_vy;
   vector<double>  *Primaries_vz;
   vector<double>  *Primaries_M;
   vector<double>  *Primaries_E;
   vector<double>  *Primaries_P;
   vector<double>  *Primaries_t;
   vector<double>  *Primaries_theta;
   vector<double>  *Primaries_phi;

   // List of branches
   TBranch        *b_ev;   //!
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
   TBranch        *b_Earm_CDET_Scint_det_esum;   //!
   TBranch        *b_Earm_CDET_Scint_hit_nhits;   //!
   TBranch        *b_Earm_CDET_Scint_hit_row;   //!
   TBranch        *b_Earm_CDET_Scint_hit_col;   //!
   TBranch        *b_Earm_CDET_Scint_hit_cell;   //!
   TBranch        *b_Earm_CDET_Scint_hit_plane;   //!
   TBranch        *b_Earm_CDET_Scint_hit_wire;   //!
   TBranch        *b_Earm_CDET_Scint_hit_xcell;   //!
   TBranch        *b_Earm_CDET_Scint_hit_ycell;   //!
   TBranch        *b_Earm_CDET_Scint_hit_zcell;   //!
   TBranch        *b_Earm_CDET_Scint_hit_xcellg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_ycellg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_zcellg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_xhit;   //!
   TBranch        *b_Earm_CDET_Scint_hit_yhit;   //!
   TBranch        *b_Earm_CDET_Scint_hit_zhit;   //!
   TBranch        *b_Earm_CDET_Scint_hit_xhitg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_yhitg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_zhitg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_sumedep;   //!
   TBranch        *b_Earm_CDET_Scint_hit_tavg;   //!
   TBranch        *b_Earm_CDET_Scint_hit_trms;   //!
   TBranch        *b_Earm_CDET_Scint_hit_tmin;   //!
   TBranch        *b_Earm_CDET_Scint_hit_tmax;   //!
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
   TBranch        *b_Earm_ECalTF1_det_esum;   //!
   TBranch        *b_Earm_ECalTF1_hit_nhits;   //!
   TBranch        *b_Earm_ECalTF1_hit_row;   //!
   TBranch        *b_Earm_ECalTF1_hit_col;   //!
   TBranch        *b_Earm_ECalTF1_hit_cell;   //!
   TBranch        *b_Earm_ECalTF1_hit_plane;   //!
   TBranch        *b_Earm_ECalTF1_hit_wire;   //!
   TBranch        *b_Earm_ECalTF1_hit_xcell;   //!
   TBranch        *b_Earm_ECalTF1_hit_ycell;   //!
   TBranch        *b_Earm_ECalTF1_hit_zcell;   //!
   TBranch        *b_Earm_ECalTF1_hit_xcellg;   //!
   TBranch        *b_Earm_ECalTF1_hit_ycellg;   //!
   TBranch        *b_Earm_ECalTF1_hit_zcellg;   //!
   TBranch        *b_Earm_ECalTF1_hit_xhit;   //!
   TBranch        *b_Earm_ECalTF1_hit_yhit;   //!
   TBranch        *b_Earm_ECalTF1_hit_zhit;   //!
   TBranch        *b_Earm_ECalTF1_hit_xhitg;   //!
   TBranch        *b_Earm_ECalTF1_hit_yhitg;   //!
   TBranch        *b_Earm_ECalTF1_hit_zhitg;   //!
   TBranch        *b_Earm_ECalTF1_hit_sumedep;   //!
   TBranch        *b_Earm_ECalTF1_hit_tavg;   //!
   TBranch        *b_Earm_ECalTF1_hit_trms;   //!
   TBranch        *b_Earm_ECalTF1_hit_tmin;   //!
   TBranch        *b_Earm_ECalTF1_hit_tmax;   //!
   TBranch        *b_Harm_FPP1_hit_nhits;   //!
   TBranch        *b_Harm_FPP1_hit_plane;   //!
   TBranch        *b_Harm_FPP1_hit_strip;   //!
   TBranch        *b_Harm_FPP1_hit_x;   //!
   TBranch        *b_Harm_FPP1_hit_y;   //!
   TBranch        *b_Harm_FPP1_hit_z;   //!
   TBranch        *b_Harm_FPP1_hit_polx;   //!
   TBranch        *b_Harm_FPP1_hit_poly;   //!
   TBranch        *b_Harm_FPP1_hit_polz;   //!
   TBranch        *b_Harm_FPP1_hit_t;   //!
   TBranch        *b_Harm_FPP1_hit_trms;   //!
   TBranch        *b_Harm_FPP1_hit_tmin;   //!
   TBranch        *b_Harm_FPP1_hit_tmax;   //!
   TBranch        *b_Harm_FPP1_hit_tx;   //!
   TBranch        *b_Harm_FPP1_hit_ty;   //!
   TBranch        *b_Harm_FPP1_hit_xin;   //!
   TBranch        *b_Harm_FPP1_hit_yin;   //!
   TBranch        *b_Harm_FPP1_hit_zin;   //!
   TBranch        *b_Harm_FPP1_hit_xout;   //!
   TBranch        *b_Harm_FPP1_hit_yout;   //!
   TBranch        *b_Harm_FPP1_hit_zout;   //!
   TBranch        *b_Harm_FPP1_hit_txp;   //!
   TBranch        *b_Harm_FPP1_hit_typ;   //!
   TBranch        *b_Harm_FPP1_hit_xg;   //!
   TBranch        *b_Harm_FPP1_hit_yg;   //!
   TBranch        *b_Harm_FPP1_hit_zg;   //!
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
   TBranch        *b_Harm_FPP1_Track_Sx;   //!
   TBranch        *b_Harm_FPP1_Track_Sy;   //!
   TBranch        *b_Harm_FPP1_Track_Sz;   //!
   TBranch        *b_Harm_FPP1_Track_Xfit;   //!
   TBranch        *b_Harm_FPP1_Track_Yfit;   //!
   TBranch        *b_Harm_FPP1_Track_Xpfit;   //!
   TBranch        *b_Harm_FPP1_Track_Ypfit;   //!
   TBranch        *b_Harm_FPP2_hit_nhits;   //!
   TBranch        *b_Harm_FPP2_hit_plane;   //!
   TBranch        *b_Harm_FPP2_hit_strip;   //!
   TBranch        *b_Harm_FPP2_hit_x;   //!
   TBranch        *b_Harm_FPP2_hit_y;   //!
   TBranch        *b_Harm_FPP2_hit_z;   //!
   TBranch        *b_Harm_FPP2_hit_polx;   //!
   TBranch        *b_Harm_FPP2_hit_poly;   //!
   TBranch        *b_Harm_FPP2_hit_polz;   //!
   TBranch        *b_Harm_FPP2_hit_t;   //!
   TBranch        *b_Harm_FPP2_hit_trms;   //!
   TBranch        *b_Harm_FPP2_hit_tmin;   //!
   TBranch        *b_Harm_FPP2_hit_tmax;   //!
   TBranch        *b_Harm_FPP2_hit_tx;   //!
   TBranch        *b_Harm_FPP2_hit_ty;   //!
   TBranch        *b_Harm_FPP2_hit_xin;   //!
   TBranch        *b_Harm_FPP2_hit_yin;   //!
   TBranch        *b_Harm_FPP2_hit_zin;   //!
   TBranch        *b_Harm_FPP2_hit_xout;   //!
   TBranch        *b_Harm_FPP2_hit_yout;   //!
   TBranch        *b_Harm_FPP2_hit_zout;   //!
   TBranch        *b_Harm_FPP2_hit_txp;   //!
   TBranch        *b_Harm_FPP2_hit_typ;   //!
   TBranch        *b_Harm_FPP2_hit_xg;   //!
   TBranch        *b_Harm_FPP2_hit_yg;   //!
   TBranch        *b_Harm_FPP2_hit_zg;   //!
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
   TBranch        *b_Harm_FPP2_Track_Sx;   //!
   TBranch        *b_Harm_FPP2_Track_Sy;   //!
   TBranch        *b_Harm_FPP2_Track_Sz;   //!
   TBranch        *b_Harm_FPP2_Track_Xfit;   //!
   TBranch        *b_Harm_FPP2_Track_Yfit;   //!
   TBranch        *b_Harm_FPP2_Track_Xpfit;   //!
   TBranch        *b_Harm_FPP2_Track_Ypfit;   //!
   TBranch        *b_Harm_FT_hit_nhits;   //!
   TBranch        *b_Harm_FT_hit_plane;   //!
   TBranch        *b_Harm_FT_hit_strip;   //!
   TBranch        *b_Harm_FT_hit_x;   //!
   TBranch        *b_Harm_FT_hit_y;   //!
   TBranch        *b_Harm_FT_hit_z;   //!
   TBranch        *b_Harm_FT_hit_polx;   //!
   TBranch        *b_Harm_FT_hit_poly;   //!
   TBranch        *b_Harm_FT_hit_polz;   //!
   TBranch        *b_Harm_FT_hit_t;   //!
   TBranch        *b_Harm_FT_hit_trms;   //!
   TBranch        *b_Harm_FT_hit_tmin;   //!
   TBranch        *b_Harm_FT_hit_tmax;   //!
   TBranch        *b_Harm_FT_hit_tx;   //!
   TBranch        *b_Harm_FT_hit_ty;   //!
   TBranch        *b_Harm_FT_hit_xin;   //!
   TBranch        *b_Harm_FT_hit_yin;   //!
   TBranch        *b_Harm_FT_hit_zin;   //!
   TBranch        *b_Harm_FT_hit_xout;   //!
   TBranch        *b_Harm_FT_hit_yout;   //!
   TBranch        *b_Harm_FT_hit_zout;   //!
   TBranch        *b_Harm_FT_hit_txp;   //!
   TBranch        *b_Harm_FT_hit_typ;   //!
   TBranch        *b_Harm_FT_hit_xg;   //!
   TBranch        *b_Harm_FT_hit_yg;   //!
   TBranch        *b_Harm_FT_hit_zg;   //!
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
   TBranch        *b_Harm_FT_Track_Sx;   //!
   TBranch        *b_Harm_FT_Track_Sy;   //!
   TBranch        *b_Harm_FT_Track_Sz;   //!
   TBranch        *b_Harm_FT_Track_Xfit;   //!
   TBranch        *b_Harm_FT_Track_Yfit;   //!
   TBranch        *b_Harm_FT_Track_Xpfit;   //!
   TBranch        *b_Harm_FT_Track_Ypfit;   //!
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
   TBranch        *b_Harm_HCalScint_det_esum;   //!
   TBranch        *b_Harm_HCalScint_hit_nhits;   //!
   TBranch        *b_Harm_HCalScint_hit_row;   //!
   TBranch        *b_Harm_HCalScint_hit_col;   //!
   TBranch        *b_Harm_HCalScint_hit_cell;   //!
   TBranch        *b_Harm_HCalScint_hit_plane;   //!
   TBranch        *b_Harm_HCalScint_hit_wire;   //!
   TBranch        *b_Harm_HCalScint_hit_xcell;   //!
   TBranch        *b_Harm_HCalScint_hit_ycell;   //!
   TBranch        *b_Harm_HCalScint_hit_zcell;   //!
   TBranch        *b_Harm_HCalScint_hit_xcellg;   //!
   TBranch        *b_Harm_HCalScint_hit_ycellg;   //!
   TBranch        *b_Harm_HCalScint_hit_zcellg;   //!
   TBranch        *b_Harm_HCalScint_hit_xhit;   //!
   TBranch        *b_Harm_HCalScint_hit_yhit;   //!
   TBranch        *b_Harm_HCalScint_hit_zhit;   //!
   TBranch        *b_Harm_HCalScint_hit_xhitg;   //!
   TBranch        *b_Harm_HCalScint_hit_yhitg;   //!
   TBranch        *b_Harm_HCalScint_hit_zhitg;   //!
   TBranch        *b_Harm_HCalScint_hit_sumedep;   //!
   TBranch        *b_Harm_HCalScint_hit_tavg;   //!
   TBranch        *b_Harm_HCalScint_hit_trms;   //!
   TBranch        *b_Harm_HCalScint_hit_tmin;   //!
   TBranch        *b_Harm_HCalScint_hit_tmax;   //!
   TBranch        *b_primaries_Sigma;   //!
   TBranch        *b_primaries_Ebeam;   //!
   TBranch        *b_primaries_Eprime;   //!
   TBranch        *b_primaries_Q2;   //!
   TBranch        *b_primaries_xbj;   //!
   TBranch        *b_primaries_y;   //!
   TBranch        *b_primaries_W2;   //!
   TBranch        *b_primaries_theta_e;   //!
   TBranch        *b_primaries_phi_e;   //!
   TBranch        *b_primaries_px_e;   //!
   TBranch        *b_primaries_py_e;   //!
   TBranch        *b_primaries_pz_e;   //!
   TBranch        *b_primaries_vx_e;   //!
   TBranch        *b_primaries_vy_e;   //!
   TBranch        *b_primaries_vz_e;   //!
   TBranch        *b_primaries_Egamma;   //!
   TBranch        *b_primaries_theta_gamma;   //!
   TBranch        *b_primaries_phi_gamma;   //!
   TBranch        *b_primaries_px_gamma;   //!
   TBranch        *b_primaries_py_gamma;   //!
   TBranch        *b_primaries_pz_gamma;   //!
   TBranch        *b_primaries_vx_gamma;   //!
   TBranch        *b_primaries_vy_gamma;   //!
   TBranch        *b_primaries_vz_gamma;   //!
   TBranch        *b_Nprimaries;   //!
   TBranch        *b_Primaries_PID;   //!
   TBranch        *b_Primaries_genflag;   //!
   TBranch        *b_Primaries_Px;   //!
   TBranch        *b_Primaries_Py;   //!
   TBranch        *b_Primaries_Pz;   //!
   TBranch        *b_Primaries_vx;   //!
   TBranch        *b_Primaries_vy;   //!
   TBranch        *b_Primaries_vz;   //!
   TBranch        *b_Primaries_M;   //!
   TBranch        *b_Primaries_E;   //!
   TBranch        *b_Primaries_P;   //!
   TBranch        *b_Primaries_t;   //!
   TBranch        *b_Primaries_theta;   //!
   TBranch        *b_Primaries_phi;   //!

   gep_pythia6_tree(TTree *tree=0);
   virtual ~gep_pythia6_tree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gep_pythia6_tree_cxx
gep_pythia6_tree::gep_pythia6_tree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("gep_12GeV2_pythia6_opticalphotons_1kevents.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("gep_12GeV2_pythia6_opticalphotons_1kevents.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

gep_pythia6_tree::~gep_pythia6_tree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gep_pythia6_tree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gep_pythia6_tree::LoadTree(Long64_t entry)
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

void gep_pythia6_tree::Init(TTree *tree)
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
   Earm_CDET_Scint_hit_cell = 0;
   Earm_CDET_Scint_hit_plane = 0;
   Earm_CDET_Scint_hit_wire = 0;
   Earm_CDET_Scint_hit_xcell = 0;
   Earm_CDET_Scint_hit_ycell = 0;
   Earm_CDET_Scint_hit_zcell = 0;
   Earm_CDET_Scint_hit_xcellg = 0;
   Earm_CDET_Scint_hit_ycellg = 0;
   Earm_CDET_Scint_hit_zcellg = 0;
   Earm_CDET_Scint_hit_xhit = 0;
   Earm_CDET_Scint_hit_yhit = 0;
   Earm_CDET_Scint_hit_zhit = 0;
   Earm_CDET_Scint_hit_xhitg = 0;
   Earm_CDET_Scint_hit_yhitg = 0;
   Earm_CDET_Scint_hit_zhitg = 0;
   Earm_CDET_Scint_hit_sumedep = 0;
   Earm_CDET_Scint_hit_tavg = 0;
   Earm_CDET_Scint_hit_trms = 0;
   Earm_CDET_Scint_hit_tmin = 0;
   Earm_CDET_Scint_hit_tmax = 0;
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
   Earm_ECalTF1_hit_cell = 0;
   Earm_ECalTF1_hit_plane = 0;
   Earm_ECalTF1_hit_wire = 0;
   Earm_ECalTF1_hit_xcell = 0;
   Earm_ECalTF1_hit_ycell = 0;
   Earm_ECalTF1_hit_zcell = 0;
   Earm_ECalTF1_hit_xcellg = 0;
   Earm_ECalTF1_hit_ycellg = 0;
   Earm_ECalTF1_hit_zcellg = 0;
   Earm_ECalTF1_hit_xhit = 0;
   Earm_ECalTF1_hit_yhit = 0;
   Earm_ECalTF1_hit_zhit = 0;
   Earm_ECalTF1_hit_xhitg = 0;
   Earm_ECalTF1_hit_yhitg = 0;
   Earm_ECalTF1_hit_zhitg = 0;
   Earm_ECalTF1_hit_sumedep = 0;
   Earm_ECalTF1_hit_tavg = 0;
   Earm_ECalTF1_hit_trms = 0;
   Earm_ECalTF1_hit_tmin = 0;
   Earm_ECalTF1_hit_tmax = 0;
   Harm_FPP1_hit_plane = 0;
   Harm_FPP1_hit_strip = 0;
   Harm_FPP1_hit_x = 0;
   Harm_FPP1_hit_y = 0;
   Harm_FPP1_hit_z = 0;
   Harm_FPP1_hit_polx = 0;
   Harm_FPP1_hit_poly = 0;
   Harm_FPP1_hit_polz = 0;
   Harm_FPP1_hit_t = 0;
   Harm_FPP1_hit_trms = 0;
   Harm_FPP1_hit_tmin = 0;
   Harm_FPP1_hit_tmax = 0;
   Harm_FPP1_hit_tx = 0;
   Harm_FPP1_hit_ty = 0;
   Harm_FPP1_hit_xin = 0;
   Harm_FPP1_hit_yin = 0;
   Harm_FPP1_hit_zin = 0;
   Harm_FPP1_hit_xout = 0;
   Harm_FPP1_hit_yout = 0;
   Harm_FPP1_hit_zout = 0;
   Harm_FPP1_hit_txp = 0;
   Harm_FPP1_hit_typ = 0;
   Harm_FPP1_hit_xg = 0;
   Harm_FPP1_hit_yg = 0;
   Harm_FPP1_hit_zg = 0;
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
   Harm_FPP1_Track_Sx = 0;
   Harm_FPP1_Track_Sy = 0;
   Harm_FPP1_Track_Sz = 0;
   Harm_FPP1_Track_Xfit = 0;
   Harm_FPP1_Track_Yfit = 0;
   Harm_FPP1_Track_Xpfit = 0;
   Harm_FPP1_Track_Ypfit = 0;
   Harm_FPP2_hit_plane = 0;
   Harm_FPP2_hit_strip = 0;
   Harm_FPP2_hit_x = 0;
   Harm_FPP2_hit_y = 0;
   Harm_FPP2_hit_z = 0;
   Harm_FPP2_hit_polx = 0;
   Harm_FPP2_hit_poly = 0;
   Harm_FPP2_hit_polz = 0;
   Harm_FPP2_hit_t = 0;
   Harm_FPP2_hit_trms = 0;
   Harm_FPP2_hit_tmin = 0;
   Harm_FPP2_hit_tmax = 0;
   Harm_FPP2_hit_tx = 0;
   Harm_FPP2_hit_ty = 0;
   Harm_FPP2_hit_xin = 0;
   Harm_FPP2_hit_yin = 0;
   Harm_FPP2_hit_zin = 0;
   Harm_FPP2_hit_xout = 0;
   Harm_FPP2_hit_yout = 0;
   Harm_FPP2_hit_zout = 0;
   Harm_FPP2_hit_txp = 0;
   Harm_FPP2_hit_typ = 0;
   Harm_FPP2_hit_xg = 0;
   Harm_FPP2_hit_yg = 0;
   Harm_FPP2_hit_zg = 0;
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
   Harm_FPP2_Track_Sx = 0;
   Harm_FPP2_Track_Sy = 0;
   Harm_FPP2_Track_Sz = 0;
   Harm_FPP2_Track_Xfit = 0;
   Harm_FPP2_Track_Yfit = 0;
   Harm_FPP2_Track_Xpfit = 0;
   Harm_FPP2_Track_Ypfit = 0;
   Harm_FT_hit_plane = 0;
   Harm_FT_hit_strip = 0;
   Harm_FT_hit_x = 0;
   Harm_FT_hit_y = 0;
   Harm_FT_hit_z = 0;
   Harm_FT_hit_polx = 0;
   Harm_FT_hit_poly = 0;
   Harm_FT_hit_polz = 0;
   Harm_FT_hit_t = 0;
   Harm_FT_hit_trms = 0;
   Harm_FT_hit_tmin = 0;
   Harm_FT_hit_tmax = 0;
   Harm_FT_hit_tx = 0;
   Harm_FT_hit_ty = 0;
   Harm_FT_hit_xin = 0;
   Harm_FT_hit_yin = 0;
   Harm_FT_hit_zin = 0;
   Harm_FT_hit_xout = 0;
   Harm_FT_hit_yout = 0;
   Harm_FT_hit_zout = 0;
   Harm_FT_hit_txp = 0;
   Harm_FT_hit_typ = 0;
   Harm_FT_hit_xg = 0;
   Harm_FT_hit_yg = 0;
   Harm_FT_hit_zg = 0;
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
   Harm_FT_Track_Sx = 0;
   Harm_FT_Track_Sy = 0;
   Harm_FT_Track_Sz = 0;
   Harm_FT_Track_Xfit = 0;
   Harm_FT_Track_Yfit = 0;
   Harm_FT_Track_Xpfit = 0;
   Harm_FT_Track_Ypfit = 0;
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
   Harm_HCalScint_hit_wire = 0;
   Harm_HCalScint_hit_xcell = 0;
   Harm_HCalScint_hit_ycell = 0;
   Harm_HCalScint_hit_zcell = 0;
   Harm_HCalScint_hit_xcellg = 0;
   Harm_HCalScint_hit_ycellg = 0;
   Harm_HCalScint_hit_zcellg = 0;
   Harm_HCalScint_hit_xhit = 0;
   Harm_HCalScint_hit_yhit = 0;
   Harm_HCalScint_hit_zhit = 0;
   Harm_HCalScint_hit_xhitg = 0;
   Harm_HCalScint_hit_yhitg = 0;
   Harm_HCalScint_hit_zhitg = 0;
   Harm_HCalScint_hit_sumedep = 0;
   Harm_HCalScint_hit_tavg = 0;
   Harm_HCalScint_hit_trms = 0;
   Harm_HCalScint_hit_tmin = 0;
   Harm_HCalScint_hit_tmax = 0;
   Primaries_PID = 0;
   Primaries_genflag = 0;
   Primaries_Px = 0;
   Primaries_Py = 0;
   Primaries_Pz = 0;
   Primaries_vx = 0;
   Primaries_vy = 0;
   Primaries_vz = 0;
   Primaries_M = 0;
   Primaries_E = 0;
   Primaries_P = 0;
   Primaries_t = 0;
   Primaries_theta = 0;
   Primaries_phi = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev", &ev_count, &b_ev);
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
   fChain->SetBranchAddress("Earm.CDET_Scint.det.esum", &Earm_CDET_Scint_det_esum, &b_Earm_CDET_Scint_det_esum);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.nhits", &Earm_CDET_Scint_hit_nhits, &b_Earm_CDET_Scint_hit_nhits);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.row", &Earm_CDET_Scint_hit_row, &b_Earm_CDET_Scint_hit_row);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.col", &Earm_CDET_Scint_hit_col, &b_Earm_CDET_Scint_hit_col);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.cell", &Earm_CDET_Scint_hit_cell, &b_Earm_CDET_Scint_hit_cell);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.plane", &Earm_CDET_Scint_hit_plane, &b_Earm_CDET_Scint_hit_plane);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.wire", &Earm_CDET_Scint_hit_wire, &b_Earm_CDET_Scint_hit_wire);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.xcell", &Earm_CDET_Scint_hit_xcell, &b_Earm_CDET_Scint_hit_xcell);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.ycell", &Earm_CDET_Scint_hit_ycell, &b_Earm_CDET_Scint_hit_ycell);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.zcell", &Earm_CDET_Scint_hit_zcell, &b_Earm_CDET_Scint_hit_zcell);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.xcellg", &Earm_CDET_Scint_hit_xcellg, &b_Earm_CDET_Scint_hit_xcellg);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.ycellg", &Earm_CDET_Scint_hit_ycellg, &b_Earm_CDET_Scint_hit_ycellg);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.zcellg", &Earm_CDET_Scint_hit_zcellg, &b_Earm_CDET_Scint_hit_zcellg);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.xhit", &Earm_CDET_Scint_hit_xhit, &b_Earm_CDET_Scint_hit_xhit);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.yhit", &Earm_CDET_Scint_hit_yhit, &b_Earm_CDET_Scint_hit_yhit);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.zhit", &Earm_CDET_Scint_hit_zhit, &b_Earm_CDET_Scint_hit_zhit);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.xhitg", &Earm_CDET_Scint_hit_xhitg, &b_Earm_CDET_Scint_hit_xhitg);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.yhitg", &Earm_CDET_Scint_hit_yhitg, &b_Earm_CDET_Scint_hit_yhitg);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.zhitg", &Earm_CDET_Scint_hit_zhitg, &b_Earm_CDET_Scint_hit_zhitg);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.sumedep", &Earm_CDET_Scint_hit_sumedep, &b_Earm_CDET_Scint_hit_sumedep);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.tavg", &Earm_CDET_Scint_hit_tavg, &b_Earm_CDET_Scint_hit_tavg);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.trms", &Earm_CDET_Scint_hit_trms, &b_Earm_CDET_Scint_hit_trms);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.tmin", &Earm_CDET_Scint_hit_tmin, &b_Earm_CDET_Scint_hit_tmin);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.tmax", &Earm_CDET_Scint_hit_tmax, &b_Earm_CDET_Scint_hit_tmax);
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
   fChain->SetBranchAddress("Earm.ECalTF1.det.esum", &Earm_ECalTF1_det_esum, &b_Earm_ECalTF1_det_esum);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.nhits", &Earm_ECalTF1_hit_nhits, &b_Earm_ECalTF1_hit_nhits);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.row", &Earm_ECalTF1_hit_row, &b_Earm_ECalTF1_hit_row);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.col", &Earm_ECalTF1_hit_col, &b_Earm_ECalTF1_hit_col);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.cell", &Earm_ECalTF1_hit_cell, &b_Earm_ECalTF1_hit_cell);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.plane", &Earm_ECalTF1_hit_plane, &b_Earm_ECalTF1_hit_plane);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.wire", &Earm_ECalTF1_hit_wire, &b_Earm_ECalTF1_hit_wire);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.xcell", &Earm_ECalTF1_hit_xcell, &b_Earm_ECalTF1_hit_xcell);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.ycell", &Earm_ECalTF1_hit_ycell, &b_Earm_ECalTF1_hit_ycell);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.zcell", &Earm_ECalTF1_hit_zcell, &b_Earm_ECalTF1_hit_zcell);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.xcellg", &Earm_ECalTF1_hit_xcellg, &b_Earm_ECalTF1_hit_xcellg);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.ycellg", &Earm_ECalTF1_hit_ycellg, &b_Earm_ECalTF1_hit_ycellg);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.zcellg", &Earm_ECalTF1_hit_zcellg, &b_Earm_ECalTF1_hit_zcellg);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.xhit", &Earm_ECalTF1_hit_xhit, &b_Earm_ECalTF1_hit_xhit);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.yhit", &Earm_ECalTF1_hit_yhit, &b_Earm_ECalTF1_hit_yhit);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.zhit", &Earm_ECalTF1_hit_zhit, &b_Earm_ECalTF1_hit_zhit);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.xhitg", &Earm_ECalTF1_hit_xhitg, &b_Earm_ECalTF1_hit_xhitg);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.yhitg", &Earm_ECalTF1_hit_yhitg, &b_Earm_ECalTF1_hit_yhitg);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.zhitg", &Earm_ECalTF1_hit_zhitg, &b_Earm_ECalTF1_hit_zhitg);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.sumedep", &Earm_ECalTF1_hit_sumedep, &b_Earm_ECalTF1_hit_sumedep);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.tavg", &Earm_ECalTF1_hit_tavg, &b_Earm_ECalTF1_hit_tavg);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.trms", &Earm_ECalTF1_hit_trms, &b_Earm_ECalTF1_hit_trms);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.tmin", &Earm_ECalTF1_hit_tmin, &b_Earm_ECalTF1_hit_tmin);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.tmax", &Earm_ECalTF1_hit_tmax, &b_Earm_ECalTF1_hit_tmax);
   fChain->SetBranchAddress("Harm.FPP1.hit.nhits", &Harm_FPP1_hit_nhits, &b_Harm_FPP1_hit_nhits);
   fChain->SetBranchAddress("Harm.FPP1.hit.plane", &Harm_FPP1_hit_plane, &b_Harm_FPP1_hit_plane);
   fChain->SetBranchAddress("Harm.FPP1.hit.strip", &Harm_FPP1_hit_strip, &b_Harm_FPP1_hit_strip);
   fChain->SetBranchAddress("Harm.FPP1.hit.x", &Harm_FPP1_hit_x, &b_Harm_FPP1_hit_x);
   fChain->SetBranchAddress("Harm.FPP1.hit.y", &Harm_FPP1_hit_y, &b_Harm_FPP1_hit_y);
   fChain->SetBranchAddress("Harm.FPP1.hit.z", &Harm_FPP1_hit_z, &b_Harm_FPP1_hit_z);
   fChain->SetBranchAddress("Harm.FPP1.hit.polx", &Harm_FPP1_hit_polx, &b_Harm_FPP1_hit_polx);
   fChain->SetBranchAddress("Harm.FPP1.hit.poly", &Harm_FPP1_hit_poly, &b_Harm_FPP1_hit_poly);
   fChain->SetBranchAddress("Harm.FPP1.hit.polz", &Harm_FPP1_hit_polz, &b_Harm_FPP1_hit_polz);
   fChain->SetBranchAddress("Harm.FPP1.hit.t", &Harm_FPP1_hit_t, &b_Harm_FPP1_hit_t);
   fChain->SetBranchAddress("Harm.FPP1.hit.trms", &Harm_FPP1_hit_trms, &b_Harm_FPP1_hit_trms);
   fChain->SetBranchAddress("Harm.FPP1.hit.tmin", &Harm_FPP1_hit_tmin, &b_Harm_FPP1_hit_tmin);
   fChain->SetBranchAddress("Harm.FPP1.hit.tmax", &Harm_FPP1_hit_tmax, &b_Harm_FPP1_hit_tmax);
   fChain->SetBranchAddress("Harm.FPP1.hit.tx", &Harm_FPP1_hit_tx, &b_Harm_FPP1_hit_tx);
   fChain->SetBranchAddress("Harm.FPP1.hit.ty", &Harm_FPP1_hit_ty, &b_Harm_FPP1_hit_ty);
   fChain->SetBranchAddress("Harm.FPP1.hit.xin", &Harm_FPP1_hit_xin, &b_Harm_FPP1_hit_xin);
   fChain->SetBranchAddress("Harm.FPP1.hit.yin", &Harm_FPP1_hit_yin, &b_Harm_FPP1_hit_yin);
   fChain->SetBranchAddress("Harm.FPP1.hit.zin", &Harm_FPP1_hit_zin, &b_Harm_FPP1_hit_zin);
   fChain->SetBranchAddress("Harm.FPP1.hit.xout", &Harm_FPP1_hit_xout, &b_Harm_FPP1_hit_xout);
   fChain->SetBranchAddress("Harm.FPP1.hit.yout", &Harm_FPP1_hit_yout, &b_Harm_FPP1_hit_yout);
   fChain->SetBranchAddress("Harm.FPP1.hit.zout", &Harm_FPP1_hit_zout, &b_Harm_FPP1_hit_zout);
   fChain->SetBranchAddress("Harm.FPP1.hit.txp", &Harm_FPP1_hit_txp, &b_Harm_FPP1_hit_txp);
   fChain->SetBranchAddress("Harm.FPP1.hit.typ", &Harm_FPP1_hit_typ, &b_Harm_FPP1_hit_typ);
   fChain->SetBranchAddress("Harm.FPP1.hit.xg", &Harm_FPP1_hit_xg, &b_Harm_FPP1_hit_xg);
   fChain->SetBranchAddress("Harm.FPP1.hit.yg", &Harm_FPP1_hit_yg, &b_Harm_FPP1_hit_yg);
   fChain->SetBranchAddress("Harm.FPP1.hit.zg", &Harm_FPP1_hit_zg, &b_Harm_FPP1_hit_zg);
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
   fChain->SetBranchAddress("Harm.FPP1.Track.Sx", &Harm_FPP1_Track_Sx, &b_Harm_FPP1_Track_Sx);
   fChain->SetBranchAddress("Harm.FPP1.Track.Sy", &Harm_FPP1_Track_Sy, &b_Harm_FPP1_Track_Sy);
   fChain->SetBranchAddress("Harm.FPP1.Track.Sz", &Harm_FPP1_Track_Sz, &b_Harm_FPP1_Track_Sz);
   fChain->SetBranchAddress("Harm.FPP1.Track.Xfit", &Harm_FPP1_Track_Xfit, &b_Harm_FPP1_Track_Xfit);
   fChain->SetBranchAddress("Harm.FPP1.Track.Yfit", &Harm_FPP1_Track_Yfit, &b_Harm_FPP1_Track_Yfit);
   fChain->SetBranchAddress("Harm.FPP1.Track.Xpfit", &Harm_FPP1_Track_Xpfit, &b_Harm_FPP1_Track_Xpfit);
   fChain->SetBranchAddress("Harm.FPP1.Track.Ypfit", &Harm_FPP1_Track_Ypfit, &b_Harm_FPP1_Track_Ypfit);
   fChain->SetBranchAddress("Harm.FPP2.hit.nhits", &Harm_FPP2_hit_nhits, &b_Harm_FPP2_hit_nhits);
   fChain->SetBranchAddress("Harm.FPP2.hit.plane", &Harm_FPP2_hit_plane, &b_Harm_FPP2_hit_plane);
   fChain->SetBranchAddress("Harm.FPP2.hit.strip", &Harm_FPP2_hit_strip, &b_Harm_FPP2_hit_strip);
   fChain->SetBranchAddress("Harm.FPP2.hit.x", &Harm_FPP2_hit_x, &b_Harm_FPP2_hit_x);
   fChain->SetBranchAddress("Harm.FPP2.hit.y", &Harm_FPP2_hit_y, &b_Harm_FPP2_hit_y);
   fChain->SetBranchAddress("Harm.FPP2.hit.z", &Harm_FPP2_hit_z, &b_Harm_FPP2_hit_z);
   fChain->SetBranchAddress("Harm.FPP2.hit.polx", &Harm_FPP2_hit_polx, &b_Harm_FPP2_hit_polx);
   fChain->SetBranchAddress("Harm.FPP2.hit.poly", &Harm_FPP2_hit_poly, &b_Harm_FPP2_hit_poly);
   fChain->SetBranchAddress("Harm.FPP2.hit.polz", &Harm_FPP2_hit_polz, &b_Harm_FPP2_hit_polz);
   fChain->SetBranchAddress("Harm.FPP2.hit.t", &Harm_FPP2_hit_t, &b_Harm_FPP2_hit_t);
   fChain->SetBranchAddress("Harm.FPP2.hit.trms", &Harm_FPP2_hit_trms, &b_Harm_FPP2_hit_trms);
   fChain->SetBranchAddress("Harm.FPP2.hit.tmin", &Harm_FPP2_hit_tmin, &b_Harm_FPP2_hit_tmin);
   fChain->SetBranchAddress("Harm.FPP2.hit.tmax", &Harm_FPP2_hit_tmax, &b_Harm_FPP2_hit_tmax);
   fChain->SetBranchAddress("Harm.FPP2.hit.tx", &Harm_FPP2_hit_tx, &b_Harm_FPP2_hit_tx);
   fChain->SetBranchAddress("Harm.FPP2.hit.ty", &Harm_FPP2_hit_ty, &b_Harm_FPP2_hit_ty);
   fChain->SetBranchAddress("Harm.FPP2.hit.xin", &Harm_FPP2_hit_xin, &b_Harm_FPP2_hit_xin);
   fChain->SetBranchAddress("Harm.FPP2.hit.yin", &Harm_FPP2_hit_yin, &b_Harm_FPP2_hit_yin);
   fChain->SetBranchAddress("Harm.FPP2.hit.zin", &Harm_FPP2_hit_zin, &b_Harm_FPP2_hit_zin);
   fChain->SetBranchAddress("Harm.FPP2.hit.xout", &Harm_FPP2_hit_xout, &b_Harm_FPP2_hit_xout);
   fChain->SetBranchAddress("Harm.FPP2.hit.yout", &Harm_FPP2_hit_yout, &b_Harm_FPP2_hit_yout);
   fChain->SetBranchAddress("Harm.FPP2.hit.zout", &Harm_FPP2_hit_zout, &b_Harm_FPP2_hit_zout);
   fChain->SetBranchAddress("Harm.FPP2.hit.txp", &Harm_FPP2_hit_txp, &b_Harm_FPP2_hit_txp);
   fChain->SetBranchAddress("Harm.FPP2.hit.typ", &Harm_FPP2_hit_typ, &b_Harm_FPP2_hit_typ);
   fChain->SetBranchAddress("Harm.FPP2.hit.xg", &Harm_FPP2_hit_xg, &b_Harm_FPP2_hit_xg);
   fChain->SetBranchAddress("Harm.FPP2.hit.yg", &Harm_FPP2_hit_yg, &b_Harm_FPP2_hit_yg);
   fChain->SetBranchAddress("Harm.FPP2.hit.zg", &Harm_FPP2_hit_zg, &b_Harm_FPP2_hit_zg);
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
   fChain->SetBranchAddress("Harm.FPP2.Track.Sx", &Harm_FPP2_Track_Sx, &b_Harm_FPP2_Track_Sx);
   fChain->SetBranchAddress("Harm.FPP2.Track.Sy", &Harm_FPP2_Track_Sy, &b_Harm_FPP2_Track_Sy);
   fChain->SetBranchAddress("Harm.FPP2.Track.Sz", &Harm_FPP2_Track_Sz, &b_Harm_FPP2_Track_Sz);
   fChain->SetBranchAddress("Harm.FPP2.Track.Xfit", &Harm_FPP2_Track_Xfit, &b_Harm_FPP2_Track_Xfit);
   fChain->SetBranchAddress("Harm.FPP2.Track.Yfit", &Harm_FPP2_Track_Yfit, &b_Harm_FPP2_Track_Yfit);
   fChain->SetBranchAddress("Harm.FPP2.Track.Xpfit", &Harm_FPP2_Track_Xpfit, &b_Harm_FPP2_Track_Xpfit);
   fChain->SetBranchAddress("Harm.FPP2.Track.Ypfit", &Harm_FPP2_Track_Ypfit, &b_Harm_FPP2_Track_Ypfit);
   fChain->SetBranchAddress("Harm.FT.hit.nhits", &Harm_FT_hit_nhits, &b_Harm_FT_hit_nhits);
   fChain->SetBranchAddress("Harm.FT.hit.plane", &Harm_FT_hit_plane, &b_Harm_FT_hit_plane);
   fChain->SetBranchAddress("Harm.FT.hit.strip", &Harm_FT_hit_strip, &b_Harm_FT_hit_strip);
   fChain->SetBranchAddress("Harm.FT.hit.x", &Harm_FT_hit_x, &b_Harm_FT_hit_x);
   fChain->SetBranchAddress("Harm.FT.hit.y", &Harm_FT_hit_y, &b_Harm_FT_hit_y);
   fChain->SetBranchAddress("Harm.FT.hit.z", &Harm_FT_hit_z, &b_Harm_FT_hit_z);
   fChain->SetBranchAddress("Harm.FT.hit.polx", &Harm_FT_hit_polx, &b_Harm_FT_hit_polx);
   fChain->SetBranchAddress("Harm.FT.hit.poly", &Harm_FT_hit_poly, &b_Harm_FT_hit_poly);
   fChain->SetBranchAddress("Harm.FT.hit.polz", &Harm_FT_hit_polz, &b_Harm_FT_hit_polz);
   fChain->SetBranchAddress("Harm.FT.hit.t", &Harm_FT_hit_t, &b_Harm_FT_hit_t);
   fChain->SetBranchAddress("Harm.FT.hit.trms", &Harm_FT_hit_trms, &b_Harm_FT_hit_trms);
   fChain->SetBranchAddress("Harm.FT.hit.tmin", &Harm_FT_hit_tmin, &b_Harm_FT_hit_tmin);
   fChain->SetBranchAddress("Harm.FT.hit.tmax", &Harm_FT_hit_tmax, &b_Harm_FT_hit_tmax);
   fChain->SetBranchAddress("Harm.FT.hit.tx", &Harm_FT_hit_tx, &b_Harm_FT_hit_tx);
   fChain->SetBranchAddress("Harm.FT.hit.ty", &Harm_FT_hit_ty, &b_Harm_FT_hit_ty);
   fChain->SetBranchAddress("Harm.FT.hit.xin", &Harm_FT_hit_xin, &b_Harm_FT_hit_xin);
   fChain->SetBranchAddress("Harm.FT.hit.yin", &Harm_FT_hit_yin, &b_Harm_FT_hit_yin);
   fChain->SetBranchAddress("Harm.FT.hit.zin", &Harm_FT_hit_zin, &b_Harm_FT_hit_zin);
   fChain->SetBranchAddress("Harm.FT.hit.xout", &Harm_FT_hit_xout, &b_Harm_FT_hit_xout);
   fChain->SetBranchAddress("Harm.FT.hit.yout", &Harm_FT_hit_yout, &b_Harm_FT_hit_yout);
   fChain->SetBranchAddress("Harm.FT.hit.zout", &Harm_FT_hit_zout, &b_Harm_FT_hit_zout);
   fChain->SetBranchAddress("Harm.FT.hit.txp", &Harm_FT_hit_txp, &b_Harm_FT_hit_txp);
   fChain->SetBranchAddress("Harm.FT.hit.typ", &Harm_FT_hit_typ, &b_Harm_FT_hit_typ);
   fChain->SetBranchAddress("Harm.FT.hit.xg", &Harm_FT_hit_xg, &b_Harm_FT_hit_xg);
   fChain->SetBranchAddress("Harm.FT.hit.yg", &Harm_FT_hit_yg, &b_Harm_FT_hit_yg);
   fChain->SetBranchAddress("Harm.FT.hit.zg", &Harm_FT_hit_zg, &b_Harm_FT_hit_zg);
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
   fChain->SetBranchAddress("Harm.FT.Track.Sx", &Harm_FT_Track_Sx, &b_Harm_FT_Track_Sx);
   fChain->SetBranchAddress("Harm.FT.Track.Sy", &Harm_FT_Track_Sy, &b_Harm_FT_Track_Sy);
   fChain->SetBranchAddress("Harm.FT.Track.Sz", &Harm_FT_Track_Sz, &b_Harm_FT_Track_Sz);
   fChain->SetBranchAddress("Harm.FT.Track.Xfit", &Harm_FT_Track_Xfit, &b_Harm_FT_Track_Xfit);
   fChain->SetBranchAddress("Harm.FT.Track.Yfit", &Harm_FT_Track_Yfit, &b_Harm_FT_Track_Yfit);
   fChain->SetBranchAddress("Harm.FT.Track.Xpfit", &Harm_FT_Track_Xpfit, &b_Harm_FT_Track_Xpfit);
   fChain->SetBranchAddress("Harm.FT.Track.Ypfit", &Harm_FT_Track_Ypfit, &b_Harm_FT_Track_Ypfit);
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
   fChain->SetBranchAddress("Harm.HCalScint.det.esum", &Harm_HCalScint_det_esum, &b_Harm_HCalScint_det_esum);
   fChain->SetBranchAddress("Harm.HCalScint.hit.nhits", &Harm_HCalScint_hit_nhits, &b_Harm_HCalScint_hit_nhits);
   fChain->SetBranchAddress("Harm.HCalScint.hit.row", &Harm_HCalScint_hit_row, &b_Harm_HCalScint_hit_row);
   fChain->SetBranchAddress("Harm.HCalScint.hit.col", &Harm_HCalScint_hit_col, &b_Harm_HCalScint_hit_col);
   fChain->SetBranchAddress("Harm.HCalScint.hit.cell", &Harm_HCalScint_hit_cell, &b_Harm_HCalScint_hit_cell);
   fChain->SetBranchAddress("Harm.HCalScint.hit.plane", &Harm_HCalScint_hit_plane, &b_Harm_HCalScint_hit_plane);
   fChain->SetBranchAddress("Harm.HCalScint.hit.wire", &Harm_HCalScint_hit_wire, &b_Harm_HCalScint_hit_wire);
   fChain->SetBranchAddress("Harm.HCalScint.hit.xcell", &Harm_HCalScint_hit_xcell, &b_Harm_HCalScint_hit_xcell);
   fChain->SetBranchAddress("Harm.HCalScint.hit.ycell", &Harm_HCalScint_hit_ycell, &b_Harm_HCalScint_hit_ycell);
   fChain->SetBranchAddress("Harm.HCalScint.hit.zcell", &Harm_HCalScint_hit_zcell, &b_Harm_HCalScint_hit_zcell);
   fChain->SetBranchAddress("Harm.HCalScint.hit.xcellg", &Harm_HCalScint_hit_xcellg, &b_Harm_HCalScint_hit_xcellg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.ycellg", &Harm_HCalScint_hit_ycellg, &b_Harm_HCalScint_hit_ycellg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.zcellg", &Harm_HCalScint_hit_zcellg, &b_Harm_HCalScint_hit_zcellg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.xhit", &Harm_HCalScint_hit_xhit, &b_Harm_HCalScint_hit_xhit);
   fChain->SetBranchAddress("Harm.HCalScint.hit.yhit", &Harm_HCalScint_hit_yhit, &b_Harm_HCalScint_hit_yhit);
   fChain->SetBranchAddress("Harm.HCalScint.hit.zhit", &Harm_HCalScint_hit_zhit, &b_Harm_HCalScint_hit_zhit);
   fChain->SetBranchAddress("Harm.HCalScint.hit.xhitg", &Harm_HCalScint_hit_xhitg, &b_Harm_HCalScint_hit_xhitg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.yhitg", &Harm_HCalScint_hit_yhitg, &b_Harm_HCalScint_hit_yhitg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.zhitg", &Harm_HCalScint_hit_zhitg, &b_Harm_HCalScint_hit_zhitg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.sumedep", &Harm_HCalScint_hit_sumedep, &b_Harm_HCalScint_hit_sumedep);
   fChain->SetBranchAddress("Harm.HCalScint.hit.tavg", &Harm_HCalScint_hit_tavg, &b_Harm_HCalScint_hit_tavg);
   fChain->SetBranchAddress("Harm.HCalScint.hit.trms", &Harm_HCalScint_hit_trms, &b_Harm_HCalScint_hit_trms);
   fChain->SetBranchAddress("Harm.HCalScint.hit.tmin", &Harm_HCalScint_hit_tmin, &b_Harm_HCalScint_hit_tmin);
   fChain->SetBranchAddress("Harm.HCalScint.hit.tmax", &Harm_HCalScint_hit_tmax, &b_Harm_HCalScint_hit_tmax);
   fChain->SetBranchAddress("primaries.Sigma", &primaries_Sigma, &b_primaries_Sigma);
   fChain->SetBranchAddress("primaries.Ebeam", &primaries_Ebeam, &b_primaries_Ebeam);
   fChain->SetBranchAddress("primaries.Eprime", &primaries_Eprime, &b_primaries_Eprime);
   fChain->SetBranchAddress("primaries.Q2", &primaries_Q2, &b_primaries_Q2);
   fChain->SetBranchAddress("primaries.xbj", &primaries_xbj, &b_primaries_xbj);
   fChain->SetBranchAddress("primaries.y", &primaries_y, &b_primaries_y);
   fChain->SetBranchAddress("primaries.W2", &primaries_W2, &b_primaries_W2);
   fChain->SetBranchAddress("primaries.theta_e", &primaries_theta_e, &b_primaries_theta_e);
   fChain->SetBranchAddress("primaries.phi_e", &primaries_phi_e, &b_primaries_phi_e);
   fChain->SetBranchAddress("primaries.px_e", &primaries_px_e, &b_primaries_px_e);
   fChain->SetBranchAddress("primaries.py_e", &primaries_py_e, &b_primaries_py_e);
   fChain->SetBranchAddress("primaries.pz_e", &primaries_pz_e, &b_primaries_pz_e);
   fChain->SetBranchAddress("primaries.vx_e", &primaries_vx_e, &b_primaries_vx_e);
   fChain->SetBranchAddress("primaries.vy_e", &primaries_vy_e, &b_primaries_vy_e);
   fChain->SetBranchAddress("primaries.vz_e", &primaries_vz_e, &b_primaries_vz_e);
   fChain->SetBranchAddress("primaries.Egamma", &primaries_Egamma, &b_primaries_Egamma);
   fChain->SetBranchAddress("primaries.theta_gamma", &primaries_theta_gamma, &b_primaries_theta_gamma);
   fChain->SetBranchAddress("primaries.phi_gamma", &primaries_phi_gamma, &b_primaries_phi_gamma);
   fChain->SetBranchAddress("primaries.px_gamma", &primaries_px_gamma, &b_primaries_px_gamma);
   fChain->SetBranchAddress("primaries.py_gamma", &primaries_py_gamma, &b_primaries_py_gamma);
   fChain->SetBranchAddress("primaries.pz_gamma", &primaries_pz_gamma, &b_primaries_pz_gamma);
   fChain->SetBranchAddress("primaries.vx_gamma", &primaries_vx_gamma, &b_primaries_vx_gamma);
   fChain->SetBranchAddress("primaries.vy_gamma", &primaries_vy_gamma, &b_primaries_vy_gamma);
   fChain->SetBranchAddress("primaries.vz_gamma", &primaries_vz_gamma, &b_primaries_vz_gamma);
   fChain->SetBranchAddress("Primaries.Nprimaries", &Primaries_Nprimaries, &b_Nprimaries);
   fChain->SetBranchAddress("Primaries.PID", &Primaries_PID, &b_Primaries_PID);
   fChain->SetBranchAddress("Primaries.genflag", &Primaries_genflag, &b_Primaries_genflag);
   fChain->SetBranchAddress("Primaries.Px", &Primaries_Px, &b_Primaries_Px);
   fChain->SetBranchAddress("Primaries.Py", &Primaries_Py, &b_Primaries_Py);
   fChain->SetBranchAddress("Primaries.Pz", &Primaries_Pz, &b_Primaries_Pz);
   fChain->SetBranchAddress("Primaries.vx", &Primaries_vx, &b_Primaries_vx);
   fChain->SetBranchAddress("Primaries.vy", &Primaries_vy, &b_Primaries_vy);
   fChain->SetBranchAddress("Primaries.vz", &Primaries_vz, &b_Primaries_vz);
   fChain->SetBranchAddress("Primaries.M", &Primaries_M, &b_Primaries_M);
   fChain->SetBranchAddress("Primaries.E", &Primaries_E, &b_Primaries_E);
   fChain->SetBranchAddress("Primaries.P", &Primaries_P, &b_Primaries_P);
   fChain->SetBranchAddress("Primaries.t", &Primaries_t, &b_Primaries_t);
   fChain->SetBranchAddress("Primaries.theta", &Primaries_theta, &b_Primaries_theta);
   fChain->SetBranchAddress("Primaries.phi", &Primaries_phi, &b_Primaries_phi);
   Notify();
}

Bool_t gep_pythia6_tree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gep_pythia6_tree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gep_pythia6_tree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gep_pythia6_tree_cxx
