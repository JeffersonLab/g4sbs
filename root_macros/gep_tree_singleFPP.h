//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Mar 22 13:33:14 2021 by ROOT version 6.14/04
// from TTree T/Geant4 SBS Simulation
// found on file: /volatile/halla/sbs/puckett/g4sbs_output/gep_12GeV2_elastic/fppoption1/thick89/gep12_elastic_fppoption1_thick89_job191.root
//////////////////////////////////////////////////////////

#ifndef gep_tree_singleFPP_h
#define gep_tree_singleFPP_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"

class gep_tree_singleFPP {
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
   Double_t        ev_phiS;
   Double_t        ev_thetaS;
   Double_t        ev_MX2;
   Double_t        ev_Sx;
   Double_t        ev_Sy;
   Double_t        ev_Sz;
   Double_t        ev_s;
   Double_t        ev_t;
   Double_t        ev_u;
   Double_t        ev_costhetaCM;
   Double_t        ev_Egamma;
   Int_t           ev_nucl;
   Int_t           ev_fnucl;
   Int_t           ev_hadr;
   Int_t           ev_earmaccept;
   Int_t           ev_harmaccept;
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
   vector<int>     *Earm_CDET_Scint_hit_otridx;
   vector<int>     *Earm_CDET_Scint_hit_ptridx;
   vector<int>     *Earm_CDET_Scint_hit_sdtridx;
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
   vector<int>     *Earm_ECalTF1_hit_otridx;
   vector<int>     *Earm_ECalTF1_hit_ptridx;
   vector<int>     *Earm_ECalTF1_hit_sdtridx;

   //Perhaps we should declare the FPP2 variables as well just so the existing macros don't crash when we use the single-FPP option

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
   vector<int>     *Harm_FPP1_hit_otridx;
   vector<int>     *Harm_FPP1_hit_ptridx;
   vector<int>     *Harm_FPP1_hit_sdtridx;
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
   vector<int>     *Harm_FPP1_Track_otridx;
   vector<int>     *Harm_FPP1_Track_ptridx;
   vector<int>     *Harm_FPP1_Track_sdtridx;

   //Declare the FPP2 variables as dummy variables just to keep compilers happy:
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
   vector<int>     *Harm_FPP2_hit_otridx;
   vector<int>     *Harm_FPP2_hit_ptridx;
   vector<int>     *Harm_FPP2_hit_sdtridx;
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
   vector<int>     *Harm_FPP2_Track_otridx;
   vector<int>     *Harm_FPP2_Track_ptridx;
   vector<int>     *Harm_FPP2_Track_sdtridx;



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
   vector<int>     *Harm_FT_hit_otridx;
   vector<int>     *Harm_FT_hit_ptridx;
   vector<int>     *Harm_FT_hit_sdtridx;
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
   vector<int>     *Harm_FT_Track_otridx;
   vector<int>     *Harm_FT_Track_ptridx;
   vector<int>     *Harm_FT_Track_sdtridx;
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
   vector<int>     *Harm_HCalScint_hit_otridx;
   vector<int>     *Harm_HCalScint_hit_ptridx;
   vector<int>     *Harm_HCalScint_hit_sdtridx;
   Int_t           OTrack_ntracks;
   vector<int>     *OTrack_TID;
   vector<int>     *OTrack_MID;
   vector<int>     *OTrack_PID;
   vector<double>  *OTrack_posx;
   vector<double>  *OTrack_posy;
   vector<double>  *OTrack_posz;
   vector<double>  *OTrack_momx;
   vector<double>  *OTrack_momy;
   vector<double>  *OTrack_momz;
   vector<double>  *OTrack_polx;
   vector<double>  *OTrack_poly;
   vector<double>  *OTrack_polz;
   vector<double>  *OTrack_Etot;
   vector<double>  *OTrack_T;
   Int_t           PTrack_ntracks;
   vector<int>     *PTrack_TID;
   vector<int>     *PTrack_PID;
   vector<double>  *PTrack_posx;
   vector<double>  *PTrack_posy;
   vector<double>  *PTrack_posz;
   vector<double>  *PTrack_momx;
   vector<double>  *PTrack_momy;
   vector<double>  *PTrack_momz;
   vector<double>  *PTrack_polx;
   vector<double>  *PTrack_poly;
   vector<double>  *PTrack_polz;
   vector<double>  *PTrack_Etot;
   vector<double>  *PTrack_T;
   Int_t           SDTrack_ntracks;
   vector<int>     *SDTrack_TID;
   vector<int>     *SDTrack_MID;
   vector<int>     *SDTrack_PID;
   vector<double>  *SDTrack_posx;
   vector<double>  *SDTrack_posy;
   vector<double>  *SDTrack_posz;
   vector<double>  *SDTrack_momx;
   vector<double>  *SDTrack_momy;
   vector<double>  *SDTrack_momz;
   vector<double>  *SDTrack_polx;
   vector<double>  *SDTrack_poly;
   vector<double>  *SDTrack_polz;
   vector<double>  *SDTrack_Etot;
   vector<double>  *SDTrack_T;
   vector<double>  *SDTrack_vx;
   vector<double>  *SDTrack_vy;
   vector<double>  *SDTrack_vz;
   vector<double>  *SDTrack_vnx;
   vector<double>  *SDTrack_vny;
   vector<double>  *SDTrack_vnz;
   vector<double>  *SDTrack_vEkin;

   // List of branches
   TBranch        *b_ev;   //!
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
   TBranch        *b_Earm_CDET_Scint_hit_otridx;   //!
   TBranch        *b_Earm_CDET_Scint_hit_ptridx;   //!
   TBranch        *b_Earm_CDET_Scint_hit_sdtridx;   //!
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
   TBranch        *b_Earm_ECalTF1_hit_otridx;   //!
   TBranch        *b_Earm_ECalTF1_hit_ptridx;   //!
   TBranch        *b_Earm_ECalTF1_hit_sdtridx;   //!
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
   TBranch        *b_Harm_FPP1_hit_otridx;   //!
   TBranch        *b_Harm_FPP1_hit_ptridx;   //!
   TBranch        *b_Harm_FPP1_hit_sdtridx;   //!
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
   TBranch        *b_Harm_FPP1_Track_otridx;   //!
   TBranch        *b_Harm_FPP1_Track_ptridx;   //!
   TBranch        *b_Harm_FPP1_Track_sdtridx;   //!
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
   TBranch        *b_Harm_FT_hit_otridx;   //!
   TBranch        *b_Harm_FT_hit_ptridx;   //!
   TBranch        *b_Harm_FT_hit_sdtridx;   //!
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
   TBranch        *b_Harm_FT_Track_otridx;   //!
   TBranch        *b_Harm_FT_Track_ptridx;   //!
   TBranch        *b_Harm_FT_Track_sdtridx;   //!
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
   TBranch        *b_Harm_HCalScint_hit_otridx;   //!
   TBranch        *b_Harm_HCalScint_hit_ptridx;   //!
   TBranch        *b_Harm_HCalScint_hit_sdtridx;   //!
   TBranch        *b_OTrack_ntracks;   //!
   TBranch        *b_OTrack_TID;   //!
   TBranch        *b_OTrack_MID;   //!
   TBranch        *b_OTrack_PID;   //!
   TBranch        *b_OTrack_posx;   //!
   TBranch        *b_OTrack_posy;   //!
   TBranch        *b_OTrack_posz;   //!
   TBranch        *b_OTrack_momx;   //!
   TBranch        *b_OTrack_momy;   //!
   TBranch        *b_OTrack_momz;   //!
   TBranch        *b_OTrack_polx;   //!
   TBranch        *b_OTrack_poly;   //!
   TBranch        *b_OTrack_polz;   //!
   TBranch        *b_OTrack_Etot;   //!
   TBranch        *b_OTrack_T;   //!
   TBranch        *b_PTrack_ntracks;   //!
   TBranch        *b_PTrack_TID;   //!
   TBranch        *b_PTrack_PID;   //!
   TBranch        *b_PTrack_posx;   //!
   TBranch        *b_PTrack_posy;   //!
   TBranch        *b_PTrack_posz;   //!
   TBranch        *b_PTrack_momx;   //!
   TBranch        *b_PTrack_momy;   //!
   TBranch        *b_PTrack_momz;   //!
   TBranch        *b_PTrack_polx;   //!
   TBranch        *b_PTrack_poly;   //!
   TBranch        *b_PTrack_polz;   //!
   TBranch        *b_PTrack_Etot;   //!
   TBranch        *b_PTrack_T;   //!
   TBranch        *b_SDTrack_ntracks;   //!
   TBranch        *b_SDTrack_TID;   //!
   TBranch        *b_SDTrack_MID;   //!
   TBranch        *b_SDTrack_PID;   //!
   TBranch        *b_SDTrack_posx;   //!
   TBranch        *b_SDTrack_posy;   //!
   TBranch        *b_SDTrack_posz;   //!
   TBranch        *b_SDTrack_momx;   //!
   TBranch        *b_SDTrack_momy;   //!
   TBranch        *b_SDTrack_momz;   //!
   TBranch        *b_SDTrack_polx;   //!
   TBranch        *b_SDTrack_poly;   //!
   TBranch        *b_SDTrack_polz;   //!
   TBranch        *b_SDTrack_Etot;   //!
   TBranch        *b_SDTrack_T;   //!
   TBranch        *b_SDTrack_vx;   //!
   TBranch        *b_SDTrack_vy;   //!
   TBranch        *b_SDTrack_vz;   //!
   TBranch        *b_SDTrack_vnx;   //!
   TBranch        *b_SDTrack_vny;   //!
   TBranch        *b_SDTrack_vnz;   //!
   TBranch        *b_SDTrack_vEkin;   //!

   gep_tree_singleFPP(TTree *tree=0);
   virtual ~gep_tree_singleFPP();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef gep_tree_singleFPP_cxx
gep_tree_singleFPP::gep_tree_singleFPP(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/volatile/halla/sbs/puckett/g4sbs_output/gep_12GeV2_elastic/fppoption1/thick89/gep12_elastic_fppoption1_thick89_job191.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/volatile/halla/sbs/puckett/g4sbs_output/gep_12GeV2_elastic/fppoption1/thick89/gep12_elastic_fppoption1_thick89_job191.root");
      }
      f->GetObject("T",tree);

   }
   Init(tree);
}

gep_tree_singleFPP::~gep_tree_singleFPP()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gep_tree_singleFPP::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gep_tree_singleFPP::LoadTree(Long64_t entry)
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

void gep_tree_singleFPP::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
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
   Earm_CDET_Scint_hit_otridx = 0;
   Earm_CDET_Scint_hit_ptridx = 0;
   Earm_CDET_Scint_hit_sdtridx = 0;
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
   Earm_ECalTF1_hit_otridx = 0;
   Earm_ECalTF1_hit_ptridx = 0;
   Earm_ECalTF1_hit_sdtridx = 0;
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
   Harm_FPP1_hit_otridx = 0;
   Harm_FPP1_hit_ptridx = 0;
   Harm_FPP1_hit_sdtridx = 0;
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
   Harm_FPP1_Track_otridx = 0;
   Harm_FPP1_Track_ptridx = 0;
   Harm_FPP1_Track_sdtridx = 0;

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
   Harm_FPP2_hit_otridx = 0;
   Harm_FPP2_hit_ptridx = 0;
   Harm_FPP2_hit_sdtridx = 0;
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
   Harm_FPP2_Track_otridx = 0;
   Harm_FPP2_Track_ptridx = 0;
   Harm_FPP2_Track_sdtridx = 0;


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
   Harm_FT_hit_otridx = 0;
   Harm_FT_hit_ptridx = 0;
   Harm_FT_hit_sdtridx = 0;
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
   Harm_FT_Track_otridx = 0;
   Harm_FT_Track_ptridx = 0;
   Harm_FT_Track_sdtridx = 0;
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
   Harm_HCalScint_hit_otridx = 0;
   Harm_HCalScint_hit_ptridx = 0;
   Harm_HCalScint_hit_sdtridx = 0;
   OTrack_TID = 0;
   OTrack_MID = 0;
   OTrack_PID = 0;
   OTrack_posx = 0;
   OTrack_posy = 0;
   OTrack_posz = 0;
   OTrack_momx = 0;
   OTrack_momy = 0;
   OTrack_momz = 0;
   OTrack_polx = 0;
   OTrack_poly = 0;
   OTrack_polz = 0;
   OTrack_Etot = 0;
   OTrack_T = 0;
   PTrack_TID = 0;
   PTrack_PID = 0;
   PTrack_posx = 0;
   PTrack_posy = 0;
   PTrack_posz = 0;
   PTrack_momx = 0;
   PTrack_momy = 0;
   PTrack_momz = 0;
   PTrack_polx = 0;
   PTrack_poly = 0;
   PTrack_polz = 0;
   PTrack_Etot = 0;
   PTrack_T = 0;
   SDTrack_TID = 0;
   SDTrack_MID = 0;
   SDTrack_PID = 0;
   SDTrack_posx = 0;
   SDTrack_posy = 0;
   SDTrack_posz = 0;
   SDTrack_momx = 0;
   SDTrack_momy = 0;
   SDTrack_momz = 0;
   SDTrack_polx = 0;
   SDTrack_poly = 0;
   SDTrack_polz = 0;
   SDTrack_Etot = 0;
   SDTrack_T = 0;
   SDTrack_vx = 0;
   SDTrack_vy = 0;
   SDTrack_vz = 0;
   SDTrack_vnx = 0;
   SDTrack_vny = 0;
   SDTrack_vnz = 0;
   SDTrack_vEkin = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("ev", &ev_count, &b_ev);
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
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.otridx", &Earm_CDET_Scint_hit_otridx, &b_Earm_CDET_Scint_hit_otridx);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.ptridx", &Earm_CDET_Scint_hit_ptridx, &b_Earm_CDET_Scint_hit_ptridx);
   fChain->SetBranchAddress("Earm.CDET_Scint.hit.sdtridx", &Earm_CDET_Scint_hit_sdtridx, &b_Earm_CDET_Scint_hit_sdtridx);
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
   fChain->SetBranchAddress("Earm.ECalTF1.hit.otridx", &Earm_ECalTF1_hit_otridx, &b_Earm_ECalTF1_hit_otridx);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.ptridx", &Earm_ECalTF1_hit_ptridx, &b_Earm_ECalTF1_hit_ptridx);
   fChain->SetBranchAddress("Earm.ECalTF1.hit.sdtridx", &Earm_ECalTF1_hit_sdtridx, &b_Earm_ECalTF1_hit_sdtridx);
   fChain->SetBranchAddress("Harm.FPP.hit.nhits", &Harm_FPP1_hit_nhits, &b_Harm_FPP1_hit_nhits);
   fChain->SetBranchAddress("Harm.FPP.hit.plane", &Harm_FPP1_hit_plane, &b_Harm_FPP1_hit_plane);
   fChain->SetBranchAddress("Harm.FPP.hit.strip", &Harm_FPP1_hit_strip, &b_Harm_FPP1_hit_strip);
   fChain->SetBranchAddress("Harm.FPP.hit.x", &Harm_FPP1_hit_x, &b_Harm_FPP1_hit_x);
   fChain->SetBranchAddress("Harm.FPP.hit.y", &Harm_FPP1_hit_y, &b_Harm_FPP1_hit_y);
   fChain->SetBranchAddress("Harm.FPP.hit.z", &Harm_FPP1_hit_z, &b_Harm_FPP1_hit_z);
   fChain->SetBranchAddress("Harm.FPP.hit.polx", &Harm_FPP1_hit_polx, &b_Harm_FPP1_hit_polx);
   fChain->SetBranchAddress("Harm.FPP.hit.poly", &Harm_FPP1_hit_poly, &b_Harm_FPP1_hit_poly);
   fChain->SetBranchAddress("Harm.FPP.hit.polz", &Harm_FPP1_hit_polz, &b_Harm_FPP1_hit_polz);
   fChain->SetBranchAddress("Harm.FPP.hit.t", &Harm_FPP1_hit_t, &b_Harm_FPP1_hit_t);
   fChain->SetBranchAddress("Harm.FPP.hit.trms", &Harm_FPP1_hit_trms, &b_Harm_FPP1_hit_trms);
   fChain->SetBranchAddress("Harm.FPP.hit.tmin", &Harm_FPP1_hit_tmin, &b_Harm_FPP1_hit_tmin);
   fChain->SetBranchAddress("Harm.FPP.hit.tmax", &Harm_FPP1_hit_tmax, &b_Harm_FPP1_hit_tmax);
   fChain->SetBranchAddress("Harm.FPP.hit.tx", &Harm_FPP1_hit_tx, &b_Harm_FPP1_hit_tx);
   fChain->SetBranchAddress("Harm.FPP.hit.ty", &Harm_FPP1_hit_ty, &b_Harm_FPP1_hit_ty);
   fChain->SetBranchAddress("Harm.FPP.hit.xin", &Harm_FPP1_hit_xin, &b_Harm_FPP1_hit_xin);
   fChain->SetBranchAddress("Harm.FPP.hit.yin", &Harm_FPP1_hit_yin, &b_Harm_FPP1_hit_yin);
   fChain->SetBranchAddress("Harm.FPP.hit.zin", &Harm_FPP1_hit_zin, &b_Harm_FPP1_hit_zin);
   fChain->SetBranchAddress("Harm.FPP.hit.xout", &Harm_FPP1_hit_xout, &b_Harm_FPP1_hit_xout);
   fChain->SetBranchAddress("Harm.FPP.hit.yout", &Harm_FPP1_hit_yout, &b_Harm_FPP1_hit_yout);
   fChain->SetBranchAddress("Harm.FPP.hit.zout", &Harm_FPP1_hit_zout, &b_Harm_FPP1_hit_zout);
   fChain->SetBranchAddress("Harm.FPP.hit.txp", &Harm_FPP1_hit_txp, &b_Harm_FPP1_hit_txp);
   fChain->SetBranchAddress("Harm.FPP.hit.typ", &Harm_FPP1_hit_typ, &b_Harm_FPP1_hit_typ);
   fChain->SetBranchAddress("Harm.FPP.hit.xg", &Harm_FPP1_hit_xg, &b_Harm_FPP1_hit_xg);
   fChain->SetBranchAddress("Harm.FPP.hit.yg", &Harm_FPP1_hit_yg, &b_Harm_FPP1_hit_yg);
   fChain->SetBranchAddress("Harm.FPP.hit.zg", &Harm_FPP1_hit_zg, &b_Harm_FPP1_hit_zg);
   fChain->SetBranchAddress("Harm.FPP.hit.trid", &Harm_FPP1_hit_trid, &b_Harm_FPP1_hit_trid);
   fChain->SetBranchAddress("Harm.FPP.hit.mid", &Harm_FPP1_hit_mid, &b_Harm_FPP1_hit_mid);
   fChain->SetBranchAddress("Harm.FPP.hit.pid", &Harm_FPP1_hit_pid, &b_Harm_FPP1_hit_pid);
   fChain->SetBranchAddress("Harm.FPP.hit.vx", &Harm_FPP1_hit_vx, &b_Harm_FPP1_hit_vx);
   fChain->SetBranchAddress("Harm.FPP.hit.vy", &Harm_FPP1_hit_vy, &b_Harm_FPP1_hit_vy);
   fChain->SetBranchAddress("Harm.FPP.hit.vz", &Harm_FPP1_hit_vz, &b_Harm_FPP1_hit_vz);
   fChain->SetBranchAddress("Harm.FPP.hit.p", &Harm_FPP1_hit_p, &b_Harm_FPP1_hit_p);
   fChain->SetBranchAddress("Harm.FPP.hit.edep", &Harm_FPP1_hit_edep, &b_Harm_FPP1_hit_edep);
   fChain->SetBranchAddress("Harm.FPP.hit.beta", &Harm_FPP1_hit_beta, &b_Harm_FPP1_hit_beta);
   fChain->SetBranchAddress("Harm.FPP.hit.otridx", &Harm_FPP1_hit_otridx, &b_Harm_FPP1_hit_otridx);
   fChain->SetBranchAddress("Harm.FPP.hit.ptridx", &Harm_FPP1_hit_ptridx, &b_Harm_FPP1_hit_ptridx);
   fChain->SetBranchAddress("Harm.FPP.hit.sdtridx", &Harm_FPP1_hit_sdtridx, &b_Harm_FPP1_hit_sdtridx);
   fChain->SetBranchAddress("Harm.FPP.Track.ntracks", &Harm_FPP1_Track_ntracks, &b_Harm_FPP1_Track_ntracks);
   fChain->SetBranchAddress("Harm.FPP.Track.TID", &Harm_FPP1_Track_TID, &b_Harm_FPP1_Track_TID);
   fChain->SetBranchAddress("Harm.FPP.Track.PID", &Harm_FPP1_Track_PID, &b_Harm_FPP1_Track_PID);
   fChain->SetBranchAddress("Harm.FPP.Track.MID", &Harm_FPP1_Track_MID, &b_Harm_FPP1_Track_MID);
   fChain->SetBranchAddress("Harm.FPP.Track.NumHits", &Harm_FPP1_Track_NumHits, &b_Harm_FPP1_Track_NumHits);
   fChain->SetBranchAddress("Harm.FPP.Track.NumPlanes", &Harm_FPP1_Track_NumPlanes, &b_Harm_FPP1_Track_NumPlanes);
   fChain->SetBranchAddress("Harm.FPP.Track.NDF", &Harm_FPP1_Track_NDF, &b_Harm_FPP1_Track_NDF);
   fChain->SetBranchAddress("Harm.FPP.Track.Chi2fit", &Harm_FPP1_Track_Chi2fit, &b_Harm_FPP1_Track_Chi2fit);
   fChain->SetBranchAddress("Harm.FPP.Track.Chi2true", &Harm_FPP1_Track_Chi2true, &b_Harm_FPP1_Track_Chi2true);
   fChain->SetBranchAddress("Harm.FPP.Track.X", &Harm_FPP1_Track_X, &b_Harm_FPP1_Track_X);
   fChain->SetBranchAddress("Harm.FPP.Track.Y", &Harm_FPP1_Track_Y, &b_Harm_FPP1_Track_Y);
   fChain->SetBranchAddress("Harm.FPP.Track.Xp", &Harm_FPP1_Track_Xp, &b_Harm_FPP1_Track_Xp);
   fChain->SetBranchAddress("Harm.FPP.Track.Yp", &Harm_FPP1_Track_Yp, &b_Harm_FPP1_Track_Yp);
   fChain->SetBranchAddress("Harm.FPP.Track.T", &Harm_FPP1_Track_T, &b_Harm_FPP1_Track_T);
   fChain->SetBranchAddress("Harm.FPP.Track.P", &Harm_FPP1_Track_P, &b_Harm_FPP1_Track_P);
   fChain->SetBranchAddress("Harm.FPP.Track.Sx", &Harm_FPP1_Track_Sx, &b_Harm_FPP1_Track_Sx);
   fChain->SetBranchAddress("Harm.FPP.Track.Sy", &Harm_FPP1_Track_Sy, &b_Harm_FPP1_Track_Sy);
   fChain->SetBranchAddress("Harm.FPP.Track.Sz", &Harm_FPP1_Track_Sz, &b_Harm_FPP1_Track_Sz);
   fChain->SetBranchAddress("Harm.FPP.Track.Xfit", &Harm_FPP1_Track_Xfit, &b_Harm_FPP1_Track_Xfit);
   fChain->SetBranchAddress("Harm.FPP.Track.Yfit", &Harm_FPP1_Track_Yfit, &b_Harm_FPP1_Track_Yfit);
   fChain->SetBranchAddress("Harm.FPP.Track.Xpfit", &Harm_FPP1_Track_Xpfit, &b_Harm_FPP1_Track_Xpfit);
   fChain->SetBranchAddress("Harm.FPP.Track.Ypfit", &Harm_FPP1_Track_Ypfit, &b_Harm_FPP1_Track_Ypfit);
   fChain->SetBranchAddress("Harm.FPP.Track.otridx", &Harm_FPP1_Track_otridx, &b_Harm_FPP1_Track_otridx);
   fChain->SetBranchAddress("Harm.FPP.Track.ptridx", &Harm_FPP1_Track_ptridx, &b_Harm_FPP1_Track_ptridx);
   fChain->SetBranchAddress("Harm.FPP.Track.sdtridx", &Harm_FPP1_Track_sdtridx, &b_Harm_FPP1_Track_sdtridx);
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
   fChain->SetBranchAddress("Harm.FT.hit.otridx", &Harm_FT_hit_otridx, &b_Harm_FT_hit_otridx);
   fChain->SetBranchAddress("Harm.FT.hit.ptridx", &Harm_FT_hit_ptridx, &b_Harm_FT_hit_ptridx);
   fChain->SetBranchAddress("Harm.FT.hit.sdtridx", &Harm_FT_hit_sdtridx, &b_Harm_FT_hit_sdtridx);
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
   fChain->SetBranchAddress("Harm.FT.Track.otridx", &Harm_FT_Track_otridx, &b_Harm_FT_Track_otridx);
   fChain->SetBranchAddress("Harm.FT.Track.ptridx", &Harm_FT_Track_ptridx, &b_Harm_FT_Track_ptridx);
   fChain->SetBranchAddress("Harm.FT.Track.sdtridx", &Harm_FT_Track_sdtridx, &b_Harm_FT_Track_sdtridx);
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
   fChain->SetBranchAddress("Harm.HCalScint.hit.otridx", &Harm_HCalScint_hit_otridx, &b_Harm_HCalScint_hit_otridx);
   fChain->SetBranchAddress("Harm.HCalScint.hit.ptridx", &Harm_HCalScint_hit_ptridx, &b_Harm_HCalScint_hit_ptridx);
   fChain->SetBranchAddress("Harm.HCalScint.hit.sdtridx", &Harm_HCalScint_hit_sdtridx, &b_Harm_HCalScint_hit_sdtridx);
   fChain->SetBranchAddress("OTrack.ntracks", &OTrack_ntracks, &b_OTrack_ntracks);
   fChain->SetBranchAddress("OTrack.TID", &OTrack_TID, &b_OTrack_TID);
   fChain->SetBranchAddress("OTrack.MID", &OTrack_MID, &b_OTrack_MID);
   fChain->SetBranchAddress("OTrack.PID", &OTrack_PID, &b_OTrack_PID);
   fChain->SetBranchAddress("OTrack.posx", &OTrack_posx, &b_OTrack_posx);
   fChain->SetBranchAddress("OTrack.posy", &OTrack_posy, &b_OTrack_posy);
   fChain->SetBranchAddress("OTrack.posz", &OTrack_posz, &b_OTrack_posz);
   fChain->SetBranchAddress("OTrack.momx", &OTrack_momx, &b_OTrack_momx);
   fChain->SetBranchAddress("OTrack.momy", &OTrack_momy, &b_OTrack_momy);
   fChain->SetBranchAddress("OTrack.momz", &OTrack_momz, &b_OTrack_momz);
   fChain->SetBranchAddress("OTrack.polx", &OTrack_polx, &b_OTrack_polx);
   fChain->SetBranchAddress("OTrack.poly", &OTrack_poly, &b_OTrack_poly);
   fChain->SetBranchAddress("OTrack.polz", &OTrack_polz, &b_OTrack_polz);
   fChain->SetBranchAddress("OTrack.Etot", &OTrack_Etot, &b_OTrack_Etot);
   fChain->SetBranchAddress("OTrack.T", &OTrack_T, &b_OTrack_T);
   fChain->SetBranchAddress("PTrack.ntracks", &PTrack_ntracks, &b_PTrack_ntracks);
   fChain->SetBranchAddress("PTrack.TID", &PTrack_TID, &b_PTrack_TID);
   fChain->SetBranchAddress("PTrack.PID", &PTrack_PID, &b_PTrack_PID);
   fChain->SetBranchAddress("PTrack.posx", &PTrack_posx, &b_PTrack_posx);
   fChain->SetBranchAddress("PTrack.posy", &PTrack_posy, &b_PTrack_posy);
   fChain->SetBranchAddress("PTrack.posz", &PTrack_posz, &b_PTrack_posz);
   fChain->SetBranchAddress("PTrack.momx", &PTrack_momx, &b_PTrack_momx);
   fChain->SetBranchAddress("PTrack.momy", &PTrack_momy, &b_PTrack_momy);
   fChain->SetBranchAddress("PTrack.momz", &PTrack_momz, &b_PTrack_momz);
   fChain->SetBranchAddress("PTrack.polx", &PTrack_polx, &b_PTrack_polx);
   fChain->SetBranchAddress("PTrack.poly", &PTrack_poly, &b_PTrack_poly);
   fChain->SetBranchAddress("PTrack.polz", &PTrack_polz, &b_PTrack_polz);
   fChain->SetBranchAddress("PTrack.Etot", &PTrack_Etot, &b_PTrack_Etot);
   fChain->SetBranchAddress("PTrack.T", &PTrack_T, &b_PTrack_T);
   fChain->SetBranchAddress("SDTrack.ntracks", &SDTrack_ntracks, &b_SDTrack_ntracks);
   fChain->SetBranchAddress("SDTrack.TID", &SDTrack_TID, &b_SDTrack_TID);
   fChain->SetBranchAddress("SDTrack.MID", &SDTrack_MID, &b_SDTrack_MID);
   fChain->SetBranchAddress("SDTrack.PID", &SDTrack_PID, &b_SDTrack_PID);
   fChain->SetBranchAddress("SDTrack.posx", &SDTrack_posx, &b_SDTrack_posx);
   fChain->SetBranchAddress("SDTrack.posy", &SDTrack_posy, &b_SDTrack_posy);
   fChain->SetBranchAddress("SDTrack.posz", &SDTrack_posz, &b_SDTrack_posz);
   fChain->SetBranchAddress("SDTrack.momx", &SDTrack_momx, &b_SDTrack_momx);
   fChain->SetBranchAddress("SDTrack.momy", &SDTrack_momy, &b_SDTrack_momy);
   fChain->SetBranchAddress("SDTrack.momz", &SDTrack_momz, &b_SDTrack_momz);
   fChain->SetBranchAddress("SDTrack.polx", &SDTrack_polx, &b_SDTrack_polx);
   fChain->SetBranchAddress("SDTrack.poly", &SDTrack_poly, &b_SDTrack_poly);
   fChain->SetBranchAddress("SDTrack.polz", &SDTrack_polz, &b_SDTrack_polz);
   fChain->SetBranchAddress("SDTrack.Etot", &SDTrack_Etot, &b_SDTrack_Etot);
   fChain->SetBranchAddress("SDTrack.T", &SDTrack_T, &b_SDTrack_T);
   fChain->SetBranchAddress("SDTrack.vx", &SDTrack_vx, &b_SDTrack_vx);
   fChain->SetBranchAddress("SDTrack.vy", &SDTrack_vy, &b_SDTrack_vy);
   fChain->SetBranchAddress("SDTrack.vz", &SDTrack_vz, &b_SDTrack_vz);
   fChain->SetBranchAddress("SDTrack.vnx", &SDTrack_vnx, &b_SDTrack_vnx);
   fChain->SetBranchAddress("SDTrack.vny", &SDTrack_vny, &b_SDTrack_vny);
   fChain->SetBranchAddress("SDTrack.vnz", &SDTrack_vnz, &b_SDTrack_vnz);
   fChain->SetBranchAddress("SDTrack.vEkin", &SDTrack_vEkin, &b_SDTrack_vEkin);
   Notify();
}

Bool_t gep_tree_singleFPP::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gep_tree_singleFPP::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gep_tree_singleFPP::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef gep_tree_singleFPP_cxx
