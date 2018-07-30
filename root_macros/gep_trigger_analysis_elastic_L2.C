#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "gep_tree_elastic.C"
#include "G4SBSRunData.hh"
#include "TRandom3.h"
#include "TChainElement.h"
#include "TObjArray.h"
#include "TString.h"
#include "TObjString.h"
//#include "TIter.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"

#include "TVector3.h"

#include "TF1.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>

using namespace std;

TF1 *gaussian = new TF1("gaussian", "[0]*exp(-0.5*pow((x-[1])/[2],2))", 0.0,5000.0);
TF1 *gaussplusexpo = new TF1("gaussplusexpo", "[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-[4]*x)",0.0,5000.0);

//const double Rcal = 4.9; //m, to back surface
const double Mp = 0.938272; //GeV
const double PI = TMath::Pi();

void gep_trigger_analysis_elastic_L2( const char *rootfilename, const char *thresholdfilename_ecal="database/ecal_trigger_thresholds_12GeV2.txt", const char *outputfilename="gep_L2_trigger_analysis_elastic_temp.root", int pheflag=0, const char *assocfilename="ECAL_HCAL_L2_correlations_nophe.txt"){

  double nominal_threshold_HCAL = 0.5;
  double nominal_threshold_ECAL = 0.85;
  
  TFile *fout = new TFile(outputfilename,"RECREATE");
  TChain *C = new TChain("T");
  C->Add(rootfilename);

  gep_tree_elastic *T = new gep_tree_elastic( C );

  G4SBSRunData *rd;

  long ngen = 0;
  int nfiles = 0;
  
  TObjArray *FileList = C->GetListOfFiles();
  TIter next(FileList);

  TChainElement *chEl = 0;

  set<TString> bad_file_list;

  map<TString,double> ECALdist_file;
  map<TString,double> ECALtheta_file;
  
  while( (chEl=(TChainElement*)next() )){
    TFile newfile(chEl->GetTitle());
    newfile.GetObject("run_data",rd);
    if( rd ){
      ngen += rd->fNtries;
      nfiles++;

      ECALdist_file[chEl->GetTitle()] = rd->fBBdist;
      ECALtheta_file[chEl->GetTitle()] = rd->fBBtheta;

    } else {
      bad_file_list.insert( chEl->GetTitle());
    }
  }

  cout << "number of generated events = " << ngen << endl;
  
  set<int> list_of_nodes_ecal;
  map<int, set<int> > cells_logic_sums_ecal; //mapping between node numbers and cell numbers
  map<int, double> logic_mean_ecal; //mean peak positions by node number
  map<int, double> logic_sigma_ecal; //peak width by node number
  map<int, double> threshold_ecal; //threshold by node number
  map<std::pair<int,int>, int > cell_rowcol_ecal; //cell numbers mapped by unique row and column pairs
  map<int,set<int> > nodes_cells_ecal; //mapping of nodes by cell number:
  map<int,int> rows_cells_ecal;
  map<int,int> cols_cells_ecal;
  map<int,double> xcells_ecal;
  map<int,double> ycells_ecal;
  
  //keep track of min and max x by row number:
  double ycellmin,ycellmax;
  map<int,double> ycell_rows;
  //map<int,double> cellsize_rows;
  map<int,double> xcellmin_rows;
  map<int,double> xcellmax_rows;
  
  int minrow=1000,maxrow=-1;
  
  set<int> rows_ecal;
  map<int,set<int> > columns_rows_ecal;
  
  map<int,double> elastic_peak_new_ecal;
  map<int,double> sigma_new_ecal;
  map<int,double> threshold_new_ecal;

  ifstream mapfile_ecal("database/ecal_gep_blockmap.txt");
  ifstream logicfile_ecal("database/GEP_ECAL_L2sums.txt");
  //ifstream thresholdfile(thresholdfilename);

  TString currentline;

  ycellmin = 1.0e9;
  ycellmax = -1.0e9;


  //Get cell mapping info:
  while( currentline.ReadLine( mapfile_ecal ) ){
    cout << currentline << endl;
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(",");
      int ntokens = tokens->GetEntries();

      if( ntokens == 5 ){
	TString scell = ( (TObjString*) (*tokens)[0] )->GetString();
	TString srow  = ( (TObjString*) (*tokens)[1] )->GetString();
	TString scol  = ( (TObjString*) (*tokens)[2] )->GetString();
	TString sxcell = ( (TObjString*) (*tokens)[3] )->GetString();
	TString sycell = ( (TObjString*) (*tokens)[4] )->GetString();

	int cellnum = scell.Atoi();
	int row = srow.Atoi();
	int col = scol.Atoi();
	double xcell = sxcell.Atoi();
	double ycell = sycell.Atoi();

	std::pair<int, int> rowcoltemp(row,col);

	cell_rowcol_ecal[rowcoltemp] = cellnum;
	rows_cells_ecal[cellnum] = row;
	cols_cells_ecal[cellnum] = col;
	xcells_ecal[cellnum] = xcell/100.0;
	ycells_ecal[cellnum] = ycell/100.0;

	if( ycell_rows.empty() || ycell/100.0 < ycellmin ) ycellmin = ycell/100.0;
	if( ycell_rows.empty() || ycell/100.0 > ycellmax ) ycellmax = ycell/100.0;
	  
	ycell_rows[row] = ycell/100.0;

	// if( row < 30 ) {
	//   cellsize_rows[row] = 4.0; //cm
	// } else {
	//   cellsize_rows[row] = 4.2; //cm
	// }
	
	if( xcellmin_rows.empty() || xcell/100.0 < xcellmin_rows[row] ){
	  xcellmin_rows[row] = xcell/100.0;
	}

	if( xcellmax_rows.empty() || xcell/100.0 > xcellmax_rows[row] ){
	  xcellmax_rows[row] = xcell/100.0;
	}
      }
    }
  }
  
  int current_node = 1;

  bool first_cell = true;

  //Get trigger logic summing info:
  
  while( currentline.ReadLine(logicfile_ecal) ){
    if( !currentline.BeginsWith( "#" ) ){
      
      TObjArray *tokens = currentline.Tokenize(" ");
      int ntokens = tokens->GetEntries();

      TString snode = ( (TObjString*) (*tokens)[0] )->GetString();
      int nodenumber = snode.Atoi() + 1;
      
      TString sncell_node = ( (TObjString*) (*tokens)[1] )->GetString();
      int ncell_node = sncell_node.Atoi();
      
      if( ntokens >= ncell_node + 2 ){
	list_of_nodes_ecal.insert( nodenumber );
	for( int itoken = 2; itoken < ncell_node+2; itoken++ ){
	  TString scell = ( (TObjString*) (*tokens)[itoken] )->GetString();
	  int cell = scell.Atoi();

	  cells_logic_sums_ecal[nodenumber].insert( cell );

	  //provide default values;
	  logic_mean_ecal[nodenumber] = 3000.0;
	  logic_sigma_ecal[nodenumber] = 0.06*3000.0;
	  threshold_ecal[nodenumber] = 0.8;
	  
	  nodes_cells_ecal[ cell ].insert( nodenumber );
	}
      }
    }
  }

  set<int> list_of_nodes_hcal;
  map<int, set<int> > cells_logic_sums_hcal; //mapping between node numbers and cell numbers
  map<int, double> logic_mean_hcal; //mean peak positions by node number
  map<int, double> logic_sigma_hcal; //peak width by node number
  map<int, double> threshold_hcal; //threshold by node number
  map<std::pair<int,int>, int > cell_rowcol_hcal; //cell numbers mapped by unique row and column pairs
  map<int,set<int> > nodes_cells_hcal; //mapping of nodes by cell number:

  current_node = 1;
  //  bool first_cell = true;

  int nrows_hcal=24;
  int ncols_hcal=12;
  for( int row=1; row<=nrows_hcal-3; row++ ){
    for( int col=1; col<=ncols_hcal-3; col++ ){
      list_of_nodes_hcal.insert( current_node );
      for( int m=col; m<=col+3; m++ ){ //Add all blocks in each 4x4 sum to the current node:
	for( int n=row; n<=row+3; n++ ){ //Add
	  int cell = (n-1) + nrows_hcal*(m-1) + 1;
	  cells_logic_sums_hcal[current_node].insert( cell );
	  nodes_cells_hcal[cell].insert(current_node);
	  std::pair<int,int> rowcol(n,m);
	  cell_rowcol_hcal[rowcol] = cell;
	  //Use some sensible default value for logic mean and logic sigma:
	}
      }

      logic_mean_hcal[current_node] = 1141.5;
      logic_sigma_hcal[current_node] = 356.25;
      
      current_node++;
    }
  }
  
  
  TH1D::SetDefaultSumw2();

  //double PI = TMath::Pi();

  //Photoelectron statistics:
  double phe_per_GeV_ECAL = 1000.0/1.33; //~ 750 pe/GeV
  double phe_per_GeV_HCAL = 1000.0/0.30; //~ 3,333 pe/GeV (but sampling fraction is small)

  //read in alternate threshold:
  ifstream thresholdfile_ecal(thresholdfilename_ecal);
  if( thresholdfile_ecal ){
    int node;
    double mean,sigma;
    while( thresholdfile_ecal >> node >> mean >> sigma ){
      if( list_of_nodes_ecal.find( node ) != list_of_nodes_ecal.end() ){
	logic_mean_ecal[ node ] = mean;
	logic_sigma_ecal[ node ] = sigma;
      }
    }
  }

  //read in alternate threshold:
  // ifstream thresholdfile_hcal(thresholdfilename_hcal);
  // if( thresholdfile_hcal ){
  //   int node;
  //   double mean,sigma;
  //   while( thresholdfile_hcal >> node >> mean >> sigma ){
  //     if( list_of_nodes_hcal.find( node ) != list_of_nodes_hcal.end() ){
  // 	logic_mean_hcal[ node ] = mean;
  // 	logic_sigma_hcal[ node ] = sigma;
  //     }
  //   }
  // }

  ifstream assocfile( assocfilename );

  bool use_ECAL_HCAL_associations=false;
  map<int, set<int> > ECAL_nodes_HCAL;
  if( assocfile ){
    while( !assocfile.eof() ){
      int hcalnode, N;
      assocfile >> hcalnode >> N;
      for( int i=0; i<N; i++ ){
	int ecalnode;
	assocfile >> ecalnode;

	ECAL_nodes_HCAL[hcalnode].insert(ecalnode);
      }
    }
  }
  
  fout->cd();

  // TH1D *hrate_vs_threshold_ECAL = new TH1D("hrate_vs_threshold_ECAL","",30,0.0,1.5);
  // //TH1D *hnum_logic_sums_fired_vs_threshold = new TH1D("hnum_logic_sums_fired_vs_threshold

  // TH1D *hrate_vs_threshold_HCAL = new TH1D("hrate_vs_threshold_HCAL","",40,0.0,2.0);

  //TH2D *htrue_coincidence_rate_vs_threshold_ECAL_HCAL = new TH2D("htrue_coincidence_rate_vs_threshold_ECAL_HCAL","",40,0,2.0,30,0,1.5);
  
  TH2D *hnphesum_vs_node_ECAL_all = new TH2D("hnphesum_vs_node_ECAL_all","",list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5,100,0.0,5000.0);
  TH2D *hnphesum_vs_node_HCAL_all = new TH2D("hnphesum_vs_node_HCAL_all","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);

  //TH2D *hnphesum_vs_node_ECAL_FTcut = new TH2D("hnphesum_vs_node_ECAL_FTcut","",list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5,100,0.0,5000.0);
  TH2D *hnphesum_vs_node_HCAL_FTcut = new TH2D("hnphesum_vs_node_HCAL_FTcut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCAL_FPP1cut = new TH2D("hnphesum_vs_node_HCAL_FPP1cut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCAL_FPP2cut = new TH2D("hnphesum_vs_node_HCAL_FPP2cut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCAL_FPPbothcut = new TH2D("hnphesum_vs_node_HCAL_FPPbothcut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCAL_FPPeithercut = new TH2D("hnphesum_vs_node_HCAL_FPPeithercut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);

  TH2D *hnphesum_vs_node_HCALmax_all = new TH2D("hnphesum_vs_node_HCALmax_all","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCALmax_FTcut = new TH2D("hnphesum_vs_node_HCALmax_FTcut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCALmax_FPP1cut = new TH2D("hnphesum_vs_node_HCALmax_FPP1cut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCALmax_FPP2cut = new TH2D("hnphesum_vs_node_HCALmax_FPP2cut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCALmax_FPPbothcut = new TH2D("hnphesum_vs_node_HCALmax_FPPbothcut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCALmax_FPPeithercut = new TH2D("hnphesum_vs_node_HCALmax_FPPeithercut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  
  
  TH2D *hmaxnode_ECAL_vs_HCAL = new TH2D("hmaxnode_ECAL_vs_HCAL","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5);
  TH2D *hallnodes_ECAL_vs_HCAL = new TH2D("hallnodes_ECAL_vs_HCAL","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5);

  TH1D *hshouldhit_vs_threshold_ECAL = new TH1D("hshouldhit_vs_threshold_ECAL","",30,0.025,1.525);
  TH1D *hefficiency_vs_threshold_ECAL = new TH1D("hefficiency_vs_threshold_ECAL","",30,0.025,1.525);
  TH1D *hshouldhit_vs_threshold_ECAL_FTcut = new TH1D("hshouldhit_vs_threshold_ECAL_FTcut","",30,0.025,1.525);
  TH1D *hefficiency_vs_threshold_ECAL_FTcut = new TH1D("hefficiency_vs_threshold_ECAL_FTcut","",30,0.025,1.525);
  
  TH1D *hefficiency_vs_threshold_HCAL_FTcut = new TH1D("hefficiency_vs_threshold_HCAL_FTcut","",30,0.025,1.525);
  TH1D *hefficiency_vs_threshold_HCAL_FPP1cut = new TH1D("hefficiency_vs_threshold_HCAL_FPP1cut","",30,0.025,1.525);
  TH1D *hefficiency_vs_threshold_HCAL_FPP2cut = new TH1D("hefficiency_vs_threshold_HCAL_FPP2cut","",30,0.025,1.525);

  TH2D *hefficiency_vs_threshold_ECAL_HCAL_coincidence_FTcut = new TH2D("hefficiency_vs_threshold_ECAL_HCAL_coincidence_FTcut","",30,0.025,1.525,30,0.025,1.525);
  TH2D *hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FTcut = new TH2D("hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FTcut","",30,0.025,1.525,30,0.025,1.525);

  TH2D *hefficiency_vs_threshold_ECAL_HCAL_coincidence_FPP1cut = new TH2D("hefficiency_vs_threshold_ECAL_HCAL_coincidence_FPP1cut","",30,0.025,1.525,30,0.025,1.525);
  TH2D *hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FPP1cut = new TH2D("hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FPP1cut","",30,0.025,1.525,30,0.025,1.525);

  TH2D *hefficiency_vs_threshold_ECAL_HCAL_coincidence_FPP2cut = new TH2D("hefficiency_vs_threshold_ECAL_HCAL_coincidence_FPP2cut","",30,0.025,1.525,30,0.025,1.525);
  TH2D *hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FPP2cut = new TH2D("hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FPP2cut","",30,0.025,1.525,30,0.025,1.525);
  
  TH1D *hshouldhit_HCAL_FTcut = new TH1D("hshouldhit_HCAL_FTcut","",30,0.025,1.525);
  TH1D *hshouldhit_HCAL_FPP1cut = new TH1D("hshouldhit_HCAL_FPP1cut","",30,0.025,1.525);
  TH1D *hshouldhit_HCAL_FPP2cut = new TH1D("hshouldhit_HCAL_FPP2cut","",30,0.025,1.525);
  
  TH2D *hnphe_vs_sum_edep_ECAL = new TH2D("hnphe_vs_sum_edep_ECAL","",125,0.0,5.0,125,0.0,5000.0 );
  TH2D *hnphe_vs_sum_edep_HCAL = new TH2D("hnphe_vs_sum_edep_HCAL","",125,0.0,0.75,125,0.0,2500.0 );

  TH1D *hthetaFPP1 = new TH1D("hthetaFPP1","",120,0.0,12.0);
  TH1D *hthetaFPP2 = new TH1D("hthetaFPP2","",120,0.0,12.0);

  TH1D *hthetaFPP1_cointrig = new TH1D("hthetaFPP1_cointrig","",120,0.0,12.0);
  TH1D *hthetaFPP2_cointrig = new TH1D("hthetaFPP2_cointrig","",120,0.0,12.0);
  
  TH1D *hthetaFPP1_HCALtrig = new TH1D("hthetaFPP1_HCALtrig","",120,0.0,12.0);
  TH1D *hthetaFPP2_HCALtrig = new TH1D("hthetaFPP2_HCALtrig","",120,0.0,12.0);
  
  TH1D *hshouldhit_vs_Q2_ECAL_FTcut = new TH1D("hshouldhit_vs_Q2_ECAL_FTcut","",100,8.0,16.0);
  TH1D *hefficiency_vs_Q2_ECAL_FTcut = new TH1D("hefficiency_vs_Q2_ECAL_FTcut","",100,8.0,16.0);

  TH1D *HCAL_hit_edep = new TH1D("HCAL_hit_edep","",200,0.0,1.0); //GeV
  TH1D *HCAL_hit_num_phe = new TH1D("HCAL_hit_num_phe","",200,0.0,5000.0);

  TH1D *HCAL_maxhit_edep = new TH1D("HCAL_maxhit_edep","",200,0.0,1.0); //GeV
  TH1D *HCAL_maxhit_num_phe = new TH1D("HCAL_maxhit_num_phe","",200,0.0,5000.0);
  
  
  double Ibeam = 75.0e-6; //Amps
  double Ltarget = 40.0; //cm
  double e = 1.602e-19; //electron charge;
  double rho_target = 0.072; //g/cm^3
  double N_A = 6.022e23; //atoms/mol:
  double Mmol_H = 1.008; //g/mol
  double Lumi = rho_target * Ltarget * N_A / Mmol_H * Ibeam/e; //~ 8e38;
  
  TRandom3 num(0);

  cout << "Entering event loop " << endl;
  
  long nevent=0;
  for( nevent=0; nevent<C->GetEntries(); ++nevent ){
    T->GetEntry(nevent);

    double weight;
    //cross section is given in mb: 1 mb = 1e-3 * 1e-24 = 1e-27 cm^2
    // if (pythia6flag != 0 ){
    //   weight = Lumi * T->primaries_Sigma * 1.0e-27/ double(ngen); //luminosity times cross section / number of events generated.
    // } else {
    weight = T->ev_rate / double(nfiles);

    TString fnametemp = C->GetFile()->GetName();

    if( bad_file_list.find( fnametemp ) == bad_file_list.end() ){

      double R = ECALdist_file[fnametemp];
      double thetacal = ECALtheta_file[fnametemp];

      // cout << "(R, thetacal)=(" << R << ", " << thetacal*57.3 << endl;
    //if( Q2cut == 0 || (T->ev_Q2 >= 10.5 && T->ev_Q2 <= 14.0) ){
    
      bool FTtrack = false;
      int itrack_FT=-1;
      for( int itrack=0; itrack<T->Harm_FT_Track_ntracks; itrack++ ){
	if( (*(T->Harm_FT_Track_MID))[itrack] == 0 &&
	    (*(T->Harm_FT_Track_PID))[itrack] == 2212 ){ //primary elastically scattered proton track in FT:
	  FTtrack = true;
	  itrack_FT = itrack;
	}
      }

      TVector3 nhat_FT, nhat_FPP1, nhat_FPP2;
      if( FTtrack ){
	nhat_FT.SetXYZ( (*(T->Harm_FT_Track_Xp))[itrack_FT],
			(*(T->Harm_FT_Track_Yp))[itrack_FT],
			1.0 );
	nhat_FT = nhat_FT.Unit();
      }

      double thetaFPP1, thetaFPP2, pFPP1, pFPP2;
      bool FPP1track = false, FPP2track = false;
      //    if( FTtrack )
      if( T->Harm_FPP1_Track_ntracks > 0 && FTtrack ){
	double thetamin;
	for( int itrack=0; itrack<T->Harm_FPP1_Track_ntracks; itrack++ ){
	  if( (*(T->Harm_FPP1_Track_MID))[itrack] == 0 ){
	    nhat_FPP1.SetXYZ( (*(T->Harm_FPP1_Track_Xp))[itrack],
			      (*(T->Harm_FPP1_Track_Yp))[itrack],
			      1.0 );
	    nhat_FPP1 = nhat_FPP1.Unit();

	    thetaFPP1 = acos( nhat_FPP1.Dot( nhat_FT ) );
	  
	    pFPP1 = (*(T->Harm_FPP1_Track_P))[itrack];

	  
	    FPP1track = true;
	  }   
	}
	if( FPP1track ) hthetaFPP1->Fill(thetaFPP1*180.0/PI,weight);
      }

      if( T->Harm_FPP2_Track_ntracks > 0 && FTtrack && FPP1track){
	for( int itrack=0; itrack<T->Harm_FPP2_Track_ntracks; itrack++ ){
	  if( (*(T->Harm_FPP2_Track_MID))[itrack] == 0 ){
	    nhat_FPP2.SetXYZ( (*(T->Harm_FPP2_Track_Xp))[itrack],
			      (*(T->Harm_FPP2_Track_Yp))[itrack],
			      1.0 );
	    nhat_FPP2 = nhat_FPP2.Unit();
	    thetaFPP2 = acos( nhat_FPP2.Dot( nhat_FPP1 ) );
	    pFPP2 = (*(T->Harm_FPP2_Track_P))[itrack];
	    //FPP2track = thetaFPP2 < 24.0*PI/180.0 && pFPP2/T->ev_np > 0.5;
	    FPP2track = thetaFPP2 < 12.0*PI/180.0 && pFPP2/T->ev_np > 0.5;
	    if( FPP2track ) hthetaFPP2->Fill(thetaFPP2*180.0/PI,weight);
	  }
	}
      }
      
      double nu = T->ev_Q2 / 2.0 / 0.938272;
      double pp_elastic = sqrt(pow(nu,2)+2.0*.938272*nu);
       
      //}

      // double R = T->gen_dbb;
      // double thetacal = T->gen_thbb;
    
      //ECAL is on beam left:
      TVector3 nhat_e( sin( T->ev_th )*cos( T->ev_ph ), sin(T->ev_th)*sin(T->ev_ph), cos(T->ev_th) );
      TVector3 vertex( T->ev_vx, T->ev_vy, T->ev_vz );
      TVector3 ecal_z( sin(thetacal), 0, cos(thetacal) );
      TVector3 ecal_y(0,1,0);
      TVector3 ecal_x = (ecal_y.Cross(ecal_z)).Unit();

      TVector3 Rcalo = R * ecal_z;
      //ecal_z dot (vertex + s * nhat_e - Rcalo ) = 0;
      double s = (Rcalo - vertex).Dot( ecal_z )/ ( nhat_e.Dot( ecal_z ) );

      TVector3 pos_calo = vertex + s * nhat_e;

      double xcalo = (pos_calo - Rcalo).Dot( ecal_x );
      double ycalo = (pos_calo - Rcalo).Dot( ecal_y ); //
    
      if( (nevent+1) % 1000 == 0 ){ cout << "Event number " << nevent+1 << ", event weight = " << weight << endl; }
    
      map<int,double> node_sums; //initialize all node sums to zero:
      for( set<int>::iterator inode = list_of_nodes_ecal.begin(); inode != list_of_nodes_ecal.end(); ++inode ){
	node_sums[ *inode ] = 0.0;
      }

      bool should_hit_ECAL = false;

      //cout << "(xcalo,ycalo,should_hit_ECAL)=("
      
      if( ycalo >= ycellmin && ycalo <= ycellmax  ){
	//make an initial guess at which row: (row runs from 1 to N):
	int closest_row = int( (ycalo - ycellmin)/(ycellmax-ycellmin)*double(ycell_rows.size()) );

	map<int,double>::iterator rowguess = ycell_rows.find( closest_row );

	if( rowguess != ycell_rows.end() ){
	//since map is sorted on row, we can grab rowmin and rowmax:
	  int rowmin = (ycell_rows.begin())->first;
	  int rowmax = (ycell_rows.rbegin())->first;
	  
	  while( closest_row < rowmax && fabs( ycalo - ycell_rows[closest_row+1] ) < fabs( ycalo - ycell_rows[closest_row] ) ){ closest_row++; }
	  while( closest_row > 0 && fabs( ycalo - ycell_rows[closest_row-1] ) < fabs( ycalo - ycell_rows[closest_row] ) ){ closest_row--; }

	  //cout << "(rowmin,rowmax,closest row) = (" << rowmin << ", " << rowmax << ", " << closest_row << ")" << endl;

	  //cout << "xcellmin, xcellmax = " << xcellmin_rows[closest_row] << ", " << xcellmax_rows[closest_row] << ")" << endl;
	  if( xcalo >= xcellmin_rows[closest_row]  &&
	      xcalo <= xcellmax_rows[closest_row]  ){

	    //cout << "should hit ECAL SHOULD be true here" << endl;
	    should_hit_ECAL = true;
	  }
	}
      }

      // cout << "(xcalo,ycalo,should_hit_ECAL) = ("
      // 	   << xcalo << ", " << ycalo << ", " << should_hit_ECAL << ")" << endl;

      int nphe = 0;
    
      if( pheflag == 0 ){
	for( int ihit = 0; ihit<T->Earm_ECalTF1_hit_nhits; ihit++ ){
	  int rowhit = ( *(T->Earm_ECalTF1_hit_row))[ihit]+1;
	  int colhit = ( *(T->Earm_ECalTF1_hit_col))[ihit]+1;
	  std::pair<int,int> rowcolhit( rowhit,colhit );
	
	  //int cellhit = cell_rowcol_ecal[rowcolhit];
	  int cellhit = (*(T->Earm_ECalTF1_hit_cell))[ihit];
	  
	  //int trigger_group = nodes_cells_ecal[cellhit];
	
	  double edep = (*(T->Earm_ECalTF1_hit_sumedep))[ihit];

	  double mean = 752.2*edep;
	  double sigma = 52.0*sqrt(edep) + 20.76*edep;

	  nphe = TMath::Max(0,TMath::Nint(num.Gaus(mean,sigma)));
	
	  for( set<int>::iterator inode = nodes_cells_ecal[cellhit].begin(); inode != nodes_cells_ecal[cellhit].end(); ++inode ){
	    node_sums[ *inode ] += double(nphe);
	  }
	
	}
      } else {
	for( int ihit = 0; ihit<T->Earm_ECAL_hit_nhits; ihit++){
	  int rowhit = ( *(T->Earm_ECAL_hit_row))[ihit]+1;
	  int colhit = ( *(T->Earm_ECAL_hit_col))[ihit]+1;
	  std::pair<int,int> rowcolhit( rowhit,colhit );
	
	  //int cellhit = cell_rowcol_ecal[rowcolhit];

	  int cellhit = (*(T->Earm_ECAL_hit_PMT))[ihit];
	  
	  //int trigger_group = nodes_cells_ecal[cellhit];
	
	  //	double edep = (*(T->Earm_ECalTF1_hit_sumedep))[ihit];

	  int nphe = (*(T->Earm_ECAL_hit_NumPhotoelectrons))[ihit];
	
	  for( set<int>::iterator inode = nodes_cells_ecal[cellhit].begin(); inode != nodes_cells_ecal[cellhit].end(); ++inode ){
	    node_sums[ *inode ] += double(nphe);
	  }
	  for( int jhit=0; jhit<T->Earm_ECalTF1_hit_nhits; jhit++ ){
	    if( (*(T->Earm_ECalTF1_hit_cell))[jhit] == cellhit && 
		fabs( (*(T->Earm_ECAL_hit_Time_avg))[ihit]-(*(T->Earm_ECalTF1_hit_tavg))[jhit]-2.5)<=10.0 ){
	      hnphe_vs_sum_edep_ECAL->Fill( (*(T->Earm_ECalTF1_hit_sumedep))[jhit], nphe );
	    }   
	  }
	}
      
	//node_sums[ trigger_group ] += double(nphe);
      }

      vector<int> trigger_nodes_fired(hefficiency_vs_threshold_ECAL->GetNbinsX());
      for( int ithr=0; ithr<hefficiency_vs_threshold_ECAL->GetNbinsX(); ithr++ ){
	trigger_nodes_fired[ithr] = 0;
      }

      int maxnode_ECAL=-1;
      int maxnode_HCAL=-1;
      double maxsum_ECAL = 0.0;
      double maxsum_HCAL = 0.0;

      bool ECALtrig_nominal = false;
    
      int nominal_threshold_bin_HCAL = hefficiency_vs_threshold_HCAL_FTcut->FindBin(nominal_threshold_HCAL);
      int nominal_threshold_bin_ECAL = hefficiency_vs_threshold_ECAL->FindBin(nominal_threshold_ECAL);
    
      for( set<int>::iterator inode = list_of_nodes_ecal.begin(); inode != list_of_nodes_ecal.end(); ++inode ){
	for( int bin=1; bin<=hefficiency_vs_threshold_ECAL->GetNbinsX(); bin++ ){
	  if( node_sums[*inode]/logic_mean_ecal[*inode] > hefficiency_vs_threshold_ECAL->GetBinCenter(bin) ){
	    //cout << "node above threshold, nphe, peak position = " << node_sums[*inode] << ", " << logic_mean_ecal[*inode] << endl;
	    trigger_nodes_fired[bin-1]++;
	    if( bin == nominal_threshold_bin_ECAL ) ECALtrig_nominal = true;
	  }
	
	}
	if( node_sums[*inode] > maxsum_ECAL ) {
	  maxsum_ECAL = node_sums[*inode];
	  maxnode_ECAL = *inode;
	}
      
	if( node_sums[*inode] > 0.0 ) hnphesum_vs_node_ECAL_all->Fill( *inode, node_sums[*inode], weight );
      }

      if( should_hit_ECAL ){
	for( int ithr=0; ithr<hefficiency_vs_threshold_ECAL->GetNbinsX(); ithr++ ){
	  hshouldhit_vs_threshold_ECAL->Fill( hefficiency_vs_threshold_ECAL->GetBinCenter(ithr+1), weight );
	  if( trigger_nodes_fired[ithr] > 0 ){
	    hefficiency_vs_threshold_ECAL->Fill( hefficiency_vs_threshold_ECAL->GetBinCenter(ithr+1), weight );
	  }
	  if( FTtrack ){
	    hshouldhit_vs_threshold_ECAL_FTcut->Fill( hefficiency_vs_threshold_ECAL->GetBinCenter(ithr+1), weight );
	    if( trigger_nodes_fired[ithr] > 0 ){
	      hefficiency_vs_threshold_ECAL_FTcut->Fill( hefficiency_vs_threshold_ECAL->GetBinCenter(ithr+1), weight );
	    }
	  }
	}
      }

      if( FTtrack ){
	hshouldhit_vs_Q2_ECAL_FTcut->Fill( T->ev_Q2, weight );
	if( ECALtrig_nominal ){
	  hefficiency_vs_Q2_ECAL_FTcut->Fill( T->ev_Q2, weight );
	}
      }
    
      map<int,double> node_sums_hcal;
      for( set<int>::iterator inode = list_of_nodes_hcal.begin(); inode != list_of_nodes_hcal.end(); ++inode ){
	node_sums_hcal[*inode] = 0.0;
      }

      //int nphe = 0;
    
      if( pheflag == 0 ){

	double maxedep = 0.0;
	double maxphe = 0.0;
	
	for( int ihit=0; ihit<T->Harm_HCalScint_hit_nhits; ihit++ ){
	  int rowhit = (*(T->Harm_HCalScint_hit_row))[ihit]+1;
	  int colhit = (*(T->Harm_HCalScint_hit_col))[ihit]+1;
	  std::pair<int,int> rowcolhit(rowhit,colhit);
	  int cellhit = cell_rowcol_hcal[rowcolhit];
	  //int trigger_group = nodes_cells_hcal[cellhit];
	  double edep = (*(T->Harm_HCalScint_hit_sumedep))[ihit];
	  //nphe = num.Poisson( phe_per_GeV_HCAL * edep );
	  double mean = 2981.0*edep;
	  double sigma = 69.54*sqrt(edep) + 155.3*edep;

	  nphe = TMath::Max(0,TMath::Nint(num.Gaus(mean,sigma)));
	
	  //cout << "HCAL hit " << ihit+1 << " node, edep, nphe = " << trigger_group << ", " << edep << ", " << nphe << endl;
	  //node_sums_hcal[trigger_group] += double(nphe);
	  for( set<int>::iterator inode = nodes_cells_hcal[cellhit].begin(); inode != nodes_cells_hcal[cellhit].end(); ++inode ){
	  
	    node_sums_hcal[*inode] += double(nphe);
	  
	  }

	  HCAL_hit_edep->Fill( edep );
	  HCAL_hit_num_phe->Fill(nphe);
	  if( edep > maxedep || ihit == 0 ){
	    maxedep = edep;
	  }
	  if( nphe > maxphe || ihit == 0 ){
	    maxphe = nphe;
	  }
	}

	HCAL_maxhit_edep->Fill(maxedep);
	HCAL_maxhit_num_phe->Fill(maxphe);
	
      } else {
	double maxedep = 0.0;
	double maxphe = 0.0;
	
	for( int jhit=0; jhit<T->Harm_HCal_hit_nhits; jhit++ ){
	  int rowhit = (*(T->Harm_HCal_hit_row))[jhit]+1;
	  int colhit = (*(T->Harm_HCal_hit_col))[jhit]+1;
	  std::pair<int,int> rowcolhit(rowhit,colhit);
	  int cellhit = cell_rowcol_hcal[rowcolhit];
	  nphe = (*(T->Harm_HCal_hit_NumPhotoelectrons))[jhit];
	  for( set<int>::iterator inode = nodes_cells_hcal[cellhit].begin(); inode != nodes_cells_hcal[cellhit].end(); ++inode ){
	  
	    node_sums_hcal[*inode] += double(nphe);
	  
	  }

	  HCAL_hit_num_phe->Fill( nphe );

	  if( jhit == 0 || nphe > maxphe ){
	    maxphe = nphe;
	  }
	  
	  for( int khit=0; khit<T->Harm_HCalScint_hit_nhits; khit++ ){
	    if( (*(T->Harm_HCalScint_hit_row))[khit]+1 == rowhit &&
		(*(T->Harm_HCalScint_hit_col))[khit]+1 == colhit &&
		fabs( (*(T->Harm_HCal_hit_Time_avg))[jhit]-(*(T->Harm_HCalScint_hit_tavg))[khit] - 8.6 )<=15.0 ){
	      hnphe_vs_sum_edep_HCAL->Fill( (*(T->Harm_HCalScint_hit_sumedep))[khit], nphe );
	      HCAL_hit_edep->Fill( (*(T->Harm_HCalScint_hit_sumedep))[khit] );
	      if( jhit == 0 || (*(T->Harm_HCalScint_hit_sumedep))[khit] > maxedep ){
		maxedep = (*(T->Harm_HCalScint_hit_sumedep))[khit];
	      }
	    }
	  }
	
	}

	HCAL_maxhit_edep->Fill(maxedep);
	HCAL_maxhit_num_phe->Fill(maxphe);
	
      }
    
      vector<int> trigger_nodes_fired_hcal(hefficiency_vs_threshold_HCAL_FTcut->GetNbinsX());
      for( int ithr=0; ithr<hefficiency_vs_threshold_HCAL_FTcut->GetNbinsX(); ithr++ ){
	trigger_nodes_fired_hcal[ithr] = 0;
      }

      vector<int> coin_trigger_fired( hefficiency_vs_threshold_HCAL_FTcut->GetNbinsX()*hefficiency_vs_threshold_ECAL->GetNbinsX() );
      for( int ithr=0; ithr<coin_trigger_fired.size(); ithr++ ){
	coin_trigger_fired[ithr] = 0;
      }

    

      bool cointrig_nominal_threshold = false;
    
      for( set<int>::iterator inode = list_of_nodes_hcal.begin(); inode != list_of_nodes_hcal.end(); ++inode ){
	for( int bin=1; bin<=hefficiency_vs_threshold_HCAL_FTcut->GetNbinsX(); bin++ ){
	  if( node_sums_hcal[*inode]/logic_mean_hcal[*inode] > hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(bin) ){ //this HCAL sum fired:
	    trigger_nodes_fired_hcal[bin-1]++;
	    for( set<int>::iterator enode = list_of_nodes_ecal.begin(); enode != list_of_nodes_ecal.end(); ++enode ){
	      if( ECAL_nodes_HCAL[*inode].find(*enode) != ECAL_nodes_HCAL[*inode].end() ){ //Check associated ECAL trigger sums:
		for( int ebin=1; ebin<=hefficiency_vs_threshold_ECAL->GetNbinsX(); ebin++ ){ //check ECAL sums:
		  if( node_sums[ *enode ]/logic_mean_ecal[*enode] > hefficiency_vs_threshold_ECAL->GetBinCenter(ebin) ){ //this ECAL sum fired:
		    coin_trigger_fired[ (ebin-1) + (bin-1)*hefficiency_vs_threshold_ECAL->GetNbinsX() ]++;
		    if( ebin == nominal_threshold_bin_ECAL && bin == nominal_threshold_bin_HCAL ){
		      cointrig_nominal_threshold = true;
		    }
		  }
		}
	      } 
	    }
	  }
	}
	if( node_sums_hcal[*inode] > maxsum_HCAL ) {
	  maxsum_HCAL = node_sums_hcal[*inode];
	  maxnode_HCAL = *inode;
	}

	hnphesum_vs_node_HCAL_all->Fill( *inode, node_sums_hcal[*inode], weight );
	if( FTtrack ){
	  hnphesum_vs_node_HCAL_FTcut->Fill( *inode, node_sums_hcal[*inode], weight );
	  if( FPP1track ) hnphesum_vs_node_HCAL_FPP1cut->Fill( *inode, node_sums_hcal[*inode], weight );
	  if( FPP2track ) hnphesum_vs_node_HCAL_FPP2cut->Fill( *inode, node_sums_hcal[*inode], weight );
	  if( FPP1track && FPP2track ) hnphesum_vs_node_HCAL_FPPbothcut->Fill( *inode, node_sums_hcal[*inode], weight );
	  if( FPP1track || FPP2track ) hnphesum_vs_node_HCAL_FPPeithercut->Fill( *inode, node_sums_hcal[*inode], weight );
	}
      }

      if( cointrig_nominal_threshold ){
	if( FPP1track ) hthetaFPP1_cointrig->Fill(thetaFPP1*180.0/PI,weight);
	if( FPP2track ) hthetaFPP2_cointrig->Fill(thetaFPP2*180.0/PI,weight);
      }
    
      for( int bin=1; bin<=hefficiency_vs_threshold_HCAL_FTcut->GetNbinsX(); bin++ ){
	if( FTtrack ) hshouldhit_HCAL_FTcut->Fill( hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(bin), weight );
	if( FTtrack && FPP1track ) hshouldhit_HCAL_FPP1cut->Fill( hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(bin), weight );
	if( FTtrack && FPP2track ) hshouldhit_HCAL_FPP2cut->Fill( hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(bin), weight );
	if( trigger_nodes_fired_hcal[bin-1] > 0 ){
	  if( FTtrack ){
	    hefficiency_vs_threshold_HCAL_FTcut->Fill( hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(bin), weight );
	    if( FPP1track ) hefficiency_vs_threshold_HCAL_FPP1cut->Fill( hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(bin), weight );
	    if( FPP2track ) hefficiency_vs_threshold_HCAL_FPP2cut->Fill( hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(bin), weight );
	  }
	}
      }

      for( int ithr=0; ithr<coin_trigger_fired.size(); ithr++ ){
	int bin_e = ithr%(hefficiency_vs_threshold_ECAL->GetNbinsX())+1;
	int bin_h = ithr/(hefficiency_vs_threshold_HCAL_FTcut->GetNbinsX())+1;
	double thr_e = hefficiency_vs_threshold_ECAL_HCAL_coincidence_FTcut->GetYaxis()->GetBinCenter(bin_e);
	double thr_h = hefficiency_vs_threshold_ECAL_HCAL_coincidence_FTcut->GetXaxis()->GetBinCenter(bin_h);
	if( FTtrack ) hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FTcut->Fill( thr_h, thr_e, weight );
	if( FTtrack && FPP1track ) hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FPP1cut->Fill( thr_h, thr_e, weight );
	if( FTtrack && FPP2track ) hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FPP2cut->Fill( thr_h, thr_e, weight );
	if( coin_trigger_fired[ithr] > 0 ){
	  if( FTtrack ) hefficiency_vs_threshold_ECAL_HCAL_coincidence_FTcut->Fill( thr_h, thr_e, weight );
	  if( FTtrack && FPP1track ) hefficiency_vs_threshold_ECAL_HCAL_coincidence_FPP1cut->Fill( thr_h, thr_e, weight );
	  if( FTtrack && FPP2track ) hefficiency_vs_threshold_ECAL_HCAL_coincidence_FPP2cut->Fill( thr_h, thr_e, weight );
	}
      }
    
      hnphesum_vs_node_HCALmax_all->Fill( maxnode_HCAL, node_sums_hcal[maxnode_HCAL], weight );
      if( FTtrack ){
	hnphesum_vs_node_HCALmax_FTcut->Fill( maxnode_HCAL, node_sums_hcal[maxnode_HCAL], weight );
	if( FPP1track ) hnphesum_vs_node_HCALmax_FPP1cut->Fill( maxnode_HCAL, node_sums_hcal[maxnode_HCAL], weight );
	if( FPP2track ) hnphesum_vs_node_HCALmax_FPP2cut->Fill( maxnode_HCAL, node_sums_hcal[maxnode_HCAL], weight );
	if( FPP1track && FPP2track ) hnphesum_vs_node_HCALmax_FPPbothcut->Fill( maxnode_HCAL, node_sums_hcal[maxnode_HCAL], weight );
	if( FPP1track || FPP2track ) hnphesum_vs_node_HCALmax_FPPeithercut->Fill( maxnode_HCAL, node_sums_hcal[maxnode_HCAL], weight );
      }
    
      // for( int ithr=0; ithr<hefficiency_vs_threshold_HCAL_FTcut->GetNbinsX(); ithr++ ){
      //   if( trigger_nodes_fired_hcal[ithr] > 0 ) hefficiency_vs_threshold_HCAL_FTcut->Fill( hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(ithr+1),weight );
      //   for( int jthr=0; jthr<hefficiency_vs_threshold_ECAL->GetNbinsX(); jthr++ ){
      // 	if( trigger_nodes_fired[jthr] > 0 && trigger_nodes_fired_hcal[ithr] > 0 ){
      // 	  //htrue_coincidence_rate_vs_threshold_ECAL_HCAL->Fill( hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(ithr+1),hefficiency_vs_threshold_ECAL->GetBinCenter(jthr+1),weight );
      // 	}
      //   }
      // }

      // for( set<int>::iterator inode = list_of_nodes_ecal.begin(); inode != list_of_nodes_ecal.end(); ++inode ){
      //   for( set<int>::iterator jnode = list_of_nodes_hcal.begin(); jnode != list_of_nodes_hcal.end(); ++jnode ){
      // 	//Fill the correlation histogram for all true coincidence events for which ECAL and HCAL node are both above threshold:
      // 	if( node_sums[*inode] >= nominal_threshold_ECAL*logic_mean_ecal[*inode] && node_sums_hcal[*jnode] >= nominal_threshold_HCAL*logic_mean_hcal[*jnode] ){
      // 	  hallnodes_ECAL_vs_HCAL->Fill( *jnode, *inode, weight );
      // 	}
      //   }
      // }
      //if( maxsum_ECAL >= nominal_threshold_ECAL*logic_mean_ecal[maxnode_ECAL] && maxsum_HCAL >= nominal_threshold_HCAL*logic_mean_hcal[maxnode_HCAL] ){
      if( FTtrack && FPP2track) {
	hmaxnode_ECAL_vs_HCAL->Fill( maxnode_HCAL, maxnode_ECAL, weight );
      }
      //}
    }
  }

  hefficiency_vs_threshold_HCAL_FTcut->Divide( hshouldhit_HCAL_FTcut );
  hefficiency_vs_threshold_HCAL_FPP1cut->Divide( hshouldhit_HCAL_FPP1cut );
  hefficiency_vs_threshold_HCAL_FPP2cut->Divide( hshouldhit_HCAL_FPP2cut );

  hefficiency_vs_threshold_ECAL_HCAL_coincidence_FTcut->Divide( hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FTcut );
  hefficiency_vs_threshold_ECAL_HCAL_coincidence_FPP1cut->Divide( hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FPP1cut );
  hefficiency_vs_threshold_ECAL_HCAL_coincidence_FPP2cut->Divide( hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FPP2cut );
  
  // TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  
  // c1->Divide(2,1);

  // c1->cd(1)->SetLogy();
  // hefficiency_vs_threshold_ECAL->SetMarkerStyle(20);
  // hefficiency_vs_threshold_ECAL->Draw();
  
  // c1->cd(2)->SetLogy();
  // hefficiency_vs_threshold_HCAL_FTcut->SetMarkerStyle(20);
  // hefficiency_vs_threshold_HCAL_FTcut->Draw();

  hefficiency_vs_threshold_ECAL->Divide( hshouldhit_vs_threshold_ECAL );
  hefficiency_vs_threshold_ECAL_FTcut->Divide(hshouldhit_vs_threshold_ECAL_FTcut );
  
  hefficiency_vs_threshold_HCAL_FTcut->SetMarkerStyle(20);
  hefficiency_vs_threshold_ECAL->SetMarkerStyle(20);

  hefficiency_vs_Q2_ECAL_FTcut->Divide(hshouldhit_vs_Q2_ECAL_FTcut);
  
  fout->Write();
  fout->Close();
}
