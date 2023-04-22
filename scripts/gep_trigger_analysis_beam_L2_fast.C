#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "gep_tree_beam.C"
#include "G4SBSRunData.hh"
#include "TRandom3.h"
#include "TChainElement.h"
#include "TObjArray.h"
#include "TString.h"
#include "TObjString.h"
//#include "TIter.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
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

void gep_trigger_analysis_beam_L2_fast( const char *rootfilename, const char *outputfilename="gep_trigrates_temp.root", double ECAL_threshold=0.8, double HCAL_threshold=0.5, const char *logicfilename_ecal="database/GEP_ECAL_L2sums.txt", const char *thresholdfilename_ecal="database/ecal_trigger_thresholds_12GeV2.txt", const char *thresholdfilename_hcal="database/hcal_trigger_thresholds_12GeV2.txt", int pheflag=0, const char *assocfilename="database/ECAL_HCAL_L2_default.txt" ){
  
  //double thetacal = thetacaldeg*PI/180.0;
  
  //The difference between this analysis and the "slow" one is that we only look at coincidence rates (real and accidental) for one choice of threshold:
  
  TFile *fout = new TFile(outputfilename,"RECREATE");
  TChain *C = new TChain("T");
  C->Add(rootfilename);

  C->SetBranchStatus("*",0);
  C->SetBranchStatus("Harm.HCalScint.hit.nhits",1);
  C->SetBranchStatus("Harm.HCalScint.hit.row",1);
  C->SetBranchStatus("Harm.HCalScint.hit.col",1);
  C->SetBranchStatus("Harm.HCalScint.hit.cell",1);
  C->SetBranchStatus("Harm.HCalScint.hit.sumedep",1);

  C->SetBranchStatus("Earm.ECalTF1.hit.nhits",1);
  C->SetBranchStatus("Earm.ECalTF1.hit.row",1);
  C->SetBranchStatus("Earm.ECalTF1.hit.col",1);
  C->SetBranchStatus("Earm.ECalTF1.hit.cell",1);
  C->SetBranchStatus("Earm.ECalTF1.hit.sumedep",1);
  //T->SetBranchStatus("Harm.HCalScint.hit.nhits",1);
  //T->SetBranchStatus("Harm.HCalScint.hit.nhits",1);

  //  T->SetBranchStatus("*",0);
  C->SetBranchStatus("Harm.HCal.hit.nhits",1);
  C->SetBranchStatus("Harm.HCal.hit.row",1);
  C->SetBranchStatus("Harm.HCal.hit.col",1);
  C->SetBranchStatus("Harm.HCal.hit.PMT",1);
  C->SetBranchStatus("Harm.HCal.hit.NumPhotoElectrons",1);

  C->SetBranchStatus("Earm.ECAL.hit.nhits",1);
  C->SetBranchStatus("Earm.ECAL.hit.row",1);
  C->SetBranchStatus("Earm.ECAL.hit.col",1);
  C->SetBranchStatus("Earm.ECAL.hit.PMT",1);
  C->SetBranchStatus("Earm.ECAL.hit.NumPhotoElectrons",1);

  // C->SetBranchStatus("primaries.Sigma",1);
  // C->SetBranchStatus("primaries.Q2",1);
  // C->SetBranchStatus("primaries.W2",1);
  // C->SetBranchStatus("primaries.y",1);
  // C->SetBranchStatus("primaries.theta_e",1);
  // C->SetBranchStatus("primaries.xbj",1);
  // //C->SetBranchStatus("primaries.phi_e",1);
  // C->SetBranchStatus("Primaries.genflag",1);
  // C->SetBranchStatus("Primaries.PID",1);
  // C->SetBranchStatus("Primaries.P",1);
  // C->SetBranchStatus("Primaries.theta",1);
  
  gep_tree_beam *T = new gep_tree_beam( C );

  
  

  G4SBSRunData *rd;

  long ngen = 0;
  int nfiles = 0;
  
  TObjArray *FileList = C->GetListOfFiles();
  TIter next(FileList);

  TChainElement *chEl = 0;

  set<TString> bad_file_list;

  map<TString,double> Rcal_file;
  map<TString,double> thetacal_file;
  map<TString,double> LumiFile;

  double Lumi = 6.0e38;

  //double sigma_default=2.2769680e-30; //cm2, 11 GeV beam on proton target, minimum bias
  
  while( (chEl=(TChainElement*)next() )){
    TFile newfile(chEl->GetTitle());

    if( !newfile.TestBit(TFile::kRecovered) ){
    
      newfile.GetObject("run_data",rd);
      TTree *Ttemp;
      newfile.GetObject("T",Ttemp);
      double sigtemp;
      if( rd && Ttemp ){
	ngen += rd->fNtries;

	Rcal_file[chEl->GetTitle()] = rd->fBBdist;
	thetacal_file[chEl->GetTitle()] = rd->fBBtheta;
	LumiFile[chEl->GetTitle()] = rd->fLuminosity;

	Lumi = rd->fLuminosity;
      
	cout << "file " << chEl->GetTitle() << ", ngen=" << rd->fNtries << endl;
      
	// Ttemp->SetBranchAddress("primaries.Sigma",&sigtemp);
      
	// if( Ttemp->GetEntries() > 0 ){
  
	//   Ttemp->GetEntry(0);

	//   if ( !isnan( sigtemp ) ){
	//     sigma_default = sigtemp;
	//     cout << "first event cross section = " << sigtemp << " cm^2" << endl;
	//   } else {
	//     cout << "bad cross section!!!!" << endl;
	//   }
	// }

	//	Ttemp->ResetBranchAddresses();

	//      Ttemp->Delete();
	//newfile.Close();
	nfiles++;
      } else {
	bad_file_list.insert( chEl->GetTitle());
      }
    } else {
      bad_file_list.insert( chEl->GetTitle());
      //newfile.Close();
    }

    newfile.Close();
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
  map<int,double> cellsize_rows;
  map<int,double> xcellmin_rows;
  map<int,double> xcellmax_rows;
  
  int minrow=1000,maxrow=-1;
  
  set<int> rows_ecal;
  map<int,set<int> > columns_rows_ecal;
  
  map<int,double> elastic_peak_new_ecal;
  map<int,double> sigma_new_ecal;
  map<int,double> threshold_new_ecal;
  
  ifstream logicfile_ecal(logicfilename_ecal);
  //ifstream thresholdfile(thresholdfilename);

  TString currentline;
  
  int current_node = 1;

  bool first_cell = true;

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
	  logic_mean_ecal[nodenumber] = 3.0;
	  logic_sigma_ecal[nodenumber] = 0.06*3.0;
	  threshold_ecal[nodenumber] = 0.8;
	  
	  nodes_cells_ecal[ cell ].insert( nodenumber );
	}
      }
    }
  }
  
  // while( currentline.ReadLine( logicfile_ecal ) ){
  //   if( !currentline.BeginsWith( "#" ) ){
      
  //     TObjArray *tokens = currentline.Tokenize(" ");
  //     int ntokens = tokens->GetEntries();
  //     if( ntokens >= 11 ){
  // 	cout << currentline.Data() << ", ntokens = " << ntokens << endl;
	
  // 	TString snode = ( (TObjString*) (*tokens)[0] )->GetString();
  // 	int nodenumber = snode.Atoi();
	
  // 	TString scell = ( (TObjString*) (*tokens)[1] )->GetString();
  // 	int cellnumber = scell.Atoi();
	
  // 	TString speakpos = ( (TObjString*) (*tokens)[8] )->GetString();
  // 	double mean = speakpos.Atof();
	
  // 	TString ssigma = ( (TObjString*) (*tokens)[9] )->GetString();
  // 	double sigma = ssigma.Atof();

  // 	TString sthreshold = ( (TObjString*) (*tokens)[10] )->GetString();
  // 	double threshold = sthreshold.Atof();

  // 	TString srow = ( (TObjString*) (*tokens)[2] )->GetString();
  // 	TString scol = ( (TObjString*) (*tokens)[3] )->GetString();

  // 	std::pair<int,int> rowcoltemp( srow.Atoi(), scol.Atoi() );

  // 	cell_rowcol_ecal[rowcoltemp] = cellnumber;
	
  // 	list_of_nodes_ecal.insert( nodenumber );

  // 	cells_logic_sums_ecal[nodenumber].insert( cellnumber );

  // 	logic_mean_ecal[nodenumber] = mean;
  // 	logic_sigma_ecal[nodenumber] = sigma;
  // 	threshold_ecal[nodenumber] = threshold;

  // 	nodes_cells_ecal[ cellnumber ].insert(nodenumber);

  // 	TString sxcell = ( (TObjString*) (*tokens)[4] )->GetString(); 
  // 	TString sycell = ( (TObjString*) (*tokens)[5] )->GetString(); 

  // 	cols_cells_ecal[cellnumber] = scol.Atoi(); 
  // 	rows_cells_ecal[cellnumber] = srow.Atoi();

  // 	xcells_ecal[cellnumber] = sxcell.Atof()/1000.0; //convert to m
  // 	ycells_ecal[cellnumber] = sycell.Atof()/1000.0; //convert to m

  // 	if( ycell_rows.empty() || sycell.Atof()/1000.0 < ycellmin ) ycellmin = sycell.Atof()/1000.0;
  // 	if( ycell_rows.empty() || sycell.Atof()/1000.0 > ycellmax ) ycellmax = sycell.Atof()/1000.0;
	
  // 	ycell_rows[srow.Atoi()] = sycell.Atof()/1000.0;
  // 	TString ssize = ( (TObjString*) (*tokens)[6] )->GetString();
  // 	double size = ssize.Atof();

  // 	cellsize_rows[srow.Atoi()] = size/1000.0; 
	
  // 	if( xcellmin_rows.empty() || sxcell.Atof()/1000.0 < xcellmin_rows[srow.Atoi()] ){
  // 	  xcellmin_rows[srow.Atoi()] = sxcell.Atof()/1000.0;
  // 	}
  // 	if( xcellmax_rows.empty() || sxcell.Atof()/1000.0 > xcellmax_rows[srow.Atoi()] ){
  // 	  xcellmax_rows[srow.Atoi()] = sxcell.Atof()/1000.0;
  // 	}
	
  //     }
  //   }
  // }

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

      // logic_mean_hcal[current_node] = 1229.0;
      // logic_sigma_hcal[current_node] = 309.0;

      logic_mean_hcal[current_node] = 0.6;
      logic_sigma_hcal[current_node] = 0.2;
      //logic_sigma_hcal[current_node] = 346.3;
      
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
  ifstream thresholdfile_hcal(thresholdfilename_hcal);
  if( thresholdfile_hcal ){
    int node;
    double mean,sigma;
    while( thresholdfile_hcal >> node >> mean >> sigma ){
      if( list_of_nodes_hcal.find( node ) != list_of_nodes_hcal.end() ){
	logic_mean_hcal[ node ] = mean;
	logic_sigma_hcal[ node ] = sigma;
      }
    }
  }

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
  
  // TH2D *hnphesum_vs_node_ECAL_all = new TH2D("hnphesum_vs_node_ECAL_all","",list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5,100,0.0,5000.0);
  // TH2D *hnphesum_vs_node_HCAL_all = new TH2D("hnphesum_vs_node_HCAL_all","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);

  // //TH2D *hnphesum_vs_node_ECAL_FTcut = new TH2D("hnphesum_vs_node_ECAL_FTcut","",list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5,100,0.0,5000.0);
  // TH2D *hnphesum_vs_node_HCAL_FTcut = new TH2D("hnphesum_vs_node_HCAL_FTcut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  // TH2D *hnphesum_vs_node_HCAL_FPP1cut = new TH2D("hnphesum_vs_node_HCAL_FPP1cut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  // TH2D *hnphesum_vs_node_HCAL_FPP2cut = new TH2D("hnphesum_vs_node_HCAL_FPP2cut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  // TH2D *hnphesum_vs_node_HCAL_FPPbothcut = new TH2D("hnphesum_vs_node_HCAL_FPPbothcut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  // TH2D *hnphesum_vs_node_HCAL_FPPeithercut = new TH2D("hnphesum_vs_node_HCAL_FPPeithercut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);

  // TH2D *hnphesum_vs_node_HCALmax_all = new TH2D("hnphesum_vs_node_HCALmax_all","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  // TH2D *hnphesum_vs_node_HCALmax_FTcut = new TH2D("hnphesum_vs_node_HCALmax_FTcut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  // TH2D *hnphesum_vs_node_HCALmax_FPP1cut = new TH2D("hnphesum_vs_node_HCALmax_FPP1cut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  // TH2D *hnphesum_vs_node_HCALmax_FPP2cut = new TH2D("hnphesum_vs_node_HCALmax_FPP2cut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  // TH2D *hnphesum_vs_node_HCALmax_FPPbothcut = new TH2D("hnphesum_vs_node_HCALmax_FPPbothcut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  // TH2D *hnphesum_vs_node_HCALmax_FPPeithercut = new TH2D("hnphesum_vs_node_HCALmax_FPPeithercut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  
  
  // TH2D *hmaxnode_ECAL_vs_HCAL = new TH2D("hmaxnode_ECAL_vs_HCAL","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5);
  // TH2D *hallnodes_ECAL_vs_HCAL = new TH2D("hallnodes_ECAL_vs_HCAL","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5);

  // TH1D *hshouldhit_vs_threshold_ECAL = new TH1D("hshouldhit_vs_threshold_ECAL","",30,0.0,1.5);
  // TH1D *hefficiency_vs_threshold_ECAL = new TH1D("hefficiency_vs_threshold_ECAL","",30,0.0,1.5);
  // TH1D *hshouldhit_vs_threshold_ECAL_FTcut = new TH1D("hshouldhit_vs_threshold_ECAL_FTcut","",30,0.0,1.5);
  // TH1D *hefficiency_vs_threshold_ECAL_FTcut = new TH1D("hefficiency_vs_threshold_ECAL_FTcut","",30,0.0,1.5);
  
  // TH1D *hefficiency_vs_threshold_HCAL_FTcut = new TH1D("hefficiency_vs_threshold_HCAL_FTcut","",30,0.0,1.5);
  // TH1D *hefficiency_vs_threshold_HCAL_FPP1cut = new TH1D("hefficiency_vs_threshold_HCAL_FPP1cut","",30,0.0,1.5);
  // TH1D *hefficiency_vs_threshold_HCAL_FPP2cut = new TH1D("hefficiency_vs_threshold_HCAL_FPP2cut","",30,0.0,1.5);

  // TH1D *hshouldhit_HCAL_FTcut = new TH1D("hshouldhit_HCAL_FTcut","",30,0.0,1.5);
  // TH1D *hshouldhit_HCAL_FPP1cut = new TH1D("hshouldhit_HCAL_FPP1cut","",30,0.0,1.5);
  // TH1D *hshouldhit_HCAL_FPP2cut = new TH1D("hshouldhit_HCAL_FPP2cut","",30,0.0,1.5);
  
  // TH2D *hnphe_vs_sum_edep_ECAL = new TH2D("hnphe_vs_sum_edep_ECAL","",125,0.0,5.0,125,0.0,5000.0 );
  // TH2D *hnphe_vs_sum_edep_HCAL = new TH2D("hnphe_vs_sum_edep_HCAL","",125,0.0,0.75,125,0.0,2500.0 );

  //These histograms are the distributions of all trigger sums normalized as a fraction of elastic peak signal:
  TH1D *hEratio_allsums_ECAL = new TH1D("hEratio_allsums_ECAL","ECAL;E/E_{elastic}; rate", 200, 0.0, 2.0 );
  TH1D *hEratio_allsums_HCAL = new TH1D("hEratio_allsums_HCAL","HCAL;E/E_{elastic}; rate", 250, 0.0, 2.5 );

  TH2D *hEratio_logicsums_ECAL = new TH2D("hEratio_logicsums_ECAL", "ECAL; logic sum ; E/E_{elastic}", list_of_nodes_ecal.size(), 0.5, list_of_nodes_ecal.size() + 0.5, 100, 0., 2. );
  TH2D *hEratio_logicsums_HCAL = new TH2D("hEratio_logicsums_HCAL", "HCAL; logic sum; E/E_{elastic}", list_of_nodes_hcal.size(), 0.5, list_of_nodes_hcal.size() + 0.5, 125, 0., 2.5 );
  
  TH1D *hrate_vs_threshold_ECAL = new TH1D("hrate_vs_threshold_ECAL","",30,0.025,1.525);
  TH1D *hrate_vs_fixed_energy_threshold_ECAL = new TH1D("hrate_vs_fixed_energy_threshold_ECAL","",35,0.0,7.0 );
  TH1D *hrate_vs_threshold_HCAL = new TH1D("hrate_vs_threshold_HCAL","",30,0.025,1.525);
  TH1D *hrate_vs_fixed_energy_threshold_HCAL = new TH1D("hrate_vs_fixed_energy_threshold_HCAL","",30,0.0,10.0 );
  TH2D *hrate_vs_threshold_logic_sums_HCAL = new TH2D("hrate_vs_threshold_logic_sums_HCAL","",list_of_nodes_hcal.size(),0.5,0.5+list_of_nodes_hcal.size(),30,.025,1.525);
  TH2D *hrate_vs_threshold_logic_sums_HCALmax = new TH2D("hrate_vs_threshold_logic_sums_HCALmax","",list_of_nodes_hcal.size(),0.5,0.5+list_of_nodes_hcal.size(),30,.025,1.525);
  TH2D *hrate_E_vs_threshold_logic_sums_HCAL = new TH2D("hrate_E_vs_threshold_logic_sums_HCAL","",list_of_nodes_hcal.size(),0.5,0.5+list_of_nodes_hcal.size(),30,.025,1.525);
  //TH2D *hrate_vs_threshold_logic_sums_HCAL_nocointrig = new TH2D("hrate_vs_threshold_logic_sums_HCAL_nocointrig","",list_of_nodes_hcal.size(),0.5,0.5+list_of_nodes_hcal.size(),30,.025,1.525);
  TH2D *hrate_vs_threshold_logic_sums_ECAL = new TH2D("hrate_vs_threshold_logic_sums_ECAL","",list_of_nodes_ecal.size(),0.5,0.5+list_of_nodes_ecal.size(),30,.025,1.525);
  //TH2D *hrate_vs_threshold_logic_sums_ECAL_nocointrig = new TH2D("hrate_vs_threshold_logic_sums_ECAL_nocointrig","",list_of_nodes_ecal.size(),0.5,0.5+list_of_nodes_ecal.size(),30,.025,1.525);
  
  //TH2D *htrue_coincidence_rate_vs_threshold_ECAL_HCAL = new TH2D("htrue_coincidence_rate_vs_threshold_ECAL_HCAL","",30,.025,1.525,30,.025,1.525);
  TH2D *htrue_coincidence_rate_vs_threshold_ECAL_HCAL = new TH2D("htrue_coincidence_rate_vs_threshold_ECAL_HCAL","",30,.025,1.525,30,.025,1.525);
  TH2D *haccidental_coincidence_rate_vs_threshold_ECAL_HCAL = new TH2D("haccidental_coincidence_rate_vs_threshold_ECAL_HCAL","",30,.025,1.525,30,.025,1.525);
  //TH2D *haccidental_coincidence_rate_vs_threshold_ECAL_HCAL = new TH2D("haccidental_coincidence_rate_vs_threshold_ECAL_HCAL","",30,.025,1.525,30,.025,1.525);

  TH1D *hrate_ECAL_vs_ECAL_logic_sum = new TH1D("hrate_ECAL_vs_ECAL_logic_sum","",list_of_nodes_ecal.size(),0.5,0.5+list_of_nodes_ecal.size() );
  TH1D *hrate_HCAL_vs_ECAL_logic_sum = new TH1D("hrate_HCAL_vs_ECAL_logic_sum","",list_of_nodes_ecal.size(),0.5,0.5+list_of_nodes_ecal.size() );
  TH1D *hrate_ECAL_vs_HCAL_logic_sum = new TH1D("hrate_ECAL_vs_HCAL_logic_sum","",list_of_nodes_hcal.size(),0.5,0.5+list_of_nodes_hcal.size() );
  TH1D *hrate_HCAL_vs_HCAL_logic_sum = new TH1D("hrate_HCAL_vs_HCAL_logic_sum","",list_of_nodes_hcal.size(),0.5,0.5+list_of_nodes_hcal.size() );

  TH1D *hrealcoin_rate_vs_HCAL_logic_sum = new TH1D("hrealcoin_rate_vs_HCAL_logic_sum","",list_of_nodes_hcal.size(),0.5,0.5+list_of_nodes_hcal.size() );
  TH1D *haccidental_rate_vs_HCAL_logic_sum = new TH1D("haccidental_rate_vs_HCAL_logic_sum","",list_of_nodes_hcal.size(),0.5,0.5+list_of_nodes_hcal.size() );

  TH1D *hNnodes_ECAL_nominal_threshold = new TH1D("hNnodes_ECAL_nominal_threshold","",list_of_nodes_ecal.size(),0.5,0.5+list_of_nodes_ecal.size() );
  TH1D *hNnodes_HCAL_nominal_threshold = new TH1D("hNnodes_HCAL_nominal_threshold","",list_of_nodes_hcal.size(),0.5,0.5+list_of_nodes_hcal.size() ); 
  
  
  TH1D *hQ2 = new TH1D("hQ2","",200,0.0,10.0);
  TH1D *hW = new TH1D("hW","",200,0.0,6.0);
  TH1D *hy = new TH1D("hy","",200,0.0,1.0);
  TH2D *hQ2_theta_e = new TH2D("hQ2_theta_e","",200,0.0,45.0,200,0.0,10.0);
  TH2D *hW_theta_e = new TH2D("hW_theta_e","",200,0.0,45.0,200,0.0,10.0);
  TH2D *hQ2_W = new TH2D("hQ2_W","",200,0.0,6.0,200,0.0,10.0);
  TH2D *hy_theta_e = new TH2D("hy_theta_e","",200,0.0,45.0,200,0.0,1.0);
  
  
  // TH2D *hpion_p_vs_theta = new TH2D("hpion_p_vs_theta","",200,0.0,45.0,200,0.0,11.0 );
  // TH2D *hproton_p_vs_theta = new TH2D("hproton_p_vs_theta","",200,0.0,45.0,200,0.0,11.0);
  // TH2D *hneutron_p_vs_theta = new TH2D("hneutron_p_vs_theta","",200,0.0,45.0,200,0.0,11.0);
  // TH2D *hgamma_p_vs_theta = new TH2D("hgamma_p_vs_theta","",200,0.0,45.0,200,0.0,11.0);
  // TH2D *helectron_p_vs_theta = new TH2D("helectron_p_vs_theta","",200,0.0,45.0,200,0.0,11.0);
  // TH2D *hother_p_vs_theta = new TH2D("hother_p_vs_theta","",200,0.0,45.0,200,0.0,11.0);

  //plot the "overlap rates" between node i and j for ECAL and HCAL:
  TH3D *hrates_overlap_HCAL = new TH3D("hprob_overlap_HCAL","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,30,0.025,1.525);
  TH3D *hrates_overlap_ECAL = new TH3D("hprob_overlap_ECAL","",list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5,list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5,30,0.025,1.525);

  TH2D *hNnodes_fired_ECAL_vs_threshold = new TH2D("hNnodes_fired_ECAL_vs_threshold","",30,0.025,1.525,list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5);
  TH2D *hNnodes_fired_HCAL_vs_threshold = new TH2D("hNnodes_fired_HCAL_vs_threshold","",30,0.025,1.525,list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5);
  
  double Ibeam = 70.0e-6; //Amps
  double Ltarget = 30.0; //cm
  double e = 1.602e-19; //electron charge;
  double rho_target = 0.072; //g/cm^3
  double N_A = 6.022e23; //atoms/mol:
  double Mmol_H = 1.008; //g/mol
  //double Lumi = rho_target * Ltarget * N_A / Mmol_H * Ibeam/e; //~ 8e38;

  cout << "Lumi = " << Lumi << endl;
  
  TRandom3 num(0);

  cout << "Entering event loop " << endl;

  cout << "ngen total = " << ngen << endl;

  int treenum=-1;
  int oldtreenum=treenum;

  double real_coin_total_rate = 0.0;
  double total_real_coin_events = 0;
  double real_coin_sum2 = 0.0;
  
  long nevent=0, ntotal = C->GetEntries();
  for( nevent=0; nevent<C->GetEntries(); ++nevent ){
    T->GetEntry(nevent);

    treenum = C->GetTreeNumber();

    TString fname = C->GetFile()->GetName();

    //bool any_real_coin = false;
    
    if( bad_file_list.find(fname.Data()) == bad_file_list.end() ){ //good file
    
      if( treenum != oldtreenum ){
	cout << "new tree event " << nevent << endl;
	oldtreenum = treenum;

	cout << "file name = " << fname << endl;
	Lumi = LumiFile[fname];

	cout << "Luminosity = " << Lumi << endl;
      }

      //For beam generator, the normalization is simply incident electrons/second / total number of generated electrons
      double weight = Ibeam/e/double(ngen);
      //cross section is given in mb: 1 mb = 1e-3 * 1e-24 = 1e-27 cm^2
      // if (pythia6flag != 0 ){
      //no longer true: cross section is now given in cm2

      // if( !isnan(T->primaries_Sigma) ){
      
      // 	weight = Lumi * T->primaries_Sigma / double(ngen); //luminosity times cross section / number of events generated.
      // } else {
      // 	weight = Lumi * sigma_default / double(ngen);
      // }
      // } else {
      //weight = T->ev_rate / double(nfiles);
      //}

      // hQ2->Fill( T->primaries_Q2, weight );
      // hW->Fill( sqrt(T->primaries_W2), weight );
      // hy->Fill( T->primaries_y, weight );
      // hQ2_theta_e->Fill( 180.0/PI*T->primaries_theta_e, T->primaries_Q2, weight  );
      // hW_theta_e->Fill( 180.0/PI*T->primaries_theta_e, sqrt(T->primaries_W2), weight  );
      // hQ2_W->Fill( sqrt(T->primaries_W2), T->primaries_Q2, weight  );
      // hy_theta_e->Fill( 180.0/PI*T->primaries_theta_e, T->primaries_y, weight  );
    
      // for( int ipart=0; ipart<T->Primaries_Nprimaries; ipart++ ){
      // 	if( (*(T->Primaries_genflag))[ipart] == 1 ){
      // 	  if( fabs( (*(T->Primaries_PID))[ipart] ) == 211 ){ //charged pion:
      // 	    hpion_p_vs_theta->Fill( (*(T->Primaries_theta))[ipart]*180.0/PI, (*(T->Primaries_P))[ipart], weight );
      // 	  } else if( (*(T->Primaries_PID))[ipart] == 2212 ){ //proton:
      // 	    hproton_p_vs_theta->Fill( (*(T->Primaries_theta))[ipart]*180.0/PI, (*(T->Primaries_P))[ipart], weight );
      // 	  } else if( (*(T->Primaries_PID))[ipart] == 2112 ){ //neutron:
      // 	    hneutron_p_vs_theta->Fill( (*(T->Primaries_theta))[ipart]*180.0/PI, (*(T->Primaries_P))[ipart], weight );
      // 	  } else if( (*(T->Primaries_PID))[ipart] == 22 ){ //gamma:
      // 	    hgamma_p_vs_theta->Fill( (*(T->Primaries_theta))[ipart]*180.0/PI, (*(T->Primaries_P))[ipart], weight );
      // 	  } else if( fabs( (*(T->Primaries_PID))[ipart] ) == 11 ){ //e+/e-
      // 	    helectron_p_vs_theta->Fill( (*(T->Primaries_theta))[ipart]*180.0/PI, (*(T->Primaries_P))[ipart], weight );
      // 	  } else {
      // 	    hother_p_vs_theta->Fill( (*(T->Primaries_theta))[ipart]*180.0/PI, (*(T->Primaries_P))[ipart], weight );
      // 	  }
      // 	}
      // }
    
      //    double R = T->gen_dbb;
      //double thetacal = T->gen_thbb;
    
      if( (nevent+1) % 1000 == 0 ){ cout << "Event number " << nevent+1 << " of " << ntotal << ", event weight = " << weight << endl; }
    
      map<int,double> node_sums; //initialize all node sums to zero:
      for( set<int>::iterator inode = list_of_nodes_ecal.begin(); inode != list_of_nodes_ecal.end(); ++inode ){
	node_sums[ *inode ] = 0.0;
      }

      pheflag = 0;
      //  int nphe = 0;
      double nphe = 0.;
      
      // if( pheflag == 0 ){
      for( int ihit = 0; ihit<T->Earm_ECalTF1_hit_nhits; ihit++ ){
	int rowhit = ( *(T->Earm_ECalTF1_hit_row))[ihit]+1;
	int colhit = ( *(T->Earm_ECalTF1_hit_col))[ihit]+1;
	int cellhit = ( *(T->Earm_ECalTF1_hit_cell))[ihit];
	std::pair<int,int> rowcolhit( rowhit,colhit );
	
	//int cellhit = cell_rowcol_ecal[rowcolhit];
	
	
	//int trigger_group = nodes_cells_ecal[cellhit];
	
	double edep = (*(T->Earm_ECalTF1_hit_sumedep))[ihit];

	double mean = 528.*edep;
	//double sigma = 52.0*sqrt(edep) + 20.76*edep;
	nphe = num.Gaus(mean,sqrt(mean));
	
	//	nphe = TMath::Max(0,TMath::Nint(num.Gaus(mean,sigma)));
	double esmear = nphe/528.;
	
	for( set<int>::iterator inode = nodes_cells_ecal[cellhit].begin(); inode != nodes_cells_ecal[cellhit].end(); ++inode ){
	  node_sums[ *inode ] += esmear;
	}
	
      }
      // } else {
      // 	for( int ihit = 0; ihit<T->Earm_ECAL_hit_nhits; ihit++){
      // 	  int rowhit = ( *(T->Earm_ECAL_hit_row))[ihit]+1;
      // 	  int colhit = ( *(T->Earm_ECAL_hit_col))[ihit]+1;
      // 	  std::pair<int,int> rowcolhit( rowhit,colhit );
	
      // 	  int cellhit = cell_rowcol_ecal[rowcolhit];
	
      // 	  //int trigger_group = nodes_cells_ecal[cellhit];
	
      // 	  //	double edep = (*(T->Earm_ECalTF1_hit_sumedep))[ihit];

      // 	  int nphe = (*(T->Earm_ECAL_hit_NumPhotoelectrons))[ihit];
	
      // 	  for( set<int>::iterator inode = nodes_cells_ecal[cellhit].begin(); inode != nodes_cells_ecal[cellhit].end(); ++inode ){
      // 	    node_sums[ *inode ] += double(nphe);
      // 	  }
      // 	}
      
      // 	//node_sums[ trigger_group ] += double(nphe);
      // }

      // vector<int> trigger_nodes_fired_vs_E(hrate_vs_fixed_energy_threshold_ECAL->GetNbinsX());
      // vector<int> trigger_nodes_fired(hrate_vs_threshold_ECAL->GetNbinsX());
      // map<int, set<int> > list_of_nodes_fired_vs_threshold_ECAL; //map of threshold bins and list of fired ECAL nodes
      
      // for( int ithr=0; ithr<hrate_vs_threshold_ECAL->GetNbinsX(); ithr++ ){
      // 	trigger_nodes_fired[ithr] = 0;
      // }

      // for( int ithr=0; ithr<trigger_nodes_fired_vs_E.size(); ithr++ ){
      // 	trigger_nodes_fired_vs_E[ithr] = 0;
      // }

      int maxnode_ECAL=-1;
      int maxnode_HCAL=-1;
      double maxsum_ECAL = 0.0;
      double maxsum_HCAL = 0.0;

      int ntrig_ECAL_nominal_threshold = 0;
      int ntrig_HCAL_nominal_threshold = 0;
      
      for( set<int>::iterator inode = list_of_nodes_ecal.begin(); inode != list_of_nodes_ecal.end(); ++inode ){
	// for( int bin=1; bin<=hrate_vs_threshold_ECAL->GetNbinsX(); bin++ ){
	//   if( node_sums[*inode]/logic_mean_ecal[*inode] > hrate_vs_threshold_ECAL->GetBinCenter(bin) ){
	//     //cout << "node above threshold, nphe, peak position = " << node_sums[*inode] << ", " << logic_mean_ecal[*inode] << endl;
	//     trigger_nodes_fired[bin-1]++;
	//     hrate_vs_threshold_logic_sums_ECAL->Fill( *inode, hrate_vs_threshold_ECAL->GetBinCenter(bin), weight );
	//     list_of_nodes_fired_vs_threshold_ECAL[bin].insert( *inode );
	//   }
	// }

	if( node_sums[*inode] > 0.01*logic_mean_ecal[*inode] ){
	
	  hEratio_allsums_ECAL->Fill( node_sums[*inode]/logic_mean_ecal[*inode], weight );
	  hEratio_logicsums_ECAL->Fill( *inode, node_sums[*inode]/logic_mean_ecal[*inode], weight );

	}
	
	if( node_sums[*inode]/logic_mean_ecal[*inode] > ECAL_threshold ){
	  ntrig_ECAL_nominal_threshold++;
	}

	// for( int bin=1; bin<=hrate_vs_fixed_energy_threshold_ECAL->GetNbinsX(); bin++ ){
	//   if( node_sums[*inode]/752.2 >= hrate_vs_fixed_energy_threshold_ECAL->GetBinCenter(bin) ){
	//     trigger_nodes_fired_vs_E[bin-1]++;
	//   }
	// }

	//Always measure trigger sums as a fraction of elastic peak energy:
	if( node_sums[*inode]/logic_mean_ecal[*inode] > maxsum_ECAL ) {
	  maxsum_ECAL = node_sums[*inode]/logic_mean_ecal[*inode];
	  maxnode_ECAL = *inode;
	}
      
	//if( node_sums[*inode] > 0.0 ) hnphesum_vs_node_ECAL_all->Fill( *inode, node_sums[*inode], weight );
      }

      //For single-arm trigger rate calculation, we ask if the ECAL logic sum with the maximum amplitude is above threshold for that bin:
      for( int ithr=1; ithr<=hrate_vs_threshold_ECAL->GetNbinsX(); ithr++ ){
	double thrtemp = hrate_vs_threshold_ECAL->GetBinCenter( ithr ); 
	//if( maxsum_ECAL/logic_mean_ecal[maxnode_ECAL] > thrtemp ){
	if( maxsum_ECAL > thrtemp ){
	  hrate_vs_threshold_ECAL->Fill( thrtemp, weight );
	}
      }
      
      // for( int ithr=0; ithr<trigger_nodes_fired.size(); ithr++ ){
      // 	if( trigger_nodes_fired[ithr] > 0 ){
      // 	  double threshold = hrate_vs_threshold_ECAL->GetBinCenter( ithr + 1 );
      // 	  hrate_vs_threshold_ECAL->Fill( threshold, weight );
      // 	  hNnodes_fired_ECAL_vs_threshold->Fill( threshold, trigger_nodes_fired[ithr], weight );
      // 	  //for( set<int>::iterator inode = list_of_nodes_fired_vs_threshold_ECAL[ithr+1].begin(); inode != list_of_nodes_fired_vs_threshold_ECAL[ithr+1].end(); ++inode ){
      // 	  //  for( set<int>::iterator jnode = inode; jnode != list_of_nodes_fired_vs_threshold_ECAL[ithr+1].end(); ++jnode ){
      // 	  //    hrates_overlap_ECAL->Fill( *inode, *jnode, threshold, weight );
      // 	  //  }
      // 	  //}
      // 	}
      // }

      // for( int ithr=0; ithr<trigger_nodes_fired_vs_E.size(); ithr++ ){
      // 	if( trigger_nodes_fired_vs_E[ithr] > 0 ){
      // 	  double threshold = hrate_vs_fixed_energy_threshold_ECAL->GetBinCenter(ithr+1);
      // 	  hrate_vs_fixed_energy_threshold_ECAL->Fill( threshold, weight );
      // 	}
      // }
    
      map<int,double> node_sums_hcal;
      for( set<int>::iterator inode = list_of_nodes_hcal.begin(); inode != list_of_nodes_hcal.end(); ++inode ){
	node_sums_hcal[*inode] = 0.0;
      }

      //int nphe = 0;
    
      //if( pheflag == 0 ){
      for( int ihit=0; ihit<T->Harm_HCalScint_hit_nhits; ihit++ ){
	int rowhit = (*(T->Harm_HCalScint_hit_row))[ihit]+1;
	int colhit = (*(T->Harm_HCalScint_hit_col))[ihit]+1;
	std::pair<int,int> rowcolhit(rowhit,colhit);
	int cellhit = cell_rowcol_hcal[rowcolhit];
	//int trigger_group = nodes_cells_hcal[cellhit];
	double edep = (*(T->Harm_HCalScint_hit_sumedep))[ihit];
	//nphe = num.Poisson( phe_per_GeV_HCAL * edep );
	double mean = 2981.0*edep;
	//	double sigma = 69.54*sqrt(edep) + 155.3*edep;

	//nphe = TMath::Max(0,TMath::Nint(num.Gaus(mean,sigma)));
	nphe = num.Gaus(mean,sqrt(mean));

	double esmear = nphe/2981.;
	
	//cout << "HCAL hit " << ihit+1 << " node, edep, nphe = " << trigger_group << ", " << edep << ", " << nphe << endl;
	//node_sums_hcal[trigger_group] += double(nphe);
	for( set<int>::iterator inode = nodes_cells_hcal[cellhit].begin(); inode != nodes_cells_hcal[cellhit].end(); ++inode ){
	  
	  node_sums_hcal[*inode] += esmear;
	  
	}

	
      }
      // } else {
      // 	for( int jhit=0; jhit<T->Harm_HCal_hit_nhits; jhit++ ){
      // 	  int rowhit = (*(T->Harm_HCal_hit_row))[jhit]+1;
      // 	  int colhit = (*(T->Harm_HCal_hit_col))[jhit]+1;
      // 	  std::pair<int,int> rowcolhit(rowhit,colhit);
      // 	  int cellhit = cell_rowcol_hcal[rowcolhit];
      // 	  nphe = (*(T->Harm_HCal_hit_NumPhotoelectrons))[jhit];
      // 	  for( set<int>::iterator inode = nodes_cells_hcal[cellhit].begin(); inode != nodes_cells_hcal[cellhit].end(); ++inode ){
	  
      // 	    node_sums_hcal[*inode] += double(nphe);
	  
      // 	  }
      // 	}
      // }
    
      // vector<int> trigger_nodes_fired_hcal(hrate_vs_threshold_HCAL->GetNbinsX());
      // vector<int> trigger_nodes_fired_vs_E_hcal(hrate_vs_threshold_HCAL->GetNbinsX());
      // map<int, set<int> > list_of_nodes_fired_vs_threshold_HCAL; //mapping between threshold bins and list of fired HCAL nodes:
    
      // for( int ithr=0; ithr<hrate_vs_threshold_HCAL->GetNbinsX(); ithr++ ){
      // 	trigger_nodes_fired_hcal[ithr] = 0;
      // 	trigger_nodes_fired_vs_E_hcal[ithr] = 0;
      // }

      // vector<int> coin_trigger_fired( hrate_vs_threshold_HCAL->GetNbinsX()*hrate_vs_threshold_ECAL->GetNbinsX() );
      // for( int ithr=0; ithr<coin_trigger_fired.size(); ithr++ ){
      // 	coin_trigger_fired[ithr] = 0;
      // }

      //map<int,vector<int> > ECAL_nodes_fired_count; //count the number of ECAL nodes associated with each HCAL node that fire in this event, regardless of
      // whether HCAL fires, as a function of ECAL threshold

      map<int,int> ECAL_nodes_fired_count;

      bool any_real_coin = false;
      
      for( set<int>::iterator inode = list_of_nodes_hcal.begin(); inode != list_of_nodes_hcal.end(); ++inode ){
	// for( int bin=1; bin<=hrate_vs_threshold_HCAL->GetNbinsX(); bin++ ){
	//   if( node_sums_hcal[*inode]/logic_mean_hcal[*inode] > hrate_vs_threshold_HCAL->GetBinCenter(bin) ){
	//     trigger_nodes_fired_hcal[bin-1]++;
	//     hrate_vs_threshold_logic_sums_HCAL->Fill( *inode, hrate_vs_threshold_HCAL->GetBinCenter(bin), weight );
	//     //bool anycoin_hcal = false;

	//     list_of_nodes_fired_vs_threshold_HCAL[bin].insert( *inode );
	//   }

	  

	//   if( node_sums_hcal[*inode]/logic_mean_hcal[*inode]*6.4 >= hrate_vs_fixed_energy_threshold_HCAL->GetBinCenter(bin) ){
	//     trigger_nodes_fired_vs_E_hcal[bin-1]++;
	//   }
	// }
	if( node_sums_hcal[*inode] > 0.01*logic_mean_hcal[*inode] ){
	
	  hEratio_allsums_HCAL->Fill( node_sums_hcal[*inode]/logic_mean_hcal[*inode], weight );
	  hEratio_logicsums_HCAL->Fill( *inode, node_sums_hcal[*inode]/logic_mean_hcal[*inode], weight );
	}
	  
	if( node_sums_hcal[*inode]/logic_mean_hcal[*inode] > HCAL_threshold ){
	  ntrig_HCAL_nominal_threshold++;
	  hrate_HCAL_vs_HCAL_logic_sum->Fill( *inode, weight );
	}

	ECAL_nodes_fired_count[*inode] = 0;
	
	for( set<int>::iterator enode = ECAL_nodes_HCAL[*inode].begin(); enode != ECAL_nodes_HCAL[*inode].end(); ++enode ){
	  
	  //for( int ebin=1; ebin<=hrate_vs_threshold_ECAL->GetNbinsX(); ebin++ ){ //check ECAL sums:
	  if( node_sums[ *enode ]/logic_mean_ecal[*enode] > ECAL_threshold ){ //this ECAL sum fired:
	    //coin_trigger_fired[ (ebin-1) + (bin-1)*hrate_vs_threshold_ECAL->GetNbinsX() ]++;
	    //anycoin_hcal = true;
	    ECAL_nodes_fired_count[*inode]++;
	  }
	 
	}

	if( node_sums_hcal[*inode] > maxsum_HCAL ) {
	  maxsum_HCAL = node_sums_hcal[*inode];
	  maxnode_HCAL = *inode;
	}

	//if( node_sums_hcal[*inode] > HCAL_threshol
	
	if( ECAL_nodes_fired_count[*inode] > 0 ){
	  hrate_ECAL_vs_HCAL_logic_sum->Fill( *inode, weight );
	  if( node_sums_hcal[*inode]/logic_mean_hcal[*inode] > HCAL_threshold ){
	    hrealcoin_rate_vs_HCAL_logic_sum->Fill( *inode, weight );
	    any_real_coin = true;

	  }
	}
	
	

	// //ECAL_nodes_fired_count[*inode].resize( hrate_vs_threshold_ECAL->GetNbinsX() );
	// for( int ebin=1; ebin<=hrate_vs_threshold_ECAL->GetNbinsX(); ebin++ ){
	//   ECAL_nodes_fired_count[*inode][ebin-1] = 0;
	//   for( set<int>::iterator enode = ECAL_nodes_HCAL[*inode].begin(); enode != ECAL_nodes_HCAL[*inode].end(); ++enode ){
	//     if( node_sums[*enode]/logic_mean_ecal[*enode] > hrate_vs_threshold_ECAL->GetBinCenter(ebin) ){
	//       ECAL_nodes_fired_count[*inode][ebin-1]++;
	//     }
	//   }
	//   //regardless of whether HCAL fired, compute rate in associated ECAL sums:
	//   if( ECAL_nodes_fired_count[*inode][ebin-1] > 0 ) hrate_E_vs_threshold_logic_sums_HCAL->Fill( *inode, hrate_vs_threshold_ECAL->GetBinCenter(ebin), weight );
	// }
      
      }

      //Plot HCAL singles rate vs threshold:
      for( int ithr=1; ithr<=hrate_vs_threshold_HCAL->GetNbinsX(); ithr++ ){
	double thrtemp = hrate_vs_threshold_HCAL->GetBinCenter( ithr );

	if( maxsum_HCAL/logic_mean_hcal[maxnode_HCAL] >= thrtemp ){
	  hrate_vs_threshold_HCAL->Fill( thrtemp, weight );
	}
      }
      

      //if( maxsum_HCAL > HCAL_threshold && maxsum_ECAL > ECAL_threshold ){
      if( any_real_coin ){
	real_coin_total_rate += weight;
	real_coin_sum2 += pow(weight,2);
	total_real_coin_events += 1.0;
      }
      
      // for( int ithr=0; ithr<coin_trigger_fired.size(); ithr++ ){
      // 	int bin_e = ithr%(hrate_vs_threshold_ECAL->GetNbinsX())+1;
      // 	int bin_h = ithr/(hrate_vs_threshold_HCAL->GetNbinsX())+1;
      // 	double thr_e = htrue_coincidence_rate_vs_threshold_ECAL_HCAL->GetYaxis()->GetBinCenter(bin_e);
      // 	double thr_h = htrue_coincidence_rate_vs_threshold_ECAL_HCAL->GetXaxis()->GetBinCenter(bin_h);
      // 	if( coin_trigger_fired[ithr] > 0 ){
      // 	  htrue_coincidence_rate_vs_threshold_ECAL_HCAL->Fill( thr_h, thr_e, weight );
      // 	}
      // }

      hNnodes_HCAL_nominal_threshold->Fill( ntrig_HCAL_nominal_threshold );
      hNnodes_ECAL_nominal_threshold->Fill( ntrig_ECAL_nominal_threshold );
      
      // for( int ithr=0; ithr<trigger_nodes_fired_hcal.size(); ithr++ ){
      // 	double threshold = hrate_vs_threshold_HCAL->GetBinCenter( ithr + 1 );
      // 	if( trigger_nodes_fired_hcal[ithr] > 0 ){
      // 	  hrate_vs_threshold_HCAL->Fill( threshold, weight );
      // 	  hNnodes_fired_HCAL_vs_threshold->Fill( threshold, trigger_nodes_fired_hcal[ithr], weight );

      // 	  //for( set<int>::iterator inode = list_of_nodes_fired_vs_threshold_HCAL[ithr+1].begin(); inode != list_of_nodes_fired_vs_threshold_HCAL[ithr+1].end(); ++inode ){
      // 	  //  for( set<int>::iterator jnode = inode; jnode != list_of_nodes_fired_vs_threshold_HCAL[ithr+1].end(); ++jnode ){
	    
      // 	  //    hrates_overlap_HCAL->Fill( *inode, *jnode, threshold, weight );
	    
      // 	  //  }
      // 	  //}
      // 	}
      // 	if( node_sums_hcal[ maxnode_HCAL ]/logic_mean_hcal[maxnode_HCAL] >= threshold ){
      // 	  hrate_vs_threshold_logic_sums_HCALmax->Fill( maxnode_HCAL, threshold, weight );
      // 	}
      // 	if( trigger_nodes_fired_vs_E_hcal[ithr] > 0 ){
      // 	  hrate_vs_fixed_energy_threshold_HCAL->Fill( hrate_vs_fixed_energy_threshold_HCAL->GetBinCenter( ithr+1 ), weight );
      // 	}
      //}
    }
  }

  double sum_accidental = 0.0;
  double sum2_dacc = 0.0;
  
  double HCAL_avg_mult = hNnodes_HCAL_nominal_threshold->GetMean();
  
  for( int ih=1; ih<=list_of_nodes_hcal.size(); ih++ ){
    double hrate = hrate_HCAL_vs_HCAL_logic_sum->GetBinContent(ih);
    double erate = hrate_ECAL_vs_HCAL_logic_sum->GetBinContent(ih);

    double dhrate = hrate_HCAL_vs_HCAL_logic_sum->GetBinError(ih);
    double derate = hrate_ECAL_vs_HCAL_logic_sum->GetBinError(ih);
    
    double dt = 50.0e-9;

    double rate_accidental = hrate*erate*dt;

    haccidental_rate_vs_HCAL_logic_sum->SetBinContent( ih, rate_accidental );
    if( rate_accidental > 0 ){
      double drate_accidental = rate_accidental * sqrt(pow(derate/erate,2)+pow(dhrate/hrate,2));
      
      sum_accidental += rate_accidental;
      
      sum2_dacc += pow(drate_accidental,2);
    }
  }

  cout << "Total accidental coincidence rate = " << sum_accidental/HCAL_avg_mult << " +/- "
       << sqrt(sum2_dacc)/HCAL_avg_mult << endl;
  
  double realcoin_avg_weight = real_coin_total_rate / total_real_coin_events;
  double realcoin_avg_weight2 = real_coin_sum2 / total_real_coin_events;

  double realcoin_sigma_weight = sqrt(realcoin_avg_weight2 - pow(realcoin_avg_weight,2));

  double sigavg_realcoin = realcoin_sigma_weight/sqrt(total_real_coin_events);

  cout << "Total real coincidence rate = " << real_coin_total_rate << " +/- "
       << real_coin_total_rate / sqrt(TMath::Max(1.0,total_real_coin_events) ) << endl;
  
  
  //Accidental rate analysis:
  //In each threshold bin, we sum over all HCAL
  
  // How to analyze accidental coincidence rate? First
  // we loop over ECAL and HCAL threshold bins.
  // then we loop over HCAL logic groups. We compute the accidental coincidence rate for a
  // single HCAL logic group as the rate in that group times the rate at which any associated ECAL logic group fires.
  // Then, we sum over all HCAL logic groups, and divide by the average HCAL trigger multiplicity per event.

  // TString nametemp;
  
  // TH1D *htemp;
  // for( int ebin=1; ebin<=hrate_vs_threshold_ECAL->GetNbinsX(); ebin++ ){
  //   double ethr = hrate_vs_threshold_ECAL->GetBinCenter(ebin);

  //   //Get average trigger multiplicity per event for ECAL:  
  //   htemp = hNnodes_fired_ECAL_vs_threshold->ProjectionY( nametemp.Format("hNnodes_fired_ECAL_vs_threshold_py_xbin%d",ebin), ebin, ebin );

  //   double ECAL_avg_multiplicity = htemp->GetMean();
  //   cout << "Bin " << ebin << ", threshold = " << ethr << ", average trigger multiplicity = " << ECAL_avg_multiplicity << endl;
    
  //   for( int hbin=1; hbin<=hrate_vs_threshold_HCAL->GetNbinsX(); hbin++ ){
  //     double hthr = hrate_vs_threshold_HCAL->GetBinCenter(hbin);

  //     htemp = hNnodes_fired_HCAL_vs_threshold->ProjectionY( nametemp.Format("hNnodes_fired_HCAL_vs_threshold_py_xbin%d",hbin), hbin, hbin );

  //     double HCAL_avg_multiplicity = htemp->GetMean();

  //     cout << "Bin " << hbin << ", threshold = " << hthr << ", average trigger multiplicity = "
  // 	   << HCAL_avg_multiplicity << endl;
      
  //     double RealCoinRate = htrue_coincidence_rate_vs_threshold_ECAL_HCAL->GetBinContent( hbin, ebin );
  //     //Now for these values of ECAL threshold and HCAL threshold, we want to sum over all HCAL and divide by
  //     //average ECAL and HCAL trigger multiplicities:
  //     double dt = 30.0e-9; //30 ns coin-time window:
  //     double SumAccidentalRate = 0.0;
  //     for( set<int>::iterator inode=list_of_nodes_hcal.begin(); inode != list_of_nodes_hcal.end(); inode++ ){
  // 	double rate_H = hrate_vs_threshold_logic_sums_HCAL->GetBinContent( *inode, hbin );
  // 	double rate_E = hrate_E_vs_threshold_logic_sums_HCAL->GetBinContent( *inode, ebin );

  // 	if( rate_H > 0 && rate_E > 0 && HCAL_avg_multiplicity > 0. ){
  // 	  SumAccidentalRate += rate_H * rate_E * dt / HCAL_avg_multiplicity;
  // 	}
  //     }
  //     //The total coincidence rate includes real coincidences
  //     haccidental_coincidence_rate_vs_threshold_ECAL_HCAL->SetBinContent( hbin, ebin, TMath::Max(0.0,SumAccidentalRate ));
  //   }
  // }
  
 
  // TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  
  // c1->Divide(2,1);

  // c1->cd(1)->SetLogy();
  // hefficiency_vs_threshold_ECAL->SetMarkerStyle(20);
  // hefficiency_vs_threshold_ECAL->Draw();
  
  // c1->cd(2)->SetLogy();
  // hefficiency_vs_threshold_HCAL->SetMarkerStyle(20);
  // hefficiency_vs_threshold_HCAL->Draw();

  TH1D *hEintegral_allsums_ECAL = new TH1D(*hEratio_allsums_ECAL);
  TH2D *hEintegral_logicsums_ECAL = new TH2D(*hEratio_logicsums_ECAL);
  
  hEintegral_allsums_ECAL->SetName("hEintegral_allsums_ECAL");
  for( int i=1; i<=hEintegral_allsums_ECAL->GetNbinsX(); i++ ){
    double total, error;
    total = hEratio_allsums_ECAL->IntegralAndError( i, hEintegral_allsums_ECAL->GetNbinsX(), error );
    
    hEintegral_allsums_ECAL->SetBinContent(i, total);
    hEintegral_allsums_ECAL->SetBinError(i, error);

    for( int j=1; j<=hEratio_logicsums_ECAL->GetNbinsX(); j++ ){
      TH1D *htemp = hEratio_logicsums_ECAL->ProjectionY("",j,j);

      total = htemp->IntegralAndError( i, htemp->GetNbinsX(), error );

      int globalbin = hEratio_logicsums_ECAL->GetBin( j, i );
      
      hEintegral_logicsums_ECAL->SetBinContent( globalbin, total );
      hEintegral_logicsums_ECAL->SetBinError( globalbin, error );
    }
  }

  TH1D *hEintegral_allsums_HCAL = new TH1D(*hEratio_allsums_HCAL);
  hEintegral_allsums_HCAL->SetName("hEintegral_allsums_HCAL");
  for( int i=1; i<=hEintegral_allsums_HCAL->GetNbinsX(); i++ ){
    double total, error;
    total = hEratio_allsums_HCAL->IntegralAndError( i, hEintegral_allsums_HCAL->GetNbinsX(), error );
    
    hEintegral_allsums_HCAL->SetBinContent(i, total);
    hEintegral_allsums_HCAL->SetBinError(i, error);
  }
  
  fout->Write();
  fout->Close();
}
