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

void gep_trigger_analysis_elastic( const char *rootfilename, const char *logicfilename_ecal, const char *logicfilename_hcal, const char *outputfilename, int pheflag=0 ){

  double thetacal = thetacaldeg*PI/180.0;
  
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
  
  while( (chEl=(TChainElement*)next() )){
    TFile newfile(chEl->GetTitle());
    newfile.GetObject("run_data",rd);
    if( rd ){
      ngen += rd->fNtries;
      nfiles++;
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
  
  while( currentline.ReadLine( logicfile_ecal ) ){
    if( !currentline.BeginsWith( "#" ) ){
      
      TObjArray *tokens = currentline.Tokenize(" ");
      int ntokens = tokens->GetEntries();
      if( ntokens >= 11 ){
	cout << currentline.Data() << ", ntokens = " << ntokens << endl;
	
	TString snode = ( (TObjString*) (*tokens)[0] )->GetString();
	int nodenumber = snode.Atoi();
	
	TString scell = ( (TObjString*) (*tokens)[1] )->GetString();
	int cellnumber = scell.Atoi();
	
	TString speakpos = ( (TObjString*) (*tokens)[8] )->GetString();
	double mean = speakpos.Atof();
	
	TString ssigma = ( (TObjString*) (*tokens)[9] )->GetString();
	double sigma = ssigma.Atof();

	TString sthreshold = ( (TObjString*) (*tokens)[10] )->GetString();
	double threshold = sthreshold.Atof();

	TString srow = ( (TObjString*) (*tokens)[2] )->GetString();
	TString scol = ( (TObjString*) (*tokens)[3] )->GetString();

	std::pair<int,int> rowcoltemp( srow.Atoi(), scol.Atoi() );

	cell_rowcol_ecal[rowcoltemp] = cellnumber;
	
	list_of_nodes_ecal.insert( nodenumber );

	cells_logic_sums_ecal[nodenumber].insert( cellnumber );

	logic_mean_ecal[nodenumber] = mean;
	logic_sigma_ecal[nodenumber] = sigma;
	threshold_ecal[nodenumber] = threshold;

	nodes_cells_ecal[ cellnumber ].insert(nodenumber);

	TString sxcell = ( (TObjString*) (*tokens)[4] )->GetString(); 
	TString sycell = ( (TObjString*) (*tokens)[5] )->GetString(); 

	cols_cells_ecal[cellnumber] = scol.Atoi(); 
	rows_cells_ecal[cellnumber] = srow.Atoi();

	xcells_ecal[cellnumber] = sxcell.Atof()/1000.0; //convert to m
	ycells_ecal[cellnumber] = sycell.Atof()/1000.0; //convert to m

	if( ycell_rows.empty() || sycell.Atof()/1000.0 < ycellmin ) ycellmin = sycell.Atof()/1000.0;
	if( ycell_rows.empty() || sycell.Atof()/1000.0 > ycellmax ) ycellmax = sycell.Atof()/1000.0;
	
	ycell_rows[srow.Atoi()] = sycell.Atof()/1000.0;
	TString ssize = ( (TObjString*) (*tokens)[6] )->GetString();
	double size = ssize.Atof();

	cellsize_rows[srow.Atoi()] = size/1000.0; 
	
	if( xcellmin_rows.empty() || sxcell.Atof()/1000.0 < xcellmin_rows[srow.Atoi()] ){
	  xcellmin_rows[srow.Atoi()] = sxcell.Atof()/1000.0;
	}
	if( xcellmax_rows.empty() || sxcell.Atof()/1000.0 > xcellmax_rows[srow.Atoi()] ){
	  xcellmax_rows[srow.Atoi()] = sxcell.Atof()/1000.0;
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

  ifstream logicfile_hcal(logicfilename_hcal);

  current_node = 1;
  //  bool first_cell = true;

  while( currentline.ReadLine(logicfile_hcal) ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");
      int ntokens = tokens->GetEntries();
      if( ntokens >= 11 ){
	cout << currentline.Data() << ", ntokens = " << ntokens << endl;

	TString snode = ( (TObjString*) (*tokens)[0] )->GetString();
	int nodenumber = snode.Atoi();
	
	TString scell = ( (TObjString*) (*tokens)[1] )->GetString();
	int cellnumber = scell.Atoi();
	
	TString speakpos = ( (TObjString*) (*tokens)[8] )->GetString();
	double mean = speakpos.Atof();
	
	TString ssigma = ( (TObjString*) (*tokens)[9] )->GetString();
	double sigma = ssigma.Atof();

	TString sthreshold = ( (TObjString*) (*tokens)[10] )->GetString();
	double threshold = sthreshold.Atof();

	TString srow = ( (TObjString*) (*tokens)[2] )->GetString();
	TString scol = ( (TObjString*) (*tokens)[3] )->GetString();

	std::pair<int,int> rowcoltemp( srow.Atoi(), scol.Atoi() );

	cell_rowcol_hcal[rowcoltemp] = cellnumber;
	
	list_of_nodes_hcal.insert( nodenumber );

	cells_logic_sums_hcal[nodenumber].insert( cellnumber );

	logic_mean_hcal[nodenumber] = mean;
	logic_sigma_hcal[nodenumber] = sigma;
	threshold_hcal[nodenumber] = threshold;

	nodes_cells_hcal[ cellnumber ].insert(nodenumber);
	
      }
    }
  }
  
  TH1D::SetDefaultSumw2();

  //double PI = TMath::Pi();

  //Photoelectron statistics:
  double phe_per_GeV_ECAL = 1000.0/1.33; //~ 750 pe/GeV
  double phe_per_GeV_HCAL = 1000.0/0.30; //~ 3,333 pe/GeV (but sampling fraction is small)

  fout->cd();

  // TH1D *hrate_vs_threshold_ECAL = new TH1D("hrate_vs_threshold_ECAL","",30,0.0,1.5);
  // //TH1D *hnum_logic_sums_fired_vs_threshold = new TH1D("hnum_logic_sums_fired_vs_threshold

  // TH1D *hrate_vs_threshold_HCAL = new TH1D("hrate_vs_threshold_HCAL","",40,0.0,2.0);

  //TH2D *htrue_coincidence_rate_vs_threshold_ECAL_HCAL = new TH2D("htrue_coincidence_rate_vs_threshold_ECAL_HCAL","",40,0,2.0,30,0,1.5);
  
  TH2D *hnphesum_vs_node_ECAL_all = new TH2D("hnphesum_vs_node_ECAL_all","",list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5,100,0.0,5000.0);
  TH2D *hnphesum_vs_node_HCAL_all = new TH2D("hnphesum_vs_node_HCAL_all","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);

  TH2D *hnphesum_vs_node_ECAL_FTcut = new TH2D("hnphesum_vs_node_ECAL_FTcut","",list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5,100,0.0,5000.0);
  TH2D *hnphesum_vs_node_HCAL_FTcut = new TH2D("hnphesum_vs_node_HCAL_FTcut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);

  TH2D *hnphesum_vs_node_HCAL_FPP1cut = new TH2D("hnphesum_vs_node_HCAL_FPP1cut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCAL_FPP2cut = new TH2D("hnphesum_vs_node_HCAL_FPP2cut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCAL_FPPbothcut = new TH2D("hnphesum_vs_node_HCAL_FPPbothcut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCAL_FPPeithercut = new TH2D("hnphesum_vs_node_HCAL_FPPeithercut","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  
  
  TH2D *hmaxnode_ECAL_vs_HCAL = new TH2D("hmaxnode_ECAL_vs_HCAL","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5);
  TH2D *hallnodes_ECAL_vs_HCAL = new TH2D("hallnodes_ECAL_vs_HCAL","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5);

  TH2D *hmaxnode_ECAL_vs_HCAL = new TH2D("hmaxnode_ECAL_vs_HCAL","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5);
  TH2D *hallnodes_ECAL_vs_HCAL = new TH2D("hallnodes_ECAL_vs_HCAL","",list_of_nodes_hcal.size(),0.5,list_of_nodes_hcal.size()+0.5,list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5);

  TH1D *hshouldhit_vs_threshold_ECAL = new TH1D("hshouldhit_vs_threshold_ECAL","",30,0.0,1.5);
  TH1D *hefficiency_vs_threshold_ECAL = new TH1D("hefficiency_vs_threshold_ECAL","",30,0.0,1.5);
  TH1D *hshouldhit_vs_threshold_ECAL_FTcut = new TH1D("hshouldhit_vs_threshold_ECAL_FTcut","",30,0.0,1.5);
  TH1D *hefficiency_vs_threshold_ECAL_FTcut = new TH1D("hefficiency_vs_threshold_ECAL_FTcut","",30,0.0,1.5);
  TH1D *hefficiency_vs_threshold_HCAL_FTcut = new TH1D("hefficiency_vs_threshold_HCAL_FTcut","",30,0.0,1.5);
  TH1D *hefficiency_vs_threshold_HCAL_FPP1cut = new TH1D("hefficiency_vs_threshold_HCAL_FPP1cut","",30,0.0,1.5);
  TH1D *hefficiency_vs_threshold_HCAL_FPP2cut = new TH1D("hefficiency_vs_threshold_HCAL_FPP2cut","",30,0.0,1.5);

  TH2D *hnphe_vs_sum_edep_ECAL = new TH2D("hnphe_vs_sum_edep_ECAL","",250,0.0,5.0,250,0.0,5000.0 );
  TH2D *hnphe_vs_sum_edep_HCAL = new TH2D("hnphe_vs_sum_edep_HCAL","",250,0.0,1.0,250,0.0,2500.0 );
  
  
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

    bool FTtrack = false;
    for( int itrack=0; itrack<T->Harm_FT_Track_ntracks; itrack++ ){
      if( (*(T->Harm_FT_Track_MID))[itrack] == 0 &&
	  (*(T->Harm_FT_Track_PID))[itrack] == 2212 ){ //primary elastically scattered proton track in FT:
	FTtrack = true;
      }
    }

    double 
    bool FPP1_singletrack = false, FPP2_singletrack = false;
    if( FTtrack )
    
    double nu = T->ev_Q2 / 2.0 / 0.938272;
    double pp_elastic = sqrt(pow(nu,2)+2.0*.938272*nu);
    
    double weight;
    //cross section is given in mb: 1 mb = 1e-3 * 1e-24 = 1e-27 cm^2
    // if (pythia6flag != 0 ){
    //   weight = Lumi * T->primaries_Sigma * 1.0e-27/ double(ngen); //luminosity times cross section / number of events generated.
    // } else {
    weight = T->ev_rate / double(nfiles);
      //}

    double R = T->gen_dbb;
    double thetacal = T->gen_thbb;
    
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

    if( ycalo >= ycellmin && ycalo <= ycellmax ){
      //make an initial guess at which row: (row runs from 1 to N):
      int closest_row = int( (ycalo - ycellmin)/4.0 ) + 1;

      map<int,double>::iterator rowguess = ycell_rows.find( closest_row );

      while( rowguess != ycell_rows.end() && ycalo > ycell_rows[rowguess->first] + 0.5*cellsize_rows[rowguess->first] ){ ++rowguess; }
      while( rowguess != ycell_rows.end() && ycalo < ycell_rows[rowguess->first] - 0.5*cellsize_rows[rowguess->first] ){ --rowguess; }

      if( rowguess != ycell_rows.end() ){
	closest_row = rowguess->first;
	if( xcalo >= xcellmin_rows[closest_row] + 0.5*cellsize_rows[closest_row] &&
	    xcalo <= xcellmax_rows[closest_row] - 0.5*cellsize_rows[closest_row] &&
	    ycalo >= ycellmin + 0.5*cellsize_rows[closest_row] && ycalo <= ycellmax - 0.5*cellsize_rows[closest_row] ){
	  should_hit_ECAL = true;
	}
      }
    }
    
    for( int ihit = 0; ihit<T->Earm_ECalTF1_hit_nhits; ihit++ ){
      int rowhit = ( *(T->Earm_ECalTF1_hit_row))[ihit]+1;
      int colhit = ( *(T->Earm_ECalTF1_hit_col))[ihit]+1;
      std::pair<int,int> rowcolhit( rowhit,colhit );

      int cellhit = cell_rowcol_ecal[rowcolhit];

      //int trigger_group = nodes_cells_ecal[cellhit];
      
      double edep = (*(T->Earm_ECalTF1_hit_sumedep))[ihit];

      if( pheflag != 0 ){
	for( int jhit = 0; jhit<T->Earm_ECAL_hit_nhits; jhit++ ){
	  if( (*(T->Earm_ECAL_hit_row))[jhit] == rowhit &&
	      (*(T->Earm_ECAL_hit_col))[jhit] == colhit &&
	      fabs( (*(T->Earm_ECAL_hit_Time_avg))[jhit] - (*(T->Earm_ECalTF1_hit_tavg))[ihit] - 2.5 ) <= 10.0 ){
	    hnphe_sum_vs_edep_ECAL->Fill( edep, (*(T->Earm_ECAL_hit_NumPhotoelectrons))[jhit] );
	    for( set<int>::iterator inode = nodes_cells_ecal[cellhit].begin(); inode != nodes_cells_ecal[cellhit].end(); ++inode ){
	      node_sums[*inode] += double( (*(T->Earm_ECAL_hit_NumPhotoelectrons))[jhit] );
	    }
	  }
	}
      } else {
      
	int nphe = num.Poisson( phe_per_GeV_ECAL * edep );
	
	for( set<int>::iterator inode = nodes_cells_ecal[cellhit].begin(); inode != nodes_cells_ecal[cellhit].end(); ++inode ){
	  node_sums[ *inode ] += double(nphe);
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
    
    for( set<int>::iterator inode = list_of_nodes_ecal.begin(); inode != list_of_nodes_ecal.end(); ++inode ){
      for( int bin=1; bin<=hefficiency_vs_threshold_ECAL->GetNbinsX(); bin++ ){
	if( node_sums[*inode]/logic_mean_ecal[*inode] > hefficiency_vs_threshold_ECAL->GetBinLowEdge(bin) ){
	  //cout << "node above threshold, nphe, peak position = " << node_sums[*inode] << ", " << logic_mean_ecal[*inode] << endl;
	  trigger_nodes_fired[bin-1]++;
	}
	
      }
      if( node_sums[*inode] > maxsum_ECAL ) {
	maxsum_ECAL = node_sums[*inode];
	maxnode_ECAL = *inode;
      }
      
      if( node_sums[*inode] > 0.0 ) hnphesum_vs_node_ECAL->Fill( *inode, node_sums[*inode], weight );
    }

    if( should_hit_ECAL ){
      for( int ithr=0; ithr<hefficiency_vs_threshold_ECAL->GetNbinsX(); ithr++ ){
	hshouldhit_vs_threshold_ECAL->Fill( hefficiency_vs_threshold_ECAL->GetBinCenter(ithr+1), weight );
	if( trigger_nodes_fired[ithr] > 0 ){
	  hefficiency_vs_threshold_ECAL->Fill( hefficiency_vs_threshold_ECAL->GetBinCenter(ithr+1), weight );
	}
      }
    }
    
    map<int,double> node_sums_hcal;
    for( set<int>::iterator inode = list_of_nodes_hcal.begin(); inode != list_of_nodes_hcal.end(); ++inode ){
      node_sums_hcal[*inode] = 0.0;
    }
    
    for( int ihit=0; ihit<T->Harm_HCalScint_hit_nhits; ihit++ ){
      int rowhit = (*(T->Harm_HCalScint_hit_row))[ihit]+1;
      int colhit = (*(T->Harm_HCalScint_hit_col))[ihit]+1;
      std::pair<int,int> rowcolhit(rowhit,colhit);
      int cellhit = cell_rowcol_hcal[rowcolhit];
      //int trigger_group = nodes_cells_hcal[cellhit];
      double edep = (*(T->Harm_HCalScint_hit_sumedep))[ihit];
      int nphe = num.Poisson( phe_per_GeV_HCAL * edep );
      //cout << "HCAL hit " << ihit+1 << " node, edep, nphe = " << trigger_group << ", " << edep << ", " << nphe << endl;
      //node_sums_hcal[trigger_group] += double(nphe);
      for( set<int>::iterator inode = nodes_cells_hcal[cellhit].begin(); inode != nodes_cells_hcal[cellhit].end(); ++inode ){
	
	node_sums_hcal[*inode] += double(nphe);
	
      }
    }
    
    vector<int> trigger_nodes_fired_hcal(hrate_vs_threshold_HCAL->GetNbinsX());
    for( int ithr=0; ithr<hrate_vs_threshold_HCAL->GetNbinsX(); ithr++ ){
      trigger_nodes_fired_hcal[ithr] = 0;
    }
    
    for( set<int>::iterator inode = list_of_nodes_hcal.begin(); inode != list_of_nodes_hcal.end(); ++inode ){
      for( int bin=1; bin<=hrate_vs_threshold_HCAL->GetNbinsX(); bin++ ){
	if( node_sums_hcal[*inode]/logic_mean_hcal[*inode] > hrate_vs_threshold_HCAL->GetBinLowEdge(bin) ){
	  trigger_nodes_fired_hcal[bin-1]++;
	}
      }
      if( node_sums_hcal[*inode] > maxsum_HCAL ) {
	maxsum_HCAL = node_sums_hcal[*inode];
	maxnode_HCAL = *inode;
      }
	
      if( node_sums_hcal[ *inode ] > 0.0 ) hnphesum_vs_node_HCAL->Fill( *inode, node_sums_hcal[*inode], weight );
    }
    for( int ithr=0; ithr<hrate_vs_threshold_HCAL->GetNbinsX(); ithr++ ){
      if( trigger_nodes_fired_hcal[ithr] > 0 ) hrate_vs_threshold_HCAL->Fill( hrate_vs_threshold_HCAL->GetBinCenter(ithr+1),weight );
      for( int jthr=0; jthr<hefficiency_vs_threshold_ECAL->GetNbinsX(); jthr++ ){
	if( trigger_nodes_fired[jthr] > 0 && trigger_nodes_fired_hcal[ithr] > 0 ){
	  htrue_coincidence_rate_vs_threshold_ECAL_HCAL->Fill( hrate_vs_threshold_HCAL->GetBinCenter(ithr+1),hefficiency_vs_threshold_ECAL->GetBinCenter(jthr+1),weight );
	}
      }
    }

    for( set<int>::iterator inode = list_of_nodes_ecal.begin(); inode != list_of_nodes_ecal.end(); ++inode ){
      for( set<int>::iterator jnode = list_of_nodes_hcal.begin(); jnode != list_of_nodes_hcal.end(); ++jnode ){
	//Fill the correlation histogram for all true coincidence events for which ECAL and HCAL node are both above threshold:
	if( node_sums[*inode] >= nominal_threshold_ECAL*logic_mean_ecal[*inode] && node_sums_hcal[*jnode] >= nominal_threshold_HCAL*logic_mean_hcal[*jnode] ){
	  hallnodes_ECAL_vs_HCAL->Fill( *jnode, *inode, weight );
	}
      }
    }
    if( maxsum_ECAL >= nominal_threshold_ECAL*logic_mean_ecal[maxnode_ECAL] && maxsum_HCAL >= nominal_threshold_HCAL*logic_mean_hcal[maxnode_HCAL] ){
      hmaxnode_ECAL_vs_HCAL->Fill( maxnode_HCAL, maxnode_ECAL, weight );
    }
    //}
  }
 
  // TCanvas *c1 = new TCanvas("c1","c1",1200,900);
  
  // c1->Divide(2,1);

  // c1->cd(1)->SetLogy();
  // hefficiency_vs_threshold_ECAL->SetMarkerStyle(20);
  // hefficiency_vs_threshold_ECAL->Draw();
  
  // c1->cd(2)->SetLogy();
  // hrate_vs_threshold_HCAL->SetMarkerStyle(20);
  // hrate_vs_threshold_HCAL->Draw();

  hrate_vs_threshold_HCAL->SetMarkerStyle(20);
  hefficiency_vs_threshold_ECAL->SetMarkerStyle(20);
  
  fout->Write();
  fout->Close();
}
