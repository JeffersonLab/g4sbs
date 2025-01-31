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
#include <algorithm>

using namespace std;

TF1 *gaussian = new TF1("gaussian", "[0]*exp(-0.5*pow((x-[1])/[2],2))", 0.0,5000.0);
TF1 *gaussplusexpo = new TF1("gaussplusexpo", "[0]*exp(-0.5*pow((x-[1])/[2],2))+[3]*exp(-[4]*x)",0.0,5000.0);

//const double Rcal = 4.9; //m, to back surface
const double Mp = 0.938272; //GeV
const double PI = TMath::Pi();

struct calblock {
  int cell;
  int row;
  int col;
  double x;
  double y;
  double E; //smeared energy
  double Etrue; //true energy
  double nphe; //estimated nphe
};

struct calcluster {
  int node; //"node" index
  int cell; //cell index of center cell;
  int bin;  //bin number of cluster = binx + biny * nbinsx
  int binx,biny;
  double Esum_true;
  double nphesum;
  double Esum; //cluster energy sum
  double Esum_norm; //cluster energy sum / elastic peak energy
};

bool CompareCalBlocks( const calblock &b1, const calblock &b2 ){
  return b1.E > b2.E;
}

bool CompareCalClusters( const calcluster &c1, const calcluster &c2 ){
  return c1.Esum > c2.Esum; 
}

void gep_trigger_analysis_elastic_L2( const char *rootfilename, const char *outputfilename="gep_L2_trigger_analysis_elastic_temp.root", const char *assocfilename="database/ECAL_HCAL_L2_default.txt", double normfac_ECAL=0.88, double normfac_HCAL=1.064, double thresh_ECAL=0.8, double thresh_HCAL=0.5, double pdeflect_default=0.7853, int pheflag=0){

  double nominal_threshold_HCAL = thresh_HCAL;
  double nominal_threshold_ECAL = thresh_ECAL;
  
  TFile *fout = new TFile(outputfilename,"RECREATE");
  TChain *C = new TChain("T");
  C->Add(rootfilename);

  gep_tree_elastic *T = new gep_tree_elastic( C );

  C->SetBranchStatus("*",0);

  C->SetBranchStatus("ev*",1);
  C->SetBranchStatus("Harm.FT.Track.*",1);
  C->SetBranchStatus("Harm.FPP1.Track.*",1);
  C->SetBranchStatus("Earm.ECalTF1.hit.*",1);
  C->SetBranchStatus("Harm.HCalScint.hit.*",1);
  
  G4SBSRunData *rd;

  long ngen = 0;
  int nfiles = 0;
  
  TObjArray *FileList = C->GetListOfFiles();
  TIter next(FileList);

  TChainElement *chEl = 0;

  set<TString> bad_file_list;

  double Ebeam_default = 10.688;
  double ECALtheta_default = 29.75*TMath::DegToRad();
  double ECALdist_default = 4.7;
  double HCALdist_default = 10.0;
  double HCALvoff_default = 0.75;
  double SBSdist_default = 1.6;
  double SBStheta_default = 16.9*TMath::DegToRad();

  //double pdeflect_default = 0.7853;
  
  map<TString,double> Ebeam_file;
  map<TString,double> ECALdist_file;
  map<TString,double> ECALtheta_file;
  map<TString,double> HCALdist_file;
  map<TString,double> HCALvoff_file;
  map<TString,double> SBStheta_file;
  map<TString,double> SBSdist_file;
  map<TString,double> Lumi_file;
  map<TString,double> GenVol_file;
  
  while( (chEl=(TChainElement*)next() )){
    TFile newfile(chEl->GetTitle());
    newfile.GetObject("run_data",rd);
    if( rd ){
      ngen += rd->fNtries;
      nfiles++;

      ECALdist_file[chEl->GetTitle()] = rd->fBBdist;
      ECALtheta_file[chEl->GetTitle()] = rd->fBBtheta;
      Ebeam_file[chEl->GetTitle()] = rd->fBeamE;
      Lumi_file[chEl->GetTitle()] = rd->fLuminosity;
      GenVol_file[chEl->GetTitle()] = rd->fGenVol;
      HCALdist_file[chEl->GetTitle()] = rd->fHCALdist;
      HCALvoff_file[chEl->GetTitle()] = rd->fHCALvoff;
      SBSdist_file[chEl->GetTitle()] = rd->fSBSdist;
      SBStheta_file[chEl->GetTitle()] = rd->fSBStheta;
      
      Ebeam_default = rd->fBeamE;
      ECALtheta_default = rd->fBBtheta;
      ECALdist_default = rd->fBBdist;
      HCALdist_default = rd->fHCALdist;
      HCALvoff_default = rd->fHCALvoff;
      SBSdist_default = rd->fSBSdist;
      SBStheta_default = rd->fSBStheta;
      
    } else {
      bad_file_list.insert( chEl->GetTitle());
    }
  }

  cout << "number of generated events = " << ngen << endl;

  TVector3 ecal_zaxis(sin(ECALtheta_default),0,cos(ECALtheta_default));
  TVector3 ecal_yaxis(0,1,0);
  TVector3 ecal_xaxis = ecal_yaxis.Cross(ecal_zaxis).Unit();
  TVector3 ecal_pos = ECALdist_default * ecal_zaxis;
  
  //Make every cell a cluster center:
  
  set<int> list_of_nodes_ecal;
  map<int, int> cluster_centers_ecal_by_node; //Record the center cell of each cluster!
  map<int, int> nodes_ecal_by_cluster_center; //key = cell, mapped value = node
  map<int, set<int> > cells_logic_sums_ecal; //mapping between node numbers and cell numbers
  map<int, double> logic_mean_ecal; //mean peak positions by node number
  map<int, double> logic_sigma_ecal; //peak width by node number
  map<int, double> threshold_ecal; //threshold by node number
  map<std::pair<int,int>, int > cell_rowcol_ecal; //cell numbers mapped by unique row and column pairs
  map<int,set<int> > nodes_cells_ecal; //mapping of nodes by cell number:

  map<int,double> Eprime_expect_by_node_ecal; 
  
  //redo everything in terms of fixed size arrays:
  int ncells_ecal;
  set<int> cells_ecal;      //If cells are somehow not in order: 
  map<int,int> rows_cells_ecal;
  map<int,int> cols_cells_ecal;
  
  map<int,double> xcell_ecal;
  map<int,double> ycell_ecal;
  map<int,int> nearestcol_up_ecal;
  map<int,int> nearestcol_down_ecal;
  //  map<int,int> bins_cells_ecal; //This is for later
  set<int> list_of_bins_ecal; //list of bins containing at least one cluster center
  
  map<int,int> binx_by_node_ecal;
  map<int,int> biny_by_node_ecal;
  map<int,int> binglobal_by_node_ecal;

  map<int,int> binx_by_cell_ecal;
  map<int,int> biny_by_cell_ecal;
  map<int,int> binglobal_by_cell_ecal;

  map<int,set<int> > cells_by_binglobal_ecal;
  map<int,set<int> > nodes_by_binglobal_ecal;


  
  // map<int,int> rows_cells_ecal;
  // map<int,int> cols_cells_ecal;
  // map<int,double> xcells_ecal;
  // map<int,double> ycells_ecal;
  
  //keep track of min and max x by row number:
  double ycellmin,ycellmax;
  double xcellmin,xcellmax;
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
  //ifstream logicfile_ecal("database/GEP_ECAL_L2sums.txt");
  //ifstream thresholdfile(thresholdfilename);

  //Cell runs from 0 to N-1
  
  TString currentline;

  ycellmin = 1.0e9;
  ycellmax = -1.0e9;
  xcellmin = 1.e9;
  xcellmax = -1.e9;


  //Get cell mapping info:
  while( currentline.ReadLine( mapfile_ecal ) ){
    //    cout << currentline << endl;
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
	double xcell = sxcell.Atof();
	double ycell = sycell.Atof();

	std::pair<int, int> rowcoltemp(row,col);

	rows_ecal.insert(row);
	columns_rows_ecal[row].insert(col); //if we're lazy we can ASSUME that column within a row runs from 0 to ncolumns - 1.
	
	cells_ecal.insert( cellnum );
	cell_rowcol_ecal[rowcoltemp] = cellnum;
	rows_cells_ecal[cellnum] = row;
	cols_cells_ecal[cellnum] = col;
	xcell_ecal[cellnum] = xcell/100.0;
	ycell_ecal[cellnum] = ycell/100.0;

	if( ycell_rows.empty() || ycell/100.0 < ycellmin ) ycellmin = ycell/100.0;
	if( ycell_rows.empty() || ycell/100.0 > ycellmax ) ycellmax = ycell/100.0;
	  
	ycell_rows[row] = ycell/100.0;

	maxrow = (row > maxrow) ? row : maxrow;
	minrow = (row < minrow) ? row : minrow;
	
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

	xcellmin = xcell/100.0 < xcellmin ? xcell/100.0 : xcellmin;
	xcellmax = xcell/100.0 > xcellmax ? xcell/100.0 : xcellmax;
	
      }
    }
  }

  std::cout << "(xcellmin, xcellmax)=(" << xcellmin << ", " << xcellmax << ") m" << std::endl;
  std::cout << "(ycellmin, ycellmax)=(" << ycellmin << ", " << ycellmax << ") m" << std::endl;

  //we want to superimpose a rectangular binning on the ECAL clusters:
  
  ncells_ecal = cells_ecal.size();
  //The SAFEST approach to forming the cluster groupings is to go by rows and columns:

  for( auto row : rows_ecal ){
    for( auto col : columns_rows_ecal[row] ){
      std::pair<int,int> rowcoltemp( row, col );
      int cell = cell_rowcol_ecal[rowcoltemp];

      //grab position, although I think we are mainly interested in x (horizontal) here:
      double xcell = xcell_ecal[cell];
      double ycell = ycell_ecal[cell];
      
      if( row+1 < rows_ecal.size() ){ //if the row above exists, find nearest column:
	int nearestcol=-1;
	double mindiff = 1000.0;
	for( auto colup : columns_rows_ecal[row+1] ){
	  std::pair<int,int> rowcolup(row+1,colup);
	  int cellup = cell_rowcol_ecal[rowcolup];
	  double xup = xcell_ecal[cellup];

	  bool goodx = fabs(xup-xcell)<0.5*0.043 && fabs(xup-xcell)<mindiff;
	  
	  if( goodx ){
	    nearestcol = colup;
	    mindiff = fabs(xup-xcell);
	  }
	}

	nearestcol_up_ecal[cell] = nearestcol;
      } else {
	nearestcol_up_ecal[cell] = -1;
      }

      if( row > 0 ){ //if the row below exists, find nearest column below:
	int nearestcol = -1;
	double mindiff = 1000.0;
	for( auto coldown : columns_rows_ecal[row-1] ){
	  std::pair<int,int> rowcoldown(row-1,coldown);
	  int celldown = cell_rowcol_ecal[rowcoldown];
	  double xdown = xcell_ecal[celldown];
	  bool goodx = fabs(xdown-xcell)<0.5*0.043 && fabs(xdown-xcell) < mindiff;
	  if( goodx ){
	    nearestcol = coldown;
	    mindiff = abs(xdown-xcell);
	  }
	}
	nearestcol_down_ecal[cell] = nearestcol;
      } else {
	nearestcol_down_ecal[cell] = -1;
      }	
    }
  }

  int node = 0;

  double xspace = 0.04292;
  double yspace = 0.04293;

  //We impose a regular rectangular binning on ECAL to account for its non-rectangular shape
  //since the minimum x and y cell positions occur for EDGE blocks; let's set the low edge of the grid at xmin + block spacing/2 etc.
  double binxmin = xcellmin + 0.5*xspace;
  double binymin = ycellmin + 0.5*yspace;

  int nbinsx = 0, nbinsy = 0;
  while( binxmin + 3.0*xspace*nbinsx < xcellmax ) nbinsx++;
  while( binymin + 3.0*yspace*nbinsy < ycellmax ) nbinsy++; 

  double binxmax = binxmin + 3.0*xspace*nbinsx;
  double binymax = binymin + 3.0*yspace*nbinsy;
  
  //  double binxmax = xcellmax + 0.5*xspace;
  //double binymax = ycellmax + 0.5*yspace;

  //this rounds to the nearest integer number of bins and adds 1
  //int nbinsx = int( (binxmax-binxmin)/(2.0*xspace) + 0.5) + 1;
  //int nbinsy = int( (binymax-binymin)/(2.0*yspace) + 0.5) + 1;

  
  
  int nbinstot = nbinsx*nbinsy;
  
  //now implement 3x5 groupings and binning:
  for( auto cell : cells_ecal ){
    int row = rows_cells_ecal[cell];
    int col = cols_cells_ecal[cell];

    double xtemp = xcell_ecal[cell];
    double ytemp = ycell_ecal[cell];

    //int binxtemp = int( (xtemp - xcellmin + 0.5*xspace)/(2.0*xspace) );
    
    if( row > 0 && row + 1 < rows_ecal.size() &&
	col > 0 && col + 1 < columns_rows_ecal[row].size() ){ //not an edge block; make a cluster!
      list_of_nodes_ecal.insert(node); //"node" starts at ZERO and goes to Nnodes - 1
      cluster_centers_ecal_by_node[node] = cell;
      nodes_ecal_by_cluster_center[cell] = node;

      //here we are calculating the global position of the cell center coordinate in the g4sbs coordinate system:
      TVector3 clusterpos_global = ecal_pos + xtemp * ecal_xaxis + ytemp * ecal_yaxis;

      //Calculate the polar scattering angle for this cluster center assuming a point target at the origin
      double clustertheta = acos( clusterpos_global.Z() / clusterpos_global.Mag() );

      //Calculate the expected scattered electron energy for this scattering angle:
      double clusterEprime = Ebeam_default/(1.0 + Ebeam_default/Mp * (1.0-cos(clustertheta)));

      //store expected scattered electron energy for this cluster center:
      Eprime_expect_by_node_ecal[node] = clusterEprime; 

      //figure out the bins:
      int binxtemp = int( (xtemp - binxmin)/(3.0*xspace) );
      int binytemp = int( (ytemp - binymin)/(3.0*yspace) );
      int binglobaltemp = binxtemp + binytemp * nbinsx;
      if( binxtemp >= nbinsx ) cout << "warning: cluster center with bin x >= nbinsx, this will screw up global bin counting" << endl;

      list_of_bins_ecal.insert( binglobaltemp );

      binx_by_node_ecal[node] = binxtemp;
      biny_by_node_ecal[node] = binytemp;
      binglobal_by_node_ecal[node] = binglobaltemp;
      binx_by_cell_ecal[cell] = binxtemp;
      biny_by_cell_ecal[cell] = binytemp;
      binglobal_by_cell_ecal[cell] = binglobaltemp;

      logic_mean_ecal[node] = Eprime_expect_by_node_ecal[node]*normfac_ECAL;
      logic_sigma_ecal[node] = 0.06*logic_mean_ecal[node];
      
      cells_by_binglobal_ecal[binglobaltemp].insert( cell );
      nodes_by_binglobal_ecal[binglobaltemp].insert( node );
      
      //cells_logic_sums_ecal[node].insert( cell );
      
      for( int rowj=row-2; rowj<=row+2; rowj++ ){
	int colmin=col-1, colmax=col+1;
	if( rowj < row ){
	  colmin = nearestcol_down_ecal[cell]-1;
	  colmax = colmin + 2;
	} else if( rowj > row ){
	  colmin = nearestcol_up_ecal[cell]-1;
	  colmax = colmin + 2;
	}

	//check valid row and column indices:
	if( rowj >= 0 && rowj < rows_ecal.size() ){ 
	  for( int colj = colmin; colj <= colmax; colj++ ){
	    if( colj >= 0 && colj < columns_rows_ecal[rowj].size() ){
	      std::pair<int,int> rowcolj(rowj,colj);
	      auto testcell = cell_rowcol_ecal.find(rowcolj);

	      if( testcell != cell_rowcol_ecal.end() ){
		cells_logic_sums_ecal[node].insert( testcell->second );
		nodes_cells_ecal[cell].insert(node);
	      }
	    }
	  }
	}
      }
      
      node++;
    }
  }

  ofstream clustermap_ecal("database/ecal_cluster_mapping.txt");

  TString header = "# Format: cluster index, cluster center cell index, x center (m), y center (m), x bin number, y bin number, global bin number, total number of cells, list of cells";

  header.Form( "# Format: %15s, %25s, %15s, %15s, %15s, %15s, %20s, %25s, %25s, %15s", "cluster index", "cluster center cell index", "x center (m)", "y center (m)", "x bin number", "y bin number", "global bin number", "Elastic e- energy (GeV)","total number of cells", "list of cells" ); 
  
  clustermap_ecal << header.Data() << endl;
  for( auto node : list_of_nodes_ecal ){
    //TString thisnode;
    // clustermap_ecal << node << ", " << cluster_centers_ecal_by_node[node] << ", "
    // 		    << binx_by_node_ecal[node] << ", " << biny_by_node_ecal[node] << ", "
    // 		    << binglobal_by_node_ecal[node] << ", " << cells_logic_sums_ecal[node].size();

    TString line;
    line.Form( "          %15d, %25d, %15.6g, %15.6g, %15d, %15d, %20d, %25g, %25d", 
	       node, cluster_centers_ecal_by_node[node], xcell_ecal[cluster_centers_ecal_by_node[node]], ycell_ecal[cluster_centers_ecal_by_node[node]],
	       binx_by_node_ecal[node], biny_by_node_ecal[node], binglobal_by_node_ecal[node], Eprime_expect_by_node_ecal[node], cells_logic_sums_ecal[node].size() );
    
    clustermap_ecal << line;
    
    for( auto cell : cells_logic_sums_ecal[node] ){
      clustermap_ecal << ", " << cell;
    }

    clustermap_ecal << endl;
    
    //  clustermap_ecal << ", " << xcell_ecal[cluster_centers_ecal_by_node[node]]
    //		    << ", " << ycell_ecal[cluster_centers_ecal_by_node[node]] << endl; 
					 

  }

  //Also output cluster bins by 
  ofstream clusterbins_ecal("database/ecal_cluster_binning.txt");

  //Here we want the global bin number, x bin number, y bin number, X low, X high, Y low, Y high (or maybe just the centers?) number of clusters in the bin, and the list of clusters contained in that global bin:
  header.Form("# Format: %20s, %15s, %15s, %20s, %20s, %20s, %15s","global bin number", "X bin number", "Y bin number", "X bin center (m)", "Y bin center (m)", "Number of clusters", "List of clusters");

  clusterbins_ecal << header.Data() << endl;

  for( auto bin : list_of_bins_ecal ){
    TString line;

    int binx = bin%nbinsx;
    int biny = bin/nbinsx;
    int nclust = nodes_by_binglobal_ecal[bin].size();
    line.Form( "          %20d, %15d, %15d, %20g, %20g, %20d", bin, binx, biny,
	       binxmin + (3.0*binx + 1.0 )*xspace,
	       binymin + (3.0*biny + 1.0 )*yspace,
	       nclust );
    
    clusterbins_ecal << line; 
    for( auto clust : nodes_by_binglobal_ecal[bin] ){
      clusterbins_ecal << ", " << clust;
    }
    clusterbins_ecal << endl;
    
  }
  
  // int current_node = 1;

  // bool first_cell = true;

  // //Get trigger logic summing info:
  
  // while( currentline.ReadLine(logicfile_ecal) ){
  //   if( !currentline.BeginsWith( "#" ) ){
      
  //     TObjArray *tokens = currentline.Tokenize(" ");
  //     int ntokens = tokens->GetEntries();

  //     TString snode = ( (TObjString*) (*tokens)[0] )->GetString();
  //     int nodenumber = snode.Atoi() + 1;
      
  //     TString sncell_node = ( (TObjString*) (*tokens)[1] )->GetString();
  //     int ncell_node = sncell_node.Atoi();
      
  //     if( ntokens >= ncell_node + 2 ){
  // 	list_of_nodes_ecal.insert( nodenumber );
  // 	for( int itoken = 2; itoken < ncell_node+2; itoken++ ){
  // 	  TString scell = ( (TObjString*) (*tokens)[itoken] )->GetString();
  // 	  int cell = scell.Atoi();

  // 	  cells_logic_sums_ecal[nodenumber].insert( cell );

  // 	  //provide default values;
  // 	  logic_mean_ecal[nodenumber] = 3.0;
  // 	  logic_sigma_ecal[nodenumber] = 0.06*3.0;
  // 	  threshold_ecal[nodenumber] = 0.8;
	  
  // 	  nodes_cells_ecal[ cell ].insert( nodenumber );
  // 	}
  //     }
  //   }
  // }

  //For HCAL we want to read in HCAL_map.txt but only for the purpose of mapping cell indices to block positions. HCAL will always have exactly 24 rows x 12 columns; however, if the mapping between cell and row/column indices in g4sbs changes in the future, then significant parts of this code will need re-writing:

  ifstream HCAL_cell_map("database/HCAL_map.txt");
  map<int,double> xcell_hcal;
  map<int,double> ycell_hcal;
  
  //For HCAL we are going to use 3x3 sums:
  set<int> list_of_nodes_hcal;
  map<int, int> cluster_centers_hcal_by_node;
  map<int, int> nodes_hcal_by_cluster_center;
  map<int, set<int> > cells_logic_sums_hcal; //mapping between node numbers and cell numbers
  map<int, double> Tp_expect_hcal; //Proton KE by cluster center
  map<int, double> logic_mean_hcal; //mean peak positions by node number
  map<int, double> logic_sigma_hcal; //peak width by node number
  map<int, double> threshold_hcal; //threshold by node number
  map<std::pair<int,int>, int > cell_rowcol_hcal; //cell numbers mapped by unique row and column pairs
  map<int,set<int> > nodes_cells_hcal; //mapping of nodes by cell number:
  //For HCAL, the bins are also 2x2
  
  set<int>  list_of_bins_hcal;
  map<int,int> binx_by_cell_hcal;
  map<int,int> biny_by_cell_hcal;
  map<int,int> binglobal_by_cell_hcal;
  
  map<int,int> binx_by_node_hcal;
  map<int,int> biny_by_node_hcal;
  map<int,int> binglobal_by_node_hcal;
  map<int,set<int> > nodes_by_binglobal_hcal;
  map<int,double> xcellavg_by_binglobal_hcal;
  map<int,double> ycellavg_by_binglobal_hcal;
  
  while( currentline.ReadLine( HCAL_cell_map ) ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(",");
      int ntokens = tokens->GetEntries();
      if( ntokens == 5 ){
	TString scell = ( (TObjString*) (*tokens)[0] )->GetString();
	TString srow  = ( (TObjString*) (*tokens)[1] )->GetString();
	TString scol  = ( (TObjString*) (*tokens)[2] )->GetString();
	TString sxcell = ( (TObjString*) (*tokens)[3] )->GetString();
	TString sycell = ( (TObjString*) (*tokens)[4] )->GetString();

	int cell = scell.Atoi();
	int row = srow.Atoi();
	int col = scol.Atoi();
	double xcell = sxcell.Atof();
	double ycell = sycell.Atof();

	std::pair<int,int> rowcoltemp(row,col);
	cell_rowcol_hcal[rowcoltemp] = cell;
	xcell_hcal[cell] = xcell/100.0; //convert to meters
	ycell_hcal[cell] = ycell/100.0;
      }
    }
  }

  //Calculate "expected" proton K.E. for this node (assumes no deflection in CH2):
  TVector3 HCAL_zaxis(-sin(SBStheta_default),0,cos(SBStheta_default));
  TVector3 HCAL_yaxis(0,1,0);
  TVector3 HCAL_xaxis = HCAL_yaxis.Cross(HCAL_zaxis).Unit();

  TVector3 HCALpos = HCALdist_default * HCAL_zaxis + HCALvoff_default * HCAL_yaxis;
  
  int current_node = 0;
  //  bool first_cell = true;

  int nrows_hcal=24;
  int ncols_hcal=12;
  for( int row=1; row<nrows_hcal-1; row++ ){ //rows 1 to 22
    for( int col=1; col<ncols_hcal-1; col++ ){ //columns 1 to 10
      list_of_nodes_hcal.insert( current_node );

      int binx = (col-1)/2;
      int biny = (row-1)/2;
      int binglobal = binx + 5 * biny;
      
      for( int m=col-1; m<=col+1; m++ ){ //Add all blocks in each 3x3 sum to the current node:
	for( int n=row-1; n<=row+1; n++ ){ //Add
	  //in HCAL, cell = col + row*12, and column runs from left to right and row runs from top to bottom
	  //int cell = m + ncols_hcal * n;
	  int cell = cell_rowcol_hcal[std::make_pair(n,m)];
	  // std::cout << "Cell index (calculated, mapped)=("
	  // 	    << m+ncols_hcal*n << ", " << cell << ")" << std::endl;

	  cells_logic_sums_hcal[current_node].insert( cell );
	  nodes_cells_hcal[cell].insert(current_node);
	  //Now the following mapping is done while reading the HCAL map file from g4sbs
	  //std::pair<int,int> rowcol(n,m);
	  //cell_rowcol_hcal[rowcol] = cell;

	  if( m == col && n == row ){
	    cluster_centers_hcal_by_node[current_node] = cell;
	    nodes_hcal_by_cluster_center[cell] = current_node;

	    binx_by_node_hcal[current_node] = binx;
	    biny_by_node_hcal[current_node] = biny;
	    binglobal_by_node_hcal[current_node] = binglobal;

	    binx_by_cell_hcal[cell] = binx;
	    biny_by_cell_hcal[cell] = biny;
	    binglobal_by_cell_hcal[cell] = binglobal;

	    nodes_by_binglobal_hcal[binglobal].insert(current_node);

	    list_of_bins_hcal.insert( binglobal );
	  }
	  
	}
      }

      //Use some sensible default value for logic mean and logic sigma:

      //double Eprime_expect = Ebeam_default/(1.0+Ebeam_default/Mp*(1.0-cos(ECALtheta_default)));
      double xcelltemp = xcell_hcal[cluster_centers_hcal_by_node[current_node]];
      double ycelltemp = ycell_hcal[cluster_centers_hcal_by_node[current_node]];

      TVector3 cellpos = HCALpos + xcelltemp * HCAL_xaxis + ycelltemp * HCAL_yaxis;

      TVector3 cellpos_corr = cellpos - pdeflect_default * HCAL_yaxis;

      double ptheta_cell = acos( cellpos_corr.Z() / cellpos_corr.Mag() );

      double pp_ptheta = 2.0*Mp*Ebeam_default * (Mp + Ebeam_default)*cos(ptheta_cell) / ( pow(Mp,2) + 2.0*Mp*Ebeam_default + pow(Ebeam_default*sin(ptheta_cell),2) ); 
      
      double Tp_expect = sqrt( pow(pp_ptheta,2)+pow(Mp,2) ) - Mp;
      
      //do this in terms of energy, not photoelectrons:
      logic_mean_hcal[current_node] = 0.065*Tp_expect*normfac_HCAL;
      logic_sigma_hcal[current_node] = 0.5*logic_mean_hcal[current_node];
      Tp_expect_hcal[current_node] = Tp_expect;
      
      current_node++;
    }
  }

  

  //TODO: write out cluster mapping and binning info for HCAL as well:

  ofstream clustermap_hcal("database/hcal_cluster_mapping.txt");

  header = "# Format: node index, cluster center cell index, x bin number, y bin number, global bin number, total number of cells, list of cells";

  header.Form( "# Format: %15s, %25s, %15s, %15s, %15s, %15s, %20s, %25s, %25s, %25s, %20s, %15s", "cluster index", "cluster center cell index", "x center (m)", "y center (m)", "x bin number", "y bin number", "global bin number", "Elastic proton KE (GeV)", "edep in HCAL scint (GeV)","Sampling Fraction","total number of cells", "list of cells" );
  clustermap_hcal << header.Data() << endl;

  for ( auto node : list_of_nodes_hcal ){
    // clustermap_hcal << node << ", " << cluster_centers_hcal_by_node[node] << ", "
    // 		    << binx_by_node_hcal[node] << ", " << biny_by_node_hcal[node] << ", "
    // 		    << binglobal_by_node_hcal[node] << ", " << cells_logic_sums_hcal[node].size();

    TString line;
    line.Form( "          %15d, %25d, %15.6g, %15.6g, %15d, %15d, %20d, %25g, %25g, %25d",
	       node, cluster_centers_hcal_by_node[node], xcell_hcal[cluster_centers_hcal_by_node[node]], ycell_hcal[cluster_centers_hcal_by_node[node]],
	       binx_by_node_hcal[node], biny_by_node_hcal[node], binglobal_by_node_hcal[node], Tp_expect_hcal[node], logic_mean_hcal[node], cells_logic_sums_hcal[node].size() );
    clustermap_hcal << line;
    for( auto cell : cells_logic_sums_hcal[node] ){
      clustermap_hcal << ", " << cell;
    }
    clustermap_hcal << endl;
  }
  
  ofstream clusterbins_hcal("database/hcal_cluster_binning.txt");

  header.Form("# Format: %20s, %15s, %15s, %20s, %20s","global bin number", "X bin number", "Y bin number", "Number of clusters", "List of clusters");

  clusterbins_hcal << header.Data() << endl;

  for( auto bin : list_of_bins_hcal ){
    int binx = bin % 5;
    int biny = bin/5;
    int nclust = nodes_by_binglobal_hcal[bin].size();

    TString line;
    line.Form( "          %20d, %15d, %15d, %20d", bin, binx, biny, nclust); 

    clusterbins_hcal << line;

    for( auto clust : nodes_by_binglobal_hcal[bin] ){
      clusterbins_hcal << ", " << clust;
    }
    clusterbins_hcal << endl;
  }
	      
  TH1D::SetDefaultSumw2();

  //double PI = TMath::Pi();

  //Photoelectron statistics:
  double phe_per_GeV_ECAL = 1000.0/1.33; //~ 750 pe/GeV
  double phe_per_GeV_HCAL = 1000.0/0.30; //~ 3,333 pe/GeV (but sampling fraction is small)

  
  
  //read in alternate threshold:
  // ifstream thresholdfile_ecal(thresholdfilename_ecal);
  // if( thresholdfile_ecal ){
  //   int node;
  //   double mean,sigma;
  //   while( thresholdfile_ecal >> node >> mean >> sigma ){
  //     if( list_of_nodes_ecal.find( node ) != list_of_nodes_ecal.end() ){
  // 	logic_mean_ecal[ node ] = mean;
  // 	logic_sigma_ecal[ node ] = sigma;
  //     }
  //   }
  // }

  //read in alternate threshold:
  // ifstream thresholdfile_hcal(thresholdfilename_hcal);
  // if( thresholdfile_hcal ){
  //   int node;
  //   double mean,sigma;
  //   while( thresholdfile_hcal >> node >> mean >> sigma ){
  //     if( list_of_nodes_hcal.find( node ) != list_of_nodes_hcal.end() ){
  // 	logic_mean_hcal[ node ] = mean;
  //  	logic_sigma_hcal[ node ] = sigma;
  //     }
  //   }
  // }

  ifstream assocfile( assocfilename );

  //the following is not used:
  // bool use_ECAL_HCAL_associations=false;
  map<int, set<int> > ECAL_nodes_HCAL;
  if( assocfile ){
    while( !assocfile.eof() ){

      TString currentline;
      while( currentline.ReadLine(assocfile) ){
	if( !currentline.BeginsWith("#") ) {
	  TObjArray *tokens = currentline.Tokenize(",");
	  if( tokens->GetEntries() > 2 ){
	    TString shcalbin = ( (TObjString*) (*tokens)[0] )->GetString();
	    TString sNecalbins = ( (TObjString*) (*tokens)[1] )->GetString();

	    int hcalbin = shcalbin.Atoi();
	    int N = sNecalbins.Atoi();

	    if( tokens->GetEntries() >= N+2 ){
	      for( int i=0; i<N; i++ ){
		TString secalbin = ( (TObjString*) (*tokens)[i+2] )->GetString();
		int ecalbin = secalbin.Atoi();

		ECAL_nodes_HCAL[hcalbin].insert(ecalbin);
	      }
	    }
	  }
	  
	}
      }
    }
  }
  
  fout->cd();

  // TH1D *hrate_vs_threshold_ECAL = new TH1D("hrate_vs_threshold_ECAL","",30,0.0,1.5);
  // //TH1D *hnum_logic_sums_fired_vs_threshold = new TH1D("hnum_logic_sums_fired_vs_threshold

  // TH1D *hrate_vs_threshold_HCAL = new TH1D("hrate_vs_threshold_HCAL","",40,0.0,2.0);

  //TH2D *htrue_coincidence_rate_vs_threshold_ECAL_HCAL = new TH2D("htrue_coincidence_rate_vs_threshold_ECAL_HCAL","",40,0,2.0,30,0,1.5);

  TH1D *hnclust_ECAL = new TH1D("hnclust_ECAL","ECAL ; number of clusters ; ", 101,-0.5,100.5);
  TH1D *hnclust_HCAL = new TH1D("hnclust_HCAL","HCAL ; number of clusters ; ", 101,-0.5,100.5);

  TH1D *hEsum_smear_ECALmax = new TH1D("hEsum_smear_ECALmax","ECAL ; smeared E of max cluster (GeV);",3000,0.0,6.0);
  TH1D *hEsum_true_ECALmax = new TH1D("hEsum_true_ECALmax","ECAL ; true E of max cluster (GeV);",3000,0.0,6.0);
  TH1D *hnphesum_ECALmax = new TH1D("hnphesum_ECALmax","ECAL ; num photoelectrons;",3000,0.0,3000);
  
  TH1D *hEsum_smear_HCALmax_all = new TH1D("hEsum_smear_HCALmax_all","HCAL (all); smeared E of max cluster (GeV);",200,0,1.0);
  TH1D *hEsum_true_HCALmax_all = new TH1D("hEsum_true_HCALmax_all", "HCAL (all); true E of max cluster (GeV);",200,0,1.0);
  TH1D *hnphesum_HCALmax_all = new TH1D("hnphesum_HCALmax_all", "HCAL (all); num photoelectrons;", 2000,0,4000);

  TH1D *hEsum_smear_HCALmax_FTcut = new TH1D("hEsum_smear_HCALmax_FTcut","HCAL (FT); smeared E of max cluster (GeV);",200,0,1.0);
  TH1D *hEsum_true_HCALmax_FTcut = new TH1D("hEsum_true_HCALmax_FTcut", "HCAL (FT); true E of max cluster (GeV);",200,0,1.0);
  TH1D *hnphesum_HCALmax_FTcut = new TH1D("hnphesum_HCALmax_FTcut", "HCAL (FT); num photoelectrons;", 2000,0,4000);

  TH1D *hEsum_smear_HCALmax_FPP1cut = new TH1D("hEsum_smear_HCALmax_FPP1cut","HCAL (FPP1); smeared E of max cluster (GeV);",200,0,1.0);
  TH1D *hEsum_true_HCALmax_FPP1cut = new TH1D("hEsum_true_HCALmax_FPP1cut", "HCAL (FPP1); true E of max cluster (GeV);",200,0,1.0);
  TH1D *hnphesum_HCALmax_FPP1cut = new TH1D("hnphesum_HCALmax_FPP1cut", "HCAL (FPP1); num photoelectrons;", 2000,0,4000);
  
  TH2D *hnphesum_vs_node_ECAL_all = new TH2D("hnphesum_vs_node_ECAL_all","",list_of_nodes_ecal.size()+1,-0.5,list_of_nodes_ecal.size()+0.5,100,0.0,5000.0);
  TH2D *hnphesum_vs_node_ECAL_max = new TH2D("hnphesum_vs_node_ECAL_max","",list_of_nodes_ecal.size()+1,-0.5,list_of_nodes_ecal.size()+0.5,100,0.0,5000.0);

  TH2D *hEsum_smear_vs_node_ECAL_all = new TH2D("hEsum_smear_vs_node_ECAL_all","",list_of_nodes_ecal.size()+1,-0.5,list_of_nodes_ecal.size()+0.5,200,0,10.0);
  TH2D *hEsum_smear_vs_node_ECAL_max = new TH2D("hEsum_smear_vs_node_ECAL_max","",list_of_nodes_ecal.size()+1,-0.5,list_of_nodes_ecal.size()+0.5,200,0,10.0);

  TH2D *hEsum_true_vs_node_ECAL_all = new TH2D("hEsum_true_vs_node_ECAL_all","",list_of_nodes_ecal.size()+1,-0.5,list_of_nodes_ecal.size()+0.5,200,0,10.0);
  TH2D *hEsum_true_vs_node_ECAL_max = new TH2D("hEsum_true_vs_node_ECAL_max","",list_of_nodes_ecal.size()+1,-0.5,list_of_nodes_ecal.size()+0.5,200,0,10.0);

  TH1D *hEsum_true_norm_ECAL_max = new TH1D("hEsum_true_norm_ECAL_max","All events; True cluster energy/expected energy;", 2000,0.0,2.0);
  TH1D *hEsum_smear_norm_ECAL_max = new TH1D("hEsum_smear_norm_ECAL_max","All events; Smeared cluster energy/expected energy;",2000,0.0,2.0);
  TH1D *hEsum_smear_norm_ECAL_max_shouldhit = new TH1D("hEsum_smear_norm_ECAL_max_shouldhit","Should hit ECAL; Smeared cluster energy/expected energy;",2000,0.0,2.0);

  TH1D *hEsum_true_norm_HCAL_max_FPP1cut = new TH1D("hEsum_true_norm_HCAL_max_FPP1cut", "Good FPP1 ; HCAL true energy/expected;",2000,0.0,2.0);
  TH1D *hEsum_smear_norm_HCAL_max_FPP1cut = new TH1D("hEsum_smear_norm_HCAL_max_FPP1cut", "Good FPP1 ; HCAL smeared energy/expected;",2000,0.0,2.0);

  TH2D *hEsum_smear_norm_vs_node_ECAL_max = new TH2D("hEsum_smear_norm_vs_node_ECAL_max","All events; ECAL node; ECAL smeared energy/expected",list_of_nodes_ecal.size()+1,-0.5,list_of_nodes_ecal.size()+0.5,200,0,2);
  TH2D *hEsum_smear_norm_vs_node_ECAL_max_shouldhit = new TH2D("hEsum_smear_norm_vs_node_ECAL_max_shouldhit","Should hit ECAL; ECAL node; ECAL smeared energy/expected",list_of_nodes_ecal.size()+1,-0.5,list_of_nodes_ecal.size()+0.5,200,0,2);

  TH2D *hEsum_smear_norm_vs_node_HCAL_max_FPP1cut = new TH2D("hEsum_smear_norm_vs_node_HCAL_max_FPP1cut","FPP1 cut; HCAL node; HCAL cluster smeared energy/expected",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,200,0,2);
  
  TH2D *hnphesum_vs_node_HCAL_all = new TH2D("hnphesum_vs_node_HCAL_all","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  
  //TH2D *hnphesum_vs_node_ECAL_FTcut = new TH2D("hnphesum_vs_node_ECAL_FTcut","",list_of_nodes_ecal.size(),0.5,list_of_nodes_ecal.size()+0.5,100,0.0,5000.0);
  TH2D *hnphesum_vs_node_HCAL_FTcut = new TH2D("hnphesum_vs_node_HCAL_FTcut","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCAL_FPP1cut = new TH2D("hnphesum_vs_node_HCAL_FPP1cut","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  //TH2D *hnphesum_vs_node_HCAL_FPP2cut = new TH2D("hnphesum_vs_node_HCAL_FPP2cut","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  //TH2D *hnphesum_vs_node_HCAL_FPPbothcut = new TH2D("hnphesum_vs_node_HCAL_FPPbothcut","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  //TH2D *hnphesum_vs_node_HCAL_FPPeithercut = new TH2D("hnphesum_vs_node_HCAL_FPPeithercut","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);

  TH2D *hnphesum_vs_node_HCALmax_all = new TH2D("hnphesum_vs_node_HCALmax_all","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCALmax_FTcut = new TH2D("hnphesum_vs_node_HCALmax_FTcut","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  TH2D *hnphesum_vs_node_HCALmax_FPP1cut = new TH2D("hnphesum_vs_node_HCALmax_FPP1cut","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  //TH2D *hnphesum_vs_node_HCALmax_FPP2cut = new TH2D("hnphesum_vs_node_HCALmax_FPP2cut","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  //TH2D *hnphesum_vs_node_HCALmax_FPPbothcut = new TH2D("hnphesum_vs_node_HCALmax_FPPbothcut","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  //TH2D *hnphesum_vs_node_HCALmax_FPPeithercut = new TH2D("hnphesum_vs_node_HCALmax_FPPeithercut","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,100,0.0,3500.0);
  

  TH2D *hEsum_true_vs_node_HCALmax_all = new TH2D("hEsum_true_vs_node_HCALmax_all","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,200,0,2.0);
  TH2D *hEsum_true_vs_node_HCALmax_FTcut = new TH2D("hEsum_true_vs_node_HCALmax_FTcut","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,200,0,2.0);
  TH2D *hEsum_true_vs_node_HCALmax_FPP1cut = new TH2D("hEsum_true_vs_node_HCALmax_FPP1cut","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,200,0,2.0);

  TH2D *hEsum_smear_vs_node_HCALmax_all = new TH2D("hEsum_smear_vs_node_HCALmax_all","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,200,0,2.0);
  TH2D *hEsum_smear_vs_node_HCALmax_FTcut = new TH2D("hEsum_smear_vs_node_HCALmax_FTcut","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,200,0,2.0);
  TH2D *hEsum_smear_vs_node_HCALmax_FPP1cut = new TH2D("hEsum_smear_vs_node_HCALmax_FPP1cut","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,200,0,2.0);
  
  
  
  TH2D *hmaxnode_ECAL_vs_HCAL = new TH2D("hmaxnode_ECAL_vs_HCAL","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,list_of_nodes_ecal.size()+1,-0.5,list_of_nodes_ecal.size()+0.5);
  TH2D *hallnodes_ECAL_vs_HCAL = new TH2D("hallnodes_ECAL_vs_HCAL","",list_of_nodes_hcal.size()+1,-0.5,list_of_nodes_hcal.size()+0.5,list_of_nodes_ecal.size()+1,-0.5,list_of_nodes_ecal.size()+0.5);

  
  auto firstbin_ecal = *(list_of_bins_ecal.begin());
  auto lastbin_ecal = *(list_of_bins_ecal.rbegin());
  auto firstbin_hcal = *(list_of_bins_hcal.begin());
  auto lastbin_hcal = *(list_of_bins_hcal.rbegin());

  int necal=lastbin_ecal-firstbin_ecal+1;
  int nhcal=lastbin_hcal-firstbin_hcal+1;
  
  TH2D *hmaxbin_ECAL_vs_HCAL = new TH2D("hmaxbin_ECAL_vs_HCAL","",nhcal,firstbin_hcal-0.5,lastbin_hcal+0.5,necal,firstbin_ecal-0.5,lastbin_ecal+0.5);
  TH2D *hallbins_ECAL_vs_HCAL = new TH2D("hallbins_ECAL_vs_HCAL","",nhcal,firstbin_hcal-0.5,lastbin_hcal+0.5, necal, firstbin_ecal-0.5,lastbin_ecal+0.5);
  
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

  TH1D *hthetaFPP1_shouldhit = new TH1D("hthetaFPP1_shouldhit","",120,0.0,12.0);
  
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
  
  
  double Ibeam = 70.0e-6; //Amps
  double Ltarget = 30.0; //cm
  double e = 1.602e-19; //electron charge;
  double rho_target = 0.072; //g/cm^3
  double N_A = 6.022e23; //atoms/mol:
  double Mmol_H = 1.008; //g/mol
  double Lumi = rho_target * Ltarget * N_A / Mmol_H * Ibeam/e; //~ 5.6e38;
  
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
    weight = T->ev_rate / double(nfiles); //wrong, but not important for our purposes here:
    
    
    TString fnametemp = C->GetFile()->GetName();

    if( bad_file_list.find( fnametemp ) == bad_file_list.end() ){

      double R = ECALdist_file[fnametemp];
      double thetacal = ECALtheta_file[fnametemp];
      double Ebeam = Ebeam_file[fnametemp];

      weight = T->ev_sigma * Lumi_file[fnametemp] * GenVol_file[fnametemp] / double(ngen);
      
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

      TVector3 nhat_FT, nhat_FPP1;
      if( FTtrack ){
	nhat_FT.SetXYZ( (*(T->Harm_FT_Track_Xp))[itrack_FT],
			(*(T->Harm_FT_Track_Yp))[itrack_FT],
			1.0 );
	nhat_FT = nhat_FT.Unit();
      }

      double thetaFPP1, pFPP1;
      bool FPP1track = false;
      //    if( FTtrack )

      int itrack_FPP1=-1;
      if( T->Harm_FPP1_Track_ntracks > 0 && FTtrack ){
	double thetamin=PI+0.001;
	for( int itrack=0; itrack<T->Harm_FPP1_Track_ntracks; itrack++ ){
	  if( (*(T->Harm_FPP1_Track_MID))[itrack] == 0 ){ //should we insist on a primary proton track or should we allow more?
	    nhat_FPP1.SetXYZ( (*(T->Harm_FPP1_Track_Xp))[itrack],
			      (*(T->Harm_FPP1_Track_Yp))[itrack],
			      1.0 );
	    nhat_FPP1 = nhat_FPP1.Unit();
	    
	    thetaFPP1 = acos( nhat_FPP1.Dot( nhat_FT ) );	  
	    
	    pFPP1 = (*(T->Harm_FPP1_Track_P))[itrack];
	    
	    if( itrack == 0 || thetaFPP1 < thetamin ){
	      thetamin = thetaFPP1;
	      itrack_FPP1 = itrack;
	    }
	    
	    FPP1track = true;
	  }
	}
	if( FPP1track ) hthetaFPP1->Fill(thetaFPP1*180.0/PI,weight);
      }

      //no FPP2 tracks in final plan:
      // if( T->Harm_FPP2_Track_ntracks > 0 && FTtrack && FPP1track){
      // 	for( int itrack=0; itrack<T->Harm_FPP2_Track_ntracks; itrack++ ){
      // 	  if( (*(T->Harm_FPP2_Track_MID))[itrack] == 0 ){
      // 	    nhat_FPP2.SetXYZ( (*(T->Harm_FPP2_Track_Xp))[itrack],
      // 			      (*(T->Harm_FPP2_Track_Yp))[itrack],
      // 			      1.0 );
      // 	    nhat_FPP2 = nhat_FPP2.Unit();
      // 	    thetaFPP2 = acos( nhat_FPP2.Dot( nhat_FPP1 ) );
      // 	    pFPP2 = (*(T->Harm_FPP2_Track_P))[itrack];
      // 	    //FPP2track = thetaFPP2 < 24.0*PI/180.0 && pFPP2/T->ev_np > 0.5;
      // 	    FPP2track = thetaFPP2 < 12.0*PI/180.0 && pFPP2/T->ev_np > 0.5;
      // 	    if( FPP2track ) hthetaFPP2->Fill(thetaFPP2*180.0/PI,weight);
      // 	  }
      // 	}
      // }
      
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
    
      // map<int,double> node_sums; //initialize all node sums to zero:
      // map<int,double> node_Esum_true; 
      // map<int,double> node_Esum_smear; 
      // for( set<int>::iterator inode = list_of_nodes_ecal.begin(); inode != list_of_nodes_ecal.end(); ++inode ){
      // 	node_sums[ *inode ] = 0.0;
      // 	node_Esum_true[ *inode ] = 0.0;
      // 	node_Esum_smear[ *inode ] = 0.0;
      // }

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

      if( should_hit_ECAL && FTtrack && FPP1track ){
	hthetaFPP1_shouldhit->Fill( thetaFPP1*180.0/PI, weight );
      }
      
      pheflag = 0;
      
      //int nphe = 0;
      double nphe = 0.0;

      int nhits_ECAL = T->Earm_ECalTF1_hit_nhits;

      vector<calblock> ECALblocks_all;
      vector<calblock> ECALblocks_unused;
      
      // if( pheflag == 0 ){
      //We assume we aren't doing full optical photon simulations here:
      for( int ihit = 0; ihit<T->Earm_ECalTF1_hit_nhits; ihit++ ){
	int rowhit = ( *(T->Earm_ECalTF1_hit_row))[ihit];
	int colhit = ( *(T->Earm_ECalTF1_hit_col))[ihit];
	std::pair<int,int> rowcolhit( rowhit,colhit );
	
	//int cellhit = cell_rowcol_ecal[rowcolhit];
	int cellhit = (*(T->Earm_ECalTF1_hit_cell))[ihit];
	  
	//int trigger_group = nodes_cells_ecal[cellhit];
	
	double edep = (*(T->Earm_ECalTF1_hit_sumedep))[ihit];

	//double mean = 752.2*edep;
	//double sigma = 52.0*sqrt(edep) + 20.76*edep;

	double mean = 528.*edep; //from C16_GEANT_report.pdf, 528 phe/GeV
	nphe = num.Gaus(mean,sqrt(mean));

	double esmear = nphe/528.;

	calblock thisblock;
	thisblock.cell = cellhit;
	thisblock.row = rowhit;
	thisblock.col = colhit;
	thisblock.x = T->Earm_ECalTF1_hit_xcell->at(ihit);
	thisblock.y = T->Earm_ECalTF1_hit_ycell->at(ihit);
	thisblock.E = esmear;
	thisblock.Etrue = edep;
	thisblock.nphe = nphe;
	
	ECALblocks_all.push_back( thisblock );
	ECALblocks_unused.push_back( thisblock );
	
	//nphe = TMath::Max(0,TMath::Nint(num.Gaus(mean,sigma)));
	
	// for( set<int>::iterator inode = nodes_cells_ecal[cellhit].begin(); inode != nodes_cells_ecal[cellhit].end(); ++inode ){
	//   node_sums[ *inode ] += nphe;
	//   node_Esum_true[ *inode ] += edep;
	//   node_Esum_smear[ *inode ] += esmear; 
	// }
	
      }

      vector<calcluster> ECALclusters;
      
      while( !ECALblocks_unused.empty() ){
	// Sort remaining unused ECAL blocks by energy:
	std::sort( ECALblocks_unused.begin(), ECALblocks_unused.end(), CompareCalBlocks ); //This should sort in descending order by energy:

	//now form the clusters around all local maxima:

	auto maxblk = ECALblocks_unused.begin();

	//calblock max = *calblock_iterator;

	//Only form a cluser if the current local maximum is in the list of possible cluster centers; i.e., if it is not an edge block:
	if( nodes_ecal_by_cluster_center.find((*maxblk).cell) != nodes_ecal_by_cluster_center.end() ){
	  
	  calcluster thiscluster;
	  thiscluster.cell = (*maxblk).cell;
	  thiscluster.node = nodes_ecal_by_cluster_center[(*maxblk).cell];
	  thiscluster.bin = binglobal_by_cell_ecal[(*maxblk).cell];
	  thiscluster.binx = binx_by_cell_ecal[(*maxblk).cell];
	  thiscluster.biny = biny_by_cell_ecal[(*maxblk).cell];
	  thiscluster.Esum = 0.0;
	  thiscluster.Esum_true = 0.0;
	  thiscluster.nphesum = 0.0;
	  
	  //Next we need to loop on all blocks, add up the energies from any blocks in this cluster, and erase unused blocks in this cluster
	  // so they won't also form clusters:
	  ECALblocks_unused.erase( maxblk );

	  set<int> blockstoerase;
	  for( int iblk=0; iblk<ECALblocks_all.size(); iblk++ ){
	    calblock blk_i = ECALblocks_all[iblk];
	    if( cells_logic_sums_ecal[thiscluster.node].find( blk_i.cell ) != cells_logic_sums_ecal[thiscluster.node].end() ){
	      blockstoerase.insert(blk_i.cell);
	      //don't double-count the energy of the seed:
	      thiscluster.Esum += blk_i.E;
	      thiscluster.Esum_true += blk_i.Etrue;
	      thiscluster.nphesum += blk_i.nphe;
	    }
	  }

	  auto itblk = ECALblocks_unused.begin();

	  while( itblk != ECALblocks_unused.end() ){
	    if( blockstoerase.find((*itblk).cell) != blockstoerase.end() ){
	      itblk = ECALblocks_unused.erase( itblk );
	    } else {
	      ++itblk;
	    }
	  }
	  thiscluster.Esum_norm = thiscluster.Esum / logic_mean_ecal[thiscluster.node];
	  ECALclusters.push_back( thiscluster );
	} else { // current maximum is an edge block; erase without forming a cluster:
	  ECALblocks_unused.erase( maxblk );
	}
      } //done with ECAL clustering loop;

      //Sort ECAL clusters by energy (existing comparison operator is based on smeared energy sum, NOT normalized):
      std::sort( ECALclusters.begin(), ECALclusters.end(), CompareCalClusters );

      vector<int> trigger_nodes_fired(hefficiency_vs_threshold_ECAL->GetNbinsX(),0);
      

      int maxnode_ECAL=-1;
      int maxnode_HCAL=-1;
      double maxsum_ECAL = 0.0;
      double maxsum_HCAL = 0.0;

      bool ECALtrig_nominal = false;
      bool HCALtrig_nominal = false;
      
      int nominal_threshold_bin_HCAL = hefficiency_vs_threshold_HCAL_FTcut->FindBin(nominal_threshold_HCAL);
      int nominal_threshold_bin_ECAL = hefficiency_vs_threshold_ECAL->FindBin(nominal_threshold_ECAL);

      hnclust_ECAL->Fill( ECALclusters.size(), weight );
      
      if( ECALclusters.size() > 0 ){
	maxnode_ECAL = ECALclusters[0].node;
	//For THIS purpose (elastic events analysis), I think using the highest smeared energy cluster is as good as (or better) than trying to use the normalized energy:
	maxsum_ECAL = ECALclusters[0].Esum_norm; 

	hEsum_smear_ECALmax->Fill(ECALclusters[0].Esum, weight );
	hEsum_true_ECALmax->Fill( ECALclusters[0].Esum_true, weight );
	hnphesum_ECALmax->Fill( ECALclusters[0].nphesum, weight );

	hnphesum_vs_node_ECAL_max->Fill( maxnode_ECAL, ECALclusters[0].nphesum, weight );
	hEsum_smear_vs_node_ECAL_max->Fill( maxnode_ECAL, ECALclusters[0].Esum, weight );
	hEsum_true_vs_node_ECAL_max->Fill( maxnode_ECAL, ECALclusters[0].Esum_true, weight );

	hEsum_true_norm_ECAL_max->Fill( ECALclusters[0].Esum_true / logic_mean_ecal[maxnode_ECAL], weight );

       
	hEsum_smear_norm_ECAL_max->Fill( maxsum_ECAL, weight );
	hEsum_smear_norm_vs_node_ECAL_max->Fill( maxnode_ECAL, maxsum_ECAL, weight );
	if( should_hit_ECAL ){
	  hEsum_smear_norm_ECAL_max_shouldhit->Fill( maxsum_ECAL, weight );
	  hEsum_smear_norm_vs_node_ECAL_max_shouldhit->Fill( maxnode_ECAL, maxsum_ECAL, weight );
	}
	// This could be made faster by making a differential (rate versus energy) histogram rather than
	// accumulating an integrated rate above threshold at this level of granularity:
	//for( int bin=1; bin<=hefficiency_vs_threshold_ECAL->GetNbinsX(); bin++ ){
	// trigger_nodes_fired[bin-1] = 0;
	for( auto clust : ECALclusters ){
	  int triggerbin = int( hefficiency_vs_threshold_ECAL->FindBin( clust.Esum_norm ) );
	  if( clust.Esum_norm < hefficiency_vs_threshold_ECAL->GetBinCenter(triggerbin) ) triggerbin--;
	   
	  //  if( triggerbin >= nominal_threshold_bin_ECAL ) ECALtrig_nominal = true;
	  if( clust.Esum_norm >= nominal_threshold_ECAL ) ECALtrig_nominal = true;
	  for( int bin=1; bin<=triggerbin; bin++ ){
	    if( bin<=hefficiency_vs_threshold_ECAL->GetNbinsX() ) trigger_nodes_fired[bin-1]++;
	  }

	  hEsum_smear_vs_node_ECAL_all->Fill( clust.node, clust.Esum, weight );
	  //if( node_Esum_smear[*inode] > 0.0 ) {
	  hnphesum_vs_node_ECAL_all->Fill( clust.node, clust.nphesum, weight );
	  hEsum_true_vs_node_ECAL_all->Fill( clust.node, clust.Esum_true, weight );
	  //hEsum_smear_vs_node_ECAL_all->Fill( *inode, node_Esum_smear[*inode], weight );
	}
	
     
      }
      
      if( should_hit_ECAL ){
	for( int ithr=1; ithr<=hefficiency_vs_threshold_ECAL->GetNbinsX(); ithr++ ){
	  hshouldhit_vs_threshold_ECAL->Fill( hefficiency_vs_threshold_ECAL->GetBinCenter(ithr), weight );
	  if( trigger_nodes_fired[ithr-1] > 0 ){
	    hefficiency_vs_threshold_ECAL->Fill( hefficiency_vs_threshold_ECAL->GetBinCenter(ithr), weight );
	  }
	  if( FTtrack ){
	    hshouldhit_vs_threshold_ECAL_FTcut->Fill( hefficiency_vs_threshold_ECAL->GetBinCenter(ithr), weight );
	    if( trigger_nodes_fired[ithr-1] > 0 ){
	      hefficiency_vs_threshold_ECAL_FTcut->Fill( hefficiency_vs_threshold_ECAL->GetBinCenter(ithr), weight );
	    }
	  }
	}
      }

      if( FTtrack && should_hit_ECAL ){
	hshouldhit_vs_Q2_ECAL_FTcut->Fill( T->ev_Q2, weight );
	if( ECALtrig_nominal ){
	  hefficiency_vs_Q2_ECAL_FTcut->Fill( T->ev_Q2, weight );
	}
      }
    
      vector<calblock> HCALblocks_all;
      vector<calblock> HCALblocks_unused;
      
      //int nphe = 0;
    
      //if( pheflag == 0 ){

      double maxedep = 0.0;
      double maxphe = 0.0;
	
      for( int ihit=0; ihit<T->Harm_HCalScint_hit_nhits; ihit++ ){
	int rowhit = (*(T->Harm_HCalScint_hit_row))[ihit]; //why are we adding 1? We don't do that when we form the clusters and bins!
	int colhit = (*(T->Harm_HCalScint_hit_col))[ihit]; 
	std::pair<int,int> rowcolhit(rowhit,colhit);
	int cellhit = cell_rowcol_hcal[rowcolhit];
	//int trigger_group = nodes_cells_hcal[cellhit];
	double edep = (*(T->Harm_HCalScint_hit_sumedep))[ihit];
	//nphe = num.Poisson( phe_per_GeV_HCAL * edep );
	double mean = 2981.0*edep;
	//double sigma = 69.54*sqrt(edep) + 155.3*edep;

	//nphe = TMath::Max(0,TMath::Nint(num.Gaus(mean,sigma)));

	nphe = num.Gaus(mean,sqrt(mean));

	double esmear = nphe/2981.0;
	
	//cout << "HCAL hit " << ihit+1 << " node, edep, nphe = " << trigger_group << ", " << edep << ", " << nphe << endl;
	//node_sums_hcal[trigger_group] += double(nphe);
	calblock thisblock;
	thisblock.cell = cellhit;
	thisblock.row = rowhit;
	thisblock.col = colhit;
	thisblock.x = T->Harm_HCalScint_hit_xcell->at(ihit);
	thisblock.y = T->Harm_HCalScint_hit_ycell->at(ihit);
	thisblock.E = esmear;
	thisblock.Etrue = edep;
	thisblock.nphe = nphe;
	
	HCALblocks_all.push_back( thisblock );
	HCALblocks_unused.push_back( thisblock );
	
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
      
      vector<calcluster> HCALclusters;

      while( !HCALblocks_unused.empty() ){
	std::sort( HCALblocks_unused.begin(), HCALblocks_unused.end(), CompareCalBlocks );

	auto maxblk = HCALblocks_unused.begin();

	//Here again we have to add a check that we aren't looking at an edge block for the current local maximum,
	// it needs to be in the defined list of cluster centers:

	if( nodes_hcal_by_cluster_center.find((*maxblk).cell) != nodes_hcal_by_cluster_center.end() ){
	  
	  calcluster thiscluster;
	  
	  thiscluster.cell = (*maxblk).cell;
	  thiscluster.node = nodes_hcal_by_cluster_center[(*maxblk).cell];
	  thiscluster.bin = binglobal_by_cell_hcal[(*maxblk).cell];
	  thiscluster.binx = binx_by_cell_hcal[(*maxblk).cell];
	  thiscluster.biny = biny_by_cell_hcal[(*maxblk).cell];
	  thiscluster.Esum = 0.0;
	  thiscluster.nphesum = 0.0;
	  thiscluster.Esum_true = 0.0;
	  
	  HCALblocks_unused.erase( maxblk );
	  
	  set<int> blockstoerase;
	  for( int iblk=0; iblk<HCALblocks_all.size(); iblk++ ){
	    calblock blk_i = HCALblocks_all[iblk];
	    if( cells_logic_sums_hcal[thiscluster.node].find( blk_i.cell ) != cells_logic_sums_hcal[thiscluster.node].end() ){
	      blockstoerase.insert( blk_i.cell );
	      thiscluster.Esum += blk_i.E;
	      thiscluster.Esum_true += blk_i.Etrue;
	      thiscluster.nphesum += blk_i.nphe;
	    }
	    
	  }
	  
	  auto itblk = HCALblocks_unused.begin();
	  
	  while( itblk != HCALblocks_unused.end() ){
	    if( blockstoerase.find((*itblk).cell) != blockstoerase.end() ){
	      itblk = HCALblocks_unused.erase( itblk );
	    } else {
	      ++itblk;
	    }
	  }
	  
	  thiscluster.Esum_norm = thiscluster.Esum/logic_mean_hcal[thiscluster.node];
	  HCALclusters.push_back( thiscluster );
	} else {
	  HCALblocks_unused.erase( maxblk );
	}
      } // done with HCAL clustering loop:

      std::sort( HCALclusters.begin(), HCALclusters.end(), CompareCalClusters );
      
      
      vector<int> trigger_nodes_fired_hcal(hefficiency_vs_threshold_HCAL_FTcut->GetNbinsX(),0);
      

      vector<int> coin_trigger_fired( hefficiency_vs_threshold_HCAL_FTcut->GetNbinsX()*hefficiency_vs_threshold_ECAL->GetNbinsX(),0 );
 

      bool cointrig_nominal_threshold = false;

      hnclust_HCAL->Fill( HCALclusters.size(), weight );
      
      if( HCALclusters.size() > 0 ){
	maxnode_HCAL = HCALclusters[0].node;
	maxsum_HCAL = HCALclusters[0].Esum;

	hEsum_smear_HCALmax_all->Fill( HCALclusters[0].Esum, weight );
	hEsum_true_HCALmax_all->Fill( HCALclusters[0].Esum_true, weight );
	hnphesum_HCALmax_all->Fill( HCALclusters[0].nphesum, weight );
	if( FTtrack ) {
	  hEsum_smear_HCALmax_FTcut->Fill( HCALclusters[0].Esum, weight );
	  hEsum_true_HCALmax_FTcut->Fill( HCALclusters[0].Esum_true, weight );
	  hnphesum_HCALmax_FTcut->Fill( HCALclusters[0].nphesum, weight );
	  if( FPP1track ){
	    hEsum_smear_HCALmax_FPP1cut->Fill( HCALclusters[0].Esum, weight );
	    hEsum_true_HCALmax_FPP1cut->Fill( HCALclusters[0].Esum_true, weight );
	    hnphesum_HCALmax_FPP1cut->Fill( HCALclusters[0].nphesum, weight );

	    hEsum_true_norm_HCAL_max_FPP1cut->Fill( HCALclusters[0].Esum_true/logic_mean_hcal[maxnode_HCAL], weight );
	    hEsum_smear_norm_HCAL_max_FPP1cut->Fill( HCALclusters[0].Esum_norm, weight );
	    hEsum_smear_norm_vs_node_HCAL_max_FPP1cut->Fill( maxnode_HCAL, HCALclusters[0].Esum_norm, weight );
	  }
	}
	
	//for( int bin=1; bin<=hefficiency_vs_threshold_HCAL_FTcut->GetNbinsX(); bin++ ){
	// trigger_nodes_fired_hcal[bin-1] = 0;
	for( auto clust : HCALclusters ){

	  int triggerbin = int( hefficiency_vs_threshold_HCAL_FTcut->FindBin( clust.Esum_norm ) );
	  if( clust.Esum_norm < hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(triggerbin) ) triggerbin--;

	  // if( triggerbin >= nominal_threshold_bin_HCAL ) HCALtrig_nominal = true;
	  if( clust.Esum_norm >= nominal_threshold_HCAL ) HCALtrig_nominal = true;
	  for( int bin=1; bin<=triggerbin; bin++ ){
	    if( bin <= hefficiency_vs_threshold_HCAL_FTcut->GetNbinsX() ) trigger_nodes_fired_hcal[bin-1]++;
	  }
	}

	calcluster clust = HCALclusters[0];
	hEsum_smear_vs_node_HCALmax_all->Fill( clust.node, clust.Esum, weight );
	hEsum_true_vs_node_HCALmax_all->Fill( clust.node, clust.Esum_true, weight );
	hnphesum_vs_node_HCALmax_all->Fill( clust.node, clust.nphesum, weight );
	if( FTtrack ){
	  hEsum_smear_vs_node_HCALmax_FTcut->Fill( clust.node, clust.Esum, weight );
	  hEsum_true_vs_node_HCALmax_FTcut->Fill( clust.node, clust.Esum_true, weight );
	  hnphesum_vs_node_HCALmax_FTcut->Fill( clust.node, clust.nphesum, weight );
	  if( FPP1track ){
	    hEsum_smear_vs_node_HCALmax_FPP1cut->Fill( clust.node, clust.Esum, weight );
	    hEsum_true_vs_node_HCALmax_FPP1cut->Fill( clust.node, clust.Esum_true, weight );
	    hnphesum_vs_node_HCALmax_FPP1cut->Fill( clust.node, clust.nphesum, weight );
	  }
	}
      }
      
      //Now decide about coincidence trigger. Loop on all combinations of ECAL and HCAL:
      for( auto eclust : ECALclusters ){
	for( auto hclust : HCALclusters ){
	  if( ECAL_nodes_HCAL[hclust.bin].find( eclust.bin ) != ECAL_nodes_HCAL[hclust.bin].end() ){
	    if( hclust.Esum_norm >= nominal_threshold_HCAL && eclust.Esum_norm >= nominal_threshold_ECAL ){
	      cointrig_nominal_threshold = true;
	    }
	    int threshbin_ECAL = hefficiency_vs_threshold_ECAL->FindBin( eclust.Esum_norm );
	    int threshbin_HCAL = hefficiency_vs_threshold_HCAL_FTcut->FindBin( hclust.Esum_norm);

	    if( eclust.Esum_norm < hefficiency_vs_threshold_ECAL->GetBinCenter(threshbin_ECAL) ) threshbin_ECAL--;
	    if( hclust.Esum_norm < hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(threshbin_HCAL) ) threshbin_HCAL--;
	    
	    //we need to set all global bins for which both ECAL and HCAL thresholds are less than the respective trigger sums:
	    for( int binh=1; binh<=threshbin_HCAL; binh++ ){
	      if( binh > hefficiency_vs_threshold_HCAL_FTcut->GetNbinsX() ){
		break;
	      }
	      for( int bine=1; bine<=threshbin_ECAL; bine++ ){
		if( bine > hefficiency_vs_threshold_ECAL->GetNbinsX() ){
		  break;
		}
		int threshbin_global = bine-1 + (binh-1)*hefficiency_vs_threshold_ECAL->GetNbinsX();
		coin_trigger_fired[threshbin_global]++;
	      }
	      
	    }   
	  } //end check for coincidence lookup table

	  //This needs to be outside the check on the coincidence lookup table!
	  // And this is done REGARDLESS of the presence of FT or FPP1 tracks! If we have both ECAL and HCAL above threshold
	  // ANYWHERE for an elastic event, we want that bin combination to go in the lookup table (I should think...)!
	  if( eclust.Esum_norm >= nominal_threshold_ECAL && hclust.Esum_norm >= nominal_threshold_HCAL && FTtrack && FPP1track ){
	    hallnodes_ECAL_vs_HCAL->Fill( hclust.node, eclust.node, weight );
	    hallbins_ECAL_vs_HCAL->Fill( hclust.bin, eclust.bin, weight );
	  }
	  
	} //end loop over HCAL clusters
      } //end loop over ECAL clusters;
      
      if( cointrig_nominal_threshold ){
	if( FTtrack && FPP1track && should_hit_ECAL ) hthetaFPP1_cointrig->Fill(thetaFPP1*180.0/PI,weight);
	//	if( FPP2track ) hthetaFPP2_cointrig->Fill(thetaFPP2*180.0/PI,weight);
      }
    
      for( int bin=1; bin<=hefficiency_vs_threshold_HCAL_FTcut->GetNbinsX(); bin++ ){
	if( FTtrack ) hshouldhit_HCAL_FTcut->Fill( hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(bin), weight );
	if( FTtrack && FPP1track ) hshouldhit_HCAL_FPP1cut->Fill( hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(bin), weight );
	//	if( FTtrack && FPP2track ) hshouldhit_HCAL_FPP2cut->Fill( hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(bin), weight );
	if( trigger_nodes_fired_hcal[bin-1] > 0 ){
	  if( FTtrack ){
	    hefficiency_vs_threshold_HCAL_FTcut->Fill( hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(bin), weight );
	    if( FPP1track ) hefficiency_vs_threshold_HCAL_FPP1cut->Fill( hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(bin), weight );
	    //if( FPP2track ) hefficiency_vs_threshold_HCAL_FPP2cut->Fill( hefficiency_vs_threshold_HCAL_FTcut->GetBinCenter(bin), weight );
	  }
	}
      }

      for( int ithr=0; ithr<coin_trigger_fired.size(); ithr++ ){
	int bin_e = ithr%(hefficiency_vs_threshold_ECAL->GetNbinsX())+1;
	int bin_h = ithr/(hefficiency_vs_threshold_ECAL->GetNbinsX())+1;
	double thr_e = hefficiency_vs_threshold_ECAL_HCAL_coincidence_FTcut->GetYaxis()->GetBinCenter(bin_e);
	double thr_h = hefficiency_vs_threshold_ECAL_HCAL_coincidence_FTcut->GetXaxis()->GetBinCenter(bin_h);
	if( should_hit_ECAL ){
	  if( FTtrack ) hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FTcut->Fill( thr_h, thr_e, weight );
	  if( FTtrack && FPP1track ) hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FPP1cut->Fill( thr_h, thr_e, weight );
	  //if( FTtrack && FPP2track ) hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FPP2cut->Fill( thr_h, thr_e, weight );
	  if( coin_trigger_fired[ithr] > 0 ){
	    if( FTtrack ) hefficiency_vs_threshold_ECAL_HCAL_coincidence_FTcut->Fill( thr_h, thr_e, weight );
	    if( FTtrack && FPP1track ) hefficiency_vs_threshold_ECAL_HCAL_coincidence_FPP1cut->Fill( thr_h, thr_e, weight );
	    // if( FTtrack && FPP2track ) hefficiency_vs_threshold_ECAL_HCAL_coincidence_FPP2cut->Fill( thr_h, thr_e, weight );
	  }
	}
      }

      if( should_hit_ECAL && FTtrack && FPP1track && HCALtrig_nominal ){
	hthetaFPP1_HCALtrig->Fill( thetaFPP1*180.0/PI, weight );
      }
      
      //the "should_hit_ECAL" flag probably redundant here
      if( should_hit_ECAL && FTtrack && FPP1track && maxsum_ECAL >= nominal_threshold_ECAL &&
	  maxsum_HCAL >= nominal_threshold_HCAL && maxnode_ECAL >= 0 && maxnode_HCAL >= 0 ) {
	hmaxnode_ECAL_vs_HCAL->Fill( maxnode_HCAL, maxnode_ECAL, weight );
	
	int binHCAL = HCALclusters[0].bin;
	int binECAL = ECALclusters[0].bin;
	hmaxbin_ECAL_vs_HCAL->Fill( binHCAL, binECAL, weight );
      }
    }
  }

  hefficiency_vs_threshold_HCAL_FTcut->Divide( hshouldhit_HCAL_FTcut );
  hefficiency_vs_threshold_HCAL_FPP1cut->Divide( hshouldhit_HCAL_FPP1cut );
  //hefficiency_vs_threshold_HCAL_FPP2cut->Divide( hshouldhit_HCAL_FPP2cut );

  hefficiency_vs_threshold_ECAL_HCAL_coincidence_FTcut->Divide( hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FTcut );
  hefficiency_vs_threshold_ECAL_HCAL_coincidence_FPP1cut->Divide( hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FPP1cut );
  //hefficiency_vs_threshold_ECAL_HCAL_coincidence_FPP2cut->Divide( hshouldhit_vs_threshold_ECAL_HCAL_coincidence_FPP2cut );
  
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
