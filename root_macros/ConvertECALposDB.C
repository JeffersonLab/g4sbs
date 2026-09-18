#include <iostream>
#include <fstream>
#include "TString.h"
#include <vector>
#include <set>
#include <map>
#include "TObjArray.h"
#include "TObjString.h"

using namespace std;

const int nrows=69;

void ConvertECALposDB( const char *inputfilename="database/ecal_gep_blockmap.txt", const char *outfilename="ECALposDB_temp.dat", double xoff=0.0, double yoff=-0.05946 ){

  int ncols_row[69] = { 12, 12, 12,
			18, 18, 18,
			21, 21, 21, 21, 21, 21,
			24, 24, 24, 24, 24, 24,
			27, 27, 27, 27, 27, 27,
			27, 27, 27, 27, 27, 27,
			27, 27, 27, 27, 27, 27,
			27, 27, 27, 27, 27, 27,
			27, 27, 27, 27, 27, 27,
			27, 27, 27, 27, 27, 27,
			24, 24, 24, 24, 24, 24,
			21, 21, 21, 21, 21, 21,
			18, 18, 18 };

  ifstream infile(inputfilename);
  
  if( !infile ) return;

  ofstream outfile(outfilename);

  std::set<int> list_cells;
  std::map<int,int> row_cell;
  std::map<int,int> col_cell;
  std::map<int,double> ypos_cell;
  std::map<int,double> xpos_cell;
  
  TString currentline;
  
  while( currentline.ReadLine(infile) ){

    if( !currentline.BeginsWith( "#" ) ){
      TObjArray *tokens = currentline.Tokenize(",");

      if( tokens->GetEntries() >= 5 ){
	TString scell = ( (TObjString*) (*tokens)[0] )->GetString();
	TString srow = ( (TObjString*) (*tokens)[1] )->GetString();
	TString scol = ( (TObjString*) (*tokens)[2] )->GetString();
	TString sycell = ( (TObjString*) (*tokens)[3] )->GetString();
	TString sxcell = ( (TObjString*) (*tokens)[4] )->GetString();

	int cell = scell.Atoi();
	int row = srow.Atoi();
	int col = scol.Atoi();
	double xpos = -sxcell.Atof()/100.0;
	double ypos = sycell.Atof()/100.0;
	
	list_cells.insert(cell);
	row_cell[cell] = row;
	col_cell[cell] = col;
	ypos_cell[cell] = ypos + yoff;
	xpos_cell[cell] = xpos + xoff;
	
      }

    }
  }

  outfile << "earm.ecal.xpos = " << endl;

  int lastrow = -1;
  
  for( auto cell : list_cells ){
    if( cell > 0 && col_cell[cell] == 0 ){ //new row:
      outfile << " # ECAL row " << row_cell[cell]-1 << endl;
    }
    outfile << Form("%10.7g ", xpos_cell[cell]);
    lastrow = row_cell[cell];
  }

  outfile << " # ECAL row " << lastrow << endl << endl;

  lastrow = -1;

  outfile << "earm.ecal.ypos = " << endl;
  for( auto cell : list_cells ){
    if( cell > 0 && col_cell[cell] == 0 ){ //new row:
      outfile << " # ECAL row " << row_cell[cell]-1 << endl;
    }
    outfile << Form("%10.7g ", ypos_cell[cell]);
    lastrow = row_cell[cell];
  }

  outfile << " # ECAL row " << lastrow << endl;
  
   
  
}
