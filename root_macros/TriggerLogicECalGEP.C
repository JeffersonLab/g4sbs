#include <iostream>
#include <fstream>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include "TCanvas.h"
#include "TBox.h"
#include "TEllipse.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"

//Based on hard-coded geometry information in G4SBSECal.cc...
void TriggerLogicECalGEP(const char *mapfilename){
  
  // Place the blocks
  int NrowsSM_40 = 10; //total rows = 30 
  int NrowsSM_42 = 15; //total rows = 45
  
  int NcolsSM_40[10] = {3, 5, 6, 7, 8, 8, 9, 9, 9, 9};// from bottom to top
  int NcolsSM_42[15] = {9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 7, 6, 5, 3};// from bottom to top

  int Nrows_40 = 3*NrowsSM_40; //blocks in this section are on a regular rectangular grid
  int Nrows_42 = 3*NrowsSM_42; //blocks in this section are on a (different) rectangular grid. 
  
  double yfp_start_40[10] = {-62.0, -62.0, -62.0, -62.0, -62.0, 
			       -58.0, -54.0, -54.0, -54.0, -54.0};//start of the block edge, cm
  double yfp_start_42[15] = {-54.0, -54.0, -54.0, -54.0, -54.0, 
			       -54.0, -54.0, -54.0, -54.0, -58.0,
			       -62.0, -62.0, -62.0, -62.0, -62.0};//start of the block edge, cm

  int coloffset40[10] = {0,0,0,0,0,
		       1,2,2,2,2};
  int coloffset42[15] = {2,2,2,2,2,
			 2,2,2,2,1,
			 0,0,0,0,0};
  
  double xfpstart = -153.8;
  
  //Basic logic unit is 4 (vertical) x 2 (horizontal) at level 1.
  //At level 2, we combine no more than 4 L1 sums, overlapping by 1 in both vertical and horizontal directions:

  //row, column, X, and Y mapped by cell number:
  set<int> cell_list;
  map<int,int> row_cell;
  map<int,int> col_cell;
  map<int,int> SMrow_cell;
  map<int,int> SMcol_cell;
  map<int,double> X_cell;
  map<int,double> Y_cell;
  map<int,int> L1sum_cell; //L1 sum ID indexed by cell number

  map<int,set<int> > cell_list_L1sum; //list of cell numbers indexed by L1 sum index.
  map<int,double>    XavgL1;
  map<int,double>    YavgL1;
  
  ifstream mapfile(mapfilename);

  TString currentline;

  gStyle->SetOptStat(0);
  
  TCanvas *c1 = new TCanvas("c1","c1",2000,1000);

  c1->Divide(2,1,.001,.001);

  c1->cd(1);

  
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.05);


  
  TH2D *hframe = new TH2D("hframe","",250,-180.0,180.,250,-180.,180.);

  hframe->GetYaxis()->SetTitle("y (cm)");
  hframe->GetXaxis()->SetTitle("x (cm)");

  hframe->Draw();

  int ncolmax40=-1;
  int ncolmax42=-1;
  
  for( int irowSM=0; irowSM<10; irowSM++ ){
    int ncoltemp = 3*NcolsSM_40[irowSM]+coloffset40[irowSM];
    ncolmax40 = (irowSM==0 || ncoltemp > ncolmax40 ) ? ncoltemp : ncolmax40;
  }

  for( int irowSM=0; irowSM<15; irowSM++ ){
    int ncoltemp = 3*NcolsSM_42[irowSM]+coloffset42[irowSM];
    ncolmax42 = (irowSM==0 || ncoltemp > ncolmax42 ) ? ncoltemp : ncolmax42;
  }

  int NcolL1_40 = ncolmax40 / 2;
  if( ncolmax40 % 2 != 0 ) NcolL1_40++;

  int NcolL1_42 = ncolmax42 / 2;
  if( ncolmax42 % 2 != 0 ) NcolL1_42++;

  int NrowL1_40 = Nrows_40/4;
  if( Nrows_40 % 4 != 0 ) NrowL1_40++;

  int NrowL1_42 = Nrows_42/4;
  if( Nrows_42 % 4 != 0 ) NrowL1_42++;
  
  bool first = true;

  TBox b;
  b.SetFillStyle(0);
  b.SetLineStyle(1);
  b.SetLineColor(1);
  b.SetLineWidth(2);
  
  while( currentline.ReadLine(mapfile) ){
    cout << currentline << endl;
    
    if( !currentline.BeginsWith( "#" ) ){ //new cell information:
      TObjArray *tokens = ( (TObjArray*) currentline.Tokenize(",") );
      int ntokens = tokens->GetEntries();
      if( ntokens == 5 ){
	int cell = ( (TObjString*) (*tokens)[0] )->GetString().Atoi();
	int row = ( (TObjString*) (*tokens)[1] )->GetString().Atoi();
	int col = ( (TObjString*) (*tokens)[2] )->GetString().Atoi();
	double xtemp = ( (TObjString*) (*tokens)[3] )->GetString().Atof();
	double ytemp = ( (TObjString*) (*tokens)[4] )->GetString().Atof();

	cell_list.insert(cell);
	row_cell[cell] = row;
	col_cell[cell] = col;
	X_cell[cell] = xtemp;
	Y_cell[cell] = ytemp;

	//c1->cd();
	
	if( row < Nrows_40 ){
	  b.DrawBox( xtemp-2.0, ytemp-2.0, xtemp+2.0,ytemp+2.0 );

	  int coltemp = col + coloffset40[ row/3 ];
	  
	  SMrow_cell[cell]  = row/3;
	  SMcol_cell[cell]  = col/3;

	  int irow_L1 = row/4;
	  int icol_L1 = coltemp/2;
	  int icell_L1 = icol_L1 + irow_L1*NcolL1_40;

	  L1sum_cell[cell] = icell_L1;
	  cell_list_L1sum[icell_L1].insert(cell);
	  
	} else {
	  b.DrawBox( xtemp-2.1, ytemp-2.1, xtemp+2.1,ytemp+2.1);

	  int coltemp = col + coloffset42[ (row-Nrows_40)/3 ];

	  SMrow_cell[cell]  = (row-Nrows_40)/3;
	  SMcol_cell[cell]  = col/3;

	  int irow_L1 = (row-Nrows_40)/4;
	  int icol_L1 = coltemp/2;
	  int icell_L1 = icol_L1 + irow_L1*NcolL1_42 + NrowL1_40*NcolL1_40;

	  L1sum_cell[cell] = icell_L1;
	  cell_list_L1sum[icell_L1].insert(cell);
	  
	}

	//c1->Update();
	
	first = false;

	
      }
    }
  }

  int colors[5] = {2,3,4,5,6};

  //c1->cd(1);
  //hframe->Draw();

  c1->cd(2);
  hframe->Draw();
  
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.05);
  
  for( map<int,set<int> >::iterator i=cell_list_L1sum.begin(); i != cell_list_L1sum.end(); ++i ){
    //Draw transparent boxes of alternating colors around L1 sums.

    int icell_L1 = i->first;
    
    set<int> cell_list_temp = i->second;

    double xsum = 0.;
    double ysum = 0.;
    
    for( set<int>::iterator j = cell_list_temp.begin(); j != cell_list_temp.end(); ++j ){
      int icelltemp = *j;

      int rowtemp = row_cell[icelltemp];
      int coltemp = col_cell[icelltemp];

      b.SetFillStyle(1001);
      b.SetLineStyle(1);
      b.SetLineColor(1);
      b.SetLineWidth(2);
      b.SetFillColorAlpha( colors[icell_L1%5], 0.3 );
      double xtemp = X_cell[icelltemp];
      double ytemp = Y_cell[icelltemp];

      c1->cd(1);
      
      if( rowtemp < Nrows_40 ){
	b.DrawBox( xtemp-2.0, ytemp-2.0, xtemp+2.0,ytemp+2.0 );
      } else {
	b.DrawBox( xtemp-2.1, ytemp-2.1, xtemp+2.1,ytemp+2.1);
      }

      c1->cd(2);
      b.SetFillStyle(0);
      if( rowtemp < Nrows_40 ){
	b.DrawBox( xtemp-2.0, ytemp-2.0, xtemp+2.0,ytemp+2.0 );
      } else {
	b.DrawBox( xtemp-2.1, ytemp-2.1, xtemp+2.1,ytemp+2.1);
      }
      
      xsum += xtemp;
      ysum += ytemp;
    }

    XavgL1[icell_L1] = xsum/double(cell_list_temp.size());
    YavgL1[icell_L1] = ysum/double(cell_list_temp.size());

  }

  c1->cd(2);
  //hframe->Draw();
  
  gPad->SetLeftMargin(0.14);
  gPad->SetRightMargin(0.05);
  gPad->SetBottomMargin(0.14);
  gPad->SetTopMargin(0.05);
  
  //Next, we want to populate the list of cells for Level-2 sums:

  map<int, set<int> > cell_list_L2sum; //list of cells by L2 sum
  map<int, set<int> > L2sum_cell; //list of L2 sum(s) by cell

  int icell_L2 = 0;
  //First: make all groups containing only blocks in the 4.0 cm section:
  for( int irow_L1=0; irow_L1<NrowL1_40-1; irow_L1++ ){
    for( int icol_L1=0; icol_L1<NcolL1_40-1; icol_L1++ ){

      int i1 = icol_L1 + irow_L1*NcolL1_40;
      int i2 = (icol_L1+1) + irow_L1*NcolL1_40; //col+1
      int i3 = icol_L1 + (irow_L1+1)*NcolL1_40; //row+1
      int i4 = (icol_L1+1) + (irow_L1+1)*NcolL1_40; //col+1, row+1;

      double xmin=1000.0,xmax=-1000.0,ymin=1000.0,ymax=-1000.0;
      
      if( cell_list_L1sum.find( i1 ) != cell_list_L1sum.end() ){
	for( set<int>::iterator k = cell_list_L1sum[i1].begin(); k != cell_list_L1sum[i1].end(); ++k ){
	  int cell = *k;
	  cell_list_L2sum[icell_L2].insert( cell );
	  L2sum_cell[cell].insert(icell_L2);

	  double xtemp = X_cell[cell];
	  double ytemp = Y_cell[cell];

	  xmin = (xtemp - 2.0 < xmin ) ? xtemp - 2.0 : xmin;
	  xmax = (xtemp + 2.0 > xmax ) ? xtemp + 2.0 : xmax;
	  ymin = (ytemp - 2.0 < ymin ) ? ytemp - 2.0 : ymin;
	  ymax = (ytemp + 2.0 > ymax ) ? ytemp + 2.0 : ymax;

	  b.SetFillStyle(1001);
	  b.SetFillColorAlpha( colors[icell_L2%5], 0.2 );
	  b.SetLineStyle(1);
	  b.SetLineColor( colors[icell_L2%5] );
	  b.SetLineWidth( 1 );
	  b.DrawBox( xtemp-1.9, ytemp-1.9, xtemp+1.9, ytemp+1.9 );
	}
      }
       
      if( cell_list_L1sum.find( i2 ) != cell_list_L1sum.end() ){
	for( set<int>::iterator k = cell_list_L1sum[i2].begin(); k != cell_list_L1sum[i2].end(); ++k ){
	  int cell = *k;
	  cell_list_L2sum[icell_L2].insert( cell );
	  L2sum_cell[cell].insert(icell_L2);

	  double xtemp = X_cell[cell];
	  double ytemp = Y_cell[cell];

	  xmin = (xtemp - 2.0 < xmin ) ? xtemp - 2.0 : xmin;
	  xmax = (xtemp + 2.0 > xmax ) ? xtemp + 2.0 : xmax;
	  ymin = (ytemp - 2.0 < ymin ) ? ytemp - 2.0 : ymin;
	  ymax = (ytemp + 2.0 > ymax ) ? ytemp + 2.0 : ymax;

	  b.SetFillStyle(1001);
	  b.SetFillColorAlpha( colors[icell_L2%5], 0.2 );
	  b.SetLineStyle(1);
	  b.SetLineColor( colors[icell_L2%5] );
	  b.SetLineWidth( 1 );
	  b.DrawBox( xtemp-1.9, ytemp-1.9, xtemp+1.9, ytemp+1.9 );
	}
      }

      if( cell_list_L1sum.find( i3 ) != cell_list_L1sum.end() ){
	for( set<int>::iterator k = cell_list_L1sum[i3].begin(); k != cell_list_L1sum[i3].end(); ++k ){
	  int cell = *k;
	  cell_list_L2sum[icell_L2].insert( cell );
	  L2sum_cell[cell].insert(icell_L2);

	  double xtemp = X_cell[cell];
	  double ytemp = Y_cell[cell];

	  xmin = (xtemp - 2.0 < xmin ) ? xtemp - 2.0 : xmin;
	  xmax = (xtemp + 2.0 > xmax ) ? xtemp + 2.0 : xmax;
	  ymin = (ytemp - 2.0 < ymin ) ? ytemp - 2.0 : ymin;
	  ymax = (ytemp + 2.0 > ymax ) ? ytemp + 2.0 : ymax;

	  b.SetFillStyle(1001);
	  b.SetFillColorAlpha( colors[icell_L2%5], 0.2 );
	  b.SetLineStyle(1);
	  b.SetLineColor( colors[icell_L2%5] );
	  b.SetLineWidth( 1 );
	  b.DrawBox( xtemp-1.9, ytemp-1.9, xtemp+1.9, ytemp+1.9 );
	}
      }

      if( cell_list_L1sum.find( i4 ) != cell_list_L1sum.end() ){
	for( set<int>::iterator k = cell_list_L1sum[i4].begin(); k != cell_list_L1sum[i4].end(); ++k ){
	  int cell = *k;
	  cell_list_L2sum[icell_L2].insert( cell );
	  L2sum_cell[cell].insert(icell_L2);

	  double xtemp = X_cell[cell];
	  double ytemp = Y_cell[cell];

	  xmin = (xtemp - 2.0 < xmin ) ? xtemp - 2.0 : xmin;
	  xmax = (xtemp + 2.0 > xmax ) ? xtemp + 2.0 : xmax;
	  ymin = (ytemp - 2.0 < ymin ) ? ytemp - 2.0 : ymin;
	  ymax = (ytemp + 2.0 > ymax ) ? ytemp + 2.0 : ymax;

	  b.SetFillStyle(1001);
	  b.SetFillColorAlpha( colors[icell_L2%5], 0.2 );
	  b.SetLineStyle(1);
	  b.SetLineColor( colors[icell_L2%5] );
	  b.SetLineWidth( 1 );
	  b.DrawBox( xtemp-1.9, ytemp-1.9, xtemp+1.9, ytemp+1.9 );
	}
      }
      
      if( cell_list_L2sum[icell_L2].size() > 0 ){	
	icell_L2++;
	c1->Update();
	cout << "Press Enter to continue";
      
	currentline.ReadLine(cin,kFALSE);
      }
	//}

      
    }
  }

  //We also need to populate the edge region between the two sections:

  int irow_top_40 = NrowL1_40-1;
  for( int icol_L1=0; icol_L1<NcolL1_40-1; icol_L1++ ){

    int icell_L1_40 = icol_L1 + irow_top_40*NcolL1_40;
    if( cell_list_L1sum.find(icell_L1_40) != cell_list_L1sum.end() ){

      double minXsep = 1000.0;
      int icol_L1_42_nearest = -1;
      for( int jcol=0; jcol<NcolL1_42-1; jcol++ ){ //find the nearest-neighbor L1 sums in X
	int icell_L1_42 = jcol + NcolL1_40*NrowL1_40;
	if( cell_list_L1sum.find(icell_L1_42) != cell_list_L1sum.end() ){
	  double Xsep = fabs( XavgL1[icell_L1_40] - XavgL1[icell_L1_42] );
	  icol_L1_42_nearest = (Xsep < minXsep || icol_L1_42_nearest==-1 ) ? jcol : icol_L1_42_nearest;
	  minXsep = (Xsep < minXsep ) ? Xsep : minXsep;
	}
      }

      cout << "Nearest neighbor column = " << icol_L1_42_nearest << endl;
      
      if( icol_L1_42_nearest >= 0 ){
	int i1 = icol_L1 + irow_top_40*NcolL1_40;
	int i2 = (icol_L1+1) + irow_top_40*NcolL1_40;
	int i3 = icol_L1_42_nearest + NcolL1_40*NrowL1_40;
	int i4 = (icol_L1_42_nearest+1) + NcolL1_40*NrowL1_40;

	if( cell_list_L1sum.find( i1 ) != cell_list_L1sum.end() ){
	  for( set<int>::iterator k = cell_list_L1sum[i1].begin(); k != cell_list_L1sum[i1].end(); ++k ){
	    int cell = *k;
	    cell_list_L2sum[icell_L2].insert( cell );
	    L2sum_cell[cell].insert(icell_L2);
	    
	    double xtemp = X_cell[cell];
	    double ytemp = Y_cell[cell];
	    
	   
	    
	    b.SetFillStyle(1001);
	    b.SetFillColorAlpha( colors[icell_L2%5], 0.2 );
	    b.SetLineStyle(1);
	    b.SetLineColor( colors[icell_L2%5] );
	    b.SetLineWidth( 1 );
	    b.DrawBox( xtemp-1.9, ytemp-1.9, xtemp+1.9, ytemp+1.9 );
	  }
	}

	if( cell_list_L1sum.find( i2 ) != cell_list_L1sum.end() ){
	  for( set<int>::iterator k = cell_list_L1sum[i2].begin(); k != cell_list_L1sum[i2].end(); ++k ){
	    int cell = *k;
	    cell_list_L2sum[icell_L2].insert( cell );
	    L2sum_cell[cell].insert(icell_L2);
	    
	    double xtemp = X_cell[cell];
	    double ytemp = Y_cell[cell];
	    
	    
	    
	    b.SetFillStyle(1001);
	    b.SetFillColorAlpha( colors[icell_L2%5], 0.2 );
	    b.SetLineStyle(1);
	    b.SetLineColor( colors[icell_L2%5] );
	    b.SetLineWidth( 1 );
	    b.DrawBox( xtemp-1.9, ytemp-1.9, xtemp+1.9, ytemp+1.9 );
	  }
	}

	if( cell_list_L1sum.find( i3 ) != cell_list_L1sum.end() ){
	  for( set<int>::iterator k = cell_list_L1sum[i3].begin(); k != cell_list_L1sum[i3].end(); ++k ){
	    int cell = *k;
	    cell_list_L2sum[icell_L2].insert( cell );
	    L2sum_cell[cell].insert(icell_L2);
	    
	    double xtemp = X_cell[cell];
	    double ytemp = Y_cell[cell];
	    
	   
	    
	    b.SetFillStyle(1001);
	    b.SetFillColorAlpha( colors[icell_L2%5], 0.2 );
	    b.SetLineStyle(1);
	    b.SetLineColor( colors[icell_L2%5] );
	    b.SetLineWidth( 1 );
	    b.DrawBox( xtemp-1.9, ytemp-1.9, xtemp+1.9, ytemp+1.9 );
	  }
	}

	if( cell_list_L1sum.find( i4 ) != cell_list_L1sum.end() ){
	  for( set<int>::iterator k = cell_list_L1sum[i4].begin(); k != cell_list_L1sum[i4].end(); ++k ){
	    int cell = *k;
	    cell_list_L2sum[icell_L2].insert( cell );
	    L2sum_cell[cell].insert(icell_L2);
	    
	    double xtemp = X_cell[cell];
	    double ytemp = Y_cell[cell];
	    
	   
	    
	    b.SetFillStyle(1001);
	    b.SetFillColorAlpha( colors[icell_L2%5], 0.2 );
	    b.SetLineStyle(1);
	    b.SetLineColor( colors[icell_L2%5] );
	    b.SetLineWidth( 1 );
	    b.DrawBox( xtemp-1.9, ytemp-1.9, xtemp+1.9, ytemp+1.9 );
	  }
	}	
      }

      if( cell_list_L2sum[icell_L2].size() > 0 ){
	
	// b.SetFillStyle(1001);
	// b.SetFillColorAlpha( colors[icell_L2%5], 0.2 );
	// b.SetLineStyle(3);
	// b.SetLineColor( colors[icell_L2%5] );
	// b.SetLineWidth( 3 );
	// b.DrawBox( xmin-1.0, ymin-1.0, xmax+1.0, ymax+1.0 );
	
	
	
	icell_L2++;
      
	c1->Update();
	cout << "Press Enter to continue";
	
	currentline.ReadLine(cin,kFALSE);
      }
    }
  }
  
  
  for( int irow_L1=0; irow_L1<NrowL1_42-1; irow_L1++ ){
    for( int icol_L1=0; icol_L1<NcolL1_42-1; icol_L1++ ){

      int i1 = NrowL1_40*NcolL1_40 + icol_L1 + irow_L1*NcolL1_42;
      int i2 = NrowL1_40*NcolL1_40 + (icol_L1+1) + irow_L1*NcolL1_42; //col+1
      int i3 = NrowL1_40*NcolL1_40 + icol_L1 + (irow_L1+1)*NcolL1_42; //row+1
      int i4 = NrowL1_40*NcolL1_40 + (icol_L1+1) + (irow_L1+1)*NcolL1_42; //col+1, row+1;

      //Make sure that each L1 sum is only used as a seed for a logic group once:

      double xmin=1000.0,xmax=-1000.0,ymin=1000.0,ymax=-1000.0;
      
      if( cell_list_L1sum.find( i1 ) != cell_list_L1sum.end() ){

	for( set<int>::iterator k = cell_list_L1sum[i1].begin(); k != cell_list_L1sum[i1].end(); ++k ){
	  int cell = *k;
	  cell_list_L2sum[icell_L2].insert( cell );
	  L2sum_cell[cell].insert(icell_L2);

	  double xtemp = X_cell[cell];
	  double ytemp = Y_cell[cell];

	  xmin = (xtemp - 2.0 < xmin ) ? xtemp - 2.0 : xmin;
	  xmax = (xtemp + 2.0 > xmax ) ? xtemp + 2.0 : xmax;
	  ymin = (ytemp - 2.0 < ymin ) ? ytemp - 2.0 : ymin;
	  ymax = (ytemp + 2.0 > ymax ) ? ytemp + 2.0 : ymax;

	  b.SetFillStyle(1001);
	  b.SetFillColorAlpha( colors[icell_L2%5], 0.2 );
	  b.SetLineStyle(1);
	  b.SetLineColor( colors[icell_L2%5] );
	  b.SetLineWidth( 1 );
	  b.DrawBox( xtemp-1.9, ytemp-1.9, xtemp+1.9, ytemp+1.9 );
	}
      }
      if( cell_list_L1sum.find( i2 ) != cell_list_L1sum.end() ){
	for( set<int>::iterator k = cell_list_L1sum[i2].begin(); k != cell_list_L1sum[i2].end(); ++k ){
	  int cell = *k;
	  cell_list_L2sum[icell_L2].insert( cell );
	  L2sum_cell[cell].insert(icell_L2);

	  double xtemp = X_cell[cell];
	  double ytemp = Y_cell[cell];

	  xmin = (xtemp - 2.0 < xmin ) ? xtemp - 2.0 : xmin;
	  xmax = (xtemp + 2.0 > xmax ) ? xtemp + 2.0 : xmax;
	  ymin = (ytemp - 2.0 < ymin ) ? ytemp - 2.0 : ymin;
	  ymax = (ytemp + 2.0 > ymax ) ? ytemp + 2.0 : ymax;

	  b.SetFillStyle(1001);
	  b.SetFillColorAlpha( colors[icell_L2%5], 0.2 );
	  b.SetLineStyle(1);
	  b.SetLineColor( colors[icell_L2%5] );
	  b.SetLineWidth( 1 );
	  b.DrawBox( xtemp-1.9, ytemp-1.9, xtemp+1.9, ytemp+1.9 );
	}
      }

      if( cell_list_L1sum.find( i3 ) != cell_list_L1sum.end() ){
	for( set<int>::iterator k = cell_list_L1sum[i3].begin(); k != cell_list_L1sum[i3].end(); ++k ){
	  int cell = *k;
	  cell_list_L2sum[icell_L2].insert( cell );
	  L2sum_cell[cell].insert(icell_L2);

	  double xtemp = X_cell[cell];
	  double ytemp = Y_cell[cell];

	  xmin = (xtemp - 2.0 < xmin ) ? xtemp - 2.0 : xmin;
	  xmax = (xtemp + 2.0 > xmax ) ? xtemp + 2.0 : xmax;
	  ymin = (ytemp - 2.0 < ymin ) ? ytemp - 2.0 : ymin;
	  ymax = (ytemp + 2.0 > ymax ) ? ytemp + 2.0 : ymax;

	  b.SetFillStyle(1001);
	  b.SetFillColorAlpha( colors[icell_L2%5], 0.2 );
	  b.SetLineStyle(1);
	  b.SetLineColor( colors[icell_L2%5] );
	  b.SetLineWidth( 1 );
	  b.DrawBox( xtemp-1.9, ytemp-1.9, xtemp+1.9, ytemp+1.9 );
	}
      }

      if( cell_list_L1sum.find( i4 ) != cell_list_L1sum.end() ){
	for( set<int>::iterator k = cell_list_L1sum[i4].begin(); k != cell_list_L1sum[i4].end(); ++k ){
	  int cell = *k;
	  cell_list_L2sum[icell_L2].insert( cell );
	  L2sum_cell[cell].insert(icell_L2);

	  double xtemp = X_cell[cell];
	  double ytemp = Y_cell[cell];

	  xmin = (xtemp - 2.0 < xmin ) ? xtemp - 2.0 : xmin;
	  xmax = (xtemp + 2.0 > xmax ) ? xtemp + 2.0 : xmax;
	  ymin = (ytemp - 2.0 < ymin ) ? ytemp - 2.0 : ymin;
	  ymax = (ytemp + 2.0 > ymax ) ? ytemp + 2.0 : ymax;

	  b.SetFillStyle(1001);
	  b.SetFillColorAlpha( colors[icell_L2%5], 0.2 );
	  b.SetLineStyle(1);
	  b.SetLineColor( colors[icell_L2%5] );
	  b.SetLineWidth( 1 );
	  b.DrawBox( xtemp-1.9, ytemp-1.9, xtemp+1.9, ytemp+1.9 );
	}
      }
 
      if( cell_list_L2sum[icell_L2].size() > 0 ){
	
	// b.SetFillStyle(1001);
	// b.SetFillColorAlpha( colors[icell_L2%5], 0.2 );
	// b.SetLineStyle(3);
	// b.SetLineColor( colors[icell_L2%5] );
	// b.SetLineWidth( 3 );
	// b.DrawBox( xmin-1.0, ymin-1.0, xmax+1.0, ymax+1.0 );
	  
    
	  
	icell_L2++;

	c1->Update();
	cout << "Press Enter to continue";
	
	currentline.ReadLine(cin,kFALSE);
      }

      
    }
  }

  ofstream L1file("GEP_ECAL_L1sums.txt"); //listing of cell numbers by level 1 sum index
  ofstream L2file("GEP_ECAL_L2sums.txt"); //listing of cell numbers by level 2 sum index
  //ofstream L1file_by_cell("GEP_ECAL_L1sums_by_cell.txt"); //listing of L1 sum (should only be 1) per cell
  //ofstream L2file_by_cell("GEP_ECAL_L2sums_by_cell.txt"); //listing of L2 sums (can be up to four? six? per cell)

  int ngood1 = 0;
  
  for( map<int,set<int> >::iterator i = cell_list_L1sum.begin(); i != cell_list_L1sum.end(); ++i ){
    //What info to include: L1 sum index, ncells;
    //Then: cell information:
    int icell_L1 = i->first;

    set<int> cell_list_temp = i->second;

    TString stemp;
    stemp.Form( " %8d  %8d  ", ngood1, cell_list_temp.size() );

    currentline = "";
    currentline += stemp;
    
    for( set<int>::iterator j = cell_list_temp.begin(); j != cell_list_temp.end(); ++j ){
      int icell = *j;

      stemp.Form( "%7d ", icell );
      currentline += stemp;
    }
    L1file << currentline << endl;

    ngood1++;
  }

  int ngood2 = 0; //Number of "good" level-2 sums
  int minblock = 4;
  
  for( map<int,set<int> >::iterator i = cell_list_L2sum.begin(); i != cell_list_L2sum.end(); ++i ){
    int icell_L2 = i->first;
    set<int> cell_list_temp = i->second;

    TString stemp;
    stemp.Form( " %16d  %16d  ", ngood2, cell_list_temp.size() );
    currentline = "";
    currentline += stemp;
    for( set<int>::iterator j = cell_list_temp.begin(); j != cell_list_temp.end(); ++j ){
      int icell = *j;
      stemp.Form( "%7d ", icell );
      currentline += stemp;
    }

    if( cell_list_temp.size() > minblock ){
      L2file << currentline << endl;
      ngood2++;
    }
  }
  
    
  
  
  
}
