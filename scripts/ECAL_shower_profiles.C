#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include <vector>
#include "TFile.h"

void get_shower_profile(TH1D *h, int &nbins, double &min, double &max, vector<double> &frac ){
  //Make acceptance-averaged profiles:

  //double xmom_min=-1.0,xmom_max=1.0
  //vector<double> frac;

  int binmax = h->GetMaximumBin();
  int binhigh=binmax;
  int binlow=binmax;

  nbins = h->GetNbinsX();
  double tot = h->Integral();

  while( h->Integral(1,binlow)/tot >= 1.e-4 && binlow > 1 ){binlow--;}
  while( h->Integral(binhigh,nbins)/tot >= 1.e-4 && binhigh < nbins ){binhigh++;}

  tot = h->Integral(binlow,binhigh);
  min = h->GetBinLowEdge(binlow);
  max = h->GetBinLowEdge(binhigh)+h->GetBinWidth(binhigh);

  //vector<double> frac;

  frac.clear();
  
  frac.push_back(0.0);
  for(int bin=binlow; bin<=binhigh; bin++ ){
    frac.push_back( h->Integral(binlow,bin)/tot );
    cout << "bin, frac = " << bin << ", " << frac[frac.size()-1] << endl;
  }

  nbins = frac.size()-1;
}

void ECAL_shower_profiles( const char *infilename, const char *outfilename="ECAL_shower_profiles.txt" ){

  ofstream outfile(outfilename);
  
  TFile *fin = new TFile(infilename,"READ");


  TH1D *hxmom_avg,*hymom_avg;
  TH2D *hxmom_x, *hymom_y;

  fin->GetObject("hxmom",hxmom_avg);
  fin->GetObject("hymom",hymom_avg);
  fin->GetObject("hxmom_x",hxmom_x);
  fin->GetObject("hymom_y",hymom_y);

  //How should output look? we need a "DEFAULT" x and y distribution, and then x and y dependent distributions:
  
  double xmin_avg, xmax_avg;
  int nbins_xavg;
  vector<double> fracx_avg;
  
  get_shower_profile( hxmom_avg, nbins_xavg, xmin_avg, xmax_avg, fracx_avg );
  
  //outfile << "ECAL_shower_profile_x_default" << endl;

  //outfile << "ECAL_shower_profile_x_default" << endl;

  
  TString currentline;
  currentline.Form("ECAL_shower_profile_x_default %d  %9.5g %9.5g", nbins_xavg, xmin_avg, xmax_avg );

  outfile << currentline << endl;

  for( int bin=0; bin<fracx_avg.size(); bin++ ){
    currentline.Form(" %15.6g ", fracx_avg[bin] );
    outfile << currentline;
    if( (bin+1)%8 == 0 ) outfile << endl;
  }

  outfile << endl << endl;;

  int nbins_yavg;
  double ymin_avg, ymax_avg;
  vector<double> fracy_avg;

  get_shower_profile( hymom_avg, nbins_yavg, ymin_avg, ymax_avg, fracy_avg );

  //outfile << "" << endl;
  
  currentline.Form( "ECAL_shower_profile_y_default %d  %9.5g %9.5g", nbins_yavg, ymin_avg, ymax_avg );

  outfile << currentline << endl;
  
  for(int bin=0; bin<fracy_avg.size(); bin++ ){
    currentline.Form(" %15.6g ", fracy_avg[bin] );
    outfile << currentline;
    if( (bin+1)%8 == 0 ) outfile << endl;
  }
  outfile << endl << endl;
  
  //  outfile << endl << endl;

  int nbins_xclust = hxmom_x->GetNbinsX();

  double xclust_min = hxmom_x->GetXaxis()->GetXmin();
  double xclust_max = hxmom_x->GetXaxis()->GetXmax();

  currentline.Form("ECAL_shower_profiles_x %d  %9.5g %9.5g", nbins_xclust, xclust_min, xclust_max );

  outfile << currentline << endl;

  TH1D *htemp;
  
  for( int bin=1; bin<=nbins_xclust; bin++ ){
    int nbins_xtemp;
    double xmintemp,xmaxtemp;
    vector<double> fracxtemp;

    TString histname;
    histname.Form( "hxmom_x_py_bin%d", bin );

    htemp = hxmom_x->ProjectionY(histname.Data(), bin,bin );
    
    get_shower_profile( htemp, nbins_xtemp, xmintemp, xmaxtemp, fracxtemp );

    currentline.Form( "ECAL_shower_profile_x_bin %d  %d  %9.5g %9.5g", bin, nbins_xtemp, xmintemp, xmaxtemp );

    outfile << currentline << endl;

    for( int j=0; j<fracxtemp.size(); j++ ){
      currentline.Form(" %15.6g ", fracxtemp[j] );
      outfile << currentline;
      if( (j+1)%8 == 0 ) outfile << endl;
    }
    outfile << endl << endl;
  }

  outfile << endl << endl;


  int nbins_yclust = hymom_y->GetNbinsX();

  double yclust_min = hymom_y->GetXaxis()->GetXmin();
  double yclust_max = hymom_y->GetXaxis()->GetXmax();

  currentline.Form("ECAL_shower_profiles_y %d  %9.5g %9.5g", nbins_yclust, yclust_min, yclust_max );

  outfile << currentline << endl;

  //TH1D *htemp;
  
  for( int bin=1; bin<=nbins_yclust; bin++ ){
    int nbins_ytemp;
    double ymintemp,ymaxtemp;
    vector<double> fracytemp;

    TString histname;
    histname.Form( "hymom_y_py_bin%d", bin );

    htemp = hymom_y->ProjectionY(histname.Data(), bin,bin );
    
    get_shower_profile( htemp, nbins_ytemp, ymintemp, ymaxtemp, fracytemp );

    currentline.Form( "ECAL_shower_profile_y_bin %d  %d  %9.5g %9.5g", bin, nbins_ytemp, ymintemp, ymaxtemp );

    outfile << currentline << endl;

    for( int j=0; j<fracytemp.size(); j++ ){
      currentline.Form(" %15.6g ", fracytemp[j] );
      outfile << currentline;
      if( (j+1)%8 == 0 ) outfile << endl;
    }
    outfile << endl << endl;
  }

  outfile << endl << endl;

  

  
}
