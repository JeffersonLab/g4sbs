#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include <iostream>
#include <fstream>
#include <set>
#include <map>



void ECAL_HCAL_Correlations(const char *rootfilename, const char *outfilename, double threshold=0.01){
  TFile *fin = new TFile( rootfilename, "READ" );
  ofstream outfile(outfilename);

  TH2D *hcorrplot;
  fin->GetObject("hmaxnode_ECAL_vs_HCAL",hcorrplot);

  TH1D *htemp;

  map<int,set<int> > ECAL_list_HCAL;
  
  for( int i=1; i<=hcorrplot->GetNbinsX(); i++ ){
    int node = i;
    htemp = hcorrplot->ProjectionY("",i,i);

    double total = htemp->Integral(1,htemp->GetNbinsX());
    
    for( int j=1; j<=htemp->GetNbinsX(); j++ ){
      if( htemp->GetBinContent(j)/total >= threshold ){
	ECAL_list_HCAL[node].insert(j);
      }
    }

    outfile << node << " " << ECAL_list_HCAL[node].size();
    for( set<int>::iterator k=ECAL_list_HCAL[node].begin(); k!=ECAL_list_HCAL[node].end(); ++k ){
      outfile << " " << *k; 
    }
    outfile << endl;
  }

}
