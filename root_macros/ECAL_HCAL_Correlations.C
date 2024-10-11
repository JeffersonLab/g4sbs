#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include <iostream>
#include <fstream>
#include <set>
#include <map>



void ECAL_HCAL_Correlations(const char *rootfilename, const char *outfilename, double threshold=0.0){
  TFile *fin = new TFile( rootfilename, "READ" );
  ofstream outfile(outfilename);

  TH2D *hcorrplot;
  fin->GetObject("hallbins_ECAL_vs_HCAL",hcorrplot);

  TH1D *htemp;

  map<int,set<int> > ECAL_list_HCAL;

  TString header;
  header.Form("#%16s, %30s, %15s", "HCAL bin number", "Number of associated ECAL bins", "list of associated ECAL bins");

  outfile << header.Data() << endl;
  
  for( int i=1; i<=hcorrplot->GetNbinsX(); i++ ){
    int nodeH = int(hcorrplot->GetXaxis()->GetBinCenter(i));

    cout << "(i, HCAL bin)=(" << i << ", " << nodeH << ")" << endl;
    htemp = hcorrplot->ProjectionY("",i,i);

    double total = htemp->Integral(1,htemp->GetNbinsX());
    
    for( int j=1; j<=htemp->GetNbinsX(); j++ ){
      if( htemp->GetBinContent(j)/total > threshold ){
	int nodeE = int(htemp->GetBinCenter(j));
	cout << "(j,ECAL bin)=(" << j << ", " << nodeE << ")" << endl;
	ECAL_list_HCAL[nodeH].insert(nodeE);
      }
    }

    TString temp;
    temp.Form( " %16d, %30d", nodeH, ECAL_list_HCAL[nodeH].size() );

    outfile << temp.Data(); 
    //outfile << nodeH << "     " << ECAL_list_HCAL[nodeH].size() << "    ";
    //for( set<int>::iterator k=ECAL_list_HCAL[nodeH].begin(); k!=ECAL_list_HCAL[nodeH].end(); ++k ){
    for( auto k : ECAL_list_HCAL[nodeH] ){
      outfile << ", " << k; 
    }
    outfile << endl;
  }

}
