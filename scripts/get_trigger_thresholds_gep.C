#include "TH2D.h"
#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include <iostream>
#include <fstream>
#include "TSystem.h"

void get_trigger_thresholds_gep( const char *rootfilename, const char *outfilename_ECAL="ecal_trigger_temp.txt", const char *outfilename_HCAL="hcal_trigger_temp.txt" ){

  TH2D *hnodesums_ECAL;
  TH2D *hnodesums_HCAL;

  TFile *fin = new TFile(rootfilename,"READ");
  fin->GetObject( "hEsum_smear_vs_node_ECAL_max", hnodesums_ECAL );
  //fin->GetObject( "hnphesum_vs_node_ECAL_all", hnodesums_ECAL );
  fin->GetObject( "hEsum_smear_vs_node_HCALmax_FPP1cut", hnodesums_HCAL );

  TF1 *gaussfunc = new TF1("gaussfunc","[0]*exp(-0.5*pow((x-[1])/[2],2))", 0.0,5000.0);
  TF1 *gaussplusexpo = new TF1("gaussplusexpo","[0]*exp(-[1]*x)+[2]*exp(-0.5*pow((x-[3])/[4],2))",0.0,5000.0 );

  double startpar_gauss[3] = {1.0,3.0,0.2};
  double startpar_gaussexpo[5] = {.1,.001,1,1500.0,300.0};

  gaussfunc->SetParameters(startpar_gauss);
  gaussfunc->SetParLimits(2,0,5000);
  gaussplusexpo->SetParameters(startpar_gaussexpo);

  //do ECAL first:
  int nnodes_ECAL = hnodesums_ECAL->GetNbinsX();

  TH1D *htemp;

  TCanvas *ctemp;
  
  ctemp = new TCanvas("c0","c0",1200,900);
  ctemp->Divide(5,5,.001,.001);

  double nodes_ECAL[nnodes_ECAL], ex[nnodes_ECAL], mean_ECAL[nnodes_ECAL], sigma_ECAL[nnodes_ECAL], dmean_ECAL[nnodes_ECAL], dsigma_ECAL[nnodes_ECAL];
  
  for( int i=1; i<=nnodes_ECAL; i++ ){

    //double startpar_gauss[3] = {1.0,3000.0,200.0};
    
    int page = (i-1)/25;
    int panel = (i-1)%25 + 1;
    
    int node = int( hnodesums_ECAL->GetXaxis()->GetBinCenter(i) );
    htemp = hnodesums_ECAL->ProjectionY("",i,i);

    gaussfunc->SetParameters(startpar_gauss);
    //gaussfunc->SetParameter(
    
    gaussfunc->SetParameter(0,0.4*htemp->GetBinContent(htemp->GetMaximumBin()));
    
    TString title;
    title.Form("ECAL node %d",node);

    htemp->SetTitle(title.Data());
    double xmin = 0.5, xmax = 7.0;
    int binmin = htemp->FindBin(xmin);
    int binmax = htemp->FindBin(xmax);

    htemp->GetXaxis()->SetRange(binmin,binmax);

    gaussfunc->SetParameter(1,htemp->GetBinCenter(htemp->GetMaximumBin()));
    gaussfunc->SetParameter(2,0.2);
    
    cout << "panel = " << panel << endl;
    
    ctemp->cd(panel);
    ctemp->Update();

    int maxbin = htemp->GetMaximumBin();
    double max = htemp->GetBinContent(maxbin);
    int binlo = maxbin, binhi = maxbin;

    while( htemp->GetBinContent(binlo--) >= 0.55*max && binlo >= binmin ){};
    while( htemp->GetBinContent(binhi++) >= 0.01*max && binhi <= binmax ){};

    double xlo = htemp->GetBinCenter( binlo );
    double xhi = htemp->GetBinCenter( binhi );

    bool okay = false;

    while( !okay ){

      //gSystem->Sleep(10);
      okay = true;
      
      htemp->Fit(gaussfunc,"SQ","",xlo,xhi );
      htemp->DrawCopy();
      ctemp->Update();

      gPad->Modified();

      gSystem->ProcessEvents();
      
      // cout << "fit okay (y/n)?" << endl;
      // TString reply;
      // reply.ReadLine(cin,kFALSE);

      // if( !reply.BeginsWith("n") ){
      // 	okay = true;
      // } else {
      // 	cout << "xmin = ";
      // 	cin >> xlo;
      // 	cout << "xmax = ";
      // 	cin >> xhi;
      // }
    }

    nodes_ECAL[i-1] = node;
    mean_ECAL[i-1] = gaussfunc->GetParameter(1);
    dmean_ECAL[i-1] = gaussfunc->GetParError(1);
    sigma_ECAL[i-1] = gaussfunc->GetParameter(2);
    dsigma_ECAL[i-1] = gaussfunc->GetParError(2);
    ex[i-1] = 0;
    
    if( panel == 25 || i == nnodes_ECAL ){

      ctemp->Update();
      
      if( page == 0 ){
	ctemp->Print( "ECAL_sums.pdf(","pdf");
      } else {
	ctemp->Print( "ECAL_sums.pdf","pdf");
      }
      
      ctemp->Clear();
      
      ctemp->Divide(5,5,.001,.001);
      //ctemp->Update();
    }
  }

  //  ctemp->Print("ECAL_sums.pdf)","pdf");
  
  ofstream ecalfile( outfilename_ECAL );

  for( int inode=0; inode<nnodes_ECAL; inode++ ){
    TString line;
    line.Form( "%5d %15.6g %15.6g", int(nodes_ECAL[inode]), mean_ECAL[inode], sigma_ECAL[inode] );
    ecalfile << line.Data() << endl;
  }
  
  
  
  TGraphErrors *gmean_ECAL = new TGraphErrors(nnodes_ECAL, nodes_ECAL, mean_ECAL, ex, dmean_ECAL );
  TGraphErrors *gsigma_ECAL = new TGraphErrors(nnodes_ECAL, nodes_ECAL, sigma_ECAL, ex, dsigma_ECAL );
  gmean_ECAL->SetMarkerStyle(20);
  gsigma_ECAL->SetMarkerStyle(20);

  ctemp = new TCanvas( "cgraphs","cgraphs",1200,900);
  ctemp->Divide(2,1,.001,.001);
  ctemp->cd(1);
  gmean_ECAL->Draw("AP");
  ctemp->cd(2);
  gsigma_ECAL->Draw("AP");
  ctemp->Print("ECAL_sums.pdf)","pdf");

  TCanvas *ch = new TCanvas("ch","ch",1200,900);
  ch->Divide(5,5,.001,.001);
  
  int nnodes_HCAL = hnodesums_HCAL->GetNbinsX();
  
  //  TH1D *hnodesums_HCAL_py = hnodesums_HCAL->ProjectionY("hnodesums_HCAL_py",1,hnodesums_HCAL->GetNbinsX() );

  double nodes_HCAL[nnodes_HCAL], exh[nnodes_HCAL], mean_HCAL[nnodes_HCAL], sigma_HCAL[nnodes_HCAL], dmean_HCAL[nnodes_HCAL], dsigma_HCAL[nnodes_HCAL];

  startpar_gauss[0] = 1.;
  startpar_gauss[1] = 0.6;
  startpar_gauss[2] = 0.2;
  gaussfunc->SetParameters(startpar_gauss);
  
  for( int i=1; i<=nnodes_HCAL; i++ ){
    gaussfunc->SetParameters(startpar_gauss);
    
    int page = (i-1)/25;
    int panel = (i-1)%25+1;
    int node = int( hnodesums_HCAL->GetXaxis()->GetBinCenter(i) );
    htemp = hnodesums_HCAL->ProjectionY("",i,i);
    
    TString title;
    title.Form("HCAL node %d;Num Photoelectrons;",node);
    htemp->SetTitle(title.Data());
    double xmin = 0.1, xmax = 2.0;
    int binmin = htemp->FindBin(xmin);
    int binmax = htemp->FindBin(xmax);
    htemp->GetXaxis()->SetRange(binmin,binmax);

    cout << "panel = " << panel << endl;

    ch->cd(panel);
    ch->Update();

    int maxbin = htemp->GetMaximumBin();
    double max = htemp->GetBinContent(maxbin);

    gaussfunc->SetParameter(0,0.4*max);
    
    int binlo = maxbin, binhi=maxbin;
    while( htemp->GetBinContent(binlo--) >= 0.4*max&&binlo >= binmin ){};
    while( htemp->GetBinContent(binhi++) >= 0.05*max&&binhi <= binmax ){};

    double xlo = htemp->GetBinCenter(binlo);
    double xhi = htemp->GetBinCenter(binhi);

    xlo = 0.25;
    xhi = 1.5;
    
    //   xlo = 750.;
    //xhi = 3000.;
    
    bool okay = false;

    while( !okay ){
      okay = true;
      htemp->Fit( gaussfunc, "SQ", "", xlo, xhi );
      htemp->DrawCopy();
      ch->Update();
      gPad->Modified();
      gSystem->ProcessEvents();
    }
    nodes_HCAL[i-1] = node;
    exh[i-1] = 0.0;
    mean_HCAL[i-1] = gaussfunc->GetParameter(1);
    dmean_HCAL[i-1] = gaussfunc->GetParError(1);
    sigma_HCAL[i-1] = gaussfunc->GetParameter(2);
    dsigma_HCAL[i-1] = gaussfunc->GetParError(2);

    if( panel == 25 || i == nnodes_HCAL ){
      ctemp->Update();
      if( page == 0 ){
	ch->Print("HCAL_sums.pdf(","pdf");
      } else {
	ch->Print("HCAL_sums.pdf","pdf");
      }
      ch->Clear();
      ch->Divide(5,5,.001,.001);
    }
  }
  
  ofstream hcalfile( outfilename_HCAL );
  
  for( int inode=0; inode<nnodes_HCAL; inode++ ){
    TString line;
    line.Form( "%5d %15.6g %15.6g", int(nodes_HCAL[inode]), mean_HCAL[inode], sigma_HCAL[inode] );
    hcalfile << line.Data() << endl;
  }
  
  // hnodesums_HCAL_py->Fit(gaussfunc,"","",500.0,2500.0);

  // TString line;
  // line.Form( "%s   %15.6g  %15.6g", "#ALL nodes", gaussfunc->GetParameter(1), gaussfunc->GetParameter(2) );

  // hcalfile << line << endl;

  TGraphErrors *gmean_HCAL = new TGraphErrors(nnodes_HCAL, nodes_HCAL, mean_HCAL, exh, dmean_HCAL );
  TGraphErrors *gsigma_HCAL = new TGraphErrors(nnodes_HCAL, nodes_HCAL, sigma_HCAL, exh, dsigma_HCAL );

  gmean_HCAL->SetMarkerStyle(20);
  gsigma_HCAL->SetMarkerStyle(20);

  ch->Clear();
  ch->Divide(2,1,.001,.001);

  ch->cd(1);
  gmean_HCAL->Draw("AP");
  ch->cd(2);
  gsigma_HCAL->Draw("AP");
  
  ch->Print("HCAL_sums.pdf)","pdf");
}
