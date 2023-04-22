#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "G4SBSRunData.hh"
#include "WAPP_pythia_tree.C"
#include "TH1D.h"
#include "TH2D.h"
#include "TCut.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TEventList.h"
#include "TF2.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set> 
#include "TObjArray.h"
#include "TChainElement.h"
#include "TRandom3.h"


void WAPP_trigrates_pythia( const char *configfilename, const char *outputfilename ){

  TRandom3 num(0);
  ifstream configfile( configfilename );

  TFile *fout = new TFile(outputfilename,"RECREATE");

  TH1D *hBBPS = new TH1D("hBBPS","BigBite preshower energy sum",150,0.0,1.0 );
  TH1D *hBBSH = new TH1D("hBBSH","BigBite shower energy sum",150,0.0,2.0 );
  TH1D *hBBSHcutPS = new TH1D("hBBSHcutPS","",150,0.0,2.0);
  TH1D *hBBSHanticutPS = new TH1D("hBBSHanticutPS","",150,0.0,2.0);
  TH1D *hBBPScutSH = new TH1D("hBBPScutSH","",150,0.0,1.0);
  TH1D *hBBPSanticutSH = new TH1D("hBBPSanticutSH","",150,0.0,1.0);
  TH1D *hBBTotalShower = new TH1D("bBBTotalShower", "BigBite PS+SH energy sum",150,0.0,2.0);
  TH2D *hBBPS_vs_SH = new TH2D("hBBPS_vs_SH","Preshower vs. Shower", 100, 0.0,2.0, 100, 0.0, 1.0 );
  
  TH1D *hHCAL = new TH1D("hHCAL","HCAL scint. energy sum", 150, 0.0, 0.5 );
  
  TChain *Cproton = new TChain("T");
  TChain *Cneutron = new TChain("T");
  
  TFile *ftemp;
  G4SBSRunData *rdtemp;

  double Ngen_proton = 0.0;
  double Ngen_neutron = 0.0;

  double sigmatot_proton, sigmatot_neutron;
  
  //assuming LD2 target

  double lumi_proton, lumi_neutron;
  
  TString currentline;

  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlistp") ){
    if( !currentline.BeginsWith("#") ) Cproton->Add(currentline);
  }

  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlistn") ){
    if( !currentline.BeginsWith("#") ) Cneutron->Add(currentline);
  }

  TObjArray *protonfiles = Cproton->GetListOfFiles();

  for( int i=0; i<protonfiles->GetEntries(); i++ ){
    TString fname = ( (TChainElement*) (*protonfiles)[i] )->GetTitle();

    ftemp = new TFile(fname,"READ");

    ftemp->GetObject("run_data", rdtemp );

    Ngen_proton += rdtemp->fNtries;
    lumi_proton = rdtemp->fLuminosity;

    ftemp->Close();
    ftemp->Delete();
  }

  TObjArray *neutronfiles = Cneutron->GetListOfFiles();

  for( int i=0; i<neutronfiles->GetEntries(); i++ ){
    TString fname = ( (TChainElement*) (*neutronfiles)[i] )->GetTitle();

    ftemp = new TFile(fname,"READ");

    ftemp->GetObject("run_data", rdtemp );

    Ngen_neutron += rdtemp->fNtries;
    lumi_neutron = rdtemp->fLuminosity;

    ftemp->Close();
    ftemp->Delete();
  }

  cout << "Ngen_p, Ngen_n, Lumi_p, Lumi_n = " << Ngen_proton << ", " << Ngen_neutron << ", " << lumi_proton << ", " << lumi_neutron
       << endl;

  TCut cut="";
  
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    cut += currentline;
  }

  double thresh_HCAL = 0.08;
  configfile >> thresh_HCAL;

  double thresh_BBSH = 0.5;
  configfile >> thresh_BBSH;

  double Emin_BBPS = 0.03;
  configfile >> Emin_BBPS;
  double Emax_BBPS = 0.1;
  configfile >> Emax_BBPS;
  
  double Ibeam=5.0e-6;
  configfile >> Ibeam;
  
  double Ltgt=15.0; //cm
  configfile >> Ltgt;

  double rho_D2 = 0.1638;
  double Mmol_D2 = 2.0141; //g/mol

  double Z_D2 = 1.0;
  double N_D2 = 1.0;
  
  int userad=0;
  configfile >> userad;

  double radthick=0.0;
  configfile >> radthick; //assumed to be Cu.

  double coin_window = 30e-9;
  configfile >> coin_window;
  
  double radlength_Cu = 1.436; //cm

  double thick_Cu = radthick * radlength_Cu; //cm

  double Z_Cu = 29.0;
  double N_Cu = 0.69*(63.-Z_Cu)+0.31*(65.-Z_Cu); //natural isotopic abundances

  double rho_Cu = 8.960; //g/cm^3

  double Mmol_Cu = 63.546; //g/mole;

  double N_A = 6.022e23; //Avogadro

  double nuclei_per_cm2_D2 = rho_D2 * Ltgt / Mmol_D2 * N_A;
  double protons_per_cm2_D2 = nuclei_per_cm2_D2 * Z_D2;
  double neutrons_per_cm2_D2 = nuclei_per_cm2_D2 * N_D2;

  double nuclei_per_cm2_Cu = rho_Cu* thick_Cu / Mmol_Cu * N_A;

  double protons_per_cm2_Cu = nuclei_per_cm2_Cu * Z_Cu;
  double neutrons_per_cm2_Cu = nuclei_per_cm2_Cu * N_Cu;

  double electrons_per_sec = Ibeam/1.602e-19;

  lumi_proton = (protons_per_cm2_D2 + protons_per_cm2_Cu )* electrons_per_sec;
  lumi_neutron = (neutrons_per_cm2_D2 + neutrons_per_cm2_Cu ) * electrons_per_sec;
  
  cout << "lumi proton = " << lumi_proton << endl;
  cout << "lumi neutron = " << lumi_neutron << endl;

  cout << "lumi total = " << lumi_proton + lumi_neutron << endl;
  
  double rate_BB_proton = 0.0;
  double rate_HCAL_proton = 0.0;

  double rate_eBB_proton = 0.0;
  double rate_ecoin_proton = 0.0;
  
  TEventList *elist_p = new TEventList("elist_p");

  Cproton->Draw(">>elist_p",cut);
  
  WAPP_pythia_tree *Tproton = new WAPP_pythia_tree(Cproton);

  long neventp = 0;

  double rate_coin_proton = 0.0;

  int count_realcoin_p = 0;
  
  while( Cproton->GetEntry(elist_p->GetEntry(neventp++)) ){
    if( neventp % 1000 == 0 ) cout << neventp << endl;

    double weight = Tproton->primaries_Sigma * lumi_proton/Ngen_proton;

    double EPS = num.Gaus(200.0*Tproton->Earm_BBPSTF1_det_esum,sqrt(fabs(200.0*Tproton->Earm_BBPSTF1_det_esum) ))/200.0;
    double ESH = num.Gaus(200.0*Tproton->Earm_BBSHTF1_det_esum,sqrt(fabs(200.0*Tproton->Earm_BBSHTF1_det_esum) ) )/200.0;
    
    if( Tproton->Earm_BBPSTF1_hit_nhits>0 ) {
      hBBPS->Fill( EPS, weight );
      if( ESH >= thresh_BBSH ){
	hBBPScutSH->Fill( EPS, weight );
      } else {
	hBBPSanticutSH->Fill( EPS, weight );
      }
    }
      
    if( Tproton->Earm_BBSHTF1_hit_nhits>0 ){
      hBBSH->Fill( ESH, weight );
      if( EPS >= Emin_BBPS && EPS <= Emax_BBPS ){
	hBBSHcutPS->Fill( ESH, weight );
      } else {
	hBBSHanticutPS->Fill( ESH, weight );
      }
    }
      
    if( EPS + ESH >0 ) {
      hBBTotalShower->Fill( EPS + ESH, weight );
      hBBPS_vs_SH->Fill( ESH, EPS, weight );
    } 
    
    if( Tproton->Harm_HCalScint_det_esum > 0 ) hHCAL->Fill( Tproton->Harm_HCalScint_det_esum, weight );

    if( EPS + ESH >= thresh_BBSH ){
      rate_eBB_proton += weight;
      if( Tproton->Harm_HCalScint_det_esum >= thresh_HCAL ){
	rate_ecoin_proton += weight;
      }
    }
    
    if( EPS >= Emin_BBPS &&
	EPS <= Emax_BBPS &&
	ESH >= thresh_BBSH ){
    
      rate_BB_proton += weight;

      if( Tproton->Harm_HCalScint_det_esum >= thresh_HCAL ){
	rate_coin_proton += weight;
	count_realcoin_p++;
      }
      
    }

    if( Tproton->Harm_HCalScint_det_esum >= thresh_HCAL ){
      rate_HCAL_proton += weight;
    }
  }

  long neventn=0;
  
  TEventList *elist_n = new TEventList("elist_n");

  Cneutron->Draw(">>elist_n",cut);
  
  WAPP_pythia_tree *Tneutron = new WAPP_pythia_tree(Cneutron);

  double rate_BB_neutron = 0.0;
  double rate_HCAL_neutron = 0.0;
  double rate_coin_neutron = 0.0;

  double rate_eBB_neutron = 0.0;
  double rate_ecoin_neutron = 0.0;
  
  int count_realcoin_n=0;
  
  while( Cneutron->GetEntry(elist_n->GetEntry(neventn++)) ){
    if( neventn % 1000 == 0 ) cout << neventn << endl;

    double weight = Tneutron->primaries_Sigma * lumi_neutron/Ngen_neutron;

    double EPS = num.Gaus(200.0*Tneutron->Earm_BBPSTF1_det_esum,sqrt(fabs(200.0*Tneutron->Earm_BBPSTF1_det_esum) ))/200.0;
    double ESH = num.Gaus(200.0*Tneutron->Earm_BBSHTF1_det_esum,sqrt(fabs(200.0*Tneutron->Earm_BBSHTF1_det_esum) ) )/200.0;
    
    if( Tneutron->Earm_BBPSTF1_hit_nhits>0 ) {
      hBBPS->Fill( EPS, weight );
      if( ESH >= thresh_BBSH ){
	hBBPScutSH->Fill( EPS, weight );
      } else {
	hBBPSanticutSH->Fill( EPS, weight );
      }
    }
      
    if( Tneutron->Earm_BBSHTF1_hit_nhits>0 ){
      hBBSH->Fill( ESH, weight );
      if( EPS >= Emin_BBPS && EPS <= Emax_BBPS ){
	hBBSHcutPS->Fill( ESH, weight );
      } else {
	hBBSHanticutPS->Fill( ESH, weight );
      }
    }
    if( EPS + ESH >0 ) {
      hBBTotalShower->Fill( EPS + ESH, weight );
      hBBPS_vs_SH->Fill( ESH, EPS, weight );
    }
    
    if( Tneutron->Harm_HCalScint_det_esum > 0 ) hHCAL->Fill( Tneutron->Harm_HCalScint_det_esum, weight );

    if( EPS + ESH >= thresh_BBSH ){
      rate_eBB_neutron += weight;
      if( Tneutron->Harm_HCalScint_det_esum >= thresh_HCAL ){
	rate_ecoin_neutron += weight;
      }
    }
    
    if( EPS >= Emin_BBPS &&
	EPS <= Emax_BBPS &&
	ESH >= thresh_BBSH ){
    
      rate_BB_neutron += Tneutron->primaries_Sigma * lumi_neutron/Ngen_neutron;

      if( Tneutron->Harm_HCalScint_det_esum >= thresh_HCAL ){
	rate_coin_neutron += weight;

	count_realcoin_n++;
      }
      
    }

    if( Tneutron->Harm_HCalScint_det_esum >= thresh_HCAL ){
      rate_HCAL_neutron += Tneutron->primaries_Sigma * lumi_neutron/Ngen_neutron;
    }
  }

  cout << "For a beam current of Ibeam = " << Ibeam << "uA and Cu radiator thickness of " << radthick << "X0, with LD2 target thickness = "
       << Ltgt << " cm:" << endl;
  cout << "Total rate in BB due to protons in target + rad. = " << rate_BB_proton << endl;
  cout << "Total rate in BB due to neutrons in target + rad. = " << rate_BB_neutron << endl;
  cout << "Total rate in BB = " << rate_BB_proton + rate_BB_neutron << endl;
  
  cout << "Total rate in HCAL due to protons in target + rad. = " << rate_HCAL_proton << endl;
  cout << "Total rate in HCAL due to neutrons in target + rad. = " << rate_HCAL_neutron << endl;
  cout << "Total rate in HCAL = " << rate_HCAL_proton + rate_HCAL_neutron << endl;

  cout << "Assuming a coincidence time window of " << coin_window*1e9 << " ns, accidental coincidence rate = "
       << (rate_HCAL_proton + rate_HCAL_neutron)*(rate_BB_proton + rate_BB_neutron)*coin_window << endl;

  cout << "Real coincidence rate due to protons in target + rad. = " << rate_coin_proton << endl;
  cout << "Real coincidence rate due to neutrons in target + rad. = " << rate_coin_neutron << endl;
  cout << "Total real coincidence rate = " << rate_coin_proton + rate_coin_neutron << endl;

  cout << "Number of real coincidence events on proton in this PYTHIA sample = " << count_realcoin_p << endl;
  cout << "Number of real coincidence events on neutron in this PYTHIA sample = " << count_realcoin_n << endl;

  cout << "Using 'electron' logic (threshold on sum of PS+SH): " << endl;
  cout << "Total rate in BB = " << rate_eBB_proton + rate_eBB_neutron << endl;
  cout << "Real coincidence rate = " << rate_ecoin_proton + rate_ecoin_neutron << endl;
  cout << "Accidental coincidence rate = " << (rate_eBB_proton + rate_eBB_neutron)*(rate_HCAL_proton + rate_HCAL_neutron)*coin_window << endl;
  
  elist_p->Delete();
  elist_n->Delete();
  fout->Write();
  

}
