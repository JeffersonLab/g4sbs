#include "G4SBSRunData.hh"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObjArray.h"
#include "TChainElement.h"
#include "TClonesArray.h"
#include <vector>
#include <iostream>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include "TGraph.h"
//#include "gmn_tree.h"
#include "gmn_tree.C"
#include "TRandom3.h"

const double me = 0.511e-3;

void background_rates_GMN( const char *setupfilename, const char *outfilename ){

  TRandom3 num(0);
  
  ifstream setupfile(setupfilename);

  TString filename;

  TChain *C = new TChain("T");
  
  while( setupfile >> filename && !filename.BeginsWith("endlist")){
    if( !filename.BeginsWith("#") ){
      C->Add(filename);
    }
  }

  int pheflag = 0;
  setupfile >> pheflag;

  double mean_npeGeV=300.0;
  double rindex_PbGl=1.68;

  setupfile >> mean_npeGeV;
  setupfile >> rindex_PbGl;

  int generator_flag = 0; //0 = beam, 1 = pythia, 2 = other:
  setupfile >> generator_flag; 

  cout << "mean npe/GeV BBCAL = " << mean_npeGeV << endl;
  cout << "rindex Pb-Gl = " << rindex_PbGl << endl;
  cout << "generator flag = " << generator_flag << endl;
  
  G4SBSRunData *rd;
  
  // C->Add(rootfilename);

  map<TString,double> Lumi_file;
  map<TString,double> genvol_file; 
  
  double Lumi_default = 3.97e36; //1 uA on 15-cm LH2
  double genvol_default = 1.0;
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
      Lumi_file[chEl->GetTitle()] = rd->fLuminosity;
      genvol_file[chEl->GetTitle()] = rd->fGenVol;
      Lumi_default = rd->fLuminosity;
      genvol_default = rd->fGenVol;
      nfiles++;
    } else {
      bad_file_list.insert( chEl->GetTitle());
    }
    //cout << chEl->GetTitle() << endl;
  }

  cout << "Total number of generated events = " << ngen << endl;

  gmn_tree *T = new gmn_tree(C);

  TFile *fout = new TFile( outfilename, "RECREATE" );

  TH2D *hnphe_vs_edep_HCAL = new TH2D("hnphe_vs_edep_HCAL","",500,0.0,0.5,501,-0.5,500.5);
  TH2D *hrate_vs_nphe_HCAL = new TH2D("hrate_vs_nphe_HCAL","",288,-0.5,287.5,500,0.5,500.5);
  TH2D *hrate_vs_edep_HCALscint = new TH2D("hrate_vs_edep_HCALscint","",288,-0.5,287.5,500,0.0,1.0); //GeV

  TH2D *hnphe_vs_edep_CDET = new TH2D("hnphe_vs_edep_CDET","",500,0.0,0.25,501,-0.5,500.5);
  TH2D *hrate_vs_nphe_CDET = new TH2D("hrate_vs_nphe_CDET","",2352,-0.5,2351.5,500,0.5,500.5);
  TH2D *hrate_vs_edep_CDETscint = new TH2D("hrate_vs_edep_CDETscint","",2352,-0.5,2351.5,500,0.0,1.0); //GeV

  //TH2D *hnphe_vs_edep_BBhodo = new TH2D("hnphe_vs_edep_BBhodo","",500,0.0,0.25,501,-0.5,500.5);
  //TH2D *hrate_vs_nphe_BBhodo = new TH2D("hrate_vs_nphe_BBhodo","",2352,-0.5,2351.5,500,0.5,500.5);
  TH2D *hrate_vs_edep_BBHodoScint = new TH2D("hrate_vs_edep_BBHodoScint","",90,-0.5,89.5,500,0.0,0.15); //GeV

  TH2D *hnphe_vs_edep_BBPS = new TH2D("hnphe_vs_edep_BBPS","",500,0.0,0.25,501,-0.5,500.5);
  TH2D *hrate_vs_nphe_BBPS = new TH2D("hrate_vs_nphe_BBPS","",54,-0.5,53.5,501,-0.5,500.5);
  TH2D *hrate_vs_edep_BBPSTF1 = new TH2D("hrate_vs_edep_BBPSTF1","",54,-0.5,53.5,500,0.0,0.25); //GeV

  TH2D *hnphe_vs_edep_BBSH = new TH2D("hnphe_vs_edep_BBSH","",500,0.0,0.5,501,-0.5,500.5);
  TH2D *hrate_vs_nphe_BBSH = new TH2D("hrate_vs_nphe_BBSH","",189,-0.5,188.5,501,-0.5,500.5);
  TH2D *hrate_vs_edep_BBSHTF1 = new TH2D("hrate_vs_edep_BBSHTF1","",189,-0.5,188.5,500,0.0,0.5); //GeV

  TH2D *hrate_vs_nphe_BBPS_from_edep = new TH2D("hrate_vs_nphe_BBPS_from_edep","",54,-0.5,53.5,501,-0.5,500.5);
  TH2D *hrate_vs_nphe_BBSH_from_edep = new TH2D("hrate_vs_nphe_BBSH_from_edep","",189,-0.5,188.5,501,-0.5,500.5);

  TH1D *hitrate_vs_layer_BBGEM = new TH1D("hitrate_vs_layer_BBGEM","",5,0.5,5.5);
  TH1D *hrate_edep_vs_layer_BBGEM = new TH1D("hrate_edep_vs_layer_BBGEM","BigBite GEM primary ionization region; GEM layer; Energy deposit rate (MeV/cm^{2}/s)", 5,0.5,5.5); 
  TH2D *hitrate_vs_X_BBGEM = new TH2D("hitrate_vs_X_BBGEM","Hit rate (Hz/cm^{2}/uA), 15-cm LD2",5,0.5,5.5,100,-1.05,1.05);
  TH2D *hitrate_vs_Y_BBGEM = new TH2D("hitrate_vs_Y_BBGEM","Hit rate (Hz/cm^{2}/uA), 15-cm LD2",5,0.5,5.5,100,-0.31,0.31);

  TH1D *hrate_vs_esum_BBCAL = new TH1D("hrate_vs_esum_BBCAL","event rate vs shower + preshower edep.",400,0.0,4.0);
  hrate_vs_esum_BBCAL->SetXTitle("BB SH + PS smeared edep (GeV)");
  hrate_vs_esum_BBCAL->SetYTitle("Event rate on 15-cm LD2 (Hz/(10 MeV)/uA)");

  TH1D *hrate_vs_esum_HCAL = new TH1D("hrate_vs_esum_HCAL","event rate vs HCAL scint total edep.",1000,0.0,1.0);
  hrate_vs_esum_HCAL->SetXTitle("HCal summed energy deposit (GeV)");
  hrate_vs_esum_HCAL->SetYTitle("Event rate on 15-cm LD2 (Hz/MeV/uA)");
  
  hitrate_vs_X_BBGEM->SetXTitle("BigBite GEM layer");
  hitrate_vs_X_BBGEM->SetYTitle("Hit X (m)");

  hitrate_vs_Y_BBGEM->SetXTitle("BigBite GEM layer");
  hitrate_vs_Y_BBGEM->SetYTitle("Hit Y (m)");
  
  double BBGEM_area_cm2[5] = {40.0*150.0, 40.0*150.0, 40.0*150.0, 40.0*150.0, 60.0*200.0 };
  double BBGEM_LX[5] = {150.,150.,150.,150.,200.};
  double BBGEM_LY[5] = {40.,40.,40.,40.,60.};

  TClonesArray *BBGEM_edep = new TClonesArray("TH1D",5);
  TClonesArray *BBGEM_logedep_keV = new TClonesArray("TH1D",5);
  TClonesArray *BBGEM_yvsx = new TClonesArray("TH2D",5);
  TClonesArray *BBGEM_dyvsdx = new TClonesArray("TH2D",5);
  

  for( int layer=0; layer<5; layer++ ){
    new( (*BBGEM_edep)[layer] ) TH1D( Form("h1_BBGEM_Edep_%d",layer), Form("BB GEM layer %d; Energy deposit (MeV)",layer), 500,0,0.1);
    new( (*BBGEM_logedep_keV)[layer] ) TH1D( Form("h1_BBGEM_Edep_log_%d",layer), Form("BB GEM layer %d;log(Energy deposit (keV))",layer), 200, -5.,5.);
    new( (*BBGEM_yvsx)[layer] ) TH2D( Form("h1_BBGEM_yVsx_%d",layer), Form("BB GEM layer %d; x hit (m); y hit (m)",layer), 105, -1.05, 1.05, 36, -0.36, 0.36 );
    new( (*BBGEM_dyvsdx)[layer] ) TH2D( Form("h1_BBGEM_dyVsdx_%d",layer), Form("BB GEM layer %d; dx/dz ; dy/dz",layer), 100,-0.03,0.03,100,-0.03,0.03);
  }
  
  
  double pmtnum[288];
  double sumedep_HCAL[288];
  // double sum2edep_HCAL[288];
  double sumnphe_HCAL[288];
  double hitrate_HCAL[288];
  
  for( int ipmt=0; ipmt<288; ipmt++ ){
    sumedep_HCAL[ipmt] = 0.0;
    pmtnum[ipmt] = ipmt;
    sumnphe_HCAL[ipmt] = 0.0;
    hitrate_HCAL[ipmt] = 0.0;
  }

  double pmt_CDET[2352];
  double sumedep_CDET[2352];
  double sumnphe_CDET[2352];
  double hitrate_CDET[2352];
  
  for( int ipmt=0; ipmt<2352; ipmt++ ){
    pmt_CDET[ipmt] = ipmt;
    sumedep_CDET[ipmt] = 0.0;
    sumnphe_CDET[ipmt] = 0.0;
    hitrate_CDET[ipmt] = 0.0;
  }

  double pmtBBHodo[90];
  double sumedep_BBHodo[90];
  double hitrate_BBHodo[90];
  
  for( int ipmt=0; ipmt<90; ipmt++ ){
    pmtBBHodo[ipmt] = ipmt;
    sumedep_BBHodo[ipmt] = 0.0;
    hitrate_BBHodo[ipmt] = 0.0;
  }

  double PMT_PS[2*27];
  double PMT_SH[189];

  double sumedep_PS[54];
  double sumedep_SH[189];
  double sumnphe_PS[54];
  double sumnphe_SH[189];
  double hitrate_PS[54];
  double hitrate_SH[189];
  
  for( int i=0; i<54; i++ ){
    PMT_PS[i] = i;
    sumedep_PS[i] = 0.0;
    sumnphe_PS[i] = 0.0;
    hitrate_PS[i] = 0.0;
  }

  for( int i=0; i<189; i++ ){
    PMT_SH[i] = i;
    sumedep_SH[i] = 0.0;
    sumnphe_SH[i] = 0.0;
    hitrate_SH[i] = 0.0;
  }

  
  double thresh_CDET = .0055; //5.5 MeV
  double thresh_BBhodo = 0.008; //3 MeV
  // double thresh_BBPS   = 0.02*0.751; //7.5 MeV
  // double thresh_BBSH   = 0.02*2.83; //28.3 MeV
  double thresh_BBPS   = 0.05;
  double thresh_BBSH   = 0.05;
  double thresh_HCAL   = 0.02*0.56; //~11 MeV
  
  long nevent=0;

  double Ibeam = 1e-6; //A
  double weight = Ibeam/double(ngen)/1.602e-19;

  double Xbinwidth = 210.0/100.0;
  double Ybinwidth = 62.0/100.0;
  
  while( T->GetEntry( nevent++ ) ){
    if( nevent % 1000 == 0 ) cout << nevent << endl;
    
    TFile *f = ( (TChain*) (T->fChain) )->GetFile();

    if( bad_file_list.find( f->GetName() ) == bad_file_list.end() ){

      if( generator_flag == 1 ){ //PYTHIA:
	weight = Lumi_file[f->GetName()] * T->primaries_Sigma / double(ngen);

	//cout << "lumi, sigma, weight = " << Lumi_file[f->GetName()] << ", " << T->primaries_Sigma << ", " << weight << endl;
      }
      if( generator_flag == 2 ){ //Generic g4sbs "built-in" generators:
	weight = Lumi_file[f->GetName()] * T->ev_sigma * genvol_file[f->GetName()] / double(ngen);	
      }
      
      //Add GEMs:
      for( int ihit=0; ihit<T->Earm_BBGEM_hit_nhits; ihit++ ){
	int plane = (*(T->Earm_BBGEM_hit_plane))[ihit];

	hitrate_vs_layer_BBGEM->Fill( plane, weight/BBGEM_area_cm2[plane-1] );
	hitrate_vs_X_BBGEM->Fill( plane, (*(T->Earm_BBGEM_hit_x))[ihit], weight/Xbinwidth/BBGEM_LY[plane-1] );
	hitrate_vs_Y_BBGEM->Fill( plane, (*(T->Earm_BBGEM_hit_y))[ihit], weight/Ybinwidth/BBGEM_LX[plane-1] );

	double edep = (*(T->Earm_BBGEM_hit_edep))[ihit];

	hrate_edep_vs_layer_BBGEM->Fill( plane, weight * edep * 1000.0 / BBGEM_area_cm2[plane-1] );

	double xhit = (*(T->Earm_BBGEM_hit_x))[ihit];
	double yhit = (*(T->Earm_BBGEM_hit_y))[ihit];
	//double dxhit = (*(T->Earm_BBGEM_hit_txp))[ihit];
	//double dyhit = (*(T->Earm_BBGEM_hit_typ))[ihit];
	double dxhit = (*(T->Earm_BBGEM_hit_xout))[ihit]-(*(T->Earm_BBGEM_hit_xin))[ihit];
	double dyhit = (*(T->Earm_BBGEM_hit_yout))[ihit]-(*(T->Earm_BBGEM_hit_yin))[ihit];
	double logedep_keV = log( edep*1.e6 );
	
	( (TH1D*) (*BBGEM_edep)[plane-1] )->Fill( edep *1000.0);
	( (TH1D*) (*BBGEM_logedep_keV)[plane-1] )->Fill( logedep_keV );
	( (TH2D*) (*BBGEM_yvsx)[plane-1] )->Fill( xhit, yhit );
	( (TH2D*) (*BBGEM_dyvsdx)[plane-1] )->Fill( dxhit, dyhit );
	
      }
      
      // if( pheflag != 0 ){ //optical photons ON: only consider hits with optical photons detected in PMTs
      // 	for( int ihit=0; ihit<T->Harm_HCal_hit_nhits; ihit++ ){
      // 	  int PMT = (*(T->Harm_HCal_hit_PMT))[ihit];
      // 	  int nphe = (*(T->Harm_HCal_hit_NumPhotoelectrons))[ihit];

      // 	  double edep = 0.0;
	
      // 	  for( int jhit=0; jhit<T->Harm_HCalScint_hit_nhits; jhit++ ){
      // 	    if( (*(T->Harm_HCalScint_hit_cell))[jhit] == PMT ){
      // 	      edep = (*(T->Harm_HCalScint_hit_sumedep))[jhit];
      // 	    }
      // 	  }
      // 	  hrate_vs_nphe_HCAL->Fill( PMT, nphe, weight );
      // 	  hrate_vs_edep_HCALscint->Fill( PMT, edep, weight );

      // 	  hnphe_vs_edep_HCAL->Fill( edep, nphe, weight );
	
      // 	  sumedep_HCAL[PMT] += edep * weight;
      // 	  sumnphe_HCAL[PMT] += nphe * weight;

      // 	  if( edep >= thresh_HCAL ) hitrate_HCAL[PMT] += weight;
      // 	}

      // 	// for( int ihit=0; ihit<T->Harm_CDET_hit_nhits; ihit++ ){
      // 	//   int PMT = (*(T->Harm_CDET_hit_PMT))[ihit];
      // 	//   int nphe = (*(T->Harm_CDET_hit_NumPhotoelectrons))[ihit];

      // 	//   double edep = 0.0;
	
      // 	//   for( int jhit=0; jhit<T->Harm_CDET_Scint_hit_nhits; jhit++ ){
      // 	//     if( (*(T->Harm_CDET_Scint_hit_cell))[jhit] == PMT ){
      // 	//       edep = (*(T->Harm_CDET_Scint_hit_sumedep))[jhit];
      // 	//     }
      // 	//   }
      // 	//   hrate_vs_nphe_CDET->Fill( PMT, nphe, weight );
      // 	//   hrate_vs_edep_CDETscint->Fill( PMT, edep, weight );

      // 	//   hnphe_vs_edep_CDET->Fill( edep, nphe, weight );
	
      // 	//   sumedep_CDET[PMT] += edep * weight;
      // 	//   sumnphe_CDET[PMT] += nphe * weight;

      // 	//   if( edep >= thresh_CDET ) hitrate_CDET[PMT] += weight;
      // 	// }

      // 	for( int ihit=0; ihit<T->Earm_BBHodoScint_hit_nhits; ihit++ ){
      // 	  double edep = (*(T->Earm_BBHodoScint_hit_sumedep))[ihit];
      // 	  int PMT = (*(T->Earm_BBHodoScint_hit_cell))[ihit];
      // 	  hrate_vs_edep_BBHodoScint->Fill( PMT, edep, weight );

      // 	  sumedep_BBHodo[PMT] += edep * weight;

      // 	  if( edep >= thresh_BBhodo ) hitrate_BBHodo[PMT] += weight;
      // 	}

      // 	//BB PS:
      // 	for( int ihit=0; ihit<T->Earm_BBPS_hit_nhits; ihit++ ){
      // 	  int PMT = (*(T->Earm_BBPS_hit_PMT))[ihit];
      // 	  int nphe = (*(T->Earm_BBPS_hit_NumPhotoelectrons))[ihit];
      // 	  double edep = 0.0;
      // 	  for( int jhit=0; jhit<T->Earm_BBPSTF1_hit_nhits; jhit++ ){
      // 	    int cell = (*(T->Earm_BBPSTF1_hit_cell))[jhit];
      // 	    if( cell == PMT ){
      // 	      edep = (*(T->Earm_BBPSTF1_hit_sumedep))[jhit];
      // 	    }
      // 	  }

      // 	  hrate_vs_nphe_BBPS->Fill( PMT, nphe, weight );
      // 	  hrate_vs_edep_BBPSTF1->Fill( PMT, edep, weight );

      // 	  hnphe_vs_edep_BBPS->Fill( edep, nphe, weight );

      // 	  sumnphe_PS[PMT] += nphe * weight;
      // 	  sumedep_PS[PMT] += edep * weight;

      // 	  if( edep >= thresh_BBPS ) hitrate_PS[PMT] += weight;
      // 	}

      // 	//BB SH:
      // 	for( int ihit=0; ihit<T->Earm_BBSH_hit_nhits; ihit++ ){
      // 	  int PMT = (*(T->Earm_BBSH_hit_PMT))[ihit];
      // 	  int nphe = (*(T->Earm_BBSH_hit_NumPhotoelectrons))[ihit];
      // 	  double edep = 0.0;
      // 	  for( int jhit=0; jhit<T->Earm_BBSHTF1_hit_nhits; jhit++ ){
      // 	    int cell = (*(T->Earm_BBSHTF1_hit_cell))[jhit];
      // 	    if( cell == PMT ){
      // 	      edep = (*(T->Earm_BBSHTF1_hit_sumedep))[jhit];
      // 	    }
      // 	  }

      // 	  hrate_vs_nphe_BBSH->Fill( PMT, nphe, weight );
      // 	  hrate_vs_edep_BBSHTF1->Fill( PMT, edep, weight );

      // 	  hnphe_vs_edep_BBSH->Fill( edep, nphe, weight );
	
      // 	  sumnphe_SH[PMT] += nphe * weight;
      // 	  sumedep_SH[PMT] += edep * weight;

      // 	  if( edep >= thresh_BBSH ) hitrate_SH[PMT] += weight;
      // 	}
      
      // } else { // no optical photons:
      //HCAL:
      for( int jhit=0; jhit<T->Harm_HCalScint_hit_nhits; jhit++ ){

	int PMT = (*(T->Harm_HCalScint_hit_cell))[jhit];
	double edep = (*(T->Harm_HCalScint_hit_sumedep))[jhit];

	if( edep >= thresh_HCAL ) hitrate_HCAL[PMT] += weight;

	hrate_vs_edep_HCALscint->Fill( PMT, edep, weight );

	sumedep_HCAL[PMT] += edep * weight;
	  
      }

      hrate_vs_esum_HCAL->Fill( T->Harm_HCalScint_det_esum, weight );

      //CDET:
      // for( int jhit=0; jhit<T->Harm_CDET_Scint_hit_nhits; jhit++ ){
      //   int PMT = (*(T->Harm_CDET_Scint_hit_cell))[jhit];
      //   double edep = (*(T->Harm_CDET_Scint_hit_sumedep))[jhit];

      //   if( edep >= thresh_CDET ) hitrate_CDET[PMT] += weight;

      //   hrate_vs_edep_CDETscint->Fill( PMT, edep, weight );

      //   sumedep_CDET[PMT] += edep * weight;
      // }

      //BB timing hodo:
      for( int ihit=0; ihit<T->Earm_BBHodoScint_hit_nhits; ihit++ ){
	double edep = (*(T->Earm_BBHodoScint_hit_sumedep))[ihit];
	int PMT = (*(T->Earm_BBHodoScint_hit_cell))[ihit];
	hrate_vs_edep_BBHodoScint->Fill( PMT, edep, weight );

	sumedep_BBHodo[PMT] += edep * weight;

	if( edep >= thresh_BBhodo ) hitrate_BBHodo[PMT] += weight;
      }

      //BB PS:
	
      for( int jhit=0; jhit<T->Earm_BBPSTF1_hit_nhits; jhit++ ){
	int PMT = (*(T->Earm_BBPSTF1_hit_cell))[jhit];
	double edep = (*(T->Earm_BBPSTF1_hit_sumedep))[jhit];

	double beta_edep = sqrt(pow(edep,2)+2.0*me*edep)/(me+edep);
	//cout << "beta_edep = " << beta_edep << endl;

	//sin^2 (thetaC) = 1 - 1/(beta*n)^2, 

	double sin2thetaC_edep = (beta_edep > 1./rindex_PbGl ) ? 1.0 - pow( beta_edep * rindex_PbGl, -2 ) : 0.0;
	//if( beta_edep > 1.0/rindex_PbGl )
	
	//double sin2thetaC_edep = 1.0-pow(TMath::Max(beta_edep,1.0/rindex_PbGl)*rindex_PbGl,-2);
	double sin2thetaC_max = 1.0-pow(rindex_PbGl,-2);

	double mean_npe = edep * mean_npeGeV * sin2thetaC_edep/sin2thetaC_max; 
	
	double npe_edep = num.Gaus(mean_npeGeV, sqrt(mean_npeGeV) );

	hrate_vs_nphe_BBPS_from_edep->Fill( PMT, npe_edep, weight );
	  
	if( edep >= thresh_BBPS ) hitrate_PS[PMT] += weight;

	hrate_vs_edep_BBPSTF1->Fill( PMT, edep, weight );

	sumedep_PS[PMT] += edep * weight;
      }

      //BB SH:
	
      for( int jhit=0; jhit<T->Earm_BBSHTF1_hit_nhits; jhit++ ){
	int PMT = (*(T->Earm_BBSHTF1_hit_cell))[jhit];
	double edep = (*(T->Earm_BBSHTF1_hit_sumedep))[jhit];

	double beta_edep = sqrt(pow(edep,2)+2.0*me*edep)/(me+edep);
	// Let edep be the kinetic energy of the electron; then Etot = me + edep
	// p^2 = (me + edep)^2 - me^2 = me^2 - me^2 + edep2 + 2medep
	//cout << "beta_edep = " << beta_edep << endl;

	//sin^2 (thetaC) = 1 - 1/(beta*n)^2, 

	double sin2thetaC_edep = (beta_edep > 1./rindex_PbGl ) ? 1.0 - pow( beta_edep * rindex_PbGl, -2 ) : 0.0;
	//if( beta_edep > 1.0/rindex_PbGl )
	
	//double sin2thetaC_edep = 1.0-pow(TMath::Max(beta_edep,1.0/rindex_PbGl)*rindex_PbGl,-2);
	double sin2thetaC_max = 1.0-pow(rindex_PbGl,-2);

	double mean_npe = edep * mean_npeGeV * sin2thetaC_edep/sin2thetaC_max; 
	
	double npe_edep = num.Gaus(mean_npeGeV, sqrt(mean_npeGeV) );

	hrate_vs_nphe_BBSH_from_edep->Fill( PMT, npe_edep, weight );
	  
	if( edep >= thresh_BBSH ) hitrate_SH[PMT] += weight;

	hrate_vs_edep_BBSHTF1->Fill( PMT, edep, weight );

	sumedep_SH[PMT] += edep * weight;
      }

      double etot_PS = T->Earm_BBPSTF1_det_esum;
      double etot_SH = T->Earm_BBSHTF1_det_esum;

      double npemean_PS = etot_PS * mean_npeGeV;
      double npemean_SH = etot_SH * mean_npeGeV;

      double npesmear_PS = num.Gaus( npemean_PS, sqrt(npemean_PS) );
      double npesmear_SH = num.Gaus( npemean_SH, sqrt(npemean_SH) );

      //Smeared energy:
      double etot_PS_smear = npesmear_PS/mean_npeGeV;
      double etot_SH_smear = npesmear_SH/mean_npeGeV;

      double esum_BBCAL_smear = etot_PS_smear + etot_SH_smear;

      //cout << "esum_BBCAL_smear = " << esum_BBCAL_smear << endl;
      
      hrate_vs_esum_BBCAL->Fill( esum_BBCAL_smear, weight );
	
      
    }
  }

  TGraph *edep_rate_HCAL = new TGraph( 288, pmtnum, sumedep_HCAL );
  TGraph *nphe_rate_HCAL = new TGraph( 288, pmtnum, sumnphe_HCAL );

  edep_rate_HCAL->SetMarkerStyle(20);
  edep_rate_HCAL->GetYaxis()->SetTitle("GeV/s (nphe >= 1)");
  edep_rate_HCAL->GetXaxis()->SetTitle("HCAL PMT number");
  nphe_rate_HCAL->SetMarkerStyle(20);
  nphe_rate_HCAL->GetYaxis()->SetTitle("phe/s");
  nphe_rate_HCAL->GetXaxis()->SetTitle("HCAL PMT number");
  
  edep_rate_HCAL->Write("edep_total_rate_vs_PMT_HCAL");
  nphe_rate_HCAL->Write("nphe_total_rate_vs_PMT_HCAL");

  TGraph *ghitrate_HCAL = new TGraph( 288, pmtnum, hitrate_HCAL );
  ghitrate_HCAL->SetMarkerStyle(20);
  ghitrate_HCAL->GetXaxis()->SetTitle("HCAL PMT number");
  ghitrate_HCAL->GetYaxis()->SetTitle("Hit rate (thresh. 11 MeV)");
  ghitrate_HCAL->Write("hitrate_HCAL");
  
  TGraph *edep_rate_CDET = new TGraph( 2352, pmt_CDET, sumedep_CDET );
  TGraph *nphe_rate_CDET = new TGraph( 2352, pmt_CDET, sumnphe_CDET );

  edep_rate_CDET->SetMarkerStyle(20);
  edep_rate_CDET->GetYaxis()->SetTitle("GeV/s (nphe >= 1)");
  edep_rate_CDET->GetXaxis()->SetTitle("CDET PMT number");
  nphe_rate_CDET->SetMarkerStyle(20);
  nphe_rate_CDET->GetYaxis()->SetTitle("phe/s");
  nphe_rate_CDET->GetXaxis()->SetTitle("CDET PMT number");

  
  edep_rate_CDET->Write("edep_total_rate_vs_PMT_CDET");
  nphe_rate_CDET->Write("nphe_total_rate_vs_PMT_CDET");

  TGraph *ghitrate_CDET = new TGraph( 2352, pmt_CDET, hitrate_CDET );
  ghitrate_CDET->SetMarkerStyle(20);
  ghitrate_CDET->GetXaxis()->SetTitle("CDET PMT number");
  ghitrate_CDET->GetYaxis()->SetTitle("Hit rate (thresh. 5.5 MeV)");
  ghitrate_CDET->Write("hitrate_CDET");
  
  TGraph *edep_rate_BBHodo = new TGraph( 90, pmtBBHodo, sumedep_BBHodo );
  edep_rate_BBHodo->SetMarkerStyle(20);
  edep_rate_BBHodo->GetYaxis()->SetTitle("GeV/s");
  edep_rate_BBHodo->GetXaxis()->SetTitle("BB hodo PMT number");

  edep_rate_BBHodo->Write("edep_rate_BBHodo");

  TGraph *ghitrate_BBHodo = new TGraph( 90, pmtBBHodo, hitrate_BBHodo );
  ghitrate_BBHodo->SetMarkerStyle(20);
  ghitrate_BBHodo->GetXaxis()->SetTitle("BBHodo PMT number");
  ghitrate_BBHodo->GetYaxis()->SetTitle("Hit rate (thresh. 3 MeV)");
  ghitrate_BBHodo->Write("hitrate_BBHodo");
  
  TGraph *edep_rate_PS = new TGraph( 54, PMT_PS, sumedep_PS );
  TGraph *nphe_rate_PS = new TGraph( 54, PMT_PS, sumnphe_PS );

  edep_rate_PS->SetMarkerStyle(20);
  edep_rate_PS->GetYaxis()->SetTitle( "GeV/s (nphe >= 1)" );
  edep_rate_PS->GetXaxis()->SetTitle( "PS PMT number" );
  edep_rate_PS->Write("edep_rate_PS");

  nphe_rate_PS->SetMarkerStyle(20);
  nphe_rate_PS->GetYaxis()->SetTitle( "nphe/s" );
  nphe_rate_PS->GetXaxis()->SetTitle( "PS PMT number" );
  nphe_rate_PS->Write("nphe_rate_PS");

  TGraph *ghitrate_PS = new TGraph( 54, PMT_PS, hitrate_PS );
  ghitrate_PS->SetMarkerStyle(20);
  ghitrate_PS->GetXaxis()->SetTitle("PS PMT number");
  ghitrate_PS->GetYaxis()->SetTitle("Hit rate (thresh. 7.5 MeV)");
  ghitrate_PS->Write("hitrate_PS");
  
  TGraph *edep_rate_SH = new TGraph( 189, PMT_SH, sumedep_SH );
  TGraph *nphe_rate_SH = new TGraph( 189, PMT_SH, sumnphe_SH );

  edep_rate_SH->SetMarkerStyle(20);
  edep_rate_SH->GetYaxis()->SetTitle( "GeV/s (nphe >= 1)" );
  edep_rate_SH->GetXaxis()->SetTitle( "SH PMT number" );
  edep_rate_SH->Write("edep_rate_SH");

  nphe_rate_SH->SetMarkerStyle(20);
  nphe_rate_SH->GetYaxis()->SetTitle( "nphe/s" );
  nphe_rate_SH->GetXaxis()->SetTitle( "SH PMT number" );
  nphe_rate_SH->Write("nphe_rate_SH");

  TGraph *ghitrate_SH = new TGraph( 189, PMT_SH, hitrate_SH );
  ghitrate_SH->SetMarkerStyle(20);
  ghitrate_SH->GetXaxis()->SetTitle("SH PMT number");
  ghitrate_SH->GetYaxis()->SetTitle("Hit rate (thresh. 28.3 MeV)");
  ghitrate_SH->Write("hitrate_SH");

  TH1D *hrate_vs_threshold_BBCAL = new TH1D( *hrate_vs_esum_BBCAL );
  hrate_vs_threshold_BBCAL->SetName("hrate_vs_threshold_BBCAL");
  hrate_vs_threshold_BBCAL->SetTitle("Integrated rate vs energy threshold, BBCAL");
  hrate_vs_threshold_BBCAL->SetYTitle("Hz/uA on 15-cm LD2");
  hrate_vs_threshold_BBCAL->SetXTitle("PS+SH threshold (GeV)");

  for( int i=1; i<=hrate_vs_threshold_BBCAL->GetNbinsX(); i++ ){
    double int_i, dint_i;
    int_i = hrate_vs_esum_BBCAL->IntegralAndError( i, hrate_vs_threshold_BBCAL->GetNbinsX(), dint_i );

    hrate_vs_threshold_BBCAL->SetBinContent(i, int_i );
    hrate_vs_threshold_BBCAL->SetBinError(i, dint_i );
    
  }

  TH1D *hrate_vs_threshold_HCAL = new TH1D( *hrate_vs_esum_HCAL );
  hrate_vs_threshold_HCAL->SetName("hrate_vs_threshold_HCAL");
  hrate_vs_threshold_HCAL->SetTitle("Integrated rate vs energy threshold, HCAL");
  hrate_vs_threshold_HCAL->SetYTitle("Hz/uA on 15-cm LD2");
  hrate_vs_threshold_HCAL->SetXTitle("HCAL threshold (GeV)");

  for( int i=1; i<=hrate_vs_threshold_HCAL->GetNbinsX(); i++ ){
    double int_i, dint_i;
    int_i = hrate_vs_esum_HCAL->IntegralAndError( i, hrate_vs_threshold_HCAL->GetNbinsX(), dint_i );

    hrate_vs_threshold_HCAL->SetBinContent(i, int_i );
    hrate_vs_threshold_HCAL->SetBinError(i, dint_i );
    
  }

  cout << "Total number of generated events = " << ngen << endl;

  BBGEM_edep->Write();
  BBGEM_logedep_keV->Write();
  BBGEM_yvsx->Write();
  BBGEM_dyvsdx->Write();
  
  fout->Write();
  
}
