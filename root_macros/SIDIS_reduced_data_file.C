#include "g4sbs_tree.C"
#include <iostream>
#include <fstream>
#include <vector>
#include "TSystem.h"
#include "TString.h"
#include "TEventList.h"
#include "TCut.h"
#include "TObject.h"
#include "TROOT.h"
#include "G4SBSRunData.hh"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TRandom.h"

const double mpi  = 0.13957018; 
const double mpi0 = 0.1349766;
const double mK   = 0.493677;
const double mp   = 0.938272046;
const double mn   = 0.939565379;

//Declare tree variables at global scope:
//Kinematics:
double ThetaBB, ThetaSBS, Ebeam;

//"True" kinematic variables:
double vx, vy, vz;
double ep, eth, eph, Eh, hp, hth, hph;
double xbj,y,z,Q2,pT,phih,W2,MX2;
int nucleon;
int hadron;
double Weight;

//Track variables:
double exfp, eyfp, expfp, eypfp;
double hxfp, hyfp, hxpfp, hypfp;

//"Reconstructed" target quantities:
double extar, eytar, exptar, eyptar, evzrecon, eprecon, ethrecon, ephrecon;
double hxtar, hytar, hxptar, hyptar, hvzrecon, hprecon, hthrecon, hphrecon;

//"Reconstructed" SIDIS kinematic quantities:
double xrecon, yrecon, zrecon, Q2recon, pTrecon, phihrecon,W2recon,MX2recon;
double mpi0recon; //only applicable for pi0 case.

TFile *fout;
TTree *Tout;

void init_tree(){
  Tout = new TTree("Tout", "SIDIS reduced data file");
  
  Tout->Branch("Weight", &Weight, "Weight/D");
  Tout->Branch("Hadron", &hadron, "Hadron/I");
  Tout->Branch("Nucleon", &nucleon, "Nucleon/I");
  Tout->Branch("ThetaBB", &ThetaBB, "ThetaBB/D");
  Tout->Branch("ThetaSBS", &ThetaSBS, "ThetaSBS/D");
  Tout->Branch("Ebeam", &Ebeam, "Ebeam/D");
  Tout->Branch("vx", &vx, "vx/D");
  Tout->Branch("vy", &vy, "vy/D");
  Tout->Branch("vz", &vz, "vz/D");
  Tout->Branch("ep", &ep, "ep/D");
  Tout->Branch("eth", &eth, "eth/D");
  Tout->Branch("eph", &eph, "eph/D");
  Tout->Branch("Eh", &Eh, "Eh/D");
  Tout->Branch("hp", &hp, "hp/D");
  Tout->Branch("hth", &hth, "hth/D");
  Tout->Branch("hph", &hph, "hph/D");
  Tout->Branch("x", &xbj, "x/D");
  Tout->Branch("y", &y, "y/D");
  Tout->Branch("z", &z, "z/D");
  Tout->Branch("Q2", &Q2, "Q2/D");
  Tout->Branch("pT", &pT, "pT/D");
  Tout->Branch("phih", &phih, "phih/D");
  Tout->Branch("W2", &W2, "W2/D");
  Tout->Branch("MX2", &MX2, "MX2/D");
  //focal plane track variables:
  Tout->Branch("exfp", &exfp, "exfp/D");
  Tout->Branch("eyfp", &eyfp, "eyfp/D");
  Tout->Branch("expfp", &exfp, "expfp/D");
  Tout->Branch("eypfp", &eyfp, "eypfp/D");
  Tout->Branch("hxfp", &hxfp, "hxfp/D");
  Tout->Branch("hyfp", &hyfp, "hyfp/D");
  Tout->Branch("hxpfp", &hxfp, "hxpfp/D");
  Tout->Branch("hypfp", &hyfp, "hypfp/D");
  //reconstructed particle kinematics:
  Tout->Branch("extar",&extar,"extar/D");
  Tout->Branch("eytar",&eytar,"eytar/D");
  Tout->Branch("exptar",&exptar,"exptar/D");
  Tout->Branch("eyptar",&eyptar,"eyptar/D");
  Tout->Branch("evzrecon",&evzrecon,"evzrecon/D");
  Tout->Branch("eprecon",&eprecon,"eprecon/D");
  Tout->Branch("ethrecon",&ethrecon,"ethrecon/D");
  Tout->Branch("ephrecon",&ephrecon,"ephrecon/D");
  
  Tout->Branch("hxtar",&hxtar,"hxtar/D");
  Tout->Branch("hytar",&hytar,"hytar/D");
  Tout->Branch("hxptar",&hxptar,"hxptar/D");
  Tout->Branch("hyptar",&hyptar,"hyptar/D");
  Tout->Branch("hvzrecon",&hvzrecon,"hvzrecon/D");
  Tout->Branch("hprecon",&hprecon,"hprecon/D");
  Tout->Branch("hthrecon",&hthrecon,"hthrecon/D");
  Tout->Branch("hphrecon",&hphrecon,"hphrecon/D");
  
  //Reconstructed SIDIS quantities:
  Tout->Branch("xrecon", &xrecon, "xrecon/D");
  Tout->Branch("yrecon", &yrecon, "yrecon/D");
  Tout->Branch("zrecon", &zrecon, "zrecon/D");
  Tout->Branch("Q2recon", &Q2recon, "Q2recon/D");
  Tout->Branch("pTrecon", &pTrecon, "pTrecon/D");
  Tout->Branch("phihrecon", &phihrecon, "phihrecon/D");
  Tout->Branch("W2recon", &W2recon, "W2recon/D");
  Tout->Branch("MX2recon", &MX2recon, "MX2recon/D");
  Tout->Branch("mpi0recon", &mpi0recon, "mpi0recon/D");
}

void change_file(const char *newfilename){
  fout->cd();
  Tout->Write();
  fout->Close();
  
  //fout->Open(newfilename,"RECREATE");
  fout = new TFile(newfilename,"RECREATE");
  init_tree();
}

void SIDIS_reduced_data_file( const char *setupfilename, const char *outputfilename ){
  
  double pixelsize = 0.15;
  double Hcal_height = 3.3;
  double Hcal_width = 1.65;

  TRandom3 num(0);

  double Mh[7] = {mp, mK, mpi, mpi0, mpi, mK, mp};

  double PI = TMath::Pi();
  //  double Mh = mpi; //default to mpi;

  ifstream setupfile(setupfilename); 

  TString currentline;

  vector<int> ngen;
  vector<int> ntries;
  //vector<TString> list_of_files;
  vector<double> file_weights;

  TChain *C = new TChain("T");

  int nfiles = 0;

  //Add files to the tree:

  G4SBSRunData *rtemp;

  while( currentline.ReadLine(setupfile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      TFile Ftemp(currentline.Data(),"READ");
      if( !Ftemp.IsZombie() ){
	Ftemp.GetObject("run_data",rtemp);
	if( rtemp ){ 
	  ngen.push_back( rtemp->fNthrown );
	  ntries.push_back( rtemp->fNtries );
	  file_weights.push_back( double( rtemp->fNthrown )/double( rtemp->fNtries ) );
	  cout << "weight = " << file_weights[nfiles] << endl;
	  C->Add(currentline.Data());
	  nfiles++;
	}
      }
    }
  }

  //Unlike previous case, we do not apply a TCut and perform a TTree::Draw() command at this stage

  g4sbs_tree *T = new g4sbs_tree(C);

  double chargesep_m = 0.25; //slope of x' vs xfp x' = b + m * xfp //Bend angle is the same at a given field strength; therefore the  
  double chargesep_b = 0.0;
  
  double ndays;
  //setupfile >> ndays;
  setupfile >> chargesep_m >> chargesep_b;

  TString fname_sbs_optics_upbend, fname_sbs_optics_downbend, fname_bb_optics;
  fname_sbs_optics_upbend.ReadLine( setupfile );
  fname_sbs_optics_downbend.ReadLine( setupfile );
  fname_bb_optics.ReadLine( setupfile );

  ifstream file_sbs_optics_upbend(fname_sbs_optics_upbend.Data());
  ifstream file_sbs_optics_downbend(fname_sbs_optics_downbend.Data());
  ifstream file_bb_optics(fname_bb_optics.Data());

  double nterms_SBS[2];
  vector<double> SBScoeff[2][4];
  vector<double> SBSexpon[2][5];

  //int nterms_SBS_up = 0;
  file_sbs_optics_upbend >> nterms_SBS[0];
  for(int term=0; term<nterms_SBS[0]; term++){
    for(int icoeff=0; icoeff<4; icoeff++){
      double Ctemp;
      file_sbs_optics_upbend >> Ctemp;
      SBScoeff[0][icoeff].push_back( Ctemp );
    }
    for(int iexpon=0; iexpon<5; iexpon++){
      int itemp;
      file_sbs_optics_upbend >> itemp;
      SBSexpon[0][iexpon].push_back( itemp );
    }
  }

  //  int nterms_SBS_ = 0;
  file_sbs_optics_downbend >> nterms_SBS[1];
  for(int term=0; term<nterms_SBS[1]; term++){
    for(int icoeff=0; icoeff<4; icoeff++){
      double Ctemp;
      file_sbs_optics_downbend >> Ctemp;
      SBScoeff[1][icoeff].push_back( Ctemp );
    }
    for(int iexpon=0; iexpon<5; iexpon++){
      int itemp;
      file_sbs_optics_downbend >> itemp;
      SBSexpon[1][iexpon].push_back(itemp);
    }
  }
  

  vector<double> BBcoeff[4];
  vector<double> BBexpon[5];

  int nterms_BB=0;
  file_bb_optics >> nterms_BB;
  for(int term=0; term<nterms_BB; term++){
    for(int icoeff=0; icoeff<4; icoeff++){
      double Ctemp;
      file_bb_optics >> Ctemp;
      BBcoeff[icoeff].push_back(Ctemp);
    }
    for(int iexpon=0; iexpon<5; iexpon++){
      int itemp;
      file_bb_optics >> itemp;
      BBexpon[iexpon].push_back( itemp );
    }
  }


  for(int i=0; i<nfiles; i++){
    file_weights[i] /= double(nfiles); //Now each event weight only represents a rate.
    //file_weights[i] *= ndays * 24.0 * 3600.0;
  }

  int max_events_per_file=1000000;

  int ifile=0;
  TString prefix = outputfilename;
  prefix.ReplaceAll(".root","");
  TString outfilename;
  outfilename.Form("%s_%d.root",prefix.Data(),ifile);

  fout = new TFile( outfilename, "RECREATE" );
  init_tree();
  
  long nevent=0;
  int ngood = 0;
  while( T->GetEntry( nevent++ ) ){
    //cout << "loaded New event, event = " << nevent << endl;
    if( nevent%10000 == 0 ) cout << "nevent = " << nevent << ", tree number = " << T->fChain->GetTreeNumber() << endl;
    
    Ebeam = T->gen_Ebeam; //GeV
    ThetaBB = T->gen_thbb; //rad
    ThetaSBS = T->gen_thhcal; //rad

    TVector3 BBxaxis(0,-1,0);
    TVector3 BBzaxis(-sin(ThetaBB),0,cos(ThetaBB));
    TVector3 BByaxis = (BBzaxis.Cross(BBxaxis)).Unit();

    TVector3 SBSxaxis(0,-1,0);
    TVector3 SBSzaxis(sin(ThetaSBS),0,cos(ThetaSBS));
    TVector3 SBSyaxis = (SBSzaxis.Cross(SBSxaxis)).Unit();

    nucleon = T->ev_nucl; //proton or neutron?
    hadron = T->ev_hadr; //pion, kaon, charge?

    Weight = T->ev_rate * file_weights[T->fChain->GetTreeNumber()];
    
    vx = T->ev_vx;
    vy = T->ev_vy;
    vz = T->ev_vz;
    ep = T->ev_ep;
    eth = T->ev_th;
    eph = T->ev_ph;

    hp = T->ev_np;
    Eh = sqrt(pow( Mh[hadron+3],2) + pow(hp,2));
    hth = T->ev_nth;
    hph = T->ev_nph;
    
    xbj = T->ev_xbj; // = Q2/2P dot q
    y = 1.0 - ep/Ebeam;
    z = T->ev_z;
    Q2 = T->ev_Q2;
    W2 = T->ev_W2;
    MX2 = T->ev_MX2;
    
    pT = T->ev_phperp;
    phih = T->ev_phih;

    //Reconstruct electron side, regardless whether this is a pi0 event:
    bool goodelectrontrack = false;
    bool goodhadrontrack = false;
    for(int track=0; track<T->ntracks; track++){
      if( (*(T->trackerid))[track] == 0 && (*(T->trackid))[track] == 1 ){ //Primary electron track!
	exfp = (*(T->trackxfit))[track];
	eyfp = (*(T->trackyfit))[track];
	expfp = (*(T->trackxpfit))[track];
	eypfp = (*(T->trackypfit))[track];
	goodelectrontrack = true;
      }
    
      if( (*(T->trackerid))[track] == 1 && (*(T->trackid))[track] == 2 ){ //Primary hadron track!
	hxfp = (*(T->trackxfit))[track];
	hyfp = (*(T->trackyfit))[track];
	hxpfp = (*(T->trackxpfit))[track];
	hypfp = (*(T->trackypfit))[track];
	goodhadrontrack = true;
      }
    }
    
    
    if( goodelectrontrack ){
     
      //Reconstruct electron side:
      extar = 0.0;
      for(int iter=0; iter<2; iter++){
	double recon_sum[4] = {0.0,0.0,0.0,0.0};
	for(int coeff=0; coeff<4; coeff++){
	  for(int term=0; term<nterms_BB; term++){
	    recon_sum[coeff] += BBcoeff[coeff][term] * pow(exfp,BBexpon[0][term]) * pow(eyfp,BBexpon[1][term]) * pow(expfp,BBexpon[2][term]) * pow(eypfp,BBexpon[3][term]) * pow(extar,BBexpon[4][term]);
	  }
	}
	exptar = recon_sum[0]; 
	eyptar = recon_sum[1];
	eytar = recon_sum[2];
	eprecon = 1.0/recon_sum[3];
	
	evzrecon = eytar / ( sin(ThetaBB) - cos(ThetaBB)*eyptar );
	extar = -vy - evzrecon * cos(ThetaBB) * exptar;
      }

      

      TVector3 ephat_spec( exptar, eyptar, 1.0 );
      ephat_spec = ephat_spec.Unit();
      
      TVector3 ephat = ephat_spec.X() * BBxaxis + ephat_spec.Y() * BByaxis + ephat_spec.Z() * BBzaxis;
      TLorentzVector kPrime( eprecon * ephat, eprecon );
      TLorentzVector kBeam( 0.0, 0.0, Ebeam, Ebeam );
      TLorentzVector qVect = kBeam - kPrime;
      TLorentzVector PVect(0,0,0,mp);
      TLorentzVector WVect = PVect + qVect;

      ethrecon = kPrime.Theta();
      ephrecon = kPrime.Phi();

      Q2recon = -qVect.M2();
      xrecon = Q2recon / (2.0*mp*(Ebeam-eprecon));
      W2recon = WVect.M2();
      yrecon = 1.0-eprecon/Ebeam;

      TVector3 lepton_plane_yaxis = (kBeam.Vect().Cross(kPrime.Vect()) ).Unit();
      TVector3 lepton_plane_xaxis = (lepton_plane_yaxis.Cross( qVect.Vect() ) ).Unit();
      TVector3 lepton_plane_zaxis = qVect.Vect().Unit();

      if( hadron == 0 ){ //pi0 case: take two highest-energy hits with mid == 2 && pid == 22
	int hit1 = -1, hit2 = -1;
	double ehitmax = 0.0, ehit1=0.0, ehit2=0.0;
	double xhit1, yhit1, xhit2, yhit2;
	int ix1, iy1, ix2, iy2;
	for( int ihit=0; ihit<T->hc_ndata; ihit++){
	  for( int jhit=ihit+1; jhit<T->hc_ndata; jhit++){
	    if( T->hc_pid[ihit] == 22 && T->hc_pid[jhit] == 22 && T->hc_mid[ihit] == 2 && T->hc_mid[jhit] == 2 ){ //Two photons from pi0 decay:
	      hit1 = ihit; 
	      hit2 = jhit; 
	      ehit1 = T->hc_e[ihit];
	      ehit2 = T->hc_e[jhit];
	      
	      xhit1 = T->hc_x[ihit];
	      yhit1 = T->hc_y[ihit];
	      
	      xhit2 = T->hc_x[jhit];
	      yhit2 = T->hc_y[jhit];

	      ix1 = int( (xhit1 + 0.5*Hcal_width)/pixelsize );
	      iy1 = int( (yhit1 + 0.5*Hcal_height)/pixelsize );
	      ix2 = int( (xhit2 + 0.5*Hcal_width)/pixelsize );
	      iy2 = int( (yhit2 + 0.5*Hcal_height)/pixelsize );

	      double x1recon = (ix1+0.5)*pixelsize - 0.5*Hcal_width;
	      double y1recon = (iy1+0.5)*pixelsize - 0.5*Hcal_height;
	      double x2recon = (ix2+0.5)*pixelsize - 0.5*Hcal_width;
	      double y2recon = (iy2+0.5)*pixelsize - 0.5*Hcal_height;
	      
	      TVector3 ray1recon( x1recon*cos(ThetaSBS)+T->gen_dhcal*sin(ThetaSBS), y1recon, -x1recon*sin(ThetaSBS)+T->gen_dhcal*cos(ThetaSBS) - evzrecon );
	      TVector3 ray2recon( x2recon*cos(ThetaSBS)+T->gen_dhcal*sin(ThetaSBS), y2recon, -x2recon*sin(ThetaSBS)+T->gen_dhcal*cos(ThetaSBS) - evzrecon );

	      double E1recon = num.Gaus( ehit1, 0.14*sqrt(ehit1) );
	      double E2recon = num.Gaus( ehit2, 0.14*sqrt(ehit2) );

	      if( (abs(ix1-ix2) > 1 || abs(iy1-iy2) > 1) && (E1recon >= 2.0 || E2recon >= 2.0) ){
		TLorentzVector Photon1( E1recon*ray1recon.Unit(), E1recon );
		TLorentzVector Photon2( E2recon*ray2recon.Unit(), E2recon );

		TLorentzVector TwoPhotons = Photon1 + Photon2;

		goodhadrontrack = true;
		mpi0recon = TwoPhotons.M();
		hprecon = TwoPhotons.P();
		hthrecon = TwoPhotons.Theta();
		hphrecon = TwoPhotons.Phi();
		
		TLorentzVector XVect = WVect - TwoPhotons;
		MX2recon = XVect.M2();

		zrecon = TwoPhotons.Dot(PVect)/qVect.Dot(PVect);

		TVector3 pTVect = TwoPhotons.Vect() - TwoPhotons.Vect().Dot(qVect.Vect())  / qVect.Vect().Mag2() * qVect.Vect();
		pTrecon = pTVect.Mag();
		phihrecon = atan2( pTVect.Dot( lepton_plane_yaxis ), pTVect.Dot( lepton_plane_xaxis ) );
	      } 
	    }
	  }
	}
      } else if( goodhadrontrack ){ // Charged hadrons: look for track and reconstruct:
	hxtar = 0.0;
	int polarity = 0; //upbending --> x' < x
	if( hxpfp > chargesep_b + chargesep_m * hxfp ){ //downbending:
	  polarity = 1;
	}


	for(int iter=0; iter<2; iter++){
	  double recon_sum[4] = {0,0,0,0};
	  for(int coeff=0; coeff<4; coeff++){
	    for(int term=0; term<nterms_SBS[polarity]; term++){
	      recon_sum[coeff] += SBScoeff[polarity][coeff][term] * 
		pow(hxfp,SBSexpon[polarity][0][term]) * 
		pow(hyfp,SBSexpon[polarity][1][term]) * 
		pow(hxpfp,SBSexpon[polarity][2][term]) * 
		pow(hypfp,SBSexpon[polarity][3][term]) * 
		pow(hxtar,SBSexpon[polarity][4][term]);
	    }
	  }
	  hxptar = recon_sum[0];
	  hyptar = recon_sum[1];
	  hytar = recon_sum[2];
	  hprecon = 1.0/recon_sum[3];
	  
	  hvzrecon = -hytar / ( sin(ThetaSBS) + cos(ThetaSBS)*hyptar );
	  hxtar = -vy - hvzrecon * cos(ThetaSBS) * hxptar;
	}
       

	TVector3 hphat_spec(hxptar,hyptar,1.0);
	hphat_spec = hphat_spec.Unit();

	TVector3 hphat = hphat_spec.X() * SBSxaxis + hphat_spec.Y() * SBSyaxis + hphat_spec.Z() * SBSzaxis;
	TLorentzVector Ph( hprecon * hphat, sqrt(pow(hprecon,2) + pow(Mh[hadron+3],2) ) );

	hthrecon = Ph.Theta();
	hphrecon = Ph.Phi();

	zrecon = Ph.Dot(PVect)/qVect.Dot(PVect);

	TLorentzVector XVect = WVect - Ph;
	MX2recon = XVect.M2();

	TVector3 pTVect = Ph.Vect() - Ph.Vect().Dot(qVect.Vect()) / qVect.Vect().Mag2() * qVect.Vect();
	pTrecon = pTVect.Mag();
	phihrecon = atan2( pTVect.Dot( lepton_plane_yaxis ), pTVect.Dot( lepton_plane_xaxis ) );


      }
    }
    

    if( goodelectrontrack && goodhadrontrack ){
      Tout->Fill();
      ngood++;
      if( ngood == max_events_per_file ){
	ifile++;
	outfilename.Form("%s_%d.root",prefix.Data(),ifile);
	change_file( outfilename.Data() );
      }
    }
  }

  Tout->Write();
  fout->Close();


}
