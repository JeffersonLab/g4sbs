#include "gep_tree_with_spin.C"
#include "TChain.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include "TVector3.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TDecompSVD.h"
#include "TRotation.h"
#include "TString.h"
#include "TMath.h"
#include "TCut.h"
#include "TEventList.h"
#include "G4SBSRunData.hh"
#include "TObjArray.h"
#include "TTreeFormula.h"

const double Mp = 0.938272046; //GeV/c^2
const double mu_p = 2.792847356; //mu_N
const double kappa_p = mu_p - 1.0; //anomalous magnetic moment
const double PI = TMath::Pi();
const double SBS_thetabend = 0.0*PI/180.0; //central bend angle

//The purpose of this macro is to create a reduced ROOT tree with a flat "ntuple" structure that contains only the information about an event required for the fitting of forward and backward optical reconstruction coefficients and spin transport matrices for SBS in the "GEP" configuration, using the g4sbs ROOT trees:
void G4SBS_spin_transport( const char *inputfilename, const char *outputfilename, int pflag=0 ){
  ifstream infile(inputfilename);
  
  TFile *fout = new TFile(outputfilename,"RECREATE");
  TTree *Tout = new TTree("Tout","SBS optics and spin transport for GEP");

  double beta, gamma; //relativistic factors for the proton.
  double xfp, yfp, xpfp, ypfp, t; //include time-of-flight!
  double xfprecon,yfprecon,xpfprecon,ypfprecon; //includes detector resolution!
  double xtar, ytar, xptar, yptar, p;
  double chi, chiphi; //usual dipole precession angles
  //Euler angles for the trajectory bend:
  double phitrack, thetatrack, psitrack;
  double Pxtg, Pytg, Pztg; //In TRANSPORT coordinates!
  double Pxfp, Pyfp, Pzfp; //In TRANSPORT coordinates!
  double Pxfpgeom, Pyfpgeom, Pzfpgeom;
  double Pxfpdipole, Pyfpdipole, Pzfpdipole;
  double phispin, thetaspin, psispin;

  double Rxx,Rxy,Rxz,Ryx,Ryy,Ryz,Rzx,Rzy,Rzz;
  
  Tout->Branch("beta",&beta);
  Tout->Branch("gamma",&gamma);
  Tout->Branch("xfp",&xfp);
  Tout->Branch("yfp",&yfp);
  Tout->Branch("xpfp",&xpfp);
  Tout->Branch("ypfp",&ypfp);
  Tout->Branch("xfprecon",&xfprecon);
  Tout->Branch("yfprecon",&yfprecon);
  Tout->Branch("xpfprecon",&xpfprecon);
  Tout->Branch("ypfprecon",&ypfprecon);
  Tout->Branch("t",&t);
  Tout->Branch("xtar",&xtar);
  Tout->Branch("ytar",&ytar);
  Tout->Branch("xptar",&xptar);
  Tout->Branch("yptar",&yptar);
  Tout->Branch("p",&p);
  Tout->Branch("chi",&chi);
  Tout->Branch("chiphi",&chiphi);
  Tout->Branch("phitrack",&phitrack);
  Tout->Branch("thetatrack",&thetatrack);
  Tout->Branch("psitrack",&psitrack);
  //Tout->Branch("phispin",&phispin);
  //Tout->Branch("thetaspin",&thetaspin);
  //Tout->Branch("psispin",&psispin);
  Tout->Branch("Pxtg",&Pxtg);
  Tout->Branch("Pytg",&Pytg);
  Tout->Branch("Pztg",&Pztg);
  Tout->Branch("Pxfp",&Pxfp);
  Tout->Branch("Pyfp",&Pyfp);
  Tout->Branch("Pzfp",&Pzfp);
  Tout->Branch("Pxfpgeom",&Pxfpgeom);
  Tout->Branch("Pyfpgeom",&Pyfpgeom);
  Tout->Branch("Pzfpgeom",&Pzfpgeom);
  // Tout->Branch("Rxx",&Rxx);
  // Tout->Branch("Rxy",&Rxy);
  // Tout->Branch("Rxz",&Rxz);
  // Tout->Branch("Ryx",&Ryx);
  // Tout->Branch("Ryy",&Ryy);
  // Tout->Branch("Ryz",&Ryz);
  // Tout->Branch("Rzx",&Rzx);
  // Tout->Branch("Rzy",&Rzy);
  // Tout->Branch("Rzz",&Rzz);
  Tout->Branch("Pxfpdipole",&Pxfpdipole);
  Tout->Branch("Pyfpdipole",&Pyfpdipole);
  Tout->Branch("Pzfpdipole",&Pzfpdipole);

  double total_precession_angle; //"big" precession angle"
  double total_precession_axis_yaw; //"yaw" of precession axis wrt y axis
  double total_precession_axis_roll; //"roll" of precession axis wrt y axis 
  //double total_precession_axis_z;
  // Tout->Branch("total_precession_angle",&total_precession_angle);
  // Tout->Branch("total_precession_axis_yaw",&total_precession_axis_yaw);
  // Tout->Branch("total_precession_axis_roll",&total_precession_axis_roll);

  
  TChain *C = new TChain("T");

  map<TString,double> SBSang_file;
  //map<TString,double> SBStrkrpitch_file;

  double SBStheta_default = 16.9*PI/180.0;
  
  TFile *ftemp;
  G4SBSRunData *rd;

  cout << "Reading list of files..." << endl;
  TString currentline;
  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){

      C->Add(currentline.Data());
	     
      //      cout << currentline << endl;
      
      //     ftemp = new TFile(currentline.Data(), "READ");
      
      //ftemp->ls();
      
    }
  }

  std::set<TString> bad_file_list;

  for( int ifile=0; ifile<C->GetNtrees(); ifile++ ){
    TObjArray *list_temp = C->GetListOfFiles();
    
    TString fname_temp = (*list_temp)[0]->GetTitle();
    cout << "file " << ifile << ", name = " << fname_temp.Data() << endl; 
    
    TFile *ftemp = new TFile(fname_temp.Data(),"READ");
    
    if( !ftemp->IsZombie() ){
      G4SBSRunData *rdtemp;
      ftemp->GetObject("run_data",rdtemp);
      if( rdtemp ){
  	SBSang_file[fname_temp] = rdtemp->fSBStheta;
  	//SBStrkrpitch_file[fname_temp] = rdtemp->f
      } else {
  	bad_file_list.insert( fname_temp );
      }
      ftemp->Close();
      delete ftemp;
    } else {
      bad_file_list.insert( fname_temp );
      ftemp->Close(); 
      delete ftemp;
    }
  }
  
  cout << "finished reading list of files" 
       << endl;
  
  if( C->GetNtrees() <= 0 ){
    cout << "Did not load any files, quitting..." << endl;
    return;
  }

  TCut global_cut = "";
  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      global_cut += currentline.Data();
    }
  }

  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", global_cut, C );
  //  TEventList *elist = new TEventList("elist");
  //C->Draw(">>elist",global_cut);
  
  gep_tree_with_spin *T = new gep_tree_with_spin(C);

  long nevent=0;

  double SBStheta = SBStheta_default;
  int treenum = 0;
  int oldtreenum = -1;
  
  while( T->GetEntry(nevent++) ){
    if( nevent%10000 == 0 ) cout << nevent << endl;

    //Check for file changeover and grab SBS angle 
   
    treenum = C->GetTreeNumber();

    bool goodfile = true;
    
    if( treenum != oldtreenum ){
      GlobalCut->UpdateFormulaLeaves();
      
      oldtreenum = treenum;
    
      TString fname = C->GetFile()->GetName();

      if( SBSang_file.find(fname) != SBSang_file.end() ){
	SBStheta = SBSang_file[fname];
      } 
      if( bad_file_list.find( fname ) != bad_file_list.end() ) goodfile = false;
    }
    
    bool passed_global_cut = GlobalCut->EvalInstance(0) != 0;

    if( goodfile && passed_global_cut ){
      if( T->Harm_FT_Track_ntracks == 1 &&
	  (*(T->Harm_FT_Track_MID))[0] == 0 ){ //good track in SBS is required regardless of globalcut
	xfp = (*(T->Harm_FT_Track_X))[0];
	yfp = (*(T->Harm_FT_Track_Y))[0];
	xpfp = (*(T->Harm_FT_Track_Xp))[0];
	ypfp = (*(T->Harm_FT_Track_Yp))[0];
	
	xfprecon = (*(T->Harm_FT_Track_Xfit))[0];
	xpfprecon = (*(T->Harm_FT_Track_Xpfit))[0];
	yfprecon = (*(T->Harm_FT_Track_Yfit))[0];
	ypfprecon = (*(T->Harm_FT_Track_Ypfit))[0];
	
	t = (*(T->Harm_FT_Track_T))[0];
	//p = T->ev_ep; //use "track" momentum (because most of the eloss will occur on the way out of the target?
	p = (*(T->Harm_FT_Track_P))[0];
	beta = p/sqrt(pow(p,2)+pow(Mp,2));
	gamma = sqrt(1.0 + pow(p/Mp,2));

	//Now need target quantities!

	//ouble px,py,pz;
	//if( pflag == 0 ){ //default is to assume gun generator! use electron momentum variables: 
	double px = T->ev_epx;
	double py = T->ev_epy;
	double pz = T->ev_epz;

	if( pflag != 0 ){ //override default behavior, use nucleon momentum variables: 
	  px = T->ev_npx;
	  py = T->ev_npy;
	  pz = T->ev_npz;
	}

	//double SBStheta = T->gen_thsbs; //on beam right!

	TVector3 phall(px,py,pz);
	TVector3 phall_unit = phall.Unit();
	TVector3 SBSzaxis( -sin(SBStheta),0,cos(SBStheta) );
	TVector3 SBSxaxis(0,-1,0);
	TVector3 SBSyaxis = (SBSzaxis.Cross(SBSxaxis)).Unit();

	TVector3 punit_TRANSPORT( phall_unit.Dot( SBSxaxis ),
				  phall_unit.Dot( SBSyaxis ),
				  phall_unit.Dot( SBSzaxis ) );

	TVector3 vertex_hall( T->ev_vx, T->ev_vy, T->ev_vz ); //in m
      
	xptar = punit_TRANSPORT.X()/punit_TRANSPORT.Z();
	yptar = punit_TRANSPORT.Y()/punit_TRANSPORT.Z();
      
	TVector3 vertex_TRANSPORT( vertex_hall.Dot( SBSxaxis ),
				   vertex_hall.Dot( SBSyaxis ),
				   vertex_hall.Dot( SBSzaxis ) );

	xtar = vertex_TRANSPORT.X() - xptar * vertex_TRANSPORT.Z();
	ytar = vertex_TRANSPORT.Y() - yptar * vertex_TRANSPORT.Z();

	double thetabend = SBS_thetabend + atan(xptar) - atan(xpfp);
	double phibend = atan(ypfp) - atan(yptar);

	chi = gamma * kappa_p * thetabend;
	chiphi = gamma * kappa_p * phibend;

	//Compute Euler angles: for this we need to compute the axes of the tgt and fp comoving coordinate systems:
	TVector3 pfp_TRANSPORT( xpfp, ypfp, 1.0 );
	pfp_TRANSPORT = pfp_TRANSPORT.Unit();

	TVector3 SBSzaxis_fp( -sin(SBS_thetabend), 0, cos(SBS_thetabend) );
	TVector3 SBSxaxis_fp( cos(SBS_thetabend), 0, sin(SBS_thetabend) );
	TVector3 SBSyaxis_fp(0,1,0);

	TVector3 pfp_bent =
	  pfp_TRANSPORT.X() * SBSxaxis_fp +
	  pfp_TRANSPORT.Y() * SBSyaxis_fp +
	  pfp_TRANSPORT.Z() * SBSzaxis_fp;

	TVector3 xaxis_comoving_fp = SBSxaxis_fp;
	TVector3 zaxis_comoving_fp = pfp_bent;
	TVector3 yaxis_comoving_fp = (zaxis_comoving_fp.Cross( xaxis_comoving_fp )).Unit();
	xaxis_comoving_fp = (yaxis_comoving_fp.Cross(zaxis_comoving_fp)).Unit();

	TVector3 xaxis_comoving_tgt(1,0,0);
	TVector3 zaxis_comoving_tgt = punit_TRANSPORT;
	TVector3 yaxis_comoving_tgt = (zaxis_comoving_tgt.Cross(xaxis_comoving_tgt)).Unit();
	xaxis_comoving_tgt = (yaxis_comoving_tgt.Cross(zaxis_comoving_tgt)).Unit();

	// phitrack = atan2( yaxis_comoving_fp.Dot( xaxis_comoving_tgt ), -zaxis_comoving_fp.Dot( xaxis_comoving_tgt ) );
	// psitrack = atan2( xaxis_comoving_fp.Dot( yaxis_comoving_tgt ), xaxis_comoving_fp.Dot( zaxis_comoving_tgt ) );
	// thetatrack = acos( xaxis_comoving_fp.Dot( xaxis_comoving_tgt ) );

	TVector3 rotation_axis_track = (punit_TRANSPORT.Cross( pfp_bent) ).Unit();
	double rotation_angle_track = acos( punit_TRANSPORT.Dot( pfp_bent ) );

	TRotation Rtrack;
	Rtrack.Rotate( rotation_angle_track, rotation_axis_track );

	phitrack = atan2( Rtrack.XY(), -Rtrack.XZ() );
	thetatrack = acos( Rtrack.XX() );
	psitrack = atan2( Rtrack.YX(), Rtrack.ZX() );
      
	//In TRANSPORT coordinates at the target:
	Pxtg = T->ev_Sx;
	Pytg = T->ev_Sy;
	Pztg = T->ev_Sz;
      
	//In TRANSPORT coordinates at the fp:
	Pxfp = (*(T->Harm_FT_Track_Sx))[0];
	Pyfp = (*(T->Harm_FT_Track_Sy))[0];
	Pzfp = (*(T->Harm_FT_Track_Sz))[0];

	TVector3 SpinTg(Pxtg,Pytg,Pztg);
	SpinTg = SpinTg.Unit();

	TVector3 SpinFp(Pxfp,Pyfp,Pzfp);
	SpinFp = SpinFp.Unit();


	TVector3 SpinFp_bent =
	  SpinFp.X() * SBSxaxis_fp +
	  SpinFp.Y() * SBSyaxis_fp +
	  SpinFp.Z() * SBSzaxis_fp;
      
	//Compute geometric approximation rotation matrix:
	TVector3 SpinTg_comoving( SpinTg.Dot( xaxis_comoving_tgt ),
				  SpinTg.Dot( yaxis_comoving_tgt ),
				  SpinTg.Dot( zaxis_comoving_tgt ) );

	TRotation R;
	R.RotateX( -gamma*kappa_p*(phitrack+psitrack)/2.0 );
	R.RotateY( -gamma*kappa_p*thetatrack );
	R.RotateX( -gamma*kappa_p*(phitrack+psitrack)/2.0 );

	TRotation Rdipole;
	Rdipole.RotateY( -gamma*kappa_p*thetatrack );
      
	TVector3 SpinFp_comoving_geom = R * SpinTg_comoving;
	TVector3 SpinFp_comoving_dipole = Rdipole * SpinTg_comoving;

	//Now this is expressed in "bent" FP coordinates or target transport coordinates:
	TVector3 SpinFp_bent_geom =
	  SpinFp_comoving_geom.X() * xaxis_comoving_fp +
	  SpinFp_comoving_geom.Y() * yaxis_comoving_fp +
	  SpinFp_comoving_geom.Z() * zaxis_comoving_fp; 

	//TVector3 SpinFp_TRANSPORT_geom( SpinFp_bent_geom.Dot( SBS

	TVector3 SpinFp_TRANSPORT_geom( SpinFp_bent_geom.Dot( SBSxaxis_fp ),
					SpinFp_bent_geom.Dot( SBSyaxis_fp ),
					SpinFp_bent_geom.Dot( SBSzaxis_fp ) );

	//Now this is expressed in "bent" FP coordinates or target transport coordinates:
	TVector3 SpinFp_bent_dipole =
	  SpinFp_comoving_dipole.X() * xaxis_comoving_fp +
	  SpinFp_comoving_dipole.Y() * yaxis_comoving_fp +
	  SpinFp_comoving_dipole.Z() * zaxis_comoving_fp; 

	//TVector3 SpinFp_TRANSPORT_dipole( SpinFp_bent_dipole.Dot( SBS

	TVector3 SpinFp_TRANSPORT_dipole( SpinFp_bent_dipole.Dot( SBSxaxis_fp ),
					  SpinFp_bent_dipole.Dot( SBSyaxis_fp ),
					  SpinFp_bent_dipole.Dot( SBSzaxis_fp ) );
      
	// cout << "(Pxfp, Pyfp, Pzfp)=(" << Pxfp << ", " << Pyfp << ", " << Pzfp << "), geometric approx = (" << SpinFp_TRANSPORT_geom.X()
	// 	   << ", " << SpinFp_TRANSPORT_geom.Y() << ", " << SpinFp_TRANSPORT_geom.Z() << ")" << endl;

	Pxfpgeom = SpinFp_TRANSPORT_geom.X();
	Pyfpgeom = SpinFp_TRANSPORT_geom.Y();
	Pzfpgeom = SpinFp_TRANSPORT_geom.Z();

	Pxfpdipole = SpinFp_TRANSPORT_dipole.X();
	Pyfpdipole = SpinFp_TRANSPORT_dipole.Y();
	Pzfpdipole = SpinFp_TRANSPORT_dipole.Z();

	//Compute the total rotation of the spin in fixed transport coordinates:
	TVector3 spin_rotation_axis = (SpinTg.Cross(SpinFp)).Unit();
	double spin_rotation_angle = acos( SpinTg.Dot(SpinFp));

	//cout << "Spin rotation axis in TRANSPORT coordinates:" << endl;
	//spin_rotation_axis.Print();
	//cout << "Total spin rotation angle in TRANSPORT coordinates = " << spin_rotation_angle*180.0/TMath::Pi() << " deg" << endl;
      
      
	TRotation Rspin;
	Rspin.Rotate( spin_rotation_angle, spin_rotation_axis );

	phispin = atan2( Rspin.XY(), -Rspin.XZ() );
	thetaspin = acos( Rspin.XX() );
	psispin = atan2( Rspin.YX(), Rspin.ZX() );

	Rxx = Rspin.XX();
	Rxy = Rspin.XY();
	Rxz = Rspin.XZ();
	Ryx = Rspin.YX();
	Ryy = Rspin.YY();
	Ryz = Rspin.YZ();
	Rzx = Rspin.ZX();
	Rzy = Rspin.ZY();
	Rzz = Rspin.ZZ();

	//Now compute full rotation a different way:
	TRotation RdipoleT;
	RdipoleT.RotateY( -(chi + thetabend - SBS_thetabend) );

	TRotation RdipoleTinv = RdipoleT.Inverse();
	//If Rtotal = Rz Rx Rd
	// And we know Rtotal and Rd,
	// Then Rz Rx = Rtotal * Rdinv
      
	TRotation RzRx = Rspin * RdipoleTinv; //

	thetaspin = chi + thetabend - SBS_thetabend;
	phispin = asin( RzRx.ZY() );
	//psispin = atan2( RzRx.YX(), RzRx.XX() ); 
	psispin = atan2( -RzRx.XY(), RzRx.YY() );

	// cout << "Rspin: " << endl;
	// cout << "| xx, xy, xz | = |" << Rspin.XX() << ", " << Rspin.XY() << ", " << Rspin.XZ() << "|" << endl
	// 	   << "| yx, yy, yz | = |" << Rspin.YX() << ", " << Rspin.YY() << ", " << Rspin.YZ() << "|" << endl
	// 	   << "| zx, zy, zz | = |" << Rspin.ZX() << ", " << Rspin.ZY() << ", " << Rspin.ZZ() << "|" << endl;

	// cout << "Rdipole: " << endl;
	// cout << "| xx, xy, xz | = |" << Rdipole.XX() << ", " << Rdipole.XY() << ", " << Rdipole.XZ() << "|" << endl
	// 	   << "| yx, yy, yz | = |" << Rdipole.YX() << ", " << Rdipole.YY() << ", " << Rdipole.YZ() << "|" << endl
	// 	   << "| zx, zy, zz | = |" << Rdipole.ZX() << ", " << Rdipole.ZY() << ", " << Rdipole.ZZ() << "|" << endl;
      
	// cout << "RzRx: " << endl;
	// cout << "| xx, xy, xz | = |" << RzRx.XX() << ", " << RzRx.XY() << ", " << RzRx.XZ() << "|" << endl
	// 	   << "| yx, yy, yz | = |" << RzRx.YX() << ", " << RzRx.YY() << ", " << RzRx.YZ() << "|" << endl
	// 	   << "| zx, zy, zz | = |" << RzRx.ZX() << ", " << RzRx.ZY() << ", " << RzRx.ZZ() << "|" << endl;
      
      
      
	//Since "RzRx" should be close to the identity, the XX and ZZ matrix elements should be "small".
      
	Tout->Fill();
      }
    }
  }

  //  elist->Delete();
  
  fout->Write();
  
}
