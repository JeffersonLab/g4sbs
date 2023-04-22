#include "G4SBSRunData.hh"
#include "gep_tree_with_spin.C"
//#include "gep_pythia6_tree.C"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TObjArray.h"
#include "TObjString.h"
#include <iostream>
#include <fstream>
#include <vector>
#include "TString.h"
#include <set>
#include <map>
#include <cstdio>
#include <string>
#include <sstream>
#include "TEventList.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TCut.h"

const double Mp = 0.938272; //GeV
const double Lx_scint_CDET = 0.51; //m
const double Ly_scint_CDET = 0.005; //5 mm
const double mu_p = 2.793; // nuclear magneton
const double PI = TMath::Pi();
const double SBS_tracker_pitch=5.0*PI/180.0; //5 degrees

const double ECAL_phe_per_GeV=300.0;
const double ECAL_max_cell_size = 0.042; //meters
const double X0_ECAL = 0.0274; //radiation length
const double Ec_ECAL = 0.015; //GeV
const double yoff_ECAL = 0.0077; //meters, average y deflection in SBS fringe field, to be SUBTRACTED from recconstructed shower coordinate!
const double sigx_ECAL = 0.008; //meters, shower x coordinate resolution
const double sigy_ECAL = 0.006; //meters, shower y coordinate resolution
const double Ltgt = 0.4; //meters
const double sigy_CDET = 0.003;

// Bin width for vertex z "filtering": 
//const double vz_bin_width = 3.0*0.0064;
//const double vz_scan_stepsize = 


const double dE_E_MPV = 0.13; //Most probable electron energy loss before reaching ECAL:

TRandom3 *random_generator;

//Parameters for optics reconstruction (forward and backward)

int nterms_optics;
vector<double> Cxptar,Cyptar,Cytar,Cpthetabend;
vector<vector<int> > Cexpon;

int nterms_foptics;
vector<double> Cxfp, Cyfp, Cxpfp, Cypfp;
vector<vector<int> > Cfexpon;


  
//Parameters for shower coordinate reconstruction:
struct shower_profile_t {
  int nbins_mom;
  double xmin;
  double xmax;
  vector<double> frac;
};

shower_profile_t profxdefault;
shower_profile_t profydefault;

int profx_nbins;
int profy_nbins;
double profx_xmin,profx_xmax;
double profy_ymin,profy_ymax;
vector<shower_profile_t> profx;
vector<shower_profile_t> profy;

void calc_shower_coordinates( double xmom, double ymom, double xmax, double ymax, double Eclust, double Rcal, double &xclust, double &yclust, double &xf, double &yf ){

  //calculate longitudinal depth of max. shower energy deposition
  double tmax = TMath::Max(0.0, X0_ECAL * (log( Eclust/Ec_ECAL ) - 0.5) );

  int binx = int( (100.0*xmax-profx_xmin)/(profx_xmax-profx_xmin)*profx_nbins );

  shower_profile_t proftemp;
  if( binx >= 0 && binx < profx_nbins ){
    proftemp = profx[binx];
  } else {
    proftemp = profxdefault;
  }
  
  int binxmom = int( (xmom - proftemp.xmin)/(proftemp.xmax-proftemp.xmin) * proftemp.nbins_mom );

  //set some sensible default value:
  xclust = xmax + xmom*ECAL_max_cell_size; 
  
  //do linear interpolation within the bin:
  if( binxmom < 0 ){
    xclust = xmax - 0.5 * ECAL_max_cell_size;
    //xf = -0.5*ECAL_max_cell_size;
  } else if( binxmom < proftemp.nbins_mom ){
    double binwidth = (proftemp.xmax-proftemp.xmin)/double(proftemp.nbins_mom);
    double flo = proftemp.frac[binxmom];
    double fhi = proftemp.frac[binxmom+1];
    
    double xlo = proftemp.xmin + binxmom * binwidth;
    //double xhi = xlo + binwidth;
    //linear interpolation within the bin:
    double f = flo + (fhi-flo) * (xmom-xlo)/binwidth;
    xclust = xmax + (f - 0.5)*ECAL_max_cell_size;

    // xf = xclust-xmax;
    //cout << "xmom, binxmom, f = " << xmom << ", " << binxmom << ", " << f << endl;
  } else {
    xclust = xmax + 0.5 * ECAL_max_cell_size;
    
  }

  //before incident-angle correction, record x position within cell:
  xf = xclust-xmax;
  
  int biny = int( (100.0*ymax-profy_ymin)/(profy_ymax-profy_ymin)*profy_nbins );

  if( biny >= 0 && biny < profy_nbins ){
    proftemp = profy[biny];
  } else {
    proftemp = profydefault;
  }

  int binymom = int( (ymom - proftemp.xmin)/(proftemp.xmax-proftemp.xmin) * proftemp.nbins_mom );

  yclust = ymax + ymom*ECAL_max_cell_size;

  if( binymom < 0 ){
    yclust = ymax - 0.5 * ECAL_max_cell_size;
  } else if( binymom < proftemp.nbins_mom ){
    double binwidth = (proftemp.xmax-proftemp.xmin)/double(proftemp.nbins_mom);
    double flo = proftemp.frac[binymom];
    double fhi = proftemp.frac[binymom+1];
    
    double ylo = proftemp.xmin + binymom * binwidth;
    //double xhi = xlo + binwidth;
    //linear interpolation within the bin:
    double f = flo + (fhi-flo) * (ymom-ylo)/binwidth;
    yclust = ymax + (f - 0.5)*ECAL_max_cell_size;

    //cout << "ymom, biny, binymom, f = " << ymom << ", " << biny << "," <<  binymom << ", " << f << endl;
  } else {
    yclust = ymax + 0.5 * ECAL_max_cell_size;
  }

  yf = yclust-ymax;
  
  //Apply incident-angle correction under the assumption that track starts at (x,y,z) = (0,0,0)
  double xptemp = xclust/Rcal;
  double yptemp = yclust/Rcal;

  double dz = tmax / sqrt(1.0 + pow(xptemp,2) + pow(yptemp,2) );

  xclust -= dz * xptemp;
  yclust -= dz * yptemp;
  
  return;
  
}

bool load_shower_profiles( const char *filename ){
  ifstream fshower(filename);

  TString currentfile;
  
  currentfile.ReadFile( fshower );
  
  //Now the whole file is loaded into the TString:

  int istart, istop, length;

  std::istringstream sstream_temp;
  
  TString stemp,substrtemp;
  stemp = "ECAL_shower_profile_x_default";

  //Read to the next newline character:
  istart = currentfile.Index(stemp);
  istop = currentfile.Index( '\n', istart ); 

  if( istop < 0 ) istop = currentfile.Length()-1;
  
  substrtemp = TString(currentfile(istart,istop-istart));

  cout << "'" << substrtemp << "'" << endl;

  substrtemp.ReplaceAll(stemp,"");

  cout << "'" << substrtemp << "'" << endl;

  //  sstream_temp = std::istringstream(substrtemp.Data());
  sstream_temp.str(std::string(substrtemp.Data()));

  //extract number of bins and min and max moment value for default x profile:
  sstream_temp >> profxdefault.nbins_mom >> profxdefault.xmin >> profxdefault.xmax;
  
  cout << "values extracted from sstream: profxdefault.nbins_mom = " << profxdefault.nbins_mom
       << ", xmin, xmax = " << profxdefault.xmin << ", " << profxdefault.xmax << endl;
  istart = istop;
  istop = currentfile.Index( "ECAL_shower_profile", istart ); //returns position of next header

  if( istop < 0 ) istop = currentfile.Length()-1;
  
  //Grab profile data:
  substrtemp = TString(currentfile(istart+1,istop-istart-1));
  
  cout << "'" << substrtemp << "'" << endl;
  //for( int bin=0; bin<=profxdefault.nbins_mom; bin++ ){
  sstream_temp.clear();
  sstream_temp.str(std::string(substrtemp.Data()));

  cout << "Current string = '" << sstream_temp.str() << "'" << endl;
  
  profxdefault.frac.resize(profxdefault.nbins_mom+1);
  
  //extract bin contents from profile data:
  for( int bin=0; bin<=profxdefault.nbins_mom; bin++ ){
    
    sstream_temp >> profxdefault.frac[bin];
    cout << "bin, frac = " << bin << ", " << profxdefault.frac[bin] << endl;
  }
  
  //Do the same for the default y profile:
  stemp = "ECAL_shower_profile_y_default";

  
								       
  istart = currentfile.Index(stemp);
  istop = currentfile.Index( '\n', istart ); 

  if( istop < 0 ) istop = currentfile.Length()-1;
  
  substrtemp = TString(currentfile(istart,istop-istart));

  cout << "'" << substrtemp << "'" << endl;

  substrtemp.ReplaceAll(stemp,"");

  cout << "'" << substrtemp << "'" << endl;

  sstream_temp.clear();
  sstream_temp.str(std::string(substrtemp.Data()));

  sstream_temp >> profydefault.nbins_mom >> profydefault.xmin >> profydefault.xmax;

  istart = istop;
  istop = currentfile.Index( "ECAL_shower_profile", istart ); //returns position of next header

  if( istop < 0 ) istop = currentfile.Length()-1;
  
  substrtemp = TString(currentfile(istart+1,istop-istart-1));

  cout << "'" << substrtemp << "'" << endl;
  //for( int bin=0; bin<=profxdefault.nbins_mom; bin++ ){
  sstream_temp.clear();
  sstream_temp.str(std::string(substrtemp.Data()));

  profydefault.frac.resize(profydefault.nbins_mom+1);
    
  for( int bin=0; bin<=profydefault.nbins_mom; bin++ ){
    sstream_temp >> profydefault.frac[bin];
    cout << "bin, frac = " << bin << ", " << profxdefault.frac[bin] << endl;
  }

  //Now read x-dependent x profiles and y-dependent y profiles:
  stemp = "ECAL_shower_profiles_x";

  istart = currentfile.Index( stemp );
  istop = currentfile.Index('\n',istart);

  if( istop < 0 ) istop = currentfile.Length()-1;
  
  substrtemp = TString(currentfile(istart,istop-istart));

  substrtemp.ReplaceAll(stemp,"");

  sstream_temp.clear();
  sstream_temp.str(std::string(substrtemp.Data()));

  sstream_temp >> profx_nbins >> profx_xmin >> profx_xmax;

  profx.resize( profx_nbins );
  
  for( int xbin=0; xbin<profx_nbins; xbin++ ){
    stemp.Form("ECAL_shower_profile_x_bin %d",xbin+1);
    istart = currentfile.Index(stemp);
    istop = currentfile.Index('\n',istart);

    if( istop < 0 ) istop = currentfile.Length()-1;
    
    substrtemp = TString(currentfile(istart,istop-istart));
    substrtemp.ReplaceAll(stemp,"");

    sstream_temp.clear();
    sstream_temp.str(std::string(substrtemp.Data()));

    sstream_temp >> profx[xbin].nbins_mom >> profx[xbin].xmin >> profx[xbin].xmax;

    istart = istop;
    istop = currentfile.Index("ECAL",istart);

    if( istop < 0 ) istop = currentfile.Length()-1;
    
    substrtemp = TString(currentfile(istart+1,istop-istart-1));

    sstream_temp.clear();
    sstream_temp.str(std::string(substrtemp.Data()));

    profx[xbin].frac.resize(profx[xbin].nbins_mom+1);

    for( int bin=0; bin<=profx[xbin].nbins_mom; bin++ ){
      sstream_temp >> profx[xbin].frac[bin];

      cout << "xbin, bin, frac = " << xbin << ", " << bin << ", " << profx[xbin].frac[bin] << endl;
    }
  }

  ///y-dependent y profiles:
  stemp = "ECAL_shower_profiles_y";

  istart = currentfile.Index( stemp );
  istop = currentfile.Index('\n',istart);

  if( istop < 0 ) istop = currentfile.Length()-1;
  
  substrtemp = TString(currentfile(istart,istop-istart));

  substrtemp.ReplaceAll(stemp,"");

  sstream_temp.clear();
  sstream_temp.str(std::string(substrtemp.Data()));

  sstream_temp >> profy_nbins >> profy_ymin >> profy_ymax;

  profy.resize( profy_nbins );
  
  for( int ybin=0; ybin<profy_nbins; ybin++ ){
    stemp.Form("ECAL_shower_profile_y_bin %d",ybin+1);
    istart = currentfile.Index(stemp);
    istop = currentfile.Index('\n',istart);

    if( istop < 0 ) istop = currentfile.Length()-1;
    
    substrtemp = TString(currentfile(istart,istop-istart));
    substrtemp.ReplaceAll(stemp,"");

    sstream_temp.clear();
    sstream_temp.str(std::string(substrtemp.Data()));

    sstream_temp >> profy[ybin].nbins_mom >> profy[ybin].xmin >> profy[ybin].xmax;

    istart = istop;
    istop = currentfile.Index("ECAL",istart);

    cout << "ybin, istop = " << ybin << ", " << istop << endl;

    if( istop < 0 ) istop = currentfile.Length()-1;
    
    substrtemp = TString(currentfile(istart+1,istop-istart-1));

    sstream_temp.clear();
    sstream_temp.str(std::string(substrtemp.Data()));

    profy[ybin].frac.resize(profy[ybin].nbins_mom+1);

    for( int bin=0; bin<=profy[ybin].nbins_mom; bin++ ){
      sstream_temp >> profy[ybin].frac[bin];

      cout << "ybin, bin, frac = " << ybin << ", " << bin << ", " << profy[ybin].frac[bin] << endl;
    }
  }
  
  return true;
}

void SBS_tgt_reconstruct( double xfp, double yfp, double xpfp, double ypfp, double xtar, double &xptar, double &yptar, double &ytar, double &p ){
  double sum_xptar=0.0;
  double sum_yptar=0.0;
  double sum_ytar=0.0;
  double sum_ptheta=0.0;

  for( int i=0; i<nterms_optics; i++ ){
    double term = pow(xfp,Cexpon[i][0]) * pow(yfp,Cexpon[i][1]) * pow(xpfp,Cexpon[i][2])
      * pow(ypfp,Cexpon[i][3])*pow(xtar,Cexpon[i][4]);
    sum_xptar += Cxptar[i]*term;
    sum_yptar += Cyptar[i]*term;
    sum_ytar += Cytar[i]*term;
    sum_ptheta += Cpthetabend[i]*term;
  }

  TVector3 zaxis_FP(-sin(SBS_tracker_pitch),0,cos(SBS_tracker_pitch));
  TVector3 yaxis_FP(0,1,0);
  TVector3 xaxis_FP = yaxis_FP.Cross(zaxis_FP).Unit();

  TVector3 nhat_fp(xpfp,ypfp,1.0);
  nhat_fp = nhat_fp.Unit();

  TVector3 nhat_fp_global = nhat_fp.X() * xaxis_FP + nhat_fp.Y() * yaxis_FP + nhat_fp.Z() * zaxis_FP;
  TVector3 nhat_tgt(xptar,yptar,1.0);
  nhat_tgt = nhat_tgt.Unit();

  double thetabend = acos( nhat_fp_global.Dot( nhat_tgt ) );

  p = sum_ptheta/thetabend;

  xptar = sum_xptar;
  yptar = sum_yptar;
  ytar = sum_ytar;
}

void SBS_fp_reconstruct( double xtar, double ytar, double xptar, double yptar, double p, double &xfp, double &yfp, double &xpfp, double &ypfp ){
  double sum_xfp = 0.0;
  double sum_yfp = 0.0;
  double sum_xpfp = 0.0;
  double sum_ypfp = 0.0;

  for( int i=0; i<nterms_foptics; i++ ){
    double term = pow(xptar, Cfexpon[i][0]) * pow(yptar,Cfexpon[i][1]) * pow(ytar,Cfexpon[i][2]) * pow(1.0/p,Cfexpon[i][3]) * pow(xtar, Cfexpon[i][4] );
    sum_xfp += Cxfp[i]*term;
    sum_yfp += Cyfp[i]*term;
    sum_xpfp += Cxpfp[i]*term;
    sum_ypfp += Cypfp[i]*term;
  }
  xfp = sum_xfp;
  yfp = sum_yfp;
  xpfp = sum_xpfp;
  ypfp = sum_ypfp;
}

void find_ECAL_clusters( gep_tree_with_spin *T, int &nclust, vector<double> &xclust, vector<double> &yclust, vector<double> &Eclust, vector<int> &nhitclust, vector<vector<int> > &hitlist_clust ){ //This method is only intended to group ECAL hits together in clusters and do "crude" coordinate reconstuction (shower center-of-gravity)

  nclust = 0;
  xclust.clear();
  yclust.clear();
  Eclust.clear();
  nhitclust.clear();
  //nxclust.clear();
  //nyclust.clear();
  hitlist_clust.clear();
  
  double Eclustmin = TMath::Max(1.0,0.5*T->ev_ep); //50% of elastic or 1 GeV, whichever is bigger

  bool foundclust = false;

  int nhitstot = T->Earm_ECalTF1_hit_nhits;

  //cout << "ECAL cluster finding, nhits = " << nhitstot << endl;

  set<int> unused_hits;
  //vector<bool> hitused(nhitstot);
  vector<double> Ehit_recon(nhitstot);
  
  for( int hit=0; hit<nhitstot; hit++ ){
    //hitused[hit] = false;
    unused_hits.insert(hit);

    double Etemp = (*(T->Earm_ECalTF1_hit_sumedep))[hit];
    double nphemean = Etemp*ECAL_phe_per_GeV;
    
    double nphe = random_generator->PoissonD(nphemean);

    Ehit_recon[hit] = nphe/ECAL_phe_per_GeV;

    //cout << "(ihit, Etrue, Erecon)=(" << hit << ", " << Etemp << ", " << Ehit_recon[hit] << ")" << endl;
    
  }
  
  //repeat cluster search until no new clusters found:
  do {
    foundclust = false;
    //Step 1: loop over all unused hits; find maximum.
    double Ehitmax = 0.0;
    int ihitmax = -1; //position in hit array of hit with largest energy.

    for( set<int>::iterator hit=unused_hits.begin(); hit != unused_hits.end(); ++hit ){
      int ihit = *hit;

      double Ehit = Ehit_recon[ihit];

      ihitmax = Ehit > Ehitmax ? ihit : ihitmax;
      Ehitmax = Ehit > Ehitmax ? Ehit : Ehitmax;
      
    }

    //cout << "(ihitmax, Ehitmax)=(" << ihitmax << ", " << Ehitmax << ")" << endl;

    if( ihitmax < 0 ) {
      foundclust = false;
    } else { //found a new maximum: start a cluster around this maximum
      //unused_hits.erase( ihitmax );

      int ncellclust_temp = 1;
      vector<int> hitlist_temp;
      hitlist_temp.push_back( ihitmax );

      double sum_logE = log(Ehitmax);
      double Eclust_temp = Ehitmax;
      double xclust_temp = (*(T->Earm_ECalTF1_hit_xcell))[ihitmax]*Ehitmax;
      double yclust_temp = (*(T->Earm_ECalTF1_hit_ycell))[ihitmax]*Ehitmax;
      //int nxclust_temp = 1;
      //int nyclust_temp = 1;
      
      int icelltemp = 0;
      //Next, we want a do while loop over all hits in the cluster to add "nearest neighbors"
      while ( icelltemp < hitlist_temp.size() ) {

	// unused_hits.erase( hitlist_temp[icelltemp] );
	// at the beginning of each iteration of finding nearest-neighbors, mark all hits
	// associated with this cluster as used, starting with icelltemp
	// everything preceding icelltemp should have already been removed!
	for( int jhit=icelltemp; jhit<hitlist_temp.size(); jhit++ ){
	  unused_hits.erase( hitlist_temp[jhit] ); 
	}
	//grab information about the current cell:
	int rowtemp = (*(T->Earm_ECalTF1_hit_row))[hitlist_temp[icelltemp]];
	int coltemp = (*(T->Earm_ECalTF1_hit_col))[hitlist_temp[icelltemp]];
	double xcelltemp = (*(T->Earm_ECalTF1_hit_xcell))[hitlist_temp[icelltemp]];
	double ycelltemp = (*(T->Earm_ECalTF1_hit_ycell))[hitlist_temp[icelltemp]];
	double Ecelltemp = Ehit_recon[hitlist_temp[icelltemp]];
	//loop over the list of unused hits, adding any nearest-neighbor hits found:
	//note that when we start this loop, the central maximum of the cluster has already been removed from the
	//set of unused hits!

	//Do NOT modify the set while iterating through it!
	for( set<int>::iterator hit=unused_hits.begin(); hit != unused_hits.end(); ++hit ){
	  int ihit = *hit;

	  //check for same row and column +/- 1, or row +/- 1 and (x - xcell)

	  double xhit = (*(T->Earm_ECalTF1_hit_xcell))[ihit];
	  double yhit = (*(T->Earm_ECalTF1_hit_ycell))[ihit];
	  double Ehit = Ehit_recon[ihit];

	  if( pow(xhit-xcelltemp,2)+pow(yhit-ycelltemp,2)<=pow(1.5*ECAL_max_cell_size,2) ){ //nearest-neighbor:
	    hitlist_temp.push_back( ihit );
	    xclust_temp += xhit*Ehit;
	    yclust_temp += yhit*Ehit;
	    Eclust_temp += Ehit;
	    sum_logE += log(Ehit);
	  }
	}

	icelltemp++; //after first iteration, which always happens, icelltemp = 1.
	// if any nearest-neighbors were found, hitlist_temp.size() > 1
	// we keep looping over all hits in the list of hits added to this cluster until
	// no more unused nearest-neighbor hits are found!
	// As long as at least one new nearest neighbor was found, we keep going!
	// As soon as we don't find any new neighbors, the condition below fails and we exit the loop!
	
      } 

      //Add this cluster to the output arrays:

      if( Eclust_temp >= Eclustmin ){
      
	xclust.push_back( xclust_temp/Eclust_temp );
	yclust.push_back( yclust_temp/Eclust_temp );
	Eclust.push_back( Eclust_temp );
	nhitclust.push_back( hitlist_temp.size() );
	hitlist_clust.push_back( hitlist_temp );
      
	foundclust = true;

	// cout << "found cluster E, x, y = " << Eclust_temp << ", " << xclust_temp/Eclust_temp
	//      << ", " << yclust_temp/Eclust_temp << endl;
      }
    }
  } while( foundclust );
  
  nclust = xclust.size();

}

void Fit_3D_track( vector<double> xpoints, vector<double> ypoints, vector<double> zpoints, vector<double> wx, vector<double> wy,
		   double &X, double &Y, double &Xp, double &Yp ){

  //Setting up the matrices for the fit:
  
  TMatrixD A(4,4);
  TVectorD b(4);

  for( int i=0; i<4; i++ ){
    for( int j=0; j<4; j++ ){
      A(i,j) = 0.0;
    }
    b(i) = 0.0;
  }

  int npoints = xpoints.size();

  if( npoints<2 ) return;

  int ndf=0;

  //For a 3D fit to a straight-line:
  // chi^2 = sum_i wxi * (xi- (X + Xp*zi))^2 + wyi*(y - (Y+Yp*zi))^2
  // dchi^2/dX = -2 * (xi - (X+Xp*zi))* wxi = 0
  // dchi^2/dY = -2 * (yi - (Y+Yp*zi))* wyi = 0
  // dchi^2/dXp = -2 * (xi - (X+Xp*zi))*zi * wxi = 0
  // dchi^2/dYp = -2 * (yi - (Y+Yp*zi))*zi * wyi = 0
  for( int i=0; i<npoints; i++ ){
    b(0) += wx[i]*xpoints[i];
    b(1) += wy[i]*ypoints[i];
    b(2) += wx[i]*xpoints[i]*zpoints[i];
    b(3) += wy[i]*ypoints[i]*zpoints[i];

    A(0,0) += wx[i];
    A(0,1) += 0.0;
    A(0,2) += wx[i]*zpoints[i];
    A(0,3) += 0.0;

    A(1,0) += 0.0;
    A(1,1) += wy[i];
    A(1,2) += 0.0;
    A(1,3) += wy[i]*zpoints[i];

    A(2,0) += wx[i]*zpoints[i];
    A(2,1) += 0.0;
    A(2,2) += wx[i]*pow(zpoints[i],2);
    A(2,3) += 0.0;

    A(3,0) += 0.0;
    A(3,1) += wy[i]*zpoints[i];
    A(3,2) += 0.0;
    A(3,3) += wy[i]*pow(zpoints[i],2);
    
  }

  TVectorD solution = A.Invert() * b;

  X = solution(0);
  Y = solution(1);
  Xp = solution(2);
  Yp = solution(3);
}

void ElasticEP_Correlation_Analysis(const char *inputfilename, const char *outputfilename="temp.root"){
  //start by reading list of files and creating the TChain:

  random_generator = new TRandom3(0);
  
  TFile *fout = new TFile(outputfilename,"RECREATE");

  TH1D *hdpp = new TH1D("hdpp","",100,-0.1,0.1); //fraction of p
  TH1D *hdptheta = new TH1D("hdptheta","",100,-0.01,0.01); //radians
  TH1D *hdpphi   = new TH1D("hdpphi","",100,-0.03,0.03);
  TH1D *hpmissp = new TH1D("hpmissp","",100,-500.0,500.0); //MeV
  TH1D *hpmissp_fractional = new TH1D("hpmissp_fractional","",100,-0.1,0.1);
  TH1D *hpmissp_fractional_cut = new TH1D("hpmissp_fractional_cut","",100,-0.1,0.1);
  TH1D *hpmissp_fractional_anticut = new TH1D("hpmissp_fractional_anticut","",100,-0.1,0.1);
  TH1D *hdvz = new TH1D("hdvz","",100,-50.,50.); //mm

  TH1D *htrackchi2ndf_fit = new TH1D("htrackchi2ndf_fit","",100,0.0,10.0);
  TH1D *htrackchi2ndf_true = new TH1D("htrackchi2ndf_true","",100,0.0,1.0);
  TH1D *htracknhits = new TH1D("htracknhits","",8,-0.5,7.5);

  htrackchi2ndf_fit->GetXaxis()->SetTitle("Track #chi^{2}/NDF, smeared hit positions");
  htrackchi2ndf_true->GetXaxis()->SetTitle("Track #chi^{2}/NDF, true hit positions");
  htracknhits->GetXaxis()->SetTitle("Number of GEM layers with hits on track");
  
  hdpp->GetXaxis()->SetTitle("p_{p}(recon)/p_{p}(true)-1");
  hdptheta->GetXaxis()->SetTitle("#theta_{p}(recon)-#theta_{p}(true) (rad)");
  hdpphi->GetXaxis()->SetTitle("#phi_{p}(recon)-#phi_{p}(true) (rad)");
  hpmissp->GetXaxis()->SetTitle("p_{p}(recon)-p_{el}(#theta_{p}(recon)) (MeV/c)");
  hpmissp_fractional->GetXaxis()->SetTitle("(p_{p}(recon)-p_{el}(#theta_{p}(recon)))/p_{true}");
  hdvz->GetXaxis()->SetTitle("vertex z - z_{true} (mm)");

  TH1D *hnclust = new TH1D("hnclust","",11,-0.5,10.5);
  TH2D *hxyclust = new TH2D("hxyclust","",100,-0.75,0.75,100,-1.75,1.75); //m
  TH1D *hEclust = new TH1D("hEclust","",100,0.0,6.0);
  TH2D *hEclust_vs_Etrue = new TH2D("hEclust_vs_Etrue","",100,0.0,6.0,100,0.0,6.0);
  TH1D *hdEclust_fractional = new TH1D("hdEclust_fractional","",100,-1.0,1.0);
  TH2D *hdEclust_fractional_vs_Eclust = new TH2D("hdEclust_fractional_vs_Eclust","",100,0.0,6.0,100,-1.0,1.0);
  TH1D *hnhitclust = new TH1D("hnhitclust","",50,0.5,50.5);
  TH1D *hdxclust = new TH1D("hdxclust","",100,-10.,10.); //cm, xclust - x "true"
  TH1D *hdyclust = new TH1D("hdyclust","",100,-10.,10.); //cm, yclust - y "true"
  TH1D *hxmom = new TH1D("hxmom","",100,-1.0,1.0);
  TH1D *hymom = new TH1D("hymom","",100,-1.0,1.0);
  TH2D *hdxx = new TH2D("hdxx","",300,-75.0,75.0,100,-10.,10.);
  TH2D *hdyy = new TH2D("hdyy","",700,-175.,175.,100,-10.,10.);

  //For shower profile mapping:
  TH2D *hxmom_x = new TH2D("hxmom_x","",5,-70.,70.,100,-1.0,1.0);
  TH2D *hymom_y = new TH2D("hymom_y","",9,-150.,150.,100,-1.0,1.0);

  TH2D *hxyclust_corrected = new TH2D("hxyclust_corrected","",100,-0.75,0.75,100,-1.75,1.75);
  TH1D *hdxclust_corrected = new TH1D("hdxclust_corrected","",100,-10.0,10.0);
  TH1D *hdyclust_corrected = new TH1D("hdyclust_corrected","",100,-10.0,10.0);

  TH2D *hdxx_corrected = new TH2D("hdxx_corrected","",300,-75.,75.,100,-10.,10.);
  TH2D *hdyy_corrected = new TH2D("hdyy_corrected","",700,-175.,175.,100,-10.,10.);

  TH1D *hxcorr_xmax = new TH1D("hxcorr_xmax","",100,-ECAL_max_cell_size,ECAL_max_cell_size);
  TH1D *hycorr_ymax = new TH1D("hycorr_ymax","",100,-ECAL_max_cell_size,ECAL_max_cell_size);
  
  //TH1D *hxcl_xmax_corr = new TH1D("hxcl_xmax_corr","",100,-ECAL_max_

  TH1D *hnplanes_CDET = new TH1D("hnplanes_CDET","",3,-0.5,2.5);
  TH1D *hnhits1_CDET = new TH1D("hnhits1_CDET","",6,-0.5,5.5);
  TH1D *hnhits2_CDET = new TH1D("hnhits2_CDET","",6,-0.5,5.5);
  TH1D *hdy1_CDET_ECAL = new TH1D("hdy1_CDET_ECAL","",100,-0.05,.05);
  TH1D *hdy2_CDET_ECAL = new TH1D("hdy2_CDET_ECAL","",100,-0.05,.05);
  TH1D *hdy12_CDET = new TH1D("hdy12_CDET","",100,-0.05,.05);

  TH1D *hEtot1_CDET = new TH1D("hEtot1_CDET","",100,0.0,0.05);
  TH1D *hEtot2_CDET = new TH1D("hEtot2_CDET","",100,0.0,0.05);

  TH1D *hdEclust_Eprime_eth = new TH1D("hdEclust_Eprime_eth","",100,-1.0,1.0);
  TH1D *hdEclust_Eprime_eth_cut = new TH1D("hdEclust_Eprime_eth_cut","",100,-1.0,1.0);
  TH1D *hdEclust_Eprime_eth_anticut = new TH1D("hdEclust_Eprime_eth_anticut","",100,-1.0,1.0);
  
  TH1D *hdethc_Eclust = new TH1D("hdethc_Eclust","",100,-0.2,0.2);
  TH1D *hdethc_thtrue = new TH1D("hdethc_thtrue","",100,-0.03,0.03);
  TH2D *hdethc_thtrue_vs_z = new TH2D("hdethc_thtrue_vs_z","",100,-0.55*Ltgt,0.55*Ltgt,100,-.03,.03);
  TH1D *hdephc_phtrue = new TH1D("hdephc_phtrue","",100,-0.03,0.03);

  TH1D *hdxfp = new TH1D("hdxfp","",100,-0.1,0.1);
  TH1D *hdyfp = new TH1D("hdyfp","",100,-0.1,0.1);
  TH1D *hdxpfp = new TH1D("hdxpfp","",100,-0.02,0.02);
  TH1D *hdypfp = new TH1D("hdypfp","",100,-0.02,0.02);

  TH2D *hdxfp_vs_ztrue = new TH2D("hdxfp_vs_ztrue","",100,-0.55*Ltgt,0.55*Ltgt, 100, -0.1,0.1);
  TH2D *hdyfp_vs_ztrue = new TH2D("hdyfp_vs_ztrue","",100,-0.55*Ltgt,0.55*Ltgt, 100, -0.1,0.1);
  TH2D *hdxpfp_vs_ztrue = new TH2D("hdxpfp_vs_ztrue","",100,-0.55*Ltgt,0.55*Ltgt, 100, -0.02,0.02);
  TH2D *hdypfp_vs_ztrue = new TH2D("hdypfp_vs_ztrue","",100,-0.55*Ltgt,0.55*Ltgt, 100, -0.02,0.02);
  
  

  TH1D *hdpp_etheta_c = new TH1D("hdpp_etheta_c","",100,-0.2,0.2);
  TH1D *hdpp_eclust = new TH1D("hdpp_eclust","",100,-0.2,0.2);
  TH1D *hdpp_etheta_vtx_known = new TH1D("hdpp_etheta_vtx_known","",100,-0.2,0.2);

  TH1D *hdxfp_vtx_known = new TH1D("hdxfp_vtx_known","",100,-0.1,0.1);
  TH1D *hdyfp_vtx_known = new TH1D("hdyfp_vtx_known","",100,-0.1,0.1);
  TH1D *hdxpfp_vtx_known = new TH1D("hdxpfp_vtx_known","",100,-0.02,0.02);
  TH1D *hdypfp_vtx_known = new TH1D("hdypfp_vtx_known","",100,-0.02,0.02);

  //Finally, let's look at exclusivity cut variables computed from reconstructed electron arm and proton arm
  //quantities:

  TH1D *hdeltax_pth = new TH1D("hdeltax_pth","",100,-30.0,30.0);
  TH1D *hdeltay_pth = new TH1D("hdeltay_pth","",100,-30.0,30.0);
  TH2D *hdxdy_pth = new TH2D("hdxdy_pth","",100,-30.0,30.0,100,-30.0,30.0);

  TH1D *hdeltax_pth_cut = new TH1D("hdeltax_pth_cut","",100,-30.0,30.0);
  TH1D *hdeltay_pth_cut = new TH1D("hdeltay_pth_cut","",100,-30.0,30.0);
  TH2D *hdxdy_pth_cut = new TH2D("hdxdy_pth_cut","",100,-30.0,30.0,100,-30.0,30.0);

  TH1D *hdeltax_pth_anticut = new TH1D("hdeltax_pth_anticut","",100,-30.0,30.0);
  TH1D *hdeltay_pth_anticut = new TH1D("hdeltay_pth_anticut","",100,-30.0,30.0);
  TH2D *hdxdy_pth_anticut = new TH2D("hdxdy_pth_anticut","",100,-30.0,30.0,100,-30.0,30.0);

  
  
  TH1D *hdeltax_pp = new TH1D("hdeltax_pp","",100,-30.0,30.0);
  TH1D *hdeltay_pp = new TH1D("hdeltay_pp","",100,-30.0,30.0);
  TH2D *hdxdy_pp = new TH2D("hdxdy_pp","",100,-30.0,30.0,100,-30.,30.);

  // TH1D *hdeltax_4v = new TH1D("hdeltax_4v","",100,-30.,30.);
  // TH1D *hdeltay_4v = new TH1D("hdeltay_4v","",100,-30.,30.);
  // TH2D *hdxdy_4v = new TH2D("hdxdy_4v","",100,-30.,30.,100,-30.,30.);

  TH1D *hdphip = new TH1D("hdphip","",100,-0.02,0.02);
  TH1D *hdthetap = new TH1D("hdthetap","",100,-0.02,0.02);
  TH1D *hpmisse = new TH1D("hpmisse","",100,-0.1,0.1);
  TH1D *hdpp_eth_pth = new TH1D("hdpp_eth_pth","",100,-0.1,0.1);
  TH1D *hdthetae_true = new TH1D("hdthetae_true","",100,-0.02,0.02);
  TH1D *hdphie_true = new TH1D("hdphie_true","",100,-0.02,0.02);
    
  TChain *C = new TChain("T");
  
  ifstream infile(inputfilename);
  
  TString currentline;

  //What information is needed to analyze the EP kinematic correlation?
  G4SBSRunData *rd;
  TFile *ftemp;

  map<TString,double> Ebeam_file;
  map<TString,double> SBStheta_file;
  map<TString,double> ECALtheta_file;
  map<TString,double> ECALdist_file;
  //map<TString,long> ngen_file;
  map<TString,double> GenVol_file;
  map<TString,double> Luminosity_file;
  
  long ngen_total=0;
  
  double SBStheta_default = 16.9*PI/180.0;
  double ECALtheta_default = 29.0*PI/180.0;
  double ECALdist_default = 4.5;
  double Ebeam_default = 11.0;
  double GenVol_default = 1.0;
  double Luminosity_default = 8e38;

  int pythia6_flag=0;
  
  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){ //not commented out:
      //Only add the file to the chain if we can load and read the rundata object:
      
      ftemp = new TFile(currentline.Data(),"READ");

      if( !ftemp->IsZombie() ){
	ftemp->GetObject("run_data",rd);
	
	if( rd ){
	  Ebeam_file[currentline] = rd->fBeamE; // GeV
	  SBStheta_file[currentline] = rd->fSBStheta; //radians
	  ECALtheta_file[currentline] = rd->fBBtheta; //radians
	  ECALdist_file[currentline] = rd->fBBdist; //meters = distance from origin to front surface of ECAL

	  //ngen_file[currentline] = rd->fNtries;
	  GenVol_file[currentline] = rd->fGenVol;
	  Luminosity_file[currentline] = rd->fLuminosity;
	  ngen_total += rd->fNtries;
	  //TODO: add SBS tracker pitch angle and/or distance to list of things we grab from the
	  //ROOT tree directly (after re-running simulation with correct default values in run_data).
	  
	  C->Add(currentline.Data());

	  Ebeam_default = Ebeam_file[currentline];
	  SBStheta_default = SBStheta_file[currentline];
	  ECALtheta_default = ECALtheta_file[currentline];
	  ECALdist_default = ECALdist_file[currentline];
	  GenVol_default = GenVol_file[currentline];
	  Luminosity_default = Luminosity_file[currentline];
	  
	  
	}
      }

    }
  }

  TCut global_cut = "";
  
  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      global_cut += currentline;
    }
  }

  if(C->GetNtrees() <= 0 ) return;

  TEventList *elist = new TEventList("elist");

  C->Draw(">>elist",global_cut);
  
  TString opticsfilename="",fopticsfilename="";

  TString ECALshowerprofile_filename="";

  double dxmin=-4.0, dxmax=4.0;
  double dymin=-3.0, dymax=3.0;
  double dppmin=-0.03, dppmax=0.03;
  double dEmin=-0.2, dEmax=0.2;
  double Eclustmin = 2.0;
  
  while( currentline.ReadLine(infile) && !currentline.BeginsWith("endconfig") ){
    if( !currentline.BeginsWith("#") ){
      TObjArray *tokens = currentline.Tokenize(" ");

      int ntokens = tokens->GetEntries();
      if( ntokens >= 2 ){
      
	TString keyword = ( (TObjString*) (*tokens)[0] )->GetString();

	if( keyword.BeginsWith("forward_optics") ){
	  fopticsfilename = ( (TObjString*) (*tokens)[1] )->GetString();
	}

	if( keyword.BeginsWith("optics") ){
	  opticsfilename = ( (TObjString*) (*tokens)[1] )->GetString();
	}

	if( keyword.BeginsWith("ECAL_shower_profiles") ){
	  ECALshowerprofile_filename = ( (TObjString*) (*tokens)[1] )->GetString();
	}

	if( keyword.BeginsWith("pythia_flag") ){
	  pythia6_flag = ( (TObjString*) (*tokens)[1] )->GetString().Atoi();
	}

	if( keyword.BeginsWith("dxcut") && ntokens >= 3 ){
	  dxmin = ( (TObjString*) (*tokens)[1] )->GetString().Atof();
	  dxmax = ( (TObjString*) (*tokens)[2] )->GetString().Atof();
	}

	if( keyword.BeginsWith("dycut") && ntokens >= 3 ){
	  dymin = ( (TObjString*) (*tokens)[1] )->GetString().Atof();
	  dymax = ( (TObjString*) (*tokens)[2] )->GetString().Atof();
	}

	if( keyword.BeginsWith("dppcut") && ntokens >= 3 ){
	  dppmin = ( (TObjString*) (*tokens)[1] )->GetString().Atof();
	  dppmax = ( (TObjString*) (*tokens)[2] )->GetString().Atof();
	}

	if( keyword.BeginsWith("dEcut") && ntokens >= 3 ){
	  dEmin = ( (TObjString*) (*tokens)[1] )->GetString().Atof();
	  dEmax = ( (TObjString*) (*tokens)[2] )->GetString().Atof();
	}

	if( keyword.BeginsWith("ECAL_threshold") ){
	  Eclustmin = ( (TObjString*) (*tokens)[1] )->GetString().Atof();
	}
      }
    }
  }

  nterms_optics = 0;
  nterms_foptics = 0;

  Cxptar.clear();
  Cyptar.clear();
  Cytar.clear();
  Cpthetabend.clear();

  Cexpon.clear();

  Cxfp.clear();
  Cyfp.clear();
  Cxpfp.clear();
  Cypfp.clear();

  Cfexpon.clear();
  
  if( opticsfilename != "" ){
    ifstream opticsfile(opticsfilename.Data());

    while(currentline.ReadLine( opticsfile ) ){ //Only concern ourselves with lines corresponding to actual terms in the expansion:
      if( !currentline.BeginsWith("#") ){

	//cout << currentline << endl;
	
	double Cxptemp, Cyptemp, Cytemp, Cpthtemp;

	vector<int> etemp(5);

	if( sscanf( currentline.Data(), " %lg %lg %lg %lg %d %d %d %d %d",
		    &Cxptemp, &Cyptemp, &Cytemp, &Cpthtemp,
		    &(etemp[0]), &(etemp[1]), &(etemp[2]), &(etemp[3]), &(etemp[4]) ) == 9 ){

	  Cxptar.push_back(Cxptemp);
	  Cyptar.push_back(Cyptemp);
	  Cytar.push_back(Cytemp);
	  Cpthetabend.push_back(Cpthtemp);

	  Cexpon.push_back(etemp);

	  nterms_optics++;
	}
      }
    }
  }

  if( fopticsfilename != "" ){
    ifstream fopticsfile(fopticsfilename.Data());

    while(currentline.ReadLine( fopticsfile ) ){ //Only concern ourselves with lines corresponding to actual terms in the expansion:
      if( !currentline.BeginsWith("#") ){

	//cout << currentline << endl;
	
	double Cxptemp, Cyptemp, Cxtemp, Cytemp;

	vector<int> etemp(5);

	if( sscanf( currentline.Data(), " %lg %lg %lg %lg %d %d %d %d %d",
		    &Cxtemp, &Cytemp, &Cxptemp, &Cyptemp,
		    &(etemp[0]), &(etemp[1]), &(etemp[2]), &(etemp[3]), &(etemp[4]) ) == 9 ){

	  Cxfp.push_back(Cxtemp);
	  Cyfp.push_back(Cytemp);
	  Cxpfp.push_back(Cxptemp);
	  Cypfp.push_back(Cyptemp);

	  Cfexpon.push_back(etemp);

	  nterms_foptics++;
	}
      }
    }
  }

  if( nterms_optics <= 0 || nterms_foptics <= 0 ) return;

  cout << "n terms, reconstruction optics = " << nterms_optics << endl;
  cout << "n terms, forward optics = " << nterms_foptics << endl;

  bool shower_initialized = false;
  
  if( ECALshowerprofile_filename != "" ){
    
    shower_initialized = load_shower_profiles( ECALshowerprofile_filename );
    
  }
  
  //C->Print();

  gep_tree_with_spin *T = new gep_tree_with_spin(C);

  long nevent = 0;

  int treenum=0;
  int oldtreenum=-1;

  //  double SBStheta, ECALtheta;

  //Vectors needed for coordinate transformations
  TVector3 SBS_xaxis_global(0,-1,0);
  TVector3 SBS_zaxis_global(-sin(SBStheta_default),0,cos(SBStheta_default));
  TVector3 SBS_yaxis_global = SBS_zaxis_global.Cross(SBS_xaxis_global).Unit();

  //Should we adopt a "transport-like" coordinate system for ECAL? For now, no, go with GEANT4-like coordinate system:
  TVector3 ECAL_zaxis(sin(ECALtheta_default),0,cos(ECALtheta_default));
  TVector3 ECAL_yaxis(0,1,0);
  TVector3 ECAL_xaxis = ECAL_yaxis.Cross(ECAL_zaxis).Unit();
  TVector3 R0_ECAL = ECALdist_default * ECAL_zaxis; //Coordinates of ECAL center in global coordinate system

  double Ebeam = Ebeam_default;
  double SBStheta = SBStheta_default;
  double ECALtheta = ECALtheta_default;
  double ECALdist = ECALdist_default;
  double GenVol = GenVol_default;
  double Luminosity = Luminosity_default;
  
  double weight = GenVol*Luminosity/double(ngen_total);
  long ngen_temp = 0;

  double nevents_all = 0.0;
  double nevents_cut = 0.0;
  
  double sum_weights_all = 0.0;
  double sum_weights_cut = 0.0;

  double sum2_weights_all = 0.0;
  double sum2_weights_cut = 0.0;
  
  while( T->GetEntry(elist->GetEntry(nevent++) ) ){
    if( nevent % 1000 == 0 ) cout << nevent << endl;
    
    treenum = C->GetTreeNumber();
    if( treenum != oldtreenum ){

      TString fnametemp = C->GetFile()->GetName();

      if( Ebeam_file.find(fnametemp) != Ebeam_file.end() ){
	Ebeam = Ebeam_file[fnametemp];
	SBStheta = SBStheta_file[fnametemp];
	ECALtheta = ECALtheta_file[fnametemp];
	ECALdist = ECALdist_file[fnametemp];

	//update three-vectors:
	SBS_zaxis_global.SetXYZ(-sin(SBStheta),0,cos(SBStheta) );
	SBS_yaxis_global = SBS_zaxis_global.Cross(SBS_xaxis_global).Unit();

	ECAL_zaxis.SetXYZ(sin(ECALtheta),0,cos(ECALtheta));
	ECAL_xaxis = ECAL_yaxis.Cross(ECAL_zaxis).Unit();
	R0_ECAL = ECALdist * ECAL_zaxis;

	GenVol = GenVol_file[fnametemp];
	Luminosity = Luminosity_file[fnametemp];
      }

      oldtreenum = treenum;
    }

    weight = T->ev_sigma * GenVol * Luminosity/double(ngen_total);

    if( pythia6_flag != 0 ) weight = T->primaries_Sigma * Luminosity/double(ngen_total);

    int itrack_best = 0;
    double ptrue_best = T->ev_np;
    
    
    if( pythia6_flag != 0 ){ // choose highest energy track
      for( int itr=0; itr<T->Harm_FT_Track_ntracks; itr++ ){
	if( itr==0 || (*(T->Harm_FT_Track_P))[itr] > ptrue_best ){
	  ptrue_best = (*(T->Harm_FT_Track_P))[itr];
	  itrack_best = itr;
	}
      }
    }


    //Get primary proton track in FT GEMs (if it's there):
    if( T->Harm_FT_Track_ntracks >= 1 && itrack_best >= 0 && 
	(*(T->Harm_FT_Track_MID))[itrack_best] == 0 ){

      int bestmatch_primary = -1;
      if( pythia6_flag != 0 ){
	for( int ipr=0; ipr<T->Primaries_Nprimaries; ipr++ ){
	  if( ipr == 0 || ( (*(T->Primaries_PID))[ipr] == (*(T->Harm_FT_Track_PID))[itrack_best] &&
			    fabs( (*(T->Primaries_P))[ipr]-ptrue_best ) <= 0.01*ptrue_best ) ){
	    bestmatch_primary = ipr;
	  }
	}
      }
      double xfptemp = (*(T->Harm_FT_Track_Xfit))[itrack_best];
      double yfptemp = (*(T->Harm_FT_Track_Yfit))[itrack_best];
      double xpfptemp = (*(T->Harm_FT_Track_Xpfit))[itrack_best];
      double ypfptemp = (*(T->Harm_FT_Track_Ypfit))[itrack_best];

      double chi2fit = (*(T->Harm_FT_Track_Chi2fit))[itrack_best];
      double chi2true = (*(T->Harm_FT_Track_Chi2true))[itrack_best];

      int NDF = (*(T->Harm_FT_Track_NDF))[itrack_best];

      htrackchi2ndf_fit->Fill( chi2fit/double(NDF), weight );
      htrackchi2ndf_true->Fill( chi2true/double(NDF), weight );
      htracknhits->Fill( double( (*(T->Harm_FT_Track_NumPlanes))[itrack_best] ), weight );
      
      double xtartemp = -T->ev_vy;

      double xptartemp,yptartemp,ytartemp,ptemp;
      int niter=1; //number of iterations of target reconstruction after first iteration for xtar correction:
      for(int iter=0; iter<=niter; iter++){
	SBS_tgt_reconstruct(xfptemp,yfptemp,xpfptemp,ypfptemp,xtartemp,
			    xptartemp,yptartemp,ytartemp,ptemp);
	
	double vztemp = ytartemp / (sin(SBStheta)-yptartemp*cos(SBStheta));

	xtartemp = -T->ev_vy - xptartemp * vztemp * cos(SBStheta);
      }

      double vzrecon = ytartemp/(sin(SBStheta)-yptartemp*cos(SBStheta));
      
      TVector3 nhat_tgt_recon_spec(xptartemp,yptartemp,1.0);
      nhat_tgt_recon_spec = nhat_tgt_recon_spec.Unit();

      TVector3 nhat_tgt_recon_global = nhat_tgt_recon_spec.X() * SBS_xaxis_global +
	nhat_tgt_recon_spec.Y() * SBS_yaxis_global +
	nhat_tgt_recon_spec.Z() * SBS_zaxis_global;

      double pthetarecon = acos( nhat_tgt_recon_global.Z() );
      //Should be centered around +/-PI for SBS on beam right in global coordinate system:
      double pphirecon = TMath::ATan2( nhat_tgt_recon_global.Y(), nhat_tgt_recon_global.X() ); 
      if( pphirecon < 0.0 ) pphirecon += 2.0*PI;

      double pphitrue = (T->ev_nph < 0 ) ? T->ev_nph + 2.0*PI : T->ev_nph;

      if( pythia6_flag != 0 && bestmatch_primary >= 0 ){
	pphitrue = (*(T->Primaries_phi))[bestmatch_primary];

	pphitrue = (pphitrue<0) ? pphitrue + 2.0*PI : pphitrue;
      }
      
      double pprecon = ptemp;

      double ptrue = T->ev_np;
      double thtrue = T->ev_nth;
      if( pythia6_flag != 0 ) ptrue = ptrue_best;
      if( pythia6_flag != 0 && bestmatch_primary >= 0 ) thtrue = (*(T->Primaries_theta))[bestmatch_primary];
      
      hdpp->Fill( pprecon/ptrue - 1.0, weight );
      hdptheta->Fill( pthetarecon - thtrue, weight );
      hdpphi->Fill( pphirecon - pphitrue, weight );

      double pp_ptheta_recon = 2.0*Mp*Ebeam*(Mp+Ebeam)*cos(pthetarecon)/(pow(Mp,2)+2.*Mp*Ebeam + pow(Ebeam*sin(pthetarecon),2));

      hpmissp->Fill( 1000.0*(pprecon-pp_ptheta_recon), weight );
      hpmissp_fractional->Fill( (pprecon-pp_ptheta_recon)/ptrue, weight );
      hdvz->Fill( 1000.0*(vzrecon - T->ev_vz), weight );

      //Next step is to do clustering of hits in ECAL, association with CDET scintillator hits,
      //and electron angle and energy reconstruction. This is the "hard" part

      int nclust;
      vector<double> xclust,yclust,Eclust,xclust_corr,yclust_corr;
      vector<int> nhitclust;
      vector<vector<int> > hitlist_clust;

      find_ECAL_clusters(T,nclust, xclust, yclust, Eclust, nhitclust, hitlist_clust );

      xclust_corr.resize(nclust);
      yclust_corr.resize(nclust);
      
      hnclust->Fill( nclust, weight );

      TVector3 ephat_true_global(T->ev_epx, T->ev_epy, T->ev_epz);
      ephat_true_global = ephat_true_global.Unit();

      TVector3 vertex_true_global(T->ev_vx, T->ev_vy, T->ev_vz );

      TVector3 vertex_true_ECAL( vertex_true_global.Dot(ECAL_xaxis),
				 vertex_true_global.Dot(ECAL_yaxis),
				 vertex_true_global.Dot(ECAL_zaxis) );
      TVector3 ephat_true_ECAL( ephat_true_global.Dot(ECAL_xaxis),
				ephat_true_global.Dot(ECAL_yaxis),
				ephat_true_global.Dot(ECAL_zaxis) );

      double xclust_true = vertex_true_ECAL.X()
	+ ephat_true_ECAL.X()/ephat_true_ECAL.Z() * (ECALdist - vertex_true_ECAL.Z() );
      double yclust_true = vertex_true_ECAL.Y()
	+ ephat_true_ECAL.Y()/ephat_true_ECAL.Z() * (ECALdist - vertex_true_ECAL.Z() );

      //nevermind, don't do that.

      int ibest=-1;

      double Eclust_max = 0.0;
      double xclust_max, yclust_max;
      
      for( int iclust=0; iclust<nclust; iclust++ ){

	ibest = Eclust[iclust] > Eclust_max ? iclust : ibest;
	Eclust_max = Eclust[iclust] > Eclust_max ? Eclust[iclust] : Eclust_max;
	
	hxyclust->Fill(xclust[iclust],yclust[iclust], weight );
	hEclust->Fill(Eclust[iclust], weight );
	hEclust_vs_Etrue->Fill( T->ev_ep, Eclust[iclust], weight );
	hdEclust_fractional->Fill( 1.0-Eclust[iclust]/T->ev_ep, weight );
	hdEclust_fractional_vs_Eclust->Fill( T->ev_ep, 1.0-Eclust[iclust]/T->ev_ep, weight  );
	hnhitclust->Fill(nhitclust[iclust], weight );
	hdxclust->Fill( (xclust[iclust] - xclust_true)*100.0, weight  );
	hdyclust->Fill( (yclust[iclust] - yclust_true)*100.0, weight  );

	double xmom = (xclust[iclust] - (*(T->Earm_ECalTF1_hit_xcell))[hitlist_clust[iclust][0]])/ECAL_max_cell_size;
	double ymom = (yclust[iclust] - (*(T->Earm_ECalTF1_hit_ycell))[hitlist_clust[iclust][0]])/ECAL_max_cell_size;

	double xmax = xclust[iclust] - xmom*ECAL_max_cell_size;
	double ymax = yclust[iclust] - ymom*ECAL_max_cell_size;

	double xf,yf;
	
	hxmom->Fill( xmom, weight  );
	hymom->Fill( ymom, weight  );

	hdxx->Fill( 100.0*xclust_true, 100.0*(xclust[iclust]-xclust_true), weight  );
	hdyy->Fill( 100.0*yclust_true, 100.0*(yclust[iclust]-yclust_true), weight  );

	hxmom_x->Fill( 100.*xclust[iclust], xmom, weight  );
	hymom_y->Fill( 100.*yclust[iclust], ymom, weight  );

	if( shower_initialized ){

	  double xcorrected, ycorrected;
	  calc_shower_coordinates( xmom, ymom, xmax, ymax, Eclust[iclust], ECALdist, xcorrected, ycorrected, xf, yf );
	  //at this point everything is in meters:
	  hxyclust_corrected->Fill( xcorrected, ycorrected, weight  );
	  hdxclust_corrected->Fill( 100.0 * (xcorrected - xclust_true ), weight  );
	  hdyclust_corrected->Fill( 100.0 * (ycorrected - yclust_true ), weight  );

	  hdxx_corrected->Fill( 100.0*xclust_true, 100.0 * (xcorrected - xclust_true ), weight  );
	  hdyy_corrected->Fill( 100.0*yclust_true, 100.0 * (ycorrected - yclust_true ), weight  );

	  hxcorr_xmax->Fill( xf, weight  );
	  hycorr_ymax->Fill( yf, weight  );

	  xclust_corr[iclust] = xcorrected;
	  yclust_corr[iclust] = ycorrected;
	  
	}
	
      }

      //cout << "ibest = " << ibest << endl;
      
      if( ibest >= 0 ){ //Then consider CDET info. We cannot assume tracking info is available, since we want to use ECAL info to aid tracking:

	//Also assume track originates from the origin:

	if( Eclust[ibest] >= Eclustmin ){
	
	  double Eprime_ECAL = Eclust[ibest]/(1.0-dE_E_MPV);
	
	  double xECAL=xclust[ibest], yECAL=yclust[ibest];

	  //cout << "xECAL,yECAL = " << xECAL << ", " << yECAL << endl;
	
	  if( shower_initialized ){
	    xECAL = xclust_corr[ibest];
	    yECAL = yclust_corr[ibest];
	  }

	  TVector3 ECALpos( xECAL, yECAL, ECALdist );
	
	  //ECALpos.Print();
	
	  TVector3 ECALnhat = ECALpos.Unit();

	  TVector3 ECALnhat_global =
	    ECALnhat.X() * ECAL_xaxis +
	    ECALnhat.Y() * ECAL_yaxis +
	    ECALnhat.Z() * ECAL_zaxis;

	  double etheta_temp = acos(ECALnhat_global.Z() );

	  double Eprime_etheta_temp = Ebeam/(1.+Ebeam/Mp*(1.-cos(etheta_temp)));

	  hdEclust_Eprime_eth->Fill( 1.0 - Eprime_ECAL/Eprime_etheta_temp, weight  );
	
	  // track slopes in ECAL local coordinate system, assuming straight line from the origin:
	  double xpECAL = ECALnhat.X()/ECALnhat.Z();
	  double ypECAL = ECALnhat.Y()/ECALnhat.Z(); 
	
	  //Find all CDET hits within some tolerance of the projection of the "ECAL track" back to both CDET planes:

	  int nplanes_fired = 0;

	  int nhits_plane[2] = {0,0};
	  double yavg_plane[2] = {0.,0.};
	  double xavg_plane[2] = {0.,0.};
	  double zplane[2] = {0.,0.};
	  double Etot_plane[2] = {0.,0.};

	  //	cout << "looping over CDET hits...." << endl;
	
	  //require both planes to fire:
	  for( int hit=0; hit<T->Earm_CDET_Scint_hit_nhits; hit++ ){
	    double yhit = (*(T->Earm_CDET_Scint_hit_ycell))[hit];
	    double xhit = (*(T->Earm_CDET_Scint_hit_xcell))[hit];
	    double zhit = (*(T->Earm_CDET_Scint_hit_zcell))[hit] + ECALdist;

	    double Ehit = (*(T->Earm_CDET_Scint_hit_sumedep))[hit];
	  
	    int plane = (*(T->Earm_CDET_Scint_hit_plane))[hit];
	  
	    TVector3 CDETpos( xhit, yhit, zhit );

	    double sigy_total = sqrt(pow(sigy_ECAL,2) + pow(Ly_scint_CDET,2)/12.0);

	    double ytest = yECAL + ypECAL * (zhit - ECALdist);
	    double xtest = xECAL + xpECAL * (zhit - ECALdist);

	    // if( pow(yhit-ytest,2) <= pow(3.0*sigy_total,2) &&
	    //     fabs( xhit - xtest ) <= Lx_scint_CDET/2.0 + 3.0*sigx_ECAL ){
	    if( pow(yhit-ytest,2) <= pow(3.0*sigy_total,2) &&
		fabs( xhit - xtest ) <= Lx_scint_CDET + 3.0*sigx_ECAL ){
	      //good hit:
	      if( nhits_plane[plane-1] == 0 ) nplanes_fired++;
	      nhits_plane[plane-1]++;
	      yavg_plane[plane-1] += yhit * Ehit;
	      xavg_plane[plane-1] += xhit * Ehit;
	      Etot_plane[plane-1] += Ehit;
	      zplane[plane-1] = zhit;
	    }
	  }

	  //	cout << "done looping over CDET hits... " << endl;
	
	  hnplanes_CDET->Fill(nplanes_fired, weight );
	  hnhits1_CDET->Fill(nhits_plane[0], weight );
	  hnhits2_CDET->Fill(nhits_plane[1], weight );
	  if( nplanes_fired == 2 ){
	    for( int plane=0; plane<2; plane++ ){
	      yavg_plane[plane] /= Etot_plane[plane];
	      xavg_plane[plane] /= Etot_plane[plane];
	    }

	    hdy1_CDET_ECAL->Fill( yavg_plane[0] - (yECAL + (zplane[0]-ECALdist)*ypECAL), weight  );
	    hdy2_CDET_ECAL->Fill( yavg_plane[1] - (yECAL + (zplane[1]-ECALdist)*ypECAL), weight  );

	    hdy12_CDET->Fill( yavg_plane[0] - (yavg_plane[1] + (zplane[0]-zplane[1])*ypECAL), weight );

	    hEtot1_CDET->Fill( Etot_plane[0], weight  );
	    hEtot2_CDET->Fill( Etot_plane[1], weight  );

	    //Next step: Fit straight line to ECAL+CDET "tracks" with one point fixed at origin, use forward optics to reconstructed expected FP track position,
	    //compare to actual. Also compute usual exclusivity cut quantities like deltax, deltay, pmisse, etc.
	    // Will probably end up with better resolution if we use proton angles as opposed to proton momentum.

	    //Matrices for 3D fit of tracks based on ECAL cluster plus CDET points:
	  
	    vector<double> xup,yup,zup;
	    vector<double> xdown,ydown,zdown;
	    vector<double> xc,yc,zc;
	    vector<double> xz,yz,zz;

	  
	    vector<double> wxup,wyup,wxdown,wydown,wxc,wyc;
	    vector<double> wxz,wyz;
	  
	    //Compute electron scattering angles for three assumed vertex points:
	    // upstream end of target, center of target, and downstream end of target:
	    TVector3 vertex_up_global(vertex_true_global.X(),vertex_true_global.Y(),-Ltgt/2.0);
	    TVector3 vertex_down_global(vertex_true_global.X(),vertex_true_global.Y(),Ltgt/2.0);
	    TVector3 vertex_center_global(vertex_true_global.X(),vertex_true_global.Y(),0);

       
	  
	    TVector3 vertex_up_ECAL( vertex_up_global.Dot(ECAL_xaxis),
				     vertex_up_global.Dot(ECAL_yaxis),
				     vertex_up_global.Dot(ECAL_zaxis) );

	    //assuming vertex to be at upstream end:
	    xup.push_back( vertex_up_ECAL.X() );
	    yup.push_back( vertex_up_ECAL.Y() );
	    zup.push_back( vertex_up_ECAL.Z() );

	    wxup.push_back( pow(0.001,-2) );
	    wyup.push_back( pow(0.001,-2) );
	  
	    xup.push_back( xavg_plane[0] );
	    yup.push_back( yavg_plane[0] - yoff_ECAL );
	    zup.push_back( zplane[0] );

	    wxup.push_back( pow( Lx_scint_CDET/2.0, -2 ) );
	    wyup.push_back( pow( sigy_CDET, -2 ) );

	    xup.push_back( xavg_plane[1] );
	    yup.push_back( yavg_plane[1] - yoff_ECAL );
	    zup.push_back( zplane[1] );
	  
	    wxup.push_back( pow( Lx_scint_CDET/2.0, -2 ) );
	    wyup.push_back( pow( sigy_CDET, -2 ) );

	    xup.push_back( xECAL );
	    yup.push_back( yECAL - yoff_ECAL );
	    zup.push_back( ECALdist );

	    wxup.push_back( pow(sigx_ECAL,-2) );
	    wyup.push_back( pow(sigy_ECAL,-2) );

	    double Xup,Yup,dXdZup,dYdZup;

	    Fit_3D_track( xup, yup, zup, wxup, wyup, Xup, Yup, dXdZup, dYdZup ); 

	    TVector3 ehat_up_ECAL( dXdZup, dYdZup, 1.0 );
	    ehat_up_ECAL = ehat_up_ECAL.Unit();

	    TVector3 ehat_up_global =
	      ehat_up_ECAL.X() * ECAL_xaxis +
	      ehat_up_ECAL.Y() * ECAL_yaxis +
	      ehat_up_ECAL.Z() * ECAL_zaxis;

	    double etheta_up = acos( ehat_up_global.Z() );

	    TVector3 vertex_down_ECAL( vertex_down_global.Dot(ECAL_xaxis),
				       vertex_down_global.Dot(ECAL_yaxis),
				       vertex_down_global.Dot(ECAL_zaxis) );
	  
	    xdown.push_back( vertex_down_ECAL.X() );
	    ydown.push_back( vertex_down_ECAL.Y() );
	    zdown.push_back( vertex_down_ECAL.Z() );

	    wxdown.push_back( pow(0.001,-2) );
	    wydown.push_back( pow(0.001,-2) );

	    xdown.push_back( xavg_plane[0] );
	    ydown.push_back( yavg_plane[0] - yoff_ECAL );
	    zdown.push_back( zplane[0] );

	    wxdown.push_back( pow( Lx_scint_CDET/2.0, -2 ) );
	    wydown.push_back( pow(sigy_CDET, -2 ) );

	    xdown.push_back( xavg_plane[1] );
	    ydown.push_back( yavg_plane[1] - yoff_ECAL );
	    zdown.push_back( zplane[1] );

	    wxdown.push_back( pow( Lx_scint_CDET/2.0, -2 ) );
	    wydown.push_back( pow( sigy_CDET, -2 ) );

	    xdown.push_back( xECAL );
	    ydown.push_back( yECAL );
	    zdown.push_back( ECALdist );

	    wxdown.push_back( pow(sigx_ECAL,-2) );
	    wydown.push_back( pow(sigy_ECAL,-2) );

	    double Xdown, Ydown, dXdZdown, dYdZdown;

	    Fit_3D_track( xdown, ydown, zdown, wxdown, wydown,
			  Xdown, Ydown, dXdZdown, dYdZdown );

	    TVector3 ehat_down_ECAL( dXdZdown, dYdZdown, 1.0 );
	    ehat_down_ECAL = ehat_down_ECAL.Unit();

	    TVector3 ehat_down_global =
	      ehat_down_ECAL.X() * ECAL_xaxis +
	      ehat_down_ECAL.Y() * ECAL_yaxis +
	      ehat_down_ECAL.Z() * ECAL_zaxis;

	    double etheta_down = acos( ehat_down_global.Z() );

	    TVector3 vertex_center_ECAL( vertex_center_global.Dot(ECAL_xaxis),
					 vertex_center_global.Dot(ECAL_yaxis),
					 vertex_center_global.Dot(ECAL_zaxis) );

	    xc.push_back( vertex_center_ECAL.X() );
	    yc.push_back( vertex_center_ECAL.Y() );
	    zc.push_back( vertex_center_ECAL.Z() );

	    wxc.push_back( pow(0.001,-2) );
	    wyc.push_back( pow(0.001,-2) );

	    xc.push_back( xavg_plane[0] );
	    yc.push_back( yavg_plane[0] - yoff_ECAL );
	    zc.push_back( zplane[0] );

	    wxc.push_back( pow( Lx_scint_CDET/2.0, -2 ) );
	    wyc.push_back( pow( sigy_CDET, -2 ) );

	    xc.push_back( xavg_plane[1] );
	    yc.push_back( yavg_plane[1] - yoff_ECAL );
	    zc.push_back( zplane[1] );

	    wxc.push_back( pow( Lx_scint_CDET/2.0, -2 ) );
	    wyc.push_back( pow( Ly_scint_CDET/2.0, -2 ) );

	    xc.push_back( xECAL );
	    yc.push_back( yECAL - yoff_ECAL );
	    zc.push_back( ECALdist );

	    wxc.push_back( pow(sigx_ECAL,-2) );
	    wyc.push_back( pow(sigy_ECAL,-2) );

	    double Xc, Yc, dXdZc, dYdZc;
	    Fit_3D_track( xc, yc, zc, wxc, wyc,
			  Xc, Yc, dXdZc, dYdZc );

	    // cout << "(xECAL, yECAL, Xc, Yc, Xpc, Ypc) = (" << xECAL << ", " << yECAL
	    //      << dXdZc << ", " << dYdZc << ")" << endl;
	  
	    TVector3 ehat_c_ECAL( dXdZc, dYdZc, 1.0 );
	    ehat_c_ECAL = ehat_c_ECAL.Unit();
	    TVector3 ehat_c_global =
	      ehat_c_ECAL.X() * ECAL_xaxis +
	      ehat_c_ECAL.Y() * ECAL_yaxis +
	      ehat_c_ECAL.Z() * ECAL_zaxis;

	    double etheta_c = acos( ehat_c_global.Z() );

	    double nu_Eclust = Ebeam - Eprime_ECAL;
	    double Q2_Eclust = 2.*Mp*nu_Eclust;

	    double pp_Eclust = sqrt( pow(nu_Eclust,2) + 2.*Mp*nu_Eclust);
	  
	    double etheta_Eclust = acos(1. - Q2_Eclust/(2.*Ebeam*Eprime_ECAL));
	  
	    // cout << "etheta c, etheta_up, etheta_down, etheta(E'), etheta true = " << etheta_c * 180.0/PI << " deg, "
	    //      << etheta_up * 180.0/PI << " deg, "
	    //      << etheta_down * 180.0/PI << " deg, "
	    //      << etheta_Eclust * 180.0/PI << " deg, "
	    //      << T->ev_th * 180.0/PI << " deg." << endl;
				 
	    double ephi_c = atan2( ehat_c_global.Y(), ehat_c_global.X() );

	    //cout << "ephi_c, ephi true = " << ephi_c * 180.0/PI << ", " << T->ev_ph * 180.0/PI << endl;

	    hdethc_Eclust->Fill( etheta_c - etheta_Eclust, weight  );
	    hdethc_thtrue->Fill( etheta_c - T->ev_th, weight  );
	    hdephc_phtrue->Fill( ephi_c - T->ev_ph, weight  );
	    hdethc_thtrue_vs_z->Fill( T->ev_vz, etheta_c - T->ev_th, weight  );

	    //Calculate proton kinematics from reconstructed electron:

	    double Eprime_etheta_c = Ebeam / (1.0 + Ebeam/Mp*(1.0-cos(etheta_c)));
	    double nu_etheta_c = Ebeam - Eprime_etheta_c;

	    double pp_etheta_c = sqrt(pow(nu_etheta_c,2) + 2.*Mp*nu_etheta_c );

	    double pphi_ephi_c = ephi_c + PI;

	    double ptheta_etheta_c = acos( (Ebeam - Eprime_etheta_c * cos(etheta_c))/pp_etheta_c );

	    TVector3 nhat_p_eth_global( sin(ptheta_etheta_c)*cos(pphi_ephi_c),
					sin(ptheta_etheta_c)*sin(pphi_ephi_c),
					cos(ptheta_etheta_c) );

	    TVector3 nhat_p_eth_SBS( nhat_p_eth_global.Dot(SBS_xaxis_global),
				     nhat_p_eth_global.Dot(SBS_yaxis_global),
				     nhat_p_eth_global.Dot(SBS_zaxis_global) );

	    double xptar_eth = nhat_p_eth_SBS.X()/nhat_p_eth_SBS.Z();
	    double yptar_eth = nhat_p_eth_SBS.Y()/nhat_p_eth_SBS.Z();
	    double ytar_eth = 0.0;
	    double xtar_eth = 0.0;

	    double xfp_eth,yfp_eth,xpfp_eth,ypfp_eth;

	    SBS_fp_reconstruct( 0.0, 0.0, xptar_eth, yptar_eth, pp_etheta_c,
				xfp_eth, yfp_eth, xpfp_eth, ypfp_eth );
	  
	  
	    hdxfp->Fill( xfp_eth - xfptemp, weight  );
	    hdyfp->Fill( yfp_eth - yfptemp, weight  );
	    hdxpfp->Fill( xpfp_eth - xpfptemp, weight  );
	    hdypfp->Fill( ypfp_eth - ypfptemp, weight  );

	    hdxfp_vs_ztrue->Fill( T->ev_vz, xfp_eth - xfptemp, weight  );
	    hdyfp_vs_ztrue->Fill( T->ev_vz, yfp_eth - yfptemp, weight  );
	    hdxpfp_vs_ztrue->Fill( T->ev_vz, xpfp_eth - xpfptemp, weight  );
	    hdypfp_vs_ztrue->Fill( T->ev_vz, ypfp_eth - ypfptemp, weight  );

	    hdpp_etheta_c->Fill( pp_etheta_c/T->ev_np - 1., weight  );
	    hdpp_eclust->Fill( pp_Eclust/T->ev_np - 1., weight  );

	    xz.push_back( vertex_true_ECAL.X() );
	    yz.push_back( vertex_true_ECAL.Y() );
	    zz.push_back( vertex_true_ECAL.Z() );

	    wxz.push_back( pow(0.001,-2) );
	    wyz.push_back( pow(0.001,-2) );

	    xz.push_back( xavg_plane[0] );
	    yz.push_back( yavg_plane[0] - yoff_ECAL );
	    zz.push_back( zplane[0] );

	    wxz.push_back( pow( Lx_scint_CDET/2.0, -2 ) );
	    wyz.push_back( pow( sigy_CDET, -2 ) );

	    xz.push_back( xavg_plane[1] );
	    yz.push_back( yavg_plane[1] - yoff_ECAL );
	    zz.push_back( zplane[1] );

	    wxz.push_back( pow( Lx_scint_CDET/2.0, -2 ) );
	    wyz.push_back( pow( sigy_CDET, -2 ) );

	    xz.push_back( xECAL );
	    yz.push_back( yECAL - yoff_ECAL );
	    zz.push_back( ECALdist );

	    wxz.push_back( pow( sigx_ECAL, -2 ) );
	    wyz.push_back( pow( sigy_ECAL, -2 ) );

	    double Xz, Yz, Xpz, Ypz;

	    Fit_3D_track( xz, yz, zz, wxz, wyz,
			  Xz, Yz, Xpz, Ypz );

	    TVector3 ehat_ztrue_ECAL( Xpz, Ypz, 1.0 );

	    ehat_ztrue_ECAL = ehat_ztrue_ECAL.Unit();

	    TVector3 ehat_ztrue_global =
	      ehat_ztrue_ECAL.X() * ECAL_xaxis +
	      ehat_ztrue_ECAL.Y() * ECAL_yaxis +
	      ehat_ztrue_ECAL.Z() * ECAL_zaxis;

	    //	  TVector3 vtx_ztrue_ECAL( Xz, Yz, 
	    double etheta_ztrue = acos( ehat_ztrue_global.Z() );
	    double ephi_ztrue = atan2( ehat_ztrue_global.Y(), ehat_ztrue_global.X() );
	  
	    double Eprime_etheta_zt = Ebeam/(1.0 + Ebeam/Mp*(1.-cos(etheta_ztrue)));

	    double nu_etheta_zt = Ebeam - Eprime_etheta_zt;
	    double pp_etheta_zt = sqrt(pow(nu_etheta_zt,2)+2.*Mp*nu_etheta_zt );

	    double ptheta_etheta_zt = acos( (Ebeam-Eprime_etheta_zt*cos(etheta_ztrue))/
					    pp_etheta_zt );
	  
	    double pphi_ephi_zt = ephi_ztrue + PI;

	    hdpp_etheta_vtx_known->Fill( pp_etheta_zt/T->ev_np - 1.0, weight  );

	    hdpp_etheta_vtx_known->GetXaxis()->SetTitle("p_{p}(e arm)/p_{p}^{true}-1 (v_{z} = v_{z}(true))");
	  
	    TVector3 vertex_true_SBS( vertex_true_global.Dot( SBS_xaxis_global ),
				      vertex_true_global.Dot( SBS_yaxis_global ),
				      vertex_true_global.Dot( SBS_zaxis_global ) );

	    TVector3 phat_eth_zt_global( sin(ptheta_etheta_zt)*cos(pphi_ephi_zt),
					 sin(ptheta_etheta_zt)*sin(pphi_ephi_zt),
					 cos(ptheta_etheta_zt) );

	    TVector3 phat_eth_zt_SBS( phat_eth_zt_global.Dot(SBS_xaxis_global),
				      phat_eth_zt_global.Dot(SBS_yaxis_global),
				      phat_eth_zt_global.Dot(SBS_zaxis_global) );

	    double xptar_eth_zt = phat_eth_zt_SBS.X()/phat_eth_zt_SBS.Z();
	    double yptar_eth_zt = phat_eth_zt_SBS.Y()/phat_eth_zt_SBS.Z();
	    double ytar_eth_zt = vertex_true_SBS.Y() - vertex_true_SBS.Z()*yptar_eth_zt;
	    double xtar_eth_zt = vertex_true_SBS.X() - vertex_true_SBS.Z()*xptar_eth_zt;

	    //Now reconstruct SBS fp quantities from electron arm quantities assuming vertex position is known:

	    double xfp_eth_zt, yfp_eth_zt, xpfp_eth_zt, ypfp_eth_zt;
	    SBS_fp_reconstruct( xtar_eth_zt, ytar_eth_zt, xptar_eth_zt, yptar_eth_zt, pp_etheta_zt,
				xfp_eth_zt, yfp_eth_zt, xpfp_eth_zt, ypfp_eth_zt );

	    hdxfp_vtx_known->Fill( xfp_eth_zt - xfptemp, weight  );
	    hdyfp_vtx_known->Fill( yfp_eth_zt - yfptemp, weight  );
	    hdxpfp_vtx_known->Fill( xpfp_eth_zt - xpfptemp, weight  );
	    hdypfp_vtx_known->Fill( ypfp_eth_zt - ypfptemp, weight  );
			     

	    TVector3 vertex_recon_global(T->ev_vx, T->ev_vy, vzrecon );

	    TVector3 vertex_recon_ECAL( vertex_recon_global.Dot(ECAL_xaxis),
					vertex_recon_global.Dot(ECAL_yaxis),
					vertex_recon_global.Dot(ECAL_zaxis) );

	    //calculate electron kinematics using ONLY measured proton scattering angles and beam energy
	    //no momentum!
	    double nu_ptheta_recon = sqrt(pow(pp_ptheta_recon,2)+pow(Mp,2))-Mp;
	    double Eprime_ptheta_recon = Ebeam - nu_ptheta_recon;
	    double Q2_ptheta_recon = 2.*Mp*nu_ptheta_recon;
	    double etheta_ptheta_recon = acos( 1.0 - Q2_ptheta_recon/(2.*Ebeam*Eprime_ptheta_recon) );
	  
	    double ephi_pphi_recon = pphirecon - PI;

	    TVector3 ehat_ptheta_pphi( sin(etheta_ptheta_recon)*cos(ephi_pphi_recon),
				       sin(etheta_ptheta_recon)*sin(ephi_pphi_recon),
				       cos(etheta_ptheta_recon) );

	    TVector3 ehat_pth_pph_ECAL( ehat_ptheta_pphi.Dot(ECAL_xaxis),
					ehat_ptheta_pphi.Dot(ECAL_yaxis),
					ehat_ptheta_pphi.Dot(ECAL_zaxis) );

	    double xclust_SBS = vertex_recon_ECAL.X() + ehat_pth_pph_ECAL.X()/ehat_pth_pph_ECAL.Z() * ( ECALdist - vertex_recon_ECAL.Z() );
	    double yclust_SBS = vertex_recon_ECAL.Y() + ehat_pth_pph_ECAL.Y()/ehat_pth_pph_ECAL.Z() * ( ECALdist - vertex_recon_ECAL.Z() );

	    hdeltax_pth->Fill( 100.*(xECAL - xclust_SBS), weight );
	    hdeltay_pth->Fill( 100.*(yECAL - yoff_ECAL - yclust_SBS), weight );
	    hdxdy_pth->Fill( 100.*(xECAL - xclust_SBS), 100.*(yECAL-yoff_ECAL - yclust_SBS), weight );

	    double nu_pp_recon = sqrt(pow(pprecon,2)+pow(Mp,2))-Mp;
	    double Eprime_pp_recon = Ebeam - nu_pp_recon;
	    double Q2_pp_recon = 2.*Mp*nu_pp_recon;

	    double etheta_pp_recon = acos( 1.0 - Q2_pp_recon/(2.*Ebeam*Eprime_pp_recon));

	    //predict electron trajectory using pp, pphi:
	    TVector3 ehat_pp_pphi( sin(etheta_pp_recon)*cos(ephi_pphi_recon),
				   sin(etheta_pp_recon)*sin(ephi_pphi_recon),
				   cos(etheta_pp_recon) );

	    TVector3 ehat_pp_pphi_ECAL( ehat_pp_pphi.Dot(ECAL_xaxis),
					ehat_pp_pphi.Dot(ECAL_yaxis),
					ehat_pp_pphi.Dot(ECAL_zaxis) );

	    double xclust_pp = vertex_recon_ECAL.X() + ehat_pp_pphi_ECAL.X()/ehat_pp_pphi_ECAL.Z() * ( ECALdist - vertex_recon_ECAL.Z() );
	    double yclust_pp = vertex_recon_ECAL.Y() + ehat_pp_pphi_ECAL.Y()/ehat_pp_pphi_ECAL.Z() * (ECALdist - vertex_recon_ECAL.Z() );

	    hdeltax_pp->Fill( 100.*(xECAL - xclust_pp ), weight  );
	    hdeltay_pp->Fill( 100.*(yECAL-yoff_ECAL - yclust_pp ), weight  );
	    hdxdy_pp->Fill( 100.*(xECAL - xclust_pp ), 100.*(yECAL-yoff_ECAL - yclust_pp ), weight  );

	    //hdphip->Fill( pphirecon
	  
	    //NOW we can do final electron scattering angle reconstruction using SBS reconstructed vertex:

	    vector<double> xfinal,yfinal,zfinal,wxfinal,wyfinal;

	    xfinal.push_back( vertex_recon_ECAL.X() );
	    yfinal.push_back( vertex_recon_ECAL.Y() );
	    zfinal.push_back( vertex_recon_ECAL.Z() );

	    TVector3 sigvtx_global(0.001,0.001,0.007);

	    TVector3 sigvtx_ECAL( sigvtx_global.Dot( ECAL_xaxis ),
				  sigvtx_global.Dot( ECAL_yaxis ),
				  sigvtx_global.Dot( ECAL_zaxis ) );

	    wxfinal.push_back( pow( sigvtx_ECAL.X(), -2 ) );
	    wyfinal.push_back( pow( sigvtx_ECAL.Y(), -2 ) );

	    xfinal.push_back( xavg_plane[0] );
	    yfinal.push_back( yavg_plane[0] - yoff_ECAL );
	    zfinal.push_back( zplane[0] );

	    wxfinal.push_back( pow( Lx_scint_CDET, -2 ) );
	    wyfinal.push_back( pow( sigy_CDET, -2 ) );

	    xfinal.push_back( xavg_plane[1] );
	    yfinal.push_back( yavg_plane[1] - yoff_ECAL );
	    zfinal.push_back( zplane[1] );

	    wxfinal.push_back( pow( Lx_scint_CDET, -2 ) );
	    wyfinal.push_back( pow( sigy_CDET, -2 ) );

	    xfinal.push_back( xECAL );
	    yfinal.push_back( yECAL - yoff_ECAL );
	    zfinal.push_back( ECALdist );
			     
	    wxfinal.push_back( pow( sigx_ECAL, -2 ) );
	    wyfinal.push_back( pow( sigy_ECAL, -2 ) );

	    double eX, eY, eXp, eYp;

	    Fit_3D_track( xfinal, yfinal, zfinal, wxfinal, wyfinal,
			  eX, eY, eXp, eYp );

	    TVector3 ehat_final_ECAL( eXp, eYp, 1.0 );
	    ehat_final_ECAL = ehat_final_ECAL.Unit();

	    TVector3 ehat_final_global =
	      ehat_final_ECAL.X() * ECAL_xaxis +
	      ehat_final_ECAL.Y() * ECAL_yaxis +
	      ehat_final_ECAL.Z() * ECAL_zaxis;

	    double etheta_recon = acos( ehat_final_global.Z() );
	    double ephi_recon = atan2( ehat_final_global.Y(), ehat_final_global.X() );

	    double Eprime_etheta_recon = Ebeam / (1.0 + Ebeam/Mp*(1.-cos(etheta_recon)));
	    double nu_etheta_recon = Ebeam - Eprime_etheta_recon;

	    double pp_etheta_recon = sqrt(pow(nu_etheta_recon,2) + 2.*Mp*nu_etheta_recon);

	    double ptheta_etheta_recon = acos( (Ebeam-Eprime_etheta_recon*cos(etheta_recon))/pp_etheta_recon);
	  
	    hdphip->Fill( pphirecon - ephi_recon - PI, weight  );
	    hdthetap->Fill( pthetarecon - ptheta_etheta_recon, weight  );
	    hpmisse->Fill( (pprecon - pp_etheta_recon)/ptrue, weight  );
	    hdpp_eth_pth->Fill( (pp_ptheta_recon - pp_etheta_recon)/ptrue, weight  );
	    hdthetae_true->Fill( etheta_recon - T->ev_th, weight  );
	    hdphie_true->Fill( ephi_recon - T->ev_ph, weight );

	    //Evaluate cuts: first, 1-Eclust/E'(thetae): apply all cuts except deltaE:
	    bool passed_dx = fabs( 100.*(xECAL - xclust_SBS) - 0.5*(dxmin+dxmax) ) <= 0.5*(dxmax-dxmin);
	    bool passed_dy = fabs( 100.*(yECAL - yoff_ECAL - yclust_SBS) - 0.5*(dxmax+dxmin) ) <= 0.5*(dymax-dymin);
	    bool passed_dE = fabs( 1.0 - Eprime_ECAL/Eprime_etheta_temp - 0.5*(dEmax+dEmin) ) <= 0.5*(dEmax-dEmin);
	    bool passed_dpp = fabs( (pprecon-pp_ptheta_recon)/ptrue - 0.5*(dppmax+dppmin) ) <= 0.5*(dppmax-dppmin);

	    nevents_all += 1.0;
	    sum_weights_all += weight;
	    sum2_weights_all += pow(weight,2);
	    if( passed_dx && passed_dy && passed_dE && passed_dpp ){
	      sum_weights_cut += weight;
	      sum2_weights_cut += pow(weight,2);
	      nevents_cut += 1.0;
	    }
	    
	    if( passed_dx && passed_dy && passed_dpp ){
	      hdEclust_Eprime_eth_cut->Fill( 1.0 - Eprime_ECAL/Eprime_etheta_temp, weight );
	    } else {
	      hdEclust_Eprime_eth_anticut->Fill( 1.0 - Eprime_ECAL/Eprime_etheta_temp, weight );
	    }

	    if( passed_dy && passed_dpp && passed_dE ){
	      hdeltax_pth_cut->Fill( 100.*(xECAL - xclust_SBS), weight ); 
	    } else {
	      hdeltax_pth_anticut->Fill( 100.*(xECAL - xclust_SBS), weight ); 
	    }

	    if( passed_dx && passed_dpp && passed_dE ){
	      hdeltay_pth_cut->Fill( 100.*(yECAL - yclust_SBS), weight ); 
	    } else {
	      hdeltay_pth_anticut->Fill( 100.*(yECAL - yclust_SBS), weight ); 
	    }

	    if( passed_dx && passed_dy && passed_dE ){
	      hpmissp_fractional_cut->Fill( (pprecon-pp_ptheta_recon)/ptrue, weight );
	    } else {
	      hpmissp_fractional_anticut->Fill( (pprecon-pp_ptheta_recon)/ptrue, weight );
	    }
	  }
	}
      }
    }
    
  }

  hnclust->GetXaxis()->SetTitle("Number of ECAL clusters found");
  hxyclust->GetXaxis()->SetTitle("Cluster x (m), before corrections");
  hxyclust->GetYaxis()->SetTitle("Cluster y (m), before corrections");
  hEclust->GetXaxis()->SetTitle("Cluster energy (GeV), uncorrected");
  hEclust_vs_Etrue->GetXaxis()->SetTitle("Incident electron energy E'_{e} (GeV)");
  hEclust_vs_Etrue->GetYaxis()->SetTitle("Cluster energy (GeV), uncorrected");
  hdEclust_fractional->GetXaxis()->SetTitle("1-E_{clust}/E'_{e}");
  hdEclust_fractional_vs_Eclust->GetXaxis()->SetTitle("E'_{e} (GeV)");
  hdEclust_fractional_vs_Eclust->GetYaxis()->SetTitle("1-E_{clust}/E'_{e}");
  hnhitclust->GetXaxis()->SetTitle("Number of hits per cluster");
  hdxclust->GetXaxis()->SetTitle("x_{clust} - x_{true} (cm), before corrections");
  hdyclust->GetXaxis()->SetTitle("y_{clust} - y_{true} (cm), before corrections");

  hxmom->GetXaxis()->SetTitle("(#bar{x}-x_{max})/L");
  hymom->GetXaxis()->SetTitle("(#bar{y}-y_{max})/L");

  hdxx->GetXaxis()->SetTitle("x_{true} (cm)");
  hdxx->GetYaxis()->SetTitle("x_{clust}-x_{true} (cm), before corrections");
  hdyy->GetXaxis()->SetTitle("y_{true} (cm)");
  hdyy->GetYaxis()->SetTitle("y_{clust}-y_{true} (cm), before corrections");

  hxmom_x->GetYaxis()->SetTitle("(#bar{x}-x_{max})/L");
  hxmom_x->GetXaxis()->SetTitle("x_{max} (cm)");

  hymom_y->GetYaxis()->SetTitle("(#bar{y}-y_{max})/L");
  hymom_y->GetXaxis()->SetTitle("y_{max} (cm)");

  hxyclust_corrected->GetXaxis()->SetTitle("Cluster x (m), corrected");
  hxyclust_corrected->GetYaxis()->SetTitle("Cluster y (m), corrected");
  
  hdxclust_corrected->GetXaxis()->SetTitle("x_{clust} - x_{true} (cm), after corrections");
  hdyclust_corrected->GetXaxis()->SetTitle("y_{clust} - y_{true} (cm), after corrections");
  hdxx_corrected->GetXaxis()->SetTitle("x_{true} (cm)");
  hdxx_corrected->GetYaxis()->SetTitle("x_{clust} - x_{true} (cm), corrected");
  hdyy_corrected->GetXaxis()->SetTitle("y_{true} (cm)");
  hdyy_corrected->GetYaxis()->SetTitle("y_{clust} - y_{true} (cm), corrected");

  hxcorr_xmax->GetXaxis()->SetTitle("x_{clust} - x_{max} (m), after corrections");
  hycorr_ymax->GetXaxis()->SetTitle("y_{clust} - y_{max} (m), after corrections");

  hnplanes_CDET->GetXaxis()->SetTitle("Number of CDET planes fired within cluster tolerance");
  hnhits1_CDET->GetXaxis()->SetTitle("Number of hits in CDET plane 1");
  hnhits2_CDET->GetXaxis()->SetTitle("Number of hits in CDET plane 2");

  hdy1_CDET_ECAL->GetXaxis()->SetTitle("y_{1}(CDET) - y(ECAL) (m)");
  hdy2_CDET_ECAL->GetXaxis()->SetTitle("y_{2}(CDET) - y(ECAL) (m)");
  hdy12_CDET->GetXaxis()->SetTitle("y_{1}(CDET) - y_{2} (CDET) (m)");
  hEtot1_CDET->GetXaxis()->SetTitle("CDET plane 1 energy deposit (GeV)");
  hEtot2_CDET->GetXaxis()->SetTitle("CDET plane 2 energy deposit (GeV)");

  hdEclust_Eprime_eth->GetXaxis()->SetTitle("1-E_{clust}(corrected)/E'(#theta_{e})");
  hdethc_Eclust->GetXaxis()->SetTitle("#theta_{e} - #theta_{e}(E_{clust}) (rad), v_{z} = 0 assumed");
  hdethc_thtrue->GetXaxis()->SetTitle("#theta_{e}(recon) - #theta_{e}(true) (rad), v_{z} = 0 assumed");
  hdethc_thtrue_vs_z->GetXaxis()->SetTitle("true vertex z (m)");
  hdethc_thtrue_vs_z->GetYaxis()->SetTitle("#theta_{e} (recon) - #theta_{e} (true) (rad), v_{z} = 0 assumed");

  hdephc_phtrue->GetXaxis()->SetTitle("#phi_{e} (recon) - #phi_{e} (true) (rad), v_{z} = 0 assumed");

  hdxfp->GetXaxis()->SetTitle("x_{fp} (e arm) - x_{fp} (true) (m), v_{z} = 0 assumed");
  hdyfp->GetXaxis()->SetTitle("y_{fp} (e arm) - y_{fp} (true) (m), v_{z} = 0 assumed");
  hdxpfp->GetXaxis()->SetTitle("x'_{fp} (e arm) - x'_{fp} (true), v_{z} = 0 assumed");
  hdypfp->GetXaxis()->SetTitle("y'_{fp} (e arm) - y'_{fp} (true), v_{z} = 0 assumed");

  hdxfp_vs_ztrue->GetXaxis()->SetTitle("v_{z}^{true} (m)");
  hdxfp_vs_ztrue->GetYaxis()->SetTitle("x_{fp} (e arm) - x_{fp} (true) (m), v_{z}=0 assumed");

  hdyfp_vs_ztrue->GetXaxis()->SetTitle("v_{z}^{true} (m)");
  hdyfp_vs_ztrue->GetYaxis()->SetTitle("y_{fp} (e arm) - y_{fp} (true) (m), v_{z}=0 assumed");

  hdxpfp_vs_ztrue->GetXaxis()->SetTitle("v_{z}^{true} (m)");
  hdxpfp_vs_ztrue->GetYaxis()->SetTitle("x'_{fp} (e arm) - x'_{fp} (true), v_{z}=0 assumed");

  hdypfp_vs_ztrue->GetXaxis()->SetTitle("v_{z}^{true} (m)");
  hdypfp_vs_ztrue->GetYaxis()->SetTitle("y'_{fp} (e arm) - y'_{fp} (true), v_{z}=0 assumed");
  
  
  
  hdpp_etheta_c->GetXaxis()->SetTitle("p_{p}(e arm)/p_{p}(true)-1, v_{z} = 0 assumed");
  hdpp_eclust->GetXaxis()->SetTitle("p_{p}(E_{clust})/p_{p}(true)-1");

  hdxfp_vtx_known->GetXaxis()->SetTitle("x_{fp} (e arm) - x_{fp} (true) (m), v_{z} = v_{z}(true)");
  hdyfp_vtx_known->GetXaxis()->SetTitle("y_{fp} (e arm) - y_{fp} (true) (m), v_{z} = v_{z}(true)");
  hdxpfp_vtx_known->GetXaxis()->SetTitle("x'_{fp} (e arm) - x'_{fp} (true), v_{z} = v_{z}(true)");
  hdypfp_vtx_known->GetXaxis()->SetTitle("y'_{fp} (e arm) - y'_{fp} (true), v_{z} = v_{z}(true)");
  
  hdeltax_pth->GetXaxis()->SetTitle("x_{clust} - x_{clust}(#theta_{p}) (cm), v_{z} = v_{z}^{true}");
  hdeltay_pth->GetXaxis()->SetTitle("y_{clust} - y_{clust}(#theta_{p}) (cm), v_{z} = v_{z}^{true}");
  hdxdy_pth->GetXaxis()->SetTitle("x_{clust} - x_{clust}(#theta_{p}) (cm), v_{z} = v_{z}^{true}");
  hdxdy_pth->GetYaxis()->SetTitle("y_{clust} - y_{clust}(#theta_{p}) (cm), v_{z} = v_{z}^{true}");

  hdeltax_pp->GetXaxis()->SetTitle("x_{clust} - x_{clust}(p_{p}) (cm), v_{z} = v_{z}^{true}");
  hdeltay_pp->GetXaxis()->SetTitle("y_{clust} - y_{clust}(p_{p}) (cm), v_{z} = v_{z}^{true}");
  hdxdy_pp->GetXaxis()->SetTitle("x_{clust} - x_{clust}(p_{p}) (cm), v_{z} = v_{z}^{true}");
  hdxdy_pp->GetYaxis()->SetTitle("y_{clust} - y_{clust}(p_{p}) (cm), v_{z} = v_{z}^{true}");
  
  hdphip->GetXaxis()->SetTitle("#phi_{p}(recon) - #phi_{e} (recon)-#pi (rad)");
  hdthetap->GetXaxis()->SetTitle("#theta_{p}(recon) - #theta_{p}(#theta_{e}) (rad)");
  hpmisse->GetXaxis()->SetTitle("(p_{p}(recon) - p_{el}(#theta_{e}))/p_{p}^{true}");
  hdpp_eth_pth->GetXaxis()->SetTitle( "(p_{p}(#theta_{p})-p_{p}(#theta_{e}))/p_{p}^{true}");
  hdthetae_true->GetXaxis()->SetTitle( "#theta_{e}(recon) - #theta_{e}(true) (rad)");
  hdphie_true->GetXaxis()->SetTitle("#phi_{e}(recon) - #phi_{e}(true) (rad)" );

  double dsum_weights_all = sum_weights_all/sqrt(nevents_all);
  double dsum_weights_cut = sum_weights_cut/sqrt(nevents_cut);
  
  cout << "Total event rate = " << sum_weights_all << "+/-" << dsum_weights_all << " Hz" << endl;
  cout << "Total event rate passing all four exclusivity cuts = " << sum_weights_cut << "+/-" << dsum_weights_cut << " Hz" << endl;

  elist->Delete();
  
  fout->Write();
}
