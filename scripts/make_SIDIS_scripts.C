#include <iostream>
#include <fstream>
#include "TString.h"
#include "TMath.h"
#include "TVector3.h"
#include "TObjArray.h"
#include "TObjString.h"

void make_SIDIS_scripts( const char *inputfilename ){

  double PI = TMath::Pi();

  double Mp = 0.938272;

  int itgt;
  int nevents=100000;

  int fieldclamp=0;
  int ckov_flag = 1;

  double Ibeam, Ltgt, Ptgt;
  double xrast, yrast;
  double Ebeam, thbb, Dbb, thsbs, Dsbs, Dhcal, Drich;
  double SBS_magfield;
  double pmin_SBS, pmin_BB;

  ifstream inputfile(inputfilename);

  inputfile >> nevents;
  inputfile >> itgt; //0 = LH2, 1 = LD2, 2 = H2 gas, 3 = 3He gas:
  inputfile >> fieldclamp;
  inputfile >> ckov_flag;
  inputfile >> Ibeam >> Ltgt >> Ptgt;
  inputfile >> xrast >> yrast;
  inputfile >> Ebeam >> thbb >> Dbb >> thsbs >> Dsbs;
  inputfile >> SBS_magfield;
  inputfile >> pmin_SBS >> pmin_BB;
  
  Dhcal = Dsbs + 4.0;
  Drich = Dsbs + 2.5;

  //The files we want to create are: pi+/-, K+/-, pi0, with upbending and downbending configurations for the charged species.
  // This means 9 files.

  //Now, let's compute angle generation limits: we know that at D = 1.5 meters, xptar acceptance is ~ +/- 0.4, yptar acceptance
  // ~+/-0.2 (conservatively). 
  double dxptar_BB = 0.35 * 1.5/Dbb;
  double dyptar_BB = 0.15 * 1.5/Dbb;
  
  double BBthetamin = thbb - atan(dyptar_BB)*180.0/PI;
  
  double pxhat_spec = dxptar_BB/sqrt(1.0+pow(dxptar_BB,2) + pow(dyptar_BB,2) );
  double pyhat_spec = -dyptar_BB/sqrt(1.0+pow(dxptar_BB,2) + pow(dyptar_BB,2) );
  double pzhat_spec = 1.0/sqrt(1.0+pow(dxptar_BB,2) + pow(dyptar_BB,2) );

  double pxhat_lab = pxhat_spec;
  double pyhat_lab = pyhat_spec*cos(thbb*PI/180.0) - pzhat_spec * sin(thbb*PI/180.0);
  double pzhat_lab = pyhat_spec*sin(thbb*PI/180.0) + pzhat_spec * cos(thbb*PI/180.0);
  
  //  cout << "pzhat_spec, pzhat_lab = " << pzhat_spec << ", " << pzhat_lab << endl;

  double BBthetamax = 180.0/PI * acos(pzhat_lab);
  
  // cout << "BigBite thetamax = " << BBthetamax << endl;

  //Phi max/min computed from smallest in-plane angle, largest out-of-plane
  pyhat_spec = dyptar_BB/sqrt(1.0+pow(dxptar_BB,2)+pow(dyptar_BB,2));
  pxhat_spec = -dxptar_BB/sqrt(1.0+pow(dxptar_BB,2)+pow(dyptar_BB,2));

  pyhat_lab = pyhat_spec*cos(thbb*PI/180.0) - pzhat_spec * sin(thbb*PI/180.0);

  double BBphimin = 180.0/PI * atan2( -pxhat_spec, pyhat_lab );
  double BBphimax = 360.0 + 180.0/PI * atan2( pxhat_spec, pyhat_lab );
  
  double dxptar_SBS = 0.3 * 2.45/Dsbs;
  double dyptar_SBS = 0.1 * 2.45/Dsbs;

  double SBSthetamin = thsbs - atan(dyptar_SBS)*180.0/PI;
  
  pxhat_spec = dxptar_SBS/sqrt(1.0+pow(dxptar_SBS,2)+pow(dyptar_SBS,2));
  pyhat_spec = dyptar_SBS/sqrt(1.0+pow(dxptar_SBS,2)+pow(dyptar_SBS,2));
  pzhat_spec = 1.0/sqrt(1.0+pow(dxptar_SBS,2)+pow(dyptar_SBS,2));

  pxhat_lab = pxhat_spec;
  pyhat_lab = pyhat_spec*cos(thsbs*PI/180.0) + pzhat_spec*sin(thsbs*PI/180.0);
  pzhat_lab = -pyhat_spec*sin(thsbs*PI/180.0) + pzhat_spec*cos(thsbs*PI/180.0);
  
  double SBSthetamax = 180.0/PI * acos( pzhat_lab );

  pyhat_spec = -dyptar_SBS/sqrt(1.0+pow(dxptar_SBS,2)+pow(dyptar_SBS,2));
  pxhat_spec = dxptar_SBS/sqrt(1.0+pow(dxptar_SBS,2)+pow(dyptar_SBS,2));
  
  pyhat_lab = pyhat_spec*cos(thsbs*PI/180.0) + pzhat_spec*sin(thsbs*PI/180.0);
  
  double SBSphimin = 180.0/PI*atan2( -pxhat_spec, pyhat_lab );
  double SBSphimax = 180.0/PI*atan2( pxhat_spec, pyhat_lab );

  TString outfilename = "sidis_";

  TString hadrons[5] = {"km_", "pim_", "pi0_", "pip_", "kp_"};
  TString hadroncmd[5] = {"K-", "pi-", "pi0", "pi+", "K+" };
  TString targets[4] = {"LH2", "LD2", "H2", "3He" };

  for(int ihadron=-2; ihadron<=2; ihadron++){
    
    for(int sbspol=-1; sbspol<=1; sbspol += 2 ){

      outfilename = "sidis_";
      outfilename += hadrons[ihadron+2];

      if( sbspol*ihadron > 0 ){
	outfilename += "upbend_";
      } else if( sbspol*ihadron < 0 ){
	outfilename += "downbend_";
      }

      char cthD[80];
      sprintf(cthD,"_ebeam%gGeV_thbb%gdeg_thsbs%gdeg_Dbb%gm_Dsbs%gm",Ebeam, thbb, thsbs, Dbb, Dsbs );
      
      outfilename += targets[itgt];
      outfilename += cthD;
      outfilename += ".mac";
      ofstream outputfile(outfilename.Data());

      ifstream template_macro("sidis_template.mac"); 
      TString currentline;
      while( currentline.ReadLine(template_macro) ){
	TObjArray *tokens = currentline.Tokenize(" ");
	TString cmd_prefix = ( (TObjString*) (*tokens)[0] )->GetString();
	// cout << "read line, ntokens = " << tokens->GetEntries() << ": ";
	// for(int i=0; i<tokens->GetEntries(); i++){
	//   TString token = ( (TObjString*) (*tokens)[i] )->GetString();
	//   cout << token.Data() << ", ";
	// }
	// cout << endl;
	
	int ntokens = tokens->GetEntries();
	
	TString value=" ", unit=" ";
	TString cmd_temp;

	if( cmd_prefix == "/g4sbs/beamcur" ){
	  //value = " ";
	  //value += Ibeam;
	  value.Form(" %5.3g",Ibeam);
	  unit = " muA";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/target" ) {
	  value = " ";
	  value += targets[itgt]; 
	  //value.Form( " %s",targets[itgt] );
	  unit = "";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	  
	} else if( cmd_prefix == "/g4sbs/targpres" ){
	  //value = " ";
	  //value += Ptgt;
	  value.Form(" %5.3g",Ptgt);
	  unit = " atmosphere";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	  
	} else if( cmd_prefix == "/g4sbs/targlen" ){
	  //value = " ";
	  //value += Ltgt;
	  value.Form( " %5.3g", Ltgt);
	  unit = " cm";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	  
	} else if( cmd_prefix == "/g4sbs/rasterx" ){
	  //	  value = " ";
	  //value += xrast;
	  value.Form(" %5.3g", xrast);
	  unit = " mm";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/rastery" ){
	  //value = " ";
	  //value += yrast;
	  value.Form(" %5.3g", yrast);
	  unit = " mm";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/hadron" ){
	  value = " ";
	  value += hadroncmd[ihadron+2];
	  //value.Form(" %s",hadroncmd[ihadron+2]);
	  unit = "";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/beamE" ){
	  //value = " ";
	  //value += Ebeam;
	  value.Form(" %5.3g",Ebeam);
	  unit = " GeV";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/bbang" ){
	  //value = " ";
	  //value += thbb;
	  value.Form(" %5.3g", thbb);
	  unit = " deg";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/bbdist" ){
	  //value = " ";
	  //value += Dbb;
	  value.Form(" %5.3g", Dbb);
	  unit = " m";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/hcalang"){
	  //value = " ";
	  //value += thsbs;
	  value.Form(" %5.3g", thsbs);
	  unit = " deg";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/hcaldist"){
	  //value = " ";
	  //value += Dhcal;
	  value.Form(" %5.3g",Dhcal);
	  unit = " m";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/48D48dist" ){
	  //value = " ";
	  //value += Dsbs;
	  value.Form(" %5.3g", Dsbs);
	  unit = " m";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/sbsmagfield" ){
	  //value = " ";
	  //value += SBS_magfield*sbspol;
	  value.Form(" %5.3g", SBS_magfield*sbspol);
	  unit = " tesla";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/sbsclampopt" ){
	  //value = " ";
	  //value += fieldclamp;
	  value.Form(" %d",fieldclamp);
	  unit = "";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/richdist" ){
	  //value = " ";
	  //value += Drich;
	  value.Form(" %5.3g", Drich);
	  unit = " m";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/thmin" ){
	  //value = " ";
	  //value += BBthetamin;
	  value.Form(" %5.3g", BBthetamin);
	  unit = " deg";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/phmin" ){
	  //value = " ";
	  //value += BBphimin;
	  value.Form(" %5.3g",BBphimin);
	  unit = " deg";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/thmax" ){
	  //value = " ";
	  //value += BBthetamax;
	  value.Form(" %5.3g",BBthetamax);
	  unit = " deg";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/phmax" ){
	  //value = " ";
	  //value += BBphimax;
	  value.Form(" %5.3g", BBphimax);
	  unit = " deg";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/hthmin" ){
	  //value = " ";
	  //value += SBSthetamin;
	  value.Form(" %5.3g", SBSthetamin);
	  unit = " deg";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/hthmax" ){
	  // value = " ";
	  // value += SBSthetamax;
	  value.Form(" %5.3g", SBSthetamax);
	  unit = " deg";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/hphmin" ){
	  // value = " ";
	  // value += SBSphimin;
	  value.Form(" %5.3g", SBSphimin);
	  unit = " deg";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/hphmax" ){
	  // value = " ";
	  // value += SBSphimax;
	  value.Form(" %5.3g", SBSphimax);
	  unit = " deg";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/eemin" ){
	  // value = " ";
	  // value += pmin_BB;
	  value.Form(" %5.3g", pmin_BB);
	  unit = " GeV";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/eemax" ){

	  double Eemax = Ebeam/(1.0 + Ebeam/Mp*(1.0-cos(BBthetamin*PI/180.0)));
	  // value = " ";
	  // value += Eemax;
	  value.Form(" %5.3g", Eemax);
	  unit = " GeV";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/ehmin" ){
	  double Ehmin = pmin_SBS;
	  // value = " ";
	  // value += Ehmin;
	  value.Form( " %5.3g", Ehmin);
	  unit = " GeV";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else if( cmd_prefix == "/g4sbs/ehmax" ){
	  double Ehmax = Ebeam - pmin_BB;
	  // value = " ";
	  // value += Ehmax;
	  value.Form(" %5.3g", Ehmax);
	  unit = " GeV";
	  cmd_temp = cmd_prefix + value + unit;
	  outputfile << cmd_temp.Data() << endl;
	} else {
	  outputfile << currentline.Data() << endl;
	}  

      }
      if( ckov_flag == 0 ){
	outputfile << "/process/inactivate Cerenkov" << endl;
      }
      
      TString rootfilename = outfilename;
      rootfilename.ReplaceAll(".mac",".root");

      outputfile << "/g4sbs/filename " << rootfilename.Data() << endl;
      outputfile << "/g4sbs/run " << nevents << endl;
    }
  }

 
  

}
