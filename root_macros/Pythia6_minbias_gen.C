#include "TPythia6.h"
#include "TFile.h"
#include "TTree.h"
#include "TParticle.h"
#include "TClonesArray.h"
#include <iostream>
#include <fstream>
#include "TString.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;

const double Mp = 0.938272046;
const double me = 0.511e-3;

void Pythia6_minbias_gen( const char *configfilename, const char *outputfilename ){
  ifstream infile(configfilename);

  double Ebeam = 11.0;
  infile >> Ebeam;

  double pbeam = sqrt(pow(Ebeam,2)-pow(me,2));
  
  TVector3 k(0,0,pbeam);
  
  TLorentzVector kbeam( k, Ebeam );
  TVector3 Ptarget(0,0,0);
  TLorentzVector Pproton(Ptarget,Mp);

  TLorentzVector Ptotal = kbeam + Pproton;
  double Ecm = Ptotal.M();

  if( Ebeam < 1.0 ) return;
  
  TPythia6 Generator;
  Generator.SetPARP(2,1.0); //Set minimum CM energy to 1 GeV instead of the default 10 GeV:

  Generator.Initialize( "FIXT", "gamma/e-", "p", 11.0 );

  
}
