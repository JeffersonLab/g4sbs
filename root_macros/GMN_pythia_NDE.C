#include "gmn_pythia_tree.C"
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include <iostream>
#include <fstream>
#include <vector>

//The basic goal of this macro is to identify gamma p --> pi+ n events in BB + HCAL, smear reconstructed quantities by the detector resolutions, and estimate the relative rates of one-pion versus multipion events as a function of various cuts (pion momentum, reconstructed photon energy, etc). 

void GMN_pythia_NDE( const char *infilename, const char *outfilename = "temp.root" ){
  
  TChain *C = new TChain("T");


  TFile *fout = new 
}
