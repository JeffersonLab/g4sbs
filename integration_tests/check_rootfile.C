#include <vector>

int kNumOfEntries = 100;

std::vector<TString> gCheckHitsOnDetector;
std::vector<int> gCheckGeneratedNucleons;

// The tree
TTree *gTree = 0;

// Checks that these detectors are non-zero
int check_hits_nonzero()
{
  for( size_t i = 0; i < gCheckHitsOnDetector.size(); i++ ) {
    gTree->Draw(TString::Format("%s.hit.nhits>>hCountCheck",
          gCheckHitsOnDetector[i].Data()),
        TString::Format("%s.hit.nhits>0",gCheckHitsOnDetector[i].Data()));
    TH1F *h = (TH1F*)gDirectory->Get("hCountCheck");
    if(!h || h->GetEntries() <= 0) {
      std::cout << "No hits in " << gCheckHitsOnDetector[i].Data() << std::endl;
      return -1;
    } else {
      std::cout << gCheckHitsOnDetector[i].Data() << " has hits." << std::endl;
    }
    delete h;
  }
  return 0;
}


// This is a very simple function that has various sanity tests on the
// produced rootfile.
int standard_tests(const char *fileName)
{
  TFile *file = new TFile(fileName,"READ");
  gTree = (TTree*)file->Get("T");

  // Check if ROOT file and Tree were generated
  if(!gTree) {
    std::cerr << "Tree not found in " << fileName << std::endl;
    return -1;
  }

  // Check the number of entries
  int num_entries = gTree->GetEntries();
  std::cout << "Entries in the tree: " << num_entries
    << " out of " << kNumOfEntries << " expected." << std::endl;
  if(num_entries != kNumOfEntries) {
    std::cerr << "Expected " << kNumOfEntries << " entries in the tree, "
      << "but got only " << num_entries << std::endl;
    return -1;
  }

  // Check that all listed detectors have non-zero counts
  if( check_hits_nonzero() != 0)
    return -1;

  // If we reached this point, then no errors were encountered. Return 0
  return 0;
}

double GetQ2FromFileName(TString fileName)
{
  // Ensure fileName is not null
  if(fileName.IsNull())
    return 0.0;

  // Get rid of the GeV2.root part
  fileName = fileName.ReplaceAll("GeV2.root","");

  // Remove the exp_ part (where exp == gmn, gep, gen, etc...)
  fileName = fileName(4,fileName.Length());
  return fileName.Atof();
}


int perform_Q2Tests(TString fileName)
{
  // Determine the expected Q2
  Double_t Q2 = GetQ2FromFileName(fileName);
  if(Q2==0) {
    std::cout << "Could not determine Q^2 from the file name" << std::endl;
    return -1;
  } else  {
    std::cout << "Checking that resulting Q^2 is " << Q2 << " GeV^2: ";
  }
  TCanvas *c = new TCanvas("c","c",600,600);
  gTree->Draw("ev.Q2>>hQ2Test");
  c->SaveAs("test.png");
  TH1F *hQ2Test = (TH1F*)gDirectory->Get("hQ2Test");
  if(!hQ2Test) {
    std::cout << "Could not determine Q^2 in tree!" << std::endl;
    return -1;
  }
  Double_t tQ2 = hQ2Test->GetMean();
  Double_t tQ2E = hQ2Test->GetMeanError();
  std::cout << " Tree Q^2 is " << tQ2 << " +/- " << tQ2E;
  if(TMath::Abs(tQ2-tQ2) <= 1.0) {
    std::cout << " good match within +/- 1.0 GeV^2" << std::endl;
  } else {
    std::cout << " too different from expected." << std::endl;
    return -1;
  }

  return 0;
}


int standard_ff_tests(TString fileName)
{
  // Perform the standard_tests
  int status = standard_tests(fileName);
  if(status != 0)
    return status;

  status = perform_Q2Tests(fileName);
  if(status != 0)
    return status;

  return 0;
}

int check_rootfile(TString fileName, TString exp, bool opticalPhotons = false)
{
  if(exp.EqualTo("GEn",TString::kIgnoreCase)) {
    // Expect to see both protons and neutrons generated in GEn
    gCheckGeneratedNucleons.push_back(0);
    gCheckGeneratedNucleons.push_back(1);

    // Ensure these detectors have hits
    gCheckHitsOnDetector.push_back("Earm.BBGEM");
    gCheckHitsOnDetector.push_back("Earm.BBPSTF1");
    gCheckHitsOnDetector.push_back("Earm.BBSHTF1");
    gCheckHitsOnDetector.push_back("Harm.HCalScint");

    if(opticalPhotons) {
      gCheckHitsOnDetector.push_back("Earm.BBPS");
      gCheckHitsOnDetector.push_back("Earm.BBSH");
      gCheckHitsOnDetector.push_back("Earm.GRINCH");
      gCheckHitsOnDetector.push_back("Harm.HCal");
    }

    // Run standard tests
    standard_ff_tests(fileName);
  } else if (exp.EqualTo("GEp",TString::kIgnoreCase)) {
    // Expect to see only generated protons in GEp
    gCheckGeneratedNucleons.push_back(1);

    // Ensure these detectors have hits
    gCheckHitsOnDetector.push_back("Harm.FT");
    gCheckHitsOnDetector.push_back("Harm.FPP1");
    gCheckHitsOnDetector.push_back("Harm.FPP2");
    gCheckHitsOnDetector.push_back("Harm.HCalScint");
    gCheckHitsOnDetector.push_back("Earm.CDET_Scint");
    gCheckHitsOnDetector.push_back("Earm.ECalTF1");

    if(opticalPhotons) {
      gCheckHitsOnDetector.push_back("Earm.CDET");
      gCheckHitsOnDetector.push_back("Earm.ECAL");
      gCheckHitsOnDetector.push_back("Harm.HCal");
    }

    // Run standard tests
    standard_ff_tests(fileName);
  } else if (exp.EqualTo("GMn",TString::kIgnoreCase)) {
    // Expect to see both protons and neutrons generated in GMn
    gCheckGeneratedNucleons.push_back(0);
    gCheckGeneratedNucleons.push_back(1);

    // Ensure these detectors have hits
    gCheckHitsOnDetector.push_back("Earm.BBGEM");
    gCheckHitsOnDetector.push_back("Earm.BBPSTF1");
    gCheckHitsOnDetector.push_back("Earm.BBSHTF1");
    gCheckHitsOnDetector.push_back("Harm.HCalScint");

    if(opticalPhotons) {
      gCheckHitsOnDetector.push_back("Earm.BBPS");
      gCheckHitsOnDetector.push_back("Earm.BBSH");
      gCheckHitsOnDetector.push_back("Earm.GRINCH");
      gCheckHitsOnDetector.push_back("Harm.HCal");
    }

    // Run standard tests
    standard_ff_tests(fileName);
  }
}
