#ifndef G4SBSTDISGen_HH
#define G4SBSTDISGen_HH

// I will try to make it as clean as possible

//#include "TFile.h"
//#include "TTree.h"
//#include "TChain.h"
#include "G4SBSEventGen.hh"
#include "G4SBSIO.hh"

#include "G4SystemOfUnits.hh"
#include "CLHEP/Units/SystemOfUnits.h"
#include "sbstypes.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"

#include "cteq/cteqpdf.h"
#include "DSS2007FF.hh"


using namespace CLHEP;
using namespace G4SBS;


class G4SBSTDISGen //should I use inheritance (CA)
{
public: 
  G4SBSTDISGen();
  ~G4SBSTDISGen();

  G4SBSIO *fIO;
   
  bool Generate(Kine_t, Nucl_t, G4LorentzVector , G4LorentzVector);

  void Initial4Vectors(G4LorentzVector, G4LorentzVector);
  G4ThreeVector GenerateElectronAngles(Kine_t);
  void GenerateScatterElectron(Kine_t, G4LorentzVector, G4LorentzVector, G4ThreeVector);
  void PhotonKinematics(G4LorentzVector);
  void GenerateFinalState(Kine_t, G4LorentzVector);

  void GenerateInitialStateTDIS(Kine_t);
  void GenerateFinalStateTDIS(Kine_t, G4LorentzVector);
  G4double GenerateTDIS(G4LorentzVector, Kine_t);

  void kbinning(G4double , G4double &, G4double &);
  G4double f2_pi(G4double, G4double , G4double );

  G4double F2_N(double x, double Q2, Kine_t Kine);

  void SelectHadron();
  void EventRejected();
  void GenerateSIDIS();
  void GenerateHadronKinematics();
  G4double dissigma(G4double);

  G4ThreeVector DirectionVector(G4double, G4double);
  G4ThreeVector FermiMomentum();
  void FermiDistribution();

  G4double PhotoD_XS(G4double E_photon); //PhotoDisintegration Cross Section
  
  G4double VXPhoton_flux(double photon_energy, double beam_energy);
 
  G4double QuasiElasticXS(G4double beam_energy, G4int z1, G4int n1, G4int partID, G4double momentum, G4double angle);
  G4double FluxCorrection(G4LorentzVector ni, G4LorentzVector ei);

  G4double ElasticXS(G4double beam_energy, G4double scatter_e_energy, G4double theta, Nucl_t nu);

  //do it correctly, full dependences
  G4double SIDISXS(G4LorentzVector ef_NrestL,G4double H1L,G4double H2L, G4double etheta_NrestL);

  G4double MottXS(G4double, G4double);
  G4double DipoleFF();
  G4double GE( Nucl_t);
  G4double GM( Nucl_t);
  G4double tau();
  void qe_epc();
  //  TGraph *gD;

  G4double alpha()
  {return fine_structure_const;} // fine structure, to shorten the definition (useless?)

  void setFinalNucleon(Nucl_t ft) //ft is a temporary variable
  {tFinalNucl = ft;}

  // to be sent to PrimaryGeneratorAction
  Nucl_t tGetFinalNucleon(){ return tFinalNucl; }

  // momentum-energy of the scattered electron to be sent to PrimaryGeneratorAction
  G4ThreeVector tGetElectronP(){ return tElectronP; }
  G4double tGetElectronE(){ return tElectronE; }

  // momentum-energy of the scattered nucleon to be sent to PrimaryGeneratorAction
  G4ThreeVector tGetNucleonP(){ return tNucleonP; }
  G4double tGetNucleonE(){ return tNucleonE; }

  G4double tGetnu(){ return tnu; }
  G4double tGetQ2(){ return tQ2; }
  G4double tGetW2(){ return W2; }
  G4double tGetMX(){ return Mx2; }
  G4double tGetxbj(){ return xbj; }

  G4double tGetKE(){ return KE; }


  G4double tGetQEsigma(){ return QEsigma; }
  G4double tGetELAsigma(){ return ELAsigma; }
  G4double tGetSIDISsigma(){ return SIDISsigma; }

  G4LorentzVector tGetelef_lab(){ return tElectronf_lab; }//to rootfile
  G4LorentzVector tGetnucf_lab(){ return tNucleonf_lab; }//to rootfile


  G4double tGetf2N(){ return f2N; }
  G4double tGetf2pi(){ return f2pi; }
  G4double tGetDISsigma(){ return DISsigma;} //DIS cross section (from TDISkine)
  G4double tGetTDISsigma(){ return TDISsigma;} //TDIS cross section

  G4LorentzVector tGetiProton(){ return iProton; }// TDIS initial proton
  G4LorentzVector tGetiNeutron(){ return iNeutron; }// TDIS initial neutron (Deuterium)(NOT USED)
  G4LorentzVector tGetfProton(){ return fProton; }// TDIS final proton
  G4LorentzVector tGetfPion(){ return fPion; }//to rootfile (NOT USED)

  G4double tGetxpi(){ return xpi; }
  G4double tGettpi(){ return tpi; }
  G4double tGetypi(){ return ypi; }

  G4double tGetpt(){ return pt; } // I should change the name of these variables
  G4double tGetz(){ return z; } // with the prefix 't' to keep the same nomenclature

  G4double tGettxa(){ return txa; }
  G4double tGettya(){ return tya; } // I should drop the 'a' of the name is confussing

  G4LorentzVector tGetHad_lab(){ return mom_had_final; }//stupid getter in order to set the nucleon momentum inside the switch (SHOULD BE A BETTER WAY)


   G4LorentzVector tGetiElectron(){ return ei_lab; }// initial 4-vector electron


  void tSetThMin(double v){tThMin = v;}
  void tSetThMax(double v){tThMax = v;}

  void tSetPhMin(double v){tPhMin = v;}
  void tSetPhMax(double v){tPhMax = v;}
  
  void tSetEeMin(double v){tEeMin = v;}
  void tSetEeMax(double v){tEeMax = v;}


  void tSetThMin_had(double v){tThMin_had = v;}
  void tSetThMax_had(double v){tThMax_had = v;}

  void tSetPhMin_had(double v){tPhMin_had = v;}
  void tSetPhMax_had(double v){tPhMax_had = v;}
  
  void tSetEMin_had(double v){tEMin_had = v;}
  void tSetEMax_had(double v){tEMax_had = v;}

  void tSetHadronType(Hadron_t h ){tHadronType = h; }
  G4int tGetSIDISHadron(){ return SIDISHadron; }

  void tSetThows(double v){NoOfThrows = v;}
  G4double tGetNoOfThrows(){ return NoOfThrows; }

  G4double GetSigma(Kine_t);
  // BAD PROGRAMMING


  void tinitcteqpdf();
  double tF2N(double x, double Q2,  Nucl_t nucl);
  double tdissigma( double ebeam, double th, double eprime, Nucl_t nucl );
  double tdissigma_p(double eb, double th, double ep);
  double tdissigma_n(double eb, double th, double ep);

  double F2N(double x, double Q2,  Nucl_t nucl);

private:

  G4int counter; 

  cteq_pdf_t *__tdis_pdf;

  //maybe we need to deleta most of these definitions

  G4double PI;       // 3.1415...
  G4double e;        // electron charge
  G4double m_deu;    // Deuterom mass (MeV)
  G4double m_pro;    // Deuterom mass (MeV)
  G4double deu_bind; // deuterium binding energy (MeV)
  G4double h_bar;    // reduced Planck Constant (MeV s)
  G4double c;        // speed of light
  G4double m_e;      // mass of the electron
  //maybe we need to deleta most of these definitions

  G4double Mp, Mn,  Mt; // proton, neutron (target p/n)

  G4double Eeprime_lab, Peprime_lab; 

  G4double KE;

  G4double NoOfThrows;

  G4ThreeVector boost_Nrest; //boost from LAB to Nucleon Rest frame
  G4LorentzVector ei_Nrest, ni_Nrest;
  G4LorentzVector Pisum_lab;
  G4double Ebeam_Nrest, Ebeam_lab;

  G4double th, ph; //the random generated electron angles

  G4LorentzVector ei_lab;

  G4ThreeVector null3D;//= G4ThreeVector(0,0,0);//is there a null vector by definition?

  G4ThreeVector kfhat_lab;
  G4LorentzVector ef_lab, q_lab;
  G4LorentzVector ef_Nrest;
  G4LorentzVector nf_Nrest;
  G4double Eprime_Nrest;

  G4double eth_Nrest;// scatter electron angle in Nucleon rest frame
  G4double thN_Nrest;//knock-out nucleon theta angle in Nucleon rest frame (QE case)

  G4LorentzVector q_Nrest;
  G4LorentzVector nf_lab;

  G4double hTheta, hPhi, hE, hP; //tSIDIS variables
  G4LorentzVector Phad_lab, Phad_Nrest;
  G4LorentzVector Pfsum_lab; // final state 4 momentum

  G4int icharge, ihadron;
  G4LorentzVector  mom_had_final;
  G4ThreeVector Phad_Nrest_vect, q_Nrest_vect, Phad_perp;
  G4double Ph_perp, tb;

  G4double pt, z; // TAKE CARE, z is defined diferently in SIDIS and TDIS

  G4double txa; //in g4sbs is xbj BUT a bit different and it is stored separetly

  G4LorentzVector iProton, iNeutron; //initial proton and neutron for the TDIS case
                                     // only proton for tTDISKinH and both for tTDISKinD
  G4LorentzVector fProton, fPion; // the final proton and pion for the TDIS case

  G4double tpi, ypi, fpi, xpi; 

  G4double f2pi; //pion structure function (parametrization, not sure!)
  G4double f2N; //nucleon structure function
  G4double DISsigma, TDISsigma;


  G4double H1, H2; //SIDIS structure functions
  DSS2007FF tFragFunc;

  G4double tThMin, tThMax, tPhMin, tPhMax; //Angular generation limits for electron arm
  G4double tEeMin, tEeMax; //Electron energy generation limits

  G4double tThMin_had, tThMax_had, tPhMin_had, tPhMax_had; //Angular generation limits for hadron
  G4double tEMin_had, tEMax_had; //Hadron energy generation limits
  G4double Mh;
  G4double nu_Nrest;
  G4double etheta_Nrest,theta_pq_Nrest;

  G4double iGE; // internal calculus G
  G4double iGM; // internal calculus GM

  G4double tQ2; // momentum transfer

  G4double tya; //

  // Nucl_t tNuclType, tFinalNucl;
 
  G4double FluxC; //flux correction from Nucleon Rest to Lab frames
  
  Nucl_t tNuclType, tFinalNucl;
 
  Hadron_t tHadronType; 
  G4int SIDISHadron;//I don't know what I can't access the hadron

  G4double W2, xbj, Mx2;
  G4ThreeVector tElectronP, tNucleonP;
  G4double tElectronE, tNucleonE;
  G4double  QEsigma, ELAsigma, SIDISsigma;
  G4double tnu;
  G4LorentzVector tElectronf_lab, tNucleonf_lab;










  // TWO FIRST COLUMNS FOR FERMI DISTRIBUTION MOMENTA 
  // FROM FILE moment_ld2b.dat
  // THE FILE CONTAINS TWO MORE COLUMS BUT THERE IS NO DOCUMENTATION ABOUT IT

  G4double f_mom[201]={0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, //0-9
		     0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, //10-19
		       1, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, //20-29
		     1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 
		       2, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 
		     2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 
		       3, 3.05, 3.1, 3.15, 3.2, 3.25, 3.3, 3.35, 3.4, 3.45, 
		     3.5, 3.55, 3.6, 3.65, 3.7, 3.75, 3.8, 3.85, 3.9, 3.95, 
		       4, 4.05, 4.1, 4.15, 4.2, 4.25, 4.3, 4.35, 4.4, 4.45, 
		     4.5, 4.55, 4.6, 4.65, 4.7, 4.75, 4.8, 4.85, 4.9, 4.95, 
		       5, 5.05, 5.1, 5.15, 5.2, 5.25, 5.3, 5.35, 5.4, 5.45, 
		     5.5, 5.55, 5.6, 5.65, 5.7, 5.75, 5.8, 5.85, 5.9, 5.95, 
		       6, 6.05, 6.1, 6.15, 6.2, 6.25, 6.3, 6.35, 6.4, 6.45, 
		     6.5, 6.55, 6.6, 6.65, 6.7, 6.75, 6.8, 6.85, 6.9, 6.95, 
		       7, 7.05, 7.1, 7.15, 7.2, 7.25, 7.3, 7.35, 7.4, 7.45, 
		     7.5, 7.55, 7.6, 7.65, 7.7, 7.75, 7.8, 7.85, 7.9, 7.95, 
		       8, 8.05, 8.1, 8.15, 8.2, 8.25, 8.3, 8.35, 8.4, 8.45,
		     8.5, 8.55, 8.6, 8.65, 8.7, 8.75, 8.8, 8.85, 8.9, 8.95, 
		       9, 9.05, 9.1, 9.15, 9.2, 9.25, 9.3, 9.35, 9.4, 9.45, 
		     9.5, 9.55, 9.6, 9.65, 9.7, 9.75, 9.8, 9.85, 9.9, 9.95, 10};


 G4double f_distr[201]={3214,  2925.7,       2257.2,     1554.5,     1007.6,      638.71,     405.07,     260.04,     169.84,     113.03,        //0-9
			76.63,     52.87,      37.071,     26.381,     19.029,     13.896,     10.263,      7.66,       5.7737,     4.3925,     //10-19
			3.3716,     2.6103,     2.038,      1.6045,     1.2738,     1.0198,     0.82343,    0.67066,    0.55103,    0.45678,    //20-29
			0.38205,    0.3224,     0.27448,    0.23571,    0.20412,    0.17819,    0.15673,    0.13883,    0.12379,    0.11103,    //30-39
			0.10012,    0.090719,   0.082546,   0.075386,   0.069066,   0.063448,   0.058422,   0.053898,   0.049803,   0.046079,   //40-49
			0.042678,   0.03956,    0.036691,   0.034046,   0.031601,   0.029336,   0.027234,   0.025283,   0.023468,   0.021779,   //50-59
			0.020206,   0.018742,   0.017377,   0.016105,   0.014921,   0.013817,   0.012788,   0.011831,   0.010939,   0.010109,   //60-69
			0.009337,   0.008619,   0.0079516,  0.0073315,  0.0067557,  0.0062213,  0.0057255,  0.005266,   0.0048401,  0.0044458,  //70-79
			0.0040809,  0.0037433,  0.0034313,  0.0031431,  0.0028771,  0.0026316,  0.0024053,  0.0021968,  0.0020048,  0.0018282,  //80-89
			0.0016658,  0.0015165,  0.0013795,  0.0012538,  0.0011385,  0.001033,   0.00093631, 0.00084789, 0.00076709, 0.00069329, //90-99
			0.00062596, 0.00056457, 0.00050865, 0.00045776, 0.00041148, 0.00036944, 0.0003313,  0.00029671, 0.00026538, 0.00023704, //100-109
			0.00021143, 0.0001883,  0.00016745, 0.00014867, 0.00013177, 0.00011659, 0.00010297, 9.0774E-05, 7.9857E-05, 7.0106E-05, 
			6.1408E-05, 5.3663E-05, 4.6778E-05, 4.0669E-05, 3.5259E-05, 3.0477E-05, 2.6259E-05, 2.2548E-05, 1.929E-05,  1.6438E-05, 
			1.3947E-05, 1.1778E-05, 9.8963E-06, 8.2688E-06, 6.8666E-06, 5.6636E-06, 4.6361E-06, 3.763E-06,  3.0254E-06, 2.4062E-06, 
			1.8902E-06, 1.4639E-06, 1.1152E-06, 8.3336E-07, 6.0894E-07, 4.335E-07,  2.9961E-07, 2.0075E-07, 1.3119E-07, 8.5916E-08, 
			6.0567E-08, 5.1347E-08, 5.4978E-08, 6.8637E-08, 8.9908E-08, 1.1673E-07, 1.4737E-07, 1.8036E-07, 2.1449E-07, 2.4877E-07,
			2.8237E-07, 3.1466E-07, 3.4513E-07, 3.7339E-07, 3.9918E-07, 4.2229E-07, 4.4261E-07, 4.601E-07,  4.7475E-07, 4.8661E-07, 
			4.9576E-07, 5.023E-07,  5.0637E-07, 5.0809E-07, 5.0765E-07, 5.0519E-07, 5.0087E-07, 4.9488E-07, 4.8738E-07, 4.7852E-07, 
			4.6848E-07, 4.574E-07,  4.4543E-07, 4.3272E-07, 4.1938E-07, 4.0555E-07, 3.9135E-07, 3.7686E-07, 3.6221E-07, 3.4746E-07, 
			3.327E-07,  3.1801E-07, 3.0345E-07, 2.8908E-07, 2.7495E-07, 2.611E-07,  2.4757E-07, 2.3439E-07, 2.216E-07,  2.092E-07, 1.9722E-07};




};

extern G4SBSTDISGen *tdishandler;

#endif//G4SBSTDISGen_HH
