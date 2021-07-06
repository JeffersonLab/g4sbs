// TDIS Generators class 
// Carlos Ayerbe (gayoso@jlab.org)


#include "TMath.h"
#include <vector>
#include "TF1.h" //test

#include "CLHEP/Random/RandGauss.h"
//#include "cteq/cteqpdf.h"
#include "globals.hh"
#include "CLHEP/Random/RandFlat.h"
#include "G4SBSTDISGen.hh"
#include "G4SBSIO.hh"
#include "G4SystemOfUnits.hh"
#include "G4SBSEventGen.hh"
#include "G4Proton.hh" 
#include "G4PionPlus.hh"

#define TH TRUE // use Tim Hobbs code
#define DEBUG FALSE


// Should be in the header? (CA)
extern"C"{
  double epc_func_(double *, int*, int*, int*, double *, double *);
}


// this one was coded as subroutine
extern"C"{
  void f2pi_sub_(double *, double*, double*, double*, double*, double * );
  // parameteters to send:
  //  energy beam, xbj, electron scattering angle, min k, max k
  // parameters to receive:
  // f2_pi (the numerator to calculate the sigma TDIS)
}


// WHY SHOULD BE DECLARED GLOBAL?

using namespace CLHEP;


G4SBSTDISGen::G4SBSTDISGen()
{
  //TAKE CARE!! The handler is initialized practically BEFORE all the classes.
  // any definition from messenger will appear wrong


  if(DEBUG)
    G4cout<<"<G4SBSTDISGen::G4SBSTDISGen()>: initializing"<<G4endl;
  // A bunch of constants (some of them should be replaced by CLHEP)




  PI = TMath::Pi();  //3.14159

  e = 1.602e-19;     // electron charge
  m_deu = 1875.6;    // Deuterom mass (MeV)
  m_pro = 938.3;     // Deuterom mass (MeV)
  deu_bind = 2.224;  // deuterium binding energy (MeV)
  h_bar = 6.582e-22; // reduced Planck Constant (MeV s)
  c = 299792458;     // speed of light
  m_e = 0.511;       // mass of the electron (MeV)
 

  // these are useless, we can get the mass of the particle
  // with ---> ni_Nrest.m()
  Mp = proton_mass_c2; // I can't use PDG definition HERE because it was not defined when this class is initiazed
  Mn = neutron_mass_c2; 
  //  iQ2 = 0;

  tQ2 = 0;

  counter = 0;

  if(DEBUG)
    G4cout<<"Initializing"<<G4endl;
  //  init DIS cteq pdf

  tinitcteqpdf();
  if(DEBUG)
    G4cout<<"crashing"<<G4endl;
}

G4SBSTDISGen::~G4SBSTDISGen()
{
  if(DEBUG)
    G4cout<<"G4SBSTDISGen::~G4SBSTDISGen()"<<G4endl;
}

G4SBSTDISGen *tdishandler=NULL;


bool G4SBSTDISGen::Generate(Kine_t tKineType, Nucl_t nucl, G4LorentzVector ei, G4LorentzVector ni)
{

  if(DEBUG)
    G4cout<<"Welcome to the generators"<<G4endl;
  if(DEBUG)
    G4cout << "GENERATETDIS" << G4endl;

  if(DEBUG)
    G4cout<<"Entering  <G4SBSTDISGen::Generate>"<<G4endl;
  if(DEBUG)
    G4cout<<"kine: "<< tKineType <<G4endl; 
  if(DEBUG)
    G4cout<<"nucl: "<< nucl <<G4endl; 
  if(DEBUG)
    G4cout<<"ei: "<< ei <<G4endl; 
  if(DEBUG)
    G4cout<<"ni: "<< ni <<G4endl; 
  //note that the units are in MeV, which I understand are the natural units in Geant4

  ei_lab = ei;

  // to keep coherence with the nomenclature used, I should change n=0, p=1
  if(nucl == 0)
    {
      if(DEBUG)
	G4cout << "nucl: "<<  nucl << " A NEUTRON" << G4endl;
      Mt = Mn;
    }
  else
    {
      if(DEBUG)
	G4cout << "nucl: "<<  nucl << " A PROTON" << G4endl;
      Mt = Mp;
    }
   

    if (counter == 0)
    {
      if(DEBUG)
	G4cout <<" tThMin = "<<tThMin/deg<<G4endl;
      if(DEBUG)
	G4cout <<" tThMax = "<<tThMax/deg <<G4endl;
      
      if(DEBUG)
	G4cout <<" tPhMin = "<<tPhMin/deg <<G4endl;
      if(DEBUG)
	G4cout <<" tPhMax = " <<tPhMax/deg <<G4endl;
      
      if(DEBUG)
	G4cout <<" tEeMin = "<<tEeMin <<G4endl;
      if(DEBUG)
	G4cout <<" tEeMax = "<<tEeMax <<G4endl;
      
      
      if (tKineType == tSIDIS)
	{
	  if(DEBUG)
	    G4cout << " tThMin_had  = "<< tThMin_had/deg << G4endl;
	  if(DEBUG)
	    G4cout << " tThMax_had  = "<< tThMax_had/deg << G4endl; 
	  
	  if(DEBUG)
	    G4cout << " tPhMin_had  = "<< tPhMin_had/deg << G4endl;
	  if(DEBUG)
	    G4cout << " tPhMax_had  = "<< tPhMax_had/deg << G4endl; 
	  
	  if(DEBUG)
	    G4cout << " tEMin_had  = "<< tEMin_had << G4endl;
	  if(DEBUG)
	    G4cout << " tEMax_had  = "<< tEMax_had << G4endl; 
	}
      
      counter ++;
    }

  // I put the initial 4-vectors treatment in a separate method
  // since it is common for all reactions

  // this method receive the 4-vector beam and 4-vector nucleon
  // WHICH could be at rest (H likewise) or Fermi smeared (CHECK!)

  // So, we have two frames, Lab system and Nucleon rest
  // both frames coincides IF the nucleon is NOT Fermi smeared

  Initial4Vectors(ei, ni); 
 
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  //GENERATE KINEMATICS (copy from G4SBSEventGen)
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  // I put it apart. Also common for many reactions
  // It produces the randomized electron angles and
  // an unitary vector in such direction

  // perhaps it is not very smart...
  kfhat_lab = GenerateElectronAngles(tKineType);

  if(DEBUG)
    G4cout<<"kfhat_lab: "<< kfhat_lab <<G4endl; 

  // ni here is just for the elastic and QE cases
  GenerateScatterElectron(tKineType, ei, ni, kfhat_lab);

  if(DEBUG)
    G4cout<<" ef_lab_momentum: "<< ef_lab.vect().unit()<<G4endl;

  // there is a not very good programming coding
  // I am using the variables global-ish.
  // ei and ni are not global.
  // Not sure if it is ok (CA)
  // The name Photon Kinematics is used in a very wide way
  PhotonKinematics(ei);

  if(tKineType != tTDISKinD &&tKineType != tTDISKinH )
    {
      GenerateFinalState(tKineType, ni);
    }
  else
    {
      GenerateInitialStateTDIS(tKineType);
      GenerateFinalStateTDIS(tKineType, ei);
    }

  // Note: we need to check the units of this correction
  // and take care when is applicable. Like this, it is only useful in the elastic case.
  // 26.05.2020 aparently it will also works for SIDIS
  FluxC = FluxCorrection(ni, ei); // calculate the flux correction for this kinematics

  // In principle, each kinematic case just need to carry the 
  // values calculated until here. Then each case will be 
  // treated differently, maybe just one line, maybe a whole 
  // method (as the pion structure generator will need)

   if(DEBUG)
    G4cout<<"kine: "<<tKineType<<G4endl; 

  // REPAIR THIS!!!
  //These values are set fixed now for Deuterium target with the epc code
  // later on, the target used should be carried here
  // NOTE: I don't know why but these line cannot be inside the switch    
  G4int z1 = 1; //atomic number (number of protons)
  G4int n1 = 1; //THIS IS NUMBER OF NEUTRONS!!! 
  G4int partID = (Nucl_t) nucl; 

  // PART OF SIDIS FOR A PROTON OR PI+

  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  //Cross sections 
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

  switch(tKineType)
    {
    case tElastic:
      xbj = 1.0; //force the x_Bj value for elastic
      if(DEBUG)
	G4cout<<Ebeam_Nrest<<" "<<Eprime_Nrest<<" "<<eth_Nrest<<" "<<nucl<<G4endl;
      
      // EPAF: I think we want epc function for both elastic and quasielastic
      ELAsigma = ElasticXS(Ebeam_Nrest, Eprime_Nrest, eth_Nrest, nucl);
      //QuasiElasticXS(Ebeam_Nrest, z1, n1, partID, Phad_Nrest .vect().mag(), thN_Nrest);
      
      if(DEBUG)
	G4cout<<"enter in TDIS elastic: "<<ELAsigma/barn<<G4endl;
      break;
      
      //------------------
      //QUASIELASTIC CASE
      //------------------
    case tQuasiElastic:
   
      xbj = 1.0; //force the x_Bj value for elastic (DD suggestion, makes sense)

      //always in Nucleon Rest Frame
      thN_Nrest = Phad_Nrest.theta() ; //nucleon theta angle in N rest frame
      
      QEsigma = QuasiElasticXS(Ebeam_Nrest, z1, n1, partID, Phad_Nrest .vect().mag(), thN_Nrest);
      
      if(DEBUG)
	G4cout<<"[TDIS] quasielastic sigma (uB/(MeV sr): "<<QEsigma<<G4endl;
      // This is from definition BUT it is NOT clear for me. Nevertheless, it coincides with the epc code
      // value 'tp' which I believe is the kinetic energy  (CA)
      if(DEBUG)
	G4cout<<" KE from definition:  "<< sqrt( pow(Phad_Nrest.vect().mag(),2) + ni_Nrest.m2()) - ni_Nrest.m()<<G4endl;   
      KE =  sqrt( pow(Phad_Nrest.vect().mag(),2) + ni_Nrest.m2()) - ni_Nrest.m();
      
      break;
    

      //-----------------------------------------
      // SIDIS CASE (protons or pion final state)
      //-----------------------------------------
    case tSIDIS:
      
      if(DEBUG)
	G4cout<<"enter in TDIS SIDIS: "<<G4endl;
      
      SelectHadron();
      if(DEBUG)
	G4cout << "(CASE) Hadron: " << tHadronType << " ihadron: " << ihadron << G4endl;
      
      // 1. the sum of outgoing electron and hadron energies cannot exceed the 
      //   incoming electron and nucleon energies (assuming the collision takes place on a single nucleon):
      // 2. the momentum fraction should be 0<=x<=1
      // 3. SIDIS variables where kinematically forbidden (z>1)
      
      if(DEBUG)
	G4cout<<"SIDIS: Pfsum_lab.m2() "<<  Pfsum_lab.m2()<< " xbj: "<<xbj<< " z: "<< z <<" Q2: "<< tQ2/(GeV*GeV) << G4endl;
      
      //the cases Kinematically forbidden --> abort:
      if( Pfsum_lab.m2() >  Pisum_lab.m2() || 
	  Pfsum_lab.e()  >  Pisum_lab.e()  ||
	  xbj >= 1.0                       || 
	  xbj <= 0.0                       ||
	  z > 1.0                          ||
	  tQ2/(GeV*GeV) < 1.0 ) 
	//extra condition. Values,  Q^2<0.25 or x_Bj<~0.04 give some problems calculating the PDFs
	//Check Andrew and Eric answer about it (Carlos logbook)
	{
	  if(DEBUG)
	    G4cout<<"SIDIS: EVENT REJECTED!! "<<G4endl;
	  EventRejected();
	  return false;
	}	
      GenerateSIDIS();
      
      SIDISsigma = SIDISXS(ef_Nrest, H1, H2, etheta_Nrest );//I prefer to put the cross section formula separated (CA)
      if(DEBUG)
	G4cout<<"[TDIS] SIDIS sigma (units?): "<<SIDISsigma/cm2*pow(GeV,2)<<G4endl;
      
      break;


      //-------------------------------------------
      // TDIS CASE (DK parametrization or TH  code)
      //-------------------------------------------

    case tTDISKinH:
    case tTDISKinD:
      if(DEBUG)
	G4cout<<"enter in TDIS case: "<<G4endl;
      
      f2pi = GenerateTDIS(ei, tKineType); // I send the kine variable, which could be H or D
      if(DEBUG)
	G4cout  <<" f2pi: "<<f2pi<<G4endl;
      if (f2pi<=0) //not sure if the calculations give negative numbers
	{
	  if(DEBUG)
	    G4cout<<"TDIS: EVENT REJECTED!! "<<G4endl;
	  return false;
	}

      if(DEBUG)
	if( Pfsum_lab.m2() >  Pisum_lab.m2() || 
	    Pfsum_lab.e()  >  Pisum_lab.e()  )
	  {
	    G4cout<<"Pfsum_lab.m2() : "<<Pfsum_lab.m2() << " Pfsum_lab.e() : "<<Pfsum_lab.e() <<G4endl;
	    G4cout<<"Pisum_lab.m2() : "<<Pisum_lab.m2() << " Pisum_lab.e() : "<<Pisum_lab.e() <<G4endl;
	  }
      
      //sometimes f2pi = 0, possible because it is out of the xbj range where it is defined
      // I think that similarly to the SIDIS case, where sometimes the case was rejected
      // here I should include a return false

      DISsigma = 0;
      TDISsigma = 0;
      
      if (xbj > 0 && xbj < 1.0) // this condition should be checked. In principle xbj should be always <0.3
	{
	  // for the neutron/proton target
	  f2N   = F2_N(xbj, tQ2, tKineType);// proton or neutron determined inside the function   

	  DISsigma = dissigma(f2N)*nanobarn; // the number calculated in the function is in nb, but it doesn't have an unit

	  if(DEBUG)
	    G4cout<<"dissigma: "<<  DISsigma <<" dissigma(cm^2): "<<  DISsigma/cm2<<G4endl;

	  // DISsigma = dissigma(f2N) *nanobarn; // in g4sbs the function needs Energy Beam, Theta_e', Energy e' 
	  // and it was sent in GeV. But also f2N was recalculated inside. Since I am using
	  // global variables, and to save code, I just send f2N.
	  // I repeated the same proceadure has G4SBSDIS.hh. Making use of internal units was creating a mess
	  // at least like this, we know that everything is consistent/
	 
	  TDISsigma = DISsigma * (f2pi/f2N);
	}
      
      if(DEBUG)
	G4cout<<"DISsigma(cm2): "<<DISsigma/cm2<<" DISsigma(nb): "<<DISsigma/nanobarn<<" f2pi: "<<f2pi<<" f2N: "<<f2N<<" TDISsigma: "<<TDISsigma/cm2<<G4endl; 

      break;
      
    case tInelastic:
      
      break;
      
    default:
      break;
    }






  // all this should be moved to a method WHERE we put the variables to be stored in the rootfile



  G4double x_Nrest = -q_Nrest.m2()/(2.0*ni_Nrest.dot( q_Nrest ) );
  

  
  //  tNucleonf_lab = nf_lab; // to rootfile


  //this is the 4-vector knock-out nucleon (QE) or the one produced in SIDIS
  // BUT it should be differentiated, since in QE could be p or n, but only protons in SIDIS
  tNucleonf_lab = Phad_lab ;// to rootfile 


  if(DEBUG)
    G4cout<<"Phad : "<< Phad_Nrest << G4endl;
  if(DEBUG)
    G4cout<<"Phad 3-mon: "<< Phad_Nrest.vect().mag() <<" Phad ene: "<< Phad_Nrest.e()<<G4endl;
  //****** more kinematic values*****

  // 

  // double xa =  Q2/(2*Mt*nu);
  // double ya =  nu/EBeam;

  //********************************


  // momentum-energy of the scattered electron to be sent to PrimaryGeneratorAction
  tElectronP = ef_lab.vect();
  tElectronE = ef_lab.e();

  // momentum-energy of the scattered nucleon to be sent to PrimaryGeneratorAction
  tNucleonP = Phad_lab.vect();
  tNucleonE = Phad_lab.e();

  tFinalNucl = nucl;
  if(DEBUG)
    G4cout<<"nucleon sent: "<<tFinalNucl<<G4endl;


  return true;
}


void G4SBSTDISGen::Initial4Vectors(G4LorentzVector ei, G4LorentzVector ni)
{
  //NOT SURE ABOUT KEEPING THE SAME NAME FOR THE VARIABLES

  //This proceadure is common in EventGen for the
  // Elastic, Inelastic, SIDIS, TDIS
  
  if(DEBUG)
    G4cout<<"ei: "<< ei <<G4endl; 
  if(DEBUG)
    G4cout<<"ni: "<< ni <<G4endl; 

  // This creates a boost to a rest Nucleon in case it is Fermi smeared
  boost_Nrest = ni.boostVector();
  // boostVector() return the spatial coordinates divided by the time component
  //if the nucleon is on rest, it will return 0

  // Total initial 4-momentum in lab
  Pisum_lab = ei + ni;

  // 4-momentum vectors in Nucleon rest frame
  // defined FIRST from the lab frame, them boosted to Nucleon rest
  ei_Nrest = ei;
  ni_Nrest = ni;

  ei_Nrest.boost( -boost_Nrest );
  ni_Nrest.boost( -boost_Nrest );

  if(DEBUG)
    G4cout<<"ei_Nrest: "<< ei_Nrest <<G4endl; 
  if(DEBUG)
    G4cout<<"ni_Nrest: "<< ni_Nrest <<G4endl; 

  // Note that, if the Nucleon is not smeared, 
  // ei_Nrest = ei; ni_Nrest = ni 

  //electron beam energy in Nucleon rest frame
  Ebeam_Nrest = ei_Nrest.e();

  //electron beam energy in Lab frame
  Ebeam_lab = ei.e();

  if(DEBUG)
    G4cout<<"ei_Nrest.e(): "<< ei_Nrest.e() <<G4endl; 
  if(DEBUG)
    G4cout<<"ei.m2(): "<< ei.m2() <<G4endl; 


  return;

}

G4ThreeVector G4SBSTDISGen::GenerateElectronAngles(Kine_t Kine)
{
  G4double ThetaMax, ThetaMin, PhiMax, PhiMin;

  if(Kine != tTDISKinD && Kine != tTDISKinH )
    {
      ThetaMax = tThMax;
      ThetaMin = tThMin;
      PhiMax   = tPhMax; //Note: phi is generated WITHOUT considering being detected by sbs (electron arm)
      PhiMin   = tPhMin; //if that is the case (detect the electron), phi should be rotated 180 deg (see TDIS case) 
   if(DEBUG)
      G4cout<<"PhiMax: "<<PhiMax<<" PhiMin: "<<PhiMin<< G4endl;
   if(DEBUG)
      G4cout<<"tPhMax: "<<tPhMax/deg<<" tPhMin: "<<tPhMin/deg<< G4endl;
      th = acos( CLHEP::RandFlat::shoot(cos(ThetaMax), cos(ThetaMin)) );
    }
  else
    {
      //SBS acceptance (hard coded) for TDIS case
      ThetaMax =  17*deg;
      ThetaMin =  12*deg;
      PhiMax   = (180-12)*deg;
      PhiMin   = (180+12)*deg;
     // redefinition for FIXED angle
      //ThetaMax = ThetaMin = 35*deg; 
      th = acos( CLHEP::RandFlat::shoot(cos(ThetaMax), cos(ThetaMin)));
      //   th = 35*deg;
    }

  //Generate both theta and phi angles in the LAB frame:

  ph = CLHEP::RandFlat::shoot(PhiMin, PhiMax); 

  if(DEBUG)
    G4cout<<"Electron th: "<<th/deg<<" Electron ph: "<<ph/deg<< G4endl;

  //  G4ThreeVector direction_vec ( sin(th)*cos(ph), sin(th)*sin(ph),  cos(th) );
  //unit vector in the direction of the scattered electron in the LAB frame:
 
  return DirectionVector(th, ph);
}

void G4SBSTDISGen::GenerateScatterElectron(Kine_t Kine,  G4LorentzVector ei, G4LorentzVector ni, G4ThreeVector direction_vec)
{
  // Outgoing energy of the scattered electron in the LAB frame accounting for 
  // the initial nucleon motion (no off-shell or binding energy corrections, just Fermi momentum)
    
    
  // I choose a switch since I am not sure if there will be different cases
  // to randomize the scattered electron energy


  //Generate electron and hadron angles and energies in the LAB frame:
  //These will then be boosted to the nucleon rest frame to compute the differential cross section.


  //NOTE: the use of ni (initial nucleon 4-vector) here is just to calculate the electron energy
  // in the elastic or quasielastic case. 


  switch (Kine)
    {
    case tElastic:
    case tQuasiElastic:
      if(DEBUG)
	G4cout<< " Ela or QE electron energy"<< G4endl;
      Eeprime_lab = (ei.e()*(ni.e()-ni.pz()))/(ei.e()*(1.-cos(th))+ni.e()-ni.vect().dot( direction_vec ));
      break;

    case tSIDIS:
    case tTDISKinH:
    case tTDISKinD:
      Eeprime_lab = CLHEP::RandFlat::shoot(tEeMin, tEeMax );
      if(DEBUG)
	G4cout<< " SIDIS/TDIS electron energy (LAB): "<< Eeprime_lab << G4endl;
      break;
      
    default:    
      if(DEBUG)
	G4cout<< " This case should never be activated"<< G4endl;
      break;
    }

  //scattered electrom momentum module (Lab System) (is that correct?)
  Peprime_lab = sqrt(pow(Eeprime_lab, 2) - ei.m2());
    
  if(DEBUG)
    G4cout<<"Eeprime_lab: "<< Eeprime_lab <<G4endl; 
  if(DEBUG)
    G4cout<<"Peprime_lab: "<< Peprime_lab <<G4endl; 
  
    
  // 3-momentum scatter electron in lab frame:
  G4ThreeVector kf_lab = Peprime_lab* direction_vec;

  if(DEBUG)
    G4cout<<"kf_lab: "<< kf_lab <<G4endl; 
  G4LorentzVector electron_final_lab( kf_lab, Eeprime_lab );

  //Four-momentum of scattered electron in the LAB frame:
  ef_lab = electron_final_lab; // there should be an easier way (CA)


  // Calculate four-momentum of scattered electron boosted to the 
  // nucleon REST frame for cross section calculation:
  ef_Nrest = ef_lab; // first define the 4-vector...
  ef_Nrest.boost( -boost_Nrest ); // then boost to the Nucleon rest frame

  //scattered electron energy in Nucleon Rest frame
  Eprime_Nrest = ef_Nrest.e();  

  tElectronf_lab = ef_lab; //to rootfile


  // this is the angle used to calculate Cross-sections 
  // at least Elastic and Mott (for the moment)

  //Compute the e- scattering angle in the nucleon rest frame:
  // it is necesary in TDIS-TH
  eth_Nrest = acos( ei_Nrest.vect().unit().dot( ef_Nrest.vect().unit()) );

    
  return; // I would love to return the two vectors created here (Lab and Nrest)
}



// I know!! terrible programming using global variables
void G4SBSTDISGen::PhotonKinematics(G4LorentzVector ei)
{
  //q vector in the LAB frame (virtual photon momentum):
  q_lab = ei - ef_lab;
  if(DEBUG)
    G4cout<<"ei: "<< ei<<" ef_lab: "<<ef_lab <<G4endl; 
  if(DEBUG)
    G4cout<<"q_lab: "<< q_lab <<G4endl; 
  // Q2 definition (Lorentz Invariant)
  tQ2 = -q_lab.m2(); // set the variable to be used along the whole class
 
  tnu =  q_lab.e(); // return Energy component
 
  tya =  tnu/ ei.e(); //the fractional energy loss of the electron

  if(DEBUG)
    G4cout<<"Photon Kinematics. Q2: "<< tQ2/(GeV*GeV)<< " (GeV^2) "<<" nu: "<< tnu<<G4endl;

  return;
}


void G4SBSTDISGen::GenerateFinalState(Kine_t Kine, G4LorentzVector ni)
{
  // x-Bjorken uses the initial nucleon 4-vector.
  // this vector is calculated differently for TDIS, thus xbj is different.
  // not much, but enough. This is the reason I put it here.

  //x-Bjorken (Lorentz Invariant) (I cross checked using both frames)
  xbj = tQ2/(2.0 * ni.dot( q_lab ) );

  if(DEBUG)
    G4cout<<"xbj(Photon): "<<xbj<<G4endl;


  //Mh is now defined by the hadron case

  // Mh = Mp; //this is a very restrictive SIDIS case to just protons
  //defining Mh as the mass of the nucleon is more general... BUT how affects the 
  // SIDIS case   


  
  // technically this is the virtual photon in the nucleon rest frame
  // and I defined all the photon properties in other place (CA)
  q_Nrest = ei_Nrest - ef_Nrest;
  
  switch(Kine)
    {
    case tElastic:
    case tQuasiElastic:
      Phad_Nrest = ni_Nrest + q_Nrest; // only true in Elastic and QE
      
      //boost back to LAB frame
      Phad_lab = Phad_Nrest; //copy the 4-vector
      Phad_lab.boost( boost_Nrest ); //then boost to LAB frame
      
      // other cases
      break;
      
    case tSIDIS:

      if(DEBUG)
	G4cout<<"GenerateFinalState: SIDIS"<<G4endl;
      GenerateHadronKinematics();

      Phad_lab = tGetHad_lab();//stupid getter in order to set the nucleon momentum inside the switch (SHOULD BE A BETTER WAY)

      Phad_Nrest = Phad_lab;  //copy the 4-vector
      Phad_Nrest.boost( -boost_Nrest ); //then boost to LAB frame
      break;
      
    default:
      if(DEBUG)
	G4cout<< " This case should never be activated"<< G4endl;
      break;
    }

  //I can't use the CLHEP vectors inside the switch case, so I need to make the declaration out of it 
  
  //Some of these variables are only useful for SIDIS, but it won't hurt putting out of the class
  // perhaps, slower the run?

  Phad_Nrest_vect = Phad_Nrest.vect();

  nu_Nrest = q_Nrest.e();
  q_Nrest_vect    = q_Nrest.vect();

  Phad_perp = Phad_Nrest_vect - ( Phad_Nrest_vect.dot(q_Nrest_vect)/q_Nrest_vect.mag2() ) * q_Nrest_vect ;
  
  Ph_perp = sqrt( Phad_perp.mag2() );

  Pfsum_lab = Phad_lab + ef_lab;

  if(DEBUG)
    G4cout<<"Phad_lab "<<Phad_lab<< " ef_lab: "<<ef_lab<<G4endl;

  if(DEBUG)
    G4cout<<"Phad_Nrest "<<Phad_Nrest<< " ni_Nrest: "<<ni_Nrest<<G4endl;

  //Compute SIDIS kinematic quantities
  // z = P dot Ph / P dot q:
  z = ni_Nrest.dot( Phad_Nrest ) / ni_Nrest.dot( q_Nrest ); //This quantity is also Lorentz-invariant
  // does it makes sense or could provides problems if the case is not SIDIS? (check with cout)


  W2 = (ni_Nrest + q_Nrest).mag2(); //definetely check this equation
  Mx2 = (ni_Nrest + q_Nrest - Phad_Nrest ).mag2();// missing mass squared

  if(DEBUG)
    G4cout<<"W2: "<<W2<<" Mx2: "<<Mx2<< G4endl;

  return;
}

// is this only valid for SIDIS? 
void G4SBSTDISGen::GenerateHadronKinematics()
{
  if(DEBUG)
    G4cout<<"GenerateHadronKinematics: SIDIS"<<G4endl;

  if(DEBUG)
    G4cout<<"Hadron Energy: min:"<< tEMin_had << " max: "<<tEMax_had<< " proton mass: "<< Mh <<G4endl;
  
  hTheta = acos( CLHEP::RandFlat::shoot( cos( tThMax_had ), cos( tThMin_had ) ) );
  hPhi   =       CLHEP::RandFlat::shoot( tPhMin_had, tPhMax_had );
  hE     =       CLHEP::RandFlat::shoot( tEMin_had, tEMax_had );
  
  //because THIS condition, the minimum energy is the Mh, IN OUR CASE, Mp (always proton)
  if(DEBUG)
    G4cout<<"mass of the hadron: "<< Mh<< G4endl;
  
  hP = sqrt( pow(hE,2)-pow(Mh,2) );
  
  //unit vector in the direction of the scattered hadron in the LAB frame (this is just for clarity):
  G4LorentzVector  HadMom_local (hP*DirectionVector(hTheta, hPhi), hE );
  
  mom_had_final = HadMom_local; //I didn't find an easier way to do it (CA)
  
  if(DEBUG)
    G4cout<<"hTheta: "<<hTheta/deg<<" hPhi: "<<hPhi/deg<<" hE: "<< hE/GeV<<" hP: "<<hP<<" mom_had_final: "<<mom_had_final<<G4endl;
  
  return;
}
//IDEA, I can make the function G4LorentzVector and return the vector, instead of a global variable (same for the others)

void G4SBSTDISGen::SelectHadron()
{
  G4cout<<"tHadronType: "<<tHadronType<<G4endl;

  switch( tHadronType )
    {
    case kPiPlus:
      Mh = G4PionPlus::PionPlusDefinition()->GetPDGMass();
      icharge = 1;
      ihadron = 0;
      SIDISHadron = 0; //same order as sbstypes.hh
      break;
      
    case kP:
      Mh = G4Proton::ProtonDefinition()->GetPDGMass();
      icharge = 1;
      ihadron = 2;
      SIDISHadron = 5; //same order as sbstypes.hh
      break;
      
    default:
      Mh = G4Proton::ProtonDefinition()->GetPDGMass();
      icharge = 1;
      ihadron = 2;
      break;
    }


} 

void G4SBSTDISGen::GenerateInitialStateTDIS(Kine_t Kine)
{
 
  G4ThreeVector proton1 = FermiMomentum(); // it generates a vector with module Fermi momentum (in GeV)

  if (Kine == tTDISKinD)// The deuterium case
    {
      // This proton is the spectator in the Deuterium case BECAUSE the process is n->p+pi-
      // thus THIS one is measured IN COINCIDENCE with the proton calculated ahead
      // Note that they are always in lab frame --> proton (in D) has only fermi momentum 
      
      txa = tQ2/(2*Mn*tnu);

      //PROTON
      G4double P_p1     = proton1.mag();//momentum of the proton
      G4double E_p1     = sqrt (pow(P_p1,2) + pow(Mp,2));//energy of the proton
      
      iProton.set(E_p1, proton1);//

      if(DEBUG)
	G4cout<<"P_p1: " << P_p1<< " Mp: " << Mp<< " E_p1: "<< E_p1 <<G4endl;      

      //NEUTRON
      G4double P_n  = proton1.mag();
      G4double E_n  = m_deu - E_p1;
      
      iNeutron.set(E_n, -proton1); // neutron kinematics from proton
 
    }
  else // The hydrogen case
    {
      txa = tQ2/(2*Mp*tnu);

      iProton.set(Mp, null3D);

      iNeutron.set(0);
    }
  
  // if(DEBUG)
    G4cout<<"iProton: " << iProton  <<G4endl;
  if(DEBUG)
    G4cout<<"iNeutron: "<< iNeutron.e()<< " "<<iNeutron.e()/GeV<<G4endl;
  
  // New statement for TDIS
  Pisum_lab = tGetiElectron() + iProton;
  if(DEBUG)
    G4cout<<"Pisum_lab.m2() : "<<Pisum_lab.m2() << " Pisum_lab.e() : "<<Pisum_lab.e() <<G4endl;
  
//NOTE: all 4-mom vectors are in Lab system. 
  return;
}
 
void G4SBSTDISGen::GenerateFinalStateTDIS(Kine_t Kine,G4LorentzVector ei )
{
  G4LorentzVector q = q_lab;// calculated in PhotonKinematics
  G4cout<<"ei:"<< ei <<G4endl;
  G4cout<<"iProton: " << iProton  <<G4endl;   

  // The proton to be detected
  //===============
  //  pt  = CLHEP::RandFlat::shoot(0.0, 0.5*GeV);

  // I changed it to MeV in order to be consistent with the next calculation (CA)
  pt  = CLHEP::RandFlat::shoot(50.0, 500*MeV);// corrected from the proposal, was 0-500 MeV
  
  // I am assuming that this z is the same as the one in SIDIS
  // but here is randomize
  z   = CLHEP::RandFlat::shoot(0.0, 1.0); 

  if(DEBUG)
    G4cout<<"GenerateFinalStateTDIS z :"<< z <<G4endl;
  
  double znq;
  
  // xbj is defined by different initial nucleons, so it is 
  // different. y was also calculated there, and althought
  // it looks different, its value is practically the same
  // up to 5 decimals--> so y is recalculated here (try to simplify it)
  if(Kine ==  tTDISKinD)// The deuterium case
    { 
      xbj = tQ2/(2*iNeutron.dot(q));
      if(DEBUG)
	G4cout<<"xbj(TDIS): "<<xbj<<" tQ2: "<< tQ2<<" iNeutron.dot(q): "<<iNeutron.dot(q) <<G4endl;

      znq = z*iNeutron.dot(q);
      //      Mx2 = (q + iNeutron).mag2();
     
      W2  = (q + iNeutron).mag2();// DK stored fW2 = Mx2; which I think it is not correct (I keep both)
      tya = (iNeutron.dot(q))/(iNeutron.dot(ei));
    }
  else
    {
      xbj = tQ2/(2*iProton.dot(q));
      znq = z*iProton.dot(q);
      Mx2 = (q + iProton).mag2();
      W2  = (q + iProton).mag2(); // DK stored fW2 = Mx2; which I think it is not correct (I keep both)
      tya = (iProton.dot(q))/(iProton.dot(ei));
    }
  
  //  if(DEBUG)
    G4cout<<"xbj(TDIS): "<<xbj<<" tya: "<< tya<<" znq: "<<znq<<" Mx2: "<< Mx2 <<G4endl;
    G4cout<<"q: "<<q<<" iProton.dot(q): "<< iProton.dot(q)<<" iProton.dot(ei): "<<iProton.dot(ei) <<G4endl;
  
  // pt is randomize, so we need to calculate the z-component and then
  // calculate the whole 4-vector
  
  G4double fProton_P_z = (-1.0 * znq*q.z() + sqrt( pow(znq*q.z(),2) + tQ2*pow(tnu,2) * (pow(Mp,2) + pow(pt,2) ) - tQ2*pow(znq,2))) /tQ2;//z-component of the final proton
  if(DEBUG)
    G4cout<<"  z: "<< z<<" iProton: "<<iProton<<"\n q: "<< q << G4endl;
  if(DEBUG)
    G4cout<<"  znq: "<< znq<< " q.z(): "<<q.z()<< " tQ2: "<<tQ2<< " tnu: "<<tnu << " Mp: "<<Mp << " pt: "<<pt << G4endl;

  if(DEBUG)
    G4cout<<"  znq*q.z(): "<< znq*q.z()<< "  pow(znq*q.z(),2) + tQ2*pow(tnu,2) * (pow(Mp,2) + pow(pt,2): "<<  pow(znq*q.z(),2) + tQ2*pow(tnu,2) * (pow(Mp,2) + pow(pt,2))<< G4endl;

  if(DEBUG)
    G4cout<<" tQ2*pow(znq,2): "<< tQ2*pow(znq,2)<<G4endl;

  G4double fProton_P  = sqrt (pow(fProton_P_z, 2) + pow(pt,2));//momentum of the final proton
 
  G4double fProton_E = sqrt ( pow(fProton_P,2) + pow(Mp,2));
  
  G4double theta_p2 = acos (fProton_P_z/fProton_P);
  G4double phi_p2 = CLHEP::RandFlat::shoot( 0.0, 360.0*deg);

  fProton.set( fProton_E, fProton_P*DirectionVector(theta_p2, phi_p2));
  
  if(DEBUG)
    G4cout<<"fProton: "<< fProton<< " fProton_P: "<<fProton_P<< " energy (proton): "<< fProton.e()<<" "<<fProton.e()/GeV<<G4endl;


  // Hadron (PION IN THIS CASE IN PARTICULAR)
  //========
  if(Kine ==  tTDISKinD)// The deuterium case
    {
      Mx2 = (q + iNeutron - fProton).mag2(); //missing mass (CA)

      fPion.setE(iNeutron.e() - fProton.e());
      fPion.setPx(iNeutron.px() - fProton.px());
      fPion.setPy(iNeutron.py() - fProton.py());
      fPion.setPz(iNeutron.pz() - fProton.pz());
    }
  
  else
    {
      Mx2 = (q + iProton - fProton).mag2(); //missing mass (CA)

      fPion.setE(Mp - fProton.e());
      fPion.setPx(- fProton.px());
      fPion.setPy(- fProton.py());
      fPion.setPz(- fProton.pz());
    }


  Pfsum_lab = fProton + ef_lab; // should I include the pion? DD says no because it will be off-shell thus it could be anything
  if(DEBUG)
    G4cout<<"Pfsum_lab.m2() : "<<Pfsum_lab.m2() << " Pfsum_lab.e() : "<<Pfsum_lab.e() <<G4endl;
  
  //energy of pion sometimes negative, what does it means? it happens in g4sbs too
  if(DEBUG)
    G4cout<<"fPion: "<< fPion<< " energy (pion): "<< fPion.e()<<" tpi: "<< fPion.m2()/(GeV*GeV)<<" fProton.e(): "<< fProton.e() <<G4endl;
  // fPion defined by its components
  
  return;
}



G4double G4SBSTDISGen::GenerateTDIS(G4LorentzVector ei, Kine_t Kine)
{
  // Calculate the TDIS variables and Pion Structure Function
  // needed later to calculate the TDIS cross-section

  if(DEBUG)
    G4cout<<"GENERATE TDIS xbj: "<<xbj<<" z: "<<z<<" fPion: "<<fPion<<G4endl;

  G4LorentzVector q = q_lab;// calculated in PhotonKinematics

  // THESE QUANTITIES WERE CALCULATED IN g4sbs BUT ARE NOT USED
  G4ThreeVector p1f_vec = iProton.vect();
  G4ThreeVector p2f_vec = fProton.vect();
  G4ThreeVector q_vec = q.vect();
  
  double theta1 = acos( p1f_vec.unit().dot( q_vec.unit() ) );
  double theta2 = acos( p2f_vec.unit().dot( q_vec.unit() ) );
  // END NOT USED


  xpi = xbj/(1 - z);
  tpi = fPion.m2(); //altough called tpi it is simply the Mandelstan t. Here, it coincides with the pion 4-mom
  ypi = fPion.dot(q)/(fPion.dot(ei));

  
  G4double P_pi = fPion.vect().mag()/GeV;

  G4double P_recoil_proton = p2f_vec.mag();

  if(DEBUG)
    G4cout<<"P_pi: "<<P_pi<<" xbj: "<<xbj<<G4endl;
  
  // if (xbj > 0.055 && xbj < 0.3)
  //   {
  //     if ( Kine == tTDISKinD )// Deuterium case
  // 	fpi = 2*f2_pi(P_pi, xbj, theta2/deg); //According to Kijun, the fit was made in deg

  //     else // Hydrogen case
  // 	fpi = f2_pi(P_pi, xbj, theta2/deg);  //According to Kijun, the fit was made in deg
  //   }


  // I changed the momentum to send to the function. DK sent the pion, but
  // I believe should be the recoil proton (CA)

  if(DEBUG)
    G4cout<<"Recoil proton momentum: "<<P_recoil_proton<<" xbj: "<<xbj<<G4endl;



  
  if(! TH) //if TH code or not. In the preprocessor.
    {
      if (xbj > 0.055 && xbj < 0.3)
	{
	  if ( Kine == tTDISKinD )// Deuterium case
	    fpi = 2*f2_pi(P_recoil_proton, xbj, theta2/deg); //According to Kijun, the fit was made in deg

	  else // Hydrogen case
	    fpi = f2_pi(P_recoil_proton, xbj, theta2/deg);  //According to Kijun, the fit was made in deg
	}
      else
	{
	  if(DEBUG)
	    G4cout<<"Out of xbj range "<<xbj<<G4endl;
	  fpi = 0.0;
	}

    }
  else
    {
      if(DEBUG)
	G4cout<<" TH code"<< G4endl;

      G4double km  = P_recoil_proton/GeV; // in order to use similar nomenclature as TH (CA)
      G4double km1 = 0;
      G4double km2 = 0;

      kbinning(km, km1, km2); //function to obtain the bound limits of the momentum, parameters called by reference
    
      //    if(DEBUG)
	G4cout<<"Recoil proton momentum: "<<P_recoil_proton<<G4endl;
	//      if(DEBUG)
	G4cout<<"km: "<< km<<" km1: "<<km1<<" km2: "<<km2<<G4endl;
	//      if(DEBUG)
	G4cout << " xbj: "<<xbj<<"z: " << z << " 1-z: "<< 1-z<<" xpi: " << xpi <<G4endl;

      double f2_pi_TH;
      // I am using the Nucleon rest frame variables (global) 
      if (xbj > 0.055 && xbj < 0.3 &&
	  xpi > 0 && xpi < 1 &&
	  km < 0.4 ) // I am not sure if I should apply any other cut here 
	{
	  if ( Kine == tTDISKinD )// Deuterium case
	    {
	      f2pi_sub_(&Ebeam_Nrest, &xbj, &eth_Nrest,  &km1, &km2,  &fpi); 
	      //  fpi = 2*fpi; // the factor 2 was included in DK code, but I think that SHOULD be divided by 2

	      if(DEBUG)
		G4cout<<"f2pi from sub: "<<fpi<<G4endl;
	    }
	  else // Hydrogen case
	    {
	      f2pi_sub_(&Ebeam_Nrest, &xbj, &eth_Nrest, &km1, &km2, &fpi);
	    }
	}
      else
	{
	  if(DEBUG)
	    G4cout<<"Out of xbj range "<<xbj<<G4endl;
	  fpi = 0.0;
	}

    }


  return fpi;
}


void G4SBSTDISGen::kbinning(G4double km_b, G4double &km_b1, G4double &km_b2)
{
  // I am using this method to simplify any possible change in the k binning
  // the sub-index _b is just to remark that they belong just here --> (b)inning
  // km_b = hadron momentum 
  if(DEBUG)
    G4cout<<"***km_b: "<<km_b<<G4endl;
  if(DEBUG)
    G4cout<<"****km_b1: "<<km_b1<<" km_b2: "<<km_b2<<G4endl;

  // const int bins = 7;
  // G4double  k_min_arr[bins]={0.06, 0.100, 0.150, 0.200, 0.250, 0.300, 0.350}; //in GeV 
  // G4double  k_max_arr[bins]={0.100, 0.150, 0.200, 0.250, 0.300, 0.350, 0.400}; 
 
  const int bins = 100;
  G4double k_min_arr[bins];
  G4double k_max_arr[bins];

  G4double step = 0.3/(bins-1); //but it implies starting in 0.05GeV and not 0.06GeV

  for(int tt = 0; tt < bins; tt++)
    {
      k_min_arr[tt] = 0.050 + step*tt;
      k_max_arr[tt] = 0.100 + step*tt;


      if (km_b >= k_min_arr[tt] && km_b <= k_max_arr[tt]+0.00001)//with this we assure the uniqueness of the bin selection
	{
	  km_b1 = k_min_arr[tt];
	  km_b2 = k_max_arr[tt];
	}
    }


  // if (km_b >= k_min_arr[0] && km_b <= k_max_arr[0])
  //   {
  //     km_b1 = k_min_arr[0];
  //     km_b2 = k_max_arr[0];
  //   }

  // if (km_b > k_min_arr[1] && km_b <= k_max_arr[1])
  //   {
  //     km_b1 = k_min_arr[1];
  //     km_b2 = k_max_arr[1];
  //   }

  // if (km_b > k_min_arr[2] && km_b <= k_max_arr[2])
  //   {
  //     km_b1 = k_min_arr[2];
  //     km_b2 = k_max_arr[2];
  //   }

  // if (km_b > k_min_arr[3] && km_b <= k_max_arr[3])
  //   {
  //     km_b1 = k_min_arr[3];
  //     km_b2 = k_max_arr[3];
  //   }

  // if (km_b > k_min_arr[4] && km_b <= k_max_arr[4])
  //   {
  //     km_b1 = k_min_arr[4];
  //     km_b2 = k_max_arr[4];
  //   }

  // if (km_b > k_min_arr[5] && km_b <= k_max_arr[5])
  //   {
  //     km_b1 = k_min_arr[5];
  //     km_b2 = k_max_arr[5];
  //   }

  // if (km_b > k_min_arr[6] && km_b <= k_max_arr[6])
  //   {
  //     km_b1 = k_min_arr[6];
  //     km_b2 = k_max_arr[6];
  //   }

  //  if(DEBUG)
    G4cout<<"km_b: "<<km_b<<G4endl;
    //  if(DEBUG)
    G4cout<<"km_b1: "<<km_b1<<" km_b2: "<<km_b2<<G4endl;

}

G4double G4SBSTDISGen::F2_N(double x, double Q2, Kine_t Kine)
{
  // NOTE 1: this is the second time I use the CTEQ function, maybe a dedicated one?
  // NOTE 2: in the original, the Q2 dependency WAS not indicated by its unit (/GeV)
  // but I think it should be sent in GeV

  G4double qu = cteq_pdf_evolvepdf(__tdis_pdf, 1, x, sqrt(Q2)/GeV );
  G4double qd = cteq_pdf_evolvepdf(__tdis_pdf, 2, x, sqrt(Q2)/GeV );
  G4double qubar = cteq_pdf_evolvepdf(__tdis_pdf, -1, x, sqrt(Q2)/GeV );
  G4double qdbar = cteq_pdf_evolvepdf(__tdis_pdf, -2, x, sqrt(Q2)/GeV );
  
  if(DEBUG)
    G4cout<<"qu: "<<qu<<" qd: "<<qd<<" qubar: "<<qubar<<" qdbar: "<<qdbar<< G4endl;

  G4double quv = qu-qubar;
  G4double qdv = qd-qdbar;
  
  G4double qs = cteq_pdf_evolvepdf(__tdis_pdf, 3, x, sqrt(Q2)/GeV );
  
  G4double F2 = 0.0; 
  G4double e_u =  2.0/3.0;
  G4double e_d = -1.0/3.0;
  
  if( Kine == tTDISKinD) 
    {
      F2 += x*( e_u*e_u*qdv + e_d*e_d*quv );  //Deuterium case
    }
  else
    {
      F2 += x*( e_u*e_u*quv + e_d*e_d*qdv );  //Hydrogen case
    }

  // Sea quarks
  F2  += x*(2.0*e_u*e_u*qubar + 2.0*e_d*e_d*(qdbar + qs));
  
  return F2;
  
}


//TDIS
// subroutine to calculate the f2pi as function of recoiled nucleon momentum, xbj, theta (Timothy J. Hobbs)
// what is "th" = theta angle
//
// This is user parametrization by fit the Wally's codes 3Var_x.f() with integration of finite momentum range
// typ = 2 !EXPONENTIAL FORM FACTOR (s-depdendent exp)
// dis = 1 !NEUTRAL EXCHANGE
// FLAG = 0  --- THE PION CONTRIBUTION      | J = 0 + 1/2
double G4SBSTDISGen::f2_pi(double p, double x, double th)
{
  if(DEBUG)
    G4cout<<"ENTERING f2_pi"<<G4endl;

  double p0, p1, p2, p3, p4, p5;
  int xflag = 0;
  double fk = 0.0;
  double fth = 0.0;

  if (p > 0.05 && p <= 0.1){
    if(DEBUG)
      G4cout<<"[0,05, 0.1] "<<p<<G4endl;
    p0 = 3.656e-5;
    p1 = -0.000402;
    p2 = -0.008886;
    p3 = 0.07359;
    p4 = 1.079;
    p5 = -8.953;

    if (x < 0.0555 || x > 0.083) xflag = 1;

    
    fk = -0.954 + 66.5*p -1632.4*p*p + 14573*p*p*p; 
    
  }

  if (p > 0.1 && p <= 0.2){
    if(DEBUG)
      G4cout<<"[0,1, 0.2] "<<p<<G4endl;
    p0 = 0.000287;
    p1 = 0.009397;
    p2 = -0.2632;
    p3 = 2.029;
    p4 = -5.878;
    p5 = 4.664;

    if (x < 0.0555 || x > 0.16) xflag = 1;
    fk = 0.464 -15.4*p + 126.5*p*p;
    
  }
  
  if (p > 0.2 && p <= 0.3){
    if(DEBUG)
      G4cout<<"[0,2, 0.5] "<<p<<G4endl;
    p0 = 0.0003662;
    p1 = 0.02811;
    p2 = -0.4566;
    p3 = 2.428;
    p4 = -5.107;
    p5 = 3.222;

    if (x < 0.0555 || x > 0.226) xflag =1;

    fk = -1.133 + 8.5354*p;
  }

  if (p > 0.3 && p <= 0.5){
    if(DEBUG)
      G4cout<<"[0,3, 0.5] "<<p<<G4endl;
    p0 = 0.0009412;
    p1 = 0.01366;
    p2 = -0.1744;
    p3 = 0.3864;
    p4 = 0.6615;
    p5 = -2.113;
    
    if (x < 0.0555 || x > 0.281) xflag = 1;

    fk = -1.345 + 9.47*p -7.91*p*p;

  }

  if (p < 0.05 || p > 0.5){
    if(DEBUG)
      G4cout<<"0,05< p  >0.5 GeV"<< p <<G4endl;
    p0 = 0.0;
    p1 = 0.0;
    p2 = 0.0;
    p3 = 0.0;
    p4 = 0.0;
    p5 = 0.0;

  }

  if (th < 1.8 || th > 74)
    {
      if(DEBUG)
	G4cout<<"Out of angle range"<<G4endl;
      fth = 0.0;
    }
  else
    fth = -0.183 + 0.0976*th -0.0024*th*th + 0.000015*th*th*th; 

  if( xflag == 1 || x < 0.0555 || x > 0.3) 
    {
      if(DEBUG)
	G4cout<<"Out of xbj range (f2_pi): "<< xbj<<G4endl;
      return 0.0;
    }
  else
    {
      double f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4) + p5*pow(x,5);
      if(DEBUG)
	G4cout<<"Returning f2_pi "<<" f2: "<< f2<< " fk: "<< fk<< " fth: "<<fth<<" "<<f2*fk*fth <<G4endl;
      return f2*fk*fth;
    }
  
  
}



G4double G4SBSTDISGen::dissigma(G4double F2)
{
  G4double F1 = F2/(2.0*xbj);
  
  //NOTE: there was a HUGE mess in units here.

  // I was trying to use the internal system of units, but emulating the code from EventGen
  // was creating more trouble than solutions. 
  // I decided to do the same as G4SBSDIS, all the units in GeV and the final result convert in nanobarns/GeV^2

  // From PDG (eq. 18.8)
  G4double ds_dxdy = 4.0*PI*pow(alpha(),2) * ( (1.0 - tya - pow(xbj*tya*Mp/GeV,2) /tQ2/(GeV*GeV)) *F2 + pow(tya,2)*xbj*F1)
    /(xbj*tya*tQ2/(GeV*GeV));

  //tya was recalculated in GenerateFinalStateTDIS

  if(DEBUG)
    G4cout<<"ef_lab.e(): "<<ef_lab.e()<< " Mp: "<<Mp<<" tnu: "<<tnu<< G4endl;
  
  //this is ALSO from PDG (eq 18.1)
  G4double ds_dOmega_dE  = ds_dxdy* (ef_lab.e()/GeV)/(2.0 * PI * Mp/GeV * tnu/GeV);

  if(DEBUG)
    G4cout<<"ds_dOmega_dE (nb/GeV²): "<<ds_dOmega_dE*0.197*0.197*1e7<< G4endl;

  if(DEBUG)//this is a cross check that the units are ok when return the values
    G4cout<<"ds_dOmega_dE (cm²/GeV²): "<<ds_dOmega_dE*0.197*0.197*1e7/cm2<< G4endl;

  return ds_dOmega_dE *0.197*0.197*1e7; // GeV^-2 * nb

}




  
//maybe in the header? I like the functions in the code

G4ThreeVector G4SBSTDISGen::DirectionVector(G4double theta_dir, G4double phi_dir)
{
  return G4ThreeVector(sin(theta_dir)*cos(phi_dir), sin(theta_dir)*sin(phi_dir), cos(theta_dir));
}


  
// TDIS (copy and paste from EventGen but documented)
// hard coded the file into the program perhaps make it faster (?)
// specially since it only access two columns
G4ThreeVector G4SBSTDISGen::FermiMomentum()
{
  const int nfermi = 201;

  G4double CDF_fermi[nfermi]={0};// cummulative distribution function
  G4double f_mom_local[nfermi]={0};

  // This loop reads the data array in the header
  // is the equivalent to the external file readou in EventGen
  // BUT take care!!! we cannot use recursevely the variables
  // because it will change its value.

  // f_mom units are fm^-1 <--- thus we multiply by hbarc and the result is GeV

  for(int i = 0; i < nfermi; i++)
    {
      f_mom_local[i] = f_mom[i] * CLHEP::hbarc/(GeV*fermi); // hbarc in GeV fm. 
      // NOTE: it only gives the correct order of magnitude, but doesn't convert the "unit"
  
      CDF_fermi[i] = f_distr[i] * pow(f_mom_local[i], 2);
      
      if (i>0)//from the original code, I feel there are a couple of missing points, but it is not important
	CDF_fermi[i] =  CDF_fermi[i] + CDF_fermi[i-1];
     }



  if (CDF_fermi[nfermi-1] > 0.0)
    {
      for (int j = 0; j < nfermi; j++)
	{
	  CDF_fermi[j] = CDF_fermi[j]/CDF_fermi[nfermi-1];
	}
    }
  else
    {
      cout << "error in fermi file" << endl;
    }
  

  
  G4double fermi = -1.0*MeV;// test. We provide an unit to this variable
  G4double xran;

  // IS THIS LEGIT? I am just using another way to have a Fermi Momentum from a distribution. I could use flat 0-350MeV

  while (fermi < 0. || fermi > 0.250) // repeat until fermi is positive (always will be) OR the case bigger than 350MeV
{
  xran = CLHEP::RandFlat::shoot(0.0,1.0);

  // original fermi value is in GeV. I gave units to the variable, BUT its original value
  // came from GeV, thus I have to multiply by 1e3 to have them with the proper magnitude order.
  for (int i = 0; i < nfermi; i++)
    {
      if (CDF_fermi[i] == xran)
	{
	  fermi = f_mom[i];
	}
      else if (xran > CDF_fermi[i] && xran < CDF_fermi[i+1])
	{
	  G4double numerator   = ( f_mom[i+1] -  f_mom[i]);
	  G4double denominator = ( CDF_fermi[i+1] -   CDF_fermi[i]);
	  G4double slope = numerator/denominator;
	  fermi = f_mom[i] + (xran - CDF_fermi[i])*slope ;

	  if(DEBUG) 
	    G4cout<<"num: "<< numerator<< " denominator: "<< denominator<< " slope: "<< slope<<G4endl;
	  if(DEBUG)
	    G4cout<<"f_mom[i]: "<<  f_mom[i]<<" f_mom[i](MeV): "<<  f_mom[i]/GeV << " CDF_fermi[i]: "<< CDF_fermi[i]<<G4endl;
	  if(DEBUG)
	    G4cout<<"f_mom[i+1]: "<<  f_mom[i+1] << " CDF_fermi[i+1]: "<< CDF_fermi[i+1]<<G4endl;
    	}
    }
 

  }//close while loop

  if(DEBUG)
    cout << fermi << "\t" << xran << endl;
  
  G4double fth = acos( CLHEP::RandFlat::shoot( cos(180*deg), cos(0) ) );
  G4double fphi = CLHEP::RandFlat::shoot(0.0,360*deg);

  if(DEBUG)
    G4cout<<"Fermi th: "<< fth<<" "<<fth/deg<<"Fermi fphi: "<<fphi<<" "<<fphi/deg<<G4endl;

  // double fermix = fermi*sin(fth)*cos(fphi);
  // double fermiy = fermi*sin(fth)*sin(fphi);
  // double fermiz = fermi*cos(fth);
   
  return (fermi*1e3)*DirectionVector(fth, fphi); // Note that the units are GeV
  // I simplified and rewrite the whole function making use of CLHEP functions
  // but it is important to control that everything runs ok
}



void G4SBSTDISGen::GenerateSIDIS()
{
  //This is a very general approach, BUT WE ONLY WANT PROTONS

  if(DEBUG)
    G4cout << "(SIDIS) Hadron: " << tHadronType << " ihadron: " << ihadron << G4endl;

  G4double x = xbj; 

  //Get PDFs: sqrt(Q2) has units of energy, we should divide by GeV as argument to CTEQ: 
  double u = cteq_pdf_evolvepdf(__tdis_pdf, 1, x, sqrt(tQ2)/GeV );
  double d = cteq_pdf_evolvepdf(__tdis_pdf, 2, x, sqrt(tQ2)/GeV );
  double ubar = cteq_pdf_evolvepdf(__tdis_pdf, -1, x, sqrt(tQ2)/GeV );
  double dbar = cteq_pdf_evolvepdf(__tdis_pdf, -2, x, sqrt(tQ2)/GeV );
  double st = cteq_pdf_evolvepdf(__tdis_pdf, 3, x, sqrt(tQ2)/GeV );
  double sbar = st;
  
  //Gaussian model for transverse momentum widths of quark distribution (kperp) and fragmentation (pperp):
  double kperp2_avg = 0.25 * GeV * GeV;
  double pperp2_avg = 0.20 * GeV * GeV;
  
  
  tb = 1.0/( pow(z,2)*kperp2_avg + pperp2_avg );
  
  double e_u = 2.0/3.0;
  double e_d = -1.0/3.0;
  double e_s = -1.0/3.0;
  
  vector<double> Dqh;
  tFragFunc.GetFFs( ihadron, icharge, z, sqrt(tQ2)/GeV, Dqh ); 
  
  //Compute SIDIS structure function for a proton:
  H2 = x * tb/twopi*exp(-tb*pow(Ph_perp,2)) * ( pow(e_u,2) * (u * Dqh[0]  + ubar * Dqh[1]) + 
	   		  	                pow(e_d,2) * (d * Dqh[2]  + dbar * Dqh[3]) + 
					        pow(e_s,2) * (st * Dqh[4] + sbar * Dqh[5]) );
  
  H1 = H2/(2.0*x); //Callan-Gross relation

  if(DEBUG)
    G4cout << "Dqh[0]: " << Dqh[0] << " Dqh[1]: " << Dqh[1] << G4endl;
  if(DEBUG)
    G4cout << "Dqh[2]: " << Dqh[2] << " Dqh[3]: " << Dqh[3] << G4endl;

  if(DEBUG)
    G4cout<<"u: "<<u<<" d: "<<d<<" ubar: "<<ubar<<" dbar: "<<dbar<<" H2: "<< H2<< " H1: "<< H1<< G4endl;
  //some is already calculated

  etheta_Nrest = eth_Nrest; //it was calculated here, but it is the same definition, let's keep it like this for the moment

  theta_pq_Nrest = acos( Phad_Nrest_vect.unit().dot( q_Nrest_vect.unit() ) );

  if(DEBUG)
    G4cout<<"Phad_Nrest_vect.unit(): "<< Phad_Nrest_vect.unit()<< " q_Nrest_vect.unit(): "<< q_Nrest_vect.unit()<<G4endl;

  if(DEBUG)
    G4cout<<"theta_pq_Nrest: "<<   theta_pq_Nrest/deg << G4endl;




  return;
}



// I should put FULL dependences 
G4double G4SBSTDISGen::SIDISXS(G4LorentzVector ef_NrestL,G4double H1L,G4double H2L, G4double etheta_NrestL)
{
  if(DEBUG)
    G4cout<<"H1L: "<< H1L<<" H2L: "<<H2L<<G4endl;

  G4double sigsemi = 4.0 * pow(alpha(),2) * pow(hbarc, 2) * pow(ef_NrestL.e(),2)/pow(tQ2,2) * ( 2.0*H1L/proton_mass_c2 * pow(sin(etheta_Nrest/2.0),2) + H2L/nu_Nrest * pow(cos(etheta_NrestL/2.0),2) );
   
   
  //  double sigsemi = 4.0 * pow(fine_structure_const,2) * hbarc_squared * pow(ef_Nrest.e(),2)/pow(Q2,2) * ( 2.0*H1/proton_mass_c2 * pow(sin(etheta_Nrest/2.0),2) + H2/nu_Nrest * pow(cos(etheta_Nrest/2.0),2) );
   
   
  //This jacobian factor converts the cross section from d5sig/dE'dOmega_e dz dPh_perp^2 dphi_h  to 
  // d5sig/dE'dOmega_e dE_h dOmega_h. 
  double jacobian = 2.0 * Phad_Nrest_vect.mag2() * cos( theta_pq_Nrest ) / nu_Nrest;
  if(DEBUG)
    G4cout<<"Jacobian: "<<jacobian<<" FluxC: "<<FluxC<<G4endl;

  sigsemi *= jacobian; 

  //SIDIS uses the same flux correction as elastic, maybe is universal
  return sigsemi*FluxC;
}



G4double G4SBSTDISGen::FluxCorrection(G4LorentzVector ni, G4LorentzVector ei)
{
  //From G4SBS this correction is the same for the elastic and SIDIS case, I guess because is the electron
  //cross-section what it is under consideration


  G4double beta            = Phad_lab.vect().mag()/Phad_lab.e(); //beta = p/e
  G4double costheta_eN_lab = (ei.vect().unit() ).dot( ni.vect().unit() );
  G4double betaN_lab       = ni.beta();
  G4double gammaN_lab      = ni.gamma();
  G4double flux_Nrest      = 4.0*ni.m()*ei_Nrest.e();
  G4double flux_lab        = 4.0*ei.e()*ni.e()*sqrt( 2.0*(1.0-betaN_lab*costheta_eN_lab) - pow(gammaN_lab,-2) );

  //The lines above already converted the cross section to GEANT4 units. Now this has dimensions of area, is expressed in the lab frame, and is differential in solid angle only! (note from EventGen)

  return flux_Nrest/flux_lab;
}




G4double G4SBSTDISGen::ElasticXS(G4double beam_energy, G4double scatter_e_energy, G4double theta, Nucl_t nu)
{
  // I redefined the variables, just for facility notation
  // I prefer to keep large names at the beggining to clear know what are we talking about
  G4double E = beam_energy;
  G4double E_prime = scatter_e_energy;


  //Mott x Recoil fraction x form factors relation

  G4double eSigma = MottXS(theta, E)* (E_prime/E)*
    ( (pow(GE(nu),2)+tau()*pow(GM(nu),2)/(1.0+tau()) + 2.0*tau()*pow(GM(nu),2)*pow(tan(theta/2.0),2) )); 
  // Dimensions of area 
  
  return eSigma*FluxC;
}



G4double G4SBSTDISGen::PhotoD_XS(G4double E_photon)
{
  // E_photon; // virtual photon energy in MeV

  // Bethe-Peirls Deuterium photodesintegration cross-section
  // http://www1.jinr.ru/Pepan_letters/panl_2013_3/16_did.pdf

  if(E_photon < deu_bind)
    {
      return 0; //to extend the X axis to 0
    }
  else
    {

      // a is the constant factor of the Bethe-Peirls cross-section 
      // READ NOTE BELOW!!
      
      //  a = (8*PI)/(3*M);//h_bar,c, e = 1 
      //  a = (8*PI*pow(e,2)*h_bar)/(3*M*c);

      // b is the factor energy-dependent of the Bethe-Peirls cross-section 
      G4double b = (sqrt(deu_bind) * sqrt(pow(E_photon - deu_bind,3)))/pow(E_photon,3);
	
      return (2.4/0.056195)*b; // READ NOTE BELOW!!
    }

  // NOTE:
  // The paper indicates the maximum is 2.4mb @4.4MeV (which is true from  
  // real data base)
  // IN PRINCIPLE XS is in fm² (actually in MeV⁻² doing the
  // dimensional analysis with e=h=c=1) 
  // BUT considering the values of the constant units in many systems
  // there is no a consensus about their values.
  // Following a suggestion from Rey Cruz, I calculated the constant value
  // knowing the maximum. 2.4mb@4.4MeV. In other words /sigma = K * f(E)
  // where f(E) is the factor energy dependent (b here).
  // thus K = 2.4 mb / 0.056195 MeV^-1
}


G4double G4SBSTDISGen::VXPhoton_flux(double E_photon, double E_beam)
{
  // virtual photon flux, check with Raffo's CLAS12 J/psi photoproduction proposal

  G4double Q2max = 0.3; // (MeV/c)^2 Cut-off value. The formula is not so sensitive 
  // to this value. It was tested with 0.3, 3, 30 and 300 (MeV/c)²
  // and the change in flux is neglegible
  
  G4double x = E_photon / E_beam;

  if (x >= 1) return 0;
  
  G4double Q2min = m_e *x*x/(1-x);

  G4double term1 = (1 - x + x*x/2);
  G4double term2 = log(Q2max/Q2min);
  G4double term3 = alpha() /(E_beam*x*TMath::Pi());
  
  return term3 * ( term1 * term2 - (1-x)) ;
}



G4double G4SBSTDISGen::QuasiElasticXS(G4double beam_energy, G4int z1, G4int n1, G4int partID, G4double momentum, G4double angle)
{
  // This function calls a fortran function based in the 
  // EPC code by Lightbody and O'Connell, 
  // J.S. O'Connell and J.W.Lightbody, Computers in Physics 2,57(1988).
  // the function used here was modified by O. Rondon 
  // and adapted as a function by C. Ayerbe

  // the arguments needed (in order):
  // beam energy (in MeV)
  // atomic number
  // number of nucleons
  // kind of particle // 1=proton, -1=neutron, 2=pi+, -2=pi-, 0=pi0
  // momentum of the emmited particle
  // angle of emission (in deg, it is converted in rad inside the function)

  angle = angle*(180./TMath::Pi());

  if(DEBUG)
    G4cout<<"QuasiElasticXS: beam_energy: "<<beam_energy<< " momentum: "<<momentum<<" angle: "<<angle<<G4endl;

  // EPC notation (just neutron/proton now)
  if (partID == 0)
    {
      partID = -1; // neutron
    }
  else
    {
      partID = 1; //proton
    }

  return  epc_func_(&beam_energy, &z1, &n1, &partID, &momentum, &angle);

}


G4double G4SBSTDISGen::MottXS(G4double theta, G4double beam_energy)
{
  //the theta angle used is the theta in the NUCLEON REST frame!!

  G4double M1 = cos(theta/2)*alpha();

  G4double M2 = 2*beam_energy*pow(sin(theta/2),2);

  //MXS: Mott Cross Section. Keep in mind, you don't care units.
  G4double MXS = pow(hbarc, 2) * pow(M1/M2, 2);

  return MXS; 

  // According to Andrew, if we don't touch the units, Geant4 works internally with them
  // in other words, BECAUSE we give some units to certain numbers, GEANT4 knows
  // how to convert to the units WE WANT. 
  // For example, the final result here will be in mm^2, if we want barns
  // just MXS/barn and the final number is in such unit.
}


G4double G4SBSTDISGen::DipoleFF()
{
  return pow(1.0 + tQ2/(0.71*GeV*GeV), -2.0);
}


G4double G4SBSTDISGen::GE( Nucl_t nucleon)
{
  //I need to find where this comes from

  if (nucleon == 0) //a proton
    {
      iGE = (1.0-0.24*tau())/
	(1.0 + 10.98*tau() + 12.82*pow(tau(),2) + 21.97*pow(tau(),3));
    }
  else
    {
      iGE = (1.520*tau() + 2.629*pow(tau(),2) + 3.055*pow(tau(),2)*DipoleFF())/
	(1.0+5.222*tau()+0.040*pow(tau(),2)+11.438*pow(tau(),3));
    }
  
  return iGE;

}

G4double G4SBSTDISGen::GM( Nucl_t nucleon)
{
  if (nucleon == 0) //a proton
    {
      iGM = 2.79*(1.0+0.12*tau())/
	(1.0 + 10.97*tau() + 18.86*pow(tau(),2) + 6.55*pow(tau(),3) );
    }
  else
    {
      iGM = -1.913*(1.0+2.33*tau())/ 
	(1.0 + 14.72*tau() + 24.20*pow(tau(),2) + 84.1*pow(tau(),3));
    }
  
  return iGM;
  
}

G4double G4SBSTDISGen::tau()
{
  return tQ2/(4.0*Mp*Mp);
}


void G4SBSTDISGen::EventRejected()
{
  //gave null values to these variables when the event produced is rejected
  //fSigma = 0.0; // check what is the proper name
  

  //I was using nucleon, originally was hadron. Due to the generalization of SIDIS, should I go back to hadron?
  tNucleonE = 0.0;
  tNucleonP = G4ThreeVector();
  
  tElectronE = 0.0;
  tElectronP = G4ThreeVector();
  
  //      fWeight = 0.0; // ???
  xbj = -1.0; // change name 
  z = -1.0;
  
  return;
}





// THIS IS A TERRIBLE SOLUTION BUT I CAN'T FIND A WAY TO AVOID THE MULTIPLE DECLARATION ERROR

// This is the G4SBSDIS.hh code. I was unable to use here without conflict with G4SBSEventGen.

// changes (mainly addind a 't'):
//  __dis_pdf     --> __tdis_pdf;
//  initcteqpdf() --> tinitcteqpdf()
//  F2N           --> tF2N
//  dissigma      --> tdissigma
//  dissigma_p    -->  tdissigma_p
//  dissigma_n    -->  tdissigma_n



//cteq_pdf_t *__tdis_pdf;


void G4SBSTDISGen::tinitcteqpdf(){
  __tdis_pdf = cteq_pdf_alloc_id(400); // mode 400 = cteq6.6?

  assert(__tdis_pdf);
}


double G4SBSTDISGen::tF2N(double x, double Q2,  Nucl_t nucl){
  
  double qu = cteq_pdf_evolvepdf(__tdis_pdf, 1, x, sqrt(Q2) );
  double qd = cteq_pdf_evolvepdf(__tdis_pdf, 2, x, sqrt(Q2) );
  double qubar = cteq_pdf_evolvepdf(__tdis_pdf, -1, x, sqrt(Q2) );
  double qdbar = cteq_pdf_evolvepdf(__tdis_pdf, -2, x, sqrt(Q2) );

  double quv = qu-qubar;
  double qdv = qd-qdbar;

  double qs = cteq_pdf_evolvepdf(__tdis_pdf, 3, x, sqrt(Q2) );

  double F2 = 0.0; 
  double e_u =  2.0/3.0;
  double e_d = -1.0/3.0;

  if( nucl == kProton ){
    F2 += x*( e_u*e_u*quv + e_d*e_d*qdv ); 
  }
  if( nucl == kNeutron){
    F2 += x*( e_u*e_u*qdv + e_d*e_d*quv ); 
  }
  // Sea quarks
  F2  += x*(2.0*e_u*e_u*qubar + 2.0*e_d*e_d*(qdbar + qs));

  return F2;

}

// tdissigma is dissigma
// ebeam must be ei_Nrest.e()
// th should be eTh_Nrest  --> all thNR (theta NucleonRest) was th
// eprime is ef_Nrest.e()


double G4SBSTDISGen::tdissigma( double ebeam, double thNR, double eprime, Nucl_t nucl )
{
  // Return in nb/(GeV*sr) //nb in this file from Adikaram it is commented return in nb only - need to check TDIS xsec units
  
  double Q2 = 2.0*eprime*ebeam*(1.0-cos(thNR));
  //  double nu = ebeam-eprime;
  G4double nu = tnu;

  // every Mp value here changed by Mh
  //  double Mp = 0.938;
 
  //  double x = Q2/(2.0*Mp*nu);
  G4double x = xbj;

  double y = nu/ebeam;
  

  if( ! (0.0 < x && x < 1.0 && 0.0 < y && y < 1.0) )
    {
      printf("WARNING %s line %d  x = %f, y = %f -> eprime = %f GeV   th = %f deg  ebeam = %f GeV\n", __FILE__,
	     __LINE__, x, y, eprime, th*180/3.14159, ebeam );
      //	exit(1);
      return 0.0;;
    }
  
  double qu = cteq_pdf_evolvepdf(__tdis_pdf, 1, x, sqrt(Q2) );
  double qd = cteq_pdf_evolvepdf(__tdis_pdf, 2, x, sqrt(Q2) );
  double qubar = cteq_pdf_evolvepdf(__tdis_pdf, -1, x, sqrt(Q2) );
  double qdbar = cteq_pdf_evolvepdf(__tdis_pdf, -2, x, sqrt(Q2) );
  
  double quv = qu-qubar;
  double qdv = qd-qdbar;
  
  double qs = cteq_pdf_evolvepdf(__tdis_pdf, 3, x, sqrt(Q2) );
  
  double F2 = 0.0; 
  double e_u =  2.0/3.0;
  double e_d = -1.0/3.0;
  
  if( nucl == kProton ){
    F2 += x*( e_u*e_u*quv + e_d*e_d*qdv ); 
  }
  if( nucl == kNeutron){
    F2 += x*( e_u*e_u*qdv + e_d*e_d*quv ); 
  }
  // Sea quarks
  F2  += x*(2.0*e_u*e_u*qubar + 2.0*e_d*e_d*(qdbar + qs));
  double F1 = F2/(2.0*x);
  
  // From PDG
  double ds_dxdy = 4.0*3.14159*((1.0-y-pow(x*y*Mh,2.0)/Q2)*F2+y*y*x*F1)
			/(x*y*Q2*137.0*137.0);
  
  // In GeV^-2
  double ds_dOmega_dE = ds_dxdy*eprime/(2.0*3.14159*Mh*nu);
  
  return ds_dOmega_dE*0.197*0.197*1e7; // GeV2 -> nb
}
/*
double G4SBSTDISGen::tdissigma_p(double eb, double thNR, double ep){
  return tdissigma( eb, thNR, ep, kProton);
}
double G4SBSTDISGen::tdissigma_n(double eb, double thNR, double ep){
  return tdissigma( eb, thNR, ep, kNeutron );
}
*/



G4double G4SBSTDISGen::GetSigma(Kine_t kin)
{
  switch(kin){
  case (tElastic):
    return ELAsigma;
    break;
  case (tQuasiElastic):
    return QEsigma;
  break;
  case (tSIDIS):
    return SIDISsigma;
    break;
  case (kTDISGen):
  case (tTDISKinH):
  case (tTDISKinD):
    return TDISsigma;
  break;
  default:
    return QEsigma;
    break;
  }
}

