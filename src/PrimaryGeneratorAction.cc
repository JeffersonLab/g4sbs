// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class PrimaryGeneratorAction
// Generation of particles
// 20/05/13 JRMA adapted from SBS equivalent, under construction
// 14/02/15 JRMA add event generation fron ROOT histogram

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4String.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "G4IonTable.hh"
#include "RootIO.hh"
#include "DIS.hh"
#include "Riostream.h"

PrimaryGeneratorAction::
PrimaryGeneratorAction(DetectorConstruction* dc)
{
  fDC = dc;
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  //create a messenger for this class
  fPGMessenger = new PrimaryGeneratorMessenger(this);
  //Get Particle table pointer
  fParticleTable=G4ParticleTable::GetParticleTable();
  fPDef=NULL;
  f4Vector=NULL;
  fPartType=NULL;
  fMass=NULL;
  fTrackThis=NULL;
 //default mode is g4 command line input
  fMode=EPGA_g4;
  fNevent=1;
  fNToBeTcount=0;
  fNTracked=0;
  fNpart=1;//for interactive use 
  fIpart = 0;
  fRe = fThe = 0.0;
  fName = fHname = NULL;
  fIsWindow = 0;        // default not calculating window effects
}



PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  delete fPGMessenger;
}



void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)

{
  //This function is called at the begining of event
  G4float Mass;
  G4float P;
  G4ThreeVector pvec;
  switch(fMode){
  case EPGA_g4:
    //standard g4 command line event generator
    fParticleGun->GeneratePrimaryVertex(anEvent);
    Mass=fParticleGun->GetParticleDefinition()->GetPDGMass();
    P=sqrt(fParticleGun->GetParticleEnergy()*fParticleGun->GetParticleEnergy()
	   - Mass*Mass);
    pvec=(fParticleGun->GetParticleMomentumDirection())*P;
    //fBeamLorentzVec->SetXYZM(pvec.x(),pvec.y(),pvec.z(),Mass);
    break;
  case EPGA_ROOT:
    if(fRootInput){
      //Get the event from input tree
      fRootInput->GetEvent();
      //Set vertex position
      fParticleGun->SetParticlePosition(G4ThreeVector(fPos[0]*cm,
						      fPos[1]*cm,
						      fPos[2]*cm));
      //Loop tracked particles, set particle gun, create vertex each particle
      for(G4int i=0;i<fNTracked;i++){
	fParticleGun->SetParticleDefinition(fPDef[fTrackThis[i]]);
	G4ThreeVector pdir = G4ThreeVector( f4Vector[fTrackThis[i]][0]*GeV,
					    f4Vector[fTrackThis[i]][1]*GeV,
					    f4Vector[fTrackThis[i]][2]*GeV
					    ).unit();
	fParticleGun->SetParticleMomentumDirection( pdir );
	fParticleGun->SetParticleEnergy(f4Vector[fTrackThis[i]][3]*GeV
					- fMass[fTrackThis[i]]);
	fParticleGun->GeneratePrimaryVertex(anEvent);
      }
    }
    else{ G4cerr<<"ROOT input mode specidied but no input file given"
		<<G4endl; exit(1);}
    break;
  case EPGA_multi:
    /*
    fSrcPos = fDC->GetSrcPos();
    fSrcSize = fDC->GetSrcSize();
    GenVertex();
    //Loop tracked particles, set particle gun, create vertex each particle
    for(G4int i=0;i<fNTracked;i++){
      fParticleGun->SetParticleDefinition(fPDef[i]);
      fParticleGun->SetParticleMomentumDirection( fPDir[i] );
      G4double energy = fPTmin[i] + fRandom->Uniform()*(fPTmax[i]-fPTmin[i]);
      fParticleGun->SetParticleEnergy( energy );
      fParticleGun->GeneratePrimaryVertex(anEvent);
    }
    */
    break;
  case EPGA_Ebeam:
    GenVertex();
    GenBrem(anEvent);
    break;
  case EPGA_2DSample:
    //GenVertex();
    Gen2DSample(anEvent);
    break;
  case EPGA_TDISp:
  case EPGA_TDISn:
    GenTDIS();
    TrackTDIS(anEvent);
    break;
  default:
    G4cerr<<"Unknown mode given to ESSNPrimiaryGeneratorAction"<<G4endl; 
    exit(1);

  }
}

void PrimaryGeneratorAction::SetUpROOTInput(G4String filename){
  // Check if Ebeam
  if(filename=="Ebeam"){
    SetUpEbeam();
    return;
  }
  else if(fMode==EPGA_2DSample){
    fName = (char*)filename.data();
    fRootInput = new RootIO(fName,E2DInputFile,fHname);
    return;
  }
  else if((fMode==EPGA_TDISp) || (fMode==EPGA_TDISn)){
    fName = (char*)filename.data();
    InitTDIS();
    //fRootInput = new RootIO(fName,E2DInputFile,fHname);
    return;
  }
  //Set ROOT input to default mode
  fMode=EPGA_ROOT;
  //Open ROOT file
  fName = (char*)filename.data();
  fRootInput = new RootIO(fName,EInputFile);
  //for(G4int i=0; i<3; i++){ fRootInput->SetH3(i); }
  fNpart = fRootInput->GetNpart();
  fPartType = fRootInput->GetPartType();
  f4Vector = fRootInput->Get4Vector();
  fPos = fRootInput->GetPos();
  fMass = new G4float[fNpart];
  if( !fPDef ) fPDef = new G4ParticleDefinition*[fNpart];
  for(G4int i=0;i<fNpart;i++){
    fPDef[i] =
	fParticleTable->FindParticle(fPartType[i]);
      // Try Ion table if not found
    if(!fPDef[i])	
      fPDef[i] = fParticleTable->GetIonTable()
	->GetIon(fPartType[i]);
    // Not found
    if(!fPDef[i]){
      printf("Undefined particle PDG ID = %d\n",fPartType[i]);
    }
    else{
      fMass[i]=fPDef[i]->GetPDGMass();
      printf("Particle PDG=%d, Mass=%g\n", fPartType[i],fMass[i]);
    }
  }
  if(fNTracked==0){//default track all
    fNTracked=fNpart-1; //-1 for beam
    fTrackThis=new G4int[fNTracked];
    for(G4int i=1;i<=fNTracked;i++) fTrackThis[i-1]=i;
  }
  for(G4int i=0; i<fNTracked; i++){
    G4int j = fTrackThis[i];
    printf("Tracking Particle PDG=%d, Mass=%g\n", fPartType[j],fMass[j]);
  }
}

void PrimaryGeneratorAction::SetUpEbeam()
{
  //  p = fParticleTable->GetIonTable()->GetIon(1000020040);
  fMode = EPGA_Ebeam;
  fPDef[0] = fParticleTable->FindParticle(11);
  fPDef[1] = fParticleTable->FindParticle(11);
  fMe = fPDef[0]->GetPDGMass();
  //  fEe = 0.855*GeV;
  fEe = fParticleGun->GetParticleEnergy();
  fPartType = new G4int[fNpart];
  fRootInput = new RootIO(fNpart,fPartType);
}

//-----------------------------------------------------------------------------
void PrimaryGeneratorAction::GenBrem(G4Event* anEvent)
{
  G4double E = fEe - fMe;
  G4double th = G4RandGauss::shoot(0.0,fThe*deg);
  G4double phi = G4UniformRand()*360*deg;
  G4double pz = E*cos(th);
  G4double px = E*sin(th)*cos(phi);
  G4double py = E*sin(th)*sin(phi);
  fPDir[0] = G4ThreeVector(px,py,pz);  // save g lab direction
  fParticleGun->SetParticleDefinition(fPDef[0]);
  fParticleGun->SetParticleMomentumDirection( fPDir[0] );
  fParticleGun->SetParticleEnergy(E);
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//-----------------------------------------------------------------------------
void PrimaryGeneratorAction::Gen2DSample(G4Event* anEvent)
{
  G4double p,costh,h;
  h = fDC->GetZst();
  h = -h + G4UniformRand()*2.0*h;
  if(fIsWindow){
    if( h < 0.0 ) h = -fDC->GetZst();
    else h = fDC->GetZst();
  }
  G4ThreeVector pos = G4ThreeVector(0.0,0.0,h);
  fParticleGun->SetParticlePosition(pos);
  fRootInput->Sample2D(p,costh);
  G4double sinth = sqrt(1.0 - costh*costh);
  G4double mp = fPDef[0]->GetPDGMass();
  G4double E = sqrt(p*p + mp*mp);
  //G4double th = G4RandGauss::shoot(0.0,fThe*deg);
  G4double phi = G4UniformRand()*360*deg;
  G4double pz = p*costh;
  G4double px = p*sinth*cos(phi);
  G4double py = p*sinth*sin(phi);
  fPDir[0] = G4ThreeVector(px,py,pz);  // save g lab direction
  fParticleGun->SetParticleDefinition(fPDef[0]);
  fParticleGun->SetParticleMomentumDirection( fPDir[0] );
  fParticleGun->SetParticleEnergy(E-mp);
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//-----------------------------------------------------------------------------
void PrimaryGeneratorAction::GenVertex()
{
  // Sample vertex position randomly over internal active part of source
  // cylinder dimensions obtained from detector construction class
  //
  G4double r = G4UniformRand()*fRe*mm;
  G4double phi = G4UniformRand()*360*deg;
  G4double h = -fDC->GetZext();
  G4ThreeVector pos = G4ThreeVector(r*cos(phi),r*sin(phi),h);
  fParticleGun->SetParticlePosition(pos);
}
/*
//-----------------------------------------------------------------------------
TH1D* PrimaryGeneratorAction::GenDsig(G4double E)
{
  // Sample vertex position randomly over internal active part of source
  // cylinder dimensions obtained from detector construction class
  //
  G4double m = fMe;
  G4double a = 0.625;
  G4double d = 27.0;
  G4double C = 9*a*a/(9 + d);
  G4double thC = m/E;
  TH1D* hdsig = new TH1D("BremAngle","GTagg",1000,0,3*thC);
  // angles out to 3 * char angle
  for(G4int i=0; i<1000; i++){
    G4double th = thC*3*i/1000;
    //G4double th = 0.25*thC*i/1000;
    G4double u = E*th/m;
    G4double f = C*(u*exp(-a*u) + d*u*exp(-3*a*u));
    hdsig->SetBinContent(i+1,f);
  }
  return hdsig;
}
*/
//-----------------------------------------------------------------------------
G4int PrimaryGeneratorAction::GetNEvents()
{
  if(fRootInput)
    return fRootInput->GetNevent();
  else
    return 0;
}  

//-----------------------------------------------------------------------------
void PrimaryGeneratorAction::InitTDIS()
{
  // Set up tagged deep inelastic scattering
  // init DIS cteq pdf
  initcteqpdf();
  // Get deuteron fermi momentum table
  InitDfermi();
  fMd     = 1.87561282*GeV;
  SetNTracked(6);
  fNpart = 6;
  fPDef[0] = fParticleTable->FindParticle(11);  //electron
  fMe = fPDef[0]->GetPDGMass();
  fEe = fParticleGun->GetParticleEnergy();
  fPDef[1] = fParticleTable->FindParticle(11);    // scattered electron
  fPDef[2] = fParticleTable->FindParticle(2212);  // struck proton
  fPDef[3] = fParticleTable->FindParticle(2212);  // spectator proton
  fPDef[4] = fParticleTable->FindParticle(2112);  // spectator neutron
  fPDef[5] = fParticleTable->FindParticle(211);   // pi+
  fMp = fPDef[2]->GetPDGMass();
  fMn = fPDef[4]->GetPDGMass();
  if(fMode == EPGA_TDISp){
    fTtype = 0;
    fNTracked = 1;
    fTrackThis[0] = 2;
  }
  else{
    fTtype = 1;
    fNTracked = 2;
    fTrackThis[0] = 3;
    fTrackThis[1] = 2;
  }
  fPartType = new G4int[fNpart];
  fRootInput = new RootIO(fNpart,fPartType);
  f4Vector = fRootInput->Get4Vector();
  fPtdis[0] = &fPei;
  fPtdis[1] = &fPef;
  fPtdis[2] = &fPp2f;
  fPtdis[3] = &fPp1f;
  fPtdis[4] = &fPnf;
  fPtdis[5] = &fPpif;
}
//-----------------------------------------------------------------------------
void PrimaryGeneratorAction::TrackTDIS(G4Event* anEvent)
{
  G4double h;
  h = fDC->GetZst();
  h = -h + G4UniformRand()*2.0*h;
  if(fIsWindow){
    if( h < 0.0 ) h = -fDC->GetZst();
    else h = fDC->GetZst();
  }
  G4ThreeVector pos(0.0,0.0,h);
  //Set vertex position
  fParticleGun->SetParticlePosition(pos);
  fRootInput->SetWgt(fSigDIS,0);
  fRootInput->SetWgt(fSigTDIS,1);
  for(G4int i=0;i<fNTracked;i++){
    G4int j = fTrackThis[i];
    SaveTrack(j);
    G4double mass = fPDef[j]->GetPDGMass();
    fParticleGun->SetParticleDefinition(fPDef[j]);
    G4ThreeVector pdir = G4ThreeVector( f4Vector[j][0],
					f4Vector[j][1],
					f4Vector[j][2]
					).unit();
    fParticleGun->SetParticleMomentumDirection( pdir );
    fParticleGun->SetParticleEnergy(f4Vector[j][3] - mass);
    fParticleGun->GeneratePrimaryVertex(anEvent);
  }
}

//-----------------------------------------------------------------------------
void PrimaryGeneratorAction::GenTDIS()
{
  // Sample vertex position randomly over internal active part of source
  // cylinder dimensions obtained from detector construction class
  //
  G4double pei = sqrt(fEe*fEe - fMe*fMe);
  fPei.set(0.0,0.0,pei,fEe);
  //The nucleon could be a proton or a neutron. It has initial 4-momentum ni:
  double Mt;
  G4LorentzVector* pf;
  G4LorentzVector q;
  if(fTtype){
    Mt = fMn;
    pf = &fPnf;
  }
  else{
    Mt = fMp;
    pf = &fPp1f;
  }
  double E_e, theta_e, P_e, phi_e;
  double Px_p1, Py_p1, Pz_p1, P_p1, E_p1;
  double Px_n, Py_n, Pz_n, E_n;
  double Px_p2, Py_p2, Pz_p2, P_p2, E_p2, theta_p2, phi_p2;
  double Px_pi, Py_pi, Pz_pi, E_pi;
  double pt,znq;
  //G4LorentzVector ef,q,p1f,nf,p2f;
  // Generate kinematics of the scattering electron and photon
  do{
    do{
      E_e      = URand( 1.0*GeV, 10.0*GeV);
      double cs  = URand( 0.96126, 0.99027);
      //theta_e  = URand( 5*deg, 45*deg);
      theta_e = acos(cs);
      P_e      = sqrt( E_e*E_e - fPei.mag2() );
      phi_e    = URand( -12*deg , 12*deg) ; 
      fPef.set( P_e*sin(theta_e)*cos(phi_e), P_e*sin(theta_e)*sin(phi_e),
	      P_e*cos(theta_e),E_e );
      q = fPei - fPef;     // virtual photon
      fQ2 = -q.mag2();
    }while( (fQ2 > 10*GeV*GeV) );
    fnu =  q.e();
    fxa =  fQ2/(2*Mt*fnu);
    fya =  fnu/fEe;
    if(!fTtype){
      E_p1  = fMp;   // proton at rest
      Px_p1 = 0.0;
      Py_p1 = 0.0;
      Pz_p1 = 0.0;
      E_n  = 0.0;   // no neutron
      Px_n = 0.0;
      Py_n = 0.0;
      Pz_n = 0.0;
    }
    else{
      // First Proton    
      GetDfermi();  // proton has fermi motion
      Px_p1  = fPfermi.x();
      Py_p1  = fPfermi.y();
      Pz_p1  = fPfermi.z();
      P_p1   = fPfermi.mag();
      E_p1   = sqrt (P_p1*P_p1 + fMp*fMp);
      // Neutron
      Px_n = -Px_p1;
      Py_n = -Py_p1;
      Pz_n = -Pz_p1;
      E_n  = fMd - E_p1;
    }    
    fPp1f.set( Px_p1, Py_p1, Pz_p1, E_p1);
    fPnf.set( Px_n, Py_n, Pz_n, E_n);
    // Second proton
    pt  = URand(0.0, 0.5*GeV);
    fz   = URand(0.0, 1.0);
    fxbj = fQ2/(2*pf->dot(q));
    znq = fz*pf->dot(q);
    fMx2 = (q + *pf).mag2();
    fy = (pf->dot(q))/(pf->dot(fPei));
  }while( fxbj > 1.0 );
  
  Pz_p2 = ( -1.0*znq*q.z() + sqrt( znq*q.z()*znq*q.z() +fQ2*(q.e()*q.e())*(fMp*fMp + pt*pt) - fQ2*znq*znq))/fQ2;
  P_p2  = sqrt (Pz_p2*Pz_p2 + pt*pt);
  theta_p2 = acos (Pz_p2/P_p2);
  E_p2 = sqrt ( P_p2*P_p2 + fMp*fMp);
  phi_p2 = URand( 0.0, 360.0*deg) ;
  Px_p2 = pt*cos(phi_p2);
  Py_p2 = pt*sin(phi_p2);
  fPp2f.set( Px_p2,Py_p2,Pz_p2,E_p2);

// Hadron 
//========
  if(fTtype){
    E_pi  = E_n - E_p2;
    Px_pi = Px_n - Px_p2; 
    Py_pi = Py_n - Py_p2;
    Pz_pi = Pz_n - Pz_p2;
  }
  else{
    E_pi  = fMp - E_p2;
    Px_pi = - Px_p2; 
    Py_pi = - Py_p2;
    Pz_pi = - Pz_p2;
  }
  double p_pi = sqrt(Px_pi*Px_pi + Py_pi*Py_pi + Pz_pi*Pz_pi);
  //G4ThreeVector pi_vec(Px_pi, Py_pi, Pz_pi);
  fPpif.set(Px_pi, Py_pi, Pz_pi, E_pi);
  fxpi = fxbj/(1 - fz);
  ftpi = (E_pi*E_pi) - (p_pi*p_pi);
  fypi = fPpif.dot(q)/(fPpif.dot(fPei));

  fSigDIS = fSigTDIS = fF2N = 0.0;
  if (fxbj > 0.055 && fxbj < 0.3){
    //nTDIS++;
    if( fTtype )
      fFpi = 2*F2pi(p_pi/GeV, fxbj, theta_p2/deg);
    else
      fFpi = F2pi(p_pi/GeV, fxbj, theta_p2/deg);
    //hF2p->Fill(fxbj,log10(fpi),p_pi);
    // From the DIS event generator
    //double minE = 0.0;
     //if (fxbj > 0 && fxbj < 1.0){
    // for the neutron/proton target
    fF2N   = F2N(fxbj, fQ2, fTtype);  // F2N is defined at G4SBSDIS.hh 
    if( fTtype )
      fSigDIS    = dissigma_n(fEe/GeV, theta_e/rad, E_e/GeV);
	//*((fEe-minE)/GeV);
    else
      fSigDIS    = dissigma_p(fEe/GeV, theta_e/rad, E_e/GeV);
	//*((fEe-minE)/GeV);
    
    fSigTDIS   = fSigDIS * (fFpi/fF2N);
  }
}

//-----------------------------------------------------------------------------
void PrimaryGeneratorAction::InitDfermi()
{ 
  const double fgev=0.1973;
  char line[256];
  // Read in fermi-momentum tabulation
  FILE* ffm = fopen("data/moment_ld2b.dat","r");
  if(!ffm){
    printf("Fatal error opening data/moment_ld2b.dat\n");
    exit(-1);
  }
  G4int i = 0; 
  for(;;){
    if(!fgets(line,255,ffm))
      break;
    if(sscanf(line,"%lf%lf",fPpfermi+i,fPfdis+i) != 2){
      printf("Fatal error reading deuteron fermi momentum file: %s\n",line);
      exit(-1);
    }
    fPpfermi[i] = fPpfermi[i]*fgev;
    fPfdis[i] = fPfdis[i]*fPpfermi[i]*fPpfermi[i];
    i++;
    if(i==nfermi) break;
  }
  if(i != nfermi){
    printf("Fatal error reading deuteron fermi momentum file\n");
    exit(-1);
  }
  for (G4int j = 1; j < nfermi; j++){
    fPfdis[j] = fPfdis[j] +  fPfdis[j-1];
  }
  for (G4int j = 0; j < nfermi; j++){
    fPfdis[j] = fPfdis[j]/fPfdis[nfermi-1];
  }
}

//-----------------------------------------------------------------------------
void PrimaryGeneratorAction::GetDfermi()
{  
  // randomly sample fermi momentum
  G4double xran = G4UniformRand();
  G4double pferm = -1.0;
  for (G4int i = 0; i < nfermi; i++){
    if (fPfdis[i] == xran){
      pferm = fPpfermi[i];
    }
    else if (xran > fPfdis[i] && xran < fPfdis[i+1]){
      G4double numerator   = ( fPpfermi[i+1] -  fPpfermi[i]);
      G4double denominator = (  fPfdis[i+1] -   fPfdis[i]);
      G4double slope = numerator/denominator;
      pferm = fPpfermi[i] + (xran - fPfdis[i])*slope;
      break;     
    }
  }
  pferm *= GeV;
  xran = G4UniformRand();
  G4double fth = acos(2*xran - 1.0);
  //G4double xran = G4UniformRand();
  G4double fphi = URand(0.0,2*PI);

  fPfermi.set( pferm*sin(fth)*cos(fphi), pferm*sin(fth)*sin(fphi),
	       pferm*cos(fth) );
}

//-----------------------------------------------------------------------------
G4double PrimaryGeneratorAction::F2pi(G4double p, G4double x, G4double th)
{
  // Empirical calculation of pion structure function
  // valid angle, x, p range
  if(th < 1.8 || th > 74)
    return 0.0;
  if((x < 0.0555) || (x > 0.3))
    return 0.0;
  if((p < 0.05) || (p > 0.5))
    return 0.0;
  //
  G4double p0, p1, p2, p3, p4, p5;
  G4double fk;
  G4double fth;
  if (p > 0.05 && p <= 0.1){
    if (x < 0.0555 || x > 0.083)
      return 0.0;
    p0 = 3.656e-5;
    p1 = -0.000402;
    p2 = -0.008886;
    p3 = 0.07359;
    p4 = 1.079;
    p5 = -8.953;
    fk = -0.954 + 66.5*p -1632.4*p*p + 14573*p*p*p;     
  }
  else if (p > 0.1 && p <= 0.2){
    if (x < 0.0555 || x > 0.16)
      return 0.0;
    p0 = 0.000287;
    p1 = 0.009397;
    p2 = -0.2632;
    p3 = 2.029;
    p4 = -5.878;
    p5 = 4.664;
    fk = 0.464 -15.4*p + 126.5*p*p;
  }
  else if (p > 0.2 && p <= 0.3){
    if (x < 0.0555 || x > 0.226)
      return 0.0;
    p0 = 0.0003662;
    p1 = 0.02811;
    p2 = -0.4566;
    p3 = 2.428;
    p4 = -5.107;
    p5 = 3.222;
    fk = -1.133 + 8.5354*p;
  }
  else if (p > 0.3 && p <= 0.5){ 
    if (x < 0.0555 || x > 0.281)
      return 0.0;
    p0 = 0.0009412;
    p1 = 0.01366;
    p2 = -0.1744;
    p3 = 0.3864;
    p4 = 0.6615;
    p5 = -2.113;
    fk = -1.345 + 9.47*p -7.91*p*p;
  }
  fth = -0.183 + 0.0976*th -0.0024*th*th + 0.000015*th*th*th; 
  G4double f2 = p0 + p1*pow(x,1) + p2*pow(x,2) + p3*pow(x,3) + p4*pow(x,4)
    + p5*pow(x,5);
  return f2*fk*fth;  
}
