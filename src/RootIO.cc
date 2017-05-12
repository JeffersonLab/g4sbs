// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class RootIO
// Write energy, time, position and ID collected in a sebsitive detector
// element hit to a root TTree
// 05/05/14 JRMA

#include "RootIO.hh"
#include "G4RunManager.hh"
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "MCNtuple.h"
#include "TH3D.h"
#include "TH2D.h"
#include "TH1D.h"
//
//----------------------------------------------------------------------------
RootIO::RootIO(const char* fname, G4int type, const char* hname)
{
  f4Vector = NULL;
  // Root Input/Output
  fIevent = 1;
  if( type == EInputFile ){
    SetInput(fname);
    fIsInput = 1;
  }
  else if( type == E2DInputFile ){
    Set2DInput(fname,hname);
    fIsInput = 1;
  }
  else{
    SetOutput(fname);
    fIsInput = 0;
  }
  for(G4int i=0; i<8; i++) {
    fH3[i] = NULL;
    fH2[i] = NULL;
    fH1[i] = NULL;
  }
}
//
//----------------------------------------------------------------------------
RootIO::RootIO(G4int npart, G4int* parttype)
{
  // Root output...take parameters from PGA
  fNpart = npart;
  fPartType = parttype;
  f4Vector = new G4float*[npart];
  for(G4int i=0; i<npart; i++) f4Vector[i] = new G4float[5];
  fIsInput = 1;
}
//
//----------------------------------------------------------------------------
RootIO::~RootIO(){
  delete fidpart;
  delete fplab;
  delete felab;
  // for(Int_t i=0;i<fNpart;i++) delete fdircos[i];
  //if(fTreeO)delete fTree;
}
//
//-----------------------------------------------------------------------------
void RootIO::SetInput(const char* fn)
{
  // Setup data read from ROOT TTree held in TFile fName
  fFileI = new TFile(fn);
  if(!fFileI){
    printf("Error opening file %s\n",fn);
    exit(1);
  }
  //Get ntuple
  fTreeI = (TTree*)(fFileI->Get("h1"));
  if(!fTreeI){
    printf("Setup RootInput:no ntuple h1\n");
    fFileI->ls();
    exit(1);
  }
  fNevent = fTreeI->GetEntries();
  //Get the number of branches
  fNbranch = fTreeI->GetNbranches();
  // Calculate number of particles 
  // assuming 3 position branches and 5 particle variables 
  fNpart = (fNbranch - 3)/5;
  f4Vector=new Float_t*[fNpart];
  for(Int_t j=0;j<fNpart;j++) f4Vector[j]=new Float_t[5];
  fPartType=new Int_t[fNpart];
  fvertex = new G4float[3];
  //
  //Loop over branches and set addresses of variables attached to branches
  TObjArray* objarray = fTreeI->GetListOfBranches();
  for(Int_t i=0; i<fNbranch; i++){
    TBranch *branch = (TBranch*)(objarray->At(i));
    G4String bname = G4String((char*)(branch->GetName()));
    //Get vertex coordinates
    if(bname=="X_vtx")branch->SetAddress(&fvertex[0]);
    else if(bname=="Y_vtx")branch->SetAddress(&fvertex[1]);
    else if(bname=="Z_vtx")branch->SetAddress(&fvertex[2]);
    //Particle 4 vectors
    //First beam
    else if(bname=="Px_bm")branch->SetAddress(&f4Vector[0][0]);
    else if(bname=="Py_bm")branch->SetAddress(&f4Vector[0][1]);
    else if(bname=="Pz_bm")branch->SetAddress(&f4Vector[0][2]);
    else if(bname=="En_bm")branch->SetAddress(&f4Vector[0][3]);
    else if(bname=="Pt_bm")branch->SetAddress(&f4Vector[0][4]);
    //Then final state particles
    else if(!bname.contains("bm")&&!bname.contains("vtx")){
      char* sbname=(char*)(bname.data());
      // ID number of particle
      G4int index;
      sscanf(sbname+4,"%2d",&index);
      // Particle type...convert G3 to PDG
      G4int ind1;
      sscanf(sbname+6,"%d",&ind1);
      fPartType[index]= GetPDGfromG3(ind1); 
      if(bname.contains("Px"))branch->SetAddress(&f4Vector[index][0]);
      else if(bname.contains("Py"))branch->SetAddress(&f4Vector[index][1]);
      else if(bname.contains("Pz"))branch->SetAddress(&f4Vector[index][2]);
      else if(bname.contains("Pt"))branch->SetAddress(&f4Vector[index][4]);
      else if(bname.contains("En"))branch->SetAddress(&f4Vector[index][3]);
    }
  }
}
//
//-----------------------------------------------------------------------------
void RootIO::Set2DInput(const char* fn, const char* hn )
{
  // Setup data sampling from ROOT 2D hist held in TFile fName
  fFileI = new TFile(fn);
  if(!fFileI){
    printf("Error opening file %s\n",fn);
    exit(1);
  }
  fFileI->ls();
  //Get 2D hist
  //fFileI->cd();
  fH2s = (TH2D*)fFileI->Get(hn);
  if(!fH2s){
    printf("Setup RootInput from file %s:no 2D sampling histogram %s\n",fn,hn);
    fFileI->ls();
    exit(1);
  }
  G4int npart = 2;
  f4Vector = new G4float*[npart];
  for(G4int i=0; i<npart; i++) f4Vector[i] = new G4float[5];
}
//
//----------------------------------------------------------------------------
void RootIO::GetEvent()
{
  // Read single event from ROOT TTree
  fTreeI->GetEvent(fIevent);
  fIevent++;
}
//
//----------------------------------------------------------------------------
void RootIO::SetOutput(const char* fn)
{
  felab=new Float_t[fNpart]; 
  fplab=new Float_t[fNpart]; 
  fidpart=new Int_t[fNpart];

  //create Tree
  fTreeO=new TTree("ESSN-Out","Detector-response");
  fTreeO->SetAutoSave();
  fFileO = new TFile(fn,"RECREATE");
  fTreeO->SetDirectory(fFileO);

  //Collection IDs for decoding hit collections
  //  fCBCollID = -1;
  //  fTAPSCollID = -1;

  //Array stuff
  //    fArrayTot=fDET->GetNArraybars();
  fArrayTot = 16;
  fArr_i=new Int_t[fArrayTot];
  fArr_e=new Float_t[fArrayTot];
  fArr_ew=new Float_t[fArrayTot];
  fArr_t=new Float_t[fArrayTot];
  fArr_te=new Float_t[fArrayTot];
  fArr_x=new Float_t[fArrayTot];
  fArr_y=new Float_t[fArrayTot];
  fArr_z=new Float_t[fArrayTot];
  fStepTot = 1000;
  fStID = new Int_t[fStepTot];
  fStTrID = new Int_t[fStepTot];
  fStParID = new Int_t[fStepTot];
  fStStepID = new Int_t[fStepTot];
  fStdE = new Float_t[fStepTot];
  fStPx = new Float_t[fStepTot];
  fStPy = new Float_t[fStepTot];
  fStPz = new Float_t[fStepTot];
  fStPE = new Float_t[fStepTot];
  fStPreX = new Float_t[fStepTot];
  fStPreY = new Float_t[fStepTot];
  fStPreZ = new Float_t[fStepTot];
  fStPostX = new Float_t[fStepTot];
  fStPostY= new Float_t[fStepTot];
  fStPostZ= new Float_t[fStepTot];
  Int_t basket =64000;

  fTreeO->Branch("nhits",&fNhits,"fNhits/I",basket);
  fTreeO->Branch("npart",&fNpart,"fNpart/I",basket);
  fTreeO->Branch("plab",fplab,"fplab[fNpart]/F",basket);
  fTreeO->Branch("vertex",fvertex,"fvertex[3]/F",basket);
  fTreeO->Branch("beam",fBeam,"fBeam[5]/F",basket);
  fTreeO->Branch("dircos",fdircos,"fdircos[fNpart][3]/F",basket);
  fTreeO->Branch("elab",felab,"felab[fNpart]/F",basket);
  fTreeO->Branch("eleak",&feleak,"feleak/F",basket);
  fTreeO->Branch("etot",&fetot,"fetot/F",basket);
  fTreeO->Branch("Wgt",fWgt,"fWgt[8]/F",basket);
  // detector stuff
  fTreeO->Branch("arrN",&fNArr,"fNArr/I",basket);
  fTreeO->Branch("arrI",fArr_i,"fArr_i[fNArr]/I",basket);
  fTreeO->Branch("arrE",fArr_e,"fArr_e[fNArr]/F",basket);
  fTreeO->Branch("arrEw",fArr_ew,"fArr_ew[fNArr]/F",basket);
  fTreeO->Branch("arrT",fArr_t,"fArr_t[fNArr]/F",basket);
  fTreeO->Branch("arrTe",fArr_te,"fArr_te[fNArr]/F",basket);
  fTreeO->Branch("arrX",fArr_x,"fArr_x[fNArr]/F",basket);
  fTreeO->Branch("arrY",fArr_y,"fArr_y[fNArr]/F",basket);
  fTreeO->Branch("arrZ",fArr_z,"fArr_z[fNArr]/F",basket);
  fTreeO->Branch("Nst",&fNSt,"fNSt/I",basket);
  fTreeO->Branch("StID",fStID,"fStID[fNSt]/I",basket);
  fTreeO->Branch("StTrID",fStTrID,"fStTrID[fNSt]/I",basket);
  fTreeO->Branch("StParID",fStParID,"fStParID[fNSt]/I",basket);
  fTreeO->Branch("StStepID",fStStepID,"fStStepID[fNSt]/I",basket);
  fTreeO->Branch("StdE",fStdE,"fStdE[fNSt]/F",basket);
  fTreeO->Branch("StPx",fStPx,"fStPx[fNSt]/F",basket);
  fTreeO->Branch("StPy",fStPy,"fStPy[fNSt]/F",basket);
  fTreeO->Branch("StPz",fStPz,"fStPz[fNSt]/F",basket);
  fTreeO->Branch("StPE",fStPE,"fStPE[fNSt]/F",basket);
  fTreeO->Branch("StPreX",fStPreX,"fStPreX[fNSt]/F",basket);
  fTreeO->Branch("StPreY",fStPreY,"fStPreY[fNSt]/F",basket);
  fTreeO->Branch("StPreZ",fStPreZ,"fStPreZ[fNSt]/F",basket);
  fTreeO->Branch("StPostX",fStPostX,"fStPostX[fNSt]/F",basket);
  fTreeO->Branch("StPostY",fStPostY,"fStPostY[fNSt]/F",basket);
  fTreeO->Branch("StPostZ",fStPostZ,"fStPostZ[fNSt]/F",basket);

}
//
//---------------------------------------------------------------------------
void RootIO::WriteHit(G4HCofThisEvent* HitsColl)
{
  // Write hist information to TTree connected variables
  fNhits = fNArr = fNSt = 0;
  fNArr = *fNHits;
  *fNHits = 0;
  fetot = 0;
  G4int nst;
  //Get the ball hit info to be written to output
  ArrayHitsCollection* hc = (ArrayHitsCollection*)(HitsColl->GetHC(0));
  G4int nent = hc->entries();
  for(Int_t i=0; i<fNArr; i++){
    G4int j = fHits[i];
    if( j > nent ){
      G4cerr<<"RootIO::WriteHits, HitID out of range"<<G4endl;
      break;
    }
    ArrayHit* hit = (ArrayHit*)(hc->GetHit(j));
    fArr_e[i]=hit->GetEdep();
    fArr_ew[i]=hit->GetEdepB();
    fArr_t[i]=hit->GetTime();
    fArr_te[i]=hit->GetTimeE();
    fArr_x[i]=hit->GetPos().x();
    fArr_y[i]=hit->GetPos().y();
    fArr_z[i]=hit->GetPos().z();
    fArr_i[i]=hit->GetID();
    nst = hit->GetNSt();
    if( nst > fStepTot ) nst = 0;
    if( nst < 0 ) nst = 0;
    for(Int_t k=0; k<nst; k++){
      Int_t k1 = k + fNSt;
      if(k1 >= fStepTot){
	fNSt = fStepTot;
	break;
      }
      fStID[k1] = hit->GetID();
      fStTrID[k1] = hit->GetHitStep(k)->trID;
      fStParID[k1] = hit->GetHitStep(k)->parID;
      fStStepID[k1] = hit->GetHitStep(k)->stepID;
      fStdE[k1] = hit->GetHitStep(k)->dE;
      fStPx[k1] = hit->GetHitStep(k)->p4.x();
      fStPy[k1] = hit->GetHitStep(k)->p4.y();
      fStPz[k1] = hit->GetHitStep(k)->p4.z();
      fStPE[k1] = hit->GetHitStep(k)->p4.e();
      fStPreX[k1] = hit->GetHitStep(k)->pre.x();
      fStPreY[k1] = hit->GetHitStep(k)->pre.y();
      fStPreZ[k1] = hit->GetHitStep(k)->pre.z();
      fStPostX[k1] = hit->GetHitStep(k)->post.x();
      fStPostY[k1] = hit->GetHitStep(k)->post.y();
      fStPostZ[k1] = hit->GetHitStep(k)->post.z();
    }
    fNSt += nst;
    fHitID[j] = -1;
    fHits[i] = -1;
    hit->Clear();
  }
  WriteGenInput();
  fTreeO->Fill();
}
//
//-----------------------------------------------------------------------------
void RootIO::WriteGenInput(){
  //Note fvertex is already the pointer to fPGA::fPos 
  //Get the generated input info to be written to output
  for(G4int i=0; i<5; i++) { fBeam[i] = f4Vector[0][i]; }
  //Loop over the input particles and write their real kinematics
  for(Int_t i=1;i<fNpart;i++){
    G4int j = i-1;
    felab[j] = f4Vector[i][3];
    fplab[j] = f4Vector[i][4];
    fidpart[j]=fPartType[i];
    fdircos[j][0]=f4Vector[i][0];
    fdircos[j][1]=f4Vector[i][1];
    fdircos[j][2]=f4Vector[i][2];
  }
}
//
//-----------------------------------------------------------------------------
void RootIO::WriteH3(G4int ih3, G4double x, G4double y, G4double z, G4double v)
{
  fH3[ih3]->Fill(x,y,z,v);
}
//
//-----------------------------------------------------------------------------
void RootIO::WriteH2(G4int ih2, G4double x, G4double y, G4double v)
{
  fH2[ih2]->Fill(x,y,v);
}
//
//-----------------------------------------------------------------------------
void RootIO::WriteH1(G4int ih1, G4double x, G4double v)
{
  fH1[ih1]->Fill(x,v);
}
//
//-----------------------------------------------------------------------------
void RootIO::SetH3(G4int ih3, Hparm* hp)
{
  char name[32];
  sprintf(name,"H3-%d",ih3);
  fH3[ih3] = new TH3D(name,name,hp->nX,hp->Xmin,hp->Xmax,
		      hp->nY,hp->Ymin,hp->Ymax,hp->nZ,hp->Zmin,hp->Zmax);
}
//
//-----------------------------------------------------------------------------
void RootIO::SetH2(G4int ih2, Hparm* hp)
{
  char name[32];
  sprintf(name,"H2_%d",ih2);
  fH2[ih2] = new TH2D(name,name,hp->nX,hp->Xmin,hp->Xmax,
		      hp->nY,hp->Ymin,hp->Ymax);
}
//
//-----------------------------------------------------------------------------
void RootIO::SetH1(G4int ih1, Hparm* hp)
{
  char name[32];
  sprintf(name,"H1_%d",ih1);
  fH1[ih1] = new TH1D(name,name,hp->nX,hp->Xmin,hp->Xmax);
}
//
//-----------------------------------------------------------------------------
void  RootIO::Close()
{
  if(fTreeO) fTreeO->Write();
  for(G4int i=0; i<8; i++){
    if( fH3[i] ) fH3[i]->Write();
    if( fH2[i] ) fH2[i]->Write();
    if( fH1[i] ) fH1[i]->Write();
  }
  if(fFileO) fFileO->Close();
  // if(fFileI) fFileI->Close();
}
//
//-----------------------------------------------------------------------------
 void RootIO::Sample2D(G4double& p, G4double& costh){
   fH2s->GetRandom2(costh,p);
 }
