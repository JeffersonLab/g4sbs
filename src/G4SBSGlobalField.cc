#include <TMD5.h>
#include "TCanvas.h"
#include "TH2F.h"
#include "TStyle.h"
#include "G4SBSGlobalField.hh"
#include "G4SBSMagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "sbstypes.hh"

#include "G4SBSToscaField.hh"
#include "G4SBSRun.hh"
#include <sys/stat.h>

#include <vector>

#define MAXBUFF 1024

G4SBSGlobalField::G4SBSGlobalField() {
  fInverted = false;
}


G4SBSGlobalField::~G4SBSGlobalField() {
}

void G4SBSGlobalField::GetFieldValue(const double Point[3],double *Bfield) const {

  unsigned int i;
  double Bfield_onemap[3];

  for( i = 0; i < 3; i++ ){ Bfield[i] = 0.0; }

  for (std::vector<G4SBSMagneticField *>::const_iterator it = fFields.begin() ; it != fFields.end(); it++){
    (*it)->GetFieldValue(Point, Bfield_onemap);
    for( i = 0; i < 3; i++ ){ Bfield[i] += Bfield_onemap[i]; }
  }


  return;
}


void G4SBSGlobalField::SetInvertField(G4bool b) {

  for (std::vector<G4SBSMagneticField *>::iterator it = fFields.begin() ; it != fFields.end(); it++){
    (*it)->InvertField(b);
  }

  return;
}



void G4SBSGlobalField::AddField( G4SBSMagneticField *f ){ 
  f->InvertField(fInverted);
  fFields.push_back(f); 

  // Rebuild chord finder now that the field sum has changed
  G4TransportationManager::GetTransportationManager()->GetFieldManager()->CreateChordFinder(this);

  return;
}


void G4SBSGlobalField::AddToscaField( const char *fn ){ 
  G4SBSToscaField *f = new G4SBSToscaField(fn);

  G4String fname = f->GetFilename();
  
  f->fArm = kHarm; //for now, a TOSCA field is always associated with "HARM". we may wish to change in the future.
    
  AddField(f);
  G4TransportationManager::GetTransportationManager()->GetFieldManager()->CreateChordFinder(this);

  G4SBSRunData *rd = G4SBSRun::GetRun()->GetData();
  
  TMD5 *md5 = TMD5::FileChecksum(fname.data());
  filedata_t fdata;

  strcpy(fdata.filename, fname.data() );
  strcpy(fdata.hashsum, md5->AsString() );

  G4cout << "MD5 checksum " << md5->AsString() << G4endl;

  delete md5;

  struct stat fs;
  stat(fname.data(), &fs);
  fdata.timestamp = TTimeStamp( fs.st_mtime );

  fdata.timestamp.Print();

  rd->AddMagData(fdata);


  return;
}

void G4SBSGlobalField::DropField( G4SBSMagneticField *f ){ 
  for (std::vector<G4SBSMagneticField *>::iterator it = fFields.begin(); it != fFields.end(); it++){
    if( (*it) == f ){ fFields.erase(it); }
  }
  G4TransportationManager::GetTransportationManager()->GetFieldManager()->CreateChordFinder(this);

  return;
}

void G4SBSGlobalField::ScaleFields( G4double scalefact, Arm_t arm ){
  for( std::vector<G4SBSMagneticField *>::iterator it = fFields.begin(); it != fFields.end(); ++it ){
    if( (*it)->fArm == arm ){
      (*it)->fScaleFactor = scalefact;
    }
  }
}

void G4SBSGlobalField::DebugField(G4double thEarm, G4double thHarm ){
  // New (added AJRP July 16, 2018): yz projections along the spectrometer axes (additional orientation check):

  G4ThreeVector Earm_zaxis( sin(thEarm), 0.0, cos(thEarm) );
  G4ThreeVector Earm_yaxis(0.0, 1.0, 0.0 );
  G4ThreeVector Earm_xaxis = Earm_yaxis.cross(Earm_zaxis).unit();

  G4ThreeVector Harm_zaxis( -sin(thHarm), 0.0, cos(thHarm) );
  G4ThreeVector Harm_yaxis(0.0, 1.0, 0.0 );
  G4ThreeVector Harm_xaxis = Harm_yaxis.cross(Harm_zaxis).unit();
  
  // Make a heatmap of the field strength in x-z plane for x direction
  int nstep = 201;

  double xmin = -4*m; double xmax =  4*m;
  double zmin = -1*m; double zmax =  7*m;

  TH2F *hx = new TH2F("field_x", "Field x component",
		      nstep, zmin/m, zmax/m, nstep, xmin/m, xmax/m );
  TH2F *hy = new TH2F("field_y", "Field y component",
		      nstep, zmin/m, zmax/m, nstep, xmin/m, xmax/m );
  TH2F *hz = new TH2F("field_z", "Field z component",
		      nstep, zmin/m, zmax/m, nstep, xmin/m, xmax/m );
  TH2F *h = new TH2F("field", "Field total magnitude",
		     nstep, zmin/m, zmax/m, nstep, xmin/m, xmax/m );

  TH2F *hBx_yzproj_Earm = new TH2F("hBx_yzproj_Earm","Bx, Earm yz projection",
				   nstep, zmin/m, zmax/m, nstep, xmin/m, xmax/m );
  TH2F *hBy_yzproj_Earm = new TH2F("hBy_yzproj_Earm","By, Earm yz projection",
				   nstep, zmin/m, zmax/m, nstep, xmin/m, xmax/m );
  TH2F *hBz_yzproj_Earm = new TH2F("hBz_yzproj_Earm","Bz, Earm yz projection",
				   nstep, zmin/m, zmax/m, nstep, xmin/m, xmax/m );
  TH2F *hBtot_yzproj_Earm = new TH2F("hBtot_yzproj_Earm","|B|, Earm yz projection", 
				     nstep, zmin/m, zmax/m, nstep, xmin/m, xmax/m );
  TH2F *hBx_yzproj_Harm = new TH2F("hBx_yzproj_Harm","Bx, Harm yz projection",
				   nstep, zmin/m, zmax/m, nstep, xmin/m, xmax/m );
  TH2F *hBy_yzproj_Harm = new TH2F("hBy_yzproj_Harm","By, Harm yz projection",
				   nstep, zmin/m, zmax/m, nstep, xmin/m, xmax/m );
  TH2F *hBz_yzproj_Harm = new TH2F("hBz_yzproj_Harm","Bz, Harm yz projection",
				   nstep, zmin/m, zmax/m, nstep, xmin/m, xmax/m );
  TH2F *hBtot_yzproj_Harm = new TH2F("hBtot_yzproj_Harm","|B|, Harm yz projection",
				     nstep, zmin/m, zmax/m, nstep, xmin/m, xmax/m );

  hx->GetXaxis()->SetTitle("z [m]");
  hx->GetXaxis()->CenterTitle();
  hx->GetYaxis()->SetTitle("x [m]");
  hx->GetYaxis()->CenterTitle();
  hy->GetXaxis()->SetTitle("z [m]");
  hy->GetXaxis()->CenterTitle();
  hy->GetYaxis()->SetTitle("x [m]");
  hy->GetYaxis()->CenterTitle();
  hz->GetXaxis()->SetTitle("z [m]");
  hz->GetXaxis()->CenterTitle();
  hz->GetYaxis()->SetTitle("x [m]");
  hz->GetYaxis()->CenterTitle();
  h->GetXaxis()->SetTitle("z [m]");
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->SetTitle("x [m]");
  h->GetYaxis()->CenterTitle();

  hBx_yzproj_Earm->GetXaxis()->SetTitle( "zBB (m)");
  hBx_yzproj_Earm->GetYaxis()->SetTitle( "yBB (m)");
  hBx_yzproj_Earm->GetXaxis()->CenterTitle();
  hBx_yzproj_Earm->GetYaxis()->CenterTitle();
  
  hBy_yzproj_Earm->GetXaxis()->SetTitle( "zBB (m)");
  hBy_yzproj_Earm->GetYaxis()->SetTitle( "yBB (m)");
  hBy_yzproj_Earm->GetXaxis()->CenterTitle();
  hBy_yzproj_Earm->GetYaxis()->CenterTitle();
  
  hBz_yzproj_Earm->GetXaxis()->SetTitle( "zBB (m)");
  hBz_yzproj_Earm->GetYaxis()->SetTitle( "yBB (m)");
  hBz_yzproj_Earm->GetXaxis()->CenterTitle();
  hBz_yzproj_Earm->GetYaxis()->CenterTitle();
  
  hBtot_yzproj_Earm->GetXaxis()->SetTitle( "zBB (m)");
  hBtot_yzproj_Earm->GetYaxis()->SetTitle( "yBB (m)");
  hBtot_yzproj_Earm->GetXaxis()->CenterTitle();
  hBtot_yzproj_Earm->GetYaxis()->CenterTitle();

  hBx_yzproj_Harm->GetXaxis()->SetTitle( "zSBS (m)");
  hBx_yzproj_Harm->GetYaxis()->SetTitle( "ySBS (m)");
  hBx_yzproj_Harm->GetXaxis()->CenterTitle();
  hBx_yzproj_Harm->GetYaxis()->CenterTitle();
  
  hBy_yzproj_Harm->GetXaxis()->SetTitle( "zSBS (m)");
  hBy_yzproj_Harm->GetYaxis()->SetTitle( "ySBS (m)");
  hBy_yzproj_Harm->GetXaxis()->CenterTitle();
  hBy_yzproj_Harm->GetYaxis()->CenterTitle();

  hBz_yzproj_Harm->GetXaxis()->SetTitle( "zSBS (m)");
  hBz_yzproj_Harm->GetYaxis()->SetTitle( "ySBS (m)");
  hBz_yzproj_Harm->GetXaxis()->CenterTitle();
  hBz_yzproj_Harm->GetYaxis()->CenterTitle();

  hBtot_yzproj_Harm->GetXaxis()->SetTitle( "zSBS (m)");
  hBtot_yzproj_Harm->GetYaxis()->SetTitle( "ySBS (m)");
  hBtot_yzproj_Harm->GetXaxis()->CenterTitle();
  hBtot_yzproj_Harm->GetYaxis()->CenterTitle();
  
  int i,j;

  double p[3];
  double B[3];

  G4ThreeVector ptemp;

  for( i = 0; i < nstep; i++ ){
    for( j = 0; j < nstep; j++ ){
      //i = "x" bin
      //j = "z" bin
      
      // G4double xtest = (xmax-xmin)*((double) i)/nstep + xmin;
      // G4double ytest = 0.0;
      // G4double ztest = (zmax-zmin)*((double) j)/nstep + zmin;

      G4double xtest = xmin + (i+0.5)*(xmax-xmin)/double(nstep);
      G4double ytest = 0.0;
      G4double ztest = zmin + (j+0.5)*(zmax-zmin)/double(nstep);
      
      p[1] = ytest;
      p[0] = xtest;
      p[2] = ztest;

      GetFieldValue(p,B);

      hx->Fill(ztest/m, xtest/m, B[0]/tesla);
      hy->Fill(ztest/m, xtest/m, B[1]/tesla);
      hz->Fill(ztest/m, xtest/m, B[2]/tesla);
      h->Fill(ztest/m, xtest/m, sqrt(B[0]*B[0]+B[1]*B[1]+B[2]*B[2])/tesla);

      //E arm yz projections:
      ptemp = xtest*Earm_yaxis + ztest*Earm_zaxis;

      p[0] = ptemp.getX();
      p[1] = ptemp.getY();
      p[2] = ptemp.getZ();

      GetFieldValue(p,B);

      hBx_yzproj_Earm->Fill( ztest/m, xtest/m, B[0]/tesla );
      hBy_yzproj_Earm->Fill( ztest/m, xtest/m, B[1]/tesla );
      hBz_yzproj_Earm->Fill( ztest/m, xtest/m, B[2]/tesla );
      hBtot_yzproj_Earm->Fill( ztest/m, xtest/m, sqrt(pow(B[0]/tesla,2)+pow(B[1]/tesla,2)+pow(B[2]/tesla,2)) );

      //H arm yz projections:
      ptemp = xtest*Harm_yaxis + ztest*Harm_zaxis;

      p[0] = ptemp.getX();
      p[1] = ptemp.getY();
      p[2] = ptemp.getZ();

      GetFieldValue(p,B);

      hBx_yzproj_Harm->Fill( ztest/m, xtest/m, B[0]/tesla );
      hBy_yzproj_Harm->Fill( ztest/m, xtest/m, B[1]/tesla );
      hBz_yzproj_Harm->Fill( ztest/m, xtest/m, B[2]/tesla );
      hBtot_yzproj_Harm->Fill( ztest/m, xtest/m, sqrt(pow(B[0]/tesla,2)+pow(B[1]/tesla,2)+pow(B[2]/tesla,2)) );
      
    }
  }

  fFieldPlots.push_back(hx);
  fFieldPlots.push_back(hy);
  fFieldPlots.push_back(hz);
  fFieldPlots.push_back(h);

  fFieldPlots.push_back( hBx_yzproj_Earm );
  fFieldPlots.push_back( hBy_yzproj_Earm );
  fFieldPlots.push_back( hBz_yzproj_Earm );
  fFieldPlots.push_back( hBtot_yzproj_Earm );

  fFieldPlots.push_back( hBx_yzproj_Harm );
  fFieldPlots.push_back( hBy_yzproj_Harm );
  fFieldPlots.push_back( hBz_yzproj_Harm );
  fFieldPlots.push_back( hBtot_yzproj_Harm );
}



















