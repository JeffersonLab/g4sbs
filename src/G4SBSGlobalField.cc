#include <TMD5.h>
#include "TCanvas.h"
#include "TH2F.h"
#include "G4SBSGlobalField.hh"
#include "G4SBSMagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

#include "sbstypes.hh"

#include "G4SBSToscaField.hh"
#include "G4SBSRun.hh"
#include <sys/stat.h>

#include <vector>

#define MAXBUFF 1024

G4SBSGlobalField::G4SBSGlobalField() {
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
    AddField(f);
    G4TransportationManager::GetTransportationManager()->GetFieldManager()->CreateChordFinder(this);

    G4SBSRunData *rd = G4SBSRun::GetRun()->GetData();
    TMD5 *md5 = TMD5::FileChecksum(fn);
    filedata_t fdata;

    strcpy(fdata.filename, fn);
    strcpy(fdata.hashsum, md5->AsString() );

    G4cout << "MD5 checksum " << md5->AsString() << G4endl;

    delete md5;

    struct stat fs;
    stat(fn, &fs);
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


void G4SBSGlobalField::DebugField(){ 
    // Make a heatmap of the field strength in x-z plane for x direction
    int nstep = 100;

    double xmin = -4*m; double xmax =  4*m;
    double zmin = -1*m; double zmax =  7*m;

    TH2F *hx = new TH2F("heatmap_x", "heatmap x", nstep, xmin, xmax, nstep, zmin, zmax );
    TH2F *hz = new TH2F("heatmap_z", "heatmap z", nstep, xmin, xmax, nstep, zmin, zmax );

    int i,j;

    double p[3];
    double B[3];

    for( i = 0; i < nstep; i++ ){
	for( j = 0; j < nstep; j++ ){
	    p[1] = 0.0;

	    p[0] = (xmax-xmin)*((double) i)/nstep + xmin;
	    p[2] = (zmax-zmin)*((double) j)/nstep + zmin;

	    GetFieldValue(p,B);

	    hx->SetBinContent(i+1, j+2, B[0]/tesla);
	    hz->SetBinContent(i+1, j+2, B[2]/tesla);

	    if( B[0] != 0 ){
		printf("B = %f at %f %f\n", B[0]/tesla, p[0]/m, p[2]/m);
	    }
	}
    }

    TCanvas *c = new TCanvas();
    hx->Draw("COLZ");
    c->Print("heatmap_x.png");
    hz->Draw("COLZ");
    c->Print("heatmap_z.png");
}



















