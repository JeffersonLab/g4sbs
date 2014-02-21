#include "G4SBSGlobalField.hh"

#define MAXBUFF 1024

G4SBSGlobalField::G4SBSGlobalField(double zoffset, G4RotationMatrix *rm) {
    printf("Creating G4SBSGlobalField\n");
}


G4SBSGlobalField::~G4SBSGlobalField() {
}

void G4SBSGlobalField::GetFieldValue(const double Point[3],double *Bfield) const {

    int i;
    double Bfield_onemap[3];

    for( i = 0; i < 3; i++ ){ Bfield[i] = 0.0; }

     for (std::vector<G4MagneticField *>::iterator it = fFields.begin() ; it != fFields.end(); it++){
	 it->GetFieldValue(Point, Bfield_onemap);

	 for( i = 0; i < 3; i++ ){ Bfield[i] += Bfield_onemap[i]; }
     }


    return;
}




void G4SBSGlobalField::InvertField(bool b) const {

     for (std::vector<G4MagneticField *>::iterator it = fFields.begin() ; it != fFields.end(); it++){
	 it->InvertField(b);
     }

     return;
}



void G4SBSGlobalField::AddField( G4SBSMagneticField *f ){ 
    if( fInverted ){
	f->Invert(b);
    }

    fFields.push_back(f); 
    return;
}
























