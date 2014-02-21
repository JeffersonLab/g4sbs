#include "G4SBSGlobalField.hh"
#include "G4SBSMagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"

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
























