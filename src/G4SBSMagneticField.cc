#include "G4SBSMagneticField.hh"

#define MAXBUFF 1024

G4SBSMagneticField::G4SBSMagneticField(G4ThreeVector off, G4RotationMatrix rm) {
    fOffset = off;
    frm = rm;

    fInverted= false;
}


G4SBSMagneticField::~G4SBSMagneticField() {
}

//////////////////////////////////////////////////////////////////////////////////////////

G4SBSMappedField::G4SBSMappedField(G4ThreeVector off, G4RotationMatrix rm, const char *fn)
   : G4SBSMagneticField(off,rm) {

    strcpy( fFilename, fn);
    //fFieldVal = NULL;
    fBfield.clear();
    fN[0] = fN[1] = fN[2] = 0;
}


G4SBSMappedField::~G4SBSMappedField() {
}

