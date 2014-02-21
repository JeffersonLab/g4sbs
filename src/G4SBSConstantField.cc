#include "G4SBSConstantField.hh"

#define MAXBUFF 1024

G4SBSConstantField::G4SBSConstantField(G4ThreeVector offset, G4RotationMatrix *rm, G4ThreeVector boxdim, G4ThreeVector fval) 
	: G4SBSMagneticField( offset, rm) {

	    fBoxDim   = boxdim;
	    fFieldVal = fval;


	    if( fBoxDim.x() < 0 ){
		fBoxDim.setX( -fBoxDim.x() );
		fprintf( stderr, "WARNING:  %s - x dimension of box negative\n", __PRETTY_FUNCTION__ );
	    }

	    if( fBoxDim.y() < 0 ){
		fBoxDim.setY( -fBoxDim.y() );
		fprintf( stderr, "WARNING:  %s - y dimension of box negative\n", __PRETTY_FUNCTION__ );
	    }

	    if( fBoxDim.z() < 0 ){
		fBoxDim.setZ( -fBoxDim.z() );
		fprintf( stderr, "WARNING:  %s - z dimension of box negative\n", __PRETTY_FUNCTION__ );
	    }

	    return;
}


G4SBSConstantField::~G4SBSConstantField() {
}

void G4SBSConstantField::GetFieldValue(const double Point[3],double *Bfield) const {
    G4ThreeVector pt(Point[0], Point[1], Point[2]);

    pt = ((*frm)*pt) - fOffset;

    int jdx;
    bool isinbox = true;

    if( fabs(pt.x()) > fBoxDim.x() || fabs(pt.y()) > fBoxDim.y() || fabs(pt.z()) > fBoxDim.z() ){
	isinbox = false;
    }

    if( !isinbox ){
	for( jdx = 0; jdx < 3; jdx++ ){
	    Bfield[jdx] = 0.0;
	}
	return;
    } 

    Bfield[0] = fFieldVal.x();
    Bfield[1] = fFieldVal.y();
    Bfield[2] = fFieldVal.z();

    G4ThreeVector newB(Bfield[0], Bfield[1], Bfield[2]);

    // Rotate to global coordinates
    newB = (frm->inverse())*newB;

    if( !fInverted ){
	Bfield[0] = newB.x();
	Bfield[1] = newB.y();
	Bfield[2] = newB.z();
    } else {
	Bfield[0] = -newB.x();
	Bfield[1] = -newB.y();
	Bfield[2] = -newB.z();
    }


    return;
}

