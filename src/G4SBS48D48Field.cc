#include "G4SBS48D48Field.hh"

#define MAXBUFF 1024

G4SBS48D48Field::G4SBS48D48Field(double zoffset, G4RotationMatrix *rm) {
    printf("Creating G4SBS48D48Field\n");
    strcpy(fFilename, "GEN-map1.table");
    ReadField();

    // Intrinsic x offset to the map
    fZOffset = zoffset;
    frm = rm;

    // This map is not translatable or rotateable  (due to the inclusion of the beamline)
    // Check to make sure that the passed parameters are correct

    fUseGeantino = false;
}


G4SBS48D48Field::~G4SBS48D48Field() {
}

void G4SBS48D48Field::GetFieldValue(const double Point[3],double *Bfield) const {
    double s[3];
    int idx, jdx;
    int i, j, k;
    double sx, sy, sz;

    double point[3];

    G4ThreeVector pt(Point[0], Point[1], Point[2]);
    pt = ((*frm)*pt) - G4ThreeVector(0.0, 0.0, fZOffset - 48*2.54*cm/2.);

//    printf("Querying point %f %f %f\n", pt[0]/m, pt[1]/m, pt[2]/m);

    // Use a trilinear interpolation
    // Convert fed point to field map coordinates (180 degree rotation about y)
    // Offset accounts for box positioning
    point[0] = -pt[0];
    point[1] =  pt[1];
    point[2] = -pt[2];

    // Calculate index values for position

    for( idx = 0; idx < 3; idx++ ){
	s[idx] = (point[idx] - fMin[idx])/(fMax[idx] - fMin[idx]);


	if( s[idx] < 0.0 || s[idx] >= 1.0 ){
	    // Out of range, return 0 field
//	    printf("Out of range\n");
	    for( jdx = 0; jdx < 3; jdx++ ){
		Bfield[jdx] = 0.0;
	    }
	    return;
	}
    }

    i = (int) floor( (fN[0]-1)*s[0] );
    j = (int) floor( (fN[1]-1)*s[1] );
    k = (int) floor( (fN[2]-1)*s[2] );

    sx = (fN[0]-1)*s[0] - (double) i;
    sy = (fN[1]-1)*s[1] - (double) j;
    sz = (fN[2]-1)*s[2] - (double) k;


    // Perform interpolation
    double interp[3];
    double c00, c10, c01, c11, c0, c1;

    for( idx = 0; idx < 3; idx++ ){
	c00 = fFieldVal[i][j][k][idx]*(1.0-sx) + fFieldVal[i+1][j][k][idx]*sx;
	c10 = fFieldVal[i][j+1][k][idx]*(1.0-sx) + fFieldVal[i+1][j+1][k][idx]*sx;
	c01 = fFieldVal[i][j][k+1][idx]*(1.0-sx) + fFieldVal[i+1][j][k+1][idx]*sx;
	c11 = fFieldVal[i][j+1][k+1][idx]*(1.0-sx) + fFieldVal[i+1][j+1][k+1][idx]*sx;

	c0 = c00*(1.0-sy) + c10*sy;
	c1 = c01*(1.0-sy) + c11*sy;

	interp[idx] = c0*(1.0-sz) + c1*sz;
    }

    ///////////////////////////////////////////////
    // Make sure to put in local coordinates
    Bfield[0] =  interp[0];
    Bfield[1] = -interp[1];
    Bfield[2] =  interp[2];


    G4ThreeVector newB(Bfield[0], Bfield[1], Bfield[2]);

//    printf("%f %f %f -> %f %f %f\n", point[0]/cm, point[1]/cm, point[2]/cm, Bfield[0]/tesla, Bfield[1]/tesla, Bfield[2]/tesla );

    // Rotate to global coordinates
    newB = (frm->inverse())*newB;

    if( !fUseGeantino ){
	Bfield[0] = newB.x();
	Bfield[1] = newB.y();
	Bfield[2] = newB.z();
    } else {
	//  Flip sign of field if geantino
	//  since it's backwards and stupid
	Bfield[0] = -newB.x();
	Bfield[1] = -newB.y();
	Bfield[2] = -newB.z();
    }
    

    //printf("Returning %f %f %f\n", Bfield[0], Bfield[1], Bfield[2]);

    return;
}

void G4SBS48D48Field::ReadField(){
    printf("G4SBS48D48Field - Reading in field\n");
    FILE *f = fopen(fFilename, "r");

    if( !f ){
	fprintf(stderr, "Error: %s Line %d, %s - File %s could not be opened\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, fFilename); 
	exit(1);
    }


    // First line should have 4 values with the size of the
    // file indices
    
    int dint, idx, i,j,k;
    double ddouble, x[3], fB[3];

    char dstring[MAXBUFF];

    fscanf(f, "%d%d%d%d", &fN[2], &fN[1], &fN[0], &dint);


    // Ensure we have enough space to read this
    if( fN[0] > MAXPT || fN[1] > MAXPT || fN[2] > MAXPT  ){
	fprintf(stderr, "Error: %s Line %d, %s - File %s is too big\nRead parameters nx = %d ny = %d nz = %d, but MAXPT = %d\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, fFilename, fN[0], fN[1], fN[2], MAXPT ); 
	exit(1);
    }


    // Next 8 lines are not useful
    int nskip = 8; 

    for( idx = 0; idx < nskip; idx++ ){
	if( !fgets(dstring, MAXBUFF, f) ){
	    fprintf(stderr, "Error: %s Line %d, %s - File %s has line too long (> %d)\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, fFilename, MAXBUFF); 
	    exit(1);
	}
    }

    for( idx = 0; idx < 3; idx++ ){
	fMin[idx] =  1e9;
	fMax[idx] = -1e9;
    }

    // Read in the field map

    for( i = 0; i < fN[0]; i++ ){
	for( j = 0; j < fN[1]; j++ ){
	    for( k = 0; k < fN[2]; k++ ){
		fscanf(f, "%lf%lf%lf%lf%lf%lf%lf", &x[0], &x[1], &x[2], &fB[0], &fB[1], &fB[2], &ddouble );

	//	printf("%f %f %f %f %f %f\n", x[0], x[1], x[2], fB[0], fB[1], fB[2]);

		// Grab limits as we go alone.  Assume this is a square grid
		for( idx = 0; idx < 3; idx++ ){
		    if( x[idx]*cm < fMin[idx] ){ 
			fMin[idx] = x[idx]*cm; 
		    }
		    if( x[idx]*cm > fMax[idx] ){ 
			fMax[idx] = x[idx]*cm;
		    }

		    fFieldVal[i][j][k][idx] = fB[idx]*gauss;
		}

	    }
	}
    }

    printf("G4SBS48D48Field - Field complete\n");

    return;
}

bool G4SBS48D48Field::GoodParams() {

    if( fabs(fZOffset - 2.8*m)>1*mm ) return false;
    // FIXME:  This needs a rotation angle check
    
    return true;
}





























