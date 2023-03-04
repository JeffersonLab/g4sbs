#include "G4SBSBigBiteField.hh"
#include "G4SystemOfUnits.hh"

#define MAXBUFF 1024

G4SBSBigBiteField::G4SBSBigBiteField(G4ThreeVector offset, G4RotationMatrix rm, G4String filename) 
  : G4SBSMappedField( offset, rm, filename ) {
  ReadField();
}


G4SBSBigBiteField::~G4SBSBigBiteField() {
  // int i,j,k;
  // for( i = 0; i < fN[0]; i++ ){
  // 	for( j = 0; j < fN[1]; j++ ){
  // 	    for( k = 0; k < fN[2]; k++ ){
  // 		delete fFieldVal[i][j][k];
  // 	    }
  // 	    delete fFieldVal[i][j];
  // 	}
  // 	delete fFieldVal[i];
  // }
  // delete fFieldVal;
  // fFieldVal = NULL;
  // return;
}

void G4SBSBigBiteField::GetFieldValue(const G4double Point[3],G4double *Bfield) const {
  double scale[3];
  G4int idx, jdx;
  G4int i, j, k;
  double sx, sy, sz;

  double point[3];

  G4ThreeVector pt(Point[0], Point[1], Point[2]);
  pt = frm*pt - fOffset;

  //    printf("Querying point %f %f %f\n", pt[0]/m, pt[1]/m, pt[2]/m);

  // Use a trilinear interpolation
  // Convert fed point to field map coordinates (180 degree rotation about y)
  // Offset accounts for box positioning
  point[0] = -pt[0];
  point[1] =  pt[1];
  point[2] = -pt[2];

  // Calculate index values for position

  for( idx = 0; idx < 3; idx++ ){
    scale[idx] = (point[idx] - fMin[idx])/(fMax[idx] - fMin[idx]);

    if( (G4int) floor( (fN[idx]-1)*scale[idx] ) < 0 || 
	(G4int) floor( (fN[idx]-1)*scale[idx] ) >= fN[idx]-1 ){

      ///	if( s[idx] < 0.0 || s[idx] >= 1.0 ){
      // Out of range, return 0 field
      //	    printf("Out of range\n");
      for( jdx = 0; jdx < 3; jdx++ ){
	Bfield[jdx] = 0.0;
      }
      return;
    }

  }

  i = G4int(floor( (fN[0]-1)*scale[0] ) );
  j = G4int(floor( (fN[1]-1)*scale[1] ) );
  k = G4int(floor( (fN[2]-1)*scale[2] ) );

  sx = (fN[0]-1)*scale[0] - double(i);
  sy = (fN[1]-1)*scale[1] - double(j);
  sz = (fN[2]-1)*scale[2] - double(k);

  G4int ix = i, iy = j, iz = k;

  // Perform interpolation
  double interp[3];
  double c00, c10, c01, c11, c0, c1;

  for( idx = 0; idx < 3; idx++ ){
    // c00 = fFieldVal[i][j][k][idx]*(1.0-sx) + fFieldVal[i+1][j][k][idx]*sx;
    // c10 = fFieldVal[i][j+1][k][idx]*(1.0-sx) + fFieldVal[i+1][j+1][k][idx]*sx;
    // c01 = fFieldVal[i][j][k+1][idx]*(1.0-sx) + fFieldVal[i+1][j][k+1][idx]*sx;
    // c11 = fFieldVal[i][j+1][k+1][idx]*(1.0-sx) + fFieldVal[i+1][j+1][k+1][idx]*sx;
    
    c00 = fBfield[GetIndex(ix,iy,iz)][idx]*(1.0-sx) + fBfield[GetIndex(ix+1,iy,iz)][idx]*sx;
    c10 = fBfield[GetIndex(i,j+1,k)][idx]*(1.0-sx) + fBfield[GetIndex(i+1,j+1,k)][idx]*sx;
    c01 = fBfield[GetIndex(i,j,k+1)][idx]*(1.0-sx) + fBfield[GetIndex(i+1,j,k+1)][idx]*sx;
    c11 = fBfield[GetIndex(i,j+1,k+1)][idx]*(1.0-sx) + fBfield[GetIndex(i+1,j+1,k+1)][idx]*sx;
									     
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
  newB = (frm.inverse())*newB;

  newB *= fScaleFactor;

  if( !fInverted ){
    Bfield[0] = newB.x();
    Bfield[1] = newB.y();
    Bfield[2] = newB.z();
  } else {
    Bfield[0] = -newB.x();
    Bfield[1] = -newB.y();
    Bfield[2] = -newB.z();
  }
    

  //printf("Returning %f %f %f\n", Bfield[0], Bfield[1], Bfield[2]);

  return;
}

void G4SBSBigBiteField::ReadField(){
  printf("G4SBSBigBiteField - Reading in field\n");

  FILE *f = fopen(fFilename.data(), "r");

  if( f == NULL ){ //try prepending field map path name:
    G4cout << "BigBite field map not found..." << G4endl;
    
    const char *G4SBS_ENV_VAR = std::getenv("G4SBS");

    G4String prefix = G4SBS_ENV_VAR;
    prefix += "/share/fieldmaps/";

    //fFilename.prepend( prefix );
    fFilename.insert( 0, prefix );
    
    G4cout << "Trying " << fFilename << G4endl;

    f = fopen( fFilename.data(), "r" );

    if( !f ){
      fprintf(stderr, "Error: %s Line %d, %s - File %s could not be opened\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, fFilename.data() ); 
      exit(1);
    }
  }
  
  fBfield.clear();


  // First line should have 4 values with the size of the
  // file indices
    
  int dint, idx, i,j,k;
  double ddouble, x[3], fB[3];

  char dstring[MAXBUFF];

  fscanf(f, "%d%d%d%d", &fN[0], &fN[1], &fN[2], &dint);

  fBfield.resize( fN[2]*fN[1]*fN[0] );

  /*
  // Ensure we have enough space to read this
  if( fN[0] > MAXPT || fN[1] > MAXPT || fN[2] > MAXPT  ){
  fprintf(stderr, "Error: %s Line %d, %s - File %s is too big\nRead parameters nx = %d ny = %d nz = %d, but MAXPT = %d\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, fFilename, fN[0], fN[1], fN[2], MAXPT ); 
  exit(1);
  }
  */

  // Dynamically allocate table
  // fFieldVal = new double *** [fN[0]];
  // for( i = 0; i < fN[0]; i++ ){
  // 	fFieldVal[i] = new double ** [fN[1]];
  // 	for( j = 0; j < fN[1]; j++ ){
  // 	    fFieldVal[i][j] = new double * [fN[2]];
  // 	    for( k = 0; k < fN[2]; k++ ){
  // 		fFieldVal[i][j][k] = new double[3];
  // 	    }
  // 	}
  // }

  // Next 9 lines are not useful
  int nskip = 9; 

  for( idx = 0; idx < nskip; idx++ ){
    if( !fgets(dstring, MAXBUFF, f) ){
      fprintf(stderr, "Error: %s Line %d, %s - File %s has line too long (> %d)\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, fFilename.data(), MAXBUFF); 
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

	//printf("%f %f %f %f %f %f\n", x[0], x[1], x[2], fB[0], fB[1], fB[2]);

	// Grab limits as we go alone.  Assume this is a square grid
	for( idx = 0; idx < 3; idx++ ){
	  if( x[idx]*cm < fMin[idx] ){ 
	    fMin[idx] = x[idx]*cm; 
	  }
	  if( x[idx]*cm > fMax[idx] ){ 
	    fMax[idx] = x[idx]*cm;
	  }
	}
	//		    fFieldVal[i][j][k][idx] = fB[idx]*gauss;
	int index = k + fN[2]*j + fN[2]*fN[1]*i;

	fBfield[ index ] = G4ThreeVector( fB[0]*gauss, fB[1]*gauss, fB[2]*gauss );

      }
    }
  }

  printf("G4SBSBigBiteField - Field complete\n");

  return;
}



G4int G4SBSBigBiteField::GetIndex( G4int ix, G4int iy, G4int iz ) const {
  return iz + fN[2]*iy + fN[2]*fN[1]*ix;
}



























