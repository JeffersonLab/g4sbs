#include "G4SBSToscaField.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#define MAXBUFF 1024

G4SBSToscaField::G4SBSToscaField( G4String filename) 
  : G4SBSMappedField( G4ThreeVector(), G4RotationMatrix(),  filename ) 
{
  ReadField();    
}


G4SBSToscaField::~G4SBSToscaField() {
  // int i,j,k;

  // for( i = 0; i < fN[0]; i++ ){
  //   for( j = 0; j < fN[1]; j++ ){
  //     for( k = 0; k < fN[2]; k++ ){
  // 	delete fFieldVal[i][j][k];
  //     }
  //     delete fFieldVal[i][j];
  //   }
  //   delete fFieldVal[i];
  // }
  // delete fFieldVal;

  // fFieldVal = NULL;

  // return;
}

void G4SBSToscaField::GetFieldValue(const double Point[3],double *Bfield) const {
  double scale[3];
  int idx, jdx;
  int i, j, k;
  double sx, sy, sz;
  int bin[3];
  //double binwidth[3];
  double lowedge[3];
  
  //double point[3];
  double binfrac[3];

  //First of all: Is this the correct way to initialize the "local" coordinates of the point?
  //pt is a global point in the global coordinate system
  G4ThreeVector pt(Point[0], Point[1], Point[2]);
  pt = frm*pt - fOffset;

  //frm.print(G4cout);

  //As long as fOffset is defined in the local coordinate system relative to the 
  //origin, this is correct.

  // printf("Querying point %f %f %f\n", pt[0]/m, pt[1]/m, pt[2]/m);

  // Use a trilinear interpolation
  // Convert fed point to field map coordinates (180 degree rotation about y)
  // Overall - is from wrong polarity in field, primarily goes down +x-axis
  // which would deflect protons into the floor
  /*
    point[0] =  pt[0];
    point[1] = -pt[1];
    point[2] =  pt[2];
  */
  //point[0] = -pt[0];
  //point[1] = -pt[1];
  //point[2] = -pt[2];

  // Calculate index values for position

  //There is a smarter and more intuitive (and human-readable/understandable) way to accomplish this (but equally fast?)
  
  //scale[idx] = 3D position within the grid of this point, as a fraction of the grid extent along each dimension:
  for( idx = 0; idx < 3; idx++ ){
    scale[idx] = (pt[idx] - fMin[idx])/(fMax[idx] - fMin[idx]);

    //binwidth[idx] = (fMax[idx]-fMin[idx])/double( fN[idx]-1 );
    // Point is outside the volume of the grid: 
    if( scale[idx] < 0.0 || scale[idx] >= 1.0 ){
      // Out of range, return 0 field
      //	    printf("Out of range\n");
      for( jdx = 0; jdx < 3; jdx++ ){
	Bfield[jdx] = 0.0;
      }
      return;
    }
    
    bin[idx] = int( floor( (pt[idx]-fMin[idx])/fBinWidth[idx] ) );
    lowedge[idx] = fMin[idx] + bin[idx]*fBinWidth[idx]; //this calculation is probably unnecessary
    binfrac[idx] = (pt[idx]-lowedge[idx])/fBinWidth[idx];
  }

  sx = binfrac[0];
  sy = binfrac[1];
  sz = binfrac[2];

  i = bin[0];
  j = bin[1];
  k = bin[2];
  
  // i, j, k = largest integer value less than or equal to (number of grid points - 1) * fraction of full extent:
  // This is the bin of the grid into which this point falls
  // i = (int) floor( (fN[0]-1)*scale[0] );
  // j = (int) floor( (fN[1]-1)*scale[1] );
  // k = (int) floor( (fN[2]-1)*scale[2] );
 
  // lowedge[0] = floor( (pt[0] - fMin[0])/binwidth[0] );
  // j = floor( (pt[1] - fMin[1])/binwidth[1] );
  // k = floot( (pt[2] - fMin[2])/binwidth[2] );

  //
  // sx = (pt[0]/binwidth[0] - i);
  // sy = (pt[1]/binwidth[1] - j);
  // sz = (pt[2]/binwidth[2] - k);
  
  // // sx, sy, sz = fraction of bin extent within 3D bin where this point lies:
  // sx = (fN[0]-1)*scale[0] - (double) i;
  // sy = (fN[1]-1)*scale[1] - (double) j;
  // sz = (fN[2]-1)*scale[2] - (double) k;


  // Perform interpolation
  double interp[3];
  double c00, c10, c01, c11, c0, c1;

  for( idx = 0; idx < 3; idx++ ){ //Do trilinear interpolation within this 3D bin of each field component:
    

    c00 = fBfield[GetIndex(i,j,k)][idx]*(1.0-sx) + fBfield[GetIndex(i+1,j,k)][idx]*sx; //interpolate along X at low y, z, edges of bin
    c10 = fBfield[GetIndex(i,j+1,k)][idx]*(1.0-sx) + fBfield[GetIndex(i+1,j+1,k)][idx]*sx; //interpolate along X at high y, low z
    c01 = fBfield[GetIndex(i,j,k+1)][idx]*(1.0-sx) + fBfield[GetIndex(i+1,j,k+1)][idx]*sx; //interpolate along X at low y, high z
    c11 = fBfield[GetIndex(i,j+1,k+1)][idx]*(1.0-sx) + fBfield[GetIndex(i+1,j+1,k+1)][idx]*sx; //interpolate along x at high y, high z

    // c00 = fFieldVal[i][j][k][idx]*(1.0-sx) + fFieldVal[i+1][j][k][idx]*sx;
    // c10 = fFieldVal[i][j+1][k][idx]*(1.0-sx) + fFieldVal[i+1][j+1][k][idx]*sx;
    // c01 = fFieldVal[i][j][k+1][idx]*(1.0-sx) + fFieldVal[i+1][j][k+1][idx]*sx;
    // c11 = fFieldVal[i][j+1][k+1][idx]*(1.0-sx) + fFieldVal[i+1][j+1][k+1][idx]*sx;

    c0 = c00*(1.0-sy) + c10*sy; //Interpolate in Y low Z
    c1 = c01*(1.0-sy) + c11*sy; //Interpolate in Y at high Z

    interp[idx] = c0*(1.0-sz) + c1*sz; //Interpolate in Z
  }

  ///////////////////////////////////////////////
  // Make sure to put in local coordinates
  Bfield[0] =  interp[0];
  Bfield[1] =  interp[1];
  Bfield[2] =  interp[2];


  G4ThreeVector newB(Bfield[0], Bfield[1], Bfield[2]);

  //printf("%f %f %f -> %f %f %f\n", point[0]/cm, point[1]/cm, point[2]/cm, Bfield[0]/tesla, Bfield[1]/tesla, Bfield[2]/tesla );

  // Rotate to global coordinates
  newB = (frm.inverse())*newB;

  newB *= fScaleFactor;
  
  if( !fInverted ){
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

void G4SBSToscaField::ReadField(){
  fBfield.clear();//initialize map to empty to avoid issues.

  printf("G4SBSToscaField - Reading in field from %s\n", fFilename.data());
  FILE *f = fopen(fFilename.data(), "r");

  if( f == NULL ){ //try prepending field map path name:
    G4cout << "TOSCA field map not found..." << G4endl;
    
    const char *G4SBS_ENV_VAR = std::getenv("G4SBS");
    
    G4String prefix = G4SBS_ENV_VAR;
    prefix += "/share/fieldmaps/";
    
    fFilename.prepend( prefix );
    
    G4cout << "Trying " << fFilename << G4endl;
    
    f = fopen( fFilename.data(), "r" );
  
    if( !f ){
      fprintf(stderr, "Error: %s Line %d, %s - File %s could not be opened\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, fFilename.data()); 
      exit(1);
    }
  }

  double x,y,z, ang_x, ang_y, ang_z; //define angles relative to x, y and z axis for more flexibility in orientation of field

  // First line is the position offset
  fscanf(f, "%lf%lf%lf", &x, &y, &z);
  fOffset = G4ThreeVector(x*cm,y*cm,z*cm);
    
  // Second line is the rotation:
  // rotate around x, then y', then z'
    
  fscanf(f, "%lf%lf%lf", &ang_x, &ang_y, &ang_z );
  frm = G4RotationMatrix();
  frm.rotateX(ang_x*deg);
  frm.rotateY(ang_y*deg);
  frm.rotateZ(ang_z*deg);


  frm.print(G4cout);

  // Third line should have 4 values with the size of the
  // file indices
    
  int dint, idx, i,j,k;
  double r[3], fB[3];

  char dstring[MAXBUFF];

  int readorder[3];

  fscanf(f, "%d%d%d%d", &fN[2], &fN[1], &fN[0], &dint);
  fscanf(f, "%d%d%d", &readorder[2], &readorder[1], &readorder[0]);

  fBfield.resize( fN[2]*fN[1]*fN[0] );

  G4cout << "Bfield nx,ny,nz,ntotal = " << fN[0] << ", " 
	 << fN[1] << ", " 
	 << fN[2] << ", " << fBfield.size() << G4endl;

  /*

  // Ensure we have enough space to read this
  if( fN[0] > MAXPT || fN[1] > MAXPT || fN[2] > MAXPT  ){
  fprintf(stderr, "Error: %s Line %d, %s - File %s is too big\nRead parameters nx = %d ny = %d nz = %d, but MAXPT = %d\n", __FILE__, __LINE__, __PRETTY_FUNCTION__, fFilename, fN[0], fN[1], fN[2], MAXPT ); 
  exit(1);
  }
  */

  //fBfield.resize( fNx*fNy*fNz );

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

  dstring[0] = 'x';
  
  while( dstring[0] != '0' ){
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

  int effidx[3];
    
  //outermost index is x index, then y, then z:
  
  for( i = 0; i < fN[0]; i++ ){
    for( j = 0; j < fN[1]; j++ ){
      for( k = 0; k < fN[2]; k++ ){
	fscanf(f, "%lf%lf%lf%lf%lf%lf", &r[0], &r[1], &r[2], &fB[0], &fB[1], &fB[2]);

	//printf("%f %f %f %f %f %f\n", r[0], r[1], r[2], fB[0], fB[1], fB[2]);

	// Grab limits as we go along. 
	for( idx = 0; idx < 3; idx++ ){
	  if( r[idx]*cm < fMin[idx] ){ 
	    fMin[idx] = r[idx]*cm; 
	  }
	  if( r[idx]*cm > fMax[idx] ){ 
	    fMax[idx] = r[idx]*cm;
	  }
	}
	readorder[0] < 0 ? effidx[0] = fN[0]-i-1: effidx[0] = i;
	readorder[1] < 0 ? effidx[1] = fN[1]-j-1: effidx[1] = j;
	readorder[2] < 0 ? effidx[2] = fN[2]-k-1: effidx[2] = k;

	//fFieldVal[effidx[0]][effidx[1]][effidx[2]][idx] = fB[idx]*gauss;
	int index = effidx[2] + fN[2]*effidx[1] + fN[2]*fN[1]*effidx[0];
	
	fBfield[index] = G4ThreeVector( fB[0]*gauss, fB[1]*gauss, fB[2]*gauss );

	// G4cout << "Bx, By, Bz = " << fBfield[index][0] << ", "
	//        << fBfield[index][1] << ", "
	//        << fBfield[index][2] << G4endl;
      
      }
    }
  }

  for( idx=0; idx<3; idx++ ){
    fBinWidth[idx] = (fMax[idx]-fMin[idx])/double( fN[idx]-1 );
  }
  
  printf("G4SBSToscaField - Field complete\n");

  return;
}

int G4SBSToscaField::GetIndex( int ix, int iy, int iz ) const {
  return iz + fN[2]*iy + fN[2]*fN[1]*ix;
}
