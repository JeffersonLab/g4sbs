// SBS a Geant-4 Based Model of Hall-A Experiments with 11 GeV
// J.R.M Annand, University of Glasgow
// Class TOSCAField3D
// TOSCA generated field map
// Original by D.J.Hamilton
// 23/02/12 JRMA magnet "aperture" dimension changeable,
//               e.g. where field extends beyond aperture
// 27/02/12 JRMA implement x,y,z-offset to coordinates..covert to cm
//               change to c-style printf, scanf i/o
//

#include "TOSCAField3D.hh"

TOSCAField3D::TOSCAField3D( const char* filename, double xoff, double yoff,
			    double zoff, double scale) 
  :invertX(false),invertY(false),invertZ(false)
{    
  double lenUnit   = cm;
  double fieldUnit = gauss; 
  fXoffset = xoff/lenUnit;
  fYoffset = yoff/lenUnit;
  fZoffset = zoff/lenUnit;
  printf(" Reading TOSCA generated field map from file %s \n", filename);
  // Open the file for reading.
  FILE* file = fopen(filename,"r");
  if(!file){
    printf(" Error opening file %s\n Abort Geant-4 session\n",filename);
    exit(-1);
  }
  char* fl;
  // Ignore first blank line, if the 1st is blank
  // and get the table dimensions
  char buffer[256];
  for(;;){
    fl = fgets(buffer,256,file);
    if(sscanf(buffer,"%d%d%d",&nx,&ny,&nz) == 3) break;
  }
  printf(" X-dimension: %d,  Y-dimension: %d,  Z-dimension: %d\n",nx,ny,nz);
  // Set up storage space for table
  xField.resize( nx );
  yField.resize( nx );
  zField.resize( nx );
  G4int ix, iy, iz;
  for (ix=0; ix<nx; ix++) {
    xField[ix].resize(ny);
    yField[ix].resize(ny);
    zField[ix].resize(ny);
    for (iy=0; iy<ny; iy++) {
      xField[ix][iy].resize(nz);
      yField[ix][iy].resize(nz);
      zField[ix][iy].resize(nz);
    }
  }
  // Skip other header information    
  // The first line whose 1st integer is 0 is considered to
  // be the last line of the header.
  for(;;){
    fl = fgets(buffer,256,file);
    G4int xx;
    sscanf(buffer,"%d",&xx);
    if( xx == 0 ) break;
  }
  // Read in the data
  double xval,yval,zval,bx,by,bz;
  double permeability;                 // Not used here
  for (ix=0; ix<nx; ix++) {
    for (iy=0; iy<ny; iy++) {
      for (iz=0; iz<nz; iz++) {
	if(!(fl = fgets(buffer,256,file)) ){
	  printf("Premature end of file %s\n",filename);
	  fclose(file);
	  exit(-1);
	}
	if( sscanf(buffer,"%lf%lf%lf%lf%lf%lf%lf",
		   &xval,&yval,&zval,&bx,&by,&bz,&permeability) < 6 ){
	  printf("File format error %s, exiting\n",filename);
	  fclose(file);
	  exit(-1);
	}
	bx = bx * scale;
	by = by * scale;
	bz = bz * scale;
	//file >> xval >> yval >> zval >> bx >> by >> bz >> permeability;
	xval += fXoffset;
	yval += fYoffset;
	zval += fZoffset;
        if ( ix==0 && iy==0 && iz==0 ) {
          minx = xval * lenUnit;
          miny = yval * lenUnit;
          minz = zval * lenUnit;
        }
        xField[ix][iy][iz] = bx * fieldUnit;
        yField[ix][iy][iz] = by * fieldUnit;
        zField[ix][iy][iz] = bz * fieldUnit;
      }
    }
  }
  fclose(file);

  maxx = xval * lenUnit;
  maxy = yval * lenUnit;
  maxz = zval * lenUnit;

  // G4cout << " Read values of field from file " << filename << endl; 
  G4cout << " ---> assumed the order:  x, y, z, Bx, By, Bz "
	 << "\n ---> Min values x,y,z: " 
	 << minx/cm << " " << miny/cm << " " << minz/cm << " cm "
	 << "\n ---> Max values x,y,z: " 
	 << maxx/cm << " " << maxy/cm << " " << maxz/cm << " cm " << endl;
  printf("Field coordinates offset by:\n  x = %f cm, y = %f cm, z = %f cm\n",
	 fXoffset, fYoffset, fZoffset);

  // Should really check that the limits are not the wrong way around.
  if (maxx < minx) {swap(maxx,minx); invertX = true;} 
  if (maxy < miny) {swap(maxy,miny); invertY = true;} 
  if (maxz < minz) {swap(maxz,minz); invertZ = true;} 
  G4cout << "\nAfter reordering if neccesary"  
	 << "\n ---> Min values x,y,z: " 
	 << minx/cm << " " << miny/cm << " " << minz/cm << " cm "
	 << " \n ---> Max values x,y,z: " 
	 << maxx/cm << " " << maxy/cm << " " << maxz/cm << " cm ";

  dx = maxx - minx;
  dy = maxy - miny;
  dz = maxz - minz;
  G4cout << "\n ---> Dif values x,y,z (range): " 
	 << dx/cm << " " << dy/cm << " " << dz/cm << " cm in z "
	 << "\n------------------------------------------------------" << endl;
}

void TOSCAField3D::GetFieldValue(const double point[4],
				      double *Bfield ) const
{

  double x = point[0];
  double y = point[1];
  double z = point[2];

  // Check that the point is within the defined region 
  if ( x>=minx && x<=maxx &&
       y>=miny && y<=maxy && 
       z>=minz && z<=maxz ) {
    
    // Position of given point within region, normalized to the range
    // [0,1]
    double xfraction = (x - minx) / dx;
    double yfraction = (y - miny) / dy; 
    double zfraction = (z - minz) / dz;

    if (invertX) { xfraction = 1 - xfraction;}
    if (invertY) { yfraction = 1 - yfraction;}
    if (invertZ) { zfraction = 1 - zfraction;}

    // Need addresses of these to pass to modf below.
    // modf uses its second argument as an OUTPUT argument.
    double xdindex, ydindex, zdindex;
    
    // Position of the point within the cuboid defined by the
    // nearest surrounding tabulated points
    double xlocal = ( std::modf(xfraction*(nx-1), &xdindex));
    double ylocal = ( std::modf(yfraction*(ny-1), &ydindex));
    double zlocal = ( std::modf(zfraction*(nz-1), &zdindex));
    
    // The indices of the nearest tabulated point whose coordinates
    // are all less than those of the given point
    int xindex = static_cast<int>(xdindex);
    int yindex = static_cast<int>(ydindex);
    int zindex = static_cast<int>(zdindex);
    

#ifdef DEBUG_INTERPOLATING_FIELD
    G4cout << "Local x,y,z: " << xlocal << " " << ylocal << " " << zlocal << endl;
    G4cout << "Index x,y,z: " << xindex << " " << yindex << " " << zindex << endl;
    double valx0z0, mulx0z0, valx1z0, mulx1z0;
    double valx0z1, mulx0z1, valx1z1, mulx1z1;
    valx0z0= table[xindex  ][0][zindex];  mulx0z0=  (1-xlocal) * (1-zlocal);
    valx1z0= table[xindex+1][0][zindex];  mulx1z0=   xlocal    * (1-zlocal);
    valx0z1= table[xindex  ][0][zindex+1]; mulx0z1= (1-xlocal) * zlocal;
    valx1z1= table[xindex+1][0][zindex+1]; mulx1z1=  xlocal    * zlocal;
#endif

        // Full 3-dimensional version
    Bfield[0] =
      xField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      xField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      xField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      xField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      xField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      xField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      xField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      xField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    Bfield[1] =
      yField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      yField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      yField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      yField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      yField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      yField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      yField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      yField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;
    Bfield[2] =
      zField[xindex  ][yindex  ][zindex  ] * (1-xlocal) * (1-ylocal) * (1-zlocal) +
      zField[xindex  ][yindex  ][zindex+1] * (1-xlocal) * (1-ylocal) *    zlocal  +
      zField[xindex  ][yindex+1][zindex  ] * (1-xlocal) *    ylocal  * (1-zlocal) +
      zField[xindex  ][yindex+1][zindex+1] * (1-xlocal) *    ylocal  *    zlocal  +
      zField[xindex+1][yindex  ][zindex  ] *    xlocal  * (1-ylocal) * (1-zlocal) +
      zField[xindex+1][yindex  ][zindex+1] *    xlocal  * (1-ylocal) *    zlocal  +
      zField[xindex+1][yindex+1][zindex  ] *    xlocal  *    ylocal  * (1-zlocal) +
      zField[xindex+1][yindex+1][zindex+1] *    xlocal  *    ylocal  *    zlocal ;

  } else {
    Bfield[0] = 0.0;
    Bfield[1] = 0.0;
    Bfield[2] = 0.0;
  }
}

