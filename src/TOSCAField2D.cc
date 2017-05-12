// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class TOSCAField2D
// TOSCA generated field map. Phi-symmetric field
// 13/05/14 JRMA

#include "TOSCAField2D.hh"

TOSCAField2D::TOSCAField2D( const char* filename, double xoff, double yoff,
			    double zoff, double scale) 
{    
  double lenUnit   = cm;
  double fieldUnit = gauss; 
  //fXoffset = xoff/lenUnit;
  //fYoffset = yoff/lenUnit;
  //fZoffset = zoff/lenUnit;
  fXoffset = xoff;
  fYoffset = yoff;
  fZoffset = zoff;
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
    if(sscanf(buffer,"%u%u%u",&nx,&ny,&nz) == 3) break;
  }
  printf(" X-dimension: %d,  Y-dimension: %d,  Z-dimension: %d\n",nx,ny,nz);
  // Set up storage space for table
  BrField.resize( nx );
  BzField.resize( nx );
  unsigned int ix,iy,iz;
  for (ix=0; ix<nx; ix++) {
    BrField[ix].resize(nz);
    BzField[ix].resize(nz);
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
	//xval += fXoffset;
	//yval += fYoffset;
	//zval += fZoffset;
        if ( ix==0 && iy==0 && iz==0 ) {
          Rmin = xval * lenUnit;
          Zmin = zval * lenUnit;
        }
	if( yval == 1 ) continue;
        BrField[ix][iz] = bx * fieldUnit;
        BzField[ix][iz] = bz * fieldUnit;
      }
    }
  }
  fclose(file);

  Rmax = xval * lenUnit;
  Zmax = zval * lenUnit;
  printf("Min r,z: %f %f (mm)   Max r,z: %f %f (mm)\n", Rmin,Zmin,Rmax,Zmax);
  printf("Field coordinates offset by:\n  x = %f mm, y = %f mm, z = %f mm\n",
	 fXoffset, fYoffset, fZoffset);

  dR = (Rmax - Rmin)/(nx-1);
  dZ = (Zmax - Zmin)/(nz-1);
  printf("Map step values dR = %g mm    dZ = %g mm\n",dR,dZ);
}

//-----------------------------------------------------------------------------
void TOSCAField2D::GetFieldValue(const double pos[4], double* b) const
{
  // 2D linear interpolation of field for given corrdinate pos
  // from field tables. pos assumed to be in mm
  //
  double x = pos[0] + fXoffset;
  double y = pos[1] + fYoffset;
  double z = pos[2] + fZoffset;
  double r = sqrt(x*x + y*y);
  double cs,sn;
  if( r == 0.0 ){
    cs = sn = 0.0;
  }
  else{
    cs = x/r;
    sn = y/r;
  }
  //r = r;
  //double z = pos[2];
  double flag = 1.0;
  if( z < 0 ){
    z = -z;
    flag = -1.0;
  }
  if( (r < Rmin) || (r >= Rmax) || (z < Zmin) || (z >= Zmax) ){
    b[0] = b[1] = b[2] = 0.0;
    return;
  }
  unsigned int ir = (r - Rmin)/dR;
  unsigned int iz = (z - Zmin)/dZ;
  double r1 = ir*dR + Rmin;
  double z1 = iz*dZ + Zmin;
  // if((ir > nx-2) || (iz > nz-2)){
  //   b[0] = b[1] = b[2] = 0.0;
  //    return;
  // }
  // printf("%d %d %g %g %g %g %g %g\n",
  //	 ir,iz,r1,z1,r,z,rField[ir][iz],BzField[ir][iz]);  
  double br11 = BrField[ir][iz];
  double br21 = BrField[ir+1][iz];
  double br12 = BrField[ir][iz+1];
  double br22 = BrField[ir+1][iz+1];
  double bra = br11 + (br21 - br11)*(r - r1)/dR;
  double brb = br12 + (br22 - br12)*(r - r1)/dR;
  double brc = bra + (brb - bra)*(z - z1)/dZ;
  //printf("%g %g %g %g %g %g %g\n",br11,br21,br12,br22,bra,brb,brc); 
  //
  double bz11 = BzField[ir][iz];
  double bz21 = BzField[ir+1][iz];
  double bz12 = BzField[ir][iz+1];
  double bz22 = BzField[ir+1][iz+1];
  double bza = bz11 + (bz21 - bz11)*(r - r1)/dR;
  double bzb = bz12 + (bz22 - bz12)*(r - r1)/dR;
  double bzc = bza + (bzb - bza)*(z - z1)/dZ;
  //printf("%g %g %g %g %g %g %g\n",bz11,bz21,bz12,bz22,bza,bzb,bzc);
  b[0] = brc*cs*flag;
  b[1] = brc*sn*flag;
  b[2] = bzc;
  return;
}


