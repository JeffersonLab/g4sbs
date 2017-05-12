// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class BonusBField
// Read Field map for Bonus solenoid (origin unknown). Phi-symmetric field
// 13/05/14 JRMA

#include "BonusBField.hh"

BonusBField::BonusBField(char* mapname,double xoff, double yoff, double zoff,
			 double scale)
{
  // Read in map. Map coords in cm.
  // Map field values in Gauss
  //
  dR = 0.5;    // step in radius 0.5 cm
  dZ = 0.5;    // step in z 0.5 cm
  Xoff = xoff;
  Yoff = yoff;
  Zoff = zoff;
  Nr = 401;
  Nz = 560;
  BrField = new double[Nr*Nz];
  BzField = new double[Nr*Nz];
  ScaleFactor = scale*gauss;
  printf(" Reading field map from file %s \n", mapname);
  // Open the file for reading.
  FILE* file = fopen(mapname,"r");
  if(!file){
    printf(" Error opening file %s\n Abort Geant-4 session\n",mapname);
    exit(-1);
  }
  char* fl;
  // Read in the data
  double xval,yval,zval,bx,by,bz;
  Rmax = Zmax = 0;
  Rmin = Zmin = 99999;
  char buffer[256];
  for (int iy=0; iy<Nr; iy++) {
    for (int iz=0; iz<Nz; iz++) {
      if(!(fl = fgets(buffer,256,file)) ){
	printf("Premature end of file %s\n",mapname);
	fclose(file);
	exit(-1);
      }
      if( sscanf(buffer,"%lf%lf%lf%lf%lf%lf",&xval,&yval,&zval,&bx,&by,&bz)
	  != 6 ){
	printf("File format error %s, exiting\n",mapname);
	fclose(file);
	exit(-1);
      }
      xval += xoff;
      yval += yoff;
      zval += zoff;
      if(Rmax < yval) Rmax = yval;
      if(Rmin > yval) Rmin = yval;
      if(Zmax < zval) Zmax = zval;
      if(Zmin > zval) Zmin = zval;
      BrField[jk(iy,iz)] = by*ScaleFactor;
      BzField[jk(iy,iz)] = bz*ScaleFactor;
    }
  }
  fclose(file);

  printf("Min r,z: %f %f   Max r,z: %f %f\n", Rmin,Zmin,Rmax,Zmax);
  printf("Field coordinates offset by:\n  x = %f cm, y = %f cm, z = %f cm\n",
	 Xoff, Yoff, Zoff);
}

BonusBField::~BonusBField()
{
  delete BrField;
  delete BzField;
}


//-----------------------------------------------------------------------------
void BonusBField::GetFieldValue(const double pos[4], double* b) const
{
  // 2D linear interpolation of field for given corrdinate pos
  // from field tables. pos from G4 comes in mm so convert to cm
  double r = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
  double cos,sin;
  if( r == 0.0 ){
    cos = sin = 0.0;
  }
  else{
    cos = pos[0]/r;
    sin = pos[1]/r;
  }
  r = r/cm;
  double z = pos[2]/cm;
  double flag = 1.0;
  if( z < 0 ){
    z = -z;
    flag = -1.0;
  }
  if( (r < Rmin) || (r >= Rmax) || (z < Zmin) || (z >= Zmax) ){
    b[0] = b[1] = b[2] = 0.0;
    return;
  }
  int ir = (r - Rmin)/dR;
  int iz = (z - Zmin)/dZ;
  double r1 = ir*dR + Rmin;
  double z1 = iz*dZ + Zmin;
  // printf("%d %d %g %g %g %g %g %g\n",
  //	 ir,iz,r1,z1,r,z,rField[ir][iz],zField[ir][iz]);  
  double br11 = BrField[jk(ir,iz)];
  double br21 = BrField[jk(ir+1,iz)];
  double br12 = BrField[jk(ir,iz+1)];
  double br22 = BrField[jk(ir+1,iz+1)];
  double bra = br11 + (br21 - br11)*(r - r1)/dR;
  double brb = br12 + (br22 - br12)*(r - r1)/dR;
  double brc = bra + (brb - bra)*(z - z1)/dZ;
  //printf("%g %g %g %g %g %g %g\n",br11,br21,br12,br22,bra,brb,brc); 
  //
  double bz11 = BzField[jk(ir,iz)];
  double bz21 = BzField[jk(ir+1,iz)];
  double bz12 = BzField[jk(ir,iz+1)];
  double bz22 = BzField[jk(ir+1,iz+1)];
  double bza = bz11 + (bz21 - bz11)*(r - r1)/dR;
  double bzb = bz12 + (bz22 - bz12)*(r - r1)/dR;
  double bzc = bza + (bzb - bza)*(z - z1)/dZ;
  //printf("%g %g %g %g %g %g %g\n",bz11,bz21,bz12,bz22,bza,bzb,bzc);
  b[0] = brc*cos*flag;
  b[1] = brc*sin*flag;
  b[2] = bzc;
  return;
}

