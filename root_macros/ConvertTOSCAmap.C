#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include "TString.h"

//Bogdan's TOSCA map is in a coordinate system where -Z is along the optical axis, the origin coincides with the Hall A
// origin, +X is to beam right, and +Y is vertically up:

void ConvertTOSCAmap(const char *mapfilename, const char *outmapfilename){
  ifstream mapfile(mapfilename);

  if( mapfile ){
    ofstream outmapfile(outmapfilename);

    TString currentline;
    while( currentline.ReadLine(mapfile) ){
      if( currentline.BeginsWith("start") ) break;
      cout << currentline.Data() << endl;
      currentline.ReplaceAll('\n',"");
      currentline.ReplaceAll('\r',"");
      outmapfile << currentline.Data() << endl;
    }

    double xtemp, ytemp, ztemp, Bxtemp, Bytemp, Bztemp;
    while ( mapfile >> xtemp >> ytemp >> ztemp >> Bxtemp >> Bytemp >> Bztemp ){
      cout << "(x,y,z, Bx, By, Bz)=(" << xtemp << ", " << ytemp << ", " << ztemp
	   << ", " << Bxtemp << ", " << Bytemp << ", " << Bztemp << endl;
      //so we need to send x, Bx to -x, -Bx, z, Bz to -Bz, -z, and keep y, By the same:
      TString outline;
      outline.Form("%12.6g %12.6g %12.6g %18.12g %18.12g %18.12g", -xtemp, ytemp, -ztemp, Bxtemp, -Bytemp, Bztemp );
      outmapfile << outline.Data() << endl;
      
    }
    
  }
}
