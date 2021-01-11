#include "G4SBSUtil.hh"
//______________________________________________________________________________
namespace G4SBS {
   namespace Util {
      //______________________________________________________________________________
      void RotateVector(std::vector<G4double> R,G4ThreeVector P,G4ThreeVector &W){
	 // rotate a vector P by angles R, producing new vector W
	 // rotation angles about (x,y,z) = (gamma,beta,alpha) 
	 G4double COS_G = cos(R[0]); G4double COS_B = cos(R[1]); G4double COS_A = cos(R[2]);
	 G4double SIN_G = sin(R[0]); G4double SIN_B = sin(R[1]); G4double SIN_A = sin(R[2]);
	 // compute new coordinates 
	 G4double xp = COS_A*COS_B*P.x() + (COS_A*COS_B*SIN_G - SIN_A*COS_G)*P.y() + (COS_A*SIN_B*COS_G + SIN_A*SIN_G)*P.z();
	 G4double yp = SIN_A*COS_B*P.x() + (SIN_A*SIN_B*SIN_G + COS_A*COS_G)*P.y() + (SIN_A*SIN_B*COS_G - COS_A*SIN_G)*P.z();
	 G4double zp =      -SIN_B*P.x() +                       COS_B*SIN_G*P.y() +                       COS_B*COS_G*P.z();
	 // set new values 
	 W.setX(xp);
	 W.setY(yp);
	 W.setZ(zp);
      }
   } //::Util
} //::G4SBS
