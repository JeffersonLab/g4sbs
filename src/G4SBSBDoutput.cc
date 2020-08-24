#include "G4SBSBDoutput.hh"
//______________________________________________________________________________
G4SBSBDoutput::G4SBSBDoutput(){
   Clear();
}
//______________________________________________________________________________
G4SBSBDoutput::~G4SBSBDoutput(){

}
//______________________________________________________________________________
void G4SBSBDoutput::Clear(){
   nhits = 0;
   plane.clear();
   trid.clear();
   t.clear();
   x.clear();
   y.clear();
   z.clear();
   xg.clear();
   yg.clear();
   zg.clear();
   p.clear();
   pid.clear();
   mid.clear();
   beta.clear();
   edep.clear();
}
