#include "G4SBSICoutput.hh"
//______________________________________________________________________________
G4SBSICoutput::G4SBSICoutput(){
   Clear();
}
//______________________________________________________________________________
G4SBSICoutput::~G4SBSICoutput(){

}
//______________________________________________________________________________
void G4SBSICoutput::Clear(){
   nhits_IC = 0;
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
   ParticleHistory.Clear();
}
