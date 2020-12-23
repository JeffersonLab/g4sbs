#include "G4SBSTargetoutput.hh"
//______________________________________________________________________________
G4SBSTargetoutput::G4SBSTargetoutput(){
   Clear();
}
//______________________________________________________________________________
G4SBSTargetoutput::~G4SBSTargetoutput(){

}
//______________________________________________________________________________
void G4SBSTargetoutput::Clear(){
   nhits_Target = 0;
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
   trackLength.clear();
   ParticleHistory.Clear();
}
