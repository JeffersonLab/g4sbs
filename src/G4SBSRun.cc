#include "G4SBSRun.hh"
#include "G4SBSRunData.hh"

G4SBSRun *G4SBSRun::gSingleton = NULL;

G4SBSRun::G4SBSRun(){
    gSingleton = this;
    fRunData = new G4SBSRunData();
    fRunData->Init();
}

G4SBSRun::~G4SBSRun(){
}

G4SBSRun *G4SBSRun::GetRun(){
    if( gSingleton == NULL ){
	gSingleton = new G4SBSRun();
    }
    return gSingleton;
}
