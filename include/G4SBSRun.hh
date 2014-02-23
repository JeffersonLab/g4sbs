#ifndef __G4SBSRUN_HH
#define __G4SBSRUN_HH

/*!
 * All the information on the run
 * The data object will get put into the output
 * stream
  
   This is implemented in the soliton model
 */

#include "G4SBSRunData.hh"

class G4SBSRun {

    private:
	static G4SBSRun *gSingleton;
	 G4SBSRun();

	G4SBSRunData *fRunData;

    public:
	 static G4SBSRun *GetRun();
	~G4SBSRun();

	G4SBSRunData *GetData(){return fRunData;}
};

#endif//__G4SBSRUN_HH
