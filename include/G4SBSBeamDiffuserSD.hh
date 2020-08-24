// BeamDiffuser (BD) Sensitive Detector (SD) class 
// - Defines what happens event-by-event for the beam diffuser using the BDHit class
// - Author: D. Flay (JLab) 

#ifndef G4SBS_BEAM_DIFFUSER_SD_HH
#define G4SBS_BEAM_DIFFUSER_SD_HH

#include "G4VSensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

#include "G4SBSBDHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

class G4SBSBeamDiffuserSD : public G4VSensitiveDetector
{
   public:
      G4SBSBeamDiffuserSD(const G4String& name,const G4String& hitsCollectionName);
      virtual ~G4SBSBeamDiffuserSD();

      // methods from base class
      virtual void   Initialize(G4HCofThisEvent* hitCollection);
      virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
      virtual void   EndOfEvent(G4HCofThisEvent* hitCollection);

   private:
      G4SBSBDHitsCollection* fHitsCollection;
};

#endif
