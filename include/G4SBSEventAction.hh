
#ifndef G4SBSEventAction_h
#define G4SBSEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#define __MAXGEM 100

class G4Event;
class G4SBSIO;
class G4SBSEventGen;

class G4SBSEventAction : public G4UserEventAction
{
  public:
    G4SBSEventAction();
    virtual ~G4SBSEventAction();

  public:
    virtual void BeginOfEventAction(const G4Event*);
    virtual void EndOfEventAction(const G4Event*);

    void SetIO( G4SBSIO *io ){ fIO = io; }
    void SetEvGen( G4SBSEventGen *g ){ fevgen = g; }
    void SetGEMRes( double r ){ fGEMres = r; }

    void LoadSigmas(const char *filename);


  private:
    G4int gemCollID, hcalCollID, bbcalCollID;

    double fGEMres;

    G4SBSIO *fIO;
    G4SBSEventGen *fevgen;

    double fGEMsigma[__MAXGEM];

  public:
};

#endif

    
