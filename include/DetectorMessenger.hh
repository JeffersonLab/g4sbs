// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class ESSNDetectorMessenger
// Online control of detector configuration via keyboard
// 20/05/13 JRMA adapted from SBS equivalent, under construction

#ifndef DetectorMessenger_h
#define DetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithADouble;

class DetectorMessenger: public G4UImessenger{
public:
  DetectorMessenger(DetectorConstruction* );
  ~DetectorMessenger();
  void SetNewValue(G4UIcommand*, G4String);
private:
  DetectorConstruction* fDetector;
  G4UIdirectory*           fSBSDir;
  G4UIdirectory*           fDataDir;
  G4UIcmdWithAnInteger*    fPlugAppCmd;    
  G4UIcmdWithAString*      fTargetMatCmd;
  G4UIcmdWithADouble*      fPlugLengthCmd;
  G4UIcmdWithADouble*      fTargetRadiusCmd;
  G4UIcmdWithoutParameter* fUseSrcPbCmd;
  G4UIcmdWithoutParameter* fUpdateCmd;
  G4UIcmdWithoutParameter* fOverlapVolCmd;
 };

#endif

