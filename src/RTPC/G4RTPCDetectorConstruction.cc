// ESSN a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class ESSNDetectorMessenger
// Online control of detector configuration via keyboard
// 20/05/13 JRMA adapted from SBS equivalent, under construction

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "DetectorConstruction.hh"
#include "G4SBSDetectorConstruction"
#include "RTPC.hh"


class G4RTPCDetectorConstruction : public DetectorConstruction, public G4SBSDetectorConstruction
{

public:
  G4RTPCDetectorConstruction();
 ~G4RTPCDetectorConstruction(); 
};

#endif

