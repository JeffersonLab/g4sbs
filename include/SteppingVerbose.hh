// RTPC a Geant-4 Based Model of the Lund ESS Neutron Test Facility
// J.R.M Annand, University of Glasgow
// Class SteppingVerbose
// Output details of stepping action
// 20/05/13 JRMA adapted from SBS equivalent, under construction
//
class SteppingVerbose;

#ifndef SteppingVerbose_h
#define SteppingVerbose_h 1

#include "G4SteppingVerbose.hh"

class SteppingVerbose : public G4SteppingVerbose
{
 public:   
   SteppingVerbose();
  ~SteppingVerbose();
   void StepInfo();
   void TrackingStarted();
};

#endif
