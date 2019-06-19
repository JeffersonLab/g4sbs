//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file runAndEvent/G4SBS/include/G4SBSTrackingAction.hh
/// \brief Definition of the G4SBSTrackingAction class
//
//
//Use example RE01 as template, modify for G4SBS as needed

#ifndef G4SBSTrackingAction_h
#define G4SBSTrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "G4SBSDetectorConstruction.hh"

class G4SBSTrackingAction : public G4UserTrackingAction 
{
public:
  G4SBSTrackingAction();
  virtual ~G4SBSTrackingAction(){};
   
  virtual void PreUserTrackingAction(const G4Track*);
  virtual void PostUserTrackingAction(const G4Track*);

  //should only be invoked at start of run, after the geometry is built!
  inline void SetTargetVolumes( set<G4String> tvlist ){ fTargetVolumes = tvlist; } 
  inline void SetAnalyzerVolumes( set<G4String> avlist ){ fAnalyzerVolumes = avlist; }
  
  void Initialize( G4SBSDetectorConstruction *fdc );
  
private:

  //G4SBSDetectorConstruction *fdetcon;
  set<G4String> fTargetVolumes; //list of logical volume names to be flagged as "TARGET"
  set<G4String> fAnalyzerVolumes; //list of logical volume names to be flagged as "ANALYZER"
  
  
};

#endif
