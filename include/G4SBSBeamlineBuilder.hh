#ifndef __G4SBSBeamlineBuilder_hh
#define __G4SBSBeamlineBuilder_hh

#include "G4SBSComponent.hh"

class G4LogicalVolume;
class GSBSDetectorConstruction;

class G4SBSBeamlineBuilder: public G4SBSComponent {
public:
  G4SBSBeamlineBuilder(G4SBSDetectorConstruction *);
  ~G4SBSBeamlineBuilder();

  void BuildComponent(G4LogicalVolume *);

private:

  void MakeGEpLead(G4LogicalVolume *);
  void MakeGMnLead(G4LogicalVolume *);
  void MakeGEnRPLead(G4LogicalVolume *); // extra shielding after magnet yoke
  void MakeGEnLead(G4LogicalVolume *);
  void MakeGEnClamp(G4LogicalVolume *);
  void MakeSIDISLead( G4LogicalVolume * );

  //void MakeEntranceBeamline(G4LogicalVolume *);
  void MakeCommonExitBeamline(G4LogicalVolume *);
  void MakeGEpBeamline(G4LogicalVolume *);
  void MakeGMnBeamline(G4LogicalVolume *);
  void Make3HeBeamline(G4LogicalVolume *);// for GEn, A1n, SIDIS
  void MakeDefaultBeamline(G4LogicalVolume *);// Old beam line...
  void MakeToyBeamline(G4LogicalVolume *); // "Toy" beam line for playing around with extreme forward angles of detectors, etc:

  // for a beam diffuser that sits in front of the dump (useful for beam steering studies)
  G4double fBDLength;
  void MakeBeamDump(G4LogicalVolume *logicMother);  // default location  
  void MakeBeamDump_Diffuser(G4LogicalVolume *logicMother,G4double z0=0);
  void MakeBeamDump_ISOWallWeldment(G4LogicalVolume *logicMother,G4double z0=0);
  void MakeBeamDump_UpstreamPipe(G4LogicalVolume *logicMother,G4double z0=0);
  void MakeBeamDump_DownstreamPipe(G4LogicalVolume *logicMother,G4double z0=0);
  
};

#endif//__G4SBSBeamlineBuilder_hh
