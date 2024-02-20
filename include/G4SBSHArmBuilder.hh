#ifndef __G4SBSHArmBuilder_hh
#define __G4SBSHArmBuilder_hh

#include "G4SBSComponent.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class G4SBSHArmBuilder: public G4SBSComponent {
public:
  G4SBSHArmBuilder(G4SBSDetectorConstruction *);
  ~G4SBSHArmBuilder();

  void BuildComponent(G4LogicalVolume *);

  void Set48D48Ang(double a){ f48D48ang = a; }
  void SetHCALDist(double a){ fHCALdist= a;   }
  void SetHCALVOffset(double a){ fHCALvertical_offset = a; }
  void SetHCALHOffset(double a){ fHCALhorizontal_offset = a; }
  void SetHCALAngOffset(double a){ fHCALangular_offset = a; }
  void Set48D48Dist(double a){ f48D48dist = a; }
  void SetLACDist( double a ){ fLACdist = a; }
  void SetLACVOffset( double a ){ fLACvertical_offset = a; }
  void SetLACHOffset( double a ){ fLAChorizontal_offset = a; }
  void SetRICHdist( double d ){ fRICHdist = d; } //Set RICH detector distance
  void SetRICHHoffset( double d ){ fRICHhorizontal_offset = d; }
  void SetRICHVoffset( double d ){ fRICHvertical_offset = d; }
  void SetRICHSnoutExtension( double d ){ fRichSnoutExtension = d; }
  void SetFieldClampConfig48D48( int option ){ f48D48_fieldclamp_config = option; }
  void SetTrackerPitch(double a){ fSBS_tracker_pitch = a; }
  void SetTrackerDist(double d){ fSBS_tracker_dist = d; }
  void SetSBSSieve(bool a){fBuildSBSSieve = a;};
  void SetRICHgas( G4String s ){ fRICHgas = s; }
  void SetRICH_use_aerogel( G4bool b ){ fRICH_use_aerogel = b; }
  void SetFPP_CH2thick(int ifpp, double thick );
  void SetGEPFPPoption(int option){ fGEPFPPoption = option; }
  void SetGENRPAnalyzerOption(int ia ){ fGEnRP_analyzer_option = ia; }
  void SetFTabsthick( G4double t ){ fFTabsthick = t; }
  void SetFTuseabsorber( G4bool b ){ fFTuseabsorber = b; }
  void SetFTabsmaterial( G4String s ){ fFTabsmaterial = s; }
  void SetCDetReady( G4bool b ){fCDetReady = b;}
  
  void Make48D48(G4LogicalVolume*, double);
  void MakeSBSFieldClamps(G4LogicalVolume*);
  void MakeHCAL(G4LogicalVolume*, G4double);
  void MakeHCALV2(G4LogicalVolume*, G4double);
  void MakeCDET(G4LogicalVolume*, G4double, G4double);
  void MakeNeutronVeto(G4LogicalVolume*, G4double);
  void MakeFPP(G4LogicalVolume*, G4RotationMatrix*, G4ThreeVector );
  void MakeRICH(G4LogicalVolume *);
  void MakeRICH_new(G4LogicalVolume *, bool extended_snout = false);
  void MakeTracker(G4LogicalVolume *, G4int nplanes = 5);
  void MakeGEpFPP(G4LogicalVolume *);
  //void MakeGEpFPP_BB( G4LogicalVolume *);
  void MakeTracker_A1n(G4LogicalVolume *);
  void MakeElectronModeSBS(G4LogicalVolume *);
  void MakeSBSSieveSlit(G4LogicalVolume *);
  void MakeLAC(G4LogicalVolume *);
  void MakePolarimeterGEnRP(G4LogicalVolume *);

  //************************
  void MakeTest(G4LogicalVolume *);
  
  double f48D48ang;
  double f48D48dist;
  int f48D48_fieldclamp_config; //Configuration of field clamp. There could be several of these.
  double fHCALdist;
  double fHCALvertical_offset;  // Vertical offset (from center) of HCAL
  double fHCALhorizontal_offset; // Horizontal offset (from SBS center line) of HCAL (by convention, +X is toward smaller angle)
  double fHCALangular_offset; // Angular offset of HCAL wrt exit beamline (+ = away from beamline).
  double fRICHdist; //Distance from target to RICH entry window (in horizontal plane)
  double fRICHvertical_offset; //Vertical offset (from center)
  double fRICHhorizontal_offset; //Horizontal offset (from SBS center line, + = toward beamline).
  double fSBS_tracker_pitch; //SBS tracker pitch angle for electron mode
  double fSBS_tracker_dist; 
  
  // Useful constants
  double f48D48depth;
  double f48D48height;
  double f48D48width;

  bool fUseLocalField;
  
  double fFieldStrength;
  
  bool fBuildSBSSieve;

  bool fCDetReady;
  bool fUseNeutronVeto;
  
  G4String fRICHgas; //String defining
  G4bool   fRICH_use_aerogel; // Flag to use or not use aerogel
  G4double fRichSnoutExtension;
  
  int    fGEnRP_analyzer_option; // configuration of GEn-RP analyzer

  double fLACdist; //Distance to CLAS Large-angle calorimeter
  double fLACvertical_offset; //vertical offset of center of LAC wr
  double fLAChorizontal_offset; //horizontal offset of center of LAC

  double fCH2thickFPP[2];

  //Parameters of possible FT absorber for GEP to reduce backgrounds:
  G4double fFTabsthick;
  G4bool   fFTuseabsorber; //flag to toggle building of an absorber
  G4String fFTabsmaterial; //string defining FT absorber material

  int fGEPFPPoption; //options for different arrangements of the GEP FPP:

  //Make the width and height of the CH2 blocks user-configurable:
  double fGEP_CH2width[2];
  double fGEP_CH2height[2]; 
  double fGEP_CH2zpos[2];
  double fGEP_CH2yoff[2];
  //double fGEP_trkryoff[3];
  
private:

};

#endif//__G4SBSHArmBuilder_hh
