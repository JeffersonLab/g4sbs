#ifndef __Beamline_hh
#define __Beamline_hh

class G4LogicalVolume;
class DetectorConstruction;

class Beamline {
public:
  Beamline(DetectorConstruction*);
  ~Beamline();
  void BuildComponent(G4LogicalVolume *);
private:
  //void MakeGEpLead(G4LogicalVolume *);
  //void MakeGEnLead(G4LogicalVolume *);
  //void MakeGEnClamp(G4LogicalVolume *);
  //void MakeSIDISLead( G4LogicalVolume * );
  bool fOvLap = false;
  DetectorConstruction* fRtag;
};
#endif
