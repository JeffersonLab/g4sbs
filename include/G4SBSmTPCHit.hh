#ifndef G4SBSmTPCHit_h 
#define G4SBSmTPCHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"

class G4SBSmTPCHit : public G4VHit
{
public:
  G4SBSmTPCHit();
  ~G4SBSmTPCHit();
  
  G4SBSmTPCHit( const G4SBSmTPCHit &hit );
  const G4SBSmTPCHit & operator=( const G4SBSmTPCHit &hit );
  
  G4int operator==(const G4SBSmTPCHit &hit ) const;

  inline void * operator new(size_t);
  inline void operator delete(void *aHit);
  
  void Draw();
  void Print();
  // std::vector<G4double> integrateDgt(G4ThreeVector lpos, G4double edep, G4ThreeVector lmom);


  //const std::map<G4String, G4AttDef> *GetAttDefs() const;
  //std::vector<G4AttValue> *CreateAttValues() const;

  //Utility functions to map rows and columns to PMT numbers
  // G4int calc_row( G4int PMT );
  // G4int calc_col( G4int PMT );

private:

  G4int fTrackID;              //track number of current track
  G4int fMotherID;             //track number of mother of current track
  G4int fTrackPID;             //particle ID of current track
  // G4int fMotherPID;            //particle ID of current track
  G4ThreeVector fpos;          //position of the step (global)
  G4ThreeVector fLpos;         //position of the step  (local to the volume)
  G4ThreeVector fdirection;    //momentum direction of the step
  G4ThreeVector fLdirection;   //Local direction of the step:
  G4ThreeVector fvertex;       //Primary vertex of the current track
  // G4ThreeVector fMothervertex; //Primary vertex of the track's mother
  // G4ThreeVector fvertexdirection; //momentum direction at the primary vertex;
  // G4ThreeVector fMothervertexdirection; //Mother particle momentum direction
  G4double fedep;              //Energy deposit
  G4double fdx;                //Step length
  G4double fenergy;            //optical photon energy in this step
  G4double ftime;              //Global time since start of event
  G4ThreeVector fmom;              //Global momentum of particle
  G4ThreeVector fLmom;              //Local momentum of particle
  G4double fCell;        // mtpc chamber no
  // G4int fPMTnumber;            //PMT number;
  // G4int frownumber;            //PMT row number
  // G4int fcolnumber;            //PMT column number
  // G4int foriginvol;            //Volume in which the track of this optical photon originated (1=Aerogel, 2=RICHbox, 3=lucite, 0=other)
  G4ThreeVector fCellCoord;       //Local coordinates of center of PMT photocathode
  G4ThreeVector fGlobalCellCoord; //Global coordinates of center of PMT photocathode
  //Need to store a pointer to logical volume in order to retrieve quantum efficiency info later:
  G4LogicalVolume *fvolume_log;

  G4double fZTravel;            // z travel
  G4int fNStrips;               // estimated number of readout pads hit
  
  //rewrite for more compact tree structure: store G4 Track IDs in the hit classes and then
  //use the existing maps in the SDtrackoutput classes:
  G4int otridx, ptridx, sdtridx; 
  //static std::map<G4String, G4AttDef> fAttDefs; //For visualization
  
public:
  //Get/set methods for sensitive detector:
  
  //Track IDs:
  inline void SetTrackID( G4int tid ){ fTrackID = tid; }
  inline G4int GetTrackID() const { return fTrackID; }
  
  inline void SetMotherID( G4int mid ){ fMotherID = mid; }
  inline G4int GetMotherID() const { return fMotherID; }
  
  //Particle IDs:
  inline void SetTrackPID( G4int pid ){ fTrackPID = pid; }
  inline G4int GetTrackPID() const { return fTrackPID; }
  
  // inline void SetMotherPID( G4int mid ){ fMotherPID = mid; }
  // inline G4int GetMotherPID() const { return fMotherPID; }

  //Track positions and vertices:
  inline void SetPos( G4ThreeVector x ){ fpos = x; }
  inline G4ThreeVector GetPos() const { return fpos; }

  inline void SetLPos( G4ThreeVector x ){ fLpos = x; }
  inline G4ThreeVector GetLPos() const { return fLpos; }
  
  inline void SetDirection( G4ThreeVector n ){ fdirection = n.unit(); }
  inline G4ThreeVector GetDirection() const { return fdirection; }
  
  inline void SetLDirection( G4ThreeVector n ){ fLdirection = n.unit(); }
  inline G4ThreeVector GetLDirection() const { return fLdirection; }

  inline void SetVertex( G4ThreeVector v ){ fvertex = v; }
  inline G4ThreeVector GetVertex() const { return fvertex; }

  // inline void SetVertexDirection( G4ThreeVector n ) { fvertexdirection = n.unit(); }
  // inline G4ThreeVector GetVertexDirection() const { return fvertexdirection; }

  // inline void SetMotherVertex( G4ThreeVector v ){ fMothervertex = v; }
  // inline G4ThreeVector GetMotherVertex() const { return fMothervertex; }

  // inline void SetMotherVertexDirection( G4ThreeVector n ) { fMothervertexdirection = n.unit(); }
  // inline G4ThreeVector GetMotherVertexDirection() const { return fMothervertexdirection; }

  //Step physical properties (energy deposit, energy, time, step length, mom):
  inline void SetEdep( G4double dE ){ fedep = dE; }
  inline G4double GetEdep() const { return fedep; }
  
  inline void Setdx( G4double d ){ fdx = d; }
  inline G4double Getdx() const { return fdx; }
  
  inline void SetEnergy( G4double E ){ fenergy = E; }
  inline G4double GetEnergy() const { return fenergy; }
  
  inline void SetTime( G4double t ){ ftime = t; }
  inline G4double GetTime() const { return ftime; }

  inline void SetMom( G4ThreeVector m ){ fmom = m; }
  inline G4ThreeVector GetMom() const { return fmom; }

  inline void SetLMom( G4ThreeVector lm ){ fLmom = lm; }
  inline G4ThreeVector GetLMom() const { return fLmom; }

  //Digital information regarding PMT (PMT number, row number, column number):
  inline void SetCell( G4int n ){ fCell = n; }
  inline G4int GetCell() const { return fCell; }

  // inline void SetPMTnumber( G4int n ){ fPMTnumber = n; }
  // inline G4int GetPMTnumber() const { return fPMTnumber; }

  // inline void Setrownumber( G4int i ){ frownumber = i; }
  // inline G4int Getrownumber() const { return frownumber; }

  // inline void Setcolnumber( G4int i ){ fcolnumber = i; }
  // inline G4int Getcolnumber() const { return fcolnumber; }
  
  // inline void SetOriginVol( G4int i ){ foriginvol = i; }
  // inline G4int GetOriginVol() const { return foriginvol; }

  inline void SetCellCoord( G4ThreeVector x ){ fCellCoord = x; }
  inline G4ThreeVector GetCellCoord() const { return fCellCoord; }
  
  inline void SetGlobalCellCoord( G4ThreeVector x ){ fGlobalCellCoord = x; }
  inline G4ThreeVector GetGlobalCellCoord() const { return fGlobalCellCoord; }

  inline void SetLogicalVolume( G4LogicalVolume *v ){  fvolume_log = v; }
  inline G4LogicalVolume *GetLogicalVolume() const { return fvolume_log; }

  inline void SetZTravel( G4double z ){ fZTravel = z; }
  inline G4double GetZTravel() const { return fZTravel; }

  inline void SetNStrips( G4double nstrips ){ fNStrips = nstrips; }
  inline G4int GetNStrips() const { return fNStrips; }

  inline void SetOTrIdx(G4int idx){ otridx = idx; }
  inline void SetPTrIdx(G4int idx){ ptridx = idx; }
  inline void SetSDTrIdx(G4int idx){ sdtridx = idx; }
  
  inline G4int GetOTrIdx() const { return otridx; }
  inline G4int GetPTrIdx() const { return ptridx; }
  inline G4int GetSDTrIdx() const { return sdtridx; }
};

typedef G4THitsCollection<G4SBSmTPCHit> G4SBSmTPCHitsCollection;

extern G4Allocator<G4SBSmTPCHit> *G4SBSmTPCHitAllocator;

inline void *G4SBSmTPCHit::operator new(size_t)
{
  if( !G4SBSmTPCHitAllocator )  G4SBSmTPCHitAllocator = new G4Allocator<G4SBSmTPCHit>;
  return (void *) G4SBSmTPCHitAllocator->MallocSingle();
}

inline void G4SBSmTPCHit::operator delete(void *aHit)
{
  G4SBSmTPCHitAllocator->FreeSingle( (G4SBSmTPCHit*) aHit );
}
				  


#endif
