#ifndef G4SBSECalHit_h 
#define G4SBSECalHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"

class G4SBSECalHit : public G4VHit
{
public:
  G4SBSECalHit();
  ~G4SBSECalHit();
  
  G4SBSECalHit( const G4SBSECalHit &hit );
  const G4SBSECalHit & operator=( const G4SBSECalHit &hit );
  
  G4int operator==(const G4SBSECalHit &hit ) const;

  inline void * operator new(size_t);
  inline void operator delete(void *aHit);
  
  virtual void Draw();
  virtual void Print();

  //Utility functions to map rows and columns to PMT numbers
  // G4int calc_row( G4int PMT );
  // G4int calc_col( G4int PMT );

private:
  
  //What are the relevant properties of a ECal hit?
  //    1) Timing
  //    2) # of photoelectrons
  //    3) Energy deposited by primary particle
  //    4) What else? Computationally hard to keep track of certain 
  //       quantities because we expect ~10^4 - 10^5 optical photons to be 
  //       created during an EM shower.

  G4int fTrackID;              //track number of current track: needed to access optical photon energy 
  G4ThreeVector fpos;          //position of the step (global) 
  G4ThreeVector fLpos;         //position of the step  (local to the volume)
  G4double fedep;              //Energy deposit
  G4double fdx;                //Step length
  G4double fenergy;            //optical photon energy in this step
  G4double ftime;              //Global time since start of event
  G4int fPMTnumber;            //PMT number;
  G4int frownumber;            //PMT row number
  G4int fcolnumber;            //PMT column number
  G4int fplanenumber;          //PMT "plane" number:
  G4ThreeVector CellCoords;    //"local" coordinate of center of cell in which hit occurs.
  G4ThreeVector GlobalCellCoords; //"global" coordinate of center of cell in which hit occurs:

  //Better than carrying around all the baggage associated with a logical volume pointer; let's just grab the
  //detection efficiency if it has been defined for the material in which this hit occurs:
  G4double fQuantumEfficiency;

  //Indices for track info:
  G4int otridx, ptridx, sdtridx;
  
public:
  //Get/set methods for sensitive detector:
  
  //Track IDs:
  inline void SetTrackID( G4int tid ){ fTrackID = tid; }
  inline G4int GetTrackID() const { return fTrackID; }

  //Track positions and vertices:
  inline void SetPos( G4ThreeVector x ){ fpos = x; }
  inline G4ThreeVector GetPos() const { return fpos; }

  inline void SetLPos( G4ThreeVector x ){ fLpos = x; }
  inline G4ThreeVector GetLPos() const { return fLpos; }

  //Step physical properties (energy deposit, energy, time, step length):
  inline void SetEdep( G4double dE ){ fedep = dE; }
  inline G4double GetEdep() const { return fedep; }
  
  inline void Setdx( G4double d ){ fdx = d; }
  inline G4double Getdx() const { return fdx; }
  
  inline void Setenergy( G4double E ){ fenergy = E; }
  inline G4double Getenergy() const { return fenergy; }
  
  inline void SetTime( G4double t ){ ftime = t; }
  inline G4double GetTime() const { return ftime; }

  //Digital information regarding PMT (PMT number, row number, column number):
  inline void SetPMTnumber( G4int n ){ fPMTnumber = n; }
  inline G4int GetPMTnumber() const { return fPMTnumber; }

  inline void Setrownumber( G4int i ){ frownumber = i; }
  inline G4int Getrownumber() const { return frownumber; }

  inline void Setcolnumber( G4int i ){ fcolnumber = i; }
  inline G4int Getcolnumber() const { return fcolnumber; }

  inline void Setplanenumber( G4int i ){ fplanenumber = i; }
  inline G4int Getplanenumber() const { return fplanenumber; }
  
  inline void SetCellCoords( G4ThreeVector x ){CellCoords = x;}
  inline G4ThreeVector GetCellCoords() const { return CellCoords; }
  
  inline void SetGlobalCellCoords( G4ThreeVector x ){ GlobalCellCoords = x; }
  inline G4ThreeVector GetGlobalCellCoords() const { return GlobalCellCoords; }
  
  
  // inline void SetOriginVol( G4int i ){ foriginvol = i; }
  // inline G4int GetOriginVol() const { return foriginvol; }

  // inline void SetLogicalVolume( G4LogicalVolume *v ){  fvolume_log = v; }
  //inline G4LogicalVolume *GetLogicalVolume() const { return fvolume_log; }
  inline void SetQuantumEfficiency( G4double QE ){ fQuantumEfficiency = QE; }
  inline G4double GetQuantumEfficiency() const { return fQuantumEfficiency; }

  inline void SetOTrIdx(G4int idx){ otridx = idx; }
  inline void SetPTrIdx(G4int idx){ ptridx = idx; }
  inline void SetSDTrIdx(G4int idx){ sdtridx = idx; }

  inline G4int GetOTrIdx() const { return otridx; }
  inline G4int GetPTrIdx() const { return ptridx; }
  inline G4int GetSDTrIdx() const { return sdtridx; }
};

typedef G4THitsCollection<G4SBSECalHit> G4SBSECalHitsCollection;

extern G4Allocator<G4SBSECalHit> *G4SBSECalHitAllocator;

inline void *G4SBSECalHit::operator new(size_t)
{
  if( !G4SBSECalHitAllocator )  G4SBSECalHitAllocator = new G4Allocator<G4SBSECalHit>;
  return (void *) G4SBSECalHitAllocator->MallocSingle();
}

inline void G4SBSECalHit::operator delete(void *aHit)
{
  G4SBSECalHitAllocator->FreeSingle( (G4SBSECalHit*) aHit );
}

#endif
