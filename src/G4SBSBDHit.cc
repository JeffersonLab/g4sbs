#include "G4SBSBDHit.hh"
G4ThreadLocal G4Allocator<G4SBSBDHit>* G4SBSBDHitAllocator = 0;
//______________________________________________________________________________
G4SBSBDHit::G4SBSBDHit()
 : G4VHit(),
   fEdep(0.),
   fTrackLength(0.),
   fEtot(0.),
   fBeta(0.),
   fHitTime(0.),
   fTrackID(-1),
   fPlane(-1),
   fPID(-1),
   fMID(-1),
   fVerbosity(0)
{
   fPos.setX(0);    fPos.setY(0);    fPos.setZ(0);
   fLabPos.setX(0); fLabPos.setY(0); fLabPos.setZ(0);
   fMom.setX(0);    fMom.setY(0);    fMom.setZ(0);
}
//______________________________________________________________________________
G4SBSBDHit::~G4SBSBDHit() 
{

}
//______________________________________________________________________________
G4SBSBDHit::G4SBSBDHit(const G4SBSBDHit& right)
  : G4VHit()
{
  fEdep        = right.fEdep;
  fTrackLength = right.fTrackLength;
  fEtot        = right.fEtot;
  fBeta        = right.fBeta;
  fHitTime     = right.fHitTime;
  fPlane       = right.fPlane;
  fPID         = right.fPID;
  fMID         = right.fMID;
  fVerbosity   = right.fVerbosity; 
  fPos         = right.fPos;
  fLabPos      = right.fLabPos;
  fMom         = right.fMom;
}
//______________________________________________________________________________
const G4SBSBDHit& G4SBSBDHit::operator=(const G4SBSBDHit& right)
{
  fEdep        = right.fEdep;
  fTrackLength = right.fTrackLength;
  fEtot        = right.fEtot;
  fBeta        = right.fBeta;
  fHitTime     = right.fHitTime;
  fPlane       = right.fPlane;
  fPID         = right.fPID;
  fMID         = right.fMID;
  fVerbosity   = right.fVerbosity; 
  fPos         = right.fPos;
  fLabPos      = right.fLabPos;
  fMom         = right.fMom;
  return *this;
}
//______________________________________________________________________________
G4bool G4SBSBDHit::operator==(const G4SBSBDHit& right) const
{
  return ( this == &right ) ? true : false;
}
//______________________________________________________________________________
void G4SBSBDHit::Print()
{
   if(fVerbosity>0){
      G4cout
	 << "Edep: " 
	 << std::setw(7) << G4BestUnit(fEdep,"Energy")
	 << " track length: " 
	 << std::setw(7) << G4BestUnit(fTrackLength,"Length")
	 << " track x: " 
	 << std::setw(7) << G4BestUnit(fPos.getX(),"Length") 
	 << " track y: " 
	 << std::setw(7) << G4BestUnit(fPos.getY(),"Length") 
	 << " track z: " 
	 << std::setw(7) << G4BestUnit(fPos.getZ(),"Length") 
	 << " plane no: " 
	 << std::setw(7) << fPlane  
	 << G4endl;
   }
}
