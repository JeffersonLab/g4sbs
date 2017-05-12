// RecTagg a Geant-4 Based Model of RTPC for recoil-proton tagging
// J.R.M Annand, University of Glasgow
// Class ArrayHit
// Define the interesting parameter from the sensitive volumes
// ie the active regions of the detectors
// 06/05/14 JRMA

#include "ArrayHit.hh"
#include "G4Color.hh"
#include "G4VisAttributes.hh"

G4Allocator<ArrayHit> ArrayHitAllocator;

ArrayHit::ArrayHit()
{
  fEdep = 0;
  fEdepB = 0;
  fPos.set(0,0,0);
  fID = -1;
  //fPDG = 0;
  fTime = 0;
  fTimeE = 0;
  fNst = 0;
}

ArrayHit::~ArrayHit()
{ ;}


ArrayHit::ArrayHit(const ArrayHit& right)
  :G4VHit()
{
  fEdep = right.fEdep;
  fEdepB = right.fEdepB;
  fPos = right.fPos;
  fID = right.fID;
  //fPDG = right.fPDG;
  fTime = right.fTime;
  fTimeE = right.fTimeE;
  //fPol = right.fPol;
  fNst = right.fNst;
  for( G4int i=0; i<fNst; i++ )
    fStepList[i] = right.fStepList[i];
}


const ArrayHit& ArrayHit::operator=(const ArrayHit& right)
{
  fEdep = right.fEdep;
  fEdepB = right.fEdepB;
  fPos = right.fPos;
  fID = right.fID;
  //fPDG = right.fPDG;
  fTime = right.fTime;
  fTimeE = right.fTimeE;
  //fPol = right.fPol;
  fNst = right.fNst;
  for( G4int i=0; i<fNst; i++ )
    fStepList[i] = right.fStepList[i];
  return *this;
}

/*
int ArrayHit::operator==(const ArrayHit& right) const
{
  return 0;
}
*/

void ArrayHit::Draw()
{
}


void ArrayHit::Print()
{;}












