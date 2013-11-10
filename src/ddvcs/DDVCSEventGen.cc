#include "DDVCSEventGen.hh"

DDVCSEventGen::DDVCSEventGen()
{
  fPairCAngle = 0 * deg;
  fPairPhiAngle = 0* deg;
  fPairDThetaAngle = 0* deg;
  fPairE= 1 * GeV;
  fToFres = 0.5*ns;
  fQPrime2 = 2.5*GeV*GeV;
}



DDVCSEventGen::~DDVCSEventGen()
{

}

void DDVCSEventGen::GeneratePair()
{
  G4ThreeVector pairax,pairm,pairp;
   G4ThreeVector rotax;
   G4RotationMatrix pairrot;
   pairm.setRThetaPhi(fPairE,(fPairCAngle+fPairDThetaAngle),fPairPhiAngle);
   pairp.setRThetaPhi(fPairE,(fPairCAngle-fPairDThetaAngle),-fPairPhiAngle);
   rotax.setRThetaPhi(fPairE,(fPairCAngle),0);
   printf("Rotate by %f \n",fPairRotAngle);
   pairrot=pairrot.rotate(fPairRotAngle,rotax);
   printf("Before rotation theta : %f Phi : %f\n",pairm.theta(),pairm.phi());
   pairm = pairrot * pairm;
   printf("After rotation theta : %f Phi : %f\n",pairm.theta(),pairm.phi());
   pairp = pairrot * pairp;
   fQP.setRThetaPhi(fPairE,pairp.theta(),pairp.phi());
   fQM.setRThetaPhi(fPairE,pairm.theta(),pairm.phi());

}
