#ifndef _DDVCSEventGen_
#define _DDVCSEventGen_
#include "G4SBSEventGen.hh"
#include "G4RotationMatrix.hh"

class DDVCSEventGen : public G4SBSEventGen
{
public:
  DDVCSEventGen();
  ~DDVCSEventGen();
  G4LorentzVector GetQP(){return fQP;}
  G4LorentzVector GetQM(){return fQM;}
   void SetPairE(double v){fPairE = v;}
        void SetQprime2(double v){fQPrime2 = v;}
        void SetPairPart(G4String v){fPairPart = v;}
        void SetPairCAngle(double v){fPairCAngle = v;}
        void SetPairRotAngle(double v){fPairRotAngle = v;}
        void SetPairDThetaAngle(double v){fPairDThetaAngle = v;}
        void SetPairPhiAngle(double v){fPairPhiAngle = v;}
        G4String GetPairPart(){return fPairPart;}
        double GetPairCAngle(){return fPairCAngle;}
        double GetPairPhiAngle(){return fPairPhiAngle;}
        double GetPairDThetaAngle(){return fPairDThetaAngle;}
        double GetPairRotAngle(){return fPairRotAngle;}
  void GeneratePair();
protected:
  double  fPairE, fPairCAngle,fPairPhiAngle,fPairRotAngle,fPairDThetaAngle,fQPrime2;
    G4LorentzVector fQM,fQP;
    G4String fPairPart; 
}
  ;
#endif
