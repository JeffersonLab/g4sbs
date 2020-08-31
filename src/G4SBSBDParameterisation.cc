#include "G4SBSBDParameterisation.hh"
//______________________________________________________________________________
G4SBSBDParameterisation::G4SBSBDParameterisation(char Hall,G4ThreeVector r0){
   fHall    = Hall; // hall label 
   fR0      = r0;   // origin of part (relative to mother volume) 
   InitParameters();
}
//______________________________________________________________________________
G4SBSBDParameterisation::~G4SBSBDParameterisation(){
   delete[] fThickness;
   delete[] fStartPhi;
   delete[] fDeltaPhi;
   delete[] fColor;
}
//______________________________________________________________________________
void G4SBSBDParameterisation::InitParameters(){
   // check which hall we are setting up for 
   if(fHall!='A' && fHall!='C'){
      std::cout << "[G4SBSBDParameterisation::InitParameters]: Invalid Hall = "
                << fHall << std::endl;
      exit(1);
   }

   // initialize materials  
   auto nistManager = G4NistManager::Instance();
   nistManager->FindOrBuildMaterial("G4_Al");

   // define some constants 
   G4double inch        = 25.4*mm;
   G4double thk_031     = 0.03125*inch;
   G4double thk_062     = 0.06250*inch;
   G4double thk_094     = 0.09375*inch;
   G4double thk_100     = 0.10000*inch;
   G4double thk_125     = 0.12500*inch;
   G4double startPhi_30 = 255.*deg;
   G4double startPhi_60 = 255.*deg;
   G4double startPhi_90 = 225.*deg;
   G4double deltaPhi_30 = 30.*deg;
   G4double deltaPhi_60 = 60.*deg;
   G4double deltaPhi_90 = 90.*deg;

   fGap        = 0.195*inch;
   fWidth      = 2.*inch;               // FIXME: Estimates! 
   fRadius_min = 5.*inch;               // FIXME: Estimates!  
   fRadius_max = fRadius_min + fWidth;

   // make arrays large enough for all the layers we may need 
   const int NP_MAX = 16;
   fThickness = new double[NP_MAX];
   fStartPhi  = new double[NP_MAX];
   fDeltaPhi  = new double[NP_MAX];
   fColor     = new int[NP_MAX];
   for(int i=0;i<NP_MAX;i++){
      fThickness[i] = 0;
      fStartPhi[i]  = 0;
      fDeltaPhi[i]  = 0;
      fColor[i]     = 0;
   }

   // set thickness array 
   if(fHall=='A'){
      fNPlanes = 15;
      // thicknesses  
      fThickness[0]  = thk_125; fThickness[1]  = thk_125; fThickness[2]  = thk_100; fThickness[3]  = thk_125;
      fThickness[4]  = thk_125; fThickness[5]  = thk_100; fThickness[6]  = thk_125; fThickness[7]  = thk_100;
      fThickness[8]  = thk_125; fThickness[9]  = thk_100; fThickness[10] = thk_125; fThickness[11] = thk_125;
      fThickness[12] = thk_100; fThickness[13] = thk_125; fThickness[14] = thk_125;
      // start angles 
      fStartPhi[0]   = startPhi_30; fStartPhi[1]  = startPhi_30; fStartPhi[2]  = startPhi_30; fStartPhi[3]  = startPhi_60;
      fStartPhi[4]   = startPhi_60; fStartPhi[5]  = startPhi_60; fStartPhi[6]  = startPhi_90; fStartPhi[7]  = startPhi_90;
      fStartPhi[8]   = startPhi_90; fStartPhi[9]  = startPhi_60; fStartPhi[10] = startPhi_60; fStartPhi[11] = startPhi_60;
      fStartPhi[12]  = startPhi_30; fStartPhi[13] = startPhi_30; fStartPhi[14] = startPhi_30;
      // step angles 
      fDeltaPhi[0]   = deltaPhi_30; fDeltaPhi[1]  = deltaPhi_30; fDeltaPhi[2]  = deltaPhi_30; fDeltaPhi[3]  = deltaPhi_60;
      fDeltaPhi[4]   = deltaPhi_60; fDeltaPhi[5]  = deltaPhi_60; fDeltaPhi[6]  = deltaPhi_90; fDeltaPhi[7]  = deltaPhi_90;
      fDeltaPhi[8]   = deltaPhi_90; fDeltaPhi[9]  = deltaPhi_60; fDeltaPhi[10] = deltaPhi_60; fDeltaPhi[11] = deltaPhi_60;
      fDeltaPhi[12]  = deltaPhi_30; fDeltaPhi[13] = deltaPhi_30; fDeltaPhi[14] = deltaPhi_30;
      // colors 
      fColor[0]      = diffuser::kBlue;    fColor[1]  = diffuser::kBlue;    fColor[2]  = diffuser::kMagenta; fColor[3]  = diffuser::kBlue;
      fColor[4]      = diffuser::kBlue;    fColor[5]  = diffuser::kMagenta; fColor[6]  = diffuser::kBlue;    fColor[7]  = diffuser::kMagenta;
      fColor[8]      = diffuser::kBlue;    fColor[9]  = diffuser::kMagenta; fColor[10] = diffuser::kBlue;    fColor[11] = diffuser::kBlue;
      fColor[12]     = diffuser::kMagenta; fColor[13] = diffuser::kBlue;    fColor[14] = diffuser::kBlue;
   }else if(fHall=='C'){
      fNPlanes = 16;
      // thicknesses  
      fThickness[0]  = thk_125; fThickness[1]  = thk_125; fThickness[2]  = thk_125; fThickness[3]  = thk_125;
      fThickness[4]  = thk_125; fThickness[5]  = thk_125; fThickness[6]  = thk_094; fThickness[7]  = thk_094;
      fThickness[8]  = thk_094; fThickness[9]  = thk_094; fThickness[10] = thk_094; fThickness[11] = thk_094;
      fThickness[12] = thk_062; fThickness[13] = thk_062; fThickness[14] = thk_031; fThickness[15] = thk_031;
      // start angles 
      fStartPhi[0]   = startPhi_30; fStartPhi[1]  = startPhi_30; fStartPhi[2]  = startPhi_30; fStartPhi[3]  = startPhi_30;
      fStartPhi[4]   = startPhi_60; fStartPhi[5]  = startPhi_60; fStartPhi[6]  = startPhi_90; fStartPhi[7]  = startPhi_90;
      fStartPhi[8]   = startPhi_90; fStartPhi[9]  = startPhi_60; fStartPhi[10] = startPhi_60; fStartPhi[11] = startPhi_60;
      fStartPhi[12]  = startPhi_30; fStartPhi[13] = startPhi_30; fStartPhi[14] = startPhi_30; fStartPhi[15] = startPhi_30;
      // step angles 
      fDeltaPhi[0]   = deltaPhi_30; fDeltaPhi[1]  = deltaPhi_30; fDeltaPhi[2]  = deltaPhi_30; fDeltaPhi[3]  = deltaPhi_30;
      fDeltaPhi[4]   = deltaPhi_60; fDeltaPhi[5]  = deltaPhi_60; fDeltaPhi[6]  = deltaPhi_90; fDeltaPhi[7]  = deltaPhi_90;
      fDeltaPhi[8]   = deltaPhi_90; fDeltaPhi[9]  = deltaPhi_60; fDeltaPhi[10] = deltaPhi_60; fDeltaPhi[11] = deltaPhi_60;
      fDeltaPhi[12]  = deltaPhi_30; fDeltaPhi[13] = deltaPhi_30; fDeltaPhi[14] = deltaPhi_30; fDeltaPhi[15] = deltaPhi_30;
      // colors 
      fColor[0]      = diffuser::kBlue;   fColor[1]  = diffuser::kBlue;   fColor[2]  = diffuser::kBlue;  fColor[3]  = diffuser::kBlue;
      fColor[4]      = diffuser::kBlue;   fColor[5]  = diffuser::kBlue;   fColor[6]  = diffuser::kGreen; fColor[7]  = diffuser::kGreen;
      fColor[8]      = diffuser::kGreen;  fColor[9]  = diffuser::kGreen;  fColor[10] = diffuser::kGreen; fColor[11] = diffuser::kGreen;
      fColor[12]     = diffuser::kYellow; fColor[13] = diffuser::kYellow; fColor[14] = diffuser::kRed;   fColor[15] = diffuser::kRed;
   }
   // total thickness (used for placement of layers) 
   fTotalThickness = fGap*( (double)fNPlanes - 1. );
   for(int i=0;i<fNPlanes;i++) fTotalThickness += fThickness[i];
}
//______________________________________________________________________________
void G4SBSBDParameterisation::ComputeTransformation(const G4int copyNo,
      G4VPhysicalVolume *physVol) const{
   // z coordinate: distance along z to *center* of the plate  
   // sum over previous layers 
   G4double Ls=0;
   for(int i=0;i<copyNo;i++) Ls += fThickness[i];
   // put it all together 
   G4double xp = fR0.x();
   G4double yp = fR0.y() + fRadius_min + 0.5*fWidth;
   G4double z  = fR0.z();
   G4double zp = -0.5*fTotalThickness + z + Ls + (double)(copyNo-1)*fGap + 0.5*fThickness[copyNo];
   // set the 3-vector
   G4ThreeVector P = G4ThreeVector(xp,yp,zp);
   physVol->SetTranslation(P); // set position 
   physVol->SetRotation(0);    // no rotation  
}
//______________________________________________________________________________
void G4SBSBDParameterisation::ComputeDimensions(G4Tubs &plate,
      const G4int copyNo,const G4VPhysicalVolume *physVol) const{
   // set dimensions
   plate.SetInnerRadius(fRadius_min);
   plate.SetOuterRadius(fRadius_max);
   plate.SetZHalfLength(fThickness[copyNo]/2.);
   plate.SetStartPhiAngle(fStartPhi[copyNo]);
   plate.SetDeltaPhiAngle(fDeltaPhi[copyNo]);
   // determine color by copy number 
   G4VisAttributes *vis = new G4VisAttributes();
   if(fColor[copyNo]==diffuser::kRed    ) vis->SetColour( G4Colour::Red()     );
   if(fColor[copyNo]==diffuser::kYellow ) vis->SetColour( G4Colour::Yellow()  );
   if(fColor[copyNo]==diffuser::kGreen  ) vis->SetColour( G4Colour::Green()   );
   if(fColor[copyNo]==diffuser::kBlue   ) vis->SetColour( G4Colour::Blue()    );
   if(fColor[copyNo]==diffuser::kMagenta) vis->SetColour( G4Colour::Magenta() );
   // set properties 
   physVol->GetLogicalVolume()->SetVisAttributes(vis);
}
