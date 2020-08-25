#ifndef SBSTYPES_HH
#define SBSTYPES_HH

#include "TTimeStamp.h"

#define MAXTARG 5
#define NNUCL   2

#define __RUNSTR_LEN 255
#define __MAXFILE_LEN 1048576 // MB file size limit

// #include "TPDGCode.h"  // kind of works, but you have to change function definitions to accept PDG_t types 

// encapsulate in a namespace so we don't step on the ROOT and GEANT4 definitions 
namespace G4SBS { 
   // particle type definitions (matches the GEANT4 standard) 
   enum Nucl_t { kProton, kNeutron };
   enum Hadron_t { kPiPlus, kPiMinus, kPi0, kKPlus, kKMinus, kP, kPbar}; //Hadron types for SIDIS event generator
   // target type; include fictional neutron target
   enum Targ_t { kH2, kD2, kLH2, kLD2, k3He, kNeutTarg, kCfoil, kOptics };
   // kinematic type 
   enum Kine_t { kElastic, kFlat, kInelastic, kDIS, kBeam, kSIDIS, kGun, kWiser, kPYTHIA6, kGMnElasticCheck, kCosmics, kPionPhoto};
   //enum Exp_t { kGEp, kNeutronExp, kSIDISExp, kC16, kA1n, kTDIS, kNDVCS, kGEnRP, kGEMHCtest};
   // experiment type
   enum Exp_t { kGEp, kGMN, kGEN, kSIDISExp, kC16, kA1n, kTDIS, kNDVCS, kGEnRP, kGEMHCtest, kGEPpositron, kWAPP};
   // detector arm type (for association of detector modules with spectrometer arms. Presently "E arm" and "H arm" are possible) 
   enum Arm_t { kEarm, kHarm };
   // sensitive detector type  
   enum SDet_t { kGEM, kCAL, kRICH, kECAL, kBD }; 

   struct filedata_t {
      char filename[__RUNSTR_LEN];
      char hashsum[__RUNSTR_LEN];
      TTimeStamp timestamp;
   };

   // switches for GEn
   // Helmholtz coils or shielding 
   // - 146  => Q2 = 1.46  (GeV/c)^2  
   // - 368  => Q2 = 3.68  (GeV/c)^2  
   // - 677  => Q2 = 6.77  (GeV/c)^2  
   // - 1018 => Q2 = 10.18 (GeV/c)^2
   // for shielding only:  
   // - full => full window cut (remove panel 1, 2, 3, and door) 
   // - new  => new design from Bert Metzger (6/2020)  
   enum GEnConfig_t {
      kSBS_GEN_146  = 1,
      kSBS_GEN_368  = 2,
      kSBS_GEN_677  = 3,
      kSBS_GEN_1018 = 4,
      kSBS_GEN_full = 5,
      kSBS_GEN_new  = 6
   };

} 

#endif//SBSTYPES_HH
