#ifndef SBSTYPES_HH
#define SBSTYPES_HH

#include "globals.hh"

#define MAXTARG 5
#define NNUCL   2

// Include fictional neutron target
enum Targ_t { kH2, kLH2, kLD2, k3He, kNeutTarg };
enum Nucl_t { kProton, kNeutron };
enum Hadron_t { kPiPlus, kPiMinus, kPi0, kKPlus, kKMinus }; //Hadron types for SIDIS event generator
enum Kine_t { kElastic, kFlat, kInelastic, kDIS, kBeam, kSIDIS };
enum Exp_t { kGEp, kNeutronExp, kSIDISExp };

#endif//SBSTYPES_HH
