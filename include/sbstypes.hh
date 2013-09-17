#ifndef SBSTYPES_HH
#define SBSTYPES_HH

#include "globals.hh"

#define MAXTARG 5
#define NNUCL   2

// Include fictional neutron target
enum Targ_t { kH2, kLH2, kLD2, k3He, kNeutTarg };
enum Nucl_t { kProton, kNeutron };
enum Pion_t { kpip, kpim, kpi0 };
enum Kine_t { kElastic, kFlat, kInelastic, kDIS, kBeam };
enum Exp_t { kGEp, kNeutronExp };

#endif//SBSTYPES_HH
