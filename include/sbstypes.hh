#ifndef SBSTYPES_HH
#define SBSTYPES_HH

#include "TTimeStamp.h"

#define MAXTARG 5
#define NNUCL   2

#define __RUNSTR_LEN 255
#define __MAXFILE_LEN 1048576 // MB file size limit

// Include fictional neutron target
enum Targ_t { kH2, kLH2, kLD2, k3He, kNeutTarg };
enum Nucl_t { kProton, kNeutron };
enum Hadron_t { kPiPlus, kPiMinus, kPi0, kKPlus, kKMinus }; //Hadron types for SIDIS event generator
enum Kine_t { kElastic, kFlat, kInelastic, kDIS, kBeam, kSIDIS };
enum Exp_t { kGEp, kNeutronExp, kSIDISExp };

struct filedata_t {
    char filename[__RUNSTR_LEN];
    char hashsum[__RUNSTR_LEN];
    TTimeStamp timestamp;
};


#endif//SBSTYPES_HH
