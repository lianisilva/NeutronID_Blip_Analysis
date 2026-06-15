// ROOT includes
#include "TH1.h"

// XSecAnalyzer includes
#include "XSecAnalyzer/Selections/EventCategoriesBlip_n_Np.hh"

std::map< int, std::pair< std::string, int > > Blip_n_Np_MAP = {
  { kUnknown, { "Unknown", kGray } },
  { kNuMuCC_prim_n_1p, { "CCmu1pNn", kGreen } },
  { kNuMuCC_prim_n_Np, { "CCmuNpNn", kGreen + 1 } }, 
  { kNuMuCC_second_n_Xp, { "secondary (n,Xp)", kOrange + 4 } },
  { kNuMuCC_prim_n_gamma, { "primary (n,gamma)", kCyan - 3 } },
  { kNuMuCC_second_n_gamma, { "secondary (n,gamma)", kCyan - 6 } },
  { kNuMuCC_ncap_gamma, { "ncapture gamma", kCyan - 4 } },
  { kNuMuCC_cosmics_rads, { "cosmics & radiogenics", kBlue - 2 } },

  { kNuMuCCOther, { "Other #nu_{#mu} CC", kAzure } },
  { kNuECC, { "#nu_{e} CC", kViolet } },
  { kNC, { "NC", kOrange } },
  { kOOFV, {"Out FV", kRed + 3 } },
  { kOther, { "Other", kRed + 1 } }
};

