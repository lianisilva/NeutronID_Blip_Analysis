#pragma once

// Standard library includes
#include <map>
#include <string>

// Enum used to label event categories of interest for analysis plots in
// the CC0pi 1p/2p/Np/Xp analyses
enum EventCategoryBlip_n_Np {

  // Unable to categorize (e.g., because the event is real data and thus
  // has no MC truth information)
  kUnknown = 0,

  // Signal events broken down by number of protons
  kNuMuCC_prim_n_1p = 1,
  kNuMuCC_prim_n_Np = 2,

  kNuMuCC_second_n_Xp = 3,
  kNuMuCC_prim_n_gamma = 4,
  kNuMuCC_second_n_gamma = 5,
  kNuMuCC_ncap_gamma = 6,
  kNuMuCC_cosmics_rads = 7,   // this might go in unknown

  // Any true numu CC event which does not satisfy the criteria for inclusion
  // in one of the other categories above
  kNuMuCCOther = 8,

  // True nue CC event
  kNuECC = 9,

  // True neutral current event for any neutrino flavor
  kNC = 10,

  // True neutrino vertex (any reaction mode and flavor combination) is outside
  // of the fiducial volume
  kOOFV = 11,

  // All events that do not fall within any of the other categories (e.g.,
  // numubar CC)
  kOther = 12,

};

extern std::map< int, std::pair< std::string, int > > Blip_n_Np_MAP;

