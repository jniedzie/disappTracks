//
//  Config.hpp
//
//  Created by Jeremi Niedziela on 09/10/2018.
//

#ifndef Config_h
#define Config_h

#include <vector>
#include <string>

using namespace std;

// Analysis configuration
const bool printHeaders = true;

// Limit number of events loaded (-1 means load all available)
const int maxNeventsBackground  = 10000;
const int maxNeventsSignal      = -1;
const int maxNeventsData        = 10000;

// turn on/off different backgrounds, signals and data samples
enum EBackground{
  kQCD,
  kZmumuJets,
  kTT,
  kVV,
  kWmunuJets,
  kZnunuJets,
  kNbackgrounds
};

enum ESignal{
  kWino_M_300_cTau_3,
  kWino_M_300_cTau_10,
  kWino_M_300_cTau_30,
  kWino_M_500_cTau_10,
  kWino_M_500_cTau_20,
  kWino_M_650_cTau_10,
  kWino_M_650_cTau_20,
  kWino_M_800_cTau_10,
  kWino_M_800_cTau_20,
  kWino_M_1000_cTau_10,
  kWino_M_1000_cTau_20,
  kNsignals
};

const bool runSignal[kNsignals] = {
  false,  // wino m=300 cτ=3
  false,  // wino m=300 cτ=10
  true,   // wino m=300 cτ=30
  false,  // wino m=500 cτ=10
  false,  // wino m=500 cτ=20
  false,  // wino m=650 cτ=10
  false,  // wino m=650 cτ=20
  false,  // wino m=800 cτ=10
  false,  // wino m=800 cτ=20
  false,  // wino m=1000 cτ=10
  false,  // wino m=1000 cτ=20
};

enum EData{
  kElectron_Run2017B,
  kNdata
};

#endif /* Config_h */
