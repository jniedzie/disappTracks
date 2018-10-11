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

// Limit number of events loaded (-1 means - load all available)
const int maxNeventsBackground  = 10000;
const int maxNeventsSignal      = -1;
const int maxNeventsData        = -1;

// turn on/off different backgrounds, signals and data samples
enum EBackground{
  kZmumuJets,
  kTT,
  kVV,
  kWmunuJets,
  kQCD,
  kZnunuJets,
  kNbackgrounds
};

enum ESignal{
  kWino_M_300_cTau_10,
  //  kWino_M_300_cTau_3,
  //  kWino_M_300_cTau_30,
  //  kWino_M_500_cTau_10,
  kNsignals
};

enum EData{
//  kElectron_Run2017B,
  kNdata
};

#endif /* Config_h */
