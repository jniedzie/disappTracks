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
const bool printYields = true;
const bool interactive = false;
const bool showLegends = false;
const int performCutsLevel = 2;

// Limit number of events loaded (-1 means load all available)
const int maxNeventsBackground  = -1;
const int maxNeventsSignal      = -1;
const int maxNeventsData        = -1;

// turn on/off different backgrounds, signals and data samples
const vector<bool> runBackground = {
  true,   // QCD
  true,   // Z->μμ + jets
  true,   // tops
  true,   // VV
  true,   // W->μν + jets
  true,   // Z->νν + jets
};

const vector<bool> runSignal = {
  true,   // wino m=300 cτ=3
  true,   // wino m=300 cτ=10
  true,   // wino m=300 cτ=30
  true,   // wino m=500 cτ=10
  true,   // wino m=500 cτ=20
  true,   // wino m=650 cτ=10
  true,   // wino m=650 cτ=20
  true,   // wino m=800 cτ=10
  true,   // wino m=800 cτ=20
  true,   // wino m=1000 cτ=10
  true,   // wino m=1000 cτ=20
};

const vector<bool> runData = {
  false,  // 2017B
};

#endif /* Config_h */
