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
const int performCutsLevel = 2;

const bool saveEvents = false;
const bool printYields = true;
const bool printBackgroundDetails = false;

const bool drawStandardPlots = false;
const bool drawPerLayerPlots = false;
const bool showLegends = false;

// Limit number of events loaded (-1 means load all available)
const int maxNeventsBackground  = -1;
const int maxNeventsSignal      = -1;
const int maxNeventsData        = -1;

// turn on/off different backgrounds, signals and data samples
const vector<bool> runBackground = {
  false,   // QCD
  true,   // Z->μμ + jets
  false,   // tops
  false,   // VV
  true,   // W->μν + jets
  true,   // Z->νν + jets
};

const vector<bool> runSignal = {
  false,   // wino m=300 cτ=3
  false,   // wino m=300 cτ=10
  false,   // wino m=300 cτ=30
  false,   // wino m=500 cτ=10
  false,   // wino m=500 cτ=20
  false,   // wino m=650 cτ=10
  false,   // wino m=650 cτ=20
  false,   // wino m=800 cτ=10
  false,   // wino m=800 cτ=20
  false,   // wino m=1000 cτ=10
  false,   // wino m=1000 cτ=20
};

const vector<bool> runData = {
  false,  // 2017B
};

#endif /* Config_h */
