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
const bool printDataDetails = false;
const bool printSignalDetails = false;

const bool drawStandardPlots = false;
const bool drawPerLayerPlots = false;
const bool showLegends = false;
const bool scanMETbinning = false;
const bool doMETbinning = false;

// Limit number of events loaded (-1 means load all available)
const int maxNeventsBackground  = -1;
const int maxNeventsSignal      = -1;
const int maxNeventsData        = -1;

//                             2015   2016    2017    2018
const double totalLuminosity = 3.81 + 37.76 + 41.37 + 63.97; // in fb^-1
//const double totalLuminosity = 41.37; // in fb^-1

enum ECategory{
  k2tracks,
  k3layers,
  k4layers
};

const ECategory category = k3layers;

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
  true,  // 2017 (B+C+D)
};

#endif /* Config_h */
