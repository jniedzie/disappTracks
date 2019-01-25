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
const int performCutsLevel = 20;

const bool saveEvents = false;
const bool printYields = true;
const bool printBackgroundDetails = false;
const bool printDataDetails = false;
const bool printSignalDetails = false;

const bool drawStandardPlots = true;
const bool drawPerLayerPlots = true;
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

const ECategory category = k4layers;

// turn on/off different backgrounds, signals and data samples
const vector<bool> runBackground = {
  false,   // QCD
  false,   // Z->μμ + jets
  false,   // tops
  false,   // VV
  false,   // W->μν + jets
  false,   // Z->νν + jets
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
  true,  // 2017 (B+C+D)
};

#endif /* Config_h */
