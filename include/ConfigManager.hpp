//
//  FitterConfig.hpp
//
//  Copyright © 2019 Jeremi Niedziela. All rights reserved.
//

#ifndef FitterConfig_hpp
#define FitterConfig_hpp

#include "Helpers.hpp"

struct ConfigManager;
extern ConfigManager config;

/// Wrapper on a config file that provides access to options from the code.
/// Reads a config file in markdown format. For all options that are not specified in the config
/// default values will be used.
struct ConfigManager {
  /// Default constructor
  /// \_path Path to the config file
  ConfigManager(string _path="");
  
  /// Prints complete information about current configuration
  void Print();
  
  map<string, double> params;
  
  //------------------------------//
  //    Analysis configuration    //
  //------------------------------//
  
  int performCutsLevel = 2;
  
  bool saveEvents;
  bool printYields;
  bool printBackgroundDetails;
  bool printDataDetails;
  bool printSignalDetails;
  
  bool drawStandardPlots;
  bool drawPerLayerPlots;
  bool showLegends;
  bool scanMETbinning;
  bool doMETbinning;
  
  int maxNeventsBackground;
  int maxNeventsSignal;
  int maxNeventsData;
  
  bool loadFriendTree;
  
  double totalLuminosity;
  
  string category;
  
  vector<bool> runBackground;
  vector<bool> runSignal;
  vector<bool> runData;
  
  //------------------------------//
  //  Helix fitter configuration  //
  //------------------------------//
  
  // General settings
  bool injectPionHits;
  int nTests;
  int nNoiseHits;
  int nTrackerLayers;
  int verbosity;
  
  // Fitter parameter
  
  double startR0, minR0, maxR0;
  double minRslope, maxRslope;
  
  double minS0, maxS0;
  double minSslope, maxSslope;
  double minX0, maxX0;
  double minY0, maxY0;
  double minZ0, maxZ0;
  
  // Random chargino's track parameters
  double maxEta;
  int nTrackHits;
  
    // Random pion's parameters
  double minPx;
  double minPy;
  double minPz;
  double maxPx;
  double maxPy;
  double maxPz;
  
  const char* outputPath;
  
  // Benchmark parameters
  double toleranceX;
  double toleranceY;
  double toleranceZ;
  double tolerancePx;
  double tolerancePy;
  double tolerancePz;
  
  //------------------------------//
  // Event Display configuration  //
  //------------------------------//
  
  bool showGeometryPixel;
  bool showGeometryStrip;
  bool showGeometryEcal;
  bool showGeometryHcal;
  
  bool drawTrackerClusters;
  bool drawMET;
  bool drawJets;
  bool drawPionSimHits;
  bool drawPionClusters;
  bool drawCharginoSimHits;
  
private:
  unique_ptr<TEnv> configFile;
};

#endif /* FitterConfig_hpp */
