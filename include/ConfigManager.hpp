//
//  FitterConfig.hpp
//
//  Copyright © 2019 Jeremi Niedziela. All rights reserved.
//

#ifndef FitterConfig_hpp
#define FitterConfig_hpp

#include "Helpers.hpp"

struct ConfigManager {
  ConfigManager(string _path="");
  void Print();
  
  // General settings
  bool injectPionHits;
  int nTests;
  int nNoiseHits;
  int nTrackerLayers;
  
  // Fitter parameters
  double helixThickness;
  double circleThickness;
  double linesToleranceForCircles;
  double linesToleranceForRegularity;
  double stepPz;
  double zRegularityTolerance;
  int minNpointsAlongZ;
  
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
  
private:
  unique_ptr<TEnv> config;
};

#endif /* FitterConfig_hpp */
