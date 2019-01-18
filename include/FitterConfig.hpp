//
//  FitterConfig.hpp
//
//  Copyright © 2019 Jeremi Niedziela. All rights reserved.
//

#ifndef FitterConfig_hpp
#define FitterConfig_hpp

#include "Helpers.hpp"

class FitterConfig {
public:
  FitterConfig(string _path);
  
  void Print();
  
  // Fitter parameters
  double GetHelixThickness(){return helixThickness;}
  double GetCircleThickness(){return circleThickness;}
  double GetLinesToleranceForCircles(){return linesToleranceForCircles;}
  double GetLinesToleranceForRegularity(){return linesToleranceForRegularity;}
  double GetStepPz(){return stepPz;}
  double GetZregularityTolerance(){return zRegularityTolerance;}
  int GetMinPointsAlongZ(){return minNpointsAlongZ;}
  
  // Random chargino's track parameters
  double GetMaxTrackEta(){return maxEta;}
  double GetMinL(){return minL;}
  double GetMaxL(){return maxL;}
  
  // Random pion's parameters
  double GetMinPx(){return minPx;}
  double GetMinPy(){return minPy;}
  double GetMinPz(){return minPz;}
  double GetMaxPx(){return maxPx;}
  double GetMaxPy(){return maxPy;}
  double GetMaxPz(){return maxPz;}
  
  // General settings
  bool GetInjectPionHits(){return injectPionHits;}
  int GetNtests(){return nTests;}
  const char* GetOutputPath(){return outputPath;}
  int GetNnoiseHits(){return nNoiseHits;}
  
  // Benchmark parameters
  double GetToleranceX(){return toleranceX;}
  double GetToleranceY(){return toleranceY;}
  double GetToleranceZ(){return toleranceZ;}
  double GetTolerancePx(){return tolerancePx;}
  double GetTolerancePy(){return tolerancePy;}
  double GetTolerancePz(){return tolerancePz;}

//private:
  unique_ptr<TEnv> config;
  
  double helixThickness;
  double circleThickness;
  double linesToleranceForCircles;
  double linesToleranceForRegularity;
  double stepPz;
  double zRegularityTolerance;
  int minNpointsAlongZ;
  double maxEta;
  double minL;
  double maxL;
  double minPx;
  double minPy;
  double minPz;
  double maxPx;
  double maxPy;
  double maxPz;
  bool injectPionHits;
  int nTests;
  const char* outputPath;
  double toleranceX;
  double toleranceY;
  double toleranceZ;
  double tolerancePx;
  double tolerancePy;
  double tolerancePz;
  int nNoiseHits;
};

#endif /* FitterConfig_hpp */
