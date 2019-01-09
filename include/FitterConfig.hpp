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
  
  // Fitter parameters
  double GetHelixThickness();
  double GetCircleThickness();
  double GetLinesTolerance();
  double GetStepPz();
  double GetZregularityTolerance();
  int GetMinPointsAlongZ();
  
  // Random chargino's track parameters
  double GetMaxTrackEta();
  double GetMinL();
  double GetMaxL();
  
  // Random pion's parameters
  double GetMinPx();
  double GetMinPy();
  double GetMinPz();
  double GetMaxPx();
  double GetMaxPy();
  double GetMaxPz();
  
  // General settings
  bool GetInjectPionHits();
  int GetNtests();
  const char* GetOutputPath();
  int GetNnoiseHits();
  
  // Benchmark parameters
  double GetToleranceX();
  double GetToleranceY();
  double GetToleranceZ();
  double GetTolerancePx();
  double GetTolerancePy();
  double GetTolerancePz();

private:
  unique_ptr<TEnv> config;
};

#endif /* FitterConfig_hpp */
