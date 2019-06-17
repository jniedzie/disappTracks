//  PerformaceMonitor.hpp
//
//  Created by Jeremi Niedziela on 17/06/2019.

#ifndef PerformanceMonitor_hpp
#define PerformanceMonitor_hpp

#include "Helpers.hpp"

class PerformanceMonitor
{
public:
  PerformanceMonitor();
  
  PerformanceMonitor(string _name, int _nBins, double min, double max,
                     int _nTests, int nEvents);
  
  void operator=(const PerformanceMonitor &pm);
  
  void SetValues(int iTest, int iEvent, double valueSignal, double valueBackground);
  
  void CalcEfficiency(function<double(int)> SetParamValue, int nAnalyzedEvents);
  
  void DrawRocGraphs();
  void DrawHists();
  
private:
  int nTests;
  
  double thresholdMin;
  double thresholdMax;
  double thresholdStep;
  int nBins;
  
  vector<vector<double>> valuesSignal;// [iTest][iEvent]
  vector<vector<double>> valuesBackground;
  
  vector<double> efficiency;
  vector<double> fakeRate;
  
  TH1D *histSignal;
  TH1D *histBackground;
  
  vector<TF1*> rocFun;
  vector<TGraph*> rocGraph;
  
  string name;
  
  
  TF1* GetRocFunction();
};

#endif /* PerformanceMonitor_hpp */
