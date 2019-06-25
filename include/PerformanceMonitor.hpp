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
  
  PerformanceMonitor(string _name, int _nBins, double min, double max, int nEvents);
  
  void operator=(const PerformanceMonitor &pm);
  
  void SetValues(int iEvent, double valueSignal, double valueBackground);
  
  void CalcEfficiency(int nAnalyzedEvents);
  
  void DrawRocGraph(bool first);
  void DrawHists();
  
  void PrintFakesEfficiency();
  void PrintParams();
  
  // trivial getters
  inline double GetAUC(){ return auc; }
  inline double GetMaxEfficiency(){ return maxEfficiency; }
  inline double GetSignificanceInitial(){ return significanceInitial; }
  inline double GetSignificanceAfterL0(){ return significanceAfterL0; }
  inline double GetSignificanceAfterL1(){ return significanceAfterL1; }
  inline double GetInvFakeAtHighestEff(){ return invFakeAtHighestEff; }
  
private:
  double thresholdMin;
  double thresholdMax;
  double thresholdStep;
  int nBins;
  
  vector<double> valuesSignal;// [iEvent]
  vector<double> valuesBackground;
  
  TH1D *histSignal;
  TH1D *histBackground;
  
  TF1* rocFun;
  TGraph* rocGraph;
  
  string name;
  
  vector<double> efficiency;  ///< efficiency for given threshold
  vector<double> fakeRate;    ///< fake rate for given threshold
  
  double auc;                 ///< area under ROC curve
  double maxEfficiency;       ///< maximum efficiency lower than 1.0
  double significanceInitial; ///< max significance assuming initial N_sig and N_bck, only when fake rate !=0
  double significanceAfterL0; ///< max significance assuming L0 N_sig and N_bck, only when fake rate !=0
  double significanceAfterL1; ///< max significance assuming L1 N_sig and N_bck, only when fake rate !=0
  double invFakeAtHighestEff; ///< 1/fake_rate for the highest efficiency lower than 1.0
  
  TF1* GetRocFunction();
};

#endif /* PerformanceMonitor_hpp */
