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
  
  PerformanceMonitor(string _name, int _nBins, double min, double max);
  
  void operator=(const PerformanceMonitor &pm);
  
  void SetValues(double valueSignal, double valueBackground);
  
  void CalcEfficiency();
  
  void DrawRocGraph(bool first);
  void DrawHists();
  
  void PrintFakesEfficiency();
  void PrintParams();
  
  // trivial getters
  
  inline double GetValueByName(string name){
    if(name=="auc") return auc;
    if(name=="sigma_init") return significanceInitial;
    if(name=="sigma_L0") return significanceAfterL0;
    if(name=="sigma_L1") return significanceAfterL1;
    if(name=="max_eff") return maxEfficiency;
    if(name=="min_fake") return invFakeAtHighestEff;
    if(name=="max_dist_fake") return maxDistToSqrtFake;
    if(name=="avg_dist_fake") return avgDistToSqrtFake;
    
    return -inf;
  }
  
  inline double GetAUC(){ return auc; }
  inline double GetMaxEfficiency(){ return maxEfficiency; }
  inline double GetSignificanceInitial(){ return significanceInitial; }
  inline double GetSignificanceAfterL0(){ return significanceAfterL0; }
  inline double GetSignificanceAfterL1(){ return significanceAfterL1; }
  inline double GetInvFakeAtHighestEff(){ return invFakeAtHighestEff; }
  inline double GetMaxDistToSqrtFake(){ return maxDistToSqrtFake; }
  inline double GetAvgDistToSqrtFake(){ return avgDistToSqrtFake; }
  
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
  double maxDistToSqrtFake;   ///< (max) Distance to √c_fake, which is a minimum to be useful
  double avgDistToSqrtFake;   ///< (avg) Distance to √c_fake, which is a minimum to be useful
  
  TF1* GetRocFunction();
};

#endif /* PerformanceMonitor_hpp */
