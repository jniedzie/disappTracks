//  PerformaceMonitor.hpp
//
//  Created by Jeremi Niedziela on 17/06/2019.

#ifndef PerformanceMonitor_hpp
#define PerformanceMonitor_hpp

#include "Helpers.hpp"

class PerformanceMonitor;
typedef map<string, PerformanceMonitor> Monitors;

/**
 Monitors performance of the looper tagger, by storing values of some discriminating
 variable (e.g. number of hits on helix) and calculating efficiency and fake rate.
 */
class PerformanceMonitor
{
public:
  /// Default constructor
  PerformanceMonitor();
  
  /**
   Constructor that creates histograms and prepares all internals
   \param _name Name of the monitor
   \param _nBins Number of histogram bins, used also as a variable step for efficiency
   \param min Minimum variable value
   \param max Maximum variable value
   */
  PerformanceMonitor(string _name, int _nBins, double min, double max);
  
  /// Assignment operator
  void operator=(const PerformanceMonitor &pm);
  
  /**
   Adds value of monitored parameter to the monitor for signal or background
   \param value Value of the parameter to be monitored
   \param signal Set to `true` for signal, `false` for background
   */
  void SetValue(double value, bool signal);
  
  /**
   Calculated values of all internal parameters: AUC, Z-score at different cut levels,
   max efficiency, min fake rate, max and average distance to √f curve.
   */
  void CalcEfficiency();
  
  /**
   Draws ROC curve in the current pad
   \param first Specify whether this is the first time a graph is drawn in this pad or not
  */
  void DrawRocGraph(bool first);
  
  /// Draws signal and background histograms in the current pad
  void DrawHists();
  
  /// Prints two columns: fake rate and efficiency, for different values of threshold
  void PrintFakesEfficiency();
  
  /// Prints values of internal parameters
  void PrintParams();
  
  /// Returns internal parameter given its name, or `-inf` if name is incorrect
  inline double GetValueByName(string name){
    if(name=="auc")           return auc;
    if(name=="sigma_init")    return significanceInitial;
    if(name=="sigma_L0")      return significanceAfterL0;
    if(name=="sigma_L1")      return significanceAfterL1;
    if(name=="max_eff")       return maxEfficiency;
    if(name=="min_fake")      return invFakeAtHighestEff;
    if(name=="max_dist_fake") return maxDistToSqrtFake;
    if(name=="avg_dist_fake") return avgDistToSqrtFake;
    return -inf;
  }
  
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
