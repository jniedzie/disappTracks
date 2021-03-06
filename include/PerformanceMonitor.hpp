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
   \param color Optionally, color of the ROC graph points and curve
   \param alternativeColors When switched on, alternative colors will be used for signal and background histograms
   */
  PerformanceMonitor(string _name, string _title, int _nBins, double min, double max, EColor color = kRed,
                     bool alternativeColors=false, int _fixUpperThresholdBin=-1);
  
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
   \param legend If legend is provided, graph will be added to it
  */
  void DrawRocGraph(bool first, TLegend *legend = nullptr);
  
  /**
   Draws signal and background histograms in the current pad
   \param first Specify whether this is the first time hists are drawn in this pad or not
   \param legend If legend is provided, hists will be added to it
   */
  void DrawHists(bool first, TLegend *legend = nullptr);
  
  /// Prints two columns: fake rate and efficiency, for different values of threshold
  void PrintFakesEfficiency();
  
  /// Prints values of internal parameters
  void PrintParams();
  
  /// Returns internal parameter given its name
  inline double GetValueByName(string name){ return params[name];}
  
  /// Returns the highest distance between efficiency and c_eff = sqrt(c_fake) curve
  double GetMaxDistanceFromSqrtFake(double &bestEff, double &bestFake,
                                    int &thresholdLowBin, int &thresholdUpBin);
  
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
  TGraphErrors* rocGraph;
  
  string name, title;
  
  vector<vector<double>> efficiency;  ///< efficiency for given lower and upper threshold
  vector<vector<double>> fakeRate;    ///< fake rate for given lower and upper threshold
  
  map<string, double> params;
  
  TF1* GetRocFunction();
  void CountEventsAboveThreshold();
  void ResetParams();
  
  int fixedUpperThresholdBin;
  
};

#endif /* PerformanceMonitor_hpp */
