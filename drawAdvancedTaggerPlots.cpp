//  drawAdvancedTaggerPlots.cpp
//
//  Created by Jeremi Niedziela on 18/09/2019.

#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "EventProcessor.hpp"
#include "EventSet.hpp"
#include "HelixProcessor.hpp"
#include "PerformanceMonitor.hpp"
#include "Logger.hpp"

string configPath = "configs/taggerPlotting.md";
string cutLevel = "afterHelixTagging/";

xtracks::EDataType dataType = xtracks::kSignal;

enum ETestParams {
  kNoBins,
  kPionPt,
  kCharginoEta,
  kCharginoNlayers,
  kMET,
  kCharginoCharge,
  kCharginoPt,
  nTestParams
};

ETestParams testParam = kCharginoEta;

map<ETestParams, vector<range<double>>> paramRanges = {
  { kNoBins, {range<double>()} },
  { kPionPt,
    {
      range<double>(0, 150),
      range<double>(150, 300),
      range<double>(300, 450),
      range<double>(450, 600),
      range<double>(600, inf) }
  },
  { kCharginoEta,
    { range<double>(0, 0.5),
      range<double>(0.5, 1.0),
      range<double>(1.0, 1.5),
      range<double>(1.5, 3.0) }
  },
  { kCharginoNlayers,
    { range<double>(3, 3),
      range<double>(4, 4),
      range<double>(5, 5),
      range<double>(6, 6),
      range<double>(7, inf) }
  },
  { kMET,
    {
      range<double>(200, 500),
      range<double>(500, 600),
      range<double>(600, 700),
      range<double>(700, inf) }
  },
  { kCharginoCharge,
    { range<double>(-inf, 0),
      range<double>(0, inf) }
  },
  { kCharginoPt,
    {
      range<double>(0, 200),
      range<double>(200, 400),
      range<double>(400, 600),
      range<double>(600, 800),
      range<double>(800, 1000),
      range<double>(1000, inf)
    }
  },
};

void GetBinsForRanges(float *binsArray, int &arraySize)
{
  vector<range<double>> ranges = paramRanges[testParam];
  vector<double> bins;
  bins.push_back(ranges[0].GetMin());
  
  for(range<double> r : ranges) bins.push_back(r.GetMax());
  
  for(int iBin=0; iBin<bins.size(); iBin++) binsArray[iBin] = bins[iBin];
  arraySize = (int)bins.size()-1;
}

bool IsEventOk(const Event &event)
{
  if(testParam == kPionPt){
    return event.GetGenPionHelices().size() == 1;
  }
  
  if(testParam == kCharginoEta || testParam == kCharginoNlayers ||
     testParam == kCharginoCharge || testParam == kCharginoPt){
    return event.GetNtracks() == 1;
  }
  
  return true;
}

double EventToParam(const Event &event)
{
  if(testParam == kPionPt)          return event.GetGenPionHelices().front().GetMomentum().GetTransverse();
  if(testParam == kCharginoEta)     return fabs(event.GetTrack(0)->GetEta());
  if(testParam == kCharginoNlayers) return event.GetTrack(0)->GetNtrackerLayers();
  if(testParam == kMET)             return event.GetMetPt();
  if(testParam == kCharginoCharge)  return event.GetTrack(0)->GetCharge();
  if(testParam == kCharginoPt)      return event.GetTrack(0)->GetPt();
  return 0;
}

//  name          title  color pad
map<string, tuple<string, int, int>> monitorTypes = {
  {"avg_hits"   , {"Average number of hits on helix", 40         , 1 }},
  {"max_hits"   , {"Maximum number of hits on helix", 46         , 2 }},
  {"avg_layers" , {"Average number of helix layers" , kMagenta   , 3 }},
  {"max_layers" , {"Maximum number of helix layers" , kViolet+1  , 4 }},
  {"avg_length" , {"Average length of helix"        , kGreen+2   , 5 }},
  {"max_length" , {"Maximum length of helix"        , 49         , 6 }},
  {"n_helices"  , {"Number of helices per event"    , kCyan+1    , 7 }},
};

/// Defines types of monitors and initializes them
vector<Monitors> CreateMonitors()
{
  vector<Monitors> monitors;
  
  for(range<double> r : paramRanges[testParam]){ // This gives a range of variable, e.g. eta
    Monitors monitorsForParamValue;
    for(auto &[name, params] : monitorTypes){
      auto [title, color, pad] = params;
      int max = 20, nBins = 20;
      if(name=="avg_length" || name=="max_length"){ max = 10; nBins = 40; }
      monitorsForParamValue[name] = PerformanceMonitor(name, title, nBins, 0, max, (EColor)color);
    }
    monitors.push_back(monitorsForParamValue);
  }
  return monitors;
}

/// Fills `monitors` with data found in `events`. Will use signal or background events, depending on `isSignal` value
void FillMonitors(vector<Monitors> &monitors, const EventSet &events, bool isSignal)
{
  for(int iEvent=0;
          iEvent<events.size(dataType, isSignal ? kTaggerSignal : kTaggerBackground);
          iEvent++){
    auto event = events.At(dataType, isSignal ? kTaggerSignal : kTaggerBackground, iEvent);
    if(!IsEventOk(*event)) continue;
    
    double value = -inf;
    
    if(testParam == kCharginoEta) value = fabs(event->GetTrack(0)->GetEta());
    
    int monitorBin = -1;
    
    for(int bin=0; bin<paramRanges[testParam].size(); bin++){
      if(paramRanges[testParam][bin].IsInside(value)){
        monitorBin = bin;
        break;
      }
    }
    
    if(monitorBin < 0){
      Log(0)<<"Could not find bin for parameter value!!\n";
      continue;
    }
    
    for(auto &[name, monitor] : monitors[monitorBin]){
      double value = helixProcessor.GetHelicesParamsByMonitorName(event->GetHelices(), name);
      monitor.SetValue(value, isSignal);
    }
  }
}

/// Calculates internal parameters of monitors and draws resulting plots
void DrawMonitors(vector<Monitors> &monitors)
{
  TCanvas *canvasPerformance  = new TCanvas("canvasPerformance", "canvasPerformance", 800, 600);
  canvasPerformance->Divide(1, 1);
  
  TLegend *legend = new TLegend(0.5, 0.6, 0.9, 0.9);
  
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  
  float histBins[100];
  int nBins;
  
  GetBinsForRanges(histBins, nBins);
  
  map<string, TH1D*> maxDistance;
  
  
  canvasPerformance->cd();
  
  int nameIter=0;
  double max=-inf, min=inf;
  
  for(auto &[name, params] : monitorTypes){
    auto [title, color, pad] = params;
    
    maxDistance[name] = new TH1D("maxDistance", "maxDistance", nBins, histBins);
    maxDistance[name]->SetBarWidth(0.1);
    maxDistance[name]->SetBarOffset(0.1 + nameIter * 0.1);
    maxDistance[name]->SetFillColor(color);
    maxDistance[name]->SetStats(0);
    
    for(int paramBin=0; paramBin<paramRanges[testParam].size(); paramBin++){
      
      PerformanceMonitor monitor = monitors[paramBin][name];
      monitor.CalcEfficiency();
      double distance = monitor.GetMaxDistanceFromSqrtFake();
      
      maxDistance[name]->SetBinContent(paramBin, distance);
      if(distance < min) min = distance;
      if(distance > max) max = distance;
    }
    maxDistance[name]->Draw(nameIter == 0 ? "b" : "b same");
    legend->AddEntry(maxDistance[name], title.c_str(), "f");
    
    nameIter++;
  }
//  maxDistance["avg_hits"]->SetMaximum(1.1*max);
//  maxDistance["avg_hits"]->SetMinimum(1.1*min);
  maxDistance["avg_hits"]->SetMaximum(2);
  maxDistance["avg_hits"]->SetMinimum(-2);
  
  
  canvasPerformance->cd();
  legend->Draw("same");
  
  canvasPerformance->Update();
  canvasPerformance->SaveAs("plots/tagger_performance.pdf");
}

/// The program execution starting point.
int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = ConfigManager(configPath);
  auto helixProcessor = HelixProcessor();
  
  EventSet events;
  events.LoadEventsFromFiles(cutLevel);
  
  vector<Monitors> monitors = CreateMonitors();
  
  FillMonitors(monitors, events, true);
  FillMonitors(monitors, events, false);
  
  DrawMonitors(monitors);
  
  theApp.Run();
  return 0;
}
