//  drawTaggerPlots.cpp
//
//  Created by Jeremi Niedziela on 16/09/2019.

#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "EventProcessor.hpp"
#include "EventSet.hpp"
#include "HelixProcessor.hpp"
#include "PerformanceMonitor.hpp"

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

ETestParams testParam = kNoBins;

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
      range<double>(1.5, inf) }
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

/// Defines types of monitors and initializes them
Monitors CreateMonitors()
{
  vector<tuple<string, string, int>> monitorTypes = {
    {"avg_hits"   , "Average number of hits on helix", 40 },
    {"max_hits"   , "Maximum number of hits on helix", 46 },
    {"avg_layers" , "Average number of helix layers" , 30 },
    {"max_layers" , "Maximum number of helix layers" , 33 },
    {"avg_length" , "Average length of helix"        , kGreen+2  },
    {"max_length" , "Maximum length of helix"        , 49 },
    {"n_helices"  , "Number of helices per event"    , kRed+1 },
  };
  
  Monitors monitors;
  
  for(auto &[name, title, color] : monitorTypes){
    int max = 20, nBins = 20;
    if(name=="avg_length" || name=="max_length"){ max = 10; nBins = 40; }
    monitors[name] = PerformanceMonitor(name, title, nBins, 0, max, (EColor)color);
  }
  return monitors;
}

/// Fills `monitors` with data found in `events`. Will use signal or background events, depending on `isSignal` value
void FillMonitors(Monitors &monitors, const EventSet &events, bool isSignal)
{
  for(int iEvent=0;
          iEvent<events.size(dataType, isSignal ? kTaggerSignal : kTaggerBackground);
          iEvent++){
    auto event = events.At(dataType, isSignal ? kTaggerSignal : kTaggerBackground, iEvent);
    
    for(auto &[name, monitor] : monitors){
      double value = helixProcessor.GetHelicesParamsByMonitorName(event->GetHelices(), name);
      monitor.SetValue(value, isSignal);
    }
  }
}

/// Calculates internal parameters of monitors and draws resulting plots
void DrawMonitors(Monitors &monitors)
{
  TCanvas *canvasROC    = new TCanvas("canvasROC", "canvasROC", 800, 600);
  TCanvas *canvasDists  = new TCanvas("canvasDists", "canvasDists", 1000, 1500);
  canvasDists->Divide(2, 4);
  
  TLegend *legend = new TLegend(0.5, 0.1, 0.9, 0.4);
  
  int iPad=1;
  bool first = true;
  for(auto &[name, monitor] : monitors){
    monitor.CalcEfficiency();
    canvasROC->cd();
    monitor.DrawRocGraph(first, legend);
    if(first) first = false;
    canvasDists->cd(iPad++);
    monitor.DrawHists();
    
    
    monitor.PrintFakesEfficiency();
  }
  
  canvasROC->cd();
  legend->Draw("same");
  
  canvasROC->Update();
  canvasDists->Update();
  
  canvasROC->SaveAs("plots/tagger_roc.pdf");
  canvasDists->SaveAs("plots/tagger_distributions.pdf");
}

/// The program execution starting point.
int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = ConfigManager(configPath);
  auto helixProcessor = HelixProcessor();
  
  EventSet events;
  events.LoadEventsFromFiles(cutLevel);
  
  map<string, PerformanceMonitor> monitors = CreateMonitors();
  
  FillMonitors(monitors, events, true);
  FillMonitors(monitors, events, false);
  
  DrawMonitors(monitors);
  
  theApp.Run();
  return 0;
}
