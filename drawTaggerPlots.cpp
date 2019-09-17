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

/// Defines types of monitors and initializes them
Monitors CreateMonitors()
{
  vector<tuple<string, string, int>> monitorTypes = {
    {"avg_hits"   , "Average number of hits on helix", 40 },
    {"max_hits"   , "Maximum number of hits on helix", 46 },
    {"avg_layers" , "Average number of helix layers" , 30 },
    {"max_layers" , "Maximum number of helix layers" , 33 },
    {"avg_length" , "Average length of helix"        , 9  },
    {"max_length" , "Maximum length of helix"        , 49 },
    {"n_helices"  , "Number of helices per event"    , 29 },
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
    auto event = events.At(dataType, kTaggerSignal, iEvent);
    
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
