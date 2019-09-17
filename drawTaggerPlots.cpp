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
  vector<string> monitorTypes = {
    "avg_hits",
    "max_hits",
    "avg_layers",
    "max_layers",
    "avg_length",
    "max_length",
    "n_helices"
  };
  
  Monitors monitors;
  
  for(auto monitorType : monitorTypes){
    int max = 20, nBins = 20;
    if(monitorType=="avg_length" || monitorType=="max_length"){ max = 10; nBins = 40; }
    monitors[monitorType] = PerformanceMonitor(monitorType, nBins, 0, max);
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
  TCanvas *canvasROC    = new TCanvas("canvasROC", "canvasROC", 1000, 1500);
  TCanvas *canvasDists  = new TCanvas("canvasDists", "canvasDists", 1000, 1500);
  canvasROC->Divide(2, 3);
  canvasDists->Divide(2, 3);
  
  int iPad=1;
  for(auto &[name, monitor] : monitors){
    monitor.CalcEfficiency();
    canvasROC->cd(iPad);
    monitor.DrawRocGraph(true);
    canvasDists->cd(iPad++);
    monitor.DrawHists();
  }
  
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
