//  drawBasicTaggerPlots.cpp
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
Monitors CreateMonitors()
{
  Monitors monitors;
  
  for(auto &[name, params] : monitorTypes){
    auto [title, color, pad] = params;
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
  TCanvas *canvasROC          = new TCanvas("canvasROC", "canvasROC", 800, 600);
  TCanvas *canvasDists        = new TCanvas("canvasDists", "canvasDists", 1000, 1500);
  canvasDists->Divide(2, 4);
  
  TLegend *legendROC   = new TLegend(0.5, 0.1, 0.9, 0.4);
  TLegend *legendDists = new TLegend(0.1, 0.1, 0.4, 0.4);
  
  gStyle->SetOptStat(0);
  
  bool first = true;
  for(auto &[name, monitor] : monitors){
    monitor.CalcEfficiency();
    canvasROC->cd();
    monitor.DrawRocGraph(first, legendROC);
    int pad = get<2>(monitorTypes[name]);
    canvasDists->cd(pad);
    monitor.DrawHists(first ? legendDists : nullptr);
    monitor.PrintFakesEfficiency();
    if(first) first = false;
  }
  
  canvasROC->cd();
  legendROC->Draw("same");
  
  canvasDists->cd(8);
  legendDists->Draw("same");
  
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
  
  Monitors monitors = CreateMonitors();
  
  FillMonitors(monitors, events, true);
  FillMonitors(monitors, events, false);
  
  DrawMonitors(monitors);
  
  theApp.Run();
  return 0;
}
