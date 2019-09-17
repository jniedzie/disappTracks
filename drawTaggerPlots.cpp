//  drawTaggerPlots.cpp
//
//  Created by Jeremi Niedziela on 16/09/2019.

#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "EventProcessor.hpp"
#include "EventSet.hpp"
#include "HelixProcessor.hpp"
#include "PerformanceMonitor.hpp"

string configPath = "configs/helixTagger.md";
string cutLevel = "afterHelixTagging/";

xtracks::EDataType dataType = xtracks::kSignal;

vector<string> monitorTypes = {
  "avg_hits",
  "max_hits",
  "avg_layers",
  "max_layers",
  "avg_length",
  "max_length",
  "n_helices"
};

/**
 The program execution starting point.
 */
int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = ConfigManager(configPath);
  auto helixProcessor = HelixProcessor();
  
  EventSet events;
  events.LoadEventsFromFiles(cutLevel);
  
  TCanvas *canvas = new TCanvas("canvas", "canvas", 2880,1800);
  canvas->Divide(4,3);
  
  TH1D *maxNhitsSignal     = new TH1D("maxNhitsSignal", "maxNhitsSignal", 50, 0, 50);
  TH1D *maxNhitsBackground = new TH1D("maxNhitsBackground", "maxNhitsBackground", 50, 0, 50);
  
  
  map<string, PerformanceMonitor> monitors;
    
  for(auto monitorType : monitorTypes){
    int max = 20, nBins = 20;
    if(monitorType=="avg_length" || monitorType=="max_length"){ max = 10; nBins = 40; }
    monitors[monitorType] = PerformanceMonitor(monitorType, nBins, 0, max);
  }
  
  for(int iEvent=0; iEvent<events.size(dataType, kTaggerSignal); iEvent++){
    auto event = events.At(dataType, kTaggerSignal, iEvent);
    maxNhitsSignal->Fill(helixProcessor.GetMaxNhits(event->GetHelices()));
    
    for(string monitorType : monitorTypes){
      double value = helixProcessor.GetHelicesParamsByMonitorName(event->GetHelices(), monitorType);
      monitors[monitorType].SetValue(value, true);
    }
  }
  
  for(int iEvent=0; iEvent<events.size(dataType, kTaggerBackground); iEvent++){
    auto event = events.At(dataType, kTaggerBackground, iEvent);
    maxNhitsBackground->Fill(helixProcessor.GetMaxNhits(event->GetHelices()));
    for(string monitorType : monitorTypes){
      double value = helixProcessor.GetHelicesParamsByMonitorName(event->GetHelices(), monitorType);
      monitors[monitorType].SetValue(value, false);
    }
  }
  
  
  int iPad=1;
  for(auto &[name, monitor] : monitors){
    monitor.CalcEfficiency();
    canvas->cd(iPad++);
    monitor.DrawRocGraph(true);
    canvas->cd(iPad++);
    monitor.DrawHists();
  }
 
  canvas->Update();
  theApp.Run();
  return 0;
}
