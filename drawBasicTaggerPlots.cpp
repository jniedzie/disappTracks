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
string cutLevel = "afterHelixTagging";
string suffix = "_withEndcaps";





xtracks::EDataType dataType = xtracks::kSignal;

//  name          title  color pad thresholdUpBin
map<string, tuple<string, int, int, int>> monitorTypes = {
  {"avg_hits"   , {"Average number of hits on helix", 40         , 1 , 19 }},
  {"max_hits"   , {"Maximum number of hits on helix", 46         , 2 , 19 }},
  {"avg_layers" , {"Average number of helix layers" , kMagenta   , 3 , 7, }},
  {"max_layers" , {"Maximum number of helix layers" , kViolet+1  , 4 , 13 }},
  {"avg_length" , {"Average length of helix"        , kGreen+2   , 5 , 19 }},
  {"max_length" , {"Maximum length of helix"        , 49         , 6 , 34 }},
  {"n_helices"  , {"Number of helices per event"    , kCyan+1    , 7 , 9  }},
};

/// Defines types of monitors and initializes them
Monitors CreateMonitors(bool withPU)
{
  Monitors monitors;
  
  for(auto &[name, params] : monitorTypes){
    auto [title, color, pad, thresholdUp] = params;
    int max = 20, nBins = 20;
    if(name=="avg_length" || name=="max_length"){ max = 10; nBins = 40; }
    monitors[name] = PerformanceMonitor(name, title, nBins, 0, max, (EColor)color, withPU, thresholdUp);
  }
  return monitors;
}

/**
 Fills monitors with data found in events.
 \param isSignal If true, will use signal events, otherwise background
 \param withPU If true, will use event with pile-up
 */
void FillMonitors(Monitors &monitors, const EventSet &events, bool isSignal, bool withPU)
{
  ESignal dataSet;
  if(isSignal && withPU)        dataSet = kTaggerSignalWithPU;
  else if(isSignal && !withPU)  dataSet = kTaggerSignalNoPU;
  else if(!isSignal && withPU)  dataSet = kTaggerBackgroundWithPU;
  else if(!isSignal && !withPU) dataSet = kTaggerBackgroundNoPU;
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    for(int iEvent=0; iEvent<events.size(dataType, dataSet, year); iEvent++){
      auto event = events.At(dataType, dataSet, year, iEvent);
      if(!event->WasTagged()) continue;
      
      for(auto &[name, monitor] : monitors){
        double value = helixProcessor.GetHelicesParamsByMonitorName(event->GetHelices(), name);
        monitor.SetValue(value, isSignal);
      }
    }
  }
}

/// Calculates internal parameters of monitors and draws resulting plots
void DrawMonitors(Monitors &monitorsNoPU, Monitors *monitorsWithPU=nullptr)
{
  TCanvas *canvasROC          = new TCanvas("canvasROC", "canvasROC", 800, 600);
  TCanvas *canvasDists        = new TCanvas("canvasDists", "canvasDists", 1000, 1500);
  canvasDists->Divide(2, 4);
  
  TLegend *legendROC   = new TLegend(0.5, 0.1, 0.9, 0.4);
  TLegend *legendDists = new TLegend(0.1, 0.1, 0.4, 0.4);
  
  gStyle->SetOptStat(0);
  
  bool first = true;
//  for(auto &[name, monitor] : monitorsNoPU){
//    monitor.CalcEfficiency();
//    canvasROC->cd();
//    monitor.DrawRocGraph(first, legendROC);
//    int pad = get<2>(monitorTypes[name]);
//    canvasDists->cd(pad);
//    monitor.DrawHists(first, first ? legendDists : nullptr);
//    monitor.PrintFakesEfficiency();
//    if(first) first = false;
//  }
  if(monitorsWithPU){
    first = true;
    for(auto &[name, monitor] : *monitorsWithPU){
      monitor.CalcEfficiency();
      canvasROC->cd();
      monitor.DrawRocGraph(first, legendROC);
      int pad = get<2>(monitorTypes[name]);
      canvasDists->cd(pad);
      monitor.DrawHists(false, first ? legendDists : nullptr);
      
      
      double bestEff, bestFake;
      int thresholdLowBin, thresholdUpBin;
      
      double dist = monitor.GetMaxDistanceFromSqrtFake(bestEff, bestFake, thresholdLowBin, thresholdUpBin);
      cout<<name<<"\tdist: "<<dist<<"\teff: "<<bestEff<<"\tfake: "<<bestFake;
      cout<<"\tthreshold: "<<thresholdLowBin<<" -- "<<thresholdUpBin<<endl;
      
//      monitor.PrintFakesEfficiency();
      if(first) first = false;
    }
  }
  
  canvasROC->cd();
  legendROC->Draw("same");
  
  canvasDists->cd(8);
  legendDists->Draw("same");
  
  canvasROC->Update();
  canvasDists->Update();
  
  canvasROC->SaveAs(("plots/tagger_roc"+suffix+".pdf").c_str());
  canvasDists->SaveAs(("plots/tagger_distributions"+suffix+".pdf").c_str());
}

/// The program execution starting point.
int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = ConfigManager(configPath);
  auto helixProcessor = HelixProcessor();
  
  EventSet events;
  events.LoadEventsFromFiles(cutLevel+suffix+"/");
  
  Monitors monitorsNoPU   = CreateMonitors(false);
  Monitors monitorsWithPU = CreateMonitors(true);
  
  FillMonitors(monitorsNoPU, events, true,  false); // signal, noPU
  FillMonitors(monitorsNoPU, events, false, false); // bkg, noPU
  FillMonitors(monitorsWithPU, events, true , true);  // signal, withPU
  FillMonitors(monitorsWithPU, events, false, true);  // bkg, withPU
  
  DrawMonitors(monitorsNoPU, &monitorsWithPU);
  
  theApp.Run();
  return 0;
}
