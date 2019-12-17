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
<<<<<<< HEAD
//string suffix = "";
=======
string suffix = "";
>>>>>>> Updating event display
//string suffix = "_paramSet3";
//string suffix = "_noHighPtHits";
//string suffix = "_default";
//string suffix = "_highMerging";
//string suffix = "_lowSeedChi";
//string suffix = "_noMissing";
string suffix = "_lowTrackChi";
//string suffix = "_removingPionHits";

xtracks::EDataType dataType = xtracks::kSignal;
ESignal signalDataset     = kChargino500_10;
ESignal backgroundDataset = kTaggerBackgroundWithPU;


/// Returns path prefix for cuts level and category selected in the config file
string getPathPrefix()
{
  string prefix = "";
   
  if(config.secondaryCategory == "Zmumu") prefix += "Zmumu/";
  if(config.secondaryCategory == "Wmunu") prefix += "Wmunu/";
  
  if(config.params["cuts_level"]==0) prefix += "after_L0/";
  if(config.params["cuts_level"]==1) prefix += "after_L1/"+config.category+"/";
  
  prefix += "afterHelixTagging"+suffix+"/";
  
  return prefix;
}

map<string, tuple<string, int, int, int>> monitorTypes = {
//  name          title                              color       pad  thresholdUpBin
  {"avg_hits"   , {"Average number of hits on helix", 40         , 1 , 19 }},
  {"max_hits"   , {"Maximum number of hits on helix", 46         , 2 , 19 }},
  {"avg_layers" , {"Average number of helix layers" , kMagenta   , 3 , 19 }},
  {"max_layers" , {"Maximum number of helix layers" , kViolet+1  , 4 , 19 }},
  {"avg_length" , {"Average length of helix"        , kGreen+2   , 5 , 39 }},
  {"max_length" , {"Maximum length of helix"        , 49         , 6 , 39 }},
  {"n_helices"  , {"Number of helices per event"    , kCyan+1    , 7 , 9  }},
};

/// Defines types of monitors and initializes them
Monitors CreateMonitors()
{
  Monitors monitors;
  
  for(auto &[name, params] : monitorTypes){
    auto [title, color, pad, thresholdUp] = params;
    int max = 20, nBins = 20;
    if(name=="avg_length" || name=="max_length"){ max = 10; nBins = 40; }
    if(name=="n_helices"){ max = 1000; nBins = 10; }
    monitors[name] = PerformanceMonitor(name, title, nBins, 0, max, (EColor)color, true, thresholdUp);
  }
  return monitors;
}

/**
 Fills monitors with data found in events.
 \param isSignal If true, will use signal events, otherwise background
 \param withPU If true, will use event with pile-up
 */
void FillMonitors(Monitors &monitors, const EventSet &events, bool isSignal)
{
  ESignal dataSet = isSignal ? signalDataset : backgroundDataset;
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    for(int iEvent=0; iEvent<events.size(dataType, dataSet, year); iEvent++){
      auto event = events.At(dataType, dataSet, year, iEvent);
      
      if(isSignal){
        bool hasHighMomnetumPion = false;
        
        for(Helix pion : event->GetGenPionHelices()){
          double pt = sqrt(pow(pion.GetMomentum().GetX(), 2) + pow(pion.GetMomentum().GetY(), 2));
          if(pt > 400){
            hasHighMomnetumPion = true;
            break;
          }
        }
        
        if(!hasHighMomnetumPion) continue;
      }
        
      bool hasLowMomentumLepton = false;
      
      for(int iLepton=0; iLepton<event->GetNleptons(); iLepton++){
        if(event->GetLepton(iLepton)->GetPt() < 1.0){
          hasLowMomentumLepton = true;
          break;
        }
      }
      
      if(hasLowMomentumLepton) continue;
      
      if(!event->WasTagged()) continue;
      
      for(auto &[name, monitor] : monitors){
        double value = helixProcessor.GetHelicesParamsByMonitorName(event->GetHelices(), name);
        monitor.SetValue(value, isSignal);
      }
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
    monitor.DrawHists(false, first ? legendDists : nullptr);
    
    
    double bestEff, bestFake;
    int thresholdLowBin, thresholdUpBin;
    
    double dist = monitor.GetMaxDistanceFromSqrtFake(bestEff, bestFake, thresholdLowBin, thresholdUpBin);
    cout<<name<<"\tdist: "<<dist<<"\teff: "<<bestEff<<"\tfake: "<<bestFake;
    cout<<"\tthreshold: "<<thresholdLowBin<<" -- "<<thresholdUpBin<<endl;
    
    //      monitor.PrintFakesEfficiency();
    if(first) first = false;
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
  
  EventSet events; events.LoadEventsFromFiles(getPathPrefix());
  
  Monitors monitors = CreateMonitors();
  
  FillMonitors(monitors, events, true);
  FillMonitors(monitors, events, false);
  
  DrawMonitors(monitors);
  
  theApp.Run();
  return 0;
}
