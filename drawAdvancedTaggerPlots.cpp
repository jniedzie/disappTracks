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
string cutLevel = "afterHelixTagging";
//string suffix = "";
string suffix = "_paramSet3";
//string suffix = "_noNoise";

xtracks::EDataType dataType = xtracks::kSignal;
ESignal signalDataset = kTaggerSignalNoPU;
ESignal backgroundDataset = kTaggerSignalNoPUpionRemoved;

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

map<ETestParams, string> paramTitles = {
  {kNoBins          , "No binning"        },
  {kPionPt          , "Pion p_{T} (MeV)"  },
  {kCharginoEta     , "Track |#eta|"      },
  {kCharginoNlayers , "Track N layers"    },
  {kMET             , "MET p_{T} (GeV)"   },
  {kCharginoCharge  , "Track charge"      },
  {kCharginoPt      , "Track p_{T} (GeV)" },
};

map<ETestParams, vector<range<double>>> paramRanges = {
  { kNoBins, {range<double>(-1, 1)} },
  { kPionPt,
    {
      range<double>(0, 150),
      range<double>(150, 300),
      range<double>(300, 450),
      range<double>(450, 600),
      range<double>(600, 750) }
  },
  { kCharginoEta,
    { range<double>(0, 0.5),
      range<double>(0.5, 1.0),
      range<double>(1.0, 1.5),
      range<double>(1.5, 2.0),
      range<double>(2.0, 2.5) }
  },
  { kCharginoNlayers,
    { range<double>(3, 3),
      range<double>(4, 4),
      range<double>(5, 5),
      range<double>(6, 6),
      range<double>(7, 7) }
  },
  { kMET,
    {
      range<double>(0  , 200),
      range<double>(200, 400),
      range<double>(400, 600),
      range<double>(600, 800),
      range<double>(800,1000) }
  },
  { kCharginoCharge,
    { range<double>(-2, 0),
      range<double>(0, 2) }
  },
  { kCharginoPt,
    {
      range<double>(0, 200),
      range<double>(200, 400),
      range<double>(400, 600),
      range<double>(600, 800),
      range<double>(800, 1000),
      range<double>(1000,1200)
    }
  },
};

void GetBinsForRanges(float *binsArray, int &arraySize, ETestParams param)
{
  vector<range<double>> ranges = paramRanges[param];
  vector<double> bins;
  bins.push_back(ranges[0].GetMin());
  
  for(range<double> r : ranges) bins.push_back(r.GetMax());
  
  for(int iBin=0; iBin<bins.size(); iBin++) binsArray[iBin] = bins[iBin];
  arraySize = (int)bins.size()-1;
}

bool IsEventOk(const Event &event, ETestParams param)
{
  if(param == kPionPt){
    return event.GetGenPionHelices().size() == 1;
  }
  
  if(param == kCharginoEta || param == kCharginoNlayers ||
     param == kCharginoCharge || param == kCharginoPt){
    return event.GetNtracks() == 1;
  }
  
  return true;
}

double EventToParam(const Event &event, ETestParams param)
{
  if(param == kPionPt)          return event.GetGenPionHelices().front().GetMomentum().GetTransverse();
  if(param == kCharginoEta)     return fabs(event.GetTrack(0)->GetEta());
  if(param == kCharginoNlayers) return event.GetTrack(0)->GetNtrackerLayers();
  if(param == kMET)             return event.GetMetPt();
  if(param == kCharginoCharge)  return event.GetTrack(0)->GetCharge();
  if(param == kCharginoPt)      return event.GetTrack(0)->GetPt();
  return 0;
}

//  name          title  color pad thresholdUpBin
map<string, tuple<string, int, int, int>> monitorTypes = {
  {"avg_hits"   , {"Average number of hits on helix", 40         , 1 , 19 }},
  {"max_hits"   , {"Maximum number of hits on helix", 46         , 2 , 19 }},
  {"avg_layers" , {"Average number of helix layers" , kMagenta   , 3 , 39, }},
  {"max_layers" , {"Maximum number of helix layers" , kViolet+1  , 4 , 19 }},
  {"avg_length" , {"Average length of helix"        , kGreen+2   , 5 , 19 }},
  {"max_length" , {"Maximum length of helix"        , 49         , 6 , 39 }},
  {"n_helices"  , {"Number of helices per event"    , kCyan+1    , 7 , 9  }},
};

/// Defines types of monitors and initializes them
vector<Monitors> CreateMonitors(ETestParams param)
{
  vector<Monitors> monitors;
  
  for(range<double> r : paramRanges[param]){ // This gives a range of variable, e.g. eta
    Monitors monitorsForParamValue;
    for(auto &[name, params] : monitorTypes){
      auto [title, color, pad, thresholdUp] = params;
      int max = 20, nBins = 20;
      if(name=="avg_length" || name=="max_length"){ max = 10; nBins = 40; }
      if(name=="n_helices"){ max = 1000; nBins = 10; }
      monitorsForParamValue[name] = PerformanceMonitor(name, title, nBins, 0, max, (EColor)color, true, thresholdUp);
    }
    monitors.push_back(monitorsForParamValue);
  }
  return monitors;
}

/**
Fills monitors with data found in events.
\param isSignal If true, will use signal events, otherwise background
*/
void FillMonitors(vector<Monitors> &monitors, const EventSet &events, bool isSignal, ETestParams param)
{
  ESignal dataSet = isSignal ? signalDataset : backgroundDataset;
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
   
    for(int iEvent=0; iEvent<events.size(dataType, dataSet, year); iEvent++){
      auto event = events.At(dataType, dataSet, year, iEvent);
      if(!event->WasTagged()) continue;
      
      if(isSignal){
        if(!IsEventOk(*event, param)) continue;
        
        double value = EventToParam(*event, param);
        int monitorBin = -1;
        
        for(int bin=0; bin<paramRanges[param].size(); bin++){
          if(paramRanges[param][bin].IsInside(value)){
            monitorBin = bin;
            break;
          }
        }
        
        if(monitorBin < 0) continue;
        
        for(auto &[name, monitor] : monitors[monitorBin]){
          double value = helixProcessor.GetHelicesParamsByMonitorName(event->GetHelices(), name);
          monitor.SetValue(value, isSignal);
        }
      }
      else{
        for(int bin=0; bin<paramRanges[param].size(); bin++){
          for(auto &[name, monitor] : monitors[bin]){
            double value = helixProcessor.GetHelicesParamsByMonitorName(event->GetHelices(), name);
            monitor.SetValue(value, isSignal);
          }
        }
      }
        
      
    }
  }
}

/// Calculates internal parameters of monitors and draws resulting plots
void DrawMonitors(vector<Monitors> &monitors, ETestParams param, bool first)
{
  TLegend *legend = new TLegend(0.5, 0.1, 0.9, 0.4);
  
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0);
  
  float histBins[100];
  int nBins;
  
  GetBinsForRanges(histBins, nBins, param);
  
  map<string, TH1D*> maxDistance;
  
  int nameIter=0;
  
  double bestEff = -inf;
  double bestFake = inf;
  double bestDist = -inf;
  
  for(auto &[name, params] : monitorTypes){
    auto [title, color, pad, thresholdUp] = params;
    
    maxDistance[name] = new TH1D(to_string(RandInt(0, inf)).c_str(), "", nBins, histBins);
    maxDistance[name]->GetXaxis()->SetTitle(paramTitles[param].c_str());
    maxDistance[name]->GetYaxis()->SetTitle("Max distance from #sqrt{c_{fake}}");
    maxDistance[name]->SetBarWidth(0.1);
    maxDistance[name]->SetBarOffset(0.1 + nameIter * 0.1);
    maxDistance[name]->SetFillColor(color);
    maxDistance[name]->SetStats(0);
    
    for(int paramBin=0; paramBin<paramRanges[param].size(); paramBin++){
      
      PerformanceMonitor monitor = monitors[paramBin][name];
      monitor.CalcEfficiency();
      double eff, fake;
      int thresholdLowBin, thresholdUpBin;
      double distance = monitor.GetMaxDistanceFromSqrtFake(eff, fake, thresholdLowBin, thresholdUpBin);
      maxDistance[name]->SetBinContent(paramBin+1, distance);
      if(distance > bestDist){
        bestDist = distance;
        bestEff = eff;
        bestFake = fake;
      }
    }
    maxDistance[name]->Draw(nameIter == 0 ? "b" : "b same");
    legend->AddEntry(maxDistance[name], title.c_str(), "f");
    
    nameIter++;
  }
  maxDistance["avg_hits"]->SetMaximum(1.3);
  maxDistance["avg_hits"]->SetMinimum(-1.3);
  
  cout<<"Best distance: "<<bestDist<<"\teff: "<<bestEff<<"\tfake: "<<bestFake<<endl;
  
  if(first) legend->Draw("same");
}

/// The program execution starting point.
int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = ConfigManager(configPath);
  auto helixProcessor = HelixProcessor();
  
  EventSet events; events.LoadEventsFromFiles(getPathPrefix());
  
  TCanvas *canvasPerformance  = new TCanvas("canvasPerformance", "canvasPerformance", 1000, 1500);
  canvasPerformance->Divide(2, 3);
  
  for(int testParam=1; testParam<nTestParams; testParam++){
    canvasPerformance->cd(testParam);
    
    vector<Monitors> monitors = CreateMonitors((ETestParams)testParam);
    
    FillMonitors(monitors, events, true, (ETestParams)testParam); // signal
    FillMonitors(monitors, events, false, (ETestParams)testParam); // bkg
    
    DrawMonitors(monitors, (ETestParams)testParam, testParam==1);
  }
    
  canvasPerformance->Update();
  canvasPerformance->SaveAs(("plots/tagger_performance"+suffix+".pdf").c_str());
  
  theApp.Run();
  return 0;
}
