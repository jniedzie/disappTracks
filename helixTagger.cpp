//  helixTagger.cpp
//
//  Created by Jeremi Niedziela on 25/01/2019.

#include "Helpers.hpp"
#include "Fitter.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "PerformanceMonitor.hpp"
#include "EventSet.hpp"

string configPath = "configs/helixTagger.md";
string cutLevel = "after_L1/all/";

bool printROCpoints = false;

enum ETestParams {
  kPionPt,
  kCharginoEta,
  kCharginoNlayers,
  kMET,
  kCharginoCharge,
  kCharginoPt,
  nTestParams
};

ETestParams testParam = kPionPt;

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;

vector<string> monitorTypes = {
  "avg_hits",
  "max_hits",
  "avg_layers",
  "max_layers",
  "avg_length",
  "max_length",
  "n_helices"
};

vector<shared_ptr<Point>> GetClustersNoEndcaps(const shared_ptr<Event> &event, bool removePionClusters);

map<ETestParams, vector<range<double>>> paramRanges = {
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
    { range<double>(200, 400),
      range<double>(400, 600),
      range<double>(600, inf) }
  },
  { kCharginoCharge,
    { range<double>(-inf, 0),
      range<double>(0, inf) }
  },
  { kCharginoPt,
    { range<double>(0, 200),
      range<double>(200, 400),
      range<double>(400, 600),
      range<double>(600, inf) }
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
  if(testParam == kPionPt)          return event.GetGenPionHelices().front().GetMomentum()->GetTransverse();
  if(testParam == kCharginoEta)     return fabs(event.GetTrack(0)->GetEta());
  if(testParam == kCharginoNlayers) return event.GetTrack(0)->GetNtrackerLayers();
  if(testParam == kMET)             return event.GetMetPt();
  if(testParam == kCharginoCharge)  return event.GetTrack(0)->GetCharge();
  if(testParam == kCharginoPt)      return event.GetTrack(0)->GetPt();
  return 0;
}

int main(int argc, char* argv[])
{
  cout.imbue(locale("de_DE"));
  // TApplication theApp("App", &argc, argv);
  
  config = ConfigManager(configPath);
  auto fitter = make_unique<Fitter>();
  TCanvas *canvas = new TCanvas("ROC", "ROC", 2880,1800);
  canvas->Divide(3,3);
  
  vector<map<string, PerformanceMonitor>> monitors; // [paramBin][monitorType]
  
  for(int iParam=0; iParam<paramRanges[testParam].size(); iParam++){
    map<string, PerformanceMonitor> monitorsForParam;
    
    for(auto monitorType : monitorTypes){
      int max = 20, nBins = 20;
      if(monitorType=="avg_length" || monitorType=="max_length"){ max = 10; nBins = 40; }
      monitorsForParam[monitorType] = PerformanceMonitor(monitorType, nBins, 0, max);
    }
    monitors.push_back(monitorsForParam);
  }
  
  EventSet events;
  events.LoadEventsFromFiles(cutLevel);
  int nEvents=events.size(dataType, setIter);
  
  auto start = now();
  
  int nAnalyzedEvents=0;
  
  for(auto iEvent=0; iEvent<nEvents; iEvent++){
    cout<<"\n\n=================================================================\n"<<endl;
    cout<<"helixTagger -- processing event "<<iEvent<<endl;
    
    auto event = events.At(dataType, setIter, iEvent);
    
    if(!IsEventOk(*event)) continue;
    double paramValue = EventToParam(*event);
    
    int iParam = -1;
    
    for(iParam=0; iParam<paramRanges[testParam].size(); iParam++){
      if(paramRanges[testParam][iParam].IsInside(paramValue)) break;
    }
    if(iParam < 0 || iParam == paramRanges[testParam].size()) continue;
    
    
    for(auto &track : event->GetTracks()){
      
      auto pointsNoEndcapsSignal = GetClustersNoEndcaps(event, false);
      auto pointsNoEndcapsBackground = GetClustersNoEndcaps(event, true);
      
      vector<Helix> fittedHelicesSignal     = fitter->FitHelices(pointsNoEndcapsSignal, *track, *event->GetVertex());
      vector<Helix> fittedHelicesBackground = fitter->FitHelices(pointsNoEndcapsBackground, *track, *event->GetVertex());
      
      // for(auto helix : fittedHelicesSignal) event->AddHelix(move(fittedHelix));
      
      for(string monitorType : monitorTypes){
        monitors[iParam][monitorType].SetValues(helixProcessor.GetHelicesParamsByMonitorName(fittedHelicesSignal, monitorType),
                                                helixProcessor.GetHelicesParamsByMonitorName(fittedHelicesBackground, monitorType));
        
      }
    }
    nAnalyzedEvents++;
  }
  
  int iPad=1;
  for(auto &monitorsForParam : monitors){
    for(auto &monitor : monitorsForParam){
      monitor.second.CalcEfficiency();
      canvas->cd(iPad++);
      monitor.second.DrawRocGraph(true);
      canvas->cd(iPad++);
      monitor.second.DrawHists();
    }
  }
  
  int iParam=0;
  for(auto &monitorsForParam : monitors){
    cout<<"\n\n============================================================"<<endl;
    cout<<"Param bin: "<<iParam++<<endl;
    
    for(auto &monitor : monitorsForParam){
      cout<<"\nMonitor: "<<monitor.first<<endl;
      if(printROCpoints) monitor.second.PrintFakesEfficiency();
      monitor.second.PrintParams();
    }
  }
  
  cout<<"N events analyzed: "<<nAnalyzedEvents<<endl;
  
  //  cout<<"helixTagger -- saving events"<<endl;
  //  events.SaveEventsToFiles("afterHelixTagging/");
  //  cout<<"helixTagger -- finished"<<endl;
  
  cout<<"Time: "<<duration(start, now());
  
  canvas->Update();
  
  TFile *outFile = new TFile("helixTagger.root", "recreate");
  outFile->cd();
  canvas->Write();
  outFile->Close();
  
  //  theApp.Run();
  return 0;
}

vector<shared_ptr<Point>> GetClustersNoEndcaps(const shared_ptr<Event> &event, bool removePionClusters)
{
  vector<shared_ptr<Point>> pointsNoEndcaps;
  
  for(auto &point : event->GetTrackerClusters()){
    if(point->GetSubDetName() == "TID" || point->GetSubDetName() == "TEC" || point->GetSubDetName() == "P1PXEC") continue;
    
    if(point->GetSubDetName() != "TIB" && point->GetSubDetName() != "TOB" && point->GetSubDetName() != "P1PXB"){
      cout<<"Weird detector:"<<point->GetSubDetName()<<endl;
    }
    
    if(removePionClusters){
      bool isPionHit = false;
      for(auto &pionCluster : event->GetPionClusters()){
        if(*pionCluster == *point){
          isPionHit = true;
          break;
        }
      }
      if(isPionHit) continue;
    }
    
    pointsNoEndcaps.push_back(point);
  }
  
  return pointsNoEndcaps;
}
