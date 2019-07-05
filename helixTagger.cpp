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
string cutLevel = "after_L1/all/";//after_L1/";

const int nEvents = 100; // max: 1287
const int eventOffset = 0;

int nAnalyzedEvents = 0;

vector<string> monitorTypes = {
  "avg_hits",
  "max_hits",
  "avg_layers",
  "max_layers",
  "avg_length",
  "max_length",
  "n_helices"
};

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;

shared_ptr<Event> GetEvent(int iEvent);
vector<shared_ptr<Point>> GetClustersNoEndcaps(const shared_ptr<Event> &event, bool removePionClusters);

int main(int argc, char* argv[])
{
  cout.imbue(locale("de_DE"));
  // TApplication theApp("App", &argc, argv);
  
  config = ConfigManager(configPath);
  auto fitter = make_unique<Fitter>();
  TCanvas *canvas = new TCanvas("ROC", "ROC", 2880,1800);
  canvas->Divide(3,3);
  
  map<string, PerformanceMonitor> monitors;
  
  for(auto monitorType : monitorTypes){
    int max = 20, nBins = 20;
    if(monitorType=="avg_length" || monitorType=="max_length"){ max = 10; nBins = 40; }
    monitors[monitorType] = PerformanceMonitor(monitorType, nBins, 0, max, nEvents);
  }
  
  vector<shared_ptr<Event>> events;
  vector<vector<shared_ptr<Point>>> pointsNoEndcapsSignal;
  vector<vector<shared_ptr<Point>>> pointsNoEndcapsBackground;
  
  for(auto iEvent=eventOffset; iEvent<eventOffset+nEvents; iEvent++){
    auto event = GetEvent(iEvent);
    events.push_back(event);
    pointsNoEndcapsSignal.push_back(GetClustersNoEndcaps(event, false));
    pointsNoEndcapsBackground.push_back(GetClustersNoEndcaps(event, true));
  }
  
  auto start = now();
  
  nAnalyzedEvents=0;
  
  for(auto iEvent=0; iEvent<nEvents; iEvent++){
    cout<<"\n\n=================================================================\n"<<endl;
    cout<<"helixTagger -- processing event "<<iEvent<<endl;
    
    auto event = events[iEvent];
    
    for(auto &track : event->GetTracks()){
      
      vector<Helix> fittedHelicesSignal     = fitter->FitHelices(pointsNoEndcapsSignal[iEvent], *track, *event->GetVertex());
      vector<Helix> fittedHelicesBackground = fitter->FitHelices(pointsNoEndcapsBackground[iEvent], *track, *event->GetVertex());
      
      // for(auto helix : fittedHelicesSignal) event->AddHelix(move(fittedHelix));
      
      for(string monitorType : monitorTypes){
        monitors[monitorType].SetValues(iEvent,
                                        helixProcessor.GetHelicesParamsByMonitorName(fittedHelicesSignal, monitorType),
                                        helixProcessor.GetHelicesParamsByMonitorName(fittedHelicesBackground, monitorType));
        
      }
    }
    nAnalyzedEvents++;
  }
  
  int iPad=1;
  for(auto &monitor : monitors){
    monitor.second.CalcEfficiency(nAnalyzedEvents);
    canvas->cd(iPad++);
    monitor.second.DrawRocGraph(true);
    canvas->cd(iPad++);
    monitor.second.DrawHists();
  }
  
  
  for(auto &monitor : monitors){
    cout<<"\n\n============================================================"<<endl;
    cout<<"Monitor: "<<monitor.first<<endl;
    monitor.second.PrintFakesEfficiency();
    monitor.second.PrintParams();
  }
  
  //  cout<<"helixTagger -- saving events"<<endl;
  //  events.SaveEventsToFiles("afterHelixTagging/");
  //  cout<<"helixTagger -- finished"<<endl;
  
  cout<<"Time: "<<duration(start, now());
  
  canvas->Update();
  //  theApp.Run();
  return 0;
}

shared_ptr<Event> GetEvent(int iEvent)
{
  EventSet events;
  events.LoadEventFromFiles(dataType, setIter, iEvent, cutLevel);
  auto event = events.At(dataType, setIter, 0);
  
  if(!event){
    cout<<"helixTagger -- event not found"<<endl;
    exit(0);
  }
  
  return event;
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
