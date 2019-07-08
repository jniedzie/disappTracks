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
  
  vector<map<string, PerformanceMonitor>> monitors; // [paramBin][monitorType]
  const int nParamBins = 4;
  
  for(int iParam=0; iParam<nParamBins; iParam++){
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
  
  nAnalyzedEvents=0;
  
  for(auto iEvent=0; iEvent<nEvents; iEvent++){
    cout<<"\n\n=================================================================\n"<<endl;
    cout<<"helixTagger -- processing event "<<iEvent<<endl;
    
//    auto event = events[iEvent];
    auto event = events.At(dataType, setIter, iEvent);
    
    if(event->GetGenPionHelices().size() != 1) continue;
    double pionPt = event->GetGenPionHelices().front().GetMomentum()->GetTransverse();
    
    int iParam = -1;
    if(pionPt < 150) iParam = 0;
    if(pionPt > 150 && pionPt < 300) iParam = 1;
    if(pionPt > 300 && pionPt < 450) iParam = 2;
    if(pionPt > 450) iParam = 3;
    
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
      monitor.second.PrintFakesEfficiency();
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
