//  helixTagger.cpp
//
//  Created by Jeremi Niedziela on 25/01/2019.

#include "Helpers.hpp"
#include "Fitter.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "PerformanceMonitor.hpp"
#include "EventSet.hpp"

string configPath = "../configs/helixTagger.md";
string cutLevel = "after_L1/all/";//after_L1/";

const int nEvents = 100; // max: 1287
const int eventOffset = 0;
const int nTests = 1;

int nAnalyzedEvents = 0;

double SetParamValue(int iTest){
  return 0;
  // here put a way to calculate param value based on the test iter:
  // then assign the value to the correct config parameter:

  double paramValue = iTest+5;
  config.doubleHitsMaxDistance = paramValue;
  
  //-------
  // seeds
  
//  double paramValue = 0.1+iTest*0.1;//pow(10, -3+iTest);
//  config.seedMaxChi2 = paramValue;

//  double paramValue = -1.0+iTest*0.2;
//  config.seedMiddleHitDeltaPhi = range<double>(paramValue, config.seedMiddleHitDeltaPhi.GetMax());
//  double paramValue = -0.5+iTest*0.2;
//  config.seedMiddleHitDeltaPhi = range<double>(config.seedMiddleHitDeltaPhi.GetMin(), paramValue);
//  double paramValue = iTest*30;
//  config.seedMiddleHitMaxDeltaZ = paramValue;
  
//  double paramValue = -1.0+iTest*0.1;
//  config.seedLastHitDeltaPhi = range<double>(paramValue, config.seedLastHitDeltaPhi.GetMax());
//  double paramValue = -0.2+iTest*0.1;
//  config.seedLastHitDeltaPhi = range<double>(config.seedLastHitDeltaPhi.GetMin(), paramValue);
//  double paramValue = iTest*30;
//  config.seedLastHitMaxDeltaZ = paramValue;

  //-------
  // tracks
  
//  double paramValue =  pow(10, -5+iTest*0.5); // 0.01+iTest*0.01;//
//  config.trackMaxChi2 = paramValue;
  
//  double paramValue = -1.5+iTest*0.2;
//  config.nextPointDeltaPhi = range<double>(paramValue, config.nextPointDeltaPhi.GetMax());
//  double paramValue = -0.5+iTest*0.2;
//  config.nextPointDeltaPhi = range<double>(config.nextPointDeltaPhi.GetMin(), paramValue);
//  double paramValue = iTest*20;
//  config.nextPointMaxDeltaZ = paramValue;
//  double paramValue = iTest*20;
//  config.nextPointMaxDeltaXY = paramValue;
//  double paramValue = iTest*0.1;
//  config.nextPointMaxDeltaT = paramValue;
  
//  double paramValue = iTest+3;
//  config.trackMinNpoints = paramValue;
  
  //-------
  // merging
  
//  double paramValue = 5+iTest;
//  config.mergingMaxDifferentPoints = paramValue;
//  double paramValue = iTest;
//  config.candidateMinNpoints = paramValue;

  //-------
  // switches
  
//  bool paramValue = iTest;
//  config.doAsymmetricConstraints = paramValue;

  return paramValue;
}

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;

shared_ptr<Event> GetEvent(int iEvent);
vector<shared_ptr<Point>> GetClustersNoEndcaps(const shared_ptr<Event> &event, bool removePionClusters);
TF1* GetRocFunction();

double GetAvgNhits(vector<Helix> helices);
int    GetMaxNhits(vector<Helix> helices);
int    GetMaxNlayers(vector<Helix> helices);
double GetAvgLength(vector<Helix> helices);
double GetMaxLength(vector<Helix> helices);
double GetMinChi2(vector<Helix> helices);
double GetMinChi2overNhits(vector<Helix> helices);
bool   DidTurnBack(vector<Helix> helices);

int main(int argc, char* argv[])
{
  cout.imbue(locale("de_DE"));
  // TApplication theApp("App", &argc, argv);
  
  config = ConfigManager(configPath);
  auto fitter = make_unique<Fitter>();
  TCanvas *canvas = new TCanvas("ROC", "ROC", 2880,1800);
  canvas->Divide(3,3);
  
  vector<map<string, PerformanceMonitor>> monitors;// [iTest][name]
  
  for(int iTest=0; iTest<nTests; iTest++){
    map<string, PerformanceMonitor> mapForTest = {
      {"n_helices" , PerformanceMonitor("N helices",  20, 0, 20 , nEvents)},
//      {"avg_hits"  , PerformanceMonitor("Avg hits",   20, 0, 20 , nEvents)},
      {"max_hits"  , PerformanceMonitor("Max hits",   20, 0, 20 , nEvents)},
      {"max_layers", PerformanceMonitor("Max layers", 20, 0, 20 , nEvents)},
//      {"avg_length", PerformanceMonitor("Avg length", 20, 0, 2  , nEvents)},
      {"max_length", PerformanceMonitor("Max length", 40, 0, 12   , nEvents)},
//      {"min_chi2"  , PerformanceMonitor("Min chi2"  , 10000, 0, 0.01 , nEvents)},
//      {"min_chi2_per_hit"  , PerformanceMonitor("Min chi2 per hit"  , 10000, 0, 0.001 , nEvents)},
      {"turned_back"  , PerformanceMonitor("Did turn back"  , 2, 0, 2 , nEvents)},
    };
    monitors.push_back(mapForTest);
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
  for(int iTest=0; iTest<nTests; iTest++){
    cout<<"\n\nparam: "<<SetParamValue(iTest)<<"\n\n"<<endl;
    nAnalyzedEvents=0;
    
    for(auto iEvent=0; iEvent<nEvents; iEvent++){
    
      if(pointsNoEndcapsSignal[iEvent].size()==0 || pointsNoEndcapsBackground[iEvent].size()==0){
        cout<<"helixTagger -- no tracker hits for event "<<iEvent<<endl;
        //      continue;
      }
      cout<<"\n\n=================================================================\n"<<endl;
      cout<<"helixTagger -- processing event "<<iEvent<<endl;
      
      auto event = events[iEvent];
      
      for(auto &track : event->GetTracks()){
        
        vector<Helix> fittedHelicesSignal     = fitter->FitHelices(pointsNoEndcapsSignal[iEvent], *track, *event->GetVertex());
        vector<Helix> fittedHelicesBackground = fitter->FitHelices(pointsNoEndcapsBackground[iEvent], *track, *event->GetVertex());

        // for(auto helix : fittedHelicesSignal) event->AddHelix(move(fittedHelix));
        
        
        monitors[iTest]["n_helices"].SetValues(iEvent,
                                               fittedHelicesSignal.size(),
                                               fittedHelicesBackground.size());
        
//        monitors[iTest]["avg_hits"].SetValues(iEvent,
//                                              GetAvgNhits(fittedHelicesSignal),
//                                              GetAvgNhits(fittedHelicesBackground));
        
        monitors[iTest]["max_hits"].SetValues(iEvent,
                                              GetMaxNhits(fittedHelicesSignal),
                                              GetMaxNhits(fittedHelicesBackground));
        
        monitors[iTest]["max_layers"].SetValues(iEvent,
                                                GetMaxNlayers(fittedHelicesSignal),
                                                GetMaxNlayers(fittedHelicesBackground));
        
//        monitors[iTest]["avg_length"].SetValues(iEvent,
//                                                GetAvgLength(fittedHelicesSignal),
//                                                GetAvgLength(fittedHelicesBackground));
        
        monitors[iTest]["max_length"].SetValues(iEvent,
                                                GetMaxLength(fittedHelicesSignal),
                                                GetMaxLength(fittedHelicesBackground));
        
//        monitors[iTest]["min_chi2"].SetValues(iEvent,
//                                              GetMinChi2(fittedHelicesSignal),
//                                              GetMinChi2(fittedHelicesBackground));
//        
//        monitors[iTest]["min_chi2_per_hit"].SetValues(iEvent,
//                                                      GetMinChi2overNhits(fittedHelicesSignal),
//                                                      GetMinChi2overNhits(fittedHelicesBackground));
        
        
        monitors[iTest]["turned_back"].SetValues(iEvent,
                                                 DidTurnBack(fittedHelicesSignal),
                                                 DidTurnBack(fittedHelicesBackground));
        
      }
      nAnalyzedEvents++;
    }
  
    int iPad=1;
    for(auto &monitor : monitors[iTest]){
      monitor.second.CalcEfficiency(nAnalyzedEvents);
      canvas->cd(iPad++);
      monitor.second.DrawRocGraph(iTest==0);
      canvas->cd(iPad++);
      monitor.second.DrawHists();
    }
    
  }
  
  
  for(auto &monitor : monitors[0]){
    cout<<"\n\n============================================================"<<endl;
    cout<<"Monitor: "<<monitor.first<<endl;

    if(monitor.first!="min_chi2" && monitor.first!="min_chi2_per_hit"){
      for(int iTest=0; iTest<nTests; iTest++){
        double param = SetParamValue(iTest);
        cout<<"Param: "<<param<<endl;
        monitors[iTest][monitor.first].PrintFakesEfficiency();
      }
    }
    for(int iTest=0; iTest<nTests; iTest++){
      double param = SetParamValue(iTest);
      cout<<"Param: "<<param<<"\t";
      monitors[iTest][monitor.first].PrintParams();
    }
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

double GetAvgNhits(vector<Helix> helices)
{
  if(helices.size()==0) return 0;
  double avgHits = 0;
  for(auto helix : helices) avgHits += helix.GetNpoints();
  avgHits /= helices.size();
  return avgHits;
}

int GetMaxNhits(vector<Helix> helices)
{
  int maxNhits = 0;
  for(auto helix : helices){
    if(helix.GetNpoints() > maxNhits) maxNhits = helix.GetNpoints();
  }
  return maxNhits;
}

int GetMaxNlayers(vector<Helix> helices)
{
  int maxNlayers = 0;
  
  for(auto helix : helices){
    uint nLayers = helix.GetNlayers();
    if(nLayers > maxNlayers) maxNlayers = nLayers;
  }
  return maxNlayers;
}

double GetAvgLength(vector<Helix> helices)
{
  if(helices.size()==0) return 0;
  double avgLength = 0;
  
  for(auto helix : helices){
    double length = fabs(helix.GetTmax() - helix.GetTmin());
    avgLength += length;
  }
  avgLength /= helices.size();
  return avgLength;
}

double GetMaxLength(vector<Helix> helices)
{
  if(helices.size()==0) return 0;
  double maxLength = -inf;
  
  for(auto helix : helices){
    double length = fabs(helix.GetTmax() - helix.GetTmin());
    if(length > maxLength) maxLength = length;
  }
  return maxLength;
}

double GetMinChi2(vector<Helix> helices)
{
  if(helices.size()==0) return 0;
  double minChi2 = inf;
  
  for(auto helix : helices){
    if(helix.GetChi2() < minChi2) minChi2 = helix.GetChi2();
  }
  return minChi2;
}

double GetMinChi2overNhits(vector<Helix> helices)
{
  if(helices.size()==0) return 0;
  double minChi2 = inf;
  
  for(auto helix : helices){
    if(helix.GetChi2()/helix.GetNpoints() < minChi2) minChi2 = helix.GetChi2()/helix.GetNpoints();
  }
  
  return minChi2;
}

bool DidTurnBack(vector<Helix> helices)
{
  if(helices.size()==0) return false;
  
  for(auto helix : helices){
    if(helix.GetFirstTurningPointIndex()>0) return true;
  }
  
  return false;
}
