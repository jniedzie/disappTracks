//
//  helixTagger.cpp
//
//  Created by Jeremi Niedziela on 25/01/2019.
//

#include "Helpers.hpp"
#include "Fitter.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "PerformanceMonitor.hpp"
#include "EventSet.hpp"

string configPath = "configs/helixTagger_maxLength.md";
string cutLevel = "after_L2/4layers/";//after_L1/";

int nEvents = 45;
const int nTests = 10;


int nAnalyzedEvents = 0;

double SetParamValue(int iTest){
//  return 0;
  // here put a way to calculate param value based on the test iter:
  // then assign the value to the correct config parameter:

//  double paramValue = iTest+5;
//  config.doubleHitsMaxDistance = paramValue;
  
  //-------
  // seeds
  
//  double paramValue = 0.01+iTest*0.01;//pow(10, 1-iTest);
//  config.seedMaxChi2 = paramValue;

//  double paramValue = -2.0+iTest*0.2;
//  config.seedMiddleHitDeltaPhi = range<double>(paramValue, config.seedMiddleHitDeltaPhi.GetMax());
//  double paramValue = -0.3+iTest*0.1;
//  config.seedMiddleHitDeltaPhi = range<double>(config.seedMiddleHitDeltaPhi.GetMin(), paramValue);
//  double paramValue = iTest*20;
//  config.seedMiddleHitMaxDeltaZ = paramValue;
  
//  double paramValue = -1.6+iTest*0.2;
//  config.seedLastHitDeltaPhi = range<double>(paramValue, config.seedLastHitDeltaPhi.GetMax());
//  double paramValue = -0.3+iTest*0.1;
//  config.seedLastHitDeltaPhi = range<double>(config.seedLastHitDeltaPhi.GetMin(), paramValue);
//  double paramValue = iTest*20;
//  config.seedLastHitMaxDeltaZ = paramValue;

  //-------
  // tracks
  
//  double paramValue = 0.01+iTest*0.01;// pow(10, -7+iTest);
//  config.trackMaxChi2 = paramValue;
  
//  double paramValue = -1.6+iTest*0.2;
//  config.nextPointDeltaPhi = range<double>(paramValue, config.nextPointDeltaPhi.GetMax());
//  double paramValue = -0.2+iTest*0.1;
//  config.nextPointDeltaPhi = range<double>(config.nextPointDeltaPhi.GetMin(), paramValue);
//  double paramValue = iTest*20;
//  config.nextPointMaxDeltaZ = paramValue;
//  double paramValue = iTest*20;
//  config.nextPointMaxDeltaXY = paramValue;
//  double paramValue = iTest+3;
//  config.trackMinNpoints = paramValue;
  
  //-------
  // merging
  
  double paramValue = iTest;
  config.mergingMaxDifferentPoints = paramValue;
//  double paramValue = iTest+3;
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

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  config = ConfigManager(configPath);
  auto fitter = make_unique<Fitter>();
  TCanvas *canvas = new TCanvas("ROC", "ROC", 2880,1800);
  canvas->Divide(4,3);
  
  map<string, PerformanceMonitor> monitors;
  
//  monitors["n_helices"]  = PerformanceMonitor("N helices",  20, 0, 20 , nTests, nEvents);
//  monitors["avg_hits"]   = PerformanceMonitor("Avg hits",   20, 0, 20 , nTests, nEvents);
//  monitors["max_hits"]   = PerformanceMonitor("Max hits",   20, 0, 20 , nTests, nEvents);
//  monitors["max_layers"] = PerformanceMonitor("Max layers", 20, 0, 20 , nTests, nEvents);
//  monitors["avg_length"] = PerformanceMonitor("Avg length", 20, 0, 2  , nTests, nEvents);
  monitors["max_length"] = PerformanceMonitor("Max length", 20, 0, 2  , nTests, nEvents);
  
  for(auto iEvent=0; iEvent<nEvents; iEvent++){
    auto event = GetEvent(iEvent);
    
    auto pointsNoEndcapsSignal     = GetClustersNoEndcaps(event, false);
    auto pointsNoEndcapsBackground = GetClustersNoEndcaps(event, true);
    
    if(pointsNoEndcapsSignal.size()==0 || pointsNoEndcapsBackground.size()==0){
      cout<<"helixTagger -- no tracker hits for event "<<iEvent<<endl;
      //      continue;
    }
    
    cout<<"\n\n=================================================================\n"<<endl;
    cout<<"helixTagger -- processing event "<<iEvent<<endl;
    
    for(int iTest=0; iTest<nTests; iTest++){
      SetParamValue(iTest);
      
      for(auto &track : event->GetTracks()){
        
        vector<Helix> fittedHelicesSignal     = fitter->FitHelices(pointsNoEndcapsSignal, *track, *event->GetVertex());
        vector<Helix> fittedHelicesBackground = fitter->FitHelices(pointsNoEndcapsBackground, *track, *event->GetVertex());
        
        
        // for(auto helix : fittedHelicesSignal) event->AddHelix(move(fittedHelix));
        
        
//        monitors["n_helices"].SetValues(iTest, iEvent,
//                                        fittedHelicesSignal.size(),
//                                        fittedHelicesBackground.size());
//
//        monitors["avg_hits"].SetValues(iTest, iEvent,
//                                       GetAvgNhits(fittedHelicesSignal),
//                                       GetAvgNhits(fittedHelicesBackground));
//
//        monitors["max_hits"].SetValues(iTest, iEvent,
//                                       GetMaxNhits(fittedHelicesSignal),
//                                       GetMaxNhits(fittedHelicesBackground));
//
//
//        monitors["max_layers"].SetValues(iTest, iEvent,
//                                         GetMaxNlayers(fittedHelicesSignal),
//                                         GetMaxNlayers(fittedHelicesBackground));
        
//        monitors["avg_length"].SetValues(iTest, iEvent,
//                                         GetAvgLength(fittedHelicesSignal),
//                                         GetAvgLength(fittedHelicesBackground));
        
        monitors["max_length"].SetValues(iTest, iEvent,
                                         GetMaxLength(fittedHelicesSignal),
                                         GetMaxLength(fittedHelicesBackground));
      }
    }
    nAnalyzedEvents++;
  }
  
  //  cout<<"helixTagger -- saving events"<<endl;
  //  events.SaveEventsToFiles("afterHelixTagging/");
  //  cout<<"helixTagger -- finished"<<endl;
  
  int iPad=1;
  for(auto &[name, monitor] : monitors){
    monitor.CalcEfficiency(SetParamValue, nAnalyzedEvents);
    canvas->cd(iPad++);
    monitor.DrawRocGraphs();
  }
  
  for(auto &[name, monitor] : monitors){
    canvas->cd(iPad++);
    monitor.DrawHists();
  }
  
  canvas->Update();
  theApp.Run();
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
    unordered_set<int> layers;
    for(int iPoint=0; iPoint<helix.GetNpoints(); iPoint++){
      int layer = helix.GetPoints()[iPoint]->GetLayer();
      if(layer > 0){
        if(iPoint < helix.GetFirstTurningPointIndex()) layers.insert( layer);
        else                                           layers.insert(-layer);
      }
    }
    if(layers.size() > maxNlayers) maxNlayers = (int)layers.size();
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
