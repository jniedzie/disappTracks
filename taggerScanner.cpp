//  taggerScanner.cpp
//
//  Created by Jeremi Niedziela on 18/06/2019.

#include "Helpers.hpp"
#include "Fitter.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "PerformanceMonitor.hpp"
#include "EventSet.hpp"

string configPath = "configs/helixTagger.md";
string cutLevel = "after_L1/all/";//after_L1/";

const int nEvents = 100;
const int eventOffset = 100;

bool removeEndcapClusters = false;

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;

auto helixFitter = make_unique<Fitter>();

shared_ptr<Event> GetEvent(int iEvent);

vector<shared_ptr<Event>> events;
vector<Points> pointsNoEndcapsSignal;
vector<Points> pointsNoEndcapsBackground;

vector<string> monitorTypes = {
  "avg_hits",
  "max_hits",
  "avg_layers",
  "max_layers",
  "avg_length",
  "max_length",
  "n_helices"
};

vector<string> optimizationVars = {
//  "auc",
//  "sigma_init",
//  "sigma_L0",
//  "sigma_L1",
//  "max_eff",
//  "min_fake",
  "max_dist_fake",
  "avg_dist_fake",
};

map<string, map<string, double>> bestMonitorValue; //[monitorType][optimizeFor]
map<string, map<string, map<string, double>>> bestParamValue; //[monitorType][optimizeFor][configParam]

//           name     start    stop   step
vector<tuple<string, double, double, double>> paramsToTest = {
//  {"double_hit_max_distance", 20, 0, -1},
//  {"seed_max_chi2", 0.01, 0.10, 0.01},
//  {"seed_middle_hit_max_delta_phi", 0, 1.0, 0.1},
//  {"seed_middle_hit_max_delta_z", 50, 300, 50},
//  {"seed_last_hit_max_delta_phi", 0.0, 1.0, 0.1},
//  {"seed_last_hit_max_delta_z", 50, 300, 50},
  {"track_max_chi2", 0.001, 0.05, 0.002},
};

void CheckParamForCurrentConfig(string paramName)
{
  map<string, PerformanceMonitor> monitors;
  
  for(auto monitorType : monitorTypes){
    int max = 20, nBins = 20;
    if(monitorType=="avg_length" || monitorType=="max_length"){ max = 10; nBins = 40; }
    monitors[monitorType] = PerformanceMonitor(monitorType, nBins, 0, max);
  }

  for(auto iEvent=0; iEvent<nEvents; iEvent++){
    auto event = events[iEvent];
    
    for(auto &track : event->GetTracks()){
      Helices fittedHelicesSignal = helixFitter->FitHelices(pointsNoEndcapsSignal[iEvent], *track, *event->GetVertex());
      Helices fittedHelicesBackground = helixFitter->FitHelices(pointsNoEndcapsBackground[iEvent], *track, *event->GetVertex());
      
      for(string monitorType : monitorTypes){
        monitors[monitorType].SetValues(helixProcessor.GetHelicesParamsByMonitorName(fittedHelicesSignal, monitorType),
                                       helixProcessor.GetHelicesParamsByMonitorName(fittedHelicesBackground, monitorType));
        
      }
    }
  }
  
  for(string monitorType : monitorTypes){
    monitors[monitorType].CalcEfficiency();
    
    for(string optimizeFor : optimizationVars){
      double value = monitors[monitorType].GetValueByName(optimizeFor);
      if(value > bestMonitorValue[monitorType][optimizeFor]){
        bestMonitorValue[monitorType][optimizeFor] = value;
        bestParamValue[monitorType][optimizeFor][paramName] = config.params[paramName];
      }
    }
  }
};

int main(int argc, char* argv[])
{
  config = ConfigManager(configPath);
  
  for(auto iEvent=eventOffset; iEvent<eventOffset+nEvents; iEvent++){
    auto event = GetEvent(iEvent);
    events.push_back(event);
    pointsNoEndcapsSignal.push_back(event->GetClusters(false, removeEndcapClusters));
    pointsNoEndcapsBackground.push_back(event->GetClusters(true, removeEndcapClusters));
  }
  
  for(string monitorType : monitorTypes){
    for(string optimizeFor : optimizationVars){
      bestMonitorValue[monitorType][optimizeFor] = -inf;
      for(auto &[paramName, min, max, step] : paramsToTest){
        bestParamValue[monitorType][optimizeFor][paramName] = -inf;
      }
    }
  }
 
  // double hit max distance
  for(auto &[name, min, max, step] : paramsToTest){
    cout<<"Testing "<<name<<endl;
    
    for(config.params[name]  = min;
        step < 0 ? config.params[name] >= max : config.params[name] <= max;
        config.params[name] += step){
      cout<<"\t"<<config.params[name]<<endl;
      CheckParamForCurrentConfig(name);
    }
    cout<<endl;
  }
  
  for(string monitorType : monitorTypes){
    cout<<"\n\nMonitor type: "<<monitorType<<endl;
    
    for(string optimizeFor : optimizationVars){
      cout<<"\tbest "<<optimizeFor<<": "<<bestMonitorValue[monitorType][optimizeFor]<<endl;

      for(auto &[paramName, min, max, step] : paramsToTest){
        cout<<"\t\t"<<paramName<<": "<<bestParamValue[monitorType][optimizeFor][paramName]<<endl;
      }
    }
  }
  
//  ofstream outFile(outPath);
//  outFile << "Init "<<optimizeFor<<": "<<initOptValue<<"\tfinal "<<optimizeFor<<":"<<maxOptValue<<endl;
//  outFile << "double_hit_max_distance: "<<bestDoubleHitsDistance<<endl;
//  outFile << "seed_max_chi2: "<<bestSeedChi2<<endl;
//  outFile << "seed_middle_hit_min_delta_phi: "<<bestMiddleMinPhi<<endl;
//  outFile << "seed_middle_hit_max_delta_phi: "<<bestMiddleMaxPhi<<endl;
//  outFile << "seed_middle_hit_max_delta_z: "<<bestMiddleZ<<endl;
//  outFile << "seed_last_hit_min_delta_phi: "<<bestLastMinPhi<<endl;
//  outFile << "seed_last_hit_max_delta_phi: "<<bestLastMaxPhi<<endl;
//  outFile << "seed_last_hit_max_delta_z: "<<bestLastZ<<endl;
//  outFile << "track_max_chi2: "<<bestTrackChi2<<endl;
//  outFile << "next_point_max_delta_z: "<<bestTrackZ<<endl;
//  outFile << "next_point_max_delta_xy: "<<bestTrackXY<<endl;
//  outFile << "next_point_max_delta_t: "<<bestTrackT<<endl;
//
//  outFile.close();
  
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
