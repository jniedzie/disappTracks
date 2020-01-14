//  taggerScanner.cpp
//
//  Created by Jeremi Niedziela on 18/06/2019.

#include "Helpers.hpp"
#include "Fitter.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "PerformanceMonitor.hpp"
#include "EventSet.hpp"

//string configPath = "configs/helixTagger.md";
string configPath = "configs/helixTaggerTunedOnTruth.md";

const int nEvents = 50;
const int eventOffset = 0;

string suffix = "";

xtracks::EDataType dataType = kSignal;

auto helixFitter = make_unique<Fitter>();

vector<shared_ptr<Event>> loadedEvents;
vector<Points> pointsNoEndcapsSignal;
vector<Points> pointsNoEndcapsBackground;


vector<string> monitorTypes = {
//  "avg_hits",
//  "max_hits",
//  "avg_layers",
//  "max_layers",
//  "avg_length",
//  "max_length",
//  "n_helices"
};

vector<string> optimizationVars = {
//  "auc",
//  "sigma_init",
//  "sigma_L0",
//  "sigma_L1",
//  "max_eff",
//  "min_fake",
//  "max_dist_fake",
//  "avg_dist_fake",
};

map<string, map<string, double>> bestMonitorValue; //[monitorType][optimizeFor]
map<string, map<string, map<string, double>>> bestParamValue; //[monitorType][optimizeFor][configParam]
double bestAvgByFraction = 0;
double bestMaxByFraction = 0;
double bestPionFraction = 0;
double latestAvgByFraction = 0;
double latestPionFraction = 0;
double latestMaxByFraction = 0;
double bestParamValueForAvg = -inf;
double bestParamValueForMax = -inf;
double bestParamValueForPionFraction = -inf;

//           name     start    stop   step   log
vector<tuple<string, double, double, double, bool>> paramsToTest = {
//  {"double_hit_max_distance"      , 20    , 20    , -1    , false },
  
//  {"seed_max_chi2"                , 1E-7  , 0.1   , 0.01  , true  },
//  {"seed_middle_hit_min_delta_phi", 0     , -1.0  , -0.1  , false },
//  {"seed_middle_hit_max_delta_phi", 0.0   , 1.0   , 0.1   , false },
//  {"seed_middle_hit_max_delta_z"  , 0     , 250   , 10    , false },
//  {"seed_last_hit_min_delta_phi"  , 0.0   , -1.0  , -0.1  , false },
  {"seed_last_hit_max_delta_phi"  , -0.3   , 0.5   , 0.1   , false },
//  {"seed_last_hit_max_delta_z"    , 0     , 200   , 10    , false },
  
//  {"track_max_chi2"               , 1E-8  , 1.00  , 0.001 , true  },
//  {"next_point_min_delta_phi"     , 0.0   , -1.0  , -0.1  , false },
//  {"next_point_max_delta_phi"     , -0.4  , 0.2   , 0.1   , false },
//  {"next_point_max_delta_z"       , 0     , 400   , 10    , false },
//  {"next_point_max_delta_xy"      , 150   , 150   , 5     , false },
//  {"next_point_max_delta_t"       , 0     , 2.0   , 0.1   , false },
};

/// Returns path prefix for cuts level and category selected in the config file
string getPathPrefix()
{
  string prefix = "";
   
  if(config.secondaryCategory == "Zmumu") prefix += "Zmumu/";
  if(config.secondaryCategory == "Wmunu") prefix += "Wmunu/";
  
  if(config.params["cuts_level"]==0) prefix += "after_L0/";
  if(config.params["cuts_level"]==1) prefix += "after_L1/"+config.category+"/";
  
  prefix += suffix+"/";
  
  return prefix;
}

void CheckParamForCurrentConfig(string paramName)
{
  map<string, PerformanceMonitor> monitors;
  
  for(auto monitorType : monitorTypes){
    int max = 20, nBins = 20;
    if(monitorType=="avg_length" || monitorType=="max_length"){ max = 10; nBins = 40; }
    monitors[monitorType] = PerformanceMonitor(monitorType, monitorType, nBins, 0, max);
  }
  cout<<"Event: ";
  
  double avgOfAvgByFraction=0;
  double avgOfMaxByFraction=0;
  
  double totalPionHitsOnHelix=0;
  double totalHitsOnHelix=0;
  
  for(auto iEvent=0; iEvent<loadedEvents.size(); iEvent++){
    auto event = loadedEvents[iEvent];
    if(iEvent%5==0) cout<<iEvent<<endl;
    
    double avgByFraction=0;
    double maxByFraction=0;
    int nHelices=0;
    
    for(auto &track : event->GetTracks()){
      Helices fittedHelicesSignal = helixFitter->FitHelices(pointsNoEndcapsSignal[iEvent], *track, *event->GetVertex());
//      Helices fittedHelicesBackground = helixFitter->FitHelices(pointsNoEndcapsBackground[iEvent], *track, *event->GetVertex());
      
//      for(string monitorType : monitorTypes){
//        monitors[monitorType].SetValue(helixProcessor.GetHelicesParamsByMonitorName(fittedHelicesSignal, monitorType), true);
//        monitors[monitorType].SetValue(helixProcessor.GetHelicesParamsByMonitorName(fittedHelicesBackground, monitorType), false);
//      }

      for(Helix helix : fittedHelicesSignal){
        double byFraction = helix.GetNtruePionHits()/((double)helix.GetNrecHits());
        avgByFraction += byFraction;
        if(byFraction > maxByFraction) maxByFraction = byFraction;
        
        totalHitsOnHelix      += helix.GetNrecHits();
        totalPionHitsOnHelix  += helix.GetNtruePionHits();
      }
      nHelices += fittedHelicesSignal.size();
      
    }
    avgByFraction /= nHelices;
    if(nHelices==0) avgByFraction=0;

    avgOfAvgByFraction += avgByFraction;
    avgOfMaxByFraction += maxByFraction;
  }
  
  double pionHitsFraction = totalPionHitsOnHelix/totalHitsOnHelix;
  
  avgOfAvgByFraction /= loadedEvents.size();
  avgOfMaxByFraction /= loadedEvents.size();
  
  if(avgOfAvgByFraction > bestAvgByFraction){
    bestAvgByFraction = avgOfAvgByFraction;
    bestParamValueForAvg = config.params[paramName];
  }
  if(avgOfMaxByFraction > bestMaxByFraction){
    bestMaxByFraction = avgOfMaxByFraction;
    bestParamValueForMax = config.params[paramName];
  }
  if(pionHitsFraction > bestPionFraction){
    bestPionFraction = pionHitsFraction;
    bestParamValueForPionFraction = config.params[paramName];
  }
  latestAvgByFraction = avgOfAvgByFraction;
  latestMaxByFraction = avgOfMaxByFraction;
  latestPionFraction = pionHitsFraction;
  
  cout<<endl;
  
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

void loadEventsAndClusters()
{
  EventSet events; events.LoadEventsFromFiles(getPathPrefix());
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    for(ESignal iSig : signals){
      if(!config.runSignal[iSig]) continue;
      
      for(auto iEvent=eventOffset; iEvent<eventOffset+nEvents; iEvent++){
        auto event = events.At(dataType, iSig, year, iEvent);
        loadedEvents.push_back(event);
        config.params["fit_noise_clusters_only"] = false;
        pointsNoEndcapsSignal.push_back(event->GetClusters());
        config.params["fit_noise_clusters_only"] = true;
        pointsNoEndcapsBackground.push_back(event->GetClusters());
      }
    }
  }
}

void initMonitors()
{
  for(string monitorType : monitorTypes){
     for(string optimizeFor : optimizationVars){
       bestMonitorValue[monitorType][optimizeFor] = -inf;
       for(auto &[paramName, min, max, step, doLog] : paramsToTest){
         bestParamValue[monitorType][optimizeFor][paramName] = -inf;
       }
     }
   }
}

int main(int argc, char* argv[])
{
  config = ConfigManager(configPath);
  
  loadEventsAndClusters();
  initMonitors();
  
  
  for(auto &[name, min, max, step, doLog] : paramsToTest){
    cout<<"Testing "<<name<<endl;
    
    for(config.params[name]  = min;
        step < 0 ? config.params[name] >= max : config.params[name] <= max;
        doLog ?  config.params[name] *= 10 : config.params[name] += step){
      cout<<"\t"<<config.params[name]<<endl;
      
      CheckParamForCurrentConfig(name);
//      cout<<"Best monitor value at the moment: "<<bestMonitorValue[monitorTypes[0]][optimizationVars[0]]<<endl;
      cout<<"latest avg: "<<latestAvgByFraction<<endl;
      cout<<"latest max: "<<latestMaxByFraction<<endl;
      cout<<"latest pion fraction: "<<latestPionFraction<<endl;
      cout<<"Best avg so far: "<<bestAvgByFraction<<endl;
      cout<<"Best max so far: "<<bestMaxByFraction<<endl;
      cout<<"Best pion fraction so far: "<<bestPionFraction<<endl;
    }
    cout<<endl;
    
//    if(bestParamValue[monitorTypes[0]][optimizationVars[0]][name] > -9999){
//      config.params[name] = bestParamValue[monitorTypes[0]][optimizationVars[0]][name];
//    }
    
    if(bestParamValueForPionFraction > -9999){
      config.params[name] = bestParamValueForPionFraction;
    }
    
    cout<<"Best avg: "<<bestAvgByFraction<<endl;
    cout<<"Best max: "<<bestMaxByFraction<<endl;
    cout<<"Best pion fraction: "<<bestPionFraction<<endl;
    cout<<"Best param value (avg): "<<bestParamValueForAvg<<endl;
    cout<<"Best param value (max): "<<bestParamValueForMax<<endl;
    cout<<"Best param value (pion fraction): "<<bestParamValueForPionFraction<<endl;
    
//    cout<<"Best monitor value:"<<bestMonitorValue[monitorTypes[0]][optimizationVars[0]]<<endl;
//    cout<<"Best param value "<<name<<": "<<bestParamValue[monitorTypes[0]][optimizationVars[0]][name]<<endl;
    
  }
  
  for(string monitorType : monitorTypes){
    cout<<"\n\nMonitor type: "<<monitorType<<endl;
    
    for(string optimizeFor : optimizationVars){
      cout<<"\tbest "<<optimizeFor<<": "<<bestMonitorValue[monitorType][optimizeFor]<<endl;

      for(auto &[paramName, min, max, step, doLog] : paramsToTest){
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
