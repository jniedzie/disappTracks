//  taggerFitter.cpp
//
//  Created by Jeremi Niedziela on 18/06/2019.

#include "Helpers.hpp"
#include "Fitter.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "PerformanceMonitor.hpp"
#include "EventSet.hpp"

string configPath = "configs/helixTagger.md";
string outFilePrefix = "L2_all";
string outPath = "taggerScannerOutput.txt";
string cutLevel = "after_L2/all/";//after_L1/";
string optimizeFor = "min_fake";
string optimizationMonitor = "max_hits";

const int nEvents = 40;
const int eventOffset = 40;

int nAnalyzedEvents = 0;

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;

string monitorType = "";
string optimizationParam = "";

vector<shared_ptr<Point>> GetClustersNoEndcaps(const shared_ptr<Event> &event, bool removePionClusters);

auto helixFitter = make_unique<Fitter>();

double GetAvgNhits(vector<Helix> helices);
int    GetMaxNhits(vector<Helix> helices);
int    GetMaxNlayers(vector<Helix> helices);
double GetAvgLength(vector<Helix> helices);
double GetMaxLength(vector<Helix> helices);

shared_ptr<Event> GetEvent(int iEvent);

vector<shared_ptr<Event>> events;
vector<vector<shared_ptr<Point>>> pointsNoEndcapsSignal;
vector<vector<shared_ptr<Point>>> pointsNoEndcapsBackground;

double GetParamForCurrentConfig(string monitorType, string optParam)
{
  nAnalyzedEvents=0;
  int max = 20;
  int nBins = 20;
  if(monitorType=="avg_length") max = 2;
  if(monitorType=="max_length"){
    nBins = 40;
    max = 12;
  }
  auto monitor = PerformanceMonitor(monitorType, 20, 0, max, nEvents);
  
  for(auto iEvent=0; iEvent<nEvents; iEvent++){
    if(pointsNoEndcapsSignal[iEvent].size()==0 || pointsNoEndcapsBackground[iEvent].size()==0){
      cout<<"helixTagger -- no tracker hits for event "<<iEvent<<endl;
      //      continue;
    }
    auto event = events[iEvent];
    
    for(auto &track : event->GetTracks()){
      vector<Helix> fittedHelicesSignal = helixFitter->FitHelices(pointsNoEndcapsSignal[iEvent], *track, *event->GetVertex());
      vector<Helix> fittedHelicesBackground = helixFitter->FitHelices(pointsNoEndcapsBackground[iEvent], *track, *event->GetVertex());
      
      if(monitorType=="avg_hits"){
        monitor.SetValues(iEvent, GetAvgNhits(fittedHelicesSignal), GetAvgNhits(fittedHelicesBackground));
      }
      else if(monitorType=="max_hits"){
        monitor.SetValues(iEvent, GetMaxNhits(fittedHelicesSignal), GetMaxNhits(fittedHelicesBackground));
      }
      else if(monitorType=="max_layers"){
        monitor.SetValues(iEvent, GetMaxNlayers(fittedHelicesSignal), GetMaxNlayers(fittedHelicesBackground));
      }
      else if(monitorType=="n_helices"){
        monitor.SetValues(iEvent, fittedHelicesSignal.size(), fittedHelicesBackground.size());
      }
      else if(monitorType=="avg_length"){
        monitor.SetValues(iEvent, GetAvgLength(fittedHelicesSignal), GetAvgLength(fittedHelicesBackground));
      }
      else if(monitorType=="max_length"){
        monitor.SetValues(iEvent, GetMaxLength(fittedHelicesSignal), GetMaxLength(fittedHelicesBackground));
      }
      else{
        cout<<"Unknown monitor type: "<<monitorType<<endl;
        cout<<"Possible options: avg_hits/max_hits/max_layers/n_helices/avg_length/max_length"<<endl;
        exit(0);
      }

    }
    nAnalyzedEvents++;
  }
  monitor.CalcEfficiency(nAnalyzedEvents);
  
  if(optParam=="auc")               return monitor.GetAUC();
  else if(optParam=="sigma_init")   return monitor.GetSignificanceInitial();
  else if(optParam=="sigma_L0")     return monitor.GetSignificanceAfterL0();
  else if(optParam=="sigma_L1")     return monitor.GetSignificanceAfterL1();
  else if(optParam=="max_eff")      return monitor.GetMaxEfficiency();
  else if(optParam=="min_fake")     return monitor.GetInvFakeAtHighestEff();
  else{
    cout<<"Uknown optimization parameter: "<<optParam<<endl;
  }
  return 0.0;
};

int main(int argc, char* argv[])
{
  config = ConfigManager(configPath);
  
  for(auto iEvent=eventOffset; iEvent<eventOffset+nEvents; iEvent++){
    auto event = GetEvent(iEvent);
    events.push_back(event);
    pointsNoEndcapsSignal.push_back(GetClustersNoEndcaps(event, false));
    pointsNoEndcapsBackground.push_back(GetClustersNoEndcaps(event, true));
  }
 
  double initOptValue = GetParamForCurrentConfig(optimizationMonitor, optimizeFor);
  
  double maxOptValue = -inf;
  
  // double hit max distance
//  double bestDoubleHitsDistance = config.doubleHitsMaxDistance;
//  for(config.doubleHitsMaxDistance = 20.0;
//      config.doubleHitsMaxDistance >= 5.0;
//      config.doubleHitsMaxDistance -= 1.0){
//    double currentAUC = GetParamForCurrentConfig(optimizationMonitor, optimizeFor);
//    cout<<"double_hit_max_distance: "<<config.doubleHitsMaxDistance<<"\t"<<optimizeFor<<": "<<currentAUC<<endl;
//
//    if(currentAUC > maxOptValue){
//      maxOptValue = currentAUC;
//      bestDoubleHitsDistance = config.doubleHitsMaxDistance;
//    }
//  }
//  config.doubleHitsMaxDistance = bestDoubleHitsDistance;
  
//   seed chi2
  double bestSeedChi2 = config.seedMaxChi2;
  for(double exponent=-3; exponent<=1; exponent+=1){
    config.seedMaxChi2 = pow(10, exponent);
    double currentAUC = GetParamForCurrentConfig(optimizationMonitor, optimizeFor);
    cout<<"seed_max_chi2: "<<config.seedMaxChi2<<"\t"<<optimizeFor<<": "<<currentAUC<<endl;

    if(currentAUC > maxOptValue){
      maxOptValue = currentAUC;
      bestSeedChi2 = config.seedMaxChi2;
    }
  }
  config.seedMaxChi2 = bestSeedChi2;
  
  // seed middle min Δφ
//  double bestMiddleMinPhi = config.seedMiddleHitDeltaPhi.GetMin();
//  for(double min = 0.0; min >= -1.0; min-=0.2){
//    config.seedMiddleHitDeltaPhi = range<double>(min, config.seedMiddleHitDeltaPhi.GetMax());
//
//    double currentAUC = GetParamForCurrentConfig(optimizationMonitor, optimizeFor);
//    cout<<"seed_middle_hit_min_delta_phi: "<<min<<"\t"<<optimizeFor<<": "<<currentAUC<<endl;
//
//    if(currentAUC > maxOptValue){
//      maxOptValue = currentAUC;
//      bestMiddleMinPhi = min;
//    }
//  }
//  config.seedMiddleHitDeltaPhi = range<double>(bestMiddleMinPhi, config.seedMiddleHitDeltaPhi.GetMax());

  // seed middle max Δφ
//  double bestMiddleMaxPhi = config.seedMiddleHitDeltaPhi.GetMax();
//  for(double max = 0.0; max <= 1.0; max+=0.2){
//    config.seedMiddleHitDeltaPhi = range<double>(config.seedMiddleHitDeltaPhi.GetMin(), max);
//
//    double currentAUC = GetParamForCurrentConfig(optimizationMonitor, optimizeFor);
//    cout<<"seed_middle_hit_max_delta_phi: "<<max<<"\t"<<optimizeFor<<": "<<currentAUC<<endl;
//
//    if(currentAUC > maxOptValue){
//      maxOptValue = currentAUC;
//      bestMiddleMaxPhi = max;
//    }
//  }
//  config.seedMiddleHitDeltaPhi = range<double>(config.seedMiddleHitDeltaPhi.GetMin(), bestMiddleMaxPhi);
  
  
  // seed middle Δz
//  double bestMiddleZ = config.seedMiddleHitMaxDeltaZ;
//  for(config.seedMiddleHitMaxDeltaZ  = 0;
//      config.seedMiddleHitMaxDeltaZ <= 200;
//      config.seedMiddleHitMaxDeltaZ += 10){
//    double currentAUC = GetParamForCurrentConfig(optimizationMonitor, optimizeFor);
//    cout<<"seed_middle_hit_max_delta_z: "<<config.seedMiddleHitMaxDeltaZ<<"\t"<<optimizeFor<<": "<<currentAUC<<endl;
//
//    if(currentAUC > maxOptValue){
//      maxOptValue = currentAUC;
//      bestMiddleZ = config.seedMiddleHitMaxDeltaZ;
//    }
//  }
//  config.seedMiddleHitMaxDeltaZ = bestMiddleZ;
  
  /*
  // seed last min Δφ
  double bestLastMinPhi = config.seedLastHitDeltaPhi.GetMin();
  for(double min = -0.4; min >= -0.7; min-=0.1){
    config.seedLastHitDeltaPhi = range<double>(min, config.seedLastHitDeltaPhi.GetMax());
    
    double currentAUC = GetParamForCurrentConfig(optimizationMonitor, optimizeFor);
    cout<<"seed_Last_hit_min_delta_phi: "<<min<<"\t"<<optimizeFor<<": "<<currentAUC<<endl;
    
    if(currentAUC > maxOptValue){
      maxOptValue = currentAUC;
      bestLastMinPhi = min;
    }
  }
  config.seedLastHitDeltaPhi = range<double>(bestLastMinPhi, config.seedLastHitDeltaPhi.GetMax());
  
  // seed last max Δφ
  double bestLastMaxPhi = config.seedLastHitDeltaPhi.GetMax();
  for(double max = 0.2; max <= 0.7; max+=0.1){
    config.seedLastHitDeltaPhi = range<double>(config.seedLastHitDeltaPhi.GetMin(), max);
    
    double currentAUC = GetParamForCurrentConfig(optimizationMonitor, optimizeFor);
    cout<<"seed_last_hit_max_delta_phi: "<<max<<"\t"<<optimizeFor<<": "<<currentAUC<<endl;
    
    if(currentAUC > maxOptValue){
      maxOptValue = currentAUC;
      bestLastMaxPhi = max;
    }
  }
  config.seedLastHitDeltaPhi = range<double>(config.seedLastHitDeltaPhi.GetMin(), bestLastMaxPhi);
   */
  
  // seed last min Δz
//  double bestLastZ = config.seedLastHitMaxDeltaZ;
//  for(config.seedLastHitMaxDeltaZ = 150;
//      config.seedLastHitMaxDeltaZ <= 300;
//      config.seedLastHitMaxDeltaZ += 10){
//    double currentAUC = GetParamForCurrentConfig(optimizationMonitor, optimizeFor);
//    cout<<"seed_last_hit_max_delta_z: "<<config.seedLastHitMaxDeltaZ<<"\t"<<optimizeFor<<": "<<currentAUC<<endl;
//
//    if(currentAUC > maxOptValue){
//      maxOptValue = currentAUC;
//      bestLastZ = config.seedLastHitMaxDeltaZ;
//    }
//  }
//  config.seedLastHitMaxDeltaZ = bestLastZ;
  
  // track chi2
//  double bestTrackChi2 = config.trackMaxChi2;
//  for(double exponent=-3.0; exponent<-1.9; exponent+=0.1){
//    config.trackMaxChi2 = pow(10, exponent);
//    double currentAUC = GetParamForCurrentConfig(optimizationMonitor, optimizeFor);
//    cout<<"track_max_chi2: "<<config.trackMaxChi2<<"\t"<<optimizeFor<<": "<<currentAUC<<endl;
//
//    if(currentAUC > maxOptValue){
//      maxOptValue = currentAUC;
//      bestTrackChi2 = config.trackMaxChi2;
//    }
//  }
//  config.trackMaxChi2 = bestTrackChi2;
  
  // next point Δz
//  double bestTrackZ = config.nextPointMaxDeltaZ;
//  for(config.nextPointMaxDeltaZ = 200;
//      config.nextPointMaxDeltaZ <= 1000;
//      config.nextPointMaxDeltaZ += 100){
//    double currentAUC = GetParamForCurrentConfig(optimizationMonitor, optimizeFor);
//    cout<<"next_point_max_delta_z: "<<config.nextPointMaxDeltaZ<<"\t"<<optimizeFor<<": "<<currentAUC<<endl;
//
//    if(currentAUC > maxOptValue){
//      maxOptValue = currentAUC;
//      bestTrackZ = config.nextPointMaxDeltaZ;
//    }
//  }
//  config.nextPointMaxDeltaZ = bestTrackZ;
  
  // next point Δxy
//  double bestTrackXY = config.nextPointMaxDeltaXY;
//  for(config.nextPointMaxDeltaXY  = 0;
//      config.nextPointMaxDeltaXY <= 100;
//      config.nextPointMaxDeltaXY += 10){
//    double currentAUC = GetParamForCurrentConfig(optimizationMonitor, optimizeFor);
//    cout<<"next_point_max_delta_xy: "<<config.nextPointMaxDeltaXY<<"\t"<<optimizeFor<<": "<<currentAUC<<endl;
//
//    if(currentAUC > maxOptValue){
//      maxOptValue = currentAUC;
//      bestTrackXY = config.nextPointMaxDeltaXY;
//    }
//  }
//  config.nextPointMaxDeltaXY = bestTrackXY;
  
  // next point Δxy
//  double bestTrackT = config.nextPointMaxDeltaT;
//  for(config.nextPointMaxDeltaT  = 0;
//      config.nextPointMaxDeltaT <= 1.5;
//      config.nextPointMaxDeltaT += 0.1){
//    double currentAUC = GetParamForCurrentConfig(optimizationMonitor, optimizeFor);
//    cout<<"next_point_max_delta_t: "<<config.nextPointMaxDeltaT<<"\t"<<optimizeFor<<": "<<currentAUC<<endl;
//    
//    if(currentAUC > maxOptValue){
//      maxOptValue = currentAUC;
//      bestTrackT = config.nextPointMaxDeltaT;
//    }
//  }
//  config.nextPointMaxDeltaT = bestTrackT;

//  SetParameter(fitter, 13, "track_min_n_points"            ,  5   ,  0    , 20  , 1  , dontFix);
//  SetParameter(fitter, 14, "merging_max_different_point"   ,  5   ,  0    , 20  , 1  , dontFix);
//  SetParameter(fitter, 15, "max_n_missing_hits"            ,  1   ,  0    , 5   , 1  , fix);
//  SetParameter(fitter, 16, "max_n_missing_hits_in_raw"     ,  1   ,  0    , 5   , 1  , fix);
//  SetParameter(fitter, 17, "merge_at_turn_back"            ,  0   ,  0    , 2   , 1  , fix);
//  SetParameter(fitter, 18, "merge_final_helices"           ,  0   ,  0    , 2   , 1  , fix);
//  SetParameter(fitter, 19, "do_asymmetric_constraints"     ,  1   ,  0    , 2   , 1  , fix);
//  SetParameter(fitter, 20, "allow_turning_back"            ,  1   ,  0    , 2   , 1  , fix);
//  SetParameter(fitter, 21, "candidate_min_n_points"        ,  3   ,  3    , 10  , 1  , dontFix);
  
  ofstream outFile(outPath);
  outFile << "Output for data: " << cutLevel << "\t optimizing for: "<<optimizeFor<<"\tusing monitor: "<<optimizationMonitor<<endl;
  outFile << "Init "<<optimizeFor<<": "<<initOptValue<<"\tfinal "<<optimizeFor<<":"<<maxOptValue<<endl;
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
  outFile.close();
  
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
