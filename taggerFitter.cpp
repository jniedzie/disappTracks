//  taggerFitter.cpp
//
//  Created by Jeremi Niedziela on 18/06/2019.

#include "Helpers.hpp"
#include "Fitter.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "PerformanceMonitor.hpp"
#include "EventSet.hpp"

string configPath = "configs/helixTagger_maxHits.md";
string cutLevel = "after_L2/4layers/";//after_L1/";

int nEvents = 10;

int nAnalyzedEvents = 0;

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

void SetParameter(ROOT::Fit::Fitter *fitter, int i, string name, double start, double min, double max, bool fix=false)
{
  fitter->Config().ParSettings(i).SetName(name);
  fitter->Config().ParSettings(i).SetValue(start);
  fitter->Config().ParSettings(i).SetLimits((min < max) ? min : max,(min < max) ? max : min);
  fitter->Config().ParSettings(i).SetStepSize(0.0001);
  if(fix) fitter->Config().ParSettings(i).Fix();
}

void SetParameter(TFitter *fitter, int i, string name, double start, double min, double max, bool fix=false)
{
  fitter->SetParameter(i, name.c_str(), start, 0.1, min, max);
  if(fix) fitter->FixParameter(i);
}

auto helixFitter = make_unique<Fitter>();
vector<shared_ptr<Event>> events;
vector<vector<shared_ptr<Point>>> pointsNoEndcapsSignal;
vector<vector<shared_ptr<Point>>> pointsNoEndcapsBackground;

//auto chi2Function = [&](const double *par) {
void chi2Function(Int_t&, Double_t*, Double_t &f, Double_t *par, Int_t)
{
  config.doubleHitsMaxDistance     = par[0];
  config.seedMaxChi2               = par[1];
  config.seedMiddleHitDeltaPhi     = range<double>(par[2], par[3]);
  config.seedMiddleHitMaxDeltaZ    = par[4];
  config.seedLastHitDeltaPhi       = range<double>(par[5], par[6]);
  config.seedLastHitMaxDeltaZ      = par[7];
  config.trackMaxChi2              = par[8];
  config.nextPointDeltaPhi         = range<double>(par[9], par[10]);
  config.nextPointMaxDeltaZ        = par[11];
  config.nextPointMaxDeltaXY       = par[12];
  config.trackMinNpoints           = par[13];
  config.mergingMaxDifferentPoints = par[14];
  config.maxNmissingHits           = par[15];
  config.maxNmissingHitsInRow      = par[16];
  config.mergeAtTurnBack           = par[17];
  config.mergeFinalHelices         = par[18];
  config.doAsymmetricConstraints   = par[19];
  config.allowTurningBack          = par[20];
  
  nAnalyzedEvents=0;
  auto monitor = PerformanceMonitor("Max hits", 20, 0, 20 , nEvents);
  
  for(auto iEvent=0; iEvent<nEvents; iEvent++){
    
    if(pointsNoEndcapsSignal[iEvent].size()==0 || pointsNoEndcapsBackground[iEvent].size()==0){
      cout<<"helixTagger -- no tracker hits for event "<<iEvent<<endl;
      //      continue;
    }
    auto event = events[iEvent];
    
    for(auto &track : event->GetTracks()){
      vector<Helix> fittedHelicesSignal = helixFitter->FitHelices(pointsNoEndcapsSignal[iEvent], *track, *event->GetVertex());
      vector<Helix> fittedHelicesBackground = helixFitter->FitHelices(pointsNoEndcapsBackground[iEvent], *track, *event->GetVertex());
      monitor.SetValues(iEvent, GetMaxNhits(fittedHelicesSignal), GetMaxNhits(fittedHelicesBackground));
    }
    nAnalyzedEvents++;
  }
  monitor.CalcEfficiency(nAnalyzedEvents);
  //    double auc = monitor.GetAUC();
  double significance = monitor.GetSignificanceAfterL0();
  double chi2 = pow(20-significance, 2)/significance;
  
  cout<<"chi2: "<<chi2<<"\tsignificance: "<<significance<<endl;
  f = chi2;
};

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = ConfigManager(configPath);
  
  for(auto iEvent=0; iEvent<nEvents; iEvent++){
    auto event = GetEvent(iEvent);
    events.push_back(event);
    pointsNoEndcapsSignal.push_back(GetClustersNoEndcaps(event, false));
    pointsNoEndcapsBackground.push_back(GetClustersNoEndcaps(event, true));
  }
  

  const int nPar = 21;
  TVirtualFitter::SetDefaultFitter("Minuit");
  TFitter *fitter = new TFitter(nPar);
  fitter->SetFCN(chi2Function);
  
  bool fix = true;
  bool dontFix = false;
  
  SetParameter(fitter, 0 , "double_hit_max_distance"       ,  16.0,  5.0  , 20  , fix);
  SetParameter(fitter, 1 , "seed_max_chi2"                 ,  3   ,  1e-10, 20  , fix);
  SetParameter(fitter, 2 , "seed_middle_hit_min_delta_phi" , -0.1 , -1.0  , 0.0 , dontFix);
  SetParameter(fitter, 3, "seed_middle_hit_max_delta_phi" ,  0.1 ,  0.0  , 1.0  , fix);
  SetParameter(fitter, 4 , "seed_middle_hit_max_delta_z"   ,  60  ,  0.0  , 300 , fix);
  SetParameter(fitter, 5 , "seed_last_hit_min_delta_phi"   , -0.8 , -1.0  ,-0.00001 , fix);
  SetParameter(fitter, 6 , "seed_last_hit_max_delta_phi"   ,  0.1 ,  0.0  , 1.0 , fix);
  SetParameter(fitter, 7 , "seed_last_hit_max_delta_z"     ,  180 ,  0.0  , 300 , fix);
  SetParameter(fitter, 8 , "track_max_chi2"                ,  1e-3,  1e-10, 2e-3, fix);
  SetParameter(fitter, 9 , "next_point_min_delta_phi"      , -0.5 , -1.0  , 0.0 , fix);
  SetParameter(fitter, 10, "next_point_max_delta_phi"      ,  0.5 ,  0.0  , 1.0 , fix);
  SetParameter(fitter, 11, "next_point_max_delta_z"        ,  220 ,  0.0  , 300 , fix);
  SetParameter(fitter, 12, "next_point_max_delta_xy"       ,  20  ,  0.0  , 300 , fix);
  SetParameter(fitter, 13, "track_min_n_points"            ,  3   ,  3    , 20  , fix);
  SetParameter(fitter, 14, "merging_max_different_point"   ,  2   ,  0    , 20  , fix);
  
  SetParameter(fitter, 15, "max_n_missing_hits"            ,  1   ,  0    , 5   , fix);
  SetParameter(fitter, 16, "max_n_missing_hits_in_raw"     ,  1   ,  0    , 5   , fix);
  SetParameter(fitter, 17, "merge_at_turn_back"            ,  0   ,  0    , 2   , fix);
  SetParameter(fitter, 18, "merge_final_helices"           ,  0   ,  0    , 2   , fix);
  SetParameter(fitter, 19, "do_asymmetric_constraints"     ,  1   ,  0    , 2   , fix);
  SetParameter(fitter, 20, "allow_turning_back"            ,  1   ,  0    , 2   , fix);
  
  double args = 0; // put to 0 for results only, or to -1 for no garbage
  fitter->ExecuteCommand( "SET PRINTOUT"  , &args, 1);
  fitter->ExecuteCommand( "SET PRINT"     , &args, 1);
  double strategyLevel[1] = {2};
  fitter->ExecuteCommand( "SET STR", strategyLevel, 1);
  double arglist[1] = {0};
  fitter->ExecuteCommand("MIGRAD", arglist, 0);

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
