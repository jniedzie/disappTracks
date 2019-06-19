//  taggerFitter.cpp
//
//  Created by Jeremi Niedziela on 18/06/2019.

#include "Helpers.hpp"
#include "Fitter.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "PerformanceMonitor.hpp"
#include "EventSet.hpp"

string configPath = "../configs/helixTagger_maxHits.md";
string outFilePrefix = "L2_4layers";
string cutLevel = "after_L2/4layers/";//after_L1/";

int nEvents = 45;

int nAnalyzedEvents = 0;

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;

string monitorType = "";
string optimizationParam = "";

vector<shared_ptr<Point>> GetClustersNoEndcaps(const shared_ptr<Event> &event, bool removePionClusters);
TF1* GetRocFunction();

double GetAvgNhits(vector<Helix> helices);
int    GetMaxNhits(vector<Helix> helices);
int    GetMaxNlayers(vector<Helix> helices);
double GetAvgLength(vector<Helix> helices);
double GetMaxLength(vector<Helix> helices);

shared_ptr<Event> GetEvent(int iEvent);

void WriteOutputToFile(const TFitter *fitter, string outPath);

void SetParameter(TFitter *fitter, int i,
                  string name, double start, double min, double max, double step,
                  bool fix=false)
{
  fitter->SetParameter(i, name.c_str(), start, step, min, max);
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
  config.seedMaxChi2               = pow(10, par[1]);
  config.seedMiddleHitDeltaPhi     = range<double>(par[2], par[3]);
  config.seedMiddleHitMaxDeltaZ    = par[4];
  config.seedLastHitDeltaPhi       = range<double>(par[5], par[6]);
  config.seedLastHitMaxDeltaZ      = par[7];
  config.trackMaxChi2              = pow(10, par[8]);
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
  config.candidateMinNpoints       = par[21];
  
  nAnalyzedEvents=0;
  int max = 20;
  if(monitorType=="avg_length") max = 2;
  if(monitorType=="max_length") max = 6;
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
  
  double optValue = inf;
  double chi2 = inf;
  
  if(optimizationParam=="auc"){
    optValue = monitor.GetAUC();
    chi2 = pow(1-optValue, 2)/optValue;
  }
  else if(optimizationParam=="max_eff"){
    optValue = monitor.GetMaxEfficiency();
    chi2 = pow(1-optValue, 2)/optValue;
  }
  else if(optimizationParam=="sigma_init"){
    optValue = monitor.GetSignificanceInitial();
    chi2 = pow(20-optValue, 2)/optValue;
  }
  else if(optimizationParam=="sigma_L0"){
    optValue = monitor.GetSignificanceAfterL0();
    chi2 = pow(20-optValue, 2)/optValue;
  }
  else{
    cout<<"Unknown optimization param: "<<optimizationParam<<endl;
    cout<<"Possible options: auc/max_eff/sigma_init/sigma_L0"<<endl;
    exit(0);
  }

  cout<<"chi2: "<<chi2<<"\t"<<optimizationParam<<": "<<optValue<<endl;
  f = chi2;
};

int main(int argc, char* argv[])
{
  if(argc!=3){
    cout<<"Usage: ./taggerFitter [avg_hits/max_hits/max_layers/n_helices/avg_length/max_length] [auc/max_eff/sigma_init/sigma_L0]"<<endl;
    exit(0);
  }
  monitorType = argv[1];
  optimizationParam = argv[2];
  
  config = ConfigManager(configPath);
  
  for(auto iEvent=0; iEvent<nEvents; iEvent++){
    auto event = GetEvent(iEvent);
    events.push_back(event);
    pointsNoEndcapsSignal.push_back(GetClustersNoEndcaps(event, false));
    pointsNoEndcapsBackground.push_back(GetClustersNoEndcaps(event, true));
  }
  

  const int nPar = 22;
  TVirtualFitter::SetDefaultFitter("Minuit");
  TFitter *fitter = new TFitter(nPar);
  fitter->SetFCN(chi2Function);
  
  bool fix = true;
  bool dontFix = false;
  
//                     i    name                            start  min     max   step  fix
  SetParameter(fitter, 0 , "double_hit_max_distance"       ,  10.0,  5.0  , 20  , 1.0, dontFix);
  SetParameter(fitter, 1 , "seed_max_chi2"                 ,  -1  ,  -10  , 1   , 1  , dontFix);
  SetParameter(fitter, 2 , "seed_middle_hit_min_delta_phi" , -0.5 , -2.0  ,-0.0 , 0.1, dontFix);
  SetParameter(fitter, 3 , "seed_middle_hit_max_delta_phi" ,  0.1 ,  0.0  , 1.0 , 0.1, dontFix);
  SetParameter(fitter, 4 , "seed_middle_hit_max_delta_z"   ,  20  ,  0.0  , 300 , 10 , dontFix);
  SetParameter(fitter, 5 , "seed_last_hit_min_delta_phi"   , -0.5 , -2.0  ,-0.1 , 0.1, dontFix);
  SetParameter(fitter, 6 , "seed_last_hit_max_delta_phi"   ,  0.0 , -0.1  , 2.0 , 0.1, dontFix);
  SetParameter(fitter, 7 , "seed_last_hit_max_delta_z"     ,  20  ,  0.0  , 300 , 10 , dontFix);
  SetParameter(fitter, 8 , "track_max_chi2"                ,  -2  ,  -6   , -3  , 1  , dontFix);
  SetParameter(fitter, 9 , "next_point_min_delta_phi"      , -0.5 , -1.0  , 0.0 , 0.1, fix);
  SetParameter(fitter, 10, "next_point_max_delta_phi"      ,  0.5 ,  0.0  , 1.0 , 0.1, fix);
  SetParameter(fitter, 11, "next_point_max_delta_z"        ,  20  ,  0.0  , 300 , 10 , dontFix);
  SetParameter(fitter, 12, "next_point_max_delta_xy"       ,  5   ,  0.0  , 200 , 10 , dontFix);
  SetParameter(fitter, 13, "track_min_n_points"            ,  5   ,  0    , 20  , 1  , dontFix);
  SetParameter(fitter, 14, "merging_max_different_point"   ,  5   ,  0    , 20  , 1  , dontFix);
  SetParameter(fitter, 15, "max_n_missing_hits"            ,  1   ,  0    , 5   , 1  , fix);
  SetParameter(fitter, 16, "max_n_missing_hits_in_raw"     ,  1   ,  0    , 5   , 1  , fix);
  SetParameter(fitter, 17, "merge_at_turn_back"            ,  0   ,  0    , 2   , 1  , fix);
  SetParameter(fitter, 18, "merge_final_helices"           ,  0   ,  0    , 2   , 1  , fix);
  SetParameter(fitter, 19, "do_asymmetric_constraints"     ,  1   ,  0    , 2   , 1  , fix);
  SetParameter(fitter, 20, "allow_turning_back"            ,  1   ,  0    , 2   , 1  , fix);
  SetParameter(fitter, 21, "candidate_min_n_points"        ,  3   ,  3    , 10  , 1  , dontFix);
  
  double args = 0; // put to 0 for results only, or to -1 for no garbage
  fitter->ExecuteCommand( "SET PRINTOUT"  , &args, 1);
  fitter->ExecuteCommand( "SET PRINT"     , &args, 1);
//  double strategyLevel[1] = {2};
//  fitter->ExecuteCommand( "SET STR", strategyLevel, 1);
  double arglist[1] = {0};
  fitter->ExecuteCommand("MIGRAD", arglist, 0);

  string outPath = "taggerFitterOutput_"+outFilePrefix+"_"+monitorType+"_"+optimizationParam+".txt";
  
  WriteOutputToFile(fitter, outPath);
  
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

void WriteOutputToFile(const TFitter *fitter, string outPath)
{
  ofstream outFile(outPath);
  outFile << "Output for data: " << cutLevel << "\t monitor: "<<monitorType<<"\toptimizing: "<<optimizationParam<<endl;
  
  double chi2, edm, errdef;
  int nvpar, nparx;
  
  fitter->GetStats(chi2, edm, errdef, nvpar, nparx);
  
  outFile << "chi2: " << chi2 <<"\tedm: "<<edm<<endl;
  
  outFile << endl;
  outFile << "double_hit_max_distance: "        << fitter->GetParameter(0)  << endl;
  outFile << endl;
  outFile << "seed_max_chi2: "                  << fitter->GetParameter(1)  << endl;
  outFile << endl;
  outFile << "seed_middle_hit_min_delta_phi: "  << fitter->GetParameter(2)  << endl;
  outFile << "seed_middle_hit_max_delta_phi: "  << fitter->GetParameter(3)  << endl;
  outFile << "seed_middle_hit_max_delta_z: "    << fitter->GetParameter(4)  << endl;
  outFile << endl;
  outFile << "seed_last_hit_min_delta_phi: "    << fitter->GetParameter(5)  << endl;
  outFile << "seed_last_hit_max_delta_phi: "    << fitter->GetParameter(6)  << endl;
  outFile << "seed_last_hit_max_delta_z: "      << fitter->GetParameter(7)  << endl;
  outFile << endl;
  outFile << "track_max_chi2: "                 << fitter->GetParameter(8)  << endl;
  outFile << endl;
  outFile << "next_point_min_delta_phi: "       << fitter->GetParameter(9)  << endl;
  outFile << "next_point_max_delta_phi: "       << fitter->GetParameter(10) << endl;
  outFile << "next_point_max_delta_z: "         << fitter->GetParameter(11) << endl;
  outFile << "next_point_max_delta_xy: "        << fitter->GetParameter(12) << endl;
  outFile << endl;
  outFile << "track_min_n_points: "             << fitter->GetParameter(13) << endl;
  outFile << endl;
  outFile << "merging_max_different_point: "    << fitter->GetParameter(14) << endl;
  outFile << "candidate_min_n_points: "         << fitter->GetParameter(21) << endl;
  outFile << "merge_at_turn_back: "             << fitter->GetParameter(17) << endl;
  outFile << "merge_final_helices: "            << fitter->GetParameter(18) << endl;
  
  outFile << "max_n_missing_hits: "             << fitter->GetParameter(15) << endl;
  outFile << "max_n_missing_hits_in_raw: "      << fitter->GetParameter(16) << endl;
  outFile << endl;
  outFile << "do_asymmetric_constraints: "      << fitter->GetParameter(19) << endl;
  outFile << endl;
  outFile << "allow_turning_back: "             << fitter->GetParameter(20) << endl;
  
  outFile.close();
}
