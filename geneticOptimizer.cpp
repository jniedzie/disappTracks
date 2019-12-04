//  geneticOptimizer
//  Created by Jeremi Niedziela on 03/12/2019.

#include "Helpers.hpp"
#include "Fitter.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "PerformanceMonitor.hpp"
#include "EventSet.hpp"

using namespace TMVA;

#define CONVSTEPS 20
#define CONVCRIT 0.0001
#define SCSTEPS 10
#define SCRATE 5
#define SCFACTOR 0.95

const int nEvents = 50;
const int eventOffset = 0;
const double maxExecTimePerEvent = 1.0; // seconds

/*
 mind: big population sizes will help in searching the domain space of the solution but you have to weight this
 out to the number of generations the extreme case of 1 generation and populationsize n is equal to a Monte Carlo
 calculation with n tries
 */
int populationSize = 40;

auto helixFitter = make_unique<Fitter>(maxExecTimePerEvent);

vector<shared_ptr<Event>> loadedEvents;
vector<Points> pointsNoEndcapsSignal;
vector<Points> pointsNoEndcapsBackground;

double bestMonitorValue;
map<string, double> bestParamValue; //[configParam]

//           name     start    stop   nBins
vector<tuple<string, double, double, int>> paramsToTest = {
  {"double_hit_max_distance"      , 10    , 20  , 6   },
  {"seed_max_chi2"                , 0.0   , 0.01, 1e5 },
  {"seed_middle_hit_min_delta_phi",-1.0   , 0.0 , 11  },
  {"seed_middle_hit_max_delta_phi", 0.0   , 1.0 , 11  },
  {"seed_middle_hit_max_delta_z"  , 50    , 300 , 6   },
  {"seed_last_hit_min_delta_phi"  ,-1.0   , 0.0 , 11  },
  {"seed_last_hit_max_delta_phi"  , 0.0   , 1.2 ,  7  },
  {"seed_last_hit_max_delta_z"    , 50    , 300 , 6   },
  {"track_max_chi2"               , 0.0   , 0.01 , 1e5},
  {"next_point_min_delta_phi"     ,-1.0   , 0.0 , 11  },
  {"next_point_max_delta_phi"     , 0.0   , 1.5 , 16  },
  {"next_point_max_delta_z"       , 0     , 1000, 11  },
  {"next_point_max_delta_xy"      , 50    , 400 , 8   },
  {"next_point_max_delta_t"       , 0     , 3.0 , 31  },
  {"do_asymmetric_constraints"    , 0     , 2.0 , 2   },
};

void loadEventsAndClusters()
{
  EventSet events; events.LoadEventsFromFiles("after_L1/all/");
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    for(auto iEvent=eventOffset; iEvent<eventOffset+nEvents; iEvent++){
      auto event = events.At(kSignal, kTaggerSignalNoPU, year, iEvent);
      loadedEvents.push_back(event);
      config.params["fit_noise_clusters_only"] = false;
      pointsNoEndcapsSignal.push_back(event->GetClusters());
      config.params["fit_noise_clusters_only"] = true;
      pointsNoEndcapsBackground.push_back(event->GetClusters());
    }
  }
}

void setDefaultConfigParams()
{
  // This are non-optmizable options:
  config.runSignal[kTaggerSignalNoPU] = true;
  config.params["verbosity_level"] = 0;
  config.params["max_N_events_signal"] = nEvents+eventOffset;
  config.params["save_events"] = 0;
  config.params["load_2017"] = 0;
  config.params["load_2018"] = 1;
  config.params["load_friend_tree"] = 1;
  config.params["load_hits"] = 1;
  
  // These will be optimized
  config.params["double_hit_max_distance"] = 20.0;
  config.params["seed_max_chi2"] = 0.002;
  config.params["seed_middle_hit_max_delta_phi"] = 0.8;
  config.params["seed_middle_hit_max_delta_z"] = 150;
  config.params["seed_last_hit_max_delta_phi"] = 0.5;
  config.params["seed_last_hit_max_delta_z"] = 200;
  config.params["track_max_chi2"] = 0.011;
  config.params["next_point_max_delta_phi"] = 1.5;
  config.params["next_point_max_delta_z"] = 500;
  config.params["next_point_max_delta_xy"] = 50;
  config.params["next_point_max_delta_t"] = 1.5;
  
  config.params["do_asymmetric_constraints"] = 0;
  config.params["seed_middle_hit_min_delta_phi"] = -0.6;
  config.params["seed_last_hit_min_delta_phi"] = -0.6;
  config.params["next_point_min_delta_phi"] = -0.6;
  
  // Not optimizing those:
  config.params["cut_noise_hits"] = 1;
  config.params["include_endcaps"] = 0;
  
  config.params["track_min_n_points"] = 3;
  config.params["track_min_n_layers"] = 2;
  config.params["min_layers_for_delta_xy"] = 5;
  
  config.params["merging_max_different_point"] = 2;
  config.params["candidate_min_n_points"] = 3;
  config.params["merge_at_turn_back"] = 0;
  config.params["merge_final_helices"] = 1;
  
  config.params["max_n_missing_hits"] = 1;
  config.params["max_n_missing_hits_in_raw"] = 1;
  
  
  config.params["allow_turning_back"] = 1;
  config.params["require_good_starting_values"] = 1;
  config.params["exp_radius_function"] = 0;
  config.params["exp_slope_function"] = 0;
  config.params["allow_one_less_layer"] = 0;
  config.params["allow_one_more_layer"] = 1;
  config.params["check_opposite_charge_below_Nlayers"] = 5;
  
  config.params["start_R0"] =  350;
  config.params["min_R0"] =  50;
  config.params["max_R0"] =  1000;
  config.params["min_Rslope"] =   0;
  config.params["max_Rslope"] =  10000;
  config.params["min_S0"] =  -10000;
  config.params["max_S0"] =  10000;
  config.params["min_Sslope"] =  -1000;
  config.params["max_Sslope"] =  0;
  config.params["min_X0"] =  -5000;
  config.params["max_X0"] =  5000;
  config.params["min_Y0"] =  -5000;
  config.params["max_Y0"] =  5000;
  config.params["min_Z0"] =  -5000;
  config.params["max_Z0"] =  5000;
}

void printParams(vector<double> &params)
{
  for(int iParam=0; iParam < params.size(); iParam++){
    string name = get<0>(paramsToTest[iParam]);
    double value = params.at(iParam);
    cout<<"Setting param: "<<name<<"\t to value: "<<value<<endl;
  }
}

class MyFitness : public IFitterTarget {
public:
  MyFitness() : IFitterTarget() {}
  
  // the fitness-function goes here
  // the factors are optimized such that the return-value of this function is minimized
  // take care!! the fitness-function must never fail, .. means: you have to prevent
  // the function from reaching undefined values (such as x=0 for 1/x or so)
  //
  // HINT: to use INTEGER variables, it is sufficient to cast the "factor" in the fitness-function
  // to (int). In this case the variable-range has to be chosen +1 ( to get 0..5, take Interval(0,6) )
  // since the introduction of "Interval" ranges can be defined with a third parameter
  // which gives the number of bins within the interval. With that technique discrete values
  // can be achieved easier. The random selection out of this discrete numbers is completly uniform.
  //
  double EstimatorFunction(vector<double> &factors){
    
    // Set values of parameters in the config
    cout<<"\n----------------------------------------"<<endl;
    for(int iParam=0; iParam < factors.size(); iParam++){
      string name = get<0>(paramsToTest[iParam]);
      config.params[name] = factors.at(iParam);
      cout<<"Setting param: "<<name<<"\t to value: "<<factors.at(iParam)<<endl;
    }
    cout<<"----------------------------------------\n"<<endl;
    
    // Create monitor
    int max = 20, nBins = 20;
    string monitorType = "avg_hits";
    
    PerformanceMonitor monitor(monitorType, monitorType, nBins, 0, max);
    
    // Run reconstruction on loaded events
    auto start = now();
    cout<<"Event: ";
    for(auto iEvent=0; iEvent<loadedEvents.size(); iEvent++){
      auto event = loadedEvents[iEvent];
      cout<<iEvent<<" ";
      
      for(auto &track : event->GetTracks()){
        Helices fittedHelicesSignal     = helixFitter->FitHelices(pointsNoEndcapsSignal[iEvent], *track, *event->GetVertex());
        Helices fittedHelicesBackground = helixFitter->FitHelices(pointsNoEndcapsBackground[iEvent], *track, *event->GetVertex());
        monitor.SetValue(helixProcessor.GetHelicesParamsByMonitorName(fittedHelicesSignal, monitorType), true);
        monitor.SetValue(helixProcessor.GetHelicesParamsByMonitorName(fittedHelicesBackground, monitorType), false);
      }
    }
    cout<<endl;
    
    cout<<"\nTime per event: "<<duration(start, now())/loadedEvents.size()<<" s"<<endl;
    
    // Get max distance to √fake curve
    monitor.CalcEfficiency();
    double value = monitor.GetValueByName("max_dist_fake");
    cout<<"Distance to √fake: "<<value<<endl;
    cout<<endl;
    
    // return something that should be minimized
    return 1-value;
  }
};

vector<Interval*> getRanges()
{
  vector<Interval*> ranges;
  
  for(auto &[name, min, max, nBins] : paramsToTest){
    ranges.push_back(new Interval(min, max, nBins));
  }
  return ranges;
}

int main(int argc, char* argv[])
{
  cout<<"Starting genetic optimizer with params:"<<endl;
  cout<<"- N events: "<<nEvents<<endl;
  cout<<"- max time per event: "<<maxExecTimePerEvent<<endl;
  cout<<"- population size: "<<populationSize<<endl;
  
  setDefaultConfigParams();
  loadEventsAndClusters();
  
  IFitterTarget *myFitness = new MyFitness();
  vector<Interval*> ranges = getRanges();
  
  GeneticAlgorithm geneticOptimizer(*myFitness, populationSize, ranges);
  
  int iIter=0;
  
  do{
    cout<<"===================================================================="<<endl;
    cout<<"\tStarting iteration "<<iIter<<endl;
    cout<<"===================================================================="<<endl;
    iIter++;
    
    geneticOptimizer.Init();              // prepares the new generation and does evolution
    geneticOptimizer.CalculateFitness();  // assess the quality of the individuals
    geneticOptimizer.GetGeneticPopulation().Print(0);
    
    geneticOptimizer.GetGeneticPopulation().TrimPopulation(); // reduce the population size to the initially defined one
    
    // tricky thing: control the speed of how fast the "solution space" is searched through
    // this function basically influences the sigma of a gaussian around the actual value
    // of the parameter where the new value will be randomly thrown.
    // when the number of improvements within the last SCSTEPS
    // A) smaller than SCRATE: divide the preset sigma by SCFACTOR
    // B) equal to SCRATE: do nothing
    // C) greater than SCRATE: multiply the preset sigma by SCFACTOR
    // if you don't know what to do, leave it unchanged or even delete this function call
    geneticOptimizer.SpreadControl(SCSTEPS, SCRATE, SCFACTOR);
    
    GeneticGenes* genes = geneticOptimizer.GetGeneticPopulation().GetGenes(0);
    vector<double> gvec = genes->GetFactors();
    
    cout<<"BEST PARAMS in this iteration:"<<endl;
    printParams(gvec);
    
  }
  while(!geneticOptimizer.HasConverged(CONVSTEPS, CONVCRIT));
  // converged if: fitness-improvement < CONVCRIT within the last CONVSTEPS loops
  
  GeneticGenes* genes = geneticOptimizer.GetGeneticPopulation().GetGenes(0);
  vector<double> gvec = genes->GetFactors();
  
  cout<<"BEST PARAMS EVER:"<<endl;
  printParams(gvec);
  
  return 0;
}
