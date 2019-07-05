#include "EventSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "CutsManager.hpp"

#include <TApplication.h>

string configPath = "configs/analysis.md";

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  // All events with initial cuts only
  config = ConfigManager(configPath);
  EventSet events;
  
  string inputPrefix;
  if(config.params["cuts_level"] == 0)       inputPrefix = "after_L0/";
  else if(config.params["cuts_level"] == 1)  inputPrefix = "after_L1/";
  else{
    cout<<"ERROR -- unknown cuts level: "<<config.params["cuts_level"]<<endl;
    exit(0);
  }
  
  events.LoadEventsFromFiles(inputPrefix);
  
  CutsManager cutsManager;
  
  EventCut eventCut;
  TrackCut trackCut;
  JetCut jetCut;
  LeptonCut leptonCut;
  
  cutsManager.GetCuts(eventCut, trackCut, jetCut, leptonCut);
  events.ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  cout<<"Drawing plots"<<endl;
  
  if(config.params["draw_standard_plots"]) events.DrawStandardPlots();
  if(config.params["draw_per_layer_plots"]) events.DrawPerLayerPlots();
 
  cout<<"Done"<<endl;
  
  theApp.Run();
  return 0;
}
