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
  if(config.performCutsLevel == 0)       inputPrefix = "after_L0/";
  else if(config.performCutsLevel == 1)  inputPrefix = "after_L1/";
  else{
    cout<<"ERROR -- unknown cuts level: "<<config.performCutsLevel<<endl;
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
  
  if(config.drawStandardPlots) events.DrawStandardPlots();
  if(config.drawPerLayerPlots) events.DrawPerLayerPlots();
 
  cout<<"Done"<<endl;
  
  theApp.Run();
  return 0;
}
