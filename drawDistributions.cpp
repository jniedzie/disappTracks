#include "EventSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "CutsManager.hpp"

#include <TApplication.h>

string configPath = "configs/analysis.md";

int main(int argc, char* argv[])
{
  TApplication *theApp = new TApplication("App", &argc, argv);
  
  // All events with initial cuts only
  config = ConfigManager(configPath);
  auto events = make_shared<EventSet>();
  
  string inputPrefix;
  if(config.performCutsLevel == 0)       inputPrefix = "after_L0/";
  else if(config.performCutsLevel == 1)  inputPrefix = "after_L1/";
  else{
    cout<<"ERROR -- unknown cuts level: "<<config.performCutsLevel<<endl;
    exit(0);
  }
  
  events->LoadEventsFromFiles(inputPrefix);
  
  CutsManager cutsManager;
  
  EventCut eventCut;
  TrackCut trackCut;
  JetCut jetCut;
  LeptonCut leptonCut;
  
  cutsManager.GetCuts(eventCut, trackCut, jetCut, leptonCut);
  
  events->ApplyCuts(eventCut,
                    trackCut,
                    make_unique<JetCut>(jetCut),
                    make_unique<LeptonCut>(leptonCut));
  
  if(config.drawStandardPlots) events->DrawStandardPlots();
  if(config.drawPerLayerPlots) events->DrawPerLayerPlots();
 
  theApp->Run();
  return 0;
}
