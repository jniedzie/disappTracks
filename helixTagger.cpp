//  helixTagger.cpp
//
//  Created by Jeremi Niedziela on 25/01/2019.

#include "Helpers.hpp"
#include "Fitter.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "EventSet.hpp"
#include "Logger.hpp"

string configPath = "configs/helixTagger.md";
string cutLevel = "after_L0/";

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kTaggerSignal;
//int setIter = kTaggerBackground;

vector<int> eventsToSkip = { };

int main(int argc, char* argv[])
{
  cout.imbue(locale("de_DE"));
  TApplication theApp("App", &argc, argv);
  
  config = ConfigManager(configPath);
  auto fitter = make_unique<Fitter>();
  
  EventSet events; events.LoadEventsFromFiles(cutLevel);
  
  auto start = now();
  int nAnalyzedEvents=0;
  
  for(auto iEvent=0; iEvent<events.size(dataType, setIter); iEvent++){
    if(find(eventsToSkip.begin(), eventsToSkip.end(), iEvent) != eventsToSkip.end()){
      Log(0)<<"\n\n=================================================================\n";
      Log(0)<<"Skipping event "<<iEvent<<"\n";
      continue;
    }
    
    Log(0)<<"\n\n=================================================================\n";
    Log(0)<<"helixTagger -- processing event "<<iEvent<<"\n";
    
    auto event = events.At(dataType, setIter, iEvent);
    
    for(auto &track : event->GetTracks()){
      Helices fittedHelices = fitter->FitHelices(event->GetClusters(), *track, *event->GetVertex());
      for(auto helix : fittedHelices) event->AddHelix(helix);
    }
    nAnalyzedEvents++;
  }
  
  Log(0)<<"N events analyzed: "<<nAnalyzedEvents<<"\n";
  Log(0)<<"Saving events...\n";
  events.SaveEventsToFiles("afterHelixTagging/");
  Log(0)<<"Time: "<<duration(start, now())<<"\n";
  
  return 0;
}
