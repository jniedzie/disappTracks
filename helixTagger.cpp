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

xtracks::EDataType dataType = xtracks::kSignal;

vector<int> eventsToSkip = { };

int main(int argc, char* argv[])
{
  if(argc != 1 && argc != 4){
    Log(0)<<"helixTagger takes no arguments or: output_path events_offset n_events \n";
    exit(0);
  }
  cout.imbue(locale("de_DE"));
  TApplication theApp("App", &argc, argv);
  
  config = ConfigManager(configPath);
  
  string cutLevel;
  if(config.params["cuts_level"]==0) cutLevel = "after_L0/";
    if(config.params["cuts_level"]==1) cutLevel = "after_L1/"+config.category+"/";
  EventSet events;
//  events.LoadEventsFromFiles(cutLevel);
  
  int eventOffset   = 0;
  string outputPath = "afterHelixTagging";
  int maxEvents     = 0;
  
  if(argc == 4){
    outputPath  = argv[1];
    eventOffset = atoi(argv[2]);
    maxEvents   = atoi(argv[3]);
  }
  
  auto fitter = make_unique<Fitter>();
  auto start = now();
  int nAnalyzedEvents=0;
  
  for(int iSig=0; iSig<kNsignals; iSig++){
    if(!config.runSignal[iSig]) continue;
    
//    if(argc == 1) maxEvents = events.size(dataType, iSig);
    
    for(auto iEvent=eventOffset; iEvent<maxEvents+eventOffset; iEvent++){
      if(find(eventsToSkip.begin(), eventsToSkip.end(), iEvent) != eventsToSkip.end()){
        Log(0)<<"\n\n=================================================================\n";
        Log(0)<<"Skipping event "<<iEvent<<"\n";
        continue;
      }
      
      Log(0)<<"\n\n=================================================================\n";
      Log(0)<<"helixTagger -- processing event "<<iEvent<<"\n";
      Log(0)<<"cut level: "<<cutLevel<<"\n";
      
      events.LoadEventFromFiles(dataType, iSig, iEvent, cutLevel);
      auto event = events.At(dataType, iSig, iEvent-eventOffset);
      
      if(!event->HasFriendData()){
        Log(0)<<"Warning -- skipping event "<<iEvent<<" as it has no friend info\n";
        event->SetWasTagged(false);
        continue;
      }
      
      for(auto &track : event->GetTracks()){
        Helices fittedHelices = fitter->FitHelices(event->GetClusters(), *track, *event->GetVertex());
        
        for(auto helix : fittedHelices){
          helix.Print();
          event->AddHelix(helix);
        }
      }
      event->SetWasTagged(true);
      nAnalyzedEvents++;
    }
    
    Log(0)<<"N events analyzed: "<<nAnalyzedEvents<<"\n";
    
    if(config.params["save_events"]){
      Log(0)<<"Saving events...\n";
      events.SaveEventsToFiles(outputPath+"/");
    }
    Log(0)<<"Time: "<<duration(start, now())<<"\n";

  }
  
  return 0;
}
