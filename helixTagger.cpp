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
  if(argc != 1 && argc != 5){
    Log(0)<<"helixTagger takes no arguments or: output_path events_offset n_events config_path\n";
    exit(0);
  }
  cout.imbue(locale("de_DE"));
//  TApplication theApp("App", &argc, argv);
  
  if(argc == 5) configPath = argv[4];
  
  cout<<"Reading config from "<<configPath<<endl;
  config = ConfigManager(configPath);
  
  string cutLevel;
  if(config.params["cuts_level"]==0) cutLevel = "after_L0/";
  if(config.params["cuts_level"]==1) cutLevel = "after_L1/"+config.category+"/";
  EventSet events;
  
  int eventOffset   = 0;
  string outputPath = cutLevel+"afterHelixTagging";
  int maxEvents     = config.params["max_N_events_signal"];
  
  if(argc == 5){
    outputPath  = argv[1];
    eventOffset = atoi(argv[2]);
    maxEvents   = atoi(argv[3]);
  }
  
  auto fitter = make_unique<Fitter>();
  auto start = now();
  int nAnalyzedEvents=0;
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    cout<<"Runnig for year: "<<year<<endl;
    
    for(ESignal iSig : signals){
      if(!config.runSignal[iSig]) continue;
      cout<<"Running for signal: "<<iSig<<endl;
      
      for(auto iEvent=eventOffset; iEvent<maxEvents+eventOffset; iEvent++){
        if(find(eventsToSkip.begin(), eventsToSkip.end(), iEvent) != eventsToSkip.end()){
          Log(2)<<"\n\n=================================================================\n";
          Log(2)<<"Skipping event "<<iEvent<<"\n";
          continue;
        }
        
        events.LoadEventsFromFiles(dataType, iSig, cutLevel, iEvent);
        auto event = events.At(dataType, iSig, year, iEvent-eventOffset);
        
        if(!event->HasFriendData()){
          Log(2)<<"Warning -- skipping event "<<iEvent<<" as it has no friend info\n";
          event->SetWasTagged(false);
          continue;
        }
        
        cout<<"\n\n=================================================================\n";
        cout<<"helixTagger -- processing event "<<iEvent<<"\n";
        cout<<"cut level: "<<cutLevel<<"\n";
        
        for(auto &track : event->GetTracks()){
          int nTrackerLayers = -1;
          if(iSig == kTaggerBackgroundWithPU || iSig == kTaggerBackgroundNoPU){
            nTrackerLayers = RandInt(3, 6);
          }
          
          Helices fittedHelices = fitter->FitHelices(event->GetClusters(), *track, *event->GetVertex(), nTrackerLayers);
          
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
  }
  
  return 0;
}
