#include "Event.hpp"
#include "EventSet.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "CutsManager.hpp"

string configPath = "configs/analysis.md";

void saveEvents(const EventSet &events, const string &suffix)
{
  if(!config.params["save_events"]) return;
  
  string prefix = "";
  
  if(config.secondaryCategory == "Zmumu") prefix += "Zmumu/";
  if(config.secondaryCategory == "Wmunu") prefix += "Wmunu/";
  
  prefix += "after_L"+to_string((int)config.params["cuts_level"])+"/";
  prefix += suffix+"/";
  
  events.SaveEventsToFiles(prefix);
}

void plotEvents(const EventSet &events)
{
  if(config.params["draw_standard_plots"])  events.DrawStandardPlots();
  if(config.params["draw_per_layer_plots"])  events.DrawPerLayerPlots();
}

void printDetails(const EventSet &events)
{
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    if(config.params["print_background_details"]){
      for(EBackground iBck : backgrounds){
        if(!config.runBackground[iBck]) continue;
        
        cout<<"Background events in "<<backgroundTitle.at(iBck)<<":"<<endl;
        for(int iEvent=0; iEvent<events.size(xtracks::kBackground, iBck, year); iEvent++){
          events.At(xtracks::kBackground, iBck, year, iEvent)->Print();
        }
      }
    }
    if(config.params["print_data_details"]){
      for(EData iData : datas){
        if(!config.runData[iData]) continue;
        
        cout<<"Data events in "<<dataTitle.at(iData)<<":"<<endl;
        for(int iEvent=0; iEvent<events.size(xtracks::kData, iData, year); iEvent++){
          events.At(xtracks::kData, iData, year, iEvent)->Print();
        }
      }
    }
    if(config.params["print_signal_details"]){
      for(ESignal iSig : signals){
        if(!config.runSignal[iSig]) continue;
        
        cout<<"Signal events in "<<signalTitle.at(iSig)<<":"<<endl;
        for(int iEvent=0;iEvent<events.size(xtracks::kSignal, iSig, year); iEvent++){
          events.At(xtracks::kSignal, iSig, year, iEvent)->Print();
        }
      }
    }
  }
}

void saveSurvivingEvents(const EventSet &events)
{
  bool updateData = false;
  
  for(EData iData : datas){
    if(config.runData[iData]){ updateData = true; break; }
  }
  
  string suffix = "L"+to_string_with_precision(config.params["cuts_level"],0);
  if(config.params["cuts_level"]==1) suffix+= "_"+config.category;
  
  if(updateData){
    ofstream dataSurvivingFile;
    dataSurvivingFile.open("results/survivingDataEventsAfter"+suffix+".txt");
    
    ofstream dataSurvivingFileByLS;
    dataSurvivingFileByLS.open("results/survivingDataEventsAfter"+suffix+"_byLS.txt");
    
    set<int> lumiSections;
    
    for(EData iData : datas){
      if(!config.runData[iData]) continue;
      for(int year : years){
        if(!config.params["load_"+to_string(year)]) continue;
        
        cout<<"Data events surviving cuts in "<<dataTitle.at(iData)<<":"<<events.size(xtracks::kData, iData, year)<<endl;
        for(int iEvent=0; iEvent<events.size(xtracks::kData, iData, year); iEvent++){
          auto event = events.At(xtracks::kData, iData, year, iEvent);
          int runNumber = event->GetRunNumber();
          int lumiSection = event->GetLumiSection();
          long long int eventNumber = event->GetEventNumber();
          
          dataSurvivingFile<<runNumber<<":"<<lumiSection<<":"<<eventNumber<<"\n";
          
          if(lumiSections.find(lumiSection) == lumiSections.end()){
            dataSurvivingFileByLS<<runNumber<<" "<<lumiSection<<"\n";
            lumiSections.insert(lumiSection);
          }
        }
      }
    }
    dataSurvivingFile.close();
  }
  
  bool updateSignals = false;
  
  for(ESignal iSig : signals){
    if(config.runSignal[iSig]){ updateSignals = true; break; }
  }
  
  if(updateSignals){
    ofstream signalSurvivingFile;
    signalSurvivingFile.open ("results/survivingSignalEventsAfter"+suffix+".txt");
    
    for(ESignal iSig : signals){
      if(!config.runSignal[iSig]) continue;
      for(int year : years){
        if(!config.params["load_"+to_string(year)]) continue;
        
        cout<<"Signal events surviving cuts in "<<signalTitle.at(iSig)<<":"<<events.size(xtracks::kSignal, iSig, year)<<endl;
        for(int iEvent=0; iEvent<events.size(xtracks::kSignal, iSig, year); iEvent++){
          auto event = events.At(xtracks::kSignal, iSig, year, iEvent);
          int runNumber = event->GetRunNumber();
          int lumiSection = event->GetLumiSection();
          long long int eventNumber = event->GetEventNumber();
          signalSurvivingFile<<runNumber<<":"<<lumiSection<<":"<<eventNumber<<"\n";
        }
      }
    }
    signalSurvivingFile.close();
  }
}

void processCuts(EventSet &events,
                 const EventCut &eventCut, const TrackCut &trackCut,
                 const JetCut &jetCut, const LeptonCut &leptonCut,
                 string suffix = "")
{
  events.ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  if(config.params["print_yields"]){
    cout<<"\n\nYields after level "<<config.params["cuts_level"]<<" cuts"<<endl;
    events.PrintYields();
  }
  saveEvents(events, suffix);
  plotEvents(events);
  printDetails(events);
  saveSurvivingEvents(events);
}

/// Returns path prefix for cuts level and category selected in the config file
string getPathPrefix()
{
  string prefix = "";
   
  if(config.secondaryCategory == "Zmumu") prefix += "Zmumu/";
  if(config.secondaryCategory == "Wmunu") prefix += "Wmunu/";
  
  if(config.params["cuts_level"]==0 || config.params["cuts_level"]==10) prefix = "";
  if(config.params["cuts_level"]==1) prefix += "after_L0/";
  if(config.params["cuts_level"]==2) prefix += "after_L1/"+config.category+"/";
  if(config.params["cuts_level"]==20) prefix += "afterHelixTagging/";
  
  return prefix;
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  config = ConfigManager(configPath);
  
  EventSet events; events.LoadEventsFromFiles(getPathPrefix());
  cout<<"\n\nInitial yields"<<endl; events.PrintYields();
  
  CutsManager cutsManager;
  EventCut eventCut; TrackCut trackCut; JetCut jetCut; LeptonCut leptonCut;
    
  if(config.secondaryCategory == "Zmumu")       cutsManager.GetZmumuCuts(eventCut, trackCut, jetCut, leptonCut);
  else if(config.secondaryCategory == "Wmunu")  cutsManager.GetWmunuCuts(eventCut, trackCut, jetCut, leptonCut);
  else                                          cutsManager.GetCuts(eventCut, trackCut, jetCut, leptonCut);
  
  if(config.params["cuts_level"] == 0){
    processCuts(events, eventCut, trackCut, jetCut, leptonCut);
  }
  else{
    processCuts(events, eventCut, trackCut, jetCut, leptonCut, config.category);
  }

  if(config.params["cuts_level"] == 20) events.DrawStandardPlots();
  cout<<"Done"<<endl;
  
  if(config.params["draw_standard_plots"]  || config.params["draw_per_layer_plots"]) theApp.Run();
  return 0;
}



