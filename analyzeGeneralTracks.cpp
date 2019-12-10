//  analyzeGeneralTracks.cpp
//  Created by Jeremi Niedziela on 06/12/2019.

#include "EventSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"

string configPath = "configs/analysis.md";

/// Returns path prefix for cuts level and category selected in the config file
string getPathPrefix()
{
  string prefix = "";
   
  if(config.secondaryCategory == "Zmumu")   prefix += "Zmumu/";
  if(config.secondaryCategory == "Wmunu")   prefix += "Wmunu/";
  if(config.secondaryCategory == "LowMET")   prefix += "LowMET/";
  
  if(config.params["cuts_level"]==0) prefix += "after_L0/";
  if(config.params["cuts_level"]==1) prefix += "after_L1/"+config.category+"/";
  if(config.params["cuts_level"]==2) prefix += "after_L1/"+config.category+"/afterHelixTagging/";
  
  return prefix;
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = ConfigManager(configPath);
  
  EventSet events; events.LoadEventsFromFiles(getPathPrefix());
  
   for(int year : years){
     if(!config.params["load_"+to_string(year)]) continue;
   
     for(ESignal iSig : signals){
       if(!config.runSignal[iSig]) continue;
       
       for(int iEvent=0; iEvent<events.size(kSignal, iSig, year); iEvent++){
         auto event = events.At(kSignal, iSig, year, iEvent);
         cout<<"iEvent: "<<iEvent<<endl;
         
         if(!event){
           cout<<"Event not found"<<endl;
           exit(0);
         }
         
         for(auto track : event->GetTracks()){
         
           for(auto generalTrack : event->GetGeneralTracks()){
             if(fabs(track->GetEta()-generalTrack.GetEta()) < 0.5 &&
                fabs(track->GetPhi()-generalTrack.GetPhi()) < 0.5){
               cout<<"Close track "<<generalTrack.GetD0()<<endl;
             }
           }
         
           
         }
       }
     }
   }
  
  
//  theApp.Run();
  return 0;
}
