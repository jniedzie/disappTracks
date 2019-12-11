//  analyzeGeneralTracks.cpp
//  Created by Jeremi Niedziela on 06/12/2019.

#include "EventSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"

string configPath = "configs/analysis.md";

TH1D *generalPionHitsFraction, *taggerPionHitsFraction;

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

void fillHists(const EventSet &events)
{
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
         
         for(auto generalTrack : event->GetGeneralTracks()){
           double trueHitsFraction = generalTrack.GetNrecPionHits()/(double)(generalTrack.GetNrecHits());
           generalPionHitsFraction->Fill(trueHitsFraction);
         }
         
         for(auto taggerHelix : event->GetHelices()){
           double trueHitsFraction = taggerHelix.GetNrecPionHits()/(double)(taggerHelix.GetNrecHits()-1);
           taggerPionHitsFraction->Fill(trueHitsFraction);
         }
       }
     }
   }
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = ConfigManager(configPath);
  
  EventSet events; events.LoadEventsFromFiles(getPathPrefix());
  
  generalPionHitsFraction = new TH1D("generalNpionHits", "generalNpionHits", 100, 0, 1);
  taggerPionHitsFraction  = new TH1D("taggerNpionHits" , "taggerNpionHits" , 100, 0, 1);
  
  fillHists(events);
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000);
  c1->Divide(2,2);
  
  c1->cd(1); generalPionHitsFraction->Draw();
  c1->cd(2); taggerPionHitsFraction->Draw();
  
  
  
  c1->Update();
  theApp.Run();
  return 0;
}
