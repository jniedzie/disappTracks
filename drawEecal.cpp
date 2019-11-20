//  drawEecal.cpp
//  Created by Jeremi Niedziela on 20/11/2019.

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
  events.LoadEventsFromFiles();
  
  EventCut eventCut; TrackCut trackCut; JetCut jetCut; LeptonCut leptonCut;
  
  // Remove bad tracks
  trackCut.SetNmissingInnerPixel(range<int>(0, 0));
  trackCut.SetNmissingMiddleTracker(range<int>(0, 0));
  trackCut.SetRelativeIsolation(range<double>(0.0, 0.5));
  trackCut.SetNlayers(range<int>(2, inf));
  trackCut.SetEta(range<double>(-2.1, 2.1));
  
  // Check MET properties
  eventCut.SetMetNoMuPt(range<double>(200,inf));
  eventCut.SetRequireMetNoMuTrigger(true);
  eventCut.SetRequirePassingAllFilters(true);
  eventCut.SetJetMetDeltaPhi(range<double>(0.5,inf));
  
  // Check jet properties
  jetCut.SetPt(range<double>(30.0, inf));
  jetCut.SetChargedHadronEnergyFraction(range<double>(0.01,0.99));
  jetCut.SetNeutralHadronEnergyFraction(range<double>(0.01,0.99));
  eventCut.SetLeadingJetPt(range<double>(100,inf));
  eventCut.SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCut.SetLeadingJetNeHEF(range<double>(0,0.8));
  eventCut.SetLeadingJetChHEF(range<double>(0.1,1.0));
  
  // Check number of objects after cuts
  eventCut.SetNtracks(range<int>(1,inf));
  eventCut.SetNjets(range<int>(1,inf));
  eventCut.SetNtaus(range<int>(0,0));
  
  // Pre-select Z->μμ events
  eventCut.SetNmuons(range<int>(2,2));
  eventCut.SetNleptons(range<int>(2,inf));
  eventCut.SetRequireMuonsFromZ(true);
  eventCut.SetRequireTwoOppositeMuons(true);
  
  events.ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  TH1D *histEecalMC   = new TH1D("histEecalMC"  , "histEecalMC"   , 1000, 0, 100);
  TH1D *histEecalData = new TH1D("histEecalData", "histEecalData" , 1000, 0, 100);
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    for(auto iBck : backgrounds){
      if(!config.runBackground[iBck]) continue;
      
      for(int iEvent=0; iEvent<events.size(kBackground, iBck, year); iEvent++){
        auto event = events.At(kBackground, iBck, year, iEvent);
        
        for(auto &track : event->GetTracks()){
          histEecalMC->Fill(track->GetCaloEmEnergy());
        }
      }
    }
    
    for(auto iData : datas){
      if(!config.runData[iData]) continue;
      
      for(int iEvent=0; iEvent<events.size(kData, iData, year); iEvent++){
        auto event = events.At(kData, iData, year, iEvent);
        
        for(auto &track : event->GetTracks()){
          histEecalData->Fill(track->GetCaloEmEnergy());
        }
      }
    }
    
  }
  
  
  histEecalMC->SetLineColor(kBlue);
  histEecalMC->SetFillColorAlpha(kBlue, 0.3);
  
  histEecalData->SetLineColor(kRed);
  histEecalData->SetMarkerColor(kRed);
  histEecalData->SetMarkerStyle(20);
  histEecalData->SetMarkerSize(1.0);
  
  TCanvas *c1 = new TCanvas("c1","c1",800, 600);
  c1->cd();
  
  histEecalMC->DrawNormalized();
  histEecalData->DrawNormalized("PEsame");
  
  c1->Update();
  
  TFile *outFile = new TFile("results/ecalEnergy.root", "recreate");
  outFile->cd();
  histEecalMC->Write();
  histEecalData->Write();
  outFile->Close();
  
  theApp.Run();
  return 0;
}

