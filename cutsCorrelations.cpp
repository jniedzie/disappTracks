//
//  cutsCorrelations.cpp
//
//  Created by Jeremi Niedziela on 22/11/2018.
//

#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "CutsManager.hpp"
#include "Logger.hpp"

/// Creates a map of 2D histograms showing correlations between different variables
map<string, TH2D*> CreateCorrelationHistograms()
{
  //  name          nBinsX minX  maxX    nBinsY minY  maxY
  map<string, tuple<int, double, double, int, double, double>> histParams = {
    {"iso_vs_dedx"                , { 100  , 0 , 20   , 100 , 0   , 0.1   }},
    {"met_vs_dedx"                , { 100  , 0 , 20   , 100 , 200 , 1200  }},
    {"trackPt_vs_missing"         , { 100  , 0 , 1000 , 20  , 0   , 20    }},
    {"deltaJetTrack_vs_missing"   , { 20   , 0 , 2    , 20  , 0   , 20    }},
  };
  
  map<string, tuple<string, string>> histAxesNames = {
    {"iso_vs_dedx"                , { "Log-likelihood per dE/dx hit (MeV/cm)" , "Relative isolation" }},
    {"met_vs_dedx"                , { "Log-likelihood per dE/dx hit (MeV/cm)" , "MET p_{T} (GeV)"    }},
    {"trackPt_vs_missing"         , { "Missing outer tracker hits"            , "Track p_{T}"        }},
    {"deltaJetTrack_vs_missing"   , { "Missing outer tracker hits"            , "#Delta R(jet,track)"}},
  };
  
  map<string, TH2D*> correlationHists;
  
  for(auto &[name, params] : histParams){
    auto &[nBinsX, minX, maxX, nBinsY, minY, maxY] = params;
    correlationHists[name] = new TH2D(name.c_str(), name.c_str(), nBinsX, minX, maxX, nBinsY, minY, maxY);
    auto &[titleX, titleY] = histAxesNames[name];
    correlationHists[name]->GetXaxis()->SetTitle(titleX.c_str());
    correlationHists[name]->GetYaxis()->SetTitle(titleY.c_str());
    correlationHists[name]->SetTitle("");
  }
  return correlationHists;
}

/// Fills correlation histograms from background events
void FillCorrelationHistograms(map<string, TH2D*> &correlationHists, const EventSet &events)
{
//  for(int iBck=0; iBck<kNbackgrounds; iBck++){
//    for(int iEvent=0;iEvent<events.size(xtracks::kBackground, iBck);iEvent++){
//      auto event = events.At(xtracks::kBackground, iBck, iEvent);
  
  for(int iEvent=0;iEvent<events.size(xtracks::kSignal, kWino_M_650_cTau_10);iEvent++){
      auto event = events.At(xtracks::kSignal, kWino_M_650_cTau_10, iEvent);
  
      
      for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
        auto track = event->GetTrack(iTrack);
        
        double avgDedx        = track->GetAverageDedx();
        double minDedx        = track->GetMinDedx();
        double likelihoodDedx = track->GetDedxLikelihood();
        
        correlationHists["iso_vs_dedx"]->Fill(likelihoodDedx, track->GetRelativeIsolation());
        correlationHists["met_vs_dedx"]->Fill(likelihoodDedx, event->GetMetNoMuPt());
        correlationHists["trackPt_vs_missing"]->Fill(track->GetPt(), track->GetNmissingOuterTrackerHits());
        
        for(int iJet=0;iJet<event->GetNjets();iJet++){
          auto jet = event->GetJet(iJet);
          
          double deltaR = sqrt(pow(track->GetPhi() - jet->GetPhi(),2)+pow(track->GetEta() - jet->GetEta(),2));
          correlationHists["deltaJetTrack_vs_missing"]->Fill(deltaR, track->GetNmissingOuterTrackerHits());
        }
      }
    }
//  }
}

/// Draws correlation histograms in a canvas
void DrawAndSaveCorrelationHistograms(const map<string, TH2D*> &correlationHists)
{
  TCanvas *c1 = new TCanvas("c1","c1",1000,1500);
  if(correlationHists.size() <= 4) c1->Divide(2,2);
  else if(correlationHists.size() <= 6) c1->Divide(2,3);
  else if(correlationHists.size() <= 9) c1->Divide(3,3);
  int iPad=1;
  for(auto &[name, hist] : correlationHists){
    c1->cd(iPad++);
    hist->Draw("colz");
    
    TCanvas *c2 = new TCanvas("c2", "c2", 800, 600);
    c2->cd();
    gPad->SetLeftMargin(0.12);
    hist->Draw("colz");
    c2->SaveAs(("plots/"+name+".pdf").c_str());
  }
  c1->Update();
}

/// Starting point of the application
int main(int argc, char* argv[])
{
  TApplication *theApp = new TApplication("App", &argc, argv);
  config = ConfigManager("configs/analysis.md");
  
  gStyle->SetOptStat(0);
  
  // All events with initial cuts only
  EventSet events;
  string prefix;
  if(config.params["cuts_level"]==0) prefix = "after_L0/";
  if(config.params["cuts_level"]==1) prefix = "after_L1/"+config.category+"/";
  events.LoadEventsFromFiles(prefix);
  
  // Draw correlation plots
  map<string, TH2D*> correlationHists = CreateCorrelationHistograms();
  FillCorrelationHistograms(correlationHists, events);
  DrawAndSaveCorrelationHistograms(correlationHists);
  
  theApp->Run();
  return 0;
}




