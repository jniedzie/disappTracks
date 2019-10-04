//  getTrackHeatMap.cpp
//
//  Created by Jeremi Niedziela on 01/10/2019.

#include "Event.hpp"
#include "EventSet.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"
#include "CutsManager.hpp"

const double Zmass = 91.1876; // GeV
const double Zwidth = 2.4952; // GeV

const double muonMass = 0.1056583745; // GeV
const double pionMass = 0.13957061; // GeV

string configPath = "configs/analysis.md";


void saveEvents(const EventSet &events)
{
  if(!config.params["save_events"]) return;
  string prefix = "after_L"+to_string((int)config.params["cuts_level"])+"/";
  events.SaveEventsToFiles(prefix);
}

/// Returns path prefix for cuts level and category selected in the config file
string getPathPrefix()
{
  string prefix;
  if(config.params["cuts_level"]==0) prefix = "";
  if(config.params["cuts_level"]==1) prefix = "after_L0/";
  return prefix;
}

void drawZline(int shift=0, int style=1)
{
  TLine *zPeakLine = new TLine(Zmass+shift*Zwidth, 0, Zmass+shift*Zwidth, 1E4);
  zPeakLine->SetLineColor(kCyan);
  zPeakLine->SetLineWidth(1);
  zPeakLine->SetLineStyle(style);
  zPeakLine->Draw();
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  config = ConfigManager(configPath);
  
  string initPrefix = getPathPrefix();

  EventSet events;
  events.LoadEventsFromFiles(initPrefix);
  cout<<"\n\nInitial yields"<<endl;
  events.PrintYields();
  
  if(config.params["cuts_level"]==0){
    EventCut eventCut; TrackCut trackCut; JetCut jetCut; LeptonCut leptonCut;
    eventCut.SetNmuons(range<int>(2,2));
    eventCut.SetNtracks(range<int>(1, inf));
    eventCut.SetRequireMuonsFromZ(true);
    eventCut.SetRequireTwoOppositeMuons(true);
    
    trackCut.SetNmissingInnerPixel(range<int>(0, 0));
    trackCut.SetNmissingMiddleTracker(range<int>(0, 0));
    trackCut.SetRelativeIsolation(range<double>(0.0, 0.5));
    trackCut.SetNlayers(range<int>(2, inf));
    trackCut.SetEta(range<double>(-2.1, 2.1));
    
    events.ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
    cout<<"\n\nEvents with >0 muon and >0 track (3 or 4 layers):"<<endl;
    events.PrintYields();
    saveEvents(events);
    return 0;
  }
  
  
  cout<<"\nCreating inv mass plot"<<endl;
  
  
  TH1D *invMass = new TH1D("invMass #mu#mu", "invMass #mu#mu", 1000, 0, 500);
  invMass->GetXaxis()->SetTitle("m_{inv} (GeV)");
  invMass->GetYaxis()->SetTitle("# entries");
  
  map<string, TH2D*> heatMaps;
  vector<string> heatMapsNames = {
    "1sigma_all", "1sigma_2layers", "1sigma_3layers", "1sigma_4layers", "1sigma_5layers", "1sigma_6layers",
    "2sigma_all", "2sigma_2layers", "2sigma_3layers", "2sigma_4layers", "2sigma_5layers", "2sigma_6layers",
    "3sigma_all", "3sigma_2layers", "3sigma_3layers", "3sigma_4layers", "3sigma_5layers", "3sigma_6layers",
  };
  
  for(string name : heatMapsNames){
    heatMaps[name] = new TH2D(name.c_str(), name.c_str(), 100, -TMath::Pi(), TMath::Pi(), 100, -2.4, 2.4);
    heatMaps[name]->GetXaxis()->SetTitle("#phi");
    heatMaps[name]->GetYaxis()->SetTitle("#eta");
  }
  
  EventCut eventCut; TrackCut trackCut; JetCut jetCut; LeptonCut leptonCut;
  
  events.ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  int nEventsWithFakeTrack = 0;
  int nEvents = 0;
  
  for(int iEvent=0; iEvent<events.size(xtracks::kData, kMET_Run2018A_CR); iEvent++){
    auto event = events.At(xtracks::kData, kMET_Run2018A_CR, iEvent);
    
    if(event->GetNleptons() != 2) continue;
    
    auto lepton1 = event->GetLepton(0);
    auto lepton2 = event->GetLepton(1);
    
    TLorentzVector muonVec1;
    muonVec1.SetPtEtaPhiM(lepton1->GetPt(), lepton1->GetEta(), lepton1->GetPhi(), muonMass);
    
    TLorentzVector muonVec2;
    muonVec2.SetPtEtaPhiM(lepton2->GetPt(), lepton2->GetEta(), lepton2->GetPhi(), muonMass);
    
    TLorentzVector dimuon = muonVec1 + muonVec2;
    
    double mass = dimuon.M();
    invMass->Fill(mass);
    
    for(int iTrack=0; iTrack<event->GetNtracks(); iTrack++){
      auto track = event->GetTrack(iTrack);
      int nLayers = track->GetNtrackerLayers();
      
      if(nLayers < 7) nEventsWithFakeTrack++;
      
      if(fabs(mass-Zmass) < Zwidth){
        heatMaps["1sigma_all"]->Fill(track->GetPhi(), track->GetEta());
        if(nLayers >= 2 && nLayers <= 6){
          string name = "1sigma_"+to_string(nLayers)+"layers";
          heatMaps[name]->Fill(track->GetPhi(), track->GetEta());
        }
      }
      if(fabs(mass-Zmass) < 2*Zwidth){
        heatMaps["2sigma_all"]->Fill(track->GetPhi(), track->GetEta());
        if(nLayers >= 2 && nLayers <= 6){
          string name = "2sigma_"+to_string(nLayers)+"layers";
          heatMaps[name]->Fill(track->GetPhi(), track->GetEta());
        }
      }
      if(fabs(mass-Zmass) < 3*Zwidth){
        heatMaps["3sigma_all"]->Fill(track->GetPhi(), track->GetEta());
        if(nLayers >= 2 && nLayers <= 6){
          string name = "3sigma_"+to_string(nLayers)+"layers";
          heatMaps[name]->Fill(track->GetPhi(), track->GetEta());
        }
      }
    }
    nEvents++;
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",1000,1500);
  c1->Divide(2,3);
  
  c1->cd(1);
  invMass->Draw();
  drawZline();
  drawZline( 1, 2);
  drawZline(-1, 2);

  c1->cd(2);
  heatMaps["1sigma_all"]->Draw("colz");
  
  c1->cd(3);
  heatMaps["1sigma_2layers"]->Draw("colz");
  
  c1->cd(4);
  heatMaps["1sigma_3layers"]->Draw("colz");
  
  c1->cd(5);
  heatMaps["1sigma_4layers"]->Draw("colz");
  
  c1->cd(6);
  heatMaps["1sigma_4layers"]->Draw("colz");
  
  c1->Update();

  TFile *outFile = new TFile("results/tracksHeatMap.root", "recreate");
  outFile->cd();
  invMass->Write();
  for(auto &[name, hist] : heatMaps) hist->Write();
  outFile->Close();
  

  cout<<"Number of events with a fake tracks below 7 layers: "<<nEventsWithFakeTrack<<endl;
  cout<<"Total number of events analyzed: "<<nEvents<<endl;
  cout<<"Fake probability: "<<(double)nEventsWithFakeTrack/nEvents<<endl;
  cout<<"2 fakes probability: "<<pow((double)nEventsWithFakeTrack/nEvents, 2)<<endl;
  
  
  theApp.Run();
  return 0;
}





