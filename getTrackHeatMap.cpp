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

void plotResults(map<string, TH2D*> &heatMaps, TH1D *invMass)
{
  TCanvas *c1 = new TCanvas("c1","c1",1000,1500);
  c1->Divide(2,3);
  
  c1->cd(1);
  invMass->GetXaxis()->SetTitle("m_{inv} (GeV)");
  invMass->GetYaxis()->SetTitle("# entries");
  invMass->Draw();
  drawZline();
  drawZline( 1, 2);
  drawZline(-1, 2);

  c1->cd(2); heatMaps["1sigma_all"]->Draw("colz");
  c1->cd(3); heatMaps["1sigma_2layers"]->Draw("colz");
  c1->cd(4); heatMaps["1sigma_3layers"]->Draw("colz");
  c1->cd(5); heatMaps["1sigma_4layers"]->Draw("colz");
  c1->cd(6); heatMaps["1sigma_4layers"]->Draw("colz");
  
  c1->Update();
}

EventSet getZmumuEvents()
{
  string initPrefix = getPathPrefix();
  
  EventSet events; events.LoadEventsFromFiles(initPrefix);

  EventCut eventCut; TrackCut trackCut; JetCut jetCut; LeptonCut leptonCut;
  CutsManager cutsManager;
  cutsManager.GetCuts(eventCut, trackCut, jetCut, leptonCut);
  
  // Pre-select Z->μμ events
  eventCut.SetNmuons(range<int>(2,2));
  eventCut.SetNleptons(range<int>(2,inf));
  eventCut.SetRequireMuonsFromZ(true);
  eventCut.SetRequireTwoOppositeMuons(true);
  
  events.ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  return events;
}

EventSet getSingleLeptonEvents()
{
  string initPrefix = getPathPrefix();
  
  EventSet events; events.LoadEventsFromFiles(initPrefix);

  EventCut eventCut; TrackCut trackCut; JetCut jetCut; LeptonCut leptonCut;
  CutsManager cutsManager;
  cutsManager.GetCuts(eventCut, trackCut, jetCut, leptonCut);
  
  // Pre-select W->μν events
  eventCut.SetNmuons(range<int>(1,1));
  eventCut.SetNleptons(range<int>(1,inf));
  
  events.ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  return events;
}

double getDimuonMass(const shared_ptr<Event> &event)
{
  if(event->GetNleptons() != 2) return -1;
  
  auto lepton1 = event->GetLepton(0);
  auto lepton2 = event->GetLepton(1);
  
  TLorentzVector muonVec1;
  muonVec1.SetPtEtaPhiM(lepton1->GetPt(), lepton1->GetEta(), lepton1->GetPhi(), muonMass);
  
  TLorentzVector muonVec2;
  muonVec2.SetPtEtaPhiM(lepton2->GetPt(), lepton2->GetEta(), lepton2->GetPhi(), muonMass);
  
  TLorentzVector dimuon = muonVec1 + muonVec2;
  
  return dimuon.M();
}

map<string, TH2D*> initHeatMaps()
{
  map<string, TH2D*> heatMaps;
  vector<string> heatMapsNames = {
    "1sigma_all", "1sigma_2layers", "1sigma_3layers", "1sigma_4layers", "1sigma_5layers", "1sigma_6layers",
    "2sigma_all", "2sigma_2layers", "2sigma_3layers", "2sigma_4layers", "2sigma_5layers", "2sigma_6layers",
    "3sigma_all", "3sigma_2layers", "3sigma_3layers", "3sigma_4layers", "3sigma_5layers", "3sigma_6layers",
    "4sigma_all", "4sigma_2layers", "4sigma_3layers", "4sigma_4layers", "4sigma_5layers", "4sigma_6layers",
  };
  
  for(string name : heatMapsNames){
    heatMaps[name] = new TH2D(name.c_str(), name.c_str(), 100, -TMath::Pi(), TMath::Pi(), 100, -2.4, 2.4);
    heatMaps[name]->GetXaxis()->SetTitle("#phi");
    heatMaps[name]->GetYaxis()->SetTitle("#eta");
  }
  return heatMaps;
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  config = ConfigManager(configPath);
  
  EventSet events = getZmumuEvents();
  
  cout<<"\nCreating inv mass plot"<<endl;
  
  TH1D *invMass = new TH1D("invMass #mu#mu", "invMass #mu#mu", 1000, 0, 500);
  
  
  map<string, TH2D*> heatMaps = initHeatMaps();
  
  EventCut eventCut; TrackCut trackCut; JetCut jetCut; LeptonCut leptonCut;
  trackCut.SetCaloEmEnergy(range<double>(0.0, 2.0));
  trackCut.SetNlayers(range<int>(2, 6));
  events.ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  vector<int> nEventsWithFakeTrack(5, 0);
  vector<int> nEvents(5, 0);
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    
    for(int iEvent=0; iEvent<events.size(xtracks::kData, kControlRegion, year); iEvent++){
      auto event = events.At(xtracks::kData, kControlRegion, year, iEvent);
      
      double mass = getDimuonMass(event);
      if(mass < 0) continue;
      
      invMass->Fill(mass);
      
      for(int i=0; i<5; i++){
        if(fabs(mass-Zmass) < i*Zwidth){
          if(event->GetNtracks() == 1) nEventsWithFakeTrack[i]++;
          nEvents[i]++;
        }
      }
      
      for(int iTrack=0; iTrack<event->GetNtracks(); iTrack++){
        auto track = event->GetTrack(iTrack);
        int nLayers = track->GetNtrackerLayers();
        
        for(int i=0; i<5; i++){
          if(fabs(mass-Zmass) < i*Zwidth){
            heatMaps[to_string(i)+"sigma_all"]->Fill(track->GetPhi(), track->GetEta());
            if(nLayers >= 2 && nLayers <= 6){
              string name = to_string(i)+"sigma_"+to_string(nLayers)+"layers";
              heatMaps[name]->Fill(track->GetPhi(), track->GetEta());
            }
          }
        }
      }
    }
  }
  plotResults(heatMaps, invMass);
  

  TFile *outFile = new TFile("results/tracksHeatMap.root", "recreate");
  outFile->cd();
  invMass->Write();
  for(auto &[name, hist] : heatMaps) hist->Write();
  outFile->Close();
  

  cout<<"Number of events with a fake track below 7 layers: "<<nEventsWithFakeTrack[4]<<endl;
  cout<<"Total number of events analyzed: "<<nEvents[4]<<endl;
  
  
  for(int i=1; i<5; i++){
    cout<<"2 fakes probability ("<<i<<"σ): "<<pow((double)nEventsWithFakeTrack[i]/nEvents[i], 2);
    cout<<" +/- "<<2*nEventsWithFakeTrack[i]/pow(nEvents[i],2)*sqrt(nEventsWithFakeTrack[i]+1/nEvents[i])<<endl;
  }
  
  
  EventSet singleLeptonEvents = getSingleLeptonEvents();
  
  int nSingleLeptonEvents =  singleLeptonEvents.size(kData, kControlRegion, 2018);
  
  trackCut.SetCaloEmEnergy(range<double>(0.0, 2.0));
  trackCut.SetNlayers(range<int>(2, 6));
  eventCut.SetNtracks(range<int>(2, 2));
  
  singleLeptonEvents.ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  int nTwoFakeTrackEvents =  singleLeptonEvents.size(kData, kControlRegion, 2018);
  
  cout<<"N single lepton events: "<<nSingleLeptonEvents<<endl;
  cout<<"N single lepton events with two fake tracks: "<<nTwoFakeTrackEvents<<endl;
  cout<<"Fraction: "<<(double)nTwoFakeTrackEvents/nSingleLeptonEvents<<endl;
  
  
  theApp.Run();
  return 0;
}





