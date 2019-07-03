//  analyzeTagger.cpp
//
//  Created by Jeremi Niedziela on 25/01/2019.

#include "Helpers.hpp"
#include "Fitter.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "PerformanceMonitor.hpp"
#include "EventSet.hpp"

string configPath = "configs/eventDisplay.md";
string cutLevel = "after_L1/all/";//after_L1/";

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;

shared_ptr<Event> GetEvent(int iEvent);
vector<shared_ptr<Point>> GetClustersNoEndcaps(const shared_ptr<Event> &event, bool removePionClusters);

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  config = ConfigManager(configPath);
  auto fitter = make_unique<Fitter>();
  TCanvas *canvas = new TCanvas("Tagger analysis", "Tagger analysis", 2880,1800);
  canvas->Divide(2,2);

  
  TH1D *nCommonPointsHistPion  = new TH1D("n rec clusters on helix pion", "n rec clusters on helix pion", 20, 0, 20);
  TH1D *nCommonPointsHistNoise = new TH1D("n rec clusters on helix noise", "n rec clusters on helix noise", 20, 0, 20);
  TH1D *nCommonPointsHistAll   = new TH1D("n rec clusters on helix all", "n rec clusters on helix all", 20, 0, 20);
  
  EventSet events;
  events.LoadEventsFromFiles(cutLevel);
  int nEvents=events.size(dataType, setIter);
  
  auto start = now();

  for(auto iEvent=0; iEvent<nEvents; iEvent++){
    auto event = events.At(dataType, setIter, iEvent);
    
    cout<<"\n\n=================================================================\n"<<endl;
    cout<<"helixTagger -- processing event "<<iEvent<<endl;
    
    if(event->GetNtracks() != 1 || event->GetGenPionHelices().size() != 1) continue;
    
    auto track     = event->GetTracks().front();
    auto pionHelix = event->GetGenPionHelices().front();
    pionHelix.SetPoints(event->GetPionClusters());
    
    vector<Helix> fittedHelicesPion = fitter->FitHelices(event->GetPionClusters(), *track, *event->GetVertex());

    bool removePionClusters = true;
    auto pointsNoEndcapsNoPion = GetClustersNoEndcaps(event, removePionClusters);
    vector<Helix> fittedHelicesNoise = fitter->FitHelices(pointsNoEndcapsNoPion, *track, *event->GetVertex());
    
    removePionClusters = false;
    auto pointsNoEndcaps = GetClustersNoEndcaps(event, removePionClusters);
    vector<Helix> fittedHelicesAll = fitter->FitHelices(pointsNoEndcaps, *track, *event->GetVertex());
    
    for(auto helix : fittedHelicesPion){
      nCommonPointsHistPion->Fill(helixProcessor.GetNcommonPoints(pionHelix, helix));
    }
    for(auto helix : fittedHelicesNoise){
      nCommonPointsHistNoise->Fill(helixProcessor.GetNcommonPoints(pionHelix, helix));
    }
    for(auto helix : fittedHelicesAll){
      nCommonPointsHistAll->Fill(helixProcessor.GetNcommonPoints(pionHelix, helix));
    }
    
  }
  
  cout<<"Time: "<<duration(start, now());
  
  canvas->cd(1);
  nCommonPointsHistPion->SetLineColor(kGreen+2);
  nCommonPointsHistPion->SetFillColorAlpha(kGreen+2, 0.3);
  nCommonPointsHistPion->Draw();
  nCommonPointsHistNoise->SetLineColor(kRed);
  nCommonPointsHistNoise->SetFillColorAlpha(kRed, 0.3);
  nCommonPointsHistNoise->Draw("same");
  nCommonPointsHistAll->SetLineColor(kBlue);
  nCommonPointsHistAll->SetFillColorAlpha(kBlue, 0.3);
  nCommonPointsHistAll->Draw("same");
  
  double legendW=0.25, legendH=0.40, legendX=0.65, legendY=0.1;
  TLegend *leg = new TLegend(legendX,legendY,legendX+legendW,legendY+legendH);
  leg->AddEntry(nCommonPointsHistPion, "Pion hits only", "elp");
  leg->AddEntry(nCommonPointsHistNoise,"Noise hits only", "elp");
  leg->AddEntry(nCommonPointsHistAll, "Full event", "elp");
  
  leg->Draw();
  
  canvas->Update();
  
  TFile *outFile = new TFile("results/analyzeTagger.root","recreate");
  outFile->cd();
  nCommonPointsHistPion->Write();
  nCommonPointsHistNoise->Write();
  nCommonPointsHistAll->Write();
  
  outFile->Close();
  
  
  theApp.Run();
  return 0;
}

shared_ptr<Event> GetEvent(int iEvent)
{
  EventSet events;
  events.LoadEventFromFiles(dataType, setIter, iEvent, cutLevel);
  auto event = events.At(dataType, setIter, 0);
  
  if(!event){
    cout<<"helixTagger -- event not found"<<endl;
    exit(0);
  }
  
  return event;
}


vector<shared_ptr<Point>> GetClustersNoEndcaps(const shared_ptr<Event> &event, bool removePionClusters)
{
  vector<shared_ptr<Point>> pointsNoEndcaps;
  
  for(auto &point : event->GetTrackerClusters()){
    if(point->GetSubDetName() == "TID" || point->GetSubDetName() == "TEC" || point->GetSubDetName() == "P1PXEC") continue;
    
    if(point->GetSubDetName() != "TIB" && point->GetSubDetName() != "TOB" && point->GetSubDetName() != "P1PXB"){
      cout<<"Weird detector:"<<point->GetSubDetName()<<endl;
    }
    
    if(removePionClusters){
      bool isPionHit = false;
      for(auto &pionCluster : event->GetPionClusters()){
        if(*pionCluster == *point){
          isPionHit = true;
          break;
        }
      }
      if(isPionHit) continue;
    }
    
    pointsNoEndcaps.push_back(point);
  }
  
  return pointsNoEndcaps;
}
