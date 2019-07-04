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
  canvas->Divide(4,3);

  
  TH1D *nCommonPointsHistPion  = new TH1D("n rec clusters on helix pion", "n rec clusters on helix pion", 20, 0, 20);
  TH1D *nCommonPointsHistNoise = new TH1D("n rec clusters on helix noise", "n rec clusters on helix noise", 20, 0, 20);
  TH1D *nCommonPointsHistAll   = new TH1D("n rec clusters on helix all", "n rec clusters on helix all", 20, 0, 20);
  
  const int nThresholds = 10;
  
  int nEventsWithHelixPion = 0;
  int nEventsWithHelixNoise = 0;
  int nEventsWithHelixAll = 0;
  
  int nEventsWithTrueHelixPion[nThresholds] = {0};
  int nEventsWithTrueHelixNoise[nThresholds] = {0};
  
  int nEventsWithTrueAndFakeAll[nThresholds] = {0};
  int nEventsWithTrueOrFakeAll[nThresholds] = {0};
  int nEventsWithFakeOnlyAll[nThresholds] = {0};
  int nEventsWithTrueOnlyAll[nThresholds] = {0};
  int nEventsWithTrueAll[nThresholds] = {0};
  
  map<string, tuple<int, double, double>> histParams = {
    {"nMaxHits", {30, 0, 30}},
    {"nMaxLayers", {30, 0, 30}},
    {"nMaxLength", {30, 0, 6}},
    {"nAvgHits", {30, 0, 30}},
    {"nAvgLayers", {30, 0, 30}},
    {"nAvgLength", {30, 0, 6}},
    {"slope", {100, -1000, 1000}},
    {"slopeFactor", {100, -4000, 0}},
    {"radius", {100, 0, 1000}},
    {"radiusFactor", {100, 0, 2000}},
  };
  
  map<string, pair<TH1D*, TH1D*>> hists;
  
  for(auto &[title, params] : histParams){
    auto &[nBins, min, max] = params;
    hists[title] = make_pair(new TH1D((title+"Pion").c_str(), (title+"Pion").c_str(), nBins, min, max),
                             new TH1D((title+"Noise").c_str(), (title+"Noise").c_str(), nBins, min, max));
  }
  
  int nEventsAnalyzed = 0;
  
  EventSet events;
  events.LoadEventsFromFiles(cutLevel);
  int nEvents=events.size(dataType, setIter);
  
  auto start = now();

  for(auto iEvent=0; iEvent<nEvents; iEvent++){
    auto event = events.At(dataType, setIter, iEvent);
    
    if(iEvent%10==0)  cout<<"|";
    else              cout<<".";
    
    if(event->GetNtracks() != 1 || event->GetGenPionHelices().size() != 1) continue;
    
    nEventsAnalyzed++;
    
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
    
    size_t maxNclusters = 0;

    for(auto helix : fittedHelicesPion){
      size_t nPionClusters = helixProcessor.GetNcommonPoints(pionHelix, helix);
      if(nPionClusters > maxNclusters) maxNclusters = nPionClusters;
    }
    
    if(fittedHelicesPion.size() > 0){
      nEventsWithHelixPion++;
      
      nCommonPointsHistPion->Fill(maxNclusters);
      hists["nMaxHits"].first->Fill(helixProcessor.GetMaxNhits(fittedHelicesPion));
      hists["nMaxLayers"].first->Fill(helixProcessor.GetMaxNlayers(fittedHelicesPion));
      hists["nMaxLength"].first->Fill(helixProcessor.GetMaxLength(fittedHelicesPion));
      hists["nAvgHits"].first->Fill(helixProcessor.GetAvgNhits(fittedHelicesPion));
      hists["nAvgLayers"].first->Fill(helixProcessor.GetAvgNlayers(fittedHelicesPion));
      hists["nAvgLength"].first->Fill(helixProcessor.GetAvgLength(fittedHelicesPion));
      
      for(auto &helix : fittedHelicesPion){
        hists["slope"].first->Fill(helix.GetSlope(helix.GetTmin()));
        hists["slopeFactor"].first->Fill(helix.GetSlopeFactor());
        hists["radius"].first->Fill(helix.GetRadius(helix.GetTmin()));
        hists["radiusFactor"].first->Fill(helix.GetRadiusFactor());
      }
      
      for(int iThreshold=0; iThreshold<nThresholds; iThreshold++){
        if(maxNclusters >= iThreshold) nEventsWithTrueHelixPion[iThreshold]++;
      }
    }
    
    maxNclusters=0;
    for(auto helix : fittedHelicesNoise){
      size_t nPionClusters = helixProcessor.GetNcommonPoints(pionHelix, helix);
      if(nPionClusters > maxNclusters) maxNclusters = nPionClusters;
    }

    if(fittedHelicesNoise.size() > 0){
      nEventsWithHelixNoise++;
      
      nCommonPointsHistNoise->Fill(maxNclusters);
      hists["nMaxHits"].second->Fill(helixProcessor.GetMaxNhits(fittedHelicesNoise));
      hists["nMaxLayers"].second->Fill(helixProcessor.GetMaxNlayers(fittedHelicesNoise));
      hists["nMaxLength"].second->Fill(helixProcessor.GetMaxLength(fittedHelicesNoise));
      hists["nAvgHits"].second->Fill(helixProcessor.GetAvgNhits(fittedHelicesNoise));
      hists["nAvgLayers"].second->Fill(helixProcessor.GetAvgNlayers(fittedHelicesNoise));
      hists["nAvgLength"].second->Fill(helixProcessor.GetAvgLength(fittedHelicesNoise));
      
      for(auto &helix : fittedHelicesNoise){
        hists["slope"].second->Fill(helix.GetSlope(helix.GetTmin()));
        hists["slopeFactor"].second->Fill(helix.GetSlopeFactor());
        hists["radius"].second->Fill(helix.GetRadius(helix.GetTmin()));
        hists["radiusFactor"].second->Fill(helix.GetRadiusFactor());
      }
      
      for(int iThreshold=0; iThreshold<nThresholds; iThreshold++){
        if(maxNclusters >= iThreshold) nEventsWithTrueHelixNoise[iThreshold]++;
      }
    }
    
    maxNclusters=0;
    int nTrueHelices[nThresholds]={0};
    int nFakeHelices[nThresholds]={0};
    for(auto helix : fittedHelicesAll){
      size_t nPionClusters = helixProcessor.GetNcommonPoints(pionHelix, helix);
      if(nPionClusters > maxNclusters) maxNclusters = nPionClusters;
      
      for(int iThreshold=0; iThreshold<nThresholds; iThreshold++){
        if(nPionClusters >= iThreshold) nTrueHelices[iThreshold]++;
        else                            nFakeHelices[iThreshold]++;
      }
    }
    
    if(fittedHelicesAll.size() > 0){
      nCommonPointsHistAll->Fill(maxNclusters);
      
      nEventsWithHelixAll++;
      
      for(int iThreshold=0; iThreshold<nThresholds; iThreshold++){
        if(nTrueHelices[iThreshold]!=0) nEventsWithTrueAll[iThreshold]++;
        if(nTrueHelices[iThreshold]!=0 && nFakeHelices[iThreshold]!=0) nEventsWithTrueAndFakeAll[iThreshold]++;
        if(nTrueHelices[iThreshold]!=0 || nFakeHelices[iThreshold]!=0) nEventsWithTrueOrFakeAll[iThreshold]++;
        
        if(nTrueHelices[iThreshold]>0 && nFakeHelices[iThreshold]==0) nEventsWithTrueOnlyAll[iThreshold]++;
        if(nTrueHelices[iThreshold]==0 && nFakeHelices[iThreshold]!=0) nEventsWithFakeOnlyAll[iThreshold]++;
      }
    }
  }
  cout<<endl;
  
  cout<<"Time: "<<duration(start, now())<<endl;
  
  cout<<"events with helix (pion): "<<100*nEventsWithHelixPion/(double)nEventsAnalyzed<<" %"<<endl;
  cout<<"events with helix (Noise): "<<100*nEventsWithHelixNoise/(double)nEventsAnalyzed<<" %"<<endl;
  cout<<"events with helix (All): "<<100*nEventsWithHelixAll/(double)nEventsAnalyzed<<" %"<<endl;
  
  for(int iThreshold=0; iThreshold<nThresholds; iThreshold++){
    cout<<"\n\n==================================="<<endl;
    cout<<"Results for threshold: "<<iThreshold<<endl;
    
    double c_eff_pure  = nEventsWithTrueHelixPion[iThreshold]/(double)nEventsAnalyzed;
    double c_fake_pure = (nEventsWithHelixNoise-nEventsWithTrueHelixNoise[iThreshold])/(double)nEventsAnalyzed;
    double c_eff_and   = nEventsWithTrueAndFakeAll[iThreshold]/(double)nEventsAnalyzed;
    double c_eff_or    = nEventsWithTrueOrFakeAll[iThreshold]/(double)nEventsAnalyzed;
    double c_fake      = nEventsWithFakeOnlyAll[iThreshold]/(double)nEventsAnalyzed;
    double c_eff_full  = nEventsWithTrueOnlyAll[iThreshold]/(double)nEventsAnalyzed;
    double c_eff       = nEventsWithTrueAll[iThreshold]/(double)nEventsAnalyzed;
    
    cout.imbue(locale("de_DE"));
    cout<<"\n"<<endl;
    cout<<"c_eff^pure: "<<c_eff_pure<<endl;
    cout<<"c_fake^pure: "<<c_fake_pure<<endl;
    cout<<"c_eff^and: "<<c_eff_and<<endl;
    cout<<"c_eff^or: "<<c_eff_or<<endl;
    cout<<"c_fake: "<<c_fake<<endl;
    cout<<"c_eff^full: "<<c_eff_full<<endl;
    cout<<"c_eff: "<<c_eff<<endl;
    cout<<setprecision(3);
    cout<<c_eff_pure<<" "<<c_fake_pure<<" "<<c_eff<<" "<<c_fake<<" "<<c_eff_full<<endl;
  }
  
  cout<<"Pure ROC curve points:"<<endl;
  for(int iThreshold=0; iThreshold<nThresholds; iThreshold++){
    double c_fake_pure = (nEventsWithHelixNoise-nEventsWithTrueHelixNoise[iThreshold])/(double)nEventsAnalyzed;
    double c_eff_pure  = nEventsWithTrueHelixPion[iThreshold]/(double)nEventsAnalyzed;
    
    cout<<c_fake_pure<<"\t"<<c_eff_pure<<endl;
  }
  
  cout<<"ROC curve points:"<<endl;
   for(int iThreshold=0; iThreshold<nThresholds; iThreshold++){
     double c_fake      = nEventsWithFakeOnlyAll[iThreshold]/(double)nEventsAnalyzed;
     double c_eff       = nEventsWithTrueAll[iThreshold]/(double)nEventsAnalyzed;
     cout<<c_fake<<"\t"<<c_eff<<endl;
   }
  
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
  
  int iPad=2;
  TFile *outFile = new TFile("results/analyzeTagger.root","recreate");
  
  for(auto &[title, histPair] : hists){
    canvas->cd(iPad++);
    histPair.first->SetLineColor(kGreen+2);
    histPair.first->SetFillColorAlpha(kGreen+2, 0.3);
    histPair.first->DrawNormalized(iPad==3 ? "" : "same");
    histPair.second->SetLineColor(kRed);
    histPair.second->SetFillColorAlpha(kRed, 0.3);
    histPair.second->DrawNormalized("same");
    
    outFile->cd();
    histPair.first->Write();
    histPair.second->Write();
  }
  
  double legendW=0.25, legendH=0.40, legendX=0.65, legendY=0.1;
  TLegend *leg = new TLegend(legendX,legendY,legendX+legendW,legendY+legendH);
  leg->AddEntry(nCommonPointsHistPion, "Pion hits only", "elp");
  leg->AddEntry(nCommonPointsHistNoise,"Noise hits only", "elp");
  leg->AddEntry(nCommonPointsHistAll, "Full event", "elp");
  
  leg->Draw();
  
  outFile->cd();
  nCommonPointsHistPion->Write();
  nCommonPointsHistNoise->Write();
  nCommonPointsHistAll->Write();
  outFile->Close();
  
  canvas->Update();
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
