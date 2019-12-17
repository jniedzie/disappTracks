//  analyzeGeneralTracks.cpp
//  Created by Jeremi Niedziela on 06/12/2019.

#include "EventSet.hpp"
#include "Helpers.hpp"
#include "ConfigManager.hpp"

string configPath = "configs/analysis.md";
string outPath = "results/plots_general_vs_tagger_test";

int nEvents;
int eventOffset;

enum ESource { kGeneral, kTagger };
constexpr initializer_list<ESource> sources = { kGeneral, kTagger };
map<ESource, int> colors = {{kGeneral, kRed}, {kTagger, kGreen}};

struct Hists {
  Hists(string title, int nBins, double min, double max){
    for(ESource source : sources){
      string fullTitle = (source == kTagger ? "tagger_" : "") + title;
      hists[source] = new TH1D(fullTitle.c_str(), fullTitle.c_str(), nBins, min, max);
      hists[source]->SetLineColor(colors[source]);
      hists[source]->SetMarkerColor(colors[source]);
      hists[source]->SetMarkerStyle(20);
      hists[source]->SetMarkerSize(1.0);
    }
  }
  
  void Draw(string options=""){
    hists[kGeneral]->Draw(options.c_str());
    hists[kTagger]->Draw(("same"+options).c_str());
  }
  
  void DrawNormalized(string options=""){
    hists[kGeneral]->DrawNormalized(options.c_str());
    hists[kTagger]->DrawNormalized(("same"+options).c_str());
  }

  void Write(){ hists[kGeneral]->Write(); hists[kTagger]->Write(); }
  
  map<ESource, TH1D*> hists;
};

struct Graphs {
  Graphs(string title){
    for(ESource source : sources){
      string fullTitle = (source == kTagger ? "tagger_" : "") + title;
      graphs[source] = new TGraphErrors();
      graphs[source]->SetTitle(fullTitle.c_str());
      graphs[source]->SetName(fullTitle.c_str());
      graphs[source]->SetLineColor(colors[source]);
      graphs[source]->SetMarkerColor(colors[source]);
      graphs[source]->SetMarkerStyle(20);
      graphs[source]->SetMarkerSize(1.0);
    }
  }
  
  void Draw(string options=""){
    graphs[kGeneral]->Draw(("A"+options).c_str());
    graphs[kTagger]->Draw(("same"+options).c_str());
  }

  void AddPoint(double x, double num, double den, ESource source){
    graphs[source]->SetPoint(iPoint[source], x, num/den);
    graphs[source]->SetPointError(iPoint[source]++, 0, num/den*sqrt(1/num + 1/den));
  }
  
  void AddPointROC(double fakeVal, double trueVal, double norm, ESource source){
    graphs[source]->SetPoint(iPoint[source]++, fakeVal/norm, trueVal/norm);
  }
  
  void Write(){ graphs[kGeneral]->Write(); graphs[kTagger]->Write(); }
  
  map<ESource, TGraphErrors*> graphs;
  map<ESource, int> iPoint;
};

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

struct Counter{
  void clearFound(){
    fakeFound[kGeneral].clear();
    fakeFound[kTagger].clear();
    trueFound[kGeneral].clear();
    trueFound[kTagger].clear();
  }
  map<ESource, map<double, int>> nFakes;
  map<ESource, map<double, int>> nTrue;
  map<ESource, map<double, bool>> fakeFound;
  map<ESource, map<double, bool>> trueFound;
};

void fillCounters(Counter &byFraction, Counter &byCompatibility, Counter &byHits,
                  double trackTruth, double trackCompatibility, int nPionHits,
                  ESource source)
{
  for(double fraction=0; fraction<1.1; fraction+=0.05){
    if(trackTruth >= fraction && !byFraction.trueFound[source][fraction]){
      byFraction.nTrue[source][fraction]++;
      byFraction.trueFound[source][fraction] = true;
    }
    if(trackTruth < fraction && !byFraction.fakeFound[source][fraction]){
      byFraction.nFakes[source][fraction]++;
      byFraction.fakeFound[source][fraction] = true;
    }
    
    if(trackCompatibility >= fraction && !byCompatibility.trueFound[source][fraction]){
      byCompatibility.nTrue[source][fraction]++;
      byCompatibility.trueFound[source][fraction] = true;
    }
    if(trackCompatibility < fraction && !byCompatibility.fakeFound[source][fraction]){
      byCompatibility.nFakes[source][fraction]++;
      byCompatibility.fakeFound[source][fraction] = true;
    }
  }
  for(double nHits=0; nHits<=21.0; nHits+=1.0){
    if(nPionHits >= nHits && !byHits.trueFound[source][nHits]){
      byHits.nTrue[source][nHits]++;
      byHits.trueFound[source][nHits] = true;
    }
    if(nPionHits < nHits && !byHits.fakeFound[source][nHits]){
      byHits.nFakes[source][nHits]++;
      byHits.fakeFound[source][nHits] = true;
    }
  }
}

void fillHists(Hists &trackTruthFraction, Hists &trackToPionSimilarity,
               Graphs &trueEventsByFraction,
               Graphs &trueEventsByHits,
               Graphs &trueEventsByCompatibility,
               Graphs &fakesEventsByFraction,
               Graphs &fakeEventsByCompatibility,
               Graphs &rocFraction,
               Graphs &rocCompatibility)
{
  Counter byFraction;
  Counter byHits;
  Counter byCompatibility;
  
  map<ESource, int> nAnalyzed;
  EventSet events; events.LoadEventsFromFiles(getPathPrefix());
  
  for(int year : years){
    if(!config.params["load_"+to_string(year)]) continue;
    cout<<"Year: "<<year<<endl;
    for(ESignal iSig : signals){
      if(!config.runSignal[iSig]) continue;
      cout<<"Sample: "<<iSig<<endl;
      cout<<"N events: "<<nEvents<<"\toffset: "<<eventOffset<<endl;
      
      for(int iEvent=0; iEvent<events.size(kSignal, iSig, year); iEvent++){
        auto event = events.At(kSignal, iSig, year, iEvent);
        
        if(!event){ cout<<"Event not found"<<endl; exit(0); }
        if(!event->WasTagged()) continue;
        
        byFraction.clearFound();
        byCompatibility.clearFound();
        byHits.clearFound();
        
        double nPionHits = event->GetPionClusters().size();
        if(nPionHits == 0) continue;
        nPionHits /= event->GetGenPionHelices().size();
        
        for(auto generalTrack : event->GetGeneralTracks()){
          if(!generalTrack.IsLooper()) continue;
          
          double trackTruth = generalTrack.GetNrecPionHits()/(double)(generalTrack.GetNrecHits());
          if(generalTrack.GetNrecHits() <= 0) trackTruth=0;
          trackTruthFraction.hists[kGeneral]->Fill(trackTruth);
          
          double trackPionCompatibility = generalTrack.GetNrecPionHits()/nPionHits;
          if(nPionHits <= 0) trackPionCompatibility=0;
          trackToPionSimilarity.hists[kGeneral]->Fill(trackPionCompatibility);
          
          fillCounters(byFraction, byCompatibility, byHits,
                       trackTruth, trackPionCompatibility, generalTrack.GetNrecPionHits(),
                       kGeneral);
        }
        
        for(auto taggerHelix : event->GetHelices()){
          
          double trackTruth = taggerHelix.GetNrecPionHits()/(double)(taggerHelix.GetNrecHits()-1);
          if(taggerHelix.GetNrecHits() <= 0) trackTruth=0;
          trackTruthFraction.hists[kTagger]->Fill(trackTruth);
          
          double trackPionCompatibility = taggerHelix.GetNrecPionHits()/nPionHits;
          if(nPionHits <= 0) trackPionCompatibility=0;
          trackToPionSimilarity.hists[kTagger]->Fill(trackPionCompatibility);
          
          fillCounters(byFraction, byCompatibility, byHits,
                       trackTruth, trackPionCompatibility, taggerHelix.GetNrecPionHits(),
                       kTagger);
        }
        nAnalyzed[kGeneral]++;
        if(event->GetHelices().size() !=0) nAnalyzed[kTagger]++;
      }
    }
  }
  
  for(ESource source : sources){
    for(double fraction=0; fraction<=1.0; fraction+=0.05){
      trueEventsByFraction.AddPoint(fraction, byFraction.nTrue[source][fraction], nAnalyzed[source], source);
      trueEventsByCompatibility.AddPoint(fraction, byCompatibility.nTrue[source][fraction], nAnalyzed[source], source);
      
      fakesEventsByFraction.AddPoint(fraction, byFraction.nFakes[source][fraction] , nAnalyzed[source], source);
      fakeEventsByCompatibility.AddPoint(fraction, byCompatibility.nFakes[source][fraction] , nAnalyzed[source], source);
      
      rocFraction.AddPointROC(byFraction.nFakes[source][fraction],
                              byFraction.nTrue[source][fraction], nAnalyzed[source], source);
      rocCompatibility.AddPointROC(byCompatibility.nFakes[source][fraction],
                                   byCompatibility.nTrue[source][fraction], nAnalyzed[source], source);
    }
    
    for(double nHits=0; nHits<=21.0; nHits+=1.0){
      trueEventsByHits.AddPoint(nHits, byHits.nTrue[source][nHits], nAnalyzed[source], source);
    }
  }
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config = ConfigManager(configPath);
  
  int part = 0;
  nEvents = 100;
  
  if(argc==3){
    nEvents = atoi(argv[1]);
    part = atoi(argv[2]);
  }
  cout<<"Running part "<<part<<" with "<<nEvents<<" events"<<endl;
  eventOffset = part*nEvents;
  
  // how many of helix hits are the true pion hits
  Hists trackTruthFraction("track_truth_fraction", 100, 0, 1);
  // track compatibility with pion helix
  Hists trackToPionSimilarity("track_compatibility", 100, 0, 1);
  
  // Events with true helix, depending on fraction
  Graphs trueEventsByFraction("true_events_by_fraction");
  // Events with true helix, depending on n hits
  Graphs trueEventsByHits("true_events_by_hits");
  // Events with true helix, depending on cmpatibility fraction
  Graphs trueEventsByCompatibility("true_events_by_compatibility");

  
  // Average fraction of fake helices, where definition of fake depends on fraction of hits on helix that are pion hits
  Graphs fakesEventsByFraction("fake_events_by_fraction");
  
  // Average fraction of fake helices, where definition of fake depends on fraction of pion hits that are on helix
  Graphs fakeEventsByCompatibility("fake_events_by_compatibility");
  
  Graphs rocFraction("roc_fraction");
  Graphs rocCompatibility("roc_compatibility");
  
  fillHists(trackTruthFraction, trackToPionSimilarity,
            trueEventsByFraction, trueEventsByHits, trueEventsByCompatibility,
            fakesEventsByFraction, fakeEventsByCompatibility,
            rocFraction, rocCompatibility);
  
  
  
  TCanvas *c1 = new TCanvas("c1", "c1", 1800, 1800);
  c1->Divide(3, 3);
  
  c1->cd(1);
  gPad->SetLogy(true);
  trackTruthFraction.DrawNormalized();
  
  c1->cd(2);
  gPad->SetLogy(true);
  trackToPionSimilarity.DrawNormalized();
  
  c1->cd(3); trueEventsByFraction.Draw("P");
  c1->cd(4); trueEventsByHits.Draw("P");
  c1->cd(5); trueEventsByCompatibility.Draw("P");
  c1->cd(6); fakesEventsByFraction.Draw("P");
  c1->cd(7); fakeEventsByCompatibility.Draw("P");
  c1->cd(8); rocFraction.Draw("P");
  c1->cd(9); rocCompatibility.Draw("P");
  
  c1->Update();
  
  TFile *outFile = new TFile((outPath+"_p"+to_string(part)+".root").c_str(), "recreate");
  outFile->cd();
  
  trackTruthFraction.Write();
  trackToPionSimilarity.Write();
  trueEventsByFraction.Write();
  trueEventsByHits.Write();
  trueEventsByCompatibility.Write();
  fakesEventsByFraction.Write();
  fakeEventsByCompatibility.Write();
  
  outFile->Close();
  
  theApp.Run();
  return 0;
}
