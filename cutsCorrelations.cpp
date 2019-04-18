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

vector<tuple<double,double>> GetDifferencesForCriticalValues(shared_ptr<EventSet> &events,
                                               double criticalMet = 500, double critialDedx = 3.5)
{
  auto eventsA = shared_ptr<EventSet>(new EventSet(*events));
  auto eventsB = shared_ptr<EventSet>(new EventSet(*events));
  auto eventsC = shared_ptr<EventSet>(new EventSet(*events));
  auto eventsD = shared_ptr<EventSet>(new EventSet(*events));

  auto eventCut = EventCut();
  auto trackCut = unique_ptr<TrackCut>(new TrackCut());
  auto jetCut   = unique_ptr<JetCut>(new JetCut());
  auto leptonCut= unique_ptr<LeptonCut>(new LeptonCut());
  
  // + standard cuts to be applied after L2 selections
  eventCut.SetNtracks(range<int>(1,1));
  eventCut.SetNjets(range<int>(1,inf));
  eventCut.SetLeadingJetPt(range<double>(100,inf));
  eventCut.SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCut.SetLeadingJetNeHEF(range<double>(-inf,0.8));
  eventCut.SetLeadingJetChHEF(range<double>(0.1,inf));

  // A: B + S
  eventCut.SetMetPt(range<double>(criticalMet,inf));           // select Background
  trackCut->SetDedxPerCluster(range<double>(critialDedx,inf));  // select Signal
  eventsA->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);

  // B: B + B
  eventCut.SetMetPt(range<double>(criticalMet,inf));           // select Background
  trackCut->SetDedxPerCluster(range<double>(0,critialDedx));    // select Background
  eventsB->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);

  // C: S + S
  eventCut.SetMetPt(range<double>(0.0,criticalMet));           // select Signal
  trackCut->SetDedxPerCluster(range<double>(critialDedx,inf));  // select Signal
  eventsC->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);

  // D: S + B
  eventCut.SetMetPt(range<double>(criticalMet,inf));           // select Signal
  trackCut->SetDedxPerCluster(range<double>(0,critialDedx));    // select Background
  eventsD->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);

  vector<tuple<double,double>> results;
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!config.runBackground[iBck]){
      results.push_back(make_tuple(inf,inf));
      continue;
    }
    double expectedBck = eventsC->weightedSize(xtracks::kBackground, iBck);
    double predictedBck = eventsA->weightedSize(xtracks::kBackground, iBck)/eventsB->weightedSize(xtracks::kBackground, iBck)*eventsD->weightedSize(xtracks::kBackground, iBck);
    
    results.push_back(make_tuple(expectedBck,predictedBck));
    }
  return results;
}

double GetFraction(shared_ptr<EventSet> &events, double criticalMet, double critialDedx, double stepX)
{
  auto eventsA = shared_ptr<EventSet>(new EventSet(*events));
  auto eventsB = shared_ptr<EventSet>(new EventSet(*events));
  auto eventsC = shared_ptr<EventSet>(new EventSet(*events));
  auto eventsD = shared_ptr<EventSet>(new EventSet(*events));
  
  EventCut eventCut;
  auto trackCut = unique_ptr<TrackCut>(new TrackCut());
  auto jetCut   = unique_ptr<JetCut>(new JetCut());
  auto leptonCut= unique_ptr<LeptonCut>(new LeptonCut());
  
  // + standard cuts to be applied after L2 selections
  eventCut.SetNtracks(range<int>(1,1));
  eventCut.SetNjets(range<int>(1,inf));
  eventCut.SetLeadingJetPt(range<double>(100,inf));
  eventCut.SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCut.SetLeadingJetNeHEF(range<double>(-inf,0.8));
  eventCut.SetLeadingJetChHEF(range<double>(0.1,inf));
  
  // A: B + S
  eventCut.SetMetPt(range<double>(criticalMet,inf));           // select Background
  trackCut->SetDedxPerCluster(range<double>(critialDedx,critialDedx+stepX));  // select Signal
  eventsA->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  // B: B + B
  eventCut.SetMetPt(range<double>(0,criticalMet));           // select Background
  trackCut->SetDedxPerCluster(range<double>(critialDedx, critialDedx+stepX));    // select Background
  eventsB->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  int nBackgroundA = 0;
  int nBackgroundB = 0;
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!config.runBackground[iBck]) continue;
    
    nBackgroundA += eventsA->weightedSize(xtracks::kBackground, iBck);
    nBackgroundB += eventsB->weightedSize(xtracks::kBackground, iBck);
  }
  
  return nBackgroundA/(double)nBackgroundB;
}

int main(int argc, char* argv[])
{
  TApplication *theApp = new TApplication("App", &argc, argv);
  
  config = ConfigManager("configs/analysis.md");
  
  // All events with initial cuts only
  shared_ptr<EventSet> events = shared_ptr<EventSet>(new EventSet());
  events->LoadEventsFromFiles("after_L0/");
  
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->Divide(2,2);
  
  TH2D *dedx_vs_isolation = new TH2D("dedx_vs_isolation","dedx_vs_isolation",100,0,10,100,0,0.1);
  TH2D *met_vs_dedx = new TH2D("met_vs_dedx","met_vs_dedx",100,0,10,100,200,1200);
  TH2D *trackPt_vs_missing = new TH2D("trackPt_vs_missing","trackPt_vs_missing",100,0,1000,20,0,20);
  TH2D *deltaJetTrack_vs_missing = new TH2D("deltaJetTrack_vs_missing","deltaJetTrack_vs_missing",20,0,2,20,0,20);
  
  for(int iEvent=0;iEvent<events->size(xtracks::kBackground, kQCD);iEvent++){
    auto event = events->At(xtracks::kBackground, kQCD, iEvent);
    
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      auto track = event->GetTrack(iTrack);
      dedx_vs_isolation->Fill(track->GetDeDxForHit(0),track->GetRelativeIsolation());
      
      double avgDedx=0;
      for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){avgDedx += track->GetAverageDedx();}
      avgDedx /= event->GetNtracks();
      
      met_vs_dedx->Fill(avgDedx, event->GetMetPt());
      
      trackPt_vs_missing->Fill(track->GetPt(),track->GetNmissingOuterTrackerHits());
      
      
      for(int iJet=0;iJet<event->GetNjets();iJet++){
        auto jet = event->GetJet(iJet);
        
        double deltaR = sqrt(pow(track->GetPhi() - jet->GetPhi(),2)+pow(track->GetEta() - jet->GetEta(),2));
        deltaJetTrack_vs_missing->Fill(deltaR,track->GetNmissingOuterTrackerHits());
      }
    }
  }
  c1->cd(1);
  dedx_vs_isolation->GetXaxis()->SetTitle("dedx in layer[0]");
  dedx_vs_isolation->GetXaxis()->SetTitle("Relative isolation");
  dedx_vs_isolation->Draw("colz");
  c1->cd(2);
  trackPt_vs_missing->GetXaxis()->SetTitle("Track p_{T}");
  trackPt_vs_missing->GetXaxis()->SetTitle("Missing outer tracker hits");
  trackPt_vs_missing->Draw("colz");
  c1->cd(3);
  deltaJetTrack_vs_missing->GetXaxis()->SetTitle("#Delta R(jet,track)");
  deltaJetTrack_vs_missing->GetXaxis()->SetTitle("Missing outer tracker hits");
  deltaJetTrack_vs_missing->Draw("colz");
//  c1->cd(4);
//  met_vs_dedx->GetXaxis()->SetTitle("Average dE/dx (MeV)");
//  met_vs_dedx->GetXaxis()->SetTitle("MET p_{T} (GeV)");
//  met_vs_dedx->Draw("colz");
  
  auto eventCut = EventCut();
  auto trackCut = unique_ptr<TrackCut>(new TrackCut());
  auto jetCut   = unique_ptr<JetCut>(new JetCut());
  auto leptonCut= unique_ptr<LeptonCut>(new LeptonCut());
  
  // pick category
  trackCut->SetNpixelLayers(range<int>(3,3));
  eventCut.SetNtracks(range<int>(1,1));
  
  // other L2 cuts
//  trackCut->SetCaloEmEnergy(range<double>(0.0,0.4));
//  trackCut->SetDedxPerCluster(range<double>(2.0,inf));
//  trackCut->SetTrackMetDeltaPhi(range<double>(-2.3,2.3));
  
//  trackCut->SetCaloEmEnergy(range<double>(0.0,3.0));
//  trackCut->SetDedxPerCluster(range<double>(2.0,inf));
//  eventCut.SetJetMetDeltaPhi(range<double>(0.7,inf));
//  trackCut->SetTrackMetDeltaPhi(range<double>(-2.3,2.3));
  
  
  // + standard cuts to be applied after L2 selections
  eventCut.SetNjets(range<int>(1,inf));
  eventCut.SetLeadingJetPt(range<double>(100,inf));
  eventCut.SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCut.SetLeadingJetNeHEF(range<double>(-inf,0.8));
  eventCut.SetLeadingJetChHEF(range<double>(0.1,inf));
  
  events->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  events->DrawStandardPlots();
  
  /*
  vector<tuple<double,double>> results = GetDifferencesForCriticalValues(events, 0.15, 4.0);
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    auto [expected, predicted] = results[iBck];
    double diff = fabs(expected - predicted)/expected;
    
    cout<<"\n"<<backgroundTitle[iBck]<<endl;
    cout<<"Expected:"<<expected<<"\tABCD:"<<predicted<<endl;
    cout<<"Relative difference:"<< diff<<endl;
  }
  */
  
  
  vector<tuple<tuple<double,double>,double,double>> bestResults;
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    bestResults.push_back(make_tuple(make_tuple(0,inf),inf,inf));
  }
  
  TH1D *ratio = new TH1D("ratio","ratio",50,1.0,6.0);
  double criticalMet = 300;
  double dedxStep = 0.1;
//  for(double criticalMet = 200;criticalMet<310;criticalMet+=10){
//    cout<<"MET pt:"<<criticalMet<<endl;
    for(double criticalDedx = 1.0;criticalDedx<6.0;criticalDedx+=dedxStep){
      cout<<"critical dedx:"<<criticalDedx<<endl;
      double fraction =  GetFraction(events, criticalMet, criticalDedx, dedxStep);
      ratio->Fill(criticalDedx, fraction);
      
      /*
      vector<tuple<double,double>> results = GetDifferencesForCriticalValues(events, criticalMet, criticalDedx);
      
      for(int iBck=0;iBck<kNbackgrounds;iBck++){
        auto [expected, predicted] = results[iBck];
        double relDiff = fabs(expected-predicted)/expected;
        
        auto [bestExpected, bestPredicted] = get<0>(bestResults[iBck]);
        double bestDiff = fabs(bestExpected - bestPredicted)/bestExpected;
        
        if(relDiff < bestDiff){
          bestResults[iBck] = make_tuple(make_tuple(expected,predicted),criticalMet,criticalDedx);
        }
      }
       */
    }
//  }
  c1->cd(4);
  ratio->Draw();
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    auto [bestExpected, bestPredicted] = get<0>(bestResults[iBck]);
    double bestCriticalMet = get<1>(bestResults[iBck]);
    double bestCriticalDedx = get<2>(bestResults[iBck]);
    double bestDiff = fabs(bestExpected - bestPredicted)/bestExpected;
    
    cout<<"\n"<<backgroundTitle[iBck]<<endl;
    cout<<"Expected:"<<bestExpected<<"\tABCD:"<<bestPredicted<<endl;
    cout<<"Relative difference:"<< bestDiff<<endl;
    cout<<"Critical MET:"<<bestCriticalMet<<endl;
    cout<<"Critical dedx:"<<bestCriticalDedx<<endl;
  }
  
  c1->Update();
  theApp->Run();
  
  return 0;
}




