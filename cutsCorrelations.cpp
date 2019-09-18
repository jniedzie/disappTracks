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

/**
 Creates ABCD histogram. In case of backgrounds, all samples are merged into the same histogram.
 \param events Events with which histogram will be filled
 \param criticalMet Position of bins border on MET axis
 \param criticalDedx Position of bins border on dE/dx axis
 \param dataType Specifies whether background or signal events should be analyzed
 \param setIter For signal, specified which samples to use
 */
TH2D* GetABCDplot(const EventSet &events, double criticalMet, double criticalDedx,
                  xtracks::EDataType dataType, int setIter=0)
{
  CutsManager cutsManager;
  EventCut eventCut; TrackCut trackCut; JetCut jetCut; LeptonCut leptonCut;
  cutsManager.GetCuts(eventCut, trackCut, jetCut, leptonCut);
  
  EventSet eventsA(events), eventsB(events), eventsC(events), eventsD(events);
  string title = "ABCD_"+to_string_with_precision(criticalMet, 0)+"_"+to_string_with_precision(criticalDedx, 1);
  if(dataType == xtracks::kSignal)  title += ("_"+signalName[setIter]);
  else                              title += "_Backgrounds";
  
  float binsMet[]   = { 0, (float)criticalMet   , 1000 };
  float binsDedx[]  = { 0, (float)criticalDedx  , 6.0 };
  
  TH2D *hist = new TH2D(title.c_str(), title.c_str(), 2, binsDedx, 2, binsMet);
  
  // A
  trackCut.SetDedxPerCluster(range<double>(0, criticalDedx));
  eventCut.SetMetNoMuPt(range<double>(criticalMet,inf));
  eventsA.ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  // B
  trackCut.SetDedxPerCluster(range<double>(0, criticalDedx));
  eventCut.SetMetNoMuPt(range<double>(0,criticalMet));
  eventsB.ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  // C
  trackCut.SetDedxPerCluster(range<double>(criticalDedx, inf));
  eventCut.SetMetNoMuPt(range<double>(criticalMet,inf));
  eventsC.ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  // D
  trackCut.SetDedxPerCluster(range<double>(criticalDedx, inf));
  eventCut.SetMetNoMuPt(range<double>(0,criticalMet));
  eventsD.ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  if(dataType == xtracks::kBackground){
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      if(!config.runBackground[iBck]) continue;
      
      hist->Fill(0.0                  , 1.1*criticalMet , eventsA.weightedSize(dataType, iBck));
      hist->Fill(0.0                  , 0.0             , eventsB.weightedSize(dataType, iBck));
      hist->Fill(1.1*criticalDedx     , 1.1*criticalMet , eventsC.weightedSize(dataType, iBck));
      hist->Fill(1.1*criticalDedx     , 0.0             , eventsD.weightedSize(dataType, iBck));
    }
  }
  else if(dataType == xtracks::kSignal){
    hist->Fill(0.0                  , 1.1*criticalMet , eventsA.weightedSize(dataType, setIter));
    hist->Fill(0.0                  , 0.0             , eventsB.weightedSize(dataType, setIter));
    hist->Fill(1.1*criticalDedx     , 1.1*criticalMet , eventsC.weightedSize(dataType, setIter));
    hist->Fill(1.1*criticalDedx     , 0.0             , eventsD.weightedSize(dataType, setIter));
  }
  
  return hist;
}

/**
 Returns number of counts in A, B, C and D regions determined by criticalMet and criticalDedx.
 A: MET > criticalMet, dE/dx < criticalDedx
 B: MET < criticalMet, dE/dx < criticalDedx
 C: MET > criticalMet, dE/dx > criticalDedx
 D: MET < criticalMet, dE/dx > criticalDedx
 \param events Events from which ABCD will be calculated
 \param criticalMet Position of bins border on MET axis
 \param criticalDedx Position of bins border on dE/dx axis
 \param dataType Specifies whether background or signal events should be analyzed
 \param setIter For signal, specify which sample to use
 */
tuple<double,double,double,double> GetABCD(const EventSet &events,
                                           double criticalMet,
                                           double criticalDedx,
                                           xtracks::EDataType dataType, int setIter=0)
{
  TH2D *abcdPlot = GetABCDplot(events, criticalMet, criticalDedx, dataType, setIter);
  
  double A = abcdPlot->GetBinContent(1, 2);
  double B = abcdPlot->GetBinContent(1, 1);
  double C = abcdPlot->GetBinContent(2, 2);
  double D = abcdPlot->GetBinContent(2, 1);
  
  return make_tuple(A,B,C,D);
}

/**
 Returns expected and predicted yield in the signal bin together with errors.
 \param events Events from which ABCD will be calculated
 \param criticalMet Position of bins border on MET axis
 \param criticalDedx Position of bins border on dE/dx axis
 \param dataType Specifies whether background or signal events should be analyzed
 \param setIter For signal, specify which sample to use
 \return tuple containing: expected yield, expected yield error, predicted yield, prediction error
 */
tuple<double,double,double,double> GetDifference(const EventSet &events,
                                                 double criticalMet,
                                                 double criticalDedx,
                                                 xtracks::EDataType dataType, int setIter=0)
{
  auto [A,B,C,D] = GetABCD(events, criticalMet, criticalDedx, dataType, setIter);
  
  double expected = C;
  double predicted = A/B*D;
  
  double expectedErr = sqrt(expected);
  double predictedErr = sqrt(1/A+1/B+1/D)*predicted;
  
  return make_tuple(expected, expectedErr, predicted, predictedErr);
}


/// Creates a map of 2D histograms showing correlations between different variables
map<string, TH2D*> CreateCorrelationHistograms()
{
  //  name          nBinsX minX  maxX    nBinsY minY  maxY
  map<string, tuple<int, double, double, int, double, double>> histParams = {
    {"iso_vs_dedx"                , { 100  , 0 , 10   , 100 , 0   , 0.1   }},
    {"met_vs_dedx"                , { 100  , 0 , 10   , 100 , 200 , 1200  }},
    {"trackPt_vs_missing"         , { 100  , 0 , 1000 , 20  , 0   , 20    }},
    {"deltaJetTrack_vs_missing"   , { 20   , 0 , 2    , 20  , 0   , 20    }},
  };
  
  map<string, tuple<string, string>> histAxesNames = {
    {"iso_vs_dedx"                , { "Average dE/dx (MeV/cm)"  , "Relative isolation"      }},
    {"met_vs_dedx"                , { "Average dE/dx (MeV/cm)"  , "MET p_{T} (GeV)"         }},
    {"trackPt_vs_missing"         , { "Missing outer tracker hits" , "Track p_{T}"          }},
    {"deltaJetTrack_vs_missing"   , { "Missing outer tracker hits" , "#Delta R(jet,track)"  }},
  };
  
  map<string, TH2D*> correlationHists;
  
  for(auto &[name, params] : histParams){
    auto &[nBinsX, minX, maxX, nBinsY, minY, maxY] = params;
    correlationHists[name] = new TH2D(name.c_str(), name.c_str(), nBinsX, minX, maxX, nBinsY, minY, maxY);
    auto &[titleX, titleY] = histAxesNames[name];
    correlationHists[name]->GetXaxis()->SetTitle(titleX.c_str());
    correlationHists[name]->GetYaxis()->SetTitle(titleY.c_str());
  }
  return correlationHists;
}

/// Fills correlation histograms from background events
void FillCorrelationHistograms(map<string, TH2D*> &correlationHists, const EventSet &events)
{
  for(int iBck=0; iBck<kNbackgrounds; iBck++){
    for(int iEvent=0;iEvent<events.size(xtracks::kBackground, iBck);iEvent++){
      auto event = events.At(xtracks::kBackground, iBck, iEvent);
      
      for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
        auto track = event->GetTrack(iTrack);
        
        double avgDedx = track->GetAverageDedx();
        
        correlationHists["iso_vs_dedx"]->Fill(avgDedx, track->GetRelativeIsolation());
        correlationHists["met_vs_dedx"]->Fill(avgDedx, event->GetMetPt());
        correlationHists["trackPt_vs_missing"]->Fill(track->GetPt(), track->GetNmissingOuterTrackerHits());
        
        for(int iJet=0;iJet<event->GetNjets();iJet++){
          auto jet = event->GetJet(iJet);
          
          double deltaR = sqrt(pow(track->GetPhi() - jet->GetPhi(),2)+pow(track->GetEta() - jet->GetEta(),2));
          correlationHists["deltaJetTrack_vs_missing"]->Fill(deltaR, track->GetNmissingOuterTrackerHits());
        }
      }
    }
  }
}

/// Draws correlation histograms in a canvas
void DrawCorrelationHistograms(const map<string, TH2D*> &correlationHists)
{
  TCanvas *c1 = new TCanvas("c1","c1",1000,1500);
  if(correlationHists.size() <= 4) c1->Divide(2,2);
  else if(correlationHists.size() <= 6) c1->Divide(2,3);
  else if(correlationHists.size() <= 9) c1->Divide(3,3);
  int iPad=1;
  for(auto &[_, hist] : correlationHists){
    c1->cd(iPad++);
    hist->Draw("colz");
  }
  c1->Update();
}

/**
 Searches for criticalMet and criticalDedx for which relative difference between predicted and expected number
 of events in the signal bin is the smallers. Uses background events to find those values.
 \param events Events to be analyzed
 \param minMet Starting value of MET pt
 \param maxMet Maximum value of MET pt
 \param stepMet Step of MET pt
 \param minDedx Starting value of dE/dx
 \param maxDedx Maximum value of dE/dx
 \param stepDedx Step of dE/dx
 \return returns tuple containing: best critical MET, best critical dE/dx
 */
tuple<double,double> GetBestMetAndDedx(const EventSet &events,
                                       double minMet, double maxMet, double stepMet,
                                       double minDedx, double maxDedx, double stepDedx)
{
  double bestRelativeDifference = inf;
  double bestMet=0, bestDedx=0, bestExpected=0, bestExpectedErr=0, bestPredicted=0, bestPredictedErr=0;
  
  for(double criticalMet = minMet; criticalMet <= maxMet; criticalMet += stepMet){
    cout<<"MET pt:"<<criticalMet<<endl;
    
    for(double criticalDedx = minDedx; criticalDedx <= maxDedx; criticalDedx += stepDedx){
      auto [expected, expectedErr, predicted, predictedErr] = GetDifference(events,
                                                                            criticalMet, criticalDedx,
                                                                            xtracks::kBackground);
      
      double relDiff = fabs(expected-predicted)/expected;
      
      if(relDiff < bestRelativeDifference){
        bestRelativeDifference = relDiff;
        bestMet = criticalMet;
        bestDedx = criticalDedx;
        bestExpected = expected;
        bestPredicted = predicted;
        bestExpectedErr = expectedErr;
        bestPredictedErr = predictedErr;
      }
    }
  }
  return make_tuple(bestMet, bestDedx);
}

/// Draws and saves ABCD plots for given values of critical MET and critical dE/dx
void DrawAndSaveABCDplots(const EventSet &events, double criticalMet, double criticalDedx)
{
  TCanvas *abcdCanvas = new TCanvas("ABCD", "ABCD", 1000, 1500);
  abcdCanvas->Divide(4,3);
  
  TFile *outFile = new TFile("results/abcd_plots.root","recreate");
  
  gStyle->SetOptStat(0);
  
  cout<<"Plotting background ABCD"<<endl;
  TH2D *abcdPlotBackgrounds = GetABCDplot(events, criticalMet, criticalDedx, xtracks::kBackground);
  int iPad=1;
  abcdCanvas->cd(iPad++);
  abcdPlotBackgrounds->SetMarkerSize(3.0);
  abcdPlotBackgrounds->Draw("colzText");
  outFile->cd();
  abcdPlotBackgrounds->Write();
  
  for(int iSig=0; iSig<kNsignals; iSig++){
    if(!config.runSignal[iSig]) continue;
    cout<<"Plotting "<<signalTitle[iSig]<<" ABCD"<<endl;
    TH2D *abcdPlot = GetABCDplot(events, criticalMet, criticalDedx, xtracks::kSignal, iSig);
    abcdCanvas->cd(iPad++);
    abcdPlot->SetMarkerSize(3.0);
    abcdPlot->Draw("colzText");
    outFile->cd();
    abcdPlot->Write();
  }
  
  abcdCanvas->Update();
  abcdCanvas->Write();
  outFile->Close();
}

/// Starting point of the application
int main(int argc, char* argv[])
{
  cout.imbue(locale("de_DE"));
  TApplication *theApp = new TApplication("App", &argc, argv);
  config = ConfigManager("configs/analysis.md");
  
  // All events with initial cuts only
  EventSet events;
  string prefix = "after_L"+to_string_with_precision(config.params["cuts_level"], 0)+"/"+config.category+"/";
  events.LoadEventsFromFiles(prefix);
  
  // Draw correlation plots
//  map<string, TH2D*> correlationHists = CreateCorrelationHistograms();
//  FillCorrelationHistograms(correlationHists, events);
//  DrawCorrelationHistograms(correlationHists);
  
  // Find the best values of critical MET and critical dE/dx
//  double minMet = 210, maxMet = 300, stepMet = 10;
//  double minDedx = 3.0, maxDedx = 4.5, stepDedx = 0.1;
//  auto [bestMet, bestDedx] = GetBestMetAndDedx(events, minMet, maxMet, stepMet, minDedx, maxDedx, stepDedx, xtracks::kBackground);
  
  double bestMet=240, bestDedx=4.3;
  
  // Print restuls
  cout<<"Best MET: "<<bestMet<<"\tdE/dx: "<<bestDedx<<endl;
  
  auto [expected, expectedErr, predicted, predictedErr] = GetDifference(events, bestMet, bestDedx, xtracks::kBackground);
  cout<<"Expected: "<<expected<<" +/- "<<expectedErr<<"\tpredicted: "<<predicted<<" +/- "<<predictedErr<<endl;
  
  double relDiff = fabs(expected-predicted)/expected;
  double error = sqrt(pow(predicted*expectedErr, 2) + pow(expected*predictedErr, 2)) / pow(expected, 2);
  double ratio = predicted/expected;
  
  cout<<"Relative difference: "<<relDiff<<" +/- "<<error<<endl;
  cout<<"Ratio: "<<ratio<<endl;
  
  
  DrawAndSaveABCDplots(events, bestMet, bestDedx);
  
  theApp->Run();
  return 0;
}




