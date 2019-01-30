//
//  helixFitter.cpp
//
//  Created by Jeremi Niedziela on 17/12/2018.
//

#include "Helpers.hpp"
#include "HelixProcessor.hpp"
#include "Circle.hpp"
#include "PointsProcessor.hpp"
#include "Event.hpp"
#include "EventSet.hpp"
#include "Fitter.hpp"
#include "Display.hpp"
#include "ConfigManager.hpp"
#include "MonitorsManager.hpp"
#include "TrackProcessor.hpp"

int verbosityLevel = 2;

string configPath = "configs/helixFitter.md";
unique_ptr<PointsProcessor> pointsProcessor;
unique_ptr<HelixProcessor> helixProcessor;
unique_ptr<MonitorsManager> monitorsManager;

void InjectPionPointsToCollectionOfPoints(const unique_ptr<Helix> &pionHelix,
                                          shared_ptr<vector<Point>> pixelPoints)
{
  shared_ptr<vector<Point>> pionPoints = pionHelix->GetPoints();
  pixelPoints->insert(pixelPoints->end(),pionPoints->begin(), pionPoints->end());
}

void PerformTests(int &nSuccess, int &nFullSuccess)
{
  int nTests = config->nTests;
  auto fitter = make_unique<Fitter>();
  auto trackProcessor = make_unique<TrackProcessor>();
  
  for(int i=0;i<nTests;i++){
    if(verbosityLevel >= 2){
      cout<<"\n========================================================"<<endl;
      cout<<"Test iter:"<<i<<endl;
    }
    
    shared_ptr<Track> track = trackProcessor->GetRandomTrack(config->nTrackHits, config->maxEta);
    
    shared_ptr<vector<Point>> pixelPoints = pointsProcessor->GetRandomPoints(config->nNoiseHits);
//    unique_ptr<Event> event = make_unique<Event>();
//    event->SetRunNumber(297100);
//    event->SetLumiSection(136);
//    event->SetEventNumber(245000232);
//    shared_ptr<vector<Point>> pixelPoints = event->GetTrackerHits();
    
    unique_ptr<Helix> pionHelix = helixProcessor->GetRandomPionHelix(track);
    if(config->injectPionHits) InjectPionPointsToCollectionOfPoints(pionHelix, pixelPoints);
    
    unique_ptr<Helix> bestHelix = fitter->GetBestFittingHelix(pixelPoints, track);
    monitorsManager->FillMonitors(bestHelix, pionHelix, track);
    
    auto successCode = monitorsManager->GetFittingStatus(bestHelix, pionHelix);
    
    if(successCode == MonitorsManager::kFail)         continue;
    if(successCode == MonitorsManager::kSuccess)      nSuccess++;
    if(successCode == MonitorsManager::kFullSuccess){
      nSuccess++;
      nFullSuccess++;
    }
      
    if(verbosityLevel >= 2){
      cout<<"Pion helix:"; pionHelix->Print();
      cout<<"Fitted helix:"; bestHelix->Print();
      
      if(bestHelix->GetCharge() != pionHelix->GetCharge()){
        cout<<"\n\nwrong charge\n\n"<<endl;
        cout<<"best charge:"<<bestHelix->GetCharge()<<endl;
        cout<<"true charge:"<<pionHelix->GetCharge()<<endl;
        cout<<"best t shift:"<<bestHelix->GetTmin()<<endl;
        cout<<"true t shift:"<<pionHelix->GetTmin()<<endl;
      }
    }
  }
}

void ScanParameter()
{
  string paramName = "n_noise_hits";
  double paramMin = 0;
  double paramMax = 1500;
  double paramStep = 100;
  
  TH1D *eff_vs_param = new TH1D(("eff_vs_"+paramName).c_str(),
                                ("eff_vs_"+paramName).c_str(),
                                (paramMax-paramMin)/paramStep+1,paramMin,paramMax);
  
  TH1D *full_eff_vs_param = new TH1D(("full_eff_vs_"+paramName).c_str(),
                                     ("full_eff_vs_"+paramName).c_str(),
                                     (paramMax-paramMin)/paramStep+1,paramMin,paramMax);
  
  for(double param=paramMin;param<paramMax;param+=paramStep){
    cout<<"param:"<<param<<endl;
    config->nNoiseHits = param;
    
    int nTests = config->nTests;
    int nFullSuccess = 0;
    int nSuccess = 0;
    PerformTests(nSuccess, nFullSuccess);
    
    eff_vs_param->SetBinContent(eff_vs_param->GetXaxis()->FindFixBin(param), nSuccess/(double)nTests);
    full_eff_vs_param->SetBinContent(full_eff_vs_param->GetXaxis()->FindFixBin(param), nFullSuccess/(double)nTests);
  }
  
  TCanvas *c2 = new TCanvas("c2","c2",1280,800);
  c2->Divide(2,2);
  c2->cd(1);
  eff_vs_param->Draw();
  c2->cd(2);
  full_eff_vs_param->Draw();
  c2->Update();
  eff_vs_param->SaveAs(("eff_vs_"+paramName+".root").c_str());
  full_eff_vs_param->SaveAs(("full_eff_vs_"+paramName+".root").c_str());
  
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  config          = make_unique<ConfigManager>(configPath);
  pointsProcessor = make_unique<PointsProcessor>();
  helixProcessor  = make_unique<HelixProcessor>();
  monitorsManager = make_unique<MonitorsManager>();
  
  auto startTime = now();
  
  int nTests = config->nTests;
  int nSuccess, nFullSuccess;
  nFullSuccess = nSuccess = 0;
  PerformTests(nSuccess, nFullSuccess);
  cout<<"Percentage of successful fits:"<<nSuccess/(double)nTests<<endl;
  cout<<"Percentage of fully successful fits:"<<nFullSuccess/(double)nTests<<endl;
  
  monitorsManager->PlotAndSaveMonitors();
  
//  ScanParameter();
  config->Print();
  
  double timeElapsed = duration(startTime, now());
  cout<<"Average time per event:"<<timeElapsed/nTests<<" seconds"<<endl;
  
  theApp.Run();
  return 0;
}
