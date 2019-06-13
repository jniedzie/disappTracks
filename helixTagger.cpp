//
//  helixTagger.cpp
//
//  Created by Jeremi Niedziela on 25/01/2019.
//

#include "Helpers.hpp"
#include "Fitter.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"
#include "EventSet.hpp"

string configPath = "configs/helixTagger.md";
string cutLevel = "after_L2/4layers/";//after_L1/";

int nEvents = 45;
int nTests = 10;

double SetParamValue(int iTest){
  // here put a way to calculate param value based on the test iter:
  // then assign the value to the correct config parameter:

//  double paramValue = iTest;      config.doubleHitsMaxDistance = paramValue;
  
//  double paramValue = pow(10, -iTest+2);  config.seedMaxChi2 = paramValue;

//    double paramValue = -iTest*0.1;  config.seedMiddleHitDeltaPhi = range<double>(paramValue, config.seedMiddleHitDeltaPhi.GetMax());
//  double paramValue = iTest*0.1;  config.seedMiddleHitDeltaPhi = range<double>(0, paramValue);
//  double paramValue = iTest*10;   config.seedMiddleHitMaxDeltaZ = paramValue;
  
  double paramValue = -iTest*0.1;  config.seedLastHitDeltaPhi = range<double>(paramValue, config.seedLastHitDeltaPhi.GetMax());
//  double paramValue = iTest*0.1;  config.seedLastHitDeltaPhi = range<double>(0, paramValue);
//  double paramValue = iTest*10;   config.seedLastHitMaxDeltaZ = paramValue;

//  double paramValue = pow(10, -iTest+2);  config.trackMaxChi2 = paramValue;
//  double paramValue = iTest*0.1;          config.nextPointDeltaPhi = range<double>(0, paramValue);
//  double paramValue = iTest*10;           config.nextPointMaxDeltaZ = paramValue;
//  config.nextPointMaxDeltaXY = paramValue;
//  double paramValue = iTest+3;            config.trackMinNpoints = paramValue;
  
//  double paramValue = iTest;    config.mergingMaxDifferentPoints = paramValue;
//  double paramValue = iTest+3;  config.candidateMinNpoints = paramValue;

  return paramValue;
}

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;

shared_ptr<Event> GetEvent(int iEvent);
vector<shared_ptr<Point>> GetClustersNoEndcaps(const shared_ptr<Event> &event, bool removePionClusters);

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  config = ConfigManager(configPath);
  auto fitter = make_unique<Fitter>();
  TCanvas *canvas = new TCanvas("ROC", "ROC", 800,600);
  canvas->cd();
  
  
  TGraph *rocGraph[nTests];
  TF1 *rocFun[nTests];
  
  vector<tuple<double, double, double, double>> aucAndMaxEffForTest;
  
  vector<vector<double>> nHelixPointsSignal;// [iTest][iEvent]
  vector<vector<double>> nHelixPointsBackground;
  
  vector<vector<int>> nMaxHelixPointsSignal;// [iTest][iEvent]
  vector<vector<int>> nMaxHelixPointsBackground;
  
  for(int iTest=0; iTest<nTests; iTest++){
    vector<double> vec;
    vector<int> vecInt;
    
    for(int iEvent=0; iEvent<nEvents; iEvent++){
      vec.push_back(0.0);
      vecInt.push_back(0);
    }
    nHelixPointsSignal.push_back(vec);
    nHelixPointsBackground.push_back(vec);
    nMaxHelixPointsSignal.push_back(vecInt);
    nMaxHelixPointsBackground.push_back(vecInt);
  }
  
  int nAnalyzedEvents=0;
  
  for(auto iEvent=0; iEvent<nEvents; iEvent++){
    auto event = GetEvent(iEvent);
    
    auto pointsNoEndcapsSignal     = GetClustersNoEndcaps(event, false);
    auto pointsNoEndcapsBackground = GetClustersNoEndcaps(event, true);
    
    if(pointsNoEndcapsSignal.size()==0 || pointsNoEndcapsBackground.size()==0){
      cout<<"helixTagger -- no tracker hits for event "<<iEvent<<endl;
      //      continue;
    }
    
    cout<<"\n\n=================================================================\n"<<endl;
    cout<<"helixTagger -- processing event "<<iEvent<<endl;
    
    for(int iTest=0; iTest<nTests; iTest++){
      SetParamValue(iTest);
      
      for(auto &track : event->GetTracks()){
        
        vector<Helix> fittedHelicesSignal     = fitter->FitHelices(pointsNoEndcapsSignal, *track, *event->GetVertex());
        vector<Helix> fittedHelicesBackground = fitter->FitHelices(pointsNoEndcapsBackground, *track, *event->GetVertex());
        
        double averageHelixNpointsSignal=0;
        double averageHelixNpointsBackground=0;
        
        int maxHelixNpointsSignal = 0;
        int maxHelixNpointsBackground = 0;
        
        for(auto helix : fittedHelicesSignal){
          averageHelixNpointsSignal += helix.GetNpoints();
          if(helix.GetNpoints() > maxHelixNpointsSignal) maxHelixNpointsSignal = helix.GetNpoints();
          //      event->AddHelix(move(fittedHelix));
        }
        
        for(auto helix : fittedHelicesBackground){
          averageHelixNpointsBackground += helix.GetNpoints();
          if(helix.GetNpoints() > maxHelixNpointsBackground) maxHelixNpointsBackground = helix.GetNpoints();
        }
        
        averageHelixNpointsSignal /= fittedHelicesSignal.size()>0 ? fittedHelicesSignal.size() : 1;
        nHelixPointsSignal[iTest][iEvent] = averageHelixNpointsSignal;
        nMaxHelixPointsSignal[iTest][iEvent] = maxHelixNpointsSignal;
        
        averageHelixNpointsBackground /= fittedHelicesBackground.size()>0 ? fittedHelicesBackground.size() : 1;
        nHelixPointsBackground[iTest][iEvent] = averageHelixNpointsBackground;
        nMaxHelixPointsBackground[iTest][iEvent] = maxHelixNpointsBackground;
      }
    }
    nAnalyzedEvents++;
  }
  
  for(int iTest=0; iTest<nTests; iTest++){
    
    vector<double> efficiency;
    vector<double> fakeRate;
    
    vector<double> efficiencyMax;
    vector<double> fakeRateMax;
    
    int maxNpoints = 20;
    for(int iThreshold=0; iThreshold<maxNpoints; iThreshold++){
      efficiency.push_back(0);
      fakeRate.push_back(0);
      efficiencyMax.push_back(0);
      fakeRateMax.push_back(0);
    }
    
    
    for(double nPoints : nHelixPointsSignal[iTest]){
      for(int iThreshold=0; iThreshold<maxNpoints; iThreshold++){
        if(nPoints >= iThreshold) efficiency[iThreshold]+=1;
      }
    }
    for(double nPoints : nHelixPointsBackground[iTest]){
      for(int iThreshold=0; iThreshold<maxNpoints; iThreshold++){
        if(nPoints >= iThreshold) fakeRate[iThreshold]+=1;
      }
    }
    
    for(int nPoints : nMaxHelixPointsSignal[iTest]){
      for(int iThreshold=0; iThreshold<maxNpoints; iThreshold++){
        if(nPoints >= iThreshold) efficiencyMax[iThreshold]+=1;
      }
    }
    for(int nPoints : nMaxHelixPointsBackground[iTest]){
      for(int iThreshold=0; iThreshold<maxNpoints; iThreshold++){
        if(nPoints >= iThreshold) fakeRateMax[iThreshold]+=1;
      }
    }
    
    rocFun[iTest] = new TF1(Form("rocFun[%i]",iTest), "[2]*exp(-pow(x-[0],2)/(2*[1]*[1]))/sqrt(2*3.1415*[1]*[1])+[3]",0,1);
    rocFun[iTest]->SetParameter(0, 0); // μ
    rocFun[iTest]->SetParameter(1, 1); // σ
    rocFun[iTest]->SetParameter(2, 1); // a
    rocFun[iTest]->SetParameter(3, 1); // b
    
    rocGraph[iTest] = new TGraph();
    
    double paramValue = SetParamValue(iTest);
    
    double maxEff = -inf;
    double maxEffMax = -inf;
    double maxEffOverFake = -inf;
    
    cout<<"\nFake-eff for test param value: "<<paramValue<<"(test "<<iTest<<")"<<endl;
    
    for(int iThreshold=0; iThreshold<maxNpoints; iThreshold++){
      efficiency[iThreshold] /= nAnalyzedEvents;
      fakeRate[iThreshold] /= nAnalyzedEvents;
      efficiencyMax[iThreshold] /= nAnalyzedEvents;
      fakeRateMax[iThreshold] /= nAnalyzedEvents;
      
      if(efficiency[iThreshold] > maxEff && efficiency[iThreshold] != 1.0)
        maxEff = efficiency[iThreshold];
      
      if(efficiencyMax[iThreshold] > maxEffMax && efficiencyMax[iThreshold] != 1.0)
        maxEffMax = efficiencyMax[iThreshold];
      
      double sigmaApprox = 100*efficiencyMax[iThreshold]/sqrt(100*efficiencyMax[iThreshold]+100000*fakeRateMax[iThreshold]);
      
      if(sigmaApprox > maxEffOverFake) maxEffOverFake = sigmaApprox;
      
      cout<<fakeRateMax[iThreshold]<<"\t"<<efficiencyMax[iThreshold]<<endl;
      rocGraph[iTest]->SetPoint(iThreshold, fakeRateMax[iThreshold], efficiencyMax[iThreshold]);

//      cout<<fakeRate[iThreshold]<<"\t"<<efficiency[iThreshold]<<endl;
//      rocGraph[iTest]->SetPoint(iThreshold, fakeRate[iThreshold], efficiency[iThreshold]);
    }
    
    rocGraph[iTest]->SetMarkerStyle(20);
    rocGraph[iTest]->SetMarkerSize(2.0);
    rocGraph[iTest]->SetMarkerColor(iTest);
    
    rocGraph[iTest]->Draw(iTest==0 ? "AP" : "Psame");
    rocGraph[iTest]->Fit(rocFun[iTest]);
    
    double auc = rocFun[iTest]->Integral(0,1);
    
//    aucAndMaxEffForTest.push_back(make_tuple(paramValue, auc, maxEff));
    aucAndMaxEffForTest.push_back(make_tuple(paramValue, auc, maxEffMax, maxEffOverFake));
  }
  
  
  for(auto &[param, auc, eff, effOverFake] : aucAndMaxEffForTest){
    cout<<"Param: "<<param<<"\tAUC: "<<auc<<"\tmax eff: "<<eff<<"\teff/fake: "<<effOverFake<<endl;
  }
  
  //  cout<<"helixTagger -- saving events"<<endl;
  //  events.SaveEventsToFiles("afterHelixTagging/");
  //  cout<<"helixTagger -- finished"<<endl;
  
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
