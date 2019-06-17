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
int nTests = 15;

double SetParamValue(int iTest){
//  return 0;
  // here put a way to calculate param value based on the test iter:
  // then assign the value to the correct config parameter:

//  double paramValue = iTest+5;
//  config.doubleHitsMaxDistance = paramValue;
  
//  double paramValue = pow(10, -iTest+2);  config.seedMaxChi2 = paramValue;

//  double paramValue = -2.0+iTest*0.1;
//  config.seedMiddleHitDeltaPhi = range<double>(paramValue, config.seedMiddleHitDeltaPhi.GetMax());
//  double paramValue = -0.9+iTest*0.1;
//  config.seedMiddleHitDeltaPhi = range<double>(config.seedMiddleHitDeltaPhi.GetMin(), paramValue);
//  double paramValue = iTest*20;   config.seedMiddleHitMaxDeltaZ = paramValue;
  
//  double paramValue = -2.0+iTest*0.1;
//  config.seedLastHitDeltaPhi = range<double>(paramValue, config.seedLastHitDeltaPhi.GetMax());
//  double paramValue = -1.0+iTest*0.1;
//  config.seedLastHitDeltaPhi = range<double>(config.seedLastHitDeltaPhi.GetMin(), paramValue);
  double paramValue = iTest*20;   config.seedLastHitMaxDeltaZ = paramValue;

//  double paramValue = (iTest+5)*pow(10, -4);
//  config.trackMaxChi2 = paramValue;
  
//  double paramValue = -1.5+iTest*0.1;
//  config.nextPointDeltaPhi = range<double>(paramValue, config.nextPointDeltaPhi.GetMax());
//  double paramValue = -0.7+iTest*0.1;
//  config.nextPointDeltaPhi = range<double>(config.nextPointDeltaPhi.GetMin(), paramValue);
//  double paramValue = iTest*20;
//  config.nextPointMaxDeltaZ = paramValue;
//  double paramValue = iTest*20;
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
TF1* GetRocFunction();

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  config = ConfigManager(configPath);
  auto fitter = make_unique<Fitter>();
  TCanvas *canvas = new TCanvas("ROC", "ROC", 2880,1800);
  canvas->Divide(4,2);
  
  
  TGraph *rocGraph[nTests];
  TGraph *rocGraphMax[nTests];
  TGraph *rocGraphLayers[nTests];
  
  TF1 *rocFun[nTests];
  TF1 *rocFunMax[nTests];
  TF1 *rocFunLayers[nTests];
  
  TH1D *nHelicesSignal     = new TH1D("n helices", "n helices", 100, 0, 100);
  TH1D *nHelicesBackground = new TH1D("n helices bck", "n helices bck", 100, 0, 100);
  
  TH1D *maxNhitsSignal     = new TH1D("max n hits", "max n hits", 20, 0, 20);
  TH1D *maxNhitsBackground = new TH1D("max n hits bck", "max n hits bck", 20, 0, 20);
  
  TH1D *maxNlayersSignal     = new TH1D("max n layers", "max n layers", 20, 0, 20);
  TH1D *maxNlayersBackground = new TH1D("max n layers bck", "max n layers bck", 20, 0, 20);
  
  TH1D *avgNhitsSignal     = new TH1D("avg n hits", "avg n hits", 20, 0, 20);
  TH1D *avgNhitsBackground = new TH1D("avg n hits bck", "avg n hits bck", 20, 0, 20);
  
  vector<tuple<double, double, double, double>> aucAndMaxEffForTest;
  vector<tuple<double, double, double, double>> aucAndMaxEffForTestMax;
  vector<tuple<double, double, double, double>> aucAndMaxEffForTestLayers;
  
  vector<vector<double>> nHelixPointsSignal;// [iTest][iEvent]
  vector<vector<double>> nHelixPointsBackground;
  
  vector<vector<int>> nMaxHelixPointsSignal;// [iTest][iEvent]
  vector<vector<int>> nMaxHelixPointsBackground;
  
  vector<vector<int>> nHelixLayersSignal;// [iTest][iEvent]
  vector<vector<int>> nHelixLayersBackground;
  
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
    nHelixLayersSignal.push_back(vecInt);
    nHelixLayersBackground.push_back(vecInt);
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
        
        int maxHelixNlayersSignal = 0;
        int maxHelixNlayersBackground = 0;
        
        for(auto helix : fittedHelicesSignal){
          averageHelixNpointsSignal += helix.GetNpoints();
          if(helix.GetNpoints() > maxHelixNpointsSignal) maxHelixNpointsSignal = helix.GetNpoints();
          
          unordered_set<int> layers;
          
          for(int iPoint=0; iPoint<helix.GetNpoints(); iPoint++){
            auto point = helix.GetPoints()[iPoint];
            int layer = point->GetLayer();
            if(layer > 0){
              if(iPoint < helix.GetFirstTurningPointIndex()) layers.insert( layer);
              else                                           layers.insert(-layer);
            }
          }
          if(layers.size() > maxHelixNlayersSignal) maxHelixNlayersSignal = (int)layers.size();
          //      event->AddHelix(move(fittedHelix));
        }
        
        for(auto helix : fittedHelicesBackground){
          averageHelixNpointsBackground += helix.GetNpoints();
          if(helix.GetNpoints() > maxHelixNpointsBackground) maxHelixNpointsBackground = helix.GetNpoints();
          
          unordered_set<int> layers;
          
          for(int iPoint=0; iPoint<helix.GetNpoints(); iPoint++){
            auto point = helix.GetPoints()[iPoint];
            int layer = point->GetLayer();
            if(layer > 0){
              if(iPoint < helix.GetFirstTurningPointIndex()) layers.insert( layer);
              else                                           layers.insert(-layer);
            }
          }
          if(layers.size() > maxHelixNlayersBackground) maxHelixNlayersBackground = (int)layers.size();
        }
        
        averageHelixNpointsSignal /= fittedHelicesSignal.size()>0 ? fittedHelicesSignal.size() : 1;
        nHelixPointsSignal[iTest][iEvent] = averageHelixNpointsSignal;
        nMaxHelixPointsSignal[iTest][iEvent] = maxHelixNpointsSignal;
        nHelixLayersSignal[iTest][iEvent] = maxHelixNlayersSignal;
        
        averageHelixNpointsBackground /= fittedHelicesBackground.size()>0 ? fittedHelicesBackground.size() : 1;
        nHelixPointsBackground[iTest][iEvent] = averageHelixNpointsBackground;
        nMaxHelixPointsBackground[iTest][iEvent] = maxHelixNpointsBackground;
        nHelixLayersBackground[iTest][iEvent] = maxHelixNlayersBackground;
        
        nHelicesSignal->Fill(fittedHelicesSignal.size());
        nHelicesBackground->Fill(fittedHelicesBackground.size());
        
        maxNhitsSignal->Fill(maxHelixNpointsSignal);
        maxNhitsBackground->Fill(maxHelixNpointsBackground);
        
        maxNlayersSignal->Fill(maxHelixNlayersSignal);
        maxNlayersBackground->Fill(maxHelixNlayersBackground);
        
        avgNhitsSignal->Fill(averageHelixNpointsSignal);
        avgNhitsBackground->Fill(averageHelixNpointsBackground);
      }
    }
    nAnalyzedEvents++;
  }
  
  for(int iTest=0; iTest<nTests; iTest++){
    
    vector<double> efficiency;
    vector<double> fakeRate;
    
    vector<double> efficiencyMax;
    vector<double> fakeRateMax;
    
    vector<double> efficiencyLayers;
    vector<double> fakeRateLayers;
    
    int maxNpoints = 20;
    for(int iThreshold=0; iThreshold<maxNpoints; iThreshold++){
      efficiency.push_back(0);
      fakeRate.push_back(0);
      efficiencyMax.push_back(0);
      fakeRateMax.push_back(0);
      efficiencyLayers.push_back(0);
      fakeRateLayers.push_back(0);
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
    
    for(int nLayers : nHelixLayersSignal[iTest]){
      for(int iThreshold=0; iThreshold<maxNpoints; iThreshold++){
        if(nLayers >= iThreshold) efficiencyLayers[iThreshold]+=1;
      }
    }
    for(int nLayers : nHelixLayersBackground[iTest]){
      for(int iThreshold=0; iThreshold<maxNpoints; iThreshold++){
        if(nLayers >= iThreshold) fakeRateLayers[iThreshold]+=1;
      }
    }
    
    rocFun[iTest] = GetRocFunction();
    rocFunMax[iTest] = GetRocFunction();
    rocFunLayers[iTest] = GetRocFunction();
    
    rocGraph[iTest] = new TGraph();
    rocGraphMax[iTest] = new TGraph();
    rocGraphLayers[iTest] = new TGraph();
    
    double paramValue = SetParamValue(iTest);
    
    double maxEff = -inf;
    double maxEffMax = -inf;
    double maxEffLayers = -inf;
    
    double maxEffOverFake = -inf;
    double maxEffMaxOverFake = -inf;
    double maxEffLayersOverFake = -inf;
    
    cout<<"\nFake-eff for test param value: "<<paramValue<<"(test "<<iTest<<")"<<endl;
    
    for(int iThreshold=0; iThreshold<maxNpoints; iThreshold++){
      efficiency[iThreshold] /= nAnalyzedEvents;
      fakeRate[iThreshold] /= nAnalyzedEvents;
      efficiencyMax[iThreshold] /= nAnalyzedEvents;
      fakeRateMax[iThreshold] /= nAnalyzedEvents;
      efficiencyLayers[iThreshold] /= nAnalyzedEvents;
      fakeRateLayers[iThreshold] /= nAnalyzedEvents;
      
      if(efficiency[iThreshold] > maxEff && efficiency[iThreshold] != 1.0)
        maxEff = efficiency[iThreshold];
      if(efficiencyMax[iThreshold] > maxEffMax && efficiencyMax[iThreshold] != 1.0)
        maxEffMax = efficiencyMax[iThreshold];
      if(efficiencyLayers[iThreshold] > maxEffLayers && efficiencyLayers[iThreshold] != 1.0)
        maxEffLayers = efficiencyLayers[iThreshold];
      
      double sigmaApprox = 100*efficiency[iThreshold]/sqrt(100*efficiency[iThreshold]+100000*fakeRate[iThreshold]);
      double sigmaApproxMax = 100*efficiencyMax[iThreshold]/sqrt(100*efficiencyMax[iThreshold]+100000*fakeRateMax[iThreshold]);
      double sigmaApproxLayers = 100*efficiencyLayers[iThreshold]/sqrt(100*efficiencyLayers[iThreshold]+100000*fakeRateLayers[iThreshold]);
      
      
      if(sigmaApprox > maxEffOverFake) maxEffOverFake = sigmaApprox;
      if(sigmaApproxMax > maxEffMaxOverFake) maxEffMaxOverFake = sigmaApproxMax;
      if(sigmaApproxLayers > maxEffLayersOverFake) maxEffLayersOverFake = sigmaApproxLayers;
      
      
      rocGraph[iTest]->SetPoint(iThreshold, fakeRate[iThreshold], efficiency[iThreshold]);
      rocGraphMax[iTest]->SetPoint(iThreshold, fakeRateMax[iThreshold], efficiencyMax[iThreshold]);
      rocGraphLayers[iTest]->SetPoint(iThreshold, fakeRateLayers[iThreshold], efficiencyLayers[iThreshold]);

//      cout<<fakeRate[iThreshold]<<"\t"<<efficiency[iThreshold]<<endl;
//      rocGraph[iTest]->SetPoint(iThreshold, fakeRate[iThreshold], efficiency[iThreshold]);
    }
    cout<<"\n\nFake-eff avg:"<<endl;
    for(int iThreshold=0; iThreshold<maxNpoints; iThreshold++){
      cout<<fakeRate[iThreshold]<<"\t"<<efficiency[iThreshold]<<endl;
    }
    cout<<"\n\nFake-eff max:"<<endl;
    for(int iThreshold=0; iThreshold<maxNpoints; iThreshold++){
      cout<<fakeRateMax[iThreshold]<<"\t"<<efficiencyMax[iThreshold]<<endl;
    }
    cout<<"\n\nFake-eff layers:"<<endl;
    for(int iThreshold=0; iThreshold<maxNpoints; iThreshold++){
      cout<<fakeRateLayers[iThreshold]<<"\t"<<efficiencyLayers[iThreshold]<<endl;
    }
    
    
    rocGraph[iTest]->SetMarkerStyle(20);
    rocGraph[iTest]->SetMarkerSize(2.0);
    rocGraph[iTest]->SetMarkerColor(iTest);
    
    rocGraphMax[iTest]->SetMarkerStyle(20);
    rocGraphMax[iTest]->SetMarkerSize(2.0);
    rocGraphMax[iTest]->SetMarkerColor(iTest);
    
    rocGraphLayers[iTest]->SetMarkerStyle(20);
    rocGraphLayers[iTest]->SetMarkerSize(2.0);
    rocGraphLayers[iTest]->SetMarkerColor(iTest);
    
    canvas->cd(1);
    rocGraph[iTest]->Draw(iTest==0 ? "AP" : "Psame");
    rocGraph[iTest]->Fit(rocFun[iTest]);
    
    canvas->cd(2);
    rocGraphMax[iTest]->Draw(iTest==0 ? "AP" : "Psame");
    rocGraphMax[iTest]->Fit(rocFunMax[iTest]);
    
    canvas->cd(3);
    rocGraphLayers[iTest]->Draw(iTest==0 ? "AP" : "Psame");
    rocGraphLayers[iTest]->Fit(rocFunLayers[iTest]);
    
    double auc = rocFun[iTest]->Integral(0,1);
    double aucMax = rocFunMax[iTest]->Integral(0,1);
    double aucLayers = rocFunLayers[iTest]->Integral(0,1);
    
//    aucAndMaxEffForTest.push_back(make_tuple(paramValue, auc, maxEff));
    aucAndMaxEffForTest.push_back(make_tuple(paramValue, auc, maxEff, maxEffOverFake));
    aucAndMaxEffForTestMax.push_back(make_tuple(paramValue, aucMax, maxEffMax, maxEffMaxOverFake));
    aucAndMaxEffForTestLayers.push_back(make_tuple(paramValue, aucLayers, maxEffLayers, maxEffLayersOverFake));
  }
  
  cout<<"\n\nAUC's for avg:"<<endl;
  for(auto &[param, auc, eff, effOverFake] : aucAndMaxEffForTest){
    cout<<"Param: "<<param<<"\tAUC: "<<auc<<"\tmax eff: "<<eff<<"\teff/fake: "<<effOverFake<<endl;
  }
  cout<<"\n\nAUC's for max:"<<endl;
  for(auto &[param, auc, eff, effOverFake] : aucAndMaxEffForTestMax){
    cout<<"Param: "<<param<<"\tAUC: "<<auc<<"\tmax eff: "<<eff<<"\teff/fake: "<<effOverFake<<endl;
  }
  cout<<"\n\nAUC's for layers:"<<endl;
  for(auto &[param, auc, eff, effOverFake] : aucAndMaxEffForTestLayers){
    cout<<"Param: "<<param<<"\tAUC: "<<auc<<"\tmax eff: "<<eff<<"\teff/fake: "<<effOverFake<<endl;
  }
  
  //  cout<<"helixTagger -- saving events"<<endl;
  //  events.SaveEventsToFiles("afterHelixTagging/");
  //  cout<<"helixTagger -- finished"<<endl;
  
  canvas->cd(4);
  nHelicesSignal->SetFillColorAlpha(kGreen+1, 0.5);
  nHelicesBackground->SetFillColorAlpha(kRed, 0.5);
  nHelicesSignal->Draw();
  nHelicesBackground->Draw("same");
  
  canvas->cd(5);
  maxNhitsSignal->SetFillColorAlpha(kGreen+1, 0.5);
  maxNhitsBackground->SetFillColorAlpha(kRed, 0.5);
  maxNhitsSignal->Draw();
  maxNhitsBackground->Draw("same");
  
  canvas->cd(6);
  avgNhitsSignal->SetFillColorAlpha(kGreen+1, 0.5);
  avgNhitsBackground->SetFillColorAlpha(kRed, 0.5);
  avgNhitsSignal->Draw();
  avgNhitsBackground->Draw("same");
  
  canvas->cd(7);
  maxNlayersSignal->SetFillColorAlpha(kGreen+1, 0.5);
  maxNlayersBackground->SetFillColorAlpha(kRed, 0.5);
  maxNlayersSignal->Draw();
  maxNlayersBackground->Draw("same");
  
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

TF1* GetRocFunction()
{
  TF1 *fun = new TF1(Form("rocFun[%i]",RandInt(0, 999999999)), "[2]*exp(-pow(x-[0],2)/(2*[1]*[1]))/sqrt(2*3.1415*[1]*[1])+[3]",0,1);
  fun->SetParameter(0, 0); // μ
  fun->SetParameter(1, 1); // σ
  fun->SetParameter(2, 1); // a
  fun->SetParameter(3, 1); // b
  
  return fun;
}
