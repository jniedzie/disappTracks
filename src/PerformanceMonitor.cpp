//  PerformaceMonitor.cpp
//
//  Created by Jeremi Niedziela on 17/06/2019.

#include "PerformanceMonitor.hpp"

PerformanceMonitor::PerformanceMonitor() :
nTests(0)
{
  
}

PerformanceMonitor::PerformanceMonitor(string _name, int _nBins, double min, double max,
                                       int _nTests, int nEvents) :
name(_name),
nTests(_nTests),
thresholdMin(min),
thresholdMax(max),
nBins(_nBins)
{
  histSignal     = new TH1D(Form("%s", name.c_str()), Form("%s", name.c_str()), nBins, min, max);
  histBackground = new TH1D(Form("%s bck", name.c_str()), Form("%s bck", name.c_str()), nBins, min, max);
  
  histSignal->SetFillColorAlpha(kGreen+1, 0.3);
  histBackground->SetFillColorAlpha(kRed, 0.3);
  
  for(int iTest=0; iTest<nTests; iTest++){
    vector<double> vec;
    
    for(int iEvent=0; iEvent<nEvents; iEvent++) vec.push_back(0.0);
    
    valuesSignal.push_back(vec);
    valuesBackground.push_back(vec);
    
    rocFun.push_back(nullptr);
    rocGraph.push_back(nullptr);
  }
  thresholdStep = (max-min)/nBins;
  
  for(int i=0; i<nBins; i++){
    efficiency.push_back(0);
    fakeRate.push_back(0);
  }
  
}

void PerformanceMonitor::operator=(const PerformanceMonitor &pm)
{
  nTests           = pm.nTests;
  thresholdMin     = pm.thresholdMin;
  thresholdMax     = pm.thresholdMax;
  thresholdStep    = pm.thresholdStep;
  nBins            = pm.nBins;
  valuesSignal     = pm.valuesSignal;
  valuesBackground = pm.valuesBackground;
  efficiency       = pm.efficiency;
  fakeRate         = pm.fakeRate;
  histSignal       = pm.histSignal;
  histBackground   = pm.histBackground;
  rocFun           = pm.rocFun;
  rocGraph         = pm.rocGraph;
  name             = pm.name;
}

void PerformanceMonitor::SetValues(int iTest, int iEvent, double valueSignal, double valueBackground)
{
  valuesSignal[iTest][iEvent] = valueSignal;
  valuesBackground[iTest][iEvent] = valueBackground;
  histSignal->Fill(valueSignal);
  histBackground->Fill(valueBackground);
}

void PerformanceMonitor::CalcEfficiency(function<double(int)> SetParamValue, int nAnalyzedEvents)
{
  cout<<"\n\n========================================================"<<endl;
  cout<<"Performance monitor: "<<name<<endl;
  
  vector<double> auc;
  vector<double> maxEfficiency;
  vector<double> significanceInitial;
  vector<double> significanceAfterL0all;
  
  for(int iTest=0; iTest<nTests; iTest++){
    
    for(int iThreshold=0; iThreshold<nBins; iThreshold++){
      efficiency[iThreshold] = 0.0;
      fakeRate[iThreshold]   = 0.0;
      
      double threshold = thresholdMin + iThreshold*thresholdStep;
      for(double nPoints : valuesSignal[iTest]){
        if(nPoints >= threshold) efficiency[iThreshold]+=1;
      }
      for(double nPoints : valuesBackground[iTest]){
        if(nPoints >= threshold) fakeRate[iThreshold]+=1;
      }
    }
    
    rocFun[iTest] = GetRocFunction();
    rocGraph[iTest] = new TGraph();
    rocGraph[iTest]->SetTitle(Form("%s %i", name.c_str(), iTest));
    rocGraph[iTest]->SetMarkerStyle(20);
    rocGraph[iTest]->SetMarkerSize(1.0);
    rocGraph[iTest]->SetMarkerColor(iTest+1);
    
    double paramValue = SetParamValue(iTest);
    
    double maxEff = -inf;
    double maxInitial = -inf;
    double maxL0all = -inf;
    
    cout<<"\nFake-eff for test param value: "<<paramValue<<"(test "<<iTest<<")"<<endl;
    
    
    for(int iThreshold=0; iThreshold<nBins; iThreshold++){
      efficiency[iThreshold]  /= nAnalyzedEvents;
      fakeRate[iThreshold]    /= nAnalyzedEvents;
      
      if(efficiency[iThreshold] > maxEff && efficiency[iThreshold] != 1.0) maxEff = efficiency[iThreshold];
      
      double sigmaApproxInitial = 3986*efficiency[iThreshold]/sqrt(3986*efficiency[iThreshold]+6E+07*fakeRate[iThreshold]);
      if(sigmaApproxInitial > maxInitial && fakeRate[iThreshold] > 0) maxInitial = sigmaApproxInitial;
      
      double sigmaApproxL0all = 1388*efficiency[iThreshold]/sqrt(1388*efficiency[iThreshold]+1E+05*fakeRate[iThreshold]);
      if(sigmaApproxL0all > maxL0all && fakeRate[iThreshold] > 0) maxL0all = sigmaApproxL0all;
      
      rocGraph[iTest]->SetPoint(iThreshold, fakeRate[iThreshold], efficiency[iThreshold]);
    }
    
    cout<<"\n\nFake-eff avg:"<<endl;
    
    for(int iThreshold=0; iThreshold<nBins; iThreshold++){
      cout<<fakeRate[iThreshold]<<"\t"<<efficiency[iThreshold]<<endl;
    }
    
    rocGraph[iTest]->Fit(rocFun[iTest],"Q");
    
    auc.push_back(rocFun[iTest]->Integral(0,1));
    maxEfficiency.push_back(maxEff);
    significanceInitial.push_back(maxInitial);
    significanceAfterL0all.push_back(maxL0all);
  }
  
  for(int iTest=0; iTest<nTests; iTest++){
    cout<<"Param: "<<SetParamValue(iTest)<<"\tAUC: "<<auc[iTest]<<"\tmax eff: "<<maxEfficiency[iTest]<<"\tsign Initial: "<<significanceInitial[iTest]<<"\tsign L0 all: "<<significanceAfterL0all[iTest]<<endl;
  }
}

TF1* PerformanceMonitor::GetRocFunction()
{
  TF1 *fun = new TF1(Form("rocFun[%i]",RandInt(0, 999999999)), "[2]*exp(-pow(x-[0],2)/(2*[1]*[1]))/sqrt(2*3.1415*[1]*[1])+[3]",0,1);
  fun->SetParameter(0, 0); // μ
  fun->SetParameter(1, 1); // σ
  fun->SetParameter(2, 1); // a
  fun->SetParameter(3, 1); // b
  
  return fun;
}

void PerformanceMonitor::DrawRocGraphs()
{
  for(int iTest=0; iTest<nTests; iTest++){
    rocGraph[iTest]->Draw(iTest==0 ? "AP" : "Psame");
  }
}

void PerformanceMonitor::DrawHists()
{
  histSignal->Draw();
  histBackground->Draw("same");
}
