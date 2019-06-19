//  PerformaceMonitor.cpp
//
//  Created by Jeremi Niedziela on 17/06/2019.

#include "PerformanceMonitor.hpp"

PerformanceMonitor::PerformanceMonitor()
{
  
}

PerformanceMonitor::PerformanceMonitor(string _name, int _nBins, double min, double max, int nEvents) :
name(_name),
thresholdMin(min),
thresholdMax(max),
nBins(_nBins),
auc(0),
maxEfficiency(0),
significanceInitial(0),
significanceAfterL0(0)
{
  histSignal     = new TH1D(Form("%s%i", name.c_str(), RandInt(0, inf)),
                            Form("%s%i", name.c_str(), RandInt(0, inf)),
                            nBins, min, max);
  histBackground = new TH1D(Form("%s%i", name.c_str(), RandInt(0, inf)),
                            Form("%s%i", name.c_str(), RandInt(0, inf)),
                            nBins, min, max);
  
  histSignal->SetFillColorAlpha(kGreen+1, 0.3);
  histBackground->SetFillColorAlpha(kRed, 0.3);
  
  for(int iEvent=0; iEvent<nEvents; iEvent++){
    valuesSignal.push_back(0.0);
    valuesBackground.push_back(0.0);
  }
  
  rocFun = GetRocFunction();
  rocGraph = new TGraph();
  rocGraph->SetTitle(name.c_str());
  rocGraph->SetMarkerStyle(20);
  rocGraph->SetMarkerSize(1.0);
  rocGraph->SetMarkerColor(kRed);
  
  thresholdStep = (max-min)/(nBins+1);
  
  for(int i=0; i<nBins; i++){
    efficiency.push_back(0);
    fakeRate.push_back(0);
  }
}

void PerformanceMonitor::operator=(const PerformanceMonitor &pm)
{
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

void PerformanceMonitor::SetValues(int iEvent, double valueSignal, double valueBackground)
{
  valuesSignal[iEvent] = valueSignal;
  valuesBackground[iEvent] = valueBackground;
  histSignal->Fill(valueSignal);
  histBackground->Fill(valueBackground);
}

void PerformanceMonitor::CalcEfficiency(int nAnalyzedEvents)
{
  for(int iThreshold=0; iThreshold<nBins; iThreshold++){
    efficiency[iThreshold] = 0.0;
    fakeRate[iThreshold]   = 0.0;
    
    double threshold = thresholdMin + iThreshold*thresholdStep;
    for(double nPoints : valuesSignal){
      if(nPoints >= threshold) efficiency[iThreshold]+=1;
    }
    for(double nPoints : valuesBackground){
      if(nPoints >= threshold) fakeRate[iThreshold]+=1;
    }
  }
  
  maxEfficiency = -inf;
  significanceInitial = -inf;
  significanceAfterL0 = -inf;
  significanceAfterL1 = -inf;
  
  set<pair<double, double>> rocXY;
  
  for(int iThreshold=0; iThreshold<nBins; iThreshold++){
    efficiency[iThreshold]  /= nAnalyzedEvents;
    fakeRate[iThreshold]    /= nAnalyzedEvents;
    
    if(efficiency[iThreshold] > maxEfficiency && efficiency[iThreshold] != 1.0)
      maxEfficiency = efficiency[iThreshold];
    
    double sigmaApproxInitial = 3986*efficiency[iThreshold]/sqrt(3986*efficiency[iThreshold]+6E+07*fakeRate[iThreshold]);
    
    if(sigmaApproxInitial > significanceInitial && fakeRate[iThreshold] > 0)
      significanceInitial = sigmaApproxInitial;
    
    double sigmaApproxL0all = 1388*efficiency[iThreshold]/sqrt(1388*efficiency[iThreshold]+1E+05*fakeRate[iThreshold]);
    
    if(sigmaApproxL0all > significanceAfterL0 && fakeRate[iThreshold] > 0)
      significanceAfterL0 = sigmaApproxL0all;
    
    double sigmaApproxL1all = 1166*efficiency[iThreshold]/sqrt(1166*efficiency[iThreshold]+6573*fakeRate[iThreshold]);
    
    if(sigmaApproxL1all > significanceAfterL1 && fakeRate[iThreshold] > 0)
      significanceAfterL1 = sigmaApproxL1all;
    
    rocGraph->SetPoint(iThreshold, fakeRate[iThreshold], efficiency[iThreshold]);
    rocXY.insert(make_pair(fakeRate[iThreshold], efficiency[iThreshold]));
  }
//  rocGraph->Fit(rocFun,"Q");
  
  double aucXY=0;

  for(auto &xy : rocXY){
    aucXY += xy.second-xy.first;
  }
  
//  auc = rocFun->Integral(0,1);
  auc = aucXY;
}


void PerformanceMonitor::PrintFakesEfficiency()
{
  cout<<"Fake-eff avg:"<<endl;
  for(int iThreshold=0; iThreshold<nBins; iThreshold++){
    cout<<fakeRate[iThreshold]<<"\t"<<efficiency[iThreshold]<<endl;
  }
  cout<<endl;
}

void PerformanceMonitor::PrintParams()
{
  cout<<"AUC: "<<auc<<"\tmax eff: "<<maxEfficiency<<"\t";
  cout<<"sign Initial: "<<significanceInitial<<"\t";
  cout<<"sign L0 all: "<<significanceAfterL0<<endl;
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

void PerformanceMonitor::DrawRocGraph(bool first)
{
  rocGraph->Draw(first ? "AP" : "Psame");
}

void PerformanceMonitor::DrawHists()
{
  histSignal->Draw();
  histBackground->Draw("same");
}
