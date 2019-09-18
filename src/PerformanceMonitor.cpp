//  PerformaceMonitor.cpp
//
//  Created by Jeremi Niedziela on 17/06/2019.

#include "PerformanceMonitor.hpp"
#include "Logger.hpp"

PerformanceMonitor::PerformanceMonitor()
{
  
}

PerformanceMonitor::PerformanceMonitor(string _name, string _title,
                                       int _nBins, double min, double max,
                                       EColor color) :
name(_name),
title(_title),
thresholdMin(min),
thresholdMax(max),
nBins(_nBins)
{
  histSignal     = new TH1D(Form("%s%i", name.c_str(), RandInt(0, inf)),
                            title.c_str(),
                            nBins, min, max);
  histBackground = new TH1D(Form("%s%i", name.c_str(), RandInt(0, inf)),
                            title.c_str(),
                            nBins, min, max);
  
  histSignal->SetFillColorAlpha(kGreen+1, 0.3);
  histSignal->SetLineColor(kGreen+1);
  histSignal->SetMarkerColor(kGreen+1);
  
  histBackground->SetFillColorAlpha(kRed, 0.3);
  histBackground->SetLineColor(kRed);
  histBackground->SetMarkerColor(kRed);
  
  rocFun = GetRocFunction();
  rocGraph = new TGraphErrors();
  rocGraph->SetTitle("Looper tagger ROC curves");
  rocGraph->GetXaxis()->SetTitle("Fake rate");
  rocGraph->GetYaxis()->SetTitle("Efficiency");
  rocGraph->SetMarkerStyle(20);
  rocGraph->SetMarkerSize(0.8);
  rocGraph->SetMarkerColor(color);
  rocGraph->SetLineColor(color);
  rocGraph->SetLineStyle(2);
  
  thresholdStep = (max-min)/(nBins+1);
  
  for(int i=0; i<nBins; i++){
    efficiency.push_back(0);
    fakeRate.push_back(0);
  }
  ResetParams();
}

void PerformanceMonitor::ResetParams()
{
  params = {
    {"auc"            , -inf}, // area under ROC curve
    {"sigma_init"     , -inf}, // maximum efficiency lower than 1.0
    {"sigma_L0"       , -inf}, // max significance assuming initial N_sig and N_bck, only when fake rate !=0
    {"sigma_L1"       , -inf}, // max significance assuming L0 N_sig and N_bck, only when fake rate !=0
    {"max_eff"        , -inf}, // max significance assuming L1 N_sig and N_bck, only when fake rate !=0
    {"min_fake"       , -inf}, // 1/fake_rate for the highest efficiency lower than 1.0
    {"max_dist_fake"  , -inf}, // (max) Distance to √c_fake, which is a minimum to be useful
    {"avg_dist_fake"  , -inf}, // (avg) Distance to √c_fake, which is a minimum to be useful
  };
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
  title            = pm.title;
  params           = pm.params;
}

void PerformanceMonitor::SetValue(double value, bool signal)
{
  if(signal){
    valuesSignal.push_back(value);
    histSignal->Fill(value);
  }
  else{
    valuesBackground.push_back(value);
    histBackground->Fill(value);
  }
}

void PerformanceMonitor::CountEventsAboveThreshold()
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
}

void PerformanceMonitor::CalcEfficiency()
{
  CountEventsAboveThreshold();
  ResetParams();
  
  bool first=true;
  
  set<pair<double, double>> rocXY;
  
  for(int iThreshold=0; iThreshold<nBins; iThreshold++){
    double effError = sqrt(1/valuesSignal.size() + 1/efficiency[iThreshold]);
    efficiency[iThreshold]  /= valuesSignal.size();
    effError *= efficiency[iThreshold];
    
    double fakeError = sqrt(1/valuesBackground.size() + 1/fakeRate[iThreshold]);
    fakeRate[iThreshold]    /= valuesBackground.size();
    fakeError *= fakeRate[iThreshold];
    
    if(efficiency[iThreshold] > params["max_eff"] && efficiency[iThreshold] != 1.0){
      params["max_eff"]   = efficiency[iThreshold];
      params["min_fake"]  = 1/fakeRate[iThreshold];
    }
    double sigmaApproxInitial = 3986*efficiency[iThreshold]/sqrt(3986*efficiency[iThreshold]+6E+07*fakeRate[iThreshold]);
    
    if(sigmaApproxInitial > params["sigma_init"] && fakeRate[iThreshold] > 0)
      params["sigma_init"] = sigmaApproxInitial;
    
    double sigmaApproxL0all = 1388*efficiency[iThreshold]/sqrt(1388*efficiency[iThreshold]+1E+05*fakeRate[iThreshold]);
    
    if(sigmaApproxL0all > params["sigma_L0"] && fakeRate[iThreshold] > 0)
      params["sigma_L0"] = sigmaApproxL0all;
    
    double sigmaApproxL1all = 1166*efficiency[iThreshold]/sqrt(1166*efficiency[iThreshold]+6573*fakeRate[iThreshold]);
    
    if(sigmaApproxL1all > params["sigma_L1"] && fakeRate[iThreshold] > 0)
      params["sigma_L1"] = sigmaApproxL1all;
    
    double distToSqrtFake = efficiency[iThreshold]-sqrt(fakeRate[iThreshold]);
    
    if(fakeRate[iThreshold] > 0.01 && fakeRate[iThreshold] < 0.99){
      if(distToSqrtFake > params["max_dist_fake"]) params["max_dist_fake"] = distToSqrtFake;
      if(first){
        params["avg_dist_fake"] = distToSqrtFake;
        first = false;
      }
      else{
        params["avg_dist_fake"] += distToSqrtFake;
      }
    }
    
    rocGraph->SetPoint(iThreshold, fakeRate[iThreshold], efficiency[iThreshold]);
    rocGraph->SetPointError(iThreshold, fakeError, effError);
    rocXY.insert(make_pair(fakeRate[iThreshold], efficiency[iThreshold]));
  }
  params["avg_dist_fake"] /= nBins;
  
//  rocGraph->Fit(rocFun,"Q");
  params["auc"] = rocFun->Integral(0,1);
}

double PerformanceMonitor::GetMaxDistanceFromSqrtFake()
{
  double maxDistance = -inf;
  
  for(int iThreshold=0; iThreshold<efficiency.size(); iThreshold++){
    if(efficiency[iThreshold]==0 || efficiency[iThreshold]==1) continue;
    double distance = efficiency[iThreshold] - sqrt(fakeRate[iThreshold]);
    if(distance > maxDistance) maxDistance = distance;
  }
  return maxDistance;
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
  for(auto &[name, value] : params) Log(1)<<name<<": "<<value<<"\t";
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

void PerformanceMonitor::DrawRocGraph(bool first, TLegend *legend)
{
  rocGraph->Draw(first ? "APLE" : "PLEsame");
  if(first){
    rocGraph->GetYaxis()->SetRangeUser(0,1);
    rocGraph->GetXaxis()->SetLimits(0,1);
    
    TLine *line = new TLine(0,0,1,1);
    line->SetLineColor(kBlack);
    line->Draw("same");
    
    TF1 *sqrtFun = new TF1("sqrt(fake)","sqrt(x)",0,1);
    sqrtFun->SetLineColor(kRed);
    sqrtFun->SetLineStyle(1);
    sqrtFun->SetLineWidth(2.0);
    sqrtFun->Draw("same");
    legend->AddEntry(sqrtFun, "c_{eff} = #sqrt{c_{fake}}", "L");
  }
  if(legend){
    legend->AddEntry(rocGraph, title.c_str(), "PL");
  }
}

void PerformanceMonitor::DrawHists(TLegend *legend)
{
  histSignal->Sumw2();
  histSignal->Scale(1/histSignal->GetEntries());
  histSignal->Draw();
  histSignal->SetMaximum(1.0);
  histBackground->Sumw2();
  histBackground->DrawNormalized("same");
  
  if(legend){
    legend->AddEntry(histSignal     , "Signal"      , "PE");
    legend->AddEntry(histBackground , "Background"  , "PE");
  }
  
}
