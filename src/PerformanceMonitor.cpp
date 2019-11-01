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
                                       EColor color, bool alternativeColors,
                                       int _fixUpperThresholdBin) :
name(_name),
title(_title),
thresholdMin(min),
thresholdMax(max),
nBins(_nBins),
fixedUpperThresholdBin(_fixUpperThresholdBin)
{
  histSignal     = new TH1D(Form("%s%i", name.c_str(), RandInt(0, inf)),
                            title.c_str(),
                            nBins, min, max);
  histBackground = new TH1D(Form("%s%i", name.c_str(), RandInt(0, inf)),
                            title.c_str(),
                            nBins, min, max);
  
  histSignal->SetFillColorAlpha(alternativeColors ? kBlue : kGreen+1, 0.3);
  histSignal->SetLineColor(alternativeColors ? kBlue : kGreen+1);
  histSignal->SetMarkerColor(alternativeColors ? kBlue : kGreen+1);
  
  histBackground->SetFillColorAlpha(alternativeColors ? kOrange+1 : kRed, 0.3);
  histBackground->SetLineColor(alternativeColors ? kOrange+1 : kRed);
  histBackground->SetMarkerColor(alternativeColors ? kOrange+1 : kRed);
  
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
  
  thresholdStep = (max-min)/(nBins);
  
  for(int i=0; i<nBins; i++){
    vector<double> tmp(nBins, 0);
    for(int j=0; j<nBins; j++){
      efficiency.push_back(tmp);
      fakeRate.push_back(tmp);
    }
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
  fixedUpperThresholdBin = pm.fixedUpperThresholdBin;
  
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
  for(int iThresholdLow=0; iThresholdLow<efficiency.size(); iThresholdLow++){
    int thresholdUpMin = iThresholdLow+1;
    int thresholdUpMax = (int)efficiency[0].size();
    if(fixedUpperThresholdBin > 0){
      thresholdUpMin = fixedUpperThresholdBin;
      thresholdUpMax = fixedUpperThresholdBin+1;
    }
    for(int iThresholdUp=thresholdUpMin; iThresholdUp<thresholdUpMax; iThresholdUp++){
      efficiency[iThresholdLow][iThresholdUp] = 0.0;
      fakeRate[iThresholdLow][iThresholdUp]   = 0.0;
    
      double thresholdLow = thresholdMin + iThresholdLow*thresholdStep;
      double thresholdUp  = thresholdMin + iThresholdUp*thresholdStep;
      for(double nPoints : valuesSignal){
        if(nPoints >= thresholdLow && nPoints <= thresholdUp) efficiency[iThresholdLow][iThresholdUp]+=1;
      }
      for(double nPoints : valuesBackground){
        if(nPoints >= thresholdLow && nPoints <= thresholdUp) fakeRate[iThresholdLow][iThresholdUp]+=1;
      }
    }
  }
}

void PerformanceMonitor::CalcEfficiency()
{
  CountEventsAboveThreshold();
  ResetParams();
  
  bool first=true;
  
  set<pair<double, double>> rocXY;
  
  double bestEff=-inf;
  double bestFake=inf;
  
  int iPoint=0;
  
  for(int iThresholdLow=0; iThresholdLow<efficiency.size(); iThresholdLow++){
    int thresholdUpMin = iThresholdLow+1;
    int thresholdUpMax = (int)efficiency[0].size();
    if(fixedUpperThresholdBin > 0){
      thresholdUpMin = fixedUpperThresholdBin;
      thresholdUpMax = fixedUpperThresholdBin+1;
    }
    for(int iThresholdUp=thresholdUpMin; iThresholdUp<thresholdUpMax; iThresholdUp++){
      double effError;
      if(valuesSignal.size()==0 || efficiency[iThresholdLow][iThresholdUp]==0){
        effError = 0;
      }
      else{
        effError = sqrt(1/valuesSignal.size() + 1/efficiency[iThresholdLow][iThresholdUp]);
        efficiency[iThresholdLow][iThresholdUp]  /= valuesSignal.size();
        effError *= efficiency[iThresholdLow][iThresholdUp];
      }
      
      double fakeError;
      if(valuesBackground.size()==0 || fakeRate[iThresholdLow][iThresholdUp]==0){
        fakeError=0;
      }
      else{
        fakeError = sqrt(1/valuesBackground.size() + 1/fakeRate[iThresholdLow][iThresholdUp]);
        fakeRate[iThresholdLow][iThresholdUp]    /= valuesBackground.size();
        fakeError *= fakeRate[iThresholdLow][iThresholdUp];
      }
      
      if(efficiency[iThresholdLow][iThresholdUp] > params["max_eff"] && efficiency[iThresholdLow][iThresholdUp] != 1.0){
        params["max_eff"]   = efficiency[iThresholdLow][iThresholdUp];
        params["min_fake"]  = 1/fakeRate[iThresholdLow][iThresholdUp];
      }
      double sigmaApproxInitial = 3986*efficiency[iThresholdLow][iThresholdUp]/sqrt(3986*efficiency[iThresholdLow][iThresholdUp]+6E+07*fakeRate[iThresholdLow][iThresholdUp]);
      
      if(sigmaApproxInitial > params["sigma_init"] && fakeRate[iThresholdLow][iThresholdUp] > 0)
        params["sigma_init"] = sigmaApproxInitial;
      
      double sigmaApproxL0all = 1388*efficiency[iThresholdLow][iThresholdUp]/sqrt(1388*efficiency[iThresholdLow][iThresholdUp]+1E+05*fakeRate[iThresholdLow][iThresholdUp]);
      
      if(sigmaApproxL0all > params["sigma_L0"] && fakeRate[iThresholdLow][iThresholdUp] > 0)
        params["sigma_L0"] = sigmaApproxL0all;
      
      double sigmaApproxL1all = 1166*efficiency[iThresholdLow][iThresholdUp]/sqrt(1166*efficiency[iThresholdLow][iThresholdUp]+6573*fakeRate[iThresholdLow][iThresholdUp]);
      
      if(sigmaApproxL1all > params["sigma_L1"] && fakeRate[iThresholdLow][iThresholdUp] > 0)
        params["sigma_L1"] = sigmaApproxL1all;
      
      double distToSqrtFake = efficiency[iThresholdLow][iThresholdUp]-sqrt(fakeRate[iThresholdLow][iThresholdUp]);
      
      if(fakeRate[iThresholdLow][iThresholdUp] > 0.01 && fakeRate[iThresholdLow][iThresholdUp] < 0.99){
        if(distToSqrtFake > params["max_dist_fake"]){
          params["max_dist_fake"] = distToSqrtFake;
          bestEff = efficiency[iThresholdLow][iThresholdUp];
          bestFake = fakeRate[iThresholdLow][iThresholdUp];
        }
        if(first){
          params["avg_dist_fake"] = distToSqrtFake;
          first = false;
        }
        else{
          params["avg_dist_fake"] += distToSqrtFake;
        }
      }
      
      rocGraph->SetPoint(iPoint, fakeRate[iThresholdLow][iThresholdUp], efficiency[iThresholdLow][iThresholdUp]);
      rocGraph->SetPointError(iPoint, fakeError, effError);
      rocXY.insert(make_pair(fakeRate[iThresholdLow][iThresholdUp], efficiency[iThresholdLow][iThresholdUp]));
      
      iPoint++;
    }
  }
  params["avg_dist_fake"] /= nBins;
  
//  cout<<"Distance: "<<params["max_dist_fake"]<<"\teff: "<<bestEff<<"\tfake: "<<bestFake<<endl;
  
  //  rocGraph->Fit(rocFun,"Q");
  params["auc"] = rocFun->Integral(0,1);
}
  
double PerformanceMonitor::GetMaxDistanceFromSqrtFake(double &bestEff, double &bestFake, int &thresholdLowBin, int &thresholdUpBin)
{
  double maxDistance = -inf;
  
  for(int iThresholdLow=0; iThresholdLow<efficiency.size(); iThresholdLow++){
    int thresholdUpMin = iThresholdLow+1;
    int thresholdUpMax = (int)efficiency[0].size();
    if(fixedUpperThresholdBin > 0){
      thresholdUpMin = fixedUpperThresholdBin;
      thresholdUpMax = fixedUpperThresholdBin+1;
    }
    for(int iThresholdUp=thresholdUpMin; iThresholdUp<thresholdUpMax; iThresholdUp++){
      
      if(fakeRate[iThresholdLow][iThresholdUp] < 0.01 || fakeRate[iThresholdLow][iThresholdUp] > 0.99) continue;
      
      double distance = efficiency[iThresholdLow][iThresholdUp] - sqrt(fakeRate[iThresholdLow][iThresholdUp]);
      if(distance > maxDistance){
        maxDistance = distance;
        bestEff = efficiency[iThresholdLow][iThresholdUp];
        bestFake = fakeRate[iThresholdLow][iThresholdUp];
        thresholdLowBin = iThresholdLow;
        thresholdUpBin  = iThresholdUp;
      }
    }
  }
  
  return maxDistance;
}

void PerformanceMonitor::PrintFakesEfficiency()
{
  cout<<"Fake-eff avg:"<<endl;
  for(int iThresholdLow=0; iThresholdLow<efficiency.size(); iThresholdLow++){
    int thresholdUpMin = iThresholdLow+1;
    int thresholdUpMax = (int)efficiency[0].size();
    if(fixedUpperThresholdBin > 0){
      thresholdUpMin = fixedUpperThresholdBin;
      thresholdUpMax = fixedUpperThresholdBin+1;
    }
    for(int iThresholdUp=thresholdUpMin; iThresholdUp<thresholdUpMax; iThresholdUp++){
      
      if(fakeRate[iThresholdLow][iThresholdUp] != 0 && efficiency[iThresholdLow][iThresholdUp] != 1){
        double thresholdLow = thresholdMin + iThresholdLow*thresholdStep;
        double thresholdUp  = thresholdMin + iThresholdUp*thresholdStep;
        
        cout<<"Thresholds: "<<thresholdLow<<" - "<<thresholdUp<<"\t";
        cout<<fakeRate[iThresholdLow][iThresholdUp]<<"\t"<<efficiency[iThresholdLow][iThresholdUp]<<endl;
      }
    }
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
    if(legend) legend->AddEntry(sqrtFun, "c_{eff} = #sqrt{c_{fake}}", "L");
  }
  if(legend){
    legend->AddEntry(rocGraph, title.c_str(), "PL");
  }
}

void PerformanceMonitor::DrawHists(bool first, TLegend *legend)
{
  histSignal->Sumw2();
  histSignal->Scale(1/histSignal->GetEntries());
  histSignal->Draw(first ? "" : "same");
  histSignal->SetMaximum(1.0);
  histBackground->Sumw2();
  histBackground->DrawNormalized("same");
  
  if(legend){
    legend->AddEntry(histSignal     , "Signal"      , "PE");
    legend->AddEntry(histBackground , "Background"  , "PE");
  }
  
}
