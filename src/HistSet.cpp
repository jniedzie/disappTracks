//
//  HistSet.cpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "HistSet.hpp"

#include <TLegend.h>

HistSet::HistSet(const char* title, int nBins, double min, double max) :
signal(nullptr),
background(nullptr),
data(nullptr)
{
  signal     = new TH1D(Form("%s (signal)",title),Form("%s (signal)",title),nBins,min,max);
  background = new TH1D(Form("%s (background)",title),Form("%s (background)",title),nBins,min,max);
  data       = new TH1D(Form("%s (data)",title),Form("%s (data)",title),nBins,min,max);
  
  signal->Sumw2();
  background->Sumw2();
  data->Sumw2();
  
}

HistSet::~HistSet()
{
  
  
}

void HistSet::Draw(TCanvas *c1, int pad)
{
  TLegend *leg = GetLegend(0.15,0.25,0.75,0.5,"Data type");
  
  double maxYvalue = -9999999;
  double minYvalue =  9999999;
  
  for(int i=0;i<signal->GetNbinsX();i++){
    if(signal->GetBinContent(i) > maxYvalue) maxYvalue = signal->GetBinContent(i);
    if(signal->GetBinContent(i) < minYvalue) minYvalue = signal->GetBinContent(i);
  }
  for(int i=0;i<background->GetNbinsX();i++){
    if(background->GetBinContent(i) > maxYvalue) maxYvalue = background->GetBinContent(i);
    if(background->GetBinContent(i) < minYvalue) minYvalue = background->GetBinContent(i);
  }
  if(minYvalue < 0) minYvalue *= 1.2;
  else              minYvalue *= 0.8;
  if(maxYvalue < 0) minYvalue *= 0.8;
  else              minYvalue *= 1.2;
  
  
  c1->cd(pad);
  signal->Draw();
  background->SetLineColor(kRed);
  background->Draw("same");
  
  signal->GetYaxis()->SetRangeUser(minYvalue, maxYvalue);
  
  leg->AddEntry(signal,"Signal","lp");
  leg->AddEntry(background,"Background","lp");
  leg->Draw();
  
  c1->Update();
}


TLegend* HistSet::GetLegend(double legendW, double legendH, double legendX, double legendY,const char* header)
{
  TLegend *leg = new TLegend(legendX,legendY,legendX+legendW,legendY+legendH);
  leg->SetHeader(header);
  return leg;
}
