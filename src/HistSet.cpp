//
//  HistSet.cpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#include "HistSet.hpp"
#include "Helpers.hpp"
#include "Track.hpp"

#include <TLegend.h>

HistSet::HistSet() :
signal(nullptr),
background(nullptr),
data(nullptr)
{
  
}

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

void HistSet::FillFromEvents(Events *signalEvents, Events *backgroundEvents, Events *dataEvents, EVar var)
{
  if(var == kDedx || var == kSizeX || var == kSizeY){
    FillFromEventsPerLayer(signalEvents, backgroundEvents, dataEvents, var);
  }
  else{
    FillFromEventsGlobal(signalEvents, backgroundEvents, dataEvents, var);
  }
}

void HistSet::FillFromEventsPerLayer(Events *signalEvents, Events *backgroundEvents, Events *dataEvents, EVar var)
{
  const char* title= GetTitle(var);
  int nBins=0;
  double min=0.0, max=10.0;
  
  if(var == kDedx){
    nBins = 50;
  }
  if(var == kSizeX || var == kSizeY){
    nBins = 10;
  }
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    TH1D *histSignal = new TH1D(Form("%s_layer[%i]_signal",title,iLayer),Form("%s_layer[%i]_signal",title, iLayer),nBins,min,max);
    TH1D *histBackground = new TH1D(Form("%s_layer[%i]_background",title,iLayer),Form("%s_layer[%i]_background",title,iLayer),nBins,min,max);
    TH1D *histData = new TH1D(Form("%s_layer[%i]_data",title,iLayer),Form("%s_layer[%i]_data",title,iLayer),nBins,min,max);
    signalPerLayer.push_back(histSignal);
    backgroundPerLayer.push_back(histBackground);
    dataPerLayer.push_back(histData);
  }
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    
    if(signalEvents){// fill signal events
      for(int iEvent=0;iEvent<signalEvents->size();iEvent++){
        for(int iTrack=0;iTrack<signalEvents->At(iEvent)->GetNtracks();iTrack++){
          Track *track = signalEvents->At(iEvent)->GetTrack(iTrack);
          double value = 0.0;
          
          if(var == kDedx)        value = track->GetDeDxInLayer(iLayer);
          else if(var == kSizeX)  value = track->GetSizeXinLayer(iLayer);
          else if(var == kSizeY)  value = track->GetSizeYinLayer(iLayer);
          
          if(value > 0.0001) signalPerLayer[iLayer]->Fill(value);
        }
      }
    }
    if(backgroundEvents){// fill background events
      for(int iEvent=0;iEvent<backgroundEvents->size();iEvent++){
        for(int iTrack=0;iTrack<backgroundEvents->At(iEvent)->GetNtracks();iTrack++){
          Track *track = backgroundEvents->At(iEvent)->GetTrack(iTrack);
          
          double value = 0.0;
          
          if(var == kDedx)        value = track->GetDeDxInLayer(iLayer);
          else if(var == kSizeX)  value = track->GetSizeXinLayer(iLayer);
          else if(var == kSizeY)  value = track->GetSizeYinLayer(iLayer);
          
          if(value > 0.0001) backgroundPerLayer[iLayer]->Fill(value);
        }
      }
    }
    if(dataEvents){// fill data events
      for(int iEvent=0;iEvent<dataEvents->size();iEvent++){
        for(int iTrack=0;iTrack<dataEvents->At(iEvent)->GetNtracks();iTrack++){
          Track *track = dataEvents->At(iEvent)->GetTrack(iTrack);
          double value = 0.0;
          
          if(var == kDedx)        value = track->GetDeDxInLayer(iLayer);
          else if(var == kSizeX)  value = track->GetSizeXinLayer(iLayer);
          else if(var == kSizeY)  value = track->GetSizeYinLayer(iLayer);
          
          if(value > 0.0001) dataPerLayer[iLayer]->Fill(value);
        }
      }
    }
  }
}

void HistSet::FillFromEventsGlobal(Events *signalEvents, Events *backgroundEvents, Events *dataEvents, EVar var)
{
  
  // fill fomr signal events
  for(int iEvent=0;iEvent<signalEvents->size();iEvent++){
    Event *event = signalEvents->At(iEvent);
    
    for(int iJet=0;iJet<event->GetNjets();iJet++){
      Jet *jet = event->GetJet(iJet);
      
      double value=0.0;
      
      if(var == kJetPt) value = jet->GetPt();
      if(var == kJetEta) value = jet->GetEta();
      if(var == kJetPhi) value = jet->GetPhi();
        
      signal->Fill(value);
    }
  }
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

void HistSet::DrawPerLayer(EVar var)
{
  TLegend *leg = GetLegend(0.15,0.25,0.75,0.5,"Data type");
  
  TCanvas *c1 = new TCanvas(GetTitle(var),GetTitle(var),2880,1800);
  
  if(var == kDedx) c1->Divide(4,4);
  if(var == kSizeX || var == kSizeY) c1->Divide(2,2);

  for(int iLayer=0;iLayer<nLayers;iLayer++){
    if(iLayer > 3 && (var==kSizeX || var==kSizeY)) continue;
    
    c1->cd(iLayer+1);
    
    signalPerLayer[iLayer]->SetLineColor(kGreen+1);
    signalPerLayer[iLayer]->Draw();
    
    backgroundPerLayer[iLayer]->SetLineColor(kRed+1);
    backgroundPerLayer[iLayer]->Draw("same");
  }
  
  leg->AddEntry(signalPerLayer[0],"Signal","lp");
  leg->AddEntry(backgroundPerLayer[0],"Background","lp");
  leg->Draw();
  
  c1->Update();
}

TLegend* HistSet::GetLegend(double legendW, double legendH, double legendX, double legendY,const char* header)
{
  TLegend *leg = new TLegend(legendX,legendY,legendX+legendW,legendY+legendH);
  leg->SetHeader(header);
  return leg;
}


const char* HistSet::GetTitle(EVar var)
{
  if(var == kDedx) return "dedx";
  if(var == kSizeX) return "sizeX";
  if(var == kSizeY) return "sizeY";
  
  return "";
}

