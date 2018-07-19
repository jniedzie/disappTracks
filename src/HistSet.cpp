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
#include <THStack.h>

HistSet::HistSet(const char* title, int nBins, double min, double max) :
signal(nullptr),
background(nullptr),
data(nullptr),
var(kCustom),
customTitle(title)
{
  signal     = new TH1D(Form("%s (signal)",title),Form("%s (signal)",title),nBins,min,max);
  background = new TH1D(Form("%s (background)",title),Form("%s (background)",title),nBins,min,max);
  data       = new TH1D(Form("%s (data)",title),Form("%s (data)",title),nBins,min,max);
  
  signal->Sumw2(DoSumw2());
  background->Sumw2(DoSumw2());
  data->Sumw2(DoSumw2());
}

HistSet::HistSet(EVar _var) :
signal(nullptr),
background(nullptr),
data(nullptr),
var(_var)
{
  const char* title = GetTitle();
  int nBins = GetNbins();
  double min = GetMin();
  double max = GetMax();
  
  signal     = new TH1D(Form("%s (signal)",title),Form("%s (signal)",title),nBins,min,max);
  background = new TH1D(Form("%s (background)",title),Form("%s (background)",title),nBins,min,max);
  data       = new TH1D(Form("%s (data)",title),Form("%s (data)",title),nBins,min,max);
  
  signal->Sumw2(DoSumw2());
  background->Sumw2(DoSumw2());
  data->Sumw2(DoSumw2());
}

HistSet::~HistSet()
{
  
  
}

void HistSet::FillFromEvents(Events *signalEvents, Events *backgroundEvents, Events *dataEvents)
{
  if(var == kDedx || var == kSizeX || var == kSizeY){
    FillFromEventsPerLayer(signalEvents, backgroundEvents, dataEvents);
  }
  else{
    FillFromEventsGlobal(signalEvents, backgroundEvents, dataEvents);
  }
}

void HistSet::FillFromEventsPerLayer(Events *signalEvents, Events *backgroundEvents, Events *dataEvents)
{
  const char* title= GetTitle();
  int nBins = GetNbins();
  double min = GetMin();
  double max = GetMax();
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    TH1D *histSignal = new TH1D(Form("%s_layer[%i]_signal",title,iLayer),Form("%s_layer[%i]_signal",title, iLayer),nBins,min,max);
    
    Fill(histSignal,signalEvents,iLayer);
    signalPerLayer.push_back(histSignal);
    
    TH1D *histBackground = new TH1D(Form("%s_layer[%i]_background",title,iLayer),Form("%s_layer[%i]_background",title,iLayer),nBins,min,max);
    Fill(histBackground,backgroundEvents,iLayer);
    backgroundPerLayer.push_back(histBackground);
    
    TH1D *histData = new TH1D(Form("%s_layer[%i]_data",title,iLayer),Form("%s_layer[%i]_data",title,iLayer),nBins,min,max);
    Fill(histData,dataEvents,iLayer);
    dataPerLayer.push_back(histData);
  }
}

void HistSet::FillFromEventsGlobal(Events *signalEvents, Events *backgroundEvents, Events *dataEvents)
{
  Fill(signal, signalEvents, var);
  Fill(background, backgroundEvents, var);
  Fill(data, dataEvents, var);
}

void HistSet::Fill(TH1D* hist, Events *events, int iLayer)
{
  if(events){
    for(int iEvent=0;iEvent<events->size();iEvent++){
      Event *event = events->At(iEvent);
      
      for(int iJet=0;iJet<event->GetNjets();iJet++){
        Jet *jet = event->GetJet(iJet);
        
        double value=0.0;
        
        if(var == kJetPt)   value = jet->GetPt();
        if(var == kJetEta)  value = jet->GetEta();
        if(var == kJetPhi)  value = jet->GetPhi();
        
        if(var == kJetPt || var == kJetEta || var == kJetPhi){
          hist->Fill(value);
        }
      }
      
      for(int iTrack=0;iTrack<events->At(iEvent)->GetNtracks();iTrack++){
        Track *track = events->At(iEvent)->GetTrack(iTrack);
        double value = 0.0;
        
        if(var == kTrackNclusters)            value = track->GetNclusters();
        else if(var == kTrackTotalDedx)       value = track->GetTotalDedx();
        else if(var == kTrackDedxPerCluster)  value = track->GetTotalDedx()/track->GetNclusters();
        else if(var == kTrackPt)              value = track->GetPt();
        else if(var == kTrackEta)             value = track->GetEta();
        else if(var == kTrackPhi)             value = track->GetPhi();
        else if(var == kTrackCaloEm)          value = track->GetCaloEmEnergy();
        else if(var == kTrackCaloHad)         value = track->GetCaloHadEnergy();
        else if(var == kTrackDxy)             value = track->GetDxy();
        else if(var == kTrackDz)              value = track->GetDz();
        else if(var == kTrackCharge)          value = track->GetCharge();
        else if(var == kTrackMass)            value = track->GetMass();
        else if(var == kTrackPid)             value = track->GetPid();
        
        else if(var == kDedx)   value = track->GetDeDxInLayer(iLayer);
        else if(var == kSizeX)  value = track->GetSizeXinLayer(iLayer);
        else if(var == kSizeY)  value = track->GetSizeYinLayer(iLayer);
        
        if(var == kDedx || var == kSizeX || var == kSizeY){
          if(value > 0.00001) hist->Fill(value);
        }
        else if(var == kTrackNclusters || var == kTrackTotalDedx || var == kTrackDedxPerCluster || var == kTrackPt
             || var == kTrackEta || var == kTrackPhi || var == kTrackCaloEm || var == kTrackCaloHad
             || var == kTrackDxy ||var == kTrackDz   || var == kTrackCharge || var == kTrackMass || var == kTrackPid){
          hist->Fill(value);
        }
      }
    }
  }
}

void HistSet::Draw(TCanvas *c1, int pad)
{
  TLegend *leg = GetLegend();
  leg->AddEntry(signal,"Signal","lp");
  leg->AddEntry(background,"Background","lp");
  leg->AddEntry(data,"Data","lp");
  
  c1->cd(pad);
  signal->SetLineColor(kBlue);
  signal->SetFillStyle(1000);
  signal->SetFillColorAlpha(kBlue, 0.2);
  if(ShouldNormalize()) signal->Scale(1/signal->Integral());
  if(!DoSumw2()) signal->Sumw2(false);
  
  background->SetLineColor(kRed);
  background->SetFillStyle(1000);
  background->SetFillColorAlpha(kRed, 0.2);
  if(ShouldNormalize())  background->Scale(1/background->Integral());
  if(!DoSumw2()) background->Sumw2(false);
  
  data->SetLineColor(kGreen);
  data->SetFillStyle(1000);
  data->SetFillColorAlpha(kGreen, 0.2);
  if(ShouldNormalize())  data->Scale(1/data->Integral());
  if(!DoSumw2()) data->Sumw2(false);
  
  THStack *stack = new THStack(GetTitle(),GetTitle());
  stack->Add(signal);
  stack->Add(background);
  if(analyzeData) stack->Add(data);
  
  stack->Draw("nostack");
  leg->Draw();
  c1->Update();
}

void HistSet::DrawPerLayer()
{
  TLegend *leg = GetLegend();
  leg->AddEntry(signalPerLayer[0],"Signal","lp");
  leg->AddEntry(backgroundPerLayer[0],"Background","lp");
  leg->AddEntry(dataPerLayer[0],"Data","lp");
  
  TCanvas *c1 = new TCanvas(GetTitle(),GetTitle(),2880,1800);
  
  if(var == kDedx) c1->Divide(4,4);
  if(var == kSizeX || var == kSizeY) c1->Divide(2,2);

  for(int iLayer=0;iLayer<nLayers;iLayer++){
    if(iLayer > 3 && (var==kSizeX || var==kSizeY)) continue;
    
    c1->cd(iLayer+1);
    
    signalPerLayer[iLayer]->SetLineColor(kBlue);
    signalPerLayer[iLayer]->SetFillStyle(1000);
    signalPerLayer[iLayer]->SetFillColorAlpha(kBlue,0.2);
    signalPerLayer[iLayer]->Scale(1/signalPerLayer[iLayer]->Integral());
    signalPerLayer[iLayer]->Sumw2(DoSumw2());
    
    backgroundPerLayer[iLayer]->SetLineColor(kRed+1);
    backgroundPerLayer[iLayer]->SetFillStyle(1000);
    backgroundPerLayer[iLayer]->SetFillColorAlpha(kRed, 0.2);
    backgroundPerLayer[iLayer]->Scale(1/backgroundPerLayer[iLayer]->Integral());
    backgroundPerLayer[iLayer]->Sumw2(DoSumw2());
    
    dataPerLayer[iLayer]->SetLineColor(kGreen+1);
    dataPerLayer[iLayer]->SetFillStyle(1000);
    dataPerLayer[iLayer]->SetFillColorAlpha(kGreen, 0.2);
    dataPerLayer[iLayer]->Scale(1/dataPerLayer[iLayer]->Integral());
    dataPerLayer[iLayer]->Sumw2(DoSumw2());

    THStack *stack = new THStack(Form("%s_layer[%i]",GetTitle(),iLayer),Form("%s_layer[%i]",GetTitle(),iLayer));
    stack->Add(signalPerLayer[iLayer]);
    stack->Add(backgroundPerLayer[iLayer]);
    if(analyzeData) stack->Add(dataPerLayer[iLayer]);
    
    stack->Draw("nostack");
    leg->Draw();
  }
  
  c1->Update();
}

TLegend* HistSet::GetLegend()
{
  double legendW=0.15, legendH=0.25, legendX=0.75, legendY=0.65;
  TLegend *leg = new TLegend(legendX,legendY,legendX+legendW,legendY+legendH);
  leg->SetHeader("Sample type:");
  return leg;
}


const char* HistSet::GetTitle()
{
  if(var == kCustom) return customTitle;
  
  if(var == kTrackNclusters) return "N clusters per track";
  if(var == kTrackTotalDedx) return "total dedx per track";
  if(var == kTrackDedxPerCluster) return "total dedx per track / n clusters";
  if(var == kTrackPt) return "pt dist";
  if(var == kTrackEta) return "eta dist";
  if(var == kTrackPhi) return "phi dist";
  if(var == kTrackCaloEm) return "EM calo energy";
  if(var == kTrackCaloHad) return "Hadron calo energy";
  if(var == kTrackDxy) return "Displacement in XY";
  if(var == kTrackDz) return "Displacement in Z";
  if(var == kTrackCharge) return "Charge dist";
  if(var == kTrackMass) return "Mass dist";
  if(var == kTrackPid) return "PDG PID";
  
  if(var == kJetPt) return "Jet pt";
  if(var == kJetEta) return "Jet eta";
  if(var == kJetPhi) return "Jet phi";
  
  if(var == kDedx) return "dedx";
  if(var == kSizeX) return "sizeX";
  if(var == kSizeY) return "sizeY";
  
  return "";
}

int HistSet::GetNbins()
{
  if(var == kTrackNclusters) return 20;
  if(var == kTrackTotalDedx) return 50;
  if(var == kTrackDedxPerCluster) return 50;
  if(var == kTrackPt) return 50;
  if(var == kTrackEta) return 50;
  if(var == kTrackPhi) return 50;
  if(var == kTrackCaloEm) return 100;
  if(var == kTrackCaloHad) return 100;
  if(var == kTrackDxy) return 100;
  if(var == kTrackDz) return 100;
  if(var == kTrackCharge) return 100;
  if(var == kTrackMass) return 500;
  if(var == kTrackPid) return 441;
  
  if(var == kJetPt) return 100;
  if(var == kJetEta) return 50;
  if(var == kJetPhi) return 50;
  
  if(var == kDedx) return 50;
  if(var == kSizeX) return 10;
  if(var == kSizeY) return 10;
  
  return 0;
}

double HistSet::GetMin()
{
  if(var == kTrackNclusters) return 0.0;
  if(var == kTrackTotalDedx) return 0.0;
  if(var == kTrackDedxPerCluster) return 0.0;
  if(var == kTrackPt) return 0.0;
  if(var == kTrackEta) return -3.0;
  if(var == kTrackPhi) return -3.5;
  if(var == kTrackCaloEm) return 0.0;
  if(var == kTrackCaloHad) return 0.0;
  if(var == kTrackDxy) return -0.02;
  if(var == kTrackDz) return -0.02;
  if(var == kTrackCharge) return -10.0;
  if(var == kTrackMass) return 0.0;
  if(var == kTrackPid) return -220;
  
  if(var == kJetPt) return 0.0;
  if(var == kJetEta) return -3.0;
  if(var == kJetPhi) return -3.5;
  
  if(var == kDedx) return 0.0;
  if(var == kSizeX) return 0.0;
  if(var == kSizeY) return 0.0;
  return 0.0;
}

double HistSet::GetMax()
{
  if(var == kTrackNclusters) return 20;
  if(var == kTrackTotalDedx) return 140;
  if(var == kTrackDedxPerCluster) return 14;
  if(var == kTrackPt) return 1000.0;
  if(var == kTrackEta) return 3.0;
  if(var == kTrackPhi) return 3.5;
  if(var == kTrackCaloEm) return 30.0;
  if(var == kTrackCaloHad) return 30.0;
  if(var == kTrackDxy) return 0.02;
  if(var == kTrackDz) return 0.02;
  if(var == kTrackCharge) return 10.0;
  if(var == kTrackMass) return 0.25;
  if(var == kTrackPid) return 220;
  
  if(var == kJetPt) return 1000.0;
  if(var == kJetEta) return 3.0;
  if(var == kJetPhi) return 3.5;
  
  if(var == kDedx) return 10.0;
  if(var == kSizeX) return 10.0;
  if(var == kSizeY) return 10.0;
  return 1000.0;
}

bool HistSet::ShouldNormalize()
{
  if(var == kCustom) return false;
  
  return true;
}

bool HistSet::DoSumw2()
{
  if(var == kCustom) return false;
  
  if(var == kTrackCaloEm)  return false;
  if(var == kTrackCaloHad) return false;
  
  if(var == kDedx) return false;
  if(var == kSizeX) return false;
  if(var == kSizeY) return false;
  
  return true;
}


