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

HistSet::HistSet(EVar var) :
signal(nullptr),
background(nullptr),
data(nullptr)
{
  const char* title = GetTitle(var);
  int nBins = GetNbins(var);
  double min = GetMin(var);
  double max = GetMax(var);
  
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
  int nBins = GetNbins(var);
  double min = GetMin(var);
  double max = GetMax(var);
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    TH1D *histSignal = new TH1D(Form("%s_layer[%i]_signal",title,iLayer),Form("%s_layer[%i]_signal",title, iLayer),nBins,min,max);
    
    Fill(histSignal,signalEvents,var,iLayer);
    signalPerLayer.push_back(histSignal);
    
    TH1D *histBackground = new TH1D(Form("%s_layer[%i]_background",title,iLayer),Form("%s_layer[%i]_background",title,iLayer),nBins,min,max);
    Fill(histBackground,backgroundEvents,var,iLayer);
    backgroundPerLayer.push_back(histBackground);
    
    TH1D *histData = new TH1D(Form("%s_layer[%i]_data",title,iLayer),Form("%s_layer[%i]_data",title,iLayer),nBins,min,max);
    Fill(histData,dataEvents,var,iLayer);
    dataPerLayer.push_back(histData);
  }
}

void HistSet::FillFromEventsGlobal(Events *signalEvents, Events *backgroundEvents, Events *dataEvents, EVar var)
{
  Fill(signal, signalEvents, var);
  Fill(background, backgroundEvents, var);
  Fill(data, dataEvents, var);
}

void HistSet::Fill(TH1D* hist, Events *events, EVar var, int iLayer)
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
        
        if(var == kJetPt || var == kJetEta || var == kJetPt){
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

int HistSet::GetNbins(EVar var)
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

double HistSet::GetMin(EVar var)
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

double HistSet::GetMax(EVar var)
{
  if(var == kTrackNclusters) return 20;
  if(var == kTrackTotalDedx) return 200;
  if(var == kTrackDedxPerCluster) return 20;
  if(var == kTrackPt) return 1000.0;
  if(var == kTrackEta) return 3.0;
  if(var == kTrackPhi) return 3.5;
  if(var == kTrackCaloEm) return 500.0;
  if(var == kTrackCaloHad) return 500.0;
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


