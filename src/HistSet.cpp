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
#include <TText.h>

HistSet::HistSet(const char* title, int nBins, double min, double max) :
data(nullptr),
var(kCustom),
customTitle(title),
showNonZeroBinPosX(false)
{
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    signal[iSig] = new TH1D(Form("%s (signal %s)",title,signalTitle[iSig]),
                            Form("%s (signal %s)",title,signalTitle[iSig]),
                            nBins,min,max);
    signal[iSig]->Sumw2(DoSumw2());
  }
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    background[iBck] = new TH1D(Form("%s (background %s)",title,backgroundTitle[iBck]),
                                Form("%s (background %s)",title,backgroundTitle[iBck]),
                                nBins,min,max);
    
    background[iBck]->Sumw2(DoSumw2());
  }
  
  
  data = new TH1D(Form("%s (data)",title),Form("%s (data)",title),nBins,min,max);
  data->Sumw2(DoSumw2());
}

HistSet::HistSet(EVar _var) :
data(nullptr),
var(_var),
showNonZeroBinPosX(false)
{
  const char* title = GetTitle();
  int nBins = GetNbins();
  double min = GetMin();
  double max = GetMax();
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    signal[iSig] = new TH1D(Form("%s (signal %s)",title,signalTitle[iSig]),
                            Form("%s (signal %s)",title,signalTitle[iSig]),
                            nBins,min,max);
    signal[iSig]->Sumw2(DoSumw2());
  }
  
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    background[iBck] = new TH1D(Form("%s (background %s)",title,backgroundTitle[iBck]),
                                Form("%s (background %s)",title,backgroundTitle[iBck]),
                                nBins,min,max);
    background[iBck]->Sumw2(DoSumw2());
  }
  
  data       = new TH1D(Form("%s (data)",title),Form("%s (data)",title),nBins,min,max);
  data->Sumw2(DoSumw2());
}

HistSet::~HistSet()
{
  
  
}

void HistSet::FillFromEvents(Events *signalEvents[kNsignals], Events *backgroundEvents[kNbackgrounds], Events *dataEvents)
{
  if(var == kDedx || var == kSizeX || var == kSizeY){
    FillFromEventsPerLayer(signalEvents, backgroundEvents, dataEvents);
  }
  else{
    for(int iSig=0;iSig<kNsignals;iSig++){
      Fill(signal[iSig], signalEvents[iSig]);
    }
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      Fill(background[iBck], backgroundEvents[iBck]);
    }
    Fill(data, dataEvents);
  }
}

void HistSet::FillFromEventsPerLayer(Events *signalEvents[kNsignals], Events *backgroundEvents[kNbackgrounds], Events *dataEvents)
{
  const char* title= GetTitle();
  int nBins = GetNbins();
  double min = GetMin();
  double max = GetMax();
  
  for(int iLayer=0;iLayer<nLayers;iLayer++){
    
    for(int iSig=0;iSig<kNsignals;iSig++){
      TH1D *histSignal = new TH1D(Form("%s_layer[%i]_signal_%s",title,iLayer,signalTitle[iSig]),
                                  Form("%s_layer[%i]_signal_%s",title,iLayer,signalTitle[iSig]),
                                  nBins,min,max);
      
      Fill(histSignal,signalEvents[iSig],iLayer);
      signalPerLayer[iSig].push_back(histSignal);
    }
      
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      TH1D *histBackground = new TH1D(Form("%s_layer[%i]_background_%s",title,iLayer,backgroundTitle[iBck]),
                                      Form("%s_layer[%i]_background_%s",title,iLayer,backgroundTitle[iBck]),
                                      nBins,min,max);
      
      Fill(histBackground,backgroundEvents[iBck],iLayer);
      backgroundPerLayer[iBck].push_back(histBackground);
    }
    TH1D *histData = new TH1D(Form("%s_layer[%i]_data",title,iLayer),Form("%s_layer[%i]_data",title,iLayer),nBins,min,max);
    Fill(histData,dataEvents,iLayer);
    dataPerLayer.push_back(histData);
  }
}

void HistSet::Fill(TH1D* hist, Events *events, int iLayer)
{
  if(events){
    for(int iEvent=0;iEvent<events->size();iEvent++){
      Event *event = events->At(iEvent);
      
      double value = 0.0;
      
      if(var == kNvertices)       value = event->GetNvertices();
      else if(var == kNisoTracks) value = event->GetNtracks();
      else if(var == kNjets)      value = event->GetNjets();
      else if(var == kNjets30)    value = event->GetNjet30();
      else if(var == kNjets30a)   value = event->GetNjet30a();
      else if(var == kMetSumEt)   value = event->GetMetSumEt();
      else if(var == kMetPt)      value = event->GetMetPt();
      else if(var == kMetMass)    value = event->GetMetMass();
      else if(var == kMetEta)     value = event->GetMetEta();
      else if(var == kMetPhi)     value = event->GetMetPhi();
      
//      if(var == kNjets) cout<<value<<endl;
      
      if(var == kNvertices || var == kNisoTracks || var == kNjets || var == kNjets30 || var == kNjets30a || var == kMetSumEt || var == kMetPt || var == kMetMass || var == kMetEta || var == kMetPhi){
        hist->Fill(value, event->GetWeight());
        continue;
      }
      
      for(int iJet=0;iJet<event->GetNjets();iJet++){
        Jet *jet = event->GetJet(iJet);
        
        double value=0.0;
        
        if(var == kJetPt)   value = jet->GetPt();
        else if(var == kJetEta)  value = jet->GetEta();
        else if(var == kJetPhi)  value = jet->GetPhi();
        else if(var == kMetJetDphi) value = event->GetMetPhi() - jet->GetPhi();
        
        if(var == kJetPt || var == kJetEta || var == kJetPhi || var == kMetJetDphi){
          hist->Fill(value,event->GetWeight());
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
          if(value > 0.00001) hist->Fill(value, event->GetWeight());
        }
        else if(var == kTrackNclusters || var == kTrackTotalDedx || var == kTrackDedxPerCluster || var == kTrackPt
             || var == kTrackEta || var == kTrackPhi || var == kTrackCaloEm || var == kTrackCaloHad
             || var == kTrackDxy ||var == kTrackDz   || var == kTrackCharge || var == kTrackMass || var == kTrackPid){
          hist->Fill(value, event->GetWeight());
        }
      }
    }
  }
}

void HistSet::Draw(TCanvas *c1, int pad)
{
  TLegend *leg = GetLegend();
  if(showNonZeroBinPosX){
    for(int iSig=0;iSig<kNsignals;iSig++){
      leg->AddEntry(signal[iSig],
                    Form("Signal %s (%.2f %%)",signalTitle[iSig], 100*GetNonZeroBinPosX(signal[iSig])),"lp");
    }
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      leg->AddEntry(background[iBck],
                    Form("Background %s (%.2f %%)",backgroundTitle[iBck], 100*GetNonZeroBinPosX(background[iBck])),
                    "lp");
    }
    leg->AddEntry(data,Form("Data (%.2f %%)",100*GetNonZeroBinPosX(data)),"lp");
  }
  else{
    for(int iSig=0;iSig<kNsignals;iSig++){
      leg->AddEntry(signal[iSig],Form("Signal %s",signalTitle[iSig]),"lp");
    }
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      leg->AddEntry(background[iBck],Form("Background %s",backgroundTitle[iBck]),"lp");
    }
    leg->AddEntry(data,"Data","lp");
  }
    
  c1->cd(pad);
  for(int iSig=0;iSig<kNsignals;iSig++){
    signal[iSig]->SetLineColor(SignalColor((ESignal)iSig));
//    signal[iSig]->SetLineStyle(3);
    signal[iSig]->SetMarkerStyle(signalMarkers[iSig]);
    signal[iSig]->SetMarkerColor(SignalColor((ESignal)iSig));
    signal[iSig]->SetFillStyle(fillStyleSignal);
    signal[iSig]->SetFillColorAlpha(SignalColor((ESignal)iSig), fillOpacity);
    if(ShouldNormalize()) signal[iSig]->Scale(1/signal[iSig]->Integral());
    if(!DoSumw2()) signal[iSig]->Sumw2(false);
  }
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    background[iBck]->SetLineColor(BackColor((EBackground)iBck));
    background[iBck]->SetFillStyle(fillStyleBack);
    background[iBck]->SetFillColorAlpha(BackColor((EBackground)iBck), fillOpacity);
    if(ShouldNormalize())  background[iBck]->Scale(1/background[iBck]->Integral());
    if(!DoSumw2()) background[iBck]->Sumw2(false);
  }
  data->SetLineColor(kGreen);
  data->SetFillStyle(1000);
  data->SetFillColorAlpha(kGreen, fillOpacity);
  if(ShouldNormalize())  data->Scale(1/data->Integral());
  if(!DoSumw2()) data->Sumw2(false);
  
  THStack *stack = new THStack(GetTitle(),GetTitle());
  for(int iSig=0;iSig<kNsignals;iSig++){
    stack->Add(signal[iSig]);
  }
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    stack->Add(background[iBck]);
  }
  if(analyzeData) stack->Add(data);
  stack->Draw("nostack");
  
  leg->Draw();
  c1->Update();
}

void HistSet::DrawPerLayer()
{
  TLegend *leg = GetLegend();
  for(int iSig=0;iSig<kNsignals;iSig++){
    leg->AddEntry(signalPerLayer[iSig][0],Form("Signal %s",signalTitle[iSig]),"lp");
  }
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    leg->AddEntry(backgroundPerLayer[iBck][0],Form("Background %s",backgroundTitle[iBck]),"lp");
  }
  leg->AddEntry(dataPerLayer[0],"Data","lp");
  
  TCanvas *c1 = new TCanvas(GetTitle(),GetTitle(),2880,1800);
  
  if(var == kDedx) c1->Divide(4,4);
  if(var == kSizeX || var == kSizeY) c1->Divide(2,2);

  for(int iLayer=0;iLayer<nLayers;iLayer++){
    if(iLayer > 3 && (var==kSizeX || var==kSizeY)) continue;
    
    c1->cd(iLayer+1);
    
    for(int iSig=0;iSig<kNsignals;iSig++){
      signalPerLayer[iSig][iLayer]->SetLineColor(SignalColor((ESignal)iSig));
      signalPerLayer[iSig][iLayer]->SetFillStyle(fillStyleSignal);
      signalPerLayer[iSig][iLayer]->SetFillColorAlpha(SignalColor((ESignal)iSig),fillOpacity);
      signalPerLayer[iSig][iLayer]->Scale(1/signalPerLayer[iSig][iLayer]->Integral());
      signalPerLayer[iSig][iLayer]->Sumw2(DoSumw2());
    }
      
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      backgroundPerLayer[iBck][iLayer]->SetLineColor(BackColor((EBackground)iBck));
      backgroundPerLayer[iBck][iLayer]->SetFillStyle(fillStyleBack);
      backgroundPerLayer[iBck][iLayer]->SetFillColorAlpha(BackColor((EBackground)iBck), fillOpacity);
      backgroundPerLayer[iBck][iLayer]->Scale(1/backgroundPerLayer[iBck][iLayer]->Integral());
      backgroundPerLayer[iBck][iLayer]->Sumw2(DoSumw2());
    }
    dataPerLayer[iLayer]->SetLineColor(kGreen+1);
    dataPerLayer[iLayer]->SetFillStyle(1000);
    dataPerLayer[iLayer]->SetFillColorAlpha(kGreen, fillOpacity);
    dataPerLayer[iLayer]->Scale(1/dataPerLayer[iLayer]->Integral());
    dataPerLayer[iLayer]->Sumw2(DoSumw2());

    THStack *stack = new THStack(Form("%s_layer[%i]",GetTitle(),iLayer),Form("%s_layer[%i]",GetTitle(),iLayer));
    for(int iSig=0;iSig<kNsignals;iSig++){
      stack->Add(signalPerLayer[iSig][iLayer]);
    }
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      stack->Add(backgroundPerLayer[iBck][iLayer]);
    }
    if(analyzeData) stack->Add(dataPerLayer[iLayer]);
    
    stack->Draw("nostack");
    leg->Draw();
  }
  
  c1->Update();
}

TLegend* HistSet::GetLegend()
{
  double legendW=0.25, legendH=0.80, legendX=0.65, legendY=0.1;
  TLegend *leg = new TLegend(legendX,legendY,legendX+legendW,legendY+legendH);
  leg->SetHeader("Sample type:");
  return leg;
}


const char* HistSet::GetTitle()
{
  if(var == kCustom) return customTitle;
  
  if(var == kNvertices) return "N good vertices";
  if(var == kNisoTracks) return "N iso tracks";
  if(var == kNjets) return "N jets";
  if(var == kNjets30) return "N jets with pt > 30, |eta|<2.4";
  if(var == kNjets30a) return "N jets with pt > 30, |eta|<4.7";
  if(var == kMetSumEt) return "MET sum Et";
  if(var == kMetPt) return "MET pT";
  if(var == kMetMass) return "MET mass";
  if(var == kMetEta) return "MET eta";
  if(var == kMetPhi) return "MET phi";
  if(var == kMetJetDphi) return "#Delta #phi (p_{T}^{jet},p_{T}^{MET})";
  
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
  if(var == kNvertices)   return 50;
  if(var == kNisoTracks)  return 20;
  if(var == kNjets)       return 20;
  if(var == kNjets30)     return 20;
  if(var == kNjets30a)    return 20;
  if(var == kMetSumEt)    return 100;
  if(var == kMetPt)       return 100;
  if(var == kMetMass)     return 100;
  if(var == kMetEta)      return 100;
  if(var == kMetPhi)      return 100;
  if(var == kMetJetDphi)  return 100;
  
  if(var == kTrackNclusters)      return 20;
  if(var == kTrackTotalDedx)      return 50;
  if(var == kTrackDedxPerCluster) return 50;
  if(var == kTrackPt)             return 50;
  if(var == kTrackEta)            return 50;
  if(var == kTrackPhi)            return 50;
  if(var == kTrackCaloEm)         return 100;
  if(var == kTrackCaloHad)        return 100;
  if(var == kTrackDxy)            return 100;
  if(var == kTrackDz)             return 100;
  if(var == kTrackCharge)         return 100;
  if(var == kTrackMass)           return 500;
  if(var == kTrackPid)            return 441;
  
  if(var == kJetPt)   return 100;
  if(var == kJetEta)  return 50;
  if(var == kJetPhi)  return 50;
  
  if(var == kDedx)    return 50;
  if(var == kSizeX)   return 10;
  if(var == kSizeY)   return 10;
  
  return 0;
}

double HistSet::GetMin()
{
  if(var == kNvertices)   return 0;
  if(var == kNisoTracks)  return 0;
  if(var == kNjets)       return 0;
  if(var == kNjets30)     return 0;
  if(var == kNjets30a)    return 0;
  if(var == kMetSumEt)    return -20;
  if(var == kMetPt)       return 0;
  if(var == kMetMass)     return -10e-6;
  if(var == kMetEta)      return -3.5;
  if(var == kMetPhi)      return -3.5;
  if(var == kMetJetDphi)  return -3.5;
  
  if(var == kTrackNclusters)      return 0.0;
  if(var == kTrackTotalDedx)      return 0.0;
  if(var == kTrackDedxPerCluster) return 0.0;
  if(var == kTrackPt)             return 0.0;
  if(var == kTrackEta)            return -3.0;
  if(var == kTrackPhi)            return -3.5;
  if(var == kTrackCaloEm)         return 0.0;
  if(var == kTrackCaloHad)        return 0.0;
  if(var == kTrackDxy)            return -0.02;
  if(var == kTrackDz)             return -0.02;
  if(var == kTrackCharge)         return -10.0;
  if(var == kTrackMass)           return 0.0;
  if(var == kTrackPid)            return -220;
  
  if(var == kJetPt)   return 0.0;
  if(var == kJetEta)  return -3.0;
  if(var == kJetPhi)  return -3.5;
  
  if(var == kDedx)    return 0.0;
  if(var == kSizeX)   return 0.0;
  if(var == kSizeY)   return 0.0;
  return 0.0;
}

double HistSet::GetMax()
{
  if(var == kNvertices)   return 100;
  if(var == kNisoTracks)  return 20;
  if(var == kNjets)       return 20;
  if(var == kNjets30)     return 20;
  if(var == kNjets30a)    return 20;
  if(var == kMetSumEt)    return 5000;
  if(var == kMetPt)       return 500;
  if(var == kMetMass)     return 10e-6;
  if(var == kMetEta)      return 3.5;
  if(var == kMetPhi)      return 3.5;
  if(var == kMetJetDphi)  return 3.5;
  
  if(var == kTrackNclusters)      return 22;
  if(var == kTrackTotalDedx)      return 140;
  if(var == kTrackDedxPerCluster) return 14;
  if(var == kTrackPt)             return 1000.0;
  if(var == kTrackEta)            return 3.0;
  if(var == kTrackPhi)            return 3.5;
  if(var == kTrackCaloEm)         return 30.0;
  if(var == kTrackCaloHad)        return 30.0;
  if(var == kTrackDxy) 	          return 0.02;
  if(var == kTrackDz)             return 0.02;
  if(var == kTrackCharge)         return 10.0;
  if(var == kTrackMass)           return 0.25;
  if(var == kTrackPid)            return 220;
  
  if(var == kJetPt)   return 1000.0;
  if(var == kJetEta)  return 3.0;
  if(var == kJetPhi)  return 3.5;
  
  if(var == kDedx)    return 13.0;
  if(var == kSizeX)   return 13.0;
  if(var == kSizeY)   return 13.0;
  return 1000.0;
}

bool HistSet::ShouldNormalize()
{
//  if(var == kCustom) return false;
  
//  if(var == kTrackTotalDedx) return false;
  
  return false;
}

bool HistSet::DoSumw2()
{
  if(var == kCustom) return false;
  
//  if(var == kNvertices)   return false;
//  if(var == kNisoTracks)  return false;
//  if(var == kNjets)       return false;
//  if(var == kNjets30)     return false;
//  if(var == kNjets30a)    return false;
//  if(var == kMetSumEt)    return false;
//  if(var == kMetPt)       return false;
//  if(var == kMetMass)     return false;
//  if(var == kMetEta)      return false;
//  if(var == kMetPhi)      return false;
  
//  if(var == kTrackNclusters)  return false;
//  if(var == kTrackPt)         return false;
//  if(var == kTrackCaloEm)     return false;
//  if(var == kTrackCaloHad)    return false;
  
//  if(var == kDedx)  return false;
//  if(var == kSizeX) return false;
//  if(var == kSizeY) return false;
  
  return true;
}

double HistSet::GetNonZeroBinPosX(TH1D *hist)
{
  for(int i=0;i<hist->GetNbinsX();i++){
    if(hist->GetBinContent(i) != 0){
      return hist->GetBinCenter(i);
    }
  }
  return -1;
}


