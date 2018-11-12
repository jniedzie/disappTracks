//
//  HistSet.cpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#include "HistSet.hpp"
#include "Helpers.hpp"
#include "Track.hpp"

#include <TLegend.h>
#include <THStack.h>
#include <TText.h>

HistSet::HistSet(EVar _var) :
var(_var)
{
  title= get<0>(settings.at(var)).c_str();
  nBins = get<1>(settings.at(var));
  min = get<2>(settings.at(var));
  max = get<3>(settings.at(var));
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    signal.push_back(new TH1D(Form("%s (signal %s)",title,signalTitle[iSig].c_str()),
                              Form("%s (signal %s)",title,signalTitle[iSig].c_str()),
                              nBins,min,max));
  }
  
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    background.push_back(new TH1D(Form("%s (background %s)",title,backgroundTitle[iBck].c_str()),
                                  Form("%s (background %s)",title,backgroundTitle[iBck].c_str()),
                                  nBins,min,max));
  }
  for(int iData=0;iData<kNdata;iData++){
    data.push_back(new TH1D(Form("%s (data %s)",title,dataTitle[iData].c_str()),
                            Form("%s (data %s)",title,dataTitle[iData].c_str()),
                            nBins,min,max));
  }
}

HistSet::~HistSet()
{
  
  
}

void HistSet::FillFromEvents(shared_ptr<EventSet> events)
{
  if(var == kDedx || var == kSizeX || var == kSizeY){
    FillFromEventsPerLayer(events);
  }
  else{
    for(int iSig=0;iSig<kNsignals;iSig++){
      if(!runSignal[iSig]) continue;
      Fill(signal[iSig], events,EventSet::kSignal, iSig);
    }
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      if(!runBackground[iBck]) continue;
      Fill(background[iBck],events,EventSet::kBackground, iBck);
    }
    for(int iData=0;iData<kNdata;iData++){
      if(!runData[iData]) continue;
      Fill(data[iData],events,EventSet::kData, iData);
    }
  }
}

void HistSet::FillFromEventsPerLayer(shared_ptr<EventSet> events)
{
  for(int iSig=0;iSig<kNsignals;iSig++){
    vector<TH1D*> signalVector;
    for(int iDetId=0;iDetId<nLayers;iDetId++){
      TH1D *histSignal = new TH1D(Form("%s_subDet[%i]_signal_%s",title,iDetId,signalTitle[iSig].c_str()),
                                  Form("%s_subDet[%i]_signal_%s",title,iDetId,signalTitle[iSig].c_str()),
                                  nBins,min,max);
      
      Fill(histSignal,events,EventSet::kSignal,iSig,iDetId);
      signalVector.push_back(histSignal);
    }
    signalPerLayer.push_back(signalVector);
  }
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    vector<TH1D*> backgroundVector;
    for(int iDetId=0;iDetId<nLayers;iDetId++){
      TH1D *histBackground = new TH1D(Form("%s_subDet[%i]_background_%s",title,iDetId,backgroundTitle[iBck].c_str()),
                                      Form("%s_subDet[%i]_background_%s",title,iDetId,backgroundTitle[iBck].c_str()),
                                      nBins,min,max);
      
      Fill(histBackground,events,EventSet::kBackground,iBck,iDetId);
      backgroundVector.push_back(histBackground);
    }
    backgroundPerLayer.push_back(backgroundVector);
  }
  for(int iData=0;iData<kNdata;iData++){
    vector<TH1D*> dataVector;
    for(int iDetId=0;iDetId<nLayers;iDetId++){
      
      TH1D *histData = new TH1D(Form("%s_subDet[%i]_data",title,iDetId),
                                Form("%s_subDet[%i]_data",title,iDetId),
                                nBins,min,max);
      Fill(histData,events,EventSet::kData,iData,iDetId);
      dataVector.push_back(histData);
    }
    dataPerLayer.push_back(dataVector);
  }
}

void HistSet::Fill(TH1D* hist, shared_ptr<EventSet> events,
                   EventSet::EDataType dataType, int setIter,
                   int iDetId)
{
  for(int iEvent=0;iEvent<events->size(dataType, setIter);iEvent++){
    shared_ptr<Event> event = events->At(dataType, setIter, iEvent);
    
    double value = 0.0;
    
    if(IsPerEventVariable(var)){
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
      
      hist->Fill(value, event->GetWeight());
    }
    
    if(IsPerJetVariable(var)){
      for(int iJet=0;iJet<event->GetNjets();iJet++){
        Jet *jet = event->GetJet(iJet);
        
        if(var == kJetPt)           value = jet->GetPt();
        else if(var == kJetEta)     value = jet->GetEta();
        else if(var == kJetPhi)     value = jet->GetPhi();
        else if(var == kMetJetDphi) value = event->GetMetPhi() - jet->GetPhi();
        else if(var == kJetCHF)     value = jet->GetChargedHadronEnergyFraction();
        else if(var == kJetNHF)     value = jet->GetNeutralHadronEnergyFraction();
        else if(var == kJetTrackDr){
          for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
            Track *track = event->GetTrack(iTrack);
            value = sqrt(pow(track->GetPhi() - jet->GetPhi(),2)+pow(track->GetEta() - jet->GetEta(),2));
            hist->Fill(value,event->GetWeight());
          }
          continue;
        }
        hist->Fill(value,event->GetWeight());
      }
    }
    
    if(IsPerTrackVariable(var)){
      for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
        Track *track = event->GetTrack(iTrack);
        
        if(var == kTrackNclusters)                    value = track->GetNdetIDs();
        else if(var == kTrackTotalDedx)               value = track->GetTotalDedx();
        else if(var == kTrackDedxPerCluster)          value = track->GetTotalDedx()/track->GetNdetIDs();
        else if(var == kTrackPt)                      value = track->GetPt();
        else if(var == kTrackEta)                     value = track->GetEta();
        else if(var == kTrackPhi)                     value = track->GetPhi();
        else if(var == kTrackCaloEm)                  value = track->GetCaloEmEnergy();
        else if(var == kTrackCaloHad)                 value = track->GetCaloHadEnergy();
        else if(var == kTrackDxy)                     value = track->GetDxy();
        else if(var == kTrackDz)                      value = track->GetDz();
        else if(var == kTrackCharge)                  value = track->GetCharge();
        else if(var == kTrackMass)                    value = track->GetMass();
        else if(var == kTrackPid)                     value = track->GetPid();
        else if(var == kTrackMissingOuterTrackerHits) value = track->GetNmissingOuterTrackerHits();
        else if(var == kTrackPixelHits)               value = track->GetNpixelHits();
        else if(var == kTrackTrackerHits)             value = track->GetNtrackerHits();
        else if(var == kTrackRelativeIsolation)       value = track->GetRelativeIsolation();
        else if(var == kTrackAbsoluteIsolation)       value = track->GetAbsoluteIsolation();
        else if(var == kTrackMetDphi)                 value = event->GetMetPhi() - track->GetPhi();
        else if(var == kDedx || var == kSizeX || var == kSizeY){
          for(int i=0;i<nLayers;i++){
            int detId = track->GetSubDetIdInLayer(i);
            if(detId == iDetId){
              if(var == kDedx)        value = track->GetDedxInSubDet(i);
              else if(var == kSizeX)  value = track->GetSizeXinLayer(i);
              else if(var == kSizeY)  value = track->GetSizeYinLayer(i);
              
              if(value > 0.00001) hist->Fill(value, event->GetWeight());
            }
          }
          continue;
        }
        else if(var == kTrackDedxPerHit){
          for(int i=0;i<nLayers;i++){
            value = track->GetDeDxInLayer(i);
            if(value > 0.00001) hist->Fill(value, event->GetWeight());
          }
          continue;
        }
        
        hist->Fill(value, event->GetWeight());
      }
    }
  }
}

void HistSet::Draw(TCanvas *c1, int pad)
{
  TLegend *leg = GetLegend();
  
  for(int iSig=0;iSig<(int)signal.size();iSig++){
    if(!runSignal[iSig]) continue;
    leg->AddEntry(signal[iSig],Form("Signal %s",signalTitle[iSig].c_str()),"lp");
  }
  for(int iBck=0;iBck<(int)background.size();iBck++){
    if(!runBackground[iBck]) continue;
    leg->AddEntry(background[iBck],Form("Background %s",backgroundTitle[iBck].c_str()),"lp");
  }
  for(int iData=0;iData<(int)data.size();iData++){
    if(!runData[iData]) continue;
    leg->AddEntry(data[iData],Form("Data  %s",dataTitle[iData].c_str()),"lp");
  }
  
  c1->cd(pad);
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    signal[iSig]->SetLineColor(SignalColor((ESignal)iSig));
    signal[iSig]->SetMarkerStyle(signalMarkers[iSig]);
    signal[iSig]->SetMarkerColor(SignalColor((ESignal)iSig));
    signal[iSig]->SetFillStyle(fillStyleSignal);
    signal[iSig]->SetFillColorAlpha(SignalColor((ESignal)iSig), fillOpacity);
    if(ShouldNormalize()) signal[iSig]->Scale(1/signal[iSig]->Integral());
  }
  double bckIntegral = 0;
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]) continue;
    bckIntegral += background[iBck]->Integral();
  }
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]) continue;
    background[iBck]->SetLineColor(BackColor((EBackground)iBck));
    background[iBck]->SetFillStyle(fillStyleBack);
    background[iBck]->SetFillColorAlpha(BackColor((EBackground)iBck), fillOpacity);
    if(ShouldNormalize()) background[iBck]->Scale(1/bckIntegral);
  }
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]) continue;
    data[iData]->SetLineColor(DataColor((EData)iData));
    data[iData]->SetMarkerColor(DataColor((EData)iData));
    data[iData]->SetMarkerStyle(20);
    data[iData]->SetMarkerSize(1.0);
    data[iData]->SetFillStyle(fillStyleData);
    data[iData]->SetFillColorAlpha(DataColor((EData)iData), fillOpacity);
    if(ShouldNormalize())  data[iData]->Scale(1/data[iData]->Integral());
  }
  
  THStack *backgroundStack = new THStack(title,title);
  THStack *signalStack = new THStack(title,title);
  THStack *dataStack = new THStack(title,title);
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]) continue;
    backgroundStack->Add(background[iBck]);
  }
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    signalStack->Add(signal[iSig]);
  }
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]) continue;
    dataStack->Add(data[iData]);
  }
  
  backgroundStack->Draw("HIST");
  signalStack->Draw("nostack,same,p");
  dataStack->Draw("nostack,same,p");
  
  if(showLegends) leg->Draw();
  c1->Update();
}

void HistSet::DrawPerLayer()
{
  TLegend *leg = GetLegend();
  for(int iSig=0;iSig<(int)signalPerLayer.size();iSig++){
    if(!runSignal[iSig]) continue;
    leg->AddEntry(signalPerLayer[iSig][1],Form("Signal %s",signalTitle[iSig].c_str()),"lp");
  }
  for(int iBck=0;iBck<(int)backgroundPerLayer.size();iBck++){
    if(!runBackground[iBck]) continue;
    leg->AddEntry(backgroundPerLayer[iBck][1],Form("Background %s",backgroundTitle[iBck].c_str()),"lp");
  }
  for(int iData=0;iData<(int)dataPerLayer.size();iData++){
    if(!runData[iData]) continue;
    leg->AddEntry(dataPerLayer[iData][1],Form("Data %s",dataTitle[iData].c_str()),"lp");
  }
  TCanvas *c1 = new TCanvas(title,title,2880,1800);
  
  if(var == kDedx) c1->Divide(2,3);
  if(var == kSizeX || var == kSizeY) c1->Divide(2,2);

  for(int iDetId=0;iDetId<nLayers;iDetId++){
    if(iDetId > 3 && (var==kSizeX || var==kSizeY)) continue;
    if((iDetId == 0 || iDetId > 6) && var==kDedx) continue;
    
    c1->cd(var==kDedx ? iDetId : iDetId+1);
    
    for(int iSig=0;iSig<(int)signalPerLayer.size();iSig++){
      signalPerLayer[iSig][iDetId]->SetLineColor(SignalColor((ESignal)iSig));
      signalPerLayer[iSig][iDetId]->SetMarkerStyle(signalMarkers[iSig]);
      signalPerLayer[iSig][iDetId]->SetMarkerColor(SignalColor((ESignal)iSig));
      signalPerLayer[iSig][iDetId]->SetFillStyle(fillStyleSignal);
      signalPerLayer[iSig][iDetId]->SetFillColorAlpha(SignalColor((ESignal)iSig),fillOpacity);
      signalPerLayer[iSig][iDetId]->Scale(1/signalPerLayer[iSig][iDetId]->Integral());
    }
      
    for(int iBck=0;iBck<(int)backgroundPerLayer.size();iBck++){
      backgroundPerLayer[iBck][iDetId]->SetLineColor(BackColor((EBackground)iBck));
      backgroundPerLayer[iBck][iDetId]->SetFillStyle(fillStyleBack);
      backgroundPerLayer[iBck][iDetId]->SetFillColorAlpha(BackColor((EBackground)iBck), fillOpacity);
      backgroundPerLayer[iBck][iDetId]->Scale(1/backgroundPerLayer[iBck][iDetId]->Integral());
    }
    for(int iData=0;iData<(int)dataPerLayer.size();iData++){
      dataPerLayer[iData][iDetId]->SetLineColor(DataColor((EData)iData));
      dataPerLayer[iData][iDetId]->SetMarkerColor(DataColor((EData)iData));
      dataPerLayer[iData][iDetId]->SetMarkerStyle(20);
      dataPerLayer[iData][iDetId]->SetMarkerSize(1.0);
      dataPerLayer[iData][iDetId]->SetFillStyle(fillStyleData);
      dataPerLayer[iData][iDetId]->SetFillColorAlpha(DataColor((EData)iData), fillOpacity);
      dataPerLayer[iData][iDetId]->Scale(1/dataPerLayer[iData][iDetId]->Integral());
    }
    
    THStack *backgroundStack  = new THStack(Form("%s_subDetId[%i]",title,iDetId),
                                            Form("%s_subDetId[%i]",title,iDetId));
    
    THStack *signalStack      = new THStack(Form("%s_subDetId[%i]",title,iDetId),
                                            Form("%s_subDetId[%i]",title,iDetId));
    
    THStack *dataStack        = new THStack(Form("%s_subDetId[%i]",title,iDetId),
                                            Form("%s_subDetId[%i]",title,iDetId));
    
    
    
    for(int iBck=0;iBck<(int)backgroundPerLayer.size();iBck++){backgroundStack->Add(backgroundPerLayer[iBck][iDetId]);}
    for(int iSig=0;iSig<(int)signalPerLayer.size();iSig++){  signalStack->Add(signalPerLayer[iSig][iDetId]);}
    for(int iData=0;iData<(int)dataPerLayer.size();iData++){dataStack->Add(dataPerLayer[iData][iDetId]);}
    
    backgroundStack->Draw("nostack,p");
    signalStack->Draw("nostack,same,p");
    dataStack->Draw("nostack,same,p");
    
    if(showLegends) leg->Draw();
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

bool HistSet::ShouldNormalize()
{
//  if(var == kCustom) return false;
//  if(var == kTrackTotalDedx) return false;
  return true;
}



