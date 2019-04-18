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
  title = get<0>(settings.at(var)).c_str();
  nBins = get<1>(settings.at(var));
  min   = get<2>(settings.at(var));
  max   = get<3>(settings.at(var));
  logy  = get<4>(settings.at(var));
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    string histTitle = title+" (signal "+signalTitle[iSig]+")";
    signal.push_back(shared_ptr<TH1D>(new TH1D(histTitle.c_str(), histTitle.c_str(), nBins,min,max)));
  }
  
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    string histTitle = title+" (background "+backgroundTitle[iBck]+")";
    background.push_back(shared_ptr<TH1D>(new TH1D(histTitle.c_str(), histTitle.c_str(), nBins,min,max)));
  }
  for(int iData=0;iData<kNdata;iData++){
    string histTitle = title+" (data "+dataTitle[iData]+")";
    data.push_back(shared_ptr<TH1D>(new TH1D(histTitle.c_str(), histTitle.c_str(), nBins,min,max)));
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
      if(!config->runSignal[iSig]) continue;
      Fill(signal[iSig], events,xtracks::kSignal, iSig);
    }
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      if(!config->runBackground[iBck]) continue;
      Fill(background[iBck],events,xtracks::kBackground, iBck);
    }
    for(int iData=0;iData<kNdata;iData++){
      if(!config->runData[iData]) continue;
      Fill(data[iData],events,xtracks::kData, iData);
    }
  }
}

void HistSet::FillFromEventsPerLayer(shared_ptr<EventSet> events)
{
  for(int iSig=0;iSig<kNsignals;iSig++){
    vector<shared_ptr<TH1D>> signalVector;
    for(int iDetId=0;iDetId<nLayers;iDetId++){
      string histTitle = title + "_subDet["+to_string(iDetId)+"]_signal_"+ signalTitle[iSig];
      auto histSignal = shared_ptr<TH1D>(new TH1D(histTitle.c_str(),histTitle.c_str(),nBins,min,max));
      
      Fill(histSignal,events,xtracks::kSignal,iSig,iDetId);
      signalVector.push_back(histSignal);
    }
    signalPerLayer.push_back(signalVector);
  }
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    vector<shared_ptr<TH1D>> backgroundVector;
    for(int iDetId=0;iDetId<nLayers;iDetId++){
      string histTitle = title + "_subDet["+to_string(iDetId)+"]_background_"+ backgroundTitle[iBck];
      auto histBackground = shared_ptr<TH1D>(new TH1D(histTitle.c_str(),histTitle.c_str(),nBins,min,max));
      
      Fill(histBackground,events,xtracks::kBackground,iBck,iDetId);
      backgroundVector.push_back(histBackground);
    }
    backgroundPerLayer.push_back(backgroundVector);
  }
  for(int iData=0;iData<kNdata;iData++){
    vector<shared_ptr<TH1D>> dataVector;
    for(int iDetId=0;iDetId<nLayers;iDetId++){
      string histTitle = title + "_subDet["+to_string(iDetId)+"]_data_"+ dataTitle[iData];
      auto histData = shared_ptr<TH1D>(new TH1D(histTitle.c_str(),histTitle.c_str(),nBins,min,max));
      
      Fill(histData,events,xtracks::kData,iData,iDetId);
      dataVector.push_back(histData);
    }
    dataPerLayer.push_back(dataVector);
  }
}

void HistSet::Fill(const shared_ptr<TH1D> &hist,
                   const shared_ptr<EventSet> &events,
                   xtracks::EDataType dataType, int setIter,
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
      else if(var == kNhelices)   value = event->GetNhelices();
      
      hist->Fill(value, event->GetWeight());
    }
    if(IsPerHelixVariable(var)){
      for(int iHelix=0;iHelix<event->GetNhelices();iHelix++){
        if(var == kHelixX)            value = event->GetHelix(iHelix)->GetOrigin().GetX();
        else if(var == kHelixY)       value = event->GetHelix(iHelix)->GetOrigin().GetY();
        else if(var == kHelixZ)       value = event->GetHelix(iHelix)->GetOrigin().GetZ();
        else if(var == kHelixPx)      value = event->GetHelix(iHelix)->GetMomentum()->GetX();
        else if(var == kHelixPy)      value = event->GetHelix(iHelix)->GetMomentum()->GetY();
        else if(var == kHelixPz)      value = event->GetHelix(iHelix)->GetMomentum()->GetZ();
        else if(var == kHelixCharge)  value = event->GetHelix(iHelix)->GetCharge();
        hist->Fill(value, event->GetWeight());
      }
    }
    if(IsPerJetVariable(var)){
      for(int iJet=0;iJet<event->GetNjets();iJet++){
        shared_ptr<Jet> jet = event->GetJet(iJet);
        
        if(var == kJetPt)           value = jet->GetPt();
        else if(var == kJetEta)     value = jet->GetEta();
        else if(var == kJetPhi)     value = jet->GetPhi();
        else if(var == kMetJetDphi){
          TLorentzVector metVector, jetVector;
          metVector.SetPtEtaPhiM(event->GetMetPt(), event->GetMetEta(), event->GetMetPhi(), event->GetMetMass());
          jetVector.SetPtEtaPhiM(jet->GetPt(), jet->GetEta(), jet->GetPhi(), jet->GetMass());
          value = metVector.DeltaPhi(jetVector);
        }
        else if(var == kJetCHF)     value = jet->GetChHEF();
        else if(var == kJetNHF)     value = jet->GetNeHEF();
        else if(var == kJetTrackDr){
          for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
            shared_ptr<Track> track = event->GetTrack(iTrack);
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
        shared_ptr<Track> track = event->GetTrack(iTrack);
        
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
        else if(var == kTrackTrackerLayers)           value = track->GetNtrackerLayers();
        else if(var == kTrackRelativeIsolation)       value = track->GetRelativeIsolation();
        else if(var == kTrackAbsoluteIsolation)       value = track->GetAbsoluteIsolation();
        else if(var == kTrackMetDphi){
          TLorentzVector metVector, trackVector;
          metVector.SetPtEtaPhiM(event->GetMetPt(), event->GetMetEta(), event->GetMetPhi(), event->GetMetMass());
          trackVector.SetPtEtaPhiM(track->GetPt(), track->GetEta(), track->GetPhi(), track->GetMass());
          value = metVector.DeltaPhi(trackVector);
        }
        else if(var == kDedx || var == kSizeX || var == kSizeY){
          for(int i=0;i<track->GetNdEdxHits();i++){
            int detId = track->GetSubDetIdForHit(i);
            if(detId == iDetId){
              if(var == kDedx)        value = track->GetDedxInSubDet(i);
              else if(var == kSizeX)  value = track->GetSizeXforHit(i);
              else if(var == kSizeY)  value = track->GetSizeYforHit(i);
              
              if(value > 0.00001) hist->Fill(value, event->GetWeight());
            }
          }
          continue;
        }
        else if(var == kTrackDedxPerHit){
          for(int i=0;i<track->GetNdEdxHits();i++){
            value = track->GetDeDxForHit(i);
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
    if(!config->runSignal[iSig]) continue;
    leg->AddEntry(&*signal[iSig],Form("Signal %s",signalTitle[iSig].c_str()),"lp");
  }
  for(int iBck=0;iBck<(int)background.size();iBck++){
    if(!config->runBackground[iBck]) continue;
    leg->AddEntry(&*background[iBck],Form("Background %s",backgroundTitle[iBck].c_str()),"lp");
  }
  for(int iData=0;iData<(int)data.size();iData++){
    if(!config->runData[iData]) continue;
    leg->AddEntry(&*data[iData],Form("Data  %s",dataTitle[iData].c_str()),"lp");
  }
  
  c1->cd(pad);
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!config->runSignal[iSig]) continue;
    signal[iSig]->SetLineColor(SignalColor((ESignal)iSig));
    signal[iSig]->SetMarkerStyle(signalMarkers[iSig]);
    signal[iSig]->SetMarkerColor(SignalColor((ESignal)iSig));
    signal[iSig]->SetFillStyle(fillStyleSignal);
    signal[iSig]->SetFillColorAlpha(SignalColor((ESignal)iSig), fillOpacity);
    if(ShouldNormalize()) signal[iSig]->Scale(1/signal[iSig]->Integral());
  }
  double bckIntegral = 0;
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!config->runBackground[iBck]) continue;
    bckIntegral += background[iBck]->Integral();
  }
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!config->runBackground[iBck]) continue;
    background[iBck]->SetLineColor(BackColor((EBackground)iBck));
    background[iBck]->SetFillStyle(fillStyleBack);
    background[iBck]->SetFillColorAlpha(BackColor((EBackground)iBck), fillOpacity);
    if(ShouldNormalize()) background[iBck]->Scale(1/bckIntegral);
  }
  for(int iData=0;iData<kNdata;iData++){
    if(!config->runData[iData]) continue;
    data[iData]->SetLineColor(DataColor((EData)iData));
    data[iData]->SetMarkerColor(DataColor((EData)iData));
    data[iData]->SetMarkerStyle(20);
    data[iData]->SetMarkerSize(1.0);
    data[iData]->SetFillStyle(fillStyleData);
    data[iData]->SetFillColorAlpha(DataColor((EData)iData), fillOpacity);
    if(ShouldNormalize())  data[iData]->Scale(1/data[iData]->Integral());
  }
  
  THStack *backgroundStack = new THStack(title.c_str(),title.c_str());
  THStack *signalStack     = new THStack(title.c_str(),title.c_str());
  THStack *dataStack       = new THStack(title.c_str(),title.c_str());
  
  int nActiveBackgrounds = 0;
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!config->runBackground[iBck]) continue;
    backgroundStack->Add(&*background[iBck]);
    nActiveBackgrounds++;
  }
  
  int nActiveSignals = 0;
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!config->runSignal[iSig]) continue;
    signalStack->Add(&*signal[iSig]);
    nActiveSignals++;
  }
  for(int iData=0;iData<kNdata;iData++){
    if(!config->runData[iData]) continue;
    dataStack->Add(&*data[iData]);
  }
  
  if(nActiveBackgrounds > 0){
    backgroundStack->Draw("HIST");
    signalStack->Draw("nostack,same,p");
    dataStack->Draw("nostack,same,p");
    
    double maxValue = backgroundStack->GetMaximum();
    if(signalStack->GetMaximum() > max) maxValue = signalStack->GetMaximum();
    if(dataStack->GetMaximum() > max)   maxValue = dataStack->GetMaximum();
    backgroundStack->GetYaxis()->SetLimits(0, 1.1*maxValue);
  }
  else if(nActiveSignals > 0){
    signalStack->Draw("nostack,p");
    dataStack->Draw("nostack,same,p");
    
    double maxValue = signalStack->GetMaximum();
    if(dataStack->GetMaximum() > max)   maxValue = dataStack->GetMaximum();
    signalStack->GetYaxis()->SetLimits(0, 1.1*maxValue);
  }
  else{
    dataStack->Draw("nostack,p");
    dataStack->GetYaxis()->SetLimits(0, 1.1*dataStack->GetMaximum());
  }
  
  if(logy) gPad->SetLogy();
  
  if(config->showLegends) leg->Draw();
  c1->Update();
}

void HistSet::DrawPerLayer()
{
  TLegend *leg = GetLegend();
  for(int iSig=0;iSig<(int)signalPerLayer.size();iSig++){
    if(!config->runSignal[iSig]) continue;
    leg->AddEntry(&*signalPerLayer[iSig][1],Form("Signal %s",signalTitle[iSig].c_str()),"lp");
  }
  for(int iBck=0;iBck<(int)backgroundPerLayer.size();iBck++){
    if(!config->runBackground[iBck]) continue;
    leg->AddEntry(&*backgroundPerLayer[iBck][1],Form("Background %s",backgroundTitle[iBck].c_str()),"lp");
  }
  for(int iData=0;iData<(int)dataPerLayer.size();iData++){
    if(!config->runData[iData]) continue;
    leg->AddEntry(&*dataPerLayer[iData][1],Form("Data %s",dataTitle[iData].c_str()),"lp");
  }
  TCanvas *c1 = new TCanvas(title.c_str(),title.c_str(),2880,1800);
  
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
    
    THStack *backgroundStack  = new THStack(Form("%s_subDetId[%i]",title.c_str(),iDetId),
                                            Form("%s_subDetId[%i]",title.c_str(),iDetId));
    
    THStack *signalStack      = new THStack(Form("%s_subDetId[%i]",title.c_str(),iDetId),
                                            Form("%s_subDetId[%i]",title.c_str(),iDetId));
    
    THStack *dataStack        = new THStack(Form("%s_subDetId[%i]",title.c_str(),iDetId),
                                            Form("%s_subDetId[%i]",title.c_str(),iDetId));
    
    
    
    for(int iBck=0;iBck<(int)backgroundPerLayer.size();iBck++){backgroundStack->Add(&*backgroundPerLayer[iBck][iDetId]);}
    for(int iSig=0;iSig<(int)signalPerLayer.size();iSig++){  signalStack->Add(&*signalPerLayer[iSig][iDetId]);}
    for(int iData=0;iData<(int)dataPerLayer.size();iData++){dataStack->Add(&*dataPerLayer[iData][iDetId]);}
    
    backgroundStack->Draw("nostack,p");
    signalStack->Draw("nostack,same,p");
    dataStack->Draw("nostack,same,p");
    
    if(config->showLegends) leg->Draw();
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



