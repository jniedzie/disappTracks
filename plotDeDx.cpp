#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"

#include <TApplication.h>

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  //---------------------------------------------------------------------------
  // Load MC and data files
  //---------------------------------------------------------------------------
  Events *eventsSignal[kNsignals];
  Events *eventsBackground[kNbackgrounds];
  Events *eventsData = analyzeData ? new Events(inFileNameData,2) : nullptr;
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    eventsSignal[iSig] = new Events(inFileNameSignal[iSig],1);
  }
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
   eventsBackground[iBck] = new Events(inFileNameBackground[iBck],0);
  }
  
  
  //---------------------------------------------------------------------------
  // Define event, track and jet cuts
  //---------------------------------------------------------------------------
  EventCut *oneTrackOneJetEventCut = new EventCut(EventCut::kOneTrack);
  EventCut *highMetEventCut = new EventCut(EventCut::kMet100GeVOneTrackOneJet);
  
  TrackCut *shortTrackCut = new TrackCut(TrackCut::kShort);
  TrackCut *shortAboveTrasholdTrackCut = new TrackCut(TrackCut::kShortAboveThreshold);
  TrackCut *shortLowDedxTrackCut = new TrackCut(TrackCut::kShortLowTotal);
  
  JetCut *highPtJetCut = new JetCut(JetCut::kHighPt);
  
  EventCut  *bestEventCut = new EventCut(EventCut::kOneTrack);
  TrackCut  *bestTrackCut = new TrackCut(TrackCut::kEmpty);
  JetCut    *bestJetCut   = new JetCut(JetCut::kEmpty);
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    eventsSignal[iSig] = eventsSignal[iSig]->ApplyCuts(bestEventCut, bestTrackCut, bestJetCut);
  }

//  for(int iBck=0;iBck<kNbackgrounds;iBck++){
//    eventsBackground[iBck] = eventsBackground[iBck]->ApplyCuts(bestEventCut, bestTrackCut, bestJetCut);
//  }
  
  //---------------------------------------------------------------------------
  // Create standard per event, per track and per jet plots
  //---------------------------------------------------------------------------
  HistSet *nVertices = new HistSet(kNvertices);
  HistSet *nIsoTrack = new HistSet(kNisoTracks);
  HistSet *nJet      = new HistSet(kNjets);
  HistSet *nJet30    = new HistSet(kNjets30);
  HistSet *nJet30a   = new HistSet(kNjets30a);
  HistSet *nMetSumEt = new HistSet(kMetSumEt);
  HistSet *nMetPt    = new HistSet(kMetPt);
  HistSet *nMetMass  = new HistSet(kMetMass);
  HistSet *nMetEta   = new HistSet(kMetEta);
  HistSet *nMetPhi   = new HistSet(kMetPhi);
  
  HistSet *nClustersPerTrack = new HistSet(kTrackNclusters);
  HistSet *totalDeDx = new HistSet(kTrackTotalDedx);
  
  HistSet *totalDeDxByNclusters = new HistSet(kTrackDedxPerCluster);
  HistSet *pt = new HistSet(kTrackPt);
  HistSet *eta = new HistSet(kTrackEta);
  HistSet *phi = new HistSet(kTrackPhi);
  HistSet *caloEm = new HistSet(kTrackCaloEm);
  HistSet *caloHad = new HistSet(kTrackCaloHad);
  
  HistSet *dxy = new HistSet(kTrackDxy);
  HistSet *dz = new HistSet(kTrackDz);
  HistSet *charge = new HistSet(kTrackCharge);
  HistSet *mass = new HistSet(kTrackMass);
  HistSet *pid = new HistSet(kTrackPid);
  
  HistSet *jet_pt = new HistSet(kJetPt);
  HistSet *jet_eta = new HistSet(kJetEta);
  HistSet *jet_phi = new HistSet(kJetPhi);
  
  nVertices->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  nIsoTrack->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  nJet->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  nJet30->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  nJet30a->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  nMetSumEt->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  nMetPt->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  nMetMass->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  nMetEta->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  nMetPhi->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  
  nClustersPerTrack->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  totalDeDx->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  
  totalDeDxByNclusters->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  pt->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  eta->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  phi->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  caloEm->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  caloHad->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  
  dxy->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  dz->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  charge->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  mass->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  pid->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  
  jet_pt->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  jet_eta->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  jet_phi->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  
  
  // Plot histograms
  TCanvas *canvasEvents = new TCanvas("Events","Events",2880,1800);
  canvasEvents->Divide(3,4);
  
  nVertices->Draw(canvasEvents,1);
  nIsoTrack->Draw(canvasEvents,2);
  nJet->Draw(canvasEvents,3);
  nJet30->Draw(canvasEvents,4);
  nJet30a->Draw(canvasEvents,5);
  nMetSumEt->Draw(canvasEvents,6);
  nMetPt->Draw(canvasEvents,7);
  nMetMass->Draw(canvasEvents,8);
  nMetEta->Draw(canvasEvents,9);
  nMetPhi->Draw(canvasEvents,10);
  
  TCanvas *canvasTrack = new TCanvas("Tracks","Tracks",2880,1800);
  canvasTrack->Divide(3,3);
  
  nClustersPerTrack->Draw(canvasTrack,1);
  totalDeDx->Draw(canvasTrack,2);
  totalDeDxByNclusters->Draw(canvasTrack,3);
  pt->Draw(canvasTrack,4);
  eta->Draw(canvasTrack,5);
  phi->Draw(canvasTrack,6);
  caloEm->Draw(canvasTrack,7);
  caloHad->Draw(canvasTrack,8);
//  charge->Draw(canvasTrack,3);
//  mass->Draw(canvasTrack,4);
//  pid->Draw(canvasTrack,5);
//  dxy->Draw(canvasTrack,1);
//  dz->Draw(canvasTrack,2);
  
  TCanvas *canvasJets = new TCanvas("Jets","Jets",2880,1800);
  canvasJets->Divide(2,2);
  jet_pt->Draw(canvasJets, 1);
  jet_eta->Draw(canvasJets, 2);
  jet_phi->Draw(canvasJets, 3);
  
  //---------------------------------------------------------------------------
  // Create per layer plots
  //---------------------------------------------------------------------------
  
  HistSet *dedxPerLayer = new HistSet(kDedx);
  dedxPerLayer->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  dedxPerLayer->DrawPerLayer();
  
  HistSet *sizeXperLayer = new HistSet(kSizeX);
  sizeXperLayer->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  sizeXperLayer->DrawPerLayer();
  
  HistSet *sizeYperLayer = new HistSet(kSizeY);
  sizeYperLayer->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  sizeYperLayer->DrawPerLayer();
  
  //---------------------------------------------------------------------------
  // Number of tracks in events passing different cuts
  //---------------------------------------------------------------------------
  
  int nSignal[kNsignals];
  int nShortSignal[kNsignals];
  int nShortAboveThresholdSignal[kNsignals];
  int nShortLowTotalSignal[kNsignals];
  int nShortLowTotalHighJetSignal[kNsignals];
  int nShortLowTotalHighMetSignal[kNsignals];
  int nBestSignal[kNsignals];
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    nSignal[iSig] =                     eventsSignal[iSig]->ApplyCuts(oneTrackOneJetEventCut,nullptr,nullptr)->size();
    nShortSignal[iSig] =                eventsSignal[iSig]->ApplyCuts(oneTrackOneJetEventCut,shortTrackCut,nullptr)->size();
    nShortAboveThresholdSignal[iSig] =  eventsSignal[iSig]->ApplyCuts(oneTrackOneJetEventCut,shortAboveTrasholdTrackCut,nullptr)->size();
    nShortLowTotalSignal[iSig] =        eventsSignal[iSig]->ApplyCuts(oneTrackOneJetEventCut,shortLowDedxTrackCut,nullptr)->size();
    nShortLowTotalHighJetSignal[iSig] = eventsSignal[iSig]->ApplyCuts(oneTrackOneJetEventCut, shortLowDedxTrackCut, highPtJetCut)->size();
    nShortLowTotalHighMetSignal[iSig] = eventsSignal[iSig]->ApplyCuts(highMetEventCut, shortLowDedxTrackCut, nullptr)->size();
    nBestSignal[iSig]                 = eventsSignal[iSig]->ApplyCuts(bestEventCut, bestTrackCut, bestJetCut)->size();
  }
  
  int nBackground[kNbackgrounds];
  int nShortBackground[kNbackgrounds];
  int nShortAboveThresholdBackground[kNbackgrounds];
  int nShortLowTotalBackground[kNbackgrounds];
  int nShortLowTotalHighJetBackground[kNbackgrounds];
  int nShortLowTotalHighMetBackground[kNbackgrounds];
  int nBestBackground[kNbackgrounds];
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    nBackground[iBck] =                     eventsBackground[iBck]->ApplyCuts(oneTrackOneJetEventCut,nullptr,nullptr)->size();
    nShortBackground[iBck] =                eventsBackground[iBck]->ApplyCuts(oneTrackOneJetEventCut,shortTrackCut,nullptr)->size();
    nShortAboveThresholdBackground[iBck] =  eventsBackground[iBck]->ApplyCuts(oneTrackOneJetEventCut,shortAboveTrasholdTrackCut,nullptr)->size();
    nShortLowTotalBackground[iBck] =        eventsBackground[iBck]->ApplyCuts(oneTrackOneJetEventCut,shortLowDedxTrackCut,nullptr)->size();
    nShortLowTotalHighJetBackground[iBck] = eventsBackground[iBck]->ApplyCuts(oneTrackOneJetEventCut, shortLowDedxTrackCut, highPtJetCut)->size();
    nShortLowTotalHighMetBackground[iBck] = eventsBackground[iBck]->ApplyCuts(highMetEventCut, shortLowDedxTrackCut, nullptr)->size();
    nBestBackground[iBck]                 = eventsBackground[iBck]->ApplyCuts(bestEventCut, bestTrackCut, bestJetCut)->size();
  }
    
  int nData=0, nShortData=0, nShortAboveThresholdData=0, nShortLowTotalData=0, nShortLowTotalHighJetData=0, nShortLowTotalHighMetData=0, nBestData=0;
  if(analyzeData){
    nData =                     eventsData->ApplyCuts(oneTrackOneJetEventCut,nullptr,nullptr)->size();
    nShortData =                eventsData->ApplyCuts(oneTrackOneJetEventCut,shortTrackCut,nullptr)->size();
    nShortAboveThresholdData =  eventsData->ApplyCuts(oneTrackOneJetEventCut,shortAboveTrasholdTrackCut,nullptr)->size();
    nShortLowTotalData =        eventsData->ApplyCuts(oneTrackOneJetEventCut,shortLowDedxTrackCut,nullptr)->size();
    nShortLowTotalHighJetData = eventsData->ApplyCuts(oneTrackOneJetEventCut, shortLowDedxTrackCut, highPtJetCut)->size();
    nShortLowTotalHighMetData = eventsData->ApplyCuts(highMetEventCut, shortLowDedxTrackCut, nullptr)->size();
    nBestData                 = eventsData->ApplyCuts(bestEventCut, bestTrackCut, bestJetCut)->size();
  }
  
  HistSet *nEvents = new HistSet("N events",10000,0,10e6);
  HistSet *nShortEvents = new HistSet("N events with short tracks (3 hits only) (\%)",100,0,1.6);
  HistSet *nShortAboveEvents = new HistSet("N short tracks (3 hits) above threshold (2.5 MeV / cluster) (\%)",100,0,1.6);
  HistSet *nShortLowTotalEvents = new HistSet("N events (cut set 1) (\%)",100,0,1.6);
  HistSet *nShortLowTotalHighJetEvents = new HistSet("N events (cut set 1 + jet cut) (\%)",100,0,1.6);
  HistSet *nShortLowTotalHighMetEvents = new HistSet("N events (cut set 1 + MET cut) (\%)",100,0,1.6);
  HistSet *nBest = new HistSet("N events (best cut set) (\%)",100,0,1.6);
  
  nShortEvents->SetShowNonZerBinPosX();
  nShortAboveEvents->SetShowNonZerBinPosX();
  nShortLowTotalEvents->SetShowNonZerBinPosX();
  nShortLowTotalHighJetEvents->SetShowNonZerBinPosX();
  nShortLowTotalHighMetEvents->SetShowNonZerBinPosX();
  nBest->SetShowNonZerBinPosX();
  
  nEvents->FillData(nData);
  nShortEvents->FillData(nShortData/(double)nData);
  nShortAboveEvents->FillData(nShortAboveThresholdData/(double)nData);
  nShortLowTotalEvents->FillData(nShortLowTotalData/(double)nData);
  nShortLowTotalHighJetEvents->FillData(nShortLowTotalHighJetData/(double)nData);
  nShortLowTotalHighMetEvents->FillData(nShortLowTotalHighMetData/(double)nData);
  nBest->FillData(nBestData/(double)nData);
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    nEvents->FillSignal((ESignal)iSig,nSignal[iSig]);
    nShortEvents->FillSignal((ESignal)iSig,nShortSignal[iSig]/(double)nSignal[iSig]);
    nShortAboveEvents->FillSignal((ESignal)iSig,nShortAboveThresholdSignal[iSig]/(double)nSignal[iSig]);
    nShortLowTotalEvents->FillSignal((ESignal)iSig,nShortLowTotalSignal[iSig]/(double)nSignal[iSig]);
    nShortLowTotalHighJetEvents->FillSignal((ESignal)iSig,nShortLowTotalHighJetSignal[iSig]/(double)nSignal[iSig]);
    nShortLowTotalHighMetEvents->FillSignal((ESignal)iSig,nShortLowTotalHighMetSignal[iSig]/(double)nSignal[iSig]);
    nBest->FillSignal((ESignal)iSig,nBestSignal[iSig]/(double)nSignal[iSig]);
  }
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    nEvents->FillBackground((EBackground)iBck, nBackground[iBck]);
    nShortEvents->FillBackground((EBackground)iBck, nShortBackground[iBck]/(double)nBackground[iBck]);
    nShortAboveEvents->FillBackground((EBackground)iBck, nShortAboveThresholdBackground[iBck]/(double)nBackground[iBck]);
    nShortLowTotalEvents->FillBackground((EBackground)iBck, nShortLowTotalBackground[iBck]/(double)nBackground[iBck]);
    nShortLowTotalHighJetEvents->FillBackground((EBackground)iBck, nShortLowTotalHighJetBackground[iBck]/(double)nBackground[iBck]);
    nShortLowTotalHighMetEvents->FillBackground((EBackground)iBck, nShortLowTotalHighMetBackground[iBck]/(double)nBackground[iBck]);
    nBest->FillBackground((EBackground)iBck, nBestBackground[iBck]/(double)nBackground[iBck]);
  }
    
  
  
  TCanvas *canvasNtracks = new TCanvas("Number of tracks","Number of tracks",2880,1800);
  canvasNtracks->Divide(2,3);
  
//  nEvents->Draw(canvasNtracks, 1);
  nShortEvents->Draw(canvasNtracks, 1);
  nShortAboveEvents->Draw(canvasNtracks, 2);
  nShortLowTotalEvents->Draw(canvasNtracks, 3);
  nShortLowTotalHighJetEvents->Draw(canvasNtracks, 4);
  nShortLowTotalHighMetEvents->Draw(canvasNtracks, 5);
  nBest->Draw(canvasNtracks, 6);
  
  
  theApp.Run();
  return 0;
}



