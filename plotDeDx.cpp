#include "Helpers.hpp"
#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"

#include <TApplication.h>

//  string inFileNameSignal = "../jniedzie/mcSignal/tree.root";
//string inFileNameBackground = "../jniedzie/mcBackground/tree.root";

string inFileNameSignal = "../adish/Signal/tree.root";
string inFileNameBackground = "../adish/Background/tree.root";
string inFileNameData = "../adish/Data/tree.root";

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  //---------------------------------------------------------------------------
  // Load MC and data files
  //---------------------------------------------------------------------------
  Events *eventsSignal = new Events(inFileNameSignal);
  Events *eventsBackground = new Events(inFileNameBackground);
  Events *eventsData = analyzeData ? new Events(inFileNameData) : nullptr;
  
  //---------------------------------------------------------------------------
  // Define event, track and jet cuts
  //---------------------------------------------------------------------------
  EventCut *highMetEventCut = new EventCut(EventCut::kMet100GeV);
  
  TrackCut *shortTrackCut = new TrackCut(TrackCut::kShort);
  TrackCut *shortAboveTrasholdTrackCut = new TrackCut(TrackCut::kShortAboveThreshold);
  TrackCut *shortLowDedxTrackCut = new TrackCut(TrackCut::kShortLowTotalDEdx);
  
  JetCut *highPtJetCut = new JetCut(JetCut::kHighPt);
  
  //---------------------------------------------------------------------------
  // Create standard per event, per track and per jet plots
  //---------------------------------------------------------------------------
  HistSet *nVertices = new HistSet(HistSet::kNvertices);
  HistSet *nIsoTrack = new HistSet(HistSet::kNisoTracks);
  HistSet *nJet      = new HistSet(HistSet::kNjets);
  HistSet *nJet30    = new HistSet(HistSet::kNjets30);
  HistSet *nJet30a   = new HistSet(HistSet::kNjets30a);
  HistSet *nMetSumEt = new HistSet(HistSet::kMetSumEt);
  HistSet *nMetPt    = new HistSet(HistSet::kMetPt);
  HistSet *nMetMass  = new HistSet(HistSet::kMetMass);
  HistSet *nMetEta   = new HistSet(HistSet::kMetEta);
  HistSet *nMetPhi   = new HistSet(HistSet::kMetPhi);
  
  HistSet *nClustersPerTrack = new HistSet(HistSet::kTrackNclusters);
  HistSet *totalDeDx = new HistSet(HistSet::kTrackTotalDedx);
  
  HistSet *totalDeDxByNclusters = new HistSet(HistSet::kTrackDedxPerCluster);
  HistSet *pt = new HistSet(HistSet::kTrackPt);
  HistSet *eta = new HistSet(HistSet::kTrackEta);
  HistSet *phi = new HistSet(HistSet::kTrackPhi);
  HistSet *caloEm = new HistSet(HistSet::kTrackCaloEm);
  HistSet *caloHad = new HistSet(HistSet::kTrackCaloHad);
  
  HistSet *dxy = new HistSet(HistSet::kTrackDxy);
  HistSet *dz = new HistSet(HistSet::kTrackDz);
  HistSet *charge = new HistSet(HistSet::kTrackCharge);
  HistSet *mass = new HistSet(HistSet::kTrackMass);
  HistSet *pid = new HistSet(HistSet::kTrackPid);
  
  HistSet *jet_pt = new HistSet(HistSet::kJetPt);
  HistSet *jet_eta = new HistSet(HistSet::kJetEta);
  HistSet *jet_phi = new HistSet(HistSet::kJetPhi);
  
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
  
  HistSet *dedxPerLayer = new HistSet(HistSet::kDedx);
  dedxPerLayer->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  dedxPerLayer->DrawPerLayer();
  
  HistSet *sizeXperLayer = new HistSet(HistSet::kSizeX);
  sizeXperLayer->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  sizeXperLayer->DrawPerLayer();
  
  HistSet *sizeYperLayer = new HistSet(HistSet::kSizeY);
  sizeYperLayer->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  sizeYperLayer->DrawPerLayer();
  
  //---------------------------------------------------------------------------
  // Number of tracks in events passing different cuts
  //---------------------------------------------------------------------------
  
  int nSignal =                     eventsSignal->SizeNonEmpty();
  int nShortSignal =                eventsSignal->ApplyTrackCut(shortTrackCut)->SizeNonEmpty();
  int nShortAboveThresholdSignal =  eventsSignal->ApplyTrackCut(shortAboveTrasholdTrackCut)->SizeNonEmpty();
  int nShortLowTotalSignal =        eventsSignal->ApplyTrackCut(shortLowDedxTrackCut)->SizeNonEmpty();
  int nShortLowTotalHighJetSignal = eventsSignal->ApplyCuts(nullptr, shortLowDedxTrackCut, highPtJetCut)->SizeNonEmpty();
  int nShortLowTotalHighMetSignal = eventsSignal->ApplyCuts(highMetEventCut, shortLowDedxTrackCut, nullptr)->SizeNonEmpty();
  
  int nBackground =                     eventsBackground->SizeNonEmpty();
  int nShortBackground = 	              eventsBackground->ApplyTrackCut(shortTrackCut)->SizeNonEmpty();
  int nShortAboveThresholdBackground =  eventsBackground->ApplyTrackCut(shortAboveTrasholdTrackCut)->SizeNonEmpty();
  int nShortLowTotalBackground =        eventsBackground->ApplyTrackCut(shortLowDedxTrackCut)->SizeNonEmpty();
  int nShortLowTotalHighJetBackground = eventsBackground->ApplyCuts(nullptr, shortLowDedxTrackCut, highPtJetCut)->SizeNonEmpty();
  int nShortLowTotalHighMetBackground = eventsBackground->ApplyCuts(highMetEventCut, shortLowDedxTrackCut, nullptr)->SizeNonEmpty();
  
  int nData=0, nShortData=0, nShortAboveThresholdData=0, nShortLowTotalData=0, nShortLowTotalHighJetData=0, nShortLowTotalHighMetData=0;
  if(analyzeData){
    nData =                     eventsData->SizeNonEmpty();
    nShortData =                eventsData->ApplyTrackCut(shortTrackCut)->SizeNonEmpty();
    nShortAboveThresholdData =  eventsData->ApplyTrackCut(shortAboveTrasholdTrackCut)->SizeNonEmpty();
    nShortLowTotalData =        eventsData->ApplyTrackCut(shortLowDedxTrackCut)->SizeNonEmpty();
    nShortLowTotalHighJetData = eventsData->ApplyCuts(nullptr, shortLowDedxTrackCut, highPtJetCut)->SizeNonEmpty();
    nShortLowTotalHighMetData = eventsData->ApplyCuts(highMetEventCut, shortLowDedxTrackCut, nullptr)->SizeNonEmpty();
  }
  
  HistSet *nEvents = new HistSet("N events",10000,0,10e6);
  HistSet *nShortEvents = new HistSet("N events with short tracks (3 hits only) (\%)",100,0,1.0);
  HistSet *nShortAboveEvents = new HistSet("N short tracks (3 hits) above threshold (2.5 MeV / cluster) (\%)",100,0,1.0);
  HistSet *nShortLowTotalEvents = new HistSet("N events (cut set 1) (\%)",100,0,1.0);
  HistSet *nShortLowTotalHighJetEvents = new HistSet("N events (cut set 1 + jet cut) (\%)",100,0,1.0);
  HistSet *nShortLowTotalHighMetEvents = new HistSet("N events (cut set 1 + MET cut) (\%)",100,0,1.0);
  
  nShortEvents->SetShowNonZerBinPosX();
  nShortAboveEvents->SetShowNonZerBinPosX();
  nShortLowTotalEvents->SetShowNonZerBinPosX();
  nShortLowTotalHighJetEvents->SetShowNonZerBinPosX();
  nShortLowTotalHighMetEvents->SetShowNonZerBinPosX();
  
  nEvents->FillSignal(nSignal);
  nEvents->FillBackground(nBackground);
  nEvents->FillData(nData);
  
  nShortEvents->FillSignal(nShortSignal/(double)nSignal);
  nShortEvents->FillBackground(nShortBackground/(double)nBackground);
  nShortEvents->FillData(nShortData/(double)nData);
  
  nShortAboveEvents->FillSignal(nShortAboveThresholdSignal/(double)nSignal);
  nShortAboveEvents->FillBackground(nShortAboveThresholdBackground/(double)nBackground);
  nShortAboveEvents->FillData(nShortAboveThresholdData/(double)nData);
  
  nShortLowTotalEvents->FillSignal(nShortLowTotalSignal/(double)nSignal);
  nShortLowTotalEvents->FillBackground(nShortLowTotalBackground/(double)nBackground);
  nShortLowTotalEvents->FillData(nShortLowTotalData/(double)nData);
  
  nShortLowTotalHighJetEvents->FillSignal(nShortLowTotalHighJetSignal/(double)nSignal);
  nShortLowTotalHighJetEvents->FillBackground(nShortLowTotalHighJetBackground/(double)nBackground);
  nShortLowTotalHighJetEvents->FillData(nShortLowTotalHighJetData/(double)nData);
  
  nShortLowTotalHighMetEvents->FillSignal(nShortLowTotalHighMetSignal/(double)nSignal);
  nShortLowTotalHighMetEvents->FillBackground(nShortLowTotalHighMetBackground/(double)nBackground);
  nShortLowTotalHighMetEvents->FillData(nShortLowTotalHighMetData/(double)nData);
  
  TCanvas *canvasNtracks = new TCanvas("Number of tracks","Number of tracks",2880,1800);
  canvasNtracks->Divide(2,3);
  
//  nEvents->Draw(canvasNtracks, 1);
  nShortEvents->Draw(canvasNtracks, 1);
  nShortAboveEvents->Draw(canvasNtracks, 2);
  nShortLowTotalEvents->Draw(canvasNtracks, 3);
  nShortLowTotalHighJetEvents->Draw(canvasNtracks, 4);
  nShortLowTotalHighMetEvents->Draw(canvasNtracks, 5);
  
  
  theApp.Run();
  return 0;
}



