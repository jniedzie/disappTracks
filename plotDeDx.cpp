#include "Helpers.hpp"
#include "Event.hpp"
#include "TrackCut.hpp"
#include "HistSet.hpp"

#include <TApplication.h>

const bool analyzeData = false;


int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
//  string inFileNameSignal = "../jniedzie/mcSignal/tree.root";
  string inFileNameSignal = "../adish/Signal/tree.root";
  
  string inFileNameBackground = "../jniedzie/mcBackground/tree.root";
//  const char *inFileNameBackground = "../adish/Background/tree.root";
  string inFileNameData = "../adish/Data/tree.root";
  
  
  TrackCut *shortTrackCut = new TrackCut(TrackCut::kShort);
  TrackCut *shortAboveTrasholdTrackCut = new TrackCut(TrackCut::kShortAboveThreshold);
  
  Events *eventsSignal = new Events(inFileNameSignal);
  int nTracksSignal = eventsSignal->GetNtracks();
  int nShortTracksSignal = eventsSignal->ApplyTrackCut(shortTrackCut)->GetNtracks();
  int nShortTracksAboveThresholdSignal = eventsSignal->ApplyTrackCut(shortAboveTrasholdTrackCut)->GetNtracks();
  
  Events *eventsBackground = new Events(inFileNameBackground);
  int nTracksBackground = eventsBackground->GetNtracks();
  int nShortTracksBackground = eventsBackground->ApplyTrackCut(shortTrackCut)->GetNtracks();
  int nShortTracksAboveThresholdBackground = eventsBackground->ApplyTrackCut(shortAboveTrasholdTrackCut)->GetNtracks();
  
  int nTracksData=0, nShortTracksData=0, nShortTracksAboveThresholdData=0;
  if(analyzeData){
    Events *eventsData = new Events(inFileNameData);
    nTracksData = eventsData->GetNtracks();
    nShortTracksData = eventsData->ApplyTrackCut(shortTrackCut)->GetNtracks();
    nShortTracksAboveThresholdData = eventsData->ApplyTrackCut(shortAboveTrasholdTrackCut)->GetNtracks();
  }
  
  HistSet *nTracks = new HistSet("N tracks",10000,0,10000);
  HistSet *nShortTracks = new HistSet("N short tracks (\%)",100,0,1.0);
  HistSet *nShortTracksAbove = new HistSet("N short tracks above threshold (\%)",100,0,1.0);
  
  nTracks->FillSignal(nTracksSignal);
  nTracks->FillBackground(nTracksBackground);
  nTracks->FillData(nTracksData);
  
  nShortTracks->FillSignal(nShortTracksSignal/(double)nTracksSignal);
  nShortTracks->FillBackground(nShortTracksBackground/(double)nTracksBackground);
  nShortTracks->FillData(nShortTracksData/(double)nTracksData);
  
  nShortTracksAbove->FillSignal(nShortTracksAboveThresholdSignal/(double)nTracksSignal);
  nShortTracksAbove->FillBackground(nShortTracksAboveThresholdBackground/(double)nTracksBackground);
  nShortTracksAbove->FillData(nShortTracksAboveThresholdData/(double)nTracksData);
  
  TCanvas *canvasNtracks = new TCanvas("Number of tracks","Number of tracks",2880,1800);
  canvasNtracks->Divide(2,2);
  
  nTracks->Draw(canvasNtracks, 1);
  nShortTracks->Draw(canvasNtracks, 2);
  nShortTracksAbove->Draw(canvasNtracks, 3);
  
  // Per layer plots
  HistSet *dedxPerLayer = new HistSet();
  dedxPerLayer->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kDedx);
  dedxPerLayer->DrawPerLayer(HistSet::kDedx);
  
  HistSet *sizeXperLayer = new HistSet();
  sizeXperLayer->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kSizeX);
  sizeXperLayer->DrawPerLayer(HistSet::kSizeX);
  
  HistSet *sizeYperLayer = new HistSet();
  sizeYperLayer->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kSizeY);
  sizeYperLayer->DrawPerLayer(HistSet::kSizeY);
  
  
  // Global plots
  HistSet *nClustersPerTrack = new HistSet("N clusters per track",20,0,20);
  HistSet *totalDeDx = new HistSet("total dedx per track",50,0,200);
  HistSet *dedxPerCluster = new HistSet("dedx per cluster",50,0,20);
  HistSet *totalDeDxByNclusters = new HistSet("total dedx per track / n clusters",50,0,20);
  HistSet *pt = new HistSet("pt dist",50,0,1000);
  HistSet *eta = new HistSet("eta dist",50,-3,3);
  HistSet *phi = new HistSet("phi dist",50,-3.5,3.5);
  HistSet *caloEm = new HistSet("EM calo energy",100,0,500);
  HistSet *caloHad = new HistSet("Hadron calo energy",100,0,500);
  
  HistSet *dxy = new HistSet("Displacement in XY",100,-0.02,0.02);
  HistSet *dz = new HistSet("Displacement in Z",100,-0.02,0.02);
  HistSet *charge = new HistSet("Charge",100,-10,10);
  HistSet *mass = new HistSet("Mass",500,0,0.25);
  HistSet *pid = new HistSet("PDG PID",441,-220,220);
  
  HistSet *jet_pt = new HistSet("Jet pt",100,0,1000);
  HistSet *jet_eta = new HistSet("Jet eta",50,-3,3);
  HistSet *jet_phi = new HistSet("Jet phi",50,-3.5,3.5);
  
  jet_pt->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kJetPt);
  jet_eta->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kJetEta);
  jet_phi->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kJetPhi);
  
  // Fill signal histograms
  for(int iEvent=0;iEvent<eventsSignal->size();iEvent++){
    Event *event = eventsSignal->At(iEvent);
    
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      Track *track = event->GetTrack(iTrack);
      
      nClustersPerTrack->FillSignal(track->GetNclusters());
      totalDeDx->FillSignal(track->GetTotalDedx());
      totalDeDxByNclusters->FillSignal(track->GetTotalDedx()/track->GetNclusters());
      pt->FillSignal(track->GetPt());
      eta->FillSignal(track->GetEta());
      phi->FillSignal(track->GetPhi());
      caloEm->FillSignal(track->GetCaloEmEnergy());
      caloHad->FillSignal(track->GetCaloHadEnergy());
      dxy->FillSignal(track->GetDxy());
      dz->FillSignal(track->GetDz());
      charge->FillSignal(track->GetCharge());
      mass->FillSignal(track->GetMass());
      pid->FillSignal(track->GetPid());
      
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        if(track->GetDeDxInLayer(iLayer) > 0.0000001){
          dedxPerCluster->FillSignal(track->GetDeDxInLayer(iLayer));
        }
      }
    }
  }
  
  // Fill background histograms
  for(int iEvent=0;iEvent<eventsBackground->size();iEvent++){
    Event *event = eventsBackground->At(iEvent);
    
    for(int iTrack=0;iTrack<event->GetNtracks();iTrack++){
      Track *track = event->GetTrack(iTrack);
      
      nClustersPerTrack->FillBackground(track->GetNclusters());
      totalDeDx->FillBackground(track->GetTotalDedx());
      totalDeDxByNclusters->FillBackground(track->GetTotalDedx()/track->GetNclusters());
      pt->FillBackground(track->GetPt());
      eta->FillBackground(track->GetEta());
      phi->FillBackground(track->GetPhi());
      caloEm->FillBackground(track->GetCaloEmEnergy());
      caloHad->FillBackground(track->GetCaloHadEnergy());
      dxy->FillBackground(track->GetDxy());
      dz->FillBackground(track->GetDz());
      charge->FillBackground(track->GetCharge());
      mass->FillBackground(track->GetMass());
      pid->FillBackground(track->GetPid());
      
      
      for(int iLayer=0;iLayer<nLayers;iLayer++){
        if(track->GetDeDxInLayer(iLayer) > 0.0000001){
          dedxPerCluster->FillBackground(track->GetDeDxInLayer(iLayer));
        }
      }
    }
  }
  
  // Plot histograms
  TCanvas *c2 = new TCanvas("c2","c2",2880,1800);
  c2->Divide(3,3);
  
  nClustersPerTrack->Draw(c2,1);
  totalDeDx->Draw(c2,2);
  dedxPerCluster->Draw(c2,3);
  totalDeDxByNclusters->Draw(c2,4);
  pt->Draw(c2,5);
  eta->Draw(c2,6);
  phi->Draw(c2,7);
  caloEm->Draw(c2,8);
  caloHad->Draw(c2,9);
  
  TCanvas *c3 = new TCanvas("c3","c3",2880,1800);
  c3->Divide(2,3);
  dxy->Draw(c3,1);
  dz->Draw(c3,2);
  charge->Draw(c3,3);
  mass->Draw(c3,4);
  pid->Draw(c3,5);
  jet_pt->Draw(c3, 6);
  
  theApp.Run();
  
  return 0;
}



