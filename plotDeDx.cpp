#include "Helpers.hpp"
#include "Event.hpp"
#include "TrackCut.hpp"
#include "HistSet.hpp"

#include <TApplication.h>

const bool analyzeData = false;

//  string inFileNameSignal = "../jniedzie/mcSignal/tree.root";
string inFileNameSignal = "../adish/Signal/tree.root";

string inFileNameBackground = "../jniedzie/mcBackground/tree.root";
//  const char *inFileNameBackground = "../adish/Background/tree.root";
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
  // Define track cuts
  //---------------------------------------------------------------------------
  TrackCut *shortTrackCut = new TrackCut(TrackCut::kShort);
  TrackCut *shortAboveTrasholdTrackCut = new TrackCut(TrackCut::kShortAboveThreshold);
  
  
  //---------------------------------------------------------------------------
  // Create standard per track and per jet plots
  //---------------------------------------------------------------------------
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
  
  
  nClustersPerTrack->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kTrackNclusters);
  totalDeDx->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kTrackTotalDedx);
  
  totalDeDxByNclusters->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kTrackDedxPerCluster);
  pt->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kTrackPt);
  eta->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kTrackEta);
  phi->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kTrackPhi);
  caloEm->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kTrackCaloEm);
  caloHad->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kTrackCaloHad);
  
  dxy->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kTrackDxy);
  dz->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kTrackDz);
  charge->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kTrackCharge);
  mass->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kTrackMass);
  pid->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kTrackPid);
  
  jet_pt->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kJetPt);
  jet_eta->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kJetEta);
  jet_phi->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kJetPhi);
  
  
  // Plot histograms
  TCanvas *canvasTrack = new TCanvas("Tracks","Tracks",2880,1800);
  canvasTrack->Divide(3,3);
  
  nClustersPerTrack->Draw(canvasTrack,1);
  totalDeDx->Draw(canvasTrack,2);
  totalDeDxByNclusters->Draw(canvasTrack,2);
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
  
  HistSet *dedxPerLayer = new HistSet();
  dedxPerLayer->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kDedx);
  dedxPerLayer->DrawPerLayer(HistSet::kDedx);
  
  HistSet *sizeXperLayer = new HistSet();
  sizeXperLayer->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kSizeX);
  sizeXperLayer->DrawPerLayer(HistSet::kSizeX);
  
  HistSet *sizeYperLayer = new HistSet();
  sizeYperLayer->FillFromEvents(eventsSignal, eventsBackground, nullptr, HistSet::kSizeY);
  sizeYperLayer->DrawPerLayer(HistSet::kSizeY);
  
  //---------------------------------------------------------------------------
  // Create other custom plots
  //---------------------------------------------------------------------------
  
  
  int nTracksSignal = eventsSignal->GetNtracks();
  int nShortTracksSignal = eventsSignal->ApplyTrackCut(shortTrackCut)->GetNtracks();
  int nShortTracksAboveThresholdSignal = eventsSignal->ApplyTrackCut(shortAboveTrasholdTrackCut)->GetNtracks();
  
  
  int nTracksBackground = eventsBackground->GetNtracks();
  int nShortTracksBackground = eventsBackground->ApplyTrackCut(shortTrackCut)->GetNtracks();
  int nShortTracksAboveThresholdBackground = eventsBackground->ApplyTrackCut(shortAboveTrasholdTrackCut)->GetNtracks();
  
  int nTracksData=0, nShortTracksData=0, nShortTracksAboveThresholdData=0;
  if(analyzeData){
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
  
  
  theApp.Run();
  return 0;
}



