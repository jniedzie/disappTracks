#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"

#include <TApplication.h>

void LoadEventsFromFiles(vector<Events*> &eventsSignal, vector<Events*> &eventsBackground, vector<Events*> &eventsData)
{
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]){
      eventsData.push_back(nullptr);
    }
    else{
      eventsData.push_back(new Events(inFileNameData[iData], Events::kData, maxNeventsData));
    }
  }
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]){
      eventsSignal.push_back(nullptr);
    }
    else{
      eventsSignal.push_back(new Events(inFileNameSignal[iSig], Events::kSignal, maxNeventsSignal,(ESignal)iSig));
    }
  }
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]){
      eventsBackground.push_back(nullptr);
    }
    else{
      eventsBackground.push_back(new Events());
      
      for(string path : inFileNameBackground[iBck]){
        eventsBackground[iBck]->AddEventsFromFile(path, Events::kBackground, maxNeventsBackground);
      }
    }
  }
}

void ApplyCuts(vector<Events*> &eventsSignal, vector<Events*> &eventsBackground, vector<Events*> &eventsData,
               EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut)
{
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    eventsSignal[iSig] = eventsSignal[iSig]->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  }
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]) continue;
    eventsBackground[iBck] = eventsBackground[iBck]->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  }
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]) continue;
    eventsData[iData] = eventsData[iData]->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  }
}

void DrawStandardPlots(vector<Events*> &eventsSignal, vector<Events*> &eventsBackground, vector<Events*> &eventsData)
{
  // Create standard per event, per track and per jet plots
  
  HistSet *nVertices    = new HistSet(kNvertices);
  HistSet *nIsoTrack    = new HistSet(kNisoTracks);
  HistSet *nJet         = new HistSet(kNjets);
  HistSet *nJet30       = new HistSet(kNjets30);
  HistSet *nJet30a      = new HistSet(kNjets30a);
  HistSet *nMetSumEt    = new HistSet(kMetSumEt);
  HistSet *nMetPt       = new HistSet(kMetPt);
  HistSet *nMetMass     = new HistSet(kMetMass);
  HistSet *nMetEta      = new HistSet(kMetEta);
  HistSet *nMetPhi      = new HistSet(kMetPhi);
  HistSet *nMetJetDphi  = new HistSet(kMetJetDphi);
  
  HistSet *nClustersPerTrack  = new HistSet(kTrackNclusters);
  HistSet *totalDeDx          = new HistSet(kTrackTotalDedx);
  
  HistSet *totalDeDxByNclusters = new HistSet(kTrackDedxPerCluster);
  HistSet *pt                   = new HistSet(kTrackPt);
  HistSet *eta                  = new HistSet(kTrackEta);
  HistSet *phi                  = new HistSet(kTrackPhi);
  HistSet *caloEm               = new HistSet(kTrackCaloEm);
  HistSet *caloHad              = new HistSet(kTrackCaloHad);
  
  HistSet *dxy    = new HistSet(kTrackDxy);
  HistSet *dz     = new HistSet(kTrackDz);
  HistSet *charge = new HistSet(kTrackCharge);
  HistSet *mass   = new HistSet(kTrackMass);
  HistSet *pid    = new HistSet(kTrackPid);
  
  HistSet *jet_pt   = new HistSet(kJetPt);
  HistSet *jet_eta  = new HistSet(kJetEta);
  HistSet *jet_phi  = new HistSet(kJetPhi);
  
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
  nMetJetDphi->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  
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
  canvasEvents->Divide(3,3);
  
  nVertices->Draw(canvasEvents,1);
  nIsoTrack->Draw(canvasEvents,2);
  nJet->Draw(canvasEvents,3);
  //  nJet30->Draw(canvasEvents,4);
  //  nJet30a->Draw(canvasEvents,5);
  jet_pt->Draw(canvasEvents, 4);
  nMetSumEt->Draw(canvasEvents,5);
  nMetPt->Draw(canvasEvents,6);
  nMetJetDphi->Draw(canvasEvents,7);
  //  nMetMass->Draw(canvasEvents,8);
  //  nMetEta->Draw(canvasEvents,9);
  //  nMetPhi->Draw(canvasEvents,10);
  
  TCanvas *canvasTrack = new TCanvas("Tracks","Tracks",2880,1800);
  canvasTrack->Divide(2,3);
  
  nClustersPerTrack->Draw(canvasTrack,1);
  totalDeDx->Draw(canvasTrack,2);
  totalDeDxByNclusters->Draw(canvasTrack,3);
  pt->Draw(canvasTrack,4);
  caloEm->Draw(canvasTrack,5);
  caloHad->Draw(canvasTrack,6);
  //  eta->Draw(canvasTrack,5);
  //  phi->Draw(canvasTrack,6);
  //  charge->Draw(canvasTrack,3);
  //  mass->Draw(canvasTrack,4);
  //  pid->Draw(canvasTrack,5);
  //  dxy->Draw(canvasTrack,1);
  //  dz->Draw(canvasTrack,2);
  
  //  TCanvas *canvasJets = new TCanvas("Jets","Jets",2880,1800);
  //  canvasJets->Divide(2,2);
  
  //  jet_eta->Draw(canvasJets, 2);
  //  jet_phi->Draw(canvasJets, 3);
}

void DrawPerLayerPlots(vector<Events*> &eventsSignal, vector<Events*> &eventsBackground, vector<Events*> &eventsData)
{
  HistSet *dedxPerLayer = new HistSet(kDedx);
  dedxPerLayer->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  dedxPerLayer->DrawPerLayer();
  
  //  HistSet *sizeXperLayer = new HistSet(kSizeX);
  //  sizeXperLayer->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  //  sizeXperLayer->DrawPerLayer();
  //
  //  HistSet *sizeYperLayer = new HistSet(kSizeY);
  //  sizeYperLayer->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  //  sizeYperLayer->DrawPerLayer();
}

void PrintSignalToBackground(
                  vector<Events*> &eventsSignal, vector<Events*> &eventsBackground, vector<Events*> &eventsData,
                  vector<int> initNsignal)
{
  int nEvents[kNsignals];
  int nSignalTotal=0, nBackgroundTotal=0;
  int nSignal[kNsignals], nBackground[kNbackgrounds];
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    nEvents[iSig] = initNsignal[iSig];
    nSignal[iSig] = eventsSignal[iSig]->size();
    nSignalTotal += nSignal[iSig];
    
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      if(!runBackground[iBck]) continue;
      nEvents[iSig] += eventsBackground[iBck]->size();
    }
    for(int iData=0;iData<kNdata;iData++){
      if(!runData[iData]) continue;
      nEvents[iSig] += eventsData[iData]->size();
    }
  }
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]) continue;
    
    if(eventsBackground[iBck]){
      nBackground[iBck] = eventsBackground[iBck]->size();
      nBackgroundTotal += nBackground[iBck];
    }
    else{
      nBackground[iBck] = 0;
    }
  }
  
  int nData=0,nDataTotal=0;
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]) continue;
    nData = eventsData[iData]->size();
  }
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    if(printHeaders) cout<<signalTitle[iSig]<<"\t";
    cout<<nSignal[iSig]/(double)nEvents[iSig]<<endl;
    
    if(printHeaders) cout<<"S/B:\t";
    cout<<nSignal[iSig]/(double)nBackgroundTotal<<endl;
  }
  
  if(printHeaders) cout<<"Bck total:\t";
  cout<<nBackgroundTotal/(double)nEvents[0]<<endl;
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]) continue;
    
    if(printHeaders) cout<<backgroundTitle[iBck]<<"\t";
    cout<<nBackground[iBck]/(double)nEvents[0]<<endl;
  }
  
  if(printHeaders) cout<<"N data/total:";
  cout<<nData/(double)nEvents[0]<<endl;
}

int main(int argc, char* argv[])
{
  TApplication theApp("App", &argc, argv);
  
  // All events with initial cuts only
  vector<Events*> eventsSignal;
  vector<Events*> eventsBackground;
  vector<Events*> eventsData;
  
  LoadEventsFromFiles(eventsSignal, eventsBackground, eventsData);
  
  //---------------------------------------------------------------------------
  // Define event, track and jet cuts
  //---------------------------------------------------------------------------
  
  unsigned int eventCutOptions =
    EventCut::kOneTrack
  | EventCut::kOneJet;
  
  unsigned int trackCutOptions =
    TrackCut::kLowCalo
    | TrackCut::kLowDEdx;
//  | TrackCut::kHighPt;
//  | TrackCut::kMedium;
  
  unsigned int jetCutOptions =
  JetCut::kPt200GeV;
  
  EventCut  *initialEventCut = new EventCut((EventCut::ECut)eventCutOptions);
  TrackCut  *initialTrackCut = new TrackCut(TrackCut::kEmpty);
  JetCut    *initialJetCut   = new JetCut(JetCut::kEmpty);
  
  ApplyCuts(eventsSignal, eventsBackground, eventsData,
            initialEventCut, initialTrackCut, initialJetCut, nullptr);
  
  
  
  DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
  DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
  
  vector<int> initNsignal;
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig])  initNsignal.push_back(0);
    else                  initNsignal.push_back(eventsSignal[iSig]->size());
  }
  
  PrintSignalToBackground(eventsSignal, eventsBackground, eventsData, initNsignal);
  
  unsigned int bestEventCutOptions =
  EventCut::kOneTrack
  | EventCut::kOneJet
  | EventCut::kMet100GeV;
  
  EventCut  *bestEventCut = new EventCut((EventCut::ECut)bestEventCutOptions);
  TrackCut  *bestTrackCut = new TrackCut((TrackCut::ECut)trackCutOptions);
  JetCut    *bestJetCut   = new JetCut((JetCut::ECut)jetCutOptions);
  
  vector<Events*> eventsSignalBestCut = eventsSignal;
  vector<Events*> eventsBackgroundBestCut = eventsBackground;
  vector<Events*> eventsDataBestCut = eventsData;
  
  ApplyCuts(eventsSignalBestCut, eventsBackgroundBestCut,eventsDataBestCut , bestEventCut, bestTrackCut, bestJetCut, nullptr);
  
  PrintSignalToBackground(eventsSignalBestCut, eventsBackgroundBestCut, eventsDataBestCut, initNsignal);
  

  
  theApp.Run();
  return 0;
}



