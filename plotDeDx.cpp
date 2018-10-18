#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"

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

void DrawStandardPlots(vector<Events*> &eventsSignal, vector<Events*> &eventsBackground, vector<Events*> &eventsData, string prefix="")
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
  HistSet *missingOuterTracker  = new HistSet(kTrackMissingOuterTrackerHits);
  HistSet *pixelHits            = new HistSet(kTrackPixelHits);
  HistSet *trackerHits          = new HistSet(kTrackTrackerHits);
  HistSet *isolation            = new HistSet(kTrackRelativeIsolation);
  HistSet *trackMetDphi         = new HistSet(kTrackMetDphi);
  HistSet *dedx                 = new HistSet(kTrackDedxPerHit);
  
  HistSet *dxy    = new HistSet(kTrackDxy);
  HistSet *dz     = new HistSet(kTrackDz);
  HistSet *charge = new HistSet(kTrackCharge);
  HistSet *mass   = new HistSet(kTrackMass);
  HistSet *pid    = new HistSet(kTrackPid);
  
  HistSet *jet_pt   = new HistSet(kJetPt);
  HistSet *jet_eta  = new HistSet(kJetEta);
  HistSet *jet_phi  = new HistSet(kJetPhi);
  HistSet *jetTrackDr  = new HistSet(kJetTrackDr);
  
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
  caloHad->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  missingOuterTracker->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  pixelHits->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  trackerHits->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  isolation->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  trackMetDphi->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  dedx->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  
  dxy->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  dz->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  charge->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  mass->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  pid->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  
  jet_pt->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  jet_eta->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  jet_phi->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  jetTrackDr->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  
  
  // Plot histograms
  TCanvas *canvasEvents = new TCanvas((prefix+"Events").c_str(),(prefix+"Events").c_str(),2880,1800);
  canvasEvents->Divide(3,2);
  
  nVertices->Draw(canvasEvents,1);
  nIsoTrack->Draw(canvasEvents,2);
  nJet->Draw(canvasEvents,3);
  jet_pt->Draw(canvasEvents, 4);
  nMetPt->Draw(canvasEvents,5);
  nMetJetDphi->Draw(canvasEvents,6);
  
  TCanvas *canvasTrack = new TCanvas((prefix+"Tracks").c_str(),(prefix+"Tracks").c_str(),2880,1800);
  canvasTrack->Divide(4,3);
  
  pt->Draw(canvasTrack,1);
  caloEm->Draw(canvasTrack,2);
  caloHad->Draw(canvasTrack,3);
  missingOuterTracker->Draw(canvasTrack,4);
  pixelHits->Draw(canvasTrack,5);
  trackerHits->Draw(canvasTrack,6);
  isolation->Draw(canvasTrack,7);
  trackMetDphi->Draw(canvasTrack,8);
  eta->Draw(canvasTrack,9);
  dedx->Draw(canvasTrack,10);
  
  TCanvas *canvasJets = new TCanvas("Jets","Jets",2880,1800);
  canvasJets->Divide(2,2);
  
  jetTrackDr->Draw(canvasJets, 1);
  jet_eta->Draw(canvasJets, 2);
  jet_phi->Draw(canvasJets, 3);
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

void PrintYields(vector<Events*> &eventsSignal, vector<Events*> &eventsBackground, vector<Events*> &eventsData)
{
  int nBackgroundTotal=0;
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck] || !eventsBackground[iBck]) continue;
    nBackgroundTotal += eventsBackground[iBck]->weightedSize();
    
    if(printHeaders) cout<<backgroundTitle[iBck]<<"\t";
    cout<<eventsBackground[iBck]->weightedSize()<<endl;
  }

  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    if(printHeaders) cout<<signalTitle[iSig]<<"\tN events:\t";
    cout<<eventsSignal[iSig]->weightedSize()<<endl;
  }
    
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    if(printHeaders) cout<<signalTitle[iSig]<<"\tS/B:\t";
    cout<<eventsSignal[iSig]->weightedSize()/(double)nBackgroundTotal<<endl;
  }
}

int main(int argc, char* argv[])
{
  TApplication *theApp = new TApplication("App", &argc, argv);
  
  // All events with initial cuts only
  vector<Events*> eventsSignal, eventsBackground, eventsData;
  LoadEventsFromFiles(eventsSignal, eventsBackground, eventsData);
  
  cout<<"\n\nInitial yields"<<endl;
  PrintYields(eventsSignal, eventsBackground, eventsData);
  
  //---------------------------------------------------------------------------
  // Level 0
  //---------------------------------------------------------------------------
  
  EventCut  *eventCut_L0 = new EventCut();
  TrackCut  *trackCut_L0 = new TrackCut();
  JetCut    *jetCut_L0   = new JetCut();
  
  eventCut_L0->SetNtracks(1, 999999);
  eventCut_L0->SetMinNjets(1);
  eventCut_L0->SetMaxNmuons(0);
  eventCut_L0->SetMaxNtau(0);
  eventCut_L0->SetMaxNlepton(0);
  eventCut_L0->SetMinMetNoMuPt(200);
  eventCut_L0->SetRequireMetNoMuTrigger(true);
  eventCut_L0->SetRequirePassingAllFilters(true);
  eventCut_L0->SetHighJetMinPt(100);
  eventCut_L0->SetHighJetMaxEta(2.4);
  eventCut_L0->SetHighJetMaxNeHEF(0.8);
  eventCut_L0->SetHighJetMinChHEF(0.1);
  
  //  trackCut_L0->SetRequireSameNpixelHitsLayers(true);
  trackCut_L0->SetNmissingInnerPixel(0, 0);
  trackCut_L0->SetNmissingMiddleTracker(0, 0);
  trackCut_L0->SetNpixelLayers(2, 999999);
  trackCut_L0->SetMaxEta(2.1);
  
  jetCut_L0->SetPtRange(30, 999999);
  
  ApplyCuts(eventsSignal, eventsBackground, eventsData,eventCut_L0, trackCut_L0, jetCut_L0, nullptr);
  cout<<"\n\nYields after level 0 cuts"<<endl;
  PrintYields(eventsSignal, eventsBackground, eventsData);
  
  if(plotAfterLevel == 0){
    DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
    //  DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
  }
  
  //---------------------------------------------------------------------------
  // Level 1
  //---------------------------------------------------------------------------
  
  EventCut  *eventCut_L1 = eventCut_L0;
  TrackCut  *trackCut_L1 = trackCut_L0;
  JetCut    *jetCut_L1   = jetCut_L0;
  
  eventCut_L1->SetNtracks(1, 999999);
  eventCut_L1->SetMinNjets(1);
  eventCut_L1->SetHighJetMinPt(100);
  eventCut_L1->SetHighJetMaxEta(2.4);
  eventCut_L1->SetHighJetMaxNeHEF(0.8);
  eventCut_L1->SetHighJetMinChHEF(0.1);
  
  trackCut_L1->SetMaxRelativeIsolation(0.15);
  
  jetCut_L1->SetMinTrackDeltaR(0.2);
  
  ApplyCuts(eventsSignal, eventsBackground, eventsData,eventCut_L1, trackCut_L1, jetCut_L1, nullptr);
  cout<<"\n\nYields after level 1 cuts"<<endl;
  PrintYields(eventsSignal, eventsBackground, eventsData);
  
  if(plotAfterLevel == 1){
    DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
    //  DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
  }
  
  //---------------------------------------------------------------------------
  // Level 2
  //---------------------------------------------------------------------------
  
  EventCut  *eventCut_L2 = eventCut_L1;
  TrackCut  *trackCut_L2 = trackCut_L1;
  JetCut    *jetCut_L2   = jetCut_L1;
  
  eventCut_L2->SetMinMetPt(230);
  eventCut_L2->SetMinJetMetPhi(0.5);

  trackCut_L2->SetMaxEmCalo(0.5);
  trackCut_L2->SetMaxHadCalo(0.5);
  trackCut_L2->SetNmissingOuterTracker(1, 999999);
  
  ApplyCuts(eventsSignal, eventsBackground, eventsData,eventCut_L2, trackCut_L2, jetCut_L2, nullptr);
  cout<<"\n\nYields after level 2 cuts"<<endl;
  PrintYields(eventsSignal, eventsBackground, eventsData);
  
  if(plotAfterLevel == 2){
    DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
    //  DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
  }
  
  //---------------------------------------------------------------------------
  // Sub-categories
  //---------------------------------------------------------------------------
  
  // here select sub-categories of events (with 1 or two tracks, with 3 or 4 layers etc.)
  //  eventCut_L0->SetNtracks(2, 2);
  
  if(plotAfterLevel == 3){
    DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
    //  DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
  }
  
  vector<Events*> currentEventsSignal, currentEventsBackground, currentEventsData;
  
  if(interactive){
    string option;
    double value;

    EventCut  *currentEventCut  = new EventCut();
    TrackCut  *currentTrackCut  = new TrackCut();
    JetCut    *currentJetCut    = new JetCut();
    
    currentEventCut->SetMinNjets(1);
    currentEventCut->SetNtracks(1,999999);
    currentEventCut->SetHighJetMaxNeHEF(0.8);
    currentEventCut->SetHighJetMinChHEF(0.1);
    currentEventCut->SetHighJetMaxEta(2.4);
    
    bool stop = false;
    
    while(!stop){
      cin >> option >> value;
      
      if(option == "high_jet_pt_min")     currentEventCut->SetHighJetMinPt(value);
      else if(option == "met_pt_min")     currentEventCut->SetMinMetPt(value);
      else if(option == "min_dets")       currentTrackCut->SetNdets(value, currentTrackCut->GetMaxDets());
      else if(option == "max_dets")       currentTrackCut->SetNdets(currentTrackCut->GetMinDets(), value);
      else if(option == "max_em_calo")    currentTrackCut->SetMaxEmCalo(value);
      else if(option == "max_had_calo")   currentTrackCut->SetMaxHadCalo(value);
      else if(option == "break"){
        stop = true;
        break;
      }
      else{
        cout<<"Unknown option"<<endl;
      }
      
      cout<<"current min dets:"<<currentTrackCut->GetMinDets()<<"\tmax dets:"<<currentTrackCut->GetMaxDets()<<endl;
      
      currentEventsSignal.clear();
      currentEventsBackground.clear();
      currentEventsData.clear();
      
      currentEventsSignal = eventsSignal;
      currentEventsBackground = eventsBackground;
      currentEventsData = eventsData;
      
      ApplyCuts(currentEventsSignal, currentEventsBackground,currentEventsData,
                currentEventCut, currentTrackCut, currentJetCut, nullptr);
      
      PrintYields(currentEventsSignal, currentEventsBackground, currentEventsData);
    }
  }
  
  if(showPlots) theApp->Run();
  
  return 0;
}



