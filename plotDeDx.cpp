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
  TCanvas *canvasEvents = new TCanvas((prefix+"Events").c_str(),(prefix+"Events").c_str(),2880,1800);
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
  
  TCanvas *canvasTrack = new TCanvas((prefix+"Tracks").c_str(),(prefix+"Tracks").c_str(),2880,1800);
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

void MakeComparison(vector<Events*> &eventsSignal1, vector<Events*> &eventsBackground1, vector<Events*> &eventsData1,
                    vector<Events*> &eventsSignal2, vector<Events*> &eventsBackground2, vector<Events*> &eventsData2)
{
  // Create standard per event, per track and per jet plots
  
  HistSet *nVertices1    = new HistSet(kNvertices);
  HistSet *nIsoTrack1    = new HistSet(kNisoTracks);
  HistSet *nJet1         = new HistSet(kNjets);
  HistSet *nMetSumEt1    = new HistSet(kMetSumEt);
  HistSet *nMetPt1       = new HistSet(kMetPt);
  HistSet *nMetJetDphi1  = new HistSet(kMetJetDphi);
  HistSet *nClustersPerTrack1  = new HistSet(kTrackNclusters);
  HistSet *totalDeDx1          = new HistSet(kTrackTotalDedx);
  HistSet *totalDeDxByNclusters1 = new HistSet(kTrackDedxPerCluster);
  HistSet *pt1                   = new HistSet(kTrackPt);
  HistSet *caloEm1               = new HistSet(kTrackCaloEm);
  HistSet *caloHad1              = new HistSet(kTrackCaloHad);
  HistSet *jet_pt1   = new HistSet(kJetPt);
  
  HistSet *nVertices2    = new HistSet(kNvertices);
  HistSet *nIsoTrack2    = new HistSet(kNisoTracks);
  HistSet *nJet2         = new HistSet(kNjets);
  HistSet *nMetSumEt2    = new HistSet(kMetSumEt);
  HistSet *nMetPt2       = new HistSet(kMetPt);
  HistSet *nMetJetDphi2  = new HistSet(kMetJetDphi);
  HistSet *nClustersPerTrack2  = new HistSet(kTrackNclusters);
  HistSet *totalDeDx2          = new HistSet(kTrackTotalDedx);
  HistSet *totalDeDxByNclusters2 = new HistSet(kTrackDedxPerCluster);
  HistSet *pt2                   = new HistSet(kTrackPt);
  HistSet *caloEm2               = new HistSet(kTrackCaloEm);
  HistSet *caloHad2              = new HistSet(kTrackCaloHad);
  HistSet *jet_pt2   = new HistSet(kJetPt);
  
  nIsoTrack1->FillFromEvents(eventsSignal1, eventsBackground1, eventsData1);
  nJet1->FillFromEvents(eventsSignal1, eventsBackground1, eventsData1);
  nMetSumEt1->FillFromEvents(eventsSignal1, eventsBackground1, eventsData1);
  nMetPt1->FillFromEvents(eventsSignal1, eventsBackground1, eventsData1);
  nMetJetDphi1->FillFromEvents(eventsSignal1, eventsBackground1, eventsData1);
  nClustersPerTrack1->FillFromEvents(eventsSignal1, eventsBackground1, eventsData1);
  pt1->FillFromEvents(eventsSignal1, eventsBackground1, eventsData1);
  caloEm1->FillFromEvents(eventsSignal1, eventsBackground1, eventsData1);
  caloHad1->FillFromEvents(eventsSignal1, eventsBackground1, eventsData1);
  
  nIsoTrack2->FillFromEvents(eventsSignal2, eventsBackground2, eventsData2);
  nJet2->FillFromEvents(eventsSignal2, eventsBackground2, eventsData2);
  nMetSumEt2->FillFromEvents(eventsSignal2, eventsBackground2, eventsData2);
  nMetPt2->FillFromEvents(eventsSignal2, eventsBackground2, eventsData2);
  nMetJetDphi2->FillFromEvents(eventsSignal2, eventsBackground2, eventsData2);
  nClustersPerTrack2->FillFromEvents(eventsSignal2, eventsBackground2, eventsData2);
  pt2->FillFromEvents(eventsSignal2, eventsBackground2, eventsData2);
  caloEm2->FillFromEvents(eventsSignal2, eventsBackground2, eventsData2);
  caloHad2->FillFromEvents(eventsSignal2, eventsBackground2, eventsData2);
  
  // Plot histograms
  TCanvas *canvasEvents = new TCanvas("Events comp.","Events comp",2880,1800);
  canvasEvents->Divide(4,4);
  
  nVertices1->Draw(canvasEvents,1);
  nVertices2->Draw(canvasEvents,2);
  
  nIsoTrack1->Draw(canvasEvents,3);
  nIsoTrack2->Draw(canvasEvents,4);
  
  nJet1->Draw(canvasEvents,5);
  nJet2->Draw(canvasEvents,6);
  
  jet_pt1->Draw(canvasEvents, 7);
  jet_pt2->Draw(canvasEvents, 8);
  
  nMetSumEt1->Draw(canvasEvents,9);
  nMetSumEt2->Draw(canvasEvents,10);
  
  nMetPt1->Draw(canvasEvents,11);
  nMetPt2->Draw(canvasEvents,12);
  
  nMetJetDphi1->Draw(canvasEvents,13);
  nMetJetDphi2->Draw(canvasEvents,14);
  
  TCanvas *canvasTrack = new TCanvas("Tracks comp.","Tracks comp.",2880,1800);
  canvasTrack->Divide(4,3);
  
  nClustersPerTrack1->Draw(canvasTrack,1);
  nClustersPerTrack2->Draw(canvasTrack,2);
  
  totalDeDx1->Draw(canvasTrack,3);
  totalDeDx2->Draw(canvasTrack,4);
  
  totalDeDxByNclusters1->Draw(canvasTrack,5);
  totalDeDxByNclusters2->Draw(canvasTrack,6);
  
  pt1->Draw(canvasTrack,7);
  pt2->Draw(canvasTrack,8);
  
  caloEm1->Draw(canvasTrack,9);
  caloEm2->Draw(canvasTrack,10);
  
  caloHad1->Draw(canvasTrack,11);
  caloHad2->Draw(canvasTrack,12);
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
                  vector<Events*> &eventsSignal, vector<Events*> &eventsBackground, vector<Events*> &eventsData)
{
  // sum up events from all background types
  int nBackgroundTotal=0;
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck] || !eventsBackground[iBck]) continue;
    nBackgroundTotal += eventsBackground[iBck]->weightedSize();
    
    if(printHeaders) cout<<backgroundTitle[iBck]<<"\t";
    cout<<eventsBackground[iBck]->weightedSize()<<endl;
  }
  
  // print S/B ratio for each signal type
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    if(printHeaders) cout<<signalTitle[iSig]<<"\tN events:\t";
    cout<<eventsSignal[iSig]->size()<<endl;
    
    if(printHeaders) cout<<signalTitle[iSig]<<"\tS/B:\t";
    cout<<eventsSignal[iSig]->size()/(double)nBackgroundTotal<<endl;
  }
}

int main(int argc, char* argv[])
{
  TApplication *theApp = new TApplication("App", &argc, argv);
  
  // All events with initial cuts only
  vector<Events*> eventsSignal, eventsBackground, eventsData;
  LoadEventsFromFiles(eventsSignal, eventsBackground, eventsData);
  
  //---------------------------------------------------------------------------
  // Define event, track and jet cuts
  //---------------------------------------------------------------------------
  
  EventCut  *initialEventCut = new EventCut();
  TrackCut  *initialTrackCut = new TrackCut();
  JetCut    *initialJetCut   = new JetCut();
  
  // obvious event cuts that we are sure we want to apply
  initialEventCut->SetNtracks(1, 999999);
  initialEventCut->SetMinNjets(1);
  initialEventCut->SetMaxNmuons(0);
  initialEventCut->SetMaxNtau(0);
  initialEventCut->SetMaxNlepton(0);
  
  // obvious track cuts
  initialTrackCut->SetRequireSameNpixelHitsLayers(true);
  initialTrackCut->SetNmissingInnerPixel(0, 0);
  initialTrackCut->SetNmissingMiddleTracker(0, 0);
  initialTrackCut->SetNpixelHits(0, 999999);
  
  // obvious jet cuts
  initialJetCut->SetMaxEta(2.4);
  initialJetCut->SetPtRange(80, 999999);
  
  
  initialEventCut->SetMinMetPt(200);
  initialEventCut->SetMinJetMetPhi(0.5);
  //  initialEventCut->SetNtracks(2, 2);
  
  ApplyCuts(eventsSignal, eventsBackground, eventsData,
            initialEventCut, initialTrackCut, initialJetCut, nullptr);
  
  DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
//  DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
  
  PrintSignalToBackground(eventsSignal, eventsBackground, eventsData);
  
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
      
      PrintSignalToBackground(currentEventsSignal, currentEventsBackground, currentEventsData);
    }
  }
  
  if(showPlots){
    if(interactive){
      MakeComparison(eventsSignal, eventsBackground, eventsData,
                     currentEventsSignal, currentEventsBackground, currentEventsData);
    }
    theApp->Run();
  }
    
  return 0;
}



