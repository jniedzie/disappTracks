#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"

#include "TGraph.h"

#include <TApplication.h>

void LoadEventsFromFiles(vector<Events*> &eventsSignal, vector<Events*> &eventsBackground, vector<Events*> &eventsData, string prefix="")
{
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]){
      eventsData.push_back(nullptr);
    }
    else{
      eventsData.push_back(new Events((inFileNameData[iData]+prefix+"tree.root"), Events::kData, maxNeventsData));
    }
  }
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]){
      eventsSignal.push_back(nullptr);
    }
    else{
      eventsSignal.push_back(new Events((inFileNameSignal[iSig]+prefix+"tree.root"), Events::kSignal, maxNeventsSignal,(ESignal)iSig));
    }
  }
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]){
      eventsBackground.push_back(nullptr);
    }
    else{
      eventsBackground.push_back(new Events());
      
      if(prefix==""){
        for(string path : inFileNameBackground[iBck]){
          eventsBackground[iBck]->AddEventsFromFile((path+prefix+"tree.root"),
                                                    Events::kBackground, maxNeventsBackground);
        }
      }
      else{
        string path = inFileNameBackground[iBck][0];
        eventsBackground[iBck]->AddEventsFromFile((path+prefix+"tree.root"),
                                                  Events::kBackground, maxNeventsBackground);
      }
    }
  }
}

void SaveEventsToFiles(vector<Events*> &eventsSignal, vector<Events*> &eventsBackground, vector<Events*> &eventsData, string prefix="after_L/")
{
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    system(("mkdir -p "+inFileNameSignal[iSig]+prefix).c_str());
    eventsSignal[iSig]->SaveToTree((inFileNameSignal[iSig]+prefix+"tree.root").c_str());
  }
  
  for(int iBck=0;iBck<kNbackgrounds;iBck++){
    if(!runBackground[iBck]) continue;
    
    // merged events will be stored in the first directory for given background
    string path = inFileNameBackground[iBck][0];
    system(("mkdir -p "+path+prefix).c_str());
    eventsBackground[iBck]->SaveToTree((path+prefix+"tree.root").c_str());
  }
  
  for(int iData=0;iData<kNdata;iData++){
    if(!runData[iData]) continue;
    system(("mkdir -p "+inFileNameData[iData]+prefix).c_str());
    eventsData[iData]->SaveToTree((inFileNameData[iData]+prefix+"tree.root").c_str());
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
  map<string, HistSet*> hists;
  
  hists["nVertices"]  = new HistSet(kNvertices);
  hists["nIsoTrack"]  = new HistSet(kNisoTracks);
  hists["nJet"]       = new HistSet(kNjets);
  hists["nJet30"]     = new HistSet(kNjets30);
  hists["nJet30a"]    = new HistSet(kNjets30a);
  hists["nMetSumEt"]  = new HistSet(kMetSumEt);
  hists["nMetPt"]     = new HistSet(kMetPt);
  hists["nMetMass"]   = new HistSet(kMetMass);
  hists["nMetEta"]    = new HistSet(kMetEta);
  hists["nMetPhi"]    = new HistSet(kMetPhi);
  hists["nMetJetDphi"]= new HistSet(kMetJetDphi);
  
  hists["nClustersPerTrack"]  = new HistSet(kTrackNclusters);
  hists["totalDeDx"]          = new HistSet(kTrackTotalDedx);
  hists["totalDeDxByNclusters"] = new HistSet(kTrackDedxPerCluster);
  hists["missingOuterTracker"]  = new HistSet(kTrackMissingOuterTrackerHits);
  
  hists["pt"]           = new HistSet(kTrackPt);
  hists["eta"]          = new HistSet(kTrackEta);
  hists["phi"]          = new HistSet(kTrackPhi);
  hists["caloEm"]       = new HistSet(kTrackCaloEm);
  hists["caloHad"]      = new HistSet(kTrackCaloHad);
  hists["pixelHits"]    = new HistSet(kTrackPixelHits);
  hists["trackerHits"]  = new HistSet(kTrackTrackerHits);
  hists["isolation"]    = new HistSet(kTrackRelativeIsolation);
  hists["trackMetDphi"] = new HistSet(kTrackMetDphi);
  hists["dedx"]         = new HistSet(kTrackDedxPerHit);
  
  hists["dxy"]    = new HistSet(kTrackDxy);
  hists["dz"]     = new HistSet(kTrackDz);
  hists["charge"] = new HistSet(kTrackCharge);
  hists["mass"]   = new HistSet(kTrackMass);
  hists["pid"]    = new HistSet(kTrackPid);
  
  hists["jet_pt"]     = new HistSet(kJetPt);
  hists["jet_eta"]    = new HistSet(kJetEta);
  hists["jet_phi"]    = new HistSet(kJetPhi);
  hists["jetTrackDr"] = new HistSet(kJetTrackDr);
  
  for(auto hist : hists){
    hist.second->FillFromEvents(eventsSignal, eventsBackground, eventsData);
  }
  
  // Plot histograms
  TCanvas *canvasEvents = new TCanvas((prefix+"Events").c_str(),(prefix+"Events").c_str(),2880,1800);
  canvasEvents->Divide(3,2);
  
  hists["nVertices"]->Draw(canvasEvents,1);
  hists["nIsoTrack"]->Draw(canvasEvents,2);
  hists["nJet"]->Draw(canvasEvents,3);
  hists["jet_pt"]->Draw(canvasEvents, 4);
  hists["nMetPt"]->Draw(canvasEvents,5);
  hists["nMetJetDphi"]->Draw(canvasEvents,6);
  
  TCanvas *canvasTrack = new TCanvas((prefix+"Tracks").c_str(),(prefix+"Tracks").c_str(),2880,1800);
  canvasTrack->Divide(4,3);
  
  hists["pt"]->Draw(canvasTrack,1);
  hists["caloEm"]->Draw(canvasTrack,2);
  hists["caloHad"]->Draw(canvasTrack,3);
  hists["missingOuterTracker"]->Draw(canvasTrack,4);
  hists["pixelHits"]->Draw(canvasTrack,5);
  hists["trackerHits"]->Draw(canvasTrack,6);
  hists["isolation"]->Draw(canvasTrack,7);
  hists["trackMetDphi"]->Draw(canvasTrack,8);
  hists["eta"]->Draw(canvasTrack,9);
  hists["dedx"]->Draw(canvasTrack,10);
  hists["dxy"]->Draw(canvasTrack,11);
  hists["dz"]->Draw(canvasTrack,12);
  
  TCanvas *canvasJets = new TCanvas("Jets","Jets",2880,1800);
  canvasJets->Divide(2,2);
  
  hists["jetTrackDr"]->Draw(canvasJets, 1);
  hists["jet_eta"]->Draw(canvasJets, 2);
  hists["jet_phi"]->Draw(canvasJets, 3);
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
    
    if(printYields){
      if(printHeaders) cout<<backgroundTitle[iBck]<<"\t";
      cout<<eventsBackground[iBck]->weightedSize()<<endl;
    }
  }

  if(printYields){
    for(int iSig=0;iSig<kNsignals;iSig++){
      if(!runSignal[iSig]) continue;
      if(printHeaders) cout<<signalTitle[iSig]<<"\tN events:\t";
      cout<<eventsSignal[iSig]->weightedSize()<<endl;
    }
  }
    
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    if(printHeaders) cout<<signalTitle[iSig]<<"\tS/sqrt(B):\t";
    cout<<eventsSignal[iSig]->weightedSize()/sqrt((double)nBackgroundTotal)<<endl;
  }
}

int main(int argc, char* argv[])
{
  TApplication *theApp = new TApplication("App", &argc, argv);
  
  // All events with initial cuts only
  vector<Events*> eventsSignal, eventsBackground, eventsData;
  
  string initPrefix;
  if(performCutsLevel==0) initPrefix = "";
  if(performCutsLevel==1) initPrefix = "after_L0/";
  if(performCutsLevel==2) initPrefix = "after_L1/";
  if(performCutsLevel==3) initPrefix = "after_L2/";
  if(performCutsLevel==10) initPrefix = "";
  
  LoadEventsFromFiles(eventsSignal, eventsBackground, eventsData, initPrefix);
  
  cout<<"\n\nInitial yields"<<endl;
  PrintYields(eventsSignal, eventsBackground, eventsData);
  
  //---------------------------------------------------------------------------
  // Level 0
  //---------------------------------------------------------------------------

  if(performCutsLevel == 0){
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
    
    SaveEventsToFiles(eventsSignal, eventsBackground, eventsData, "after_L0/");
    DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
    //  DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
  }
  
  //---------------------------------------------------------------------------
  // Level 1
  //---------------------------------------------------------------------------
  
  if(performCutsLevel == 1){
    EventCut  *eventCut_L1 = new EventCut();
    TrackCut  *trackCut_L1 = new TrackCut();
    JetCut    *jetCut_L1   = new JetCut();
    
    // L1 cuts
    trackCut_L1->SetMaxRelativeIsolation(0.15);
    jetCut_L1->SetMinTrackDeltaR(0.4);
    
    // + standard cuts to be applied after L2 selections
    eventCut_L1->SetNtracks(1, 9999999);
    eventCut_L1->SetMinNjets(1);
    eventCut_L1->SetHighJetMinPt(100);
    eventCut_L1->SetHighJetMaxEta(2.4);
    eventCut_L1->SetHighJetMaxNeHEF(0.8);
    eventCut_L1->SetHighJetMinChHEF(0.1);
    
    ApplyCuts(eventsSignal, eventsBackground, eventsData,eventCut_L1, trackCut_L1, jetCut_L1, nullptr);
    cout<<"\n\nYields after level 1 cuts"<<endl;
    PrintYields(eventsSignal, eventsBackground, eventsData);
    SaveEventsToFiles(eventsSignal, eventsBackground, eventsData, "after_L1/");
//    DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
//  DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
  }
    
  //---------------------------------------------------------------------------
  // Level 2
  //---------------------------------------------------------------------------
  if(performCutsLevel == 2){
    EventCut  *eventCut_L2 = new EventCut();
    TrackCut  *trackCut_L2 = new TrackCut();
    JetCut    *jetCut_L2   = new JetCut();
    
    // L2 cuts
    eventCut_L2->SetMinMetPt(230);
    eventCut_L2->SetMinJetMetPhi(0.5);

    trackCut_L2->SetMaxEmCalo(0.5);
    trackCut_L2->SetMaxHadCalo(0.5);
    trackCut_L2->SetNmissingOuterTracker(1, 999999);
    trackCut_L2->SetMinDedxPerCluster(2.2);
    
    // + standard cuts to be applied after L2 selections
    eventCut_L2->SetNtracks(1, 999999);
    eventCut_L2->SetMinNjets(1);
    eventCut_L2->SetHighJetMinPt(100);
    eventCut_L2->SetHighJetMaxEta(2.4);
    eventCut_L2->SetHighJetMaxNeHEF(0.8);
    eventCut_L2->SetHighJetMinChHEF(0.1);
    
    ApplyCuts(eventsSignal, eventsBackground, eventsData, eventCut_L2, trackCut_L2, jetCut_L2, nullptr);
    cout<<"\n\nYields after level 2 cuts"<<endl;
    PrintYields(eventsSignal, eventsBackground, eventsData);
    SaveEventsToFiles(eventsSignal, eventsBackground, eventsData, "after_L2/");
//    DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
//  DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
    
  }
  
  //---------------------------------------------------------------------------
  // Level 3
  //---------------------------------------------------------------------------
  if(performCutsLevel == 3){
    EventCut  *eventCut_L3 = new EventCut();
    TrackCut  *trackCut_L3 = new TrackCut();
    JetCut    *jetCut_L3   = new JetCut();
    
    // L2 cuts
//    eventCut_L3->SetMinMetPt(230);
//    eventCut_L3->SetMinJetMetPhi(0.5);
    
//    trackCut_L3->SetMaxEmCalo(0.5);
//    trackCut_L3->SetMaxHadCalo(0.5);
    trackCut_L3->SetNmissingOuterTracker(3, 999999);
    
    // + standard cuts to be applied after L2 selections
    eventCut_L3->SetNtracks(1, 999999);
    eventCut_L3->SetMinNjets(1);
    eventCut_L3->SetHighJetMinPt(100);
    eventCut_L3->SetHighJetMaxEta(2.4);
    eventCut_L3->SetHighJetMaxNeHEF(0.8);
    eventCut_L3->SetHighJetMinChHEF(0.1);
    
    ApplyCuts(eventsSignal, eventsBackground, eventsData, eventCut_L3, trackCut_L3, jetCut_L3, nullptr);
    cout<<"\n\nYields after level 3 cuts"<<endl;
    PrintYields(eventsSignal, eventsBackground, eventsData);
    SaveEventsToFiles(eventsSignal, eventsBackground, eventsData, "after_L3/");
    DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
    //  DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
  }
  
  //---------------------------------------------------------------------------
  // Adish cuts
  //---------------------------------------------------------------------------
  if(performCutsLevel == 10){
    EventCut  *eventCut_adish = new EventCut();
    TrackCut  *trackCut_adish = new TrackCut();
    JetCut    *jetCut_adish   = new JetCut();
    
    // adish cuts
    eventCut_adish->SetRequireMetNoMuTrigger(true);
    eventCut_adish->SetMinMetNoMuPt(200);
    
    eventCut_adish->SetMinNjets(1);

    eventCut_adish->SetHighJetMinPt(100);
    eventCut_adish->SetHighJetMaxEta(2.4);
    eventCut_adish->SetHighJetMaxNeHEF(0.8);
    eventCut_adish->SetHighJetMinChHEF(0.1);

    eventCut_adish->SetMinJetMetPhi(0.5);
    jetCut_adish->SetPtRange(30, 999999);
    
    eventCut_adish->SetMaxNmuons(0);
    eventCut_adish->SetMaxNtau(0);
    eventCut_adish->SetMaxNlepton(0);
    
    ApplyCuts(eventsSignal, eventsBackground, eventsData,eventCut_adish, trackCut_adish, jetCut_adish, nullptr);
    cout<<"\n\nYields after adish cuts"<<endl;
    PrintYields(eventsSignal, eventsBackground, eventsData);
//    SaveEventsToFiles(eventsSignal, eventsBackground, eventsData, "adish_cuts/");
//    DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
    //  DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
  }
  
    /*
  
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
    JetCut    *currentJetCut"]    = new JetCut();
    
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
*/
  theApp->Run();
  
  return 0;
}



