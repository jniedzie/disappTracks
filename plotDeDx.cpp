#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"

#include "TGraph.h"

#include <TApplication.h>

int main(int argc, char* argv[])
{
  TApplication *theApp = new TApplication("App", &argc, argv);
  
  // All events with initial cuts only
  vector<shared_ptr<Events>> eventsSignal, eventsBackground, eventsData;
  
  string initPrefix;
  if(performCutsLevel==0) initPrefix = "";
  if(performCutsLevel==1) initPrefix = "after_L0/";
  if(performCutsLevel==2) initPrefix = "after_L1/";
  if(performCutsLevel==3) initPrefix = "after_L2/";
  if(performCutsLevel==10) initPrefix = "";
  
  Events::LoadEventsFromFiles(eventsSignal, eventsBackground, eventsData, initPrefix);
  
  cout<<"\n\nInitial yields"<<endl;
  Events::PrintYields(eventsSignal, eventsBackground, eventsData);
  
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
    
    Events::ApplyCuts(eventsSignal, eventsBackground, eventsData,eventCut_L0, trackCut_L0, jetCut_L0, nullptr);
    cout<<"\n\nYields after level 0 cuts"<<endl;
    Events::PrintYields(eventsSignal, eventsBackground, eventsData);
    
    Events::SaveEventsToFiles(eventsSignal, eventsBackground, eventsData, "after_L0/");
    HistSet::DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
    //  HistSet::DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
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
    jetCut_L1->SetMinTrackDeltaR(0.2);
    
    // + standard cuts to be applied after L2 selections
    eventCut_L1->SetNtracks(1, 9999999);
    eventCut_L1->SetMinNjets(1);
    eventCut_L1->SetHighJetMinPt(100);
    eventCut_L1->SetHighJetMaxEta(2.4);
    eventCut_L1->SetHighJetMaxNeHEF(0.8);
    eventCut_L1->SetHighJetMinChHEF(0.1);

    Events::ApplyCuts(eventsSignal, eventsBackground, eventsData,eventCut_L1, trackCut_L1, jetCut_L1, nullptr);
    cout<<"\n\nYields after level 1 cuts"<<endl;
    Events::PrintYields(eventsSignal, eventsBackground, eventsData);
    Events::SaveEventsToFiles(eventsSignal, eventsBackground, eventsData, "after_L1/");
    HistSet::DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
//  HistSet::DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
  }
    
  //---------------------------------------------------------------------------
  // Level 2
  //---------------------------------------------------------------------------
  if(performCutsLevel == 2){
    EventCut  *eventCut_L2 = new EventCut();
    TrackCut  *trackCut_L2 = new TrackCut();
    JetCut    *jetCut_L2   = new JetCut();

    // pick category
//    trackCut_L2->SetNpixelLayers(4, 4);
    eventCut_L2->SetNtracks(2,2);
    
    // play with these cuts
    trackCut_L2->SetNmissingOuterTracker(8, 999999);
    trackCut_L2->SetMinDedxPerCluster(3.5);
//    eventCut_L2->SetMinMetPt(250);
//    trackCut_L2->SetMaxEmCalo(0.1);
//    trackCut_L2->SetMaxHadCalo(0.1);
    
    // cuts not to be optimized
    eventCut_L2->SetMinJetMetPhi(0.5);
    
    // + standard cuts to be applied after L2 selections
    eventCut_L2->SetMinNjets(1);
    eventCut_L2->SetHighJetMinPt(100);
    eventCut_L2->SetHighJetMaxEta(2.4);
    eventCut_L2->SetHighJetMaxNeHEF(0.8);
    eventCut_L2->SetHighJetMinChHEF(0.1);
    
    Events::ApplyCuts(eventsSignal, eventsBackground, eventsData, eventCut_L2, trackCut_L2, jetCut_L2, nullptr);
    cout<<"\n\nYields after level 2 cuts"<<endl;
    Events::PrintYields(eventsSignal, eventsBackground, eventsData);
    
//    for(int iBck=0;iBck<kNbackgrounds;iBck++){
//      if(!runBackground[iBck]) continue;
//      cout<<"Background events in "<<backgroundTitle[iBck]<<":"<<endl;
//      for(int iEvent=0;iEvent<eventsBackground[iBck]->size();iEvent++){
//        eventsBackground[iBck]->At(iEvent)->Print();
//      }
//    }
    
    //    Events::SaveEventsToFiles(eventsSignal, eventsBackground, eventsData, "after_L2/");
//    HistSet::DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
//  HistSet::DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
    
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
    
    Events::ApplyCuts(eventsSignal, eventsBackground, eventsData, eventCut_L3, trackCut_L3, jetCut_L3, nullptr);
    cout<<"\n\nYields after level 3 cuts"<<endl;
    Events::PrintYields(eventsSignal, eventsBackground, eventsData);
    Events::SaveEventsToFiles(eventsSignal, eventsBackground, eventsData, "after_L3/");
    HistSet::DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
    //  HistSet::DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
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
    
    Events::ApplyCuts(eventsSignal, eventsBackground, eventsData,eventCut_adish, trackCut_adish, jetCut_adish, nullptr);
    cout<<"\n\nYields after adish cuts"<<endl;
    Events::PrintYields(eventsSignal, eventsBackground, eventsData);
    //    Events::SaveEventsToFiles(eventsSignal, eventsBackground, eventsData, "adish_cuts/");
//    HistSet::DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
    //  HistSet::DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
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
      
      Events::ApplyCuts(currentEventsSignal, currentEventsBackground,currentEventsData,
                currentEventCut, currentTrackCut, currentJetCut, nullptr);
      
      Events::PrintYields(currentEventsSignal, currentEventsBackground, currentEventsData);
    }
  }
*/
  theApp->Run();
  
  return 0;
}



