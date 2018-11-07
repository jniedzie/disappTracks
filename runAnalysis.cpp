#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"

#include "TGraph.h"

#include <TApplication.h>

void ProcessCuts(vector<shared_ptr<Events>> eventsSignal,
                 vector<shared_ptr<Events>> eventsBackground,
                 vector<shared_ptr<Events>> eventsData,
                 EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut)
{
  Events::ApplyCuts(eventsSignal, eventsBackground, eventsData,eventCut, trackCut, jetCut, leptonCut);
  
  if(printYields){
    cout<<"\n\nYields after level "<<performCutsLevel<<" cuts"<<endl;
    Events::PrintYields(eventsSignal, eventsBackground, eventsData);
  }
  if(saveEvents){
    string prefix = "after_L"+to_string(performCutsLevel)+"/";
    if(performCutsLevel==10) prefix = "adish_cuts";
    Events::SaveEventsToFiles(eventsSignal, eventsBackground, eventsData, prefix);
  }
  if(drawStandardPlots){
    HistSet::DrawStandardPlots(eventsSignal, eventsBackground, eventsData);
  }
  if(drawPerLayerPlots){
    HistSet::DrawPerLayerPlots(eventsSignal, eventsBackground, eventsData);
  }
}

int main(int argc, char* argv[])
{
  TApplication *theApp = new TApplication("App", &argc, argv);
  
  // All events with initial cuts only
  vector<shared_ptr<Events>> eventsSignal, eventsBackground, eventsData;
  
  string initPrefix = "after_L"+to_string(performCutsLevel-1)+"/";
  if(performCutsLevel==0 || performCutsLevel==10) initPrefix = "";
  
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
    trackCut_L0->SetNmissingInnerPixel(range<int>(0, 0));
    trackCut_L0->SetNmissingMiddleTracker(range<int>(0, 0));
    trackCut_L0->SetNpixelLayers(range<int>(2, 999999));
    trackCut_L0->SetEta(range<double>(-2.1, 2.1));

    jetCut_L0->SetPt(range<double>(30.0, 999999));
    
    ProcessCuts(eventsSignal, eventsBackground, eventsData,eventCut_L0, trackCut_L0, jetCut_L0, nullptr);
  }
  
  //---------------------------------------------------------------------------
  // Level 1
  //---------------------------------------------------------------------------
  
  if(performCutsLevel == 1){
    EventCut  *eventCut_L1 = new EventCut();
    TrackCut  *trackCut_L1 = new TrackCut();
    JetCut    *jetCut_L1   = new JetCut();
    
    // L1 cuts
    trackCut_L1->SetRelativeIsolation(range<double>(0.0, 0.15));
    jetCut_L1->SetTrackDeltaR(range<double>(0.2,999999));
    
    // + standard cuts to be applied after L2 selections
    eventCut_L1->SetNtracks(1, 9999999);
    eventCut_L1->SetMinNjets(1);
    eventCut_L1->SetHighJetMinPt(100);
    eventCut_L1->SetHighJetMaxEta(2.4);
    eventCut_L1->SetHighJetMaxNeHEF(0.8);
    eventCut_L1->SetHighJetMinChHEF(0.1);

    ProcessCuts(eventsSignal, eventsBackground, eventsData,eventCut_L1, trackCut_L1, jetCut_L1, nullptr);
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
    trackCut_L2->SetNmissingOuterTracker(range<int>(8, 999999));
    trackCut_L2->SetDedxPerCluster(range<double>(3.5,999999));
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
    
    ProcessCuts(eventsSignal, eventsBackground, eventsData, eventCut_L2, trackCut_L2, jetCut_L2, nullptr);
    
    if(printBackgroundDetails){
      for(int iBck=0;iBck<kNbackgrounds;iBck++){
        if(!runBackground[iBck]) continue;
        cout<<"Background events in "<<backgroundTitle[iBck]<<":"<<endl;
        for(int iEvent=0;iEvent<eventsBackground[iBck]->size();iEvent++){
          eventsBackground[iBck]->At(iEvent)->Print();
        }
      }
    }
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
    trackCut_L3->SetNmissingOuterTracker(range<int>(3, 999999));
    
    // + standard cuts to be applied after L2 selections
    eventCut_L3->SetNtracks(1, 999999);
    eventCut_L3->SetMinNjets(1);
    eventCut_L3->SetHighJetMinPt(100);
    eventCut_L3->SetHighJetMaxEta(2.4);
    eventCut_L3->SetHighJetMaxNeHEF(0.8);
    eventCut_L3->SetHighJetMinChHEF(0.1);
    
    ProcessCuts(eventsSignal, eventsBackground, eventsData, eventCut_L3, trackCut_L3, jetCut_L3, nullptr);
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
    jetCut_adish->SetPt(range<double>(30, 999999));
    
    eventCut_adish->SetMaxNmuons(0);
    eventCut_adish->SetMaxNtau(0);
    eventCut_adish->SetMaxNlepton(0);
    
    ProcessCuts(eventsSignal, eventsBackground, eventsData,eventCut_adish, trackCut_adish, jetCut_adish, nullptr);
  }
  
  if(drawStandardPlots || drawPerLayerPlots)  theApp->Run();
  return 0;
}



