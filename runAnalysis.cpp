#include "Event.hpp"
#include "EventSet.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"

#include "TGraph.h"

#include <TApplication.h>

void ProcessCuts(shared_ptr<EventSet> events,
                 const unique_ptr<EventCut> &eventCut,const  unique_ptr<TrackCut> &trackCut,
                 const unique_ptr<JetCut> &jetCut,const unique_ptr<LeptonCut> &leptonCut)
{
  events->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  if(printYields){
    cout<<"\n\nYields after level "<<performCutsLevel<<" cuts"<<endl;
    events->PrintYields();
  }
  if(saveEvents){
    string prefix = "after_L"+to_string(performCutsLevel)+"/";
    if(performCutsLevel==10) prefix = "adish_cuts";
    events->SaveEventsToFiles(prefix);
  }
  if(drawStandardPlots){
    events->DrawStandardPlots();
  }
  if(drawPerLayerPlots){
    events->DrawPerLayerPlots();
  }
  if(printBackgroundDetails){
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      if(!runBackground[iBck]) continue;
      cout<<"Background events in "<<backgroundTitle[iBck]<<":"<<endl;
      for(int iEvent=0;iEvent<events->size(EventSet::kBackground,iBck);iEvent++){
        events->At(EventSet::kBackground,iBck,iEvent)->Print();
      }
    }
  }
}

int main(int argc, char* argv[])
{
  TApplication *theApp = new TApplication("App", &argc, argv);
  
  // All events with initial cuts only
  shared_ptr<EventSet> events = shared_ptr<EventSet>(new EventSet);
  
  string initPrefix = "after_L"+to_string(performCutsLevel-1)+"/";
  if(performCutsLevel==0 || performCutsLevel==10) initPrefix = "";
  
  events->LoadEventsFromFiles(initPrefix);
  cout<<"\n\nInitial yields"<<endl;
  events->PrintYields();
  
  //---------------------------------------------------------------------------
  // Level 0
  //---------------------------------------------------------------------------

  if(performCutsLevel == 0){
    unique_ptr<EventCut>  eventCut_L0 = unique_ptr<EventCut>(new EventCut());
    unique_ptr<TrackCut>  trackCut_L0 = unique_ptr<TrackCut>(new TrackCut());
    unique_ptr<JetCut>    jetCut_L0   = unique_ptr<JetCut>(new JetCut());
    unique_ptr<LeptonCut> leptonCut_L0= unique_ptr<LeptonCut>(new LeptonCut());
    
    eventCut_L0->SetNtracks(range<int>(1, inf));
    eventCut_L0->SetNjets(range<int>(1,inf));
    eventCut_L0->SetNmuons(range<int>(0,0));
    eventCut_L0->SetNtaus(range<int>(0,0));
    eventCut_L0->SetNleptons(range<int>(0,0));

    eventCut_L0->SetMetNoMuPt(range<double>(200,inf));
    eventCut_L0->SetRequireMetNoMuTrigger(true);

    eventCut_L0->SetRequirePassingAllFilters(true);

    eventCut_L0->SetLeadingJetPt(range<double>(100,inf));
    eventCut_L0->SetLeadingJetEta(range<double>(-2.4,2.4));
    eventCut_L0->SetLeadingJetNeHEF(range<double>(-inf,0.8));
    eventCut_L0->SetLeadingJetChHEF(range<double>(0.1,inf));
    
    //  trackCut_L0->SetRequireSameNpixelHitsLayers(true);
    trackCut_L0->SetNmissingInnerPixel(range<int>(0, 0));
    trackCut_L0->SetNmissingMiddleTracker(range<int>(0, 0));
    trackCut_L0->SetNpixelLayers(range<int>(2, inf));
    trackCut_L0->SetEta(range<double>(-2.1, 2.1));

    jetCut_L0->SetPt(range<double>(30.0, inf));
    
    ProcessCuts(events,eventCut_L0, trackCut_L0, jetCut_L0, leptonCut_L0);
  }
  
  //---------------------------------------------------------------------------
  // Level 1
  //---------------------------------------------------------------------------
  
  if(performCutsLevel == 1){
    unique_ptr<EventCut>  eventCut_L1 = unique_ptr<EventCut>(new EventCut());
    unique_ptr<TrackCut>  trackCut_L1 = unique_ptr<TrackCut>(new TrackCut());
    unique_ptr<JetCut>    jetCut_L1   = unique_ptr<JetCut>(new JetCut());
    unique_ptr<LeptonCut> leptonCut_L1= unique_ptr<LeptonCut>(new LeptonCut());
    
    // L1 cuts
    trackCut_L1->SetRelativeIsolation(range<double>(0.0, 0.15));
    
    jetCut_L1->SetTrackDeltaR(range<double>(0.2,inf));
    jetCut_L1->SetChargedHadronEnergyFraction(range<double>(0.01,0.99));
    jetCut_L1->SetNeutralHadronEnergyFraction(range<double>(0.01,1.00));
    
    eventCut_L1->SetJetMetDeltaPhi(range<double>(0.5,inf));
    
    // + standard cuts to be applied after L2 selections
    eventCut_L1->SetNtracks(range<int>(1, inf));
    eventCut_L1->SetNjets(range<int>(1,inf));
    eventCut_L1->SetLeadingJetPt(range<double>(100,inf));
    eventCut_L1->SetLeadingJetEta(range<double>(-2.4,2.4));
    eventCut_L1->SetLeadingJetNeHEF(range<double>(-inf,0.8));
    eventCut_L1->SetLeadingJetChHEF(range<double>(0.1,inf));

    ProcessCuts(events,eventCut_L1, trackCut_L1, jetCut_L1, leptonCut_L1);
  }
    
  //---------------------------------------------------------------------------
  // Level 2
  //---------------------------------------------------------------------------
  if(performCutsLevel == 2){
    unique_ptr<EventCut>  eventCut_L2 = unique_ptr<EventCut>(new EventCut());
    unique_ptr<TrackCut>  trackCut_L2 = unique_ptr<TrackCut>(new TrackCut());
    unique_ptr<JetCut>    jetCut_L2   = unique_ptr<JetCut>(new JetCut());
    unique_ptr<LeptonCut> leptonCut_L2= unique_ptr<LeptonCut>(new LeptonCut());

    // pick category
//    trackCut_L2->SetNpixelLayers(range<int>(3, 3));
    eventCut_L2->SetNtracks(range<int>(2,2));
    
    // play with these cuts
    trackCut_L2->SetNmissingOuterTracker(range<int>(1, inf));
    trackCut_L2->SetCaloEmEnergy(range<double>(0.0,8.0));
//    trackCut_L2->SetCaloHadEnergy(range<double>(0.0,10.0));
//    trackCut_L2->SetDedxPerCluster(range<double>(2.0,inf));
//    eventCut_L2->SetMetPt(range<double>(230,inf));
//    trackCut_L2->SetPt(range<double>(100,inf));
//    trackCut_L2->SetTrackMetDeltaPhi(range<double>(-2.3,2.3));
//    eventCut_L2->SetJetMetDeltaPhi(range<double>(0.7,inf));
    
    
    // + standard cuts to be applied after L2 selections
    eventCut_L2->SetNjets(range<int>(1,inf));
    eventCut_L2->SetLeadingJetPt(range<double>(100,inf));
    eventCut_L2->SetLeadingJetEta(range<double>(-2.4,2.4));
    eventCut_L2->SetLeadingJetNeHEF(range<double>(-inf,0.8));
    eventCut_L2->SetLeadingJetChHEF(range<double>(0.1,inf));
    
    ProcessCuts(events, eventCut_L2, trackCut_L2, jetCut_L2, leptonCut_L2);
  }
  
  //---------------------------------------------------------------------------
  // Adish cuts
  //---------------------------------------------------------------------------
  if(performCutsLevel == 10){
    unique_ptr<EventCut>  eventCut_adish = unique_ptr<EventCut>(new EventCut());
    unique_ptr<TrackCut>  trackCut_adish = unique_ptr<TrackCut>(new TrackCut());
    unique_ptr<JetCut>    jetCut_adish = unique_ptr<JetCut>(new JetCut());
    unique_ptr<LeptonCut> leptonCut_adish= unique_ptr<LeptonCut>(new LeptonCut());
    
    // adish cuts
    eventCut_adish->SetRequireMetNoMuTrigger(true);
    eventCut_adish->SetMetNoMuPt(range<double>(200,inf));
    
    eventCut_adish->SetNjets(range<int>(1,inf));

    eventCut_adish->SetLeadingJetPt(range<double>(100,inf));
    eventCut_adish->SetLeadingJetEta(range<double>(-2.4,2.4));
    eventCut_adish->SetLeadingJetNeHEF(range<double>(-inf,0.8));
    eventCut_adish->SetLeadingJetChHEF(range<double>(0.1,inf));

    eventCut_adish->SetJetMetDeltaPhi(range<double>(0.5,inf));
    jetCut_adish->SetPt(range<double>(30, inf));
    
    eventCut_adish->SetNmuons(range<int>(0,0));
    eventCut_adish->SetNtaus(range<int>(0,0));
    eventCut_adish->SetNleptons(range<int>(0,0));
    
    ProcessCuts(events,eventCut_adish, trackCut_adish, jetCut_adish, leptonCut_adish);
  }
  
  if(drawStandardPlots || drawPerLayerPlots)  theApp->Run();
  return 0;
}



