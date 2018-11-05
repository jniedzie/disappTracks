#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"

#include "TGraph.h"

#include <TApplication.h>

ESignal testSignal = kWino_M_300_cTau_30;

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

void ApplyCuts(vector<shared_ptr<Events>> &eventsSignal,
               vector<shared_ptr<Events>> &eventsBackground,
               vector<shared_ptr<Events>> &eventsData,
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

int main(int argc, char* argv[])
{
  TApplication *theApp = new TApplication("App", &argc, argv);
  
  // All events with initial cuts only
  vector<Events*> eventsSignal, eventsBackground, eventsData;
  LoadEventsFromFiles(eventsSignal, eventsBackground, eventsData, "after_L1/");
  
  EventCut  *eventCut = new EventCut();
  TrackCut  *trackCut = new TrackCut();
  JetCut    *jetCut   = new JetCut();
  
  // + standard cuts to be applied after L2 selections
  eventCut->SetNtracks(1, 999999);
  eventCut->SetMinNjets(1);
  eventCut->SetHighJetMinPt(100);
  eventCut->SetHighJetMaxEta(2.4);
  eventCut->SetHighJetMaxNeHEF(0.8);
  eventCut->SetHighJetMinChHEF(0.1);
  
  TGraph *sb[kNsignals];
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    sb[iSig] = new TGraph();
  }
  
  int iPoint=0;
  double cutMin=0, cutMax=20, cutStep=1;
  
  shared_ptr<Events> eventsAfterCuts[kNsignals];
  shared_ptr<Events> backAfterCuts[kNbackgrounds];
  
  for(double cut=cutMin;cut<cutMax; cut += cutStep){
    cout<<"cut:"<<cut<<endl;
    
//    trackCut->SetMinDedxPerCluster(cut);
//    eventCut->SetMinMetPt(cut);
//    eventCut->SetMinJetMetPhi(cut);
//    trackCut->SetMaxEmCalo(cut);
//    trackCut->SetMaxHadCalo(cut);
//    trackCut->SetNmissingOuterTracker(cut,999999);
    
    for(int iSig=0;iSig<kNsignals;iSig++){
      if(!runSignal[iSig]) continue;
      eventsAfterCuts[iSig] = eventsSignal[iSig]->ApplyCuts(eventCut, trackCut, jetCut, nullptr);
    }
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      if(!runBackground[iBck]) continue;
      backAfterCuts[iBck] = eventsBackground[iBck]->ApplyCuts(eventCut, trackCut, jetCut, nullptr);
    }
    
    double nBackgroundTotal=0;
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      if(!runBackground[iBck] || !backAfterCuts[iBck]) continue;
      nBackgroundTotal += backAfterCuts[iBck]->weightedSize();
    }
    nBackgroundTotal = sqrt(nBackgroundTotal);
    
    for(int iSig=0;iSig<kNsignals;iSig++){
      if(!runSignal[iSig]) continue;
      double val = eventsAfterCuts[iSig]->weightedSize()/nBackgroundTotal;
      sb[iSig]->SetPoint(iPoint,cut,val);
    }
    iPoint++;
  }
  cout<<"\ndone\n"<<endl;
  
  
  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  c2->cd();
  
  bool first=true;
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    
    sb[iSig]->SetMarkerStyle(signalMarkers[iSig]);
    sb[iSig]->SetMarkerSize(1.0);
    sb[iSig]->SetMarkerColor(SignalColor((ESignal)iSig));
    if(first){
      sb[iSig]->Draw("AP");
      sb[iSig]->GetXaxis()->SetTitle("cut value");
      sb[iSig]->GetYaxis()->SetTitle("S/#sqrt{B}");
      sb[iSig]->SetMinimum(0);
      first = false;
    }
    else{
      sb[iSig]->Draw("Psame");
    }
  }
  
  c2->Update();
  
  theApp->Run();
  
  return 0;
}



