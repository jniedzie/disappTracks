#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"

#include "TGraph.h"

#include <TApplication.h>

ESignal testSignal = kWino_M_300_cTau_30;

int main(int argc, char* argv[])
{
  TApplication *theApp = new TApplication("App", &argc, argv);
  
  // All events with initial cuts only
  vector<shared_ptr<Events>> eventsSignal, eventsBackground, eventsData;
  Events::LoadEventsFromFiles(eventsSignal, eventsBackground, eventsData, "after_L1/");
  
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



