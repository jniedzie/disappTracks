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
  shared_ptr<EventSet> events;
  events->LoadEventsFromFiles("after_L1/");
  
  auto eventCut = unique_ptr<EventCut>(new EventCut());
  auto trackCut = unique_ptr<TrackCut>(new TrackCut());
  auto jetCut   = unique_ptr<JetCut>(new JetCut());
  auto leptonCut= unique_ptr<LeptonCut>(new LeptonCut());
  
  // + standard cuts to be applied after L2 selections
  eventCut->SetNtracks(range<int>(1, inf));
  eventCut->SetNjets(range<int>(1,inf));
  eventCut->SetLeadingJetPt(range<double>(100,inf));
  eventCut->SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCut->SetLeadingJetNeHEF(range<double>(-inf,0.8));
  eventCut->SetLeadingJetChHEF(range<double>(0.1,inf));
  
  TGraph *sb[kNsignals];
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    sb[iSig] = new TGraph();
  }
  
  int iPoint=0;
  double cutMin=0, cutMax=20, cutStep=1;

  shared_ptr<EventSet> eventsAfterCuts;
  
  for(double cut=cutMin;cut<cutMax; cut += cutStep){
    cout<<"cut:"<<cut<<endl;
    
//    trackCut->SetMinDedxPerCluster(cut);
//    eventCut->SetMinMetPt(cut);
//    eventCut->SetMinJetMetPhi(cut);
//    trackCut->SetMaxEmCalo(cut);
//    trackCut->SetMaxHadCalo(cut);
//    trackCut->SetNmissingOuterTracker(cut,inf);
    
    eventsAfterCuts = events;
    events->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
    
    double nBackgroundTotal=0;
    for(int iBck=0;iBck<kNbackgrounds;iBck++){
      if(!runBackground[iBck]) continue;
      nBackgroundTotal += eventsAfterCuts->weightedSize(EventSet::kBackground,iBck);
    }
    nBackgroundTotal = sqrt(nBackgroundTotal);
    
    for(int iSig=0;iSig<kNsignals;iSig++){
      if(!runSignal[iSig]) continue;
      double val = eventsAfterCuts->weightedSize(EventSet::kSignal,iSig)/nBackgroundTotal;
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



