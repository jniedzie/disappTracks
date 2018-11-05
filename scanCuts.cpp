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

int main(int argc, char* argv[])
{
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
  
  double dedxMin=3.0, dedxMax=6.1, dedxStep=0.1;
  double caloEmMin=0.1, caloEmMax=1.2,caloEmStep=0.5;
  double caloHadMin=0.1, caloHadMax=1.2,caloHadStep=0.5;
  double metPtMin=200, metPtMax=401, metPtStep=25;
  double jetMetPhiMin=0.1, jetMetPhiMax=1.5, jetMetPhiStep=0.1;
  double minNmissingMin=0, minNmissingMax=12, minNmissingStep=1;
  
  double bestSb[kNsignals];
  for(int iSig=0;iSig<kNsignals;iSig++){
    bestSb[iSig] = 0;
  }

  Events* eventsAfterCuts[kNsignals];
  Events* backAfterCuts[kNbackgrounds];
  
  for(double dedxCut=dedxMin;dedxCut<dedxMax; dedxCut += dedxStep){
//    for(double caloEmCut=caloEmMin;caloEmCut<caloEmMax; caloEmCut += caloEmStep){
//      for(double caloHadCut=caloHadMin;caloHadCut<caloHadMax; caloHadCut += caloHadStep){
//        for(double metPtCut=metPtMin;metPtCut<metPtMax; metPtCut += metPtStep){
          for(double jetMetPhiCut=jetMetPhiMin;jetMetPhiCut<jetMetPhiMax; jetMetPhiCut += jetMetPhiStep){
            for(double minNmissingCut=minNmissingMin;minNmissingCut<minNmissingMax; minNmissingCut += minNmissingStep){
              
              trackCut->SetMinDedxPerCluster(dedxCut);
//              trackCut->SetMaxEmCalo(caloEmCut);
//              trackCut->SetMaxHadCalo(caloHadCut);
//              eventCut->SetMinMetPt(metPtCut);
              eventCut->SetMinJetMetPhi(jetMetPhiCut);
              trackCut->SetNmissingOuterTracker(minNmissingCut, 999999);
              
              
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
                double sb = eventsAfterCuts[iSig]->weightedSize()/nBackgroundTotal;
              
                if(sb > bestSb[iSig]){
                  bestSb[iSig] = sb;
                  cout<<"S/B better than ever for signal "<<signalTitle[iSig]<<":"<<sb<<endl;
                  cout<<"min dedx:"<<dedxCut<<endl;
//                  cout<<"max em calo:"<<caloEmCut<<endl;
//                  cout<<"max had calo:"<<caloHadCut<<endl;
//                  cout<<"min MET pt:"<<metPtCut<<endl;
                  cout<<"min jet,met phi:"<<jetMetPhiCut<<endl;
                  cout<<"min N missing hits:"<<minNmissingCut<<endl;
                }
              }
                
              // clean up
              for(int iSig=0;iSig<kNsignals;iSig++){
                if(!runSignal[iSig]) continue;
                  if(eventsAfterCuts[iSig]){
                    delete eventsAfterCuts[iSig];
                    eventsAfterCuts[iSig] = nullptr;
                  }
              }
              for(int iBck=0;iBck<kNbackgrounds;iBck++){
                if(!runBackground[iBck]) continue;
                if(backAfterCuts[iBck]){
                  delete backAfterCuts[iBck];
                  backAfterCuts[iBck] = nullptr;
                }
              }
              
            }
          }
//        }
//      }
//    }
  }
  cout<<"\ndone\n"<<endl;
  
  return 0;
}



