#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"

int main(int argc, char* argv[])
{
  // All events with initial cuts only
  vector<shared_ptr<Events>> eventsSignal, eventsBackground, eventsData;
  Events::LoadEventsFromFiles(eventsSignal, eventsBackground, eventsData, "after_L1/");
  
  EventCut  *eventCut = new EventCut();
  TrackCut  *trackCut = new TrackCut();
  JetCut    *jetCut   = new JetCut();
  
  // + standard cuts to be applied after L2 selections
  eventCut->SetNtracks(range<int>(2,2));
  eventCut->SetNjets(range<int>(1,999999));
  eventCut->SetLeadingJetPt(range<double>(100,999999));
  eventCut->SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCut->SetLeadingJetNeHEF(range<double>(-999999,0.8));
  eventCut->SetLeadingJetChHEF(range<double>(0.1,999999));
  
  double dedxMin=2.0, dedxMax=4.3, dedxStep=0.1;
  double caloEmMin=0.1, caloEmMax=1.2,caloEmStep=0.5;
  double caloHadMin=0.1, caloHadMax=1.2,caloHadStep=0.5;
//  double metPtMin=200, metPtMax=401, metPtStep=25;
//  double jetMetPhiMin=0.1, jetMetPhiMax=1.4, jetMetPhiStep=0.1;
  double minNmissingMin=0, minNmissingMax=12, minNmissingStep=1;
  double nPixelHitsMin=0, nPixelHitsMax=8, nPixelHitsStep=1;
  
  double bestSb[kNsignals];
  for(int iSig=0;iSig<kNsignals;iSig++){
    bestSb[iSig] = 0;
  }

  shared_ptr<Events> eventsAfterCuts[kNsignals];
  shared_ptr<Events> backAfterCuts[kNbackgrounds];
  
  eventCut->SetJetMetDeltaPhi(range<double>(0.5,999999));
  
  double sb_sum_best = 0;
  
  for(double dedxCut=dedxMin;dedxCut<dedxMax; dedxCut += dedxStep){
//    for(double caloEmCut=caloEmMin;caloEmCut<caloEmMax; caloEmCut += caloEmStep){
//      for(double caloHadCut=caloHadMin;caloHadCut<caloHadMax; caloHadCut += caloHadStep){
//        for(double metPtCut=metPtMin;metPtCut<metPtMax; metPtCut += metPtStep){
//          for(double jetMetPhiCut=jetMetPhiMin;jetMetPhiCut<jetMetPhiMax; jetMetPhiCut += jetMetPhiStep){
            for(double minNmissingCut=minNmissingMin;minNmissingCut<minNmissingMax; minNmissingCut += minNmissingStep){
//            for(double nPixelHitsCut=nPixelHitsMin;nPixelHitsCut<nPixelHitsMax; nPixelHitsCut += nPixelHitsStep){
              
              trackCut->SetDedxPerCluster(range<double>(dedxCut,999999));
//              trackCut->SetMaxEmCalo(caloEmCut);
//              trackCut->SetMaxHadCalo(caloHadCut);
              //              eventCut->SetMinMetPt(metPtCut);
              trackCut->SetNmissingOuterTracker(range<int>(minNmissingCut, 999999));
//              trackCut->SetNpixelHits(nPixelHitsCut, nPixelHitsCut);
              
              
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
              
              double sb_sum=0;
              
              for(int iSig=0;iSig<kNsignals;iSig++){
                if(!runSignal[iSig]) continue;
                double sb = eventsAfterCuts[iSig]->weightedSize()/sqrt(nBackgroundTotal+eventsAfterCuts[iSig]->weightedSize());
              
                sb_sum += sb*sb;
                
//                if(sb > bestSb[iSig]){
//                  bestSb[iSig] = sb;
//                  cout<<"S/B better than ever for signal "<<signalTitle[iSig]<<":"<<sb<<endl;
//                  cout<<"min dedx:"<<dedxCut<<endl;
////                  cout<<"max em calo:"<<caloEmCut<<endl;
////                  cout<<"max had calo:"<<caloHadCut<<endl;
////                  cout<<"min MET pt:"<<metPtCut<<endl;
////                  cout<<"min jet,met phi:"<<jetMetPhiCut<<endl;
//                  cout<<"min N missing hits:"<<minNmissingCut<<endl;
////                  cout<<"N pixel hits:"<<nPixelHitsCut<<endl;
//                }
              }
              sb_sum = sqrt(sb_sum);
              
              if(sb_sum > sb_sum_best){
                sb_sum_best = sb_sum;
                cout<<"S/B better than ever:"<<sb_sum<<endl;
                cout<<"min dedx:"<<dedxCut<<endl;
                //                  cout<<"max em calo:"<<caloEmCut<<endl;
                //                  cout<<"max had calo:"<<caloHadCut<<endl;
                //                  cout<<"min MET pt:"<<metPtCut<<endl;
                //                  cout<<"min jet,met phi:"<<jetMetPhiCut<<endl;
                cout<<"min N missing hits:"<<minNmissingCut<<endl;
                //                  cout<<"N pixel hits:"<<nPixelHitsCut<<endl;
              }
              
            }
//          }
//        }
//      }
//    }
  }
  cout<<"\ndone\n"<<endl;
  
  return 0;
}



