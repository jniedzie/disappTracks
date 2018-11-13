#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"

int main()
{
  // All events with initial cuts only
  shared_ptr<EventSet> events = shared_ptr<EventSet>(new EventSet());
  events->LoadEventsFromFiles("after_L1/");
  
  EventCut  *eventCut = new EventCut();
  TrackCut  *trackCut = new TrackCut();
  JetCut    *jetCut   = new JetCut();
  
  eventCut->SetJetMetDeltaPhi(range<double>(0.5,inf));
  
  trackCut->SetNpixelLayers(range<int>(4,4));
  eventCut->SetNtracks(range<int>(1,1));
  
  // + standard cuts to be applied after L2 selections
  eventCut->SetNjets(range<int>(1,inf));
  eventCut->SetLeadingJetPt(range<double>(100,inf));
  eventCut->SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCut->SetLeadingJetNeHEF(range<double>(-inf,0.8));
  eventCut->SetLeadingJetChHEF(range<double>(0.1,inf));
  
  double dedxMin=0.0, dedxMax=5.1, dedxStep=0.1;
  double caloEmMin=0.0, caloEmMax=20.0,caloEmStep=1.0;
  double caloHadMin=0.0, caloHadMax=15.1,caloHadStep=2.0;
//  double metPtMin=200, metPtMax=401, metPtStep=25;
//  double jetMetPhiMin=0.1, jetMetPhiMax=1.4, jetMetPhiStep=0.1;
  double minNmissingMin=0, minNmissingMax=10, minNmissingStep=1;
//  double nPixelHitsMin=0, nPixelHitsMax=8, nPixelHitsStep=1;
  
  double bestSb[kNsignals] = {0};
  
  shared_ptr<EventSet> eventsAfterCuts;
  
  //     S/B,   dedx,   em cal, had cal, n_miss
  tuple<double, double, double, double, int> bestResults[kNsignals];
  
//  double sb_sum_best = 0;
  
  for(double dedxCut=dedxMin;dedxCut<dedxMax; dedxCut += dedxStep){
    cout<<'\r'<< setw(2) << setfill('0') << (dedxCut-dedxMin)/(dedxMax-dedxMin)<<"%"<<flush;
    
    for(double caloEmCut=caloEmMax;caloEmCut>caloEmMin; caloEmCut -= caloEmStep){
      eventsAfterCuts = shared_ptr<EventSet>(new EventSet(*events));
      
//      for(double caloHadCut=caloHadMax;caloHadCut>caloHadMin; caloHadCut -= caloHadStep){
        for(double minNmissingCut=minNmissingMin;minNmissingCut<minNmissingMax; minNmissingCut += minNmissingStep){
          cout<<"dE/dx:"<<dedxCut<<"\tcalo EM:"<<caloEmCut<<"\tn miss:"<<minNmissingCut<<endl;
//        for(double metPtCut=metPtMin;metPtCut<metPtMax; metPtCut += metPtStep){
//            for(double nPixelHitsCut=nPixelHitsMin;nPixelHitsCut<nPixelHitsMax; nPixelHitsCut += nPixelHitsStep){
              
              trackCut->SetDedxPerCluster(range<double>(dedxCut,inf));
              trackCut->SetCaloEmEnergy(range<double>(0.0,caloEmCut));
//              trackCut->SetCaloHadEnergy(range<double>(0.0,caloHadCut));
              trackCut->SetNmissingOuterTracker(range<int>(minNmissingCut, inf));
          
//              eventCut->SetMetPt(range<double>(metPtCut,inf));
//              trackCut->SetNpixelHits(range<int>(nPixelHitsCut,nPixelHitsCut));
              eventsAfterCuts->ApplyCuts(eventCut, trackCut, jetCut, nullptr);
              
              double nBackgroundTotal=0;
              for(int iBck=0;iBck<kNbackgrounds;iBck++){
                if(!runBackground[iBck]) continue;
                nBackgroundTotal += eventsAfterCuts->weightedSize(EventSet::kBackground, iBck);
              }
              double sb_sum=0;
              
              for(int iSig=0;iSig<kNsignals;iSig++){
                if(!runSignal[iSig]) continue;
                double sb = eventsAfterCuts->weightedSize(EventSet::kSignal,iSig)/sqrt(nBackgroundTotal+eventsAfterCuts->weightedSize(EventSet::kSignal,iSig));
              
                sb_sum += sb*sb;
                
                if(sb > bestSb[iSig]){
                  bestSb[iSig] = sb;
                  bestResults[iSig] = make_tuple(sb,dedxCut,caloEmCut,inf,minNmissingCut);
                }
              }
            }
          }
//        }
//      }
//    }
  }
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    
    cout<<"Best result for "<<signalTitle[iSig]<<":"<<endl;
    cout<<"S/sqrt(S+B):"<<get<0>(bestResults[iSig])<<endl;
    cout<<"min dedx:"<<get<1>(bestResults[iSig])<<endl;
    cout<<"max em calo:"<<get<2>(bestResults[iSig])<<endl;
    cout<<"max had calo:"<<get<3>(bestResults[iSig])<<endl;
    //                  cout<<"min MET pt:"<<metPtCut<<endl;
    //                  cout<<"min jet,met phi:"<<jetMetPhiCut<<endl;
    cout<<"min N missing hits:"<<get<4>(bestResults[iSig])<<endl;
    //                  cout<<"N pixel hits:"<<nPixelHitsCut<<endl;
  }
  
  cout<<"\ndone\n"<<endl;
  
  return 0;
}



