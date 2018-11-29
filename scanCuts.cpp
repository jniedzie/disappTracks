#include "Event.hpp"
#include "TrackCut.hpp"
#include "JetCut.hpp"
#include "HistSet.hpp"
#include "Helpers.hpp"

struct ForRange {
  ForRange(){}
  ForRange(double _min, double _max, double _step) : min(_min), max(_max), step(_step) {}
  double min, max, step;
};

int main()
{
  // All events with initial cuts only
  shared_ptr<EventSet> events = shared_ptr<EventSet>(new EventSet());
  events->LoadEventsFromFiles("after_L1/");
  
  auto eventCut = unique_ptr<EventCut>(new EventCut());
  auto trackCut = unique_ptr<TrackCut>(new TrackCut());
  auto jetCut   = unique_ptr<JetCut>(new JetCut());
  auto leptonCut= unique_ptr<LeptonCut>(new LeptonCut());
  
  eventCut->SetJetMetDeltaPhi(range<double>(0.5,inf));
  
  trackCut->SetNpixelLayers(range<int>(4,4));
  eventCut->SetNtracks(range<int>(1,1));
  
  // + standard cuts to be applied after L2 selections
  eventCut->SetNjets(range<int>(1,inf));
  eventCut->SetLeadingJetPt(range<double>(100,inf));
  eventCut->SetLeadingJetEta(range<double>(-2.4,2.4));
  eventCut->SetLeadingJetNeHEF(range<double>(-inf,0.8));
  eventCut->SetLeadingJetChHEF(range<double>(0.1,inf));
  
  events->ApplyCuts(eventCut, trackCut, jetCut, leptonCut);
  
  map<string,ForRange> ranges;
  
  ranges["dedx"]          = ForRange(0.0,   6.0,    0.5);
  ranges["missingOuter"]  = ForRange(0,     13,     1);
  ranges["trackPt"]       = ForRange(50,    400,    25);
  ranges["trackMetPhi"]   = ForRange(1.9,   3.2,    0.1); // goes from max to min
  ranges["relIso"]        = ForRange(0.000, 0.151,  0.01);// goes from max to min
  
  double bestSb[kNsignals] = {0};
  
  shared_ptr<EventSet> eventsAfterCuts;
  map<string,double> bestResults[kNsignals];
  
//  double sb_sum_best = 0;
  
  for(double dedxCut=ranges["dedx"].min;dedxCut<ranges["dedx"].max; dedxCut += ranges["dedx"].step){
    for(double missingCut=ranges["missingOuter"].min;missingCut<ranges["missingOuter"].max; missingCut += ranges["missingOuter"].step){
      for(double trackPtCut=ranges["trackPt"].min;trackPtCut<ranges["trackPt"].max; trackPtCut += ranges["trackPt"].step){
        for(double trackMetPhiCut=ranges["trackMetPhi"].max;trackMetPhiCut>ranges["trackMetPhi"].min; trackMetPhiCut -= ranges["trackMetPhi"].step){
          
          // we can make a copy of the original events here, as in the next loop cuts will become tighter and tighter,
          // so events rejected in i-th iteration would be also rejected in (i+1)-th iteration
          eventsAfterCuts = shared_ptr<EventSet>(new EventSet(*events));
          
          for(double relIsoCut=ranges["relIso"].max;relIsoCut>ranges["relIso"].min; relIsoCut -= ranges["relIso"].step){
            cout<<"dE/dx:"<<dedxCut<<"\tmissing outer:"<<missingCut<<"\tn track pt:"<<trackPtCut<<"\treliso:"<<relIsoCut<<endl;
            
            trackCut->SetDedxPerCluster(range<double>(dedxCut,inf));
//            trackCut->SetCaloEmEnergy(range<double>(0.0,caloEmCut));
//            trackCut->SetCaloHadEnergy(range<double>(0.0,caloHadCut));
            trackCut->SetNmissingOuterTracker(range<int>(missingCut,inf));
            trackCut->SetPt(range<double>(trackPtCut,inf));
            trackCut->SetTrackMetDeltaPhi(range<double>(-trackMetPhiCut,trackMetPhiCut));
            trackCut->SetRelativeIsolation(range<double>(0,relIsoCut));
//            eventCut->SetMetPt(range<double>(230,inf));
            
            eventsAfterCuts->ApplyCuts(eventCut, trackCut, jetCut, nullptr);
            
            double nBackgroundTotal=0;
            for(int iBck=0;iBck<kNbackgrounds;iBck++){
              if(!runBackground[iBck]) continue;
              nBackgroundTotal += eventsAfterCuts->weightedSize(EventSet::kBackground, iBck);
            }

            for(int iSig=0;iSig<kNsignals;iSig++){
              if(!runSignal[iSig]) continue;
              double sb = eventsAfterCuts->weightedSize(EventSet::kSignal,iSig)/sqrt(nBackgroundTotal+eventsAfterCuts->weightedSize(EventSet::kSignal,iSig));
              
              if(sb > bestSb[iSig]){
                bestSb[iSig] = sb;
                bestResults[iSig]["dedx"] = dedxCut;
                bestResults[iSig]["missingOuter"] = missingCut;
                bestResults[iSig]["trackPt"] = trackPtCut;
                bestResults[iSig]["trackMetPhi"] = trackMetPhiCut;
                bestResults[iSig]["relIso"] = relIsoCut;
              }
            }
          }
        }
      }
    }
  }
  
  for(int iSig=0;iSig<kNsignals;iSig++){
    if(!runSignal[iSig]) continue;
    cout<<"Best result for "<<signalTitle[iSig]<<":"<<bestSb[iSig]<<endl;
    for(auto const& [title, val] : bestResults[iSig]){
      cout<<title<<":"<<val<<endl;
    }
  }
  
  return 0;
}



