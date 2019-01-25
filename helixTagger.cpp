//
//  helixTagger.cpp
//
//  Created by Jeremi Niedziela on 25/01/2019.
//

#include "Helpers.hpp"
#include "Fitter.hpp"
#include "FitterConfig.hpp"
#include "HelixProcessor.hpp"

#include "EventSet.hpp"

bool injectPion = true;
string configPath = "configs/helixFitter.md";

uint searchRun = 297100;
uint searchLumi = 136;
unsigned long long searchEvent = 245000232;

unique_ptr<HelixProcessor> helixProcessor;

void InjectPion(shared_ptr<vector<Point>> trackerPoints,
                shared_ptr<Track> track){
  double theta = track->GetTheta();
  double phi = track->GetPhi();
  
  int nHits = track->GetNtrackerLayers();
  double minL = layerR[nHits-1]/sin(theta);
  double maxL = layerR[nHits]/sin(theta);
  double decayR = RandDouble(minL, maxL);
  
  double x = decayR*sin(theta)*cos(phi);
  double y = decayR*sin(theta)*sin(phi);
  double z = decayR*cos(theta);
  
  track->SetDecayPoint(make_unique<Point>(x,y,z));
  unique_ptr<Helix> pionHelix = helixProcessor->GetRandomPionHelix(track);
  pionHelix->Print();
  trackerPoints->insert(trackerPoints->end(),pionHelix->GetPoints()->begin(), pionHelix->GetPoints()->end());
}

int main(int argc, char* argv[])
{
  auto config          = make_shared<FitterConfig>(configPath);
  auto pointsProcessor = make_unique<PointsProcessor>(config);
  auto fitter          = make_unique<Fitter>(config);
  helixProcessor       = make_unique<HelixProcessor>(config);
  
  auto events = make_shared<EventSet>();
  events->LoadEventsFromFiles("after_L2/3layers/");
  
  cout<<"helixTagger -- events loaded"<<endl;
  
  for(int iEvent=0; iEvent<events->size(EventSet::kData, kElectron_Run2017B); iEvent++){
    cout<<"\n\n=================================================================\n"<<endl;
    cout<<"helixTagger -- processing event "<<iEvent<<endl;
  
//    auto event = events->GetEvent(EventSet::kData, searchRun, searchLumi, searchEvent);
    auto event = events->At(EventSet::kData, kElectron_Run2017B, iEvent);
    
    shared_ptr<vector<Point>> trackerPoints = event->GetTrackerHits();
    
    if(!trackerPoints || trackerPoints->size()==0){
      cout<<"helixTagger -- no tracker hits for event "<<iEvent<<endl;
//      continue;
    }
    cout<<"helixTagger -- tracker points loaded"<<endl;
    
    for(int iTrack=0; iTrack<event->GetNtracks(); iTrack++){
      cout<<"helixTagger -- fitting helix for track:"<<iTrack<<endl;
      auto track = event->GetTrack(iTrack);
      
      if(injectPion) InjectPion(trackerPoints, track);
      
      unique_ptr<Helix> fittedHelix = fitter->GetBestFittingHelix(trackerPoints,track);
      if(fittedHelix){
        fittedHelix->Print();
        event->AddHelix(move(fittedHelix));
        cout<<"helixTagger -- fitted helix added to the event"<<endl;
      }
      else{
        cout<<"helixTagger -- could not fit a helix"<<endl;
      }
    }
  }
  
  cout<<"helixTagger -- saving events"<<endl;
  events->SaveEventsToFiles("afterHelixTagging/");
  cout<<"helixTagger -- finished"<<endl;
  return 0;
}
