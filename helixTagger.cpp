//
//  helixTagger.cpp
//
//  Created by Jeremi Niedziela on 25/01/2019.
//

#include "Helpers.hpp"
#include "Fitter.hpp"
#include "ConfigManager.hpp"
#include "HelixProcessor.hpp"

#include "EventSet.hpp"

string configPath = "configs/helixTagger.md";

uint searchRun = 297100;
uint searchLumi = 136;
unsigned long long searchEvent = 245000232;

unique_ptr<HelixProcessor> helixProcessor;

void InjectPion(vector<shared_ptr<Point>> trackerPoints,
                shared_ptr<Track> track){
  double theta = track->GetTheta();
  double phi = track->GetPhi();

  int lastBarrelLayer = track->GetLastBarrelLayer();
  
  double minR = layerR[lastBarrelLayer-1];
  double maxR = layerR[lastBarrelLayer];
  double decayR = RandDouble(minR, maxR);
  
  double x = decayR*cos(phi);
  double y = decayR*sin(phi);
  double z = decayR/sin(theta)*cos(theta);
  
  track->SetDecayPoint(Point(x,y,z));
  unique_ptr<Helix> pionHelix = helixProcessor->GetRandomPionHelix(track);
  pionHelix->Print();
  trackerPoints.insert(trackerPoints.end(),pionHelix->GetPoints().begin(), pionHelix->GetPoints().end());
}

int main(int argc, char* argv[])
{
  config = ConfigManager(configPath);
  auto pointsProcessor = make_unique<PointsProcessor>();
  auto fitter          = make_unique<Fitter>();
  helixProcessor       = make_unique<HelixProcessor>();
  
  EventSet events;
  events.LoadEventsFromFiles("after_L1/");
  
  cout<<"helixTagger -- events loaded"<<endl;
  
  for(int iEvent=0; iEvent<events.size(xtracks::kSignal, kWino_M_300_cTau_10); iEvent++){
    cout<<"\n\n=================================================================\n"<<endl;
    cout<<"helixTagger -- processing event "<<iEvent<<endl;
  
//    auto event = events.GetEvent(xtracks::kData, searchRun, searchLumi, searchEvent);
    auto event = events.At(xtracks::kSignal, kWino_M_300_cTau_10, iEvent);
    
//    shared_ptr<vector<Point>> trackerPoints = event->GetTrackerHits();
    vector<shared_ptr<Point>> trackerPoints = pointsProcessor->GetRandomPoints(config.nNoiseHits);
    
    if(trackerPoints.size()==0){
      cout<<"helixTagger -- no tracker hits for event "<<iEvent<<endl;
//      continue;
    }
    cout<<"helixTagger -- tracker points loaded"<<endl;
    
    for(int iTrack=0; iTrack<event->GetNtracks(); iTrack++){
      cout<<"helixTagger -- fitting helix for track:"<<iTrack<<endl;
      auto track = event->GetTrack(iTrack);
      
      if(config.injectPionHits) InjectPion(trackerPoints, track);
      
      unique_ptr<Helix> fittedHelix = fitter->GetBestFittingHelix(trackerPoints, *track, *event->GetVertex());
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
  events.SaveEventsToFiles("afterHelixTagging/");
  cout<<"helixTagger -- finished"<<endl;
  return 0;
}
