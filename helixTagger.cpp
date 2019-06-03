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
string cutLevel = "after_L1/3layers/";//after_L1/";

int nEvents = 50;

xtracks::EDataType dataType = xtracks::kSignal;
int setIter = kWino_M_300_cTau_10;

bool removePionClusters = false;

//vector<int> eventIndices = { 3, 7, 28, 30, 37, 43, };
vector<int> eventIndices = { 3,7,11,13,18,19,22,26,27,28,30,31,33,34,37,38,39,41,43,44, };

map<int, double> correctFractions;
map<int, double> fakeFractions;
map<int, double> nHelixPoints;

shared_ptr<Event> GetEvent(int iEvent);
vector<shared_ptr<Point>> GetClustersNoEndcaps(const shared_ptr<Event> &event, bool removePionClusters);

int main(int argc, char* argv[])
{
  config = ConfigManager(configPath);
  auto fitter = make_unique<Fitter>();
  
  cout<<"helixTagger -- events loaded"<<endl;
  
  int nAnalyzedEvents=0;
  
//  for(auto iEvent : eventIndices){
  for(auto iEvent=0; iEvent<nEvents; iEvent++){
    if(iEvent==30 || iEvent==156) continue;

    auto event = GetEvent(iEvent);
    
    cout<<"\n\n=================================================================\n"<<endl;
    cout<<"helixTagger -- processing event "<<iEvent<<endl;
    
    auto pointsNoEndcaps = GetClustersNoEndcaps(event, removePionClusters);
    
    if(pointsNoEndcaps.size()==0){
      cout<<"helixTagger -- no tracker hits for event "<<iEvent<<endl;
      //      continue;
    }
    cout<<"helixTagger -- tracker points loaded"<<endl;
    
    for(auto &track : event->GetTracks()){
      
      vector<Helix> fittedHelices = fitter->FitHelices(pointsNoEndcaps, *track, *event->GetVertex());
      
      double averageCorrectFraction=0;
      double averageFakeFraction=0;
      double averageHelixNpoints=0;
      
      for(auto helix : fittedHelices){
        helix.Print();
        
        auto helixPoints = helix.GetPoints();
        averageHelixNpoints += helixPoints.size();
        
        int nPionPoints = 0;
        for(auto &pionPoint : event->GetPionClusters()){
          for(auto &helixPoint : helixPoints){
            if(*pionPoint == *helixPoint) nPionPoints++;
          }
        }
        
        int nFakePoints = (int)helixPoints.size()-1 - nPionPoints;

        double fractionCorrectPoints = ((double)nPionPoints) / (helixPoints.size()-1);
        double fractionFakePoints    = ((double)nFakePoints) / (helixPoints.size()-1);
        
        averageCorrectFraction += fractionCorrectPoints;
        averageFakeFraction += fractionFakePoints;
        cout<<"Fraction of correctly assigned points:"<<fractionCorrectPoints<<endl;
        cout<<"Fraction of fake points:"<<fractionFakePoints<<endl;
        
        //      event->AddHelix(move(fittedHelix));
      }
      
      averageCorrectFraction  /= fittedHelices.size()>0 ? fittedHelices.size() : 1;
      averageFakeFraction     /= fittedHelices.size()>0 ? fittedHelices.size() : 1;
      averageHelixNpoints     /= fittedHelices.size()>0 ? fittedHelices.size() : 1;
      
      correctFractions[iEvent] = averageCorrectFraction;
      fakeFractions[iEvent]    = averageFakeFraction;
      nHelixPoints[iEvent]     = averageHelixNpoints;
    }
    
    nAnalyzedEvents++;
  }
  
//  for(auto &[iEvent, correctFraction] : correctFractions){
//    cout<<"Event: "<<iEvent<<"\tcorrect hits fraction:"<<correctFraction<<endl;
//  }
//
//  for(auto &[iEvent, fakeFraction] : fakeFractions){
//    cout<<"Event: "<<iEvent<<"\tfake hits fraction:"<<fakeFraction<<endl;
//  }
  vector<double> efficiency;
  for(int iThreshold=0; iThreshold<15; iThreshold++){efficiency.push_back(0);}
  
  for(auto &[iEvent, nPoints] : nHelixPoints){
//    cout<<"Event: "<<iEvent<<"\tn helix points:"<<nPoints<<endl;
    for(int iThreshold=0; iThreshold<15; iThreshold++){
      if(nPoints >= iThreshold) efficiency[iThreshold]+=1;
    }
  }

  for(int iThreshold=0; iThreshold<15; iThreshold++){
    efficiency[iThreshold] /= nAnalyzedEvents;
//    cout<<"Efficiency (N >= "<<iThreshold<<"): "<<efficiency[iThreshold]<<endl;
    cout<<efficiency[iThreshold]<<endl;
  }
  
  //  cout<<"helixTagger -- saving events"<<endl;
  //  events.SaveEventsToFiles("afterHelixTagging/");
  //  cout<<"helixTagger -- finished"<<endl;
  return 0;
}

shared_ptr<Event> GetEvent(int iEvent)
{
  EventSet events;
  events.LoadEventFromFiles(dataType, setIter, iEvent, cutLevel);
  auto event = events.At(dataType, setIter, 0);
  
  if(!event){
    cout<<"helixTagger -- event not found"<<endl;
    exit(0);
  }
  
  return event;
}


vector<shared_ptr<Point>> GetClustersNoEndcaps(const shared_ptr<Event> &event, bool removePionClusters)
{
  vector<shared_ptr<Point>> pointsNoEndcaps;
  
  for(auto &point : event->GetTrackerClusters()){
    if(point->GetSubDetName() == "TID" || point->GetSubDetName() == "TEC" || point->GetSubDetName() == "P1PXEC") continue;
    
    if(point->GetSubDetName() != "TIB" && point->GetSubDetName() != "TOB" && point->GetSubDetName() != "P1PXB"){
      cout<<"Weird detector:"<<point->GetSubDetName()<<endl;
    }
    
    if(removePionClusters){
      bool isPionHit = false;
      for(auto &pionCluster : event->GetPionClusters()){
        if(*pionCluster == *point){
          isPionHit = true;
          break;
        }
      }
      if(isPionHit) continue;
    }
    
    pointsNoEndcaps.push_back(point);
  }
  
  return pointsNoEndcaps;
}
