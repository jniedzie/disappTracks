//
//  Event.hpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#ifndef Event_hpp
#define Event_hpp

#include "Helpers.hpp"
#include "Track.hpp"
#include "TrackCut.hpp"
#include "Jet.hpp"

class Event;

class Events {
public:
  Events(std::string fileName);
  Events();
  ~Events();
  
  void AddEvent(Event *event){events.push_back(event);}
  Events* ApplyTrackCut(TrackCut *cut);
  
  inline unsigned long size(){return events.size();}
  Event* At(int index){return events[index];}
  int GetNtracks();
  
private:
  std::vector<Event*> events;
};

//---------------------------------------------------------------------------------------
// Single event class
//---------------------------------------------------------------------------------------

class Event{
public:
  Event(){};
  ~Event(){};
  
  void AddTrack(Track *track){tracks.push_back(track);}
  unsigned long GetNtracks(){return tracks.size(); }
  Track* GetTrack(int i){return tracks[i];}
  
  void Print();
//  
//  Event* FilterShortTracksAboveThreshold(double threshold);
//  Event* FilterShortTracks();
  
  Event* ApplyTrackCut(TrackCut *cut);
  
  vector<Track*> GetTracksPassingCut(TrackCut *cut);
  
private:
  vector<Track*> tracks;  // vector of isolated tracks
  vector<Jet*>   jets;    // vector of jets
  
};

#endif /* Event_hpp */
