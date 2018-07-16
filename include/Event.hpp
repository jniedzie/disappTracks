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

class Event;

class Events {
public:
  Events(std::string fileName);
  ~Events();
  
  inline unsigned long size(){return events.size();}
  Event* operator[] (const int index);
  
private:
  std::vector<Event*> events;
  std::map<unsigned long long,Event*> GetEventsFromFile(std::string fileName);
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
  
  Event* FilterShortTracksAboveThreshold(double threshold);
  Event* FilterShortTracks();
  
private:
  vector<Track*> tracks; // vector of isolated tracks
  
  
};

#endif /* Event_hpp */
