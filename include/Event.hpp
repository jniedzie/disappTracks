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

class Event{
public:
  Event(){};
  ~Event(){};
  
  void AddTrack(Track *track){tracks.push_back(track);}
  int GetNtracks(){return tracks.size(); }
  Track* GetTrack(int i){return tracks[i];}
  
  void Print();
  
  Event* FilterShortTracksAboveThreshold(double threshold);
  Event* FilterShortTracks();
  
  static map<unsigned long long,Event*> GetEventsFromFile(const char *fileName);
  static vector<Event*> GetEventsVectorFromFile(const char *fileName);

private:
  vector<Track*> tracks; // vector of isolated tracks
  
};

#endif /* Event_hpp */
