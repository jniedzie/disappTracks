//
//  Event.hpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#ifndef Event_hpp
#define Event_hpp

#include "Helpers.hpp"
#include "EventCut.hpp"
#include "Track.hpp"
#include "TrackCut.hpp"
#include "Jet.hpp"
#include "JetCut.hpp"

class Event;

class Events {
public:
  Events(std::string fileName);
  Events();
  ~Events();
  
  void AddEvent(Event *event){events.push_back(event);}
  
  Events* ApplyCuts(EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut);
  Events* ApplyEventCut(EventCut *cut);
  Events* ApplyTrackCut(TrackCut *cut);
  Events* ApplyJetCut(JetCut *cut);
  
  inline int size(){return (int)events.size();}
  Event* At(int index){return events[index];}

private:
  Events(const Events &event);
  
  std::vector<Event*> events;
};

//---------------------------------------------------------------------------------------
// Single event class
//---------------------------------------------------------------------------------------

class Event{
public:
  Event();
  ~Event();
  
  // setters
  inline void AddTrack(Track *track){tracks.push_back(track);}
  inline void AddJet(Jet *jet){jets.push_back(jet);}
  
  inline void SetNvertices(int n){nVertices = n;}
  inline void SetNjet30(int n){nJet30 = n;}
  inline void SetNjet30a(int n){nJet30a = n;}
  inline void SetMetSumEt(double val){metSumEt = val;}
  inline void SetMetPt(double val){metPt = val;}
  inline void SetMetMass(double val){metMass = val;}
  inline void SetMetPhi(double val){metPhi = val;}
  inline void SetMetEta(double val){metEta = val;}
  
  // getters
  inline unsigned long GetNtracks(){return tracks.size(); }
  inline unsigned long GetNjets(){return jets.size(); }
  inline int GetNvertices(){return nVertices;}
  inline int GetNjet30(){return nJet30;}
  inline int GetNjet30a(){return nJet30a;}
  inline double GetMetSumEt(){return metSumEt;}
  inline double GetMetPt(){return metPt;}
  inline double GetMetMass(){return metMass;}
  inline double GetMetPhi(){return metPhi;}
  inline double GetMetEta(){return metEta;}
  
  inline Track*  GetTrack(int i){return tracks[i];}
  inline Jet*    GetJet(int i){return jets[i];}
  
  // other methods
  void Print();
  
  /// Returns a new event with only tracks passing the cut
  Event* ApplyTrackCut(TrackCut *cut);
  
  /// Returns a new event with only jets passing the cut
  Event* ApplyJetCut(JetCut *cut);
  
  /// Check if event passes a cut
  bool IsPassingCut(EventCut *cut);
private:
  vector<Track*> tracks;  // vector of isolated tracks
  vector<Jet*>   jets;    // vector of jets
  
  int nVertices;          // Number of good verices
  int nJet30;             // Number of jets with pt > 30, |eta|<2.4
  int nJet30a;            // Number of jets with pt > 30, |eta|<4.7
  double metSumEt;
  double metPt;
  double metMass;
  double metPhi;
  double metEta;
  
};

#endif /* Event_hpp */
