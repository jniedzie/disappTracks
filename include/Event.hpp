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
#include "Lepton.hpp"
#include "LeptonCut.hpp"

class Event;

class Events {
public:
  Events(std::string fileName, int dataType=0); // 0 - background, 1 - signal, 2 - data
  Events();
  ~Events();
  
  void AddEvent(Event *event){events.push_back(event);}
  
  Events* ApplyCuts(EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut);
  Events* ApplyEventCut(EventCut *cut);
  Events* ApplyTrackCut(TrackCut *cut);
  Events* ApplyJetCut(JetCut *cut);
  Events* ApplyLeptonCut(LeptonCut *cut);
  
  inline int size(){return (int)events.size();}
  double WeightedSize();
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
  inline void AddLepton(Lepton *lepton){leptons.push_back(lepton);}
  
  inline void SetWeight(double val){weight = val;}
  
  inline void SetNvertices(int n){nVertices = n;}
  inline void SetNjet30(int n){nJet30 = n;}
  inline void SetNjet30a(int n){nJet30a = n;}
  inline void SetNlepton(int n){nLepton = n;}
  inline void SetNtau(int n){nTau = n;}
  
  inline void SetMetSumEt(double val){metSumEt = val;}
  inline void SetMetPt(double val){metPt = val;}
  inline void SetMetMass(double val){metMass = val;}
  inline void SetMetPhi(double val){metPhi = val;}
  inline void SetMetEta(double val){metEta = val;}
  
  inline void SetMetNoMuPt(double val){metNoMuPt = val;}
  inline void SetMetNoMuMass(double val){metNoMuMass = val;}
  inline void SetMetNoMuPhi(double val){metNoMuPhi = val;}
  inline void SetMetNoMuEta(double val){metNoMuEta = val;}
  inline void SetHasNoMuTrigger(bool val){metNoMuTrigger = val;}
  
  // getters
  inline double GetWeight(){return weight;}
  
  inline unsigned long GetNtracks(){return tracks.size(); }
  inline unsigned long GetNjets(){return jets.size(); }
  inline int GetNvertices(){return nVertices;}
  inline int GetNjet30(){return nJet30;}
  inline int GetNjet30a(){return nJet30a;}
  inline int GetNlepton(){return nLepton;}
  inline int GetNtau(){return nTau;}
  
  inline double GetMetSumEt(){return metSumEt;}
  inline double GetMetPt(){return metPt;}
  inline double GetMetMass(){return metMass;}
  inline double GetMetPhi(){return metPhi;}
  inline double GetMetEta(){return metEta;}
  
  inline double GetMetNoMuPt(){return metNoMuPt;}
  inline double GetMetNoMuMass(){return metNoMuMass;}
  inline double GetMetNoMuPhi(){return metNoMuPhi;}
  inline double GetMetNoMuEta(){return metNoMuEta;}
  inline bool HetMetNoMuTrigger(){return metNoMuTrigger;}
  
  inline Track*  GetTrack(int i){return tracks[i];}
  inline Jet*    GetJet(int i){return jets[i];}
  inline Lepton* GetLepton(int i){return leptons[i];}
  
  // other methods
  void Print();
  
  /// Returns a new event with only tracks passing the cut
  Event* ApplyTrackCut(TrackCut *cut);
  
  /// Returns a new event with only jets passing the cut
  Event* ApplyJetCut(JetCut *cut);
  
  /// Returns a new event with only leptons passing the cut
  Event* ApplyLeptonCut(LeptonCut *cut);
  
  /// Check if event passes a cut
  bool IsPassingCut(EventCut *cut);
  
  /// Returns an event with global params copied (but no colletions of tracks, jets, leptons
  Event* CopyThisEventProperties();
  
private:
  vector<Track*>  tracks;   // vector of isolated tracks
  vector<Jet*>    jets;     // vector of jets
  vector<Lepton*> leptons;  // vector of leptons
  
  double weight;          // Weight for this event resulting from lumi, xsec and number of events generated
  int nVertices;          // Number of good verices
  int nJet30;             // Number of jets with pt > 30, |eta|<2.4
  int nJet30a;            // Number of jets with pt > 30, |eta|<4.7
  int nLepton;
  int nTau;
  
  double metSumEt;
  double metPt;
  double metMass;
  double metPhi;
  double metEta;
  
  bool metNoMuTrigger;
  double metNoMuPt;
  double metNoMuMass;
  double metNoMuPhi;
  double metNoMuEta;
  
  
};

#endif /* Event_hpp */
