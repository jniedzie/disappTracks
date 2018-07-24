//
//  EventCut.hpp
//
//  Created by Jeremi Niedziela on 23/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#ifndef EventCut_hpp
#define EventCut_hpp

class EventCut {
public:
  enum ECut {
    kEmpty,
    kOneTrack,                ///< require at least one track
    kOneJet,                  ///< require at least one jet
    kOneTrackOneJet,          ///< require at least one track and one jet
    kMet100GeV,               ///< require MET above 100 GeV
    kMet100GeVOneJet,         ///< require MET above 100 GeV and at least one jet
    kMet100GeVOneTrack,       ///< require MET above 100 GeV and at least one track
    kMet100GeVOneTrackOneJet, ///< require MET above 100 GeV, at least one track and at least one jet
  };
  
  EventCut(ECut cutType=kEmpty);
  ~EventCut();
  
  // setters
  inline void SetMinMetPt(double min){minMetPt = min;}
  inline void SetMinNjets(int min){minNjets = min;}
  inline void SetMinNtracks(int min){minNtracks = min;}
  
  // getters
  inline double GetMinMetPt(){return minMetPt;}
  inline int    GetMinNjets(){return minNjets;}
  inline int    GetMinNtracks(){return minNtracks;}
  
private:
  double minMetPt;      ///< min MET pT
  int minNjets;         ///< min number of jets
  int minNtracks;       ///< min number of tracks
  
};

#endif /* EventCut_hpp */
