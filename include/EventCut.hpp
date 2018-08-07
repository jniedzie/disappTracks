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
    kEmpty          = 1,
    kOneTrack       = 1 << 1, ///< at least one track
    kOneJet         = 1 << 2, ///< at least one jet
    kMet100GeV      = 1 << 3, ///< MET pt ≥ 100 GeV
    kMetNoMu100GeV  = 1 << 4, ///< MET no mu pt ≥ 100 GeV
    kMetNoMuTrigger = 1 << 5, ///< require MET no mu trigger to fire
    kNoLepton       = 1 << 6, ///< requires no leptons in the event
    kNoTau          = 1 << 7, ///< requires no tau in the event
  };
  
  EventCut(ECut cutType=kEmpty);
  ~EventCut();
  
  // setters
  inline void SetMinMetPt(double min){minMetPt = min;}
  inline void SetMinMetNoMuPt(double min){minMetNoMuPt = min;}
  inline void SetMinNjets(int min){minNjets = min;}
  inline void SetMinNtracks(int min){minNtracks = min;}
  inline void SetMaxNlepton(int max){maxNleptons = max;}
  inline void SetMaxNtau(int max){maxNtau = max;}
  inline void SetRequireMetNoMuTrigger(bool val){metNoMuTrigger = val;}
  
  // getters
  inline double GetMinMetPt(){return minMetPt;}
  inline double GetMinMetNoMuPt(){return minMetNoMuPt;}
  inline int    GetMinNjets(){return minNjets;}
  inline int    GetMinNtracks(){return minNtracks;}
  inline int    GetMaxNlepton(){return maxNleptons;}
  inline int    GetMaxNtau(){return maxNtau;}
  inline bool   RequiresMetNoMuTrigger(){return metNoMuTrigger;}
  
private:
  double minMetPt;      ///< min MET pT
  double minMetNoMuPt;  ///< min MET no mu pT
  int minNjets;         ///< min number of jets
  int minNtracks;       ///< min number of tracks
  int maxNleptons;      ///< max number of leptons
  int maxNtau;         ///< max number of tau
  bool metNoMuTrigger;  ///< should require MET no mu trigger
};

#endif /* EventCut_hpp */
