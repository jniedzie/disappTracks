//  EventCut.hpp
//
//  Created by Jeremi Niedziela on 23/07/2018.

#ifndef EventCut_hpp
#define EventCut_hpp

#include "Helpers.hpp"

/// Class containing definition of the event selection criteria. User can define ranges
/// of allowed parameters and required flags.
class EventCut {
public:
  /// Default constructor
  EventCut();
  
  /// Copy constructor
  EventCut(const EventCut &c);
  
  /// Default destructor
  ~EventCut();
  
  // Setters
  inline void SetNtracks(range<int> val){nTracks=val;}
  inline void SetNjets(range<int> val){nJets=val;}
  inline void SetNleptons(range<int> val){nLeptons=val;}
  inline void SetNmuons(range<int> val){nMuons=val;}
  inline void SetNtaus(range<int> val){nTaus=val;}
  
  inline void SetMetPt(range<double> val){metPt=val;}
  inline void SetMetNoMuPt(range<double> val){metNoMuPt=val;}
  inline void SetJetMetDeltaPhi(range<double> val){jetMetDeltaPhi=val;}
  inline void SetJetMuonDeltaPhi(range<double> val){jetMuonDeltaPhi=val;}
  inline void SetTrackMuonDeltaPhi(range<double> val){trackMuonDeltaPhi=val;}
  
  inline void SetLeadingJetPt(range<double> val){leadingJetPt=val;}
  inline void SetLeadingJetEta(range<double> val){leadingJetEta=val;}
  inline void SetLeadingJetChHEF(range<double> val){leadingJetChHEF=val;}
  inline void SetLeadingJetNeHEF(range<double> val){leadingJetNeHEF=val;}
  
  inline void SetRequireMetNoMuTrigger(bool val){requiresMetNoMuTrigger = val;}
  inline void SetRequireMuonsFromZ(bool val){requiresMuonsFromZ = val;}
  inline void SetRequireTightMuon(bool val){requiresTightMuon = val;}
  inline void SetRequireTwoOppositeMuons(bool val){requiresTwoOpositeMuons = val;}
  inline void SetRequirePassingAllFilters(bool val){requiresPassingAllFilters = val;}
  
private:
  range<int> nTracks;   ///< allowed number of tracks
  range<int> nJets;     ///< allowed number of jets
  range<int> nLeptons;  ///< allowed number of leptons
  range<int> nMuons;    ///< allowed number if muons
  range<int> nTaus;     ///< allowed number of taus
  
  range<double> metPt;              ///< allowed MET transverse momentum
  range<double> metNoMuPt;          ///< allowed MET no mu transverse momentum
  range<double> jetMetDeltaPhi;     ///< allowed angle between jet and MET
  range<double> jetMuonDeltaPhi;    ///< allowed angle between jet and muon
  range<double> trackMuonDeltaPhi;  ///< allowed angle between track and muon
  range<double> jetTrackDeltaR;     ///< allowed angle between track and muon
  
  range<double> leadingJetPt;     ///< allowed pt of the leading jet
  range<double> leadingJetEta;    ///< allowed pseudorapidity of the leading jet
  range<double> leadingJetChHEF;  ///< allowed charged hadron energy fraction of the leading jet
  range<double> leadingJetNeHEF;  ///< allowed neutral hadron energy fraction of the leading jet
  
  bool requiresMetNoMuTrigger;    ///< should require MET no mu trigger
  bool requiresMuonsFromZ;        ///< should require two muons in the event with invariant mass close to Z mass
  bool requiresTightMuon;         ///< should require at least one muon to pass tight id
  bool requiresTwoOpositeMuons;   ///< should require exactly two muons with opposite signs
  bool requiresPassingAllFilters; ///< should event pass all MET filters
  
  friend class EventProcessor;
};

#endif /* EventCut_hpp */
