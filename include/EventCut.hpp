//
//  EventCut.hpp
//
//  Created by Jeremi Niedziela on 23/07/2018.
//

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
  
  inline void SetLeadingJetPt(range<double> val){leadingJetPt=val;}
  inline void SetLeadingJetEta(range<double> val){leadingJetEta=val;}
  inline void SetLeadingJetChHEF(range<double> val){leadingJetChHEF=val;}
  inline void SetLeadingJetNeHEF(range<double> val){leadingJetNeHEF=val;}
  
  inline void SetRequireMetNoMuTrigger(bool val){metNoMuTrigger = val;}
  inline void SetRequireMuonsFromZ(bool val){muonsFromZ = val;}
  inline void SetRequireMetJetPhi0p5(bool val){metJetPhi = val;}
  inline void SetRequireMetNoMuJetPhi0p5(bool val){metNoMuJetPhi = val;}
  inline void SetRequireMuJetR0p4(bool val){muJetR0p4 = val;}
  inline void SetRequireMuTrackR0p4(bool val){muTrackR0p4 = val;}
  inline void SetRequireHighJet(bool val){highJet = val;}
  inline void SetRequireTightMuon(bool val){tightMuon = val;}
  inline void SetRequireTwoOppositeMuons(bool val){twoOpositeMuons = val;}
  
  inline void SetRequirePassingAllFilters(bool val){requirePassAllFilters = val;}
  
  
  // Getters
  inline range<int>    GetNtracks(){return nTracks;}
  inline range<int>    GetNjets(){return nJets;}
  inline range<int>    GetNleptons(){return nLeptons;}
  inline range<int>    GetNmuons(){return nMuons;}
  inline range<int>    GetNtaus(){return nTaus;}

  inline range<double> GetMetPt(){return metPt;}
  inline range<double> GetMetNoMuPt(){return metNoMuPt;}
  inline range<double> GetJetMetDeltaPhi(){return jetMetDeltaPhi;}
  
  inline range<double> GetLeadingJetPt(){return leadingJetPt;}
  inline range<double> GetLeadingJetEta(){return leadingJetEta;}
  inline range<double> GetLeadingJetChHEF(){return leadingJetChHEF;}
  inline range<double> GetLeadingJetNeHEF(){return leadingJetNeHEF;}
  
  inline bool RequiresMetNoMuTrigger() const {return metNoMuTrigger;}
  inline bool RequiresMuonsFromZ() const {return muonsFromZ;}
  inline bool RequiresMetJetPhi0p5() const {return metJetPhi;}
  inline bool RequiresMetNoMuJetPhi0p5() const {return metNoMuJetPhi;}
  inline bool RequiresMuJetR0p4() const {return muJetR0p4;}
  inline bool RequiresMuTrackR0p4() const {return muTrackR0p4;}
  inline bool RequiresHighJet() const {return highJet;}
  inline bool RequiresTightMuon() const {return tightMuon;}
  inline bool RequiresTwoOppositeMuons() const {return twoOpositeMuons;}
 
  inline bool GetRequiresPassingAllFilters(){return requirePassAllFilters;}
  
private:
  range<int> nTracks;         ///< allowed number of tracks
  range<int> nJets;           ///< allowed number of jets
  range<int> nLeptons;        ///< allowed number of leptons
  range<int> nMuons;          ///< allowed number if muons
  range<int> nTaus;            ///< allowed number of taus
  
  range<double> metPt;           ///< allowed MET transverse momentum
  range<double> metNoMuPt;       ///< allowed MET no mu transverse momentum
  range<double> jetMetDeltaPhi;  ///< allowed angle between jet and MET
  
  range<double> leadingJetPt;    ///< allowed pt of the leading jet
  range<double> leadingJetEta;   ///< allowed pseudorapidity of the leading jet
  range<double> leadingJetChHEF; ///< allowed charged hadron energy fraction of the leading jet
  range<double> leadingJetNeHEF; ///< allowed neutral hadron energy fraction of the leading jet
  
  bool metNoMuTrigger;  ///< should require MET no mu trigger
  bool muonsFromZ;      ///< should require two muons in the event with invariant mass close to Z mass
  bool metJetPhi;       ///< should require Δφ(MET pt,jet pt) ≥ 0.5 for each jet
  bool metNoMuJetPhi;   ///< should require Δφ(MET pt,jet pt) ≥ 0.5 for each jet
  bool muJetR0p4;       ///< should require ΔR(mu,jet) ≥ 0.4 for each jet and each muon
  bool muTrackR0p4;     ///< should require ΔR(mu,track) ≥ 0.4 for each track and each muon
  bool highJet;         ///< should require at least one jet above 100 GeV
  bool tightMuon;       ///< should require at least one muon to pass tight id
  bool twoOpositeMuons; ///< should require exactly two muons with opposite signs
  
  bool requirePassAllFilters; ///< should event pass all filters
  
  friend class EventProcessor;
};

#endif /* EventCut_hpp */
