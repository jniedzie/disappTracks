//
//  EventCut.hpp
//
//  Created by Jeremi Niedziela on 23/07/2018.
//

#ifndef EventCut_hpp
#define EventCut_hpp

/// Keeps representation of cuts to be applied to events
class EventCut {
public:
  /// List of possible event cuts. Can be binary combined
  /// (e.g. kOneTrack | kOne Jet means require both one track and one jet).
  enum ECut {
    kEmpty            = 1,        ///< if this option is specified, no other options will be set!
    kOneTrack         = 1 << 1,   ///< at least one track
    kOneJet           = 1 << 2,   ///< at least one jet
    kOneLepton        = 1 << 3,   ///< at least one lepton
    kNoLepton         = 1 << 4,   ///< requires no leptons in the event
    kNoTau            = 1 << 5,   ///< requires no tau in the event
    kOneMuon          = 1 << 6,   ///< among leptons at least one is a muon
    kTwoMuon          = 1 << 7,   ///< there are exactly two muons in the event
    kMuonsFromZ       = 1 << 8,   ///< require two muons in the event with invariant mass close to Z mass
    kTightMuon        = 1 << 9,   ///< require at least one muon to pass tight id
    kTwoOppositeMuons = 1 << 10,  ///< require exactly two muons with opposite signs
    kMuTrackR0p4      = 1 << 11,  ///< ΔR(mu,track) ≥ 0.4 for each track and each muon
    kMetNoMuTrigger   = 1 << 12,  ///< require MET no mu trigger to fire
    kMet100GeV        = 1 << 13,  ///< MET pt ≥ 100 GeV
    kMet200GeV        = 1 << 14,  ///< MET pt ≥ 200 GeV
    kMetNoMu100GeV    = 1 << 15,  ///< MET no mu pt ≥ 100 GeV
    kMetNoMu200GeV    = 1 << 16,  ///< MET no mu pt ≥ 200 GeV
    kMetJetPhi0p5     = 1 << 17,  ///< Δφ(MET pt,jet pt) ≥ 0.5 for each jet
    kMetNoMuJetPhi0p5 = 1 << 18,  ///< Δφ(MET no mu pt,jet pt) ≥ 0.5 for each jet
    kMuJetR0p4        = 1 << 19,  ///< ΔR(mu,jet) ≥ 0.4 for each jet and each muon
    kHighJetPt100GeV  = 1 << 20,  ///< require highest pt jet to have pt ≥ 100 GeV
    kHighJetChHEF0p1  = 1 << 21,  ///< require highest pt jet to have charged hadron energy fraction above 0.1
    kHighJetNeHEF0p8  = 1 << 22,  ///< require highest pt jet to have neutral hadron energy fraction below 0.8
    kHighJetEta2p4    = 1 << 23,  ///< require highest pt jet to have |eta| ≤ 2.4
    kHighJet          = 1 << 24,  ///< require at least one jet above 100 GeV
  };
  
  /// Default constructor. Creates event cut with specified options.
  EventCut(ECut cutType=kEmpty);
  
  /// Default destructor
  ~EventCut();
  
  /// Setters - allow to modify event cut options
  inline void SetMinNjets(int min){minNjets = min;}
  inline void SetNtracks(int min, int max){minNtracks = min;maxNtracks = max;}
  inline void SetMinNleptons(int min){minNleptons = min;}
  inline void SetMaxNlepton(int max){maxNleptons = max;}
  inline void SetMinNmuons(int min){minNmuons = min;}
  inline void SetMaxNmuons(int max){maxNmuons = max;}
  inline void SetMaxNtau(int max){maxNtau = max;}
  inline void SetRequireMetNoMuTrigger(bool val){metNoMuTrigger = val;}
  inline void SetRequireMuonsFromZ(bool val){muonsFromZ = val;}
  inline void SetRequireMetJetPhi0p5(bool val){metJetPhi = val;}
  inline void SetRequireMetNoMuJetPhi0p5(bool val){metNoMuJetPhi = val;}
  inline void SetRequireMuJetR0p4(bool val){muJetR0p4 = val;}
  inline void SetRequireMuTrackR0p4(bool val){muTrackR0p4 = val;}
  inline void SetRequireHighJet(bool val){highJet = val;}
  inline void SetRequireTightMuon(bool val){tightMuon = val;}
  inline void SetRequireTwoOppositeMuons(bool val){twoOpositeMuons = val;}
  inline void SetMinMetPt(double min){minMetPt = min;}
  inline void SetMinMetNoMuPt(double min){minMetNoMuPt = min;}
  inline void SetHighJetMinPt(double val){highJetMinPt = val;}
  inline void SetHighJetMinChHEF(double val){highJetMinChHEF = val;}
  inline void SetHighJetMaxNeHEF(double val){highJetMaxNeHEF = val;}
  inline void SetHighJetMaxEta(double val){highJetMaxEta = val;}
  inline void SetMinJetMetPhi(double val){minJetMetPhi = val;}
  
  // Getters - give access to options of this event cut
  inline int    GetMinNjets(){return minNjets;}
  inline int    GetMinNtracks(){return minNtracks;}
  inline int    GetMaxNtracks(){return maxNtracks;}
  inline int    GetMinNleptons(){return minNleptons;}
  inline int    GetMaxNleptons(){return maxNleptons;}
  inline int    GetMinNmuons(){return minNmuons;}
  inline int    GetMaxNmuons(){return maxNmuons;}
  inline int    GetMaxNtau(){return maxNtau;}
  inline bool   RequiresMetNoMuTrigger(){return metNoMuTrigger;}
  inline bool   RequiresMuonsFromZ(){return muonsFromZ;}
  inline bool   RequiresMetJetPhi0p5(){return metJetPhi;}
  inline bool   RequiresMetNoMuJetPhi0p5(){return metNoMuJetPhi;}
  inline bool   RequiresMuJetR0p4(){return muJetR0p4;}
  inline bool   RequiresMuTrackR0p4(){return muTrackR0p4;}
  inline bool   RequiresHighJet(){return highJet;}
  inline bool   RequiresTightMuon(){return tightMuon;}
  inline bool   RequiresTwoOppositeMuons(){return twoOpositeMuons;}
  inline double GetMinMetPt(){return minMetPt;}
  inline double GetMinMetNoMuPt(){return minMetNoMuPt;}
  inline double GetHighJetMinPt(){return highJetMinPt;}
  inline double GetHighJetMinChHEF(){return highJetMinChHEF;}
  inline double GetHighJetMaxNeHEF(){return highJetMaxNeHEF;}
  inline double GetHighJetMaxEta(){return highJetMaxEta;}
  inline double GetMinJetMetPhi(){return minJetMetPhi;}
  
private:
  int minNjets;           ///< min number of jets
  int minNtracks;         ///< min number of tracks
  int maxNtracks;         ///< max number of tracks
  int minNleptons;        ///< min number of leptons
  int maxNleptons;        ///< max number of leptons
  int minNmuons;          ///< min number if muons
  int maxNmuons;          ///< max number if muons
  int maxNtau;            ///< max number of tau
  bool metNoMuTrigger;    ///< should require MET no mu trigger
  bool muonsFromZ;        ///< should require two muons in the event with invariant mass close to Z mass
  bool metJetPhi;         ///< should require Δφ(MET pt,jet pt) ≥ 0.5 for each jet
  bool metNoMuJetPhi;     ///< should require Δφ(MET pt,jet pt) ≥ 0.5 for each jet
  bool muJetR0p4;         ///< should require ΔR(mu,jet) ≥ 0.4 for each jet and each muon
  bool muTrackR0p4;       ///< should require ΔR(mu,track) ≥ 0.4 for each track and each muon
  bool highJet;           ///< should require at least one jet above 100 GeV
  bool tightMuon;         ///< should require at least one muon to pass tight id
  bool twoOpositeMuons;   ///< should require exactly two muons with opposite signs
  double minMetPt;        ///< min MET pT
  double minMetNoMuPt;    ///< min MET no mu pT
  double highJetMinPt;    ///< min pt of the highest pt jet
  double highJetMinChHEF; ///< min charged hadron energy fraction of the highest pt jet
  double highJetMaxNeHEF; ///< max neutral hadron energy fraction of the highest pt jet
  double highJetMaxEta;   ///< max |eta| of the highest pt jet
  
  double minJetMetPhi;    ///< min angle between jet and MET
};

#endif /* EventCut_hpp */
