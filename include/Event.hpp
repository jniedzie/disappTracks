//
//  Event.hpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
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

#include <utility>

/// This class represents a single event. It contains global information about the event, as well as
/// collection of tracks, jets and leptons belonging to this event.
/// It also checks if this event passes some set of selection criteria represented by EventCut object.
class Event{
public:
  /// Default constructor
  Event();
  
  /// Copy constructor
  Event(const Event &event);
  
  /// Default destructor
  ~Event();
  
  /// Prints basic information about the event
  void Print();
  
  /// Removes tracks that don't pass the cut from the tracks collection
  void ApplyTrackCut(const unique_ptr<TrackCut> &cut);
  
  /// Removes jets that don't pass the cut from the jets collection
  void ApplyJetCut(const unique_ptr<JetCut> &cut);
  
  /// Removes leptons that don't pass the cut from the leptons collection
  void ApplyLeptonCut(const unique_ptr<LeptonCut> &cut);
  
  /// Check if event passes a cut
  bool IsPassingCut(const unique_ptr<EventCut> &cut);
  
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
  
  inline void SetGoodVerticesFlag(bool val){flag_goodVertices = val;}
  inline void SetBadPFmuonFlag(bool val){flag_badPFmuon = val;}
  inline void SetHBHEnoiseFlag(bool val){flag_HBHEnoise = val;}
  inline void SetHBHEnoiseIsoFlag(bool val){flag_HBHEnoiseIso = val;}
  inline void SetEcalDeadCellFlag(bool val){flag_EcalDeadCell = val;}
  inline void SetEeBadScFlag(bool val){flag_eeBadSc = val;}
  inline void SetBadChargedCandidateFlag(bool val){flag_badChargedCandidate = val;}
  inline void SetEcalBadCalibFlag(bool val){flag_ecalBadCalib = val;}
  inline void SetGlobalTightHalo2016Flag(bool val){flag_globalTightHalo2016 = val;}
  
  inline void SetNgenChargino(int val){nGenChargino = val;}
  inline void SetXsec(double val){xsec = val;}
  inline void SetWgtSum(double val){wgtsum = val;}
  inline void SetGenWeight(double val){genWeight = val;}
  
  // getters
  inline double GetWeight(){return weight;}
  
  inline int GetNtracks(){return (int)tracks.size(); }
  inline int GetNjets(){return (int)jets.size(); }
  inline int GetNcentralJets(){
    int n=0;
    for(auto j : jets){
      if(!j->IsForward()) n++;
    }
    return n;
  }
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
  
  inline bool GetGoodVerticesFlag(){return flag_goodVertices;}
  inline bool GetBadPFmuonFlag(){return flag_badPFmuon;}
  inline bool GetHBHEnoiseFlag(){return flag_HBHEnoise;}
  inline bool GetHBHEnoiseIsoFlag(){return flag_HBHEnoiseIso;}
  inline bool GetEcalDeadCellFlag(){return flag_EcalDeadCell;}
  inline bool GetEeBadScFlag(){return flag_eeBadSc;}
  inline bool GetBadChargedCandidateFlag(){return flag_badChargedCandidate;}
  inline bool GetEcalBadCalibFlag(){return flag_ecalBadCalib;}
  inline bool GetGlobalTightHalo2016Flag(){return flag_globalTightHalo2016;}
  
  inline int GetNgenChargino(){return nGenChargino;}
  inline double GetXsec(){return xsec;}
  inline double GetWgtSum(){return wgtsum;}
  inline double GetGenWeight(){return genWeight;}
  
  inline Track*  GetTrack(int i){return tracks[i];}
  inline Jet*    GetJet(int i){return jets[i];}
  inline Lepton* GetLepton(int i){return leptons[i];}
  
private:
  vector<Track*>  tracks;   ///< Vector of isolated tracks
  vector<Jet*>    jets;     ///< Vector of jets
  vector<Lepton*> leptons;  ///< Vector of leptons
  
  double weight;          ///< Weight for this event resulting from lumi, cross section and generator weights
  
  int nVertices;          ///< Number of good verices
  int nJet30;             ///< Number of jets with pt > 30, |eta|<2.4
  int nJet30a;            ///< Number of jets with pt > 30, |eta|<4.7
  int nLepton;            ///< Number of leptons
  int nTau;               ///< Number of tau leptons
  
  double  metSumEt;       ///< Sum of missing transverse energy
  double  metPt;          ///< MET transverse momentum
  double  metEta;         ///< MET pseudorapidity
  double  metPhi;         ///< MET polar angle
  double  metMass;        ///< MET mass
  
  bool    metNoMuTrigger; ///< Did HLT_BIT_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight trigger fired for this event?
  double  metNoMuPt;      ///< MET transverse momentum (calculated without muons)
  double  metNoMuEta;     ///< MET pseudorapidity (calculated without muons)
  double  metNoMuPhi;     ///< MET polar angle (calculated without muons)
  double  metNoMuMass;    ///< MET mass (calculated without muons)
  
  // flags for different filters (true means that event passes the filter)
  bool flag_goodVertices;
  bool flag_badPFmuon;
  bool flag_HBHEnoise;
  bool flag_HBHEnoiseIso;
  bool flag_EcalDeadCell;
  bool flag_eeBadSc;
  bool flag_badChargedCandidate; // an exception, this is zero when event passes filter
  bool flag_ecalBadCalib;
  bool flag_globalTightHalo2016;
  
  int nGenChargino;
  double xsec;
  double wgtsum;
  double genWeight;
};

#endif /* Event_hpp */
