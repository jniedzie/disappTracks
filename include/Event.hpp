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

class Event;

class Events {
public:
  enum EDataType{
    kBackground,
    kSignal,
    kData
  };
  
  /// Default constructor. Loads events from ROOT tree
  /// \param fileName Path to the ROOT file with ntuples from which events will be loaded
  /// \param dataType Event weigth will be calculated differently background, signal and data
  /// \param maxNevents Load just maxNevents from file and then stop
  /// \param iSig Signal type for correct cross section assignment
  Events(std::string fileName, EDataType dataType=kBackground, int maxNevents=-1,
         ESignal iSig=kNsignals);
  
  /// Empty constructor. Creates an empty class with no events.
  Events();
  
  /// Default destructor
  ~Events();
  
  /// Adds event to the collection of events
  /// \param event Event object to be added to the collection
  void AddEvent(Event *event){events.push_back(event);}
  
  /// Adds events from specified path to the existing events collection
  /// \param fileName Path to the ROOT file with ntuples from which events will be loaded
  /// \param dataType Event weigth will be calculated differently for background, signal and data
  /// \param maxNevents Load just maxNevents from file and then stop
  /// \param iSig Signal type for correct cross section assignment
  void AddEventsFromFile(std::string fileName, EDataType dataType=kBackground, int maxNevents=-1,
                         ESignal iSig=kNsignals);
  
  void SaveToTree(std::string fileName);
  
  /// Applies cuts in this order: track, jet, lepton, event
  /// \param eventCut   Cuts to be applied to events
  /// \param trackCut   Cuts to be applied to tracks
  /// \param jetCut     Cuts to be applied to jets
  /// \param leptonCut  Cuts to be applied to leptons
  /// \return Returns collection of events passing cuts, containing only tracks, jets and leptons that passed all cuts
  Events* ApplyCuts(EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut);
  
  /// Applies cuts on events
  /// \param cut   Cuts to be applied to events
  /// \return Returns collection of events passing all cuts
  Events* ApplyEventCut(EventCut *cut);
  
  /// Applies cuts on tracks
  /// \param cut   Cuts to be applied to tracks
  /// \return Returns collection of events, containing only those tracks that passed all cuts
  Events* ApplyTrackCut(TrackCut *cut);
  
  /// Applies cuts on jets
  /// \param cut   Cuts to be applied to jets
  /// \return Returns collection of events, containing only those jets that passed all cuts
  Events* ApplyJetCut(JetCut *cut);
  
  /// Applies cuts on leptons
  /// \param cut   Cuts to be applied to leptons
  /// \return Returns collection of events, containing only those leptons that passed all cuts
  Events* ApplyLeptonCut(LeptonCut *cut);
  
  /// Returns number of events in this collection
  inline int size(){return (int)events.size();}
  
  double weightedSize();
  
  /// Returns the event with given index
  Event* At(int index){return events[index];}

private:
  /// Copy constructor. Copies all events in the collection.
  Events(const Events &event);
  
  std::vector<Event*> events; ///< Vector of events
};

//---------------------------------------------------------------------------------------
// Single event class
//---------------------------------------------------------------------------------------

/// Representation of a single event.
///
/// This class contains global information about the event, as well as
/// collection of tracks, jets and leptons belonging to this event.
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
  
  inline void SetGoodVerticesFlag(bool val){flag_goodVertices = val;}
  inline void SetBadPFmuonFlag(bool val){flag_badPFmuon = val;}
  inline void SetHBHEnoiseFlag(bool val){flag_HBHEnoise = val;}
  inline void SetHBHEnoiseIsoFlag(bool val){flag_HBHEnoiseIso = val;}
  inline void SetEcalDeadCellFlag(bool val){flag_EcalDeadCell = val;}
  inline void SetEeBadScFlag(bool val){flag_eeBadSc = val;}
  inline void SetBadChargedCandidateFlag(bool val){flag_badChargedCandidate = val;}
  inline void SetEcalBadCalibFlag(bool val){flag_ecalBadCalib = val;}
  inline void SetGlobalTightHalo2016Flag(bool val){flag_globalTightHalo2016 = val;}
  
  // getters
  inline double GetWeight(){return weight;}
  
  inline unsigned long GetNtracks(){return tracks.size(); }
  inline unsigned long GetNjets(){return jets.size(); }
  inline unsigned long GetNcentralJets(){
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
};

#endif /* Event_hpp */
