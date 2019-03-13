//
//  Event.hpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#ifndef Event_hpp
#define Event_hpp

#include "Helpers.hpp"
#include "EventCut.hpp"
#include "TrackProcessor.hpp"
#include "Jet.hpp"
#include "JetCut.hpp"
#include "Lepton.hpp"
#include "LeptonCut.hpp"
#include "HelixProcessor.hpp"

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
  
  /// Look for a friend tree and load gen-level information + all tracker clusters
  void LoadAdditionalInfo();
  
  /// Returns vector of tracker clusters not assigned to any tracks
  inline shared_ptr<vector<Point>> GetTrackerClusters(){return trackerClusters;}
  
  /// Returns helices of generated pion(s)
  inline shared_ptr<vector<unique_ptr<Helix>>> GetGenPionHelices(){return genPionHelices;}
  
  /// Returns vector of pion(s) sim hits
  inline shared_ptr<vector<Point>> GetPionSimHits(){return pionSimHits;}
  
  /// Returns vector of chargino(s) sim hits
  inline shared_ptr<vector<Point>> GetCharginoSimHits(){return charginoSimHits;}
  
  // setters
  inline void AddTrack(shared_ptr<Track> track){tracks.push_back(track);}
  inline void AddJet(shared_ptr<Jet> jet){jets.push_back(jet);}
  inline void AddLepton(shared_ptr<Lepton> lepton){leptons.push_back(lepton);}
  inline void AddHelix(shared_ptr<Helix> helix){helices.push_back(helix);}
  
  inline void SetLumiSection(uint val){lumiSection = val;}
  inline void SetRunNumber(uint val){runNumber = val;}
  inline void SetEventNumber(unsigned long long val){eventNumber = val;}
  
  inline void SetWeight(double val){weight = val;}
  
  inline void SetVertex(unique_ptr<Point> val){vertex = make_unique<Point>(val);}
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
  
  inline void SetDataType(xtracks::EDataType val){dataType = val;}
  inline void SetSetIter(int val){setIter = val;}
  
  // getters
  inline uint  GetLumiSection(){return lumiSection;}
  inline uint  GetRunNumber(){return runNumber;}
  inline unsigned long long GetEventNumber(){return eventNumber;}
  
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
  inline unique_ptr<Point> GetVertex(){return make_unique<Point>(vertex);}
  inline int GetNjet30(){return nJet30;}
  inline int GetNjet30a(){return nJet30a;}
  inline int GetNlepton(){return nLepton;}
  inline int GetNtau(){return nTau;}
  inline int GetNhelices(){return (int)helices.size();}
  
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
  
  inline shared_ptr<Track>  GetTrack(int i){return tracks[i];}
  inline shared_ptr<Jet>    GetJet(int i){return jets[i];}
  inline shared_ptr<Lepton> GetLepton(int i){return leptons[i];}
  inline shared_ptr<Helix>  GetHelix(int i){return helices[i];}
	
	inline vector<shared_ptr<Track>>   GetTracks(){return tracks;}
  inline vector<shared_ptr<Jet>>     GetJets(){return jets;}
  inline vector<shared_ptr<Lepton>>  GetLeptons(){return leptons;}
  inline vector<shared_ptr<Helix>>   GetHelices(){return helices;}
private:
  vector<shared_ptr<Track>>  tracks;   ///< Vector of isolated tracks
  vector<shared_ptr<Jet>>    jets;     ///< Vector of jets
  vector<shared_ptr<Lepton>> leptons;  ///< Vector of leptons
  vector<shared_ptr<Helix>>  helices;  ///< Parameters of the fitted helices (one per track)
  
  xtracks::EDataType dataType;     ///< Type of the event (signal/background/data)
  int setIter;            ///< Iterator of the dataset (e.g. which type of background it is)
  
  uint lumiSection;        ///< ID of lumi section
  uint runNumber;          ///< ID of the run
  unsigned long long eventNumber;        ///< ID of the event
  
  double weight;          ///< Weight for this event resulting from lumi, cross section and generator weights
  
  int nVertices;          ///< Number of good verices
  unique_ptr<Point> vertex; ///< Primary vertex of the event
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
  
  // Additional info from RAW-RECO files
  shared_ptr<vector<Point>> trackerClusters;  ///< All reconstructed tracker clusters not assigned to any track
  shared_ptr<vector<Point>> pionSimHits;  ///< Sim hits associated with generated pion(s) coming from chargino decay vertex
  shared_ptr<vector<Point>> charginoSimHits; ///< Sim hits associated with generated chargino(s)
  shared_ptr<vector<unique_ptr<Helix>>> genPionHelices; ///< Helix representing gen-level pion(s)
  
  unique_ptr<TrackProcessor> trackProcessor;
  
  friend class EventProcessor;
};

#endif /* Event_hpp */
