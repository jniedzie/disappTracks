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
  
  /// Returns vector of tracker clusters not assigned to any tracks
  inline Points GetTrackerClusters(){return trackerClusters;}
  
  /// Returns vector of rec pion clusters
  inline Points GetPionClusters(){return pionClusters;}
  
  /// Returns helices of generated pion(s)
  inline Helices GetGenPionHelices() const {return genPionHelices;}
  
  /// Returns vector of pion(s) sim hits
  inline Points GetPionSimHits(){return pionSimHits;}
  
  /// Returns vector of chargino(s) sim hits
  inline Points GetCharginoSimHits(){return charginoSimHits;}
  
  /// Returns tracks of generated chargino(s)
  inline vector<Track> GetGenCharginoTracks() const {return genCharginoTrack;}
  
  // setters
  inline void AddTrack(shared_ptr<Track> track){tracks.push_back(track);}
  inline void AddJet(shared_ptr<Jet> jet){jets.push_back(jet);}
  inline void AddLepton(shared_ptr<Lepton> lepton){leptons.push_back(lepton);}
  inline void AddHelix(Helix helix){helices.push_back(helix);}
  
  inline void SetLumiSection(uint val){lumiSection = val;}
  inline void SetRunNumber(uint val){runNumber = val;}
  inline void SetEventNumber(unsigned long long val){eventNumber = val;}
  
  inline void SetWeight(double val){weight = val;}
  
  inline void SetVertex(unique_ptr<Point> val){vertex = make_unique<Point>(*val);}
  inline void SetNvertices(int n){nVertices = n;}
  inline void SetNjet30(int n){nJet30 = n;}
  inline void SetNjet30a(int n){nJet30a = n;}
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
  
  inline void SetWasTagged(bool val){wasTagged = val;}
  
  // getters
  inline uint  GetLumiSection() const {return lumiSection;}
  inline uint  GetRunNumber() const {return runNumber;}
  inline unsigned long long GetEventNumber() const {return eventNumber;}
  
  inline double GetWeight() const {return weight;}
  
  inline int GetNtracks() const {return (int)tracks.size(); }
  inline int GetNjets() const {return (int)jets.size(); }
  inline int GetNleptons() const {return (int)leptons.size(); }
  inline int GetNcentralJets() const {
    int n=0;
    for(auto j : jets){
      if(!j->IsForward()) n++;
    }
    return n;
  }
  inline int GetNvertices() const {return nVertices;}
  inline unique_ptr<Point> GetVertex() const {return make_unique<Point>(*vertex);}
  inline int GetNjet30() const {return nJet30;}
  inline int GetNjet30a() const {return nJet30a;}
  inline int GetNtau() const {return nTau;}
  inline int GetNhelices() const {return (int)helices.size();}
  
  inline double GetMetSumEt() const {return metSumEt;}
  inline double GetMetPt() const {return metPt;}
  inline double GetMetMass() const {return metMass;}
  inline double GetMetPhi() const {return metPhi;}
  inline double GetMetEta() const {return metEta;}
  
  inline double GetMetNoMuPt() const {return metNoMuPt;}
  inline double GetMetNoMuMass() const {return metNoMuMass;}
  inline double GetMetNoMuPhi() const {return metNoMuPhi;}
  inline double GetMetNoMuEta() const {return metNoMuEta;}
  inline bool HetMetNoMuTrigger(){return metNoMuTrigger;}
  
  inline bool GetGoodVerticesFlag() const {return flag_goodVertices;}
  inline bool GetBadPFmuonFlag() const {return flag_badPFmuon;}
  inline bool GetHBHEnoiseFlag() const {return flag_HBHEnoise;}
  inline bool GetHBHEnoiseIsoFlag() const {return flag_HBHEnoiseIso;}
  inline bool GetEcalDeadCellFlag() const {return flag_EcalDeadCell;}
  inline bool GetEeBadScFlag() const {return flag_eeBadSc;}
  inline bool GetBadChargedCandidateFlag() const {return flag_badChargedCandidate;}
  inline bool GetEcalBadCalibFlag() const {return flag_ecalBadCalib;}
  inline bool GetGlobalTightHalo2016Flag() const {return flag_globalTightHalo2016;}
  
  inline int GetNgenChargino() const {return nGenChargino;}
  inline double GetXsec() const {return xsec;}
  inline double GetWgtSum() const {return wgtsum;}
  inline double GetGenWeight() const {return genWeight;}
  
  inline bool HasFriendData() const {return hasFriendData;}
  inline bool WasTagged() const {return wasTagged;}
  
  inline shared_ptr<Track>  GetTrack(int i) const {return tracks[i];}
  inline shared_ptr<Jet>    GetJet(int i) const {return jets[i];}
  inline shared_ptr<Lepton> GetLepton(int i) const {return leptons[i];}
  inline Helix              GetHelix(int i) const {return helices[i];}
	
	inline vector<shared_ptr<Track>>    GetTracks() const {return tracks;}
  inline vector<shared_ptr<Jet>>      GetJets() const {return jets;}
  inline vector<shared_ptr<Lepton>>   GetLeptons() const {return leptons;}
  inline Helices                      GetHelices() const {return helices;}

  /**
   Returns tracker clusters. End-caps are included or not based on `include_endcaps` option
   in the config. Pion clusters are removed from the collection or not, according to
   `fit_noise_clusters_only` option in the config.
   */
  Points GetClusters();
  
private:
  vector<shared_ptr<Track>>  tracks;   ///< Vector of isolated tracks
  vector<shared_ptr<Jet>>    jets;     ///< Vector of jets
  vector<shared_ptr<Lepton>> leptons;  ///< Vector of leptons
  Helices  helices;  ///< Parameters of the fitted helices (one per track)
  
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
  bool hasFriendData;       ///< Indicates whether of not additional info is available for this event
  bool wasTagged;           ///< Indicates whether of not tagger was ran on this event
  Points trackerClusters;   ///< All reconstructed tracker clusters not assigned to any track
  Points pionClusters;      ///< Reconstructed clusters associated with generated pion(s)
  Points pionSimHits;       ///< Sim hits associated with generated pion(s) coming from chargino decay vertex
  Points charginoSimHits;   ///< Sim hits associated with generated chargino(s)
  Helices genPionHelices;   ///< Helix representing gen-level pion(s)
  
  vector<Track> genCharginoTrack; ///< Gen-level information about charginos
  
  TTree *friendTree;
  
  friend class EventProcessor;
};

#endif /* Event_hpp */
