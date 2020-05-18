//
//  EventProcessor.hpp
//
//  Created by Jeremi Niedziela on 30/01/2019.
//

#ifndef EventProcessor_hpp
#define EventProcessor_hpp

#include "Helpers.hpp"
#include "Event.hpp"
#include "TrackProcessor.hpp"
#include "JetProcessor.hpp"
#include "LeptonProcessor.hpp"

class EventProcessor;
extern EventProcessor eventProcessor;

class EventProcessor {
public:
  EventProcessor();
  
  ~EventProcessor();
  
  /// Removes tracks that don't pass the cut from the tracks collection
  void ApplyTrackCut(shared_ptr<Event> event, const TrackCut &cut, vector<int> *trackCutReasons=nullptr);
  
  /// Removes jets that don't pass the cut from the jets collection
  void ApplyJetCut(shared_ptr<Event> event, const JetCut &cut);
  
  /// Removes leptons that don't pass the cut from the leptons collection
  void ApplyLeptonCut(shared_ptr<Event> event, const LeptonCut &cut);
  
  /// Check if event passes a cut
  bool IsPassingCut(const shared_ptr<Event> event, const EventCut &cut, vector<int> *cutReasons = nullptr);
  
  /// Link class variables to branches of a specified tree
  /// \param tree Tree from which event parameters will be read
  void SetupBranchesForReading(TTree *tree, TTree *friendTree=nullptr, TTree *prefireTree=nullptr);
  
  /// Returns an event with parameters and objects read from tree previously set with SetupBranchesForReading(..)
  shared_ptr<Event> GetEventFromTree(xtracks::EDataType dataType, int setIter, int year,
                                     TTree *friendTree=nullptr, TTree *prefireTree=nullptr, TH1D *metWeights=nullptr);
  
  /// Link class variables to branches of a specified tree
  /// \param tree Tree to which evnet parameters and objects will be saved
  void SetupBranchesForWriting(TTree *tree);
  
  /// Writes an event to the tree previously set with SetupBranchesForWriting(...)
  void SaveEventToTree(shared_ptr<Event> event);
  
  /// Read friend tree for given dataType and setIter
  TTree* GetFriendTree(xtracks::EDataType dataType, int setIter);
  
  ///
  
  vector<shared_ptr<Event>> survivingEvents;
  
private:
  
  map<string, float>                singleValuesFloat;  ///< Float per-event variables in the current entry
  
  static const int nMax = 10000;   ///< Maximum supported number of tracks per event
  map<string, float[nMax]>       arrayValuesFloat;  ///< Float per-event variables in the current entry
  
  map<string, int>                  singleValuesInt;     ///< Int per-event variables in the current entry
  map<string, uint>                 singleValuesUint;    ///< uint per-event variables in the current entry
  map<string, unsigned long long>   singleValuesUlonglong;///<long per-event variables in the current entry
  
  vector<string> singleNamesFloat;      ///< Names or float per-event variables
  vector<string> arrayNamesFloat;      ///< Names or float per-event variables
  vector<string> singleNamesInt;        ///< Names or int per-event variables
  vector<string> singleNamesUint;       ///< Names or uint per-event variables
  vector<string> singleNamesUlongLong;  ///< Names or long per-event variables
  
  
  
  
  map<string, vector<double>* > arrayValuesFriendFloat;  ///< Float variables in the current entry
  map<string, vector<int>* >   arrayValuesFriendInt;    ///< Int variables in the current entry
  
  vector<string> arrayNamesFriendFloat;     ///< Names or float variables
  vector<string> arrayNamesFriendInt;       ///< Names or int variables
  
  int nGenChargino;
  float prefireWeight;
};

#endif /* EventProcessor_hpp */
