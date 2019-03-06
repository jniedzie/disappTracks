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

class EventProcessor {
public:
  EventProcessor();
  
  ~EventProcessor();
  
  /// Removes tracks that don't pass the cut from the tracks collection
  void ApplyTrackCut(shared_ptr<Event> event, const unique_ptr<TrackCut> &cut);
  
  /// Removes jets that don't pass the cut from the jets collection
  void ApplyJetCut(shared_ptr<Event> event, const unique_ptr<JetCut> &cut);
  
  /// Removes leptons that don't pass the cut from the leptons collection
  void ApplyLeptonCut(shared_ptr<Event> event, const unique_ptr<LeptonCut> &cut);
  
  /// Check if event passes a cut
  bool IsPassingCut(const shared_ptr<Event> event, const unique_ptr<EventCut> &cut);
  
  /// Link class variables to branches of a specified tree
  /// \param tree Tree from which event parameters will be read
  void SetupBranchesForReading(TTree *tree);
  
  /// Returns an event with parameters and objects read from tree previously set with SetupBranchesForReading(..)
  shared_ptr<Event> GetEventFromTree(xtracks::EDataType dataType, int setIter);
  
  /// Link class variables to branches of a specified tree
  /// \param tree Tree to which evnet parameters and objects will be saved
  void SetupBranchesForWriting(TTree *tree);
  
  /// Writes an event to the tree previously set with SetupBranchesForWriting(...)
  void SaveEventToTree(shared_ptr<Event> event);
  
private:
  unique_ptr<TrackProcessor>  trackProcessor;
  unique_ptr<JetProcessor>    jetProcessor;
  unique_ptr<LeptonProcessor> leptonProcessor;
  
  map<string, float>                singleValuesFloat;  ///< Float per-event variables in the current entry
  map<string, int>                  singleValuesInt;     ///< Int per-event variables in the current entry
  map<string, uint>                 singleValuesUint;    ///< uint per-event variables in the current entry
  map<string, unsigned long long>   singleValuesUlonglong;///<long per-event variables in the current entry
  
  vector<string> singleNamesFloat;     ///< Names or float per-event variables
  vector<string> singleNamesInt;       ///< Names or int per-event variables
  vector<string> singleNamesUint;
  vector<string> singleNamesUlongLong;
  
};

#endif /* EventProcessor_hpp */
