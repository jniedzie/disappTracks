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
  
private:
  unique_ptr<TrackProcessor>  trackProcessor;
  unique_ptr<JetProcessor>    jetProcessor;
  unique_ptr<LeptonProcessor> leptonProcessor;
};

#endif /* EventProcessor_hpp */
