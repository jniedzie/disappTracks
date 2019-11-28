//
//  CutsManager.hpp
//
//  Created by Jeremi Niedziela on 11/04/2019.
//

#ifndef CutsManager_hpp
#define CutsManager_hpp

#include "Helpers.hpp"

#include "EventProcessor.hpp"
#include "TrackProcessor.hpp"
#include "JetProcessor.hpp"
#include "LeptonProcessor.hpp"

class CutsManager {
public:
  /// Default constructor
  CutsManager();
  
  /// Ddefault destructor
  ~CutsManager();
  
  void GetCuts(EventCut &eventCut, TrackCut &trackCut, JetCut &jetCut, LeptonCut &leptonCut);
  void GetZmumuCuts(EventCut &eventCut, TrackCut &trackCut, JetCut &jetCut, LeptonCut &leptonCut);
  void GetWmunuCuts(EventCut &eventCut, TrackCut &trackCut, JetCut &jetCut, LeptonCut &leptonCut);
  void GetLowMetCuts(EventCut &eventCut, TrackCut &trackCut, JetCut &jetCut, LeptonCut &leptonCut);
  
private:
  
};

#endif /* CutsManager_hpp */
