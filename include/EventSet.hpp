//
//  EventSet.hpp
//
//  Created by Jeremi Niedziela on 08/11/2018.
//

#ifndef EventSet_hpp
#define EventSet_hpp

#include "Helpers.hpp"
#include "Event.hpp"

#include <optional>

class EventSet {
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
  EventSet(string fileName, EDataType dataType=kBackground, int maxNevents=-1, ESignal iSig=kNsignals);
  
  /// Empty constructor. Creates an empty class with no events.
  EventSet();
  
  /// Copy constructor. Copies all events in the collection.
  EventSet(const EventSet &event);
  
  /// Default destructor
  ~EventSet();
  
  static void LoadEventsFromFiles(vector<shared_ptr<EventSet>> &eventsSignal,
                                  vector<shared_ptr<EventSet>> &eventsBackground,
                                  vector<shared_ptr<EventSet>> &eventsData,
                                  string prefix="");
  
  static void SaveEventsToFiles(vector<shared_ptr<EventSet>> &eventsSignal,
                                vector<shared_ptr<EventSet>> &eventsBackground,
                                vector<shared_ptr<EventSet>> &eventsData,
                                string prefix="after_L/");
  
 
  
  static void PrintYields(vector<shared_ptr<EventSet>> &eventsSignal,
                          vector<shared_ptr<EventSet>> &eventsBackground,
                          vector<shared_ptr<EventSet>> &eventsData);
  
  /// Applies cuts in this order: track, jet, lepton, event to three sets of events: signal, background and data.
  /// \param eventsSignal Reference to vector of signal event sets
  /// \param eventsBackground Reference to vector of background event sets
  /// \param eventsData Reference to vector of data event sets
  /// \param eventCut   Cuts to be applied to events
  /// \param trackCut   Cuts to be applied to tracks
  /// \param jetCut     Cuts to be applied to jets
  /// \param leptonCut  Cuts to be applied to leptons
  static void ApplyCuts(vector<shared_ptr<EventSet>> &eventsSignal,
                        vector<shared_ptr<EventSet>> &eventsBackground,
                        vector<shared_ptr<EventSet>> &eventsData,
                        EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut);
  
 
  
  /// Returns number of events in this collection
  inline int size(){return (int)events.size();}
  
  double weightedSize();
  
  /// Returns the event with given index
  shared_ptr<Event> At(int index){return events[index];}
  
private:
  vector<std::shared_ptr<Event>> events; ///< Vector of events
  
  /// Adds event to the collection of events
  /// \param event Event object to be added to the collection
  void AddEvent(shared_ptr<Event> event){events.push_back(event);}
  
  /// Adds events from specified path to the existing events collection
  /// \param fileName Path to the ROOT file with ntuples from which events will be loaded
  /// \param dataType Event weigth will be calculated differently for background, signal and data
  /// \param maxNevents Load just maxNevents from file and then stop
  /// \param iSig Signal type for correct cross section assignment
  void AddEventsFromFile(string fileName, EDataType dataType=kBackground, int maxNevents=-1,
                         ESignal iSig=kNsignals);
  
  /// Applies cuts in this order: track, jet, lepton, event
  /// \param eventCut   Cuts to be applied to events
  /// \param trackCut   Cuts to be applied to tracks
  /// \param jetCut     Cuts to be applied to jets
  /// \param leptonCut  Cuts to be applied to leptons
  /// \return Returns collection of events passing cuts, containing only tracks, jets and leptons that passed all cuts
  shared_ptr<EventSet> ApplyCuts(EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut);
  
  /// Applies cuts on events
  /// \param cut   Cuts to be applied to events
  /// \return Returns collection of events passing all cuts
  shared_ptr<EventSet> ApplyEventCut(EventCut *cut);
  
  /// Applies cuts on tracks
  /// \param cut   Cuts to be applied to tracks
  /// \return Returns collection of events, containing only those tracks that passed all cuts
  shared_ptr<EventSet> ApplyTrackCut(TrackCut *cut);
  
  /// Applies cuts on jets
  /// \param cut   Cuts to be applied to jets
  /// \return Returns collection of events, containing only those jets that passed all cuts
  shared_ptr<EventSet> ApplyJetCut(JetCut *cut);
  
  /// Applies cuts on leptons
  /// \param cut   Cuts to be applied to leptons
  /// \return Returns collection of events, containing only those leptons that passed all cuts
  shared_ptr<EventSet> ApplyLeptonCut(LeptonCut *cut);
  
  void SaveToTree(string fileName);
};


#endif /* EventSet_hpp */
