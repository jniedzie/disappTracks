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
  
  /// Empty constructor. Creates an empty class with no events.
  EventSet();
  
  /// Copy constructor. Copies all events in the collection.
  EventSet(const EventSet &event);
  
  /// Default destructor
  ~EventSet();
  
  void LoadEventsFromFiles(string prefix="");
  
  void SaveEventsToFiles(string prefix="after_L/");
  
 
  void PrintYields();
  
  /// Applies cuts in this order: track, jet, lepton, event to three sets of events: signal, background and data.
  /// \param eventCut   Cuts to be applied to events
  /// \param trackCut   Cuts to be applied to tracks
  /// \param jetCut     Cuts to be applied to jets
  /// \param leptonCut  Cuts to be applied to leptons
  void ApplyCuts(EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut);
  
 
  
  /// Returns number of events in this collection
  int size(EDataType dataType, int setIter){
    if(dataType == kSignal){
      return (int)eventsSignal[(ESignal)setIter].size();
    }
    else if(dataType == kBackground){
      return (int)eventsBackground[(EBackground)setIter].size();
    }
    else if(dataType == kData){
      return (int)eventsData[(EData)setIter].size();
    }
    else{
      throw out_of_range("Unknown data type provided");
    }
    return 1;
  }
  
  double weightedSize(EDataType type, int setIter);
  
  /// Returns the event with given index
  shared_ptr<Event> At(EDataType dataType, int setIter, int index){
    if(dataType == kSignal){
      return eventsSignal[(ESignal)setIter][index];
    }
    else if(dataType == kBackground){
      return eventsBackground[(EBackground)setIter][index];
    }
    else if(dataType == kData){
      return eventsData[(EData)setIter][index];
    }
    else{
      throw out_of_range("Unknown data type provided");
    }
  }
  
private:
  /// Default constructor. Loads events from ROOT tree
  /// \param fileName Path to the ROOT file with ntuples from which events will be loaded
  /// \param dataType Event weigth will be calculated differently background, signal and data
  /// \param maxNevents Load just maxNevents from file and then stop
  /// \param iSig Signal type for correct cross section assignment
  EventSet(string fileName, EDataType dataType=kBackground, int maxNevents=-1, ESignal iSig=kNsignals);
  
  vector<vector<shared_ptr<Event>>> eventsSignal;     ///< Vector of signal events - [ESignal][iEvent]
  vector<vector<shared_ptr<Event>>> eventsBackground; ///< Vector of backgrnd events - [EBackkground][iEvent]
  vector<vector<shared_ptr<Event>>> eventsData;       ///< Vector of data events - [EData][iEvent]
  
  /// Adds event to the collection of events
  /// \param event Event object to be added to the collection
  void AddEvent(shared_ptr<Event> event, EDataType dataType, int setIter){
    
    if(dataType == kSignal){
      eventsSignal[(ESignal)setIter].push_back(event);
    }
    else if(dataType == kBackground){
      eventsBackground[(EBackground)setIter].push_back(event);
    }
    else if(dataType == kData){
      eventsData[(EData)setIter].push_back(event);
    }
    else{
      throw out_of_range("Unknown data type provided");
    }
  }
  
  /// Adds events from specified path to the existing events collection
  /// \param fileName Path to the ROOT file with ntuples from which events will be loaded
  /// \param dataType Event weigth will be calculated differently for background, signal and data
  /// \param maxNevents Load just maxNevents from file and then stop
  /// \param setIter Data set type (i.a. for correct cross section assignment). Will be casted to ESignal, EBackground or EData, depending on dataType parameter provided
  void AddEventsFromFile(string fileName, EDataType dataType=kBackground, int maxNevents=-1, int setIter=-1);
  
  /// Applies cuts in this order: track, jet, lepton, event
  /// \param eventCut   Cuts to be applied to events
  /// \param trackCut   Cuts to be applied to tracks
  /// \param jetCut     Cuts to be applied to jets
  /// \param leptonCut  Cuts to be applied to leptons
  /// \return Returns collection of events passing cuts, containing only tracks, jets and leptons that passed all cuts
  vector<shared_ptr<Event>> ApplyCuts(EventCut *eventCut, TrackCut *trackCut, JetCut *jetCut, LeptonCut *leptonCut, EDataType dataType, int setIter);
  
  /// Applies cuts on events
  /// \param cut   Cuts to be applied to events
  /// \return Returns collection of events passing all cuts
  vector<shared_ptr<Event>> ApplyEventCut(vector<shared_ptr<Event>> events, EventCut *cut);
  
  /// Applies cuts on tracks
  /// \param cut   Cuts to be applied to tracks
  /// \return Returns collection of events, containing only those tracks that passed all cuts
  vector<shared_ptr<Event>> ApplyTrackCut(vector<shared_ptr<Event>> events, TrackCut *cut);
  
  /// Applies cuts on jets
  /// \param cut   Cuts to be applied to jets
  /// \return Returns collection of events, containing only those jets that passed all cuts
  vector<shared_ptr<Event>> ApplyJetCut(vector<shared_ptr<Event>> events, JetCut *cut);
  
  /// Applies cuts on leptons
  /// \param cut   Cuts to be applied to leptons
  /// \return Returns collection of events, containing only those leptons that passed all cuts
  vector<shared_ptr<Event>> ApplyLeptonCut(vector<shared_ptr<Event>> events, LeptonCut *cut);
  
  void SaveToTree(string fileName, EDataType type, int setIter);
};


#endif /* EventSet_hpp */
