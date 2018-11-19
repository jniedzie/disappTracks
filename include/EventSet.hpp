//
//  EventSet.hpp
//
//  Created by Jeremi Niedziela on 08/11/2018.
//

#ifndef EventSet_hpp
#define EventSet_hpp

#include "Helpers.hpp"
#include "Event.hpp"

/// This class contains three vectors, for signal, background and data events. In each of those vectors,
/// there's a vector of events of given type (for instance, in the background vector there will be a
/// vector of events for Z->mumu, another vector with events for Z->nunu etc.).
/// The class allows to load those events from files, apply cuts, get basic information about each of
/// the collections and store filtered events in ROOT files.
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
  
  /// Loads signal, background and data events according to settings in the config.
  /// \param prefix If specified, prefix will be appended at the end of the path, before "tree.root"
  void LoadEventsFromFiles(string prefix="");
  
  /// Loads events from one tree only.
  /// \param dataType Specifies whether signal, background or data events should be loaded.
  /// \param setIter Specifies which set should be loaded (e.g. kZmumuJets for Z->mumu)
  /// \param prefix If specified, prefix will be appended at the end of the path, before "tree.root"
  void LoadEventsFromFiles(EDataType dataType, int setIter, string prefix="");
  
  /// Saves signal, background and data events in the ROOT files, according to settings in the config.
  /// \param prefix If specified, prefix will be appended at the end of the path, before "tree.root"
  void SaveEventsToFiles(string prefix="after_L/");
  
  /// Prints yields of signal, background and data events, as well as S/sqrt(S+B) ratio.
  void PrintYields();
  
  vector<double> GetSignificance(bool inData=false);
  
  /// Applies cuts in this order: track, jet, lepton, event to three sets of events: signal, background and data.
  /// \param eventCut   Cuts to be applied to events
  /// \param trackCut   Cuts to be applied to tracks
  /// \param jetCut     Cuts to be applied to jets
  /// \param leptonCut  Cuts to be applied to leptons
  void ApplyCuts(const unique_ptr<EventCut> &eventCut,const unique_ptr<TrackCut> &trackCut,
                 const unique_ptr<JetCut> &jetCut,const unique_ptr<LeptonCut> &leptonCut);
  
  void DrawStandardPlots(string prefix="");
  void DrawPerLayerPlots();
  
  /// Returns number of events in this collection
  /// \param dataType Specifies which data type should be taken into concideration: signal, background or data
  /// \param setIter Specifies which set for given dataType to look at (e.g. kZmumuJets for Z->mumu)
  int size(EDataType dataType, int setIter);
  
  /// Returns weighted size of a collection (including luminosity, cross section and generator weights.
  /// \param dataType Specifies which data type should be taken into concideration: signal, background or data
  /// \param setIter Specifies which set for given dataType to look at (e.g. kZmumuJets for Z->mumu)
  double weightedSize(EDataType dataType, int setIter);
  
  /// Returns the event with given index
  /// \param dataType Specifies from which data type: signal, background or data the events should be taken
  /// \param setIter Specifies which set for given dataType to look at (e.g. kZmumuJets for Z->mumu)
  shared_ptr<Event> At(EDataType dataType, int setIter, int index);
  
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
  void AddEvent(shared_ptr<Event> event, EDataType dataType, int setIter);
  
  /// Adds events from specified path to the existing events collection
  /// \param fileName Path to the ROOT file with ntuples from which events will be loaded
  /// \param dataType Event weigth will be calculated differently for background, signal and data
  /// \param maxNevents Load just maxNevents from file and then stop
  /// \param setIter Data set type (i.a. for correct cross section assignment). Will be casted to ESignal, EBackground or EData, depending on dataType parameter provided
  void AddEventsFromFile(string fileName, EDataType dataType=kBackground, int maxNevents=-1, int setIter=-1);
 
  void SaveToTree(string fileName, EDataType type, int setIter);
};


#endif /* EventSet_hpp */
