//
//  TrackProcessor.hpp
//
//  Created by Jeremi Niedziela on 30/01/2019.
//

#ifndef TrackProcessor_hpp
#define TrackProcessor_hpp

#include "Helpers.hpp"
#include "Track.hpp"

/// Class description
class TrackProcessor {
public:
  /// Default constructor
  TrackProcessor();
  
  /// Default destructor
  ~TrackProcessor();
  
  /// Returns a random track with basic properties filled in (eta, phi, decay point given number of layers it went through)
  /// \param nLayers Number of layers which track passed before decaying
  /// \param maxEta Pseudorapidity limit for the track
  shared_ptr<Track> GetRandomTrack(int nLayers, double maxEta);
  
  /// Check if track passes selection criteria
  /// \param track Tracks to which cuts should be applied
  /// \param cut Tracks selection criteria to be checked
  bool IsPassingCut(const shared_ptr<Track> track,
                    const unique_ptr<TrackCut> &cut);
  
  /// Link class variables to branches of a specified tree
  /// \param tree Tree from which track parameters will be read
  void SetupBranches(TTree *tree);
  
  /// Returns a vector of tracks with parameters read from tree previously set with SetupBranches(..)
  vector<shared_ptr<Track>> GetTracksFromTree();
  
private:
  static const int maxNtracks = 1000;   ///< Maximum supported number of tracks per event
  int nTracks;                          ///< Number of tracks in the current tree entry
  
  map<string, float[maxNtracks] >  arrayValuesFloat;  ///< Float per-track variables in the current entry
  map<string, int[maxNtracks] >    arrayValuesInt;    ///< Int per-track variables in the current entry
};

#endif /* TrackProcessor_hpp */
