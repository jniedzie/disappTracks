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
  unique_ptr<Track> GetRandomTrack(int nLayers, double maxEta);
  
  /// Check if track passes selection criteria
  /// \param track Tracks to which cuts should be applied
  /// \param cut Tracks selection criteria to be checked
  bool IsPassingCut(const shared_ptr<Track> track,
                    const unique_ptr<TrackCut> &cut);
  
private:
  
};

#endif /* TrackProcessor_hpp */
