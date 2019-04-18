//  TrackCut.hpp
//
//  Created by Jeremi Niedziela on 16/07/2018.

#ifndef TrackCut_hpp
#define TrackCut_hpp

#include "Helpers.hpp"

/// Class containing definition of the track selection criteria.
/// User can define ranges of allowed parameters and required flags.
class TrackCut {
public:
  /// Default constructor
  TrackCut();
  
  /// Copy constructor
  TrackCut(const TrackCut &c);
  
  /// Default desctructor
  ~TrackCut();
  
  // Setters
  inline void SetPt(range<double> val){pt = val;}
  inline void SetEta(range<double> val){eta = val;}
  inline void SetCaloEmEnergy(range<double> val){caloEmEnergy = val;}
  inline void SetCaloHadEnergy(range<double> val){caloHadEnergy = val;}
  inline void SetRelativeIsolation(range<double> val){relativeIsolation = val;}
  
  inline void SetRequireMcMatch(bool val){requiresMcMatch=val;}
  
  inline void SetNlayers(range<int> val){nLayers=val;}
  inline void SetNpixelLayers(range<int> val){nPixelLayers=val;}
  inline void SetNpixelHits(range<int> val){nPixelHits=val;}
  inline void SetNmissingInnerPixel(range<int> val){nMissingInnerPixel=val;}
  inline void SetNmissingOuterTracker(range<int> val){nMissingOuterTracker=val;}
  inline void SetNmissingMiddleTracker(range<int> val){nMissingMiddleTracker=val;}
  
  inline void SetDedxPerCluster(range<double> val){dedxPerCluster=val;}
  inline void SetTotalDedx(range<double> val){totalDeDx=val;}
  inline void SetNdetIDs(range<int> val){nDetIDs=val;}
  inline void SetNdedxClusters(range<int> val){nDedxClusters = val;}
  
  inline void SetTrackMetDeltaPhi(range<double> val){trackMetDeltaPhi = val;}

private:
  range<double> pt;                 ///< allowed transverse momentum of the track
  range<double> eta;                ///< allowed pseudorapidity
  range<double> caloEmEnergy;       ///< allowed energy deposit in EM calorimeter
  range<double> caloHadEnergy;      ///< allowed energy deposit in Hadronic calorimeter
  range<double> relativeIsolation;  ///< allowed relative isolation in dR=0.3
  range<double> trackMetDeltaPhi;   ///< allowed angle between track and MET
  
  bool requiresMcMatch;              ///< require mcMatch flag to be different than 0
  
  range<int> nLayers;               ///< allowed number of tracker layers
  range<int> nPixelLayers;          ///< allowed number of pixel layers
  range<int> nPixelHits;            ///< allowed number of pixel hits
  range<int> nMissingInnerPixel;    ///< allowed number of missing inner pixel hits
  range<int> nMissingOuterTracker;  ///< allowed number of missing outer tracker hits
  range<int> nMissingMiddleTracker; ///< allowed number of missing middle tracker hits
  
  range<double> dedxPerCluster;     ///< allowed dedx at each track's cluster
  range<double> totalDeDx;          ///< allowed total dedx along the track
  range<int> nDetIDs;               ///< allowed number of subdetectors
  range<int> nDedxClusters;         ///< allowed number of dedx clusters along the track
  
  friend class TrackProcessor;
};

#endif /* TrackCut_hpp */
