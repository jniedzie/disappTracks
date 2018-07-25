//
//  TrackCut.hpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#ifndef TrackCut_hpp
#define TrackCut_hpp

class TrackCut {
public:
  
  
  enum ECut {
    kEmpty,
    kShort,                 ///< 3 to 4 dE/dx clusters
    kShortAboveThreshold,   ///< 3 to 4 dE/dx clusters, each cluster ≥ 2.5 MeV
    kShortLowTotalDEdx      ///< max 10 dE/dx clusters, ΣdE/dx ≥ 38 MeV
  };
  
  TrackCut(ECut cutType=kEmpty);
  ~TrackCut();
  
  inline int GetMinDedxClusters(){return minDedxClusters;}
  inline int GetMaxDedxClusters(){return maxDedxClusters;}
  inline double GetMinDedxPerCluster(){return minDedxPerCluster;}
  inline double GetMinTotalDedx(){return minTotalDeDx;}
  inline double GetMaxTotalDedx(){return maxTotalDeDx;}
  
  inline void SetNdedxClusters(int min, int max){minDedxClusters=min;maxDedxClusters=max;}
  inline void SetMinDedxPerCluster(double min){minDedxPerCluster=min;}
  inline void SetTotalDedx(double min, double max){minTotalDeDx=min; maxTotalDeDx=max;}
  
private:
  int minDedxClusters;      ///< min number of dedx clusters along the track
  int maxDedxClusters;      ///< max number of dedx clusters along the track
  double minDedxPerCluster; ///< min dedx at each track's cluster
  double minTotalDeDx;      ///< min total dedx along the track
  double maxTotalDeDx;      ///< max total dedx along the track
  
};

#endif /* TrackCut_hpp */
