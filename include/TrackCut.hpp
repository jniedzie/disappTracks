//
//  TrackCut.hpp
//  disappTracksTarget
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
    kShort,
    kShortAboveThreshold
  };
  
  TrackCut(ECut cutType=kEmpty);
  ~TrackCut();
  
  int GetMinDedxClusters(){return minDedxClusters;}
  int GetMaxDedxClusters(){return minDedxClusters;}
  double GetMinDedxPerCluster(){return minDedxPerCluster;}
  double GetMinTotalDedx(){return minTotalDeDx;}
  
  void SetNdedxClusters(int min, int max){minDedxClusters=min;maxDedxClusters=max;}
  void SetMinDedxPerCluster(double min){minDedxPerCluster=min;}
  void SetMinTotalDedx(double min){minTotalDeDx=min;}
  
private:
  int minDedxClusters;      ///< min number of dedx clusters along the track
  int maxDedxClusters;      ///< max number of dedx clusters along the track
  double minDedxPerCluster; ///< min dedx at each track's cluster
  double minTotalDeDx;      ///< min total dedx along the track
  
};

#endif /* TrackCut_hpp */
