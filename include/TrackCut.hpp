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
    kEmpty    = 1,
    kHighPt   = 1 << 1, ///< pT ≥ 200 GeV
    kLowCalo  = 1 << 2, ///< energy deposit in both EM and Had calo ≤ 0.5 GeV
    kLowDEdx  = 1 << 3, ///< ΣdE/dx ≤ 38 MeV
    kShort    = 1 << 4, ///< 3-4 clusters
    kMedium   = 1 << 5, ///< 3-8 clusters
  };
  
  TrackCut(ECut cutType=kEmpty);
  ~TrackCut();
  
  inline int GetMinDedxClusters(){return minDedxClusters;}
  inline int GetMaxDedxClusters(){return maxDedxClusters;}
  inline double GetMinDedxPerCluster(){return minDedxPerCluster;}
  inline double GetMinTotalDedx(){return minTotalDeDx;}
  inline double GetMaxTotalDedx(){return maxTotalDeDx;}
  inline double GetMinPt(){return minPt;}
  inline double GetMaxEmCalo(){return maxEmCalo;}
  inline double GetMaxHadCalo(){return maxHadCalo;}
  
  inline void SetNdedxClusters(int min, int max){minDedxClusters=min;maxDedxClusters=max;}
  inline void SetMinDedxPerCluster(double min){minDedxPerCluster=min;}
  inline void SetTotalDedx(double min, double max){minTotalDeDx=min; maxTotalDeDx=max;}
  inline void SetMinPt(double min){minPt = min;}
  inline void SetMaxEmCalo(double max){maxEmCalo = max;}
  inline void SetMaxHadCalo(double max){maxHadCalo = max;}
  
private:
  int minDedxClusters;      ///< min number of dedx clusters along the track
  int maxDedxClusters;      ///< max number of dedx clusters along the track
  double minDedxPerCluster; ///< min dedx at each track's cluster
  double minTotalDeDx;      ///< min total dedx along the track
  double maxTotalDeDx;      ///< max total dedx along the track
  double minPt;             ///< min transverse momentum of the track
  double maxEmCalo;         ///< max energy deposit in EM calorimeter
  double maxHadCalo;        ///< max energy deposit in hadronic calorimeter
};

#endif /* TrackCut_hpp */
