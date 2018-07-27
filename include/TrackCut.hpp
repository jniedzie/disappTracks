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
    kHighPt,                ///< pT ≥ 200 GeV
    kLowCalo,               ///< energy deposit in both EM and Had calo ≤ 0.5 GeV
    kLowTotal,              ///< ΣdE/dx ≤ 38 MeV
    kLowCaloLowDEdx,        ///< calo ≤ 0.5 GeV, ΣdE/dx ≤ 38 MeV
    kShort,                 ///< 3-4 clusters
    kShortLowCalo,          ///< 3-4 clusters, calo ≤ 0.5 GeV
    kShortHighPt,           ///< 3-4 clusters, pT ≥ 200 GeV
    kShortLowTotal,         ///< 3-4 clusters, ΣdE/dx ≤ 38 MeV
    kShortLowTotalHighPt,   ///< 3-4 clusters, ΣdE/dx ≤ 38 MeV, pT ≥ 200 GeV
    kShortAboveThreshold,   ///< 3 to 4 dE/dx clusters, each cluster ≥ 2.5 MeV
    kMedium,                ///< 3-8 clusters
    kMediumLowCalo,         ///< 3-8 clusters, calo ≤ 0.5 GeV
    kMediumHighPt,          ///< 3-8 clusters, pT ≥ 200 GeV
    kMediumLowTotal,        ///< 3-8 clusters, ΣdE/dx ≤ 38 MeV
    kMediumLowTotalHighPt,  ///< 3-8 clusters, ΣdE/dx ≤ 38 MeV, pT ≥ 200 GeV
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
