//
//  TrackCut.hpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#ifndef TrackCut_hpp
#define TrackCut_hpp

class TrackCut {
public:
  enum ECut {
    kEmpty    = 1,
    kPt200GeV = 1 << 1, ///< pT ≥ 200 GeV
    kLowCalo  = 1 << 2, ///< energy deposit in both EM and Had calo ≤ 0.5 GeV
    kLowDEdx  = 1 << 3, ///< ΣdE/dx ≤ 38 MeV
    kShort    = 1 << 4, ///< 3-4 clusters
    kMedium   = 1 << 5, ///< 3-8 clusters
    kPt50GeV  = 1 << 6, ///< pT ≥ 50 GeV
    kEta2p4   = 1 << 7, ///< |eta| ≤ 2.4
  };
  
  TrackCut(ECut cutType=kEmpty);
  ~TrackCut();
  
  inline int GetMinDedxClusters(){return minDedxClusters;}
  inline int GetMaxDedxClusters(){return maxDedxClusters;}
  inline int GetMinDets(){return minDets;}
  inline int GetMaxDets(){return maxDets;}
  inline double GetMinDedxPerCluster(){return minDedxPerCluster;}
  inline double GetMinTotalDedx(){return minTotalDeDx;}
  inline double GetMaxTotalDedx(){return maxTotalDeDx;}
  inline double GetMinPt(){return minPt;}
  inline double GetMaxEmCalo(){return maxEmCalo;}
  inline double GetMaxHadCalo(){return maxHadCalo;}
  inline double GetMaxEta(){return maxEta;}
  
  inline void SetNdedxClusters(int min, int max){minDedxClusters=min;maxDedxClusters=max;}
  inline void SetNdets(int min, int max){minDets=min;maxDets=max;}
  inline void SetMinDedxPerCluster(double min){minDedxPerCluster=min;}
  inline void SetTotalDedx(double min, double max){minTotalDeDx=min; maxTotalDeDx=max;}
  inline void SetMinPt(double min){minPt = min;}
  inline void SetMaxEmCalo(double max){maxEmCalo = max;}
  inline void SetMaxHadCalo(double max){maxHadCalo = max;}
  inline void SetMaxEta(double max){maxEta = max;}
  
private:
  int minDedxClusters;      ///< min number of dedx clusters along the track
  int maxDedxClusters;      ///< max number of dedx clusters along the track
  int minDets;              ///< min number of subdetectors
  int maxDets;              ///< max number of subdetectors
  double minDedxPerCluster; ///< min dedx at each track's cluster
  double minTotalDeDx;      ///< min total dedx along the track
  double maxTotalDeDx;      ///< max total dedx along the track
  double minPt;             ///< min transverse momentum of the track
  double maxEmCalo;         ///< max energy deposit in EM calorimeter
  double maxHadCalo;        ///< max energy deposit in hadronic calorimeter
  double maxEta;            ///< maximum pseudorapidity
};

#endif /* TrackCut_hpp */
