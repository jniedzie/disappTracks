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
  
  // Getters
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
  inline double GetMaxRelativeIsolation(){return maxRelIso;}
  
  inline bool GetRequireSameNpixelHitsLayers(){return sameNpixelHitsLayers;}
  inline bool GetRequireSameNtrackerHitsLayers(){return sameNtrackerHitsLayers;}
  
  inline int GetMinNpixelHits(){return minNpixelHits;}
  inline int GetMaxNpixelHits(){return maxNpixelHits;}
  
  inline int GetMinNpixelLayers(){return minNpixelLayers;}
  inline int GetMaxNpixelLayers(){return maxNpixelLayers;}
  
  inline int GetMinNmissingInnerPixel(){return minMissingInnerPixel;}
  inline int GetMaxNmissingInnerPixel(){return maxMissingInnerPixel;}
  
  inline int GetMinNmissingMiddleTracker(){return minMissingMiddleTracker;}
  inline int GetMaxNmissingMiddleTracker(){return maxMissingMiddleTracker;}
  
  inline int GetMinNmissingOuterTracker(){return minMissingOuterTracker;}
  inline int GetMaxNmissingOuterTracker(){return maxMissingOuterTracker;}
  
  // Setters
  inline void SetNdedxClusters(int min, int max){minDedxClusters=min;maxDedxClusters=max;}
  inline void SetNdets(int min, int max){minDets=min;maxDets=max;}
  inline void SetMinDedxPerCluster(double min){minDedxPerCluster=min;}
  inline void SetTotalDedx(double min, double max){minTotalDeDx=min; maxTotalDeDx=max;}
  inline void SetMinPt(double min){minPt = min;}
  inline void SetMaxEmCalo(double max){maxEmCalo = max;}
  inline void SetMaxHadCalo(double max){maxHadCalo = max;}
  inline void SetMaxEta(double max){maxEta = max;}
  inline void SetMaxRelativeIsolation(double val){maxRelIso = val;}
  
  inline void SetRequireSameNpixelHitsLayers(bool val){sameNpixelHitsLayers=val;}
  inline void SetRequireSameNtrackerHitsLayers(bool val){sameNtrackerHitsLayers=val;}
  
  inline void SetNpixelHits(int min, int max){minNpixelHits=min;maxNpixelHits=max;}
  inline void SetNpixelLayers(int min, int max){minNpixelLayers=min;maxNpixelLayers=max;}
  
  inline void SetNmissingInnerPixel(int min, int max){minMissingInnerPixel=min;maxMissingInnerPixel=max;}
  inline void SetNmissingMiddleTracker(int min, int max){
    minMissingMiddleTracker=min;maxMissingMiddleTracker=max;
  }
  inline void SetNmissingOuterTracker(int min, int max){
    minMissingOuterTracker=min;maxMissingOuterTracker=max;
  }
  
  
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
  double maxRelIso;         ///< max relative isolation in dR=0.3
  
  bool sameNpixelHitsLayers;    ///< require the same number of hits and layers in the pixel
  bool sameNtrackerHitsLayers;  ///< require the same number of hits and layers in the tracker
  
  int minNpixelHits;            ///< min number of pixel hits
  int maxNpixelHits;            ///< max number of pixel hits
  
  int minNpixelLayers;            ///< min number of pixel layers
  int maxNpixelLayers;            ///< max number of pixel layers
  
  int minMissingInnerPixel;     ///< min number of missing inner pixel hits
  int maxMissingInnerPixel;     ///< max number of missing inner pixel hits
  int minMissingMiddleTracker;  ///< min number of missing middle tracker hits
  int maxMissingMiddleTracker;  ///< max number of missing middle tracker hits
  int minMissingOuterTracker;   ///< min number of missing outer tracker hits
  int maxMissingOuterTracker;   ///< max number of missing outer tracker hits
};

#endif /* TrackCut_hpp */
