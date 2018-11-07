//
//  TrackCut.hpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#ifndef TrackCut_hpp
#define TrackCut_hpp

#include "Helpers.hpp"

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
  
  
  // Setters
  inline void SetPt(range<double> val){pt = val;}
  inline void SetEta(range<double> val){eta = val;}
  inline void SetCaloEmEnergy(range<double> val){caloEmEnergy = val;}
  inline void SetCaloHadEnergy(range<double> val){caloHadEnergy = val;}
  inline void SetRelativeIsolation(range<double> val){relativeIsolation = val;}
  
  inline void SetRequireSameNpixelHitsLayers(bool val){sameNpixelHitsLayers=val;}
  inline void SetRequireSameNtrackerHitsLayers(bool val){sameNtrackerHitsLayers=val;}
  
  inline void SetNpixelLayers(range<int> val){nPixelLayers=val;}
  inline void SetNpixelHits(range<int> val){nPixelHits=val;}
  inline void SetNmissingInnerPixel(range<int> val){nMissingInnerPixel=val;}
  inline void SetNmissingOuterTracker(range<int> val){nMissingOuterTracker=val;}
  inline void SetNmissingMiddleTracker(range<int> val){nMissingMiddleTracker=val;}
  
  inline void SetDedxPerCluster(range<double> val){dedxPerCluster=val;}
  inline void SetTotalDedx(range<double> val){totalDeDx=val;}
  inline void SetNdetIDs(range<int> val){nDetIDs=val;}
  inline void SetNdedxClusters(range<int> val){nDedxClusters = val;}
  
  // Getters
  inline range<double> GetPt(){return pt;}
  inline range<double> GetEta(){return eta;}
  inline range<double> GetCaloEmEnergy(){return caloEmEnergy;}
  inline range<double> GetCaloHadEnergy(){return caloHadEnergy;}
  inline range<double> GetRelativeIsolation(){return relativeIsolation;}
  
  inline bool GetRequireSameNpixelHitsLayers(){return sameNpixelHitsLayers;}
  inline bool GetRequireSameNtrackerHitsLayers(){return sameNtrackerHitsLayers;}
  
  inline range<int> GetNpixelLayers(){return nPixelLayers;}
  inline range<int> GetNpixelHits(){return nPixelHits;}
  inline range<int> GetNmissingInnerPixel(){return nMissingInnerPixel;}
  inline range<int> GetNmissingOuterTracker(){return nMissingOuterTracker;}
  inline range<int> GetNmissingMiddleTracker(){return nMissingMiddleTracker;}
  
  inline range<double>  GetDedxPerCluster(){return dedxPerCluster;}
  inline range<double>  GetTotalDedx(){return totalDeDx;}
  inline range<int>     GetNdetIDs(){return nDetIDs;}
  inline range<int>     GetNdedxClusters(){return nDedxClusters;}
  
private:
  range<double> pt;                 ///< allowed transverse momentum of the track
  range<double> eta;                ///< allowed pseudorapidity
  range<double> caloEmEnergy;       ///< allowed energy deposit in EM calorimeter
  range<double> caloHadEnergy;      ///< allowed energy deposit in Hadronic calorimeter
  range<double> relativeIsolation;  ///< max relative isolation in dR=0.3
  
  bool sameNpixelHitsLayers;        ///< require the same number of hits and layers in the pixel
  bool sameNtrackerHitsLayers;      ///< require the same number of hits and layers in the tracker
  
  range<int> nPixelLayers;          ///< allowed number of pixel layers
  range<int> nPixelHits;            ///< allowed number of pixel hits
  range<int> nMissingInnerPixel;    ///< allowed number of missing inner pixel hits
  range<int> nMissingOuterTracker;  ///< allowed number of missing outer tracker hits
  range<int> nMissingMiddleTracker; ///< allowed number of missing middle tracker hits
  
  range<double> dedxPerCluster;     ///< allowed dedx at each track's cluster
  range<double> totalDeDx;          ///< allowed total dedx along the track
  range<int> nDetIDs;               ///< allowed number of subdetectors
  range<int> nDedxClusters;         ///< allowed number of dedx clusters along the track
};

#endif /* TrackCut_hpp */
