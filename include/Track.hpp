//
//  Track.hpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#ifndef Track_hpp
#define Track_hpp

#include "Helpers.hpp"
#include "TrackCut.hpp"

/// This class represents a single track. It contains information about its kinematical
/// properties, number of hits in different parts of the detector, energy loss in the tracker etc.
/// It also checks if this track passes some set of selection criteria represented by TrackCut object.
class Track{
public:
  /// Default constructor
  Track();
  
  /// Default destructor
  ~Track(){};
  
  /// Print basic information about the track
  void Print();
  
  /// Check if track passes selection criteria
  /// \param cut Tracks selection criteria to be checked
  bool IsPassingCut(TrackCut *cut);
  
  
  // Setters
  void SetDeDxInLayer(int layer, float value);
  void SetSubDetIdInLayer(int layer, int id);
  inline void SetSizeXinLayer(int layer, int size){sizeX[layer] = size;}
  inline void SetSizeYinLayer(int layer, int size){sizeY[layer] = size;}
  
  inline void SetPt(double val){pt = val;}
  inline void SetEta(double val){eta=val;}
  inline void SetPhi(double val){phi=val;}
  inline void SetMass(double val){mass = val;}
  
  inline void SetCaloEmEnergy(double val) {caloEmEnergy = val;}
  inline void SetCaloHadEnergy(double val){caloHadEnergy = val;}
  inline void SetRelativeIsolation(double val){relativeIsolation = val;}
  inline void SetDxy(double val,double err) {dxy=val;dxyErr=err;}
  inline void SetDz(double val,double err)  {dz =val;dzErr =err;}
  
  inline void SetCharge(int val){charge = val;}
  inline void SetPid(int val){pid = val;}
  
  inline void SetNtrackerLayers(int n){nTrackerLayers = n;}
  inline void SetNpixelLayers(int n){nPixelLayers = n;}
  inline void SetNtrackerHits(int n){nTrackerHits = n;}
  inline void SetNpixelHits(int n){nPixelHits = n;}
  inline void SetNmissingInnerPixelHits(int n){nMissingInnerPixelHits = n;}
  inline void SetNmissingOuterPixelHits(int n){nMissingOuterPixelHits = n;}
  inline void SetNmissingInnerStripHits(int n){nMissingInnerStripHits = n;}
  inline void SetNmissingOuterStripHits(int n){nMissingOuterStripHits = n;}
  inline void SetNmissingInnerTrackerHits(int n){nMissingInnerTrackerHits = n;}
  inline void SetNmissingOuterTrackerHits(int n){nMissingOuterTrackerHits = n;}
  inline void SetNmissingMiddleTrackerHits(int n){nMissingMiddleTrackerHits = n;}
  
  // Getters
  inline double  GetDeDxInLayer(int layer){return dedx[layer];}
  inline double  GetTotalDedx(){return accumulate(dedx.begin(),dedx.end(),0.0);}
  double  GetDedxInSubDet(int det);
  
  inline int     GetSubDetIdInLayer(int layer){return subDetId[layer];}
  inline int     GetSizeXinLayer(int layer){return sizeX[layer];}
  inline int     GetSizeYinLayer(int layer){return sizeY[layer];}
  
  inline double  GetPt(){return pt;}
  inline double  GetEta(){return eta;}
  inline double  GetPhi(){return phi;}
  inline double  GetMass(){return mass;}
  
  inline double  GetCaloEmEnergy(){return caloEmEnergy;}
  inline double  GetCaloHadEnergy(){return caloHadEnergy;}
  inline double  GetRelativeIsolation(){return relativeIsolation;}
  inline double  GetAbsoluteIsolation(){return relativeIsolation*pt;}
  
  inline double  GetDxy(){return dxy;}
  inline double  GetDxyErr(){return dxyErr;}
  inline double  GetDz(){return dz;}
  inline double  GetDzErr(){return dzErr;}
  
  inline int     GetCharge(){return charge;}
  inline int     GetPid(){return pid;}
  
  inline int GetNtrackerLayers(){return nTrackerLayers;}
  inline int GetNpixelLayers(){return nPixelLayers;}
  inline int GetNtrackerHits(){return nTrackerHits;}
  inline int GetNpixelHits(){return nPixelHits;}
  inline int GetNmissingInnerPixelHits(){return nMissingInnerPixelHits;}
  inline int GetNmissingOuterPixelHits(){return nMissingOuterPixelHits;}
  inline int GetNmissingInnerStripHits(){return nMissingInnerStripHits;}
  inline int GetNmissingOuterStripHits(){return nMissingOuterStripHits;}
  inline int GetNmissingInnerTrackerHits(){return nMissingInnerTrackerHits;}
  inline int GetNmissingOuterTrackerHits(){return nMissingOuterTrackerHits;}
  inline int GetNmissingMiddleTrackerHits(){return nMissingMiddleTrackerHits;}
  
  inline int GetNdetIDs(){return nDetIDs;}
  inline int GetNclusters(){return nClusters;}
  
private:
  std::vector<float> dedx;    ///< dE/dx in consecutive layers
  std::vector<int> subDetId;  ///< Sub-detector IDs for each layer
  std::vector<int> sizeX;     ///< Cluster size X in each layer
  std::vector<int> sizeY;     ///< Cluster size Y in each layer
  double pt;                  ///< Transverse momentum (GeV)
  double eta;                 ///< Pseudorapidity
  double phi;                 ///< Polar angle
  double mass;                ///< Particle's mass (GeV)
  double caloEmEnergy;        ///< Energy deposit in EM calorimeter (GeV)
  double caloHadEnergy;       ///< Energy deposit in Hadronic calorimeter (GeV)
  double relativeIsolation;   ///< Relative track isolation in cone dR=0.3
  double dxy;                 ///< Distance from primary vertex in XY
  double dxyErr;              ///< Uncertainty on the distance from primary vertex in XY
  double dz;                  ///< Distance from primary vertex in Z
  double dzErr;               ///< Uncertainty on the distance from primary vertex in Z
  int charge;                 ///< Particle's charge
  int pid;                    ///< Particle's PDG PID code
  
  
  int nTrackerLayers;             ///< Number of tracker layers
  int nPixelLayers;               ///< Number of pixel layers
  int nTrackerHits;               ///< Number of tracker hits
  int nPixelHits;                 ///< Number of pixel hits
  int nMissingInnerPixelHits;     ///< Number of missing inner pixel hits
  int nMissingOuterPixelHits;     ///< Number of missing outer pixel hits
  int nMissingInnerStripHits;     ///< Number of missing inner strip hits
  int nMissingOuterStripHits;     ///< Number of missing outer strip hits
  int nMissingInnerTrackerHits;   ///< Number of missing inner tracker hits
  int nMissingOuterTrackerHits;   ///< Number of missing outer tracker hits
  int nMissingMiddleTrackerHits;  ///< Number of missing middle tracker hits
  
  
  int nDetIDs;    ///< Total number of sub-detectors hit by this track
  int nClusters;  ///< Total number of clusters belonging to this track
  
  /// Calculates additional properties that are not read directly from ntuples.
  /// Should be called after adding setting some of the track's properties
  void CalculateInternals();
};

#endif /* Track_hpp */
