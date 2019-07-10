//  Track.hpp
//
//  Created by Jeremi Niedziela on 16/07/2018.

#ifndef Track_hpp
#define Track_hpp

#include "Helpers.hpp"
#include "TrackCut.hpp"
#include "Point.hpp"

/// This class represents a single track. It contains information about its kinematical
/// properties, number of hits in different parts of the detector, energy loss in the tracker etc.
/// It also checks if this track passes some set of selection criteria represented by TrackCut object.
class Track{
public:
  /// Default constructor
  Track();
  
  /// Simple track constructor
  Track(double _eta, double _phi, int _charge, int _nTrackerLayers, double _pt);
  
  /// Copy constructor
  Track(const Track &t);
  
  /// Default destructor
  ~Track(){};
  
  /// Print basic information about the track
  void Print();
  
  // Setters
  inline void SetEventMetPt(double val){eventMetPt = val;}
  inline void SetEventMetEta(double val){eventMetEta = val;}
  inline void SetEventMetPhi(double val){eventMetPhi = val;}
  inline void SetEventMetMass(double val){eventMetMass = val;}
  
  inline void SetDecayPoint(const Point &val){decayPoint = val;}
  
  // Getters
  int GetLastBarrelLayer();
  
  inline int     GetNdEdxHits(){return (int)dedx.size();}
  inline int     GetNnotEmptyDedxHits(){
    int maxHit=0;
    for(int iHit=0; iHit<dedx.size(); iHit++){
      if(dedx[iHit] > 0) maxHit = iHit;
    }
    return maxHit+1;
  }
	inline double  GetDeDxForHit(int iHit){return iHit < dedx.size() ? dedx[iHit] : 0;}
  inline double  GetTotalDedx(){return accumulate(dedx.begin(),dedx.end(),0.0);}
  double         GetAverageDedx();
  double         GetDedxInSubDet(int det);
  double         GetDedxInBarrelLayer(int iLayer);
  
  inline int     GetSubDetIdForHit(int iHit){return iHit < subDetId.size() ? subDetId[iHit] : 0;}
	inline int     GetSizeXforHit(int iHit){return iHit < sizeX.size() ? sizeX[iHit] : 0;}
	inline int     GetSizeYforHit(int iHit){return iHit < sizeY.size() ? sizeY[iHit] : 0;}
	inline int     GetDetTypeForHit(int iHit){return iHit < detType.size() ? detType[iHit] : 0;}
	inline int     GetLayerForHit(int iHit){return iHit < layer.size() ? layer[iHit] : 0;}
	inline int     GetLadderForHit(int iHit){return iHit < ladder.size() ? ladder[iHit] : 0;}
  
  inline double  GetPt(){return pt;}
  inline double  GetEta(){return eta;}
  inline double  GetTheta() const {return 2*atan(exp(-eta));}
  inline double  GetPhi() const {return phi;}
  inline double  GetMass(){return mass;}
  
  inline double  GetCaloEmEnergy(){return caloEmEnergy;}
  inline double  GetCaloHadEnergy(){return caloHadEnergy;}
  inline double  GetRelativeIsolation(){return relativeIsolation;}
  inline double  GetAbsoluteIsolation(){return relativeIsolation*pt;}
  
  inline double  GetDxy(){return dxy;}
  inline double  GetDxyErr(){return dxyErr;}
  inline double  GetDz(){return dz;}
  inline double  GetDzErr(){return dzErr;}
  
  inline int     GetCharge() const {return charge;}
  inline int     GetPid(){return pid;}
  inline int     GetMcMatch(){return mcMatch;}
  
  inline int GetNtrackerLayers() const {return nTrackerLayers;}
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
  inline int GetNdedxClusters(){return nDedxClusters;}
  
  inline Point& GetDecayPoint(){return decayPoint;}
  
private:
  vector<float> dedx;         ///< dE/dx in consecutive hits
  vector<int> subDetId;       ///< Sub-detector IDs for each hit
  vector<int> sizeX;          ///< Cluster size X of each hit
  vector<int> sizeY;          ///< Cluster size Y of each hit
  vector<int> detType;        ///< Type of detector from which the hit comes (0 = strips, 1 = bpix, 2 = fpix)
  vector<int> layer;          ///< Layer for given dE/dx hit
  vector<int> ladder;         ///< Ladder for gien dE/dx hit
  
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
  int mcMatch;                ///< Does the reconstructed particle match an mc generated one
  
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
  int nDedxClusters;  ///< Total number of clusters belonging to this track
  
  double eventMetPt;   ///< MET transverse momentum of the event that contains this track
  double eventMetEta;  ///< MET pseudorapidity of the event that contains this track
  double eventMetPhi;  ///< MET polar angle of the event that contains this track
  double eventMetMass; ///< MET mass of the event that contains this track
  
  /// Calculates additional properties that are not read directly from ntuples.
  /// Should be called after adding setting some of the track's properties
  void CalculateInternals();
  
  Point decayPoint;
  
  friend class TrackProcessor;
};

#endif /* Track_hpp */
