//
//  Track.hpp
//
//  Created by Jeremi Niedziela on 16/07/2018.
//

#ifndef Track_hpp
#define Track_hpp

#include "Helpers.hpp"
#include "TrackCut.hpp"

#include <vector>

class Track{
public:
  Track();
  ~Track(){};
  
  void CalculateInternals(); // Call this method when everything was already added to a track
  
  // Setters
  void SetDeDxInLayer(int layer, float value){dedx[layer] = value;}
  void SetSubDetIdInLayer(int layer, int id){subDetId[layer] = id;}
  void SetSizeXinLayer(int layer, int size){sizeX[layer] = size;}
  void SetSizeYinLayer(int layer, int size){sizeY[layer] = size;}
  void SetEta(double _eta){eta=_eta;}
  void SetPhi(double _phi){phi=_phi;}
  void SetCaloEmEnergy(double _energy){caloEmEnergy = _energy;}
  void SetCaloHadEnergy(double _energy){caloHadEnergy = _energy;}
  void SetDxy(double _d,double _dErr){dxy=_d;dxyErr=_dErr;}
  void SetDz(double _d,double _dErr){dz=_d;dzErr=_dErr;}
  void SetCharge(int _charge){charge = _charge;}
  void SetMass(double _mass){mass = _mass;}
  void SetPt(double _pt){pt = _pt;}
  void SetPid(int _pid){pid = _pid;}
  void SetRelativeIsolation(double val){relIso03 = val;}
  
  void SetNtrackerLayers(int n){nTrackerLayers = n;}
  void SetNpixelLayers(int n){nPixelLayers = n;}
  void SetNtrackerHits(int n){nTrackerHits = n;}
  void SetNpixelHits(int n){nPixelHits = n;}
  void SetNmissingInnerPixelHits(int n){nMissingInnerPixelHits = n;}
  void SetNmissingOuterPixelHits(int n){nMissingOuterPixelHits = n;}
  void SetNmissingInnerStripHits(int n){nMissingInnerStripHits = n;}
  void SetNmissingOuterStripHits(int n){nMissingOuterStripHits = n;}
  void SetNmissingInnerTrackerHits(int n){nMissingInnerTrackerHits = n;}
  void SetNmissingOuterTrackerHits(int n){nMissingOuterTrackerHits = n;}
  void SetNmissingMiddleTrackerHits(int n){nMissingMiddleTrackerHits = n;}
  
  // Getters
  float   GetDeDxInLayer(int layer){return dedx[layer];}
  int     GetSubDetIdInLayer(int layer){return subDetId[layer];}
  int     GetSizeXinLayer(int layer){return sizeX[layer];}
  int     GetSizeYinLayer(int layer){return sizeY[layer];}
  double  GetEta(){return eta;}
  double  GetPhi(){return phi;}
  
  float   GetTotalDedx(){return accumulate(dedx.begin(),dedx.end(),0.0);}
  double  GetCaloEmEnergy(){return caloEmEnergy;}
  double  GetCaloHadEnergy(){return caloHadEnergy;}
  double  GetDxy(){return dxy;}
  double  GetDxyErr(){return dxyErr;}
  double  GetDz(){return dz;}
  double  GetDzErr(){return dzErr;}
  int     GetCharge(){return charge;}
  double  GetMass(){return mass;}
  double  GetPt(){return pt;}
  int     GetPid(){return pid;}
  double  GetRelativeIsolation(){return relIso03;}
  double  GetAbsoluteIsolation(){return relIso03*pt;}
  double  GetDedxInSubDet(int det);
  
  int GetNtrackerLayers(){return nTrackerLayers;}
  int GetNpixelLayers(){return nPixelLayers;}
  int GetNtrackerHits(){return nTrackerHits;}
  int GetNpixelHits(){return nPixelHits;}
  int GetNmissingInnerPixelHits(){return nMissingInnerPixelHits;}
  int GetNmissingOuterPixelHits(){return nMissingOuterPixelHits;}
  int GetNmissingInnerStripHits(){return nMissingInnerStripHits;}
  int GetNmissingOuterStripHits(){return nMissingOuterStripHits;}
  int GetNmissingInnerTrackerHits(){return nMissingInnerTrackerHits;}
  int GetNmissingOuterTrackerHits(){return nMissingOuterTrackerHits;}
  int GetNmissingMiddleTrackerHits(){return nMissingMiddleTrackerHits;}
  
  int GetNclusters(){return nClusters;}
  int GetNdetIDs(){return nDetIDs;}
  
  // Other methods
  bool IsPassingCut(TrackCut *cut);
  void Print();
private:
  std::vector<float> dedx;   // dedx in consecutive layers
  std::vector<int> subDetId; // sub-detector IDs for each layer
  std::vector<int> sizeX;    // cluster size X in each layer
  std::vector<int> sizeY;    // cluster size Y in each layer
  double eta;
  double phi;
  double caloEmEnergy;
  double caloHadEnergy;
  double dxy;
  double dxyErr;
  double dz;
  double dzErr;
  int charge;
  double mass;
  double pt;
  int pid;
  double relIso03;  ////< Relative track isolation in cone dR=0.3
  
  int nTrackerLayers;           // Number of tracker layers
  int nPixelLayers;             // Number of pixel layers
  int nTrackerHits;             // Number of tracker hits
  int nPixelHits;               // Number of pixel hits
  int nMissingInnerPixelHits;   // Number of missing inner pixel hits
  int nMissingOuterPixelHits;   // Number of missing outer pixel hits
  int nMissingInnerStripHits;   // Number of missing inner strip hits
  int nMissingOuterStripHits;   // Number of missing outer strip hits
  int nMissingInnerTrackerHits; // Number of missing inner tracker hits
  int nMissingOuterTrackerHits; // Number of missing outer tracker hits
  int nMissingMiddleTrackerHits;// Number of missing middle tracker hits
  
  int nDetIDs;  // total number of sub-detectors hit by this track
  int nClusters;// total number of clusters belonging to this track
};

#endif /* Track_hpp */
