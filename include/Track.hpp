//
//  Track.hpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 16/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
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
  void SetIsShort(bool _isShort){isShort = _isShort;}
  
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
  
  bool    GetIsShort(){return isShort;}
  int     GetNclusters();
  
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
  
  bool isShort; // track is short if it has max 3 dedx points
};

#endif /* Track_hpp */
