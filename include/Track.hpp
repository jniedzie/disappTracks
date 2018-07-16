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

#include <vector>

class Track{
public:
  Track();
  ~Track(){};
  
  void SetDeDxInLayer(int layer, float value){dedx[layer] = value;}
  void SetSubDetIdInLayer(int layer, int id){subDetId[layer] = id;}
  void SetEta(double _eta){eta=_eta;}
  void SetPhi(double _phi){phi=_phi;}
  void SetIsShort(bool _isShort){isShort = _isShort;}
  
  float GetDeDxInLayer(int layer){return dedx[layer];}
  int GetSubDetIdInLayer(int layer){return subDetId[layer];}
  double GetEta(){return eta;}
  double GetPhi(){return phi;}
  bool GetIsShort(){return isShort;}
  float GetTotalDedx(){return accumulate(dedx.begin(),dedx.end(),0.0);}
  
  int GetNclusters();
  
  void Print();
private:
  std::vector<float> dedx;   // dedx in consecutive layers
  std::vector<int> subDetId; // sub-detector IDs for each layer
  double eta;
  double phi;
  bool isShort; // track is short if it has max 3 dedx points
};

#endif /* Track_hpp */
