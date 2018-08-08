//
//  Lepton.hpp
//
//  Created by Jeremi Niedziela on 07/08/2018.
//

#ifndef Lepton_hpp
#define Lepton_hpp


#include "Helpers.hpp"
#include "LeptonCut.hpp"

#include <vector>

class Lepton{
public:
  Lepton();
  ~Lepton(){};
  
  // Setters
  void SetEta(double _eta){eta=_eta;}
  void SetPhi(double _phi){phi=_phi;}
  void SetPt(double _pt){pt = _pt;}
  void SetPid(int _pid){pid = _pid;}
  void SetIsolation(double _isolation){isolation = _isolation;}
  void SetTightID(bool id){tightID = id;}
  
  // Getters
  double  GetEta(){return eta;}
  double  GetPhi(){return phi;}
  double  GetPt(){return pt;}
  int     GetPid(){return pid;}
  double  GetIsolation(){return isolation;}
  bool    GetTightID(){return tightID;}
  
  // Other methods
  bool IsPassingCut(LeptonCut *cut);
  
private:
  double eta;
  double phi;
  double pt;
  bool tightID;
  int pid;
  double isolation;
};


#endif /* Lepton_hpp */
