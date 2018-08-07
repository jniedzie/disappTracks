//
//  Lepton.hpp
//  xDisappLeptons
//
//  Created by Jeremi Niedziela on 07/08/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#ifndef Lepton_hpp
#define Lepton_hpp


#include "Helpers.hpp"
//#include "LeptonCut.hpp"

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
  
  // Getters
  double  GetEta(){return eta;}
  double  GetPhi(){return phi;}
  double  GetPt(){return pt;}
  int     GetPid(){return pid;}
  
  // Other methods
  bool IsPassingCut(/*LeptonCut *cut*/);
  
private:
  double eta;
  double phi;
  double pt;
  int pid;
};


#endif /* Lepton_hpp */
