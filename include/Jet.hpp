//
//  Jet.hpp
//  disappTracksTarget
//
//  Created by Jeremi Niedziela on 17/07/2018.
//  Copyright © 2018 Jeremi Niedziela. All rights reserved.
//

#ifndef Jet_hpp
#define Jet_hpp

class Jet{
public:
  Jet();
  ~Jet();

  inline void SetPt(double _pt){pt = _pt;}
  inline void SetEta(double _eta){eta = _eta;}
  inline void SetPhi(double _phi){phi = _phi;}
  
  inline double GetPt(){return pt;}
  inline double GetEta(){return eta;}
  inline double GetPhi(){return phi;}
  
  bool IsPassingCut(/*JetCut *cut*/);
  void Print();
  
private:
  double pt;     ///< Transverse momentum
  double eta;    ///< Pseudorapidity
  double phi;    ///< Polar angle
};

#endif /* Jet_hpp */
