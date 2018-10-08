//
//  Jet.hpp
//
//  Created by Jeremi Niedziela on 17/07/2018.
//

#ifndef Jet_hpp
#define Jet_hpp

#include "JetCut.hpp"

class Jet{
public:
  Jet();
  ~Jet();

  inline void SetPt(double _pt){pt = _pt;}
  inline void SetEta(double _eta){eta = _eta;}
  inline void SetPhi(double _phi){phi = _phi;}
  inline void SetMass(double _mass){mass = _mass;}
  inline void SetChargedHadronEnergyFraction(double frac){chargedHadronEnergyFraction = frac;}
  inline void SetNeutralHadronEnergyFraction(double frac){neutralHadronEnergyFraction = frac;}
  inline void SetIsForward(double val){isForward = val;}
  
  inline double GetPt(){return pt;}
  inline double GetEta(){return eta;}
  inline double GetPhi(){return phi;}
  inline double GetMass(){return mass;}
  inline double GetChargedHadronEnergyFraction(){return chargedHadronEnergyFraction;}
  inline double GetNeutralHadronEnergyFraction(){return neutralHadronEnergyFraction;}
  inline bool   IsForward(){return isForward;}
  
  bool IsPassingCut(JetCut *cut);
  void Print();
  
private:
  double pt;     ///< Transverse momentum
  double eta;    ///< Pseudorapidity
  double phi;    ///< Polar angle
  double mass;   ///< Mass
  double chargedHadronEnergyFraction; ///< Energy fraction carried by charged hadrons
  double neutralHadronEnergyFraction; ///< Energy fraction carried by neutral hadrons
  bool isForward;   ///< is this a forward jet, or a regular one
};

#endif /* Jet_hpp */
