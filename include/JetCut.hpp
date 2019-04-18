//  JetCut.hpp
//
//  Created by Jeremi Niedziela on 20/07/2018.

#ifndef JetCut_hpp
#define JetCut_hpp

#include "Helpers.hpp"

/// Class containing definition of the jet selection criteria.
/// User can define ranges of allowed parameters and required flags.
class JetCut {
public:
  /// Default constructor
  JetCut();
  
  /// Copy constructor
  JetCut(const JetCut &c);
  
  /// Default desctructor
  ~JetCut();
  
  // Setters
  inline void SetPt(range<double> val){pt=val;}
  inline void SetEta(range<double> val){eta=val;}
  inline void SetEtaForward(range<double> val){etaForward=val;}
  inline void SetChargedHadronEnergyFraction(range<double> val){ChHEF=val;}
  inline void SetNeutralHadronEnergyFraction(range<double> val){NeHEF=val;}
  
private:
  range<double> pt;         ///< allowed transverse momentum of the jet
  range<double> eta;        ///< allowed pseudorapidity
  range<double> etaForward; ///< allowed pseudorapidity for forward jets
  range<double> ChHEF;      ///< allowed charged hadron energy fraction
  range<double> NeHEF;      ///< allowed neutral hadron energy fraction
  
  friend class JetProcessor;
};

#endif /* JetCut_hpp */
