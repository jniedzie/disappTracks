//
//  JetCut.hpp
//
//  Created by Jeremi Niedziela on 20/07/2018.
//

#ifndef JetCut_hpp
#define JetCut_hpp

#include "Helpers.hpp"

/// Class containing definition of the jet selection criteria.
/// User can define ranges of allowed parameters and required flags.
class JetCut {
public:
  /// Default constructor
  JetCut();
  
  /// Default desctructor
  ~JetCut();
  
  // Setters
  inline void SetPt(range<double> val){pt=val;}
  inline void SetEta(range<double> val){eta=val;}
  inline void SetEtaForward(range<double> val){etaForward=val;}
  inline void SetChargedHadronEnergyFraction(range<double> val){chargedHadronEnergyFraction=val;}
  inline void SetNeutralHadronEnergyFraction(range<double> val){neutralHadronEnergyFraction=val;}
  inline void SetTrackDeltaR(range<double> val){trackDeltaR=val;}
  
  // Getters
  inline range<double> GetPt(){return pt;}
  inline range<double> GetEta(){return eta;}
  inline range<double> GetEtaForward(){return etaForward;}
  inline range<double> GetChargedHadronEnergyFraction(){return chargedHadronEnergyFraction;}
  inline range<double> GetNeutralHadronEnergyFraction(){return neutralHadronEnergyFraction;}
  inline range<double> GetTrackDeltaR(){return trackDeltaR;}
  
private:
  range<double> pt;                           ///< allowed transverse momentum of the jet
  range<double> eta;                          ///< allowed pseudorapidity
  range<double> etaForward;                   ///< allowed pseudorapidity for forward jets
  range<double> chargedHadronEnergyFraction;  ///< allowed charged hadron energy fraction
  range<double> neutralHadronEnergyFraction;  ///< allowed neutral hadron energy fraction
  range<double> trackDeltaR;                  ///< allowed separation with any of the tracks
};

#endif /* JetCut_hpp */
