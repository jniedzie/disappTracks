//
//  Jet.hpp
//
//  Created by Jeremi Niedziela on 17/07/2018.
//

#ifndef Jet_hpp
#define Jet_hpp

#include "JetCut.hpp"

/// This class represents a single jet. It contains information about its kinematical
/// properties, charged and neutral hadron fraction etc. It also checks if this jet passes some set
/// of selection criteria represented by JetCut object.
class Jet{
public:
  /// Default constructor
  Jet();
  
  /// Default destructor
  ~Jet();

  /// Print basic information about the jet
  void Print();
  
  /// Checks if jet passes selection criteria
  /// \param cut Jets selection criteria to be checked
  bool IsPassingCut(JetCut *cut);
  
  
  // Setters
  inline void SetPt(double val){pt = val;}
  inline void SetEta(double val){eta = val;}
  inline void SetPhi(double val){phi = val;}
  inline void SetMass(double val){mass = val;}
  inline void SetChargedHadronEnergyFraction(double val){chargedHadronEnergyFraction = val;}
  inline void SetNeutralHadronEnergyFraction(double val){neutralHadronEnergyFraction = val;}
  inline void SetIsForward(bool val){isForward = val;}
  
  // Getters
  inline double GetPt(){return pt;}
  inline double GetEta(){return eta;}
  inline double GetPhi(){return phi;}
  inline double GetMass(){return mass;}
  inline double GetChargedHadronEnergyFraction(){return chargedHadronEnergyFraction;}
  inline double GetNeutralHadronEnergyFraction(){return neutralHadronEnergyFraction;}
  inline bool   IsForward(){return isForward;}
  
private:
  double pt;                          ///< Transverse momentum
  double eta;                         ///< Pseudorapidity
  double phi;                         ///< Polar angle
  double mass;                        ///< Mass
  double chargedHadronEnergyFraction; ///< Energy fraction carried by charged hadrons
  double neutralHadronEnergyFraction; ///< Energy fraction carried by neutral hadrons
  bool   isForward;                   ///< is this a forward jet, or a regular one
};

#endif /* Jet_hpp */
