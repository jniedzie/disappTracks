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
  
  // Getters
  inline double GetPt(){return pt;}
  inline double GetEta(){return eta;}
  inline double GetPhi(){return phi;}
  inline double GetMass(){return mass;}
  inline double GetChHEF(){return chHEF;}
  inline double GetNeHEF(){return neHEF;}
  inline bool   IsForward(){return isForward;}
  
private:
  double pt;        ///< Transverse momentum
  double eta;       ///< Pseudorapidity
  double phi;       ///< Polar angle
  double mass;      ///< Mass
  double chHEF;     ///< Energy fraction carried by charged hadrons
  double neHEF;     ///< Energy fraction carried by neutral hadrons
  bool   isForward; ///< is this a forward jet, or a regular one
  
  friend class JetProcessor;
};

#endif /* Jet_hpp */
