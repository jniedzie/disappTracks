//  LeptonCut.hpp
//
//  Created by Jeremi Niedziela on 07/08/2018.

#ifndef LeptonCut_hpp
#define LeptonCut_hpp

#include "Helpers.hpp"

/// Class containing definition of the lepton selection criteria.
/// User can define ranges of allowed parameters and required flags.
class LeptonCut {
public:
  /// Default constructor
  LeptonCut();
  
  /// Copy constructor
  LeptonCut(const LeptonCut &c);
  
  /// Default constructor
  ~LeptonCut();
  
  // setters
  inline void SetPt(range<double> val){pt=val;}
  inline void SetRelativeIsolation(range<double> val){relativeIsolation=val;}
  inline void SetRequireTightID(bool id){requiresTightID = id;}
  
private:
  range<double> pt;                 ///< allowed transverse momentum of the lepton
  range<double> relativeIsolation;  ///< allowed relative isolation in dR=0.4
  bool requiresTightID;             ///< should tight ID be required
  
  friend class LeptonProcessor;
};

#endif /* LeptonCut_hpp */
