//
//  Lepton.hpp
//
//  Created by Jeremi Niedziela on 07/08/2018.
//

#ifndef Lepton_hpp
#define Lepton_hpp

#include "Helpers.hpp"
#include "LeptonCut.hpp"

/// This class represents a single lepton. It contains information about its kinematical
/// properties, isolation, tight ID flag etc. It also checks if this track passes some set of
/// selection criteria represented by TrackCut object.
class Lepton{
public:
  /// Default constructor
  Lepton();
  
  /// Default destructor
  ~Lepton(){};
  
  // Getters
  double  GetPt(){return pt;}
  double  GetEta(){return eta;}
  double  GetPhi(){return phi;}
  double  GetRelativeIsolation(){return relativeIsolation;}
  bool    GetTightID(){return tightID;}
  int     GetPid(){return pid;}
  
private:
  double pt;                  ///< Transverse momentum (GeV)
  double eta;                 ///< Pseudorapidity
  double phi;                 ///< Polar angle
  double relativeIsolation;   ///< Relative track isolation in cone dR=0.4
  bool tightID;               ///< Flag: does it pass tight ID criteria
  int pid;                    ///< Particle's PDG PID code
  
  friend class LeptonProcessor;
};


#endif /* Lepton_hpp */
