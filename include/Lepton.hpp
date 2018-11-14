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
  
  /// Check if lepton passes selection criteria
  /// \param cut Lepton selection criteria to be checked
  bool IsPassingCut(const unique_ptr<LeptonCut> &cut);
  
  
  // Setters
  void SetPt(double val){pt = val;}
  void SetEta(double val){eta=val;}
  void SetPhi(double val){phi=val;}
  void SetRelativeIsolation(double val){relativeIsolation = val;}
  void SetTightID(bool val){tightID = val;}
  void SetPid(int val){pid = val;}
  
  
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
};


#endif /* Lepton_hpp */
